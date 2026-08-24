#!/usr/bin/env python
"""Phase 2: one (family, protonation state, connectivity) ladder cell.

Runs the proven Phase-1 stage machinery on a crystallographic cluster
cut from the kaolinite deck cell (HANDOFF.md §2b):

    uv run python scripts/phase2_ladder.py --family oss --gpu
    uv run python scripts/phase2_ladder.py --family oss --n-intact 2 --gpu
    uv run python scripts/phase2_ladder.py --family oaa --dry-run

Families: oss (300s Si-O-Si), osa (400s Si-O-Al2), oaa (500s Al-OH-Al).
Each run is one ladder cell: hydrolysis of the center bridge by the
attacking species, with the attacked metal's connectivity set by
``--n-intact`` (Q-style: intact bridges, center included). Stages are
checkpointed as XYZ under runs/phase2/<cell>/ — a crashed campaign
resumes where it stopped. Frozen-shell clusters get partial-Hessian
thermochemistry automatically (rot/trans cancel between complex and TS).

``--dry-run`` builds the cluster + attack complex, writes geometry and
metadata, and exits before any DFT — the sizing step.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    bootstrap_cli(
        "phase2_ladder",
        default_run_root=Path(__file__).resolve().parent.parent / "runs",
    )

import numpy as np  # noqa: E402

from quarry.clusters import Cluster, hydronium, water  # noqa: E402
from quarry.crystal import attack_complex, from_deck_cell  # noqa: E402
from quarry.pipeline import (  # noqa: E402
    HARTREE_TO_KJ,
    DftSettings,
    energy,
    frequencies,
    optimize,
)
from quarry.rates import rate_from_thermo, thermo_from_frequencies  # noqa: E402
from quarry.store import Store, geometry_hash  # noqa: E402
from quarry.ts import (  # noqa: E402
    ScanNoMaximumError,
    constrained_scan,
    find_ts,
    first_interior_maximum,
    quick_irc,
    scan_to_maximum,
    scan_ts_guess,
)

FAMILIES = {"oss": "Oss", "osa": "Osa", "oaa": "Oaa"}
STATES = {"neutral": (water, 0), "acid": (hydronium, +1)}
KCAL = 4.184
# Imaginary modes below this are PHVA numerical noise, not reaction modes.
NOISE_FLOOR_CM = 30.0
# A "barrier" below this is a trivial-rearrangement saddle (learnings).
MIN_PLAUSIBLE_DE_KJ = 20.0
# The Phase-1 free-dimer anchor for the lattice-resistance comparison.
SI_NEUTRAL_FREE_DIMER_DG_KJ = 113.05
APPROACH_SEED_VERSION = 1

# Per-element approach parameters (Angstrom).
APPROACH = {
    "Si": {
        "distances": [2.8, 2.6, 2.4, 2.2, 2.1, 2.0, 1.9],
        "pin": 1.90,
        "limit": 2.6,
    },
    "Al": {
        "distances": [3.0, 2.8, 2.6, 2.4, 2.3, 2.2, 2.1, 2.0],
        "pin": 2.00,
        "limit": 2.8,
    },
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def preload_cutensor() -> None:
    """Load the cutensor-cu12 wheel's libraries with RTLD_GLOBAL.

    The wheel lands in site-packages/cutensor/lib, which is not on the
    dynamic loader's path; cupy then silently falls back to its einsum
    contraction engine, which materializes temporaries that OOM a 24 GB
    card at ~60 atoms/def2-svp (seen live, twice). Preloading by
    absolute path makes cupy's dlopen find them regardless of
    LD_LIBRARY_PATH. The core library is mandatory for GPU runs; the
    multi-GPU companion is optional.
    """
    import contextlib
    import ctypes
    import sysconfig

    libdir = Path(sysconfig.get_paths()["purelib"]) / "cutensor" / "lib"
    core = libdir / "libcutensor.so.2"
    try:
        ctypes.CDLL(str(core), mode=ctypes.RTLD_GLOBAL)
    except OSError as exc:
        raise RuntimeError(
            f"failed to load required cuTENSOR core library {core}: {exc}; "
            "install cutensor-cu12 in the active environment before using --gpu"
        ) from exc

    with contextlib.suppress(OSError):
        ctypes.CDLL(str(libdir / "libcutensorMg.so.2"), mode=ctypes.RTLD_GLOBAL)


def trim_gpu_pool() -> None:
    """Return CuPy pool blocks to the CUDA driver between stages.

    Multi-hour campaigns otherwise fragment the pool into
    cudaErrorMemoryAllocation (seen live: stage-2 DF-K gradient OOM
    right after a 2 h stage-1 optimization on a 24 GB card with ~7 GB
    resident). No-op without a GPU.
    """
    try:
        import cupy

        cupy.get_default_memory_pool().free_all_blocks()
        cupy.get_default_pinned_memory_pool().free_all_blocks()
    except Exception:
        pass


def save_xyz(cluster: Cluster, path: Path) -> None:
    path.write_text(cluster.to_xyz())


def load_xyz(path: Path, template: Cluster) -> Cluster:
    from dataclasses import replace

    lines = path.read_text().splitlines()
    n = int(lines[0])
    coords = np.array(
        [[float(x) for x in line.split()[1:4]] for line in lines[2 : 2 + n]]
    )
    return replace(template, coords=coords)


def approach_seed_signature(
    complex_guess: Cluster,
    *,
    m_index: int,
    br_index: int,
    ow_index: int,
    pin_a: float,
) -> dict[str, object]:
    """Parameters that make an approach checkpoint safe to resume."""
    return {
        "version": APPROACH_SEED_VERSION,
        "complex_guess_geometry_hash": geometry_hash(complex_guess.to_xyz()),
        "m_index": m_index,
        "br_index": br_index,
        "ow_index": ow_index,
        "pin_a": pin_a,
    }


def save_approach_seed(seed: Cluster, path: Path, signature: dict[str, object]) -> None:
    """Write an approach checkpoint and its compatibility data."""
    save_xyz(seed, path)
    path.with_suffix(".json").write_text(json.dumps(signature, indent=2))


def load_compatible_approach_seed(
    path: Path, template: Cluster, expected: dict[str, object]
) -> Cluster | None:
    """Load a seed only when its recorded inputs match this run exactly."""
    signature_path = path.with_suffix(".json")
    if not path.exists():
        return None
    if not signature_path.exists():
        log("  ignoring approach_seed.xyz: compatibility metadata is missing")
        return None
    try:
        actual = json.loads(signature_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        log(f"  ignoring approach_seed.xyz: invalid compatibility metadata ({exc})")
        return None
    if actual != expected:
        log("  ignoring approach_seed.xyz: parameters or input geometry changed")
        return None
    try:
        return load_xyz(path, template)
    except (OSError, ValueError) as exc:
        log(f"  ignoring approach_seed.xyz: checkpoint is unreadable ({exc})")
        return None


def checkpointed(path: Path, template: Cluster, compute) -> Cluster:
    if path.exists():
        log(f"  resume: {path.name} exists, skipping recompute")
        return load_xyz(path, template)
    result = compute()
    save_xyz(result, path)
    return result


def channel_escape_reason(
    ts_guess: Cluster, ts: Cluster, m_index: int, ow_index: int, limit_a: float
) -> str | None:
    r_guess = float(
        np.linalg.norm(ts_guess.coords[m_index] - ts_guess.coords[ow_index])
    )
    r_ts = float(np.linalg.norm(ts.coords[m_index] - ts.coords[ow_index]))
    if r_ts > limit_a:
        return f"saddle r(M-Ow)={r_ts:.2f} A is outside the bonding channel"
    if r_ts > r_guess + 0.5:
        return f"saddle escaped the approach channel by {r_ts - r_guess:.2f} A"
    return None


def proton_neb_guess(
    approach_seed: Cluster,
    complex_opt: Cluster,
    settings: DftSettings,
    run_dir: Path,
    *,
    m_index: int,
    br_index: int,
    ow_index: int,
    pin_a: float,
) -> Cluster:
    """Concerted proton-transfer/product/CI-NEB guess (Phase-1 route,
    parameterized for arbitrary cluster indices)."""
    product_path = run_dir / "product.xyz"
    product_seed: Cluster | None = None
    h_candidates = _attacker_h_indices(approach_seed, ow_index)

    if product_path.exists():
        product = load_xyz(product_path, approach_seed)
        r_prod_ow = float(
            np.linalg.norm(product.coords[m_index] - product.coords[ow_index])
        )
        r_prod_br = float(
            np.linalg.norm(product.coords[m_index] - product.coords[br_index])
        )
        log(
            "  resume: product.xyz exists; "
            f"r(M-Ow)={r_prod_ow:.2f} A, r(M-Obr)={r_prod_br:.2f} A"
        )
        if r_prod_ow <= pin_a and r_prod_br >= 2.2:
            return crest_from_product(
                product,
                settings,
                m_index=m_index,
                br_index=br_index,
                ow_index=ow_index,
                pin_a=pin_a,
            )
        log("  saved product is not hydrolyzed; extending M-Obr cleavage")
        product_path.replace(run_dir / "product.rejected-rollback.xyz")
        product_seed = product

    if product_seed is None:
        seed_r = float(
            np.linalg.norm(
                approach_seed.coords[m_index] - approach_seed.coords[ow_index]
            )
        )
        if abs(seed_r - pin_a) > 0.02:
            log(
                f"  pinning resumed approach seed: "
                f"r(M-Ow) {seed_r:.2f} -> {pin_a:.2f} A"
            )
            approach_seed = constrained_scan(
                approach_seed,
                settings,
                atom_i=m_index,
                atom_j=ow_index,
                distances_a=[pin_a],
            )[0][2]

        log(f"  stage 2b: pin r(M-Ow)={pin_a:.2f} A, drive H(water) -> O(bridge)")
        h_idx = min(
            h_candidates,
            key=lambda i: np.linalg.norm(
                approach_seed.coords[i] - approach_seed.coords[br_index]
            ),
        )
        r0 = float(
            np.linalg.norm(approach_seed.coords[h_idx] - approach_seed.coords[br_index])
        )
        log(f"  driving H{h_idx} from r(Obr-H)={r0:.2f} A inward")
        first = max(1.05, round(r0 - 0.15, 2))
        distances = [round(float(x), 2) for x in np.arange(first, 1.04, -0.15)]
        try:
            pscan = scan_to_maximum(
                approach_seed,
                settings,
                atom_i=br_index,
                atom_j=h_idx,
                distances_a=distances or [first],
                fixed_distances=[(m_index, ow_index, pin_a)],
                extend_step_a=0.06,
                min_distance_a=0.95,
                progress=lambda r, e: (
                    log(f"  r(Obr-H)={r:.2f} A  E={e:.6f} Ha"),
                    trim_gpu_pool(),
                )[0],
            )
        except ScanNoMaximumError as exc:
            # Embedded clusters: proton transfer can be coupled to bridge
            # rupture, so no crest exists along the proton coordinate
            # alone (seen live on oss-neutral-n4-s2). The floor endpoint
            # has the proton delivered — a valid product seed; CI-NEB is
            # pinned to both basins and finds the concerted col anyway.
            log(
                "  proton scan monotonic to floor — transfer is coupled "
                "to bridge rupture; using delivered-proton endpoint"
            )
            pscan = exc.scan
            product_seed = pscan[-1][2]
        else:
            crest_idx = first_interior_maximum(pscan)
            if crest_idx is None:
                raise RuntimeError("proton scan produced no interior crest")
            product_seed = pscan[min(crest_idx + 1, len(pscan) - 1)][2]

    h_idx = min(
        h_candidates,
        key=lambda i: np.linalg.norm(
            product_seed.coords[i] - product_seed.coords[br_index]
        ),
    )
    current_br = float(
        np.linalg.norm(product_seed.coords[m_index] - product_seed.coords[br_index])
    )
    break_targets = [
        r for r in [1.85, 2.05, 2.30, 2.60, 3.00, 3.40, 3.60] if r > current_br + 0.05
    ]
    if not break_targets:
        break_targets = [round(current_br + 0.30, 2)]

    log("  stage 2c: breaking M-Obr with the new bonds held")
    broken = constrained_scan(
        product_seed,
        settings,
        atom_i=m_index,
        atom_j=br_index,
        distances_a=break_targets,
        fixed_distances=[
            (m_index, ow_index, pin_a - 0.20),
            (br_index, h_idx, 0.98),
        ],
    )
    for r, e, _ in broken:
        log(f"  r(M-Obr)={r:.2f} A  E={e:.6f} Ha")
    product = optimize(broken[-1][2], settings)
    save_xyz(product, product_path)

    r_prod_ow = float(
        np.linalg.norm(product.coords[m_index] - product.coords[ow_index])
    )
    r_prod_br = float(
        np.linalg.norm(product.coords[m_index] - product.coords[br_index])
    )
    log(f"  product r(M-Ow)={r_prod_ow:.2f} A, r(M-Obr)={r_prod_br:.2f} A")
    if r_prod_ow > pin_a or r_prod_br < 2.2:
        raise RuntimeError(
            "product rolled back toward reactants; expected M-Ow bonded "
            "and M-Obr broken"
        )
    return crest_from_product(
        product,
        settings,
        m_index=m_index,
        br_index=br_index,
        ow_index=ow_index,
        pin_a=pin_a,
    )


def crest_from_product(
    product: Cluster,
    settings: DftSettings,
    *,
    m_index: int,
    br_index: int,
    ow_index: int,
    pin_a: float,
) -> Cluster:
    """TS guess: drive the broken bridge back inward to the interior crest.

    A CI-NEB on a 60+-atom embedded cluster costs ~5 min per band step
    (measured live: 80 MDMin steps = 7 h, still 0.43 eV/A) — but the
    stage-2c rupture scan already proves an interior maximum exists
    along r(M-Obr) with the new bonds held. Re-cross the col from the
    product side with the same constrained-scan machinery; Sella and
    the in-channel/verify gates then judge the crest exactly as they
    would a NEB peak. neb_ts_guess remains the escaped-saddle fallback.
    """
    r_br = float(np.linalg.norm(product.coords[m_index] - product.coords[br_index]))
    first = round(r_br - 0.20, 2)
    if first < 1.70:
        raise ValueError(
            f"product M-Obr bond ({r_br:.2f} A) is too short for a reverse crest scan"
        )
    distances = [round(r, 2) for r in np.arange(r_br - 0.20, 1.99, -0.15)]
    if not distances:
        distances = [first]
    log("  stage 2d: reverse crest scan r(M-Obr) from the product side")
    scan = scan_to_maximum(
        product,
        settings,
        atom_i=m_index,
        atom_j=br_index,
        distances_a=distances,
        fixed_distances=[
            (m_index, ow_index, pin_a - 0.20),
            (br_index, _nearest_h(product, br_index, ow_index), 0.98),
        ],
        extend_step_a=0.10,
        min_distance_a=1.70,
        progress=lambda r, e: (
            log(f"  r(M-Obr)={r:.2f} A  E={e:.6f} Ha"),
            trim_gpu_pool(),
        )[0],
    )
    return scan_ts_guess(scan)


def _attacker_h_indices(cluster: Cluster, ow_index: int) -> list[int]:
    """Return appended attacker H indices, failing fast if ordering changed.

    ``attack_complex`` appends water/hydronium with its oxygen first. This
    builder contract preserves proton identity after transfer, when the
    delivered proton is intentionally no longer bonded to the water oxygen.
    """
    if not 0 <= ow_index < len(cluster.symbols) or cluster.symbols[ow_index] != "O":
        raise ValueError(f"attacker oxygen index {ow_index} is not an O atom")
    h_indices = list(range(ow_index + 1, len(cluster.symbols)))
    if not h_indices or any(cluster.symbols[i] != "H" for i in h_indices):
        raise ValueError(
            "attacker atom ordering changed: expected only H atoms after "
            f"oxygen index {ow_index}"
        )
    return h_indices


def _nearest_h(cluster: Cluster, br_index: int, ow_index: int) -> int:
    """The attacker proton closest to the bridge oxygen."""
    return min(
        _attacker_h_indices(cluster, ow_index),
        key=lambda i: np.linalg.norm(cluster.coords[i] - cluster.coords[br_index]),
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--family", choices=sorted(FAMILIES), required=True)
    ap.add_argument("--state", choices=sorted(STATES), default="neutral")
    ap.add_argument("--n-intact", type=int, default=None)
    ap.add_argument("--metal-shells", type=int, default=2)
    ap.add_argument("--center-index", type=int, default=None)
    ap.add_argument(
        "--deck",
        default=None,
        help="petra deck carrying the [cell] (default: repo kaolinite.toml)",
    )
    ap.add_argument("--xc", default="b3lyp")
    ap.add_argument("--basis", default="def2-svp")
    ap.add_argument("--gpu", action="store_true")
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--nice", type=int, default=10, help="process niceness increment")
    ap.add_argument("--log", help="tee output to this path (default: qm/runs)")
    ap.add_argument("--temperature", type=float, default=298.15)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument(
        "--run-root",
        default=None,
        help="parent of the run dir (default: qm/runs) — tests redirect this",
    )
    args = ap.parse_args()

    os.environ.setdefault("OMP_NUM_THREADS", str(args.threads))
    if args.gpu:
        preload_cutensor()
    qm_root = Path(__file__).resolve().parent.parent
    deck = (
        Path(args.deck)
        if args.deck
        else (qm_root.parent / "petra" / "examples" / "kaolinite.toml")
    )
    settings = DftSettings(
        xc=args.xc, basis=args.basis, density_fit=True, use_gpu=args.gpu
    )

    attacker_factory, charge_offset = STATES[args.state]
    site_kind = FAMILIES[args.family]

    cc = from_deck_cell(
        deck,
        site_kind,
        center_index=args.center_index,
        metal_shells=args.metal_shells,
        n_intact=args.n_intact,
        target_charge=0,
    )
    attacker = attacker_factory()
    complex_guess, ow_index = attack_complex(cc, attacker)
    m_index, br_index = cc.attacked_index, cc.bridge_index
    assert br_index is not None  # ladder families always attack a bridge oxygen
    metal = cc.cluster.symbols[m_index]
    approach = APPROACH[metal]

    cell_name = f"{args.family}-{args.state}-n{cc.n_intact}-s{cc.metal_shells}"
    run_root = Path(args.run_root) if args.run_root else qm_root / "runs"
    run_dir = run_root / "phase2" / f"{cell_name}-{args.xc}-{args.basis}"
    run_dir.mkdir(parents=True, exist_ok=True)
    log(f"run dir: {run_dir}")
    log(f"settings: {settings}")
    log(
        f"cell: {cell_name}  cluster {cc.cluster.formula} "
        f"({len(cc.cluster.symbols)} atoms, {len(cc.cluster.frozen_indices)} frozen, "
        f"charge {cc.cluster.charge + attacker.charge}), attacked {metal}@{m_index}, "
        f"bridge O@{br_index}"
    )
    for line in cc.termination_log:
        log(f"  termination: {line}")

    save_xyz(cc.cluster, run_dir / "cluster.xyz")
    save_xyz(complex_guess, run_dir / "complex_guess.xyz")
    metadata = cc.metadata()
    metadata["state"] = args.state
    metadata["attacker"] = attacker.name
    metadata["charge_offset"] = charge_offset
    metadata["deck"] = str(deck)
    metadata["method"] = f"{args.xc}/{args.basis}/df"
    metadata["gpu"] = args.gpu
    metadata["temperature_k"] = args.temperature
    metadata["written_at"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    try:
        import subprocess

        metadata["driver_git_commit"] = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=Path(__file__).parent,
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    except Exception:
        metadata["driver_git_commit"] = None
    (run_dir / "metadata.json").write_text(json.dumps(metadata, indent=2))
    if args.dry_run:
        log("dry run: geometry + metadata written, no DFT")
        return 0

    # Stage 1 — optimize reactant complex and separated fragments.
    # Crystallographic positions are sane starting points; no cheap
    # pre-opt stage (the crude legacy cell is a logged systematic).
    log("stage 1: optimizing reactant complex + fragments")
    complex_opt = checkpointed(
        run_dir / "complex.xyz",
        complex_guess,
        lambda: optimize(complex_guess, settings),
    )
    cluster_opt = checkpointed(
        run_dir / "cluster_opt.xyz",
        cc.cluster,
        lambda: optimize(cc.cluster, settings),
    )
    attacker_opt = checkpointed(
        run_dir / "attacker.xyz", attacker, lambda: optimize(attacker, settings)
    )
    trim_gpu_pool()

    # Stage 2 — relaxed approach scan, then the concerted proton route
    # if the approach coordinate alone never crosses the ridge.
    log(f"stage 2: relaxed scan r({metal}-Ow), auto-extending to interior maximum")
    ts_guess_path = run_dir / "ts_guess.xyz"
    route_path = run_dir / "ts_guess.route"
    approach_seed_path = run_dir / "approach_seed.xyz"
    seed_signature = approach_seed_signature(
        complex_guess,
        m_index=m_index,
        br_index=br_index,
        ow_index=ow_index,
        pin_a=approach["pin"],
    )
    route = route_path.read_text().strip() if route_path.exists() else "direct"
    if ts_guess_path.exists():
        ts_guess = load_xyz(ts_guess_path, complex_opt)
        log(f"  resume: ts_guess.xyz exists ({route} route)")
    else:
        route = "direct"
        route_path.unlink(missing_ok=True)
        base: Cluster | None = None
        base = load_compatible_approach_seed(
            approach_seed_path, complex_opt, seed_signature
        )
        if base is not None:
            # A crashed proton-route attempt already proved the direct
            # route dead — skip the approach scan entirely.
            log("  resume: approach_seed.xyz exists, skipping approach scan")

        else:
            try:
                scan = scan_to_maximum(
                    complex_opt,
                    settings,
                    atom_i=m_index,
                    atom_j=ow_index,
                    distances_a=approach["distances"],
                    progress=lambda r, e: (
                        log(f"  r={r:.2f} A  E={e:.6f} Ha"),
                        trim_gpu_pool(),
                    )[0],
                )
                ts_guess = scan_ts_guess(scan)
            except ScanNoMaximumError as exc:
                log("  approach coordinate alone does not cross the ridge;")
                base = min(exc.scan, key=lambda p: abs(p[0] - approach["pin"]))[2]
                save_approach_seed(base, approach_seed_path, seed_signature)
        if base is not None:
            route = "proton-neb"
            route_path.write_text(route)
            ts_guess = proton_neb_guess(
                base,
                complex_opt,
                settings,
                run_dir,
                m_index=m_index,
                br_index=br_index,
                ow_index=ow_index,
                pin_a=approach["pin"],
            )
        save_xyz(ts_guess, ts_guess_path)
        if route == "proton-neb":
            route_path.write_text(route)

    # Stage 3 — Sella saddle search with the escaped-channel gates.
    trim_gpu_pool()
    log("stage 3: Sella saddle search")
    ts_path = run_dir / "ts.xyz"
    trajectory_path = run_dir / "sella.traj"
    ts = checkpointed(
        ts_path,
        ts_guess,
        lambda: find_ts(ts_guess, settings, trajectory=str(trajectory_path)),
    )
    escape = channel_escape_reason(ts_guess, ts, m_index, ow_index, approach["limit"])
    if escape and route == "direct":
        log(f"  !! {escape}; pivoting once to proton-drive + product + CI-NEB")
        neb_guess = proton_neb_guess(
            ts_guess,
            complex_opt,
            settings,
            run_dir,
            m_index=m_index,
            br_index=br_index,
            ow_index=ow_index,
            pin_a=approach["pin"],
        )
        if ts_path.exists():
            ts_path.replace(run_dir / "ts.rejected-direct.xyz")
        if trajectory_path.exists():
            trajectory_path.replace(run_dir / "sella.rejected-direct.traj")
        ts_guess = neb_guess
        save_xyz(ts_guess, ts_guess_path)
        route = "proton-neb"
        route_path.write_text(route)
        log("stage 3b: Sella saddle search from CI-NEB guess")
        ts = checkpointed(
            ts_path,
            ts_guess,
            lambda: find_ts(ts_guess, settings, trajectory=str(trajectory_path)),
        )
        escape = channel_escape_reason(
            ts_guess, ts, m_index, ow_index, approach["limit"]
        )
    if escape:
        log(f"  !! {escape}; {route} route exhausted — aborting before thermochemistry")
        return 1

    # Stage 4 — verify (PHVA-aware) + quick IRC.
    trim_gpu_pool()
    log("stage 4: frequencies + quick-IRC")
    ts_freq = frequencies(ts, settings)
    imag = ts_freq.imaginary_cm
    significant = imag[imag > NOISE_FLOOR_CM]
    log(f"  TS imaginary modes: {np.round(imag, 1).tolist()}")
    if significant.size != 1:
        log(
            f"  !! expected exactly 1 imaginary mode above {NOISE_FLOOR_CM:.0f} "
            f"cm^-1, found {significant.size} — inspect ts.xyz; aborting"
        )
        return 1
    back, fwd = quick_irc(ts, settings, noise_floor_cm=NOISE_FLOOR_CM)
    save_xyz(back, run_dir / "irc_back.xyz")
    save_xyz(fwd, run_dir / "irc_fwd.xyz")

    # Stage 5 — thermochemistry (PHVA on frozen clusters; trans/rot
    # cancel between complex and TS of identical composition).
    trim_gpu_pool()
    log("stage 5: thermochemistry")
    cx_freq = frequencies(complex_opt, settings)
    t = args.temperature

    def thermo(freq):
        return thermo_from_frequencies(
            freq.electronic_hartree * HARTREE_TO_KJ,
            freq.frequencies_cm,
            t,
            molar_mass_kg=freq.molar_mass_kg,
            rotational_temperatures_k=(
                list(freq.rotational_temperatures_k)
                if freq.rotational_temperatures_k
                else None
            ),
            linear=freq.linear,
        )

    rate = rate_from_thermo(
        thermo(cx_freq),
        thermo(ts_freq),
        imag_nu_cm=float(significant[0]),
        tunneling="wigner",
    )
    de_kj = (ts_freq.electronic_hartree - cx_freq.electronic_hartree) * HARTREE_TO_KJ
    if de_kj < MIN_PLAUSIBLE_DE_KJ:
        log(
            f"  !! dE(elec)={de_kj:.1f} kJ/mol is below the trivial-saddle "
            f"floor ({MIN_PLAUSIBLE_DE_KJ:.0f}) — chemically wrong saddle; aborting"
        )
        return 1
    frag_e = (
        ts_freq.electronic_hartree
        - energy(cluster_opt, settings)
        - energy(attacker_opt, settings)
    ) * HARTREE_TO_KJ

    results = {
        "cell": cell_name,
        "family": args.family,
        "state": args.state,
        "n_intact": cc.n_intact,
        "metal_shells": cc.metal_shells,
        "attacked_metal": metal,
        "method": f"{args.xc}/{args.basis}/df",
        "temperature_k": t,
        "route": route,
        "dE_elec_vs_complex_kj": de_kj,
        "dE_elec_vs_fragments_kj": frag_e,
        "dH_kj": rate.dh_kj,
        "dG_kj": rate.dg_kj,
        "dG_kcal": rate.dg_kj / KCAL,
        "wigner_kappa": rate.kappa,
        "k_per_s": rate.k,
        "ts_imaginary_cm": float(significant[0]),
        "ts_imaginary_all_cm": [float(x) for x in imag],
        "cluster": metadata,
        "geometry_hash": {
            "complex": geometry_hash(complex_opt.to_xyz()),
            "ts": geometry_hash(ts.to_xyz()),
        },
    }
    if args.family == "oss" and args.state == "neutral":
        shift = rate.dg_kj - SI_NEUTRAL_FREE_DIMER_DG_KJ
        results["lattice_resistance_shift_kj"] = shift
        log(
            f"  LATTICE RESISTANCE: dG‡ {rate.dg_kj:.1f} vs free-dimer "
            f"{SI_NEUTRAL_FREE_DIMER_DG_KJ:.1f} kJ/mol -> shift {shift:+.1f} kJ/mol "
            "(publishable observable — Pelmenschikov predicts +20 to +70)"
        )
    (run_dir / "results.json").write_text(json.dumps(results, indent=2))

    with Store(run_dir / "store.sqlite") as store:
        for name, cl, freq in (
            ("complex", complex_opt, cx_freq),
            ("ts", ts, ts_freq),
        ):
            sid = store.add_structure(
                f"{cell_name}-{name}",
                cl.formula,
                cl.to_xyz(),
                charge=cl.charge,
                spin=cl.spin,
            )
            jid = store.add_job(
                sid, "freq", results["method"], "gpu4pyscf" if args.gpu else "pyscf"
            )
            store.set_job_status(jid, "done")
            store.add_result(jid, "electronic", freq.electronic_hartree, "hartree")

    log("")
    log(f"=== {cell_name} @ {results['method']} ===")
    log(f"  dE(elec, TS - complex)   = {de_kj:8.1f} kJ/mol")
    log(f"  dH‡({t:.0f})              = {rate.dh_kj:8.1f} kJ/mol")
    log(
        f"  dG‡({t:.0f})              = {rate.dg_kj:8.1f} kJ/mol "
        f"({rate.dg_kj / KCAL:.1f} kcal/mol)"
    )
    log(f"  wigner kappa             = {rate.kappa:8.2f}")
    log(f"  k({t:.0f})                = {rate.k:8.3e} 1/s")
    log(f"  TS imag mode             = {float(significant[0]):8.1f}i cm^-1")
    log(f"  dE(elec, TS - fragments) = {frag_e:6.1f} kJ/mol")
    log(f"results.json + store.sqlite written to {run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
