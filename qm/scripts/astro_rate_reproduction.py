#!/usr/bin/env python
"""D2a: reproduce low-temperature astrochemical tunneling rates.

The campaign uses tiny gas-phase models as an intentionally cheap calibration
of quarry's stationary-point and asymmetric-Eckart machinery.  Every expensive
stage is checkpointed under ``qm/runs/D2a-astro-rate-reproduction``.

Run on the workstation only, with the command recorded in ``run-command.txt``::

    OMP_NUM_THREADS=16 uv run python scripts/astro_rate_reproduction.py --gpu

The published comparisons are Langmuir-Hinshelwood (unimolecular) rate fits for
surface reactions.  Gas-phase quarry models are therefore a method-transfer
stress test, not a claim to reproduce the authors' QM/MM potential exactly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import time
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

from quarry.clusters import Cluster
from quarry.pipeline import (
    HARTREE_TO_KJ,
    DftSettings,
    FrequencyResult,
    frequencies,
)
from quarry.rates import eckart_kappa, thermo_from_frequencies
from quarry.ts import (
    ScanNoMaximumError,
    find_ts,
    quick_irc,
    scan_to_maximum,
    scan_ts_guess,
)

TEMPERATURES = (50.0, 60.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0)
RUN_ROOT = (
    Path(__file__).resolve().parent.parent / "runs" / "D2a-astro-rate-reproduction"
)


@dataclass(frozen=True)
class LiteratureFit:
    alpha_s: float
    beta: float
    gamma_k: float
    t0_k: float
    source: str
    valid_floor_k: float

    def rate(self, temperature_k: float) -> float:
        t = temperature_k
        return (
            self.alpha_s
            * (t / 300.0) ** self.beta
            * np.exp(-self.gamma_k * (t + self.t0_k) / (t * t + self.t0_k * self.t0_k))
        )


@dataclass(frozen=True)
class Reaction:
    key: str
    label: str
    cluster: Cluster
    scan_i: int
    scan_j: int
    scan_distances_a: tuple[float, ...]
    scan_floor_a: float
    literature_barrier_k: float | None
    literature_imag_cm: float | None
    literature_fit: LiteratureFit | None
    method: DftSettings
    barrierless: bool = False


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def _formaldehyde(incoming: np.ndarray, *, key: str) -> Cluster:
    # Planar H2CO seed: C, O, Ha, Hb, incoming H.
    return Cluster(
        name=key,
        symbols=["C", "O", "H", "H", "H"],
        coords=np.array(
            [
                [0.0, 0.0, 0.0],
                [1.205, 0.0, 0.0],
                [-0.575, 0.935, 0.0],
                [-0.575, -0.935, 0.0],
                incoming,
            ]
        ),
        spin=1,
    )


def reactions(*, gpu: bool, basis: str) -> dict[str, Reaction]:
    pwb6k = DftSettings(
        xc="pwb6k", basis=basis, dispersion="d3bj", density_fit=True, use_gpu=gpu
    )
    bhlyp = DftSettings(xc="bhandhlyp", basis=basis, density_fit=True, use_gpu=gpu)
    song = "Song & Kaestner 2020, arXiv:2009.05442, Table 2/Table 4"
    meisner = "Meisner, Lamberts & Kaestner 2017, arXiv:1708.05559 + SI"
    return {
        "h-h2co-ch3o": Reaction(
            key="h-h2co-ch3o",
            label="H + H2CO -> CH3O",
            cluster=_formaldehyde(np.array([0.0, 0.0, 3.0]), key="h-h2co-ch3o"),
            scan_i=0,
            scan_j=4,
            scan_distances_a=(2.8, 2.5, 2.2, 2.0, 1.85, 1.7, 1.55, 1.4, 1.25),
            scan_floor_a=1.1,
            literature_barrier_k=1900.0,
            literature_imag_cm=831.0,
            literature_fit=LiteratureFit(3146e10, 1.0, 830.0, 119.6, song, 59.0),
            method=pwb6k,
        ),
        "h-h2co-ch2oh": Reaction(
            key="h-h2co-ch2oh",
            label="H + H2CO -> CH2OH",
            cluster=_formaldehyde(np.array([3.0, 0.0, 0.0]), key="h-h2co-ch2oh"),
            scan_i=1,
            scan_j=4,
            scan_distances_a=(2.8, 2.5, 2.2, 1.95, 1.75, 1.55, 1.4, 1.25, 1.1),
            scan_floor_a=0.95,
            literature_barrier_k=5670.0,
            literature_imag_cm=1637.0,
            literature_fit=LiteratureFit(1.14e10, 1.0, 3146.0, 219.3, song, 75.0),
            method=pwb6k,
        ),
        "h-h2co-h2-hco": Reaction(
            key="h-h2co-h2-hco",
            label="H + H2CO -> H2 + HCO",
            cluster=_formaldehyde(np.array([-0.575, 0.935, 2.8]), key="h-h2co-h2-hco"),
            scan_i=2,
            scan_j=4,
            scan_distances_a=(2.6, 2.3, 2.0, 1.75, 1.5, 1.3, 1.15, 1.02, 0.9),
            scan_floor_a=0.78,
            literature_barrier_k=3030.0,
            literature_imag_cm=1663.0,
            literature_fit=LiteratureFit(4.13e10, 1.0, 1222.0, 147.7, song, 59.0),
            method=pwb6k,
        ),
        "oh-h2": Reaction(
            key="oh-h2",
            label="OH + H2 -> H2O + H",
            cluster=Cluster(
                name="oh-h2",
                symbols=["O", "H", "H", "H"],
                coords=np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [-0.97, 0.0, 0.0],
                        [3.0, 0.0, 0.0],
                        [3.74, 0.0, 0.0],
                    ]
                ),
                spin=1,
            ),
            scan_i=0,
            scan_j=2,
            scan_distances_a=(2.8, 2.5, 2.2, 1.9, 1.65, 1.45, 1.33, 1.2, 1.05),
            scan_floor_a=0.9,
            literature_barrier_k=2900.0,
            literature_imag_cm=1259.6,
            literature_fit=LiteratureFit(7.64e10, 1.0, 1339.9, 153.2, meisner, 60.0),
            method=bhlyp,
        ),
        "h-h2-exchange": Reaction(
            key="h-h2-exchange",
            label="H + H2 -> H2 + H",
            cluster=Cluster(
                name="h-h2-exchange",
                symbols=["H", "H", "H"],
                coords=np.array([[0.0, 0.0, 0.0], [0.74, 0.0, 0.0], [-3.0, 0.0, 0.0]]),
                spin=1,
            ),
            scan_i=0,
            scan_j=2,
            scan_distances_a=(2.8, 2.5, 2.2, 1.9, 1.6, 1.35, 1.15, 1.0, 0.9),
            scan_floor_a=0.75,
            literature_barrier_k=5030.0,
            literature_imag_cm=1510.0,
            literature_fit=None,
            method=bhlyp,
        ),
        "h-oh-control": Reaction(
            key="h-oh-control",
            label="H + OH -> H2O (barrierless control)",
            cluster=Cluster(
                name="h-oh-control",
                symbols=["O", "H", "H"],
                coords=np.array([[0.0, 0.0, 0.0], [-0.97, 0.0, 0.0], [3.0, 0.0, 0.0]]),
                spin=0,
            ),
            scan_i=0,
            scan_j=2,
            scan_distances_a=(2.8, 2.5, 2.2, 1.9, 1.65, 1.4, 1.2, 1.0),
            scan_floor_a=0.9,
            literature_barrier_k=None,
            literature_imag_cm=None,
            literature_fit=None,
            method=bhlyp,
            barrierless=True,
        ),
    }


def save_xyz(cluster: Cluster, path: Path) -> None:
    path.write_text(cluster.to_xyz())


def load_xyz(path: Path, template: Cluster) -> Cluster:
    lines = path.read_text().splitlines()
    n = int(lines[0])
    coords = np.array(
        [[float(v) for v in line.split()[1:4]] for line in lines[2 : 2 + n]]
    )
    return replace(template, coords=coords)


def geometry_hash(cluster: Cluster) -> str:
    return hashlib.sha256(cluster.to_xyz().encode()).hexdigest()


def frequency_payload(freq: FrequencyResult) -> dict:
    return {
        "frequencies_cm": freq.frequencies_cm.tolist(),
        "imaginary_cm": freq.imaginary_cm.tolist(),
        "electronic_hartree": freq.electronic_hartree,
        "molar_mass_kg": freq.molar_mass_kg,
        "rotational_temperatures_k": list(freq.rotational_temperatures_k or []),
        "linear": freq.linear,
    }


def thermo(freq: FrequencyResult, temperature: float):
    return thermo_from_frequencies(
        freq.electronic_hartree * HARTREE_TO_KJ,
        freq.frequencies_cm,
        temperature,
        molar_mass_kg=freq.molar_mass_kg,
        rotational_temperatures_k=(
            list(freq.rotational_temperatures_k)
            if freq.rotational_temperatures_k
            else None
        ),
        linear=freq.linear,
    )


def run_reaction(reaction: Reaction, *, force: bool = False) -> dict:
    run_dir = RUN_ROOT / reaction.key
    run_dir.mkdir(parents=True, exist_ok=True)
    result_path = run_dir / "results.json"
    if result_path.exists() and not force:
        log(f"{reaction.key}: results.json exists; resume complete")
        return json.loads(result_path.read_text())

    ts_guess_path = run_dir / "ts_guess.xyz"
    scan_path = run_dir / "scan.json"
    if ts_guess_path.exists() and not force and not reaction.barrierless:
        log(f"{reaction.key}: ts_guess.xyz exists; skipping completed scan")
        ts_guess = load_xyz(ts_guess_path, reaction.cluster)
    else:
        log(f"{reaction.key}: relaxed scan")
        try:
            scan = scan_to_maximum(
                reaction.cluster,
                reaction.method,
                atom_i=reaction.scan_i,
                atom_j=reaction.scan_j,
                distances_a=list(reaction.scan_distances_a),
                min_distance_a=reaction.scan_floor_a,
                progress=lambda r, e: log(f"{reaction.key}: r={r:.3f} A E={e:.10f} Ha"),
            )
            scan_path.write_text(
                json.dumps(
                    [{"distance_a": r, "energy_hartree": e} for r, e, _ in scan],
                    indent=2,
                )
            )
        except ScanNoMaximumError as exc:
            scan = exc.scan
            scan_path.write_text(
                json.dumps(
                    [{"distance_a": r, "energy_hartree": e} for r, e, _ in scan],
                    indent=2,
                )
            )
            if not reaction.barrierless:
                raise
            energies = np.array([e for _, e, _ in scan])
            monotonic = bool(np.all(np.diff(energies) <= 1e-5))
            result = {
                "reaction": reaction.label,
                "key": reaction.key,
                "classification": "barrierless-control",
                "scan_monotonic_downhill": monotonic,
                "barrier_zpe_kj_mol": 0.0,
                "imaginary_frequency_cm": None,
                "rates": [
                    {"temperature_k": t, "eckart_kappa": 1.0} for t in TEMPERATURES
                ],
                "method": vars(reaction.method),
                "provenance": provenance(reaction, [reaction.cluster]),
            }
            result_path.write_text(json.dumps(result, indent=2))
            return result

        if reaction.barrierless:
            raise RuntimeError("barrierless control produced an interior maximum")

        ts_guess = scan_ts_guess(scan)
        save_xyz(ts_guess, ts_guess_path)

    ts_path = run_dir / "ts.xyz"
    ts = (
        load_xyz(ts_path, ts_guess)
        if ts_path.exists()
        else find_ts(ts_guess, reaction.method, trajectory=str(run_dir / "sella.traj"))
    )
    save_xyz(ts, ts_path)
    ts_freq = frequencies(ts, reaction.method)
    if ts_freq.n_imaginary != 1 or ts_freq.imaginary_cm[0] < 200.0:
        raise RuntimeError(
            f"{reaction.key}: expected one chemical imaginary mode, "
            f"got {ts_freq.imaginary_cm}"
        )

    back_path, fwd_path = run_dir / "irc_back.xyz", run_dir / "irc_fwd.xyz"
    if back_path.exists() and fwd_path.exists():
        back, fwd = load_xyz(back_path, ts), load_xyz(fwd_path, ts)
    else:
        back, fwd = quick_irc(ts, reaction.method)
        save_xyz(back, back_path)
        save_xyz(fwd, fwd_path)

    # The reactant basin has the longer scanned bond; product has the shorter one.
    def distance(cluster: Cluster) -> float:
        return float(
            np.linalg.norm(
                cluster.coords[reaction.scan_i] - cluster.coords[reaction.scan_j]
            )
        )

    reactant, product = sorted((back, fwd), key=distance, reverse=True)
    reactant_freq, product_freq = (
        frequencies(reactant, reaction.method),
        frequencies(product, reaction.method),
    )
    if reactant_freq.n_imaginary or product_freq.n_imaginary:
        raise RuntimeError(f"{reaction.key}: quick-IRC endpoint is not a minimum")

    barrier_zpe = (
        ts_freq.electronic_hartree * HARTREE_TO_KJ
        + thermo(ts_freq, 1.0).zpe_kj
        - reactant_freq.electronic_hartree * HARTREE_TO_KJ
        - thermo(reactant_freq, 1.0).zpe_kj
    )
    reverse_zpe = (
        ts_freq.electronic_hartree * HARTREE_TO_KJ
        + thermo(ts_freq, 1.0).zpe_kj
        - product_freq.electronic_hartree * HARTREE_TO_KJ
        - thermo(product_freq, 1.0).zpe_kj
    )
    if barrier_zpe <= 0 or reverse_zpe <= 0:
        detail = f"forward={barrier_zpe}, reverse={reverse_zpe}"
        raise RuntimeError(f"{reaction.key}: non-positive Eckart barrier ({detail})")

    rates = []
    imag = float(ts_freq.imaginary_cm[0])
    for temperature in TEMPERATURES:
        tr, tt = thermo(reactant_freq, temperature), thermo(ts_freq, temperature)
        dg = tt.gibbs - tr.gibbs
        kappa = eckart_kappa(imag, barrier_zpe, reverse_zpe, temperature)
        classical = (
            1.380649e-23
            * temperature
            / 6.62607015e-34
            * np.exp(-dg / (0.00831446261815324 * temperature))
        )
        lit = (
            reaction.literature_fit.rate(temperature)
            if reaction.literature_fit
            else None
        )
        rates.append(
            {
                "temperature_k": temperature,
                "delta_g_kj_mol": dg,
                "eckart_kappa": kappa,
                "quarry_k_s-1": kappa * classical,
                "literature_k_s-1": lit,
                "quarry_over_literature": (kappa * classical / lit if lit else None),
                "literature_extrapolated": bool(
                    reaction.literature_fit
                    and temperature < reaction.literature_fit.valid_floor_k
                ),
            }
        )

    result = {
        "reaction": reaction.label,
        "key": reaction.key,
        "classification": "first-order-saddle",
        "barrier_zpe_kj_mol": barrier_zpe,
        "reverse_barrier_zpe_kj_mol": reverse_zpe,
        "imaginary_frequency_cm": imag,
        "literature_barrier_k": reaction.literature_barrier_k,
        "literature_imaginary_cm": reaction.literature_imag_cm,
        "rates": rates,
        "method": vars(reaction.method),
        "frequency_data": {
            "reactant": frequency_payload(reactant_freq),
            "ts": frequency_payload(ts_freq),
            "product": frequency_payload(product_freq),
        },
        "provenance": provenance(reaction, [reactant, ts, product]),
    }
    result_path.write_text(json.dumps(result, indent=2))
    return result


def provenance(reaction: Reaction, structures: list[Cluster]) -> dict:
    try:
        git_sha = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parent, text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        git_sha = "unknown"
    return {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "git_sha": git_sha,
        "python": platform.python_version(),
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
        "gpu_requested": reaction.method.use_gpu,
        "literature_source": reaction.literature_fit.source
        if reaction.literature_fit
        else None,
        "geometry_sha256": {
            cluster.name: geometry_hash(cluster) for cluster in structures
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reaction", action="append", help="reaction key (repeatable; default all)"
    )
    parser.add_argument(
        "--gpu", action="store_true", help="use GPU4PySCF (required workstation path)"
    )
    parser.add_argument(
        "--basis", default="def2-svp", help="DFT basis for calibration campaign"
    )
    parser.add_argument(
        "--force", action="store_true", help="recompute completed results"
    )
    args = parser.parse_args()
    os.environ.setdefault("OMP_NUM_THREADS", "16")
    if int(os.environ["OMP_NUM_THREADS"]) > 16:
        raise SystemExit("OMP_NUM_THREADS must be <=16")
    all_reactions = reactions(gpu=args.gpu, basis=args.basis)
    selected = args.reaction or list(all_reactions)
    unknown = sorted(set(selected) - set(all_reactions))
    if unknown:
        raise SystemExit(f"unknown reactions: {unknown}")
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    results = []
    failures = []
    for key in selected:
        try:
            results.append(run_reaction(all_reactions[key], force=args.force))
        except Exception as exc:  # campaign isolation: preserve and continue
            failure = {
                "key": key,
                "error_type": type(exc).__name__,
                "error": str(exc),
                "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            }
            failures.append(failure)
            failure_path = RUN_ROOT / key / "failure.json"
            failure_path.write_text(json.dumps(failure, indent=2))
            log(f"{key}: FAILED but campaign continues: {exc}")
    (RUN_ROOT / "results.json").write_text(json.dumps(results, indent=2))
    (RUN_ROOT / "failures.json").write_text(json.dumps(failures, indent=2))
    log(
        f"complete: {len(results)} results, {len(failures)} failures -> "
        f"{RUN_ROOT / 'results.json'}"
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
