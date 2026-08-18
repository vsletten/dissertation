#!/usr/bin/env python
"""Phase 1: reproduce the Xiao & Lasaga hydrolysis barriers.

Runs one reaction end-to-end (HANDOFF.md §4, SURVEY.md §7.1):

    uv run python scripts/phase1_xiao_lasaga.py --reaction si-neutral
    uv run python scripts/phase1_xiao_lasaga.py --reaction si-acid --gpu

Reactions: si-neutral / si-acid attack the (HO)3Si-O-Si(OH)3 dimer with
H2O / H3O+; al-neutral / al-acid do the same to (HO)3Si-O-Al(OH)3^-.
Literature anchors: X&L 1994/96 got ~29 kcal/mol (neutral), ~24 (acid,
H3O+), ~19 (base); modern AIMD puts Q1-ish Si-O-T hydrolysis at 54-71
kJ/mol (SURVEY.md §7.2).

Protocol (each stage checkpointed as XYZ under runs/phase1/<reaction>/,
so a crashed campaign resumes where it stopped):
  1. optimize the reactant complex and the separated fragments
  2. relaxed scan of r(Si-Ow) toward the pentacoordinate geometry
  3. Sella saddle search from the scan maximum
  4. verify (exactly 1 imaginary mode) + quick-IRC basin check
  5. frequencies on complex and TS -> quasi-RRHO dG‡(298.15 K)
  6. barriers vs both references (complex and separated fragments),
     results.json + SQLite store with provenance

Every number this prints is also written to disk. Nothing is lost if a
terminal dies.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import numpy as np  # noqa: E402

from quarry.clusters import (  # noqa: E402
    Cluster,
    aluminosilicate_dimer,
    disilicate,
    hydrolysis_complex,
    hydronium,
    water,
)
from quarry.pipeline import (  # noqa: E402
    HARTREE_TO_KJ,
    DftSettings,
    energy,
    frequencies,
    optimize,
)
from quarry.rates import rate_from_thermo, thermo_from_frequencies  # noqa: E402
from quarry.store import Store  # noqa: E402
from quarry.ts import find_ts, quick_irc, scan_maximum, scan_to_maximum  # noqa: E402

REACTIONS = {
    "si-neutral": (disilicate, water),
    "si-acid": (disilicate, hydronium),
    "al-neutral": (aluminosilicate_dimer, water),
    "al-acid": (aluminosilicate_dimer, hydronium),
}
# Indices from the builders: 0 = bridging O, 1 = Si under attack;
# the attacker's O is the first atom appended after the dimer.
SI_INDEX = 1
KCAL = 4.184


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


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


def checkpointed(path: Path, template: Cluster, compute) -> Cluster:
    """Load a stage result from disk or compute-and-save it."""
    if path.exists():
        log(f"  resume: {path.name} exists, skipping recompute")
        return load_xyz(path, template)
    result = compute()
    save_xyz(result, path)
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reaction", choices=sorted(REACTIONS), required=True)
    ap.add_argument("--xc", default="b3lyp")
    ap.add_argument(
        "--basis",
        default="def2-svp",
        help="geometry/frequency basis (def2-svp first pass)",
    )
    ap.add_argument("--gpu", action="store_true")
    ap.add_argument(
        "--threads",
        type=int,
        default=16,
        help="OMP thread cap for CPU stages (default 16 — leave the box usable)",
    )
    ap.add_argument("--temperature", type=float, default=298.15)
    args = ap.parse_args()

    os.environ.setdefault("OMP_NUM_THREADS", str(args.threads))
    settings = DftSettings(
        xc=args.xc, basis=args.basis, density_fit=True, use_gpu=args.gpu
    )

    run_dir = (
        Path(__file__).resolve().parent.parent
        / "runs"
        / "phase1"
        / (f"{args.reaction}-{args.xc}-{args.basis}")
    )
    run_dir.mkdir(parents=True, exist_ok=True)
    log(f"run dir: {run_dir}")
    log(f"settings: {settings}")

    dimer_factory, attacker_factory = REACTIONS[args.reaction]
    dimer, attacker = dimer_factory(), attacker_factory()
    complex_guess = hydrolysis_complex(dimer, attacker)
    ow_index = len(dimer.symbols)  # attacker O

    # Stage 0 — cheap, robust pre-optimization of the hand-built guess.
    # HF/STO-3G converges where a hybrid at a strained guess may not;
    # the production method then starts from a relaxed structure.
    log("stage 0: HF/STO-3G pre-optimization of the complex guess")
    preopt_settings = DftSettings(xc="hf", basis="sto-3g")
    complex_pre = checkpointed(
        run_dir / "complex_preopt.xyz",
        complex_guess,
        lambda: optimize(complex_guess, preopt_settings),
    )

    # Stage 1 — optimize reactant complex and separated fragments.
    log("stage 1: optimizing reactant complex + fragments")
    complex_opt = checkpointed(
        run_dir / "complex.xyz",
        complex_pre,
        lambda: optimize(complex_pre, settings),
    )
    dimer_opt = checkpointed(
        run_dir / "dimer.xyz", dimer, lambda: optimize(dimer, settings)
    )
    attacker_opt = checkpointed(
        run_dir / "attacker.xyz",
        attacker,
        lambda: optimize(attacker, settings),
    )

    # Stage 2 — relaxed scan pulling the attacker O onto Si, extended
    # until the energy maximum is interior (ridge actually crossed).
    log("stage 2: relaxed scan r(Si-Ow), auto-extending to interior maximum")
    ts_guess_path = run_dir / "ts_guess.xyz"
    if ts_guess_path.exists():
        ts_guess = load_xyz(ts_guess_path, complex_opt)
        log("  resume: ts_guess.xyz exists")
    else:
        scan = scan_to_maximum(
            complex_opt,
            settings,
            atom_i=SI_INDEX,
            atom_j=ow_index,
            distances_a=[2.8, 2.6, 2.4, 2.2, 2.1, 2.0, 1.9],
            progress=lambda r, e: log(f"  r={r:.2f} A  E={e:.6f} Ha"),
        )
        ts_guess = scan_maximum(scan)
        save_xyz(ts_guess, ts_guess_path)

    # Stage 3 — saddle search.
    log("stage 3: Sella saddle search")
    ts = checkpointed(
        run_dir / "ts.xyz",
        ts_guess,
        lambda: find_ts(ts_guess, settings, trajectory=str(run_dir / "sella.traj")),
    )

    # Guard: a saddle that drifted back out of the approach channel is a
    # trivial complex-rearrangement saddle, not the hydrolysis TS
    # (observed live: 208i water wag at r(Si-Ow)=3.25 from a bad guess).
    r_guess = float(
        np.linalg.norm(ts_guess.coords[SI_INDEX] - ts_guess.coords[ow_index])
    )
    r_ts = float(np.linalg.norm(ts.coords[SI_INDEX] - ts.coords[ow_index]))
    log(f"  r(Si-Ow): guess {r_guess:.2f} A -> saddle {r_ts:.2f} A")
    if r_ts > r_guess + 0.5:
        log(
            "  !! saddle escaped the approach channel (r(Si-Ow) grew by "
            f"{r_ts - r_guess:.2f} A) — this is not the hydrolysis TS. "
            "Delete ts.xyz/ts_guess.xyz and rerun with a tighter scan, or "
            "drive a combined coordinate. Aborting before thermochemistry."
        )
        return 1

    # Stage 4 — verify + quick IRC.
    log("stage 4: frequencies + quick-IRC")
    ts_freq = frequencies(ts, settings)
    log(f"  TS imaginary modes: {np.round(ts_freq.imaginary_cm, 1).tolist()}")
    if ts_freq.n_imaginary != 1:
        log("  !! not a clean first-order saddle — inspect ts.xyz; aborting")
        return 1
    back, fwd = quick_irc(ts, settings)
    save_xyz(back, run_dir / "irc_back.xyz")
    save_xyz(fwd, run_dir / "irc_fwd.xyz")

    # Stage 5 — thermochemistry.
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

    th_cx = thermo(cx_freq)
    th_ts = thermo(ts_freq)
    rate = rate_from_thermo(
        th_cx,
        th_ts,
        imag_nu_cm=float(ts_freq.imaginary_cm[0]),
        tunneling="wigner",
    )

    de_kj = (ts_freq.electronic_hartree - cx_freq.electronic_hartree) * HARTREE_TO_KJ
    # Electronic barrier vs separated, relaxed fragments (X&L's reference).
    frag_e = (
        ts_freq.electronic_hartree
        - energy(dimer_opt, settings)
        - energy(attacker_opt, settings)
    ) * HARTREE_TO_KJ
    results = {
        "reaction": args.reaction,
        "method": f"{args.xc}/{args.basis}/df",
        "temperature_k": t,
        "dE_elec_vs_complex_kj": de_kj,
        "dH_kj": rate.dh_kj,
        "dG_kj": rate.dg_kj,
        "dG_kcal": rate.dg_kj / KCAL,
        "wigner_kappa": rate.kappa,
        "k_per_s": rate.k,
        "ts_imaginary_cm": float(ts_freq.imaginary_cm[0]),
        "literature": {
            "xiao_lasaga_acid_kcal": 24.0,
            "xiao_lasaga_base_kcal": 19.0,
            "xiao_lasaga_neutral_kcal": 29.0,
            "modern_aimd_q1_q3_kj": [54.0, 81.0],
        },
    }
    results["dE_elec_vs_fragments_kj"] = frag_e

    (run_dir / "results.json").write_text(json.dumps(results, indent=2))

    with Store(run_dir / "store.sqlite") as store:
        for name, cl, freq in (
            ("complex", complex_opt, cx_freq),
            ("ts", ts, ts_freq),
        ):
            sid = store.add_structure(
                f"{args.reaction}-{name}",
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
        sid = store.add_structure(
            f"{args.reaction}-attacker",
            attacker_opt.formula,
            attacker_opt.to_xyz(),
            charge=attacker_opt.charge,
            spin=attacker_opt.spin,
        )

    log("")
    log(f"=== {args.reaction} @ {results['method']} ===")
    log(f"  dE(elec, TS - complex) = {de_kj:8.1f} kJ/mol")
    log(f"  dH‡(298)               = {rate.dh_kj:8.1f} kJ/mol")
    log(
        f"  dG‡(298)               = {rate.dg_kj:8.1f} kJ/mol "
        f"({rate.dg_kj / KCAL:.1f} kcal/mol)"
    )
    log(f"  wigner kappa           = {rate.kappa:8.2f}")
    log(f"  k(298)                 = {rate.k:8.3e} 1/s")
    log(f"  TS imag mode           = {ts_freq.imaginary_cm[0]:8.1f}i cm^-1")
    log(f"  dE(elec, TS - fragments) = {frag_e:6.1f} kJ/mol")
    log("  X&L anchors: neutral ~29, acid ~24, base ~19 kcal/mol")
    log(f"results.json + store.sqlite written to {run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
