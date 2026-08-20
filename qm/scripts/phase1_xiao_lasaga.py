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
    FrequencyResult,
    energy,
    frequencies,
    optimize,
)
from quarry.rates import Thermo, rate_from_thermo, thermo_from_frequencies  # noqa: E402
from quarry.store import Store  # noqa: E402
from quarry.ts import (  # noqa: E402
    ScanNoMaximumError,
    constrained_scan,
    find_ts,
    first_interior_maximum,
    neb_ts_guess,
    quick_irc,
    reaction_path_vector,
    scan_to_maximum,
    scan_ts_guess,
)

REACTIONS = {
    "si-neutral": (disilicate, water),
    "si-acid": (disilicate, hydronium),
    "al-neutral": (aluminosilicate_dimer, water),
    "al-acid": (aluminosilicate_dimer, hydronium),
}
# Indices from the builders: 0 = bridging O, 1 = Si under attack;
# the attacker's O is the first atom appended after the dimer.
BR_INDEX = 0
SI_INDEX = 1
KCAL = 4.184
SI_NEUTRAL_DG_KJ = 113.048262
SI_NEUTRAL_DG_KCAL = 27.019
REACTANT_BASIN = (False, True, False)
ASSOCIATIVE_BASIN = (True, True, True)
HYDROLYZED_BASIN = (True, False, True)


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


def channel_escape_reason(ts_guess: Cluster, ts: Cluster, ow_index: int) -> str | None:
    """Explain why a refined saddle left the hydrolysis channel, if it did."""
    r_guess = float(
        np.linalg.norm(ts_guess.coords[SI_INDEX] - ts_guess.coords[ow_index])
    )
    r_ts = float(np.linalg.norm(ts.coords[SI_INDEX] - ts.coords[ow_index]))
    if r_ts > 2.6:
        return f"saddle r(Si-Ow)={r_ts:.2f} A is outside the bonding channel"
    if r_ts > r_guess + 0.5:
        return f"saddle escaped the approach channel by {r_ts - r_guess:.2f} A"
    return None


def hydrolysis_basin_signature(
    cluster: Cluster, ow_index: int
) -> tuple[bool, bool, bool]:
    """Return (Si-Ow bonded, Si-Obr bonded, proton on Obr) for IRC gating."""
    h_candidates = [ow_index + 1, ow_index + 2]
    return (
        float(np.linalg.norm(cluster.coords[SI_INDEX] - cluster.coords[ow_index]))
        < 2.3,
        float(np.linalg.norm(cluster.coords[SI_INDEX] - cluster.coords[BR_INDEX]))
        < 2.3,
        min(
            float(np.linalg.norm(cluster.coords[BR_INDEX] - cluster.coords[h_idx]))
            for h_idx in h_candidates
        )
        < 1.25,
    )


def quick_irc_channel_reason(back: Cluster, fwd: Cluster, ow_index: int) -> str | None:
    """Explain why quick-IRC did not connect reactant and hydrolyzed basins."""
    return quick_irc_basin_reason(
        back,
        fwd,
        ow_index,
        expected={REACTANT_BASIN, HYDROLYZED_BASIN},
    )


def quick_irc_basin_reason(
    back: Cluster,
    fwd: Cluster,
    ow_index: int,
    *,
    expected: set[tuple[bool, bool, bool]],
) -> str | None:
    """Explain why two quick-IRC endpoints missed the expected basins."""
    actual = {
        hydrolysis_basin_signature(back, ow_index),
        hydrolysis_basin_signature(fwd, ow_index),
    }
    if actual != expected:
        return f"quick-IRC basin signatures {sorted(actual)} != {sorted(expected)}"
    return None


def quick_irc_addition_reason(back: Cluster, fwd: Cluster, ow_index: int) -> str | None:
    """Require reactant ↔ associative-intermediate connectivity."""
    return quick_irc_basin_reason(
        back,
        fwd,
        ow_index,
        expected={REACTANT_BASIN, ASSOCIATIVE_BASIN},
    )


def quick_irc_cleavage_reason(back: Cluster, fwd: Cluster, ow_index: int) -> str | None:
    """Require associative-intermediate ↔ hydrolyzed-product connectivity."""
    return quick_irc_basin_reason(
        back,
        fwd,
        ow_index,
        expected={ASSOCIATIVE_BASIN, HYDROLYZED_BASIN},
    )


def endpoint_in_basin(
    back: Cluster,
    fwd: Cluster,
    ow_index: int,
    signature: tuple[bool, bool, bool],
) -> Cluster:
    """Return the unique quick-IRC endpoint in ``signature`` or fail closed."""
    matches = [
        endpoint
        for endpoint in (back, fwd)
        if hydrolysis_basin_signature(endpoint, ow_index) == signature
    ]
    if len(matches) != 1:
        detail = f"basin {signature}, found {len(matches)}"
        raise RuntimeError(f"expected one quick-IRC endpoint in {detail}")
    return matches[0]


def sequential_barrier_metrics(
    reactant: Thermo,
    intermediate: Thermo,
    addition_ts: Thermo,
    cleavage_ts: Thermo,
) -> dict[str, float | str]:
    """Free-energy profile metrics for R → I → P through two saddles.

    Local forward barriers identify the slow elementary step.  The overall
    profile barrier, used for the like-for-like comparison with a concerted
    reactant-referenced barrier, is the highest saddle above the initial
    reactant complex.  Both are reported because a stabilized intermediate can
    make those two notions differ substantially.
    """
    addition_local = addition_ts.gibbs - reactant.gibbs
    cleavage_local = cleavage_ts.gibbs - intermediate.gibbs
    addition_from_reactant = addition_local
    cleavage_from_reactant = cleavage_ts.gibbs - reactant.gibbs
    local = {"addition": addition_local, "cleavage": cleavage_local}
    profile = {
        "addition": addition_from_reactant,
        "cleavage": cleavage_from_reactant,
    }
    return {
        "intermediate_dG_vs_reactant_kj": intermediate.gibbs - reactant.gibbs,
        "addition_dG_dagger_kj": addition_local,
        "cleavage_dG_dagger_kj": cleavage_local,
        "addition_ts_dG_vs_reactant_kj": addition_from_reactant,
        "cleavage_ts_dG_vs_reactant_kj": cleavage_from_reactant,
        "overall_profile_dG_dagger_kj": max(profile.values()),
        "rate_limiting_local_step": max(local, key=local.__getitem__),
        "highest_profile_ts": max(profile, key=profile.__getitem__),
    }


def thermo_result(freq: FrequencyResult, temperature: float) -> Thermo:
    """Convert one stationary point's projected frequencies to quasi-RRHO."""
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


def proton_neb_guess(
    approach_seed: Cluster,
    complex_opt: Cluster,
    settings: DftSettings,
    run_dir: Path,
    ow_index: int,
) -> Cluster:
    """Build the concerted proton-transfer/product/CI-NEB TS guess.

    ``approach_seed`` may be the r=1.90 A point from a completed scan or a
    resumed direct-crest checkpoint. Existing product checkpoints are
    validated before reuse; a rolled-back product is extended farther along
    Si-Obr cleavage rather than repeating the expensive proton scan.
    """
    product_path = run_dir / "product.xyz"
    product_seed: Cluster | None = None
    h_candidates = [ow_index + 1, ow_index + 2]

    if product_path.exists():
        product = load_xyz(product_path, approach_seed)
        r_prod_ow = float(
            np.linalg.norm(product.coords[SI_INDEX] - product.coords[ow_index])
        )
        r_prod_br = float(
            np.linalg.norm(product.coords[SI_INDEX] - product.coords[BR_INDEX])
        )
        log(
            "  resume: product.xyz exists; "
            f"r(Si-Ow)={r_prod_ow:.2f} A, r(Si-Obr)={r_prod_br:.2f} A"
        )
        if r_prod_ow <= 1.9 and r_prod_br >= 2.2:
            # Endpoint alignment in ``neb_ts_guess`` removes the rigid-frame
            # mismatch that made the full hydrolysis product fail its first
            # image SCF.  Use the actual product basin here: the older
            # rolled-back r(Si-Obr)≈2.07 checkpoint yielded a product-like
            # four-imaginary-mode NEB peak, not the concerted saddle.
            log(
                "  stage 2d: staged CI-NEB complex -> aligned full product "
                "(5 images; MDMin 1.5 pre-relax -> 0.2 climb eV/A)"
            )
            return neb_ts_guess(
                complex_opt,
                product,
                settings,
                n_images=5,
                fmax_ev_a=0.2,
                max_steps=240,
                pre_relax_fmax_ev_a=1.5,
                pre_relax_steps=120,
            )
        log("  saved product is not hydrolyzed; extending Si-Obr cleavage")
        rejected = run_dir / "product.rejected-rollback.xyz"
        product_path.replace(rejected)
        product_seed = product

    if product_seed is None:
        seed_r = float(
            np.linalg.norm(
                approach_seed.coords[SI_INDEX] - approach_seed.coords[ow_index]
            )
        )
        if abs(seed_r - 1.90) > 0.02:
            log(f"  pinning resumed approach seed: r(Si-Ow) {seed_r:.2f} -> 1.90 A")
            approach_seed = constrained_scan(
                approach_seed,
                settings,
                atom_i=SI_INDEX,
                atom_j=ow_index,
                distances_a=[1.90],
            )[0][2]

        log("  stage 2b: pin r(Si-Ow)=1.90 A, drive H(water) -> O(bridge)")
        h_idx = min(
            h_candidates,
            key=lambda i: np.linalg.norm(
                approach_seed.coords[i] - approach_seed.coords[BR_INDEX]
            ),
        )
        r0 = float(
            np.linalg.norm(approach_seed.coords[h_idx] - approach_seed.coords[BR_INDEX])
        )
        log(f"  driving H{h_idx} from r(Obr-H)={r0:.2f} A inward")
        first = max(1.05, round(r0 - 0.15, 2))
        distances = [round(float(x), 2) for x in np.arange(first, 1.04, -0.15)]
        pscan = scan_to_maximum(
            approach_seed,
            settings,
            atom_i=BR_INDEX,
            atom_j=h_idx,
            distances_a=distances or [first],
            fixed_distances=[(SI_INDEX, ow_index, 1.90)],
            extend_step_a=0.06,
            min_distance_a=0.95,
            progress=lambda r, e: log(f"  r(Obr-H)={r:.2f} A  E={e:.6f} Ha"),
        )
        crest_idx = first_interior_maximum(pscan)
        if crest_idx is None:
            raise RuntimeError("proton scan produced no interior crest")
        product_seed = pscan[min(crest_idx + 1, len(pscan) - 1)][2]

    h_idx = min(
        h_candidates,
        key=lambda i: np.linalg.norm(
            product_seed.coords[i] - product_seed.coords[BR_INDEX]
        ),
    )
    current_br = float(
        np.linalg.norm(product_seed.coords[SI_INDEX] - product_seed.coords[BR_INDEX])
    )
    break_targets = [
        r for r in [1.85, 2.05, 2.30, 2.60, 3.00, 3.40, 3.60] if r > current_br + 0.05
    ]
    if not break_targets:
        break_targets = [round(current_br + 0.30, 2)]

    log("  stage 2c: breaking Si-Obr with the new bonds held")
    broken = constrained_scan(
        product_seed,
        settings,
        atom_i=SI_INDEX,
        atom_j=BR_INDEX,
        distances_a=break_targets,
        fixed_distances=[
            (SI_INDEX, ow_index, 1.70),
            (BR_INDEX, h_idx, 0.98),
        ],
    )
    for r, e, _ in broken:
        log(f"  r(Si-Obr)={r:.2f} A  E={e:.6f} Ha")
    product = optimize(broken[-1][2], settings)
    save_xyz(product, product_path)

    r_prod_ow = float(
        np.linalg.norm(product.coords[SI_INDEX] - product.coords[ow_index])
    )
    r_prod_br = float(
        np.linalg.norm(product.coords[SI_INDEX] - product.coords[BR_INDEX])
    )
    log(f"  product r(Si-Ow)={r_prod_ow:.2f} A, r(Si-Obr)={r_prod_br:.2f} A")
    if r_prod_ow > 1.9 or r_prod_br < 2.2:
        raise RuntimeError(
            "product rolled back toward reactants; expected Si-Ow bonded "
            "and Si-Obr broken"
        )
    log(
        "  stage 2d: staged CI-NEB complex -> aligned full product "
        "(5 images; MDMin 1.5 pre-relax -> 0.2 climb eV/A)"
    )
    return neb_ts_guess(
        complex_opt,
        product,
        settings,
        n_images=5,
        fmax_ev_a=0.2,
        max_steps=240,
        pre_relax_fmax_ev_a=1.5,
        pre_relax_steps=120,
    )


def finish_al_neutral_sequential(
    *,
    complex_opt: Cluster,
    cleavage_ts: Cluster,
    cleavage_freq: FrequencyResult,
    cleavage_back: Cluster,
    cleavage_fwd: Cluster,
    settings: DftSettings,
    run_dir: Path,
    ow_index: int,
    temperature: float,
    dimer_opt: Cluster,
    attacker_opt: Cluster,
) -> int:
    """Close the proven R → associative I → hydrolyzed P mechanism."""
    cleavage_failure = quick_irc_cleavage_reason(cleavage_back, cleavage_fwd, ow_index)
    if cleavage_failure:
        raise RuntimeError(f"invalid cleavage handoff: {cleavage_failure}")

    intermediate = endpoint_in_basin(
        cleavage_back, cleavage_fwd, ow_index, ASSOCIATIVE_BASIN
    )
    hydrolyzed = endpoint_in_basin(
        cleavage_back, cleavage_fwd, ow_index, HYDROLYZED_BASIN
    )
    save_xyz(intermediate, run_dir / "intermediate.xyz")
    save_xyz(hydrolyzed, run_dir / "hydrolyzed_product.xyz")
    save_xyz(cleavage_ts, run_dir / "cleavage_ts.xyz")
    save_xyz(cleavage_back, run_dir / "cleavage_irc_back.xyz")
    save_xyz(cleavage_fwd, run_dir / "cleavage_irc_fwd.xyz")

    log("stage 4b: directed reactant -> associative-intermediate saddle")
    addition_guess = checkpointed(
        run_dir / "addition_ts_guess.xyz",
        complex_opt,
        lambda: neb_ts_guess(
            complex_opt,
            intermediate,
            settings,
            n_images=5,
            fmax_ev_a=0.20,
            max_steps=240,
            pre_relax_fmax_ev_a=1.0,
            pre_relax_steps=120,
        ),
    )
    addition_ts = checkpointed(
        run_dir / "addition_directed_ts.xyz",
        addition_guess,
        lambda: find_ts(
            addition_guess,
            settings,
            trajectory=str(run_dir / "addition_directed_sella.traj"),
            initial_mode=reaction_path_vector(complex_opt, intermediate),
            internal=False,
        ),
    )
    addition_freq = frequencies(addition_ts, settings)
    log(
        "  addition TS imaginary modes: "
        f"{np.round(addition_freq.imaginary_cm, 1).tolist()}"
    )
    if addition_freq.n_imaginary != 1:
        log("  !! addition saddle is not first order; aborting")
        return 1

    accepted_addition_irc = None
    for displacement in (0.50, 1.00):
        back, fwd = quick_irc(
            addition_ts,
            settings,
            displacement_a=displacement,
            frequency=addition_freq,
        )
        tag = f"{displacement:.2f}".replace(".", "p")
        save_xyz(back, run_dir / f"addition_irc_back_{tag}.xyz")
        save_xyz(fwd, run_dir / f"addition_irc_fwd_{tag}.xyz")
        failure = quick_irc_addition_reason(back, fwd, ow_index)
        if failure is None:
            accepted_addition_irc = (back, fwd, displacement)
            break
        log(f"  displacement {displacement:.2f} A: {failure}")
    if accepted_addition_irc is None:
        log("  !! addition saddle does not connect reactant and intermediate")
        return 1
    addition_back, addition_fwd, accepted_displacement = accepted_addition_irc
    save_xyz(addition_back, run_dir / "addition_irc_back.xyz")
    save_xyz(addition_fwd, run_dir / "addition_irc_fwd.xyz")

    log("stage 5: sequential thermochemistry")
    complex_freq = frequencies(complex_opt, settings)
    intermediate_freq = frequencies(intermediate, settings)
    for label, freq in (
        ("reactant", complex_freq),
        ("intermediate", intermediate_freq),
    ):
        if freq.n_imaginary:
            log(f"  !! {label} minimum has {freq.n_imaginary} imaginary modes")
            return 1

    th_complex = thermo_result(complex_freq, temperature)
    th_intermediate = thermo_result(intermediate_freq, temperature)
    th_addition = thermo_result(addition_freq, temperature)
    th_cleavage = thermo_result(cleavage_freq, temperature)
    addition_rate = rate_from_thermo(
        th_complex,
        th_addition,
        imag_nu_cm=float(addition_freq.imaginary_cm[0]),
        tunneling="wigner",
    )
    cleavage_rate = rate_from_thermo(
        th_intermediate,
        th_cleavage,
        imag_nu_cm=float(cleavage_freq.imaginary_cm[0]),
        tunneling="wigner",
    )
    metrics = sequential_barrier_metrics(
        th_complex, th_intermediate, th_addition, th_cleavage
    )
    highest = str(metrics["highest_profile_ts"])
    highest_thermo = th_addition if highest == "addition" else th_cleavage
    highest_freq = addition_freq if highest == "addition" else cleavage_freq
    profile_rate = rate_from_thermo(
        th_complex,
        highest_thermo,
        imag_nu_cm=float(highest_freq.imaginary_cm[0]),
        tunneling="wigner",
    )
    overall_dg = float(metrics["overall_profile_dG_dagger_kj"])
    si_neutral_dg = SI_NEUTRAL_DG_KJ
    al_easier = overall_dg < si_neutral_dg
    de_kj = (
        highest_freq.electronic_hartree - complex_freq.electronic_hartree
    ) * HARTREE_TO_KJ
    frag_e = (
        highest_freq.electronic_hartree
        - energy(dimer_opt, settings)
        - energy(attacker_opt, settings)
    ) * HARTREE_TO_KJ

    results = {
        "reaction": "al-neutral",
        "mechanism": "sequential-associative",
        "method": f"{settings.xc}/{settings.basis}/df",
        "temperature_k": temperature,
        "dE_elec_vs_complex_kj": de_kj,
        "dH_kj": profile_rate.dh_kj,
        "dG_kj": overall_dg,
        "dG_kcal": overall_dg / KCAL,
        "wigner_kappa": profile_rate.kappa,
        "k_per_s": profile_rate.k,
        "ts_imaginary_cm": float(highest_freq.imaginary_cm[0]),
        "dE_elec_vs_fragments_kj": frag_e,
        "profile": metrics,
        "steps": {
            "addition": {
                "dH_kj": addition_rate.dh_kj,
                "dG_kj": addition_rate.dg_kj,
                "dG_kcal": addition_rate.dg_kj / KCAL,
                "wigner_kappa": addition_rate.kappa,
                "k_per_s": addition_rate.k,
                "ts_imaginary_cm": float(addition_freq.imaginary_cm[0]),
                "accepted_quick_irc_displacement_a": accepted_displacement,
            },
            "cleavage": {
                "dH_kj": cleavage_rate.dh_kj,
                "dG_kj": cleavage_rate.dg_kj,
                "dG_kcal": cleavage_rate.dg_kj / KCAL,
                "wigner_kappa": cleavage_rate.kappa,
                "k_per_s": cleavage_rate.k,
                "ts_imaginary_cm": float(cleavage_freq.imaginary_cm[0]),
                "accepted_quick_irc_displacement_a": 0.50,
            },
        },
        "comparison": {
            "si_neutral_dG_kj": si_neutral_dg,
            "si_neutral_dG_kcal": SI_NEUTRAL_DG_KCAL,
            "al_easier_than_si": al_easier,
            "xiao_lasaga_easier_si_o_al_ordering_confirmed": al_easier,
        },
        "literature": {
            "xiao_lasaga_acid_kcal": 24.0,
            "xiao_lasaga_base_kcal": 19.0,
            "xiao_lasaga_neutral_kcal": 29.0,
            "modern_aimd_q1_q3_kj": [54.0, 81.0],
        },
    }
    results_tmp = run_dir / "results.json.tmp"
    results_tmp.write_text(json.dumps(results, indent=2))

    store_tmp = run_dir / "store.sequential.tmp.sqlite"
    store_tmp.unlink(missing_ok=True)
    with Store(store_tmp) as store:
        for name, cluster, freq in (
            ("complex", complex_opt, complex_freq),
            ("intermediate", intermediate, intermediate_freq),
            ("addition-ts", addition_ts, addition_freq),
            ("cleavage-ts", cleavage_ts, cleavage_freq),
        ):
            sid = store.add_structure(
                f"al-neutral-{name}",
                cluster.formula,
                cluster.to_xyz(),
                charge=cluster.charge,
                spin=cluster.spin,
            )
            jid = store.add_job(
                sid,
                "freq",
                str(results["method"]),
                "gpu4pyscf" if settings.use_gpu else "pyscf",
            )
            store.set_job_status(jid, "done")
            store.add_result(jid, "electronic", freq.electronic_hartree, "hartree")
    store_tmp.replace(run_dir / "store.sqlite")
    results_tmp.replace(run_dir / "results.json")

    log("")
    log(f"=== al-neutral sequential @ {results['method']} ===")
    log(
        f"  addition dG‡(298) = {addition_rate.dg_kj:8.3f} kJ/mol "
        f"({addition_rate.dg_kj / KCAL:.3f} kcal/mol)"
    )
    log(
        f"  cleavage dG‡(298) = {cleavage_rate.dg_kj:8.3f} kJ/mol "
        f"({cleavage_rate.dg_kj / KCAL:.3f} kcal/mol)"
    )
    log(
        f"  overall profile dG‡(298) = {overall_dg:8.3f} kJ/mol "
        f"({overall_dg / KCAL:.3f} kcal/mol)"
    )
    ordering = "CONFIRMED" if al_easier else "NOT CONFIRMED"
    log(f"  X&L easier Si-O-Al ordering: {ordering} vs si-neutral {si_neutral_dg:.3f}")
    log(f"results.json + store.sqlite written to {run_dir}")
    return 0


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
    ap.add_argument(
        "--approach",
        choices=["flank", "backside"],
        default="flank",
        help="attack geometry: flank = X&L 4-center (proton can reach the "
        "bridging O), backside = SN2-like (no concerted transfer)",
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
        / (f"{args.reaction}-{args.xc}-{args.basis}-{args.approach}")
    )
    run_dir.mkdir(parents=True, exist_ok=True)
    log(f"run dir: {run_dir}")
    log(f"settings: {settings}")

    dimer_factory, attacker_factory = REACTIONS[args.reaction]
    dimer, attacker = dimer_factory(), attacker_factory()
    complex_guess = hydrolysis_complex(dimer, attacker, mode=args.approach)
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
    # The neutral 4-center mechanism needs a second driven coordinate:
    # if approach alone never peaks (observed live: monotonic to 1.66 A,
    # +115 kJ/mol), pin r(Si-Ow) and drive the water proton onto the
    # bridging O — the concerted proton transfer of the X&L TS.
    log("stage 2: relaxed scan r(Si-Ow), auto-extending to interior maximum")
    ts_guess_path = run_dir / "ts_guess.xyz"
    route_path = run_dir / "ts_guess.route"
    route = route_path.read_text().strip() if route_path.exists() else "direct"
    if ts_guess_path.exists():
        ts_guess = load_xyz(ts_guess_path, complex_opt)
        log(f"  resume: ts_guess.xyz exists ({route} route)")
        if route == "proton-neb" and channel_escape_reason(
            ts_guess, ts_guess, ow_index
        ):
            r_bad = float(
                np.linalg.norm(ts_guess.coords[SI_INDEX] - ts_guess.coords[ow_index])
            )
            log(
                "  saved NEB peak is outside the channel "
                f"(r(Si-Ow)={r_bad:.2f} A); rebuilding from the closer endpoint"
            )
            ts_guess_path.replace(run_dir / "ts_guess.rejected-neb-escape.xyz")
            stale_traj = run_dir / "sella.traj"
            if stale_traj.exists():
                stale_traj.replace(run_dir / "sella.rejected-neb-escape.traj")
            ts_guess = proton_neb_guess(
                ts_guess, complex_opt, settings, run_dir, ow_index
            )
            save_xyz(ts_guess, ts_guess_path)
    elif route == "proton-neb" and (run_dir / "product.xyz").exists():
        # A failed bounded NEB leaves the route and validated product durable
        # but no ts_guess.xyz.  Resume that exact stage instead of repeating
        # the already-completed direct/proton scans.
        log("  resume: proton-NEB route/product exist; rebuilding only the band")
        ts_guess = proton_neb_guess(
            complex_opt, complex_opt, settings, run_dir, ow_index
        )
        save_xyz(ts_guess, ts_guess_path)
    else:
        # No reusable product exists. Reconstruct from the direct route.
        route = "direct"
        route_path.unlink(missing_ok=True)
        try:
            scan = scan_to_maximum(
                complex_opt,
                settings,
                atom_i=SI_INDEX,
                atom_j=ow_index,
                distances_a=[2.8, 2.6, 2.4, 2.2, 2.1, 2.0, 1.9],
                progress=lambda r, e: log(f"  r={r:.2f} A  E={e:.6f} Ha"),
            )
            ts_guess = scan_ts_guess(scan)
        except ScanNoMaximumError as exc:
            log("  approach coordinate alone does not cross the ridge;")
            base = min(exc.scan, key=lambda p: abs(p[0] - 1.90))[2]
            ts_guess = proton_neb_guess(base, complex_opt, settings, run_dir, ow_index)
            route = "proton-neb"
        save_xyz(ts_guess, ts_guess_path)
        if route == "proton-neb":
            route_path.write_text(route)

    # Stage 3 — saddle search. A direct approach crest gets one bounded
    # fallback through the concerted proton/product/CI-NEB route if Sella
    # escapes. A proton-NEB saddle that escapes is terminal.
    log("stage 3: Sella saddle search")
    ts_path = run_dir / "ts.xyz"
    trajectory_path = run_dir / "sella.traj"
    ts = checkpointed(
        ts_path,
        ts_guess,
        lambda: find_ts(ts_guess, settings, trajectory=str(trajectory_path)),
    )

    r_guess = float(
        np.linalg.norm(ts_guess.coords[SI_INDEX] - ts_guess.coords[ow_index])
    )
    r_ts = float(np.linalg.norm(ts.coords[SI_INDEX] - ts.coords[ow_index]))
    log(f"  r(Si-Ow): guess {r_guess:.2f} A -> saddle {r_ts:.2f} A")
    escape = channel_escape_reason(ts_guess, ts, ow_index)
    if escape and route == "direct":
        log(f"  !! {escape}; pivoting once to proton-drive + product + CI-NEB")
        neb_guess = proton_neb_guess(ts_guess, complex_opt, settings, run_dir, ow_index)

        # Preserve the rejected direct-search receipt, then replace all
        # resumable checkpoints with the new route before refining again.
        rejected_ts = run_dir / "ts.rejected-direct.xyz"
        if ts_path.exists():
            ts_path.replace(rejected_ts)
        rejected_traj = run_dir / "sella.rejected-direct.traj"
        if trajectory_path.exists():
            trajectory_path.replace(rejected_traj)
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
        r_guess = float(
            np.linalg.norm(ts_guess.coords[SI_INDEX] - ts_guess.coords[ow_index])
        )
        r_ts = float(np.linalg.norm(ts.coords[SI_INDEX] - ts.coords[ow_index]))
        log(f"  r(Si-Ow): NEB guess {r_guess:.2f} A -> saddle {r_ts:.2f} A")
        escape = channel_escape_reason(ts_guess, ts, ow_index)
    if escape:
        log(f"  !! {escape}; {route} route exhausted — aborting before thermochemistry")
        return 1

    # Stage 4 — verify + quick IRC.
    log("stage 4: frequencies + quick-IRC")
    ts_freq = frequencies(ts, settings)
    log(f"  TS imaginary modes: {np.round(ts_freq.imaginary_cm, 1).tolist()}")
    if ts_freq.n_imaginary != 1:
        log("  !! not a clean first-order saddle — inspect ts.xyz; aborting")
        return 1
    back, fwd = quick_irc(
        ts,
        settings,
        displacement_a=0.50,
        frequency=ts_freq,
    )
    save_xyz(back, run_dir / "irc_back.xyz")
    save_xyz(fwd, run_dir / "irc_fwd.xyz")
    irc_failure = quick_irc_channel_reason(back, fwd, ow_index)
    if irc_failure:
        cleavage_failure = quick_irc_cleavage_reason(back, fwd, ow_index)
        if args.reaction == "al-neutral" and cleavage_failure is None:
            log("  full path resolved as sequential associative hydrolysis")
            return finish_al_neutral_sequential(
                complex_opt=complex_opt,
                cleavage_ts=ts,
                cleavage_freq=ts_freq,
                cleavage_back=back,
                cleavage_fwd=fwd,
                settings=settings,
                run_dir=run_dir,
                ow_index=ow_index,
                temperature=args.temperature,
                dimer_opt=dimer_opt,
                attacker_opt=attacker_opt,
            )
        log(f"  !! {irc_failure}; aborting before thermochemistry")
        return 1

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
