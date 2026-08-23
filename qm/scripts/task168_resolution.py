#!/usr/bin/env python
"""Resolve TASK-168 after the bounded associative-intermediate sweep.

If the original high-lying intermediate survives the diagnostic sweep, this
executes exactly one densified, checkpointed BFGS climbing band and one Sella
refinement. Failure to produce an in-channel first-order R<->I saddle is the
pre-approved mechanistic finding: a quasi-barrierless uphill addition slide,
with the proven I<->P cleavage saddle closing the reactant-referenced profile.
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
    hydrolysis_complex,
    water,
)
from quarry.pipeline import (  # noqa: E402
    DftSettings,
    FrequencyResult,
    frequencies,
    frequency_geometry_fingerprint,
    frequency_settings_fingerprint,
)
from quarry.rates import rate_from_thermo  # noqa: E402
from quarry.store import Store  # noqa: E402
from quarry.ts import (  # noqa: E402
    find_ts,
    neb_ts_guess,
    quick_irc,
    reaction_path_vector,
)
from scripts.phase1_xiao_lasaga import (  # noqa: E402
    AL_INDEX,
    ASSOCIATIVE_BASIN,
    BR_INDEX,
    HYDROLYZED_BASIN,
    SI_INDEX,
    SI_NEUTRAL_DG_KCAL,
    SI_NEUTRAL_DG_KJ,
    addition_channel_escape_reason,
    begin_sequential_run,
    endpoint_in_basin,
    hydrolysis_basin_signature,
    load_xyz,
    minimum_frequency_reason,
    quick_irc_addition_reason,
    quick_irc_cleavage_reason,
    save_xyz,
    thermo_result,
    write_sequential_run_status,
)

DEFAULT_RUN_DIR = (
    Path(__file__).resolve().parent.parent
    / "runs"
    / "phase1"
    / "al-neutral-b3lyp-def2-svp-flank"
)
KCAL = 4.184
OBSERVED_ADDITION_CREST_ELECTRONIC_KJ = 150.471


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def densified_piecewise_band(
    reactant: Cluster,
    crest: Cluster,
    intermediate: Cluster,
) -> list[np.ndarray]:
    """Nine images with the final five progressively densified toward I."""
    if not (reactant.symbols == crest.symbols == intermediate.symbols):
        raise ValueError("piecewise band atom mappings differ")
    coords = []
    for fraction in (0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0):
        coords.append(reactant.coords + fraction * (crest.coords - reactant.coords))
    for fraction in (0.35, 0.60, 0.78, 0.91, 1.0):
        coords.append(crest.coords + fraction * (intermediate.coords - crest.coords))
    return [np.asarray(image, dtype=float).copy() for image in coords]


def save_frequency(path: Path, result: FrequencyResult) -> None:
    temporary = path.with_name(f".{path.name}.{time.time_ns()}.tmp.npz")
    np.savez(
        temporary,
        frequencies_cm=result.frequencies_cm,
        imaginary_cm=result.imaginary_cm,
        electronic_hartree=np.array(result.electronic_hartree),
        molar_mass_kg=np.array(result.molar_mass_kg),
        rotational_temperatures_k=np.array(result.rotational_temperatures_k or ()),
        linear=np.array(result.linear),
        imaginary_mode=(
            result.imaginary_mode
            if result.imaginary_mode is not None
            else np.empty((0, 3))
        ),
        imaginary_modes=(
            result.imaginary_modes
            if result.imaginary_modes is not None
            else np.empty((0, 0, 3))
        ),
        geometry_fingerprint=np.array(result.geometry_fingerprint),
        settings_fingerprint=np.array(result.settings_fingerprint),
    )
    temporary.replace(path)


def load_frequency(
    path: Path, cluster: Cluster, settings: DftSettings
) -> FrequencyResult | None:
    if not path.is_file():
        return None
    with np.load(path, allow_pickle=False) as data:
        geometry_fingerprint = str(data["geometry_fingerprint"].item())
        settings_fingerprint = str(data["settings_fingerprint"].item())
        if geometry_fingerprint != frequency_geometry_fingerprint(cluster):
            return None
        if settings_fingerprint != frequency_settings_fingerprint(settings):
            return None
        rotational = tuple(float(value) for value in data["rotational_temperatures_k"])
        imaginary_mode = np.asarray(data["imaginary_mode"])
        imaginary_modes = (
            np.asarray(data["imaginary_modes"])
            if "imaginary_modes" in data.files
            else np.empty((0, 0, 3))
        )
        return FrequencyResult(
            frequencies_cm=np.asarray(data["frequencies_cm"]),
            imaginary_cm=np.asarray(data["imaginary_cm"]),
            electronic_hartree=float(data["electronic_hartree"]),
            molar_mass_kg=float(data["molar_mass_kg"]),
            rotational_temperatures_k=rotational or None,
            linear=bool(data["linear"]),
            imaginary_mode=imaginary_mode if imaginary_mode.size else None,
            imaginary_modes=imaginary_modes if imaginary_modes.size else None,
            geometry_fingerprint=geometry_fingerprint,
            settings_fingerprint=settings_fingerprint,
        )


def cached_frequency(
    path: Path, cluster: Cluster, settings: DftSettings
) -> FrequencyResult:
    cached = load_frequency(path, cluster, settings)
    if cached is not None:
        print(f"TASK168_FREQ_RESUME path={path.name}", flush=True)
        return cached
    if path.exists():
        path.replace(
            path.with_name(f"{path.stem}.rejected-{time.time_ns()}{path.suffix}")
        )
    result = frequencies(cluster, settings)
    save_frequency(path, result)
    print(
        f"TASK168_FREQ_DONE path={path.name} n_imaginary={result.n_imaginary} "
        f"electronic_hartree={result.electronic_hartree:.12f}",
        flush=True,
    )
    return result


def newest_pre_relaxed_checkpoint(checkpoint_root: Path) -> Path | None:
    candidates = []
    if not checkpoint_root.is_dir():
        return None
    for path in checkpoint_root.iterdir():
        manifest_path = path / "manifest.json"
        if not manifest_path.is_file():
            continue
        manifest = json.loads(manifest_path.read_text())
        if (
            manifest.get("stage") == "pre-relax-final"
            and manifest.get("converged") is True
        ):
            candidates.append(path)
    return max(candidates, key=lambda path: path.stat().st_mtime_ns, default=None)


def attempt_addition_saddle(
    run_dir: Path,
    reactant: Cluster,
    intermediate: Cluster,
    settings: DftSettings,
    ow_index: int,
) -> tuple[Cluster, FrequencyResult, Cluster, Cluster] | None:
    receipt_path = run_dir / "task168-addition-attempt.json"
    if receipt_path.is_file():
        receipt = json.loads(receipt_path.read_text())
        if receipt.get("status") == "no-saddle":
            print("TASK168_ADDITION_ATTEMPT_ALREADY_EXHAUSTED", flush=True)
            return None
        if receipt.get("status") == "accepted":
            ts = load_xyz(run_dir / receipt["ts"], reactant)
            back = load_xyz(run_dir / receipt["irc_back"], reactant)
            fwd = load_xyz(run_dir / receipt["irc_fwd"], reactant)
            freq = cached_frequency(run_dir / receipt["frequency"], ts, settings)
            return ts, freq, back, fwd

    crest = load_xyz(run_dir / "addition_pinned_relaxed.xyz", reactant)
    checkpoint_root = run_dir / "task168-addition-neb-checkpoints"
    pre_relaxed = newest_pre_relaxed_checkpoint(checkpoint_root)
    try:
        if pre_relaxed is not None:
            print(f"TASK168_NEB_RESUME checkpoint={pre_relaxed}", flush=True)
            guess = neb_ts_guess(
                reactant,
                intermediate,
                settings,
                n_images=9,
                fmax_ev_a=0.05,
                max_steps=180,
                pre_relax_fmax_ev_a=0.20,
                pre_relax_steps=180,
                optimizer_dt=0.02,
                optimizer_maxstep=0.03,
                climb_optimizer="bfgs",
                checkpoint_dir=checkpoint_root,
                checkpoint_interval=5,
                resume_from=pre_relaxed,
            )
        else:
            guess = neb_ts_guess(
                reactant,
                intermediate,
                settings,
                n_images=9,
                fmax_ev_a=0.05,
                max_steps=180,
                pre_relax_fmax_ev_a=0.20,
                pre_relax_steps=180,
                optimizer_dt=0.02,
                optimizer_maxstep=0.03,
                climb_optimizer="bfgs",
                checkpoint_dir=checkpoint_root,
                checkpoint_interval=5,
                initial_image_coords=densified_piecewise_band(
                    reactant, crest, intermediate
                ),
            )
        save_xyz(guess, run_dir / "addition_task168_neb_guess.xyz")
        reactive_core = sorted(
            {BR_INDEX, SI_INDEX, AL_INDEX, ow_index, ow_index + 1, ow_index + 2}
        )
        ts = find_ts(
            guess,
            settings,
            trajectory=str(run_dir / "addition_task168_sella.traj"),
            initial_mode=reaction_path_vector(
                reactant, intermediate, active_indices=reactive_core
            ),
            internal=False,
            max_steps=300,
        )
        save_xyz(ts, run_dir / "addition_task168_ts.xyz")
        if reason := addition_channel_escape_reason(
            guess, ts, ow_index, acid_path=False
        ):
            raise RuntimeError(reason)
        freq = cached_frequency(
            run_dir / "addition_task168_ts.frequency.npz", ts, settings
        )
        if freq.n_imaginary != 1:
            raise RuntimeError(
                f"TASK-168 addition saddle has {freq.n_imaginary} imaginary modes"
            )
        accepted = None
        for displacement in (0.50, 1.00):
            back, fwd = quick_irc(
                ts,
                settings,
                displacement_a=displacement,
                frequency=freq,
            )
            if quick_irc_addition_reason(back, fwd, ow_index) is None:
                accepted = (back, fwd, displacement)
                break
        if accepted is None:
            raise RuntimeError("TASK-168 addition saddle does not connect R and I")
        back, fwd, displacement = accepted
        save_xyz(back, run_dir / "addition_task168_irc_back.xyz")
        save_xyz(fwd, run_dir / "addition_task168_irc_fwd.xyz")
        receipt = {
            "status": "accepted",
            "ts": "addition_task168_ts.xyz",
            "frequency": "addition_task168_ts.frequency.npz",
            "irc_back": "addition_task168_irc_back.xyz",
            "irc_fwd": "addition_task168_irc_fwd.xyz",
            "quick_irc_displacement_a": displacement,
        }
        atomic_json(receipt_path, receipt)
        return ts, freq, back, fwd
    except Exception as exc:
        receipt = {
            "status": "no-saddle",
            "reason": f"{type(exc).__name__}: {exc}",
            "checkpoint_root": checkpoint_root.name,
            "kill_criterion": "no attempt #5; report sequential mechanism finding",
        }
        atomic_json(receipt_path, receipt)
        print(
            "TASK168_ADDITION_NO_SADDLE " + json.dumps(receipt, sort_keys=True),
            flush=True,
        )
        return None


def finalize(
    run_dir: Path,
    settings: DftSettings,
    temperature: float,
    addition: tuple[Cluster, FrequencyResult, Cluster, Cluster] | None,
) -> dict:
    template = hydrolysis_complex(aluminosilicate_dimer(), water(), mode="flank")
    ow_index = len(aluminosilicate_dimer().symbols)
    reactant = load_xyz(run_dir / "complex.xyz", template)
    intermediate = load_xyz(run_dir / "intermediate.task168-selected.xyz", template)
    cleavage_ts = load_xyz(run_dir / "cleavage_ts.xyz", template)
    cleavage_back = load_xyz(run_dir / "cleavage_irc_back.xyz", template)
    cleavage_fwd = load_xyz(run_dir / "cleavage_irc_fwd.xyz", template)
    if quick_irc_cleavage_reason(cleavage_back, cleavage_fwd, ow_index):
        raise RuntimeError("preserved cleavage saddle no longer passes I<->P gate")
    endpoint_in_basin(cleavage_back, cleavage_fwd, ow_index, ASSOCIATIVE_BASIN)
    endpoint_in_basin(cleavage_back, cleavage_fwd, ow_index, HYDROLYZED_BASIN)

    begin_sequential_run(run_dir)
    reactant_freq = cached_frequency(
        run_dir / "complex.task168.frequency.npz", reactant, settings
    )
    intermediate_freq = cached_frequency(
        run_dir / "intermediate.task168.frequency.npz", intermediate, settings
    )
    cleavage_freq = cached_frequency(
        run_dir / "cleavage_ts.task168.frequency.npz", cleavage_ts, settings
    )
    for label, freq in (
        ("reactant", reactant_freq),
        ("intermediate", intermediate_freq),
    ):
        if reason := minimum_frequency_reason(label, freq):
            write_sequential_run_status(run_dir, "failed", detail=reason)
            raise RuntimeError(reason)
    if cleavage_freq.n_imaginary != 1:
        reason = f"cleavage saddle has {cleavage_freq.n_imaginary} imaginary modes"
        write_sequential_run_status(run_dir, "failed", detail=reason)
        raise RuntimeError(reason)

    th_reactant = thermo_result(reactant_freq, temperature)
    th_intermediate = thermo_result(intermediate_freq, temperature)
    th_cleavage = thermo_result(cleavage_freq, temperature)
    cleavage_rate = rate_from_thermo(
        th_intermediate,
        th_cleavage,
        imag_nu_cm=float(cleavage_freq.imaginary_cm[0]),
        tunneling="wigner",
    )
    profile_rate = rate_from_thermo(
        th_reactant,
        th_cleavage,
        imag_nu_cm=float(cleavage_freq.imaginary_cm[0]),
        tunneling="wigner",
    )
    intermediate_dg = th_intermediate.gibbs - th_reactant.gibbs
    overall_dg = th_cleavage.gibbs - th_reactant.gibbs
    addition_attempt = json.loads(
        (run_dir / "task168-addition-attempt.json").read_text()
    )
    sweep = json.loads(
        (run_dir / "task168-intermediate-sweep" / "manifest.json").read_text()
    )["summary"]

    results: dict = {
        "reaction": "al-neutral",
        "method": f"{settings.xc}/{settings.basis}/df",
        "temperature_k": temperature,
        "mechanism": (
            "sequential-associative"
            if addition is not None
            else "sequential-associative-uphill-slide-no-resolved-addition-saddle"
        ),
        "mechanism_finding": addition is None,
        "intermediate_sweep": sweep,
        "profile": {
            "intermediate_dG_vs_reactant_kj": intermediate_dg,
            "cleavage_dG_dagger_from_intermediate_kj": cleavage_rate.dg_kj,
            "cleavage_ts_dG_vs_reactant_kj": overall_dg,
            "overall_profile_dG_dagger_kj": overall_dg,
            "highest_profile_ts": "cleavage",
        },
        "dH_kj": profile_rate.dh_kj,
        "dG_kj": overall_dg,
        "dG_kcal": overall_dg / KCAL,
        "wigner_kappa": profile_rate.kappa,
        "k_per_s": profile_rate.k,
        "ts_imaginary_cm": float(cleavage_freq.imaginary_cm[0]),
        "steps": {
            "cleavage": {
                "dH_kj": cleavage_rate.dh_kj,
                "dG_kj": cleavage_rate.dg_kj,
                "dG_kcal": cleavage_rate.dg_kj / KCAL,
                "wigner_kappa": cleavage_rate.kappa,
                "k_per_s": cleavage_rate.k,
                "ts_imaginary_cm": float(cleavage_freq.imaginary_cm[0]),
                "quick_irc": "associative-intermediate<->hydrolyzed-product",
            }
        },
        "addition_attempt": addition_attempt,
        "comparison": {
            "si_neutral_dG_kj": SI_NEUTRAL_DG_KJ,
            "si_neutral_dG_kcal": SI_NEUTRAL_DG_KCAL,
            "al_easier_than_si": overall_dg < SI_NEUTRAL_DG_KJ,
            "xiao_lasaga_easier_si_o_al_ordering_confirmed": overall_dg
            < SI_NEUTRAL_DG_KJ,
        },
        "absolute_numbers_require_a2_retier": True,
    }
    if addition is None:
        results["steps"]["addition"] = {
            "classification": "quasi-barrierless-uphill-slide",
            "accepted_first_order_saddle": False,
            "observed_electronic_crest_lower_bound_kj": (
                OBSERVED_ADDITION_CREST_ELECTRONIC_KJ
            ),
            "note": "No attempt #5 under the Victor-approved TASK-168 kill criterion.",
        }
    else:
        addition_ts, addition_freq, addition_back, addition_fwd = addition
        if quick_irc_addition_reason(addition_back, addition_fwd, ow_index):
            raise RuntimeError("accepted addition handoff failed R<->I recheck")
        th_addition = thermo_result(addition_freq, temperature)
        addition_rate = rate_from_thermo(
            th_reactant,
            th_addition,
            imag_nu_cm=float(addition_freq.imaginary_cm[0]),
            tunneling="wigner",
        )
        results["steps"]["addition"] = {
            "dH_kj": addition_rate.dh_kj,
            "dG_kj": addition_rate.dg_kj,
            "dG_kcal": addition_rate.dg_kj / KCAL,
            "wigner_kappa": addition_rate.kappa,
            "k_per_s": addition_rate.k,
            "ts_imaginary_cm": float(addition_freq.imaginary_cm[0]),
        }
        save_xyz(addition_ts, run_dir / "addition_task168_ts.xyz")

    results_tmp = run_dir / "results.json.tmp"
    results_tmp.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n")
    store_tmp = run_dir / "store.task168.tmp.sqlite"
    store_tmp.unlink(missing_ok=True)
    structures = [
        ("complex", reactant, reactant_freq),
        ("intermediate", intermediate, intermediate_freq),
        ("cleavage-ts", cleavage_ts, cleavage_freq),
    ]
    if addition is not None:
        structures.append(("addition-ts", addition[0], addition[1]))
    with Store(store_tmp) as store:
        for name, cluster, freq in structures:
            structure_id = store.add_structure(
                f"al-neutral-{name}",
                cluster.formula,
                cluster.to_xyz(),
                charge=cluster.charge,
                spin=cluster.spin,
            )
            job_id = store.add_job(
                structure_id,
                "freq",
                str(results["method"]),
                "gpu4pyscf" if settings.use_gpu else "pyscf",
            )
            store.set_job_status(job_id, "done")
            store.add_result(job_id, "electronic", freq.electronic_hartree, "hartree")
    store_tmp.replace(run_dir / "store.sqlite")
    results_tmp.replace(run_dir / "results.json")
    write_sequential_run_status(
        run_dir,
        "completed",
        detail=(
            "accepted addition+cleavage saddles"
            if addition is not None
            else (
                "Victor-approved mechanistic finding after bounded addition "
                "kill criterion"
            )
        ),
    )
    print("TASK168_FINAL " + json.dumps(results, sort_keys=True), flush=True)
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--temperature", type=float, default=298.15)
    args = parser.parse_args()
    sweep_manifest = json.loads(
        (args.run_dir / "task168-intermediate-sweep" / "manifest.json").read_text()
    )
    summary = sweep_manifest["summary"]
    if summary["lower_intermediate_found"]:
        raise RuntimeError(
            "a significantly lower intermediate was found; rebuild both elementary "
            "bands through the standard ladder before using the old cleavage saddle"
        )
    template = hydrolysis_complex(aluminosilicate_dimer(), water(), mode="flank")
    reactant = load_xyz(args.run_dir / "complex.xyz", template)
    intermediate = load_xyz(
        args.run_dir / "intermediate.task168-selected.xyz", template
    )
    ow_index = len(aluminosilicate_dimer().symbols)
    if hydrolysis_basin_signature(intermediate, ow_index) != ASSOCIATIVE_BASIN:
        raise RuntimeError("selected intermediate is outside the associative basin")
    settings = DftSettings(
        xc="b3lyp", basis="def2-svp", density_fit=True, use_gpu=args.gpu
    )
    addition = attempt_addition_saddle(
        args.run_dir, reactant, intermediate, settings, ow_index
    )
    finalize(args.run_dir, settings, args.temperature, addition)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
