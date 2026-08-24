#!/usr/bin/env python
"""Bounded al-neutral associative-intermediate conformer/proton sweep.

TASK-168 stage 2 asks whether the high-lying associative intermediate used by
its sequential mechanism is merely a bad proton/hydroxyl conformer. This
runner generates a finite, deterministic set of one-coordinate hydrogen
rotamers, optimizes each at the requested DFT tier, rejects candidates that
leave the associative basin, and atomically journals every result.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from dataclasses import asdict, replace
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    bootstrap_cli(
        "task168_intermediate_sweep",
        default_run_root=Path(__file__).resolve().parent.parent / "runs",
    )

import numpy as np  # noqa: E402

from quarry.clusters import (  # noqa: E402
    Cluster,
    aluminosilicate_dimer,
    hydrolysis_complex,
    water,
)
from quarry.pipeline import HARTREE_TO_KJ, DftSettings, energy, optimize  # noqa: E402
from scripts.phase1_xiao_lasaga import (  # noqa: E402
    AL_INDEX,
    ASSOCIATIVE_BASIN,
    BR_INDEX,
    SI_INDEX,
    hydrolysis_basin_signature,
    load_xyz,
    save_xyz,
)

DEFAULT_RUN_DIR = (
    Path(__file__).resolve().parent.parent
    / "runs"
    / "phase1"
    / "al-neutral-b3lyp-def2-svp-flank"
)
TERMINAL_HYDROXYLS = (
    (SI_INDEX, 3, 4),
    (SI_INDEX, 5, 6),
    (SI_INDEX, 7, 8),
    (AL_INDEX, 9, 10),
    (AL_INDEX, 11, 12),
    (AL_INDEX, 13, 14),
)
SIGNIFICANT_LOWERING_KJ = 2.0
MIN_PAIR_DISTANCE_A = 0.55


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def rotate_point_about_axis(
    point: np.ndarray,
    axis_start: np.ndarray,
    axis_end: np.ndarray,
    angle_degrees: float,
) -> np.ndarray:
    """Rotate one Cartesian point around an oriented axis with Rodrigues' rule."""
    axis = np.asarray(axis_end, dtype=float) - np.asarray(axis_start, dtype=float)
    norm = float(np.linalg.norm(axis))
    if not np.isfinite(norm) or norm < 1e-12:
        raise ValueError("rotation axis is non-finite or degenerate")
    unit = axis / norm
    vector = np.asarray(point, dtype=float) - np.asarray(axis_end, dtype=float)
    angle = np.deg2rad(angle_degrees)
    rotated = (
        vector * np.cos(angle)
        + np.cross(unit, vector) * np.sin(angle)
        + unit * np.dot(unit, vector) * (1.0 - np.cos(angle))
    )
    return np.asarray(axis_end, dtype=float) + rotated


def minimum_pair_distance(cluster: Cluster) -> float:
    delta = cluster.coords[:, None, :] - cluster.coords[None, :, :]
    distances = np.linalg.norm(delta, axis=2)
    distances[np.diag_indices_from(distances)] = np.inf
    return float(np.min(distances))


def intermediate_conformer_seeds(
    intermediate: Cluster,
    ow_index: int,
) -> list[tuple[str, Cluster]]:
    """Return the finite TASK-168 hydrogen-rotamer diagnostic set.

    Each candidate changes one interpretable coordinate except the deliberate
    all-terminal inversion. Heavy atoms and electronic state are byte-identical.
    This is a diagnostic sweep, not an unconstrained combinatorial search.
    """

    def rotated(label: str, rotations: list[tuple[int, int, int, float]]) -> Cluster:
        coords = intermediate.coords.copy()
        for center, oxygen, hydrogen, angle in rotations:
            coords[hydrogen] = rotate_point_about_axis(
                coords[hydrogen], coords[center], coords[oxygen], angle
            )
        return replace(intermediate, name=f"{intermediate.name}-{label}", coords=coords)

    seeds: list[tuple[str, Cluster]] = [("baseline", intermediate)]
    for _center, oxygen, hydrogen in TERMINAL_HYDROXYLS:
        label = f"terminal-oh-{oxygen:02d}-180"
        seeds.append((label, rotated(label, [(_center, oxygen, hydrogen, 180.0)])))
    all_terminal = [(*triple, 180.0) for triple in TERMINAL_HYDROXYLS]
    seeds.append(("all-terminal-oh-180", rotated("all-terminal-oh-180", all_terminal)))

    water_h = (ow_index + 1, ow_index + 2)
    transfer_h = min(
        water_h,
        key=lambda index: np.linalg.norm(
            intermediate.coords[index] - intermediate.coords[BR_INDEX]
        ),
    )
    residual_h = next(index for index in water_h if index != transfer_h)
    for angle in (-120.0, 120.0):
        label = f"bridge-proton-{int(angle):+d}"
        seeds.append(
            (
                label,
                rotated(label, [(SI_INDEX, BR_INDEX, transfer_h, angle)]),
            )
        )
    seeds.append(
        (
            "water-proton-180",
            rotated("water-proton-180", [(SI_INDEX, ow_index, residual_h, 180.0)]),
        )
    )
    return seeds


def run_sweep(run_dir: Path, settings: DftSettings, *, max_steps: int) -> dict:
    template = hydrolysis_complex(aluminosilicate_dimer(), water(), mode="flank")
    intermediate_path = run_dir / "intermediate.xyz"
    reactant_path = run_dir / "complex.xyz"
    if not intermediate_path.is_file() or not reactant_path.is_file():
        raise FileNotFoundError(
            "TASK-168 requires existing complex.xyz and intermediate.xyz"
        )
    intermediate = load_xyz(intermediate_path, template)
    reactant = load_xyz(reactant_path, template)
    ow_index = len(aluminosilicate_dimer().symbols)
    if hydrolysis_basin_signature(intermediate, ow_index) != ASSOCIATIVE_BASIN:
        raise RuntimeError("saved intermediate.xyz is not in the associative basin")

    sweep_dir = run_dir / "task168-intermediate-sweep"
    sweep_dir.mkdir(exist_ok=True)
    manifest_path = sweep_dir / "manifest.json"
    manifest = (
        json.loads(manifest_path.read_text())
        if manifest_path.is_file()
        else {
            "schema_version": 1,
            "settings": asdict(settings),
            "max_steps": max_steps,
            "significant_lowering_kj": SIGNIFICANT_LOWERING_KJ,
            "candidates": {},
        }
    )
    if manifest["settings"] != asdict(settings) or manifest["max_steps"] != max_steps:
        raise RuntimeError("existing sweep manifest belongs to different settings")

    reactant_energy = energy(reactant, settings)
    seeds = intermediate_conformer_seeds(intermediate, ow_index)
    for label, seed in seeds:
        existing = manifest["candidates"].get(label)
        if existing and existing.get("status") in {"accepted", "rejected", "failed"}:
            print(
                f"TASK168_SWEEP_RESUME label={label} status={existing['status']}",
                flush=True,
            )
            continue
        seed_path = sweep_dir / f"{label}.seed.xyz"
        result_path = sweep_dir / f"{label}.optimized.xyz"
        save_xyz(seed, seed_path)
        pair_distance = minimum_pair_distance(seed)
        if pair_distance < MIN_PAIR_DISTANCE_A:
            manifest["candidates"][label] = {
                "status": "rejected",
                "reason": f"seed minimum pair distance {pair_distance:.6f} A",
                "seed": seed_path.name,
            }
            atomic_json(manifest_path, manifest)
            continue
        print(f"TASK168_SWEEP_START label={label}", flush=True)
        try:
            optimized = optimize(seed, settings, max_steps=max_steps)
            save_xyz(optimized, result_path)
            electronic = energy(optimized, settings)
            signature = hydrolysis_basin_signature(optimized, ow_index)
            accepted = signature == ASSOCIATIVE_BASIN
            record = {
                "status": "accepted" if accepted else "rejected",
                "reason": None if accepted else f"basin signature {signature}",
                "seed": seed_path.name,
                "optimized": result_path.name,
                "electronic_hartree": electronic,
                "electronic_kj_vs_reactant": (electronic - reactant_energy)
                * HARTREE_TO_KJ,
                "minimum_pair_distance_a": minimum_pair_distance(optimized),
                "basin_signature": list(signature),
            }
        except Exception as exc:  # each bounded seed must not erase earlier evidence
            record = {
                "status": "failed",
                "reason": f"{type(exc).__name__}: {exc}",
                "seed": seed_path.name,
            }
        manifest["candidates"][label] = record
        atomic_json(manifest_path, manifest)
        print(
            f"TASK168_SWEEP_DONE label={label} status={record['status']} "
            f"reason={record.get('reason')}",
            flush=True,
        )

    accepted = {
        label: record
        for label, record in manifest["candidates"].items()
        if record["status"] == "accepted"
    }
    if "baseline" not in accepted:
        raise RuntimeError(
            "baseline intermediate did not survive its own minimum/basin gate"
        )
    best_label, best = min(
        accepted.items(), key=lambda item: item[1]["electronic_hartree"]
    )
    baseline = accepted["baseline"]
    lowering_kj = (
        baseline["electronic_hartree"] - best["electronic_hartree"]
    ) * HARTREE_TO_KJ
    selected_label = best_label if lowering_kj > SIGNIFICANT_LOWERING_KJ else "baseline"
    selected_record = accepted[selected_label]
    selected = load_xyz(sweep_dir / selected_record["optimized"], template)
    save_xyz(selected, run_dir / "intermediate.task168-selected.xyz")
    manifest["summary"] = {
        "reactant_electronic_hartree": reactant_energy,
        "baseline_label": "baseline",
        "best_label": best_label,
        "best_lowering_vs_baseline_kj": lowering_kj,
        "selected_label": selected_label,
        "selected_electronic_hartree": selected_record["electronic_hartree"],
        "selected_electronic_kj_vs_reactant": selected_record[
            "electronic_kj_vs_reactant"
        ],
        "lower_intermediate_found": selected_label != "baseline",
        "selected_structure": "../intermediate.task168-selected.xyz",
    }
    atomic_json(manifest_path, manifest)
    print(
        "TASK168_SWEEP_SUMMARY " + json.dumps(manifest["summary"], sort_keys=True),
        flush=True,
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--nice", type=int, default=10)
    parser.add_argument("--log")
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    parser.add_argument("--xc", default="b3lyp")
    parser.add_argument("--basis", default="def2-svp")
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--max-steps", type=int, default=180)
    args = parser.parse_args()
    if args.max_steps <= 0 or args.max_steps > 240:
        parser.error("--max-steps must be in 1..240")
    settings = DftSettings(
        xc=args.xc,
        basis=args.basis,
        density_fit=True,
        use_gpu=args.gpu,
    )
    run_sweep(args.run_dir, settings, max_steps=args.max_steps)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
