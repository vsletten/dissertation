#!/usr/bin/env python3
"""Localize the A2a cleavage saddle with one bounded full-system dimer run."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    _ETIQUETTE = bootstrap_cli(
        "a2a-cleavage-localize",
        default_run_root="/mnt/data/vsletten/dissertation-data/task208-a2a-path-rebuild",
    )

import numpy as np

from quarry.pipeline import (
    frequency_geometry_fingerprint,
    frequency_settings_fingerprint,
)
from quarry.ts import find_ts_dimer
from scripts import a2a_path_rebuild as a2a
from scripts import production_energetics as a2

LOCALIZATION_VERSION = "a2a-cleavage-full-system-local-dimer-v2"


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    a2.atomic_json(path, payload)


def localization_inputs(
    checkpoint_root: Path,
    neb_crest,
    conditioned_crest,
    local_indices: list[int],
) -> tuple[np.ndarray, dict[str, Any]]:
    """Bind the mode and envelope to the exact persisted adjacent images."""
    mode, tangent = a2a.final_neb_tangent(
        checkpoint_root,
        neb_crest,
        conditioned_crest,
        local_indices,
        strategy=LOCALIZATION_VERSION,
    )
    checkpoint = Path(tangent["checkpoint"])
    left = a2a.load_reference(
        checkpoint / tangent["left_image"],
        role="left-neighbor",
        template=conditioned_crest,
    )
    right = a2a.load_reference(
        checkpoint / tangent["right_image"],
        role="right-neighbor",
        template=conditioned_crest,
    )
    local = np.asarray(local_indices, dtype=int)
    conditioned_radii = [
        float(np.linalg.norm(neighbor.coords[local] - conditioned_crest.coords[local]))
        for neighbor in (left, right)
    ]
    peak_radii = [
        float(np.linalg.norm(neighbor.coords[local] - neb_crest.coords[local]))
        for neighbor in (left, right)
    ]
    if any(not np.isfinite(radius) or radius <= 0.0 for radius in peak_radii):
        raise ValueError("NEB crest has an invalid adjacent-image radius")
    if any(not np.isfinite(radius) or radius <= 0.0 for radius in conditioned_radii):
        raise ValueError("conditioned crest has an invalid adjacent-image radius")
    # The conditioner relaxes spectators and can sit much farther from the band
    # images than the exact climb peak.  Measuring from it widened v1's envelope
    # from 0.344/0.546 A to 1.533/1.643 A.  Keep the exact persisted peak-neighbor
    # radii as the contract; the conditioned crest is only the search origin.
    trust_radius = min(peak_radii)
    guard_radius = max(max(peak_radii), 1.5 * trust_radius)
    receipt = {
        **tangent,
        "peak_left_local_radius_a": peak_radii[0],
        "peak_right_local_radius_a": peak_radii[1],
        "conditioned_left_local_radius_a": conditioned_radii[0],
        "conditioned_right_local_radius_a": conditioned_radii[1],
        "local_trust_radius_a": trust_radius,
        "local_guard_radius_a": guard_radius,
        "envelope_origin": "exact-persisted-climb-peak-adjacent-images",
        "checkpoint_manifest_sha256": a2.sha256_path(checkpoint / "manifest.json"),
        "left_image_sha256": a2.sha256_path(checkpoint / tangent["left_image"]),
        "right_image_sha256": a2.sha256_path(checkpoint / tangent["right_image"]),
    }
    return mode, receipt


def run(args: argparse.Namespace) -> int:
    run_dir = args.run_dir.resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    a2a.preflight_gpu_contraction_engine(args.gpu)
    conditioned_path = args.conditioned_crest.resolve()
    neb_crest_path = args.neb_crest.resolve()
    checkpoint_root = args.checkpoint_root.resolve()
    conditioned = a2a.load_reference(conditioned_path, role="conditioned-crest")
    neb_crest = a2a.load_reference(
        neb_crest_path,
        role="neb-crest",
        template=conditioned,
    )
    local_indices = a2a.active_indices(args.local_indices, len(conditioned.symbols))
    if any(index in conditioned.frozen_indices for index in local_indices):
        raise ValueError("local-indices contains a frozen atom")
    mode, neighborhood = localization_inputs(
        checkpoint_root,
        neb_crest,
        conditioned,
        local_indices,
    )
    r2scan3c, _, _ = a2a.a2.settings(use_gpu=args.gpu)
    identity = {
        "localization_version": LOCALIZATION_VERSION,
        "conditioned_crest_sha256": a2.sha256_path(conditioned_path),
        "conditioned_crest_geometry": frequency_geometry_fingerprint(conditioned),
        "neb_crest_sha256": a2.sha256_path(neb_crest_path),
        "settings": frequency_settings_fingerprint(r2scan3c),
        "neighborhood": neighborhood,
        "parameters": {
            "local_indices": local_indices,
            "fmax_ev_a": args.fmax,
            "max_steps": args.max_steps,
            "maximum_translation_a": args.maximum_translation,
            "dimer_separation_a": args.dimer_separation,
            "local_restraint_k_ev_a2": args.restraint_k,
            "f_rot_min_ev_a": args.f_rot_min,
            "f_rot_max_ev_a": args.f_rot_max,
            "max_num_rot": args.max_num_rot,
        },
        "production_single_points_run": False,
        "finite_difference_hessian_run": False,
        "full_irc_run": False,
    }
    candidate_path = run_dir / "full-system-local-dimer-transition-state.xyz"
    algorithm_receipt: dict[str, Any] = {}

    def compute():
        candidate, receipt = find_ts_dimer(
            conditioned,
            r2scan3c,
            initial_mode=mode,
            local_indices=local_indices,
            local_trust_radius_a=neighborhood["local_trust_radius_a"],
            local_guard_radius_a=neighborhood["local_guard_radius_a"],
            local_restraint_k_ev_a2=args.restraint_k,
            fmax_ev_a=args.fmax,
            max_steps=args.max_steps,
            maximum_translation_a=args.maximum_translation,
            dimer_separation_a=args.dimer_separation,
            f_rot_min_ev_a=args.f_rot_min,
            f_rot_max_ev_a=args.f_rot_max,
            max_num_rot=args.max_num_rot,
            trajectory=str(run_dir / "full-system-local-dimer.traj"),
            logfile=str(run_dir / "full-system-local-dimer.log"),
            eigenmode_logfile=str(run_dir / "full-system-local-dimer-mode.log"),
        )
        algorithm_receipt.update(receipt)
        return candidate

    candidate = a2.checkpoint_cluster(
        candidate_path,
        conditioned,
        compute,
        identity=identity,
    )
    receipt_path = run_dir / "full-system-local-dimer.receipt.json"
    if not algorithm_receipt:
        if not receipt_path.is_file():
            raise ValueError("cached dimer candidate lacks its algorithm receipt")
        cached = json.loads(receipt_path.read_text())
        if cached.get("identity") != identity:
            raise ValueError("cached dimer algorithm receipt identity drift")
        algorithm_receipt = dict(cached["algorithm"])
    receipt = {
        "identity": identity,
        "algorithm": algorithm_receipt,
        "candidate_path": str(candidate_path),
        "candidate_sha256": a2.sha256_path(candidate_path),
        "candidate_geometry": frequency_geometry_fingerprint(candidate),
        "candidate_persisted_only_after_physical_convergence": True,
    }
    atomic_json(receipt_path, receipt)
    atomic_json(
        run_dir / "localization_status.json",
        {
            "status": "completed",
            "localization_version": LOCALIZATION_VERSION,
            "completed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "candidate_sha256": receipt["candidate_sha256"],
            "receipt_sha256": a2.sha256_path(receipt_path),
        },
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


def execute_with_status(args: argparse.Namespace) -> int:
    status_path = args.run_dir.resolve() / "localization_status.json"
    atomic_json(
        status_path,
        {
            "status": "running",
            "localization_version": LOCALIZATION_VERSION,
            "started_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        },
    )
    try:
        return run(args)
    except Exception as exc:
        atomic_json(
            status_path,
            {
                "status": "failed",
                "localization_version": LOCALIZATION_VERSION,
                "error_type": type(exc).__name__,
                "error": str(exc),
                "failed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            },
        )
        raise


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--conditioned-crest", type=Path, required=True)
    result.add_argument("--neb-crest", type=Path, required=True)
    result.add_argument("--checkpoint-root", type=Path, required=True)
    result.add_argument("--run-dir", type=Path, required=True)
    result.add_argument("--local-indices", default="0,1,2,15,16,17")
    result.add_argument("--gpu", action="store_true")
    result.add_argument("--fmax", type=float, default=0.02)
    result.add_argument("--max-steps", type=int, default=120)
    result.add_argument("--maximum-translation", type=float, default=0.03)
    result.add_argument("--dimer-separation", type=float, default=0.01)
    result.add_argument("--restraint-k", type=float, default=50.0)
    result.add_argument("--f-rot-min", type=float, default=0.01)
    result.add_argument("--f-rot-max", type=float, default=0.05)
    result.add_argument("--max-num-rot", type=int, default=4)
    result.add_argument("--threads", type=int, default=16)
    result.add_argument("--nice", type=int, default=10)
    result.add_argument("--log")
    return result


def main() -> int:
    args = parser().parse_args()
    return execute_with_status(args)


if __name__ == "__main__":
    raise SystemExit(main())
