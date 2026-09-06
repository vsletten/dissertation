#!/usr/bin/env python3
"""Cold optimizer-free verifier for the bounded A3e H35 owner-basin run."""

from __future__ import annotations

import argparse
import ast
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    bootstrap_cli(
        "a3e_verify",
        default_run_root=Path(
            "/mnt/data/vsletten/dissertation-data/"
            "a3e-oss-neutral-n1-proton-microstate-stability/logs"
        ),
        gpu_owner="a3e_verify",
    )

from quarry.pipeline import (  # noqa: E402
    energy,
    frequencies,
    frequency_geometry_fingerprint,
    frequency_settings_fingerprint,
)
from scripts import a3b_proton_microstate_stability as a3b  # noqa: E402
from scripts import a3e_oss_owner_basin as a3e  # noqa: E402


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def _read_stage(
    output_root: Path,
    source: a3e.SourceEvidence,
    spec: a3e.StageSpec,
    code_revision: str,
    parent: dict[str, Any],
    *,
    recompute_calculators: bool,
) -> tuple[dict[str, Any], Any | None]:
    stage_dir = output_root / "stages" / spec.directory
    reservation_path = stage_dir / "reservation.json"
    receipt_path = stage_dir / "receipt.json"
    reservation = json.loads(reservation_path.read_text())
    receipt = json.loads(receipt_path.read_text())
    signature = a3e.stage_signature(source, spec, code_revision)
    if (
        reservation.get("schema") != "a3e-stage-reservation-v1"
        or reservation.get("status") != "reserved"
        or reservation.get("signature") != signature
        or reservation.get("parent") != parent
        or receipt.get("schema") != "a3e-stage-receipt-v1"
        or receipt.get("signature") != signature
        or receipt.get("parent") != parent
    ):
        raise RuntimeError(f"{spec.stage_id} reservation/receipt identity mismatch")
    if not (
        receipt.get("optimizer", {}).get("fresh_instance") is True
        and receipt.get("optimizer", {}).get("geometric_default_fresh_hessian") is True
    ):
        raise RuntimeError(f"{spec.stage_id} optimizer freshness evidence mismatch")
    optimizer = receipt.get("optimizer") or {}
    if not (
        optimizer.get("observed_calls") == 1
        and optimizer.get("observed_retries") == 0
        and optimizer.get("observed_max_steps") == a3e.MAX_STEPS
    ):
        raise RuntimeError(f"{spec.stage_id} observed optimizer budget mismatch")
    raw_path = stage_dir / "raw-endpoint.xyz"
    if receipt.get("status") == "optimizer-failed":
        if raw_path.exists() or receipt.get("raw_endpoint") is not None:
            raise RuntimeError(
                f"{spec.stage_id} optimizer failure has inconsistent raw evidence"
            )
        return receipt, None
    if not raw_path.is_file() or receipt.get("raw_endpoint", {}).get(
        "sha256"
    ) != a3e.sha256_path(raw_path):
        raise RuntimeError(f"{spec.stage_id} raw endpoint missing or hash-mismatched")
    raw_cluster = a3b.read_cluster(raw_path, source.cluster)
    if receipt.get("raw_endpoint") != a3e._artifact(raw_cluster, raw_path):
        raise RuntimeError(f"{spec.stage_id} raw endpoint artifact mismatch")
    if receipt.get("status") != "complete":
        return receipt, None
    endpoint_path = stage_dir / "endpoint.xyz"
    if not endpoint_path.is_file() or receipt.get("endpoint", {}).get(
        "sha256"
    ) != a3e.sha256_path(endpoint_path):
        raise RuntimeError(f"{spec.stage_id} endpoint missing or hash-mismatched")
    endpoint = a3b.read_cluster(endpoint_path, source.cluster)
    if receipt.get("endpoint", {}).get(
        "geometry_fingerprint"
    ) != frequency_geometry_fingerprint(endpoint):
        raise RuntimeError(f"{spec.stage_id} endpoint geometry fingerprint mismatch")
    active = a3e.constraints(source, spec)
    raw_structure = a3b.structural_gate(
        raw_cluster, source.cluster, active, spec.stage_id
    )[1]
    structure = a3b.structural_gate(endpoint, source.cluster, active, spec.stage_id)[1]
    owners = a3e._owner_labels(endpoint)
    if receipt.get("observed_owners") != owners:
        raise RuntimeError(f"{spec.stage_id} owner receipt mismatch")
    if receipt.get("structure", {}).get("constraint_residuals") != structure.get(
        "constraint_residuals"
    ):
        raise RuntimeError(f"{spec.stage_id} constraint residual receipt mismatch")
    if receipt.get("structure", {}).get("owner_changes") != structure.get(
        "owner_changes"
    ):
        raise RuntimeError(f"{spec.stage_id} owner-change receipt mismatch")
    if receipt.get("structure", {}).get(
        "raw_maximum_frozen_coordinate_drift_a"
    ) != raw_structure.get("raw_maximum_frozen_coordinate_drift_a"):
        raise RuntimeError(f"{spec.stage_id} raw frozen-shell receipt mismatch")
    verified = dict(receipt)
    verified["structure"] = {
        **structure,
        "raw_maximum_frozen_coordinate_drift_a": raw_structure[
            "raw_maximum_frozen_coordinate_drift_a"
        ],
    }
    verified["owner_changes"] = structure["owner_changes"]
    verified["owner_retaining"] = structure["owner_changes"] == []
    if recompute_calculators:
        observed_energy = float(energy(endpoint, spec.settings))
        if not math.isfinite(observed_energy):
            raise RuntimeError(f"{spec.stage_id} recomputed energy is non-finite")
        metrics = a3b.gradient_metrics(endpoint, spec.settings, active)
        verified["evidence"] = {"energy_hartree": observed_energy, **metrics}
        verified["stationary"] = bool(
            receipt.get("optimizer", {}).get("converged") is True
            and metrics["passed"] is True
        )
        if (
            spec.fully_released
            and verified["owner_retaining"]
            and verified["stationary"]
        ):
            result = frequencies(endpoint, spec.settings)
            imaginary = [float(value) for value in result.imaginary_cm]
            if not math.isfinite(float(result.electronic_hartree)) or not all(
                math.isfinite(value) for value in imaginary
            ):
                raise RuntimeError("released PHVA recomputation is non-finite")
            verified["phva"] = {
                "status": (
                    "passed"
                    if not any(value > a3e.NOISE_FLOOR_CM for value in imaginary)
                    else "failed"
                ),
                "imaginary_cm": imaginary,
                "noise_floor_cm": a3e.NOISE_FLOOR_CM,
                "fresh_hessian": True,
                "geometry_fingerprint": frequency_geometry_fingerprint(endpoint),
                "settings_fingerprint": frequency_settings_fingerprint(spec.settings),
                "electronic_hartree": float(result.electronic_hartree),
            }
    return verified, endpoint


def verify_experiment(
    output_root: Path,
    *,
    verifier_identity: str,
    source_root: Path = a3e.DEFAULT_SOURCE_ROOT,
    source_override: a3e.SourceEvidence | None = None,
    recompute_calculators: bool = True,
) -> dict[str, Any]:
    """Rehash source and independently derive the A3e terminal classification."""
    output_root = output_root.resolve()
    candidate: dict[str, Any] = {}
    with a3e.exclusive_run(output_root):
        try:
            candidate_path = output_root / "candidate-terminal.json"
            candidate = json.loads(candidate_path.read_text())
            if candidate.get("schema") != "a3e-candidate-terminal-v1":
                raise RuntimeError("candidate terminal schema mismatch")
            executor_identity = candidate.get("executor_identity")
            if (
                not isinstance(executor_identity, str)
                or not verifier_identity
                or verifier_identity == executor_identity
            ):
                raise RuntimeError(
                    "verifier identity must differ from executor identity"
                )
            if not recompute_calculators:
                raise RuntimeError("verification requires calculator recomputation")
            source = source_override or a3e.validate_source(source_root)
            if candidate.get("source") != a3e.source_map(source):
                raise RuntimeError("candidate source binding mismatch")
            code_revision = candidate.get("code_revision")
            if not isinstance(code_revision, str) or len(code_revision) != 40:
                raise RuntimeError("candidate code revision is invalid")

            receipts: dict[str, dict[str, Any]] = {}
            endpoints: dict[str, Any | None] = {}
            parent = {"stage": None, "receipt_sha256": None}
            for spec in a3e.STAGES:
                stage_state = candidate.get("stages", {}).get(spec.stage_id, {})
                if stage_state.get("status") == "not-run":
                    stage_dir = output_root / "stages" / spec.directory
                    if stage_dir.exists():
                        raise RuntimeError(
                            f"{spec.stage_id} is marked not-run but has artifacts"
                        )
                    receipts[spec.stage_id] = {
                        "schema": "a3e-stage-receipt-v1",
                        "status": "not-run",
                        "stage": spec.stage_id,
                    }
                    endpoints[spec.stage_id] = None
                else:
                    receipt_path = (
                        output_root / "stages" / spec.directory / "receipt.json"
                    )
                    if stage_state.get("receipt_sha256") != a3e.sha256_path(
                        receipt_path
                    ):
                        raise RuntimeError(
                            f"{spec.stage_id} candidate receipt hash mismatch"
                        )
                    receipt, endpoint = _read_stage(
                        output_root,
                        source,
                        spec,
                        code_revision,
                        parent,
                        recompute_calculators=recompute_calculators,
                    )
                    receipts[spec.stage_id] = receipt
                    endpoints[spec.stage_id] = endpoint
                if spec.parent is None:
                    parent = {
                        "stage": spec.stage_id,
                        "receipt_sha256": (
                            a3e.sha256_path(
                                output_root / "stages" / spec.directory / "receipt.json"
                            )
                            if (
                                output_root / "stages" / spec.directory / "receipt.json"
                            ).is_file()
                            else None
                        ),
                    }
                elif (
                    output_root / "stages" / spec.directory / "receipt.json"
                ).is_file():
                    parent = {
                        "stage": spec.stage_id,
                        "receipt_sha256": a3e.sha256_path(
                            output_root / "stages" / spec.directory / "receipt.json"
                        ),
                    }

            conditioning = receipts[a3e.STAGES[0].stage_id]
            constrained = receipts[a3e.STAGES[1].stage_id]
            released = receipts[a3e.STAGES[2].stage_id]
            if candidate.get("experiment_budget") != a3e._experiment_budget(receipts):
                raise RuntimeError("candidate experiment_budget mismatch")
            if constrained.get("status") != "not-run" and not (
                conditioning.get("status") == "complete"
                and conditioning.get("owner_retaining") is True
                and conditioning.get("owner_changes") == []
                and a3e._finite_energy(conditioning) is not None
            ):
                raise RuntimeError(
                    "constrained stage ran without valid conditioning parent"
                )
            if released.get(
                "status"
            ) != "not-run" and not a3e.constrained_release_allowed(constrained):
                raise RuntimeError(
                    "released stage ran without a valid constrained release gate"
                )
            classification = a3e.classify(constrained, released)
            if classification != candidate.get("classification"):
                raise RuntimeError(
                    "candidate classification mismatch: "
                    f"claimed {candidate.get('classification')!r}, "
                    f"derived {classification!r}"
                )
            verified_classification = {
                a3e.ACCEPTED_CANDIDATE: a3e.VERIFIED_ACCEPTED,
                a3e.NO_BASIN_CANDIDATE: a3e.VERIFIED_NO_BASIN,
                a3e.INCONCLUSIVE_CANDIDATE: a3e.VERIFIED_INCONCLUSIVE,
            }[classification]
            result = {
                "schema": "a3e-verified-terminal-v1",
                "written_at": now(),
                "status": "verified",
                "classification": verified_classification,
                "executor_identity": executor_identity,
                "verifier_identity": verifier_identity,
                "candidate_terminal_sha256": a3e.sha256_path(candidate_path),
                "source_rehashed": source_override is None,
                "calculator_evidence_recomputed": True,
                "code_revision": code_revision,
                "stage_receipt_sha256": {
                    spec.stage_id: candidate.get("stages", {})
                    .get(spec.stage_id, {})
                    .get("receipt_sha256")
                    for spec in a3e.STAGES
                },
            }
        except Exception as exc:
            result = {
                "schema": "a3e-verified-terminal-v1",
                "written_at": now(),
                "status": "rejected",
                "classification": a3e.VERIFIED_INCONCLUSIVE,
                "executor_identity": candidate.get("executor_identity"),
                "verifier_identity": verifier_identity,
                "detail": f"{type(exc).__name__}: {exc}",
            }
        a3e.atomic_json(output_root / "verified-terminal.json", result)
        return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=a3e.DEFAULT_SOURCE_ROOT)
    parser.add_argument("--output-root", type=Path, default=a3e.DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--verifier-identity", required=True)
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--gpu-mem-gb", type=float, default=16.0)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--nice", type=int, default=10)
    parser.add_argument("--log", type=Path)
    args = parser.parse_args()
    if args.threads < 1 or args.threads > 16:
        parser.error("--threads must be within 1..16")
    if args.nice < 10:
        parser.error("--nice must be >=10")
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.environ["MKL_NUM_THREADS"] = str(args.threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(args.threads)
    candidate = json.loads((args.output_root / "candidate-terminal.json").read_text())
    production_ran = (
        candidate.get("stages", {}).get("constrained-production", {}).get("status")
        != "not-run"
    )
    if production_ran and not args.gpu:
        parser.error("--gpu is required to recompute production calculator evidence")
    if args.gpu:
        a3e.preload_cutensor()
    result = verify_experiment(
        args.output_root,
        source_root=args.source_root,
        verifier_identity=args.verifier_identity,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "verified" else 1


if __name__ == "__main__":
    tree = ast.parse(Path(__file__).read_text())
    raise SystemExit(main())
