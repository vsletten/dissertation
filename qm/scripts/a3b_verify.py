#!/usr/bin/env python3
"""Independent, optimizer-free verifier for the bounded A3b experiment."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

DEFAULT_SOURCE_ROOT = Path(
    "/mnt/data/vsletten/dissertation-data/a3a-reactant-minimum-recovery"
)
DEFAULT_OUTPUT_ROOT = Path(
    "/mnt/data/vsletten/dissertation-data/a3b-osa-neutral-n1-proton-microstate-stability"
)

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    bootstrap_cli(
        "a3b_verify",
        default_run_root=DEFAULT_OUTPUT_ROOT / "logs",
        gpu_owner="a3b_verify",
    )

import numpy as np  # noqa: E402

from quarry.clusters import Cluster  # noqa: E402
from quarry.pipeline import (  # noqa: E402
    DftSettings,
    energy,
    frequencies,
    frequency_geometry_fingerprint,
    gradient,
)
from quarry.store import geometry_hash  # noqa: E402
from scripts import a3b_proton_microstate_stability as a3b  # noqa: E402
from scripts.phase2_ladder import oxygen_proton_owners  # noqa: E402

DEFAULT_VERIFIER_IDENTITY = "hermes-a3b-independent-verifier"


@dataclass(frozen=True)
class StaticSource:
    cluster: Cluster
    source_root: Path
    closeout_manifest_sha256: str
    conditioned_xyz_sha256: str
    conditioned_receipt_sha256: str
    conditioned_geometry_hash: str

    @classmethod
    def from_executor_source(cls, source: a3b.SourceEvidence) -> StaticSource:
        return cls(
            cluster=source.cluster,
            source_root=source.source_root,
            closeout_manifest_sha256=source.closeout_manifest_sha256,
            conditioned_xyz_sha256=source.conditioned_xyz_sha256,
            conditioned_receipt_sha256=source.conditioned_receipt_sha256,
            conditioned_geometry_hash=source.conditioned_geometry_hash,
        )


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def _load_source(source_root: Path) -> StaticSource:
    """Independently rehash fixed A3a paths before reading the candidate."""
    source_root = source_root.resolve()
    manifest_path = source_root / "closeout-manifest.json"
    xyz_path = source_root / a3b.SOURCE_RUN_RELATIVE / "complex_conditioned.xyz"
    receipt_path = source_root / a3b.SOURCE_RUN_RELATIVE / "complex_conditioned.json"
    observed = {
        "manifest": a3b.sha256_path(manifest_path),
        "xyz": a3b.sha256_path(xyz_path),
        "receipt": a3b.sha256_path(receipt_path),
    }
    expected = {
        "manifest": a3b.EXPECTED_CLOSEOUT_MANIFEST_SHA256,
        "xyz": a3b.EXPECTED_CONDITIONED_XYZ_SHA256,
        "receipt": a3b.EXPECTED_CONDITIONED_RECEIPT_SHA256,
    }
    for key in expected:
        if observed[key] != expected[key]:
            raise RuntimeError(f"source {key} SHA-256 mismatch: {observed[key]}")
    manifest = json.loads(manifest_path.read_text())
    receipt = json.loads(receipt_path.read_text())
    if (
        manifest.get("schema") != "a3a-closeout-manifest-v1"
        or manifest.get("status") != "failed"
    ):
        raise RuntimeError("source closeout identity/status mismatch")
    if receipt.get("schema") != "phase2-reactant-conditioning-v1":
        raise RuntimeError("source conditioned receipt schema mismatch")
    signature = receipt.get("signature")
    if not isinstance(signature, dict):
        raise RuntimeError("source conditioned receipt omitted signature")
    symbols = signature.get("symbols")
    frozen = signature.get("frozen_indices")
    if not isinstance(symbols, list) or not isinstance(frozen, list):
        raise RuntimeError("source conditioned identity is malformed")
    lines = xyz_path.read_text().splitlines()
    template = Cluster(
        lines[1],
        list(symbols),
        np.zeros((len(symbols), 3)),
        charge=int(signature.get("charge", 0)),
        spin=int(signature.get("spin", 0)),
        frozen_indices=list(frozen),
    )
    cluster = a3b.read_cluster(xyz_path, template)
    fingerprint = geometry_hash(cluster.to_xyz())
    if receipt.get("endpoint_geometry_hash") != fingerprint:
        raise RuntimeError("source conditioned geometry hash mismatch")
    owners = oxygen_proton_owners(cluster)
    if owners.get(50) != 26 or owners.get(57) != 29:
        raise RuntimeError("source required owner assignment changed")
    if signature.get("oxygen_proton_owners") != [
        f"H{h}:O{o}" for h, o in sorted(owners.items())
    ]:
        raise RuntimeError("source complete owner map changed")
    return StaticSource(
        cluster=cluster,
        source_root=source_root,
        closeout_manifest_sha256=observed["manifest"],
        conditioned_xyz_sha256=observed["xyz"],
        conditioned_receipt_sha256=observed["receipt"],
        conditioned_geometry_hash=fingerprint,
    )


def _constraints(
    source: StaticSource, spec: a3b.StageSpec
) -> list[tuple[int, int, float]]:
    return [
        (
            oxygen,
            hydrogen,
            float(
                np.linalg.norm(
                    source.cluster.coords[oxygen] - source.cluster.coords[hydrogen]
                )
            ),
        )
        for oxygen, hydrogen in spec.constraint_pairs
    ]


def _structure(
    endpoint: Cluster,
    source: StaticSource,
    constraints: list[tuple[int, int, float]],
    stage: str,
) -> dict[str, Any]:
    reference = source.cluster
    if endpoint.symbols != reference.symbols:
        raise RuntimeError(f"{stage} endpoint identity/order changed")
    if endpoint.charge != reference.charge or endpoint.spin != reference.spin:
        raise RuntimeError(f"{stage} endpoint charge/spin changed")
    if endpoint.frozen_indices != reference.frozen_indices:
        raise RuntimeError(f"{stage} endpoint frozen set changed")
    if endpoint.coords.shape != reference.coords.shape or not np.all(
        np.isfinite(endpoint.coords)
    ):
        raise RuntimeError(f"{stage} endpoint coordinates are invalid")
    frozen = sorted(reference.frozen_indices)
    frozen_drift = (
        float(np.max(np.abs(endpoint.coords[frozen] - reference.coords[frozen])))
        if frozen
        else 0.0
    )
    if frozen_drift > a3b.ADVISORY_PREOPT_MAX_PROJECTED_FROZEN_DRIFT_A:
        raise RuntimeError(f"{stage} endpoint frozen shell drifted")
    delta = endpoint.coords[:, None, :] - endpoint.coords[None, :, :]
    distances = np.linalg.norm(delta, axis=2)
    distances[np.diag_indices_from(distances)] = np.inf
    minimum_pair = float(np.min(distances))
    if minimum_pair <= a3b.ADVISORY_PREOPT_MIN_PAIR_A:
        raise RuntimeError(f"{stage} endpoint collision")
    owners = oxygen_proton_owners(endpoint)
    reference_owners = oxygen_proton_owners(reference)
    residuals = []
    for oxygen, hydrogen, target in constraints:
        observed = float(
            np.linalg.norm(endpoint.coords[oxygen] - endpoint.coords[hydrogen])
        )
        residual = abs(observed - target)
        if residual > a3b.CONSTRAINT_RESIDUAL_MAX_A:
            raise RuntimeError(f"{stage} active constraint residual exceeds gate")
        residuals.append(residual)
    return {
        "minimum_pair_distance_a": minimum_pair,
        "maximum_frozen_coordinate_drift_a": frozen_drift,
        "constraint_residuals_a": residuals,
        "owners": [f"H{h}:O{o}" for h, o in sorted(owners.items())],
        "owner_retaining": owners == reference_owners,
        "owner_changes": [
            f"H{hydrogen}:O{reference_owners.get(hydrogen)}->O{owners.get(hydrogen)}"
            for hydrogen in sorted(reference_owners.keys() | owners.keys())
            if reference_owners.get(hydrogen) != owners.get(hydrogen)
        ],
    }


def _project_gradient(
    endpoint: Cluster,
    values: np.ndarray,
    constraints: list[tuple[int, int, float]],
) -> np.ndarray:
    projected = np.asarray(values, dtype=float).copy()
    if projected.shape != endpoint.coords.shape or not np.all(np.isfinite(projected)):
        raise RuntimeError("recomputed gradient is malformed or non-finite")
    frozen = set(endpoint.frozen_indices)
    for index in frozen:
        projected[index] = 0.0
    for oxygen, hydrogen, _target in constraints:
        vector = endpoint.coords[hydrogen] - endpoint.coords[oxygen]
        unit = vector / float(np.linalg.norm(vector))
        normal = np.zeros_like(projected)
        if oxygen not in frozen:
            normal[oxygen] = -unit
        if hydrogen not in frozen:
            normal[hydrogen] = unit
        norm2 = float(np.sum(normal * normal))
        if norm2:
            projected -= float(np.sum(projected * normal)) / norm2 * normal
    return projected


def _gradient_metrics(
    endpoint: Cluster,
    settings: DftSettings,
    constraints: list[tuple[int, int, float]],
) -> dict[str, Any]:
    projected = _project_gradient(endpoint, gradient(endpoint, settings), constraints)
    free = sorted(set(range(len(endpoint.symbols))) - set(endpoint.frozen_indices))
    values = projected[free]
    rms = float(np.sqrt(np.mean(values**2))) if values.size else 0.0
    maximum = float(np.max(np.abs(values))) if values.size else 0.0
    return {
        "rms": rms,
        "maximum": maximum,
        "passed": (
            math.isfinite(rms)
            and math.isfinite(maximum)
            and rms <= a3b.GRADIENT_RMS_MAX_HARTREE_PER_BOHR
            and maximum <= a3b.GRADIENT_MAX_HARTREE_PER_BOHR
        ),
    }


def _finite_energy(receipt: dict[str, Any]) -> float | None:
    value = receipt.get("evidence", {}).get("energy_hartree")
    if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
        return None
    return float(value)


def _is_recovery(receipt: dict[str, Any]) -> bool:
    imaginary = receipt.get("phva", {}).get("imaginary_cm", [])
    return bool(
        receipt.get("status") == "complete"
        and receipt.get("fully_unconstrained") is True
        and receipt.get("owner_retaining") is True
        and receipt.get("stationary") is True
        and receipt.get("phva", {}).get("status") == "passed"
        and isinstance(imaginary, list)
        and all(
            isinstance(value, (int, float)) and math.isfinite(float(value))
            for value in imaginary
        )
        and not any(float(value) > a3b.NOISE_FLOOR_CM for value in imaginary)
        and _finite_energy(receipt) is not None
    )


def _classify(
    common: dict[str, Any], final_a: dict[str, Any], final_b: dict[str, Any]
) -> str:
    if _is_recovery(final_a) or _is_recovery(final_b):
        return a3b.RECOVERY_CANDIDATE
    common_energy = _finite_energy(common)
    energies = [_finite_energy(final_a), _finite_energy(final_b)]
    if common_energy is None or any(value is None for value in energies):
        return a3b.INCONCLUSIVE
    if all(
        receipt.get("status") == "complete"
        and receipt.get("fully_unconstrained") is True
        and receipt.get("owner_retaining") is False
        and receipt.get("stationary") is True
        for receipt in (final_a, final_b)
    ) and all(
        common_energy - float(value) > a3b.DOWNHILL_MIN_HARTREE
        for value in energies
        if value is not None
    ):
        return a3b.NO_BASIN_CANDIDATE
    return a3b.INCONCLUSIVE


def _verified_receipt(
    output_root: Path,
    source: StaticSource,
    spec: a3b.StageSpec,
    candidate_settings: DftSettings,
    *,
    recompute_calculators: bool,
) -> dict[str, Any]:
    stage_dir = output_root / "stages" / spec.directory
    receipt_path = stage_dir / "receipt.json"
    if not receipt_path.is_file():
        if spec.parent is None:
            raise RuntimeError(f"{spec.stage_id} receipt is required")
        return {"status": "not-run", "stage": spec.stage_id}
    receipt = json.loads(receipt_path.read_text())
    if receipt.get("stage") != spec.stage_id:
        raise RuntimeError(f"{spec.stage_id} stage identity mismatch")
    expected_source = {
        "closeout_manifest_sha256": source.closeout_manifest_sha256,
        "conditioned_xyz_sha256": source.conditioned_xyz_sha256,
        "conditioned_receipt_sha256": source.conditioned_receipt_sha256,
        "conditioned_geometry_hash": source.conditioned_geometry_hash,
    }
    expected_constraints = _constraints(source, spec)
    expected_signature = {
        "schema": "a3b-stage-signature-v1",
        "stage": spec.stage_id,
        "route": spec.route,
        "ordinal": spec.ordinal,
        "parent": spec.parent,
        "source": expected_source,
        "settings": asdict(candidate_settings),
        "constraints": [list(item) for item in expected_constraints],
        "released": spec.released,
        "budget": {"max_steps": 100, "continuation_allowed": False},
        "fresh_optimizer": True,
        "geometric_default_fresh_hessian": True,
    }
    if receipt.get("signature") != expected_signature:
        raise RuntimeError(f"{spec.stage_id} signature mismatch")
    if receipt.get("schema") != "a3b-stage-receipt-v1":
        raise RuntimeError(f"{spec.stage_id} receipt schema mismatch")
    if receipt.get("source") != expected_source:
        raise RuntimeError(f"{spec.stage_id} source binding mismatch")
    if receipt.get("budget") != {"max_steps": 100, "continuation_allowed": False}:
        raise RuntimeError(f"{spec.stage_id} budget mismatch")
    optimizer = receipt.get("optimizer", {})
    if not (
        optimizer.get("fresh_instance") is True
        and optimizer.get("geometric_default_fresh_hessian") is True
    ):
        raise RuntimeError(f"{spec.stage_id} optimizer freshness mismatch")
    expected_parent = {
        "stage": spec.parent,
        "receipt_sha256": (
            a3b.sha256_path(
                output_root
                / "stages"
                / a3b.STAGE_BY_ID[spec.parent].directory
                / "receipt.json"
            )
            if spec.parent
            else None
        ),
    }
    if receipt.get("parent") != expected_parent:
        raise RuntimeError(f"{spec.stage_id} parent DAG mismatch")
    if receipt.get("settings") != asdict(candidate_settings):
        raise RuntimeError(f"{spec.stage_id} settings mismatch")
    if receipt.get("constraints") != [list(item) for item in expected_constraints]:
        raise RuntimeError(f"{spec.stage_id} constraint set mismatch")
    seed_cluster = source.cluster
    if spec.parent:
        seed_cluster = a3b.read_cluster(
            output_root
            / "stages"
            / a3b.STAGE_BY_ID[spec.parent].directory
            / "endpoint.xyz",
            source.cluster,
        )
    expected_seed_fingerprint = frequency_geometry_fingerprint(seed_cluster)
    if receipt.get("seed", {}).get("geometry_fingerprint") != expected_seed_fingerprint:
        raise RuntimeError(f"{spec.stage_id} seed geometry mismatch")
    reservation_path = stage_dir / "reservation.json"
    reservation = json.loads(reservation_path.read_text())
    if (
        reservation.get("schema") != "a3b-stage-reservation-v1"
        or reservation.get("status") != "reserved"
        or reservation.get("signature") != expected_signature
        or reservation.get("parent") != expected_parent
        or reservation.get("seed", {}).get("geometry_fingerprint")
        != expected_seed_fingerprint
        or reservation.get("seed", {}).get("sha256")
        != hashlib.sha256(seed_cluster.to_xyz().encode()).hexdigest()
    ):
        raise RuntimeError(f"{spec.stage_id} reservation mismatch")
    raw = receipt.get("raw_endpoint")
    if receipt.get("status") == "complete" and not isinstance(raw, dict):
        raise RuntimeError(f"{spec.stage_id} raw endpoint metadata missing")
    if isinstance(raw, dict):
        raw_path = stage_dir / "raw-endpoint.xyz"
        if not raw_path.is_file():
            raise RuntimeError(f"{spec.stage_id} raw endpoint artifact missing")
        if raw.get("sha256") != a3b.sha256_path(raw_path):
            raise RuntimeError(f"{spec.stage_id} raw endpoint byte hash mismatch")
    if receipt.get("status") != "complete":
        return receipt
    endpoint_path = stage_dir / "endpoint.xyz"
    if receipt.get("endpoint", {}).get("sha256") != a3b.sha256_path(endpoint_path):
        raise RuntimeError(f"{spec.stage_id} endpoint byte hash mismatch")
    endpoint = a3b.read_cluster(endpoint_path, source.cluster)
    if receipt.get("endpoint", {}).get(
        "geometry_fingerprint"
    ) != frequency_geometry_fingerprint(endpoint):
        raise RuntimeError(f"{spec.stage_id} endpoint geometry fingerprint mismatch")
    structure = _structure(endpoint, source, expected_constraints, spec.stage_id)
    if receipt.get("observed_owners") != structure["owners"]:
        raise RuntimeError(f"{spec.stage_id} owner receipt mismatch")
    if receipt.get("owner_changes") != structure["owner_changes"]:
        raise RuntimeError(f"{spec.stage_id} owner-change receipt mismatch")

    verified = dict(receipt)
    verified["owner_retaining"] = structure["owner_retaining"]
    if recompute_calculators:
        observed_energy = float(energy(endpoint, candidate_settings))
        if not math.isfinite(observed_energy):
            raise RuntimeError(f"{spec.stage_id} recomputed non-finite energy")
        metrics = _gradient_metrics(endpoint, candidate_settings, expected_constraints)
        verified["evidence"] = {
            "energy_hartree": observed_energy,
            "projected_gradient_rms_hartree_per_bohr": metrics["rms"],
            "projected_gradient_max_hartree_per_bohr": metrics["maximum"],
            "passed": metrics["passed"],
        }
        verified["stationary"] = bool(
            receipt.get("optimizer", {}).get("converged")
        ) and bool(metrics["passed"])
        verified["phva"] = {"status": "not-required", "imaginary_cm": []}
        if (
            spec.fully_unconstrained
            and verified["owner_retaining"]
            and verified["stationary"]
        ):
            frequency = frequencies(endpoint, candidate_settings)
            imaginary = [float(value) for value in frequency.imaginary_cm]
            if not math.isfinite(float(frequency.electronic_hartree)) or not all(
                math.isfinite(value) for value in imaginary
            ):
                raise RuntimeError(f"{spec.stage_id} recomputed PHVA is non-finite")
            verified["phva"] = {
                "status": (
                    "passed"
                    if not any(value > a3b.NOISE_FLOOR_CM for value in imaginary)
                    else "failed"
                ),
                "imaginary_cm": imaginary,
                "noise_floor_cm": a3b.NOISE_FLOOR_CM,
                "fresh_hessian": True,
            }
    return verified


def verify_experiment(
    output_root: Path,
    *,
    source_root: Path = DEFAULT_SOURCE_ROOT,
    verifier_identity: str = DEFAULT_VERIFIER_IDENTITY,
    source_override: StaticSource | None = None,
) -> dict[str, Any]:
    """Re-derive the candidate from static artifacts; never run an optimizer."""
    output_root = output_root.resolve()
    candidate_path = output_root / "candidate-terminal.json"
    candidate: dict[str, Any] = {}
    try:
        candidate = json.loads(candidate_path.read_text())
        if candidate.get("schema") != "a3b-candidate-terminal-v1":
            raise RuntimeError("candidate schema mismatch")
        executor_identity = candidate.get("executor_identity")
        if (
            not isinstance(executor_identity, str)
            or executor_identity == verifier_identity
        ):
            raise RuntimeError("verifier identity must differ from executor identity")
        source = source_override or _load_source(source_root)
        expected_source = {
            "closeout_manifest_sha256": source.closeout_manifest_sha256,
            "conditioned_xyz_sha256": source.conditioned_xyz_sha256,
            "conditioned_receipt_sha256": source.conditioned_receipt_sha256,
            "conditioned_geometry_hash": source.conditioned_geometry_hash,
        }
        if candidate.get("source") != expected_source:
            raise RuntimeError("candidate source binding mismatch")
        settings_data = candidate.get("settings")
        if not isinstance(settings_data, dict):
            raise RuntimeError("candidate settings missing")
        settings = DftSettings(**settings_data)
        receipts = {
            spec.stage_id: _verified_receipt(
                output_root,
                source,
                spec,
                settings,
                recompute_calculators=source_override is None,
            )
            for spec in a3b.STAGES
        }
        classification = _classify(
            receipts["common-dual"],
            receipts["h50-then-h57-unconstrained"],
            receipts["h57-then-h50-unconstrained"],
        )
        if candidate.get("classification") != classification:
            claimed = candidate.get("classification")
            raise RuntimeError(
                "candidate classification mismatch: "
                f"claimed {claimed!r}, derived {classification!r}"
            )
        result = {
            "schema": "a3b-verified-terminal-v1",
            "written_at": now(),
            "status": "verified",
            "classification": classification,
            "executor_identity": executor_identity,
            "verifier_identity": verifier_identity,
            "candidate_terminal_sha256": a3b.sha256_path(candidate_path),
            "stage_receipt_sha256": {
                spec.stage_id: a3b.sha256_path(
                    output_root / "stages" / spec.directory / "receipt.json"
                )
                for spec in a3b.STAGES
                if (output_root / "stages" / spec.directory / "receipt.json").is_file()
            },
            "calculator_evidence_recomputed": source_override is None,
            "source_rehashed": source_override is None,
        }
    except Exception as exc:
        result = {
            "schema": "a3b-verified-terminal-v1",
            "written_at": now(),
            "status": "rejected",
            "classification": a3b.INCONCLUSIVE,
            "executor_identity": candidate.get("executor_identity"),
            "verifier_identity": verifier_identity,
            "detail": f"{type(exc).__name__}: {exc}",
        }
    atomic_json(output_root / "verified-terminal.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--gpu-mem-gb", type=float, default=16.0)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--nice", type=int, default=10)
    parser.add_argument("--log", type=Path)
    parser.add_argument("--verifier-identity", default=DEFAULT_VERIFIER_IDENTITY)
    args = parser.parse_args()
    candidate = json.loads((args.output_root / "candidate-terminal.json").read_text())
    if bool(candidate.get("settings", {}).get("use_gpu")) != args.gpu:
        parser.error("--gpu must match the candidate's recorded use_gpu setting")
    result = verify_experiment(
        args.output_root,
        source_root=args.source_root,
        verifier_identity=args.verifier_identity,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "verified" else 1


if __name__ == "__main__":
    raise SystemExit(main())
