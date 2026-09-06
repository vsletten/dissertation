#!/usr/bin/env python3
"""Bounded A3b two-order release experiment for the Osa-neutral n=1 basin."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from dataclasses import asdict, dataclass, replace
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
        "a3b_proton_microstate_stability",
        default_run_root=DEFAULT_OUTPUT_ROOT / "logs",
        gpu_owner="a3b_proton_microstate_stability",
    )

import numpy as np  # noqa: E402

from quarry.clusters import Cluster, water  # noqa: E402
from quarry.crystal import attack_complex, from_deck_cell  # noqa: E402
from quarry.pipeline import (  # noqa: E402
    DftSettings,
    energy,
    frequencies,
    frequency_geometry_fingerprint,
    frequency_settings_fingerprint,
    gradient,
    optimize_bounded,
)
from quarry.store import geometry_hash  # noqa: E402
from scripts.phase2_ladder import (  # noqa: E402
    ADVISORY_PREOPT_MAX_PROJECTED_FROZEN_DRIFT_A,
    ADVISORY_PREOPT_MAX_RAW_FROZEN_DRIFT_A,
    ADVISORY_PREOPT_MIN_PAIR_A,
    NOISE_FLOOR_CM,
    oxygen_proton_owners,
    preload_cutensor,
)

EXPECTED_CLOSEOUT_MANIFEST_SHA256 = (
    "88039a9cb5f1db921a4bf49d93b195a07df023275199c4d987592b327943ecec"
)
EXPECTED_CONDITIONED_XYZ_SHA256 = (
    "279d168cab2f604c5bb9c3eef83e04adbd0d66bff15a2a35237f750c1d087a47"
)
EXPECTED_CONDITIONED_RECEIPT_SHA256 = (
    "e42635ffd41a72623fe3c13aceb1dc0f5b51926e2da279cd7a577efe4127e8d0"
)
SOURCE_RUN_RELATIVE = Path("runs/phase2/osa-neutral-n1-s2-b3lyp-def2-svp")
MAX_STEPS = 100
OWNER_CONSTRAINTS = ((26, 50), (29, 57))
CONSTRAINT_RESIDUAL_MAX_A = 1.0e-4
GRADIENT_RMS_MAX_HARTREE_PER_BOHR = 3.0e-4
GRADIENT_MAX_HARTREE_PER_BOHR = 4.5e-4
DOWNHILL_MIN_HARTREE = 1.0e-6
EXECUTOR_IDENTITY = "hermes-custom-build-001"
RECOVERY_CANDIDATE = "unverified: recovery candidate"
NO_BASIN_CANDIDATE = "unverified: no-basin candidate"
INCONCLUSIVE = "inconclusive"


@dataclass(frozen=True)
class SourceEvidence:
    cluster: Cluster
    source_root: Path
    closeout_manifest_path: Path
    closeout_manifest_sha256: str
    conditioned_xyz_path: Path
    conditioned_xyz_sha256: str
    conditioned_receipt_path: Path
    conditioned_receipt_sha256: str
    conditioned_geometry_hash: str


@dataclass(frozen=True)
class StageSpec:
    stage_id: str
    directory: str
    parent: str | None
    route: str
    ordinal: int
    constraint_pairs: tuple[tuple[int, int], ...]
    released: str | None
    fully_unconstrained: bool


STAGES = (
    StageSpec(
        "common-dual",
        "common-dual",
        None,
        "common",
        0,
        OWNER_CONSTRAINTS,
        None,
        False,
    ),
    StageSpec(
        "h50-then-h57-single",
        "h50-then-h57/01-release-h50",
        "common-dual",
        "h50-then-h57",
        1,
        ((29, 57),),
        "H50-O26",
        False,
    ),
    StageSpec(
        "h50-then-h57-unconstrained",
        "h50-then-h57/02-release-h57",
        "h50-then-h57-single",
        "h50-then-h57",
        2,
        (),
        "H57-O29",
        True,
    ),
    StageSpec(
        "h57-then-h50-single",
        "h57-then-h50/01-release-h57",
        "common-dual",
        "h57-then-h50",
        1,
        ((26, 50),),
        "H57-O29",
        False,
    ),
    StageSpec(
        "h57-then-h50-unconstrained",
        "h57-then-h50/02-release-h50",
        "h57-then-h50-single",
        "h57-then-h50",
        2,
        (),
        "H50-O26",
        True,
    ),
)
STAGE_BY_ID = {stage.stage_id: stage for stage in STAGES}


@dataclass
class StageOutcome:
    receipt: dict[str, Any]
    endpoint: Cluster | None


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def atomic_xyz(path: Path, cluster: Cluster) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(cluster.to_xyz())
    temporary.replace(path)


def read_cluster(path: Path, template: Cluster | None = None) -> Cluster:
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        raise RuntimeError(f"malformed XYZ: {path}")
    count = int(lines[0])
    records = [line.split() for line in lines[2:]]
    if len(records) != count or any(len(record) != 4 for record in records):
        raise RuntimeError(f"malformed XYZ atom records: {path}")
    symbols = [record[0] for record in records]
    coords = np.asarray([[float(value) for value in record[1:]] for record in records])
    if coords.shape != (count, 3) or not np.all(np.isfinite(coords)):
        raise RuntimeError(f"non-finite or malformed XYZ coordinates: {path}")
    if template is None:
        return Cluster(lines[1], symbols, coords)
    if symbols != template.symbols:
        raise RuntimeError(f"XYZ identity/order mismatch: {path}")
    return replace(template, name=lines[1], coords=coords)


def _canonical_template(repo_root: Path) -> Cluster:
    cell = from_deck_cell(
        repo_root / "petra/examples/kaolinite.toml",
        "Osa",
        metal_shells=2,
        n_intact=1,
        target_charge=0,
    )
    complex_guess, _ = attack_complex(cell, water())
    return complex_guess


def _source_cluster_from_receipt(xyz_path: Path, receipt: dict[str, Any]) -> Cluster:
    signature = receipt.get("signature")
    if not isinstance(signature, dict):
        raise RuntimeError("conditioned receipt omitted signature")
    symbols = signature.get("symbols")
    frozen = signature.get("frozen_indices")
    if not isinstance(symbols, list) or not all(
        isinstance(item, str) for item in symbols
    ):
        raise RuntimeError("conditioned receipt symbols are invalid")
    if not isinstance(frozen, list) or not all(
        isinstance(item, int) for item in frozen
    ):
        raise RuntimeError("conditioned receipt frozen shell is invalid")
    lines = xyz_path.read_text().splitlines()
    template = Cluster(
        lines[1] if len(lines) > 1 else "conditioned",
        list(symbols),
        np.zeros((len(symbols), 3)),
        charge=int(signature.get("charge", 0)),
        spin=int(signature.get("spin", 0)),
        frozen_indices=list(frozen),
    )
    return read_cluster(xyz_path, template)


def validate_source(
    source_root: Path,
    *,
    expected_closeout_sha256: str = EXPECTED_CLOSEOUT_MANIFEST_SHA256,
    expected_conditioned_xyz_sha256: str = EXPECTED_CONDITIONED_XYZ_SHA256,
    expected_conditioned_receipt_sha256: str = EXPECTED_CONDITIONED_RECEIPT_SHA256,
    repo_root: Path | None = None,
) -> SourceEvidence:
    """Rehash and validate A3a's three pinned source artifacts before DFT."""
    source_root = source_root.resolve()
    manifest_path = source_root / "closeout-manifest.json"
    xyz_path = source_root / SOURCE_RUN_RELATIVE / "complex_conditioned.xyz"
    receipt_path = source_root / SOURCE_RUN_RELATIVE / "complex_conditioned.json"
    observed = {
        "closeout": sha256_path(manifest_path),
        "xyz": sha256_path(xyz_path),
        "receipt": sha256_path(receipt_path),
    }
    if observed["closeout"] != expected_closeout_sha256:
        raise RuntimeError(
            f"A3a closeout manifest SHA-256 mismatch: {observed['closeout']}"
        )
    if observed["xyz"] != expected_conditioned_xyz_sha256:
        raise RuntimeError(f"conditioned XYZ SHA-256 mismatch: {observed['xyz']}")
    if observed["receipt"] != expected_conditioned_receipt_sha256:
        raise RuntimeError(
            f"conditioned receipt SHA-256 mismatch: {observed['receipt']}"
        )

    manifest = json.loads(manifest_path.read_text())
    if (
        manifest.get("schema") != "a3a-closeout-manifest-v1"
        or manifest.get("status") != "failed"
    ):
        raise RuntimeError("A3a closeout manifest identity/status mismatch")
    preserved = manifest.get("preserved_evidence")
    if not isinstance(preserved, list):
        raise RuntimeError("A3a closeout manifest omitted preserved evidence")
    preserved_hashes = {
        item.get("sha256") for item in preserved if isinstance(item, dict)
    }
    if not {observed["xyz"], observed["receipt"]} <= preserved_hashes:
        raise RuntimeError(
            "A3a closeout manifest does not bind both conditioned artifacts"
        )

    receipt = json.loads(receipt_path.read_text())
    if receipt.get("schema") != "phase2-reactant-conditioning-v1":
        raise RuntimeError("conditioned receipt schema mismatch")
    cluster = _source_cluster_from_receipt(xyz_path, receipt)
    conditioned_geometry_hash = geometry_hash(cluster.to_xyz())
    if receipt.get("endpoint_geometry_hash") != conditioned_geometry_hash:
        raise RuntimeError("conditioned receipt endpoint geometry hash mismatch")
    owners = oxygen_proton_owners(cluster)
    if owners.get(50) != 26 or owners.get(57) != 29:
        raise RuntimeError(
            "conditioned source lacks required H50:O26 and H57:O29 owners"
        )
    expected_owners = receipt.get("signature", {}).get("oxygen_proton_owners")
    observed_owners = [f"H{h}:O{o}" for h, o in sorted(owners.items())]
    if expected_owners != observed_owners:
        raise RuntimeError("conditioned source owner map drifted from its receipt")

    if (
        expected_closeout_sha256 == EXPECTED_CLOSEOUT_MANIFEST_SHA256
        and expected_conditioned_xyz_sha256 == EXPECTED_CONDITIONED_XYZ_SHA256
        and expected_conditioned_receipt_sha256 == EXPECTED_CONDITIONED_RECEIPT_SHA256
    ):
        root = repo_root or Path(__file__).resolve().parents[2]
        canonical = _canonical_template(root)
        if (
            canonical.symbols != cluster.symbols
            or canonical.charge != cluster.charge
            or canonical.spin != cluster.spin
            or canonical.frozen_indices != cluster.frozen_indices
        ):
            raise RuntimeError(
                "conditioned source identity differs from canonical Osa-neutral n=1"
            )

    return SourceEvidence(
        cluster=cluster,
        source_root=source_root,
        closeout_manifest_path=manifest_path,
        closeout_manifest_sha256=observed["closeout"],
        conditioned_xyz_path=xyz_path,
        conditioned_xyz_sha256=observed["xyz"],
        conditioned_receipt_path=receipt_path,
        conditioned_receipt_sha256=observed["receipt"],
        conditioned_geometry_hash=conditioned_geometry_hash,
    )


def _constraints(source: Cluster, spec: StageSpec) -> list[tuple[int, int, float]]:
    return [
        (
            oxygen,
            hydrogen,
            float(np.linalg.norm(source.coords[oxygen] - source.coords[hydrogen])),
        )
        for oxygen, hydrogen in spec.constraint_pairs
    ]


def _owner_labels(cluster: Cluster) -> list[str]:
    return [f"H{h}:O{o}" for h, o in sorted(oxygen_proton_owners(cluster).items())]


def structural_gate(
    raw: Cluster,
    reference: Cluster,
    constraints: list[tuple[int, int, float]],
    stage: str,
) -> tuple[Cluster, dict[str, Any]]:
    """Project the shell exactly while retaining proton transfer as evidence."""
    if raw.symbols != reference.symbols:
        raise RuntimeError(f"{stage} changed atom identity/order")
    if raw.charge != reference.charge or raw.spin != reference.spin:
        raise RuntimeError(f"{stage} changed charge or spin")
    if raw.frozen_indices != reference.frozen_indices:
        raise RuntimeError(f"{stage} changed the frozen shell")
    if raw.coords.shape != reference.coords.shape or not np.all(
        np.isfinite(raw.coords)
    ):
        raise RuntimeError(f"{stage} produced malformed or non-finite coordinates")
    frozen = sorted(reference.frozen_indices)
    raw_drift = (
        float(np.max(np.abs(raw.coords[frozen] - reference.coords[frozen])))
        if frozen
        else 0.0
    )
    if raw_drift > ADVISORY_PREOPT_MAX_RAW_FROZEN_DRIFT_A:
        raise RuntimeError(
            f"{stage} exceeded raw frozen-shell bound ({raw_drift:.6f} A)"
        )
    coords = raw.coords.copy()
    coords[frozen] = reference.coords[frozen]
    endpoint = replace(raw, coords=coords)
    projected_drift = (
        float(np.max(np.abs(endpoint.coords[frozen] - reference.coords[frozen])))
        if frozen
        else 0.0
    )
    if projected_drift > ADVISORY_PREOPT_MAX_PROJECTED_FROZEN_DRIFT_A or (
        frozen and not np.array_equal(endpoint.coords[frozen], reference.coords[frozen])
    ):
        raise RuntimeError(f"{stage} violated exact frozen shell")
    delta = endpoint.coords[:, None, :] - endpoint.coords[None, :, :]
    distances = np.linalg.norm(delta, axis=2)
    distances[np.diag_indices_from(distances)] = np.inf
    minimum_pair = float(np.min(distances))
    if minimum_pair <= ADVISORY_PREOPT_MIN_PAIR_A:
        raise RuntimeError(f"{stage} produced collision ({minimum_pair:.6f} A)")
    owners = oxygen_proton_owners(endpoint)
    residuals = []
    for oxygen, hydrogen, target in constraints:
        observed = float(
            np.linalg.norm(endpoint.coords[oxygen] - endpoint.coords[hydrogen])
        )
        residual = abs(observed - target)
        residuals.append(
            {
                "oxygen": oxygen,
                "hydrogen": hydrogen,
                "target_a": target,
                "observed_a": observed,
                "absolute_residual_a": residual,
            }
        )
        if residual > CONSTRAINT_RESIDUAL_MAX_A:
            raise RuntimeError(
                f"{stage} active H{hydrogen}-O{oxygen} residual "
                f"{residual:.6g} A exceeds gate"
            )
    reference_owners = oxygen_proton_owners(reference)
    owner_changes = [
        f"H{hydrogen}:O{reference_owners.get(hydrogen)}->O{owners.get(hydrogen)}"
        for hydrogen in sorted(reference_owners.keys() | owners.keys())
        if reference_owners.get(hydrogen) != owners.get(hydrogen)
    ]
    return endpoint, {
        "minimum_pair_distance_a": minimum_pair,
        "raw_maximum_frozen_coordinate_drift_a": raw_drift,
        "projected_maximum_frozen_coordinate_drift_a": projected_drift,
        "constraint_residual_max_a": CONSTRAINT_RESIDUAL_MAX_A,
        "constraint_residuals": residuals,
        "oxygen_proton_owners": [f"H{h}:O{o}" for h, o in sorted(owners.items())],
        "owner_changes": owner_changes,
    }


def _project_gradient(
    cluster: Cluster,
    raw_gradient: np.ndarray,
    constraints: list[tuple[int, int, float]],
) -> np.ndarray:
    projected = np.asarray(raw_gradient, dtype=float).copy()
    if projected.shape != cluster.coords.shape or not np.all(np.isfinite(projected)):
        raise RuntimeError("independent production gradient is malformed or non-finite")
    frozen = set(cluster.frozen_indices)
    for index in frozen:
        projected[index] = 0.0
    for oxygen, hydrogen, _target in constraints:
        direction = cluster.coords[hydrogen] - cluster.coords[oxygen]
        norm = float(np.linalg.norm(direction))
        if not math.isfinite(norm) or norm <= 0.0:
            raise RuntimeError("active constraint has invalid direction")
        unit = direction / norm
        normal = np.zeros_like(projected)
        if oxygen not in frozen:
            normal[oxygen] = -unit
        if hydrogen not in frozen:
            normal[hydrogen] = unit
        denominator = float(np.sum(normal * normal))
        if denominator <= 0.0:
            continue
        projected -= float(np.sum(projected * normal)) / denominator * normal
    return projected


def gradient_metrics(
    cluster: Cluster,
    settings: DftSettings,
    constraints: list[tuple[int, int, float]],
) -> dict[str, Any]:
    projected = _project_gradient(cluster, gradient(cluster, settings), constraints)
    free = sorted(set(range(len(cluster.symbols))) - set(cluster.frozen_indices))
    values = projected[free]
    rms = float(np.sqrt(np.mean(values**2))) if values.size else 0.0
    maximum = float(np.max(np.abs(values))) if values.size else 0.0
    return {
        "projected_gradient_hartree_per_bohr": projected.tolist(),
        "projected_free_atom_indices": free,
        "projected_gradient_rms_hartree_per_bohr": rms,
        "projected_gradient_max_hartree_per_bohr": maximum,
        "gradient_rms_threshold_hartree_per_bohr": GRADIENT_RMS_MAX_HARTREE_PER_BOHR,
        "gradient_max_threshold_hartree_per_bohr": GRADIENT_MAX_HARTREE_PER_BOHR,
        "passed": (
            math.isfinite(rms)
            and math.isfinite(maximum)
            and rms <= GRADIENT_RMS_MAX_HARTREE_PER_BOHR
            and maximum <= GRADIENT_MAX_HARTREE_PER_BOHR
        ),
    }


def _stage_signature(
    source: SourceEvidence,
    settings: DftSettings,
    spec: StageSpec,
    constraints: list[tuple[int, int, float]],
) -> dict[str, Any]:
    return {
        "schema": "a3b-stage-signature-v1",
        "stage": spec.stage_id,
        "route": spec.route,
        "ordinal": spec.ordinal,
        "parent": spec.parent,
        "source": {
            "closeout_manifest_sha256": source.closeout_manifest_sha256,
            "conditioned_xyz_sha256": source.conditioned_xyz_sha256,
            "conditioned_receipt_sha256": source.conditioned_receipt_sha256,
            "conditioned_geometry_hash": source.conditioned_geometry_hash,
        },
        "settings": asdict(settings),
        "constraints": [[i, j, target] for i, j, target in constraints],
        "released": spec.released,
        "budget": {"max_steps": MAX_STEPS, "continuation_allowed": False},
        "fresh_optimizer": True,
        "geometric_default_fresh_hessian": True,
    }


def _resume_stage(
    stage_dir: Path,
    spec: StageSpec,
    template: Cluster,
    signature: dict[str, Any],
    seed: Cluster,
    parent: dict[str, Any],
) -> StageOutcome | None:
    reservation_path = stage_dir / "reservation.json"
    receipt_path = stage_dir / "receipt.json"
    if not reservation_path.exists() and not receipt_path.exists():
        return None
    if not reservation_path.exists():
        return StageOutcome(
            {"status": "orphaned-receipt", "stage": spec.stage_id}, None
        )
    if not receipt_path.exists():
        return StageOutcome(
            {"status": "orphaned-reservation", "stage": spec.stage_id}, None
        )
    try:
        reservation = json.loads(reservation_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        return StageOutcome(
            {
                "status": "reservation-invalid",
                "stage": spec.stage_id,
                "detail": f"{type(exc).__name__}: {exc}",
            },
            None,
        )
    expected_seed = frequency_geometry_fingerprint(seed)
    if (
        reservation.get("status") != "reserved"
        or reservation.get("signature") != signature
        or reservation.get("parent") != parent
        or reservation.get("seed", {}).get("geometry_fingerprint") != expected_seed
    ):
        return StageOutcome(
            {"status": "reservation-identity-mismatch", "stage": spec.stage_id},
            None,
        )
    try:
        receipt = json.loads(receipt_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        return StageOutcome(
            {
                "status": "receipt-invalid",
                "stage": spec.stage_id,
                "detail": f"{type(exc).__name__}: {exc}",
            },
            None,
        )
    if (
        receipt.get("signature") != signature
        or receipt.get("parent") != parent
        or receipt.get("seed", {}).get("geometry_fingerprint") != expected_seed
    ):
        return StageOutcome(
            {"status": "receipt-identity-mismatch", "stage": spec.stage_id}, None
        )
    if receipt.get("status") != "complete":
        status = receipt.get("status", "orphaned-reservation")
        return StageOutcome({**receipt, "status": status}, None)
    endpoint_path = stage_dir / "endpoint.xyz"
    if not endpoint_path.is_file():
        return StageOutcome({**receipt, "status": "orphaned-complete-receipt"}, None)
    if receipt.get("endpoint", {}).get("sha256") != sha256_path(endpoint_path):
        return StageOutcome({**receipt, "status": "endpoint-hash-mismatch"}, None)
    endpoint = read_cluster(endpoint_path, template)
    if receipt.get("endpoint", {}).get(
        "geometry_fingerprint"
    ) != frequency_geometry_fingerprint(endpoint):
        return StageOutcome({**receipt, "status": "endpoint-geometry-mismatch"}, None)
    return StageOutcome(receipt, endpoint)


def run_stage(
    output_root: Path,
    source: SourceEvidence,
    settings: DftSettings,
    spec: StageSpec,
    seed: Cluster,
    parent_receipt: dict[str, Any] | None,
) -> StageOutcome:
    stage_dir = output_root / "stages" / spec.directory
    constraints = _constraints(source.cluster, spec)
    signature = _stage_signature(source, settings, spec, constraints)
    parent = {
        "stage": spec.parent,
        "receipt_sha256": (
            sha256_path(
                output_root
                / "stages"
                / STAGE_BY_ID[spec.parent].directory
                / "receipt.json"
            )
            if spec.parent
            else None
        ),
    }
    resumed = _resume_stage(stage_dir, spec, source.cluster, signature, seed, parent)
    if resumed is not None:
        return resumed
    stage_dir.mkdir(parents=True, exist_ok=True)
    reservation = {
        "schema": "a3b-stage-reservation-v1",
        "status": "reserved",
        "reserved_at": now(),
        "signature": signature,
        "seed": {
            "sha256": hashlib.sha256(seed.to_xyz().encode()).hexdigest(),
            "geometry_fingerprint": frequency_geometry_fingerprint(seed),
        },
        "parent": parent,
    }
    atomic_json(stage_dir / "reservation.json", reservation)
    try:
        kwargs: dict[str, Any] = {"max_steps": MAX_STEPS}
        if constraints:
            kwargs["fixed_distances"] = constraints
        optimized = optimize_bounded(seed, settings, **kwargs)
    except Exception as exc:
        receipt = {
            "schema": "a3b-stage-receipt-v1",
            "status": "optimizer-failed",
            "terminal_stage": "optimizer-call",
            "stage": spec.stage_id,
            "route": spec.route,
            "ordinal": spec.ordinal,
            "parent": parent,
            "signature": signature,
            "source": signature["source"],
            "settings": asdict(settings),
            "seed": {
                "geometry_hash": geometry_hash(seed.to_xyz()),
                "geometry_fingerprint": frequency_geometry_fingerprint(seed),
            },
            "constraints": signature["constraints"],
            "released": spec.released,
            "fully_unconstrained": spec.fully_unconstrained,
            "detail": f"{type(exc).__name__}: {exc}",
            "budget": signature["budget"],
            "optimizer": {
                "converged": False,
                "fresh_instance": True,
                "geometric_default_fresh_hessian": True,
            },
        }
        atomic_json(stage_dir / "receipt.json", receipt)
        return StageOutcome(receipt, None)

    atomic_xyz(stage_dir / "raw-endpoint.xyz", optimized.cluster)
    receipt = {
        "schema": "a3b-stage-receipt-v1",
        "status": "pending-structural-gates",
        "terminal_stage": "structural-gates",
        "stage": spec.stage_id,
        "route": spec.route,
        "ordinal": spec.ordinal,
        "parent": parent,
        "signature": signature,
        "source": signature["source"],
        "settings": asdict(settings),
        "seed": {
            "geometry_hash": geometry_hash(seed.to_xyz()),
            "geometry_fingerprint": frequency_geometry_fingerprint(seed),
        },
        "constraints": signature["constraints"],
        "released": spec.released,
        "fully_unconstrained": spec.fully_unconstrained,
        "budget": signature["budget"],
        "optimizer": {
            "converged": bool(optimized.converged),
            "fresh_instance": True,
            "geometric_default_fresh_hessian": True,
        },
        "raw_endpoint": {
            "path": "raw-endpoint.xyz",
            "sha256": sha256_path(stage_dir / "raw-endpoint.xyz"),
            "geometry_hash": geometry_hash(optimized.cluster.to_xyz()),
            "geometry_fingerprint": frequency_geometry_fingerprint(optimized.cluster),
            "charge": optimized.cluster.charge,
            "spin": optimized.cluster.spin,
            "frozen_indices": list(optimized.cluster.frozen_indices),
        },
    }
    atomic_json(stage_dir / "receipt.json", receipt)

    try:
        endpoint, raw_structure = structural_gate(
            optimized.cluster, source.cluster, constraints, spec.stage_id
        )
        atomic_xyz(stage_dir / "endpoint.xyz", endpoint)
        endpoint = read_cluster(stage_dir / "endpoint.xyz", source.cluster)
        structure = structural_gate(
            endpoint, source.cluster, constraints, spec.stage_id
        )[1]
        structure["raw_maximum_frozen_coordinate_drift_a"] = raw_structure[
            "raw_maximum_frozen_coordinate_drift_a"
        ]
        receipt["terminal_stage"] = "energy-evidence"
        observed_energy = float(energy(endpoint, settings))
        if not math.isfinite(observed_energy):
            raise RuntimeError("production energy is non-finite")
        receipt["terminal_stage"] = "gradient-evidence"
        metrics = gradient_metrics(endpoint, settings, constraints)
        owners = _owner_labels(endpoint)
        owner_retaining = owners == _owner_labels(source.cluster)
        stationary = bool(optimized.converged) and bool(metrics["passed"])
        phva: dict[str, Any] = {"status": "not-required", "imaginary_cm": []}
        if spec.fully_unconstrained and owner_retaining and stationary:
            receipt["terminal_stage"] = "phva-evidence"
            frequency = frequencies(endpoint, settings)
            imaginary = [float(value) for value in frequency.imaginary_cm]
            if not math.isfinite(float(frequency.electronic_hartree)) or not all(
                math.isfinite(value) for value in imaginary
            ):
                raise RuntimeError("PHVA returned non-finite energy/frequencies")
            phva = {
                "status": (
                    "passed"
                    if not any(value > NOISE_FLOOR_CM for value in imaginary)
                    else "failed"
                ),
                "imaginary_cm": imaginary,
                "noise_floor_cm": NOISE_FLOOR_CM,
                "fresh_hessian": True,
                "geometry_fingerprint": frequency_geometry_fingerprint(endpoint),
                "settings_fingerprint": frequency_settings_fingerprint(settings),
                "electronic_hartree": float(frequency.electronic_hartree),
            }
        receipt.update(
            {
                "status": "complete",
                "terminal_stage": "stage-complete",
                "completed_at": now(),
                "structure": structure,
                "endpoint": {
                    "path": "endpoint.xyz",
                    "sha256": sha256_path(stage_dir / "endpoint.xyz"),
                    "geometry_hash": geometry_hash(endpoint.to_xyz()),
                    "geometry_fingerprint": frequency_geometry_fingerprint(endpoint),
                    "charge": endpoint.charge,
                    "spin": endpoint.spin,
                    "frozen_indices": list(endpoint.frozen_indices),
                },
                "reference_owners": _owner_labels(source.cluster),
                "observed_owners": owners,
                "owner_changes": structure["owner_changes"],
                "owner_retaining": owner_retaining,
                "stationary": stationary,
                "evidence": {"energy_hartree": observed_energy, **metrics},
                "phva": phva,
            }
        )
    except Exception as exc:
        receipt.update(
            {
                "status": "gate-rejected",
                "terminal_stage": receipt["terminal_stage"],
                "completed_at": now(),
                "detail": f"{type(exc).__name__}: {exc}",
                "stationary": False,
            }
        )
        atomic_json(stage_dir / "receipt.json", receipt)
        return StageOutcome(receipt, None)

    atomic_json(stage_dir / "receipt.json", receipt)
    return StageOutcome(receipt, endpoint)


def _finite_energy(receipt: dict[str, Any]) -> float | None:
    value = receipt.get("evidence", {}).get("energy_hartree")
    if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
        return None
    return float(value)


def _recovery(receipt: dict[str, Any]) -> bool:
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
        and not any(float(value) > NOISE_FLOOR_CM for value in imaginary)
        and _finite_energy(receipt) is not None
    )


def classify(
    common: dict[str, Any],
    final_a: dict[str, Any],
    final_b: dict[str, Any],
) -> str:
    """Classify only a candidate; independent verification publishes the verdict."""
    if _recovery(final_a) or _recovery(final_b):
        return RECOVERY_CANDIDATE
    common_energy = _finite_energy(common)
    energies = [_finite_energy(final_a), _finite_energy(final_b)]
    if common_energy is None or any(value is None for value in energies):
        return INCONCLUSIVE
    finals = (final_a, final_b)
    no_basin = all(
        receipt.get("status") == "complete"
        and receipt.get("fully_unconstrained") is True
        and receipt.get("owner_retaining") is False
        and receipt.get("stationary") is True
        for receipt in finals
    ) and all(common_energy - float(value) > DOWNHILL_MIN_HARTREE for value in energies)
    return NO_BASIN_CANDIDATE if no_basin else INCONCLUSIVE


def _missing(stage_id: str, detail: str) -> StageOutcome:
    return StageOutcome(
        {"status": "not-run", "stage": stage_id, "detail": detail}, None
    )


def run_experiment(
    output_root: Path,
    source: SourceEvidence,
    *,
    use_gpu: bool,
    executor_identity: str = EXECUTOR_IDENTITY,
) -> dict[str, Any]:
    """Execute each predeclared optimizer budget at most once."""
    output_root.mkdir(parents=True, exist_ok=True)
    settings = DftSettings(
        xc="b3lyp", basis="def2-svp", density_fit=True, use_gpu=use_gpu
    )
    outcomes: dict[str, StageOutcome] = {}
    common_spec = STAGES[0]
    common = run_stage(output_root, source, settings, common_spec, source.cluster, None)
    outcomes[common_spec.stage_id] = common

    if common.endpoint is not None and common.receipt.get("stationary") is True:
        first_a = run_stage(
            output_root, source, settings, STAGES[1], common.endpoint, common.receipt
        )
        outcomes[STAGES[1].stage_id] = first_a
        outcomes[STAGES[2].stage_id] = (
            run_stage(
                output_root,
                source,
                settings,
                STAGES[2],
                first_a.endpoint,
                first_a.receipt,
            )
            if first_a.endpoint is not None
            else _missing(STAGES[2].stage_id, "parent stage produced no endpoint")
        )
        first_b = run_stage(
            output_root, source, settings, STAGES[3], common.endpoint, common.receipt
        )
        outcomes[STAGES[3].stage_id] = first_b
        outcomes[STAGES[4].stage_id] = (
            run_stage(
                output_root,
                source,
                settings,
                STAGES[4],
                first_b.endpoint,
                first_b.receipt,
            )
            if first_b.endpoint is not None
            else _missing(STAGES[4].stage_id, "parent stage produced no endpoint")
        )
    else:
        for spec in STAGES[1:]:
            outcomes[spec.stage_id] = _missing(
                spec.stage_id, "common dual-constrained stationary seed unavailable"
            )

    final_a = outcomes[STAGES[2].stage_id].receipt
    final_b = outcomes[STAGES[4].stage_id].receipt
    classification = classify(common.receipt, final_a, final_b)
    terminal = {
        "schema": "a3b-candidate-terminal-v1",
        "written_at": now(),
        "executor_identity": executor_identity,
        "classification": classification,
        "source": {
            "closeout_manifest_sha256": source.closeout_manifest_sha256,
            "conditioned_xyz_sha256": source.conditioned_xyz_sha256,
            "conditioned_receipt_sha256": source.conditioned_receipt_sha256,
            "conditioned_geometry_hash": source.conditioned_geometry_hash,
        },
        "settings": asdict(settings),
        "stages": {
            stage_id: outcome.receipt.get("status", "unknown")
            for stage_id, outcome in outcomes.items()
        },
        "independent_verification_required": True,
        "forbidden_outputs_emitted": False,
    }
    if classification == INCONCLUSIVE:
        terminal["detail"] = (
            "bounded release experiment was inconclusive; no stage may be replayed "
            "without a new executable mechanism card"
        )
    atomic_json(output_root / "candidate-terminal.json", terminal)
    return terminal


def dry_run(source: SourceEvidence, *, use_gpu: bool) -> dict[str, Any]:
    settings = DftSettings(
        xc="b3lyp", basis="def2-svp", density_fit=True, use_gpu=use_gpu
    )
    return {
        "schema": "a3b-dry-run-v1",
        "source": {
            "conditioned_xyz_sha256": source.conditioned_xyz_sha256,
            "conditioned_geometry_hash": source.conditioned_geometry_hash,
            "owners": _owner_labels(source.cluster),
        },
        "settings": asdict(settings),
        "stages": [
            _stage_signature(source, settings, spec, _constraints(source.cluster, spec))
            for spec in STAGES
        ],
        "production_artifacts_written": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--gpu-mem-gb", type=float, default=16.0)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--nice", type=int, default=10)
    parser.add_argument("--log")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--executor-identity", default=EXECUTOR_IDENTITY)
    args = parser.parse_args()
    if args.threads > 16:
        parser.error("--threads must be <=16")
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.environ["MKL_NUM_THREADS"] = str(args.threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(args.threads)
    if args.gpu:
        preload_cutensor()
    source = validate_source(args.source_root)
    if args.dry_run:
        print(json.dumps(dry_run(source, use_gpu=args.gpu), indent=2, sort_keys=True))
        return 0
    terminal = run_experiment(
        args.output_root,
        source,
        use_gpu=args.gpu,
        executor_identity=args.executor_identity,
    )
    print(json.dumps(terminal, indent=2, sort_keys=True))
    return 0 if terminal["classification"] != INCONCLUSIVE else 2


if __name__ == "__main__":
    raise SystemExit(main())
