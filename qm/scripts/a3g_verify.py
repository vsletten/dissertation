#!/usr/bin/env python3
"""Cold optimizer-free verifier for the bounded A3g H52 owner-basin run."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import os
import shutil
import sys
import time
from contextlib import contextmanager
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    bootstrap_cli(
        "a3g_verify",
        default_run_root=Path(
            "/mnt/data/vsletten/dissertation-data/"
            "a3g-oaa-neutral-n2-proton-microstate-stability/logs"
        ),
        gpu_owner="a3g_verify",
    )

import numpy as np  # noqa: E402

from quarry.clusters import Cluster, water  # noqa: E402
from quarry.crystal import attack_complex, from_deck_cell  # noqa: E402
from quarry.pipeline import (  # noqa: E402
    energy,
    frequencies,
    frequency_geometry_fingerprint,
    frequency_settings_fingerprint,
    gradient,
)
from quarry.store import geometry_hash  # noqa: E402
from scripts import a3g_oaa_owner_basin as a3g  # noqa: E402
from scripts.phase2_ladder import preload_cutensor  # noqa: E402


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
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    temporary.replace(path)


def canonical_template(repo_root: Path) -> Cluster:
    cell = from_deck_cell(
        repo_root / "petra/examples/kaolinite.toml",
        "Oaa",
        center_index=18,
        metal_shells=2,
        n_intact=2,
        target_charge=0,
    )
    cluster, _product = attack_complex(cell, water())
    return cluster


def validate_source(
    source_root: Path,
    *,
    expected_hashes: dict[str, str] | None = None,
    repo_root: Path | None = None,
) -> a3g.SourceEvidence:
    """Rehash and semantically bind every artifact in the A3g evidence contract."""
    source_root = source_root.resolve()
    expected = expected_hashes or a3g.EXPECTED_SOURCE_SHA256
    if set(expected) != set(a3g.SOURCE_PATHS):
        raise RuntimeError("source expected-hash set is incomplete")
    observed = {
        name: sha256_path(source_root / relative)
        for name, relative in a3g.SOURCE_PATHS.items()
    }
    for name, digest in expected.items():
        if observed[name] != digest:
            raise RuntimeError(f"source {name} SHA-256 mismatch: {observed[name]}")

    family = json.loads((source_root / a3g.SOURCE_PATHS["family_receipt"]).read_text())
    if (
        family.get("schema") != "a3-family-campaign-terminal-v1"
        or family.get("success") is not False
        or family.get("family") != "oaa"
        or family.get("state") != "neutral"
        or family.get("current_n_intact") != 2
        or family.get("completed") != []
        or family.get("expected_git_sha") != a3g.EXPECTED_EXECUTION_SOURCE
        or family.get("observed_git_sha") != a3g.EXPECTED_EXECUTION_SOURCE
    ):
        raise RuntimeError("family terminal identity/outcome/source mismatch")

    metadata = json.loads((source_root / a3g.SOURCE_PATHS["metadata"]).read_text())
    if (
        metadata.get("site_kind") != "Oaa"
        or metadata.get("state") != "neutral"
        or metadata.get("n_intact") != 2
        or metadata.get("center_site") != 18
        or metadata.get("metal_shells") != 2
        or metadata.get("charge") != 0
        or metadata.get("method") != "b3lyp/def2-svp/df"
        or metadata.get("driver_git_commit") != a3g.EXPECTED_EXECUTION_SOURCE
    ):
        raise RuntimeError("source metadata identity/settings mismatch")

    root = (repo_root or Path(__file__).resolve().parents[2]).resolve()
    template = canonical_template(root)
    seed = read_cluster(source_root / a3g.SOURCE_PATHS["seed"], template)
    if (
        seed.symbols != template.symbols
        or seed.charge != template.charge
        or seed.spin != template.spin
        or seed.frozen_indices != template.frozen_indices
        or len(seed.symbols) != 58
    ):
        raise RuntimeError("source seed identity/order/state/frozen-shell mismatch")
    owners = oxygen_proton_owners(seed)
    if owners.get(a3g.H52) != a3g.O15:
        raise RuntimeError(f"source seed H52 owner is O{owners.get(a3g.H52)}, not O15")
    if metadata.get("n_atoms") + 3 != len(seed.symbols):
        raise RuntimeError("metadata cell atom count does not bind attacked complex")
    return a3g.SourceEvidence(seed, source_root, observed, metadata)


def source_map(source: a3g.SourceEvidence) -> dict[str, Any]:
    return {
        "root": str(source.source_root),
        "hashes": dict(sorted(source.hashes.items())),
        "execution_source": a3g.EXPECTED_EXECUTION_SOURCE,
        "branch": a3g.EXPECTED_BRANCH,
        "seed_geometry_hash": geometry_hash(source.cluster.to_xyz()),
        "seed_geometry_fingerprint": frequency_geometry_fingerprint(source.cluster),
        "symbols": list(source.cluster.symbols),
        "charge": source.cluster.charge,
        "spin": source.cluster.spin,
        "frozen_indices": list(source.cluster.frozen_indices),
        "owners": _owner_labels(source.cluster),
    }


@contextmanager
def exclusive_run(output_root: Path):
    lock_path = output_root / "run.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("A3g experiment is already active") from exc
        try:
            yield
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def revoke_terminals(output_root: Path) -> None:
    paths = [
        output_root / name
        for name in ("candidate-terminal.json", "verified-terminal.json")
        if (output_root / name).exists()
    ]
    if not paths:
        return
    revoked = output_root / "revoked" / f"{time.time_ns()}-{os.getpid()}"
    revoked.mkdir(parents=True, exist_ok=False)
    for path in paths:
        shutil.move(path, revoked / path.name)


def stage_signature(
    source: a3g.SourceEvidence,
    spec: a3g.StageSpec,
    code_revision: str,
) -> dict[str, Any]:
    active = constraints(source, spec)
    return {
        "schema": "a3g-stage-signature-v1",
        "stage": spec.stage_id,
        "method": spec.method,
        "settings": asdict(spec.settings),
        "source": source_map(source),
        "constraints": [list(item) for item in active],
        "budget": {"max_steps": a3g.MAX_STEPS, "retry_allowed": False},
        "parent": spec.parent,
        "fresh_optimizer": True,
        "geometric_default_fresh_hessian": True,
        "fully_released": spec.fully_released,
        "code_revision": code_revision,
    }


def _artifact(cluster: Cluster, path: Path) -> dict[str, Any]:
    return {
        "path": path.name,
        "sha256": sha256_path(path),
        "geometry_hash": geometry_hash(cluster.to_xyz()),
        "geometry_fingerprint": frequency_geometry_fingerprint(cluster),
        "charge": cluster.charge,
        "spin": cluster.spin,
        "frozen_indices": list(cluster.frozen_indices),
    }


def _experiment_budget(receipts: dict[str, dict[str, Any]]) -> dict[str, int]:
    mapping = {
        "owner-conditioning": "owner_conditioning",
        "constrained-production": "constrained_production",
        "released-production": "released_production",
    }
    budget = {key: 0 for key in mapping.values()}
    budget["retries"] = 0
    for stage_id, key in mapping.items():
        optimizer = (receipts.get(stage_id) or {}).get("optimizer") or {}
        calls = optimizer.get("observed_calls")
        retries = optimizer.get("observed_retries")
        if isinstance(calls, int):
            budget[key] = calls
        if isinstance(retries, int):
            budget["retries"] += retries
    return budget


def _finite_energy(receipt: dict[str, Any]) -> float | None:
    value = receipt.get("evidence", {}).get("energy_hartree")
    if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
        return None
    return float(value)


def read_cluster(path: Path, template: Cluster) -> Cluster:
    lines = path.read_text().splitlines()
    count = int(lines[0])
    rows = [line.split() for line in lines[2:]]
    if (
        count != len(template.symbols)
        or len(rows) != count
        or any(len(row) != 4 for row in rows)
    ):
        raise RuntimeError("invalid XYZ atom records")
    if [row[0] for row in rows] != template.symbols:
        raise RuntimeError("XYZ atom order mismatch")
    coords = np.array([[float(x) for x in row[1:]] for row in rows])
    if not np.isfinite(coords).all():
        raise RuntimeError("non-finite XYZ")
    return replace(template, name=lines[1], coords=coords)


def oxygen_proton_owners(cluster: Cluster) -> dict[int, int]:
    # Distance matrix and unique membership rather than executor's sorted pairs.
    hydrogens = np.flatnonzero(np.array(cluster.symbols) == "H")
    oxygens = np.flatnonzero(np.array(cluster.symbols) == "O")
    if len(hydrogens) and not len(oxygens):
        raise RuntimeError("no oxygen owners")
    owners = {}
    for h in hydrogens:
        distances = np.sqrt(
            np.sum((cluster.coords[oxygens] - cluster.coords[h]) ** 2, axis=1)
        )
        bonded = np.flatnonzero(distances <= 1.25)
        if len(bonded) != 1:
            raise RuntimeError("proton must have exactly one oxygen within 1.25 A")
        owner = int(bonded[0])
        others = np.delete(distances, owner)
        if len(others) and np.min(others) - distances[owner] < 0.15:
            raise RuntimeError("ambiguous proton owner margin")
        owners[int(h)] = int(oxygens[owner])
    return owners


def _owner_labels(cluster: Cluster) -> list[str]:
    return [f"H{h}:O{o}" for h, o in oxygen_proton_owners(cluster).items()]


def constraints(source: a3g.SourceEvidence, spec: a3g.StageSpec) -> list:
    if spec.fully_released:
        return []
    return [
        (
            15,
            52,
            float(
                np.linalg.norm(source.cluster.coords[52] - source.cluster.coords[15])
            ),
        )
    ]


def inspect_structure(
    raw: Cluster, endpoint: Cluster, source: a3g.SourceEvidence, active: list
) -> dict:
    reference = source.cluster
    for cluster in (raw, endpoint):
        if (
            cluster.symbols != reference.symbols
            or cluster.charge != reference.charge
            or cluster.spin != reference.spin
            or cluster.frozen_indices != reference.frozen_indices
            or cluster.coords.shape != reference.coords.shape
            or not np.isfinite(cluster.coords).all()
        ):
            raise RuntimeError("structural identity mismatch")
    frozen = reference.frozen_indices
    drift = max(
        (
            abs(float(raw.coords[i, j] - reference.coords[i, j]))
            for i in frozen
            for j in range(3)
        ),
        default=0.0,
    )
    if drift > 0.020:
        raise RuntimeError("raw frozen-shell drift")
    expected = raw.coords.copy()
    expected[frozen] = reference.coords[frozen]
    if not np.array_equal(expected, endpoint.coords):
        raise RuntimeError("raw/projected endpoint relationship mismatch")
    if not np.array_equal(endpoint.coords[frozen], reference.coords[frozen]):
        raise RuntimeError("projected frozen shell not exact")
    minimum = min(
        math.dist(endpoint.coords[i], endpoint.coords[j])
        for i in range(len(endpoint.symbols))
        for j in range(i)
    )
    if minimum <= 0.60:
        raise RuntimeError("endpoint collision")
    owners = oxygen_proton_owners(endpoint)
    original = oxygen_proton_owners(reference)
    residuals = []
    for o, h, target in active:
        observed = float(np.linalg.norm(endpoint.coords[o] - endpoint.coords[h]))
        residual = abs(observed - target)
        if residual > 1.0e-4:
            raise RuntimeError("constraint residual")
        residuals.append(
            {
                "oxygen": o,
                "hydrogen": h,
                "target_a": target,
                "observed_a": observed,
                "absolute_residual_a": residual,
            }
        )
    return {
        "minimum_pair_distance_a": minimum,
        "raw_maximum_frozen_coordinate_drift_a": drift,
        "projected_maximum_frozen_coordinate_drift_a": 0.0,
        "constraint_residual_max_a": 1.0e-4,
        "constraint_residuals": residuals,
        "oxygen_proton_owners": [f"H{h}:O{o}" for h, o in owners.items()],
        "owner_changes": [
            f"H{h}:O{original[h]}->O{owners[h]}"
            for h in original
            if original[h] != owners[h]
        ],
    }


def independent_metrics(cluster: Cluster, raw: Any, active: list) -> dict:
    values = np.asarray(raw, dtype=float)
    if values.shape != cluster.coords.shape or not np.isfinite(values).all():
        raise RuntimeError("non-finite full gradient")
    free = [i for i in range(len(cluster.symbols)) if i not in cluster.frozen_indices]
    flat = values[free].reshape(-1)
    # Orthogonal projection via a least-squares constraint Jacobian in free DOFs.
    jacobian = np.zeros((len(active), len(flat)))
    for row, (o, h, _) in enumerate(active):
        direction = cluster.coords[h] - cluster.coords[o]
        direction /= np.linalg.norm(direction)
        for index, atom in enumerate(free):
            if atom == o:
                jacobian[row, index * 3 : index * 3 + 3] = -direction
            elif atom == h:
                jacobian[row, index * 3 : index * 3 + 3] = direction
    if len(active):
        flat = flat - jacobian.T @ np.linalg.lstsq(jacobian.T, flat, rcond=None)[0]
    projected = np.zeros_like(values)
    projected[free] = flat.reshape(-1, 3)
    rms = math.sqrt(float(flat @ flat) / len(flat)) if len(flat) else 0.0
    maximum = max(map(abs, flat), default=0.0)
    return {
        "projected_gradient_hartree_per_bohr": projected.tolist(),
        "projected_free_atom_indices": free,
        "projected_gradient_rms_hartree_per_bohr": rms,
        "projected_gradient_max_hartree_per_bohr": float(maximum),
        "gradient_rms_threshold_hartree_per_bohr": 3.0e-4,
        "gradient_max_threshold_hartree_per_bohr": 4.5e-4,
        "passed": rms <= 3.0e-4 and maximum <= 4.5e-4,
    }


def reproduce(cluster: Cluster, settings: Any, saved: dict) -> dict:
    value = float(energy(cluster, settings))
    full = np.asarray(gradient(cluster, settings), dtype=float)
    prior = np.asarray(saved.get("gradient_hartree_per_bohr"), dtype=float)
    if (
        not math.isfinite(value)
        or full.shape != cluster.coords.shape
        or prior.shape != full.shape
        or not np.isfinite(full).all()
        or not np.isfinite(prior).all()
        or not isinstance(saved.get("energy_hartree"), (int, float))
        or not math.isfinite(saved["energy_hartree"])
        or abs(value - saved["energy_hartree"]) > a3g.ENERGY_REPRO_TOL
        or not np.allclose(full, prior, atol=a3g.GRADIENT_REPRO_TOL, rtol=0)
        or saved.get("settings") != asdict(settings)
        or saved.get("geometry_fingerprint") != frequency_geometry_fingerprint(cluster)
    ):
        raise RuntimeError(
            "raw/endpoint energy, full gradient or settings not reproducible"
        )
    return {"energy_hartree": value, "gradient_hartree_per_bohr": full.tolist()}


def phva(cluster: Cluster, settings: Any, value: float) -> dict:
    result = frequencies(cluster, settings)
    imaginary = [float(v) for v in result.imaginary_cm]
    if (
        not math.isfinite(float(result.electronic_hartree))
        or not all(math.isfinite(v) and v >= 0 for v in imaginary)
        or abs(float(result.electronic_hartree) - value) > a3g.ENERGY_REPRO_TOL
        or result.geometry_fingerprint != frequency_geometry_fingerprint(cluster)
        or result.settings_fingerprint != frequency_settings_fingerprint(settings)
    ):
        raise RuntimeError("PHVA evidence mismatch")
    return {
        "status": "passed"
        if max(imaginary, default=0) <= a3g.NOISE_FLOOR_CM
        else "failed",
        "imaginary_cm": imaginary,
        "electronic_hartree": float(result.electronic_hartree),
        "fresh_hessian": True,
        "noise_floor_cm": a3g.NOISE_FLOOR_CM,
        "geometry_fingerprint": frequency_geometry_fingerprint(cluster),
        "settings_fingerprint": frequency_settings_fingerprint(settings),
    }


def _read_stage(
    output_root,
    source,
    spec,
    code_revision,
    parent,
    *,
    recompute_calculators,
    expected_seed,
):
    stage_dir = output_root / "stages" / spec.directory
    reservation = json.loads((stage_dir / "reservation.json").read_text())
    receipt = json.loads((stage_dir / "receipt.json").read_text())
    signature = stage_signature(source, spec, code_revision)
    if (
        reservation.get("schema") != "a3g-stage-reservation-v1"
        or reservation.get("status") != "reserved"
        or receipt.get("schema") != "a3g-stage-receipt-v1"
        or any(
            item.get("signature") != signature or item.get("parent") != parent
            for item in (reservation, receipt)
        )
    ):
        raise RuntimeError("reservation/receipt identity mismatch")
    expected_input = {
        "sha256": hashlib.sha256(expected_seed.to_xyz().encode()).hexdigest(),
        "geometry_fingerprint": frequency_geometry_fingerprint(expected_seed),
    }
    if (
        reservation.get("seed") != expected_input
        or receipt.get("seed") != expected_input
    ):
        raise RuntimeError("stage input seed/parent mismatch")
    optimizer = receipt.get("optimizer", {})
    if (
        optimizer.get("fresh_instance") is not True
        or optimizer.get("geometric_default_fresh_hessian") is not True
        or optimizer.get("observed_calls") != 1
        or optimizer.get("observed_retries") != 0
        or optimizer.get("observed_max_steps") != 100
    ):
        raise RuntimeError("optimizer budget/freshness mismatch")
    if receipt.get("status") == "optimizer-failed":
        if (stage_dir / "raw-endpoint.xyz").exists():
            raise RuntimeError("optimizer failure has raw endpoint")
        return receipt, None
    clusters = {}
    for key, filename in (
        ("raw_endpoint", "raw-endpoint.xyz"),
        ("endpoint", "endpoint.xyz"),
    ):
        path = stage_dir / filename
        if key == "endpoint" and receipt.get("status") != "complete":
            break
        if not path.is_file() or receipt.get(key, {}).get("sha256") != sha256_path(
            path
        ):
            raise RuntimeError(f"{key} hash-mismatched")
        cluster = read_cluster(path, source.cluster)
        if receipt.get(key) != _artifact(cluster, path):
            raise RuntimeError(f"{key} artifact mismatch")
        clusters[key] = cluster
    raw = clusters["raw_endpoint"]
    reproduce(raw, spec.settings, receipt.get("raw_evidence", {}))
    if receipt.get("status") != "complete":
        projected_coords = raw.coords.copy()
        projected_coords[source.cluster.frozen_indices] = source.cluster.coords[
            source.cluster.frozen_indices
        ]
        projected = replace(raw, coords=projected_coords)
        try:
            inspect_structure(raw, projected, source, constraints(source, spec))
        except RuntimeError:
            # A fresh, independently computed integrity failure closes the route.
            return {**receipt, "stationary": False, "owner_retaining": False}, None
        raise RuntimeError("incomplete stage has no independently established outcome")
    endpoint = clusters["endpoint"]
    active = constraints(source, spec)
    structure = inspect_structure(raw, endpoint, source, active)
    saved_structure = receipt.get("structure", {})
    for key, value in structure.items():
        if key in ("minimum_pair_distance_a",):
            if not math.isclose(
                value, saved_structure.get(key, math.inf), abs_tol=1e-12
            ):
                raise RuntimeError("structural receipt mismatch")
        elif saved_structure.get(key) != value:
            raise RuntimeError(f"raw/structural receipt mismatch: {key}")
    owners = _owner_labels(endpoint)
    if receipt.get("observed_owners") != owners or receipt.get(
        "reference_owners"
    ) != _owner_labels(source.cluster):
        raise RuntimeError("owner receipt mismatch")
    fresh = reproduce(endpoint, spec.settings, receipt.get("endpoint_evidence", {}))
    metrics = independent_metrics(endpoint, fresh["gradient_hartree_per_bohr"], active)
    saved = receipt.get("evidence", {})
    for key, value in metrics.items():
        if not np.allclose(value, saved.get(key), atol=1e-12, rtol=0):
            raise RuntimeError("projected gradient receipt mismatch")
    if (
        abs(fresh["energy_hartree"] - saved.get("energy_hartree", math.inf))
        > a3g.ENERGY_REPRO_TOL
    ):
        raise RuntimeError("endpoint energy receipt mismatch")
    verified = {
        **receipt,
        "structure": structure,
        "owner_changes": structure["owner_changes"],
        "observed_owners": owners,
        "owner_retaining": structure["owner_changes"] == [],
        "evidence": {**metrics, "energy_hartree": fresh["energy_hartree"]},
        "stationary": optimizer.get("converged") is True and metrics["passed"],
    }
    if spec.fully_released and verified["stationary"]:
        verified["phva"] = phva(endpoint, spec.settings, fresh["energy_hartree"])
        saved_phva = receipt.get("phva", {})
        for key, value in verified["phva"].items():
            prior = saved_phva.get(key)
            if key == "imaginary_cm":
                if (
                    not isinstance(prior, list)
                    or len(prior) != len(value)
                    or not np.allclose(prior, value, atol=1e-3, rtol=0)
                ):
                    raise RuntimeError("PHVA imaginary modes mismatch")
            elif key == "electronic_hartree":
                if (
                    not isinstance(prior, (int, float))
                    or not math.isfinite(prior)
                    or abs(prior - value) > a3g.ENERGY_REPRO_TOL
                ):
                    raise RuntimeError("PHVA energy mismatch")
            elif prior != value:
                raise RuntimeError(f"PHVA receipt mismatch: {key}")
    return verified, endpoint


def constrained_release_allowed(receipt: dict) -> bool:
    if (
        receipt.get("status") != "complete"
        or not receipt.get("stationary")
        or receipt.get("owner_changes") != []
    ):
        return False
    residuals = receipt.get("structure", {}).get("constraint_residuals", [])
    return (
        len(residuals) == 1
        and residuals[0]["oxygen"] == 15
        and residuals[0]["hydrogen"] == 52
        and residuals[0]["absolute_residual_a"] <= 1e-4
        and receipt["evidence"]["passed"]
        and _finite_energy(receipt) is not None
    )


def classify(constrained: dict, released: dict) -> str:
    if not constrained_release_allowed(constrained):
        return a3g.INCONCLUSIVE_CANDIDATE
    if (
        released.get("status") != "complete"
        or not released.get("stationary")
        or not released.get("evidence", {}).get("passed")
        or released.get("phva", {}).get("status") != "passed"
        or _finite_energy(released) is None
    ):
        return a3g.INCONCLUSIVE_CANDIDATE
    before = constrained["observed_owners"]
    after = released["observed_owners"]
    if before == after and released["owner_changes"] == []:
        return a3g.ACCEPTED_CANDIDATE
    expected = ["H52:O9" if owner == "H52:O15" else owner for owner in before]
    if (
        after == expected
        and released["owner_changes"] == ["H52:O15->O9"]
        and released["evidence"]["energy_hartree"]
        < constrained["evidence"]["energy_hartree"] - a3g.DOWNHILL_MIN_HARTREE
    ):
        return a3g.ALTERNATIVE_CANDIDATE
    return a3g.INCONCLUSIVE_CANDIDATE


def verify_experiment(
    output_root: Path,
    *,
    verifier_identity: str,
    source_root: Path = a3g.DEFAULT_SOURCE_ROOT,
    source_override: a3g.SourceEvidence | None = None,
    recompute_calculators: bool = True,
) -> dict[str, Any]:
    """Rehash source and independently derive the A3g terminal classification."""
    output_root = output_root.resolve()
    candidate: dict[str, Any] = {}
    with exclusive_run(output_root):
        try:
            candidate_path = output_root / "candidate-terminal.json"
            candidate = json.loads(candidate_path.read_text())
            if candidate.get("schema") != "a3g-candidate-terminal-v1":
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
            source = source_override or validate_source(source_root)
            if candidate.get("source") != source_map(source):
                raise RuntimeError("candidate source binding mismatch")
            code_revision = candidate.get("code_revision")
            if not isinstance(code_revision, str) or len(code_revision) != 40:
                raise RuntimeError("candidate code revision is invalid")

            receipts: dict[str, dict[str, Any]] = {}
            endpoints: dict[str, Any | None] = {}
            parent = {"stage": None, "receipt_sha256": None}
            for spec in a3g.STAGES:
                stage_state = candidate.get("stages", {}).get(spec.stage_id, {})
                if stage_state.get("status") == "not-run":
                    stage_dir = output_root / "stages" / spec.directory
                    if stage_dir.exists():
                        raise RuntimeError(
                            f"{spec.stage_id} is marked not-run but has artifacts"
                        )
                    receipts[spec.stage_id] = {
                        "schema": "a3g-stage-receipt-v1",
                        "status": "not-run",
                        "stage": spec.stage_id,
                    }
                    endpoints[spec.stage_id] = None
                else:
                    receipt_path = (
                        output_root / "stages" / spec.directory / "receipt.json"
                    )
                    if stage_state.get("receipt_sha256") != sha256_path(receipt_path):
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
                        expected_seed=(
                            source.cluster
                            if spec.parent is None
                            else endpoints[spec.parent]
                        ),
                    )
                    if stage_state.get("status") != receipt.get("status"):
                        raise RuntimeError("candidate stage status mismatch")
                    receipts[spec.stage_id] = receipt
                    endpoints[spec.stage_id] = endpoint
                if spec.parent is None:
                    parent = {
                        "stage": spec.stage_id,
                        "receipt_sha256": (
                            sha256_path(
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
                        "receipt_sha256": sha256_path(
                            output_root / "stages" / spec.directory / "receipt.json"
                        ),
                    }

            conditioning = receipts[a3g.STAGES[0].stage_id]
            constrained = receipts[a3g.STAGES[1].stage_id]
            released = receipts[a3g.STAGES[2].stage_id]
            if candidate.get("experiment_budget") != _experiment_budget(receipts):
                raise RuntimeError("candidate experiment_budget mismatch")
            if constrained.get("status") != "not-run" and not (
                conditioning.get("status") == "complete"
                and conditioning.get("owner_retaining") is True
                and conditioning.get("owner_changes") == []
                and _finite_energy(conditioning) is not None
            ):
                raise RuntimeError(
                    "constrained stage ran without valid conditioning parent"
                )
            if released.get("status") != "not-run" and not constrained_release_allowed(
                constrained
            ):
                raise RuntimeError(
                    "released stage ran without a valid constrained release gate"
                )
            classification = classify(constrained, released)
            if classification != candidate.get("classification"):
                raise RuntimeError(
                    "candidate classification mismatch: "
                    f"claimed {candidate.get('classification')!r}, "
                    f"derived {classification!r}"
                )
            verified_classification = {
                a3g.ACCEPTED_CANDIDATE: a3g.VERIFIED_ACCEPTED,
                a3g.ALTERNATIVE_CANDIDATE: a3g.VERIFIED_ALTERNATIVE,
                a3g.INCONCLUSIVE_CANDIDATE: a3g.VERIFIED_INCONCLUSIVE,
            }[classification]
            result = {
                "schema": "a3g-verified-terminal-v1",
                "written_at": now(),
                "status": "verified",
                "classification": verified_classification,
                "executor_identity": executor_identity,
                "verifier_identity": verifier_identity,
                "candidate_terminal_sha256": sha256_path(candidate_path),
                "source_rehashed": source_override is None,
                "calculator_evidence_recomputed": True,
                "code_revision": code_revision,
                "stage_receipt_sha256": {
                    spec.stage_id: candidate.get("stages", {})
                    .get(spec.stage_id, {})
                    .get("receipt_sha256")
                    for spec in a3g.STAGES
                },
            }
        except Exception as exc:
            revoke_terminals(output_root)
            result = {
                "schema": "a3g-verified-terminal-v1",
                "written_at": now(),
                "status": "rejected",
                "classification": a3g.VERIFIED_INCONCLUSIVE,
                "executor_identity": candidate.get("executor_identity"),
                "verifier_identity": verifier_identity,
                "detail": f"{type(exc).__name__}: {exc}",
            }
        atomic_json(output_root / "verified-terminal.json", result)
        return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=a3g.DEFAULT_SOURCE_ROOT)
    parser.add_argument("--output-root", type=Path, default=a3g.DEFAULT_OUTPUT_ROOT)
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
        preload_cutensor()
    result = verify_experiment(
        args.output_root,
        source_root=args.source_root,
        verifier_identity=args.verifier_identity,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "verified" else 1


if __name__ == "__main__":
    raise SystemExit(main())
