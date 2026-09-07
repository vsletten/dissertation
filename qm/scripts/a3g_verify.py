#!/usr/bin/env python3
"""Cold optimizer-free verifier for the bounded A3g H52 owner-basin run."""

from __future__ import annotations

import argparse
import ast
import fcntl
import hashlib
import json
import math
import os
import re
import socket
import subprocess
import sys
import time
import uuid
from contextlib import contextmanager
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from types import CodeType
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
    DftSettings,
    energy,
    frequencies,
    frequency_geometry_fingerprint,
    frequency_settings_fingerprint,
    gradient,
)
from quarry.store import geometry_hash  # noqa: E402
from scripts import a3g_oaa_owner_basin as a3g  # noqa: E402
from scripts.phase2_ladder import preload_cutensor  # noqa: E402

AUDITED_EXECUTOR_SOURCES = {
    "qm/scripts/a3g_oaa_owner_basin.py": (
        "ab594073b47caa5b35985b7899b9a8541260d73241222b00b817daa5e3919261"
    ),
    "qm/quarry/pipeline.py": (
        "6d92a01e1d275c2f97926731cc61e295eba360ce2e8850ac48a0a2e0ded6297f"
    ),
}

# Independently maintained card contract; never read executor runtime constants.
DEFAULT_SOURCE_ROOT = Path(
    "/mnt/data/vsletten/dissertation-data/task274-a3-oaa-neutral-family-20260906"
)
DEFAULT_OUTPUT_ROOT = Path(
    "/mnt/data/vsletten/dissertation-data/a3g-oaa-neutral-n2-proton-microstate-stability"
)
SOURCE_CELL = Path("runs/phase2/oaa-neutral-n2-s2-b3lyp-def2-svp")
SOURCE_PATHS = {
    "family_receipt": Path("terminal-receipt.json"),
    "family_progress": Path("family-progress.json"),
    "child_log": Path("logs/oaa-neutral-n2.log"),
    "seed": SOURCE_CELL / "complex_guess.xyz",
    "metadata": SOURCE_CELL / "metadata.json",
    "launch": Path("launch-receipt.txt"),
    "restoration": Path("restoration-receipt.txt"),
}
EXPECTED_SOURCE_SHA256 = {
    "family_receipt": (
        "e6cfb7063d1dc6bdfe9ce88324183f54838560aa379673a3c1ca2ce20fb531e6"
    ),
    "family_progress": (
        "2d75c859be5fe2496432537bfb7e9892b5ccf228377c762caad98d050457a92b"
    ),
    "child_log": "ff607898882e7f16c1fd6d52f1283c01646ab32e9dccda603bc80ef54af746a2",
    "seed": "bcf597e78dede88afe0f2b7af360c71ebb90b1fbacbb476eeb02ba6365dfbedd",
    "metadata": "6fdd24d73fb1d2cd5dd6145345562c5803abcf13e7286d57ae5bf6865a2a6190",
    "launch": "48bb6e018372caf56828d37aeb18915b2a579a58eafd62a6e75bb2bd6d5821a2",
    "restoration": "ad5a78f88fc84be7ba794b274af431d35b1d97da8cab855fccc054d9755d9646",
}
EXPECTED_EXECUTION_SOURCE = "97cea5b585c95f18e499f53317ca8864a4c542d2"
EXPECTED_BRANCH = "agents/A3g-oaa-neutral-n2-proton-microstate-stability"
H52 = 52
O15 = 15
O9 = 9
MAX_STEPS = 100
DOWNHILL_MIN_HARTREE = 1.0e-6
EXECUTOR_IDENTITY = "hermes-custom-build-001"
ACCEPTED_CANDIDATE = "unverified: accepted reactant minimum"
ALTERNATIVE_CANDIDATE = "unverified: lower-energy alternative microstate candidate"
INCONCLUSIVE_CANDIDATE = "unverified: inconclusive terminal failure"
VERIFIED_ACCEPTED = "accepted reactant minimum"
VERIFIED_ALTERNATIVE = "verified lower-energy alternative microstate"
VERIFIED_INCONCLUSIVE = "inconclusive terminal failure"


NOISE_FLOOR_CM = 30.0
ENERGY_REPRO_TOL = 1e-7
GRADIENT_REPRO_TOL = 1e-7


@dataclass(frozen=True)
class StageSpec:
    stage_id: str
    directory: str
    method: str
    settings: DftSettings
    parent: str | None
    constrain_h52: bool
    require_production_release_gate: bool
    fully_released: bool


STAGES = (
    StageSpec(
        "owner-conditioning",
        "01-owner-conditioning",
        "hf/sto-3g",
        DftSettings(xc="hf", basis="sto-3g", use_gpu=False),
        None,
        True,
        False,
        False,
    ),
    StageSpec(
        "constrained-production",
        "02-constrained-production",
        "b3lyp/def2-svp/df",
        DftSettings(xc="b3lyp", basis="def2-svp", density_fit=True, use_gpu=True),
        "owner-conditioning",
        True,
        True,
        False,
    ),
    StageSpec(
        "released-production",
        "03-released-production",
        "b3lyp/def2-svp/df",
        DftSettings(xc="b3lyp", basis="def2-svp", density_fit=True, use_gpu=True),
        "constrained-production",
        False,
        False,
        True,
    ),
)
STAGE_BY_ID = {stage.stage_id: stage for stage in STAGES}


from quarry.durable_receipts import durable_mkdir, durable_replace  # noqa: E402

# Card-scoped downstream artifacts, including family publication and store sidecars.
FORBIDDEN_ARTIFACT_NAMES = frozenset(
    {
        "results.json",
        "store.sqlite",
        "store.sqlite-wal",
        "store.sqlite-shm",
        "store.sqlite-journal",
        "store.task168.tmp.sqlite",
        "store.sequential.tmp.sqlite",
        "ts.xyz",
        "barrier.json",
        "petra.toml",
        "family-progress.json",
        "terminal-receipt.json",
        "terminal.json",
    }
)


def forbidden_inventory(output_root: Path) -> list[str]:
    # No subtree exemptions: evidence/attempt directories cannot hide publication.
    return sorted(
        str(path.relative_to(output_root))
        for path in output_root.rglob("*")
        if path.name in FORBIDDEN_ARTIFACT_NAMES
    )


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    durable_mkdir(path.parent)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    with temporary.open("rb") as stream:
        os.fsync(stream.fileno())
    durable_replace(temporary, path)


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
    expected = expected_hashes or EXPECTED_SOURCE_SHA256
    if set(expected) != set(SOURCE_PATHS):
        raise RuntimeError("source expected-hash set is incomplete")
    observed = {
        name: sha256_path(source_root / relative)
        for name, relative in SOURCE_PATHS.items()
    }
    for name, digest in expected.items():
        if observed[name] != digest:
            raise RuntimeError(f"source {name} SHA-256 mismatch: {observed[name]}")

    family = json.loads((source_root / SOURCE_PATHS["family_receipt"]).read_text())
    if (
        family.get("schema") != "a3-family-campaign-terminal-v1"
        or family.get("success") is not False
        or family.get("family") != "oaa"
        or family.get("state") != "neutral"
        or family.get("current_n_intact") != 2
        or family.get("completed") != []
        or family.get("expected_git_sha") != EXPECTED_EXECUTION_SOURCE
        or family.get("observed_git_sha") != EXPECTED_EXECUTION_SOURCE
    ):
        raise RuntimeError("family terminal identity/outcome/source mismatch")

    metadata = json.loads((source_root / SOURCE_PATHS["metadata"]).read_text())
    if (
        metadata.get("site_kind") != "Oaa"
        or metadata.get("state") != "neutral"
        or metadata.get("n_intact") != 2
        or metadata.get("center_site") != 18
        or metadata.get("metal_shells") != 2
        or metadata.get("charge") != 0
        or metadata.get("method") != "b3lyp/def2-svp/df"
        or metadata.get("driver_git_commit") != EXPECTED_EXECUTION_SOURCE
    ):
        raise RuntimeError("source metadata identity/settings mismatch")

    root = (repo_root or Path(__file__).resolve().parents[2]).resolve()
    template = canonical_template(root)
    seed = read_cluster(source_root / SOURCE_PATHS["seed"], template)
    if (
        seed.symbols != template.symbols
        or seed.charge != template.charge
        or seed.spin != template.spin
        or seed.frozen_indices != template.frozen_indices
        or len(seed.symbols) != 58
    ):
        raise RuntimeError("source seed identity/order/state/frozen-shell mismatch")
    owners = oxygen_proton_owners(seed)
    if owners.get(H52) != O15:
        raise RuntimeError(f"source seed H52 owner is O{owners.get(H52)}, not O15")
    if metadata.get("n_atoms") + 3 != len(seed.symbols):
        raise RuntimeError("metadata cell atom count does not bind attacked complex")
    return a3g.SourceEvidence(seed, source_root, observed, metadata)


def source_map(source: a3g.SourceEvidence) -> dict[str, Any]:
    return {
        "root": str(source.source_root),
        "hashes": dict(sorted(source.hashes.items())),
        "execution_source": EXPECTED_EXECUTION_SOURCE,
        "branch": EXPECTED_BRANCH,
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
    durable_mkdir(lock_path.parent)
    with lock_path.open("a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("A3g experiment is already active") from exc
        try:
            yield
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def stage_signature(
    source: a3g.SourceEvidence,
    spec: StageSpec,
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
        "budget": {"max_steps": MAX_STEPS, "retry_allowed": False},
        "parent": spec.parent,
        "fresh_optimizer_requested": True,
        "default_hessian_requested": True,
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


def constraints(source: a3g.SourceEvidence, spec: StageSpec) -> list:
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


def reproduce(cluster: Cluster, settings: Any, saved: dict, observations: list) -> dict:
    value = float(energy(cluster, settings))
    observations.append("energy")
    full = np.asarray(gradient(cluster, settings), dtype=float)
    observations.append("gradient")
    prior = np.asarray(saved.get("gradient_hartree_per_bohr"), dtype=float)
    if (
        not math.isfinite(value)
        or full.shape != cluster.coords.shape
        or prior.shape != full.shape
        or not np.isfinite(full).all()
        or not np.isfinite(prior).all()
        or not isinstance(saved.get("energy_hartree"), (int, float))
        or not math.isfinite(saved["energy_hartree"])
        or abs(value - saved["energy_hartree"]) > ENERGY_REPRO_TOL
        or not np.allclose(full, prior, atol=GRADIENT_REPRO_TOL, rtol=0)
        or saved.get("settings") != asdict(settings)
        or saved.get("geometry_fingerprint") != frequency_geometry_fingerprint(cluster)
    ):
        raise RuntimeError(
            "raw/endpoint energy, full gradient or settings not reproducible"
        )
    return {"energy_hartree": value, "gradient_hartree_per_bohr": full.tolist()}


def phva(cluster: Cluster, settings: Any, value: float, observations: list) -> dict:
    result = frequencies(cluster, settings)
    observations.append("phva")
    imaginary = [float(v) for v in result.imaginary_cm]
    if (
        not math.isfinite(float(result.electronic_hartree))
        or not all(math.isfinite(v) and v >= 0 for v in imaginary)
        or abs(float(result.electronic_hartree) - value) > ENERGY_REPRO_TOL
        or result.geometry_fingerprint != frequency_geometry_fingerprint(cluster)
        or result.settings_fingerprint != frequency_settings_fingerprint(settings)
    ):
        raise RuntimeError("PHVA evidence mismatch")
    return {
        "status": "passed" if max(imaginary, default=0) <= NOISE_FLOOR_CM else "failed",
        "imaginary_cm": imaginary,
        "electronic_hartree": float(result.electronic_hartree),
        "fresh_hessian_requested": True,
        "noise_floor_cm": NOISE_FLOOR_CM,
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
    process,
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
    observation = json.loads((stage_dir / "optimizer-observation.json").read_text())
    if (
        optimizer.get("fresh_optimizer_requested") is not True
        or optimizer.get("default_hessian_requested") is not True
        or optimizer.get("observed_calls") != 1
        or optimizer.get("observed_retries") != 0
        or optimizer.get("requested_max_steps") != 100
        or optimizer.get("runtime_observation") != observation
        or observation.get("calls") != 1
        or observation.get("kernel_calls") not in (0, 1)
        or observation.get("pid") != process.get("pid")
        or observation.get("hostname") != process.get("hostname")
    ):
        raise RuntimeError("optimizer budget/process observation mismatch")
    if receipt.get("status") == "optimizer-failed":
        if (stage_dir / "raw-endpoint.xyz").exists():
            raise RuntimeError("optimizer failure has raw endpoint")
        return receipt, {}
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
    for prefix, key in (("raw", "raw_evidence"), ("endpoint", "endpoint_evidence")):
        if key not in receipt:
            continue
        saved = receipt[key]
        for quantity in ("energy_hartree", "gradient_hartree_per_bohr"):
            suffix = "energy" if quantity == "energy_hartree" else "gradient"
            item = json.loads((stage_dir / f"{prefix}-{suffix}.json").read_text())
            if item != {
                k: saved[k] for k in ("settings", "geometry_fingerprint", quantity)
            }:
                raise RuntimeError("individual raw calculator receipt mismatch")
    return receipt, clusters


def _recompute_stage(output_root, source, spec, receipt, clusters, observations):
    stage_dir = output_root / "stages" / spec.directory
    if receipt.get("status") == "optimizer-failed":
        return receipt, None
    raw = clusters["raw_endpoint"]
    reproduce(raw, spec.settings, receipt.get("raw_evidence", {}), observations)
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
    fresh = reproduce(
        endpoint, spec.settings, receipt.get("endpoint_evidence", {}), observations
    )
    metrics = independent_metrics(endpoint, fresh["gradient_hartree_per_bohr"], active)
    saved = receipt.get("evidence", {})
    for key, value in metrics.items():
        if not np.allclose(value, saved.get(key), atol=1e-12, rtol=0):
            raise RuntimeError("projected gradient receipt mismatch")
    if (
        abs(fresh["energy_hartree"] - saved.get("energy_hartree", math.inf))
        > ENERGY_REPRO_TOL
    ):
        raise RuntimeError("endpoint energy receipt mismatch")
    verified = {
        **receipt,
        "structure": structure,
        "owner_changes": structure["owner_changes"],
        "observed_owners": owners,
        "owner_retaining": structure["owner_changes"] == [],
        "evidence": {**metrics, "energy_hartree": fresh["energy_hartree"]},
        "stationary": receipt.get("optimizer", {}).get("converged") is True
        and metrics["passed"],
    }
    if spec.fully_released and verified["stationary"]:
        verified["phva"] = phva(
            endpoint, spec.settings, fresh["energy_hartree"], observations
        )
        saved_phva = receipt.get("phva", {})
        raw_phva = json.loads((stage_dir / "raw-phva.json").read_text())
        if receipt.get("raw_phva_sha256") != sha256_path(stage_dir / "raw-phva.json"):
            raise RuntimeError("raw PHVA hash mismatch")
        output = raw_phva.pop("output", None)
        if not isinstance(output, dict) or any(
            output.get(key) != raw_phva.get(key)
            for key in (
                "imaginary_cm",
                "electronic_hartree",
                "geometry_fingerprint",
                "settings_fingerprint",
            )
        ):
            raise RuntimeError("raw PHVA output mismatch")
        if raw_phva != {
            "settings": asdict(spec.settings),
            "input_geometry_fingerprint": frequency_geometry_fingerprint(endpoint),
            **{
                key: saved_phva.get(key)
                for key in (
                    "imaginary_cm",
                    "electronic_hartree",
                    "geometry_fingerprint",
                    "settings_fingerprint",
                )
            },
        }:
            raise RuntimeError("raw PHVA receipt mismatch")
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
                    or abs(prior - value) > ENERGY_REPRO_TOL
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
        return INCONCLUSIVE_CANDIDATE
    if (
        released.get("status") != "complete"
        or not released.get("stationary")
        or not released.get("evidence", {}).get("passed")
        or released.get("phva", {}).get("status") != "passed"
        or _finite_energy(released) is None
    ):
        return INCONCLUSIVE_CANDIDATE
    before = constrained["observed_owners"]
    after = released["observed_owners"]
    if before == after and released["owner_changes"] == []:
        return ACCEPTED_CANDIDATE
    expected = ["H52:O9" if owner == "H52:O15" else owner for owner in before]
    if (
        after == expected
        and released["owner_changes"] == ["H52:O15->O9"]
        and released["evidence"]["energy_hartree"]
        < constrained["evidence"]["energy_hartree"] - DOWNHILL_MIN_HARTREE
    ):
        return ALTERNATIVE_CANDIDATE
    return INCONCLUSIVE_CANDIDATE


def runtime_provenance():
    """Observe the verifier process locally; caller labels cannot supply these facts."""
    return {
        "pid": os.getpid(),
        "hostname": socket.gethostname(),
        "boot_id": Path("/proc/sys/kernel/random/boot_id").read_text().strip(),
        "process_start_ticks": Path("/proc/self/stat")
        .read_text()
        .rsplit(")", 1)[1]
        .split()[19],
        "executable_sha256": sha256_path(Path("/proc/self/exe")),
        "implementation_sha256": sha256_path(Path(__file__)),
    }


def validate_identity(identity, executor, runtime=None):
    runtime = runtime if runtime is not None else runtime_provenance()
    if not isinstance(executor, str) or not executor.strip():
        raise RuntimeError("missing executor identity")
    allowed = {"worker_id", "profile", "hostname", "implementation_sha256"}
    if (
        not isinstance(identity, dict)
        or set(identity) - allowed
        or any(
            not isinstance(identity.get(key), str) or not identity[key].strip()
            for key in ("worker_id", "profile")
        )
    ):
        raise RuntimeError(
            "structured verifier identity receipt required; runtime facts are internal"
        )
    if identity["worker_id"] == executor:
        raise RuntimeError("verifier identity must differ from executor identity")
    for key in ("hostname", "implementation_sha256"):
        if key in identity and identity[key] != runtime[key]:
            raise RuntimeError("verifier identity host/implementation mismatch")


def validate_process_independence(runtime, executor):
    keys = ("hostname", "boot_id", "pid", "process_start_ticks")
    if any(not executor.get(key) for key in keys):
        raise RuntimeError("missing executor process identity")
    if all(str(runtime[key]) == str(executor[key]) for key in keys):
        raise RuntimeError("same-process verification is forbidden")


def rejection_attempt(output_root, detail, observations=()):
    result = {
        "schema": "a3g-verification-attempt-v1",
        "written_at": now(),
        "status": "rejected",
        "classification": VERIFIED_INCONCLUSIVE,
        "detail": detail,
        "calculator_evidence_recomputed": bool(observations),
        "calculator_evidence_recomputed_count": len(observations),
    }
    atomic_json(
        output_root / "verification-attempts" / f"{uuid.uuid4().hex}.json", result
    )
    return result


def verify_experiment(output_root: Path, **kwargs) -> dict[str, Any]:
    """Reject invalid/contending attempts without modifying canonical evidence."""
    try:
        return _verify_experiment(output_root, **kwargs)
    except Exception as exc:
        return rejection_attempt(output_root.resolve(), f"{type(exc).__name__}: {exc}")


def compiled_code_digest(code):
    # Independent normalization of executable instructions; paths/line locations
    # are excluded so a separate worktree can check the recorded loaded code.
    fields = {
        "bytecode": code.co_code.hex(),
        "constants": [
            compiled_code_digest(c) if isinstance(c, CodeType) else repr(c)
            for c in code.co_consts
        ],
        "names": code.co_names,
        "variables": code.co_varnames,
        "freevars": code.co_freevars,
        "cellvars": code.co_cellvars,
        "flags": code.co_flags,
        "arguments": [
            code.co_argcount,
            code.co_posonlyargcount,
            code.co_kwonlyargcount,
        ],
    }
    return hashlib.sha256(json.dumps(fields, sort_keys=True).encode()).hexdigest()


def validate_revision(worktree, revision, provenance):
    def git(*args):
        return subprocess.run(
            ["git", "-C", str(worktree), *args],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()

    if (
        git("symbolic-ref", "--short", "HEAD") != EXPECTED_BRANCH
        or git("rev-parse", "HEAD") != revision
        or git("rev-parse", f"origin/{EXPECTED_BRANCH}") != revision
        or git("status", "--porcelain")
    ):
        raise RuntimeError(
            "verification requires exact clean local/remote worktree revision"
        )
    verifier_path = worktree / "qm/scripts/a3g_verify.py"
    if sha256_path(verifier_path) != sha256_path(Path(__file__)):
        raise RuntimeError(
            "running verifier is not the selected worktree implementation"
        )
    if (
        not isinstance(provenance.get("pid"), int)
        or provenance["pid"] <= 0
        or not provenance.get("hostname")
    ):
        raise RuntimeError("missing executor process provenance")
    if (
        provenance.get("executable_sha256") != sha256_path(Path("/proc/self/exe"))
        or not provenance.get("boot_id")
        or not str(provenance.get("process_start_ticks", "")).isdigit()
    ):
        raise RuntimeError("missing executable/process-start provenance")
    # Audited complete sources: one direct call in run_stage; pipeline constructs
    # a new SCF and calls one mutually exclusive kernel branch, with no retry or
    # Hessian input. Hash pins deliberately require a new audit after any change.
    for function, filename in (
        ("optimize_bounded", "qm/quarry/pipeline.py"),
        ("run_stage", "qm/scripts/a3g_oaa_owner_basin.py"),
    ):
        code = (worktree / filename).read_text()
        node = next(
            n
            for n in ast.parse(code).body
            if isinstance(n, ast.FunctionDef) and n.name == function
        )
        digest = hashlib.sha256(
            (ast.get_source_segment(code, node) + "\n").encode()
        ).hexdigest()
        compiled = compile(code, str(worktree / filename), "exec", dont_inherit=True)
        function_code = next(
            c
            for c in compiled.co_consts
            if isinstance(c, CodeType) and c.co_name == function
        )
        if provenance.get("loaded_code", {}).get(function) != compiled_code_digest(
            function_code
        ):
            raise RuntimeError("loaded executable code differs from audited source")
        if provenance.get("callables", {}).get(function) != digest:
            raise RuntimeError("executed callable differs from audited source")
    for name, digest in AUDITED_EXECUTOR_SOURCES.items():
        path = worktree / name
        if (
            sha256_path(path) != digest
            or hashlib.sha256(
                subprocess.run(
                    ["git", "-C", str(worktree), "show", f"{revision}:{name}"],
                    check=True,
                    capture_output=True,
                ).stdout
            ).hexdigest()
            != digest
            or provenance.get("files", {}).get(name) != digest
        ):
            raise RuntimeError("executor differs from audited one-shot implementation")


def _prepare_stages(output_root, source, candidate, code_revision, reservation):
    """Read and bind every stage receipt and endpoint before any calculator call."""
    prepared = {}
    endpoints: dict[str, Any | None] = {}
    parent = {"stage": None, "receipt_sha256": None}
    for spec in STAGES:
        stage_state = candidate.get("stages", {}).get(spec.stage_id, {})
        if spec.parent is None and stage_state.get("status") == "not-run":
            raise RuntimeError("required first stage was not attempted")
        if stage_state.get("status") == "not-run":
            stage_dir = output_root / "stages" / spec.directory
            if stage_dir.exists():
                raise RuntimeError(
                    f"{spec.stage_id} is marked not-run but has artifacts"
                )
            receipt = {
                "schema": "a3g-stage-receipt-v1",
                "status": "not-run",
                "stage": spec.stage_id,
            }
            clusters = {}
            endpoints[spec.stage_id] = None
        else:
            receipt_path = output_root / "stages" / spec.directory / "receipt.json"
            if stage_state.get("receipt_sha256") != sha256_path(receipt_path):
                raise RuntimeError(f"{spec.stage_id} candidate receipt hash mismatch")
            receipt, clusters = _read_stage(
                output_root,
                source,
                spec,
                code_revision,
                parent,
                process=reservation["provenance"],
                expected_seed=(
                    source.cluster if spec.parent is None else endpoints[spec.parent]
                ),
            )
            if stage_state.get("status") != receipt.get("status"):
                raise RuntimeError("candidate stage status mismatch")
            endpoints[spec.stage_id] = clusters.get("endpoint")
        prepared[spec.stage_id] = (receipt, clusters)
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
        elif (output_root / "stages" / spec.directory / "receipt.json").is_file():
            parent = {
                "stage": spec.stage_id,
                "receipt_sha256": sha256_path(
                    output_root / "stages" / spec.directory / "receipt.json"
                ),
            }

    if candidate.get("experiment_budget") != _experiment_budget(
        {stage_id: receipt for stage_id, (receipt, _) in prepared.items()}
    ):
        raise RuntimeError("candidate experiment_budget mismatch")
    return prepared


@dataclass(frozen=True)
class VerificationPreflight:
    candidate: dict[str, Any]
    verifier_identity: dict[str, Any]
    runtime: dict[str, Any]
    source: a3g.SourceEvidence
    source_rehashed: bool
    code_revision: str
    reservation: dict[str, Any]
    stages: dict[str, tuple[dict[str, Any], dict[str, Cluster]]]


def verification_preflight(
    output_root: Path,
    *,
    verifier_identity: dict | Path,
    worktree: Path,
    source_root: Path,
    source_override: a3g.SourceEvidence | None = None,
    recompute_calculators: bool = True,
) -> VerificationPreflight:
    """Validate calculator-free prerequisites under the caller's run lock."""
    if (output_root / "verified-terminal.json").exists():
        raise RuntimeError("replay attempt: verified terminal already exists")
    if isinstance(verifier_identity, Path):
        verifier_identity = json.loads(verifier_identity.read_text())
    runtime = runtime_provenance()
    inventory = forbidden_inventory(output_root)
    if inventory:
        raise RuntimeError(f"forbidden artifacts present: {inventory}")
    candidate_path = output_root / "candidate-terminal.json"
    candidate = json.loads(candidate_path.read_text())
    if candidate.get("schema") != "a3g-candidate-terminal-v1":
        raise RuntimeError("candidate terminal schema mismatch")
    executor_identity = candidate.get("executor_identity")
    validate_identity(verifier_identity, executor_identity, runtime)
    if (
        candidate.get("forbidden_artifacts") != []
        or candidate.get("forbidden_outputs_emitted") is not False
    ):
        raise RuntimeError("candidate forbidden artifact inventory mismatch")
    if not recompute_calculators:
        raise RuntimeError("verification requires calculator recomputation")
    source = source_override or validate_source(source_root, repo_root=worktree)
    if candidate.get("source") != source_map(source):
        raise RuntimeError("candidate source binding mismatch")
    code_revision = candidate.get("code_revision")
    if (
        not isinstance(code_revision, str)
        or re.fullmatch("[0-9a-f]{40}", code_revision) is None
    ):
        raise RuntimeError("candidate code revision is invalid")

    marker = output_root / "experiment-reservation.json"
    reservation = json.loads(marker.read_text())
    if (
        candidate.get("experiment_reservation_sha256") != sha256_path(marker)
        or reservation.get("schema") != "a3g-experiment-reservation-v2"
        or reservation.get("source") != source_map(source)
        or reservation.get("code_revision") != code_revision
        or reservation.get("executor_identity") != executor_identity
    ):
        raise RuntimeError("experiment reservation mismatch")
    validate_process_independence(runtime, reservation.get("provenance", {}))
    validate_revision(worktree, code_revision, reservation.get("provenance", {}))
    return VerificationPreflight(
        candidate,
        verifier_identity,
        runtime,
        source,
        source_override is None,
        code_revision,
        reservation,
        _prepare_stages(output_root, source, candidate, code_revision, reservation),
    )


def _verify_experiment(
    output_root: Path,
    *,
    verifier_identity: dict | Path,
    worktree: Path | None = None,
    source_root: Path = DEFAULT_SOURCE_ROOT,
    source_override: a3g.SourceEvidence | None = None,
    recompute_calculators: bool = True,
    initialize_gpu: bool | None = None,
) -> dict[str, Any]:
    """Share CLI/library preflight and retain its snapshot through recomputation."""
    output_root = output_root.resolve()
    worktree = (worktree or Path(__file__).resolve().parents[2]).resolve()
    observations = []
    with exclusive_run(output_root):
        try:
            preflight = verification_preflight(
                output_root,
                verifier_identity=verifier_identity,
                worktree=worktree,
                source_root=source_root,
                source_override=source_override,
                recompute_calculators=recompute_calculators,
            )
            production_ran = (
                preflight.candidate.get("stages", {})
                .get("constrained-production", {})
                .get("status")
                != "not-run"
            )
            # None is the library path, whose caller manages calculator setup.
            if initialize_gpu is False and production_ran:
                raise RuntimeError(
                    "--gpu is required to recompute production calculator evidence"
                )
            if initialize_gpu:
                preload_cutensor()
            result = _recompute_experiment(output_root, preflight, observations)
        except Exception as exc:
            return rejection_attempt(
                output_root, f"{type(exc).__name__}: {exc}", observations
            )
        atomic_json(output_root / "verified-terminal.json", result)
        return result


def _recompute_experiment(
    output_root: Path, preflight: VerificationPreflight, observations: list
) -> dict[str, Any]:
    candidate = preflight.candidate
    candidate_path = output_root / "candidate-terminal.json"
    verifier_identity = preflight.verifier_identity
    executor_identity = candidate["executor_identity"]
    runtime = preflight.runtime
    source = preflight.source
    code_revision = preflight.code_revision
    receipts = {}
    for spec in STAGES:
        receipt, clusters = preflight.stages[spec.stage_id]
        if receipt["status"] != "not-run":
            receipt, _endpoint = _recompute_stage(
                output_root, source, spec, receipt, clusters, observations
            )
        receipts[spec.stage_id] = receipt

    conditioning = receipts[STAGES[0].stage_id]
    constrained = receipts[STAGES[1].stage_id]
    released = receipts[STAGES[2].stage_id]
    if constrained.get("status") != "not-run" and not (
        conditioning.get("status") == "complete"
        and conditioning.get("owner_retaining") is True
        and conditioning.get("owner_changes") == []
        and _finite_energy(conditioning) is not None
    ):
        raise RuntimeError("constrained stage ran without valid conditioning parent")
    if released.get("status") != "not-run" and not constrained_release_allowed(
        constrained
    ):
        raise RuntimeError(
            "released stage ran without a valid constrained release gate"
        )
    if constrained.get("status") == "not-run" and (
        conditioning.get("status") == "complete"
        and conditioning.get("owner_retaining") is True
        and conditioning.get("owner_changes") == []
        and _finite_energy(conditioning) is not None
    ):
        raise RuntimeError("required constrained stage was skipped")
    if released.get("status") == "not-run" and constrained_release_allowed(constrained):
        raise RuntimeError("required released stage was skipped")
    classification = classify(constrained, released)
    if classification != candidate.get("classification"):
        raise RuntimeError(
            "candidate classification mismatch: "
            f"claimed {candidate.get('classification')!r}, "
            f"derived {classification!r}"
        )
    verified_classification = {
        ACCEPTED_CANDIDATE: VERIFIED_ACCEPTED,
        ALTERNATIVE_CANDIDATE: VERIFIED_ALTERNATIVE,
        INCONCLUSIVE_CANDIDATE: VERIFIED_INCONCLUSIVE,
    }[classification]
    if forbidden_inventory(output_root):
        raise RuntimeError("forbidden artifacts appeared during verification")
    result = {
        "schema": "a3g-verified-terminal-v1",
        "written_at": now(),
        "status": "verified",
        "classification": verified_classification,
        "executor_identity": executor_identity,
        "verifier_identity": verifier_identity,
        "verifier_provenance": runtime,
        "forbidden_artifacts": [],
        "candidate_terminal_sha256": sha256_path(candidate_path),
        "source_rehashed": preflight.source_rehashed,
        "calculator_evidence_recomputed": bool(observations),
        "calculator_evidence_recomputed_count": len(observations),
        "code_revision": code_revision,
        "stage_receipt_sha256": {
            spec.stage_id: candidate.get("stages", {})
            .get(spec.stage_id, {})
            .get("receipt_sha256")
            for spec in STAGES
        },
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--verifier-identity",
        type=Path,
        required=True,
        help="JSON worker identity receipt",
    )
    parser.add_argument(
        "--worktree", type=Path, default=Path(__file__).resolve().parents[2]
    )
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
    result = verify_experiment(
        args.output_root,
        source_root=args.source_root,
        verifier_identity=args.verifier_identity,
        worktree=args.worktree,
        initialize_gpu=args.gpu,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "verified" else 1


if __name__ == "__main__":
    raise SystemExit(main())
