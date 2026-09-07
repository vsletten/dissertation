#!/usr/bin/env python3
"""Bounded A3g H52 owner-constrained to released production-basin experiment."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

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

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    bootstrap_cli(
        "a3g_oaa_owner_basin",
        default_run_root=DEFAULT_OUTPUT_ROOT / "logs",
        gpu_owner="a3g_oaa_owner_basin",
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
from scripts import a3b_proton_microstate_stability as a3b  # noqa: E402
from scripts.phase2_ladder import (  # noqa: E402
    NOISE_FLOOR_CM,
    preload_cutensor,
)


@dataclass(frozen=True)
class SourceEvidence:
    cluster: Cluster
    source_root: Path
    hashes: dict[str, str]
    metadata: dict[str, Any]


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


def atomic_xyz(path: Path, cluster: Cluster) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        "\n".join(
            [str(len(cluster.symbols)), cluster.name]
            + [
                f"{symbol} {x:.17g} {y:.17g} {z:.17g}"
                for symbol, (x, y, z) in zip(
                    cluster.symbols, cluster.coords, strict=True
                )
            ]
        )
        + "\n"
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


def oxygen_proton_owners(cluster: Cluster) -> dict[int, int]:
    owners = {}
    oxygens = [i for i, symbol in enumerate(cluster.symbols) if symbol == "O"]
    for h, symbol in enumerate(cluster.symbols):
        if symbol != "H":
            continue
        distances = sorted(
            (float(np.linalg.norm(cluster.coords[h] - cluster.coords[o])), o)
            for o in oxygens
        )
        if (
            not distances
            or not all(math.isfinite(d) for d, _ in distances)
            or sum(d <= 1.25 for d, _ in distances) != 1
            or (len(distances) > 1 and distances[1][0] - distances[0][0] < 0.15)
        ):
            raise RuntimeError(f"H{h} lacks exactly one unambiguous oxygen owner")
        owners[h] = distances[0][1]
    return owners


def _owner_labels(cluster: Cluster) -> list[str]:
    return [f"H{h}:O{o}" for h, o in sorted(oxygen_proton_owners(cluster).items())]


def validate_source(
    source_root: Path,
    *,
    expected_hashes: dict[str, str] | None = None,
    repo_root: Path | None = None,
) -> SourceEvidence:
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
    seed = a3b.read_cluster(source_root / SOURCE_PATHS["seed"], template)
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
    return SourceEvidence(seed, source_root, observed, metadata)


def source_map(source: SourceEvidence) -> dict[str, Any]:
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


def resolve_code_revision(worktree: Path) -> str:
    worktree = worktree.resolve()
    branch = subprocess.run(
        ["git", "-C", str(worktree), "symbolic-ref", "--short", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    head = subprocess.run(
        ["git", "-C", str(worktree), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    remote = subprocess.run(
        ["git", "-C", str(worktree), "rev-parse", f"origin/{branch}"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    dirty = subprocess.run(
        ["git", "-C", str(worktree), "status", "--porcelain"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    if branch != EXPECTED_BRANCH:
        raise RuntimeError(f"wrong branch for A3g execution: {branch}")
    if dirty:
        raise RuntimeError("worktree must be fully clean before production execution")
    if head != remote:
        raise RuntimeError("local and remote implementation revisions differ")
    return head


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


def h52_constraint(reference: Cluster) -> list[tuple[int, int, float]]:
    return [
        (
            O15,
            H52,
            float(np.linalg.norm(reference.coords[O15] - reference.coords[H52])),
        )
    ]


def constraints(
    source: SourceEvidence, spec: StageSpec
) -> list[tuple[int, int, float]]:
    return h52_constraint(source.cluster) if spec.constrain_h52 else []


def stage_signature(
    source: SourceEvidence,
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


def _optimizer_record(
    *,
    converged: bool,
    observed_calls: int,
    observed_retries: int,
    observed_max_steps: int,
) -> dict[str, Any]:
    return {
        "converged": converged,
        "fresh_instance": True,
        "geometric_default_fresh_hessian": True,
        "observed_calls": observed_calls,
        "observed_retries": observed_retries,
        "observed_max_steps": observed_max_steps,
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


def _resume_stage(
    stage_dir: Path,
    expected_signature: dict[str, Any],
    template: Cluster,
    parent: dict[str, Any],
) -> tuple[dict[str, Any], Cluster | None] | None:
    if stage_dir.exists():
        return ({"status": "orphaned-stage", "signature": expected_signature}, None)
    return None


def run_stage(
    output_root: Path,
    source: SourceEvidence,
    spec: StageSpec,
    seed: Cluster,
    code_revision: str,
    parent_receipt: dict[str, Any] | None,
) -> tuple[dict[str, Any], Cluster | None]:
    stage_dir = output_root / "stages" / spec.directory
    signature = stage_signature(source, spec, code_revision)
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
    resumed = _resume_stage(stage_dir, signature, source.cluster, parent)
    if resumed is not None:
        return resumed

    stage_dir.mkdir(parents=True, exist_ok=True)
    reservation = {
        "schema": "a3g-stage-reservation-v1",
        "status": "reserved",
        "reserved_at": now(),
        "signature": signature,
        "parent": parent,
        "seed": {
            "sha256": hashlib.sha256(seed.to_xyz().encode()).hexdigest(),
            "geometry_fingerprint": frequency_geometry_fingerprint(seed),
        },
    }
    atomic_json(stage_dir / "reservation.json", reservation)
    active = constraints(source, spec)
    kwargs: dict[str, Any] = {"max_steps": MAX_STEPS}
    if active:
        kwargs["fixed_distances"] = active
    observed_calls = 0
    observed_retries = 0
    observed_max_steps = int(kwargs["max_steps"])
    try:
        optimized = optimize_bounded(seed, spec.settings, **kwargs)
        observed_calls = 1
    except Exception as exc:
        observed_calls = 1
        receipt = {
            "schema": "a3g-stage-receipt-v1",
            "status": "optimizer-failed",
            "terminal_stage": "optimizer-call",
            "completed_at": now(),
            "stage": spec.stage_id,
            "signature": signature,
            "parent": parent,
            "seed": reservation["seed"],
            "optimizer": _optimizer_record(
                converged=False,
                observed_calls=observed_calls,
                observed_retries=observed_retries,
                observed_max_steps=observed_max_steps,
            ),
            "detail": f"{type(exc).__name__}: {exc}",
        }
        atomic_json(stage_dir / "receipt.json", receipt)
        return receipt, None

    raw_path = stage_dir / "raw-endpoint.xyz"
    atomic_xyz(raw_path, optimized.cluster)
    receipt = {
        "schema": "a3g-stage-receipt-v1",
        "status": "pending-structural-gates",
        "terminal_stage": "structural-gates",
        "stage": spec.stage_id,
        "signature": signature,
        "parent": parent,
        "seed": reservation["seed"],
        "optimizer": _optimizer_record(
            converged=bool(optimized.converged),
            observed_calls=observed_calls,
            observed_retries=observed_retries,
            observed_max_steps=observed_max_steps,
        ),
        "raw_endpoint": _artifact(optimized.cluster, raw_path),
    }
    atomic_json(stage_dir / "receipt.json", receipt)

    try:
        # Persist raw calculator evidence before *any* scientific gate. Rejected
        # geometries remain reproducible; projected geometry gets its own evidence.
        receipt["terminal_stage"] = "raw-calculator-evidence"
        receipt["raw_evidence"] = calculator_evidence(optimized.cluster, spec.settings)
        atomic_json(stage_dir / "receipt.json", receipt)
        receipt["terminal_stage"] = "structural-gates"
        endpoint, raw_structure = a3b.structural_gate(
            optimized.cluster,
            source.cluster,
            active,
            spec.stage_id,
        )
        endpoint_path = stage_dir / "endpoint.xyz"
        atomic_xyz(endpoint_path, endpoint)
        endpoint = a3b.read_cluster(endpoint_path, source.cluster)
        receipt["endpoint"] = _artifact(endpoint, endpoint_path)
        atomic_json(stage_dir / "receipt.json", receipt)
        structure = a3b.structural_gate(
            endpoint,
            source.cluster,
            active,
            spec.stage_id,
        )[1]
        structure["raw_maximum_frozen_coordinate_drift_a"] = raw_structure[
            "raw_maximum_frozen_coordinate_drift_a"
        ]
        receipt["terminal_stage"] = "energy-evidence"
        endpoint_evidence = calculator_evidence(endpoint, spec.settings)
        receipt["endpoint_evidence"] = endpoint_evidence
        atomic_json(stage_dir / "receipt.json", receipt)
        observed_energy = endpoint_evidence["energy_hartree"]
        metrics = gradient_metrics(
            endpoint, endpoint_evidence["gradient_hartree_per_bohr"], active
        )
        oxygen_proton_owners(endpoint)
        owner_retaining = structure["owner_changes"] == []
        stationary = bool(optimized.converged) and bool(metrics["passed"])
        phva: dict[str, Any] = {"status": "not-required", "imaginary_cm": []}
        if spec.fully_released and stationary:
            receipt["terminal_stage"] = "phva-evidence"
            result = frequencies(endpoint, spec.settings)
            imaginary = [float(value) for value in result.imaginary_cm]
            if not math.isfinite(float(result.electronic_hartree)) or not all(
                math.isfinite(value) and value >= 0 for value in imaginary
            ):
                raise RuntimeError("PHVA returned non-finite evidence")
            if (
                result.geometry_fingerprint != frequency_geometry_fingerprint(endpoint)
                or result.settings_fingerprint
                != frequency_settings_fingerprint(spec.settings)
                or abs(float(result.electronic_hartree) - observed_energy)
                > ENERGY_REPRO_TOL
            ):
                raise RuntimeError("PHVA geometry/settings/energy mismatch")
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
                "settings_fingerprint": frequency_settings_fingerprint(spec.settings),
                "electronic_hartree": float(result.electronic_hartree),
            }
        receipt.update(
            {
                "status": "complete",
                "terminal_stage": "stage-complete",
                "completed_at": now(),
                "structure": structure,
                "endpoint": _artifact(endpoint, endpoint_path),
                "reference_owners": _owner_labels(source.cluster),
                "observed_owners": _owner_labels(endpoint),
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
                "completed_at": now(),
                "detail": f"{type(exc).__name__}: {exc}",
                "stationary": False,
            }
        )
        atomic_json(stage_dir / "receipt.json", receipt)
        return receipt, None

    atomic_json(stage_dir / "receipt.json", receipt)
    return receipt, endpoint


ENERGY_REPRO_TOL = 1.0e-7
GRADIENT_REPRO_TOL = 1.0e-7


def calculator_evidence(cluster: Cluster, settings: DftSettings) -> dict[str, Any]:
    value = float(energy(cluster, settings))
    forces = np.asarray(gradient(cluster, settings), dtype=float)
    if (
        not math.isfinite(value)
        or forces.shape != cluster.coords.shape
        or not np.all(np.isfinite(forces))
    ):
        raise RuntimeError("non-finite or malformed energy/full gradient evidence")
    return {
        "energy_hartree": value,
        "gradient_hartree_per_bohr": forces.tolist(),
        "settings": asdict(settings),
        "geometry_fingerprint": frequency_geometry_fingerprint(cluster),
    }


def gradient_metrics(cluster: Cluster, raw: Any, active: list) -> dict[str, Any]:
    projected = a3b._project_gradient(cluster, np.asarray(raw), active)
    free = sorted(set(range(len(cluster.symbols))) - set(cluster.frozen_indices))
    values = projected[free]
    rms = float(np.sqrt(np.mean(values**2))) if values.size else 0.0
    maximum = float(np.max(np.abs(values))) if values.size else 0.0
    return {
        "projected_gradient_hartree_per_bohr": projected.tolist(),
        "projected_free_atom_indices": free,
        "projected_gradient_rms_hartree_per_bohr": rms,
        "projected_gradient_max_hartree_per_bohr": maximum,
        "gradient_rms_threshold_hartree_per_bohr": 3.0e-4,
        "gradient_max_threshold_hartree_per_bohr": 4.5e-4,
        "passed": rms <= 3.0e-4 and maximum <= 4.5e-4,
    }


def _finite_energy(receipt: dict[str, Any]) -> float | None:
    value = receipt.get("evidence", {}).get("energy_hartree")
    if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
        return None
    return float(value)


def constrained_release_allowed(receipt: dict[str, Any]) -> bool:
    residuals = receipt.get("structure", {}).get("constraint_residuals", [])
    return bool(
        receipt.get("status") == "complete"
        and receipt.get("optimizer", {}).get("converged") is True
        and receipt.get("stationary") is True
        and receipt.get("owner_retaining") is True
        and receipt.get("owner_changes") == []
        and receipt.get("evidence", {}).get("passed") is True
        and len(residuals) == 1
        and residuals[0].get("oxygen") == O15
        and residuals[0].get("hydrogen") == H52
        and residuals[0].get("absolute_residual_a", math.inf)
        <= a3b.CONSTRAINT_RESIDUAL_MAX_A
        and _finite_energy(receipt) is not None
    )


def classify(
    constrained: dict[str, Any],
    released: dict[str, Any],
) -> str:
    if not constrained_release_allowed(constrained):
        return INCONCLUSIVE_CANDIDATE
    released_energy = _finite_energy(released)
    constrained_energy = _finite_energy(constrained)
    if (
        released.get("status") == "complete"
        and released.get("stationary") is True
        and released.get("owner_retaining") is True
        and released.get("owner_changes") == []
        and released.get("observed_owners") == constrained.get("observed_owners")
        and released.get("evidence", {}).get("passed") is True
        and released.get("phva", {}).get("status") == "passed"
        and released_energy is not None
    ):
        return ACCEPTED_CANDIDATE
    if (
        released.get("status") == "complete"
        and released.get("stationary") is True
        and released.get("evidence", {}).get("passed") is True
        and released.get("owner_changes") == [f"H{H52}:O{O15}->O{O9}"]
        and released.get("observed_owners")
        == [
            (f"H{H52}:O{O9}" if item == f"H{H52}:O{O15}" else item)
            for item in constrained.get("observed_owners", [])
        ]
        and constrained_energy is not None
        and released_energy is not None
        and released_energy < constrained_energy - DOWNHILL_MIN_HARTREE
        and released.get("phva", {}).get("status") == "passed"
    ):
        return ALTERNATIVE_CANDIDATE
    return INCONCLUSIVE_CANDIDATE


def _not_run(stage: StageSpec, detail: str) -> dict[str, Any]:
    return {
        "schema": "a3g-stage-receipt-v1",
        "status": "not-run",
        "stage": stage.stage_id,
        "detail": detail,
    }


def _run_locked(
    output_root: Path,
    source: SourceEvidence,
    *,
    code_revision: str,
    executor_identity: str,
) -> dict[str, Any]:
    output_root.mkdir(parents=True, exist_ok=True)
    revoke_terminals(output_root)
    receipts: dict[str, dict[str, Any]] = {}
    # The reservation is written before any stage. Even a crash before stage one
    # spends this route; no continuation or orphan replay is authorized.
    marker = output_root / "experiment-reservation.json"
    if marker.exists():
        terminal = {
            "schema": "a3g-candidate-terminal-v1",
            "classification": INCONCLUSIVE_CANDIDATE,
            "detail": "experiment already reserved; zero retries/orphan replay",
            "executor_identity": executor_identity,
            "independent_verification_required": True,
            "forbidden_outputs_emitted": False,
        }
        atomic_json(output_root / "candidate-terminal.json", terminal)
        return terminal
    atomic_json(marker, {"source": source_map(source), "code_revision": code_revision})

    conditioning, conditioned = run_stage(
        output_root, source, STAGES[0], source.cluster, code_revision, None
    )
    receipts[STAGES[0].stage_id] = conditioning
    conditioning_ok = bool(
        conditioned is not None
        and conditioning.get("status") == "complete"
        and conditioning.get("owner_retaining") is True
        and conditioning.get("owner_changes") == []
        and _finite_energy(conditioning) is not None
    )
    if conditioning_ok:
        assert conditioned is not None
        constrained, constrained_endpoint = run_stage(
            output_root,
            source,
            STAGES[1],
            conditioned,
            code_revision,
            conditioning,
        )
    else:
        constrained = _not_run(STAGES[1], "owner-preserving conditioning unavailable")
        constrained_endpoint = None
    receipts[STAGES[1].stage_id] = constrained

    if constrained_endpoint is not None and constrained_release_allowed(constrained):
        released, _released_endpoint = run_stage(
            output_root,
            source,
            STAGES[2],
            constrained_endpoint,
            code_revision,
            constrained,
        )
    else:
        released = _not_run(
            STAGES[2],
            "constrained production endpoint failed the independent release gate",
        )
    receipts[STAGES[2].stage_id] = released
    classification = classify(constrained, released)
    terminal = {
        "schema": "a3g-candidate-terminal-v1",
        "written_at": now(),
        "executor_identity": executor_identity,
        "classification": classification,
        "source": source_map(source),
        "code_revision": code_revision,
        "stages": {
            stage.stage_id: {
                "status": receipts[stage.stage_id].get("status", "unknown"),
                "receipt_sha256": (
                    sha256_path(
                        output_root / "stages" / stage.directory / "receipt.json"
                    )
                    if (
                        output_root / "stages" / stage.directory / "receipt.json"
                    ).is_file()
                    else None
                ),
            }
            for stage in STAGES
        },
        "independent_verification_required": True,
        "forbidden_outputs_emitted": False,
        "experiment_budget": _experiment_budget(receipts),
    }
    if classification == INCONCLUSIVE_CANDIDATE:
        terminal["detail"] = (
            "the single finite constrained-to-released route was inconclusive; "
            "this route and the Oaa-neutral n=2/4/6 family are closed; no replay"
        )
    atomic_json(output_root / "candidate-terminal.json", terminal)
    return terminal


def run_experiment(
    output_root: Path,
    source: SourceEvidence,
    *,
    code_revision: str,
    executor_identity: str = EXECUTOR_IDENTITY,
) -> dict[str, Any]:
    output_root.mkdir(parents=True, exist_ok=True)
    with exclusive_run(output_root):
        return _run_locked(
            output_root,
            source,
            code_revision=code_revision,
            executor_identity=executor_identity,
        )


def dry_run(source: SourceEvidence, *, code_revision: str) -> dict[str, Any]:
    return {
        "schema": "a3g-dry-run-v1",
        "source": source_map(source),
        "code_revision": code_revision,
        "stages": [stage_signature(source, spec, code_revision) for spec in STAGES],
        "production_artifacts_written": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--worktree", type=Path, default=Path(__file__).resolve().parents[2]
    )
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--gpu-mem-gb", type=float, default=16.0)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--nice", type=int, default=10)
    parser.add_argument("--log", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--executor-identity", default=EXECUTOR_IDENTITY)
    args = parser.parse_args()
    if args.threads < 1 or args.threads > 16:
        parser.error("--threads must be within 1..16")
    if args.nice < 10:
        parser.error("--nice must be >=10")
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.environ["MKL_NUM_THREADS"] = str(args.threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(args.threads)
    try:
        source = validate_source(args.source_root, repo_root=args.worktree)
        revision = resolve_code_revision(args.worktree)
        if args.dry_run:
            print(
                json.dumps(
                    dry_run(source, code_revision=revision), indent=2, sort_keys=True
                )
            )
            return 0
        if not args.gpu:
            parser.error("--gpu is required for non-dry runs")
        preload_cutensor()
        terminal = run_experiment(
            args.output_root,
            source,
            code_revision=revision,
            executor_identity=args.executor_identity,
        )
    except Exception as exc:
        args.output_root.mkdir(parents=True, exist_ok=True)
        terminal = {
            "schema": "a3g-candidate-terminal-v1",
            "written_at": now(),
            "executor_identity": args.executor_identity,
            "classification": INCONCLUSIVE_CANDIDATE,
            "terminal_stage": "preflight-or-runner-error",
            "detail": f"{type(exc).__name__}: {exc}",
            "independent_verification_required": True,
            "forbidden_outputs_emitted": False,
        }
        with exclusive_run(args.output_root):
            revoke_terminals(args.output_root)
            atomic_json(args.output_root / "candidate-terminal.json", terminal)
        print(json.dumps(terminal, indent=2, sort_keys=True))
        return 1
    print(json.dumps(terminal, indent=2, sort_keys=True))
    return 0 if terminal["classification"] != INCONCLUSIVE_CANDIDATE else 2


if __name__ == "__main__":
    raise SystemExit(main())
