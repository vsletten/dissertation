#!/usr/bin/env python3
"""Bounded A3c triad-conditioned stationary-seed experiment."""

from __future__ import annotations

import argparse
import fcntl
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

DEFAULT_A3A_ROOT = Path(
    "/mnt/data/vsletten/dissertation-data/a3a-reactant-minimum-recovery"
)
DEFAULT_A3B_ROOT = Path(
    "/mnt/data/vsletten/dissertation-data/a3b-osa-neutral-n1-proton-microstate-stability"
)
DEFAULT_OUTPUT_ROOT = Path(
    "/mnt/data/vsletten/dissertation-data/a3c-osa-neutral-n1-mobile-proton-triad-conditioning"
)

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    bootstrap_cli(
        "a3c_triad_conditioning",
        default_run_root=DEFAULT_OUTPUT_ROOT / "logs",
        gpu_owner="a3c_triad_conditioning",
    )

import numpy as np  # noqa: E402

from quarry.clusters import Cluster  # noqa: E402
from quarry.pipeline import (  # noqa: E402
    DftSettings,
    energy,
    optimize_bounded,
)
from quarry.store import geometry_hash  # noqa: E402
from scripts import a3b_proton_microstate_stability as a3b  # noqa: E402
from scripts.phase2_ladder import oxygen_proton_owners, preload_cutensor  # noqa: E402

MAX_STEPS = 100
OWNER_CONSTRAINTS = ((26, 50), (27, 52), (29, 57))
EXPECTED_TRIAD_OWNERS = {50: 26, 52: 27, 57: 29}
EXECUTOR_IDENTITY = "hermes-custom-build-001"
CANDIDATE_SEED = "unverified: triad-conditioned stationary seed"
CANDIDATE_FAILURE = "unverified: triad-conditioning failure"
VERIFIED_SEED = "triad-conditioned stationary seed"
VERIFIED_FAILURE = "triad-conditioning failure"

A3B_RELATIVE_PATHS = {
    "candidate_terminal": Path("candidate-terminal.json"),
    "verified_terminal": Path("verified-terminal.json"),
    "stage_receipt": Path("stages/common-dual/receipt.json"),
    "raw_endpoint": Path("stages/common-dual/raw-endpoint.xyz"),
    "endpoint": Path("stages/common-dual/endpoint.xyz"),
}
EXPECTED_A3B_SHA256 = {
    "candidate_terminal": (
        "a575400e36e6d6dc9addb5d176eb80520233e68e6500f7d86e7ad547c792df30"
    ),
    "verified_terminal": (
        "7734ebee4fc03546ccd7fc5470934d91474d6f7a2308567e63e5188f5649e113"
    ),
    "stage_receipt": "a0226c4b260a198934523a7e02b4ce7807b7d9b77137874517bf0c0c57d93fa2",
    "raw_endpoint": "cb61c135e5f0e4b2fcb877eee6d9020c67890937fb1eb7ac6f08cc0a93ae3105",
    "endpoint": "d0ef248967812ec7461f5581769a0eddc54f17e5d1e1ef1b8d2c0be8dd899bdb",
}


@dataclass(frozen=True)
class SourceEvidence:
    a3a: a3b.SourceEvidence
    a3b_root: Path
    a3b_hashes: dict[str, str]


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def _source_map(source: SourceEvidence) -> dict[str, Any]:
    return {
        "a3a": {
            "closeout_manifest_sha256": source.a3a.closeout_manifest_sha256,
            "conditioned_xyz_sha256": source.a3a.conditioned_xyz_sha256,
            "conditioned_receipt_sha256": source.a3a.conditioned_receipt_sha256,
            "conditioned_geometry_hash": source.a3a.conditioned_geometry_hash,
        },
        "a3b": dict(sorted(source.a3b_hashes.items())),
    }


def _triad_retained(cluster: Cluster) -> bool:
    owners = oxygen_proton_owners(cluster)
    return all(
        owners.get(hydrogen) == oxygen
        for hydrogen, oxygen in EXPECTED_TRIAD_OWNERS.items()
    )


def validate_sources(
    a3a_root: Path,
    a3b_root: Path,
    *,
    expected_a3a_closeout_sha256: str = a3b.EXPECTED_CLOSEOUT_MANIFEST_SHA256,
    expected_a3a_xyz_sha256: str = a3b.EXPECTED_CONDITIONED_XYZ_SHA256,
    expected_a3a_receipt_sha256: str = a3b.EXPECTED_CONDITIONED_RECEIPT_SHA256,
    expected_a3b_sha256: dict[str, str] | None = None,
    repo_root: Path | None = None,
) -> SourceEvidence:
    """Rehash and semantically validate the fixed A3a and A3b evidence chain."""
    a3a = a3b.validate_source(
        a3a_root,
        expected_closeout_sha256=expected_a3a_closeout_sha256,
        expected_conditioned_xyz_sha256=expected_a3a_xyz_sha256,
        expected_conditioned_receipt_sha256=expected_a3a_receipt_sha256,
        repo_root=repo_root,
    )
    if not _triad_retained(a3a.cluster):
        raise RuntimeError(
            "A3a conditioned source lacks the exact H50/H52/H57 owner triad"
        )

    a3b_root = a3b_root.resolve()
    expected = expected_a3b_sha256 or EXPECTED_A3B_SHA256
    if set(expected) != set(A3B_RELATIVE_PATHS):
        raise RuntimeError("A3b expected-hash set is incomplete")
    observed = {
        name: a3b.sha256_path(a3b_root / relative)
        for name, relative in A3B_RELATIVE_PATHS.items()
    }
    for name, digest in expected.items():
        if observed[name] != digest:
            raise RuntimeError(f"A3b {name} SHA-256 mismatch: {observed[name]}")

    candidate = json.loads(
        (a3b_root / A3B_RELATIVE_PATHS["candidate_terminal"]).read_text()
    )
    verified = json.loads(
        (a3b_root / A3B_RELATIVE_PATHS["verified_terminal"]).read_text()
    )
    receipt = json.loads((a3b_root / A3B_RELATIVE_PATHS["stage_receipt"]).read_text())
    expected_a3a = _source_map(SourceEvidence(a3a, a3b_root, observed))["a3a"]
    if (
        candidate.get("schema") != "a3b-candidate-terminal-v1"
        or candidate.get("classification") != a3b.INCONCLUSIVE
        or candidate.get("source") != expected_a3a
    ):
        raise RuntimeError(
            "A3b candidate terminal identity/classification/source mismatch"
        )
    if (
        verified.get("schema") != "a3b-verified-terminal-v1"
        or verified.get("status") != "verified"
        or verified.get("classification") != a3b.INCONCLUSIVE
        or verified.get("calculator_evidence_recomputed") is not True
        or verified.get("source_rehashed") is not True
        or verified.get("candidate_terminal_sha256") != observed["candidate_terminal"]
        or verified.get("stage_receipt_sha256", {}).get("common-dual")
        != observed["stage_receipt"]
    ):
        raise RuntimeError("A3b independent verified terminal is not authoritative")
    expected_settings = asdict(
        DftSettings(xc="b3lyp", basis="def2-svp", density_fit=True, use_gpu=True)
    )
    if (
        receipt.get("schema") != "a3b-stage-receipt-v1"
        or receipt.get("status") != "complete"
        or receipt.get("stage") != "common-dual"
        or receipt.get("source") != expected_a3a
        or receipt.get("settings") != expected_settings
        or receipt.get("budget")
        != {"max_steps": MAX_STEPS, "continuation_allowed": False}
        or receipt.get("optimizer", {}).get("fresh_instance") is not True
        or receipt.get("optimizer", {}).get("geometric_default_fresh_hessian")
        is not True
        or receipt.get("stationary") is not False
        or receipt.get("owner_changes") != ["H52:O27->O32"]
        or receipt.get("raw_endpoint", {}).get("sha256") != observed["raw_endpoint"]
        or receipt.get("endpoint", {}).get("sha256") != observed["endpoint"]
    ):
        raise RuntimeError("A3b common-stage terminal evidence mismatch")
    return SourceEvidence(a3a=a3a, a3b_root=a3b_root, a3b_hashes=observed)


def resolve_code_revision(worktree: Path) -> str:
    worktree = worktree.resolve()
    head = subprocess.run(
        ["git", "-C", str(worktree), "rev-parse", "HEAD"],
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
    if dirty:
        raise RuntimeError("worktree must be fully clean before production execution")
    branch = subprocess.run(
        ["git", "-C", str(worktree), "symbolic-ref", "--short", "HEAD"],
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
    if remote != head:
        raise RuntimeError("local and remote implementation revisions differ")
    return head


@contextmanager
def exclusive_run(output_root: Path):
    """Prevent concurrent processes from spending the single optimizer budget."""
    lock_path = output_root / "run.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("A3c experiment is already active") from exc
        try:
            yield
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def revoke_canonical_terminals(output_root: Path) -> None:
    """Move stale candidate/verifier terminals out of canonical read paths."""
    existing = [
        path
        for path in (
            output_root / "candidate-terminal.json",
            output_root / "verified-terminal.json",
        )
        if path.exists()
    ]
    if not existing:
        return
    revoked = output_root / "revoked" / f"{time.time_ns()}-{os.getpid()}"
    revoked.mkdir(parents=True, exist_ok=False)
    for path in existing:
        shutil.move(path, revoked / path.name)


def constraints(source: SourceEvidence) -> list[tuple[int, int, float]]:
    cluster = source.a3a.cluster
    return [
        (
            oxygen,
            hydrogen,
            float(np.linalg.norm(cluster.coords[oxygen] - cluster.coords[hydrogen])),
        )
        for oxygen, hydrogen in OWNER_CONSTRAINTS
    ]


def signature(
    source: SourceEvidence,
    settings: DftSettings,
    code_revision: str,
) -> dict[str, Any]:
    return {
        "schema": "a3c-stage-signature-v1",
        "stage": "triad-conditioned",
        "source": _source_map(source),
        "settings": asdict(settings),
        "seed": {
            "geometry_hash": geometry_hash(source.a3a.cluster.to_xyz()),
            "geometry_fingerprint": a3b.frequency_geometry_fingerprint(
                source.a3a.cluster
            ),
        },
        "constraints": [list(item) for item in constraints(source)],
        "budget": {
            "max_steps": MAX_STEPS,
            "continuation_allowed": False,
            "retry_allowed": False,
        },
        "fresh_optimizer": True,
        "geometric_default_fresh_hessian": True,
        "code_revision": code_revision,
    }


def _classify_receipt(receipt: dict[str, Any]) -> str:
    accepted = bool(
        receipt.get("status") == "complete"
        and receipt.get("optimizer", {}).get("converged") is True
        and receipt.get("stationary") is True
        and receipt.get("owner_retaining") is True
        and receipt.get("structure", {}).get("owner_changes") == []
        and receipt.get("evidence", {}).get("passed") is True
    )
    return CANDIDATE_SEED if accepted else CANDIDATE_FAILURE


def _resume(
    output_root: Path, expected_signature: dict[str, Any]
) -> dict[str, Any] | None:
    reservation_path = output_root / "reservation.json"
    receipt_path = output_root / "receipt.json"
    if not reservation_path.exists() and not receipt_path.exists():
        return None
    if not reservation_path.is_file() or not receipt_path.is_file():
        return {
            "schema": "a3c-stage-receipt-v1",
            "status": "stale-evidence-rejected",
            "terminal_stage": "reservation-readback",
            "detail": (
                "reservation/receipt pair is incomplete; budget will not be replayed"
            ),
            "signature": expected_signature,
        }
    try:
        reservation = json.loads(reservation_path.read_text())
        receipt = json.loads(receipt_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        return {
            "schema": "a3c-stage-receipt-v1",
            "status": "stale-evidence-rejected",
            "terminal_stage": "reservation-readback",
            "detail": f"{type(exc).__name__}: {exc}",
            "signature": expected_signature,
        }
    if (
        reservation.get("schema") != "a3c-stage-reservation-v1"
        or reservation.get("status") != "reserved"
        or reservation.get("signature") != expected_signature
        or receipt.get("signature") != expected_signature
    ):
        return {
            "schema": "a3c-stage-receipt-v1",
            "status": "stale-evidence-rejected",
            "terminal_stage": "reservation-identity",
            "detail": (
                "existing reservation or receipt does not match the exact run identity"
            ),
            "signature": expected_signature,
        }
    if receipt.get("status") == "complete":
        endpoint_path = output_root / "endpoint.xyz"
        if not endpoint_path.is_file() or receipt.get("endpoint", {}).get(
            "sha256"
        ) != a3b.sha256_path(endpoint_path):
            return {
                "schema": "a3c-stage-receipt-v1",
                "status": "stale-evidence-rejected",
                "terminal_stage": "endpoint-readback",
                "detail": "completed receipt endpoint is absent or hash-mismatched",
                "signature": expected_signature,
            }
    return receipt


def _run_experiment_locked(
    output_root: Path,
    source: SourceEvidence,
    *,
    use_gpu: bool,
    code_revision: str,
    executor_identity: str = EXECUTOR_IDENTITY,
) -> dict[str, Any]:
    """Spend exactly one fresh 100-step triad-conditioned optimizer budget."""
    output_root.mkdir(parents=True, exist_ok=True)
    settings = DftSettings(
        xc="b3lyp", basis="def2-svp", density_fit=True, use_gpu=use_gpu
    )
    run_signature = signature(source, settings, code_revision)
    receipt = _resume(output_root, run_signature)
    if receipt is None:
        revoke_canonical_terminals(output_root)
        reservation = {
            "schema": "a3c-stage-reservation-v1",
            "status": "reserved",
            "reserved_at": now(),
            "signature": run_signature,
        }
        a3b.atomic_json(output_root / "reservation.json", reservation)
        try:
            optimized = optimize_bounded(
                source.a3a.cluster,
                settings,
                max_steps=MAX_STEPS,
                fixed_distances=constraints(source),
            )
        except Exception as exc:
            receipt = {
                "schema": "a3c-stage-receipt-v1",
                "status": "optimizer-failed",
                "terminal_stage": "optimizer-call",
                "completed_at": now(),
                "signature": run_signature,
                "source": run_signature["source"],
                "settings": run_signature["settings"],
                "seed": run_signature["seed"],
                "constraints": run_signature["constraints"],
                "budget": run_signature["budget"],
                "code_revision": code_revision,
                "optimizer": {
                    "converged": False,
                    "fresh_instance": True,
                    "geometric_default_fresh_hessian": True,
                },
                "detail": f"{type(exc).__name__}: {exc}",
                "stationary": False,
                "owner_retaining": False,
            }
            a3b.atomic_json(output_root / "receipt.json", receipt)
        else:
            a3b.atomic_xyz(output_root / "raw-endpoint.xyz", optimized.cluster)
            receipt = {
                "schema": "a3c-stage-receipt-v1",
                "status": "pending-structural-gates",
                "terminal_stage": "structural-gates",
                "signature": run_signature,
                "source": run_signature["source"],
                "settings": run_signature["settings"],
                "seed": run_signature["seed"],
                "constraints": run_signature["constraints"],
                "budget": run_signature["budget"],
                "code_revision": code_revision,
                "optimizer": {
                    "converged": bool(optimized.converged),
                    "fresh_instance": True,
                    "geometric_default_fresh_hessian": True,
                },
                "raw_endpoint": {
                    "path": "raw-endpoint.xyz",
                    "sha256": a3b.sha256_path(output_root / "raw-endpoint.xyz"),
                    "geometry_hash": geometry_hash(optimized.cluster.to_xyz()),
                    "geometry_fingerprint": a3b.frequency_geometry_fingerprint(
                        optimized.cluster
                    ),
                    "charge": optimized.cluster.charge,
                    "spin": optimized.cluster.spin,
                    "frozen_indices": list(optimized.cluster.frozen_indices),
                },
            }
            a3b.atomic_json(output_root / "receipt.json", receipt)
            try:
                endpoint, raw_structure = a3b.structural_gate(
                    optimized.cluster,
                    source.a3a.cluster,
                    constraints(source),
                    "triad-conditioned",
                )
                a3b.atomic_xyz(output_root / "endpoint.xyz", endpoint)
                endpoint = a3b.read_cluster(
                    output_root / "endpoint.xyz", source.a3a.cluster
                )
                structure = a3b.structural_gate(
                    endpoint,
                    source.a3a.cluster,
                    constraints(source),
                    "triad-conditioned",
                )[1]
                structure["raw_maximum_frozen_coordinate_drift_a"] = raw_structure[
                    "raw_maximum_frozen_coordinate_drift_a"
                ]
                receipt["terminal_stage"] = "energy-evidence"
                observed_energy = float(energy(endpoint, settings))
                if not math.isfinite(observed_energy):
                    raise RuntimeError("production energy is non-finite")
                receipt["terminal_stage"] = "gradient-evidence"
                metrics = a3b.gradient_metrics(endpoint, settings, constraints(source))
                owner_retaining = (
                    _triad_retained(endpoint) and structure["owner_changes"] == []
                )
                stationary = bool(optimized.converged) and bool(metrics["passed"])
                receipt.update(
                    {
                        "status": "complete",
                        "terminal_stage": "stage-complete",
                        "completed_at": now(),
                        "structure": structure,
                        "endpoint": {
                            "path": "endpoint.xyz",
                            "sha256": a3b.sha256_path(output_root / "endpoint.xyz"),
                            "geometry_hash": geometry_hash(endpoint.to_xyz()),
                            "geometry_fingerprint": a3b.frequency_geometry_fingerprint(
                                endpoint
                            ),
                            "charge": endpoint.charge,
                            "spin": endpoint.spin,
                            "frozen_indices": list(endpoint.frozen_indices),
                        },
                        "observed_owners": a3b._owner_labels(endpoint),
                        "owner_retaining": owner_retaining,
                        "stationary": stationary,
                        "evidence": {"energy_hartree": observed_energy, **metrics},
                    }
                )
            except Exception as exc:
                receipt.update(
                    {
                        "status": "gate-rejected",
                        "completed_at": now(),
                        "detail": f"{type(exc).__name__}: {exc}",
                        "stationary": False,
                        "owner_retaining": False,
                    }
                )
            a3b.atomic_json(output_root / "receipt.json", receipt)

    if receipt.get("status") == "stale-evidence-rejected":
        revoke_canonical_terminals(output_root)

    classification = _classify_receipt(receipt)
    terminal = {
        "schema": "a3c-candidate-terminal-v1",
        "written_at": now(),
        "executor_identity": executor_identity,
        "classification": classification,
        "terminal_stage": receipt.get("terminal_stage", "unknown"),
        "receipt_status": receipt.get("status", "unknown"),
        "receipt_sha256": a3b.sha256_path(output_root / "receipt.json"),
        "source": run_signature["source"],
        "settings": run_signature["settings"],
        "seed": run_signature["seed"],
        "constraints": run_signature["constraints"],
        "budget": run_signature["budget"],
        "code_revision": code_revision,
        "independent_verification_required": True,
        "forbidden_outputs_emitted": False,
    }
    candidate_path = output_root / "candidate-terminal.json"
    if candidate_path.is_file():
        existing = json.loads(candidate_path.read_text())
        stable_keys = (
            "schema",
            "executor_identity",
            "classification",
            "terminal_stage",
            "receipt_status",
            "receipt_sha256",
            "source",
            "settings",
            "seed",
            "constraints",
            "budget",
            "code_revision",
        )
        if all(existing.get(key) == terminal.get(key) for key in stable_keys):
            return existing
        revoke_canonical_terminals(output_root)
    a3b.atomic_json(candidate_path, terminal)
    return terminal


def run_experiment(
    output_root: Path,
    source: SourceEvidence,
    *,
    use_gpu: bool,
    code_revision: str,
    executor_identity: str = EXECUTOR_IDENTITY,
) -> dict[str, Any]:
    output_root.mkdir(parents=True, exist_ok=True)
    with exclusive_run(output_root):
        return _run_experiment_locked(
            output_root,
            source,
            use_gpu=use_gpu,
            code_revision=code_revision,
            executor_identity=executor_identity,
        )


def dry_run(
    source: SourceEvidence, *, use_gpu: bool, code_revision: str
) -> dict[str, Any]:
    settings = DftSettings(
        xc="b3lyp", basis="def2-svp", density_fit=True, use_gpu=use_gpu
    )
    return {
        "schema": "a3c-dry-run-v1",
        "signature": signature(source, settings, code_revision),
        "source_owners": a3b._owner_labels(source.a3a.cluster),
        "production_artifacts_written": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a3a-root", type=Path, default=DEFAULT_A3A_ROOT)
    parser.add_argument("--a3b-root", type=Path, default=DEFAULT_A3B_ROOT)
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
    if args.threads > 16:
        parser.error("--threads must be <=16")
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.environ["MKL_NUM_THREADS"] = str(args.threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(args.threads)
    if args.gpu:
        preload_cutensor()
    source = validate_sources(args.a3a_root, args.a3b_root, repo_root=args.worktree)
    revision = resolve_code_revision(args.worktree)
    if args.dry_run:
        print(
            json.dumps(
                dry_run(source, use_gpu=args.gpu, code_revision=revision),
                indent=2,
                sort_keys=True,
            )
        )
        return 0
    if not args.gpu:
        parser.error("--gpu is required for non-dry runs")
    terminal = run_experiment(
        args.output_root,
        source,
        use_gpu=args.gpu,
        code_revision=revision,
        executor_identity=args.executor_identity,
    )
    print(json.dumps(terminal, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
