#!/usr/bin/env python3
"""Independent optimizer-free verifier for the bounded A3c experiment."""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from quarry.pipeline import DftSettings, energy  # noqa: E402
from scripts import a3b_proton_microstate_stability as a3b  # noqa: E402
from scripts import a3c_triad_conditioning as a3c  # noqa: E402

DEFAULT_VERIFIER_IDENTITY = "hermes-a3c-independent-verifier"

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    bootstrap_cli(
        "a3c_verify",
        default_run_root=a3c.DEFAULT_OUTPUT_ROOT / "logs",
        gpu_owner="a3c_verify",
    )


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def _expected_signature(
    source: a3c.SourceEvidence,
    settings: DftSettings,
    code_revision: str,
) -> dict[str, Any]:
    return a3c.signature(source, settings, code_revision)


def _verify_receipt(
    output_root: Path,
    source: a3c.SourceEvidence,
    settings: DftSettings,
    code_revision: str,
    *,
    recompute_calculators: bool,
) -> dict[str, Any]:
    receipt_path = output_root / "receipt.json"
    reservation_path = output_root / "reservation.json"
    receipt = json.loads(receipt_path.read_text())
    reservation = json.loads(reservation_path.read_text())
    expected_signature = _expected_signature(source, settings, code_revision)
    if (
        reservation.get("schema") != "a3c-stage-reservation-v1"
        or reservation.get("status") != "reserved"
        or reservation.get("signature") != expected_signature
    ):
        raise RuntimeError("reservation identity mismatch")
    if (
        receipt.get("schema") != "a3c-stage-receipt-v1"
        or receipt.get("signature") != expected_signature
        or receipt.get("source") != expected_signature["source"]
        or receipt.get("settings") != expected_signature["settings"]
        or receipt.get("seed") != expected_signature["seed"]
        or receipt.get("constraints") != expected_signature["constraints"]
        or receipt.get("budget") != expected_signature["budget"]
        or receipt.get("code_revision") != code_revision
    ):
        raise RuntimeError("receipt identity/source/settings/seed/budget mismatch")
    optimizer = receipt.get("optimizer", {})
    if not (
        optimizer.get("fresh_instance") is True
        and optimizer.get("geometric_default_fresh_hessian") is True
    ):
        raise RuntimeError("optimizer freshness evidence mismatch")

    raw = receipt.get("raw_endpoint")
    if isinstance(raw, dict):
        raw_path = output_root / "raw-endpoint.xyz"
        if not raw_path.is_file() or raw.get("sha256") != a3b.sha256_path(raw_path):
            raise RuntimeError("raw endpoint missing or hash-mismatched")
    elif receipt.get("status") in {
        "complete",
        "gate-rejected",
        "pending-structural-gates",
    }:
        raise RuntimeError("post-optimizer receipt omitted raw endpoint evidence")

    verified = dict(receipt)
    if receipt.get("status") == "complete":
        endpoint_path = output_root / "endpoint.xyz"
        if not endpoint_path.is_file() or receipt.get("endpoint", {}).get(
            "sha256"
        ) != a3b.sha256_path(endpoint_path):
            raise RuntimeError("endpoint missing or hash-mismatched")
        endpoint = a3b.read_cluster(endpoint_path, source.a3a.cluster)
        expected_fingerprint = a3b.frequency_geometry_fingerprint(endpoint)
        if (
            receipt.get("endpoint", {}).get("geometry_fingerprint")
            != expected_fingerprint
        ):
            raise RuntimeError("endpoint geometry fingerprint mismatch")
        structure = a3b.structural_gate(
            endpoint,
            source.a3a.cluster,
            a3c.constraints(source),
            "triad-conditioned",
        )[1]
        if receipt.get("observed_owners") != a3b._owner_labels(endpoint):
            raise RuntimeError("endpoint owner map receipt mismatch")
        if receipt.get("structure", {}).get("constraint_residuals") != structure.get(
            "constraint_residuals"
        ):
            raise RuntimeError("constraint residual receipt mismatch")
        verified["structure"] = structure
        verified["owner_retaining"] = (
            a3c._triad_retained(endpoint) and structure["owner_changes"] == []
        )
        if recompute_calculators:
            observed_energy = float(energy(endpoint, settings))
            if not math.isfinite(observed_energy):
                raise RuntimeError("recomputed energy is non-finite")
            metrics = a3b.gradient_metrics(endpoint, settings, a3c.constraints(source))
            verified["evidence"] = {"energy_hartree": observed_energy, **metrics}
            verified["stationary"] = bool(optimizer.get("converged")) and bool(
                metrics["passed"]
            )
    return verified


def verify_experiment(
    output_root: Path,
    *,
    a3a_root: Path = a3c.DEFAULT_A3A_ROOT,
    a3b_root: Path = a3c.DEFAULT_A3B_ROOT,
    verifier_identity: str = DEFAULT_VERIFIER_IDENTITY,
    source_override: a3c.SourceEvidence | None = None,
    recompute_calculators: bool | None = None,
) -> dict[str, Any]:
    """Re-derive A3c from static evidence and independently recompute final gates."""
    output_root = output_root.resolve()
    candidate: dict[str, Any] = {}
    try:
        candidate_path = output_root / "candidate-terminal.json"
        candidate = json.loads(candidate_path.read_text())
        if candidate.get("schema") != "a3c-candidate-terminal-v1":
            raise RuntimeError("candidate terminal schema mismatch")
        executor_identity = candidate.get("executor_identity")
        if (
            not isinstance(executor_identity, str)
            or executor_identity == verifier_identity
        ):
            raise RuntimeError("verifier identity must differ from executor identity")
        source = source_override or a3c.validate_sources(a3a_root, a3b_root)
        settings_data = candidate.get("settings")
        if not isinstance(settings_data, dict):
            raise RuntimeError("candidate settings missing")
        required_settings = asdict(
            DftSettings(
                xc="b3lyp",
                basis="def2-svp",
                density_fit=True,
                use_gpu=True,
            )
        )
        if settings_data != required_settings:
            raise RuntimeError(
                "candidate is not the required B3LYP/def2-SVP/DF GPU run"
            )
        settings = DftSettings(**settings_data)
        expected_settings = asdict(settings)
        code_revision = candidate.get("code_revision")
        if not isinstance(code_revision, str) or len(code_revision) != 40:
            raise RuntimeError("candidate code revision is invalid")
        if candidate.get("source") != a3c._source_map(source):
            raise RuntimeError("candidate source binding mismatch")
        if candidate.get("settings") != expected_settings:
            raise RuntimeError("candidate settings binding mismatch")
        if candidate.get("receipt_sha256") != a3b.sha256_path(
            output_root / "receipt.json"
        ):
            raise RuntimeError("candidate receipt hash mismatch")
        should_recompute = (
            source_override is None
            if recompute_calculators is None
            else recompute_calculators
        )
        receipt = _verify_receipt(
            output_root,
            source,
            settings,
            code_revision,
            recompute_calculators=should_recompute,
        )
        classification = a3c._classify_receipt(receipt)
        if classification != candidate.get("classification"):
            raise RuntimeError(
                "candidate classification mismatch: "
                f"claimed {candidate.get('classification')!r}, "
                f"derived {classification!r}"
            )
        verified_classification = (
            a3c.VERIFIED_SEED
            if classification == a3c.CANDIDATE_SEED
            else a3c.VERIFIED_FAILURE
        )
        result = {
            "schema": "a3c-verified-terminal-v1",
            "written_at": now(),
            "status": "verified",
            "classification": verified_classification,
            "executor_identity": executor_identity,
            "verifier_identity": verifier_identity,
            "candidate_terminal_sha256": a3b.sha256_path(candidate_path),
            "receipt_sha256": a3b.sha256_path(output_root / "receipt.json"),
            "source_rehashed": source_override is None,
            "calculator_evidence_recomputed": should_recompute,
            "code_revision": code_revision,
            "terminal_stage": receipt.get("terminal_stage"),
        }
    except Exception as exc:
        result = {
            "schema": "a3c-verified-terminal-v1",
            "written_at": now(),
            "status": "rejected",
            "classification": a3c.VERIFIED_FAILURE,
            "executor_identity": candidate.get("executor_identity"),
            "verifier_identity": verifier_identity,
            "detail": f"{type(exc).__name__}: {exc}",
        }
    a3b.atomic_json(output_root / "verified-terminal.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a3a-root", type=Path, default=a3c.DEFAULT_A3A_ROOT)
    parser.add_argument("--a3b-root", type=Path, default=a3c.DEFAULT_A3B_ROOT)
    parser.add_argument("--output-root", type=Path, default=a3c.DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--gpu-mem-gb", type=float, default=16.0)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--nice", type=int, default=10)
    parser.add_argument("--log", type=Path)
    parser.add_argument("--verifier-identity", default=DEFAULT_VERIFIER_IDENTITY)
    args = parser.parse_args()
    if args.threads > 16:
        parser.error("--threads must be <=16")
    os.environ["OMP_NUM_THREADS"] = str(args.threads)
    os.environ["MKL_NUM_THREADS"] = str(args.threads)
    os.environ["OPENBLAS_NUM_THREADS"] = str(args.threads)
    candidate = json.loads((args.output_root / "candidate-terminal.json").read_text())
    if bool(candidate.get("settings", {}).get("use_gpu")) != args.gpu:
        parser.error("--gpu must match the candidate's recorded use_gpu setting")
    result = verify_experiment(
        args.output_root,
        a3a_root=args.a3a_root,
        a3b_root=args.a3b_root,
        verifier_identity=args.verifier_identity,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "verified" else 1


if __name__ == "__main__":
    raise SystemExit(main())
