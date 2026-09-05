#!/usr/bin/env python3
"""Hash-gated launcher for the single A3a Osa-neutral reactant recovery."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from collections.abc import Sequence
from datetime import datetime
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    bootstrap_cli(
        "a3a_reactant_recovery",
        default_run_root=Path(
            "/mnt/data/vsletten/dissertation-data/a3a-reactant-minimum-recovery/logs"
        ),
        argv=["--nice", "0"],
    )

from quarry.clusters import Cluster, water  # noqa: E402
from quarry.crystal import attack_complex, from_deck_cell  # noqa: E402
from quarry.store import geometry_hash  # noqa: E402
from scripts import phase2_ladder as phase2  # noqa: E402

DEFAULT_SOURCE_ROOT = Path(
    "/mnt/data/vsletten/dissertation-data/task274-a3-barrier-ladder-20260905"
)
DEFAULT_OUTPUT_ROOT = Path(
    "/mnt/data/vsletten/dissertation-data/a3a-reactant-minimum-recovery"
)
SOURCE_ATTEMPTS = (
    (
        Path("."),
        "6e2bb923b461515295eae8f0bc357a6c13e9614904a3c0bc2833fe17272c4071",
    ),
    (
        Path("recovery-hf-preopt-bc73b8e"),
        "a25ba7611cd0d6d793e1a5fe291ce319c8a8cd825a14ad7e9bff23ccc48d9a66",
    ),
    (
        Path("recovery-advisory-1338c6f"),
        "982ff963d9b2e421e0102881687771c085708e1ef2a938e7bff80bfa6b6092d9",
    ),
)
CELL_RELATIVE = Path("runs/phase2/osa-neutral-n1-s2-b3lyp-def2-svp")
EXPECTED_ADVISORY_OWNER_CHANGES = (
    "H50:O26->O31",
    "H52:O27->O32",
    "H57:O29->O20",
)


def now() -> str:
    return datetime.now().astimezone().isoformat()


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def canonical_complex_guess(repo_root: Path) -> Cluster:
    cell = from_deck_cell(
        repo_root / "petra/examples/kaolinite.toml",
        "Osa",
        metal_shells=2,
        n_intact=1,
        target_charge=0,
    )
    complex_guess, _ = attack_complex(cell, water())
    return complex_guess


def git_head(repo_root: Path) -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def outcome_artifact_hashes(output_root: Path) -> list[dict[str, object]]:
    excluded = {"source-evidence-manifest.json", "terminal-receipt.json"}
    return [
        {
            "path": str(path.relative_to(output_root)),
            "size": path.stat().st_size,
            "sha256": phase2.sha256_path(path),
        }
        for path in sorted(output_root.rglob("*"))
        if path.is_file()
        and path.name not in excluded
        and not path.name.endswith(".tmp")
    ]


def validate_source_evidence(
    source_root: Path,
    complex_guess: Cluster,
    *,
    attempts: Sequence[tuple[Path, str]] = SOURCE_ATTEMPTS,
    expected_owner_changes: Sequence[str] = EXPECTED_ADVISORY_OWNER_CHANGES,
) -> dict[str, object]:
    """Rehash every failed attempt and reject the drifted advisory seed."""
    records: list[dict[str, object]] = []
    expected_geometry_hash = geometry_hash(complex_guess.to_xyz())
    latest_attempt_root: Path | None = None
    for attempt_index, (relative_root, expected_terminal_sha) in enumerate(
        attempts, start=1
    ):
        attempt_root = (source_root / relative_root).resolve()
        terminal = attempt_root / "terminal-receipt.json"
        observed_terminal_sha = phase2.sha256_path(terminal)
        if observed_terminal_sha != expected_terminal_sha:
            raise RuntimeError(
                f"attempt {attempt_index} terminal receipt SHA-256 mismatch: "
                f"{observed_terminal_sha}"
            )
        terminal_payload = json.loads(terminal.read_text())
        if (
            terminal_payload.get("success") is not False
            or terminal_payload.get("family") != "osa"
            or terminal_payload.get("state") != "neutral"
            or 1 not in terminal_payload.get("requested_n_intact", [])
            or terminal_payload.get("current_n_intact") != 1
            or terminal_payload.get("observed_git_sha")
            != terminal_payload.get("expected_git_sha")
        ):
            raise RuntimeError(f"attempt {attempt_index} terminal identity mismatch")
        source_guess_path = attempt_root / CELL_RELATIVE / "complex_guess.xyz"
        source_guess = phase2.load_xyz(source_guess_path, complex_guess)
        source_geometry_hash = geometry_hash(source_guess.to_xyz())
        if source_geometry_hash != expected_geometry_hash:
            raise RuntimeError(
                f"attempt {attempt_index} source complex geometry drifted"
            )
        records.append(
            {
                "attempt": attempt_index,
                "root": str(attempt_root),
                "terminal_receipt": str(terminal),
                "terminal_receipt_sha256": observed_terminal_sha,
                "terminal_expected_git_sha": terminal_payload.get("expected_git_sha"),
                "terminal_reason": terminal_payload.get("terminal_reason"),
                "complex_guess": str(source_guess_path),
                "complex_guess_sha256": phase2.sha256_path(source_guess_path),
                "complex_guess_geometry_hash": source_geometry_hash,
            }
        )
        latest_attempt_root = attempt_root

    if latest_attempt_root is None:
        raise RuntimeError("no source attempts were supplied")
    advisory_path = latest_attempt_root / CELL_RELATIVE / "complex_preopt.xyz"
    advisory = phase2.load_xyz(advisory_path, complex_guess)
    observed_changes = phase2.proton_owner_changes(complex_guess, advisory)
    if observed_changes != list(expected_owner_changes):
        raise RuntimeError(
            "latest advisory seed did not reproduce the pinned owner changes: "
            f"{observed_changes}"
        )
    try:
        phase2.reactant_geometry_gate(
            advisory, complex_guess, stage="latest advisory seed"
        )
    except RuntimeError as exc:
        rejection = str(exc)
    else:
        raise RuntimeError(
            "latest advisory seed unexpectedly passed the microstate gate"
        )

    return {
        "schema": "a3a-source-evidence-manifest-v1",
        "validated_at": now(),
        "canonical_complex_guess_geometry_hash": expected_geometry_hash,
        "attempts": records,
        "latest_advisory_seed": {
            "path": str(advisory_path),
            "sha256": phase2.sha256_path(advisory_path),
            "owner_changes": observed_changes,
            "rejection": rejection,
            "production_calculator_called": False,
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument(
        "--worktree", type=Path, default=Path(__file__).resolve().parents[2]
    )
    args = parser.parse_args()

    started_at = now()
    started_monotonic = time.monotonic()
    args.output_root.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_root / "source-evidence-manifest.json"
    terminal_path = args.output_root / "terminal-receipt.json"
    terminal_path.unlink(missing_ok=True)
    driver_run_dir = args.output_root / CELL_RELATIVE
    status = "failed"
    stage = "source-evidence-validation"
    detail = "runner did not reach the driver"
    driver_returncode: int | None = None
    manifest: dict[str, object] | None = None

    try:
        complex_guess = canonical_complex_guess(args.worktree)
        manifest = validate_source_evidence(args.source_root, complex_guess)
        implementation_sha = git_head(args.worktree)
        manifest["implementation"] = {
            "git_sha": implementation_sha,
            "driver": str(args.worktree / "qm/scripts/phase2_ladder.py"),
            "driver_sha256": phase2.sha256_path(
                args.worktree / "qm/scripts/phase2_ladder.py"
            ),
            "runner": str(Path(__file__).resolve()),
            "runner_sha256": phase2.sha256_path(Path(__file__).resolve()),
        }
        atomic_json(manifest_path, manifest)

        stage = "reactant-recovery-driver"
        log_path = args.output_root / "logs/a3a-reactant-recovery.log"
        log_path.parent.mkdir(parents=True, exist_ok=True)
        command = [
            sys.executable,
            str(args.worktree / "qm/scripts/phase2_ladder.py"),
            "--family",
            "osa",
            "--state",
            "neutral",
            "--n-intact",
            "1",
            "--xc",
            "b3lyp",
            "--basis",
            "def2-svp",
            "--gpu",
            "--gpu-mem-gb",
            "16",
            "--threads",
            "16",
            "--nice",
            "10",
            "--reactant-recovery-only",
            "--run-root",
            str(args.output_root / "runs"),
            "--log",
            str(log_path),
        ]
        child_env = os.environ.copy()
        worktree_pythonpath = str(args.worktree / "qm")
        inherited_pythonpath = child_env.get("PYTHONPATH")
        child_env["PYTHONPATH"] = (
            f"{worktree_pythonpath}{os.pathsep}{inherited_pythonpath}"
            if inherited_pythonpath
            else worktree_pythonpath
        )
        completed = subprocess.run(
            command, cwd=args.worktree, env=child_env, check=False
        )
        driver_returncode = completed.returncode
        if completed.returncode != 0:
            production_terminal = driver_run_dir / "production-terminal.json"
            production_detail = (
                json.loads(production_terminal.read_text())
                if production_terminal.exists()
                else None
            )
            detail = f"driver exited {completed.returncode}: {production_detail}"
            raise RuntimeError(detail)

        stage = "result-verification"
        required = [
            driver_run_dir / "complex.xyz",
            driver_run_dir / "complex.json",
            driver_run_dir / "production-terminal.json",
            driver_run_dir / "reactant-result.json",
        ]
        missing = [str(path) for path in required if not path.is_file()]
        if missing:
            raise RuntimeError(f"driver success omitted required artifacts: {missing}")
        forbidden = [
            driver_run_dir / "results.json",
            driver_run_dir / "store.sqlite",
            driver_run_dir / "ts.xyz",
            driver_run_dir / "barrier.json",
        ]
        present = [str(path) for path in forbidden if path.exists()]
        if present:
            raise RuntimeError(
                f"reactant-only run emitted forbidden outputs: {present}"
            )
        result = json.loads((driver_run_dir / "reactant-result.json").read_text())
        if result.get("stopped_before_saddle_or_barrier") is not True:
            raise RuntimeError("reactant-only stop receipt is missing")
        status = "success"
        detail = "Osa-neutral n=1 reactant minimum accepted; no saddle/barrier work ran"
    except Exception as exc:
        if not detail.startswith("driver exited"):
            detail = f"{type(exc).__name__}: {exc}"
    finally:
        if manifest is not None:
            try:
                manifest["outcome"] = {
                    "status": status,
                    "stage": stage,
                    "artifacts": outcome_artifact_hashes(args.output_root),
                }
                atomic_json(manifest_path, manifest)
            except Exception as exc:
                status = "failed"
                stage = "hash-manifest-finalization"
                detail = f"{type(exc).__name__}: {exc}"
        payload: dict[str, object] = {
            "schema": "a3a-reactant-minimum-terminal-v1",
            "started_at": started_at,
            "completed_at": now(),
            "elapsed_seconds": time.monotonic() - started_monotonic,
            "status": status,
            "stage": stage,
            "detail": detail,
            "driver_returncode": driver_returncode,
            "source_manifest": str(manifest_path),
            "source_manifest_sha256": (
                phase2.sha256_path(manifest_path) if manifest_path.exists() else None
            ),
            "driver_run_dir": str(driver_run_dir),
        }
        atomic_json(terminal_path, payload)

    return 0 if status == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
