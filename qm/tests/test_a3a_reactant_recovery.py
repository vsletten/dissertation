from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

from quarry.clusters import Cluster
from scripts import a3a_reactant_recovery as runner
from scripts import phase2_ladder as phase2


def _source_fixture(tmp_path: Path):
    guess = Cluster(
        "source",
        ["O", "O", "H", "H"],
        np.array(
            [
                [0.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
                [0.96, 0.0, 0.0],
                [3.96, 0.0, 0.0],
            ]
        ),
    )
    attempts = []
    for attempt_number in range(1, 4):
        attempt_root = tmp_path / f"attempt-{attempt_number}"
        run_dir = attempt_root / runner.CELL_RELATIVE
        run_dir.mkdir(parents=True)
        phase2.save_xyz(guess, run_dir / "complex_guess.xyz")
        terminal = attempt_root / "terminal-receipt.json"
        terminal.write_text(
            json.dumps(
                {
                    "schema": "task274-a3-family-terminal-v1",
                    "success": False,
                    "family": "osa",
                    "state": "neutral",
                    "requested_n_intact": [1],
                    "current_n_intact": 1,
                    "expected_git_sha": f"sha-{attempt_number}",
                    "observed_git_sha": f"sha-{attempt_number}",
                    "terminal_reason": "cell-n1-failed",
                }
            )
        )
        attempts.append(
            (
                Path(f"attempt-{attempt_number}"),
                phase2.sha256_path(terminal),
            )
        )
    advisory = Cluster("advisory", guess.symbols, guess.coords.copy())
    advisory.coords[2] = np.array([2.04, 0.0, 0.0])
    phase2.save_xyz(
        advisory,
        tmp_path / "attempt-3" / runner.CELL_RELATIVE / "complex_preopt.xyz",
    )
    return guess, attempts


def test_source_evidence_is_rehashed_and_advisory_seed_rejected(tmp_path):
    guess, attempts = _source_fixture(tmp_path)

    manifest = runner.validate_source_evidence(
        tmp_path,
        guess,
        attempts=attempts,
        expected_owner_changes=["H2:O0->O1"],
    )

    assert manifest["schema"] == "a3a-source-evidence-manifest-v1"
    manifest_attempts = manifest["attempts"]
    advisory = manifest["latest_advisory_seed"]
    assert isinstance(manifest_attempts, list)
    assert isinstance(advisory, dict)
    assert len(manifest_attempts) == 3
    assert advisory["owner_changes"] == ["H2:O0->O1"]
    assert advisory["production_calculator_called"] is False
    assert "changed the reactant proton microstate" in advisory["rejection"]


def test_source_evidence_refuses_terminal_receipt_hash_drift(tmp_path):
    guess, attempts = _source_fixture(tmp_path)
    attempts[1] = (attempts[1][0], "0" * 64)

    with pytest.raises(
        RuntimeError, match="attempt 2 terminal receipt SHA-256 mismatch"
    ):
        runner.validate_source_evidence(
            tmp_path,
            guess,
            attempts=attempts,
            expected_owner_changes=["H2:O0->O1"],
        )


def test_source_evidence_refuses_input_geometry_drift(tmp_path):
    guess, attempts = _source_fixture(tmp_path)
    path = tmp_path / "attempt-2" / runner.CELL_RELATIVE / "complex_guess.xyz"
    drifted = Cluster("drifted", guess.symbols, guess.coords.copy())
    drifted.coords[0, 1] = 0.1
    phase2.save_xyz(drifted, path)
    terminal = tmp_path / "attempt-2" / "terminal-receipt.json"
    attempts[1] = (attempts[1][0], phase2.sha256_path(terminal))

    with pytest.raises(RuntimeError, match="attempt 2 source complex geometry drifted"):
        runner.validate_source_evidence(
            tmp_path,
            guess,
            attempts=attempts,
            expected_owner_changes=["H2:O0->O1"],
        )


def test_outcome_manifest_excludes_outer_wrapper_files(tmp_path):
    (tmp_path / "stable.json").write_text("{}")
    (tmp_path / "launcher.log").write_text("still growing")
    (tmp_path / "restoration-receipt.txt").write_text("written later")

    records = runner.outcome_artifact_hashes(tmp_path)

    assert [record["path"] for record in records] == ["stable.json"]


def test_previous_manifest_and_terminal_are_revoked_together(tmp_path):
    for name in ("source-evidence-manifest.json", "terminal-receipt.json"):
        (tmp_path / name).write_text("stale")

    runner.revoke_previous_receipts(tmp_path)

    assert not (tmp_path / "source-evidence-manifest.json").exists()
    assert not (tmp_path / "terminal-receipt.json").exists()


def test_driver_failure_propagates_exact_child_stage(tmp_path):
    (tmp_path / "production-terminal.json").write_text(
        json.dumps(
            {
                "stage": "production-endpoint-geometry-gate",
                "detail": "RuntimeError: proton owner changed",
            }
        )
    )

    stage, detail = runner.driver_failure(tmp_path, 1)

    assert stage == "production-endpoint-geometry-gate"
    assert detail == "driver exited 1: RuntimeError: proton owner changed"


def test_outcome_manifest_excludes_configured_recovery_log(tmp_path):
    (tmp_path / "stable.json").write_text("{}")
    log_path = tmp_path / "logs/a3a-reactant-recovery.log"
    log_path.parent.mkdir(parents=True)
    log_path.write_text("still growing")

    records = runner.outcome_artifact_hashes(tmp_path)

    assert [record["path"] for record in records] == ["stable.json"]


def test_driver_command_uses_static_relative_output_paths():
    command = runner.driver_command()

    assert command[0] == sys.executable
    assert command[1] == str(
        Path(runner.__file__).resolve().parent / "phase2_ladder.py"
    )
    assert command[command.index("--run-root") + 1] == "runs"
    assert command[command.index("--log") + 1] == "logs/a3a-reactant-recovery.log"
    assert "--reactant-recovery-only" in command
