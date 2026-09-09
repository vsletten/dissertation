"""Optimizer-free adjudication tests for terminal CALC-005 failures."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

from scripts import calc005_si_attachment as driver
from scripts import calc005_si_attachment_verify as verifier

REPO = Path(__file__).resolve().parents[2]
DECK = REPO / "petra/examples/kaolinite.toml"
EXPECTED_UNIT = {
    "MainPID": "0",
    "Result": "exit-code",
    "ExecMainStatus": "2",
    "ActiveState": "failed",
    "SubState": "failed",
}


class CallbackFailureBackend:
    production_backend = False

    def optimize(self, role, cluster, settings, *, max_steps, on_optimizer_enter):
        del role, cluster, settings, max_steps
        on_optimizer_enter()
        raise RuntimeError("Nuclear gradients of <test DFRKS_Scanner> not converged")

    def gradient(self, *args, **kwargs):
        raise AssertionError("failure adjudication must not calculate gradients")

    def frequencies(self, *args, **kwargs):
        raise AssertionError("failure adjudication must not calculate frequencies")

    def energy(self, *args, **kwargs):
        raise AssertionError("failure adjudication must not calculate energies")


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _failure_root(tmp_path: Path) -> tuple[Path, str]:
    commit = subprocess.run(
        ["git", "-C", str(REPO), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    source_provenance = {
        "repository": str(REPO),
        "branch": "agents/test-failure-verifier",
        "head": commit,
        "origin_head": commit,
        "clean": True,
        "tracked_source_files": driver._required_production_sources(DECK),
    }
    terminal = driver.run_pilot(
        tmp_path,
        DECK,
        CallbackFailureBackend(),
        executor_identity="executor-a",
        source_commit=commit,
        use_gpu=True,
        execution_envelope={"mode": "in-process-test"},
        source_provenance=source_provenance,
    )
    assert terminal["verdict"] == "incomplete-computational-failure"
    log = tmp_path / "logs/a3i-executor.log"
    log.parent.mkdir(parents=True)
    log.write_text(
        "geomeTRIC started.\n"
        "maxiter                   150\n"
        "Step    0 : Gradient = 1.0e-1/2.0e-1 (rms/max)\n"
        + json.dumps(terminal, indent=2, sort_keys=True)
        + "\n"
    )
    prior_hash = hashlib.sha256(
        (json.dumps(terminal, indent=2, sort_keys=True) + "\n").encode()
    ).hexdigest()
    terminal["restoration"] = {
        "status": "verified-restored",
        "checked_at": "2099-09-09T02:00:00-07:00",
        "bootstrap_session_closed": True,
        "qi2_lease_released": True,
        "lease_path": str(tmp_path / "lease.json"),
        "prior_terminal_sha256": prior_hash,
        "unit": "test-calc005.service",
        "unit_active_state": "failed",
        "unit_result": "exit-code",
        "unit_exec_main_status": 2,
        "unit_main_pid": 0,
    }
    _write_json(tmp_path / "terminal-receipt.json", terminal)
    return tmp_path, commit


def _current_blob(repo: Path, _commit: str, relative: str) -> bytes:
    return (repo / relative).read_bytes()


def _restored(_restoration: dict) -> dict:
    return {"lease_present": False, "unit": EXPECTED_UNIT}


def _verify(root: Path) -> dict:
    return verifier.verify_failed_pilot(
        root,
        DECK,
        verifier_identity="different-worker",
        repo=REPO,
        blob_reader=_current_blob,
        restoration_probe=_restored,
    )


def test_failure_verifier_passes_terminal_rejection_without_calculators(tmp_path):
    root, commit = _failure_root(tmp_path)
    before = {
        str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in root.rglob("*")
        if path.is_file()
    }

    result = _verify(root)

    assert result["status"] == "verified-pass"
    assert result["pilot_disposition"] == "terminally-rejected"
    assert result["source_commit"] == commit
    assert result["optimizer_calls"] == 0
    assert result["calculator_calls"] == 0
    assert result["failed_edges"] == []
    assert result["optimizer_budget"] == {
        "C": {"calls": 1, "max_steps": 150, "retries": 0},
        "V": {"calls": 0, "max_steps": 150, "retries": 0},
        "SiOH4": {"calls": 0, "max_steps": 150, "retries": 0},
    }
    assert result["physical_state"]["condensed_state"] == 204
    assert result["physical_state"]["cycle"] == ("Al6H38O30Si -> Al6H34O26 + H4O4Si")
    assert result["petra_boundary"]["production_value_emitted"] is False
    after = {
        str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in root.rglob("*")
        if path.is_file() and path.name != "verified-terminal.json"
    }
    assert after == before
    assert json.loads((root / "verified-terminal.json").read_text()) == result


def test_failure_verifier_rejects_an_endpoint_after_callback_failure(tmp_path):
    root, _commit = _failure_root(tmp_path)
    (root / "components/C/raw-endpoint.xyz").write_text("forbidden\n")

    result = _verify(root)

    assert result["status"] == "verified-fail"
    assert "checkpoint namespace" in result["detail"]
    assert result["optimizer_calls"] == 0
    assert result["calculator_calls"] == 0


def test_failure_verifier_rejects_cross_receipt_source_drift(tmp_path):
    root, _commit = _failure_root(tmp_path)
    terminal_path = root / "terminal-receipt.json"
    terminal = json.loads(terminal_path.read_text())
    terminal["source"]["deck_sha256"] = "0" * 64
    _write_json(terminal_path, terminal)

    result = _verify(root)

    assert result["status"] == "verified-fail"
    assert "source provenance copies differ" in result["detail"]


def test_failure_verifier_rejects_live_restoration_drift(tmp_path):
    root, _commit = _failure_root(tmp_path)
    result = verifier.verify_failed_pilot(
        root,
        DECK,
        verifier_identity="different-worker",
        repo=REPO,
        blob_reader=_current_blob,
        restoration_probe=lambda _restoration: {
            "lease_present": True,
            "unit": EXPECTED_UNIT,
        },
    )

    assert result["status"] == "verified-fail"
    assert "live executor restoration state differs" in result["detail"]
