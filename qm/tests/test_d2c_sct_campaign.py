from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest

from scripts import d2c_input_bundle as bundle
from scripts import d2c_sct_campaign as campaign

FIXED_GIT_SHA = "a" * 40
FIXED_DEPENDENCIES = {
    "ase": "3.29.0",
    "geometric": "1.1",
    "numpy": "2.3.2",
    "pyscf": "2.10.0",
    "sella": "2.5.1",
}
BUNDLE_ROOT = Path(__file__).parents[1] / "data" / "D2c-instanton-tier" / "d2b-inputs"


def _preflight(run_root: Path, bundle_root: Path = BUNDLE_ROOT) -> dict:
    return campaign.create_preflight_receipt(
        bundle_root,
        run_root,
        git_sha=FIXED_GIT_SHA,
        dependency_versions=FIXED_DEPENDENCIES,
        created_utc="2026-09-07T12:00:00Z",
    )


def test_dry_run_identity_and_exact_four_route_inventory(tmp_path: Path):
    first = _preflight(tmp_path / "first")
    second = _preflight(tmp_path / "second")

    assert first["state"] == "pending"
    assert first["identity"] == second["identity"]
    assert first["campaign"]["git_sha"] == FIXED_GIT_SHA
    assert first["campaign"]["bundle_manifest_sha256"] == bundle.sha256_path(
        BUNDLE_ROOT / "manifest.json"
    )
    assert list(first["routes"]) == list(bundle.ROUTES)
    assert len(first["routes"]) == 4
    for route_receipt in first["routes"].values():
        assert route_receipt["symbols"]
        assert len(route_receipt["masses_amu"]) == len(route_receipt["symbols"])
        assert all(
            stage["state"] == "pending" for stage in route_receipt["stages"].values()
        )
    assert first["campaign_stages"]["branching_common_reference_gate"]["required"]
    assert first["campaign_stages"]["final_freeze"]["required"]
    assert first["accepted_result"] is None


def test_dry_run_refuses_bundle_drift_without_writing_receipt(tmp_path: Path):
    copied = tmp_path / "bundle"
    shutil.copytree(BUNDLE_ROOT, copied)
    target = copied / bundle.ROUTES[0] / "ts.xyz"
    target.write_text(target.read_text() + "\n")
    run_root = tmp_path / "run"

    with pytest.raises(ValueError, match="bundled input drifted"):
        _preflight(run_root, copied)

    assert not (run_root / campaign.PREFLIGHT_RECEIPT).exists()


def test_preflight_receipt_is_non_overwriting(tmp_path: Path):
    run_root = tmp_path / "run"
    first = _preflight(run_root)
    receipt_path = run_root / campaign.PREFLIGHT_RECEIPT
    original_bytes = receipt_path.read_bytes()

    with pytest.raises(FileExistsError, match="already exists"):
        _preflight(run_root)

    assert receipt_path.read_bytes() == original_bytes
    assert json.loads(original_bytes)["identity"] == first["identity"]


def test_preflight_refuses_and_preserves_stale_accepted_result(tmp_path: Path):
    run_root = tmp_path / "run"
    run_root.mkdir()
    stale = run_root / "accepted-result.json"
    stale.write_text('{"state":"accepted","identity":"stale"}\n')
    original_bytes = stale.read_bytes()

    with pytest.raises(FileExistsError, match="stale accepted/final result"):
        _preflight(run_root)

    assert stale.read_bytes() == original_bytes
    assert not (run_root / campaign.PREFLIGHT_RECEIPT).exists()


def test_preflight_refuses_and_preserves_final_freeze_receipt(tmp_path: Path):
    run_root = tmp_path / "run"
    run_root.mkdir()
    stale = run_root / "final-freeze.json"
    stale.write_bytes(b"untrusted final freeze receipt\n")
    original_bytes = stale.read_bytes()

    with pytest.raises(FileExistsError, match="stale accepted/final result"):
        _preflight(run_root)

    assert stale.read_bytes() == original_bytes
    assert not (run_root / campaign.PREFLIGHT_RECEIPT).exists()


def test_preflight_refuses_and_preserves_declared_route_stage_receipt(
    tmp_path: Path,
):
    run_root = tmp_path / "run"
    stale = run_root / bundle.ROUTES[0] / "irc-forward" / "receipt.json"
    stale.parent.mkdir(parents=True)
    stale.write_bytes(b"not JSON and must not be parsed or replaced\n")
    original_bytes = stale.read_bytes()

    with pytest.raises(FileExistsError, match="stale stage receipt"):
        _preflight(run_root)

    assert stale.read_bytes() == original_bytes
    assert not (run_root / campaign.PREFLIGHT_RECEIPT).exists()


def test_cli_refuses_execution_without_dry_run(tmp_path: Path, capsys):
    with pytest.raises(SystemExit) as excinfo:
        campaign.main(
            ["--run-root", str(tmp_path / "run")],
            git_sha=FIXED_GIT_SHA,
            dependency_versions=FIXED_DEPENDENCIES,
        )

    assert excinfo.value.code == 2
    assert "--dry-run is required" in capsys.readouterr().err
    assert not (tmp_path / "run").exists()


def test_cli_enforces_compute_etiquette_bounds(tmp_path: Path, capsys):
    with pytest.raises(SystemExit):
        campaign.main(
            [
                "--dry-run",
                "--run-root",
                str(tmp_path / "run"),
                "--threads",
                "17",
            ],
            git_sha=FIXED_GIT_SHA,
            dependency_versions=FIXED_DEPENDENCIES,
        )
    assert "--threads must be <= 16" in capsys.readouterr().err

    with pytest.raises(SystemExit):
        campaign.main(
            [
                "--dry-run",
                "--run-root",
                str(tmp_path / "run"),
                "--nice",
                "9",
            ],
            git_sha=FIXED_GIT_SHA,
            dependency_versions=FIXED_DEPENDENCIES,
        )
    assert "--nice must be >= 10" in capsys.readouterr().err


def test_git_identity_refuses_dirty_worktree(monkeypatch, tmp_path: Path):
    calls = []

    def fake_run(command, **kwargs):
        calls.append(command)
        return subprocess.CompletedProcess(
            command, 0, stdout="?? untracked.py\n", stderr=""
        )

    monkeypatch.setattr(campaign.subprocess, "run", fake_run)
    with pytest.raises(RuntimeError, match="requires a clean tracked and untracked"):
        campaign._git_sha(tmp_path)

    assert calls == [["git", "status", "--porcelain", "--untracked-files=all"]]
