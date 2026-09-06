"""Unit gates for the bounded A3 family campaign supervisor."""

from __future__ import annotations

import hashlib
import json
import signal
import sqlite3
import subprocess

import pytest

from scripts import a3_family_campaign as campaign


def write_valid_result(tmp_path, **changes):
    cell = "oss-neutral-n2-s2"
    method = "b3lyp/def2-svp/df"
    xyz = {
        "complex": "1\ncomplex\nH 0.0 0.0 0.0\n",
        "ts": "1\nts\nH 0.1 0.0 0.0\n",
    }
    payload = {
        "cell": cell,
        "family": "oss",
        "state": "neutral",
        "n_intact": 2,
        "method": method,
        "geometry_hash": {
            role: hashlib.sha256(text.encode("utf-8")).hexdigest()
            for role, text in xyz.items()
        },
        "dG_kj": 123.4,
        "dH_kj": 111.0,
        "ts_imaginary_cm": 98.0,
        "route": "proton-neb",
    }
    payload.update(changes)
    result = tmp_path / "results.json"
    result.write_text(json.dumps(payload))
    store = tmp_path / "store.sqlite"
    with sqlite3.connect(store) as connection:
        connection.executescript(
            "CREATE TABLE structures "
            "(id INTEGER PRIMARY KEY, name TEXT, xyz TEXT, geometry_hash TEXT);"
            "CREATE TABLE jobs "
            "(id INTEGER PRIMARY KEY, structure_id INTEGER, kind TEXT, "
            "method TEXT, engine TEXT, status TEXT);"
            "CREATE TABLE results (job_id INTEGER, key TEXT, value REAL, units TEXT);"
        )
        for index, role in enumerate(("complex", "ts"), start=1):
            connection.execute(
                "INSERT INTO structures VALUES (?, ?, ?, ?)",
                (
                    index,
                    f"{payload['cell']}-{role}",
                    xyz[role],
                    payload["geometry_hash"][role],
                ),
            )
            connection.execute(
                "INSERT INTO jobs VALUES (?, ?, 'freq', ?, 'gpu4pyscf', 'done')",
                (index, index, payload["method"]),
            )
            connection.execute(
                "INSERT INTO results VALUES (?, 'electronic', ?, 'hartree')",
                (index, -100.0 + index),
            )
    return result, store


def test_wait_for_gpu_retries_until_available():
    responses = iter([(False, "busy pid=7"), (False, "still busy"), (True, "free")])
    clock = iter([0.0, 0.0, 1.0, 2.0])
    sleeps: list[float] = []
    heartbeats: list[str] = []

    message = campaign.wait_for_gpu(
        timeout_seconds=10,
        poll_seconds=1,
        heartbeat=heartbeats.append,
        lane_probe=lambda: next(responses),
        monotonic=lambda: next(clock),
        sleep=sleeps.append,
    )

    assert message == "free"
    assert heartbeats == ["busy pid=7", "still busy", "free"]
    assert sleeps == [1, 1]


def test_wait_for_gpu_times_out_with_last_owner():
    clock = iter([0.0, 0.0, 2.0])
    with pytest.raises(TimeoutError, match="live owner"):
        campaign.wait_for_gpu(
            timeout_seconds=1,
            poll_seconds=1,
            heartbeat=lambda _message: None,
            lane_probe=lambda: (False, "live owner"),
            monotonic=lambda: next(clock),
            sleep=lambda _seconds: None,
        )


def test_validate_result_binds_identity_finite_values_and_hashes(tmp_path):
    result, store = write_valid_result(tmp_path)

    record = campaign.validate_result(result, family="oss", state="neutral", n_intact=2)

    assert record["dG_kj"] == 123.4
    assert record["route"] == "proton-neb"
    assert record["result_sha256"] == campaign.sha256_path(result)
    assert record["store_sha256"] == campaign.sha256_path(store)


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"n_intact": 3}, "identity mismatch"),
        ({"dG_kj": float("nan")}, "invalid dG_kj"),
        ({"route": ""}, "no route provenance"),
    ],
)
def test_validate_result_rejects_false_green_payloads(tmp_path, change, message):
    changes = {"route": "direct", **change}
    result, _store = write_valid_result(tmp_path, **changes)

    with pytest.raises(RuntimeError, match=message):
        campaign.validate_result(result, family="oss", state="neutral", n_intact=2)


def test_validate_result_rejects_unrelated_sqlite_store(tmp_path):
    result, store = write_valid_result(tmp_path)
    store.unlink()
    with sqlite3.connect(store) as connection:
        connection.execute("CREATE TABLE unrelated (value TEXT)")

    with pytest.raises(RuntimeError, match="store schema is incomplete"):
        campaign.validate_result(result, family="oss", state="neutral", n_intact=2)


@pytest.mark.parametrize(
    "change",
    [
        {"cell": "not-oss-neutral-n2-s2"},
        {"method": "hf/sto-3g/df"},
    ],
)
def test_validate_result_rejects_self_consistent_wrong_identity(tmp_path, change):
    result, _store = write_valid_result(tmp_path, **change)

    with pytest.raises(RuntimeError, match="wrong (cell identity|method)"):
        campaign.validate_result(result, family="oss", state="neutral", n_intact=2)


def test_validate_result_rejects_non_gpu_engine(tmp_path):
    result, store = write_valid_result(tmp_path)
    with sqlite3.connect(store) as connection:
        connection.execute("UPDATE jobs SET engine = 'pyscf'")

    with pytest.raises(RuntimeError, match="store provenance mismatch"):
        campaign.validate_result(result, family="oss", state="neutral", n_intact=2)


def test_validate_result_rejects_wrong_oaa_center(tmp_path):
    result, _store = write_valid_result(
        tmp_path,
        family="oaa",
        state="neutral",
        n_intact=4,
        cell="oaa-neutral-n4-s2",
        cluster={"center_site": 18},
    )

    with pytest.raises(RuntimeError, match="wrong crystallographic center"):
        campaign.validate_result(result, family="oaa", state="neutral", n_intact=4)


def test_oaa_campaign_rejects_unconstructable_odd_connectivity(tmp_path):
    with pytest.raises(SystemExit, match="drawn from 2,4,6"):
        campaign.main(
            [
                "--worktree",
                str(tmp_path),
                "--run-root",
                str(tmp_path / "run"),
                "--expected-git-sha",
                "deadbeef",
                "--family",
                "oaa",
                "--state",
                "neutral",
                "--cells",
                "1",
                "2",
                "4",
                "6",
            ]
        )


def test_oaa_n4_campaign_pins_exact_crystallographic_center(tmp_path, monkeypatch):
    worktree = tmp_path / "worktree"
    worktree.mkdir()
    run_root = tmp_path / "run"
    run_root.mkdir()
    commands = []
    monkeypatch.setattr(campaign, "verify_revision", lambda *_args: ("branch", "sha"))
    monkeypatch.setattr(campaign, "observed_git_sha", lambda *_args: "sha")
    monkeypatch.setattr(campaign, "wait_for_gpu", lambda **_kwargs: "free")
    monkeypatch.setattr(
        campaign,
        "validate_result",
        lambda *_args, **_kwargs: {"n_intact": 4},
    )

    def run(command, **_kwargs):
        commands.append(command)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(campaign.subprocess, "run", run)
    assert (
        campaign.main(
            [
                "--worktree",
                str(worktree),
                "--run-root",
                str(run_root),
                "--expected-git-sha",
                "sha",
                "--family",
                "oaa",
                "--state",
                "neutral",
                "--cells",
                "4",
            ]
        )
        == 0
    )
    center_flag = commands[0].index("--center-index")
    assert commands[0][center_flag + 1] == "23"


def test_atomic_json_never_emits_nan(tmp_path):
    path = tmp_path / "receipt.json"
    with pytest.raises(ValueError):
        campaign.atomic_json(path, {"bad": float("nan")})
    assert not path.exists()


def test_main_writes_receipt_when_observed_git_sha_lookup_fails(tmp_path):
    missing = tmp_path / "missing-worktree"
    run_root = tmp_path / "run"
    run_root.mkdir()

    exit_code = campaign.main(
        [
            "--worktree",
            str(missing),
            "--run-root",
            str(run_root),
            "--expected-git-sha",
            "deadbeef",
            "--family",
            "oss",
            "--state",
            "neutral",
            "--cells",
            "1",
        ]
    )

    receipt_path = run_root / "terminal-receipt.json"
    assert receipt_path.is_file()
    receipt = json.loads(receipt_path.read_text())
    assert exit_code == 1
    assert receipt["success"] is False
    assert receipt["expected_git_sha"] == "deadbeef"
    assert receipt["error"]
    assert receipt["observed_git_sha"] is None


def test_main_terminal_receipt_records_final_progress_failure(tmp_path, monkeypatch):
    worktree = tmp_path / "worktree"
    worktree.mkdir()
    run_root = tmp_path / "run"
    run_root.mkdir()
    monkeypatch.setattr(campaign, "verify_revision", lambda *_args: ("branch", "sha"))
    monkeypatch.setattr(campaign, "observed_git_sha", lambda *_args: "sha")
    monkeypatch.setattr(campaign, "wait_for_gpu", lambda **_kwargs: "free")
    monkeypatch.setattr(
        campaign,
        "validate_result",
        lambda *_args, **_kwargs: {"n_intact": 1},
    )
    monkeypatch.setattr(
        campaign.subprocess,
        "run",
        lambda *_args, **_kwargs: subprocess.CompletedProcess([], 0),
    )
    real_atomic_json = campaign.atomic_json

    def fail_only_final_progress(path, payload):
        if path.name == "family-progress.json" and payload.get("status") == "complete":
            raise OSError("simulated final progress failure")
        real_atomic_json(path, payload)

    monkeypatch.setattr(campaign, "atomic_json", fail_only_final_progress)
    exit_code = campaign.main(
        [
            "--worktree",
            str(worktree),
            "--run-root",
            str(run_root),
            "--expected-git-sha",
            "sha",
            "--family",
            "oss",
            "--state",
            "neutral",
            "--cells",
            "1",
        ]
    )

    receipt = json.loads((run_root / "terminal-receipt.json").read_text())
    assert exit_code == 1
    assert receipt["success"] is False
    assert receipt["terminal_reason"] == "progress-finalization-failed"
    assert "simulated final progress failure" in receipt["error"]


def test_main_defers_signal_during_terminal_receipt_commit(tmp_path, monkeypatch):
    worktree = tmp_path / "worktree"
    worktree.mkdir()
    run_root = tmp_path / "run"
    run_root.mkdir()
    monkeypatch.setattr(campaign, "verify_revision", lambda *_args: ("branch", "sha"))
    monkeypatch.setattr(campaign, "observed_git_sha", lambda *_args: "sha")
    monkeypatch.setattr(campaign, "wait_for_gpu", lambda **_kwargs: "free")
    monkeypatch.setattr(
        campaign,
        "validate_result",
        lambda *_args, **_kwargs: {"n_intact": 1},
    )
    monkeypatch.setattr(
        campaign.subprocess,
        "run",
        lambda *_args, **_kwargs: subprocess.CompletedProcess([], 0),
    )
    real_atomic_json = campaign.atomic_json
    injected = False

    def inject_signal_during_terminal(path, payload):
        nonlocal injected
        if path.name == "terminal-receipt.json" and not injected:
            injected = True
            handler = signal.getsignal(signal.SIGTERM)
            assert callable(handler)
            handler(signal.SIGTERM, None)
        real_atomic_json(path, payload)

    monkeypatch.setattr(campaign, "atomic_json", inject_signal_during_terminal)
    exit_code = campaign.main(
        [
            "--worktree",
            str(worktree),
            "--run-root",
            str(run_root),
            "--expected-git-sha",
            "sha",
            "--family",
            "oss",
            "--state",
            "neutral",
            "--cells",
            "1",
        ]
    )

    receipt = json.loads((run_root / "terminal-receipt.json").read_text())
    assert injected is True
    assert exit_code == 0
    assert receipt["success"] is True
