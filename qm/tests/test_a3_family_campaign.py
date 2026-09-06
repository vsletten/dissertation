"""Unit gates for the bounded A3 family campaign supervisor."""

from __future__ import annotations

import hashlib
import json
import signal
import sqlite3
import subprocess
from pathlib import Path

import pytest

from quarry.crystal import from_deck_cell
from quarry.pipeline import HARTREE_TO_KJ
from quarry.store import Store, geometry_hash
from scripts import a3_family_campaign as campaign

REPO = Path(__file__).resolve().parents[2]
DECK = REPO / "petra" / "examples" / "kaolinite.toml"


def write_valid_result(tmp_path, **changes):
    complex_energy = -100.0
    ts_energy = -99.95
    xyz = {
        "complex": (
            "4\ncomplex\nSi  0.0 0.0 0.0\nO 1.6 0.0 0.0\n"
            "H   2.5  0.0 0.0\nH 1.6 0.9 0.0\n"
        ),
        "ts": (
            "4\nts\nSi 0.1  0.0 0.0\nO 1.6 0.0 0.0\nH  2.5 0.0  0.0\nH 1.6 0.9 0.0\n"
        ),
    }
    payload = {
        "cell": "oss-neutral-n2-s2",
        "family": "oss",
        "state": "neutral",
        "n_intact": 2,
        "method": "b3lyp/def2-svp/df",
        "geometry_hash": {role: geometry_hash(text) for role, text in xyz.items()},
        "cluster": {
            "site_kind": "Oss",
            "center_site": 10,
            "metal_shells": 2,
            "n_intact_requested": 2,
            "n_intact": 2,
            "n_atoms": 1,
            "n_frozen": 1,
            "formula": "Si",
            "charge": 0,
            "state": "neutral",
            "charge_offset": 0,
            "method": "b3lyp/def2-svp/df",
            "gpu": True,
        },
        "dE_elec_vs_complex_kj": (ts_energy - complex_energy) * HARTREE_TO_KJ,
        "metal_shells": 2,
        "attacked_metal": "Si",
        "dG_kj": 123.4,
        "dH_kj": 111.0,
        "ts_imaginary_cm": 98.0,
        "route": "proton-neb",
    }
    cluster_changes = changes.pop("cluster", {})
    payload.update(changes)
    family = payload["family"]
    state = payload["state"]
    n_intact = payload["n_intact"]
    payload["cluster"].update(
        {
            "site_kind": campaign.FAMILY_SITE_KINDS.get(family, "unknown"),
            "center_site": campaign.FAMILY_CELL_CENTERS.get(family, {}).get(
                n_intact, 10
            )
            or 10,
            "n_intact": n_intact,
            "n_intact_requested": n_intact,
            "state": state,
            "charge_offset": 1 if state == "acid" else 0,
        }
    )
    payload["attacked_metal"] = "Al" if family == "oaa" else "Si"
    payload["cluster"].update(cluster_changes)
    result = tmp_path / "results.json"
    result.write_text(json.dumps(payload))
    store = tmp_path / "store.sqlite"
    charge = 1 if state == "acid" else 0
    with Store(store) as evidence:
        for role, energy in (("complex", complex_energy), ("ts", ts_energy)):
            structure_id = evidence.add_structure(
                f"{payload['cell']}-{role}",
                "H2OSi",
                xyz[role],
                charge=charge,
                spin=0,
            )
            job_id = evidence.add_job(
                structure_id,
                "freq",
                payload["method"],
                "gpu4pyscf",
            )
            evidence.set_job_status(job_id, "done")
            evidence.add_result(job_id, "electronic", energy, "hartree")
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
        ({"dG_kj": True}, "invalid dG_kj"),
        ({"route": ""}, "invalid route provenance"),
        ({"route": "invented-route"}, "invalid route provenance"),
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


def test_validate_result_rejects_reduced_lookalike_store_schema(tmp_path):
    result, store = write_valid_result(tmp_path)
    store.unlink()
    with sqlite3.connect(store) as connection:
        connection.executescript(
            "CREATE TABLE structures "
            "(id INTEGER PRIMARY KEY, name TEXT, xyz TEXT, geometry_hash TEXT);"
            "CREATE TABLE jobs "
            "(id INTEGER PRIMARY KEY, structure_id INTEGER, kind TEXT, "
            "method TEXT, engine TEXT, status TEXT);"
            "CREATE TABLE results (job_id INTEGER, key TEXT, value REAL, units TEXT);"
        )

    with pytest.raises(RuntimeError, match="store schema mismatch"):
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


@pytest.mark.parametrize(("column", "value"), [("charge", 1), ("spin", 1)])
def test_validate_result_rejects_wrong_electronic_identity(tmp_path, column, value):
    result, store = write_valid_result(tmp_path)
    with sqlite3.connect(store) as connection:
        connection.execute(f"UPDATE structures SET {column} = ?", (value,))

    with pytest.raises(RuntimeError, match="store provenance mismatch"):
        campaign.validate_result(result, family="oss", state="neutral", n_intact=2)


def test_validate_result_rejects_invalid_provenance_timestamp(tmp_path):
    result, store = write_valid_result(tmp_path)
    with sqlite3.connect(store) as connection:
        connection.execute("UPDATE jobs SET created_at = 'not-an-iso-timestamp'")

    with pytest.raises(RuntimeError, match="timestamp"):
        campaign.validate_result(result, family="oss", state="neutral", n_intact=2)


def test_validate_result_rejects_foreign_key_violation(tmp_path):
    result, store = write_valid_result(tmp_path)
    with sqlite3.connect(store) as connection:
        connection.execute("UPDATE jobs SET structure_id = 999 WHERE id = 1")

    with pytest.raises(RuntimeError, match="foreign-key check failed"):
        campaign.validate_result(result, family="oss", state="neutral", n_intact=2)


def test_validate_result_rejects_energy_inconsistent_with_summary(tmp_path):
    result, store = write_valid_result(tmp_path)
    with sqlite3.connect(store) as connection:
        connection.execute(
            "UPDATE results SET value = value + 0.01 "
            "WHERE job_id = (SELECT id FROM jobs WHERE structure_id = 2)"
        )

    with pytest.raises(RuntimeError, match="electronic barrier does not match"):
        campaign.validate_result(result, family="oss", state="neutral", n_intact=2)


def test_validate_result_uses_store_canonical_geometry_hash(tmp_path):
    result, store = write_valid_result(tmp_path)
    with sqlite3.connect(store) as connection:
        xyz, stored_hash = connection.execute(
            "SELECT xyz, geometry_hash FROM structures ORDER BY id LIMIT 1"
        ).fetchone()

    assert stored_hash == geometry_hash(xyz)
    assert stored_hash != hashlib.sha256(xyz.encode()).hexdigest()
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


@pytest.mark.parametrize(
    ("n_intact", "center_site", "n_atoms", "n_frozen"),
    [(2, 18, 55, 17), (4, 23, 68, 25), (6, 18, 81, 33)],
)
def test_oaa_center_map_builds_every_exact_even_connectivity(
    n_intact, center_site, n_atoms, n_frozen
):
    assert campaign.FAMILY_CELL_CENTERS["oaa"][n_intact] == center_site

    cluster = from_deck_cell(
        DECK,
        "Oaa",
        center_index=center_site,
        metal_shells=2,
        n_intact=n_intact,
        target_charge=0,
    )

    assert cluster.center_site == center_site
    assert cluster.n_intact == n_intact
    assert len(cluster.cluster.symbols) == n_atoms
    assert len(cluster.cluster.frozen_indices) == n_frozen


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
