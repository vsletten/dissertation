from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pytest

from quarry.clusters import Cluster
from quarry.pipeline import DftSettings
from scripts import a2b1_reactant_scf_recovery as recovery
from scripts import production_energetics as a2


def cluster() -> Cluster:
    return Cluster(
        name="reactant.r2scan3c",
        symbols=["H", "H"],
        coords=np.array([[0.0, 0.0, 0.0], [0.75, 0.0, 0.0]]),
        charge=-1,
        spin=0,
    )


def settings() -> DftSettings:
    return DftSettings(
        xc="wb97m-v",
        basis="def2-tzvpd",
        solvent="smd",
        density_fit=True,
    )


def evidence() -> recovery.InputEvidence:
    item = cluster()
    method = settings()
    evidence_path = Path(recovery.__file__).resolve()
    evidence_reference = {
        "path": str(evidence_path),
        "sha256": a2.sha256_path(evidence_path),
    }
    return recovery.InputEvidence(
        cluster=item,
        settings=method,
        bindings={
            "recovery_contract": recovery.RECOVERY_CONTRACT,
            "legacy_terminal": evidence_reference,
            "legacy_source_receipt": evidence_reference,
            "legacy_settings": evidence_reference,
            "reactant_geometry": {
                **evidence_reference,
                "fingerprint": a2.frequency_geometry_fingerprint(item),
                "charge": -1,
                "spin": 0,
            },
            "reactant_frequency": {
                **evidence_reference,
                "imaginary_count": 0,
            },
            "production_method": a2.PRODUCTION_METHOD,
            "production_settings_fingerprint": a2.frequency_settings_fingerprint(
                method
            ),
            "driver": str(evidence_path),
            "driver_sha256": a2.sha256_path(evidence_path),
        },
    )


def args(run_dir: Path, target: Path) -> argparse.Namespace:
    return argparse.Namespace(
        run_dir=run_dir,
        target_energy=target,
        legacy_run_root=run_dir / "legacy",
        gpu=False,
        gpu_mem_gb=18.0,
        threads=2,
        nice=0,
        wall_seconds=60,
        finalize_if_running=False,
        log=None,
    )


class FakeScf:
    def __init__(self, *, converged: bool, energy: float, density: float) -> None:
        self.converged = converged
        self.energy = energy
        self.density = np.full((2, 2), density)
        self.cycles = 0
        self.callback = None
        self.max_cycle = 0
        self.diis = True
        self.diis_start_cycle = 1
        self.damp = 0.0
        self.level_shift = 0.0
        self.received_dm0 = None

    def get_init_guess(self, *, key: str):
        assert key == "huckel"
        return np.eye(2)

    def kernel(self, *, dm0=None):
        self.received_dm0 = np.asarray(dm0)
        self.cycles = self.max_cycle
        if self.callback is not None:
            self.callback(
                {
                    "cycle": self.max_cycle - 1,
                    "e_tot": self.energy,
                    "norm_gorb": 1.5e-4,
                    "norm_ddm": 2.5e-4,
                }
            )
        return self.energy

    def make_rdm1(self):
        return self.density


def install_fakes(
    monkeypatch: pytest.MonkeyPatch, scfs: list[FakeScf]
) -> list[FakeScf]:
    queue = iter(scfs)
    monkeypatch.setattr(recovery.pipeline, "build_mol", lambda *_args: object())
    monkeypatch.setattr(recovery.pipeline, "_make_scf", lambda *_args: next(queue))
    monkeypatch.setattr(
        recovery, "validate_inputs", lambda *_args, **_kwargs: evidence()
    )
    monkeypatch.setattr(recovery, "_install_alarm", lambda _seconds: None)
    return scfs


def test_both_attempts_fail_atomically_and_publish_nothing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    target = tmp_path / "reactant.wb97m-v.energy.json"
    target.write_text("stale energy\n")
    (run_dir / "results.json").write_text("stale results\n")
    (run_dir / "store.sqlite").write_text("stale store\n")
    install_fakes(
        monkeypatch,
        [
            FakeScf(converged=False, energy=-10.0, density=1.0),
            FakeScf(converged=False, energy=-10.1, density=2.0),
        ],
    )

    with pytest.raises(RuntimeError, match="finalization did not converge"):
        recovery.execute(args(run_dir, target))

    assert not target.exists()
    assert not (run_dir / "results.json").exists()
    assert not (run_dir / "store.sqlite").exists()
    attempts = [
        json.loads((run_dir / "attempt-001-huckel-roothaan.json").read_text()),
        json.loads((run_dir / "attempt-002-cdiis-finalization.json").read_text()),
    ]
    assert [item["converged"] for item in attempts] == [False, False]
    assert attempts[0]["solver"]["algorithm"] == "damped-level-shifted-roothaan"
    assert attempts[0]["solver"]["diis"] is False
    assert attempts[1]["solver"]["algorithm"] == "unshifted-cdiis-finalization"
    assert (
        attempts[1]["initial_density"]["sha256"]
        == attempts[0]["final_density"]["sha256"]
    )
    terminal = json.loads((run_dir / "terminal-receipt.json").read_text())
    status = json.loads((run_dir / "run_status.json").read_text())
    assert terminal["outcome"] == "incomplete-computational-failure"
    assert terminal["running_record_present"] is False
    assert len(terminal["attempts"]) == 2
    assert status["status"] == "failed"
    assert len(list((run_dir / "quarantine").glob("*/*"))) == 3


def test_success_publishes_only_hash_bound_reactant_energy(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run_dir = tmp_path / "run"
    target = tmp_path / "reactant.wb97m-v.energy.json"
    scfs = install_fakes(
        monkeypatch,
        [
            FakeScf(converged=False, energy=-10.0, density=1.0),
            FakeScf(converged=True, energy=-10.25, density=2.0),
        ],
    )

    assert recovery.execute(args(run_dir, target)) == 0
    payload = json.loads(target.read_text())
    assert payload["converged"] is True
    assert payload["scf_contract"] == recovery.RECOVERY_CONTRACT
    assert payload["electronic_hartree"] == pytest.approx(-10.25)
    assert scfs[1].received_dm0 is not None
    assert np.array_equal(scfs[1].received_dm0, scfs[0].density)
    assert scfs[0].diis is False
    assert scfs[0].diis_start_cycle == 51
    assert scfs[0].damp == pytest.approx(0.5)
    assert scfs[0].level_shift == pytest.approx(0.5)
    assert scfs[1].diis is True
    assert scfs[1].diis_start_cycle == 1
    assert scfs[1].damp == pytest.approx(0.0)
    assert scfs[1].level_shift == pytest.approx(0.0)
    assert not (run_dir / "results.json").exists()
    assert not (run_dir / "store.sqlite").exists()
    terminal = json.loads((run_dir / "terminal-receipt.json").read_text())
    assert terminal["outcome"] == "success"
    assert terminal["canonical_energy"]["sha256"] == a2.sha256_path(target)

    item = evidence()
    assert a2.checkpoint_converged_energy(
        target, item.cluster, item.settings, a2.PRODUCTION_METHOD
    ) == pytest.approx(-10.25)


def test_recovery_cache_rejects_attempt_receipt_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run_dir = tmp_path / "run"
    target = tmp_path / "reactant.wb97m-v.energy.json"
    install_fakes(
        monkeypatch,
        [
            FakeScf(converged=False, energy=-10.0, density=1.0),
            FakeScf(converged=True, energy=-10.25, density=2.0),
        ],
    )
    recovery.execute(args(run_dir, target))
    attempt = run_dir / "attempt-002-cdiis-finalization.json"
    attempt.write_text(attempt.read_text() + "\n")
    item = evidence()
    with pytest.raises(ValueError, match="attempt receipt hash drift"):
        a2.checkpoint_converged_energy(
            target, item.cluster, item.settings, a2.PRODUCTION_METHOD
        )


def test_recovery_cache_rejects_rehashed_solver_semantic_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run_dir = tmp_path / "run"
    target = tmp_path / "reactant.wb97m-v.energy.json"
    install_fakes(
        monkeypatch,
        [
            FakeScf(converged=False, energy=-10.0, density=1.0),
            FakeScf(converged=True, energy=-10.25, density=2.0),
        ],
    )
    recovery.execute(args(run_dir, target))
    attempt = run_dir / "attempt-002-cdiis-finalization.json"
    attempt_payload = json.loads(attempt.read_text())
    attempt_payload["solver"]["algorithm"] = "wrong-solver"
    a2.atomic_json(attempt, attempt_payload)
    target_payload = json.loads(target.read_text())
    target_payload["recovery_attempt_receipts"][1]["sha256"] = a2.sha256_path(attempt)
    a2.atomic_json(target, target_payload)
    item = evidence()
    with pytest.raises(ValueError, match="solver contract drift"):
        a2.checkpoint_converged_energy(
            target, item.cluster, item.settings, a2.PRODUCTION_METHOD
        )


def test_recovery_cache_rejects_rehashed_non_solver_semantic_drift(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run_dir = tmp_path / "run"
    target = tmp_path / "reactant.wb97m-v.energy.json"
    install_fakes(
        monkeypatch,
        [
            FakeScf(converged=False, energy=-10.0, density=1.0),
            FakeScf(converged=True, energy=-10.25, density=2.0),
        ],
    )
    recovery.execute(args(run_dir, target))
    attempt = run_dir / "attempt-001-huckel-roothaan.json"
    attempt_payload = json.loads(attempt.read_text())
    attempt_payload["initial_density"]["kind"] = "not-huckel"
    attempt_payload["bindings"]["legacy_terminal"] = {
        "path": "/tampered",
        "sha256": "0" * 64,
    }
    a2.atomic_json(attempt, attempt_payload)
    target_payload = json.loads(target.read_text())
    target_payload["recovery_attempt_receipts"][0]["sha256"] = a2.sha256_path(attempt)
    a2.atomic_json(target, target_payload)
    item = evidence()
    with pytest.raises(ValueError, match="binding drift|evidence hash drift"):
        a2.checkpoint_converged_energy(
            target, item.cluster, item.settings, a2.PRODUCTION_METHOD
        )


def test_systemd_backstop_finalizes_running_state_and_quarantines_outputs(
    tmp_path: Path,
) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    target = tmp_path / "reactant.wb97m-v.energy.json"
    target.write_text("partial\n")
    a2.atomic_json(
        run_dir / "run_status.json",
        {
            "recovery_version": recovery.RECOVERY_VERSION,
            "recovery_contract": recovery.RECOVERY_CONTRACT,
            "status": "running",
        },
    )
    invocation = args(run_dir, target)

    assert recovery.finalize_if_running(invocation) == 0
    assert not target.exists()
    terminal = json.loads((run_dir / "terminal-receipt.json").read_text())
    status = json.loads((run_dir / "run_status.json").read_text())
    assert terminal["outcome"] == "incomplete-computational-failure"
    assert terminal["error_type"] == "RecoveryTimeout"
    assert status["status"] == "failed"


def test_systemd_backstop_fails_closed_before_running_status_exists(
    tmp_path: Path,
) -> None:
    run_dir = tmp_path / "run"
    target = tmp_path / "reactant.wb97m-v.energy.json"
    target.write_text("partial\n")
    invocation = args(run_dir, target)

    assert recovery.finalize_if_running(invocation) == 0
    assert not target.exists()
    terminal = json.loads((run_dir / "terminal-receipt.json").read_text())
    assert terminal["outcome"] == "incomplete-computational-failure"
    assert (run_dir / "run_status.json").is_file()


def test_main_rejects_wall_clock_override(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        "sys.argv",
        ["a2b1", "--wall-seconds", str(recovery.DRIVER_WALL_SECONDS + 1)],
    )
    with pytest.raises(ValueError, match="fixed at 7200"):
        recovery.main()


def test_terminal_failure_cannot_be_replayed(
    tmp_path: Path,
) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    target = tmp_path / "reactant.wb97m-v.energy.json"
    a2.atomic_json(
        run_dir / "terminal-receipt.json",
        {"outcome": "incomplete-computational-failure"},
    )
    target.write_text("stale energy\n")
    (run_dir / "results.json").write_text("stale result\n")
    (run_dir / "store.sqlite").write_text("stale store\n")
    invocation = args(run_dir, target)
    with pytest.raises(RuntimeError, match="no retry is authorized"):
        recovery.execute(invocation)
    assert not target.exists()
    assert not (run_dir / "results.json").exists()
    assert not (run_dir / "store.sqlite").exists()

    target.write_text("stale energy again\n")
    (run_dir / "results.json").write_text("stale result again\n")
    (run_dir / "store.sqlite").write_text("stale store again\n")
    assert recovery.finalize_if_running(invocation) == 0
    assert not target.exists()
    assert not (run_dir / "results.json").exists()
    assert not (run_dir / "store.sqlite").exists()
