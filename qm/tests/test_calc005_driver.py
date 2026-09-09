"""Fail-closed, calculator-free tests for the CALC-005 production driver."""

from __future__ import annotations

import ast
import json
import math
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from quarry.calc005 import compose_thermochemistry
from quarry.pipeline import (
    FrequencyResult,
    OptimizationResult,
    frequency_geometry_fingerprint,
)
from scripts import calc005_si_attachment as driver
from scripts import calc005_si_attachment_verify as verifier

REPO = Path(__file__).resolve().parents[2]
DECK = REPO / "petra/examples/kaolinite.toml"


def _calculation(receipt):
    return json.loads(
        Path(receipt["artifacts"]["calculation_receipt"]["path"]).read_text()
    )


def _production_envelope(pid=4242):
    return {
        "measured": True,
        "systemd_unit": "calc005.service",
        "cgroup_path": "/sys/fs/cgroup/system.slice/calc005.service",
        "runtime_max_seconds": 43200,
        "memory_max_bytes": 32 * 1024**3,
        "memory_swap_max_bytes": 4 * 1024**3,
        "cpu_quota_percent": 1600,
        "nice": 10,
        "thread_environment": {
            "OMP_NUM_THREADS": 16,
            "MKL_NUM_THREADS": 16,
            "OPENBLAS_NUM_THREADS": 16,
        },
        "qi2_lease": {
            "path": "/tmp/lease.json",
            "owner": driver.EXPECTED_GPU_OWNER,
            "pid": pid,
            "ttl_hours": 13.0,
            "expected_gb": 16.0,
            "maximum_gb": 18.0,
        },
        "shared_service_mutation": False,
        "restoration": "pending-process-exit",
    }


class FakeBackend:
    def __init__(
        self, *, converged: bool = True, imaginary: dict[str, list[float]] | None = None
    ):
        self.converged = converged
        self.imaginary = imaginary or {}
        self.optimize_calls: list[tuple[str, int]] = []
        self.energy_calls: list[tuple[str, str]] = []

    def optimize(self, role, cluster, settings, *, max_steps):
        self.optimize_calls.append((role, max_steps))
        return OptimizationResult(
            cluster=cluster, converged=self.converged, max_steps=max_steps
        )

    def gradient(self, role, cluster, settings):
        return np.zeros_like(cluster.coords)

    def frequencies(self, role, cluster, settings):
        imaginary = np.asarray(self.imaginary.get(role, []), dtype=float)
        mode_count = driver.expected_mode_count(role, cluster)
        return FrequencyResult(
            frequencies_cm=np.linspace(100.0, 1700.0, mode_count - len(imaginary)),
            imaginary_cm=imaginary,
            electronic_hartree=-100.0,
            molar_mass_kg=0.1,
            rotational_temperatures_k=(1.0, 2.0, 3.0) if role == "SiOH4" else None,
            linear=False,
            geometry_fingerprint=frequency_geometry_fingerprint(cluster),
            settings_fingerprint=driver.frequency_settings_fingerprint(settings),
        )

    def energy(self, role, cluster, settings):
        self.energy_calls.append((role, settings.xc))
        return {"C": -300.0, "V": -200.0, "SiOH4": -100.02}[role]


def test_exact_settings_and_fixed_scope(tmp_path):
    geometry, production = driver.calc005_settings(use_gpu=True)
    assert driver.MAX_STEPS == 150
    assert geometry.xc == "r2scan"
    assert geometry.basis == "def2-mtzvpp"
    assert geometry.composite == "r2scan3c"
    assert geometry.solvent is None
    assert production.xc == "wb97m-v"
    assert production.basis == "def2-tzvpd"
    assert production.solvent == "smd"
    assert production.dispersion is None
    refused = tmp_path / "unused"
    with pytest.raises(ValueError, match="only.*i=1"):
        driver.run_pilot(refused, DECK, FakeBackend(), environment_index=2)
    assert not refused.exists()


def test_happy_path_persists_exact_components_terms_and_no_water(tmp_path):
    backend = FakeBackend()
    receipt = driver.run_pilot(
        tmp_path,
        DECK,
        backend,
        executor_identity="executor-a",
        source_commit="a" * 40,
    )

    assert receipt["verdict"] == "passed-protocol-pilot"
    assert receipt["classification"] == "thermodynamic/non-kinetic/non-emittable"
    calculation = _calculation(receipt)
    assert calculation["stoichiometry"] == {"C": -1, "V": 1, "SiOH4": 1}
    assert "H2O" not in json.dumps(calculation["thermochemistry"])
    assert calculation["thermochemistry"]["standard_state_kj_mol"] == pytest.approx(
        7.958500693927389, abs=1e-15
    )
    assert calculation["thermochemistry"]["recomposition_error_kj_mol"] <= 1e-8
    assert [role for role, steps in backend.optimize_calls] == ["C", "V", "SiOH4"]
    assert all(steps == 150 for _role, steps in backend.optimize_calls)
    assert calculation["optimizer_budget"] == {
        "C": {"calls": 1, "max_steps": 150, "retries": 0},
        "V": {"calls": 1, "max_steps": 150, "retries": 0},
        "SiOH4": {"calls": 1, "max_steps": 150, "retries": 0},
    }
    assert calculation["components"]["C"]["frequency"]["hessian"] == "PHVA"
    assert calculation["components"]["V"]["frequency"]["hessian"] == "PHVA"
    assert calculation["components"]["SiOH4"]["frequency"]["hessian"] == "full"
    assert (tmp_path / "terminal-receipt.json").is_file()
    assert not (tmp_path / "calc005-result.json").exists()
    assert not (tmp_path / "store.sqlite").exists()
    generation = tmp_path / "generations" / receipt["generation"]
    assert sorted(path.name for path in generation.iterdir()) == [
        "calc005-result.json",
        "calculation-receipt.json",
        "store.sqlite",
    ]
    for role in ("C", "V", "SiOH4"):
        component = tmp_path / "components" / role
        assert (component / "raw-endpoint.xyz").is_file()
        assert (component / "endpoint.xyz").is_file()
        assert (
            json.loads((component / "receipt.json").read_text())["status"] == "accepted"
        )


def test_run_and_validate_cli_surfaces_remain_backend_mockable(
    tmp_path, monkeypatch, capsys
):
    backend = FakeBackend()
    monkeypatch.setattr(driver, "PipelineBackend", lambda: backend)
    provenance = {
        "repository": str(REPO),
        "branch": "test",
        "head": "4" * 40,
        "origin_head": "4" * 40,
        "clean": True,
        "tracked_source_files": [],
    }
    monkeypatch.setattr(driver, "verify_production_source", lambda *_args: provenance)
    monkeypatch.setattr(driver, "measure_production_envelope", _production_envelope)

    run_status = driver.main(
        [
            "run",
            "--deck",
            str(DECK),
            "--output-root",
            str(tmp_path),
            "--gpu",
        ]
    )
    run_payload = json.loads(capsys.readouterr().out)
    validate_status = driver.main(
        [
            "validate",
            "--deck",
            str(DECK),
            "--output-root",
            str(tmp_path),
        ]
    )
    validate_payload = json.loads(capsys.readouterr().out)

    assert run_status == 0
    assert run_payload["verdict"] == "passed-protocol-pilot"
    assert validate_status == 0
    assert validate_payload["status"] == "valid"
    assert validate_payload["optimizer_calls"] == 0
    assert backend.optimize_calls == [("C", 150), ("V", 150), ("SiOH4", 150)]


def test_optimizer_exhaustion_is_terminal_zero_retry_and_preserves_raw(tmp_path):
    backend = FakeBackend(converged=False)
    receipt = driver.run_pilot(
        tmp_path, DECK, backend, executor_identity="executor-a", source_commit="b" * 40
    )

    assert receipt["verdict"] == "incomplete-computational-failure"
    assert len(backend.optimize_calls) == 1
    assert receipt["optimizer_budget"]["C"] == {
        "calls": 1,
        "max_steps": 150,
        "retries": 0,
    }
    assert (tmp_path / "components/C/raw-endpoint.xyz").is_file()
    assert not (tmp_path / "calc005-result.json").exists()
    assert not (tmp_path / "store.sqlite").exists()


def test_post_optimizer_calculator_failure_is_terminal_and_preserves_raw(tmp_path):
    class BrokenGradient(FakeBackend):
        def gradient(self, role, cluster, settings):
            raise RuntimeError("simulated gradient failure")

    backend = BrokenGradient()
    receipt = driver.run_pilot(
        tmp_path, DECK, backend, executor_identity="executor-a", source_commit="9" * 40
    )

    assert receipt["verdict"] == "incomplete-computational-failure"
    assert receipt["failed_component"] == "C"
    assert "simulated gradient failure" in receipt["detail"]
    assert backend.optimize_calls == [("C", 150)]
    assert (tmp_path / "components/C/raw-endpoint.xyz").is_file()
    assert not (tmp_path / "calc005-result.json").exists()
    assert not (tmp_path / "store.sqlite").exists()


@pytest.mark.parametrize("defect", ["gradient", "frequency-binding", "energy"])
def test_independent_numerical_evidence_fails_closed(tmp_path, defect):
    class DefectiveEvidence(FakeBackend):
        def gradient(self, role, cluster, settings):
            value = super().gradient(role, cluster, settings)
            if role == "C" and defect == "gradient":
                value[0, 0] = 1.0
            return value

        def frequencies(self, role, cluster, settings):
            value = super().frequencies(role, cluster, settings)
            if role == "C" and defect == "frequency-binding":
                return replace(value, geometry_fingerprint="wrong-endpoint")
            return value

        def energy(self, role, cluster, settings):
            if role == "C" and defect == "energy":
                return float("nan")
            return super().energy(role, cluster, settings)

    backend = DefectiveEvidence()
    receipt = driver.run_pilot(
        tmp_path, DECK, backend, executor_identity="executor-a", source_commit="8" * 40
    )

    assert receipt["verdict"] == "rejected-physical-state"
    assert backend.optimize_calls == [("C", 150)]
    assert (tmp_path / "components/C/raw-endpoint.xyz").is_file()
    assert not (tmp_path / "calc005-result.json").exists()
    assert not (tmp_path / "store.sqlite").exists()


@pytest.mark.parametrize("defect", ["owner", "nonfinite", "nonminimum"])
def test_physical_and_numerical_gates_reject_after_raw_persistence(tmp_path, defect):
    class Defective(FakeBackend):
        def optimize(self, role, cluster, settings, *, max_steps):
            result = super().optimize(role, cluster, settings, max_steps=max_steps)
            if role != "C":
                return result
            coords = cluster.coords.copy()
            if defect == "owner":
                hydrogen = cluster.symbols.index("H")
                oxygen = [
                    i for i, symbol in enumerate(cluster.symbols) if symbol == "O"
                ][-1]
                coords[hydrogen] = coords[oxygen] + np.array([0.96, 0.0, 0.0])
            elif defect == "nonfinite":
                coords[0, 0] = np.nan
            return replace(result, cluster=replace(cluster, coords=coords))

    backend = Defective(imaginary={"C": [45.0]} if defect == "nonminimum" else None)
    receipt = driver.run_pilot(
        tmp_path, DECK, backend, executor_identity="executor-a", source_commit="c" * 40
    )

    assert receipt["verdict"] == "rejected-physical-state"
    assert (tmp_path / "components/C/raw-endpoint.xyz").is_file()
    assert not (tmp_path / "calc005-result.json").exists()
    assert not (tmp_path / "store.sqlite").exists()


def test_corrupt_or_spent_checkpoint_refuses_without_another_optimizer(tmp_path):
    component = tmp_path / "components/C"
    component.mkdir(parents=True)
    (component / "reservation.json").write_text("{not-json")
    backend = FakeBackend()

    receipt = driver.run_pilot(
        tmp_path, DECK, backend, executor_identity="executor-a", source_commit="d" * 40
    )

    assert receipt["verdict"] == "incomplete-computational-failure"
    assert backend.optimize_calls == []
    assert receipt["failed_component"] == "C"
    assert receipt["quarantine"]


def test_forged_success_terminal_is_quarantined_without_optimization(tmp_path):
    tmp_path.joinpath("terminal-receipt.json").write_text(
        json.dumps(
            {
                "schema": "calc005-terminal-receipt-v1",
                "verdict": "passed-protocol-pilot",
                "source_commit": "5" * 40,
            }
        )
    )
    backend = FakeBackend()

    receipt = driver.run_pilot(
        tmp_path, DECK, backend, executor_identity="executor-a", source_commit="5" * 40
    )

    assert receipt["verdict"] == "incomplete-computational-failure"
    assert receipt["failed_component"] == "preflight"
    assert "terminal receipt refused" in receipt["detail"]
    assert backend.optimize_calls == []
    assert Path(receipt["quarantine"]).is_dir()


def test_stale_accepted_checkpoint_is_not_reused_across_source_commits(tmp_path):
    first = FakeBackend()
    driver.run_pilot(
        tmp_path,
        DECK,
        first,
        executor_identity="executor-a",
        source_commit="a" * 40,
    )
    (tmp_path / "terminal-receipt.json").unlink()
    second = FakeBackend()

    receipt = driver.run_pilot(
        tmp_path,
        DECK,
        second,
        executor_identity="executor-a",
        source_commit="b" * 40,
    )

    assert receipt["verdict"] == "incomplete-computational-failure"
    assert "identity drifted" in receipt["detail"]
    assert second.optimize_calls == []
    assert receipt["quarantine"] == "spent-checkpoint-retained-in-place"
    assert (tmp_path / "components/C/optimizer-entered.json").is_file()


def test_tampered_accepted_checkpoint_is_terminal_without_reoptimization(tmp_path):
    first = FakeBackend()
    driver.run_pilot(
        tmp_path,
        DECK,
        first,
        executor_identity="executor-a",
        source_commit="7" * 40,
    )
    (tmp_path / "terminal-receipt.json").unlink()
    checkpoint = tmp_path / "components/C/receipt.json"
    payload = json.loads(checkpoint.read_text())
    payload["thermochemistry"]["zpe_kj_mol"] += 1.0
    checkpoint.write_text(json.dumps(payload))
    second = FakeBackend()

    receipt = driver.run_pilot(
        tmp_path, DECK, second, executor_identity="executor-a", source_commit="7" * 40
    )

    assert receipt["verdict"] == "incomplete-computational-failure"
    assert "payload drifted" in receipt["detail"]
    assert second.optimize_calls == []


def test_optimizer_interruption_is_recorded_and_never_replayed(tmp_path):
    class Interrupted(FakeBackend):
        def optimize(self, role, cluster, settings, *, max_steps):
            self.optimize_calls.append((role, max_steps))
            raise KeyboardInterrupt("simulated stop")

    first = Interrupted()
    with pytest.raises(KeyboardInterrupt, match="simulated stop"):
        driver.run_pilot(
            tmp_path,
            DECK,
            first,
            executor_identity="executor-a",
            source_commit="e" * 40,
        )
    reservation = tmp_path / "components/C/reservation.json"
    assert reservation.is_file()

    second = FakeBackend()
    receipt = driver.run_pilot(
        tmp_path, DECK, second, executor_identity="executor-a", source_commit="e" * 40
    )
    assert receipt["verdict"] == "incomplete-computational-failure"
    assert second.optimize_calls == []


def test_graph_publication_failure_never_promotes_a_canonical_value(
    tmp_path, monkeypatch
):
    def reject_store(*_args, **_kwargs):
        raise RuntimeError("simulated Store failure")

    monkeypatch.setattr(driver, "write_calc005_store", reject_store)
    backend = FakeBackend()

    receipt = driver.run_pilot(
        tmp_path, DECK, backend, executor_identity="executor-a", source_commit="6" * 40
    )

    assert receipt["verdict"] == "incomplete-computational-failure"
    assert receipt["failed_component"] == "evidence-graph-publication"
    assert len(backend.optimize_calls) == 3
    assert not (tmp_path / "calc005-result.json").exists()
    assert not (tmp_path / "store.sqlite").exists()


def test_cold_verifier_recomputes_every_gate_without_optimizing(tmp_path):
    backend = FakeBackend()
    driver.run_pilot(
        tmp_path,
        DECK,
        backend,
        executor_identity="executor-a",
        source_commit="1" * 40,
    )
    before = list(backend.optimize_calls)

    result = verifier.verify_pilot(
        tmp_path,
        DECK,
        backend,
        verifier_identity="different-worker",
    )

    assert result["status"] == "verified-pass"
    assert result["optimizer_calls"] == 0
    assert backend.optimize_calls == before
    assert result["petra_boundary"]["production_value_emitted"] is False
    assert result["detailed_balance"]["equation"].endswith("exp(S_10/RT)")


def test_cold_verifier_fails_closed_on_endpoint_tampering(tmp_path):
    backend = FakeBackend()
    driver.run_pilot(
        tmp_path,
        DECK,
        backend,
        executor_identity="executor-a",
        source_commit="2" * 40,
    )
    endpoint = tmp_path / "components/C/endpoint.xyz"
    endpoint.write_text(endpoint.read_text().replace("Al ", "Si ", 1))

    result = verifier.verify_pilot(
        tmp_path,
        DECK,
        backend,
        verifier_identity="different-worker",
    )

    assert result["status"] == "verified-fail"
    assert "hash" in result["detail"].lower()
    assert backend.optimize_calls == [("C", 150), ("V", 150), ("SiOH4", 150)]


def test_structural_gate_rejects_moved_center_si_and_destroyed_heavy_topology(tmp_path):
    class BrokenTopology(FakeBackend):
        def optimize(self, role, cluster, settings, *, max_steps):
            result = super().optimize(role, cluster, settings, max_steps=max_steps)
            if role == "C":
                coords = cluster.coords.copy()
                center = cluster.symbols.index("Si")
                coords[center] += np.array([8.0, 8.0, 8.0])
                return replace(result, cluster=replace(cluster, coords=coords))
            return result

    receipt = driver.run_pilot(
        tmp_path,
        DECK,
        BrokenTopology(),
        executor_identity="executor-a",
        source_commit="3" * 40,
    )
    assert receipt["verdict"] == "rejected-physical-state"
    assert "heavy connectivity" in receipt["detail"]


@pytest.mark.parametrize("missing", [1, 21])
def test_frequency_gate_rejects_empty_or_incomplete_mode_vectors(tmp_path, missing):
    class IncompleteModes(FakeBackend):
        def frequencies(self, role, cluster, settings):
            result = super().frequencies(role, cluster, settings)
            if role == "C":
                return replace(result, frequencies_cm=result.frequencies_cm[:-missing])
            return result

    receipt = driver.run_pilot(
        tmp_path,
        DECK,
        IncompleteModes(),
        executor_identity="executor-a",
        source_commit="4" * 40,
    )
    assert receipt["verdict"] == "rejected-physical-state"
    assert "mode count" in receipt["detail"]


def test_frequency_mode_count_accounts_for_below_noise_imaginary_modes(tmp_path):
    backend = FakeBackend(imaginary={"C": [12.0]})
    receipt = driver.run_pilot(
        tmp_path,
        DECK,
        backend,
        executor_identity="executor-a",
        source_commit="5" * 40,
    )
    frequency = _calculation(receipt)["components"]["C"]["frequency"]
    assert frequency["imaginary_mode_count"] == 1
    assert frequency["real_mode_count"] + frequency["imaginary_mode_count"] == 165


def test_compose_thermochemistry_refuses_nonprotocol_temperature():
    terms = {
        role: {
            "electronic_kj_mol": 0.0,
            "zpe_kj_mol": 0.0,
            "thermal_kj_mol": 0.0,
            "entropy_kj_mol_k": 0.0,
        }
        for role in ("C", "V", "SiOH4")
    }
    with pytest.raises(ValueError, match="exactly 298.15"):
        compose_thermochemistry(terms, temperature_k=298.15000000000003)


def test_precall_reservation_crash_is_recoverable_without_spending(
    tmp_path, monkeypatch
):
    original = driver.atomic_json
    interrupted = {"done": False}

    def crash_before_entered(path, payload):
        if path.name == "optimizer-entered.json" and not interrupted["done"]:
            interrupted["done"] = True
            raise KeyboardInterrupt("before entered marker promotion")
        original(path, payload)

    monkeypatch.setattr(driver, "atomic_json", crash_before_entered)
    first = FakeBackend()
    with pytest.raises(KeyboardInterrupt, match="before entered"):
        driver.run_pilot(
            tmp_path,
            DECK,
            first,
            executor_identity="executor-a",
            source_commit="6" * 40,
        )
    assert first.optimize_calls == []
    assert (tmp_path / "components/C/reservation.json").is_file()
    assert not (tmp_path / "components/C/optimizer-entered.json").exists()

    second = FakeBackend()
    receipt = driver.run_pilot(
        tmp_path,
        DECK,
        second,
        executor_identity="executor-a",
        source_commit="6" * 40,
    )
    assert receipt["verdict"] == "passed-protocol-pilot"
    assert second.optimize_calls == [("C", 150), ("V", 150), ("SiOH4", 150)]


def test_crash_after_optimizer_entered_marker_spends_budget(tmp_path):
    class CrashAfterEntered(FakeBackend):
        def optimize(self, role, cluster, settings, *, max_steps):
            self.optimize_calls.append((role, max_steps))
            raise KeyboardInterrupt("process died after optimizer entry")

    first = CrashAfterEntered()
    with pytest.raises(KeyboardInterrupt, match="after optimizer entry"):
        driver.run_pilot(
            tmp_path,
            DECK,
            first,
            executor_identity="executor-a",
            source_commit="6" * 40,
        )
    assert first.optimize_calls == [("C", 150)]
    assert (tmp_path / "components/C/optimizer-entered.json").is_file()

    second = FakeBackend()
    receipt = driver.run_pilot(
        tmp_path,
        DECK,
        second,
        executor_identity="executor-a",
        source_commit="6" * 40,
    )
    assert receipt["verdict"] == "incomplete-computational-failure"
    assert "already spent" in receipt["detail"]
    assert second.optimize_calls == []


def test_validate_cli_returns_nonzero_for_rejected_and_malformed_receipts(
    tmp_path, capsys
):
    rejected = driver.run_pilot(
        tmp_path,
        DECK,
        FakeBackend(converged=False),
        executor_identity="executor-a",
        source_commit="7" * 40,
    )
    assert rejected["verdict"] == "incomplete-computational-failure"
    assert (
        driver.main(["validate", "--deck", str(DECK), "--output-root", str(tmp_path)])
        == 1
    )
    assert json.loads(capsys.readouterr().out)["status"] == "invalid"
    (tmp_path / "terminal-receipt.json").write_text("{malformed")
    assert (
        driver.main(["validate", "--deck", str(DECK), "--output-root", str(tmp_path)])
        == 1
    )
    assert json.loads(capsys.readouterr().out)["status"] == "invalid"


def test_production_source_refuses_unrelated_repo_dirty_state_and_origin_drift(
    tmp_path, monkeypatch
):
    with pytest.raises(RuntimeError, match="this driver's repository"):
        driver.verify_production_source(tmp_path, DECK)

    monkeypatch.setattr(
        driver, "_git", lambda _repo, *args: "dirty" if args[0] == "status" else ""
    )
    with pytest.raises(RuntimeError, match="dirty"):
        driver.verify_production_source(REPO, DECK)

    def divergent(_repo, *args):
        if args[:2] == ("status", "--porcelain=v1"):
            return ""
        if args[:2] == ("branch", "--show-current"):
            return "feature"
        if args[-1] == "HEAD":
            return "a" * 40
        if args[-1] == "origin/feature":
            return "b" * 40
        return ""

    monkeypatch.setattr(driver, "_git", divergent)
    with pytest.raises(RuntimeError, match="origin"):
        driver.verify_production_source(REPO, DECK)


def test_verifier_has_no_executor_import_or_forbidden_helper_calls():
    source = Path(verifier.__file__).read_text()
    tree = ast.parse(source)
    forbidden = {
        "validate_run",
        "read_xyz",
        "structural_gate",
        "gradient_gate",
        "_frequency_payload",
        "_thermochemistry",
        "compose_thermochemistry",
    }
    imported = {
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, (ast.Import, ast.ImportFrom))
        for alias in node.names
    }
    called = {
        node.func.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    } | {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    assert "calc005_si_attachment" not in imported
    assert forbidden.isdisjoint(called)


def test_cold_verifier_uses_tolerances_for_one_ulp_numerical_differences(tmp_path):
    executor_backend = FakeBackend()
    driver.run_pilot(
        tmp_path,
        DECK,
        executor_backend,
        executor_identity="executor-a",
        source_commit="8" * 40,
    )

    class OneUlp(FakeBackend):
        def gradient(self, role, cluster, settings):
            value = super().gradient(role, cluster, settings)
            value[0, 0] = np.nextafter(0.0, 1.0)
            return value

        def frequencies(self, role, cluster, settings):
            value = super().frequencies(role, cluster, settings)
            frequencies = value.frequencies_cm.copy()
            frequencies[0] = np.nextafter(frequencies[0], math.inf)
            return replace(value, frequencies_cm=frequencies)

        def energy(self, role, cluster, settings):
            return np.nextafter(super().energy(role, cluster, settings), math.inf)

    result = verifier.verify_pilot(
        tmp_path, DECK, OneUlp(), verifier_identity="different-worker"
    )
    assert result["status"] == "verified-pass", result


def test_atomic_publication_exposes_no_partial_canonical_generation(
    tmp_path, monkeypatch
):
    monkeypatch.setattr(
        driver,
        "validate_staged_calc005_store",
        lambda *_args: (_ for _ in ()).throw(RuntimeError("staged validation failed")),
    )
    receipt = driver.run_pilot(
        tmp_path,
        DECK,
        FakeBackend(),
        executor_identity="executor-a",
        source_commit="9" * 40,
    )
    assert receipt["verdict"] == "incomplete-computational-failure"
    assert receipt["failed_component"] == "evidence-graph-publication"
    generations = tmp_path / "generations"
    assert generations.is_dir()
    assert list(generations.iterdir()) == []
    assert not (tmp_path / "calc005-result.json").exists()
    assert not (tmp_path / "calculation-receipt.json").exists()
    assert not (tmp_path / "store.sqlite").exists()
    assert json.loads((tmp_path / "terminal-receipt.json").read_text())["verdict"] == (
        "incomplete-computational-failure"
    )


def test_production_envelope_records_actual_exact_limits_and_lease(
    tmp_path, monkeypatch
):
    lease_path = tmp_path / "lease.json"
    lease = {
        "owner": driver.EXPECTED_GPU_OWNER,
        "pid": 4242,
        "started": "2026-09-08T00:00:00Z",
        "ttl": 13.0,
        "expected_gb": 16.0,
    }
    real_read_text = Path.read_text

    def measured_read_text(path, *args, **kwargs):
        value = str(path)
        if value == "/proc/self/cgroup":
            return "0::/system.slice/calc005.service\n"
        if value.endswith("/memory.max"):
            return f"{driver.EXPECTED_MEMORY_MAX_BYTES}\n"
        if value.endswith("/memory.swap.max"):
            return f"{driver.EXPECTED_MEMORY_SWAP_MAX_BYTES}\n"
        if value.endswith("/cpu.max"):
            return "1600000 100000\n"
        if path == lease_path:
            return json.dumps(lease)
        return real_read_text(path, *args, **kwargs)

    monkeypatch.setattr(Path, "read_text", measured_read_text)
    monkeypatch.setattr(
        driver.subprocess,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(stdout="12h\n"),
    )
    monkeypatch.setattr(driver.os, "getpriority", lambda *_args: 10)
    monkeypatch.setattr(driver.os, "getpid", lambda: 4242)
    monkeypatch.setenv("GPU_LEASE_PATH", str(lease_path))
    for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        monkeypatch.setenv(name, "16")

    measured = driver.measure_production_envelope()
    assert measured["runtime_max_seconds"] == 43200
    assert measured["memory_max_bytes"] == 32 * 1024**3
    assert measured["memory_swap_max_bytes"] == 4 * 1024**3
    assert measured["cpu_quota_percent"] == 1600
    assert measured["nice"] == 10
    assert measured["thread_environment"] == {
        "OMP_NUM_THREADS": 16,
        "MKL_NUM_THREADS": 16,
        "OPENBLAS_NUM_THREADS": 16,
    }
    assert measured["qi2_lease"]["owner"] == driver.EXPECTED_GPU_OWNER
    assert measured["qi2_lease"]["ttl_hours"] == 13.0

    lease["ttl"] = 12.999
    with pytest.raises(RuntimeError, match="lease_ttl"):
        driver.measure_production_envelope()


def test_bootstrap_source_pins_exact_thirteen_hour_gpu_ttl():
    source = Path(driver.__file__).read_text()
    tree = ast.parse(source)
    assignments = {
        target.id: node.value.value
        for node in tree.body
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Constant)
        for target in node.targets
        if isinstance(target, ast.Name)
    }
    assert assignments["GPU_TTL_HOURS"] == 13.0
    bootstrap = next(
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "bootstrap_cli"
    )
    ttl = next(
        keyword.value
        for keyword in bootstrap.keywords
        if keyword.arg == "gpu_ttl_hours"
    )
    assert isinstance(ttl, ast.Name)
    assert ttl.id == "GPU_TTL_HOURS"


def test_success_manifest_has_no_store_or_terminal_back_edge(tmp_path):
    receipt = driver.run_pilot(
        tmp_path,
        DECK,
        FakeBackend(),
        executor_identity="executor-a",
        source_commit="a" * 40,
    )
    calculation = _calculation(receipt)
    assert "artifacts" not in calculation
    assert "store_sha256" not in calculation
    assert "terminal_sha256" not in calculation
    assert receipt["artifacts"]["calculation_receipt"]["canonical_sha256"]
    assert receipt["artifacts"]["store"]["sha256"]
    assert receipt["receipt_payload_sha256"] == driver.receipt_payload_sha256(receipt)
