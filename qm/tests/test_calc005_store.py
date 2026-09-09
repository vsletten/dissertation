"""CALC-005 temporary Store evidence-graph tests."""

from __future__ import annotations

import hashlib
import json
import sqlite3
from pathlib import Path

import pytest

from quarry.calc005 import build_calc005_pair
from quarry.calc005_store import (
    calculation_receipt_sha256,
    receipt_payload_sha256,
    validate_calc005_store,
    write_calc005_store,
)
from scripts import calc005_si_attachment as driver

REPO = Path(__file__).resolve().parents[2]
DECK = REPO / "petra/examples/kaolinite.toml"


def _accepted(tmp_path):
    backend = driver.AnalyticalBackend()
    receipt = driver.run_pilot(
        tmp_path, DECK, backend, executor_identity="store-test", source_commit="f" * 40
    )
    assert receipt["verdict"] == "passed-protocol-pilot"
    return receipt


def _store_path(receipt):
    return Path(receipt["artifacts"]["store"]["path"])


def _file_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _coordinated_receipt_tamper(receipt, mutator, *, repair_checkpoint=True):
    calculation_path = Path(receipt["artifacts"]["calculation_receipt"]["path"])
    calculation = json.loads(calculation_path.read_text())
    mutator(calculation)
    component = calculation["components"]["C"]
    if repair_checkpoint:
        component["checkpoint_sha256"] = driver._checkpoint_sha256(component)
    calculation_path.write_text(
        json.dumps(calculation, indent=2, sort_keys=True) + "\n"
    )

    store_path = _store_path(receipt)
    with sqlite3.connect(store_path) as connection:
        connection.execute(
            "UPDATE calc005_cycles SET calculation_receipt_sha256=?",
            (calculation_receipt_sha256(calculation),),
        )
        if repair_checkpoint:
            connection.execute(
                "UPDATE calc005_components SET checkpoint_sha256=? WHERE role='C'",
                (component["checkpoint_sha256"],),
            )
        connection.commit()

    receipt["artifacts"]["calculation_receipt"]["sha256"] = _file_sha256(
        calculation_path
    )
    receipt["artifacts"]["calculation_receipt"]["canonical_sha256"] = (
        calculation_receipt_sha256(calculation)
    )
    receipt["artifacts"]["store"]["sha256"] = _file_sha256(store_path)
    receipt["receipt_payload_sha256"] = receipt_payload_sha256(receipt)
    return receipt


def test_store_has_exact_cycle_component_construction_and_job_graph(tmp_path):
    receipt = _accepted(tmp_path)
    report = validate_calc005_store(
        _store_path(receipt),
        receipt,
        build_calc005_pair(DECK, environment_index=1),
    )

    assert report["roles"] == ["C", "SiOH4", "V"]
    assert report["coefficients"] == {"C": -1, "V": 1, "SiOH4": 1}
    assert report["job_count"] == 9
    assert report["construction_water_atoms"] == 3
    assert report["status"] == "passed-protocol-pilot"


@pytest.mark.parametrize(
    "sql",
    [
        "DELETE FROM calc005_components WHERE role='V'",
        "UPDATE calc005_components SET coefficient=9 WHERE role='SiOH4'",
        "UPDATE jobs SET status='failed' WHERE id="
        "(SELECT opt_job_id FROM calc005_components WHERE role='C')",
        "UPDATE results SET units='eV' WHERE key='electronic_hartree'",
        "UPDATE results SET value=value+1 WHERE key='zpe' AND job_id="
        "(SELECT freq_job_id FROM calc005_components WHERE role='C')",
        "UPDATE calc005_cycles SET label='kinetic barrier'",
        "DELETE FROM calc005_construction WHERE origin_kind LIKE 'hydrolysis-water-%'",
        "INSERT INTO jobs(structure_id,kind,method,engine,status,created_at) "
        "VALUES(1,'sp','orphan','fake','done','2026-01-01T00:00:00+00:00')",
    ],
)
def test_read_only_validator_refuses_sabotaged_graph(tmp_path, sql):
    receipt = _accepted(tmp_path)
    path = _store_path(receipt)
    with sqlite3.connect(path) as connection:
        connection.execute(sql)
        connection.commit()
    receipt["artifacts"]["store"]["sha256"] = _file_sha256(path)
    receipt["receipt_payload_sha256"] = receipt_payload_sha256(receipt)

    with pytest.raises(RuntimeError):
        validate_calc005_store(
            path, receipt, build_calc005_pair(DECK, environment_index=1)
        )


def test_writer_refuses_nonterminal_or_kinetic_payload(tmp_path):
    pair = build_calc005_pair(DECK, environment_index=1)
    with pytest.raises(RuntimeError, match="passed-protocol-pilot"):
        write_calc005_store(tmp_path / "bad.sqlite", pair, {}, {"verdict": "failed"})


@pytest.mark.parametrize(
    ("defect", "message"),
    [
        ("gate-rejected", "status is not accepted"),
        ("optimizer-unconverged", "optimizer budget/convergence"),
        ("optimizer-retry", "optimizer budget/convergence"),
        ("gradient-threshold", "gradient gate failed"),
        ("frequency-incomplete", "frequency gate/count"),
        ("frequency-imaginary", "frequency gate/count"),
        ("sp-unconverged", "single point failed"),
        ("checkpoint-hash", "checkpoint hash drifted"),
    ],
)
def test_store_rederives_acceptance_gates_despite_coordinated_manifest(
    tmp_path, defect, message
):
    receipt = _accepted(tmp_path)

    def mutate(calculation):
        component = calculation["components"]["C"]
        if defect == "gate-rejected":
            component["status"] = "gate-rejected"
        elif defect == "optimizer-unconverged":
            component["optimizer"]["converged"] = False
        elif defect == "optimizer-retry":
            component["optimizer"]["observed_calls"] = 2
            component["optimizer"]["observed_retries"] = 1
        elif defect == "gradient-threshold":
            component["gradient"]["rms_hartree_per_bohr"] = 3.1e-4
        elif defect == "frequency-incomplete":
            component["frequency"]["frequencies_cm"].pop()
            component["frequency"]["real_mode_count"] -= 1
            component["frequency"]["observed_mode_count"] -= 1
        elif defect == "frequency-imaginary":
            component["frequency"]["frequencies_cm"].pop()
            component["frequency"]["real_mode_count"] -= 1
            component["frequency"]["imaginary_cm"].append(30.0001)
            component["frequency"]["imaginary_mode_count"] += 1
        elif defect == "sp-unconverged":
            component["production_single_point"]["converged"] = False
        elif defect == "checkpoint-hash":
            component["gradient"]["max_hartree_per_bohr"] = 1.0e-12
        else:  # pragma: no cover - guarded by the parameter list
            raise AssertionError(defect)

    coordinated = _coordinated_receipt_tamper(
        receipt,
        mutate,
        repair_checkpoint=defect != "checkpoint-hash",
    )
    with pytest.raises(RuntimeError, match=message):
        validate_calc005_store(
            _store_path(coordinated),
            coordinated,
            build_calc005_pair(DECK, environment_index=1),
        )


def test_store_recomputes_thermochemistry_from_raw_modes(tmp_path):
    receipt = _accepted(tmp_path)

    def coordinated_thermo_lie(calculation):
        calculation["components"]["C"]["thermochemistry"]["zpe_kj_mol"] += 1.0

    coordinated = _coordinated_receipt_tamper(receipt, coordinated_thermo_lie)
    with pytest.raises(RuntimeError, match="recompute from raw modes"):
        validate_calc005_store(
            _store_path(coordinated),
            coordinated,
            build_calc005_pair(DECK, environment_index=1),
        )


def test_store_rejects_forged_production_envelope_even_when_rehashed(tmp_path):
    receipt = _accepted(tmp_path)

    def forged_envelope(calculation):
        calculation["execution_envelope"] = {
            "measured": True,
            "runtime_max_seconds": 43200,
            "qi2_lease": {"owner": "somebody-else", "ttl_hours": 13.0},
        }

    coordinated = _coordinated_receipt_tamper(receipt, forged_envelope)
    with pytest.raises(RuntimeError, match="execution envelope"):
        validate_calc005_store(
            _store_path(coordinated),
            coordinated,
            build_calc005_pair(DECK, environment_index=1),
        )


def test_store_rejects_production_backend_labeled_in_process(tmp_path):
    receipt = _accepted(tmp_path)

    def production_backend_lie(calculation):
        calculation["backend_kind"] = "production"

    coordinated = _coordinated_receipt_tamper(receipt, production_backend_lie)
    with pytest.raises(RuntimeError, match="execution envelope"):
        validate_calc005_store(
            _store_path(coordinated),
            coordinated,
            build_calc005_pair(DECK, environment_index=1),
        )
