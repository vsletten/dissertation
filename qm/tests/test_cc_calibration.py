"""Cheap contracts for the open-source A2 coupled-cluster layer."""

from __future__ import annotations

import argparse
import json
import math
from types import SimpleNamespace

import pytest

from scripts import cc_calibration as cc


def receipt(engine: str, basis: str, scf: float, correlation: float) -> dict:
    return {
        "engine": engine,
        "basis": basis,
        "identity": {"engine": engine, "basis": basis},
        "input_sha256": "abc",
        "charge": 0,
        "spin_2s": 0,
        "scf_hartree": scf,
        "correlation_hartree": correlation,
        "total_hartree": scf + correlation,
    }


def test_frozen_core_orbitals_for_si_neutral_cluster():
    # H8O8Si2: 8 oxygen 1s orbitals + 2 silicon Ne cores (5 each).
    assert cc.frozen_core_orbitals(["H"] * 8 + ["O"] * 8 + ["Si"] * 2) == 18


def test_byteqc_frozen_core_df_uses_one_complete_auxiliary_block():
    with_df = SimpleNamespace(blockdim=240, get_naoaux=lambda: 1337)
    coupled_cluster = SimpleNamespace(with_df=with_df)

    assert cc.configure_byteqc_frozen_core_blocking(coupled_cluster) == 1337
    assert with_df.blockdim == 1337


def test_two_point_extrapolations_recover_synthetic_limits():
    hf_limit = -100.0
    hf_amplitude = 0.2
    tz_hf = hf_limit + hf_amplitude * math.exp(-cc.HF_CBS_ALPHA * 3)
    qz_hf = hf_limit + hf_amplitude * math.exp(-cc.HF_CBS_ALPHA * 4)
    assert cc.extrapolate_hf(tz_hf, qz_hf) == pytest.approx(hf_limit)

    corr_limit = -1.0
    corr_amplitude = 0.5
    tz_corr = corr_limit + corr_amplitude / 3**3
    qz_corr = corr_limit + corr_amplitude / 4**3
    assert cc.extrapolate_correlation(tz_corr, qz_corr) == pytest.approx(corr_limit)


def test_receipt_grid_requires_complete_role_basis_matrix():
    values = [
        "reactant=cc-pVTZ=r-tz.json",
        "reactant=cc-pVQZ=r-qz.json",
        "ts=cc-pVTZ=t-tz.json",
        "ts=cc-pVQZ=t-qz.json",
    ]
    grid = cc.receipt_grid(values)
    assert set(grid) == {"reactant", "ts"}
    assert set(grid["reactant"]) == {"tz", "qz"}
    with pytest.raises(ValueError, match="cover reactant/ts"):
        cc.receipt_grid(values[:-1])


def test_existing_receipt_fails_closed_on_identity_drift(tmp_path):
    path = tmp_path / "receipt.json"
    payload = receipt("byteqc-canonical-ccsd(t)", "cc-pVTZ", -10.0, -1.0)
    path.write_text(json.dumps(payload))
    assert cc.existing_receipt(
        path,
        engine="byteqc-canonical-ccsd(t)",
        basis="cc-pVTZ",
        input_sha256="abc",
        identity=payload["identity"],
    )
    with pytest.raises(ValueError, match="identity drift"):
        cc.existing_receipt(
            path,
            engine="byteqc-canonical-ccsd(t)",
            basis="cc-pVQZ",
            input_sha256="abc",
            identity=payload["identity"],
        )
    with pytest.raises(ValueError, match="expected basis cc-pVQZ"):
        cc.load_receipt(path, "byteqc-canonical-ccsd(t)", "cc-pVQZ")


def test_calibration_engines_refuse_open_shell_before_import():
    with pytest.raises(ValueError, match="closed-shell spin=0 only"):
        cc.byteqc_job(argparse.Namespace(spin=1))
    with pytest.raises(ValueError, match="supports spin=0 only"):
        cc.psi4_job(argparse.Namespace(spin=1))


def test_calibration_resource_contract_requires_exact_byteqc_gpu_lease():
    byteqc = argparse.Namespace(
        command="byteqc", gpu=False, gpu_mem_gb=16.0, gpu_memory_gb=16
    )
    with pytest.raises(ValueError, match="requires --gpu"):
        cc.validate_resource_contract(byteqc)

    byteqc.gpu = True
    byteqc.gpu_mem_gb = 12.0
    with pytest.raises(ValueError, match="must be at least"):
        cc.validate_resource_contract(byteqc)

    byteqc.gpu_mem_gb = 18.0
    cc.validate_resource_contract(byteqc)

    with pytest.raises(ValueError, match="only for ByteQC"):
        cc.validate_resource_contract(argparse.Namespace(command="psi4", gpu=True))


def test_summarize_reports_barrier_delta_and_gate(tmp_path):
    paths = {engine: {} for engine in ("canonical", "dlpno")}
    engine_names = {
        "canonical": "byteqc-canonical-ccsd(t)",
        "dlpno": "psi4-dlpno-ccsd(t)",
    }
    for engine, engine_name in engine_names.items():
        for role, offset in (("reactant", 0.0), ("ts", 0.04)):
            for basis, label in (("cc-pVTZ", "tz"), ("cc-pVQZ", "qz")):
                # DLPNO TS is 0.0001 Eh above canonical (~0.26 kJ/mol delta).
                local_delta = 0.0001 if engine == "dlpno" and role == "ts" else 0.0
                payload = receipt(
                    engine_name,
                    basis,
                    -100.0 + offset + local_delta,
                    -1.0,
                )
                path = tmp_path / f"{engine}-{role}-{label}.json"
                path.write_text(json.dumps(payload))
                paths[engine].setdefault(role, {})[label] = path
    output = tmp_path / "summary.json"
    result = cc.summarize(
        argparse.Namespace(
            canonical=paths["canonical"],
            dlpno=paths["dlpno"],
            output=output,
        )
    )
    assert result["canonical_barrier_kj"] == pytest.approx(0.04 * cc.HARTREE_TO_KJ)
    assert result["dlpno_minus_canonical_kj"] == pytest.approx(
        0.0001 * cc.HARTREE_TO_KJ
    )
    assert result["gate_pass"] is True
    assert json.loads(output.read_text()) == result

    drift_path = paths["dlpno"]["ts"]["qz"]
    drift = json.loads(drift_path.read_text())
    drift["input_sha256"] = "different-geometry"
    drift_path.write_text(json.dumps(drift))
    with pytest.raises(ValueError, match="do not share one exact geometry"):
        cc.summarize(
            argparse.Namespace(
                canonical=paths["canonical"],
                dlpno=paths["dlpno"],
                output=output,
            )
        )
