from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import pytest

SCRIPTS = Path(__file__).parents[1] / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
MODULE_PATH = SCRIPTS / "e3b_periodic_dft.py"
SPEC = importlib.util.spec_from_file_location("e3b_periodic_dft", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
e3b = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = e3b
SPEC.loader.exec_module(e3b)

WORKBOOK = Path(
    "/mnt/data/vsletten/dissertation-data/e3a-classical-neb-barriers-20260902/"
    "sources/nteme-2022-supplement/Results Of MD Simulations.xlsx"
)
CAMPAIGN = Path(
    "/mnt/data/vsletten/dissertation-data/e3a-classical-neb-campaign-grid-20260903"
)


def method_files(tmp_path: Path) -> Any:
    basis = tmp_path / "BASIS_MOLOPT"
    potential = tmp_path / "GTH_POTENTIALS"
    d3 = tmp_path / "dftd3.dat"
    basis.write_text("test basis data\n", encoding="utf-8")
    potential.write_text("test potential data\n", encoding="utf-8")
    d3.write_text("test D3 parameter data\n", encoding="utf-8")
    return e3b.MethodProvenance(
        cp2k_version="2026.1",
        cp2k_image_digest="sha256:" + "a" * 64,
        basis_set_file=basis,
        potential_file=potential,
        d3_parameter_file=d3,
    )


def test_require_sha256_fails_closed(tmp_path: Path) -> None:
    path = tmp_path / "input.json"
    path.write_text("{}\n", encoding="utf-8")
    with pytest.raises(ValueError, match="hash mismatch"):
        e3b.require_sha256(path, "0" * 64, label="campaign result")


def test_cp2k_smoke_status_is_fail_closed() -> None:
    complete = e3b.classify_cp2k_output(
        "SCF run converged in 12 steps\n"
        "ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: -123.5\n"
        "PROGRAM ENDED AT 2026-09-14\n",
        returncode=0,
        timed_out=False,
    )
    assert complete["status"] == "converged"
    assert complete["energy_hartree"] == -123.5

    assert (
        e3b.classify_cp2k_output(
            "SCF run NOT converged\n", returncode=0, timed_out=False
        )["status"]
        == "incomplete-scf"
    )
    assert (
        e3b.classify_cp2k_output("", returncode=None, timed_out=True)["status"]
        == "incomplete-timeout"
    )
    assert (
        e3b.classify_cp2k_output(
            "SCF run converged in 2 steps\n", returncode=1, timed_out=False
        )["status"]
        == "incomplete-execution"
    )


def test_smoke_runner_writes_timing_receipt(tmp_path: Path) -> None:
    fake = tmp_path / "fake-cp2k"
    fake.write_text(
        "#!/usr/bin/env python3\n"
        "print('SCF run converged in 3 steps')\n"
        "print('ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: -42.25')\n"
        "print('PROGRAM ENDED AT 2026-09-14')\n",
        encoding="utf-8",
    )
    fake.chmod(0o755)
    input_path = tmp_path / "smoke.inp"
    input_path.write_text("&GLOBAL\n  RUN_TYPE ENERGY_FORCE\n&END GLOBAL\n")
    receipt_path = tmp_path / "smoke-result.json"

    receipt = e3b.run_smoke(
        executable=str(fake),
        input_path=input_path,
        receipt_path=receipt_path,
        timeout_seconds=5.0,
    )

    assert receipt["status"] == "converged"
    assert receipt["energy_hartree"] == -42.25
    assert receipt["elapsed_seconds"] >= 0.0
    assert receipt["input_sha256"] == e3b.sha256(input_path)
    assert receipt["stdout_sha256"] == e3b.sha256(tmp_path / "smoke.stdout.log")
    assert json.loads(receipt_path.read_text(encoding="utf-8")) == receipt


def _prepared_analysis_fixture() -> dict[str, Any]:
    return {
        "schema": "e3b-cp2k-preparation-v1",
        "models": {
            "model": {
                "atom_identity_sha256": "identity",
                "classical_reference": {
                    "barrier_kcal_mol": 100.0,
                    "status": "converged",
                },
            }
        },
    }


def _observation_fixture() -> dict[str, Any]:
    return {
        "schema": "e3b-cp2k-observations-v1",
        "models": {
            "model": {
                "atom_identity_sha256": "identity",
                "cp2k": {
                    "status": "converged",
                    "timed_out": False,
                    "scf_converged_all_images": True,
                    "endpoints_converged": True,
                    "neb_converged": True,
                    "barrier_kcal_mol": 112.0,
                },
                "matched_classical": {
                    "status": "converged",
                    "timed_out": False,
                    "neb_converged": True,
                    "barrier_kcal_mol": 102.0,
                },
            }
        },
    }


def test_analysis_refuses_unbound_manual_barriers() -> None:
    result = e3b.analyze_models(
        _prepared_analysis_fixture(),
        _observation_fixture(),
        transfer_tolerance_kcal_mol=5.0,
        endorsement_tolerance_kcal_mol=5.0,
    )
    model = result["models"]["model"]
    assert model["typed_outcome"] == "incomplete-unverified-observation"
    assert model["matched_cell_correction_kcal_mol"] is None
    assert model["calibrated_barrier_kcal_mol"] is None
    assert model["verdict"] == "no-correction"


def test_analysis_refuses_incomplete_source_classical_barrier() -> None:
    preparation = _prepared_analysis_fixture()
    preparation["models"]["model"]["classical_reference"]["status"] = (
        "incomplete-convergence"
    )
    result = e3b.analyze_models(preparation, _observation_fixture())
    model = result["models"]["model"]
    assert model["typed_outcome"] == "incomplete-source-classical"
    assert model["calibrated_barrier_kcal_mol"] is None


@pytest.mark.parametrize(
    ("mutation", "outcome"),
    [
        (("cp2k", "status", "incomplete-scf"), "incomplete-scf"),
        (("cp2k", "status", "incomplete-timeout"), "incomplete-timeout"),
        (("cp2k", "status", "incomplete-convergence"), "incomplete-convergence"),
        (
            ("matched_classical", "status", "incomplete-convergence"),
            "incomplete-matched-classical",
        ),
    ],
)
def test_analysis_refuses_incomplete_results(
    mutation: tuple[str, str, object], outcome: str
) -> None:
    observations = _observation_fixture()
    section, field, value = mutation
    observations["models"]["model"][section][field] = value
    result = e3b.analyze_models(
        _prepared_analysis_fixture(), observations, transfer_tolerance_kcal_mol=5.0
    )
    model = result["models"]["model"]
    assert model["typed_outcome"] == outcome
    assert model["matched_cell_correction_kcal_mol"] is None
    assert model["calibrated_barrier_kcal_mol"] is None


def test_analysis_does_not_trust_manual_transfer_numbers() -> None:
    observations = _observation_fixture()
    observations["models"]["model"]["matched_classical"]["barrier_kcal_mol"] = 106.0
    result = e3b.analyze_models(
        _prepared_analysis_fixture(), observations, transfer_tolerance_kcal_mol=5.0
    )
    model = result["models"]["model"]
    assert model["typed_outcome"] == "incomplete-unverified-observation"
    assert model["matched_cell_transfer_pass"] is False
    assert model["matched_cell_correction_kcal_mol"] is None


def test_analysis_rejects_atom_identity_drift() -> None:
    observations = _observation_fixture()
    observations["models"]["model"]["atom_identity_sha256"] = "different"
    result = e3b.analyze_models(_prepared_analysis_fixture(), observations)
    model = result["models"]["model"]
    assert model["typed_outcome"] == "incomplete-atom-identity"
    assert model["matched_cell_correction_kcal_mol"] is None


def test_analyze_command_generates_result_and_manifest(tmp_path: Path) -> None:
    prepared_root = tmp_path / "prepared"
    prepared_root.mkdir()
    preparation_path = prepared_root / "preparation.json"
    preparation = _prepared_analysis_fixture()
    preparation_path.write_text(
        json.dumps(preparation, sort_keys=True) + "\n", encoding="utf-8"
    )
    e3b.write_manifest(prepared_root)
    observations = _observation_fixture()
    observations["preparation_sha256"] = e3b.sha256(preparation_path)
    observations_path = tmp_path / "observations.json"
    observations_path.write_text(
        json.dumps(observations, sort_keys=True) + "\n", encoding="utf-8"
    )
    out = tmp_path / "analysis"
    args = e3b.argparse.Namespace(
        prepared_root=prepared_root,
        observations=observations_path,
        out=out,
        transfer_tolerance=5.0,
        endorsement_tolerance=5.0,
    )

    assert e3b.command_analyze(args) == 2
    result = json.loads((out / "analysis.json").read_text(encoding="utf-8"))
    assert (
        result["models"]["model"]["typed_outcome"]
        == "incomplete-unverified-observation"
    )
    e3b.verify_manifest(out)


@pytest.mark.skipif(
    not WORKBOOK.exists() or not CAMPAIGN.exists(),
    reason="hash-pinned E3a evidence unavailable",
)
def test_real_e3a_inputs_prepare_three_periodic_spot_checks(tmp_path: Path) -> None:
    pytest.importorskip("openpyxl")
    out = tmp_path / "prepared"
    preparation = e3b.prepare_campaign(
        workbook=WORKBOOK,
        campaign_root=CAMPAIGN,
        out=out,
        method=method_files(tmp_path),
    )

    assert preparation["schema"] == "e3b-cp2k-preparation-v1"
    assert tuple(preparation["models"]) == e3b.MODEL_ORDER
    assert preparation["source"]["workbook_sha256"] == e3b.WORKBOOK_SHA256
    assert preparation["source"]["campaign_result_sha256"] == (
        e3b.CAMPAIGN_RESULT_SHA256
    )

    expected_atoms = {
        "reconstructed-replication": 334,
        "dehydroxylate-lattice": 331,
        "xenon-divacancy": 334,
    }
    for name, record in preparation["models"].items():
        assert record["atom_count"] == expected_atoms[name]
        assert record["supercell"] == [2, 2, 1]
        assert record["gates"]["all_pass"] is True
        assert abs(record["gates"]["net_charge_e3a_e"]) <= 1.0e-8
        assert record["gates"]["route_pass"] is True
        assert record["gates"]["identity_pass"] is True
        model_dir = out / "models" / name
        neb_text = (model_dir / "neb.inp").read_text(encoding="utf-8")
        assert neb_text.count("&REPLICA") == 8
        assert "BAND_TYPE CI-NEB" in neb_text
        assert "TYPE DFTD3" in neb_text
        assert "PARAMETER_FILE_NAME ../../cp2k-data/dftd3.dat" in neb_text
        assert len(list((model_dir / "images").glob("replica-*.xyz"))) == 8
        atom_map = json.loads((model_dir / "atom-map.json").read_text())
        assert len(atom_map["atoms"]) == expected_atoms[name]
        assert len({atom["source_id"] for atom in atom_map["atoms"]}) == len(
            atom_map["atoms"]
        )
        assert (model_dir / "matched-classical" / "in.neb").is_file()

    dehyd = preparation["models"]["dehydroxylate-lattice"]
    assert dehyd["gates"]["five_coordinate_al_ids"] == [10, 11]
    assert dehyd["gates"]["coordination_pass"] is True
    assert dehyd["scientific_scope"] == (
        "local post-dehydroxylation Ar hop; not a dehydroxylation reaction path"
    )

    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    assert "manifest.json" not in {entry["path"] for entry in manifest["files"]}
    for entry in manifest["files"]:
        artifact = out / entry["path"]
        assert entry["bytes"] == artifact.stat().st_size
        assert entry["sha256"] == e3b.sha256(artifact)
