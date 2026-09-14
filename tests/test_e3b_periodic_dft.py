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
        "from pathlib import Path\n"
        "Path('runtime-side-effect.restart').write_text('runtime only')\n"
        "print('SCF run converged in 3 steps')\n"
        "print('ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: -42.25')\n"
        "print('PROGRAM ENDED AT 2026-09-14')\n",
        encoding="utf-8",
    )
    fake.chmod(0o755)
    prepared_root = tmp_path / "prepared"
    model_dir = prepared_root / "models" / "reconstructed-replication"
    data_dir = prepared_root / "cp2k-data"
    model_dir.mkdir(parents=True)
    data_dir.mkdir()
    input_path = model_dir / "smoke.inp"
    input_path.write_text(
        "&GLOBAL\n  RUN_TYPE ENERGY_FORCE\n&END GLOBAL\n"
        "BASIS_SET_FILE_NAME ../../cp2k-data/BASIS_SET\n"
        "COORD_FILE_NAME initial.xyz\n",
        encoding="utf-8",
    )
    (model_dir / "initial.xyz").write_text("0\nfixture\n", encoding="utf-8")
    (data_dir / "BASIS_SET").write_text("fixture basis\n", encoding="utf-8")
    (prepared_root / "preparation.json").write_text(
        json.dumps(
            {
                "schema": "e3b-cp2k-preparation-v1",
                "method": {
                    "cp2k_version": "2026.1",
                    "cp2k_image_digest": "sha256:" + "a" * 64,
                },
                "models": {
                    "reconstructed-replication": {"atom_identity_sha256": "identity"}
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )
    e3b.write_manifest(prepared_root)
    manifest_before = (prepared_root / "manifest.json").read_bytes()
    runtime_dir = tmp_path / "runtime"
    receipt_path = runtime_dir / "smoke-result.json"

    receipt = e3b.run_smoke(
        executable=str(fake),
        prepared_root=prepared_root,
        model="reconstructed-replication",
        runtime_dir=runtime_dir,
        receipt_path=receipt_path,
        timeout_seconds=5.0,
    )

    assert receipt["status"] == "converged"
    assert receipt["energy_hartree"] == -42.25
    assert receipt["elapsed_seconds"] >= 0.0
    assert receipt["schema"] == "e3b-cp2k-smoke-v2"
    assert receipt["input"]["sha256"] == e3b.sha256(input_path)
    assert receipt["output"]["stdout"]["sha256"] == e3b.sha256(
        runtime_dir / "smoke.stdout.log"
    )
    assert receipt["method"]["cp2k_version"] == "2026.1"
    assert receipt["resources"]["omp_threads"] == 1
    assert {item["prepared_path"] for item in receipt["dependencies"]} == {
        "cp2k-data/BASIS_SET",
        "models/reconstructed-replication/initial.xyz",
        "models/reconstructed-replication/smoke.inp",
    }
    assert receipt["cleanup"]["work_directory_removed"] is True
    assert not (runtime_dir / "work").exists()
    assert receipt["post_run_state"]["prepared_manifest_valid"] is True
    assert (prepared_root / "manifest.json").read_bytes() == manifest_before
    e3b.verify_manifest(prepared_root)
    assert json.loads(receipt_path.read_text(encoding="utf-8")) == receipt
    second_runtime = tmp_path / "runtime-second"
    second = e3b.run_smoke(
        executable=str(fake),
        prepared_root=prepared_root,
        model="reconstructed-replication",
        runtime_dir=second_runtime,
        receipt_path=second_runtime / "smoke-result.json",
        timeout_seconds=5.0,
    )
    assert second["status"] == "converged"
    e3b.verify_manifest(prepared_root)


def test_smoke_rejects_runtime_inside_prepared_tree(tmp_path: Path) -> None:
    prepared_root = tmp_path / "prepared"
    model_dir = prepared_root / "models" / "reconstructed-replication"
    model_dir.mkdir(parents=True)
    (model_dir / "smoke.inp").write_text("&GLOBAL\n&END GLOBAL\n")
    (prepared_root / "preparation.json").write_text(
        json.dumps(
            {
                "schema": "e3b-cp2k-preparation-v1",
                "method": {},
                "models": {
                    "reconstructed-replication": {"atom_identity_sha256": "identity"}
                },
            }
        )
        + "\n"
    )
    e3b.write_manifest(prepared_root)

    with pytest.raises(ValueError, match="outside the prepared tree"):
        e3b.run_smoke(
            executable="cp2k",
            prepared_root=prepared_root,
            model="reconstructed-replication",
            runtime_dir=prepared_root / "runtime",
            receipt_path=prepared_root / "runtime" / "smoke-result.json",
            timeout_seconds=5.0,
        )
    e3b.verify_manifest(prepared_root)


def _source_reference(
    barrier: float = 100.0, *, converged: bool = True
) -> dict[str, Any]:
    return {
        "schema": "e3a-source-barrier-reference-v1",
        "barrier_kcal_mol": barrier,
        "status": "converged" if converged else "incomplete-convergence",
        "source_cell": [6, 3, 1],
        "campaign_result": {
            "path": "campaign-result.json",
            "sha256": "1" * 64,
            "bytes": 123,
        },
        "campaign_model_record_sha256": "2" * 64,
        "raw_neb_output": {
            "path": "model/screen.log",
            "sha256": "3" * 64,
            "bytes": 456,
        },
        "raw_endpoint_outputs": {
            "initial": {
                "path": "model/initial.min.stdout.log",
                "sha256": "4" * 64,
                "bytes": 100,
            },
            "endpoint": {
                "path": "model/endpoint.min.stdout.log",
                "sha256": "5" * 64,
                "bytes": 101,
            },
        },
        "convergence": {
            "independent_parser": "e3a_campaign.parse_last_thermo+parse_neb",
            "requested_ftol_kcal_mol_angstrom": 0.01,
            "max_replica_force_kcal_mol_angstrom": 0.005,
            "converged_to_requested_ftol": converged,
            "final_step": 1000,
            "expected_final_step": 1000,
            "complete_step_history": True,
            "parsed_forward_barrier_kcal_mol": barrier,
            "initial_endpoint_force_pass": True,
            "final_endpoint_force_pass": True,
            "initial_endpoint_normal_termination": True,
            "final_endpoint_normal_termination": True,
        },
    }


def _prepared_analysis_fixture() -> dict[str, Any]:
    return {
        "schema": "e3b-cp2k-preparation-v1",
        "models": {
            "model": {
                "atom_identity_sha256": "identity",
                "classical_reference": _source_reference(),
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


def _write_raw_fixture(root: Path, *, replica_count: int = 8) -> dict[str, Path]:
    root.mkdir(parents=True)
    paths: dict[str, Path] = {}
    for role, energy in (("cp2k_initial", -100.0), ("cp2k_endpoint", -100.01)):
        path = root / f"{role}.log"
        path.write_text(
            "SCF run converged in 3 steps\n"
            f"ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: {energy}\n"
            "GEOMETRY OPTIMIZATION COMPLETED\n"
            "PROGRAM ENDED AT 2026-09-14\n",
            encoding="utf-8",
        )
        paths[role] = path
    band = root / "band.log"
    band.write_text(
        "".join("SCF run converged in 4 steps\n" for _ in range(replica_count))
        + f"BAND| Number of replicas {replica_count}\n"
        + "BAND OPTIMIZATION CONVERGED\nPROGRAM ENDED AT 2026-09-14\n",
        encoding="utf-8",
    )
    paths["cp2k_band"] = band
    energies = root / "replica-energies.log"
    barrier_hartree = 112.0 / e3b.HARTREE_TO_KCAL_MOL
    profile = [
        -100.0,
        -100.0 + barrier_hartree * 0.3,
        -100.0 + barrier_hartree * 0.7,
        -100.0 + barrier_hartree,
        -100.0 + barrier_hartree * 0.8,
        -100.0 + barrier_hartree * 0.5,
        -100.0 + barrier_hartree * 0.2,
        -100.01,
    ][:replica_count]
    energies.write_text(
        "".join(
            f"Replica {index + 1} energy [a.u.] = {energy:.16f}\n"
            for index, energy in enumerate(profile)
        ),
        encoding="utf-8",
    )
    paths["cp2k_replica_energies"] = energies
    for role, energy in (("lammps_initial", -500.0), ("lammps_endpoint", -501.0)):
        path = root / f"{role}.log"
        path.write_text(
            f"100 {energy} 0.005 0.001 0.0\nTotal wall time: 0:00:01\n",
            encoding="utf-8",
        )
        paths[role] = path
    neb = root / "neb.screen"
    row = [1000.0, 0.005, 0.001, 0.0, 0.0, 0.0, 102.0, 103.0, 0.0]
    row.extend(float(index) for index in range(2 * replica_count))
    neb.write_text(
        " ".join(str(value) for value in row) + "\nTotal wall time: 0:00:01\n",
        encoding="utf-8",
    )
    paths["lammps_neb"] = neb
    return paths


def _write_prepared_fixture(tmp_path: Path) -> tuple[Path, dict[str, Any]]:
    prepared_root = tmp_path / "prepared-evidence"
    prepared_root.mkdir()
    preparation = _prepared_analysis_fixture()
    preparation_path = prepared_root / "preparation.json"
    e3b.write_json(preparation_path, preparation)
    e3b.write_manifest(prepared_root)
    return prepared_root, preparation


def test_hash_bound_raw_fixtures_produce_calibration(tmp_path: Path) -> None:
    prepared_root, preparation = _write_prepared_fixture(tmp_path)
    raw = _write_raw_fixture(tmp_path / "raw")
    evidence_root = tmp_path / "evidence"
    evidence = e3b.collect_evidence(
        prepared_root=prepared_root,
        raw_outputs={"model": raw},
        out=evidence_root,
    )

    result = e3b.analyze_models(
        preparation,
        evidence,
        evidence_root=evidence_root,
        preparation_sha256=e3b.sha256(prepared_root / "preparation.json"),
    )
    model = result["models"]["model"]
    assert model["typed_outcome"] == "complete-endorsed"
    assert model["matched_cell_transfer_pass"] is True
    assert model["matched_cell_transfer_delta_kcal_mol"] == pytest.approx(2.0)
    assert model["matched_cell_correction_kcal_mol"] == pytest.approx(10.0)
    assert model["calibrated_barrier_kcal_mol"] == pytest.approx(110.0)
    assert model["endorsement_delta_kcal_mol"] == pytest.approx(2.0)
    assert model["verdict"] == "endorsed"
    e3b.verify_manifest(evidence_root)


def test_collect_and_analyze_cli_handlers_use_raw_spec(tmp_path: Path) -> None:
    prepared_root, _ = _write_prepared_fixture(tmp_path)
    raw = _write_raw_fixture(tmp_path / "raw-with-arbitrary-names")
    raw_spec = tmp_path / "raw-spec.json"
    e3b.write_json(
        raw_spec,
        {
            "schema": "e3b-raw-output-spec-v1",
            "models": {
                "model": {role: str(path) for role, path in raw.items()},
            },
        },
    )
    evidence_root = tmp_path / "evidence-cli"
    collect_args = e3b.build_parser().parse_args(
        [
            "collect-evidence",
            "--prepared-root",
            str(prepared_root),
            "--raw-spec",
            str(raw_spec),
            "--out",
            str(evidence_root),
        ]
    )
    assert collect_args.func(collect_args) == 0

    analysis_root = tmp_path / "analysis-cli"
    analyze_args = e3b.build_parser().parse_args(
        [
            "analyze",
            "--prepared-root",
            str(prepared_root),
            "--evidence-root",
            str(evidence_root),
            "--out",
            str(analysis_root),
        ]
    )
    assert analyze_args.func(analyze_args) == 0
    result = e3b.read_json(analysis_root / "analysis.json")
    assert result["models"]["model"]["typed_outcome"] == "complete-endorsed"
    e3b.verify_manifest(analysis_root)


def test_analysis_reparses_raw_files_and_rejects_tampered_evidence(
    tmp_path: Path,
) -> None:
    prepared_root, preparation = _write_prepared_fixture(tmp_path)
    evidence_root = tmp_path / "evidence"
    evidence = e3b.collect_evidence(
        prepared_root=prepared_root,
        raw_outputs={"model": _write_raw_fixture(tmp_path / "raw")},
        out=evidence_root,
    )
    evidence["models"]["model"]["cp2k"]["barrier_kcal_mol"] = 1.0

    result = e3b.analyze_models(
        preparation,
        evidence,
        evidence_root=evidence_root,
        preparation_sha256=e3b.sha256(prepared_root / "preparation.json"),
    )
    assert result["models"]["model"]["typed_outcome"] == "incomplete-evidence-integrity"
    assert result["models"]["model"]["calibrated_barrier_kcal_mol"] is None


def test_analysis_rejects_raw_output_changed_after_collection(tmp_path: Path) -> None:
    prepared_root, preparation = _write_prepared_fixture(tmp_path)
    evidence_root = tmp_path / "evidence"
    evidence = e3b.collect_evidence(
        prepared_root=prepared_root,
        raw_outputs={"model": _write_raw_fixture(tmp_path / "raw")},
        out=evidence_root,
    )
    raw_path = (
        evidence_root / evidence["models"]["model"]["raw_outputs"]["cp2k_band"]["path"]
    )
    raw_path.write_text(raw_path.read_text() + "tampered\n")

    result = e3b.analyze_models(
        preparation,
        evidence,
        evidence_root=evidence_root,
        preparation_sha256=e3b.sha256(prepared_root / "preparation.json"),
    )
    model = result["models"]["model"]
    assert model["typed_outcome"] == "incomplete-evidence-integrity"
    assert model["matched_cell_correction_kcal_mol"] is None


def test_analysis_rejects_nominally_converged_unbound_source(tmp_path: Path) -> None:
    prepared_root, preparation = _write_prepared_fixture(tmp_path)
    del preparation["models"]["model"]["classical_reference"]["raw_endpoint_outputs"]
    evidence_root = tmp_path / "evidence"
    evidence = e3b.collect_evidence(
        prepared_root=prepared_root,
        raw_outputs={"model": _write_raw_fixture(tmp_path / "raw")},
        out=evidence_root,
    )

    result = e3b.analyze_models(
        preparation,
        evidence,
        evidence_root=evidence_root,
        preparation_sha256=e3b.sha256(prepared_root / "preparation.json"),
    )
    model = result["models"]["model"]
    assert model["typed_outcome"] == "incomplete-evidence-integrity"
    assert model["calibrated_barrier_kcal_mol"] is None


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("missing-normal-end", "normal termination"),
        ("missing-scf", "SCF"),
        ("missing-band-convergence", "BAND convergence"),
        ("wrong-replica-count", "replica count"),
        ("nonfinite-energy", "finite"),
    ],
)
def test_cp2k_raw_parser_fails_closed(
    tmp_path: Path, mutation: str, message: str
) -> None:
    raw = _write_raw_fixture(
        tmp_path / "raw", replica_count=7 if mutation == "wrong-replica-count" else 8
    )
    if mutation == "missing-normal-end":
        raw["cp2k_initial"].write_text(
            raw["cp2k_initial"].read_text().replace("PROGRAM ENDED AT", "STOPPED AT")
        )
    elif mutation == "missing-scf":
        raw["cp2k_initial"].write_text(
            raw["cp2k_initial"].read_text().replace("SCF run converged", "SCF run")
        )
    elif mutation == "missing-band-convergence":
        raw["cp2k_band"].write_text(
            raw["cp2k_band"]
            .read_text()
            .replace("BAND OPTIMIZATION CONVERGED", "BAND OPTIMIZATION STOPPED")
        )
    elif mutation == "nonfinite-energy":
        lines = raw["cp2k_replica_energies"].read_text().splitlines()
        lines[3] = "Replica 4 energy [a.u.] = NaN"
        raw["cp2k_replica_energies"].write_text("\n".join(lines) + "\n")

    with pytest.raises(ValueError, match=message):
        e3b.parse_cp2k_barrier(
            initial_log=raw["cp2k_initial"],
            endpoint_log=raw["cp2k_endpoint"],
            band_log=raw["cp2k_band"],
            replica_energies=raw["cp2k_replica_energies"],
            expected_replicas=8,
        )


def test_cp2k_parser_accepts_replica_block_energy_output(tmp_path: Path) -> None:
    raw = _write_raw_fixture(tmp_path / "raw")
    blocks = []
    for line in raw["cp2k_replica_energies"].read_text().splitlines():
        match = e3b.re.match(r"Replica (\d+) energy \[a.u.\] = (\S+)", line)
        assert match is not None
        blocks.append(
            f"----- REPLICA Nr. {match.group(1)} -----\n"
            "ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: "
            f"{match.group(2)}\n"
        )
    raw["cp2k_replica_energies"].write_text("".join(blocks))

    parsed = e3b.parse_cp2k_barrier(
        initial_log=raw["cp2k_initial"],
        endpoint_log=raw["cp2k_endpoint"],
        band_log=raw["cp2k_band"],
        replica_energies=raw["cp2k_replica_energies"],
    )
    assert parsed["barrier_kcal_mol"] == pytest.approx(112.0)


@pytest.mark.parametrize(
    "mutation",
    ["missing-normal-end", "unconverged-endpoint", "wrong-replica-count"],
)
def test_matched_classical_raw_parser_fails_closed(
    tmp_path: Path, mutation: str
) -> None:
    raw = _write_raw_fixture(tmp_path / "raw")
    if mutation == "missing-normal-end":
        raw["lammps_initial"].write_text(
            raw["lammps_initial"].read_text().replace("Total wall time", "Stopped")
        )
    elif mutation == "unconverged-endpoint":
        raw["lammps_initial"].write_text(
            raw["lammps_initial"].read_text().replace("0.005", "0.05")
        )
    else:
        fields = raw["lammps_neb"].read_text().splitlines()[0].split()
        raw["lammps_neb"].write_text(
            " ".join(fields[:-2]) + "\nTotal wall time: 0:00:01\n"
        )

    with pytest.raises(ValueError):
        e3b.parse_matched_classical_barrier(
            initial_log=raw["lammps_initial"],
            endpoint_log=raw["lammps_endpoint"],
            neb_log=raw["lammps_neb"],
        )


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
    assert tuple(preparation["source"]["reduction"]["window_start"]) == (5, 2)

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
    gate = dehyd["gates"]
    assert gate["five_coordinate_al_ids"] == [10, 11]
    assert gate["coordination_pass"] is True
    criterion = gate["five_coordinate_al_criterion"]
    assert criterion["schema"] == "e3b-five-coordinate-al-criterion-v1"
    assert criterion["parent_model"] == "reconstructed-replication"
    assert criterion["al_o_cutoff_angstrom"] == 2.3
    assert criterion["route_interaction_cutoff_angstrom"] == 6.0
    assert criterion["removed_hydroxyl_oxygen_ids"] == [76]
    assert criterion["candidate_al_source_ids"] == [10, 11]
    for candidate in criterion["candidates"]:
        assert candidate["parent_o_coordination"] == 6
        assert candidate["transformed_o_coordination"] == 5
        assert candidate["removed_hydroxyl_o_in_parent_shell"] is True
        assert candidate["removed_hydroxyl_o_declared_removed"] is True
        assert candidate["removed_hydroxyl_o_absent_from_transformed_model"] is True
        assert candidate["route_interaction_pass"] is True
        assert len(candidate["route_image_distances_angstrom"]) == 8
        assert candidate["minimum_route_distance_angstrom"] <= 6.0
    assert dehyd["scientific_scope"] == (
        "local post-dehydroxylation Ar hop; not a dehydroxylation reaction path"
    )

    for record in preparation["models"].values():
        source = record["classical_reference"]
        assert source["schema"] == "e3a-source-barrier-reference-v1"
        assert source["status"] == "incomplete-convergence"
        assert source["campaign_result"]["sha256"] == e3b.CAMPAIGN_RESULT_SHA256
        assert len(source["campaign_model_record_sha256"]) == 64
        assert len(source["raw_neb_output"]["sha256"]) == 64
        assert set(source["raw_endpoint_outputs"]) == {"initial", "endpoint"}
        assert all(
            len(item["sha256"]) == 64
            for item in source["raw_endpoint_outputs"].values()
        )
        assert source["convergence"]["independent_parser"] == (
            "e3a_campaign.parse_last_thermo+parse_neb"
        )
        assert source["convergence"]["converged_to_requested_ftol"] is False
        assert (
            source["convergence"]["max_replica_force_kcal_mol_angstrom"]
            > (source["convergence"]["requested_ftol_kcal_mol_angstrom"])
        )

    source_atoms, source_cell, barriers = e3b.contract.nteme_neb.read_workbook(WORKBOOK)
    route = e3b.contract.nteme_neb.select_route(barriers, "divacancy", 1)
    pristine, target_cell, _ = e3b.reduce_pristine(source_atoms, source_cell)
    built = {
        model.name: model
        for model in e3b.build_spot_models(pristine, target_cell, route)
    }
    parent_model = built["reconstructed-replication"]
    dehyd_model = built["dehydroxylate-lattice"]
    bad_transformation = dict(dehyd_model.metadata["transformation"])
    bad_transformation["removed_site_ids"] = [81, 84]
    bad_model = e3b.replace(
        dehyd_model,
        metadata={**dehyd_model.metadata, "transformation": bad_transformation},
    )
    provenance_gate = e3b._coordination_gate(bad_model, parent_model)
    assert provenance_gate["five_coordinate_al_ids"] == [10, 11]
    assert provenance_gate["coordination_pass"] is False
    assert all(
        candidate["route_interaction_pass"] is True
        and candidate["removed_hydroxyl_o_declared_removed"] is False
        for candidate in provenance_gate["five_coordinate_al_criterion"]["candidates"]
    )

    manifest = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    assert "manifest.json" not in {entry["path"] for entry in manifest["files"]}
    for entry in manifest["files"]:
        artifact = out / entry["path"]
        assert entry["bytes"] == artifact.stat().st_size
        assert entry["sha256"] == e3b.sha256(artifact)
