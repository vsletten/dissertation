from __future__ import annotations

import importlib.util
import sys
import time
from dataclasses import replace
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

MODEL = "model"
CELL = {"a": 10.0, "b": 11.0, "c": 12.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0}


def _write_e3a_source(
    root: Path, barrier: float = 100.0
) -> tuple[dict[str, Any], dict[str, Any]]:
    workbook = root / "source-workbook.fixture"
    workbook.write_text("deterministic E3a workbook fixture\n", encoding="utf-8")
    campaign_root = root / "e3a-campaign"
    model_root = campaign_root / MODEL
    model_root.mkdir(parents=True)

    endpoint_paths = {}
    for label, energy in (("initial", -700.0), ("endpoint", -701.0)):
        path = model_root / f"{label}.min.stdout.log"
        path.write_text(
            f"100 {energy} 0.005 0.001 0.0\nTotal wall time: 0:00:01\n",
            encoding="utf-8",
        )
        endpoint_paths[label] = path

    neb_path = model_root / "screen.log"
    row = [1000.0, 0.005, 0.001, 0.0, 0.0, 0.0, barrier, barrier + 1.0, 0.0]
    row.extend(float(index) for index in range(16))
    neb_path.write_text(
        " ".join(str(value) for value in row) + "\nTotal wall time: 0:00:01\n",
        encoding="utf-8",
    )
    model_record = {
        "neb": {
            "computed_forward_barrier_kcal_mol": barrier,
            "max_replica_force_kcal_mol_angstrom": 0.005,
            "converged_to_requested_ftol": True,
            "screen_sha256": e3b.sha256(neb_path),
        }
    }
    campaign = {
        "schema": "e3a-classical-neb-campaign-v1",
        "source": {"source_sha256": e3b.sha256(workbook)},
        "settings": {
            "neb_ftol_kcal_mol_angstrom": 0.01,
            "min_ftol_kcal_mol_angstrom": 0.01,
            "neb_relax_steps": 400,
            "neb_climb_steps": 600,
        },
        "models": {MODEL: model_record},
    }
    result_path = campaign_root / "campaign-result.json"
    e3b.write_json(result_path, campaign)
    manifest_entries = []
    for path in (result_path, neb_path, *endpoint_paths.values()):
        manifest_entries.append(
            {
                "path": str(path.relative_to(campaign_root)),
                "bytes": path.stat().st_size,
                "sha256": e3b.sha256(path),
            }
        )
    manifest_path = campaign_root / "manifest.json"
    e3b.write_json(
        manifest_path,
        {"schema": "sha256-manifest-v1", "files": manifest_entries},
    )
    source = {
        "workbook_path": str(workbook.resolve()),
        "workbook_sha256": e3b.sha256(workbook),
        "campaign_root": str(campaign_root.resolve()),
        "campaign_manifest_sha256": e3b.sha256(manifest_path),
        "campaign_result_sha256": e3b.sha256(result_path),
    }
    reference = e3b._classical_reference(
        campaign, MODEL, campaign_root, e3b.read_json(manifest_path)
    )
    assert reference["status"] == "converged"
    return source, reference


def _write_prepared_fixture(
    root: Path, *, source_barrier: float = 100.0
) -> tuple[Path, dict[str, Any]]:
    prepared = root / "prepared"
    model_root = prepared / "models" / MODEL
    matched = model_root / "matched-classical"
    images = model_root / "images"
    data = prepared / "cp2k-data"
    matched.mkdir(parents=True)
    images.mkdir()
    data.mkdir(parents=True)

    for name in ("BASIS_SET", "GTH_POTENTIALS", "dftd3.dat"):
        (data / name).write_text(f"fixture {name}\n", encoding="utf-8")
    for name in ("initial.xyz", "endpoint.xyz"):
        (model_root / name).write_text("3\nfixture\n", encoding="utf-8")
    for name in ("initial-opt.inp", "endpoint-opt.inp", "neb.inp", "smoke.inp"):
        (model_root / name).write_text(
            f"&GLOBAL\n PROJECT {name}\n&END GLOBAL\n", encoding="utf-8"
        )
    for index in range(e3b.IMAGE_COUNT):
        (images / f"replica-{index:02d}.xyz").write_text(
            "3\nfixture image\n", encoding="utf-8"
        )
    for name in (
        "in.min.initial",
        "initial.seed.data",
        "in.min.endpoint",
        "endpoint.seed.data",
        "in.neb",
    ):
        (matched / name).write_text(f"fixture {name}\n", encoding="utf-8")

    source, reference = _write_e3a_source(root, source_barrier)
    method = {
        "cp2k_version": "2024.3",
        "cp2k_image_digest": "sha256:" + "a" * 64,
        "files": {
            name: {
                "artifact_path": f"cp2k-data/{name}",
                "sha256": e3b.sha256(data / name),
            }
            for name in ("BASIS_SET", "GTH_POTENTIALS", "dftd3.dat")
        },
    }
    preparation = {
        "schema": "e3b-cp2k-preparation-v1",
        "source": source,
        "method": method,
        "models": {
            MODEL: {
                "atom_count": 3,
                "atom_identity_sha256": "identity-model",
                "periodic_cell": CELL,
                "periodic_cell_sha256": e3b.canonical_sha256(CELL),
                "classical_reference": reference,
            }
        },
    }
    e3b.write_json(prepared / "preparation.json", preparation)
    e3b.write_manifest(prepared)
    return prepared, preparation


def _identity(project: str, atom_count: int = 3, *, a: float = 10.0) -> str:
    return (
        f" GLOBAL| Project name {project}\n"
        f" - Atoms: {atom_count}\n"
        f" CELL| Vector a [angstrom]: {a:.8f} 0.00000000 0.00000000\n"
        " CELL| Vector b [angstrom]: 0.00000000 11.00000000 0.00000000\n"
        " CELL| Vector c [angstrom]: 0.00000000 0.00000000 12.00000000\n"
    )


def _cp2k_profile(barrier: float) -> list[float]:
    delta = barrier / e3b.HARTREE_TO_KCAL_MOL
    return [
        -100.0,
        -100.0 + 0.3 * delta,
        -100.0 + 0.7 * delta,
        -100.0 + delta,
        -100.0 + 0.8 * delta,
        -100.0 + 0.5 * delta,
        -100.0 + 0.2 * delta,
        -100.01,
    ]


def _lammps_neb_row(barrier: float, max_force: float) -> list[float]:
    row = [1000.0, max_force, 0.001, 0.0, 0.0, 0.0, barrier, barrier + 1.0, 0.0]
    energies = [-500.0, -490.0, -480.0, -470.0, -475.0, -485.0, -495.0, -501.0]
    for index, energy in enumerate(energies):
        row.extend((index / 7.0, energy))
    return row


def _write_raw_fixture(
    root: Path,
    *,
    model: str = MODEL,
    cp2k_barrier: float = 112.0,
    matched_barrier: float = 102.0,
) -> dict[str, Path]:
    root.mkdir(parents=True)
    paths: dict[str, Path] = {}
    for role, suffix, energy in (
        ("cp2k_initial", "initial", -100.0),
        ("cp2k_endpoint", "endpoint", -100.01),
    ):
        path = root / f"{role}.out"
        path.write_text(
            _identity(f"e3b-{model}-{suffix}")
            + " SCF run converged in 3 steps\n"
            + f" ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: {energy}\n"
            + " GEOMETRY OPTIMIZATION COMPLETED\n"
            + " PROGRAM ENDED AT 2026-09-14 12:00:00\n",
            encoding="utf-8",
        )
        paths[role] = path

    profile = _cp2k_profile(cp2k_barrier)
    band = root / "e3b-band.out"
    band.write_text(
        _identity(f"e3b-{model}-neb")
        + " BAND TYPE = CI-NEB\n"
        + " BAND TYPE OPTIMIZATION = DIIS\n"
        + " STEP NUMBER = 44\n"
        + f" NUMBER OF NEB REPLICA = {e3b.IMAGE_COUNT}\n"
        + " ENERGIES [au] = "
        + " ".join(f"{value:.16f}" for value in profile)
        + "\n BAND TOTAL ENERGY [au] = -799.0\n"
        + " RMS DISPLACEMENT = 0.00010 (YES)\n"
        + " MAX DISPLACEMENT = 0.00020 (YES)\n"
        + " RMS FORCE = 0.00010 (YES)\n"
        + " MAX FORCE = 0.00020 (YES)\n"
        + " PROGRAM ENDED AT 2026-09-14 12:00:00\n",
        encoding="utf-8",
    )
    paths["cp2k_band"] = band
    profile_path = root / "e3b-band-final.profile"
    profile_path.write_text(
        "ENERGIES [au] = " + " ".join(f"{value:.16f}" for value in profile) + "\n",
        encoding="utf-8",
    )
    paths["cp2k_band_profile"] = profile_path

    for index, role in enumerate(e3b.CP2K_REPLICA_OUTPUT_ROLES, 1):
        path = root / f"e3b-{model}-neb-BAND{index:02d}.out"
        path.write_text(
            _identity(f"e3b-{model}-neb-BAND{index:02d}")
            + " SCF run converged in 4 steps\n"
            + " ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: -100.0\n"
            + " PROGRAM ENDED AT 2026-09-14 12:00:00\n",
            encoding="utf-8",
        )
        paths[role] = path

    for role, energy in (("lammps_initial", -500.0), ("lammps_endpoint", -501.0)):
        path = root / f"{role}.log"
        path.write_text(
            f"100 {energy} 0.005 0.001 0.0\nTotal wall time: 0:00:01\n",
            encoding="utf-8",
        )
        paths[role] = path

    header = (
        "Step MaxReplicaForce MaxAtomForce GradV0 GradV1 GradVc EBF EBR RDT "
        + " ".join(f"RD{index} PE{index}" for index in range(e3b.IMAGE_COUNT))
    )
    preclimb = _lammps_neb_row(90.0, 0.2)
    climb = _lammps_neb_row(matched_barrier, 0.005)
    neb = root / "neb.screen"
    neb.write_text(
        header
        + "\n"
        + " ".join(str(value) for value in preclimb)
        + "\nSetting up climbing ...\nClimbing replica = 4\n"
        + header
        + "\n"
        + " ".join(str(value) for value in climb)
        + "\nLoop time of 1.0 on 8 procs for 1000 steps\nTotal wall time: 0:00:01\n",
        encoding="utf-8",
    )
    paths["lammps_neb"] = neb
    assert set(paths) == set(e3b.RAW_OUTPUT_ROLES)
    return paths


def _fake_executable(path: Path, body: str = "print('fixture executable')\n") -> Path:
    path.write_text("#!/usr/bin/env python3\n" + body, encoding="utf-8")
    path.chmod(0o755)
    return path


def _create_receipt(
    prepared: Path, raw: dict[str, Path], path: Path, executable: Path
) -> dict[str, Any]:
    return e3b.create_execution_receipt(
        prepared_root=prepared,
        model=MODEL,
        raw_outputs=raw,
        executables={"cp2k": executable, "lammps": executable},
        executable_identities={"cp2k": "test/fake", "lammps": "test/fake"},
        commands={unit: [str(executable), unit] for unit in e3b.SOLVER_UNITS},
        resources={unit: {"ranks": 1, "threads": 1} for unit in e3b.SOLVER_UNITS},
        out=path,
    )


def _collect_bundle(
    tmp_path: Path,
    *,
    source_barrier: float = 100.0,
    cp2k_barrier: float = 112.0,
    matched_barrier: float = 102.0,
) -> tuple[Path, dict[str, Any], Path, dict[str, Any]]:
    prepared, preparation = _write_prepared_fixture(
        tmp_path, source_barrier=source_barrier
    )
    raw = _write_raw_fixture(
        tmp_path / "raw",
        cp2k_barrier=cp2k_barrier,
        matched_barrier=matched_barrier,
    )
    fake = _fake_executable(tmp_path / "solver-fixture")
    receipt_path = tmp_path / "execution-receipt.json"
    _create_receipt(prepared, raw, receipt_path, fake)
    evidence_root = tmp_path / "evidence"
    evidence = e3b.collect_evidence(
        prepared_root=prepared,
        raw_outputs={MODEL: raw},
        execution_receipts={MODEL: receipt_path},
        out=evidence_root,
    )
    return prepared, preparation, evidence_root, evidence


def _parse_cp2k(raw: dict[str, Path], **overrides: Any) -> dict[str, Any]:
    arguments = {
        "initial_log": raw["cp2k_initial"],
        "endpoint_log": raw["cp2k_endpoint"],
        "band_log": raw["cp2k_band"],
        "profile_log": raw["cp2k_band_profile"],
        "replica_logs": [raw[role] for role in e3b.CP2K_REPLICA_OUTPUT_ROLES],
        "expected_model": MODEL,
        "expected_atom_count": 3,
        "expected_cell": CELL,
    }
    arguments.update(overrides)
    return e3b.parse_cp2k_barrier(**arguments)


def _parse_lammps(raw: dict[str, Path]) -> dict[str, Any]:
    return e3b.parse_matched_classical_barrier(
        initial_log=raw["lammps_initial"],
        endpoint_log=raw["lammps_endpoint"],
        neb_log=raw["lammps_neb"],
    )


def test_cp2k_2024_3_native_band_and_replica_outputs_parse(tmp_path: Path) -> None:
    raw = _write_raw_fixture(tmp_path / "raw")
    parsed = _parse_cp2k(raw)
    assert parsed["barrier_kcal_mol"] == pytest.approx(112.0)
    assert parsed["final_convergence_criteria"] == {
        "rms_displacement": True,
        "max_displacement": True,
        "rms_force": True,
        "max_force": True,
    }
    assert len(parsed["replica_outputs"]) == 8
    assert all(item["scf_converged"] for item in parsed["replica_outputs"])
    assert all(item["normal_termination"] for item in parsed["replica_outputs"])


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("wrong-count", "replica count"),
        ("incomplete-criteria", "criteria incomplete"),
        ("criterion-no", "did not all pass"),
        ("missing-replica-log", "replica log count"),
        ("replica-scf-failure", "SCF did not converge"),
        ("replica-no-footer", "normal termination"),
        ("nonfinite-profile", "not finite"),
        ("malformed-profile", "wrong replica count"),
    ],
)
def test_cp2k_native_parser_fails_closed(
    tmp_path: Path, mutation: str, message: str
) -> None:
    raw = _write_raw_fixture(tmp_path / "raw")
    replica_logs = [raw[role] for role in e3b.CP2K_REPLICA_OUTPUT_ROLES]
    if mutation == "wrong-count":
        raw["cp2k_band"].write_text(
            raw["cp2k_band"]
            .read_text()
            .replace("NUMBER OF NEB REPLICA = 8", "NUMBER OF NEB REPLICA = 7")
        )
    elif mutation == "incomplete-criteria":
        raw["cp2k_band"].write_text(
            raw["cp2k_band"].read_text().replace(" RMS FORCE = 0.00010 (YES)\n", "")
        )
    elif mutation == "criterion-no":
        raw["cp2k_band"].write_text(
            raw["cp2k_band"]
            .read_text()
            .replace("MAX FORCE = 0.00020 (YES)", "MAX FORCE = 0.00020 (NO)")
        )
    elif mutation == "missing-replica-log":
        replica_logs.pop()
    elif mutation == "replica-scf-failure":
        replica_logs[3].write_text(
            replica_logs[3]
            .read_text()
            .replace("SCF run converged in 4 steps", "SCF run NOT converged")
        )
    elif mutation == "replica-no-footer":
        replica_logs[3].write_text(
            replica_logs[3].read_text().replace("PROGRAM ENDED AT", "STOPPED AT")
        )
    elif mutation == "nonfinite-profile":
        fields = raw["cp2k_band_profile"].read_text().split()
        fields[6] = "NaN"
        raw["cp2k_band_profile"].write_text(" ".join(fields) + "\n")
    else:
        raw["cp2k_band_profile"].write_text("ENERGIES [au] = -100 -99\n")

    with pytest.raises(ValueError, match=message):
        _parse_cp2k(raw, replica_logs=replica_logs)


def test_cp2k_parser_accepts_numeric_ener_rows(tmp_path: Path) -> None:
    raw = _write_raw_fixture(tmp_path / "raw")
    profile = _cp2k_profile(112.0)
    raw["cp2k_band_profile"].write_text(
        "# Step Nr. Time[fs] Kin.[a.u.] Temp[K] Pot.[a.u.]\n"
        + "".join(
            f"{index} 0.0 0.0 0.0 {energy:.16f} 0.0\n"
            for index, energy in enumerate(profile)
        )
    )
    assert _parse_cp2k(raw)["barrier_kcal_mol"] == pytest.approx(112.0)


def test_cp2k_native_identity_is_cross_checked(tmp_path: Path) -> None:
    raw = _write_raw_fixture(tmp_path / "raw", model="other")
    with pytest.raises(ValueError, match="project"):
        _parse_cp2k(raw)


def test_lammps_parser_requires_completed_climbing_phase(tmp_path: Path) -> None:
    raw = _write_raw_fixture(tmp_path / "raw")
    parsed = _parse_lammps(raw)
    assert parsed["barrier_kcal_mol"] == 102.0
    assert parsed["climbing_phase"] is True
    assert parsed["endpoint_energy_deltas_kcal_mol"] == {
        "initial": 0.0,
        "endpoint": 0.0,
    }


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("preclimb-only", "climbing phase"),
        ("missing-second-header", "climbing progress header"),
        ("aborted", "normal termination after climb"),
        ("fatal", "fatal/error"),
        ("climb-unconverged", "did not converge"),
        ("endpoint-mismatch", "endpoint energies"),
    ],
)
def test_lammps_regular_aborted_or_mismatched_neb_cannot_promote(
    tmp_path: Path, mutation: str, message: str
) -> None:
    raw = _write_raw_fixture(tmp_path / "raw")
    text = raw["lammps_neb"].read_text()
    if mutation == "preclimb-only":
        text = text.split("Setting up climbing")[0] + "Total wall time: 0:00:01\n"
    elif mutation == "missing-second-header":
        before, after = text.split("Climbing replica = 4\n", 1)
        after = after.split("\n", 1)[1]
        text = before + "Climbing replica = 4\n" + after
    elif mutation == "aborted":
        text = text.replace("Loop time of 1.0 on 8 procs for 1000 steps\n", "")
        text = text.replace("Total wall time: 0:00:01\n", "")
    elif mutation == "fatal":
        text += "ERROR: Lost atoms\n"
    elif mutation == "climb-unconverged":
        lines = text.splitlines()
        fields = lines[-3].split()
        fields[1] = "0.05"
        lines[-3] = " ".join(fields)
        text = "\n".join(lines) + "\n"
    else:
        lines = text.splitlines()
        fields = lines[-3].split()
        fields[10] = "-499.0"
        lines[-3] = " ".join(fields)
        text = "\n".join(lines) + "\n"
    raw["lammps_neb"].write_text(text)
    with pytest.raises(ValueError, match=message):
        _parse_lammps(raw)


def test_execution_receipt_binds_every_prepared_input_and_raw_output(
    tmp_path: Path,
) -> None:
    prepared, _ = _write_prepared_fixture(tmp_path)
    raw = _write_raw_fixture(tmp_path / "raw")
    fake = _fake_executable(tmp_path / "solver")
    receipt_path = tmp_path / "receipt.json"
    receipt = _create_receipt(prepared, raw, receipt_path, fake)
    assert receipt["preparation"]["sha256"] == e3b.sha256(prepared / "preparation.json")
    assert receipt["periodic_cell_record_sha256"] == e3b.canonical_sha256(CELL)
    assert set(receipt["raw_outputs"]) == set(e3b.RAW_OUTPUT_ROLES)
    assert set(receipt["units"]) == set(e3b.SOLVER_UNITS)
    assert all(unit["prepared_inputs"] for unit in receipt["units"].values())
    assert receipt["executables"]["cp2k"]["identity"] == "test/fake"
    assert receipt["executables"]["cp2k"]["identity_verification"] == "test-fixture"


def test_collect_requires_receipt_and_rejects_changed_output(tmp_path: Path) -> None:
    prepared, _ = _write_prepared_fixture(tmp_path)
    raw = _write_raw_fixture(tmp_path / "raw")
    with pytest.raises(ValueError, match="execution-receipt model set"):
        e3b.collect_evidence(
            prepared_root=prepared,
            raw_outputs={MODEL: raw},
            execution_receipts={},
            out=tmp_path / "missing-receipt",
        )

    fake = _fake_executable(tmp_path / "solver")
    receipt_path = tmp_path / "receipt.json"
    _create_receipt(prepared, raw, receipt_path, fake)
    raw["cp2k_band"].write_text(raw["cp2k_band"].read_text() + "changed\n")
    with pytest.raises(ValueError, match="hash mismatch"):
        e3b.collect_evidence(
            prepared_root=prepared,
            raw_outputs={MODEL: raw},
            execution_receipts={MODEL: receipt_path},
            out=tmp_path / "changed-output",
        )


def test_collect_rejects_cross_model_swapped_logs_even_with_fresh_receipt(
    tmp_path: Path,
) -> None:
    prepared, _ = _write_prepared_fixture(tmp_path)
    swapped = _write_raw_fixture(tmp_path / "raw", model="other")
    fake = _fake_executable(tmp_path / "solver")
    receipt_path = tmp_path / "receipt.json"
    _create_receipt(prepared, swapped, receipt_path, fake)
    with pytest.raises(ValueError, match="project"):
        e3b.collect_evidence(
            prepared_root=prepared,
            raw_outputs={MODEL: swapped},
            execution_receipts={MODEL: receipt_path},
            out=tmp_path / "swapped",
        )


def test_analysis_reopens_source_files_and_applies_corrected_outcome(
    tmp_path: Path,
) -> None:
    prepared, preparation, evidence_root, evidence = _collect_bundle(tmp_path)
    result = e3b.analyze_models(
        preparation,
        evidence,
        prepared_root=prepared,
        evidence_root=evidence_root,
        preparation_sha256=e3b.sha256(prepared / "preparation.json"),
        transfer_tolerance_kcal_mol=5.0,
        endorsement_tolerance_kcal_mol=5.0,
    )
    model = result["models"][MODEL]
    assert model["typed_outcome"] == "complete-corrected"
    assert model["matched_cell_transfer_delta_kcal_mol"] == pytest.approx(2.0)
    assert model["matched_cell_correction_kcal_mol"] == pytest.approx(10.0)
    assert model["calibrated_barrier_kcal_mol"] == pytest.approx(110.0)
    assert model["endorsement_delta_kcal_mol"] == pytest.approx(10.0)
    assert model["verdict"] == "corrected"


def test_small_cp2k_correction_is_endorsed(tmp_path: Path) -> None:
    prepared, preparation, evidence_root, evidence = _collect_bundle(
        tmp_path, cp2k_barrier=105.0
    )
    result = e3b.analyze_models(
        preparation,
        evidence,
        prepared_root=prepared,
        evidence_root=evidence_root,
        preparation_sha256=e3b.sha256(prepared / "preparation.json"),
        endorsement_tolerance_kcal_mol=5.0,
    )
    model = result["models"][MODEL]
    assert model["typed_outcome"] == "complete-endorsed"
    assert model["matched_cell_correction_kcal_mol"] == pytest.approx(3.0)
    assert model["calibrated_barrier_kcal_mol"] == pytest.approx(103.0)
    assert model["endorsement_delta_kcal_mol"] == pytest.approx(3.0)
    assert model["verdict"] == "endorsed"


def test_analysis_rehashes_and_reparses_e3a_source_artifacts(tmp_path: Path) -> None:
    prepared, preparation, evidence_root, evidence = _collect_bundle(tmp_path)
    source_neb = Path(preparation["source"]["campaign_root"]) / MODEL / "screen.log"
    source_neb.write_text(source_neb.read_text() + "tampered\n")
    result = e3b.analyze_models(
        preparation,
        evidence,
        prepared_root=prepared,
        evidence_root=evidence_root,
        preparation_sha256=e3b.sha256(prepared / "preparation.json"),
    )
    model = result["models"][MODEL]
    assert model["typed_outcome"] == "incomplete-evidence-integrity"
    assert "hash mismatch" in model["reason"]
    assert model["calibrated_barrier_kcal_mol"] is None


def test_analysis_rejects_tampered_stored_parser_fields(tmp_path: Path) -> None:
    prepared, preparation, evidence_root, evidence = _collect_bundle(tmp_path)
    evidence["models"][MODEL]["cp2k"]["barrier_kcal_mol"] = 1.0
    result = e3b.analyze_models(
        preparation,
        evidence,
        prepared_root=prepared,
        evidence_root=evidence_root,
        preparation_sha256=e3b.sha256(prepared / "preparation.json"),
    )
    assert result["models"][MODEL]["typed_outcome"] == "incomplete-evidence-integrity"


def _write_smoke_prepared(root: Path) -> Path:
    prepared = root / "smoke-prepared"
    model_root = prepared / "models" / "reconstructed-replication"
    data = prepared / "cp2k-data"
    model_root.mkdir(parents=True)
    data.mkdir()
    (model_root / "smoke.inp").write_text(
        "&GLOBAL\n RUN_TYPE ENERGY_FORCE\n&END GLOBAL\n"
        "BASIS_SET_FILE_NAME ../../cp2k-data/BASIS_SET\n"
        "COORD_FILE_NAME initial.xyz\n"
    )
    (model_root / "initial.xyz").write_text("0\nfixture\n")
    (data / "BASIS_SET").write_text("fixture basis\n")
    e3b.write_json(
        prepared / "preparation.json",
        {
            "schema": "e3b-cp2k-preparation-v1",
            "method": {
                "cp2k_version": "2024.3",
                "cp2k_image_digest": "sha256:" + "a" * 64,
            },
            "models": {
                "reconstructed-replication": {"atom_identity_sha256": "identity"}
            },
        },
    )
    e3b.write_manifest(prepared)
    return prepared


def test_smoke_receipt_labels_fake_and_unenforced_declarations(tmp_path: Path) -> None:
    prepared = _write_smoke_prepared(tmp_path)
    fake = _fake_executable(
        tmp_path / "fake-cp2k",
        "print('SCF run converged in 3 steps')\n"
        "print('ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: -42.25')\n"
        "print('PROGRAM ENDED AT 2026-09-14')\n",
    )
    executable_hash_before = e3b.sha256(fake)
    runtime = tmp_path / "runtime"
    receipt = e3b.run_smoke(
        executable=str(fake),
        executable_identity="test/fake",
        prepared_root=prepared,
        model="reconstructed-replication",
        runtime_dir=runtime,
        receipt_path=runtime / "smoke-result.json",
        timeout_seconds=5.0,
        mpi_ranks=4,
        omp_threads=2,
        memory_limit_mib=512,
        execution_method="fixture direct launch",
    )
    assert receipt["schema"] == "e3b-cp2k-smoke-v3"
    assert receipt["status"] == "converged"
    executable = receipt["method"]["resolved_executable_before_launch"]
    assert executable["identity"] == "test/fake"
    assert executable["sha256"] == executable_hash_before
    assert receipt["resource_declarations"]["mpi_ranks"] == {
        "value": 4,
        "verification": "declared-not-enforced-by-generic-runner",
    }
    assert receipt["resource_declarations"]["memory_limit_mib"]["verification"] == (
        "declared-not-enforced-by-generic-runner"
    )
    assert receipt["environment_controls"]["omp_threads"]["mechanism"] == (
        "OMP_NUM_THREADS"
    )
    assert receipt["external_controller_evidence"] is None
    assert receipt["cleanup"]["internal_process_group"]["new_session"] is True
    assert receipt["cleanup"]["internal_process_group"]["subprocess_reaped"] is True
    assert receipt["cleanup"]["work_directory_removed"] is True


def test_smoke_timeout_terminates_descendant_process_group(tmp_path: Path) -> None:
    prepared = _write_smoke_prepared(tmp_path)
    marker = tmp_path / "descendant-survived"
    child_code = (
        "import signal,time,pathlib; "
        "signal.signal(signal.SIGTERM, signal.SIG_IGN); "
        "time.sleep(1.2); "
        f"pathlib.Path({str(marker)!r}).write_text('survived')"
    )
    fake = _fake_executable(
        tmp_path / "timeout-cp2k",
        "import subprocess, sys, time\n"
        f"subprocess.Popen([sys.executable, '-c', {child_code!r}])\n"
        "print('started descendant', flush=True)\n"
        "time.sleep(10)\n",
    )
    runtime = tmp_path / "timeout-runtime"
    receipt = e3b.run_smoke(
        executable=str(fake),
        executable_identity="test/fake",
        prepared_root=prepared,
        model="reconstructed-replication",
        runtime_dir=runtime,
        receipt_path=runtime / "smoke-result.json",
        timeout_seconds=0.1,
    )
    assert receipt["status"] == "incomplete-timeout"
    cleanup = receipt["cleanup"]["internal_process_group"]
    assert cleanup["timeout_term_sent"] is True
    assert cleanup["timeout_kill_sent"] is True
    assert cleanup["subprocess_reaped"] is True
    time.sleep(1.3)
    assert not marker.exists()


def test_typed_controller_readback_is_recorded_only_when_verified(
    tmp_path: Path,
) -> None:
    prepared = _write_smoke_prepared(tmp_path)
    fake = _fake_executable(
        tmp_path / "fake-cp2k",
        "print('SCF run converged in 1 steps')\n"
        "print('ENERGY| Total FORCE_EVAL ( QS ) energy [a.u.]: -1')\n"
        "print('PROGRAM ENDED AT now')\n",
    )
    readback = {"unit": "e3b-smoke.scope", "MemoryMax": 536870912}
    controller = {
        "schema": "e3b-controller-readback-v1",
        "controller": "systemd",
        "verified": True,
        "readback": readback,
        "readback_sha256": e3b.canonical_sha256(readback),
    }
    runtime = tmp_path / "controller-runtime"
    receipt = e3b.run_smoke(
        executable=str(fake),
        executable_identity="test/fake",
        prepared_root=prepared,
        model="reconstructed-replication",
        runtime_dir=runtime,
        receipt_path=runtime / "smoke-result.json",
        timeout_seconds=2.0,
        controller_evidence=controller,
    )
    assert receipt["external_controller_evidence"] == controller


def test_al_causality_rejects_lost_other_oxygen_and_gained_replacement() -> None:
    atom = e3b.contract.nteme_neb.Atom
    cell = e3b.contract.nteme_neb.Cell(30.0, 30.0, 30.0, 90.0, 90.0, 90.0)
    moving = atom(3, "Ar", 0.0, 0.0, 3.0)
    al = atom(10, "Al1", 0.0, 0.0, 0.0)
    parent_oxygen = [
        atom(76, "O1", 1.0, 0.0, 0.0),
        atom(75, "O1", -1.0, 0.0, 0.0),
        atom(1, "O2", 0.0, 1.0, 0.0),
        atom(2, "O2", 0.0, -1.0, 0.0),
        atom(4, "O3", 0.0, 0.0, 1.0),
        atom(5, "O3", 0.0, 0.0, -1.0),
    ]
    parent = e3b.contract.Model(
        name="reconstructed-replication",
        atoms=tuple([moving, al, *parent_oxygen]),
        cell=cell,
        metadata={},
        endpoint_atoms=None,
    )
    transformed_oxygen = [
        replace(item, element="O2") if item.id == 75 else item
        for item in parent_oxygen
        if item.id not in {76, 5}
    ]
    transformed_oxygen.append(atom(99, "O2", 0.5, 0.5, 0.0))
    endpoint_moving = replace(moving, z=2.0)
    transformed_atoms = tuple([moving, al, *transformed_oxygen])
    transformed = e3b.contract.Model(
        name="dehydroxylate-lattice",
        atoms=transformed_atoms,
        cell=cell,
        metadata={
            "transformation": {
                "selected_oh_oxygen_ids": [75, 76],
                "residual_oxygen_site_id": 75,
                "removed_site_ids": [76],
            },
            "route": {"moving_site_id": 3},
        },
        endpoint_atoms=tuple(
            endpoint_moving if item.id == 3 else item for item in transformed_atoms
        ),
    )
    gate = e3b._coordination_gate(transformed, parent)
    candidate = gate["five_coordinate_al_criterion"]["candidates"][0]
    assert candidate["parent_o_coordination"] == 6
    assert candidate["transformed_o_coordination"] == 5
    assert candidate["removed_hydroxyl_o_in_parent_shell"] is True
    assert candidate["lost_o_neighbor_source_ids"] == [5, 76]
    assert candidate["gained_o_neighbor_source_ids"] == [99]
    assert candidate["exact_removed_o_neighbor_change"] is False
    assert candidate["coordination_change_caused_by_transformation"] is False
    assert gate["coordination_pass"] is False


def test_collect_and_analyze_cli_use_typed_receipt_spec(tmp_path: Path) -> None:
    prepared, _ = _write_prepared_fixture(tmp_path)
    raw = _write_raw_fixture(tmp_path / "raw")
    fake = _fake_executable(tmp_path / "solver")
    receipt_path = tmp_path / "receipt.json"
    raw_spec = tmp_path / "raw-spec.json"
    e3b.write_json(
        raw_spec,
        {
            "schema": "e3b-raw-output-spec-v1",
            "models": {
                MODEL: {
                    "execution_receipt": str(receipt_path),
                    "outputs": {role: str(path) for role, path in raw.items()},
                }
            },
        },
    )
    execution_spec = tmp_path / "execution-spec.json"
    e3b.write_json(
        execution_spec,
        {
            "schema": "e3b-execution-declarations-v1",
            "executables": {
                solver: {"path": str(fake), "identity": "test/fake"}
                for solver in ("cp2k", "lammps")
            },
            "units": {
                unit: {
                    "command": [str(fake), unit],
                    "resources": {"ranks": 1, "threads": 1},
                }
                for unit in e3b.SOLVER_UNITS
            },
        },
    )
    receipt_args = e3b.build_parser().parse_args(
        [
            "create-execution-receipt",
            "--prepared-root",
            str(prepared),
            "--model",
            MODEL,
            "--raw-spec",
            str(raw_spec),
            "--execution-spec",
            str(execution_spec),
            "--out",
            str(receipt_path),
        ]
    )
    assert receipt_args.func(receipt_args) == 0
    assert receipt_path.is_file()
    evidence_root = tmp_path / "evidence-cli"
    collect_args = e3b.build_parser().parse_args(
        [
            "collect-evidence",
            "--prepared-root",
            str(prepared),
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
            str(prepared),
            "--evidence-root",
            str(evidence_root),
            "--out",
            str(analysis_root),
        ]
    )
    assert analyze_args.func(analyze_args) == 0
    result = e3b.read_json(analysis_root / "analysis.json")
    assert result["models"][MODEL]["typed_outcome"] == "complete-corrected"
    e3b.verify_manifest(analysis_root)
