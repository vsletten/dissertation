"""CPU-only orchestration and adjudication tests for A3b."""

from __future__ import annotations

import ast
import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from quarry.clusters import Cluster
from quarry.pipeline import DftSettings
from quarry.store import geometry_hash
from scripts import a3b_proton_microstate_stability as a3b
from scripts import a3b_verify as verifier


def _cluster() -> Cluster:
    symbols = ["He"] * 58
    for index in (20, 26, 29, 31):
        symbols[index] = "O"
    symbols[50] = "H"
    symbols[57] = "H"
    coords = np.array([[float(index * 5), 0.0, 0.0] for index in range(58)])
    coords[50] = coords[26] + np.array([0.96, 0.0, 0.0])
    coords[57] = coords[29] + np.array([0.96, 0.0, 0.0])
    return Cluster(
        "conditioned",
        symbols,
        coords,
        charge=0,
        spin=0,
        frozen_indices=[0, 20, 26, 29, 31],
    )


def _source(cluster: Cluster | None = None) -> a3b.SourceEvidence:
    cluster = cluster or _cluster()
    return a3b.SourceEvidence(
        cluster=cluster,
        source_root=Path("/fixture/a3a"),
        closeout_manifest_path=Path("/fixture/a3a/closeout-manifest.json"),
        closeout_manifest_sha256=a3b.EXPECTED_CLOSEOUT_MANIFEST_SHA256,
        conditioned_xyz_path=Path("/fixture/a3a/conditioned.xyz"),
        conditioned_xyz_sha256=a3b.EXPECTED_CONDITIONED_XYZ_SHA256,
        conditioned_receipt_path=Path("/fixture/a3a/conditioned.json"),
        conditioned_receipt_sha256=a3b.EXPECTED_CONDITIONED_RECEIPT_SHA256,
        conditioned_geometry_hash=geometry_hash(cluster.to_xyz()),
    )


def _source_tree(tmp_path: Path) -> tuple[Path, dict[str, str]]:
    source_root = tmp_path / "a3a"
    run_dir = source_root / a3b.SOURCE_RUN_RELATIVE
    run_dir.mkdir(parents=True)
    cluster = _cluster()
    xyz_path = run_dir / "complex_conditioned.xyz"
    xyz_path.write_text(cluster.to_xyz())
    receipt_path = run_dir / "complex_conditioned.json"
    receipt_path.write_text(
        json.dumps(
            {
                "schema": "phase2-reactant-conditioning-v1",
                "endpoint_geometry_hash": geometry_hash(cluster.to_xyz()),
                "signature": {
                    "symbols": cluster.symbols,
                    "charge": cluster.charge,
                    "spin": cluster.spin,
                    "frozen_indices": cluster.frozen_indices,
                    "oxygen_proton_owners": ["H50:O26", "H57:O29"],
                },
            },
            sort_keys=True,
        )
    )
    hashes = {
        "xyz": a3b.sha256_path(xyz_path),
        "receipt": a3b.sha256_path(receipt_path),
    }
    manifest_path = source_root / "closeout-manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "schema": "a3a-closeout-manifest-v1",
                "status": "failed",
                "preserved_evidence": [
                    {"path": str(xyz_path), "sha256": hashes["xyz"]},
                    {"path": str(receipt_path), "sha256": hashes["receipt"]},
                ],
            },
            sort_keys=True,
        )
    )
    hashes["manifest"] = a3b.sha256_path(manifest_path)
    return source_root, hashes


def _patch_calculators(
    monkeypatch, endpoints, converged=None, energies=None, imag=None
):
    calls = []
    converged = converged or [True] * len(endpoints)
    results = iter(zip(endpoints, converged, strict=True))

    def fake_optimize(cluster, settings, **kwargs):
        calls.append((cluster, settings, kwargs))
        endpoint, did_converge = next(results)
        return SimpleNamespace(cluster=endpoint, converged=did_converge, max_steps=100)

    monkeypatch.setattr(a3b, "optimize_bounded", fake_optimize)
    monkeypatch.setattr(
        a3b, "gradient", lambda cluster, _settings: np.zeros_like(cluster.coords)
    )
    energy_map = energies or {endpoint.name: -10.0 for endpoint in endpoints}
    monkeypatch.setattr(
        a3b, "energy", lambda cluster, _settings: energy_map[cluster.name]
    )
    imaginary = np.asarray([] if imag is None else imag, dtype=float)
    monkeypatch.setattr(
        a3b,
        "frequencies",
        lambda cluster, _settings: SimpleNamespace(
            imaginary_cm=imaginary,
            electronic_hartree=energy_map[cluster.name],
            geometry_fingerprint=a3b.frequency_geometry_fingerprint(cluster),
            settings_fingerprint=a3b.frequency_settings_fingerprint(
                DftSettings(xc="b3lyp", basis="def2-svp", density_fit=True)
            ),
        ),
    )
    return calls


def _named_endpoints(source: Cluster, *, transfer_finals: bool = False):
    endpoints = []
    for index in range(5):
        coords = source.coords.copy()
        coords[1, 1] = 0.01 * (index + 1)
        endpoint = replace(source, name=f"stage-{index + 1}", coords=coords)
        endpoints.append(endpoint)
    if transfer_finals:
        for index in (2, 4):
            coords = endpoints[index].coords.copy()
            coords[50] = source.coords[31] + np.array([0.96, 0.0, 0.0])
            coords[57] = source.coords[20] + np.array([0.96, 0.0, 0.0])
            endpoints[index] = replace(endpoints[index], coords=coords)
    return endpoints


def test_validate_source_rehashes_all_three_pinned_artifacts(tmp_path):
    source_root, hashes = _source_tree(tmp_path)
    evidence = a3b.validate_source(
        source_root,
        expected_closeout_sha256=hashes["manifest"],
        expected_conditioned_xyz_sha256=hashes["xyz"],
        expected_conditioned_receipt_sha256=hashes["receipt"],
    )
    assert evidence.cluster.symbols == _cluster().symbols
    assert evidence.conditioned_geometry_hash == geometry_hash(_cluster().to_xyz())

    conditioned = source_root / a3b.SOURCE_RUN_RELATIVE / "complex_conditioned.xyz"
    conditioned.write_text(
        conditioned.read_text().replace("130.96000000", "130.95000000")
    )
    with pytest.raises(RuntimeError, match="conditioned XYZ SHA-256 mismatch"):
        a3b.validate_source(
            source_root,
            expected_closeout_sha256=hashes["manifest"],
            expected_conditioned_xyz_sha256=hashes["xyz"],
            expected_conditioned_receipt_sha256=hashes["receipt"],
        )


def test_five_stage_dag_has_exact_calls_constraints_and_fresh_budgets(
    tmp_path, monkeypatch
):
    source = _source()
    endpoints = _named_endpoints(source.cluster)
    calls = _patch_calculators(monkeypatch, endpoints)

    terminal = a3b.run_experiment(tmp_path, source, use_gpu=False)

    assert terminal["classification"] == a3b.RECOVERY_CANDIDATE
    assert len(calls) == 5
    h50 = (26, 50, pytest.approx(0.96))
    h57 = (29, 57, pytest.approx(0.96))
    assert calls[0][2] == {"max_steps": 100, "fixed_distances": [h50, h57]}
    assert calls[1][2] == {"max_steps": 100, "fixed_distances": [h57]}
    assert calls[2][2] == {"max_steps": 100}
    assert calls[3][2] == {"max_steps": 100, "fixed_distances": [h50]}
    assert calls[4][2] == {"max_steps": 100}
    assert calls[1][0].name == "stage-1"
    assert calls[2][0].name == "stage-2"
    assert calls[3][0].name == "stage-1"
    assert calls[4][0].name == "stage-4"
    for spec in a3b.STAGES:
        receipt = json.loads(
            (tmp_path / "stages" / spec.directory / "receipt.json").read_text()
        )
        assert receipt["budget"] == {"max_steps": 100, "continuation_allowed": False}
        assert receipt["optimizer"]["fresh_instance"] is True
        assert receipt["optimizer"]["geometric_default_fresh_hessian"] is True


def test_nonconverged_first_route_does_not_suppress_sibling(tmp_path, monkeypatch):
    source = _source()
    endpoints = _named_endpoints(source.cluster)
    calls = _patch_calculators(
        monkeypatch, endpoints, converged=[True, False, False, True, True]
    )
    terminal = a3b.run_experiment(tmp_path, source, use_gpu=False)
    assert len(calls) == 5
    assert terminal["classification"] == a3b.RECOVERY_CANDIDATE
    assert terminal["stages"]["h50-then-h57-unconstrained"] == "complete"
    assert terminal["stages"]["h57-then-h50-unconstrained"] == "complete"


def test_orphan_reservation_is_never_retried(tmp_path, monkeypatch):
    source = _source()
    stage_dir = tmp_path / "stages" / a3b.STAGES[0].directory
    stage_dir.mkdir(parents=True)
    (stage_dir / "reservation.json").write_text(json.dumps({"status": "reserved"}))
    monkeypatch.setattr(
        a3b,
        "optimize_bounded",
        lambda *_args, **_kwargs: pytest.fail("orphaned budget was retried"),
    )
    terminal = a3b.run_experiment(tmp_path, source, use_gpu=False)
    assert terminal["classification"] == a3b.INCONCLUSIVE
    assert terminal["stages"]["common-dual"] == "orphaned-reservation"


def test_raw_endpoint_and_pending_receipt_exist_before_structural_gate(
    tmp_path, monkeypatch
):
    source = _source()
    endpoints = _named_endpoints(source.cluster)
    _patch_calculators(monkeypatch, endpoints)
    real_gate = a3b.structural_gate
    observed = []

    def checking_gate(raw, reference, constraints, stage):
        stage_dir = tmp_path / "stages" / a3b.STAGE_BY_ID[stage].directory
        observed.append(
            (
                (stage_dir / "raw-endpoint.xyz").is_file(),
                json.loads((stage_dir / "receipt.json").read_text())["status"],
            )
        )
        return real_gate(raw, reference, constraints, stage)

    monkeypatch.setattr(a3b, "structural_gate", checking_gate)
    a3b.run_experiment(tmp_path, source, use_gpu=False)
    # Each endpoint is gated once as the raw optimizer object and once after
    # the exact persisted/projected XYZ round trip.
    assert observed == [(True, "pending-structural-gates")] * 10


@pytest.mark.parametrize(
    (
        "retain_a",
        "retain_b",
        "stationary_a",
        "stationary_b",
        "imag_a",
        "imag_b",
        "ea",
        "eb",
        "expected",
    ),
    [
        (True, False, True, True, [], [], -10.0, -10.0, a3b.RECOVERY_CANDIDATE),
        (True, True, True, True, [30.0], [99.0], -10.0, -10.0, a3b.RECOVERY_CANDIDATE),
        (
            False,
            False,
            True,
            True,
            [],
            [],
            -10.000002,
            -10.000002,
            a3b.NO_BASIN_CANDIDATE,
        ),
        (False, False, True, True, [], [], -10.000001, -10.000002, a3b.INCONCLUSIVE),
        (False, False, True, False, [], [], -10.000002, -10.000002, a3b.INCONCLUSIVE),
        (True, False, True, True, [30.0000001], [], -10.0, -10.0, a3b.INCONCLUSIVE),
    ],
)
def test_classification_matrix_and_strict_boundaries(
    retain_a,
    retain_b,
    stationary_a,
    stationary_b,
    imag_a,
    imag_b,
    ea,
    eb,
    expected,
):
    def receipt(retain, stationary, imaginary, energy):
        return {
            "status": "complete",
            "fully_unconstrained": True,
            "owner_retaining": retain,
            "stationary": stationary,
            "evidence": {"energy_hartree": energy},
            "phva": {"status": "passed", "imaginary_cm": imaginary},
        }

    common = {"status": "complete", "evidence": {"energy_hartree": -10.0}}
    assert (
        a3b.classify(
            common,
            receipt(retain_a, stationary_a, imag_a, ea),
            receipt(retain_b, stationary_b, imag_b, eb),
        )
        == expected
    )


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_classification_rejects_nonfinite_energy(value):
    final = {
        "status": "complete",
        "fully_unconstrained": True,
        "owner_retaining": False,
        "stationary": True,
        "evidence": {"energy_hartree": value},
        "phva": {"status": "not-required", "imaginary_cm": []},
    }
    common = {"status": "complete", "evidence": {"energy_hartree": -10.0}}
    assert a3b.classify(common, final, final) == a3b.INCONCLUSIVE


def test_verifier_has_no_optimizer_import_and_accepts_untampered_run(
    tmp_path, monkeypatch
):
    tree = ast.parse(Path(verifier.__file__).read_text())
    imported = {
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, (ast.Import, ast.ImportFrom))
        for alias in node.names
    }
    assert not any("optim" in name for name in imported)

    source = _source()
    endpoints = _named_endpoints(source.cluster)
    _patch_calculators(monkeypatch, endpoints)
    a3b.run_experiment(tmp_path, source, use_gpu=False)
    result = verifier.verify_experiment(
        tmp_path,
        source_override=verifier.StaticSource.from_executor_source(source),
    )
    assert result["status"] == "verified"
    assert result["classification"] == a3b.RECOVERY_CANDIDATE
    assert result["verifier_identity"] != result["executor_identity"]
    assert (tmp_path / "verified-terminal.json").is_file()


def test_verifier_detects_endpoint_receipt_and_parent_dag_tampering(
    tmp_path, monkeypatch
):
    source = _source()
    endpoints = _named_endpoints(source.cluster)
    _patch_calculators(monkeypatch, endpoints)
    a3b.run_experiment(tmp_path, source, use_gpu=False)
    receipt_path = tmp_path / "stages" / a3b.STAGES[2].directory / "receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["parent"]["stage"] = "common-dual"
    receipt_path.write_text(json.dumps(receipt))

    result = verifier.verify_experiment(
        tmp_path,
        source_override=verifier.StaticSource.from_executor_source(source),
    )
    assert result["status"] == "rejected"
    assert "parent" in result["detail"]


def test_resume_rejects_signature_drift_without_spending_new_budget(
    tmp_path, monkeypatch
):
    source = _source()
    endpoints = _named_endpoints(source.cluster)
    _patch_calculators(monkeypatch, endpoints)
    first = a3b.run_experiment(tmp_path, source, use_gpu=False)
    assert first["classification"] == a3b.RECOVERY_CANDIDATE

    monkeypatch.setattr(
        a3b,
        "optimize_bounded",
        lambda *_args, **_kwargs: pytest.fail("signature drift spent a new budget"),
    )
    second = a3b.run_experiment(tmp_path, source, use_gpu=True)
    assert second["classification"] == a3b.INCONCLUSIVE
    assert second["stages"]["common-dual"] == "reservation-identity-mismatch"


def test_verifier_rejects_candidate_source_tampering_with_static_fixture(
    tmp_path, monkeypatch
):
    source = _source()
    endpoints = _named_endpoints(source.cluster)
    _patch_calculators(monkeypatch, endpoints)
    a3b.run_experiment(tmp_path, source, use_gpu=False)
    candidate_path = tmp_path / "candidate-terminal.json"
    candidate = json.loads(candidate_path.read_text())
    candidate["source"]["conditioned_geometry_hash"] = "0" * 64
    candidate_path.write_text(json.dumps(candidate))

    result = verifier.verify_experiment(
        tmp_path,
        source_override=verifier.StaticSource.from_executor_source(source),
    )
    assert result["status"] == "rejected"
    assert "source binding" in result["detail"]


@pytest.mark.parametrize("tamper", ["parent-hash", "reservation"])
def test_verifier_rejects_parent_hash_and_reservation_tampering(
    tmp_path, monkeypatch, tamper
):
    source = _source()
    endpoints = _named_endpoints(source.cluster)
    _patch_calculators(monkeypatch, endpoints)
    a3b.run_experiment(tmp_path, source, use_gpu=False)
    spec = a3b.STAGES[1] if tamper == "parent-hash" else a3b.STAGES[0]
    if tamper == "parent-hash":
        path = tmp_path / "stages" / spec.directory / "receipt.json"
        payload = json.loads(path.read_text())
        payload["parent"]["receipt_sha256"] = "0" * 64
    else:
        path = tmp_path / "stages" / spec.directory / "reservation.json"
        payload = json.loads(path.read_text())
        payload["status"] = "spent"
    path.write_text(json.dumps(payload))

    result = verifier.verify_experiment(
        tmp_path,
        source_override=verifier.StaticSource.from_executor_source(source),
    )
    assert result["status"] == "rejected"
    expected = "parent" if tamper == "parent-hash" else "reservation"
    assert expected in result["detail"]


def test_no_forbidden_scientific_outputs_are_emitted(tmp_path, monkeypatch):
    source = _source()
    endpoints = _named_endpoints(source.cluster, transfer_finals=True)
    energy_values = {
        endpoint.name: -10.0 - index * 0.001 for index, endpoint in enumerate(endpoints)
    }
    _patch_calculators(monkeypatch, endpoints, energies=energy_values)
    terminal = a3b.run_experiment(tmp_path, source, use_gpu=False)
    assert terminal["classification"] == a3b.NO_BASIN_CANDIDATE
    forbidden = {"results.json", "store.sqlite", "ts.xyz", "barrier.json", "petra.toml"}
    assert not any(path.name in forbidden for path in tmp_path.rglob("*"))
