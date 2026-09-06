"""CPU-only orchestration and independent-verifier tests for A3e."""

from __future__ import annotations

import ast
import json
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from quarry.clusters import Cluster
from quarry.store import geometry_hash
from scripts import a3b_proton_microstate_stability as a3b
from scripts import a3e_oss_owner_basin as a3e
from scripts import a3e_verify as verifier


def _cluster() -> Cluster:
    symbols = ["He"] * 40
    symbols[a3e.O14] = "O"
    symbols[a3e.O21] = "O"
    symbols[a3e.H35] = "H"
    coords = np.array([[float(index * 5), 0.0, 0.0] for index in range(40)])
    coords[a3e.H35] = coords[a3e.O21] + np.array([0.96, 0.0, 0.0])
    return Cluster(
        "oss-neutral-n1",
        symbols,
        coords,
        charge=0,
        spin=0,
        frozen_indices=[0, a3e.O14, a3e.O21],
    )


def _source(cluster: Cluster | None = None) -> a3e.SourceEvidence:
    cluster = cluster or _cluster()
    return a3e.SourceEvidence(
        cluster=cluster,
        source_root=Path("/fixture/source"),
        hashes={
            name: str(index) * 64 for index, name in enumerate(a3e.SOURCE_PATHS, 1)
        },
        metadata={"site_kind": "Oss", "state": "neutral", "n_intact": 1},
    )


def _source_tree(tmp_path: Path) -> tuple[Path, dict[str, str], Cluster]:
    root = tmp_path / "source"
    for relative in a3e.SOURCE_PATHS.values():
        (root / relative).parent.mkdir(parents=True, exist_ok=True)
    cluster = _cluster()
    (root / a3e.SOURCE_PATHS["seed"]).write_text(cluster.to_xyz())
    (root / a3e.SOURCE_PATHS["child_log"]).write_text("failed owner H35:O21->O14\n")
    (root / a3e.SOURCE_PATHS["launch"]).write_text("launch\n")
    (root / a3e.SOURCE_PATHS["restoration"]).write_text("restored\n")
    (root / a3e.SOURCE_PATHS["family_receipt"]).write_text(
        json.dumps(
            {
                "schema": "a3-family-campaign-terminal-v1",
                "success": False,
                "family": "oss",
                "state": "neutral",
                "current_n_intact": 1,
                "completed": [],
                "expected_git_sha": a3e.EXPECTED_EXECUTION_SOURCE,
                "observed_git_sha": a3e.EXPECTED_EXECUTION_SOURCE,
            },
            sort_keys=True,
        )
    )
    (root / a3e.SOURCE_PATHS["metadata"]).write_text(
        json.dumps(
            {
                "site_kind": "Oss",
                "state": "neutral",
                "n_intact": 1,
                "charge": 0,
                "method": "b3lyp/def2-svp/df",
                "driver_git_commit": a3e.EXPECTED_EXECUTION_SOURCE,
                "n_atoms": 37,
            },
            sort_keys=True,
        )
    )
    hashes = {
        name: a3e.sha256_path(root / relative)
        for name, relative in a3e.SOURCE_PATHS.items()
    }
    return root, hashes, cluster


def _patch_calculators(
    monkeypatch, endpoints, energies=None, converged=None, imag=None
):
    calls = []
    converged = converged or [True] * len(endpoints)
    results = iter(zip(endpoints, converged, strict=True))

    def fake_optimize(cluster, settings, **kwargs):
        calls.append((cluster, settings, kwargs))
        endpoint, did_converge = next(results)
        return SimpleNamespace(cluster=endpoint, converged=did_converge, max_steps=100)

    energy_values = energies or {endpoint.name: -10.0 for endpoint in endpoints}
    monkeypatch.setattr(a3e, "optimize_bounded", fake_optimize)
    monkeypatch.setattr(
        a3e, "energy", lambda cluster, _settings: energy_values[cluster.name]
    )
    monkeypatch.setattr(
        a3b, "gradient", lambda cluster, _settings: np.zeros_like(cluster.coords)
    )
    imaginary = np.asarray([] if imag is None else imag, dtype=float)
    monkeypatch.setattr(
        a3e,
        "frequencies",
        lambda cluster, settings: SimpleNamespace(
            imaginary_cm=imaginary,
            electronic_hartree=energy_values[cluster.name],
            geometry_fingerprint=a3e.frequency_geometry_fingerprint(cluster),
            settings_fingerprint=a3e.frequency_settings_fingerprint(settings),
        ),
    )
    return calls


def _endpoints(source: Cluster, *, transfer: bool = False):
    values = []
    for index, name in enumerate(("conditioned", "constrained", "released"), 1):
        coords = source.coords.copy()
        coords[1, 1] = 0.01 * index
        values.append(replace(source, name=name, coords=coords))
    if transfer:
        coords = values[-1].coords.copy()
        coords[a3e.H35] = source.coords[a3e.O14] + np.array([0.96, 0.0, 0.0])
        values[-1] = replace(values[-1], coords=coords)
    return values


def test_validate_source_rehashes_full_contract_and_binds_seed(tmp_path, monkeypatch):
    root, hashes, cluster = _source_tree(tmp_path)
    monkeypatch.setattr(a3e, "canonical_template", lambda _root: cluster)
    source = a3e.validate_source(root, expected_hashes=hashes, repo_root=tmp_path)
    assert source.hashes == hashes
    assert a3e._owner_labels(source.cluster) == ["H35:O21"]
    assert geometry_hash(source.cluster.to_xyz()) == geometry_hash(cluster.to_xyz())

    (root / a3e.SOURCE_PATHS["child_log"]).write_text("tampered\n")
    with pytest.raises(RuntimeError, match="child_log SHA-256 mismatch"):
        a3e.validate_source(root, expected_hashes=hashes, repo_root=tmp_path)


def test_exact_three_stage_budget_and_fresh_release(tmp_path, monkeypatch):
    source = _source()
    endpoints = _endpoints(source.cluster)
    calls = _patch_calculators(monkeypatch, endpoints)
    terminal = a3e.run_experiment(
        tmp_path,
        source,
        code_revision="a" * 40,
    )
    assert terminal["classification"] == a3e.ACCEPTED_CANDIDATE
    assert len(calls) == 3
    expected_constraint = [(a3e.O21, a3e.H35, pytest.approx(0.96))]
    assert calls[0][2] == {"max_steps": 100, "fixed_distances": expected_constraint}
    assert calls[0][1].xc == "hf"
    assert calls[0][1].basis == "sto-3g"
    assert calls[1][2] == {"max_steps": 100, "fixed_distances": expected_constraint}
    assert calls[1][1].xc == "b3lyp"
    assert calls[1][1].density_fit is True
    assert calls[2][2] == {"max_steps": 100}
    for spec in a3e.STAGES:
        receipt = json.loads(
            (tmp_path / "stages" / spec.directory / "receipt.json").read_text()
        )
        assert receipt["optimizer"]["fresh_instance"] is True
        assert receipt["optimizer"]["geometric_default_fresh_hessian"] is True
        assert receipt["signature"]["budget"] == {
            "max_steps": 100,
            "retry_allowed": False,
        }
        assert (tmp_path / "stages" / spec.directory / "raw-endpoint.xyz").is_file()


def test_constrained_gate_blocks_release_without_second_production_budget(
    tmp_path, monkeypatch
):
    source = _source()
    endpoints = _endpoints(source.cluster)[:2]
    calls = _patch_calculators(monkeypatch, endpoints, converged=[True, False])
    terminal = a3e.run_experiment(tmp_path, source, code_revision="b" * 40)
    assert len(calls) == 2
    assert terminal["classification"] == a3e.INCONCLUSIVE_CANDIDATE
    assert terminal["stages"]["released-production"]["status"] == "not-run"
    assert not (tmp_path / "stages" / a3e.STAGES[2].directory).exists()


def test_production_no_basin_requires_exact_h35_transfer_and_downhill_energy(
    tmp_path, monkeypatch
):
    source = _source()
    endpoints = _endpoints(source.cluster, transfer=True)
    calls = _patch_calculators(
        monkeypatch,
        endpoints,
        energies={"conditioned": -9.0, "constrained": -10.0, "released": -10.01},
    )
    terminal = a3e.run_experiment(tmp_path, source, code_revision="c" * 40)
    assert len(calls) == 3
    assert terminal["classification"] == a3e.NO_BASIN_CANDIDATE
    released = json.loads(
        (tmp_path / "stages" / a3e.STAGES[2].directory / "receipt.json").read_text()
    )
    assert released["owner_changes"] == ["H35:O21->O14"]
    assert released["phva"]["status"] == "not-required"


def test_no_basin_strict_downhill_boundary_is_inconclusive():
    constrained = {
        "status": "complete",
        "optimizer": {"converged": True},
        "stationary": True,
        "owner_retaining": True,
        "owner_changes": [],
        "observed_owners": ["H35:O21"],
        "structure": {
            "constraint_residuals": [
                {
                    "oxygen": a3e.O21,
                    "hydrogen": a3e.H35,
                    "absolute_residual_a": 0.0,
                }
            ]
        },
        "evidence": {"energy_hartree": -10.0, "passed": True},
    }
    released = {
        "status": "complete",
        "stationary": True,
        "owner_retaining": False,
        "owner_changes": ["H35:O21->O14"],
        "observed_owners": ["H35:O14"],
        "evidence": {
            "energy_hartree": -10.0 - a3e.DOWNHILL_MIN_HARTREE,
            "passed": True,
        },
        "phva": {"status": "not-required"},
    }
    assert a3e.classify(constrained, released) == a3e.INCONCLUSIVE_CANDIDATE


def test_raw_endpoint_and_pending_receipt_precede_structural_gate(
    tmp_path, monkeypatch
):
    source = _source()
    endpoints = _endpoints(source.cluster)
    _patch_calculators(monkeypatch, endpoints)
    real_gate = a3b.structural_gate
    observations = []

    def checking_gate(raw, reference, active, stage):
        spec = a3e.STAGE_BY_ID[stage]
        stage_dir = tmp_path / "stages" / spec.directory
        observations.append(
            (
                (stage_dir / "raw-endpoint.xyz").is_file(),
                json.loads((stage_dir / "receipt.json").read_text())["status"],
            )
        )
        return real_gate(raw, reference, active, stage)

    monkeypatch.setattr(a3b, "structural_gate", checking_gate)
    a3e.run_experiment(tmp_path, source, code_revision="d" * 40)
    assert observations == [(True, "pending-structural-gates")] * 6


def test_reservation_prevents_replay_after_orphaned_budget(tmp_path, monkeypatch):
    source = _source()
    stage_dir = tmp_path / "stages" / a3e.STAGES[0].directory
    stage_dir.mkdir(parents=True)
    (stage_dir / "reservation.json").write_text(json.dumps({"status": "reserved"}))
    monkeypatch.setattr(
        a3e,
        "optimize_bounded",
        lambda *_args, **_kwargs: pytest.fail("orphaned budget was replayed"),
    )
    terminal = a3e.run_experiment(tmp_path, source, code_revision="e" * 40)
    assert terminal["classification"] == a3e.INCONCLUSIVE_CANDIDATE
    assert terminal["stages"]["owner-conditioning"]["status"] == "orphaned-stage"


def test_independent_verifier_recomputes_and_detects_tampering(tmp_path, monkeypatch):
    tree = ast.parse(Path(verifier.__file__).read_text())
    imported = {
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, (ast.Import, ast.ImportFrom))
        for alias in node.names
    }
    assert not any("optim" in name for name in imported)

    source = _source()
    endpoints = _endpoints(source.cluster)
    energy_values = {endpoint.name: -10.0 for endpoint in endpoints}
    _patch_calculators(monkeypatch, endpoints, energies=energy_values)
    a3e.run_experiment(tmp_path, source, code_revision="f" * 40)
    monkeypatch.setattr(
        verifier, "energy", lambda cluster, _settings: energy_values[cluster.name]
    )
    monkeypatch.setattr(
        verifier,
        "frequencies",
        lambda cluster, settings: SimpleNamespace(
            imaginary_cm=np.asarray([], dtype=float),
            electronic_hartree=energy_values[cluster.name],
            geometry_fingerprint=a3e.frequency_geometry_fingerprint(cluster),
            settings_fingerprint=a3e.frequency_settings_fingerprint(settings),
        ),
    )
    result = verifier.verify_experiment(
        tmp_path,
        source_override=source,
        verifier_identity="cold-worker-2",
    )
    assert result["status"] == "verified"
    assert result["classification"] == a3e.VERIFIED_ACCEPTED
    assert result["verifier_identity"] != result["executor_identity"]

    released = tmp_path / "stages" / a3e.STAGES[2].directory / "endpoint.xyz"
    released.write_text(released.read_text().replace("released", "tampered"))
    rejected = verifier.verify_experiment(
        tmp_path,
        source_override=source,
        verifier_identity="cold-worker-2",
    )
    assert rejected["status"] == "rejected"
    assert "hash-mismatched" in rejected["detail"]


def test_verifier_refuses_executor_identity_and_no_recompute(tmp_path, monkeypatch):
    source = _source()
    endpoints = _endpoints(source.cluster)
    _patch_calculators(monkeypatch, endpoints)
    a3e.run_experiment(tmp_path, source, code_revision="1" * 40)
    same = verifier.verify_experiment(
        tmp_path,
        source_override=source,
        verifier_identity=a3e.EXECUTOR_IDENTITY,
    )
    assert same["status"] == "rejected"
    refused = verifier.verify_experiment(
        tmp_path,
        source_override=source,
        verifier_identity="cold-worker-2",
        recompute_calculators=False,
    )
    assert refused["status"] == "rejected"
    assert "calculator recomputation" in refused["detail"]


def test_verifier_rejects_coordinated_raw_endpoint_and_receipt_tampering(
    tmp_path, monkeypatch
):
    source = _source()
    endpoints = _endpoints(source.cluster)
    energy_values = {endpoint.name: -10.0 for endpoint in endpoints}
    _patch_calculators(monkeypatch, endpoints, energies=energy_values)
    a3e.run_experiment(tmp_path, source, code_revision="g" * 40)
    monkeypatch.setattr(
        verifier, "energy", lambda cluster, _settings: energy_values[cluster.name]
    )
    monkeypatch.setattr(
        verifier,
        "frequencies",
        lambda cluster, settings: SimpleNamespace(
            imaginary_cm=np.asarray([], dtype=float),
            electronic_hartree=energy_values[cluster.name],
            geometry_fingerprint=a3e.frequency_geometry_fingerprint(cluster),
            settings_fingerprint=a3e.frequency_settings_fingerprint(settings),
        ),
    )

    spec = a3e.STAGES[-1]
    stage_dir = tmp_path / "stages" / spec.directory
    raw_path = stage_dir / "raw-endpoint.xyz"
    tampered = a3b.read_cluster(raw_path, source.cluster)
    coords = tampered.coords.copy()
    coords[0] = coords[0] + np.array([0.05, 0.0, 0.0])
    tampered = replace(tampered, coords=coords)
    a3e.atomic_xyz(raw_path, tampered)
    receipt_path = stage_dir / "receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["raw_endpoint"] = a3e._artifact(tampered, raw_path)
    a3e.atomic_json(receipt_path, receipt)
    candidate_path = tmp_path / "candidate-terminal.json"
    candidate = json.loads(candidate_path.read_text())
    candidate["stages"][spec.stage_id]["receipt_sha256"] = a3e.sha256_path(receipt_path)
    a3e.atomic_json(candidate_path, candidate)

    result = verifier.verify_experiment(
        tmp_path,
        source_override=source,
        verifier_identity="cold-worker-2",
    )
    assert result["status"] == "rejected"
    assert "raw" in result["detail"].casefold()


def test_no_forbidden_downstream_outputs(tmp_path, monkeypatch):
    source = _source()
    endpoints = _endpoints(source.cluster)
    _patch_calculators(monkeypatch, endpoints)
    a3e.run_experiment(tmp_path, source, code_revision="2" * 40)
    forbidden = {"results.json", "store.sqlite", "ts.xyz", "barrier.json", "petra.toml"}
    assert not any(path.name in forbidden for path in tmp_path.rglob("*"))
