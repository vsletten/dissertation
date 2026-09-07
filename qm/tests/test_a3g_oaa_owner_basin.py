"""CPU-only orchestration and independent-verifier tests for A3g."""

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
from scripts import a3g_oaa_owner_basin as a3g
from scripts import a3g_verify as verifier


def _cluster() -> Cluster:
    symbols = ["He"] * 58
    symbols[a3g.O9] = "O"
    symbols[a3g.O15] = "O"
    symbols[a3g.H52] = "H"
    coords = np.array([[float(index * 5), 0.0, 0.0] for index in range(58)])
    coords[a3g.H52] = coords[a3g.O15] + np.array([0.96, 0.0, 0.0])
    return Cluster(
        "oaa-neutral-n2",
        symbols,
        coords,
        charge=0,
        spin=0,
        frozen_indices=[0, a3g.O9, a3g.O15],
    )


def _source(cluster: Cluster | None = None) -> a3g.SourceEvidence:
    cluster = cluster or _cluster()
    return a3g.SourceEvidence(
        cluster=cluster,
        source_root=Path("/fixture/source"),
        hashes={
            name: str(index) * 64 for index, name in enumerate(a3g.SOURCE_PATHS, 1)
        },
        metadata={"site_kind": "Oaa", "state": "neutral", "n_intact": 2},
    )


def _source_tree(tmp_path: Path) -> tuple[Path, dict[str, str], Cluster]:
    root = tmp_path / "source"
    for relative in a3g.SOURCE_PATHS.values():
        (root / relative).parent.mkdir(parents=True, exist_ok=True)
    cluster = _cluster()
    (root / a3g.SOURCE_PATHS["family_progress"]).write_text("{}\n")
    (root / a3g.SOURCE_PATHS["seed"]).write_text(cluster.to_xyz())
    (root / a3g.SOURCE_PATHS["child_log"]).write_text("failed owner H52:O15->O9\n")
    (root / a3g.SOURCE_PATHS["launch"]).write_text("launch\n")
    (root / a3g.SOURCE_PATHS["restoration"]).write_text("restored\n")
    (root / a3g.SOURCE_PATHS["family_receipt"]).write_text(
        json.dumps(
            {
                "schema": "a3-family-campaign-terminal-v1",
                "success": False,
                "family": "oaa",
                "state": "neutral",
                "current_n_intact": 2,
                "completed": [],
                "expected_git_sha": a3g.EXPECTED_EXECUTION_SOURCE,
                "observed_git_sha": a3g.EXPECTED_EXECUTION_SOURCE,
            },
            sort_keys=True,
        )
    )
    (root / a3g.SOURCE_PATHS["metadata"]).write_text(
        json.dumps(
            {
                "site_kind": "Oaa",
                "state": "neutral",
                "n_intact": 2,
                "charge": 0,
                "method": "b3lyp/def2-svp/df",
                "driver_git_commit": a3g.EXPECTED_EXECUTION_SOURCE,
                "n_atoms": 55,
                "center_site": 18,
                "metal_shells": 2,
            },
            sort_keys=True,
        )
    )
    hashes = {
        name: a3g.sha256_path(root / relative)
        for name, relative in a3g.SOURCE_PATHS.items()
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
    monkeypatch.setattr(a3g, "optimize_bounded", fake_optimize)
    monkeypatch.setattr(
        a3g, "energy", lambda cluster, _settings: energy_values[cluster.name]
    )
    monkeypatch.setattr(
        a3g, "gradient", lambda cluster, _settings: np.zeros_like(cluster.coords)
    )
    imaginary = np.asarray([] if imag is None else imag, dtype=float)
    monkeypatch.setattr(
        a3g,
        "frequencies",
        lambda cluster, settings: SimpleNamespace(
            imaginary_cm=imaginary,
            electronic_hartree=energy_values[cluster.name],
            geometry_fingerprint=a3g.frequency_geometry_fingerprint(cluster),
            settings_fingerprint=a3g.frequency_settings_fingerprint(settings),
        ),
    )
    monkeypatch.setattr(
        verifier, "gradient", lambda cluster, _settings: np.zeros_like(cluster.coords)
    )
    monkeypatch.setattr(verifier, "energy", a3g.energy)
    monkeypatch.setattr(verifier, "frequencies", a3g.frequencies)
    return calls


def _endpoints(source: Cluster, *, transfer: bool = False):
    values = []
    for index, name in enumerate(("conditioned", "constrained", "released"), 1):
        coords = source.coords.copy()
        coords[1, 1] = 0.01 * index
        values.append(replace(source, name=name, coords=coords))
    if transfer:
        coords = values[-1].coords.copy()
        coords[a3g.H52] = source.coords[a3g.O9] + np.array([0.96, 0.0, 0.0])
        values[-1] = replace(values[-1], coords=coords)
    return values


def test_validate_source_rehashes_full_contract_and_binds_seed(tmp_path, monkeypatch):
    root, hashes, cluster = _source_tree(tmp_path)
    monkeypatch.setattr(a3g, "canonical_template", lambda _root: cluster)
    source = a3g.validate_source(root, expected_hashes=hashes, repo_root=tmp_path)
    assert source.hashes == hashes
    assert a3g._owner_labels(source.cluster) == ["H52:O15"]
    assert geometry_hash(source.cluster.to_xyz()) == geometry_hash(cluster.to_xyz())

    (root / a3g.SOURCE_PATHS["child_log"]).write_text("tampered\n")
    with pytest.raises(RuntimeError, match="child_log SHA-256 mismatch"):
        a3g.validate_source(root, expected_hashes=hashes, repo_root=tmp_path)


def test_exact_three_stage_budget_and_fresh_release(tmp_path, monkeypatch):
    source = _source()
    endpoints = _endpoints(source.cluster)
    calls = _patch_calculators(monkeypatch, endpoints)
    terminal = a3g.run_experiment(
        tmp_path,
        source,
        code_revision="a" * 40,
    )
    assert terminal["classification"] == a3g.ACCEPTED_CANDIDATE
    assert len(calls) == 3
    expected_constraint = [(a3g.O15, a3g.H52, pytest.approx(0.96))]
    assert calls[0][2] == {"max_steps": 100, "fixed_distances": expected_constraint}
    assert calls[0][1].xc == "hf"
    assert calls[0][1].basis == "sto-3g"
    assert calls[1][2] == {"max_steps": 100, "fixed_distances": expected_constraint}
    assert calls[1][1].xc == "b3lyp"
    assert calls[1][1].density_fit is True
    assert calls[2][2] == {"max_steps": 100}
    for spec in a3g.STAGES:
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
    terminal = a3g.run_experiment(tmp_path, source, code_revision="b" * 40)
    assert len(calls) == 2
    assert terminal["classification"] == a3g.INCONCLUSIVE_CANDIDATE
    assert terminal["stages"]["released-production"]["status"] == "not-run"
    assert not (tmp_path / "stages" / a3g.STAGES[2].directory).exists()


def test_production_alternative_requires_exact_h52_transfer_and_downhill_energy(
    tmp_path, monkeypatch
):
    source = _source()
    endpoints = _endpoints(source.cluster, transfer=True)
    calls = _patch_calculators(
        monkeypatch,
        endpoints,
        energies={"conditioned": -9.0, "constrained": -10.0, "released": -10.01},
    )
    terminal = a3g.run_experiment(tmp_path, source, code_revision="c" * 40)
    assert len(calls) == 3
    assert terminal["classification"] == a3g.ALTERNATIVE_CANDIDATE
    released = json.loads(
        (tmp_path / "stages" / a3g.STAGES[2].directory / "receipt.json").read_text()
    )
    assert released["owner_changes"] == ["H52:O15->O9"]
    assert released["phva"]["status"] == "passed"


def test_alternative_strict_downhill_boundary_is_inconclusive():
    constrained = {
        "status": "complete",
        "optimizer": {"converged": True},
        "stationary": True,
        "owner_retaining": True,
        "owner_changes": [],
        "observed_owners": ["H52:O15"],
        "structure": {
            "constraint_residuals": [
                {
                    "oxygen": a3g.O15,
                    "hydrogen": a3g.H52,
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
        "owner_changes": ["H52:O15->O9"],
        "observed_owners": ["H52:O9"],
        "evidence": {
            "energy_hartree": -10.0 - a3g.DOWNHILL_MIN_HARTREE,
            "passed": True,
        },
        "phva": {"status": "not-required"},
    }
    assert a3g.classify(constrained, released) == a3g.INCONCLUSIVE_CANDIDATE


def test_raw_endpoint_and_pending_receipt_precede_structural_gate(
    tmp_path, monkeypatch
):
    source = _source()
    endpoints = _endpoints(source.cluster)
    _patch_calculators(monkeypatch, endpoints)
    real_gate = a3b.structural_gate
    observations = []

    def checking_gate(raw, reference, active, stage):
        spec = a3g.STAGE_BY_ID[stage]
        stage_dir = tmp_path / "stages" / spec.directory
        observations.append(
            (
                (stage_dir / "raw-endpoint.xyz").is_file(),
                json.loads((stage_dir / "receipt.json").read_text())["status"],
            )
        )
        return real_gate(raw, reference, active, stage)

    monkeypatch.setattr(a3b, "structural_gate", checking_gate)
    a3g.run_experiment(tmp_path, source, code_revision="d" * 40)
    assert observations == [(True, "pending-structural-gates")] * 6


def test_reservation_prevents_replay_after_orphaned_budget(tmp_path, monkeypatch):
    source = _source()
    stage_dir = tmp_path / "stages" / a3g.STAGES[0].directory
    stage_dir.mkdir(parents=True)
    (stage_dir / "reservation.json").write_text(json.dumps({"status": "reserved"}))
    monkeypatch.setattr(
        a3g,
        "optimize_bounded",
        lambda *_args, **_kwargs: pytest.fail("orphaned budget was replayed"),
    )
    terminal = a3g.run_experiment(tmp_path, source, code_revision="e" * 40)
    assert terminal["classification"] == a3g.INCONCLUSIVE_CANDIDATE
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
    a3g.run_experiment(tmp_path, source, code_revision="f" * 40)
    monkeypatch.setattr(
        verifier, "energy", lambda cluster, _settings: energy_values[cluster.name]
    )
    monkeypatch.setattr(
        verifier,
        "frequencies",
        lambda cluster, settings: SimpleNamespace(
            imaginary_cm=np.asarray([], dtype=float),
            electronic_hartree=energy_values[cluster.name],
            geometry_fingerprint=a3g.frequency_geometry_fingerprint(cluster),
            settings_fingerprint=a3g.frequency_settings_fingerprint(settings),
        ),
    )
    result = verifier.verify_experiment(
        tmp_path,
        source_override=source,
        verifier_identity="cold-worker-2",
    )
    assert result["status"] == "verified"
    assert result["classification"] == a3g.VERIFIED_ACCEPTED
    assert result["verifier_identity"] != result["executor_identity"]

    released = tmp_path / "stages" / a3g.STAGES[2].directory / "endpoint.xyz"
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
    a3g.run_experiment(tmp_path, source, code_revision="1" * 40)
    same = verifier.verify_experiment(
        tmp_path,
        source_override=source,
        verifier_identity=a3g.EXECUTOR_IDENTITY,
    )
    assert same["status"] == "rejected"
    a3g.atomic_json(
        tmp_path / "candidate-terminal.json",
        json.loads(
            next((tmp_path / "revoked").rglob("candidate-terminal.json")).read_text()
        ),
    )
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
    a3g.run_experiment(tmp_path, source, code_revision="g" * 40)
    monkeypatch.setattr(
        verifier, "energy", lambda cluster, _settings: energy_values[cluster.name]
    )
    monkeypatch.setattr(
        verifier,
        "frequencies",
        lambda cluster, settings: SimpleNamespace(
            imaginary_cm=np.asarray([], dtype=float),
            electronic_hartree=energy_values[cluster.name],
            geometry_fingerprint=a3g.frequency_geometry_fingerprint(cluster),
            settings_fingerprint=a3g.frequency_settings_fingerprint(settings),
        ),
    )

    spec = a3g.STAGES[-1]
    stage_dir = tmp_path / "stages" / spec.directory
    raw_path = stage_dir / "raw-endpoint.xyz"
    tampered = a3b.read_cluster(raw_path, source.cluster)
    coords = tampered.coords.copy()
    coords[0] = coords[0] + np.array([0.05, 0.0, 0.0])
    tampered = replace(tampered, coords=coords)
    a3g.atomic_xyz(raw_path, tampered)
    receipt_path = stage_dir / "receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["raw_endpoint"] = a3g._artifact(tampered, raw_path)
    a3g.atomic_json(receipt_path, receipt)
    candidate_path = tmp_path / "candidate-terminal.json"
    candidate = json.loads(candidate_path.read_text())
    candidate["stages"][spec.stage_id]["receipt_sha256"] = a3g.sha256_path(receipt_path)
    a3g.atomic_json(candidate_path, candidate)

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
    a3g.run_experiment(tmp_path, source, code_revision="2" * 40)
    forbidden = {"results.json", "store.sqlite", "ts.xyz", "barrier.json", "petra.toml"}
    assert not any(path.name in forbidden for path in tmp_path.rglob("*"))


def test_spent_stage_receipts_persist_one_call_zero_retry_100_step_budget(
    tmp_path, monkeypatch
):
    source = _source()
    endpoints = _endpoints(source.cluster)
    _patch_calculators(monkeypatch, endpoints)
    terminal = a3g.run_experiment(tmp_path, source, code_revision="3" * 40)
    assert terminal["experiment_budget"] == {
        "owner_conditioning": 1,
        "constrained_production": 1,
        "released_production": 1,
        "retries": 0,
    }
    for spec in a3g.STAGES:
        receipt = json.loads(
            (tmp_path / "stages" / spec.directory / "receipt.json").read_text()
        )
        assert receipt["optimizer"]["observed_calls"] == 1
        assert receipt["optimizer"]["observed_retries"] == 0
        assert receipt["optimizer"]["observed_max_steps"] == 100


def test_verifier_rejects_spent_stage_budget_that_is_not_one_call_zero_retry_100_steps(
    tmp_path, monkeypatch
):
    source = _source()
    endpoints = _endpoints(source.cluster)
    energy_values = {endpoint.name: -10.0 for endpoint in endpoints}
    _patch_calculators(monkeypatch, endpoints, energies=energy_values)
    a3g.run_experiment(tmp_path, source, code_revision="4" * 40)
    monkeypatch.setattr(
        verifier, "energy", lambda cluster, _settings: energy_values[cluster.name]
    )
    monkeypatch.setattr(
        verifier,
        "frequencies",
        lambda cluster, settings: SimpleNamespace(
            imaginary_cm=np.asarray([], dtype=float),
            electronic_hartree=energy_values[cluster.name],
            geometry_fingerprint=a3g.frequency_geometry_fingerprint(cluster),
            settings_fingerprint=a3g.frequency_settings_fingerprint(settings),
        ),
    )

    spec = a3g.STAGES[2]
    receipt_path = tmp_path / "stages" / spec.directory / "receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["optimizer"]["observed_calls"] = 2
    receipt["optimizer"]["observed_retries"] = 1
    receipt["optimizer"]["observed_max_steps"] = 200
    a3g.atomic_json(receipt_path, receipt)
    candidate_path = tmp_path / "candidate-terminal.json"
    candidate = json.loads(candidate_path.read_text())
    candidate["stages"][spec.stage_id]["receipt_sha256"] = a3g.sha256_path(receipt_path)
    candidate["experiment_budget"]["released_production"] = 2
    candidate["experiment_budget"]["retries"] = 1
    a3g.atomic_json(candidate_path, candidate)

    result = verifier.verify_experiment(
        tmp_path,
        source_override=source,
        verifier_identity="cold-worker-2",
    )
    assert result["status"] == "rejected"
    assert "budget" in result["detail"].casefold()


@pytest.fixture(autouse=True)
def forbid_real_calculators(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("a real electronic-structure calculator was called")

    from quarry import pipeline

    for module in (pipeline, a3g, verifier, a3b):
        for name in ("energy", "gradient", "frequencies", "optimize_bounded"):
            if hasattr(module, name):
                monkeypatch.setattr(module, name, forbidden)


@pytest.mark.parametrize("artifact", list(a3g.SOURCE_PATHS))
def test_every_source_hash_tamper_stops_before_calculator(
    tmp_path, monkeypatch, artifact
):
    root, hashes, cluster = _source_tree(tmp_path)
    for module in (a3g, verifier):
        monkeypatch.setattr(module, "canonical_template", lambda _root: cluster)
    path = root / a3g.SOURCE_PATHS[artifact]
    path.write_bytes(path.read_bytes() + b"tamper")
    for module in (a3g, verifier):
        with pytest.raises(RuntimeError, match="SHA-256 mismatch"):
            module.validate_source(root, expected_hashes=hashes)


def test_exact_seed_template_and_card_constants():
    template = a3g.canonical_template(Path(__file__).resolve().parents[2])
    independent = verifier.canonical_template(Path(__file__).resolve().parents[2])
    assert len(template.symbols) == 58
    assert template.charge == template.spin == 0
    assert "oaa-c18-s2-n2" in template.name
    assert len(template.frozen_indices) == 17
    assert np.array_equal(template.coords, independent.coords)
    assert a3g.oxygen_proton_owners(template)[52] == 15
    assert a3g.h52_constraint(template) == [
        (15, 52, np.linalg.norm(template.coords[52] - template.coords[15]))
    ]
    source = _source(template)
    for spec in a3g.STAGES:
        assert verifier.stage_signature(source, spec, "a" * 40) == a3g.stage_signature(
            source, spec, "a" * 40
        )
    assert len(a3g.EXPECTED_SOURCE_SHA256) == 7
    card = (
        Path(__file__).resolve().parents[2]
        / "docs/program/cards/A3g-oaa-neutral-n2-proton-microstate-stability.md"
    ).read_text()
    assert all(value in card for value in a3g.EXPECTED_SOURCE_SHA256.values())
    assert a3g.EXPECTED_EXECUTION_SOURCE in card


@pytest.mark.parametrize("transfer", [False, True])
@pytest.mark.parametrize("imag", [[], [50.0], [100.01]])
def test_both_positive_routes_require_fresh_phva(tmp_path, monkeypatch, transfer, imag):
    source = _source()
    endpoints = _endpoints(source.cluster, transfer=transfer)
    _patch_calculators(
        monkeypatch,
        endpoints,
        energies={"conditioned": -9.0, "constrained": -10.0, "released": -10.01},
        imag=imag,
    )
    terminal = a3g.run_experiment(tmp_path, source, code_revision="a" * 40)
    expected = (
        (a3g.ALTERNATIVE_CANDIDATE if transfer else a3g.ACCEPTED_CANDIDATE)
        if not imag or max(imag) <= a3g.NOISE_FLOOR_CM
        else a3g.INCONCLUSIVE_CANDIDATE
    )
    assert terminal["classification"] == expected
    result = verifier.verify_experiment(
        tmp_path, source_override=source, verifier_identity="independent"
    )
    assert result["status"] == "verified", result
    assert result["classification"] == (
        a3g.VERIFIED_ALTERNATIVE
        if expected == a3g.ALTERNATIVE_CANDIDATE
        else expected.removeprefix("unverified: ")
    )


@pytest.mark.parametrize(
    "fault",
    [
        "rms",
        "max",
        "residual",
        "owner",
        "shell",
        "collision",
        "energy",
        "gradient",
        "unconverged",
    ],
)
def test_constrained_gates_prevent_release(tmp_path, monkeypatch, fault):
    source = _source()
    endpoints = _endpoints(source.cluster)[:2]
    coords = endpoints[1].coords.copy()
    if fault == "residual":
        coords[52, 0] += 0.0002
    elif fault == "owner":
        coords[52] = coords[9] + [0.96, 0, 0]
    elif fault == "shell":
        coords[0, 0] += 0.0201
    elif fault == "collision":
        coords[2] = coords[1] + [0.59, 0, 0]
    endpoints[1] = replace(endpoints[1], coords=coords)
    calls = _patch_calculators(
        monkeypatch, endpoints, converged=[True, fault != "unconverged"]
    )

    def forces(cluster, settings):
        values = np.zeros_like(cluster.coords)
        if cluster.name == "constrained":
            if fault == "rms":
                values[:] = 0.00031
            elif fault == "max":
                values[1, 0] = 0.000451
            elif fault == "gradient":
                values[1, 0] = np.nan
        return values

    monkeypatch.setattr(a3g, "gradient", forces)
    if fault == "energy":
        monkeypatch.setattr(
            a3g, "energy", lambda c, s: np.nan if c.name == "constrained" else -10.0
        )
    result = a3g.run_experiment(tmp_path, source, code_revision="a" * 40)
    assert len(calls) == 2
    assert result["classification"] == a3g.INCONCLUSIVE_CANDIDATE
    assert result["stages"]["released-production"]["status"] == "not-run"
    receipt = json.loads(
        (tmp_path / "stages" / a3g.STAGES[1].directory / "receipt.json").read_text()
    )
    assert "raw_endpoint" in receipt
    if fault not in ("energy", "gradient"):
        assert len(receipt["raw_evidence"]["gradient_hartree_per_bohr"]) == 58


@pytest.mark.parametrize(
    "transfer", ["wrong-proton", "wrong-target", "extra-transfer", "exact"]
)
def test_only_exact_transfer_with_every_other_owner_unchanged(
    tmp_path, monkeypatch, transfer
):
    original = _cluster()
    symbols = original.symbols.copy()
    symbols[51] = "H"
    symbols[8] = "O"
    coords = original.coords.copy()
    coords[51] = coords[8] + [0.96, 0, 0]
    source = _source(replace(original, symbols=symbols, coords=coords))
    endpoints = _endpoints(
        source.cluster, transfer=transfer in ("extra-transfer", "exact")
    )
    coords = endpoints[2].coords.copy()
    if transfer == "wrong-target":
        coords[52] = coords[8] + [-0.96, 0, 0]
    if transfer in ("wrong-proton", "extra-transfer"):
        coords[51] = coords[9] + [-0.96, 0, 0]
    endpoints[2] = replace(endpoints[2], coords=coords)
    _patch_calculators(
        monkeypatch,
        endpoints,
        energies={"conditioned": -9.0, "constrained": -10.0, "released": -11.0},
    )
    terminal = a3g.run_experiment(tmp_path, source, code_revision="a" * 40)
    assert terminal["classification"] == (
        a3g.ALTERNATIVE_CANDIDATE if transfer == "exact" else a3g.INCONCLUSIVE_CANDIDATE
    )
    verified = verifier.verify_experiment(
        tmp_path, source_override=source, verifier_identity="cold"
    )
    assert verified["status"] == "verified", verified
    assert verified["classification"] == (
        a3g.VERIFIED_ALTERNATIVE
        if terminal["classification"] == a3g.ALTERNATIVE_CANDIDATE
        else terminal["classification"].removeprefix("unverified: ")
    )


@pytest.mark.parametrize("delta", [0.0, -0.01, 1e-6, 1.0001e-6])
def test_strict_predeclared_downhill_tolerance(tmp_path, monkeypatch, delta):
    source = _source()
    endpoints = _endpoints(source.cluster, transfer=True)
    _patch_calculators(
        monkeypatch,
        endpoints,
        energies={"conditioned": 0.0, "constrained": 0.0, "released": -delta},
    )
    terminal = a3g.run_experiment(tmp_path, source, code_revision="a" * 40)
    assert terminal["classification"] == (
        a3g.ALTERNATIVE_CANDIDATE if delta > 1e-6 else a3g.INCONCLUSIVE_CANDIDATE
    )


def test_raw_full_precision_and_calculators_persist_before_all_gates(
    tmp_path, monkeypatch
):
    source = _source()
    endpoints = _endpoints(source.cluster)
    coords = endpoints[0].coords.copy()
    coords[1, 1] = np.pi / 100
    endpoints[0] = replace(endpoints[0], coords=coords)
    _patch_calculators(monkeypatch, endpoints)
    real = a3b.structural_gate

    def gate(raw, reference, active, stage):
        directory = tmp_path / "stages" / a3g.STAGE_BY_ID[stage].directory
        receipt = json.loads((directory / "receipt.json").read_text())
        evidence = receipt["raw_evidence"]
        persisted = verifier.read_cluster(
            directory / "raw-endpoint.xyz", source.cluster
        )
        assert evidence["geometry_fingerprint"] == a3g.frequency_geometry_fingerprint(
            persisted
        )
        assert evidence["settings"] == receipt["signature"]["settings"]
        assert np.asarray(evidence["gradient_hartree_per_bohr"]).shape == (58, 3)
        assert np.isfinite(evidence["energy_hartree"])
        return real(raw, reference, active, stage)

    monkeypatch.setattr(a3b, "structural_gate", gate)
    terminal = a3g.run_experiment(tmp_path, source, code_revision="a" * 40)
    assert terminal["classification"] == a3g.ACCEPTED_CANDIDATE
    result = verifier.verify_experiment(
        tmp_path, source_override=source, verifier_identity="cold"
    )
    assert result["status"] == "verified", result


@pytest.mark.parametrize(
    "artifact", ["raw-endpoint.xyz", "receipt.json", "reservation.json"]
)
def test_every_orphan_artifact_blocks_stage_replay(tmp_path, artifact):
    source = _source()
    stage_dir = tmp_path / "stages" / a3g.STAGES[0].directory
    stage_dir.mkdir(parents=True)
    (stage_dir / artifact).write_text("orphan")
    result = a3g.run_experiment(tmp_path, source, code_revision="a" * 40)
    assert result["classification"] == a3g.INCONCLUSIVE_CANDIDATE
    assert result["experiment_budget"]["owner_conditioning"] == 0


def test_completed_experiment_cannot_be_replayed(tmp_path, monkeypatch):
    source = _source()
    calls = _patch_calculators(monkeypatch, _endpoints(source.cluster))
    a3g.run_experiment(tmp_path, source, code_revision="a" * 40)
    result = a3g.run_experiment(tmp_path, source, code_revision="a" * 40)
    assert len(calls) == 3
    assert result["classification"] == a3g.INCONCLUSIVE_CANDIDATE
    assert "zero retries" in result["detail"]


@pytest.mark.parametrize("entrypoint", ["executor", "verifier"])
def test_source_hash_failure_quarantines_both_stale_terminals(
    tmp_path, monkeypatch, entrypoint
):
    source = _source()
    _patch_calculators(monkeypatch, _endpoints(source.cluster))
    a3g.run_experiment(tmp_path, source, code_revision="a" * 40)
    old_candidate = (tmp_path / "candidate-terminal.json").read_bytes()
    old_verified = b'{"status":"verified","classification":"accepted reactant minimum"}'
    (tmp_path / "verified-terminal.json").write_bytes(old_verified)

    def fail(*args, **kwargs):
        raise RuntimeError("source family_progress SHA-256 mismatch")

    if entrypoint == "executor":
        monkeypatch.setattr(a3g, "validate_source", fail)
        monkeypatch.setattr(
            a3g.sys, "argv", ["a3g", "--output-root", str(tmp_path), "--dry-run"]
        )
        assert a3g.main() == 1
        assert not (tmp_path / "verified-terminal.json").exists()
    else:
        monkeypatch.setattr(verifier, "validate_source", fail)
        result = verifier.verify_experiment(tmp_path, verifier_identity="cold")
        assert result["status"] == "rejected"
        assert not (tmp_path / "candidate-terminal.json").exists()
    assert any(
        p.read_bytes() == old_candidate
        for p in (tmp_path / "revoked").rglob("candidate-terminal.json")
    )
    assert any(
        p.read_bytes() == old_verified
        for p in (tmp_path / "revoked").rglob("verified-terminal.json")
    )


def test_verifier_has_no_runtime_executor_or_shared_science_calls(
    tmp_path, monkeypatch
):
    source = _source()
    _patch_calculators(monkeypatch, _endpoints(source.cluster))
    a3g.run_experiment(tmp_path, source, code_revision="a" * 40)

    def fail(*args, **kwargs):
        pytest.fail("verifier delegated scientific verification")

    for module, names in (
        (
            a3g,
            [
                "classify",
                "constrained_release_allowed",
                "gradient_metrics",
                "oxygen_proton_owners",
                "stage_signature",
                "source_map",
                "run_stage",
            ],
        ),
        (
            a3b,
            [
                "structural_gate",
                "gradient_metrics",
                "_project_gradient",
                "read_cluster",
            ],
        ),
    ):
        for name in names:
            monkeypatch.setattr(module, name, fail)
    result = verifier.verify_experiment(
        tmp_path, source_override=source, verifier_identity="cold"
    )
    assert result["status"] == "verified", result
    tree = ast.parse(Path(verifier.__file__).read_text())
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and isinstance(node.func.value, ast.Name)
        ):
            assert node.func.value.id != "a3b"
            if node.func.value.id == "a3g":
                assert node.func.attr == "SourceEvidence"


@pytest.mark.parametrize("frozen", [[], [15], [15, 52]])
def test_independent_projection_matches_orthogonal_geometry(frozen):
    cluster = replace(_cluster(), frozen_indices=frozen)
    active = [(15, 52, 0.96)]
    raw = np.random.default_rng(42).normal(size=(58, 3)) * 0.001
    first = a3g.gradient_metrics(cluster, raw, active)
    second = verifier.independent_metrics(cluster, raw, active)
    assert np.allclose(
        first["projected_gradient_hartree_per_bohr"],
        second["projected_gradient_hartree_per_bohr"],
        atol=1e-15,
    )
    projected = np.array(second["projected_gradient_hartree_per_bohr"])
    assert np.all(projected[frozen] == 0)
    if 52 not in frozen:
        assert projected[52, 0] == pytest.approx(projected[15, 0], abs=1e-15)
    assert first["passed"] == second["passed"]


@pytest.mark.parametrize("quantity", ["energy", "gradient", "phva"])
def test_alternative_requires_independently_reproducible_evidence(
    tmp_path, monkeypatch, quantity
):
    source = _source()
    _patch_calculators(
        monkeypatch,
        _endpoints(source.cluster, transfer=True),
        energies={"conditioned": -9.0, "constrained": -10.0, "released": -11.0},
    )
    a3g.run_experiment(tmp_path, source, code_revision="a" * 40)
    if quantity == "energy":
        monkeypatch.setattr(verifier, "energy", lambda c, s: 999.0)
    elif quantity == "gradient":
        monkeypatch.setattr(verifier, "gradient", lambda c, s: np.ones_like(c.coords))
    else:
        monkeypatch.setattr(
            verifier,
            "frequencies",
            lambda c, s: SimpleNamespace(
                imaginary_cm=[999.0],
                electronic_hartree=-11.0,
                geometry_fingerprint=a3g.frequency_geometry_fingerprint(c),
                settings_fingerprint=a3g.frequency_settings_fingerprint(s),
            ),
        )
    result = verifier.verify_experiment(
        tmp_path, source_override=source, verifier_identity="cold"
    )
    assert result["status"] == "rejected"


@pytest.mark.parametrize(
    "field",
    [
        "raw_energy",
        "raw_gradient",
        "raw_settings",
        "projection",
        "metric",
        "owners",
        "phva_settings",
        "phva_modes",
        "seed",
        "parent",
        "fresh_hessian",
    ],
)
def test_coordinated_receipt_tampering_is_detected(tmp_path, monkeypatch, field):
    source = _source()
    _patch_calculators(monkeypatch, _endpoints(source.cluster))
    a3g.run_experiment(tmp_path, source, code_revision="a" * 40)
    stage = a3g.STAGES[-1]
    path = tmp_path / "stages" / stage.directory / "receipt.json"
    receipt = json.loads(path.read_text())
    if field == "raw_energy":
        receipt["raw_evidence"]["energy_hartree"] -= 1
    elif field == "raw_gradient":
        receipt["raw_evidence"]["gradient_hartree_per_bohr"][0][0] = 1
    elif field == "raw_settings":
        receipt["raw_evidence"]["settings"]["basis"] = "sto-3g"
    elif field == "projection":
        receipt["evidence"]["projected_gradient_hartree_per_bohr"][1][0] = 1
    elif field == "metric":
        receipt["evidence"]["projected_gradient_max_hartree_per_bohr"] = 1
    elif field == "owners":
        receipt["observed_owners"] = []
    elif field == "phva_settings":
        receipt["phva"]["settings_fingerprint"] = "tampered"
    elif field == "phva_modes":
        receipt["phva"]["imaginary_cm"] = [50.0]
    elif field == "seed":
        receipt["seed"]["geometry_fingerprint"] = "tampered"
    elif field == "parent":
        receipt["parent"]["receipt_sha256"] = "tampered"
    else:
        receipt["optimizer"]["geometric_default_fresh_hessian"] = False
    a3g.atomic_json(path, receipt)
    candidate_path = tmp_path / "candidate-terminal.json"
    candidate = json.loads(candidate_path.read_text())
    candidate["stages"][stage.stage_id]["receipt_sha256"] = a3g.sha256_path(path)
    a3g.atomic_json(candidate_path, candidate)
    result = verifier.verify_experiment(
        tmp_path, source_override=source, verifier_identity="cold"
    )
    assert result["status"] == "rejected", result


@pytest.mark.parametrize("failure", ["owner", "collision", "shell", "residual"])
def test_independently_verified_structural_failure(tmp_path, monkeypatch, failure):
    source = _source()
    endpoints = _endpoints(source.cluster)[:2]
    coords = endpoints[-1].coords.copy()
    if failure == "owner":
        coords[52] = (coords[9] + coords[15]) / 2
    elif failure == "collision":
        coords[2] = coords[1] + [0.5, 0, 0]
    elif failure == "shell":
        coords[0, 0] += 0.021
    else:
        coords[52, 0] += 0.001
    endpoints[-1] = replace(endpoints[-1], coords=coords)
    _patch_calculators(monkeypatch, endpoints)
    terminal = a3g.run_experiment(tmp_path, source, code_revision="a" * 40)
    assert terminal["classification"] == a3g.INCONCLUSIVE_CANDIDATE
    result = verifier.verify_experiment(
        tmp_path, source_override=source, verifier_identity="cold"
    )
    assert result["status"] == "verified", result
    assert result["classification"] == a3g.VERIFIED_INCONCLUSIVE
