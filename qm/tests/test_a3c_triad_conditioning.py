"""CPU-only orchestration and independent-verifier tests for A3c."""

from __future__ import annotations

import ast
import json
import subprocess
import threading
from dataclasses import asdict, replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from quarry.clusters import Cluster
from quarry.pipeline import DftSettings
from quarry.store import geometry_hash
from scripts import a3b_proton_microstate_stability as a3b
from scripts import a3c_triad_conditioning as a3c
from scripts import a3c_verify as verifier


def _cluster() -> Cluster:
    symbols = ["He"] * 75
    for index in (20, 26, 27, 29, 31, 32):
        symbols[index] = "O"
    for index in (50, 52, 57):
        symbols[index] = "H"
    coords = np.array([[float(index * 5), 0.0, 0.0] for index in range(75)])
    coords[50] = coords[26] + np.array([0.96, 0.0, 0.0])
    coords[52] = coords[27] + np.array([0.97, 0.0, 0.0])
    coords[57] = coords[29] + np.array([0.98, 0.0, 0.0])
    return Cluster(
        "conditioned",
        symbols,
        coords,
        charge=0,
        spin=0,
        frozen_indices=[0, 20, 26, 27, 29, 31, 32],
    )


def _source(cluster: Cluster | None = None) -> a3c.SourceEvidence:
    cluster = cluster or _cluster()
    a3a = a3b.SourceEvidence(
        cluster=cluster,
        source_root=Path("/fixture/a3a"),
        closeout_manifest_path=Path("/fixture/a3a/closeout-manifest.json"),
        closeout_manifest_sha256="1" * 64,
        conditioned_xyz_path=Path("/fixture/a3a/conditioned.xyz"),
        conditioned_xyz_sha256="2" * 64,
        conditioned_receipt_path=Path("/fixture/a3a/conditioned.json"),
        conditioned_receipt_sha256="3" * 64,
        conditioned_geometry_hash=geometry_hash(cluster.to_xyz()),
    )
    return a3c.SourceEvidence(
        a3a=a3a,
        a3b_root=Path("/fixture/a3b"),
        a3b_hashes={
            name: str(index) * 64
            for index, name in enumerate(a3c.A3B_RELATIVE_PATHS, 4)
        },
    )


def _write_a3b_tree(root: Path, source: a3c.SourceEvidence) -> dict[str, str]:
    for relative in a3c.A3B_RELATIVE_PATHS.values():
        (root / relative).parent.mkdir(parents=True, exist_ok=True)
    raw_path = root / a3c.A3B_RELATIVE_PATHS["raw_endpoint"]
    endpoint_path = root / a3c.A3B_RELATIVE_PATHS["endpoint"]
    raw_path.write_text(source.a3a.cluster.to_xyz())
    endpoint = replace(source.a3a.cluster, name="a3b-endpoint")
    endpoint_path.write_text(endpoint.to_xyz())
    raw_hash = a3b.sha256_path(raw_path)
    endpoint_hash = a3b.sha256_path(endpoint_path)
    expected_settings = asdict(
        DftSettings(xc="b3lyp", basis="def2-svp", density_fit=True, use_gpu=True)
    )
    receipt = {
        "schema": "a3b-stage-receipt-v1",
        "status": "complete",
        "stage": "common-dual",
        "source": a3c._source_map(source)["a3a"],
        "settings": expected_settings,
        "budget": {"max_steps": 100, "continuation_allowed": False},
        "optimizer": {
            "fresh_instance": True,
            "geometric_default_fresh_hessian": True,
        },
        "stationary": False,
        "owner_changes": ["H52:O27->O32"],
        "raw_endpoint": {"sha256": raw_hash},
        "endpoint": {"sha256": endpoint_hash},
    }
    receipt_path = root / a3c.A3B_RELATIVE_PATHS["stage_receipt"]
    receipt_path.write_text(json.dumps(receipt, sort_keys=True))
    receipt_hash = a3b.sha256_path(receipt_path)
    candidate = {
        "schema": "a3b-candidate-terminal-v1",
        "classification": a3b.INCONCLUSIVE,
        "source": a3c._source_map(source)["a3a"],
    }
    candidate_path = root / a3c.A3B_RELATIVE_PATHS["candidate_terminal"]
    candidate_path.write_text(json.dumps(candidate, sort_keys=True))
    candidate_hash = a3b.sha256_path(candidate_path)
    verified = {
        "schema": "a3b-verified-terminal-v1",
        "status": "verified",
        "classification": a3b.INCONCLUSIVE,
        "calculator_evidence_recomputed": True,
        "source_rehashed": True,
        "candidate_terminal_sha256": candidate_hash,
        "stage_receipt_sha256": {"common-dual": receipt_hash},
    }
    verified_path = root / a3c.A3B_RELATIVE_PATHS["verified_terminal"]
    verified_path.write_text(json.dumps(verified, sort_keys=True))
    return {
        name: a3b.sha256_path(root / relative)
        for name, relative in a3c.A3B_RELATIVE_PATHS.items()
    }


def _patch_calculators(monkeypatch, endpoint: Cluster, *, converged: bool = True):
    calls = []

    def fake_optimize(cluster, settings, **kwargs):
        calls.append((cluster, settings, kwargs))
        return SimpleNamespace(cluster=endpoint, converged=converged, max_steps=100)

    monkeypatch.setattr(a3c, "optimize_bounded", fake_optimize)
    monkeypatch.setattr(a3c, "energy", lambda *_args: -4022.0)
    monkeypatch.setattr(
        a3b, "gradient", lambda cluster, _settings: np.zeros_like(cluster.coords)
    )
    return calls


def test_validate_sources_rehashes_a3a_and_verified_a3b_chain(tmp_path, monkeypatch):
    source = _source()
    root = tmp_path / "a3b"
    hashes = _write_a3b_tree(root, source)
    monkeypatch.setattr(a3b, "validate_source", lambda *_args, **_kwargs: source.a3a)

    observed = a3c.validate_sources(
        Path("/fixture/a3a"),
        root,
        expected_a3b_sha256=hashes,
    )
    assert observed.a3b_hashes == hashes
    assert a3c._triad_retained(observed.a3a.cluster)

    endpoint_path = root / a3c.A3B_RELATIVE_PATHS["endpoint"]
    endpoint_path.write_text(
        endpoint_path.read_text().replace("a3b-endpoint", "tampered")
    )
    with pytest.raises(RuntimeError, match="A3b endpoint SHA-256 mismatch"):
        a3c.validate_sources(
            Path("/fixture/a3a"),
            root,
            expected_a3b_sha256=hashes,
        )


def test_single_fresh_budget_uses_exact_triad_constraints_and_binds_code(
    tmp_path, monkeypatch
):
    source = _source()
    endpoint = replace(source.a3a.cluster, name="endpoint")
    calls = _patch_calculators(monkeypatch, endpoint)

    terminal = a3c.run_experiment(
        tmp_path,
        source,
        use_gpu=False,
        code_revision="a" * 40,
    )

    assert terminal["classification"] == a3c.CANDIDATE_SEED
    assert len(calls) == 1
    assert calls[0][2] == {
        "max_steps": 100,
        "fixed_distances": [
            (26, 50, pytest.approx(0.96)),
            (27, 52, pytest.approx(0.97)),
            (29, 57, pytest.approx(0.98)),
        ],
    }
    receipt = json.loads((tmp_path / "receipt.json").read_text())
    assert receipt["budget"] == {
        "max_steps": 100,
        "continuation_allowed": False,
        "retry_allowed": False,
    }
    assert receipt["code_revision"] == "a" * 40
    assert receipt["optimizer"]["fresh_instance"] is True
    assert receipt["optimizer"]["geometric_default_fresh_hessian"] is True
    assert receipt["owner_retaining"] is True
    assert receipt["stationary"] is True


def test_raw_endpoint_and_pending_receipt_precede_structural_gate(
    tmp_path, monkeypatch
):
    source = _source()
    endpoint = replace(source.a3a.cluster, name="endpoint")
    _patch_calculators(monkeypatch, endpoint)
    real_gate = a3b.structural_gate
    observations = []

    def checking_gate(raw, reference, active_constraints, stage):
        observations.append(
            (
                (tmp_path / "raw-endpoint.xyz").is_file(),
                json.loads((tmp_path / "receipt.json").read_text())["status"],
            )
        )
        return real_gate(raw, reference, active_constraints, stage)

    monkeypatch.setattr(a3b, "structural_gate", checking_gate)
    a3c.run_experiment(tmp_path, source, use_gpu=False, code_revision="b" * 40)
    assert observations == [
        (True, "pending-structural-gates"),
        (True, "pending-structural-gates"),
    ]


def test_stale_signature_is_rejected_without_retry(tmp_path, monkeypatch):
    source = _source()
    endpoint = replace(source.a3a.cluster, name="endpoint")
    calls = _patch_calculators(monkeypatch, endpoint)
    a3c.run_experiment(tmp_path, source, use_gpu=False, code_revision="c" * 40)
    assert len(calls) == 1
    monkeypatch.setattr(
        a3c,
        "optimize_bounded",
        lambda *_args, **_kwargs: pytest.fail("stale run spent a second budget"),
    )
    terminal = a3c.run_experiment(
        tmp_path,
        source,
        use_gpu=True,
        code_revision="d" * 40,
    )
    assert terminal["classification"] == a3c.CANDIDATE_FAILURE
    assert terminal["receipt_status"] == "stale-evidence-rejected"


def test_concurrent_launch_cannot_spend_a_second_budget(tmp_path, monkeypatch):
    source = _source()
    endpoint = replace(source.a3a.cluster, name="endpoint")
    entered = threading.Event()
    release = threading.Event()
    calls = []

    def slow_optimize(cluster, settings, **kwargs):
        calls.append((cluster, settings, kwargs))
        entered.set()
        assert release.wait(timeout=5)
        return SimpleNamespace(cluster=endpoint, converged=True, max_steps=100)

    monkeypatch.setattr(a3c, "optimize_bounded", slow_optimize)
    monkeypatch.setattr(a3c, "energy", lambda *_args: -4022.0)
    monkeypatch.setattr(
        a3b,
        "gradient",
        lambda cluster, _settings: np.zeros_like(cluster.coords),
    )
    errors = []

    def first_run():
        try:
            a3c.run_experiment(
                tmp_path,
                source,
                use_gpu=True,
                code_revision="7" * 40,
            )
        except Exception as exc:  # pragma: no cover - surfaced by assertion below
            errors.append(exc)

    thread = threading.Thread(target=first_run)
    thread.start()
    assert entered.wait(timeout=5)
    try:
        with pytest.raises(RuntimeError, match="already active"):
            a3c.run_experiment(
                tmp_path,
                source,
                use_gpu=True,
                code_revision="7" * 40,
            )
    finally:
        release.set()
        thread.join(timeout=5)
    assert not thread.is_alive()
    assert errors == []
    assert len(calls) == 1


def test_resolve_code_revision_requires_clean_pushed_tree(tmp_path):
    repo = tmp_path / "repo"
    remote = tmp_path / "remote.git"

    def git(*args):
        return subprocess.run(
            ["git", *args],
            check=True,
            capture_output=True,
            text=True,
        )

    git("init", "--bare", str(remote))
    git("init", "-b", "main", str(repo))
    (repo / "tracked.py").write_text("VALUE = 1\n")
    git("-C", str(repo), "add", "tracked.py")
    git(
        "-C",
        str(repo),
        "-c",
        "user.name=Test",
        "-c",
        "user.email=test@example.invalid",
        "commit",
        "-m",
        "fixture",
    )
    git("-C", str(repo), "remote", "add", "origin", str(remote))
    git("-C", str(repo), "push", "-u", "origin", "main")
    assert len(a3c.resolve_code_revision(repo)) == 40
    (repo / "untracked.py").write_text("VALUE = 2\n")
    with pytest.raises(RuntimeError, match="fully clean"):
        a3c.resolve_code_revision(repo)


def test_nonconverged_or_transferred_endpoint_fails_closed(tmp_path, monkeypatch):
    source = _source()
    transferred_coords = source.a3a.cluster.coords.copy()
    transferred_coords[52] = source.a3a.cluster.coords[32] + np.array([0.96, 0.0, 0.0])
    endpoint = replace(
        source.a3a.cluster, name="transferred", coords=transferred_coords
    )
    _patch_calculators(monkeypatch, endpoint, converged=False)
    terminal = a3c.run_experiment(
        tmp_path, source, use_gpu=False, code_revision="e" * 40
    )
    assert terminal["classification"] == a3c.CANDIDATE_FAILURE
    receipt = json.loads((tmp_path / "receipt.json").read_text())
    assert receipt["status"] in {"complete", "gate-rejected"}
    assert receipt["stationary"] is False
    assert (tmp_path / "raw-endpoint.xyz").is_file()


def test_independent_verifier_accepts_and_detects_tampering(tmp_path, monkeypatch):
    tree = ast.parse(Path(verifier.__file__).read_text())
    imported = {
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, (ast.Import, ast.ImportFrom))
        for alias in node.names
    }
    assert not any("optim" in name for name in imported)

    source = _source()
    endpoint = replace(source.a3a.cluster, name="endpoint")
    _patch_calculators(monkeypatch, endpoint)
    a3c.run_experiment(tmp_path, source, use_gpu=True, code_revision="f" * 40)
    result = verifier.verify_experiment(
        tmp_path,
        source_override=source,
        recompute_calculators=False,
    )
    assert result["status"] == "verified"
    assert result["classification"] == a3c.VERIFIED_SEED
    assert result["verifier_identity"] != result["executor_identity"]

    endpoint_path = tmp_path / "endpoint.xyz"
    endpoint_path.write_text(endpoint_path.read_text().replace("endpoint", "tampered"))
    rejected = verifier.verify_experiment(
        tmp_path,
        source_override=source,
        recompute_calculators=False,
    )
    assert rejected["status"] == "rejected"
    assert "hash-mismatched" in rejected["detail"]


def test_verifier_rejects_mutually_consistent_wrong_method(tmp_path, monkeypatch):
    source = _source()
    endpoint = replace(source.a3a.cluster, name="endpoint")
    _patch_calculators(monkeypatch, endpoint)
    a3c.run_experiment(tmp_path, source, use_gpu=True, code_revision="8" * 40)

    reservation_path = tmp_path / "reservation.json"
    receipt_path = tmp_path / "receipt.json"
    candidate_path = tmp_path / "candidate-terminal.json"
    reservation = json.loads(reservation_path.read_text())
    receipt = json.loads(receipt_path.read_text())
    candidate = json.loads(candidate_path.read_text())
    for payload in (reservation["signature"], receipt["signature"]):
        payload["settings"]["xc"] = "hf"
    receipt["settings"]["xc"] = "hf"
    candidate["settings"]["xc"] = "hf"
    reservation_path.write_text(json.dumps(reservation))
    receipt_path.write_text(json.dumps(receipt))
    candidate["receipt_sha256"] = a3b.sha256_path(receipt_path)
    candidate_path.write_text(json.dumps(candidate))

    result = verifier.verify_experiment(
        tmp_path,
        source_override=source,
        recompute_calculators=False,
    )
    assert result["status"] == "rejected"
    assert "required B3LYP" in result["detail"]


def test_stale_rerun_revokes_prior_verified_terminal(tmp_path, monkeypatch):
    source = _source()
    endpoint = replace(source.a3a.cluster, name="endpoint")
    _patch_calculators(monkeypatch, endpoint)
    a3c.run_experiment(tmp_path, source, use_gpu=True, code_revision="9" * 40)
    result = verifier.verify_experiment(
        tmp_path,
        source_override=source,
        recompute_calculators=False,
    )
    assert result["status"] == "verified"
    assert (tmp_path / "verified-terminal.json").is_file()

    terminal = a3c.run_experiment(
        tmp_path,
        source,
        use_gpu=True,
        code_revision="a" * 40,
    )
    assert terminal["classification"] == a3c.CANDIDATE_FAILURE
    assert not (tmp_path / "verified-terminal.json").exists()
    assert list((tmp_path / "revoked").rglob("verified-terminal.json"))


def test_no_downstream_scientific_outputs_are_emitted(tmp_path, monkeypatch):
    source = _source()
    endpoint = replace(source.a3a.cluster, name="endpoint")
    _patch_calculators(monkeypatch, endpoint)
    a3c.run_experiment(tmp_path, source, use_gpu=False, code_revision="1" * 40)
    forbidden = {"results.json", "store.sqlite", "ts.xyz", "barrier.json", "petra.toml"}
    assert not any(path.name in forbidden for path in tmp_path.rglob("*"))
