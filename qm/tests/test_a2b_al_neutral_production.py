from __future__ import annotations

import argparse
import json
import sqlite3
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from quarry.clusters import Cluster
from quarry.pipeline import DftSettings, FrequencyResult
from quarry.store import Store
from scripts import a2b_al_neutral_production as a2b
from scripts import production_energetics as a2


def neutral(name: str, *, state: str) -> Cluster:
    if state == "reactant":
        si_x, attacker_x, transferred_h_x, transferred_h_y = 1.60, 4.20, 4.20, 0.95
    elif state == "intermediate":
        si_x, attacker_x, transferred_h_x, transferred_h_y = 1.60, 3.45, 0.90, 0.0
    elif state == "product":
        si_x, attacker_x, transferred_h_x, transferred_h_y = 3.10, 4.65, 0.90, 0.0
    else:
        si_x, attacker_x, transferred_h_x, transferred_h_y = 2.30, 4.05, 0.90, 0.0
    return Cluster(
        name=name,
        symbols=["O", "Si", "O", "H", "H", "O", "H"],
        coords=np.array(
            [
                [0.0, 0.0, 0.0],
                [si_x, 0.0, 0.0],
                [attacker_x, 0.0, 0.0],
                [transferred_h_x, transferred_h_y, 0.0],
                [attacker_x, -0.95, 0.0],
                [0.0, 2.8, 0.0],
                [-0.95, 0.0, 0.0],
            ],
            dtype=float,
        ),
        charge=-1,
        spin=0,
    )


def fixture_source(root: Path) -> tuple[a2b.SourceContract, dict[str, Cluster]]:
    clusters = {
        "reactant": neutral("al-neutral-complex", state="reactant"),
        "intermediate": neutral("al-neutral-intermediate", state="intermediate"),
        "transition-state": neutral("al-neutral-cleavage-ts", state="ts"),
        "product": neutral("al-neutral-product", state="product"),
    }
    root.mkdir(parents=True)
    paths = {
        "reactant": "complex.xyz",
        "intermediate": "intermediate.xyz",
        "transition-state": "ts.xyz",
        "product": "hydrolyzed_product.xyz",
    }
    for role, relative in paths.items():
        (root / relative).write_text(a2.exact_xyz(clusters[role]))
    results = {
        "reaction": "al-neutral",
        "method": a2b.SVP_METHOD,
        "mechanism": "sequential-associative-uphill-slide-no-resolved-addition-saddle",
        "addition_attempt": {"status": "no-saddle"},
        "steps": {
            "addition": {
                "accepted_first_order_saddle": False,
                "classification": "quasi-barrierless-uphill-slide",
            }
        },
        "profile": {"highest_profile_ts": "cleavage"},
    }
    (root / "results.json").write_text(json.dumps(results, sort_keys=True) + "\n")
    with Store(root / "store.sqlite") as store:
        for role in ("reactant", "intermediate", "transition-state"):
            cluster = clusters[role]
            structure_id = store.add_structure(
                cluster.name,
                cluster.formula,
                a2.exact_xyz(cluster),
                charge=cluster.charge,
                spin=cluster.spin,
            )
            job = store.add_job(
                structure_id,
                "freq",
                a2b.SVP_METHOD,
                "fixture",
            )
            store.set_job_status(job, "done")
            value = {
                "reactant": -10.0,
                "intermediate": -9.97,
                "transition-state": -9.95,
            }[role]
            store.add_result(job, "electronic", value, "hartree")
    artifact_sha = {
        relative: a2.sha256_path(root / relative)
        for relative in (
            "store.sqlite",
            "results.json",
            "complex.xyz",
            "intermediate.xyz",
            "ts.xyz",
            "hydrolyzed_product.xyz",
        )
    }
    contract = a2b.SourceContract(
        artifact_sha256=artifact_sha,
        role_ids=dict(a2b.ROLE_IDS),
        role_names=dict(a2b.ROLE_NAMES),
        expected_basins={
            role: a2.si_neutral_signature(clusters[role], 2)
            for role in ("reactant", "intermediate", "product")
        },
    )
    return contract, clusters


def fake_frequency(cluster: Cluster, *, transition_state: bool) -> FrequencyResult:
    return FrequencyResult(
        frequencies_cm=np.array([100.0, 200.0, 300.0]),
        imaginary_cm=np.array([84.0]) if transition_state else np.array([]),
        electronic_hartree=-9.95 if transition_state else -10.0,
        molar_mass_kg=0.018,
        rotational_temperatures_k=(1.0, 2.0, 3.0),
        linear=False,
        geometry_fingerprint=a2.frequency_geometry_fingerprint(cluster),
        settings_fingerprint="fixture-settings",
    )


def args(run_dir: Path, source_root: Path) -> argparse.Namespace:
    return argparse.Namespace(
        source_root=source_root,
        run_dir=run_dir,
        attacker_index=2,
        gpu=False,
        gpu_mem_gb=18.0,
        minimum_steps=5,
        saddle_steps=5,
        irc_steps=5,
        imaginary_floor=30.0,
        threads=2,
        nice=0,
        log=None,
    )


def test_source_contract_binds_hashes_jobs_topology_and_rejected_addition(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source"
    contract, clusters = fixture_source(source)
    evidence = a2b.load_source_evidence(
        source,
        attacker_index=2,
        contract=contract,
    )
    assert evidence.svp_energies == {
        "reactant": -10.0,
        "intermediate": -9.97,
        "transition-state": -9.95,
    }
    assert evidence.receipt["rejected_addition_saddle_used"] is False
    assert evidence.clusters["product"].symbols == clusters["reactant"].symbols

    (source / "results.json").write_text(
        (source / "results.json")
        .read_text()
        .replace(
            '"accepted_first_order_saddle": false',
            '"accepted_first_order_saddle": true',
        )
    )
    changed = replace(
        contract,
        artifact_sha256={
            **contract.artifact_sha256,
            "results.json": a2.sha256_path(source / "results.json"),
        },
    )
    with pytest.raises(ValueError, match="two-imaginary-mode addition candidate"):
        a2b.load_source_evidence(source, attacker_index=2, contract=changed)


def test_source_contract_rejects_byte_and_store_job_drift(tmp_path: Path) -> None:
    source = tmp_path / "source"
    contract, _ = fixture_source(source)
    with (source / "complex.xyz").open("a") as stream:
        stream.write("\n")
    with pytest.raises(ValueError, match="source artifact drift"):
        a2b.load_source_evidence(source, attacker_index=2, contract=contract)

    contract, _ = fixture_source(tmp_path / "source2")
    source2 = tmp_path / "source2"
    connection = sqlite3.connect(source2 / "store.sqlite")
    connection.execute("UPDATE jobs SET method = 'drifted' WHERE id = 1")
    connection.commit()
    connection.close()
    changed = replace(
        contract,
        artifact_sha256={
            **contract.artifact_sha256,
            "store.sqlite": a2.sha256_path(source2 / "store.sqlite"),
        },
    )
    with pytest.raises(ValueError, match="source job drift"):
        a2b.load_source_evidence(source2, attacker_index=2, contract=changed)


def test_full_irc_requires_exact_intermediate_and_product_typed_identity() -> None:
    intermediate = neutral("intermediate", state="intermediate")
    product = neutral("product", state="product")
    receipt = a2b.require_sequential_irc(
        (intermediate, product), intermediate, product, attacker_index=2
    )
    assert receipt["typed_identity"].startswith("basin+")
    with pytest.raises(RuntimeError, match="exact accepted intermediate"):
        a2b.require_sequential_irc(
            (intermediate, intermediate), intermediate, product, attacker_index=2
        )


def test_driver_is_resume_safe_and_publishes_complete_store(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_root = tmp_path / "source"
    contract, source_clusters = fixture_source(source_root)
    evidence = a2b.SourceEvidence(
        root=source_root,
        clusters=source_clusters,
        svp_energies={
            "reactant": -10.0,
            "intermediate": -9.97,
            "transition-state": -9.95,
        },
        receipt={"fixture": True},
        contract=contract,
    )
    monkeypatch.setattr(a2b, "load_source_evidence", lambda *_args, **_kwargs: evidence)
    monkeypatch.setattr(a2b, "assert_source_unchanged", lambda _source: None)
    monkeypatch.setattr(
        a2b, "require_role_identities", lambda *_args, **_kwargs: {"ok": True}
    )
    monkeypatch.setattr(
        a2b,
        "require_sequential_irc",
        lambda *_args, **_kwargs: {"connected": True},
    )

    calls = {"minimum": 0, "ts": 0, "frequency": 0, "irc": 0, "energy": 0}

    def minimum(cluster: Cluster, *_args, **_kwargs) -> Cluster:
        calls["minimum"] += 1
        return cluster

    def transition_state(cluster: Cluster, *_args, **_kwargs) -> Cluster:
        calls["ts"] += 1
        return cluster

    def frequency(
        path: Path,
        cluster: Cluster,
        _settings: DftSettings,
        *,
        finite_difference: bool,
    ) -> FrequencyResult:
        assert finite_difference is True
        if path.exists():
            return a2.frequency_from_payload(json.loads(path.read_text()))
        calls["frequency"] += 1
        result = fake_frequency(
            cluster, transition_state="transition-state" in path.name
        )
        payload = a2.frequency_payload(result)
        payload["hessian_method"] = "finite-difference-gradient"
        a2.atomic_json(path, payload)
        return result

    def irc(_cluster: Cluster, *_args, **_kwargs) -> tuple[Cluster, Cluster]:
        calls["irc"] += 1
        return source_clusters["intermediate"], source_clusters["product"]

    method_offset = {
        a2.PRODUCTION_METHOD: -20.0,
        a2.B3LYP_D4_METHOD: -30.0,
    }
    role_delta = {
        "reactant": 0.0,
        "intermediate": 0.03,
        "transition-state": 0.05,
        "product": -0.01,
    }

    def single_point(
        path: Path, cluster: Cluster, settings: DftSettings, method: str
    ) -> float:
        if path.exists():
            return a2.checkpoint_energy(path, cluster, settings, method)
        calls["energy"] += 1
        role = next(role for role in role_delta if role in path.name)
        value = method_offset[method] + role_delta[role]
        a2.atomic_json(
            path,
            {
                "method": method,
                "electronic_hartree": value,
                "geometry_fingerprint": a2.frequency_geometry_fingerprint(cluster),
                "settings_fingerprint": a2.frequency_settings_fingerprint(settings),
            },
        )
        return value

    monkeypatch.setattr(a2, "optimize_minimum", minimum)
    monkeypatch.setattr(a2b, "find_ts", transition_state)
    monkeypatch.setattr(a2, "checkpoint_frequency", frequency)
    monkeypatch.setattr(a2b, "full_irc", irc)
    monkeypatch.setattr(a2, "checkpoint_converged_energy", single_point)

    run_dir = tmp_path / "run"
    invocation = args(run_dir, source_root)
    assert a2b.execute_with_status(invocation) == 0
    first_calls = dict(calls)
    assert first_calls == {
        "minimum": 3,
        "ts": 1,
        "frequency": 4,
        "irc": 1,
        "energy": 8,
    }
    assert a2b.execute_with_status(invocation) == 0
    assert calls == first_calls

    result = json.loads((run_dir / "results.json").read_text())
    assert result["rejected_addition_saddle_used"] is False
    assert result["production_profile_maximum_role"] == "transition-state"
    assert set(result["electronic_barriers_kj_mol"]) == {
        a2b.SVP_METHOD,
        a2.R2SCAN3C_METHOD,
        a2.PRODUCTION_METHOD,
        a2.B3LYP_D4_METHOD,
    }
    assert (
        result["stationary_points"]["transition-state"]["significant_imaginary_count"]
        == 1
    )
    terminal = json.loads((run_dir / "terminal-receipt.json").read_text())
    assert terminal["outcome"] == "success"
    assert terminal["running_record_present"] is False
    assert (run_dir / "quarantine").is_dir()

    connection = sqlite3.connect(run_dir / "store.sqlite")
    try:
        assert connection.execute("PRAGMA integrity_check").fetchone()[0] == "ok"
        assert connection.execute("SELECT COUNT(*) FROM structures").fetchone()[0] == 7
        assert connection.execute("SELECT COUNT(*) FROM jobs").fetchone()[0] == 16
    finally:
        connection.close()


def test_converged_energy_uses_bounded_newton_and_reuses_receipt(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cluster = neutral("reactant", state="reactant")
    settings = DftSettings(
        xc="wb97m-v",
        basis="def2-tzvpd",
        solvent="smd",
        density_fit=True,
    )

    class FakeScf:
        converged = True
        max_cycle = 20
        newton_calls = 0
        kernel_calls = 0

        def newton(self):
            self.newton_calls += 1
            return self

        def kernel(self):
            self.kernel_calls += 1
            return -42.5

    fake = FakeScf()
    monkeypatch.setattr(a2.pipeline, "build_mol", lambda *_args: object())
    monkeypatch.setattr(a2.pipeline, "_make_scf", lambda *_args: fake)
    path = tmp_path / "energy.json"
    assert a2.checkpoint_converged_energy(
        path, cluster, settings, a2.PRODUCTION_METHOD
    ) == pytest.approx(-42.5)
    assert fake.newton_calls == 1
    assert fake.kernel_calls == 1
    assert fake.max_cycle == 100
    receipt = json.loads(path.read_text())
    assert receipt["convergence_route"] == "newton-first"
    assert receipt["converged"] is True
    assert a2.checkpoint_converged_energy(
        path, cluster, settings, a2.PRODUCTION_METHOD
    ) == pytest.approx(-42.5)
    assert fake.kernel_calls == 1


def test_failure_quarantines_stale_outputs_and_publishes_terminal_receipt(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "results.json").write_text("stale\n")
    (run_dir / "store.sqlite").write_text("stale\n")
    invocation = args(run_dir, tmp_path / "source")
    monkeypatch.setattr(
        a2b, "run", lambda _args: (_ for _ in ()).throw(RuntimeError("finite failure"))
    )
    with pytest.raises(RuntimeError, match="finite failure"):
        a2b.execute_with_status(invocation)
    status = json.loads((run_dir / "run_status.json").read_text())
    terminal = json.loads((run_dir / "terminal-receipt.json").read_text())
    assert status["outcome"] == "incomplete-computational-failure"
    assert terminal["outcome"] == "incomplete-computational-failure"
    assert not (run_dir / "results.json").exists()
    assert not (run_dir / "store.sqlite").exists()
    quarantined = list((run_dir / "quarantine").glob("*/results.json"))
    assert len(quarantined) == 1
