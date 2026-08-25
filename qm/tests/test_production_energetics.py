"""CPU-only contracts for the A2 production re-tiering driver."""

from __future__ import annotations

import argparse
import json
import sqlite3
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from quarry.clusters import Cluster, water
from quarry.pipeline import (
    DftSettings,
    FrequencyResult,
    energy,
    frequency_geometry_fingerprint,
)
from quarry.store import Store
from scripts import production_energetics as a2


def neutral_pair(name: str, *, product: bool = False) -> Cluster:
    if product:
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
                [4.8, 0.0, 0.0],
                [0.9, 0.0, 0.0],
                [5.4, 0.7, 0.0],
                [0.0, 5.0, 0.0],
                [0.9, 5.0, 0.0],
            ]
        )
    else:
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.6, 0.0, 0.0],
                [5.0, 0.0, 0.0],
                [5.9, 0.0, 0.0],
                [4.7, 0.9, 0.0],
                [0.0, 5.0, 0.0],
                [0.9, 5.0, 0.0],
            ]
        )
    return Cluster(name, ["O", "Si", "O", "H", "H", "O", "H"], coords)


def fake_frequency(cluster: Cluster, settings, *, transition_state: bool):
    return FrequencyResult(
        frequencies_cm=np.array([100.0, 200.0, 300.0]),
        imaginary_cm=np.array([120.0]) if transition_state else np.array([]),
        electronic_hartree=-9.9 if transition_state else -10.0,
        molar_mass_kg=0.1,
        rotational_temperatures_k=(1.0, 2.0, 3.0),
        linear=False,
        geometry_fingerprint=frequency_geometry_fingerprint(cluster),
        settings_fingerprint=a2.frequency_settings_fingerprint(settings),
    )


def source_store(path: Path) -> tuple[int, int, Cluster, Cluster]:
    reactant = neutral_pair("si-neutral-reactant")
    transition_state = neutral_pair("si-neutral-transition-state")
    transition_state.coords[2] = np.array([3.2, 0.0, 0.0])
    with Store(path) as store:
        reactant_id = store.add_structure(
            reactant.name, reactant.formula, reactant.to_xyz()
        )
        ts_id = store.add_structure(
            transition_state.name,
            transition_state.formula,
            transition_state.to_xyz(),
        )
    return reactant_id, ts_id, reactant, transition_state


def test_composite_and_production_settings_are_exact():
    r2scan3c, production, b3lyp_d4 = a2.settings(use_gpu=True)
    assert r2scan3c.composite == "r2scan3c"
    assert r2scan3c.basis == "def2-mtzvpp"
    assert production.xc == "wb97m-v"
    assert production.basis == "def2-tzvpd"
    assert production.solvent == "smd"
    assert production.dispersion is None  # VV10 is intrinsic; no double count.
    assert b3lyp_d4.dispersion == "d4"
    assert all(item.use_gpu for item in (r2scan3c, production, b3lyp_d4))


def test_ase_minimum_converges_on_cheap_water(tmp_path):
    initial = water()
    cheap = DftSettings(xc="hf", basis="sto-3g")
    optimized = a2.optimize_minimum(
        initial,
        cheap,
        max_steps=50,
        trajectory=tmp_path / "water.traj",
    )
    assert energy(optimized, cheap) < energy(initial, cheap)
    assert (tmp_path / "water.traj").exists()


def test_source_store_loader_checks_geometry_hash(tmp_path):
    path = tmp_path / "source.sqlite"
    reactant_id, _, reactant, _ = source_store(path)
    loaded, receipt = a2.load_store_structure(path, reactant_id)
    assert loaded.symbols == reactant.symbols
    assert np.array_equal(loaded.coords, reactant.coords)
    assert receipt["geometry_hash"]

    with sqlite3.connect(path) as connection:
        connection.execute(
            "UPDATE structures SET geometry_hash = 'bad' WHERE id = ?", (reactant_id,)
        )
        connection.commit()
    with pytest.raises(ValueError, match="hash mismatch"):
        a2.load_store_structure(path, reactant_id)


def test_cluster_checkpoint_is_exact_and_bound_to_source(tmp_path):
    source = neutral_pair("source")
    computed = replace(source, coords=source.coords.copy())
    computed.coords[1, 0] = 1.6000000010879999
    path = tmp_path / "minimum.xyz"
    identity = {"stage": "minimum", "settings": "exact"}

    first = a2.checkpoint_cluster(
        path,
        source,
        lambda: computed,
        identity=identity,
    )
    second = a2.checkpoint_cluster(
        path,
        source,
        lambda: pytest.fail("checkpoint recomputed"),
        identity=identity,
    )
    assert np.array_equal(first.coords, computed.coords)
    assert np.array_equal(second.coords, computed.coords)
    assert a2.frequency_geometry_fingerprint(
        first
    ) == a2.frequency_geometry_fingerprint(second)

    changed_source = replace(source, coords=source.coords.copy())
    changed_source.coords[0, 0] += 0.01
    with pytest.raises(ValueError, match="source geometry drift"):
        a2.checkpoint_cluster(
            path,
            changed_source,
            lambda: computed,
            identity=identity,
        )


def test_irc_checkpoint_binds_ts_settings_and_endpoints(tmp_path):
    transition_state = neutral_pair("transition-state")
    forward = replace(transition_state, coords=transition_state.coords.copy())
    reverse = replace(transition_state, coords=transition_state.coords.copy())
    forward.coords[2, 0] = 4.9
    reverse.coords[2, 0] = 1.7
    identity = {
        "stage": "full-irc",
        "algorithm": "sella-gonzalez-schlegel-full-irc-v1",
        "settings": "exact",
        "max_steps": 7,
    }

    first = a2.checkpoint_irc_endpoints(
        tmp_path,
        transition_state,
        lambda: (forward, reverse),
        identity=identity,
    )
    second = a2.checkpoint_irc_endpoints(
        tmp_path,
        transition_state,
        lambda: pytest.fail("IRC recomputed"),
        identity=identity,
    )
    assert np.array_equal(first[0].coords, forward.coords)
    assert np.array_equal(second[1].coords, reverse.coords)
    receipt = json.loads((tmp_path / "irc.r2scan3c.receipt.json").read_text())
    assert receipt["identity"] == identity
    assert receipt["directions"]["forward"]["path"] == "irc-forward.r2scan3c.xyz"
    assert receipt["directions"]["reverse"]["path"] == "irc-reverse.r2scan3c.xyz"

    changed_ts = replace(transition_state, coords=transition_state.coords.copy())
    changed_ts.coords[0, 0] += 0.01
    with pytest.raises(ValueError, match="transition-state geometry drift"):
        a2.checkpoint_irc_endpoints(
            tmp_path,
            changed_ts,
            lambda: (forward, reverse),
            identity=identity,
        )

    with pytest.raises(ValueError, match="identity drift"):
        a2.checkpoint_irc_endpoints(
            tmp_path,
            transition_state,
            lambda: (forward, reverse),
            identity={**identity, "max_steps": 8},
        )

    tampered = replace(forward, coords=forward.coords.copy())
    tampered.coords[1, 0] += 0.2
    (tmp_path / "irc-forward.r2scan3c.xyz").write_text(a2.exact_xyz(tampered))
    with pytest.raises(ValueError, match="forward geometry drift"):
        a2.checkpoint_irc_endpoints(
            tmp_path,
            transition_state,
            lambda: pytest.fail("IRC recomputed"),
            identity=identity,
        )

    (tmp_path / "irc-forward.r2scan3c.xyz").write_text(a2.exact_xyz(forward))
    (tmp_path / "irc.r2scan3c.receipt.json").unlink()
    with pytest.raises(ValueError, match="incomplete IRC checkpoint set"):
        a2.checkpoint_irc_endpoints(
            tmp_path,
            transition_state,
            lambda: pytest.fail("IRC recomputed"),
            identity=identity,
        )


def test_irc_checkpoint_cleans_partial_publish(monkeypatch, tmp_path):
    transition_state = neutral_pair("transition-state")
    forward = replace(transition_state, coords=transition_state.coords.copy())
    reverse = replace(transition_state, coords=transition_state.coords.copy())
    forward.coords[2, 0] = 4.9
    reverse.coords[2, 0] = 1.7
    identity = {
        "stage": "full-irc",
        "algorithm": "sella-gonzalez-schlegel-full-irc-v1",
        "settings": "exact",
        "max_steps": 7,
    }
    real_atomic = a2.atomic_json

    def boom(path, payload):
        raise OSError("simulated persist failure")

    monkeypatch.setattr(a2, "atomic_json", boom)
    with pytest.raises(OSError, match="simulated persist failure"):
        a2.checkpoint_irc_endpoints(
            tmp_path,
            transition_state,
            lambda: (forward, reverse),
            identity=identity,
        )
    assert not (tmp_path / "irc-forward.r2scan3c.xyz").exists()
    assert not (tmp_path / "irc-reverse.r2scan3c.xyz").exists()
    assert not (tmp_path / "irc.r2scan3c.receipt.json").exists()

    monkeypatch.setattr(a2, "atomic_json", real_atomic)
    loaded = a2.checkpoint_irc_endpoints(
        tmp_path,
        transition_state,
        lambda: (forward, reverse),
        identity=identity,
    )
    assert np.array_equal(loaded[0].coords, forward.coords)
    assert (tmp_path / "irc.r2scan3c.receipt.json").exists()


def test_full_irc_gate_requires_both_neutral_basins():
    reactant = neutral_pair("reactant")
    product = neutral_pair("product", product=True)
    receipt = a2.require_si_neutral_irc(
        (reactant, product), reactant, product, attacker_index=2
    )
    assert len(receipt["actual"]) == 2

    with pytest.raises(RuntimeError, match="do not match"):
        a2.require_si_neutral_irc(
            (reactant, reactant), reactant, product, attacker_index=2
        )

    wrong_ownership = replace(product, coords=product.coords.copy())
    wrong_ownership.coords[6] = np.array([5.0, -0.9, 0.0])
    assert a2.si_neutral_signature(wrong_ownership, 2) == a2.si_neutral_signature(
        product, 2
    )
    with pytest.raises(RuntimeError, match="exact H ownership"):
        a2.require_si_neutral_irc(
            (reactant, wrong_ownership),
            reactant,
            product,
            attacker_index=2,
        )


def test_checkpoint_frequency_rejects_geometry_drift(monkeypatch, tmp_path):
    cluster = neutral_pair("reactant")
    settings = a2.settings(use_gpu=False)[0]
    monkeypatch.setattr(
        a2,
        "frequencies",
        lambda current, current_settings: fake_frequency(
            current, current_settings, transition_state=False
        ),
    )
    path = tmp_path / "freq.json"
    a2.checkpoint_frequency(path, cluster, settings)
    shifted = Cluster(
        cluster.name,
        cluster.symbols,
        cluster.coords + np.array([1.0, 0.0, 0.0]),
    )
    with pytest.raises(ValueError, match="geometry fingerprint drift"):
        a2.checkpoint_frequency(path, shifted, settings)


def test_end_to_end_driver_journals_results_without_recompute(monkeypatch, tmp_path):
    source = tmp_path / "source.sqlite"
    reactant_id, ts_id, reactant, transition_state = source_store(source)
    product = neutral_pair("product", product=True)
    product_path = tmp_path / "product.xyz"
    product_path.write_text(product.to_xyz())

    monkeypatch.setattr(
        a2,
        "optimize_minimum",
        lambda current, settings, max_steps, trajectory: current,
    )
    ts_calls = []

    def find_transition_state(
        current, settings, max_steps, trajectory, initial_mode, internal
    ):
        ts_calls.append((initial_mode.copy(), internal))
        return current

    monkeypatch.setattr(a2, "find_ts", find_transition_state)

    frequency_calls: list[str] = []

    def frequencies(current, current_settings):
        frequency_calls.append(current.name)
        return fake_frequency(
            current,
            current_settings,
            transition_state="transition-state" in current.name,
        )

    monkeypatch.setattr(a2, "frequencies", frequencies)
    monkeypatch.setattr(a2, "frequencies_finite_difference", frequencies)
    monkeypatch.setattr(
        a2,
        "full_irc",
        lambda current, settings, **kwargs: (reactant, product),
    )

    def energy(current, current_settings):
        base = -20.0 if current_settings.xc == "wb97m-v" else -15.0
        return base + (0.05 if "transition-state" in current.name else 0.0)

    monkeypatch.setattr(a2, "energy", energy)
    run_dir = tmp_path / "run"
    args = argparse.Namespace(
        reaction="si-neutral",
        source_store=source,
        reactant_id=reactant_id,
        ts_id=ts_id,
        product_reference=product_path,
        attacker_index=2,
        run_dir=run_dir,
        gpu=False,
        minimum_steps=5,
        saddle_steps=5,
        irc_steps=5,
        imaginary_floor=30.0,
    )
    assert a2.execute_with_status(args) == 0
    result = json.loads((run_dir / "results.json").read_text())
    assert result["reaction"] == "si-neutral"
    assert result["full_irc"]["actual"] == result["full_irc"]["expected"]
    assert result["barrier_electronic_kj"][a2.PRODUCTION_METHOD] > 0.0
    assert (run_dir / "store.sqlite").exists()
    status = json.loads((run_dir / "run_status.json").read_text())
    assert status["status"] == "completed"
    assert status["results_sha256"] == a2.sha256_path(run_dir / "results.json")
    assert status["store_sha256"] == a2.sha256_path(run_dir / "store.sqlite")
    assert len(frequency_calls) == 3
    assert len(ts_calls) == 1
    assert ts_calls[0][0].shape == transition_state.coords.shape
    assert ts_calls[0][1] is False

    # The second run must consume exact checkpoints rather than invoke heavy work.
    monkeypatch.setattr(
        a2,
        "frequencies",
        lambda *unused, **kwargs: pytest.fail("frequency recomputed"),
    )
    monkeypatch.setattr(
        a2,
        "frequencies_finite_difference",
        lambda *unused, **kwargs: pytest.fail("frequency recomputed"),
    )
    monkeypatch.setattr(
        a2,
        "energy",
        lambda *unused, **kwargs: pytest.fail("single point recomputed"),
    )
    monkeypatch.setattr(
        a2,
        "full_irc",
        lambda *unused, **kwargs: pytest.fail("IRC recomputed"),
    )
    assert a2.execute_with_status(args) == 0
    irc_receipt = json.loads((run_dir / "irc.r2scan3c.receipt.json").read_text())
    assert irc_receipt["identity"]["algorithm"] == "sella-gonzalez-schlegel-full-irc-v1"
    assert irc_receipt["identity"]["max_steps"] == 5


def test_execute_with_status_records_failure(monkeypatch, tmp_path):
    args = argparse.Namespace(run_dir=tmp_path / "failed", reaction="si-neutral")

    def fail(current):
        raise RuntimeError("scientific gate red")

    monkeypatch.setattr(a2, "run", fail)
    with pytest.raises(RuntimeError, match="scientific gate red"):
        a2.execute_with_status(args)
    status = json.loads((args.run_dir / "run_status.json").read_text())
    assert status["status"] == "failed"
    assert status["error_type"] == "RuntimeError"
    assert status["error"] == "scientific gate red"
