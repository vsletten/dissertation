"""CPU-only contracts for the A2 production re-tiering driver."""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path

import numpy as np
import pytest

from quarry.clusters import Cluster
from quarry.pipeline import FrequencyResult, frequency_geometry_fingerprint
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
            ]
        )
    return Cluster(name, ["O", "Si", "O", "H", "H"], coords)


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

    monkeypatch.setattr(a2, "optimize", lambda current, settings, max_steps: current)
    monkeypatch.setattr(
        a2,
        "find_ts",
        lambda current, settings, max_steps, trajectory: current,
    )

    frequency_calls: list[str] = []

    def frequencies(current, current_settings):
        frequency_calls.append(current.name)
        return fake_frequency(
            current,
            current_settings,
            transition_state="transition-state" in current.name,
        )

    monkeypatch.setattr(a2, "frequencies", frequencies)
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
    assert a2.run(args) == 0
    result = json.loads((run_dir / "results.json").read_text())
    assert result["reaction"] == "si-neutral"
    assert result["full_irc"]["actual"] == result["full_irc"]["expected"]
    assert result["barrier_electronic_kj"][a2.PRODUCTION_METHOD] > 0.0
    assert (run_dir / "store.sqlite").exists()
    assert len(frequency_calls) == 3

    # The second run must consume exact checkpoints rather than invoke heavy work.
    monkeypatch.setattr(
        a2,
        "frequencies",
        lambda *unused, **kwargs: pytest.fail("frequency recomputed"),
    )
    monkeypatch.setattr(
        a2,
        "energy",
        lambda *unused, **kwargs: pytest.fail("single point recomputed"),
    )
    assert a2.run(args) == 0
