"""Cheap deterministic gates for the bounded TASK-168 intermediate sweep."""

from __future__ import annotations

import json
from dataclasses import replace

import numpy as np
import pytest

from quarry.clusters import aluminosilicate_dimer, hydrolysis_complex, water
from quarry.pipeline import FrequencyResult
from scripts import phase1_xiao_lasaga as phase1
from scripts import task168_intermediate_sweep as sweep


def associative_fixture():
    cluster = hydrolysis_complex(aluminosilicate_dimer(), water(), mode="flank")
    ow_index = len(aluminosilicate_dimer().symbols)
    coords = cluster.coords.copy()
    coords[ow_index] = coords[phase1.SI_INDEX] + np.array([0.0, 1.80, 0.0])
    coords[ow_index + 1] = coords[phase1.BR_INDEX] + np.array([0.0, 0.98, 0.0])
    coords[ow_index + 2] = coords[ow_index] + np.array([0.0, 0.96, 0.0])
    return replace(cluster, name="associative", coords=coords), ow_index


def test_rotation_preserves_axis_projection_and_radius():
    point = np.array([1.0, 0.0, 1.0])
    start = np.array([0.0, 0.0, 0.0])
    end = np.array([1.0, 0.0, 0.0])

    rotated = sweep.rotate_point_about_axis(point, start, end, 90.0)

    assert rotated == pytest.approx([1.0, -1.0, 0.0], abs=1e-12)
    assert np.linalg.norm(rotated - end) == pytest.approx(np.linalg.norm(point - end))


def test_conformer_set_is_finite_and_changes_only_hydrogens():
    intermediate, ow_index = associative_fixture()

    seeds = sweep.intermediate_conformer_seeds(intermediate, ow_index)

    assert [label for label, _ in seeds] == [
        "baseline",
        "terminal-oh-03-180",
        "terminal-oh-05-180",
        "terminal-oh-07-180",
        "terminal-oh-09-180",
        "terminal-oh-11-180",
        "terminal-oh-13-180",
        "all-terminal-oh-180",
        "bridge-proton--120",
        "bridge-proton-+120",
        "water-proton-180",
    ]
    heavy = [
        index for index, symbol in enumerate(intermediate.symbols) if symbol != "H"
    ]
    for _label, seed in seeds:
        assert seed.symbols == intermediate.symbols
        assert seed.charge == intermediate.charge
        assert seed.spin == intermediate.spin
        assert np.all(np.isfinite(seed.coords))
        assert seed.coords[heavy] == pytest.approx(intermediate.coords[heavy])


def test_sweep_selects_only_significantly_lower_associative_candidate(
    monkeypatch, tmp_path
):
    intermediate, ow_index = associative_fixture()
    reactant = replace(intermediate, name="reactant", coords=intermediate.coords.copy())
    reactant.coords[ow_index] = reactant.coords[phase1.SI_INDEX] + np.array(
        [0.0, 3.20, 0.0]
    )
    sweep.save_xyz(intermediate, tmp_path / "intermediate.xyz")
    sweep.save_xyz(reactant, tmp_path / "complex.xyz")

    monkeypatch.setattr(sweep, "optimize", lambda cluster, settings, max_steps: cluster)

    def fake_frequencies(cluster, settings):
        return FrequencyResult(
            frequencies_cm=np.array([300.0, 700.0]),
            imaginary_cm=np.array([]),
            electronic_hartree=-99.95,
            molar_mass_kg=0.1,
            rotational_temperatures_k=(1.0, 2.0, 3.0),
            linear=False,
        )

    monkeypatch.setattr(sweep, "frequencies", fake_frequencies)

    def fake_energy(cluster, settings):
        if cluster.name.endswith("terminal-oh-03-180"):
            return -99.96
        si_ow = np.linalg.norm(
            cluster.coords[phase1.SI_INDEX] - cluster.coords[ow_index]
        )
        return -99.95 if si_ow < 2.3 else -100.0

    monkeypatch.setattr(sweep, "energy", fake_energy)
    manifest = sweep.run_sweep(tmp_path, sweep.DftSettings(), max_steps=10)

    summary = manifest["summary"]
    assert summary["selected_label"] == "terminal-oh-03-180"
    assert summary["lower_intermediate_found"] is True
    assert summary["best_lowering_vs_baseline_kj"] == pytest.approx(
        0.01 * sweep.HARTREE_TO_KJ
    )
    assert (tmp_path / "intermediate.task168-selected.xyz").is_file()
    persisted = json.loads(
        (tmp_path / "task168-intermediate-sweep" / "manifest.json").read_text()
    )
    assert persisted["summary"] == summary
    assert persisted["candidates"]["terminal-oh-03-180"]["n_imaginary"] == 0
    assert persisted["candidates"]["baseline"]["n_imaginary"] == 0


def test_sweep_rejects_associative_stationary_points_with_imaginary_modes(
    monkeypatch, tmp_path
):
    intermediate, ow_index = associative_fixture()
    reactant = replace(intermediate, name="reactant", coords=intermediate.coords.copy())
    reactant.coords[ow_index] = reactant.coords[phase1.SI_INDEX] + np.array(
        [0.0, 3.20, 0.0]
    )
    sweep.save_xyz(intermediate, tmp_path / "intermediate.xyz")
    sweep.save_xyz(reactant, tmp_path / "complex.xyz")

    monkeypatch.setattr(sweep, "optimize", lambda cluster, settings, max_steps: cluster)

    def fake_energy(cluster, settings):
        if cluster.name.endswith("terminal-oh-03-180"):
            return -99.96
        si_ow = np.linalg.norm(
            cluster.coords[phase1.SI_INDEX] - cluster.coords[ow_index]
        )
        return -99.95 if si_ow < 2.3 else -100.0

    def fake_frequencies(cluster, settings):
        imaginary = (80.0,) if cluster.name.endswith("terminal-oh-03-180") else ()
        return FrequencyResult(
            frequencies_cm=np.array([300.0, 700.0]),
            imaginary_cm=np.array(imaginary),
            electronic_hartree=-99.95,
            molar_mass_kg=0.1,
            rotational_temperatures_k=(1.0, 2.0, 3.0),
            linear=False,
        )

    monkeypatch.setattr(sweep, "energy", fake_energy)
    monkeypatch.setattr(sweep, "frequencies", fake_frequencies)
    manifest = sweep.run_sweep(tmp_path, sweep.DftSettings(), max_steps=10)

    assert manifest["summary"]["selected_label"] == "baseline"
    assert manifest["summary"]["lower_intermediate_found"] is False
    saddle = manifest["candidates"]["terminal-oh-03-180"]
    assert saddle["status"] == "rejected"
    assert saddle["n_imaginary"] == 1
    assert "imaginary modes" in saddle["reason"]
