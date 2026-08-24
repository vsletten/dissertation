"""CPU-only contracts for the final TASK-168 resolution runner."""

from __future__ import annotations

import json
from dataclasses import replace

import numpy as np
import pytest

from quarry.clusters import Cluster
from quarry.pipeline import DftSettings, FrequencyResult
from quarry.rates import Thermo
from scripts import task168_resolution as resolution


def cluster(name: str, offset: float) -> Cluster:
    return Cluster(
        name,
        ["H", "H"],
        np.array([[offset, 0.0, 0.0], [offset + 0.8, 0.0, 0.0]]),
    )


def test_piecewise_band_has_exact_endpoints_and_densifies_toward_intermediate():
    reactant = cluster("r", 0.0)
    crest = cluster("crest", 3.0)
    intermediate = cluster("i", 4.0)

    images = resolution.densified_piecewise_band(reactant, crest, intermediate)

    assert len(images) == 9
    assert images[0] == pytest.approx(reactant.coords)
    assert images[-1] == pytest.approx(intermediate.coords)
    final_steps = [
        float(np.linalg.norm(images[index + 1] - images[index]))
        for index in range(4, 8)
    ]
    assert final_steps == sorted(final_steps, reverse=True)


def test_frequency_cache_round_trip_is_geometry_and_settings_bound(tmp_path):
    geometry = cluster("minimum", 0.0)
    settings = DftSettings(xc="hf", basis="sto-3g")
    result = FrequencyResult(
        frequencies_cm=np.array([100.0, 200.0]),
        imaginary_cm=np.array([]),
        electronic_hartree=-1.0,
        molar_mass_kg=0.002,
        rotational_temperatures_k=(1.0, 2.0, 3.0),
        linear=False,
        geometry_fingerprint=resolution.frequency_geometry_fingerprint(geometry),
        settings_fingerprint=resolution.frequency_settings_fingerprint(settings),
    )
    path = tmp_path / "frequency.npz"

    resolution.save_frequency(path, result)

    loaded = resolution.load_frequency(path, geometry, settings)
    assert loaded is not None
    assert loaded.electronic_hartree == -1.0
    assert loaded.n_imaginary == 0
    assert resolution.load_frequency(path, cluster("moved", 0.1), settings) is None
    assert (
        resolution.load_frequency(
            path, geometry, DftSettings(xc="b3lyp", basis="sto-3g")
        )
        is None
    )


def test_newest_pre_relaxed_checkpoint_ignores_failed_and_partial(tmp_path):
    old = tmp_path / "old"
    old.mkdir()
    (old / "manifest.json").write_text(
        json.dumps({"stage": "pre-relax-final", "converged": True})
    )
    failed = tmp_path / "failed"
    failed.mkdir()
    (failed / "manifest.json").write_text(
        json.dumps({"stage": "pre-relax-final", "converged": False})
    )
    partial = tmp_path / "partial"
    partial.mkdir()
    (partial / "manifest.json").write_text(
        json.dumps({"stage": "pre-relax-step-000005", "converged": None})
    )

    assert resolution.newest_pre_relaxed_checkpoint(tmp_path) == old


def test_mechanism_finding_publishes_only_with_verified_cleavage(monkeypatch, tmp_path):
    template = resolution.hydrolysis_complex(
        resolution.aluminosilicate_dimer(), resolution.water(), mode="flank"
    )
    ow_index = len(resolution.aluminosilicate_dimer().symbols)
    reactant_coords = template.coords.copy()
    reactant_coords[ow_index] = reactant_coords[resolution.SI_INDEX] + np.array(
        [0.0, 3.2, 0.0]
    )
    reactant_coords[ow_index + 1] = reactant_coords[ow_index] + np.array(
        [0.0, 0.96, 0.0]
    )
    reactant_coords[ow_index + 2] = reactant_coords[ow_index] + np.array(
        [0.0, -0.96, 0.0]
    )
    reactant = replace(template, name="reactant", coords=reactant_coords)
    intermediate_coords = template.coords.copy()
    intermediate_coords[ow_index] = intermediate_coords[resolution.SI_INDEX] + np.array(
        [0.0, 1.8, 0.0]
    )
    intermediate_coords[ow_index + 1] = intermediate_coords[
        resolution.BR_INDEX
    ] + np.array([0.0, 0.98, 0.0])
    intermediate_coords[ow_index + 2] = intermediate_coords[ow_index] + np.array(
        [0.0, 0.96, 0.0]
    )
    intermediate = replace(template, name="intermediate", coords=intermediate_coords)
    product_coords = intermediate.coords.copy()
    product_coords[resolution.BR_INDEX] = product_coords[
        resolution.SI_INDEX
    ] + np.array([3.5, 0.0, 0.0])
    product_coords[ow_index + 1] = product_coords[resolution.BR_INDEX] + np.array(
        [0.0, 0.98, 0.0]
    )
    product = replace(template, name="product", coords=product_coords)
    cleavage_ts = replace(intermediate, name="cleavage-ts")

    resolution.save_xyz(reactant, tmp_path / "complex.xyz")
    resolution.save_xyz(intermediate, tmp_path / "intermediate.task168-selected.xyz")
    (tmp_path / "task168-addition-lower-i-attempt.json").write_text(
        json.dumps({"status": "no-saddle", "reason": "two soft modes"})
    )
    sweep_dir = tmp_path / "task168-intermediate-sweep"
    sweep_dir.mkdir()
    (sweep_dir / "manifest.json").write_text(
        json.dumps({"summary": {"selected_label": "lower"}})
    )

    def frequency(electronic, imaginary=()):
        return FrequencyResult(
            frequencies_cm=np.array([300.0, 700.0]),
            imaginary_cm=np.array(imaginary),
            electronic_hartree=electronic,
            molar_mass_kg=0.1,
            rotational_temperatures_k=(1.0, 2.0, 3.0),
            linear=False,
        )

    reactant_freq = frequency(-100.0)
    intermediate_freq = frequency(-99.96)
    cleavage_freq = frequency(-99.95, (84.0,))

    def fake_cached(path, cluster, settings):
        return reactant_freq if path.name.startswith("complex") else intermediate_freq

    monkeypatch.setattr(resolution, "cached_frequency", fake_cached)
    monkeypatch.setattr(
        resolution,
        "thermo_result",
        lambda freq, temperature: Thermo(
            {-100.0: 0.0, -99.96: 100.0, -99.95: 130.0}[freq.electronic_hartree],
            0.0,
            0.0,
            0.0,
            temperature,
        ),
    )

    class FakeStore:
        structures = []

        def __init__(self, path):
            self.path = path
            path.touch()

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def add_structure(self, name, *args, **kwargs):
            self.structures.append(name)
            return len(self.structures)

        def add_job(self, *args, **kwargs):
            return 1

        def set_job_status(self, *args, **kwargs):
            return None

        def add_result(self, *args, **kwargs):
            return None

    monkeypatch.setattr(resolution, "Store", FakeStore)
    results = resolution.finalize(
        tmp_path,
        DftSettings(),
        298.15,
        addition=None,
        cleavage=(cleavage_ts, cleavage_freq, intermediate, product),
    )

    assert results["mechanism_finding"] is True
    assert results["mechanism"].endswith("no-resolved-addition-saddle")
    assert results["profile"]["overall_profile_dG_dagger_kj"] == 130.0
    assert results["steps"]["addition"]["accepted_first_order_saddle"] is False
    assert FakeStore.structures == [
        "al-neutral-complex",
        "al-neutral-intermediate",
        "al-neutral-cleavage-ts",
    ]
    assert json.loads((tmp_path / "run_status.json").read_text())["status"] == (
        "completed"
    )


def _mechanism_finding_closeout_inputs(tmp_path):
    template = resolution.hydrolysis_complex(
        resolution.aluminosilicate_dimer(), resolution.water(), mode="flank"
    )
    ow_index = len(resolution.aluminosilicate_dimer().symbols)
    reactant_coords = template.coords.copy()
    reactant_coords[ow_index] = reactant_coords[resolution.SI_INDEX] + np.array(
        [0.0, 3.2, 0.0]
    )
    reactant_coords[ow_index + 1] = reactant_coords[ow_index] + np.array(
        [0.0, 0.96, 0.0]
    )
    reactant_coords[ow_index + 2] = reactant_coords[ow_index] + np.array(
        [0.0, -0.96, 0.0]
    )
    reactant = replace(template, name="reactant", coords=reactant_coords)
    intermediate_coords = template.coords.copy()
    intermediate_coords[ow_index] = intermediate_coords[resolution.SI_INDEX] + np.array(
        [0.0, 1.8, 0.0]
    )
    intermediate_coords[ow_index + 1] = intermediate_coords[
        resolution.BR_INDEX
    ] + np.array([0.0, 0.98, 0.0])
    intermediate_coords[ow_index + 2] = intermediate_coords[ow_index] + np.array(
        [0.0, 0.96, 0.0]
    )
    intermediate = replace(template, name="intermediate", coords=intermediate_coords)
    product_coords = intermediate.coords.copy()
    product_coords[resolution.BR_INDEX] = product_coords[
        resolution.SI_INDEX
    ] + np.array([3.5, 0.0, 0.0])
    product_coords[ow_index + 1] = product_coords[resolution.BR_INDEX] + np.array(
        [0.0, 0.98, 0.0]
    )
    product = replace(template, name="product", coords=product_coords)
    cleavage_ts = replace(intermediate, name="cleavage-ts")

    resolution.save_xyz(reactant, tmp_path / "complex.xyz")
    resolution.save_xyz(intermediate, tmp_path / "intermediate.task168-selected.xyz")
    (tmp_path / "task168-addition-lower-i-attempt.json").write_text(
        json.dumps({"status": "no-saddle", "reason": "two soft modes"})
    )
    sweep_dir = tmp_path / "task168-intermediate-sweep"
    sweep_dir.mkdir()
    (sweep_dir / "manifest.json").write_text(
        json.dumps({"summary": {"selected_label": "lower"}})
    )
    cleavage_freq = FrequencyResult(
        frequencies_cm=np.array([300.0, 700.0]),
        imaginary_cm=np.array((84.0,)),
        electronic_hartree=-99.95,
        molar_mass_kg=0.1,
        rotational_temperatures_k=(1.0, 2.0, 3.0),
        linear=False,
    )
    return cleavage_ts, cleavage_freq, intermediate, product


def test_finalize_marks_failed_when_closeout_after_begin_raises(monkeypatch, tmp_path):
    cleavage = _mechanism_finding_closeout_inputs(tmp_path)

    def boom(path, cluster, settings):
        raise RuntimeError("frequency cache exploded")

    monkeypatch.setattr(resolution, "cached_frequency", boom)

    with pytest.raises(RuntimeError, match="frequency cache exploded"):
        resolution.finalize(
            tmp_path,
            DftSettings(),
            298.15,
            addition=None,
            cleavage=cleavage,
        )

    status = json.loads((tmp_path / "run_status.json").read_text())
    assert status["status"] == "failed"
    assert "frequency cache exploded" in status["detail"]
