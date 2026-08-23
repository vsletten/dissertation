"""CPU-only contracts for the final TASK-168 resolution runner."""

from __future__ import annotations

import json

import numpy as np
import pytest

from quarry.clusters import Cluster
from quarry.pipeline import DftSettings, FrequencyResult
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
