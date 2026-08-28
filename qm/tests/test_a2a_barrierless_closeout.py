"""CPU-only contracts for the A2a barrierless production closeout."""

from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pytest

from quarry.clusters import Cluster
from quarry.pipeline import FrequencyResult
from scripts import a2a_barrierless_closeout as closeout


def cluster(name: str) -> Cluster:
    return Cluster(
        name,
        ["H", "H"],
        np.array([[0.0, 0.0, 0.0], [0.8, 0.0, 0.0]]),
    )


def frequency(electronic_hartree: float) -> FrequencyResult:
    return FrequencyResult(
        frequencies_cm=np.array([100.0]),
        imaginary_cm=np.array([]),
        electronic_hartree=electronic_hartree,
        molar_mass_kg=0.1,
        rotational_temperatures_k=(1.0, 2.0, 3.0),
        linear=False,
        geometry_fingerprint="geometry",
        settings_fingerprint="settings",
    )


def write_classification(root, *, outcome="barrierless-shelf"):
    classification = root / "cleavage/coupled-scan-v1/complete-grid-classification.json"
    release = root / "cleavage/barrierless-downhill-release-v2/release-receipt.json"
    classification.parent.mkdir(parents=True, exist_ok=True)
    release.parent.mkdir(parents=True, exist_ok=True)
    classification.write_text(
        json.dumps(
            {
                "cell_count": 81,
                "classification": {
                    "outcome": outcome,
                    "barrier_threshold_kj_mol": 2.0,
                    "product_cell": [7, 7],
                    "candidate_kind": (
                        "complete-grid-valley-no-interior-crest-above-threshold"
                    ),
                },
            }
        )
    )
    release.write_text(
        json.dumps(
            {
                "release_version": "a2a-barrierless-downhill-release-v2",
                "status": "completed",
                "seed_cell": [7, 7],
                "zero_index_minimum": True,
                "typed_product_identity_matches": True,
                "imaginary_mode_count": 0,
                "electronic_delta_kj_mol": -1.0,
                "cartesian_displacement_a": 0.2,
                "minimum_cartesian_displacement_a": 0.001,
            }
        )
    )


def test_barrierless_classification_and_release_are_fail_closed(tmp_path):
    write_classification(tmp_path)
    classification, release = closeout._classification_and_release(tmp_path)
    assert classification["classification"]["outcome"] == "barrierless-shelf"
    assert release["zero_index_minimum"] is True

    write_classification(tmp_path, outcome="interior-saddle")
    with pytest.raises(RuntimeError, match="not barrierless-shelf"):
        closeout._classification_and_release(tmp_path)


def test_closeout_method_matrix_is_exact():
    settings = closeout.method_settings(use_gpu=True)
    assert set(settings) == {
        closeout.SVP_METHOD,
        closeout.production.PRODUCTION_METHOD,
        closeout.production.B3LYP_D4_METHOD,
    }
    assert settings[closeout.SVP_METHOD].basis == "def2-svp"
    assert settings[closeout.production.PRODUCTION_METHOD].basis == "def2-tzvpd"
    assert settings[closeout.production.B3LYP_D4_METHOD].dispersion == "d4"
    assert all(value.use_gpu for value in settings.values())


def test_dft_runs_four_by_three_matrix_and_uses_only_addition_ts(monkeypatch, tmp_path):
    roles = (
        "reactant",
        "intermediate",
        "addition-transition-state",
        "released-product",
    )
    route = closeout.AcceptedRoute(
        evidence_root=tmp_path,
        clusters={role: cluster(role) for role in roles},
        frequencies={
            "reactant": frequency(-10.0),
            "intermediate": frequency(-10.01),
            "addition-transition-state": frequency(-9.95),
            "released-product": frequency(-10.02),
        },
        validation={
            "route": [
                "reactant",
                "addition-transition-state",
                "intermediate",
                "barrierless-shelf",
                "released-product",
            ]
        },
    )
    calls = []

    def checkpoint(path, current, settings, method):
        calls.append((path.name, current.name, method))
        role_offset = 0.05 if current.name == "addition-transition-state" else 0.0
        method_offset = {
            closeout.SVP_METHOD: -20.0,
            closeout.production.PRODUCTION_METHOD: -30.0,
            closeout.production.B3LYP_D4_METHOD: -40.0,
        }[method]
        payload = {
            "method": method,
            "electronic_hartree": method_offset + role_offset,
            "geometry_fingerprint": closeout.production.frequency_geometry_fingerprint(
                current
            ),
            "settings_fingerprint": closeout.production.frequency_settings_fingerprint(
                settings
            ),
        }
        path.write_text(json.dumps(payload))
        return payload["electronic_hartree"]

    monkeypatch.setattr(closeout.production, "checkpoint_energy", checkpoint)
    monkeypatch.setattr(
        closeout.production,
        "thermo",
        lambda current_frequency, electronic: SimpleNamespace(
            gibbs=electronic * closeout.HARTREE_TO_KJ
        ),
    )

    result = closeout.run_dft(route, use_gpu=False)

    assert len(calls) == 12
    assert {name for _, name, _ in calls} == set(roles)
    assert "cleavage-transition-state" not in {name for _, name, _ in calls}
    assert result["cleavage_verified_saddle"] is False
    assert result["rejected_banked_transition_state_used"] is False
    assert (tmp_path / "production-closeout/dft-summary.json").is_file()
