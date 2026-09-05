"""CPU-only contracts for the A2a barrierless production closeout."""

from __future__ import annotations

import json
import sqlite3
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
    with pytest.raises(ValueError, match="finite energy matrix"):
        closeout._classification_and_release(tmp_path)

    write_classification(tmp_path, outcome="interior-saddle")
    with pytest.raises(ValueError, match="finite energy matrix"):
        closeout._classification_and_release(tmp_path)


def test_frequency_receipt_is_bound_to_geometry_settings_and_hessian(tmp_path):
    current = cluster("reactant")
    settings = closeout.production.settings(use_gpu=False)[0]
    payload = {
        "electronic_hartree": -10.0,
        "frequencies_cm": [100.0],
        "imaginary_cm": [],
        "molar_mass_kg": 0.1,
        "rotational_temperatures_k": [1.0, 2.0, 3.0],
        "linear": False,
        "geometry_fingerprint": closeout.production.frequency_geometry_fingerprint(
            current
        ),
        "settings_fingerprint": closeout.production.frequency_settings_fingerprint(
            settings
        ),
        "hessian_method": "finite-difference-gradient",
    }
    path = tmp_path / "frequency.json"
    path.write_text(json.dumps(payload))
    assert closeout._load_frequency(path, current).n_imaginary == 0
    for field, bad_value, match in (
        ("geometry_fingerprint", "wrong", "geometry drift"),
        ("settings_fingerprint", "wrong", "settings drift"),
        ("hessian_method", "analytic", "Hessian method drift"),
    ):
        corrupted = dict(payload)
        corrupted[field] = bad_value
        path.write_text(json.dumps(corrupted))
        with pytest.raises(ValueError, match=match):
            closeout._load_frequency(path, current)


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


@pytest.mark.parametrize(
    "method",
    [
        closeout.production.PRODUCTION_METHOD,
        closeout.production.B3LYP_D4_METHOD,
    ],
)
def test_production_scf_uses_receipted_bounded_newton(monkeypatch, tmp_path, method):
    current = cluster("reactant")
    settings = closeout.method_settings(use_gpu=False)[method]

    class FakeScf:
        converged = False
        max_cycle = 150

        def newton(self):
            return self

        def kernel(self, **kwargs):
            assert kwargs == {}
            self.converged = True
            return -30.5

    monkeypatch.setattr(closeout.pipeline, "build_mol", lambda *args: object())
    monkeypatch.setattr(closeout.pipeline, "_make_scf", lambda *args: FakeScf())
    path = tmp_path / "energy.json"

    value = closeout.checkpoint_closeout_energy(
        path,
        current,
        settings,
        method,
    )

    assert value == -30.5
    payload = json.loads(path.read_text())
    assert payload["converged"] is True
    assert payload["convergence_route"] == "newton-first"
    assert payload["scf_contract"] == "bounded-direct-diis-or-newton-v1"


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

    monkeypatch.setattr(closeout, "checkpoint_closeout_energy", checkpoint)
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
    store_path = tmp_path / "production-closeout/store.sqlite"
    assert store_path.is_file()
    with sqlite3.connect(store_path) as connection:
        assert connection.execute("SELECT count(*) FROM structures").fetchone()[0] == 4
        assert connection.execute("SELECT count(*) FROM jobs").fetchone()[0] == 17
        assert connection.execute(
            "SELECT value FROM results WHERE key = 'production_dg_dagger'"
        ).fetchone()[0] == pytest.approx(result["production_dg_dagger_kj_mol"])
        analysis = connection.execute(
            "SELECT status, detail FROM jobs WHERE kind = 'analysis'"
        ).fetchone()
        assert analysis[0] == "done"
        assert json.loads(analysis[1])["rejected_banked_transition_state_used"] is False


def test_dft_resume_with_complete_receipts_skips_gpu_import_preflight(
    monkeypatch, tmp_path
):
    roles = (
        "reactant",
        "intermediate",
        "addition-transition-state",
        "released-product",
    )
    route = closeout.AcceptedRoute(
        evidence_root=tmp_path,
        clusters={role: cluster(role) for role in roles},
        frequencies={role: frequency(-10.0) for role in roles},
        validation={"route": list(roles)},
    )
    energy_root = tmp_path / "production-closeout/energies"
    energy_root.mkdir(parents=True)
    for method, settings in closeout.method_settings(use_gpu=True).items():
        slug = method.split("/")[0].replace("(", "-").replace(")", "")
        for role, current in route.clusters.items():
            (energy_root / f"{role}.{slug}.energy.json").write_text(
                json.dumps(
                    {
                        "method": method,
                        "electronic_hartree": -10.0,
                        "geometry_fingerprint": (
                            closeout.production.frequency_geometry_fingerprint(current)
                        ),
                        "settings_fingerprint": (
                            closeout.production.frequency_settings_fingerprint(settings)
                        ),
                    }
                )
            )
    monkeypatch.setattr(
        closeout.path_driver,
        "preflight_gpu_contraction_engine",
        lambda *_: pytest.fail("complete receipt resume must not import GPU engine"),
    )
    monkeypatch.setattr(
        closeout.production,
        "thermo",
        lambda current_frequency, electronic: SimpleNamespace(
            gibbs=electronic * closeout.HARTREE_TO_KJ
        ),
    )

    result = closeout.run_dft(route, use_gpu=True)

    assert result["status"] == "dft-completed"
    assert (tmp_path / "production-closeout/store.sqlite").is_file()


def test_release_seed_energy_must_match_classified_product_cell(tmp_path):
    cell_dir = tmp_path / "cleavage/coupled-scan-v1/r07-c07"
    cell_dir.mkdir(parents=True)
    (cell_dir / "electronic-energy.json").write_text(
        json.dumps({"electronic_hartree": -10.0})
    )

    with pytest.raises(ValueError, match="classified product cell"):
        closeout._authoritative_seed_electronic_hartree(
            tmp_path,
            product_cell=[7, 7],
            claimed_seed_hartree=-9.5,
        )

    assert closeout._authoritative_seed_electronic_hartree(
        tmp_path,
        product_cell=[7, 7],
        claimed_seed_hartree=-10.0,
    ) == pytest.approx(-10.0)
