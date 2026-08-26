"""CPU-only contracts for the A2a sequential production-path rebuild."""

from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from quarry.clusters import Cluster
from quarry.pipeline import DftSettings, FrequencyResult, frequency_geometry_fingerprint
from scripts import a2a_path_rebuild as a2a


def state(name: str, basin: str) -> Cluster:
    if basin == "reactant":
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
    elif basin == "intermediate":
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.6, 0.0, 0.0],
                [3.1, 0.0, 0.0],
                [0.9, 0.0, 0.0],
                [3.7, 0.7, 0.0],
                [0.0, 5.0, 0.0],
                [0.9, 5.0, 0.0],
            ]
        )
    elif basin == "product":
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
        raise ValueError(basin)
    return Cluster(name, ["O", "Si", "O", "H", "H", "O", "H"], coords)


def frequency(
    cluster: Cluster,
    settings: DftSettings,
    *,
    imaginary_cm: float | None,
    mode: np.ndarray | None = None,
    electronic_hartree: float = -10.0,
) -> FrequencyResult:
    imaginary = np.array([]) if imaginary_cm is None else np.array([imaginary_cm])
    modes = None if mode is None else np.asarray([mode], dtype=float)
    return FrequencyResult(
        frequencies_cm=np.array([100.0, 200.0, 300.0]),
        imaginary_cm=imaginary,
        electronic_hartree=electronic_hartree,
        molar_mass_kg=0.1,
        rotational_temperatures_k=(1.0, 2.0, 3.0),
        linear=False,
        imaginary_mode=None if modes is None else modes[0],
        imaginary_modes=modes,
        geometry_fingerprint=frequency_geometry_fingerprint(cluster),
        settings_fingerprint=a2a.frequency_settings_fingerprint(settings),
    )


def write_xyz(path: Path, cluster: Cluster) -> None:
    path.write_text(a2a.a2.exact_xyz(cluster))


def test_predeclared_hessian_step_contract_accepts_stable_reaction_mode():
    cluster = state("ts", "intermediate")
    settings = DftSettings(xc="hf", basis="sto-3g")
    mode = np.zeros_like(cluster.coords)
    mode[0, 0] = 1.0
    vector = mode * 2.0
    primary = frequency(cluster, settings, imaginary_cm=100.0, mode=mode)
    secondary = frequency(cluster, settings, imaginary_cm=108.0, mode=-mode)

    receipt = a2a.require_hessian_step_stability(primary, secondary, vector)

    assert receipt["mode_cosine"] == pytest.approx(1.0)
    assert receipt["frequency_drift_cm"] == pytest.approx(8.0)
    assert receipt["primary_reaction_overlap"] == pytest.approx(1.0)


def test_hessian_step_contract_rejects_frequency_and_mode_drift():
    cluster = state("ts", "intermediate")
    settings = DftSettings(xc="hf", basis="sto-3g")
    x_mode = np.zeros_like(cluster.coords)
    x_mode[0, 0] = 1.0
    y_mode = np.zeros_like(cluster.coords)
    y_mode[0, 1] = 1.0
    primary = frequency(cluster, settings, imaginary_cm=100.0, mode=x_mode)

    with pytest.raises(RuntimeError, match="frequency changed"):
        a2a.require_hessian_step_stability(
            primary,
            frequency(cluster, settings, imaginary_cm=140.0, mode=x_mode),
            x_mode,
        )
    with pytest.raises(RuntimeError, match="mode cosine"):
        a2a.require_hessian_step_stability(
            primary,
            frequency(cluster, settings, imaginary_cm=105.0, mode=y_mode),
            x_mode,
        )


def test_fd_frequency_checkpoint_binds_step_and_mode_vectors(monkeypatch, tmp_path):
    cluster = state("ts", "intermediate")
    settings = DftSettings(xc="hf", basis="sto-3g")
    mode = np.zeros_like(cluster.coords)
    mode[0, 0] = 1.0
    expected = frequency(cluster, settings, imaginary_cm=90.0, mode=mode)
    calls = []

    def compute(current, current_settings, *, step_bohr):
        calls.append(step_bohr)
        return expected

    monkeypatch.setattr(a2a, "frequencies_finite_difference", compute)
    path = tmp_path / "frequency.json"
    first = a2a.checkpoint_fd_frequency(path, cluster, settings, step_bohr=1e-3)
    second = a2a.checkpoint_fd_frequency(path, cluster, settings, step_bohr=1e-3)
    assert calls == [1e-3]
    assert first.imaginary_modes is not None
    assert second.imaginary_modes is not None
    assert np.array_equal(first.imaginary_modes, second.imaginary_modes)
    payload = json.loads(path.read_text())
    assert payload["hessian_step_bohr"] == 1e-3

    with pytest.raises(ValueError, match="identity drift"):
        a2a.checkpoint_fd_frequency(path, cluster, settings, step_bohr=2e-3)


def test_exact_typed_minimum_basin_contract():
    a2a.require_basin(state("R", "reactant"), 2, a2a.REACTANT_BASIN)
    a2a.require_basin(state("I", "intermediate"), 2, a2a.ASSOCIATIVE_BASIN)
    a2a.require_basin(state("P", "product"), 2, a2a.PRODUCT_BASIN)
    with pytest.raises(RuntimeError, match="expected basin"):
        a2a.require_basin(state("R", "reactant"), 2, a2a.PRODUCT_BASIN)


def test_full_route_records_two_independent_segments_and_rate_limiting_cc_roles(
    monkeypatch, tmp_path
):
    references = {
        "reactant": state("reactant", "reactant"),
        "intermediate": state("intermediate", "intermediate"),
        "product": state("product", "product"),
    }
    paths = {}
    for role, cluster in references.items():
        path = tmp_path / f"{role}.xyz"
        write_xyz(path, cluster)
        paths[role] = path

    monkeypatch.setattr(
        a2a.a2,
        "optimize_minimum",
        lambda current, settings, max_steps, trajectory: current,
    )

    def fake_fd(path, cluster, settings, *, step_bohr):
        return frequency(cluster, settings, imaginary_cm=None)

    monkeypatch.setattr(a2a, "checkpoint_fd_frequency", fake_fd)

    segment_index = {"addition": 0, "cleavage": 1}

    def fake_segment(run_dir, spec, start, end, settings, **kwargs):
        coords = start.coords.copy()
        coords[1, 2] += 0.2 + segment_index[spec.slug] * 0.1
        ts = replace(
            start,
            name=f"{spec.slug}-transition-state",
            coords=coords,
        )
        mode = end.coords - start.coords
        return {
            "spec": asdict_segment(spec),
            "transition_state": ts,
            "frequency": frequency(
                ts,
                settings,
                imaginary_cm=100.0 + 10.0 * segment_index[spec.slug],
                mode=mode,
                electronic_hartree=-9.95 + 0.01 * segment_index[spec.slug],
            ),
            "hessian_step_stability": {"mode_cosine": 1.0},
            "full_irc": {
                "actual": [
                    list(a2a.a2.si_neutral_signature(start, 2)),
                    list(a2a.a2.si_neutral_signature(end, 2)),
                ]
            },
            "artifacts": {},
        }

    monkeypatch.setattr(a2a, "run_segment", fake_segment)

    role_offsets = {
        "reactant": 0.0,
        "intermediate": 0.01,
        "product": -0.01,
        "addition-transition-state": 0.04,
        "cleavage-transition-state": 0.06,
    }

    def fake_energy(path, cluster, settings, method):
        base = -20.0 if settings.xc == "wb97m-v" else -15.0
        role = next(name for name in role_offsets if path.name.startswith(name + "."))
        return base + role_offsets[role]

    monkeypatch.setattr(a2a.a2, "checkpoint_energy", fake_energy)
    run_dir = tmp_path / "run"
    args = argparse.Namespace(
        reactant_reference=paths["reactant"],
        intermediate_reference=paths["intermediate"],
        product_reference=paths["product"],
        run_dir=run_dir,
        attacker_index=2,
        active_indices="0,1,2,3,4",
        gpu=False,
        minimum_steps=5,
        neb_images=5,
        neb_pre_steps=5,
        neb_steps=5,
        saddle_steps=5,
        irc_steps=5,
    )

    assert a2a.execute_with_status(args) == 0
    result = json.loads((run_dir / "results.json").read_text())
    assert result["mechanism"] == "sequential-associative"
    assert set(result["segments"]) == {"addition", "cleavage"}
    assert result["route"] == [
        "reactant",
        "addition-transition-state",
        "intermediate",
        "cleavage-transition-state",
        "product",
    ]
    assert result["rate_limiting_segment"] == "cleavage"
    assert result["cc_calibration_roles"]["ts"].endswith(
        "cleavage/transition-state.xyz"
    )
    assert result["rejected_banked_transition_state_used"] is False
    assert (run_dir / "store.sqlite").is_file()
    status = json.loads((run_dir / "run_status.json").read_text())
    assert status["status"] == "completed"


def test_local_neb_tangent_comes_from_neighbors_of_exact_peak(tmp_path):
    template = state("reactant", "reactant")
    checkpoint_root = tmp_path / "checkpoints"
    checkpoint = checkpoint_root / "climb-final"
    checkpoint.mkdir(parents=True)
    images = []
    for index in range(5):
        current = replace(template, coords=template.coords.copy())
        current.coords[0, 0] += index * 0.1
        current.coords[1, 1] += index * index * 0.05
        path = checkpoint / f"image-{index:02d}.xyz"
        write_xyz(path, current)
        images.append(current)
    (checkpoint / "manifest.json").write_text(
        json.dumps(
            {
                "stage": "climb-final",
                "converged": True,
                "images": [f"image-{index:02d}.xyz" for index in range(5)],
            }
        )
    )
    (checkpoint_root / "latest.json").write_text(
        json.dumps({"checkpoint": checkpoint.name})
    )

    tangent, receipt = a2a.final_neb_tangent(
        checkpoint_root,
        images[2],
        template,
        [0, 1],
    )

    expected = images[3].coords - images[1].coords
    masked = np.zeros_like(expected)
    masked[[0, 1]] = expected[[0, 1]]
    expected = masked / np.linalg.norm(masked)
    assert np.allclose(tangent, expected)
    assert receipt["peak_index"] == 2
    assert receipt["strategy"] == a2a.SELLA_MODE_STRATEGY
    left_radius = np.linalg.norm(images[1].coords[[0, 1]] - images[2].coords[[0, 1]])
    right_radius = np.linalg.norm(images[3].coords[[0, 1]] - images[2].coords[[0, 1]])
    assert receipt["left_active_radius_a"] == pytest.approx(left_radius)
    assert receipt["right_active_radius_a"] == pytest.approx(right_radius)
    assert receipt["local_trust_radius_a"] == pytest.approx(
        min(left_radius, right_radius)
    )
    assert receipt["local_guard_radius_a"] == pytest.approx(
        max(left_radius, right_radius, 1.5 * min(left_radius, right_radius))
    )
    assert receipt["sella_delta_max_a"] <= a2a.LOCAL_SELLA_DELTA_MAX_A


def asdict_segment(spec: a2a.SegmentSpec) -> dict[str, str]:
    return {
        "slug": spec.slug,
        "start_role": spec.start_role,
        "end_role": spec.end_role,
    }


def test_active_indices_fail_closed():
    assert a2a.active_indices("0,2,4", 5) == [0, 2, 4]
    with pytest.raises(ValueError, match="non-empty and unique"):
        a2a.active_indices("0,0", 5)
    with pytest.raises(ValueError, match="out-of-range"):
        a2a.active_indices("0,5", 5)


def test_fixed_distance_checkpoint_identity_survives_json_roundtrip():
    targets = [(1, 15, 1.75), (0, 17, 0.99)]
    serialized = a2a.serialized_distance_targets(targets)
    assert json.loads(json.dumps(serialized)) == serialized
    assert serialized == [[1, 15, 1.75], [0, 17, 0.99]]


def test_sella_strategy_records_active_subspace_and_reconditioning():
    assert "active-subspace" in a2a.SELLA_MODE_STRATEGY
    assert "internal-conditioned" in a2a.SELLA_MODE_STRATEGY
    assert "cartesian-trust" in a2a.LOCAL_TRUST_SELLA_MODE_STRATEGY
