from __future__ import annotations

import math

import numpy as np
import pytest

from quarry.tunneling import (
    corrected_rate_from_reference,
    harmonic_turning_length_angstrom,
    sct_kappa,
    small_curvature_effective_mass,
)


def test_small_curvature_mass_is_bounded_and_reduces_at_bend():
    mass = small_curvature_effective_mass(
        1.0,
        [0.0, 0.5, 0.0],
        harmonic_turning_length_angstrom(1000.0),
    )
    assert mass[[0, 2]].tolist() == pytest.approx([1.0, 1.0])
    assert 0.0 < mass[1] < 1.0


def test_curvature_increases_explicit_path_tunnelling():
    coordinate = [0.0, 0.5, 1.0, 1.5, 2.0]
    potential = [0.0, 20.0, 30.0, 18.0, -5.0]
    zct = sct_kappa(coordinate, potential, 1.0, 20.0, quadrature_order=64)
    sct_mass = small_curvature_effective_mass(1.0, [0.0, 0.4, 0.6, 0.3, 0.0], 0.15)
    sct = sct_kappa(coordinate, potential, sct_mass, 20.0, quadrature_order=64)
    assert sct.log_kappa > zct.log_kappa
    assert sct.ground_probability > zct.ground_probability
    assert sct.barrier_kj_mol == 30.0
    assert sct.energy_floor_kj_mol == 0.0


def test_path_reversal_is_invariant():
    coordinate = np.array([0.0, 0.4, 1.1, 1.8])
    potential = np.array([2.0, 16.0, 25.0, -3.0])
    mass = np.array([1.0, 0.9, 0.95, 1.0])
    forward = sct_kappa(coordinate, potential, mass, 50.0, quadrature_order=64)
    reverse = sct_kappa(
        coordinate[::-1], potential[::-1], mass[::-1], 50.0, quadrature_order=64
    )
    assert reverse.log_kappa == pytest.approx(forward.log_kappa, rel=2.0e-12)
    assert reverse.ground_probability == pytest.approx(
        forward.ground_probability, rel=2.0e-12
    )


def test_action_and_quadrature_converge():
    coordinate = [0.0, 0.35, 0.9, 1.4, 2.0]
    potential = [0.0, 12.0, 28.0, 14.0, -2.0]
    coarse = sct_kappa(
        coordinate, potential, 1.0, 75.0, quadrature_order=64, path_grid_size=2049
    )
    fine = sct_kappa(
        coordinate, potential, 1.0, 75.0, quadrature_order=128, path_grid_size=4097
    )
    assert fine.log_kappa == pytest.approx(coarse.log_kappa, rel=2.0e-4)


def test_corrected_rate_replacement_stays_in_log_space():
    correction = sct_kappa(
        [0.0, 0.5, 1.0, 1.5],
        [0.0, 40.0, 30.0, -10.0],
        1.0,
        12.0,
        quadrature_order=64,
    )
    result = corrected_rate_from_reference(2.0e3, math.log(1.0e150), correction)
    assert math.isfinite(result)
    assert result > 0.0


def test_ground_action_matches_square_barrier():
    width_a = 1.25
    barrier_kj_mol = 18.0
    mass_amu = 1.0
    edge = 1.0e-6
    result = sct_kappa(
        [-edge, 0.0, width_a, width_a + edge],
        [0.0, barrier_kj_mol, barrier_kj_mol, 0.0],
        mass_amu,
        50.0,
        quadrature_order=64,
        path_grid_size=8193,
    )
    expected = (
        width_a
        * 1.0e-10
        * math.sqrt(
            2.0 * mass_amu * 1.66053906660e-27 * barrier_kj_mol * 1000.0 / 6.02214076e23
        )
        / 1.054571817e-34
    )
    assert result.ground_action == pytest.approx(expected, rel=3.0e-4)


@pytest.mark.parametrize(
    "coordinate,potential,mass",
    [
        ([0.0, 1.0], [0.0, 1.0], 1.0),
        ([0.0, 1.0, 0.5], [0.0, 2.0, 0.0], 1.0),
        ([0.0, 1.0, 2.0], [0.0, 2.0, 0.0], -1.0),
    ],
)
def test_invalid_paths_fail_closed(coordinate, potential, mass):
    with pytest.raises(ValueError):
        sct_kappa(coordinate, potential, mass, 20.0)
