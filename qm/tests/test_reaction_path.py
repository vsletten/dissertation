from __future__ import annotations

import math

import numpy as np
import pytest

from quarry.reaction_path import (
    build_mass_scaled_path,
    curvature_effective_mass_profile,
    path_tangents_and_curvature,
    project_transverse_hessian,
    vibrationally_adiabatic_potential,
)


def _rotation_z(angle: float) -> np.ndarray:
    return np.array(
        [
            [math.cos(angle), -math.sin(angle), 0.0],
            [math.sin(angle), math.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )


def _mass_weighted_rigid_motions(
    coordinates: np.ndarray, masses: np.ndarray
) -> np.ndarray:
    centred = coordinates - np.average(coordinates, axis=0, weights=masses)
    roots = np.sqrt(masses)
    translations = [(roots[:, None] * axis).reshape(-1) for axis in np.eye(3)]
    rotations = [
        (roots[:, None] * np.cross(axis, centred)).reshape(-1) for axis in np.eye(3)
    ]
    return np.column_stack([*translations, *rotations])


def _distinct_transverse_hessian(
    coordinates: np.ndarray, masses: np.ndarray, tangent: np.ndarray
) -> np.ndarray:
    rigid = _mass_weighted_rigid_motions(coordinates, masses)
    rigid_basis = np.linalg.qr(rigid)[0]
    tangent_vector = tangent.reshape(-1)
    tangent_vector -= rigid_basis @ (rigid_basis.T @ tangent_vector)
    tangent_vector /= np.linalg.norm(tangent_vector)
    removed = np.column_stack([rigid_basis, tangent_vector])
    transverse = np.linalg.svd(removed.T, full_matrices=True)[2][7:].T
    mass_weighted_hessian = transverse @ np.diag([2.0, 5.0]) @ transverse.T
    factors = np.repeat(np.sqrt(masses), 3)
    return factors[:, None] * mass_weighted_hessian * factors[None, :]


def test_straight_path_has_mass_scaled_coordinate_and_zero_curvature():
    bond_lengths = np.array([1.0, 1.2, 1.5, 1.9, 2.4])
    coordinates = np.zeros((len(bond_lengths), 2, 3))
    coordinates[:, 1, 0] = bond_lengths
    masses = np.array([1.0, 3.0])

    path = build_mass_scaled_path(
        coordinates, masses, transition_state_index=2, reference_mass_amu=1.0
    )
    expected = math.sqrt(0.75) * (bond_lengths - bond_lengths[2])
    assert path.coordinate_angstrom == pytest.approx(expected, abs=1.0e-12)

    geometry = path_tangents_and_curvature(
        path.coordinate_angstrom, path.mass_scaled_coordinates
    )
    assert geometry.curvature_magnitudes_per_angstrom == pytest.approx(
        np.zeros(len(bond_lengths)), abs=2.0e-12
    )
    assert np.linalg.norm(geometry.unit_tangents, axis=1) == pytest.approx(1.0)


def test_nonuniform_circle_recovers_curvature():
    radius = 2.5
    theta = np.array([-0.8, -0.55, -0.2, 0.0, 0.17, 0.48, 0.9])
    coordinate = radius * theta
    points = np.column_stack(
        [radius * np.cos(theta), radius * np.sin(theta), np.zeros_like(theta)]
    )

    result = path_tangents_and_curvature(coordinate, points)

    assert result.curvature_magnitudes_per_angstrom[1:-1] == pytest.approx(
        np.full(len(theta) - 2, 1.0 / radius), rel=0.09
    )
    assert np.einsum(
        "ij,ij->i", result.unit_tangents, result.curvature_vectors_per_angstrom
    ) == pytest.approx(np.zeros(len(theta)), abs=1.0e-12)


def test_path_is_invariant_to_rigid_motion_and_reversal():
    base = np.array([[0.0, 0.0, 0.0], [1.1, 0.2, 0.0], [-0.2, 0.9, 0.4]])
    displacement = np.array(
        [[-0.08, 0.03, 0.02], [0.05, -0.04, 0.01], [0.03, 0.01, -0.03]]
    )
    parameter = np.array([-1.0, -0.4, 0.0, 0.35, 1.0])
    coordinates = np.asarray([base + value * displacement for value in parameter])
    masses = np.array([12.0, 16.0, 1.0])
    reference = build_mass_scaled_path(
        coordinates, masses, transition_state_index=2, reference_mass_amu=2.0
    )

    moved = np.asarray(
        [
            point @ _rotation_z(0.31 * (index + 1)).T
            + np.array([2.0 * index, -index, 0.4 * index])
            for index, point in enumerate(coordinates)
        ]
    )
    transformed = build_mass_scaled_path(
        moved, masses, transition_state_index=2, reference_mass_amu=2.0
    )
    reversed_path = build_mass_scaled_path(
        coordinates[::-1], masses, transition_state_index=2, reference_mass_amu=2.0
    )

    assert transformed.coordinate_angstrom == pytest.approx(
        reference.coordinate_angstrom, abs=2.0e-12
    )
    assert reversed_path.coordinate_angstrom == pytest.approx(
        -reference.coordinate_angstrom[::-1], abs=2.0e-12
    )


def test_nonlinear_transverse_projection_recovers_distinct_modes_and_constraints():
    coordinates = np.array([[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [0.1, 0.9, 0.4]])
    masses = np.array([12.0, 16.0, 1.0])
    tangent = np.array(
        [[0.12, -0.03, 0.05], [-0.08, 0.07, -0.01], [-0.04, -0.04, -0.04]]
    )
    positive_hessian = _distinct_transverse_hessian(coordinates, masses, tangent.copy())

    result = project_transverse_hessian(
        coordinates,
        masses,
        positive_hessian,
        tangent,
        negative_eigenvalue_tolerance=1.0e-10,
    )

    assert result.eigenvalues.shape == (2,)
    assert result.mass_weighted_eigenvectors.shape == (2, 3, 3)
    assert result.eigenvalues == pytest.approx([2.0, 5.0], abs=2.0e-12)

    eigenvectors = result.mass_weighted_eigenvectors.reshape(2, -1)
    rigid = _mass_weighted_rigid_motions(coordinates, masses)
    assert eigenvectors @ rigid == pytest.approx(np.zeros((2, 6)), abs=2.0e-12)
    assert eigenvectors @ tangent.reshape(-1) == pytest.approx(np.zeros(2), abs=2.0e-12)
    assert eigenvectors @ eigenvectors.T == pytest.approx(np.eye(2), abs=2.0e-12)
    factors = np.repeat(np.sqrt(masses), 3)
    mass_weighted_hessian = positive_hessian / (factors[:, None] * factors[None, :])
    for eigenvalue, eigenvector in zip(result.eigenvalues, eigenvectors, strict=True):
        assert mass_weighted_hessian @ eigenvector == pytest.approx(
            eigenvalue * eigenvector, abs=2.0e-12
        )


def test_nonlinear_transverse_projection_rejects_significant_negative_mode():
    coordinates = np.array([[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [0.1, 0.9, 0.4]])
    masses = np.array([12.0, 16.0, 1.0])
    tangent = np.array(
        [[0.12, -0.03, 0.05], [-0.08, 0.07, -0.01], [-0.04, -0.04, -0.04]]
    )
    positive = project_transverse_hessian(
        coordinates,
        masses,
        _distinct_transverse_hessian(coordinates, masses, tangent.copy()),
        tangent,
        negative_eigenvalue_tolerance=1.0e-10,
    )

    mass_factors = np.repeat(np.sqrt(masses), 3)
    mode = positive.mass_weighted_eigenvectors[0].reshape(-1)
    negative_mass_weighted = np.eye(9) - 2.0 * np.outer(mode, mode)
    negative_cartesian = (
        mass_factors[:, None] * negative_mass_weighted * mass_factors[None, :]
    )
    with pytest.raises(ValueError, match="negative transverse eigenvalue"):
        project_transverse_hessian(
            coordinates,
            masses,
            negative_cartesian,
            tangent,
            negative_eigenvalue_tolerance=1.0e-10,
        )


def test_transverse_projection_is_rigid_rotation_invariant():
    coordinates = np.array([[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [0.1, 0.9, 0.4]])
    masses = np.array([12.0, 16.0, 1.0])
    tangent = np.array(
        [[0.12, -0.03, 0.05], [-0.08, 0.07, -0.01], [-0.04, -0.04, -0.04]]
    )
    rng = np.random.default_rng(20260907)
    raw = rng.normal(size=(9, 9))
    mass_weighted_hessian = raw.T @ raw + np.eye(9)
    factors = np.repeat(np.sqrt(masses), 3)
    hessian = factors[:, None] * mass_weighted_hessian * factors[None, :]
    original = project_transverse_hessian(coordinates, masses, hessian, tangent)

    rotation = _rotation_z(0.73)
    cartesian_rotation = np.kron(np.eye(3), rotation)
    rotated = project_transverse_hessian(
        coordinates @ rotation.T + np.array([4.0, -3.0, 1.5]),
        masses,
        cartesian_rotation @ hessian @ cartesian_rotation.T,
        tangent @ rotation.T,
    )
    assert rotated.eigenvalues == pytest.approx(original.eigenvalues, rel=2.0e-12)


def test_vibrationally_adiabatic_potential_uses_exact_wavenumber_conversion():
    electronic = np.array([1.0, 2.0])
    frequencies = np.array([[100.0, 300.0], [200.0, 500.0]])
    correction = np.array([0.25, -0.5])

    result = vibrationally_adiabatic_potential(
        electronic, frequencies, high_level_correction_kj_mol=correction
    )

    quantum_kj_mol_per_cm = (
        6.62607015e-34 * 299792458.0 * 100.0 * 6.02214076e23 / 1000.0
    )
    expected_zpe = 0.5 * quantum_kj_mol_per_cm * frequencies.sum(axis=1)
    assert result.zero_point_energy_kj_mol == pytest.approx(expected_zpe, rel=1.0e-15)
    assert result.potential_kj_mol == pytest.approx(
        electronic + expected_zpe + correction, rel=1.0e-15
    )


def test_mode_resolved_curvature_profile_is_sign_and_permutation_invariant():
    coordinate = np.array([-0.6, -0.2, 0.0, 0.35, 0.8])
    curvature_vectors = np.array(
        [[0.3, 0.4, 0.0, 0.0], [0.2, 0.5, 0.0, 0.0]] * 2 + [[0.35, 0.45, 0.0, 0.0]]
    )
    modes = np.broadcast_to(
        np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]]),
        (len(coordinate), 2, 4),
    ).copy()
    frequencies = np.column_stack(
        [np.linspace(800.0, 1000.0, len(coordinate)), np.full(len(coordinate), 1200.0)]
    )

    original = curvature_effective_mass_profile(
        coordinate,
        curvature_vectors,
        modes,
        frequencies,
        reference_mass_amu=1.0,
    )
    changed = curvature_effective_mass_profile(
        coordinate,
        curvature_vectors,
        -modes[:, ::-1],
        frequencies[:, ::-1],
        reference_mass_amu=1.0,
    )

    assert changed.turning_length_angstrom == pytest.approx(
        original.turning_length_angstrom, rel=2.0e-14
    )
    assert changed.effective_mass_amu == pytest.approx(
        original.effective_mass_amu, rel=2.0e-14
    )
    assert np.all(original.effective_mass_amu > 0.0)
    assert np.all(original.effective_mass_amu <= 1.0)

    degenerate_frequencies = np.full_like(frequencies, 900.0)
    degenerate = curvature_effective_mass_profile(
        coordinate, curvature_vectors, modes, degenerate_frequencies
    )
    angle = 0.37
    mode_rotation = np.array(
        [[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]]
    )
    rotated_modes = np.einsum("mn,pni->pmi", mode_rotation, modes)
    rotated = curvature_effective_mass_profile(
        coordinate, curvature_vectors, rotated_modes, degenerate_frequencies
    )
    assert rotated.turning_length_angstrom == pytest.approx(
        degenerate.turning_length_angstrom, rel=2.0e-14
    )
    assert rotated.effective_mass_amu == pytest.approx(
        degenerate.effective_mass_amu, rel=2.0e-14
    )


def test_pilgrim_liu_two_mode_turning_length_golden_equation():
    coordinate = np.array([-1.0, 0.0, 1.0])
    components = np.array([0.3, 0.4])
    curvature_vectors = np.tile([*components, 0.0], (3, 1))
    modes = np.broadcast_to(np.eye(3)[:2], (3, 2, 3)).copy()
    frequencies = np.tile([400.0, 100.0], (3, 1))

    result = curvature_effective_mass_profile(
        coordinate, curvature_vectors, modes, frequencies
    )

    hbar_j_s = 6.62607015e-34 / (2.0 * math.pi)
    omega = 2.0 * math.pi * 2.99792458e10 * frequencies[0]
    amplitudes = np.sqrt(hbar_j_s / (1.66053906660e-27 * omega)) / 1.0e-10
    kappa = np.linalg.norm(components)
    expected = math.sqrt(kappa) * np.sum((components / amplitudes**2) ** 2) ** (-0.25)
    assert expected == pytest.approx(0.36506336503450354, rel=2.0e-15)
    assert result.turning_length_angstrom == pytest.approx(
        np.full(3, expected), rel=2.0e-15
    )


def test_exact_straight_point_uses_continuous_neighboring_turning_length():
    coordinate = np.array([-1.0, 0.0, 1.0])
    curvature_vectors = np.array([[0.3, 0.4, 0.0], [0.0, 0.0, 0.0], [0.3, 0.4, 0.0]])
    modes = np.broadcast_to(np.eye(3)[:2], (3, 2, 3)).copy()
    frequencies = np.tile([400.0, 100.0], (3, 1))

    result = curvature_effective_mass_profile(
        coordinate, curvature_vectors, modes, frequencies
    )

    assert result.turning_length_angstrom[1] == pytest.approx(
        0.5 * (result.turning_length_angstrom[0] + result.turning_length_angstrom[2])
    )
    assert result.turning_length_derivative[1] == pytest.approx(0.0, abs=1.0e-14)
    assert result.mode_curvature_components_per_angstrom[1] == pytest.approx([0.0, 0.0])
    assert result.effective_mass_amu[1] == pytest.approx(1.0)
    assert np.all(np.isfinite(result.turning_length_angstrom))


def test_duplicate_nonfinite_and_nonpositive_inputs_are_rejected():
    coordinates = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            [[1.0, 2.0, 3.0], [2.0, 2.0, 3.0]],
            [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]],
        ]
    )
    with pytest.raises(ValueError, match="duplicate"):
        build_mass_scaled_path(coordinates, [1.0, 1.0], transition_state_index=1)

    bad = coordinates.copy()
    bad[1, 0, 0] = np.nan
    with pytest.raises(ValueError, match="finite"):
        build_mass_scaled_path(bad, [1.0, 1.0], transition_state_index=1)

    with pytest.raises(ValueError, match="positive"):
        vibrationally_adiabatic_potential([0.0], [[100.0, -2.0]])
