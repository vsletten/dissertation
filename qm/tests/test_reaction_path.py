from __future__ import annotations

import math

import numpy as np
import pytest

from quarry.reaction_path import (
    build_mass_scaled_path,
    curvature_effective_mass_profile,
    molecular_path_tangents_and_curvature,
    path_tangents_and_curvature,
    project_transverse_hessian,
    project_vibrational_hessian,
    rotate_cartesian_hessian,
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


def _axis_angle_rotation(axis: np.ndarray | list[float], angle: float) -> np.ndarray:
    axis = np.asarray(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    cross_matrix = np.array(
        [
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ]
    )
    return (
        math.cos(angle) * np.eye(3)
        + (1.0 - math.cos(angle)) * np.outer(axis, axis)
        + math.sin(angle) * cross_matrix
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
    assert result.tangent_rigid_residual_magnitudes is None
    assert result.curvature_rigid_residual_magnitudes_per_angstrom is None


def test_propagated_alignment_makes_every_outward_segment_rigid_horizontal():
    parameter = np.linspace(-1.0, 1.0, 9)
    base = np.array(
        [[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [0.1, 0.9, 0.4], [0.3, -0.4, 0.8]]
    )
    coordinates = np.asarray(
        [
            base
            + np.array(
                [
                    [0.15 * value**2, 0.08 * value**3, -0.04 * value],
                    [0.12 * value, 0.18 * value**2, 0.07 * value**3],
                    [-0.09 * value, 0.05 * value**3, -0.12 * value**2],
                    [0.04 * value**2, -0.11 * value, 0.09 * value],
                ]
            )
            for value in parameter
        ]
    )
    masses = np.array([12.0, 16.0, 1.0, 14.0])
    transition_state_index = 4
    path = build_mass_scaled_path(
        coordinates, masses, transition_state_index=transition_state_index
    )

    outward_indices = [*range(transition_state_index - 1, -1, -1), *range(5, 9)]
    for point_index in outward_indices:
        predecessor_index = (
            point_index + 1 if point_index < transition_state_index else point_index - 1
        )
        displacement = (
            np.sqrt(masses)[:, None]
            * (
                path.aligned_coordinates_angstrom[point_index]
                - path.aligned_coordinates_angstrom[predecessor_index]
            )
        ).reshape(-1)
        local_rotations = _mass_weighted_rigid_motions(
            path.aligned_coordinates_angstrom[point_index], masses
        )[:, 3:]
        assert local_rotations.T @ displacement == pytest.approx(
            np.zeros(3), abs=2.0e-11
        )


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

    imposed_rotations = [
        _axis_angle_rotation(axis, angle)
        for axis, angle in zip(
            (
                [1.0, 2.0, 3.0],
                [-2.0, 1.0, 0.5],
                [0.3, -1.0, 2.0],
                [1.5, 0.2, -0.7],
                [-0.4, 2.0, 1.0],
            ),
            (0.31, -0.73, 1.11, -0.47, 0.88),
            strict=True,
        )
    ]
    moved = np.asarray(
        [
            point @ rotation.T + np.array([2.0 * index, -index, 0.4 * index])
            for index, (point, rotation) in enumerate(
                zip(coordinates, imposed_rotations, strict=True)
            )
        ]
    )
    transformed = build_mass_scaled_path(
        moved, masses, transition_state_index=2, reference_mass_amu=2.0
    )
    reversed_path = build_mass_scaled_path(
        coordinates[::-1], masses, transition_state_index=2, reference_mass_amu=2.0
    )

    reference_geometry = molecular_path_tangents_and_curvature(reference)
    transformed_geometry = molecular_path_tangents_and_curvature(transformed)
    reversed_geometry = molecular_path_tangents_and_curvature(reversed_path)
    transition_state_rotation = imposed_rotations[2]
    transformed_tangents_in_reference_frame = (
        transformed_geometry.unit_tangents.reshape(-1, 3, 3) @ transition_state_rotation
    ).reshape(len(parameter), -1)
    transformed_curvature_in_reference_frame = (
        transformed_geometry.curvature_vectors_per_angstrom.reshape(-1, 3, 3)
        @ transition_state_rotation
    ).reshape(len(parameter), -1)

    assert transformed.coordinate_angstrom == pytest.approx(
        reference.coordinate_angstrom, abs=2.0e-12
    )
    assert reversed_path.coordinate_angstrom == pytest.approx(
        -reference.coordinate_angstrom[::-1], abs=2.0e-12
    )
    assert transformed_tangents_in_reference_frame == pytest.approx(
        reference_geometry.unit_tangents, abs=2.0e-10
    )
    assert transformed_curvature_in_reference_frame == pytest.approx(
        reference_geometry.curvature_vectors_per_angstrom, abs=2.0e-8
    )
    assert reversed_path.aligned_coordinates_angstrom == pytest.approx(
        reference.aligned_coordinates_angstrom[::-1], abs=2.0e-12
    )
    assert reversed_geometry.unit_tangents == pytest.approx(
        -reference_geometry.unit_tangents[::-1], abs=2.0e-10
    )
    assert reversed_geometry.curvature_vectors_per_angstrom == pytest.approx(
        reference_geometry.curvature_vectors_per_angstrom[::-1], abs=2.0e-8
    )

    assert reference_geometry.tangent_rigid_residual_magnitudes is not None
    assert (
        reference_geometry.curvature_rigid_residual_magnitudes_per_angstrom is not None
    )
    assert np.all(reference_geometry.tangent_rigid_residual_magnitudes >= 0.0)
    assert np.all(
        reference_geometry.curvature_rigid_residual_magnitudes_per_angstrom >= 0.0
    )
    for point, tangent, curvature in zip(
        reference.aligned_coordinates_angstrom,
        reference_geometry.unit_tangents,
        reference_geometry.curvature_vectors_per_angstrom,
        strict=True,
    ):
        rigid = _mass_weighted_rigid_motions(point, masses)
        assert rigid.T @ tangent == pytest.approx(np.zeros(6), abs=2.0e-12)
        assert rigid.T @ curvature == pytest.approx(np.zeros(6), abs=2.0e-11)
        assert tangent @ curvature == pytest.approx(0.0, abs=2.0e-12)


def test_alignment_rotations_map_input_rows_and_coherently_rotate_hessians():
    base = np.array([[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [0.1, 0.9, 0.4]])
    displacement = np.array(
        [[0.03, -0.02, 0.01], [-0.04, 0.05, 0.02], [0.01, -0.03, 0.04]]
    )
    internal = np.asarray([base + value * displacement for value in (-0.4, 0.0, 0.5)])
    imposed_rotations = [_rotation_z(angle) for angle in (0.8, -0.35, 0.45)]
    coordinates = np.asarray(
        [
            point @ rotation.T + np.array([3.0 * index, -2.0, 0.7])
            for index, (point, rotation) in enumerate(
                zip(internal, imposed_rotations, strict=True)
            )
        ]
    )
    masses = np.array([12.0, 16.0, 1.0])
    path = build_mass_scaled_path(coordinates, masses, transition_state_index=1)

    for point, aligned, rotation in zip(
        coordinates,
        path.aligned_coordinates_angstrom,
        path.input_to_aligned_rotations,
        strict=True,
    ):
        centred = point - np.average(point, axis=0, weights=masses)
        assert centred @ rotation == pytest.approx(aligned, abs=2.0e-12)
        assert rotation.T @ rotation == pytest.approx(np.eye(3), abs=2.0e-12)
        assert np.linalg.det(rotation) == pytest.approx(1.0, abs=2.0e-12)

    point_index = 0
    rotation = path.input_to_aligned_rotations[point_index]
    aligned_tangent = (
        molecular_path_tangents_and_curvature(path)
        .unit_tangents[point_index]
        .reshape(3, 3)
    )
    input_tangent = aligned_tangent @ rotation.T
    rng = np.random.default_rng(20260907)
    raw = rng.normal(size=(9, 9))
    mass_weighted_hessian = raw.T @ raw + np.diag(np.arange(1.0, 10.0))
    factors = np.repeat(np.sqrt(masses), 3)
    input_hessian = factors[:, None] * mass_weighted_hessian * factors[None, :]

    input_modes = project_transverse_hessian(
        coordinates[point_index], masses, input_hessian, input_tangent
    )
    aligned_modes = project_transverse_hessian(
        path.aligned_coordinates_angstrom[point_index],
        masses,
        rotate_cartesian_hessian(input_hessian, rotation),
        aligned_tangent,
    )
    wrong_frame_modes = project_transverse_hessian(
        path.aligned_coordinates_angstrom[point_index],
        masses,
        input_hessian,
        aligned_tangent,
    )

    assert aligned_modes.eigenvalues == pytest.approx(
        input_modes.eigenvalues, rel=2.0e-12
    )
    assert not np.allclose(wrong_frame_modes.eigenvalues, input_modes.eigenvalues)

    tensor_hessian = input_hessian.reshape(3, 3, 3, 3)
    assert rotate_cartesian_hessian(tensor_hessian, rotation).reshape(
        9, 9
    ) == pytest.approx(rotate_cartesian_hessian(input_hessian, rotation), abs=2.0e-12)
    with pytest.raises(ValueError, match="proper orthogonal"):
        rotate_cartesian_hessian(input_hessian, np.diag([-1.0, 1.0, 1.0]))


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
    assert (
        result.negative_eigenvalue_tolerance_units == "same as cartesian_hessian / amu"
    )

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


def test_full_vibrational_projection_preserves_negative_index_and_rigid_invariance():
    coordinates = np.array([[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [0.1, 0.9, 0.4]])
    masses = np.array([12.0, 16.0, 1.0])
    rigid = np.linalg.qr(_mass_weighted_rigid_motions(coordinates, masses))[0]
    vibrational = np.linalg.svd(rigid.T, full_matrices=True)[2][6:].T
    expected = np.array([-2.0, 3.0, 5.0])
    mass_weighted = vibrational @ np.diag(expected) @ vibrational.T
    factors = np.repeat(np.sqrt(masses), 3)
    hessian = factors[:, None] * mass_weighted * factors[None, :]

    result = project_vibrational_hessian(coordinates, masses, hessian)
    assert result.eigenvalues == pytest.approx(expected, abs=2.0e-12)
    assert result.mass_weighted_eigenvectors.shape == (3, 3, 3)
    eigenvectors = result.mass_weighted_eigenvectors.reshape(3, -1)
    assert eigenvectors @ rigid == pytest.approx(np.zeros((3, 6)), abs=2.0e-12)

    rotation = _rotation_z(0.41)
    cartesian_rotation = np.kron(np.eye(3), rotation)
    rotated = project_vibrational_hessian(
        coordinates @ rotation.T + np.array([3.0, -2.0, 0.5]),
        masses,
        cartesian_rotation @ hessian @ cartesian_rotation.T,
    )
    assert rotated.eigenvalues == pytest.approx(expected, abs=2.0e-12)


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


def test_exact_straight_point_uses_finite_direction_neutral_baseline():
    coordinate = np.array([-1.0, 0.0, 1.0])
    curvature_vectors = np.array([[0.3, 0.4, 0.0], [0.0, 0.0, 0.0], [0.3, 0.4, 0.0]])
    modes = np.broadcast_to(np.eye(3)[:2], (3, 2, 3)).copy()
    frequencies = np.tile([400.0, 100.0], (3, 1))

    result = curvature_effective_mass_profile(
        coordinate, curvature_vectors, modes, frequencies
    )
    unregularized = curvature_effective_mass_profile(
        coordinate,
        curvature_vectors,
        modes,
        frequencies,
        straightness_tolerance_per_angstrom=0.0,
    )

    assert result.turning_length_angstrom[1] == pytest.approx(
        np.sqrt(np.mean(result.mode_turning_amplitudes_angstrom[1] ** 2)),
        rel=2.0e-15,
    )
    assert result.turning_length_derivative[1] == pytest.approx(0.0, abs=1.0e-14)
    assert result.mode_curvature_components_per_angstrom[1] == pytest.approx([0.0, 0.0])
    assert result.effective_mass_amu[1] == pytest.approx(1.0)
    assert np.all(np.isfinite(result.turning_length_angstrom))
    assert np.all(np.isfinite(unregularized.turning_length_angstrom))
    assert unregularized.turning_length_angstrom[1] == pytest.approx(
        result.turning_length_angstrom[1], rel=2.0e-15
    )


def test_near_zero_curvature_uses_same_straight_limit_as_exact_zero():
    coordinate = np.array([-1.0, 0.0, 1.0])
    exact = np.array([[0.3, 0.4, 0.0], [0.0, 0.0, 0.0], [0.3, 0.4, 0.0]])
    near = exact.copy()
    near[1, 0] = 1.0e-15
    modes = np.broadcast_to(np.eye(3)[:2], (3, 2, 3)).copy()
    frequencies = np.tile([400.0, 100.0], (3, 1))

    exact_result = curvature_effective_mass_profile(
        coordinate, exact, modes, frequencies
    )
    near_result = curvature_effective_mass_profile(coordinate, near, modes, frequencies)

    assert near_result.straightness_tolerance_per_angstrom == pytest.approx(1.0e-12)
    assert near_result.curvature_magnitudes_per_angstrom[1] == pytest.approx(1.0e-15)
    assert near_result.turning_length_angstrom == pytest.approx(
        exact_result.turning_length_angstrom, rel=2.0e-6
    )
    assert near_result.effective_mass_amu == pytest.approx(
        exact_result.effective_mass_amu, rel=2.0e-6
    )
    with pytest.raises(ValueError, match="straightness_tolerance"):
        curvature_effective_mass_profile(
            coordinate,
            exact,
            modes,
            frequencies,
            straightness_tolerance_per_angstrom=-1.0,
        )


def test_turning_length_and_mass_are_continuous_across_regularization_scale():
    coordinate = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
    tolerance = 0.2
    modes = np.broadcast_to(np.eye(3)[:2], (len(coordinate), 2, 3)).copy()
    frequencies = np.tile([400.0, 100.0], (len(coordinate), 1))

    def profile(centre_scale: float):
        scales = np.array([0.4, 0.8, centre_scale, 1.4, 2.0])
        curvature = np.column_stack(
            [scales * tolerance, np.zeros((len(coordinate), 2))]
        )
        return curvature_effective_mass_profile(
            coordinate,
            curvature,
            modes,
            frequencies,
            straightness_tolerance_per_angstrom=tolerance,
        )

    below = profile(0.999)
    above = profile(1.001)

    assert below.curvature_magnitudes_per_angstrom[2] == pytest.approx(
        0.999 * tolerance, rel=2.0e-15
    )
    assert above.curvature_magnitudes_per_angstrom[2] == pytest.approx(
        1.001 * tolerance, rel=2.0e-15
    )
    assert above.turning_length_angstrom == pytest.approx(
        below.turning_length_angstrom, rel=5.0e-4
    )
    assert above.effective_mass_amu == pytest.approx(
        below.effective_mass_amu, rel=3.0e-4
    )


def test_wholly_straight_profile_has_reference_mass_even_with_varying_modes():
    coordinate = np.array([-1.0, -0.3, 0.0, 0.4, 1.0])
    curvature = np.zeros((len(coordinate), 3))
    modes = np.broadcast_to(np.eye(3)[:2], (len(coordinate), 2, 3)).copy()
    frequencies = np.column_stack(
        [
            np.linspace(300.0, 700.0, len(coordinate)),
            np.linspace(900.0, 500.0, len(coordinate)),
        ]
    )

    result = curvature_effective_mass_profile(
        coordinate,
        curvature,
        modes,
        frequencies,
        reference_mass_amu=2.5,
    )

    assert np.all(np.isfinite(result.turning_length_angstrom))
    assert result.effective_mass_amu == pytest.approx(np.full(len(coordinate), 2.5))


def test_curved_molecular_path_composes_through_transverse_curvature_mass():
    parameter = np.linspace(-0.6, 0.6, 7)
    base = np.array([[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [0.1, 0.9, 0.4]])
    coordinates = np.asarray(
        [
            base
            + np.array(
                [
                    [0.03 * value**2, 0.0, 0.0],
                    [0.08 * value, 0.04 * value**2, 0.02 * value],
                    [-0.05 * value, 0.06 * value, -0.04 * value**2],
                ]
            )
            for value in parameter
        ]
    )
    masses = np.array([12.0, 16.0, 1.0])
    path = build_mass_scaled_path(coordinates, masses, transition_state_index=3)
    geometry = molecular_path_tangents_and_curvature(path)
    cartesian_hessian = np.diag(np.repeat(masses, 3))
    mode_results = [
        project_transverse_hessian(point, masses, cartesian_hessian, tangent)
        for point, tangent in zip(
            path.aligned_coordinates_angstrom,
            geometry.unit_tangents,
            strict=True,
        )
    ]
    eigenvectors = np.asarray(
        [result.mass_weighted_eigenvectors for result in mode_results]
    )
    flattened_modes = eigenvectors.reshape(len(parameter), 2, -1)
    projected = np.einsum(
        "pm,pmi->pi",
        np.einsum(
            "pmi,pi->pm", flattened_modes, geometry.curvature_vectors_per_angstrom
        ),
        flattened_modes,
    )
    residual = np.linalg.norm(
        geometry.curvature_vectors_per_angstrom - projected, axis=1
    )
    assert residual == pytest.approx(np.zeros(len(parameter)), abs=2.0e-11)

    profile = curvature_effective_mass_profile(
        path.coordinate_angstrom,
        geometry.curvature_vectors_per_angstrom,
        eigenvectors,
        np.full((len(parameter), 2), 500.0),
    )

    assert profile.curvature_magnitudes_per_angstrom == pytest.approx(
        geometry.curvature_magnitudes_per_angstrom, rel=2.0e-12, abs=2.0e-12
    )
    assert np.all(np.isfinite(profile.turning_length_angstrom))
    assert np.all(np.isfinite(profile.effective_mass_amu))
    assert np.all(
        (profile.effective_mass_amu > 0.0) & (profile.effective_mass_amu <= 1.0)
    )


def test_curvature_mass_rejects_incomplete_molecular_transverse_basis():
    coordinate = np.array([-1.0, 0.0, 1.0])
    curvature = np.zeros((3, 9))
    modes = np.broadcast_to(np.eye(9)[:1], (3, 1, 9)).copy()
    frequencies = np.full((3, 1), 500.0)

    with pytest.raises(ValueError, match="complete 3N-7"):
        curvature_effective_mass_profile(coordinate, curvature, modes, frequencies)


def test_duplicate_nonfinite_and_nonpositive_inputs_are_rejected():
    coordinates = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            [[1.0, 2.0, 3.0], [2.2, 2.0, 3.0]],
            [[4.0, -1.0, 2.0], [5.0, -1.0, 2.0]],
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
