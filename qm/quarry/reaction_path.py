"""Fail-closed numerical primitives for a bounded molecular reaction path.

The coordinate convention is the one used by reaction-path Hamiltonians.  Cartesian
coordinates are in Angstrom and atomic masses are in unified atomic mass units (amu).
A mass-scaled Cartesian displacement is
``dq_i = sqrt(m_i / mu) dx_i`` for an explicit reference mass ``mu``; consequently
both ``q`` and the signed path coordinate ``s`` are in Angstrom.  Hessian eigenvalues
remain in the caller's Cartesian-Hessian energy/length-squared units divided by amu.

These routines do not run electronic-structure calculations.  They validate and reduce
already-computed path geometries, Hessians, energies, and frequencies.  In particular,
negative transverse modes are never made positive with ``abs`` or silently discarded.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from quarry.tunneling import small_curvature_effective_mass

_AMU_KG = 1.66053906660e-27
_ANGSTROM_M = 1.0e-10
_C_CM_S = 2.99792458e10
_H_J_S = 6.62607015e-34
_HBAR_J_S = _H_J_S / (2.0 * math.pi)
_AVOGADRO_MOL = 6.02214076e23
_WAVENUMBER_TO_KJ_MOL = _H_J_S * 299792458.0 * 100.0 * _AVOGADRO_MOL / 1000.0
_DUPLICATE_TOLERANCE_ANGSTROM = 1.0e-12


@dataclass(frozen=True)
class MassScaledPath:
    """A rigid-motion-free path in one transition-state-centred frame."""

    coordinate_angstrom: np.ndarray
    aligned_coordinates_angstrom: np.ndarray
    mass_scaled_coordinates: np.ndarray
    masses_amu: np.ndarray
    transition_state_index: int
    reference_mass_amu: float


@dataclass(frozen=True)
class PathDifferentialGeometry:
    """Unit tangents and curvature ``d tangent / ds`` on a nonuniform path."""

    coordinate_angstrom: np.ndarray
    unit_tangents: np.ndarray
    curvature_vectors_per_angstrom: np.ndarray
    curvature_magnitudes_per_angstrom: np.ndarray


@dataclass(frozen=True)
class TransverseHessianModes:
    """The exactly ``3N-7`` non-rigid, non-reaction Hessian eigenpairs.

    Eigenvectors are orthonormal in mass-weighted Cartesian coordinates and are
    returned with shape ``(3N-7, N, 3)``.  Eigenvalues use the input Hessian's
    energy/length-squared unit divided by amu.
    """

    eigenvalues: np.ndarray
    mass_weighted_eigenvectors: np.ndarray
    transverse_basis: np.ndarray
    negative_eigenvalue_tolerance: float


@dataclass(frozen=True)
class VibrationallyAdiabaticPotential:
    """Electronic, transverse-ZPE, high-level, and summed path energies."""

    electronic_energy_kj_mol: np.ndarray
    zero_point_energy_kj_mol: np.ndarray
    high_level_correction_kj_mol: np.ndarray
    potential_kj_mol: np.ndarray
    wavenumber_quantum_kj_mol_per_cm: float


@dataclass(frozen=True)
class CurvatureMassProfile:
    """Mode-resolved Liu/Pilgrim turning length and SCT effective mass."""

    curvature_magnitudes_per_angstrom: np.ndarray
    mode_curvature_components_per_angstrom: np.ndarray
    mode_turning_amplitudes_angstrom: np.ndarray
    turning_length_angstrom: np.ndarray
    turning_length_derivative: np.ndarray
    effective_mass_amu: np.ndarray
    reference_mass_amu: float


def _positive_masses(
    masses_amu: Sequence[float] | np.ndarray, atom_count: int
) -> np.ndarray:
    masses = np.asarray(masses_amu, dtype=float)
    if masses.shape != (atom_count,):
        raise ValueError(f"masses_amu must have shape ({atom_count},)")
    if not np.all(np.isfinite(masses)) or np.any(masses <= 0.0):
        raise ValueError("masses_amu must be finite and positive")
    return masses


def _mass_weighted_kabsch(
    moving: np.ndarray, reference: np.ndarray, masses_amu: np.ndarray
) -> np.ndarray:
    """Centre and optimally rotate ``moving`` onto centred ``reference``."""

    centre = np.average(moving, axis=0, weights=masses_amu)
    centred = moving - centre
    covariance = centred.T @ (masses_amu[:, None] * reference)
    left, _, right_transpose = np.linalg.svd(covariance)
    if np.linalg.det(left @ right_transpose) < 0.0:
        left[:, -1] *= -1.0
    return centred @ (left @ right_transpose)


def build_mass_scaled_path(
    coordinates_angstrom: Sequence[Sequence[Sequence[float]]] | np.ndarray,
    masses_amu: Sequence[float] | np.ndarray,
    *,
    transition_state_index: int,
    reference_mass_amu: float = 1.0,
) -> MassScaledPath:
    """Align an ordered molecular path and construct signed mass-scaled ``s``.

    Every geometry is mass-centred and fitted by a proper mass-weighted Kabsch
    rotation to the explicit transition-state geometry.  This removes arbitrary
    per-point translations and rotations without changing atom order.  ``s`` is
    strictly increasing in input order and exactly zero at ``transition_state_index``.
    Rigidly equivalent duplicate points, including non-adjacent duplicates, are
    rejected because their tangent and Hessian association would be ambiguous.
    """

    coordinates = np.asarray(coordinates_angstrom, dtype=float)
    if coordinates.ndim != 3 or coordinates.shape[2] != 3:
        raise ValueError("coordinates_angstrom must have shape (n_points, n_atoms, 3)")
    point_count, atom_count, _ = coordinates.shape
    if point_count < 3 or atom_count < 2:
        raise ValueError(
            "a molecular path requires at least three points and two atoms"
        )
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("coordinates_angstrom must be finite")
    if not 0 < transition_state_index < point_count - 1:
        raise ValueError("transition_state_index must identify an interior path point")
    if not math.isfinite(reference_mass_amu) or reference_mass_amu <= 0.0:
        raise ValueError("reference_mass_amu must be finite and positive")
    masses = _positive_masses(masses_amu, atom_count)

    transition_state = coordinates[transition_state_index]
    transition_state = transition_state - np.average(
        transition_state, axis=0, weights=masses
    )
    aligned = np.asarray(
        [
            _mass_weighted_kabsch(point, transition_state, masses)
            for point in coordinates
        ]
    )
    mass_scale = np.sqrt(masses / reference_mass_amu)[None, :, None]
    mass_scaled = (aligned * mass_scale).reshape(point_count, -1)

    pair_differences = mass_scaled[:, None, :] - mass_scaled[None, :, :]
    pair_distances = np.linalg.norm(pair_differences, axis=2)
    duplicate_pairs = np.argwhere(
        np.triu(pair_distances <= _DUPLICATE_TOLERANCE_ANGSTROM, k=1)
    )
    if duplicate_pairs.size:
        first, second = duplicate_pairs[0]
        raise ValueError(
            "path contains rigidly equivalent duplicate points "
            f"at indices {first} and {second}"
        )

    segment_lengths = np.linalg.norm(np.diff(mass_scaled, axis=0), axis=1)
    cumulative = np.concatenate(([0.0], np.cumsum(segment_lengths)))
    coordinate = cumulative - cumulative[transition_state_index]
    return MassScaledPath(
        coordinate_angstrom=coordinate,
        aligned_coordinates_angstrom=aligned,
        mass_scaled_coordinates=mass_scaled,
        masses_amu=masses.copy(),
        transition_state_index=transition_state_index,
        reference_mass_amu=float(reference_mass_amu),
    )


def _strict_coordinate(
    coordinate_angstrom: Sequence[float] | np.ndarray, *, minimum_points: int = 3
) -> np.ndarray:
    coordinate = np.asarray(coordinate_angstrom, dtype=float)
    if coordinate.ndim != 1 or coordinate.size < minimum_points:
        raise ValueError(
            "coordinate_angstrom must be one-dimensional with at least "
            f"{minimum_points} points"
        )
    if not np.all(np.isfinite(coordinate)):
        raise ValueError("coordinate_angstrom must be finite")
    if not np.all(np.diff(coordinate) > 0.0):
        raise ValueError(
            "coordinate_angstrom must be strictly increasing without duplicates"
        )
    return coordinate


def path_tangents_and_curvature(
    coordinate_angstrom: Sequence[float] | np.ndarray,
    mass_scaled_coordinates: Sequence[Sequence[float]] | np.ndarray,
) -> PathDifferentialGeometry:
    """Differentiate a mass-scaled path on its nonuniform ``s`` grid.

    Second-order nonuniform finite differences are used at both interior and edge
    points.  Curvature is explicitly projected normal to the local normalized
    tangent, preventing finite-difference normalization noise from masquerading as
    physical curvature.
    """

    coordinate = _strict_coordinate(coordinate_angstrom)
    points = np.asarray(mass_scaled_coordinates, dtype=float)
    if points.ndim != 2 or points.shape[0] != coordinate.size or points.shape[1] < 2:
        raise ValueError(
            "mass_scaled_coordinates must have shape (n_points, n_dimensions)"
        )
    if not np.all(np.isfinite(points)):
        raise ValueError("mass_scaled_coordinates must be finite")
    if np.any(np.linalg.norm(np.diff(points, axis=0), axis=1) <= 1.0e-14):
        raise ValueError("mass_scaled_coordinates contains adjacent duplicate points")

    derivatives = np.gradient(points, coordinate, axis=0, edge_order=2)
    derivative_norms = np.linalg.norm(derivatives, axis=1)
    if np.any(~np.isfinite(derivative_norms)) or np.any(derivative_norms <= 1.0e-14):
        raise ValueError("path has an undefined local tangent")
    tangents = derivatives / derivative_norms[:, None]
    curvature = np.gradient(tangents, coordinate, axis=0, edge_order=2)
    curvature -= np.einsum("ij,ij->i", curvature, tangents)[:, None] * tangents
    magnitudes = np.linalg.norm(curvature, axis=1)
    return PathDifferentialGeometry(
        coordinate_angstrom=coordinate.copy(),
        unit_tangents=tangents,
        curvature_vectors_per_angstrom=curvature,
        curvature_magnitudes_per_angstrom=magnitudes,
    )


def _rigid_motion_basis(
    coordinates_angstrom: np.ndarray, masses_amu: np.ndarray
) -> np.ndarray:
    centred = coordinates_angstrom - np.average(
        coordinates_angstrom, axis=0, weights=masses_amu
    )
    roots = np.sqrt(masses_amu)
    columns: list[np.ndarray] = []
    for axis in np.eye(3):
        columns.append((roots[:, None] * axis).reshape(-1))
    for axis in np.eye(3):
        columns.append((roots[:, None] * np.cross(axis, centred)).reshape(-1))
    rigid = np.column_stack(columns)
    left, singular_values, _ = np.linalg.svd(rigid, full_matrices=False)
    tolerance = max(rigid.shape) * np.finfo(float).eps * singular_values[0]
    rank = int(np.count_nonzero(singular_values > tolerance))
    if rank != 6:
        raise ValueError(
            "coordinates must describe a nonlinear molecule with six rigid modes"
        )
    return left[:, :6]


def project_transverse_hessian(
    coordinates_angstrom: Sequence[Sequence[float]] | np.ndarray,
    masses_amu: Sequence[float] | np.ndarray,
    cartesian_hessian: Sequence[Sequence[float]] | np.ndarray,
    reaction_tangent_mass_scaled: Sequence[Sequence[float]] | np.ndarray,
    *,
    negative_eigenvalue_tolerance: float = 1.0e-8,
) -> TransverseHessianModes:
    """Project a Cartesian Hessian into a nonlinear path's transverse space.

    The supplied reaction tangent must use the mass-scaled Cartesian convention
    of :func:`build_mass_scaled_path`.  Translations, rotations, and the rigid-free
    tangent are removed by an orthonormal null-space projection before diagonalizing
    ``M^-1/2 H M^-1/2``.  A nonlinear ``N``-atom molecule therefore returns exactly
    ``3N-7`` eigenpairs.  Any eigenvalue below ``-negative_eigenvalue_tolerance`` is
    a hard error; small values inside the declared tolerance are returned unchanged.
    """

    coordinates = np.asarray(coordinates_angstrom, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3 or coordinates.shape[0] < 3:
        raise ValueError(
            "coordinates_angstrom must have shape (n_atoms, 3), n_atoms >= 3"
        )
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("coordinates_angstrom must be finite")
    atom_count = coordinates.shape[0]
    masses = _positive_masses(masses_amu, atom_count)
    if (
        not math.isfinite(negative_eigenvalue_tolerance)
        or negative_eigenvalue_tolerance < 0.0
    ):
        raise ValueError(
            "negative_eigenvalue_tolerance must be finite and non-negative"
        )

    hessian = np.asarray(cartesian_hessian, dtype=float)
    if hessian.shape == (atom_count, 3, atom_count, 3):
        hessian = hessian.reshape(3 * atom_count, 3 * atom_count)
    if hessian.shape != (3 * atom_count, 3 * atom_count):
        raise ValueError("cartesian_hessian must have shape (3N, 3N) or (N, 3, N, 3)")
    if not np.all(np.isfinite(hessian)):
        raise ValueError("cartesian_hessian must be finite")
    symmetry_scale = max(1.0, float(np.max(np.abs(hessian))))
    if not np.allclose(hessian, hessian.T, rtol=1.0e-10, atol=1.0e-12 * symmetry_scale):
        raise ValueError("cartesian_hessian must be symmetric")
    hessian = 0.5 * (hessian + hessian.T)

    rigid = _rigid_motion_basis(coordinates, masses)
    tangent = np.asarray(reaction_tangent_mass_scaled, dtype=float)
    if tangent.shape == (atom_count, 3):
        tangent = tangent.reshape(-1)
    if tangent.shape != (3 * atom_count,) or not np.all(np.isfinite(tangent)):
        raise ValueError(
            "reaction_tangent_mass_scaled must be a finite (N, 3) or (3N,) vector"
        )
    tangent = tangent - rigid @ (rigid.T @ tangent)
    tangent_norm = float(np.linalg.norm(tangent))
    if tangent_norm <= 1.0e-12:
        raise ValueError("reaction tangent has no non-rigid component")
    tangent /= tangent_norm

    removed = np.column_stack([rigid, tangent])
    _, singular_values, right_transpose = np.linalg.svd(removed.T, full_matrices=True)
    tolerance = max(removed.T.shape) * np.finfo(float).eps * singular_values[0]
    if int(np.count_nonzero(singular_values > tolerance)) != 7:
        raise ValueError("reaction tangent is not independent of rigid motion")
    transverse_basis = right_transpose[7:].T
    expected_modes = 3 * atom_count - 7
    if transverse_basis.shape != (3 * atom_count, expected_modes):
        raise RuntimeError("failed to construct the exact 3N-7 transverse space")

    mass_factors = np.repeat(np.sqrt(masses), 3)
    mass_weighted_hessian = hessian / (mass_factors[:, None] * mass_factors[None, :])
    projected = transverse_basis.T @ mass_weighted_hessian @ transverse_basis
    eigenvalues, eigenvectors_in_basis = np.linalg.eigh(0.5 * (projected + projected.T))
    if eigenvalues[0] < -negative_eigenvalue_tolerance:
        raise ValueError(
            "significant negative transverse eigenvalue: "
            f"{eigenvalues[0]:.12g} < -{negative_eigenvalue_tolerance:.12g}"
        )
    eigenvectors = (transverse_basis @ eigenvectors_in_basis).T.reshape(
        expected_modes, atom_count, 3
    )
    return TransverseHessianModes(
        eigenvalues=eigenvalues,
        mass_weighted_eigenvectors=eigenvectors,
        transverse_basis=transverse_basis,
        negative_eigenvalue_tolerance=float(negative_eigenvalue_tolerance),
    )


def vibrationally_adiabatic_potential(
    electronic_energy_kj_mol: Sequence[float] | np.ndarray,
    transverse_frequencies_cm: Sequence[Sequence[float]] | np.ndarray,
    *,
    high_level_correction_kj_mol: Sequence[float] | np.ndarray | None = None,
) -> VibrationallyAdiabaticPotential:
    """Add transverse harmonic ZPE and an optional same-length energy correction.

    Frequencies are positive wavenumbers in cm^-1.  The conversion uses the exact
    SI definitions of Planck's constant, the speed of light, and Avogadro's number:
    one wavenumber quantum is ``h c N_A`` after cm^-1 -> m^-1 and J -> kJ.
    """

    electronic = np.asarray(electronic_energy_kj_mol, dtype=float)
    frequencies = np.asarray(transverse_frequencies_cm, dtype=float)
    if electronic.ndim != 1 or electronic.size == 0:
        raise ValueError("electronic_energy_kj_mol must be a non-empty 1D array")
    if frequencies.ndim != 2 or frequencies.shape[0] != electronic.size:
        raise ValueError(
            "transverse_frequencies_cm must have shape (n_points, n_modes)"
        )
    if frequencies.shape[1] == 0:
        raise ValueError("at least one transverse frequency is required")
    if not np.all(np.isfinite(electronic)):
        raise ValueError("electronic energies must be finite")
    if not np.all(np.isfinite(frequencies)) or np.any(frequencies <= 0.0):
        raise ValueError("transverse frequencies must be finite and positive")

    if high_level_correction_kj_mol is None:
        correction = np.zeros_like(electronic)
    else:
        correction = np.asarray(high_level_correction_kj_mol, dtype=float)
        if correction.shape != electronic.shape:
            raise ValueError(
                "high-level correction must have the same length as energies"
            )
        if not np.all(np.isfinite(correction)):
            raise ValueError("high-level correction must be finite")
    zero_point = 0.5 * _WAVENUMBER_TO_KJ_MOL * frequencies.sum(axis=1)
    return VibrationallyAdiabaticPotential(
        electronic_energy_kj_mol=electronic.copy(),
        zero_point_energy_kj_mol=zero_point,
        high_level_correction_kj_mol=correction.copy(),
        potential_kj_mol=electronic + zero_point + correction,
        wavenumber_quantum_kj_mol_per_cm=_WAVENUMBER_TO_KJ_MOL,
    )


def curvature_effective_mass_profile(
    coordinate_angstrom: Sequence[float] | np.ndarray,
    curvature_vectors_per_angstrom: Sequence[Sequence[float]] | np.ndarray,
    transverse_eigenvectors: np.ndarray,
    transverse_frequencies_cm: Sequence[Sequence[float]] | np.ndarray,
    *,
    reference_mass_amu: float = 1.0,
) -> CurvatureMassProfile:
    """Derive the mode-resolved Liu/Pilgrim SCT mass profile.

    At each point, the curvature vector is resolved into the complete orthonormal
    transverse normal-mode basis.  With ``B_k = e_k . curvature``,
    ``kappa = sqrt(sum_k B_k**2)``, and mode turning amplitudes
    ``t_k = sqrt(hbar/(mu*omega_k))``, the Pilgrim/Liu turning length is
    ``tbar = sqrt(kappa) * sum_k((B_k/t_k**2)**2)**(-1/4)``.  This expression is
    invariant to normal-mode signs and permutations, and to rotations within an
    exactly degenerate mode subspace.  The nonuniform derivative ``dtbar/ds`` and
    curvature are then passed unchanged to
    :func:`small_curvature_effective_mass`.

    At an exactly straight point (``kappa == 0``) the formula's directional limit
    is undefined.  Isolated straight points therefore receive a coordinate-linear
    interpolation of neighboring bent-point values.  A wholly straight path uses
    one constant positive diagnostic amplitude, so ``dtbar/ds`` is exactly zero and
    cannot manufacture an effective-mass correction.
    """

    coordinate = _strict_coordinate(coordinate_angstrom)
    curvature = np.asarray(curvature_vectors_per_angstrom, dtype=float)
    if curvature.ndim != 2 or curvature.shape[0] != coordinate.size:
        raise ValueError("curvature vectors must have shape (n_points, n_dimensions)")
    if not np.all(np.isfinite(curvature)):
        raise ValueError("curvature vectors must be finite")
    if not math.isfinite(reference_mass_amu) or reference_mass_amu <= 0.0:
        raise ValueError("reference_mass_amu must be finite and positive")

    modes = np.asarray(transverse_eigenvectors, dtype=float)
    if modes.ndim < 3 or modes.shape[0] != coordinate.size:
        raise ValueError("transverse_eigenvectors must start with (n_points, n_modes)")
    modes = modes.reshape(modes.shape[0], modes.shape[1], -1)
    if modes.shape[2] != curvature.shape[1] or modes.shape[1] == 0:
        raise ValueError("transverse eigenvectors and curvature dimensions differ")
    if not np.all(np.isfinite(modes)):
        raise ValueError("transverse eigenvectors must be finite")
    gram = np.einsum("pmi,pni->pmn", modes, modes)
    identity = np.broadcast_to(np.eye(modes.shape[1]), gram.shape)
    if not np.allclose(gram, identity, rtol=1.0e-9, atol=1.0e-10):
        raise ValueError("transverse eigenvectors must be orthonormal at every point")

    frequencies = np.asarray(transverse_frequencies_cm, dtype=float)
    if frequencies.shape != modes.shape[:2]:
        raise ValueError("transverse frequencies must match path points and modes")
    if not np.all(np.isfinite(frequencies)) or np.any(frequencies <= 0.0):
        raise ValueError("transverse frequencies must be finite and positive")

    vector_magnitudes = np.linalg.norm(curvature, axis=1)
    components = np.einsum("pmi,pi->pm", modes, curvature)
    kappa = np.linalg.norm(components, axis=1)
    bent = vector_magnitudes > 0.0
    projected_curvature = np.einsum("pm,pmi->pi", components, modes)
    projection_residuals = np.linalg.norm(curvature - projected_curvature, axis=1)
    if np.any(projection_residuals[bent] > 1.0e-7 * vector_magnitudes[bent]):
        raise ValueError(
            "curvature vector is not contained in the transverse mode space"
        )

    omega = 2.0 * math.pi * _C_CM_S * frequencies
    amplitudes = (
        np.sqrt(_HBAR_J_S / (reference_mass_amu * _AMU_KG * omega)) / _ANGSTROM_M
    )
    turning = np.empty(coordinate.size, dtype=float)
    if np.any(bent):
        inverse_quartic_sum = np.sum(
            (components[bent] / amplitudes[bent] ** 2) ** 2, axis=1
        )
        if np.any(~np.isfinite(inverse_quartic_sum)) or np.any(
            inverse_quartic_sum <= 0.0
        ):
            raise ValueError("bent-point turning-length denominator must be positive")
        turning[bent] = np.sqrt(kappa[bent]) * inverse_quartic_sum ** (-0.25)
        turning[~bent] = np.interp(coordinate[~bent], coordinate[bent], turning[bent])
    else:
        turning.fill(float(np.sqrt(np.mean(amplitudes**2))))
    derivative = np.gradient(turning, coordinate, edge_order=2)
    effective_mass = small_curvature_effective_mass(
        reference_mass_amu,
        kappa,
        turning,
        derivative,
    )
    return CurvatureMassProfile(
        curvature_magnitudes_per_angstrom=kappa,
        mode_curvature_components_per_angstrom=components,
        mode_turning_amplitudes_angstrom=amplitudes,
        turning_length_angstrom=turning,
        turning_length_derivative=derivative,
        effective_mass_amu=effective_mass,
        reference_mass_amu=float(reference_mass_amu),
    )
