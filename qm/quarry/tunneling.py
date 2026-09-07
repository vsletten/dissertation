"""Semiclassical tunnelling corrections on an explicit reaction path.

The Eckart correction in :mod:`quarry.rates` reconstructs an analytic barrier
from two heights and one imaginary frequency.  This module instead integrates
the action through a supplied vibrationally adiabatic path.  A path-dependent
effective mass implements the small-curvature correction of Liu et al.
(JACS 115, 2408, 1993), in the formulation used by Pilgrim 2021.5.

Coordinates are Angstrom, molar energies are kJ/mol, and masses are amu.  The
thermal correction is evaluated in log space so 10--20 K calculations do not
overflow even when kappa is enormous.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

_AMU_KG = 1.66053906660e-27
_ANGSTROM_M = 1.0e-10
_HBAR_J_S = 1.054571817e-34
_KJ_MOL_TO_J = 1000.0 / 6.02214076e23
_R_KJ = 8.31446261815324e-3
_C_CM_S = 2.99792458e10


@dataclass(frozen=True)
class SCTResult:
    """Small-curvature transmission receipt for one temperature."""

    temperature_k: float
    kappa: float
    log_kappa: float
    barrier_kj_mol: float
    energy_floor_kj_mol: float
    ground_probability: float
    ground_action: float
    representative_energy_kj_mol: float
    quadrature_order: int


def harmonic_turning_length_angstrom(
    imaginary_frequency_cm: float, reference_mass_amu: float = 1.0
) -> float:
    """Ground-state harmonic turning length ``sqrt(hbar/(mu*omega))``."""

    if not math.isfinite(imaginary_frequency_cm) or imaginary_frequency_cm <= 0.0:
        raise ValueError("imaginary_frequency_cm must be finite and positive")
    if not math.isfinite(reference_mass_amu) or reference_mass_amu <= 0.0:
        raise ValueError("reference_mass_amu must be finite and positive")
    omega = 2.0 * math.pi * _C_CM_S * imaginary_frequency_cm
    length_m = math.sqrt(_HBAR_J_S / (reference_mass_amu * _AMU_KG * omega))
    return length_m / _ANGSTROM_M


def small_curvature_effective_mass(
    reference_mass_amu: float,
    curvature_per_angstrom: Sequence[float] | np.ndarray,
    turning_length_angstrom: Sequence[float] | np.ndarray | float,
    turning_length_derivative: Sequence[float] | np.ndarray | float = 0.0,
) -> np.ndarray:
    """Return the Pilgrim/Liu small-curvature effective mass profile.

    This is Eq. 14 of Liu et al. (1993), ``mu_eff = mu exp[-2*k*t -
    (k*t)^2 + (dt/ds)^2]``, capped at the zero-curvature mass.  Curvature and
    turning length must use reciprocal units.  The derivative is dimensionless.
    """

    if not math.isfinite(reference_mass_amu) or reference_mass_amu <= 0.0:
        raise ValueError("reference_mass_amu must be finite and positive")
    curvature = np.asarray(curvature_per_angstrom, dtype=float)
    turning = np.broadcast_to(
        np.asarray(turning_length_angstrom, dtype=float), curvature.shape
    )
    derivative = np.broadcast_to(
        np.asarray(turning_length_derivative, dtype=float), curvature.shape
    )
    if curvature.ndim != 1 or curvature.size < 2:
        raise ValueError(
            "curvature profile must be one-dimensional with at least two points"
        )
    if not np.all(np.isfinite(curvature)) or np.any(curvature < 0.0):
        raise ValueError("curvature profile must be finite and non-negative")
    if not np.all(np.isfinite(turning)) or np.any(turning <= 0.0):
        raise ValueError("turning lengths must be finite and positive")
    if not np.all(np.isfinite(derivative)):
        raise ValueError("turning-length derivatives must be finite")
    kt = curvature * turning
    exponent = -2.0 * kt - kt * kt + derivative * derivative
    factor = np.minimum(np.exp(np.clip(exponent, -745.0, 0.0)), 1.0)
    return reference_mass_amu * factor


def _validated_path(
    path_coordinate_a: Sequence[float] | np.ndarray,
    potential_kj_mol: Sequence[float] | np.ndarray,
    effective_mass_amu: Sequence[float] | np.ndarray | float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    coordinate = np.asarray(path_coordinate_a, dtype=float)
    potential = np.asarray(potential_kj_mol, dtype=float)
    if coordinate.ndim != 1 or potential.ndim != 1 or coordinate.size != potential.size:
        raise ValueError("path coordinate and potential must be equal-length 1D arrays")
    if coordinate.size < 3:
        raise ValueError("explicit path requires at least three points")
    if not np.all(np.isfinite(coordinate)) or not np.all(np.isfinite(potential)):
        raise ValueError("path coordinate and potential must be finite")
    delta = np.diff(coordinate)
    if not (np.all(delta > 0.0) or np.all(delta < 0.0)):
        raise ValueError("path coordinate must be strictly monotonic")
    if np.all(delta < 0.0):
        coordinate = coordinate[::-1]
        potential = potential[::-1]
    mass = np.broadcast_to(
        np.asarray(effective_mass_amu, dtype=float), coordinate.shape
    ).copy()
    if not np.all(np.isfinite(mass)) or np.any(mass <= 0.0):
        raise ValueError("effective mass must be finite and positive")
    if np.all(delta < 0.0):
        mass = mass[::-1]
    return coordinate, potential, mass


def _logsumexp(values: np.ndarray) -> float:
    maximum = float(np.max(values))
    if not math.isfinite(maximum):
        return maximum
    return maximum + math.log(float(np.exp(values - maximum).sum()))


def sct_kappa(
    path_coordinate_a: Sequence[float] | np.ndarray,
    potential_kj_mol: Sequence[float] | np.ndarray,
    effective_mass_amu: Sequence[float] | np.ndarray | float,
    temperature_k: float,
    *,
    energy_floor_kj_mol: float | None = None,
    quadrature_order: int = 96,
    path_grid_size: int = 4097,
) -> SCTResult:
    """Canonical SCT transmission coefficient for an explicit path.

    The sub-barrier transmission is ``1/(1 + exp(2*theta))`` with the WKB
    action integrated over every forbidden segment.  The thermal integral is
    normalized to the classical step-barrier integral; an over-barrier term of
    one is retained.  ``energy_floor_kj_mol`` defaults to the higher endpoint,
    Pilgrim's SCT lower-energy convention.
    """

    if not math.isfinite(temperature_k) or temperature_k <= 0.0:
        raise ValueError("temperature_k must be finite and positive")
    if quadrature_order < 16:
        raise ValueError("quadrature_order must be at least 16")
    if path_grid_size < 257:
        raise ValueError("path_grid_size must be at least 257")
    coordinate, potential, mass = _validated_path(
        path_coordinate_a, potential_kj_mol, effective_mass_amu
    )
    barrier = float(np.max(potential))
    floor = (
        max(float(potential[0]), float(potential[-1]))
        if energy_floor_kj_mol is None
        else float(energy_floor_kj_mol)
    )
    if not math.isfinite(floor) or floor >= barrier:
        raise ValueError("energy floor must be finite and below the path barrier")

    dense_x = np.linspace(float(coordinate[0]), float(coordinate[-1]), path_grid_size)
    dense_v = np.interp(dense_x, coordinate, potential)
    dense_mu = np.interp(dense_x, coordinate, mass)

    def action(energy_kj_mol: float) -> float:
        deficit_j = np.maximum(dense_v - energy_kj_mol, 0.0) * _KJ_MOL_TO_J
        integrand = np.sqrt(2.0 * dense_mu * _AMU_KG * deficit_j)
        return float(np.trapezoid(integrand, dense_x * _ANGSTROM_M) / _HBAR_J_S)

    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
    energies = floor + 0.5 * (nodes + 1.0) * (barrier - floor)
    energy_weights = 0.5 * (barrier - floor) * weights
    actions = np.asarray([action(float(energy)) for energy in energies])
    log_probabilities = -np.logaddexp(0.0, 2.0 * actions)
    log_reflections = -np.logaddexp(0.0, -2.0 * actions)
    beta = 1.0 / (_R_KJ * temperature_k)
    log_i1_terms = (
        np.log(beta * energy_weights) + log_probabilities - beta * (energies - barrier)
    )
    # Pilgrim's I2 maps the same quadrature points above the barrier and
    # accounts for non-classical reflection. I3 is the analytic high-energy
    # tail. Keeping all three terms prevents a low-temperature approximation
    # from leaking into this generic API.
    log_i2_terms = (
        np.log(beta * energy_weights) + log_reflections - beta * (barrier - energies)
    )
    log_i3 = -beta * (barrier - floor)
    combined_logs = np.concatenate([log_i1_terms, log_i2_terms, np.asarray([log_i3])])
    log_kappa = _logsumexp(combined_logs)
    kappa = (
        math.exp(log_kappa)
        if log_kappa < math.log(float("1.7976931348623157e308"))
        else math.inf
    )

    ground_action = action(floor)
    ground_probability = math.exp(-float(np.logaddexp(0.0, 2.0 * ground_action)))
    representative_energies = np.concatenate(
        [energies, 2.0 * barrier - energies, np.asarray([2.0 * barrier - floor])]
    )
    representative = float(representative_energies[int(np.argmax(combined_logs))])
    return SCTResult(
        temperature_k=float(temperature_k),
        kappa=kappa,
        log_kappa=log_kappa,
        barrier_kj_mol=barrier,
        energy_floor_kj_mol=floor,
        ground_probability=ground_probability,
        ground_action=ground_action,
        representative_energy_kj_mol=representative,
        quadrature_order=quadrature_order,
    )


def corrected_rate_from_reference(
    reference_rate_s: float,
    reference_log_kappa: float,
    correction: SCTResult,
) -> float:
    """Replace a reference tunnelling correction without under/overflow."""

    if not math.isfinite(reference_rate_s) or reference_rate_s <= 0.0:
        raise ValueError("reference_rate_s must be finite and positive")
    if not math.isfinite(reference_log_kappa):
        raise ValueError("reference_log_kappa must be finite")
    log_rate = math.log(reference_rate_s) - reference_log_kappa + correction.log_kappa
    if log_rate > math.log(float("1.7976931348623157e308")):
        return math.inf
    return math.exp(log_rate)
