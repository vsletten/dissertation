"""Fresh native PySCF Hessians for reaction-path qualification.

This module is deliberately separate from :mod:`quarry.pipeline`: historical
one-shot verifiers hash-pin that module byte-for-byte. PySCF's native Hessian
layout is converted once here before reaction-path code sees it.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, replace
from typing import Any

import numpy as np

from quarry.clusters import Cluster
from quarry.pipeline import (
    BOHR_TO_ANGSTROM,
    DftSettings,
    _gradient_method,
    _hessian_method,
    _is_gpu_hessian_contiguity_assertion,
    _make_scf,
    build_mol,
    frequency_geometry_fingerprint,
    frequency_settings_fingerprint,
)

HARTREE_TO_EV = 27.211386245988
HESSIAN_SYMMETRY_ABSOLUTE_MAX = 1.0e-8
HESSIAN_SYMMETRY_SPECTRAL_RELATIVE_MAX = 1.0e-8


def _immutable_f64(array: np.ndarray) -> np.ndarray:
    """Return an irrevocably read-only float64 array backed by immutable bytes."""

    contiguous = np.ascontiguousarray(array, dtype="<f8")
    return np.frombuffer(contiguous.tobytes(), dtype="<f8").reshape(contiguous.shape)


@dataclass(frozen=True)
class HessianSymmetryMetrics:
    """Raw Cartesian-Hessian skew metrics evaluated before symmetrization."""

    maximum_absolute_asymmetry: float
    matrix_absolute_maximum: float
    symmetric_spectral_norm: float
    asymmetry_spectral_norm: float
    spectral_relative_asymmetry: float
    absolute_limit: float
    spectral_relative_limit: float
    accepted: bool


def cartesian_hessian_symmetry_metrics(hessian: Any) -> HessianSymmetryMetrics:
    """Measure global absolute and scale-relative skew without hiding low scales."""

    matrix = np.asarray(hessian, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] < 1 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Cartesian Hessian must be a non-empty square matrix")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("Cartesian Hessian contains non-finite values")
    symmetric = 0.5 * (matrix + matrix.T)
    asymmetry = matrix - matrix.T
    maximum_absolute_asymmetry = float(np.max(np.abs(asymmetry)))
    matrix_absolute_maximum = float(np.max(np.abs(matrix)))
    symmetric_spectral_norm = float(np.linalg.norm(symmetric, ord=2))
    asymmetry_spectral_norm = float(np.linalg.norm(asymmetry, ord=2))
    if symmetric_spectral_norm == 0.0:
        spectral_relative_asymmetry = (
            0.0 if asymmetry_spectral_norm == 0.0 else math.inf
        )
    else:
        spectral_relative_asymmetry = asymmetry_spectral_norm / symmetric_spectral_norm
    accepted = (
        maximum_absolute_asymmetry < HESSIAN_SYMMETRY_ABSOLUTE_MAX
        and spectral_relative_asymmetry < HESSIAN_SYMMETRY_SPECTRAL_RELATIVE_MAX
    )
    return HessianSymmetryMetrics(
        maximum_absolute_asymmetry=maximum_absolute_asymmetry,
        matrix_absolute_maximum=matrix_absolute_maximum,
        symmetric_spectral_norm=symmetric_spectral_norm,
        asymmetry_spectral_norm=asymmetry_spectral_norm,
        spectral_relative_asymmetry=spectral_relative_asymmetry,
        absolute_limit=HESSIAN_SYMMETRY_ABSOLUTE_MAX,
        spectral_relative_limit=HESSIAN_SYMMETRY_SPECTRAL_RELATIVE_MAX,
        accepted=accepted,
    )


def _symmetry_error(metrics: HessianSymmetryMetrics) -> str:
    return (
        "Cartesian Hessian is not symmetric: "
        f"max_abs_asymmetry={metrics.maximum_absolute_asymmetry:.12g}, "
        f"matrix_abs_max={metrics.matrix_absolute_maximum:.12g}, "
        f"absolute_limit={metrics.absolute_limit:.12g}, "
        f"spectral_relative_asymmetry={metrics.spectral_relative_asymmetry:.12g}, "
        f"spectral_relative_limit={metrics.spectral_relative_limit:.12g}"
    )


@dataclass(frozen=True)
class NativeHessianResult:
    """One fresh SCF/gradient/Hessian evaluation in canonical Cartesian layout."""

    electronic_hartree: float
    gradient_hartree_per_bohr: np.ndarray
    physical_fmax_ev_per_angstrom: float
    cartesian_hessian_hartree_per_bohr2: np.ndarray
    requested_backend: str
    actual_backend: str
    gpu_fallback_used: bool
    geometry_fingerprint: str
    settings_fingerprint: str

    def __post_init__(self) -> None:
        gradient = np.asarray(self.gradient_hartree_per_bohr, dtype=float)
        hessian = np.asarray(self.cartesian_hessian_hartree_per_bohr2, dtype=float)
        if gradient.ndim != 2 or gradient.shape[1] != 3 or gradient.shape[0] < 1:
            raise ValueError("physical gradient must have shape (N, 3)")
        if not np.all(np.isfinite(gradient)):
            raise ValueError("physical gradient contains non-finite values")
        if hessian.shape != (3 * gradient.shape[0], 3 * gradient.shape[0]):
            raise ValueError(
                "canonical Cartesian Hessian shape does not match gradient"
            )
        if not np.all(np.isfinite(hessian)):
            raise ValueError("canonical Cartesian Hessian contains non-finite values")
        symmetry = cartesian_hessian_symmetry_metrics(hessian)
        if not symmetry.accepted:
            raise ValueError(_symmetry_error(symmetry))
        expected_fmax = float(np.max(np.linalg.norm(gradient, axis=1)))
        expected_fmax *= HARTREE_TO_EV / BOHR_TO_ANGSTROM
        if (
            not math.isfinite(self.electronic_hartree)
            or not math.isfinite(self.physical_fmax_ev_per_angstrom)
            or self.physical_fmax_ev_per_angstrom < 0.0
            or not math.isclose(
                self.physical_fmax_ev_per_angstrom,
                expected_fmax,
                rel_tol=1.0e-12,
                abs_tol=1.0e-12,
            )
        ):
            raise ValueError("physical fmax is invalid or inconsistent with gradient")
        allowed_backends = {"pyscf", "gpu4pyscf"}
        if (
            self.requested_backend not in allowed_backends
            or self.actual_backend not in allowed_backends
        ):
            raise ValueError("native Hessian backend identity is invalid")
        valid_provenance = {
            ("pyscf", "pyscf", False),
            ("gpu4pyscf", "gpu4pyscf", False),
            ("gpu4pyscf", "pyscf", True),
        }
        if (
            self.requested_backend,
            self.actual_backend,
            self.gpu_fallback_used,
        ) not in valid_provenance:
            raise ValueError("native Hessian fallback provenance is inconsistent")
        if not self.geometry_fingerprint or not self.settings_fingerprint:
            raise ValueError("native Hessian fingerprints must be non-empty")
        gradient = _immutable_f64(gradient)
        canonical_hessian = _immutable_f64(0.5 * (hessian + hessian.T))
        object.__setattr__(self, "gradient_hartree_per_bohr", gradient)
        object.__setattr__(
            self,
            "cartesian_hessian_hartree_per_bohr2",
            canonical_hessian,
        )


def canonicalize_pyscf_hessian(hessian: Any, atom_count: int) -> np.ndarray:
    """Convert ``(atom_i, atom_j, xyz_i, xyz_j)`` to canonical ``(3N,3N)``."""

    native = np.asarray(
        hessian.get() if hasattr(hessian, "get") else hessian,
        dtype=float,
    )
    if native.shape != (atom_count, atom_count, 3, 3):
        raise ValueError(
            "PySCF Hessian must have shape "
            f"({atom_count}, {atom_count}, 3, 3), got {native.shape}"
        )
    if not np.all(np.isfinite(native)):
        raise ValueError("PySCF Hessian contains non-finite values")
    canonical = native.transpose(0, 2, 1, 3).reshape(3 * atom_count, 3 * atom_count)
    symmetry = cartesian_hessian_symmetry_metrics(canonical)
    if not symmetry.accepted:
        raise ValueError(
            "PySCF Hessian is not symmetric in canonical Cartesian layout: "
            + _symmetry_error(symmetry)
        )
    return 0.5 * (canonical + canonical.T)


def _execute_native_hessian(
    cluster: Cluster, settings: DftSettings
) -> tuple[Any, float, Any, DftSettings, bool]:
    mf = _make_scf(build_mol(cluster, settings), settings)
    electronic = mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"SCF did not converge for {cluster.name}")
    try:
        return (
            mf,
            float(electronic),
            _hessian_method(mf, settings).kernel(),
            settings,
            False,
        )
    except AssertionError as gpu_error:
        if not settings.use_gpu or not _is_gpu_hessian_contiguity_assertion(gpu_error):
            raise
        warnings.warn(
            "GPU4PySCF Hessian assertion; retrying this Hessian on CPU",
            RuntimeWarning,
            stacklevel=2,
        )
        actual_settings = replace(settings, use_gpu=False)
        mf = _make_scf(build_mol(cluster, actual_settings), actual_settings)
        electronic = mf.kernel()
        if not mf.converged:
            raise RuntimeError(
                f"CPU fallback SCF did not converge for {cluster.name}"
            ) from gpu_error
        return (
            mf,
            float(electronic),
            _hessian_method(mf, actual_settings).kernel(),
            actual_settings,
            True,
        )


def native_cartesian_hessian(
    cluster: Cluster, settings: DftSettings
) -> NativeHessianResult:
    """Evaluate a fresh physical gradient and canonical analytic Hessian."""

    mf, electronic, native, actual_settings, fallback_used = _execute_native_hessian(
        cluster, settings
    )
    gradient_value = _gradient_method(mf, actual_settings).kernel()
    if hasattr(gradient_value, "get"):
        gradient_value = gradient_value.get()
    gradient = np.asarray(gradient_value, dtype=float)
    if gradient.shape != (len(cluster.symbols), 3):
        raise ValueError(
            "physical gradient must have shape "
            f"({len(cluster.symbols)}, 3), got {gradient.shape}"
        )
    if not np.all(np.isfinite(gradient)):
        raise ValueError("physical gradient contains non-finite values")
    fmax = float(np.max(np.linalg.norm(gradient, axis=1)))
    fmax *= HARTREE_TO_EV / BOHR_TO_ANGSTROM
    return NativeHessianResult(
        electronic_hartree=electronic,
        gradient_hartree_per_bohr=gradient,
        physical_fmax_ev_per_angstrom=fmax,
        cartesian_hessian_hartree_per_bohr2=canonicalize_pyscf_hessian(
            native, len(cluster.symbols)
        ),
        requested_backend="gpu4pyscf" if settings.use_gpu else "pyscf",
        actual_backend="gpu4pyscf" if actual_settings.use_gpu else "pyscf",
        gpu_fallback_used=fallback_used,
        geometry_fingerprint=frequency_geometry_fingerprint(cluster),
        settings_fingerprint=frequency_settings_fingerprint(settings),
    )
