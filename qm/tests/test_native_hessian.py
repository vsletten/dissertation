"""Deterministic gates for canonical native Hessian evaluation."""

from __future__ import annotations

import numpy as np
import pytest

from quarry import native_hessian
from quarry.clusters import water
from quarry.pipeline import BOHR_TO_ANGSTROM, DftSettings


def test_native_pyscf_hessian_layout_n6_is_canonicalized_by_axis_semantics():
    atom_count = 6
    raw = np.arange((3 * atom_count) ** 2, dtype=float).reshape(3 * atom_count, -1)
    canonical = raw + raw.T
    native = canonical.reshape(atom_count, 3, atom_count, 3).transpose(0, 2, 1, 3)

    observed = native_hessian.canonicalize_pyscf_hessian(native, atom_count)

    assert observed == pytest.approx(canonical)
    assert not np.array_equal(native.reshape(3 * atom_count, 3 * atom_count), canonical)


def test_native_pyscf_hessian_accepts_observed_backend_symmetry_noise():
    atom_count = 3
    canonical = np.eye(3 * atom_count)
    canonical[0, 1] = 1.1e-9
    native = canonical.reshape(atom_count, 3, atom_count, 3).transpose(0, 2, 1, 3)

    observed = native_hessian.canonicalize_pyscf_hessian(native, atom_count)

    assert observed == pytest.approx(0.5 * (canonical + canonical.T))


def test_native_pyscf_hessian_rejects_shape_nonfinite_and_material_asymmetry():
    with pytest.raises(ValueError, match="shape"):
        native_hessian.canonicalize_pyscf_hessian(np.eye(6), 2)
    native = np.zeros((2, 2, 3, 3))
    native[0, 0, 0, 0] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        native_hessian.canonicalize_pyscf_hessian(native, 2)
    native = np.zeros((2, 2, 3, 3))
    native[0, 1, 0, 1] = 1.0
    with pytest.raises(ValueError, match="not symmetric"):
        native_hessian.canonicalize_pyscf_hessian(native, 2)


def test_native_cartesian_hessian_records_gradient_fmax_and_backend(monkeypatch):
    cluster = water()
    native = np.zeros((3, 3, 3, 3))
    gradient_value = np.array([[0.01, 0.0, 0.0], [0.0, -0.02, 0.0], [0.0, 0.0, 0.0]])

    class FakeGradient:
        def kernel(self):
            return gradient_value

    class FakeScf:
        def nuc_grad_method(self):
            return FakeGradient()

    actual = DftSettings(xc="b3lyp", basis="sto-3g", use_gpu=False)
    monkeypatch.setattr(
        native_hessian,
        "_execute_native_hessian",
        lambda cluster, settings: (FakeScf(), -75.0, native, actual, True),
    )
    requested = DftSettings(xc="b3lyp", basis="sto-3g", use_gpu=True)

    result = native_hessian.native_cartesian_hessian(cluster, requested)

    assert result.gradient_hartree_per_bohr == pytest.approx(gradient_value)
    assert result.physical_fmax_ev_per_angstrom == pytest.approx(
        0.02 * native_hessian.HARTREE_TO_EV / BOHR_TO_ANGSTROM
    )
    assert result.requested_backend == "gpu4pyscf"
    assert result.actual_backend == "pyscf"
    assert result.gpu_fallback_used is True


@pytest.mark.parametrize("reported", [float("nan"), float("inf"), -1.0, 0.0])
def test_native_hessian_result_rejects_invalid_or_inconsistent_fmax(reported: float):
    gradient = np.full((1, 3), 0.01)
    with pytest.raises(ValueError, match="fmax"):
        native_hessian.NativeHessianResult(
            electronic_hartree=-1.0,
            gradient_hartree_per_bohr=gradient,
            physical_fmax_ev_per_angstrom=reported,
            cartesian_hessian_hartree_per_bohr2=np.eye(3),
            requested_backend="pyscf",
            actual_backend="pyscf",
            gpu_fallback_used=False,
            geometry_fingerprint="geometry",
            settings_fingerprint="settings",
        )
