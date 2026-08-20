"""Gates for the PySCF-backed pipeline rungs.

Everything here is CPU-only and tiny (H2O at HF/STO-3G) so the default
gate stays green in any session. The thermo cross-check gates quarry's
own formulas against pyscf's audited harmonic implementation.
"""

import numpy as np
import pytest

from quarry.clusters import water
from quarry.pipeline import (
    DftSettings,
    FrequencyResult,
    build_mol,
    energy,
    frequencies,
    gradient,
    optimize,
)
from quarry.rates import thermo_from_frequencies

CHEAP = DftSettings(xc="hf", basis="sto-3g")


@pytest.fixture(scope="module")
def opt_water():
    return optimize(water(), CHEAP)


@pytest.fixture(scope="module")
def freq_water(opt_water) -> FrequencyResult:
    return frequencies(opt_water, CHEAP)


class TestElectronicStructure:
    def test_water_hf_sto3g_energy(self):
        # HF/STO-3G water sits near -74.96 Eh for any sane geometry.
        e = energy(water(), CHEAP)
        assert -75.1 < e < -74.8

    def test_gradient_shape_and_finiteness(self):
        g = gradient(water(), CHEAP)
        assert g.shape == (3, 3)
        assert np.all(np.isfinite(g))

    def test_mol_carries_charge_and_spin(self):
        from quarry.clusters import hydronium

        mol = build_mol(hydronium(), CHEAP)
        assert mol.charge == 1
        assert mol.spin == 0

    def test_density_fit_energy_close_to_exact(self):
        from dataclasses import replace

        e_exact = energy(water(), CHEAP)
        e_df = energy(water(), replace(CHEAP, density_fit=True))
        # RI error is well under 0.1 mHa/atom for sane aux bases.
        assert e_df == pytest.approx(e_exact, abs=1e-3)
        assert e_df != e_exact  # DF path actually taken


class TestOptimization:
    def test_energy_decreases(self, opt_water):
        assert energy(opt_water, CHEAP) < energy(water(), CHEAP)

    def test_oh_bond_at_hf_sto3g_value(self, opt_water):
        # Literature HF/STO-3G r(O-H) = 0.989 A.
        o, h = opt_water.coords[0], opt_water.coords[1]
        assert np.linalg.norm(o - h) == pytest.approx(0.989, abs=0.02)

    def test_gradient_near_zero_at_minimum(self, opt_water):
        g = gradient(opt_water, CHEAP)
        assert np.abs(g).max() < 1e-3  # Hartree/Bohr, geomeTRIC default conv


class TestFrequencies:
    def test_three_real_modes_no_imaginary(self, freq_water):
        assert freq_water.n_imaginary == 0
        assert freq_water.frequencies_cm.shape == (3,)

    def test_mode_range_is_physical(self, freq_water):
        # HF/STO-3G water: bend ~2170, stretches ~4140/4390 cm^-1.
        assert freq_water.frequencies_cm[0] == pytest.approx(2170.0, abs=150.0)
        assert freq_water.frequencies_cm[-1] == pytest.approx(4390.0, abs=250.0)

    def test_rotational_temperatures_nonlinear(self, freq_water):
        assert freq_water.linear is False
        assert len(freq_water.rotational_temperatures_k) == 3
        # Water's rotational temperatures are tens of K.
        assert all(5.0 < t < 60.0 for t in freq_water.rotational_temperatures_k)


class TestThermoCrossCheck:
    """quarry.rates vs pyscf.hessian.thermo on the same frequencies."""

    @staticmethod
    def _pyscf_thermo(opt_water, freq_water):
        from pyscf.hessian import thermo as pt

        from quarry.pipeline import _make_scf

        mf = _make_scf(build_mol(opt_water, CHEAP), CHEAP)
        mf.kernel()
        hess = mf.Hessian().kernel()
        info = pt.harmonic_analysis(mf.mol, hess)
        return pt.thermo(mf, info["freq_au"], 298.15, 101325.0)

    def test_zpe_and_entropy_match_pyscf(self, opt_water, freq_water):
        ref = self._pyscf_thermo(opt_water, freq_water)

        def val(key):
            v = ref[key]
            return float(v[0]) if isinstance(v, tuple | list | np.ndarray) else float(v)

        ours = thermo_from_frequencies(
            0.0,
            freq_water.frequencies_cm,
            298.15,
            molar_mass_kg=freq_water.molar_mass_kg,
            rotational_temperatures_k=list(freq_water.rotational_temperatures_k),
            symmetry_number=2,  # pyscf detects water's C2v sigma=2
            pressure_pa=101325.0,
            qrrho_cutoff_cm=1e-9,  # pure harmonic for the comparison
        )
        hartree_to_kj = 2625.4996394798254
        assert ours.zpe_kj == pytest.approx(val("ZPE") * hartree_to_kj, rel=1e-6)
        assert ours.entropy_kj_per_k == pytest.approx(
            val("S_tot") * hartree_to_kj, rel=1e-3
        )


class TestPartialHessian:
    """PHVA path: frequencies() on a cluster with frozen_indices."""

    def test_phva_mode_count_and_stretches(self, opt_water):
        from dataclasses import replace

        # Freeze the oxygen: 2 free H -> 6 PHVA modes, no trans/rot
        # projection, no rotational temperatures.
        frozen = replace(opt_water, frozen_indices=[0])
        freq = frequencies(frozen, CHEAP)
        total = freq.frequencies_cm.size + freq.imaginary_cm.size
        assert total == 6
        assert freq.rotational_temperatures_k is None
        # The OH stretches and the bend survive the freeze.
        assert np.sum(freq.frequencies_cm > 1500) >= 3
        # Residual imaginary modes are pure rotation-about-O noise.
        assert np.all(freq.imaginary_cm < 100)

    def test_phva_energy_matches_full(self, opt_water):
        from dataclasses import replace

        e_full = frequencies(opt_water, CHEAP).electronic_hartree
        frozen = replace(opt_water, frozen_indices=[0])
        e_phva = frequencies(frozen, CHEAP).electronic_hartree
        assert e_phva == pytest.approx(e_full, abs=1e-8)

    def test_phva_imaginary_mode_zero_on_frozen_atoms(self, opt_water):
        from dataclasses import replace

        # Displace one H off the minimum so a real imaginary mode exists,
        # then check the returned mode vector is zero on the frozen atom.
        bent = replace(
            opt_water,
            coords=opt_water.coords
            + np.array([[0, 0, 0], [0.4, -0.3, 0.2], [0, 0, 0]]),
            frozen_indices=[0],
        )
        freq = frequencies(bent, CHEAP)
        if freq.imaginary_mode is not None:
            assert np.allclose(freq.imaginary_mode[0], 0.0)
