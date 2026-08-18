"""Gates for the in-house rate module.

Closed-form pieces are held to 1e-12 (relative); physical sanity checks
bound the rest. Mirrors petra's own eyring gates (petra-core rate.rs).
"""

import math

import numpy as np
import pytest

from quarry.rates import (
    KB,
    NA,
    R_KJ,
    H,
    Thermo,
    arrhenius_fit,
    eckart_kappa,
    eyring_k,
    eyring_k_dh_ds,
    rate_from_thermo,
    thermo_from_frequencies,
    wigner_kappa,
)

T298 = 298.15


class TestEyring:
    def test_zero_barrier_is_kbt_over_h(self):
        assert eyring_k(0.0, T298) == pytest.approx(KB * T298 / H, rel=1e-12)

    def test_hand_computed_barrier(self):
        # ΔG‡ = 50 kJ/mol at 298.15 K: k = (kB T / h) exp(-50/(R T)).
        expected = (KB * T298 / H) * math.exp(-50.0 / (R_KJ * T298))
        assert eyring_k(50.0, T298) == pytest.approx(expected, rel=1e-12)

    def test_dh_ds_form_matches_dg_form(self):
        dh, ds = 60.0, -0.05  # kJ/mol, kJ/(mol K)
        assert eyring_k_dh_ds(dh, ds, T298) == pytest.approx(
            eyring_k(dh - T298 * ds, T298), rel=1e-12
        )

    def test_petra_gate_value(self):
        # petra-core rate.rs computes with the dissertation's truncated
        # gas constant (1.987e-3 kcal/mol/K); quarry uses CODATA R.
        # Document the systematic offset: ~0.17% in k at 10 kcal/mol,
        # 298 K — negligible against DFT barrier error, but real, so the
        # bridge (emit.py) sends barriers, never pre-exponentiated rates.
        ea_kj = 10.0 * 4.184
        k_quarry = 1.0e13 * math.exp(-ea_kj / (R_KJ * T298))
        k_petra = 1.0e13 * math.exp(-10.0 / (1.987e-3 * T298))
        assert k_quarry == pytest.approx(k_petra, rel=3e-3)
        assert k_quarry != pytest.approx(k_petra, rel=1e-4)


class TestTunneling:
    def test_wigner_closed_form(self):
        # 1000i cm^-1 at 298.15 K.
        u = H * 2.99792458e10 * 1000.0 / (KB * T298)
        expected = 1.0 + u * u / 24.0
        assert wigner_kappa(1000.0, T298) == pytest.approx(expected, rel=1e-12)

    def test_wigner_grows_as_t_drops(self):
        assert wigner_kappa(1500.0, 250.0) > wigner_kappa(1500.0, 350.0)

    def test_eckart_exceeds_unity(self):
        assert eckart_kappa(1200.0, 80.0, 90.0, T298) > 1.0

    def test_eckart_approaches_wigner_for_small_frequency(self):
        # For a broad/low-tunneling barrier both corrections -> the same
        # leading-order series.
        k_e = eckart_kappa(150.0, 60.0, 60.0, T298)
        k_w = wigner_kappa(150.0, T298)
        assert k_e == pytest.approx(k_w, rel=2e-2)

    def test_eckart_asymmetry_orientation(self):
        # Exothermic direction (deep product well) tunnels at least as
        # much as the symmetric barrier of the same forward height.
        sym = eckart_kappa(1000.0, 60.0, 60.0, T298)
        exo = eckart_kappa(1000.0, 60.0, 160.0, T298)
        assert exo > 1.0
        assert sym > 1.0

    def test_eckart_high_temperature_limit(self):
        assert eckart_kappa(800.0, 70.0, 70.0, 2000.0) == pytest.approx(1.0, abs=0.1)


class TestThermo:
    def test_monatomic_ideal_gas_sackur_tetrode(self):
        # Argon at 298.15 K, 1 bar: S = 154.846 J/(mol K) (literature).
        th = thermo_from_frequencies(
            0.0,
            [],
            T298,
            molar_mass_kg=39.948e-3,
            rotational_temperatures_k=None,
        )
        assert th.entropy_kj_per_k * 1000.0 == pytest.approx(154.846, abs=0.01)
        # H = 5/2 RT for a monatomic ideal gas.
        assert th.enthalpy == pytest.approx(2.5 * R_KJ * T298, rel=1e-12)

    def test_zpe_closed_form(self):
        # One mode at 1000 cm^-1: ZPE = NA h c nu / 2.
        th = thermo_from_frequencies(
            0.0,
            [1000.0],
            T298,
            molar_mass_kg=18.0e-3,
            rotational_temperatures_k=[10.0, 15.0, 25.0],
        )
        zpe = 0.5 * NA * H * 2.99792458e10 * 1000.0 / 1000.0
        assert th.zpe_kj == pytest.approx(zpe, rel=1e-12)

    def test_gibbs_identity(self):
        th = Thermo(
            electronic_kj=-100.0,
            zpe_kj=50.0,
            thermal_kj=10.0,
            entropy_kj_per_k=0.2,
            temperature=300.0,
        )
        assert th.gibbs == pytest.approx(th.enthalpy - 300.0 * 0.2, rel=1e-12)

    def test_qrrho_lowers_low_frequency_entropy(self):
        common = dict(molar_mass_kg=0.1, rotational_temperatures_k=[1.0, 2.0, 3.0])
        harmonic = thermo_from_frequencies(
            0.0, [20.0, 1500.0], T298, qrrho_cutoff_cm=1e-6, **common
        )
        qrrho = thermo_from_frequencies(
            0.0, [20.0, 1500.0], T298, qrrho_cutoff_cm=100.0, **common
        )
        # A 20 cm^-1 harmonic mode carries absurd entropy; the damped
        # rotor must reduce it.
        assert qrrho.entropy_kj_per_k < harmonic.entropy_kj_per_k

    def test_imaginary_mode_rejected(self):
        with pytest.raises(ValueError):
            thermo_from_frequencies(
                0.0,
                [-500.0, 1000.0],
                T298,
                molar_mass_kg=0.018,
                rotational_temperatures_k=[10.0, 15.0, 25.0],
            )


class TestRateFromThermo:
    @staticmethod
    def _pair(dg_target: float) -> tuple[Thermo, Thermo]:
        reactant = Thermo(0.0, 0.0, 0.0, 0.0, T298)
        ts = Thermo(dg_target, 0.0, 0.0, 0.0, T298)
        return reactant, ts

    def test_matches_direct_eyring(self):
        reactant, ts = self._pair(75.0)
        r = rate_from_thermo(reactant, ts)
        assert r.k == pytest.approx(eyring_k(75.0, T298), rel=1e-12)

    def test_wigner_multiplies(self):
        reactant, ts = self._pair(75.0)
        r = rate_from_thermo(reactant, ts, imag_nu_cm=1100.0, tunneling="wigner")
        assert r.k == pytest.approx(
            wigner_kappa(1100.0, T298) * eyring_k(75.0, T298), rel=1e-12
        )

    def test_temperature_mismatch_rejected(self):
        reactant = Thermo(0.0, 0.0, 0.0, 0.0, 298.15)
        ts = Thermo(50.0, 0.0, 0.0, 0.0, 300.0)
        with pytest.raises(ValueError):
            rate_from_thermo(reactant, ts)


class TestArrheniusFit:
    def test_recovers_exact_parameters(self):
        a, ea = 2.5e12, 85.0
        temps = np.array([280.0, 300.0, 320.0, 340.0])
        ks = a * np.exp(-ea / (R_KJ * temps))
        a_fit, ea_fit = arrhenius_fit(temps, ks)
        assert a_fit == pytest.approx(a, rel=1e-9)
        assert ea_fit == pytest.approx(ea, rel=1e-9)
