"""Gates for the PySCF-backed pipeline rungs.

Everything here is CPU-only and tiny (H2O at HF/STO-3G) so the default
gate stays green in any session. The thermo cross-check gates quarry's
own formulas against pyscf's audited harmonic implementation.
"""

from pathlib import Path

import numpy as np
import pytest

from quarry import pipeline
from quarry.clusters import water
from quarry.pipeline import (
    DftSettings,
    FrequencyResult,
    build_mol,
    energy,
    frequencies,
    frequencies_finite_difference,
    frequency_geometry_fingerprint,
    frequency_settings_fingerprint,
    gradient,
    optimize,
)
from quarry.rates import thermo_from_frequencies

CHEAP = DftSettings(xc="hf", basis="sto-3g")
R2SCAN3C = DftSettings(
    xc="r2scan",
    basis="def2-mtzvpp",
    composite="r2scan3c",
    density_fit=True,
)


def test_r2scan3c_settings_refuse_partial_or_double_counted_variants():
    assert pipeline._normalized_composite(R2SCAN3C) == "r2scan3c"
    with pytest.raises(ValueError, match="requires xc='r2scan'"):
        pipeline._normalized_composite(
            DftSettings(xc="b3lyp", basis="def2-mtzvpp", composite="r2scan3c")
        )
    with pytest.raises(ValueError, match="supplies its own D4"):
        pipeline._normalized_composite(
            DftSettings(
                xc="r2scan",
                basis="def2-mtzvpp",
                composite="r2scan3c",
                dispersion="d4",
            )
        )


def test_r2scan3c_d4_set_param_uses_official_keyword_tuple(monkeypatch):
    """Pin the real pyscf-dispersion set_param contract, not (s6, s8, a1, a2).

    Official r2SCAN-3c D4: s6=1.0, s8=0.0, a1=0.42, a2=5.65, s9=2.0
    (Grimme et al., J. Chem. Phys. 154, 064103 (2021)).
    """
    import sys
    import types

    init_kwargs: dict[str, object] = {}
    set_args: tuple[object, ...] = ()
    set_kwargs: dict[str, object] = {}

    class FakeD4:
        def __init__(self, *args, **kwargs):
            init_kwargs.update(kwargs)

        def set_param(self, *args, **kwargs):
            nonlocal set_args, set_kwargs
            set_args = args
            set_kwargs = kwargs

        def get_dispersion(self, grad=False):
            return {"energy": 0.0}

    class FakeGcp:
        def __init__(self, *args, **kwargs):
            pass

        def get_counterpoise(self, grad=False):
            return {"energy": 0.0}

    fake_dftd4 = types.SimpleNamespace(DFTD4Dispersion=FakeD4)
    fake_gcp = types.SimpleNamespace(GCP=FakeGcp)
    fake_disp = types.SimpleNamespace(dftd4=fake_dftd4, gcp=fake_gcp)
    monkeypatch.setitem(sys.modules, "pyscf.dispersion", fake_disp)
    monkeypatch.setitem(sys.modules, "pyscf.dispersion.dftd4", fake_dftd4)
    monkeypatch.setitem(sys.modules, "pyscf.dispersion.gcp", fake_gcp)

    pipeline._r2scan3c_correction(object(), gradient=False, use_gpu=False)

    assert set_args == ()
    assert set_kwargs == {
        "s6": 1.0,
        "s8": 0.0,
        "a1": 0.42,
        "a2": 5.65,
        "s9": 2.0,
    }
    assert init_kwargs["xc"] == "r2scan-3c"


def test_r2scan3c_energy_includes_exact_d4_and_gcp_terms():
    pytest.importorskip("pyscf.dispersion.gcp")
    composite = energy(water(), R2SCAN3C)
    plain_settings = DftSettings(xc="r2scan", basis="def2-mtzvpp", density_fit=True)
    plain = energy(water(), plain_settings)
    correction = pipeline._r2scan3c_correction(
        build_mol(water(), R2SCAN3C), gradient=False, use_gpu=False
    )["energy"]
    assert composite - plain == pytest.approx(float(correction), abs=1e-9)


def test_r2scan3c_patches_direct_geometric_gradient_factory(monkeypatch):
    expected = np.full((3, 3), 0.125)

    class Gradient:
        def __init__(self, base):
            self.base = base

    class MeanField:
        mol = object()
        scf_summary = {}

        def nuc_grad_method(self):
            return Gradient(self)

    monkeypatch.setattr(
        pipeline,
        "_r2scan3c_correction",
        lambda mol, *, gradient, use_gpu: {
            "energy": -0.01,
            "gradient": expected,
        },
    )
    mean_field = pipeline._attach_composite_energy(MeanField(), R2SCAN3C)
    derivative = mean_field.nuc_grad_method()
    assert np.array_equal(derivative.get_dispersion(), expected)


def test_frequency_settings_fingerprint_ignores_execution_backend():
    cpu = DftSettings(xc="b3lyp", basis="def2-svp", density_fit=True)
    gpu = DftSettings(xc="b3lyp", basis="def2-svp", density_fit=True, use_gpu=True)

    assert frequency_settings_fingerprint(cpu) == frequency_settings_fingerprint(gpu)


def test_frequency_geometry_fingerprint_is_exact_identity():
    """The geometry fingerprint is an exact cache key, not a descriptor.

    It hashes the raw coordinate bytes (plus symbols/charge/spin) so a
    reusable Hessian is only ever reused for the *identical* geometry.
    Any coordinate change — including a rigid translation or rotation —
    must therefore change the fingerprint; that is the point, not a bug.
    """
    from dataclasses import replace

    base = water()
    base_fp = frequency_geometry_fingerprint(base)

    # Deterministic: an identical cluster yields the identical fingerprint.
    assert frequency_geometry_fingerprint(water()) == base_fp

    # A rigid translation changes the raw coordinates, so the exact-identity
    # key changes (the cached Hessian is not reused for a shifted geometry).
    translated = replace(base, coords=base.coords + np.array([1.0, 2.0, 3.0]))
    assert frequency_geometry_fingerprint(translated) != base_fp

    # Electronic-state changes also change the key.
    charged = replace(base, charge=1)
    assert frequency_geometry_fingerprint(charged) != base_fp


def test_frequency_geometry_fingerprint_includes_frozen_indices():
    """The frozen-atom set is part of the geometry identity.

    PHVA (partial Hessian) and TS verification branch on ``frozen_indices``,
    so a Hessian computed with a different frozen set must never be reused
    for the same coordinates.
    """
    from dataclasses import replace

    base = water()
    frozen = replace(base, frozen_indices=[0])

    assert frequency_geometry_fingerprint(frozen) != frequency_geometry_fingerprint(
        base
    )
    # Order-independent: the frozen set is sorted before hashing.
    assert frequency_geometry_fingerprint(
        replace(base, frozen_indices=[0, 1])
    ) == frequency_geometry_fingerprint(replace(base, frozen_indices=[1, 0]))


def test_gpu_hessian_assertion_retries_only_hessian_on_cpu(monkeypatch):
    calls = []

    class FakeHessian:
        def __init__(self, gpu):
            self.gpu = gpu

        def kernel(self):
            if self.gpu:
                raise AssertionError("gpu4pyscf non-contiguous UKS Hessian")
            return np.eye(3)

    class FakeScf:
        converged = True

        def __init__(self, gpu):
            self.gpu = gpu

        def kernel(self):
            calls.append(self.gpu)
            return -1.0

        def Hessian(self):
            return FakeHessian(self.gpu)

    monkeypatch.setattr(pipeline, "build_mol", lambda cluster, settings: object())
    monkeypatch.setattr(
        pipeline,
        "_make_scf",
        lambda mol, settings: FakeScf(settings.use_gpu),
    )
    gpu = DftSettings(xc="b3lyp", basis="sto-3g", use_gpu=True)

    with pytest.warns(RuntimeWarning, match="retrying this Hessian on CPU"):
        mf, energy_value, hessian = pipeline._scf_hessian(water(), gpu)

    assert calls == [True, False]
    assert mf.gpu is False
    assert energy_value == -1.0
    assert np.array_equal(hessian, np.eye(3))


def test_gpu_hessian_unrelated_assertion_is_not_masked(monkeypatch):
    """A non-contiguity AssertionError must propagate, not trigger CPU fallback."""

    class FakeHessian:
        def kernel(self):
            raise AssertionError("some unrelated programming bug")

    class FakeScf:
        converged = True

        def kernel(self):
            return -1.0

        def Hessian(self):
            return FakeHessian()

    monkeypatch.setattr(pipeline, "build_mol", lambda cluster, settings: object())
    monkeypatch.setattr(
        pipeline,
        "_make_scf",
        lambda mol, settings: FakeScf(),
    )
    gpu = DftSettings(xc="b3lyp", basis="sto-3g", use_gpu=True)

    with pytest.raises(AssertionError, match="unrelated programming bug"):
        pipeline._scf_hessian(water(), gpu)


@pytest.fixture(scope="module")
def opt_water():
    return optimize(water(), CHEAP)


@pytest.fixture(scope="module")
def freq_water(opt_water) -> FrequencyResult:
    return frequencies(opt_water, CHEAP)


def test_finite_difference_hessian_matches_analytic_water(opt_water, freq_water):
    numerical = frequencies_finite_difference(opt_water, CHEAP, step_bohr=1e-3)
    assert numerical.n_imaginary == freq_water.n_imaginary == 0
    assert numerical.frequencies_cm == pytest.approx(
        freq_water.frequencies_cm,
        abs=0.2,
    )
    assert numerical.electronic_hartree == pytest.approx(
        freq_water.electronic_hartree,
        abs=1e-10,
    )


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
    def test_exhausted_step_bound_is_rejected(self, monkeypatch):
        from pyscf.geomopt import geometric_solver

        monkeypatch.setattr(pipeline, "build_mol", lambda cluster, settings: object())
        monkeypatch.setattr(pipeline, "_make_scf", lambda mol, settings: object())

        class FakeMol:
            natm = len(water().symbols)

            def atom_symbol(self, index):
                return water().symbols[index]

            def atom_coords(self):
                return water().coords / pipeline.BOHR_TO_ANGSTROM

        monkeypatch.setattr(
            geometric_solver,
            "kernel",
            lambda method, **kwargs: (False, FakeMol()),
        )

        with pytest.raises(RuntimeError, match="did not converge within 7 steps"):
            pipeline.optimize(water(), CHEAP, max_steps=7)

    def test_bounded_optimization_combines_distance_and_frozen_constraints(
        self, monkeypatch
    ):
        from dataclasses import replace

        from pyscf.geomopt import geometric_solver

        cluster = replace(water(), frozen_indices=[2])
        monkeypatch.setattr(pipeline, "build_mol", lambda cluster, settings: object())
        monkeypatch.setattr(pipeline, "_make_scf", lambda mol, settings: object())
        captured = {}

        class FakeMol:
            natm = len(cluster.symbols)

            def atom_symbol(self, index):
                return cluster.symbols[index]

            def atom_coords(self):
                return cluster.coords / pipeline.BOHR_TO_ANGSTROM

        def fake_kernel(method, **kwargs):
            captured["constraint"] = Path(kwargs["constraints"]).read_text()
            return False, FakeMol()

        monkeypatch.setattr(geometric_solver, "kernel", fake_kernel)

        result = pipeline.optimize_bounded(
            cluster,
            CHEAP,
            max_steps=9,
            fixed_distances=[(0, 1, 0.9572)],
        )

        assert result.converged is False
        assert captured["constraint"] == (
            "$set\ndistance 1 2 0.95720000\n$freeze\nxyz 3\n"
        )

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

    def test_result_carries_geometry_and_settings_fingerprints(
        self, opt_water, freq_water
    ):
        # frequencies() must stamp the result with the same fingerprints the
        # downstream reuse logic (quick_irc, PHVA) recomputes to detect stale
        # or mismatched Hessians.
        assert freq_water.geometry_fingerprint == frequency_geometry_fingerprint(
            opt_water
        )
        assert freq_water.settings_fingerprint == frequency_settings_fingerprint(CHEAP)


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
