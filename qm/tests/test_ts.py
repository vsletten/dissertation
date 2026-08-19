"""Gates for the TS machinery, all CPU at HF/STO-3G.

HCN <-> HNC isomerization is the canonical cheap saddle: 3 atoms, a
well-conditioned single imaginary mode, seconds per gradient. It
exercises the whole ladder — ASE adapter, Sella saddle search, one-
imaginary-mode verification, quick-IRC basin check — without touching
anything silicate-sized.
"""

import numpy as np
import pytest

from quarry.clusters import (
    Cluster,
    aluminosilicate_dimer,
    disilicate,
    hydrolysis_complex,
    hydronium,
    merge,
    water,
)
from quarry.pipeline import DftSettings, energy, frequencies, optimize
from quarry.ts import find_ts, make_ase_calculator, quick_irc, verify_ts

CHEAP = DftSettings(xc="hf", basis="sto-3g")


def hcn() -> Cluster:
    return Cluster(
        name="hcn",
        symbols=["C", "N", "H"],
        coords=np.array([[0.0, 0.0, 0.0], [1.16, 0.0, 0.0], [-1.06, 0.0, 0.0]]),
    )


def hcn_ts_guess() -> Cluster:
    # H halfway around the C-N axis — near the known cyclic TS.
    return Cluster(
        name="hcn-ts-guess",
        symbols=["C", "N", "H"],
        coords=np.array([[0.0, 0.0, 0.0], [1.20, 0.0, 0.0], [0.45, 1.15, 0.0]]),
    )


class TestAseAdapter:
    def test_energy_matches_pipeline(self):
        from ase import Atoms

        w = water()
        atoms = Atoms(symbols=w.symbols, positions=w.coords)
        atoms.calc = make_ase_calculator(CHEAP, w.charge, w.spin)
        e_ase = atoms.get_potential_energy() / 27.211386245988
        assert e_ase == pytest.approx(energy(w, CHEAP), rel=1e-9)

    def test_forces_shape_and_direction(self):
        from ase import Atoms

        w = water()
        atoms = Atoms(symbols=w.symbols, positions=w.coords)
        atoms.calc = make_ase_calculator(CHEAP, w.charge, w.spin)
        f = atoms.get_forces()
        assert f.shape == (3, 3)
        assert np.all(np.isfinite(f))

    def test_forces_match_finite_difference(self):
        # Physics gate on signs/units/coordinate handling: F_x(H1) must
        # equal -dE/dx to finite-difference accuracy.
        from ase import Atoms

        w = water()
        atoms = Atoms(symbols=w.symbols, positions=w.coords)
        atoms.calc = make_ase_calculator(CHEAP, w.charge, w.spin)
        f_analytic = atoms.get_forces()[1, 0]  # H1, x-component, eV/A

        h = 1e-3  # Angstrom
        e = []
        for sign in (+1.0, -1.0):
            shifted = Atoms(symbols=w.symbols, positions=w.coords)
            shifted.positions[1, 0] += sign * h
            shifted.calc = make_ase_calculator(CHEAP, w.charge, w.spin)
            e.append(shifted.get_potential_energy())
        f_numeric = -(e[0] - e[1]) / (2.0 * h)
        assert f_analytic == pytest.approx(f_numeric, abs=1e-3)


@pytest.mark.slow
class TestSaddleSearch:
    @pytest.fixture(scope="class")
    def ts(self):
        return find_ts(hcn_ts_guess(), CHEAP, max_steps=200)

    def test_exactly_one_imaginary_mode(self, ts):
        freq = verify_ts(ts, CHEAP)
        # HF/STO-3G HCN<->HNC saddle: one large imaginary frequency.
        assert 500.0 < freq.imaginary_cm[0] < 3500.0
        assert freq.imaginary_mode is not None
        assert freq.imaginary_mode.shape == (3, 3)

    def test_saddle_above_both_minima(self, ts):
        e_ts = energy(ts, CHEAP)
        e_hcn = energy(optimize(hcn(), CHEAP), CHEAP)
        assert e_ts > e_hcn

    def test_quick_irc_reaches_two_distinct_minima(self, ts):
        back, fwd = quick_irc(ts, CHEAP)
        e_back, e_fwd = energy(back, CHEAP), energy(fwd, CHEAP)
        e_ts = energy(ts, CHEAP)
        assert e_back < e_ts and e_fwd < e_ts
        # HCN and HNC differ: the two basins are not the same structure.
        ch_back = np.linalg.norm(back.coords[2] - back.coords[0])
        ch_fwd = np.linalg.norm(fwd.coords[2] - fwd.coords[0])
        assert abs(ch_back - ch_fwd) > 0.3

    def test_neb_peak_approximates_the_saddle(self, ts):
        from quarry.ts import neb_ts_guess

        back, fwd = quick_irc(ts, CHEAP)
        guess = neb_ts_guess(back, fwd, CHEAP, n_images=5, max_steps=80)
        e_guess = energy(guess, CHEAP)
        e_ts = energy(ts, CHEAP)
        # The climbing image should land near the true saddle energy,
        # from above or below within the loose NEB convergence.
        assert abs(e_guess - e_ts) < 0.02  # Hartree (~50 kJ/mol slack)
        assert e_guess > energy(back, CHEAP)


class TestNebConvergence:
    @staticmethod
    def _endpoints():
        reactant = Cluster("r", ["H"], np.array([[0.0, 0.0, 0.0]]))
        product = Cluster("p", ["H"], np.array([[1.0, 0.0, 0.0]]))
        return reactant, product

    @pytest.mark.parametrize(
        ("stage_results", "message"),
        [
            ([False], "pre-relaxation did not converge"),
            ([True, False], "climbing-image"),
        ],
    )
    def test_unconverged_band_is_rejected(self, monkeypatch, stage_results, message):
        import ase.mep
        import ase.optimize

        neb_calls = []
        aligned_endpoint_rms = []
        fire_calls = []
        results = iter(stage_results)

        class FakeNeb:
            def __init__(self, images, **kwargs):
                self.images = images
                self.climb = kwargs["climb"]
                neb_calls.append(kwargs)
                delta = images[-1].positions - images[0].positions
                aligned_endpoint_rms.append(float(np.sqrt(np.mean(delta**2))))

            def interpolate(self, *, method):
                assert method == "idpp"

        class FakeFire:
            def __init__(self, neb, **kwargs):
                self.neb = neb
                fire_calls.append((neb.climb, kwargs))

            def run(self, *, fmax, steps):
                fire_calls[-1] += (fmax, steps)
                return next(results)

        monkeypatch.setattr(ase.mep, "NEB", FakeNeb)
        monkeypatch.setattr(ase.optimize, "FIRE", FakeFire)
        reactant, product = self._endpoints()

        with pytest.raises(RuntimeError, match=message):
            from quarry.ts import neb_ts_guess

            neb_ts_guess(
                reactant,
                product,
                CHEAP,
                n_images=3,
                pre_relax_steps=4,
                max_steps=5,
            )

        assert neb_calls == [
            {
                "climb": False,
                "method": "improvedtangent",
                "remove_rotation_and_translation": True,
            }
        ]
        assert aligned_endpoint_rms == pytest.approx([0.0], abs=1e-12)
        assert fire_calls[0] == (
            False,
            {"dt": 0.05, "dtmax": 0.2, "maxstep": 0.05},
            0.3,
            4,
        )
        if len(stage_results) == 2:
            assert fire_calls[1] == (
                True,
                {"dt": 0.05, "dtmax": 0.2, "maxstep": 0.05},
                0.1,
                5,
            )


class TestConstraints:
    """These run real geomeTRIC constraint plumbing — the inline-text
    regression (FileNotFoundError on the constraint string) died here."""

    def test_constrained_scan_holds_distance_and_raises_energy(self):
        from quarry.ts import constrained_scan, scan_maximum

        w = optimize(water(), CHEAP)
        scan = constrained_scan(w, CHEAP, atom_i=0, atom_j=1, distances_a=[1.10, 1.30])
        for r_target, _e, cl in scan:
            r_actual = np.linalg.norm(cl.coords[0] - cl.coords[1])
            assert r_actual == pytest.approx(r_target, abs=0.01)
        # Stretching an O-H bond off equilibrium must cost energy.
        assert scan[1][1] > scan[0][1]
        assert scan_maximum(scan).name.endswith("r1.30")

    def test_optimize_respects_frozen_atoms(self):
        from dataclasses import replace

        w = replace(water(), frozen_indices=[0, 1])
        opt = optimize(w, CHEAP)
        # Frozen atoms stay put to numerical drift (~1e-5 A observed).
        assert np.allclose(opt.coords[0], w.coords[0], atol=1e-4)
        assert np.allclose(opt.coords[1], w.coords[1], atol=1e-4)
        # The free H must have moved off the deliberately-bad guess.
        assert not np.allclose(opt.coords[2], w.coords[2], atol=1e-3)


class TestScanExtension:
    @staticmethod
    def _fake_point(r, e):
        cl = Cluster(f"p-r{r:.2f}", ["H"], np.zeros((1, 3)))
        return (r, e, cl)

    def test_interior_maximum_detection(self):
        from quarry.ts import has_interior_maximum

        rising = [
            self._fake_point(r, e) for r, e in [(2.8, 0.0), (2.4, 1.0), (2.0, 2.0)]
        ]
        peaked = [
            self._fake_point(r, e) for r, e in [(2.8, 0.0), (2.4, 2.0), (2.0, 1.0)]
        ]
        assert not has_interior_maximum(rising)
        assert has_interior_maximum(peaked)
        assert not has_interior_maximum(peaked[:2])  # too short to say
        # A plateau shared with an endpoint is not a crossed ridge.
        plateau = [
            self._fake_point(r, e) for r, e in [(2.8, 0.0), (2.4, 2.0), (2.0, 2.0)]
        ]
        assert not has_interior_maximum(plateau)

    def test_first_crest_beats_compression_wall(self):
        # The observed proton-transfer profile: crest, bond formation
        # dip, then a monotonic compression wall that out-climbs the
        # crest at the endpoint. The guess must be the first crest.
        from quarry.ts import first_interior_maximum, scan_ts_guess

        profile = [
            (2.26, 0.001),
            (1.96, 0.007),  # the real crest
            (1.81, -0.001),  # proton snapped over
            (1.21, 0.010),
            (1.06, 0.015),  # compression wall, global max at endpoint
        ]
        scan = [self._fake_point(r, e) for r, e in profile]
        assert first_interior_maximum(scan) == 1
        assert scan_ts_guess(scan).name == "p-r1.96"

    def test_scan_ts_guess_requires_interior_crest(self):
        from quarry.ts import scan_ts_guess

        rising = [
            self._fake_point(r, e) for r, e in [(2.8, 0.0), (2.4, 1.0), (2.0, 2.0)]
        ]
        with pytest.raises(ValueError, match="no interior"):
            scan_ts_guess(rising)

    def test_scan_extends_until_peak(self, monkeypatch):
        import quarry.ts as ts_mod

        # Energy profile: rises to a peak at r=1.92 then falls.
        profile = {2.0: 1.0, 1.92: 2.0, 1.84: 1.5}

        def fake_scan(cluster, settings, *, atom_i, atom_j, distances_a, **kwargs):
            return [self._fake_point(r, profile[round(r, 2)]) for r in distances_a]

        monkeypatch.setattr(ts_mod, "constrained_scan", fake_scan)
        scan = ts_mod.scan_to_maximum(
            Cluster("x", ["H"], np.zeros((1, 3))),
            CHEAP,
            atom_i=0,
            atom_j=0,
            distances_a=[2.0, 1.92],
        )
        assert [round(r, 2) for r, _, _ in scan] == [2.0, 1.92, 1.84]
        assert ts_mod.scan_maximum(scan).name == "p-r1.92"

    def test_spike_is_retried_from_both_neighbors_before_crest_detection(
        self, monkeypatch
    ):
        import quarry.ts as ts_mod

        # The r=2.20 point mimics the live al-neutral optimizer hop: it is
        # >40 kJ/mol above the interpolation of its neighbors and would be
        # mistaken for the first crest without a bidirectional repair.
        profile = {
            2.60: 0.000,
            2.40: 0.010,
            2.20: 0.060,
            2.10: 0.023,
            2.00: 0.015,
        }
        retry_energies = {"p-r2.40": 0.020, "p-r2.10": 0.014}
        retry_seeds = []

        def fake_scan(cluster, settings, *, atom_i, atom_j, distances_a, **kwargs):
            if len(distances_a) > 1:
                return [self._fake_point(r, profile[round(r, 2)]) for r in distances_a]
            r = distances_a[0]
            retry_seeds.append(cluster.name)
            return [self._fake_point(r, retry_energies[cluster.name])]

        monkeypatch.setattr(ts_mod, "constrained_scan", fake_scan)
        scan = ts_mod.scan_to_maximum(
            Cluster("x", ["H"], np.zeros((1, 3))),
            CHEAP,
            atom_i=0,
            atom_j=0,
            distances_a=[2.60, 2.40, 2.20, 2.10, 2.00],
        )

        assert retry_seeds == ["p-r2.40", "p-r2.10"]
        assert scan[2][1] == pytest.approx(0.014)
        assert scan[2][2].name == "p-r2.20"
        assert ts_mod.first_interior_maximum(scan) == 3
        assert ts_mod.scan_ts_guess(scan).name == "p-r2.10"

    def test_scan_raises_at_floor_without_peak(self, monkeypatch):
        import quarry.ts as ts_mod

        def fake_scan(cluster, settings, *, atom_i, atom_j, distances_a, **kwargs):
            return [self._fake_point(r, 3.0 - r) for r in distances_a]

        monkeypatch.setattr(ts_mod, "constrained_scan", fake_scan)
        with pytest.raises(ts_mod.ScanNoMaximumError, match="interior") as exc_info:
            ts_mod.scan_to_maximum(
                Cluster("x", ["H"], np.zeros((1, 3))),
                CHEAP,
                atom_i=0,
                atom_j=0,
                distances_a=[1.8, 1.75],
            )
        # The exception carries the scan so callers can pivot to a
        # second coordinate from the best structure already computed.
        assert len(exc_info.value.scan) >= 2
        assert exc_info.value.scan[0][0] == 1.8

    def test_fixed_distances_held_during_scan(self):
        from quarry.ts import constrained_scan

        w = optimize(water(), CHEAP)
        scan = constrained_scan(
            w,
            CHEAP,
            atom_i=0,
            atom_j=1,
            distances_a=[1.15],
            fixed_distances=[(0, 2, 1.05)],
        )
        _, _, cl = scan[0]
        assert np.linalg.norm(cl.coords[0] - cl.coords[1]) == pytest.approx(
            1.15, abs=0.01
        )
        assert np.linalg.norm(cl.coords[0] - cl.coords[2]) == pytest.approx(
            1.05, abs=0.01
        )


class TestVerifyTs:
    def test_minimum_rejected_as_ts(self):
        w = optimize(water(), CHEAP)
        with pytest.raises(RuntimeError, match="expected exactly 1"):
            verify_ts(w, CHEAP)

    def test_minimum_has_no_imaginary_mode_vector(self):
        w = optimize(water(), CHEAP)
        freq = frequencies(w, CHEAP)
        assert freq.imaginary_mode is None


class TestComplexBuilders:
    def test_merge_offsets_frozen_indices_and_sums_charge(self):
        a = disilicate()
        a.frozen_indices = [1]
        b = hydronium()
        b.frozen_indices = [0]
        m = merge(a, b)
        assert m.charge == 1
        assert m.frozen_indices == [1, len(a.symbols)]
        assert len(m.symbols) == len(a.symbols) + len(b.symbols)

    def test_backside_attack_complex(self):
        c = hydrolysis_complex(disilicate(), water(), mode="backside")
        assert c.charge == 0
        assert len(c.symbols) == 15 + 3
        # Attacker O sits ~approach distance from the Si under attack,
        # farther from the bridging O than the Si is (backside).
        si, o_attack = c.coords[1], c.coords[15]
        assert np.linalg.norm(o_attack - si) == pytest.approx(3.2, abs=0.05)
        assert np.linalg.norm(o_attack - c.coords[0]) > np.linalg.norm(si - c.coords[0])

    def test_flank_attack_complex_geometry(self):
        c = hydrolysis_complex(disilicate(), water())  # flank is default
        obr, si, ow = c.coords[0], c.coords[1], c.coords[15]
        assert np.linalg.norm(ow - si) == pytest.approx(3.2, abs=0.05)
        # Ow-Si-Obr angle near the 4-center arrangement (80 deg built in).
        v1, v2 = ow - si, obr - si
        cos = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        assert 60.0 < np.degrees(np.arccos(cos)) < 100.0
        # The aimed water H must actually be able to reach the bridging O
        # (backside had it at 3.9 A — the measured failure).
        h_to_obr = min(
            np.linalg.norm(c.coords[16] - obr), np.linalg.norm(c.coords[17] - obr)
        )
        assert h_to_obr < 2.6

    def test_unknown_mode_rejected(self):
        with pytest.raises(ValueError, match="unknown attack mode"):
            hydrolysis_complex(disilicate(), water(), mode="sideways")

    def test_acid_attack_complex_charge(self):
        c = hydrolysis_complex(aluminosilicate_dimer(), hydronium())
        assert c.charge == 0  # -1 dimer + +1 hydronium
        assert c.name == "aluminosilicate-dimer+hydronium"
