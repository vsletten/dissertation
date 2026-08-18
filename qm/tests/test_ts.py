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

    def test_neutral_attack_complex(self):
        c = hydrolysis_complex(disilicate(), water())
        assert c.charge == 0
        assert len(c.symbols) == 15 + 3
        # Attacker O sits ~approach distance from the Si under attack,
        # farther from the bridging O than the Si is (backside).
        si, o_attack = c.coords[1], c.coords[15]
        assert np.linalg.norm(o_attack - si) == pytest.approx(3.2, abs=0.05)
        assert np.linalg.norm(o_attack - c.coords[0]) > np.linalg.norm(si - c.coords[0])

    def test_acid_attack_complex_charge(self):
        c = hydrolysis_complex(aluminosilicate_dimer(), hydronium())
        assert c.charge == 0  # -1 dimer + +1 hydronium
        assert c.name == "aluminosilicate-dimer+hydronium"
