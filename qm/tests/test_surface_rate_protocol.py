"""Fast gates for the D2b surface-rate-protocol campaign driver."""

import math

import numpy as np
import pytest

from quarry.rates import (
    surface_thermo_from_frequencies,
    thermo_from_frequencies,
)
from scripts import surface_rate_protocol as surf


def test_surface_thermo_drops_translation_rotation_and_pv():
    freqs = [120.0, 800.0, 1600.0, 3200.0]
    gas = thermo_from_frequencies(
        -1000.0,
        freqs,
        100.0,
        molar_mass_kg=0.031,
        rotational_temperatures_k=[4.0, 1.0, 0.9],
    )
    surface = surface_thermo_from_frequencies(-1000.0, freqs, 100.0)

    assert surface.zpe_kj == pytest.approx(gas.zpe_kj, rel=1e-12)
    # Vibration-only: strictly less thermal energy and entropy than the gas.
    assert surface.thermal_kj < gas.thermal_kj
    assert surface.entropy_kj_per_k < gas.entropy_kj_per_k
    assert surface.temperature == 100.0


def test_surface_thermo_is_mass_independent():
    freqs = [500.0, 1500.0]
    a = surface_thermo_from_frequencies(0.0, freqs, 80.0)
    b = surface_thermo_from_frequencies(0.0, freqs, 80.0)
    assert a.gibbs == pytest.approx(b.gibbs, rel=1e-14)


def test_surface_thermo_atom_has_zero_vibrational_terms():
    atom = surface_thermo_from_frequencies(-10.0, [], 50.0)
    assert atom.zpe_kj == 0.0
    assert atom.thermal_kj == 0.0
    assert atom.entropy_kj_per_k == 0.0
    assert atom.gibbs == pytest.approx(-10.0)


def test_reaction_set_covers_families_sites_and_doublets():
    targets = surf.reactions(gpu=True, basis="def2-svp")

    assert set(targets) == {
        "h-co",
        "h-co-1w-oside",
        "h-co-1w-cside",
        "h-co-2w",
        "h-h2co-ch3o",
        "h-h2co-ch3o-1w",
        "h-h2co-ch3o-2w",
        "h-h2co-h2-hco",
        "h-h2co-h2-hco-1w",
        "h-h2co-h2-hco-2w",
    }
    for reaction in targets.values():
        assert reaction.cluster.spin == 1
        assert reaction.method.use_gpu
        # Substrate order is invariant: scanned atoms must be substrate atoms.
        n_substrate = len(reaction.cluster.symbols) - 3 * reaction.n_water
        assert reaction.scan_i < n_substrate
        assert reaction.scan_j < n_substrate

    families = {reaction.family for reaction in targets.values()}
    assert families == {"h-co", "h-h2co-ch3o", "h-h2co-h2-hco"}
    # Site sampling: every family has a gas reference and >= 2 wet sites.
    for family in families:
        sites = {r.site for r in targets.values() if r.family == family}
        assert "gas" in sites
        assert len(sites - {"gas"}) >= 2


def test_water_blocks_are_appended_after_substrate():
    targets = surf.reactions(gpu=False, basis="def2-svp")
    wet = targets["h-h2co-ch3o-2w"]
    dry = targets["h-h2co-ch3o"]

    assert wet.cluster.symbols[: len(dry.cluster.symbols)] == dry.cluster.symbols
    np.testing.assert_allclose(
        wet.cluster.coords[: len(dry.cluster.symbols)], dry.cluster.coords
    )
    assert wet.cluster.symbols[len(dry.cluster.symbols) :] == ["O", "H", "H"] * 2


def test_cluster_seeds_have_no_atom_collisions():
    targets = surf.reactions(gpu=False, basis="def2-svp")
    for reaction in targets.values():
        coords = reaction.cluster.coords
        deltas = coords[:, None, :] - coords[None, :, :]
        distances = np.linalg.norm(deltas, axis=2)
        n = len(coords)
        off_diagonal = distances[~np.eye(n, dtype=bool)]
        assert off_diagonal.min() > 0.85, reaction.key


def test_literature_fit_matches_modified_arrhenius_form():
    fit = surf.LiteratureFit(3146e10, 1.0, 830.0, 119.6, "fixture", 59.0)
    expected = (
        3146e10
        * (75.0 / 300.0)
        * np.exp(-830.0 * (75.0 + 119.6) / (75.0**2 + 119.6**2))
    )
    assert fit.rate(75.0) == pytest.approx(expected, rel=1e-14)


def test_fuchs_tables_are_the_published_rows():
    assert [row.temperature_k for row in surf.FUCHS_H_CO] == [12.0, 13.5, 15.0, 16.5]
    assert surf.FUCHS_H_CO[0].barrier_k == 390.0
    assert surf.FUCHS_H_H2CO[-1].barrier_k == 500.0
    for row in (*surf.FUCHS_H_CO, *surf.FUCHS_H_H2CO):
        assert row.rate_s > 0.0
        assert row.barrier_err_k > 0.0


def test_fuchs_effective_barrier_roundtrip():
    # k = nu exp(-E/T)  <=>  E = -T ln(k/nu)
    barrier_k = 450.0
    t = 15.0
    k = surf.FUCHS_PREFACTOR_S * math.exp(-barrier_k / t)
    assert surf.fuchs_effective_barrier_k(k, t) == pytest.approx(barrier_k, rel=1e-12)


def test_crossover_temperature_matches_hbar_omega_over_2pi_kb():
    # 831 cm^-1 (Song & Kaestner's CH3O channel): T_c = 1.4388*831/(2 pi).
    assert surf.crossover_temperature_k(831.0) == pytest.approx(190.3, abs=0.2)
    assert surf.crossover_temperature_k(-831.0) == surf.crossover_temperature_k(831.0)


def test_rate_table_cc_correction_shifts_barrier_consistently():
    class FakeFreq:
        def __init__(self, elec, freqs):
            self.electronic_hartree = elec
            self.frequencies_cm = np.array(freqs)
            self.imaginary_cm = np.array([])
            self.molar_mass_kg = 0.031
            self.rotational_temperatures_k = [4.0, 1.0, 0.9]
            self.linear = False

    reactant = FakeFreq(-100.0, [100.0, 900.0, 2000.0])
    ts = FakeFreq(-99.995, [850.0, 1900.0])
    rows = surf.rate_table(
        reactant,
        ts,
        imag_cm=900.0,
        barrier_zpe_dft_kj=15.0,
        reverse_zpe_dft_kj=120.0,
        cc_delta_kj=3.0,
        literature_fit=None,
        temperatures=(100.0, 300.0),
        extrapolated=False,
    )
    for row in rows:
        assert row["k_headline_s"] == row["k_surface_cc_s"]
        assert row["k_surface_cc_s"] < row["k_surface_dft_s"]  # higher barrier
        assert row["eckart_extrapolated"] is False
        assert row["fuchs_effective_barrier_k"] > 0.0


def test_rate_table_without_cc_falls_back_to_dft_headline():
    class FakeFreq:
        def __init__(self, elec, freqs):
            self.electronic_hartree = elec
            self.frequencies_cm = np.array(freqs)
            self.imaginary_cm = np.array([])
            self.molar_mass_kg = 0.031
            self.rotational_temperatures_k = [4.0, 1.0, 0.9]
            self.linear = False

    rows = surf.rate_table(
        FakeFreq(-100.0, [100.0, 900.0, 2000.0]),
        FakeFreq(-99.995, [850.0, 1900.0]),
        imag_cm=900.0,
        barrier_zpe_dft_kj=15.0,
        reverse_zpe_dft_kj=120.0,
        cc_delta_kj=None,
        literature_fit=surf.LiteratureFit(1e12, 1.0, 800.0, 100.0, "fixture", 59.0),
        temperatures=(50.0,),
        extrapolated=True,
    )
    (row,) = rows
    assert row["k_headline_s"] == row["k_surface_dft_s"]
    assert "k_surface_cc_s" not in row
    assert row["fit_extrapolated"] is True
    assert row["eckart_extrapolated"] is True


def test_geometry_hash_ignores_cluster_name():
    cluster = surf.reactions(gpu=False, basis="sto-3g")["h-co"].cluster
    renamed = type(cluster)(
        name=f"{cluster.name}-ts-back",
        symbols=list(cluster.symbols),
        coords=cluster.coords.copy(),
        spin=cluster.spin,
    )
    assert surf.geometry_hash(cluster) == surf.geometry_hash(renamed)


def test_geometry_hash_changes_with_coordinates():
    cluster = surf.reactions(gpu=False, basis="sto-3g")["h-co"].cluster
    moved = cluster.coords.copy()
    moved[-1, 0] += 0.1
    assert surf.geometry_hash(cluster) != surf.geometry_hash(
        type(cluster)(
            name=cluster.name,
            symbols=list(cluster.symbols),
            coords=moved,
            spin=cluster.spin,
        )
    )


def _fake_scan(energies):
    payload = surf.reactions(gpu=False, basis="sto-3g")["h-co"].cluster
    radii = np.linspace(2.8, 1.4, len(energies))
    return [(float(r), float(e), payload) for r, e in zip(radii, energies, strict=True)]


def test_interior_global_maximum_accepts_shallow_interior_crest():
    # The live h-co failure: crest clears one neighbor by only 2e-5 Ha but
    # clears both endpoints by miles.
    scan = _fake_scan([-113.7714, -113.7688, -113.76779, -113.76781, -113.789])
    guess = surf.interior_global_maximum(scan)
    assert guess is scan[2][2]


def test_interior_global_maximum_rejects_endpoint_maximum():
    with pytest.raises(ValueError, match="endpoint"):
        surf.interior_global_maximum(_fake_scan([-1.0, -1.1, -1.2, -1.3]))


def test_interior_global_maximum_rejects_noise_level_rise():
    flat = [-1.0, -1.0 + 5e-5, -1.0 + 2e-5]
    with pytest.raises(ValueError, match="noise"):
        surf.interior_global_maximum(_fake_scan(flat))


def test_saddle_gate_rejects_soft_and_multiple_imaginary_modes():
    class FakeFreq:
        def __init__(self, imaginary):
            self.imaginary_cm = np.array(imaginary)

    with pytest.raises(RuntimeError, match="trivial rearrangement"):
        surf.require_chemical_saddle(FakeFreq([120.0]), name="soft")
    with pytest.raises(RuntimeError, match="expected exactly 1"):
        surf.require_chemical_saddle(FakeFreq([900.0, 300.0]), name="double")
    assert surf.require_chemical_saddle(FakeFreq([900.0, 20.0]), name="ok") == 900.0


def test_minimum_gate_tolerates_only_numerical_noise():
    class FakeFreq:
        def __init__(self, imaginary):
            self.imaginary_cm = np.array(imaginary)

    surf.require_minimum(FakeFreq([12.0, 25.0]), name="soft-librations")
    with pytest.raises(RuntimeError, match="not a minimum"):
        surf.require_minimum(FakeFreq([80.0]), name="bad")
