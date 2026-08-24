"""Gates for the cluster builders: topology, stoichiometry, sane bonds."""

import numpy as np
import pytest

from quarry.clusters import (
    ACID_MICROSOLVATION_FAMILIES,
    BENCHMARKS,
    Cluster,
    SiteFamily,
    acid_microsolvation_hbond_edges,
    aluminosilicate_dimer,
    disilicate,
    hydronium,
    hydronium_hydrate,
    protonated_bridge_complex,
    silicic_acid,
    silicic_acid_hydrate,
    water,
)


def _min_interatomic(c: Cluster) -> float:
    d = np.linalg.norm(c.coords[:, None, :] - c.coords[None, :, :], axis=-1)
    return float(d[np.triu_indices(len(c.symbols), k=1)].min())


class TestBasics:
    def test_water(self):
        w = water()
        assert w.formula == "H2O"
        assert w.charge == 0

    def test_hydronium_charge(self):
        assert hydronium().charge == 1
        assert hydronium().formula == "H3O"

    def test_silicic_acid(self):
        s = silicic_acid()
        assert s.formula == "H4O4Si"
        assert s.site_family is SiteFamily.SI
        # Four Si-O bonds at the guess length.
        si = s.coords[0]
        o_dists = [
            np.linalg.norm(s.coords[i] - si)
            for i, sym in enumerate(s.symbols)
            if sym == "O"
        ]
        assert len(o_dists) == 4
        assert all(d == pytest.approx(1.63, abs=0.01) for d in o_dists)

    def test_disilicate(self):
        d = disilicate()
        assert d.formula == "H6O7Si2"
        assert d.charge == 0
        assert d.site_family is SiteFamily.SI_O_SI

    def test_aluminosilicate_dimer(self):
        d = aluminosilicate_dimer()
        assert d.formula == "AlH6O7Si"
        assert d.charge == -1
        assert d.site_family is SiteFamily.SI_O_AL2

    @pytest.mark.parametrize(
        ("dimer_factory", "formula", "charge", "family"),
        [
            (disilicate, "H9O8Si2", 1, SiteFamily.SI_O_SI),
            (aluminosilicate_dimer, "AlH9O8Si", 0, SiteFamily.SI_O_AL2),
        ],
    )
    @pytest.mark.parametrize("mode", ["flank", "backside"])
    def test_protonated_bridge_complex_is_pre_equilibrated_water_attack(
        self, dimer_factory, formula, charge, family, mode
    ):
        dimer = dimer_factory()
        complex_ = protonated_bridge_complex(dimer, mode=mode)
        ow_index = len(dimer.symbols)
        bridge_proton = ow_index + 1

        assert complex_.formula == formula
        assert complex_.charge == charge
        assert complex_.site_family is family
        assert complex_.symbols[ow_index : ow_index + 4] == ["O", "H", "H", "H"]
        assert np.linalg.norm(
            complex_.coords[0] - complex_.coords[bridge_proton]
        ) == pytest.approx(0.96, abs=0.01)
        assert (
            np.linalg.norm(complex_.coords[ow_index] - complex_.coords[bridge_proton])
            > 1.25
        )
        assert all(
            np.linalg.norm(complex_.coords[ow_index] - complex_.coords[index])
            == pytest.approx(0.96, abs=0.01)
            for index in (ow_index + 2, ow_index + 3)
        )
        assert _min_interatomic(complex_) > 0.85

    @pytest.mark.parametrize("n_water", [3, 4, 5, 6])
    @pytest.mark.parametrize(
        ("dimer_factory", "charge", "family"),
        [
            (disilicate, 1, SiteFamily.SI_O_SI),
            (aluminosilicate_dimer, 0, SiteFamily.SI_O_AL2),
        ],
    )
    def test_microsolvated_protonated_bridge_has_sane_deterministic_shell(
        self, n_water, dimer_factory, charge, family
    ):
        dimer = dimer_factory()
        complex_ = protonated_bridge_complex(dimer, n_water=n_water)
        repeated = protonated_bridge_complex(dimer, n_water=n_water)
        ow_index = len(dimer.symbols)

        assert len(complex_.symbols) == len(dimer.symbols) + 4 + 3 * (n_water - 1)
        assert complex_.charge == charge
        assert complex_.site_family is family
        assert complex_.symbols[ow_index : ow_index + 4] == ["O", "H", "H", "H"]
        assert np.array_equal(complex_.coords, repeated.coords)
        assert _min_interatomic(complex_) > 0.85
        assert f"{n_water} explicit waters" in complex_.note

    @pytest.mark.parametrize("n_water", [3, 4, 5, 6])
    @pytest.mark.parametrize("family", ACID_MICROSOLVATION_FAMILIES)
    @pytest.mark.parametrize("dimer_factory", [disilicate, aluminosilicate_dimer])
    def test_acid_conformer_families_preserve_atom_contract_without_collisions(
        self, n_water, family, dimer_factory
    ):
        dimer = dimer_factory()
        complex_ = protonated_bridge_complex(
            dimer, n_water=n_water, conformer_family=family
        )
        repeated = protonated_bridge_complex(
            dimer, n_water=n_water, conformer_family=family
        )
        ow_index = len(dimer.symbols)

        assert complex_.symbols == repeated.symbols
        assert np.array_equal(complex_.coords, repeated.coords)
        assert complex_.symbols[ow_index : ow_index + 4] == ["O", "H", "H", "H"]
        assert all(
            complex_.symbols[ow_index + 4 + 3 * i : ow_index + 7 + 3 * i]
            == ["O", "H", "H"]
            for i in range(n_water - 1)
        )
        assert _min_interatomic(complex_) > 0.85
        assert family in complex_.name
        assert family in complex_.note

        edges = acid_microsolvation_hbond_edges(ow_index, n_water, family)
        nodes = {0, ow_index, *(ow_index + 4 + 3 * i for i in range(n_water - 1))}
        adjacency = {node: set() for node in nodes}
        for donor, acceptor in edges:
            adjacency[donor].add(acceptor)
            adjacency[acceptor].add(donor)
            donor_hydrogens = [
                index
                for index, symbol in enumerate(complex_.symbols)
                if symbol == "H"
                and np.linalg.norm(complex_.coords[index] - complex_.coords[donor])
                < 1.25
            ]
            assert donor_hydrogens
            assert (
                min(
                    np.linalg.norm(complex_.coords[index] - complex_.coords[acceptor])
                    for index in donor_hydrogens
                )
                < 2.45
            )
            assert (
                np.linalg.norm(complex_.coords[donor] - complex_.coords[acceptor])
                < 3.40
            )
        visited = {next(iter(nodes))}
        frontier = list(visited)
        while frontier:
            frontier.extend(adjacency[frontier.pop()] - visited)
            visited.update(frontier)
        assert visited == nodes

    @pytest.mark.parametrize("n_water", [3, 4, 5, 6])
    def test_acid_conformer_families_are_distinct_at_each_water_count(self, n_water):
        dimer = disilicate()
        ow_index = len(dimer.symbols)
        extra_oxygen_indices = [ow_index + 4 + 3 * i for i in range(n_water - 1)]
        solvent_shapes = {
            family: protonated_bridge_complex(
                dimer, n_water=n_water, conformer_family=family
            ).coords[extra_oxygen_indices]
            for family in ACID_MICROSOLVATION_FAMILIES
        }

        for index, family in enumerate(ACID_MICROSOLVATION_FAMILIES):
            for other in ACID_MICROSOLVATION_FAMILIES[index + 1 :]:
                assert not np.allclose(solvent_shapes[family], solvent_shapes[other])

    def test_acid_conformer_family_validation(self):
        with pytest.raises(ValueError, match="conformer_family"):
            protonated_bridge_complex(
                disilicate(), n_water=4, conformer_family="random-shell"
            )
        with pytest.raises(ValueError, match="n_water=1"):
            protonated_bridge_complex(
                disilicate(), n_water=1, conformer_family="bridge-donor-chain"
            )

    @pytest.mark.parametrize("n_water", [3, 4, 5, 6])
    def test_hydronium_hydrate_matches_complex_fragment_stoichiometry(self, n_water):
        reagent = hydronium_hydrate(n_water)

        assert len(reagent.symbols) == 4 + 3 * (n_water - 1)
        assert reagent.formula == f"H{2 * n_water + 1}O{n_water}"
        assert reagent.charge == 1
        assert _min_interatomic(reagent) > 0.85

    @pytest.mark.parametrize("n_water", [0, 2, 7])
    def test_microsolvation_rejects_unsupported_shell_sizes(self, n_water):
        with pytest.raises(ValueError, match="n_water"):
            hydronium_hydrate(n_water)
        with pytest.raises(ValueError, match="n_water"):
            protonated_bridge_complex(disilicate(), n_water=n_water)

    def test_no_atom_collisions_in_any_benchmark(self):
        for factory in BENCHMARKS.values():
            c = factory()
            assert _min_interatomic(c) > 0.85, c.name  # > shortest real bond

    def test_hydrate_series(self):
        for n in (0, 1, 3, 6, 12, 24):
            c = silicic_acid_hydrate(n)
            assert len(c.symbols) == 9 + 3 * n
        with pytest.raises(ValueError):
            silicic_acid_hydrate(25)
        with pytest.raises(ValueError):
            silicic_acid_hydrate(-1)

    def test_hydrate_shells_never_collide(self):
        # Both shells populated at the maximum: waters must not overlap
        # each other or the core.
        c = silicic_acid_hydrate(24)
        assert _min_interatomic(c) > 0.85


class TestClusterContract:
    def test_xyz_roundtrip_shape(self):
        s = silicic_acid()
        lines = s.to_xyz().splitlines()
        assert int(lines[0]) == len(s.symbols)
        assert len(lines) == len(s.symbols) + 2

    def test_pyscf_atom_string(self):
        atoms = water().to_pyscf_atom().split("; ")
        assert len(atoms) == 3
        assert atoms[0].startswith("O ")

    def test_frozen_indices_validated(self):
        with pytest.raises(ValueError):
            Cluster("bad", ["O"], np.zeros((1, 3)), frozen_indices=[5])

    def test_symbols_coords_mismatch_rejected(self):
        with pytest.raises(ValueError):
            Cluster("bad", ["O", "H"], np.zeros((1, 3)))
