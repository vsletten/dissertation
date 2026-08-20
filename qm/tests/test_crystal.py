"""Gates for the crystallographic cluster builder (HANDOFF.md §2a).

All CPU-only, no DFT: taxonomy match against the deck cell, termination
and charge bookkeeping, frozen shell peripheral, no atom collisions,
cross-check of the deck's crystallography against legacy data.cell.
"""

from pathlib import Path

import numpy as np
import pytest

from quarry.clusters import SiteFamily, water
from quarry.crystal import (
    KIND_TO_FAMILY,
    DeckCell,
    attack_complex,
    from_deck_cell,
    min_interatomic_distance,
)

REPO = Path(__file__).resolve().parents[2]
DECK = REPO / "petra" / "examples" / "kaolinite.toml"
LEGACY_CELL = REPO / "legacy" / "cpp-model" / "data.cell"

LEGACY_TYPE_TO_KIND = {1: "Al", 2: "Si", 3: "Oss", 4: "Osa", 5: "Oaa"}


@pytest.fixture(scope="module")
def cell() -> DeckCell:
    return DeckCell.from_deck(DECK)


class TestDeckCell:
    def test_site_census(self, cell):
        counts: dict[str, int] = {}
        for s in cell.sites:
            counts[s.kind] = counts.get(s.kind, 0) + 1
        assert counts == {"Al": 4, "Si": 4, "Oss": 6, "Osa": 4, "Oaa": 8}
        assert len(cell.bonds) == 40

    def test_taxonomy_covers_all_kinds(self, cell):
        assert {s.kind for s in cell.sites} <= set(KIND_TO_FAMILY)
        assert KIND_TO_FAMILY["Oss"] is SiteFamily.SI_O_SI
        assert KIND_TO_FAMILY["Osa"] is SiteFamily.SI_O_AL2
        assert KIND_TO_FAMILY["Oaa"] is SiteFamily.AL_OH_AL

    def test_graph_degrees_are_full_coordination(self, cell):
        expected = {"Al": 6, "Si": 4, "Oss": 2, "Osa": 3, "Oaa": 2}
        for i, s in enumerate(cell.sites):
            assert cell.degree(i) == expected[s.kind], f"site {i} ({s.kind})"

    def test_bond_lengths_sane(self, cell):
        # Every bonded M-O pair sits at a bondable distance. Bounds are
        # loose because the dissertation-era cell is crude (1.38-2.73 A
        # measured — identical in legacy data.cell, so the deck import is
        # faithful); the pipeline re-optimizes the free region, and the
        # frozen shell carries this as a documented systematic.
        for b in cell.bonds:
            ri = cell.cart((b.i, (0, 0, 0)))
            rj = cell.cart((b.j, b.dcell))
            d = float(np.linalg.norm(ri - rj))
            assert 1.3 < d < 3.0, f"bond {b.i}-{b.j} at {d:.2f} A"


def _parse_legacy_cell(path: Path):
    """(kinds, adjacency) from legacy data.cell; adjacency is per-site
    sets of (neighbor, dcell). -1 entries are padding."""
    lines = path.read_text().splitlines()
    n_sites = int(lines[2].split()[0])
    kinds: dict[int, str] = {}
    adj: dict[int, set] = {}
    current = None
    for line in lines[3:]:
        parts = line.split("#")[0].split()
        if not parts:
            continue
        if len(parts) >= 9 or (len(parts) >= 5 and "." in parts[2]):
            current = int(parts[0])
            kinds[current] = LEGACY_TYPE_TO_KIND[int(parts[1]) // 100]
            adj[current] = set()
            nbr = parts[5:9]
        else:
            nbr = parts[:4]
        if len(nbr) == 4 and current is not None:
            j = int(nbr[0])
            if j >= 0:
                adj[current].add((j, tuple(int(x) for x in nbr[1:4])))
    assert len(kinds) == n_sites
    return kinds, adj


class TestLegacyCrossCheck:
    def test_kinds_match_data_cell(self, cell):
        kinds, _ = _parse_legacy_cell(LEGACY_CELL)
        assert [cell.sites[i].kind for i in range(len(cell.sites))] == [
            kinds[i] for i in range(len(kinds))
        ]

    def test_bond_graph_matches_data_cell(self, cell):
        _, legacy_adj = _parse_legacy_cell(LEGACY_CELL)
        deck_adj = {
            i: {(j, d) for j, d in cell._adjacency[i]} for i in range(len(cell.sites))
        }
        assert deck_adj == legacy_adj


class TestBuilder:
    def test_oss_reference_dimer(self):
        cc = from_deck_cell(DECK, "Oss", metal_shells=1)
        c = cc.cluster
        # The crystallographic analog of the X&L disilicate benchmark.
        assert c.formula == "H6O7Si2"
        assert c.charge == 0
        assert cc.bridge_index == 0
        assert c.symbols[cc.attacked_index] == "Si"
        assert cc.n_intact == 4  # all four bridges of the attacked Si intact
        assert not c.frozen_indices  # metal_shells < 2: free cluster
        assert min_interatomic_distance(c.coords) > 0.75
        # Bridge O actually bridges: both Si within bonding distance.
        br = c.coords[cc.bridge_index]
        si_dists = [
            np.linalg.norm(c.coords[i] - br)
            for i, s in enumerate(c.symbols)
            if s == "Si"
        ]
        assert len(si_dists) == 2
        assert all(1.4 < d < 2.0 for d in si_dists)

    def test_oaa_structural_hydroxyl(self):
        cc = from_deck_cell(DECK, "Oaa", metal_shells=1)
        c = cc.cluster
        assert c.charge == 0
        br = c.coords[cc.bridge_index]
        h_dists = [
            np.linalg.norm(c.coords[i] - br)
            for i, s in enumerate(c.symbols)
            if s == "H"
        ]
        assert min(h_dists) < 1.1  # the bridging OH proton

    @pytest.mark.parametrize("kind", ["Al", "Si", "Oss", "Osa", "Oaa"])
    def test_charge_bookkeeping_all_families(self, kind):
        cc = from_deck_cell(DECK, kind, metal_shells=2)
        c = cc.cluster
        q = (
            4 * sum(1 for s in c.symbols if s == "Si")
            + 3 * sum(1 for s in c.symbols if s == "Al")
            - 2 * sum(1 for s in c.symbols if s == "O")
            + sum(1 for s in c.symbols if s == "H")
        )
        assert q == 0 == c.charge
        assert min_interatomic_distance(c.coords) > 0.75

    def test_target_charge_deprotonation(self):
        cc = from_deck_cell(DECK, "Oss", metal_shells=1, target_charge=-1)
        assert cc.cluster.charge == -1
        assert any("deprotonated" in line for line in cc.termination_log)

    def test_frozen_shell_peripheral(self):
        cc = from_deck_cell(DECK, "Oss", metal_shells=2)
        c = cc.cluster
        assert c.frozen_indices  # lattice resistance present
        frozen = set(c.frozen_indices)
        # Reactive center free, and everything bonded to the attacked metal free.
        assert cc.attacked_index not in frozen
        assert cc.bridge_index not in frozen
        m = c.coords[cc.attacked_index]
        for i, s in enumerate(c.symbols):
            if s == "O" and np.linalg.norm(c.coords[i] - m) < 2.3:
                assert i not in frozen
        # Frozen atoms are exactly the outermost-shell heavy atoms;
        # protons (constructed guesses, not crystallography) never freeze.
        assert all(cc.atom_shell[i] == cc.metal_shells for i in frozen)
        assert all(c.symbols[i] != "H" for i in frozen)
        # And they are geometrically peripheral.
        center = c.coords[cc.bridge_index]
        d_frozen = [np.linalg.norm(c.coords[i] - center) for i in frozen]
        d_free_metal = [
            np.linalg.norm(c.coords[i] - center)
            for i in range(len(c.symbols))
            if i not in frozen and c.symbols[i] in ("Si", "Al")
        ]
        assert np.mean(d_frozen) > np.mean(d_free_metal)

    def test_n_intact_ladder_on_si(self):
        # Q-style connectivity rungs on the Si under attack (Oss center).
        for n in (1, 2, 3, 4):
            cc = from_deck_cell(DECK, "Oss", metal_shells=2, n_intact=n)
            c = cc.cluster
            assert cc.n_intact == n
            assert c.charge == 0
            assert min_interatomic_distance(c.coords) > 0.75
            # Attacked Si keeps full O coordination: pruned bridges
            # became terminal hydroxyls, not vacancies.
            # Thresholds sized to the crude legacy cell (bonds 1.4-2.8 A,
            # next-nearest O at 2.8+).
            m = c.coords[cc.attacked_index]
            o_near = [
                i
                for i, s in enumerate(c.symbols)
                if s == "O" and np.linalg.norm(c.coords[i] - m) < 2.3
            ]
            assert len(o_near) == 4
            # Exactly n of those oxygens still bridge to another metal.
            bridging = 0
            for i in o_near:
                other = [
                    j
                    for j, s in enumerate(c.symbols)
                    if s in ("Si", "Al")
                    and j != cc.attacked_index
                    and np.linalg.norm(c.coords[j] - c.coords[i]) < 2.8
                ]
                bridging += bool(other)
            assert bridging == n

    def test_pruning_shrinks_cluster(self):
        full = from_deck_cell(DECK, "Oss", metal_shells=2)
        pruned = from_deck_cell(DECK, "Oss", metal_shells=2, n_intact=1)
        n_si = lambda cc: sum(1 for s in cc.cluster.symbols if s == "Si")  # noqa: E731
        assert n_si(pruned) < n_si(full)

    @pytest.mark.parametrize("kind", ["Oss", "Osa", "Oaa"])
    def test_default_size_in_survey_window(self, kind):
        cc = from_deck_cell(DECK, kind, metal_shells=2)
        n_centers = sum(1 for s in cc.cluster.symbols if s in ("Si", "Al"))
        assert 4 <= n_centers <= 30
        assert 20 <= len(cc.cluster.symbols) <= 160

    def test_metal_center_family(self):
        cc = from_deck_cell(DECK, "Si", metal_shells=2)
        c = cc.cluster
        assert cc.bridge_index is None
        assert c.symbols[cc.attacked_index] == "Si"
        assert cc.atom_shell[cc.attacked_index] == 0
        assert c.charge == 0


class TestAttackComplex:
    def test_water_attack_on_oss(self):
        cc = from_deck_cell(DECK, "Oss", metal_shells=2)
        complex_, ow = attack_complex(cc, water())
        n0 = len(cc.cluster.symbols)
        assert ow == n0
        assert len(complex_.symbols) == n0 + 3
        assert complex_.symbols[ow] == "O"
        # Attacker O at the approach distance from the attacked metal.
        d = np.linalg.norm(complex_.coords[ow] - complex_.coords[cc.attacked_index])
        assert d == pytest.approx(3.2, abs=0.05)
        # Frozen shell carried through the merge untouched.
        assert complex_.frozen_indices == cc.cluster.frozen_indices
        assert min_interatomic_distance(complex_.coords) > 0.75
