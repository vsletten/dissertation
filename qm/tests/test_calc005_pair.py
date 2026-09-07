"""CPU-only physical and provenance gates for the CALC-005 live states."""

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

import quarry.calc005 as calc005
from quarry.calc005 import build_calc005_pair, probe_calc005_pairs
from quarry.crystal import AtomOrigin

REPO = Path(__file__).resolve().parents[2]
DECK = REPO / "petra" / "examples" / "kaolinite.toml"

EXPECTED = {
    1: ("Al6H36O29Si", "Al6H38O30Si", "Al6H34O26"),
    2: ("Al6H38O36Si4", "Al6H42O38Si4", "Al6H38O34Si3"),
    3: ("Al7H43O46Si7", "Al7H49O49Si7", "Al7H45O45Si6"),
    4: ("Al7H43O52Si10", "Al7H51O56Si10", "Al7H47O52Si9"),
}


@pytest.mark.parametrize("index", [1, 2, 3, 4])
def test_live_pair_has_exact_state_expansion_and_balanced_formula(index):
    pair = build_calc005_pair(DECK, environment_index=index)
    seed_formula, occupied_formula, vacancy_formula = EXPECTED[index]

    assert pair.condensed.cluster.formula == seed_formula
    assert pair.occupied.formula == occupied_formula
    assert pair.vacancy.formula == vacancy_formula
    assert pair.silicic_acid.formula == "H4O4Si"
    assert pair.occupied.charge == pair.vacancy.charge == 0
    assert pair.silicic_acid.charge == 0
    assert pair.occupied.spin == pair.vacancy.spin == pair.silicic_acid.spin == 0
    assert pair.hydrolysis_water_count == index
    assert pair.validate()["atom_conservation"] is True


def test_pilot_states_match_live_adsorb_desorb_transition():
    pair = build_calc005_pair(DECK, environment_index=1)

    assert pair.x == 1
    assert pair.y == 0
    assert pair.occupied_states == {
        "center": 205,
        "osa": 402,
        "oss": [303, 303, 303],
    }
    assert pair.vacancy_states == {
        "center": 200,
        "osa": 404,
        "oss": [300, 300, 300],
    }


def test_pair_preserves_support_and_frozen_identity_without_distance_deletion():
    pair = build_calc005_pair(DECK, environment_index=1)
    report = pair.validate()

    assert report["support_origin_match"] is True
    assert report["frozen_origin_match"] is True
    assert report["frozen_coordinate_match"] is True
    assert report["four_center_oh_groups"] is True
    assert report["hydrolysis_origin_count"] == 3
    assert report["min_interatomic_distance_a"] >= 0.75
    assert report["max_nearest_oh_a"] <= 1.25
    assert report["min_owner_margin_a"] >= 0.15
    assert report["frozen_heavy_only"] is True
    assert report["frozen_peripheral"] is True


def test_validator_rejects_state_charge_spin_formula_and_nonfinite_corruption():
    pair = build_calc005_pair(DECK, environment_index=1)
    wrong_states = replace(pair, occupied_states={**pair.occupied_states, "osa": 403})
    wrong_charge = replace(pair, occupied=replace(pair.occupied, charge=7))
    wrong_spin = replace(pair, vacancy=replace(pair.vacancy, spin=4))
    wrong_symbols = list(pair.occupied.symbols)
    wrong_symbols[-1] = "O"
    wrong_formula = replace(
        pair, occupied=replace(pair.occupied, symbols=wrong_symbols)
    )
    bad_coords = pair.occupied.coords.copy()
    bad_coords[0, 0] = np.nan
    nonfinite = replace(pair, occupied=replace(pair.occupied, coords=bad_coords))

    assert wrong_states.validate()["exact_state_match"] is False
    assert wrong_charge.validate()["charge_spin_match"] is False
    assert wrong_spin.validate()["charge_spin_match"] is False
    assert wrong_formula.validate()["exact_formula_match"] is False
    assert nonfinite.validate()["finite_coordinates"] is False


def test_validator_rejects_hydrolysis_provenance_and_intended_owner_corruption():
    pair = build_calc005_pair(DECK, environment_index=1)
    origins = list(pair.occupied_origins)
    index = next(
        i
        for i, origin in enumerate(origins)
        if origin.kind == "hydrolysis-water-h-support"
    )
    origins[index] = AtomOrigin(
        "hydrolysis-water-h-support",
        (10, (0, 0, 0)),
        2,
    )
    corrupted = replace(pair, occupied_origins=tuple(origins))

    report = corrupted.validate()
    assert report["hydrolysis_origin_match"] is False
    assert report["intended_owner_match"] is False


def test_probe_is_deterministic_and_covers_every_rung():
    first = probe_calc005_pairs(DECK)
    second = probe_calc005_pairs(DECK)

    assert first == second
    assert [row["environment_index"] for row in first["pairs"]] == [1, 2, 3, 4]
    assert first["two_pass_deterministic"] is True
    assert len({row["occupied_geometry_hash"] for row in first["pairs"]}) == 4
    assert len({row["vacancy_geometry_hash"] for row in first["pairs"]}) == 4


def test_probe_refuses_nondeterministic_or_nondistinct_output(monkeypatch):
    original = calc005.build_calc005_pair
    calls = 0

    def unstable(*args, **kwargs):
        nonlocal calls
        calls += 1
        pair = original(*args, **kwargs)
        if calls > 4:
            return replace(pair, occupied=replace(pair.occupied, name="drifted"))
        return pair

    monkeypatch.setattr(calc005, "build_calc005_pair", unstable)
    with pytest.raises(RuntimeError, match="deterministic"):
        calc005.probe_calc005_pairs(DECK)

    monkeypatch.setattr(calc005, "build_calc005_pair", original)
    monkeypatch.setattr(calc005, "geometry_hash", lambda _xyz: "same")
    with pytest.raises(RuntimeError, match="distinct"):
        calc005.probe_calc005_pairs(DECK)


@pytest.mark.parametrize("bad", [0, 5, True, 1.0])
def test_pair_rejects_unsupported_or_aliased_environment_index(bad):
    with pytest.raises((TypeError, ValueError), match="environment_index"):
        build_calc005_pair(DECK, environment_index=bad)


def test_validator_refuses_coordinated_and_condensed_tampering():
    pair = build_calc005_pair(DECK, environment_index=2)

    seed_origins = list(pair.condensed.atom_origins)
    seed_origins[0] = replace(seed_origins[0], ordinal=99)
    bad_seed_origin = replace(
        pair,
        condensed=replace(pair.condensed, atom_origins=tuple(seed_origins)),
    )
    assert bad_seed_origin.validate()["canonical_rebuild_match"] is False

    changed = replace(pair.vacancy_origins[0], ordinal=88)
    vacancy_origins = list(pair.vacancy_origins)
    occupied_origins = list(pair.occupied_origins)
    vacancy_origins[0] = occupied_origins[0] = changed
    bad_provenance = replace(
        pair,
        vacancy_origins=tuple(vacancy_origins),
        occupied_origins=tuple(occupied_origins),
    )
    assert bad_provenance.validate()["canonical_rebuild_match"] is False

    wrong_vacancy_coords = pair.vacancy.coords.copy()
    h_index = next(
        i
        for i, origin in enumerate(pair.vacancy_origins)
        if origin.kind == "hydrolysis-water-h-support"
    )
    wrong_oxygen = next(
        i
        for i, (symbol, origin) in enumerate(
            zip(pair.vacancy.symbols, pair.vacancy_origins, strict=True)
        )
        if symbol == "O" and origin.node != pair.vacancy_origins[h_index].node
    )
    wrong_vacancy_coords[h_index] = pair.vacancy.coords[wrong_oxygen] + [0.96, 0, 0]
    bad_owner = replace(
        pair,
        vacancy=replace(pair.vacancy, coords=wrong_vacancy_coords),
    )
    assert bad_owner.validate()["intended_owner_match"] is False

    v_symbols = list(pair.vacancy.symbols)
    al_index = v_symbols.index("Al")
    si_index = v_symbols.index("Si")
    v_symbols[al_index], v_symbols[si_index] = v_symbols[si_index], v_symbols[al_index]
    c_symbols = list(pair.occupied.symbols)
    c_symbols[al_index], c_symbols[si_index] = c_symbols[si_index], c_symbols[al_index]
    bad_elements = replace(
        pair,
        vacancy=replace(pair.vacancy, symbols=v_symbols),
        occupied=replace(pair.occupied, symbols=c_symbols),
    )
    assert bad_elements.validate()["canonical_rebuild_match"] is False

    bad_seed_charge = replace(
        pair,
        condensed=replace(
            pair.condensed,
            cluster=replace(pair.condensed.cluster, charge=1),
        ),
    )
    assert bad_seed_charge.validate()["charge_spin_match"] is False

    seed_coords = pair.condensed.cluster.coords.copy()
    seed_coords[0, 0] = np.nan
    bad_seed_nan = replace(
        pair,
        condensed=replace(
            pair.condensed,
            cluster=replace(pair.condensed.cluster, coords=seed_coords),
        ),
    )
    assert bad_seed_nan.validate()["finite_coordinates"] is False
