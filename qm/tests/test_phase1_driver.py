"""Fast orchestration gates for the Xiao & Lasaga Phase-1 driver."""

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from quarry.clusters import (
    Cluster,
    aluminosilicate_dimer,
    disilicate,
    protonated_bridge_complex,
)
from quarry.pipeline import DftSettings, FrequencyResult
from quarry.rates import Thermo
from scripts import phase1_xiao_lasaga as phase1

CHEAP = DftSettings(xc="hf", basis="sto-3g")


def geometry(name: str, si_ow: float, si_obr: float = 1.6) -> Cluster:
    return Cluster(
        name=name,
        symbols=["O", "Si", "O", "H", "H"],
        coords=np.array(
            [
                [si_obr, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [si_ow, 0.0, 0.0],
                [si_ow, 0.9, 0.0],
                [si_ow, -0.9, 0.0],
            ]
        ),
    )


def _place_attacking_water(cluster: Cluster, ow_index: int, distance: float) -> Cluster:
    """Move the residual water as a rigid group without moving the bridge proton."""
    coords = cluster.coords.copy()
    direction = coords[ow_index] - coords[phase1.SI_INDEX]
    direction /= np.linalg.norm(direction)
    target = coords[phase1.SI_INDEX] + distance * direction
    delta = target - coords[ow_index]
    for index in (ow_index, ow_index + 2, ow_index + 3):
        coords[index] += delta
    return replace(cluster, coords=coords)


def acid_geometries(
    dimer_factory=disilicate,
) -> tuple[Cluster, Cluster, int, tuple[int, ...]]:
    """Full benchmark acid reactant/product with the opposite bridge intact."""
    dimer = dimer_factory()
    reactant = protonated_bridge_complex(dimer)
    ow_index = len(dimer.symbols)
    proton_indices = (ow_index + 1, ow_index + 2, ow_index + 3)
    product = _place_attacking_water(reactant, ow_index, 1.70)
    coords = product.coords.copy()

    bridge = coords[phase1.BR_INDEX]
    silicon = coords[phase1.SI_INDEX]
    direction = bridge - silicon
    direction /= np.linalg.norm(direction)
    target_bridge = silicon + 3.50 * direction
    leaving_delta = target_bridge - bridge
    # Atom 2 is the center across the bridge; 9:15 are its three OH groups.
    for index in (phase1.BR_INDEX, phase1.AL_INDEX, *range(9, 15), ow_index + 1):
        coords[index] += leaving_delta
    # Hydrolysis leaves Si-OH plus a doubly protonated bridge oxygen.  Move one
    # of the attacker's two hydrogens onto the departing bridge group.
    coords[ow_index + 2] = coords[phase1.BR_INDEX] + np.array([0.0, 0.98, 0.0])
    product = replace(product, name="acid-product", coords=coords)
    return reactant, product, ow_index, proton_indices


def test_channel_escape_guard_is_absolute_and_relative():
    guess = geometry("guess", 2.2)
    outside = phase1.channel_escape_reason(guess, geometry("far", 2.7), 2)
    relative = phase1.channel_escape_reason(
        geometry("tight-guess", 1.8), geometry("drift", 2.4), 2
    )

    assert outside is not None and "outside" in outside
    assert relative is not None and "escaped" in relative
    assert phase1.channel_escape_reason(guess, geometry("channel", 2.4), 2) is None


def test_acid_addition_guard_allows_loose_association_crest():
    guess = geometry("acid-addition-guess", 3.43)
    loose = geometry("acid-addition-ts", 3.63)
    escaped = geometry("escaped", 4.30)

    assert (
        phase1.addition_channel_escape_reason(guess, loose, 2, acid_path=True) is None
    )
    assert phase1.addition_channel_escape_reason(guess, escaped, 2, acid_path=True)


def test_quick_irc_requires_reactant_and_hydrolyzed_basins():
    reactant = geometry("reactant", 3.2, si_obr=1.6)
    reactant.coords[3:] += np.array([0.0, 3.0, 0.0])
    product = geometry("product", 1.7, si_obr=3.5)
    product.coords[3] = np.array([3.5, 0.98, 0.0])
    trapped = geometry("trapped", 1.8, si_obr=1.8)
    trapped.coords[3] = np.array([1.8, 0.98, 0.0])

    assert phase1.quick_irc_channel_reason(reactant, product, 2) is None
    reason = phase1.quick_irc_channel_reason(trapped, trapped, 2)
    assert reason is not None and "basin signatures" in reason


@pytest.mark.parametrize("dimer_factory", [disilicate, aluminosilicate_dimer])
def test_acid_quick_irc_requires_full_pre_equilibrium_and_product_speciation(
    dimer_factory,
):
    reactant, product, ow_index, proton_indices = acid_geometries(dimer_factory)

    assert (
        phase1.acid_basin_signature(reactant, ow_index, proton_indices)
        == phase1.ACID_REACTANT_BASIN
    )
    assert (
        phase1.acid_basin_signature(product, ow_index, proton_indices)
        == phase1.ACID_PRODUCT_BASIN
    )
    assert (
        phase1.acid_quick_irc_channel_reason(
            reactant,
            product,
            ow_index,
            proton_indices,
        )
        is None
    )

    intermediate = _place_attacking_water(reactant, ow_index, 1.80)
    terminal_product = replace(product, coords=product.coords.copy())
    terminal_product.coords[ow_index + 2] = terminal_product.coords[3] + np.array(
        [0.0, 0.0, 0.98]
    )
    assert (
        phase1.acid_basin_signature(terminal_product, ow_index, proton_indices)
        in phase1.ACID_HYDROLYZED_BASINS
    )
    assert (
        phase1.acid_quick_irc_addition_reason(
            reactant, intermediate, ow_index, proton_indices
        )
        is None
    )
    assert (
        phase1.acid_quick_irc_cleavage_reason(
            intermediate, terminal_product, ow_index, proton_indices
        )
        is None
    )
    assert (
        phase1.acid_endpoint_in_basins(
            intermediate,
            terminal_product,
            ow_index,
            proton_indices,
            phase1.ACID_HYDROLYZED_BASINS,
        )
        is terminal_product
    )

    dissociated_product = replace(
        terminal_product, coords=terminal_product.coords.copy()
    )
    dissociated_product.coords[ow_index + 2] = np.array([100.0, 100.0, 100.0])
    assert (
        phase1.acid_basin_signature(dissociated_product, ow_index, proton_indices)
        not in phase1.ACID_HYDROLYZED_BASINS
    )
    assert phase1.acid_quick_irc_cleavage_reason(
        intermediate, dissociated_product, ow_index, proton_indices
    )


@pytest.mark.parametrize("dimer_factory", [disilicate, aluminosilicate_dimer])
def test_acid_speciation_guards_reject_lost_proton_bridge_or_water(dimer_factory):
    reactant, product, ow_index, proton_indices = acid_geometries(dimer_factory)

    lost_proton = replace(reactant, coords=reactant.coords.copy())
    lost_proton.coords[proton_indices[0]] += np.array([0.0, 2.0, 0.0])
    assert phase1.protonated_bridge_reason(lost_proton, ow_index, proton_indices)

    broken_opposite = replace(reactant, coords=reactant.coords.copy())
    broken_opposite.coords[phase1.AL_INDEX] += np.array([3.0, 0.0, 0.0])
    assert phase1.protonated_bridge_reason(broken_opposite, ow_index, proton_indices)

    hydroxide_attacker = replace(reactant, coords=reactant.coords.copy())
    hydroxide_attacker.coords[ow_index + 3] += np.array([0.0, 0.0, 2.0])
    assert phase1.protonated_bridge_reason(hydroxide_attacker, ow_index, proton_indices)

    deprotonated_product = replace(product, coords=product.coords.copy())
    deprotonated_product.coords[proton_indices[0]] += np.array([0.0, 2.0, 0.0])
    assert phase1.acid_product_reason(deprotonated_product, ow_index, proton_indices)


def test_sequential_quick_irc_gates_each_elementary_step():
    reactant = geometry("reactant", 3.2, si_obr=1.6)
    reactant.coords[3:] += np.array([0.0, 3.0, 0.0])
    intermediate = geometry("intermediate", 1.8, si_obr=1.8)
    intermediate.coords[3] = np.array([1.8, 0.98, 0.0])
    product = geometry("product", 1.7, si_obr=3.5)
    product.coords[3] = np.array([3.5, 0.98, 0.0])

    assert phase1.quick_irc_addition_reason(reactant, intermediate, 2) is None
    assert phase1.quick_irc_cleavage_reason(intermediate, product, 2) is None
    assert phase1.quick_irc_cleavage_reason(reactant, product, 2) is not None
    assert phase1.quick_irc_addition_reason(intermediate, product, 2) is not None
    assert phase1.quick_irc_cleavage_reason(reactant, intermediate, 2) is not None


def test_endpoint_in_basin_requires_exactly_one_match():
    reactant = geometry("reactant", 3.2, si_obr=1.6)
    reactant.coords[3:] += np.array([0.0, 3.0, 0.0])
    intermediate = geometry("intermediate", 1.8, si_obr=1.8)
    intermediate.coords[3] = np.array([1.8, 0.98, 0.0])

    assert (
        phase1.endpoint_in_basin(
            reactant,
            intermediate,
            2,
            phase1.ASSOCIATIVE_BASIN,
        )
        is intermediate
    )
    with np.testing.assert_raises_regex(RuntimeError, "found 0"):
        phase1.endpoint_in_basin(
            reactant,
            intermediate,
            2,
            phase1.HYDROLYZED_BASIN,
        )


def test_sequential_barrier_metrics_separate_local_and_profile_barriers():
    def thermo(gibbs_kj: float) -> Thermo:
        return Thermo(gibbs_kj, 0.0, 0.0, 0.0, 298.15)

    metrics = phase1.sequential_barrier_metrics(
        reactant=thermo(0.0),
        intermediate=thermo(-20.0),
        addition_ts=thermo(70.0),
        cleavage_ts=thermo(100.0),
    )

    assert metrics["addition_dG_dagger_kj"] == 70.0
    assert metrics["cleavage_dG_dagger_kj"] == 120.0
    assert metrics["overall_profile_dG_dagger_kj"] == 100.0
    assert metrics["rate_limiting_local_step"] == "cleavage"
    assert metrics["highest_profile_ts"] == "cleavage"


def test_minimum_and_comparison_helpers_fail_closed():
    minimum = FrequencyResult(
        frequencies_cm=np.array([300.0]),
        imaginary_cm=np.array([]),
        electronic_hartree=-1.0,
        molar_mass_kg=0.1,
        rotational_temperatures_k=(1.0, 2.0, 3.0),
        linear=False,
    )
    unstable = replace(minimum, imaginary_cm=np.array([50.0]))

    assert phase1.minimum_frequency_reason("reactant", minimum) is None
    reason = phase1.minimum_frequency_reason("reactant", unstable)
    assert reason is not None and "1 imaginary" in reason
    comparison = phase1.comparison_against_si_neutral(
        "al-acid", acid_path=True, barrier_kj=82.193
    )
    assert comparison["acid_lower_than_si_neutral"] is True


def test_basin_equivalence_rejects_same_signature_different_conformer():
    from dataclasses import replace

    canonical = geometry("intermediate", 1.8, si_obr=1.8)
    canonical.coords[3] = np.array([1.8, 0.98, 0.0])
    changed_coords = canonical.coords.copy()
    changed_coords[4] += np.array([0.0, 0.0, 1.5])
    candidate = replace(canonical, coords=changed_coords)

    assert phase1.hydrolysis_basin_signature(candidate, 2) == (phase1.ASSOCIATIVE_BASIN)
    reason = phase1.basin_equivalence_reason(
        candidate,
        canonical,
        [0, 1, 2, 3, 4],
        candidate_energy_hartree=-100.0,
        canonical_energy_hartree=-100.0,
    )
    assert reason is not None and "RMSD" in reason


def test_begin_sequential_run_quarantines_stale_final_outputs(tmp_path):
    (tmp_path / "results.json").write_text("stale")
    (tmp_path / "store.sqlite").write_text("stale")

    phase1.begin_sequential_run(tmp_path)

    assert not (tmp_path / "results.json").exists()
    assert not (tmp_path / "store.sqlite").exists()
    assert len(list(tmp_path.glob("results.superseded-*.json"))) == 1
    assert len(list(tmp_path.glob("store.superseded-*.sqlite"))) == 1
    status = json.loads((tmp_path / "run_status.json").read_text())
    assert status["status"] == "running"


def test_block_run_quarantines_stale_success_and_records_reason(tmp_path):
    (tmp_path / "results.json").write_text("stale")
    (tmp_path / "store.sqlite").write_text("stale")

    phase1.block_run(tmp_path, "reactant is not a minimum")

    assert not (tmp_path / "results.json").exists()
    assert not (tmp_path / "store.sqlite").exists()
    assert len(list(tmp_path.glob("results.blocked-*.json"))) == 1
    assert len(list(tmp_path.glob("store.blocked-*.sqlite"))) == 1
    status = json.loads((tmp_path / "run_status.json").read_text())
    assert status == {
        "status": "blocked",
        "updated_at": status["updated_at"],
        "detail": "reactant is not a minimum",
    }


@pytest.mark.parametrize(
    ("terminal_failure", "expected_detail", "expected_gate_calls"),
    [
        pytest.param(
            "channel-escape",
            "saddle r(Si-Ow)=3.00 A is outside the bonding channel; "
            "acid-neb route exhausted before thermochemistry",
            (0, 0),
            id="channel-escape",
        ),
        pytest.param(
            "saddle",
            "not a clean first-order saddle — inspect ts.xyz",
            (1, 0),
            id="non-first-order-saddle",
        ),
        pytest.param(
            "quick-irc",
            "quick IRC missed the product basin; aborting before thermochemistry",
            (1, 1),
            id="quick-irc",
        ),
    ],
)
def test_acid_main_records_terminal_stage_failures(
    monkeypatch, tmp_path, terminal_failure, expected_detail, expected_gate_calls
):
    dimer = disilicate()
    reactant = protonated_bridge_complex(dimer)
    ow_index = len(dimer.symbols)
    approach = _place_attacking_water(reactant, ow_index, 1.90)
    saddle = _place_attacking_water(reactant, ow_index, 2.00)
    escaped = _place_attacking_water(reactant, ow_index, 3.00)

    monkeypatch.setattr(
        phase1.sys,
        "argv",
        ["phase1_xiao_lasaga.py", "--reaction", "si-acid"],
    )
    monkeypatch.setattr(phase1, "reaction_run_slug", lambda *args: tmp_path)
    monkeypatch.setattr(phase1, "optimize", lambda cluster, settings: cluster)
    monkeypatch.setattr(
        phase1,
        "scan_to_maximum",
        lambda *args, **kwargs: [(1.90, -1.0, approach)],
    )
    monkeypatch.setattr(phase1, "acid_neb_guess", lambda *args, **kwargs: saddle)
    monkeypatch.setattr(
        phase1,
        "refine_phase1_saddle",
        lambda *args, **kwargs: (
            escaped if terminal_failure == "channel-escape" else saddle
        ),
    )
    frequency = FrequencyResult(
        frequencies_cm=np.array([300.0, 700.0, 1100.0]),
        imaginary_cm=np.array([500.0] if terminal_failure == "quick-irc" else []),
        electronic_hartree=-100.0,
        molar_mass_kg=0.1,
        rotational_temperatures_k=(1.0, 2.0, 3.0),
        linear=False,
        imaginary_mode=(
            np.ones_like(saddle.coords) if terminal_failure == "quick-irc" else None
        ),
    )
    gate_calls = {"frequencies": 0, "quick_irc": 0}

    def mocked_frequencies(cluster, settings):
        gate_calls["frequencies"] += 1
        return frequency

    def mocked_quick_irc(*args, **kwargs):
        gate_calls["quick_irc"] += 1
        return reactant, reactant

    monkeypatch.setattr(phase1, "frequencies", mocked_frequencies)
    monkeypatch.setattr(phase1, "quick_irc", mocked_quick_irc)
    monkeypatch.setattr(
        phase1,
        "acid_quick_irc_channel_reason",
        lambda *args, **kwargs: "quick IRC missed the product basin",
    )
    monkeypatch.setattr(
        phase1,
        "acid_quick_irc_cleavage_reason",
        lambda *args, **kwargs: "quick IRC missed the product basin",
    )

    assert phase1.main() == 1
    status = json.loads((tmp_path / "run_status.json").read_text())
    assert status["status"] == "blocked"
    assert status["detail"] == expected_detail
    assert (gate_calls["frequencies"], gate_calls["quick_irc"]) == expected_gate_calls


def test_sequential_closeout_writes_gated_profile(monkeypatch, tmp_path):
    reactant = geometry("complex", 3.2, si_obr=1.6)
    reactant.coords[3:] += np.array([0.0, 3.0, 0.0])
    intermediate = geometry("intermediate", 1.8, si_obr=1.8)
    intermediate.coords[3] = np.array([1.8, 0.98, 0.0])
    product = geometry("product", 1.7, si_obr=3.5)
    product.coords[3] = np.array([3.5, 0.98, 0.0])
    addition_guess = geometry("addition-guess", 2.0, si_obr=1.7)
    addition_ts = geometry("addition-ts", 2.0, si_obr=1.7)
    cleavage_ts = geometry("cleavage-ts", 1.8, si_obr=2.0)

    def frequency(cluster, electronic, imaginary=()):
        return FrequencyResult(
            frequencies_cm=np.array([300.0, 700.0, 1100.0]),
            imaginary_cm=np.array(imaginary),
            electronic_hartree=electronic,
            molar_mass_kg=0.1,
            rotational_temperatures_k=(1.0, 2.0, 3.0),
            linear=False,
            imaginary_mode=(
                np.ones_like(cluster.coords) if len(imaginary) == 1 else None
            ),
        )

    frequencies_by_name = {
        "complex": frequency(reactant, -100.00),
        "intermediate": frequency(intermediate, -100.01),
        "addition-ts": frequency(addition_ts, -99.96, (500.0,)),
        "cleavage-ts": frequency(cleavage_ts, -99.95, (300.0,)),
    }

    monkeypatch.setattr(phase1, "neb_ts_guess", lambda *args, **kwargs: addition_guess)
    relaxed_calls = []

    def fake_relax(cluster, settings, **kwargs):
        relaxed_calls.append((cluster, kwargs))
        return addition_guess

    monkeypatch.setattr(phase1, "relax_at_fixed_distances", fake_relax)

    def fake_find_ts(cluster, settings, **kwargs):
        assert cluster is addition_guess
        assert kwargs["internal"] is False
        assert kwargs["initial_mode"].shape == addition_guess.coords.shape
        return addition_ts

    monkeypatch.setattr(phase1, "find_ts", fake_find_ts)
    monkeypatch.setattr(
        phase1,
        "quick_irc",
        lambda *args, **kwargs: (reactant, intermediate),
    )
    monkeypatch.setattr(
        phase1,
        "frequencies",
        lambda cluster, settings: frequencies_by_name[cluster.name],
    )
    monkeypatch.setattr(
        phase1,
        "energy",
        lambda cluster, settings: frequencies_by_name[cluster.name].electronic_hartree,
    )

    class FakeStore:
        def __init__(self, path):
            self.path = Path(path)
            self.path.touch()

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def add_structure(self, *args, **kwargs):
            return 1

        def add_job(self, *args, **kwargs):
            return 1

        def set_job_status(self, *args, **kwargs):
            return None

        def add_result(self, *args, **kwargs):
            return None

    monkeypatch.setattr(phase1, "Store", FakeStore)

    (tmp_path / "results.json").write_text("stale")
    (tmp_path / "store.sqlite").write_text("stale")
    status = phase1.finish_al_neutral_sequential(
        complex_opt=reactant,
        cleavage_ts=cleavage_ts,
        cleavage_freq=frequencies_by_name["cleavage-ts"],
        cleavage_back=intermediate,
        cleavage_fwd=product,
        settings=CHEAP,
        run_dir=tmp_path,
        ow_index=2,
        temperature=298.15,
        dimer_opt=reactant,
        attacker_opt=reactant,
    )

    results = json.loads((tmp_path / "results.json").read_text())
    assert status == 0
    assert len(relaxed_calls) == 1
    assert len(relaxed_calls[0][1]["fixed_distances"]) == 2
    assert results["mechanism"] == "sequential-associative"
    assert set(results["steps"]) == {"addition", "cleavage"}
    assert results["profile"]["highest_profile_ts"] == "cleavage"
    assert (tmp_path / "store.sqlite").exists()
    run_status = json.loads((tmp_path / "run_status.json").read_text())
    assert run_status["status"] == "completed"
    assert len(list(tmp_path.glob("results.superseded-*.json"))) == 1
    assert len(list(tmp_path.glob("store.superseded-*.sqlite"))) == 1


def test_proton_neb_fallback_reconstructs_resumed_approach_seed(monkeypatch, tmp_path):
    direct_crest = geometry("direct-crest", 2.2)
    complex_opt = geometry("complex", 3.2)
    calls = []

    def fake_constrained_scan(
        cluster,
        settings,
        *,
        atom_i,
        atom_j,
        distances_a,
        fixed_distances=(),
        **kwargs,
    ):
        calls.append((cluster.name, atom_i, atom_j, distances_a, fixed_distances))
        if (atom_i, atom_j) == (phase1.SI_INDEX, 2):
            return [(1.90, -1.0, geometry("approach-r1.90", 1.90))]
        broken = geometry("broken", 1.70, si_obr=2.60)
        return [(r, -2.0, broken) for r in distances_a]

    def fake_proton_scan(cluster, settings, **kwargs):
        points = [
            (1.40, 0.0, geometry("proton-start", 1.90)),
            (1.25, 2.0, geometry("proton-crest", 1.90)),
            (1.10, 1.0, geometry("proton-product-seed", 1.90)),
        ]
        return points

    neb_inputs = []

    def fake_neb(reactant, product, settings, **kwargs):
        neb_inputs.append((reactant, product))
        assert kwargs == {
            "n_images": 5,
            "fmax_ev_a": 0.2,
            "max_steps": 240,
            "pre_relax_fmax_ev_a": 1.5,
            "pre_relax_steps": 120,
        }
        return geometry("neb-guess", 2.0)

    monkeypatch.setattr(phase1, "constrained_scan", fake_constrained_scan)
    monkeypatch.setattr(phase1, "scan_to_maximum", fake_proton_scan)
    monkeypatch.setattr(phase1, "optimize", lambda cluster, settings: cluster)
    monkeypatch.setattr(phase1, "neb_ts_guess", fake_neb)

    result = phase1.proton_neb_guess(direct_crest, complex_opt, CHEAP, tmp_path, 2)

    assert result.name == "neb-guess"
    assert calls[0][1:4] == (phase1.SI_INDEX, 2, [1.90])
    assert neb_inputs[0][0] is complex_opt
    assert (
        np.linalg.norm(neb_inputs[0][1].coords[1] - neb_inputs[0][1].coords[2]) == 1.70
    )
    assert (tmp_path / "product.xyz").exists()


def test_rolled_back_product_extends_cleavage_without_repeating_proton_scan(
    monkeypatch, tmp_path
):
    direct_crest = geometry("direct-crest", 2.2)
    complex_opt = geometry("complex", 3.2)
    phase1.save_xyz(
        geometry("rolled-back", 1.70, si_obr=2.07), tmp_path / "product.xyz"
    )
    break_targets = []

    def fake_break_scan(cluster, settings, *, distances_a, **kwargs):
        break_targets.extend(distances_a)
        broken = geometry("fully-broken", 1.70, si_obr=3.60)
        return [(r, -2.0, broken) for r in distances_a]

    monkeypatch.setattr(phase1, "constrained_scan", fake_break_scan)
    monkeypatch.setattr(
        phase1,
        "scan_to_maximum",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("proton scan must not repeat")
        ),
    )
    monkeypatch.setattr(phase1, "optimize", lambda cluster, settings: cluster)
    monkeypatch.setattr(
        phase1,
        "neb_ts_guess",
        lambda reactant, product, settings, **kwargs: geometry("neb-guess", 2.0),
    )

    result = phase1.proton_neb_guess(direct_crest, complex_opt, CHEAP, tmp_path, 2)

    assert result.name == "neb-guess"
    assert break_targets == [2.3, 2.6, 3.0, 3.4, 3.6]
    assert (tmp_path / "product.rejected-rollback.xyz").exists()


def test_completed_product_uses_full_aligned_neb_endpoint(monkeypatch, tmp_path):
    approach = geometry("approach", 2.2)
    complex_opt = geometry("complex", 3.2)
    phase1.save_xyz(geometry("full-product", 1.66, 3.50), tmp_path / "product.xyz")
    near = geometry("near-product", 1.67, 2.07)
    near.coords[3] = np.array([2.07, 0.96, 0.0])
    phase1.save_xyz(near, tmp_path / "product.rejected-rollback.xyz")
    endpoints = []

    def fake_neb(reactant, product, settings, **kwargs):
        endpoints.append((product, kwargs))
        return geometry("neb-guess", 2.0)

    monkeypatch.setattr(phase1, "neb_ts_guess", fake_neb)

    result = phase1.proton_neb_guess(approach, complex_opt, CHEAP, tmp_path, 2)

    assert result.name == "neb-guess"
    assert len(endpoints) == 1
    product, kwargs = endpoints[0]
    assert np.linalg.norm(product.coords[1] - product.coords[0]) == 3.50
    assert kwargs == {
        "n_images": 5,
        "fmax_ev_a": 0.2,
        "max_steps": 240,
        "pre_relax_fmax_ev_a": 1.5,
        "pre_relax_steps": 120,
    }


def test_acid_neb_builds_product_without_repeating_proton_transfer(
    monkeypatch, tmp_path
):
    reactant, product, ow_index, proton_indices = acid_geometries()
    approach = _place_attacking_water(reactant, ow_index, 1.90)
    break_targets = []
    fixed = []

    def fake_break_scan(
        cluster, settings, *, distances_a, fixed_distances=(), **kwargs
    ):
        break_targets.extend(distances_a)
        fixed.extend(fixed_distances)
        return [(distance, -2.0, product) for distance in distances_a]

    endpoints = []

    def fake_neb(reactant_endpoint, product_endpoint, settings, **kwargs):
        endpoints.append((reactant_endpoint, product_endpoint, kwargs))
        return replace(approach, name="acid-neb-guess")

    monkeypatch.setattr(phase1, "constrained_scan", fake_break_scan)
    monkeypatch.setattr(phase1, "optimize", lambda cluster, settings: cluster)
    monkeypatch.setattr(phase1, "neb_ts_guess", fake_neb)
    monkeypatch.setattr(
        phase1,
        "scan_to_maximum",
        lambda *args, **kwargs: (_ for _ in ()).throw(
            AssertionError("pre-equilibrated acid path must not scan proton transfer")
        ),
    )

    result = phase1.acid_neb_guess(
        approach,
        reactant,
        CHEAP,
        tmp_path,
        ow_index=ow_index,
        proton_indices=proton_indices,
    )

    assert result.name == "acid-neb-guess"
    assert break_targets == [1.85, 2.05, 2.3, 2.6, 3.0, 3.4, 3.6]
    assert (phase1.SI_INDEX, ow_index, 1.70) in fixed
    assert (phase1.BR_INDEX, proton_indices[0], 0.96) in fixed
    assert (tmp_path / "product.xyz").exists()
    assert len(endpoints) == 1
    assert endpoints[0][0] is reactant
    assert endpoints[0][2]["climb_optimizer"] == "ode"


def test_acid_neb_resumes_saved_product_after_band_failure(monkeypatch, tmp_path):
    reactant, product, ow_index, proton_indices = acid_geometries()
    approach = _place_attacking_water(reactant, ow_index, 1.90)
    break_calls = []
    neb_calls = []

    def fake_break_scan(cluster, settings, *, distances_a, **kwargs):
        break_calls.append(tuple(distances_a))
        return [(distance, -2.0, product) for distance in distances_a]

    def flaky_neb(reactant_endpoint, product_endpoint, settings, **kwargs):
        neb_calls.append(product_endpoint)
        if len(neb_calls) == 1:
            raise RuntimeError("bounded CI-NEB failed")
        return replace(approach, name="resumed-acid-neb")

    monkeypatch.setattr(phase1, "constrained_scan", fake_break_scan)
    monkeypatch.setattr(phase1, "optimize", lambda cluster, settings: cluster)
    monkeypatch.setattr(phase1, "neb_ts_guess", flaky_neb)

    with pytest.raises(RuntimeError, match="CI-NEB failed"):
        phase1.acid_neb_guess(
            approach,
            reactant,
            CHEAP,
            tmp_path,
            ow_index,
            proton_indices,
        )
    assert (tmp_path / "product.xyz").exists()

    resumed = phase1.acid_neb_guess(
        reactant,
        reactant,
        CHEAP,
        tmp_path,
        ow_index,
        proton_indices,
    )
    assert resumed.name == "resumed-acid-neb"
    assert len(break_calls) == 1
    assert len(neb_calls) == 2


def test_acid_neb_rejects_cached_product_without_reusing_it_as_seed(
    monkeypatch, tmp_path
):
    reactant, product, ow_index, proton_indices = acid_geometries()
    approach = _place_attacking_water(reactant, ow_index, 1.90)
    invalid_cached = replace(approach, coords=approach.coords.copy())
    invalid_cached.coords[proton_indices[0]] = np.array([100.0, 100.0, 100.0])
    phase1.save_xyz(invalid_cached, tmp_path / "product.xyz")
    seen_seeds = []

    def fake_break_scan(cluster, settings, *, distances_a, **kwargs):
        seen_seeds.append(cluster)
        return [(distance, -2.0, product) for distance in distances_a]

    monkeypatch.setattr(phase1, "constrained_scan", fake_break_scan)
    monkeypatch.setattr(phase1, "optimize", lambda cluster, settings: cluster)
    monkeypatch.setattr(
        phase1,
        "neb_ts_guess",
        lambda *args, **kwargs: replace(approach, name="recovered-acid-neb"),
    )

    result = phase1.acid_neb_guess(
        approach,
        reactant,
        CHEAP,
        tmp_path,
        ow_index,
        proton_indices,
    )

    assert result.name == "recovered-acid-neb"
    assert seen_seeds[0] is approach
    rebuilt = phase1.load_xyz(tmp_path / "product.xyz", product)
    assert phase1.acid_product_reason(rebuilt, ow_index, proton_indices) is None
    assert list(tmp_path.glob("product.rejected-acid-rollback-*.xyz"))


def test_acid_saddle_refinement_uses_endpoint_directed_cartesian_mode(
    monkeypatch, tmp_path
):
    reactant, product, ow_index, proton_indices = acid_geometries()
    guess = _place_attacking_water(reactant, ow_index, 2.0)
    phase1.save_xyz(product, tmp_path / "product.xyz")
    calls = []

    def fake_find_ts(cluster, settings, **kwargs):
        calls.append(kwargs)
        return cluster

    monkeypatch.setattr(phase1, "find_ts", fake_find_ts)
    result = phase1.refine_phase1_saddle(
        guess,
        reactant,
        CHEAP,
        trajectory_path=tmp_path / "sella.traj",
        acid_path=True,
        route="acid-neb",
        run_dir=tmp_path,
        ow_index=ow_index,
        proton_indices=proton_indices,
    )

    assert result is guess
    assert len(calls) == 1
    assert calls[0]["internal"] is False
    assert calls[0]["initial_mode"].shape == reactant.coords.shape
    assert np.linalg.norm(calls[0]["initial_mode"]) == pytest.approx(1.0)
    assert calls[0]["trajectory"] == str(tmp_path / "sella.traj")


def test_acid_addition_scan_guess_persists_structures_and_manifest(
    monkeypatch, tmp_path
):
    reactant = geometry("reactant", 3.2)
    crest = geometry("crest", 2.6)
    product_side = geometry("product-side", 2.4)
    points = [
        (2.8, -2.0, geometry("loose", 2.8)),
        (2.6, -1.0, crest),
        (2.4, -2.0, product_side),
    ]
    monkeypatch.setattr(phase1, "scan_to_maximum", lambda *args, **kwargs: points)

    result = phase1.acid_addition_scan_guess(reactant, CHEAP, tmp_path, 2)

    assert result is crest
    manifest = json.loads((tmp_path / "addition_scan" / "scan.json").read_text())
    assert [point["distance_a"] for point in manifest] == [2.8, 2.6, 2.4]
    assert (tmp_path / "addition_scan" / "r-2.60.xyz").exists()


def test_acid_checkpoint_namespace_is_mechanism_qualified():
    assert phase1.reaction_run_slug("si-neutral", "b3lyp", "def2-svp", "flank") == (
        "si-neutral-b3lyp-def2-svp-flank"
    )
    assert phase1.reaction_run_slug("si-acid", "b3lyp", "def2-svp", "flank") == (
        "si-acid-preprotonated-v2-b3lyp-def2-svp-flank"
    )
