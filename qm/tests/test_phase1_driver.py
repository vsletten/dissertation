"""Fast orchestration gates for the Xiao & Lasaga Phase-1 driver."""

import numpy as np

from quarry.clusters import Cluster
from quarry.pipeline import DftSettings
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


def test_channel_escape_guard_is_absolute_and_relative():
    guess = geometry("guess", 2.2)
    outside = phase1.channel_escape_reason(guess, geometry("far", 2.7), 2)
    relative = phase1.channel_escape_reason(
        geometry("tight-guess", 1.8), geometry("drift", 2.4), 2
    )

    assert outside is not None and "outside" in outside
    assert relative is not None and "escaped" in relative
    assert phase1.channel_escape_reason(guess, geometry("channel", 2.4), 2) is None


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


def test_sequential_quick_irc_gates_each_elementary_step():
    reactant = geometry("reactant", 3.2, si_obr=1.6)
    reactant.coords[3:] += np.array([0.0, 3.0, 0.0])
    intermediate = geometry("intermediate", 1.8, si_obr=1.8)
    intermediate.coords[3] = np.array([1.8, 0.98, 0.0])
    product = geometry("product", 1.7, si_obr=3.5)
    product.coords[3] = np.array([3.5, 0.98, 0.0])

    assert phase1.quick_irc_addition_reason(reactant, intermediate, 2) is None
    assert phase1.quick_irc_cleavage_reason(intermediate, product, 2) is None
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
