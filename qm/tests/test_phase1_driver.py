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

    def fake_neb(reactant, product, settings):
        neb_inputs.append((reactant, product))
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
