"""Fast orchestration gates for the Phase-2 ladder driver (no DFT)."""

import json
import sys

import numpy as np
import pytest

from quarry.clusters import Cluster
from scripts import phase2_ladder as phase2


def geometry(name: str, m_ow: float, m_obr: float = 1.6) -> Cluster:
    return Cluster(
        name=name,
        symbols=["O", "Si", "O", "H", "H"],
        coords=np.array(
            [
                [m_obr, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [m_ow, 0.0, 0.0],
                [m_ow, 0.9, 0.0],
                [m_ow, -0.9, 0.0],
            ]
        ),
    )


def test_channel_escape_guard_uses_metal_limit():
    guess = geometry("guess", 2.2)
    # Outside the Si channel (limit 2.6) but inside the Al channel (2.8).
    drifted = geometry("drift", 2.7)
    assert phase2.channel_escape_reason(guess, drifted, 1, 2, 2.6) is not None
    assert phase2.channel_escape_reason(guess, drifted, 1, 2, 2.8) is None
    # Relative escape trips regardless of the absolute limit.
    rel = phase2.channel_escape_reason(
        geometry("g", 1.8), geometry("d", 2.4), 1, 2, 2.8
    )
    assert rel is not None and "escaped" in rel


def test_nearest_h_does_not_depend_on_water_atom_order():
    cluster = Cluster(
        name="reordered-water",
        symbols=["O", "H", "Si", "H", "O", "C", "H"],
        coords=np.array(
            [
                [0.0, 0.0, 0.0],  # bridge oxygen
                [0.1, 0.0, 0.0],  # unrelated hydrogen, far from water oxygen
                [3.0, 0.0, 0.0],
                [1.1, 0.0, 0.0],  # delivered water hydrogen
                [2.0, 0.0, 0.0],  # water oxygen
                [0.2, 0.0, 0.0],  # non-H at the old positional candidate
                [2.0, 0.9, 0.0],  # other water hydrogen
            ]
        ),
    )

    assert phase2._nearest_h(cluster, br_index=0, ow_index=4) == 3


@pytest.mark.parametrize("family", ["oss", "osa", "oaa"])
def test_dry_run_writes_geometry_and_metadata(family, tmp_path, monkeypatch):
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "phase2_ladder.py",
            "--family",
            family,
            "--dry-run",
            "--run-root",
            str(tmp_path),
        ],
    )
    assert phase2.main() == 0
    candidates = sorted((tmp_path / "phase2").glob(f"{family}-neutral-*"))
    assert candidates, "dry run produced no run dir"
    run_dir = candidates[-1]
    meta = json.loads((run_dir / "metadata.json").read_text())
    assert meta["state"] == "neutral"
    assert meta["charge"] == 0
    assert meta["n_frozen"] > 0
    assert (run_dir / "cluster.xyz").exists()
    assert (run_dir / "complex_guess.xyz").exists()


def test_approach_parameters_cover_both_metals():
    assert phase2.APPROACH["Si"]["pin"] == pytest.approx(1.90)
    assert phase2.APPROACH["Al"]["pin"] == pytest.approx(2.00)
    for p in phase2.APPROACH.values():
        assert p["limit"] > p["pin"]
        assert min(p["distances"]) >= p["pin"]
