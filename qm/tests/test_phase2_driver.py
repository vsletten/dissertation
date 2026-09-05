"""Fast orchestration gates for the Phase-2 ladder driver (no DFT)."""

import ctypes
import json
import sys
import sysconfig
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from quarry.clusters import Cluster
from quarry.pipeline import DftSettings
from scripts import phase2_ladder as phase2

CHEAP = DftSettings(xc="hf", basis="sto-3g")


def test_preload_cutensor_fails_loudly_when_core_library_cannot_load(
    monkeypatch, tmp_path
):
    monkeypatch.setattr(sysconfig, "get_paths", lambda: {"purelib": str(tmp_path)})

    attempted = []

    def fail_load(path, *, mode):
        attempted.append((Path(path).name, mode))
        raise OSError("shared object not found")

    monkeypatch.setattr(ctypes, "CDLL", fail_load)

    with pytest.raises(RuntimeError, match="required cuTENSOR core library"):
        phase2.preload_cutensor()

    assert attempted[0] == ("libcutensor.so.2", ctypes.RTLD_GLOBAL)


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


def test_nearest_h_tracks_transferred_attacker_proton():
    cluster = Cluster(
        name="transferred-proton",
        symbols=["O", "H", "Si", "O", "H", "H"],
        coords=np.array(
            [
                [0.0, 0.0, 0.0],  # bridge oxygen
                [0.1, 0.0, 0.0],  # unrelated cluster hydrogen
                [3.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],  # appended attacker oxygen
                [0.98, 0.0, 0.0],  # delivered attacker H, no longer O-H bonded
                [2.0, 0.9, 0.0],  # other attacker H
            ]
        ),
    )

    assert phase2._nearest_h(cluster, br_index=0, ow_index=3) == 4


def test_attacker_h_indices_fail_fast_if_builder_order_changes():
    reordered = Cluster(
        "reordered",
        ["O", "Si", "O", "H", "O", "H"],
        np.zeros((6, 3)),
    )
    with pytest.raises(ValueError, match="ordering changed"):
        phase2._attacker_h_indices(reordered, 2)


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


def test_reverse_crest_scan_has_a_target_for_short_product_bond(monkeypatch):
    product = geometry("short-product", 1.8, m_obr=2.1)
    captured = {}

    def fake_scan(cluster, settings, *, distances_a, **kwargs):
        captured["distances"] = distances_a
        return [(distances_a[0], 0.0, cluster)]

    monkeypatch.setattr(phase2, "scan_to_maximum", fake_scan)
    monkeypatch.setattr(phase2, "scan_ts_guess", lambda scan: scan[0][2])

    assert (
        phase2.crest_from_product(
            product,
            CHEAP,
            m_index=1,
            br_index=0,
            ow_index=2,
            pin_a=1.9,
        )
        is product
    )
    assert captured["distances"] == [1.9]


def test_approach_seed_requires_matching_signature(tmp_path):
    seed = geometry("seed", 1.9)
    template = geometry("template", 3.2)
    signature = phase2.approach_seed_signature(
        template, m_index=1, br_index=0, ow_index=2, pin_a=1.9
    )
    path = tmp_path / "approach_seed.xyz"
    phase2.save_approach_seed(seed, path, signature)

    loaded = phase2.load_compatible_approach_seed(path, template, signature)
    assert loaded is not None
    assert np.array_equal(loaded.coords, seed.coords)

    incompatible = {**signature, "pin_a": 2.0}
    assert phase2.load_compatible_approach_seed(path, template, incompatible) is None


def test_reactant_complex_uses_checkpointed_hf_preoptimization(tmp_path, monkeypatch):
    guess = replace(geometry("guess", 3.2), frozen_indices=[0])
    preoptimized = replace(guess, coords=guess.coords + 0.01)
    production = replace(guess, coords=guess.coords + 0.02)
    production_settings = DftSettings(
        xc="b3lyp", basis="def2-svp", density_fit=True, use_gpu=True
    )
    calls = []

    def fake_optimize(cluster, settings):
        calls.append((cluster, settings))
        return preoptimized if settings.xc == "hf" else production

    monkeypatch.setattr(phase2, "optimize", fake_optimize)

    result = phase2.optimize_reactant_complex(tmp_path, guess, production_settings)

    assert [settings for _, settings in calls] == [
        DftSettings(xc="hf", basis="sto-3g"),
        production_settings,
    ]
    assert calls[1][0] is preoptimized
    assert result is production
    assert (tmp_path / "complex_preopt.xyz").exists()
    assert (tmp_path / "complex.xyz").exists()


def test_reactant_complex_preoptimization_failure_never_promotes_production(
    tmp_path, monkeypatch
):
    guess = geometry("guess", 3.2)
    calls = []

    def fail_preopt(cluster, settings):
        calls.append(settings)
        raise RuntimeError("preoptimization failed")

    monkeypatch.setattr(phase2, "optimize", fail_preopt)

    with pytest.raises(RuntimeError, match="preoptimization failed"):
        phase2.optimize_reactant_complex(tmp_path, guess, CHEAP)

    assert calls == [DftSettings(xc="hf", basis="sto-3g")]
    assert not (tmp_path / "complex_preopt.xyz").exists()
    assert not (tmp_path / "complex.xyz").exists()


def test_reactant_complex_resumes_preoptimization_before_production(
    tmp_path, monkeypatch
):
    guess = replace(geometry("guess", 3.2), frozen_indices=[0])
    preoptimized = replace(guess, coords=guess.coords + 0.01)
    phase2.save_xyz(preoptimized, tmp_path / "complex_preopt.xyz")
    production_settings = DftSettings(
        xc="b3lyp", basis="def2-svp", density_fit=True, use_gpu=True
    )
    calls = []

    def fake_optimize(cluster, settings):
        calls.append((cluster, settings))
        return replace(cluster, coords=cluster.coords + 0.01)

    monkeypatch.setattr(phase2, "optimize", fake_optimize)

    phase2.optimize_reactant_complex(tmp_path, guess, production_settings)

    assert len(calls) == 1
    resumed, settings = calls[0]
    assert settings == production_settings
    assert resumed.charge == guess.charge
    assert resumed.spin == guess.spin
    assert resumed.frozen_indices == guess.frozen_indices
    assert resumed.coords == pytest.approx(preoptimized.coords, abs=1e-9)


def test_existing_reactant_complex_skips_both_optimization_rungs(tmp_path, monkeypatch):
    guess = geometry("guess", 3.2)
    completed = replace(guess, coords=guess.coords + 0.02)
    phase2.save_xyz(completed, tmp_path / "complex.xyz")
    monkeypatch.setattr(
        phase2,
        "optimize",
        lambda *_args, **_kwargs: pytest.fail("completed complex was recomputed"),
    )

    resumed = phase2.optimize_reactant_complex(tmp_path, guess, CHEAP)

    assert resumed.coords == pytest.approx(completed.coords, abs=1e-9)
    assert not (tmp_path / "complex_preopt.xyz").exists()


def test_resume_persists_proton_route_before_reentering_proton_stage(
    monkeypatch, tmp_path
):
    cluster = geometry("cluster", 3.0)
    cc = SimpleNamespace(
        cluster=cluster,
        attacked_index=1,
        bridge_index=0,
        n_intact=4,
        metal_shells=2,
        termination_log=[],
        metadata=lambda: {"charge": 0},
    )
    run_dir = tmp_path / "phase2" / "oss-neutral-n4-s2-b3lyp-def2-svp"
    run_dir.mkdir(parents=True)
    complex_guess = geometry("complex", 3.0)
    phase2.save_approach_seed(
        geometry("approach-seed", 1.9),
        run_dir / "approach_seed.xyz",
        phase2.approach_seed_signature(
            complex_guess, m_index=1, br_index=0, ow_index=2, pin_a=1.9
        ),
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "phase2_ladder.py",
            "--family",
            "oss",
            "--run-root",
            str(tmp_path),
        ],
    )
    monkeypatch.setattr(phase2, "from_deck_cell", lambda *args, **kwargs: cc)
    monkeypatch.setattr(
        phase2, "attack_complex", lambda cc, attacker: (complex_guess, 2)
    )
    monkeypatch.setattr(phase2, "optimize", lambda cluster, settings: cluster)
    monkeypatch.setattr(phase2, "trim_gpu_pool", lambda: None)

    def stop_after_route_checkpoint(*args, **kwargs):
        assert (run_dir / "ts_guess.route").read_text() == "proton-neb"
        raise RuntimeError("route checkpoint observed")

    monkeypatch.setattr(phase2, "proton_neb_guess", stop_after_route_checkpoint)

    with pytest.raises(RuntimeError, match="route checkpoint observed"):
        phase2.main()
