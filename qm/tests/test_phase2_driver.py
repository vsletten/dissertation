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
    preoptimized = replace(guess, coords=guess.coords.copy())
    preoptimized.coords[1:] += 0.01
    production = replace(guess, coords=guess.coords.copy())
    production.coords[1:] += 0.02
    production_settings = DftSettings(
        xc="b3lyp", basis="def2-svp", density_fit=True, use_gpu=True
    )
    bounded_calls = []
    production_calls = []

    def fake_bounded(cluster, settings, *, max_steps):
        bounded_calls.append((cluster, settings, max_steps))
        return SimpleNamespace(cluster=preoptimized, converged=False)

    def fake_optimize(cluster, settings):
        production_calls.append((cluster, settings))
        return production

    monkeypatch.setattr(phase2, "optimize_bounded", fake_bounded)
    monkeypatch.setattr(phase2, "gradient", lambda *_args: np.full((5, 3), 0.01))
    monkeypatch.setattr(phase2, "optimize", fake_optimize)

    result = phase2.optimize_reactant_complex(tmp_path, guess, production_settings)

    assert bounded_calls == [(guess, DftSettings(xc="hf", basis="sto-3g"), 100)]
    assert len(production_calls) == 1
    production_seed, observed_settings = production_calls[0]
    assert observed_settings == production_settings
    assert production_seed.coords == pytest.approx(preoptimized.coords, abs=1e-9)
    assert production_seed.frozen_indices == guess.frozen_indices
    assert result.coords == pytest.approx(production.coords, abs=1e-8)
    assert (tmp_path / "complex_preopt.xyz").exists()
    assert (tmp_path / "complex.xyz").exists()
    assert (tmp_path / "complex.json").exists()
    receipt = json.loads((tmp_path / "complex_preopt.json").read_text())
    assert receipt["optimizer_converged"] is False
    assert receipt["signature"]["convergence_is_advisory"] is True
    assert receipt["signature"]["max_steps"] == 100
    assert receipt["production_qualification"]["status"] == "passed"
    assert receipt["geometry_gate"]["oxygen_proton_owners"] == ["H3:O2", "H4:O2"]
    assert receipt["endpoint_geometry_hash"] == phase2.geometry_hash(
        preoptimized.to_xyz()
    )


def test_reactant_complex_preoptimization_failure_never_promotes_production(
    tmp_path, monkeypatch
):
    guess = replace(geometry("guess", 3.2), frozen_indices=[0])
    preoptimized = replace(guess, coords=guess.coords.copy())
    monkeypatch.setattr(
        phase2,
        "optimize_bounded",
        lambda *_args, **_kwargs: SimpleNamespace(
            cluster=preoptimized, converged=False
        ),
    )
    monkeypatch.setattr(
        phase2,
        "gradient",
        lambda *_args: (_ for _ in ()).throw(RuntimeError("qualification failed")),
    )
    monkeypatch.setattr(
        phase2,
        "optimize",
        lambda *_args: pytest.fail("production optimization was reached"),
    )

    with pytest.raises(RuntimeError, match="qualification failed"):
        phase2.optimize_reactant_complex(tmp_path, guess, CHEAP)

    assert (tmp_path / "complex_preopt.xyz").exists()
    assert not (tmp_path / "complex.xyz").exists()
    receipt = json.loads((tmp_path / "complex_preopt.json").read_text())
    assert receipt["production_qualification"]["status"] == "failed"
    assert receipt["production_qualification"]["error_type"] == "RuntimeError"


def test_reactant_complex_resumes_preoptimization_before_production(
    tmp_path, monkeypatch
):
    guess = replace(geometry("guess", 3.2), frozen_indices=[0])
    preoptimized = replace(guess, coords=guess.coords.copy())
    preoptimized.coords[1:] += 0.01
    production_settings = DftSettings(
        xc="b3lyp", basis="def2-svp", density_fit=True, use_gpu=True
    )
    monkeypatch.setattr(
        phase2,
        "optimize_bounded",
        lambda *_args, **_kwargs: SimpleNamespace(
            cluster=preoptimized, converged=False
        ),
    )
    monkeypatch.setattr(phase2, "gradient", lambda *_args: np.full((5, 3), 0.01))
    monkeypatch.setattr(phase2, "optimize", lambda cluster, _settings: cluster)
    phase2.optimize_reactant_complex(tmp_path, guess, production_settings)
    (tmp_path / "complex.xyz").unlink()
    calls = []

    def fake_optimize(cluster, settings):
        calls.append((cluster, settings))
        coords = cluster.coords.copy()
        coords[1:] += 0.01
        return replace(cluster, coords=coords)

    monkeypatch.setattr(
        phase2,
        "optimize_bounded",
        lambda *_args, **_kwargs: pytest.fail("qualified preopt was recomputed"),
    )
    monkeypatch.setattr(
        phase2,
        "gradient",
        lambda *_args: pytest.fail("production qualification was recomputed"),
    )
    monkeypatch.setattr(phase2, "optimize", fake_optimize)

    phase2.optimize_reactant_complex(tmp_path, guess, production_settings)

    assert len(calls) == 1
    resumed, settings = calls[0]
    assert settings == production_settings
    assert resumed.charge == guess.charge
    assert resumed.spin == guess.spin
    assert resumed.frozen_indices == guess.frozen_indices
    assert resumed.coords == pytest.approx(preoptimized.coords, abs=1e-9)


def test_advisory_preoptimization_rejects_frozen_shell_drift(tmp_path, monkeypatch):
    guess = replace(geometry("guess", 3.2), frozen_indices=[0])
    drifted = replace(guess, coords=guess.coords.copy())
    drifted.coords[0, 0] += 0.0201
    monkeypatch.setattr(
        phase2,
        "optimize_bounded",
        lambda *_args, **_kwargs: SimpleNamespace(cluster=drifted, converged=False),
    )
    monkeypatch.setattr(
        phase2,
        "gradient",
        lambda *_args: pytest.fail("raw frozen-shell guard was bypassed"),
    )
    monkeypatch.setattr(
        phase2,
        "optimize",
        lambda *_args: pytest.fail("raw frozen-shell guard was bypassed"),
    )

    with pytest.raises(RuntimeError, match="exceeded its raw frozen-shell bound"):
        phase2.optimize_reactant_complex(tmp_path, guess, CHEAP)

    assert not (tmp_path / "complex_preopt.xyz").exists()


def test_advisory_preoptimization_projects_small_frozen_drift(tmp_path, monkeypatch):
    guess = replace(geometry("guess", 3.2), frozen_indices=[0])
    drifted = replace(guess, coords=guess.coords.copy())
    drifted.coords[0, 0] += 0.0115
    monkeypatch.setattr(
        phase2,
        "optimize_bounded",
        lambda *_args, **_kwargs: SimpleNamespace(cluster=drifted, converged=False),
    )
    monkeypatch.setattr(phase2, "gradient", lambda *_args: np.full((5, 3), 0.01))
    monkeypatch.setattr(phase2, "optimize", lambda cluster, _settings: cluster)

    endpoint = phase2.optimize_reactant_complex(tmp_path, guess, CHEAP)

    assert endpoint.coords[0] == pytest.approx(guess.coords[0], abs=1e-10)
    receipt = json.loads((tmp_path / "complex_preopt.json").read_text())
    assert receipt["unprojected_maximum_frozen_coordinate_drift_a"] == pytest.approx(
        0.0115
    )
    assert receipt["geometry_gate"]["maximum_frozen_coordinate_drift_a"] == 0.0


def test_advisory_preoptimization_rejects_changed_atom_order(tmp_path, monkeypatch):
    guess = replace(geometry("guess", 3.2), frozen_indices=[0])
    reordered = replace(guess, symbols=["Si", "O", "O", "H", "H"])
    monkeypatch.setattr(
        phase2,
        "optimize_bounded",
        lambda *_args, **_kwargs: SimpleNamespace(cluster=reordered, converged=False),
    )

    with pytest.raises(RuntimeError, match="changed atom identity/order"):
        phase2.optimize_reactant_complex(tmp_path, guess, CHEAP)


def test_advisory_preoptimization_rejects_changed_proton_owner(tmp_path, monkeypatch):
    guess = Cluster(
        "two-waters",
        ["O", "O", "H", "H"],
        np.array(
            [
                [0.0, 0.0, 0.0],
                [3.0, 0.0, 0.0],
                [0.96, 0.0, 0.0],
                [3.96, 0.0, 0.0],
            ]
        ),
    )
    transferred = replace(guess, coords=guess.coords.copy())
    transferred.coords[2, 0] = 2.04
    monkeypatch.setattr(
        phase2,
        "optimize_bounded",
        lambda *_args, **_kwargs: SimpleNamespace(cluster=transferred, converged=False),
    )
    monkeypatch.setattr(
        phase2,
        "gradient",
        lambda *_args: pytest.fail("microstate-changing seed reached production"),
    )
    monkeypatch.setattr(
        phase2,
        "optimize",
        lambda *_args: pytest.fail("microstate-changing seed reached production"),
    )

    with pytest.raises(RuntimeError, match="changed the reactant proton microstate"):
        phase2.optimize_reactant_complex(tmp_path, guess, CHEAP)

    assert not (tmp_path / "complex_preopt.xyz").exists()


@pytest.mark.parametrize(
    ("oxygen_x", "hydrogen_x", "message"),
    [
        ([0.0], 1.50, "is unassigned"),
        ([0.0, 2.0], 1.00, "has ambiguous owners"),
    ],
)
def test_reactant_proton_owner_must_be_bonded_and_unambiguous(
    oxygen_x, hydrogen_x, message
):
    cluster = Cluster(
        "invalid-owner",
        [*("O" for _ in oxygen_x), "H"],
        np.array([[x, 0.0, 0.0] for x in [*oxygen_x, hydrogen_x]]),
    )

    with pytest.raises(RuntimeError, match=message):
        phase2.oxygen_proton_owners(cluster)


def test_advisory_preoptimization_resume_is_bound_to_production_settings(
    tmp_path, monkeypatch
):
    guess = replace(geometry("guess", 3.2), frozen_indices=[0])
    preoptimized = replace(guess, coords=guess.coords.copy())
    monkeypatch.setattr(
        phase2,
        "optimize_bounded",
        lambda *_args, **_kwargs: SimpleNamespace(
            cluster=preoptimized, converged=False
        ),
    )
    monkeypatch.setattr(phase2, "gradient", lambda *_args: np.full((5, 3), 0.01))
    monkeypatch.setattr(phase2, "optimize", lambda cluster, _settings: cluster)
    phase2.optimize_reactant_complex(tmp_path, guess, CHEAP)
    (tmp_path / "complex.xyz").unlink()

    with pytest.raises(RuntimeError, match="receipt signature mismatch"):
        phase2.optimize_reactant_complex(
            tmp_path, guess, replace(CHEAP, density_fit=True)
        )


@pytest.mark.parametrize("receipt_state", [None, "pending"])
def test_incomplete_advisory_checkpoint_is_preserved_and_recomputed(
    tmp_path, monkeypatch, receipt_state
):
    guess = replace(geometry("guess", 3.2), frozen_indices=[0])
    phase2.save_xyz(guess, tmp_path / "complex_preopt.xyz")
    if receipt_state is not None:
        (tmp_path / "complex_preopt.json").write_text(
            json.dumps({"production_qualification": {"status": receipt_state}})
        )
    calls = []

    def fake_bounded(*_args, **_kwargs):
        calls.append("bounded")
        return SimpleNamespace(cluster=guess, converged=False)

    monkeypatch.setattr(phase2, "optimize_bounded", fake_bounded)
    monkeypatch.setattr(phase2, "gradient", lambda *_args: np.full((5, 3), 0.01))
    monkeypatch.setattr(phase2, "optimize", lambda cluster, _settings: cluster)

    phase2.optimize_reactant_complex(tmp_path, guess, CHEAP)

    assert calls == ["bounded"]
    assert list(tmp_path.glob("complex_preopt.incomplete-*.xyz"))
    assert (tmp_path / "complex_preopt.json").exists()


def test_existing_reactant_complex_without_receipt_is_rejected(tmp_path):
    guess = geometry("guess", 3.2)
    completed = replace(guess, coords=guess.coords + 0.02)
    phase2.save_xyz(completed, tmp_path / "complex.xyz")

    with pytest.raises(RuntimeError, match="refusing unchecked complex.xyz"):
        phase2.optimize_reactant_complex(tmp_path, guess, CHEAP)


def test_verified_reactant_complex_skips_both_optimization_rungs(tmp_path, monkeypatch):
    guess = geometry("guess", 3.2)
    monkeypatch.setattr(
        phase2,
        "optimize_bounded",
        lambda *_args, **_kwargs: SimpleNamespace(cluster=guess, converged=True),
    )
    monkeypatch.setattr(phase2, "gradient", lambda *_args: np.full((5, 3), 0.01))
    monkeypatch.setattr(phase2, "optimize", lambda cluster, _settings: cluster)
    completed = phase2.optimize_reactant_complex(tmp_path, guess, CHEAP)
    monkeypatch.setattr(
        phase2,
        "optimize",
        lambda *_args, **_kwargs: pytest.fail("completed complex was recomputed"),
    )

    resumed = phase2.optimize_reactant_complex(tmp_path, guess, CHEAP)

    assert resumed.coords == pytest.approx(completed.coords, abs=1e-9)


def test_production_endpoint_rejects_changed_proton_owner(tmp_path, monkeypatch):
    guess = Cluster(
        "two-waters",
        ["O", "O", "H", "H"],
        np.array(
            [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.96, 0.0, 0.0], [3.96, 0.0, 0.0]]
        ),
    )
    transferred = replace(guess, coords=guess.coords.copy())
    transferred.coords[2, 0] = 2.04
    monkeypatch.setattr(
        phase2,
        "optimize_bounded",
        lambda *_args, **_kwargs: SimpleNamespace(cluster=guess, converged=True),
    )
    monkeypatch.setattr(phase2, "gradient", lambda *_args: np.full((4, 3), 0.01))
    monkeypatch.setattr(phase2, "optimize", lambda *_args: transferred)

    with pytest.raises(RuntimeError, match="production optimization changed"):
        phase2.optimize_reactant_complex(tmp_path, guess, CHEAP)

    assert not (tmp_path / "complex.xyz").exists()


def test_load_xyz_rejects_changed_atom_order(tmp_path):
    template = geometry("template", 3.2)
    path = tmp_path / "bad.xyz"
    lines = template.to_xyz().splitlines()
    lines[2], lines[3] = lines[3], lines[2]
    path.write_text("\n".join(lines))

    with pytest.raises(ValueError, match="atom-order or symbol mismatch"):
        phase2.load_xyz(path, template)


def test_reactant_minimum_rejects_significant_imaginary_mode():
    assert phase2.reactant_minimum_reason(np.array([12.0, 35.0])) is not None
    assert phase2.reactant_minimum_reason(np.array([12.0, 29.0])) is None


def test_quick_irc_must_span_intact_bridge_and_hydrolyzed_product():
    intact = geometry("intact", 3.0, m_obr=1.6)
    product = geometry("product", 1.9, m_obr=3.0)
    product.coords[3] = np.array([3.98, 0.0, 0.0])
    assert (
        phase2.quick_irc_acceptance_reason(
            intact, product, m_index=1, br_index=0, ow_index=2
        )
        is None
    )
    assert (
        phase2.quick_irc_acceptance_reason(
            intact, intact, m_index=1, br_index=0, ow_index=2
        )
        is not None
    )


def test_new_attempt_quarantines_stale_canonical_outputs(tmp_path):
    for name in ("results.json", "store.sqlite", "store.sqlite-wal"):
        (tmp_path / name).write_text(name)

    quarantine = phase2.quarantine_canonical_outputs(tmp_path)

    assert quarantine is not None
    assert {path.name for path in quarantine.iterdir()} == {
        "results.json",
        "store.sqlite",
        "store.sqlite-wal",
    }
    assert not (tmp_path / "results.json").exists()


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
    monkeypatch.setattr(
        phase2,
        "optimize_reactant_complex",
        lambda _run_dir, _guess, _settings: complex_guess,
    )
    monkeypatch.setattr(
        phase2,
        "frequencies",
        lambda *_args, **_kwargs: SimpleNamespace(imaginary_cm=np.array([])),
    )
    monkeypatch.setattr(phase2, "optimize", lambda cluster, settings: cluster)
    monkeypatch.setattr(phase2, "trim_gpu_pool", lambda: None)

    def stop_after_route_checkpoint(*args, **kwargs):
        assert (run_dir / "ts_guess.route").read_text() == "proton-neb"
        raise RuntimeError("route checkpoint observed")

    monkeypatch.setattr(phase2, "proton_neb_guess", stop_after_route_checkpoint)

    with pytest.raises(RuntimeError, match="route checkpoint observed"):
        phase2.main()


def test_load_xyz_rejects_extra_atom_fields(tmp_path):
    template = geometry("template", 3.2)
    path = tmp_path / "bad.xyz"
    lines = template.to_xyz().splitlines()
    lines[2] = f"{lines[2]} extra"
    path.write_text("\n".join(lines))

    with pytest.raises(ValueError, match="atom records are malformed"):
        phase2.load_xyz(path, template)


def test_load_xyz_rejects_trailing_atom_records(tmp_path):
    template = geometry("template", 3.2)
    path = tmp_path / "bad.xyz"
    lines = template.to_xyz().splitlines()
    lines.append("H 0.00000000 0.00000000 0.00000000")
    path.write_text("\n".join(lines))

    with pytest.raises(ValueError, match="line-count mismatch"):
        phase2.load_xyz(path, template)


def test_reactant_recovery_fixes_owner_bonds_then_releases_once(tmp_path, monkeypatch):
    guess = replace(geometry("guess", 3.2), frozen_indices=[0])
    conditioned = replace(guess, coords=guess.coords.copy())
    attempt_one = replace(guess, coords=guess.coords.copy())
    attempt_one.coords[1, 1] = 0.01
    attempt_two = replace(guess, coords=guess.coords.copy())
    attempt_two.coords[1, 1] = 0.02
    results = iter(
        [
            SimpleNamespace(cluster=conditioned, converged=False),
            SimpleNamespace(cluster=attempt_one, converged=False),
            SimpleNamespace(cluster=attempt_two, converged=True),
        ]
    )
    calls = []

    def fake_bounded(cluster, settings, **kwargs):
        calls.append((cluster, settings, kwargs))
        return next(results)

    monkeypatch.setattr(phase2, "optimize_bounded", fake_bounded)
    monkeypatch.setattr(phase2, "gradient", lambda *_args: np.zeros((5, 3)))
    monkeypatch.setattr(
        phase2,
        "frequencies",
        lambda *_args: SimpleNamespace(
            imaginary_cm=np.array([]), electronic_hartree=-10.0
        ),
    )

    recovered = phase2.recover_reactant_minimum(tmp_path, guess, CHEAP)

    expected_constraints = phase2.owner_bond_constraints(guess)
    assert calls[0][2] == {
        "max_steps": phase2.REACTANT_CONDITIONING_MAX_STEPS,
        "fixed_distances": expected_constraints,
    }
    assert calls[0][1] == DftSettings(xc="hf", basis="sto-3g")
    assert calls[1][2] == {"max_steps": phase2.REACTANT_PRODUCTION_MAX_STEPS}
    assert calls[2][2] == {"max_steps": phase2.REACTANT_PRODUCTION_MAX_STEPS}
    assert calls[2][0].coords == pytest.approx(attempt_one.coords, abs=1e-8)
    assert recovered.coords == pytest.approx(attempt_two.coords, abs=1e-8)
    assert (tmp_path / "complex_production_attempt_1.xyz").exists()
    assert (tmp_path / "complex_production_attempt_2.xyz").exists()
    assert (tmp_path / "complex.xyz").exists()
    receipt = json.loads((tmp_path / "complex.json").read_text())
    assert receipt["accepted_attempt"] == 2
    assert receipt["minimum_gate"]["status"] == "passed"
    terminal = json.loads((tmp_path / "production-terminal.json").read_text())
    assert terminal["status"] == "success"


def test_reactant_recovery_exhaustion_preserves_both_endpoints(tmp_path, monkeypatch):
    guess = replace(geometry("guess", 3.2), frozen_indices=[0])
    calls = []

    def fake_bounded(cluster, settings, **kwargs):
        calls.append((cluster, settings, kwargs))
        moved = replace(cluster, coords=cluster.coords.copy())
        moved.coords[1, 1] += 0.01
        return SimpleNamespace(cluster=moved, converged=False)

    monkeypatch.setattr(phase2, "optimize_bounded", fake_bounded)
    monkeypatch.setattr(
        phase2,
        "gradient",
        lambda *_args: pytest.fail("gradient ran after production exhaustion"),
    )
    monkeypatch.setattr(
        phase2,
        "frequencies",
        lambda *_args: pytest.fail("Hessian ran after production exhaustion"),
    )

    with pytest.raises(RuntimeError, match="2 x 100 steps"):
        phase2.recover_reactant_minimum(tmp_path, guess, CHEAP)

    assert len(calls) == 3  # one conditioning plus exactly two production optimizers
    for attempt in (1, 2):
        endpoint = tmp_path / f"complex_production_attempt_{attempt}.xyz"
        receipt = tmp_path / f"complex_production_attempt_{attempt}.json"
        assert endpoint.exists() and receipt.exists()
        assert json.loads(receipt.read_text())["converged"] is False
    assert not (tmp_path / "complex.xyz").exists()
    terminal = json.loads((tmp_path / "production-terminal.json").read_text())
    assert terminal["stage"] == "production-step-exhaustion"
    assert terminal["attempt"] == 2


@pytest.mark.parametrize(
    ("gradient_value", "imaginary", "message"),
    [
        (0.001, np.array([]), "gradient exceeds"),
        (0.0, np.array([45.0]), "imaginary mode"),
    ],
)
def test_reactant_recovery_never_promotes_failed_minimum(
    tmp_path, monkeypatch, gradient_value, imaginary, message
):
    guess = replace(geometry("guess", 3.2), frozen_indices=[0])
    results = iter(
        [
            SimpleNamespace(cluster=guess, converged=False),
            SimpleNamespace(cluster=guess, converged=True),
        ]
    )
    monkeypatch.setattr(
        phase2, "optimize_bounded", lambda *_args, **_kwargs: next(results)
    )
    monkeypatch.setattr(
        phase2, "gradient", lambda *_args: np.full((5, 3), gradient_value)
    )
    monkeypatch.setattr(
        phase2,
        "frequencies",
        lambda *_args: SimpleNamespace(
            imaginary_cm=imaginary, electronic_hartree=-10.0
        ),
    )

    with pytest.raises(RuntimeError, match=message):
        phase2.recover_reactant_minimum(tmp_path, guess, CHEAP)

    assert not (tmp_path / "complex.xyz").exists()
    terminal = json.loads((tmp_path / "production-terminal.json").read_text())
    assert terminal["stage"] == "independent-reactant-minimum-gate"


def test_reactant_recovery_refuses_ambiguous_owner_before_calculator(
    tmp_path, monkeypatch
):
    guess = Cluster(
        "ambiguous",
        ["O", "O", "H"],
        np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
    )
    monkeypatch.setattr(
        phase2,
        "optimize_bounded",
        lambda *_args, **_kwargs: pytest.fail(
            "calculator reached with ambiguous owner"
        ),
    )

    with pytest.raises(RuntimeError, match="ambiguous owners"):
        phase2.recover_reactant_minimum(tmp_path, guess, CHEAP)


def test_reactant_recovery_resume_rejects_endpoint_hash_drift(tmp_path, monkeypatch):
    guess = replace(geometry("guess", 3.2), frozen_indices=[0])
    results = iter(
        [
            SimpleNamespace(cluster=guess, converged=False),
            SimpleNamespace(cluster=guess, converged=False),
            SimpleNamespace(cluster=guess, converged=False),
        ]
    )
    monkeypatch.setattr(
        phase2, "optimize_bounded", lambda *_args, **_kwargs: next(results)
    )
    with pytest.raises(RuntimeError, match="2 x 100 steps"):
        phase2.recover_reactant_minimum(tmp_path, guess, CHEAP)
    receipt_path = tmp_path / "complex_production_attempt_1.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["endpoint_geometry_hash"] = "0" * 64
    receipt_path.write_text(json.dumps(receipt))

    with pytest.raises(RuntimeError, match="attempt 1 endpoint hash mismatch"):
        phase2.recover_reactant_minimum(tmp_path, guess, CHEAP)


def test_reactant_recovery_resume_rejects_production_settings_drift(
    tmp_path, monkeypatch
):
    guess = replace(geometry("guess", 3.2), frozen_indices=[0])
    results = iter(
        [
            SimpleNamespace(cluster=guess, converged=False),
            SimpleNamespace(cluster=guess, converged=False),
            SimpleNamespace(cluster=guess, converged=False),
        ]
    )
    monkeypatch.setattr(
        phase2, "optimize_bounded", lambda *_args, **_kwargs: next(results)
    )
    with pytest.raises(RuntimeError, match="2 x 100 steps"):
        phase2.recover_reactant_minimum(tmp_path, guess, CHEAP)
    receipt_path = tmp_path / "complex_production_attempt_1.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["signature"]["production"]["settings"]["basis"] = "def2-svp"
    receipt_path.write_text(json.dumps(receipt))

    with pytest.raises(RuntimeError, match="attempt 1 receipt signature mismatch"):
        phase2.recover_reactant_minimum(tmp_path, guess, CHEAP)
