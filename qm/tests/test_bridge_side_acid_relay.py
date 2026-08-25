from dataclasses import asdict, replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from quarry.clusters import (
    BRIDGE_SIDE_ACID_FAMILY,
    BRIDGE_SIDE_ACID_WATER_COUNT,
    BridgeSideHydroniumEndpoints,
    Cluster,
    aluminosilicate_dimer,
    bridge_side_hydronium_endpoints,
    disilicate,
)
from quarry.pipeline import DftSettings, FrequencyResult
from scripts import bridge_side_acid_relay as relay
from scripts import concerted_acid_relay as shared
from scripts import phase1_xiao_lasaga as phase1
from scripts import separated_acid_relay as a1f

CHEAP = DftSettings(xc="hf", basis="sto-3g")


def _frequency(cluster: Cluster, mode: np.ndarray, *, n_imaginary: int = 1):
    return FrequencyResult(
        frequencies_cm=np.array([1000.0]),
        imaginary_cm=np.full(n_imaginary, 250.0),
        electronic_hartree=-100.0,
        molar_mass_kg=0.1,
        rotational_temperatures_k=(1.0, 2.0, 3.0),
        linear=False,
        imaginary_mode=mode if n_imaginary else None,
        imaginary_modes=(
            np.repeat(mode[None, ...], n_imaginary, axis=0) if n_imaginary else None
        ),
    )


@pytest.mark.parametrize("model", [disilicate, aluminosilicate_dimer])
def test_bridge_side_endpoint_seeds_have_exact_roles_basins_and_no_collisions(model):
    ends = bridge_side_hydronium_endpoints(model())

    assert ends.hydronium_oxygen_index != ends.ow_index
    assert ends.solvent_oxygen_indices == (
        ends.outer_solvator_oxygen_index,
        ends.ow_index,
        ends.hydronium_oxygen_index,
        ends.spectator_oxygen_index,
    )
    assert len(ends.solvent_oxygen_indices) == BRIDGE_SIDE_ACID_WATER_COUNT
    assert len(ends.solvent_h_indices) == 2 * BRIDGE_SIDE_ACID_WATER_COUNT + 1
    assert ends.relay_h_indices == ends.hydronium_h_indices
    assert ends.reactant.symbols == ends.product.symbols
    assert ends.reactant.frozen_indices == ends.product.frozen_indices == []
    for endpoint, cluster in (("reactant", ends.reactant), ("product", ends.product)):
        assert np.all(np.isfinite(cluster.coords))
        assert shared.minimum_pair_distance(cluster) >= relay.MIN_PAIR_DISTANCE_A
        assert (
            relay.endpoint_gate_reason(
                cluster,
                ends,
                n_water=BRIDGE_SIDE_ACID_WATER_COUNT,
                endpoint=endpoint,
            )
            is None
        )
        assert phase1.acid_basin_signature(
            cluster,
            ends.ow_index,
            ends.solvent_h_indices,
            solvent_oxygen_indices=ends.solvent_oxygen_indices,
        ) == relay.expected_basin(BRIDGE_SIDE_ACID_WATER_COUNT, endpoint)


def test_physical_h_identity_follows_migration_and_direct_bridge_transfer():
    ends = bridge_side_hydronium_endpoints(disilicate())
    all_h = tuple(
        index for index, symbol in enumerate(ends.reactant.symbols) if symbol == "H"
    )
    reactant_owners = dict(
        zip(
            all_h,
            phase1._acid_mobile_assignments(ends.reactant, all_h),
            strict=True,
        )
    )
    product_owners = dict(
        zip(
            all_h,
            phase1._acid_mobile_assignments(ends.product, all_h),
            strict=True,
        )
    )

    assert None not in reactant_owners.values()
    assert None not in product_owners.values()
    assert reactant_owners[ends.migrated_h_index] == ends.hydronium_oxygen_index
    assert product_owners[ends.migrated_h_index] == ends.hydronium_oxygen_index
    assert reactant_owners[ends.transferred_h_index] == ends.hydronium_oxygen_index
    assert product_owners[ends.transferred_h_index] == phase1.BR_INDEX
    attacker_h = (ends.ow_index + 1, ends.ow_index + 2)
    assert all(reactant_owners[index] == ends.ow_index for index in attacker_h)
    assert all(product_owners[index] == ends.ow_index for index in attacker_h)
    assert all(
        reactant_owners[index] == ends.hydronium_oxygen_index
        for index in ends.hydronium_h_indices
    )


@pytest.mark.parametrize("n_water", [1, 3, 5, 6])
def test_bridge_side_builder_refuses_every_water_count_except_four(n_water):
    with pytest.raises(ValueError, match="exactly 4 waters"):
        bridge_side_hydronium_endpoints(disilicate(), n_water=n_water)


def test_bridge_side_builder_refuses_prior_families():
    for family in (
        "bridge-donor-chain",
        "compact-cyclic-relay",
        "separated-donor-neutral-attacker",
        "other",
    ):
        with pytest.raises(ValueError, match="family must be"):
            bridge_side_hydronium_endpoints(disilicate(), family=family)


def test_endpoint_reaction_vector_passes_direct_transfer_mode_components():
    ends = bridge_side_hydronium_endpoints(disilicate())
    mode = ends.product.coords - ends.reactant.coords

    result = relay.coupled_mode_components(
        ends.reactant,
        _frequency(ends.reactant, mode),
        ends,
    )

    assert result["accepted"] is True
    assert result["components"]["si_obr_cleavage"] > 0.0
    assert result["components"]["si_ow_attack"] < 0.0
    assert result["components"]["hydronium_h_release"] > 0.0
    assert result["components"]["hydronium_to_obr"] < 0.0


def test_direct_transfer_mode_gate_rejects_a_water_wag():
    ends = bridge_side_hydronium_endpoints(disilicate())
    mode = np.zeros_like(ends.reactant.coords)
    mode[ends.ow_index + 1, 2] = 1.0

    result = relay.coupled_mode_components(
        ends.reactant,
        _frequency(ends.reactant, mode),
        ends,
    )

    assert result["accepted"] is False
    assert "direct hydronium transfer" in result["reason"]


def test_endpoint_gate_rejects_unassigned_and_double_owned_hydrogen():
    ends = bridge_side_hydronium_endpoints(disilicate())
    unassigned = ends.reactant.coords.copy()
    unassigned[ends.transferred_h_index] += np.array([20.0, 0.0, 0.0])
    reason = relay.endpoint_gate_reason(
        replace(ends.reactant, coords=unassigned),
        ends,
        n_water=BRIDGE_SIDE_ACID_WATER_COUNT,
        endpoint="reactant",
    )
    assert reason is not None

    ambiguous = ends.reactant.coords.copy()
    ambiguous[ends.transferred_h_index] = 0.5 * (
        ambiguous[ends.hydronium_oxygen_index] + ambiguous[phase1.BR_INDEX]
    )
    reason = relay.endpoint_gate_reason(
        replace(ends.reactant, coords=ambiguous),
        ends,
        n_water=BRIDGE_SIDE_ACID_WATER_COUNT,
        endpoint="reactant",
    )
    assert reason is not None


def test_full_irc_requires_exact_typed_endpoints():
    ends = bridge_side_hydronium_endpoints(disilicate())

    assert (
        relay.irc_channel_reason(
            ends.reactant,
            ends.product,
            ends,
            n_water=BRIDGE_SIDE_ACID_WATER_COUNT,
        )
        is None
    )
    reason = relay.irc_channel_reason(
        ends.reactant,
        ends.reactant,
        ends,
        n_water=BRIDGE_SIDE_ACID_WATER_COUNT,
    )
    assert reason is not None
    assert "full IRC endpoint basins" in reason


@pytest.mark.parametrize(
    ("mechanism_version", "gate_version"),
    [
        (shared.MECHANISM_VERSION, shared.GATE_VERSION),
        (a1f.MECHANISM_VERSION, a1f.GATE_VERSION),
    ],
)
def test_a1g_terminal_refuses_a1e_and_a1f_identity(
    tmp_path, mechanism_version, gate_version
):
    bounds = relay.Bounds()
    path_dir = tmp_path / "path"
    path_dir.mkdir()
    artifact = path_dir / "endpoint.xyz"
    artifact.write_text("proof\n")
    payload = {
        "schema_version": shared.SCHEMA_VERSION,
        "mechanism_version": mechanism_version,
        "gate_version": gate_version,
        "key": f"si:4w:{BRIDGE_SIDE_ACID_FAMILY}",
        "settings": asdict(CHEAP),
        "bounds": asdict(bounds),
        "source_reactant_sha256": relay.A1F_SOURCE_REACTANT_SHA256,
        "status": "rejected",
        "artifacts": {"endpoint.xyz": shared.sha256_path(artifact)},
    }
    shared.atomic_json(path_dir / "terminal.json", payload)

    assert not shared.terminal_record_is_reusable(
        path_dir,
        key=f"si:4w:{BRIDGE_SIDE_ACID_FAMILY}",
        settings=CHEAP,
        bounds=bounds,
        mechanism_version=relay.MECHANISM_VERSION,
        gate_version=relay.GATE_VERSION,
        identity_extra={"source_reactant_sha256": relay.A1F_SOURCE_REACTANT_SHA256},
    )


def test_source_reactant_is_hash_checked_and_keeps_a1g_identity(tmp_path):
    ends = bridge_side_hydronium_endpoints(disilicate())
    source = tmp_path / "a1f-reactant.xyz"
    shared.save_xyz(ends.reactant, source)
    digest = shared.sha256_path(source)

    loaded = relay.load_source_reactant(source, digest, ends)

    assert isinstance(loaded, BridgeSideHydroniumEndpoints)
    assert np.allclose(loaded.reactant.coords, ends.reactant.coords, atol=1.0e-8)
    with pytest.raises(ValueError, match="SHA-256"):
        relay.load_source_reactant(source, "0" * 64, ends)


def test_run_path_forwards_a1g_identity_builder_gates_and_source(monkeypatch, tmp_path):
    source_ends = bridge_side_hydronium_endpoints(disilicate())
    source = tmp_path / "a1f-reactant.xyz"
    shared.save_xyz(source_ends.reactant, source)
    digest = shared.sha256_path(source)
    captured = {}

    def fake_run(*args, **kwargs):
        captured.update(kwargs)
        return {"status": "rejected"}

    monkeypatch.setattr(shared, "run_path", fake_run)
    relay.run_path(
        tmp_path,
        CHEAP,
        relay.Bounds(),
        model="si",
        source_reactant=source,
        source_reactant_sha256=digest,
    )

    assert captured["n_water"] == BRIDGE_SIDE_ACID_WATER_COUNT
    assert captured["family"] == BRIDGE_SIDE_ACID_FAMILY
    assert captured["mechanism_version"] == relay.MECHANISM_VERSION
    assert captured["gate_version"] == relay.GATE_VERSION
    assert captured["endpoint_gate"] is relay.endpoint_gate_reason
    assert captured["mode_gate"] is relay.coupled_mode_components
    assert captured["irc_gate"] is relay.irc_channel_reason
    assert captured["identity_extra"] == {
        "source_reactant_sha256": digest,
        "parent_si_terminal_sha256": None,
    }
    built = captured["endpoint_builder"](disilicate())
    assert np.allclose(built.reactant.coords, source_ends.reactant.coords, atol=1.0e-8)


def test_run_path_fake_acceptance_reaches_persisted_neb_and_saddle(
    monkeypatch, tmp_path
):
    ends = bridge_side_hydronium_endpoints(disilicate())
    optimized = iter((ends.reactant, ends.product))
    mode = ends.product.coords - ends.reactant.coords
    frequency_calls = 0
    observed = {}

    def fake_optimize(cluster, settings, *, max_steps):
        return SimpleNamespace(cluster=next(optimized), converged=True)

    def fake_frequency(path, cluster, settings):
        nonlocal frequency_calls
        frequency_calls += 1
        if frequency_calls < 3:
            return _frequency(cluster, np.zeros_like(cluster.coords), n_imaginary=0)
        return _frequency(cluster, mode)

    def fake_neb(reactant, product, settings, **kwargs):
        checkpoint_dir = kwargs["checkpoint_dir"]
        checkpoint_dir.mkdir(parents=True)
        (checkpoint_dir / "pre-relaxed-images.npz").write_bytes(b"checkpoint")
        observed["checkpoint_dir"] = checkpoint_dir
        observed["pre_relax_steps"] = kwargs["pre_relax_steps"]
        return replace(reactant, coords=0.5 * (reactant.coords + product.coords))

    def fake_find_ts(crest, settings, **kwargs):
        observed["active_mode_norm"] = float(np.linalg.norm(kwargs["initial_mode"]))
        return crest

    monkeypatch.setattr(shared, "optimize_bounded", fake_optimize)
    monkeypatch.setattr(shared, "cached_frequency", fake_frequency)
    monkeypatch.setattr(shared, "neb_ts_guess", fake_neb)
    monkeypatch.setattr(shared, "find_ts", fake_find_ts)
    monkeypatch.setattr(
        shared,
        "full_irc",
        lambda saddle, settings, **kwargs: (ends.reactant, ends.product),
    )
    monkeypatch.setattr(
        shared.phase1,
        "thermo_result",
        lambda frequency, temperature: SimpleNamespace(
            gibbs=100.0 if frequency.n_imaginary else 0.0
        ),
    )

    record = relay.run_path(tmp_path, CHEAP, relay.Bounds(), model="si")

    assert record["status"] == "accepted"
    assert record["stage"] == "completed"
    assert record["barrier_kj_mol"] == 100.0
    assert observed["pre_relax_steps"] == relay.Bounds().neb_pre_steps
    assert observed["active_mode_norm"] > 0.0
    assert (observed["checkpoint_dir"] / "pre-relaxed-images.npz").is_file()


def _write_fake_terminal(run_dir: Path, model: str, status: str):
    path_dir = shared.path_directory(
        run_dir,
        model,
        BRIDGE_SIDE_ACID_WATER_COUNT,
        BRIDGE_SIDE_ACID_FAMILY,
    )
    path_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "key": shared.path_key(
            model,
            BRIDGE_SIDE_ACID_WATER_COUNT,
            BRIDGE_SIDE_ACID_FAMILY,
        ),
        "model": model,
        "water_count": BRIDGE_SIDE_ACID_WATER_COUNT,
        "family": BRIDGE_SIDE_ACID_FAMILY,
        "status": status,
        "barrier_kj_mol": 100.0 if model == "si" else 80.0,
    }
    shared.atomic_json(path_dir / "terminal.json", payload)
    return payload


def test_campaign_runs_one_si_path_and_never_unmatched_al(monkeypatch, tmp_path):
    calls = []

    def fake_run(run_dir, settings, bounds, *, model, **kwargs):
        calls.append(model)
        return _write_fake_terminal(run_dir, model, "rejected")

    monkeypatch.setattr(relay, "run_path", fake_run)
    manifest = relay.run_campaign(
        tmp_path,
        CHEAP,
        relay.Bounds(),
        log_path=str(tmp_path / "run.log"),
        source_reactant=None,
        source_reactant_sha256=None,
    )

    assert calls == ["si"]
    assert manifest["summary"]["verdict"] == "si-path-conclusive-rejection"
    assert manifest["summary"]["al_status"] is None


def test_campaign_triggers_only_atom_matched_al_after_accepted_si(
    monkeypatch, tmp_path
):
    calls = []

    def fake_run(run_dir, settings, bounds, *, model, **kwargs):
        calls.append(
            (
                model,
                kwargs.get("source_reactant"),
                kwargs.get("parent_si_terminal_sha256"),
            )
        )
        return _write_fake_terminal(run_dir, model, "accepted")

    monkeypatch.setattr(relay, "run_path", fake_run)
    manifest = relay.run_campaign(
        tmp_path,
        CHEAP,
        relay.Bounds(),
        log_path=str(tmp_path / "run.log"),
        source_reactant=tmp_path / "source.xyz",
        source_reactant_sha256="a" * 64,
    )

    si_terminal = (
        shared.path_directory(
            tmp_path, "si", BRIDGE_SIDE_ACID_WATER_COUNT, BRIDGE_SIDE_ACID_FAMILY
        )
        / "terminal.json"
    )
    assert calls == [
        ("si", tmp_path / "source.xyz", None),
        ("al", None, shared.sha256_path(si_terminal)),
    ]
    assert manifest["summary"]["verdict"] == "matched-si-al-barrier-ready"
    assert manifest["summary"]["matched"]["ordering"] == "Si-O-Al < Si-O-Si"


def test_al_run_path_requires_exact_parent_si_terminal(tmp_path):
    with pytest.raises(ValueError, match="parent Si terminal"):
        relay.run_path(tmp_path, CHEAP, relay.Bounds(), model="al")


def test_recomputed_si_rejection_revokes_stale_al_terminal(monkeypatch, tmp_path):
    statuses = {"si": "accepted", "al": "accepted"}

    def fake_run(run_dir, settings, bounds, *, model, **kwargs):
        return _write_fake_terminal(run_dir, model, statuses[model])

    monkeypatch.setattr(relay, "run_path", fake_run)
    first = relay.run_campaign(
        tmp_path,
        CHEAP,
        relay.Bounds(),
        log_path=str(tmp_path / "run.log"),
        source_reactant=None,
        source_reactant_sha256=None,
    )
    assert sorted(first["paths"]) == [
        f"al:4w:{BRIDGE_SIDE_ACID_FAMILY}",
        f"si:4w:{BRIDGE_SIDE_ACID_FAMILY}",
    ]

    statuses["si"] = "rejected"
    second = relay.run_campaign(
        tmp_path,
        CHEAP,
        relay.Bounds(),
        log_path=str(tmp_path / "run.log"),
        source_reactant=None,
        source_reactant_sha256=None,
    )

    al_path = shared.path_directory(
        tmp_path, "al", BRIDGE_SIDE_ACID_WATER_COUNT, BRIDGE_SIDE_ACID_FAMILY
    )
    assert list(second["paths"]) == [f"si:4w:{BRIDGE_SIDE_ACID_FAMILY}"]
    assert second["summary"]["al_status"] is None
    assert not (al_path / "terminal.json").exists()
    assert len(list(al_path.glob("terminal.invalid-parent-*.json"))) == 1


@pytest.mark.parametrize("status", ["failed", "blocked"])
def test_computational_failure_remains_incomplete(monkeypatch, tmp_path, status):
    monkeypatch.setattr(
        relay,
        "run_path",
        lambda run_dir, settings, bounds, *, model, **kwargs: _write_fake_terminal(
            run_dir, model, status
        ),
    )

    manifest = relay.run_campaign(
        tmp_path,
        CHEAP,
        relay.Bounds(),
        log_path=str(tmp_path / "run.log"),
        source_reactant=None,
        source_reactant_sha256=None,
    )

    assert manifest["summary"]["verdict"].startswith("incomplete-")


def test_existing_a1f_manifest_and_source_hash_drift_are_refused(tmp_path):
    payload = a1f._expected_manifest(
        CHEAP,
        a1f.Bounds(),
        str(tmp_path / "run.log"),
    )
    shared.atomic_json(tmp_path / "manifest.json", payload)

    with pytest.raises(RuntimeError, match="different mechanism_version"):
        relay.run_campaign(
            tmp_path,
            CHEAP,
            relay.Bounds(),
            log_path=str(tmp_path / "run.log"),
            source_reactant=None,
            source_reactant_sha256=None,
        )

    tmp_path.joinpath("manifest.json").unlink()
    expected = relay._expected_manifest(
        CHEAP,
        relay.Bounds(),
        str(tmp_path / "run.log"),
        source_reactant=tmp_path / "source.xyz",
        source_reactant_sha256="a" * 64,
    )
    shared.atomic_json(tmp_path / "manifest.json", expected)
    with pytest.raises(RuntimeError, match="different source_reactant_sha256"):
        relay.run_campaign(
            tmp_path,
            CHEAP,
            relay.Bounds(),
            log_path=str(tmp_path / "run.log"),
            source_reactant=tmp_path / "source.xyz",
            source_reactant_sha256="b" * 64,
        )
