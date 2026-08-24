from dataclasses import asdict, replace

import numpy as np
import pytest

from quarry.clusters import (
    SEPARATED_ACID_FAMILY,
    SEPARATED_ACID_WATER_COUNT,
    Cluster,
    aluminosilicate_dimer,
    disilicate,
    separated_acid_relay_endpoints,
)
from quarry.pipeline import DftSettings, FrequencyResult
from scripts import concerted_acid_relay as shared
from scripts import phase1_xiao_lasaga as phase1
from scripts import separated_acid_relay as relay

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
def test_separated_endpoint_seeds_have_exact_roles_basins_and_no_collisions(model):
    ends = separated_acid_relay_endpoints(model())

    assert ends.donor_oxygen_index != ends.ow_index
    assert ends.solvent_oxygen_indices == (
        ends.donor_oxygen_index,
        ends.ow_index,
        ends.relay_oxygen_index,
        ends.spectator_oxygen_index,
    )
    assert len(ends.solvent_oxygen_indices) == SEPARATED_ACID_WATER_COUNT
    assert len(ends.solvent_h_indices) == 2 * SEPARATED_ACID_WATER_COUNT + 1
    assert ends.reactant.symbols == ends.product.symbols
    assert ends.reactant.frozen_indices == ends.product.frozen_indices == []
    for endpoint, cluster in (("reactant", ends.reactant), ("product", ends.product)):
        assert np.all(np.isfinite(cluster.coords))
        assert shared.minimum_pair_distance(cluster) >= relay.MIN_PAIR_DISTANCE_A
        assert (
            relay.endpoint_gate_reason(
                cluster,
                ends,
                n_water=SEPARATED_ACID_WATER_COUNT,
                endpoint=endpoint,
            )
            is None
        )
        assert phase1.acid_basin_signature(
            cluster,
            ends.ow_index,
            ends.solvent_h_indices,
            solvent_oxygen_indices=ends.solvent_oxygen_indices,
        ) == relay.expected_basin(SEPARATED_ACID_WATER_COUNT, endpoint)


def test_physical_h_relay_changes_exact_owners_and_keeps_attacker_neutral():
    ends = separated_acid_relay_endpoints(disilicate())
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
    assert reactant_owners[ends.donor_h_index] == ends.donor_oxygen_index
    assert product_owners[ends.donor_h_index] == ends.relay_oxygen_index
    assert reactant_owners[ends.transferred_h_index] == ends.relay_oxygen_index
    assert product_owners[ends.transferred_h_index] == phase1.BR_INDEX
    attacker_h = (ends.ow_index + 1, ends.ow_index + 2)
    assert all(reactant_owners[index] == ends.ow_index for index in attacker_h)
    assert all(product_owners[index] == ends.ow_index for index in attacker_h)
    assert len(ends.relay_h_indices) == 2


@pytest.mark.parametrize("n_water", [1, 3, 5, 6])
def test_separated_builder_refuses_every_water_count_except_four(n_water):
    with pytest.raises(ValueError, match="exactly 4 waters"):
        separated_acid_relay_endpoints(disilicate(), n_water=n_water)


def test_separated_builder_refuses_a1e_and_unscoped_families():
    for family in ("bridge-donor-chain", "compact-cyclic-relay", "other"):
        with pytest.raises(ValueError, match="family must be"):
            separated_acid_relay_endpoints(disilicate(), family=family)


def test_endpoint_reaction_vector_passes_all_six_coupled_mode_components():
    ends = separated_acid_relay_endpoints(disilicate())
    mode = ends.product.coords - ends.reactant.coords

    result = relay.coupled_mode_components(
        ends.reactant,
        _frequency(ends.reactant, mode),
        ends,
    )

    assert result["accepted"] is True
    assert result["components"]["si_obr_cleavage"] > 0.0
    assert result["components"]["si_ow_attack"] < 0.0
    assert result["components"]["donor_h_release"] > 0.0
    assert result["components"]["donor_h_to_relay"] < 0.0
    assert result["components"]["relay_h_release"] > 0.0
    assert result["components"]["relay_to_obr"] < 0.0


def test_coupled_mode_gate_rejects_a_water_wag():
    ends = separated_acid_relay_endpoints(disilicate())
    mode = np.zeros_like(ends.reactant.coords)
    mode[ends.ow_index + 1, 2] = 1.0

    result = relay.coupled_mode_components(
        ends.reactant,
        _frequency(ends.reactant, mode),
        ends,
    )

    assert result["accepted"] is False
    assert "does not couple" in result["reason"]


def test_endpoint_gate_rejects_unassigned_and_double_owned_hydrogen():
    ends = separated_acid_relay_endpoints(disilicate())
    unassigned = ends.reactant.coords.copy()
    unassigned[ends.donor_h_index] += np.array([20.0, 0.0, 0.0])
    reason = relay.endpoint_gate_reason(
        replace(ends.reactant, coords=unassigned),
        ends,
        n_water=SEPARATED_ACID_WATER_COUNT,
        endpoint="reactant",
    )
    assert reason is not None

    ambiguous = ends.reactant.coords.copy()
    ambiguous[ends.donor_h_index] = 0.5 * (
        ambiguous[ends.donor_oxygen_index] + ambiguous[ends.relay_oxygen_index]
    )
    reason = relay.endpoint_gate_reason(
        replace(ends.reactant, coords=ambiguous),
        ends,
        n_water=SEPARATED_ACID_WATER_COUNT,
        endpoint="reactant",
    )
    assert reason is not None


def test_full_irc_requires_exact_typed_endpoints():
    ends = separated_acid_relay_endpoints(disilicate())

    assert (
        relay.irc_channel_reason(
            ends.reactant,
            ends.product,
            ends,
            n_water=SEPARATED_ACID_WATER_COUNT,
        )
        is None
    )
    reason = relay.irc_channel_reason(
        ends.reactant,
        ends.reactant,
        ends,
        n_water=SEPARATED_ACID_WATER_COUNT,
    )
    assert reason is not None
    assert "full IRC endpoint basins" in reason


def test_a1f_terminal_refuses_a1e_mechanism_and_gate_identity(tmp_path):
    bounds = relay.Bounds()
    path_dir = tmp_path / "path"
    path_dir.mkdir()
    artifact = path_dir / "endpoint.xyz"
    artifact.write_text("proof\n")
    payload = {
        "schema_version": shared.SCHEMA_VERSION,
        "mechanism_version": shared.MECHANISM_VERSION,
        "gate_version": shared.GATE_VERSION,
        "key": f"si:4w:{SEPARATED_ACID_FAMILY}",
        "settings": asdict(CHEAP),
        "bounds": asdict(bounds),
        "status": "rejected",
        "artifacts": {"endpoint.xyz": shared.sha256_path(artifact)},
    }
    shared.atomic_json(path_dir / "terminal.json", payload)

    assert not shared.terminal_record_is_reusable(
        path_dir,
        key=f"si:4w:{SEPARATED_ACID_FAMILY}",
        settings=CHEAP,
        bounds=bounds,
        mechanism_version=relay.MECHANISM_VERSION,
        gate_version=relay.GATE_VERSION,
    )


def test_run_path_forwards_a1f_identity_builder_and_gates(monkeypatch, tmp_path):
    captured = {}

    def fake_run(*args, **kwargs):
        captured.update(kwargs)
        return {"status": "rejected"}

    monkeypatch.setattr(shared, "run_path", fake_run)
    relay.run_path(tmp_path, CHEAP, relay.Bounds(), model="si")

    assert captured["n_water"] == SEPARATED_ACID_WATER_COUNT
    assert captured["family"] == SEPARATED_ACID_FAMILY
    assert captured["mechanism_version"] == relay.MECHANISM_VERSION
    assert captured["gate_version"] == relay.GATE_VERSION
    assert captured["endpoint_builder"] is separated_acid_relay_endpoints
    assert captured["endpoint_gate"] is relay.endpoint_gate_reason
    assert captured["mode_gate"] is relay.coupled_mode_components
    assert captured["irc_gate"] is relay.irc_channel_reason


def _write_fake_terminal(run_dir, model, status):
    path_dir = shared.path_directory(
        run_dir,
        model,
        SEPARATED_ACID_WATER_COUNT,
        SEPARATED_ACID_FAMILY,
    )
    path_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "key": shared.path_key(
            model,
            SEPARATED_ACID_WATER_COUNT,
            SEPARATED_ACID_FAMILY,
        ),
        "model": model,
        "water_count": SEPARATED_ACID_WATER_COUNT,
        "family": SEPARATED_ACID_FAMILY,
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
    )

    assert calls == ["si"]
    assert manifest["summary"]["verdict"] == "si-path-conclusive-rejection"
    assert manifest["summary"]["al_status"] is None


def test_campaign_triggers_only_atom_matched_al_after_accepted_si(
    monkeypatch, tmp_path
):
    calls = []

    def fake_run(run_dir, settings, bounds, *, model, **kwargs):
        calls.append(model)
        return _write_fake_terminal(run_dir, model, "accepted")

    monkeypatch.setattr(relay, "run_path", fake_run)
    manifest = relay.run_campaign(
        tmp_path,
        CHEAP,
        relay.Bounds(),
        log_path=str(tmp_path / "run.log"),
    )

    assert calls == ["si", "al"]
    assert manifest["summary"]["verdict"] == "matched-si-al-barrier-ready"
    assert manifest["summary"]["matched"]["ordering"] == "Si-O-Al < Si-O-Si"


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
    )

    assert manifest["summary"]["verdict"].startswith("incomplete-")


def test_existing_a1e_manifest_is_refused(tmp_path):
    payload = shared._expected_manifest(
        CHEAP,
        shared.Bounds(),
        str(tmp_path / "run.log"),
    )
    shared.atomic_json(tmp_path / "manifest.json", payload)

    with pytest.raises(RuntimeError, match="different mechanism_version"):
        relay.run_campaign(
            tmp_path,
            CHEAP,
            relay.Bounds(),
            log_path=str(tmp_path / "run.log"),
        )
