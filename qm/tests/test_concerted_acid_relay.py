from dataclasses import asdict

import numpy as np
import pytest

from quarry.clusters import (
    CONCERTED_ACID_FAMILIES,
    CONCERTED_ACID_WATER_COUNTS,
    Cluster,
    aluminosilicate_dimer,
    concerted_acid_relay_endpoints,
    disilicate,
    water,
)
from quarry.pipeline import DftSettings, FrequencyResult
from scripts import concerted_acid_relay as relay
from scripts import phase1_xiao_lasaga as phase1

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
@pytest.mark.parametrize("n_water", CONCERTED_ACID_WATER_COUNTS)
@pytest.mark.parametrize("family", CONCERTED_ACID_FAMILIES)
def test_concerted_endpoint_seeds_have_exact_typed_basins_and_no_collisions(
    model, n_water, family
):
    ends = concerted_acid_relay_endpoints(model(), n_water=n_water, family=family)
    protons = relay.proton_indices(ends, n_water)

    assert ends.reactant.symbols == ends.product.symbols
    assert ends.reactant.frozen_indices == ends.product.frozen_indices == []
    assert np.all(np.isfinite(ends.reactant.coords))
    assert np.all(np.isfinite(ends.product.coords))
    assert relay.minimum_pair_distance(ends.reactant) >= relay.MIN_PAIR_DISTANCE_A
    assert relay.minimum_pair_distance(ends.product) >= relay.MIN_PAIR_DISTANCE_A
    assert (
        relay.endpoint_gate_reason(
            ends.reactant, ends, n_water=n_water, endpoint="reactant"
        )
        is None
    )
    assert (
        relay.endpoint_gate_reason(
            ends.product, ends, n_water=n_water, endpoint="product"
        )
        is None
    )
    assert phase1.acid_basin_signature(
        ends.reactant,
        ends.ow_index,
        protons,
        solvent_oxygen_indices=ends.solvent_oxygen_indices,
    ) == relay.expected_basin(n_water, "reactant")
    assert phase1.acid_basin_signature(
        ends.product,
        ends.ow_index,
        protons,
        solvent_oxygen_indices=ends.solvent_oxygen_indices,
    ) == relay.expected_basin(n_water, "product")


@pytest.mark.parametrize("n_water", CONCERTED_ACID_WATER_COUNTS)
def test_bridge_chain_moves_every_physical_relay_h_without_losing_ownership(n_water):
    ends = concerted_acid_relay_endpoints(
        disilicate(), n_water=n_water, family="bridge-donor-chain"
    )
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
    assert len(ends.relay_h_indices) == n_water
    assert all(
        reactant_owners[index] != product_owners[index]
        for index in ends.relay_h_indices
    )
    assert product_owners[ends.transferred_h_index] == 0
    assert list(reactant_owners.values()).count(ends.ow_index) == 3
    assert list(product_owners.values()).count(ends.ow_index) == 2


@pytest.mark.parametrize("n_water", [1, 2, 5, 6])
def test_concerted_builder_refuses_water_counts_outside_exact_scope(n_water):
    with pytest.raises(ValueError, match="exactly 3 or 4"):
        concerted_acid_relay_endpoints(
            disilicate(), n_water=n_water, family="bridge-donor-chain"
        )


def test_concerted_builder_refuses_retired_or_unscoped_families():
    with pytest.raises(ValueError, match="family must be one of"):
        concerted_acid_relay_endpoints(
            disilicate(), n_water=3, family="attacker-centered-ring"
        )


@pytest.mark.parametrize("n_water", CONCERTED_ACID_WATER_COUNTS)
@pytest.mark.parametrize("family", CONCERTED_ACID_FAMILIES)
def test_endpoint_reaction_vector_passes_coupled_mode_gate(n_water, family):
    ends = concerted_acid_relay_endpoints(disilicate(), n_water=n_water, family=family)
    mode = ends.product.coords - ends.reactant.coords
    result = relay.coupled_mode_components(
        ends.reactant, _frequency(ends.reactant, mode), ends
    )

    assert result["accepted"] is True
    assert result["components"]["si_obr_cleavage"] > 0.0
    assert result["components"]["si_ow_attack"] < 0.0
    assert result["components"]["relay_to_obr"] < 0.0


def test_coupled_mode_gate_rejects_a_water_wag():
    ends = concerted_acid_relay_endpoints(
        disilicate(), n_water=3, family="bridge-donor-chain"
    )
    mode = np.zeros_like(ends.reactant.coords)
    mode[ends.ow_index + 2, 2] = 1.0

    result = relay.coupled_mode_components(
        ends.reactant, _frequency(ends.reactant, mode), ends
    )

    assert result["accepted"] is False
    assert "does not couple" in result["reason"]


def test_full_irc_reuses_one_optimizer_for_both_directions(monkeypatch, tmp_path):
    import sella

    calls = []
    instances = []

    class FakeIRC:
        def __init__(self, atoms, **kwargs):
            self.atoms = atoms
            self.kwargs = kwargs
            instances.append(self)

        def run(self, **kwargs):
            calls.append(kwargs)
            sign = 1.0 if kwargs["direction"] == "forward" else -1.0
            self.atoms.positions[:] = water().coords + sign * 0.05
            return True

    monkeypatch.setattr(sella, "IRC", FakeIRC)
    monkeypatch.setattr(relay, "full_irc", relay.full_irc)
    import quarry.ts as ts_module

    monkeypatch.setattr(ts_module, "make_ase_calculator", lambda *args: object())
    backward, forward = ts_module.full_irc(
        water(),
        CHEAP,
        max_steps=11,
        trajectory=tmp_path / "irc.traj",
        logfile=tmp_path / "irc.log",
    )

    assert len(instances) == 1
    assert [call["direction"] for call in calls] == ["forward", "reverse"]
    assert all(call["steps"] == 11 for call in calls)
    assert not np.allclose(backward.coords, forward.coords)


def test_full_irc_fails_closed_when_either_direction_exhausts(monkeypatch):
    import sella

    import quarry.ts as ts_module

    class FakeIRC:
        def __init__(self, atoms, **kwargs):
            self.atoms = atoms

        def run(self, **kwargs):
            return kwargs["direction"] == "forward"

    monkeypatch.setattr(sella, "IRC", FakeIRC)
    monkeypatch.setattr(ts_module, "make_ase_calculator", lambda *args: object())
    with pytest.raises(RuntimeError, match="reverse direction did not converge"):
        ts_module.full_irc(water(), CHEAP, max_steps=7)


def test_terminal_resume_binds_versions_settings_bounds_and_artifact_hash(tmp_path):
    bounds = relay.Bounds()
    path_dir = tmp_path / "path"
    path_dir.mkdir()
    artifact = path_dir / "endpoint.xyz"
    artifact.write_text("proof\n")
    base = {
        "schema_version": relay.SCHEMA_VERSION,
        "mechanism_version": relay.MECHANISM_VERSION,
        "gate_version": relay.GATE_VERSION,
        "key": "si:3w:bridge-donor-chain",
        "settings": asdict(CHEAP),
        "bounds": asdict(bounds),
        "status": "rejected",
        "artifacts": {"endpoint.xyz": relay.sha256_path(artifact)},
    }
    relay.atomic_json(path_dir / "terminal.json", base)

    assert relay.terminal_record_is_reusable(
        path_dir,
        key="si:3w:bridge-donor-chain",
        settings=CHEAP,
        bounds=bounds,
    )
    artifact.write_text("drift\n")
    assert not relay.terminal_record_is_reusable(
        path_dir,
        key="si:3w:bridge-donor-chain",
        settings=CHEAP,
        bounds=bounds,
    )
    artifact.write_text("proof\n")
    base["mechanism_version"] += 1
    relay.atomic_json(path_dir / "terminal.json", base)
    assert not relay.terminal_record_is_reusable(
        path_dir,
        key="si:3w:bridge-donor-chain",
        settings=CHEAP,
        bounds=bounds,
    )


def test_optimizer_failure_is_incomplete_not_scientific_rejection(
    monkeypatch, tmp_path
):
    def explode(*args, **kwargs):
        raise RuntimeError("optimizer detonated")

    monkeypatch.setattr(relay, "optimize_bounded", explode)
    record = relay.run_path(
        tmp_path,
        CHEAP,
        relay.Bounds(),
        model="si",
        n_water=3,
        family="bridge-donor-chain",
    )

    assert record["status"] == "failed"
    assert record["stage"] == "incomplete-optimize-reactant"
    assert "optimizer detonated" in record["reason"]
    assert (tmp_path / "si/3w/bridge-donor-chain/terminal.json").is_file()


def _write_fake_terminal(run_dir, model, n_water, family, status):
    path_dir = relay.path_directory(run_dir, model, n_water, family)
    path_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "key": relay.path_key(model, n_water, family),
        "status": status,
        "barrier_kj_mol": 100.0 if model == "si" else 80.0,
    }
    relay.atomic_json(path_dir / "terminal.json", payload)
    return payload


def test_campaign_runs_exact_four_si_paths_and_never_unmatched_al(
    monkeypatch, tmp_path
):
    calls = []

    def fake_run(run_dir, settings, bounds, *, model, n_water, family):
        calls.append((model, n_water, family))
        return _write_fake_terminal(run_dir, model, n_water, family, "rejected")

    monkeypatch.setattr(relay, "run_path", fake_run)
    manifest = relay.run_campaign(
        tmp_path, CHEAP, relay.Bounds(), log_path=str(tmp_path / "run.log")
    )

    assert calls == [
        ("si", 3, "bridge-donor-chain"),
        ("si", 3, "compact-cyclic-relay"),
        ("si", 4, "bridge-donor-chain"),
        ("si", 4, "compact-cyclic-relay"),
    ]
    assert manifest["summary"]["verdict"] == ("all-four-si-paths-conclusive-rejection")
    assert manifest["summary"]["next_mechanism_card_required"] is True


def test_campaign_triggers_only_the_exact_al_match_after_accepted_si(
    monkeypatch, tmp_path
):
    calls = []

    def fake_run(run_dir, settings, bounds, *, model, n_water, family):
        calls.append((model, n_water, family))
        accepted = n_water == 3 and family == "compact-cyclic-relay"
        status = "accepted" if accepted else "rejected"
        return _write_fake_terminal(run_dir, model, n_water, family, status)

    monkeypatch.setattr(relay, "run_path", fake_run)
    manifest = relay.run_campaign(
        tmp_path, CHEAP, relay.Bounds(), log_path=str(tmp_path / "run.log")
    )

    assert calls == [
        ("si", 3, "bridge-donor-chain"),
        ("si", 3, "compact-cyclic-relay"),
        ("al", 3, "compact-cyclic-relay"),
        ("si", 4, "bridge-donor-chain"),
        ("si", 4, "compact-cyclic-relay"),
    ]
    assert manifest["summary"]["verdict"] == "matched-si-al-barrier-ready"
    assert manifest["summary"]["matched"] == [
        {
            "water_count": 3,
            "family": "compact-cyclic-relay",
            "si_terminal": manifest["paths"]["si:3w:compact-cyclic-relay"],
            "al_terminal": manifest["paths"]["al:3w:compact-cyclic-relay"],
            "si_barrier_kj_mol": 100.0,
            "al_barrier_kj_mol": 80.0,
            "ordering": "Si-O-Al < Si-O-Si",
        }
    ]


@pytest.mark.parametrize("al_status", ["failed", "blocked"])
def test_campaign_is_incomplete_when_required_matched_al_fails_or_is_blocked(
    monkeypatch, tmp_path, al_status
):
    def fake_run(run_dir, settings, bounds, *, model, n_water, family):
        accepted_si = model == "si" and n_water == 3 and family == "compact-cyclic-relay"
        if accepted_si:
            status = "accepted"
        elif model == "al":
            status = al_status
        else:
            status = "rejected"
        return _write_fake_terminal(run_dir, model, n_water, family, status)

    monkeypatch.setattr(relay, "run_path", fake_run)
    manifest = relay.run_campaign(
        tmp_path, CHEAP, relay.Bounds(), log_path=str(tmp_path / "run.log")
    )

    assert manifest["summary"]["verdict"] == "incomplete-matched-al-campaign"
    assert manifest["summary"]["al_status_counts"][al_status] == 1
    assert manifest["summary"]["matched"] == []
    assert not manifest["summary"]["verdict"].startswith("si-path-accepted")


def test_existing_manifest_refuses_mechanism_version_drift(tmp_path):
    bounds = relay.Bounds()
    payload = relay._expected_manifest(CHEAP, bounds, str(tmp_path / "run.log"))
    payload["mechanism_version"] += 1
    relay.atomic_json(tmp_path / "manifest.json", payload)

    with pytest.raises(RuntimeError, match="different mechanism_version"):
        relay.run_campaign(tmp_path, CHEAP, bounds, log_path=str(tmp_path / "run.log"))
