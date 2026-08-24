"""Finite matched proton-relay ensemble orchestration gates."""

from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from quarry.clusters import disilicate, protonated_bridge_complex
from quarry.pipeline import DftSettings, FrequencyResult
from scripts import acid_microsolvation_ensemble as ensemble


def test_deduplicate_basins_uses_heavy_atom_rmsd_and_energy():
    first = protonated_bridge_complex(
        disilicate(), n_water=4, conformer_family="bridge-donor-chain"
    )
    duplicate = replace(first, name="duplicate", coords=first.coords.copy())
    duplicate.coords[-1] += np.array([0.0, 0.0, 0.5])
    distinct = replace(first, name="distinct", coords=first.coords.copy())
    shell_oxygen = len(disilicate().symbols) + 4
    distinct.coords[shell_oxygen] += np.array([2.0, 0.0, 0.0])

    candidates = [
        ensemble.BasinCandidate("first", 4, first, -100.0),
        ensemble.BasinCandidate("duplicate", 4, duplicate, -100.0001),
        ensemble.BasinCandidate("distinct", 4, distinct, -100.0),
    ]

    assignments = ensemble.deduplicate_basins(candidates)

    assert assignments == {
        "first": "duplicate",
        "duplicate": "duplicate",
        "distinct": "distinct",
    }


def test_deduplicate_basins_never_compares_different_water_counts():
    three = protonated_bridge_complex(
        disilicate(), n_water=3, conformer_family="compact-cyclic-relay"
    )
    four = protonated_bridge_complex(
        disilicate(), n_water=4, conformer_family="compact-cyclic-relay"
    )

    assignments = ensemble.deduplicate_basins(
        [
            ensemble.BasinCandidate("three", 3, three, -100.0),
            ensemble.BasinCandidate("four", 4, four, -100.0),
        ]
    )

    assert assignments == {"three": "three", "four": "four"}


def test_terminal_record_requires_hash_bound_artifacts(tmp_path):
    artifact = tmp_path / "optimized.xyz"
    artifact.write_text("durable evidence")
    record = {
        "status": "rejected",
        "artifacts": {
            "optimized.xyz": ensemble.sha256_path(artifact),
        },
    }

    assert ensemble.terminal_record_is_reusable(tmp_path, record)
    artifact.write_text("tampered")
    assert not ensemble.terminal_record_is_reusable(tmp_path, record)


def test_accepted_terminal_record_requires_minimum_receipt_and_finite_energy(tmp_path):
    optimized = tmp_path / "optimized.xyz"
    optimized.write_text("geometry")
    record = {
        "status": "accepted",
        "n_imaginary": 0,
        "electronic_hartree": -100.0,
        "artifacts": {"optimized.xyz": ensemble.sha256_path(optimized)},
    }

    assert not ensemble.terminal_record_is_reusable(tmp_path, record)
    minimum = tmp_path / "minimum.json"
    minimum.write_text("receipt")
    for name in ("seed.xyz", "preoptimized.xyz", "minimum.json"):
        path = tmp_path / name
        if not path.exists():
            path.write_text(name)
        record["artifacts"][name] = ensemble.sha256_path(path)
    assert ensemble.terminal_record_is_reusable(tmp_path, record)
    record["electronic_hartree"] = float("nan")
    assert not ensemble.terminal_record_is_reusable(tmp_path, record)


def test_screen_seed_rejects_exhausted_optimization_without_frequency(
    tmp_path, monkeypatch
):
    settings = DftSettings(xc="b3lyp", basis="def2-svp", density_fit=True)
    manifest_path = tmp_path / "manifest.json"
    manifest = {"seeds": {}}
    frequency_called = False

    def forbidden_frequency(*_args, **_kwargs):
        nonlocal frequency_called
        frequency_called = True
        raise AssertionError("frequency must not run after exhausted optimization")

    monkeypatch.setattr(
        ensemble,
        "optimize_bounded",
        lambda cluster, _settings, max_steps: SimpleNamespace(
            cluster=cluster, converged=False
        ),
    )
    monkeypatch.setattr(ensemble, "frequencies", forbidden_frequency)

    record = ensemble.screen_seed(
        manifest,
        manifest_path,
        reaction="si-acid",
        n_water=3,
        family="bridge-donor-chain",
        settings=settings,
        max_steps=1,
    )

    assert record["status"] == "failed"
    assert "did not converge" in record["reason"]
    assert not frequency_called


def test_screen_seed_rejects_malformed_frequency_receipt(tmp_path, monkeypatch):
    settings = DftSettings(xc="b3lyp", basis="def2-svp", density_fit=True)
    manifest_path = tmp_path / "manifest.json"
    manifest = {"seeds": {}}
    monkeypatch.setattr(
        ensemble,
        "optimize_bounded",
        lambda cluster, _settings, max_steps: SimpleNamespace(
            cluster=cluster, converged=True
        ),
    )
    monkeypatch.setattr(
        ensemble,
        "frequencies",
        lambda _cluster, _settings: FrequencyResult(
            frequencies_cm=np.array([100.0]),
            imaginary_cm=np.array([]),
            electronic_hartree=float("nan"),
            molar_mass_kg=0.1,
            rotational_temperatures_k=(1.0, 2.0, 3.0),
            linear=False,
        ),
    )

    record = ensemble.screen_seed(
        manifest,
        manifest_path,
        reaction="si-acid",
        n_water=3,
        family="bridge-donor-chain",
        settings=settings,
        max_steps=1,
    )

    assert record["status"] == "failed"
    assert record.get("electronic_hartree") is None


def test_run_ensemble_exhausts_all_si_seeds_before_no_go(tmp_path, monkeypatch):
    calls = []

    def fake_screen(manifest, manifest_path, **kwargs):
        calls.append((kwargs["reaction"], kwargs["n_water"], kwargs["family"]))
        key = ensemble._seed_key(
            kwargs["reaction"], kwargs["n_water"], kwargs["family"]
        )
        record = {
            "status": "rejected",
            "reaction": kwargs["reaction"],
            "water_count": kwargs["n_water"],
            "family": kwargs["family"],
            "reason": "bridge proton transferred",
        }
        manifest["seeds"][key] = record
        ensemble.atomic_json(manifest_path, manifest)
        return record

    monkeypatch.setattr(ensemble, "screen_seed", fake_screen)

    manifest = ensemble.run_ensemble(
        tmp_path,
        DftSettings(xc="b3lyp", basis="def2-svp", density_fit=True),
        max_steps=100,
        log_path="campaign.log",
    )

    assert len(calls) == 16
    assert {reaction for reaction, _waters, _family in calls} == {"si-acid"}
    assert manifest["summary"]["si_seed_count"] == 16
    assert manifest["summary"]["si_status_counts"]["rejected"] == 16
    assert manifest["summary"]["verdict"] == (
        "model-valid-no-go-pre-equilibrated-bridge"
    )


@pytest.mark.parametrize(
    ("al_status", "expected_verdict"),
    [
        ("accepted", "matched-minimum-ready-for-barrier-ladder"),
        ("failed", "incomplete-matched-screen"),
    ],
)
def test_run_ensemble_screens_unique_si_minimum_on_matched_al(
    tmp_path, monkeypatch, al_status, expected_verdict
):
    accepted_family = "bridge-donor-chain"

    def fake_screen(manifest, manifest_path, **kwargs):
        key = ensemble._seed_key(
            kwargs["reaction"], kwargs["n_water"], kwargs["family"]
        )
        eligible = kwargs["family"] == accepted_family and kwargs["n_water"] == 3
        status = "accepted" if eligible else "rejected"
        if eligible and kwargs["reaction"] == "al-acid":
            status = al_status
        record = {
            "status": status,
            "reaction": kwargs["reaction"],
            "water_count": kwargs["n_water"],
            "family": kwargs["family"],
            "reason": None if status == "accepted" else "bounded evaluation failed",
            "electronic_hartree": -100.0 if status == "accepted" else None,
        }
        manifest["seeds"][key] = record
        ensemble.atomic_json(manifest_path, manifest)
        return record

    def fake_load(_run_dir, *, reaction, n_water, family, record, settings):
        assert settings.xc == "b3lyp"
        return ensemble.BasinCandidate(
            ensemble._seed_key(reaction, n_water, family),
            n_water,
            protonated_bridge_complex(
                disilicate(), n_water=n_water, conformer_family=family
            ),
            record["electronic_hartree"],
        )

    monkeypatch.setattr(ensemble, "screen_seed", fake_screen)
    monkeypatch.setattr(ensemble, "_load_candidate", fake_load)
    monkeypatch.setattr(
        ensemble,
        "write_matched_candidate_receipt",
        lambda run_dir, *_args: run_dir / "matched-candidates.json",
    )

    manifest = ensemble.run_ensemble(
        tmp_path,
        DftSettings(xc="b3lyp", basis="def2-svp", density_fit=True),
        max_steps=100,
        log_path="campaign.log",
    )

    assert manifest["summary"]["verdict"] == expected_verdict
    if al_status != "accepted":
        assert manifest["summary"]["matched_candidates"] == []
        return
    assert len(manifest["summary"]["matched_candidates"]) == 1
    matched = manifest["summary"]["matched_candidates"][0]
    assert matched["water_count"] == 3
    assert matched["family"] == accepted_family
    assert matched["si_seed"] == f"si-acid:3w:{accepted_family}"
    assert matched["al_seed"] == f"al-acid:3w:{accepted_family}"
    assert len(matched["barrier_commands"]) == 2


def test_matched_receipt_binds_both_models_and_refuses_tampering(tmp_path):
    settings = DftSettings(xc="b3lyp", basis="def2-svp", density_fit=True)
    family = "bridge-donor-chain"
    n_water = 3
    manifest = {
        "seeds": {},
        "summary": {"matched_candidates": [{"water_count": n_water, "family": family}]},
    }
    for reaction, dimer, energy in (
        ("si-acid", disilicate(), -100.0),
        ("al-acid", ensemble.aluminosilicate_dimer(), -101.0),
    ):
        cluster = protonated_bridge_complex(
            dimer, n_water=n_water, conformer_family=family
        )
        seed_dir = ensemble._seed_dir(tmp_path, reaction, n_water, family)
        seed_dir.mkdir(parents=True)
        for name in ("seed.xyz", "preoptimized.xyz", "optimized.xyz"):
            ensemble.phase1.save_xyz(cluster, seed_dir / name)
        cluster = ensemble.phase1.load_xyz(seed_dir / "optimized.xyz", cluster)
        frequency = FrequencyResult(
            frequencies_cm=np.array([100.0, 200.0]),
            imaginary_cm=np.array([]),
            electronic_hartree=energy,
            molar_mass_kg=0.1,
            rotational_temperatures_k=(1.0, 2.0, 3.0),
            linear=False,
        )
        ensemble.phase1.write_reactant_minimum_receipt(
            seed_dir / "minimum.json", cluster, settings, frequency
        )
        key = ensemble._seed_key(reaction, n_water, family)
        manifest["seeds"][key] = {
            "status": "accepted",
            "reaction": reaction,
            "water_count": n_water,
            "family": family,
            "n_imaginary": 0,
            "electronic_hartree": energy,
            "artifacts": ensemble._artifact_hashes(
                seed_dir,
                ("seed.xyz", "preoptimized.xyz", "optimized.xyz", "minimum.json"),
            ),
        }
    manifest_path = tmp_path / "manifest.json"
    ensemble.atomic_json(manifest_path, manifest)
    receipt = ensemble.write_matched_candidate_receipt(
        tmp_path, manifest_path, manifest, settings
    )

    selected, minimum = ensemble.phase1.validate_matched_ensemble_receipt(
        receipt,
        reaction="si-acid",
        n_water=n_water,
        family=family,
        settings=settings,
    )

    assert (
        selected.symbols
        == protonated_bridge_complex(
            disilicate(), n_water=n_water, conformer_family=family
        ).symbols
    )
    assert minimum.name == "minimum.json"
    optimized = ensemble._seed_dir(tmp_path, "si-acid", n_water, family) / (
        "optimized.xyz"
    )
    optimized.write_text(optimized.read_text() + "\n")
    with pytest.raises(ValueError, match="geometry hash mismatch"):
        ensemble.phase1.validate_matched_ensemble_receipt(
            receipt,
            reaction="si-acid",
            n_water=n_water,
            family=family,
            settings=settings,
        )
