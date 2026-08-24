"""Finite matched proton-relay ensemble orchestration gates."""

import json
import shutil
from dataclasses import asdict, replace
from types import SimpleNamespace

import numpy as np
import pytest

from quarry.clusters import disilicate, protonated_bridge_complex
from quarry.pipeline import DftSettings, FrequencyResult
from scripts import acid_microsolvation_ensemble as ensemble

TARGET_SEED = "si-acid:3w:bridge-donor-chain"


def _build_refinement_source(tmp_path):
    root = tmp_path / "source"
    root.mkdir()
    settings = DftSettings(xc="b3lyp", basis="def2-svp", density_fit=True, use_gpu=True)
    manifest = ensemble._expected_manifest(settings, 160, "/retired/ensemble.log")
    for n_water in ensemble.WATER_COUNTS:
        for family in ensemble.ACID_MICROSOLVATION_FAMILIES:
            key = ensemble._seed_key("si-acid", n_water, family)
            seed_dir = ensemble._seed_dir(root, "si-acid", n_water, family)
            seed_dir.mkdir(parents=True)
            cluster = ensemble._template("si-acid", n_water, family)
            for name in ("seed.xyz", "preoptimized.xyz", "optimized.xyz"):
                ensemble.phase1.save_xyz(cluster, seed_dir / name)
            is_target = key == TARGET_SEED
            manifest["seeds"][key] = {
                "status": "failed" if is_target else "rejected",
                "reaction": "si-acid",
                "water_count": n_water,
                "family": family,
                "reason": (
                    ensemble.OPTIMIZER_EXHAUSTION_REASON
                    if is_target
                    else "bridge proton transferred"
                ),
                "preoptimization_converged": False,
                "production_converged": not is_target,
                "artifacts": ensemble._artifact_hashes(
                    seed_dir, ("seed.xyz", "preoptimized.xyz", "optimized.xyz")
                ),
            }
    manifest["summary"] = {
        "verdict": "incomplete-si-screen",
        "si_seed_count": 16,
        "si_status_counts": {
            "accepted": 0,
            "blocked": 0,
            "failed": 1,
            "rejected": 15,
        },
    }
    manifest_path = root / "manifest.json"
    ensemble.atomic_json(manifest_path, manifest)
    (root / "ensemble.log").write_text("immutable source log\n")
    return manifest_path, ensemble.sha256_path(manifest_path)


def _load_refinement_source(tmp_path):
    manifest_path, manifest_sha256 = _build_refinement_source(tmp_path)
    source = ensemble.validate_refinement_source(
        manifest_path,
        expected_manifest_sha256=manifest_sha256,
        expected_log_sha256=ensemble.sha256_path(manifest_path.parent / "ensemble.log"),
        target_key=TARGET_SEED,
        attempt_dir=tmp_path / "attempt",
        max_steps=160,
    )
    return source, manifest_path


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


@pytest.mark.parametrize("bound", [True, 0, -1, 161, 1.5])
def test_refinement_bound_is_strictly_one_through_160(bound):
    with pytest.raises(ValueError, match="integer in 1..160"):
        ensemble._require_refinement_bound(bound)


def test_refinement_source_requires_pinned_hash_and_exact_artifacts(tmp_path):
    manifest_path, manifest_sha256 = _build_refinement_source(tmp_path)
    with pytest.raises(ValueError, match="manifest SHA-256 mismatch"):
        ensemble.validate_refinement_source(
            manifest_path,
            expected_manifest_sha256="0" * 64,
            expected_log_sha256=ensemble.sha256_path(
                manifest_path.parent / "ensemble.log"
            ),
            target_key=TARGET_SEED,
            attempt_dir=tmp_path / "attempt",
            max_steps=160,
        )

    target = (
        ensemble._seed_dir(manifest_path.parent, "si-acid", 3, "bridge-donor-chain")
        / "optimized.xyz"
    )
    target.write_text(target.read_text() + "\n")
    with pytest.raises(ValueError, match="artifact SHA-256 mismatch"):
        ensemble.validate_refinement_source(
            manifest_path,
            expected_manifest_sha256=manifest_sha256,
            expected_log_sha256=ensemble.sha256_path(
                manifest_path.parent / "ensemble.log"
            ),
            target_key=TARGET_SEED,
            attempt_dir=tmp_path / "attempt",
            max_steps=160,
        )


def test_refinement_source_refuses_settings_and_population_drift(tmp_path):
    manifest_path, _manifest_sha256 = _build_refinement_source(tmp_path)
    manifest = json.loads(manifest_path.read_text())
    manifest["settings"] = {**asdict(DftSettings()), "xc": "pbe"}
    ensemble.atomic_json(manifest_path, manifest)
    with pytest.raises(ValueError, match="settings mismatch"):
        ensemble.validate_refinement_source(
            manifest_path,
            expected_manifest_sha256=ensemble.sha256_path(manifest_path),
            expected_log_sha256=ensemble.sha256_path(
                manifest_path.parent / "ensemble.log"
            ),
            target_key=TARGET_SEED,
            attempt_dir=tmp_path / "attempt",
            max_steps=160,
        )

    source_root = manifest_path.parent
    shutil.rmtree(source_root)
    manifest_path, _manifest_sha256 = _build_refinement_source(tmp_path)
    manifest = json.loads(manifest_path.read_text())
    other_key = "si-acid:4w:bridge-donor-chain"
    manifest["seeds"][other_key]["status"] = "failed"
    ensemble.atomic_json(manifest_path, manifest)
    with pytest.raises(ValueError, match="sole failed"):
        ensemble.validate_refinement_source(
            manifest_path,
            expected_manifest_sha256=ensemble.sha256_path(manifest_path),
            expected_log_sha256=ensemble.sha256_path(
                manifest_path.parent / "ensemble.log"
            ),
            target_key=TARGET_SEED,
            attempt_dir=tmp_path / "attempt",
            max_steps=160,
        )


def test_refinement_restarts_exact_endpoint_once_and_preserves_source(
    tmp_path, monkeypatch
):
    source, manifest_path = _load_refinement_source(tmp_path)
    before = {
        path.relative_to(manifest_path.parent): ensemble.sha256_path(path)
        for path in manifest_path.parent.rglob("*")
        if path.is_file()
    }
    source_endpoint = (
        ensemble._seed_dir(manifest_path.parent, "si-acid", 3, "bridge-donor-chain")
        / "optimized.xyz"
    )
    expected = ensemble.phase1.load_xyz(
        source_endpoint, ensemble._template("si-acid", 3, "bridge-donor-chain")
    )
    optimizer_calls = []
    gate_calls = []

    def fake_optimize(cluster, settings, *, max_steps):
        optimizer_calls.append((cluster.coords.copy(), settings, max_steps))
        return SimpleNamespace(cluster=cluster, converged=True)

    def fake_gate(seed_dir, optimized, settings, **identity):
        gate_calls.append((seed_dir, optimized, settings, identity))
        return {
            "status": "rejected",
            "reason": "bridge proton transferred",
            "basin_signature": [False],
            "n_imaginary": None,
            "electronic_hartree": None,
            "minimum_pair_distance_a": ensemble.minimum_pair_distance(optimized),
        }

    monkeypatch.setattr(ensemble, "optimize_bounded", fake_optimize)
    monkeypatch.setattr(ensemble, "evaluate_optimized_endpoint", fake_gate)
    attempt_dir = tmp_path / "attempt"
    resolution = ensemble.run_failed_endpoint_refinement(
        source,
        attempt_dir=attempt_dir,
        max_steps=160,
        log_path="refinement.log",
    )

    assert len(optimizer_calls) == 1
    np.testing.assert_allclose(optimizer_calls[0][0], expected.coords)
    assert optimizer_calls[0][1] == source.settings
    assert optimizer_calls[0][2] == 160
    assert len(gate_calls) == 1
    assert gate_calls[0][3] == {
        "reaction": "si-acid",
        "n_water": 3,
        "family": "bridge-donor-chain",
    }
    assert resolution["summary"]["si_seed_count"] == 16
    assert resolution["summary"]["si_status_counts"]["rejected"] == 16
    assert resolution["summary"]["verdict"] == (
        "model-valid-no-go-pre-equilibrated-bridge"
    )
    assert len(resolution["effective_seed_refs"]) == 16
    assert resolution["effective_seed_refs"][TARGET_SEED]["source"] == (
        "a1d-refinement"
    )
    receipt = attempt_dir / resolution["attempt"]["receipt_path"]
    assert ensemble.sha256_path(receipt) == resolution["attempt"]["receipt_sha256"]
    after = {
        path.relative_to(manifest_path.parent): ensemble.sha256_path(path)
        for path in manifest_path.parent.rglob("*")
        if path.is_file()
    }
    assert after == before
    with pytest.raises(RuntimeError, match="second step budget"):
        ensemble.run_failed_endpoint_refinement(
            source,
            attempt_dir=attempt_dir,
            max_steps=160,
            log_path="refinement.log",
        )
    assert len(optimizer_calls) == 1


def test_refinement_exhaustion_remains_incomplete(tmp_path, monkeypatch):
    source, _manifest_path = _load_refinement_source(tmp_path)
    monkeypatch.setattr(
        ensemble,
        "optimize_bounded",
        lambda cluster, _settings, *, max_steps: SimpleNamespace(
            cluster=cluster, converged=False, max_steps=max_steps
        ),
    )

    resolution = ensemble.run_failed_endpoint_refinement(
        source,
        attempt_dir=tmp_path / "attempt",
        max_steps=160,
        log_path="refinement.log",
    )

    assert resolution["attempt"]["status"] == "failed"
    assert resolution["summary"]["si_status_counts"]["failed"] == 1
    assert resolution["summary"]["verdict"] == "incomplete-si-screen"


def test_refinement_gate_failure_does_not_relabel_optimizer_exhaustion(
    tmp_path, monkeypatch
):
    source, _manifest_path = _load_refinement_source(tmp_path)
    monkeypatch.setattr(
        ensemble,
        "optimize_bounded",
        lambda cluster, _settings, *, max_steps: SimpleNamespace(
            cluster=cluster, converged=True, max_steps=max_steps
        ),
    )
    monkeypatch.setattr(
        ensemble,
        "evaluate_optimized_endpoint",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            RuntimeError("frequency receipt invalid")
        ),
    )

    resolution = ensemble.run_failed_endpoint_refinement(
        source,
        attempt_dir=tmp_path / "attempt",
        max_steps=160,
        log_path="refinement.log",
    )

    receipt = json.loads(
        (tmp_path / "attempt" / resolution["attempt"]["receipt_path"]).read_text()
    )
    record = receipt["record"]
    assert record["production_converged"] is True
    assert record["failure_kind"] == "endpoint-evaluation-error"
    assert record["status"] == "failed"


def test_accepted_refinement_screens_only_exact_matched_al(tmp_path, monkeypatch):
    source, _manifest_path = _load_refinement_source(tmp_path)
    monkeypatch.setattr(
        ensemble,
        "optimize_bounded",
        lambda cluster, _settings, *, max_steps: SimpleNamespace(
            cluster=cluster, converged=True, max_steps=max_steps
        ),
    )
    monkeypatch.setattr(
        ensemble,
        "evaluate_optimized_endpoint",
        lambda *_args, **_kwargs: {
            "status": "accepted",
            "reason": None,
            "basin_signature": [False, True, True, 1, 6, 0],
            "n_imaginary": 0,
            "electronic_hartree": -100.0,
            "minimum_pair_distance_a": 1.0,
        },
    )
    al_calls = []

    def fake_al_screen(manifest, manifest_path, **kwargs):
        al_calls.append(kwargs)
        record = {
            "status": "rejected",
            "reaction": kwargs["reaction"],
            "water_count": kwargs["n_water"],
            "family": kwargs["family"],
            "reason": "matched Al basin rejected",
        }
        manifest["seeds"][
            ensemble._seed_key(kwargs["reaction"], kwargs["n_water"], kwargs["family"])
        ] = record
        ensemble.atomic_json(manifest_path, manifest)
        return record

    monkeypatch.setattr(ensemble, "screen_seed", fake_al_screen)
    resolution = ensemble.run_failed_endpoint_refinement(
        source,
        attempt_dir=tmp_path / "attempt",
        max_steps=160,
        log_path="refinement.log",
    )

    assert len(al_calls) == 1
    assert al_calls[0]["reaction"] == "al-acid"
    assert al_calls[0]["n_water"] == 3
    assert al_calls[0]["family"] == "bridge-donor-chain"
    assert resolution["summary"]["verdict"] == ("si-minimum-without-matched-al-minimum")
    assert resolution["summary"]["matched_receipt"] is None
