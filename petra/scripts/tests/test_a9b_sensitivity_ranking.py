from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import tomllib

MODULE_PATH = Path(__file__).resolve().parents[1] / "a9b_sensitivity_ranking.py"
SPEC = importlib.util.spec_from_file_location("a9b_sensitivity_ranking", MODULE_PATH)
assert SPEC and SPEC.loader
ranking = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ranking
SPEC.loader.exec_module(ranking)

VERIFIER_PATH = MODULE_PATH.with_name("verify_a9b_sensitivity_ranking.py")
VERIFIER_SPEC = importlib.util.spec_from_file_location(
    "verify_a9b_sensitivity_ranking", VERIFIER_PATH
)
assert VERIFIER_SPEC and VERIFIER_SPEC.loader
verifier = importlib.util.module_from_spec(VERIFIER_SPEC)
sys.modules[VERIFIER_SPEC.name] = verifier
VERIFIER_SPEC.loader.exec_module(verifier)

GUARD_PATH = MODULE_PATH.with_name("a9b_bounded_launch_guard.py")
GUARD_SPEC = importlib.util.spec_from_file_location(
    "a9b_bounded_launch_guard", GUARD_PATH
)
assert GUARD_SPEC and GUARD_SPEC.loader
launch_guard = importlib.util.module_from_spec(GUARD_SPEC)
sys.modules[GUARD_SPEC.name] = launch_guard
GUARD_SPEC.loader.exec_module(launch_guard)


def _json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode()


def _rehash_package(root: Path, scenario: str, seed: int) -> None:
    package = root / "runs" / scenario / f"seed-{seed}"
    internal_names = {
        "metadata.json",
        "initial.pgif.json",
        "events.jsonl",
        "likelihood.jsonl",
        "checkpoints.jsonl",
        "final.pgif.json",
    }
    internal = {
        "schema": "a9b-replica-receipt-v2",
        "status": "complete",
        "seed": seed,
        "files_sha256": {
            name: hashlib.sha256((package / name).read_bytes()).hexdigest()
            for name in sorted(internal_names)
        },
    }
    (package / "receipt.json").write_bytes(_json_bytes(internal))
    manifest = json.loads((root / "manifest.json").read_text())
    metadata = json.loads((package / "metadata.json").read_text())
    log = root / "logs" / scenario / f"seed-{seed}-attempt-1.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    log.write_text("synthetic test package\n", encoding="utf-8")
    deck_sha = next(
        row["deck_sha256"] for row in manifest["scenarios"] if row["name"] == scenario
    )
    external = {
        "schema": "a9b-job-receipt-v2",
        "scenario": scenario,
        "seed": seed,
        "attempt": 1,
        "command": [
            "nice",
            "-n",
            "10",
            metadata["inputs"]["runner_executable"],
            *metadata["argv"],
        ],
        "package": str(package.relative_to(root)),
        "log": str(log.relative_to(root)),
        "elapsed_seconds": 0.01,
        "runner_sha256": manifest["runner_sha256"],
        "deck_sha256": deck_sha,
        "log_sha256": hashlib.sha256(log.read_bytes()).hexdigest(),
        "package_files_sha256": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(package.iterdir())
        },
    }
    receipt = root / "receipts" / scenario / f"seed-{seed}.json"
    receipt.parent.mkdir(parents=True, exist_ok=True)
    receipt.write_bytes(_json_bytes(external))


def _write_jsonl(path: Path, rows: list[object]) -> None:
    path.write_text(
        "".join(
            json.dumps(row, sort_keys=True, allow_nan=False) + "\n" for row in rows
        ),
        encoding="utf-8",
    )


def _build_synthetic_package(root: Path, scenario: str, seed: int) -> Path:
    package = root / "runs" / scenario / f"seed-{seed}"
    package.mkdir(parents=True)
    manifest = json.loads((root / "manifest.json").read_text())
    contract = manifest["run_contract"]
    horizon = float(contract["horizon_s"])
    checkpoint_count = int(contract["checkpoints"])
    runner_path = root / "inputs" / "synthetic-a9b-runner"
    runner_path.write_bytes(b"synthetic runner fixture\n")
    runner_sha = hashlib.sha256(runner_path.read_bytes()).hexdigest()
    manifest["runner_sha256"] = runner_sha
    ranking.write_json_atomic(root / "manifest.json", manifest)
    scenario_row = next(row for row in manifest["scenarios"] if row["name"] == scenario)
    snapshot = json.loads((root / "inputs/initial-snapshot.pgif.json").read_text())
    initial = copy.deepcopy(snapshot)
    states = initial["nodes"]["columns"]["state"]["dict"]
    initial_data = initial["nodes"]["columns"]["state"]["data"]
    targets = json.loads((root / "inputs/path-report.json").read_text())["targets"]
    target = next(row for row in targets if row["kind"] == "Si")
    site = target["site"]
    old = initial_data[site]
    empty = states.index("Si.empty")
    final = copy.deepcopy(initial)
    final["nodes"]["columns"]["state"]["data"][site] = empty
    state_types = final["meta"]["petra"]["state_types"]
    type_dict = final["nodes"]["columns"]["type"]["dict"]
    final["nodes"]["columns"]["type"]["data"][site] = type_dict.index(
        state_types[empty]
    )
    final["meta"]["petra"]["step"] = 3
    final["meta"]["petra"]["time"] = horizon
    (package / "initial.pgif.json").write_text(
        json.dumps(initial, separators=(",", ":")), encoding="utf-8"
    )
    (package / "final.pgif.json").write_text(
        json.dumps(final, separators=(",", ":")), encoding="utf-8"
    )
    deck = tomllib.loads((root / "decks" / f"{scenario}.toml").read_text())
    reactions = [row["name"] for row in deck["reactions"]]
    reaction_ids = {name: index for index, name in enumerate(reactions)}
    mapping = [
        {
            "id": index,
            "name": name,
            "bias_factor": (
                1e6 if "hydrolysis" in name or name.startswith("desorb-") else 1.0
            ),
        }
        for index, name in enumerate(reactions)
    ]
    desorb = reaction_ids["desorb-si"]
    adsorb = reaction_ids["adsorb-si"]
    event_times = [horizon * 0.125, horizon * 0.275, horizon * 0.425]
    event_defs: list[tuple[int, int, int, str | None]] = [
        (desorb, old, empty, "Si"),
        (adsorb, empty, old, None),
        (desorb, old, empty, None),
    ]
    events: list[object] = [
        {
            "schema": "a9b-event-stream-v2",
            "seed": seed,
            "states": states,
            "state_occupants": state_types,
            "reactions": reactions,
        }
    ]
    likelihood: list[object] = [
        {"schema": "a9b-likelihood-stream-v2", "seed": seed, "formula": "fixture"}
    ]
    cumulative = 0.0
    start = 0.0
    events_by_time: dict[float, tuple[int, int, int, str | None]] = dict(
        zip(event_times, event_defs, strict=True)
    )
    checkpoint_times = [
        horizon * index / (checkpoint_count - 1) for index in range(checkpoint_count)
    ]
    segment_ends = sorted({*event_times, *checkpoint_times[1:]})
    for end in segment_ends:
        event = events_by_time.get(end)
        physical_rate = 1.0e12 if manifest["smoke"] else 1.0e6
        if event is None:
            reaction_id = None
            factor = None
            event_penalty = 0.0
            biased_rate = physical_rate * (1.0 + (seed % 1000) / 100.0)
        else:
            reaction_id, before, after, release = event
            factor = mapping[reaction_id]["bias_factor"]
            event_penalty = math.log(factor)
            biased_rate = physical_rate
        increment = (biased_rate - physical_rate) * (end - start) - event_penalty
        cumulative += increment
        segment_index = len(likelihood) - 1
        likelihood.append(
            {
                "segment": segment_index,
                "start_time_s": start,
                "end_time_s": end,
                "physical_total_rate_s-1": physical_rate,
                "biased_total_rate_s-1": biased_rate,
                "fired_reaction_id": reaction_id,
                "fired_bias_factor": factor,
                "log_likelihood_increment": increment,
                "cumulative_log_likelihood": cumulative,
            }
        )
        if event is not None:
            event_reaction_id, before, after, release = event
            event_index = len(events) - 1
            events.append(
                {
                    "event": event_index,
                    "step": event_index + 1,
                    "time_s": end,
                    "reaction_id": event_reaction_id,
                    "reaction": reactions[event_reaction_id],
                    "center_site": site,
                    "changes": [
                        {
                            "site": site,
                            "old_state_id": before,
                            "new_state_id": after,
                            "old_state": states[before],
                            "new_state": states[after],
                        }
                    ],
                    "likelihood_segment": segment_index,
                    "lattice_release_kind": release,
                }
            )
        start = end
    _write_jsonl(package / "events.jsonl", events)
    _write_jsonl(package / "likelihood.jsonl", likelihood)

    def counts(data: list[int]) -> dict[str, int]:
        return {name: data.count(index) for index, name in enumerate(states)}

    initial_counts = counts(initial_data)
    kind_totals = {
        kind: sum(
            count
            for name, count in initial_counts.items()
            if name.startswith(f"{kind}.")
        )
        for kind in {name.partition(".")[0] for name in states}
    }
    checkpoints = []
    for checkpoint_index, checkpoint_time in enumerate(checkpoint_times):
        checkpoint_state = list(initial_data)
        checkpoint_releases = {"Al": 0, "Si": 0}
        checkpoint_log = next(
            (
                row["cumulative_log_likelihood"]
                for row in reversed(likelihood[1:])
                if row["end_time_s"] <= checkpoint_time
            ),
            0.0,
        )
        checkpoint_step = 0
        for event_time, (reaction_id, _before, after, release) in zip(
            event_times, event_defs, strict=True
        ):
            if event_time > checkpoint_time:
                break
            checkpoint_state[site] = after
            checkpoint_step += 1
            if release is not None:
                checkpoint_releases[release] += 1
        area = (
            100.0
            if manifest["smoke"]
            else 100.00000000000001
            + checkpoint_index * 0.12345678901234567
            + (checkpoint_index % 3) * 1e-13
        )
        checkpoints.append(
            {
                "checkpoint": checkpoint_index,
                "checkpoint_role": "fixed_physical_time",
                "target_time_s": checkpoint_time,
                "actual_time_s": checkpoint_time,
                "step": checkpoint_step,
                "geometric_area_a2": area,
                "state_counts": counts(checkpoint_state),
                "kind_totals": kind_totals,
                "original_lattice_releases": checkpoint_releases,
                "cumulative_log_likelihood": checkpoint_log,
            }
        )
    _write_jsonl(package / "checkpoints.jsonl", checkpoints)
    deck_path = root / scenario_row["deck"]
    path_report = root / manifest["inputs"]["path_report"]
    lineage_snapshot = root / manifest["inputs"]["lineage_snapshot"]
    argv = [
        str(deck_path),
        str(path_report),
        str(lineage_snapshot),
        str(package),
        "--seed",
        str(seed),
        "--horizon",
        str(contract["horizon_s"]),
        "--checkpoints",
        str(contract["checkpoints"]),
        "--bias-factor",
        str(contract["bias_factor"]),
        "--max-events",
        str(contract["max_events"]),
    ]
    metadata = {
        "schema": "a9b-replica-package-v2",
        "survey_tier_platform_test": True,
        "not_production": True,
        "seed": seed,
        "temperature_k": 298.0,
        "horizon_s": horizon,
        "fixed_checkpoint_count": checkpoint_count,
        "recorded_checkpoint_count": checkpoint_count,
        "bias_factor": contract["bias_factor"],
        "max_events": contract["max_events"],
        "complete": True,
        "stop_reason": "physical_horizon",
        "event_count": 3,
        "likelihood_segment_count": len(likelihood) - 1,
        "final_log_likelihood": cumulative,
        "final_time_s": horizon,
        "final_step": 3,
        "original_lattice_releases": {"Al": 0, "Si": 1},
        "remaining_original_ion_ledger": 319,
        "reaction_mapping": mapping,
        "census": {
            "whole_finite_deck_original_centers": 320,
            "reachable_original_centers": 260,
            "topology_no_go_original_centers": 60,
            "whole_deck_estimand_includes_topology_no_go": True,
        },
        "inputs": {
            "deck": str(deck_path),
            "deck_sha256": scenario_row["deck_sha256"],
            "path_report": str(path_report),
            "path_report_sha256": manifest["inputs"]["path_report_sha256"],
            "lineage_snapshot": str(lineage_snapshot),
            "lineage_snapshot_sha256": manifest["inputs"]["lineage_snapshot_sha256"],
            "runner_executable": str(runner_path),
            "runner_executable_sha256": runner_sha,
        },
        "argv": argv,
    }
    (package / "metadata.json").write_bytes(_json_bytes(metadata))
    _rehash_package(root, scenario, seed)
    return package


def _synthetic_campaign(root: Path, *, smoke: bool = True) -> None:
    manifest = ranking.prepare_campaign(root, smoke=smoke)
    for row in manifest["required_runs"]:
        _build_synthetic_package(root, row["scenario"], row["seed"])
    manifest = json.loads((root / "manifest.json").read_text())
    manifest["status"] = "raw_complete"
    manifest["completed_run_count"] = len(manifest["required_runs"])
    ranking.write_json_atomic(root / "manifest.json", manifest)


class ScenarioContractTests(unittest.TestCase):
    def test_exact_29_scenarios_and_fixed_29x8_identity_matrix(self) -> None:
        scenarios = ranking.scenario_specs()
        self.assertEqual(len(scenarios), 29)
        self.assertEqual(scenarios[0].name, "nominal")
        self.assertEqual(len({item.name for item in scenarios}), 29)
        identities = ranking.required_identities(ranking.FIXED_SEEDS)
        self.assertEqual(len(identities), 232)
        self.assertEqual({seed for _, seed in identities}, set(ranking.FIXED_SEEDS))
        ranking.validate_identity_set(identities, smoke=False)
        with self.assertRaisesRegex(ValueError, "duplicate"):
            ranking.validate_identity_set([*identities, identities[0]], smoke=False)
        with self.assertRaisesRegex(ValueError, "missing"):
            ranking.validate_identity_set(identities[:-1], smoke=False)

    def test_existing_prepared_campaign_is_revalidated_before_resume(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = ranking.prepare_campaign(root, smoke=True)
            self.assertEqual(
                ranking.prepare_campaign(root, smoke=True)["scenario_count"], 29
            )
            deck = root / manifest["scenarios"][0]["deck"]
            deck.write_text(
                deck.read_text(encoding="utf-8") + "\n# drift\n", encoding="utf-8"
            )
            with self.assertRaisesRegex(ValueError, "artifact drift"):
                ranking.prepare_campaign(root, smoke=True)

    def test_resume_rejects_contradictory_tier_flags_and_temperature(self) -> None:
        for key, value in (
            ("survey_tier_platform_test", False),
            ("not_production", False),
            ("temperature_k", 999.0),
            ("production_kinetics", True),
        ):
            with self.subTest(key=key), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                manifest = ranking.prepare_campaign(root, smoke=True)
                manifest[key] = value
                ranking.write_json_atomic(root / "manifest.json", manifest)
                expected_error = (
                    "key set" if key == "production_kinetics" else f"drift in {key}"
                )
                with self.assertRaisesRegex(ValueError, expected_error):
                    ranking.prepare_campaign(root, smoke=True)

    def test_failed_preparation_never_publishes_partial_campaign_root(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "campaign"
            real_load = ranking._load_script

            def fail_during_staging(name: str) -> object:
                if name == "approximate_rate_closure":
                    raise RuntimeError("injected preparation failure")
                return real_load(name)

            with (
                mock.patch.object(
                    ranking, "_load_script", side_effect=fail_during_staging
                ),
                self.assertRaisesRegex(RuntimeError, "injected preparation failure"),
            ):
                ranking.prepare_campaign(root, smoke=True)
            self.assertFalse(root.exists())
            self.assertEqual(list(Path(directory).glob(".campaign.prepare-*")), [])

    def test_full_manifest_preregisters_exact_232_runs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            manifest = ranking.prepare_campaign(Path(directory), smoke=False)
            self.assertEqual(manifest["scenario_count"], 29)
            self.assertEqual(manifest["replicas_per_scenario"], 8)
            self.assertEqual(manifest["required_run_count"], 232)
            self.assertEqual(len(manifest["required_runs"]), 232)

    def test_only_ph4_is_simulated_and_ph3_ph5_are_comparison_bounds(self) -> None:
        scope = ranking.lab_comparison_scope({"Si": None, "Al": None})
        self.assertEqual(scope["simulated_ph"], 4)
        self.assertEqual(scope["comparison_bound_ph"], [3, 5])
        self.assertNotIn("prediction", str(scope).lower())

    def test_laboratory_comparison_reports_all_three_explicit_gaps(self) -> None:
        scope = ranking.lab_comparison_scope({"Si": 4e-12, "Al": 2e-12})
        gaps = scope["log10_gaps_vs_laboratory"]
        self.assertEqual(
            set(gaps),
            {
                "ph3_bound_not_simulated",
                "ph4_simulated",
                "ph5_bound_not_simulated",
            },
        )
        self.assertTrue(all(math.isfinite(value) for value in gaps.values()))
        self.assertEqual(scope, verifier._lab({"Si": 4e-12, "Al": 2e-12}))
        self.assertEqual(scope["simulated_kaolinite_formula_unit_rate_mol_m2_s"], 1e-12)
        self.assertEqual(
            scope["formula_unit_conversion"], "min(Si_mol_flux/2, Al_mol_flux/2)"
        )

    def test_laboratory_formula_unit_rate_is_censored_with_either_species(self) -> None:
        for rates in ({"Si": 1e-12, "Al": None}, {"Si": None, "Al": 1e-12}):
            with self.subTest(rates=rates):
                scope = ranking.lab_comparison_scope(rates)
                self.assertIsNone(
                    scope["simulated_kaolinite_formula_unit_rate_mol_m2_s"]
                )
                self.assertEqual(
                    scope["formula_unit_rate_status"],
                    "censored_species_rate_unavailable",
                )
                self.assertTrue(
                    all(
                        value is None
                        for value in scope["log10_gaps_vs_laboratory"].values()
                    )
                )
                self.assertEqual(scope, verifier._lab(rates))


class EstimatorGateTests(unittest.TestCase):
    def test_unnormalized_log_domain_estimator_and_ess(self) -> None:
        result = ranking.estimate_species([0.0] * 8, [1.0] * 8)
        self.assertEqual(result["status"], "estimated")
        self.assertEqual(result["mean"], 1.0)
        self.assertEqual(result["weight_ess"], 8.0)
        self.assertEqual(result["contribution_ess"], 8.0)

    def test_one_weight_domination_fails_weight_ess(self) -> None:
        result = ranking.estimate_species([100.0, *([0.0] * 15)], [1.0] * 16)
        self.assertEqual(result["status"], "censored_low_weight_ess")
        self.assertLess(result["weight_ess_fraction"], 0.1)

    def test_contribution_gate_is_separate_from_merged_sampling_contract(self) -> None:
        result = ranking.estimate_species([100.0, *([0.0] * 7)], [1.0] * 8)
        self.assertGreaterEqual(result["weight_ess_fraction"], 0.1)
        self.assertEqual(result["status"], "censored_low_contribution_ess")
        self.assertEqual(
            result["contribution_ess_gate_role"],
            "analysis_precision_ranking_gate_not_sampling_contract",
        )

    def test_low_contribution_ess_fails_positive_estimate(self) -> None:
        result = ranking.estimate_species([0.0] * 8, [1.0, *([0.0] * 7)])
        self.assertEqual(result["status"], "censored_low_contribution_ess")
        self.assertEqual(result["contribution_ess"], 1.0)

    def test_all_zero_is_typed_censor_without_poisson_bound(self) -> None:
        result = ranking.estimate_species([0.0] * 8, [0.0] * 8)
        self.assertEqual(result["status"], "censored_zero_biased_observations")
        self.assertIsNone(result["mean"])
        self.assertIsNone(result["upper_95"])
        self.assertEqual(result["poisson_bound"], "unavailable")
        self.assertFalse(result["rankable"])

    def test_invalid_nonfinite_weights_fail_closed(self) -> None:
        for bad in (math.nan, math.inf, -math.inf):
            with self.subTest(bad=bad), self.assertRaisesRegex(ValueError, "finite"):
                ranking.estimate_species([bad, *([0.0] * 7)], [1.0] * 8)

    def test_finite_but_unrepresentable_mean_is_typed_underflow(self) -> None:
        result = ranking.estimate_species([-1000.0] * 8, [1.0] * 8)
        self.assertEqual(result["status"], "censored_numerical_underflow")
        self.assertFalse(result["rankable"])
        self.assertIsNone(result["mean"])
        self.assertEqual(result["log_mean"], -1000.0)

    def test_finite_but_unrepresentable_mean_is_typed_overflow(self) -> None:
        for result in (
            ranking.estimate_species([1000.0] * 8, [1.0] * 8),
            verifier._estimate([1000.0] * 8, [1.0] * 8, 918273),
        ):
            self.assertEqual(result["status"], "censored_numerical_overflow")
            self.assertFalse(result["rankable"])
            self.assertIsNone(result["mean"])
            self.assertEqual(result["log_mean"], 1000.0)
            self.assertEqual(result["log_ci95"], [1000.0, 1000.0])

    def test_state_evolution_extreme_logs_are_typed_without_overflow(self) -> None:
        for point in (
            ranking._state_population_point([1000.0] * 8, [3.0] * 8),
            verifier._state_population_point([1000.0] * 8, [3.0] * 8),
        ):
            self.assertEqual(point["status"], "censored_numerical_overflow")
            self.assertIsNone(point["value"])
            self.assertAlmostEqual(point["log_mean"], 1000.0 + math.log(3.0))

    def test_finite_extreme_log_weights_are_overflow_safe_and_independent(self) -> None:
        cases = (
            (
                [1e308] * 8,
                "censored_numerical_overflow",
                8.0,
            ),
            (
                [-1e308] * 8,
                "censored_numerical_underflow",
                8.0,
            ),
            (
                [1e308, *([0.0] * 7)],
                "censored_low_contribution_ess",
                1.0,
            ),
            (
                [-1e308, *([0.0] * 7)],
                "estimated",
                7.0,
            ),
            (
                [1e308, -1e308, *([0.0] * 6)],
                "censored_low_contribution_ess",
                1.0,
            ),
        )
        for logs, status, expected_ess in cases:
            with self.subTest(logs=logs, status=status):
                analyzed = ranking.estimate_species(logs, [1.0] * 8, seed=13579)
                independently_regenerated = verifier._estimate(logs, [1.0] * 8, 13579)
                self.assertEqual(analyzed, independently_regenerated)
                self.assertEqual(analyzed["status"], status)
                self.assertEqual(analyzed["weight_ess"], expected_ess)
                self.assertTrue(math.isfinite(analyzed["weight_ess"]))


class ResponseAndStationarityTests(unittest.TestCase):
    def test_zero_mean_stays_unranked_and_has_no_log_response(self) -> None:
        response = ranking.paired_response(
            {seed: (0.0, 0.0) for seed in ranking.FIXED_SEEDS},
            {seed: (0.0, 1.0) for seed in ranking.FIXED_SEEDS},
            nominal_rankable=False,
            perturbed_rankable=True,
        )
        self.assertEqual(response["status"], "unranked_censored")
        self.assertIsNone(response["delta_log10"])
        self.assertEqual(len(response["same_seed_linear_deltas"]), 8)

    def test_paired_means_are_aggregated_in_log_domain(self) -> None:
        nominal = {seed: (1000.0, 1e-300) for seed in ranking.FIXED_SEEDS}
        changed = {seed: (1000.0, 2e-300) for seed in ranking.FIXED_SEEDS}
        response = ranking.paired_response(
            nominal,
            changed,
            nominal_rankable=True,
            perturbed_rankable=True,
        )
        self.assertEqual(response["status"], "rankable")
        self.assertAlmostEqual(response["delta_log10"], math.log10(2.0))
        self.assertTrue(math.isfinite(response["mean_linear_delta"]))

    def test_log_weight_1000_never_overflows_paired_response(self) -> None:
        nominal = {seed: (1000.0, 1.0) for seed in ranking.FIXED_SEEDS}
        changed = {seed: (1000.0, 2.0) for seed in ranking.FIXED_SEEDS}
        response = ranking.paired_response(
            nominal,
            changed,
            nominal_rankable=True,
            perturbed_rankable=True,
        )
        independent = verifier._paired(nominal, changed, True, True, 192837)
        for result in (response, independent):
            self.assertEqual(result["status"], "rankable")
            self.assertAlmostEqual(result["delta_log10"], math.log10(2.0))
            self.assertIsNone(result["mean_linear_delta"])
            self.assertEqual(result["mean_delta_sign"], 1)
            self.assertTrue(math.isfinite(result["log_abs_mean_delta_mol_m2_s"]))

    def test_extreme_log_weight_is_typed_censor_before_linear_materialization(
        self,
    ) -> None:
        nominal = {seed: (1000.0, 1.0) for seed in ranking.FIXED_SEEDS}
        changed = {seed: (1000.0, 0.0) for seed in ranking.FIXED_SEEDS}
        for result in (
            ranking.paired_response(
                nominal,
                changed,
                nominal_rankable=False,
                perturbed_rankable=False,
            ),
            verifier._paired(nominal, changed, False, False, 192837),
        ):
            self.assertEqual(result["status"], "unranked_censored")
            self.assertIsNone(result["mean_linear_delta"])
            self.assertEqual(result["mean_delta_sign"], -1)

    def test_incomplete_and_nonstationary_results_are_typed_censors(self) -> None:
        estimate = ranking.estimate_species([0.0] * 8, [1.0] * 8)
        incomplete = ranking._apply_run_and_stationarity_gates(
            estimate, runs_complete=False, stationarity_passed=True
        )
        self.assertEqual(incomplete["status"], "censored_incomplete_replicas")
        self.assertIsNone(incomplete["mean"])
        nonstationary = ranking._apply_run_and_stationarity_gates(
            estimate, runs_complete=True, stationarity_passed=False
        )
        self.assertEqual(
            nonstationary["status"],
            "censored_nonstationary_state_distribution",
        )
        self.assertFalse(nonstationary["rankable"])

    def test_stationarity_requires_every_kind_state_distribution(self) -> None:
        stable = [
            {"A.x": 5.0, "A.y": 5.0, "B.z": 10.0},
            {"A.x": 5.1, "A.y": 4.9, "B.z": 10.0},
            {"A.x": 5.0, "A.y": 5.0, "B.z": 10.0},
            {"A.x": 5.1, "A.y": 4.9, "B.z": 10.0},
            {"A.x": 5.0, "A.y": 5.0, "B.z": 10.0},
        ]
        self.assertTrue(ranking.population_stationarity(stable)["passed"])
        unstable = [dict(row) for row in stable]
        unstable[-1]["A.x"] = 8.0
        unstable[-1]["A.y"] = 2.0
        result = ranking.population_stationarity(unstable)
        self.assertFalse(result["passed"])
        self.assertEqual(result["kinds"]["A"]["status"], "nonstationary")

    def test_irrelevance_needs_complete_interval_inside_equivalence_band(self) -> None:
        self.assertTrue(ranking.interval_demonstrates_irrelevance([-0.05, 0.04]))
        self.assertFalse(ranking.interval_demonstrates_irrelevance([-0.2, 0.01]))
        self.assertFalse(ranking.interval_demonstrates_irrelevance(None))


class AtomicWriteTests(unittest.TestCase):
    def test_nonoverwrite_atomic_json_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "receipt.json"
            ranking.write_json_atomic_new(path, {"status": "complete"})
            with self.assertRaises(FileExistsError):
                ranking.write_json_atomic_new(path, {"status": "replaced"})
            self.assertIn("complete", path.read_text(encoding="utf-8"))


class IndependentRawReplayTests(unittest.TestCase):
    def _fixture(self, root: Path) -> Path:
        ranking.prepare_campaign(root, smoke=True)
        return _build_synthetic_package(root, "nominal", ranking.FIXED_SEEDS[0])

    def test_readsorption_and_duplicate_release_consumption_count_once(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._fixture(root)
            audited = verifier._audit_replica(root, "nominal", ranking.FIXED_SEEDS[0])
            self.assertEqual(audited["release_counts"], {"Al": 0, "Si": 1})
            self.assertEqual(audited["event_count"], 3)

    def test_zero_hazard_fired_event_is_impossible_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            package = self._fixture(root)
            rows = [
                json.loads(line)
                for line in (package / "likelihood.jsonl").read_text().splitlines()
            ]
            fired = next(
                row for row in rows[1:] if row["fired_reaction_id"] is not None
            )
            self.assertGreater(fired["physical_total_rate_s-1"], 0)
            self.assertGreater(fired["biased_total_rate_s-1"], 0)
            fired["physical_total_rate_s-1"] = 0.0
            fired["biased_total_rate_s-1"] = 0.0
            _write_jsonl(package / "likelihood.jsonl", rows)
            _rehash_package(root, "nominal", ranking.FIXED_SEEDS[0])
            manifest = json.loads((root / "manifest.json").read_text())
            with self.assertRaisesRegex(ValueError, "require positive"):
                ranking._verify_package(
                    root, "nominal", ranking.FIXED_SEEDS[0], manifest
                )
            with self.assertRaisesRegex(ValueError, "require positive"):
                verifier._audit_replica(root, "nominal", ranking.FIXED_SEEDS[0])

    def test_likelihood_event_pgif_timestamp_lineage_and_checkpoint_tampering(
        self,
    ) -> None:
        mutations = {
            "likelihood": (
                "likelihood.jsonl",
                lambda rows: rows[1].__setitem__("log_likelihood_increment", 0.0),
            ),
            "event_state": (
                "events.jsonl",
                lambda rows: rows[1]["changes"][0].__setitem__("old_state_id", 999),
            ),
            "timestamp": (
                "events.jsonl",
                lambda rows: rows[1].__setitem__("time_s", 9e-13),
            ),
            "lineage": (
                "events.jsonl",
                lambda rows: rows[1].__setitem__("lattice_release_kind", None),
            ),
            "checkpoint": (
                "checkpoints.jsonl",
                lambda rows: rows[1]["state_counts"].__setitem__("Si.empty", 999),
            ),
            "checkpoint_target": (
                "checkpoints.jsonl",
                lambda rows: rows[1].__setitem__("target_time_s", 7e-13),
            ),
        }
        for label, (filename, mutate) in mutations.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                package = self._fixture(root)
                rows = [
                    json.loads(line)
                    for line in (package / filename).read_text().splitlines()
                ]
                mutate(rows)
                _write_jsonl(package / filename, rows)
                _rehash_package(root, "nominal", ranking.FIXED_SEEDS[0])
                with self.assertRaises((ValueError, IndexError)):
                    verifier._audit_replica(root, "nominal", ranking.FIXED_SEEDS[0])

    def test_initial_and_final_pgif_tampering_is_detected_after_rehash(self) -> None:
        for filename in ("initial.pgif.json", "final.pgif.json"):
            with (
                self.subTest(filename=filename),
                tempfile.TemporaryDirectory() as directory,
            ):
                root = Path(directory)
                package = self._fixture(root)
                pgif = json.loads((package / filename).read_text())
                pgif["nodes"]["columns"]["state"]["data"][0] = 999
                (package / filename).write_text(json.dumps(pgif), encoding="utf-8")
                _rehash_package(root, "nominal", ranking.FIXED_SEEDS[0])
                with self.assertRaises((ValueError, IndexError)):
                    verifier._audit_replica(root, "nominal", ranking.FIXED_SEEDS[0])

    def test_edge_type_coordinate_metadata_and_timestamp_tampering_after_rehash(
        self,
    ) -> None:
        mutations = {
            "edge": lambda pgif: pgif["edges"]["dst"].__setitem__(
                0, pgif["edges"]["dst"][0] + 1
            ),
            "type": lambda pgif: pgif["nodes"]["columns"]["type"]["data"].__setitem__(
                0, 0
            ),
            "coordinate": lambda pgif: pgif["nodes"]["columns"]["x"][
                "data"
            ].__setitem__(0, pgif["nodes"]["columns"]["x"]["data"][0] + 0.25),
            "graph_metadata": lambda pgif: pgif["meta"].__setitem__("kind", "tampered"),
            "timestamp": lambda pgif: pgif["meta"]["petra"].__setitem__("time", 1e-15),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                package = self._fixture(root)
                path = package / "initial.pgif.json"
                pgif = json.loads(path.read_text())
                mutate(pgif)
                path.write_text(json.dumps(pgif), encoding="utf-8")
                _rehash_package(root, "nominal", ranking.FIXED_SEEDS[0])
                with self.assertRaises((ValueError, IndexError)):
                    verifier._audit_replica(root, "nominal", ranking.FIXED_SEEDS[0])

    def test_runner_compares_full_snapshot_to_canonical_deck_pgif(self) -> None:
        source = (
            ranking.PETRA / "crates/petra-deck/examples/a9b_importance_sampling.rs"
        ).read_text(encoding="utf-8")
        self.assertIn("petra_io::snapshot_json(deck, engine)", source)
        self.assertIn("if snapshot != &canonical", source)
        self.assertIn("metadata, nodes, coordinates, or edges", source)

    def test_exact_command_receipt_tampering_is_detected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._fixture(root)
            receipt_path = (
                root / "receipts" / "nominal" / f"seed-{ranking.FIXED_SEEDS[0]}.json"
            )
            receipt = json.loads(receipt_path.read_text())
            receipt["command"][-1] = "999"
            receipt_path.write_bytes(_json_bytes(receipt))
            with self.assertRaisesRegex(ValueError, "binding"):
                verifier._audit_replica(root, "nominal", ranking.FIXED_SEEDS[0])

    def test_reservoir_and_path_hash_bindings_are_independently_checked(self) -> None:
        for input_name in ("reservoir_contract", "path_report", "ph4_evidence"):
            with (
                self.subTest(input_name=input_name),
                tempfile.TemporaryDirectory() as directory,
            ):
                root = Path(directory)
                manifest = ranking.prepare_campaign(root, smoke=True)
                path = root / manifest["inputs"][input_name]
                path.write_bytes(path.read_bytes() + b"\n")
                with self.assertRaisesRegex(ValueError, "hash mismatch"):
                    verifier._verify_and_regenerate_decks(root, manifest)


class IndependentByteRegenerationTests(unittest.TestCase):
    def test_verifier_rejects_preregistered_contract_tampering(self) -> None:
        mutations = {
            "seeds": [999],
            "top_k_preregistered": 2,
            "equivalence_band_log10": [-1.0, 1.0],
            "weight_ess_fraction_min": 0.25,
            "species_contribution_ess_min": 1.0,
            "state_population_stationarity": {"required": "nothing"},
            "run_contract": {
                "horizon_s": 1e-12,
                "checkpoints": 3,
                "bias_factor": 1e6,
                "max_events": 999,
            },
        }
        with tempfile.TemporaryDirectory() as directory:
            manifest = ranking.prepare_campaign(Path(directory), smoke=True)
            for key, value in mutations.items():
                with self.subTest(key=key):
                    tampered = copy.deepcopy(manifest)
                    tampered[key] = value
                    with self.assertRaisesRegex(ValueError, f"mismatch in {key}"):
                        verifier._validate_manifest_contract(tampered)

    def test_independent_verifier_regenerates_all_derived_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _synthetic_campaign(root)
            analysis = ranking.analyze_campaign(root)
            receipt = verifier.verify_campaign(root, write_verification=True)
            self.assertTrue(receipt["derived_bytes_equal"])
            self.assertEqual(receipt["scenario_decks_regenerated"], 29)
            self.assertEqual(receipt["raw_replicas_replayed"], 29)
            self.assertFalse(analysis["artifact_verified"])
            self.assertFalse(analysis["ordinal_ranking_complete"])
            self.assertEqual(
                analysis["verification_status"], "pending_independent_verification"
            )
            self.assertEqual(
                analysis["scientific_acceptance"],
                "unverified_pending_independent_verification",
            )
            self.assertEqual(
                analysis["preliminary_conclusions"]["status"],
                "unverified_preliminary",
            )
            finalized = receipt["finalized_analysis"]
            self.assertTrue(finalized["artifact_verified"])
            self.assertEqual(finalized["verification_status"], "independently_verified")
            self.assertFalse(finalized["ordinal_ranking_complete"])
            self.assertEqual(finalized["ordinal_ranking_status"], "not_evaluated_smoke")
            self.assertEqual(finalized["scientific_acceptance"], "not_evaluated_smoke")
            evolution = analysis["scenarios"]["nominal"][
                "state_population_evolution_unnormalized_is"
            ]
            self.assertEqual([row["time_s"] for row in evolution], [0.0, 5e-13, 1e-12])
            initial_counts = json.loads(
                (root / "runs/nominal/seed-90401/checkpoints.jsonl")
                .read_text()
                .splitlines()[0]
            )["state_counts"]
            self.assertEqual(evolution[0]["state_counts"], initial_counts)
            self.assertIn("Si.empty", evolution[-1]["state_counts"])
            self.assertTrue((root / "verification.json").is_file())

            reanalysis = ranking.analyze_campaign(root)
            self.assertFalse(reanalysis["artifact_verified"])
            self.assertFalse((root / "verification.json").exists())
            verifier.verify_campaign(root, write_verification=True)
            self.assertTrue((root / "verification.json").is_file())
            (root / "analysis.json").write_bytes(b"{}\n")
            with self.assertRaisesRegex(ValueError, "analysis bytes differ"):
                verifier.verify_campaign(root, write_verification=True)
            self.assertFalse((root / "verification.json").exists())
            self.assertFalse(analysis["artifact_verified"])

    def test_full_232_run_bundle_is_exactly_independently_regenerated(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _synthetic_campaign(root, smoke=False)
            analysis = ranking.analyze_campaign(root)
            manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))

            verifier._verify_and_regenerate_decks(root, manifest)
            independently_derived = verifier.independently_derive(root, manifest)
            independently_regenerated = verifier._pending_analysis(
                independently_derived
            )

            self.assertEqual(analysis["run_count"], 232)
            self.assertEqual(independently_derived["run_count"], 232)
            self.assertEqual(
                (root / "analysis.json").read_bytes(),
                _json_bytes(independently_regenerated),
            )
            evolution = analysis["scenarios"]["nominal"][
                "state_population_evolution_unnormalized_is"
            ]
            checkpoints = []
            final_logs = []
            for seed in ranking.FIXED_SEEDS:
                package = root / "runs" / "nominal" / f"seed-{seed}"
                checkpoint_rows = [
                    json.loads(line)
                    for line in (package / "checkpoints.jsonl").read_text().splitlines()
                ]
                checkpoints.append(checkpoint_rows[1])
                metadata = json.loads((package / "metadata.json").read_text())
                final_logs.append(metadata["final_log_likelihood"])
            self.assertGreater(
                len({row["cumulative_log_likelihood"] for row in checkpoints}), 1
            )
            state = "Si.empty"
            expected = math.fsum(
                math.exp(row["cumulative_log_likelihood"]) * row["state_counts"][state]
                for row in checkpoints
            ) / len(checkpoints)
            old_final_weight_result = math.fsum(
                math.exp(log_weight) * row["state_counts"][state]
                for log_weight, row in zip(final_logs, checkpoints, strict=True)
            ) / len(checkpoints)
            self.assertAlmostEqual(evolution[1]["state_counts"][state], expected)
            self.assertNotAlmostEqual(expected, old_final_weight_result)

    def test_verifier_rejects_tier_and_temperature_manifest_tampering(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _synthetic_campaign(root)
            ranking.analyze_campaign(root)
            original = json.loads((root / "manifest.json").read_text())
            for key, value in (
                ("survey_tier_platform_test", False),
                ("not_production", False),
                ("temperature_k", 999.0),
                ("production_kinetics", True),
            ):
                with self.subTest(key=key):
                    tampered = copy.deepcopy(original)
                    tampered[key] = value
                    ranking.write_json_atomic(root / "manifest.json", tampered)
                    expected_error = (
                        "key set"
                        if key == "production_kinetics"
                        else f"mismatch in {key}"
                    )
                    with self.assertRaisesRegex(ValueError, expected_error):
                        verifier.verify_campaign(root, write_verification=False)
            ranking.write_json_atomic(root / "manifest.json", original)


class LauncherContractTests(unittest.TestCase):
    def test_public_launcher_self_dispatches_to_fixed_bounded_systemd_unit(
        self,
    ) -> None:
        launcher = MODULE_PATH.with_name("launch_a9b_sensitivity_campaign.sh")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fake_bin = root / "bin"
            fake_bin.mkdir()
            capture = root / "args.json"
            fake = fake_bin / "systemd-run"
            fake.write_text(
                "#!/usr/bin/env python3\n"
                "import json, os, pathlib, sys\n"
                "pathlib.Path(os.environ['A9B_TEST_CAPTURE']).write_text(json.dumps(sys.argv[1:]))\n",
                encoding="utf-8",
            )
            fake.chmod(0o755)
            environment = {
                **os.environ,
                "PATH": f"{fake_bin}:{os.environ['PATH']}",
                "A9B_TEST_CAPTURE": str(capture),
            }
            completed = subprocess.run(
                [f"./{launcher.name}", str(root / "campaign"), "3"],
                check=False,
                capture_output=True,
                text=True,
                cwd=launcher.parent,
                env=environment,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            argv = json.loads(capture.read_text())
            self.assertIn("--unit=a9b-sensitivity-campaign", argv)
            self.assertIn("--property=RuntimeMaxSec=86400", argv)
            launcher_index = argv.index(str(launcher))
            self.assertEqual(
                argv[launcher_index - 3 : launcher_index],
                ["/usr/bin/nice", "-n", "10"],
            )
            self.assertFalse(any("A9B_SYSTEMD_INNER" in value for value in argv))
            runtime = int(
                next(
                    value.rsplit("=", 1)[1]
                    for value in argv
                    if "RuntimeMaxSec=" in value
                )
            )
            self.assertLessEqual(runtime, 86400)
            source = launcher.read_text(encoding="utf-8")
            self.assertLess(
                source.index("trap finish EXIT"), source.index("/proc/loadavg")
            )
            self.assertIn('"a9b-campaign-launch-receipt-v2"', source)
            self.assertIn('"runtime_max_sec"', source)

    def test_forged_private_launcher_environment_cannot_reach_inner_workload(
        self,
    ) -> None:
        launcher = MODULE_PATH.with_name("launch_a9b_sensitivity_campaign.sh")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fake_bin = root / "bin"
            fake_bin.mkdir()
            capture = root / "args.json"
            fake = fake_bin / "systemd-run"
            fake.write_text(
                "#!/usr/bin/env python3\n"
                "import json, os, pathlib, sys\n"
                "pathlib.Path(os.environ['A9B_TEST_CAPTURE']).write_text(json.dumps(sys.argv[1:]))\n",
                encoding="utf-8",
            )
            fake.chmod(0o755)
            completed = subprocess.run(
                [str(launcher), str(root / "campaign")],
                check=False,
                capture_output=True,
                text=True,
                env={
                    **os.environ,
                    "PATH": f"{fake_bin}:{os.environ['PATH']}",
                    "A9B_TEST_CAPTURE": str(capture),
                    "A9B_SYSTEMD_INNER": "1",
                    "A9B_SYSTEMD_UNIT": "a9b-sensitivity-campaign",
                    "A9B_NICED": "1",
                    "INVOCATION_ID": "forged",
                },
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue(capture.is_file())
            self.assertFalse((root / "campaign-operator").exists())

    def test_process_niceness_is_read_from_proc_stat_and_not_environment(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            stat = Path(directory) / "stat"
            fields_from_state_through_nice = [
                "S",
                *[str(value) for value in range(4, 18)],
                "20",
                "13",
            ]
            stat.write_text(
                f"123 (worker name with spaces) {' '.join(fields_from_state_through_nice)}\n",
                encoding="utf-8",
            )
            self.assertEqual(launch_guard.process_niceness(stat), 13)
            self.assertEqual(
                launch_guard.process_niceness(), os.getpriority(os.PRIO_PROCESS, 0)
            )
            with self.assertRaisesRegex(ValueError, "niceness"):
                launch_guard.require_minimum_niceness(9, 10)

        launcher = MODULE_PATH.with_name("launch_a9b_sensitivity_campaign.sh")
        source = launcher.read_text(encoding="utf-8")
        self.assertNotIn("${A9B_NICED", source)
        self.assertIn('int(os.environ["A9B_LIVE_NICENESS"])', source)

    def test_live_systemd_proof_validation_is_exact_and_bounded(self) -> None:
        expected = "a9b-sensitivity-campaign.service"
        cgroup = "/user.slice/user-1000.slice/user@1000.service/app.slice/a9b.service"
        proof = {
            "unit_id": expected,
            "active_state": "active",
            "control_group": cgroup,
            "runtime_max_usec": 86_400_000_000,
            "process_cgroup": cgroup,
        }
        self.assertEqual(
            launch_guard.validate_live_proof(proof, expected, 86_400),
            86_400_000_000,
        )
        for field, forged in (
            ("unit_id", "attacker.service"),
            ("active_state", "inactive"),
            ("control_group", "/attacker.service"),
            ("runtime_max_usec", 86_400_000_001),
        ):
            with self.subTest(field=field):
                altered = {**proof, field: forged}
                with self.assertRaises(ValueError):
                    launch_guard.validate_live_proof(altered, expected, 86_400)


if __name__ == "__main__":
    unittest.main()
