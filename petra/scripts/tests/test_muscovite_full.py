from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import tomllib

SCRIPTS = Path(__file__).resolve().parents[1]
ANALYSIS_PATH = SCRIPTS / "muscovite_full_analysis.py"
BUILD_PATH = SCRIPTS / "build_muscovite_full_deck.py"
SWEEP_PATH = SCRIPTS / "muscovite_grain_size_sweep.py"
DECK_PATH = SCRIPTS.parents[0] / "decks" / "muscovite-full-mechanism.toml"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


analysis = load_module("muscovite_full_analysis", ANALYSIS_PATH)
builder = load_module("build_muscovite_full_deck", BUILD_PATH)


class FullMechanismDeckTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.text = DECK_PATH.read_text(encoding="utf-8")
        cls.deck = tomllib.loads(cls.text)

    def test_tracked_deck_is_reproducibly_generated(self) -> None:
        self.assertEqual(builder.render_deck(), self.text)

    def test_sized_deck_updates_volume_surface_and_name_deterministically(self) -> None:
        rendered = builder.render_deck(dims=(6, 8, 10))
        sized = tomllib.loads(rendered)
        self.assertEqual(sized["structure"]["lattice"]["dims"], [6, 8, 10])
        self.assertEqual(
            sized["structure"]["init"][-1]["region"],
            {"axis": 2, "min": 9, "max": 9},
        )
        self.assertEqual(sized["deck"]["name"], "muscovite-full-6x8x10")
        self.assertEqual(rendered, builder.render_deck(dims=(6, 8, 10)))

    def test_schedule_and_reservoirs_are_explicit(self) -> None:
        schedule = self.deck["execution"]["schedule"]
        self.assertGreaterEqual(len(schedule), 5)
        self.assertEqual(schedule[0]["temperature"], 773.15)
        self.assertEqual(schedule[-1]["temperature"], 973.15)
        kinds = {kind["name"]: kind for kind in self.deck["structure"]["kinds"]}
        self.assertIn("Defect_zone", kinds)
        self.assertIn("Interface", kinds)
        self.assertIn("Octahedral_trap", kinds)
        self.assertIn(
            "Ar39_trapped", {s["name"] for s in kinds["Octahedral_trap"]["states"]}
        )
        species = {item["name"] for item in self.deck["structure"]["species"]}
        self.assertTrue({"Ar40", "Ar39", "Ar36"} <= species)

    def test_delamination_is_local_dehydroxylation_driven(self) -> None:
        rules = {rule["name"]: rule for rule in self.deck["dynamics"]["rules"]}
        delamination = rules["delaminate_interface"]
        self.assertEqual(delamination["center"]["kind"], "Interface")
        self.assertTrue(
            any(
                guard.get("kind") == "OH_acceptor"
                and guard.get("state") == ["O_residual"]
                and guard.get("label") == "delam_driver"
                for guard in delamination["guards"]
            )
        )
        self.assertNotIn("modifier", delamination)

    def test_species_specific_release_and_zone_specific_hops_exist(self) -> None:
        names = {rule["name"] for rule in self.deck["dynamics"]["rules"]}
        for species in ("Ar40", "Ar39", "Ar36"):
            self.assertTrue(
                any(name.startswith(f"release_{species}_") for name in names)
            )
            self.assertTrue(
                any(name.startswith(f"hop_pristine_{species}_") for name in names)
            )
            self.assertTrue(
                any(name.startswith(f"hop_extended_{species}_") for name in names)
            )


class GrainSizeSweepDriverTests(unittest.TestCase):
    def test_size_ladder_and_event_receipt_parsing_are_fail_closed(self) -> None:
        sweep = load_module("muscovite_grain_size_sweep", SWEEP_PATH)
        self.assertEqual(sweep.parse_sizes("4x4x6,8x8x12"), ((4, 4, 6), (8, 8, 12)))
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            sweep.parse_sizes("8x8x12,4x4x6")
        with tempfile.TemporaryDirectory() as directory:
            ensemble = Path(directory) / "ensemble.csv"
            ensemble.write_text(
                "seed,steps,time,state\n11,10,1.0,0\n12,15,1.0,0\n",
                encoding="utf-8",
            )
            self.assertEqual(sweep.total_events(ensemble), 25)
            receipt = analysis.SizeRunReceipt(
                dims=(4, 4, 6),
                sites=768,
                replicas=2,
                elapsed_seconds=1.0,
                total_events=25,
                replay_verified=True,
            )
            campaign = Path(directory) / "campaign.json"
            sweep.write_campaign_receipts(campaign, [receipt])
            self.assertEqual(
                json.loads(campaign.read_text(encoding="utf-8"))[0]["sites"], 768
            )
            self.assertFalse(campaign.with_suffix(".json.tmp").exists())
            matching = analysis.EnsembleAnalysis(
                dims=(4, 4, 6), sites=768, replicas=2, steps=()
            )
            sweep.validate_result_receipt(matching, receipt)
            with self.assertRaisesRegex(ValueError, "replica count"):
                sweep.validate_result_receipt(replace(matching, replicas=1), receipt)
            malformed = campaign.with_name("malformed.json")
            malformed.write_text(
                json.dumps(
                    [
                        {
                            **json.loads(campaign.read_text(encoding="utf-8"))[0],
                            "replay_verified": "false",
                        }
                    ]
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "replay_verified"):
                sweep.load_campaign_receipts(malformed)
            ensemble.write_text(
                "seed,steps,time,state\n11,10,1.0,0\n11,15,1.0,0\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "duplicate seed"):
                sweep.total_events(ensemble, expected_replicas=2)
            duplicate = campaign.with_name("duplicate.json")
            payload = json.loads(campaign.read_text(encoding="utf-8"))
            duplicate.write_text(json.dumps(payload + payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "strictly increasing"):
                sweep.load_campaign_receipts(duplicate)

    def test_nonnegative_integer_preserves_identity_beyond_float_precision(self) -> None:
        sweep = load_module("muscovite_grain_size_sweep", SWEEP_PATH)
        huge = 2**53 + 1
        self.assertEqual(sweep._nonnegative_integer(str(huge), "seed"), huge)
        self.assertEqual(sweep._nonnegative_integer("0", "steps"), 0)
        with self.assertRaisesRegex(ValueError, "non-negative integer"):
            sweep._nonnegative_integer("-1", "seed")
        with self.assertRaisesRegex(ValueError, "non-negative integer"):
            sweep._nonnegative_integer("1.0", "seed")
        with self.assertRaisesRegex(ValueError, "non-negative integer"):
            sweep._nonnegative_integer("1e2", "steps")
        with tempfile.TemporaryDirectory() as directory:
            ensemble = Path(directory) / "ensemble.csv"
            ensemble.write_text(
                f"seed,steps,time,state\n{huge},10,1.0,0\n{huge + 1},15,1.0,0\n",
                encoding="utf-8",
            )
            self.assertEqual(sweep.total_events(ensemble), 25)
            rows = sweep.ensemble_rows(ensemble)
            self.assertEqual(
                [sweep._nonnegative_integer(row["seed"], "seed") for row in rows],
                [huge, huge + 1],
            )


class FullMechanismAnalysisTests(unittest.TestCase):
    def _write_fixture(self, root: Path, increase: bool = False) -> tuple[Path, Path]:
        deck = root / "deck.toml"
        deck.write_text(
            """[execution]\n[[execution.schedule]]\ntemperature=773.15\nduration=10.0\n[[execution.schedule]]\ntemperature=873.15\nduration=20.0\n""",
            encoding="utf-8",
        )
        populations = root / "populations.csv"
        rows = [
            [0, 0.0, 40, 10, 30, 20],
            [5, 10.0, 30, 8, 24, 15],
            [9, 30.0, 20, 4, 15, 8],
        ]
        if increase:
            rows[1][2] = 41
        with populations.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                [
                    "step",
                    "time",
                    "K_gallery_A.Ar40",
                    "K_gallery_A.Ar39",
                    "Octahedral_trap.Ar39_trapped",
                    "K_gallery_A.Ar36",
                ]
            )
            writer.writerows(rows)
        return deck, populations

    def test_species_release_age_and_ratio_are_segment_resolved(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            deck, populations = self._write_fixture(root)
            summary = analysis.analyze(deck, populations, j_factor=0.01)
        self.assertEqual([row.released_ar40 for row in summary.steps], [10, 10])
        self.assertEqual([row.released_ar39 for row in summary.steps], [8, 13])
        self.assertEqual([row.released_ar36 for row in summary.steps], [5, 7])
        self.assertAlmostEqual(summary.steps[0].ar36_ar40, 0.5)
        self.assertGreater(summary.steps[0].apparent_age_ma, 0.0)
        self.assertEqual(summary.initial["Ar39"], 40)

    def test_inventory_increase_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            deck, populations = self._write_fixture(root, increase=True)
            with self.assertRaisesRegex(ValueError, "inventory increased"):
                analysis.analyze(deck, populations, j_factor=0.01)

    def test_uncomputed_barrier_brackets_are_monotonic_and_labeled_proxy(self) -> None:
        rows = analysis.sensitivity_brackets(
            mechanism="delamination",
            temperatures_k=[773.15, 973.15],
            nominal_kcal=58.0,
            delta_kcal=5.0,
        )
        self.assertEqual({row.mechanism for row in rows}, {"delamination"})
        self.assertEqual(
            {row.provenance for row in rows},
            {"proxy ±5 kcal/mol; not computed kinetics"},
        )
        for temperature in (773.15, 973.15):
            values = [
                row.rate_multiplier for row in rows if row.temperature_k == temperature
            ]
            self.assertGreater(values[0], values[1])
            self.assertGreater(values[1], values[2])

    def test_json_output_is_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            deck, populations = self._write_fixture(root)
            summary = analysis.analyze(deck, populations, j_factor=0.01)
            first = analysis.summary_json(summary)
            second = analysis.summary_json(summary)
        self.assertEqual(first, second)
        self.assertEqual(json.loads(first)["initial"]["Ar40"], 40)

    def test_ensemble_analysis_preserves_replica_distributions_and_step_bands(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            deck = root / "deck.toml"
            deck.write_text(
                """[structure]
[[structure.kinds]]
name = "Gallery"
initial = "K"
[[structure.kinds.states]]
name = "K"
occupant = "K"
[[structure.kinds.states]]
name = "Ar40"
occupant = "Ar40"
[[structure.kinds.states]]
name = "Ar39"
occupant = "Ar39"
[[structure.kinds.states]]
name = "Ar36"
occupant = "Ar36"
[[structure.kinds]]
name = "Trap"
initial = "vacant"
[[structure.kinds.states]]
name = "vacant"
occupant = "vacant"
[[structure.kinds.states]]
name = "Ar39_trapped"
occupant = "Ar39"
[execution]
[[execution.schedule]]
temperature = 773.15
duration = 10.0
[[execution.schedule]]
temperature = 873.15
duration = 20.0
""",
                encoding="utf-8",
            )
            observables = root / "observables.csv"
            with observables.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(
                    ["replica", "seed", "step", "time", "kind", "index", "value"]
                )
                # Flattened state order: K, Ar40, Ar39, Ar36, vacant, Ar39_trapped.
                for replica, seed, scale, ar36_scale in (
                    (0, 11, 1, 0),
                    (1, 12, 2, 2),
                ):
                    for step, time_s, counts in (
                        (
                            0,
                            0.0,
                            [0, 40 * scale, 10 * scale, 20 * ar36_scale, 0, 30 * scale],
                        ),
                        (
                            5,
                            10.0,
                            [0, 30 * scale, 8 * scale, 15 * ar36_scale, 0, 24 * scale],
                        ),
                        (
                            9,
                            30.0,
                            [0, 20 * scale, 4 * scale, 8 * ar36_scale, 0, 15 * scale],
                        ),
                    ):
                        for index, value in enumerate(counts):
                            writer.writerow(
                                [
                                    replica,
                                    seed,
                                    step,
                                    time_s,
                                    "state_counts",
                                    index,
                                    value,
                                ]
                            )
            result = analysis.analyze_ensemble(
                deck, observables, dims=(4, 4, 6), j_factor=0.01
            )
            lines = observables.read_text(encoding="utf-8").splitlines()
            observables.write_text("\n".join(lines[:-1]) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "incomplete state-count sample"):
                analysis.analyze_ensemble(
                    deck, observables, dims=(4, 4, 6), j_factor=0.01
                )
            broken = lines.copy()
            fields = broken[-1].split(",")
            fields[-1] = "1.5"
            broken[-1] = ",".join(fields)
            observables.write_text("\n".join(broken) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "non-negative integer"):
                analysis.analyze_ensemble(
                    deck, observables, dims=(4, 4, 6), j_factor=0.01
                )
        self.assertEqual(result.replicas, 2)
        self.assertEqual(result.sites, 768)
        self.assertEqual(result.steps[0].released_ar40.values, (10.0, 20.0))
        self.assertEqual(result.steps[0].released_ar40.mean, 15.0)
        self.assertEqual(result.steps[1].released_ar39.values, (13.0, 26.0))
        self.assertEqual(result.steps[0].released_ar36.values, (0.0, 10.0))
        self.assertEqual(result.steps[0].cumulative_ar36_fraction.values, (0.25,))
        self.assertTrue(
            all(value > 0.0 for value in result.steps[0].apparent_age_ma.values)
        )
        self.assertEqual(result.steps[0].ar36_ar40.values, (0.0, 0.5))
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)
            analysis.write_ensemble_products(
                [result],
                [
                    analysis.SizeRunReceipt(
                        dims=(4, 4, 6),
                        sites=768,
                        replicas=2,
                        elapsed_seconds=2.0,
                        total_events=100,
                        replay_verified=True,
                    )
                ],
                out,
            )
            bands = (out / "spectrum-bands.csv").read_text(encoding="utf-8")
            scaling = (out / "size-scaling.csv").read_text(encoding="utf-8")
        self.assertIn("apparent_age_ma_ci95_low", bands.splitlines()[0])
        self.assertIn("apparent_age_ma_defined_replicas", bands.splitlines()[0])
        self.assertIn("10.0;20.0", bands)
        self.assertIn("events_per_second", scaling.splitlines()[0])
        self.assertIn(",50.0,", scaling)
        larger = analysis.EnsembleAnalysis(
            dims=(6, 6, 9), sites=2592, replicas=2, steps=result.steps
        )
        stability = analysis.assess_stability([result, larger])
        self.assertEqual(stability.stabilized_at_sites, 2592)
        self.assertEqual(stability.max_release_fraction_delta, 0.0)
        self.assertEqual(stability.max_age_relative_delta, 0.0)
        self.assertEqual(stability.max_age_defined_fraction_delta, 0.0)
        partial_age = replace(
            larger,
            steps=(
                replace(
                    larger.steps[0],
                    apparent_age_ma=analysis.distribution_band(
                        [larger.steps[0].apparent_age_ma.mean], seed=1
                    ),
                ),
                *larger.steps[1:],
            ),
        )
        incomplete_support = analysis.assess_stability([result, partial_age])
        self.assertFalse(incomplete_support.comparisons[0].stable)
        self.assertEqual(
            incomplete_support.comparisons[0].max_age_defined_fraction_delta, 0.5
        )
        missing_age = replace(
            larger,
            steps=(
                replace(
                    larger.steps[0],
                    apparent_age_ma=analysis.distribution_band([], seed=1),
                ),
                *larger.steps[1:],
            ),
        )
        not_comparable = analysis.assess_stability([result, missing_age])
        self.assertIsNone(not_comparable.comparisons[0].max_age_relative_delta)
        self.assertEqual(
            not_comparable.comparisons[0].max_age_defined_fraction_delta, 1.0
        )
        self.assertIsNone(not_comparable.stabilized_at_sites)
        mismatched_schedule = replace(
            larger,
            steps=(
                replace(larger.steps[0], temperature_k=999.0),
                *larger.steps[1:],
            ),
        )
        with self.assertRaisesRegex(ValueError, "same schedule"):
            analysis.assess_stability([result, mismatched_schedule])


if __name__ == "__main__":
    unittest.main()
