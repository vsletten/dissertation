from __future__ import annotations

import csv
import importlib.util
import json
import sys
import tempfile
import tomllib
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1]
ANALYSIS_PATH = SCRIPTS / "muscovite_full_analysis.py"
BUILD_PATH = SCRIPTS / "build_muscovite_full_deck.py"
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


if __name__ == "__main__":
    unittest.main()
