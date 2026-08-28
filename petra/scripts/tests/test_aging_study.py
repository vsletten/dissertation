from __future__ import annotations

import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1]
ANALYSIS_PATH = SCRIPTS / "aging_study.py"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


aging = load_module("aging_study", ANALYSIS_PATH)


class AgingStudyTests(unittest.TestCase):
    def test_render_deck_changes_only_name_and_defect_probability(self) -> None:
        template = """[deck]
name = \"kaolinite-aging-smoke\"
[[init]]
name = \"ground-surface-defects\"
probability = 0.25
[simulation]
steps = 600
"""
        rendered = aging.render_deck(template, 0.5)
        self.assertIn('name = "kaolinite-aging-density-0p500"', rendered)
        self.assertIn("probability = 0.500000", rendered)
        self.assertIn("steps = 600", rendered)
        with self.assertRaisesRegex(ValueError, "between zero and one"):
            aging.render_deck(template, 1.1)

    def _write_fixture(self, root: Path) -> Path:
        path = root / "observables.csv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["replica", "seed", "step", "time", "kind", "index", "value"])
            for replica, seed, scale in ((0, 10, 1.0), (1, 11, 2.0)):
                for step, time, defects, rate, ages in (
                    (0, 0.0, 4.0, 100.0 * scale, [0.0, 0.0]),
                    (140, 10.0 * scale, 0.0, 1.0 * scale, [8.0, 12.0]),
                ):
                    for index, value in enumerate((10.0, defects, 20.0)):
                        writer.writerow([replica, seed, step, time, "state_counts", index, value])
                    writer.writerow([replica, seed, step, time, "event_rates", 0, rate])
                    writer.writerow([replica, seed, step, time, "rate_spectra", 0, rate])
                    writer.writerow([replica, seed, step, time, "rate_spectra", 1, rate / 10.0])
                    for index, value in enumerate((50.0, 12.0, 12.0)):
                        writer.writerow([replica, seed, step, time, "surface_area", index, value])
                    for index, value in enumerate(ages):
                        writer.writerow([replica, seed, step, time, "exposure_age", index, value])
        return path

    def test_analyze_density_reports_ratio_bands_and_distributions(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = self._write_fixture(Path(directory))
            result = aging.analyze_density(path, density=0.25, aged_step=140)
        self.assertEqual(result.replicas, 2)
        self.assertEqual(result.initial_defects.mean, 4.0)
        self.assertEqual(result.drop_ratio.mean, 100.0)
        self.assertEqual(result.drop_ratio.ci95, (100.0, 100.0))
        self.assertEqual(result.geometric_normalized_drop.mean, 100.0)
        self.assertAlmostEqual(result.bet_normalized_drop.mean, 100.0)
        aged = result.curves[-1]
        self.assertEqual(aged.step, 140)
        self.assertEqual(aged.defects.mean, 0.0)
        self.assertEqual(aged.exposure_age_mean.mean, 10.0)
        self.assertGreater(aged.rate_log10_width.mean, 0.0)
        with tempfile.TemporaryDirectory() as directory:
            plot = Path(directory) / "plot.svg"
            aging.write_svg([result], plot)
            self.assertIn("2 replicas per density", plot.read_text(encoding="utf-8"))

    def test_analysis_fails_closed_if_aged_sample_retains_defects(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = self._write_fixture(Path(directory))
            rows = path.read_text(encoding="utf-8").replace(
                "0,10,140,10.0,state_counts,1,0.0",
                "0,10,140,10.0,state_counts,1,1.0",
            )
            path.write_text(rows, encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "still has defects"):
                aging.analyze_density(path, density=0.25, aged_step=140)


if __name__ == "__main__":
    unittest.main()
