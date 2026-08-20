from __future__ import annotations

import csv
import importlib.util
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "muscovite_release.py"
SPEC = importlib.util.spec_from_file_location("muscovite_release", MODULE_PATH)
assert SPEC and SPEC.loader
muscovite_release = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = muscovite_release
SPEC.loader.exec_module(muscovite_release)


class CylinderInversionTests(unittest.TestCase):
    def test_round_trip_across_release_range(self) -> None:
        for fraction in (0.01, 0.10, 0.30, 0.60, 0.90, 0.99):
            tau = muscovite_release.invert_cylinder_fraction(fraction)
            self.assertAlmostEqual(
                muscovite_release.cylinder_fraction(tau), fraction, places=10
            )

    def test_short_time_limit(self) -> None:
        fraction = 0.02
        expected = (
            2.0 / (3.141592653589793**0.5)
            - (4.0 / 3.141592653589793 - fraction) ** 0.5
        ) ** 2
        self.assertAlmostEqual(
            muscovite_release.invert_cylinder_fraction(fraction),
            expected,
            places=12,
        )

    def test_crank_series_reference_values(self) -> None:
        self.assertAlmostEqual(
            muscovite_release.cylinder_fraction(0.01),
            0.21547393817949312,
            places=12,
        )
        self.assertAlmostEqual(
            muscovite_release.cylinder_fraction(0.1),
            0.6058241939666917,
            places=12,
        )


class EventParsingTests(unittest.TestCase):
    def test_run_analysis_uses_release_reactions_only(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            run = Path(directory)
            with (run / "populations.csv").open(
                "w", newline="", encoding="utf-8"
            ) as handle:
                writer = csv.writer(handle)
                writer.writerow(
                    ["step", "time", "K_gallery_A.Ar40", "K_gallery_B.Ar40"]
                )
                writer.writerow([0, 0.0, 60, 40])
            header = {
                "petra_traj": 1,
                "reactions": ["hop", "release_A_surface", "release_B_surface"],
            }
            with (run / "events.jsonl").open("w", encoding="utf-8") as handle:
                handle.write(json.dumps(header) + "\n")
                for index in range(1, 81):
                    reaction = 1 if index % 2 else 2
                    handle.write(json.dumps([index, float(index), reaction, []]) + "\n")
                handle.write(json.dumps([81, 81.0, 0, []]) + "\n")
            initial = muscovite_release._initial_argon(run / "populations.csv")
            times, counts = muscovite_release._release_events(run / "events.jsonl")
            self.assertEqual(initial, 100)
            self.assertEqual(len(times), 80)
            self.assertEqual(counts["hop"], 1)
            self.assertEqual(counts["release_A_surface"], 40)
            self.assertEqual(counts["release_B_surface"], 40)
            points = muscovite_release.quantile_points(times, initial)
            self.assertEqual(points[-1].fraction, 0.80)


class QualitativeGateTests(unittest.TestCase):
    def test_gate_requires_rise_then_order_of_magnitude_fall(self) -> None:
        values = [1.0, 2.0, 3.0, 2.5, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05]
        points = tuple(
            muscovite_release.ReleasePoint(
                fraction=0.02 * (index + 1),
                released=index + 1,
                time_s=float(index + 1),
                sqrt_time_sqrt_s=(index + 1) ** 0.5,
                cylinder_tau=0.01 * (index + 1),
                apparent_da2_per_s=value,
            )
            for index, value in enumerate(values)
        )
        gate = muscovite_release.evaluate_gate(points, 100, 80)
        self.assertTrue(gate.early_rise)
        self.assertTrue(gate.later_fall)
        self.assertTrue(gate.pass_all)

    def test_flat_curve_fails(self) -> None:
        points = tuple(
            muscovite_release.ReleasePoint(
                fraction=0.02 * (index + 1),
                released=index + 1,
                time_s=float(index + 1),
                sqrt_time_sqrt_s=(index + 1) ** 0.5,
                cylinder_tau=0.01 * (index + 1),
                apparent_da2_per_s=1.0,
            )
            for index in range(10)
        )
        gate = muscovite_release.evaluate_gate(points, 100, 80)
        self.assertFalse(gate.pass_all)

    def test_small_population_cannot_pass_on_quantile_noise(self) -> None:
        values = [1.0, 2.0, 3.0, 2.5, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05]
        points = tuple(
            muscovite_release.ReleasePoint(
                fraction=0.02 * (index + 1),
                released=index + 1,
                time_s=float(index + 1),
                sqrt_time_sqrt_s=(index + 1) ** 0.5,
                cylinder_tau=0.01 * (index + 1),
                apparent_da2_per_s=value,
            )
            for index, value in enumerate(values)
        )
        gate = muscovite_release.evaluate_gate(points, 50, 40)
        self.assertFalse(gate.enough_release)
        self.assertFalse(gate.pass_all)

    def test_fewer_than_six_crossings_returns_failing_gate_not_raise(self) -> None:
        points = tuple(
            muscovite_release.ReleasePoint(
                fraction=0.02 * (index + 1),
                released=index + 1,
                time_s=float(index + 1),
                sqrt_time_sqrt_s=(index + 1) ** 0.5,
                cylinder_tau=0.01 * (index + 1),
                apparent_da2_per_s=1.0,
            )
            for index in range(3)
        )
        gate = muscovite_release.evaluate_gate(points, 100, 6)
        self.assertFalse(gate.pass_all)
        self.assertFalse(gate.early_rise)
        self.assertFalse(gate.later_fall)
        self.assertTrue(math.isnan(gate.rise_ratio))
        self.assertTrue(math.isnan(gate.fall_ratio))


class SvgPanelTests(unittest.TestCase):
    def _result(self, pass_all: bool):
        values = [1.0, 2.0, 3.0, 2.5, 2.0, 1.0, 0.5, 0.2, 0.1, 0.05]
        points = tuple(
            muscovite_release.ReleasePoint(
                fraction=0.02 * (index + 1),
                released=index + 1,
                time_s=float(index + 1),
                sqrt_time_sqrt_s=(index + 1) ** 0.5,
                cylinder_tau=0.01 * (index + 1),
                apparent_da2_per_s=value,
            )
            for index, value in enumerate(values)
        )
        gate = muscovite_release.evaluate_gate(points, 100, 80)
        if not pass_all:
            gate = muscovite_release.evaluate_gate(points, 50, 40)
        return muscovite_release.RunResult(
            label="run",
            temperature_c=500.0,
            run_dir=".",
            reaction_counts={},
            points=points,
            gate=gate,
        )

    def test_svg_verdict_text_matches_gate(self) -> None:
        passing = muscovite_release._svg_panels((self._result(True),), "release")
        self.assertIn("gate PASS", passing)
        self.assertNotIn("gate FAIL", passing)
        failing = muscovite_release._svg_panels((self._result(False),), "release")
        self.assertIn("gate FAIL", failing)
        self.assertNotIn("gate PASS", failing)

    def test_svg_verdict_color_matches_gate(self) -> None:
        passing = muscovite_release._svg_panels((self._result(True),), "release")
        failing = muscovite_release._svg_panels((self._result(False),), "release")
        self.assertIn('fill="#166534"', passing)
        self.assertIn('fill="#b42318"', failing)


if __name__ == "__main__":
    unittest.main()
