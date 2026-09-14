from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

import tomllib

SCRIPTS = Path(__file__).resolve().parents[1]
PETRA = SCRIPTS.parent
DATA = PETRA / "data" / "muscovite-1998"
DECKS = PETRA / "decks" / "muscovite-1998"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


load_module("build_muscovite_full_deck", SCRIPTS / "build_muscovite_full_deck.py")
load_module("muscovite_full_analysis", SCRIPTS / "muscovite_full_analysis.py")
load_module("muscovite_grain_size_sweep", SCRIPTS / "muscovite_grain_size_sweep.py")
comparison = load_module(
    "muscovite_1998_comparison", SCRIPTS / "muscovite_1998_comparison.py"
)


class ComparisonDeckTests(unittest.TestCase):
    def test_replay_evidence_records_matching_sizes_and_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            roots = [Path(temporary) / suffix for suffix in ("a", "b")]
            for root in roots:
                root.mkdir()
                for filename in comparison.REPLAY_ARTIFACTS:
                    (root / filename).write_bytes(b"same-seed-output\n")

            evidence = comparison._replay_evidence(roots)

        self.assertEqual(set(evidence), set(comparison.REPLAY_ARTIFACTS))
        for pair in evidence.values():
            self.assertEqual(pair["replay_a"], pair["replay_b"])
            self.assertEqual(pair["replay_a"]["bytes"], 17)
            self.assertEqual(len(pair["replay_a"]["sha256"]), 64)

    def test_replay_evidence_rejects_divergence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            roots = [Path(temporary) / suffix for suffix in ("a", "b")]
            for root in roots:
                root.mkdir()
                for filename in comparison.REPLAY_ARTIFACTS:
                    (root / filename).write_bytes(b"same-seed-output\n")
            (roots[1] / "observables.csv").write_bytes(b"different\n")

            with self.assertRaisesRegex(RuntimeError, "observables.csv"):
                comparison._replay_evidence(roots)

    def test_replica_prefix_evidence_ignores_later_ensemble_members(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = [Path(temporary) / name for name in ("primary.csv", "replay.csv")]
            paths[0].write_text(
                "replica,seed,value\n0,10,a\n1,11,b\n2,12,c\n", encoding="utf-8"
            )
            paths[1].write_text(
                "replica,seed,value\n0,10,a\n1,11,b\n", encoding="utf-8"
            )

            self.assertEqual(
                comparison._replica_prefix_evidence(paths[0]),
                comparison._replica_prefix_evidence(paths[1]),
            )

    def test_artifact_verification_rejects_post_run_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "observables.csv"
            path.write_bytes(b"campaign-output\n")
            evidence = comparison._artifact_evidence(path)
            path.write_bytes(b"tampered-output\n")

            with self.assertRaisesRegex(RuntimeError, "artifact drift"):
                comparison._verify_artifact(path, evidence)

    def test_all_six_tracked_decks_are_deterministically_generated(self) -> None:
        expected = []
        for barrier_label in comparison.DELAMINATION_SENSITIVITY_KCAL_MOL:
            for dims in comparison.COMPARISON_SIZES:
                slug = "x".join(str(value) for value in dims)
                path = DECKS / f"muscovite-1998-{barrier_label}-{slug}.toml"
                expected.append(path)
                self.assertEqual(
                    comparison.render_comparison_deck(dims, barrier_label),
                    path.read_text(encoding="utf-8"),
                )
        self.assertEqual(len(expected), 6)

    def test_decks_wire_e3a_barriers_and_comparison_schedule(self) -> None:
        for (
            barrier_label,
            barrier,
        ) in comparison.DELAMINATION_SENSITIVITY_KCAL_MOL.items():
            path = DECKS / f"muscovite-1998-{barrier_label}-8x8x6.toml"
            deck = tomllib.loads(path.read_text(encoding="utf-8"))
            rules = {rule["name"]: rule for rule in deck["dynamics"]["rules"]}
            self.assertEqual(
                rules["delaminate_interface"]["rate"]["eyring"]["dh"], barrier
            )
            self.assertEqual(
                rules["hop_pristine_Ar40_A_to_B_intra"]["rate"]["eyring"]["dh"],
                comparison.E3A_LOCAL_DEHYDROXYLATE_HOP_KCAL_MOL,
            )
            self.assertEqual(
                rules["hop_extended_Ar40_A_to_B_intra"]["rate"]["eyring"]["dh"],
                comparison.E3A_EXTENDED_ZONE_HOP_KCAL_MOL,
            )
            self.assertIn("incomplete-input-contract", path.read_text(encoding="utf-8"))
            self.assertIn("incomplete-convergence", path.read_text(encoding="utf-8"))
            self.assertEqual(deck["structure"]["lattice"]["dims"], [8, 8, 6])
            schedule = deck["execution"]["schedule"]
            self.assertEqual(len(schedule), 15)
            self.assertEqual(schedule[0], {"temperature": 773.15, "duration": 600.0})
            self.assertEqual(schedule[-1], {"temperature": 1473.15, "duration": 600.0})

    def test_generate_refuses_no_content_drift(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = comparison.generate_decks(Path(temporary))
            self.assertEqual(len(paths), 6)
            self.assertEqual(len({path.read_bytes() for path in paths}), 6)


class DigitizedDataTests(unittest.TestCase):
    def test_release_digitization_has_three_grain_sizes_and_uncertainty(self) -> None:
        rows = comparison._read_release_data(DATA / "figure4-release-rates.csv")
        sizes = {
            int(row["grain_size_um"])
            for row in rows
            if row["experiment"].startswith("500C_isothermal_")
        }
        self.assertEqual(sizes, {58, 165, 3000})
        self.assertTrue(
            all(float(row["estimated_uncertainty_x1e6_per_c"]) > 0 for row in rows)
        )
        peaks = {
            name: float(
                max(
                    (row for row in rows if row["experiment"] == name),
                    key=lambda row: float(row["release_rate_x1e6_per_c"]),
                )["temperature_c"]
            )
            for name in (
                "500C_isothermal_small",
                "500C_isothermal_medium",
                "500C_isothermal_large",
            )
        }
        self.assertLess(peaks["500C_isothermal_small"], peaks["500C_isothermal_large"])

    def test_age_digitization_covers_staircase_and_old_initial_steps(self) -> None:
        rows = comparison._read_age_data(DATA / "figure3-and-5-age-spectra.csv")
        experiments = {row["experiment"] for row in rows}
        self.assertIn("hydrothermal_700C_2kbar_45d", experiments)
        self.assertIn("500C_isothermal_old_steps", experiments)
        old = [row for row in rows if row["experiment"] == "500C_isothermal_old_steps"]
        self.assertGreater(max(float(row["age_ma"]) for row in old), 1000.0)
        self.assertTrue(all(float(row["estimated_uncertainty_ma"]) > 0 for row in rows))


class EvaluationTests(unittest.TestCase):
    @staticmethod
    def synthetic_rows() -> list[dict[str, object]]:
        rows = []
        for barrier in comparison.DELAMINATION_SENSITIVITY_KCAL_MOL:
            for size_index, dims in enumerate(comparison.COMPARISON_SIZES):
                cumulative = 0.0
                for segment, temperature in enumerate((500.0, 700.0, 900.0, 1100.0), 1):
                    release = (4 - size_index) if segment == size_index + 1 else 0.2
                    cumulative = min(1.0, cumulative + 0.25)
                    rows.append(
                        {
                            "barrier_label": barrier,
                            "dims": "x".join(str(value) for value in dims),
                            "temperature_c": temperature,
                            "released_ar40_mean": release,
                            "released_ar40_ci95_low": max(0.0, release - 0.1),
                            "released_ar40_ci95_high": release + 0.1,
                            "released_ar39_mean": 1.0 if segment == 1 else 5.0,
                            "cumulative_ar40_fraction_mean": cumulative,
                            "cumulative_ar39_fraction_mean": cumulative,
                            "apparent_age_ma_mean": 1050.0 if segment == 1 else 900.0,
                            "apparent_age_ma_ci95_low": 1030.0
                            if segment == 1
                            else 880.0,
                            "apparent_age_ma_ci95_high": 1070.0
                            if segment == 1
                            else 920.0,
                        }
                    )
        return rows

    def test_verdict_reports_observed_and_each_e3a_bracket(self) -> None:
        release = comparison._read_release_data(DATA / "figure4-release-rates.csv")
        verdict = comparison.evaluate(self.synthetic_rows(), release)
        self.assertTrue(verdict["observed"]["grain_size_crossover"])
        self.assertEqual(set(verdict["synthetic"]), {"low", "high"})
        self.assertEqual(len(verdict["section_5_claims"]), 7)
        self.assertEqual(
            verdict["section_5_claims"]["6_grain_size_delamination_fraction"][
                "verdict"
            ],
            "reproduced",
        )
        self.assertEqual(
            verdict["section_5_claims"]["5_recoil_old_initial_steps"]["verdict"],
            "reproduced",
        )
        self.assertIn(
            "6 of 6 volume/sensitivity families",
            verdict["section_5_claims"]["5_recoil_old_initial_steps"]["mechanism"],
        )
        self.assertIn(
            "surface-connected lateral",
            verdict["section_5_claims"]["6_grain_size_delamination_fraction"][
                "mechanism"
            ],
        )
        self.assertIn(
            "surface-connected lateral",
            verdict["section_5_claims"]["1_two_stage_non_fickian_loss"]["mechanism"],
        )
        self.assertTrue(verdict["synthetic"]["low"]["grain_size_crossover_reproduced"])
        self.assertTrue(
            all(
                verdict["synthetic"]["low"][
                    "recoil_driven_distortion_reproduced"
                ].values()
            )
        )
        self.assertTrue(
            verdict["synthetic"]["high"]["high_temperature_merge_reproduced"]
        )

    def test_crossover_requires_three_point_volume_ordering(self) -> None:
        release = comparison._read_release_data(DATA / "figure4-release-rates.csv")
        rows = self.synthetic_rows()
        for row in rows:
            if row["dims"] == "8x8x6":
                row["released_ar40_mean"] = (
                    10.0 if row["temperature_c"] == 500.0 else 0.2
                )
        verdict = comparison.evaluate(rows, release)
        self.assertFalse(verdict["synthetic"]["low"]["grain_size_crossover_reproduced"])
        self.assertTrue(verdict["observed"]["grain_size_crossover"])
        self.assertTrue(
            comparison._monotonic_volume_crossover(
                {"4x4x6": 700.0, "8x8x6": 800.0, "12x12x6": 1025.0}
            )
        )
        self.assertFalse(
            comparison._monotonic_volume_crossover(
                {"4x4x6": 700.0, "8x8x6": 600.0, "12x12x6": 1025.0}
            )
        )

    def test_svg_contains_observed_and_synthetic_layers(self) -> None:
        release = comparison._read_release_data(DATA / "figure4-release-rates.csv")
        ages = comparison._read_age_data(DATA / "figure3-and-5-age-spectra.csv")
        svg = comparison.comparison_svg(self.synthetic_rows(), release, ages)
        ET.fromstring(svg)
        self.assertTrue(svg.startswith("<svg"))
        self.assertIn("solid: Fig. 4", svg)
        self.assertIn("black: Fig. 3d", svg)
        self.assertGreaterEqual(svg.count("<polygon"), 5)
        self.assertIn("12x12x6", str(self.synthetic_rows()))
        self.assertGreaterEqual(svg.count("<polyline"), 8)


if __name__ == "__main__":
    unittest.main()
