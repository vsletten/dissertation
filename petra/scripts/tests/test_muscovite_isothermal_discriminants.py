from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

import tomllib

SCRIPTS = Path(__file__).resolve().parents[1]
PETRA = SCRIPTS.parent
DATA = PETRA / "data" / "muscovite-1998"


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


load("build_muscovite_full_deck", SCRIPTS / "build_muscovite_full_deck.py")
release = load("muscovite_release", SCRIPTS / "muscovite_release.py")
e4b = load(
    "muscovite_isothermal_discriminants",
    SCRIPTS / "muscovite_isothermal_discriminants.py",
)


class E4bDeckTests(unittest.TestCase):
    def test_twelve_decks_generate_and_parse_deterministically(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            paths = e4b.generate_decks(Path(temporary))
            self.assertEqual(len(paths), 12)
            for path in paths:
                first = path.read_text(encoding="utf-8")
                deck = tomllib.loads(first)
                self.assertEqual(len(deck["execution"]["schedule"]), 11)
                self.assertEqual(
                    sum(step["duration"] for step in deck["execution"]["schedule"]),
                    216000.0,
                )
                self.assertEqual(deck["execution"]["schedule"][-1]["duration"], 66000.0)
                self.assertIn("E4b isothermal discriminator", deck["deck"]["comment"])
                self.assertEqual(
                    first,
                    e4b.render_e4b_deck(
                        tuple(deck["structure"]["lattice"]["dims"]),
                        round(deck["dynamics"]["thermo"]["temperature"] - 273.15),
                        "hydrothermal-2kbar-screen"
                        if deck["dynamics"]["thermo"]["activity"]["H2O_vacancy"] < 1
                        else "vacuum",
                    ),
                )

    def test_reservoir_xe_and_h2o_hypotheses_are_explicit_and_labelled(self) -> None:
        text = e4b.render_e4b_deck((4, 4, 6), 700, "hydrothermal-2kbar-screen")
        deck = tomllib.loads(text)
        rules = {rule["name"]: rule for rule in deck["dynamics"]["rules"]}
        self.assertEqual(deck["dynamics"]["thermo"]["activity"]["H2O_vacancy"], 5.0e-4)
        self.assertEqual(rules["dehydroxylate_pair"]["consumes"], ["H2O_vacancy"])
        self.assertEqual(
            rules["hop_extended_Ar36_A_to_B_intra"]["rate"]["eyring"]["dh"], 40.0
        )
        self.assertEqual(
            rules["hop_extended_Ar40_A_to_B_intra"]["rate"]["eyring"]["dh"],
            e4b.E3A_AR_EXTENDED_KCAL_MOL,
        )
        self.assertEqual(
            rules["hop_pristine_Ar40_A_to_B_intra"]["rate"]["eyring"]["dh"],
            e4b.E3A_AR_LOCAL_KCAL_MOL,
        )
        self.assertEqual(
            rules["hop_extended_Xe_A_to_B_intra"]["rate"]["eyring"]["dh"],
            e4b.E3A_XE_SCREEN_KCAL_MOL,
        )
        self.assertIn("incomplete-convergence", text)
        self.assertIn("never a production calibration", text)

    def test_figure_digitization_covers_6_through_11_with_uncertainty(self) -> None:
        import csv

        with (DATA / "figures6-11-time-series.csv").open(
            newline="", encoding="utf-8"
        ) as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual({int(row["figure"]) for row in rows}, set(range(6, 12)))
        self.assertTrue(all(float(row["estimated_x_uncertainty"]) > 0 for row in rows))
        self.assertTrue(all(float(row["estimated_y_uncertainty"]) > 0 for row in rows))
        self.assertEqual({int(row["grain_size_um"]) for row in rows}, {58, 165, 3000})


class AggregationTests(unittest.TestCase):
    def test_bootstrap_mean_interval_is_deterministic(self) -> None:
        values = [0.0, 0.25, 0.5, 0.75, 1.0]
        first = e4b.mean_ci95(values, seed=42)
        second = e4b.mean_ci95(values, seed=42)
        self.assertEqual(first, second)
        self.assertEqual(first[0], 0.5)
        self.assertLess(first[1], first[0])
        self.assertGreater(first[2], first[0])

    def test_single_replica_interval_is_not_fabricated(self) -> None:
        self.assertEqual(e4b.mean_ci95([3.5], seed=42), (3.5, 3.5, 3.5))
        self.assertEqual(e4b.mean_ci95([], seed=42), (None, None, None))


class CylinderInversionTests(unittest.TestCase):
    def test_hand_checkable_short_time_forward_inverse_vector(self) -> None:
        # Infinite-cylinder short-time branch: F=4*sqrt(tau/pi)-tau.
        # For tau=1e-4 the hand-computed fraction is 0.022467583341910253.
        tau = 1.0e-4
        fraction = 0.022467583341910253
        self.assertAlmostEqual(release.cylinder_fraction(tau), fraction, places=15)
        self.assertAlmostEqual(
            release.invert_cylinder_fraction(fraction), tau, places=12
        )
        # The apparent D/a² branch is Δtau/Δt; 1e-4 over 100 s = 1e-6 s⁻¹.
        self.assertAlmostEqual(tau / 100.0, 1.0e-6, places=15)


if __name__ == "__main__":
    unittest.main()
