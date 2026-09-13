from __future__ import annotations

import importlib.util
import json
import math
import os
import sys
import tempfile
import tomllib
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "approximate_rate_closure.py"
SPEC = importlib.util.spec_from_file_location("approximate_rate_closure", MODULE_PATH)
assert SPEC and SPEC.loader
closure = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = closure
SPEC.loader.exec_module(closure)
DECK = Path(__file__).resolve().parents[2] / "examples" / "kaolinite-approx.toml"
LEGACY_DECK = DECK.with_name("kaolinite.toml")


def _write_sabotage(root: Path, replacement: tuple[str, str]) -> Path:
    path = root / "sabotaged.toml"
    source = DECK.read_text(encoding="utf-8")
    old, new = replacement
    assert source.count(old) == 1
    path.write_text(source.replace(old, new), encoding="utf-8")
    return path


def _points(increments: list[int], *, area_final: float = 100.0) -> list:
    cumulative = 0
    points = [closure.SteadyPoint(0, 0.0, 0, 100.0, 100)]
    for index, increment in enumerate(increments, start=1):
        cumulative += increment
        area = area_final if index == len(increments) else 100.0
        solid = 0 if area == 0 else 100
        points.append(
            closure.SteadyPoint(index * 1000, float(index), cumulative, area, solid)
        )
    return points


class DeckContractTests(unittest.TestCase):
    def test_canonical_deck_and_family_partition(self) -> None:
        contract = closure.validate_deck(DECK)
        self.assertEqual(len(contract.barriers), 22)
        self.assertEqual(set(contract.barriers), set(closure.EXPECTED_REACTIONS))
        self.assertEqual(
            {family for family, _ in contract.annotations.values()},
            set(closure.FAMILY_REACTIONS),
        )
        for family, names in closure.FAMILY_REACTIONS.items():
            self.assertEqual(
                {
                    name
                    for name, annotation in contract.annotations.items()
                    if annotation[0] == family
                },
                set(names),
            )

    def test_topology_init_and_reaction_semantics_match_legacy(self) -> None:
        old = tomllib.loads(LEGACY_DECK.read_text(encoding="utf-8"))
        new = tomllib.loads(DECK.read_text(encoding="utf-8"))
        for key in ("species", "kinds", "aliases", "cell", "lattice", "init"):
            self.assertEqual(new[key], old[key], key)
        self.assertEqual(len(new["reactions"]), 22)
        for old_reaction, new_reaction in zip(
            old["reactions"], new["reactions"], strict=True
        ):
            self.assertEqual(
                {key: value for key, value in new_reaction.items() if key != "rate"},
                {key: value for key, value in old_reaction.items() if key != "rate"},
            )

    def test_exact_temperature_sabotage_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = _write_sabotage(
                Path(directory), ("temperature = 298.0", "temperature = 298.15")
            )
            with self.assertRaisesRegex(ValueError, "exactly.*298.0"):
                closure.validate_deck(path)

    def test_schedule_and_wrong_prefactor_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            scheduled = _write_sabotage(
                root,
                (
                    "[simulation]\n",
                    "[[execution.schedule]]\ntemperature = 298.0\nduration = 1.0\n\n[simulation]\n",
                ),
            )
            with self.assertRaisesRegex(ValueError, "schedules are forbidden"):
                closure.validate_deck(scheduled)
            wrong = _write_sabotage(
                root,
                (
                    'center = { kind = "Oss", state = ["br"] }\nguards = [{ distance = 2, state = ["@hydrolyzed"], frozen = false, min = 1 }]\nrate = { arrhenius = { prefactor = 1.0e13, ea = 27.019 } }',
                    'center = { kind = "Oss", state = ["br"] }\nguards = [{ distance = 2, state = ["@hydrolyzed"], frozen = false, min = 1 }]\nrate = { arrhenius = { prefactor = 1.0e12, ea = 27.019 } }',
                ),
            )
            with self.assertRaisesRegex(ValueError, "prefactor must be exactly"):
                closure.validate_deck(wrong)

    def test_missing_observable_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = _write_sabotage(
                Path(directory),
                (
                    '[[observables.series]]\nkind = "exposure_age"\nstate = "@solid"\naxis = 1',
                    "# exposure age observable sabotaged",
                ),
            )
            with self.assertRaisesRegex(ValueError, "missing required observables"):
                closure.validate_deck(path)

    def test_seven_replicas_and_duplicate_seeds_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "at least 8"):
            closure.validate_seeds(range(7))
        with self.assertRaisesRegex(ValueError, "duplicate"):
            closure.validate_seeds((1, 2, 3, 4, 5, 6, 7, 7))

    def test_sensitivity_changes_only_the_named_family_rates(self) -> None:
        contract = closure.validate_deck(DECK)
        nominal = {
            reaction["name"]: reaction for reaction in contract.parsed["reactions"]
        }
        for family in closure.FAMILY_REACTIONS:
            for perturbation in closure.PERTURBATIONS:
                first = closure.perturb_deck(contract, family, perturbation)
                self.assertEqual(
                    first, closure.perturb_deck(contract, family, perturbation)
                )
                changed = {
                    reaction["name"]
                    for reaction in tomllib.loads(first)["reactions"]
                    if reaction != nominal[reaction["name"]]
                }
                self.assertEqual(changed, set(closure.FAMILY_REACTIONS[family]))

    def test_ea_minus_three_fails_before_a_negative_barrier(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = _write_sabotage(
                Path(directory),
                (
                    "prefactor = 1.0e13, ea = 6.400",
                    "prefactor = 1.0e13, ea = 2.400",
                ),
            )
            contract = closure.validate_deck(path)
            with self.assertRaisesRegex(ValueError, "would make.*negative"):
                closure.perturb_deck(contract, "adsorption", "ea-minus-3")


class EventAndUnitTests(unittest.TestCase):
    def test_event_gross_and_net_accounting_by_reaction_name(self) -> None:
        contract = closure.validate_deck(DECK)
        states = [
            f"{kind['name']}.{state['name']}"
            for kind in contract.parsed["kinds"]
            for state in kind["states"]
        ]
        state_id = {name: index for index, name in enumerate(states)}
        reactions = [reaction["name"] for reaction in contract.parsed["reactions"]]
        reaction_id = {name: index for index, name in enumerate(reactions)}
        rows = [
            [
                1,
                1.0,
                reaction_id["desorb-si"],
                [[4, state_id["Si.oh4"], state_id["Si.empty"]]],
            ],
            [
                2,
                2.0,
                reaction_id["adsorb-si"],
                [[4, state_id["Si.empty"], state_id["Si.oh4"]]],
            ],
            [
                3,
                3.0,
                reaction_id["desorb-al"],
                [[1, state_id["Al.l6"], state_id["Al.empty"]]],
            ],
            [
                4,
                4.0,
                reaction_id["adsorb-al"],
                [[1, state_id["Al.empty"], state_id["Al.l6"]]],
            ],
        ]
        header = {
            "petra_traj": 1,
            "deck": contract.parsed["deck"]["name"],
            "seed": 99,
            "n_sites": 10,
            "states": states,
            "reactions": reactions,
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "events.jsonl"
            path.write_text(
                "\n".join([json.dumps(header), *(json.dumps(row) for row in rows)])
                + "\n",
                encoding="utf-8",
            )
            events = closure.parse_events(path, contract, expected_seed=99)
        self.assertEqual(
            closure.dissolution_event_counts(events, 0.0, 4.0),
            {
                "gross_si": 1,
                "adsorb_si": 1,
                "net_si": 0,
                "gross_al": 1,
                "adsorb_al": 1,
                "net_al": 0,
            },
        )

    def test_exact_unit_conversion_oracle(self) -> None:
        self.assertEqual(closure.AVOGADRO_EXACT, 6.02214076e23)
        self.assertEqual(closure.ANGSTROM2_TO_M2, 1.0e-20)
        # 602,214,076 events = 1e-15 mol exactly by the SI-defined N_A;
        # 1e5 A^2 s = 1e-15 m^2 s, hence exactly 1 mol m^-2 s^-1.
        self.assertAlmostEqual(
            closure.event_count_to_flux(602_214_076, 1.0e5), 1.0, places=15
        )

    def test_evolving_area_uses_trapezoidal_physical_time_integral(self) -> None:
        self.assertEqual(
            closure.integrate_area([0.0, 1.0, 3.0], [10.0, 20.0, 40.0]),
            75.0,
        )


class SteadyStateGateTests(unittest.TestCase):
    def test_constant_per_replica_tail_passes(self) -> None:
        gate = closure.assess_steady_state(_points([2, 2, 2, 2, 2, 2, 2]), 7000)
        self.assertTrue(gate.passed, gate.reasons)
        self.assertEqual(gate.gross_events, 12)

    def test_early_no_event_insufficient_cadence_and_absorbing_fail(self) -> None:
        self.assertFalse(
            closure.assess_steady_state(_points([2] * 7), 8000).passed,
            "early stop must fail",
        )
        self.assertFalse(closure.assess_steady_state(_points([0] * 7), 7000).passed)
        self.assertFalse(closure.assess_steady_state(_points([2] * 5), 5000).passed)
        self.assertFalse(
            closure.assess_steady_state(_points([2] * 7, area_final=0.0), 7000).passed
        )

    def test_monotonic_trend_fails(self) -> None:
        gate = closure.assess_steady_state(_points([1, 2, 4, 8, 16, 32, 64]), 7000)
        self.assertFalse(gate.passed)
        self.assertTrue(
            any("trend" in reason or "shift" in reason for reason in gate.reasons)
        )

    def test_opposite_replica_trends_cannot_cancel_in_mean(self) -> None:
        increasing = closure.assess_steady_state(
            _points([1, 2, 4, 8, 16, 32, 64]), 7000
        )
        decreasing_increments = [64, 63, 61, 57, 49, 33, 1]
        decreasing = closure.assess_steady_state(_points(decreasing_increments), 7000)
        self.assertFalse(increasing.passed)
        self.assertFalse(decreasing.passed)
        self.assertEqual(
            [
                a + b
                for a, b in zip(
                    [1, 2, 4, 8, 16, 32, 64], decreasing_increments, strict=True
                )
            ],
            [65] * 7,
            "ensemble averaging would hide the opposed trends",
        )


class DeterminismTests(unittest.TestCase):
    def test_bootstrap_and_verification_receipts_are_byte_deterministic(self) -> None:
        values = [1.0, 2.0, 3.0, 5.0, 8.0, 13.0, 21.0, 34.0]
        first = closure.bootstrap_summary(values, "fixed-derived-output")
        self.assertEqual(
            first, closure.bootstrap_summary(values, "fixed-derived-output")
        )
        self.assertTrue(all(math.isfinite(value) for value in first))
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "receipt.json"
            closure.write_json_atomic(
                path, {"summary": first, "seeds": list(closure.DEFAULT_SEEDS)}
            )
            payload = path.read_bytes()
            closure.write_json_atomic(
                path, {"summary": first, "seeds": list(closure.DEFAULT_SEEDS)}
            )
            self.assertEqual(path.read_bytes(), payload)

    def test_direct_runner_uses_required_flags_caps_and_hash_receipt(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fake = root / "fake-petra"
            fake.write_text(
                """#!/usr/bin/env python3
import json
import os
import pathlib
import sys
args = sys.argv[1:]
out = pathlib.Path(args[args.index('--out') + 1])
out.mkdir()
payload = {'args': args, 'threads': {name: os.environ[name] for name in (
    'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'RAYON_NUM_THREADS')}}
(out / 'events.jsonl').write_text(json.dumps(payload) + '\\n')
for name in ('populations.csv', 'observables.csv', 'snapshot.pgif.json'):
    (out / name).write_text(name + '\\n')
""",
                encoding="utf-8",
            )
            fake.chmod(0o755)
            deck = root / "deck.toml"
            deck.write_text("fixture", encoding="utf-8")
            output = root / "output"
            log = root / "run.log"
            env = os.environ.copy()
            for variable in (
                "OMP_NUM_THREADS",
                "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "RAYON_NUM_THREADS",
            ):
                env[variable] = "16"
            receipt = closure._run_one(fake, deck, output, log, 90401, 30, root, env)
            event_receipt = json.loads((output / "events.jsonl").read_text())
            self.assertEqual(
                event_receipt["threads"],
                {name: "16" for name in event_receipt["threads"]},
            )
            args = event_receipt["args"]
            self.assertIn("--viz", args)
            self.assertIn("--paranoid", args)
            self.assertEqual(args[args.index("--ensemble") + 1], "1")
            self.assertEqual(args[args.index("--seed") + 1], "90401")
            self.assertEqual(
                receipt["sha256"]["events.jsonl"],
                closure.sha256_file(output / "events.jsonl"),
            )
            with self.assertRaisesRegex(ValueError, "refusing to overwrite"):
                closure._run_one(fake, deck, output, log, 90401, 30, root, env)


if __name__ == "__main__":
    unittest.main()
