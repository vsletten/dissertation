from __future__ import annotations

import csv
import importlib.util
import json
import math
import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

import tomllib

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


def _points(
    si_increments: list[int],
    *,
    al_increments: list[int] | None = None,
    area_final: float = 100.0,
    solid_si_final: int = 50,
    solid_al_final: int = 50,
    solid_si_samples: list[int] | None = None,
    solid_al_samples: list[int] | None = None,
    si_propensities: list[float] | None = None,
    al_propensities: list[float] | None = None,
) -> list:
    al_increments = si_increments if al_increments is None else al_increments
    sample_count = len(si_increments) + 1
    solid_si_samples = (
        [50] * (sample_count - 1) + [solid_si_final]
        if solid_si_samples is None
        else solid_si_samples
    )
    solid_al_samples = (
        [50] * (sample_count - 1) + [solid_al_final]
        if solid_al_samples is None
        else solid_al_samples
    )
    assert len(solid_si_samples) == sample_count
    assert len(solid_al_samples) == sample_count
    si_propensities = (
        [1.0] * sample_count if si_propensities is None else si_propensities
    )
    al_propensities = (
        [1.0] * sample_count if al_propensities is None else al_propensities
    )
    assert len(si_propensities) == sample_count
    assert len(al_propensities) == sample_count
    cumulative_si = 0
    cumulative_al = 0
    points = [
        closure.SteadyPoint(
            0,
            0.0,
            0,
            0,
            100.0,
            solid_si_samples[0],
            solid_al_samples[0],
            si_propensities[0],
            al_propensities[0],
        )
    ]
    for index, (si_increment, al_increment) in enumerate(
        zip(si_increments, al_increments, strict=True), start=1
    ):
        cumulative_si += si_increment
        cumulative_al += al_increment
        area = area_final if index == len(si_increments) else 100.0
        solid_si = solid_si_samples[index]
        solid_al = solid_al_samples[index]
        if area == 0:
            solid_si = solid_al = 0
        points.append(
            closure.SteadyPoint(
                index * 1000,
                float(index),
                cumulative_si,
                cumulative_al,
                area,
                solid_si,
                solid_al,
                si_propensities[index],
                al_propensities[index],
            )
        )
    return points


def _build_synthetic_campaign(root: Path) -> tuple[Path, Path]:
    raw = root / "raw"
    decks = raw / "decks"
    logs = raw / "logs"
    runs = raw / "runs"
    decks.mkdir(parents=True)
    logs.mkdir()
    runs.mkdir()
    binary = root / "fake-petra"
    binary.write_text("synthetic binary fixture\n", encoding="utf-8")
    binary.chmod(0o755)
    contract = closure.validate_deck(DECK)
    states = [
        f"{kind['name']}.{state['name']}"
        for kind in contract.parsed["kinds"]
        for state in kind["states"]
    ]
    state_ids = {name: index for index, name in enumerate(states)}
    reaction_names = [entry.name for entry in closure.REACTION_REGISTRY]
    reaction_ids = {name: index for index, name in enumerate(reaction_names)}
    scenario_records = []
    receipts = []

    for scenario in closure.scenarios():
        family_index = (
            list(closure.FAMILY_REACTIONS).index(scenario.family) + 1
            if scenario.family is not None
            else 0
        )
        propensity_multiplier = {
            None: 1.0,
            "ea-minus-3": 1.0 + 0.10 * family_index,
            "ea-plus-3": 1.0 / (1.0 + 0.05 * family_index),
            "prefactor-x0.1": 1.0 / (1.0 + 0.20 * family_index),
            "prefactor-x10": 1.0 + 0.20 * family_index,
        }[scenario.perturbation]
        text = (
            contract.text
            if scenario.family is None
            else closure.perturb_deck(contract, scenario.family, scenario.perturbation)
        )
        deck_path = decks / f"{scenario.name}.toml"
        deck_path.write_text(text, encoding="utf-8")
        scenario_records.append(
            {
                **closure.asdict(scenario),
                "deck": str(deck_path.resolve()),
                "deck_sha256": closure.sha256_file(deck_path),
            }
        )
        scenario_logs = logs / scenario.name
        scenario_runs = runs / scenario.name
        scenario_logs.mkdir()
        scenario_runs.mkdir()
        for replica, seed in enumerate(closure.DEFAULT_SEEDS):
            run_dir = scenario_runs / f"replica-{replica:02d}-seed-{seed}"
            run_dir.mkdir()
            log_path = scenario_logs / f"replica-{replica:02d}-seed-{seed}.log"
            log_path.write_text("synthetic completed run\n", encoding="utf-8")
            header = {
                "petra_traj": 1,
                "deck": contract.parsed["deck"]["name"],
                "seed": seed,
                "n_sites": 21,
                "states": states,
                "state_types": [
                    state["occupant"]
                    for kind in contract.parsed["kinds"]
                    for state in kind["states"]
                ],
                "reactions": reaction_names,
            }
            event_rows = []
            for index in range(1, 21):
                forward = index % 2 == 1
                event_rows.append(
                    [
                        index * 10_000,
                        float(index),
                        reaction_ids[
                            "R14-alohal-hydrolysis"
                            if forward
                            else "R15-alohal-condensation"
                        ],
                        [
                            [
                                0,
                                state_ids["Oaa.br" if forward else "Oaa.hy"],
                                state_ids["Oaa.hy" if forward else "Oaa.br"],
                            ]
                        ],
                    ]
                )
            (run_dir / "events.jsonl").write_text(
                "\n".join(
                    [json.dumps(header), *(json.dumps(row) for row in event_rows)]
                )
                + "\n",
                encoding="utf-8",
            )
            population_rows = []
            observable_rows = []
            for index in range(21):
                counts = {state: 0 for state in states}
                counts["Si.oh4"] = 10
                counts["Al.l6"] = 10
                counts["Oaa.br"] = int(index % 2 == 0)
                counts["Oaa.hy"] = int(index % 2 == 1)
                step = index * 10_000
                population_rows.append({"step": step, "time": float(index), **counts})
                observable_time = 0.0 if index == 0 else float(index) + 1.0e-7
                values = {
                    "state_counts": [counts[state] for state in states],
                    "event_rates": [0.0] * len(reaction_names),
                    "rate_spectra": [1.0],
                    "surface_area": [100.0, 1.0, 1.0],
                    "exposure_age": [observable_time],
                }
                values["event_rates"][reaction_ids["desorb-si"]] = (
                    2.0 + 0.1 * replica
                ) * propensity_multiplier
                values["event_rates"][reaction_ids["desorb-al"]] = (
                    1.0 + 0.05 * replica
                ) * propensity_multiplier
                for kind, kind_values in values.items():
                    for value_index, value in enumerate(kind_values):
                        observable_rows.append(
                            {
                                "replica": 0,
                                "seed": seed,
                                "step": step,
                                "time": observable_time,
                                "kind": kind,
                                "index": value_index,
                                "value": value,
                            }
                        )
            closure._write_csv_atomic(
                run_dir / "populations.csv",
                ["step", "time", *states],
                population_rows,
            )
            closure._write_csv_atomic(
                run_dir / "observables.csv",
                ["replica", "seed", "step", "time", "kind", "index", "value"],
                observable_rows,
            )
            (run_dir / "snapshot.pgif.json").write_text("{}\n", encoding="utf-8")
            artifact_hashes = {
                name: closure.sha256_file(run_dir / name)
                for name in (
                    "events.jsonl",
                    "populations.csv",
                    "observables.csv",
                    "snapshot.pgif.json",
                )
            }
            artifact_hashes["log"] = closure.sha256_file(log_path)
            receipts.append(
                {
                    "schema": "a9-run-receipt-v2",
                    "seed": seed,
                    "elapsed_seconds": 0.01,
                    "command": [
                        "nice",
                        "-n",
                        "10",
                        str(binary.resolve()),
                        str(deck_path.resolve()),
                        "--seed",
                        str(seed),
                        "--ensemble",
                        "1",
                        "--out",
                        str(run_dir.resolve()),
                        "--viz",
                        "--paranoid",
                    ],
                    "output": str(run_dir.resolve()),
                    "log": str(log_path.resolve()),
                    "sha256": artifact_hashes,
                    "scenario": scenario.name,
                    "replica": replica,
                }
            )
    checkpoint = raw / "checkpoint.json"
    closure.write_json_atomic(
        checkpoint,
        {
            "schema": closure.CHECKPOINT_SCHEMA,
            "status": "complete",
            "receipts": receipts,
        },
    )
    closure.write_json_atomic(
        raw / "manifest.json",
        {
            "schema": closure.RAW_SCHEMA,
            "status": "complete",
            "survey_tier": True,
            "temperature_k": 298.0,
            "units": "kcal/mol",
            "seeds": list(closure.DEFAULT_SEEDS),
            "replicas": closure.REPLICA_COUNT,
            "workers": 1,
            "timeout_seconds": 30,
            "source_deck": str(DECK.resolve()),
            "source_deck_sha256": closure.sha256_file(DECK),
            "petra_binary": str(binary.resolve()),
            "petra_binary_sha256": closure.sha256_file(binary),
            "scenarios": scenario_records,
            "completed_runs": len(receipts),
            "checkpoint_sha256": closure.sha256_file(checkpoint),
        },
    )
    return raw, root / "derived"


def _contaminate_one_run_with_adsorption(raw: Path, scenario: str) -> None:
    run = raw / "runs" / scenario / "replica-00-seed-90401" / "events.jsonl"
    rows = [json.loads(line) for line in run.read_text(encoding="utf-8").splitlines()]
    header = rows[0]
    state_ids = {name: index for index, name in enumerate(header["states"])}
    rows[-2] = [
        190_000,
        19.0,
        header["reactions"].index("adsorb-si"),
        [[4, state_ids["Si.empty"], state_ids["Si.oh4"]]],
    ]
    rows[-1] = [
        200_000,
        20.0,
        header["reactions"].index("R14-alohal-hydrolysis"),
        [[0, state_ids["Oaa.br"], state_ids["Oaa.hy"]]],
    ]
    run.write_text("\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8")
    checkpoint_path = raw / "checkpoint.json"
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    receipt = next(
        receipt
        for receipt in checkpoint["receipts"]
        if receipt["scenario"] == scenario and receipt["replica"] == 0
    )
    receipt["sha256"]["events.jsonl"] = closure.sha256_file(run)
    closure.write_json_atomic(checkpoint_path, checkpoint)
    manifest_path = raw / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["checkpoint_sha256"] = closure.sha256_file(checkpoint_path)
    closure.write_json_atomic(manifest_path, manifest)


def _make_one_run_propensity_nonstationary(raw: Path, scenario: str) -> None:
    path = raw / "runs" / scenario / "replica-00-seed-90401" / "observables.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        rows = list(reader)
    assert fields is not None
    desorb_indices = {
        str(index)
        for index, entry in enumerate(closure.REACTION_REGISTRY)
        if entry.name in {"desorb-si", "desorb-al"}
    }
    for row in rows:
        if row["kind"] == "event_rates" and row["index"] in desorb_indices:
            block = min(int(row["step"]) // 40_000, 4)
            row["value"] = str(float(row["value"]) * (2**block))
    closure._write_csv_atomic(path, fields, rows)
    checkpoint_path = raw / "checkpoint.json"
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    receipt = next(
        receipt
        for receipt in checkpoint["receipts"]
        if receipt["scenario"] == scenario and receipt["replica"] == 0
    )
    receipt["sha256"]["observables.csv"] = closure.sha256_file(path)
    closure.write_json_atomic(checkpoint_path, checkpoint)
    manifest_path = raw / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["checkpoint_sha256"] = closure.sha256_file(checkpoint_path)
    closure.write_json_atomic(manifest_path, manifest)


def _zero_all_desorption_propensities(raw: Path) -> None:
    checkpoint_path = raw / "checkpoint.json"
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    desorb_indices = {
        str(index)
        for index, entry in enumerate(closure.REACTION_REGISTRY)
        if entry.name in {"desorb-si", "desorb-al"}
    }
    for receipt in checkpoint["receipts"]:
        path = Path(receipt["output"]) / "observables.csv"
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            fields = reader.fieldnames
            rows = list(reader)
        assert fields is not None
        for row in rows:
            if row["kind"] == "event_rates" and row["index"] in desorb_indices:
                row["value"] = "0.0"
        closure._write_csv_atomic(path, fields, rows)
        receipt["sha256"]["observables.csv"] = closure.sha256_file(path)
    closure.write_json_atomic(checkpoint_path, checkpoint)
    manifest_path = raw / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["checkpoint_sha256"] = closure.sha256_file(checkpoint_path)
    closure.write_json_atomic(manifest_path, manifest)


class DeckContractTests(unittest.TestCase):
    def test_canonical_deck_and_family_partition(self) -> None:
        contract = closure.validate_deck(DECK)
        self.assertEqual(contract.step_limit, 200_000)
        self.assertEqual(contract.report_every, 10_000)
        self.assertIn("numerical open-flow sink", contract.text)
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
        self.assertEqual(len(closure.REACTION_REGISTRY), 22)
        self.assertEqual(
            closure.REGISTRY_BY_NAME["desorb-al"].provenance_class, "heuristic"
        )
        self.assertEqual(
            closure.REGISTRY_BY_NAME["desorb-si"].observable_type,
            "mixed_desorption_proxy",
        )
        self.assertIn("electronic", closure.REGISTRY_BY_NAME["desorb-si"].rationale)

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

    def test_dilute_reservoir_sabotage_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            activity = _write_sabotage(root, ("Al = 1.0e-30", "Al = 1.0"))
            with self.assertRaisesRegex(ValueError, "open-flow sink"):
                closure.validate_deck(activity)
            chemical_potential = _write_sabotage(
                root, ("Al = -1.0\nSi = -1.0", "Al = 0.0\nSi = -1.0")
            )
            with self.assertRaisesRegex(ValueError, "far-from-equilibrium"):
                closure.validate_deck(chemical_potential)

    def test_provenance_records_exact_thermodynamic_boundary(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "provenance.csv"
            closure._write_provenance(path)
            with path.open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
        conditions = {
            row["constant_name"]: row
            for row in rows
            if row["record_type"] == "condition"
        }
        self.assertEqual(
            set(conditions),
            {
                "temperature",
                "dissolved_cation_activity",
                "dissolved_cation_mu",
                "effective_consumed_cation_factor_298k",
                "simulation_steps",
                "report_every_steps",
            },
        )
        self.assertEqual(conditions["temperature"]["constant_value"], "298.0")
        self.assertEqual(
            conditions["dissolved_cation_activity"]["constant_value"], "1e-30"
        )
        self.assertEqual(conditions["dissolved_cation_mu"]["constant_value"], "-1.0")
        self.assertIn(
            "numerical open-flow sink",
            conditions["dissolved_cation_activity"]["rationale"],
        )
        self.assertEqual(conditions["simulation_steps"]["constant_value"], "200000")
        self.assertEqual(conditions["report_every_steps"]["constant_value"], "10000")

    def test_canonical_run_length_and_cadence_are_required(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            short = _write_sabotage(root, ("steps = 200000", "steps = 20000"))
            with self.assertRaisesRegex(ValueError, "exactly 200000"):
                closure.validate_deck(short)
            sparse = _write_sabotage(
                root,
                (
                    "[observables]\nreport_every = 10000",
                    "[observables]\nreport_every = 20000",
                ),
            )
            with self.assertRaisesRegex(ValueError, "exactly 10000"):
                closure.validate_deck(sparse)
            simulation_cadence = _write_sabotage(
                root,
                (
                    "seed = 42\nreport_every = 10000",
                    "seed = 42\nreport_every = 20000",
                ),
            )
            with self.assertRaisesRegex(ValueError, "simulation.report_every"):
                closure.validate_deck(simulation_cadence)

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

    def test_non_eight_replicas_and_duplicate_seeds_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "exactly 8"):
            closure.validate_seeds(range(7))
        with self.assertRaisesRegex(ValueError, "exactly 8"):
            closure.validate_seeds(range(9))
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

    def test_registry_barrier_sabotage_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = _write_sabotage(
                Path(directory),
                (
                    "prefactor = 1.0e13, ea = 6.400",
                    "prefactor = 1.0e13, ea = 2.400",
                ),
            )
            with self.assertRaisesRegex(ValueError, "exactly match registry"):
                closure.validate_deck(path)


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
                1.0,
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
            "state_types": [
                state["occupant"]
                for kind in contract.parsed["kinds"]
                for state in kind["states"]
            ],
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
        self.assertEqual(events.steps, (1, 2, 3, 4))
        self.assertEqual(events.times[:2], (1.0, 1.0))
        self.assertEqual(
            closure.dissolution_event_counts(events, 0, 4),
            {
                "lattice_si": 1,
                "gross_si": 1,
                "adsorb_si": 1,
                "net_si": 0,
                "lattice_al": 1,
                "gross_al": 1,
                "adsorb_al": 1,
                "net_al": 0,
            },
        )

    def test_event_lineage_counts_only_original_lattice_cations(self) -> None:
        contract = closure.validate_deck(DECK)
        states = [
            f"{kind['name']}.{state['name']}"
            for kind in contract.parsed["kinds"]
            for state in kind["states"]
        ]
        state_id = {name: index for index, name in enumerate(states)}
        reactions = [reaction["name"] for reaction in contract.parsed["reactions"]]
        reaction_id = {name: index for index, name in enumerate(reactions)}
        header = {
            "petra_traj": 1,
            "deck": contract.parsed["deck"]["name"],
            "seed": 99,
            "n_sites": 10,
            "states": states,
            "state_types": [
                state["occupant"]
                for kind in contract.parsed["kinds"]
                for state in kind["states"]
            ],
            "reactions": reactions,
        }
        rows = [
            [
                1,
                1.0,
                reaction_id["R14-alohal-hydrolysis"],
                [
                    [8, state_id["Oaa.br"], state_id["Oaa.hy"]],
                    [4, state_id["Si.oh4"], state_id["Si.oh3"]],
                ],
            ],
            [
                2,
                2.0,
                reaction_id["R15-alohal-condensation"],
                [
                    [8, state_id["Oaa.hy"], state_id["Oaa.br"]],
                    [4, state_id["Si.oh3"], state_id["Si.oh4"]],
                ],
            ],
            [
                3,
                3.0,
                reaction_id["desorb-si"],
                [[4, state_id["Si.oh4"], state_id["Si.empty"]]],
            ],
            [
                4,
                4.0,
                reaction_id["adsorb-si"],
                [[4, state_id["Si.empty"], state_id["Si.oh4"]]],
            ],
            [
                5,
                5.0,
                reaction_id["desorb-si"],
                [[4, state_id["Si.oh4"], state_id["Si.empty"]]],
            ],
            [
                6,
                6.0,
                reaction_id["adsorb-si"],
                [[5, state_id["Si.empty"], state_id["Si.oh4"]]],
            ],
            [
                7,
                7.0,
                reaction_id["desorb-si"],
                [[5, state_id["Si.oh4"], state_id["Si.empty"]]],
            ],
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "events.jsonl"
            path.write_text(
                "\n".join([json.dumps(header), *(json.dumps(row) for row in rows)])
                + "\n",
                encoding="utf-8",
            )
            events = closure.parse_events(path, contract, expected_seed=99)

        self.assertEqual(
            closure.dissolution_event_counts(events, 0, 7),
            {
                "lattice_si": 1,
                "gross_si": 3,
                "adsorb_si": 2,
                "net_si": 1,
                "lattice_al": 0,
                "gross_al": 0,
                "adsorb_al": 0,
                "net_al": 0,
            },
        )
        self.assertFalse(closure.propensity_origin_safe(events))
        self.assertEqual(closure._lattice_releases_through_step(events, 7, "si"), 1)

    def test_event_step_regression_and_time_regression_are_rejected(self) -> None:
        contract = closure.validate_deck(DECK)
        states = [
            f"{kind['name']}.{state['name']}"
            for kind in contract.parsed["kinds"]
            for state in kind["states"]
        ]
        reactions = [reaction["name"] for reaction in contract.parsed["reactions"]]
        header = {
            "petra_traj": 1,
            "deck": contract.parsed["deck"]["name"],
            "seed": 1,
            "n_sites": 1,
            "states": states,
            "state_types": [
                state["occupant"]
                for kind in contract.parsed["kinds"]
                for state in kind["states"]
            ],
            "reactions": reactions,
        }
        state_id = {name: index for index, name in enumerate(states)}
        reaction_id = reactions.index("R14-alohal-hydrolysis")
        base = [
            1,
            1.0,
            reaction_id,
            [[0, state_id["Oaa.br"], state_id["Oaa.hy"]]],
        ]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "events.jsonl"
            for second, message in (
                ([1, 2.0, reaction_id, base[3]], "strictly increasing"),
                ([2, 0.5, reaction_id, base[3]], "nondecreasing"),
            ):
                path.write_text(
                    "\n".join(map(json.dumps, (header, base, second))) + "\n",
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(ValueError, message):
                    closure.parse_events(path, contract, expected_seed=1)

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

    def test_propensity_integral_uses_trapezoids_and_exact_species_mapping(
        self,
    ) -> None:
        names = tuple(entry.name for entry in closure.REACTION_REGISTRY)
        si_index = names.index("desorb-si")
        al_index = names.index("desorb-al")
        samples = []
        for step, sample_time, area, si, al in (
            (100, 0.0, 10.0, 2.0, 1.0),
            (200, 1.0, 20.0, 4.0, 1.0),
            (300, 3.0, 40.0, 8.0, 3.0),
        ):
            rates = [0.0] * len(names)
            rates[si_index] = si
            rates[al_index] = al
            samples.append(
                closure.ObservableSample(
                    step,
                    sample_time,
                    {"event_rates": rates, "surface_area": [area, 0.0, 0.0]},
                )
            )
        estimate = closure.integrate_expected_gross_dissolution_from_propensity(
            samples, names, [100, 200, 300]
        )
        self.assertEqual(estimate.area_time_a2_s, 75.0)
        self.assertEqual(estimate.expected_gross_si_events_from_propensity, 15.0)
        self.assertEqual(estimate.expected_gross_al_events_from_propensity, 5.0)
        self.assertAlmostEqual(
            estimate.expected_gross_si_flux_from_propensity_mol_m2_s,
            closure._amount_to_flux(15.0, 75.0),
        )
        self.assertAlmostEqual(
            estimate.expected_gross_al_flux_from_propensity_mol_m2_s,
            closure._amount_to_flux(5.0, 75.0),
        )

    def test_invalid_propensity_and_reaction_mapping_are_rejected(self) -> None:
        names = tuple(entry.name for entry in closure.REACTION_REGISTRY)

        def samples(value: float) -> list:
            result = []
            for step, sample_time in ((100, 0.0), (200, 1.0)):
                rates = [0.0] * len(names)
                rates[names.index("desorb-si")] = value
                result.append(
                    closure.ObservableSample(
                        step,
                        sample_time,
                        {"event_rates": rates, "surface_area": [10.0, 0.0, 0.0]},
                    )
                )
            return result

        for invalid in (-1.0, float("nan")):
            with (
                self.subTest(invalid=invalid),
                self.assertRaisesRegex(ValueError, "propensit.*finite and nonnegative"),
            ):
                closure.integrate_expected_gross_dissolution_from_propensity(
                    samples(invalid), names, [100, 200]
                )
        swapped = list(names)
        si_index = swapped.index("desorb-si")
        al_index = swapped.index("desorb-al")
        swapped[si_index], swapped[al_index] = swapped[al_index], swapped[si_index]
        with self.assertRaisesRegex(ValueError, "reaction index mapping"):
            closure.integrate_expected_gross_dissolution_from_propensity(
                samples(1.0), swapped, [100, 200]
            )
        with self.assertRaisesRegex(ValueError, "sample step alignment"):
            closure.integrate_expected_gross_dissolution_from_propensity(
                samples(1.0), names, [100, 300]
            )


class SteadyStateGateTests(unittest.TestCase):
    def test_constant_two_species_four_block_tail_passes(self) -> None:
        gate = closure.assess_steady_state(_points([2] * 20), 20_000)
        self.assertTrue(gate.acceptance_passed, gate.reasons)
        self.assertEqual(gate.status, "steady-positive")
        self.assertEqual(gate.si.status, "stationary-positive")
        self.assertEqual(gate.al.status, "stationary-positive")
        self.assertEqual(gate.si.lattice_release_events, 32)
        self.assertEqual(gate.al.lattice_release_events, 32)

    def test_zero_is_stationary_complete_accepted_with_poisson_bounds(self) -> None:
        gate = closure.assess_steady_state(_points([0] * 20), 20_000)
        self.assertEqual(gate.status, "steady-zero")
        self.assertEqual(gate.dissolution_outcome, "no-dissolution")
        self.assertEqual(gate.si.status, "stationary-zero")
        self.assertEqual(gate.al.status, "stationary-zero")
        self.assertTrue(gate.evidence_complete)
        self.assertTrue(gate.acceptance_passed)
        self.assertGreater(gate.si.upper_95_mol_m2_s, 0.0)
        self.assertGreater(gate.al.upper_95_mol_m2_s, 0.0)

    def test_mixed_zero_and_positive_response_is_unresolved(self) -> None:
        gate = closure.assess_steady_state(
            _points([1] * 4 + [0] * 4 + [1] * 12), 20_000
        )
        self.assertEqual(gate.status, "nonsteady")
        self.assertEqual(gate.si.status, "unresolved-mixed-response")
        self.assertFalse(gate.evidence_complete)
        self.assertFalse(gate.acceptance_passed)

    def test_stable_depleted_tail_does_not_require_initial_inventory(self) -> None:
        populations = [100, 75, 50, 50, *([50] * 17)]
        gate = closure.assess_steady_state(
            _points(
                [2] * 20,
                solid_si_samples=populations,
                solid_al_samples=populations,
            ),
            20_000,
        )
        self.assertTrue(gate.acceptance_passed, gate.reasons)
        self.assertEqual(gate.si_population.final_fraction, 0.5)
        self.assertEqual(gate.si_population.stability_status, "stable")

    def test_evolving_tail_population_fails_acceptance(self) -> None:
        populations = [100] * 4 + list(range(100, 83, -1))
        gate = closure.assess_steady_state(
            _points(
                [2] * 20,
                solid_si_samples=populations,
                solid_al_samples=populations,
            ),
            20_000,
        )
        self.assertEqual(gate.status, "nonsteady")
        self.assertEqual(gate.population_stability_status, "evolving")
        self.assertFalse(gate.acceptance_passed)

    def test_species_population_floor_is_not_hidden_by_combined_inventory(self) -> None:
        gate = closure.assess_steady_state(
            _points([2] * 20, solid_si_final=50, solid_al_final=5), 20_000
        )
        self.assertEqual(gate.status, "absorbed")
        self.assertEqual(gate.si_population.final_fraction, 1.0)
        self.assertEqual(gate.al_population.final_fraction, 0.1)

    def test_requires_twenty_one_samples_and_complete_run(self) -> None:
        self.assertEqual(
            closure.assess_steady_state(_points([2] * 20), 21_000).status,
            "incomplete",
        )
        self.assertEqual(
            closure.assess_steady_state(_points([2] * 19), 19_000).status,
            "incomplete",
        )
        self.assertEqual(
            closure.assess_steady_state(
                _points([2] * 20, area_final=0.0), 20_000
            ).status,
            "absorbed",
        )

    def test_four_block_monotonic_trend_fails(self) -> None:
        increments = [1] * 4 + [1] * 4 + [2] * 4 + [4] * 4 + [8] * 4
        gate = closure.assess_steady_state(_points(increments), 20_000)
        self.assertEqual(gate.status, "nonsteady")
        self.assertTrue(
            any("trend" in reason or "shift" in reason for reason in gate.reasons)
        )

    def test_evolving_sampled_propensity_fails_even_with_stationary_event_counts(
        self,
    ) -> None:
        propensities = [1.0] * 5 + [2.0] * 4 + [4.0] * 4 + [8.0] * 4 + [16.0] * 4
        gate = closure.assess_steady_state(
            _points(
                [2] * 20,
                si_propensities=propensities,
                al_propensities=propensities,
            ),
            20_000,
        )
        self.assertEqual(gate.status, "nonsteady")
        self.assertGreater(abs(gate.si.propensity_trend_decades), 0.35)
        self.assertTrue(any("propensity" in reason for reason in gate.reasons))

    def test_opposite_replica_trends_cannot_cancel_in_mean(self) -> None:
        increasing_increments = [1] * 4 + [1] * 4 + [2] * 4 + [4] * 4 + [8] * 4
        decreasing_increments = [8] * 4 + [8] * 4 + [7] * 4 + [5] * 4 + [1] * 4
        increasing = closure.assess_steady_state(_points(increasing_increments), 20_000)
        decreasing = closure.assess_steady_state(_points(decreasing_increments), 20_000)
        self.assertFalse(increasing.acceptance_passed)
        self.assertFalse(decreasing.acceptance_passed)
        self.assertEqual(
            [
                a + b
                for a, b in zip(
                    increasing_increments, decreasing_increments, strict=True
                )
            ],
            [9] * 20,
            "ensemble averaging would hide the opposed trends",
        )


class CampaignEndToEndTests(unittest.TestCase):
    def test_censored_families_have_undefined_rank_but_complete_classification(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw, derived = _build_synthetic_campaign(root)
            _contaminate_one_run_with_adsorption(raw, "adsorption__ea-minus-3")
            _make_one_run_propensity_nonstationary(raw, "cation-desorption__ea-minus-3")
            analysis = closure.analyze_campaign(raw, derived)
            self.assertTrue(analysis["sensitivity_complete"])
            self.assertFalse(analysis["ordinal_ranking_complete"])
            self.assertFalse(analysis["acceptance_passed"])
            verified = closure.verify_campaign(raw, derived)
            self.assertTrue(verified["sensitivity_complete"])
            self.assertFalse(verified["ordinal_ranking_complete"])
            with (derived / "sensitivity-ranking.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                rows = list(csv.DictReader(handle))
            censored = next(row for row in rows if row["family"] == "adsorption")
            self.assertEqual(censored["status"], "censored")
            self.assertEqual(censored["rank"], "undefined")
            self.assertEqual(censored["ea-minus-3_status"], "censored")
            nonstationary = next(
                row for row in rows if row["family"] == "cation-desorption"
            )
            self.assertEqual(nonstationary["status"], "censored")
            self.assertEqual(nonstationary["rank"], "undefined")
            estimated = [row for row in rows if row["status"] == "estimated"]
            self.assertEqual({int(row["rank"]) for row in estimated}, set(range(1, 6)))
            with (derived / "per-replica-rates.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                rates = list(csv.DictReader(handle))
            contaminated = next(
                row
                for row in rates
                if row["scenario"] == "adsorption__ea-minus-3" and row["replica"] == "0"
            )
            self.assertEqual(
                contaminated["propensity_origin_status"], "origin-contaminated"
            )
            self.assertEqual(
                contaminated[
                    "expected_lattice_origin_si_flux_from_propensity_mol_m2_s"
                ],
                "undefined",
            )

    def test_all_censored_campaign_is_not_card_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw, derived = _build_synthetic_campaign(root)
            _zero_all_desorption_propensities(raw)

            analysis = closure.analyze_campaign(raw, derived)

            self.assertEqual(analysis["campaign_outcome"], "no-dissolution")
            self.assertTrue(analysis["sensitivity_complete"])
            self.assertFalse(analysis["ordinal_ranking_complete"])
            self.assertFalse(analysis["acceptance_passed"])

    def test_zero_campaign_analyzes_verifies_and_tampering_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw, derived = _build_synthetic_campaign(root)
            analysis = closure.analyze_campaign(raw, derived)
            self.assertEqual(analysis["campaign_outcome"], "no-dissolution")
            self.assertTrue(analysis["acceptance_passed"])
            self.assertTrue(analysis["ordinal_ranking_complete"])
            self.assertEqual(
                analysis["status_counts"],
                {
                    "steady-positive": 0,
                    "steady-zero": 29 * 8,
                    "nonsteady": 0,
                    "absorbed": 0,
                    "incomplete": 0,
                },
            )
            verified = closure.verify_campaign(raw, derived)
            self.assertEqual(verified, analysis)
            self.assertEqual(
                analysis["propensity_estimator_basis"], "integrated_ctmc_hazard"
            )
            self.assertIn(
                "not observed event release", analysis["propensity_interpretation"]
            )
            with (derived / "per-replica-rates.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                replica_rates = list(csv.DictReader(handle))
            self.assertTrue(
                all(
                    row["steady_state_status"] == "steady-zero" for row in replica_rates
                )
            )
            self.assertTrue(
                all(row["acceptance_passed"] == "True" for row in replica_rates)
            )
            for row in replica_rates:
                self.assertEqual(row["lattice_si_release_events"], "0")
                self.assertEqual(row["lattice_al_release_events"], "0")
                self.assertEqual(row["propensity_origin_status"], "origin-safe")
                self.assertNotIn("gross_si_flux_mol_m2_s", row)
                self.assertGreater(
                    float(
                        row["expected_lattice_origin_si_flux_from_propensity_mol_m2_s"]
                    ),
                    0.0,
                )
                self.assertGreater(
                    float(
                        row["expected_lattice_origin_al_flux_from_propensity_mol_m2_s"]
                    ),
                    0.0,
                )
                self.assertNotEqual(
                    row[
                        "log10_expected_lattice_origin_si_flux_from_propensity_mol_m2_s"
                    ],
                    "undefined",
                )
            with (derived / "stoichiometry.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                stoichiometry = list(csv.DictReader(handle))
            self.assertIn(
                "si_al_lattice_release_ratio_dimensionless_mean", stoichiometry[0]
            )
            self.assertNotIn("si_al_net_ratio_dimensionless_mean", stoichiometry[0])
            with (derived / "ensemble-rates.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                ensemble = list(csv.DictReader(handle))
            propensity_rows = [
                row
                for row in ensemble
                if row["rate_basis"] == closure.PROPENSITY_RATE_BASIS
            ]
            self.assertEqual(len(propensity_rows), 29 * 2)
            self.assertTrue(
                all(row["status"] == "estimated-positive" for row in propensity_rows)
            )
            with (derived / "sensitivity-ranking.csv").open(
                newline="", encoding="utf-8"
            ) as handle:
                sensitivity = list(csv.DictReader(handle))
            self.assertEqual(len(sensitivity), 7)
            self.assertEqual(
                {row["family"] for row in sensitivity},
                set(closure.FAMILY_REACTIONS),
            )
            self.assertEqual(
                {row["declared_response"] for row in sensitivity},
                {closure.SENSITIVITY_RESPONSE},
            )
            self.assertEqual({row["status"] for row in sensitivity}, {"estimated"})
            self.assertEqual(
                {int(row["rank"]) for row in sensitivity}, set(range(1, 8))
            )
            for row in sensitivity:
                for perturbation in closure.PERTURBATIONS:
                    self.assertEqual(row[f"{perturbation}_status"], "estimated")
                    self.assertEqual(
                        int(row[f"{perturbation}_paired_replica_count"]),
                        closure.REPLICA_COUNT,
                    )
                    self.assertTrue(
                        math.isfinite(float(row[f"{perturbation}_paired_delta_log10"]))
                    )
                    self.assertTrue(
                        math.isfinite(
                            float(row[f"{perturbation}_paired_mean_delta_mol_m2_s"])
                        )
                    )

            changed = root / "derived-changed"
            shutil.copytree(derived, changed)
            with (changed / "per-replica-rates.csv").open(
                "a", encoding="utf-8"
            ) as handle:
                handle.write("forged\n")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                closure.verify_campaign(raw, changed)

            changed_propensity = root / "derived-changed-propensity"
            shutil.copytree(derived, changed_propensity)
            propensity_path = changed_propensity / "per-replica-rates.csv"
            with propensity_path.open(newline="", encoding="utf-8") as handle:
                reader = csv.DictReader(handle)
                fields = reader.fieldnames
                propensity_rows = list(reader)
            assert fields is not None
            field = "expected_lattice_origin_si_flux_from_propensity_mol_m2_s"
            changed_value = float(propensity_rows[0][field]) * 2.0
            propensity_rows[0][field] = str(changed_value)
            propensity_rows[0][
                "log10_expected_lattice_origin_si_flux_from_propensity_mol_m2_s"
            ] = str(math.log10(changed_value))
            closure._write_csv_atomic(propensity_path, fields, propensity_rows)
            propensity_receipt_path = changed_propensity / "verification.json"
            propensity_receipt = json.loads(
                propensity_receipt_path.read_text(encoding="utf-8")
            )
            propensity_receipt["output_sha256"]["per-replica-rates.csv"] = (
                closure.sha256_file(propensity_path)
            )
            closure.write_json_atomic(propensity_receipt_path, propensity_receipt)
            with self.assertRaisesRegex(ValueError, "not reproducible"):
                closure.verify_campaign(raw, changed_propensity)

            changed_response = root / "derived-changed-response"
            shutil.copytree(derived, changed_response)
            response_path = changed_response / "sensitivity-ranking.csv"
            response_text = response_path.read_text(encoding="utf-8")
            self.assertEqual(
                response_text.count(closure.SENSITIVITY_RESPONSE),
                len(closure.FAMILY_REACTIONS),
            )
            response_path.write_text(
                response_text.replace(closure.SENSITIVITY_RESPONSE, "forged_response"),
                encoding="utf-8",
            )
            response_receipt_path = changed_response / "verification.json"
            response_receipt = json.loads(
                response_receipt_path.read_text(encoding="utf-8")
            )
            response_receipt["output_sha256"]["sensitivity-ranking.csv"] = (
                closure.sha256_file(response_path)
            )
            closure.write_json_atomic(response_receipt_path, response_receipt)
            with self.assertRaisesRegex(ValueError, "declared response"):
                closure.verify_campaign(raw, changed_response)

            missing_family = root / "derived-missing-family"
            shutil.copytree(derived, missing_family)
            sensitivity_path = missing_family / "sensitivity-ranking.csv"
            lines = sensitivity_path.read_text(encoding="utf-8").splitlines()
            sensitivity_path.write_text("\n".join(lines[:-1]) + "\n", encoding="utf-8")
            receipt_path = missing_family / "verification.json"
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            receipt["output_sha256"]["sensitivity-ranking.csv"] = closure.sha256_file(
                sensitivity_path
            )
            closure.write_json_atomic(receipt_path, receipt)
            with self.assertRaisesRegex(ValueError, "row count mismatch"):
                closure.verify_campaign(raw, missing_family)

            manifest_path = raw / "manifest.json"
            manifest_bytes = manifest_path.read_bytes()
            manifest = json.loads(manifest_bytes)
            manifest["source_deck_sha256"] = ""
            closure.write_json_atomic(manifest_path, manifest)
            with self.assertRaisesRegex(ValueError, "source deck binding"):
                closure.verify_campaign(raw, derived)
            manifest_path.write_bytes(manifest_bytes)

            missing_raw = (
                raw
                / "runs"
                / "nominal"
                / "replica-00-seed-90401"
                / "snapshot.pgif.json"
            )
            missing_raw.unlink()
            with self.assertRaisesRegex(ValueError, "missing/empty run artifact"):
                closure.verify_campaign(raw, derived)


class DeterminismTests(unittest.TestCase):
    def test_derived_hash_inventory_ignores_appledouble_sidecars(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            expected = closure.DERIVED_FILES - {"verification.json"}
            for name in expected:
                (root / name).write_text(name, encoding="utf-8")
                (root / f"._{name}").write_text("appledouble", encoding="utf-8")

            hashes = closure._derived_output_hashes(root)

            self.assertEqual(set(hashes), expected)

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
