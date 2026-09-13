#!/usr/bin/env python3
"""Run and analyze the A9 approximate kaolinite rate-closure campaign.

This is survey-tier platform testing, not calibrated or production kinetics.

Steady-state criterion (identical for nominal and every sensitivity scenario):
for each replica, require at least eight cadence samples, completion of the
configured step limit, a strictly increasing finite clock, and a final solid
cation inventory and geometric area each above 10% of their initial values.
Over the final six cadence intervals, require at least 12 gross cation-
desorption events, at least four nonzero interval fluxes, an absolute fitted
end-to-end log10-flux trend <= 0.35 decade, and an absolute first-half versus
second-half log10 mean-flux change <= 0.35 decade.  The gate is evaluated per
replica, never on an ensemble average, so opposite trends cannot cancel.  The
inventory/area floors reject an absorbing fully dissolved finite slab.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import json
import math
import os
import random
import re
import statistics
import subprocess
import threading
import time
import tomllib
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

AVOGADRO_EXACT = 6.02214076e23
ANGSTROM2_TO_M2 = 1.0e-20
ARRHENIUS_PREFACTOR = 1.0e13
MIN_REPLICAS = 8
MAX_WORKERS = 4
MAX_TIMEOUT_SECONDS = 1_800
BOOTSTRAP_RESAMPLES = 2_000
DEFAULT_SEEDS = (90401, 90403, 90407, 90409, 90413, 90419, 90421, 90427)
REQUIRED_OBSERVABLES = frozenset(
    {"state_counts", "event_rates", "rate_spectra", "surface_area", "exposure_age"}
)
PROVENANCE_CLASSES = frozenset({"computed", "literature", "heuristic"})
PERTURBATIONS = ("ea-minus-3", "ea-plus-3", "prefactor-x0.1", "prefactor-x10")
SENSITIVITY_RESPONSE = "net_release_flux"

FAMILY_REACTIONS: dict[str, tuple[str, ...]] = {
    "siloxane-neutral": (
        "R0-sio-si-hydrolysis",
        "R1-sio-si-condensation",
    ),
    "sioal-si-neutral": (
        "R2-sioal2-si-hydrolysis",
        "R3-sioal2-si-condensation",
        "R10-sial-hydrolysis",
        "R11-sial-condensation",
    ),
    "sioal-al-neutral": (
        "R4-sioal2-al-hydrolysis",
        "R5a-sioal2-al-condensation",
        "R5b-sioal2-al-condensation",
        "R12-albr-hydrolysis",
        "R13-albr-condensation",
    ),
    "connectivity-ladder": (
        "R6-second-stage-hydrolysis",
        "R7-second-stage-condensation",
        "R8a-sioha-hydrolysis",
        "R8b-sioha-hydrolysis",
        "R9-sioha-condensation",
    ),
    "al-o-al-analogue": (
        "R14-alohal-hydrolysis",
        "R15-alohal-condensation",
    ),
    "adsorption": ("adsorb-al", "adsorb-si"),
    "cation-desorption": ("desorb-al", "desorb-si"),
}
EXPECTED_REACTIONS = tuple(
    name for names in FAMILY_REACTIONS.values() for name in names
)
ANNOTATION_RE = re.compile(
    r'^name = "(?P<name>[^"]+)"\n'
    r'# a9-family = "(?P<family>[^"]+)"\n'
    r'# a9-provenance = "(?P<provenance>[^"]+)"$',
    re.MULTILINE,
)
REACTION_BLOCK_RE = re.compile(
    r"(?ms)^\[\[reactions\]\]\n(?P<body>.*?)(?=^\[\[reactions\]\]\n|\Z)"
)
RATE_RE = re.compile(
    r"rate = \{ arrhenius = \{ prefactor = (?P<prefactor>[^,]+), ea = (?P<ea>[^}]+) \} \}"
)


@dataclass(frozen=True)
class DeckContract:
    path: Path
    text: str
    parsed: dict
    annotations: dict[str, tuple[str, str]]
    barriers: dict[str, float]
    prefactors: dict[str, float]
    step_limit: int
    report_every: int


@dataclass(frozen=True)
class Scenario:
    name: str
    family: str | None
    perturbation: str | None


@dataclass(frozen=True)
class EventData:
    seed: int
    reaction_names: tuple[str, ...]
    times: tuple[float, ...]
    names: tuple[str, ...]
    counts: dict[str, int]


@dataclass(frozen=True)
class PopulationRow:
    step: int
    time: float
    states: dict[str, int]


@dataclass(frozen=True)
class SteadyPoint:
    step: int
    time: float
    gross_events: int
    area_a2: float
    solid_cations: int


@dataclass(frozen=True)
class SteadyStateGate:
    passed: bool
    reasons: tuple[str, ...]
    window_start_step: int | None
    window_end_step: int | None
    gross_events: int
    positive_intervals: int
    trend_decades: float | None
    half_change_decades: float | None


@dataclass(frozen=True)
class ReplicaRate:
    scenario: str
    replica: int
    seed: int
    window_start_step: int
    window_end_step: int
    window_start_time_s: float
    window_end_time_s: float
    area_time_a2_s: float
    gross_si_events: int
    adsorb_si_events: int
    net_si_events: int
    gross_al_events: int
    adsorb_al_events: int
    net_al_events: int
    gross_si_flux: float
    net_si_flux: float
    gross_al_flux: float
    net_al_flux: float
    si_al_net_ratio: float | None
    gate_passed: bool


def _finite_number(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be numeric")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json_atomic(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    payload = json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    with temporary.open("w", encoding="utf-8") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def _write_csv_atomic(
    path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def parse_seeds(text: str) -> tuple[int, ...]:
    try:
        seeds = tuple(int(item) for item in text.split(",") if item != "")
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "seeds must be comma-separated integers"
        ) from exc
    try:
        validate_seeds(seeds)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    return seeds


def validate_seeds(seeds: Sequence[int]) -> tuple[int, ...]:
    normalized = tuple(seeds)
    if len(normalized) < MIN_REPLICAS:
        raise ValueError(f"campaign requires at least {MIN_REPLICAS} replicas")
    if any(type(seed) is not int or seed < 0 for seed in normalized):
        raise ValueError("every seed must be a non-negative integer")
    if len(set(normalized)) != len(normalized):
        raise ValueError("campaign seeds must be unique (duplicate seed found)")
    return normalized


def _schedule_present(value: object) -> bool:
    if isinstance(value, dict):
        for key, child in value.items():
            if key == "schedule" and child not in (None, [], {}):
                return True
            if _schedule_present(child):
                return True
    elif isinstance(value, list):
        return any(_schedule_present(child) for child in value)
    return False


def validate_deck(path: Path, seeds: Sequence[int] = DEFAULT_SEEDS) -> DeckContract:
    validate_seeds(seeds)
    text = path.read_text(encoding="utf-8")
    parsed = tomllib.loads(text)
    deck_meta = parsed.get("deck", {})
    if deck_meta.get("units") != "kcal/mol":
        raise ValueError("deck energy units must be explicit kcal/mol")
    temperature = parsed.get("thermo", {}).get("temperature")
    if type(temperature) is not float or temperature != 298.0:
        raise ValueError("thermo.temperature must be exactly the TOML float 298.0")
    if _schedule_present(parsed) or re.search(r"(?m)^\s*\[+[^]]*schedule", text):
        raise ValueError("execution schedules are forbidden for A9")

    reactions = parsed.get("reactions")
    if not isinstance(reactions, list):
        raise ValueError("deck must declare reactions")
    names = [reaction.get("name") for reaction in reactions]
    if len(names) != len(set(names)):
        raise ValueError("reaction names must be unique")
    if len(names) != 22 or set(names) != set(EXPECTED_REACTIONS):
        missing = sorted(set(EXPECTED_REACTIONS) - set(names))
        extra = sorted(set(names) - set(EXPECTED_REACTIONS))
        raise ValueError(
            f"reaction coverage must be exactly 22 names; missing={missing}, extra={extra}"
        )

    matches = list(ANNOTATION_RE.finditer(text))
    annotations = {
        match.group("name"): (match.group("family"), match.group("provenance"))
        for match in matches
    }
    if len(matches) != 22 or set(annotations) != set(names):
        raise ValueError(
            "every reaction needs adjacent a9-family/a9-provenance comments"
        )
    if any(
        provenance not in PROVENANCE_CLASSES for _, provenance in annotations.values()
    ):
        raise ValueError("provenance must be computed, literature, or heuristic")
    for family, expected in FAMILY_REACTIONS.items():
        actual = {
            name for name, annotation in annotations.items() if annotation[0] == family
        }
        if actual != set(expected):
            raise ValueError(
                f"family {family!r} is incomplete or contains unrelated reactions"
            )
    unknown_families = {family for family, _ in annotations.values()} - set(
        FAMILY_REACTIONS
    )
    if unknown_families:
        raise ValueError(f"unknown barrier families: {sorted(unknown_families)}")

    barriers: dict[str, float] = {}
    prefactors: dict[str, float] = {}
    for reaction in reactions:
        name = reaction["name"]
        try:
            arrhenius = reaction["rate"]["arrhenius"]
            prefactor = _finite_number(arrhenius["prefactor"], f"{name} prefactor")
            barrier = _finite_number(arrhenius["ea"], f"{name} barrier")
        except (KeyError, TypeError) as exc:
            raise ValueError(f"{name} must have an Arrhenius rate") from exc
        if prefactor != ARRHENIUS_PREFACTOR:
            raise ValueError(f"{name} prefactor must be exactly 1e13 s^-1")
        if barrier < 0.0:
            raise ValueError(f"{name} barrier must be nonnegative")
        barriers[name] = barrier
        prefactors[name] = prefactor

    observables = parsed.get("observables", {}).get("series", [])
    observed = {item.get("kind") for item in observables if isinstance(item, dict)}
    missing_observables = REQUIRED_OBSERVABLES - observed
    if missing_observables:
        raise ValueError(f"missing required observables: {sorted(missing_observables)}")
    simulation = parsed.get("simulation", {})
    steps = simulation.get("steps")
    report_every = parsed.get("observables", {}).get("report_every")
    if type(steps) is not int or steps <= 0:
        raise ValueError("simulation.steps must be a positive integer")
    if type(report_every) is not int or report_every <= 0:
        raise ValueError("observables.report_every must be a positive integer")
    if steps // report_every < 7:
        raise ValueError("deck cadence must provide at least eight samples")
    return DeckContract(
        path=path,
        text=text,
        parsed=parsed,
        annotations=annotations,
        barriers=barriers,
        prefactors=prefactors,
        step_limit=steps,
        report_every=report_every,
    )


def scenarios() -> tuple[Scenario, ...]:
    values = [Scenario("nominal", None, None)]
    for family in FAMILY_REACTIONS:
        values.extend(
            Scenario(f"{family}__{perturbation}", family, perturbation)
            for perturbation in PERTURBATIONS
        )
    return tuple(values)


def perturb_deck(contract: DeckContract, family: str, perturbation: str) -> str:
    if family not in FAMILY_REACTIONS or perturbation not in PERTURBATIONS:
        raise ValueError("unknown sensitivity family or perturbation")
    targets = set(FAMILY_REACTIONS[family])
    changed: set[str] = set()

    def replace_block(match: re.Match[str]) -> str:
        block = match.group(0)
        name_match = re.search(r'(?m)^name = "([^"]+)"$', block)
        if name_match is None or name_match.group(1) not in targets:
            return block
        name = name_match.group(1)
        rate_match = RATE_RE.search(block)
        if rate_match is None:
            raise ValueError(f"{name}: no canonical Arrhenius rate line")
        prefactor = float(rate_match.group("prefactor"))
        barrier = float(rate_match.group("ea"))
        if perturbation == "ea-minus-3":
            barrier -= 3.0
            if barrier < 0.0:
                raise ValueError(f"{family}: Ea -3 would make {name} negative")
        elif perturbation == "ea-plus-3":
            barrier += 3.0
        elif perturbation == "prefactor-x0.1":
            prefactor *= 0.1
        else:
            prefactor *= 10.0
        replacement = (
            "rate = { arrhenius = { "
            f"prefactor = {prefactor:.12g}, ea = {barrier:.12g}"
            " } }"
        )
        changed.add(name)
        return block[: rate_match.start()] + replacement + block[rate_match.end() :]

    result = REACTION_BLOCK_RE.sub(replace_block, contract.text)
    if changed != targets:
        raise ValueError(
            f"sensitivity did not change exactly family {family}: {sorted(changed)}"
        )
    _validate_sensitivity_isolation(contract, result, family, perturbation)
    return result


def _validate_sensitivity_isolation(
    nominal: DeckContract, candidate_text: str, family: str, perturbation: str
) -> None:
    candidate = tomllib.loads(candidate_text)
    candidate_reactions = {
        reaction["name"]: reaction for reaction in candidate["reactions"]
    }
    nominal_reactions = {
        reaction["name"]: reaction for reaction in nominal.parsed["reactions"]
    }
    targets = set(FAMILY_REACTIONS[family])
    for name in EXPECTED_REACTIONS:
        before = nominal_reactions[name]
        after = candidate_reactions[name]
        if name not in targets:
            if after != before:
                raise ValueError(f"sensitivity changed unrelated reaction {name}")
            continue
        before_copy = dict(before)
        after_copy = dict(after)
        before_rate = before_copy.pop("rate")
        after_rate = after_copy.pop("rate")
        if before_copy != after_copy:
            raise ValueError(f"sensitivity changed non-rate semantics for {name}")
        b = before_rate["arrhenius"]
        a = after_rate["arrhenius"]
        expected_prefactor = b["prefactor"]
        expected_ea = b["ea"]
        if perturbation == "ea-minus-3":
            expected_ea -= 3.0
        elif perturbation == "ea-plus-3":
            expected_ea += 3.0
        elif perturbation == "prefactor-x0.1":
            expected_prefactor *= 0.1
        else:
            expected_prefactor *= 10.0
        if (
            set(a) != {"prefactor", "ea"}
            or not math.isclose(a["prefactor"], expected_prefactor, rel_tol=1e-14)
            or not math.isclose(a["ea"], expected_ea, rel_tol=1e-14, abs_tol=1e-14)
        ):
            raise ValueError(f"incorrect sensitivity transform for {name}")


def _run_one(
    petra_bin: Path,
    deck: Path,
    output: Path,
    log_path: Path,
    seed: int,
    timeout: int,
    cwd: Path,
    env: dict[str, str],
) -> dict:
    if output.exists() or log_path.exists():
        raise ValueError(f"refusing to overwrite replica output/log for seed {seed}")
    command = [
        "nice",
        "-n",
        "10",
        str(petra_bin),
        str(deck),
        "--seed",
        str(seed),
        "--ensemble",
        "1",
        "--out",
        str(output),
        "--viz",
        "--paranoid",
    ]
    started = time.monotonic()
    with log_path.open("x", encoding="utf-8") as log:
        # argv sequence, shell=False; the executable and deck are explicit local paths.
        subprocess.run(
            command,
            cwd=cwd,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
            timeout=timeout,
            shell=False,
        )
        log.flush()
        os.fsync(log.fileno())
    elapsed = time.monotonic() - started
    required = (
        "events.jsonl",
        "populations.csv",
        "observables.csv",
        "snapshot.pgif.json",
    )
    missing = [name for name in required if not (output / name).is_file()]
    if missing:
        raise RuntimeError(f"seed {seed} missing Petra outputs: {missing}")
    hashes = {name: sha256_file(output / name) for name in required}
    hashes["log"] = sha256_file(log_path)
    return {
        "seed": seed,
        "elapsed_seconds": elapsed,
        "command": command,
        "output": str(output),
        "log": str(log_path),
        "sha256": hashes,
    }


def run_campaign(
    deck_path: Path,
    petra_bin: Path,
    raw_root: Path,
    seeds: Sequence[int] = DEFAULT_SEEDS,
    workers: int = MAX_WORKERS,
    timeout: int = MAX_TIMEOUT_SECONDS,
) -> None:
    seeds = validate_seeds(seeds)
    deck_path = deck_path.resolve()
    petra_bin = petra_bin.resolve()
    raw_root = raw_root.resolve()
    if not 1 <= workers <= MAX_WORKERS:
        raise ValueError(f"workers must be between 1 and {MAX_WORKERS}")
    if not 1 <= timeout <= MAX_TIMEOUT_SECONDS:
        raise ValueError(f"timeout must be in [1,{MAX_TIMEOUT_SECONDS}]")
    if raw_root.exists():
        raise ValueError(f"refusing to overwrite campaign root: {raw_root}")
    if not petra_bin.is_file():
        raise ValueError(f"Petra executable not found: {petra_bin}")
    contract = validate_deck(deck_path, seeds)

    decks_dir = raw_root / "decks"
    logs_dir = raw_root / "logs"
    runs_dir = raw_root / "runs"
    decks_dir.mkdir(parents=True)
    logs_dir.mkdir()
    runs_dir.mkdir()
    scenario_records = []
    for scenario in scenarios():
        text = (
            contract.text
            if scenario.family is None
            else perturb_deck(contract, scenario.family, scenario.perturbation or "")
        )
        scenario_deck = decks_dir / f"{scenario.name}.toml"
        scenario_deck.write_text(text, encoding="utf-8")
        (logs_dir / scenario.name).mkdir()
        (runs_dir / scenario.name).mkdir()
        scenario_records.append(
            {
                **asdict(scenario),
                "deck": str(scenario_deck),
                "deck_sha256": sha256_file(scenario_deck),
            }
        )
    manifest = {
        "schema": 1,
        "status": "running",
        "survey_tier": True,
        "temperature_k": 298.0,
        "units": "kcal/mol",
        "seeds": list(seeds),
        "replicas": len(seeds),
        "workers": workers,
        "timeout_seconds": timeout,
        "source_deck": str(deck_path.resolve()),
        "source_deck_sha256": sha256_file(deck_path),
        "scenarios": scenario_records,
    }
    write_json_atomic(raw_root / "manifest.json", manifest)
    receipts: list[dict] = []
    write_json_atomic(
        raw_root / "checkpoint.json", {"status": "running", "receipts": receipts}
    )

    env = os.environ.copy()
    for variable in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "RAYON_NUM_THREADS",
    ):
        env[variable] = "16"
    lock = threading.Lock()
    try:
        for scenario_record in scenario_records:
            scenario_name = scenario_record["name"]
            scenario_deck = Path(scenario_record["deck"])
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = {
                    pool.submit(
                        _run_one,
                        petra_bin.resolve(),
                        scenario_deck.resolve(),
                        runs_dir / scenario_name / f"replica-{index:02d}-seed-{seed}",
                        logs_dir
                        / scenario_name
                        / f"replica-{index:02d}-seed-{seed}.log",
                        seed,
                        timeout,
                        deck_path.parent.parent.resolve(),
                        env,
                    ): (index, seed)
                    for index, seed in enumerate(seeds)
                }
                for future in as_completed(futures):
                    index, seed = futures[future]
                    receipt = future.result()
                    receipt.update({"scenario": scenario_name, "replica": index})
                    with lock:
                        receipts.append(receipt)
                        receipts.sort(
                            key=lambda item: (item["scenario"], item["replica"])
                        )
                        write_json_atomic(
                            raw_root / "checkpoint.json",
                            {"status": "running", "receipts": receipts},
                        )
    except BaseException:
        write_json_atomic(
            raw_root / "checkpoint.json", {"status": "failed", "receipts": receipts}
        )
        raise
    manifest["status"] = "complete"
    manifest["completed_runs"] = len(receipts)
    manifest["checkpoint_sha256"] = sha256_file(raw_root / "checkpoint.json")
    write_json_atomic(raw_root / "manifest.json", manifest)
    write_json_atomic(
        raw_root / "checkpoint.json", {"status": "complete", "receipts": receipts}
    )
    # Rewrite after the final checkpoint so the manifest binds its final bytes.
    manifest["checkpoint_sha256"] = sha256_file(raw_root / "checkpoint.json")
    write_json_atomic(raw_root / "manifest.json", manifest)


def _qualified_center_states(reaction: dict) -> set[str]:
    kind = reaction["center"]["kind"]
    states = reaction["center"]["state"]
    return {
        state if "." in state else f"{kind}.{state}"
        for state in states
        if not state.startswith("@")
    }


def _center_destinations(reaction: dict) -> set[str]:
    kind = reaction["center"]["kind"]
    effects = list(reaction.get("effects", []))
    for branch in reaction.get("branches", []):
        effects.extend(branch.get("effects", []))
    destinations = set()
    for effect in effects:
        if effect.get("target") == "center" and "set" in effect:
            state = effect["set"]
            destinations.add(state if "." in state else f"{kind}.{state}")
    return destinations


def parse_events(
    path: Path, contract: DeckContract, expected_seed: int | None = None
) -> EventData:
    with path.open(encoding="utf-8") as handle:
        first = handle.readline()
        if not first:
            raise ValueError(f"{path}: empty event log")
        header = json.loads(first)
        reaction_names = header.get("reactions")
        states = header.get("states")
        seed = header.get("seed")
        n_sites = header.get("n_sites")
        if reaction_names != [
            reaction["name"] for reaction in contract.parsed["reactions"]
        ]:
            raise ValueError(f"{path}: reaction table does not match deck")
        expected_states = [
            f"{kind['name']}.{state['name']}"
            for kind in contract.parsed["kinds"]
            for state in kind["states"]
        ]
        if states != expected_states:
            raise ValueError(f"{path}: state table does not match deck")
        if (
            header.get("petra_traj") != 1
            or header.get("deck") != contract.parsed["deck"]["name"]
        ):
            raise ValueError(f"{path}: trajectory header does not match deck")
        if type(seed) is not int or (
            expected_seed is not None and seed != expected_seed
        ):
            raise ValueError(f"{path}: event seed mismatch")
        if type(n_sites) is not int or n_sites <= 0:
            raise ValueError(f"{path}: invalid site count")
        state_ids = {name: index for index, name in enumerate(states)}
        reaction_by_name = {
            reaction["name"]: reaction for reaction in contract.parsed["reactions"]
        }
        transitions: dict[str, tuple[set[int], set[int]]] = {}
        for name, reaction in reaction_by_name.items():
            old_names = _qualified_center_states(reaction)
            new_names = _center_destinations(reaction)
            if (
                old_names
                and new_names
                and old_names <= state_ids.keys()
                and new_names <= state_ids.keys()
            ):
                transitions[name] = (
                    {state_ids[value] for value in old_names},
                    {state_ids[value] for value in new_names},
                )

        times: list[float] = []
        names: list[str] = []
        counts = {name: 0 for name in reaction_names}
        previous_time = -math.inf
        previous_step = 0
        last_seen_state: dict[int, int] = {}
        for line_number, line in enumerate(handle, start=2):
            row = json.loads(line)
            if not isinstance(row, list) or len(row) != 4:
                raise ValueError(f"{path}:{line_number}: invalid event row")
            step, event_time, reaction_id, changes = row
            if type(step) is not int or step != previous_step + 1:
                raise ValueError(
                    f"{path}:{line_number}: event steps must be contiguous"
                )
            event_time = _finite_number(event_time, f"{path}:{line_number} time")
            if event_time <= previous_time:
                raise ValueError(f"{path}:{line_number}: event time must increase")
            if type(reaction_id) is not int or not 0 <= reaction_id < len(
                reaction_names
            ):
                raise ValueError(f"{path}:{line_number}: invalid reaction id")
            if not isinstance(changes, list) or not changes:
                raise ValueError(f"{path}:{line_number}: event has no state changes")
            normalized_changes = []
            for change in changes:
                if not isinstance(change, list) or len(change) != 3:
                    raise ValueError(f"{path}:{line_number}: invalid state change")
                site, old, new = change
                if (
                    type(site) is not int
                    or not 0 <= site < n_sites
                    or type(old) is not int
                    or type(new) is not int
                    or not 0 <= old < len(states)
                    or not 0 <= new < len(states)
                    or old == new
                ):
                    raise ValueError(f"{path}:{line_number}: invalid state transition")
                if site in last_seen_state and last_seen_state[site] != old:
                    raise ValueError(
                        f"{path}:{line_number}: state-change replay mismatch"
                    )
                normalized_changes.append((site, old, new))
            name = reaction_names[reaction_id]
            if name in transitions:
                old_ids, new_ids = transitions[name]
                centers = [
                    change
                    for change in normalized_changes
                    if change[1] in old_ids and change[2] in new_ids
                ]
                if len(centers) != 1:
                    raise ValueError(
                        f"{path}:{line_number}: expected center transition missing for {name}"
                    )
            for site, _, new in normalized_changes:
                last_seen_state[site] = new
            previous_step = step
            previous_time = event_time
            times.append(event_time)
            names.append(name)
            counts[name] += 1
    return EventData(seed, tuple(reaction_names), tuple(times), tuple(names), counts)


def parse_populations(
    path: Path, expected_states: Sequence[str]
) -> list[PopulationRow]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != ["step", "time", *expected_states]:
            raise ValueError(f"{path}: population columns do not match deck states")
        result = []
        previous_step = -1
        previous_time = -math.inf
        for line_number, row in enumerate(reader, start=2):
            try:
                step = int(row["step"])
                sample_time = _finite_number(
                    float(row["time"]), f"{path}:{line_number} time"
                )
                states = {state: int(row[state]) for state in expected_states}
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"{path}:{line_number}: invalid population row"
                ) from exc
            if step <= previous_step or sample_time < previous_time:
                raise ValueError(
                    f"{path}:{line_number}: nonmonotonic population sample"
                )
            if any(value < 0 for value in states.values()):
                raise ValueError(f"{path}:{line_number}: negative population")
            result.append(PopulationRow(step, sample_time, states))
            previous_step, previous_time = step, sample_time
    if not result:
        raise ValueError(f"{path}: no population samples")
    return result


def parse_observables(
    path: Path, expected_seed: int | None = None
) -> dict[tuple[int, float], dict[str, list[float]]]:
    grouped: dict[tuple[int, float], dict[str, dict[int, float]]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required_columns = {"replica", "seed", "step", "time", "kind", "index", "value"}
        if reader.fieldnames is None or set(reader.fieldnames) != required_columns:
            raise ValueError(f"{path}: invalid observables columns")
        for line_number, row in enumerate(reader, start=2):
            try:
                step = int(row["step"])
                sample_time = _finite_number(
                    float(row["time"]), f"{path}:{line_number} time"
                )
                index = int(row["index"])
                value = _finite_number(
                    float(row["value"]), f"{path}:{line_number} value"
                )
                replica = int(row["replica"])
                seed = int(row["seed"])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"{path}:{line_number}: invalid observable row"
                ) from exc
            if (
                replica != 0
                or (expected_seed is not None and seed != expected_seed)
                or step < 0
                or sample_time < 0
                or index < 0
            ):
                raise ValueError(
                    f"{path}:{line_number}: invalid direct-run observable identity"
                )
            kind = row["kind"]
            indices = grouped.setdefault((step, sample_time), {}).setdefault(kind, {})
            if index in indices:
                raise ValueError(f"{path}:{line_number}: duplicate observable index")
            indices[index] = value
    if not grouped:
        raise ValueError(f"{path}: no observable samples")
    result = {}
    previous = (-1, -math.inf)
    for key in sorted(grouped):
        if key[0] <= previous[0] or key[1] < previous[1]:
            raise ValueError(f"{path}: nonmonotonic observable cadence")
        kinds = grouped[key]
        missing = REQUIRED_OBSERVABLES - kinds.keys()
        if missing:
            raise ValueError(
                f"{path}: sample {key} missing observables {sorted(missing)}"
            )
        result[key] = {}
        for kind, indexed in kinds.items():
            if set(indexed) != set(range(len(indexed))):
                raise ValueError(f"{path}: sparse {kind} indices at sample {key}")
            result[key][kind] = [indexed[index] for index in range(len(indexed))]
        if not result[key]["surface_area"] or result[key]["surface_area"][0] < 0:
            raise ValueError(f"{path}: invalid geometric surface area")
        previous = key
    return result


def integrate_area(times: Sequence[float], areas: Sequence[float]) -> float:
    if len(times) != len(areas) or len(times) < 2:
        raise ValueError("area integration requires at least two paired samples")
    integral = 0.0
    for left, right, area_left, area_right in zip(times, times[1:], areas, areas[1:]):
        left = _finite_number(left, "area time")
        right = _finite_number(right, "area time")
        area_left = _finite_number(area_left, "geometric area")
        area_right = _finite_number(area_right, "geometric area")
        if right <= left or area_left < 0.0 or area_right < 0.0:
            raise ValueError(
                "area samples require increasing time and nonnegative area"
            )
        integral += 0.5 * (area_left + area_right) * (right - left)
    if not math.isfinite(integral) or integral <= 0.0:
        raise ValueError("integrated geometric area must be positive and finite")
    return integral


def event_count_to_flux(events: int, area_time_a2_s: float) -> float:
    if type(events) is not int:
        raise ValueError("event count must be an integer")
    area_time = _finite_number(area_time_a2_s, "area-time denominator")
    if area_time <= 0.0:
        raise ValueError("area-time denominator must be positive")
    return (events / AVOGADRO_EXACT) / (area_time * ANGSTROM2_TO_M2)


def _linear_trend(values: Sequence[float]) -> float:
    xs = list(range(len(values)))
    center_x = statistics.fmean(xs)
    center_y = statistics.fmean(values)
    denominator = sum((x - center_x) ** 2 for x in xs)
    slope = (
        sum((x - center_x) * (y - center_y) for x, y in zip(xs, values, strict=True))
        / denominator
    )
    return slope * (len(values) - 1)


def assess_steady_state(
    points: Sequence[SteadyPoint], expected_steps: int
) -> SteadyStateGate:
    reasons: list[str] = []
    if len(points) < 8:
        return SteadyStateGate(
            False,
            ("insufficient cadence: need at least eight samples",),
            None,
            None,
            0,
            0,
            None,
            None,
        )
    if points[-1].step < expected_steps:
        reasons.append("run stopped early before configured step limit")
    for left, right in zip(points, points[1:]):
        if (
            right.step <= left.step
            or not math.isfinite(right.time)
            or right.time <= left.time
        ):
            reasons.append("sample steps/times are not strictly increasing and finite")
            break
    initial = points[0]
    final = points[-1]
    if (
        initial.solid_cations <= 0
        or final.solid_cations <= 0.10 * initial.solid_cations
    ):
        reasons.append("solid inventory reached the absorbing/dissolved floor")
    if initial.area_a2 <= 0.0 or final.area_a2 <= 0.10 * initial.area_a2:
        reasons.append("geometric area reached the absorbing/dissolved floor")

    tail = points[-7:]
    interval_fluxes = []
    interval_events = []
    for left, right in zip(tail, tail[1:]):
        count = right.gross_events - left.gross_events
        interval_events.append(count)
        try:
            denominator = integrate_area(
                [left.time, right.time], [left.area_a2, right.area_a2]
            )
        except ValueError:
            denominator = math.nan
        interval_fluxes.append(
            count / denominator if count > 0 and denominator > 0 else 0.0
        )
    gross_events = sum(interval_events)
    positive = [
        value for value in interval_fluxes if value > 0.0 and math.isfinite(value)
    ]
    if gross_events < 12:
        reasons.append("too few gross dissolution events in the steady-state window")
    if len(positive) < 4:
        reasons.append("insufficient nonzero dissolution-event cadence")
    trend = None
    half_change = None
    if len(positive) >= 4:
        logs = [math.log10(value) for value in positive]
        trend = _linear_trend(logs)
        if abs(trend) > 0.35:
            reasons.append("per-replica fitted dissolution trend exceeds 0.35 decade")
        first = statistics.fmean(positive[: len(positive) // 2])
        second = statistics.fmean(positive[len(positive) // 2 :])
        half_change = math.log10(second / first)
        if abs(half_change) > 0.35:
            reasons.append(
                "per-replica half-window dissolution shift exceeds 0.35 decade"
            )
    return SteadyStateGate(
        not reasons,
        tuple(reasons),
        tail[0].step,
        tail[-1].step,
        gross_events,
        len(positive),
        trend,
        half_change,
    )


def _counts_through(events: EventData, time_s: float, names: set[str]) -> int:
    end = bisect.bisect_right(events.times, time_s)
    return sum(1 for name in events.names[:end] if name in names)


def _counts_between(events: EventData, start: float, end: float, name: str) -> int:
    left = bisect.bisect_right(events.times, start)
    right = bisect.bisect_right(events.times, end)
    return sum(1 for value in events.names[left:right] if value == name)


def dissolution_event_counts(
    events: EventData, start: float, end: float
) -> dict[str, int]:
    """Return gross desorption, adsorption, and net release for Si and Al."""

    if not math.isfinite(start) or not math.isfinite(end) or end <= start:
        raise ValueError("event-accounting window must be finite and increasing")
    result: dict[str, int] = {}
    for species in ("si", "al"):
        gross = _counts_between(events, start, end, f"desorb-{species}")
        adsorbed = _counts_between(events, start, end, f"adsorb-{species}")
        result[f"gross_{species}"] = gross
        result[f"adsorb_{species}"] = adsorbed
        result[f"net_{species}"] = gross - adsorbed
    return result


def _solid_cations(row: PopulationRow) -> int:
    return sum(
        count
        for state, count in row.states.items()
        if (state.startswith("Si.") or state.startswith("Al."))
        and not state.endswith(".empty")
    )


def _analyze_replica(
    scenario: str,
    replica: int,
    seed: int,
    run_dir: Path,
    contract: DeckContract,
) -> tuple[ReplicaRate, SteadyStateGate, list[PopulationRow]]:
    state_names = [
        state
        for kind in contract.parsed["kinds"]
        for state in (f"{kind['name']}.{entry['name']}" for entry in kind["states"])
    ]
    events = parse_events(run_dir / "events.jsonl", contract, seed)
    populations = parse_populations(run_dir / "populations.csv", state_names)
    observables = parse_observables(run_dir / "observables.csv", seed)
    if populations[-1].step != len(events.times):
        raise ValueError(f"{run_dir}: event count disagrees with final population step")
    expected_sample_steps = list(
        range(0, populations[-1].step + 1, contract.report_every)
    )
    if not expected_sample_steps or expected_sample_steps[-1] != populations[-1].step:
        expected_sample_steps.append(populations[-1].step)
    if [row.step for row in populations] != expected_sample_steps:
        raise ValueError(f"{run_dir}: population cadence does not match deck")
    observable_keys = list(observables)
    population_keys = [(row.step, row.time) for row in populations]
    if len(observable_keys) != len(population_keys):
        raise ValueError(f"{run_dir}: population/observable cadence length mismatch")
    aligned = []
    for row, key in zip(populations, observable_keys, strict=True):
        if row.step != key[0] or not math.isclose(
            row.time, key[1], rel_tol=2e-6, abs_tol=1e-15
        ):
            raise ValueError(f"{run_dir}: population/observable cadence mismatch")
        sample = observables[key]
        state_counts = sample["state_counts"]
        if len(state_counts) != len(state_names) or any(
            int(value) != row.states[name]
            for name, value in zip(state_names, state_counts, strict=True)
        ):
            raise ValueError(
                f"{run_dir}: state-count observable disagrees with populations"
            )
        aligned.append((row, sample["surface_area"][0]))
    gross_names = {"desorb-si", "desorb-al"}
    points = [
        SteadyPoint(
            row.step,
            row.time,
            _counts_through(events, row.time, gross_names),
            area,
            _solid_cations(row),
        )
        for row, area in aligned
    ]
    gate = assess_steady_state(points, contract.step_limit)
    start_index = len(points) - 7
    start, end = points[start_index], points[-1]
    area_time = integrate_area(
        [point.time for point in points[start_index:]],
        [point.area_a2 for point in points[start_index:]],
    )
    accounting = dissolution_event_counts(events, start.time, end.time)
    gross_si = accounting["gross_si"]
    adsorb_si = accounting["adsorb_si"]
    net_si = accounting["net_si"]
    gross_al = accounting["gross_al"]
    adsorb_al = accounting["adsorb_al"]
    net_al = accounting["net_al"]
    ratio = net_si / net_al if net_si > 0 and net_al > 0 else None
    rate = ReplicaRate(
        scenario,
        replica,
        seed,
        start.step,
        end.step,
        start.time,
        end.time,
        area_time,
        gross_si,
        adsorb_si,
        net_si,
        gross_al,
        adsorb_al,
        net_al,
        event_count_to_flux(gross_si, area_time),
        event_count_to_flux(net_si, area_time),
        event_count_to_flux(gross_al, area_time),
        event_count_to_flux(net_al, area_time),
        ratio,
        gate.passed,
    )
    return rate, gate, populations


def _percentile(values: Sequence[float], fraction: float) -> float:
    ordered = sorted(values)
    index = round((len(ordered) - 1) * fraction)
    return ordered[index]


def bootstrap_summary(
    values: Sequence[float], label: str
) -> tuple[float, float, float]:
    if not values or any(not math.isfinite(value) for value in values):
        raise ValueError("bootstrap requires a nonempty finite distribution")
    mean = statistics.fmean(values)
    if len(values) == 1:
        return mean, mean, mean
    seed = int.from_bytes(hashlib.sha256(label.encode("utf-8")).digest()[:8], "big")
    rng = random.Random(seed)
    boot = [
        statistics.fmean(values[rng.randrange(len(values))] for _ in values)
        for _ in range(BOOTSTRAP_RESAMPLES)
    ]
    return mean, _percentile(boot, 0.025), _percentile(boot, 0.975)


def _display(value: float | None, *, positive_for_log: bool = False) -> float | str:
    if value is None or not math.isfinite(value) or (positive_for_log and value <= 0.0):
        return "undefined"
    return value


def _scenario_deck_contract(path: Path, nominal: DeckContract) -> DeckContract:
    text = path.read_text(encoding="utf-8")
    parsed = tomllib.loads(text)
    return DeckContract(
        path,
        text,
        parsed,
        nominal.annotations,
        {
            reaction["name"]: float(reaction["rate"]["arrhenius"]["ea"])
            for reaction in parsed["reactions"]
        },
        {
            reaction["name"]: float(reaction["rate"]["arrhenius"]["prefactor"])
            for reaction in parsed["reactions"]
        },
        nominal.step_limit,
        nominal.report_every,
    )


def analyze_campaign(raw_root: Path, out_dir: Path) -> dict:
    if out_dir.exists():
        raise ValueError(f"refusing to overwrite analysis directory: {out_dir}")
    manifest = json.loads((raw_root / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("status") != "complete":
        raise ValueError("campaign manifest is not complete")
    seeds = validate_seeds(manifest.get("seeds", []))
    scenario_records = manifest.get("scenarios", [])
    expected_scenarios = [asdict(value) for value in scenarios()]
    if [
        {key: record.get(key) for key in ("name", "family", "perturbation")}
        for record in scenario_records
    ] != expected_scenarios:
        raise ValueError("campaign scenario coverage/order is incomplete")
    checkpoint_path = raw_root / "checkpoint.json"
    if sha256_file(checkpoint_path) != manifest.get("checkpoint_sha256"):
        raise ValueError("campaign checkpoint hash mismatch")
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    if checkpoint.get("status") != "complete":
        raise ValueError("campaign checkpoint is not complete")
    receipts = checkpoint.get("receipts")
    expected_identities = {
        (scenario["name"], replica, seed)
        for scenario in scenario_records
        for replica, seed in enumerate(seeds)
    }
    if (
        not isinstance(receipts, list)
        or {
            (receipt.get("scenario"), receipt.get("replica"), receipt.get("seed"))
            for receipt in receipts
            if isinstance(receipt, dict)
        }
        != expected_identities
    ):
        raise ValueError("checkpoint has duplicate or missing scenario/replica seeds")
    if len(receipts) != len(expected_identities):
        raise ValueError("checkpoint receipt count includes duplicate runs")
    for receipt in receipts:
        run_dir = Path(receipt["output"])
        for filename, expected_hash in receipt["sha256"].items():
            artifact = Path(receipt["log"]) if filename == "log" else run_dir / filename
            if sha256_file(artifact) != expected_hash:
                raise ValueError(f"checkpoint artifact hash mismatch: {artifact}")
    nominal_path = Path(scenario_records[0]["deck"])
    nominal = validate_deck(nominal_path, seeds)
    out_dir.mkdir(parents=True)

    rates: list[ReplicaRate] = []
    gates: list[dict[str, object]] = []
    population_rows: list[dict[str, object]] = []
    legacy_rows: list[dict[str, object]] = []
    raw_hashes: dict[str, str] = {}
    for scenario_record in scenario_records:
        scenario_name = scenario_record["name"]
        deck_path = Path(scenario_record["deck"])
        if sha256_file(deck_path) != scenario_record["deck_sha256"]:
            raise ValueError(f"scenario deck hash mismatch: {scenario_name}")
        contract = _scenario_deck_contract(deck_path, nominal)
        if scenario_record["family"] is not None:
            expected = perturb_deck(
                nominal, scenario_record["family"], scenario_record["perturbation"]
            )
            if expected != contract.text:
                raise ValueError(f"scenario isolation mismatch: {scenario_name}")
        for replica, seed in enumerate(seeds):
            run_dir = (
                raw_root / "runs" / scenario_name / f"replica-{replica:02d}-seed-{seed}"
            )
            rate, gate, populations = _analyze_replica(
                scenario_name, replica, seed, run_dir, contract
            )
            rates.append(rate)
            gates.append(
                {
                    "scenario": scenario_name,
                    "replica": replica,
                    "seed": seed,
                    **asdict(gate),
                }
            )
            for population in populations:
                si_bins = [population.states[f"Si.oh{index}"] for index in range(5)]
                al_bins = [population.states[f"Al.l{index}"] for index in range(7)]
                population_rows.append(
                    {
                        "scenario": scenario_name,
                        "replica": replica,
                        "seed": seed,
                        "step": population.step,
                        "time": population.time,
                        **population.states,
                    }
                )
                legacy_rows.append(
                    {
                        "scenario": scenario_name,
                        "replica": replica,
                        "seed": seed,
                        "step": population.step,
                        "time": population.time,
                        "Si_total": sum(si_bins),
                        **{
                            f"Si_oh{index}": value
                            for index, value in enumerate(si_bins)
                        },
                        "Al_total": sum(al_bins),
                        **{
                            f"Al_l{index}": value for index, value in enumerate(al_bins)
                        },
                    }
                )
            for filename in (
                "events.jsonl",
                "populations.csv",
                "observables.csv",
                "snapshot.pgif.json",
            ):
                relative = (run_dir / filename).relative_to(raw_root).as_posix()
                raw_hashes[relative] = sha256_file(run_dir / filename)

    rate_fields = list(ReplicaRate.__dataclass_fields__)
    rate_rows = []
    for rate in rates:
        row = asdict(rate)
        row["si_al_net_ratio"] = _display(rate.si_al_net_ratio)
        row["log10_gross_si_flux"] = _display(
            math.log10(rate.gross_si_flux) if rate.gross_si_flux > 0 else None
        )
        row["log10_net_si_flux"] = _display(
            math.log10(rate.net_si_flux) if rate.net_si_flux > 0 else None
        )
        row["log10_gross_al_flux"] = _display(
            math.log10(rate.gross_al_flux) if rate.gross_al_flux > 0 else None
        )
        row["log10_net_al_flux"] = _display(
            math.log10(rate.net_al_flux) if rate.net_al_flux > 0 else None
        )
        rate_rows.append(row)
    _write_csv_atomic(
        out_dir / "per-replica-rates.csv",
        [
            *rate_fields,
            "log10_gross_si_flux",
            "log10_net_si_flux",
            "log10_gross_al_flux",
            "log10_net_al_flux",
        ],
        rate_rows,
    )

    ensemble_rows = []
    by_scenario = {
        scenario.name: [rate for rate in rates if rate.scenario == scenario.name]
        for scenario in scenarios()
    }
    metrics = ("gross_si_flux", "net_si_flux", "gross_al_flux", "net_al_flux")
    for scenario_name, members in by_scenario.items():
        for metric in metrics:
            values = [getattr(member, metric) for member in members]
            mean, low, high = bootstrap_summary(values, f"{scenario_name}:{metric}")
            ensemble_rows.append(
                {
                    "scenario": scenario_name,
                    "metric": metric,
                    "mean": mean,
                    "ci95_low": low,
                    "ci95_high": high,
                    "log10_mean": _display(math.log10(mean) if mean > 0 else None),
                }
            )
        ratios = [member.si_al_net_ratio for member in members]
        defined_ratios = [value for value in ratios if value is not None]
        if len(defined_ratios) == len(ratios):
            mean, low, high = bootstrap_summary(
                defined_ratios, f"{scenario_name}:si_al_net_ratio"
            )
            ensemble_rows.append(
                {
                    "scenario": scenario_name,
                    "metric": "si_al_net_ratio",
                    "mean": mean,
                    "ci95_low": low,
                    "ci95_high": high,
                    "log10_mean": _display(math.log10(mean) if mean > 0 else None),
                }
            )
        else:
            ensemble_rows.append(
                {
                    "scenario": scenario_name,
                    "metric": "si_al_net_ratio",
                    "mean": "undefined",
                    "ci95_low": "undefined",
                    "ci95_high": "undefined",
                    "log10_mean": "undefined",
                }
            )
    _write_csv_atomic(
        out_dir / "ensemble-rates.csv",
        ["scenario", "metric", "mean", "ci95_low", "ci95_high", "log10_mean"],
        ensemble_rows,
    )

    state_names = [
        f"{kind['name']}.{state['name']}"
        for kind in nominal.parsed["kinds"]
        for state in kind["states"]
    ]
    _write_csv_atomic(
        out_dir / "state-populations.csv",
        ["scenario", "replica", "seed", "step", "time", *state_names],
        population_rows,
    )
    legacy_fields = [
        "scenario",
        "replica",
        "seed",
        "step",
        "time",
        "Si_total",
        "Si_oh0",
        "Si_oh1",
        "Si_oh2",
        "Si_oh3",
        "Si_oh4",
        "Al_total",
        "Al_l0",
        "Al_l1",
        "Al_l2",
        "Al_l3",
        "Al_l4",
        "Al_l5",
        "Al_l6",
    ]
    if len(legacy_fields) != 19:
        raise AssertionError("legacy series must be 15 fields plus four identities")
    _write_csv_atomic(
        out_dir / "results-dat-equivalent.csv", legacy_fields, legacy_rows
    )

    nominal_by_seed = {rate.seed: rate for rate in by_scenario["nominal"]}
    sensitivity_rows = []
    for family in FAMILY_REACTIONS:
        row: dict[str, object] = {
            "family": family,
            "declared_response": SENSITIVITY_RESPONSE,
        }
        responses = []
        undefined = False
        for perturbation in PERTURBATIONS:
            scenario_name = f"{family}__{perturbation}"
            perturbed = {rate.seed: rate for rate in by_scenario[scenario_name]}
            for species in ("si", "al"):
                deltas = []
                for seed in seeds:
                    before = getattr(nominal_by_seed[seed], f"net_{species}_flux")
                    after = getattr(perturbed[seed], f"net_{species}_flux")
                    if before <= 0.0 or after <= 0.0:
                        undefined = True
                        deltas = []
                        break
                    deltas.append(math.log10(after) - math.log10(before))
                key = f"{perturbation}_{species}_paired_mean_delta_log10"
                if deltas:
                    value = statistics.fmean(deltas)
                    row[key] = value
                    responses.append(abs(value))
                else:
                    row[key] = "undefined"
        row["max_abs_paired_mean_delta_log10"] = (
            "undefined"
            if undefined or len(responses) != len(PERTURBATIONS) * 2
            else max(responses)
        )
        sensitivity_rows.append(row)
    sensitivity_rows.sort(
        key=lambda row: (
            row["max_abs_paired_mean_delta_log10"] == "undefined",
            -float(row["max_abs_paired_mean_delta_log10"])
            if row["max_abs_paired_mean_delta_log10"] != "undefined"
            else 0.0,
            row["family"],
        )
    )
    for rank, row in enumerate(sensitivity_rows, start=1):
        row["rank"] = (
            rank
            if row["max_abs_paired_mean_delta_log10"] != "undefined"
            else "undefined"
        )
    sensitivity_fields = [
        "rank",
        "family",
        "declared_response",
        "max_abs_paired_mean_delta_log10",
        *[
            f"{perturbation}_{species}_paired_mean_delta_log10"
            for perturbation in PERTURBATIONS
            for species in ("si", "al")
        ],
    ]
    _write_csv_atomic(
        out_dir / "sensitivity-ranking.csv", sensitivity_fields, sensitivity_rows
    )
    write_json_atomic(out_dir / "steady-state-gates.json", gates)

    output_hashes = {
        path.name: sha256_file(path)
        for path in sorted(out_dir.iterdir())
        if path.is_file() and path.name != "verification.json"
    }
    verification = {
        "schema": 1,
        "passed": all(gate["passed"] for gate in gates),
        "survey_tier": True,
        "temperature_k": 298.0,
        "avogadro_mol_inverse_exact": AVOGADRO_EXACT,
        "angstrom2_to_m2": ANGSTROM2_TO_M2,
        "flux_denominator": "trapezoidal integral of emitted geometric area over physical time",
        "sensitivity_response": SENSITIVITY_RESPONSE,
        "steady_state_criterion": (__doc__ or "")
        .split("Steady-state criterion", 1)[1]
        .strip(),
        "scenario_count": len(scenario_records),
        "replicas_per_scenario": len(seeds),
        "seeds": list(seeds),
        "raw_sha256": raw_hashes,
        "output_sha256": output_hashes,
        "failed_gates": [gate for gate in gates if not gate["passed"]],
    }
    write_json_atomic(out_dir / "verification.json", verification)
    return verification


def verify_campaign(raw_root: Path, out_dir: Path) -> dict:
    verification_path = out_dir / "verification.json"
    verification = json.loads(verification_path.read_text(encoding="utf-8"))
    seeds = validate_seeds(verification.get("seeds", []))
    if verification.get("scenario_count") != len(scenarios()):
        raise ValueError("verification scenario count mismatch")
    if verification.get("replicas_per_scenario") != len(seeds):
        raise ValueError("verification replica count mismatch")
    for relative, expected in verification.get("raw_sha256", {}).items():
        if sha256_file(raw_root / relative) != expected:
            raise ValueError(f"raw artifact hash mismatch: {relative}")
    for filename, expected in verification.get("output_sha256", {}).items():
        if sha256_file(out_dir / filename) != expected:
            raise ValueError(f"derived artifact hash mismatch: {filename}")
    if not verification.get("passed") or verification.get("failed_gates"):
        raise ValueError("one or more per-replica steady-state gates failed")
    return verification


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate = subparsers.add_parser("validate", help="validate the canonical A9 deck")
    validate.add_argument("deck", type=Path)
    validate.add_argument("--seeds", type=parse_seeds, default=DEFAULT_SEEDS)

    run = subparsers.add_parser(
        "run", help="run nominal plus all sensitivity ensembles"
    )
    run.add_argument("deck", type=Path)
    run.add_argument("petra_bin", type=Path)
    run.add_argument("raw_root", type=Path)
    run.add_argument("--seeds", type=parse_seeds, default=DEFAULT_SEEDS)
    run.add_argument("--workers", type=int, default=MAX_WORKERS)
    run.add_argument("--timeout", type=int, default=MAX_TIMEOUT_SECONDS)

    analyze = subparsers.add_parser("analyze", help="derive A9 campaign products")
    analyze.add_argument("raw_root", type=Path)
    analyze.add_argument("out_dir", type=Path)

    verify = subparsers.add_parser("verify", help="verify hashes, coverage, and gates")
    verify.add_argument("raw_root", type=Path)
    verify.add_argument("out_dir", type=Path)

    args = parser.parse_args(argv)
    if args.command == "validate":
        contract = validate_deck(args.deck, args.seeds)
        print(
            f"valid: reactions={len(contract.barriers)} families={len(FAMILY_REACTIONS)} "
            f"temperature=298.0 replicas={len(args.seeds)}"
        )
    elif args.command == "run":
        run_campaign(
            args.deck,
            args.petra_bin,
            args.raw_root,
            args.seeds,
            args.workers,
            args.timeout,
        )
    elif args.command == "analyze":
        verification = analyze_campaign(args.raw_root, args.out_dir)
        print(
            f"analyzed: scenarios={verification['scenario_count']} "
            f"replicas={verification['replicas_per_scenario']} "
            f"steady_state_passed={verification['passed']}"
        )
    else:
        verification = verify_campaign(args.raw_root, args.out_dir)
        print(
            f"verified: scenarios={verification['scenario_count']} "
            f"replicas={verification['replicas_per_scenario']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
