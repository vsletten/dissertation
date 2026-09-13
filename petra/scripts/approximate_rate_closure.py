#!/usr/bin/env python3
"""Run and analyze the A9 approximate kaolinite rate-closure campaign.

This is survey-tier platform testing, not calibrated or production kinetics.

Steady-state criterion (identical for nominal and every sensitivity scenario):
for each replica, require at least eight cadence samples, completion of the
configured step limit, strictly increasing sample steps and finite physical
exposure times, and final solid-cation inventory and geometric area each above
10% of their initial values.  Over the final six cadence intervals, evaluate
Si and Al separately.  A positive species needs at least 12 gross-desorption
events, four nonzero intervals, an absolute fitted end-to-end log10-flux trend
<= 0.35 decade, and an absolute first-half versus second-half log10 mean-flux
change <= 0.35 decade.  Zero observed events produce a one-sided Poisson 95%
upper rate bound from integrated area-time, never a fabricated log rate.
Replica verdicts are exactly steady-positive, steady-zero, nonsteady, absorbed,
or incomplete.  Only steady-positive is acceptance; steady-zero is complete
evidence for a typed no-dissolution outcome.  Population inventory trend/range
is emitted as a diagnostic but is not required to be static during genuine
steady dissolution.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import itertools
import json
import math
import os
import random
import re
import statistics
import subprocess
import tempfile
import threading
import time
from collections.abc import Iterable, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path

import tomllib

AVOGADRO_EXACT = 6.02214076e23
ANGSTROM2_TO_M2 = 1.0e-20
ARRHENIUS_PREFACTOR = 1.0e13
REPLICA_COUNT = 8
MAX_WORKERS = 4
MAX_TIMEOUT_SECONDS = 1_800
BOOTSTRAP_RESAMPLES = 2_000
DEFAULT_SEEDS = (90401, 90403, 90407, 90409, 90413, 90419, 90421, 90427)
RAW_SCHEMA = "a9-raw-campaign-v2"
CHECKPOINT_SCHEMA = "a9-checkpoint-v2"
VERIFICATION_SCHEMA = "a9-verification-v2"
POISSON_ZERO_EVENT_UPPER_COUNT_95 = -math.log(0.05)
REQUIRED_OBSERVABLES = frozenset(
    {"state_counts", "event_rates", "rate_spectra", "surface_area", "exposure_age"}
)
PROVENANCE_CLASSES = frozenset({"computed", "literature", "heuristic"})
PERTURBATIONS = ("ea-minus-3", "ea-plus-3", "prefactor-x0.1", "prefactor-x10")
SENSITIVITY_RESPONSE = "combined_gross_dissolution_flux_mol_m2_s"

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
REACTION_BLOCK_RE = re.compile(
    r"(?ms)^\[\[reactions\]\]\n(?P<body>.*?)(?=^\[\[reactions\]\]\n|\Z)"
)
RATE_RE = re.compile(
    r"rate = \{ arrhenius = \{ prefactor = (?P<prefactor>[^,]+), ea = (?P<ea>[^}]+) \} \}"
)


@dataclass(frozen=True)
class ReactionRegistryEntry:
    name: str
    ea_kcal_mol: float
    family: str
    provenance_class: str
    source: str
    method: str
    observable_type: str
    rationale: str


REACTION_REGISTRY = (
    ReactionRegistryEntry(
        "R0-sio-si-hydrolysis",
        27.019,
        "siloxane-neutral",
        "computed",
        "qm/ si-neutral CALC-001",
        "B3LYP/def2-SVP/DF",
        "activation_free_energy",
        "Direct Si-neutral survey free-energy barrier; not r2SCAN-3c.",
    ),
    ReactionRegistryEntry(
        "R1-sio-si-condensation",
        24.419,
        "siloxane-neutral",
        "heuristic",
        "R0 plus legacy Si reverse gap",
        "anchor minus 2.6 kcal/mol",
        "heuristic_reverse_barrier",
        "Preserves only the legacy Si forward/reverse gap.",
    ),
    ReactionRegistryEntry(
        "R2-sioal2-si-hydrolysis",
        27.019,
        "sioal-si-neutral",
        "computed",
        "qm/ si-neutral CALC-001",
        "B3LYP/def2-SVP/DF analogue",
        "activation_free_energy_analogue",
        "Explicit Si-neutral analogue; not r2SCAN-3c.",
    ),
    ReactionRegistryEntry(
        "R3-sioal2-si-condensation",
        24.419,
        "sioal-si-neutral",
        "heuristic",
        "R2 plus legacy Si reverse gap",
        "anchor minus 2.6 kcal/mol",
        "heuristic_reverse_barrier",
        "Preserves only the legacy Si forward/reverse gap.",
    ),
    ReactionRegistryEntry(
        "R4-sioal2-al-hydrolysis",
        32.221,
        "sioal-al-neutral",
        "computed",
        "qm/ al-neutral CALC-002",
        "B3LYP/def2-SVP/DF",
        "activation_free_energy",
        "Direct Al-neutral survey free-energy barrier; not r2SCAN-3c.",
    ),
    ReactionRegistryEntry(
        "R5a-sioal2-al-condensation",
        27.421,
        "sioal-al-neutral",
        "heuristic",
        "R4 plus legacy Al reverse gap",
        "anchor minus 4.8 kcal/mol",
        "heuristic_reverse_barrier",
        "Preserves only the legacy Al forward/reverse gap.",
    ),
    ReactionRegistryEntry(
        "R5b-sioal2-al-condensation",
        27.421,
        "sioal-al-neutral",
        "heuristic",
        "R4 plus legacy Al reverse gap",
        "anchor minus 4.8 kcal/mol",
        "heuristic_reverse_barrier",
        "Preserves only the legacy Al forward/reverse gap.",
    ),
    ReactionRegistryEntry(
        "R6-second-stage-hydrolysis",
        16.969,
        "connectivity-ladder",
        "literature",
        "Liu & Ruiz Pestana 2024 Q2",
        "71 kJ/mol divided by 4.184",
        "electronic_barrier_analogue",
        "Literature Q2 connectivity analogue; pH catalysis is not encoded.",
    ),
    ReactionRegistryEntry(
        "R7-second-stage-condensation",
        12.169,
        "connectivity-ladder",
        "heuristic",
        "R6 plus legacy Al reverse gap",
        "anchor minus 4.8 kcal/mol",
        "heuristic_reverse_barrier",
        "Preserves only the legacy Al forward/reverse gap.",
    ),
    ReactionRegistryEntry(
        "R8a-sioha-hydrolysis",
        19.359,
        "connectivity-ladder",
        "literature",
        "Liu & Ruiz Pestana 2024 Q3",
        "81 kJ/mol divided by 4.184",
        "electronic_barrier_analogue",
        "Literature Q3 connectivity analogue; pH catalysis is not encoded.",
    ),
    ReactionRegistryEntry(
        "R8b-sioha-hydrolysis",
        19.359,
        "connectivity-ladder",
        "literature",
        "Liu & Ruiz Pestana 2024 Q3",
        "81 kJ/mol divided by 4.184",
        "electronic_barrier_analogue",
        "Literature Q3 connectivity analogue; pH catalysis is not encoded.",
    ),
    ReactionRegistryEntry(
        "R9-sioha-condensation",
        16.759,
        "connectivity-ladder",
        "heuristic",
        "R8 plus legacy Si reverse gap",
        "anchor minus 2.6 kcal/mol",
        "heuristic_reverse_barrier",
        "Preserves only the legacy Si forward/reverse gap.",
    ),
    ReactionRegistryEntry(
        "R10-sial-hydrolysis",
        27.019,
        "sioal-si-neutral",
        "computed",
        "qm/ si-neutral CALC-001",
        "B3LYP/def2-SVP/DF analogue",
        "activation_free_energy_analogue",
        "Explicit Si-neutral analogue; not r2SCAN-3c.",
    ),
    ReactionRegistryEntry(
        "R11-sial-condensation",
        24.419,
        "sioal-si-neutral",
        "heuristic",
        "R10 plus legacy Si reverse gap",
        "anchor minus 2.6 kcal/mol",
        "heuristic_reverse_barrier",
        "Preserves only the legacy Si forward/reverse gap.",
    ),
    ReactionRegistryEntry(
        "R12-albr-hydrolysis",
        32.221,
        "sioal-al-neutral",
        "computed",
        "qm/ al-neutral CALC-002",
        "B3LYP/def2-SVP/DF analogue",
        "activation_free_energy_analogue",
        "Explicit Al-neutral analogue; not r2SCAN-3c.",
    ),
    ReactionRegistryEntry(
        "R13-albr-condensation",
        27.421,
        "sioal-al-neutral",
        "heuristic",
        "R12 plus legacy Al reverse gap",
        "anchor minus 4.8 kcal/mol",
        "heuristic_reverse_barrier",
        "Preserves only the legacy Al forward/reverse gap.",
    ),
    ReactionRegistryEntry(
        "R14-alohal-hydrolysis",
        19.0,
        "al-o-al-analogue",
        "literature",
        "Xiao & Lasaga 1994/1996",
        "base-scale literature analogue",
        "activation_energy_analogue",
        "Scale analogue only; pH catalysis is not encoded.",
    ),
    ReactionRegistryEntry(
        "R15-alohal-condensation",
        14.2,
        "al-o-al-analogue",
        "heuristic",
        "R14 plus legacy Al reverse gap",
        "anchor minus 4.8 kcal/mol",
        "heuristic_reverse_barrier",
        "Preserves only the legacy Al forward/reverse gap.",
    ),
    ReactionRegistryEntry(
        "adsorb-al",
        14.5,
        "adsorption",
        "heuristic",
        "legacy kaolinite attachment bucket",
        "legacy sign converted to Arrhenius Ea",
        "heuristic_attachment_barrier",
        "Retains the legacy 14.5 kcal/mol Al attachment bucket.",
    ),
    ReactionRegistryEntry(
        "adsorb-si",
        6.4,
        "adsorption",
        "heuristic",
        "legacy kaolinite attachment bucket",
        "legacy sign converted to Arrhenius Ea",
        "heuristic_attachment_barrier",
        "Retains the legacy 6.4 kcal/mol Si attachment bucket.",
    ),
    ReactionRegistryEntry(
        "desorb-al",
        32.221,
        "cation-desorption",
        "heuristic",
        "CALC-002 anchor plus legacy environment ladder",
        "mixed computed anchor and heuristic modifiers",
        "mixed_desorption_proxy",
        "CALC-005 is unavailable; this is not a computed desorption barrier.",
    ),
    ReactionRegistryEntry(
        "desorb-si",
        30.053,
        "cation-desorption",
        "heuristic",
        "r2SCAN-3c Si-neutral electronic value plus legacy environment ladder",
        "mixed electronic anchor and heuristic modifiers",
        "mixed_desorption_proxy",
        "30.053 is electronic and CALC-005 is unavailable; this is not a computed desorption barrier.",
    ),
)
REGISTRY_BY_NAME = {entry.name: entry for entry in REACTION_REGISTRY}

PROVENANCE_FIELDS = (
    "record_type",
    "reaction",
    "ea_kcal_mol",
    "family",
    "provenance_class",
    "source",
    "method",
    "observable_type",
    "rationale",
    "constant_name",
    "constant_value",
    "constant_unit",
    "conversion_expression",
)
DERIVED_FILES = frozenset(
    {
        "per-replica-rates.csv",
        "ensemble-rates.csv",
        "stoichiometry.csv",
        "state-populations.csv",
        "results-dat-equivalent.csv",
        "sensitivity-ranking.csv",
        "provenance-conversions.csv",
        "steady-state-gates.json",
        "verification.json",
    }
)


def _registry_comment(entry: ReactionRegistryEntry) -> str:
    values = (
        ("family", entry.family),
        ("provenance", entry.provenance_class),
        ("source", entry.source),
        ("method", entry.method),
        ("observable", entry.observable_type),
        ("rationale", entry.rationale),
    )
    return "\n".join(f"# a9-{key} = {json.dumps(value)}" for key, value in values)


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
    n_sites: int
    reaction_names: tuple[str, ...]
    steps: tuple[int, ...]
    times: tuple[float, ...]
    names: tuple[str, ...]
    counts: dict[str, int]


@dataclass(frozen=True)
class PopulationRow:
    step: int
    time_s: float
    states: dict[str, int]


@dataclass(frozen=True)
class ObservableSample:
    step: int
    time_s: float
    values: dict[str, list[float]]


@dataclass(frozen=True)
class SteadyPoint:
    step: int
    time_s: float
    gross_si_events: int
    gross_al_events: int
    area_a2: float
    solid_si_cations: int
    solid_al_cations: int

    @property
    def solid_cations(self) -> int:
        return self.solid_si_cations + self.solid_al_cations


@dataclass(frozen=True)
class SpeciesSteadyDiagnostic:
    status: str
    gross_events: int
    positive_intervals: int
    trend_decades: float | None
    half_change_decades: float | None
    upper_95_mol_m2_s: float | None


@dataclass(frozen=True)
class PopulationSteadyDiagnostic:
    initial_solid_cations: int
    final_solid_cations: int
    final_fraction: float | None
    tail_relative_range: float | None
    tail_fractional_trend: float | None
    stability_status: str


@dataclass(frozen=True)
class SteadyStateGate:
    status: str
    acceptance_passed: bool
    evidence_complete: bool
    dissolution_outcome: str
    reasons: tuple[str, ...]
    window_start_step: int | None
    window_end_step: int | None
    window_start_time_s: float | None
    window_end_time_s: float | None
    area_time_a2_s: float | None
    si: SpeciesSteadyDiagnostic | None
    al: SpeciesSteadyDiagnostic | None
    si_population: PopulationSteadyDiagnostic
    al_population: PopulationSteadyDiagnostic
    population_initial_solid_cations: int
    population_final_solid_cations: int
    population_final_fraction: float | None
    population_tail_relative_range: float | None
    population_tail_fractional_trend: float | None
    population_stability_status: str


@dataclass(frozen=True)
class ReplicaRate:
    scenario: str
    replica: int
    seed: int
    window_start_step: int | None
    window_end_step: int | None
    window_start_time_s: float | None
    window_end_time_s: float | None
    area_time_a2_s: float | None
    gross_si_events: int
    adsorb_si_events: int
    net_si_events: int
    gross_al_events: int
    adsorb_al_events: int
    net_al_events: int
    gross_si_flux_mol_m2_s: float | None
    net_si_flux_mol_m2_s: float | None
    gross_al_flux_mol_m2_s: float | None
    net_al_flux_mol_m2_s: float | None
    gross_si_upper_95_mol_m2_s: float | None
    gross_al_upper_95_mol_m2_s: float | None
    si_al_net_ratio_dimensionless: float | None
    steady_state_status: str
    acceptance_passed: bool


def _finite_number(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{label} must be numeric")
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
    if len(normalized) != REPLICA_COUNT:
        raise ValueError(f"campaign requires exactly {REPLICA_COUNT} replicas")
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
        raise TypeError("deck must declare reactions as a list")
    names = [reaction.get("name") for reaction in reactions]
    registry_names = [entry.name for entry in REACTION_REGISTRY]
    if names != registry_names:
        raise ValueError(
            "reaction coverage/order must exactly match the authoritative 22-reaction registry"
        )

    blocks = list(REACTION_BLOCK_RE.finditer(text))
    if len(blocks) != len(REACTION_REGISTRY):
        raise ValueError("deck must contain exactly 22 reaction blocks")
    annotations: dict[str, tuple[str, str]] = {}
    for block, entry in zip(blocks, REACTION_REGISTRY, strict=True):
        body = block.group("body")
        prefix = f'name = "{entry.name}"\n{_registry_comment(entry)}\n'
        if not body.startswith(prefix):
            raise ValueError(f"{entry.name}: provenance comments do not match registry")
        annotations[entry.name] = (entry.family, entry.provenance_class)
    if set(REGISTRY_BY_NAME) != set(EXPECTED_REACTIONS):
        raise AssertionError("family partition and reaction registry disagree")
    if any(
        entry.provenance_class not in PROVENANCE_CLASSES for entry in REACTION_REGISTRY
    ):
        raise AssertionError("reaction registry has an invalid provenance class")
    for family, expected in FAMILY_REACTIONS.items():
        actual = {entry.name for entry in REACTION_REGISTRY if entry.family == family}
        if actual != set(expected):
            raise AssertionError(f"registry family {family!r} is incomplete")

    barriers: dict[str, float] = {}
    prefactors: dict[str, float] = {}
    for reaction, entry in zip(reactions, REACTION_REGISTRY, strict=True):
        name = reaction["name"]
        try:
            arrhenius = reaction["rate"]["arrhenius"]
            prefactor = _finite_number(arrhenius["prefactor"], f"{name} prefactor")
            barrier = _finite_number(arrhenius["ea"], f"{name} barrier")
        except (KeyError, TypeError) as exc:
            raise ValueError(f"{name} must have an Arrhenius rate") from exc
        if set(arrhenius) != {"prefactor", "ea"}:
            raise ValueError(
                f"{name} Arrhenius schema must contain only prefactor and ea"
            )
        if prefactor != ARRHENIUS_PREFACTOR:
            raise ValueError(f"{name} prefactor must be exactly 1e13 s^-1")
        if barrier != entry.ea_kcal_mol:
            raise ValueError(
                f"{name} Ea must exactly match registry value {entry.ea_kcal_mol} kcal/mol"
            )
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
    missing = [
        name
        for name in required
        if not (output / name).is_file() or (output / name).stat().st_size == 0
    ]
    if missing:
        raise RuntimeError(f"seed {seed} missing Petra outputs: {missing}")
    hashes = {name: sha256_file(output / name) for name in required}
    hashes["log"] = sha256_file(log_path)
    return {
        "schema": "a9-run-receipt-v2",
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
        "schema": RAW_SCHEMA,
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
        "petra_binary": str(petra_bin),
        "petra_binary_sha256": sha256_file(petra_bin),
        "scenarios": scenario_records,
    }
    write_json_atomic(raw_root / "manifest.json", manifest)
    receipts: list[dict] = []
    write_json_atomic(
        raw_root / "checkpoint.json",
        {"schema": CHECKPOINT_SCHEMA, "status": "running", "receipts": receipts},
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
                    index, _seed = futures[future]
                    receipt = future.result()
                    receipt.update({"scenario": scenario_name, "replica": index})
                    with lock:
                        receipts.append(receipt)
                        receipts.sort(
                            key=lambda item: (item["scenario"], item["replica"])
                        )
                        write_json_atomic(
                            raw_root / "checkpoint.json",
                            {
                                "schema": CHECKPOINT_SCHEMA,
                                "status": "running",
                                "receipts": receipts,
                            },
                        )
    except BaseException:
        write_json_atomic(
            raw_root / "checkpoint.json",
            {"schema": CHECKPOINT_SCHEMA, "status": "failed", "receipts": receipts},
        )
        raise
    manifest["status"] = "complete"
    manifest["completed_runs"] = len(receipts)
    manifest["checkpoint_sha256"] = sha256_file(raw_root / "checkpoint.json")
    write_json_atomic(raw_root / "manifest.json", manifest)
    write_json_atomic(
        raw_root / "checkpoint.json",
        {"schema": CHECKPOINT_SCHEMA, "status": "complete", "receipts": receipts},
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
        _exact_keys(
            header,
            {
                "petra_traj",
                "deck",
                "seed",
                "n_sites",
                "states",
                "state_types",
                "reactions",
            },
            f"{path} event header",
        )
        reaction_names = header.get("reactions")
        states = header.get("states")
        state_types = header.get("state_types")
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
        expected_state_types = [
            state["occupant"]
            for kind in contract.parsed["kinds"]
            for state in kind["states"]
        ]
        if state_types != expected_state_types:
            raise ValueError(f"{path}: state-type table does not match deck")
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

        steps: list[int] = []
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
            if type(step) is not int or step <= previous_step:
                raise ValueError(
                    f"{path}:{line_number}: event steps must be strictly increasing"
                )
            event_time = _finite_number(event_time, f"{path}:{line_number} time")
            if event_time < previous_time:
                raise ValueError(
                    f"{path}:{line_number}: event times must be finite and nondecreasing"
                )
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
            steps.append(step)
            times.append(event_time)
            names.append(name)
            counts[name] += 1
    return EventData(
        seed,
        n_sites,
        tuple(reaction_names),
        tuple(steps),
        tuple(times),
        tuple(names),
        counts,
    )


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
                sample_time_s = _finite_number(
                    float(row["time"]), f"{path}:{line_number} time"
                )
                states = {state: int(row[state]) for state in expected_states}
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    f"{path}:{line_number}: invalid population row"
                ) from exc
            if step <= previous_step or sample_time_s < previous_time:
                raise ValueError(
                    f"{path}:{line_number}: nonmonotonic population sample"
                )
            if any(value < 0 for value in states.values()):
                raise ValueError(f"{path}:{line_number}: negative population")
            result.append(PopulationRow(step, sample_time_s, states))
            previous_step, previous_time = step, sample_time_s
    if not result:
        raise ValueError(f"{path}: no population samples")
    return result


def parse_observables(
    path: Path,
    expected_state_count: int,
    expected_reaction_count: int,
    expected_seed: int | None = None,
) -> dict[int, ObservableSample]:
    columns = ["replica", "seed", "step", "time", "kind", "index", "value"]
    grouped: dict[int, dict[str, object]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != columns:
            raise ValueError(f"{path}: invalid observables columns")
        for line_number, row in enumerate(reader, start=2):
            try:
                step = int(row["step"])
                sample_time_s = _finite_number(
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
                or sample_time_s < 0
                or index < 0
                or row["kind"] not in REQUIRED_OBSERVABLES
            ):
                raise ValueError(
                    f"{path}:{line_number}: invalid direct-run observable identity"
                )
            sample = grouped.setdefault(step, {"time_s": sample_time_s, "kinds": {}})
            if sample["time_s"] != sample_time_s:
                raise ValueError(
                    f"{path}:{line_number}: inconsistent time within sample"
                )
            kinds = sample["kinds"]
            assert isinstance(kinds, dict)
            indices = kinds.setdefault(row["kind"], {})
            if index in indices:
                raise ValueError(f"{path}:{line_number}: duplicate observable index")
            indices[index] = value
    if not grouped:
        raise ValueError(f"{path}: no observable samples")
    result: dict[int, ObservableSample] = {}
    previous_step = -1
    previous_time = -math.inf
    for step in sorted(grouped):
        raw_sample = grouped[step]
        sample_time_s = float(raw_sample["time_s"])
        if step <= previous_step or sample_time_s < previous_time:
            raise ValueError(f"{path}: nonmonotonic observable cadence")
        raw_kinds = raw_sample["kinds"]
        assert isinstance(raw_kinds, dict)
        if set(raw_kinds) != REQUIRED_OBSERVABLES:
            missing = REQUIRED_OBSERVABLES - raw_kinds.keys()
            extra = raw_kinds.keys() - REQUIRED_OBSERVABLES
            raise ValueError(
                f"{path}: sample {step} observable mismatch; missing={sorted(missing)}, extra={sorted(extra)}"
            )
        values: dict[str, list[float]] = {}
        for kind, indexed in raw_kinds.items():
            if set(indexed) != set(range(len(indexed))):
                raise ValueError(f"{path}: sparse {kind} indices at step {step}")
            values[kind] = [indexed[index] for index in range(len(indexed))]
        cardinalities = {
            "state_counts": expected_state_count,
            "event_rates": expected_reaction_count,
            "surface_area": 3,
        }
        for kind, expected in cardinalities.items():
            if len(values[kind]) != expected:
                raise ValueError(
                    f"{path}: {kind} cardinality at step {step} must be {expected}"
                )
        if not values["rate_spectra"] or not values["exposure_age"]:
            raise ValueError(
                f"{path}: variable observables must be nonempty at step {step}"
            )
        if values["surface_area"][0] < 0:
            raise ValueError(f"{path}: invalid geometric surface area")
        result[step] = ObservableSample(step, sample_time_s, values)
        previous_step, previous_time = step, sample_time_s
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


def _amount_to_flux(amount: float, area_time_a2_s: float) -> float:
    amount = _finite_number(amount, "event-equivalent amount")
    area_time = _finite_number(area_time_a2_s, "area-time denominator")
    if area_time <= 0.0:
        raise ValueError("area-time denominator must be positive")
    return (amount / AVOGADRO_EXACT) / (area_time * ANGSTROM2_TO_M2)


def event_count_to_flux(events: int, area_time_a2_s: float) -> float:
    if type(events) is not int:
        raise TypeError("event count must be an integer")
    return _amount_to_flux(float(events), area_time_a2_s)


def poisson_zero_upper_flux_95(area_time_a2_s: float) -> float:
    return _amount_to_flux(POISSON_ZERO_EVENT_UPPER_COUNT_95, area_time_a2_s)


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


def _population_diagnostics(
    points: Sequence[SteadyPoint], tail: Sequence[SteadyPoint], attribute: str
) -> PopulationSteadyDiagnostic:
    initial = getattr(points[0], attribute) if points else 0
    final = getattr(points[-1], attribute) if points else 0
    final_fraction = final / initial if initial > 0 else None
    tail_values = [getattr(point, attribute) for point in tail]
    tail_range = (
        (max(tail_values) - min(tail_values)) / tail_values[0]
        if tail_values and tail_values[0] > 0
        else None
    )
    tail_trend = (
        _linear_trend([value / tail_values[0] for value in tail_values])
        if len(tail_values) >= 2 and tail_values[0] > 0
        else None
    )
    stable = (
        final_fraction is not None
        and tail_range is not None
        and tail_trend is not None
        and abs(1.0 - final_fraction) <= 0.05
        and tail_range <= 0.05
        and abs(tail_trend) <= 0.05
    )
    return PopulationSteadyDiagnostic(
        initial,
        final,
        final_fraction,
        tail_range,
        tail_trend,
        ("stable" if stable else "evolving"),
    )


def _species_diagnostic(
    species: str,
    tail: Sequence[SteadyPoint],
    area_time_a2_s: float,
    reasons: list[str],
) -> SpeciesSteadyDiagnostic:
    attribute = f"gross_{species}_events"
    interval_fluxes: list[float] = []
    interval_events: list[int] = []
    for left, right in itertools.pairwise(tail):
        count = getattr(right, attribute) - getattr(left, attribute)
        if count < 0:
            reasons.append(f"{species} gross event counter decreased")
            count = 0
        interval_events.append(count)
        denominator = integrate_area(
            [left.time_s, right.time_s], [left.area_a2, right.area_a2]
        )
        interval_fluxes.append(
            event_count_to_flux(count, denominator) if count > 0 else 0.0
        )
    gross_events = sum(interval_events)
    positive = [value for value in interval_fluxes if value > 0.0]
    if gross_events == 0:
        return SpeciesSteadyDiagnostic(
            "zero-upper-bound",
            0,
            0,
            None,
            None,
            poisson_zero_upper_flux_95(area_time_a2_s),
        )
    if gross_events < 12 or len(positive) < 4:
        reasons.append(
            f"{species} has insufficient gross dissolution evidence "
            f"({gross_events} events, {len(positive)} positive intervals)"
        )
        return SpeciesSteadyDiagnostic(
            "censored-insufficient-events",
            gross_events,
            len(positive),
            None,
            None,
            None,
        )
    logs = [math.log10(value) for value in positive]
    trend = _linear_trend(logs)
    first = statistics.fmean(positive[: len(positive) // 2])
    second = statistics.fmean(positive[len(positive) // 2 :])
    half_change = math.log10(second / first)
    status = "steady-positive"
    if abs(trend) > 0.35:
        reasons.append(f"{species} fitted dissolution trend exceeds 0.35 decade")
        status = "nonsteady"
    if abs(half_change) > 0.35:
        reasons.append(f"{species} half-window dissolution shift exceeds 0.35 decade")
        status = "nonsteady"
    return SpeciesSteadyDiagnostic(
        status, gross_events, len(positive), trend, half_change, None
    )


def assess_steady_state(
    points: Sequence[SteadyPoint], expected_steps: int
) -> SteadyStateGate:
    reasons: list[str] = []
    tail = points[-7:] if len(points) >= 7 else points
    si_population = _population_diagnostics(points, tail, "solid_si_cations")
    al_population = _population_diagnostics(points, tail, "solid_al_cations")
    population = _population_diagnostics(points, tail, "solid_cations")

    def early(status: str, reason: str) -> SteadyStateGate:
        return SteadyStateGate(
            status,
            False,
            False,
            "unresolved",
            (reason,),
            tail[0].step if tail else None,
            tail[-1].step if tail else None,
            tail[0].time_s if tail else None,
            tail[-1].time_s if tail else None,
            None,
            None,
            None,
            si_population,
            al_population,
            population.initial_solid_cations,
            population.final_solid_cations,
            population.final_fraction,
            population.tail_relative_range,
            population.tail_fractional_trend,
            population.stability_status,
        )

    if len(points) < 8:
        return early("incomplete", "insufficient cadence: need at least eight samples")
    if points[-1].step < expected_steps:
        return early("incomplete", "run stopped early before configured step limit")
    for left, right in itertools.pairwise(points):
        if (
            right.step <= left.step
            or not math.isfinite(left.time_s)
            or not math.isfinite(right.time_s)
            or right.time_s <= left.time_s
        ):
            return early(
                "incomplete",
                "sample steps and exposure times must be strictly increasing and finite",
            )
    initial = points[0]
    final = points[-1]
    if (
        initial.solid_si_cations <= 0
        or initial.solid_al_cations <= 0
        or final.solid_si_cations <= 0.10 * initial.solid_si_cations
        or final.solid_al_cations <= 0.10 * initial.solid_al_cations
        or initial.area_a2 <= 0.0
        or final.area_a2 <= 0.10 * initial.area_a2
    ):
        return early("absorbed", "solid inventory or area reached the absorbing floor")

    area_time = integrate_area(
        [point.time_s for point in tail], [point.area_a2 for point in tail]
    )
    si = _species_diagnostic("si", tail, area_time, reasons)
    al = _species_diagnostic("al", tail, area_time, reasons)
    species_statuses = {si.status, al.status}
    if (
        "nonsteady" in species_statuses
        or "censored-insufficient-events" in species_statuses
    ):
        status = "nonsteady"
        outcome = "unresolved"
    elif species_statuses == {"zero-upper-bound"}:
        status = "steady-zero"
        outcome = "no-dissolution"
    elif "zero-upper-bound" in species_statuses:
        status = "steady-zero"
        outcome = "species-zero-upper-bound"
    else:
        status = "steady-positive"
        outcome = "positive-dissolution"
    return SteadyStateGate(
        status,
        status == "steady-positive",
        status in {"steady-positive", "steady-zero"},
        outcome,
        tuple(reasons),
        tail[0].step,
        tail[-1].step,
        tail[0].time_s,
        tail[-1].time_s,
        area_time,
        si,
        al,
        si_population,
        al_population,
        population.initial_solid_cations,
        population.final_solid_cations,
        population.final_fraction,
        population.tail_relative_range,
        population.tail_fractional_trend,
        population.stability_status,
    )


def _counts_through_step(events: EventData, step: int, names: set[str]) -> int:
    end = bisect.bisect_right(events.steps, step)
    return sum(1 for name in events.names[:end] if name in names)


def _counts_between_steps(
    events: EventData, start_step: int, end_step: int, name: str
) -> int:
    left = bisect.bisect_right(events.steps, start_step)
    right = bisect.bisect_right(events.steps, end_step)
    return sum(1 for value in events.names[left:right] if value == name)


def dissolution_event_counts(
    events: EventData, start_step: int, end_step: int
) -> dict[str, int]:
    """Count Si/Al events in the exact event-step window (start, end]."""

    if (
        type(start_step) is not int
        or type(end_step) is not int
        or start_step < 0
        or end_step <= start_step
    ):
        raise ValueError("event-accounting steps must define an increasing window")
    result: dict[str, int] = {}
    for species in ("si", "al"):
        gross = _counts_between_steps(events, start_step, end_step, f"desorb-{species}")
        adsorbed = _counts_between_steps(
            events, start_step, end_step, f"adsorb-{species}"
        )
        result[f"gross_{species}"] = gross
        result[f"adsorb_{species}"] = adsorbed
        result[f"net_{species}"] = gross - adsorbed
    return result


def _solid_cations(row: PopulationRow, species: str) -> int:
    return sum(
        count
        for state, count in row.states.items()
        if state.startswith(f"{species}.") and not state.endswith(".empty")
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
    if any(sum(row.states.values()) != events.n_sites for row in populations):
        raise ValueError(
            f"{run_dir}: population cardinality disagrees with event header"
        )
    observables = parse_observables(
        run_dir / "observables.csv", len(state_names), len(REACTION_REGISTRY), seed
    )
    final_event_step = events.steps[-1] if events.steps else 0
    if populations[-1].step != final_event_step:
        raise ValueError(f"{run_dir}: final event step disagrees with population step")
    expected_sample_steps = list(
        range(0, populations[-1].step + 1, contract.report_every)
    )
    if not expected_sample_steps or expected_sample_steps[-1] != populations[-1].step:
        expected_sample_steps.append(populations[-1].step)
    population_steps = [row.step for row in populations]
    if population_steps != expected_sample_steps:
        raise ValueError(f"{run_dir}: population cadence does not match deck")
    if list(observables) != population_steps:
        raise ValueError(
            f"{run_dir}: observables/populations must align exactly by step"
        )

    points: list[SteadyPoint] = []
    for row in populations:
        sample = observables[row.step]
        state_counts = sample.values["state_counts"]
        if any(
            not value.is_integer() or int(value) != row.states[name]
            for name, value in zip(state_names, state_counts, strict=True)
        ):
            raise ValueError(
                f"{run_dir}: state-count observable disagrees with populations"
            )
        points.append(
            SteadyPoint(
                row.step,
                sample.time_s,
                _counts_through_step(events, row.step, {"desorb-si"}),
                _counts_through_step(events, row.step, {"desorb-al"}),
                sample.values["surface_area"][0],
                _solid_cations(row, "Si"),
                _solid_cations(row, "Al"),
            )
        )
    gate = assess_steady_state(points, contract.step_limit)
    start = points[-7] if len(points) >= 7 else points[0]
    end = points[-1]
    area_time: float | None = None
    accounting = {
        "gross_si": 0,
        "adsorb_si": 0,
        "net_si": 0,
        "gross_al": 0,
        "adsorb_al": 0,
        "net_al": 0,
    }
    if len(points) >= 2 and end.step > start.step:
        try:
            area_time = integrate_area(
                [point.time_s for point in points[-7:]],
                [point.area_a2 for point in points[-7:]],
            )
        except ValueError:
            area_time = None
        accounting = dissolution_event_counts(events, start.step, end.step)
    gross_si = accounting["gross_si"]
    adsorb_si = accounting["adsorb_si"]
    net_si = accounting["net_si"]
    gross_al = accounting["gross_al"]
    adsorb_al = accounting["adsorb_al"]
    net_al = accounting["net_al"]
    ratio = net_si / net_al if net_si > 0 and net_al > 0 else None

    def flux(count: int) -> float | None:
        return event_count_to_flux(count, area_time) if area_time is not None else None

    rate = ReplicaRate(
        scenario,
        replica,
        seed,
        start.step if area_time is not None else None,
        end.step if area_time is not None else None,
        start.time_s if area_time is not None else None,
        end.time_s if area_time is not None else None,
        area_time,
        gross_si,
        adsorb_si,
        net_si,
        gross_al,
        adsorb_al,
        net_al,
        flux(gross_si),
        flux(net_si),
        flux(gross_al),
        flux(net_al),
        gate.si.upper_95_mol_m2_s if gate.si is not None else None,
        gate.al.upper_95_mol_m2_s if gate.al is not None else None,
        ratio,
        gate.status,
        gate.acceptance_passed,
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


def _exact_keys(value: object, expected: set[str], label: str) -> dict:
    if not isinstance(value, dict) or set(value) != expected:
        actual = sorted(value) if isinstance(value, dict) else type(value).__name__
        raise ValueError(f"{label} schema mismatch: {actual}")
    return value


def _valid_sha256(value: object) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def _campaign_outcome(gates: Sequence[dict[str, object]]) -> str:
    statuses = {str(gate["status"]) for gate in gates}
    for status in ("incomplete", "absorbed", "nonsteady"):
        if status in statuses:
            return status
    if "steady-zero" in statuses:
        outcomes = {str(gate["dissolution_outcome"]) for gate in gates}
        return (
            "no-dissolution"
            if outcomes == {"no-dissolution"}
            else "species-or-replica-zero-upper-bound"
        )
    return "steady-positive"


def _validate_raw_campaign(
    raw_root: Path,
) -> tuple[dict, tuple[int, ...], list[dict], DeckContract, dict[str, str]]:
    raw_root = raw_root.resolve()
    manifest_path = raw_root / "manifest.json"
    checkpoint_path = raw_root / "checkpoint.json"
    manifest = _exact_keys(
        json.loads(manifest_path.read_text(encoding="utf-8")),
        {
            "schema",
            "status",
            "survey_tier",
            "temperature_k",
            "units",
            "seeds",
            "replicas",
            "workers",
            "timeout_seconds",
            "source_deck",
            "source_deck_sha256",
            "petra_binary",
            "petra_binary_sha256",
            "scenarios",
            "completed_runs",
            "checkpoint_sha256",
        },
        "campaign manifest",
    )
    if (
        manifest["schema"] != RAW_SCHEMA
        or manifest["status"] != "complete"
        or manifest["survey_tier"] is not True
        or manifest["temperature_k"] != 298.0
        or manifest["units"] != "kcal/mol"
    ):
        raise ValueError("campaign manifest constants/status mismatch")
    seeds = validate_seeds(manifest["seeds"])
    if manifest["replicas"] != REPLICA_COUNT:
        raise ValueError("campaign manifest replica count mismatch")
    if (
        type(manifest["workers"]) is not int
        or not 1 <= manifest["workers"] <= MAX_WORKERS
        or type(manifest["timeout_seconds"]) is not int
        or not 1 <= manifest["timeout_seconds"] <= MAX_TIMEOUT_SECONDS
    ):
        raise ValueError("campaign execution bounds are invalid")

    source_deck = Path(manifest["source_deck"])
    petra_binary = Path(manifest["petra_binary"])
    for path, expected_hash, label in (
        (source_deck, manifest["source_deck_sha256"], "source deck"),
        (petra_binary, manifest["petra_binary_sha256"], "Petra binary"),
    ):
        if not _valid_sha256(expected_hash) or not path.is_file():
            raise ValueError(f"{label} binding is missing")
        if sha256_file(path) != expected_hash:
            raise ValueError(f"{label} hash mismatch")
    source_contract = validate_deck(source_deck, seeds)

    scenario_records = manifest["scenarios"]
    if not isinstance(scenario_records, list) or len(scenario_records) != len(
        scenarios()
    ):
        raise ValueError("campaign must contain exactly 29 scenarios")
    expected_scenarios = [asdict(value) for value in scenarios()]
    for record, expected in zip(scenario_records, expected_scenarios, strict=True):
        _exact_keys(
            record,
            {"name", "family", "perturbation", "deck", "deck_sha256"},
            "scenario record",
        )
        if {key: record[key] for key in expected} != expected:
            raise ValueError("campaign scenario coverage/order is incomplete")
        deck_path = raw_root / "decks" / f"{record['name']}.toml"
        if Path(record["deck"]).resolve() != deck_path or not _valid_sha256(
            record["deck_sha256"]
        ):
            raise ValueError(f"scenario deck binding mismatch: {record['name']}")
        if sha256_file(deck_path) != record["deck_sha256"]:
            raise ValueError(f"scenario deck hash mismatch: {record['name']}")
        expected_text = (
            source_contract.text
            if record["family"] is None
            else perturb_deck(source_contract, record["family"], record["perturbation"])
        )
        if deck_path.read_text(encoding="utf-8") != expected_text:
            raise ValueError(f"scenario deck bytes mismatch: {record['name']}")

    if not _valid_sha256(manifest["checkpoint_sha256"]):
        raise ValueError("checkpoint hash is empty or malformed")
    if sha256_file(checkpoint_path) != manifest["checkpoint_sha256"]:
        raise ValueError("campaign checkpoint hash mismatch")
    checkpoint = _exact_keys(
        json.loads(checkpoint_path.read_text(encoding="utf-8")),
        {"schema", "status", "receipts"},
        "campaign checkpoint",
    )
    if checkpoint["schema"] != CHECKPOINT_SCHEMA or checkpoint["status"] != "complete":
        raise ValueError("campaign checkpoint constants/status mismatch")
    receipts = checkpoint["receipts"]
    expected_identities = {
        (scenario.name, replica, seed)
        for scenario in scenarios()
        for replica, seed in enumerate(seeds)
    }
    if not isinstance(receipts, list) or len(receipts) != len(expected_identities):
        raise ValueError("checkpoint receipt count mismatch")
    identities: set[tuple[object, object, object]] = set()
    required_run_files = {
        "events.jsonl",
        "populations.csv",
        "observables.csv",
        "snapshot.pgif.json",
    }
    for receipt in receipts:
        _exact_keys(
            receipt,
            {
                "schema",
                "seed",
                "elapsed_seconds",
                "command",
                "output",
                "log",
                "sha256",
                "scenario",
                "replica",
            },
            "run receipt",
        )
        if receipt["schema"] != "a9-run-receipt-v2":
            raise ValueError("run receipt schema mismatch")
        identity = (receipt["scenario"], receipt["replica"], receipt["seed"])
        identities.add(identity)
        if identity not in expected_identities:
            raise ValueError(f"unexpected run identity: {identity}")
        scenario_name, replica, seed = identity
        expected_run = (
            raw_root
            / "runs"
            / str(scenario_name)
            / f"replica-{replica:02d}-seed-{seed}"
        )
        expected_log = (
            raw_root
            / "logs"
            / str(scenario_name)
            / f"replica-{replica:02d}-seed-{seed}.log"
        )
        expected_deck = raw_root / "decks" / f"{scenario_name}.toml"
        if (
            Path(receipt["output"]).resolve() != expected_run
            or Path(receipt["log"]).resolve() != expected_log
        ):
            raise ValueError(f"run path mismatch: {identity}")
        expected_command = [
            "nice",
            "-n",
            "10",
            str(petra_binary.resolve()),
            str(expected_deck.resolve()),
            "--seed",
            str(seed),
            "--ensemble",
            "1",
            "--out",
            str(expected_run),
            "--viz",
            "--paranoid",
        ]
        if receipt["command"] != expected_command:
            raise ValueError(f"run command mismatch: {identity}")
        elapsed = _finite_number(receipt["elapsed_seconds"], "elapsed_seconds")
        if elapsed < 0:
            raise ValueError("run elapsed_seconds must be nonnegative")
        hashes = _exact_keys(
            receipt["sha256"], required_run_files | {"log"}, "run hashes"
        )
        for filename, expected_hash in hashes.items():
            if not _valid_sha256(expected_hash):
                raise ValueError(f"empty/malformed run hash: {identity} {filename}")
            artifact = expected_log if filename == "log" else expected_run / filename
            if not artifact.is_file() or artifact.stat().st_size == 0:
                raise ValueError(f"missing/empty run artifact: {artifact}")
            if sha256_file(artifact) != expected_hash:
                raise ValueError(f"run artifact hash mismatch: {artifact}")
    if identities != expected_identities:
        raise ValueError("checkpoint has duplicate or missing scenario/replica seeds")
    if manifest["completed_runs"] != len(expected_identities):
        raise ValueError("manifest completed-run count mismatch")

    expected_files = {"manifest.json", "checkpoint.json"}
    expected_files.update(f"decks/{scenario.name}.toml" for scenario in scenarios())
    for scenario in scenarios():
        for replica, seed in enumerate(seeds):
            stem = f"{scenario.name}/replica-{replica:02d}-seed-{seed}"
            expected_files.add(f"logs/{stem}.log")
            expected_files.update(
                f"runs/{stem}/{filename}" for filename in required_run_files
            )
    actual_files = {
        path.relative_to(raw_root).as_posix()
        for path in raw_root.rglob("*")
        if path.is_file()
    }
    if actual_files != expected_files:
        raise ValueError(
            f"raw file inventory mismatch: missing={sorted(expected_files - actual_files)}, "
            f"extra={sorted(actual_files - expected_files)}"
        )
    if any(path.is_symlink() for path in raw_root.rglob("*")):
        raise ValueError("raw campaign may not contain symlinks")
    raw_hashes = {
        relative: sha256_file(raw_root / relative)
        for relative in sorted(expected_files)
    }
    nominal = validate_deck(raw_root / "decks" / "nominal.toml", seeds)
    return manifest, seeds, scenario_records, nominal, raw_hashes


def _write_provenance(path: Path) -> None:
    rows: list[dict[str, object]] = []
    for entry in REACTION_REGISTRY:
        row = {field: "" for field in PROVENANCE_FIELDS}
        row.update(
            {
                "record_type": "reaction",
                "reaction": entry.name,
                "ea_kcal_mol": entry.ea_kcal_mol,
                "family": entry.family,
                "provenance_class": entry.provenance_class,
                "source": entry.source,
                "method": entry.method,
                "observable_type": entry.observable_type,
                "rationale": entry.rationale,
            }
        )
        rows.append(row)
    for name, value, unit, expression, rationale in (
        (
            "avogadro_constant",
            AVOGADRO_EXACT,
            "mol^-1",
            "mol = dissolved_cation_events / N_A",
            "One desorbed Si or Al cation is one atom-equivalent event.",
        ),
        (
            "angstrom2_to_m2",
            ANGSTROM2_TO_M2,
            "m2/A2",
            "area_m2 = geometric_area_A2 * 1e-20",
            "Petra geometric surface area uses squared deck-cell length units (A2).",
        ),
    ):
        row = {field: "" for field in PROVENANCE_FIELDS}
        row.update(
            {
                "record_type": "conversion",
                "rationale": rationale,
                "constant_name": name,
                "constant_value": value,
                "constant_unit": unit,
                "conversion_expression": expression,
            }
        )
        rows.append(row)
    _write_csv_atomic(path, PROVENANCE_FIELDS, rows)


def _aggregate_zero_upper(members: Sequence[ReplicaRate]) -> float | None:
    if any(member.area_time_a2_s is None for member in members):
        return None
    if sum(member.gross_si_events + member.gross_al_events for member in members) != 0:
        return None
    exposure = sum(float(member.area_time_a2_s) for member in members)
    return poisson_zero_upper_flux_95(exposure)


def analyze_campaign(raw_root: Path, out_dir: Path) -> dict:
    if out_dir.exists():
        raise ValueError(f"refusing to overwrite analysis directory: {out_dir}")
    _, seeds, scenario_records, nominal, raw_hashes = _validate_raw_campaign(raw_root)
    raw_root = raw_root.resolve()
    out_dir.mkdir(parents=True)

    rates: list[ReplicaRate] = []
    gates: list[dict[str, object]] = []
    population_rows: list[dict[str, object]] = []
    legacy_rows: list[dict[str, object]] = []
    for scenario_record in scenario_records:
        scenario_name = scenario_record["name"]
        deck_path = Path(scenario_record["deck"])
        contract = _scenario_deck_contract(deck_path, nominal)
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
                        "time_s": population.time_s,
                        **population.states,
                    }
                )
                legacy_rows.append(
                    {
                        "scenario": scenario_name,
                        "replica": replica,
                        "seed": seed,
                        "step": population.step,
                        "time_s": population.time_s,
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

    rate_fields = list(ReplicaRate.__dataclass_fields__)
    log_fields = [
        "log10_gross_si_flux_mol_m2_s",
        "log10_net_si_flux_mol_m2_s",
        "log10_gross_al_flux_mol_m2_s",
        "log10_net_al_flux_mol_m2_s",
    ]
    rate_rows: list[dict[str, object]] = []
    for rate in rates:
        row = {
            key: _display(value) if value is None else value
            for key, value in asdict(rate).items()
        }
        for species in ("si", "al"):
            for basis in ("gross", "net"):
                value = getattr(rate, f"{basis}_{species}_flux_mol_m2_s")
                row[f"log10_{basis}_{species}_flux_mol_m2_s"] = _display(
                    math.log10(value) if value is not None and value > 0 else None
                )
        rate_rows.append(row)
    _write_csv_atomic(
        out_dir / "per-replica-rates.csv", [*rate_fields, *log_fields], rate_rows
    )

    by_scenario = {
        scenario.name: [rate for rate in rates if rate.scenario == scenario.name]
        for scenario in scenarios()
    }
    ensemble_fields = [
        "scenario",
        "species",
        "rate_basis",
        "status",
        "mean_mol_m2_s",
        "ci95_low_mol_m2_s",
        "ci95_high_mol_m2_s",
        "poisson_zero_upper_95_mol_m2_s",
        "log10_mean_mol_m2_s",
    ]
    ensemble_rows: list[dict[str, object]] = []
    for scenario_name, members in by_scenario.items():
        for species in ("si", "al"):
            for basis in ("gross", "net"):
                values = [
                    getattr(member, f"{basis}_{species}_flux_mol_m2_s")
                    for member in members
                ]
                defined = [value for value in values if value is not None]
                if len(defined) == len(values):
                    mean, low, high = bootstrap_summary(
                        defined, f"{scenario_name}:{species}:{basis}"
                    )
                    status = (
                        "zero-upper-bound"
                        if basis == "gross" and mean == 0
                        else "estimated"
                    )
                    upper = None
                    if status == "zero-upper-bound":
                        exposure = sum(
                            float(member.area_time_a2_s) for member in members
                        )
                        upper = poisson_zero_upper_flux_95(exposure)
                else:
                    mean = low = high = upper = None
                    status = "incomplete"
                ensemble_rows.append(
                    {
                        "scenario": scenario_name,
                        "species": species,
                        "rate_basis": basis,
                        "status": status,
                        "mean_mol_m2_s": _display(mean),
                        "ci95_low_mol_m2_s": _display(low),
                        "ci95_high_mol_m2_s": _display(high),
                        "poisson_zero_upper_95_mol_m2_s": _display(upper),
                        "log10_mean_mol_m2_s": _display(
                            math.log10(mean) if mean is not None and mean > 0 else None
                        ),
                    }
                )
    _write_csv_atomic(out_dir / "ensemble-rates.csv", ensemble_fields, ensemble_rows)

    stoichiometry_fields = [
        "scenario",
        "status",
        "si_al_net_ratio_dimensionless_mean",
        "si_al_net_ratio_dimensionless_ci95_low",
        "si_al_net_ratio_dimensionless_ci95_high",
        "net_si_events",
        "net_al_events",
    ]
    stoichiometry_rows: list[dict[str, object]] = []
    for scenario_name, members in by_scenario.items():
        ratios = [member.si_al_net_ratio_dimensionless for member in members]
        defined = [value for value in ratios if value is not None]
        if len(defined) == len(ratios):
            mean, low, high = bootstrap_summary(
                defined, f"{scenario_name}:si_al_net_ratio_dimensionless"
            )
            status = "estimated"
        else:
            mean = low = high = None
            status = "undefined-nonpositive-net"
        stoichiometry_rows.append(
            {
                "scenario": scenario_name,
                "status": status,
                "si_al_net_ratio_dimensionless_mean": _display(mean),
                "si_al_net_ratio_dimensionless_ci95_low": _display(low),
                "si_al_net_ratio_dimensionless_ci95_high": _display(high),
                "net_si_events": sum(member.net_si_events for member in members),
                "net_al_events": sum(member.net_al_events for member in members),
            }
        )
    _write_csv_atomic(
        out_dir / "stoichiometry.csv", stoichiometry_fields, stoichiometry_rows
    )

    state_names = [
        f"{kind['name']}.{state['name']}"
        for kind in nominal.parsed["kinds"]
        for state in kind["states"]
    ]
    _write_csv_atomic(
        out_dir / "state-populations.csv",
        ["scenario", "replica", "seed", "step", "time_s", *state_names],
        population_rows,
    )
    legacy_fields = [
        "scenario",
        "replica",
        "seed",
        "step",
        "time_s",
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
    sensitivity_fields = [
        "rank",
        "family",
        "declared_response",
        "status",
        "max_abs_paired_mean_delta_log10",
        "nominal_zero_upper_95_mol_m2_s",
        *[
            field
            for perturbation in PERTURBATIONS
            for field in (
                f"{perturbation}_status",
                f"{perturbation}_paired_mean_delta_log10",
                f"{perturbation}_zero_upper_95_mol_m2_s",
            )
        ],
    ]
    sensitivity_rows: list[dict[str, object]] = []
    for family in FAMILY_REACTIONS:
        row: dict[str, object] = {
            "rank": "undefined",
            "family": family,
            "declared_response": SENSITIVITY_RESPONSE,
            "status": "estimated",
            "nominal_zero_upper_95_mol_m2_s": _display(
                _aggregate_zero_upper(by_scenario["nominal"])
            ),
        }
        responses: list[float] = []
        for perturbation in PERTURBATIONS:
            scenario_name = f"{family}__{perturbation}"
            perturbed = {rate.seed: rate for rate in by_scenario[scenario_name]}
            deltas: list[float] = []
            censored_zero = False
            for seed in seeds:
                before_rate = nominal_by_seed[seed]
                after_rate = perturbed[seed]
                before = (
                    before_rate.gross_si_flux_mol_m2_s
                    + before_rate.gross_al_flux_mol_m2_s
                    if before_rate.gross_si_flux_mol_m2_s is not None
                    and before_rate.gross_al_flux_mol_m2_s is not None
                    else None
                )
                after = (
                    after_rate.gross_si_flux_mol_m2_s
                    + after_rate.gross_al_flux_mol_m2_s
                    if after_rate.gross_si_flux_mol_m2_s is not None
                    and after_rate.gross_al_flux_mol_m2_s is not None
                    else None
                )
                if before == 0 or after == 0:
                    censored_zero = True
                if before is None or after is None or before <= 0 or after <= 0:
                    deltas = []
                    break
                deltas.append(math.log10(after) - math.log10(before))
            if deltas:
                value = statistics.fmean(deltas)
                row[f"{perturbation}_status"] = "estimated"
                row[f"{perturbation}_paired_mean_delta_log10"] = value
                responses.append(abs(value))
            else:
                status = "censored-zero-gross" if censored_zero else "incomplete"
                row[f"{perturbation}_status"] = status
                row[f"{perturbation}_paired_mean_delta_log10"] = "undefined"
                row["status"] = status
            row[f"{perturbation}_zero_upper_95_mol_m2_s"] = _display(
                _aggregate_zero_upper(by_scenario[scenario_name])
            )
        row["max_abs_paired_mean_delta_log10"] = (
            max(responses) if len(responses) == len(PERTURBATIONS) else "undefined"
        )
        sensitivity_rows.append(row)
    sensitivity_rows.sort(
        key=lambda row: (
            row["max_abs_paired_mean_delta_log10"] == "undefined",
            -float(row["max_abs_paired_mean_delta_log10"])
            if row["max_abs_paired_mean_delta_log10"] != "undefined"
            else 0.0,
            str(row["family"]),
        )
    )
    next_rank = 1
    for row in sensitivity_rows:
        if row["max_abs_paired_mean_delta_log10"] != "undefined":
            row["rank"] = next_rank
            next_rank += 1
    _write_csv_atomic(
        out_dir / "sensitivity-ranking.csv", sensitivity_fields, sensitivity_rows
    )
    write_json_atomic(
        out_dir / "steady-state-gates.json",
        {"schema": "a9-steady-state-gates-v2", "gates": gates},
    )
    _write_provenance(out_dir / "provenance-conversions.csv")

    status_counts = {
        status: sum(gate["status"] == status for gate in gates)
        for status in (
            "steady-positive",
            "steady-zero",
            "nonsteady",
            "absorbed",
            "incomplete",
        )
    }
    outcome = _campaign_outcome(gates)
    output_hashes = {
        path.name: sha256_file(path)
        for path in sorted(out_dir.iterdir())
        if path.is_file() and path.name != "verification.json"
    }
    verification = {
        "schema": VERIFICATION_SCHEMA,
        "acceptance_passed": outcome == "steady-positive",
        "campaign_outcome": outcome,
        "survey_tier": True,
        "temperature_k": 298.0,
        "energy_unit": "kcal/mol",
        "flux_unit": "mol_m2_s",
        "avogadro_mol_inverse_exact": AVOGADRO_EXACT,
        "angstrom2_to_m2": ANGSTROM2_TO_M2,
        "poisson_zero_event_upper_count_95": POISSON_ZERO_EVENT_UPPER_COUNT_95,
        "flux_denominator": "trapezoidal integral of emitted geometric area_A2 over observable time_s",
        "event_window": "(window_start_step,window_end_step]",
        "sensitivity_response": SENSITIVITY_RESPONSE,
        "sensitivity_complete": len(sensitivity_rows) == len(FAMILY_REACTIONS)
        and {row["family"] for row in sensitivity_rows} == set(FAMILY_REACTIONS),
        "sensitivity_families": list(FAMILY_REACTIONS),
        "scenario_count": len(scenario_records),
        "expected_scenarios": [asdict(value) for value in scenarios()],
        "replicas_per_scenario": len(seeds),
        "seeds": list(seeds),
        "status_counts": status_counts,
        "raw_file_inventory": sorted(raw_hashes),
        "derived_file_inventory": sorted(DERIVED_FILES),
        "raw_sha256": raw_hashes,
        "output_sha256": output_hashes,
        "nonpositive_gates": [
            gate for gate in gates if gate["status"] != "steady-positive"
        ],
    }
    write_json_atomic(out_dir / "verification.json", verification)
    return json.loads((out_dir / "verification.json").read_text(encoding="utf-8"))


def _read_csv_exact(
    path: Path, fields: Sequence[str], expected_rows: int | None = None
) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != list(fields):
            raise ValueError(f"{path.name} CSV schema mismatch")
        rows = list(reader)
    if expected_rows is not None and len(rows) != expected_rows:
        raise ValueError(f"{path.name} row count mismatch")
    return rows


def _validate_derived_schemas(
    out_dir: Path,
    expected_runs: int,
    nominal: DeckContract,
    seeds: Sequence[int],
) -> None:
    expected_identities = {
        (scenario.name, str(replica), str(seed))
        for scenario in scenarios()
        for replica, seed in enumerate(seeds)
    }
    rate_fields = [
        *ReplicaRate.__dataclass_fields__,
        "log10_gross_si_flux_mol_m2_s",
        "log10_net_si_flux_mol_m2_s",
        "log10_gross_al_flux_mol_m2_s",
        "log10_net_al_flux_mol_m2_s",
    ]
    rate_rows = _read_csv_exact(
        out_dir / "per-replica-rates.csv", rate_fields, expected_runs
    )
    if {
        (row["scenario"], row["replica"], row["seed"]) for row in rate_rows
    } != expected_identities:
        raise ValueError("per-replica rate identity coverage mismatch")
    ensemble_rows = _read_csv_exact(
        out_dir / "ensemble-rates.csv",
        [
            "scenario",
            "species",
            "rate_basis",
            "status",
            "mean_mol_m2_s",
            "ci95_low_mol_m2_s",
            "ci95_high_mol_m2_s",
            "poisson_zero_upper_95_mol_m2_s",
            "log10_mean_mol_m2_s",
        ],
        len(scenarios()) * 4,
    )
    if {
        (row["scenario"], row["species"], row["rate_basis"]) for row in ensemble_rows
    } != {
        (scenario.name, species, basis)
        for scenario in scenarios()
        for species in ("si", "al")
        for basis in ("gross", "net")
    }:
        raise ValueError("ensemble rate identity coverage mismatch")
    stoichiometry_rows = _read_csv_exact(
        out_dir / "stoichiometry.csv",
        [
            "scenario",
            "status",
            "si_al_net_ratio_dimensionless_mean",
            "si_al_net_ratio_dimensionless_ci95_low",
            "si_al_net_ratio_dimensionless_ci95_high",
            "net_si_events",
            "net_al_events",
        ],
        len(scenarios()),
    )
    if {row["scenario"] for row in stoichiometry_rows} != {
        scenario.name for scenario in scenarios()
    }:
        raise ValueError("stoichiometry scenario coverage mismatch")
    sensitivity = _read_csv_exact(
        out_dir / "sensitivity-ranking.csv",
        [
            "rank",
            "family",
            "declared_response",
            "status",
            "max_abs_paired_mean_delta_log10",
            "nominal_zero_upper_95_mol_m2_s",
            *[
                field
                for perturbation in PERTURBATIONS
                for field in (
                    f"{perturbation}_status",
                    f"{perturbation}_paired_mean_delta_log10",
                    f"{perturbation}_zero_upper_95_mol_m2_s",
                )
            ],
        ],
        len(FAMILY_REACTIONS),
    )
    if {row["family"] for row in sensitivity} != set(FAMILY_REACTIONS):
        raise ValueError("sensitivity family coverage mismatch")
    provenance = _read_csv_exact(
        out_dir / "provenance-conversions.csv",
        PROVENANCE_FIELDS,
        len(REACTION_REGISTRY) + 2,
    )
    if [row["reaction"] for row in provenance[:22]] != [
        entry.name for entry in REACTION_REGISTRY
    ]:
        raise ValueError("provenance reaction coverage/order mismatch")
    state_names = [
        f"{kind['name']}.{state['name']}"
        for kind in nominal.parsed["kinds"]
        for state in kind["states"]
    ]
    population_rows = _read_csv_exact(
        out_dir / "state-populations.csv",
        ["scenario", "replica", "seed", "step", "time_s", *state_names],
    )
    legacy_rows = _read_csv_exact(
        out_dir / "results-dat-equivalent.csv",
        [
            "scenario",
            "replica",
            "seed",
            "step",
            "time_s",
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
        ],
    )
    if not population_rows or len(population_rows) != len(legacy_rows):
        raise ValueError("population/legacy series cardinality mismatch")
    population_keys = [
        (row["scenario"], row["replica"], row["seed"], row["step"])
        for row in population_rows
    ]
    legacy_keys = [
        (row["scenario"], row["replica"], row["seed"], row["step"])
        for row in legacy_rows
    ]
    if (
        population_keys != legacy_keys
        or {(scenario, replica, seed) for scenario, replica, seed, _ in population_keys}
        != expected_identities
    ):
        raise ValueError("population/legacy identity coverage mismatch")
    gate_payload = _exact_keys(
        json.loads((out_dir / "steady-state-gates.json").read_text()),
        {"schema", "gates"},
        "steady-state gates",
    )
    if (
        gate_payload["schema"] != "a9-steady-state-gates-v2"
        or len(gate_payload["gates"]) != expected_runs
    ):
        raise ValueError("steady-state gate schema/cardinality mismatch")
    gate_keys = {"scenario", "replica", "seed", *SteadyStateGate.__dataclass_fields__}
    gate_identities = set()
    for gate in gate_payload["gates"]:
        _exact_keys(gate, gate_keys, "steady-state gate")
        if gate["status"] not in {
            "steady-positive",
            "steady-zero",
            "nonsteady",
            "absorbed",
            "incomplete",
        }:
            raise ValueError("unknown steady-state status")
        gate_identities.add((gate["scenario"], str(gate["replica"]), str(gate["seed"])))
    if gate_identities != expected_identities:
        raise ValueError("steady-state gate identity coverage mismatch")


def verify_campaign(raw_root: Path, out_dir: Path) -> dict:
    _, seeds, _, nominal, raw_hashes = _validate_raw_campaign(raw_root)
    actual_files = {
        path.relative_to(out_dir).as_posix()
        for path in out_dir.rglob("*")
        if path.is_file()
    }
    if actual_files != DERIVED_FILES or any(
        path.is_symlink() for path in out_dir.rglob("*")
    ):
        raise ValueError("derived file inventory mismatch")
    verification = _exact_keys(
        json.loads((out_dir / "verification.json").read_text(encoding="utf-8")),
        {
            "schema",
            "acceptance_passed",
            "campaign_outcome",
            "survey_tier",
            "temperature_k",
            "energy_unit",
            "flux_unit",
            "avogadro_mol_inverse_exact",
            "angstrom2_to_m2",
            "poisson_zero_event_upper_count_95",
            "flux_denominator",
            "event_window",
            "sensitivity_response",
            "sensitivity_complete",
            "sensitivity_families",
            "scenario_count",
            "expected_scenarios",
            "replicas_per_scenario",
            "seeds",
            "status_counts",
            "raw_file_inventory",
            "derived_file_inventory",
            "raw_sha256",
            "output_sha256",
            "nonpositive_gates",
        },
        "verification",
    )
    if (
        verification["schema"] != VERIFICATION_SCHEMA
        or verification["survey_tier"] is not True
        or verification["temperature_k"] != 298.0
        or verification["energy_unit"] != "kcal/mol"
        or verification["flux_unit"] != "mol_m2_s"
        or verification["avogadro_mol_inverse_exact"] != AVOGADRO_EXACT
        or verification["angstrom2_to_m2"] != ANGSTROM2_TO_M2
        or verification["poisson_zero_event_upper_count_95"]
        != POISSON_ZERO_EVENT_UPPER_COUNT_95
        or verification["event_window"] != "(window_start_step,window_end_step]"
        or verification["sensitivity_response"] != SENSITIVITY_RESPONSE
        or verification["sensitivity_complete"] is not True
        or verification["sensitivity_families"] != list(FAMILY_REACTIONS)
        or verification["scenario_count"] != len(scenarios())
        or verification["expected_scenarios"]
        != [asdict(value) for value in scenarios()]
        or verification["replicas_per_scenario"] != REPLICA_COUNT
        or verification["seeds"] != list(seeds)
        or verification["raw_file_inventory"] != sorted(raw_hashes)
        or verification["derived_file_inventory"] != sorted(DERIVED_FILES)
        or verification["raw_sha256"] != raw_hashes
        or set(verification["output_sha256"]) != DERIVED_FILES - {"verification.json"}
        or any(
            not _valid_sha256(value) for value in verification["output_sha256"].values()
        )
    ):
        raise ValueError("verification constants, coverage, or hashes are invalid")
    for filename, expected_hash in verification["output_sha256"].items():
        if sha256_file(out_dir / filename) != expected_hash:
            raise ValueError(f"derived artifact hash mismatch: {filename}")
    _validate_derived_schemas(out_dir, len(scenarios()) * REPLICA_COUNT, nominal, seeds)
    with tempfile.TemporaryDirectory(prefix="a9-verify-") as directory:
        regenerated = Path(directory) / "derived"
        analyze_campaign(raw_root, regenerated)
        for filename in DERIVED_FILES:
            if (out_dir / filename).read_bytes() != (
                regenerated / filename
            ).read_bytes():
                raise ValueError(f"derived artifact is not reproducible: {filename}")
    return verification


def analyze_single_run(
    deck_path: Path, run_dir: Path, out_path: Path, seed: int
) -> dict[str, object]:
    if out_path.exists():
        raise ValueError(f"refusing to overwrite single-run analysis: {out_path}")
    contract = validate_deck(deck_path)
    rate, gate, _ = _analyze_replica("nominal", 0, seed, run_dir, contract)
    payload = {
        "schema": "a9-single-run-analysis-v2",
        "survey_tier": True,
        "acceptance_passed": gate.acceptance_passed,
        "outcome": gate.dissolution_outcome,
        "steady_state_status": gate.status,
        "rate": asdict(rate),
        "gate": asdict(gate),
    }
    write_json_atomic(out_path, payload)
    return payload


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

    analyze_run = subparsers.add_parser(
        "analyze-run",
        help="analyze one bounded direct Petra run without campaign acceptance",
    )
    analyze_run.add_argument("deck", type=Path)
    analyze_run.add_argument("run_dir", type=Path)
    analyze_run.add_argument("out_path", type=Path)
    analyze_run.add_argument("--seed", type=int, required=True)

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
            f"outcome={verification['campaign_outcome']} "
            f"acceptance_passed={verification['acceptance_passed']}"
        )
    elif args.command == "analyze-run":
        payload = analyze_single_run(args.deck, args.run_dir, args.out_path, args.seed)
        print(
            f"analyzed-run: status={payload['steady_state_status']} "
            f"outcome={payload['outcome']} "
            f"acceptance_passed={payload['acceptance_passed']}"
        )
    else:
        verification = verify_campaign(args.raw_root, args.out_dir)
        print(
            f"verified: scenarios={verification['scenario_count']} "
            f"replicas={verification['replicas_per_scenario']} "
            f"outcome={verification['campaign_outcome']} "
            f"acceptance_passed={verification['acceptance_passed']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
