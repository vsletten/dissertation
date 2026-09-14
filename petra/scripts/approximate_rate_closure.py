#!/usr/bin/env python3
"""Run and analyze the A9 approximate kaolinite rate-closure campaign.

This is survey-tier platform testing, not calibrated or production kinetics.

Each replica must complete 200,000 events and provide at least 21 equal-cadence
samples. The final 17 samples form four equal four-interval blocks. Lattice
cation lineage is replayed by site from JSONL transitions, so observed physical
dissolution rates count only original-lattice Si/Al releases, never gross or net
desorption. Sampled desorption propensity is an expected lattice-origin response
only when the entire trajectory has zero adsorption events; otherwise it is
origin-contaminated and undefined.

Si and Al response stationarity independently require all-zero blocks (typed
stationary-zero with one-sided Poisson 95% bounds) or all-positive blocks with
finite fitted log trend and half-window shift both <= 0.35 decade. Mixed
zero/positive blocks are unresolved. Tail solid Si and Al populations must each
have <= 5% relative range and <= 5% absolute fractional trend; final inventory
need not equal initial inventory. Passing stationary-zero and stationary-positive
outcomes are both scientifically complete and accepted at this survey tier.

Sensitivity ranks use the origin-safe expected propensity response. Any family
with a censored perturbation has undefined rank; ordinal ranks apply only to
fully estimated families.
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
VERIFICATION_SCHEMA = "a9-verification-v3"
POISSON_ZERO_EVENT_UPPER_COUNT_95 = -math.log(0.05)
REQUIRED_OBSERVABLES = frozenset(
    {"state_counts", "event_rates", "rate_spectra", "surface_area", "exposure_age"}
)
PROVENANCE_CLASSES = frozenset({"computed", "literature", "heuristic"})
PERTURBATIONS = ("ea-minus-3", "ea-plus-3", "prefactor-x0.1", "prefactor-x10")
SENSITIVITY_RESPONSE = (
    "combined_si_plus_al_expected_lattice_origin_flux_from_propensity_mol_m2_s"
)
PROPENSITY_RATE_BASIS = "expected_lattice_origin_from_propensity"
OBSERVED_RATE_BASIS = "lattice_origin_observed_release"
PROPENSITY_ESTIMATOR_BASIS = "integrated_ctmc_hazard"
PROPENSITY_INTERPRETATION = (
    "Trapezoidally integrated instantaneous total CTMC desorption propensity, "
    "not observed event release; used as an expected lattice-origin response only "
    "for trajectories with zero "
    "adsorption events; adsorption makes the response origin-contaminated and unresolved."
)
PROPENSITY_INTEGRAL = (
    "trapezoidal integral of instantaneous total CTMC desorb-si/desorb-al "
    "propensity over the selected observable tail"
)
OBSERVED_EVENT_ACCEPTANCE_BASIS = (
    "site-replayed original-lattice release counts plus origin-safe sampled-propensity "
    "and tail-population stationarity gates"
)
SENSITIVITY_STATISTIC = (
    "same-seed paired linear deltas; reported delta_log10 is log10(arithmetic mean "
    "perturbed response) minus log10(arithmetic mean nominal response)"
)

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
        "qm/ si-neutral CALC-002",
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
        "qm/ si-neutral CALC-002",
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
        "qm/ al-neutral CALC-003",
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
        "qm/ si-neutral CALC-002",
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
        "qm/ al-neutral CALC-003",
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
        "CALC-003 anchor plus legacy environment ladder",
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
    lattice_origin_releases: tuple[str | None, ...]
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
class PropensityExpectedGross:
    area_time_a2_s: float
    expected_gross_si_events_from_propensity: float
    expected_gross_al_events_from_propensity: float
    expected_gross_si_flux_from_propensity_mol_m2_s: float
    expected_gross_al_flux_from_propensity_mol_m2_s: float


@dataclass(frozen=True)
class SteadyPoint:
    step: int
    time_s: float
    lattice_si_releases: int
    lattice_al_releases: int
    area_a2: float
    solid_si_cations: int
    solid_al_cations: int
    si_desorb_propensity: float | None
    al_desorb_propensity: float | None

    @property
    def solid_cations(self) -> int:
        return self.solid_si_cations + self.solid_al_cations


@dataclass(frozen=True)
class SpeciesSteadyDiagnostic:
    status: str
    lattice_release_events: int
    positive_blocks: int
    observed_trend_decades: float | None
    observed_half_change_decades: float | None
    propensity_trend_decades: float | None
    propensity_half_change_decades: float | None
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
    lattice_si_release_events: int
    gross_si_events: int
    adsorb_si_events: int
    net_si_events: int
    lattice_al_release_events: int
    gross_al_events: int
    adsorb_al_events: int
    net_al_events: int
    lattice_si_release_flux_mol_m2_s: float | None
    lattice_al_release_flux_mol_m2_s: float | None
    expected_lattice_origin_si_flux_from_propensity_mol_m2_s: float | None
    expected_lattice_origin_al_flux_from_propensity_mol_m2_s: float | None
    lattice_si_upper_95_mol_m2_s: float | None
    lattice_al_upper_95_mol_m2_s: float | None
    si_al_lattice_release_ratio_dimensionless: float | None
    propensity_origin_status: str
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
    thermo = parsed.get("thermo", {})
    temperature = thermo.get("temperature")
    if type(temperature) is not float or temperature != 298.0:
        raise ValueError("thermo.temperature must be exactly the TOML float 298.0")
    if thermo.get("activity") != {"Al": 1.0e-30, "Si": 1.0e-30}:
        raise ValueError(
            "thermo.activity must fix the numerical open-flow sink at exactly 1.0e-30"
        )
    if thermo.get("mu") != {"Al": -1.0, "Si": -1.0}:
        raise ValueError(
            "thermo.mu must fix the far-from-equilibrium Al/Si reservoir at exactly -1.0 kcal/mol"
        )
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
    simulation_report_every = simulation.get("report_every")
    report_every = parsed.get("observables", {}).get("report_every")
    if type(steps) is not int or steps != 200_000:
        raise ValueError("simulation.steps must be exactly 200000")
    if type(report_every) is not int or report_every != 10_000:
        raise ValueError("observables.report_every must be exactly 10000")
    if type(simulation_report_every) is not int or simulation_report_every != 10_000:
        raise ValueError("simulation.report_every must be exactly 10000")
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
        lattice_origin_releases: list[str | None] = []
        counts = {name: 0 for name in reaction_names}
        previous_time = -math.inf
        previous_step = 0
        last_seen_state: dict[int, int] = {}
        cation_lineage: dict[int, str] = {}
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
            centers: list[tuple[int, int, int]] = []
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
            for site, old, new in normalized_changes:
                old_type = state_types[old]
                new_type = state_types[new]
                if old_type in {"Si", "Al"}:
                    cation_lineage.setdefault(site, "lattice")
                elif old_type == "vacant" and new_type in {"Si", "Al"}:
                    cation_lineage.setdefault(site, "empty")
            lattice_release: str | None = None
            if name.startswith("adsorb-"):
                cation_lineage[centers[0][0]] = "reservoir"
            elif name.startswith("desorb-"):
                site = centers[0][0]
                if cation_lineage.get(site) == "lattice":
                    lattice_release = name.removeprefix("desorb-")
                cation_lineage[site] = "empty"
            for site, _, new in normalized_changes:
                last_seen_state[site] = new
            previous_step = step
            previous_time = event_time
            steps.append(step)
            times.append(event_time)
            names.append(name)
            lattice_origin_releases.append(lattice_release)
            counts[name] += 1
    return EventData(
        seed,
        n_sites,
        tuple(reaction_names),
        tuple(steps),
        tuple(times),
        tuple(names),
        tuple(lattice_origin_releases),
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


def integrate_expected_gross_dissolution_from_propensity(
    samples: Sequence[ObservableSample],
    reaction_names: Sequence[str],
    expected_steps: Sequence[int],
) -> PropensityExpectedGross:
    """Integrate sampled Si/Al desorption hazards over one selected tail.

    ``event_rates`` is Petra's instantaneous total CTMC propensity per reaction,
    ordered exactly as the deck reaction table.  Its physical-time integral is an
    expected event count, not evidence that any event fired.
    """

    canonical_reactions = tuple(entry.name for entry in REACTION_REGISTRY)
    if tuple(reaction_names) != canonical_reactions:
        raise ValueError(
            "propensity reaction index mapping must exactly match the deck"
        )
    if len(samples) < 2:
        raise ValueError("propensity integration requires at least two samples")
    sample_steps = [sample.step for sample in samples]
    if sample_steps != list(expected_steps):
        raise ValueError("propensity sample step alignment mismatch")
    if any(type(step) is not int or step < 0 for step in sample_steps):
        raise ValueError("propensity sample steps must be nonnegative integers")

    reaction_indices = {
        species: canonical_reactions.index(f"desorb-{species}")
        for species in ("si", "al")
    }
    if len(set(reaction_indices.values())) != 2:
        raise AssertionError("Si and Al desorption reaction indices must be distinct")

    times: list[float] = []
    areas: list[float] = []
    species_propensities: dict[str, list[float]] = {"si": [], "al": []}
    previous_step = -1
    previous_time = -math.inf
    for sample in samples:
        try:
            sample_time = _finite_number(sample.time_s, "propensity sample time")
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "propensity samples require finite, strictly increasing times"
            ) from exc
        if sample.step <= previous_step or sample_time <= previous_time:
            raise ValueError(
                "propensity samples require aligned, strictly increasing steps and finite times"
            )
        raw_rates = sample.values.get("event_rates")
        if not isinstance(raw_rates, list) or len(raw_rates) != len(
            canonical_reactions
        ):
            raise ValueError(
                "propensity vector cardinality does not match reaction mapping"
            )
        rates: list[float] = []
        for raw_rate in raw_rates:
            try:
                rate = _finite_number(raw_rate, "CTMC propensity")
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    "CTMC propensities must be finite and nonnegative"
                ) from exc
            if rate < 0.0:
                raise ValueError("CTMC propensities must be finite and nonnegative")
            rates.append(rate)
        raw_area = sample.values.get("surface_area")
        if not isinstance(raw_area, list) or not raw_area:
            raise ValueError("propensity samples require geometric surface area")
        try:
            area = _finite_number(raw_area[0], "propensity geometric area")
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "propensity geometric areas must be finite and nonnegative"
            ) from exc
        if area < 0.0:
            raise ValueError(
                "propensity geometric areas must be finite and nonnegative"
            )
        times.append(sample_time)
        areas.append(area)
        for species, reaction_index in reaction_indices.items():
            species_propensities[species].append(rates[reaction_index])
        previous_step = sample.step
        previous_time = sample_time

    area_time = integrate_area(times, areas)
    expected_events: dict[str, float] = {}
    for species, propensities in species_propensities.items():
        integral = sum(
            0.5 * (left + right) * (right_time - left_time)
            for left_time, right_time, left, right in zip(
                times[:-1], times[1:], propensities[:-1], propensities[1:], strict=True
            )
        )
        if not math.isfinite(integral) or integral < 0.0:
            raise ValueError(
                "integrated CTMC propensity must be finite and nonnegative"
            )
        expected_events[species] = integral
    return PropensityExpectedGross(
        area_time,
        expected_events["si"],
        expected_events["al"],
        _amount_to_flux(expected_events["si"], area_time),
        _amount_to_flux(expected_events["al"], area_time),
    )


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
        tail_range is not None
        and tail_trend is not None
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
    attribute = f"lattice_{species}_releases"
    interval_events: list[int] = []
    for left, right in itertools.pairwise(tail):
        count = getattr(right, attribute) - getattr(left, attribute)
        if count < 0:
            reasons.append(f"{species} lattice-release counter decreased")
            count = 0
        interval_events.append(count)
    lattice_events = sum(interval_events)
    block_events = [
        sum(interval_events[index : index + 4]) for index in range(0, 16, 4)
    ]
    block_area_times = []
    for index in range(0, 16, 4):
        block = tail[index : index + 5]
        block_area_times.append(
            integrate_area(
                [point.time_s for point in block], [point.area_a2 for point in block]
            )
        )
    block_fluxes = [
        event_count_to_flux(count, denominator) if count > 0 else 0.0
        for count, denominator in zip(block_events, block_area_times, strict=True)
    ]
    propensity_attribute = f"{species}_desorb_propensity"
    propensity_values = [getattr(point, propensity_attribute) for point in tail]
    if any(value is None for value in propensity_values):
        reasons.append(f"{species} propensity response is origin-contaminated")
        return SpeciesSteadyDiagnostic(
            "origin-contaminated", lattice_events, 0, None, None, None, None, None
        )
    expected_block_fluxes = []
    for index, area_time in zip(range(0, 16, 4), block_area_times, strict=True):
        block = tail[index : index + 5]
        values = [float(getattr(point, propensity_attribute)) for point in block]
        expected_events = sum(
            0.5 * (left + right) * (right_point.time_s - left_point.time_s)
            for left_point, right_point, left, right in zip(
                block[:-1], block[1:], values[:-1], values[1:], strict=True
            )
        )
        expected_block_fluxes.append(_amount_to_flux(expected_events, area_time))
    expected_positive = [value for value in expected_block_fluxes if value > 0.0]
    if expected_positive and len(expected_positive) != 4:
        reasons.append(f"{species} propensity has mixed zero and positive blocks")
        return SpeciesSteadyDiagnostic(
            "unresolved-mixed-response", lattice_events, 0, None, None, None, None, None
        )
    expected_trend = None
    expected_half_change = None
    if expected_positive:
        expected_logs = [math.log10(value) for value in expected_block_fluxes]
        expected_trend = _linear_trend(expected_logs)
        expected_first = statistics.fmean(expected_block_fluxes[:2])
        expected_second = statistics.fmean(expected_block_fluxes[2:])
        expected_half_change = math.log10(expected_second / expected_first)
        if (
            not math.isfinite(expected_trend)
            or abs(expected_trend) > 0.35
            or not math.isfinite(expected_half_change)
            or abs(expected_half_change) > 0.35
        ):
            reasons.append(
                f"{species} sampled desorb propensity trend/shift exceeds 0.35 decade"
            )
            return SpeciesSteadyDiagnostic(
                "nonsteady",
                lattice_events,
                4,
                None,
                None,
                expected_trend,
                expected_half_change,
                None,
            )
    elif lattice_events > 0:
        reasons.append(f"{species} positive releases conflict with zero propensity")
        return SpeciesSteadyDiagnostic(
            "unresolved-mixed-response",
            lattice_events,
            0,
            None,
            None,
            expected_trend,
            expected_half_change,
            None,
        )
    positive = [value for value in block_fluxes if value > 0.0]
    if not positive:
        return SpeciesSteadyDiagnostic(
            "stationary-zero",
            0,
            0,
            None,
            None,
            expected_trend,
            expected_half_change,
            poisson_zero_upper_flux_95(area_time_a2_s),
        )
    if len(positive) != 4:
        reasons.append(f"{species} has mixed zero and positive four-block response")
        return SpeciesSteadyDiagnostic(
            "unresolved-mixed-response",
            lattice_events,
            len(positive),
            None,
            None,
            expected_trend,
            expected_half_change,
            None,
        )
    logs = [math.log10(value) for value in block_fluxes]
    trend = _linear_trend(logs)
    first = statistics.fmean(block_fluxes[:2])
    second = statistics.fmean(block_fluxes[2:])
    half_change = math.log10(second / first)
    status = "stationary-positive"
    if not math.isfinite(trend) or abs(trend) > 0.35:
        reasons.append(f"{species} fitted dissolution trend exceeds 0.35 decade")
        status = "nonsteady"
    if not math.isfinite(half_change) or abs(half_change) > 0.35:
        reasons.append(f"{species} half-window dissolution shift exceeds 0.35 decade")
        status = "nonsteady"
    return SpeciesSteadyDiagnostic(
        status,
        lattice_events,
        len(positive),
        trend,
        half_change,
        expected_trend,
        expected_half_change,
        None,
    )


def assess_steady_state(
    points: Sequence[SteadyPoint], expected_steps: int
) -> SteadyStateGate:
    reasons: list[str] = []
    tail = points[-17:] if len(points) >= 17 else points
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

    if len(points) < 21:
        return early("incomplete", "insufficient cadence: need at least 21 samples")
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
    if len({right.step - left.step for left, right in itertools.pairwise(tail)}) != 1:
        return early("incomplete", "final response blocks require equal step cadence")
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
    populations_stable = (
        si_population.stability_status == "stable"
        and al_population.stability_status == "stable"
    )
    if not populations_stable:
        reasons.append("tail solid Si/Al populations are not stationary")
    if species_statuses & {
        "nonsteady",
        "unresolved-mixed-response",
        "origin-contaminated",
    }:
        status = "nonsteady"
        outcome = "unresolved"
    elif len(species_statuses) != 1:
        reasons.append("Si/Al responses mix stationary zero and positive outcomes")
        status = "nonsteady"
        outcome = "unresolved"
    elif species_statuses == {"stationary-zero"}:
        status = "steady-zero"
        outcome = "no-dissolution"
    else:
        status = "steady-positive"
        outcome = "positive-dissolution"
    if not populations_stable:
        status = "nonsteady"
        outcome = "unresolved"
    accepted = status in {"steady-positive", "steady-zero"}
    return SteadyStateGate(
        status,
        accepted,
        accepted,
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


def _lattice_releases_through_step(events: EventData, step: int, species: str) -> int:
    end = bisect.bisect_right(events.steps, step)
    return sum(1 for value in events.lattice_origin_releases[:end] if value == species)


def _counts_between_steps(
    events: EventData, start_step: int, end_step: int, name: str
) -> int:
    left = bisect.bisect_right(events.steps, start_step)
    right = bisect.bisect_right(events.steps, end_step)
    return sum(1 for value in events.names[left:right] if value == name)


def propensity_origin_safe(events: EventData) -> bool:
    """Whether sampled desorption propensity can only represent lattice cations."""

    return events.counts["adsorb-si"] == 0 and events.counts["adsorb-al"] == 0


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
        left = bisect.bisect_right(events.steps, start_step)
        right = bisect.bisect_right(events.steps, end_step)
        lattice = sum(
            1
            for value in events.lattice_origin_releases[left:right]
            if value == species
        )
        gross = _counts_between_steps(events, start_step, end_step, f"desorb-{species}")
        adsorbed = _counts_between_steps(
            events, start_step, end_step, f"adsorb-{species}"
        )
        result[f"lattice_{species}"] = lattice
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

    origin_safe = propensity_origin_safe(events)
    si_desorb_index = events.reaction_names.index("desorb-si")
    al_desorb_index = events.reaction_names.index("desorb-al")
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
                _lattice_releases_through_step(events, row.step, "si"),
                _lattice_releases_through_step(events, row.step, "al"),
                sample.values["surface_area"][0],
                _solid_cations(row, "Si"),
                _solid_cations(row, "Al"),
                sample.values["event_rates"][si_desorb_index] if origin_safe else None,
                sample.values["event_rates"][al_desorb_index] if origin_safe else None,
            )
        )
    gate = assess_steady_state(points, contract.step_limit)
    start = points[-17] if len(points) >= 17 else points[0]
    end = points[-1]
    area_time: float | None = None
    propensity_estimate: PropensityExpectedGross | None = None
    accounting = {
        "lattice_si": 0,
        "gross_si": 0,
        "adsorb_si": 0,
        "net_si": 0,
        "lattice_al": 0,
        "gross_al": 0,
        "adsorb_al": 0,
        "net_al": 0,
    }
    if len(points) >= 2 and end.step > start.step:
        selected_points = points[-17:]
        selected_samples = [observables[point.step] for point in selected_points]
        propensity_estimate = integrate_expected_gross_dissolution_from_propensity(
            selected_samples,
            events.reaction_names,
            [point.step for point in selected_points],
        )
        area_time = propensity_estimate.area_time_a2_s
        if gate.area_time_a2_s is not None and area_time != gate.area_time_a2_s:
            raise ValueError(f"{run_dir}: propensity and observed area-time disagree")
        accounting = dissolution_event_counts(events, start.step, end.step)
    lattice_si = accounting["lattice_si"]
    gross_si = accounting["gross_si"]
    adsorb_si = accounting["adsorb_si"]
    net_si = accounting["net_si"]
    lattice_al = accounting["lattice_al"]
    gross_al = accounting["gross_al"]
    adsorb_al = accounting["adsorb_al"]
    net_al = accounting["net_al"]
    ratio = lattice_si / lattice_al if lattice_si > 0 and lattice_al > 0 else None

    def flux(count: int) -> float | None:
        return event_count_to_flux(count, area_time) if area_time is not None else None

    rate = ReplicaRate(
        scenario=scenario,
        replica=replica,
        seed=seed,
        window_start_step=start.step if area_time is not None else None,
        window_end_step=end.step if area_time is not None else None,
        window_start_time_s=start.time_s if area_time is not None else None,
        window_end_time_s=end.time_s if area_time is not None else None,
        area_time_a2_s=area_time,
        lattice_si_release_events=lattice_si,
        gross_si_events=gross_si,
        adsorb_si_events=adsorb_si,
        net_si_events=net_si,
        lattice_al_release_events=lattice_al,
        gross_al_events=gross_al,
        adsorb_al_events=adsorb_al,
        net_al_events=net_al,
        lattice_si_release_flux_mol_m2_s=flux(lattice_si),
        lattice_al_release_flux_mol_m2_s=flux(lattice_al),
        expected_lattice_origin_si_flux_from_propensity_mol_m2_s=(
            propensity_estimate.expected_gross_si_flux_from_propensity_mol_m2_s
            if propensity_estimate is not None and origin_safe
            else None
        ),
        expected_lattice_origin_al_flux_from_propensity_mol_m2_s=(
            propensity_estimate.expected_gross_al_flux_from_propensity_mol_m2_s
            if propensity_estimate is not None and origin_safe
            else None
        ),
        lattice_si_upper_95_mol_m2_s=(
            gate.si.upper_95_mol_m2_s if gate.si is not None else None
        ),
        lattice_al_upper_95_mol_m2_s=(
            gate.al.upper_95_mol_m2_s if gate.al is not None else None
        ),
        si_al_lattice_release_ratio_dimensionless=ratio,
        propensity_origin_status=(
            "origin-safe" if origin_safe else "origin-contaminated"
        ),
        steady_state_status=gate.status,
        acceptance_passed=gate.acceptance_passed,
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
        row: dict[str, object] = {field: "" for field in PROVENANCE_FIELDS}
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
    for name, value, unit, source, expression, rationale in (
        (
            "temperature",
            298.0,
            "K",
            "A9 program-card ambient-temperature execution contract",
            "T = 298.0 K exactly",
            "Ambient-temperature closure condition; validate, run, analyze, and verify reject any other value.",
        ),
        (
            "dissolved_cation_activity",
            1.0e-30,
            "dimensionless",
            "A9 far-from-equilibrium open-flow boundary condition",
            "activity(Al) = activity(Si) = 1e-30",
            "Explicit numerical open-flow sink for dissolved products in a continuously refreshed, far-from-equilibrium reservoir; not a measured pH-dependent activity.",
        ),
        (
            "dissolved_cation_mu",
            -1.0,
            "kcal/mol",
            "legacy low-chemical-potential endpoint",
            "mu(Al) = mu(Si) = -1 kcal/mol",
            "Preserves the legacy low-product-potential endpoint while the explicit activity supplies the dilute open-flow sink; Petra does not encode H+/OH- catalysis in this deck.",
        ),
        (
            "effective_consumed_cation_factor_298k",
            1.8476765674964432e-31,
            "dimensionless",
            "derived from the declared A9 reservoir",
            "activity * exp(mu / (R * T)) with R = 0.00198720425864083 kcal mol^-1 K^-1",
            "Effective mass-action factor applied to cation-consuming adsorption rules at 298 K; records the exact proxy rather than claiming a calibrated pH 3-5 chemical potential.",
        ),
        (
            "simulation_steps",
            200_000,
            "events",
            "A9 scientific-closure sampling contract",
            "simulation.steps = 200000",
            "Fixed survey-tier trajectory length providing 21 cadence samples; not a production convergence claim.",
        ),
        (
            "report_every_steps",
            10_000,
            "events",
            "A9 scientific-closure sampling contract",
            "observables.report_every = 10000",
            "Fixed cadence used for the final four equal response blocks and population stationarity gates.",
        ),
    ):
        row = {field: "" for field in PROVENANCE_FIELDS}
        row.update(
            {
                "record_type": "condition",
                "source": source,
                "observable_type": "thermodynamic_boundary_condition",
                "rationale": rationale,
                "constant_name": name,
                "constant_value": value,
                "constant_unit": unit,
                "conversion_expression": expression,
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
    for record_type, name, value, unit, expression, rationale in (
        (
            "estimator",
            "expected_lattice_origin_flux_from_propensity",
            "trapezoidal",
            "mol m^-2 s^-1",
            "expected_flux_species = integral(lambda_desorb_species(t) dt) / N_A / (integral(A_geometric(t) dt) * 1e-20)",
            PROPENSITY_INTERPRETATION,
        ),
        (
            "sensitivity-response",
            "sensitivity_response",
            SENSITIVITY_RESPONSE,
            "mol m^-2 s^-1",
            "combined response = expected_lattice_origin_si_flux_from_propensity + expected_lattice_origin_al_flux_from_propensity; delta_log10 = log10(mean paired perturbed response) - log10(mean paired nominal response)",
            "Same-seed pairs are used for every perturbation; adsorption-contaminated, incomplete, or nonpositive responses are censored and no Poisson bound enters the ranking.",
        ),
    ):
        row = {field: "" for field in PROVENANCE_FIELDS}
        row.update(
            {
                "record_type": record_type,
                "observable_type": PROPENSITY_ESTIMATOR_BASIS,
                "rationale": rationale,
                "constant_name": name,
                "constant_value": value,
                "constant_unit": unit,
                "conversion_expression": expression,
            }
        )
        rows.append(row)
    _write_csv_atomic(path, PROVENANCE_FIELDS, rows)


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
        "log10_lattice_si_release_flux_mol_m2_s",
        "log10_lattice_al_release_flux_mol_m2_s",
        "log10_expected_lattice_origin_si_flux_from_propensity_mol_m2_s",
        "log10_expected_lattice_origin_al_flux_from_propensity_mol_m2_s",
    ]
    rate_rows: list[dict[str, object]] = []
    for rate in rates:
        row = {
            key: _display(value) if value is None else value
            for key, value in asdict(rate).items()
        }
        for species in ("si", "al"):
            value = getattr(rate, f"lattice_{species}_release_flux_mol_m2_s")
            row[f"log10_lattice_{species}_release_flux_mol_m2_s"] = _display(
                math.log10(value) if value is not None and value > 0 else None
            )
            propensity_value = getattr(
                rate,
                f"expected_lattice_origin_{species}_flux_from_propensity_mol_m2_s",
            )
            row[
                f"log10_expected_lattice_origin_{species}_flux_from_propensity_mol_m2_s"
            ] = _display(
                math.log10(propensity_value)
                if propensity_value is not None and propensity_value > 0
                else None
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
            values = [
                getattr(member, f"lattice_{species}_release_flux_mol_m2_s")
                for member in members
            ]
            defined = [value for value in values if value is not None]
            if len(defined) == len(values):
                mean, low, high = bootstrap_summary(
                    defined, f"{scenario_name}:{species}:{OBSERVED_RATE_BASIS}"
                )
                status = "zero-upper-bound" if mean == 0 else "estimated"
                upper = None
                if status == "zero-upper-bound":
                    exposure = sum(float(member.area_time_a2_s) for member in members)
                    upper = poisson_zero_upper_flux_95(exposure)
            else:
                mean = low = high = upper = None
                status = "incomplete"
            ensemble_rows.append(
                {
                    "scenario": scenario_name,
                    "species": species,
                    "rate_basis": OBSERVED_RATE_BASIS,
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
            propensity_values = [
                getattr(
                    member,
                    f"expected_lattice_origin_{species}_flux_from_propensity_mol_m2_s",
                )
                for member in members
            ]
            defined_propensities = [
                value for value in propensity_values if value is not None
            ]
            if len(defined_propensities) == len(propensity_values):
                mean, low, high = bootstrap_summary(
                    defined_propensities,
                    f"{scenario_name}:{species}:{PROPENSITY_RATE_BASIS}",
                )
                status = (
                    "estimated-positive"
                    if mean > 0.0
                    else "censored-nonpositive-propensity"
                )
            else:
                mean = low = high = None
                status = "incomplete"
            ensemble_rows.append(
                {
                    "scenario": scenario_name,
                    "species": species,
                    "rate_basis": PROPENSITY_RATE_BASIS,
                    "status": status,
                    "mean_mol_m2_s": _display(mean),
                    "ci95_low_mol_m2_s": _display(low),
                    "ci95_high_mol_m2_s": _display(high),
                    "poisson_zero_upper_95_mol_m2_s": "undefined",
                    "log10_mean_mol_m2_s": _display(
                        math.log10(mean) if mean is not None and mean > 0 else None
                    ),
                }
            )
    _write_csv_atomic(out_dir / "ensemble-rates.csv", ensemble_fields, ensemble_rows)

    stoichiometry_fields = [
        "scenario",
        "status",
        "si_al_lattice_release_ratio_dimensionless_mean",
        "si_al_lattice_release_ratio_dimensionless_ci95_low",
        "si_al_lattice_release_ratio_dimensionless_ci95_high",
        "lattice_si_release_events",
        "lattice_al_release_events",
    ]
    stoichiometry_rows: list[dict[str, object]] = []
    for scenario_name, members in by_scenario.items():
        ratios = [
            member.si_al_lattice_release_ratio_dimensionless for member in members
        ]
        defined = [value for value in ratios if value is not None]
        if len(defined) == len(ratios):
            mean, low, high = bootstrap_summary(
                defined, f"{scenario_name}:si_al_lattice_release_ratio_dimensionless"
            )
            status = "estimated"
        else:
            mean = low = high = None
            status = "undefined-nonpositive-lattice-release"
        stoichiometry_rows.append(
            {
                "scenario": scenario_name,
                "status": status,
                "si_al_lattice_release_ratio_dimensionless_mean": _display(mean),
                "si_al_lattice_release_ratio_dimensionless_ci95_low": _display(low),
                "si_al_lattice_release_ratio_dimensionless_ci95_high": _display(high),
                "lattice_si_release_events": sum(
                    member.lattice_si_release_events for member in members
                ),
                "lattice_al_release_events": sum(
                    member.lattice_al_release_events for member in members
                ),
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

    def combined_propensity_response(rate: ReplicaRate) -> float | None:
        if not rate.acceptance_passed or rate.propensity_origin_status != "origin-safe":
            return None
        si = rate.expected_lattice_origin_si_flux_from_propensity_mol_m2_s
        al = rate.expected_lattice_origin_al_flux_from_propensity_mol_m2_s
        if si is None or al is None:
            return None
        value = si + al
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(
                "combined expected lattice-origin propensity flux must be finite and nonnegative"
            )
        return value

    nominal_responses = {
        seed: combined_propensity_response(nominal_by_seed[seed]) for seed in seeds
    }
    nominal_defined = [
        value for value in nominal_responses.values() if value is not None
    ]
    nominal_mean = (
        statistics.fmean(nominal_defined)
        if len(nominal_defined) == len(nominal_responses)
        else None
    )
    nominal_status = (
        "estimated-positive"
        if nominal_mean is not None and nominal_mean > 0.0
        else "censored"
    )
    sensitivity_fields = [
        "rank",
        "family",
        "declared_response",
        "status",
        "max_abs_paired_delta_log10",
        "nominal_status",
        "nominal_mean_mol_m2_s",
        *[
            field
            for perturbation in PERTURBATIONS
            for field in (
                f"{perturbation}_status",
                f"{perturbation}_paired_replica_count",
                f"{perturbation}_perturbed_mean_mol_m2_s",
                f"{perturbation}_paired_mean_delta_mol_m2_s",
                f"{perturbation}_paired_delta_log10",
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
            "nominal_status": nominal_status,
            "nominal_mean_mol_m2_s": _display(nominal_mean),
        }
        responses: list[float] = []
        perturbation_statuses: list[str] = []
        for perturbation in PERTURBATIONS:
            scenario_name = f"{family}__{perturbation}"
            perturbed = {rate.seed: rate for rate in by_scenario[scenario_name]}
            pairs = [
                (nominal_responses[seed], combined_propensity_response(perturbed[seed]))
                for seed in seeds
            ]
            complete_pairs = [
                (before, after)
                for before, after in pairs
                if before is not None and after is not None
            ]
            pair_count = len(complete_pairs)
            row[f"{perturbation}_paired_replica_count"] = pair_count
            if pair_count != len(seeds):
                status = "censored"
                perturbed_mean = None
                paired_mean_delta = None
                delta_log10 = None
            else:
                before_values = [float(before) for before, _ in complete_pairs]
                after_values = [float(after) for _, after in complete_pairs]
                paired_linear_deltas = [
                    after - before
                    for before, after in zip(before_values, after_values, strict=True)
                ]
                before_mean = statistics.fmean(before_values)
                perturbed_mean = statistics.fmean(after_values)
                paired_mean_delta = statistics.fmean(paired_linear_deltas)
                reconstructed_after_mean = before_mean + paired_mean_delta
                if not math.isclose(
                    perturbed_mean,
                    reconstructed_after_mean,
                    rel_tol=1.0e-14,
                    abs_tol=0.0,
                ):
                    raise AssertionError(
                        "paired sensitivity arithmetic is inconsistent"
                    )
                if before_mean <= 0.0 or perturbed_mean <= 0.0:
                    status = "censored"
                    delta_log10 = None
                else:
                    status = "estimated"
                    delta_log10 = math.log10(perturbed_mean) - math.log10(before_mean)
                    responses.append(abs(delta_log10))
            row[f"{perturbation}_status"] = status
            row[f"{perturbation}_perturbed_mean_mol_m2_s"] = _display(perturbed_mean)
            row[f"{perturbation}_paired_mean_delta_mol_m2_s"] = _display(
                paired_mean_delta
            )
            row[f"{perturbation}_paired_delta_log10"] = _display(delta_log10)
            perturbation_statuses.append(status)
        if any(status != "estimated" for status in perturbation_statuses):
            row["status"] = "censored"
            row["max_abs_paired_delta_log10"] = "undefined"
        else:
            row["max_abs_paired_delta_log10"] = max(responses)
        sensitivity_rows.append(row)
    sensitivity_rows.sort(
        key=lambda row: (
            row["max_abs_paired_delta_log10"] == "undefined",
            -float(row["max_abs_paired_delta_log10"])
            if row["max_abs_paired_delta_log10"] != "undefined"
            else 0.0,
            str(row["family"]),
        )
    )
    next_rank = 1
    for row in sensitivity_rows:
        if row["status"] == "estimated":
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
        "acceptance_passed": outcome in {"steady-positive", "no-dissolution"},
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
        "observed_event_acceptance_basis": OBSERVED_EVENT_ACCEPTANCE_BASIS,
        "propensity_estimator_basis": PROPENSITY_ESTIMATOR_BASIS,
        "propensity_integral": PROPENSITY_INTEGRAL,
        "propensity_interpretation": PROPENSITY_INTERPRETATION,
        "propensity_reaction_index_mapping": {
            species: list(EXPECTED_REACTIONS).index(f"desorb-{species}")
            for species in ("si", "al")
        },
        "sensitivity_response": SENSITIVITY_RESPONSE,
        "sensitivity_statistic": SENSITIVITY_STATISTIC,
        "sensitivity_complete": len(sensitivity_rows) == len(FAMILY_REACTIONS)
        and {row["family"] for row in sensitivity_rows} == set(FAMILY_REACTIONS)
        and {row["declared_response"] for row in sensitivity_rows}
        == {SENSITIVITY_RESPONSE}
        and all(row["status"] in {"estimated", "censored"} for row in sensitivity_rows),
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
        "log10_lattice_si_release_flux_mol_m2_s",
        "log10_lattice_al_release_flux_mol_m2_s",
        "log10_expected_lattice_origin_si_flux_from_propensity_mol_m2_s",
        "log10_expected_lattice_origin_al_flux_from_propensity_mol_m2_s",
    ]
    rate_rows = _read_csv_exact(
        out_dir / "per-replica-rates.csv", rate_fields, expected_runs
    )
    if {
        (row["scenario"], row["replica"], row["seed"]) for row in rate_rows
    } != expected_identities:
        raise ValueError("per-replica rate identity coverage mismatch")
    for row in rate_rows:
        for species in ("si", "al"):
            value_field = (
                f"expected_lattice_origin_{species}_flux_from_propensity_mol_m2_s"
            )
            log_field = (
                f"log10_expected_lattice_origin_{species}_flux_from_propensity_mol_m2_s"
            )
            if row[value_field] == "undefined":
                if (
                    row[log_field] != "undefined"
                    or row["propensity_origin_status"] != "origin-contaminated"
                ):
                    raise ValueError(
                        "undefined propensity flux must have undefined log10"
                    )
                continue
            if row["propensity_origin_status"] != "origin-safe":
                raise ValueError("defined propensity flux must be origin-safe")
            value = float(row[value_field])
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(
                    "per-replica propensity flux must be finite and nonnegative"
                )
            if value == 0.0:
                if row[log_field] != "undefined":
                    raise ValueError("zero propensity flux must have undefined log10")
            elif not math.isclose(
                float(row[log_field]), math.log10(value), rel_tol=1.0e-14
            ):
                raise ValueError("per-replica propensity log10 mismatch")
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
        for basis in (OBSERVED_RATE_BASIS, PROPENSITY_RATE_BASIS)
    }:
        raise ValueError("ensemble rate identity coverage mismatch")
    for row in ensemble_rows:
        if row["rate_basis"] != PROPENSITY_RATE_BASIS:
            continue
        if (
            row["status"]
            not in {
                "estimated-positive",
                "censored-nonpositive-propensity",
                "incomplete",
            }
            or row["poisson_zero_upper_95_mol_m2_s"] != "undefined"
        ):
            raise ValueError("propensity ensemble status/schema mismatch")
        summary_fields = (
            "mean_mol_m2_s",
            "ci95_low_mol_m2_s",
            "ci95_high_mol_m2_s",
        )
        if row["status"] == "incomplete":
            if any(row[field] != "undefined" for field in summary_fields):
                raise ValueError("incomplete propensity ensemble must be undefined")
            continue
        summary = [float(row[field]) for field in summary_fields]
        if any(not math.isfinite(value) or value < 0.0 for value in summary):
            raise ValueError(
                "propensity ensemble values must be finite and nonnegative"
            )
        mean, low, high = summary
        if low > high:
            raise ValueError("propensity ensemble bootstrap band is reversed")
        if row["status"] == "estimated-positive":
            if mean <= 0.0 or not math.isclose(
                float(row["log10_mean_mol_m2_s"]),
                math.log10(mean),
                rel_tol=1.0e-14,
            ):
                raise ValueError("positive propensity ensemble log10 mismatch")
        elif mean != 0.0 or row["log10_mean_mol_m2_s"] != "undefined":
            raise ValueError(
                "nonpositive propensity ensemble must be zero and censored"
            )
    stoichiometry_rows = _read_csv_exact(
        out_dir / "stoichiometry.csv",
        [
            "scenario",
            "status",
            "si_al_lattice_release_ratio_dimensionless_mean",
            "si_al_lattice_release_ratio_dimensionless_ci95_low",
            "si_al_lattice_release_ratio_dimensionless_ci95_high",
            "lattice_si_release_events",
            "lattice_al_release_events",
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
            "max_abs_paired_delta_log10",
            "nominal_status",
            "nominal_mean_mol_m2_s",
            *[
                field
                for perturbation in PERTURBATIONS
                for field in (
                    f"{perturbation}_status",
                    f"{perturbation}_paired_replica_count",
                    f"{perturbation}_perturbed_mean_mol_m2_s",
                    f"{perturbation}_paired_mean_delta_mol_m2_s",
                    f"{perturbation}_paired_delta_log10",
                )
            ],
        ],
        len(FAMILY_REACTIONS),
    )
    if {row["family"] for row in sensitivity} != set(FAMILY_REACTIONS):
        raise ValueError("sensitivity family coverage mismatch")
    if {row["declared_response"] for row in sensitivity} != {SENSITIVITY_RESPONSE}:
        raise ValueError("sensitivity declared response mismatch")
    estimated_rows = [row for row in sensitivity if row["status"] == "estimated"]
    if {row["rank"] for row in estimated_rows} != {
        str(rank) for rank in range(1, len(estimated_rows) + 1)
    }:
        raise ValueError("estimated sensitivity rank coverage mismatch")
    for row in sensitivity:
        if row["status"] not in {"estimated", "censored"}:
            raise ValueError("sensitivity status mismatch")
        if row["status"] == "censored" and row["rank"] != "undefined":
            raise ValueError("censored sensitivity rank must be undefined")
        nominal_status = row["nominal_status"]
        nominal_value = row["nominal_mean_mol_m2_s"]
        if nominal_status == "estimated-positive":
            if not math.isfinite(float(nominal_value)) or float(nominal_value) <= 0.0:
                raise ValueError(
                    "nominal propensity sensitivity response must be positive"
                )
        elif nominal_status != "censored":
            raise ValueError("nominal propensity sensitivity status mismatch")
        finite_deltas: list[float] = []
        perturbation_statuses = []
        for perturbation in PERTURBATIONS:
            pair_count = int(row[f"{perturbation}_paired_replica_count"])
            if not 0 <= pair_count <= REPLICA_COUNT:
                raise ValueError("sensitivity paired-replica coverage mismatch")
            status = row[f"{perturbation}_status"]
            perturbation_statuses.append(status)
            delta = row[f"{perturbation}_paired_delta_log10"]
            paired_mean_delta = row[f"{perturbation}_paired_mean_delta_mol_m2_s"]
            perturbed_value = row[f"{perturbation}_perturbed_mean_mol_m2_s"]
            if status == "estimated":
                if pair_count != REPLICA_COUNT:
                    raise ValueError("sensitivity paired-replica coverage mismatch")
                perturbed_mean = float(perturbed_value)
                linear_delta = float(paired_mean_delta)
                if (
                    perturbed_mean <= 0.0
                    or not math.isfinite(perturbed_mean)
                    or not math.isfinite(linear_delta)
                    or not math.isfinite(float(delta))
                    or not math.isclose(
                        perturbed_mean,
                        float(nominal_value) + linear_delta,
                        rel_tol=1.0e-14,
                        abs_tol=0.0,
                    )
                ):
                    raise ValueError("paired sensitivity arithmetic mismatch")
                finite_deltas.append(abs(float(delta)))
            elif status == "censored":
                if delta != "undefined":
                    raise ValueError("censored sensitivity delta must be undefined")
                if pair_count != REPLICA_COUNT and (
                    perturbed_value != "undefined" or paired_mean_delta != "undefined"
                ):
                    raise ValueError("censored incomplete response must be undefined")
            else:
                raise ValueError("sensitivity perturbation status mismatch")
        fully_estimated = all(status == "estimated" for status in perturbation_statuses)
        if (row["status"] == "estimated") != fully_estimated:
            raise ValueError("sensitivity family classification mismatch")
        if fully_estimated:
            if not finite_deltas or not math.isclose(
                float(row["max_abs_paired_delta_log10"]),
                max(finite_deltas),
                rel_tol=1.0e-14,
            ):
                raise ValueError("sensitivity ranking magnitude mismatch")
        elif row["max_abs_paired_delta_log10"] != "undefined":
            raise ValueError("censored sensitivity magnitude must be undefined")
    provenance = _read_csv_exact(
        out_dir / "provenance-conversions.csv",
        PROVENANCE_FIELDS,
        len(REACTION_REGISTRY) + 10,
    )
    if [row["reaction"] for row in provenance[:22]] != [
        entry.name for entry in REACTION_REGISTRY
    ]:
        raise ValueError("provenance reaction coverage/order mismatch")
    conditions = {
        row["constant_name"]: row
        for row in provenance
        if row["record_type"] == "condition"
    }
    if (
        set(conditions)
        != {
            "temperature",
            "dissolved_cation_activity",
            "dissolved_cation_mu",
            "effective_consumed_cation_factor_298k",
            "simulation_steps",
            "report_every_steps",
        }
        or conditions["temperature"]["constant_value"] != "298.0"
        or conditions["dissolved_cation_activity"]["constant_value"] != "1e-30"
        or conditions["dissolved_cation_mu"]["constant_value"] != "-1.0"
        or "numerical open-flow sink"
        not in conditions["dissolved_cation_activity"]["rationale"]
        or conditions["simulation_steps"]["constant_value"] != "200000"
        or conditions["report_every_steps"]["constant_value"] != "10000"
    ):
        raise ValueError("thermodynamic condition provenance mismatch")
    if (
        provenance[-2]["record_type"] != "estimator"
        or provenance[-2]["observable_type"] != PROPENSITY_ESTIMATOR_BASIS
        or "not observed event release" not in provenance[-2]["rationale"]
        or provenance[-1]["record_type"] != "sensitivity-response"
        or provenance[-1]["constant_value"] != SENSITIVITY_RESPONSE
        or "no Poisson bound" not in provenance[-1]["rationale"]
    ):
        raise ValueError("propensity provenance records mismatch")
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
            "observed_event_acceptance_basis",
            "propensity_estimator_basis",
            "propensity_integral",
            "propensity_interpretation",
            "propensity_reaction_index_mapping",
            "sensitivity_response",
            "sensitivity_statistic",
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
        or verification["observed_event_acceptance_basis"]
        != OBSERVED_EVENT_ACCEPTANCE_BASIS
        or verification["propensity_estimator_basis"] != PROPENSITY_ESTIMATOR_BASIS
        or verification["propensity_integral"] != PROPENSITY_INTEGRAL
        or verification["propensity_interpretation"] != PROPENSITY_INTERPRETATION
        or verification["propensity_reaction_index_mapping"]
        != {
            species: list(EXPECTED_REACTIONS).index(f"desorb-{species}")
            for species in ("si", "al")
        }
        or verification["sensitivity_response"] != SENSITIVITY_RESPONSE
        or verification["sensitivity_statistic"] != SENSITIVITY_STATISTIC
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
        "schema": "a9-single-run-analysis-v3",
        "survey_tier": True,
        "acceptance_passed": gate.acceptance_passed,
        "outcome": gate.dissolution_outcome,
        "steady_state_status": gate.status,
        "observed_event_acceptance_basis": OBSERVED_EVENT_ACCEPTANCE_BASIS,
        "propensity_estimator_basis": PROPENSITY_ESTIMATOR_BASIS,
        "propensity_interpretation": PROPENSITY_INTERPRETATION,
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
