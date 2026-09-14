#!/usr/bin/env python3
"""Prepare, run, analyze, and smoke the A9b sensitivity campaign.

This is survey-tier platform testing, never production kinetics. The simulated
scope is exactly the pH 4 reservoir anchor: pH 3 and pH 5 appear only as
laboratory-comparison bounds.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from collections import defaultdict
from collections.abc import Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import tomllib

try:
    from a9b_artifacts import (
        canonical_json_bytes,
        read_jsonl,
        remove_file_durable,
        sha256_file,
        write_json_atomic,
        write_json_atomic_new,
    )
except ModuleNotFoundError:  # importlib-based unit tests do not add this directory
    _artifacts_spec = importlib.util.spec_from_file_location(
        "a9b_artifacts", Path(__file__).with_name("a9b_artifacts.py")
    )
    assert _artifacts_spec and _artifacts_spec.loader
    _artifacts = importlib.util.module_from_spec(_artifacts_spec)
    _artifacts_spec.loader.exec_module(_artifacts)
    canonical_json_bytes = _artifacts.canonical_json_bytes
    read_jsonl = _artifacts.read_jsonl
    remove_file_durable = _artifacts.remove_file_durable
    sha256_file = _artifacts.sha256_file
    write_json_atomic = _artifacts.write_json_atomic
    write_json_atomic_new = _artifacts.write_json_atomic_new

FIXED_SEEDS = (90401, 90403, 90407, 90409, 90413, 90419, 90421, 90427)
FAMILIES = (
    "siloxane-neutral",
    "sioal-si-neutral",
    "sioal-al-neutral",
    "connectivity-ladder",
    "al-o-al-analogue",
    "adsorption",
    "cation-desorption",
)
PERTURBATIONS = ("ea-minus-3", "ea-plus-3", "prefactor-x0.1", "prefactor-x10")
MAX_WORKERS = 4
THREAD_CAP = 16
WEIGHT_ESS_FRACTION_MIN = 0.1
CONTRIBUTION_ESS_MIN = 2.0
CONTRIBUTION_ESS_ROLE = "analysis_precision_ranking_gate_not_sampling_contract"
EQUIVALENCE_BAND_LOG10 = (-0.1, 0.1)
TOP_K = 3
AVOGADRO = 6.02214076e23
A2_TO_M2 = 1.0e-20
POPULATION_TOLERANCE = 0.05
BOOTSTRAPS = 2000
MANIFEST_SCHEMA = "a9b-sensitivity-campaign-v2"
ANALYSIS_SCHEMA = "a9b-sensitivity-analysis-v2"
PACKAGE_FILES = {
    "metadata.json",
    "initial.pgif.json",
    "events.jsonl",
    "likelihood.jsonl",
    "checkpoints.jsonl",
    "final.pgif.json",
    "receipt.json",
}
PETRA = Path(__file__).resolve().parents[1]
REPO = PETRA.parent
BASE_DECK = PETRA / "examples" / "kaolinite-approx.toml"
RESERVOIR_CONTRACT = PETRA / "examples" / "kaolinite-reservoirs.toml"
COMMITTED_RESERVOIR = (
    REPO / "docs" / "program" / "results" / "a9b-reservoir-origin-contract"
)
PATH_REPORT = (
    REPO
    / "docs"
    / "program"
    / "results"
    / "a9b-mechanism-reachability"
    / "path-report.json"
)
LINEAGE_SNAPSHOT = PATH_REPORT.with_name("initial-snapshot.pgif.json")
KINETICS_LEDGER = REPO / "kinetics-db" / "minerals" / "kaolinite.toml"

REACTION_FAMILIES: dict[str, tuple[str, ...]] = {
    "siloxane-neutral": ("R0-sio-si-hydrolysis", "R1-sio-si-condensation"),
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
    "al-o-al-analogue": ("R14-alohal-hydrolysis", "R15-alohal-condensation"),
    "adsorption": ("adsorb-al", "adsorb-si"),
    "cation-desorption": ("desorb-al", "desorb-si"),
}


@dataclass(frozen=True)
class ScenarioSpec:
    name: str
    family: str | None
    perturbation: str | None


def scenario_specs() -> tuple[ScenarioSpec, ...]:
    values = [ScenarioSpec("nominal", None, None)]
    for family in FAMILIES:
        for perturbation in PERTURBATIONS:
            values.append(
                ScenarioSpec(f"{family}__{perturbation}", family, perturbation)
            )
    if len(values) != 29 or len({value.name for value in values}) != 29:
        raise AssertionError("A9b scenario preregistration drifted")
    return tuple(values)


def required_identities(seeds: Sequence[int]) -> list[tuple[str, int]]:
    return [(scenario.name, seed) for scenario in scenario_specs() for seed in seeds]


def validate_identity_set(
    identities: Sequence[tuple[str, int]], *, smoke: bool
) -> None:
    if len(set(identities)) != len(identities):
        raise ValueError("duplicate campaign run identity")
    expected_seeds = FIXED_SEEDS[:1] if smoke else FIXED_SEEDS
    expected = set(required_identities(expected_seeds))
    actual = set(identities)
    if actual != expected:
        raise ValueError(
            f"missing/extra campaign identities: missing={sorted(expected - actual)}, "
            f"extra={sorted(actual - expected)}"
        )


def _load_script(name: str) -> Any:
    path = Path(__file__).with_name(f"{name}.py")
    spec = importlib.util.spec_from_file_location(f"a9b_dependency_{name}", path)
    if not spec or not spec.loader:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _finite(value: object, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{label} must be numeric")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def _logsumexp(values: Sequence[float]) -> float:
    if not values:
        raise ValueError("logsumexp requires data")
    maximum = max(values)
    if not math.isfinite(maximum):
        raise ValueError("log values must be finite")
    return maximum + math.log(sum(math.exp(value - maximum) for value in values))


def _ess_from_logs(log_values: Sequence[float]) -> float:
    if not log_values:
        raise ValueError("ESS requires data")
    if any(not math.isfinite(value) for value in log_values):
        raise ValueError("ESS log values must be finite")
    if all(value == log_values[0] for value in log_values):
        return float(len(log_values))
    maximum = max(log_values)
    scaled_weights = [math.exp(value - maximum) for value in log_values]
    scaled_total = sum(scaled_weights)
    scaled_square_total = sum(weight * weight for weight in scaled_weights)
    result = scaled_total * scaled_total / scaled_square_total
    return min(float(len(log_values)), result)


def _unnormalized_log_mean(
    log_weights: Sequence[float], values: Sequence[float]
) -> float | None:
    if len(log_weights) != len(values) or not values:
        raise ValueError("unnormalized mean requires paired nonempty data")
    terms = []
    for log_weight, value in zip(log_weights, values, strict=True):
        _finite(log_weight, "log weight")
        value = _finite(value, "observation")
        if value < 0:
            raise ValueError("observations must be nonnegative")
        if value > 0:
            terms.append(log_weight + math.log(value))
    if not terms:
        return None
    return _logsumexp(terms) - math.log(len(values))


def _positive_linear_from_log(log_value: float) -> tuple[float | None, str]:
    try:
        value = math.exp(log_value)
    except OverflowError:
        return None, "censored_numerical_overflow"
    if not math.isfinite(value):
        return None, "censored_numerical_overflow"
    if value == 0.0:
        return None, "censored_numerical_underflow"
    return value, "represented"


def _state_population_point(
    log_weights: Sequence[float], values: Sequence[float]
) -> dict[str, object]:
    log_mean = _unnormalized_log_mean(log_weights, values)
    if log_mean is None:
        return {"value": 0.0, "log_mean": None, "status": "represented"}
    value, status = _positive_linear_from_log(log_mean)
    return {"value": value, "log_mean": log_mean, "status": status}


def _bootstrap_log_interval(
    log_weights: Sequence[float], values: Sequence[float], seed: int
) -> list[float | None]:
    rng = random.Random(seed)
    count = len(values)
    samples = []
    for _ in range(BOOTSTRAPS):
        chosen = [rng.randrange(count) for _ in range(count)]
        samples.append(
            _unnormalized_log_mean(
                [log_weights[index] for index in chosen],
                [values[index] for index in chosen],
            )
        )
    samples.sort(key=lambda value: -math.inf if value is None else value)
    return [
        samples[round((BOOTSTRAPS - 1) * 0.025)],
        samples[round((BOOTSTRAPS - 1) * 0.975)],
    ]


def estimate_species(
    log_weights: Sequence[float], observations: Sequence[float], *, seed: int = 918273
) -> dict[str, object]:
    if len(log_weights) != len(observations) or not observations:
        raise ValueError("estimator requires paired nonempty replicas")
    logs = [_finite(value, "log weight") for value in log_weights]
    obs = [_finite(value, "observation") for value in observations]
    if any(value < 0 for value in obs):
        raise ValueError("observations must be nonnegative")
    weight_ess = _ess_from_logs(logs)
    weight_fraction = weight_ess / len(logs)
    base = {
        "estimator": "unnormalized_importance_sampling_log_domain",
        "replicas": len(logs),
        "weight_ess": weight_ess,
        "weight_ess_fraction": weight_fraction,
        "contribution_ess": None,
        "contribution_ess_gate_role": CONTRIBUTION_ESS_ROLE,
        "log_mean": None,
        "log_ci95": None,
        "mean": None,
        "ci95": None,
        "upper_95": None,
        "poisson_bound": "unavailable",
        "rankable": False,
    }
    if weight_fraction < WEIGHT_ESS_FRACTION_MIN:
        return {**base, "status": "censored_low_weight_ess"}
    if not any(value > 0 for value in obs):
        return {**base, "status": "censored_zero_biased_observations"}
    contribution_logs = [
        log_weight + math.log(value)
        for log_weight, value in zip(logs, obs, strict=True)
        if value > 0
    ]
    contribution_ess = _ess_from_logs(contribution_logs)
    base["contribution_ess"] = contribution_ess
    if contribution_ess < CONTRIBUTION_ESS_MIN:
        return {**base, "status": "censored_low_contribution_ess"}
    log_mean = _unnormalized_log_mean(logs, obs)
    assert log_mean is not None
    log_interval = _bootstrap_log_interval(logs, obs, seed)
    base["log_mean"] = log_mean
    base["log_ci95"] = log_interval
    mean, mean_status = _positive_linear_from_log(log_mean)
    interval: list[float] = []
    interval_statuses = []
    for endpoint in log_interval:
        if endpoint is None:
            interval.append(0.0)
            continue
        value, status = _positive_linear_from_log(endpoint)
        interval_statuses.append(status)
        if value is not None:
            interval.append(value)
    numerical_statuses = [mean_status, *interval_statuses]
    if "censored_numerical_overflow" in numerical_statuses:
        return {**base, "status": "censored_numerical_overflow"}
    if "censored_numerical_underflow" in numerical_statuses:
        return {**base, "status": "censored_numerical_underflow"}
    assert mean is not None and len(interval) == 2
    return {
        **base,
        "status": "estimated",
        "mean": mean,
        "ci95": interval,
        "rankable": True,
    }


def _apply_run_and_stationarity_gates(
    estimate: dict[str, object],
    *,
    runs_complete: bool,
    stationarity_passed: bool,
) -> dict[str, object]:
    if not runs_complete:
        return {
            **estimate,
            "status": "censored_incomplete_replicas",
            "mean": None,
            "ci95": None,
            "upper_95": None,
            "rankable": False,
        }
    if estimate["rankable"] and not stationarity_passed:
        return {
            **estimate,
            "status": "censored_nonstationary_state_distribution",
            "mean": None,
            "ci95": None,
            "upper_95": None,
            "rankable": False,
        }
    return estimate


def paired_response(
    nominal: Mapping[int, tuple[float, float]],
    perturbed: Mapping[int, tuple[float, float]],
    *,
    nominal_rankable: bool,
    perturbed_rankable: bool,
    seed: int = 192837,
) -> dict[str, object]:
    if set(nominal) != set(perturbed):
        return {
            "status": "incomplete_pairs",
            "same_seed_linear_deltas": [],
            "mean_linear_delta": None,
            "mean_delta_sign": None,
            "log_abs_mean_delta_mol_m2_s": None,
            "delta_log10": None,
            "delta_log10_ci95": None,
        }
    seeds = sorted(nominal)

    def weighted_log_value(pair: tuple[float, float]) -> float | None:
        log_weight = _finite(pair[0], "paired log weight")
        observation = _finite(pair[1], "paired observation")
        if observation < 0:
            raise ValueError("paired observations must be nonnegative")
        if observation == 0:
            return None
        return log_weight + math.log(observation)

    def signed_log_delta(
        left: float | None, right: float | None
    ) -> tuple[int, float | None]:
        if left is None and right is None:
            return 0, None
        if left is None:
            return 1, right
        if right is None:
            return -1, left
        if left == right:
            return 0, None
        high, low = (right, left) if right > left else (left, right)
        return (1 if right > left else -1), high + math.log(-math.expm1(low - high))

    def representable_signed_value(sign: int, log_abs: float | None) -> float | None:
        if sign == 0:
            return 0.0
        assert log_abs is not None
        try:
            magnitude = math.exp(log_abs)
        except OverflowError:
            return None
        return sign * magnitude if math.isfinite(magnitude) else None

    nominal_logs = [weighted_log_value(nominal[item]) for item in seeds]
    perturbed_logs = [weighted_log_value(perturbed[item]) for item in seeds]
    signed_deltas = [
        signed_log_delta(left, right)
        for left, right in zip(nominal_logs, perturbed_logs, strict=True)
    ]
    nonzero_delta_logs = [
        value for sign, value in signed_deltas if sign and value is not None
    ]
    if not nonzero_delta_logs:
        mean_delta = 0.0
        mean_delta_sign = 0
        mean_delta_log_abs = None
    else:
        scale = max(nonzero_delta_logs)
        scaled_sum = sum(
            sign * math.exp(log_abs - scale)
            for sign, log_abs in signed_deltas
            if sign and log_abs is not None
        )
        if scaled_sum == 0.0:
            mean_delta = 0.0
            mean_delta_sign = 0
            mean_delta_log_abs = None
        else:
            mean_delta_sign = 1 if scaled_sum > 0 else -1
            mean_delta_log_abs = (
                scale + math.log(abs(scaled_sum)) - math.log(len(seeds))
            )
            mean_delta = representable_signed_value(mean_delta_sign, mean_delta_log_abs)
    result: dict[str, object] = {
        "status": "unranked_censored",
        "same_seed_linear_deltas": [
            {
                "seed": item,
                "delta_mol_m2_s": representable_signed_value(sign, log_abs),
                "delta_sign": sign,
                "log_abs_delta_mol_m2_s": log_abs,
            }
            for item, (sign, log_abs) in zip(seeds, signed_deltas, strict=True)
        ],
        "mean_linear_delta": mean_delta,
        "mean_delta_sign": mean_delta_sign,
        "log_abs_mean_delta_mol_m2_s": mean_delta_log_abs,
        "delta_log10": None,
        "delta_log10_ci95": None,
    }

    def log_mean(
        values: Mapping[int, tuple[float, float]], chosen: Sequence[int]
    ) -> float | None:
        terms = [weighted_log_value(values[item]) for item in chosen]
        positive = [term for term in terms if term is not None]
        return None if not positive else _logsumexp(positive) - math.log(len(chosen))

    nominal_log_mean = log_mean(nominal, seeds)
    perturbed_log_mean = log_mean(perturbed, seeds)
    if (
        not nominal_rankable
        or not perturbed_rankable
        or nominal_log_mean is None
        or perturbed_log_mean is None
    ):
        return result
    result["status"] = "rankable"
    result["delta_log10"] = (perturbed_log_mean - nominal_log_mean) / math.log(10.0)
    rng = random.Random(seed)
    boot = []
    for _ in range(BOOTSTRAPS):
        chosen = [seeds[rng.randrange(len(seeds))] for _ in seeds]
        n_log_mean = log_mean(nominal, chosen)
        p_log_mean = log_mean(perturbed, chosen)
        if n_log_mean is not None and p_log_mean is not None:
            boot.append((p_log_mean - n_log_mean) / math.log(10.0))
    if len(boot) == BOOTSTRAPS:
        boot.sort()
        result["delta_log10_ci95"] = [
            boot[round((BOOTSTRAPS - 1) * 0.025)],
            boot[round((BOOTSTRAPS - 1) * 0.975)],
        ]
    return result


def interval_demonstrates_irrelevance(interval: Sequence[float] | None) -> bool:
    return bool(
        interval is not None
        and len(interval) == 2
        and all(math.isfinite(value) for value in interval)
        and EQUIVALENCE_BAND_LOG10[0]
        <= interval[0]
        <= interval[1]
        <= EQUIVALENCE_BAND_LOG10[1]
    )


def population_stationarity(
    rows: Sequence[Mapping[str, float]],
    *,
    fixed_kind_totals: Mapping[str, float] | None = None,
) -> dict[str, object]:
    if len(rows) < 5:
        return {"passed": False, "status": "insufficient_checkpoints", "kinds": {}}
    states = set(rows[0])
    if any(set(row) != states for row in rows):
        raise ValueError("state checkpoint tables disagree")
    kinds: dict[str, list[str]] = defaultdict(list)
    for state in sorted(states):
        if "." not in state:
            raise ValueError("state names must be kind-qualified")
        kinds[state.partition(".")[0]].append(state)
    diagnostics: dict[str, object] = {}
    passed = True
    for kind, names in kinds.items():
        if fixed_kind_totals is None:
            totals = [sum(_finite(row[name], name) for name in names) for row in rows]
        else:
            declared = _finite(fixed_kind_totals[kind], f"{kind} fixed kind total")
            totals = [declared] * len(rows)
        if any(total <= 0 for total in totals):
            diagnostics[kind] = {"status": "invalid_zero_kind_total", "states": {}}
            passed = False
            continue
        state_diagnostics = {}
        kind_passed = True
        for name in names:
            fractions = [
                row[name] / total for row, total in zip(rows, totals, strict=True)
            ]
            tail = fractions[-5:]
            span = max(tail) - min(tail)
            trend = tail[-1] - tail[0]
            state_passed = (
                span <= POPULATION_TOLERANCE and abs(trend) <= POPULATION_TOLERANCE
            )
            kind_passed &= state_passed
            state_diagnostics[name] = {
                "tail_fraction_range": span,
                "tail_endpoint_change": trend,
                "status": "stationary" if state_passed else "nonstationary",
            }
        diagnostics[kind] = {
            "status": "stationary" if kind_passed else "nonstationary",
            "states": state_diagnostics,
        }
        passed &= kind_passed
    return {
        "passed": passed,
        "status": "stationary" if passed else "nonstationary",
        "tolerance_fraction": POPULATION_TOLERANCE,
        "kinds": diagnostics,
    }


def lab_comparison_scope(
    simulated_species_rates: Mapping[str, float | None],
) -> dict[str, object]:
    # P&K 25 C acid + neutral terms. pH 3/5 are comparison bounds only.
    rates = {}
    for ph in (3, 4, 5):
        rates[str(ph)] = 10**-11.31 * (10 ** (-ph)) ** 0.777 + 10**-13.18
    comparison: dict[str, object] = {
        "simulated_ph": 4,
        "comparison_bound_ph": [3, 5],
        "ph3_ph5_role": "laboratory_comparison_bounds_only_not_simulated",
        "laboratory_reference": "Palandri-Kharaka acid_plus_neutral_298K_approximation",
        "laboratory_formula_unit_rates_mol_m2_s": rates,
        "laboratory_rate_basis": "kaolinite_formula_units_mol_m2_s",
        "kaolinite_stoichiometry": "Al2Si2O5(OH)4",
        "formula_unit_conversion": "min(Si_mol_flux/2, Al_mol_flux/2)",
        "simulated_species_rates_mol_m2_s": {
            "Si": simulated_species_rates.get("Si"),
            "Al": simulated_species_rates.get("Al"),
        },
        "simulated_kaolinite_formula_unit_rate_mol_m2_s": None,
        "formula_unit_rate_status": "censored_species_rate_unavailable",
        "log10_gaps_vs_laboratory": {
            "ph3_bound_not_simulated": None,
            "ph4_simulated": None,
            "ph5_bound_not_simulated": None,
        },
    }
    si_rate = simulated_species_rates.get("Si")
    al_rate = simulated_species_rates.get("Al")
    if (
        isinstance(si_rate, (int, float))
        and not isinstance(si_rate, bool)
        and isinstance(al_rate, (int, float))
        and not isinstance(al_rate, bool)
        and si_rate > 0
        and al_rate > 0
        and math.isfinite(si_rate)
        and math.isfinite(al_rate)
    ):
        formula_unit_rate = min(si_rate / 2.0, al_rate / 2.0)
        comparison["simulated_kaolinite_formula_unit_rate_mol_m2_s"] = formula_unit_rate
        comparison["formula_unit_rate_status"] = "represented"
        comparison["log10_gaps_vs_laboratory"] = {
            "ph3_bound_not_simulated": math.log10(formula_unit_rate)
            - math.log10(rates["3"]),
            "ph4_simulated": math.log10(formula_unit_rate) - math.log10(rates["4"]),
            "ph5_bound_not_simulated": math.log10(formula_unit_rate)
            - math.log10(rates["5"]),
        }
    return comparison


def _pending_analysis(finalized: dict[str, Any]) -> dict[str, Any]:
    """Build the analyzer-visible report without publishing verified conclusions."""
    pending = copy.deepcopy(finalized)
    pending["preliminary_conclusions"] = {
        "status": "unverified_preliminary",
        "accepted_conclusions_withheld_until": "independent_verification_receipt",
    }
    pending.update(
        {
            "artifact_verified": False,
            "verification_status": "pending_independent_verification",
            "conclusions_status": "unverified_preliminary",
            "scientific_acceptance": "unverified_pending_independent_verification",
            "ordinal_ranking_complete": False,
            "ordinal_ranking_status": "unverified_pending_independent_verification",
            "top3_among_rankable": [],
            "families_demonstrated_irrelevant": [],
            "irrelevance_statement": "unverified pending independent verification",
        }
    )
    for scenario in pending["scenarios"].values():
        scenario["verification_status"] = "unverified_preliminary"
        scenario["combined_si_al"]["verification_status"] = "unverified_preliminary"
        for estimate in scenario["species"].values():
            estimate["verification_status"] = "unverified_preliminary"
    for family in pending["family_ranking"]:
        family["verification_status"] = "unverified_preliminary"
        family["demonstrated_irrelevant"] = False
        if "ordinal_rank_among_rankable" in family:
            family.pop("ordinal_rank_among_rankable")
        for response in family["responses"].values():
            response["verification_status"] = "unverified_preliminary"
    pending["laboratory_comparison"]["verification_status"] = "unverified_preliminary"
    return pending


def _validate_upstream_inputs() -> dict[str, object]:
    report = json.loads(PATH_REPORT.read_text(encoding="utf-8"))
    if (
        report.get("schema") != "a9b-mechanism-reachability-v2"
        or report.get("target_count") != 320
        or report.get("reachable_count") != 260
        or report.get("unreachable_count") != 60
        or len(report.get("targets", [])) != 320
        or sum(bool(row.get("reachable")) for row in report.get("targets", [])) != 260
        or report.get("snapshot_sha256") != sha256_file(LINEAGE_SNAPSHOT)
    ):
        raise ValueError("mechanism reachability census/hash contract failed")
    evidence_path = COMMITTED_RESERVOIR / "kaolinite-ph4-evidence.json"
    deck_path = COMMITTED_RESERVOIR / "kaolinite-ph4.toml"
    reservoir = _load_script("reservoir_origin_contract")
    evidence = reservoir.verify_evidence(deck_path, evidence_path)
    if evidence["selected_profile"]["ph"] != 4:
        raise ValueError("committed reservoir evidence is not pH 4")
    return {
        "path_report_sha256": sha256_file(PATH_REPORT),
        "lineage_snapshot_sha256": sha256_file(LINEAGE_SNAPSHOT),
        "committed_ph4_deck_sha256": sha256_file(deck_path),
        "committed_ph4_evidence_sha256": sha256_file(evidence_path),
    }


def _validate_existing_campaign(
    root: Path, manifest: Mapping[str, Any], *, smoke: bool
) -> None:
    immutable_keys = (
        "schema",
        "survey_tier_platform_test",
        "not_production",
        "smoke",
        "simulated_ph",
        "ph3_ph5_scope",
        "temperature_k",
        "seeds",
        "production_fixed_seeds",
        "replicas_per_scenario",
        "scenario_count",
        "required_run_count",
        "top_k_preregistered",
        "equivalence_band_log10",
        "weight_ess_fraction_min",
        "species_contribution_ess_min",
        "species_contribution_ess_role",
        "state_population_stationarity",
        "estimand",
        "run_contract",
        "inputs",
        "scenarios",
        "required_runs",
    )
    optional_keys = {"workers", "runner_sha256", "completed_run_count"}
    required_keys = set(immutable_keys) | {"status"}
    if set(manifest) != required_keys | (set(manifest) & optional_keys):
        raise ValueError("existing campaign manifest key set is not canonical")
    if manifest["status"] not in {"prepared", "running", "failed", "raw_complete"}:
        raise ValueError("existing campaign manifest status is not canonical")
    with tempfile.TemporaryDirectory(prefix=".a9b-validate-", dir=root.parent) as tmp:
        fresh_root = Path(tmp)
        fresh = prepare_campaign(fresh_root, smoke=smoke)
        for key in immutable_keys:
            if manifest.get(key) != fresh.get(key):
                raise ValueError(f"existing campaign manifest drift in {key}")
        relative_files = {
            manifest["inputs"][name]
            for name in (
                "base_deck",
                "reservoir_contract",
                "path_report",
                "lineage_snapshot",
                "ph4_evidence",
                "ph4_deck",
            )
        }
        relative_files.update(row["deck"] for row in manifest["scenarios"])
        for relative in sorted(relative_files):
            existing = root / relative
            regenerated = fresh_root / relative
            if (
                not existing.is_file()
                or existing.read_bytes() != regenerated.read_bytes()
            ):
                raise ValueError(f"prepared campaign artifact drift at {relative}")


def _prepare_campaign_in_place(
    root: Path, *, smoke: bool, upstream: Mapping[str, object]
) -> dict[str, object]:
    manifest_path = root / "manifest.json"
    inputs = root / "inputs"
    decks = root / "decks"
    for directory in (inputs, decks, root / "runs", root / "logs", root / "receipts"):
        directory.mkdir(parents=True, exist_ok=True)
    base_copy = inputs / BASE_DECK.name
    contract_copy = inputs / RESERVOIR_CONTRACT.name
    path_copy = inputs / "path-report.json"
    snapshot_copy = inputs / "initial-snapshot.pgif.json"
    for source, destination in (
        (BASE_DECK, base_copy),
        (RESERVOIR_CONTRACT, contract_copy),
        (PATH_REPORT, path_copy),
        (LINEAGE_SNAPSHOT, snapshot_copy),
    ):
        destination.write_bytes(source.read_bytes())
    reservoir = _load_script("reservoir_origin_contract")
    ph4_deck = inputs / "kaolinite-ph4.toml"
    ph4_evidence = inputs / "kaolinite-ph4-evidence.json"
    reservoir.materialize(base_copy, contract_copy, 4, ph4_deck, ph4_evidence)
    reservoir.verify_evidence(ph4_deck, ph4_evidence)
    if (
        sha256_file(ph4_deck) != upstream["committed_ph4_deck_sha256"]
        or sha256_file(ph4_evidence) != upstream["committed_ph4_evidence_sha256"]
    ):
        raise ValueError("fresh pH4 materialization differs from merged evidence")

    closure = _load_script("approximate_rate_closure")
    contract = closure.validate_deck(ph4_deck)
    closure_scenarios = closure.scenarios()
    if [item.name for item in closure_scenarios] != [
        item.name for item in scenario_specs()
    ]:
        raise ValueError(
            "reused A9 scenario generator disagrees with A9b preregistration"
        )
    scenario_rows = []
    for scenario in scenario_specs():
        text = (
            contract.text
            if scenario.family is None
            else closure.perturb_deck(contract, scenario.family, scenario.perturbation)
        )
        path = decks / f"{scenario.name}.toml"
        path.write_text(text, encoding="utf-8")
        scenario_rows.append(
            {
                **asdict(scenario),
                "deck": str(path.relative_to(root)),
                "deck_sha256": sha256_file(path),
            }
        )
    seeds = FIXED_SEEDS[:1] if smoke else FIXED_SEEDS
    identities = required_identities(seeds)
    validate_identity_set(identities, smoke=smoke)
    manifest: dict[str, object] = {
        "schema": MANIFEST_SCHEMA,
        "status": "prepared",
        "survey_tier_platform_test": True,
        "not_production": True,
        "smoke": smoke,
        "temperature_k": 298.0,
        "simulated_ph": 4,
        "ph3_ph5_scope": "laboratory_comparison_bounds_only_not_simulated",
        "seeds": list(seeds),
        "production_fixed_seeds": list(FIXED_SEEDS),
        "scenario_count": 29,
        "replicas_per_scenario": len(seeds),
        "required_run_count": len(identities),
        "top_k_preregistered": TOP_K,
        "equivalence_band_log10": list(EQUIVALENCE_BAND_LOG10),
        "weight_ess_fraction_min": WEIGHT_ESS_FRACTION_MIN,
        "species_contribution_ess_min": CONTRIBUTION_ESS_MIN,
        "species_contribution_ess_role": CONTRIBUTION_ESS_ROLE,
        "state_population_stationarity": {
            "required": "all_states_within_every_kind_for_ensemble_and_each_replica",
            "tail_checkpoint_count": 5,
            "max_fractional_span_and_endpoint_change": POPULATION_TOLERANCE,
        },
        "estimand": {
            "whole_finite_deck_centers": 320,
            "reachable_centers": 260,
            "topology_no_go_centers": 60,
            "primary": "whole_finite_deck_including_topology_no_go_in_denominator_and_area",
            "secondary": "reachable_center_census",
        },
        "inputs": {
            "base_deck": str(base_copy.relative_to(root)),
            "base_deck_sha256": sha256_file(base_copy),
            "reservoir_contract": str(contract_copy.relative_to(root)),
            "reservoir_contract_sha256": sha256_file(contract_copy),
            "ph4_deck": str(ph4_deck.relative_to(root)),
            "ph4_deck_sha256": sha256_file(ph4_deck),
            "ph4_evidence": str(ph4_evidence.relative_to(root)),
            "ph4_evidence_sha256": sha256_file(ph4_evidence),
            "path_report": str(path_copy.relative_to(root)),
            "path_report_sha256": sha256_file(path_copy),
            "lineage_snapshot": str(snapshot_copy.relative_to(root)),
            "lineage_snapshot_sha256": sha256_file(snapshot_copy),
        },
        "scenarios": scenario_rows,
        "required_runs": [
            {"scenario": scenario, "seed": seed} for scenario, seed in identities
        ],
        "run_contract": {
            "horizon_s": 1.0e-12 if smoke else 1.0e-6,
            "checkpoints": 3 if smoke else 21,
            "bias_factor": 1.0e6,
            "max_events": 10_000 if smoke else 1_000_000,
        },
    }
    write_json_atomic_new(manifest_path, manifest)
    return manifest


def prepare_campaign(root: Path, *, smoke: bool) -> dict[str, object]:
    root = root.resolve()
    manifest_path = root / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if (
            manifest.get("schema") != MANIFEST_SCHEMA
            or manifest.get("smoke") is not smoke
        ):
            raise ValueError("existing campaign manifest is incompatible")
        _validate_existing_campaign(root, manifest, smoke=smoke)
        return manifest
    if root.exists() and any(root.iterdir()):
        raise ValueError("new campaign root must be absent or empty")
    root.parent.mkdir(parents=True, exist_ok=True)
    upstream = _validate_upstream_inputs()
    staging = Path(
        tempfile.mkdtemp(prefix=f".{root.name}.prepare-", dir=root.parent)
    ).resolve()
    try:
        manifest = _prepare_campaign_in_place(staging, smoke=smoke, upstream=upstream)
        if root.exists():
            root.rmdir()
        os.replace(staging, root)
        parent_descriptor = os.open(root.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(parent_descriptor)
        finally:
            os.close(parent_descriptor)
        return manifest
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def _worker_environment(workers: int) -> dict[str, str]:
    if not 1 <= workers <= MAX_WORKERS:
        raise ValueError(f"workers must be in [1,{MAX_WORKERS}]")
    threads = max(1, THREAD_CAP // workers)
    env = os.environ.copy()
    for variable in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "RAYON_NUM_THREADS",
    ):
        env[variable] = str(threads)
    return env


def _receipt_valid(root: Path, receipt_path: Path) -> bool:
    try:
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        package = root / receipt["package"]
        if receipt.get("schema") != "a9b-job-receipt-v2" or not package.is_dir():
            return False
        for name, digest in receipt["package_files_sha256"].items():
            if sha256_file(package / name) != digest:
                return False
        return set(receipt["package_files_sha256"]) == PACKAGE_FILES
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
        return False


def _run_job(
    root: Path,
    runner: Path,
    scenario: str,
    seed: int,
    deck: Path,
    inputs: Mapping[str, str],
    contract: Mapping[str, object],
    env: Mapping[str, str],
) -> dict[str, object]:
    package = root / "runs" / scenario / f"seed-{seed}"
    receipt_path = root / "receipts" / scenario / f"seed-{seed}.json"
    if receipt_path.exists() and _receipt_valid(root, receipt_path):
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if (
            receipt.get("scenario") != scenario
            or receipt.get("seed") != seed
            or Path(receipt.get("package", ""))
            != Path("runs") / scenario / f"seed-{seed}"
            or receipt.get("runner_sha256") != sha256_file(runner)
            or receipt.get("deck_sha256") != sha256_file(deck)
        ):
            raise ValueError(
                f"existing receipt identity/input mismatch: {receipt_path}"
            )
        return receipt
    if package.exists():
        raise ValueError(f"unreceipted/nonmatching package blocks resume: {package}")
    attempt = 1
    log_dir = root / "logs" / scenario
    log_dir.mkdir(parents=True, exist_ok=True)
    while (log_dir / f"seed-{seed}-attempt-{attempt}.log").exists():
        attempt += 1
    log_path = log_dir / f"seed-{seed}-attempt-{attempt}.log"
    command = [
        "nice",
        "-n",
        "10",
        str(runner),
        str(deck),
        str(root / inputs["path_report"]),
        str(root / inputs["lineage_snapshot"]),
        str(package),
        "--seed",
        str(seed),
        "--horizon",
        str(contract["horizon_s"]),
        "--checkpoints",
        str(contract["checkpoints"]),
        "--bias-factor",
        str(contract["bias_factor"]),
        "--max-events",
        str(contract["max_events"]),
    ]
    started = time.monotonic()
    with log_path.open("x", encoding="utf-8") as log:
        subprocess.run(
            command,
            cwd=PETRA,
            env=dict(env),
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
            shell=False,
        )
        log.flush()
        os.fsync(log.fileno())
    if (
        not package.is_dir()
        or {path.name for path in package.iterdir()} != PACKAGE_FILES
    ):
        raise ValueError(f"runner package is incomplete: {package}")
    receipt: dict[str, object] = {
        "schema": "a9b-job-receipt-v2",
        "scenario": scenario,
        "seed": seed,
        "attempt": attempt,
        "command": command,
        "package": str(package.relative_to(root)),
        "log": str(log_path.relative_to(root)),
        "elapsed_seconds": time.monotonic() - started,
        "runner_sha256": sha256_file(runner),
        "deck_sha256": sha256_file(deck),
        "log_sha256": sha256_file(log_path),
        "package_files_sha256": {
            path.name: sha256_file(path) for path in sorted(package.iterdir())
        },
    }
    write_json_atomic_new(receipt_path, receipt)
    return receipt


def run_campaign(root: Path, runner: Path, *, workers: int) -> dict[str, object]:
    root = root.resolve()
    runner = runner.resolve()
    if not runner.is_file():
        raise ValueError(f"runner does not exist: {runner}")
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("not an A9b campaign")
    smoke = bool(manifest["smoke"])
    _validate_existing_campaign(root, manifest, smoke=smoke)
    identities = [(row["scenario"], row["seed"]) for row in manifest["required_runs"]]
    validate_identity_set(identities, smoke=smoke)
    scenario_decks = {row["name"]: root / row["deck"] for row in manifest["scenarios"]}
    runner_digest = sha256_file(runner)
    existing_receipts = any((root / "receipts").glob("*/*.json"))
    if (
        existing_receipts
        and manifest.get("runner_sha256") is not None
        and manifest.get("runner_sha256") != runner_digest
    ):
        raise ValueError(
            "cannot resume existing receipts with a different runner binary"
        )
    manifest["status"] = "running"
    manifest["workers"] = workers
    manifest["runner_sha256"] = runner_digest
    write_json_atomic(manifest_path, manifest)
    env = _worker_environment(workers)
    receipts = []
    lock = threading.Lock()
    try:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(
                    _run_job,
                    root,
                    runner,
                    scenario,
                    seed,
                    scenario_decks[scenario],
                    manifest["inputs"],
                    manifest["run_contract"],
                    env,
                ): (scenario, seed)
                for scenario, seed in identities
            }
            for future in as_completed(futures):
                receipt = future.result()
                with lock:
                    receipts.append(receipt)
    except BaseException:
        manifest["status"] = "failed"
        write_json_atomic(manifest_path, manifest)
        raise
    if len(receipts) != len(identities):
        raise AssertionError("executor lost a run receipt")
    manifest["status"] = "raw_complete"
    manifest["completed_run_count"] = len(receipts)
    write_json_atomic(manifest_path, manifest)
    return manifest


def _close(left: float, right: float, *, tolerance: float = 2e-12) -> bool:
    return math.isclose(left, right, rel_tol=tolerance, abs_tol=1e-18)


def _pgif_document(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if set(payload) != {"pgif", "meta", "nodes", "edges"} or payload.get("pgif") != 1:
        raise ValueError(f"{path}: invalid canonical PGIF v1 envelope")
    meta = payload["meta"]
    if (
        set(meta) != {"directed", "kind", "petra", "producer"}
        or meta["directed"] is not False
        or meta["kind"] != "kmc-lattice"
        or meta["producer"] != "petra"
        or set(meta["petra"])
        != {"deck", "state_types", "states", "step", "temperature", "time"}
    ):
        raise ValueError(f"{path}: PGIF graph metadata is not canonical Petra metadata")
    petra = meta["petra"]
    states = petra["states"]
    state_types = petra["state_types"]
    step = petra["step"]
    time = _finite(petra["time"], "PGIF time")
    _finite(petra["temperature"], "PGIF temperature")
    if (
        not isinstance(petra["deck"], str)
        or not petra["deck"]
        or not isinstance(states, list)
        or not states
        or any(not isinstance(value, str) for value in states)
        or len(states) != len(set(states))
        or not isinstance(state_types, list)
        or len(state_types) != len(states)
        or any(not isinstance(value, str) for value in state_types)
        or type(step) is not int
        or step < 0
        or time < 0
    ):
        raise ValueError(f"{path}: PGIF Petra metadata values are malformed")

    nodes = payload["nodes"]
    if set(nodes) != {"count", "columns"} or type(nodes["count"]) is not int:
        raise ValueError(f"{path}: PGIF node table is malformed")
    count = nodes["count"]
    columns = nodes["columns"]
    expected_node_types = {
        "frozen": "bool",
        "kind": "categorical",
        "state": "categorical",
        "type": "categorical",
        "x": "f32",
        "y": "f32",
        "z": "f32",
    }
    if count < 0 or set(columns) != set(expected_node_types):
        raise ValueError(f"{path}: PGIF node columns are not canonical")
    for name, expected_type in expected_node_types.items():
        column = columns[name]
        expected_keys = (
            {"type", "data", "dict"}
            if expected_type == "categorical"
            else {"type", "data"}
        )
        if (
            not isinstance(column, dict)
            or set(column) != expected_keys
            or column.get("type") != expected_type
            or not isinstance(column.get("data"), list)
            or len(column["data"]) != count
        ):
            raise ValueError(f"{path}: PGIF {name} column is malformed")
        if expected_type == "categorical":
            dictionary = column["dict"]
            if (
                not isinstance(dictionary, list)
                or not dictionary
                or any(not isinstance(value, str) for value in dictionary)
                or len(dictionary) != len(set(dictionary))
                or any(
                    type(value) is not int or not 0 <= value < len(dictionary)
                    for value in column["data"]
                )
            ):
                raise ValueError(
                    f"{path}: PGIF {name} categorical values are malformed"
                )
    if (
        columns["state"]["dict"] != states
        or columns["frozen"]["data"]
        != [bool(value) for value in columns["frozen"]["data"]]
        or any(type(value) is not bool for value in columns["frozen"]["data"])
        or any(
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(float(value))
            for name in ("x", "y", "z")
            for value in columns[name]["data"]
        )
    ):
        raise ValueError(f"{path}: PGIF state/frozen/coordinate values are malformed")
    state_data = columns["state"]["data"]
    type_dict = columns["type"]["dict"]
    kind_dict = columns["kind"]["dict"]
    for site, state_id in enumerate(state_data):
        state_name = states[state_id]
        if (
            state_types[state_id] not in type_dict
            or columns["type"]["data"][site] != type_dict.index(state_types[state_id])
            or "." not in state_name
            or state_name.partition(".")[0] not in kind_dict
            or columns["kind"]["data"][site]
            != kind_dict.index(state_name.partition(".")[0])
        ):
            raise ValueError(
                f"{path}: PGIF type/kind disagrees with state at node {site}"
            )

    edges = payload["edges"]
    if (
        set(edges) != {"count", "src", "dst", "columns"}
        or type(edges["count"]) is not int
    ):
        raise ValueError(f"{path}: PGIF edge table is malformed")
    edge_count = edges["count"]
    seam = edges["columns"].get("seam") if isinstance(edges["columns"], dict) else None
    if (
        edge_count < 0
        or set(edges["columns"]) != {"seam"}
        or not isinstance(seam, dict)
        or set(seam) != {"type", "data"}
        or seam.get("type") != "bool"
        or any(
            not isinstance(values, list) or len(values) != edge_count
            for values in (edges["src"], edges["dst"], seam["data"])
        )
        or any(type(value) is not bool for value in seam["data"])
        or any(
            type(value) is not int or not 0 <= value < count
            for values in (edges["src"], edges["dst"])
            for value in values
        )
        or any(
            source == destination
            for source, destination in zip(edges["src"], edges["dst"], strict=True)
        )
    ):
        raise ValueError(f"{path}: PGIF edges are malformed")
    return {
        "payload": payload,
        "states": states,
        "state_data": state_data,
        "frozen": columns["frozen"]["data"],
        "step": step,
        "time": time,
    }


def _pgif_static_identity(document: Mapping[str, Any]) -> dict[str, Any]:
    payload = document["payload"]
    petra = payload["meta"]["petra"]
    columns = payload["nodes"]["columns"]
    return {
        "meta": {
            **payload["meta"],
            "petra": {
                key: value
                for key, value in petra.items()
                if key not in {"step", "time"}
            },
        },
        "node_count": payload["nodes"]["count"],
        "static_node_columns": {
            name: column
            for name, column in columns.items()
            if name not in {"state", "type"}
        },
        "state_dict_and_schema": {
            key: value for key, value in columns["state"].items() if key != "data"
        },
        "type_dict_and_schema": {
            key: value for key, value in columns["type"].items() if key != "data"
        },
        "edges": payload["edges"],
    }


def _verify_package(
    root: Path,
    scenario: str,
    seed: int,
    manifest: Mapping[str, Any],
) -> dict[str, object]:
    receipt_path = root / "receipts" / scenario / f"seed-{seed}.json"
    if not _receipt_valid(root, receipt_path):
        raise ValueError(f"invalid external receipt for {scenario}/{seed}")
    external = json.loads(receipt_path.read_text(encoding="utf-8"))
    expected_package = Path("runs") / scenario / f"seed-{seed}"
    if (
        external.get("scenario") != scenario
        or external.get("seed") != seed
        or Path(external["package"]) != expected_package
        or external.get("deck_sha256")
        != next(row for row in manifest["scenarios"] if row["name"] == scenario)[
            "deck_sha256"
        ]
        or sha256_file(root / external["log"]) != external.get("log_sha256")
    ):
        raise ValueError("external package receipt identity/input mismatch")
    package = root / external["package"]
    internal = json.loads((package / "receipt.json").read_text(encoding="utf-8"))
    internal_files = PACKAGE_FILES - {"receipt.json"}
    if (
        internal.get("schema") != "a9b-replica-receipt-v2"
        or internal.get("seed") != seed
        or set(internal.get("files_sha256", {})) != internal_files
    ):
        raise ValueError("internal package receipt identity/inventory mismatch")
    for name, digest in internal["files_sha256"].items():
        if sha256_file(package / name) != digest:
            raise ValueError(f"internal hash mismatch for {name}")
    metadata = json.loads((package / "metadata.json").read_text(encoding="utf-8"))
    scenario_row = next(row for row in manifest["scenarios"] if row["name"] == scenario)
    inputs = manifest["inputs"]
    contract = manifest["run_contract"]
    metadata_inputs = metadata.get("inputs", {})
    command = external.get("command")
    expected_argv = [
        str((root / scenario_row["deck"]).resolve()),
        str((root / inputs["path_report"]).resolve()),
        str((root / inputs["lineage_snapshot"]).resolve()),
        str(package.resolve()),
        "--seed",
        str(seed),
        "--horizon",
        str(contract["horizon_s"]),
        "--checkpoints",
        str(contract["checkpoints"]),
        "--bias-factor",
        str(contract["bias_factor"]),
        "--max-events",
        str(contract["max_events"]),
    ]
    if (
        metadata.get("seed") != seed
        or metadata.get("schema") != "a9b-replica-package-v2"
        or metadata.get("survey_tier_platform_test") is not True
        or metadata.get("not_production") is not True
        or metadata.get("temperature_k") != 298.0
        or metadata.get("horizon_s") != contract["horizon_s"]
        or metadata.get("fixed_checkpoint_count") != contract["checkpoints"]
        or metadata.get("bias_factor") != contract["bias_factor"]
        or metadata.get("max_events") != contract["max_events"]
        or metadata_inputs.get("deck_sha256") != scenario_row["deck_sha256"]
        or metadata_inputs.get("path_report_sha256") != inputs["path_report_sha256"]
        or metadata_inputs.get("lineage_snapshot_sha256")
        != inputs["lineage_snapshot_sha256"]
        or metadata_inputs.get("runner_executable_sha256")
        != external.get("runner_sha256")
        or external.get("runner_sha256") != manifest.get("runner_sha256")
        or not isinstance(command, list)
        or len(command) != 18
        or command[:3] != ["nice", "-n", "10"]
        or command[3] != metadata_inputs.get("runner_executable")
        or metadata.get("argv") != expected_argv
        or command[4:] != expected_argv
    ):
        raise ValueError("metadata/run command/input binding mismatch")
    if metadata.get("census") != {
        "whole_finite_deck_original_centers": 320,
        "reachable_original_centers": 260,
        "topology_no_go_original_centers": 60,
        "whole_deck_estimand_includes_topology_no_go": True,
    }:
        raise ValueError("metadata census contract mismatch")
    mapping = metadata["reaction_mapping"]
    deck_reactions = [
        row["name"]
        for row in tomllib.loads(
            (root / scenario_row["deck"]).read_text(encoding="utf-8")
        )["reactions"]
    ]
    if len(mapping) != len(deck_reactions) or any(
        row.get("id") != index
        or row.get("name") != deck_reactions[index]
        or _finite(row.get("bias_factor"), "reaction bias factor") <= 0
        or not _close(
            float(row["bias_factor"]),
            contract["bias_factor"]
            if "hydrolysis" in deck_reactions[index]
            or deck_reactions[index].startswith("desorb-")
            else 1.0,
        )
        for index, row in enumerate(mapping)
    ):
        raise ValueError("reaction mapping/bias does not match bound scenario deck")
    reactions = [row["name"] for row in mapping]
    factors = [float(row["bias_factor"]) for row in mapping]

    likelihood_rows = read_jsonl(package / "likelihood.jsonl")
    if (
        not likelihood_rows
        or likelihood_rows[0].get("schema") != "a9b-likelihood-stream-v2"
        or likelihood_rows[0].get("seed") != seed
    ):
        raise ValueError("likelihood stream header mismatch")
    segments = likelihood_rows[1:]
    if len(segments) != metadata.get("likelihood_segment_count"):
        raise ValueError("likelihood segment count mismatch")
    cumulative = 0.0
    prior_end = 0.0
    for index, row in enumerate(segments):
        if row.get("segment") != index:
            raise ValueError("likelihood segment indices are not contiguous")
        start = _finite(row["start_time_s"], "segment start")
        end = _finite(row["end_time_s"], "segment end")
        physical = _finite(row["physical_total_rate_s-1"], "physical hazard")
        biased = _finite(row["biased_total_rate_s-1"], "biased hazard")
        if not _close(start, prior_end) or end < start or physical < 0 or biased < 0:
            raise ValueError("likelihood segment time/hazard contract failed")
        fired = row["fired_reaction_id"]
        event_term = 0.0
        if fired is not None:
            if physical <= 0 or biased <= 0:
                raise ValueError(
                    "fired likelihood segments require positive physical and biased hazards"
                )
            if type(fired) is not int or not 0 <= fired < len(reactions):
                raise ValueError("invalid fired reaction in likelihood segment")
            factor = _finite(row["fired_bias_factor"], "fired bias factor")
            if factor <= 0 or not _close(factor, factors[fired]):
                raise ValueError("likelihood fired bias factor mismatch")
            event_term = math.log(factor)
        elif row["fired_bias_factor"] is not None:
            raise ValueError("deadline segment cannot have a fired bias factor")
        increment = (biased - physical) * (end - start) - event_term
        if not _close(increment, row["log_likelihood_increment"]):
            raise ValueError("likelihood increment recomputation mismatch")
        cumulative += increment
        if not _close(cumulative, row["cumulative_log_likelihood"]):
            raise ValueError("likelihood cumulative recomputation mismatch")
        prior_end = end
    if not _close(cumulative, metadata["final_log_likelihood"]):
        raise ValueError("metadata final likelihood mismatch")

    initial_pgif = _pgif_document(package / "initial.pgif.json")
    final_pgif = _pgif_document(package / "final.pgif.json")
    state_names = initial_pgif["states"]
    initial = initial_pgif["state_data"]
    final_states = final_pgif["state_data"]
    final_step = final_pgif["step"]
    final_time = final_pgif["time"]
    if (
        initial_pgif["step"] != 0
        or initial_pgif["time"] != 0.0
        or _pgif_static_identity(initial_pgif) != _pgif_static_identity(final_pgif)
    ):
        raise ValueError("initial/final PGIF canonical identity mismatch")
    source_pgif = _pgif_document(root / inputs["lineage_snapshot"])
    if initial_pgif["payload"] != source_pgif["payload"]:
        raise ValueError(
            "replica initial PGIF differs from full bound lineage snapshot"
        )
    path_report = json.loads((root / inputs["path_report"]).read_text())
    ion_ledger = {row["site"]: row["kind"] for row in path_report["targets"]}
    events_rows = read_jsonl(package / "events.jsonl")
    if (
        not events_rows
        or events_rows[0].get("schema") != "a9b-event-stream-v2"
        or events_rows[0].get("seed") != seed
    ):
        raise ValueError("event stream header mismatch")
    if (
        events_rows[0]["states"] != state_names
        or events_rows[0]["reactions"] != reactions
        or len(events_rows) - 1 != metadata.get("event_count")
    ):
        raise ValueError("event header tables/count mismatch package")
    occupants = events_rows[0]["state_occupants"]
    if len(occupants) != len(state_names):
        raise ValueError("event state occupant table is malformed")
    current = list(initial)
    releases = {"Al": 0, "Si": 0}
    prior_time = 0.0
    event_segment_ids = set()
    checkpoints = read_jsonl(package / "checkpoints.jsonl")
    fixed_checkpoint_count = metadata["fixed_checkpoint_count"]
    horizon = _finite(metadata["horizon_s"], "horizon")
    if (
        not checkpoints
        or checkpoints[0]["checkpoint"] != 0
        or len(checkpoints) != metadata.get("recorded_checkpoint_count")
        or (metadata["complete"] and len(checkpoints) != fixed_checkpoint_count)
    ):
        raise ValueError("checkpoint stream count/start contract failed")
    event_index = 0
    prior_checkpoint_time = -1.0
    fixed_kind_totals = {
        kind: sum(
            count
            for state, count in checkpoints[0]["state_counts"].items()
            if state.startswith(f"{kind}.")
        )
        for kind in {state.partition(".")[0] for state in state_names}
    }
    for checkpoint_index, checkpoint in enumerate(checkpoints):
        checkpoint_time = _finite(checkpoint["actual_time_s"], "checkpoint time")
        target_time = _finite(checkpoint["target_time_s"], "checkpoint target time")
        role = checkpoint.get("checkpoint_role")
        if role == "fixed_physical_time":
            expected_target = horizon * checkpoint_index / (fixed_checkpoint_count - 1)
            if checkpoint_index >= fixed_checkpoint_count or not _close(
                target_time, expected_target
            ):
                raise ValueError("fixed checkpoint target cadence mismatch")
        elif role == "terminal_censor_snapshot":
            if metadata["complete"] or checkpoint_index != len(checkpoints) - 1:
                raise ValueError("terminal censor checkpoint is misplaced")
        else:
            raise ValueError("checkpoint role is invalid")
        if (
            checkpoint.get("checkpoint") != checkpoint_index
            or not _close(checkpoint_time, target_time)
            or checkpoint_time <= prior_checkpoint_time
            or checkpoint.get("kind_totals") != fixed_kind_totals
        ):
            raise ValueError("checkpoint identity/time/kind-total contract failed")
        expected_checkpoint_log = 0.0
        matching_segment = checkpoint_time == 0.0
        for segment in segments:
            segment_end = _finite(segment["end_time_s"], "segment end")
            if segment_end <= checkpoint_time or _close(segment_end, checkpoint_time):
                expected_checkpoint_log = _finite(
                    segment["cumulative_log_likelihood"], "segment cumulative"
                )
            if _close(segment_end, checkpoint_time):
                matching_segment = True
        if not matching_segment or not _close(
            expected_checkpoint_log, checkpoint["cumulative_log_likelihood"]
        ):
            raise ValueError("checkpoint likelihood/timestamp cross-check failed")
        while event_index < len(events_rows) - 1:
            row = events_rows[event_index + 1]
            event_time = _finite(row["time_s"], "event time")
            if event_time > checkpoint_time:
                break
            if row.get("event") != event_index or row.get("step") != event_index + 1:
                raise ValueError("event/step indices are not contiguous")
            if event_time < prior_time:
                raise ValueError("event timestamps are nonmonotonic")
            reaction_id = row["reaction_id"]
            if reactions[reaction_id] != row["reaction"]:
                raise ValueError("event reaction mapping mismatch")
            segment_id = row["likelihood_segment"]
            if (
                segment_id in event_segment_ids
                or segments[segment_id]["fired_reaction_id"] != reaction_id
                or not _close(segments[segment_id]["end_time_s"], event_time)
            ):
                raise ValueError("event/likelihood timestamp mapping mismatch")
            event_segment_ids.add(segment_id)
            center_site = row["center_site"]
            if type(center_site) is not int or not 0 <= center_site < len(current):
                raise ValueError("event center site is invalid")
            center_old_state = current[center_site]
            for change in row["changes"]:
                site = change["site"]
                old = change["old_state_id"]
                new = change["new_state_id"]
                if (
                    current[site] != old
                    or state_names[old] != change["old_state"]
                    or state_names[new] != change["new_state"]
                ):
                    raise ValueError("event state replay mismatch")
                current[site] = new
            expected_release = None
            if row["reaction"] in {"desorb-al", "desorb-si"}:
                kind = row["reaction"].removeprefix("desorb-").capitalize()
                if occupants[center_old_state] != kind or center_site not in {
                    change["site"] for change in row["changes"]
                }:
                    raise ValueError("desorption center/occupant semantics mismatch")
                if ion_ledger.get(center_site) == kind:
                    expected_release = kind
                    del ion_ledger[center_site]
                    releases[kind] += 1
            if row["lattice_release_kind"] != expected_release:
                raise ValueError("event lineage release annotation mismatch")
            prior_time = event_time
            event_index += 1
        counts = {name: 0 for name in state_names}
        for state in current:
            counts[state_names[state]] += 1
        if (
            checkpoint.get("checkpoint") != checkpoint_index
            or checkpoint["state_counts"] != counts
            or checkpoint.get("step") != event_index
        ):
            raise ValueError(
                "checkpoint state distribution/step disagrees with event replay"
            )
        if checkpoint["original_lattice_releases"] != releases:
            raise ValueError("checkpoint lineage counts disagree with event replay")
        prior_checkpoint_time = checkpoint_time
    if event_index != len(events_rows) - 1 or current != final_states:
        raise ValueError("final PGIF disagrees with complete event replay")
    if final_step != len(events_rows) - 1 or not _close(
        final_time, metadata["final_time_s"]
    ):
        raise ValueError("final PGIF step/time disagrees with metadata")
    if (
        not _close(prior_end, final_time)
        or not _close(prior_checkpoint_time, final_time)
        or metadata.get("remaining_original_ion_ledger") != len(ion_ledger)
    ):
        raise ValueError("likelihood/checkpoint/lineage terminal cross-check failed")
    if metadata["complete"]:
        if metadata["stop_reason"] != "physical_horizon" or not _close(
            final_time, horizon
        ):
            raise ValueError("complete replica did not reach its physical horizon")
    elif metadata["stop_reason"] == "physical_horizon" or final_time > horizon:
        raise ValueError("censored replica stop reason/time is invalid")
    if releases != metadata["original_lattice_releases"]:
        raise ValueError("metadata release counts disagree with replay")
    fired_segments = {
        index
        for index, row in enumerate(segments)
        if row["fired_reaction_id"] is not None
    }
    if fired_segments != event_segment_ids:
        raise ValueError("event stream does not cover every fired likelihood segment")

    times = [_finite(row["actual_time_s"], "checkpoint time") for row in checkpoints]
    areas = [
        _finite(row["geometric_area_a2"], "checkpoint area") for row in checkpoints
    ]
    area_time = math.fsum(
        (right - left) * (area_left + area_right) / 2.0
        for left, right, area_left, area_right in zip(
            times, times[1:], areas, areas[1:]
        )
    )
    if area_time <= 0 or not math.isfinite(area_time):
        raise ValueError("integrated area-time must be positive and finite")
    denominator = AVOGADRO * area_time * A2_TO_M2
    raw_rates = {kind: releases[kind] / denominator for kind in ("Si", "Al")}
    reachable_denominator = denominator * 260 / 320
    reachable_rates = {
        kind: releases[kind] / reachable_denominator for kind in ("Si", "Al")
    }
    fixed_checkpoints = [
        row for row in checkpoints if row["checkpoint_role"] == "fixed_physical_time"
    ]
    checkpoint_states = [row["state_counts"] for row in fixed_checkpoints]
    checkpoint_log_weights = [
        _finite(row["cumulative_log_likelihood"], "fixed checkpoint log weight")
        for row in fixed_checkpoints
    ]
    checkpoint_times = [
        _finite(row["actual_time_s"], "fixed checkpoint time")
        for row in fixed_checkpoints
    ]
    return {
        "scenario": scenario,
        "seed": seed,
        "complete": metadata["complete"],
        "stop_reason": metadata["stop_reason"],
        "log_weight": cumulative,
        "release_counts": releases,
        "whole_deck_rates": raw_rates,
        "reachable_census_rates": reachable_rates,
        "area_time_a2_s": area_time,
        "checkpoint_times_s": checkpoint_times,
        "checkpoint_state_counts": checkpoint_states,
        "checkpoint_log_weights": checkpoint_log_weights,
        "stationarity": population_stationarity(checkpoint_states),
        "event_count": len(events_rows) - 1,
        "segment_count": len(segments),
    }


def _stable_seed(*parts: str) -> int:
    digest = hashlib.sha256("\0".join(parts).encode()).digest()
    return int.from_bytes(digest[:8], "big")


def derive_analysis(root: Path) -> dict[str, object]:
    root = root.resolve()
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    if (
        manifest.get("schema") != MANIFEST_SCHEMA
        or manifest.get("status") != "raw_complete"
    ):
        raise ValueError("campaign raw data is not complete")
    smoke = bool(manifest["smoke"])
    _validate_existing_campaign(root, manifest, smoke=smoke)
    identities = [(row["scenario"], row["seed"]) for row in manifest["required_runs"]]
    validate_identity_set(identities, smoke=smoke)
    replicas = [
        _verify_package(root, scenario, seed, manifest) for scenario, seed in identities
    ]
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for replica in replicas:
        grouped[replica["scenario"]].append(replica)
    scenario_results: dict[str, dict[str, object]] = {}
    per_replica = []
    contributions: dict[str, dict[int, tuple[float, float]]] = {}
    for scenario_spec in scenario_specs():
        rows = sorted(grouped[scenario_spec.name], key=lambda row: row["seed"])
        log_weights = [row["log_weight"] for row in rows]
        all_replicas_complete = all(row["complete"] for row in rows)
        checkpoint_count = min(len(row["checkpoint_state_counts"]) for row in rows)
        checkpoint_times = rows[0]["checkpoint_times_s"][:checkpoint_count]
        if (
            checkpoint_count == 0
            or len(checkpoint_times) != checkpoint_count
            or (
                all_replicas_complete
                and any(
                    len(row["checkpoint_state_counts"]) != checkpoint_count
                    for row in rows
                )
            )
            or any(
                len(row["checkpoint_times_s"]) < checkpoint_count
                or len(row["checkpoint_log_weights"]) < checkpoint_count
                or any(
                    not _close(left, right)
                    for left, right in zip(
                        row["checkpoint_times_s"][:checkpoint_count],
                        checkpoint_times,
                        strict=True,
                    )
                )
                for row in rows
            )
        ):
            raise ValueError("replicas disagree on common fixed checkpoint cadence")
        evolution_points = []
        for checkpoint in range(checkpoint_count):
            states = rows[0]["checkpoint_state_counts"][checkpoint]
            checkpoint_log_weights = [
                row["checkpoint_log_weights"][checkpoint] for row in rows
            ]
            if checkpoint == 0:
                initial_counts = dict(states)
                if any(
                    not _close(weight, 0.0) for weight in checkpoint_log_weights
                ) or any(
                    row["checkpoint_state_counts"][0] != initial_counts for row in rows
                ):
                    raise ValueError(
                        "time-zero checkpoint weights/states must reproduce initial census"
                    )
                evolution_points.append(
                    {
                        state: {
                            "value": float(count),
                            "log_mean": math.log(count) if count > 0 else None,
                            "status": "represented",
                        }
                        for state, count in initial_counts.items()
                    }
                )
                continue
            evolution_points.append(
                {
                    state: _state_population_point(
                        checkpoint_log_weights,
                        [
                            row["checkpoint_state_counts"][checkpoint][state]
                            for row in rows
                        ],
                    )
                    for state in states
                }
            )
        evolution_counts = [
            {state: point["value"] for state, point in points.items()}
            for points in evolution_points
        ]
        evolution_statuses = {
            point["status"]
            for points in evolution_points
            for point in points.values()
            if point["status"] != "represented"
        }
        evolution = [
            {
                "time_s": time_s,
                "state_counts": counts,
                "log_state_counts": {
                    state: point["log_mean"] for state, point in points.items()
                },
                "numerical_status": (
                    "represented"
                    if not evolution_statuses
                    else next(iter(evolution_statuses))
                    if len(evolution_statuses) == 1
                    else "censored_numerical_range"
                ),
            }
            for time_s, counts, points in zip(
                checkpoint_times, evolution_counts, evolution_points, strict=True
            )
        ]
        fixed_kind_totals: dict[str, float] = defaultdict(float)
        for state, count in rows[0]["checkpoint_state_counts"][0].items():
            fixed_kind_totals[state.partition(".")[0]] += count
        if evolution_statuses:
            stationarity = {
                "passed": False,
                "status": "censored_numerical_state_evolution",
                "numerical_statuses": sorted(evolution_statuses),
                "tolerance_fraction": POPULATION_TOLERANCE,
                "kinds": {},
            }
        else:
            stationarity = population_stationarity(
                evolution_counts, fixed_kind_totals=fixed_kind_totals
            )
        all_replica_stationarity_passed = all(
            row["stationarity"]["passed"] for row in rows
        )
        stationarity_passed = bool(
            stationarity["passed"] and all_replica_stationarity_passed
        )
        species = {}
        for kind in ("Si", "Al"):
            observations = [row["whole_deck_rates"][kind] for row in rows]
            estimate = estimate_species(
                log_weights,
                observations,
                seed=_stable_seed(scenario_spec.name, kind),
            )
            species[kind] = _apply_run_and_stationarity_gates(
                estimate,
                runs_complete=all_replicas_complete,
                stationarity_passed=stationarity_passed,
            )
        combined_observations = [
            row["whole_deck_rates"]["Si"] + row["whole_deck_rates"]["Al"]
            for row in rows
        ]
        combined = estimate_species(
            log_weights,
            combined_observations,
            seed=_stable_seed(scenario_spec.name, "combined"),
        )
        combined = _apply_run_and_stationarity_gates(
            combined,
            runs_complete=all_replicas_complete,
            stationarity_passed=stationarity_passed,
        )
        raw_contributions = {
            row["seed"]: (
                row["log_weight"],
                row["whole_deck_rates"]["Si"] + row["whole_deck_rates"]["Al"],
            )
            for row in rows
        }
        contributions[scenario_spec.name] = raw_contributions
        si_mean = species["Si"]["mean"]
        al_mean = species["Al"]["mean"]
        stoichiometry = (
            si_mean / al_mean
            if isinstance(si_mean, float) and isinstance(al_mean, float) and al_mean > 0
            else None
        )
        scenario_results[scenario_spec.name] = {
            "family": scenario_spec.family,
            "perturbation": scenario_spec.perturbation,
            "species": species,
            "combined_si_al": combined,
            "stoichiometry_si_to_al_from_ensemble_means": stoichiometry,
            "state_population_evolution_unnormalized_is": evolution,
            "stationarity": stationarity,
            "all_replica_stationarity_passed": all_replica_stationarity_passed,
            "all_replicas_complete": all_replicas_complete,
            "stop_reasons": sorted({row["stop_reason"] for row in rows}),
        }
        for row in rows:
            per_replica.append(
                {
                    "scenario": row["scenario"],
                    "seed": row["seed"],
                    "log_weight": row["log_weight"],
                    "original_lattice_release_counts": row["release_counts"],
                    "whole_finite_deck_rates_mol_m2_s": row["whole_deck_rates"],
                    "reachable_census_rates_mol_m2_s": row["reachable_census_rates"],
                    "area_time_a2_s": row["area_time_a2_s"],
                    "event_count": row["event_count"],
                    "likelihood_segment_count": row["segment_count"],
                    "stationarity": row["stationarity"],
                    "stop_reason": row["stop_reason"],
                }
            )

    responses: dict[str, dict[str, object]] = {}
    family_rows = []
    for family in FAMILIES:
        perturbation_results = {}
        for perturbation in PERTURBATIONS:
            name = f"{family}__{perturbation}"
            nominal_result = scenario_results["nominal"]["combined_si_al"]
            perturbed_result = scenario_results[name]["combined_si_al"]
            response = paired_response(
                contributions["nominal"],
                contributions[name],
                nominal_rankable=bool(nominal_result["rankable"]),
                perturbed_rankable=bool(perturbed_result["rankable"]),
                seed=_stable_seed(family, perturbation, "response"),
            )
            perturbation_results[perturbation] = response
        responses[family] = perturbation_results
        rankable = all(
            result["status"] == "rankable" and result["delta_log10_ci95"] is not None
            for result in perturbation_results.values()
        )
        score = (
            max(abs(result["delta_log10"]) for result in perturbation_results.values())
            if rankable
            else None
        )
        irrelevant = bool(
            rankable
            and all(
                interval_demonstrates_irrelevance(result["delta_log10_ci95"])
                for result in perturbation_results.values()
            )
        )
        family_rows.append(
            {
                "family": family,
                "rankable": rankable,
                "score_max_abs_delta_log10": score,
                "demonstrated_irrelevant": irrelevant,
                "responses": perturbation_results,
            }
        )
    rankable_rows = sorted(
        (row for row in family_rows if row["rankable"]),
        key=lambda row: (-row["score_max_abs_delta_log10"], row["family"]),
    )
    for index, row in enumerate(rankable_rows, start=1):
        row["ordinal_rank_among_rankable"] = index
    top3 = [row["family"] for row in rankable_rows[:TOP_K]]
    irrelevant = [
        row["family"] for row in family_rows if row["demonstrated_irrelevant"]
    ]
    nominal_species = scenario_results["nominal"]["species"]
    lab_species_rates = {
        kind: estimate["mean"] if estimate["rankable"] else None
        for kind, estimate in nominal_species.items()
    }
    raw_complete = all(row["complete"] for row in replicas)
    stationary = all(
        result["stationarity"]["passed"] and result["all_replica_stationarity_passed"]
        for result in scenario_results.values()
    )
    if smoke:
        scientific_acceptance = "not_evaluated_smoke"
        ordinal_ranking_status = "not_evaluated_smoke"
    elif raw_complete and stationary:
        if not rankable_rows:
            scientific_acceptance = "accepted_typed_all_censored_no_go"
            ordinal_ranking_status = "no_rankable_families_no_go"
        elif len(rankable_rows) == len(FAMILIES):
            scientific_acceptance = "accepted_typed_outcomes"
            ordinal_ranking_status = "complete"
        else:
            scientific_acceptance = "accepted_typed_partial_ranking"
            ordinal_ranking_status = "partial_rankable_families"
    else:
        scientific_acceptance = "rejected_incomplete_or_nonstationary"
        ordinal_ranking_status = "not_scientifically_accepted"
    return _pending_analysis(
        {
            "schema": ANALYSIS_SCHEMA,
            "survey_tier_platform_test": True,
            "not_production": True,
            "artifact_verified": False,
            "scientific_acceptance": scientific_acceptance,
            "ordinal_ranking_complete": len(rankable_rows) == len(FAMILIES),
            "ordinal_ranking_status": ordinal_ranking_status,
            "simulated_scope": {
                "temperature_k": 298.0,
                "ph": 4,
                "ph3_ph5": "laboratory_comparison_bounds_only_not_simulated",
            },
            "estimator_contract": {
                "name": "unnormalized_importance_sampling",
                "aggregation": "log_domain",
                "self_normalized": False,
                "weight_ess_fraction_min": WEIGHT_ESS_FRACTION_MIN,
                "positive_species_contribution_ess_min": CONTRIBUTION_ESS_MIN,
                "positive_species_contribution_ess_role": CONTRIBUTION_ESS_ROLE,
                "zero_biased_observation_poisson_bound": "unavailable",
            },
            "census": manifest["estimand"],
            "run_count": len(replicas),
            "expected_run_count": 29 if smoke else 232,
            "per_replica": sorted(
                per_replica, key=lambda row: (row["scenario"], row["seed"])
            ),
            "scenarios": scenario_results,
            "family_ranking": family_rows,
            "top_k_preregistered": TOP_K,
            "top3_among_rankable": top3,
            "families_demonstrated_irrelevant": irrelevant,
            "irrelevance_statement": (
                "none demonstrated irrelevant"
                if not irrelevant
                else "only families whose complete response intervals lie inside the documented equivalence band"
            ),
            "equivalence_band_log10": list(EQUIVALENCE_BAND_LOG10),
            "laboratory_comparison": lab_comparison_scope(lab_species_rates),
        }
    )


def analyze_campaign(root: Path) -> dict[str, object]:
    root = root.resolve()
    remove_file_durable(root / "verification.json")
    analysis = derive_analysis(root)
    write_json_atomic(root / "analysis.json", analysis)
    return analysis


def invoke_independent_verifier(root: Path) -> dict[str, object]:
    subprocess.run(
        [
            sys.executable,
            str(Path(__file__).with_name("verify_a9b_sensitivity_ranking.py")),
            str(root),
            "--write-verification",
        ],
        check=True,
        shell=False,
    )
    return json.loads((root / "verification.json").read_text(encoding="utf-8"))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    validate = subparsers.add_parser("validate")
    validate.add_argument("output_root", type=Path)
    validate.add_argument("--smoke", action="store_true")
    run = subparsers.add_parser("run")
    run.add_argument("output_root", type=Path)
    run.add_argument("--runner", type=Path, required=True)
    run.add_argument("--workers", type=int, default=MAX_WORKERS)
    analyze = subparsers.add_parser("analyze")
    analyze.add_argument("output_root", type=Path)
    verify = subparsers.add_parser("verify")
    verify.add_argument("output_root", type=Path)
    smoke = subparsers.add_parser("smoke")
    smoke.add_argument("output_root", type=Path)
    smoke.add_argument("--runner", type=Path, required=True)
    smoke.add_argument("--workers", type=int, default=MAX_WORKERS)
    args = parser.parse_args(argv)

    if args.command == "validate":
        manifest = prepare_campaign(args.output_root, smoke=args.smoke)
        print(
            json.dumps(
                {"status": manifest["status"], "runs": manifest["required_run_count"]}
            )
        )
        return 0
    if args.command == "run":
        manifest = run_campaign(args.output_root, args.runner, workers=args.workers)
        print(
            json.dumps(
                {"status": manifest["status"], "runs": manifest["completed_run_count"]}
            )
        )
        return 0
    if args.command == "analyze":
        analysis = analyze_campaign(args.output_root)
        print(
            json.dumps(
                {
                    "artifact_verified": analysis["artifact_verified"],
                    "scientific_acceptance": analysis["scientific_acceptance"],
                    "ranking_complete": analysis["ordinal_ranking_complete"],
                }
            )
        )
        return 0
    if args.command == "verify":
        invoke_independent_verifier(args.output_root)
        return 0
    prepare_campaign(args.output_root, smoke=True)
    run_campaign(args.output_root, args.runner, workers=args.workers)
    analysis = analyze_campaign(args.output_root)
    verification = invoke_independent_verifier(args.output_root)
    print(
        json.dumps(
            {
                "artifact_verified": verification["artifact_verified"],
                "scientific_acceptance": verification["scientific_acceptance"],
                "runs": analysis["run_count"],
            }
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
