#!/usr/bin/env python3
"""Analyze the species-resolved scheduled muscovite mechanism run.

The Petra schedule runner reports populations at every segment boundary. Gallery
hops and trap escape conserve isotope inventories, so boundary-to-boundary loss
is exactly species release even though scheduled trajectory JSONL is not yet a
supported format. The apparent ages are synthetic observables, not geological
ages: they use the user-supplied irradiation J factor and the total 40K decay
constant solely to transform each step's 40Ar/39Ar release ratio.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
from dataclasses import asdict, dataclass
from itertools import pairwise
from pathlib import Path

import tomllib

K40_TOTAL_DECAY_PER_YEAR = 5.543e-10
YEARS_PER_MA = 1.0e6
GAS_CONSTANT_KCAL = 1.98720425864083e-3
PROXY_LABEL = "proxy ±5 kcal/mol; not computed kinetics"
MUSCOVITE_SITES_PER_CELL = 8
BOOTSTRAP_RESAMPLES = 2_000
MATERIAL_AR39_RELEASE_FRACTION = 0.001


@dataclass(frozen=True)
class SpectrumStep:
    segment: int
    temperature_k: float
    start_s: float
    end_s: float
    released_ar40: int
    released_ar39: int
    released_ar36: int
    cumulative_ar40: int
    cumulative_ar39: int
    cumulative_ar36: int
    apparent_age_ma: float | None
    ar36_ar40: float | None


@dataclass(frozen=True)
class AnalysisSummary:
    deck: str
    populations: str
    j_factor: float
    initial: dict[str, int]
    steps: tuple[SpectrumStep, ...]


@dataclass(frozen=True)
class SensitivityBracket:
    mechanism: str
    temperature_k: float
    barrier_kcal_mol: float
    offset_kcal_mol: float
    rate_multiplier: float
    provenance: str


@dataclass(frozen=True)
class DistributionBand:
    mean: float
    ci95: tuple[float, float]
    values: tuple[float, ...]


@dataclass(frozen=True)
class EnsembleSpectrumStep:
    segment: int
    temperature_k: float
    duration_s: float
    released_ar40: DistributionBand
    released_ar39: DistributionBand
    released_ar36: DistributionBand
    cumulative_ar40_fraction: DistributionBand
    cumulative_ar39_fraction: DistributionBand
    cumulative_ar36_fraction: DistributionBand
    apparent_age_ma: DistributionBand
    ar36_ar40: DistributionBand


@dataclass(frozen=True)
class EnsembleAnalysis:
    dims: tuple[int, int, int]
    sites: int
    replicas: int
    steps: tuple[EnsembleSpectrumStep, ...]


@dataclass(frozen=True)
class SizeRunReceipt:
    dims: tuple[int, int, int]
    sites: int
    replicas: int
    elapsed_seconds: float
    total_events: int
    replay_verified: bool


@dataclass(frozen=True)
class SizeStabilityComparison:
    previous_sites: int
    current_sites: int
    max_release_fraction_delta: float
    max_age_relative_delta: float | None
    max_age_defined_fraction_delta: float | None
    stable: bool


@dataclass(frozen=True)
class StabilityAssessment:
    stabilized_at_sites: int | None
    max_release_fraction_delta: float
    max_age_relative_delta: float | None
    max_age_defined_fraction_delta: float | None
    comparisons: tuple[SizeStabilityComparison, ...]


def _percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    return ordered[round((len(ordered) - 1) * fraction)]


def distribution_band(values: list[float], seed: int) -> DistributionBand:
    if not values:
        return DistributionBand(float("nan"), (float("nan"), float("nan")), ())
    data = tuple(float(value) for value in values)
    mean = statistics.fmean(data)
    if len(data) == 1:
        return DistributionBand(mean, (mean, mean), data)
    rng = random.Random(seed)
    bootstrap = [
        statistics.fmean(data[rng.randrange(len(data))] for _ in data)
        for _ in range(BOOTSTRAP_RESAMPLES)
    ]
    return DistributionBand(
        mean,
        (_percentile(bootstrap, 0.025), _percentile(bootstrap, 0.975)),
        data,
    )


def _schedule(deck_path: Path) -> list[dict[str, float]]:
    with deck_path.open("rb") as handle:
        deck = tomllib.load(handle)
    raw = deck.get("execution", {}).get("schedule", [])
    if not raw:
        raise ValueError(f"{deck_path}: no execution schedule")
    schedule: list[dict[str, float]] = []
    for index, item in enumerate(raw, start=1):
        temperature = float(item["temperature"])
        duration = float(item["duration"])
        if not math.isfinite(temperature) or temperature <= 0.0:
            raise ValueError(f"schedule segment {index}: invalid temperature")
        if not math.isfinite(duration) or duration <= 0.0:
            raise ValueError(f"schedule segment {index}: invalid duration")
        schedule.append({"temperature": temperature, "duration": duration})
    return schedule


def _isotope_for_column(name: str) -> str | None:
    state = name.rsplit(".", 1)[-1]
    if state == "Ar40":
        return "Ar40"
    if state in {"Ar39", "Ar39_trapped"}:
        return "Ar39"
    if state == "Ar36":
        return "Ar36"
    return None


def _population_rows(
    path: Path,
) -> tuple[dict[str, int], list[tuple[float, dict[str, int]]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "time" not in reader.fieldnames:
            raise ValueError(f"{path}: populations CSV has no time column")
        isotope_columns = {
            name: isotope
            for name in reader.fieldnames
            if (isotope := _isotope_for_column(name)) is not None
        }
        if set(isotope_columns.values()) != {"Ar40", "Ar39", "Ar36"}:
            raise ValueError(f"{path}: populations lack Ar40/Ar39/Ar36 states")
        rows: list[tuple[float, dict[str, int]]] = []
        previous_time = -math.inf
        for line_number, row in enumerate(reader, start=2):
            time_s = float(row["time"])
            if time_s < previous_time:
                raise ValueError(f"{path}:{line_number}: time moved backward")
            previous_time = time_s
            totals = {"Ar40": 0, "Ar39": 0, "Ar36": 0}
            for column, isotope in isotope_columns.items():
                totals[isotope] += int(row[column])
            rows.append((time_s, totals))
    if not rows:
        raise ValueError(f"{path}: no population rows")
    return rows[0][1], rows


def _boundary_row(
    rows: list[tuple[float, dict[str, int]]], target_s: float
) -> dict[str, int]:
    tolerance = max(1.0e-9, abs(target_s) * 1.0e-9)
    matches = [totals for time_s, totals in rows if abs(time_s - target_s) <= tolerance]
    if not matches:
        observed = ", ".join(f"{time_s:.6g}" for time_s, _ in rows[-8:])
        raise ValueError(
            f"no population report at schedule boundary {target_s:g}; tail={observed}"
        )
    return matches[-1]


def _ensemble_state_rows(
    deck_path: Path,
    observables_path: Path,
    expected_seeds: tuple[int, ...] | None = None,
) -> dict[int, list[tuple[float, dict[str, int]]]]:
    with deck_path.open("rb") as handle:
        deck = tomllib.load(handle)
    isotope_by_index: dict[int, str] = {}
    index = 0
    for kind in deck.get("structure", {}).get("kinds", []):
        for state in kind.get("states", []):
            isotope = _isotope_for_column(str(state["name"]))
            if isotope is not None:
                isotope_by_index[index] = isotope
            index += 1

    grouped: dict[tuple[int, float], dict[str, int]] = {}
    seen_indices: dict[tuple[int, float], set[int]] = {}
    replica_seeds: dict[int, int] = {}
    sample_steps: dict[tuple[int, float], int] = {}
    with observables_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["kind"] != "state_counts":
                continue
            replica = int(row["replica"])
            seed = int(row["seed"])
            step = int(row["step"])
            time_s = float(row["time"])
            sample = (replica, time_s)
            if replica in replica_seeds and replica_seeds[replica] != seed:
                raise ValueError(f"{observables_path}: replica {replica} changed seed")
            replica_seeds[replica] = seed
            if sample in sample_steps and sample_steps[sample] != step:
                raise ValueError(f"{observables_path}: sample {sample} mixed steps")
            sample_steps[sample] = step
            totals = grouped.setdefault(sample, {"Ar40": 0, "Ar39": 0, "Ar36": 0})
            state_index = int(row["index"])
            if state_index < 0 or state_index >= index:
                raise ValueError(
                    f"{observables_path}: state index {state_index} is outside deck range"
                )
            seen = seen_indices.setdefault(sample, set())
            if state_index in seen:
                raise ValueError(
                    f"{observables_path}: duplicate state index {state_index} in sample {sample}"
                )
            seen.add(state_index)
            try:
                value = float(row["value"])
            except ValueError as exc:
                raise ValueError(
                    f"{observables_path}: state count must be a non-negative integer"
                ) from exc
            if not math.isfinite(value) or value < 0 or not value.is_integer():
                raise ValueError(
                    f"{observables_path}: state count must be a non-negative integer"
                )
            isotope = isotope_by_index.get(state_index)
            if isotope is not None:
                totals[isotope] += int(value)
    expected_indices = set(range(index))
    for sample in grouped:
        if seen_indices.get(sample, set()) != expected_indices:
            raise ValueError(
                f"{observables_path}: incomplete state-count sample {sample}"
            )
    replicas: dict[int, list[tuple[float, dict[str, int]]]] = {}
    for (replica, time_s), totals in sorted(grouped.items()):
        replicas.setdefault(replica, []).append((time_s, totals))
    if not replicas:
        raise ValueError(f"{observables_path}: no ensemble state-count samples")
    if sorted(replicas) != list(range(len(replicas))):
        raise ValueError(
            f"{observables_path}: replica IDs must be contiguous from zero"
        )
    ordered_seeds = tuple(replica_seeds[index] for index in range(len(replicas)))
    if expected_seeds is not None and ordered_seeds != expected_seeds:
        raise ValueError(
            f"{observables_path}: replica seeds do not match ensemble receipt"
        )
    return replicas


def analyze_ensemble(
    deck_path: Path,
    observables_path: Path,
    dims: tuple[int, int, int],
    j_factor: float,
    expected_seeds: tuple[int, ...] | None = None,
) -> EnsembleAnalysis:
    if len(dims) != 3 or any(value <= 0 for value in dims):
        raise ValueError("lattice dimensions must be three positive integers")
    if not math.isfinite(j_factor) or j_factor <= 0.0:
        raise ValueError("J factor must be finite and positive")
    with deck_path.open("rb") as handle:
        deck = tomllib.load(handle)
    deck_dims = tuple(deck.get("structure", {}).get("lattice", {}).get("dims", ()))
    if deck_dims and deck_dims != dims:
        raise ValueError(f"deck dimensions {deck_dims} do not match requested {dims}")
    schedule = _schedule(deck_path)
    replicas = _ensemble_state_rows(deck_path, observables_path, expected_seeds)
    per_segment: list[dict[str, list[float]]] = [
        {
            "released_ar40": [],
            "released_ar39": [],
            "released_ar36": [],
            "cumulative_ar40_fraction": [],
            "cumulative_ar39_fraction": [],
            "cumulative_ar36_fraction": [],
            "apparent_age_ma": [],
            "ar36_ar40": [],
        }
        for _ in schedule
    ]
    for replica, rows in sorted(replicas.items()):
        initial = _boundary_row(rows, 0.0)
        previous = initial
        cumulative = {isotope: 0 for isotope in initial}
        boundary = 0.0
        for segment_index, segment in enumerate(schedule):
            boundary += segment["duration"]
            remaining = _boundary_row(rows, boundary)
            released = {
                isotope: previous[isotope] - remaining[isotope] for isotope in initial
            }
            if any(value < 0 for value in released.values()):
                raise ValueError(
                    f"replica {replica}: isotope inventory increased in segment {segment_index + 1}"
                )
            metrics = per_segment[segment_index]
            for isotope in ("Ar40", "Ar39", "Ar36"):
                cumulative[isotope] += released[isotope]
                metrics[f"released_{isotope.lower()}"].append(float(released[isotope]))
                if initial[isotope] > 0:
                    metrics[f"cumulative_{isotope.lower()}_fraction"].append(
                        cumulative[isotope] / initial[isotope]
                    )
            if released["Ar39"] > 0:
                metrics["apparent_age_ma"].append(
                    math.log1p(j_factor * released["Ar40"] / released["Ar39"])
                    / (K40_TOTAL_DECAY_PER_YEAR * YEARS_PER_MA)
                )
            if released["Ar40"] > 0:
                metrics["ar36_ar40"].append(released["Ar36"] / released["Ar40"])
            previous = remaining

    steps = []
    for segment_index, (segment, metrics) in enumerate(
        zip(schedule, per_segment, strict=True), start=1
    ):
        bands = {
            name: distribution_band(
                values, seed=0xE2B000 + segment_index * 100 + offset
            )
            for offset, (name, values) in enumerate(metrics.items())
        }
        steps.append(
            EnsembleSpectrumStep(
                segment=segment_index,
                temperature_k=segment["temperature"],
                duration_s=segment["duration"],
                **bands,
            )
        )
    return EnsembleAnalysis(
        dims=dims,
        sites=math.prod(dims) * MUSCOVITE_SITES_PER_CELL,
        replicas=len(replicas),
        steps=tuple(steps),
    )


def _band_cells(band: DistributionBand, replicas: int) -> list[float | int | str]:
    def finite_or_blank(value: float) -> float | str:
        return value if math.isfinite(value) else ""

    return [
        finite_or_blank(band.mean),
        finite_or_blank(band.ci95[0]),
        finite_or_blank(band.ci95[1]),
        len(band.values),
        replicas - len(band.values),
        ";".join(str(value) for value in band.values),
    ]


def write_ensemble_products(
    results: list[EnsembleAnalysis], receipts: list[SizeRunReceipt], out_dir: Path
) -> None:
    if not results or len(results) != len(receipts):
        raise ValueError("results and receipts must be non-empty and have equal length")
    out_dir.mkdir(parents=True, exist_ok=True)
    metrics = (
        "released_ar40",
        "released_ar39",
        "released_ar36",
        "cumulative_ar40_fraction",
        "cumulative_ar39_fraction",
        "cumulative_ar36_fraction",
        "apparent_age_ma",
        "ar36_ar40",
    )
    with (out_dir / "spectrum-bands.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "dims",
                "sites",
                "replicas",
                "segment",
                "temperature_k",
                "duration_s",
            ]
            + [
                field
                for metric in metrics
                for field in (
                    f"{metric}_mean",
                    f"{metric}_ci95_low",
                    f"{metric}_ci95_high",
                    f"{metric}_defined_replicas",
                    f"{metric}_missing_replicas",
                    f"{metric}_distribution",
                )
            ]
        )
        for result in results:
            dims = "x".join(str(value) for value in result.dims)
            for step in result.steps:
                row: list[float | int | str] = [
                    dims,
                    result.sites,
                    result.replicas,
                    step.segment,
                    step.temperature_k,
                    step.duration_s,
                ]
                for metric in metrics:
                    row.extend(_band_cells(getattr(step, metric), result.replicas))
                writer.writerow(row)

    with (out_dir / "size-scaling.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "dims",
                "sites",
                "replicas",
                "elapsed_seconds",
                "total_events",
                "events_per_second",
                "replay_verified",
            ]
        )
        for receipt in receipts:
            writer.writerow(
                [
                    "x".join(str(value) for value in receipt.dims),
                    receipt.sites,
                    receipt.replicas,
                    receipt.elapsed_seconds,
                    receipt.total_events,
                    receipt.total_events / receipt.elapsed_seconds,
                    str(receipt.replay_verified).lower(),
                ]
            )


def assess_stability(
    results: list[EnsembleAnalysis],
    release_fraction_tolerance: float = 0.05,
    age_relative_tolerance: float = 0.10,
    age_defined_fraction_tolerance: float = 0.05,
    material_ar39_release_fraction: float = MATERIAL_AR39_RELEASE_FRACTION,
) -> StabilityAssessment:
    if len(results) < 2:
        raise ValueError("stability assessment requires at least two lattice sizes")
    ordered = sorted(results, key=lambda result: result.sites)
    comparisons = []
    release_metrics = (
        "cumulative_ar40_fraction",
        "cumulative_ar39_fraction",
        "cumulative_ar36_fraction",
    )
    for previous, current in pairwise(ordered):
        if len(previous.steps) != len(current.steps):
            raise ValueError("all lattice sizes must have the same schedule")
        if any(
            (old.segment, old.temperature_k, old.duration_s)
            != (new.segment, new.temperature_k, new.duration_s)
            for old, new in zip(previous.steps, current.steps, strict=True)
        ):
            raise ValueError("all lattice sizes must have the same schedule")
        release_deltas = []
        age_deltas = []
        age_defined_deltas = []
        age_comparable = True
        previous_ar39_cumulative = 0.0
        current_ar39_cumulative = 0.0
        for old_step, new_step in zip(previous.steps, current.steps, strict=True):
            for metric in release_metrics:
                old_band = getattr(old_step, metric)
                new_band = getattr(new_step, metric)
                old = old_band.mean
                new = new_band.mean
                release_deltas.append(abs(new - old))
            old_ar39_segment = (
                old_step.cumulative_ar39_fraction.mean - previous_ar39_cumulative
            )
            new_ar39_segment = (
                new_step.cumulative_ar39_fraction.mean - current_ar39_cumulative
            )
            previous_ar39_cumulative = old_step.cumulative_ar39_fraction.mean
            current_ar39_cumulative = new_step.cumulative_ar39_fraction.mean
            if max(old_ar39_segment, new_ar39_segment) < material_ar39_release_fraction:
                continue
            old_age = old_step.apparent_age_ma.mean
            new_age = new_step.apparent_age_ma.mean
            age_defined_deltas.append(
                abs(
                    len(old_step.apparent_age_ma.values) / previous.replicas
                    - len(new_step.apparent_age_ma.values) / current.replicas
                )
            )
            if math.isfinite(old_age) and math.isfinite(new_age):
                age_deltas.append(abs(new_age - old_age) / max(abs(old_age), 1.0))
            elif math.isfinite(old_age) != math.isfinite(new_age):
                age_comparable = False
        max_release = max(release_deltas, default=float("inf"))
        max_age = max(age_deltas) if age_comparable and age_deltas else None
        max_age_defined = max(age_defined_deltas) if age_defined_deltas else None
        comparisons.append(
            SizeStabilityComparison(
                previous_sites=previous.sites,
                current_sites=current.sites,
                max_release_fraction_delta=max_release,
                max_age_relative_delta=max_age,
                max_age_defined_fraction_delta=max_age_defined,
                stable=max_release <= release_fraction_tolerance
                and max_age is not None
                and max_age <= age_relative_tolerance
                and max_age_defined is not None
                and max_age_defined <= age_defined_fraction_tolerance,
            )
        )
    stabilized_at = None
    for index, comparison in enumerate(comparisons):
        if all(item.stable for item in comparisons[index:]):
            stabilized_at = comparison.current_sites
            break
    return StabilityAssessment(
        stabilized_at_sites=stabilized_at,
        max_release_fraction_delta=max(
            item.max_release_fraction_delta for item in comparisons
        ),
        max_age_relative_delta=(
            max(
                item.max_age_relative_delta
                for item in comparisons
                if item.max_age_relative_delta is not None
            )
            if all(item.max_age_relative_delta is not None for item in comparisons)
            else None
        ),
        max_age_defined_fraction_delta=(
            max(
                item.max_age_defined_fraction_delta
                for item in comparisons
                if item.max_age_defined_fraction_delta is not None
            )
            if all(
                item.max_age_defined_fraction_delta is not None for item in comparisons
            )
            else None
        ),
        comparisons=tuple(comparisons),
    )


def analyze(
    deck_path: Path, populations_path: Path, j_factor: float
) -> AnalysisSummary:
    if not math.isfinite(j_factor) or j_factor <= 0.0:
        raise ValueError("J factor must be finite and positive")
    schedule = _schedule(deck_path)
    initial, rows = _population_rows(populations_path)
    previous_remaining = initial.copy()
    cumulative_time = 0.0
    steps: list[SpectrumStep] = []
    cumulative_release = {isotope: 0 for isotope in initial}
    for index, segment in enumerate(schedule, start=1):
        start_s = cumulative_time
        cumulative_time += segment["duration"]
        remaining = _boundary_row(rows, cumulative_time)
        released: dict[str, int] = {}
        for isotope in initial:
            delta = previous_remaining[isotope] - remaining[isotope]
            if delta < 0:
                raise ValueError(
                    f"{isotope} inventory increased in segment {index}: "
                    f"{previous_remaining[isotope]} -> {remaining[isotope]}"
                )
            released[isotope] = delta
            cumulative_release[isotope] += delta
        previous_remaining = remaining
        if released["Ar39"] > 0:
            ratio = released["Ar40"] / released["Ar39"]
            age_ma: float | None = math.log1p(j_factor * ratio) / (
                K40_TOTAL_DECAY_PER_YEAR * YEARS_PER_MA
            )
        else:
            age_ma = None
        contamination = (
            released["Ar36"] / released["Ar40"] if released["Ar40"] > 0 else None
        )
        steps.append(
            SpectrumStep(
                segment=index,
                temperature_k=segment["temperature"],
                start_s=start_s,
                end_s=cumulative_time,
                released_ar40=released["Ar40"],
                released_ar39=released["Ar39"],
                released_ar36=released["Ar36"],
                cumulative_ar40=cumulative_release["Ar40"],
                cumulative_ar39=cumulative_release["Ar39"],
                cumulative_ar36=cumulative_release["Ar36"],
                apparent_age_ma=age_ma,
                ar36_ar40=contamination,
            )
        )
    return AnalysisSummary(
        deck=str(deck_path),
        populations=str(populations_path),
        j_factor=j_factor,
        initial=initial,
        steps=tuple(steps),
    )


def sensitivity_brackets(
    mechanism: str,
    temperatures_k: list[float],
    nominal_kcal: float,
    delta_kcal: float,
) -> tuple[SensitivityBracket, ...]:
    if not mechanism.strip():
        raise ValueError("mechanism name must not be empty")
    if nominal_kcal <= 0.0 or delta_kcal <= 0.0 or delta_kcal >= nominal_kcal:
        raise ValueError("barrier bracket must satisfy 0 < delta < nominal")
    rows: list[SensitivityBracket] = []
    for temperature in temperatures_k:
        if not math.isfinite(temperature) or temperature <= 0.0:
            raise ValueError("temperatures must be finite and positive")
        for offset in (-delta_kcal, 0.0, delta_kcal):
            # k(E + offset) / k(E), with the common Eyring prefactor canceled.
            multiplier = math.exp(-offset / (GAS_CONSTANT_KCAL * temperature))
            rows.append(
                SensitivityBracket(
                    mechanism=mechanism,
                    temperature_k=temperature,
                    barrier_kcal_mol=nominal_kcal + offset,
                    offset_kcal_mol=offset,
                    rate_multiplier=multiplier,
                    provenance=PROXY_LABEL,
                )
            )
    return tuple(rows)


def summary_json(summary: AnalysisSummary) -> str:
    payload = asdict(summary)
    payload["interpretation"] = (
        "Synthetic step observables from a proxy-bracket mechanism; not a fitted calibration."
    )
    return json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"


def write_products(
    out_dir: Path, summary: AnalysisSummary, nominal_barrier: float
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "summary.json").write_text(summary_json(summary), encoding="utf-8")
    step_rows = [asdict(step) for step in summary.steps]
    with (out_dir / "spectrum.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(step_rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(step_rows)
    temperatures = sorted({step.temperature_k for step in summary.steps})
    sensitivity = tuple(
        row
        for mechanism, barrier in (
            ("extended-zone-hop", 40.0),
            ("delamination", nominal_barrier),
            ("octahedral-Ar39-escape", 64.0),
        )
        for row in sensitivity_brackets(mechanism, temperatures, barrier, 5.0)
    )
    with (out_dir / "sensitivity.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(asdict(sensitivity[0])),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(asdict(row) for row in sensitivity)
    (out_dir / "spectrum.svg").write_text(_spectrum_svg(summary), encoding="utf-8")


def _spectrum_svg(summary: AnalysisSummary) -> str:
    width, height = 900, 520
    left, top, plot_w, plot_h = 80, 55, 760, 360
    ages = [step.apparent_age_ma for step in summary.steps]
    finite = [age for age in ages if age is not None]
    ymax = max(finite, default=1.0) * 1.15 or 1.0
    bar_w = plot_w / max(1, len(summary.steps))
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfaf7"/>',
        '<text x="450" y="30" text-anchor="middle" font-family="sans-serif" font-size="21" font-weight="700">Synthetic 40Ar/39Ar step spectrum</text>',
        f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="white" stroke="#334155"/>',
    ]
    for index, (step, age) in enumerate(zip(summary.steps, ages, strict=True)):
        x = left + index * bar_w + 4
        if age is not None:
            h = age / ymax * plot_h
            y = top + plot_h - h
            parts.append(
                f'<rect x="{x:.2f}" y="{y:.2f}" width="{max(1.0, bar_w - 8):.2f}" height="{h:.2f}" fill="#b45309"/>'
            )
        parts.append(
            f'<text x="{x + (bar_w - 8) / 2:.2f}" y="{top + plot_h + 23}" text-anchor="middle" font-family="monospace" font-size="11">{step.temperature_k - 273.15:.0f} C</text>'
        )
    parts.extend(
        [
            f'<text x="{left + plot_w / 2}" y="{height - 35}" text-anchor="middle" font-family="sans-serif" font-size="14">incremental-heating step</text>',
            f'<text x="22" y="{top + plot_h / 2}" transform="rotate(-90 22 {top + plot_h / 2})" text-anchor="middle" font-family="sans-serif" font-size="14">apparent age (Ma)</text>',
            '<text x="450" y="500" text-anchor="middle" font-family="sans-serif" font-size="11" fill="#475569">Proxy mechanism and user-selected J factor; not a calibrated geological age</text>',
            "</svg>",
        ]
    )
    return "\n".join(parts) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deck", type=Path, required=True)
    parser.add_argument("--populations", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--j-factor", type=float, default=0.01)
    parser.add_argument("--nominal-proxy-barrier", type=float, default=58.0)
    args = parser.parse_args()
    summary = analyze(args.deck, args.populations, args.j_factor)
    write_products(args.out_dir, summary, args.nominal_proxy_barrier)
    print(summary_json(summary), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
