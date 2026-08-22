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
import tomllib
from dataclasses import asdict, dataclass
from pathlib import Path

K40_TOTAL_DECAY_PER_YEAR = 5.543e-10
YEARS_PER_MA = 1.0e6
GAS_CONSTANT_KCAL = 1.98720425864083e-3
PROXY_LABEL = "proxy ±5 kcal/mol; not computed kinetics"


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
