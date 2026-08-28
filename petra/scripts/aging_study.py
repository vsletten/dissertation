#!/usr/bin/env python3
"""Run and summarize the A5 Phase-1 finite-defect aging ensembles.

The campaign intentionally uses the reduced Phase-0 kaolinite surface-unit deck.
It measures a mechanism class; it is not a calibrated kaolinite prediction.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import random
import re
import statistics
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


DEFAULT_DENSITIES = (0.05, 0.15, 0.25, 0.50)
DEFAULT_REPLICAS = 64
DEFAULT_STEPS = 600
DEFAULT_AGED_STEP = 140
BOOTSTRAP_RESAMPLES = 2_000


@dataclass(frozen=True)
class Summary:
    mean: float
    ci95: tuple[float, float]


@dataclass(frozen=True)
class CurvePoint:
    step: int
    replicas: int
    time: Summary
    defects: Summary
    total_rate: Summary
    rate_log10_mean: Summary
    rate_log10_width: Summary
    exposure_age_mean: Summary
    exposure_age_p95: Summary
    geometric_area: Summary
    bet_proxy: Summary


@dataclass(frozen=True)
class DensityResult:
    density: float
    replicas: int
    initial_defects: Summary
    fresh_rate: Summary
    aged_rate: Summary
    drop_ratio: Summary
    geometric_normalized_drop: Summary
    bet_normalized_drop: Summary
    log10_drop: Summary
    curves: tuple[CurvePoint, ...]


@dataclass
class Sample:
    replica: int
    seed: int
    step: int
    time: float
    values: dict[str, list[float]]


def percentile(values: Iterable[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return 0.0
    return ordered[round((len(ordered) - 1) * fraction)]


def summarize(values: Iterable[float], seed: int) -> Summary:
    data = list(values)
    if not data:
        raise ValueError("cannot summarize an empty distribution")
    mean = statistics.fmean(data)
    if len(data) == 1:
        return Summary(mean, (mean, mean))
    rng = random.Random(seed)
    boot = []
    for _ in range(BOOTSTRAP_RESAMPLES):
        boot.append(statistics.fmean(data[rng.randrange(len(data))] for _ in data))
    return Summary(mean, (percentile(boot, 0.025), percentile(boot, 0.975)))


def density_slug(density: float) -> str:
    return f"{density:.3f}".replace(".", "p")


def render_deck(template: str, density: float) -> str:
    if not 0.0 < density < 1.0:
        raise ValueError("defect density must be between zero and one")
    rendered, names = re.subn(
        r'(?m)^name = "kaolinite-aging-smoke"$',
        f'name = "kaolinite-aging-density-{density_slug(density)}"',
        template,
    )
    rendered, probabilities = re.subn(
        r"(?m)^probability = 0\.25$", f"probability = {density:.6f}", rendered
    )
    if names != 1 or probabilities != 1:
        raise ValueError(
            "aging deck template must contain exactly one canonical name and probability"
        )
    return rendered


def load_samples(path: Path) -> dict[int, list[Sample]]:
    grouped: dict[tuple[int, int, int, float], dict[str, list[float]]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (
                int(row["replica"]),
                int(row["seed"]),
                int(row["step"]),
                float(row["time"]),
            )
            grouped.setdefault(key, {}).setdefault(row["kind"], []).append(
                float(row["value"])
            )
    replicas: dict[int, list[Sample]] = {}
    for (replica, seed, step, time), values in sorted(grouped.items()):
        replicas.setdefault(replica, []).append(
            Sample(replica, seed, step, time, values)
        )
    if not replicas:
        raise ValueError(f"no observable samples in {path}")
    return replicas


def sample_metric(sample: Sample, kind: str, index: int | None = None) -> float:
    values = sample.values.get(kind, [])
    if index is None:
        return sum(values)
    if len(values) <= index:
        raise ValueError(
            f"replica {sample.replica} step {sample.step} lacks {kind}[{index}]"
        )
    return values[index]


def per_sample_metrics(sample: Sample) -> dict[str, float]:
    spectrum = [value for value in sample.values.get("rate_spectra", []) if value > 0]
    rate_logs = [math.log10(value) for value in spectrum]
    ages = sample.values.get("exposure_age", [])
    surface = sample.values.get("surface_area", [])
    if len(surface) < 2:
        raise ValueError(
            f"replica {sample.replica} step {sample.step} lacks surface-area values"
        )
    return {
        "time": sample.time,
        "defects": sample_metric(sample, "state_counts", 1),
        "total_rate": sample_metric(sample, "event_rates"),
        "rate_log10_mean": statistics.fmean(rate_logs) if rate_logs else float("-inf"),
        "rate_log10_width": statistics.pstdev(rate_logs) if len(rate_logs) > 1 else 0.0,
        "exposure_age_mean": statistics.fmean(ages) if ages else 0.0,
        "exposure_age_p95": percentile(ages, 0.95),
        "geometric_area": surface[0],
        "bet_proxy": surface[1],
    }


def analyze_density(path: Path, density: float, aged_step: int) -> DensityResult:
    replicas = load_samples(path)
    fresh: list[Sample] = []
    aged: list[Sample] = []
    by_step: dict[int, list[Sample]] = {}
    for replica, samples in sorted(replicas.items()):
        del replica
        samples.sort(key=lambda sample: sample.step)
        fresh.append(samples[0])
        candidate = next((sample for sample in samples if sample.step >= aged_step), None)
        if candidate is None:
            raise ValueError(f"replica {samples[0].replica} never reaches step {aged_step}")
        defects = sample_metric(candidate, "state_counts", 1)
        if defects != 0.0:
            raise ValueError(
                f"replica {candidate.replica} aged sample still has defects ({defects:g})"
            )
        aged.append(candidate)
        for sample in samples:
            by_step.setdefault(sample.step, []).append(sample)

    def seeded(metric: str, values: Iterable[float], offset: int = 0) -> Summary:
        stable = sum(ord(char) for char in metric)
        return summarize(values, seed=0xA5_0000 + stable + offset)

    curves = []
    metric_names = (
        "time",
        "defects",
        "total_rate",
        "rate_log10_mean",
        "rate_log10_width",
        "exposure_age_mean",
        "exposure_age_p95",
        "geometric_area",
        "bet_proxy",
    )
    for step, samples in sorted(by_step.items()):
        rows = [per_sample_metrics(sample) for sample in samples]
        summaries = {
            name: seeded(name, (row[name] for row in rows), offset=step)
            for name in metric_names
        }
        curves.append(
            CurvePoint(step=step, replicas=len(samples), **summaries)
        )

    fresh_rates = [sample_metric(sample, "event_rates") for sample in fresh]
    aged_rates = [sample_metric(sample, "event_rates") for sample in aged]
    if any(value <= 0.0 for value in fresh_rates + aged_rates):
        raise ValueError("fresh and aged apparent rates must remain positive")
    ratios = [a / b for a, b in zip(fresh_rates, aged_rates, strict=True)]
    fresh_geometric = [sample_metric(sample, "surface_area", 0) for sample in fresh]
    aged_geometric = [sample_metric(sample, "surface_area", 0) for sample in aged]
    fresh_bet = [sample_metric(sample, "surface_area", 1) for sample in fresh]
    aged_bet = [sample_metric(sample, "surface_area", 1) for sample in aged]
    if any(
        value <= 0.0
        for value in fresh_geometric + aged_geometric + fresh_bet + aged_bet
    ):
        raise ValueError("fresh and aged surface measures must remain positive")
    geometric_ratios = [
        (fresh_rate / fresh_area) / (aged_rate / aged_area)
        for fresh_rate, fresh_area, aged_rate, aged_area in zip(
            fresh_rates, fresh_geometric, aged_rates, aged_geometric, strict=True
        )
    ]
    bet_ratios = [
        (fresh_rate / fresh_area) / (aged_rate / aged_area)
        for fresh_rate, fresh_area, aged_rate, aged_area in zip(
            fresh_rates, fresh_bet, aged_rates, aged_bet, strict=True
        )
    ]
    return DensityResult(
        density=density,
        replicas=len(replicas),
        initial_defects=seeded(
            "initial_defects",
            (sample_metric(sample, "state_counts", 1) for sample in fresh),
        ),
        fresh_rate=seeded("fresh_rate", fresh_rates),
        aged_rate=seeded("aged_rate", aged_rates),
        drop_ratio=seeded("drop_ratio", ratios),
        geometric_normalized_drop=seeded(
            "geometric_normalized_drop", geometric_ratios
        ),
        bet_normalized_drop=seeded("bet_normalized_drop", bet_ratios),
        log10_drop=seeded("log10_drop", (math.log10(value) for value in ratios)),
        curves=tuple(curves),
    )


def write_results(results: list[DensityResult], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    ratios_path = out_dir / "fresh-aged-ratios.csv"
    with ratios_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "defect_density",
                "replicas",
                "initial_defects_mean",
                "initial_defects_ci95_low",
                "initial_defects_ci95_high",
                "fresh_rate_mean",
                "fresh_rate_ci95_low",
                "fresh_rate_ci95_high",
                "aged_rate_mean",
                "aged_rate_ci95_low",
                "aged_rate_ci95_high",
                "drop_ratio_mean",
                "drop_ratio_ci95_low",
                "drop_ratio_ci95_high",
                "geometric_normalized_drop_mean",
                "geometric_normalized_drop_ci95_low",
                "geometric_normalized_drop_ci95_high",
                "bet_normalized_drop_mean",
                "bet_normalized_drop_ci95_low",
                "bet_normalized_drop_ci95_high",
                "log10_drop_mean",
                "log10_drop_ci95_low",
                "log10_drop_ci95_high",
            ]
        )
        for result in results:
            writer.writerow(
                [
                    result.density,
                    result.replicas,
                    result.initial_defects.mean,
                    *result.initial_defects.ci95,
                    result.fresh_rate.mean,
                    *result.fresh_rate.ci95,
                    result.aged_rate.mean,
                    *result.aged_rate.ci95,
                    result.drop_ratio.mean,
                    *result.drop_ratio.ci95,
                    result.geometric_normalized_drop.mean,
                    *result.geometric_normalized_drop.ci95,
                    result.bet_normalized_drop.mean,
                    *result.bet_normalized_drop.ci95,
                    result.log10_drop.mean,
                    *result.log10_drop.ci95,
                ]
            )

    curves_path = out_dir / "aging-curves.csv"
    with curves_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        fields = [
            "time",
            "defects",
            "total_rate",
            "rate_log10_mean",
            "rate_log10_width",
            "exposure_age_mean",
            "exposure_age_p95",
            "geometric_area",
            "bet_proxy",
        ]
        writer.writerow(
            ["defect_density", "step", "replicas"]
            + [item for field in fields for item in (f"{field}_mean", f"{field}_ci95_low", f"{field}_ci95_high")]
        )
        for result in results:
            for point in result.curves:
                row: list[float | int] = [result.density, point.step, point.replicas]
                for field in fields:
                    summary = getattr(point, field)
                    row.extend([summary.mean, *summary.ci95])
                writer.writerow(row)
    write_svg(results, out_dir / "aging-study.svg")


def _polyline(
    points: list[tuple[float, float]],
    bounds: tuple[float, float, float, float],
    frame: tuple[float, float, float, float],
) -> str:
    xmin, xmax, ymin, ymax = bounds
    x0, y0, width, height = frame
    if xmax == xmin:
        xmax = xmin + 1.0
    if ymax == ymin:
        ymax = ymin + 1.0
    return " ".join(
        f"{x0 + width * (x - xmin) / (xmax - xmin):.2f},"
        f"{y0 + height * (1.0 - (y - ymin) / (ymax - ymin)):.2f}"
        for x, y in points
    )


def write_svg(results: list[DensityResult], path: Path) -> None:
    if not results:
        raise ValueError("cannot plot an empty aging study")
    replica_count = results[0].replicas
    if any(result.replicas != replica_count for result in results):
        raise ValueError("all plotted densities must use the same replica count")
    colors = ("#38bdf8", "#34d399", "#facc15", "#fb7185", "#c084fc")
    series = []
    all_rate_points = []
    for result in results:
        points = [
            (point.exposure_age_mean.mean, math.log10(point.total_rate.mean))
            for point in result.curves
            if point.total_rate.mean > 0.0
        ]
        series.append(points)
        all_rate_points.extend(points)
    xmin = min(x for x, _ in all_rate_points)
    xmax = max(x for x, _ in all_rate_points)
    ymin = min(y for _, y in all_rate_points)
    ymax = max(y for _, y in all_rate_points)
    ratio_points = [(result.density, result.log10_drop.mean) for result in results]
    rxmin, rxmax = min(x for x, _ in ratio_points), max(x for x, _ in ratio_points)
    rymin, rymax = 0.0, max(y for _, y in ratio_points) * 1.1
    lines = []
    legend = []
    for index, (result, points) in enumerate(zip(results, series, strict=True)):
        color = colors[index % len(colors)]
        lines.append(
            f'<polyline fill="none" stroke="{color}" stroke-width="2.5" points="'
            + _polyline(points, (xmin, xmax, ymin, ymax), (90, 110, 820, 210))
            + '"/>'
        )
        legend.append(
            f'<text x="{100 + index * 190}" y="88" font-size="13" fill="{color}">p={result.density:.2f}</text>'
        )
    ratio_line = _polyline(ratio_points, (rxmin, rxmax, rymin, rymax), (90, 410, 820, 170))
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="660" viewBox="0 0 1000 660">
<rect width="100%" height="100%" fill="#0b1020"/>
<style>text{{font-family:ui-monospace,monospace;fill:#dbeafe}} .axis{{stroke:#64748b;stroke-width:1}}</style>
<text x="50" y="38" font-size="23" font-weight="bold">A5 Phase 1 — finite-defect kaolinite aging ensembles</text>
<text x="50" y="63" font-size="13" fill="#94a3b8">{replica_count} replicas per density; reduced mechanism-class deck, not calibrated kinetics</text>
{''.join(legend)}
<line class="axis" x1="90" y1="320" x2="910" y2="320"/><line class="axis" x1="90" y1="110" x2="90" y2="320"/>
{''.join(lines)}
<text x="96" y="128" font-size="14">log10 apparent propensity vs mean exposure age</text>
<text x="370" y="350" font-size="13">mean exposure age (simulation time)</text>
<line class="axis" x1="90" y1="580" x2="910" y2="580"/><line class="axis" x1="90" y1="410" x2="90" y2="580"/>
<polyline fill="none" stroke="#fb7185" stroke-width="3" points="{ratio_line}"/>
<text x="96" y="430" font-size="14">mean log10 fresh/aged apparent-rate drop</text>
<text x="390" y="615" font-size="13">initial defect probability</text>
<text x="50" y="646" font-size="12" fill="#94a3b8">Aged sample = step 140, after every replica's finite fast-site inventory is exhausted.</text>
</svg>
"""
    path.write_text(svg, encoding="utf-8")


def run_campaign(
    template_path: Path,
    petra_bin: Path,
    raw_root: Path,
    densities: tuple[float, ...],
    replicas: int,
    steps: int,
    base_seed: int,
) -> None:
    template = template_path.read_text(encoding="utf-8")
    decks = raw_root / "decks"
    logs = raw_root / "logs"
    decks.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "RAYON_NUM_THREADS"):
        env[name] = "16"
    for index, density in enumerate(densities):
        slug = density_slug(density)
        deck_path = decks / f"aging-{slug}.toml"
        deck_path.write_text(render_deck(template, density), encoding="utf-8")
        output = raw_root / f"density-{slug}"
        if output.exists():
            raise ValueError(f"refusing to overwrite existing campaign output: {output}")
        command = [
            "nice",
            "-n",
            "10",
            str(petra_bin),
            str(deck_path),
            "--steps",
            str(steps),
            "--seed",
            str(base_seed + index * 10_000),
            "--ensemble",
            str(replicas),
            "--out",
            str(output),
        ]
        with (logs / f"aging-{slug}.log").open("w", encoding="utf-8") as log:
            subprocess.run(
                command,
                check=True,
                timeout=600,
                cwd=template_path.parent.parent,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
            )


def parse_densities(text: str) -> tuple[float, ...]:
    densities = tuple(float(item) for item in text.split(","))
    if not densities or any(not 0.0 < value < 1.0 for value in densities):
        raise argparse.ArgumentTypeError("densities must be comma-separated values in (0,1)")
    return densities


def main() -> None:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="run bounded Petra ensembles")
    run.add_argument("template", type=Path)
    run.add_argument("petra_bin", type=Path)
    run.add_argument("raw_root", type=Path)
    run.add_argument("--densities", type=parse_densities, default=DEFAULT_DENSITIES)
    run.add_argument("--replicas", type=int, default=DEFAULT_REPLICAS)
    run.add_argument("--steps", type=int, default=DEFAULT_STEPS)
    run.add_argument("--base-seed", type=int, default=42_000)
    analyze = sub.add_parser("analyze", help="aggregate campaign observables")
    analyze.add_argument("raw_root", type=Path)
    analyze.add_argument("out_dir", type=Path)
    analyze.add_argument("--densities", type=parse_densities, default=DEFAULT_DENSITIES)
    analyze.add_argument("--aged-step", type=int, default=DEFAULT_AGED_STEP)
    args = parser.parse_args()
    if args.command == "run":
        if args.replicas < 2 or args.steps <= 0:
            parser.error("run requires at least two replicas and positive steps")
        run_campaign(
            args.template,
            args.petra_bin,
            args.raw_root,
            args.densities,
            args.replicas,
            args.steps,
            args.base_seed,
        )
    else:
        results = [
            analyze_density(
                args.raw_root / f"density-{density_slug(density)}" / "observables.csv",
                density,
                args.aged_step,
            )
            for density in args.densities
        ]
        write_results(results, args.out_dir)
        for result in results:
            print(
                f"density={result.density:.2f} n={result.replicas} "
                f"drop={result.drop_ratio.mean:.3f} "
                f"ci95=[{result.drop_ratio.ci95[0]:.3f},{result.drop_ratio.ci95[1]:.3f}]"
            )


if __name__ == "__main__":
    main()
