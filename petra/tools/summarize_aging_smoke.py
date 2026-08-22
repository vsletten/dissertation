#!/usr/bin/env python3
"""Condense Petra's long-form aging-smoke output and draw a dependency-free SVG."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path


def percentile(values: list[float], fraction: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    return ordered[round((len(ordered) - 1) * fraction)]


def load_samples(raw_dir: Path) -> list[dict[str, float]]:
    populations: dict[int, dict[str, float]] = {}
    with (raw_dir / "populations.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            step = int(row["step"])
            populations[step] = {
                "terrace": float(row["Kaolinite_site.terrace"]),
                "defect": float(row["Kaolinite_site.defect"]),
                "empty": float(row["Kaolinite_site.empty"]),
            }

    grouped: dict[tuple[int, float], dict[str, list[float]]] = {}
    with (raw_dir / "observables.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = (int(row["step"]), float(row["time"]))
            grouped.setdefault(key, {}).setdefault(row["kind"], []).append(
                float(row["value"])
            )

    samples = []
    for (step, time), values in sorted(grouped.items()):
        rates = values.get("event_rates", [])
        spectrum = values.get("rate_spectra", [])
        logs = [math.log10(rate) for rate in spectrum if rate > 0.0]
        surface = values.get("surface_area", [0.0, 0.0, 0.0])
        ages = values.get("exposure_age", [])
        sample = {
            "step": float(step),
            "time": time,
            **populations[step],
            "total_rate": sum(rates),
            "spectrum_mean": statistics.fmean(spectrum) if spectrum else 0.0,
            "spectrum_log10_std": statistics.pstdev(logs) if len(logs) > 1 else 0.0,
            "geometric_area": surface[0],
            "bet_proxy": surface[1],
            "exposed_sites": surface[2],
            "mean_exposure_age": statistics.fmean(ages) if ages else 0.0,
            "p95_exposure_age": percentile(ages, 0.95),
        }
        samples.append(sample)
    return samples


def write_csv(samples: list[dict[str, float]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(samples[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(samples)


def points(
    samples: list[dict[str, float]],
    key: str,
    x0: float,
    y0: float,
    width: float,
    height: float,
    *,
    logarithmic: bool = False,
    maximum: float | None = None,
) -> str:
    xs = [sample["time"] for sample in samples]
    ys = [sample[key] for sample in samples]
    if logarithmic:
        positive = [value for value in ys if value > 0.0]
        floor = min(positive) / 10.0
        ys = [math.log10(max(value, floor)) for value in ys]
    xmin, xmax = min(xs), max(xs)
    ymin = min(ys)
    ymax = maximum if maximum is not None else max(ys)
    if logarithmic and maximum is not None:
        ymax = math.log10(maximum)
    if ymax == ymin:
        ymax = ymin + 1.0
    return " ".join(
        f"{x0 + width * (x - xmin) / (xmax - xmin):.2f},"
        f"{y0 + height * (1.0 - (y - ymin) / (ymax - ymin)):.2f}"
        for x, y in zip(xs, ys, strict=True)
    )


def write_svg(samples: list[dict[str, float]], path: Path) -> None:
    width, height = 1040, 680
    rate_points = points(samples, "total_rate", 90, 100, 870, 210, logarithmic=True)
    defect_points = points(samples, "defect", 90, 400, 870, 190, maximum=576.0)
    terrace_points = points(samples, "terrace", 90, 400, 870, 190, maximum=576.0)
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
<rect width="100%" height="100%" fill="#0b1020"/>
<style>text{{font-family:ui-monospace,monospace;fill:#dbeafe}} .axis{{stroke:#64748b;stroke-width:1}} .grid{{stroke:#1e293b;stroke-width:1}}</style>
<text x="52" y="42" font-size="24" font-weight="bold">A5 Phase-0 — finite-defect surface-aging smoke (seed 42)</text>
<text x="52" y="68" font-size="14" fill="#94a3b8">Reduced kaolinite surface-unit deck; qualitative H1 observable gate, not calibrated kinetics</text>
<line class="axis" x1="90" y1="310" x2="960" y2="310"/><line class="axis" x1="90" y1="100" x2="90" y2="310"/>
<polyline fill="none" stroke="#fb7185" stroke-width="3" points="{rate_points}"/>
<text x="100" y="122" font-size="16">bulk propensity (log scale)</text>
<line class="axis" x1="90" y1="590" x2="960" y2="590"/><line class="axis" x1="90" y1="400" x2="90" y2="590"/>
<polyline fill="none" stroke="#facc15" stroke-width="3" points="{defect_points}"/>
<polyline fill="none" stroke="#38bdf8" stroke-width="3" points="{terrace_points}"/>
<text x="100" y="422" font-size="16"><tspan fill="#facc15">defect inventory</tspan><tspan> / </tspan><tspan fill="#38bdf8">terrace inventory</tspan></text>
<text x="430" y="635" font-size="15">simulation time (arbitrary smoke units)</text>
<text x="52" y="662" font-size="13" fill="#94a3b8">33 seeded fast sites exhaust by step 40; the remaining terrace channel is 1000× slower.</text>
</svg>
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(svg)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("raw_dir", type=Path)
    parser.add_argument("summary_csv", type=Path)
    parser.add_argument("plot_svg", type=Path)
    args = parser.parse_args()
    samples = load_samples(args.raw_dir)
    if not samples:
        raise SystemExit("no observable samples")
    write_csv(samples, args.summary_csv)
    write_svg(samples, args.plot_svg)
    first, aged = samples[0], samples[min(4, len(samples) - 1)]
    print(
        f"samples={len(samples)} initial_defects={first['defect']:.0f} "
        f"step40_defects={aged['defect']:.0f} "
        f"rate_ratio={aged['total_rate'] / first['total_rate']:.9g}"
    )


if __name__ == "__main__":
    main()
