#!/usr/bin/env python3
"""Post-process Petra muscovite trajectories into release and D/a^2 products.

Each --run is LABEL:TEMPERATURE_C:RUN_DIR, where RUN_DIR contains the
`populations.csv` and `events.jsonl` emitted by `petra-cli --viz`.

The apparent D/a^2 calculation deliberately asks the same question as the
Sletten & Onstott (1998) treatment: what diffusion coefficient would the
Hames-Bowring cylindrical muscovite geometry infer from this non-Fickian
release curve? We invert Crank's infinite-cylinder cumulative-loss series
and difference the resulting dimensionless time at 2%-release crossings.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import math
import re
import statistics
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class ReleasePoint:
    fraction: float
    released: int
    time_s: float
    sqrt_time_sqrt_s: float
    cylinder_tau: float
    apparent_da2_per_s: float


@dataclass(frozen=True)
class GateMetrics:
    initial_ar: int
    released_ar: int
    final_fraction: float
    early_peak_fraction: float
    early_peak_da2_per_s: float
    first_da2_per_s: float
    tail_median_da2_per_s: float
    rise_ratio: float
    fall_ratio: float
    enough_release: bool
    early_rise: bool
    later_fall: bool
    pass_all: bool


@dataclass(frozen=True)
class RunResult:
    label: str
    temperature_c: float
    run_dir: str
    reaction_counts: dict[str, int]
    points: tuple[ReleasePoint, ...]
    gate: GateMetrics


# First 64 positive roots of J0(x), sufficient once tau >= 1e-3.
# The short-time heat-content expansion handles smaller tau without needing
# hundreds of roots. Values are fixed mathematical constants, not model data.
CYLINDER_J0_ZEROS = (
    2.40482555769577,
    5.52007811028631,
    8.65372791291101,
    11.7915344390143,
    14.9309177084878,
    18.0710639679109,
    21.2116366298793,
    24.3524715307493,
    27.4934791320403,
    30.634606468432,
    33.7758202135736,
    36.917098353664,
    40.0584257646282,
    43.1997917131767,
    46.3411883716618,
    49.4826098973978,
    52.624051841115,
    55.76551075502,
    58.9069839260809,
    62.0484691902272,
    65.1899648002069,
    68.3314693298568,
    71.4729816035937,
    74.6145006437018,
    77.7560256303881,
    80.8975558711376,
    84.0390907769382,
    87.1806298436412,
    90.3221726372105,
    93.4637187819448,
    96.6052679509963,
    99.7468198586806,
    102.888374254195,
    106.029930916452,
    109.171489649805,
    112.313050280495,
    115.454612653667,
    118.596176630873,
    121.737742087951,
    124.879308913233,
    128.020877006008,
    131.162446275214,
    134.304016638305,
    137.445588020284,
    140.587160352854,
    143.72873357369,
    146.870307625797,
    150.011882456955,
    153.153458019228,
    156.295034268534,
    159.436611164263,
    162.578188668947,
    165.719766747955,
    168.861345369236,
    172.002924503078,
    175.144504121903,
    178.286084200074,
    181.427664713731,
    184.569245640639,
    187.710826960049,
    190.852408652582,
    193.993990700109,
    197.135573085661,
    200.277155793332,
)


def cylinder_fraction(tau: float) -> float:
    """Cumulative radial release F from an infinite cylinder.

    ``tau = D*t/a^2``. Below 1e-3 the cylindrical heat-content expansion
    ``F = 4*sqrt(tau/pi) - tau`` avoids slow convergence; its error at the
    switch is below 1e-5 in F. Above the switch, the J0-root series is used.
    """

    if tau <= 0.0:
        return 0.0
    if tau < 1.0e-3:
        return min(1.0, 4.0 * math.sqrt(tau / math.pi) - tau)
    retained = 4.0 * sum(
        math.exp(-(root * root) * tau) / (root * root) for root in CYLINDER_J0_ZEROS
    )
    return max(0.0, min(1.0, 1.0 - retained))


def invert_cylinder_fraction(fraction: float) -> float:
    """Numerically invert :func:`cylinder_fraction` for 0 <= F < 1."""

    if fraction <= 0.0:
        return 0.0
    if fraction >= 1.0:
        raise ValueError("cylinder inversion requires fraction < 1")
    low, high = 0.0, 1.0
    while cylinder_fraction(high) < fraction:
        high *= 2.0
    for _ in range(80):
        middle = (low + high) / 2.0
        if cylinder_fraction(middle) < fraction:
            low = middle
        else:
            high = middle
    return (low + high) / 2.0


def _initial_argon(populations: Path) -> int:
    with populations.open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle)
        try:
            first = next(rows)
        except StopIteration as exc:
            raise ValueError(f"{populations}: no population rows") from exc
    columns = [name for name in first if name and name.endswith(".Ar40")]
    if not columns:
        raise ValueError(f"{populations}: no .Ar40 state columns")
    return sum(int(first[name]) for name in columns)


def _release_events(events: Path) -> tuple[list[float], dict[str, int]]:
    with events.open(encoding="utf-8") as handle:
        try:
            header = json.loads(next(handle))
        except StopIteration as exc:
            raise ValueError(f"{events}: empty event log") from exc
        reactions = header.get("reactions")
        if not isinstance(reactions, list):
            raise ValueError(f"{events}: header has no reaction table")
        release_ids = {
            index
            for index, name in enumerate(reactions)
            if isinstance(name, str) and name.startswith("release_")
        }
        if not release_ids:
            raise ValueError(f"{events}: no release_* reactions")
        counts = [0] * len(reactions)
        release_times: list[float] = []
        previous_time = -math.inf
        for line_number, line in enumerate(handle, start=2):
            row = json.loads(line)
            if not isinstance(row, list) or len(row) != 4:
                raise ValueError(f"{events}:{line_number}: invalid event row")
            _, time_s, reaction_id, _ = row
            if not isinstance(reaction_id, int) or not 0 <= reaction_id < len(
                reactions
            ):
                raise ValueError(f"{events}:{line_number}: invalid reaction id")
            time_s = float(time_s)
            if time_s < previous_time:
                raise ValueError(f"{events}:{line_number}: event time moved backward")
            previous_time = time_s
            counts[reaction_id] += 1
            if reaction_id in release_ids:
                release_times.append(time_s)
    return release_times, dict(zip(reactions, counts, strict=True))


def quantile_points(
    release_times: list[float], initial_ar: int, increment: float = 0.02
) -> tuple[ReleasePoint, ...]:
    if initial_ar <= 0:
        raise ValueError("initial Ar population must be positive")
    if len(release_times) > initial_ar:
        raise ValueError("release events exceed the initial Ar population")
    points: list[ReleasePoint] = []
    previous_time = 0.0
    previous_tau = 0.0
    fraction = increment
    final_fraction = len(release_times) / initial_ar
    while fraction <= final_fraction + 1.0e-12 and fraction < 1.0:
        released = math.ceil(fraction * initial_ar)
        time_s = release_times[released - 1]
        tau = invert_cylinder_fraction(fraction)
        elapsed = time_s - previous_time
        if elapsed <= 0.0:
            raise ValueError("release crossing times must increase")
        apparent = (tau - previous_tau) / elapsed
        points.append(
            ReleasePoint(
                fraction=fraction,
                released=released,
                time_s=time_s,
                sqrt_time_sqrt_s=math.sqrt(time_s),
                cylinder_tau=tau,
                apparent_da2_per_s=apparent,
            )
        )
        previous_time = time_s
        previous_tau = tau
        fraction = round(fraction + increment, 12)
    return tuple(points)


def evaluate_gate(
    points: tuple[ReleasePoint, ...], initial_ar: int, released_ar: int
) -> GateMetrics:
    final_fraction = released_ar / initial_ar
    enough_release = initial_ar >= 100 and final_fraction >= 0.50
    if len(points) < 6:
        # Marginal runs (low initial Ar or short trajectories) may not reach six
        # 2%-release crossings. Surface a clear gate failure instead of aborting
        # the whole analysis, so callers/CLI can still emit outputs.
        return GateMetrics(
            initial_ar=initial_ar,
            released_ar=released_ar,
            final_fraction=final_fraction,
            early_peak_fraction=math.nan,
            early_peak_da2_per_s=math.nan,
            first_da2_per_s=math.nan,
            tail_median_da2_per_s=math.nan,
            rise_ratio=math.nan,
            fall_ratio=math.nan,
            enough_release=enough_release,
            early_rise=False,
            later_fall=False,
            pass_all=False,
        )
    early = [point for point in points if point.fraction <= 0.12 + 1.0e-12]
    peak = max(early, key=lambda point: point.apparent_da2_per_s)
    tail_values = [point.apparent_da2_per_s for point in points[-5:]]
    tail_median = statistics.median(tail_values)
    first = points[0].apparent_da2_per_s
    rise_ratio = peak.apparent_da2_per_s / first
    fall_ratio = peak.apparent_da2_per_s / tail_median
    early_rise = peak.fraction > points[0].fraction and rise_ratio >= 1.5
    later_fall = fall_ratio >= 10.0
    return GateMetrics(
        initial_ar=initial_ar,
        released_ar=released_ar,
        final_fraction=final_fraction,
        early_peak_fraction=peak.fraction,
        early_peak_da2_per_s=peak.apparent_da2_per_s,
        first_da2_per_s=first,
        tail_median_da2_per_s=tail_median,
        rise_ratio=rise_ratio,
        fall_ratio=fall_ratio,
        enough_release=enough_release,
        early_rise=early_rise,
        later_fall=later_fall,
        pass_all=enough_release and early_rise and later_fall,
    )


def analyze_run(label: str, temperature_c: float, run_dir: Path) -> RunResult:
    initial_ar = _initial_argon(run_dir / "populations.csv")
    release_times, reaction_counts = _release_events(run_dir / "events.jsonl")
    points = quantile_points(release_times, initial_ar)
    gate = evaluate_gate(points, initial_ar, len(release_times))
    return RunResult(
        label=label,
        temperature_c=temperature_c,
        run_dir=str(run_dir),
        reaction_counts=reaction_counts,
        points=points,
        gate=gate,
    )


def write_csv(path: Path, result: RunResult) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(asdict(result.points[0])),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(asdict(point) for point in result.points)


def _ticks(low: float, high: float, count: int = 5) -> list[float]:
    if count <= 1 or high <= low:
        return [low]
    return [low + index * (high - low) / (count - 1) for index in range(count)]


def _span(low: float, high: float) -> float:
    """Non-zero plot span, so a collapsed range renders as a degenerate plot.

    When every x (or y) value in a panel is identical, ``high - low`` is zero
    and the ``sx``/``sy`` scale functions would divide by zero. Return a unit
    span instead so all points map to the same coordinate rather than raising.
    """

    span = high - low
    return span if span > 0.0 else 1.0


def _svg_panels(results: tuple[RunResult, ...], metric: str) -> str:
    width = 1100
    panel_height = 340
    height = 80 + panel_height * len(results)
    left, right, top, bottom = 100, 55, 45, 65
    plot_width = width - left - right
    plot_height = panel_height - top - bottom
    title = (
        "Cumulative 40Ar release vs square-root time"
        if metric == "release"
        else "Apparent D/a^2 from Hames-Bowring cylinder inversion"
    )
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fbfaf7"/>',
        f'<text x="{width / 2}" y="30" text-anchor="middle" font-family="sans-serif" font-size="22" font-weight="700">{html.escape(title)}</text>',
    ]
    for panel, result in enumerate(results):
        origin_y = 55 + panel * panel_height
        x0, y0 = left, origin_y + top
        if metric == "release":
            xs = [point.sqrt_time_sqrt_s for point in result.points]
            ys = [point.fraction for point in result.points]
            x_label, y_label = "sqrt(time / s)", "fraction released"
            x_low, x_high = 0.0, max(xs) * 1.03
            y_low, y_high = 0.0, max(ys) * 1.05
            y_precision = 2
        else:
            xs = [point.fraction for point in result.points]
            raw_y = [point.apparent_da2_per_s for point in result.points]
            ys = [math.log10(value) for value in raw_y]
            x_label, y_label = "fraction released", "log10(D/a^2 / s^-1)"
            x_low, x_high = 0.0, max(xs) * 1.03
            y_low = math.floor(min(ys))
            y_high = math.ceil(max(ys))
            y_precision = 1

        # Guard against a degenerate panel (all values identical) collapsing
        # the range to zero, which would divide by zero in sx/sy below.
        if x_high <= x_low:
            x_high = x_low + 1.0
        if y_high <= y_low:
            y_high = y_low + 1.0

        def display_y(value: float) -> str:
            return f"{value:.{y_precision}f}"

        def sx(value: float) -> float:
            return x0 + (value - x_low) * plot_width / _span(x_low, x_high)

        def sy(value: float) -> float:
            return (
                y0 + plot_height - (value - y_low) * plot_height / _span(y_low, y_high)
            )

        parts.extend(
            [
                f'<rect x="{x0}" y="{y0}" width="{plot_width}" height="{plot_height}" fill="#ffffff" stroke="#4b5563"/>',
                f'<text x="{x0}" y="{origin_y + 25}" font-family="sans-serif" font-size="18" font-weight="700">{html.escape(result.label)}</text>',
            ]
        )
        for value in _ticks(x_low, x_high):
            x = sx(value)
            parts.append(
                f'<line x1="{x:.2f}" y1="{y0}" x2="{x:.2f}" y2="{y0 + plot_height}" stroke="#e5e7eb"/>'
            )
            parts.append(
                f'<text x="{x:.2f}" y="{y0 + plot_height + 22}" text-anchor="middle" font-family="monospace" font-size="12">{value:.2g}</text>'
            )
        for value in _ticks(y_low, y_high):
            y = sy(value)
            parts.append(
                f'<line x1="{x0}" y1="{y:.2f}" x2="{x0 + plot_width}" y2="{y:.2f}" stroke="#e5e7eb"/>'
            )
            parts.append(
                f'<text x="{x0 - 12}" y="{y + 4:.2f}" text-anchor="end" font-family="monospace" font-size="12">{display_y(value)}</text>'
            )
        coordinates = " ".join(
            f"{sx(x):.2f},{sy(y):.2f}" for x, y in zip(xs, ys, strict=True)
        )
        parts.append(
            f'<polyline points="{coordinates}" fill="none" stroke="#b42318" stroke-width="3"/>'
        )
        parts.append(
            f'<text x="{x0 + plot_width / 2}" y="{y0 + plot_height + 48}" text-anchor="middle" font-family="sans-serif" font-size="14">{html.escape(x_label)}</text>'
        )
        parts.append(
            f'<text x="{x0 + 4}" y="{origin_y + 42}" font-family="sans-serif" font-size="12" fill="#374151">{html.escape(y_label)}</text>'
        )
        verdict = result.gate
        verdict_text = "gate PASS" if verdict.pass_all else "gate FAIL"
        verdict_fill = "#166534" if verdict.pass_all else "#b42318"
        if math.isfinite(verdict.rise_ratio) and math.isfinite(verdict.fall_ratio):
            verdict_detail = (
                f"rise x{verdict.rise_ratio:.1f}; fall x{verdict.fall_ratio:.0f}; "
            )
        else:
            verdict_detail = ""
        parts.append(
            f'<text x="{x0 + plot_width - 5}" y="{origin_y + 25}" text-anchor="end" font-family="monospace" font-size="13" fill="{verdict_fill}">{verdict_detail}{verdict_text}</text>'
        )
    parts.append("</svg>")
    return "\n".join(parts) + "\n"


def _slug(label: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", label.lower()).strip("-")
    if not slug:
        raise ValueError(f"label {label!r} has no filename-safe characters")
    return slug


def parse_run(value: str) -> tuple[str, float, Path]:
    try:
        label, temperature, directory = value.split(":", 2)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "--run must be LABEL:TEMPERATURE_C:RUN_DIR"
        ) from exc
    path = Path(directory).resolve()
    if not path.is_dir():
        raise argparse.ArgumentTypeError(f"run directory does not exist: {path}")
    return label, float(temperature), path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", type=parse_run, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    results = tuple(analyze_run(*spec) for spec in args.run)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for result in results:
        write_csv(args.out_dir / f"{_slug(result.label)}.csv", result)
    (args.out_dir / "release-curves.svg").write_text(
        _svg_panels(results, "release"), encoding="utf-8"
    )
    (args.out_dir / "apparent-da2.svg").write_text(
        _svg_panels(results, "da2"), encoding="utf-8"
    )
    summary = {
        "inversion": "Hames-Bowring infinite-cylinder cumulative-loss series; 2%-release increments",
        "runs": [
            {
                "label": result.label,
                "temperature_c": result.temperature_c,
                "reaction_counts": result.reaction_counts,
                "gate": asdict(result.gate),
            }
            for result in results
        ],
        "all_gates_pass": all(result.gate.pass_all for result in results),
    }
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    for result in results:
        gate = result.gate
        print(
            f"{result.label}: released {gate.released_ar}/{gate.initial_ar} "
            f"({gate.final_fraction:.1%}); D/a^2 rise x{gate.rise_ratio:.2f}, "
            f"fall x{gate.fall_ratio:.1f}; "
            f"{'PASS' if gate.pass_all else 'FAIL'}"
        )
    print(f"wrote {args.out_dir}")
    return 0 if summary["all_gates_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
