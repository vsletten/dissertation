#!/usr/bin/env python3
"""Run and analyze the E4b isothermal reservoir-discriminator campaign.

This is a hypothesis screen, not a fit.  It combines the E4a surface-gated deck
with (1) explicit fast-access 36Ar reservoir rules, (2) E3a's Xe barrier as a
labelled screen, and (3) a forward dehydroxylation proxy: H2O activity scales
only the forward dehydroxylate_pair rate; there is no reverse dehydroxylation
rule.  The analysis applies the same infinite-cylinder cumulative-loss inversion
used by Sletten & Onstott (1998).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import random
import re
import statistics
import subprocess
import time
from collections import defaultdict
from pathlib import Path

import tomllib
from build_muscovite_full_deck import render_deck
from muscovite_release import invert_cylinder_fraction

SIZES = ((4, 4, 6), (8, 8, 6), (12, 12, 6))
TEMPERATURES_C = (500, 700)
# The engine's reservoir factor multiplies a forward rate.  H2O therefore enters
# as an explicitly named vacancy activity, a_vac=1 in vacuum and 1/P_bar for the
# ideal 2-kbar screen.  This is a directionality test, not a fugacity model.
CONDITIONS = {"vacuum": 1.0, "hydrothermal-2kbar-screen": 5.0e-4}
SQRT_TIME_MIN_HALF = (2, 4, 6, 8, 10, 15, 20, 30, 40, 50, 60)
REPLICAS = 8
BASE_SEED = 19_981
SITES_PER_CELL = 8
E3A_AR_LOCAL_KCAL_MOL = 64.095991
E3A_AR_EXTENDED_KCAL_MOL = 68.410690
E3A_XE_SCREEN_KCAL_MOL = 90.744700


def slug_dims(dims: tuple[int, int, int]) -> str:
    return "x".join(map(str, dims))


def schedule_durations() -> tuple[float, ...]:
    elapsed = [value * value * 60.0 for value in SQRT_TIME_MIN_HALF]
    return tuple(
        current - previous for previous, current in zip((0.0, *elapsed), elapsed)
    )


def _replace_rule_rate(text: str, prefix: str, barrier: float) -> str:
    pattern = re.compile(
        rf'(name = "{re.escape(prefix)}[^\"]*"\n.*?rate = \{{ eyring = \{{ dh = )[-+0-9.eE]+(, ds = 0\.0 \}} \}})',
        re.DOTALL,
    )
    text, count = pattern.subn(rf"\g<1>{barrier:.6f}\g<2>", text)
    if count != 4:
        raise ValueError(f"expected four {prefix} rules, found {count}")
    return text


def render_e4b_deck(
    dims: tuple[int, int, int], temperature_c: int, condition: str
) -> str:
    if (
        dims not in SIZES
        or temperature_c not in TEMPERATURES_C
        or condition not in CONDITIONS
    ):
        raise ValueError("unknown E4b deck coordinate")
    text = render_deck(dims)
    name = f"muscovite-e4b-{temperature_c}c-{condition}-{slug_dims(dims)}"
    text = re.sub(
        r'name = "muscovite-full(?:-[^\"]+)?-mechanism"',
        f'name = "{name}"',
        text,
        count=1,
    )
    text = re.sub(r'name = "muscovite-full-[^\"]+"', f'name = "{name}"', text, count=1)
    text = text.replace(
        'comment = "mechanism comparison; proxy brackets are not computed kinetics"',
        'comment = "E4b isothermal discriminator; Ar/Xe and H2O terms are labelled screens, not fits"',
        1,
    )

    species_anchor = '[[structure.species]]\nname = "Ar36"\nmass = 35.967545\n'
    text = text.replace(
        species_anchor,
        species_anchor
        + '[[structure.species]]\nname = "Xe"\nmass = 131.293\n'
        + '[[structure.species]]\nname = "H2O_vacancy"\n',
        1,
    )
    state_anchor = '[[structure.kinds.states]]\nname = "Ar36"\noccupant = "Ar36"\n'
    text = text.replace(
        state_anchor,
        state_anchor + '[[structure.kinds.states]]\nname = "Xe"\noccupant = "Xe"\n',
    )
    if text.count('name = "Xe"\noccupant = "Xe"') != 2:
        raise ValueError("Xe must be available in both gallery kinds")

    init_anchor = '[[structure.init]]\nname = "extended-defect-zones"'
    xe_init = """[[structure.init]]
name = "xe-screen-A"
center = { kind = "K_gallery_A", state = ["K"] }
probability = 0.05
set = "Xe"
[[structure.init]]
name = "xe-screen-B"
center = { kind = "K_gallery_B", state = ["K"] }
probability = 0.05
set = "Xe"
"""
    text = text.replace(init_anchor, xe_init + init_anchor, 1)

    # E3a values remain screens because their NEBs were incomplete-convergence.
    for isotope in ("Ar40", "Ar39", "Ar36"):
        text = _replace_rule_rate(
            text, f"hop_pristine_{isotope}_", E3A_AR_LOCAL_KCAL_MOL
        )
    for isotope in ("Ar40", "Ar39"):
        text = _replace_rule_rate(
            text, f"hop_extended_{isotope}_", E3A_AR_EXTENDED_KCAL_MOL
        )
    # Atmospheric 36Ar alone retains the E2 fast extended-defect proxy (40 kcal/mol),
    # making kinetic access explicit instead of relying on initialization alone.

    thermo_anchor = "[dynamics.thermo]\ntemperature = 773.15\n"
    activity = CONDITIONS[condition]
    text = text.replace(
        thermo_anchor,
        f"[dynamics.thermo]\ntemperature = {temperature_c + 273.15:.2f}\n"
        f"[dynamics.thermo.activity]\nH2O_vacancy = {activity:.12g}\n"
        "[dynamics.thermo.mu]\nH2O_vacancy = 0.0\n",
        1,
    )

    text = text.replace(
        "rate = { eyring = { dh = 57.0, ds = 0.0 } }\n[[dynamics.rules.modifiers]]",
        "rate = { eyring = { dh = 57.0, ds = 0.0 } }\n"
        'consumes = ["H2O_vacancy"]\n'
        "# H2O-vacancy activity is 1 in vacuum and 1/2000 in the ideal 2-kbar screen.\n"
        "[[dynamics.rules.modifiers]]",
        1,
    )

    directions = (
        ("A_to_B_intra", "K_gallery_A", "K_gallery_B", "gallery_intra"),
        ("A_to_B_inter", "K_gallery_A", "K_gallery_B", "gallery_inter"),
        ("B_to_A_intra", "K_gallery_B", "K_gallery_A", "gallery_intra"),
        ("B_to_A_inter", "K_gallery_B", "K_gallery_A", "gallery_inter"),
    )
    xe_rules: list[str] = [
        "# E4b Ar/Xe paired screen: E3a Xe = 90.744700 kcal/mol,",
        "# state incomplete-convergence; never a production calibration.",
    ]
    for zone in ("pristine", "extended"):
        for direction, center, neighbor, label in directions:
            xe_rules.extend(
                (
                    "[[dynamics.rules]]",
                    f'name = "hop_{zone}_Xe_{direction}"',
                    f'center = {{ kind = "{center}", state = ["Xe"] }}',
                    "guards = [",
                    f'  {{ kind = "{neighbor}", label = "{label}", state = ["vacant"], min = 1 }},',
                    f'  {{ kind = "{neighbor}", state = ["vacant"], min = 2 }},',
                    '  { kind = "OH_acceptor", label = "local_structure", state = ["O_residual"], min = 1 },',
                    f'  {{ kind = "Defect_zone", label = "local_defect", state = ["{zone}"], min = 1 }},',
                    "]",
                    f"rate = {{ eyring = {{ dh = {E3A_XE_SCREEN_KCAL_MOL:.6f}, ds = 0.0 }} }}",
                    "[[dynamics.rules.effects]]",
                    'target = "center"',
                    'set = "vacant"',
                    "[[dynamics.rules.effects]]",
                    'target = "neighbor"',
                    f'select = {{ kind = "{neighbor}", label = "{label}", state = ["vacant"] }}',
                    'set = "Xe"',
                )
            )
    for gallery in ("K_gallery_A", "K_gallery_B"):
        xe_rules.extend(
            (
                "[[dynamics.rules]]",
                f'name = "release_Xe_surface_{gallery[-1]}"',
                f'center = {{ kind = "{gallery}", state = ["Xe"] }}',
                'guards = [{ kind = "Surface_gate", label = "surface_gate", state = ["open"], min = 1 }]',
                "rate = { constant = 1.0 }",
                "[[dynamics.rules.effects]]",
                'target = "center"',
                'set = "vacant"',
            )
        )
    text = text.replace(
        '[execution]\nstrategy = "ctmc"',
        "\n".join(xe_rules) + '\n\n[execution]\nstrategy = "ctmc"',
        1,
    )

    start = text.index("[execution]\n")
    stop = text.index("[execution.stop]\n", start)
    schedule = ["[execution]", 'strategy = "ctmc"']
    for duration in schedule_durations():
        schedule.extend(
            (
                "[[execution.schedule]]",
                f"temperature = {temperature_c + 273.15:.2f}",
                f"duration = {duration:.1f}",
            )
        )
    text = text[:start] + "\n".join(schedule) + "\n" + text[stop:]
    return text


def generate_decks(deck_dir: Path) -> list[Path]:
    deck_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for temperature in TEMPERATURES_C:
        for condition in CONDITIONS:
            for dims in SIZES:
                path = (
                    deck_dir / f"e4b-{temperature}c-{condition}-{slug_dims(dims)}.toml"
                )
                path.write_text(
                    render_e4b_deck(dims, temperature, condition), encoding="utf-8"
                )
                paths.append(path)
    return paths


def petra_command(
    petra_bin: Path, deck: Path, out: Path, replicas: int, seed: int
) -> list[str]:
    return [
        "nice",
        "-n",
        "10",
        str(petra_bin),
        str(deck),
        "--seed",
        str(seed),
        "--ensemble",
        str(replicas),
        "--out",
        str(out),
        "--paranoid",
    ]


def run_one(command: list[str], cwd: Path, log: Path) -> float:
    env = os.environ.copy()
    for key in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "RAYON_NUM_THREADS",
    ):
        env[key] = "8"
    start = time.monotonic()
    with log.open("w", encoding="utf-8") as handle:
        # argv sequence, shell=False: trusted `nice` plus a local petra binary
        # and static thread-limit env keys for isothermal discriminants.
        # nosemgrep: python.lang.security.audit.dangerous-subprocess-use-audit
        subprocess.run(
            command,
            cwd=cwd,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=True,
            timeout=1800,
            shell=False,
        )
    return time.monotonic() - start


def run_campaign(
    petra_root: Path,
    petra_bin: Path,
    deck_dir: Path,
    raw: Path,
    replicas: int,
    seed: int,
) -> None:
    petra_root = petra_root.resolve()
    petra_bin = petra_bin.resolve()
    deck_dir = deck_dir.resolve()
    raw = raw.resolve()
    if raw.exists():
        raise ValueError(f"refusing to overwrite {raw}")
    raw.mkdir(parents=True)
    logs = raw / "logs"
    logs.mkdir()
    receipts = []
    for temperature in TEMPERATURES_C:
        for condition in CONDITIONS:
            for dims in SIZES:
                stem = f"{temperature}c-{condition}-{slug_dims(dims)}"
                deck = deck_dir / f"e4b-{stem}.toml"
                if deck.read_text(encoding="utf-8") != render_e4b_deck(
                    dims, temperature, condition
                ):
                    raise ValueError(f"tracked deck drift: {deck}")
                out = raw / stem
                elapsed = run_one(
                    petra_command(petra_bin, deck, out, replicas, seed),
                    petra_root,
                    logs / f"{stem}.log",
                )
                replays = []
                for suffix in ("a", "b"):
                    replay = raw / f"replay-{stem}-{suffix}"
                    run_one(
                        petra_command(petra_bin, deck, replay, 2, seed),
                        petra_root,
                        logs / f"replay-{stem}-{suffix}.log",
                    )
                    replays.append(replay)
                files = ("ensemble.csv", "ensemble-summary.csv", "observables.csv")
                replay_ok = all(
                    (replays[0] / name).read_bytes() == (replays[1] / name).read_bytes()
                    for name in files
                )
                if not replay_ok:
                    raise RuntimeError(f"same-seed replay diverged: {stem}")
                rows = list(
                    csv.DictReader((out / "ensemble.csv").open(encoding="utf-8"))
                )
                receipts.append(
                    {
                        "temperature_c": temperature,
                        "condition": condition,
                        "dims": list(dims),
                        "replicas": replicas,
                        "seeds": [int(row["seed"]) for row in rows],
                        "elapsed_seconds": elapsed,
                        "total_events": sum(int(row["steps"]) for row in rows),
                        "replay_verified": True,
                    }
                )
    (raw / "campaign.json").write_text(
        json.dumps(receipts, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def percentile(values: list[float], q: float) -> float:
    values = sorted(values)
    if not values:
        return math.nan
    return values[round((len(values) - 1) * q)]


def mean_ci95(
    values: list[float], seed: int
) -> tuple[float | None, float | None, float | None]:
    """Return a deterministic nonparametric bootstrap CI for the replica mean."""
    if not values:
        return None, None, None
    mean = statistics.fmean(values)
    if len(values) == 1:
        return mean, mean, mean
    rng = random.Random(seed)
    means = [
        statistics.fmean(values[rng.randrange(len(values))] for _ in values)
        for _ in range(2_000)
    ]
    return mean, percentile(means, 0.025), percentile(means, 0.975)


def state_map(deck: dict) -> dict[int, str]:
    result: dict[int, str] = {}
    index = 0
    for kind in deck["structure"]["kinds"]:
        for state in kind["states"]:
            name = str(state["name"])
            if name in {"Ar36", "Ar40", "Xe"}:
                result[index] = name
            index += 1
    return result


def population_samples(
    deck_path: Path, observables: Path
) -> dict[int, list[tuple[float, dict[str, int]]]]:
    with deck_path.open("rb") as handle:
        deck = tomllib.load(handle)
    mapping = state_map(deck)
    grouped: dict[tuple[int, float], dict[str, int]] = {}
    with observables.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["kind"] != "state_counts":
                continue
            key = (int(row["replica"]), float(row["time"]))
            totals = grouped.setdefault(key, {"Ar36": 0, "Ar40": 0, "Xe": 0})
            species = mapping.get(int(row["index"]))
            if species:
                totals[species] += int(float(row["value"]))
    result: dict[int, list[tuple[float, dict[str, int]]]] = defaultdict(list)
    for (replica, sample_time), totals in sorted(grouped.items()):
        result[replica].append((sample_time, totals))
    return dict(result)


def boundary(rows: list[tuple[float, dict[str, int]]], target: float) -> dict[str, int]:
    tolerance = max(1e-9, target * 1e-9)
    matches = [value for when, value in rows if abs(when - target) <= tolerance]
    if not matches:
        raise ValueError(f"missing schedule boundary {target}")
    return matches[-1]


def replica_series(deck: Path, observables: Path) -> list[dict[str, float | int | str]]:
    samples = population_samples(deck, observables)
    durations = schedule_durations()
    rows: list[dict[str, float | int | str]] = []
    for replica, points in samples.items():
        initial = boundary(points, 0.0)
        previous_tau = {species: 0.0 for species in initial}
        previous_time = 0.0
        elapsed = 0.0
        for sqrt_time, duration in zip(SQRT_TIME_MIN_HALF, durations, strict=True):
            elapsed += duration
            current = boundary(points, elapsed)
            metrics: dict[str, float] = {}
            for species in ("Ar36", "Ar40", "Xe"):
                if initial[species] <= 0:
                    fraction = math.nan
                    da2 = math.nan
                else:
                    fraction = (initial[species] - current[species]) / initial[species]
                    if fraction <= 0.0:
                        tau = 0.0
                    elif fraction >= 1.0:
                        tau = math.nan
                    else:
                        tau = invert_cylinder_fraction(fraction)
                    if math.isfinite(tau):
                        da2 = max(
                            0.0,
                            (tau - previous_tau[species]) / (elapsed - previous_time),
                        )
                        previous_tau[species] = tau
                    else:
                        da2 = math.nan
                metrics[f"{species.lower()}_fraction"] = fraction
                metrics[f"{species.lower()}_da2_per_s"] = da2
            a36, a40 = metrics["ar36_da2_per_s"], metrics["ar40_da2_per_s"]
            metrics["ar36_ar40_da2_ratio"] = (
                a36 / a40 if a36 > 0 and a40 > 0 else math.nan
            )
            rows.append(
                {
                    "replica": replica,
                    "sqrt_time_min_half": sqrt_time,
                    "time_s": elapsed,
                    **metrics,
                }
            )
            previous_time = elapsed
    return rows


def aggregate(raw: Path, deck_dir: Path) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for temperature_index, temperature in enumerate(TEMPERATURES_C):
        for condition_index, condition in enumerate(CONDITIONS):
            for dims_index, dims in enumerate(SIZES):
                stem = f"{temperature}c-{condition}-{slug_dims(dims)}"
                rows = replica_series(
                    deck_dir / f"e4b-{stem}.toml", raw / stem / "observables.csv"
                )
                for time_index, sqrt_time in enumerate(SQRT_TIME_MIN_HALF):
                    selected = [
                        row for row in rows if row["sqrt_time_min_half"] == sqrt_time
                    ]
                    item: dict[str, object] = {
                        "temperature_c": temperature,
                        "condition": condition,
                        "dims": slug_dims(dims),
                        "sites": math.prod(dims) * SITES_PER_CELL,
                        "replicas": len(selected),
                        "sqrt_time_min_half": sqrt_time,
                        "time_s": selected[0]["time_s"],
                    }
                    for metric_index, metric in enumerate(
                        (
                            "ar40_fraction",
                            "ar40_da2_per_s",
                            "ar36_fraction",
                            "ar36_da2_per_s",
                            "ar36_ar40_da2_ratio",
                            "xe_fraction",
                            "xe_da2_per_s",
                        )
                    ):
                        values = [
                            float(row[metric])
                            for row in selected
                            if math.isfinite(float(row[metric]))
                        ]
                        seed = (
                            0xE4B000
                            + temperature_index * 100_000
                            + condition_index * 10_000
                            + dims_index * 1_000
                            + time_index * 100
                            + metric_index
                        )
                        mean, low, high = mean_ci95(values, seed)
                        item[f"{metric}_n"] = len(values)
                        item[f"{metric}_mean"] = mean
                        item[f"{metric}_ci95_low"] = low
                        item[f"{metric}_ci95_high"] = high
                    output.append(item)
    return output


def _rise_fall(rows: list[dict[str, object]]) -> bool:
    values = [
        float(row["ar40_da2_per_s_mean"])
        for row in rows
        if row["ar40_da2_per_s_mean"] not in (None, 0.0)
        and math.isfinite(float(row["ar40_da2_per_s_mean"]))
    ]
    if len(values) < 4:
        return False
    peak = max(values[: max(2, len(values) // 2)])
    return peak >= 1.5 * values[0] and peak >= 5.0 * statistics.median(values[-3:])


def verdict(rows: list[dict[str, object]]) -> dict[str, object]:
    vacuum = [row for row in rows if row["condition"] == "vacuum"]
    claim1_cells = {}
    claim2_cells = {}
    for temperature in TEMPERATURES_C:
        for dims in SIZES:
            key = f"{temperature}c-{slug_dims(dims)}"
            selected = [
                row
                for row in vacuum
                if row["temperature_c"] == temperature
                and row["dims"] == slug_dims(dims)
            ]
            claim1_cells[key] = _rise_fall(selected)
            ratio_points = [
                {
                    "sqrt_time_min_half": int(row["sqrt_time_min_half"]),
                    "mean": float(row["ar36_ar40_da2_ratio_mean"]),
                    "effective_n": int(row["ar36_ar40_da2_ratio_n"]),
                }
                for row in selected
                if row["ar36_ar40_da2_ratio_mean"] is not None
                and float(row["ar36_ar40_da2_ratio_mean"]) > 0
                and int(row["ar36_ar40_da2_ratio_n"]) >= 4
            ]
            claim2_cells[key] = {
                "early": ratio_points[0] if ratio_points else None,
                "late": ratio_points[-1] if ratio_points else None,
                "supported_points": ratio_points,
            }
    reproduced1 = sum(claim1_cells.values())
    ratio_matches = sum(
        1
        for cell in claim2_cells.values()
        if cell["early"] is not None
        and cell["early"]["sqrt_time_min_half"] <= 10
        and 10 <= cell["early"]["mean"] <= 1e5
        and cell["late"] is not None
        and cell["late"]["sqrt_time_min_half"] >= 30
        and 0.1 <= cell["late"]["mean"] <= 10
    )
    largest = slug_dims(SIZES[-1])
    v700 = [
        row
        for row in rows
        if row["temperature_c"] == 700
        and row["dims"] == largest
        and row["sqrt_time_min_half"] == 60
    ]
    by_condition = {str(row["condition"]): row for row in v700}
    vacuum_fraction = float(by_condition["vacuum"]["ar40_fraction_mean"] or 0)
    hydro_fraction = float(
        by_condition["hydrothermal-2kbar-screen"]["ar40_fraction_mean"] or 0
    )
    xe_pairs = [
        row
        for row in vacuum
        if row["xe_fraction_mean"] is not None and row["ar40_fraction_mean"] is not None
    ]
    decoupled = any(
        abs(float(row["xe_fraction_mean"]) - float(row["ar40_fraction_mean"])) >= 0.10
        for row in xe_pairs
    )
    return {
        "scope": "provenance-labelled mechanism screen; no fitted parameters",
        "section_5_claims": {
            "1_two_stage_non_fickian_loss": {
                "verdict": "reproduced"
                if reproduced1 == 6
                else "partially_reproduced"
                if reproduced1
                else "not_reproduced",
                "passing_cells": reproduced1,
                "of": 6,
                "details": claim1_cells,
            },
            "2_distinct_reservoir_diffusivity_ratio": {
                "verdict": "reproduced"
                if ratio_matches == 6
                else "partially_reproduced"
                if ratio_matches
                else "not_reproduced",
                "matching_cells": ratio_matches,
                "of": 6,
                "details": claim2_cells,
            },
            "4_ar_xe_decoupling": {
                "verdict": "partially_reproduced" if decoupled else "not_reproduced",
                "mechanism": "E3a Xe screen creates a paired mechanistic prediction, but no source-backed Villa muscovite Xe time series exists here; this cannot be a calibrated reproduction.",
            },
            "7_hydrothermal_contrast": {
                "verdict": "partially_reproduced"
                if hydro_fraction < vacuum_fraction
                else "not_reproduced",
                "vacuum_final_ar40_fraction": vacuum_fraction,
                "hydrothermal_screen_final_ar40_fraction": hydro_fraction,
                "mechanism": "ideal H2O activity proxy tests direction only; it is not a 2-kbar calibration.",
            },
        },
    }


def write_svg(
    path: Path, rows: list[dict[str, object]], observed: list[dict[str, str]]
) -> None:
    width, height = 1200, 920
    colors = {"4x4x6": "#2563eb", "8x8x6": "#16a34a", "12x12x6": "#dc2626"}
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="600" y="28" text-anchor="middle" font-family="sans-serif" font-size="21" font-weight="700">E4b isothermal cumulative release and cylinder D/a²</text>',
    ]
    panels = (
        (500, "ar40_fraction_mean", False),
        (700, "ar40_fraction_mean", False),
        (500, "ar40_da2_per_s_mean", True),
        (700, "ar40_da2_per_s_mean", True),
    )
    for index, (temperature, metric, logarithmic) in enumerate(panels):
        col, row_index = index % 2, index // 2
        x0, y0, w, h = 75 + col * 585, 70 + row_index * 410, 500, 320
        parts.append(
            f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="#fafafa" stroke="#475569"/>'
        )
        parts.append(
            f'<text x="{x0 + 250}" y="{y0 - 12}" text-anchor="middle" font-family="sans-serif" font-size="16">{temperature} °C — {"log10 D/a² (s⁻¹)" if logarithmic else "cumulative ⁴⁰Ar fraction"}</text>'
        )
        selected_all = [
            r
            for r in rows
            if r["temperature_c"] == temperature and r["condition"] == "vacuum"
        ]
        vals = [float(r[metric]) for r in selected_all if r[metric] not in (None, 0.0)]
        observed_metric = (
            "cumulative_fraction" if not logarithmic else "log10_da2_per_s"
        )
        observed_values = [
            float(source[observed_metric])
            for source in observed
            if int(source["temperature_c"]) == temperature
            and source[observed_metric] != ""
        ]
        transformed = [
            math.log10(v) if logarithmic else v for v in vals
        ] + observed_values
        ymin, ymax = (
            (min(transformed) - 0.25, max(transformed) + 0.25)
            if transformed
            else (-10, 0)
        )
        if not logarithmic:
            ymin, ymax = 0.0, max(1.0, ymax)

        def sx(
            value: float | str, x_origin: float = x0, panel_width: float = w
        ) -> float:
            return x_origin + float(value) / 60.0 * panel_width

        def sy(
            value: float,
            y_origin: float = y0,
            panel_height: float = h,
            y_minimum: float = ymin,
            y_maximum: float = ymax,
        ) -> float:
            return (
                y_origin
                + panel_height
                - (value - y_minimum) / (y_maximum - y_minimum) * panel_height
            )

        for xtick in (0, 15, 30, 45, 60):
            xpos = sx(xtick)
            parts.append(
                f'<line x1="{xpos:.1f}" y1="{y0}" x2="{xpos:.1f}" y2="{y0 + h}" stroke="#e2e8f0"/>'
            )
            parts.append(
                f'<text x="{xpos:.1f}" y="{y0 + h + 14}" text-anchor="middle" font-family="sans-serif" font-size="10">{xtick}</text>'
            )
        for tick_index in range(5):
            yvalue = ymin + (ymax - ymin) * tick_index / 4
            ypos = sy(yvalue)
            parts.append(
                f'<line x1="{x0}" y1="{ypos:.1f}" x2="{x0 + w}" y2="{ypos:.1f}" stroke="#e2e8f0"/>'
            )
            parts.append(
                f'<text x="{x0 - 7}" y="{ypos + 3:.1f}" text-anchor="end" font-family="sans-serif" font-size="10">{yvalue:.2f}</text>'
            )
        for dims, color in colors.items():
            selected = [
                r
                for r in selected_all
                if r["dims"] == dims and r[metric] not in (None, 0.0)
            ]
            points = []
            lower_points = []
            upper_points = []
            low_key = metric.replace("_mean", "_ci95_low")
            high_key = metric.replace("_mean", "_ci95_high")
            for r in selected:
                value = float(r[metric])
                value = math.log10(value) if logarithmic else value
                sqrt_time = float(r["sqrt_time_min_half"])
                points.append(f"{sx(sqrt_time):.1f},{sy(value):.1f}")
                low = r[low_key]
                high = r[high_key]
                if low not in (None, 0.0) and high not in (None, 0.0):
                    low_value = math.log10(float(low)) if logarithmic else float(low)
                    high_value = math.log10(float(high)) if logarithmic else float(high)
                    lower_points.append(f"{sx(sqrt_time):.1f},{sy(low_value):.1f}")
                    upper_points.append(f"{sx(sqrt_time):.1f},{sy(high_value):.1f}")
            if lower_points and upper_points:
                polygon = " ".join(lower_points + list(reversed(upper_points)))
                parts.append(
                    f'<polygon points="{polygon}" fill="{color}" fill-opacity="0.13" stroke="none"/>'
                )
            if points:
                parts.append(
                    f'<polyline points="{" ".join(points)}" fill="none" stroke="{color}" stroke-width="2.5"/>'
                )
        for source in observed:
            if int(source["temperature_c"]) != temperature:
                continue
            metric_name = (
                "cumulative_fraction" if not logarithmic else "log10_da2_per_s"
            )
            if source[metric_name] == "":
                continue
            x = sx(float(source["sqrt_time_min_half"]))
            yvalue = float(source[metric_name])
            y = sy(yvalue)
            parts.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.4" fill="#111827"/>')
        parts.append(
            f'<text x="{x0 + w / 2}" y="{y0 + h + 30}" text-anchor="middle" font-family="sans-serif" font-size="13">√time (min½); black = 1998 digitization ± source CSV uncertainty</text>'
        )
    parts.append(
        '<text x="600" y="900" text-anchor="middle" font-family="sans-serif" font-size="13">colored: 4×4×6 / 8×8×6 / 12×12×6 model volumes; no physical grain-size equivalence is claimed</text>'
    )
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def analyze(raw: Path, deck_dir: Path, data_csv: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    rows = aggregate(raw, deck_dir)
    with (out / "isothermal-series.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with data_csv.open(newline="", encoding="utf-8") as handle:
        observed = list(csv.DictReader(handle))
    result = verdict(rows)
    result["campaign"] = json.loads((raw / "campaign.json").read_text(encoding="utf-8"))
    (out / "verdict.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    write_svg(out / "isothermal-overlays.svg", rows, observed)
    receipts = result["campaign"]
    (out / "campaign-receipts.json").write_text(
        json.dumps(receipts, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("generate")
    p.add_argument("deck_dir", type=Path)
    p = sub.add_parser("run")
    p.add_argument("petra_root", type=Path)
    p.add_argument("petra_bin", type=Path)
    p.add_argument("deck_dir", type=Path)
    p.add_argument("raw", type=Path)
    p.add_argument("--replicas", type=int, default=REPLICAS)
    p.add_argument("--base-seed", type=int, default=BASE_SEED)
    p = sub.add_parser("analyze")
    p.add_argument("raw", type=Path)
    p.add_argument("deck_dir", type=Path)
    p.add_argument("data_csv", type=Path)
    p.add_argument("out", type=Path)
    args = parser.parse_args()
    if args.command == "generate":
        generate_decks(args.deck_dir)
    elif args.command == "run":
        run_campaign(
            args.petra_root,
            args.petra_bin,
            args.deck_dir,
            args.raw,
            args.replicas,
            args.base_seed,
        )
    else:
        analyze(args.raw, args.deck_dir, args.data_csv, args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
