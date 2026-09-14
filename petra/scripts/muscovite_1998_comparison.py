#!/usr/bin/env python3
"""Build and evaluate the E4 Sletten-Onstott (1998) comparison campaign.

The campaign integrates E3a's bounded local-dehydroxylate and extended-zone Ar
hop estimates, retains E2's 58 ± 5 kcal/mol delamination sensitivity because
E3a emitted no delamination barrier, executes deterministic replica ensembles at
the three E2b comparison volumes, and overlays the resulting synthetic spectra
on transparent image digitizations of the 1998 paper. It is a qualitative
discrimination run, not a parameter fit.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import statistics
import subprocess
import time
from pathlib import Path

from build_muscovite_full_deck import render_deck
from muscovite_full_analysis import analyze_ensemble
from muscovite_grain_size_sweep import ensemble_rows

DELAMINATION_SENSITIVITY_KCAL_MOL = {"low": 53.0, "high": 63.0}
E3A_LOCAL_DEHYDROXYLATE_HOP_KCAL_MOL = 64.095991
E3A_EXTENDED_ZONE_HOP_KCAL_MOL = 68.410690
# E2b periodic comparison volumes: a finite-size ladder, not recovered physical
# grain diameters. Overlay against the 1998 58/165/3000 µm series is qualitative;
# E4a2 owns the surface-connected lateral release front on this same ladder.
COMPARISON_SIZES = ((4, 4, 6), (8, 8, 6), (12, 12, 6))
COMPARISON_SCHEDULE = tuple(
    (temperature + 273.15, 600.0) for temperature in range(500, 1201, 50)
)
DEFAULT_REPLICAS = 8
DEFAULT_BASE_SEED = 19_980
SITES_PER_CELL = 8
REPLAY_ARTIFACTS = ("ensemble.csv", "ensemble-summary.csv", "observables.csv")
REPLICA_PREFIX_ARTIFACTS = ("ensemble.csv", "observables.csv")


def _slug(dims: tuple[int, int, int]) -> str:
    return "x".join(str(value) for value in dims)


def render_comparison_deck(dims: tuple[int, int, int], barrier_label: str) -> str:
    """Render one tracked E4 deck from the reviewed E2 mechanism."""
    try:
        barrier = DELAMINATION_SENSITIVITY_KCAL_MOL[barrier_label]
    except KeyError as exc:
        raise ValueError(
            f"unknown delamination sensitivity label: {barrier_label}"
        ) from exc
    text = render_deck(dims=dims)
    text = text.replace(
        "# E2 full muscovite mechanism: scheduled, species- and defect-resolved.\n"
        "# Published anchors and proxy barriers are documented in\n"
        "# docs/program/results/E2-muscovite-full-mechanism.md.",
        "# E4 Sletten-Onstott 1998 comparison deck.\n"
        "# E3a incomplete-convergence Ar-hop estimates integrated explicitly;\n"
        "# E2 delamination/trap proxies retained where E3a emitted no value;\n"
        "# schedule/data limits: petra/data/muscovite-1998/README.md.",
        1,
    )
    text = text.replace(
        'name = "muscovite-full-mechanism"',
        f'name = "muscovite-1998-{barrier_label}-{_slug(dims)}"',
        1,
    )
    text = text.replace(
        f'name = "muscovite-full-{_slug(dims)}"',
        f'name = "muscovite-1998-{barrier_label}-{_slug(dims)}"',
        1,
    )
    text = text.replace(
        'comment = "mechanism comparison; proxy brackets are not computed kinetics"',
        (
            f'comment = "E4 qualitative comparison; retained delamination proxy '
            f'{barrier:.3f} kcal/mol; E3a Ar hops incomplete-convergence"'
        ),
        1,
    )
    text = text.replace(
        "# Explicit interface transition. 58 kcal/mol is a sensitivity proxy,\n"
        "# not a computed barrier; analysis brackets it by ±5 kcal/mol.",
        (
            f"# E2 delamination sensitivity proxy: {barrier:.3f} kcal/mol "
            f"({barrier_label}, 58 ± 5).\n"
            "# E3a state: incomplete-input-contract; no delamination barrier emitted."
        ),
        1,
    )
    text = text.replace(
        "rate = { eyring = { dh = 58.0, ds = 0.0 } }",
        f"rate = {{ eyring = {{ dh = {barrier:.3f}, ds = 0.0 }} }}",
        1,
    )
    pristine_rate = "rate = { eyring = { dh = 66.0, ds = 0.0 } }"
    extended_rate = "rate = { eyring = { dh = 40.0, ds = 0.0 } }"
    if text.count(pristine_rate) != 12 or text.count(extended_rate) != 12:
        raise ValueError("unexpected E2 Ar-hop rule inventory")
    text = text.replace(
        "# pristine reservoir hops: 66 kcal/mol. Extended-zone value is\n"
        "# the Nteme-2023 six-order enhancement mapped to a proxy barrier.",
        "# E3a local-dehydroxylate Ar hop: 64.095991 kcal/mol.\n"
        "# State: incomplete-convergence; bounded classical sensitivity value.",
        1,
    )
    text = text.replace(
        "# extended reservoir hops: 40 kcal/mol. Extended-zone value is\n"
        "# the Nteme-2023 six-order enhancement mapped to a proxy barrier.",
        "# E3a extended-zone Ar hop: 68.410690 kcal/mol (lowest geometry).\n"
        "# State: incomplete-convergence; bounded classical sensitivity value.",
        1,
    )
    text = text.replace(
        pristine_rate,
        f"rate = {{ eyring = {{ dh = {E3A_LOCAL_DEHYDROXYLATE_HOP_KCAL_MOL}, ds = 0.0 }} }}",
    )
    text = text.replace(
        extended_rate,
        f"rate = {{ eyring = {{ dh = {E3A_EXTENDED_ZONE_HOP_KCAL_MOL}, ds = 0.0 }} }}",
    )
    schedule_start = text.index("[execution]\n")
    stop_start = text.index("[execution.stop]\n", schedule_start)
    schedule = ["[execution]", 'strategy = "ctmc"']
    for temperature_k, duration_s in COMPARISON_SCHEDULE:
        schedule.extend(
            (
                "[[execution.schedule]]",
                f"temperature = {temperature_k:.2f}",
                f"duration = {duration_s:.1f}",
            )
        )
    text = text[:schedule_start] + "\n".join(schedule) + "\n" + text[stop_start:]
    return text


def generate_decks(deck_dir: Path) -> list[Path]:
    deck_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for barrier_label in DELAMINATION_SENSITIVITY_KCAL_MOL:
        for dims in COMPARISON_SIZES:
            path = deck_dir / f"muscovite-1998-{barrier_label}-{_slug(dims)}.toml"
            path.write_text(
                render_comparison_deck(dims, barrier_label), encoding="utf-8"
            )
            paths.append(path)
    return paths


def _run(command: list[str], cwd: Path, log_path: Path) -> float:
    env = os.environ.copy()
    for name in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "RAYON_NUM_THREADS",
    ):
        env[name] = "8"
    started = time.monotonic()
    with log_path.open("w", encoding="utf-8") as log:
        # argv sequence, shell=False: trusted `nice` plus a local petra binary path.
        # nosemgrep: python.lang.security.audit.dangerous-subprocess-use-audit
        subprocess.run(
            command,
            cwd=cwd,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
            timeout=1_800,
            shell=False,
        )
    return time.monotonic() - started


def _petra_command(
    petra_bin: Path, deck: Path, output: Path, replicas: int, seed: int
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
        str(output),
        "--paranoid",
    ]


def _artifact_evidence(path: Path) -> dict[str, int | str]:
    """Return independently re-checkable size and SHA-256 evidence for one file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {"bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def _verify_artifact(path: Path, expected: dict[str, int | str]) -> None:
    actual = _artifact_evidence(path)
    if actual != expected:
        raise RuntimeError(
            f"artifact drift for {path}: expected {expected}, observed {actual}"
        )


def _verify_prefix_replay_artifacts(
    raw_root: Path, stem: str, receipt: dict[str, object]
) -> None:
    recorded = receipt["primary_prefix_replay_artifacts"]
    computed = {
        filename: {
            "primary": _replica_prefix_evidence(raw_root / stem / filename),
            "replay_a": _replica_prefix_evidence(
                raw_root / f"replay-{stem}-a" / filename
            ),
            "replay_b": _replica_prefix_evidence(
                raw_root / f"replay-{stem}-b" / filename
            ),
        }
        for filename in REPLICA_PREFIX_ARTIFACTS
    }
    if computed != recorded:
        raise RuntimeError(
            f"primary_prefix_replay_artifacts drift for {stem}: "
            f"expected {recorded}, observed {computed}"
        )


def _replica_prefix_evidence(path: Path, replicas: int = 2) -> dict[str, int | str]:
    """Hash canonical CSV rows for the first replicas, independent of ensemble size."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"CSV header missing from {path}")
        rows = list(reader)
    if "replica" in reader.fieldnames:
        rows = [row for row in rows if int(row["replica"]) < replicas]
    elif "seed" in reader.fieldnames:
        rows = rows[:replicas]
    else:
        raise ValueError(f"replica and seed columns missing from {path}")
    encoded = json.dumps(rows, sort_keys=True, separators=(",", ":")).encode()
    return {"rows": len(rows), "sha256": hashlib.sha256(encoded).hexdigest()}


def _replay_evidence(
    replay_paths: list[Path],
) -> dict[str, dict[str, dict[str, int | str]]]:
    """Prove both same-seed replay directories contain byte-identical artifacts."""
    if len(replay_paths) != 2:
        raise ValueError("same-seed verification requires exactly two replay paths")
    evidence = {
        filename: {
            "replay_a": _artifact_evidence(replay_paths[0] / filename),
            "replay_b": _artifact_evidence(replay_paths[1] / filename),
        }
        for filename in REPLAY_ARTIFACTS
    }
    mismatches = [
        filename
        for filename, pair in evidence.items()
        if pair["replay_a"] != pair["replay_b"]
    ]
    if mismatches:
        raise RuntimeError(
            "same-seed replay diverged for artifacts: " + ", ".join(mismatches)
        )
    return evidence


def run_campaign(
    petra_root: Path,
    petra_bin: Path,
    deck_dir: Path,
    raw_root: Path,
    replicas: int = DEFAULT_REPLICAS,
    base_seed: int = DEFAULT_BASE_SEED,
) -> None:
    if replicas < 2:
        raise ValueError("E4 campaign requires at least two replicas")
    if raw_root.exists():
        raise ValueError(f"refusing to overwrite campaign root: {raw_root}")
    raw_root.mkdir(parents=True)
    logs = raw_root / "logs"
    logs.mkdir()
    receipts = []
    for barrier_label, barrier in DELAMINATION_SENSITIVITY_KCAL_MOL.items():
        for dims in COMPARISON_SIZES:
            stem = f"{barrier_label}-{_slug(dims)}"
            deck = deck_dir / f"muscovite-1998-{stem}.toml"
            if deck.read_text(encoding="utf-8") != render_comparison_deck(
                dims, barrier_label
            ):
                raise ValueError(f"tracked deck drift: {deck}")
            output = raw_root / stem
            elapsed = _run(
                _petra_command(petra_bin, deck, output, replicas, base_seed),
                petra_root,
                logs / f"{stem}.log",
            )
            replay_paths = []
            for suffix in ("a", "b"):
                replay = raw_root / f"replay-{stem}-{suffix}"
                _run(
                    _petra_command(petra_bin, deck, replay, 2, base_seed),
                    petra_root,
                    logs / f"replay-{stem}-{suffix}.log",
                )
                replay_paths.append(replay)
            replay_artifacts = _replay_evidence(replay_paths)
            primary_prefix_replay_artifacts = {}
            for filename in REPLICA_PREFIX_ARTIFACTS:
                primary_prefix = _replica_prefix_evidence(output / filename)
                replay_a_prefix = _replica_prefix_evidence(replay_paths[0] / filename)
                replay_b_prefix = _replica_prefix_evidence(replay_paths[1] / filename)
                if not primary_prefix == replay_a_prefix == replay_b_prefix:
                    raise RuntimeError(
                        f"primary/replay replica-prefix divergence for {stem}/{filename}"
                    )
                primary_prefix_replay_artifacts[filename] = {
                    "primary": primary_prefix,
                    "replay_a": replay_a_prefix,
                    "replay_b": replay_b_prefix,
                }
            rows = ensemble_rows(output / "ensemble.csv", replicas)
            primary_command = _petra_command(
                petra_bin, deck, output, replicas, base_seed
            )
            receipt = {
                "barrier_label": barrier_label,
                "barrier_kcal_mol": barrier,
                "dims": list(dims),
                "sites": math.prod(dims) * SITES_PER_CELL,
                "replicas": replicas,
                "seeds": [int(row["seed"]) for row in rows],
                "elapsed_seconds": elapsed,
                "total_events": sum(int(row["steps"]) for row in rows),
                "replay_verified": True,
                "command": primary_command,
                "deck_artifact": _artifact_evidence(deck),
                "petra_binary_artifact": _artifact_evidence(petra_bin),
                "primary_artifacts": {
                    filename: _artifact_evidence(output / filename)
                    for filename in REPLAY_ARTIFACTS
                },
                "replay_artifacts": replay_artifacts,
                "primary_prefix_replay_artifacts": primary_prefix_replay_artifacts,
            }
            receipts.append(receipt)
            (output / "receipt.json").write_text(
                json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
            )
    (raw_root / "campaign.json").write_text(
        json.dumps(receipts, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _finite(value: float) -> float | None:
    return value if math.isfinite(value) else None


def _read_release_data(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required = {
        "experiment",
        "grain_size_um",
        "temperature_c",
        "release_rate_x1e6_per_c",
    }
    if not rows or not required <= set(rows[0]):
        raise ValueError(f"invalid release digitization: {path}")
    return rows


def _read_age_data(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required = {"experiment", "fraction_ar39_start", "fraction_ar39_end", "age_ma"}
    if not rows or not required <= set(rows[0]):
        raise ValueError(f"invalid age digitization: {path}")
    return rows


def _synthetic_rows(
    results: dict[tuple[str, tuple[int, int, int]], object],
) -> list[dict[str, object]]:
    rows = []
    for (barrier_label, dims), result in results.items():
        previous = 0.0
        for step in result.steps:
            cumulative = step.cumulative_ar39_fraction.mean
            rows.append(
                {
                    "barrier_label": barrier_label,
                    "barrier_kcal_mol": DELAMINATION_SENSITIVITY_KCAL_MOL[
                        barrier_label
                    ],
                    "dims": _slug(dims),
                    "sites": result.sites,
                    "replicas": result.replicas,
                    "segment": step.segment,
                    "temperature_c": step.temperature_k - 273.15,
                    "duration_s": step.duration_s,
                    "released_ar40_mean": step.released_ar40.mean,
                    "released_ar40_ci95_low": step.released_ar40.ci95[0],
                    "released_ar40_ci95_high": step.released_ar40.ci95[1],
                    "released_ar39_mean": step.released_ar39.mean,
                    "released_ar36_mean": step.released_ar36.mean,
                    "cumulative_ar40_fraction_mean": step.cumulative_ar40_fraction.mean,
                    "cumulative_ar39_fraction_mean": cumulative,
                    "ar39_fraction_start": previous,
                    "apparent_age_ma_mean": _finite(step.apparent_age_ma.mean),
                    "apparent_age_ma_ci95_low": _finite(step.apparent_age_ma.ci95[0]),
                    "apparent_age_ma_ci95_high": _finite(step.apparent_age_ma.ci95[1]),
                    "apparent_age_defined_replicas": len(step.apparent_age_ma.values),
                    "ar36_ar40_mean": _finite(step.ar36_ar40.mean),
                }
            )
            previous = cumulative
    return rows


def _monotonic_volume_crossover(peaks: dict[str, float]) -> bool:
    small, medium, large = (_slug(dims) for dims in COMPARISON_SIZES)
    return peaks[small] < peaks[medium] < peaks[large]


def _peak_temperature(
    rows: list[dict[str, object]], dims: tuple[int, int, int]
) -> float:
    selected = [row for row in rows if row["dims"] == _slug(dims)]
    return float(
        max(selected, key=lambda row: float(row["released_ar40_mean"]))["temperature_c"]
    )


def _old_step_signature(
    rows: list[dict[str, object]], dims: tuple[int, int, int]
) -> tuple[bool, float | None]:
    ages = [
        float(row["apparent_age_ma_mean"])
        for row in rows
        if row["dims"] == _slug(dims) and row["apparent_age_ma_mean"] is not None
    ]
    if len(ages) < 3:
        return False, None
    tail = statistics.median(ages[len(ages) // 2 :])
    excess = max(ages[: max(1, len(ages) // 3)]) - tail
    return excess >= 50.0, excess


def _recoil_signature(
    rows: list[dict[str, object]], dims: tuple[int, int, int]
) -> bool:
    selected = [row for row in rows if row["dims"] == _slug(dims)]
    defined = [row for row in selected if row["apparent_age_ma_mean"] is not None]
    if len(defined) < 3:
        return False
    early_age_drop = float(defined[0]["apparent_age_ma_mean"]) > float(
        defined[1]["apparent_age_ma_mean"]
    )
    early_ar39 = float(defined[0]["released_ar39_mean"])
    later_ar39 = max(float(row["released_ar39_mean"]) for row in defined[1:])
    return early_age_drop and early_ar39 < later_ar39


def evaluate(
    synthetic: list[dict[str, object]], release_data: list[dict[str, str]]
) -> dict[str, object]:
    by_barrier: dict[str, dict[str, object]] = {}
    observed_groups: dict[str, list[dict[str, str]]] = {}
    for row in release_data:
        if row["experiment"].startswith("500C_isothermal_"):
            observed_groups.setdefault(row["experiment"], []).append(row)
    observed_peaks = {
        name: float(
            max(rows, key=lambda row: float(row["release_rate_x1e6_per_c"]))[
                "temperature_c"
            ]
        )
        for name, rows in observed_groups.items()
    }
    old_step_count = 0
    for barrier_label in DELAMINATION_SENSITIVITY_KCAL_MOL:
        rows = [row for row in synthetic if row["barrier_label"] == barrier_label]
        peaks = {
            _slug(dims): _peak_temperature(rows, dims) for dims in COMPARISON_SIZES
        }
        finals = {
            _slug(dims): float(
                [row for row in rows if row["dims"] == _slug(dims)][-1][
                    "cumulative_ar40_fraction_mean"
                ]
            )
            for dims in COMPARISON_SIZES
        }
        old_steps = {
            _slug(dims): _old_step_signature(rows, dims) for dims in COMPARISON_SIZES
        }
        old_step_count += sum(bool(value[0]) for value in old_steps.values())
        recoil = {
            _slug(dims): _recoil_signature(rows, dims) for dims in COMPARISON_SIZES
        }
        by_barrier[barrier_label] = {
            "peak_temperature_c": peaks,
            "old_initial_step": {
                key: {"reproduced": value[0], "excess_ma": value[1]}
                for key, value in old_steps.items()
            },
            "grain_size_crossover_reproduced": _monotonic_volume_crossover(peaks),
            "recoil_driven_distortion_reproduced": recoil,
            "high_temperature_merge_reproduced": max(finals.values())
            - min(finals.values())
            <= 0.05,
            "final_cumulative_ar40_fraction": finals,
        }
    crossover_count = sum(
        bool(item["grain_size_crossover_reproduced"]) for item in by_barrier.values()
    )
    if crossover_count == len(by_barrier):
        crossover_verdict = "reproduced"
        crossover_mechanism = (
            "surface-connected lateral release recovers the monotonic "
            "volume crossover in both retained proxy brackets"
        )
    elif crossover_count:
        crossover_verdict = "partially_reproduced"
        crossover_mechanism = (
            "surface-connected lateral release recovers the monotonic "
            "volume crossover in only one retained proxy bracket"
        )
    else:
        crossover_verdict = "not_reproduced"
        crossover_mechanism = (
            "surface-connected lateral release leaves the volume peaks unordered; "
            "next discriminate isothermal reservoir kinetics rather than tune the front"
        )
    family_count = len(by_barrier) * len(COMPARISON_SIZES)
    if old_step_count == family_count:
        old_step_verdict = "reproduced"
    elif old_step_count:
        old_step_verdict = "partially_reproduced"
    else:
        old_step_verdict = "not_reproduced"
    return {
        "scope": "qualitative discrimination, not fit",
        "observed": {
            "peak_temperature_c": observed_peaks,
            "grain_size_crossover": (
                observed_peaks["500C_isothermal_small"]
                < observed_peaks["500C_isothermal_medium"]
                < observed_peaks["500C_isothermal_large"]
            ),
            "old_initial_steps": True,
            "high_temperature_merge": True,
        },
        "section_5_claims": {
            "1_two_stage_non_fickian_loss": {
                "verdict": "not_reproduced",
                "mechanism": "surface-connected lateral release still does not recover a resolvable late stage",
            },
            "2_distinct_reservoir_diffusivity_ratio": {
                "verdict": "not_reproduced",
                "mechanism": "initialization alone does not create the measured kinetic separation",
            },
            "3_homogeneous_age_staircase": {
                "verdict": "partially_reproduced",
                "mechanism": "ordered steps emerge but magnitude and late structure do not",
            },
            "4_ar_xe_decoupling": {
                "verdict": "not_reproduced",
                "mechanism": "E2 has no complete Xe state and release mechanism",
            },
            "5_recoil_old_initial_steps": {
                "verdict": old_step_verdict,
                "mechanism": (
                    f"recoil distortion emerges; {old_step_count} of {family_count} "
                    "volume/sensitivity families cross the old-step gate"
                ),
            },
            "6_grain_size_delamination_fraction": {
                "verdict": crossover_verdict,
                "mechanism": crossover_mechanism,
            },
            "7_hydrothermal_contrast": {
                "verdict": "not_reproduced",
                "mechanism": "no H2O chemical-potential control exists in the current deck",
            },
        },
        "synthetic": by_barrier,
    }


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("cannot write empty synthetic table")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _polyline(points: list[tuple[float, float]], color: str, dash: str = "") -> str:
    encoded = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
    dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<polyline points="{encoded}" fill="none" stroke="{color}" stroke-width="2"{dash_attr}/>'


def _band_polygon(points: list[tuple[float, float]], color: str) -> str:
    encoded = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
    return f'<polygon points="{encoded}" fill="{color}" fill-opacity="0.13" stroke="none"/>'


def comparison_svg(
    synthetic: list[dict[str, object]],
    release_data: list[dict[str, str]],
    age_data: list[dict[str, str]],
) -> str:
    width, height = 1100, 820
    colors = {"4x4x6": "#2563eb", "8x8x6": "#d97706", "12x12x6": "#059669"}
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        f'<rect width="{width}" height="{height}" fill="#fbfaf7"/>',
        '<text x="550" y="30" text-anchor="middle" font-family="sans-serif" font-size="21" font-weight="700">E4a2 lateral-front comparison: Sletten–Onstott 1998 vs Petra</text>',
    ]
    # Panel A: release profile, normalized because raster digitization is in rate units.
    left, top, plot_w, plot_h = 75, 70, 950, 280
    parts += [
        f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="white" stroke="#334155"/>',
        '<text x="550" y="58" text-anchor="middle" font-family="sans-serif" font-size="15">Normalized 40Ar* release (solid: Fig. 4; dashed: 53 kcal delam proxy + E3a Ar hops)</text>',
    ]
    obs_map = {
        "500C_isothermal_small": "4x4x6",
        "500C_isothermal_medium": "8x8x6",
        "500C_isothermal_large": "12x12x6",
    }
    for experiment, dims_slug in obs_map.items():
        rows = [row for row in release_data if row["experiment"] == experiment]
        maximum = max(float(row["release_rate_x1e6_per_c"]) for row in rows) or 1.0
        points = [
            (
                left + (float(row["temperature_c"]) - 500.0) / 700.0 * plot_w,
                top + plot_h * (1.0 - float(row["release_rate_x1e6_per_c"]) / maximum),
            )
            for row in rows
        ]
        parts.append(_polyline(points, colors[dims_slug]))
    for dims_slug, color in colors.items():
        rows = [
            row
            for row in synthetic
            if row["barrier_label"] == "low" and row["dims"] == dims_slug
        ]
        maximum = max(float(row["released_ar40_mean"]) for row in rows) or 1.0
        upper = []
        lower = []
        points = [
            (
                left + (float(row["temperature_c"]) - 500.0) / 700.0 * plot_w,
                top + plot_h * (1.0 - float(row["released_ar40_mean"]) / maximum),
            )
            for row in rows
        ]
        for row in rows:
            x = left + (float(row["temperature_c"]) - 500.0) / 700.0 * plot_w
            upper.append(
                (
                    x,
                    top
                    + plot_h
                    * (
                        1.0
                        - min(
                            1.0,
                            max(0.0, float(row["released_ar40_ci95_high"]) / maximum),
                        )
                    ),
                )
            )
            lower.append(
                (
                    x,
                    top
                    + plot_h
                    * (
                        1.0
                        - min(
                            1.0,
                            max(0.0, float(row["released_ar40_ci95_low"]) / maximum),
                        )
                    ),
                )
            )
        parts.append(_band_polygon(upper + list(reversed(lower)), color))
        parts.append(_polyline(points, color, "7 5"))
    # Panel B: hydrothermal staircase plus largest synthetic age spectrum.
    top2 = 450
    parts += [
        f'<rect x="{left}" y="{top2}" width="{plot_w}" height="{plot_h}" fill="white" stroke="#334155"/>',
        '<text x="550" y="435" text-anchor="middle" font-family="sans-serif" font-size="15">Apparent age spectrum (black: Fig. 3d; colored: 58 ± 5 delam proxy, E3a Ar hops)</text>',
    ]
    observed = [
        row for row in age_data if row["experiment"] == "hydrothermal_700C_2kbar_45d"
    ]
    observed_points = []
    for row in observed:
        x1 = left + float(row["fraction_ar39_start"]) * plot_w
        x2 = left + float(row["fraction_ar39_end"]) * plot_w
        y = top2 + plot_h * (1.0 - float(row["age_ma"]) / 1100.0)
        observed_points.extend(((x1, y), (x2, y)))
    parts.append(_polyline(observed_points, "#111827"))
    for barrier_label, dash in (("low", ""), ("high", "6 4")):
        rows = [
            row
            for row in synthetic
            if row["barrier_label"] == barrier_label
            and row["dims"] == "12x12x6"
            and row["apparent_age_ma_mean"] is not None
        ]
        upper = []
        lower = []
        points = []
        for row in rows:
            x = left + float(row["cumulative_ar39_fraction_mean"]) * plot_w
            age = min(1100.0, max(0.0, float(row["apparent_age_ma_mean"])))
            y = top2 + plot_h * (1.0 - age / 1100.0)
            points.append((x, y))
            upper_age = min(1100.0, max(0.0, float(row["apparent_age_ma_ci95_high"])))
            lower_age = min(1100.0, max(0.0, float(row["apparent_age_ma_ci95_low"])))
            upper.append((x, top2 + plot_h * (1.0 - upper_age / 1100.0)))
            lower.append((x, top2 + plot_h * (1.0 - lower_age / 1100.0)))
        parts.append(
            _band_polygon(
                upper + list(reversed(lower)),
                "#7c3aed" if barrier_label == "low" else "#dc2626",
            )
        )
        parts.append(
            _polyline(points, "#7c3aed" if barrier_label == "low" else "#dc2626", dash)
        )
    for tick in (0, 200, 400, 600, 800, 1000, 1200):
        x = left + (tick - 500) / 700 * plot_w
        if 500 <= tick <= 1200:
            parts.append(
                f'<text x="{x:.1f}" y="{top + plot_h + 20}" text-anchor="middle" font-family="sans-serif" font-size="11">{tick}</text>'
            )
    parts += [
        f'<text x="{left + plot_w / 2}" y="{top + plot_h + 40}" text-anchor="middle" font-family="sans-serif" font-size="13">temperature (°C)</text>',
        f'<text x="{left + plot_w / 2}" y="{top2 + plot_h + 40}" text-anchor="middle" font-family="sans-serif" font-size="13">cumulative fraction 39Ar released</text>',
        '<text x="550" y="810" text-anchor="middle" font-family="sans-serif" font-size="11">Digitized points carry plot-reading uncertainty; comparison schedule is 600 s per 50 °C increment, not a recovered furnace log.</text>',
        "</svg>",
    ]
    return "\n".join(parts) + "\n"


def analyze_campaign(
    raw_root: Path, deck_dir: Path, data_dir: Path, out_dir: Path
) -> None:
    receipts = json.loads((raw_root / "campaign.json").read_text(encoding="utf-8"))
    results = {}
    for receipt in receipts:
        barrier_label = receipt["barrier_label"]
        dims = tuple(receipt["dims"])
        stem = f"{barrier_label}-{_slug(dims)}"
        deck = deck_dir / f"muscovite-1998-{stem}.toml"
        _verify_artifact(deck, receipt["deck_artifact"])
        for filename, expected in receipt["primary_artifacts"].items():
            _verify_artifact(raw_root / stem / filename, expected)
        for filename, pair in receipt["replay_artifacts"].items():
            _verify_artifact(raw_root / f"replay-{stem}-a" / filename, pair["replay_a"])
            _verify_artifact(raw_root / f"replay-{stem}-b" / filename, pair["replay_b"])
        _verify_prefix_replay_artifacts(raw_root, stem, receipt)
        if receipt["seeds"] != list(
            range(receipt["seeds"][0], receipt["seeds"][0] + receipt["replicas"])
        ):
            raise ValueError(f"non-contiguous seeds in {stem}")
        result = analyze_ensemble(
            deck,
            raw_root / stem / "observables.csv",
            dims=dims,
            j_factor=0.01,
            expected_seeds=tuple(receipt["seeds"]),
        )
        results[(barrier_label, dims)] = result
    synthetic = _synthetic_rows(results)
    release_data = _read_release_data(data_dir / "figure4-release-rates.csv")
    age_data = _read_age_data(data_dir / "figure3-and-5-age-spectra.csv")
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(out_dir / "synthetic-spectra.csv", synthetic)
    verdict = evaluate(synthetic, release_data)
    (out_dir / "comparison.json").write_text(
        json.dumps(verdict, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    (out_dir / "comparison.svg").write_text(
        comparison_svg(synthetic, release_data, age_data), encoding="utf-8"
    )
    (out_dir / "campaign-receipts.json").write_text(
        json.dumps(receipts, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    result_artifacts = {
        filename: _artifact_evidence(out_dir / filename)
        for filename in (
            "campaign-receipts.json",
            "comparison.json",
            "comparison.svg",
            "synthetic-spectra.csv",
        )
    }
    analysis_receipt = {
        "campaign_input": _artifact_evidence(raw_root / "campaign.json"),
        "data_inputs": {
            filename: _artifact_evidence(data_dir / filename)
            for filename in (
                "figure3-and-5-age-spectra.csv",
                "figure4-release-rates.csv",
            )
        },
        "script_inputs": {
            "build_muscovite_full_deck.py": _artifact_evidence(
                Path(__file__).with_name("build_muscovite_full_deck.py")
            ),
            "muscovite_1998_comparison.py": _artifact_evidence(Path(__file__)),
            "muscovite_full_analysis.py": _artifact_evidence(
                Path(__file__).with_name("muscovite_full_analysis.py")
            ),
        },
        "result_artifacts": result_artifacts,
    }
    (out_dir / "analysis-receipt.json").write_text(
        json.dumps(analysis_receipt, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(verdict, indent=2, sort_keys=True))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    generate = sub.add_parser("generate")
    generate.add_argument("deck_dir", type=Path)
    run = sub.add_parser("run")
    run.add_argument("petra_root", type=Path)
    run.add_argument("petra_bin", type=Path)
    run.add_argument("deck_dir", type=Path)
    run.add_argument("raw_root", type=Path)
    run.add_argument("--replicas", type=int, default=DEFAULT_REPLICAS)
    run.add_argument("--base-seed", type=int, default=DEFAULT_BASE_SEED)
    analyze = sub.add_parser("analyze")
    analyze.add_argument("raw_root", type=Path)
    analyze.add_argument("deck_dir", type=Path)
    analyze.add_argument("data_dir", type=Path)
    analyze.add_argument("out_dir", type=Path)
    args = parser.parse_args()
    if args.command == "generate":
        for path in generate_decks(args.deck_dir):
            print(path)
    elif args.command == "run":
        run_campaign(
            args.petra_root,
            args.petra_bin,
            args.deck_dir,
            args.raw_root,
            args.replicas,
            args.base_seed,
        )
    else:
        analyze_campaign(args.raw_root, args.deck_dir, args.data_dir, args.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
