#!/usr/bin/env python3
"""Run deterministic replica ensembles across a muscovite lattice-size ladder."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import time
from dataclasses import asdict
from itertools import pairwise
from pathlib import Path

from build_muscovite_full_deck import render_deck
from muscovite_full_analysis import (
    MUSCOVITE_SITES_PER_CELL,
    SizeRunReceipt,
    analyze_ensemble,
    assess_stability,
    write_ensemble_products,
)

DEFAULT_SIZES = (
    (4, 4, 6),
    (6, 6, 9),
    (8, 8, 12),
    (12, 12, 18),
    (16, 16, 24),
    (24, 24, 36),
    (36, 36, 54),
)
DEFAULT_REPLICAS = 32
DEFAULT_BASE_SEED = 19_980


def parse_sizes(text: str) -> tuple[tuple[int, int, int], ...]:
    sizes = []
    for item in text.split(","):
        try:
            dims = tuple(int(value) for value in item.lower().split("x"))
        except ValueError as exc:
            raise ValueError("sizes must use naxnbxnlayers syntax") from exc
        if len(dims) != 3 or any(value <= 0 for value in dims):
            raise ValueError("every size must contain three positive dimensions")
        sizes.append(dims)
    volumes = [math.prod(dims) for dims in sizes]
    if not sizes or any(current <= previous for previous, current in pairwise(volumes)):
        raise ValueError("size ladder volumes must be strictly increasing")
    return tuple(sizes)


def _nonnegative_integer(text: str, label: str) -> int:
    try:
        value = int(text)
    except ValueError as exc:
        raise ValueError(f"{label} must be a non-negative integer") from exc
    if value < 0:
        raise ValueError(f"{label} must be a non-negative integer")
    return value


def ensemble_rows(
    ensemble_csv: Path, expected_replicas: int | None = None
) -> list[dict[str, str]]:
    with ensemble_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not {"seed", "steps"} <= set(reader.fieldnames):
            raise ValueError(f"{ensemble_csv}: missing seed/steps columns")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{ensemble_csv}: missing ensemble step receipts")
    seeds = [_nonnegative_integer(row["seed"], "seed") for row in rows]
    if len(set(seeds)) != len(seeds):
        raise ValueError(f"{ensemble_csv}: duplicate seed")
    if expected_replicas is not None and len(rows) != expected_replicas:
        raise ValueError(f"{ensemble_csv}: row count does not match replica count")
    if seeds != list(range(seeds[0], seeds[0] + len(seeds))):
        raise ValueError(f"{ensemble_csv}: seeds must be ordered and contiguous")
    return rows


def total_events(ensemble_csv: Path, expected_replicas: int | None = None) -> int:
    rows = ensemble_rows(ensemble_csv, expected_replicas)
    return sum(_nonnegative_integer(row["steps"], "steps") for row in rows)


def write_json_atomic(path: Path, value) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    payload = json.dumps(value, indent=2, sort_keys=True) + "\n"
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


def write_campaign_receipts(path: Path, receipts: list[SizeRunReceipt]) -> None:
    write_json_atomic(path, [asdict(receipt) for receipt in receipts])


def load_campaign_receipts(path: Path) -> list[SizeRunReceipt]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, list) or not raw:
        raise ValueError("campaign receipts must be a non-empty list")
    receipts = []
    for index, item in enumerate(raw):
        if not isinstance(item, dict):
            raise TypeError(f"receipt {index} must be an object")
        dims = item.get("dims")
        if (
            not isinstance(dims, list)
            or len(dims) != 3
            or any(type(value) is not int or value <= 0 for value in dims)
        ):
            raise ValueError(f"receipt {index} has invalid dims")
        sites = item.get("sites")
        replicas = item.get("replicas")
        elapsed = item.get("elapsed_seconds")
        events = item.get("total_events")
        replay = item.get("replay_verified")
        if (
            type(sites) is not int
            or sites != math.prod(dims) * MUSCOVITE_SITES_PER_CELL
        ):
            raise ValueError(f"receipt {index} has invalid sites")
        if type(replicas) is not int or replicas < 2:
            raise ValueError(f"receipt {index} has invalid replicas")
        if (
            not isinstance(elapsed, (int, float))
            or isinstance(elapsed, bool)
            or not math.isfinite(elapsed)
            or elapsed <= 0
        ):
            raise ValueError(f"receipt {index} has invalid elapsed_seconds")
        if type(events) is not int or events < 0:
            raise ValueError(f"receipt {index} has invalid total_events")
        if type(replay) is not bool:
            raise ValueError(f"receipt {index} has invalid replay_verified")
        receipts.append(
            SizeRunReceipt(
                dims=tuple(dims),
                sites=sites,
                replicas=replicas,
                elapsed_seconds=float(elapsed),
                total_events=events,
                replay_verified=replay,
            )
        )
    volumes = [math.prod(receipt.dims) for receipt in receipts]
    if any(current <= previous for previous, current in pairwise(volumes)):
        raise ValueError("campaign receipt dimensions must be strictly increasing")
    if len({receipt.replicas for receipt in receipts}) != 1:
        raise ValueError("campaign receipts must use one consistent replica count")
    return receipts


def validate_result_receipt(
    result, receipt: SizeRunReceipt, ensemble_csv: Path | None = None
) -> None:
    if result.dims != receipt.dims or result.sites != receipt.sites:
        raise ValueError("analysis dimensions/site count do not match campaign receipt")
    if result.replicas != receipt.replicas:
        raise ValueError("analysis replica count does not match campaign receipt")
    if not receipt.replay_verified:
        raise ValueError("campaign receipt lacks exact-seed replay verification")
    if (
        ensemble_csv is not None
        and total_events(ensemble_csv, expected_replicas=receipt.replicas)
        != receipt.total_events
    ):
        raise ValueError("ensemble event count does not match campaign receipt")


def _run(command: list[str], log_path: Path, cwd: Path, env: dict[str, str]) -> float:
    started = time.monotonic()
    with log_path.open("w", encoding="utf-8") as log:
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
    petra_bin: Path,
    deck: Path,
    out: Path,
    replicas: int,
    base_seed: int,
) -> list[str]:
    return [
        "nice",
        "-n",
        "10",
        str(petra_bin),
        str(deck),
        "--seed",
        str(base_seed),
        "--ensemble",
        str(replicas),
        "--out",
        str(out),
    ]


def run_campaign(
    petra_root: Path,
    petra_bin: Path,
    raw_root: Path,
    sizes: tuple[tuple[int, int, int], ...],
    replicas: int,
    base_seed: int,
) -> None:
    if replicas < 2:
        raise ValueError("campaign requires at least two replicas")
    if raw_root.exists():
        raise ValueError(f"refusing to overwrite campaign root: {raw_root}")
    decks = raw_root / "decks"
    logs = raw_root / "logs"
    decks.mkdir(parents=True)
    logs.mkdir()
    env = os.environ.copy()
    for name in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "RAYON_NUM_THREADS",
    ):
        env[name] = "16"

    receipts = []
    for dims in sizes:
        slug = "x".join(str(value) for value in dims)
        deck = decks / f"muscovite-{slug}.toml"
        deck.write_text(render_deck(dims=dims), encoding="utf-8")
        output = raw_root / f"size-{slug}"
        elapsed = _run(
            _petra_command(petra_bin, deck, output, replicas, base_seed),
            logs / f"size-{slug}.log",
            petra_root,
            env,
        )

        replay_outputs = []
        for label in ("a", "b"):
            replay = raw_root / f"replay-{slug}-{label}"
            _run(
                _petra_command(petra_bin, deck, replay, 2, base_seed),
                logs / f"replay-{slug}-{label}.log",
                petra_root,
                env,
            )
            replay_outputs.append(replay)
        replay_verified = all(
            (replay_outputs[0] / filename).read_bytes()
            == (replay_outputs[1] / filename).read_bytes()
            for filename in ("ensemble.csv", "ensemble-summary.csv", "observables.csv")
        )
        if not replay_verified:
            raise RuntimeError(f"same-seed replay diverged at lattice size {slug}")

        receipt = SizeRunReceipt(
            dims=dims,
            sites=math.prod(dims) * MUSCOVITE_SITES_PER_CELL,
            replicas=replicas,
            elapsed_seconds=elapsed,
            total_events=total_events(
                output / "ensemble.csv", expected_replicas=replicas
            ),
            replay_verified=True,
        )
        receipts.append(receipt)
        write_json_atomic(output / "receipt.json", asdict(receipt))
        # Persist after every completed size so a larger-size resource ceiling
        # cannot erase the bounded campaign's already-verified receipts.
        write_campaign_receipts(raw_root / "campaign.json", receipts)


def analyze_campaign(raw_root: Path, out_dir: Path, j_factor: float) -> None:
    receipts = load_campaign_receipts(raw_root / "campaign.json")
    results = []
    campaign_seeds = None
    for receipt in receipts:
        slug = "x".join(str(value) for value in receipt.dims)
        size_receipt = json.loads(
            (raw_root / f"size-{slug}" / "receipt.json").read_text(encoding="utf-8")
        )
        expected_receipt = json.loads(json.dumps(asdict(receipt)))
        if size_receipt != expected_receipt:
            raise ValueError(f"size-{slug} receipt does not match campaign receipt")
        rows = ensemble_rows(
            raw_root / f"size-{slug}" / "ensemble.csv", receipt.replicas
        )
        seeds = tuple(_nonnegative_integer(row["seed"], "seed") for row in rows)
        if campaign_seeds is None:
            campaign_seeds = seeds
        elif seeds != campaign_seeds:
            raise ValueError("ensemble seed identities differ across lattice sizes")
        result = analyze_ensemble(
            raw_root / "decks" / f"muscovite-{slug}.toml",
            raw_root / f"size-{slug}" / "observables.csv",
            dims=receipt.dims,
            j_factor=j_factor,
            expected_seeds=seeds,
        )
        validate_result_receipt(
            result, receipt, raw_root / f"size-{slug}" / "ensemble.csv"
        )
        results.append(result)
    write_ensemble_products(results, receipts, out_dir)
    stability = assess_stability(results)
    (out_dir / "stability.json").write_text(
        json.dumps(asdict(stability), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(
        f"sizes={len(results)} largest_sites={results[-1].sites} "
        f"stabilized_at_sites={stability.stabilized_at_sites}"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run")
    run.add_argument("petra_root", type=Path)
    run.add_argument("petra_bin", type=Path)
    run.add_argument("raw_root", type=Path)
    run.add_argument(
        "--sizes",
        type=parse_sizes,
        default=DEFAULT_SIZES,
        help="comma-separated naxnbxnlayers ladder",
    )
    run.add_argument("--replicas", type=int, default=DEFAULT_REPLICAS)
    run.add_argument("--base-seed", type=int, default=DEFAULT_BASE_SEED)
    analyze = sub.add_parser("analyze")
    analyze.add_argument("raw_root", type=Path)
    analyze.add_argument("out_dir", type=Path)
    analyze.add_argument("--j-factor", type=float, default=0.01)
    args = parser.parse_args()
    if args.command == "run":
        run_campaign(
            args.petra_root,
            args.petra_bin,
            args.raw_root,
            args.sizes,
            args.replicas,
            args.base_seed,
        )
    else:
        analyze_campaign(args.raw_root, args.out_dir, args.j_factor)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
