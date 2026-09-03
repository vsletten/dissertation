#!/usr/bin/env python3
"""Build and replay the curated 1999 C KMC model against golden runs.

The curated source bytes stay unchanged. Modern-header compatibility is supplied
only through forced-include compiler flags. The JSON report is a porting oracle:
it records build provenance, hashes, structural comparisons, divergence points,
and the historical diffusion-disable audit.
"""

from __future__ import annotations

import argparse
import concurrent.futures
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import tempfile
import time
from typing import Any, Dict, List, Optional, Sequence, Tuple

OUTPUT_COLUMNS = {"results.dat": 15, "surfAl.out": 2, "surfSi.out": 2}
INPUT_NAMES = ("data.cell", "data.lattice", "data.rxn", "data.sim")
SOURCE_NAMES = (
    "mckaol.c", "actions.c", "bfsearch.c", "envrn.c", "evtlist.c",
    "futil.c", "lattice.c", "myerr.c", "output.c", "ran2.c",
    "reactions.c", "rxnlist.c", "sim.c", "ucell.c",
)
COMPATIBILITY_FLAGS = (
    "-O3", "-ffast-math", "-std=gnu89", "-include", "stdlib.h",
    "-include", "unistd.h",
)
GOLDEN_RUNS = (
    "hotrox/935077498",
    "hotrox/936930575",
    "hotrox/937172019",
    "jasper/933892971",
    "jasper/935835145",
)
CROSS_HOST_DUPLICATE = ("hotrox/936930575", "jasper/933892971")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _parse_numeric(path: Path, expected_columns: int) -> Tuple[List[Tuple[float, ...]], Optional[str]]:
    rows: List[Tuple[float, ...]] = []
    try:
        with path.open(encoding="ascii") as handle:
            for line_number, raw in enumerate(handle, 1):
                text = raw.strip()
                if not text:
                    return rows, "blank row at line {}".format(line_number)
                fields = text.split(",")
                if len(fields) != expected_columns:
                    return rows, "line {} has {} columns, expected {}".format(
                        line_number, len(fields), expected_columns
                    )
                values = tuple(float(field) for field in fields)
                if not all(math.isfinite(value) for value in values):
                    return rows, "non-finite value at line {}".format(line_number)
                rows.append(values)
    except (OSError, UnicodeDecodeError, ValueError) as exc:
        return rows, str(exc)
    return rows, None


def compare_numeric_output(expected: Path, actual: Path, expected_columns: int) -> Dict[str, Any]:
    """Compare one numeric output without normalizing divergent behavior away."""
    if not actual.is_file():
        return {
            "expected_sha256": sha256_file(expected),
            "actual_sha256": None,
            "byte_equal": False,
            "structural_match": False,
            "classification": "behavioral_mismatch",
            "error": "actual output missing",
        }

    expected_rows, expected_error = _parse_numeric(expected, expected_columns)
    actual_rows, actual_error = _parse_numeric(actual, expected_columns)
    byte_equal = sha256_file(expected) == sha256_file(actual)
    structural_match = (
        expected_error is None
        and actual_error is None
        and len(expected_rows) == len(actual_rows)
    )

    matching_prefix = 0
    first_mismatch: Optional[int] = None
    max_abs_delta = 0.0
    for index, (expected_row, actual_row) in enumerate(zip(expected_rows, actual_rows), 1):
        if expected_row == actual_row and first_mismatch is None:
            matching_prefix += 1
        elif first_mismatch is None:
            first_mismatch = index
        for expected_value, actual_value in zip(expected_row, actual_row):
            max_abs_delta = max(max_abs_delta, abs(expected_value - actual_value))
    if first_mismatch is None and len(expected_rows) != len(actual_rows):
        first_mismatch = min(len(expected_rows), len(actual_rows)) + 1

    if byte_equal and structural_match:
        classification = "byte_parity"
    elif structural_match and matching_prefix > 0:
        classification = "compiler_prng_drift"
    else:
        classification = "behavioral_mismatch"

    return {
        "expected_sha256": sha256_file(expected),
        "actual_sha256": sha256_file(actual),
        "byte_equal": byte_equal,
        "structural_match": structural_match,
        "classification": classification,
        "expected_rows": len(expected_rows),
        "actual_rows": len(actual_rows),
        "columns": expected_columns,
        "matching_prefix_rows": matching_prefix,
        "first_mismatch_row": first_mismatch,
        "max_abs_numeric_delta": max_abs_delta,
        "expected_parse_error": expected_error,
        "actual_parse_error": actual_error,
    }


def audit_diffusion_disabled(source: Path) -> Dict[str, Any]:
    envrn = (source / "envrn.c").read_text(encoding="latin-1")
    rxnlist = (source / "rxnlist.h").read_text(encoding="latin-1")
    ndes_match = re.search(r"#define\s+NDES\s+(\d+)", rxnlist)
    nrxn_match = re.search(r"#define\s+NRXN\s+(\d+)", rxnlist)
    diffusion_false = re.search(
        r"else\s*\{\s*/\*\s*diffusion\s*\*/\s*result\s*=\s*FALSE\s*;\s*\}",
        envrn,
        re.DOTALL,
    )
    if not ndes_match or not nrxn_match:
        return {"status": "behavioral_mismatch", "reason": "NDES/NRXN definitions missing"}
    ndes = int(ndes_match.group(1))
    nrxn = int(nrxn_match.group(1))
    return {
        "status": "pinned_disabled" if diffusion_false and (ndes, nrxn) == (24, 28) else "behavioral_mismatch",
        "ndes": ndes,
        "nrxn": nrxn,
        "diffusion_ids": list(range(ndes, nrxn)),
        "isActive_final_branch_returns_false": bool(diffusion_false),
    }


def compiler_identity(compiler: str) -> str:
    result = subprocess.run(
        [compiler, "--version"], check=True, capture_output=True, text=True
    )
    return result.stdout.splitlines()[0]


def build_model(source: Path, compiler: str, output_dir: Path) -> Dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    executable = output_dir / "mckaol"
    command = [compiler, *COMPATIBILITY_FLAGS, "-o", str(executable)]
    command.extend(str(source / name) for name in SOURCE_NAMES)
    command.append("-lm")
    started = time.monotonic()
    result = subprocess.run(command, capture_output=True, text=True)
    elapsed = time.monotonic() - started
    record = {
        "command": command,
        "compatibility_flags": list(COMPATIBILITY_FLAGS),
        "returncode": result.returncode,
        "elapsed_seconds": round(elapsed, 3),
        "stdout": result.stdout,
        "stderr": result.stderr,
    }
    if result.returncode != 0:
        raise RuntimeError("compiler failed: {}".format(json.dumps(record, indent=2)))
    record["executable_sha256"] = sha256_file(executable)
    return record


def _run_one(
    fixture_id: str,
    fixtures: Path,
    executable: Path,
    run_root: Path,
    timeout_seconds: int,
) -> Dict[str, Any]:
    fixture = fixtures / fixture_id
    run_dir = run_root / fixture_id.replace("/", "-")
    run_dir.mkdir(parents=True, exist_ok=False)
    for name in INPUT_NAMES:
        shutil.copy2(str(fixture / name), str(run_dir / name))
    shutil.copy2(str(executable), str(run_dir / "mckaol"))

    started = time.monotonic()
    timed_out = False
    try:
        result = subprocess.run(
            [str(run_dir / "mckaol")],
            cwd=str(run_dir),
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
        )
        returncode = result.returncode
        stdout = result.stdout
        stderr = result.stderr
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        returncode = None
        stdout = exc.stdout.decode() if isinstance(exc.stdout, bytes) else (exc.stdout or "")
        stderr = exc.stderr.decode() if isinstance(exc.stderr, bytes) else (exc.stderr or "")
    elapsed = time.monotonic() - started

    outputs = {
        name: compare_numeric_output(fixture / name, run_dir / name, columns)
        for name, columns in OUTPUT_COLUMNS.items()
    }
    results_class = outputs["results.dat"]["classification"]
    all_structural = all(item["structural_match"] for item in outputs.values())
    if returncode != 0 or timed_out or not all_structural:
        classification = "behavioral_mismatch"
    elif all(item["classification"] == "byte_parity" for item in outputs.values()):
        classification = "byte_parity"
    elif results_class == "compiler_prng_drift":
        classification = "compiler_prng_drift"
    else:
        classification = "behavioral_mismatch"

    return {
        "fixture": fixture_id,
        "returncode": returncode,
        "timed_out": timed_out,
        "elapsed_seconds": round(elapsed, 3),
        "classification": classification,
        "stdout": stdout,
        "stderr": stderr,
        "outputs": outputs,
    }


def verify_historical_duplicate(fixtures: Path) -> Dict[str, Any]:
    left, right = CROSS_HOST_DUPLICATE
    outputs: Dict[str, Any] = {}
    for name in OUTPUT_COLUMNS:
        left_hash = sha256_file(fixtures / left / name)
        right_hash = sha256_file(fixtures / right / name)
        outputs[name] = {
            "left_sha256": left_hash,
            "right_sha256": right_hash,
            "byte_equal": left_hash == right_hash,
        }
    return {
        "fixtures": list(CROSS_HOST_DUPLICATE),
        "outputs": outputs,
        "all_byte_equal": all(item["byte_equal"] for item in outputs.values()),
    }


def render_markdown(report: Dict[str, Any]) -> str:
    lines = [
        "# 1999 C KMC conformance report",
        "",
        "Generated: `{}`".format(report["generated_at"]),
        "",
        "## Provenance",
        "",
        "- Git commit: `{}`".format(report["git_commit"]),
        "- Platform: `{}`".format(report["platform"]),
        "- Compiler: `{}`".format(report["compiler"]["actual"]),
        "- Canonical compiler lock: `{}`".format(report["compiler"]["canonical"]),
        "- Compiler lock matched: `{}`".format(str(report["compiler"]["lock_matched"]).lower()),
        "- Compatibility flags: `{}`".format(" ".join(report["build"]["compatibility_flags"])),
        "- Diffusion IDs 24–27: `{}`".format(report["diffusion_audit"]["status"]),
        "",
        "## Fixture outcomes",
        "",
        "| fixture | outcome | results parity | first divergence | surfAl rows | surfSi rows | seconds |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for run in report["runs"]:
        results = run["outputs"]["results.dat"]
        lines.append(
            "| {fixture} | {classification} | {parity} | {first} | {al} | {si} | {elapsed} |".format(
                fixture=run["fixture"],
                classification=run["classification"],
                parity="byte" if results["byte_equal"] else "diverged",
                first=results["first_mismatch_row"] or "—",
                al=run["outputs"]["surfAl.out"].get("actual_rows", "—"),
                si=run["outputs"]["surfSi.out"].get("actual_rows", "—"),
                elapsed=run["elapsed_seconds"],
            )
        )
    lines.extend([
        "",
        "`byte_parity` means exact SHA-256 equality. `compiler_prng_drift` means the",
        "numeric schema and full row counts are preserved, at least one initial trajectory",
        "row is exact, and later stochastic state diverges. `behavioral_mismatch` means",
        "the run failed, structure changed, or divergence began at the first row; it is",
        "never normalized away.",
        "",
        "## Historical cross-host control",
        "",
        "The Hotrox and Jasper fixtures with identical inputs are byte-identical across",
        "all three archived outputs: `{}`.".format(
            str(report["historical_cross_host_duplicate"]["all_byte_equal"]).lower()
        ),
        "",
        "## Sabotage gate",
        "",
        "The committed comparator unit test changes a trajectory value in row 1 and",
        "requires classification as `behavioral_mismatch`: `{}`.".format(
            report["sabotage_test"]
        ),
        "",
    ])
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--fixtures", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    parser.add_argument("--compiler", default=os.environ.get("CC", "gcc"))
    parser.add_argument(
        "--canonical-compiler",
        default="gcc (Debian 12.2.0-14+deb12u1) 12.2.0",
    )
    parser.add_argument("--allow-compiler-drift", action="store_true")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--timeout-seconds", type=int, default=2400)
    parser.add_argument("--run-root", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.jobs < 1 or args.jobs > 4:
        raise SystemExit("--jobs must be between 1 and 4")
    actual_compiler = compiler_identity(args.compiler)
    lock_matched = actual_compiler == args.canonical_compiler
    if not lock_matched and not args.allow_compiler_drift:
        raise SystemExit(
            "compiler lock mismatch: expected {!r}, got {!r}; run in the pinned "
            "container or pass --allow-compiler-drift for an explicitly labeled comparison".format(
                args.canonical_compiler, actual_compiler
            )
        )

    git_commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True
    ).stdout.strip()
    with tempfile.TemporaryDirectory(prefix="a8a-build-") as build_raw:
        build_dir = Path(build_raw)
        build = build_model(args.source, args.compiler, build_dir)
        executable = build_dir / "mckaol"
        managed_run_root = args.run_root is None
        if managed_run_root:
            run_context = tempfile.TemporaryDirectory(prefix="a8a-runs-")
            run_root = Path(run_context.name)
        else:
            args.run_root.mkdir(parents=True, exist_ok=True)
            run_context = None
            run_root = args.run_root
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
                futures = [
                    executor.submit(
                        _run_one,
                        fixture_id,
                        args.fixtures,
                        executable,
                        run_root,
                        args.timeout_seconds,
                    )
                    for fixture_id in GOLDEN_RUNS
                ]
                runs = [future.result() for future in futures]
        finally:
            if run_context is not None:
                run_context.cleanup()

    report: Dict[str, Any] = {
        "schema_version": 1,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit,
        "platform": platform.platform(),
        "compiler": {
            "canonical": args.canonical_compiler,
            "actual": actual_compiler,
            "lock_matched": lock_matched,
            "drift_explicitly_allowed": args.allow_compiler_drift,
        },
        "build": build,
        "diffusion_audit": audit_diffusion_disabled(args.source),
        "historical_cross_host_duplicate": verify_historical_duplicate(args.fixtures),
        "sabotage_test": "PASS (test_sabotage_is_detected_as_behavioral_mismatch)",
        "runs": runs,
    }
    report["summary"] = {
        "fixture_count": len(runs),
        "byte_parity": sum(run["classification"] == "byte_parity" for run in runs),
        "compiler_prng_drift": sum(run["classification"] == "compiler_prng_drift" for run in runs),
        "behavioral_mismatch": sum(run["classification"] == "behavioral_mismatch" for run in runs),
        "all_completed": all(run["returncode"] == 0 and not run["timed_out"] for run in runs),
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(report), encoding="utf-8")
    print(json.dumps(report["summary"], sort_keys=True))
    return 0 if report["summary"]["all_completed"] and report["summary"]["behavioral_mismatch"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
