#!/usr/bin/env python3
"""Replay the curated 1999 C KMC model against five golden runs.

The curated C bytes remain unchanged. Modern-header compatibility is supplied
only by forced-include compiler flags. Reports retain every mismatch and include
content-addressed source, input, toolchain, and runtime provenance.
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
import signal
import subprocess
import tempfile
import time
from typing import Any, Dict, List, Optional, Sequence, Tuple
import uuid

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


def manifest(paths: Sequence[Path], base: Path) -> Dict[str, Any]:
    rows = []
    aggregate = hashlib.sha256()
    for path in sorted(paths, key=lambda item: item.relative_to(base).as_posix()):
        relative = path.relative_to(base).as_posix()
        digest = sha256_file(path)
        rows.append({"path": relative, "bytes": path.stat().st_size, "sha256": digest})
        aggregate.update(relative.encode("utf-8") + b"\0" + digest.encode("ascii") + b"\0")
    return {"sha256": aggregate.hexdigest(), "files": rows}


def source_manifest(source: Path) -> Dict[str, Any]:
    return manifest([path for path in source.iterdir() if path.is_file()], source)


def fixture_manifest(fixtures: Path) -> Dict[str, Any]:
    paths = [fixtures / fixture / name for fixture in GOLDEN_RUNS for name in INPUT_NAMES]
    return manifest(paths, fixtures)


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
    """Compare one output; never infer the cause of a numeric divergence."""
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
    expected_hash = sha256_file(expected)
    actual_hash = sha256_file(actual)
    byte_equal = expected_hash == actual_hash
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
    elif structural_match:
        classification = "numeric_divergence"
    else:
        classification = "behavioral_mismatch"

    return {
        "expected_sha256": expected_hash,
        "actual_sha256": actual_hash,
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


def _c_function_body(text: str, name: str) -> Optional[str]:
    match = re.search(r"\b{}\s*\([^)]*\)\s*\{{".format(re.escape(name)), text)
    if not match:
        return None
    start = match.end() - 1
    depth = 0
    for index in range(start, len(text)):
        if text[index] == "{":
            depth += 1
        elif text[index] == "}":
            depth -= 1
            if depth == 0:
                return text[start + 1:index]
    return None


def audit_diffusion_disabled(source: Path) -> Dict[str, Any]:
    envrn = (source / "envrn.c").read_text(encoding="latin-1")
    evtlist_path = source / "evtlist.c"
    evtlist = evtlist_path.read_text(encoding="latin-1") if evtlist_path.exists() else ""
    rxnlist = (source / "rxnlist.h").read_text(encoding="latin-1")
    ndes_match = re.search(r"#define\s+NDES\s+(\d+)", rxnlist)
    nrxn_match = re.search(r"#define\s+NRXN\s+(\d+)", rxnlist)
    body = _c_function_body(envrn, "isActive")
    diffusion_false = bool(body and re.search(
        r"else\s*\{\s*/\*\s*diffusion\s*\*/\s*result\s*=\s*FALSE\s*;\s*\}",
        body,
        re.DOTALL,
    ))
    scheduler_uses_guard = bool(re.search(r"isActive\s*\(\s*s\s*,\s*l\s*,\s*i\s*\)", evtlist))
    if not ndes_match or not nrxn_match:
        return {"status": "behavioral_mismatch", "reason": "NDES/NRXN definitions missing"}
    ndes = int(ndes_match.group(1))
    nrxn = int(nrxn_match.group(1))
    pinned = (ndes, nrxn) == (24, 28) and diffusion_false and scheduler_uses_guard
    return {
        "status": "pinned_disabled" if pinned else "behavioral_mismatch",
        "ndes": ndes,
        "nrxn": nrxn,
        "diffusion_ids": list(range(ndes, nrxn)),
        "isActive_final_branch_returns_false": diffusion_false,
        "new_evtList_calls_isActive": scheduler_uses_guard,
    }


def conformance_gate_passes(
    runs: Sequence[Dict[str, Any]],
    diffusion: Dict[str, Any],
    duplicate: Dict[str, Any],
    sabotage: Dict[str, Any],
    source_matches: bool,
    fixtures_match: bool,
    setup_error: Optional[str],
) -> bool:
    """Canonical success is exact parity; drift candidates remain non-passing."""
    return (
        setup_error is None
        and len(runs) == len(GOLDEN_RUNS)
        and [run.get("fixture") for run in runs] == list(GOLDEN_RUNS)
        and all(run.get("returncode") == 0 and not run.get("timed_out") for run in runs)
        and all(run.get("classification") == "byte_parity" for run in runs)
        and diffusion.get("status") == "pinned_disabled"
        and bool(duplicate.get("all_byte_equal"))
        and bool(sabotage.get("passed"))
        and source_matches
        and fixtures_match
    )


def run_sabotage_gate() -> Dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="a8a-sabotage-") as raw:
        root = Path(raw)
        rows = ["{:.1f},{},{}".format(index / 10.0, index + 10, index + 20) for index in range(20)]
        expected = root / "expected.dat"
        actual = root / "actual.dat"
        expected.write_text("\n".join(rows) + "\n", encoding="ascii")
        sabotaged = list(rows)
        sabotaged[11] = "1.1,22,32"
        actual.write_text("\n".join(sabotaged) + "\n", encoding="ascii")
        comparison = compare_numeric_output(expected, actual, 3)
    passed = (
        comparison["classification"] == "numeric_divergence"
        and comparison["first_mismatch_row"] == 12
        and not comparison["byte_equal"]
    )
    synthetic_runs = [
        {
            "fixture": fixture,
            "classification": "compiler_prng_drift_candidate",
            "returncode": 0,
            "timed_out": False,
        }
        for fixture in GOLDEN_RUNS
    ]
    synthetic_gate = conformance_gate_passes(
        runs=synthetic_runs,
        diffusion={"status": "pinned_disabled"},
        duplicate={"all_byte_equal": True},
        sabotage={"passed": passed},
        source_matches=True,
        fixtures_match=True,
        setup_error=None,
    )
    return {
        "passed": passed,
        "perturbed_row": 12,
        "observed_classification": comparison["classification"],
        "first_mismatch_row": comparison["first_mismatch_row"],
        "conformance_gate_passed": synthetic_gate,
    }


def _run_process(command: Sequence[str], cwd: Optional[Path], timeout: int) -> Dict[str, Any]:
    started = time.monotonic()
    kwargs: Dict[str, Any] = {
        "cwd": str(cwd) if cwd else None,
        "stdout": subprocess.PIPE,
        "stderr": subprocess.PIPE,
        "text": True,
    }
    if os.name == "posix":
        kwargs["start_new_session"] = True
    process = subprocess.Popen(list(command), **kwargs)
    timed_out = False
    try:
        stdout, stderr = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        if os.name == "posix":
            os.killpg(process.pid, signal.SIGKILL)
        else:
            process.kill()
        stdout, stderr = process.communicate()
    return {
        "command": list(command),
        "returncode": process.returncode,
        "timed_out": timed_out,
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "stdout": stdout,
        "stderr": stderr,
    }


def compiler_provenance(compiler: str, timeout: int) -> Dict[str, str]:
    version = _run_process([compiler, "--version"], None, timeout)
    target = _run_process([compiler, "-dumpmachine"], None, timeout)
    if version["timed_out"] or version["returncode"] != 0:
        raise RuntimeError("compiler identity failed: {}".format(version))
    if target["timed_out"] or target["returncode"] != 0:
        raise RuntimeError("compiler target query failed: {}".format(target))
    return {
        "identity": version["stdout"].splitlines()[0],
        "target": target["stdout"].strip(),
    }


def build_model(source: Path, compiler: str, output_dir: Path, timeout: int) -> Dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    executable = output_dir / "mckaol"
    command = [compiler, *COMPATIBILITY_FLAGS, "-o", str(executable)]
    command.extend(str(source / name) for name in SOURCE_NAMES)
    command.append("-lm")
    record = _run_process(command, None, timeout)
    record["compatibility_flags"] = list(COMPATIBILITY_FLAGS)
    if record["timed_out"] or record["returncode"] != 0:
        raise RuntimeError("compiler failed: {}".format(json.dumps(record, indent=2)))
    record["executable_sha256"] = sha256_file(executable)
    return record


def _fixture_metadata(fixture: Path) -> Dict[str, Any]:
    fields = []
    for raw in (fixture / "data.sim").read_text(encoding="latin-1").splitlines():
        fields.append(raw.split("#", 1)[0].strip())
    return {
        "steps": int(fields[0]),
        "write_every_steps": int(fields[1]),
        "movie_every_steps": int(fields[2]),
        "seed": int(fields[3]),
        "draw_bonds": int(fields[4]),
        "inputs": {name: sha256_file(fixture / name) for name in INPUT_NAMES},
    }


def _failed_run(fixture_id: str, error: str) -> Dict[str, Any]:
    return {
        "fixture": fixture_id,
        "returncode": None,
        "timed_out": False,
        "elapsed_seconds": 0.0,
        "classification": "behavioral_mismatch",
        "error": error,
        "outputs": {},
    }


def _run_one(
    fixture_id: str,
    fixtures: Path,
    executable: Path,
    execution_root: Path,
    timeout_seconds: int,
    drift_candidate_allowed: bool,
) -> Dict[str, Any]:
    fixture = fixtures / fixture_id
    run_dir = execution_root / fixture_id.replace("/", "-")
    try:
        run_dir.mkdir(parents=True, exist_ok=False)
        for name in INPUT_NAMES:
            shutil.copy2(str(fixture / name), str(run_dir / name))
        shutil.copy2(str(executable), str(run_dir / "mckaol"))
        process = _run_process([str(run_dir / "mckaol")], run_dir, timeout_seconds)
        outputs = {
            name: compare_numeric_output(fixture / name, run_dir / name, columns)
            for name, columns in OUTPUT_COLUMNS.items()
        }
        all_structural = all(item.get("structural_match", False) for item in outputs.values())
        all_byte_equal = all(item.get("byte_equal", False) for item in outputs.values())
        results_prefix = outputs["results.dat"].get("matching_prefix_rows", 0)
        if process["returncode"] != 0 or process["timed_out"] or not all_structural:
            classification = "behavioral_mismatch"
        elif all_byte_equal:
            classification = "byte_parity"
        elif drift_candidate_allowed and results_prefix >= 10:
            classification = "compiler_prng_drift_candidate"
        else:
            classification = "behavioral_mismatch"
        return {
            "fixture": fixture_id,
            "fixture_metadata": _fixture_metadata(fixture),
            "returncode": process["returncode"],
            "timed_out": process["timed_out"],
            "elapsed_seconds": process["elapsed_seconds"],
            "classification": classification,
            "drift_candidate_basis": (
                "content manifests match; compiler/architecture differs from canonical; "
                "all schemas and row counts match; >=10 initial results rows are exact"
                if classification == "compiler_prng_drift_candidate" else None
            ),
            "stdout": process["stdout"],
            "stderr": process["stderr"],
            "outputs": outputs,
        }
    except Exception as exc:
        return _failed_run(fixture_id, "{}: {}".format(type(exc).__name__, exc))


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


def _git_commit() -> Optional[str]:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"], check=True, capture_output=True,
            text=True, timeout=10,
        )
        return result.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None


def render_markdown(report: Dict[str, Any]) -> str:
    compiler = report.get("compiler", {})
    build = report.get("build", {})
    lines = [
        "# 1999 C KMC conformance report", "",
        "Generated: `{}`".format(report["generated_at"]), "",
        "## Provenance", "",
        "- Git commit: `{}`".format(report.get("git_commit") or "unavailable; content manifests govern"),
        "- Source manifest SHA-256: `{}`".format(report["source_manifest"]["sha256"]),
        "- Fixture-input manifest SHA-256: `{}`".format(report["fixture_input_manifest"]["sha256"]),
        "- Platform: `{}`".format(report["runtime"]["platform"]),
        "- Architecture: `{}`".format(report["runtime"]["architecture"]),
        "- libc: `{}`".format(report["runtime"]["libc"]),
        "- Compiler: `{}`".format(compiler.get("actual", "unavailable")),
        "- Compiler target: `{}`".format(compiler.get("target", "unavailable")),
        "- Canonical container: `{}`".format(report["toolchain_lock"].get("container_image")),
        "- Compiler lock matched: `{}`".format(str(compiler.get("lock_matched", False)).lower()),
        "- Compatibility flags: `{}`".format(" ".join(build.get("compatibility_flags", []))),
        "- Diffusion IDs 24–27: `{}`".format(report["diffusion_audit"]["status"]),
        "", "## Fixture outcomes", "",
        "| fixture | outcome | results parity | first divergence | surfAl rows | surfSi rows | seconds |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for run in report.get("runs", []):
        outputs = run.get("outputs", {})
        results = outputs.get("results.dat", {})
        al = outputs.get("surfAl.out", {})
        si = outputs.get("surfSi.out", {})
        lines.append(
            "| {fixture} | {classification} | {parity} | {first} | {al} | {si} | {elapsed} |".format(
                fixture=run["fixture"], classification=run["classification"],
                parity="byte" if results.get("byte_equal") else "diverged",
                first=results.get("first_mismatch_row") or "—",
                al=al.get("actual_rows", "—"), si=si.get("actual_rows", "—"),
                elapsed=run.get("elapsed_seconds", 0),
            )
        )
    lines.extend([
        "", "`byte_parity` is exact SHA-256 equality. A",
        "`compiler_prng_drift_candidate` is explicitly non-canonical evidence: content",
        "manifests match, the compiler/architecture differs, every schema and row count",
        "matches, and at least ten initial trajectory rows are exact. It remains a",
        "candidate rather than a normalized fact and never passes the canonical gate.",
        "Anything else is a",
        "`behavioral_mismatch` and fails the gate.", "",
        "## Independent controls", "",
        "- Historical identical-input Hotrox/Jasper outputs byte-equal: `{}`.".format(
            str(report["historical_cross_host_duplicate"]["all_byte_equal"]).lower()),
        "- Late-row sabotage detected without being labeled drift: `{}`.".format(
            str(report["sabotage_test"]["passed"]).lower()), "",
    ])
    if report.get("setup_error"):
        lines.extend(["## Setup failure", "", "```text", report["setup_error"], "```", ""])
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--fixtures", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--markdown", type=Path)
    parser.add_argument("--toolchain-lock", type=Path)
    parser.add_argument("--compiler", default=os.environ.get("CC", "gcc"))
    parser.add_argument("--allow-compiler-drift", action="store_true")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--timeout-seconds", type=int, default=2400)
    parser.add_argument("--build-timeout-seconds", type=int, default=120)
    parser.add_argument("--run-root", type=Path)
    parser.add_argument("--print-manifests", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.jobs < 1 or args.jobs > 4:
        raise SystemExit("--jobs must be between 1 and 4")
    if args.print_manifests:
        source_info = source_manifest(args.source)
        fixture_info = fixture_manifest(args.fixtures)
        print(json.dumps({
            "source_manifest_sha256": source_info["sha256"],
            "fixture_inputs_manifest_sha256": fixture_info["sha256"],
        }, sort_keys=True))
        return 0
    if not args.report:
        raise SystemExit("--report is required unless --print-manifests is used")

    runtime = {
        "platform": platform.platform(),
        "architecture": platform.machine(),
        "libc": " ".join(part for part in platform.libc_ver() if part) or "unknown",
    }
    report: Dict[str, Any] = {
        "schema_version": 2,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "git_commit": _git_commit(),
        "source_manifest": {"sha256": None, "files": []},
        "fixture_input_manifest": {"sha256": None, "files": []},
        "content_lock": {"source_matches": False, "fixture_inputs_match": False},
        "toolchain_lock": {},
        "runtime": runtime,
        "diffusion_audit": {"status": "behavioral_mismatch", "reason": "not run"},
        "historical_cross_host_duplicate": {"all_byte_equal": False, "reason": "not run"},
        "sabotage_test": {"passed": False, "reason": "not run"},
        "runs": [],
    }
    setup_error: Optional[str] = None
    source_matches = False
    fixtures_match = False
    try:
        lock_path = args.toolchain_lock or (args.source.parent / "conformance-toolchain.json")
        lock = json.loads(lock_path.read_text(encoding="utf-8"))
        report["toolchain_lock"] = lock
        source_info = source_manifest(args.source)
        fixture_info = fixture_manifest(args.fixtures)
        report["source_manifest"] = source_info
        report["fixture_input_manifest"] = fixture_info
        source_matches = source_info["sha256"] == lock.get("source_manifest_sha256")
        fixtures_match = fixture_info["sha256"] == lock.get("fixture_inputs_manifest_sha256")
        report["content_lock"] = {
            "source_matches": source_matches,
            "fixture_inputs_match": fixtures_match,
        }
        diffusion = audit_diffusion_disabled(args.source)
        duplicate = verify_historical_duplicate(args.fixtures)
        sabotage = run_sabotage_gate()
        report["diffusion_audit"] = diffusion
        report["historical_cross_host_duplicate"] = duplicate
        report["sabotage_test"] = sabotage
        if not source_matches or not fixtures_match:
            raise RuntimeError("source or fixture-input content manifest does not match toolchain lock")
        if diffusion.get("status") != "pinned_disabled":
            raise RuntimeError("historical diffusion-disable audit failed")
        if not duplicate.get("all_byte_equal"):
            raise RuntimeError("historical cross-host duplicate control failed")
        if not sabotage.get("passed") or sabotage.get("conformance_gate_passed"):
            raise RuntimeError("end-to-end sabotage gate failed")

        compiler = compiler_provenance(args.compiler, args.build_timeout_seconds)
        compiler_matches = compiler["identity"] == lock.get("compiler_identity")
        architecture_matches = runtime["architecture"] == lock.get("canonical_architecture")
        report["compiler"] = {
            "actual": compiler["identity"],
            "target": compiler["target"],
            "canonical": lock.get("compiler_identity"),
            "lock_matched": compiler_matches,
            "architecture_lock_matched": architecture_matches,
            "drift_explicitly_allowed": args.allow_compiler_drift,
        }
        if (not compiler_matches or not architecture_matches) and not args.allow_compiler_drift:
            raise RuntimeError("canonical compiler/architecture lock mismatch")
        drift_candidate_allowed = args.allow_compiler_drift and (not compiler_matches or not architecture_matches)
        with tempfile.TemporaryDirectory(prefix="a8a-build-") as build_raw:
            build_dir = Path(build_raw)
            report["build"] = build_model(
                args.source, args.compiler, build_dir, args.build_timeout_seconds
            )
            base_root = args.run_root or Path(tempfile.mkdtemp(prefix="a8a-runs-"))
            base_root.mkdir(parents=True, exist_ok=True)
            execution_root = base_root / ("execution-" + uuid.uuid4().hex)
            execution_root.mkdir(parents=True, exist_ok=False)
            report["execution_root"] = str(execution_root)
            future_to_fixture: Dict[concurrent.futures.Future, str] = {}
            with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
                for fixture_id in GOLDEN_RUNS:
                    future = executor.submit(
                        _run_one,
                        fixture_id,
                        args.fixtures,
                        build_dir / "mckaol",
                        execution_root,
                        args.timeout_seconds,
                        drift_candidate_allowed,
                    )
                    future_to_fixture[future] = fixture_id
                run_map: Dict[str, Dict[str, Any]] = {}
                for future, fixture_id in future_to_fixture.items():
                    try:
                        run_map[fixture_id] = future.result()
                    except Exception as exc:
                        run_map[fixture_id] = _failed_run(
                            fixture_id, "future {}: {}".format(type(exc).__name__, exc)
                        )
            report["runs"] = [run_map[fixture] for fixture in GOLDEN_RUNS]
    except Exception as exc:
        setup_error = "{}: {}".format(type(exc).__name__, exc)
        report["setup_error"] = setup_error

    if not report["runs"]:
        failure = setup_error or "run stage produced no outcomes"
        report["runs"] = [_failed_run(fixture, failure) for fixture in GOLDEN_RUNS]
    classifications = [run["classification"] for run in report["runs"]]
    gate_passed = conformance_gate_passes(
        runs=report["runs"],
        diffusion=report["diffusion_audit"],
        duplicate=report["historical_cross_host_duplicate"],
        sabotage=report["sabotage_test"],
        source_matches=source_matches,
        fixtures_match=fixtures_match,
        setup_error=setup_error,
    )
    report["summary"] = {
        "fixture_count": len(report["runs"]),
        "byte_parity": classifications.count("byte_parity"),
        "compiler_prng_drift_candidate": classifications.count("compiler_prng_drift_candidate"),
        "behavioral_mismatch": classifications.count("behavioral_mismatch"),
        "all_completed": all(
            run.get("returncode") == 0 and not run.get("timed_out")
            for run in report["runs"]
        ),
        "gate_passed": gate_passed,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(report), encoding="utf-8")
    print(json.dumps(report["summary"], sort_keys=True))
    return 0 if gate_passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
