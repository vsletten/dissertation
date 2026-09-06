#!/usr/bin/env python3
"""Bounded, receipt-backed supervisor for one Phase-2 ladder family.

The supervisor waits for the shared GPU lease, runs explicitly named connectivity
cells serially through ``phase2_ladder.py``, and writes an atomic terminal receipt
on every exit path.  It never treats a zero child exit as sufficient: each cell's
identity, finite scientific summary, and provenance store are checked before the
cell enters the completed ledger.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import signal
import subprocess
import sys
import time
from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    _ETIQUETTE = bootstrap_cli(
        "a3-family-campaign",
        default_run_root=Path.cwd() / "qm/runs",
    )

from quarry import etiquette

FAMILY_MAX_CONNECTIVITY = {"oss": 4, "osa": 4, "oaa": 6}
SCHEMA = "a3-family-campaign-terminal-v1"
PROGRESS_SCHEMA = "a3-family-campaign-progress-v1"


class StopRequested(RuntimeError):
    """Raised after SIGINT/SIGTERM so the terminal receipt is still written."""


def now() -> str:
    return datetime.now(UTC).isoformat().replace("+00:00", "Z")


def atomic_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    temporary.replace(path)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_output(worktree: Path, *args: str) -> str:
    return subprocess.run(
        ["git", *args],
        cwd=worktree,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def verify_revision(worktree: Path, expected_sha: str) -> str:
    branch = git_output(worktree, "rev-parse", "--abbrev-ref", "HEAD")
    if branch == "HEAD":
        raise RuntimeError("campaign worktree is detached")
    head = git_output(worktree, "rev-parse", "HEAD")
    remote = git_output(worktree, "rev-parse", f"origin/{branch}")
    dirty = git_output(worktree, "status", "--porcelain")
    if head != expected_sha or remote != expected_sha:
        raise RuntimeError(
            "source revision drifted: "
            f"expected={expected_sha} head={head} remote={remote}"
        )
    if dirty:
        raise RuntimeError(f"campaign worktree is dirty: {dirty}")
    return branch


def wait_for_gpu(
    *,
    timeout_seconds: float,
    poll_seconds: float,
    heartbeat: Callable[[str], None],
    lane_probe: Callable[[], tuple[bool, str]] = etiquette.gpu_lane_available,
    monotonic: Callable[[], float] = time.monotonic,
    sleep: Callable[[float], None] = time.sleep,
) -> str:
    deadline = monotonic() + timeout_seconds
    while True:
        try:
            available, message = lane_probe()
        except etiquette.GpuLeaseBusy as exc:
            available, message = False, str(exc)
        heartbeat(message)
        if available:
            return message
        remaining = deadline - monotonic()
        if remaining <= 0:
            raise TimeoutError(f"GPU wait exhausted: {message}")
        sleep(min(poll_seconds, remaining))


def validate_result(
    result_path: Path,
    *,
    family: str,
    state: str,
    n_intact: int,
) -> dict[str, object]:
    store_path = result_path.with_name("store.sqlite")
    if not result_path.is_file() or not store_path.is_file():
        raise RuntimeError(
            f"n={n_intact} exited zero without results.json + store.sqlite"
        )
    payload = json.loads(result_path.read_text())
    expected = {"family": family, "state": state, "n_intact": n_intact}
    observed = {key: payload.get(key) for key in expected}
    if observed != expected:
        raise RuntimeError(
            f"n={n_intact} result identity mismatch: "
            f"expected={expected} observed={observed}"
        )
    required_finite = ("dG_kj", "dH_kj", "ts_imaginary_cm")
    for key in required_finite:
        value = payload.get(key)
        if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
            raise RuntimeError(f"n={n_intact} result has invalid {key}={value!r}")
    route = payload.get("route")
    if not isinstance(route, str) or not route:
        raise RuntimeError(f"n={n_intact} result has no route provenance")
    return {
        "n_intact": n_intact,
        "result_path": str(result_path),
        "result_sha256": sha256_path(result_path),
        "store_path": str(store_path),
        "store_sha256": sha256_path(store_path),
        "dG_kj": float(payload["dG_kj"]),
        "dH_kj": float(payload["dH_kj"]),
        "ts_imaginary_cm": float(payload["ts_imaginary_cm"]),
        "route": route,
    }


def parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--worktree", type=Path, required=True)
    ap.add_argument("--run-root", type=Path, required=True)
    ap.add_argument("--expected-git-sha", required=True)
    ap.add_argument("--family", choices=sorted(FAMILY_MAX_CONNECTIVITY), required=True)
    ap.add_argument("--state", choices=("acid", "neutral"), required=True)
    ap.add_argument("--cells", type=int, nargs="+", required=True)
    ap.add_argument("--wait-for-gpu-seconds", type=float, default=36 * 3600)
    ap.add_argument("--gpu-poll-seconds", type=float, default=60)
    ap.add_argument("--gpu-mem-gb", type=float, default=16)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--nice", type=int, default=10)
    ap.add_argument("--log", type=Path, help="supervisor tee path (owned by etiquette)")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    cells = tuple(args.cells)
    if tuple(sorted(set(cells))) != cells:
        raise SystemExit("--cells must be unique and strictly increasing")
    maximum = FAMILY_MAX_CONNECTIVITY[args.family]
    if not cells or cells[0] < 1 or cells[-1] > maximum:
        raise SystemExit(f"--cells must be within 1..{maximum} for {args.family}")
    if args.threads < 1 or args.threads > 16:
        raise SystemExit("--threads must be within 1..16")
    if args.wait_for_gpu_seconds < 0 or args.gpu_poll_seconds <= 0:
        raise SystemExit("GPU wait/poll bounds must be non-negative/positive")

    root = args.run_root.resolve()
    progress_path = root / "family-progress.json"
    terminal_path = root / "terminal-receipt.json"
    if terminal_path.exists():
        raise SystemExit(f"refusing to overwrite terminal receipt: {terminal_path}")

    started = now()
    started_monotonic = time.monotonic()
    completed: list[dict[str, object]] = []
    current_cell: int | None = None
    terminal_reason = "runner-error"
    exit_code = 1
    error: str | None = None

    def heartbeat(status: str) -> None:
        atomic_json(
            progress_path,
            {
                "schema": PROGRESS_SCHEMA,
                "started_at": started,
                "updated_at": now(),
                "expected_git_sha": args.expected_git_sha,
                "family": args.family,
                "state": args.state,
                "requested_n_intact": list(cells),
                "current_n_intact": current_cell,
                "completed": completed,
                "status": status,
            },
        )

    def stop(signum: int, _frame: object) -> None:
        raise StopRequested(signal.Signals(signum).name)

    previous_handlers = {
        sig: signal.signal(sig, stop) for sig in (signal.SIGTERM, signal.SIGINT)
    }
    try:
        branch = verify_revision(args.worktree, args.expected_git_sha)
        root.mkdir(parents=True, exist_ok=True)
        (root / "logs").mkdir(exist_ok=True)
        (root / "runs").mkdir(exist_ok=True)
        lane_message = wait_for_gpu(
            timeout_seconds=args.wait_for_gpu_seconds,
            poll_seconds=args.gpu_poll_seconds,
            heartbeat=lambda message: heartbeat(f"waiting-for-gpu: {message}"),
        )
        heartbeat(f"gpu-available: {lane_message}")
        driver = args.worktree / "qm/scripts/phase2_ladder.py"
        for n_intact in cells:
            current_cell = n_intact
            verify_revision(args.worktree, args.expected_git_sha)
            heartbeat("running")
            command = [
                sys.executable,
                str(driver),
                "--family",
                args.family,
                "--state",
                args.state,
                "--n-intact",
                str(n_intact),
                "--gpu",
                "--gpu-mem-gb",
                str(args.gpu_mem_gb),
                "--threads",
                str(args.threads),
                "--nice",
                str(args.nice),
                "--run-root",
                str(root / "runs"),
                "--log",
                str(root / "logs" / f"{args.family}-{args.state}-n{n_intact}.log"),
            ]
            cell_started = time.monotonic()
            result = subprocess.run(command, cwd=args.worktree, check=False)
            if result.returncode != 0:
                terminal_reason = f"cell-n{n_intact}-failed"
                exit_code = result.returncode
                raise RuntimeError(
                    f"phase2 ladder cell n={n_intact} exited {result.returncode}"
                )
            result_path = (
                root
                / "runs/phase2"
                / f"{args.family}-{args.state}-n{n_intact}-s2-b3lyp-def2-svp"
                / "results.json"
            )
            record = validate_result(
                result_path,
                family=args.family,
                state=args.state,
                n_intact=n_intact,
            )
            record["elapsed_seconds"] = time.monotonic() - cell_started
            completed.append(record)
            heartbeat("running")
        terminal_reason = "family-complete"
        exit_code = 0
    except BaseException as exc:
        error = f"{type(exc).__name__}: {exc}"
        if isinstance(exc, StopRequested):
            terminal_reason = f"signal-{exc}"
            exit_code = 128 + getattr(signal, f"SIG{exc}")
        elif isinstance(exc, TimeoutError):
            terminal_reason = "gpu-wait-timeout"
            exit_code = 75
        elif terminal_reason == "runner-error":
            exit_code = 1
    finally:
        for sig, previous in previous_handlers.items():
            signal.signal(sig, previous)
        receipt: dict[str, object] = {
            "schema": SCHEMA,
            "started_at": started,
            "completed_at": now(),
            "elapsed_seconds": time.monotonic() - started_monotonic,
            "expected_git_sha": args.expected_git_sha,
            "observed_git_sha": git_output(args.worktree, "rev-parse", "HEAD"),
            "family": args.family,
            "state": args.state,
            "requested_n_intact": list(cells),
            "current_n_intact": current_cell,
            "completed": completed,
            "success": exit_code == 0 and len(completed) == len(cells),
            "terminal_reason": terminal_reason,
            "exit_code": exit_code,
            "error": error,
            "worktree": str(args.worktree),
            "branch": locals().get("branch"),
        }
        atomic_json(terminal_path, receipt)
        heartbeat("complete" if receipt["success"] else "failed")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
