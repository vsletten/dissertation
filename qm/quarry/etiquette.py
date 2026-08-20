"""Runtime enforcement for polite quarry campaign resource use."""

from __future__ import annotations

import argparse
import atexit
import os
import sys
import time
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
)


def enforce(threads: int = 16, niceness: int = 10) -> None:
    """Cap numerical-library threads and best-effort lower process priority.

    Unset thread variables are initialized to ``threads``. Inherited values
    above the cap (or malformed values) are clamped; valid lower values are
    preserved so enforcement never inflates an existing limit.
    """
    if threads < 1:
        raise ValueError("threads must be at least 1")
    if niceness < 0:
        geteuid = getattr(os, "geteuid", None)
        if geteuid is None or geteuid() != 0:
            raise ValueError("negative niceness requires root")

    for name in THREAD_ENV_VARS:
        inherited = os.environ.get(name)
        try:
            within_cap = inherited is not None and 1 <= int(inherited) <= threads
        except ValueError:
            within_cap = False
        if not within_cap:
            os.environ[name] = str(threads)

    try:
        nice = os.nice
    except AttributeError:
        warnings.warn(
            "os.nice() unavailable; continuing without niceness", stacklevel=2
        )
        return
    try:
        nice(niceness)
    except OSError as exc:
        warnings.warn(
            f"os.nice({niceness}) failed; continuing: {exc}",
            stacklevel=2,
        )


class _Tee:
    """Minimal text stream that duplicates writes to a durable log."""

    def __init__(self, console: TextIO, log: TextIO) -> None:
        self.console = console
        self.log = log

    def write(self, text: str) -> int:
        written = self.console.write(text)
        self.log.write(text)
        return written

    def flush(self) -> None:
        self.console.flush()
        self.log.flush()

    def __getattr__(self, name: str):
        return getattr(self.console, name)


@dataclass
class EtiquetteSession:
    """Applied campaign limits and the active tee, restorable for tests."""

    threads: int
    niceness: int
    log_path: Path
    _stdout: TextIO
    _stderr: TextIO
    _log: TextIO
    _closed: bool = False

    def close(self) -> None:
        if self._closed:
            return
        if isinstance(sys.stdout, _Tee) and sys.stdout.log is self._log:
            sys.stdout = self._stdout
        if isinstance(sys.stderr, _Tee) and sys.stderr.log is self._log:
            sys.stderr = self._stderr
        self._log.flush()
        self._log.close()
        self._closed = True


def bootstrap_cli(
    script_name: str,
    *,
    default_run_root: str | Path,
    argv: list[str] | None = None,
) -> EtiquetteSession:
    """Pre-parse etiquette flags, enforce them, and tee campaign output.

    Entry-point scripts call this before importing NumPy, PySCF, or quarry
    modules that may import them. Unknown arguments are deliberately retained
    for each script's full parser.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--nice", type=int, default=10)
    parser.add_argument("--log")
    args, _ = parser.parse_known_args(argv)

    if args.threads < 1:
        parser.error("--threads must be at least 1")
    geteuid = getattr(os, "geteuid", None)
    if args.nice < 0 and (geteuid is None or geteuid() != 0):
        parser.error("negative --nice requires root")

    if args.log:
        log_path = Path(args.log).expanduser()
    else:
        stamp = time.strftime("%Y%m%dT%H%M%S")
        log_path = Path(default_run_root) / f"{script_name}-{stamp}-{os.getpid()}.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)

    original_stdout, original_stderr = sys.stdout, sys.stderr
    print(f"campaign log: {log_path}", file=original_stdout, flush=True)
    log = log_path.open("a", encoding="utf-8", buffering=1)
    sys.stdout = _Tee(original_stdout, log)
    sys.stderr = _Tee(original_stderr, log)
    session = EtiquetteSession(
        threads=args.threads,
        niceness=args.nice,
        log_path=log_path,
        _stdout=original_stdout,
        _stderr=original_stderr,
        _log=log,
    )
    atexit.register(session.close)
    enforce(threads=args.threads, niceness=args.nice)
    print(f"campaign log: {log_path}", flush=True)
    limits = " ".join(f"{name}={os.environ[name]}" for name in THREAD_ENV_VARS)
    print(f"compute etiquette: {limits}; nice increment={args.nice}", flush=True)
    return session
