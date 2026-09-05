"""Runtime enforcement for polite quarry campaign resource use."""

from __future__ import annotations

import argparse
import atexit
import fcntl
import json
import math
import os
import signal
import stat
import sys
import time
import warnings
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path
from typing import TextIO

THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
)
GPU_TOTAL_GB = 24.0
GPU_OLLAMA_FLOOR_GB = 6.0
DEFAULT_GPU_MEM_GB = 16.0
DEFAULT_GPU_TTL_HOURS = 24.0


class GpuLeaseBusy(RuntimeError):
    """Raised when another live process owns the shared GPU lane."""


def _utc_now() -> datetime:
    return datetime.now(UTC)


def _iso_utc(value: datetime) -> str:
    return value.astimezone(UTC).isoformat().replace("+00:00", "Z")


def _default_lease_path() -> Path:
    configured = os.environ.get("GPU_LEASE_PATH")
    if configured:
        return Path(configured).expanduser()
    return Path.home() / ".local" / "state" / "gpu-lease" / "lease.json"


def _pid_is_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _validate_gpu_mem_gb(value: float) -> None:
    maximum = GPU_TOTAL_GB - GPU_OLLAMA_FLOOR_GB
    if not math.isfinite(value) or value <= 0 or value > maximum:
        raise ValueError(
            f"GPU memory limit must be > 0 and <= {maximum:g} GB "
            f"to preserve the {GPU_OLLAMA_FLOOR_GB:g} GB ollama floor"
        )


def _open_regular(path: Path, flags: int, mode: int = 0o600) -> int:
    """Open one state file without following links or accepting special files."""
    nofollow = getattr(os, "O_NOFOLLOW", 0)
    cloexec = getattr(os, "O_CLOEXEC", 0)
    fd = os.open(path, flags | nofollow | cloexec, mode)
    if not stat.S_ISREG(os.fstat(fd).st_mode):
        os.close(fd)
        raise OSError(f"GPU lease state path is not a regular file: {path}")
    return fd


def _read_lease(path: Path) -> dict[str, object]:
    try:
        fd = _open_regular(path, os.O_RDONLY)
        with os.fdopen(fd, encoding="utf-8") as stream:
            payload = json.load(stream)
    except (OSError, json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise GpuLeaseBusy(
            f"GPU lane busy (unreadable lease at {path}: {exc})"
        ) from exc
    required = {"owner", "pid", "started", "ttl", "expected_gb"}
    if not isinstance(payload, dict) or not required.issubset(payload):
        raise GpuLeaseBusy(f"GPU lane busy (invalid lease at {path})")
    if not isinstance(payload["owner"], str) or not payload["owner"]:
        raise GpuLeaseBusy(f"GPU lane busy (invalid owner in {path})")
    if not isinstance(payload["pid"], int) or payload["pid"] <= 0:
        raise GpuLeaseBusy(f"GPU lane busy (invalid pid in {path})")
    return payload


@contextmanager
def _lease_lock(lease_path: Path) -> Iterator[None]:
    lease_path.parent.mkdir(parents=True, exist_ok=True)
    lock_path = lease_path.with_name(f"{lease_path.name}.lock")
    try:
        fd = _open_regular(lock_path, os.O_CREAT | os.O_RDWR)
    except OSError as exc:
        raise GpuLeaseBusy(
            f"GPU lane busy (unsafe lock path {lock_path}: {exc})"
        ) from exc
    try:
        os.fchmod(fd, 0o600)
        fcntl.flock(fd, fcntl.LOCK_EX)
        yield
    finally:
        fcntl.flock(fd, fcntl.LOCK_UN)
        os.close(fd)


def _write_lease(path: Path, payload: dict[str, object]) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        try:
            fd = _open_regular(temporary, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            try:
                stale = _open_regular(temporary, os.O_RDONLY)
                os.close(stale)
                temporary.unlink()
                fd = _open_regular(temporary, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            except OSError as exc:
                raise GpuLeaseBusy(
                    "GPU lane busy (unrecoverable leftover lease temp "
                    f"{temporary}: {exc})"
                ) from exc
        except OSError as exc:
            raise GpuLeaseBusy(
                f"GPU lane busy (unsafe temporary lease path {temporary}: {exc})"
            ) from exc
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            json.dump(payload, stream, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
            os.fchmod(stream.fileno(), 0o600)
        try:
            os.replace(temporary, path)
        except OSError as exc:
            raise GpuLeaseBusy(
                f"GPU lane busy (failed to install lease at {path}: {exc})"
            ) from exc
    finally:
        temporary.unlink(missing_ok=True)


def _record_stale_break(
    lease_path: Path,
    payload: dict[str, object],
    *,
    now: datetime,
) -> None:
    note_path = lease_path.with_name("stale-breaks.log")
    line = (
        f"{_iso_utc(now)} broke stale GPU lease: "
        f"owner={payload['owner']} pid={payload['pid']} reason=pid-dead\n"
    )
    try:
        fd = _open_regular(note_path, os.O_CREAT | os.O_APPEND | os.O_WRONLY)
    except OSError as exc:
        raise GpuLeaseBusy(
            f"GPU lane busy (unsafe stale-break path {note_path}: {exc})"
        ) from exc
    with os.fdopen(fd, "a", encoding="utf-8") as stream:
        os.fchmod(stream.fileno(), 0o600)
        stream.write(line)
        stream.flush()
        os.fsync(stream.fileno())


@dataclass
class GpuLease:
    """One acquired process-bound GPU lease."""

    owner: str
    pid: int
    started: str
    ttl_hours: float
    expected_gb: float
    path: Path
    _acquirer_pid: int = field(default_factory=os.getpid)
    _released: bool = False
    _releasing: bool = False
    _signal_handlers: dict[int, object] = field(default_factory=dict)

    def install_signal_handlers(self) -> None:
        """Release before the default SIGTERM/SIGINT process exit."""
        for signum in (signal.SIGTERM, signal.SIGINT):
            try:
                previous = signal.getsignal(signum)
                signal.signal(signum, self._handle_signal)
            except ValueError:
                # Signal handlers can only be installed from the main thread.
                continue
            self._signal_handlers[signum] = previous

    def _handle_signal(self, signum: int, frame: object) -> None:
        previous = self._signal_handlers.get(signum, signal.SIG_DFL)
        self.release()
        if previous is signal.SIG_IGN:
            return
        if callable(previous):
            previous(signum, frame)
            return
        signal.signal(signum, signal.SIG_DFL)
        os.kill(os.getpid(), signum)

    def release(self) -> None:
        """Remove only the exact lease record this process acquired."""
        if os.getpid() != self._acquirer_pid or self._released or self._releasing:
            return
        self._releasing = True
        previous_mask = None
        pthread_sigmask = getattr(signal, "pthread_sigmask", None)
        if pthread_sigmask is not None:
            previous_mask = pthread_sigmask(
                signal.SIG_BLOCK, {signal.SIGTERM, signal.SIGINT}
            )
        try:
            with _lease_lock(self.path):
                if self.path.exists():
                    try:
                        payload = _read_lease(self.path)
                    except GpuLeaseBusy as exc:
                        warnings.warn(
                            f"refusing to remove changed GPU lease: {exc}", stacklevel=2
                        )
                    else:
                        if (
                            payload.get("owner") == self.owner
                            and payload.get("pid") == self.pid
                            and payload.get("started") == self.started
                        ):
                            self.path.unlink()
        finally:
            for signum, previous in self._signal_handlers.items():
                if signal.getsignal(signum) == self._handle_signal:
                    signal.signal(signum, previous)
            self._signal_handlers.clear()
            self._released = True
            self._releasing = False
            if previous_mask is not None:
                assert pthread_sigmask is not None
                pthread_sigmask(signal.SIG_SETMASK, previous_mask)

    def __enter__(self) -> GpuLease:
        return self

    def __exit__(self, exc_type: object, exc: object, traceback: object) -> None:
        self.release()


def acquire_gpu(
    owner: str,
    expected_gb: float,
    ttl_hours: float,
    *,
    lease_path: Path | None = None,
    pid: int | None = None,
    now: datetime | None = None,
    pid_alive: Callable[[int], bool] = _pid_is_alive,
    install_signal_handlers: bool = True,
) -> GpuLease:
    """Atomically acquire the shared GPU lane or raise :class:`GpuLeaseBusy`.

    A lease is broken only when its recorded PID is deterministically dead. A
    live process remains authoritative even beyond its declared TTL; the TTL is
    recorded for operator diagnosis, never used to steal a running campaign.
    """
    if not owner.strip():
        raise ValueError("GPU lease owner must be non-empty")
    if not math.isfinite(ttl_hours) or ttl_hours <= 0:
        raise ValueError("GPU lease ttl_hours must be finite and positive")
    _validate_gpu_mem_gb(expected_gb)
    path = lease_path or _default_lease_path()
    process_id = pid or os.getpid()
    started_at = _iso_utc(now or _utc_now())
    payload: dict[str, object] = {
        "owner": owner,
        "pid": process_id,
        "started": started_at,
        "ttl": ttl_hours,
        "expected_gb": expected_gb,
    }

    with _lease_lock(path):
        if path.exists():
            current = _read_lease(path)
            current_pid = int(current["pid"])
            if pid_alive(current_pid):
                raise GpuLeaseBusy(
                    "GPU lane busy "
                    f"(owner={current['owner']} pid={current_pid} "
                    f"expected_gb={current['expected_gb']})"
                )
            _record_stale_break(path, current, now=now or _utc_now())
            path.unlink()
        _write_lease(path, payload)

    lease = GpuLease(
        owner=owner,
        pid=process_id,
        started=started_at,
        ttl_hours=ttl_hours,
        expected_gb=expected_gb,
        path=path,
    )
    if install_signal_handlers:
        lease.install_signal_handlers()
    atexit.register(lease.release)
    return lease


def gpu_lane_available(
    *,
    lease_path: Path | None = None,
    pid_alive: Callable[[int], bool] = _pid_is_alive,
    now: datetime | None = None,
) -> tuple[bool, str]:
    """Return a preflight verdict, cleaning only deterministically dead leases."""
    path = lease_path or _default_lease_path()
    with _lease_lock(path):
        if not path.exists():
            return True, "GPU lane free (no lease)"
        current = _read_lease(path)
        current_pid = int(current["pid"])
        if pid_alive(current_pid):
            return (
                False,
                "GPU lane busy "
                f"(owner={current['owner']} pid={current_pid} "
                f"expected_gb={current['expected_gb']})",
            )
        _record_stale_break(path, current, now=now or _utc_now())
        path.unlink()
        return (
            True,
            f"GPU lane free (broke stale owner={current['owner']} pid={current_pid})",
        )


def apply_cupy_memory_limit(gpu_mem_gb: float) -> int:
    """Apply the process CuPy device-pool ceiling and return its byte value."""
    _validate_gpu_mem_gb(gpu_mem_gb)
    try:
        import cupy
    except ImportError as exc:
        raise RuntimeError(
            "--gpu requires the quarry GPU extra with CuPy installed"
        ) from exc
    limit_bytes = int(gpu_mem_gb * 1024**3)
    cupy.get_default_memory_pool().set_limit(size=limit_bytes)
    return limit_bytes


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
    gpu_lease: GpuLease | None = None
    _closed: bool = False

    def close(self) -> None:
        if self._closed:
            return
        if self.gpu_lease is not None:
            self.gpu_lease.release()
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
    gpu_owner: str | None = None,
    gpu_ttl_hours: float = DEFAULT_GPU_TTL_HOURS,
) -> EtiquetteSession:
    """Pre-parse resource flags, enforce them, and tee campaign output.

    Entry-point scripts call this before importing NumPy, PySCF, or quarry
    modules that may import them. Unknown arguments are deliberately retained
    for each script's full parser. Drivers that pass ``gpu_owner`` acquire the
    single GPU lane and apply the CuPy pool limit whenever ``--gpu`` is present.
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--nice", type=int, default=10)
    parser.add_argument("--log")
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--gpu-mem-gb", type=float, default=DEFAULT_GPU_MEM_GB)
    args, _ = parser.parse_known_args(argv)

    if args.threads < 1:
        parser.error("--threads must be at least 1")
    geteuid = getattr(os, "geteuid", None)
    if args.nice < 0 and (geteuid is None or geteuid() != 0):
        parser.error("negative --nice requires root")
    if args.gpu and gpu_owner is not None:
        try:
            _validate_gpu_mem_gb(args.gpu_mem_gb)
        except ValueError as exc:
            parser.error(str(exc))

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
    try:
        enforce(threads=args.threads, niceness=args.nice)
        print(f"campaign log: {log_path}", flush=True)
        limits = " ".join(f"{name}={os.environ[name]}" for name in THREAD_ENV_VARS)
        print(f"compute etiquette: {limits}; nice increment={args.nice}", flush=True)
        if args.gpu and gpu_owner is not None:
            session.gpu_lease = acquire_gpu(
                gpu_owner,
                args.gpu_mem_gb,
                gpu_ttl_hours,
            )
            limit_bytes = apply_cupy_memory_limit(args.gpu_mem_gb)
            print(
                "GPU lease acquired: "
                f"owner={gpu_owner} pid={os.getpid()} "
                f"pool_limit={limit_bytes} bytes "
                f"ollama_floor={GPU_OLLAMA_FLOOR_GB:g} GB",
                flush=True,
            )
    except GpuLeaseBusy as exc:
        session.close()
        print(str(exc), file=sys.stderr)
        raise SystemExit(1) from None
    except BaseException:
        session.close()
        raise
    return session


def _preflight_main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Check the shared quarry GPU lane")
    parser.parse_args(argv)
    try:
        available, message = gpu_lane_available()
    except GpuLeaseBusy as exc:
        print(str(exc), file=sys.stderr)
        return 1
    print(message)
    return 0 if available else 1


if __name__ == "__main__":
    raise SystemExit(_preflight_main())
