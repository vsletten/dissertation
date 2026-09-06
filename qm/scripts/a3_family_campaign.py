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
import sqlite3
import subprocess
import sys
import time
from collections.abc import Callable
from datetime import UTC, datetime, timedelta
from pathlib import Path

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    _ETIQUETTE = bootstrap_cli(
        "a3-family-campaign",
        default_run_root=Path.cwd() / "qm/runs",
    )

from quarry import etiquette
from quarry.pipeline import HARTREE_TO_KJ
from quarry.store import geometry_hash

FAMILY_SITE_KINDS = {"oss": "Oss", "osa": "Osa", "oaa": "Oaa"}
FAMILY_CELL_CENTERS = {
    "oss": {1: None, 2: None, 3: None, 4: None},
    "osa": {1: None, 2: None, 3: None, 4: None},
    # The first Oaa center realizes exact n=2/6, while center 23 is the first
    # crystallographic center that realizes n=4 without silently aliasing n=2.
    "oaa": {2: 18, 4: 23, 6: 18},
}
FAMILY_CONNECTIVITIES = {
    family: tuple(cells) for family, cells in FAMILY_CELL_CENTERS.items()
}
OAA_EXPECTED_BUILDS = {
    2: {
        "center_site": 18,
        "site_kind": "Oaa",
        "formula": "Al4H26O23Si2",
        "n_atoms": 55,
        "n_frozen": 17,
        "attacked_metal": "Al",
    },
    4: {
        "center_site": 23,
        "site_kind": "Oaa",
        "formula": "Al5H31O29Si3",
        "n_atoms": 68,
        "n_frozen": 25,
        "attacked_metal": "Al",
    },
    6: {
        "center_site": 18,
        "site_kind": "Oaa",
        "formula": "Al6H36O35Si4",
        "n_atoms": 81,
        "n_frozen": 33,
        "attacked_metal": "Al",
    },
}
SCHEMA = "a3-family-campaign-terminal-v1"
PROGRESS_SCHEMA = "a3-family-campaign-progress-v1"
TS_GUESS_ROUTES = frozenset({"direct", "proton-neb"})

_STORE_COLUMNS = {
    "structures": (
        ("id", "INTEGER", False, None, 1),
        ("name", "TEXT", True, None, 0),
        ("formula", "TEXT", True, None, 0),
        ("charge", "INTEGER", True, "0", 0),
        ("spin", "INTEGER", True, "0", 0),
        ("xyz", "TEXT", True, None, 0),
        ("geometry_hash", "TEXT", True, None, 0),
        ("created_at", "TEXT", True, None, 0),
    ),
    "jobs": (
        ("id", "INTEGER", False, None, 1),
        ("structure_id", "INTEGER", True, None, 0),
        ("kind", "TEXT", True, None, 0),
        ("method", "TEXT", True, None, 0),
        ("engine", "TEXT", True, None, 0),
        ("status", "TEXT", True, "'pending'", 0),
        ("detail", "TEXT", False, None, 0),
        ("created_at", "TEXT", True, None, 0),
    ),
    "results": (
        ("id", "INTEGER", False, None, 1),
        ("job_id", "INTEGER", True, None, 0),
        ("key", "TEXT", True, None, 0),
        ("value", "REAL", True, None, 0),
        ("units", "TEXT", True, None, 0),
        ("created_at", "TEXT", True, None, 0),
    ),
}
_STORE_UNIQUE_KEYS = {
    "structures": {("geometry_hash",)},
    "jobs": set(),
    "results": {("job_id", "key")},
}
_STORE_FOREIGN_KEYS = {
    "structures": set(),
    "jobs": {("structures", "structure_id", "id", "NO ACTION", "NO ACTION")},
    "results": {("jobs", "job_id", "id", "NO ACTION", "NO ACTION")},
}


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


def observed_git_sha(worktree: Path) -> str | None:
    try:
        return git_output(worktree, "rev-parse", "HEAD")
    except (OSError, subprocess.SubprocessError):
        return None


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


def _finite_scientific_value(value: object) -> bool:
    return (
        not isinstance(value, bool)
        and isinstance(value, (int, float))
        and math.isfinite(float(value))
    )


def _validate_store_schema(connection: sqlite3.Connection) -> None:
    tables = {
        row[0]
        for row in connection.execute(
            "SELECT name FROM sqlite_master "
            "WHERE type = 'table' AND name NOT LIKE 'sqlite_%'"
        )
    }
    required_tables = set(_STORE_COLUMNS)
    if not required_tables <= tables:
        missing = sorted(required_tables - tables)
        raise RuntimeError(f"store schema is incomplete: missing {missing}")
    if tables != required_tables:
        raise RuntimeError(f"store schema has unexpected tables: {sorted(tables)}")

    for table, expected_columns in _STORE_COLUMNS.items():
        observed_columns = tuple(
            (
                str(row["name"]),
                str(row["type"]).upper(),
                bool(row["notnull"]),
                row["dflt_value"],
                int(row["pk"]),
            )
            # table is a _STORE_COLUMNS key; bind it into pragma_table_info
            for row in connection.execute(
                "SELECT * FROM pragma_table_info(?)",
                (table,),
            )
        )
        if observed_columns != expected_columns:
            raise RuntimeError(f"store schema mismatch for {table}")

        unique_keys = set()
        # table is a _STORE_COLUMNS key; bind it into pragma_index_list
        for index in connection.execute(
            "SELECT * FROM pragma_index_list(?)",
            (table,),
        ):
            if not index["unique"]:
                continue
            columns = tuple(
                str(row["name"])
                # index name comes from the untrusted store; bind it into pragma_index_info
                for row in connection.execute(
                    "SELECT * FROM pragma_index_info(?)",
                    (index["name"],),
                )
            )
            unique_keys.add(columns)
        if unique_keys != _STORE_UNIQUE_KEYS[table]:
            raise RuntimeError(f"store unique-key mismatch for {table}")

        foreign_keys = {
            (
                str(row["table"]),
                str(row["from"]),
                str(row["to"]),
                str(row["on_update"]),
                str(row["on_delete"]),
            )
            # table is a _STORE_COLUMNS key; bind it into pragma_foreign_key_list
            for row in connection.execute(
                "SELECT * FROM pragma_foreign_key_list(?)",
                (table,),
            )
        }
        if foreign_keys != _STORE_FOREIGN_KEYS[table]:
            raise RuntimeError(f"store foreign-key schema mismatch for {table}")

    violations = connection.execute("PRAGMA foreign_key_check").fetchall()
    if violations:
        raise RuntimeError(
            f"store foreign-key check failed: {len(violations)} violation(s)"
        )


def _parse_store_timestamp(value: object, *, label: str) -> datetime:
    if not isinstance(value, str):
        raise RuntimeError(f"store {label} timestamp is not text")
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError as exc:
        raise RuntimeError(f"store {label} timestamp is invalid: {value!r}") from exc
    if parsed.tzinfo is None or parsed.utcoffset() != timedelta(0):
        raise RuntimeError(f"store {label} timestamp is not UTC-aware: {value!r}")
    return parsed


def _xyz_symbols(xyz: object, *, label: str) -> tuple[str, ...]:
    if not isinstance(xyz, str):
        raise RuntimeError(f"store {label} XYZ is not text")
    lines = xyz.strip().splitlines()
    try:
        atom_count = int(lines[0])
    except (IndexError, ValueError) as exc:
        raise RuntimeError(f"store {label} XYZ atom count is invalid") from exc
    if atom_count < 1 or len(lines) != atom_count + 2:
        raise RuntimeError(f"store {label} XYZ line count is invalid")
    symbols: list[str] = []
    for line in lines[2:]:
        fields = line.split()
        if len(fields) != 4 or not fields[0]:
            raise RuntimeError(f"store {label} XYZ atom record is invalid")
        try:
            coordinates = tuple(float(value) for value in fields[1:])
        except ValueError as exc:
            raise RuntimeError(f"store {label} XYZ coordinate is invalid") from exc
        if not all(math.isfinite(value) for value in coordinates):
            raise RuntimeError(f"store {label} XYZ coordinate is non-finite")
        symbols.append(fields[0])
    return tuple(symbols)


def _formula(symbols: tuple[str, ...]) -> str:
    counts: dict[str, int] = {}
    for symbol in symbols:
        counts[symbol] = counts.get(symbol, 0) + 1
    return "".join(
        f"{symbol}{count if count > 1 else ''}"
        for symbol, count in sorted(counts.items())
    )


def validate_result(
    result_path: Path,
    *,
    family: str,
    state: str,
    n_intact: int,
) -> dict[str, object]:
    store_path = result_path.with_name("store.sqlite")
    geometry_paths = {
        role: result_path.with_name(f"{role}.xyz") for role in ("complex", "ts")
    }
    if (
        not result_path.is_file()
        or not store_path.is_file()
        or any(not path.is_file() for path in geometry_paths.values())
    ):
        raise RuntimeError(
            f"n={n_intact} exited zero without results.json + store.sqlite + "
            "complex.xyz + ts.xyz"
        )
    payload = json.loads(result_path.read_text())
    if not isinstance(payload, dict):
        raise RuntimeError(f"n={n_intact} result payload is not an object")
    expected = {"family": family, "state": state, "n_intact": n_intact}
    observed = {key: payload.get(key) for key in expected}
    if type(payload.get("n_intact")) is not int or observed != expected:
        raise RuntimeError(
            f"n={n_intact} result identity mismatch: "
            f"expected={expected} observed={observed}"
        )

    required_finite = (
        "dE_elec_vs_complex_kj",
        "dG_kj",
        "dH_kj",
        "ts_imaginary_cm",
    )
    for key in required_finite:
        value = payload.get(key)
        if not _finite_scientific_value(value):
            raise RuntimeError(f"n={n_intact} result has invalid {key}={value!r}")
    route = payload.get("route")
    if route not in TS_GUESS_ROUTES:
        raise RuntimeError(
            f"n={n_intact} result has invalid route provenance {route!r}"
        )

    cell = payload.get("cell")
    method = payload.get("method")
    expected_cell = f"{family}-{state}-n{n_intact}-s2"
    expected_method = "b3lyp/def2-svp/df"
    if cell != expected_cell:
        raise RuntimeError(
            f"n={n_intact} result has wrong cell identity (expected {expected_cell})"
        )
    if method != expected_method:
        raise RuntimeError(
            f"n={n_intact} result has wrong method (expected {expected_method})"
        )

    geometry_hashes = payload.get("geometry_hash")
    if not isinstance(geometry_hashes, dict) or set(geometry_hashes) != {
        "complex",
        "ts",
    }:
        raise RuntimeError(f"n={n_intact} result has invalid geometry provenance")
    if any(
        not isinstance(value, str)
        or len(value) != 64
        or value != value.lower()
        or any(character not in "0123456789abcdef" for character in value)
        for value in geometry_hashes.values()
    ):
        raise RuntimeError(f"n={n_intact} result has invalid geometry hashes")
    durable_geometries: dict[str, str] = {}
    for role, path in geometry_paths.items():
        try:
            xyz = path.read_text()
            _xyz_symbols(xyz, label=f"durable {role}")
        except (OSError, RuntimeError) as exc:
            raise RuntimeError(
                f"n={n_intact} durable {role} geometry is invalid: {exc}"
            ) from exc
        if geometry_hash(xyz) != geometry_hashes[role]:
            raise RuntimeError(
                f"n={n_intact} durable {role} geometry hash mismatches results.json"
            )
        durable_geometries[role] = xyz

    cluster = payload.get("cluster")
    if not isinstance(cluster, dict):
        raise RuntimeError(f"n={n_intact} result has no cluster provenance")
    expected_cluster = {
        "site_kind": FAMILY_SITE_KINDS[family],
        "metal_shells": 2,
        "n_intact_requested": n_intact,
        "n_intact": n_intact,
        "charge": 0,
        "state": state,
        "charge_offset": 1 if state == "acid" else 0,
        "method": method,
        "gpu": True,
    }
    observed_cluster = {key: cluster.get(key) for key in expected_cluster}
    if observed_cluster != expected_cluster or any(
        type(observed_cluster[key]) is not type(expected_value)
        for key, expected_value in expected_cluster.items()
    ):
        raise RuntimeError(
            f"n={n_intact} result cluster identity mismatch: "
            f"expected={expected_cluster} observed={observed_cluster}"
        )
    center_site = cluster.get("center_site")
    n_atoms = cluster.get("n_atoms")
    n_frozen = cluster.get("n_frozen")
    if type(center_site) is not int or center_site < 0:
        raise RuntimeError(f"n={n_intact} result has invalid crystallographic center")
    if type(n_atoms) is not int or n_atoms < 1:
        raise RuntimeError(f"n={n_intact} result has invalid cluster atom count")
    if type(n_frozen) is not int or not 0 < n_frozen <= n_atoms:
        raise RuntimeError(f"n={n_intact} result has invalid frozen atom count")
    expected_center = FAMILY_CELL_CENTERS[family][n_intact]
    if expected_center is not None and center_site != expected_center:
        raise RuntimeError(
            f"n={n_intact} result has wrong crystallographic center "
            f"(expected {expected_center})"
        )
    expected_metal = "Al" if family == "oaa" else "Si"
    if (
        type(payload.get("metal_shells")) is not int
        or payload.get("metal_shells") != 2
        or payload.get("attacked_metal") != expected_metal
    ):
        raise RuntimeError(f"n={n_intact} result has wrong attacked-site identity")
    if family == "oaa":
        expected_build = OAA_EXPECTED_BUILDS[n_intact]
        observed_build = {
            "center_site": center_site,
            "site_kind": cluster.get("site_kind"),
            "formula": cluster.get("formula"),
            "n_atoms": n_atoms,
            "n_frozen": n_frozen,
            "attacked_metal": payload.get("attacked_metal"),
        }
        if observed_build != expected_build or any(
            type(observed_build[key]) is not type(expected_value)
            for key, expected_value in expected_build.items()
        ):
            raise RuntimeError(
                f"n={n_intact} Oaa physical build identity mismatch: "
                f"expected={expected_build} observed={observed_build}"
            )

    connection: sqlite3.Connection | None = None
    try:
        connection = sqlite3.connect(
            f"{store_path.resolve().as_uri()}?mode=ro", uri=True
        )
        connection.row_factory = sqlite3.Row
        with connection:
            integrity = connection.execute("PRAGMA integrity_check").fetchone()[0]
            if integrity != "ok":
                raise RuntimeError(f"store integrity check failed: {integrity}")
            _validate_store_schema(connection)
            # Identifiers cannot be bound; table names are frozen _STORE_COLUMNS keys.
            counts = {
                "structures": int(
                    connection.execute("SELECT count(*) FROM structures").fetchone()[0]
                ),
                "jobs": int(
                    connection.execute("SELECT count(*) FROM jobs").fetchone()[0]
                ),
                "results": int(
                    connection.execute("SELECT count(*) FROM results").fetchone()[0]
                ),
            }
            if counts != {"structures": 2, "jobs": 2, "results": 2}:
                raise RuntimeError(f"store row cardinality mismatch: {counts}")
            rows = connection.execute(
                "SELECT "
                "s.id AS structure_id, s.name, s.formula, s.charge, s.spin, "
                "s.xyz, s.geometry_hash, s.created_at AS structure_created_at, "
                "typeof(s.charge) AS charge_type, typeof(s.spin) AS spin_type, "
                "j.id AS job_id, j.structure_id AS job_structure_id, j.kind, "
                "j.method, j.engine, j.status, j.detail, "
                "j.created_at AS job_created_at, "
                "r.id AS result_id, r.job_id AS result_job_id, r.key, r.value, "
                "r.units, r.created_at AS result_created_at, "
                "typeof(r.value) AS value_type "
                "FROM structures s JOIN jobs j ON j.structure_id = s.id "
                "JOIN results r ON r.job_id = j.id ORDER BY s.name"
            ).fetchall()
    except (sqlite3.DatabaseError, OSError, RuntimeError) as exc:
        raise RuntimeError(f"n={n_intact} provenance store is invalid: {exc}") from exc
    finally:
        if connection is not None:
            connection.close()

    expected_names = {f"{cell}-complex", f"{cell}-ts"}
    if len(rows) != 2 or {row["name"] for row in rows} != expected_names:
        raise RuntimeError(f"n={n_intact} store lacks exact complex/TS provenance rows")
    expected_hashes = {
        f"{cell}-complex": geometry_hashes["complex"],
        f"{cell}-ts": geometry_hashes["ts"],
    }
    expected_charge = 1 if state == "acid" else 0
    expected_atom_count = n_atoms + (4 if state == "acid" else 3)
    symbols_by_name: dict[str, tuple[str, ...]] = {}
    energy_by_name: dict[str, float] = {}
    for row in rows:
        name = str(row["name"])
        role = name.removeprefix(f"{cell}-")
        symbols = _xyz_symbols(row["xyz"], label=name)
        symbols_by_name[name] = symbols
        structure_time = _parse_store_timestamp(
            row["structure_created_at"], label=f"{name} structure"
        )
        job_time = _parse_store_timestamp(row["job_created_at"], label=f"{name} job")
        result_time = _parse_store_timestamp(
            row["result_created_at"], label=f"{name} result"
        )
        if not structure_time <= job_time <= result_time:
            raise RuntimeError(
                f"n={n_intact} store timestamp order is invalid for {name}"
            )
        if (
            row["job_structure_id"] != row["structure_id"]
            or row["result_job_id"] != row["job_id"]
            or row["charge_type"] != "integer"
            or row["spin_type"] != "integer"
            or row["value_type"] != "real"
            or row["charge"] != expected_charge
            or row["spin"] != 0
            or row["formula"] != _formula(symbols)
            or len(symbols) != expected_atom_count
            or row["xyz"] != durable_geometries[role]
            or row["geometry_hash"] != expected_hashes[name]
            or row["geometry_hash"] != geometry_hash(str(row["xyz"]))
            or row["kind"] != "freq"
            or row["method"] != method
            or row["engine"] != "gpu4pyscf"
            or row["status"] != "done"
            or row["detail"] is not None
            or row["key"] != "electronic"
            or row["units"] != "hartree"
            or not _finite_scientific_value(row["value"])
        ):
            raise RuntimeError(f"n={n_intact} store provenance mismatch for {name}")
        energy_by_name[name] = float(row["value"])

    if symbols_by_name[f"{cell}-complex"] != symbols_by_name[f"{cell}-ts"]:
        raise RuntimeError(f"n={n_intact} store complex/TS atom identity mismatch")
    substrate_counts: dict[str, int] = {}
    for symbol in symbols_by_name[f"{cell}-complex"]:
        substrate_counts[symbol] = substrate_counts.get(symbol, 0) + 1
    substrate_counts["O"] = substrate_counts.get("O", 0) - 1
    substrate_counts["H"] = substrate_counts.get("H", 0) - (3 if state == "acid" else 2)
    if any(count < 0 for count in substrate_counts.values()):
        raise RuntimeError(f"n={n_intact} store is missing the declared attacker atoms")
    substrate_symbols = tuple(
        symbol for symbol, count in substrate_counts.items() for _ in range(count)
    )
    if cluster.get("formula") != _formula(substrate_symbols):
        raise RuntimeError(f"n={n_intact} store composition mismatches cluster formula")
    store_barrier_kj = (
        energy_by_name[f"{cell}-ts"] - energy_by_name[f"{cell}-complex"]
    ) * HARTREE_TO_KJ
    if not math.isclose(
        store_barrier_kj,
        float(payload["dE_elec_vs_complex_kj"]),
        rel_tol=1e-12,
        abs_tol=1e-8,
    ):
        raise RuntimeError(
            f"n={n_intact} store electronic barrier does not match results.json"
        )
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
    ap.add_argument("--family", choices=sorted(FAMILY_CONNECTIVITIES), required=True)
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
    allowed = FAMILY_CONNECTIVITIES[args.family]
    if not cells or any(cell not in allowed for cell in cells):
        allowed_text = ",".join(str(cell) for cell in allowed)
        raise SystemExit(f"--cells must be drawn from {allowed_text} for {args.family}")
    if args.family == "oaa" and cells != allowed[: len(cells)]:
        allowed_text = ",".join(str(cell) for cell in allowed)
        raise SystemExit(
            f"Oaa --cells must be a serial prefix of {allowed_text} starting at n=2"
        )
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
            center_index = FAMILY_CELL_CENTERS.get(args.family, {}).get(n_intact)
            if center_index is not None:
                command.extend(["--center-index", str(center_index)])
            cell_started = time.monotonic()
            # argv sequence, shell=False: trusted interpreter plus in-repo
            # phase2_ladder.py driver and static campaign flags. Not a shell string.
            # nosemgrep: python.lang.security.audit.dangerous-subprocess-use-audit
            result = subprocess.run(
                command, cwd=args.worktree, check=False, shell=False
            )
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
        deferred_signals: list[int] = []

        def defer_stop(signum: int, _frame: object) -> None:
            # The work is already terminal. Do not permit a signal to split
            # final progress from the receipt commit point.
            deferred_signals.append(signum)

        for sig in previous_handlers:
            signal.signal(sig, defer_stop)
        receipt: dict[str, object] = {
            "schema": SCHEMA,
            "started_at": started,
            "completed_at": now(),
            "elapsed_seconds": time.monotonic() - started_monotonic,
            "expected_git_sha": args.expected_git_sha,
            "observed_git_sha": observed_git_sha(args.worktree),
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
        try:
            heartbeat("complete" if receipt["success"] else "failed")
        except Exception as exc:
            progress_error = f"{type(exc).__name__}: {exc}"
            receipt["success"] = False
            receipt["terminal_reason"] = "progress-finalization-failed"
            receipt["exit_code"] = exit_code = exit_code if exit_code != 0 else 1
            receipt["error"] = (
                f"{error}; final progress write failed: {progress_error}"
                if error
                else f"final progress write failed: {progress_error}"
            )
        try:
            # This is the terminal commit point: after it exists, progress is
            # already final or the receipt itself records why that failed.
            atomic_json(terminal_path, receipt)
        finally:
            for sig, previous in previous_handlers.items():
                signal.signal(sig, previous)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
