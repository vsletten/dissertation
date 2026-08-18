"""SQLite store of structures, jobs, and results (SURVEY.md §4.3).

One user, one box: plain sqlite3, no ORM, no services. Structures are
kept as light XYZ text with a content hash; every derived number carries
its provenance (method string, geometry hash, timestamp). Big artifacts
(checkpoints, cubes) stay on disk outside git and outside this store.
"""

from __future__ import annotations

import hashlib
import sqlite3
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

_SCHEMA = """
CREATE TABLE IF NOT EXISTS structures (
    id INTEGER PRIMARY KEY,
    name TEXT NOT NULL,
    formula TEXT NOT NULL,
    charge INTEGER NOT NULL DEFAULT 0,
    spin INTEGER NOT NULL DEFAULT 0,          -- 2S (pyscf convention)
    xyz TEXT NOT NULL,                        -- plain XYZ block, Angstrom
    geometry_hash TEXT NOT NULL UNIQUE,       -- sha256 of the canonical XYZ
    created_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS jobs (
    id INTEGER PRIMARY KEY,
    structure_id INTEGER NOT NULL REFERENCES structures(id),
    kind TEXT NOT NULL,                       -- opt | freq | ts | irc | sp | screen
    method TEXT NOT NULL,                     -- e.g. "wb97m-v/def2-tzvpd/smd"
    engine TEXT NOT NULL,                     -- pyscf | gpu4pyscf | orca | xtb
    status TEXT NOT NULL DEFAULT 'pending',   -- pending | running | done | failed
    detail TEXT,                              -- free-form JSON blob for settings
    created_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS results (
    id INTEGER PRIMARY KEY,
    job_id INTEGER NOT NULL REFERENCES jobs(id),
    key TEXT NOT NULL,                        -- e.g. "energy", "gibbs", "imag_freq"
    value REAL NOT NULL,
    units TEXT NOT NULL,
    created_at TEXT NOT NULL,
    UNIQUE (job_id, key)
);
"""


def geometry_hash(xyz: str) -> str:
    """sha256 of a whitespace-canonicalized XYZ block."""
    canonical = "\n".join(" ".join(line.split()) for line in xyz.strip().splitlines())
    return hashlib.sha256(canonical.encode()).hexdigest()


def _now() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


@dataclass(frozen=True)
class Provenance:
    """What a number in an emitted deck traces back to."""

    method: str
    engine: str
    geometry_hash: str
    date: str


class Store:
    def __init__(self, path: str | Path) -> None:
        self.conn = sqlite3.connect(str(path))
        self.conn.row_factory = sqlite3.Row
        self.conn.executescript(_SCHEMA)

    def close(self) -> None:
        self.conn.close()

    def __enter__(self) -> Store:
        return self

    def __exit__(self, *exc: object) -> None:
        self.conn.commit()
        self.close()

    # -- structures --------------------------------------------------------
    def add_structure(
        self,
        name: str,
        formula: str,
        xyz: str,
        *,
        charge: int = 0,
        spin: int = 0,
    ) -> int:
        """Insert a structure; returns its id. Idempotent on geometry hash."""
        ghash = geometry_hash(xyz)
        row = self.conn.execute(
            "SELECT id FROM structures WHERE geometry_hash = ?", (ghash,)
        ).fetchone()
        if row:
            return int(row["id"])
        cur = self.conn.execute(
            "INSERT INTO structures (name, formula, charge, spin, xyz,"
            " geometry_hash, created_at) VALUES (?, ?, ?, ?, ?, ?, ?)",
            (name, formula, charge, spin, xyz, ghash, _now()),
        )
        self.conn.commit()
        return int(cur.lastrowid)

    def get_structure(self, structure_id: int) -> sqlite3.Row:
        row = self.conn.execute(
            "SELECT * FROM structures WHERE id = ?", (structure_id,)
        ).fetchone()
        if row is None:
            raise KeyError(f"no structure with id {structure_id}")
        return row

    # -- jobs --------------------------------------------------------------
    def add_job(
        self,
        structure_id: int,
        kind: str,
        method: str,
        engine: str,
        *,
        detail: str | None = None,
    ) -> int:
        cur = self.conn.execute(
            "INSERT INTO jobs (structure_id, kind, method, engine, detail,"
            " created_at) VALUES (?, ?, ?, ?, ?, ?)",
            (structure_id, kind, method, engine, detail, _now()),
        )
        self.conn.commit()
        return int(cur.lastrowid)

    def set_job_status(self, job_id: int, status: str) -> None:
        if status not in ("pending", "running", "done", "failed"):
            raise ValueError(f"unknown job status '{status}'")
        self.conn.execute("UPDATE jobs SET status = ? WHERE id = ?", (status, job_id))
        self.conn.commit()

    # -- results -----------------------------------------------------------
    def add_result(self, job_id: int, key: str, value: float, units: str) -> None:
        self.conn.execute(
            "INSERT INTO results (job_id, key, value, units, created_at)"
            " VALUES (?, ?, ?, ?, ?)"
            " ON CONFLICT (job_id, key) DO UPDATE SET"
            " value = excluded.value, units = excluded.units,"
            " created_at = excluded.created_at",
            (job_id, key, value, units, _now()),
        )
        self.conn.commit()

    def get_result(self, job_id: int, key: str) -> tuple[float, str]:
        row = self.conn.execute(
            "SELECT value, units FROM results WHERE job_id = ? AND key = ?",
            (job_id, key),
        ).fetchone()
        if row is None:
            raise KeyError(f"no result '{key}' for job {job_id}")
        return float(row["value"]), str(row["units"])

    def provenance(self, job_id: int) -> Provenance:
        """The provenance stamp for numbers produced by a job."""
        row = self.conn.execute(
            "SELECT j.method, j.engine, j.created_at, s.geometry_hash"
            " FROM jobs j JOIN structures s ON s.id = j.structure_id"
            " WHERE j.id = ?",
            (job_id,),
        ).fetchone()
        if row is None:
            raise KeyError(f"no job with id {job_id}")
        return Provenance(
            method=row["method"],
            engine=row["engine"],
            geometry_hash=row["geometry_hash"],
            date=row["created_at"],
        )
