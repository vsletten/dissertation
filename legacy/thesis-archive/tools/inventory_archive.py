#!/usr/bin/env python3
"""Generate deterministic, read-only inventories of the 1999 thesis archive."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import tempfile

ARCHIVE_FIELDS = ["path", "category", "bytes", "mtime_utc", "sha256", "annotation"]
RUN_FIELDS = [
    "host", "run_id", "run_started_utc", "temperature_K", "delta_mu_si",
    "delta_mu_al", "steps", "seed", "done_sentinel", "required_files_present",
    "results_lines", "surfAl_lines", "surfSi_lines", "bytes",
]
REQUIRED_RUN_FILES = {
    "data.cell", "data.lattice", "data.rxn", "data.sim", "results.dat",
    "surfAl.out", "surfSi.out",
}


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def annotation(path: Path) -> str:
    suffix = path.suffix.lower()
    if path.name.endswith(",v"):
        return "RCS source archive"
    if suffix == ".chk":
        return "Gaussian binary checkpoint; inventoried, not curated"
    if suffix == ".log" and "dft/results" in path.as_posix():
        return "Gaussian output; parsed into DFT-SOURCE-LOGS.csv"
    if suffix == ".com":
        return "Gaussian input deck"
    if suffix in {".wfn", ".grd"}:
        return "quantum-chemistry field/wavefunction output"
    if suffix in {".ps", ".eps", ".jpg", ".pdf", ".dvi"}:
        return "rendered/binary document or figure"
    if suffix in {".c", ".h"} or path.name == "Makefile":
        return "C model source"
    if path.name.startswith("data."):
        return "KMC run input"
    if path.name in {"results.dat", "surfAl.out", "surfSi.out"}:
        return "KMC run output"
    if suffix in {".tex", ".bib", ".sty"}:
        return "thesis/document source"
    return "archive file"


def archive_rows(source: Path) -> list[dict[str, str]]:
    rows = []
    for path in sorted((item for item in source.rglob("*") if item.is_file()), key=lambda p: p.relative_to(source).as_posix()):
        rel = path.relative_to(source)
        stat = path.stat()
        rows.append(
            {
                "path": rel.as_posix(),
                "category": rel.parts[0],
                "bytes": str(stat.st_size),
                "mtime_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(timespec="microseconds"),
                "sha256": digest(path),
                "annotation": annotation(rel),
            }
        )
    return rows


def input_value(path: Path, line_index: int) -> str:
    lines = path.read_text(encoding="latin-1").splitlines()
    if line_index >= len(lines):
        return ""
    return lines[line_index].split("#", 1)[0].strip()


def line_count(path: Path) -> str:
    if not path.exists():
        return ""
    with path.open("rb") as handle:
        return str(sum(1 for _ in handle))


def run_rows(source: Path) -> list[dict[str, str]]:
    rows = []
    results = source / "mc" / "results"
    for host_dir in sorted(item for item in results.iterdir() if item.is_dir()):
        for run_dir in sorted(item for item in host_dir.iterdir() if item.is_dir()):
            files = {item.name for item in run_dir.iterdir() if item.is_file()}
            try:
                started = datetime.fromtimestamp(int(run_dir.name), timezone.utc).isoformat()
            except (ValueError, OSError):
                started = ""
            rxn = run_dir / "data.rxn"
            sim = run_dir / "data.sim"
            rows.append(
                {
                    "host": host_dir.name,
                    "run_id": run_dir.name,
                    "run_started_utc": started,
                    "temperature_K": input_value(rxn, 0) if rxn.exists() else "",
                    "delta_mu_si": input_value(rxn, 1) if rxn.exists() else "",
                    "delta_mu_al": input_value(rxn, 2) if rxn.exists() else "",
                    "steps": input_value(sim, 0) if sim.exists() else "",
                    "seed": input_value(sim, 3) if sim.exists() else "",
                    "done_sentinel": str((run_dir / "done").exists()).lower(),
                    "required_files_present": str(REQUIRED_RUN_FILES <= files).lower(),
                    "results_lines": line_count(run_dir / "results.dat"),
                    "surfAl_lines": line_count(run_dir / "surfAl.out"),
                    "surfSi_lines": line_count(run_dir / "surfSi.out"),
                    "bytes": str(sum(item.stat().st_size for item in run_dir.iterdir() if item.is_file())),
                }
            )
    return rows


def write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    archive = archive_rows(args.source)
    runs = run_rows(args.source)
    assert len(archive) == sum(1 for item in args.source.rglob("*") if item.is_file())
    assert all(set(row) == set(ARCHIVE_FIELDS) for row in archive)
    assert all(set(row) == set(RUN_FIELDS) for row in runs)

    outputs = [
        ("ARCHIVE-INVENTORY.csv", ARCHIVE_FIELDS, archive),
        ("MC-RUN-INVENTORY.csv", RUN_FIELDS, runs),
    ]
    if args.check:
        with tempfile.TemporaryDirectory(prefix="thesis-inventory-check-") as temp:
            for name, fields, rows in outputs:
                candidate = Path(temp) / name
                write_csv(candidate, fields, rows)
                assert (args.out_dir / name).read_bytes() == candidate.read_bytes(), f"stale generated file: {name}"
    else:
        for name, fields, rows in outputs:
            write_csv(args.out_dir / name, fields, rows)
    print(f"validated {len(archive)} archive files and {len(runs)} KMC runs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
