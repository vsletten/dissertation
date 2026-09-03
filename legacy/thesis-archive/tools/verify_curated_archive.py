#!/usr/bin/env python3
"""Independent acceptance checks for the curated thesis archive.

This intentionally does not import either generator. It derives expectations from
the source tree and validates committed artifacts through a separate code path.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import re
import subprocess

EXPECTED_MANIFEST_SHA256 = "edf7dea5a690823c9a131d0e28cac5ca68634a0538b3df424dd8c91235198fdb"
LEDGER_FIELDS = [
    "system_id", "reaction_id", "ligand_class", "ref_log", "ts_log", "product_log",
    "method", "basis", "energy_treatment", "e_ref_hartree", "e_ts_hartree",
    "e_product_hartree", "barrier_forward_kj_mol", "barrier_reverse_kj_mol",
    "reaction_energy_kj_mol", "imaginary_frequencies_cm_1", "normal_termination",
    "parse_status", "notes",
]
HARTREE_TO_KJ_MOL = 2625.4996394799
GOLDEN = {
    "hotrox/935077498", "hotrox/936930575", "hotrox/937172019",
    "jasper/933892971", "jasper/935835145",
}


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def metadata_manifest(source: Path) -> bytes:
    result = subprocess.run(
        ["find", str(source), "-printf", "%P\t%y\t%s\t%T@\\n"],
        check=True,
        capture_output=True,
    )
    return b"".join(sorted(result.stdout.splitlines(keepends=True)))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def verify_inventory(source: Path, curated: Path) -> list[str]:
    rows = read_csv(curated / "ARCHIVE-INVENTORY.csv")
    source_files = {path.relative_to(source).as_posix(): path for path in source.rglob("*") if path.is_file()}
    assert len(rows) == len(source_files) == 1426
    assert {row["path"] for row in rows} == set(source_files)
    for row in rows:
        path = source_files[row["path"]]
        stat = path.stat()
        assert int(row["bytes"]) == stat.st_size
        assert row["sha256"] == sha256(path)
        expected_time = datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(timespec="microseconds")
        assert row["mtime_utc"] == expected_time
    return [f"inventory: {len(rows)} source files independently rehashed"]


def verify_dft(source: Path, curated: Path) -> list[str]:
    coverage = read_csv(curated / "DFT-SOURCE-LOGS.csv")
    source_logs = {
        path.relative_to(source).as_posix()
        for path in (source / "dft" / "results").rglob("*.log")
    }
    assert len(coverage) == len(source_logs) == 83
    assert len({row["log"] for row in coverage}) == len(coverage)
    assert {row["log"] for row in coverage} == source_logs

    ledger_path = curated / "DFT-LEDGER.csv"
    with ledger_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames == LEDGER_FIELDS
        ledger = list(reader)
    assert len(ledger) == 89
    reactions = [row for row in ledger if row["reaction_id"].startswith("ch4-")]
    sources = [row for row in ledger if row["reaction_id"] == "source-calculation"]
    assert len(reactions) == 45 and len(sources) == 44
    for row in ledger:
        assert row["parse_status"] in {"ok", "unmatched", "truncated", "ambiguous", "parse_error"}
        assert row["normal_termination"] in {"true", "false"}
    for row in reactions:
        assert not row["ts_log"]
        assert not row["e_ts_hartree"]
        assert not row["barrier_forward_kj_mol"]
        assert not row["barrier_reverse_kj_mol"]
        assert row["parse_status"] == "unmatched"
        expected = (float(row["e_product_hartree"]) - float(row["e_ref_hartree"])) * HARTREE_TO_KJ_MOL
        assert abs(expected - float(row["reaction_energy_kj_mol"])) <= 0.002

    # Independent premise check: no Gaussian input route requests a TS/QST job.
    ts_routes = []
    for path in (source / "dft" / "results").rglob("*.com"):
        text = path.read_text(encoding="latin-1").lower()
        if re.search(r"\b(?:qst2|qst3)\b|opt\s*=\s*\([^)]*\bts\b", text):
            ts_routes.append(path)
    assert not ts_routes, ts_routes
    thesis = (source / "text" / "ch4.tex").read_text(encoding="latin-1")
    assert "focuses on total reaction energies rather than transition" in thesis
    return [
        "dft: 83 source logs covered exactly once",
        "dft: fixed 19-column schema; 44 source rows + 45 reaction rows",
        "dft: reaction arithmetic independently recomputed; no TS/QST input route found",
    ]


def verify_copies(source: Path, curated: Path) -> list[str]:
    pairs: list[tuple[Path, Path]] = []
    for path in (curated / "c-model" / "source").iterdir():
        pairs.append((source / "mc" / "model" / path.name, path))
    for path in (curated / "c-model" / "appendix-source").iterdir():
        pairs.append((source / "text" / "src" / path.name, path))
    for path in (curated / "c-model" / "RCS").iterdir():
        pairs.append((source / "mc" / "model" / "RCS" / path.name, path))
    for path in (curated / "c-model" / "input").iterdir():
        pairs.append((source / "mc" / "model" / "input" / path.name, path))
    assert len(pairs) == 96
    for original, copy in pairs:
        assert original.read_bytes() == copy.read_bytes(), copy

    golden_root = curated / "golden-runs"
    actual_runs = {
        path.relative_to(golden_root).as_posix()
        for path in golden_root.glob("*/*")
        if path.is_dir()
    }
    assert actual_runs == GOLDEN
    copied_run_files = 0
    for run in GOLDEN:
        original = source / "mc" / "results" / run
        copy = golden_root / run
        assert not (copy / "plot.ps").exists()
        for path in copy.iterdir():
            assert path.read_bytes() == (original / path.name).read_bytes(), path
            copied_run_files += 1
        assert all((copy / name).exists() for name in ("data.cell", "data.lattice", "data.rxn", "data.sim", "results.dat", "surfAl.out", "surfSi.out"))
    return [f"copies: {len(pairs)} model/input files and {copied_run_files} golden-run files match source bytes"]


def pdf_pages(path: Path) -> int:
    result = subprocess.run(["pdfinfo", str(path)], check=True, capture_output=True, text=True)
    match = re.search(r"^Pages:\s+(\d+)$", result.stdout, re.MULTILINE)
    assert match
    return int(match.group(1))


def verify_pdfs(source: Path, curated: Path) -> list[str]:
    archived = curated / "text" / "Playing-Dice-with-the-Universe.archived-2007.pdf"
    rebuilt = curated / "text" / "Playing-Dice-with-the-Universe.rebuilt-2026.pdf"
    assert archived.read_bytes() == (source / "text" / "thesis.pdf").read_bytes()
    assert pdf_pages(archived) == 134
    assert pdf_pages(rebuilt) == 133
    for pdf in (archived, rebuilt):
        result = subprocess.run(["pdftotext", str(pdf), "-"], check=True, capture_output=True)
        assert b"Monte Carlo" in result.stdout and b"Organic Acids" in result.stdout
    return ["pdf: archived 134-page and rebuilt 133-page documents parse and contain thesis text"]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--curated", type=Path, required=True)
    args = parser.parse_args()

    manifest = metadata_manifest(args.source)
    actual_manifest = hashlib.sha256(manifest).hexdigest()
    assert actual_manifest == EXPECTED_MANIFEST_SHA256, actual_manifest
    checks = ["source metadata manifest unchanged"]
    checks += verify_inventory(args.source, args.curated)
    checks += verify_dft(args.source, args.curated)
    checks += verify_copies(args.source, args.curated)
    checks += verify_pdfs(args.source, args.curated)
    for check in checks:
        print(f"PASS: {check}")
    print(f"PASS: {len(checks)} independent acceptance checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
