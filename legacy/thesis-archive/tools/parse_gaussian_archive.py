#!/usr/bin/env python3
"""Build the curated 1999 Gaussian DFT ledger.

The archive is read-only. This program reads every ``dft/results/**/*.log`` file,
selects the last normally terminated Link1 job with an energy, and writes a
schema-fixed CSV plus a human-readable ledger. The reaction map is deliberately
small and evidence-backed: it reproduces the reaction-energy equations in
Chapter 4 of the thesis. No transition-state calculation is present, so no
activation barrier is invented.
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, cast

HARTREE_TO_KJ_MOL = 2625.4996394799
FIELDS = [
    "system_id",
    "reaction_id",
    "ligand_class",
    "ref_log",
    "ts_log",
    "product_log",
    "method",
    "basis",
    "energy_treatment",
    "e_ref_hartree",
    "e_ts_hartree",
    "e_product_hartree",
    "barrier_forward_kj_mol",
    "barrier_reverse_kj_mol",
    "reaction_energy_kj_mol",
    "imaginary_frequencies_cm_1",
    "normal_termination",
    "parse_status",
    "notes",
]
VALID_STATUS = {"ok", "unmatched", "truncated", "ambiguous", "parse_error"}
VALID_TREATMENT = {"electronic", "electronic+ZPE"}

SCF_RE = re.compile(r"SCF Done:\s+E\([^)]*\)\s*=\s*([-+]?\d+\.\d+(?:[DEde][-+]?\d+)?)")
ZPE_RE = re.compile(
    r"Sum of electronic and zero-point Energies=\s*([-+]?\d+\.\d+(?:[DEde][-+]?\d+)?)"
)
FREQ_RE = re.compile(r"Frequencies --\s+(.+)")
ROUTE_RE = re.compile(r"^\s*#\s*(.+)$", re.MULTILINE)
METHOD_BASIS_RE = re.compile(r"(?i)\b([ur]?(?:b3lyp|bp86|hf|mp2))/([^\s/]+)")


@dataclass(frozen=True)
class Job:
    index: int
    normal: bool
    electronic: float | None
    electronic_zpe: float | None
    method: str
    basis: str
    imaginary: tuple[float, ...]

    @property
    def energy(self) -> float | None:
        return self.electronic_zpe if self.electronic_zpe is not None else self.electronic

    @property
    def treatment(self) -> str:
        return "electronic+ZPE" if self.electronic_zpe is not None else "electronic"


@dataclass(frozen=True)
class LogResult:
    relpath: str
    jobs: tuple[Job, ...]
    selected: Job | None


@dataclass(frozen=True)
class Species:
    system_id: str
    primary: LogResult | None
    logs: tuple[LogResult, ...]


@dataclass(frozen=True)
class Reaction:
    reaction_id: str
    ligand_class: str
    refs: tuple[str, ...]
    products: tuple[str, ...]
    thesis_kcal_mol: float


# Equations and comparison values are transcribed from text/ch4.tex, Tables 4.4,
# 4.5, 4.7, 4.8, and 4.10. Species IDs are archive directory names.
REACTIONS = (
    Reaction("ch4-water-a", "water", ("h6sialo7", "h4o2"), ("h6sialo7+h2o", "h2o"), -5.9),
    Reaction("ch4-water-b", "water", ("h6sialo7", "h5o2"), ("h7sialo7", "h4o2"), -108.7),
    Reaction("ch4-water-c", "water", ("h7sialo7", "h4o2"), ("h7sialo7+h2o", "h2o"), -11.7),
    Reaction("ch4-oxalate-a", "oxalate", ("h6sialo7", "c2o4"), ("h6c2sialo11",), 75.5),
    Reaction("ch4-oxalate-b", "oxalate", ("h6sialo7+h2o", "c2o4"), ("h6c2sialo11+h2o",), 61.7),
    Reaction("ch4-oxalate-c", "oxalate", ("h6c2sialo11", "h4o2"), ("h6c2sialo11+h2o", "h2o"), -19.6),
    Reaction("ch4-oxalate-d", "oxalate", ("h7sialo7", "c2o4"), ("h7c2sialo11",), -129.0),
    Reaction("ch4-oxalate-e", "oxalate", ("h7sialo7+h2o", "c2o4"), ("h7c2sialo11+h2o",), -127.0),
    Reaction("ch4-oxalate-f", "oxalate", ("h7c2sialo11", "h4o2"), ("h7c2sialo11+h2o", "h2o"), -9.7),
    Reaction("ch4-malonate-a", "malonate", ("h6sialo7", "h2c3o4"), ("h8c3sialo11",), 91.6),
    Reaction("ch4-malonate-b", "malonate", ("h6sialo7+h2o", "h2c3o4"), ("h8c3sialo11+h2o",), 77.8),
    Reaction("ch4-malonate-c", "malonate", ("h8c3sialo11", "h4o2"), ("h8c3sialo11+h2o", "h2o"), -19.7),
    Reaction("ch4-malonate-d", "malonate", ("h7sialo7", "h2c3o4"), ("h9c3sialo11",), -109.3),
    Reaction("ch4-malonate-e", "malonate", ("h7sialo7+h2o", "h2c3o4"), ("h9c3sialo11+h2o",), -102.3),
    Reaction("ch4-malonate-f", "malonate", ("h9c3sialo11", "h4o2"), ("h9c3sialo11+h2o", "h2o"), -4.7),
    Reaction("ch4-succinate-a", "succinate", ("h6sialo7", "h4c4o4"), ("h10c4sialo11",), 108.6),
    Reaction("ch4-succinate-b", "succinate", ("h6sialo7+h2o", "h4c4o4"), ("h10c4sialo11+h2o",), 98.3),
    Reaction("ch4-succinate-c", "succinate", ("h10c4sialo11", "h4o2"), ("h10c4sialo11+h2o", "h2o"), -16.1),
    Reaction("ch4-succinate-d", "succinate", ("h7sialo7", "h4c4o4"), ("h11c4sialo11",), -72.8),
    Reaction("ch4-succinate-e", "succinate", ("h7sialo7+h2o", "h4c4o4"), ("h11c4sialo11+h2o",), -85.9),
    Reaction("ch4-succinate-f", "succinate", ("h11c4sialo11", "h4o2"), ("h11c4sialo11+h2o", "h2o"), -24.9),
    Reaction("ch4-dimer-hydrolysis-a", "unliganded", ("h6sialo7+h2o",), ("h4sio4", "h4alo4"), 35.8),
    Reaction("ch4-dimer-hydrolysis-b", "unliganded", ("h6sialo7+h2o", "h2o"), ("h4sio4", "h3alo3", "h3o2"), 109.4),
    Reaction("ch4-dimer-hydrolysis-c", "unliganded", ("h7sialo7+h2o", "h2o"), ("h4sio4", "h4alo4", "h3o"), 183.6),
    Reaction("ch4-dimer-hydrolysis-d", "unliganded", ("h7sialo7+h2o", "h2o"), ("h4sio4", "h3alo3", "h4o2"), 42.6),
    Reaction("ch4-oxalate-hydrolysis-a", "oxalate", ("h6c2sialo11+h2o",), ("h4sio4", "h4c2alo8"), 99.4),
    Reaction("ch4-oxalate-hydrolysis-b", "oxalate", ("h6c2sialo11+h2o", "h2o"), ("h4sio4", "h3c2alo7", "h3o2"), -59.2),
    Reaction("ch4-oxalate-hydrolysis-c", "oxalate", ("h7c2sialo11+h2o", "h2o"), ("h4sio4", "h3c2alo7", "h4o2"), 62.7),
    Reaction("ch4-oxalate-hydrolysis-d", "oxalate", ("h7c2sialo11+h2o", "h2o"), ("h4sio4", "h4c2alo8", "h3o"), 438.0),
    Reaction("ch4-malonate-hydrolysis-a", "malonate", ("h8c3sialo11+h2o",), ("h4sio4", "h6c3alo8"), 38.6),
    Reaction("ch4-malonate-hydrolysis-b", "malonate", ("h8c3sialo11+h2o", "h2o"), ("h4sio4", "h5c3alo7", "h3o2"), -56.0),
    Reaction("ch4-malonate-hydrolysis-c", "malonate", ("h9c3sialo11+h2o", "h2o"), ("h4sio4", "h5c3alo7", "h4o2"), 57.3),
    Reaction("ch4-malonate-hydrolysis-d", "malonate", ("h9c3sialo11+h2o", "h2o"), ("h4sio4", "h6c3alo8", "h3o"), 368.5),
    Reaction("ch4-succinate-hydrolysis-a", "succinate", ("h10c4sialo11+h2o",), ("h4sio4", "h8c4alo8"), -60.6),
    Reaction("ch4-succinate-hydrolysis-b", "succinate", ("h10c4sialo11+h2o", "h2o"), ("h4sio4", "h7c4alo7", "h3o2"), -66.6),
    Reaction("ch4-succinate-hydrolysis-c", "succinate", ("h11c4sialo11+h2o", "h2o"), ("h4sio4", "h7c4alo7", "h4o2"), 50.9),
    Reaction("ch4-succinate-hydrolysis-d", "succinate", ("h11c4sialo11+h2o", "h2o"), ("h4sio4", "h8c4alo8", "h3o"), 273.5),
    Reaction("ch4-overall-water", "unliganded", ("h6sialo7", "h4o2"), ("h4sio4", "h4alo4", "h2o"), 29.9),
    Reaction("ch4-overall-hydronium", "unliganded", ("h6sialo7", "h5o2"), ("h4sio4", "h3alo3", "h4o2"), -77.8),
    Reaction("ch4-overall-oxalate-water", "oxalate", ("h6sialo7", "c2o4", "h4o2"), ("h4sio4", "h3c2alo7", "h3o2"), -3.4),
    Reaction("ch4-overall-oxalate-hydronium", "oxalate", ("h6sialo7", "c2o4", "h5o2"), ("h4sio4", "h3c2alo7", "h4o2"), -184.7),
    Reaction("ch4-overall-malonate-water", "malonate", ("h6sialo7", "h2c3o4", "h4o2"), ("h4sio4", "h5c3alo7", "h3o2"), 16.0),
    Reaction("ch4-overall-malonate-hydronium", "malonate", ("h6sialo7", "h2c3o4", "h5o2"), ("h4sio4", "h5c3alo7", "h4o2"), -165.4),
    Reaction("ch4-overall-succinate-water", "succinate", ("h6sialo7", "h4c4o4", "h4o2"), ("h4sio4", "h7c4alo7", "h3o2"), 25.9),
    Reaction("ch4-overall-succinate-hydronium", "succinate", ("h6sialo7", "h4c4o4", "h5o2"), ("h4sio4", "h7c4alo7", "h4o2"), -155.5),
)


def _number(value: str) -> float:
    return float(value.replace("D", "E").replace("d", "e"))


def split_link1(text: str) -> list[str]:
    """Split explicit Link1 jobs while retaining single-job Gaussian 94 logs."""
    markers = list(re.finditer(r"(?m)^\s*(?:--Link1--|Entering Link 1\s*=.*)$", text))
    if not markers:
        return [text]
    chunks: list[str] = []
    start = 0
    for marker in markers:
        chunk = text[start : marker.start()]
        if chunk.strip():
            chunks.append(chunk)
        start = marker.start()
    if text[start:].strip():
        chunks.append(text[start:])
    return chunks


def parse_route(section: str) -> tuple[str, str]:
    match = ROUTE_RE.search(section)
    if not match:
        return "", ""
    route = match.group(1)
    method_basis = METHOD_BASIS_RE.search(route)
    if not method_basis:
        return "", ""
    return method_basis.group(1).upper(), method_basis.group(2)


def parse_job(index: int, section: str) -> Job:
    scf = [_number(value) for value in SCF_RE.findall(section)]
    zpe = [_number(value) for value in ZPE_RE.findall(section)]
    frequencies: list[float] = []
    for raw in FREQ_RE.findall(section):
        for token in raw.split():
            try:
                value = _number(token)
            except ValueError:
                continue
            if value < 0:
                frequencies.append(value)
    method, basis = parse_route(section)
    return Job(
        index=index,
        normal="Normal termination" in section,
        electronic=scf[-1] if scf else None,
        electronic_zpe=zpe[-1] if zpe else None,
        method=method,
        basis=basis,
        imaginary=tuple(sorted(frequencies)),
    )


def parse_log(path: Path, root: Path) -> LogResult:
    text = path.read_text(encoding="latin-1")
    jobs = tuple(parse_job(index, section) for index, section in enumerate(split_link1(text), 1))
    eligible = [job for job in jobs if job.normal and job.energy is not None]
    return LogResult(path.relative_to(root).as_posix(), jobs, eligible[-1] if eligible else None)


def choose_primary(logs: Iterable[LogResult]) -> LogResult | None:
    candidates = [log for log in logs if log.selected]
    if not candidates:
        return None

    def rank(log: LogResult) -> tuple[int, int, str]:
        name = Path(log.relpath).name.lower()
        selected = cast(Job, log.selected)
        return (int("freq" in name), int(selected.electronic_zpe is not None), name)

    return max(candidates, key=rank)


def load_species(source: Path) -> dict[str, Species]:
    results_root = source / "dft" / "results"
    grouped: dict[str, list[LogResult]] = {}
    for path in sorted(results_root.rglob("*.log")):
        rel = path.relative_to(results_root)
        if len(rel.parts) < 2:
            continue
        grouped.setdefault(rel.parts[0], []).append(parse_log(path, source))
    return {
        system_id: Species(system_id, choose_primary(logs), tuple(logs))
        for system_id, logs in sorted(grouped.items())
    }


def ligand_for(system_id: str) -> str:
    if "c2" in system_id:
        return "oxalate"
    if "c3" in system_id:
        return "malonate"
    if "c4" in system_id:
        return "succinate"
    if system_id.startswith("h") and "sial" not in system_id and "al" not in system_id:
        return "water"
    return "unliganded"


def fmt(value: float | None, places: int = 12) -> str:
    return "" if value is None else f"{value:.{places}f}"


def source_row(species: Species) -> dict[str, str]:
    primary = species.primary
    if primary is None or primary.selected is None:
        statuses = ["normal" if any(job.normal for job in log.jobs) else "truncated" for log in species.logs]
        return dict.fromkeys(FIELDS, "") | {
            "system_id": species.system_id,
            "reaction_id": "source-calculation",
            "ligand_class": ligand_for(species.system_id),
            "normal_termination": "false",
            "parse_status": "parse_error" if "normal" in statuses else "truncated",
            "notes": "No normally terminated job with a parseable energy; see source-log coverage.",
        }
    job = primary.selected
    all_logs = ";".join(log.relpath for log in species.logs)
    return dict.fromkeys(FIELDS, "") | {
        "system_id": species.system_id,
        "reaction_id": "source-calculation",
        "ligand_class": ligand_for(species.system_id),
        "ref_log": primary.relpath,
        "method": job.method,
        "basis": job.basis,
        "energy_treatment": job.treatment,
        "e_ref_hartree": fmt(job.energy),
        "imaginary_frequencies_cm_1": ";".join(fmt(v, 3) for v in job.imaginary),
        "normal_termination": "true",
        "parse_status": "unmatched",
        "notes": f"Unmatched source calculation; {len(species.logs)} log(s) parsed: {all_logs}.",
    }


def reaction_row(reaction: Reaction, species: dict[str, Species]) -> dict[str, str]:
    ids = reaction.refs + reaction.products
    missing = [sid for sid in ids if sid not in species or not species[sid].primary]
    base = dict.fromkeys(FIELDS, "") | {
        "system_id": reaction.reaction_id,
        "reaction_id": reaction.reaction_id,
        "ligand_class": reaction.ligand_class,
        "normal_termination": "false" if missing else "true",
    }
    if missing:
        base["parse_status"] = "parse_error"
        base["notes"] = "Missing normally terminated source species: " + ", ".join(missing)
        return base

    refs = [cast(LogResult, species[sid].primary) for sid in reaction.refs]
    products = [cast(LogResult, species[sid].primary) for sid in reaction.products]
    jobs = [cast(Job, log.selected) for log in refs + products]
    assert all(job.energy is not None for job in jobs)
    e_ref = sum(cast(float, cast(Job, log.selected).energy) for log in refs)
    e_product = sum(cast(float, cast(Job, log.selected).energy) for log in products)
    methods = sorted({job.method for job in jobs if job.method})
    bases = sorted({job.basis for job in jobs if job.basis})
    treatments = {job.treatment for job in jobs}
    if len(treatments) != 1:
        base.update(
            {
                "ref_log": ";".join(log.relpath for log in refs),
                "product_log": ";".join(log.relpath for log in products),
                "method": ";".join(methods),
                "basis": ";".join(bases),
                "parse_status": "ambiguous",
                "notes": "Mixed electronic and electronic+ZPE constituents; derived energies suppressed.",
            }
        )
        return base
    treatment = treatments.pop()
    delta = (e_product - e_ref) * HARTREE_TO_KJ_MOL
    observed_kcal = delta / 4.184
    discrepancy = observed_kcal - reaction.thesis_kcal_mol
    base.update(
        {
            "ref_log": ";".join(log.relpath for log in refs),
            "product_log": ";".join(log.relpath for log in products),
            "method": ";".join(methods),
            "basis": ";".join(bases),
            "energy_treatment": treatment,
            "e_ref_hartree": fmt(e_ref),
            "e_product_hartree": fmt(e_product),
            "reaction_energy_kj_mol": fmt(delta, 3),
            "parse_status": "unmatched",
            "notes": (
                "Chapter 4 reaction-energy candidate; no TS log exists, so barriers are blank. "
                f"Thesis={reaction.thesis_kcal_mol:.1f} kcal/mol; parsed={observed_kcal:.1f} "
                f"kcal/mol; delta={discrepancy:+.2f} kcal/mol."
            ),
        }
    )
    return base


def coverage_rows(species: dict[str, Species]) -> list[dict[str, str]]:
    rows = []
    for item in species.values():
        for log in item.logs:
            selected = log.selected
            rows.append(
                {
                    "log": log.relpath,
                    "system_id": item.system_id,
                    "link1_jobs": str(len(log.jobs)),
                    "normal_jobs": str(sum(job.normal for job in log.jobs)),
                    "selected_job": "" if selected is None else str(selected.index),
                    "selected_energy_hartree": "" if selected is None else fmt(selected.energy),
                    "energy_treatment": "" if selected is None else selected.treatment,
                    "method": "" if selected is None else selected.method,
                    "basis": "" if selected is None else selected.basis,
                    "imaginary_frequencies_cm_1": "" if selected is None else ";".join(fmt(v, 3) for v in selected.imaginary),
                    "status": "selected" if selected else "no-normal-energy",
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path: Path, species: dict[str, Species], rows: list[dict[str, str]], coverage: list[dict[str, str]]) -> None:
    reactions = [row for row in rows if row["reaction_id"] != "source-calculation"]
    missing = [row for row in coverage if row["status"] != "selected"]
    imaginary = [row for row in coverage if row["imaginary_frequencies_cm_1"]]
    discrepancies = []
    for row in reactions:
        match = re.search(r"delta=([+-]?\d+\.\d+)", row["notes"])
        if match:
            discrepancies.append(abs(float(match.group(1))))
    max_delta = max(discrepancies)
    lines = [
        "# 1999 Gaussian DFT ledger",
        "",
        "> Generated by `tools/parse_gaussian_archive.py`; do not hand-edit the CSVs.",
        "",
        "## Result",
        "",
        f"- Parsed **{len(coverage)}** Gaussian `.log` files in **{len(species)}** system directories.",
        f"- Built **{len(reactions)}** Chapter 4 reaction-energy candidates and **{len(species)}** source-calculation rows.",
        f"- **No transition-state calculation is identifiable in the archive.** Therefore every barrier field is intentionally blank and reaction rows are `unmatched`, not `ok`.",
        f"- {len(missing)} logs have no normally terminated job with a parseable energy; they remain visible in `DFT-SOURCE-LOGS.csv`.",
        f"- {len(imaginary)} selected jobs report one or more negative frequencies. These are recorded as diagnostics, never promoted to transition states by filename guesswork.",
        f"- The {len(reactions)} reproduced Chapter 4 energies agree with the rounded thesis tables within **{max_delta:.2f} kcal/mol**.",
        "",
        "## Files and semantics",
        "",
        "- `DFT-LEDGER.csv` has the fixed board-card schema. Semicolon-delimited `ref_log` or `product_log` values represent stoichiometric sums of species calculations.",
        "- `DFT-SOURCE-LOGS.csv` is the exhaustive audit edge: one row per source log, including failed/truncated calculations.",
        "- `e_ref_hartree` and `e_product_hartree` are stoichiometric sums. `energy_treatment=electronic+ZPE` only when every selected constituent has a Gaussian `Sum of electronic and zero-point Energies` value.",
        "- `normal_termination=true` on a reaction row means every selected constituent job terminated normally; it does **not** imply that a TS exists.",
        "- A log is split at explicit `--Link1--`/`Entering Link 1` boundaries. The last normally terminated section containing an energy is selected. Other sections stay counted in the source-log CSV.",
        "",
        "## Status meanings",
        "",
        "- `ok`: complete reactant/TS/product triplet (none found).",
        "- `unmatched`: parseable source or reaction-energy mapping with no defensible TS.",
        "- `truncated`: no normal termination in the source log.",
        "- `ambiguous`: multiple defensible mappings without a unique choice.",
        "- `parse_error`: expected numerical evidence could not be extracted.",
        "",
        "## Reaction-energy mapping",
        "",
        "The reaction equations and comparison values come from `text/ch4.tex`, Tables 4.4, 4.5, 4.7, 4.8, and 4.10–4.14 in the read-only archive. They cover water adsorption, oxalate/malonate/succinate adsorption, ligand-assisted and unliganded hydrolysis, and overall reactions. This is a reaction-*energy* ledger, not a kinetic barrier data set: Chapter 4 explicitly says transition states were not computed.",
        "",
        "| reaction_id | class | ΔE (kJ/mol) | status |",
        "|---|---:|---:|---|",
    ]
    for row in reactions:
        lines.append(f"| {row['reaction_id']} | {row['ligand_class']} | {row['reaction_energy_kj_mol']} | {row['parse_status']} |")
    lines.extend(
        [
            "",
            "## Promotion rule",
            "",
            "A future row may become `ok` only after a human-auditable mapping names reactant, transition-state, and product logs; the selected TS job terminates normally; the TS has exactly one imaginary frequency; and the displacement follows the intended bond-making/breaking coordinate. A negative frequency alone is not a TS proof.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def validate(rows: list[dict[str, str]], coverage: list[dict[str, str]], source: Path) -> None:
    assert rows, "ledger is empty"
    assert coverage, "coverage ledger is empty"
    expected_logs = len(list((source / "dft" / "results").rglob("*.log")))
    assert len(coverage) == expected_logs, (len(coverage), expected_logs)
    assert [*rows[0]] == FIELDS
    for row in rows:
        assert set(row) == set(FIELDS)
        assert row["parse_status"] in VALID_STATUS
        if row["energy_treatment"]:
            assert row["energy_treatment"] in VALID_TREATMENT
        assert row["normal_termination"] in {"true", "false"}
        imaginary = [float(value) for value in row["imaginary_frequencies_cm_1"].split(";") if value]
        assert imaginary == sorted(imaginary)
        if row["parse_status"] == "ok":
            assert row["ts_log"] and row["barrier_forward_kj_mol"]
    reaction_rows = [row for row in rows if row["reaction_id"].startswith("ch4-")]
    assert len(reaction_rows) == len(REACTIONS)
    assert all(not row["ts_log"] for row in reaction_rows)
    # Independent numerical check against the rounded thesis tables.
    for row, reaction in zip(reaction_rows, REACTIONS, strict=True):
        parsed_kcal = float(row["reaction_energy_kj_mol"]) / 4.184
        # Most rows round-trip to 0.1 kcal/mol. The archived files reproduce
        # Table 4.10(c) as 185.6 rather than the printed 183.6; retain and flag
        # that evidence instead of forcing the published number.
        assert abs(parsed_kcal - reaction.thesis_kcal_mol) <= 2.1, (reaction.reaction_id, parsed_kcal)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True, help="read-only thesis archive root")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--check", action="store_true", help="validate and compare existing outputs")
    args = parser.parse_args()

    species = load_species(args.source)
    rows = [source_row(item) for item in species.values()]
    rows.extend(reaction_row(reaction, species) for reaction in REACTIONS)
    coverage = coverage_rows(species)
    validate(rows, coverage, args.source)

    ledger_path = args.out_dir / "DFT-LEDGER.csv"
    coverage_path = args.out_dir / "DFT-SOURCE-LOGS.csv"
    markdown_path = args.out_dir / "DFT-LEDGER.md"
    if args.check:
        import tempfile

        with tempfile.TemporaryDirectory(prefix="dft-ledger-check-") as temp:
            tempdir = Path(temp)
            write_csv(tempdir / ledger_path.name, rows, FIELDS)
            write_csv(tempdir / coverage_path.name, coverage, list(coverage[0]))
            write_markdown(tempdir / markdown_path.name, species, rows, coverage)
            for existing in (ledger_path, coverage_path, markdown_path):
                candidate = tempdir / existing.name
                assert existing.read_bytes() == candidate.read_bytes(), f"stale generated file: {existing}"
    else:
        write_csv(ledger_path, rows, FIELDS)
        write_csv(coverage_path, coverage, list(coverage[0]))
        write_markdown(markdown_path, species, rows, coverage)

    print(f"validated {len(rows)} ledger rows and {len(coverage)} source logs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
