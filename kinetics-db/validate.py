#!/usr/bin/env python3
"""Validate kinetics-db v0 mineral TOML records using only the stdlib."""

from __future__ import annotations

import math
import sys
import tomllib
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
MINERALS = ROOT / "minerals"
AREA_BASES = {"BET", "geometric", "unspecified"}
MECHANISM_KINDS = {"acid", "neutral", "base", "other"}
MIN_INVENTORY_SIZE = 25
TOP_LEVEL = {
    "schema_version": int,
    "name": str,
    "formula": str,
    "mineral_group": str,
    "rate_law": str,
    "surface_area_basis": str,
    "caveats": list,
    "source": dict,
    "mechanism": list,
}
SOURCE_FIELDS = {
    "id": str,
    "citation": str,
    "url": str,
    "year": int,
    "page": int,
    "table": int,
    "notes": str,
}
MECHANISM_FIELDS = {
    "kind": str,
    "log_k25_mol_m2_s": (int, float),
    "activation_energy_kj_mol": (int, float),
    "reaction_orders": dict,
    "source_ref": str,
    "uncertainty_status": str,
    "uncertainty_log_k25": list,
    "uncertainty_activation_energy_kj_mol": list,
    "uncertainty_reaction_orders": list,
    "uncertainty_notes": list,
}


def finite_number(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(value)
    )


def nonempty_string(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def validate_file(path: Path) -> tuple[list[str], dict[str, Any] | None]:
    errors: list[str] = []
    try:
        with path.open("rb") as handle:
            record = tomllib.load(handle)
    except (OSError, tomllib.TOMLDecodeError) as exc:
        return [f"{path.name}: cannot parse: {exc}"], None

    for field, expected in TOP_LEVEL.items():
        if field not in record:
            errors.append(f"{path.name}: missing top-level field {field}")
        elif not isinstance(record[field], expected):
            errors.append(f"{path.name}: {field} has wrong type")

    if errors:
        return errors, record

    if record["schema_version"] != 0:
        errors.append(f"{path.name}: schema_version must be 0")
    for field in ("name", "formula", "mineral_group", "rate_law"):
        if not nonempty_string(record[field]):
            errors.append(f"{path.name}: {field} must be non-empty")
    if record["surface_area_basis"] not in AREA_BASES:
        errors.append(
            f"{path.name}: unsupported surface_area_basis {record['surface_area_basis']!r}"
        )
    if not all(isinstance(item, str) and item.strip() for item in record["caveats"]):
        errors.append(f"{path.name}: caveats must contain only non-empty strings")

    source = record["source"]
    for field, expected in SOURCE_FIELDS.items():
        if field not in source:
            errors.append(f"{path.name}: missing source.{field}")
        elif not isinstance(source[field], expected):
            errors.append(f"{path.name}: source.{field} has wrong type")
    if all(field in source for field in SOURCE_FIELDS):
        for field in ("id", "citation", "url", "notes"):
            if not nonempty_string(source[field]):
                errors.append(f"{path.name}: source.{field} must be non-empty")
        if not source["url"].startswith("https://"):
            errors.append(f"{path.name}: source.url must use https")
        for field in ("year", "page", "table"):
            if source[field] <= 0:
                errors.append(f"{path.name}: source.{field} must be positive")

    mechanisms = record["mechanism"]
    if not mechanisms:
        errors.append(f"{path.name}: at least one mechanism is required")
    seen_kinds: set[str] = set()
    for index, mechanism in enumerate(mechanisms):
        label = f"{path.name}: mechanism[{index}]"
        if not isinstance(mechanism, dict):
            errors.append(f"{label} must be a table")
            continue
        for field, expected in MECHANISM_FIELDS.items():
            if field not in mechanism:
                errors.append(f"{label} missing {field}")
            elif not isinstance(mechanism[field], expected):
                errors.append(f"{label}.{field} has wrong type")
        if any(field not in mechanism for field in MECHANISM_FIELDS):
            continue
        kind = mechanism["kind"]
        if kind not in MECHANISM_KINDS:
            errors.append(f"{label}.kind is unsupported: {kind!r}")
        if kind in seen_kinds:
            errors.append(f"{label}.kind duplicates {kind!r}")
        seen_kinds.add(kind)
        for field in ("log_k25_mol_m2_s", "activation_energy_kj_mol"):
            if not finite_number(mechanism[field]):
                errors.append(f"{label}.{field} must be a finite number")
        if isinstance(mechanism["reaction_orders"], dict):
            for order, value in mechanism["reaction_orders"].items():
                if not nonempty_string(order) or not finite_number(value):
                    errors.append(
                        f"{label}.reaction_orders must map name {order!r} to a finite number (got {value!r})"
                    )
        if mechanism["source_ref"] != source.get("id"):
            errors.append(f"{label}.source_ref does not match source.id")
        if mechanism["uncertainty_status"] not in {"not_reported", "reported"}:
            errors.append(f"{label}.uncertainty_status is unsupported")
        for field in (
            "uncertainty_log_k25",
            "uncertainty_activation_energy_kj_mol",
            "uncertainty_reaction_orders",
        ):
            values = mechanism[field]
            if len(values) > 1 or not all(finite_number(value) for value in values):
                errors.append(
                    f"{label}.{field} must be empty or contain one finite number"
                )
        if not all(nonempty_string(note) for note in mechanism["uncertainty_notes"]):
            errors.append(f"{label}.uncertainty_notes must contain non-empty strings")
        has_uncertainty = any(
            mechanism[field]
            for field in (
                "uncertainty_log_k25",
                "uncertainty_activation_energy_kj_mol",
                "uncertainty_reaction_orders",
                "uncertainty_notes",
            )
        )
        if mechanism["uncertainty_status"] == "reported" and not has_uncertainty:
            errors.append(f"{label}: reported uncertainty has no values or notes")
        if mechanism["uncertainty_status"] == "not_reported" and has_uncertainty:
            errors.append(
                f"{label}: not_reported uncertainty must use empty placeholders"
            )
    return errors, record


def main() -> int:
    paths = sorted(MINERALS.glob("*.toml"))
    errors: list[str] = []
    if len(paths) < MIN_INVENTORY_SIZE:
        errors.append(
            f"inventory contains {len(paths)} mineral files; "
            f"at least {MIN_INVENTORY_SIZE} required"
        )
    names: dict[str, Path] = {}
    records: list[dict[str, Any]] = []
    for path in paths:
        file_errors, record = validate_file(path)
        errors.extend(file_errors)
        if record is None:
            continue
        records.append(record)
        name = record.get("name")
        if isinstance(name, str):
            key = name.casefold()
            if key in names:
                errors.append(
                    f"duplicate mineral name {name!r}: {names[key].name}, {path.name}"
                )
            names[key] = path
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        print(
            f"FAILED: {len(errors)} error(s) across {len(paths)} file(s)",
            file=sys.stderr,
        )
        return 1
    mechanism_count = sum(len(record["mechanism"]) for record in records)
    print(f"OK: {len(paths)} mineral files, {mechanism_count} mechanisms, schema v0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
