#!/usr/bin/env python3
"""Materialize and verify the A9b pH 3-5 reservoir/origin contract.

The A9 ``kaolinite-approx.toml`` deck is retained as historical evidence. This
runner turns it into an executable deck for one declared open-flow profile and
writes a hash-bound JSON receipt naming the complete three-profile contract.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

import tomllib

PH_VALUES = (3, 4, 5)
SCHEMA_VERSION = 1
LEGACY_SINK = 1.0e-30
SAFE_NAME = re.compile(r"^[a-z0-9][a-z0-9._-]{0,63}$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
CONTRACT_KEYS = {
    "schema_version",
    "name",
    "base_deck_sha256",
    "temperature_k",
    "boundary",
    "product_activity",
    "product_activity_unit",
    "product_activity_basis",
    "acid_rate_order",
    "acid_rate_anchor_ph",
    "acid_scaled_reactions",
    "al_species",
    "si_species",
    "approximations",
    "sources",
    "profiles",
}
THERMO_BLOCK = re.compile(
    r"(?ms)^# All temperatures are Kelvin.*?^# --- build-time passes"
)


@dataclass(frozen=True)
class ReservoirProfile:
    ph: int
    activities: dict[str, float]
    al_species: str
    si_species: str


@dataclass(frozen=True)
class ReservoirContract:
    name: str
    base_deck_sha256: str
    temperature_k: float
    boundary: str
    product_activity: float
    product_activity_unit: str
    product_activity_basis: str
    acid_rate_order: float
    acid_rate_anchor_ph: int
    acid_scaled_reactions: tuple[str, ...]
    al_species: str
    si_species: str
    approximations: tuple[str, ...]
    sources: tuple[dict[str, str], ...]
    profiles: dict[int, ReservoirProfile]


def _finite_positive(value: object, name: str) -> float:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise TypeError(f"{name} must be numeric")
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return number


def load_contract(path: Path) -> ReservoirContract:
    parsed = tomllib.loads(path.read_text(encoding="utf-8"))
    if set(parsed) != CONTRACT_KEYS:
        raise ValueError(
            f"reservoir contract keys must be exactly {sorted(CONTRACT_KEYS)}"
        )
    if parsed.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"reservoir schema_version must be {SCHEMA_VERSION}")
    name = parsed.get("name")
    if not isinstance(name, str) or SAFE_NAME.fullmatch(name) is None:
        raise ValueError("reservoir name must be a safe lower-case identifier")
    base_deck_sha256 = parsed.get("base_deck_sha256")
    if (
        not isinstance(base_deck_sha256, str)
        or SHA256.fullmatch(base_deck_sha256) is None
    ):
        raise ValueError("base_deck_sha256 must be a lower-case SHA-256 digest")
    temperature = _finite_positive(parsed.get("temperature_k"), "temperature_k")
    if temperature != 298.0:
        raise ValueError("reservoir temperature_k must be exactly 298.0")
    boundary = parsed.get("boundary")
    if boundary != "open_flow_constant_activity":
        raise ValueError("reservoir boundary must be open_flow_constant_activity")
    product = _finite_positive(parsed.get("product_activity"), "product_activity")
    if product <= 1.0e-20:
        raise ValueError(
            "a realistic reservoir contract cannot reuse the 1e-30 numerical sink"
        )
    if product != 1.0e-6:
        raise ValueError(
            "product_activity must be the declared 1e-6 trace-inflow value"
        )
    product_unit = parsed.get("product_activity_unit")
    product_basis = parsed.get("product_activity_basis")
    if not isinstance(product_unit, str) or not product_unit.strip():
        raise ValueError("product_activity_unit must be a non-empty string")
    if not isinstance(product_basis, str) or not product_basis.strip():
        raise ValueError("product_activity_basis must be a non-empty string")
    acid_rate_order = _finite_positive(parsed.get("acid_rate_order"), "acid_rate_order")
    acid_rate_anchor_ph = parsed.get("acid_rate_anchor_ph")
    acid_scaled_reactions = parsed.get("acid_scaled_reactions")
    if acid_rate_order != 0.777 or acid_rate_anchor_ph != 4:
        raise ValueError("acid rate contract must use P&K n_H=0.777 anchored at pH 4")
    if not isinstance(acid_scaled_reactions, list) or acid_scaled_reactions != [
        "desorb-al",
        "desorb-si",
    ]:
        raise ValueError("acid scaling must target exactly the terminal Al/Si releases")

    al_species = parsed.get("al_species")
    si_species = parsed.get("si_species")
    if (
        al_species != "Al3+ (total-Al proxy)"
        or si_species != "H4SiO4(aq) (total-Si proxy)"
    ):
        raise ValueError(
            "reservoir aqueous species labels do not match the v1 contract"
        )
    approximations = tuple(parsed.get("approximations", ()))
    if len(approximations) < 6 or not all(
        isinstance(item, str) and item for item in approximations
    ):
        raise ValueError("reservoir contract must state every v1 approximation")
    raw_sources = parsed.get("sources")
    if not isinstance(raw_sources, list) or len(raw_sources) < 3:
        raise ValueError(
            "reservoir contract requires the three cited laboratory/ledger sources"
        )
    sources: list[dict[str, str]] = []
    for source in raw_sources:
        if not isinstance(source, dict) or set(source) != {
            "id",
            "citation",
            "url",
            "use",
        }:
            raise ValueError("each source needs exactly id, citation, url, and use")
        if not all(isinstance(value, str) and value for value in source.values()):
            raise ValueError("source fields must be non-empty strings")
        sources.append(dict(source))

    raw_profiles = parsed.get("profiles")
    if not isinstance(raw_profiles, list):
        raise TypeError("reservoir profiles must be a list")
    profiles: dict[int, ReservoirProfile] = {}
    for raw in raw_profiles:
        if not isinstance(raw, dict) or set(raw) != {"ph", "H_plus", "Al", "Si"}:
            raise ValueError("each profile needs exactly ph, H_plus, Al, and Si")
        ph = raw["ph"]
        if type(ph) is not int or ph not in PH_VALUES or ph in profiles:
            raise ValueError(
                "reservoir profiles must uniquely cover integer pH 3, 4, and 5"
            )
        h_plus = _finite_positive(raw["H_plus"], f"pH {ph} H_plus")
        al = _finite_positive(raw["Al"], f"pH {ph} Al")
        si = _finite_positive(raw["Si"], f"pH {ph} Si")
        if h_plus != 10.0 ** (-ph):
            raise ValueError(f"pH {ph} must satisfy a_H+ = 10^(-pH) exactly")
        if al != product or si != product:
            raise ValueError(f"pH {ph} Al/Si activities must equal product_activity")
        profiles[ph] = ReservoirProfile(
            ph=ph,
            activities={"H_plus": h_plus, "Al": al, "Si": si},
            al_species=al_species,
            si_species=si_species,
        )
    if set(profiles) != set(PH_VALUES):
        raise ValueError("reservoir profiles must cover pH 3, 4, and 5")

    return ReservoirContract(
        name=name,
        base_deck_sha256=base_deck_sha256,
        temperature_k=temperature,
        boundary=boundary,
        product_activity=product,
        product_activity_unit=product_unit,
        product_activity_basis=product_basis,
        acid_rate_order=acid_rate_order,
        acid_rate_anchor_ph=acid_rate_anchor_ph,
        acid_scaled_reactions=tuple(acid_scaled_reactions),
        al_species=al_species,
        si_species=si_species,
        approximations=approximations,
        sources=tuple(sources),
        profiles=profiles,
    )


def _write_atomic(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def _profile_payload(profile: ReservoirProfile) -> dict[str, object]:
    return {
        "ph": profile.ph,
        "activities": profile.activities,
        "aqueous_species": {
            "Al": profile.al_species,
            "Si": profile.si_species,
            "H_plus": "H+(aq)",
        },
    }


def verify_evidence(deck_path: Path, evidence_path: Path) -> dict[str, object]:
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    if not isinstance(evidence, dict):
        raise TypeError("reservoir evidence must be a JSON object")
    actual_hash = hashlib.sha256(deck_path.read_bytes()).hexdigest()
    if evidence.get("generated_deck_sha256") != actual_hash:
        raise ValueError("reservoir evidence does not match the generated deck")
    if evidence.get("generated_deck") != deck_path.name:
        raise ValueError("reservoir evidence generated_deck name does not match")
    return evidence


def materialize(
    base_deck: Path,
    contract_path: Path,
    ph: int,
    out_deck: Path,
    evidence_path: Path,
) -> dict[str, object]:
    reservoir = load_contract(contract_path)
    if ph not in reservoir.profiles:
        raise ValueError("selected pH must be one of 3, 4, or 5")
    if out_deck.resolve() == evidence_path.resolve():
        raise ValueError("generated deck and evidence paths must differ")
    profile = reservoir.profiles[ph]
    source = base_deck.read_text(encoding="utf-8")
    source_hash = hashlib.sha256(source.encode("utf-8")).hexdigest()
    if source_hash != reservoir.base_deck_sha256:
        raise ValueError(
            "base deck SHA-256 does not match the contract; review and rebind explicitly"
        )
    parsed_source = tomllib.loads(source)
    if parsed_source.get("thermo", {}).get("temperature") != reservoir.temperature_k:
        raise ValueError("base deck temperature does not match the reservoir contract")
    old_activity = parsed_source.get("thermo", {}).get("activity")
    if old_activity != {"Al": LEGACY_SINK, "Si": LEGACY_SINK}:
        raise ValueError("base deck must be the preserved A9 1e-30 historical sink")

    replacement = f'''# --- A9b reservoir contract (generated; do not hand-edit) -----------------
# contract = "{reservoir.name}"; pH = {ph}; boundary = "{reservoir.boundary}"
# Al means {reservoir.al_species}; Si means {reservoir.si_species}.
# H_plus records pH; P&K n_H scales terminal Si/Al release barriers from pH 4.
[thermo]
temperature = {reservoir.temperature_k:.1f}

[thermo.mu]
Al = -1.0
Si = -1.0

[thermo.activity]
Al = {profile.activities["Al"]:.1e}
Si = {profile.activities["Si"]:.1e}
H_plus = {profile.activities["H_plus"]:.1e}

# --- build-time passes'''
    generated, count = THERMO_BLOCK.subn(replacement, source, count=1)
    if count != 1:
        raise ValueError("base deck does not contain the canonical A9 thermo block")
    generated = generated.replace(
        "# deck has no H+/OH- coupling and therefore does not encode pH catalysis.",
        "# generated A9b decks apply the documented P&K H+ order to terminal release.",
    )
    barrier_delta = (
        0.00198720425864083
        * reservoir.temperature_k
        * math.log(10.0)
        * reservoir.acid_rate_order
        * (ph - reservoir.acid_rate_anchor_ph)
    )
    scaled_barriers: dict[str, dict[str, float]] = {}
    for reaction_name in reservoir.acid_scaled_reactions:
        reaction = next(
            (
                item
                for item in parsed_source["reactions"]
                if item.get("name") == reaction_name
            ),
            None,
        )
        if reaction is None:
            raise ValueError(
                f"base deck is missing acid-scaled reaction {reaction_name}"
            )
        base_ea = float(reaction["rate"]["arrhenius"]["ea"])
        selected_ea = base_ea + barrier_delta
        pattern = re.compile(
            rf'(?ms)(\[\[reactions\]\]\nname = "{re.escape(reaction_name)}"\n'
            rf"(?:(?!\[\[reactions\]\]).)*?rate = \{{ arrhenius = \{{ prefactor = 1\.0e13, ea = )"
            rf"{re.escape(str(base_ea))}( \}} \}})"
        )
        generated, replacements = pattern.subn(
            rf"\g<1>{selected_ea:.9f}\g<2>", generated, count=1
        )
        if replacements != 1:
            raise ValueError(f"could not bind acid scaling to {reaction_name}")
        scaled_barriers[reaction_name] = {
            "anchor_ea_kcal_mol": base_ea,
            "selected_ea_kcal_mol": selected_ea,
        }
    parsed_generated = tomllib.loads(generated)
    if parsed_generated["thermo"]["activity"] != {
        "Al": profile.activities["Al"],
        "Si": profile.activities["Si"],
        "H_plus": profile.activities["H_plus"],
    }:
        raise AssertionError("generated deck activity map drifted")
    if "1.0e-30" in generated or "1e-30" in generated:
        raise AssertionError(
            "generated realistic reservoir deck retained the numerical sink"
        )
    _write_atomic(out_deck, generated)
    deck_hash = hashlib.sha256(out_deck.read_bytes()).hexdigest()
    evidence: dict[str, object] = {
        "schema": "petra-a9b-reservoir-origin-v1",
        "contract": reservoir.name,
        "temperature_k": reservoir.temperature_k,
        "boundary": reservoir.boundary,
        "selected_profile": _profile_payload(profile),
        "profiles": [
            _profile_payload(reservoir.profiles[value]) for value in PH_VALUES
        ],
        "product_activity_basis": reservoir.product_activity_basis,
        "product_activity_unit": reservoir.product_activity_unit,
        "acid_rate_scaling": {
            "source": "Palandri and Kharaka 2004, kaolinite acid mechanism",
            "order_H_plus": reservoir.acid_rate_order,
            "anchor_ph": reservoir.acid_rate_anchor_ph,
            "barrier_delta_kcal_mol": barrier_delta,
            "reactions": scaled_barriers,
        },
        "approximations": list(reservoir.approximations),
        "sources": list(reservoir.sources),
        "origin_accounting": {
            "initial_occupied_cation": "original_lattice",
            "adsorbed_cation": "reservoir",
            "surface_state_transition": "preserve origin",
            "desorption_numerator": "increment only for original_lattice",
            "reservoir_desorption": "gross release only; never physical dissolution",
        },
        "source_deck": base_deck.name,
        "source_deck_sha256": source_hash,
        "contract_path": contract_path.name,
        "contract_sha256": hashlib.sha256(contract_path.read_bytes()).hexdigest(),
        "generated_deck": out_deck.name,
        "generated_deck_sha256": deck_hash,
    }
    _write_atomic(evidence_path, json.dumps(evidence, indent=2, sort_keys=True) + "\n")
    verify_evidence(out_deck, evidence_path)
    return evidence


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    validate = subparsers.add_parser("validate")
    validate.add_argument("contract", type=Path)
    build = subparsers.add_parser("materialize")
    build.add_argument("base_deck", type=Path)
    build.add_argument("contract", type=Path)
    build.add_argument("ph", type=int, choices=PH_VALUES)
    build.add_argument("out_deck", type=Path)
    build.add_argument("evidence", type=Path)
    verify = subparsers.add_parser("verify")
    verify.add_argument("deck", type=Path)
    verify.add_argument("evidence", type=Path)
    args = parser.parse_args(argv)
    if args.command == "validate":
        loaded = load_contract(args.contract)
        print(json.dumps({"contract": loaded.name, "profiles": list(loaded.profiles)}))
        return 0
    if args.command == "verify":
        evidence = verify_evidence(args.deck, args.evidence)
        print(json.dumps(evidence["selected_profile"], sort_keys=True))
        return 0
    evidence = materialize(
        args.base_deck, args.contract, args.ph, args.out_deck, args.evidence
    )
    print(json.dumps(evidence["selected_profile"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
