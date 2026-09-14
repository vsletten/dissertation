#!/usr/bin/env python3
"""Prepare and analyze bounded E3b periodic CP2K spot checks.

This harness performs no DFT during ``prepare`` or ``analyze``. It derives three
3x2x1 models from the hash-pinned E3a source, emits PBE-D3 CP2K endpoint and
8-image CI-NEB inputs plus matched-cell classical inputs, and fail-closes any
calibration whose execution, SCF, NEB, atom-identity, or cell-transfer evidence
is incomplete. The ``smoke`` command is the only execution path and is limited
to a caller-bounded one-image ENERGY_FORCE timing run.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import pathlib
import re
import shutil
import subprocess
import sys
import time
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, replace
from typing import Any

import e3a_campaign
import e3a_input_contract as contract

WORKBOOK_SHA256 = contract.nteme_neb.WORKBOOK_SHA256
CAMPAIGN_MANIFEST_SHA256 = (
    "81effdf35f1ce3600ee415c85b979da0b80669ab5d8156fb9aa56239fb1c1621"
)
CAMPAIGN_RESULT_SHA256 = (
    "2d4c13acad030b12955684af66e6abe91db6cc9a821d2ff7beec4da7eaa1f5a7"
)
MODEL_ORDER = (
    "reconstructed-replication",
    "dehydroxylate-lattice",
    "xenon-divacancy",
)
SOURCE_SUPERCELL = (6, 3, 1)
TARGET_SUPERCELL = (3, 2, 1)
SOURCE_ATOMS_PER_UNIT_CELL = 84
TARGET_PRISTINE_COUNTS = Counter(
    {
        "K": 24,
        "Al1": 48,
        "Al2": 24,
        "Si": 72,
        "O1": 48,
        "O2": 144,
        "O3": 96,
        "H": 48,
    }
)
EXPECTED_MODEL_COUNTS = {
    "reconstructed-replication": Counter(
        {
            "K": 21,
            "Al1": 48,
            "Al2": 21,
            "Si": 75,
            "O1": 48,
            "O2": 156,
            "O3": 84,
            "H": 48,
            "Ar": 1,
        }
    ),
    "dehydroxylate-lattice": Counter(
        {
            "K": 21,
            "Al1": 48,
            "Al2": 21,
            "Si": 75,
            "O1": 46,
            "O2": 157,
            "O3": 84,
            "H": 46,
            "Ar": 1,
        }
    ),
    "xenon-divacancy": Counter(
        {
            "K": 21,
            "Al1": 48,
            "Al2": 21,
            "Si": 75,
            "O1": 48,
            "O2": 156,
            "O3": 84,
            "H": 48,
            "Xe": 1,
        }
    ),
}
CHEMICAL_ELEMENT = {
    "K": "K",
    "Al1": "Al",
    "Al2": "Al",
    "Si": "Si",
    "O1": "O",
    "O2": "O",
    "O3": "O",
    "H": "H",
    "Ar": "Ar",
    "Xe": "Xe",
}
POTENTIAL = {
    "K": "GTH-PBE-q9",
    "Al": "GTH-PBE-q3",
    "Si": "GTH-PBE-q4",
    "O": "GTH-PBE-q6",
    "H": "GTH-PBE-q1",
    "Ar": "GTH-PBE-q8",
    "Xe": "GTH-PBE-q8",
}
VALENCE_ELECTRONS = {
    "K": 9,
    "Al": 3,
    "Si": 4,
    "O": 6,
    "H": 1,
    "Ar": 8,
    "Xe": 8,
}
BASIS_SET = "DZVP-MOLOPT-SR-GTH"
IMAGE_COUNT = 8
MINIMUM_DISTANCE_ANGSTROM = 0.70
AL_O_CUTOFF_ANGSTROM = 2.30
LOCAL_ROUTE_CUTOFF_ANGSTROM = 6.0
DEFAULT_TRANSFER_TOLERANCE_KCAL_MOL = 5.0
DEFAULT_ENDORSEMENT_TOLERANCE_KCAL_MOL = 5.0
ALLOWED_OBSERVATION_STATUS = {
    "converged",
    "incomplete-scf",
    "incomplete-timeout",
    "incomplete-convergence",
    "incomplete-execution",
}
SCIENTIFIC_SCOPE = {
    "reconstructed-replication": (
        "reconstructed neutral Ar route; not exact unpublished-input replication"
    ),
    "dehydroxylate-lattice": (
        "local post-dehydroxylation Ar hop; not a dehydroxylation reaction path"
    ),
    "xenon-divacancy": (
        "Xe divacancy calibration; dry-muscovite Xe remains a DFT spot check"
    ),
}


@dataclass(frozen=True)
class MethodProvenance:
    cp2k_version: str
    cp2k_image_digest: str
    basis_set_file: pathlib.Path
    potential_file: pathlib.Path
    d3_parameter_file: pathlib.Path


@dataclass(frozen=True)
class Reduction:
    grid_origin_fractional: tuple[float, float]
    window_start: tuple[int, int]
    source_ids: tuple[int, ...]


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_sha256(path: pathlib.Path, expected: str, *, label: str) -> str:
    if not path.is_file():
        raise ValueError(f"{label} is not a file: {path}")
    actual = sha256(path)
    if actual != expected:
        raise ValueError(f"{label} hash mismatch: expected {expected}, got {actual}")
    return actual


def read_json(path: pathlib.Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise TypeError(f"{path}: expected a JSON object")
    return value


def write_json(path: pathlib.Path, value: object) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def canonical_sha256(value: object) -> str:
    encoded = json.dumps(value, separators=(",", ":"), sort_keys=True).encode()
    return hashlib.sha256(encoded).hexdigest()


def prepare_output_dir(path: pathlib.Path) -> None:
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise ValueError(f"output path must be a new or empty directory: {path}")
    path.mkdir(parents=True, exist_ok=True)


def verify_e3a_inputs(
    workbook: pathlib.Path, campaign_root: pathlib.Path
) -> dict[str, Any]:
    require_sha256(workbook, WORKBOOK_SHA256, label="E3a source workbook")
    manifest_path = campaign_root / "manifest.json"
    require_sha256(
        manifest_path,
        CAMPAIGN_MANIFEST_SHA256,
        label="E3a campaign manifest",
    )
    manifest = read_json(manifest_path)
    files = manifest.get("files")
    if not isinstance(files, list):
        raise TypeError("E3a campaign manifest has no typed files list")
    result_entries = [
        entry
        for entry in files
        if isinstance(entry, dict) and entry.get("path") == "campaign-result.json"
    ]
    if len(result_entries) != 1:
        raise ValueError("E3a campaign manifest must bind campaign-result.json once")
    entry_hash = result_entries[0].get("sha256")
    if entry_hash != CAMPAIGN_RESULT_SHA256:
        raise ValueError("E3a campaign manifest binds an unexpected campaign result")
    result_path = campaign_root / "campaign-result.json"
    require_sha256(
        result_path,
        CAMPAIGN_RESULT_SHA256,
        label="E3a campaign result",
    )
    campaign = read_json(result_path)
    if campaign.get("schema") != "e3a-classical-neb-campaign-v1":
        raise ValueError("unexpected E3a campaign schema")
    source = campaign.get("source")
    if not isinstance(source, dict) or source.get("source_sha256") != WORKBOOK_SHA256:
        raise ValueError("E3a campaign is not bound to the pinned source workbook")
    models = campaign.get("models")
    if not isinstance(models, dict) or any(name not in models for name in MODEL_ORDER):
        raise ValueError("E3a campaign is missing a required E3b comparison model")
    return campaign


def _fractional_bin(value: float, origin: float, count: int) -> int:
    return math.floor(((value - origin) % 1.0) * count) % count


def _candidate_origins(values: Iterable[float], count: int) -> list[float]:
    width = 1.0 / count
    residues = sorted({value % width for value in values})
    candidates = []
    for index, lower in enumerate(residues):
        upper = residues[(index + 1) % len(residues)]
        if index == len(residues) - 1:
            upper += width
        candidates.append(((lower + upper) / 2.0) % width)
    return candidates


def find_grid_origin(
    atoms: Sequence[contract.nteme_neb.Atom], cell: contract.nteme_neb.Cell
) -> tuple[float, float]:
    fractional = [contract.fractional(atom, cell) for atom in atoms]
    x_candidates = []
    for origin in _candidate_origins((item[0] for item in fractional), 6):
        counts = Counter(_fractional_bin(item[0], origin, 6) for item in fractional)
        if set(counts.values()) == {SOURCE_ATOMS_PER_UNIT_CELL * 3}:
            x_candidates.append(origin)
    y_candidates = []
    for origin in _candidate_origins((item[1] for item in fractional), 3):
        counts = Counter(_fractional_bin(item[1], origin, 3) for item in fractional)
        if set(counts.values()) == {SOURCE_ATOMS_PER_UNIT_CELL * 6}:
            y_candidates.append(origin)
    for x_origin in sorted(x_candidates):
        for y_origin in sorted(y_candidates):
            counts = Counter(
                (
                    _fractional_bin(item[0], x_origin, 6),
                    _fractional_bin(item[1], y_origin, 3),
                )
                for item in fractional
            )
            if len(counts) == 18 and set(counts.values()) == {
                SOURCE_ATOMS_PER_UNIT_CELL
            }:
                return x_origin, y_origin
    raise ValueError("could not recover the source 6x3 unit-cell tiling")


def reduce_pristine(
    atoms: Sequence[contract.nteme_neb.Atom], cell: contract.nteme_neb.Cell
) -> tuple[list[contract.nteme_neb.Atom], contract.nteme_neb.Cell, Reduction]:
    grid_x, grid_y = find_grid_origin(atoms, cell)
    fractional = {atom.id: contract.fractional(atom, cell) for atom in atoms}
    by_tile: dict[tuple[int, int], list[contract.nteme_neb.Atom]] = {}
    for atom in atoms:
        item = fractional[atom.id]
        tile = (
            _fractional_bin(item[0], grid_x, 6),
            _fractional_bin(item[1], grid_y, 3),
        )
        by_tile.setdefault(tile, []).append(atom)

    route_ids = {3, 4, 88}
    candidates: list[tuple[int, int, list[contract.nteme_neb.Atom]]] = []
    for start_x in range(6):
        for start_y in range(3):
            tiles = {
                ((start_x + dx) % 6, (start_y + dy) % 3)
                for dx in range(3)
                for dy in range(2)
            }
            selected = [atom for tile in tiles for atom in by_tile[tile]]
            ids = {atom.id for atom in selected}
            counts = Counter(atom.element for atom in selected)
            if (
                len(selected) == SOURCE_ATOMS_PER_UNIT_CELL * 6
                and route_ids <= ids
                and counts == TARGET_PRISTINE_COUNTS
                and math.isclose(contract.net_charge(selected), 0.0, abs_tol=1.0e-8)
            ):
                candidates.append((start_x, start_y, selected))
    if len(candidates) != 1:
        windows = [(item[0], item[1]) for item in candidates]
        raise ValueError(f"expected one neutral 3x2 route window, found {windows}")
    start_x, start_y, selected = candidates[0]
    target_cell = contract.nteme_neb.Cell(
        a=cell.a * TARGET_SUPERCELL[0] / SOURCE_SUPERCELL[0],
        b=cell.b * TARGET_SUPERCELL[1] / SOURCE_SUPERCELL[1],
        c=cell.c,
        alpha=cell.alpha,
        beta=cell.beta,
        gamma=cell.gamma,
    )
    origin_x = (grid_x + start_x / 6.0) % 1.0
    origin_y = (grid_y + start_y / 3.0) % 1.0
    reduced = []
    for atom in selected:
        source_fractional = fractional[atom.id]
        target_fractional = (
            ((source_fractional[0] - origin_x) % 1.0) / 0.5,
            ((source_fractional[1] - origin_y) % 1.0) / (2.0 / 3.0),
            source_fractional[2] % 1.0,
        )
        if target_fractional[0] >= 1.0 + 1.0e-10 or target_fractional[1] >= (
            1.0 + 1.0e-10
        ):
            raise ValueError(f"atom {atom.id} lies outside the reduced cell")
        x, y, z = contract._vector(target_cell, target_fractional)
        reduced.append(replace(atom, x=x, y=y, z=z))
    reduced.sort(key=lambda atom: atom.id)
    minimum, _ids = contract.minimum_pair(reduced, target_cell)
    if minimum < MINIMUM_DISTANCE_ANGSTROM:
        raise ValueError("3x2 remap created an atomic collision")
    return (
        reduced,
        target_cell,
        Reduction(
            grid_origin_fractional=(grid_x, grid_y),
            window_start=(start_x, start_y),
            source_ids=tuple(atom.id for atom in reduced),
        ),
    )


def compensate_reduced_route(
    atoms: Sequence[contract.nteme_neb.Atom],
    cell: contract.nteme_neb.Cell,
    route_center: contract.nteme_neb.Atom,
    vacancy_count: int,
) -> tuple[list[contract.nteme_neb.Atom], list[dict[str, Any]]]:
    by_id = {atom.id: atom for atom in atoms}
    candidates = []
    for atom in atoms:
        if atom.element != "Al2":
            continue
        oxygen_ids = contract.tetrahedral_oxygen_ids(atom, atoms, cell)
        if not all(by_id[atom_id].element == "O3" for atom_id in oxygen_ids):
            continue
        distance = contract.periodic_distance(atom, route_center, cell)
        if distance >= 8.0:
            candidates.append((distance, atom.id, oxygen_ids, atom))
    selected = []
    selected_oxygen: set[int] = set()
    for candidate in sorted(candidates, key=lambda item: (-item[0], item[1])):
        distance, _atom_id, oxygen_ids, atom = candidate
        if selected_oxygen.intersection(oxygen_ids):
            continue
        if any(
            contract.periodic_distance(atom, prior[3], cell) < 4.5 for prior in selected
        ):
            continue
        selected.append(candidate)
        selected_oxygen.update(oxygen_ids)
        if len(selected) == vacancy_count:
            break
    if len(selected) != vacancy_count:
        raise ValueError(
            f"could not select {vacancy_count} reduced-cell remote compensators"
        )
    selected_ids = {item[1] for item in selected}
    compensated = []
    for atom in atoms:
        if atom.id in selected_ids:
            compensated.append(replace(atom, element="Si"))
        elif atom.id in selected_oxygen:
            compensated.append(replace(atom, element="O2"))
        else:
            compensated.append(atom)
    if not math.isclose(contract.net_charge(compensated), 0.0, abs_tol=1.0e-8):
        raise ValueError("reduced-cell charge compensation is not neutral")
    changes = [
        {
            "al2_site_id": item[1],
            "oxygen_site_ids": list(item[2]),
            "distance_from_route_angstrom": item[0],
        }
        for item in selected
    ]
    return compensated, changes


def _route_base(
    pristine: Sequence[contract.nteme_neb.Atom],
    cell: contract.nteme_neb.Cell,
    route: contract.nteme_neb.Barrier,
    species: str,
) -> tuple[
    list[contract.nteme_neb.Atom],
    tuple[contract.nteme_neb.Atom, ...],
    dict[str, Any],
]:
    by_id = {atom.id: atom for atom in pristine}
    if {route.moving_site, route.vacancy_1, route.vacancy_2} - set(by_id):
        raise ValueError("selected 3x2 cell does not contain the pinned E3a route")
    moving = by_id[route.moving_site]
    destination = by_id[route.vacancy_1]
    removed_ids = {route.vacancy_1, route.vacancy_2}
    atoms = [
        replace(atom, element=species) if atom.id == route.moving_site else atom
        for atom in pristine
        if atom.id not in removed_ids
    ]
    center = contract.periodic_midpoint(moving, destination, cell)
    atoms, compensation = compensate_reduced_route(
        atoms, cell, center, len(removed_ids) + 1
    )
    dx, dy, dz = contract.nteme_neb.minimum_image_delta(moving, destination, cell)
    endpoint_atom = contract.nteme_neb.Atom(
        route.moving_site,
        species,
        moving.x + dx,
        moving.y + dy,
        moving.z + dz,
    )
    endpoint = tuple(
        replace(atom, x=endpoint_atom.x, y=endpoint_atom.y, z=endpoint_atom.z)
        if atom.id == route.moving_site
        else atom
        for atom in atoms
    )
    metadata = {
        "outcome": "prepared-periodic-dft-spot-check",
        "route": {
            "moving_site_id": route.moving_site,
            "vacancy_site_ids": sorted(removed_ids),
            "final_coordinate_angstrom": [
                endpoint_atom.x,
                endpoint_atom.y,
                endpoint_atom.z,
            ],
            "remote_charge_compensation": compensation,
            "charge_compensation_provenance": (
                "reconstructed-charge-compensation in the neutral 3x2 route cell; "
                "not Nteme's undisclosed production compensator placement"
            ),
        },
    }
    return atoms, endpoint, metadata


def build_spot_models(
    pristine: Sequence[contract.nteme_neb.Atom],
    cell: contract.nteme_neb.Cell,
    route: contract.nteme_neb.Barrier,
) -> list[contract.Model]:
    ar_atoms, ar_endpoint, ar_metadata = _route_base(pristine, cell, route, "Ar")
    replication = contract.Model(
        name="reconstructed-replication",
        atoms=tuple(ar_atoms),
        cell=cell,
        metadata={
            **ar_metadata,
            "published_barrier_kcal_mol": route.barrier_kcal_mol,
            "limitation": (
                "This is a neutral reconstructed 3x2 calibration cell, not the "
                "unpublished Nteme production structure or the E3a 6x3 cell."
            ),
        },
        endpoint_atoms=ar_endpoint,
    )

    by_id = {atom.id: atom for atom in pristine}
    route_center = contract.periodic_midpoint(
        by_id[route.moving_site], by_id[route.vacancy_1], cell
    )
    first_o, second_o = contract.choose_oh_pair(ar_atoms, cell, route_center)
    keep_o, remove_o = sorted((first_o, second_o), key=lambda atom: atom.id)
    hydrogen_by_oxygen = dict(contract.remaining_oh_bonds(pristine, cell))
    removed = {
        remove_o.id,
        hydrogen_by_oxygen[keep_o.id],
        hydrogen_by_oxygen[remove_o.id],
    }
    dehyd_atoms = tuple(
        replace(atom, element="O2") if atom.id == keep_o.id else atom
        for atom in ar_atoms
        if atom.id not in removed
    )
    dehyd_endpoint = tuple(
        replace(atom, element="O2") if atom.id == keep_o.id else atom
        for atom in ar_endpoint
        if atom.id not in removed
    )
    dehyd_metadata = {
        **ar_metadata,
        "transformation": {
            "reaction": "2 OH(lattice) -> H2O(removed) + O_residual(lattice)",
            "selected_oh_oxygen_ids": [keep_o.id, remove_o.id],
            "selected_hydrogen_ids": [
                hydrogen_by_oxygen[keep_o.id],
                hydrogen_by_oxygen[remove_o.id],
            ],
            "residual_oxygen_site_id": keep_o.id,
            "removed_site_ids": sorted(removed),
            "residual_oxygen_type": "O2",
        },
        "limitation": (
            "The DFT path is only noble-gas migration in the local post-reaction "
            "topology; it is not a dehydroxylation path."
        ),
    }
    dehyd = contract.Model(
        name="dehydroxylate-lattice",
        atoms=dehyd_atoms,
        cell=cell,
        metadata=dehyd_metadata,
        endpoint_atoms=dehyd_endpoint,
    )

    xe_atoms, xe_endpoint, xe_metadata = _route_base(pristine, cell, route, "Xe")
    xenon = contract.Model(
        name="xenon-divacancy",
        atoms=tuple(xe_atoms),
        cell=cell,
        metadata={
            **xe_metadata,
            "limitation": (
                "The DFT result can test the E3a Xe screen, but one 3x2 route "
                "does not validate all dry-muscovite Xe environments."
            ),
        },
        endpoint_atoms=xe_endpoint,
    )
    return [replication, dehyd, xenon]


def atom_map(model: contract.Model) -> dict[str, Any]:
    atoms = [
        {
            "cp2k_index": index,
            "source_id": atom.id,
            "source_type": atom.element,
            "chemical_element": CHEMICAL_ELEMENT[atom.element],
        }
        for index, atom in enumerate(sorted(model.atoms, key=lambda item: item.id), 1)
    ]
    return {"schema": "e3b-atom-map-v1", "model": model.name, "atoms": atoms}


def _coordination_gate(model: contract.Model) -> dict[str, Any]:
    if model.name != "dehydroxylate-lattice":
        return {"coordination_pass": True, "five_coordinate_al_ids": []}
    oxygen = [atom for atom in model.atoms if atom.element.startswith("O")]
    al_coordination = {
        atom.id: sum(
            contract.periodic_distance(atom, other, model.cell) <= AL_O_CUTOFF_ANGSTROM
            for other in oxygen
        )
        for atom in model.atoms
        if atom.element == "Al1"
    }
    five_coordinate = sorted(
        atom_id for atom_id, count in al_coordination.items() if count == 5
    )
    unexpected = {
        atom_id: count
        for atom_id, count in al_coordination.items()
        if count not in {5, 6}
    }
    route = model.metadata["route"]
    assert isinstance(route, dict)
    route_center = contract.periodic_midpoint(
        next(atom for atom in model.atoms if atom.id == int(route["moving_site_id"])),
        contract.nteme_neb.Atom(
            0,
            "X",
            *[float(value) for value in route["final_coordinate_angstrom"]],
        ),
        model.cell,
    )
    by_id = {atom.id: atom for atom in model.atoms}
    route_distances = {
        str(atom_id): contract.periodic_distance(
            route_center, by_id[atom_id], model.cell
        )
        for atom_id in five_coordinate
    }
    passed = (
        five_coordinate == [10, 11]
        and not unexpected
        and all(
            value <= LOCAL_ROUTE_CUTOFF_ANGSTROM for value in route_distances.values()
        )
    )
    return {
        "coordination_pass": passed,
        "five_coordinate_al_ids": five_coordinate,
        "al_o_cutoff_angstrom": AL_O_CUTOFF_ANGSTROM,
        "distance_to_route_midpoint_angstrom": route_distances,
        "unexpected_al_coordination": unexpected,
    }


def validate_spot_model(
    model: contract.Model, source_cell: contract.nteme_neb.Cell
) -> dict[str, Any]:
    if model.endpoint_atoms is None:
        raise ValueError(f"{model.name}: missing endpoint")
    validation = contract.validate_model(model)
    endpoint_model = replace(model, atoms=model.endpoint_atoms, endpoint_atoms=None)
    endpoint_validation = contract.validate_model(endpoint_model)
    counts = Counter(atom.element for atom in model.atoms)
    composition_pass = counts == EXPECTED_MODEL_COUNTS[model.name]
    net_charge = contract.net_charge(model.atoms)
    endpoint_charge = contract.net_charge(model.endpoint_atoms)
    neutrality_pass = math.isclose(net_charge, 0.0, abs_tol=1.0e-8) and math.isclose(
        endpoint_charge, 0.0, abs_tol=1.0e-8
    )
    cell_pass = all(
        (
            math.isclose(model.cell.a, source_cell.a / 2.0, abs_tol=1.0e-10),
            math.isclose(model.cell.b, source_cell.b * 2.0 / 3.0, abs_tol=1.0e-10),
            math.isclose(model.cell.c, source_cell.c, abs_tol=1.0e-10),
            math.isclose(model.cell.alpha, source_cell.alpha, abs_tol=1.0e-10),
            math.isclose(model.cell.beta, source_cell.beta, abs_tol=1.0e-10),
            math.isclose(model.cell.gamma, source_cell.gamma, abs_tol=1.0e-10),
        )
    )
    initial = sorted(model.atoms, key=lambda atom: atom.id)
    endpoint = sorted(model.endpoint_atoms, key=lambda atom: atom.id)
    identity_pass = [(atom.id, atom.element) for atom in initial] == [
        (atom.id, atom.element) for atom in endpoint
    ]
    route = model.metadata.get("route")
    if not isinstance(route, dict):
        raise TypeError(f"{model.name}: missing typed route metadata")
    moving_id = int(route["moving_site_id"])
    vacancy_ids = {int(value) for value in route["vacancy_site_ids"]}
    initial_by_id = {atom.id: atom for atom in initial}
    endpoint_by_id = {atom.id: atom for atom in endpoint}
    fixed_identity = all(
        atom_id == moving_id
        or all(
            math.isclose(
                getattr(initial_by_id[atom_id], axis),
                getattr(endpoint_by_id[atom_id], axis),
                abs_tol=1.0e-10,
            )
            for axis in ("x", "y", "z")
        )
        for atom_id in initial_by_id
    )
    final_coordinate = tuple(
        float(value) for value in route["final_coordinate_angstrom"]
    )
    moving_endpoint = endpoint_by_id[moving_id]
    route_pass = (
        moving_id in initial_by_id
        and not vacancy_ids.intersection(initial_by_id)
        and fixed_identity
        and all(
            math.isclose(getattr(moving_endpoint, axis), value, abs_tol=1.0e-10)
            for axis, value in zip(("x", "y", "z"), final_coordinate, strict=True)
        )
    )
    coordination = _coordination_gate(model)
    chemical_counts = Counter(CHEMICAL_ELEMENT[atom.element] for atom in model.atoms)
    valence_electrons = sum(
        count * VALENCE_ELECTRONS[element] for element, count in chemical_counts.items()
    )
    initial_minimum = validation["minimum_pair_distance_angstrom"]
    endpoint_minimum = endpoint_validation["minimum_pair_distance_angstrom"]
    minimum_distance_pass = (
        isinstance(initial_minimum, (int, float))
        and isinstance(endpoint_minimum, (int, float))
        and float(initial_minimum) >= MINIMUM_DISTANCE_ANGSTROM
        and float(endpoint_minimum) >= MINIMUM_DISTANCE_ANGSTROM
    )
    gates = {
        "composition_pass": composition_pass,
        "neutrality_pass": neutrality_pass,
        "net_charge_e3a_e": net_charge,
        "endpoint_net_charge_e3a_e": endpoint_charge,
        "cp2k_total_charge_e": 0,
        "cp2k_valence_electrons": valence_electrons,
        "singlet_electron_parity_pass": valence_electrons % 2 == 0,
        "cell_pass": cell_pass,
        "route_pass": route_pass,
        "identity_pass": identity_pass,
        "minimum_distance_pass": minimum_distance_pass,
        **coordination,
    }
    gates["all_pass"] = all(
        value is True for key, value in gates.items() if key.endswith("_pass")
    )
    if not gates["all_pass"]:
        raise ValueError(f"{model.name}: preparation gate failed: {gates}")
    return gates


def validate_method(method: MethodProvenance) -> None:
    if not re.fullmatch(r"[0-9][0-9A-Za-z.+-]*", method.cp2k_version):
        raise ValueError("CP2K version must be explicit")
    if not re.fullmatch(r"sha256:[0-9a-f]{64}", method.cp2k_image_digest):
        raise ValueError("CP2K image digest must be an immutable sha256 digest")
    for label, path in (
        ("basis set file", method.basis_set_file),
        ("potential file", method.potential_file),
        ("D3 parameter file", method.d3_parameter_file),
    ):
        if not path.is_file():
            raise ValueError(f"{label} is not a file: {path}")


def materialize_method_files(
    out: pathlib.Path, method: MethodProvenance
) -> dict[str, Any]:
    validate_method(method)
    target = out / "cp2k-data"
    target.mkdir()
    files = (
        (method.basis_set_file, target / "BASIS_SET"),
        (method.potential_file, target / "GTH_POTENTIALS"),
        (method.d3_parameter_file, target / "dftd3.dat"),
    )
    records = {}
    for source, destination in files:
        shutil.copyfile(source, destination)
        records[destination.name] = {
            "source_path": str(source.resolve()),
            "artifact_path": str(destination.relative_to(out)),
            "sha256": sha256(destination),
            "bytes": destination.stat().st_size,
        }
    return {
        "electronic_structure": "periodic GPW PBE-D3 (zero damping)",
        "cp2k_version": method.cp2k_version,
        "cp2k_image_digest": method.cp2k_image_digest,
        "xc_functional": "PBE",
        "dispersion": {
            "type": "DFTD3",
            "damping": "zero",
            "reference_functional": "PBE",
            "parameter_file": files_record(records, "dftd3.dat"),
            "r_cutoff_angstrom": 15.0,
        },
        "basis": BASIS_SET,
        "potential_family": "GTH-PBE with explicit qN labels",
        "plane_wave_cutoff_ry": 600,
        "relative_cutoff_ry": 60,
        "k_points": "GAMMA",
        "eps_scf_hartree": 1.0e-6,
        "neb_endpoint_policy": (
            "endpoint replicas are optimized in-band; separate fixed-cell endpoint "
            "GEO_OPT runs must also pass the analysis convergence gate"
        ),
        "files": records,
    }


def files_record(records: Mapping[str, Any], name: str) -> dict[str, Any]:
    value = records[name]
    assert isinstance(value, dict)
    return dict(value)


def _kind_sections(model: contract.Model) -> list[str]:
    sections = []
    for source_type in sorted({atom.element for atom in model.atoms}):
        element = CHEMICAL_ELEMENT[source_type]
        sections.extend(
            (
                f"    &KIND {source_type}",
                f"      ELEMENT {element}",
                f"      BASIS_SET {BASIS_SET}",
                f"      POTENTIAL {POTENTIAL[element]}",
                "    &END KIND",
            )
        )
    return sections


def _force_eval_lines(model: contract.Model, coord_file_name: str) -> list[str]:
    lx, ly, lz, xy, xz, yz = model.cell.restricted()
    return [
        "&FORCE_EVAL",
        "  METHOD Quickstep",
        "  &DFT",
        "    BASIS_SET_FILE_NAME ../../cp2k-data/BASIS_SET",
        "    POTENTIAL_FILE_NAME ../../cp2k-data/GTH_POTENTIALS",
        "    CHARGE 0",
        "    MULTIPLICITY 1",
        "    &MGRID",
        "      CUTOFF 600",
        "      REL_CUTOFF 60",
        "    &END MGRID",
        "    &SCF",
        "      EPS_SCF 1.0E-6",
        "      MAX_SCF 100",
        "      SCF_GUESS ATOMIC",
        "      &OT",
        "        MINIMIZER DIIS",
        "        PRECONDITIONER FULL_SINGLE_INVERSE",
        "      &END OT",
        "      &OUTER_SCF",
        "        EPS_SCF 1.0E-6",
        "        MAX_SCF 10",
        "      &END OUTER_SCF",
        "    &END SCF",
        "    &XC",
        "      &XC_FUNCTIONAL PBE",
        "      &END XC_FUNCTIONAL",
        "      &VDW_POTENTIAL",
        "        POTENTIAL_TYPE PAIR_POTENTIAL",
        "        &PAIR_POTENTIAL",
        "          TYPE DFTD3",
        "          PARAMETER_FILE_NAME ../../cp2k-data/dftd3.dat",
        "          REFERENCE_FUNCTIONAL PBE",
        "          R_CUTOFF 15.0",
        "        &END PAIR_POTENTIAL",
        "      &END VDW_POTENTIAL",
        "    &END XC",
        "  &END DFT",
        "  &SUBSYS",
        "    &CELL",
        f"      A {lx:.12f} 0.000000000000 0.000000000000",
        f"      B {xy:.12f} {ly:.12f} 0.000000000000",
        f"      C {xz:.12f} {yz:.12f} {lz:.12f}",
        "      PERIODIC XYZ",
        "    &END CELL",
        "    &TOPOLOGY",
        f"      COORD_FILE_NAME {coord_file_name}",
        "      COORD_FILE_FORMAT XYZ",
        "    &END TOPOLOGY",
        *_kind_sections(model),
        "  &END SUBSYS",
        "&END FORCE_EVAL",
    ]


def write_cp2k_endpoint_input(
    path: pathlib.Path, model: contract.Model, label: str, coord_file_name: str
) -> None:
    lines = [
        "&GLOBAL",
        f"  PROJECT e3b-{model.name}-{label}",
        "  RUN_TYPE GEO_OPT",
        "  PRINT_LEVEL MEDIUM",
        "&END GLOBAL",
        *_force_eval_lines(model, coord_file_name),
        "&MOTION",
        "  &GEO_OPT",
        "    OPTIMIZER BFGS",
        "    MAX_ITER 200",
        "    MAX_FORCE 4.5E-4",
        "    RMS_FORCE 3.0E-4",
        "  &END GEO_OPT",
        "  &PRINT",
        "    &TRAJECTORY",
        "      FORMAT XYZ",
        "    &END TRAJECTORY",
        "  &END PRINT",
        "&END MOTION",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_cp2k_smoke_input(
    path: pathlib.Path, model: contract.Model, coord_file_name: str
) -> None:
    lines = [
        "&GLOBAL",
        f"  PROJECT e3b-{model.name}-smoke",
        "  RUN_TYPE ENERGY_FORCE",
        "  PRINT_LEVEL MEDIUM",
        "&END GLOBAL",
        *_force_eval_lines(model, coord_file_name),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_cp2k_neb_input(path: pathlib.Path, model: contract.Model) -> None:
    lines = [
        "&GLOBAL",
        f"  PROJECT e3b-{model.name}-neb",
        "  RUN_TYPE BAND",
        "  PRINT_LEVEL MEDIUM",
        "&END GLOBAL",
        *_force_eval_lines(model, "images/replica-00.xyz"),
        "&MOTION",
        "  &BAND",
        "    BAND_TYPE CI-NEB",
        "    NUMBER_OF_REPLICA 8",
        "    K_SPRING 0.05",
        "    ALIGN_FRAMES F",
        "    ROTATE_FRAMES F",
        "    &CONVERGENCE_CONTROL",
        "      MAX_FORCE 4.5E-4",
        "      RMS_FORCE 3.0E-4",
        "    &END CONVERGENCE_CONTROL",
        "    &OPTIMIZE_BAND",
        "      OPT_TYPE DIIS",
        "      OPTIMIZE_END_POINTS T",
        "    &END OPTIMIZE_BAND",
    ]
    for index in range(IMAGE_COUNT):
        lines.extend(
            (
                "    &REPLICA",
                f"      COORD_FILE_NAME images/replica-{index:02d}.xyz",
                "    &END REPLICA",
            )
        )
    lines.extend(("  &END BAND", "&END MOTION"))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_xyz(
    path: pathlib.Path,
    atoms: Sequence[contract.nteme_neb.Atom],
    *,
    comment: str,
) -> None:
    ordered = sorted(atoms, key=lambda atom: atom.id)
    lines = [str(len(ordered)), comment]
    lines.extend(
        f"{atom.element:<4} {atom.x:.12f} {atom.y:.12f} {atom.z:.12f}"
        for atom in ordered
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_images(model_dir: pathlib.Path, model: contract.Model) -> list[pathlib.Path]:
    if model.endpoint_atoms is None:
        raise ValueError(f"{model.name}: missing endpoint")
    initial = sorted(model.atoms, key=lambda atom: atom.id)
    endpoint = sorted(model.endpoint_atoms, key=lambda atom: atom.id)
    if [(atom.id, atom.element) for atom in initial] != [
        (atom.id, atom.element) for atom in endpoint
    ]:
        raise ValueError(f"{model.name}: endpoint atom identity drift")
    image_dir = model_dir / "images"
    image_dir.mkdir()
    paths = []
    for index in range(IMAGE_COUNT):
        fraction = index / (IMAGE_COUNT - 1)
        image = [
            replace(
                first,
                x=first.x + fraction * (last.x - first.x),
                y=first.y + fraction * (last.y - first.y),
                z=first.z + fraction * (last.z - first.z),
            )
            for first, last in zip(initial, endpoint, strict=True)
        ]
        path = image_dir / f"replica-{index:02d}.xyz"
        write_xyz(
            path,
            image,
            comment=(
                f"E3b {model.name} CI-NEB image {index + 1}/{IMAGE_COUNT}; "
                "identity in ../atom-map.json"
            ),
        )
        paths.append(path)
    return paths


def write_matched_classical_inputs(
    model_dir: pathlib.Path, model: contract.Model
) -> None:
    if model.endpoint_atoms is None:
        raise ValueError(f"{model.name}: missing endpoint")
    target = model_dir / "matched-classical"
    target.mkdir()
    contract.write_lammps_data(target / "initial.seed.data", model)
    endpoint_model = replace(model, atoms=model.endpoint_atoms, endpoint_atoms=None)
    contract.write_lammps_data(target / "endpoint.seed.data", endpoint_model)
    for label in ("initial", "endpoint"):
        e3a_campaign.write_minimize_input(
            target / f"in.min.{label}",
            data_name=f"{label}.seed.data",
            output_prefix=f"{label}.relaxed",
            ftol=0.01,
            max_steps=5000,
        )
    e3a_campaign.write_neb_input(
        target / "in.neb",
        data_name="initial.relaxed.data",
        final_name="final.coords",
        ftol=0.01,
        relax_steps=3000,
        climb_steps=10000,
    )


def _classical_reference(campaign: Mapping[str, Any], name: str) -> dict[str, Any]:
    models = campaign["models"]
    if not isinstance(models, dict) or not isinstance(models[name], dict):
        raise TypeError(f"E3a campaign has malformed model {name}")
    neb = models[name].get("neb")
    if not isinstance(neb, dict):
        raise TypeError(f"E3a campaign model {name} has no NEB result")
    barrier = neb.get("computed_forward_barrier_kcal_mol")
    if not isinstance(barrier, (int, float)) or not math.isfinite(float(barrier)):
        raise ValueError(f"E3a campaign model {name} has no finite barrier")
    converged = neb.get("converged_to_requested_ftol") is True
    return {
        "barrier_kcal_mol": float(barrier),
        "status": "converged" if converged else "incomplete-convergence",
        "source_cell": [6, 3, 1],
        "result_sha256": CAMPAIGN_RESULT_SHA256,
    }


def write_manifest(out_dir: pathlib.Path) -> list[dict[str, Any]]:
    manifest_path = out_dir / "manifest.json"
    files = []
    for path in sorted(
        item
        for item in out_dir.rglob("*")
        if item.is_file() and item.resolve() != manifest_path.resolve()
    ):
        files.append(
            {
                "path": str(path.relative_to(out_dir)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    write_json(manifest_path, {"schema": "sha256-manifest-v1", "files": files})
    return files


def verify_manifest(out_dir: pathlib.Path) -> None:
    manifest = read_json(out_dir / "manifest.json")
    files = manifest.get("files")
    if not isinstance(files, list):
        raise TypeError("prepared manifest has no files list")
    expected_paths = set()
    for entry in files:
        if not isinstance(entry, dict):
            raise TypeError("prepared manifest has a malformed entry")
        relative = entry.get("path")
        expected = entry.get("sha256")
        if not isinstance(relative, str) or not isinstance(expected, str):
            raise TypeError("prepared manifest entry is not typed")
        if pathlib.Path(relative).is_absolute() or ".." in pathlib.Path(relative).parts:
            raise ValueError("prepared manifest contains an unsafe path")
        path = out_dir / relative
        require_sha256(path, expected, label=f"prepared artifact {relative}")
        if path.stat().st_size != entry.get("bytes"):
            raise ValueError(f"prepared artifact size mismatch: {relative}")
        expected_paths.add(relative)
    manifest_path = (out_dir / "manifest.json").resolve()
    actual_paths = {
        str(path.relative_to(out_dir))
        for path in out_dir.rglob("*")
        if path.is_file() and path.resolve() != manifest_path
    }
    if actual_paths != expected_paths:
        raise ValueError("prepared artifact path set does not match its manifest")


def prepare_campaign(
    *,
    workbook: pathlib.Path,
    campaign_root: pathlib.Path,
    out: pathlib.Path,
    method: MethodProvenance,
) -> dict[str, Any]:
    campaign = verify_e3a_inputs(workbook, campaign_root)
    prepare_output_dir(out)
    method_record = materialize_method_files(out, method)
    source_atoms, source_cell, barriers = contract.nteme_neb.read_workbook(workbook)
    route = contract.nteme_neb.select_route(barriers, "divacancy", 1)
    pristine, target_cell, reduction = reduce_pristine(source_atoms, source_cell)
    models = build_spot_models(pristine, target_cell, route)
    model_records: dict[str, Any] = {}
    models_root = out / "models"
    models_root.mkdir()
    for model in models:
        gates = validate_spot_model(model, source_cell)
        mapping = atom_map(model)
        identity_sha = canonical_sha256(mapping["atoms"])
        model_dir = models_root / model.name
        model_dir.mkdir()
        write_json(model_dir / "atom-map.json", mapping)
        assert model.endpoint_atoms is not None
        write_xyz(
            model_dir / "initial.xyz",
            model.atoms,
            comment=f"E3b {model.name} initial; identity in atom-map.json",
        )
        write_xyz(
            model_dir / "endpoint.xyz",
            model.endpoint_atoms,
            comment=f"E3b {model.name} endpoint; identity in atom-map.json",
        )
        write_images(model_dir, model)
        write_cp2k_endpoint_input(
            model_dir / "initial-opt.inp", model, "initial", "initial.xyz"
        )
        write_cp2k_endpoint_input(
            model_dir / "endpoint-opt.inp", model, "endpoint", "endpoint.xyz"
        )
        write_cp2k_smoke_input(model_dir / "smoke.inp", model, "initial.xyz")
        write_cp2k_neb_input(model_dir / "neb.inp", model)
        write_matched_classical_inputs(model_dir, model)
        source_counts = Counter(atom.element for atom in model.atoms)
        chemical_counts = Counter(
            CHEMICAL_ELEMENT[atom.element] for atom in model.atoms
        )
        artifacts = {
            name: sha256(model_dir / name)
            for name in (
                "atom-map.json",
                "initial.xyz",
                "endpoint.xyz",
                "initial-opt.inp",
                "endpoint-opt.inp",
                "smoke.inp",
                "neb.inp",
            )
        }
        model_records[model.name] = {
            "scientific_scope": SCIENTIFIC_SCOPE[model.name],
            "atom_count": len(model.atoms),
            "source_type_counts": dict(sorted(source_counts.items())),
            "chemical_element_counts": dict(sorted(chemical_counts.items())),
            "supercell": list(TARGET_SUPERCELL),
            "periodic_cell": asdict(model.cell),
            "route": model.metadata["route"],
            "transformation": model.metadata.get("transformation"),
            "limitation": model.metadata["limitation"],
            "gates": gates,
            "atom_identity_sha256": identity_sha,
            "classical_reference": _classical_reference(campaign, model.name),
            "matched_classical": {
                "status": "prepared-not-run",
                "cell": list(TARGET_SUPERCELL),
                "images": IMAGE_COUNT,
                "required_for_transfer_gate": True,
            },
            "artifacts": artifacts,
        }
    preparation: dict[str, Any] = {
        "schema": "e3b-cp2k-preparation-v1",
        "source": {
            "workbook_path": str(workbook.resolve()),
            "workbook_sha256": WORKBOOK_SHA256,
            "campaign_root": str(campaign_root.resolve()),
            "campaign_manifest_sha256": CAMPAIGN_MANIFEST_SHA256,
            "campaign_result_sha256": CAMPAIGN_RESULT_SHA256,
            "source_supercell": list(SOURCE_SUPERCELL),
            "reduction": asdict(reduction),
        },
        "method": method_record,
        "execution_envelope": {
            "classification": "platform-test survey-tier spot checks",
            "per_qm_unit_wall_seconds_max": 14400,
            "endpoint_wall_seconds_each": 1800,
            "neb_wall_seconds_each": 14400,
            "smoke_wall_seconds_each": 600,
            "planned_total_wall_seconds_max": 55800,
            "images": IMAGE_COUNT,
            "heavy_compute_run_by_prepare": False,
        },
        "models": model_records,
    }
    preparation_path = out / "preparation.json"
    write_json(preparation_path, preparation)
    observation_template = {
        "schema": "e3b-cp2k-observations-v1",
        "preparation_sha256": sha256(preparation_path),
        "models": {
            name: {
                "atom_identity_sha256": record["atom_identity_sha256"],
                "cp2k": {
                    "status": "incomplete-execution",
                    "timed_out": None,
                    "scf_converged_all_images": None,
                    "endpoints_converged": None,
                    "neb_converged": None,
                    "barrier_kcal_mol": None,
                },
                "matched_classical": {
                    "status": "incomplete-execution",
                    "timed_out": None,
                    "neb_converged": None,
                    "barrier_kcal_mol": None,
                },
            }
            for name, record in model_records.items()
        },
    }
    write_json(out / "observations.template.json", observation_template)
    write_manifest(out)
    return preparation


def classify_cp2k_output(
    output: str, *, returncode: int | None, timed_out: bool
) -> dict[str, Any]:
    lower = output.lower()
    scf_failed = "scf run not converged" in lower or "scf run did not converge" in lower
    scf_converged = bool(re.search(r"scf run converged in\s+\d+\s+steps", lower))
    normal_end = "program ended at" in lower
    energies = re.findall(
        r"ENERGY\|\s+Total FORCE_EVAL.*?energy\s+\[a\.u\.\]:\s+"
        r"([-+0-9.eEdD]+)",
        output,
    )
    energy = (
        float(energies[-1].replace("D", "E").replace("d", "e")) if energies else None
    )
    if timed_out:
        status = "incomplete-timeout"
    elif scf_failed or not scf_converged:
        status = "incomplete-scf"
    elif returncode != 0 or not normal_end or energy is None:
        status = "incomplete-execution"
    else:
        status = "converged"
    return {
        "status": status,
        "timed_out": timed_out,
        "returncode": returncode,
        "scf_converged": scf_converged and not scf_failed,
        "normal_termination": normal_end,
        "energy_hartree": energy,
    }


def run_smoke(
    *,
    executable: str,
    input_path: pathlib.Path,
    receipt_path: pathlib.Path,
    timeout_seconds: float,
) -> dict[str, Any]:
    if timeout_seconds <= 0 or timeout_seconds > 14400:
        raise ValueError("smoke timeout must be in (0, 14400] seconds")
    if not input_path.is_file():
        raise ValueError(f"smoke input is not a file: {input_path}")
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    stdout_path = receipt_path.parent / "smoke.stdout.log"
    started = time.monotonic()
    timed_out = False
    returncode: int | None = None
    output = ""
    try:
        environment = os.environ.copy()
        environment.update(
            {
                "OMP_NUM_THREADS": "1",
                "OPENBLAS_NUM_THREADS": "1",
                "MKL_NUM_THREADS": "1",
            }
        )
        # Trusted local executable and fixed CP2K -i argument; never a shell string.
        # nosemgrep: python.lang.security.audit.dangerous-subprocess-use-audit
        completed = subprocess.run(
            [executable, "-i", str(input_path.resolve())],
            cwd=input_path.parent,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
            timeout=timeout_seconds,
            shell=False,
        )
        returncode = completed.returncode
        output = completed.stdout
    except subprocess.TimeoutExpired as exc:
        timed_out = True
        if isinstance(exc.stdout, bytes):
            output = exc.stdout.decode(errors="replace")
        else:
            output = exc.stdout or ""
    elapsed = time.monotonic() - started
    stdout_path.write_text(output, encoding="utf-8")
    receipt = {
        **classify_cp2k_output(
            output,
            returncode=returncode,
            timed_out=timed_out,
        ),
        "schema": "e3b-cp2k-smoke-v1",
        "command": [executable, "-i", str(input_path.resolve())],
        "timeout_seconds": timeout_seconds,
        "elapsed_seconds": elapsed,
        "input_sha256": sha256(input_path),
        "stdout_path": str(stdout_path.resolve()),
        "stdout_sha256": sha256(stdout_path),
    }
    write_json(receipt_path, receipt)
    return receipt


def _finite_number(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def _failure_result(outcome: str, reason: str) -> dict[str, Any]:
    return {
        "typed_outcome": outcome,
        "reason": reason,
        "matched_cell_transfer_pass": False,
        "matched_cell_transfer_delta_kcal_mol": None,
        "matched_cell_correction_kcal_mol": None,
        "calibrated_barrier_kcal_mol": None,
        "verdict": "no-correction",
    }


def analyze_models(
    preparation: Mapping[str, Any],
    observations: Mapping[str, Any],
    *,
    transfer_tolerance_kcal_mol: float = DEFAULT_TRANSFER_TOLERANCE_KCAL_MOL,
    endorsement_tolerance_kcal_mol: float = DEFAULT_ENDORSEMENT_TOLERANCE_KCAL_MOL,
) -> dict[str, Any]:
    if preparation.get("schema") != "e3b-cp2k-preparation-v1":
        raise ValueError("unexpected E3b preparation schema")
    if observations.get("schema") != "e3b-cp2k-observations-v1":
        raise ValueError("unexpected E3b observations schema")
    if (
        not math.isfinite(transfer_tolerance_kcal_mol)
        or not math.isfinite(endorsement_tolerance_kcal_mol)
        or transfer_tolerance_kcal_mol < 0
        or endorsement_tolerance_kcal_mol < 0
    ):
        raise ValueError("analysis tolerances must be non-negative")
    prepared_models = preparation.get("models")
    observed_models = observations.get("models")
    if not isinstance(prepared_models, Mapping) or not isinstance(
        observed_models, Mapping
    ):
        raise TypeError("preparation and observations require typed model maps")
    if set(prepared_models) != set(observed_models):
        raise ValueError("observed model set does not exactly match preparation")
    results: dict[str, Any] = {}
    for name, prepared_value in prepared_models.items():
        observed_value = observed_models[name]
        if not isinstance(prepared_value, Mapping) or not isinstance(
            observed_value, Mapping
        ):
            raise TypeError(f"{name}: malformed model record")
        if observed_value.get("atom_identity_sha256") != prepared_value.get(
            "atom_identity_sha256"
        ):
            results[name] = _failure_result(
                "incomplete-atom-identity",
                "observed atom identity does not match the prepared path",
            )
            continue
        cp2k = observed_value.get("cp2k")
        matched = observed_value.get("matched_classical")
        reference = prepared_value.get("classical_reference")
        if not isinstance(cp2k, Mapping) or not isinstance(matched, Mapping):
            raise TypeError(f"{name}: missing typed CP2K or classical observation")
        if not isinstance(reference, Mapping):
            raise TypeError(f"{name}: missing classical reference")
        cp2k_status = cp2k.get("status")
        matched_status = matched.get("status")
        if cp2k_status not in ALLOWED_OBSERVATION_STATUS:
            raise ValueError(f"{name}: unknown CP2K status {cp2k_status!r}")
        if matched_status not in ALLOWED_OBSERVATION_STATUS:
            raise ValueError(
                f"{name}: unknown matched-classical status {matched_status!r}"
            )
        if cp2k_status != "converged":
            results[name] = _failure_result(
                str(cp2k_status), "CP2K result is not typed converged"
            )
            continue
        if cp2k.get("timed_out") is not False:
            results[name] = _failure_result(
                "incomplete-timeout", "CP2K timeout state is not explicitly false"
            )
            continue
        if cp2k.get("scf_converged_all_images") is not True:
            results[name] = _failure_result(
                "incomplete-scf", "not every endpoint/image has converged SCF evidence"
            )
            continue
        if (
            cp2k.get("endpoints_converged") is not True
            or cp2k.get("neb_converged") is not True
        ):
            results[name] = _failure_result(
                "incomplete-convergence", "CP2K endpoints or CI-NEB are not converged"
            )
            continue
        cp2k_barrier = _finite_number(cp2k.get("barrier_kcal_mol"))
        if cp2k_barrier is None or cp2k_barrier < 0:
            results[name] = _failure_result(
                "incomplete-convergence", "CP2K barrier is missing or non-finite"
            )
            continue
        if matched_status != "converged":
            results[name] = _failure_result(
                "incomplete-matched-classical",
                "matched 3x2 classical CI-NEB is not typed converged",
            )
            continue
        if (
            matched.get("timed_out") is not False
            or matched.get("neb_converged") is not True
        ):
            results[name] = _failure_result(
                "incomplete-matched-classical",
                "matched classical timeout/convergence evidence is incomplete",
            )
            continue
        matched_barrier = _finite_number(matched.get("barrier_kcal_mol"))
        reference_barrier = _finite_number(reference.get("barrier_kcal_mol"))
        if (
            matched_barrier is None
            or matched_barrier < 0
            or reference_barrier is None
            or reference_barrier < 0
        ):
            results[name] = _failure_result(
                "incomplete-matched-classical",
                "matched or E3a classical barrier is missing/non-finite",
            )
            continue
        transfer_delta = matched_barrier - reference_barrier
        transfer_pass = abs(transfer_delta) <= transfer_tolerance_kcal_mol
        if not transfer_pass:
            result = _failure_result(
                "incomplete-transfer-gate",
                "3x2 matched classical barrier does not transfer to the E3a 6x3 cell",
            )
            result["matched_cell_transfer_delta_kcal_mol"] = transfer_delta
            result["transfer_tolerance_kcal_mol"] = transfer_tolerance_kcal_mol
            results[name] = result
            continue
        correction = cp2k_barrier - matched_barrier
        calibrated = reference_barrier + correction
        endorsement = abs(correction) <= endorsement_tolerance_kcal_mol
        results[name] = {
            "typed_outcome": (
                "complete-endorsement" if endorsement else "complete-correction"
            ),
            "reason": None,
            "source_classical_status": reference.get("status"),
            "source_classical_barrier_kcal_mol": reference_barrier,
            "matched_classical_barrier_kcal_mol": matched_barrier,
            "cp2k_barrier_kcal_mol": cp2k_barrier,
            "matched_cell_transfer_delta_kcal_mol": transfer_delta,
            "transfer_tolerance_kcal_mol": transfer_tolerance_kcal_mol,
            "matched_cell_transfer_pass": True,
            "matched_cell_correction_kcal_mol": correction,
            "calibrated_barrier_kcal_mol": calibrated,
            "endorsement_tolerance_kcal_mol": endorsement_tolerance_kcal_mol,
            "verdict": (
                "endorse-classical-within-tolerance"
                if endorsement
                else "apply-dft-minus-matched-classical-correction"
            ),
        }
    return {
        "schema": "e3b-cp2k-analysis-v1",
        "transfer_tolerance_kcal_mol": transfer_tolerance_kcal_mol,
        "endorsement_tolerance_kcal_mol": endorsement_tolerance_kcal_mol,
        "models": results,
    }


def command_prepare(args: argparse.Namespace) -> int:
    preparation = prepare_campaign(
        workbook=args.workbook,
        campaign_root=args.campaign_root,
        out=args.out,
        method=MethodProvenance(
            cp2k_version=args.cp2k_version,
            cp2k_image_digest=args.cp2k_image_digest,
            basis_set_file=args.basis_set_file,
            potential_file=args.potential_file,
            d3_parameter_file=args.d3_parameter_file,
        ),
    )
    print(json.dumps(preparation, indent=2, sort_keys=True))
    return 0


def command_smoke(args: argparse.Namespace) -> int:
    verify_manifest(args.prepared_root)
    input_path = args.prepared_root / "models" / args.model / "smoke.inp"
    receipt = run_smoke(
        executable=args.cp2k,
        input_path=input_path,
        receipt_path=args.receipt,
        timeout_seconds=args.timeout_seconds,
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0 if receipt["status"] == "converged" else 2


def command_analyze(args: argparse.Namespace) -> int:
    verify_manifest(args.prepared_root)
    preparation_path = args.prepared_root / "preparation.json"
    preparation = read_json(preparation_path)
    observations = read_json(args.observations)
    if observations.get("preparation_sha256") != sha256(preparation_path):
        raise ValueError("observations are not hash-bound to this preparation")
    prepare_output_dir(args.out)
    result = analyze_models(
        preparation,
        observations,
        transfer_tolerance_kcal_mol=args.transfer_tolerance,
        endorsement_tolerance_kcal_mol=args.endorsement_tolerance,
    )
    result["preparation_sha256"] = sha256(preparation_path)
    result["observations_sha256"] = sha256(args.observations)
    write_json(args.out / "analysis.json", result)
    write_manifest(args.out)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(required=True)

    prepare = sub.add_parser("prepare", help="emit the three bounded E3b inputs")
    prepare.add_argument("--workbook", type=pathlib.Path, required=True)
    prepare.add_argument("--campaign-root", type=pathlib.Path, required=True)
    prepare.add_argument("--out", type=pathlib.Path, required=True)
    prepare.add_argument("--cp2k-version", required=True)
    prepare.add_argument("--cp2k-image-digest", required=True)
    prepare.add_argument("--basis-set-file", type=pathlib.Path, required=True)
    prepare.add_argument("--potential-file", type=pathlib.Path, required=True)
    prepare.add_argument("--d3-parameter-file", type=pathlib.Path, required=True)
    prepare.set_defaults(func=command_prepare)

    smoke = sub.add_parser(
        "smoke", help="run one bounded ENERGY_FORCE timing/SCF smoke"
    )
    smoke.add_argument("--prepared-root", type=pathlib.Path, required=True)
    smoke.add_argument("--model", choices=MODEL_ORDER, required=True)
    smoke.add_argument("--cp2k", required=True)
    smoke.add_argument("--receipt", type=pathlib.Path, required=True)
    smoke.add_argument("--timeout-seconds", type=float, default=600.0)
    smoke.set_defaults(func=command_smoke)

    analyze = sub.add_parser(
        "analyze", help="apply only converged, matched-cell DFT corrections"
    )
    analyze.add_argument("--prepared-root", type=pathlib.Path, required=True)
    analyze.add_argument("--observations", type=pathlib.Path, required=True)
    analyze.add_argument("--out", type=pathlib.Path, required=True)
    analyze.add_argument(
        "--transfer-tolerance",
        type=float,
        default=DEFAULT_TRANSFER_TOLERANCE_KCAL_MOL,
    )
    analyze.add_argument(
        "--endorsement-tolerance",
        type=float,
        default=DEFAULT_ENDORSEMENT_TOLERANCE_KCAL_MOL,
    )
    analyze.set_defaults(func=command_analyze)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
