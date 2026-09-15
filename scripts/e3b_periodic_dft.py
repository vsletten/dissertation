#!/usr/bin/env python3
"""Prepare and analyze bounded E3b periodic CP2K spot checks.

This harness performs no DFT during ``prepare`` or ``analyze``. It derives three
2x2x1 models from the hash-pinned E3a source, emits PBE-D3 CP2K endpoint and
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
import signal
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
TARGET_SUPERCELL = (2, 2, 1)
SOURCE_ATOMS_PER_UNIT_CELL = 84
TARGET_PRISTINE_COUNTS = Counter(
    {
        "K": 16,
        "Al1": 32,
        "Al2": 16,
        "Si": 48,
        "O1": 32,
        "O2": 96,
        "O3": 64,
        "H": 32,
    }
)
EXPECTED_MODEL_COUNTS = {
    "reconstructed-replication": Counter(
        {
            "K": 13,
            "Al1": 32,
            "Al2": 13,
            "Si": 51,
            "O1": 32,
            "O2": 108,
            "O3": 52,
            "H": 32,
            "Ar": 1,
        }
    ),
    "dehydroxylate-lattice": Counter(
        {
            "K": 13,
            "Al1": 32,
            "Al2": 13,
            "Si": 51,
            "O1": 30,
            "O2": 109,
            "O3": 52,
            "H": 30,
            "Ar": 1,
        }
    ),
    "xenon-divacancy": Counter(
        {
            "K": 13,
            "Al1": 32,
            "Al2": 13,
            "Si": 51,
            "O1": 32,
            "O2": 108,
            "O3": 52,
            "H": 32,
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
HARTREE_TO_KCAL_MOL = 627.5094740631
MINIMUM_DISTANCE_ANGSTROM = 0.70
AL_O_CUTOFF_ANGSTROM = 2.30
LOCAL_ROUTE_CUTOFF_ANGSTROM = 6.0
MATCHED_CLASSICAL_FTOL_KCAL_MOL_ANGSTROM = 0.01
LAMMPS_ENDPOINT_ENERGY_TOLERANCE_KCAL_MOL = 1.0e-6
CP2K_REPLICA_OUTPUT_ROLES = tuple(
    f"cp2k_replica_{index:02d}" for index in range(1, IMAGE_COUNT + 1)
)
RAW_OUTPUT_ROLES = (
    "cp2k_initial",
    "cp2k_endpoint",
    "cp2k_band",
    "cp2k_band_profile",
    *CP2K_REPLICA_OUTPUT_ROLES,
    "lammps_initial",
    "lammps_endpoint",
    "lammps_neb",
)
SOLVER_UNITS = (
    "cp2k_initial",
    "cp2k_endpoint",
    "cp2k_band",
    "lammps_initial",
    "lammps_endpoint",
    "lammps_neb",
)
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


def require_outside_tree(
    target: pathlib.Path, protected_root: pathlib.Path, *, label: str
) -> None:
    target_resolved = target.resolve()
    protected_resolved = protected_root.resolve()
    if (
        target_resolved == protected_resolved
        or protected_resolved in target_resolved.parents
    ):
        raise ValueError(f"{label} must live outside the prepared tree")


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
    for start_x in range(SOURCE_SUPERCELL[0]):
        for start_y in range(SOURCE_SUPERCELL[1]):
            tiles = {
                (
                    (start_x + dx) % SOURCE_SUPERCELL[0],
                    (start_y + dy) % SOURCE_SUPERCELL[1],
                )
                for dx in range(TARGET_SUPERCELL[0])
                for dy in range(TARGET_SUPERCELL[1])
            }
            selected = [atom for tile in tiles for atom in by_tile[tile]]
            ids = {atom.id for atom in selected}
            counts = Counter(atom.element for atom in selected)
            if (
                len(selected)
                == SOURCE_ATOMS_PER_UNIT_CELL
                * TARGET_SUPERCELL[0]
                * TARGET_SUPERCELL[1]
                and route_ids <= ids
                and counts == TARGET_PRISTINE_COUNTS
                and math.isclose(contract.net_charge(selected), 0.0, abs_tol=1.0e-8)
            ):
                candidates.append((start_x, start_y, selected))
    if len(candidates) != 1:
        windows = [(item[0], item[1]) for item in candidates]
        raise ValueError(
            f"expected one neutral {TARGET_SUPERCELL[0]}x{TARGET_SUPERCELL[1]} "
            f"route window, found {windows}"
        )
    start_x, start_y, selected = candidates[0]
    target_cell = contract.nteme_neb.Cell(
        a=cell.a * TARGET_SUPERCELL[0] / SOURCE_SUPERCELL[0],
        b=cell.b * TARGET_SUPERCELL[1] / SOURCE_SUPERCELL[1],
        c=cell.c,
        alpha=cell.alpha,
        beta=cell.beta,
        gamma=cell.gamma,
    )
    origin_x = (grid_x + start_x / SOURCE_SUPERCELL[0]) % 1.0
    origin_y = (grid_y + start_y / SOURCE_SUPERCELL[1]) % 1.0
    reduced = []
    for atom in selected:
        source_fractional = fractional[atom.id]
        target_fractional = (
            ((source_fractional[0] - origin_x) % 1.0)
            / (TARGET_SUPERCELL[0] / SOURCE_SUPERCELL[0]),
            ((source_fractional[1] - origin_y) % 1.0)
            / (TARGET_SUPERCELL[1] / SOURCE_SUPERCELL[1]),
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
        raise ValueError("reduced-cell remap created an atomic collision")
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
        raise ValueError(
            f"selected {TARGET_SUPERCELL[0]}x{TARGET_SUPERCELL[1]} cell "
            "does not contain the pinned E3a route"
        )
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
                "reconstructed-charge-compensation in the neutral "
                f"{TARGET_SUPERCELL[0]}x{TARGET_SUPERCELL[1]} route cell; "
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
                "This is a neutral reconstructed "
                f"{TARGET_SUPERCELL[0]}x{TARGET_SUPERCELL[1]} calibration cell, not the "
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
                "The DFT result can test the E3a Xe screen, but one "
                f"{TARGET_SUPERCELL[0]}x{TARGET_SUPERCELL[1]} route "
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


def _coordination_gate(
    model: contract.Model, parent_model: contract.Model | None = None
) -> dict[str, Any]:
    if model.name != "dehydroxylate-lattice":
        return {"coordination_pass": True, "five_coordinate_al_ids": []}
    if parent_model is None or parent_model.name != "reconstructed-replication":
        raise ValueError(
            "dehydroxylate coordination gate requires its parent route model"
        )
    transformation = model.metadata.get("transformation")
    if not isinstance(transformation, Mapping):
        raise TypeError("dehydroxylate model has no typed transformation")
    selected_oxygen = {
        int(value) for value in transformation.get("selected_oh_oxygen_ids", [])
    }
    residual_oxygen = int(transformation["residual_oxygen_site_id"])
    removed_oxygen = sorted(selected_oxygen - {residual_oxygen})
    if len(removed_oxygen) != 1:
        raise ValueError("dehydroxylation must identify one removed hydroxyl oxygen")
    declared_removed_sites = {
        int(value) for value in transformation.get("removed_site_ids", [])
    }
    removed_oxygen_declared = all(
        oxygen_id in declared_removed_sites for oxygen_id in removed_oxygen
    )

    parent_by_id = {atom.id: atom for atom in parent_model.atoms}
    transformed_by_id = {atom.id: atom for atom in model.atoms}
    parent_oxygen = [
        atom for atom in parent_model.atoms if atom.element.startswith("O")
    ]
    transformed_oxygen = [atom for atom in model.atoms if atom.element.startswith("O")]
    parent_coordination = {
        atom.id: sum(
            contract.periodic_distance(atom, oxygen, parent_model.cell)
            <= AL_O_CUTOFF_ANGSTROM
            for oxygen in parent_oxygen
        )
        for atom in parent_model.atoms
        if atom.element == "Al1"
    }
    transformed_coordination = {
        atom.id: sum(
            contract.periodic_distance(atom, oxygen, model.cell) <= AL_O_CUTOFF_ANGSTROM
            for oxygen in transformed_oxygen
        )
        for atom in model.atoms
        if atom.element == "Al1"
    }
    candidates = sorted(
        atom_id for atom_id, count in transformed_coordination.items() if count == 5
    )
    unexpected = {
        atom_id: count
        for atom_id, count in transformed_coordination.items()
        if count not in {5, 6}
    }

    route = model.metadata.get("route")
    if not isinstance(route, Mapping):
        raise TypeError("dehydroxylate model has no typed route")
    moving_id = int(route["moving_site_id"])
    if model.endpoint_atoms is None:
        raise ValueError("dehydroxylate route has no endpoint")
    initial_moving = transformed_by_id[moving_id]
    final_moving = next(atom for atom in model.endpoint_atoms if atom.id == moving_id)
    route_positions = [
        contract.nteme_neb.Atom(
            moving_id,
            initial_moving.element,
            *(
                getattr(initial_moving, axis)
                + index
                / (IMAGE_COUNT - 1)
                * (getattr(final_moving, axis) - getattr(initial_moving, axis))
                for axis in ("x", "y", "z")
            ),
        )
        for index in range(IMAGE_COUNT)
    ]
    candidate_records = []
    for atom_id in candidates:
        parent_al = parent_by_id[atom_id]
        transformed_al = transformed_by_id[atom_id]
        removed_distances = {
            str(oxygen_id): contract.periodic_distance(
                parent_al, parent_by_id[oxygen_id], parent_model.cell
            )
            for oxygen_id in removed_oxygen
        }
        route_distances = [
            contract.periodic_distance(position, transformed_al, model.cell)
            for position in route_positions
        ]
        interaction_images = [
            index
            for index, distance in enumerate(route_distances)
            if distance <= LOCAL_ROUTE_CUTOFF_ANGSTROM
        ]
        parent_count = parent_coordination.get(atom_id)
        transformed_count = transformed_coordination[atom_id]
        parent_neighbor_ids = {
            oxygen.id
            for oxygen in parent_oxygen
            if contract.periodic_distance(parent_al, oxygen, parent_model.cell)
            <= AL_O_CUTOFF_ANGSTROM
        }
        transformed_neighbor_ids = {
            oxygen.id
            for oxygen in transformed_oxygen
            if contract.periodic_distance(transformed_al, oxygen, model.cell)
            <= AL_O_CUTOFF_ANGSTROM
        }
        lost_neighbor_ids = parent_neighbor_ids - transformed_neighbor_ids
        gained_neighbor_ids = transformed_neighbor_ids - parent_neighbor_ids
        exact_neighbor_change = lost_neighbor_ids == set(removed_oxygen) and not (
            gained_neighbor_ids
        )
        removed_in_shell = any(
            distance <= AL_O_CUTOFF_ANGSTROM for distance in removed_distances.values()
        )
        removed_absent = all(
            oxygen_id not in transformed_by_id for oxygen_id in removed_oxygen
        )
        candidate_records.append(
            {
                "al_source_id": atom_id,
                "parent_o_coordination": parent_count,
                "transformed_o_coordination": transformed_count,
                "removed_hydroxyl_o_distances_angstrom": removed_distances,
                "removed_hydroxyl_o_in_parent_shell": removed_in_shell,
                "removed_hydroxyl_o_declared_removed": removed_oxygen_declared,
                "removed_hydroxyl_o_absent_from_transformed_model": removed_absent,
                "parent_o_neighbor_source_ids": sorted(parent_neighbor_ids),
                "transformed_o_neighbor_source_ids": sorted(transformed_neighbor_ids),
                "lost_o_neighbor_source_ids": sorted(lost_neighbor_ids),
                "gained_o_neighbor_source_ids": sorted(gained_neighbor_ids),
                "exact_removed_o_neighbor_change": exact_neighbor_change,
                "coordination_change_caused_by_transformation": (
                    parent_count == 6
                    and transformed_count == 5
                    and removed_in_shell
                    and removed_oxygen_declared
                    and removed_absent
                    and exact_neighbor_change
                ),
                "route_image_distances_angstrom": route_distances,
                "minimum_route_distance_angstrom": min(route_distances),
                "interaction_image_indices_zero_based": interaction_images,
                "route_interaction_pass": bool(interaction_images),
            }
        )
    criterion = {
        "schema": "e3b-five-coordinate-al-criterion-v1",
        "definition": (
            "An Al involves the dehydroxylation route only when it changes from six "
            "parent O neighbors to five transformed O neighbors with the removed "
            "hydroxyl O as the only lost neighbor and no gained neighbor, and at least "
            "one explicit noble-gas route image enters the declared Al interaction shell."
        ),
        "parent_model": parent_model.name,
        "transformed_model": model.name,
        "al_o_cutoff_angstrom": AL_O_CUTOFF_ANGSTROM,
        "route_interaction_cutoff_angstrom": LOCAL_ROUTE_CUTOFF_ANGSTROM,
        "route_position_source": (
            "moving atom coordinates from initial plus six linearly interpolated "
            "replicas plus final endpoint"
        ),
        "route_moving_atom_source_id": moving_id,
        "removed_hydroxyl_oxygen_ids": removed_oxygen,
        "candidate_al_source_ids": candidates,
        "candidates": candidate_records,
    }
    passed = (
        candidates == [10, 11]
        and not unexpected
        and all(
            record["coordination_change_caused_by_transformation"] is True
            and record["route_interaction_pass"] is True
            for record in candidate_records
        )
    )
    return {
        "coordination_pass": passed,
        "five_coordinate_al_ids": candidates,
        "unexpected_al_coordination": unexpected,
        "five_coordinate_al_criterion": criterion,
    }


def validate_spot_model(
    model: contract.Model,
    source_cell: contract.nteme_neb.Cell,
    parent_model: contract.Model | None = None,
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
            math.isclose(
                model.cell.a,
                source_cell.a * TARGET_SUPERCELL[0] / SOURCE_SUPERCELL[0],
                abs_tol=1.0e-10,
            ),
            math.isclose(
                model.cell.b,
                source_cell.b * TARGET_SUPERCELL[1] / SOURCE_SUPERCELL[1],
                abs_tol=1.0e-10,
            ),
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
    coordination = _coordination_gate(model, parent_model)
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


def _classical_reference(
    campaign: Mapping[str, Any],
    name: str,
    campaign_root: pathlib.Path,
    campaign_manifest: Mapping[str, Any],
) -> dict[str, Any]:
    models = campaign["models"]
    if not isinstance(models, dict) or not isinstance(models[name], dict):
        raise TypeError(f"E3a campaign has malformed model {name}")
    model_record = models[name]
    neb = model_record.get("neb")
    if not isinstance(neb, dict):
        raise TypeError(f"E3a campaign model {name} has no NEB result")
    barrier = neb.get("computed_forward_barrier_kcal_mol")
    if not isinstance(barrier, (int, float)) or not math.isfinite(float(barrier)):
        raise ValueError(f"E3a campaign model {name} has no finite barrier")
    settings = campaign.get("settings")
    if not isinstance(settings, Mapping):
        raise TypeError("E3a campaign has no typed settings")
    requested_ftol = settings.get("neb_ftol_kcal_mol_angstrom")
    maximum_force = neb.get("max_replica_force_kcal_mol_angstrom")
    if not isinstance(requested_ftol, (int, float)) or not isinstance(
        maximum_force, (int, float)
    ):
        raise TypeError(f"E3a campaign model {name} has malformed convergence fields")
    if not math.isfinite(float(requested_ftol)) or not math.isfinite(
        float(maximum_force)
    ):
        raise ValueError(f"E3a campaign model {name} has non-finite convergence fields")
    manifest_files = campaign_manifest.get("files")
    if not isinstance(manifest_files, list):
        raise TypeError("E3a campaign manifest has no typed files")

    def unique_manifest_entry(
        *, path: str | None = None, digest: str | None = None
    ) -> dict[str, Any]:
        matches = [
            entry
            for entry in manifest_files
            if isinstance(entry, dict)
            and (path is None or entry.get("path") == path)
            and (digest is None or entry.get("sha256") == digest)
        ]
        if len(matches) != 1:
            raise ValueError(
                f"E3a manifest does not bind exactly one requested artifact: {path or digest}"
            )
        return dict(matches[0])

    campaign_result = unique_manifest_entry(path="campaign-result.json")
    raw_neb = unique_manifest_entry(digest=str(neb.get("screen_sha256")))
    raw_endpoints = {
        label: unique_manifest_entry(path=f"{name}/{label}.min.stdout.log")
        for label in ("initial", "endpoint")
    }
    expected_prefix = f"{name}/"
    if not str(raw_neb.get("path", "")).startswith(expected_prefix):
        raise ValueError(f"E3a model {name} NEB output belongs to another model")
    bound_artifacts = [
        ("campaign result", campaign_result),
        ("raw NEB", raw_neb),
        *((f"raw {label} endpoint", entry) for label, entry in raw_endpoints.items()),
    ]
    for label, entry in bound_artifacts:
        relative = entry.get("path")
        digest = entry.get("sha256")
        size = entry.get("bytes")
        if (
            not isinstance(relative, str)
            or not isinstance(digest, str)
            or not isinstance(size, int)
        ):
            raise TypeError(f"E3a {label} manifest record is malformed")
        artifact = campaign_root / relative
        require_sha256(artifact, digest, label=f"E3a {label}")
        if artifact.stat().st_size != size:
            raise ValueError(f"E3a {label} size mismatch")
    raw_neb_path = campaign_root / str(raw_neb["path"])
    independently_parsed_neb = e3a_campaign.parse_neb(
        raw_neb_path, float(requested_ftol)
    )
    parsed_barrier = independently_parsed_neb.get("computed_forward_barrier_kcal_mol")
    if not isinstance(parsed_barrier, (int, float)) or not math.isclose(
        float(parsed_barrier), float(barrier), rel_tol=0.0, abs_tol=1.0e-9
    ):
        raise ValueError(f"E3a model {name} campaign barrier does not match raw NEB")
    endpoint_ftol = settings.get("min_ftol_kcal_mol_angstrom")
    if not isinstance(endpoint_ftol, (int, float)) or not math.isfinite(
        float(endpoint_ftol)
    ):
        raise TypeError("E3a campaign endpoint force tolerance is malformed")
    independently_parsed_endpoints = {
        label: e3a_campaign.parse_last_thermo(campaign_root / str(entry["path"]))
        for label, entry in raw_endpoints.items()
    }
    endpoint_normal_termination = {
        label: (
            "total wall time:"
            in (campaign_root / str(entry["path"])).read_text(encoding="utf-8").lower()
            or "loop time of"
            in (campaign_root / str(entry["path"])).read_text(encoding="utf-8").lower()
        )
        for label, entry in raw_endpoints.items()
    }
    endpoint_force_pass = {
        label: math.isfinite(parsed["force_norm_kcal_mol_angstrom"])
        and parsed["force_norm_kcal_mol_angstrom"] <= float(endpoint_ftol)
        and endpoint_normal_termination[label]
        for label, parsed in independently_parsed_endpoints.items()
    }
    independent_maximum = independently_parsed_neb.get(
        "max_replica_force_kcal_mol_angstrom"
    )
    if not isinstance(independent_maximum, (int, float)) or not math.isclose(
        float(independent_maximum), float(maximum_force), rel_tol=0.0, abs_tol=1.0e-12
    ):
        raise ValueError(f"E3a model {name} campaign force does not match raw NEB")
    if (independently_parsed_neb.get("converged_to_requested_ftol") is True) != (
        neb.get("converged_to_requested_ftol") is True
    ):
        raise ValueError(
            f"E3a model {name} campaign convergence does not match raw NEB"
        )
    expected_final_step = int(settings.get("neb_relax_steps", -1)) + int(
        settings.get("neb_climb_steps", -1)
    )
    final_step = independently_parsed_neb.get("final_step")
    complete_step_history = (
        isinstance(final_step, int) and final_step == expected_final_step
    )
    independently_converged = (
        independently_parsed_neb.get("converged_to_requested_ftol") is True
        and math.isfinite(float(independent_maximum))
        and float(independent_maximum) <= float(requested_ftol)
        and complete_step_history
        and all(endpoint_force_pass.values())
    )
    return {
        "schema": "e3a-source-barrier-reference-v1",
        "barrier_kcal_mol": float(barrier),
        "status": "converged" if independently_converged else "incomplete-convergence",
        "source_cell": [6, 3, 1],
        "campaign_result": campaign_result,
        "campaign_model_record_sha256": canonical_sha256(model_record),
        "raw_neb_output": raw_neb,
        "raw_endpoint_outputs": raw_endpoints,
        "convergence": {
            "independent_parser": "e3a_campaign.parse_last_thermo+parse_neb",
            "requested_ftol_kcal_mol_angstrom": float(requested_ftol),
            "max_replica_force_kcal_mol_angstrom": float(independent_maximum),
            "converged_to_requested_ftol": independently_parsed_neb.get(
                "converged_to_requested_ftol"
            )
            is True,
            "final_step": independently_parsed_neb.get("final_step"),
            "expected_final_step": expected_final_step,
            "complete_step_history": complete_step_history,
            "parsed_forward_barrier_kcal_mol": float(parsed_barrier),
            "endpoint_requested_ftol_kcal_mol_angstrom": float(endpoint_ftol),
            "initial_endpoint_force_pass": endpoint_force_pass["initial"],
            "final_endpoint_force_pass": endpoint_force_pass["endpoint"],
            "initial_endpoint_normal_termination": endpoint_normal_termination[
                "initial"
            ],
            "final_endpoint_normal_termination": endpoint_normal_termination[
                "endpoint"
            ],
            "initial_endpoint_force_norm_kcal_mol_angstrom": (
                independently_parsed_endpoints["initial"][
                    "force_norm_kcal_mol_angstrom"
                ]
            ),
            "final_endpoint_force_norm_kcal_mol_angstrom": (
                independently_parsed_endpoints["endpoint"][
                    "force_norm_kcal_mol_angstrom"
                ]
            ),
        },
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
    campaign_manifest = read_json(campaign_root / "manifest.json")
    prepare_output_dir(out)
    method_record = materialize_method_files(out, method)
    source_atoms, source_cell, barriers = contract.nteme_neb.read_workbook(workbook)
    route = contract.nteme_neb.select_route(barriers, "divacancy", 1)
    pristine, target_cell, reduction = reduce_pristine(source_atoms, source_cell)
    models = build_spot_models(pristine, target_cell, route)
    models_by_name = {model.name: model for model in models}
    model_records: dict[str, Any] = {}
    models_root = out / "models"
    models_root.mkdir()
    for model in models:
        gates = validate_spot_model(
            model,
            source_cell,
            models_by_name["reconstructed-replication"],
        )
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
            str(path.relative_to(model_dir)): sha256(path)
            for path in sorted(item for item in model_dir.rglob("*") if item.is_file())
        }
        model_records[model.name] = {
            "scientific_scope": SCIENTIFIC_SCOPE[model.name],
            "atom_count": len(model.atoms),
            "source_type_counts": dict(sorted(source_counts.items())),
            "chemical_element_counts": dict(sorted(chemical_counts.items())),
            "supercell": list(TARGET_SUPERCELL),
            "periodic_cell": asdict(model.cell),
            "periodic_cell_sha256": canonical_sha256(asdict(model.cell)),
            "route": model.metadata["route"],
            "transformation": model.metadata.get("transformation"),
            "limitation": model.metadata["limitation"],
            "gates": gates,
            "atom_identity_sha256": identity_sha,
            "classical_reference": _classical_reference(
                campaign,
                model.name,
                campaign_root,
                campaign_manifest,
            ),
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
    raw_output_template = {
        "schema": "e3b-raw-output-spec-v1",
        "models": {
            name: {
                "execution_receipt": None,
                "outputs": {role: None for role in RAW_OUTPUT_ROLES},
            }
            for name in model_records
        },
    }
    write_json(out / "raw-outputs.template.json", raw_output_template)
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
    if energy is not None and not math.isfinite(energy):
        energy = None
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


def _cp2k_float(token: str, *, label: str) -> float:
    try:
        value = float(token.replace("D", "E").replace("d", "e"))
    except ValueError as exc:
        raise ValueError(f"{label} is malformed") from exc
    if not math.isfinite(value):
        raise ValueError(f"{label} is not finite")
    return value


def _cp2k_expected_vectors(cell: Mapping[str, Any]) -> dict[str, tuple[float, ...]]:
    try:
        value = contract.nteme_neb.Cell(
            **{
                key: float(cell[key])
                for key in ("a", "b", "c", "alpha", "beta", "gamma")
            }
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("prepared periodic cell record is malformed") from exc
    lx, ly, lz, xy, xz, yz = value.restricted()
    return {"a": (lx, 0.0, 0.0), "b": (xy, ly, 0.0), "c": (xz, yz, lz)}


def _validate_cp2k_identity(
    text: str,
    *,
    label: str,
    expected_project: str | None,
    expected_atom_count: int | None,
    expected_cell: Mapping[str, Any] | None,
    project_prefix: bool = False,
) -> dict[str, Any]:
    fatal = re.search(r"(?im)^\s*(?:\*+\s*)?(?:ERROR|FATAL|ABORT)\b", text)
    if fatal:
        raise ValueError(f"{label}: CP2K fatal/error output")
    projects = re.findall(r"(?im)^\s*GLOBAL\|\s*Project name\s+(\S+)\s*$", text)
    if expected_project is not None:
        if not projects:
            raise ValueError(f"{label}: missing CP2K project identity")
        actual_project = projects[-1]
        project_matches = (
            actual_project.startswith(expected_project + "-BAND")
            if project_prefix
            else actual_project == expected_project
        )
        if not project_matches:
            raise ValueError(
                f"{label}: CP2K project {actual_project!r} does not match {expected_project!r}"
            )
    else:
        actual_project = projects[-1] if projects else None

    atom_matches = re.findall(r"(?im)^\s*(?:GLOBAL\|\s*)?-\s*Atoms:\s*(\d+)\s*$", text)
    if expected_atom_count is not None:
        if not atom_matches:
            raise ValueError(f"{label}: missing CP2K atom count")
        if int(atom_matches[-1]) != expected_atom_count:
            raise ValueError(f"{label}: CP2K atom count does not match preparation")

    parsed_vectors: dict[str, tuple[float, float, float]] = {}
    vector_pattern = re.compile(
        r"(?im)^\s*CELL\|\s*Vector\s+([abc])\s+\[angstrom\]:\s*"
        r"([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)\s+([-+0-9.eEdD]+)"
    )
    for match in vector_pattern.finditer(text):
        parsed_vectors[match.group(1).lower()] = tuple(
            _cp2k_float(match.group(index), label=f"{label} cell vector")
            for index in range(2, 5)
        )
    if expected_cell is not None:
        if set(parsed_vectors) != {"a", "b", "c"}:
            raise ValueError(f"{label}: incomplete CP2K periodic cell identity")
        expected_vectors = _cp2k_expected_vectors(expected_cell)
        if any(
            not math.isclose(actual, expected, rel_tol=0.0, abs_tol=1.0e-8)
            for axis in ("a", "b", "c")
            for actual, expected in zip(
                parsed_vectors[axis], expected_vectors[axis], strict=True
            )
        ):
            raise ValueError(f"{label}: CP2K periodic cell does not match preparation")
    return {
        "project": actual_project,
        "atom_count": int(atom_matches[-1]) if atom_matches else None,
        "cell_vectors_angstrom": parsed_vectors or None,
    }


def _parse_cp2k_profile(path: pathlib.Path, expected_replicas: int) -> list[float]:
    text = path.read_text(encoding="utf-8")
    native_profiles: list[list[float]] = []
    marker = re.compile(r"ENERGIES\s*\[au\]\s*=", re.IGNORECASE)
    number = re.compile(
        r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][-+]?\d+)?|[-+]?nan|[-+]?inf",
        re.IGNORECASE,
    )
    for match in marker.finditer(text):
        values: list[float] = []
        remainder = text[match.end() :]
        for line in remainder.splitlines():
            tokens = line.split()
            if not tokens:
                continue
            if not all(number.fullmatch(token) for token in tokens):
                break
            values.extend(
                _cp2k_float(token, label="CP2K BAND profile energy") for token in tokens
            )
            if len(values) >= expected_replicas:
                break
        native_profiles.append(values)
    if native_profiles:
        profile = native_profiles[-1]
    else:
        numeric_rows: list[list[float]] = []
        for line in text.splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            tokens = stripped.split()
            if not all(number.fullmatch(token) for token in tokens):
                continue
            numeric_rows.append(
                [
                    _cp2k_float(token, label="CP2K .ener profile value")
                    for token in tokens
                ]
            )
        if len(numeric_rows) != expected_replicas:
            raise ValueError(
                "CP2K final native profile has wrong replica count: "
                f"expected {expected_replicas}, parsed {len(numeric_rows)}"
            )
        profile = [row[4] if len(row) >= 5 else row[-1] for row in numeric_rows]
    if len(profile) != expected_replicas:
        raise ValueError(
            "CP2K final native profile has wrong replica count: "
            f"expected {expected_replicas}, parsed {len(profile)}"
        )
    if not all(math.isfinite(value) for value in profile):
        raise ValueError("CP2K final native profile contains non-finite energy")
    return profile


def parse_cp2k_barrier(
    *,
    initial_log: pathlib.Path,
    endpoint_log: pathlib.Path,
    band_log: pathlib.Path,
    profile_log: pathlib.Path,
    replica_logs: Sequence[pathlib.Path],
    expected_replicas: int = IMAGE_COUNT,
    expected_model: str | None = None,
    expected_atom_count: int | None = None,
    expected_cell: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Parse CP2K-2024.3 endpoint and native CI-NEB output fail-closed."""

    project_base = f"e3b-{expected_model}" if expected_model is not None else None

    def parse_endpoint(path: pathlib.Path, label: str, suffix: str) -> dict[str, Any]:
        text = path.read_text(encoding="utf-8")
        lower = text.lower()
        identity = _validate_cp2k_identity(
            text,
            label=label,
            expected_project=f"{project_base}-{suffix}" if project_base else None,
            expected_atom_count=expected_atom_count,
            expected_cell=expected_cell,
        )
        if "program ended at" not in lower:
            raise ValueError(f"{label}: missing CP2K normal termination")
        if "scf run not converged" in lower or "scf run did not converge" in lower:
            raise ValueError(f"{label}: SCF did not converge")
        if not re.search(r"scf run converged in\s+\d+\s+steps", lower):
            raise ValueError(f"{label}: missing converged SCF evidence")
        if not re.search(r"geometry optimization (?:completed|converged)", lower):
            raise ValueError(f"{label}: missing endpoint geometry convergence")
        classified = classify_cp2k_output(text, returncode=0, timed_out=False)
        energy = classified.get("energy_hartree")
        if not isinstance(energy, (int, float)) or not math.isfinite(float(energy)):
            raise ValueError(f"{label}: endpoint energy is not finite")
        return {
            "energy_hartree": float(energy),
            "scf_converged": True,
            "geometry_converged": True,
            "normal_termination": True,
            "identity": identity,
        }

    initial = parse_endpoint(initial_log, "initial endpoint", "initial")
    endpoint = parse_endpoint(endpoint_log, "final endpoint", "endpoint")
    band_text = band_log.read_text(encoding="utf-8")
    band_lower = band_text.lower()
    band_identity = _validate_cp2k_identity(
        band_text,
        label="BAND output",
        expected_project=f"{project_base}-neb" if project_base else None,
        expected_atom_count=expected_atom_count,
        expected_cell=expected_cell,
    )
    if "program ended at" not in band_lower:
        raise ValueError("BAND output is missing CP2K normal termination")
    count_matches = re.findall(
        r"NUMBER OF NEB REPLICA\s*=\s*(\d+)", band_text, re.IGNORECASE
    )
    if not count_matches or int(count_matches[-1]) != expected_replicas:
        raise ValueError("CP2K BAND replica count does not match preparation")

    criteria: dict[str, bool] = {}
    for key, pattern in (
        ("rms_displacement", r"RMS\s+DISPLACEMENT"),
        ("max_displacement", r"MAX(?:IMUM)?\s+DISPLACEMENT"),
        ("rms_force", r"RMS\s+FORCE"),
        ("max_force", r"MAX(?:IMUM)?\s+FORCE"),
    ):
        matches = re.findall(rf"(?im)^.*{pattern}.*\((YES|NO)\)\s*$", band_text)
        if not matches:
            raise ValueError(f"CP2K BAND final convergence criteria incomplete: {key}")
        criteria[key] = matches[-1].upper() == "YES"
    if not all(criteria.values()):
        raise ValueError("CP2K BAND final convergence criteria did not all pass")

    if len(replica_logs) != expected_replicas:
        raise ValueError(
            "CP2K BAND replica log count does not match preparation: "
            f"expected {expected_replicas}, supplied {len(replica_logs)}"
        )
    replica_records = []
    for index, path in enumerate(replica_logs, 1):
        text = path.read_text(encoding="utf-8")
        lower = text.lower()
        identity = _validate_cp2k_identity(
            text,
            label=f"replica {index}",
            expected_project=f"{project_base}-neb" if project_base else None,
            expected_atom_count=expected_atom_count,
            expected_cell=expected_cell,
            project_prefix=project_base is not None,
        )
        if "scf run not converged" in lower or "scf run did not converge" in lower:
            raise ValueError(f"replica {index}: SCF did not converge")
        converged_cycles = len(re.findall(r"scf run converged in\s+\d+\s+steps", lower))
        if converged_cycles < 1:
            raise ValueError(f"replica {index}: missing final converged SCF evidence")
        normal_termination = "program ended at" in lower
        if not normal_termination:
            raise ValueError(f"replica {index}: missing CP2K normal termination")
        replica_records.append(
            {
                "replica": index,
                "scf_converged": True,
                "scf_converged_cycles": converged_cycles,
                "normal_termination": True,
                "identity": identity,
            }
        )

    profile = _parse_cp2k_profile(profile_log, expected_replicas)
    if not math.isclose(
        profile[0], initial["energy_hartree"], rel_tol=0.0, abs_tol=1.0e-6
    ) or not math.isclose(
        profile[-1], endpoint["energy_hartree"], rel_tol=0.0, abs_tol=1.0e-6
    ):
        raise ValueError("CP2K BAND endpoint energies do not match endpoint logs")
    barrier_hartree = max(profile) - profile[0]
    barrier_kcal_mol = barrier_hartree * HARTREE_TO_KCAL_MOL
    if barrier_hartree < 0 or not math.isfinite(barrier_kcal_mol):
        raise ValueError("CP2K barrier is not finite and non-negative")
    return {
        "status": "converged",
        "expected_replicas": expected_replicas,
        "parsed_replicas": len(profile),
        "replica_outputs": replica_records,
        "initial_endpoint": initial,
        "final_endpoint": endpoint,
        "band_identity": band_identity,
        "final_convergence_criteria": criteria,
        "replica_energies_hartree": profile,
        "barrier_hartree": barrier_hartree,
        "hartree_to_kcal_mol": HARTREE_TO_KCAL_MOL,
        "barrier_kcal_mol": barrier_kcal_mol,
        "band_converged": True,
        "normal_termination": True,
    }


def parse_matched_classical_barrier(
    *,
    initial_log: pathlib.Path,
    endpoint_log: pathlib.Path,
    neb_log: pathlib.Path,
    expected_replicas: int = IMAGE_COUNT,
    ftol: float = MATCHED_CLASSICAL_FTOL_KCAL_MOL_ANGSTROM,
) -> dict[str, Any]:
    """Parse converged matched-cell LAMMPS endpoints and NEB screen output."""

    def parse_endpoint(path: pathlib.Path, label: str) -> dict[str, float]:
        text = path.read_text(encoding="utf-8")
        if re.search(r"(?im)^\s*(?:ERROR|FATAL)\b", text):
            raise ValueError(f"{label}: LAMMPS fatal/error output")
        if (
            "total wall time:" not in text.lower()
            and "loop time of" not in text.lower()
        ):
            raise ValueError(f"{label}: missing LAMMPS normal termination")
        parsed = e3a_campaign.parse_last_thermo(path)
        numeric = [float(value) for value in parsed.values()]
        if not all(math.isfinite(value) for value in numeric):
            raise ValueError(f"{label}: non-finite endpoint thermo value")
        if parsed["force_norm_kcal_mol_angstrom"] > ftol:
            raise ValueError(f"{label}: endpoint force did not converge")
        return parsed

    initial = parse_endpoint(initial_log, "initial endpoint")
    endpoint = parse_endpoint(endpoint_log, "final endpoint")
    neb_text = neb_log.read_text(encoding="utf-8")
    if re.search(r"(?im)^\s*(?:ERROR|FATAL)\b", neb_text):
        raise ValueError("LAMMPS NEB output contains fatal/error output")
    lines = neb_text.splitlines()
    climb_markers = [
        index
        for index, line in enumerate(lines)
        if re.search(r"Setting up climbing|Climbing replica\s*=", line, re.IGNORECASE)
    ]
    if not climb_markers:
        raise ValueError("LAMMPS NEB output never entered climbing phase")
    climb_start = climb_markers[0]
    header_indices = [
        index
        for index, line in enumerate(lines)
        if re.search(r"^\s*Step\s+MaxReplicaForce\s+MaxAtomForce\b", line)
    ]
    climb_headers = [index for index in header_indices if index > climb_start]
    if len(header_indices) < 2 or not climb_headers:
        raise ValueError("LAMMPS NEB output lacks the climbing progress header")
    second_header = climb_headers[0]
    climb_rows: list[tuple[int, list[float]]] = []
    for index, line in enumerate(lines[second_header + 1 :], second_header + 1):
        fields = line.split()
        if len(fields) < 9:
            continue
        try:
            values = [float(field) for field in fields]
        except ValueError:
            continue
        climb_rows.append((index, values))
    if not climb_rows:
        raise ValueError("LAMMPS NEB climbing phase has no numeric progress row")
    final_row_index, final_row = climb_rows[-1]
    field_count = len(final_row)
    expected_fields = 9 + 2 * expected_replicas
    if field_count != expected_fields:
        parsed_replicas = (field_count - 9) // 2 if field_count >= 9 else 0
        raise ValueError(
            "LAMMPS NEB replica count does not match preparation: "
            f"expected {expected_replicas}, parsed {parsed_replicas}"
        )
    if not all(math.isfinite(value) for value in final_row):
        raise ValueError("LAMMPS NEB contains a non-finite value")
    footer_after_final = "\n".join(lines[final_row_index + 1 :]).lower()
    if (
        "total wall time:" not in footer_after_final
        and "loop time of" not in footer_after_final
    ):
        raise ValueError("LAMMPS NEB output is missing normal termination after climb")
    if final_row[1] > ftol:
        raise ValueError("LAMMPS NEB did not converge to the requested force tolerance")
    barrier = final_row[6]
    if not math.isfinite(barrier):
        raise ValueError("LAMMPS NEB barrier is not finite")
    neb_initial_energy = final_row[10]
    neb_endpoint_energy = final_row[8 + 2 * expected_replicas]
    endpoint_deltas = {
        "initial": abs(neb_initial_energy - initial["potential_energy_kcal_mol"]),
        "endpoint": abs(neb_endpoint_energy - endpoint["potential_energy_kcal_mol"]),
    }
    if any(
        delta > LAMMPS_ENDPOINT_ENERGY_TOLERANCE_KCAL_MOL
        for delta in endpoint_deltas.values()
    ):
        raise ValueError("LAMMPS NEB endpoint energies do not match minimization logs")
    neb = {
        "final_step": int(final_row[0]),
        "max_replica_force_kcal_mol_angstrom": final_row[1],
        "max_atom_force_kcal_mol_angstrom": final_row[2],
        "computed_forward_barrier_kcal_mol": barrier,
        "computed_reverse_barrier_kcal_mol": final_row[7],
        "converged_to_requested_ftol": True,
        "screen_sha256": sha256(neb_log),
    }
    return {
        "status": "converged",
        "expected_replicas": expected_replicas,
        "parsed_replicas": expected_replicas,
        "ftol_kcal_mol_angstrom": ftol,
        "initial_endpoint": initial,
        "final_endpoint": endpoint,
        "neb": neb,
        "barrier_kcal_mol": float(barrier),
        "climbing_phase": True,
        "climbing_progress_header_line": second_header + 1,
        "endpoint_energy_tolerance_kcal_mol": (
            LAMMPS_ENDPOINT_ENERGY_TOLERANCE_KCAL_MOL
        ),
        "endpoint_energy_deltas_kcal_mol": endpoint_deltas,
        "normal_termination": True,
    }


def _file_record(path: pathlib.Path, *, relative_to: pathlib.Path) -> dict[str, Any]:
    return {
        "path": str(path.relative_to(relative_to)),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def _absolute_file_record(path: pathlib.Path) -> dict[str, Any]:
    resolved = path.resolve()
    if not resolved.is_file():
        raise ValueError(f"execution artifact is not a file: {path}")
    return {
        "path": str(resolved),
        "bytes": resolved.stat().st_size,
        "sha256": sha256(resolved),
    }


def _prepared_unit_paths(
    prepared_root: pathlib.Path, model: str
) -> dict[str, list[pathlib.Path]]:
    model_root = prepared_root / "models" / model
    cp2k_data = sorted(
        path for path in (prepared_root / "cp2k-data").glob("*") if path.is_file()
    )
    images = sorted(
        path for path in (model_root / "images").glob("*") if path.is_file()
    )
    matched = model_root / "matched-classical"
    paths = {
        "cp2k_initial": [
            model_root / "initial-opt.inp",
            model_root / "initial.xyz",
            *cp2k_data,
        ],
        "cp2k_endpoint": [
            model_root / "endpoint-opt.inp",
            model_root / "endpoint.xyz",
            *cp2k_data,
        ],
        "cp2k_band": [model_root / "neb.inp", *images, *cp2k_data],
        "lammps_initial": [matched / "in.min.initial", matched / "initial.seed.data"],
        "lammps_endpoint": [
            matched / "in.min.endpoint",
            matched / "endpoint.seed.data",
        ],
        "lammps_neb": [
            matched / "in.neb",
            matched / "initial.seed.data",
            matched / "endpoint.seed.data",
        ],
    }
    for unit, unit_paths in paths.items():
        if not unit_paths or any(not path.is_file() for path in unit_paths):
            missing = [str(path) for path in unit_paths if not path.is_file()]
            raise ValueError(
                f"{model}/{unit}: prepared inputs are incomplete: {missing}"
            )
    return paths


def _resolved_executable_record(
    requested: str | pathlib.Path, *, identity: str
) -> dict[str, Any]:
    if identity not in {"cp2k", "lammps", "wrapper", "test/fake", "unverified"}:
        raise ValueError(f"unsupported executable identity: {identity!r}")
    requested_text = str(requested)
    resolved = pathlib.Path(shutil.which(requested_text) or requested_text).resolve()
    record = _absolute_file_record(resolved)
    return {
        "requested": requested_text,
        "identity": identity,
        "identity_verification": (
            "test-fixture" if identity == "test/fake" else "declared-not-verified"
        ),
        **record,
    }


def create_execution_receipt(
    *,
    prepared_root: pathlib.Path,
    model: str,
    raw_outputs: Mapping[str, pathlib.Path],
    executables: Mapping[str, str | pathlib.Path],
    executable_identities: Mapping[str, str],
    commands: Mapping[str, Sequence[str]],
    resources: Mapping[str, Mapping[str, Any]],
    out: pathlib.Path,
) -> dict[str, Any]:
    """Create a post-execution sidecar binding prepared inputs to raw outputs."""
    verify_manifest(prepared_root)
    preparation_path = prepared_root / "preparation.json"
    preparation = read_json(preparation_path)
    models = preparation.get("models")
    method = preparation.get("method")
    if not isinstance(models, Mapping) or model not in models:
        raise ValueError(f"unknown prepared execution model: {model}")
    prepared_model = models[model]
    if not isinstance(prepared_model, Mapping) or not isinstance(method, Mapping):
        raise TypeError("preparation model/method record is malformed")
    if set(raw_outputs) != set(RAW_OUTPUT_ROLES):
        raise ValueError("execution receipt requires the complete raw-output role set")
    if set(executables) != {"cp2k", "lammps"} or set(executable_identities) != {
        "cp2k",
        "lammps",
    }:
        raise ValueError(
            "execution receipt requires CP2K and LAMMPS executable identities"
        )
    if set(commands) != set(SOLVER_UNITS) or set(resources) != set(SOLVER_UNITS):
        raise ValueError(
            "execution receipt requires command/resources for every solver unit"
        )
    unit_paths = _prepared_unit_paths(prepared_root, model)
    units: dict[str, Any] = {}
    for unit in SOLVER_UNITS:
        command = commands[unit]
        if (
            not isinstance(command, Sequence)
            or isinstance(command, (str, bytes))
            or not command
            or any(not isinstance(token, str) or not token for token in command)
        ):
            raise TypeError(f"{model}/{unit}: command must be a non-empty token list")
        unit_resources = resources[unit]
        if not isinstance(unit_resources, Mapping):
            raise TypeError(f"{model}/{unit}: resources must be a typed mapping")
        units[unit] = {
            "solver": "cp2k" if unit.startswith("cp2k_") else "lammps",
            "command": list(command),
            "resources": dict(unit_resources),
            "prepared_inputs": [
                _file_record(path, relative_to=prepared_root)
                for path in unit_paths[unit]
            ],
        }
    cell = prepared_model.get("periodic_cell")
    if not isinstance(cell, Mapping):
        raise TypeError(f"{model}: prepared periodic cell is missing")
    receipt = {
        "schema": "e3b-model-execution-receipt-v1",
        "preparation": {
            "path": str(preparation_path.resolve()),
            "sha256": sha256(preparation_path),
            "manifest_sha256": sha256(prepared_root / "manifest.json"),
        },
        "model": model,
        "atom_identity_sha256": prepared_model.get("atom_identity_sha256"),
        "periodic_cell_record_sha256": canonical_sha256(cell),
        "method": {
            "record_sha256": canonical_sha256(method),
            "cp2k_version": method.get("cp2k_version"),
            "cp2k_image_digest": method.get("cp2k_image_digest"),
        },
        "executables": {
            solver: _resolved_executable_record(
                executables[solver], identity=executable_identities[solver]
            )
            for solver in ("cp2k", "lammps")
        },
        "units": units,
        "raw_outputs": {
            role: _absolute_file_record(pathlib.Path(raw_outputs[role]))
            for role in RAW_OUTPUT_ROLES
        },
    }
    write_json(out, receipt)
    return receipt


def verify_execution_receipt(
    *,
    receipt_path: pathlib.Path,
    prepared_root: pathlib.Path,
    model: str,
    raw_outputs: Mapping[str, pathlib.Path],
) -> dict[str, Any]:
    receipt = read_json(receipt_path)
    if receipt.get("schema") != "e3b-model-execution-receipt-v1":
        raise ValueError(f"{model}: unexpected execution receipt schema")
    preparation = read_json(prepared_root / "preparation.json")
    models = preparation.get("models")
    method = preparation.get("method")
    if not isinstance(models, Mapping) or not isinstance(method, Mapping):
        raise TypeError("preparation model/method record is malformed")
    prepared_model = models.get(model)
    if not isinstance(prepared_model, Mapping):
        raise ValueError(f"{model}: model is absent from preparation")
    binding = receipt.get("preparation")
    if not isinstance(binding, Mapping) or (
        binding.get("sha256") != sha256(prepared_root / "preparation.json")
        or binding.get("manifest_sha256") != sha256(prepared_root / "manifest.json")
    ):
        raise ValueError(f"{model}: execution receipt preparation binding mismatch")
    cell = prepared_model.get("periodic_cell")
    if (
        receipt.get("model") != model
        or receipt.get("atom_identity_sha256")
        != prepared_model.get("atom_identity_sha256")
        or not isinstance(cell, Mapping)
        or receipt.get("periodic_cell_record_sha256") != canonical_sha256(cell)
    ):
        raise ValueError(f"{model}: execution receipt model identity mismatch")
    receipt_method = receipt.get("method")
    if not isinstance(receipt_method, Mapping) or (
        receipt_method.get("record_sha256") != canonical_sha256(method)
        or receipt_method.get("cp2k_version") != method.get("cp2k_version")
        or receipt_method.get("cp2k_image_digest") != method.get("cp2k_image_digest")
    ):
        raise ValueError(f"{model}: execution receipt method binding mismatch")
    executables = receipt.get("executables")
    if not isinstance(executables, Mapping) or set(executables) != {"cp2k", "lammps"}:
        raise ValueError(f"{model}: executable identity records are incomplete")
    for solver, record in executables.items():
        if not isinstance(record, Mapping):
            raise TypeError(f"{model}/{solver}: malformed executable record")
        executable_path = pathlib.Path(str(record.get("path", "")))
        require_sha256(
            executable_path,
            str(record.get("sha256", "")),
            label=f"{model}/{solver} executable",
        )
        if executable_path.stat().st_size != record.get("bytes"):
            raise ValueError(f"{model}/{solver}: executable size mismatch")
    expected_unit_paths = _prepared_unit_paths(prepared_root, model)
    units = receipt.get("units")
    if not isinstance(units, Mapping) or set(units) != set(SOLVER_UNITS):
        raise ValueError(f"{model}: execution unit records are incomplete")
    for unit, expected_paths in expected_unit_paths.items():
        record = units[unit]
        if not isinstance(record, Mapping):
            raise TypeError(f"{model}/{unit}: malformed execution unit")
        command = record.get("command")
        resources = record.get("resources")
        inputs = record.get("prepared_inputs")
        if (
            not isinstance(command, list)
            or not command
            or not all(isinstance(token, str) and token for token in command)
            or not isinstance(resources, Mapping)
            or not isinstance(inputs, list)
        ):
            raise ValueError(
                f"{model}/{unit}: command/resources/input binding is malformed"
            )
        expected_records = [
            _file_record(path, relative_to=prepared_root) for path in expected_paths
        ]
        if inputs != expected_records:
            raise ValueError(
                f"{model}/{unit}: prepared input/dependency binding mismatch"
            )
    output_records = receipt.get("raw_outputs")
    if not isinstance(output_records, Mapping) or set(output_records) != set(
        RAW_OUTPUT_ROLES
    ):
        raise ValueError(f"{model}: execution receipt raw outputs are incomplete")
    for role in RAW_OUTPUT_ROLES:
        record = output_records[role]
        if not isinstance(record, Mapping):
            raise TypeError(f"{model}/{role}: malformed receipt output record")
        supplied = pathlib.Path(raw_outputs[role]).resolve()
        if pathlib.Path(str(record.get("path", ""))).resolve() != supplied:
            raise ValueError(f"{model}/{role}: receipt output path mismatch")
        require_sha256(supplied, str(record.get("sha256", "")), label=f"{model}/{role}")
        if supplied.stat().st_size != record.get("bytes"):
            raise ValueError(f"{model}/{role}: receipt output size mismatch")
    return receipt


def collect_evidence(
    *,
    prepared_root: pathlib.Path,
    raw_outputs: Mapping[str, Mapping[str, pathlib.Path]],
    execution_receipts: Mapping[str, pathlib.Path],
    out: pathlib.Path,
) -> dict[str, Any]:
    """Copy raw outputs into a hash-manifested bundle and parse them."""
    verify_manifest(prepared_root)
    require_outside_tree(out, prepared_root, label="evidence output")
    preparation_path = prepared_root / "preparation.json"
    preparation = read_json(preparation_path)
    prepared_models = preparation.get("models")
    if not isinstance(prepared_models, Mapping):
        raise TypeError("preparation has no typed model map")
    if set(raw_outputs) != set(prepared_models):
        raise ValueError("raw-output model set does not exactly match preparation")
    if set(execution_receipts) != set(prepared_models):
        raise ValueError(
            "execution-receipt model set does not exactly match preparation"
        )
    prepare_output_dir(out)
    evidence_models: dict[str, Any] = {}
    for model, supplied in raw_outputs.items():
        if pathlib.PurePath(model).name != model:
            raise ValueError(f"unsafe model name in preparation: {model!r}")
        if not isinstance(supplied, Mapping) or set(supplied) != set(RAW_OUTPUT_ROLES):
            raise ValueError(
                f"{model}: raw outputs must provide exactly {RAW_OUTPUT_ROLES}"
            )
        receipt = verify_execution_receipt(
            receipt_path=pathlib.Path(execution_receipts[model]),
            prepared_root=prepared_root,
            model=model,
            raw_outputs=supplied,
        )
        target_dir = out / "raw" / model
        target_dir.mkdir(parents=True)
        records: dict[str, Any] = {}
        copied: dict[str, pathlib.Path] = {}
        for role in RAW_OUTPUT_ROLES:
            source = pathlib.Path(supplied[role])
            if not source.is_file():
                raise ValueError(f"{model}/{role}: raw output is not a file: {source}")
            target = target_dir / f"{role}.log"
            shutil.copyfile(source, target)
            copied[role] = target
            records[role] = _file_record(target, relative_to=out)
        receipt_target = target_dir / "execution-receipt.json"
        shutil.copyfile(execution_receipts[model], receipt_target)
        receipt_record = _file_record(receipt_target, relative_to=out)
        prepared_model = prepared_models[model]
        if not isinstance(prepared_model, Mapping):
            raise TypeError(f"{model}: malformed prepared model")
        expected_atom_count = prepared_model.get("atom_count")
        expected_cell = prepared_model.get("periodic_cell")
        if not isinstance(expected_atom_count, int) or not isinstance(
            expected_cell, Mapping
        ):
            raise TypeError(f"{model}: prepared native identity fields are missing")
        cp2k = parse_cp2k_barrier(
            initial_log=copied["cp2k_initial"],
            endpoint_log=copied["cp2k_endpoint"],
            band_log=copied["cp2k_band"],
            profile_log=copied["cp2k_band_profile"],
            replica_logs=[copied[role] for role in CP2K_REPLICA_OUTPUT_ROLES],
            expected_model=model,
            expected_atom_count=expected_atom_count,
            expected_cell=expected_cell,
        )
        matched = parse_matched_classical_barrier(
            initial_log=copied["lammps_initial"],
            endpoint_log=copied["lammps_endpoint"],
            neb_log=copied["lammps_neb"],
        )
        evidence_models[model] = {
            "atom_identity_sha256": prepared_model.get("atom_identity_sha256"),
            "raw_outputs": records,
            "execution_receipt": receipt_record,
            "execution_receipt_sha256": canonical_sha256(receipt),
            "cp2k": cp2k,
            "matched_classical": matched,
        }
    evidence = {
        "schema": "e3b-parsed-evidence-v1",
        "preparation_sha256": sha256(preparation_path),
        "parser_contract": {
            "cp2k": "endpoint-geo-opt-and-band-v1",
            "matched_classical": "lammps-endpoint-and-neb-v1",
            "expected_replicas": IMAGE_COUNT,
            "hartree_to_kcal_mol": HARTREE_TO_KCAL_MOL,
        },
        "models": evidence_models,
    }
    write_json(out / "evidence.json", evidence)
    write_manifest(out)
    return evidence


def _prepared_dependency_paths(
    prepared_root: pathlib.Path, input_path: pathlib.Path
) -> list[pathlib.Path]:
    paths = {input_path.resolve()}
    for line in input_path.read_text(encoding="utf-8").splitlines():
        match = re.match(r"\s*[A-Z0-9_]*FILE_NAME\s+(\S+)", line, re.IGNORECASE)
        if not match:
            continue
        value = match.group(1).strip("'\"")
        candidate = (input_path.parent / value).resolve()
        try:
            candidate.relative_to(prepared_root.resolve())
        except ValueError as exc:
            raise ValueError(
                f"smoke input references file outside prepared root: {value}"
            ) from exc
        if not candidate.is_file():
            raise ValueError(f"smoke dependency is not a file: {value}")
        paths.add(candidate)
    return sorted(paths)


def _validated_controller_evidence(
    value: Mapping[str, Any] | None,
) -> dict[str, Any] | None:
    if value is None:
        return None
    if value.get("schema") != "e3b-controller-readback-v1":
        raise ValueError("unexpected controller evidence schema")
    if value.get("controller") not in {"systemd", "docker"}:
        raise ValueError("controller evidence must be systemd or docker")
    readback = value.get("readback")
    if (
        value.get("verified") is not True
        or not isinstance(readback, Mapping)
        or value.get("readback_sha256") != canonical_sha256(readback)
    ):
        raise ValueError("controller evidence lacks verified typed readback")
    return dict(value)


def run_smoke(
    *,
    executable: str,
    prepared_root: pathlib.Path,
    model: str,
    runtime_dir: pathlib.Path,
    receipt_path: pathlib.Path,
    timeout_seconds: float,
    mpi_ranks: int = 1,
    omp_threads: int = 1,
    memory_limit_mib: int | None = None,
    execution_method: str = "direct",
    executable_identity: str = "unverified",
    controller_evidence: Mapping[str, Any] | None = None,
    timeout_ready_file: pathlib.Path | None = None,
) -> dict[str, Any]:
    """Run a planning-only smoke probe from an isolated disposable copy."""
    if timeout_seconds <= 0 or timeout_seconds > 14400:
        raise ValueError("smoke timeout must be in (0, 14400] seconds")
    if mpi_ranks <= 0 or omp_threads <= 0 or mpi_ranks * omp_threads > 16:
        raise ValueError("smoke CPU resources must be positive and total at most 16")
    if memory_limit_mib is not None and memory_limit_mib <= 0:
        raise ValueError("smoke memory limit must be positive when declared")
    if not execution_method.strip():
        raise ValueError("smoke execution method must be explicit")
    if timeout_ready_file is not None and executable_identity != "test/fake":
        raise ValueError("timeout readiness synchronization is test-fixture only")
    resolved_timeout_ready_file = (
        timeout_ready_file.resolve() if timeout_ready_file is not None else None
    )
    if resolved_timeout_ready_file is not None and resolved_timeout_ready_file.exists():
        raise ValueError("timeout readiness file must not exist before launch")
    controller = _validated_controller_evidence(controller_evidence)
    verify_manifest(prepared_root)
    preparation_path = prepared_root / "preparation.json"
    preparation = read_json(preparation_path)
    prepared_models = preparation.get("models")
    if not isinstance(prepared_models, Mapping) or model not in prepared_models:
        raise ValueError(f"unknown prepared smoke model: {model}")
    method = preparation.get("method")
    if not isinstance(method, Mapping):
        raise TypeError("preparation has no typed method record")
    executable_record = _resolved_executable_record(
        executable, identity=executable_identity
    )
    resolved_executable = str(executable_record["path"])
    input_path = prepared_root / "models" / model / "smoke.inp"
    dependencies = _prepared_dependency_paths(prepared_root, input_path)
    require_outside_tree(runtime_dir, prepared_root, label="smoke runtime directory")
    prepare_output_dir(runtime_dir)
    if receipt_path.parent.resolve() != runtime_dir.resolve():
        raise ValueError("smoke receipt must live directly in the runtime directory")
    work_root = runtime_dir / "work"
    work_model_dir = work_root / "models" / model
    dependency_records = []
    for source in dependencies:
        relative = source.relative_to(prepared_root.resolve())
        target = work_root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        dependency_records.append(
            {
                "prepared_path": str(relative),
                "prepared_sha256": sha256(source),
                "runtime_sha256": sha256(target),
                "bytes": target.stat().st_size,
            }
        )
    runtime_input = work_model_dir / "smoke.inp"
    stdout_path = runtime_dir / "smoke.stdout.log"
    timed_out = False
    returncode: int | None = None
    output = ""
    process: subprocess.Popen[str] | None = None
    process_group_id: int | None = None
    term_sent = False
    kill_sent = False
    process_reaped = False
    environment = os.environ.copy()
    environment.update(
        {
            "OMP_NUM_THREADS": str(omp_threads),
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    command = [resolved_executable, "-i", runtime_input.name]
    started = time.monotonic()
    try:
        # Trusted caller-selected executable; arguments are fixed and shell=False.
        # nosemgrep: python.lang.security.audit.dangerous-subprocess-use-audit
        process = subprocess.Popen(
            command,
            cwd=work_model_dir,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=False,
            start_new_session=True,
        )
        process_group_id = process.pid
        if resolved_timeout_ready_file is not None:
            readiness_deadline = time.monotonic() + 5.0
            while not resolved_timeout_ready_file.exists():
                if process.poll() is not None:
                    raise RuntimeError(
                        "test fixture exited before advertising readiness"
                    )
                if time.monotonic() >= readiness_deadline:
                    raise TimeoutError("test fixture did not advertise readiness")
                time.sleep(0.001)
        try:
            output, _unused = process.communicate(timeout=timeout_seconds)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(process_group_id, signal.SIGTERM)
            term_sent = True
            try:
                output, _unused = process.communicate(timeout=0.5)
            except subprocess.TimeoutExpired:
                os.killpg(process_group_id, signal.SIGKILL)
                kill_sent = True
                output, _unused = process.communicate(timeout=1.0)
        returncode = process.returncode
        process_reaped = process.poll() is not None
    finally:
        if process is not None and process.poll() is None:
            if process_group_id is not None:
                os.killpg(process_group_id, signal.SIGKILL)
                kill_sent = True
            process.wait(timeout=1.0)
            process_reaped = True
    elapsed = time.monotonic() - started
    stdout_path.write_text(output or "", encoding="utf-8")
    copied_runtime_paths = {
        str((work_root / pathlib.Path(record["prepared_path"])).resolve())
        for record in dependency_records
    }
    generated_outputs = [
        _file_record(path, relative_to=work_root)
        for path in sorted(item for item in work_root.rglob("*") if item.is_file())
        if str(path.resolve()) not in copied_runtime_paths
    ]
    shutil.rmtree(work_root)
    prepared_valid = True
    try:
        verify_manifest(prepared_root)
    except (OSError, TypeError, ValueError):
        prepared_valid = False
    receipt = {
        **classify_cp2k_output(
            output,
            returncode=returncode,
            timed_out=timed_out,
        ),
        "schema": "e3b-cp2k-smoke-v3",
        "purpose": "planning-only ENERGY_FORCE timing evidence; not a BAND bound",
        "preparation": {
            "path": str(preparation_path.resolve()),
            "sha256": sha256(preparation_path),
            "manifest_sha256": sha256(prepared_root / "manifest.json"),
        },
        "model": model,
        "atom_identity_sha256": prepared_models[model].get("atom_identity_sha256"),
        "method": {
            "prepared_declaration": dict(method),
            "prepared_record_sha256": canonical_sha256(method),
            "resolved_executable_before_launch": executable_record,
        },
        "command": command,
        "execution_declaration": {
            "value": execution_method,
            "verification": "declared-not-enforced-by-generic-runner",
        },
        "external_controller_evidence": controller,
        "working_directory_policy": "isolated copied tree outside prepared manifest",
        "resource_declarations": {
            "mpi_ranks": {
                "value": mpi_ranks,
                "verification": "declared-not-enforced-by-generic-runner",
            },
            "memory_limit_mib": {
                "value": memory_limit_mib,
                "verification": "declared-not-enforced-by-generic-runner",
            },
        },
        "environment_controls": {
            "omp_threads": {"value": omp_threads, "mechanism": "OMP_NUM_THREADS"},
            "openblas_threads": {"value": 1, "mechanism": "OPENBLAS_NUM_THREADS"},
            "mkl_threads": {"value": 1, "mechanism": "MKL_NUM_THREADS"},
            "numexpr_threads": {"value": 1, "mechanism": "NUMEXPR_NUM_THREADS"},
        },
        "timeout_seconds": timeout_seconds,
        "elapsed_seconds": elapsed,
        "input": _file_record(input_path, relative_to=prepared_root),
        "dependencies": dependency_records,
        "output": {
            "stdout": _file_record(stdout_path, relative_to=runtime_dir),
            "generated_before_cleanup": generated_outputs,
        },
        "cleanup": {
            "policy": "remove disposable work tree after hashing generated files",
            "work_directory_removed": not work_root.exists(),
            "internal_process_group": {
                "new_session": True,
                "process_group_id": process_group_id,
                "timeout_term_sent": term_sent,
                "timeout_kill_sent": kill_sent,
                "subprocess_reaped": process_reaped,
            },
        },
        "post_run_state": {
            "prepared_manifest_valid": prepared_valid,
            "prepared_tree_modified": not prepared_valid,
        },
    }
    if not prepared_valid:
        receipt["status"] = "incomplete-prepared-integrity"
        receipt["energy_hartree"] = None
    write_json(receipt_path, receipt)
    return receipt


def _failure_result(outcome: str, reason: str) -> dict[str, Any]:
    return {
        "typed_outcome": outcome,
        "reason": reason,
        "matched_cell_transfer_pass": False,
        "matched_cell_transfer_delta_kcal_mol": None,
        "matched_cell_correction_kcal_mol": None,
        "calibrated_barrier_kcal_mol": None,
        "endorsement_delta_kcal_mol": None,
        "verdict": "no-correction",
    }


def _validated_source_barrier(
    reference: Mapping[str, Any], name: str, source: Mapping[str, Any]
) -> float:
    """Re-open, hash, and independently parse the pinned E3a source bundle."""
    required = {
        "workbook_path",
        "workbook_sha256",
        "campaign_root",
        "campaign_manifest_sha256",
        "campaign_result_sha256",
    }
    if not required <= set(source):
        raise ValueError(f"{name}: preparation source paths/hashes are incomplete")
    workbook = pathlib.Path(str(source["workbook_path"]))
    campaign_root = pathlib.Path(str(source["campaign_root"]))
    require_sha256(
        workbook, str(source["workbook_sha256"]), label="E3a source workbook"
    )
    manifest_path = campaign_root / "manifest.json"
    require_sha256(
        manifest_path,
        str(source["campaign_manifest_sha256"]),
        label="E3a campaign manifest",
    )
    result_path = campaign_root / "campaign-result.json"
    require_sha256(
        result_path,
        str(source["campaign_result_sha256"]),
        label="E3a campaign result",
    )
    manifest = read_json(manifest_path)
    campaign = read_json(result_path)
    if campaign.get("schema") != "e3a-classical-neb-campaign-v1":
        raise ValueError(f"{name}: unexpected E3a campaign schema")
    campaign_source = campaign.get("source")
    if not isinstance(campaign_source, Mapping) or campaign_source.get(
        "source_sha256"
    ) != source.get("workbook_sha256"):
        raise ValueError(f"{name}: E3a campaign workbook binding mismatch")
    files = manifest.get("files")
    if not isinstance(files, list):
        raise TypeError(f"{name}: E3a campaign manifest has no typed files")
    result_entries = [
        item
        for item in files
        if isinstance(item, Mapping) and item.get("path") == "campaign-result.json"
    ]
    if len(result_entries) != 1 or result_entries[0].get("sha256") != source.get(
        "campaign_result_sha256"
    ):
        raise ValueError(f"{name}: E3a manifest campaign-result binding mismatch")
    fresh = _classical_reference(campaign, name, campaign_root, manifest)
    if canonical_sha256(fresh) != canonical_sha256(reference):
        raise ValueError(
            f"{name}: stored source reference does not match reparsed E3a files"
        )
    if fresh.get("status") != "converged":
        raise ValueError(f"{name}: source barrier is not independently converged")
    barrier = fresh.get("barrier_kcal_mol")
    if not isinstance(barrier, (int, float)) or not math.isfinite(float(barrier)):
        raise ValueError(f"{name}: reparsed source barrier is not finite")
    return float(barrier)


def _verified_evidence_paths(
    evidence_root: pathlib.Path, records: Mapping[str, Any]
) -> dict[str, pathlib.Path]:
    if set(records) != set(RAW_OUTPUT_ROLES):
        raise ValueError("raw evidence role set is incomplete")
    resolved: dict[str, pathlib.Path] = {}
    root = evidence_root.resolve()
    for role in RAW_OUTPUT_ROLES:
        record = records[role]
        if not isinstance(record, Mapping):
            raise ValueError(f"{role}: malformed raw evidence record")
        relative = record.get("path")
        expected_hash = record.get("sha256")
        expected_size = record.get("bytes")
        if not isinstance(relative, str) or not isinstance(expected_hash, str):
            raise ValueError(f"{role}: untyped raw evidence record")
        path = (root / relative).resolve()
        try:
            path.relative_to(root)
        except ValueError as exc:
            raise ValueError(f"{role}: unsafe raw evidence path") from exc
        require_sha256(path, expected_hash, label=f"raw evidence {role}")
        if path.stat().st_size != expected_size:
            raise ValueError(f"{role}: raw evidence size mismatch")
        resolved[role] = path
    return resolved


def _verify_bundled_execution_receipt(
    *,
    evidence_root: pathlib.Path,
    observed: Mapping[str, Any],
    prepared_root: pathlib.Path,
    preparation: Mapping[str, Any],
    model: str,
    raw_records: Mapping[str, Any],
) -> None:
    record = observed.get("execution_receipt")
    if not isinstance(record, Mapping):
        raise ValueError(f"{model}: bundled execution receipt is missing")
    relative = record.get("path")
    expected_hash = record.get("sha256")
    if not isinstance(relative, str) or not isinstance(expected_hash, str):
        raise TypeError(f"{model}: bundled execution receipt record is malformed")
    root = evidence_root.resolve()
    path = (root / relative).resolve()
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise ValueError(f"{model}: unsafe bundled execution receipt path") from exc
    require_sha256(path, expected_hash, label=f"{model} bundled execution receipt")
    if path.stat().st_size != record.get("bytes"):
        raise ValueError(f"{model}: bundled execution receipt size mismatch")
    receipt = read_json(path)
    if canonical_sha256(receipt) != observed.get("execution_receipt_sha256"):
        raise ValueError(f"{model}: bundled execution receipt content mismatch")
    models = preparation.get("models")
    method = preparation.get("method")
    if not isinstance(models, Mapping) or not isinstance(method, Mapping):
        raise TypeError("preparation model/method record is malformed")
    prepared_model = models.get(model)
    if not isinstance(prepared_model, Mapping):
        raise ValueError(f"{model}: prepared model is missing")
    binding = receipt.get("preparation")
    if not isinstance(binding, Mapping) or binding.get("sha256") != sha256(
        prepared_root / "preparation.json"
    ):
        raise ValueError(f"{model}: bundled receipt preparation mismatch")
    cell = prepared_model.get("periodic_cell")
    receipt_method = receipt.get("method")
    if (
        receipt.get("model") != model
        or receipt.get("atom_identity_sha256")
        != prepared_model.get("atom_identity_sha256")
        or not isinstance(cell, Mapping)
        or receipt.get("periodic_cell_record_sha256") != canonical_sha256(cell)
        or not isinstance(receipt_method, Mapping)
        or receipt_method.get("record_sha256") != canonical_sha256(method)
    ):
        raise ValueError(f"{model}: bundled receipt model/method mismatch")
    expected_unit_paths = _prepared_unit_paths(prepared_root, model)
    units = receipt.get("units")
    if not isinstance(units, Mapping) or set(units) != set(SOLVER_UNITS):
        raise ValueError(f"{model}: bundled receipt units are incomplete")
    for unit, expected_paths in expected_unit_paths.items():
        unit_record = units[unit]
        if not isinstance(unit_record, Mapping):
            raise TypeError(f"{model}/{unit}: malformed bundled receipt unit")
        if unit_record.get("prepared_inputs") != [
            _file_record(item, relative_to=prepared_root) for item in expected_paths
        ]:
            raise ValueError(f"{model}/{unit}: bundled prepared inputs changed")
        if not isinstance(unit_record.get("command"), list) or not isinstance(
            unit_record.get("resources"), Mapping
        ):
            raise ValueError(f"{model}/{unit}: command/resources binding is missing")
    receipt_outputs = receipt.get("raw_outputs")
    if not isinstance(receipt_outputs, Mapping) or set(receipt_outputs) != set(
        RAW_OUTPUT_ROLES
    ):
        raise ValueError(f"{model}: bundled receipt raw outputs are incomplete")
    for role in RAW_OUTPUT_ROLES:
        receipt_output = receipt_outputs[role]
        evidence_output = raw_records[role]
        if not isinstance(receipt_output, Mapping) or not isinstance(
            evidence_output, Mapping
        ):
            raise TypeError(f"{model}/{role}: malformed output binding")
        if receipt_output.get("sha256") != evidence_output.get(
            "sha256"
        ) or receipt_output.get("bytes") != evidence_output.get("bytes"):
            raise ValueError(f"{model}/{role}: receipt/evidence output mismatch")


def analyze_models(
    preparation: Mapping[str, Any],
    observations: Mapping[str, Any],
    *,
    prepared_root: pathlib.Path | None = None,
    evidence_root: pathlib.Path | None = None,
    preparation_sha256: str | None = None,
    transfer_tolerance_kcal_mol: float = DEFAULT_TRANSFER_TOLERANCE_KCAL_MOL,
    endorsement_tolerance_kcal_mol: float = DEFAULT_ENDORSEMENT_TOLERANCE_KCAL_MOL,
) -> dict[str, Any]:
    if preparation.get("schema") != "e3b-cp2k-preparation-v1":
        raise ValueError("unexpected E3b preparation schema")
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
        raise TypeError("preparation and evidence require typed model maps")
    if set(prepared_models) != set(observed_models):
        raise ValueError("evidence model set does not exactly match preparation")
    results: dict[str, Any] = {}
    is_parsed_evidence = observations.get("schema") == "e3b-parsed-evidence-v1"
    if not is_parsed_evidence and observations.get("schema") != (
        "e3b-cp2k-observations-v1"
    ):
        raise ValueError("unexpected E3b evidence schema")
    if is_parsed_evidence:
        if prepared_root is None or evidence_root is None or preparation_sha256 is None:
            raise ValueError(
                "parsed evidence analysis requires prepared/evidence roots and preparation hash"
            )
        if observations.get("preparation_sha256") != preparation_sha256:
            raise ValueError("parsed evidence is not bound to this preparation")

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
                "evidence atom identity does not match the prepared path",
            )
            continue
        reference = prepared_value.get("classical_reference")
        if not isinstance(reference, Mapping):
            raise TypeError(f"{name}: missing classical reference")
        if reference.get("status") != "converged":
            results[name] = _failure_result(
                "incomplete-source-classical",
                "E3a source classical barrier is not independently converged",
            )
            continue

        if not is_parsed_evidence:
            cp2k = observed_value.get("cp2k")
            matched = observed_value.get("matched_classical")
            if not isinstance(cp2k, Mapping) or not isinstance(matched, Mapping):
                raise TypeError(f"{name}: missing typed CP2K or classical observation")
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
            elif matched_status != "converged":
                results[name] = _failure_result(
                    "incomplete-matched-classical",
                    "matched 2x2 classical CI-NEB is not typed converged",
                )
            else:
                results[name] = _failure_result(
                    "incomplete-unverified-observation",
                    "manual observations are not accepted as barrier evidence; "
                    "hash-bound raw-output parsers are required",
                )
            continue

        try:
            source = preparation.get("source")
            if not isinstance(source, Mapping):
                raise ValueError("preparation source binding is missing")
            source_barrier = _validated_source_barrier(reference, name, source)
            raw_records = observed_value.get("raw_outputs")
            if not isinstance(raw_records, Mapping):
                raise ValueError("raw evidence records are missing")
            assert evidence_root is not None
            assert prepared_root is not None
            paths = _verified_evidence_paths(evidence_root, raw_records)
            _verify_bundled_execution_receipt(
                evidence_root=evidence_root,
                observed=observed_value,
                prepared_root=prepared_root,
                preparation=preparation,
                model=name,
                raw_records=raw_records,
            )
            expected_atom_count = prepared_value.get("atom_count")
            expected_cell = prepared_value.get("periodic_cell")
            if not isinstance(expected_atom_count, int) or not isinstance(
                expected_cell, Mapping
            ):
                raise ValueError(f"{name}: prepared native identity is missing")
            reparsed_cp2k = parse_cp2k_barrier(
                initial_log=paths["cp2k_initial"],
                endpoint_log=paths["cp2k_endpoint"],
                band_log=paths["cp2k_band"],
                profile_log=paths["cp2k_band_profile"],
                replica_logs=[paths[role] for role in CP2K_REPLICA_OUTPUT_ROLES],
                expected_model=name,
                expected_atom_count=expected_atom_count,
                expected_cell=expected_cell,
            )
            reparsed_matched = parse_matched_classical_barrier(
                initial_log=paths["lammps_initial"],
                endpoint_log=paths["lammps_endpoint"],
                neb_log=paths["lammps_neb"],
            )
            if canonical_sha256(observed_value.get("cp2k")) != canonical_sha256(
                reparsed_cp2k
            ) or canonical_sha256(
                observed_value.get("matched_classical")
            ) != canonical_sha256(reparsed_matched):
                raise ValueError(
                    "stored parser fields do not match reparsed raw outputs"
                )
        except (OSError, TypeError, ValueError) as exc:
            results[name] = _failure_result("incomplete-evidence-integrity", str(exc))
            continue

        cp2k_barrier = float(reparsed_cp2k["barrier_kcal_mol"])
        matched_barrier = float(reparsed_matched["barrier_kcal_mol"])
        transfer_delta = abs(matched_barrier - source_barrier)
        transfer_pass = transfer_delta <= transfer_tolerance_kcal_mol
        if not transfer_pass:
            result = _failure_result(
                "complete-transfer-rejected",
                "matched-cell classical barrier does not transfer from the E3a cell",
            )
            result["matched_cell_transfer_delta_kcal_mol"] = transfer_delta
            results[name] = result
            continue
        correction = cp2k_barrier - matched_barrier
        calibrated = source_barrier + correction
        endorsement_delta = abs(correction)
        endorsed = endorsement_delta <= endorsement_tolerance_kcal_mol
        results[name] = {
            "typed_outcome": "complete-endorsed" if endorsed else "complete-corrected",
            "reason": (
                "hash-bound raw outputs passed all convergence and transfer gates"
            ),
            "source_classical_barrier_kcal_mol": source_barrier,
            "matched_classical_barrier_kcal_mol": matched_barrier,
            "cp2k_barrier_kcal_mol": cp2k_barrier,
            "matched_cell_transfer_pass": True,
            "matched_cell_transfer_delta_kcal_mol": transfer_delta,
            "matched_cell_correction_kcal_mol": correction,
            "calibrated_barrier_kcal_mol": calibrated,
            "endorsement_delta_kcal_mol": endorsement_delta,
            "verdict": "endorsed" if endorsed else "corrected",
            "raw_output_sha256": {
                role: raw_records[role]["sha256"] for role in RAW_OUTPUT_ROLES
            },
        }
    return {
        "schema": "e3b-cp2k-analysis-v1",
        "evidence_contract": (
            "hash-bound raw outputs reparsed during analysis"
            if is_parsed_evidence
            else "manual observations rejected"
        ),
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
    controller = (
        read_json(args.controller_evidence)
        if args.controller_evidence is not None
        else None
    )
    receipt = run_smoke(
        executable=args.cp2k,
        prepared_root=args.prepared_root,
        model=args.model,
        runtime_dir=args.runtime_dir,
        receipt_path=args.runtime_dir / "smoke-result.json",
        timeout_seconds=args.timeout_seconds,
        mpi_ranks=args.mpi_ranks,
        omp_threads=args.omp_threads,
        memory_limit_mib=args.memory_limit_mib,
        execution_method=args.execution_method,
        executable_identity=args.executable_identity,
        controller_evidence=controller,
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0 if receipt["status"] == "converged" else 2


def _spec_path(spec_path: pathlib.Path, value: str) -> pathlib.Path:
    path = pathlib.Path(value)
    return path if path.is_absolute() else spec_path.parent / path


def command_create_receipt(args: argparse.Namespace) -> int:
    raw_spec = read_json(args.raw_spec)
    execution_spec = read_json(args.execution_spec)
    models = raw_spec.get("models")
    if raw_spec.get("schema") != "e3b-raw-output-spec-v1" or not isinstance(
        models, Mapping
    ):
        raise ValueError("unexpected raw-output specification")
    raw_model = models.get(args.model)
    if not isinstance(raw_model, Mapping):
        raise ValueError(f"raw-output specification lacks model {args.model}")
    outputs = raw_model.get("outputs")
    if not isinstance(outputs, Mapping):
        raise TypeError(f"{args.model}: raw-output specification lacks outputs")
    if execution_spec.get("schema") != "e3b-execution-declarations-v1":
        raise ValueError("unexpected execution declaration schema")
    executables = execution_spec.get("executables")
    units = execution_spec.get("units")
    if not isinstance(executables, Mapping) or not isinstance(units, Mapping):
        raise TypeError("execution declarations lack executables/units")
    executable_paths: dict[str, pathlib.Path] = {}
    executable_identities: dict[str, str] = {}
    for solver in ("cp2k", "lammps"):
        record = executables.get(solver)
        if (
            not isinstance(record, Mapping)
            or not isinstance(record.get("path"), str)
            or not isinstance(record.get("identity"), str)
        ):
            raise TypeError(f"execution declaration lacks {solver} identity")
        executable_paths[solver] = _spec_path(args.execution_spec, str(record["path"]))
        executable_identities[solver] = str(record["identity"])
    commands: dict[str, Sequence[str]] = {}
    resources: dict[str, Mapping[str, Any]] = {}
    for unit in SOLVER_UNITS:
        record = units.get(unit)
        if not isinstance(record, Mapping):
            raise TypeError(f"execution declaration lacks unit {unit}")
        command = record.get("command")
        unit_resources = record.get("resources")
        if not isinstance(command, list) or not isinstance(unit_resources, Mapping):
            raise TypeError(f"execution declaration unit {unit} is malformed")
        commands[unit] = command
        resources[unit] = unit_resources
    raw_outputs = {
        role: _spec_path(args.raw_spec, str(outputs[role])) for role in RAW_OUTPUT_ROLES
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    receipt = create_execution_receipt(
        prepared_root=args.prepared_root,
        model=args.model,
        raw_outputs=raw_outputs,
        executables=executable_paths,
        executable_identities=executable_identities,
        commands=commands,
        resources=resources,
        out=args.out,
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


def command_collect(args: argparse.Namespace) -> int:
    raw_spec = read_json(args.raw_spec)
    if raw_spec.get("schema") != "e3b-raw-output-spec-v1":
        raise ValueError("unexpected raw-output specification schema")
    models = raw_spec.get("models")
    if not isinstance(models, Mapping):
        raise TypeError("raw-output specification has no typed model map")
    raw_outputs: dict[str, dict[str, pathlib.Path]] = {}
    execution_receipts: dict[str, pathlib.Path] = {}
    for model, value in models.items():
        if not isinstance(model, str) or not isinstance(value, Mapping):
            raise TypeError("raw-output specification model record is malformed")
        outputs = value.get("outputs")
        receipt_value = value.get("execution_receipt")
        if not isinstance(outputs, Mapping) or not isinstance(receipt_value, str):
            raise TypeError(f"{model}: outputs/receipt specification is malformed")
        raw_outputs[model] = {}
        for role, supplied in outputs.items():
            if not isinstance(role, str) or not isinstance(supplied, str):
                raise TypeError(f"{model}: raw-output paths must be strings")
            raw_outputs[model][role] = _spec_path(args.raw_spec, supplied)
        execution_receipts[model] = _spec_path(args.raw_spec, receipt_value)
    evidence = collect_evidence(
        prepared_root=args.prepared_root,
        raw_outputs=raw_outputs,
        execution_receipts=execution_receipts,
        out=args.out,
    )
    print(json.dumps(evidence, indent=2, sort_keys=True))
    return 0


def command_analyze(args: argparse.Namespace) -> int:
    verify_manifest(args.prepared_root)
    preparation_path = args.prepared_root / "preparation.json"
    preparation = read_json(preparation_path)
    evidence_root = getattr(args, "evidence_root", None)
    observations_path = getattr(args, "observations", None)
    if evidence_root is not None:
        verify_manifest(evidence_root)
        evidence_path = evidence_root / "evidence.json"
        observations = read_json(evidence_path)
    elif observations_path is not None:
        evidence_path = observations_path
        observations = read_json(evidence_path)
        if observations.get("preparation_sha256") != sha256(preparation_path):
            raise ValueError("observations are not hash-bound to this preparation")
    else:
        raise ValueError("analysis requires a parsed evidence bundle")
    require_outside_tree(args.out, args.prepared_root, label="analysis output")
    prepare_output_dir(args.out)
    result = analyze_models(
        preparation,
        observations,
        prepared_root=args.prepared_root,
        evidence_root=evidence_root,
        preparation_sha256=sha256(preparation_path),
        transfer_tolerance_kcal_mol=args.transfer_tolerance,
        endorsement_tolerance_kcal_mol=args.endorsement_tolerance,
    )
    result["preparation_sha256"] = sha256(preparation_path)
    result["evidence_sha256"] = sha256(evidence_path)
    write_json(args.out / "analysis.json", result)
    write_manifest(args.out)
    print(json.dumps(result, indent=2, sort_keys=True))
    outcomes = [model["typed_outcome"] for model in result["models"].values()]
    return (
        0 if outcomes and all(item.startswith("complete-") for item in outcomes) else 2
    )


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
    smoke.add_argument(
        "--runtime-dir",
        type=pathlib.Path,
        required=True,
        help="new directory outside the prepared tree for outputs and receipt",
    )
    smoke.add_argument("--timeout-seconds", type=float, default=600.0)
    smoke.add_argument("--mpi-ranks", type=int, default=1)
    smoke.add_argument("--omp-threads", type=int, default=1)
    smoke.add_argument("--memory-limit-mib", type=int)
    smoke.add_argument(
        "--execution-method",
        default="direct",
        help="unverified execution declaration recorded in the timing receipt",
    )
    smoke.add_argument(
        "--executable-identity",
        choices=("cp2k", "wrapper", "test/fake", "unverified"),
        default="unverified",
    )
    smoke.add_argument(
        "--controller-evidence",
        type=pathlib.Path,
        help="typed systemd/Docker readback produced by an external controller",
    )
    smoke.set_defaults(func=command_smoke)

    receipt = sub.add_parser(
        "create-execution-receipt",
        help="bind prepared inputs, commands/resources, executables, and raw outputs",
    )
    receipt.add_argument("--prepared-root", type=pathlib.Path, required=True)
    receipt.add_argument("--model", required=True)
    receipt.add_argument("--raw-spec", type=pathlib.Path, required=True)
    receipt.add_argument("--execution-spec", type=pathlib.Path, required=True)
    receipt.add_argument("--out", type=pathlib.Path, required=True)
    receipt.set_defaults(func=command_create_receipt)

    collect = sub.add_parser(
        "collect-evidence", help="parse and hash-bind CP2K/LAMMPS raw outputs"
    )
    collect.add_argument("--prepared-root", type=pathlib.Path, required=True)
    collect.add_argument("--raw-spec", type=pathlib.Path, required=True)
    collect.add_argument("--out", type=pathlib.Path, required=True)
    collect.set_defaults(func=command_collect)

    analyze = sub.add_parser(
        "analyze", help="apply only converged, matched-cell DFT corrections"
    )
    analyze.add_argument("--prepared-root", type=pathlib.Path, required=True)
    analyze.add_argument("--evidence-root", type=pathlib.Path, required=True)
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
