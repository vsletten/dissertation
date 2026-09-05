#!/usr/bin/env python3
"""Build the typed atomistic input families required by E3a.

The generator starts from the hash-pinned Nteme et al. (2022) supplementary
workbook.  It does not compute barriers.  It makes every published input,
derived transformation, and sensitivity assumption explicit in a JSON
manifest and emits LAMMPS run-0 inputs for cheap structural validation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import pathlib
import subprocess
import sys
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import asdict, dataclass, replace

import nteme_neb

XENON_MASS = 131.293
# Bourg & Sposito (2008), as used with ClayFF for noble gases in clay
# interlayers by Gadikota et al.: neutral Xe 12-6 LJ.  Epsilon is
# converted from 0.0185900 eV with 1 eV = 23.060547830619 kcal/mol.
XENON_SIGMA = 4.00924
XENON_EPSILON = 0.4286955841712072
MIN_DISTANCE_ANGSTROM = 0.70
STATIC_FTOL = 1.0e-6


@dataclass(frozen=True)
class Model:
    name: str
    atoms: tuple[nteme_neb.Atom, ...]
    cell: nteme_neb.Cell
    metadata: dict[str, object]
    endpoint_atoms: tuple[nteme_neb.Atom, ...] | None = None


def _vector(
    cell: nteme_neb.Cell, fractional: Sequence[float]
) -> tuple[float, float, float]:
    lx, ly, lz, xy, xz, yz = cell.restricted()
    sx, sy, sz = fractional
    return lx * sx + xy * sy + xz * sz, ly * sy + yz * sz, lz * sz


def fractional(
    atom: nteme_neb.Atom, cell: nteme_neb.Cell
) -> tuple[float, float, float]:
    inv = nteme_neb.matrix_inverse(cell)
    return (
        inv[0][0] * atom.x + inv[0][1] * atom.y + inv[0][2] * atom.z,
        inv[1][1] * atom.y + inv[1][2] * atom.z,
        inv[2][2] * atom.z,
    )


def periodic_distance(
    a: nteme_neb.Atom, b: nteme_neb.Atom, cell: nteme_neb.Cell
) -> float:
    return math.sqrt(sum(v * v for v in nteme_neb.minimum_image_delta(a, b, cell)))


def periodic_midpoint(
    a: nteme_neb.Atom, b: nteme_neb.Atom, cell: nteme_neb.Cell
) -> nteme_neb.Atom:
    dx, dy, dz = nteme_neb.minimum_image_delta(a, b, cell)
    return nteme_neb.Atom(0, "X", a.x + dx / 2, a.y + dy / 2, a.z + dz / 2)


def atom_charge(element: str) -> float:
    return 0.0 if element == "Xe" else nteme_neb.CHARGE[element]


def net_charge(atoms: Iterable[nteme_neb.Atom]) -> float:
    return sum(atom_charge(atom.element) for atom in atoms)


def nearest(
    atoms: Iterable[nteme_neb.Atom], center: nteme_neb.Atom, cell: nteme_neb.Cell
):
    return sorted(
        atoms, key=lambda atom: (periodic_distance(atom, center, cell), atom.id)
    )


def tetrahedral_oxygen_ids(
    atom: nteme_neb.Atom,
    atoms: Sequence[nteme_neb.Atom],
    cell: nteme_neb.Cell,
) -> tuple[int, ...]:
    oxygen = [candidate for candidate in atoms if candidate.element in {"O2", "O3"}]
    neighbors = [
        candidate
        for candidate in nearest(oxygen, atom, cell)
        if periodic_distance(atom, candidate, cell) < 2.25
    ]
    if len(neighbors) != 4:
        raise ValueError(
            f"tetrahedral site {atom.id} has {len(neighbors)} O2/O3 neighbors, expected 4"
        )
    return tuple(candidate.id for candidate in neighbors)


def compensate_k_vacancies(
    atoms: Sequence[nteme_neb.Atom],
    cell: nteme_neb.Cell,
    route_center: nteme_neb.Atom,
    vacancy_count: int,
) -> tuple[list[nteme_neb.Atom], list[dict[str, object]]]:
    """Neutralize K vacancies using remote Al2->Si + O3->O2 substitutions.

    With the published ClayFF charges, Al2->Si adds +0.525 e and retyping its
    four Al-adjacent O3 atoms to O2 adds 4*(+0.11875 e), exactly +1 e.
    This is the charge-algebraically consistent interpretation of Nteme's
    otherwise site-ambiguous remote tetrahedral substitutions.
    """

    candidates = sorted(
        (atom for atom in atoms if atom.element == "Al2"),
        key=lambda atom: (-periodic_distance(atom, route_center, cell), atom.id),
    )
    selected: list[nteme_neb.Atom] = []
    selected_oxygen: set[int] = set()
    changes: list[dict[str, object]] = []
    for candidate in candidates:
        if periodic_distance(candidate, route_center, cell) < 8.0:
            continue
        if any(periodic_distance(candidate, prior, cell) < 4.5 for prior in selected):
            continue
        oxygen_ids = tetrahedral_oxygen_ids(candidate, atoms, cell)
        if any(oxygen_id in selected_oxygen for oxygen_id in oxygen_ids):
            continue
        selected.append(candidate)
        selected_oxygen.update(oxygen_ids)
        changes.append(
            {
                "al2_site_id": candidate.id,
                "oxygen_site_ids": list(oxygen_ids),
                "distance_from_route_angstrom": periodic_distance(
                    candidate, route_center, cell
                ),
            }
        )
        if len(selected) == vacancy_count:
            break
    if len(selected) != vacancy_count:
        raise ValueError(
            f"could not select {vacancy_count} separated remote Al2 sites; got {len(selected)}"
        )
    selected_ids = {atom.id for atom in selected}
    compensated = []
    for atom in atoms:
        if atom.id in selected_ids:
            compensated.append(replace(atom, element="Si"))
        elif atom.id in selected_oxygen and atom.element == "O3":
            compensated.append(replace(atom, element="O2"))
        else:
            compensated.append(atom)
    if not math.isclose(net_charge(compensated), 0.0, abs_tol=1.0e-9):
        raise ValueError(
            f"charge compensation failed: q={net_charge(compensated):.12g}"
        )
    return compensated, changes


def vacancy_route(
    atoms: Sequence[nteme_neb.Atom],
    cell: nteme_neb.Cell,
    route: nteme_neb.Barrier,
    species: str,
) -> tuple[list[nteme_neb.Atom], nteme_neb.Atom, nteme_neb.Atom, dict[str, object]]:
    by_id = {atom.id: atom for atom in atoms}
    moving = by_id[route.moving_site]
    destination = by_id[route.vacancy_1]
    removed = {route.vacancy_1}
    if route.vacancy_2 is not None:
        removed.add(route.vacancy_2)
    prepared = [
        replace(atom, element=species) if atom.id == route.moving_site else atom
        for atom in atoms
        if atom.id not in removed
    ]
    center = periodic_midpoint(moving, destination, cell)
    vacancy_count = len(removed) + 1  # moving K is also replaced by neutral gas
    prepared, substitutions = compensate_k_vacancies(
        prepared, cell, center, vacancy_count
    )
    dx, dy, dz = nteme_neb.minimum_image_delta(moving, destination, cell)
    final = nteme_neb.Atom(
        moving.id, species, moving.x + dx, moving.y + dy, moving.z + dz
    )
    return (
        prepared,
        center,
        final,
        {
            "moving_site_id": moving.id,
            "vacancy_site_ids": sorted(removed),
            "final_coordinate_angstrom": [final.x, final.y, final.z],
            "remote_charge_compensation": substitutions,
            "charge_compensation_provenance": "reconstructed-charge-compensation: charge-correct inverse of Nteme prose using published ClayFF atom types; not the undisclosed production structure",
        },
    )


def choose_oh_pair(
    atoms: Sequence[nteme_neb.Atom],
    cell: nteme_neb.Cell,
    route_center: nteme_neb.Atom,
) -> tuple[nteme_neb.Atom, nteme_neb.Atom]:
    oxygens = [atom for atom in atoms if atom.element == "O1"]
    pairs: list[tuple[float, float, int, int, nteme_neb.Atom, nteme_neb.Atom]] = []
    for index, first in enumerate(oxygens):
        for second in oxygens[index + 1 :]:
            separation = periodic_distance(first, second, cell)
            if separation > 3.6:
                continue
            midpoint = periodic_midpoint(first, second, cell)
            route_distance = periodic_distance(midpoint, route_center, cell)
            pairs.append(
                (route_distance, separation, first.id, second.id, first, second)
            )
    if not pairs:
        raise ValueError("no hydroxyl O pair found within 3.6 A")
    return min(pairs)[-2:]


def reconstructed_replication_model(
    pristine: Sequence[nteme_neb.Atom],
    cell: nteme_neb.Cell,
    route: nteme_neb.Barrier,
) -> Model:
    """Build the published route with the contract's neutral reconstruction."""

    atoms, _route_center, final, route_meta = vacancy_route(pristine, cell, route, "Ar")
    endpoint = tuple(
        replace(atom, x=final.x, y=final.y, z=final.z)
        if atom.id == route.moving_site
        else atom
        for atom in atoms
    )
    return Model(
        "reconstructed-replication",
        tuple(atoms),
        cell,
        {
            "outcome": "executable-reconstructed-replication",
            "route": route_meta,
            "published_barrier_kcal_mol": route.barrier_kcal_mol,
            "acceptance_tolerance_kcal_mol": 5.0,
            "limitation": (
                "Nteme et al. do not disclose their route-specific compensator "
                "sites. This neutral model tests the source-constrained inverse "
                "operator and is not an exact reconstruction of unpublished inputs."
            ),
            "neb_endpoint": asdict(final),
        },
        endpoint,
    )


def dehydroxylated_model(
    pristine: Sequence[nteme_neb.Atom],
    cell: nteme_neb.Cell,
    route: nteme_neb.Barrier,
) -> Model:
    atoms, route_center, final, route_meta = vacancy_route(pristine, cell, route, "Ar")
    first_o, second_o = choose_oh_pair(atoms, cell, route_center)
    keep_o, remove_o = sorted((first_o, second_o), key=lambda atom: atom.id)
    bonds = nteme_neb.oh_bonds(list(pristine), cell)
    h_by_o = {oxygen_id: hydrogen_id for oxygen_id, hydrogen_id in bonds}
    removed_ids = {remove_o.id, h_by_o[keep_o.id], h_by_o[remove_o.id]}
    transformed = []
    for atom in atoms:
        if atom.id in removed_ids:
            continue
        if atom.id == keep_o.id:
            transformed.append(replace(atom, element="O2"))
        else:
            transformed.append(atom)
    endpoint = tuple(
        replace(atom, x=final.x, y=final.y, z=final.z)
        if atom.id == route.moving_site
        else atom
        for atom in transformed
    )
    return Model(
        "dehydroxylate-lattice",
        tuple(transformed),
        cell,
        {
            "outcome": "executable-sensitivity-structure",
            "route": route_meta,
            "transformation": {
                "reaction": "2 OH(lattice) -> H2O(removed) + O_residual(lattice)",
                "selected_oh_oxygen_ids": [keep_o.id, remove_o.id],
                "selected_hydrogen_ids": [h_by_o[keep_o.id], h_by_o[remove_o.id]],
                "residual_oxygen_site_id": keep_o.id,
                "removed_site_ids": sorted(removed_ids),
                "residual_oxygen_type": "O2",
                "charge_delta_e": net_charge(transformed) - net_charge(atoms),
                "selection": "nearest O1-O1 pair (<=3.6 A) to route midpoint; ties by IDs",
            },
            "relaxation": {
                "stage_1": "minimize fixed cell at 0 K",
                "stage_2": "NPT 0.5 K, 0 kbar with a,b fixed and c/tilts relaxed",
                "acceptance": "finite energy; no pair <0.70 A; residual O coordinated to two octahedral Al; no H2O retained",
            },
            "limitation": "ClayFF can compare post-reaction endpoint energetics but cannot model O-H bond breaking or H2O formation; no dehydroxylation path/barrier may be inferred.",
            "neb_endpoint": asdict(final),
        },
        endpoint,
    )


def extended_zone_model(
    pristine: Sequence[nteme_neb.Atom],
    cell: nteme_neb.Cell,
    route: nteme_neb.Barrier,
    *,
    radius_a: float = 8.0,
    radius_b: float = 5.0,
    opening: float = 1.0,
) -> Model:
    atoms, route_center, final, route_meta = vacancy_route(pristine, cell, route, "Ar")
    center_fractional = fractional(route_center, cell)

    def deform(atom: nteme_neb.Atom) -> tuple[nteme_neb.Atom, bool]:
        sx, sy, sz = fractional(atom, cell)
        dx = sx - center_fractional[0]
        dy = sy - center_fractional[1]
        dz = sz - center_fractional[2]
        dx -= round(dx)
        dy -= round(dy)
        dz -= round(dz)
        cart_xy = _vector(cell, (dx, dy, 0.0))
        radial = (cart_xy[0] / radius_a) ** 2 + (cart_xy[1] / radius_b) ** 2
        if radial >= 1.0 or math.isclose(dz, 0.0, abs_tol=1.0e-12):
            return atom, False
        weight = math.cos(math.pi * math.sqrt(radial) / 2.0) ** 2
        shift = math.copysign(opening * weight / 2.0, dz)
        # Move along the crystallographic c vector, preserving a/b coordinates.
        cx, cy, cz = _vector(cell, (0.0, 0.0, shift / cell.c))
        return replace(atom, x=atom.x + cx, y=atom.y + cy, z=atom.z + cz), True

    transformed: list[nteme_neb.Atom] = []
    moved = 0
    for atom in atoms:
        deformed, did_move = deform(atom)
        transformed.append(deformed)
        moved += int(did_move)
    deformed_final, _ = deform(final)
    endpoint = tuple(
        replace(atom, x=deformed_final.x, y=deformed_final.y, z=deformed_final.z)
        if atom.id == route.moving_site
        else atom
        for atom in transformed
    )
    return Model(
        "extended-defect-zone",
        tuple(transformed),
        cell,
        {
            "outcome": "executable-sensitivity-family-member",
            "route": route_meta,
            "transformation": {
                "shape": "periodic lenticular basal opening with cosine-squared taper",
                "center": [route_center.x, route_center.y, route_center.z],
                "radius_a_angstrom": radius_a,
                "radius_b_angstrom": radius_b,
                "peak_opening_angstrom": opening,
                "moved_atoms": moved,
                "family": {
                    "radius_a_angstrom": [6.0, 8.0, 10.0],
                    "radius_b_angstrom": [3.0, 5.0, 7.0],
                    "peak_opening_angstrom": [0.5, 1.0, 1.5],
                },
            },
            "relaxation": {
                "boundary": "periodic p p p; fixed cell for each seeded family member",
                "protocol": "restrained shell minimization, then release interior; reject if the void heals before NEB",
            },
            "limitation": "Microscopy supports lenticular voids/basal partings but does not publish one atom-resolved geometry; every dimension is a sensitivity parameter, never a measured Nteme structure.",
            "neb_endpoint": asdict(deformed_final),
        },
        endpoint,
    )


def octahedral_trap_model(
    pristine: Sequence[nteme_neb.Atom],
    cell: nteme_neb.Cell,
    route: nteme_neb.Barrier,
) -> Model:
    del route
    by_id = {atom.id: atom for atom in pristine}
    cage_atom_ids = (743, 741, 708, 711, 628, 707)
    cage_atoms = [by_id[atom_id] for atom_id in cage_atom_ids]
    if not all(atom.element.startswith("O") for atom in cage_atoms):
        raise ValueError("vacant-M1 cage IDs must all be oxygen atom types")
    seed_xyz = tuple(
        sum(getattr(atom, axis) for atom in cage_atoms) / len(cage_atoms)
        for axis in ("x", "y", "z")
    )
    center = nteme_neb.Atom(0, "Ar", *seed_xyz)
    oxygen_shell = sorted(
        (periodic_distance(center, atom, cell), atom.id)
        for atom in pristine
        if atom.element.startswith("O")
    )
    if {atom_id for _distance, atom_id in oxygen_shell[:6]} != set(cage_atom_ids):
        raise ValueError("pinned vacant-M1 cage is not the six-nearest-oxygen shell")
    if oxygen_shell[5][0] > 2.30 or oxygen_shell[6][0] < 2.50:
        raise ValueError("pinned vacant-M1 oxygen shell is not geometrically isolated")

    # Szczerba's illite recoil model supplies a 20-degree analog, not a
    # muscovite coordinate.  Select the K endpoint deterministically from the
    # positive-(001)-normal gallery: smallest angular departure from 20 deg,
    # then shortest displacement and atom ID.
    endpoint_candidates = []
    for atom in pristine:
        if atom.element != "K":
            continue
        dx, dy, dz = nteme_neb.minimum_image_delta(center, atom, cell)
        if dz <= 0.0:
            continue
        distance = math.sqrt(dx * dx + dy * dy + dz * dz)
        angle = math.degrees(math.acos(dz / distance))
        endpoint_candidates.append((abs(angle - 20.0), distance, atom.id, angle, atom))
    _angle_error, _distance, _atom_id, endpoint_angle, parent = min(endpoint_candidates)
    atoms = [
        replace(atom, element="Ar", x=seed_xyz[0], y=seed_xyz[1], z=seed_xyz[2])
        if atom.id == parent.id
        else atom
        for atom in pristine
    ]
    atoms, substitutions = compensate_k_vacancies(atoms, cell, center, 1)
    separation = periodic_distance(center, parent, cell)
    return Model(
        "octahedral-trap-escape",
        tuple(atoms),
        cell,
        {
            "outcome": "incomplete-no-minimum",
            "parent_k_site_id": parent.id,
            "recoil_seed": {
                "coordinate_angstrom": list(seed_xyz),
                "geometric_cage_atom_ids": list(cage_atom_ids),
                "oxygen_shell_distances_angstrom": [
                    distance for distance, _atom_id in oxygen_shell[:6]
                ],
                "selection": "arithmetic centroid of the six pinned oxygen atoms around one vacant-M1 cage; geometric seed for the recoil-derived reservoir, not a recoil trajectory",
                "minimum_image_distance_to_exit_angstrom": separation,
            },
            "remote_charge_compensation": substitutions,
            "charge_compensation_provenance": "reconstructed-charge-compensation: charge-correct inverse of Nteme prose using published ClayFF atom types; not the undisclosed production structure",
            "minimum_search": {
                "seed_family": "centroid plus Cartesian offsets {0,+/-0.25 A} on each axis",
                "relaxation": "freeze atoms >4.25 A from Ar; minimize to 0.05 then 0.01 kcal mol^-1 A^-1; release and minimize entire fixed cell",
                "accept": "finite energy; max force <=0.01 kcal mol^-1 A^-1; Hessian index zero; Ar remains in the same octahedral cage; all distances >=0.70 A",
                "typed_failure": "incomplete-no-minimum",
            },
            "exit_endpoint": {
                "kind": "original K675 gallery site",
                "site_id": parent.id,
                "coordinate_angstrom": [parent.x, parent.y, parent.z],
                "selection": "positive-(001)-normal K site whose cage displacement is closest to the 20-degree illite recoil analog; ties by distance then atom ID",
                "angle_from_layer_normal_degrees": endpoint_angle,
                "construction": "relaxed octahedral minimum -> original parent vacancy; CI-NEB only after the trap passes the minimum gate",
            },
            "narrower_alternative": "periodic-DFT endpoint-only optimization of this exact atom-ID transformation under E3b; the 55 kcal/mol illite value may be an analog only",
            "sources": [
                "Szczerba et al. 2015, DOI 10.1016/j.gca.2015.03.005",
                "Sletten and Onstott 1998, DOI 10.1016/S0016-7037(97)00323-2",
                "Nteme et al. 2022, DOI 10.1016/j.gca.2022.05.004",
            ],
        },
        tuple(
            replace(atom, x=parent.x, y=parent.y, z=parent.z)
            if atom.id == parent.id
            else atom
            for atom in atoms
        ),
    )


def xenon_model(
    pristine: Sequence[nteme_neb.Atom],
    cell: nteme_neb.Cell,
    route: nteme_neb.Barrier,
) -> Model:
    atoms, _center, final, route_meta = vacancy_route(pristine, cell, route, "Xe")
    return Model(
        "xenon-divacancy",
        tuple(atoms),
        cell,
        {
            "outcome": "executable-screening-structure",
            "route": route_meta,
            "parameters": {
                "charge_e": 0.0,
                "sigma_angstrom": XENON_SIGMA,
                "epsilon_kcal_mol": XENON_EPSILON,
                "mixing": "Lorentz-Berthelot with ClayFF",
                "source": "Bourg & Sposito 2008, DOI 10.1016/j.gca.2008.02.012; Gadikota et al. ClayFF clay-interlayer application",
            },
            "validation": "The Xe well is source-backed and has been used with ClayFF in hydrated clay interlayers; dry-muscovite cross interactions remain unvalidated. Treat barriers as a sensitivity screen until E3b or adsorption/solubility data calibrate them.",
            "neb_endpoint": asdict(final),
        },
        tuple(
            replace(atom, x=final.x, y=final.y, z=final.z)
            if atom.id == route.moving_site
            else atom
            for atom in atoms
        ),
    )


def type_tables() -> tuple[
    dict[str, int], dict[str, float], dict[str, float], dict[str, float]
]:
    names = (*nteme_neb.TYPE_ORDER, "Xe")
    type_id = {name: index for index, name in enumerate(names, 1)}
    mass = dict(nteme_neb.MASS)
    charge = dict(nteme_neb.CHARGE)
    sigma = dict(nteme_neb.SIGMA)
    epsilon = dict(nteme_neb.EPSILON)
    mass["Xe"] = XENON_MASS
    charge["Xe"] = 0.0
    sigma["Xe"] = XENON_SIGMA
    epsilon["Xe"] = XENON_EPSILON
    return type_id, mass, sigma, epsilon


def remaining_oh_bonds(
    atoms: Sequence[nteme_neb.Atom], cell: nteme_neb.Cell
) -> list[tuple[int, int]]:
    oxygens = [atom for atom in atoms if atom.element == "O1"]
    hydrogens = [atom for atom in atoms if atom.element == "H"]
    if len(oxygens) != len(hydrogens):
        raise ValueError(f"O1/H mismatch: {len(oxygens)} != {len(hydrogens)}")
    candidates = sorted(
        (
            periodic_distance(oxygen, hydrogen, cell),
            oxygen.id,
            hydrogen.id,
        )
        for oxygen in oxygens
        for hydrogen in hydrogens
        if periodic_distance(oxygen, hydrogen, cell) < 1.35
    )
    used_o: set[int] = set()
    used_h: set[int] = set()
    bonds = []
    for _distance, oxygen_id, hydrogen_id in candidates:
        if oxygen_id not in used_o and hydrogen_id not in used_h:
            used_o.add(oxygen_id)
            used_h.add(hydrogen_id)
            bonds.append((oxygen_id, hydrogen_id))
    if len(bonds) != len(oxygens):
        raise ValueError(f"could bind only {len(bonds)} of {len(oxygens)} OH groups")
    return sorted(bonds)


def write_lammps_data(path: pathlib.Path, model: Model) -> None:
    type_id, mass, _sigma, _epsilon = type_tables()
    bonds = remaining_oh_bonds(model.atoms, model.cell)
    lx, ly, lz, xy, xz, yz = model.cell.restricted()
    names = (*nteme_neb.TYPE_ORDER, "Xe")
    lines = [
        f"E3a {model.name}; source Nteme 2022 workbook {nteme_neb.WORKBOOK_SHA256}",
        "",
        f"{len(model.atoms)} atoms",
        f"{len(bonds)} bonds",
        "",
        f"{len(names)} atom types",
        "1 bond types",
        "",
        f"0.0 {lx:.12f} xlo xhi",
        f"0.0 {ly:.12f} ylo yhi",
        f"0.0 {lz:.12f} zlo zhi",
        f"{xy:.12f} {xz:.12f} {yz:.12f} xy xz yz",
        "",
        "Masses",
        "",
    ]
    for name in names:
        lines.append(f"{type_id[name]} {mass[name]:.10f} # {name}")
    lines.extend(("", "Atoms # full", ""))
    for atom in sorted(model.atoms, key=lambda item: item.id):
        lines.append(
            f"{atom.id} 1 {type_id[atom.element]} {atom_charge(atom.element):.8f} "
            f"{atom.x:.12f} {atom.y:.12f} {atom.z:.12f} # {atom.element}"
        )
    lines.extend(("", "Bonds", ""))
    for bond_id, (oxygen_id, hydrogen_id) in enumerate(bonds, 1):
        lines.append(f"{bond_id} 1 {oxygen_id} {hydrogen_id}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_static_input(path: pathlib.Path, data_name: str) -> None:
    type_id, _mass, sigma, epsilon = type_tables()
    names = (*nteme_neb.TYPE_ORDER, "Xe")
    lines = [
        "units           real",
        "atom_style      full",
        "atom_modify     map array sort 0 0.0",
        "boundary        p p p",
        f"read_data       {data_name}",
        "pair_style      lj/cut/coul/long 10.0",
    ]
    for name in names:
        lines.append(
            f"pair_coeff      {type_id[name]} {type_id[name]} {epsilon[name]:.10g} {sigma[name]:.10g}"
        )
    for index, name_i in enumerate(names):
        for name_j in names[index + 1 :]:
            mixed_eps = math.sqrt(epsilon[name_i] * epsilon[name_j])
            mixed_sig = 0.5 * (sigma[name_i] + sigma[name_j])
            lines.append(
                f"pair_coeff      {type_id[name_i]} {type_id[name_j]} "
                f"{mixed_eps:.10g} {mixed_sig:.10g}"
            )
    lines.extend(
        (
            "bond_style      harmonic",
            "bond_coeff      1 554.1349 1.0",
            "special_bonds   lj/coul 0.0 0.0 0.0",
            "kspace_style    ewald 1.0e-4",
            "neighbor        2.0 bin",
            "neigh_modify    delay 0 every 1 check yes",
            "thermo_style    custom step atoms bonds pe ebond evdwl ecoul elong press fmax",
            "thermo          1",
            "run             0",
        )
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def minimum_pair(
    atoms: Sequence[nteme_neb.Atom], cell: nteme_neb.Cell
) -> tuple[float, tuple[int, int]]:
    best = (math.inf, (-1, -1))
    for index, first in enumerate(atoms):
        for second in atoms[index + 1 :]:
            distance = periodic_distance(first, second, cell)
            if distance < best[0]:
                best = (distance, (first.id, second.id))
    return best


def validate_model(model: Model) -> dict[str, object]:
    ids = [atom.id for atom in model.atoms]
    if len(ids) != len(set(ids)):
        raise ValueError(f"{model.name}: duplicate atom IDs")
    if not all(
        math.isfinite(value)
        for atom in model.atoms
        for value in (atom.x, atom.y, atom.z)
    ):
        raise ValueError(f"{model.name}: non-finite coordinate")
    charge = net_charge(model.atoms)
    if not math.isclose(charge, 0.0, abs_tol=1.0e-8):
        raise ValueError(f"{model.name}: non-neutral structure q={charge:.12g}")
    restricted = model.cell.restricted()
    if not all(math.isfinite(value) and value > 0 for value in restricted[:3]):
        raise ValueError(f"{model.name}: invalid periodic cell")
    pair_distance, pair_ids = minimum_pair(model.atoms, model.cell)
    if pair_distance < MIN_DISTANCE_ANGSTROM:
        raise ValueError(
            f"{model.name}: minimum pair {pair_ids} is {pair_distance:.6f} A"
        )
    extra: dict[str, object] = {}
    if model.name == "dehydroxylate-lattice":
        transformation = model.metadata["transformation"]
        assert isinstance(transformation, dict)
        residual_id = int(transformation["residual_oxygen_site_id"])
        residual = next(atom for atom in model.atoms if atom.id == residual_id)
        octahedral_al = sorted(
            (periodic_distance(residual, atom, model.cell), atom.id)
            for atom in model.atoms
            if atom.element == "Al1"
            and periodic_distance(residual, atom, model.cell) <= 2.30
        )
        if len(octahedral_al) != 2:
            raise ValueError(
                f"{model.name}: residual O{residual_id} has {len(octahedral_al)} octahedral-Al neighbors"
            )
        extra["residual_oxygen"] = {
            "site_id": residual_id,
            "octahedral_al_neighbors": [
                atom_id for _distance, atom_id in octahedral_al
            ],
            "distances_angstrom": [distance for distance, _atom_id in octahedral_al],
        }
    return {
        "atom_count": len(model.atoms),
        "atom_counts": dict(
            sorted(Counter(atom.element for atom in model.atoms).items())
        ),
        "net_charge_e": charge,
        "minimum_pair_distance_angstrom": pair_distance,
        "minimum_pair_ids": list(pair_ids),
        "periodic_cell": asdict(model.cell),
        "cell_volume_angstrom3": restricted[0] * restricted[1] * restricted[2],
        **extra,
    }


def run_static(
    model_dir: pathlib.Path, lammps: str, input_name: str, label: str
) -> dict[str, object]:
    # argv sequence, shell=False: trusted local LAMMPS executable plus
    # static -log/-in flags and internally chosen input/label tokens.
    # nosemgrep: python.lang.security.audit.dangerous-subprocess-use-audit
    completed = subprocess.run(
        [lammps, "-log", f"{label}.lammps", "-in", input_name],
        cwd=model_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
        shell=False,
    )
    (model_dir / f"{label}.stdout.log").write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"LAMMPS static smoke {label} failed in {model_dir} with {completed.returncode}"
        )
    numeric: list[list[float]] = []
    for line in completed.stdout.splitlines():
        fields = line.split()
        if len(fields) != 10:
            continue
        try:
            numeric.append([float(field) for field in fields])
        except ValueError:
            pass
    if not numeric:
        raise RuntimeError(f"no static thermo row parsed in {model_dir}")
    row = numeric[-1]
    return {
        "exit_code": completed.returncode,
        "potential_energy_kcal_mol": row[3],
        "max_force_kcal_mol_angstrom": row[9],
        "log_sha256": hashlib.sha256(completed.stdout.encode()).hexdigest(),
    }


def build_models(workbook: pathlib.Path) -> tuple[list[Model], dict[str, object]]:
    atoms, cell, barriers = nteme_neb.read_workbook(workbook)
    route = nteme_neb.select_route(barriers, "divacancy", 1)
    models = [
        reconstructed_replication_model(atoms, cell, route),
        dehydroxylated_model(atoms, cell, route),
        extended_zone_model(atoms, cell, route),
        octahedral_trap_model(atoms, cell, route),
        xenon_model(atoms, cell, route),
    ]
    return models, {
        "source_workbook": str(workbook.resolve()),
        "source_sha256": nteme_neb.sha256(workbook),
        "source_cell": asdict(cell),
        "source_atoms": len(atoms),
        "published_route": asdict(route),
    }


def command_build(args: argparse.Namespace) -> None:
    models, source = build_models(args.workbook)
    args.out.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, object] = {
        "schema": "e3a-atomistic-input-contract-v1",
        "source": source,
        "models": {},
    }
    for model in models:
        model_dir = args.out / model.name
        model_dir.mkdir(parents=True, exist_ok=True)
        validation = validate_model(model)
        write_lammps_data(model_dir / "structure.data", model)
        write_static_input(model_dir / "in.static", "structure.data")
        if model.endpoint_atoms is None:
            raise ValueError(f"{model.name}: missing deterministic endpoint structure")
        endpoint_model = replace(model, atoms=model.endpoint_atoms, endpoint_atoms=None)
        endpoint_validation = validate_model(endpoint_model)
        write_lammps_data(model_dir / "endpoint.data", endpoint_model)
        write_static_input(model_dir / "in.endpoint.static", "endpoint.data")
        if args.no_lammps:
            static = None
            endpoint_static = None
        else:
            static = run_static(model_dir, args.lammps, "in.static", "initial-static")
            endpoint_static = run_static(
                model_dir, args.lammps, "in.endpoint.static", "endpoint-static"
            )
        record = {
            "metadata": model.metadata,
            "validation": validation,
            "static": static,
            "endpoint_validation": endpoint_validation,
            "endpoint_static": endpoint_static,
        }
        (model_dir / "manifest.json").write_text(
            json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        manifest["models"][model.name] = record
    (args.out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(required=True)
    build = sub.add_parser("build", help="build and validate all contract models")
    build.add_argument("--workbook", type=pathlib.Path, required=True)
    build.add_argument("--out", type=pathlib.Path, required=True)
    build.add_argument("--lammps", default="lmp")
    build.add_argument("--no-lammps", action="store_true")
    build.set_defaults(func=command_build)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
