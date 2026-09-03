#!/usr/bin/env python3
"""Run the bounded classical-potential E3a sensitivity campaign.

The runner consumes the hash-pinned atomistic models emitted by
``e3a_input_contract.py``.  It relaxes both endpoints, checks model-specific
minimum gates, runs an eight-image LAMMPS CI-NEB only for eligible models, and
writes one machine-readable receipt per model plus an aggregate receipt.

These calculations are classical sensitivity results.  They do not erase the
input contract's reconstructed-charge-compensation or Xe-validation limits.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import pathlib
import shutil
import subprocess
import sys
from collections.abc import Iterable, Sequence

import e3a_input_contract as contract

MODEL_ORDER = (
    "reconstructed-replication",
    "dehydroxylate-lattice",
    "extended-defect-zone-small",
    "extended-defect-zone",
    "extended-defect-zone-large",
    "octahedral-trap-escape",
    "xenon-divacancy",
)

EXTENDED_ZONE_GRID = (
    ("extended-defect-zone-small", 6.0, 3.0, 0.5),
    ("extended-defect-zone", 8.0, 5.0, 1.0),
    ("extended-defect-zone-large", 10.0, 7.0, 1.5),
)


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def potential_lines() -> list[str]:
    type_id, _mass, sigma, epsilon = contract.type_tables()
    names = (*contract.nteme_neb.TYPE_ORDER, "Xe")
    lines = ["pair_style      lj/cut/coul/long 10.0"]
    for name in names:
        lines.append(
            f"pair_coeff      {type_id[name]} {type_id[name]} "
            f"{epsilon[name]:.10g} {sigma[name]:.10g}"
        )
    for index, name_i in enumerate(names):
        for name_j in names[index + 1 :]:
            mixed_epsilon = math.sqrt(epsilon[name_i] * epsilon[name_j])
            mixed_sigma = 0.5 * (sigma[name_i] + sigma[name_j])
            lines.append(
                f"pair_coeff      {type_id[name_i]} {type_id[name_j]} "
                f"{mixed_epsilon:.10g} {mixed_sigma:.10g}"
            )
    lines.extend(
        (
            "bond_style      harmonic",
            "bond_coeff      1 554.1349 1.0",
            "special_bonds   lj/coul 0.0 0.0 0.0",
            "kspace_style    ewald 1.0e-4",
            "neighbor        2.0 bin",
            "neigh_modify    delay 0 every 1 check yes",
        )
    )
    return lines


def _id_chunks(ids: Iterable[int], width: int = 120) -> list[str]:
    tokens = [str(atom_id) for atom_id in sorted(set(ids))]
    if not tokens:
        return []
    chunks: list[str] = []
    current = ""
    for token in tokens:
        candidate = f"{current} {token}".strip()
        if current and len(candidate) > width:
            chunks.append(current)
            current = token
        else:
            current = candidate
    chunks.append(current)
    return chunks


def write_minimize_input(
    path: pathlib.Path,
    *,
    data_name: str,
    output_prefix: str,
    ftol: float,
    max_steps: int,
    local_mobile_ids: Sequence[int] = (),
) -> None:
    lines = [
        "units           real",
        "atom_style      full",
        "atom_modify     map array sort 0 0.0",
        "boundary        p p p",
        "pair_style      lj/cut/coul/long 10.0",
        "bond_style      harmonic",
        f"read_data       {data_name}",
        *(
            line
            for line in potential_lines()
            if not line.startswith(("pair_style", "bond_style"))
        ),
        "thermo_style    custom step pe fnorm fmax press",
        "thermo          100",
        "timestep        0.1",
        "min_style       fire",
    ]
    chunks = _id_chunks(local_mobile_ids)
    if chunks:
        lines.append(f"group           local-mobile id {chunks[0]} &")
        for index, chunk in enumerate(chunks[1:], 1):
            suffix = " &" if index < len(chunks) - 1 else ""
            lines.append(f"                {chunk}{suffix}")
        lines.extend(
            (
                "group           local-frozen subtract all local-mobile",
                "fix             local_hold local-frozen setforce 0.0 0.0 0.0",
                f"minimize        0.0 0.05 {max_steps} {max_steps * 10}",
                "unfix           local_hold",
            )
        )
    lines.extend(
        (
            f"minimize        0.0 {ftol:.8g} {max_steps} {max_steps * 10}",
            f"write_data      {output_prefix}.data",
            f"write_dump      all custom {output_prefix}.dump id type x y z modify sort id",
        )
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_neb_input(
    path: pathlib.Path,
    *,
    data_name: str,
    final_name: str,
    ftol: float,
    relax_steps: int,
    climb_steps: int,
) -> None:
    lines = [
        "units           real",
        "atom_style      full",
        "atom_modify     map array sort 0 0.0",
        "boundary        p p p",
        "pair_style      lj/cut/coul/long 10.0",
        "bond_style      harmonic",
        f"read_data       {data_name}",
        *(
            line
            for line in potential_lines()
            if not line.startswith(("pair_style", "bond_style"))
        ),
        "thermo_style    custom step pe ebond evdwl ecoul elong fnorm fmax",
        "thermo          100",
        "timestep        0.1",
        "min_style       quickmin",
        "variable        u uloop 8",
        "fix             path all neb 23.0",
        "dump            path all custom 100 dump.neb.$u id type x y z",
        f"neb             1.0e-10 {ftol:.8g} {relax_steps} {climb_steps} 100 final {final_name}",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_dump(path: pathlib.Path) -> dict[int, tuple[float, float, float]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    try:
        atom_header = lines.index("ITEM: ATOMS id type x y z")
    except ValueError as exc:
        raise ValueError(f"{path}: expected sorted custom dump columns") from exc
    atoms: dict[int, tuple[float, float, float]] = {}
    for line in lines[atom_header + 1 :]:
        fields = line.split()
        if len(fields) != 5:
            continue
        atom_id = int(fields[0])
        x, y, z = (float(value) for value in fields[2:5])
        atoms[atom_id] = (x, y, z)
    expected = int(lines[lines.index("ITEM: NUMBER OF ATOMS") + 1])
    if len(atoms) != expected:
        raise ValueError(f"{path}: parsed {len(atoms)} atoms, expected {expected}")
    return atoms


def write_final_coords(
    path: pathlib.Path, coords: dict[int, tuple[float, float, float]]
) -> None:
    lines = [str(len(coords))]
    for atom_id, (x, y, z) in sorted(coords.items()):
        lines.append(f"{atom_id} {x:.16g} {y:.16g} {z:.16g}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_last_thermo(path: pathlib.Path) -> dict[str, float]:
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) != 5:
            continue
        try:
            rows.append([float(field) for field in fields])
        except ValueError:
            continue
    if not rows:
        raise ValueError(f"{path}: no five-column minimization thermo rows")
    row = rows[-1]
    return {
        "step": row[0],
        "potential_energy_kcal_mol": row[1],
        "force_norm_kcal_mol_angstrom": row[2],
        "max_force_kcal_mol_angstrom": row[3],
        "pressure_atm": row[4],
    }


def parse_neb(path: pathlib.Path, ftol: float) -> dict[str, object]:
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) != 25:
            continue
        try:
            rows.append([float(field) for field in fields])
        except ValueError:
            continue
    if not rows:
        raise ValueError(f"{path}: no 25-column NEB progress rows")
    row = rows[-1]
    final_step = int(row[0])
    recent = [item[6] for item in rows if int(item[0]) >= final_step - 500]
    return {
        "final_step": final_step,
        "max_replica_force_kcal_mol_angstrom": row[1],
        "max_atom_force_kcal_mol_angstrom": row[2],
        "computed_forward_barrier_kcal_mol": row[6],
        "computed_reverse_barrier_kcal_mol": row[7],
        "barrier_span_last_500_steps_kcal_mol": max(recent) - min(recent),
        "converged_to_requested_ftol": row[1] <= ftol,
        "screen_sha256": sha256(path),
    }


def _run(argv: Sequence[str], *, cwd: pathlib.Path, log: pathlib.Path) -> None:
    env = os.environ.copy()
    env.update(
        {
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    with log.open("w", encoding="utf-8") as handle:
        # argv sequence, shell=False: trusted lmp/mpirun plus internally
        # built input-deck tokens. Not a shell string.
        # nosemgrep: python.lang.security.audit.dangerous-subprocess-use-audit
        completed = subprocess.run(
            list(argv),
            cwd=cwd,
            env=env,
            text=True,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
            shell=False,
        )
    if completed.returncode != 0:
        raise RuntimeError(f"command exited {completed.returncode}: {argv}; see {log}")


def model_atom(model: contract.Model, atom_id: int) -> contract.nteme_neb.Atom:
    return next(atom for atom in model.atoms if atom.id == atom_id)


def local_mobile_ids(model: contract.Model, radius: float = 4.25) -> list[int]:
    parent_id = int(str(model.metadata["parent_k_site_id"]))
    center = model_atom(model, parent_id)
    return [
        atom.id
        for atom in model.atoms
        if contract.periodic_distance(center, atom, model.cell) <= radius
    ]


def trap_gate(
    model: contract.Model, coords: dict[int, tuple[float, float, float]]
) -> dict[str, object]:
    parent_id = int(str(model.metadata["parent_k_site_id"]))
    original = model_atom(model, parent_id)
    relaxed = contract.nteme_neb.Atom(parent_id, "Ar", *coords[parent_id])
    displacement = contract.periodic_distance(original, relaxed, model.cell)
    return {
        "ar_displacement_from_cage_seed_angstrom": displacement,
        "ar_remains_in_seed_cage": displacement <= 1.5,
        "index_zero": None,
        "index_gate": "not-run-full-system; NEB forbidden until independently proven",
    }


def build_campaign_models(
    workbook: pathlib.Path,
) -> tuple[list[contract.Model], dict[str, object]]:
    models, source = contract.build_models(workbook)
    by_name = {model.name: model for model in models}
    pristine, cell, barriers = contract.nteme_neb.read_workbook(workbook)
    route = contract.nteme_neb.select_route(barriers, "divacancy", 1)
    for name, radius_a, radius_b, opening in EXTENDED_ZONE_GRID:
        if name == "extended-defect-zone":
            continue
        model = contract.extended_zone_model(
            pristine,
            cell,
            route,
            radius_a=radius_a,
            radius_b=radius_b,
            opening=opening,
        )
        by_name[name] = contract.replace(model, name=name)
    return [by_name[name] for name in MODEL_ORDER], source


def command_run(args: argparse.Namespace) -> None:
    if args.neb_relax_steps % 100 or args.neb_climb_steps % 100:
        raise ValueError("NEB relax/climb steps must be multiples of thermo every=100")
    models, source = build_campaign_models(args.workbook)
    by_name = {model.name: model for model in models}
    selected = MODEL_ORDER if not args.model else tuple(args.model)
    unknown = sorted(set(selected) - set(by_name))
    if unknown:
        raise ValueError(f"unknown model(s): {unknown}")
    args.out.mkdir(parents=True, exist_ok=True)
    model_results: dict[str, object] = {}
    aggregate: dict[str, object] = {
        "schema": "e3a-classical-neb-campaign-v1",
        "source": source,
        "settings": {
            "min_ftol_kcal_mol_angstrom": args.min_ftol,
            "min_steps": args.min_steps,
            "neb_ftol_kcal_mol_angstrom": args.neb_ftol,
            "neb_relax_steps": args.neb_relax_steps,
            "neb_climb_steps": args.neb_climb_steps,
            "images": 8,
            "spring_kcal_mol_angstrom2": 23.0,
            "omp_threads_per_rank": 1,
        },
        "models": model_results,
    }
    for name in selected:
        model = by_name[name]
        if model.endpoint_atoms is None:
            raise ValueError(f"{name}: missing endpoint atoms")
        model_dir = args.out / name
        model_dir.mkdir(parents=True, exist_ok=True)
        contract.write_lammps_data(model_dir / "initial.seed.data", model)
        endpoint_model = contract.replace(
            model, atoms=model.endpoint_atoms, endpoint_atoms=None
        )
        contract.write_lammps_data(model_dir / "endpoint.seed.data", endpoint_model)
        mobile = local_mobile_ids(model) if name == "octahedral-trap-escape" else []
        for label in ("initial", "endpoint"):
            write_minimize_input(
                model_dir / f"in.min.{label}",
                data_name=f"{label}.seed.data",
                output_prefix=f"{label}.relaxed",
                ftol=args.min_ftol,
                max_steps=args.min_steps,
                local_mobile_ids=mobile,
            )
            _run(
                (args.lammps, "-in", f"in.min.{label}"),
                cwd=model_dir,
                log=model_dir / f"{label}.min.stdout.log",
            )
        initial_min = parse_last_thermo(model_dir / "initial.min.stdout.log")
        endpoint_min = parse_last_thermo(model_dir / "endpoint.min.stdout.log")
        initial_coords = parse_dump(model_dir / "initial.relaxed.dump")
        endpoint_coords = parse_dump(model_dir / "endpoint.relaxed.dump")
        gates: dict[str, object] = {
            "initial_minimum_force_pass": initial_min["max_force_kcal_mol_angstrom"]
            <= args.min_ftol,
            "endpoint_minimum_force_pass": endpoint_min["max_force_kcal_mol_angstrom"]
            <= args.min_ftol,
        }
        if name == "octahedral-trap-escape":
            gates.update(trap_gate(model, initial_coords))
        eligible = bool(
            gates["initial_minimum_force_pass"]
            and gates["endpoint_minimum_force_pass"]
            and name != "octahedral-trap-escape"
        )
        result: dict[str, object] = {
            "provenance": model.metadata,
            "initial_minimum": initial_min,
            "endpoint_minimum": endpoint_min,
            "gates": gates,
            "neb": None,
        }
        if eligible:
            write_final_coords(model_dir / "final.coords", endpoint_coords)
            write_neb_input(
                model_dir / "in.neb",
                data_name="initial.relaxed.data",
                final_name="final.coords",
                ftol=args.neb_ftol,
                relax_steps=args.neb_relax_steps,
                climb_steps=args.neb_climb_steps,
            )
            _run(
                (
                    args.mpirun,
                    "-np",
                    "8",
                    args.lammps,
                    "-partition",
                    "8x1",
                    "-in",
                    "in.neb",
                ),
                cwd=model_dir,
                log=model_dir / "screen.log",
            )
            result["neb"] = parse_neb(model_dir / "screen.log", args.neb_ftol)
        else:
            result["typed_outcome"] = (
                "incomplete-index-gate"
                if name == "octahedral-trap-escape"
                else "incomplete-minimum"
            )
        receipt = model_dir / "result.json"
        receipt.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        model_results[name] = result
    aggregate_path = args.out / "campaign-result.json"
    aggregate_path.write_text(
        json.dumps(aggregate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest = []
    for path in sorted(item for item in args.out.rglob("*") if item.is_file()):
        manifest.append(
            {
                "path": str(path.relative_to(args.out)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    (args.out / "manifest.json").write_text(
        json.dumps({"files": manifest}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(aggregate, indent=2, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(required=True)
    run = sub.add_parser(
        "run", help="relax eligible models and execute bounded CI-NEBs"
    )
    run.add_argument("--workbook", type=pathlib.Path, required=True)
    run.add_argument("--out", type=pathlib.Path, required=True)
    run.add_argument("--model", action="append", choices=MODEL_ORDER)
    run.add_argument("--lammps", default=shutil.which("lmp") or "lmp")
    run.add_argument("--mpirun", default=shutil.which("mpirun") or "mpirun")
    run.add_argument("--min-ftol", type=float, default=0.01)
    run.add_argument("--min-steps", type=int, default=5000)
    run.add_argument("--neb-ftol", type=float, default=0.01)
    run.add_argument("--neb-relax-steps", type=int, default=3000)
    run.add_argument("--neb-climb-steps", type=int, default=10000)
    run.set_defaults(func=command_run)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
