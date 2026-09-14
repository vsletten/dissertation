#!/usr/bin/env python3
"""Run matched-cell LAMMPS energies on an immutable eight-image path.

This is the classical half of the E3b survey comparison. It evaluates the same
coordinates used by the PBE-D3 frozen-path profile without relaxing them. The
result is a like-for-like path energy rise, not a minimum-energy-path barrier.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
import pathlib
import re
import shutil
import subprocess
import sys
import time
from typing import Any

import e3b_periodic_dft as e3b

IMAGE_COUNT = 8
MODELS = (
    "reconstructed-replication",
    "dehydroxylate-lattice",
    "xenon-divacancy",
)


def _atomic_json(path: pathlib.Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def _file_record(path: pathlib.Path, root: pathlib.Path) -> dict[str, Any]:
    resolved = path.resolve()
    return {
        "path": str(resolved.relative_to(root.resolve())),
        "bytes": resolved.stat().st_size,
        "sha256": e3b.sha256(resolved),
    }


def parse_xyz(
    path: pathlib.Path, atom_map: list[dict[str, Any]]
) -> dict[int, tuple[float, float, float]]:
    """Bind ordered XYZ coordinates to explicit source atom IDs."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"{path}: empty XYZ")
    try:
        expected = int(lines[0])
    except ValueError as exc:
        raise ValueError(f"{path}: malformed XYZ atom count") from exc
    rows = lines[2:]
    if expected != len(rows) or expected != len(atom_map):
        raise ValueError(f"{path}: XYZ/identity atom count mismatch")
    coordinates: dict[int, tuple[float, float, float]] = {}
    for index, (line, identity) in enumerate(zip(rows, atom_map, strict=True), 1):
        fields = line.split()
        if len(fields) != 4:
            raise ValueError(f"{path}: malformed XYZ row {index}")
        cp2k_index = identity.get("cp2k_index")
        source_id = identity.get("source_id")
        source_type = identity.get("source_type")
        if cp2k_index != index or not isinstance(source_id, int):
            raise ValueError(f"{path}: malformed atom identity row {index}")
        if fields[0] != source_type:
            raise ValueError(f"{path}: XYZ label mismatch at row {index}")
        try:
            xyz = tuple(float(value) for value in fields[1:])
        except ValueError as exc:
            raise ValueError(f"{path}: non-numeric XYZ row {index}") from exc
        if len(xyz) != 3 or not all(math.isfinite(value) for value in xyz):
            raise ValueError(f"{path}: non-finite XYZ row {index}")
        if source_id in coordinates:
            raise ValueError(f"{path}: duplicate source atom id {source_id}")
        coordinates[source_id] = xyz
    return coordinates


def render_lammps_data(
    template: str, coordinates: dict[int, tuple[float, float, float]]
) -> str:
    """Replace only atom coordinates in a prepared LAMMPS data document."""
    lines = template.splitlines()
    headers = [
        index for index, line in enumerate(lines) if line.strip() == "Atoms # full"
    ]
    if len(headers) != 1:
        raise ValueError("LAMMPS template must contain one 'Atoms # full' section")
    start = headers[0] + 1
    replaced: set[int] = set()
    in_rows = False
    for index in range(start, len(lines)):
        stripped = lines[index].strip()
        if not stripped:
            if in_rows:
                break
            continue
        fields = lines[index].split()
        if not fields[0].isdigit():
            if in_rows:
                break
            continue
        in_rows = True
        if len(fields) < 7:
            raise ValueError(f"malformed LAMMPS atom row: {lines[index]}")
        atom_id = int(fields[0])
        if atom_id not in coordinates or atom_id in replaced:
            raise ValueError(f"LAMMPS atom identity mismatch: {atom_id}")
        x, y, z = coordinates[atom_id]
        fields[4:7] = [f"{x:.12f}", f"{y:.12f}", f"{z:.12f}"]
        lines[index] = " ".join(fields)
        replaced.add(atom_id)
    if replaced != set(coordinates):
        missing = sorted(set(coordinates) - replaced)
        extra = sorted(replaced - set(coordinates))
        raise ValueError(
            f"LAMMPS/XYZ atom identity mismatch: missing={missing}, extra={extra}"
        )
    return "\n".join(lines) + "\n"


def render_static_input(template: str, data_name: str) -> str:
    """Turn the prepared minimization deck into a zero-step energy deck."""
    lines = template.splitlines()
    read_indices = [
        index
        for index, line in enumerate(lines)
        if line.strip().startswith("read_data ")
    ]
    if len(read_indices) != 1:
        raise ValueError("LAMMPS template must contain one read_data command")
    output: list[str] = []
    for index, line in enumerate(lines):
        stripped = line.strip()
        if index == read_indices[0]:
            output.append(f"read_data       {data_name}")
        elif stripped.startswith(("minimize ", "write_data ", "write_dump ")):
            continue
        elif stripped.startswith("thermo "):
            output.append("thermo          1")
        else:
            output.append(line)
    output.extend(("run             0", ""))
    return "\n".join(output)


def parse_energy(path: pathlib.Path) -> float:
    text = path.read_text(encoding="utf-8", errors="replace")
    if re.search(r"(?im)^\s*(?:ERROR|FATAL)\b", text):
        raise ValueError(f"{path}: LAMMPS fatal/error output")
    if "total wall time:" not in text.lower() and "loop time of" not in text.lower():
        raise ValueError(f"{path}: missing LAMMPS normal termination")
    lines = text.splitlines()
    expected_header = ("Step", "PotEng", "Fnorm", "Fmax", "Press")
    header_indices = [
        index
        for index, line in enumerate(lines)
        if tuple(line.split()) == expected_header
    ]
    if len(header_indices) != 1:
        raise ValueError(f"{path}: expected exactly one declared thermo header")
    row_index = next(
        (
            index
            for index in range(header_indices[0] + 1, len(lines))
            if lines[index].strip()
        ),
        None,
    )
    if row_index is None:
        raise ValueError(f"{path}: declared thermo header has no row")
    fields = lines[row_index].split()
    if len(fields) != len(expected_header):
        raise ValueError(f"{path}: malformed declared thermo row")
    try:
        row = [float(field) for field in fields]
    except ValueError as exc:
        raise ValueError(f"{path}: non-numeric declared thermo row") from exc
    if row[0] != 0.0:
        raise ValueError(f"{path}: static energy row is not step zero")
    energy = row[1]
    if not math.isfinite(energy):
        raise ValueError(f"{path}: non-finite potential energy")
    return energy


def run_profile(args: argparse.Namespace) -> dict[str, Any]:
    prepared = args.prepared_root.resolve()
    runtime = args.runtime_root.resolve()
    receipt_path = args.receipt.resolve()
    if runtime.exists() or receipt_path.exists():
        raise FileExistsError("runtime root and receipt must both be new")
    if args.per_image_timeout <= 0 or args.per_image_timeout > 300:
        raise ValueError("per-image timeout must be in (0, 300] seconds")
    if args.model not in MODELS:
        raise ValueError(f"unsupported model: {args.model}")

    e3b.verify_manifest(prepared)
    preparation = e3b.read_json(prepared / "preparation.json")
    models = preparation.get("models")
    if not isinstance(models, dict) or args.model not in models:
        raise ValueError("model is absent from prepared evidence")
    model_record = models[args.model]
    if not isinstance(model_record, dict):
        raise TypeError("prepared model record is malformed")

    model_root = prepared / "models" / args.model
    atom_map_value = e3b.read_json(model_root / "atom-map.json").get("atoms")
    if not isinstance(atom_map_value, list) or not all(
        isinstance(item, dict) for item in atom_map_value
    ):
        raise TypeError("prepared atom map is malformed")
    atom_map: list[dict[str, Any]] = atom_map_value
    atom_identity_sha256 = e3b.canonical_sha256(atom_map)
    if atom_identity_sha256 != model_record.get("atom_identity_sha256"):
        raise ValueError("prepared atom map does not match declared atom identity")
    images = [
        model_root / "images" / f"replica-{index:02d}.xyz"
        for index in range(IMAGE_COUNT)
    ]
    if any(not path.is_file() for path in images):
        raise ValueError("prepared frozen path does not contain exactly eight images")

    matched = model_root / "matched-classical"
    data_template_path = matched / "initial.seed.data"
    input_template_path = matched / "in.min.initial"
    data_template = data_template_path.read_text(encoding="utf-8")
    input_template = input_template_path.read_text(encoding="utf-8")
    lammps = pathlib.Path(shutil.which(args.lammps) or args.lammps).resolve()
    if not lammps.is_file():
        raise FileNotFoundError(args.lammps)

    runtime.mkdir(parents=True)
    records: list[dict[str, Any]] = []
    status = "complete"
    failure: str | None = None
    started_at = dt.datetime.now(dt.timezone.utc)
    started = time.monotonic()
    env = os.environ.copy()
    env.update(
        {
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }
    )
    for index, coordinate_path in enumerate(images):
        data_path = runtime / f"image-{index:02d}.data"
        input_path = runtime / f"image-{index:02d}.in"
        stdout_path = runtime / f"image-{index:02d}.stdout.log"
        coordinates = parse_xyz(coordinate_path, atom_map)
        data_path.write_text(
            render_lammps_data(data_template, coordinates), encoding="utf-8"
        )
        input_path.write_text(
            render_static_input(input_template, data_path.name), encoding="utf-8"
        )
        command = [str(lammps), "-log", "none", "-in", input_path.name]
        image_started = time.monotonic()
        timed_out = False
        with stdout_path.open("wb") as output:
            try:
                completed = subprocess.run(
                    command,
                    cwd=runtime,
                    env=env,
                    stdout=output,
                    stderr=subprocess.STDOUT,
                    timeout=args.per_image_timeout,
                    check=False,
                )
                returncode = completed.returncode
            except subprocess.TimeoutExpired:
                timed_out = True
                returncode = -9
            output.flush()
            os.fsync(output.fileno())
        energy: float | None = None
        image_status = "incomplete"
        try:
            if returncode != 0 or timed_out:
                raise ValueError(f"LAMMPS exited {returncode}; timed_out={timed_out}")
            energy = parse_energy(stdout_path)
            image_status = "complete"
        except ValueError as exc:
            status = "incomplete-image"
            failure = f"image {index:02d}: {exc}"
        records.append(
            {
                "index": index,
                "status": image_status,
                "returncode": returncode,
                "timed_out": timed_out,
                "energy_kcal_mol": energy,
                "elapsed_seconds": time.monotonic() - image_started,
                "coordinate": _file_record(coordinate_path, prepared),
                "data": _file_record(data_path, runtime),
                "input": _file_record(input_path, runtime),
                "stdout": _file_record(stdout_path, runtime),
                "command": command,
            }
        )
        if status != "complete":
            break

    complete = (
        status == "complete"
        and len(records) == IMAGE_COUNT
        and all(
            isinstance(item["energy_kcal_mol"], (int, float))
            and math.isfinite(float(item["energy_kcal_mol"]))
            for item in records
        )
    )
    energies = (
        [float(item["energy_kcal_mol"]) for item in records] if complete else None
    )
    profile_rise = max(energies) - energies[0] if energies is not None else None
    if complete and (
        profile_rise is None or not math.isfinite(profile_rise) or profile_rise < 0
    ):
        complete = False
        status = "incomplete-numeric"
        failure = "profile rise is not finite and non-negative"
        energies = None
        profile_rise = None

    e3b.verify_manifest(prepared)
    runtime_files = [
        _file_record(path, runtime)
        for path in sorted(runtime.rglob("*"))
        if path.is_file() and path.resolve() != receipt_path
    ]
    receipt = {
        "schema": "e3b-classical-frozen-path-profile-v1",
        "operator": "(hermes-custom-build-001; profile=workstation)",
        "purpose": (
            "matched-cell LAMMPS energies on the immutable PBE-D3 survey path; "
            "not endpoint relaxation, classical NEB, or a production barrier"
        ),
        "status": "complete" if complete else status,
        "failure": None if complete else failure,
        "started_at": started_at.isoformat(),
        "completed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "elapsed_seconds": time.monotonic() - started,
        "model": args.model,
        "atom_identity_sha256": atom_identity_sha256,
        "prepared": {
            "path": str(prepared),
            "preparation_sha256": e3b.sha256(prepared / "preparation.json"),
            "manifest_sha256": e3b.sha256(prepared / "manifest.json"),
            "data_template_sha256": e3b.sha256(data_template_path),
            "input_template_sha256": e3b.sha256(input_template_path),
        },
        "method": {
            "basis": "matched-cell-classical-frozen-path-single-points",
            "lammps": str(lammps),
            "lammps_sha256": e3b.sha256(lammps),
            "units": "real",
            "image_count": IMAGE_COUNT,
        },
        "controller": {
            "type": "subprocess-timeout-per-image",
            "per_image_timeout_seconds": args.per_image_timeout,
            "omp_threads": 1,
        },
        "images": records,
        "profile_energies_kcal_mol": energies,
        "profile_rise_kcal_mol": profile_rise,
        "claim_boundary": {
            "supports": "classical energy rise along the exact tested immutable path",
            "does_not_support": [
                "relaxed classical endpoint energies",
                "classical minimum-energy path",
                "a production activation barrier",
            ],
        },
        "runtime_root": str(runtime),
        "runtime_manifest": {
            "files": runtime_files,
            "sha256": e3b.canonical_sha256(runtime_files),
        },
    }
    _atomic_json(receipt_path, receipt)
    return receipt


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared-root", type=pathlib.Path, required=True)
    parser.add_argument("--runtime-root", type=pathlib.Path, required=True)
    parser.add_argument("--receipt", type=pathlib.Path, required=True)
    parser.add_argument("--model", choices=MODELS, required=True)
    parser.add_argument("--lammps", default="lmp")
    parser.add_argument("--per-image-timeout", type=float, default=120.0)
    return parser


def main() -> int:
    receipt = run_profile(build_parser().parse_args())
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0 if receipt["status"] == "complete" else 2


if __name__ == "__main__":
    sys.exit(main())
