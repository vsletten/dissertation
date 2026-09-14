#!/usr/bin/env python3
"""Run one bounded eight-image frozen-path CP2K profile.

This is a survey calibration route for the one-workstation platform phase.  It
computes single-point PBE-D3 energies on the immutable prepared classical path;
it does not optimize endpoints or claim a DFT CI-NEB barrier.
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
import signal
import subprocess
import sys
import time
from typing import Any

import e3b_periodic_dft as e3b

IMAGE_COUNT = 8
HARTREE_TO_KCAL_MOL = e3b.HARTREE_TO_KCAL_MOL
MODELS = (
    "reconstructed-replication",
    "dehydroxylate-lattice",
    "xenon-divacancy",
)


class Interrupted(RuntimeError):
    """Raised when the bounded controller asks the runner to stop."""


def _signal_handler(signum: int, _frame: object) -> None:
    raise Interrupted(f"received signal {signum}")


def _file_record(path: pathlib.Path, root: pathlib.Path) -> dict[str, Any]:
    resolved = path.resolve()
    return {
        "path": str(resolved.relative_to(root.resolve())),
        "bytes": resolved.stat().st_size,
        "sha256": e3b.sha256(resolved),
    }


def _atomic_json(path: pathlib.Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def render_input(template: str, model: str, index: int) -> str:
    """Render an immutable classical image as a CP2K ENERGY_FORCE input."""
    if model not in MODELS:
        raise ValueError(f"unsupported model: {model}")
    if index not in range(IMAGE_COUNT):
        raise ValueError(f"image index out of range: {index}")
    project_pattern = re.compile(r"(?m)^\s*PROJECT\s+\S+\s*$")
    coord_pattern = re.compile(r"(?m)^\s*COORD_FILE_NAME\s+\S+\s*$")
    guess_pattern = re.compile(r"(?m)^\s*SCF_GUESS\s+\S+\s*$")
    if len(project_pattern.findall(template)) != 1:
        raise ValueError("smoke template must contain one PROJECT declaration")
    if len(coord_pattern.findall(template)) != 1:
        raise ValueError("smoke template must contain one COORD_FILE_NAME declaration")
    if len(guess_pattern.findall(template)) != 1:
        raise ValueError("smoke template must contain one SCF_GUESS declaration")
    text = project_pattern.sub(f"  PROJECT e3b-{model}-frozen-profile", template)
    text = coord_pattern.sub(
        f"      COORD_FILE_NAME images/replica-{index:02d}.xyz", text
    )
    guess = "ATOMIC" if index == 0 else "RESTART"
    return guess_pattern.sub(f"      SCF_GUESS {guess}", text)


def systemd_readback(unit: str) -> dict[str, str]:
    keys = (
        "Id",
        "ActiveState",
        "SubState",
        "RuntimeMaxUSec",
        "MemoryMax",
        "Nice",
        "CPUQuotaPerSecUSec",
        "ControlGroup",
        "InvocationID",
    )
    command = [
        "systemctl",
        "--user",
        "show",
        unit,
        *[f"--property={key}" for key in keys],
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    values: dict[str, str] = {}
    for line in result.stdout.splitlines():
        key, separator, value = line.partition("=")
        if separator:
            values[key] = value
    if values.get("Id") != unit or values.get("ActiveState") != "active":
        raise ValueError("exact systemd unit is not active")
    if values.get("Nice") != "10" or values.get("CPUQuotaPerSecUSec") != "16s":
        raise ValueError("systemd CPU/nice envelope is not exact")
    if not values.get("ControlGroup") or not values.get("InvocationID"):
        raise ValueError("systemd controller identity is incomplete")
    return values


def _stop_process(
    process: subprocess.Popen[Any], grace: float = 20.0
) -> dict[str, bool]:
    term_sent = False
    kill_sent = False
    if process.poll() is None:
        os.killpg(process.pid, signal.SIGTERM)
        term_sent = True
        try:
            process.wait(timeout=grace)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            kill_sent = True
            process.wait(timeout=10)
    return {"process_group_term_sent": term_sent, "process_group_kill_sent": kill_sent}


def _remove_container(docker: str, cidfile: pathlib.Path) -> dict[str, Any]:
    container_id = (
        cidfile.read_text(encoding="utf-8").strip() if cidfile.is_file() else ""
    )
    if container_id:
        subprocess.run(
            [docker, "rm", "-f", container_id],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    cidfile.unlink(missing_ok=True)
    absent = True
    if container_id:
        absent = (
            subprocess.run(
                [docker, "inspect", container_id],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=False,
            ).returncode
            != 0
        )
    return {
        "container_id": container_id or None,
        "container_absent": absent,
        "cidfile_absent": not cidfile.exists(),
    }


def run_profile(args: argparse.Namespace) -> dict[str, Any]:
    prepared = args.prepared_root.resolve()
    runtime = args.runtime_root.resolve()
    receipt_path = args.receipt.resolve()
    if runtime.exists() or receipt_path.exists():
        raise FileExistsError("runtime root and receipt must both be new")
    if args.timeout_seconds <= 0 or args.timeout_seconds > 13800:
        raise ValueError("total timeout must be in (0, 13800] seconds")
    if args.per_image_timeout <= 0 or args.per_image_timeout > 3600:
        raise ValueError("per-image timeout must be in (0, 3600] seconds")
    if args.mpi_ranks * args.omp_threads not in range(1, 17):
        raise ValueError("MPI ranks times OpenMP threads must be in [1, 16]")
    if args.memory_gib not in range(1, 49):
        raise ValueError("container memory must be in [1, 48] GiB")

    e3b.verify_manifest(prepared)
    preparation = e3b.read_json(prepared / "preparation.json")
    models = preparation.get("models")
    if not isinstance(models, dict) or args.model not in models:
        raise ValueError("model is absent from prepared evidence")
    model_record = models[args.model]
    if not isinstance(model_record, dict):
        raise TypeError("prepared model record is malformed")
    image_dir = prepared / "models" / args.model / "images"
    prepared_images = [
        image_dir / f"replica-{index:02d}.xyz" for index in range(IMAGE_COUNT)
    ]
    if any(not path.is_file() for path in prepared_images):
        raise ValueError("prepared frozen path does not contain exactly eight images")

    controller = systemd_readback(args.systemd_unit)
    docker = shutil.which(args.docker)
    if docker is None:
        raise FileNotFoundError(args.docker)
    runtime.mkdir(parents=True)
    work = runtime / "work"
    shutil.copytree(prepared, work, copy_function=shutil.copy2)
    model_dir = work / "models" / args.model
    template_path = model_dir / "smoke.inp"
    template = template_path.read_text(encoding="utf-8")
    image_records: list[dict[str, Any]] = []
    status = "complete"
    failure: str | None = None
    started_at = dt.datetime.now(dt.timezone.utc)
    started = time.monotonic()
    active_process: subprocess.Popen[Any] | None = None
    active_cidfile: pathlib.Path | None = None

    signal.signal(signal.SIGTERM, _signal_handler)
    signal.signal(signal.SIGINT, _signal_handler)
    try:
        for index in range(IMAGE_COUNT):
            elapsed_total = time.monotonic() - started
            remaining = args.timeout_seconds - elapsed_total
            if remaining <= 60:
                status = "incomplete-total-timeout"
                failure = "total profile envelope exhausted before next image"
                break
            input_path = model_dir / f"frozen-profile-{index:02d}.inp"
            input_path.write_text(
                render_input(template, args.model, index), encoding="utf-8"
            )
            stdout_path = runtime / f"image-{index:02d}.stdout.log"
            cidfile = runtime / f"image-{index:02d}.cid"
            active_cidfile = cidfile
            command = [
                docker,
                "run",
                "--rm",
                "--cidfile",
                str(cidfile),
                "--cpus=16",
                f"--memory={args.memory_gib}g",
                "-e",
                f"OMP_NUM_THREADS={args.omp_threads}",
                "-e",
                "OPENBLAS_NUM_THREADS=1",
                "-e",
                "MKL_NUM_THREADS=1",
                "-e",
                "NUMEXPR_NUM_THREADS=1",
                "-v",
                f"{runtime}:{runtime}",
                "-w",
                str(model_dir),
                args.image,
                "mpirun",
                "-np",
                str(args.mpi_ranks),
                args.cp2k,
                "-i",
                input_path.name,
            ]
            image_started = time.monotonic()
            timed_out = False
            cleanup: dict[str, Any] = {}
            with stdout_path.open("wb") as output:
                active_process = subprocess.Popen(
                    command,
                    cwd=model_dir,
                    stdout=output,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                )
                try:
                    active_process.wait(
                        timeout=min(args.per_image_timeout, remaining - 30)
                    )
                except subprocess.TimeoutExpired:
                    timed_out = True
                    cleanup.update(_stop_process(active_process))
                finally:
                    output.flush()
                    os.fsync(output.fileno())
            returncode = active_process.returncode
            active_process = None
            cleanup.update(_remove_container(docker, cidfile))
            active_cidfile = None
            text = stdout_path.read_text(encoding="utf-8", errors="replace")
            classification = e3b.classify_cp2k_output(
                text, returncode=returncode, timed_out=timed_out
            )
            record = {
                "index": index,
                "classification": classification,
                "elapsed_seconds": time.monotonic() - image_started,
                "coordinate": _file_record(prepared_images[index], prepared),
                "input": _file_record(input_path, runtime),
                "stdout": _file_record(stdout_path, runtime),
                "command": command,
                "cleanup": cleanup,
            }
            image_records.append(record)
            if classification.get("status") != "converged":
                status = "incomplete-image"
                failure = f"image {index:02d} did not converge"
                break
    except Interrupted as exc:
        status = "incomplete-interrupted"
        failure = str(exc)
    finally:
        if active_process is not None:
            cleanup = _stop_process(active_process)
            if image_records:
                image_records[-1].setdefault("cleanup", {}).update(cleanup)
        if active_cidfile is not None:
            _remove_container(docker, active_cidfile)

    e3b.verify_manifest(prepared)
    energies = [item["classification"].get("energy_hartree") for item in image_records]
    complete = (
        status == "complete"
        and len(energies) == IMAGE_COUNT
        and all(
            isinstance(value, (int, float)) and math.isfinite(value)
            for value in energies
        )
    )
    profile_rise = None
    if complete:
        numeric = [float(value) for value in energies]
        profile_rise = (max(numeric) - numeric[0]) * HARTREE_TO_KCAL_MOL
        if not math.isfinite(profile_rise) or profile_rise < 0:
            complete = False
            status = "incomplete-numeric"
            failure = "profile rise is not finite and non-negative"
    else:
        status = status if status != "complete" else "incomplete-profile"

    runtime_files = [
        _file_record(path, runtime)
        for path in sorted(runtime.rglob("*"))
        if path.is_file() and path.resolve() != receipt_path
    ]
    receipt = {
        "schema": "e3b-frozen-path-profile-v1",
        "operator": "(hermes-custom-build-001; profile=workstation)",
        "purpose": (
            "survey PBE-D3 single-point energies on immutable classical path images; "
            "not endpoint optimization, DFT CI-NEB, or a production barrier"
        ),
        "status": "complete" if complete else status,
        "failure": None if complete else failure,
        "started_at": started_at.isoformat(),
        "completed_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "elapsed_seconds": time.monotonic() - started,
        "model": args.model,
        "atom_identity_sha256": model_record.get("atom_identity_sha256"),
        "prepared": {
            "path": str(prepared),
            "preparation_sha256": e3b.sha256(prepared / "preparation.json"),
            "manifest_sha256": e3b.sha256(prepared / "manifest.json"),
        },
        "method": {
            "basis": "frozen-classical-path-single-points",
            "image": args.image,
            "cp2k": args.cp2k,
            "image_count": IMAGE_COUNT,
            "hartree_to_kcal_mol": HARTREE_TO_KCAL_MOL,
        },
        "resources": {
            "mpi_ranks": args.mpi_ranks,
            "omp_threads": args.omp_threads,
            "total_threads": args.mpi_ranks * args.omp_threads,
            "memory_gib": args.memory_gib,
            "total_timeout_seconds": args.timeout_seconds,
            "per_image_timeout_seconds": args.per_image_timeout,
        },
        "controller": controller,
        "images": image_records,
        "profile_energies_hartree": [float(value) for value in energies]
        if complete
        else None,
        "profile_rise_kcal_mol": profile_rise if complete else None,
        "claim_boundary": {
            "supports": "energy rise along the tested immutable classical path",
            "does_not_support": [
                "DFT-relaxed endpoint energies",
                "DFT minimum-energy path",
                "DFT CI-NEB convergence",
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
    parser.add_argument("--image", required=True)
    parser.add_argument("--cp2k", default="/opt/cp2k/exe/local/cp2k.psmp")
    parser.add_argument("--docker", default="docker")
    parser.add_argument("--mpi-ranks", type=int, default=4)
    parser.add_argument("--omp-threads", type=int, default=4)
    parser.add_argument("--memory-gib", type=int, default=48)
    parser.add_argument("--timeout-seconds", type=float, required=True)
    parser.add_argument("--per-image-timeout", type=float, default=1800)
    parser.add_argument("--systemd-unit", required=True)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    receipt = run_profile(args)
    return 0 if receipt["status"] == "complete" else 2


if __name__ == "__main__":
    sys.exit(main())
