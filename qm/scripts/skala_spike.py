#!/usr/bin/env python3
"""Run the bounded A2e Skala single-point comparison.

The source A2a store is opened only through an immutable read-only SQLite URI.
All new jobs and results are written to a fresh run directory.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import sqlite3
import time
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    _ETIQUETTE = bootstrap_cli(
        "a2e-skala-spike",
        default_run_root="/mnt/data/vsletten/dissertation-data/a2e-skala-functional-spike",
        gpu_owner="a2e_skala_functional_spike",
        gpu_ttl_hours=4.0,
    )

from quarry.pipeline import HARTREE_TO_KJ
from quarry.store import Store, geometry_hash

SOURCE_ROOT = Path(
    "/mnt/data/vsletten/dissertation-data/"
    "task208-a2a-path-rebuild-20260825/production-closeout"
)
SOURCE_STORE_SHA256 = "480cc244cf06bf7fd8edce86e4d57c3a8466e24e2793daf39be3fe27207dcd89"
FOCAL_BARRIER_KJ_MOL = 132.960133
EXPECTED_STRUCTURES = {
    1: (
        "reactant",
        "si-neutral-reactant",
        "f976b1f07cfe0db48a593854fc5d9dd794b24a243d005c5c79c5368f9c635d40",
    ),
    2: (
        "intermediate",
        "si-neutral-intermediate",
        "76ed3d0652cc93e6b024f7649102e7b53711d8054a0f0e6a3296b5540097df1e",
    ),
    3: (
        "addition-transition-state",
        "si-neutral-addition-transition-state",
        "b94fce89da1f03976e877447f154f527e7da73bd19659b7fc96a70733a849e92",
    ),
    4: (
        "released-product",
        "si-neutral-released-product",
        "dbbf58a07a595c0adfb97a3a36d29ca925b83abcf8b3a388dc5b427d3f9708e2",
    ),
}
METHODS = {
    "skala-1.1-def2-tzvp": {
        "method": "skala-1.1/def2-tzvp/gas/no-dispersion",
        "family": "skala",
        "basis": "def2-tzvp",
        "d4": False,
    },
    "skala-1.1-def2-tzvpd": {
        "method": "skala-1.1/def2-tzvpd/gas/no-dispersion",
        "family": "skala",
        "basis": "def2-tzvpd",
        "d4": False,
    },
    "skala-1.1-def2-tzvp-plus-d4-r2scan-damping": {
        "method": (
            "skala-1.1/def2-tzvp/gas + "
            "D4[r2scan-damping sensitivity; not Skala-parameterized]"
        ),
        "family": "skala",
        "basis": "def2-tzvp",
        "d4": True,
    },
    "wb97m-v-def2-tzvpd": {
        "method": "wb97m-v/def2-tzvpd/gas/native-vv10",
        "family": "wb97m-v",
        "basis": "def2-tzvpd",
        "d4": False,
    },
}
GAS_COMPARATORS_KJ_MOL = {
    "canonical-TZ-plus-TightPNO-CBS-focal": 132.960133,
    "canonical-CCSD(T)/cc-pVTZ": 128.085164,
    "TightPNO-DLPNO-CCSD(T)/cc-pVTZ": 128.330864,
    "r2SCAN-3c": 125.742278,
    "B3LYP/def2-SVP/DF": 114.173360,
}
SMD_CONTEXT_KJ_MOL = {
    "wb97m-v/def2-tzvpd/SMD(water)": 134.504248,
    "B3LYP-D4/def2-tzvpd/SMD(water)": 148.293616,
}


@dataclass(frozen=True)
class SourceStructure:
    source_id: int
    role: str
    name: str
    formula: str
    charge: int
    spin: int
    xyz: str
    geometry_hash: str


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    encoded = json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    try:
        with temporary.open("w", encoding="utf-8") as stream:
            stream.write(encoded)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def utc_now() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


def package_versions() -> dict[str, str]:
    versions = {}
    for distribution in (
        "skala",
        "skala-cuda12x",
        "torch",
        "pyscf",
        "gpu4pyscf-cuda12x",
        "pyscf-dispersion",
        "cupy-cuda12x",
    ):
        try:
            versions[distribution] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            versions[distribution] = "not-installed"
    return versions


def load_source_structures(source_store: Path) -> list[SourceStructure]:
    actual = sha256_path(source_store)
    if actual != SOURCE_STORE_SHA256:
        raise RuntimeError(
            f"A2a store SHA-256 drift: expected {SOURCE_STORE_SHA256}, found {actual}"
        )
    connection = sqlite3.connect(
        f"{source_store.resolve().as_uri()}?mode=ro&immutable=1", uri=True
    )
    connection.row_factory = sqlite3.Row
    try:
        rows = connection.execute(
            "SELECT id, name, formula, charge, spin, xyz, geometry_hash "
            "FROM structures WHERE id IN (1, 2, 3, 4) ORDER BY id"
        ).fetchall()
    finally:
        connection.close()
    if [int(row["id"]) for row in rows] != [1, 2, 3, 4]:
        raise RuntimeError("A2a store does not contain exact structure ids 1-4")
    structures = []
    for row in rows:
        source_id = int(row["id"])
        role, expected_name, expected_hash = EXPECTED_STRUCTURES[source_id]
        observed_hash = str(row["geometry_hash"])
        xyz = str(row["xyz"])
        if (
            str(row["name"]) != expected_name
            or observed_hash != expected_hash
            or geometry_hash(xyz) != expected_hash
            or str(row["formula"]) != "H8O8Si2"
            or int(row["charge"]) != 0
            or int(row["spin"]) != 0
        ):
            raise RuntimeError(f"A2a structure {source_id} identity drift")
        structures.append(
            SourceStructure(
                source_id=source_id,
                role=role,
                name=expected_name,
                formula=str(row["formula"]),
                charge=int(row["charge"]),
                spin=int(row["spin"]),
                xyz=xyz,
                geometry_hash=observed_hash,
            )
        )
    return structures


def _atom_spec(xyz: str) -> str:
    lines = xyz.strip().splitlines()
    if len(lines) < 3:
        raise ValueError("malformed XYZ: missing header")
    try:
        atom_count = int(lines[0])
    except ValueError as exc:
        raise ValueError("malformed XYZ atom count") from exc
    atom_lines = lines[2:]
    if len(atom_lines) != atom_count:
        raise ValueError("malformed XYZ atom count/body mismatch")
    for line in atom_lines:
        if len(line.split()) != 4:
            raise ValueError("malformed XYZ atom row")
    return "; ".join(atom_lines)


def _make_mol(structure: SourceStructure, basis: str) -> Any:
    from pyscf import gto

    return gto.M(
        atom=_atom_spec(structure.xyz),
        unit="Angstrom",
        basis=basis,
        charge=structure.charge,
        spin=structure.spin,
        verbose=0,
    )


def _finite(value: Any) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise RuntimeError(f"non-finite electronic energy: {number}")
    return number


def run_skala_scf(
    structure: SourceStructure,
    *,
    basis: str,
    use_gpu: bool,
) -> tuple[float, bool, int | None, dict[str, Any]]:
    mol = _make_mol(structure, basis)
    config = {"verbose": 0, "max_cycle": 150}
    if use_gpu:
        from skala.gpu4pyscf import SkalaKS

        mean_field = SkalaKS(
            mol,
            xc="skala-1.1",
            with_density_fit=False,
            with_dftd3=False,
            ks_config=config,
        )
    else:
        import torch
        from skala.pyscf import SkalaKS

        mean_field = SkalaKS(
            mol,
            xc="skala-1.1",
            with_density_fit=False,
            with_dftd3=False,
            ks_config=config,
            device=torch.device("cpu"),
        )
    value = _finite(mean_field.kernel())
    cycles = getattr(mean_field, "cycles", None)
    if cycles is not None:
        cycles = int(cycles)
    return (
        value,
        bool(mean_field.converged),
        cycles,
        {
            "density_fit": False,
            "dftd3": False,
            "grid": "SkalaKS default",
            "max_cycle": 150,
        },
    )


def run_wb97mv_scf(
    structure: SourceStructure,
    *,
    basis: str,
    use_gpu: bool,
) -> tuple[float, bool, int | None, dict[str, Any]]:
    from pyscf import dft

    mol = _make_mol(structure, basis)
    mean_field = dft.RKS(mol)
    mean_field.xc = "wb97m-v"
    mean_field.max_cycle = 150
    mean_field.verbose = 0
    if use_gpu:
        mean_field = mean_field.to_gpu()
    value = _finite(mean_field.kernel())
    cycles = getattr(mean_field, "cycles", None)
    if cycles is not None:
        cycles = int(cycles)
    return (
        value,
        bool(mean_field.converged),
        cycles,
        {
            "density_fit": False,
            "dispersion": "native VV10 in wb97m-v",
            "grid": "PySCF default",
            "max_cycle": 150,
        },
    )


def r2scan_d4_sensitivity(structure: SourceStructure, basis: str) -> float:
    """Return a separately labelled D4 sensitivity correction.

    DFT-D4 has no Skala damping parameters. The requested +D4 row therefore
    uses the library's r2SCAN damping as an explicit sensitivity only; it is
    never represented as a parameterized Skala-D4 method.
    """
    from pyscf.dispersion.dftd4 import DFTD4Dispersion

    mol = _make_mol(structure, basis)
    correction = DFTD4Dispersion(mol, xc="r2scan", atm=True).get_dispersion()
    return _finite(correction["energy"])


def classify_verdict(barrier_kj_mol: float) -> dict[str, Any]:
    delta = abs(float(barrier_kj_mol) - FOCAL_BARRIER_KJ_MOL)
    if delta <= 4.2:
        verdict = "adopt as a survey-tier single-point functional candidate"
    elif delta <= 8.4:
        verdict = "usable with stated uncertainty"
    else:
        verdict = "reject for Si-O-Si hydrolysis"
    return {
        "absolute_error_kj_mol": delta,
        "focal_barrier_kj_mol": FOCAL_BARRIER_KJ_MOL,
        "verdict": verdict,
    }


def derive_rows(jobs: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    for key, method in METHODS.items():
        by_role = {
            str(job["role"]): job
            for job in jobs
            if job["method_key"] == key and job["status"] == "done"
        }
        failures = [
            job for job in jobs if job["method_key"] == key and job["status"] != "done"
        ]
        row: dict[str, Any] = {
            "method": method["method"],
            "status": "complete"
            if len(by_role) == 4 and not failures
            else "incomplete",
            "failed_roles": [job["role"] for job in failures],
        }
        if row["status"] == "complete":
            reference = float(by_role["reactant"]["energy_hartree"])
            row["delta_kj_mol"] = {
                role: (float(record["energy_hartree"]) - reference) * HARTREE_TO_KJ
                for role, record in by_role.items()
                if role != "reactant"
            }
            barrier = row["delta_kj_mol"]["addition-transition-state"]
            row["barrier_kj_mol"] = barrier
            row["comparison_to_focal"] = classify_verdict(barrier)
        rows[key] = row
    return rows


def _record_job(
    store: Store,
    structure: SourceStructure,
    structure_id: int,
    method_key: str,
    method: dict[str, Any],
    *,
    use_gpu: bool,
    versions: dict[str, str],
    receipt_dir: Path,
) -> dict[str, Any]:
    engine = "gpu4pyscf" if use_gpu else "pyscf"
    detail = {
        "basis": method["basis"],
        "device": "cuda:0" if use_gpu else "cpu",
        "dispersion_mode": (
            "D4 r2SCAN damping sensitivity" if method["d4"] else "none/native"
        ),
        "geometry_hash": structure.geometry_hash,
        "package_versions": versions,
    }
    job_id = store.add_job(
        structure_id,
        "sp",
        str(method["method"]),
        engine,
        detail=json.dumps(detail, separators=(",", ":"), sort_keys=True),
    )
    store.set_job_status(job_id, "running")
    started = utc_now()
    started_clock = time.monotonic()
    receipt: dict[str, Any] = {
        "schema": "a2e-skala-single-point-v1",
        "job_id": job_id,
        "method_key": method_key,
        "method": method["method"],
        "role": structure.role,
        "source_structure_id": structure.source_id,
        "geometry_hash": structure.geometry_hash,
        "device": detail["device"],
        "package_versions": versions,
        "started_at": started,
    }
    try:
        if method["family"] == "skala":
            scf_energy, converged, cycles, settings = run_skala_scf(
                structure, basis=str(method["basis"]), use_gpu=use_gpu
            )
        else:
            scf_energy, converged, cycles, settings = run_wb97mv_scf(
                structure, basis=str(method["basis"]), use_gpu=use_gpu
            )
        if not converged:
            raise RuntimeError("SCF did not converge")
        dispersion = 0.0
        if bool(method["d4"]):
            dispersion = r2scan_d4_sensitivity(structure, basis=str(method["basis"]))
        energy = _finite(scf_energy + dispersion)
        wall = time.monotonic() - started_clock
        receipt.update(
            {
                "status": "done",
                "converged": True,
                "cycles": cycles,
                "wall_seconds": wall,
                "electronic_hartree": energy,
                "scf_hartree": scf_energy,
                "dispersion_hartree": dispersion,
                "settings": settings,
                "finished_at": utc_now(),
            }
        )
        store.add_result(job_id, "energy", energy, "hartree")
        store.add_result(job_id, "scf_energy", scf_energy, "hartree")
        store.add_result(job_id, "dispersion", dispersion, "hartree")
        store.add_result(job_id, "scf_converged", 1.0, "boolean")
        store.add_result(job_id, "wall_time", wall, "seconds")
        if cycles is not None:
            store.add_result(job_id, "scf_cycles", float(cycles), "count")
        store.set_job_status(job_id, "done")
    except BaseException as exc:
        receipt.update(
            {
                "status": "failed",
                "converged": False,
                "wall_seconds": time.monotonic() - started_clock,
                "error_type": type(exc).__name__,
                "error": str(exc),
                "finished_at": utc_now(),
            }
        )
        store.set_job_status(job_id, "failed")
    receipt_path = receipt_dir / f"{method_key}--{structure.role}.json"
    atomic_json(receipt_path, receipt)
    receipt["receipt_path"] = str(receipt_path)
    receipt["receipt_sha256"] = sha256_path(receipt_path)
    store.conn.execute(
        "UPDATE jobs SET detail = ? WHERE id = ?",
        (
            json.dumps(
                {**detail, "receipt_sha256": receipt["receipt_sha256"]},
                separators=(",", ":"),
                sort_keys=True,
            ),
            job_id,
        ),
    )
    store.conn.commit()
    print(
        f"job {job_id:02d} {method_key} {structure.role}: {receipt['status']} "
        f"({receipt['wall_seconds']:.3f}s)",
        flush=True,
    )
    return receipt


def run(args: argparse.Namespace) -> int:
    run_dir = args.run_dir.resolve()
    run_dir.mkdir(parents=True, exist_ok=False)
    source_store = args.source_root.resolve() / "store.sqlite"
    source_before = sha256_path(source_store)
    structures = load_source_structures(source_store)
    versions = package_versions()
    jobs: list[dict[str, Any]] = []
    store_path = run_dir / "store.sqlite"
    receipt_dir = run_dir / "receipts"
    with Store(store_path) as store:
        structure_ids = {
            structure.role: store.add_structure(
                structure.name,
                structure.formula,
                structure.xyz,
                charge=structure.charge,
                spin=structure.spin,
            )
            for structure in structures
        }
        for method_key, method in METHODS.items():
            for structure in structures:
                jobs.append(
                    _record_job(
                        store,
                        structure,
                        structure_ids[structure.role],
                        method_key,
                        method,
                        use_gpu=args.gpu,
                        versions=versions,
                        receipt_dir=receipt_dir,
                    )
                )
    source_after = sha256_path(source_store)
    if source_before != SOURCE_STORE_SHA256 or source_after != source_before:
        raise RuntimeError("A2a evidence store changed during the spike")
    results = {
        "schema": "a2e-skala-functional-spike-v1",
        "generated_at": utc_now(),
        "wall_seconds": sum(float(job["wall_seconds"]) for job in jobs),
        "execution": {
            "device": "cuda:0" if args.gpu else "cpu",
            "hostname": platform.node(),
            "python": platform.python_version(),
            "thread_environment": {
                name: os.environ.get(name)
                for name in (
                    "OMP_NUM_THREADS",
                    "MKL_NUM_THREADS",
                    "OPENBLAS_NUM_THREADS",
                )
            },
            "package_versions": versions,
        },
        "source": {
            "store_path": str(source_store),
            "store_sha256_before": source_before,
            "store_sha256_after": source_after,
            "immutable_sqlite_uri": True,
        },
        "jobs": jobs,
        "rows": derive_rows(jobs),
        "gas_comparators_kj_mol": GAS_COMPARATORS_KJ_MOL,
        "smd_context_not_gas_comparable_kj_mol": SMD_CONTEXT_KJ_MOL,
        "limits": {
            "one_reaction_is_not_validation": True,
            "geometry_optimization": False,
            "frequencies": False,
            "solvent": False,
            "coupled_cluster": False,
            "d4_sensitivity": (
                "DFT-D4 has no Skala damping parameters; the separately labelled "
                "+D4 row uses r2SCAN damping only as a sensitivity and is not a "
                "parameterized Skala-D4 method."
            ),
        },
    }
    results_path = run_dir / "results.json"
    atomic_json(results_path, results)
    manifest = {
        "schema": "a2e-skala-artifact-manifest-v1",
        "generated_at": utc_now(),
        "results_json": {
            "path": str(results_path),
            "sha256": sha256_path(results_path),
        },
        "store_sqlite": {
            "path": str(store_path),
            "sha256": sha256_path(store_path),
        },
        "receipt_count": len(jobs),
        "receipt_hashes": {
            Path(str(job["receipt_path"])).name: job["receipt_sha256"] for job in jobs
        },
        "source_store_sha256_before": source_before,
        "source_store_sha256_after": source_after,
    }
    atomic_json(run_dir / "manifest.json", manifest)
    failed = [job for job in jobs if job["status"] != "done"]
    print(f"results: {results_path}", flush=True)
    print(f"completed jobs: {len(jobs) - len(failed)}/{len(jobs)}", flush=True)
    return 0 if not failed else 2


def verify(args: argparse.Namespace) -> int:
    run_dir = args.run_dir.resolve()
    results_path = run_dir / "results.json"
    store_path = run_dir / "store.sqlite"
    manifest_path = run_dir / "manifest.json"
    payload = json.loads(results_path.read_text())
    manifest = json.loads(manifest_path.read_text())
    if payload.get("schema") != "a2e-skala-functional-spike-v1":
        raise RuntimeError("results schema drift")
    if len(payload.get("jobs", [])) != 16:
        raise RuntimeError("expected exactly 16 job receipts")
    if set(payload.get("rows", {})) != set(METHODS):
        raise RuntimeError("method row set drift")
    for job in payload["jobs"]:
        receipt_path = Path(str(job["receipt_path"]))
        if sha256_path(receipt_path) != job["receipt_sha256"]:
            raise RuntimeError(f"receipt hash drift: {receipt_path}")
    if manifest["results_json"]["sha256"] != sha256_path(results_path):
        raise RuntimeError("results.json hash drift")
    if manifest["store_sqlite"]["sha256"] != sha256_path(store_path):
        raise RuntimeError("store.sqlite hash drift")
    if sha256_path(args.source_root.resolve() / "store.sqlite") != SOURCE_STORE_SHA256:
        raise RuntimeError("A2a source store drift after run")
    connection = sqlite3.connect(f"{store_path.as_uri()}?mode=ro&immutable=1", uri=True)
    try:
        structure_count = connection.execute(
            "SELECT COUNT(*) FROM structures"
        ).fetchone()[0]
        job_count = connection.execute("SELECT COUNT(*) FROM jobs").fetchone()[0]
        done_count = connection.execute(
            "SELECT COUNT(*) FROM jobs WHERE status = 'done'"
        ).fetchone()[0]
        result_count = connection.execute("SELECT COUNT(*) FROM results").fetchone()[0]
    finally:
        connection.close()
    if structure_count != 4 or job_count != 16:
        raise RuntimeError("store row counts drift")
    recomputed = derive_rows(payload["jobs"])
    if recomputed != payload["rows"]:
        raise RuntimeError("derived comparison rows drift")
    print(
        json.dumps(
            {
                "verified": True,
                "structures": structure_count,
                "jobs": job_count,
                "done_jobs": done_count,
                "results": result_count,
                "results_sha256": sha256_path(results_path),
                "store_sha256": sha256_path(store_path),
            },
            sort_keys=True,
        )
    )
    return 0


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    root.add_argument("--source-root", type=Path, default=SOURCE_ROOT)
    # Accepted here because bootstrap_cli pre-parses and enforces them before
    # heavy imports; the full parser must still consume the same argv.
    root.add_argument("--threads", type=int, default=16)
    root.add_argument("--nice", type=int, default=10)
    root.add_argument("--log")
    root.add_argument("--gpu-mem-gb", type=float, default=16.0)
    subparsers = root.add_subparsers(dest="command", required=True)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--run-dir", type=Path, required=True)
    run_parser.add_argument("--gpu", action="store_true")
    run_parser.set_defaults(handler=run)
    verify_parser = subparsers.add_parser("verify")
    verify_parser.add_argument("--run-dir", type=Path, required=True)
    verify_parser.set_defaults(handler=verify)
    return root


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    return int(args.handler(args))


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    finally:
        _ETIQUETTE.close()
