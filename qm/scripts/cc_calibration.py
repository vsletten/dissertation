#!/usr/bin/env python3
"""Open-source coupled-cluster calibration jobs for A2.

Run ``psi4`` jobs under conda env ``calib-psi4`` and ``byteqc`` jobs under
``~/venvs/byteqc``.  Each invocation handles one immutable geometry/basis and
writes one atomic JSON receipt; ``summarize`` performs the two-point CBS barrier
extrapolation and DLPNO-vs-canonical gate without importing either engine.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path
from typing import Any

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    _ETIQUETTE = bootstrap_cli(
        "cc-calibration",
        default_run_root="/mnt/data/vsletten/dissertation-data/task207-a2-production",
    )

HARTREE_TO_KJ = 2625.4996394798254
HF_CBS_ALPHA = 1.63


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def xyz_rows(path: Path) -> list[tuple[str, float, float, float]]:
    lines = path.read_text().splitlines()
    count = int(lines[0])
    if len(lines) != count + 2:
        raise ValueError(f"{path}: XYZ line count drift")
    rows = []
    for line in lines[2:]:
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"{path}: malformed XYZ row")
        rows.append((fields[0], float(fields[1]), float(fields[2]), float(fields[3])))
    return rows


def frozen_core_orbitals(symbols: list[str]) -> int:
    """Conventional frozen-core occupied orbitals for H--Ar."""
    from pyscf.data import elements

    count = 0
    for symbol in symbols:
        charge = int(elements.charge(symbol))
        if charge <= 2:
            continue
        if charge <= 10:
            count += 1
        elif charge <= 18:
            count += 5
        else:
            raise ValueError(f"no frozen-core convention encoded for {symbol}")
    return count


def existing_receipt(
    path: Path,
    *,
    engine: str,
    basis: str,
    input_sha256: str,
) -> dict[str, Any] | None:
    if not path.exists():
        return None
    payload = json.loads(path.read_text())
    expected = (engine, basis.lower(), input_sha256)
    actual = (
        payload.get("engine"),
        str(payload.get("basis", "")).lower(),
        payload.get("input_sha256"),
    )
    if actual != expected:
        raise ValueError(f"{path}: existing receipt identity drift")
    return payload


def psi4_job(args: argparse.Namespace) -> dict[str, Any]:
    import psi4

    source = args.xyz.resolve()
    source_sha = sha256_path(source)
    if cached := existing_receipt(
        args.output,
        engine="psi4-dlpno-ccsd(t)",
        basis=args.basis,
        input_sha256=source_sha,
    ):
        return cached
    rows = xyz_rows(source)
    multiplicity = args.spin + 1
    geometry = [f"{args.charge} {multiplicity}", "no_reorient", "no_com"]
    geometry.extend(f"{s} {x:.12f} {y:.12f} {z:.12f}" for s, x, y, z in rows)
    molecule = psi4.geometry("\n".join(geometry))
    psi4.set_num_threads(args.threads)
    psi4.set_memory(f"{args.memory_gb} GB")
    psi4.core.set_output_file(str(args.log.resolve()), False)
    psi4.set_options(
        {
            "basis": args.basis,
            "freeze_core": True,
            "scf_type": "df",
            "pno_convergence": "tight",
        }
    )
    started = time.time()
    total = float(psi4.energy("dlpno-ccsd(t)", molecule=molecule))
    elapsed = time.time() - started
    scf = float(psi4.core.variable("SCF TOTAL ENERGY"))
    payload = {
        "engine": "psi4-dlpno-ccsd(t)",
        "engine_version": psi4.__version__,
        "basis": args.basis,
        "pno_convergence": "tight",
        "freeze_core": True,
        "charge": args.charge,
        "spin_2s": args.spin,
        "input": str(source),
        "input_sha256": source_sha,
        "scf_hartree": scf,
        "correlation_hartree": total - scf,
        "total_hartree": total,
        "elapsed_seconds": elapsed,
        "threads": args.threads,
        "memory_gb": args.memory_gb,
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    atomic_json(args.output, payload)
    return payload


def byteqc_job(args: argparse.Namespace) -> dict[str, Any]:
    import cupy
    import pyscf
    from byteqc import cucc
    from byteqc.cucc.ccsd_t import kernel as triples_kernel
    from pyscf import gto, lib, scf

    source = args.xyz.resolve()
    source_sha = sha256_path(source)
    if cached := existing_receipt(
        args.output,
        engine="byteqc-canonical-ccsd(t)",
        basis=args.basis,
        input_sha256=source_sha,
    ):
        return cached
    rows = xyz_rows(source)
    symbols = [row[0] for row in rows]
    lib.num_threads(args.threads)
    molecule = gto.M(
        atom=[(symbol, (x, y, z)) for symbol, x, y, z in rows],
        basis=args.basis,
        charge=args.charge,
        spin=args.spin,
        verbose=4,
        output=str(args.log.resolve()),
        max_memory=args.memory_gb * 1024,
    )
    started = time.time()
    mean_field = scf.RHF(molecule).density_fit()
    mean_field.max_cycle = 150
    mean_field.conv_tol = 1e-10
    scf_energy = float(mean_field.kernel())
    if not mean_field.converged:
        raise RuntimeError("ByteQC calibration RHF did not converge")
    frozen = frozen_core_orbitals(symbols)
    coupled_cluster = cucc.CCSD(
        mean_field,
        frozen=frozen,
        gpulim=args.gpu_memory_gb << 30,
    )
    coupled_cluster.max_cycle = 100
    coupled_cluster.conv_tol = 1e-8
    ccsd_correlation, _, _ = coupled_cluster.kernel()
    if not coupled_cluster.converged:
        raise RuntimeError("ByteQC CCSD did not converge")
    eris = coupled_cluster.ao2mo()
    triples = float(triples_kernel(coupled_cluster, eris))
    ccsd_correlation = float(ccsd_correlation)
    total = scf_energy + ccsd_correlation + triples
    elapsed = time.time() - started
    payload = {
        "engine": "byteqc-canonical-ccsd(t)",
        "engine_version": "source checkout",
        "pyscf_version": pyscf.__version__,
        "cupy_version": cupy.__version__,
        "basis": args.basis,
        "freeze_core": True,
        "frozen_core_orbitals": frozen,
        "charge": args.charge,
        "spin_2s": args.spin,
        "input": str(source),
        "input_sha256": source_sha,
        "scf_hartree": scf_energy,
        "ccsd_correlation_hartree": ccsd_correlation,
        "triples_hartree": triples,
        "correlation_hartree": ccsd_correlation + triples,
        "total_hartree": total,
        "elapsed_seconds": elapsed,
        "threads": args.threads,
        "memory_gb": args.memory_gb,
        "gpu_memory_gb": args.gpu_memory_gb,
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    atomic_json(args.output, payload)
    return payload


def extrapolate_hf(tz: float, qz: float, alpha: float = HF_CBS_ALPHA) -> float:
    ratio = pow(2.718281828459045, -alpha)
    return (qz - ratio * tz) / (1.0 - ratio)


def extrapolate_correlation(tz: float, qz: float) -> float:
    return (4.0**3 * qz - 3.0**3 * tz) / (4.0**3 - 3.0**3)


def cbs_total(tz: dict[str, Any], qz: dict[str, Any]) -> float:
    return extrapolate_hf(
        float(tz["scf_hartree"]), float(qz["scf_hartree"])
    ) + extrapolate_correlation(
        float(tz["correlation_hartree"]), float(qz["correlation_hartree"])
    )


def load_receipt(path: Path, engine: str) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if payload.get("engine") != engine:
        raise ValueError(f"{path}: expected engine {engine}")
    if str(payload.get("basis", "")).lower() not in {"cc-pvtz", "cc-pvqz"}:
        raise ValueError(f"{path}: expected cc-pVTZ or cc-pVQZ")
    return payload


def summarize(args: argparse.Namespace) -> dict[str, Any]:
    canonical_engine = "byteqc-canonical-ccsd(t)"
    dlpno_engine = "psi4-dlpno-ccsd(t)"
    canonical = {
        role: {
            basis: load_receipt(path, canonical_engine) for basis, path in paths.items()
        }
        for role, paths in args.canonical.items()
    }
    dlpno = {
        role: {basis: load_receipt(path, dlpno_engine) for basis, path in paths.items()}
        for role, paths in args.dlpno.items()
    }
    canonical_barrier = (
        cbs_total(canonical["ts"]["tz"], canonical["ts"]["qz"])
        - cbs_total(canonical["reactant"]["tz"], canonical["reactant"]["qz"])
    ) * HARTREE_TO_KJ
    dlpno_barrier = (
        cbs_total(dlpno["ts"]["tz"], dlpno["ts"]["qz"])
        - cbs_total(dlpno["reactant"]["tz"], dlpno["reactant"]["qz"])
    ) * HARTREE_TO_KJ
    delta = dlpno_barrier - canonical_barrier
    payload = {
        "method": "frozen-core CCSD(T)/CBS, cc-pVTZ/cc-pVQZ",
        "hf_extrapolation_alpha": HF_CBS_ALPHA,
        "correlation_extrapolation_power": 3,
        "canonical_barrier_kj": canonical_barrier,
        "dlpno_barrier_kj": dlpno_barrier,
        "dlpno_minus_canonical_kj": delta,
        "gate_limit_kj": 2.0,
        "gate_pass": abs(delta) <= 2.0,
        "inputs": {
            "canonical": {
                role: {basis: str(path) for basis, path in paths.items()}
                for role, paths in args.canonical.items()
            },
            "dlpno": {
                role: {basis: str(path) for basis, path in paths.items()}
                for role, paths in args.dlpno.items()
            },
        },
    }
    atomic_json(args.output, payload)
    return payload


def receipt_grid(values: list[str]) -> dict[str, dict[str, Path]]:
    result: dict[str, dict[str, Path]] = {}
    for value in values:
        fields = value.split("=", 2)
        if len(fields) != 3:
            raise ValueError("receipt must be ROLE=BASIS=PATH")
        role, basis, path = fields
        normalized = basis.lower().replace("cc-pv", "").replace("z", "")
        if role not in {"reactant", "ts"} or normalized not in {"t", "q"}:
            raise ValueError(f"invalid receipt selector {value}")
        result.setdefault(role, {})["tz" if normalized == "t" else "qz"] = Path(path)
    if set(result) != {"reactant", "ts"} or any(
        set(paths) != {"tz", "qz"} for paths in result.values()
    ):
        raise ValueError("receipts must cover reactant/ts at TZ/QZ")
    return result


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--threads", type=int, default=16)
    result.add_argument("--nice", type=int, default=10)
    result.add_argument("--log")
    subparsers = result.add_subparsers(dest="command", required=True)
    for command in ("psi4", "byteqc"):
        job = subparsers.add_parser(command)
        job.add_argument("--xyz", type=Path, required=True)
        job.add_argument("--basis", choices=("cc-pVTZ", "cc-pVQZ"), required=True)
        job.add_argument("--charge", type=int, required=True)
        job.add_argument("--spin", type=int, default=0)
        job.add_argument("--memory-gb", type=int, default=48)
        job.add_argument("--output", type=Path, required=True)
        job.add_argument("--engine-log", dest="engine_log", type=Path, required=True)
        if command == "byteqc":
            job.add_argument("--gpu-memory-gb", type=int, default=16)
    summary = subparsers.add_parser("summarize")
    summary.add_argument("--canonical", action="append", default=[], required=True)
    summary.add_argument("--dlpno", action="append", default=[], required=True)
    summary.add_argument("--output", type=Path, required=True)
    return result


def main() -> int:
    args = parser().parse_args()
    if args.command in {"psi4", "byteqc"}:
        args.log = args.engine_log
        if args.threads < 1 or args.threads > 16:
            raise ValueError("threads must be in 1..16")
        if args.memory_gb < 1:
            raise ValueError("memory-gb must be positive")
        payload = psi4_job(args) if args.command == "psi4" else byteqc_job(args)
    else:
        args.canonical = receipt_grid(args.canonical)
        args.dlpno = receipt_grid(args.dlpno)
        payload = summarize(args)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
