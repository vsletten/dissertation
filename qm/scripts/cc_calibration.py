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
import subprocess
import time
from pathlib import Path
from typing import Any

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    _ETIQUETTE = bootstrap_cli(
        "cc-calibration",
        default_run_root="/mnt/data/vsletten/dissertation-data/task207-a2-production",
        gpu_owner="cc_calibration",
    )

import numpy as np

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


def clean_git_head(path: Path) -> str:
    status = subprocess.run(
        ["git", "-C", str(path), "status", "--porcelain", "--untracked-files=no"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if status:
        raise RuntimeError(f"engine source checkout has tracked changes: {path}")
    return subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


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


def configure_byteqc_frozen_core_blocking(coupled_cluster: Any) -> int:
    """Avoid ByteQC 2.5's undersized reused DF buffer with frozen orbitals."""
    naux = int(coupled_cluster.with_df.get_naoaux())
    coupled_cluster.with_df.blockdim = naux
    return naux


def existing_receipt(
    path: Path,
    *,
    engine: str,
    basis: str,
    input_sha256: str,
    identity: dict[str, Any],
) -> dict[str, Any] | None:
    if not path.exists():
        return None
    payload = json.loads(path.read_text())
    expected = (engine, basis.lower(), input_sha256, identity)
    actual = (
        payload.get("engine"),
        str(payload.get("basis", "")).lower(),
        payload.get("input_sha256"),
        payload.get("identity"),
    )
    if actual != expected:
        raise ValueError(f"{path}: existing receipt identity drift")
    return payload


def psi4_job(args: argparse.Namespace) -> dict[str, Any]:
    if args.spin != 0:
        raise ValueError("Psi4 DLPNO calibration currently supports spin=0 only")

    import psi4

    source = args.xyz.resolve()
    source_sha = sha256_path(source)
    identity = {
        "engine_version": psi4.__version__,
        "method": "dlpno-ccsd(t)",
        "basis": args.basis,
        "pno_convergence": "tight",
        "freeze_core": True,
        "scf_type": "df",
        "charge": args.charge,
        "spin_2s": args.spin,
    }
    if cached := existing_receipt(
        args.output,
        engine="psi4-dlpno-ccsd(t)",
        basis=args.basis,
        input_sha256=source_sha,
        identity=identity,
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
        "identity": identity,
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
    if args.spin != 0:
        raise ValueError(
            "ByteQC calibration currently supports closed-shell spin=0 only"
        )

    import byteqc
    import cupy
    import pyscf
    from byteqc import cucc
    from byteqc.cucc.ccsd_t import kernel as triples_kernel
    from pyscf import gto, lib, scf

    source = args.xyz.resolve()
    source_sha = sha256_path(source)
    rows = xyz_rows(source)
    symbols = [row[0] for row in rows]
    frozen = frozen_core_orbitals(symbols)
    if byteqc.__file__ is None:
        raise RuntimeError("ByteQC module has no source path")
    byteqc_root = Path(byteqc.__file__).resolve().parent
    byteqc_commit = clean_git_head(byteqc_root)
    identity = {
        "engine_commit": byteqc_commit,
        "pyscf_version": pyscf.__version__,
        "method": "canonical-ccsd(t)",
        "basis": args.basis,
        "freeze_core_orbitals": frozen,
        "density_fit": True,
        "df_aux_blocking": "single-block-frozen-core-v1",
        "charge": args.charge,
        "spin_2s": args.spin,
    }
    if cached := existing_receipt(
        args.output,
        engine="byteqc-canonical-ccsd(t)",
        basis=args.basis,
        input_sha256=source_sha,
        identity=identity,
    ):
        return cached
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
    coupled_cluster = cucc.CCSD(
        mean_field,
        frozen=frozen,
        gpulim=args.gpu_memory_gb << 30,
    )
    # ByteQC 2.5 reuses the MO-sized transform output as the next AO-sized
    # scratch buffer. That buffer is too small whenever frozen core makes
    # nmo < nao. One complete auxiliary block avoids the invalid reuse while
    # preserving the exact frozen-core canonical calculation.
    configure_byteqc_frozen_core_blocking(coupled_cluster)
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
        "engine_version": byteqc_commit,
        "identity": identity,
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
    ratio = float(np.exp(-alpha))
    return (qz - ratio * tz) / (1.0 - ratio)


def extrapolate_correlation(tz: float, qz: float) -> float:
    return (4.0**3 * qz - 3.0**3 * tz) / (4.0**3 - 3.0**3)


def cbs_total(tz: dict[str, Any], qz: dict[str, Any]) -> float:
    return extrapolate_hf(
        float(tz["scf_hartree"]), float(qz["scf_hartree"])
    ) + extrapolate_correlation(
        float(tz["correlation_hartree"]), float(qz["correlation_hartree"])
    )


def load_receipt(
    path: Path,
    engine: str,
    expected_basis: str,
) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if payload.get("engine") != engine:
        raise ValueError(f"{path}: expected engine {engine}")
    if str(payload.get("basis", "")).lower() != expected_basis.lower():
        raise ValueError(f"{path}: expected basis {expected_basis}")
    if not isinstance(payload.get("identity"), dict):
        raise ValueError(f"{path}: missing engine/settings identity")
    return payload


def summarize(args: argparse.Namespace) -> dict[str, Any]:
    canonical_engine = "byteqc-canonical-ccsd(t)"
    dlpno_engine = "psi4-dlpno-ccsd(t)"
    basis_names = {"tz": "cc-pVTZ", "qz": "cc-pVQZ"}
    canonical = {
        role: {
            basis: load_receipt(path, canonical_engine, basis_names[basis])
            for basis, path in paths.items()
        }
        for role, paths in args.canonical.items()
    }
    dlpno = {
        role: {
            basis: load_receipt(path, dlpno_engine, basis_names[basis])
            for basis, path in paths.items()
        }
        for role, paths in args.dlpno.items()
    }
    for role in ("reactant", "ts"):
        receipts = [*canonical[role].values(), *dlpno[role].values()]
        input_hashes = {receipt.get("input_sha256") for receipt in receipts}
        states = {
            (receipt.get("charge"), receipt.get("spin_2s")) for receipt in receipts
        }
        if None in input_hashes or len(input_hashes) != 1:
            raise ValueError(f"{role} receipts do not share one exact geometry")
        if len(states) != 1:
            raise ValueError(f"{role} receipts do not share one electronic state")
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
    result.add_argument("--gpu", action="store_true")
    result.add_argument("--gpu-mem-gb", type=float, default=16.0)
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


def validate_resource_contract(args: argparse.Namespace) -> None:
    if args.command == "byteqc":
        if not args.gpu:
            raise ValueError(
                "ByteQC calibration requires --gpu and the shared GPU lease"
            )
        if float(args.gpu_mem_gb) < float(args.gpu_memory_gb):
            raise ValueError("--gpu-mem-gb must be at least --gpu-memory-gb")
    elif args.gpu:
        raise ValueError("--gpu is valid only for ByteQC calibration jobs")


def main() -> int:
    args = parser().parse_args()
    validate_resource_contract(args)
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
