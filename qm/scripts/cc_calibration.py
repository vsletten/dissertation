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
import os
import shutil
import subprocess
import time
from contextlib import contextmanager
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
CCSD_RESTART_SCHEMA = "byteqc-ccsd-amplitudes-v1"
CCSD_CONVERGED_SCHEMA = "byteqc-ccsd-converged-v1"


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{time.time_ns()}.tmp")
    try:
        with temporary.open("w") as stream:
            stream.write(json.dumps(payload, indent=2, sort_keys=True) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        directory = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        temporary.unlink(missing_ok=True)


def canonical_json_sha256(payload: dict[str, Any]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def ccsd_converged_marker_path(checkpoint_path: Path) -> Path:
    return checkpoint_path.with_suffix(checkpoint_path.suffix + ".converged.json")


def write_ccsd_converged_marker(
    checkpoint_path: Path,
    checkpoint_metadata: dict[str, Any],
) -> dict[str, Any]:
    marker = {
        "schema": CCSD_CONVERGED_SCHEMA,
        "checkpoint_manifest_sha256": canonical_json_sha256(checkpoint_metadata),
        "identity_sha256": checkpoint_metadata["identity_sha256"],
        "completed_cycle": checkpoint_metadata["completed_cycle"],
        "ccsd_correlation_hartree": checkpoint_metadata["ccsd_correlation_hartree"],
        "written_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    atomic_json(ccsd_converged_marker_path(checkpoint_path), marker)
    return marker


def load_ccsd_converged_marker(
    checkpoint_path: Path,
    checkpoint_metadata: dict[str, Any],
) -> dict[str, Any] | None:
    path = ccsd_converged_marker_path(checkpoint_path)
    if not path.exists():
        return None
    try:
        marker = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError(f"{path}: malformed CCSD convergence marker") from exc
    expected = (
        CCSD_CONVERGED_SCHEMA,
        canonical_json_sha256(checkpoint_metadata),
        checkpoint_metadata["identity_sha256"],
        checkpoint_metadata["completed_cycle"],
        checkpoint_metadata["ccsd_correlation_hartree"],
    )
    actual = (
        marker.get("schema"),
        marker.get("checkpoint_manifest_sha256"),
        marker.get("identity_sha256"),
        marker.get("completed_cycle"),
        marker.get("ccsd_correlation_hartree"),
    )
    if actual != expected:
        raise ValueError(f"{path}: CCSD convergence marker identity drift")
    return marker


def numpy_sha256(value: Any) -> str:
    array = np.ascontiguousarray(value)
    digest = hashlib.sha256()
    digest.update(str(array.dtype).encode())
    digest.update(json.dumps(array.shape).encode())
    digest.update(memoryview(array).cast("B"))
    return digest.hexdigest()


def _numpy_chunk(value: Any, start: int, stop: int) -> np.ndarray:
    chunk = value[start:stop]
    if hasattr(chunk, "asnumpy"):
        chunk = chunk.asnumpy()
    elif hasattr(chunk, "get"):
        chunk = chunk.get()
    return np.ascontiguousarray(chunk)


def _checkpoint_dataset(store: Any, name: str, value: Any) -> str:
    shape = tuple(int(item) for item in value.shape)
    if not shape or shape[0] < 1:
        raise ValueError(f"{name}: amplitude tensor must be non-empty")
    dtype = np.dtype(value.dtype)
    dataset = store.create_dataset(
        name, shape=shape, dtype=dtype, chunks=(1, *shape[1:])
    )
    digest = hashlib.sha256()
    digest.update(str(dtype).encode())
    digest.update(json.dumps(shape).encode())
    for start in range(shape[0]):
        chunk = _numpy_chunk(value, start, start + 1)
        if chunk.shape != (1, *shape[1:]):
            raise ValueError(f"{name}: amplitude chunk shape drift")
        if not np.isfinite(chunk).all():
            raise ValueError(f"{name}: amplitude checkpoint contains non-finite values")
        digest.update(memoryview(chunk).cast("B"))
        dataset[start : start + 1] = chunk
    return digest.hexdigest()


def write_ccsd_checkpoint(
    path: Path,
    *,
    identity: dict[str, Any],
    completed_cycle: int,
    ccsd_correlation_hartree: float,
    t1: Any,
    t2: Any,
) -> dict[str, Any]:
    """Atomically persist one fully evaluated ByteQC CCSD amplitude state."""
    import h5py

    if completed_cycle < 1:
        raise ValueError("CCSD checkpoint cycle must be positive")
    energy = float(ccsd_correlation_hartree)
    if not np.isfinite(energy):
        raise ValueError("CCSD checkpoint energy must be finite")
    path.parent.mkdir(parents=True, exist_ok=True)
    checkpoint_bytes = int(t1.nbytes) + int(t2.nbytes)
    required_free = checkpoint_bytes * 2 + (1 << 30)
    available_free = shutil.disk_usage(path.parent).free
    if available_free < required_free:
        raise RuntimeError(
            "insufficient disk for atomic CCSD checkpoint: "
            f"need {required_free} bytes, have {available_free}"
        )
    temporary = path.with_name(f".{path.name}.{time.time_ns()}.tmp")
    try:
        with h5py.File(temporary, "w") as store:
            tensors = {}
            for name, value in (("t1", t1), ("t2", t2)):
                tensors[name] = {
                    "shape": list(value.shape),
                    "dtype": str(np.dtype(value.dtype)),
                    "sha256": _checkpoint_dataset(store, name, value),
                }
            metadata = {
                "schema": CCSD_RESTART_SCHEMA,
                "identity": identity,
                "identity_sha256": canonical_json_sha256(identity),
                "completed_cycle": completed_cycle,
                "ccsd_correlation_hartree": energy,
                "checkpoint_bytes": checkpoint_bytes,
                "tensors": tensors,
                "written_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            }
            store.attrs["manifest_json"] = json.dumps(metadata, sort_keys=True)
            store.attrs["manifest_sha256"] = canonical_json_sha256(metadata)
            store.flush()
        descriptor = os.open(temporary, os.O_RDONLY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
        os.replace(temporary, path)
        ccsd_converged_marker_path(path).unlink(missing_ok=True)
        directory = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        temporary.unlink(missing_ok=True)
    return metadata


def load_ccsd_checkpoint(
    path: Path,
    *,
    identity: dict[str, Any],
    coupled_cluster: Any,
) -> tuple[dict[str, Any], Any, Any] | None:
    """Load a validated restart directly into ByteQC's selected buffer backends."""
    import h5py

    if not path.exists():
        return None
    with h5py.File(path, "r") as store:
        try:
            metadata_raw = store.attrs["manifest_json"]
            if isinstance(metadata_raw, bytes):
                metadata_raw = metadata_raw.decode()
            if not isinstance(metadata_raw, str):
                raise TypeError("manifest_json is not text")
            metadata = json.loads(metadata_raw)
        except (KeyError, TypeError, json.JSONDecodeError) as exc:
            raise ValueError(f"{path}: malformed CCSD checkpoint metadata") from exc
        if canonical_json_sha256(metadata) != store.attrs.get("manifest_sha256"):
            raise ValueError(f"{path}: CCSD checkpoint manifest checksum mismatch")
        expected_identity_sha = canonical_json_sha256(identity)
        if (
            metadata.get("schema") != CCSD_RESTART_SCHEMA
            or metadata.get("identity") != identity
            or metadata.get("identity_sha256") != expected_identity_sha
        ):
            raise ValueError(f"{path}: CCSD checkpoint identity drift")
        completed_cycle = metadata.get("completed_cycle")
        energy = metadata.get("ccsd_correlation_hartree")
        if (
            not isinstance(completed_cycle, int)
            or completed_cycle < 1
            or not isinstance(energy, (int, float))
            or not np.isfinite(float(energy))
        ):
            raise ValueError(f"{path}: invalid CCSD checkpoint state")

        nocc = int(coupled_cluster.nocc)
        nvir = int(coupled_cluster.nmo - nocc)
        expected_shapes = {"t1": (nocc, nvir), "t2": (nocc, nocc, nvir, nvir)}
        loaded = {}
        for name, shape in expected_shapes.items():
            if name not in store or tuple(store[name].shape) != shape:
                raise ValueError(f"{path}: {name} amplitude shape drift")
            dataset: Any = store[name]
            tensor_metadata = metadata.get("tensors", {}).get(name, {})
            if tensor_metadata.get("shape") != list(shape) or tensor_metadata.get(
                "dtype"
            ) != str(np.dtype(dataset.dtype)):
                raise ValueError(f"{path}: {name} amplitude manifest drift")
            target = coupled_cluster.pool.new(name, shape, dataset.dtype)
            digest = hashlib.sha256()
            digest.update(str(np.dtype(dataset.dtype)).encode())
            digest.update(json.dumps(shape).encode())
            for start in range(shape[0]):
                chunk = np.ascontiguousarray(dataset[start : start + 1])
                if not np.isfinite(chunk).all():
                    raise ValueError(
                        f"{path}: {name} amplitude checkpoint is non-finite"
                    )
                digest.update(memoryview(chunk).cast("B"))
                target[start : start + 1] = chunk
            if digest.hexdigest() != tensor_metadata.get("sha256"):
                raise ValueError(f"{path}: {name} amplitude checksum mismatch")
            loaded[name] = target
    return metadata, loaded["t1"], loaded["t2"]


class ByteQCEnergyCheckpoint:
    """Wrap ByteQC energy evaluation and persist every completed CCSD cycle."""

    def __init__(
        self,
        path: Path,
        identity: dict[str, Any],
        *,
        base_cycle: int = 0,
    ) -> None:
        self.path = path
        self.identity = identity
        self.base_cycle = base_cycle
        self.energy_calls = 0
        self.last_checkpoint: dict[str, Any] | None = None

    def wrap(self, original):
        def checkpointing_energy(t1=None, t2=None, eris=None, with_em=False):
            result = original(t1, t2, eris, with_em)
            # ByteQC evaluates the initial MP2 guess once before its loop. Every
            # later call occurs after damping/DIIS and is the completed cycle
            # printed to the engine log. Persist before the next update starts.
            if self.energy_calls > 0:
                self.last_checkpoint = write_ccsd_checkpoint(
                    self.path,
                    identity=self.identity,
                    completed_cycle=self.base_cycle + self.energy_calls,
                    ccsd_correlation_hartree=float(result[0]),
                    t1=t1,
                    t2=t2,
                )
            self.energy_calls += 1
            return result

        return checkpointing_energy


@contextmanager
def byteqc_energy_checkpoint(
    coupled_cluster: Any,
    path: Path,
    identity: dict[str, Any],
    *,
    base_cycle: int = 0,
):
    checkpoint = ByteQCEnergyCheckpoint(path, identity, base_cycle=base_cycle)
    original = coupled_cluster.energy
    coupled_cluster.energy = checkpoint.wrap(original)
    try:
        yield checkpoint
    finally:
        coupled_cluster.energy = original


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


def byteqc_frozen_core_nr_e2(
    eri: Any,
    mo_coeff: Any,
    orbs_slice: tuple[int, int, int, int],
    aosym: str = "s1",
    mosym: str = "s1",
    out: Any = None,
    ao_loc: Any = None,
    *,
    array_module: Any = None,
    lib_module: Any = None,
    contraction_module: Any = None,
) -> Any:
    """Transform one DF block without reusing an MO-sized buffer as AO scratch.

    ByteQC 2.5 passes the previous MO-transformed block back as ``out``.  Its
    stock ``nr_e2`` also passes that buffer to ``unpack_tril`` even though the
    AO matrix is larger whenever frozen core makes ``nmo < nao``.  Keep the
    reusable buffer for the final MO result only and allocate the AO scratch
    independently.  This preserves bounded DF blocks for cc-pVQZ instead of
    forcing the entire auxiliary basis onto the GPU at once.
    """
    if array_module is None or lib_module is None or contraction_module is None:
        import cupy
        from byteqc import lib as byteqc_lib
        from byteqc.cucc import culib

        array_module = cupy
        lib_module = byteqc_lib
        contraction_module = culib

    assert eri.flags.c_contiguous
    assert mo_coeff.dtype == np.double
    k0, k1, l0, l1 = orbs_slice
    kc = k1 - k0
    lc = l1 - l0
    kl_count = kc * lc
    nrow = eri.shape[0]
    if nrow * kl_count == 0:
        return array_module.empty((nrow, kl_count))
    assert aosym == "s2" and mosym == "s1" and kc == lc and ao_loc is None

    # Never expose the smaller, prior MO result to the larger AO unpack.
    mat = lib_module.unpack_tril(eri)
    tmp = contraction_module.contraction("ij", mo_coeff[k0:], "mjk", mat, "mik")
    mat = None
    result = contraction_module.gemm(
        "N", "T", tmp.reshape(-1, tmp.shape[-1]), mo_coeff[l0:], buf=out
    )
    return result.reshape((tmp.shape[0], tmp.shape[1], -1))


@contextmanager
def byteqc_frozen_core_transform(coupled_cluster: Any, *, dfccsd_module: Any = None):
    """Install the frozen-core-safe DF transform for one bounded CC job."""
    if dfccsd_module is None:
        from byteqc.cucc import dfccsd as dfccsd_module

    naux = int(coupled_cluster.with_df.get_naoaux())
    blockdim = min(int(coupled_cluster.with_df.blockdim), naux)
    if blockdim < 1:
        raise ValueError("ByteQC DF block size must be positive")
    coupled_cluster.with_df.blockdim = blockdim
    original = dfccsd_module.nr_e2
    dfccsd_module.nr_e2 = byteqc_frozen_core_nr_e2
    try:
        yield blockdim
    finally:
        dfccsd_module.nr_e2 = original


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


def restart_checkpoint_path(args: argparse.Namespace, source: Path) -> Path:
    path = (
        args.restart_checkpoint or args.output.with_suffix(".ccsd-restart.h5")
    ).resolve()
    protected_paths = {source.resolve(), args.output.resolve(), args.log.resolve()}
    if path in protected_paths:
        raise ValueError(
            "CCSD restart checkpoint must differ from XYZ, receipt, and engine log"
        )
    return path


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
    driver_source_sha256 = sha256_path(Path(__file__).resolve())
    identity = {
        "engine_commit": byteqc_commit,
        "driver_source_sha256": driver_source_sha256,
        "pyscf_version": pyscf.__version__,
        "method": "canonical-ccsd(t)",
        "basis": args.basis,
        "freeze_core_orbitals": frozen,
        "density_fit": True,
        "df_aux_blocking": "bounded-frozen-core-safe-transform-v2",
        "charge": args.charge,
        "spin_2s": args.spin,
    }
    restart_path = restart_checkpoint_path(args, source)
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
    coupled_cluster.max_cycle = 100
    coupled_cluster.conv_tol = 1e-8
    restart_identity = {
        **identity,
        "checkpoint_schema": CCSD_RESTART_SCHEMA,
        "input_sha256": source_sha,
        "scf_contract": {
            "method": "RHF-density-fit",
            "max_cycle": mean_field.max_cycle,
            "conv_tol": mean_field.conv_tol,
            "energy_hartree": scf_energy,
            "mo_coeff_sha256": numpy_sha256(mean_field.mo_coeff),
            "mo_energy_sha256": numpy_sha256(mean_field.mo_energy),
            "mo_occ_sha256": numpy_sha256(mean_field.mo_occ),
        },
        "ccsd_convergence_contract": {
            "max_cycle_semantics": "100-per-launch; completed_cycle-is-cumulative",
            "max_cycle": coupled_cluster.max_cycle,
            "conv_tol": coupled_cluster.conv_tol,
            "conv_tol_normt": coupled_cluster.conv_tol_normt,
            "diis_space": coupled_cluster.diis_space,
            "diis_start_cycle": coupled_cluster.diis_start_cycle,
            "diis_start_energy_diff": coupled_cluster.diis_start_energy_diff,
            "iterative_damping": coupled_cluster.iterative_damping,
            "level_shift": coupled_cluster.level_shift,
        },
    }
    # ByteQC 2.5 reuses an MO-sized DF result as larger AO-unpack scratch when
    # frozen core makes nmo < nao. Keep bounded auxiliary blocks and replace
    # only that transform route for both CCSD and the later triples ERIs.
    with byteqc_frozen_core_transform(coupled_cluster):
        eris = coupled_cluster.ao2mo()
        restart = load_ccsd_checkpoint(
            restart_path,
            identity=restart_identity,
            coupled_cluster=coupled_cluster,
        )
        base_cycle = 0
        initial_t1 = initial_t2 = None
        restart_metadata = None
        convergence_marker = None
        if restart is not None:
            restart_metadata, initial_t1, initial_t2 = restart
            base_cycle = int(restart_metadata["completed_cycle"])
            convergence_marker = load_ccsd_converged_marker(
                restart_path, restart_metadata
            )
        if convergence_marker is not None:
            assert restart_metadata is not None
            coupled_cluster.t1 = initial_t1
            coupled_cluster.t2 = initial_t2
            coupled_cluster.e_corr = float(restart_metadata["ccsd_correlation_hartree"])
            coupled_cluster.converged = True
            ccsd_correlation = coupled_cluster.e_corr
            final_checkpoint_metadata = restart_metadata
        else:
            with byteqc_energy_checkpoint(
                coupled_cluster,
                restart_path,
                restart_identity,
                base_cycle=base_cycle,
            ) as energy_checkpoint:
                ccsd_correlation, _, _ = coupled_cluster.kernel(
                    t1=initial_t1,
                    t2=initial_t2,
                    eris=eris,
                )
            final_checkpoint_metadata = energy_checkpoint.last_checkpoint
            if not coupled_cluster.converged:
                raise RuntimeError("ByteQC CCSD did not converge")
            if final_checkpoint_metadata is None:
                raise RuntimeError("ByteQC CCSD converged without a durable checkpoint")
            convergence_marker = write_ccsd_converged_marker(
                restart_path, final_checkpoint_metadata
            )
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
        "ccsd_restart": {
            "schema": CCSD_RESTART_SCHEMA,
            "used": restart_metadata is not None,
            "resumed_manifest": restart_metadata,
            "final_ccsd_manifest": final_checkpoint_metadata,
            "convergence_marker": convergence_marker,
            "identity_sha256": canonical_json_sha256(restart_identity),
        },
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    atomic_json(args.output, payload)
    restart_path.unlink(missing_ok=True)
    ccsd_converged_marker_path(restart_path).unlink(missing_ok=True)
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
            job.add_argument("--restart-checkpoint", type=Path)
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
