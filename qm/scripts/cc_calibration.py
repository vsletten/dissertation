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
import inspect
import json
import os
import platform
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
CCSD_RESTART_SCHEMA_V1 = "byteqc-ccsd-amplitudes-v1"
CCSD_RESTART_SCHEMA = "byteqc-ccsd-amplitudes-v2"
CCSD_CONVERGED_SCHEMA = "byteqc-ccsd-converged-v2"
CHECKPOINT_ARRAY_NAMES = (
    "mo_coeff",
    "mo_energy",
    "mo_occ",
    "scf_energy",
    "t1",
    "t2",
)
CHECKPOINT_DTYPE = np.dtype("float64")
SCF_ENERGY_VALIDATION_ATOL = 1e-8
SCF_ORBITAL_VALIDATION_ATOL = 1e-5
BYTEQC_SCF_SETTINGS = {
    "method": "RHF",
    "density_fit": True,
    "max_cycle": 150,
    "conv_tol": 1e-10,
}
BYTEQC_CCSD_SETTINGS = {
    "max_cycle": 100,
    "conv_tol": 1e-8,
    "conv_tol_normt": 1e-5,
    "diis": True,
    "diis_space": 6,
    "diis_start_cycle": 0,
    "diis_start_energy_diff": 1e9,
    "iterative_damping": 1.0,
    "level_shift": 0.0,
}
BYTEQC_RESTART_CONTRACT = {
    "schema": CCSD_RESTART_SCHEMA,
    "generation": "single-atomic-hdf5-scf-and-amplitudes",
    "resume": "restore-scf-before-ccsd-and-ao2mo;rerun-ccsd-kernel",
    "cycle_numbering": "completed-cycle-is-cumulative;max-cycle-is-per-launch",
    "convergence_marker": "validated-provenance-only;never-skips-kernel",
    "array_dtype": str(CHECKPOINT_DTYPE),
}
BYTEQC_TRANSFORM_CONTRACT = {
    "density_fit": True,
    "implementation": "bounded-frozen-core-safe-transform-v2",
    "aux_blocking": "min(configured-blockdim,naux)",
}


class CCSDLaunchCheckpointed(RuntimeError):
    """Stop a non-converged launch after one crash-durable CCSD cycle."""

    def __init__(self, checkpoint: dict[str, Any]) -> None:
        self.checkpoint = checkpoint
        cycle = checkpoint.get("completed_cycle")
        super().__init__(
            f"ByteQC launch stopped after durable CCSD cycle {cycle}; resume required"
        )


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


def byteqc_restart_identity(
    *,
    input_sha256: str,
    engine_commit: str,
    pyscf_version: str,
    cupy_version: str,
    numpy_version: str,
    python_version: str,
    basis: str,
    charge: int,
    spin_2s: int,
    frozen_core_orbitals: int,
) -> dict[str, Any]:
    """Build the stable numerical/contract identity for a ByteQC restart.

    SCF results and resource-only controls are deliberately absent. A resume
    restores the exact stored SCF state instead of requiring a repeat SCF to
    reproduce orbital bytes.
    """
    return {
        "input_sha256": input_sha256,
        "runtime": {
            "byteqc_commit": engine_commit,
            "pyscf_version": pyscf_version,
            "cupy_version": cupy_version,
            "numpy_version": numpy_version,
            "python_version": python_version,
        },
        "electronic_structure": {
            "method": "canonical-ccsd(t)",
            "basis": basis,
            "charge": charge,
            "spin_2s": spin_2s,
            "frozen_core_orbitals": frozen_core_orbitals,
        },
        "restart_contract": dict(BYTEQC_RESTART_CONTRACT),
        "transform_contract": dict(BYTEQC_TRANSFORM_CONTRACT),
        "scf_settings": dict(BYTEQC_SCF_SETTINGS),
        "ccsd_settings": dict(BYTEQC_CCSD_SETTINGS),
    }


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


def _checkpoint_dataset(store: Any, name: str, value: Any) -> dict[str, Any]:
    shape = tuple(int(item) for item in value.shape)
    dtype = np.dtype(value.dtype)
    if dtype != CHECKPOINT_DTYPE:
        raise ValueError(f"{name}: checkpoint dtype must be {CHECKPOINT_DTYPE}")
    if shape:
        if shape[0] < 1:
            raise ValueError(f"{name}: checkpoint array must be non-empty")
        chunks = (1, *shape[1:])
    else:
        chunks = None
    dataset = store.create_dataset(name, shape=shape, dtype=dtype, chunks=chunks)
    digest = hashlib.sha256()
    digest.update(str(dtype).encode())
    digest.update(json.dumps(shape).encode())
    if not shape:
        chunk = np.ascontiguousarray(value).reshape(())
        if not np.isfinite(chunk).all():
            raise ValueError(f"{name}: checkpoint contains non-finite values")
        digest.update(memoryview(chunk.reshape(1)).cast("B"))
        dataset[()] = chunk
    else:
        for start in range(shape[0]):
            chunk = _numpy_chunk(value, start, start + 1)
            if chunk.shape != (1, *shape[1:]):
                raise ValueError(f"{name}: checkpoint chunk shape drift")
            if not np.isfinite(chunk).all():
                raise ValueError(f"{name}: checkpoint contains non-finite values")
            digest.update(memoryview(chunk).cast("B"))
            dataset[start : start + 1] = chunk
    return {
        "shape": list(shape),
        "dtype": str(dtype),
        "nbytes": int(value.nbytes),
        "sha256": digest.hexdigest(),
    }


def _checkpoint_manifest(store: Any, path: Path) -> dict[str, Any]:
    try:
        metadata_raw = store.attrs["manifest_json"]
        if isinstance(metadata_raw, bytes):
            metadata_raw = metadata_raw.decode()
        if not isinstance(metadata_raw, str):
            raise TypeError("manifest_json is not text")
        metadata = json.loads(metadata_raw)
    except (KeyError, TypeError, json.JSONDecodeError) as exc:
        raise ValueError(f"{path}: malformed CCSD checkpoint metadata") from exc
    schema = metadata.get("schema")
    if schema == CCSD_RESTART_SCHEMA_V1:
        raise ValueError(
            f"{path}: legacy v1 CCSD checkpoint is unsupported; "
            "it will not be migrated or overwritten"
        )
    if schema != CCSD_RESTART_SCHEMA:
        raise ValueError(f"{path}: unsupported CCSD checkpoint schema {schema!r}")
    if canonical_json_sha256(metadata) != store.attrs.get("manifest_sha256"):
        raise ValueError(f"{path}: CCSD checkpoint manifest checksum mismatch")
    if set(metadata.get("arrays", {})) != set(CHECKPOINT_ARRAY_NAMES):
        raise ValueError(f"{path}: CCSD checkpoint array manifest drift")
    return metadata


def _validate_checkpoint_identity(
    path: Path, metadata: dict[str, Any], identity: dict[str, Any]
) -> None:
    expected_identity_sha = canonical_json_sha256(identity)
    if (
        metadata.get("identity") != identity
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


def _validate_dataset_header(
    store: Any,
    path: Path,
    metadata: dict[str, Any],
    name: str,
    *,
    expected_shape: tuple[int, ...] | None = None,
) -> Any:
    if name not in store:
        raise ValueError(f"{path}: missing {name} checkpoint dataset")
    dataset = store[name]
    shape = tuple(int(item) for item in dataset.shape)
    dtype = np.dtype(dataset.dtype)
    array_metadata = metadata["arrays"].get(name, {})
    if expected_shape is not None and shape != expected_shape:
        raise ValueError(f"{path}: {name} checkpoint shape drift")
    if dtype != CHECKPOINT_DTYPE:
        raise ValueError(f"{path}: {name} checkpoint dtype drift")
    if (
        array_metadata.get("shape") != list(shape)
        or array_metadata.get("dtype") != str(dtype)
        or array_metadata.get("nbytes") != int(dataset.nbytes)
        or not isinstance(array_metadata.get("sha256"), str)
    ):
        raise ValueError(f"{path}: {name} checkpoint manifest drift")
    return dataset


def _read_validated_dataset(
    store: Any,
    path: Path,
    metadata: dict[str, Any],
    name: str,
    *,
    expected_shape: tuple[int, ...] | None = None,
    target: Any = None,
) -> Any:
    dataset = _validate_dataset_header(
        store, path, metadata, name, expected_shape=expected_shape
    )
    shape = tuple(int(item) for item in dataset.shape)
    digest = hashlib.sha256()
    digest.update(str(CHECKPOINT_DTYPE).encode())
    digest.update(json.dumps(shape).encode())
    if not shape:
        value = np.asarray(dataset[()], dtype=CHECKPOINT_DTYPE).reshape(())
        if not np.isfinite(value).all():
            raise ValueError(f"{path}: {name} checkpoint is non-finite")
        digest.update(memoryview(value.reshape(1)).cast("B"))
        result = value
    else:
        result = np.empty(shape, dtype=CHECKPOINT_DTYPE) if target is None else target
        for start in range(shape[0]):
            chunk = np.ascontiguousarray(dataset[start : start + 1])
            if not np.isfinite(chunk).all():
                raise ValueError(f"{path}: {name} checkpoint is non-finite")
            digest.update(memoryview(chunk).cast("B"))
            result[start : start + 1] = chunk
    if digest.hexdigest() != metadata["arrays"][name]["sha256"]:
        raise ValueError(f"{path}: {name} checkpoint checksum mismatch")
    return result


def _validate_dataset_payload(
    store: Any,
    path: Path,
    metadata: dict[str, Any],
    name: str,
    *,
    expected_shape: tuple[int, ...],
) -> None:
    """Stream-validate a large dataset without retaining a second host copy."""
    dataset = _validate_dataset_header(
        store, path, metadata, name, expected_shape=expected_shape
    )
    digest = hashlib.sha256()
    digest.update(str(CHECKPOINT_DTYPE).encode())
    digest.update(json.dumps(expected_shape).encode())
    for start in range(expected_shape[0]):
        chunk = np.ascontiguousarray(dataset[start : start + 1])
        if not np.isfinite(chunk).all():
            raise ValueError(f"{path}: {name} checkpoint is non-finite")
        digest.update(memoryview(chunk).cast("B"))  # type: ignore[arg-type]
    if digest.hexdigest() != metadata["arrays"][name]["sha256"]:
        raise ValueError(f"{path}: {name} checkpoint checksum mismatch")


def write_ccsd_checkpoint(
    path: Path,
    *,
    identity: dict[str, Any],
    completed_cycle: int,
    ccsd_correlation_hartree: float,
    mo_coeff: Any,
    mo_energy: Any,
    mo_occ: Any,
    scf_energy: float,
    t1: Any,
    t2: Any,
) -> dict[str, Any]:
    """Atomically persist one exact SCF and evaluated ByteQC CCSD state."""
    import h5py

    if not isinstance(completed_cycle, int) or isinstance(completed_cycle, bool):
        raise ValueError("CCSD checkpoint cycle must be an integer")
    if completed_cycle < 1:
        raise ValueError("CCSD checkpoint cycle must be positive")
    correlation_energy = float(ccsd_correlation_hartree)
    scf_energy_array = np.asarray(float(scf_energy), dtype=CHECKPOINT_DTYPE)
    if not np.isfinite(correlation_energy) or not np.isfinite(scf_energy_array).all():
        raise ValueError("CCSD and SCF checkpoint energies must be finite")
    arrays = {
        "mo_coeff": mo_coeff,
        "mo_energy": mo_energy,
        "mo_occ": mo_occ,
        "scf_energy": scf_energy_array,
        "t1": t1,
        "t2": t2,
    }
    shapes = {
        name: tuple(int(item) for item in value.shape) for name, value in arrays.items()
    }
    nmo = shapes["mo_coeff"][1] if len(shapes["mo_coeff"]) == 2 else -1
    nocc = shapes["t1"][0] if len(shapes["t1"]) == 2 else -1
    nvir = shapes["t1"][1] if len(shapes["t1"]) == 2 else -1
    if (
        nmo < 1
        or shapes["mo_energy"] != (nmo,)
        or shapes["mo_occ"] != (nmo,)
        or shapes["scf_energy"] != ()
        or nocc < 1
        or nvir < 1
        or shapes["t2"] != (nocc, nocc, nvir, nvir)
    ):
        raise ValueError("CCSD checkpoint SCF/amplitude shape contract drift")
    for name, value in arrays.items():
        if np.dtype(value.dtype) != CHECKPOINT_DTYPE:
            raise ValueError(f"{name}: checkpoint dtype must be {CHECKPOINT_DTYPE}")

    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        try:
            with h5py.File(path, "r") as existing:
                _checkpoint_manifest(existing, path)
        except OSError as exc:
            raise ValueError(
                f"{path}: existing file is not a supported CCSD checkpoint; "
                "refusing to overwrite it"
            ) from exc
    checkpoint_bytes = sum(int(value.nbytes) for value in arrays.values())
    # disk_usage().free already excludes the old checkpoint. Atomic replacement
    # needs one complete temporary generation plus a 1 GiB safety reserve.
    required_free = checkpoint_bytes + (1 << 30)
    available_free = shutil.disk_usage(path.parent).free
    if available_free < required_free:
        raise RuntimeError(
            "insufficient disk for atomic CCSD checkpoint: "
            f"need {required_free} bytes, have {available_free}"
        )
    temporary = path.with_name(f".{path.name}.{time.time_ns()}.tmp")
    try:
        with h5py.File(temporary, "w") as store:
            array_manifest = {
                name: _checkpoint_dataset(store, name, value)
                for name, value in arrays.items()
            }
            metadata = {
                "schema": CCSD_RESTART_SCHEMA,
                "identity": identity,
                "identity_sha256": canonical_json_sha256(identity),
                "completed_cycle": completed_cycle,
                "ccsd_correlation_hartree": correlation_energy,
                "checkpoint_bytes": checkpoint_bytes,
                "arrays": array_manifest,
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
        marker_path = ccsd_converged_marker_path(path)
        marker_path.unlink(missing_ok=True)
        directory = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
        os.replace(temporary, path)
        directory = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        temporary.unlink(missing_ok=True)
    return metadata


def _validate_restored_scf_state(
    path: Path,
    *,
    identity: dict[str, Any],
    mean_field: Any,
    mo_coeff: np.ndarray,
    mo_energy: np.ndarray,
    mo_occ: np.ndarray,
    scf_energy: float,
    t1_shape: tuple[int, ...],
) -> None:
    """Reject physically inconsistent orbitals before mutating the SCF object."""
    electron_count = int(mean_field.mol.nelectron)
    if electron_count < 1 or electron_count % 2:
        raise ValueError(f"{path}: restored RHF electron count is invalid")
    occupation_is_closed_shell = np.logical_or(
        np.isclose(mo_occ, 0.0, atol=1e-12, rtol=0.0),
        np.isclose(mo_occ, 2.0, atol=1e-12, rtol=0.0),
    )
    if not occupation_is_closed_shell.all() or not np.isclose(
        mo_occ.sum(), electron_count, atol=1e-10, rtol=0.0
    ):
        raise ValueError(f"{path}: restored RHF occupations are invalid")
    if np.any(np.diff(mo_energy) < -1e-10):
        raise ValueError(f"{path}: restored RHF orbital energies are not ordered")

    overlap = np.asarray(mean_field.get_ovlp(), dtype=CHECKPOINT_DTYPE)
    if overlap.shape != (mo_coeff.shape[0], mo_coeff.shape[0]):
        raise ValueError(f"{path}: restored RHF overlap shape drift")
    metric = mo_coeff.T @ overlap @ mo_coeff
    if not np.allclose(
        metric,
        np.eye(mo_coeff.shape[1], dtype=CHECKPOINT_DTYPE),
        atol=1e-8,
        rtol=1e-8,
    ):
        raise ValueError(f"{path}: restored RHF orbitals are not orthonormal")

    density = mean_field.make_rdm1(mo_coeff, mo_occ)
    hcore = mean_field.get_hcore()
    veff = mean_field.get_veff(mean_field.mol, density)
    recomputed_energy = float(mean_field.energy_tot(dm=density, h1e=hcore, vhf=veff))
    if not np.isclose(
        recomputed_energy, scf_energy, atol=SCF_ENERGY_VALIDATION_ATOL, rtol=0.0
    ):
        raise ValueError(f"{path}: restored RHF total energy is inconsistent")
    fock = np.asarray(
        mean_field.get_fock(h1e=hcore, s1e=overlap, vhf=veff, dm=density),
        dtype=CHECKPOINT_DTYPE,
    )
    if fock.shape != overlap.shape or not np.isfinite(fock).all():
        raise ValueError(f"{path}: restored RHF Fock matrix is invalid")
    mo_fock = mo_coeff.T @ fock @ mo_coeff
    off_diagonal = mo_fock - np.diag(np.diag(mo_fock))
    if not np.allclose(
        off_diagonal, 0.0, atol=SCF_ORBITAL_VALIDATION_ATOL, rtol=0.0
    ) or not np.allclose(
        np.diag(mo_fock),
        mo_energy,
        atol=SCF_ORBITAL_VALIDATION_ATOL,
        rtol=1e-9,
    ):
        raise ValueError(f"{path}: restored RHF orbitals do not diagonalize Fock")

    electronic_structure = identity.get("electronic_structure")
    if not isinstance(electronic_structure, dict):
        raise ValueError(f"{path}: restart identity lacks electronic structure")
    frozen = electronic_structure.get("frozen_core_orbitals")
    if not isinstance(frozen, int) or isinstance(frozen, bool) or frozen < 0:
        raise ValueError(f"{path}: restart identity has invalid frozen core")
    occupied = electron_count // 2
    expected_t1_shape = (occupied - frozen, mo_coeff.shape[1] - occupied)
    if min(expected_t1_shape) < 1 or t1_shape != expected_t1_shape:
        raise ValueError(f"{path}: amplitudes disagree with restored RHF occupations")


def restore_ccsd_checkpoint_scf(
    path: Path,
    *,
    identity: dict[str, Any],
    mean_field: Any,
) -> dict[str, Any] | None:
    """Validate and restore exact v2 SCF state before CCSD construction."""
    import h5py

    if not path.exists():
        return None
    with h5py.File(path, "r") as store:
        metadata = _checkpoint_manifest(store, path)
        _validate_checkpoint_identity(path, metadata, identity)
        t1_shape = tuple(metadata["arrays"]["t1"].get("shape", ()))
        t2_shape = tuple(metadata["arrays"]["t2"].get("shape", ()))
        if len(t1_shape) != 2 or t2_shape != (
            t1_shape[0],
            t1_shape[0],
            t1_shape[1],
            t1_shape[1],
        ):
            raise ValueError(f"{path}: amplitude checkpoint shape drift")
        _validate_dataset_header(store, path, metadata, "t1", expected_shape=t1_shape)
        _validate_dataset_header(store, path, metadata, "t2", expected_shape=t2_shape)

        coeff_dataset = _validate_dataset_header(store, path, metadata, "mo_coeff")
        coeff_shape = tuple(int(item) for item in coeff_dataset.shape)
        if len(coeff_shape) != 2 or min(coeff_shape) < 1:
            raise ValueError(f"{path}: mo_coeff checkpoint shape drift")
        nmo = coeff_shape[1]
        expected_nao = int(mean_field.mol.nao_nr())
        if coeff_shape[0] != expected_nao:
            raise ValueError(f"{path}: mo_coeff checkpoint shape drift")
        mo_coeff = _read_validated_dataset(
            store, path, metadata, "mo_coeff", expected_shape=coeff_shape
        )
        mo_energy = _read_validated_dataset(
            store, path, metadata, "mo_energy", expected_shape=(nmo,)
        )
        mo_occ = _read_validated_dataset(
            store, path, metadata, "mo_occ", expected_shape=(nmo,)
        )
        scf_energy = _read_validated_dataset(
            store, path, metadata, "scf_energy", expected_shape=()
        )

    _validate_restored_scf_state(
        path,
        identity=identity,
        mean_field=mean_field,
        mo_coeff=mo_coeff,
        mo_energy=mo_energy,
        mo_occ=mo_occ,
        scf_energy=float(scf_energy),
        t1_shape=t1_shape,
    )

    # Validate the multi-gigabyte amplitude payload before expensive AO→MO work,
    # but only after the cheap SCF state has passed its physical checks.
    with h5py.File(path, "r") as store:
        current_metadata = _checkpoint_manifest(store, path)
        _validate_checkpoint_identity(path, current_metadata, identity)
        if current_metadata != metadata:
            raise ValueError(
                f"{path}: CCSD checkpoint generation changed during resume"
            )
        _validate_dataset_payload(store, path, metadata, "t1", expected_shape=t1_shape)
        _validate_dataset_payload(store, path, metadata, "t2", expected_shape=t2_shape)

    mean_field.mo_coeff = mo_coeff
    mean_field.mo_energy = mo_energy
    mean_field.mo_occ = mo_occ
    mean_field.e_tot = float(scf_energy)
    mean_field.converged = True
    return metadata


def load_ccsd_checkpoint_amplitudes(
    path: Path,
    *,
    identity: dict[str, Any],
    coupled_cluster: Any,
    expected_metadata: dict[str, Any],
) -> tuple[Any, Any]:
    """Validate/load amplitudes into ByteQC-owned selected buffer backends."""
    import h5py

    with h5py.File(path, "r") as store:
        metadata = _checkpoint_manifest(store, path)
        _validate_checkpoint_identity(path, metadata, identity)
        if metadata != expected_metadata:
            raise ValueError(
                f"{path}: CCSD checkpoint generation changed during resume"
            )
        nocc = int(coupled_cluster.nocc)
        nvir = int(coupled_cluster.nmo - nocc)
        expected_shapes = {"t1": (nocc, nvir), "t2": (nocc, nocc, nvir, nvir)}
        loaded = {}
        for name, shape in expected_shapes.items():
            target = coupled_cluster.pool.new(name, shape, CHECKPOINT_DTYPE)
            loaded[name] = _read_validated_dataset(
                store,
                path,
                metadata,
                name,
                expected_shape=shape,
                target=target,
            )
    return loaded["t1"], loaded["t2"]


def scf_checkpoint_state(mean_field: Any, scf_energy: float) -> dict[str, Any]:
    """Snapshot the exact finite host SCF state used to construct ByteQC."""
    state = {
        "mo_coeff": np.ascontiguousarray(mean_field.mo_coeff),
        "mo_energy": np.ascontiguousarray(mean_field.mo_energy),
        "mo_occ": np.ascontiguousarray(mean_field.mo_occ),
        "scf_energy": float(scf_energy),
    }
    nmo = state["mo_coeff"].shape[1] if state["mo_coeff"].ndim == 2 else -1
    if (
        nmo < 1
        or state["mo_energy"].shape != (nmo,)
        or state["mo_occ"].shape != (nmo,)
        or any(
            np.dtype(state[name].dtype) != CHECKPOINT_DTYPE
            for name in ("mo_coeff", "mo_energy", "mo_occ")
        )
        or not all(
            np.isfinite(state[name]).all()
            for name in ("mo_coeff", "mo_energy", "mo_occ")
        )
        or not np.isfinite(state["scf_energy"])
    ):
        raise ValueError("SCF checkpoint state shape/dtype/finiteness drift")
    return state


def _constructor_supports_orbitals(constructor: Any) -> bool:
    parameters = inspect.signature(constructor).parameters.values()
    names = {parameter.name for parameter in parameters}
    return {"mo_coeff", "mo_occ"} <= names or any(
        parameter.kind is inspect.Parameter.VAR_KEYWORD for parameter in parameters
    )


def prepare_byteqc_coupled_cluster(
    *,
    mean_field: Any,
    cucc_module: Any,
    frozen: int,
    gpulim: int,
    restart_path: Path,
    restart_identity: dict[str, Any],
) -> tuple[Any, float, dict[str, Any], dict[str, Any] | None]:
    """Restore/run SCF, then construct CCSD with the exact orbital state."""
    for name, value in BYTEQC_SCF_SETTINGS.items():
        if name not in {"method", "density_fit"}:
            setattr(mean_field, name, value)
    restart_metadata = restore_ccsd_checkpoint_scf(
        restart_path, identity=restart_identity, mean_field=mean_field
    )
    if restart_metadata is None:
        scf_energy = float(mean_field.kernel())
        if not mean_field.converged:
            raise RuntimeError("ByteQC calibration RHF did not converge")
    else:
        scf_energy = float(mean_field.e_tot)
    scf_state = scf_checkpoint_state(mean_field, scf_energy)

    constructor_kwargs = {"frozen": frozen, "gpulim": gpulim}
    if _constructor_supports_orbitals(cucc_module.CCSD):
        constructor_kwargs.update(
            mo_coeff=mean_field.mo_coeff,
            mo_occ=mean_field.mo_occ,
        )
    coupled_cluster = cucc_module.CCSD(mean_field, **constructor_kwargs)
    for name, value in BYTEQC_CCSD_SETTINGS.items():
        setattr(coupled_cluster, name, value)
    return coupled_cluster, scf_energy, scf_state, restart_metadata


class ByteQCEnergyCheckpoint:
    """Wrap ByteQC energy evaluation and persist every completed CCSD cycle."""

    def __init__(
        self,
        path: Path,
        identity: dict[str, Any],
        scf_state: dict[str, Any],
        *,
        base_cycle: int = 0,
    ) -> None:
        self.path = path
        self.identity = identity
        self.scf_state = scf_state
        self.base_cycle = base_cycle
        self.energy_calls = 0
        self.max_checkpoints = 1
        self.checkpoints_written = 0
        self.last_checkpoint: dict[str, Any] | None = None

    def wrap(self, original):
        def checkpointing_energy(t1=None, t2=None, eris=None, with_em=False):
            result = original(t1, t2, eris, with_em)
            # ByteQC evaluates the initial MP2 guess once before its loop. Every
            # later call occurs after damping/DIIS and is the completed cycle
            # printed to the engine log. Persist before the next update starts.
            if self.energy_calls > 0:
                checkpoint = write_ccsd_checkpoint(
                    self.path,
                    identity=self.identity,
                    completed_cycle=self.base_cycle + self.energy_calls,
                    ccsd_correlation_hartree=float(result[0]),
                    **self.scf_state,
                    t1=t1,
                    t2=t2,
                )
                self.last_checkpoint = checkpoint
                self.checkpoints_written += 1
            self.energy_calls += 1
            return result

        return checkpointing_energy

    def wrap_update_amps(self, original):
        def bounded_update(*args, **kwargs):
            if self.checkpoints_written >= self.max_checkpoints:
                if self.last_checkpoint is None:
                    raise RuntimeError("CCSD checkpoint count has no durable state")
                raise CCSDLaunchCheckpointed(self.last_checkpoint)
            return original(*args, **kwargs)

        return bounded_update


@contextmanager
def byteqc_energy_checkpoint(
    coupled_cluster: Any,
    path: Path,
    identity: dict[str, Any],
    scf_state: dict[str, Any],
    *,
    base_cycle: int = 0,
):
    checkpoint = ByteQCEnergyCheckpoint(
        path, identity, scf_state, base_cycle=base_cycle
    )
    original_energy = coupled_cluster.energy
    original_update_amps = coupled_cluster.update_amps
    coupled_cluster.energy = checkpoint.wrap(original_energy)
    coupled_cluster.update_amps = checkpoint.wrap_update_amps(original_update_amps)
    try:
        yield checkpoint
    finally:
        coupled_cluster.energy = original_energy
        coupled_cluster.update_amps = original_update_amps


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
    restart_identity = byteqc_restart_identity(
        input_sha256=source_sha,
        engine_commit=byteqc_commit,
        pyscf_version=pyscf.__version__,
        cupy_version=cupy.__version__,
        numpy_version=np.__version__,
        python_version=platform.python_version(),
        basis=args.basis,
        charge=args.charge,
        spin_2s=args.spin,
        frozen_core_orbitals=frozen,
    )
    mean_field = scf.RHF(molecule).density_fit()
    coupled_cluster, scf_energy, scf_state, restart_metadata = (
        prepare_byteqc_coupled_cluster(
            mean_field=mean_field,
            cucc_module=cucc,
            frozen=frozen,
            gpulim=args.gpu_memory_gb << 30,
            restart_path=restart_path,
            restart_identity=restart_identity,
        )
    )
    # ByteQC 2.5 reuses an MO-sized DF result as larger AO-unpack scratch when
    # frozen core makes nmo < nao. Keep bounded auxiliary blocks and replace
    # only that transform route for both CCSD and the later triples ERIs.
    with byteqc_frozen_core_transform(coupled_cluster):
        eris = coupled_cluster.ao2mo()
        base_cycle = 0
        initial_t1 = initial_t2 = None
        resumed_convergence_marker = None
        if restart_metadata is not None:
            base_cycle = int(restart_metadata["completed_cycle"])
            resumed_convergence_marker = load_ccsd_converged_marker(
                restart_path, restart_metadata
            )
            initial_t1, initial_t2 = load_ccsd_checkpoint_amplitudes(
                restart_path,
                identity=restart_identity,
                coupled_cluster=coupled_cluster,
                expected_metadata=restart_metadata,
            )
        # A marker is provenance only. Re-enter ByteQC's convergence contract
        # from the stored amplitudes instead of trusting it as final evidence.
        with byteqc_energy_checkpoint(
            coupled_cluster,
            restart_path,
            restart_identity,
            scf_state,
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
            "scf_state_source": "checkpoint"
            if restart_metadata is not None
            else "fresh",
            "resumed_manifest": restart_metadata,
            "resumed_convergence_marker": resumed_convergence_marker,
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
    if not isinstance(payload.get("input_sha256"), str):
        raise ValueError(f"{path}: missing geometry SHA-256")
    for field in ("scf_hartree", "correlation_hartree", "total_hartree"):
        try:
            value = float(payload[field])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{path}: missing finite {field}") from exc
        if not np.isfinite(value):
            raise ValueError(f"{path}: non-finite {field}")
    expected_total = float(payload["scf_hartree"]) + float(
        payload["correlation_hartree"]
    )
    if not np.isclose(
        float(payload["total_hartree"]), expected_total, rtol=0.0, atol=1e-9
    ):
        raise ValueError(f"{path}: total energy arithmetic mismatch")
    return payload


def validate_role_receipts(
    matrices: dict[str, dict[str, dict[str, dict[str, Any]]]],
) -> None:
    """Require one exact geometry and electronic state per reaction role."""
    for role in ("reactant", "ts"):
        receipts = [
            receipt for matrix in matrices.values() for receipt in matrix[role].values()
        ]
        input_hashes = {receipt.get("input_sha256") for receipt in receipts}
        states = {
            (receipt.get("charge"), receipt.get("spin_2s")) for receipt in receipts
        }
        if None in input_hashes or len(input_hashes) != 1:
            raise ValueError(f"{role} receipts do not share one exact geometry")
        if len(states) != 1:
            raise ValueError(f"{role} receipts do not share one electronic state")


def barrier_kj(reactant: dict[str, Any], ts: dict[str, Any]) -> float:
    return (
        float(ts["total_hartree"]) - float(reactant["total_hartree"])
    ) * HARTREE_TO_KJ


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
    validate_role_receipts({"canonical": canonical, "dlpno": dlpno})
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


def summarize_focal(args: argparse.Namespace) -> dict[str, Any]:
    """Build a canonical-TZ focal barrier without inventing canonical TS/QZ.

    The expensive canonical QZ calculation exists only for the reactant.  The
    complete TightPNO TZ/QZ pair supplies the CBS basis correction, while the
    independently complete canonical/DLPNO TZ barriers supply the method gate.
    """
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
    validate_role_receipts({"canonical": canonical, "dlpno": dlpno})

    canonical_tz_barrier = barrier_kj(
        canonical["reactant"]["tz"], canonical["ts"]["tz"]
    )
    dlpno_tz_barrier = barrier_kj(dlpno["reactant"]["tz"], dlpno["ts"]["tz"])
    dlpno_qz_barrier = barrier_kj(dlpno["reactant"]["qz"], dlpno["ts"]["qz"])
    dlpno_cbs_barrier = (
        cbs_total(dlpno["ts"]["tz"], dlpno["ts"]["qz"])
        - cbs_total(dlpno["reactant"]["tz"], dlpno["reactant"]["qz"])
    ) * HARTREE_TO_KJ
    agreement_delta = dlpno_tz_barrier - canonical_tz_barrier
    basis_correction = dlpno_cbs_barrier - dlpno_tz_barrier
    reactant_qz_delta_hartree = float(dlpno["reactant"]["qz"]["total_hartree"]) - float(
        canonical["reactant"]["qz"]["total_hartree"]
    )
    payload = {
        "method": (
            "canonical CCSD(T)/cc-pVTZ barrier + TightPNO "
            "DLPNO-CCSD(T) TZ/QZ CBS basis correction"
        ),
        "route": "canonical-tz-plus-dlpno-cbs-minus-dlpno-tz",
        "canonical_ts_qz_computed": False,
        "canonical_ts_qz_reason": (
            "directive-adjusted route: no second canonical QZ slow-mode calculation"
        ),
        "hf_extrapolation_alpha": HF_CBS_ALPHA,
        "correlation_extrapolation_power": 3,
        "canonical_tz_barrier_kj": canonical_tz_barrier,
        "dlpno_tz_barrier_kj": dlpno_tz_barrier,
        "dlpno_qz_barrier_kj": dlpno_qz_barrier,
        "dlpno_cbs_barrier_kj": dlpno_cbs_barrier,
        "dlpno_cbs_minus_tz_basis_correction_kj": basis_correction,
        "focal_point_barrier_kj": canonical_tz_barrier + basis_correction,
        "agreement_basis": "cc-pVTZ",
        "dlpno_minus_canonical_kj": agreement_delta,
        "gate_limit_kj": 2.0,
        "gate_pass": abs(agreement_delta) <= 2.0,
        "reactant_qz_anchor": {
            "canonical_total_hartree": float(
                canonical["reactant"]["qz"]["total_hartree"]
            ),
            "dlpno_total_hartree": float(dlpno["reactant"]["qz"]["total_hartree"]),
            "dlpno_minus_canonical_hartree": reactant_qz_delta_hartree,
            "dlpno_minus_canonical_kj": reactant_qz_delta_hartree * HARTREE_TO_KJ,
        },
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


def _parse_receipt_grid(values: list[str]) -> dict[str, dict[str, Path]]:
    result: dict[str, dict[str, Path]] = {}
    for value in values:
        fields = value.split("=", 2)
        if len(fields) != 3:
            raise ValueError("receipt must be ROLE=BASIS=PATH")
        role, basis, path = fields
        normalized = basis.lower().replace("cc-pv", "").replace("z", "")
        if role not in {"reactant", "ts"} or normalized not in {"t", "q"}:
            raise ValueError(f"invalid receipt selector {value}")
        basis_key = "tz" if normalized == "t" else "qz"
        if basis_key in result.setdefault(role, {}):
            raise ValueError(f"duplicate receipt selector {role}/{basis_key}")
        result[role][basis_key] = Path(path)
    return result


def receipt_grid(values: list[str]) -> dict[str, dict[str, Path]]:
    result = _parse_receipt_grid(values)
    if set(result) != {"reactant", "ts"} or any(
        set(paths) != {"tz", "qz"} for paths in result.values()
    ):
        raise ValueError("receipts must cover reactant/ts at TZ/QZ")
    return result


def focal_canonical_receipt_grid(values: list[str]) -> dict[str, dict[str, Path]]:
    result = _parse_receipt_grid(values)
    expected = {"reactant": {"tz", "qz"}, "ts": {"tz"}}
    if {role: set(paths) for role, paths in result.items()} != expected:
        raise ValueError(
            "focal canonical receipts must cover reactant TZ/QZ and TS TZ exactly"
        )
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
    for command in ("summarize", "summarize-focal"):
        summary = subparsers.add_parser(command)
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
    elif args.command == "summarize":
        args.canonical = receipt_grid(args.canonical)
        args.dlpno = receipt_grid(args.dlpno)
        payload = summarize(args)
    else:
        args.canonical = focal_canonical_receipt_grid(args.canonical)
        args.dlpno = receipt_grid(args.dlpno)
        payload = summarize_focal(args)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
