#!/usr/bin/env python
"""D2c SCT campaign preflight and immutable resume contract.

Only ``--dry-run`` is implemented.  The command verifies the frozen D2b bundle,
validates all four transition-state atom mappings against the existing D2b reaction
templates, and writes one non-overwriting pending receipt.  It deliberately performs
no IRC, Hessian, high-level, VAG, or tunnelling calculation.
"""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import subprocess
import sys
import time
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

if __name__ == "__main__":
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    from quarry.etiquette import bootstrap_cli

    _guard = argparse.ArgumentParser(add_help=False)
    _guard.add_argument("--threads", type=int, default=16)
    _guard.add_argument("--nice", type=int, default=10)
    _guard_args, _ = _guard.parse_known_args()
    if not 1 <= _guard_args.threads <= 16:
        _guard.error("--threads must be <= 16 and at least 1")
    if _guard_args.nice < 10:
        _guard.error("--nice must be >= 10")
    _ETIQUETTE = bootstrap_cli(
        "d2c_sct_campaign",
        default_run_root=Path(__file__).resolve().parent.parent / "runs",
    )

import numpy as np  # noqa: E402

from quarry.clusters import Cluster  # noqa: E402
from quarry.native_hessian import NativeHessianResult  # noqa: E402
from quarry.pipeline import frequency_geometry_fingerprint  # noqa: E402
from quarry.reaction_path import (  # noqa: E402
    hessian_eigenvalues_to_wavenumbers_cm,
    project_vibrational_hessian,
)
from scripts import d2c_input_bundle
from scripts.production_energetics import load_xyz_like
from scripts.surface_rate_protocol import reactions

SCHEMA = "d2c-sct-campaign-preflight-v2"
PREFLIGHT_RECEIPT = "preflight.json"
DEFAULT_BUNDLE_ROOT = (
    Path(__file__).resolve().parent.parent
    / "data"
    / "D2c-instanton-tier"
    / "d2b-inputs"
)
DFT_SETTINGS: dict[str, Any] = {
    "xc": "pwb6k",
    "basis": "def2-svp",
    "solvent": None,
    "dispersion": "d3bj",
    "composite": None,
    "grid_level": None,
    "density_fit": True,
    "use_gpu": True,
}
DEPENDENCY_DISTRIBUTIONS = ("numpy", "pyscf", "geometric", "sella", "ase")
GPU4PYSCF_DISTRIBUTIONS = (
    "gpu4pyscf-cuda12x",
    "gpu4pyscf-cuda11x",
    "gpu4pyscf",
)
CUPY_DISTRIBUTIONS = ("cupy-cuda13x", "cupy-cuda12x", "cupy-cuda11x", "cupy")
DEPENDENCY_VERSION_KEYS = frozenset(
    (
        *DEPENDENCY_DISTRIBUTIONS,
        "gpu_backend",
        "gpu4pyscf_distribution",
        "gpu4pyscf",
        "cupy_distribution",
        "cupy",
        "cuda_runtime",
        "cuda_driver",
        "cuda_device_count",
    )
)
# These geometry hashes are trusted source constants, deliberately independent of
# the mutable manifest being verified.  They bind atom row identity even if a copied
# bundle and all of its ordinary byte receipts are coherently rewritten.
TRUSTED_CANONICAL_TS_GEOMETRY_SHA256 = {
    "h-co-1w-oside": "7c121ddaddfb47932b95301593524498fdb7ed98faaf741ee383e04c49880106",
    "h-co-1w-cside": "4265cd1ee5dfcc082027bafa8322587aa505f200589adcaf9cbb97e5b7eb7f2b",
    "h-h2co-ch3o-1w": (
        "c5602e65b972f2703d467781cca1c8e1299eb06d68d13fd3889cbf75c294d6bf"
    ),
    "h-h2co-h2-hco-1w": (
        "a1b804833d74cd77b7382cd516b84af7981fd37f18417ae7a858a7898a8225bf"
    ),
}
# Ground-state isotopic masses.  Carbon-12 is exact by definition; the others are
# the neutral-atom masses used by this bounded H/C/O campaign.
ISOTOPIC_MASSES_AMU = {
    "H": 1.00782503223,
    "C": 12.0,
    "O": 15.99491461957,
}
REFERENCE_MASS_AMU = 1.0
BOUNDS: dict[str, Any] = {
    "transition_state_qualification": {
        "physical_fmax_ev_per_angstrom_exclusive_maximum": 0.02,
        "negative_eigenvalue_tolerance_hartree_per_bohr2_amu": 1.0e-8,
        "required_significant_imaginary_mode_count": 1,
        "minimum_reaction_imaginary_wavenumber_cm": 200.0,
        "minimum_mapped_reaction_vector_overlap": 0.20,
        "minimum_irc_tangent_overlap": 0.80,
        "spectator_imaginary_modes_forbidden": True,
        "full_index_gate_precedes_transverse_projection": True,
    },
    "irc": {
        "directions": ["forward", "reverse"],
        "algorithm": "sella-gonzalez-schlegel",
        "step_size_angstrom": 0.05,
        "maximum_steps_per_direction": 200,
        "outer_fmax_ev_per_angstrom": 0.05,
        "inner_fmax_ev_per_angstrom": 0.01,
    },
    "hessian": {
        "coverage": "every retained IRC point including TS and endpoints",
        "maximum_retained_points_per_direction": 201,
        "cartesian_hessian_units": "hartree / bohr^2",
        "transverse_negative_eigenvalue_tolerance": 1.0e-8,
        "transverse_negative_eigenvalue_tolerance_units": "hartree / bohr^2 / amu",
        "transverse_eigenvalue_units": "hartree / bohr^2 / amu",
        "required_transverse_mode_count": "3N-7",
    },
    "sct": {
        "temperature_kelvin": [
            12.0,
            13.5,
            15.0,
            16.5,
            20.0,
            30.0,
            40.0,
            50.0,
            60.0,
            75.0,
            100.0,
            150.0,
            200.0,
            250.0,
            300.0,
        ],
        "quadrature_order": 96,
        "path_grid_size": 4097,
        "reference_mass_amu": REFERENCE_MASS_AMU,
        "straightness_tolerance_per_angstrom": 1.0e-12,
    },
}
ROUTE_STAGE_CONTRACT: dict[str, dict[str, Any]] = {
    "transition_state_qualification": {
        "required": True,
        "receipt": "ts-qualification/receipt.json",
        "requirement": (
            "fresh physical gradient and strict full 3N-6 first-order-saddle gate"
        ),
    },
    "irc_forward": {
        "required": True,
        "receipt": "irc-forward/receipt.json",
        "requirement": "bounded IRC in the forward direction",
    },
    "irc_reverse": {
        "required": True,
        "receipt": "irc-reverse/receipt.json",
        "requirement": "bounded IRC in the reverse direction",
    },
    "typed_irc_path": {
        "required": True,
        "receipt": "path/receipt.json",
        "requirement": "typed endpoints and one oriented reactant-to-product path",
    },
    "hessian_every_path_point": {
        "required": True,
        "receipt": "hessians/receipt.json",
        "requirement": "Cartesian Hessian at every retained path point",
    },
    "vibrationally_adiabatic_potential": {
        "required": True,
        "receipt": "vag/receipt.json",
        "requirement": "positive 3N-7 transverse modes and exact ZPE conversion",
    },
    "high_level_correction": {
        "required": True,
        "receipt": "high-level/receipt.json",
        "requirement": "same-length path correction on hash-bound geometries",
    },
    "sct": {
        "required": True,
        "receipt": "sct/receipt.json",
        "requirement": "mode-resolved curvature mass and bounded SCT quadrature",
    },
}
CAMPAIGN_STAGE_CONTRACT: dict[str, dict[str, Any]] = {
    "branching_common_reference_gate": {
        "required": True,
        "receipt": "branching-common-reference.json",
        "requirement": (
            "both competing H2CO channels use one exact common-reactant receipt"
        ),
    },
    "final_freeze": {
        "required": True,
        "receipt": "final-freeze.json",
        "requirement": "verify every stage and freeze final immutable result hashes",
    },
}
_FORBIDDEN_ACCEPTED_RESULTS = (
    "accepted-result.json",
    "final-freeze.json",
    "final-result.json",
    "results.json",
)


def validate_transition_state_gate(
    cluster: Cluster,
    masses_amu: Any,
    native_hessian: NativeHessianResult,
    *,
    expected_settings_fingerprint: str,
    reaction_vector_mass_scaled: Any,
) -> dict[str, Any]:
    """Apply the strict fresh first-order-saddle gate before any IRC call."""

    coordinates = np.asarray(cluster.coords, dtype=float)
    if native_hessian.geometry_fingerprint != frequency_geometry_fingerprint(cluster):
        raise ValueError(
            "native Hessian geometry fingerprint does not match TS geometry"
        )
    if (
        not expected_settings_fingerprint
        or native_hessian.settings_fingerprint != expected_settings_fingerprint
    ):
        raise ValueError("native Hessian settings fingerprint does not match campaign")
    masses = np.asarray(masses_amu, dtype=float)
    expected_masses = np.asarray(
        [ISOTOPIC_MASSES_AMU[symbol] for symbol in cluster.symbols], dtype=float
    )
    if masses.shape != expected_masses.shape or not np.array_equal(
        masses, expected_masses
    ):
        raise ValueError(
            "mass vector does not match the receipt-bound isotopic standard"
        )
    bounds = BOUNDS["transition_state_qualification"]
    fmax_limit = bounds["physical_fmax_ev_per_angstrom_exclusive_maximum"]
    if native_hessian.physical_fmax_ev_per_angstrom >= fmax_limit:
        raise ValueError(
            "transition state is not stationary: physical fmax "
            f"{native_hessian.physical_fmax_ev_per_angstrom:.12g} >= "
            f"{fmax_limit:.12g} eV/A"
        )
    modes = project_vibrational_hessian(
        coordinates,
        masses,
        native_hessian.cartesian_hessian_hartree_per_bohr2,
    )
    tolerance = bounds["negative_eigenvalue_tolerance_hartree_per_bohr2_amu"]
    negative = modes.eigenvalues < -tolerance
    near_zero = np.abs(modes.eigenvalues) <= tolerance
    if int(np.count_nonzero(negative)) != 1 or bool(np.any(near_zero)):
        raise ValueError(
            "transition state must have exactly one significant negative "
            "full-vibrational mode and no zero/noise-floor modes"
        )
    if bool(np.any(modes.eigenvalues[~negative] <= 0.0)):
        raise ValueError(
            "every non-reaction full-vibrational mode must be strictly positive"
        )
    wavenumbers = hessian_eigenvalues_to_wavenumbers_cm(modes.eigenvalues)
    imaginary = float(abs(wavenumbers[negative][0]))
    minimum_imaginary = bounds["minimum_reaction_imaginary_wavenumber_cm"]
    if imaginary < minimum_imaginary:
        raise ValueError(
            f"reaction imaginary mode {imaginary:.12g} cm^-1 is below "
            f"{minimum_imaginary:.12g} cm^-1"
        )
    reaction_vector = np.asarray(reaction_vector_mass_scaled, dtype=float)
    if reaction_vector.shape == coordinates.shape:
        reaction_vector = reaction_vector.reshape(-1)
    if reaction_vector.shape != (3 * len(cluster.symbols),) or not np.all(
        np.isfinite(reaction_vector)
    ):
        raise ValueError(
            "reaction vector must be a finite mass-scaled (N,3) or (3N,) vector"
        )
    reaction_vector = modes.vibrational_basis @ (
        modes.vibrational_basis.T @ reaction_vector
    )
    reaction_norm = float(np.linalg.norm(reaction_vector))
    if reaction_norm <= 1.0e-12:
        raise ValueError("reaction vector has no non-rigid vibrational component")
    reaction_vector /= reaction_norm
    unstable_mode = modes.mass_weighted_eigenvectors[
        np.flatnonzero(negative)[0]
    ].reshape(-1)
    reaction_overlap = float(abs(np.dot(unstable_mode, reaction_vector)))
    minimum_overlap = bounds["minimum_mapped_reaction_vector_overlap"]
    if reaction_overlap < minimum_overlap:
        raise ValueError(
            f"imaginary mode overlap {reaction_overlap:.12g} is below mapped-reaction "
            f"minimum {minimum_overlap:.12g}"
        )
    gradient_bytes = np.ascontiguousarray(
        native_hessian.gradient_hartree_per_bohr, dtype="<f8"
    ).tobytes()
    hessian_bytes = np.ascontiguousarray(
        native_hessian.cartesian_hessian_hartree_per_bohr2, dtype="<f8"
    ).tobytes()
    return {
        "accepted": True,
        "physical_fmax_ev_per_angstrom": (native_hessian.physical_fmax_ev_per_angstrom),
        "physical_fmax_exclusive_limit_ev_per_angstrom": fmax_limit,
        "vibrational_mode_count": int(modes.eigenvalues.size),
        "imaginary_mode_count": 1,
        "imaginary_wavenumber_cm": imaginary,
        "minimum_imaginary_wavenumber_cm": minimum_imaginary,
        "mapped_reaction_vector_overlap": reaction_overlap,
        "minimum_mapped_reaction_vector_overlap": minimum_overlap,
        "eigenvalues_hartree_per_bohr2_amu": modes.eigenvalues.tolist(),
        "masses_amu": masses.tolist(),
        "gradient_sha256": hashlib.sha256(gradient_bytes).hexdigest(),
        "canonical_hessian_sha256": hashlib.sha256(hessian_bytes).hexdigest(),
        "requested_backend": native_hessian.requested_backend,
        "actual_backend": native_hessian.actual_backend,
        "gpu_fallback_used": native_hessian.gpu_fallback_used,
        "geometry_fingerprint": native_hessian.geometry_fingerprint,
        "settings_fingerprint": native_hessian.settings_fingerprint,
    }


def _require_sha(value: str, *, length: int, label: str) -> str:
    if len(value) != length or any(
        character not in "0123456789abcdef" for character in value
    ):
        raise ValueError(f"{label} must be {length} lowercase hexadecimal digits")
    return value


def _git_sha(repository_root: Path) -> str:
    status = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=all"],
        cwd=repository_root,
        check=True,
        capture_output=True,
        text=True,
    )
    if status.stdout:
        raise RuntimeError(
            "D2c campaign execution requires a clean tracked and untracked Git worktree"
        )
    completed = subprocess.run(
        ["git", "rev-parse", "--verify", "HEAD"],
        cwd=repository_root,
        check=True,
        capture_output=True,
        text=True,
    )
    return _require_sha(completed.stdout.strip(), length=40, label="Git SHA")


def _installed_distribution(label: str, candidates: tuple[str, ...]) -> tuple[str, str]:
    installed: list[tuple[str, str]] = []
    for distribution in candidates:
        try:
            installed.append((distribution, importlib.metadata.version(distribution)))
        except importlib.metadata.PackageNotFoundError:
            continue
    if not installed:
        raise RuntimeError(
            f"required D2c {label} distribution is not installed; expected one of "
            + ", ".join(candidates)
        )
    if len(installed) != 1:
        names = ", ".join(distribution for distribution, _ in installed)
        raise RuntimeError(
            f"ambiguous D2c {label} distributions are installed: {names}"
        )
    return installed[0]


def _dependency_versions() -> dict[str, str]:
    """Inventory the exact CPU packages, GPU backend, and live CUDA runtime."""

    versions: dict[str, str] = {}
    for distribution in DEPENDENCY_DISTRIBUTIONS:
        try:
            versions[distribution] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError as exc:
            raise RuntimeError(
                f"required D2c dependency is not installed: {distribution}"
            ) from exc

    gpu_distribution, gpu_version = _installed_distribution(
        "GPU4PySCF", GPU4PYSCF_DISTRIBUTIONS
    )
    cupy_distribution, cupy_version = _installed_distribution(
        "CuPy", CUPY_DISTRIBUTIONS
    )
    try:
        __import__("gpu4pyscf")
        cupy = __import__("cupy")
        runtime = cupy.cuda.runtime
        runtime_version = runtime.runtimeGetVersion()
        driver_version = runtime.driverGetVersion()
        device_count = runtime.getDeviceCount()
    except Exception as exc:
        raise RuntimeError(
            "required D2c GPU backend or CUDA runtime is unavailable"
        ) from exc
    if not isinstance(device_count, int) or device_count < 1:
        raise RuntimeError("required D2c CUDA runtime exposes no GPU devices")

    versions.update(
        {
            "gpu_backend": "gpu4pyscf",
            "gpu4pyscf_distribution": gpu_distribution,
            "gpu4pyscf": gpu_version,
            "cupy_distribution": cupy_distribution,
            "cupy": cupy_version,
            "cuda_runtime": str(runtime_version),
            "cuda_driver": str(driver_version),
            "cuda_device_count": str(device_count),
        }
    )
    return versions


def _canonical_hash(payload: dict[str, Any]) -> str:
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _pending_stages(contract: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    return {
        name: {
            **stage,
            "state": "pending",
            "accepted": False,
            "resume_requires_campaign_identity": True,
        }
        for name, stage in contract.items()
    }


def _declared_stage_receipt_paths(run_root: Path) -> tuple[Path, ...]:
    route_receipts = (
        run_root / route / stage["receipt"]
        for route in d2c_input_bundle.ROUTES
        for stage in ROUTE_STAGE_CONTRACT.values()
    )
    campaign_receipts = (
        run_root / stage["receipt"] for stage in CAMPAIGN_STAGE_CONTRACT.values()
    )
    return (*route_receipts, *campaign_receipts)


def _semantic_atom_identity_labels(route: str, template: Any) -> tuple[str, ...]:
    """Bind the index semantics defined by the D2b reaction-template builders."""

    if template.n_water != 1:
        raise ValueError(f"D2c route must use exactly one template water: {route}")
    if template.family == "h-co":
        labels = (
            "carbon_monoxide_carbon",
            "carbon_monoxide_oxygen",
            "incoming_hydrogen",
            "water_1_oxygen",
            "water_1_hydrogen_a",
            "water_1_hydrogen_b",
        )
        expected_symbols = ("C", "O", "H", "O", "H", "H")
        expected_scan_pair = (0, 2)
    elif template.family == "h-h2co-ch3o":
        labels = (
            "formaldehyde_carbon",
            "formaldehyde_oxygen",
            "formaldehyde_hydrogen_a",
            "formaldehyde_hydrogen_b",
            "incoming_hydrogen",
            "water_1_oxygen",
            "water_1_hydrogen_a",
            "water_1_hydrogen_b",
        )
        expected_symbols = ("C", "O", "H", "H", "H", "O", "H", "H")
        expected_scan_pair = (0, 4)
    elif template.family == "h-h2co-h2-hco":
        labels = (
            "formaldehyde_carbon",
            "formaldehyde_oxygen",
            "abstracted_formaldehyde_hydrogen",
            "retained_formaldehyde_hydrogen",
            "incoming_hydrogen",
            "water_1_oxygen",
            "water_1_hydrogen_a",
            "water_1_hydrogen_b",
        )
        expected_symbols = ("C", "O", "H", "H", "H", "O", "H", "H")
        expected_scan_pair = (2, 4)
    else:
        raise ValueError(f"unsupported D2c reaction-template family: {route}")
    if (
        tuple(template.cluster.symbols) != expected_symbols
        or (template.scan_i, template.scan_j) != expected_scan_pair
    ):
        raise ValueError(f"semantic atom template order drifted: {route}")
    return labels


def _route_inventory(bundle_root: Path, manifest: dict[str, Any]) -> dict[str, Any]:
    templates = reactions(gpu=True, basis="def2-svp")
    inventory: dict[str, Any] = {}
    for route in d2c_input_bundle.ROUTES:
        template = templates.get(route)
        if template is None:
            raise ValueError(f"missing existing reaction template: {route}")
        manifest_route = manifest["routes"][route]
        if manifest_route.get("method") != DFT_SETTINGS:
            raise ValueError(f"frozen DFT settings drifted: {route}")
        expected_geometry_hash = TRUSTED_CANONICAL_TS_GEOMETRY_SHA256[route]
        declared_canonical_hashes = manifest_route.get(
            "canonical_checkpoint_geometry_sha256"
        )
        observed_geometry_hash = d2c_input_bundle.geometry_hash_xyz(
            bundle_root / route / "ts.xyz"
        )
        if (
            not isinstance(declared_canonical_hashes, dict)
            or declared_canonical_hashes.get("ts.xyz") != expected_geometry_hash
            or observed_geometry_hash != expected_geometry_hash
        ):
            raise ValueError(
                f"canonical transition-state geometry/mapping identity drifted: {route}"
            )
        transition_state = load_xyz_like(
            bundle_root / route / "ts.xyz",
            template.cluster,
            name=f"{route}-ts",
        )
        try:
            masses = [
                ISOTOPIC_MASSES_AMU[symbol] for symbol in transition_state.symbols
            ]
        except KeyError as exc:
            raise ValueError(
                f"no frozen isotopic mass for element {exc.args[0]}"
            ) from exc
        if not all(math.isfinite(mass) and mass > 0.0 for mass in masses):
            raise RuntimeError("internal isotopic mass table is invalid")
        atom_identity_labels = _semantic_atom_identity_labels(route, template)
        if tuple(transition_state.symbols) != tuple(template.cluster.symbols):
            raise ValueError(f"transition-state atom symbols drifted: {route}")
        atom_mapping_sha256 = _canonical_hash(
            {
                "route": route,
                "canonical_transition_state_geometry_sha256": expected_geometry_hash,
                "atoms": [
                    {"index": index, "symbol": symbol, "semantic_identity": label}
                    for index, (symbol, label) in enumerate(
                        zip(
                            transition_state.symbols,
                            atom_identity_labels,
                            strict=True,
                        )
                    )
                ],
            }
        )
        inventory[route] = {
            "symbols": list(transition_state.symbols),
            "atom_identity_labels": list(atom_identity_labels),
            "atom_mapping_sha256": atom_mapping_sha256,
            "canonical_transition_state_geometry_sha256": expected_geometry_hash,
            "masses_amu": masses,
            "transition_state_sha256": d2c_input_bundle.sha256_path(
                bundle_root / route / "ts.xyz"
            ),
            "stages": _pending_stages(ROUTE_STAGE_CONTRACT),
        }
    if tuple(inventory) != d2c_input_bundle.ROUTES or len(inventory) != 4:
        raise RuntimeError(
            "D2c preflight must enumerate exactly the four frozen routes"
        )
    return inventory


@contextmanager
def _exclusive_run_claim(run_root: Path) -> Iterator[None]:
    """Hold a crash-recoverable exclusive claim on one campaign run root."""

    run_root.mkdir(parents=True, exist_ok=True)
    claim_path = run_root / ".preflight.lock"
    flags = os.O_RDWR | os.O_CREAT | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(claim_path, flags, 0o600)
    locked = False
    try:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(
                f"D2c run root is already claimed by another preflight: {run_root}"
            ) from exc
        locked = True
        claim = json.dumps(
            {"pid": os.getpid(), "claimed_unix_ns": time.time_ns()},
            sort_keys=True,
        ).encode()
        os.ftruncate(descriptor, 0)
        os.write(descriptor, claim + b"\n")
        os.fsync(descriptor)
        yield
    finally:
        if locked:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        os.close(descriptor)


def _atomic_non_overwriting_json(path: Path, payload: dict[str, Any]) -> None:
    """Publish complete JSON atomically while refusing an existing target."""

    if path.exists() or path.is_symlink():
        raise FileExistsError(f"preflight receipt already exists: {path}")
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0)
    data = (
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode()
    descriptor = os.open(temporary, flags, 0o600)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            descriptor = -1
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        try:
            os.link(temporary, path, follow_symlinks=False)
        except FileExistsError as exc:
            raise FileExistsError(f"preflight receipt already exists: {path}") from exc
        directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        temporary.unlink(missing_ok=True)


def create_preflight_receipt(
    bundle_root: Path,
    run_root: Path,
    *,
    git_sha: str | None = None,
    dependency_versions: dict[str, str] | None = None,
    created_utc: str | None = None,
) -> dict[str, Any]:
    """Validate all immutable inputs and publish one pending campaign receipt.

    ``git_sha`` and ``dependency_versions`` are injectable to keep unit tests free
    of live Git and environment discovery.  Normal CLI use always resolves both.
    No existing receipt or accepted/final result is deleted, reused, or overwritten.
    """

    bundle_root = bundle_root.resolve()
    run_root = run_root.resolve()
    manifest = d2c_input_bundle.verify_bundle(bundle_root)
    manifest_sha = d2c_input_bundle.sha256_path(bundle_root / "manifest.json")
    revision = _require_sha(
        _git_sha(Path(__file__).resolve().parents[2]) if git_sha is None else git_sha,
        length=40,
        label="Git SHA",
    )
    dependencies = dict(
        _dependency_versions() if dependency_versions is None else dependency_versions
    )
    if set(dependencies) != DEPENDENCY_VERSION_KEYS:
        raise ValueError("dependency version inventory must be exact and complete")
    if any(not isinstance(value, str) or not value for value in dependencies.values()):
        raise ValueError("dependency versions must be non-empty strings")
    if (
        dependencies["gpu_backend"] != "gpu4pyscf"
        or dependencies["gpu4pyscf_distribution"] not in GPU4PYSCF_DISTRIBUTIONS
        or dependencies["cupy_distribution"] not in CUPY_DISTRIBUTIONS
        or not dependencies["cuda_runtime"].isdigit()
        or not dependencies["cuda_driver"].isdigit()
    ):
        raise ValueError("dependency inventory has an invalid GPU backend identity")
    try:
        cuda_device_count = int(dependencies["cuda_device_count"])
    except ValueError as exc:
        raise ValueError("CUDA device count must be a positive integer") from exc
    if cuda_device_count < 1:
        raise ValueError("CUDA device count must be a positive integer")
    routes = _route_inventory(bundle_root, manifest)

    identity_payload = {
        "schema": SCHEMA,
        "bundle_manifest_sha256": manifest_sha,
        "git_sha": revision,
        "dft_settings": DFT_SETTINGS,
        "dependencies": dependencies,
        "python": platform.python_version(),
        "mass_standard": "ground-state neutral isotopic masses",
        "reference_mass_amu": REFERENCE_MASS_AMU,
        "routes": {
            route: {
                "symbols": record["symbols"],
                "atom_identity_labels": record["atom_identity_labels"],
                "atom_mapping_sha256": record["atom_mapping_sha256"],
                "canonical_transition_state_geometry_sha256": record[
                    "canonical_transition_state_geometry_sha256"
                ],
                "masses_amu": record["masses_amu"],
                "transition_state_sha256": record["transition_state_sha256"],
            }
            for route, record in routes.items()
        },
        "bounds": BOUNDS,
        "required_route_stages": ROUTE_STAGE_CONTRACT,
        "required_campaign_stages": CAMPAIGN_STAGE_CONTRACT,
    }
    identity = _canonical_hash(identity_payload)
    receipt = {
        "schema": SCHEMA,
        "state": "pending",
        "accepted_result": None,
        "dry_run": True,
        "created_utc": created_utc
        or time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "identity": identity,
        "campaign": identity_payload,
        "routes": routes,
        "campaign_stages": _pending_stages(CAMPAIGN_STAGE_CONTRACT),
        "resume_policy": {
            "identity_match_required": True,
            "stage_receipt_hash_required": True,
            "partial_or_unbound_stage_reuse_forbidden": True,
            "accepted_results_before_final_freeze_forbidden": True,
        },
    }

    with _exclusive_run_claim(run_root):
        for name in _FORBIDDEN_ACCEPTED_RESULTS:
            existing = run_root / name
            if existing.exists() or existing.is_symlink():
                raise FileExistsError(
                    f"refusing stale accepted/final result in run root: {existing}"
                )
        for existing in _declared_stage_receipt_paths(run_root):
            if existing.exists() or existing.is_symlink():
                raise FileExistsError(
                    f"refusing stale stage receipt in run root: {existing}"
                )
        _atomic_non_overwriting_json(run_root / PREFLIGHT_RECEIPT, receipt)
    return receipt


def main(
    argv: list[str] | None = None,
    *,
    git_sha: str | None = None,
    dependency_versions: dict[str, str] | None = None,
) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--bundle-root", type=Path, default=DEFAULT_BUNDLE_ROOT)
    parser.add_argument("--run-root", type=Path, required=True)
    # Pre-parsed by bootstrap_cli for an executable invocation; retained here so
    # direct tests and the complete parser enforce the same campaign limits.
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--nice", type=int, default=10)
    parser.add_argument("--log")
    args = parser.parse_args(argv)
    if not args.dry_run:
        parser.error("production execution is not implemented; --dry-run is required")
    if not 1 <= args.threads <= 16:
        parser.error("--threads must be <= 16 and at least 1")
    if args.nice < 10:
        parser.error("--nice must be >= 10")

    receipt = create_preflight_receipt(
        args.bundle_root,
        args.run_root,
        git_sha=git_sha,
        dependency_versions=dependency_versions,
    )
    print(
        json.dumps(
            {
                "status": "pending",
                "dry_run": True,
                "identity": receipt["identity"],
                "routes": list(receipt["routes"]),
                "receipt": str(args.run_root.resolve() / PREFLIGHT_RECEIPT),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
