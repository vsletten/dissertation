#!/usr/bin/env python
"""Bounded scalar-energy curvatures along two receipt-pinned D2c modes.

This diagnostic is evidence-only forever.  It evaluates one shared center and
q=+/-0.05,+/-0.10 Bohr*sqrt(amu) along the independently verified Richardson
unstable and lowest-positive modes.  It never evaluates a gradient or Hessian
and cannot qualify a transition state, launch IRC/SCT, or authorize D3b.
"""

from __future__ import annotations

# Import order is deliberate: CLI etiquette is established before scientific imports.
# ruff: noqa: E402
import argparse
import hashlib
import importlib.metadata
import inspect
import json
import math
import os
import platform
import stat
import subprocess
import sys
from contextlib import suppress
from dataclasses import replace
from datetime import UTC, datetime
from pathlib import Path
from types import CodeType
from typing import Any

_THIS_MODULE_IMPORT_CODE = sys._getframe().f_code
_QM_ROOT = Path(__file__).resolve().parent.parent
if str(_QM_ROOT) not in sys.path:
    sys.path.insert(0, str(_QM_ROOT))

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    _guard = argparse.ArgumentParser(add_help=False)
    _guard.add_argument("--threads", type=int, default=16)
    _guard.add_argument("--nice", type=int, default=10)
    _guard_args, _ = _guard.parse_known_args()
    if any(
        arg in {"--gpu", "--gpu-mem-gb", "--cuda"}
        or arg.startswith(("--gpu=", "--gpu-mem-gb=", "--cuda="))
        for arg in sys.argv[1:]
    ):
        _guard.error("GPU/CUDA flags are forbidden for this CPU-only diagnostic")
    if not 1 <= _guard_args.threads <= 16:
        _guard.error("--threads must be between 1 and 16")
    if _guard_args.nice < 10:
        _guard.error("--nice must be at least 10")
    _ETIQUETTE = bootstrap_cli(
        "d2c_mode_projected_energy_diagnostic",
        default_run_root=Path("/mnt/data/vsletten/dissertation-data"),
    )

import numpy as np
from pyscf import df as _pyscf_df
from pyscf import dft as _pyscf_dft
from pyscf import dispersion as _pyscf_dispersion
from pyscf import gto as _pyscf_gto
from pyscf import scf as _pyscf_scf

_PRELOADED_CPU_MODULES = (
    _pyscf_df,
    _pyscf_dft,
    _pyscf_dispersion,
    _pyscf_gto,
    _pyscf_scf,
)

from quarry.pipeline import (
    BOHR_TO_ANGSTROM,
    _make_scf,
    build_mol,
    frequency_geometry_fingerprint,
)
from quarry.reaction_path import (
    hessian_eigenvalues_to_wavenumbers_cm,
    project_vibrational_hessian,
)
from scripts import d2c_hessian_symmetry_diagnostic as hessian_diagnostic
from scripts import d2c_input_bundle
from scripts import d2c_sct_campaign as campaign

# These resident references are part of the production boundary. Tests replace only
# the external PySCF constructor seam; production identity construction rejects
# replaced modules, aliases, helpers, or local orchestration code.
_RESIDENT_MODULES = {
    "quarry.pipeline": sys.modules[build_mol.__module__],
    "quarry.reaction_path": sys.modules[project_vibrational_hessian.__module__],
    hessian_diagnostic.__name__: hessian_diagnostic,
    d2c_input_bundle.__name__: d2c_input_bundle,
    campaign.__name__: campaign,
}
_RESIDENT_ALIASES = {
    "_make_scf": _make_scf,
    "build_mol": build_mol,
    "frequency_geometry_fingerprint": frequency_geometry_fingerprint,
    "hessian_eigenvalues_to_wavenumbers_cm": hessian_eigenvalues_to_wavenumbers_cm,
    "project_vibrational_hessian": project_vibrational_hessian,
}
_RESIDENT_ALIAS_CODE = {
    name: value.__code__ for name, value in _RESIDENT_ALIASES.items()
}
_RESIDENT_HELPERS = {
    (campaign.__name__, name): getattr(campaign, name)
    for name in (
        "_canonical_hash",
        "_code_digest",
        "_read_bounded_regular_snapshot",
        "_read_json_object",
        "_require_json_string",
        "_require_sha",
        "_safe_absolute_root",
        "_strict_json_equal",
    )
} | {
    (hessian_diagnostic.__name__, name): getattr(hessian_diagnostic, name)
    for name in (
        "_claimed_json_object",
        "_claimed_root_entry_exists",
        "_claimed_write",
        "_claimed_write_noreplace",
        "_exclusive_output_claim",
        "_load_half_step_source",
        "_read_bounded_regular_snapshot_at",
        "_terminal_write_noreplace",
        "_validated_completed_half_step_receipt",
    )
}
_RESIDENT_HELPER_CODE = {
    key: value.__code__ for key, value in _RESIDENT_HELPERS.items()
}
campaign._OBSERVED_MODULE_CODE.setdefault(__name__, _THIS_MODULE_IMPORT_CODE)

SCHEMA = "d2c-mode-projected-scalar-energy-diagnostic-v1"
ROUTE = "h-co-1w-oside"
SOURCE_ROOT = Path(
    "/mnt/data/vsletten/dissertation-data/task300-d2c-fd-halfstep-83d372d-r2"
)
SOURCE_RECEIPT_SHA256 = (
    "dbde8c707bddbb8f1d4a226755d4e52efbd15a437bd612255bc0e9c969d5f84a"
)
SOURCE_TERMINAL_SHA256 = (
    "39ef7b67ea4fc549e223eb2e1ac883aeda7c580726eb5d438d082b6c7e047a58"
)
SOURCE_STATUS_SHA256 = (
    "e5b827fc462021ef085df2e0b00a6e22bdd7678f11960e146662b544c898b64e"
)
SOURCE_MATRIX_RECEIPT_SHA256 = (
    "94ba29b302c09b62048213e6bfa327bca1ddb98f3abbfe8260af8458edfaef15"
)
SOURCE_MATRIX_SHA256 = (
    "65878bff5729f2c133efd5da705e6400689221e5e0754774d3f030d9a420619b"
)
SOURCE_GIT_SHA = "83d372d778431f2e733fec8b8ce64e5a5ae19b99"
BUNDLE_MANIFEST_SHA256 = (
    "2c0825ff79b9debc8240fd9cb15d27c91a6a4bd192302b857ede70308fa7f99d"
)
BUNDLE_FILE_COUNT = 29
BUNDLE_TOTAL_BYTES = 71934
SOURCE_CENTER_ENERGY_HARTREE = -190.15534901824873
SOURCE_CENTER_S2 = 0.7593962152748404
SOURCE_CENTER_RECEIPT_SHA256 = (
    "1aab510ab505c894f8c85ee3584d32285294ca16d79fc84c26cf488ab3d38cdf"
)
SOURCE_CENTER_DENSITY_SHA256 = (
    "e1d62f88e77e30074b093405d76439d2e2d04667c0a0c937777996d4de554ece"
)
MODE_HASHES = {
    "unstable": "20b2f5449e3410f8bd3b1dc6864cbf1d779a228211fb8b6281bc329ac9b3c5a1",
    "lowest_positive": (
        "cd5168f557f0eab8d1c54079b075801cd94036acc41d461ab9bbe17eda59ee2a"
    ),
}
MODE_SOURCE_EIGENVALUES = {
    "unstable": -0.017127293821654257,
    "lowest_positive": 2.5235155145404263e-06,
}
MODE_SOURCE_FREQUENCIES_CM = {
    "unstable": -672.7421694955478,
    "lowest_positive": 8.165960350269595,
}
MODE_ORDER = ("unstable", "lowest_positive")
Q_STEPS = (-0.10, -0.05, 0.05, 0.10)
POINT_ORDER = (
    "center",
    "unstable-m0.10",
    "unstable-m0.05",
    "unstable-p0.05",
    "unstable-p0.10",
    "lowest-positive-m0.10",
    "lowest-positive-m0.05",
    "lowest-positive-p0.05",
    "lowest-positive-p0.10",
)
CONTRACT: dict[str, Any] = {
    "route": ROUTE,
    "method": {
        "backend": "pyscf-cpu",
        "xc": "pwb6k",
        "dispersion": "d3bj",
        "basis": "def2-svp",
        "density_fit": True,
        "grid_level": 5,
        "required_grid_point_count": 199560,
        "scf_tolerance": 1.0e-12,
        "scf_max_cycle": 150,
        "scf_retry_count": 0,
        "charge": 0,
        "spin_2s": 1,
        "density_initial_guess": (
            "one newly converged center density used unchanged as dm0 for all eight "
            "independent displaced SCFs"
        ),
    },
    "sampling": {
        "coordinate_units": "bohr * sqrt(amu)",
        "mode_convention": (
            "unit Euclidean mass-weighted vector; x(q)=x(0)+q*v/sqrt(mass_amu)"
        ),
        "shared_center_count": 1,
        "modes": list(MODE_ORDER),
        "q_steps": list(Q_STEPS),
        "total_scf_single_point_energies": 9,
        "maximum_absolute_q": 0.10,
        "maximum_per_atom_cartesian_displacement_bohr_exclusive": 0.10,
        "optional_rotated_control": "omitted to preserve the exact nine-point bound",
    },
    "analysis": {
        "curvature_units": "hartree / (bohr^2 * amu)",
        "first_derivative_units": "hartree / (bohr * sqrt(amu))",
        "kappa_h": "(E(+h)-2*E(0)+E(-h))/h^2",
        "gradient_h": "(E(+h)-E(-h))/(2*h)",
        "richardson_curvature": "(4*kappa_0.05-kappa_0.10)/3",
        "richardson_first_derivative": "(4*g_0.05-g_0.10)/3",
        "estimated_stationary_q": "-g_R/kappa_R",
        "estimated_stationary_q_bound_rationale": (
            "abs(q*) < 0.05 Bohr*sqrt(amu) keeps the quadratic stationary estimate "
            "strictly inside the inner sampled stencil; it resolves only each "
            "projected direction"
        ),
        "stationarity_scope": (
            "two sampled directions only; full molecular stationarity unresolved"
        ),
    },
    "thresholds": {
        "spin_square_inclusive_minimum": 0.74,
        "spin_square_inclusive_maximum": 0.80,
        "displaced_vs_center_spin_square_exclusive_maximum_delta": 0.01,
        (
            "displaced_final_density_vs_center_exclusive_maximum_"
            "normalized_frobenius_delta"
        ): 0.05,
        "center_energy_vs_source_exclusive_maximum_delta_hartree": 1.0e-8,
        "center_s2_vs_source_exclusive_maximum_delta": 1.0e-4,
        "curvature_convergence": (
            "abs(kappa_R-kappa_0.05) < max(1e-6,0.005*abs(kappa_R))"
        ),
        "curvature_convergence_absolute_floor_exclusive": 1.0e-6,
        "curvature_convergence_relative_fraction_exclusive": 0.005,
        "unstable_curvature_exclusive_maximum": -1.0e-8,
        "unstable_absolute_frequency_cm_exclusive_minimum": 200.0,
        "lowest_positive_curvature_exclusive_minimum": 5.0e-6,
        "estimated_stationary_absolute_q_exclusive_maximum": 0.05,
    },
    "source": {
        "root": str(SOURCE_ROOT),
        "receipt_sha256": SOURCE_RECEIPT_SHA256,
        "terminal_sha256": SOURCE_TERMINAL_SHA256,
        "status_sha256": SOURCE_STATUS_SHA256,
        "matrix_receipt_sha256": SOURCE_MATRIX_RECEIPT_SHA256,
        "richardson_symmetric_matrix_sha256": SOURCE_MATRIX_SHA256,
        "git_sha": SOURCE_GIT_SHA,
        "mode_hashes": MODE_HASHES,
        "mode_source_eigenvalues": MODE_SOURCE_EIGENVALUES,
        "mode_source_frequencies_cm": MODE_SOURCE_FREQUENCIES_CM,
        "bundle_manifest_sha256": BUNDLE_MANIFEST_SHA256,
        "bundle_file_count": BUNDLE_FILE_COUNT,
        "bundle_total_bytes": BUNDLE_TOTAL_BYTES,
        "center_energy_hartree": SOURCE_CENTER_ENERGY_HARTREE,
        "center_s2": SOURCE_CENTER_S2,
        "center_receipt_sha256": SOURCE_CENTER_RECEIPT_SHA256,
        "center_density_sha256": SOURCE_CENTER_DENSITY_SHA256,
    },
    "policy": {
        "evidence_only_forever": True,
        "accepted_campaign_result": False,
        "transition_state_qualification": False,
        "irc_authority": False,
        "sct_authority": False,
        "d3b_authority": False,
        "full_stationarity_resolved": False,
        "center_gradient_computed": False,
        "center_fmax_gate": "omitted because no gradient is evaluated",
        "failed_or_partial_root_resume": "forbidden; create a fresh output root",
        "robust_root_continuity": (
            "unavailable beyond predeclared S2 and density continuity proxies; each "
            "displaced robust-root gate therefore fails closed"
        ),
    },
}


class ScientificRejection(ValueError):
    """The complete scalar evidence failed one or more predeclared gates."""


def _json_bytes(payload: dict[str, Any]) -> bytes:
    return (
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode()


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _finite_float(value: Any, *, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be a finite number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def _now_utc() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _utc_timestamp(value: Any, *, label: str) -> datetime:
    if type(value) is not str or not value.endswith("Z"):
        raise ValueError(f"{label} must be a canonical UTC timestamp")
    try:
        parsed = datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=UTC)
    except ValueError as exc:
        raise ValueError(f"{label} must be a canonical UTC timestamp") from exc
    if parsed.isoformat().replace("+00:00", "Z") != value:
        raise ValueError(f"{label} must be a canonical UTC timestamp")
    return parsed


def _canonical_mode(vector: Any) -> np.ndarray:
    mode = np.asarray(vector, dtype=float).reshape(-1)
    if mode.shape != (18,) or not np.all(np.isfinite(mode)):
        raise ValueError("source mode must contain exactly 18 finite values")
    norm = float(np.linalg.norm(mode))
    if not math.isfinite(norm) or abs(norm - 1.0) > 2.0e-14:
        raise ValueError("source mass-weighted mode is not normalized")
    pivot = int(np.argmax(np.abs(mode)))
    if mode[pivot] < 0.0:
        mode = -mode
    return np.frombuffer(np.ascontiguousarray(mode, dtype="<f8").tobytes(), dtype="<f8")


def _mode_sha(mode: np.ndarray) -> str:
    return _sha(np.ascontiguousarray(mode, dtype="<f8").tobytes())


def _regular_tree_inventory(root: Path) -> tuple[int, int]:
    """Inventory the pinned bundle without unbounded traversal or reads."""

    nofollow = getattr(os, "O_NOFOLLOW", None)
    if nofollow is None:
        raise RuntimeError("O_NOFOLLOW is required for the tracked bundle inventory")
    root_flags = os.O_RDONLY | os.O_DIRECTORY | nofollow | getattr(os, "O_CLOEXEC", 0)
    root_fd = os.open(root, root_flags)
    count = 0
    total = 0
    directory_count = 0
    pending: list[int] = [root_fd]
    try:
        while pending:
            directory_fd = pending.pop()
            try:
                directory_count += 1
                if directory_count > BUNDLE_FILE_COUNT + 1:
                    raise ValueError(
                        "tracked D2c route bundle has too many directories"
                    )
                with os.scandir(directory_fd) as entries:
                    for entry in entries:
                        if entry.is_symlink():
                            raise ValueError(
                                "tracked D2c route bundle cannot contain symlinks"
                            )
                        if entry.is_dir(follow_symlinks=False):
                            pending.append(
                                os.open(entry.name, root_flags, dir_fd=directory_fd)
                            )
                            continue
                        if not entry.is_file(follow_symlinks=False):
                            raise ValueError(
                                "tracked D2c route bundle member is not regular"
                            )
                        if count >= BUNDLE_FILE_COUNT:
                            raise ValueError(
                                "tracked D2c route bundle exceeds 29 files"
                            )
                        descriptor = os.open(
                            entry.name,
                            os.O_RDONLY | nofollow | getattr(os, "O_CLOEXEC", 0),
                            dir_fd=directory_fd,
                        )
                        try:
                            before = os.fstat(descriptor)
                            remaining_budget = BUNDLE_TOTAL_BYTES - total
                            if (
                                not stat.S_ISREG(before.st_mode)
                                or before.st_size <= 0
                                or before.st_size > remaining_budget
                            ):
                                raise ValueError(
                                    "tracked D2c route bundle exceeds its byte budget"
                                )
                            digest = hashlib.sha256()
                            remaining = before.st_size
                            while remaining:
                                chunk = os.read(descriptor, min(remaining, 64 * 1024))
                                if not chunk:
                                    raise ValueError(
                                        "tracked D2c route bundle member changed "
                                        "while read"
                                    )
                                digest.update(chunk)
                                remaining -= len(chunk)
                            if os.read(descriptor, 1):
                                raise ValueError(
                                    "tracked D2c route bundle member grew while read"
                                )
                            after = os.fstat(descriptor)
                            if (
                                before.st_dev,
                                before.st_ino,
                                before.st_size,
                                before.st_mtime_ns,
                                before.st_ctime_ns,
                            ) != (
                                after.st_dev,
                                after.st_ino,
                                after.st_size,
                                after.st_mtime_ns,
                                after.st_ctime_ns,
                            ):
                                raise ValueError(
                                    "tracked D2c route bundle member changed while read"
                                )
                        finally:
                            os.close(descriptor)
                        count += 1
                        total += before.st_size
            finally:
                os.close(directory_fd)
        return count, total
    except BaseException:
        for descriptor in pending:
            os.close(descriptor)
        raise


def _canonical_json_at(
    path: Path, *, label: str, maximum_bytes: int
) -> tuple[dict[str, Any], bytes]:
    payload, raw = campaign._read_json_object(
        path, label=label, max_bytes=maximum_bytes
    )
    if _json_bytes(payload) != raw:
        raise ValueError(f"{label} is not strict canonical JSON")
    return payload, raw


def _read_bundle_relative_at(
    root_fd: int, relative: Path, *, label: str, maximum_bytes: int
) -> bytes:
    nofollow = getattr(os, "O_NOFOLLOW", None)
    if nofollow is None:
        raise RuntimeError("O_NOFOLLOW is required for bundle verification")
    if (
        relative.is_absolute()
        or not relative.parts
        or any(part in {"", ".", ".."} for part in relative.parts)
    ):
        raise ValueError(f"{label} has an unsafe relative path")
    parent_fd = os.dup(root_fd)
    try:
        for part in relative.parts[:-1]:
            child_fd = os.open(
                part,
                os.O_RDONLY | os.O_DIRECTORY | nofollow | getattr(os, "O_CLOEXEC", 0),
                dir_fd=parent_fd,
            )
            os.close(parent_fd)
            parent_fd = child_fd
        return hessian_diagnostic._read_bounded_regular_snapshot_at(
            parent_fd,
            relative.parts[-1],
            label=label,
            maximum_bytes=maximum_bytes,
        )
    except OSError as exc:
        raise ValueError(f"{label} must remain below the bundle root") from exc
    finally:
        os.close(parent_fd)


def _verify_bounded_bundle(root: Path) -> tuple[dict[str, Any], bytes]:
    """Verify the exact manifest inventory through one pinned root descriptor."""

    nofollow = getattr(os, "O_NOFOLLOW", None)
    if nofollow is None:
        raise RuntimeError("O_NOFOLLOW is required for bundle verification")
    root_fd = os.open(
        root,
        os.O_RDONLY | os.O_DIRECTORY | nofollow | getattr(os, "O_CLOEXEC", 0),
    )
    before = os.fstat(root_fd)
    try:
        manifest_raw = _read_bundle_relative_at(
            root_fd,
            Path("manifest.json"),
            label="tracked D2c route bundle manifest",
            maximum_bytes=min(BUNDLE_TOTAL_BYTES, 1024 * 1024),
        )
        try:
            manifest = json.loads(manifest_raw)
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise ValueError(
                "tracked D2c route bundle manifest is invalid JSON"
            ) from exc
        if type(manifest) is not dict or _json_bytes(manifest) != manifest_raw:
            raise ValueError(
                "tracked D2c route bundle manifest is not strict canonical JSON"
            )
        if manifest.get("schema") != d2c_input_bundle.SCHEMA or manifest.get(
            "selection", {}
        ).get("routes") != list(d2c_input_bundle.ROUTES):
            raise ValueError("tracked D2c route bundle manifest contract drifted")
        routes = manifest.get("routes")
        if type(routes) is not dict or set(routes) != set(d2c_input_bundle.ROUTES):
            raise ValueError("tracked D2c route bundle route inventory drifted")
        expected_paths = {Path("manifest.json")}
        for route in d2c_input_bundle.ROUTES:
            route_record = routes[route]
            files = route_record.get("files") if type(route_record) is dict else None
            if type(files) is not dict or set(files) != set(d2c_input_bundle.FILES):
                raise ValueError(f"tracked D2c route file inventory drifted: {route}")
            for name in d2c_input_bundle.FILES:
                relative = Path(route) / name
                expected_paths.add(relative)
                record = files[name]
                if (
                    type(record) is not dict
                    or set(record) != {"bytes", "sha256"}
                    or type(record["bytes"]) is not int
                    or not 1 <= record["bytes"] <= BUNDLE_TOTAL_BYTES
                ):
                    raise ValueError(f"tracked D2c file receipt is invalid: {relative}")
                campaign._require_sha(
                    campaign._require_json_string(
                        record["sha256"], label=f"{relative} SHA-256"
                    ),
                    length=64,
                    label=f"{relative} SHA-256",
                )
                raw = _read_bundle_relative_at(
                    root_fd,
                    relative,
                    label=f"tracked D2c bundled input {relative}",
                    maximum_bytes=record["bytes"],
                )
                if len(raw) != record["bytes"] or _sha(raw) != record["sha256"]:
                    raise ValueError(f"tracked D2c bundled input drifted: {relative}")
        if len(expected_paths) != BUNDLE_FILE_COUNT:
            raise ValueError(
                "tracked D2c expected bundle inventory is not exactly 29 files"
            )
        after = os.fstat(root_fd)
        published = os.stat(root, follow_symlinks=False)
        if (
            not stat.S_ISDIR(published.st_mode)
            or (
                before.st_dev,
                before.st_ino,
                before.st_mtime_ns,
                before.st_ctime_ns,
            )
            != (after.st_dev, after.st_ino, after.st_mtime_ns, after.st_ctime_ns)
            or (
                after.st_dev,
                after.st_ino,
                after.st_mtime_ns,
                after.st_ctime_ns,
            )
            != (
                published.st_dev,
                published.st_ino,
                published.st_mtime_ns,
                published.st_ctime_ns,
            )
        ):
            raise RuntimeError(
                "tracked D2c route bundle root changed during verification"
            )
        return manifest, manifest_raw
    finally:
        os.close(root_fd)


def _load_source_evidence(source_root: Path = SOURCE_ROOT) -> dict[str, Any]:
    """Validate the complete negative analytic/FD/half-step ancestry and modes."""

    root = campaign._safe_absolute_root(source_root)
    if root != SOURCE_ROOT:
        raise ValueError("scalar diagnostic requires the exact pinned half-step root")
    receipt, receipt_sha = hessian_diagnostic._validated_completed_half_step_receipt(
        root, verify_resident_identity=False
    )
    if receipt_sha != SOURCE_RECEIPT_SHA256:
        raise ValueError("half-step source receipt SHA-256 mismatch")
    if (
        receipt.get("git_sha") != SOURCE_GIT_SHA
        or receipt.get("confirmation_passed") is not False
        or receipt.get("accepted_campaign_result") is not False
        or receipt.get("route") != ROUTE
    ):
        raise ValueError("half-step source is not the pinned negative internal verdict")
    terminal, terminal_raw = _canonical_json_at(
        root / "terminal.json", label="half-step terminal", maximum_bytes=1024 * 1024
    )
    status, status_raw = _canonical_json_at(
        root / "status.json", label="half-step status", maximum_bytes=4 * 1024 * 1024
    )
    if (
        _sha(terminal_raw) != SOURCE_TERMINAL_SHA256
        or _sha(status_raw) != SOURCE_STATUS_SHA256
    ):
        raise ValueError("half-step terminal/status SHA-256 mismatch")
    for label, payload in (("terminal", terminal), ("status", status)):
        if (
            payload.get("state") != "completed"
            or payload.get("accepted_campaign_result") is not False
            or payload.get("confirmation_passed") is not False
        ):
            raise ValueError(f"half-step {label} lost its negative internal verdict")

    source = hessian_diagnostic._load_half_step_source(
        Path(hessian_diagnostic.PRIOR_FD_ROOT)
    )
    prior_fd = source["prior_receipt"]
    analytic = source["reference_receipt"]
    if (
        prior_fd.get("confirmation_passed") is not False
        or prior_fd.get("accepted_campaign_result") is not False
    ):
        raise ValueError("prior FD ancestry is not the required negative verdict")
    center_source, center_source_raw = _canonical_json_at(
        Path(hessian_diagnostic.PRIOR_FD_ROOT) / "points" / "center" / "receipt.json",
        label="authoritative prior center receipt",
        maximum_bytes=64 * 1024,
    )
    if (
        _sha(center_source_raw) != SOURCE_CENTER_RECEIPT_SHA256
        or center_source.get("point_key") != "center"
        or center_source.get("scf_converged") is not True
        or center_source.get("backend") != "pyscf-cpu"
        or center_source.get("electronic_hartree") != SOURCE_CENTER_ENERGY_HARTREE
        or center_source.get("spin_square", [None])[0] != SOURCE_CENTER_S2
        or center_source.get("density_final_sha256") != SOURCE_CENTER_DENSITY_SHA256
    ):
        raise ValueError("authoritative prior center energy/S2/density receipt drifted")
    if analytic.get("accepted_campaign_result") is not False:
        raise ValueError("analytic ancestry is not evidence-only")
    analytic_cases = analytic.get("cases")
    if type(analytic_cases) is not list or len(analytic_cases) != 4:
        raise ValueError("analytic ancestry case inventory drifted")
    analytic_internal_verdicts = []
    for case in analytic_cases:
        if type(case) is not dict or type(case.get("components")) is not dict:
            raise ValueError("analytic ancestry case verdict is malformed")
        electronic = case["components"].get("electronic", {}).get("symmetry", {})
        total = case["components"].get("total", {}).get("symmetry", {})
        if (
            electronic.get("accepted") is not False
            or total.get("accepted") is not False
        ):
            raise ValueError(
                "analytic ancestry lost its negative internal symmetry verdict"
            )
        analytic_internal_verdicts.append(
            {
                "case": case.get("name"),
                "electronic_symmetry_accepted": False,
                "total_symmetry_accepted": False,
            }
        )

    matrix_receipt, matrix_receipt_raw = _canonical_json_at(
        root / "matrices" / "receipt.json",
        label="half-step matrix receipt",
        maximum_bytes=4 * 1024 * 1024,
    )
    if _sha(matrix_receipt_raw) != SOURCE_MATRIX_RECEIPT_SHA256:
        raise ValueError("half-step matrix receipt SHA-256 mismatch")
    artifact = matrix_receipt.get("artifacts", {}).get("richardson-symmetric.f64")
    expected_artifact = {
        "dtype": "little-endian float64",
        "path": "matrices/richardson-symmetric.f64",
        "shape": [18, 18],
        "sha256": SOURCE_MATRIX_SHA256,
        "units": "hartree / bohr^2",
    }
    if artifact != expected_artifact:
        raise ValueError("Richardson matrix artifact metadata drifted")
    matrix_raw = campaign._read_bounded_regular_snapshot(
        root / artifact["path"],
        label="half-step Richardson matrix",
        maximum_bytes=18 * 18 * 8,
    )
    if len(matrix_raw) != 18 * 18 * 8 or _sha(matrix_raw) != SOURCE_MATRIX_SHA256:
        raise ValueError("half-step Richardson matrix bytes drifted")
    matrix = np.frombuffer(matrix_raw, dtype="<f8").reshape(18, 18)
    if not np.all(np.isfinite(matrix)) or not np.array_equal(matrix, matrix.T):
        raise ValueError("half-step Richardson matrix is not finite and symmetric")

    bundle_root = campaign._safe_absolute_root(campaign.DEFAULT_BUNDLE_ROOT)
    bundle_count, bundle_bytes = _regular_tree_inventory(bundle_root)
    _, manifest_raw = _verify_bounded_bundle(bundle_root)
    final_bundle_count, final_bundle_bytes = _regular_tree_inventory(bundle_root)
    if (
        _sha(manifest_raw) != BUNDLE_MANIFEST_SHA256
        or bundle_count != BUNDLE_FILE_COUNT
        or bundle_bytes != BUNDLE_TOTAL_BYTES
        or final_bundle_count != bundle_count
        or final_bundle_bytes != bundle_bytes
    ):
        raise ValueError("tracked D2c route bundle byte inventory drifted")

    modes = project_vibrational_hessian(
        source["transition_state"].coords, source["masses"], matrix
    )
    frequencies = hessian_eigenvalues_to_wavenumbers_cm(modes.eigenvalues)
    selected = {}
    for index, name in enumerate(MODE_ORDER):
        mode = _canonical_mode(modes.mass_weighted_eigenvectors[index])
        eigenvalue = _finite_float(modes.eigenvalues[index], label=f"{name} eigenvalue")
        frequency = _finite_float(frequencies[index], label=f"{name} frequency")
        if (
            _mode_sha(mode) != MODE_HASHES[name]
            or eigenvalue != MODE_SOURCE_EIGENVALUES[name]
            or frequency != MODE_SOURCE_FREQUENCIES_CM[name]
        ):
            raise ValueError(f"pinned {name} mode vector or spectrum drifted")
        selected[name] = {
            "index": index,
            "mass_weighted_vector": mode.tolist(),
            "mass_weighted_vector_sha256": MODE_HASHES[name],
            "source_eigenvalue_hartree_per_bohr2_amu": eigenvalue,
            "source_signed_wavenumber_cm": frequency,
        }

    masses = np.asarray(source["masses"], dtype=float)
    if masses.shape != (6,) or not np.all(np.isfinite(masses)) or np.any(masses <= 0.0):
        raise ValueError("source isotopic masses are invalid")
    masses_raw = np.ascontiguousarray(masses, dtype="<f8").tobytes()
    return {
        "ancestry": {
            "analytic": {
                "receipt_sha256": hessian_diagnostic.REFERENCE_RECEIPT_SHA256,
                "accepted_campaign_result": False,
                "internal_symmetry_policy_passed": False,
                "internal_verdicts": analytic_internal_verdicts,
            },
            "finite_difference": {
                "receipt_sha256": hessian_diagnostic.PRIOR_FD_RECEIPT_SHA256,
                "terminal_sha256": hessian_diagnostic.PRIOR_FD_TERMINAL_SHA256,
                "status_sha256": hessian_diagnostic.PRIOR_FD_STATUS_SHA256,
                "confirmation_passed": False,
                "accepted_campaign_result": False,
            },
            "half_step": {
                "receipt_sha256": SOURCE_RECEIPT_SHA256,
                "terminal_sha256": SOURCE_TERMINAL_SHA256,
                "status_sha256": SOURCE_STATUS_SHA256,
                "matrix_receipt_sha256": SOURCE_MATRIX_RECEIPT_SHA256,
                "richardson_symmetric_matrix_sha256": SOURCE_MATRIX_SHA256,
                "confirmation_passed": False,
                "accepted_campaign_result": False,
            },
        },
        "bundle": {
            "root": str(bundle_root),
            "manifest_sha256": BUNDLE_MANIFEST_SHA256,
            "file_count": bundle_count,
            "total_bytes": bundle_bytes,
        },
        "reference_binding": receipt["reference_binding"],
        "reference_binding_sha256": receipt["reference_binding_sha256"],
        "authoritative_center": {
            "receipt_sha256": SOURCE_CENTER_RECEIPT_SHA256,
            "electronic_hartree": SOURCE_CENTER_ENERGY_HARTREE,
            "spin_square": SOURCE_CENTER_S2,
            "density_sha256": SOURCE_CENTER_DENSITY_SHA256,
            "energy_exclusive_tolerance_hartree": CONTRACT["thresholds"][
                "center_energy_vs_source_exclusive_maximum_delta_hartree"
            ],
            "s2_exclusive_tolerance": CONTRACT["thresholds"][
                "center_s2_vs_source_exclusive_maximum_delta"
            ],
        },
        "transition_state_geometry_fingerprint": frequency_geometry_fingerprint(
            source["transition_state"]
        ),
        "masses_amu": masses.tolist(),
        "masses_sha256": _sha(masses_raw),
        "modes": selected,
        "transition_state": source["transition_state"],
        "settings": source["settings"],
    }


def _code_sha(code: CodeType) -> str:
    return campaign._code_digest(code)


def _nested_code_sha256(code: CodeType) -> set[str]:
    values = {_code_sha(code)}
    for constant in code.co_consts:
        if isinstance(constant, CodeType):
            values.update(_nested_code_sha256(constant))
    return values


def _assert_resident_bindings() -> None:
    for module_name, module in _RESIDENT_MODULES.items():
        if sys.modules.get(module_name) is not module:
            raise RuntimeError(f"critical resident module was replaced: {module_name}")
    for name, expected in _RESIDENT_ALIASES.items():
        observed = globals().get(name)
        if (
            observed is not expected
            or observed.__code__ is not _RESIDENT_ALIAS_CODE[name]
        ):
            raise RuntimeError(f"critical resident callable was replaced: {name}")
    for (module_name, name), expected in _RESIDENT_HELPERS.items():
        module = _RESIDENT_MODULES[module_name]
        observed = getattr(module, name, None)
        if (
            observed is not expected
            or observed.__code__ is not _RESIDENT_HELPER_CODE[(module_name, name)]
        ):
            raise RuntimeError(
                f"critical resident helper callable was replaced: {module_name}.{name}"
            )
    if (
        globals().get("hessian_diagnostic")
        is not _RESIDENT_MODULES[hessian_diagnostic.__name__]
        or globals().get("d2c_input_bundle")
        is not _RESIDENT_MODULES[d2c_input_bundle.__name__]
        or globals().get("campaign") is not _RESIDENT_MODULES[campaign.__name__]
    ):
        raise RuntimeError("critical resident module alias was replaced")
    for name, expected in globals().get("_RESIDENT_LOCAL_CALLABLES", {}).items():
        observed = globals().get(name)
        if (
            observed is not expected
            or observed.__code__ is not _RESIDENT_LOCAL_CODE[name]
        ):
            raise RuntimeError(f"critical local callable was replaced: {name}")


def _current_module_execution_record() -> dict[str, Any]:
    source = Path(__file__)
    if source.is_symlink():
        raise RuntimeError("resident scalar diagnostic source cannot be a symlink")
    source = source.resolve(strict=True)
    raw = campaign._read_bounded_regular_snapshot(
        source, label="resident scalar diagnostic source", maximum_bytes=4 * 1024 * 1024
    )
    compiled = compile(
        raw.decode("utf-8"),
        str(source),
        "exec",
        dont_inherit=True,
        optimize=sys.flags.optimize,
    )
    if campaign._OBSERVED_MODULE_CODE.get(
        __name__
    ) is not _THIS_MODULE_IMPORT_CODE or _code_sha(
        _THIS_MODULE_IMPORT_CODE
    ) != _code_sha(compiled):
        raise RuntimeError(
            "scalar diagnostic resident module code disagrees with source"
        )
    source_code = _nested_code_sha256(compiled)
    callables = {}
    for name, value in sorted(_RESIDENT_LOCAL_CALLABLES.items()):
        code_sha = _code_sha(value.__code__)
        if code_sha not in source_code:
            raise RuntimeError(
                f"resident scalar callable disagrees with source: {name}"
            )
        callables[name] = {
            "qualname": value.__qualname__,
            "code_sha256": code_sha,
            "signature": str(inspect.signature(value)),
        }
    return {
        "origin": str(source),
        "sha256": _sha(raw),
        "byte_count": len(raw),
        "loaded_module_code_sha256": _code_sha(_THIS_MODULE_IMPORT_CODE),
        "callables": callables,
    }


def _git_source_snapshot() -> dict[str, Any]:
    source = Path(__file__)
    if source.is_symlink():
        raise RuntimeError("scalar diagnostic source cannot be a symlink")
    source = source.resolve(strict=True)
    repository = source.parents[2]
    status = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=all"],
        cwd=repository,
        check=True,
        capture_output=True,
        text=True,
    )
    if status.stdout:
        raise RuntimeError("scalar diagnostic requires a clean Git worktree")
    git_sha = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repository,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    git_tree = subprocess.run(
        ["git", "rev-parse", "HEAD^{tree}"],
        cwd=repository,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    campaign._require_sha(git_sha, length=40, label="scalar source Git commit")
    campaign._require_sha(git_tree, length=40, label="scalar source Git tree")
    source_paths = {source}
    for module in _RESIDENT_MODULES.values():
        module_path = Path(module.__file__)
        if module_path.is_symlink():
            raise RuntimeError(
                f"critical source cannot be a symlink: {module.__name__}"
            )
        source_paths.add(module_path.resolve(strict=True))
    source_files: dict[str, str] = {}
    for path in sorted(source_paths, key=str):
        relative = path.relative_to(repository)
        raw = campaign._read_bounded_regular_snapshot(
            path, label=f"critical source {relative}", maximum_bytes=4 * 1024 * 1024
        )
        committed = subprocess.run(
            ["git", "show", f"{git_sha}:{relative.as_posix()}"],
            cwd=repository,
            check=True,
            capture_output=True,
        ).stdout
        if committed != raw:
            raise RuntimeError(f"critical source is not committed: {relative}")
        source_files[relative.as_posix()] = _sha(raw)
    return {
        "git_commit": git_sha,
        "git_tree": git_tree,
        "source_files": source_files,
    }


def _module_execution_manifest() -> dict[str, Any]:
    records: dict[str, Any] = {}
    callables_by_module: dict[str, dict[str, Any]] = {}
    for name, value in _RESIDENT_ALIASES.items():
        callables_by_module.setdefault(value.__module__, {})[name] = value
    for (module_name, name), value in _RESIDENT_HELPERS.items():
        callables_by_module.setdefault(module_name, {})[name] = value
    for module_name, module in _RESIDENT_MODULES.items():
        path = Path(module.__file__)
        if path.is_symlink():
            raise RuntimeError(
                f"critical executable module is a symlink: {module_name}"
            )
        path = path.resolve(strict=True)
        raw = campaign._read_bounded_regular_snapshot(
            path,
            label=f"critical executable module {module_name}",
            maximum_bytes=4 * 1024 * 1024,
        )
        compiled = compile(
            raw.decode("utf-8"),
            str(path),
            "exec",
            dont_inherit=True,
            optimize=sys.flags.optimize,
        )
        source_code = _nested_code_sha256(compiled)
        resident_module_code = campaign._OBSERVED_MODULE_CODE.get(module_name)
        if resident_module_code is None:
            candidate = getattr(module, "_THIS_MODULE_IMPORT_CODE", None)
            if isinstance(candidate, CodeType):
                resident_module_code = candidate
        loaded_module_code_sha = None
        if resident_module_code is not None:
            loaded_module_code_sha = _code_sha(resident_module_code)
            if loaded_module_code_sha != _code_sha(compiled):
                raise RuntimeError(
                    f"resident module-level code disagrees with source: {module_name}"
                )
        callable_records = {}
        source_cache = {path: source_code}
        source_hashes = {path: _sha(raw)}
        for alias, value in sorted(callables_by_module.get(module_name, {}).items()):
            candidates = [("resident", value)]
            wrapped = getattr(value, "__wrapped__", None)
            if inspect.isfunction(wrapped):
                candidates.append(("wrapped", wrapped))
            code_records = {}
            for kind, candidate in candidates:
                code_origin = Path(candidate.__code__.co_filename)
                if code_origin.is_symlink():
                    raise RuntimeError(
                        f"resident callable source is a symlink: {module_name}.{alias}"
                    )
                code_origin = code_origin.resolve(strict=True)
                if code_origin not in source_cache:
                    code_raw = campaign._read_bounded_regular_snapshot(
                        code_origin,
                        label=f"callable source {module_name}.{alias}",
                        maximum_bytes=4 * 1024 * 1024,
                    )
                    code_compiled = compile(
                        code_raw.decode("utf-8"),
                        str(code_origin),
                        "exec",
                        dont_inherit=True,
                        optimize=sys.flags.optimize,
                    )
                    source_cache[code_origin] = _nested_code_sha256(code_compiled)
                    source_hashes[code_origin] = _sha(code_raw)
                code_sha = _code_sha(candidate.__code__)
                if code_sha not in source_cache[code_origin]:
                    raise RuntimeError(
                        "resident callable disagrees with on-disk source: "
                        f"{module_name}.{alias}"
                    )
                code_records[kind] = {
                    "qualname": candidate.__qualname__,
                    "code_sha256": code_sha,
                    "origin": str(code_origin),
                    "source_sha256": source_hashes[code_origin],
                    "signature": str(inspect.signature(candidate)),
                }
            callable_records[alias] = code_records
        records[module_name] = {
            "origin": str(path),
            "sha256": _sha(raw),
            "byte_count": len(raw),
            "loaded_module_code_sha256": loaded_module_code_sha,
            "callables": callable_records,
        }
    records[__name__] = _current_module_execution_record()
    return dict(sorted(records.items()))


def _native_payload_manifest() -> dict[str, Any]:
    for module in _PRELOADED_CPU_MODULES:
        if sys.modules.get(module.__name__) is not module:
            raise RuntimeError(
                f"preloaded CPU dependency module was replaced: {module.__name__}"
            )
    distributions = tuple(
        importlib.metadata.distribution(name)
        for name in ("numpy", "pyscf", "pyscf-dispersion", "scipy")
    )
    allowed: set[Path] = set()
    for distribution in distributions:
        for relative in distribution.files or ():
            name = str(relative).lower()
            if not (name.endswith((".so", ".dylib", ".dll", ".pyd")) or ".so." in name):
                continue
            path = Path(distribution.locate_file(relative))
            if path.exists():
                allowed.add(path.resolve(strict=True))
    maps_fd = os.open(
        "/proc/self/maps",
        os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        maps_chunks: list[bytes] = []
        maps_bytes = 0
        while chunk := os.read(maps_fd, 64 * 1024):
            maps_bytes += len(chunk)
            if maps_bytes > 16 * 1024 * 1024:
                raise RuntimeError(
                    "CPU native executable mapping inventory is oversized"
                )
            maps_chunks.append(chunk)
        maps_raw = b"".join(maps_chunks)
    finally:
        os.close(maps_fd)
    if not maps_raw:
        raise RuntimeError("CPU native executable mapping inventory is empty")
    mapped: set[Path] = set()
    for line in maps_raw.decode("utf-8").splitlines():
        fields = line.split(maxsplit=5)
        if len(fields) != 6 or "x" not in fields[1] or not fields[5].startswith("/"):
            continue
        raw_path = fields[5]
        if raw_path.endswith(" (deleted)"):
            candidate = Path(raw_path.removesuffix(" (deleted)"))
            if candidate in allowed:
                raise RuntimeError(f"CPU native payload was deleted: {candidate}")
            continue
        path = Path(raw_path).resolve(strict=True)
        if path in allowed:
            mapped.add(path)
    records: dict[str, Any] = {}
    environment = Path(sys.prefix).resolve(strict=True)
    for path in sorted(mapped, key=str):
        if path.is_symlink():
            raise RuntimeError(f"CPU native payload cannot be a symlink: {path}")
        raw = campaign._read_bounded_regular_snapshot(
            path, label=f"CPU native payload {path}", maximum_bytes=128 * 1024 * 1024
        )
        try:
            key = path.relative_to(environment).as_posix()
        except ValueError:
            key = str(path)
        records[key] = {
            "origin": str(path),
            "sha256": _sha(raw),
            "byte_count": len(raw),
        }
    if not records:
        raise RuntimeError("CPU scientific stack exposed no native payloads")
    return records


def _cpu_dependency_versions() -> dict[str, str]:
    return {
        distribution: importlib.metadata.version(distribution)
        for distribution in ("numpy", "pyscf", "pyscf-dispersion", "scipy")
    }


def _python_executable_manifest() -> dict[str, Any]:
    path = Path(sys.executable).resolve(strict=True)
    raw = campaign._read_bounded_regular_snapshot(
        path, label="CPU Python executable", maximum_bytes=128 * 1024 * 1024
    )
    return {"origin": str(path), "sha256": _sha(raw), "byte_count": len(raw)}


def _current_identity() -> dict[str, Any]:
    """Attest the exact CPU-only resident execution boundary."""

    _assert_resident_bindings()
    source_snapshot = _git_source_snapshot()
    identity = {
        "source_snapshot": source_snapshot,
        "python": platform.python_version(),
        "dependencies": _cpu_dependency_versions(),
        "python_executable": _python_executable_manifest(),
        "executable_modules": _module_execution_manifest(),
        "native_payloads": _native_payload_manifest(),
        "accelerator": "none",
        "cpu_only": True,
    }
    encoded = json.dumps(identity, sort_keys=True).lower()
    if "cuda" in encoded or "gpu" in encoded:
        raise RuntimeError(
            "CPU execution identity contains a forbidden accelerator flag"
        )
    _assert_resident_bindings()
    if _git_source_snapshot() != source_snapshot:
        raise RuntimeError("Git source snapshot changed during execution attestation")
    return identity


def _public_evidence(evidence: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in evidence.items()
        if key not in {"transition_state", "settings"}
    }


def _preflight_payload(
    identity: dict[str, Any], evidence: dict[str, Any]
) -> dict[str, Any]:
    payload = {
        "schema": SCHEMA,
        "state": "preflight-completed",
        "accepted_campaign_result": False,
        "created_utc": _now_utc(),
        "contract": CONTRACT,
        "contract_sha256": campaign._canonical_hash(CONTRACT),
        "execution_identity": identity,
        "execution_identity_sha256": campaign._canonical_hash(identity),
        "source_evidence": _public_evidence(evidence),
    }
    payload["preflight_identity"] = campaign._canonical_hash(payload)
    return payload


def create_preflight(source_root: Path, output_root: Path) -> dict[str, Any]:
    if (
        _assert_resident_bindings is not _ATTESTED_RESIDENT_ASSERT
        or _assert_resident_bindings.__code__ is not _ATTESTED_RESIDENT_ASSERT_CODE
    ):
        raise RuntimeError("resident-boundary validator was replaced")
    _ATTESTED_RESIDENT_ASSERT()
    source = campaign._safe_absolute_root(source_root)
    output = campaign._safe_absolute_root(output_root)
    with hessian_diagnostic._exclusive_output_claim(output) as claim:
        if not claim.root_created:
            raise FileExistsError(
                "scalar diagnostic preflight requires a fresh output root"
            )
        identity = _current_identity()
        evidence = _load_source_evidence(source)
        if _current_identity() != identity:
            raise RuntimeError("execution identity changed during scalar preflight")
        preflight = _preflight_payload(identity, evidence)
        hessian_diagnostic._claimed_write_noreplace(
            claim, output / "preflight.json", _json_bytes(preflight)
        )
        return preflight


def _validated_preflight(
    claim: hessian_diagnostic._OutputRootClaim,
) -> tuple[dict[str, Any], bytes]:
    preflight, raw = hessian_diagnostic._claimed_json_object(
        claim,
        claim.path / "preflight.json",
        label="scalar diagnostic preflight",
        max_bytes=4 * 1024 * 1024,
    )
    expected_keys = {
        "schema",
        "state",
        "accepted_campaign_result",
        "created_utc",
        "contract",
        "contract_sha256",
        "execution_identity",
        "execution_identity_sha256",
        "source_evidence",
        "preflight_identity",
    }
    identity = preflight.get("execution_identity")
    unsigned = dict(preflight)
    observed_preflight_identity = unsigned.pop("preflight_identity", None)
    if (
        set(preflight) != expected_keys
        or _json_bytes(preflight) != raw
        or preflight.get("schema") != SCHEMA
        or preflight.get("state") != "preflight-completed"
        or preflight.get("accepted_campaign_result") is not False
        or preflight.get("contract") != CONTRACT
        or preflight.get("contract_sha256") != campaign._canonical_hash(CONTRACT)
        or type(identity) is not dict
        or preflight.get("execution_identity_sha256")
        != campaign._canonical_hash(identity)
        or observed_preflight_identity != campaign._canonical_hash(unsigned)
    ):
        raise ValueError("scalar diagnostic preflight contract is invalid")
    _utc_timestamp(preflight["created_utc"], label="preflight created_utc")
    encoded_identity = json.dumps(identity, sort_keys=True).lower()
    if (
        identity.get("cpu_only") is not True
        or identity.get("accelerator") != "none"
        or "cuda" in encoded_identity
        or "gpu" in encoded_identity
    ):
        raise ValueError("scalar diagnostic preflight is not CPU-only")
    _validate_source_snapshot(identity.get("source_snapshot"))
    return preflight, raw


def _point_key(mode: str | None, q: float) -> str:
    if mode is None:
        return "center"
    prefix = mode.replace("_", "-")
    sign = "m" if q < 0.0 else "p"
    return f"{prefix}-{sign}{abs(q):.2f}"


def _point_plan(evidence: dict[str, Any]) -> tuple[dict[str, Any], ...]:
    transition_state = evidence["transition_state"]
    center = np.asarray(transition_state.coords, dtype=float)
    masses = np.asarray(evidence["masses_amu"], dtype=float)
    plan = [
        {
            "point_key": "center",
            "mode": None,
            "q_bohr_sqrt_amu": 0.0,
            "coords": center.copy(),
        }
    ]
    for mode_name in MODE_ORDER:
        mode = np.asarray(
            evidence["modes"][mode_name]["mass_weighted_vector"], dtype=float
        ).reshape(6, 3)
        for q in Q_STEPS:
            displacement = q * mode / np.sqrt(masses)[:, None]
            maximum = float(np.max(np.linalg.norm(displacement, axis=1)))
            if (
                not math.isfinite(maximum)
                or maximum
                >= CONTRACT["sampling"][
                    "maximum_per_atom_cartesian_displacement_bohr_exclusive"
                ]
            ):
                raise ValueError("mode displacement exceeds its Cartesian safety bound")
            plan.append(
                {
                    "point_key": _point_key(mode_name, q),
                    "mode": mode_name,
                    "q_bohr_sqrt_amu": q,
                    "maximum_atom_displacement_bohr": maximum,
                    "coords": center + displacement * BOHR_TO_ANGSTROM,
                }
            )
    keys = [point["point_key"] for point in plan]
    coordinate_hashes = [
        _sha(np.ascontiguousarray(point["coords"], dtype="<f8").tobytes())
        for point in plan
    ]
    if (
        tuple(keys) != POINT_ORDER
        or len(set(keys)) != 9
        or len(set(coordinate_hashes)) != 9
    ):
        raise RuntimeError(
            "scalar diagnostic point plan is not exactly nine unique geometries"
        )
    return tuple(plan)


def _density_snapshot(value: Any, *, label: str) -> tuple[np.ndarray, str, float]:
    density = np.asarray(value, dtype=float)
    if density.size < 1 or not np.all(np.isfinite(density)):
        raise ValueError(f"{label} must be nonempty and finite")
    density = np.frombuffer(
        np.ascontiguousarray(density, dtype="<f8").tobytes(), dtype="<f8"
    ).reshape(density.shape)
    return density, _sha(density.tobytes()), float(np.linalg.norm(density))


def _evaluate_scalar_energy(
    cluster: Any, settings: Any, *, dm0: np.ndarray | None
) -> tuple[dict[str, Any], np.ndarray]:
    """Evaluate one fresh CPU SCF energy; deliberately no gradient or Hessian."""

    fixed = replace(settings, use_gpu=False, grid_level=5)
    settings_identity = {
        "xc": fixed.xc,
        "basis": fixed.basis,
        "solvent": fixed.solvent,
        "dispersion": fixed.dispersion,
        "composite": fixed.composite,
        "grid_level": fixed.grid_level,
        "density_fit": fixed.density_fit,
        "use_gpu": fixed.use_gpu,
    }
    expected_settings_identity = {
        "xc": "pwb6k",
        "basis": "def2-svp",
        "solvent": None,
        "dispersion": "d3bj",
        "composite": None,
        "grid_level": 5,
        "density_fit": True,
        "use_gpu": False,
    }
    if settings_identity != expected_settings_identity:
        raise ValueError(
            "scalar diagnostic calculator settings are not exact CPU settings"
        )
    if cluster.charge != 0 or cluster.spin != 1:
        raise ValueError("scalar diagnostic requires the neutral doublet")
    mf = _make_scf(build_mol(cluster, fixed), fixed)
    mf.conv_tol = 1.0e-12
    mf.max_cycle = 150
    if getattr(mf.grids, "level", None) != 5:
        raise RuntimeError("PySCF grid level 5 was not established")
    initial_sha = None
    initial_density = None
    if dm0 is None:
        energy = float(mf.kernel())
    else:
        initial_density, initial_sha, _ = _density_snapshot(
            dm0, label="shared center density guess"
        )
        energy = float(mf.kernel(dm0=initial_density))
        if (
            _density_snapshot(initial_density, label="shared center density guess")[1]
            != initial_sha
        ):
            raise RuntimeError("SCF mutated the shared center density-matrix guess")
    spin = mf.spin_square()
    final_density, final_sha, final_norm = _density_snapshot(
        mf.make_rdm1(), label="final SCF density"
    )
    if initial_density is None:
        center_norm = final_norm
        density_delta = 0.0
        root_metric = "center-self"
        root_verified = True
    else:
        center_norm = float(np.linalg.norm(initial_density))
        denominator = max(final_norm, center_norm)
        if not math.isfinite(denominator) or denominator <= 0.0:
            raise ValueError("density continuity scale must be finite and positive")
        density_delta = float(
            np.linalg.norm(final_density - initial_density) / denominator
        )
        # No exact orbital/root-following invariant is available at this boundary.
        # S2 and density continuity remain explicit proxies, never proof of the root.
        root_metric = "unavailable-beyond-s2-and-density-proxies"
        root_verified = False
    return (
        {
            "backend": "pyscf-cpu",
            "electronic_hartree": energy,
            "scf_converged": bool(getattr(mf, "converged", False)),
            "grid_point_count": int(np.asarray(mf.grids.coords).shape[0]),
            "spin_2s": int(cluster.spin),
            "spin_square": [float(spin[0]), float(spin[1])],
            "scf_tolerance": float(mf.conv_tol),
            "scf_max_cycle": int(mf.max_cycle),
            "density_initial_guess_sha256": initial_sha,
            "density_final_sha256": final_sha,
            "density_shape": list(final_density.shape),
            "final_frobenius_norm": final_norm,
            "center_frobenius_norm": center_norm,
            "final_vs_center_normalized_frobenius_delta": density_delta,
            "robust_root_continuity_metric": root_metric,
            "robust_root_continuity_verified": root_verified,
            "gradient_computed": False,
            "hessian_computed": False,
        },
        final_density,
    )


_ATTESTED_EVALUATOR = _evaluate_scalar_energy
_ATTESTED_EVALUATOR_CODE = _evaluate_scalar_energy.__code__


def _validated_observation(observation: Any) -> dict[str, Any]:
    if type(observation) is not dict:
        raise ValueError("scalar evaluator result must be an object")
    expected = {
        "backend",
        "electronic_hartree",
        "scf_converged",
        "grid_point_count",
        "spin_2s",
        "spin_square",
        "scf_tolerance",
        "scf_max_cycle",
        "density_initial_guess_sha256",
        "density_final_sha256",
        "density_shape",
        "final_frobenius_norm",
        "center_frobenius_norm",
        "final_vs_center_normalized_frobenius_delta",
        "robust_root_continuity_metric",
        "robust_root_continuity_verified",
        "gradient_computed",
        "hessian_computed",
    }
    if set(observation) != expected:
        raise ValueError("scalar evaluator result schema is not exact")
    if (
        observation["backend"] != "pyscf-cpu"
        or type(observation["scf_converged"]) is not bool
        or type(observation["grid_point_count"]) is not int
        or type(observation["spin_2s"]) is not int
        or type(observation["scf_max_cycle"]) is not int
        or type(observation["robust_root_continuity_verified"]) is not bool
        or observation["gradient_computed"] is not False
        or observation["hessian_computed"] is not False
    ):
        raise ValueError("scalar evaluator identity is invalid")
    spin = observation["spin_square"]
    shape = observation["density_shape"]
    if type(spin) is not list or len(spin) != 2:
        raise ValueError("scalar evaluator spin result is invalid")
    if (
        type(shape) is not list
        or not shape
        or any(type(value) is not int or value <= 0 for value in shape)
    ):
        raise ValueError("scalar evaluator density shape is invalid")
    for field in ("density_final_sha256",):
        campaign._require_sha(
            campaign._require_json_string(observation[field], label=field),
            length=64,
            label=field,
        )
    initial_sha = observation["density_initial_guess_sha256"]
    if initial_sha is not None:
        campaign._require_sha(
            campaign._require_json_string(initial_sha, label="density initial SHA-256"),
            length=64,
            label="density initial SHA-256",
        )
    metric = observation["robust_root_continuity_metric"]
    if metric not in {"center-self", "unavailable-beyond-s2-and-density-proxies"}:
        raise ValueError("scalar evaluator root-continuity record is invalid")
    if (metric == "center-self") is not observation["robust_root_continuity_verified"]:
        raise ValueError("scalar evaluator root-continuity verdict is invalid")
    return {
        **observation,
        "electronic_hartree": _finite_float(
            observation["electronic_hartree"], label="scalar energy"
        ),
        "scf_tolerance": _finite_float(
            observation["scf_tolerance"], label="SCF tolerance"
        ),
        "spin_square": [
            _finite_float(spin[0], label="S2"),
            _finite_float(spin[1], label="spin multiplicity"),
        ],
        "final_frobenius_norm": _finite_float(
            observation["final_frobenius_norm"], label="final density norm"
        ),
        "center_frobenius_norm": _finite_float(
            observation["center_frobenius_norm"], label="center density norm"
        ),
        "final_vs_center_normalized_frobenius_delta": _finite_float(
            observation["final_vs_center_normalized_frobenius_delta"],
            label="density continuity delta",
        ),
    }


def _gate(name: str, value: Any, comparison: str, limit: Any) -> dict[str, Any]:
    if comparison == "equals":
        passed = type(value) is type(limit) and value == limit
        observed = value
    elif comparison == "inclusive_range":
        observed = _finite_float(value, label=name)
        low, high = limit
        passed = low <= observed <= high
    elif comparison == "strictly_less_than":
        observed = _finite_float(value, label=name)
        passed = observed < float(limit)
    elif comparison == "strictly_greater_than":
        observed = _finite_float(value, label=name)
        passed = observed > float(limit)
    else:
        raise ValueError(f"unsupported gate comparison: {comparison}")
    return {
        "name": name,
        "value": observed,
        "comparison": comparison,
        "limit": limit,
        "passed": bool(passed),
    }


def _center_source_state_is_accepted(observation: dict[str, Any]) -> bool:
    thresholds = CONTRACT["thresholds"]
    return (
        observation["scf_converged"] is True
        and observation["spin_2s"] == 1
        and thresholds["spin_square_inclusive_minimum"]
        <= observation["spin_square"][0]
        <= thresholds["spin_square_inclusive_maximum"]
        and abs(observation["electronic_hartree"] - SOURCE_CENTER_ENERGY_HARTREE)
        < thresholds["center_energy_vs_source_exclusive_maximum_delta_hartree"]
        and abs(observation["spin_square"][0] - SOURCE_CENTER_S2)
        < thresholds["center_s2_vs_source_exclusive_maximum_delta"]
        and observation["robust_root_continuity_verified"] is True
    )


def _analyze(points: list[dict[str, Any]]) -> dict[str, Any]:
    if len(points) != 9 or [point["point_key"] for point in points] != list(
        POINT_ORDER
    ):
        raise ValueError("analysis requires the exact ordered nine-point set")
    by_key = {point["point_key"]: point for point in points}
    energies = {
        key: _finite_float(point["observation"]["electronic_hartree"], label=key)
        for key, point in by_key.items()
    }
    center = energies["center"]
    gates = []
    for point in points:
        observation = point["observation"]
        key = point["point_key"]
        point_gates = [
            _gate(f"{key}:scf_converged", observation["scf_converged"], "equals", True),
            _gate(f"{key}:spin_2s", observation["spin_2s"], "equals", 1),
            _gate(
                f"{key}:grid_point_count",
                observation["grid_point_count"],
                "equals",
                CONTRACT["method"]["required_grid_point_count"],
            ),
            _gate(
                f"{key}:scf_tolerance",
                observation["scf_tolerance"],
                "equals",
                CONTRACT["method"]["scf_tolerance"],
            ),
            _gate(
                f"{key}:scf_max_cycle",
                observation["scf_max_cycle"],
                "equals",
                CONTRACT["method"]["scf_max_cycle"],
            ),
            _gate(
                f"{key}:S2",
                observation["spin_square"][0],
                "inclusive_range",
                [
                    CONTRACT["thresholds"]["spin_square_inclusive_minimum"],
                    CONTRACT["thresholds"]["spin_square_inclusive_maximum"],
                ],
            ),
        ]
        if key == "center":
            point_gates.extend(
                (
                    _gate(
                        "center:source_energy_delta",
                        abs(
                            observation["electronic_hartree"]
                            - SOURCE_CENTER_ENERGY_HARTREE
                        ),
                        "strictly_less_than",
                        CONTRACT["thresholds"][
                            "center_energy_vs_source_exclusive_maximum_delta_hartree"
                        ],
                    ),
                    _gate(
                        "center:source_S2_delta",
                        abs(observation["spin_square"][0] - SOURCE_CENTER_S2),
                        "strictly_less_than",
                        CONTRACT["thresholds"][
                            "center_s2_vs_source_exclusive_maximum_delta"
                        ],
                    ),
                    _gate(
                        "center:root_continuity",
                        observation["robust_root_continuity_verified"],
                        "equals",
                        True,
                    ),
                )
            )
        else:
            center_observation = by_key["center"]["observation"]
            point_gates.extend(
                (
                    _gate(
                        f"{key}:S2_delta_from_center",
                        abs(
                            observation["spin_square"][0]
                            - center_observation["spin_square"][0]
                        ),
                        "strictly_less_than",
                        CONTRACT["thresholds"][
                            "displaced_vs_center_spin_square_exclusive_maximum_delta"
                        ],
                    ),
                    _gate(
                        f"{key}:density_delta_from_center",
                        observation["final_vs_center_normalized_frobenius_delta"],
                        "strictly_less_than",
                        CONTRACT["thresholds"][
                            "displaced_final_density_vs_center_exclusive_maximum_"
                            "normalized_frobenius_delta"
                        ],
                    ),
                    _gate(
                        f"{key}:robust_root_continuity",
                        observation["robust_root_continuity_verified"],
                        "equals",
                        True,
                    ),
                )
            )
        gates.extend(point_gates)
    mode_results = {}
    for mode_name in MODE_ORDER:
        prefix = mode_name.replace("_", "-")
        e_m10 = energies[f"{prefix}-m0.10"]
        e_m05 = energies[f"{prefix}-m0.05"]
        e_p05 = energies[f"{prefix}-p0.05"]
        e_p10 = energies[f"{prefix}-p0.10"]
        k05 = (e_p05 - 2.0 * center + e_m05) / (0.05**2)
        k10 = (e_p10 - 2.0 * center + e_m10) / (0.10**2)
        g05 = (e_p05 - e_m05) / (2.0 * 0.05)
        g10 = (e_p10 - e_m10) / (2.0 * 0.10)
        k_r = (4.0 * k05 - k10) / 3.0
        g_r = (4.0 * g05 - g10) / 3.0
        for label, value in (
            ("kappa_0.05", k05),
            ("kappa_0.10", k10),
            ("kappa_R", k_r),
            ("g_0.05", g05),
            ("g_0.10", g10),
            ("g_R", g_r),
        ):
            _finite_float(value, label=f"{mode_name} {label}")
        q_star = (
            None
            if k_r == 0.0
            else _finite_float(-g_r / k_r, label=f"{mode_name} estimated stationary q")
        )
        convergence_limit = max(
            CONTRACT["thresholds"]["curvature_convergence_absolute_floor_exclusive"],
            CONTRACT["thresholds"]["curvature_convergence_relative_fraction_exclusive"]
            * abs(k_r),
        )
        frequency = float(hessian_eigenvalues_to_wavenumbers_cm(np.array([k_r]))[0])
        q_star_gate = (
            {
                "name": f"{mode_name}:estimated_stationary_q",
                "value": None,
                "comparison": "strictly_less_than",
                "limit": CONTRACT["thresholds"][
                    "estimated_stationary_absolute_q_exclusive_maximum"
                ],
                "passed": False,
            }
            if q_star is None
            else _gate(
                f"{mode_name}:estimated_stationary_q",
                abs(q_star),
                "strictly_less_than",
                CONTRACT["thresholds"][
                    "estimated_stationary_absolute_q_exclusive_maximum"
                ],
            )
        )
        mode_gates = [
            _gate(
                f"{mode_name}:curvature_convergence",
                abs(k_r - k05),
                "strictly_less_than",
                convergence_limit,
            ),
            q_star_gate,
        ]
        if mode_name == "unstable":
            mode_gates.extend(
                (
                    _gate(
                        "unstable:kappa_R",
                        k_r,
                        "strictly_less_than",
                        CONTRACT["thresholds"]["unstable_curvature_exclusive_maximum"],
                    ),
                    _gate(
                        "unstable:absolute_frequency_cm",
                        abs(frequency),
                        "strictly_greater_than",
                        CONTRACT["thresholds"][
                            "unstable_absolute_frequency_cm_exclusive_minimum"
                        ],
                    ),
                )
            )
        else:
            mode_gates.append(
                _gate(
                    "lowest_positive:kappa_R",
                    k_r,
                    "strictly_greater_than",
                    CONTRACT["thresholds"][
                        "lowest_positive_curvature_exclusive_minimum"
                    ],
                )
            )
        gates.extend(mode_gates)
        mode_results[mode_name] = {
            "energy_deltas_hartree": {
                "m0.10": e_m10 - center,
                "m0.05": e_m05 - center,
                "p0.05": e_p05 - center,
                "p0.10": e_p10 - center,
            },
            "kappa_0.05_hartree_per_bohr2_amu": k05,
            "kappa_0.10_hartree_per_bohr2_amu": k10,
            "kappa_richardson_hartree_per_bohr2_amu": k_r,
            "curvature_convergence_absolute_delta": abs(k_r - k05),
            "curvature_convergence_exclusive_limit": convergence_limit,
            "g_0.05_hartree_per_bohr_sqrt_amu": g05,
            "g_0.10_hartree_per_bohr_sqrt_amu": g10,
            "g_richardson_hartree_per_bohr_sqrt_amu": g_r,
            "estimated_stationary_q_bohr_sqrt_amu": q_star,
            "signed_wavenumber_cm": frequency,
            "gates": mode_gates,
        }
    passed = all(gate["passed"] is True for gate in gates)
    return {
        "diagnostic_passed": passed,
        "accepted_campaign_result": False,
        "full_stationarity_resolved": False,
        "stationarity_statement": (
            "No center gradient was computed; only two projected first derivatives "
            "were estimated from scalar energies, so full stationarity is unresolved."
        ),
        "point_count": len(points),
        "unique_geometry_count": len({point["geometry_sha256"] for point in points}),
        "center_energy_hartree": center,
        "center_density_sha256": by_key["center"]["observation"][
            "density_final_sha256"
        ],
        "robust_root_continuity_available": False,
        "same_electronic_state_verified": False,
        "root_continuity_statement": (
            "All displaced SCFs used the exact center density and passed explicit "
            "S2/density proxy gates only when recorded by those gates. No exact "
            "robust root-following metric is available, so every displaced "
            "robust-root gate fails closed."
        ),
        "modes": mode_results,
        "gates": gates,
    }


def _initial_status(preflight: dict[str, Any], preflight_raw: bytes) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "state": "running",
        "accepted_campaign_result": False,
        "diagnostic_passed": False,
        "route": ROUTE,
        "started_utc": _now_utc(),
        "preflight_sha256": _sha(preflight_raw),
        "preflight_identity": preflight["preflight_identity"],
        "execution_identity_sha256": preflight["execution_identity_sha256"],
        "contract_sha256": preflight["contract_sha256"],
        "source_receipt_sha256": SOURCE_RECEIPT_SHA256,
        "source_evidence_sha256": campaign._canonical_hash(
            preflight["source_evidence"]
        ),
        "source_snapshot": preflight["execution_identity"]["source_snapshot"],
        "center_density": None,
        "completed_points": [],
        "current_point": None,
    }


def _point_receipt(
    definition: dict[str, Any],
    observation: dict[str, Any],
    status: dict[str, Any],
    geometry_sha: str,
) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "kind": "scalar-energy-single-point",
        "accepted_campaign_result": False,
        "route": ROUTE,
        "point_key": definition["point_key"],
        "mode": definition["mode"],
        "q_bohr_sqrt_amu": definition["q_bohr_sqrt_amu"],
        "maximum_atom_displacement_bohr": definition.get(
            "maximum_atom_displacement_bohr", 0.0
        ),
        "geometry_sha256": geometry_sha,
        "preflight_sha256": status["preflight_sha256"],
        "preflight_identity": status["preflight_identity"],
        "execution_identity_sha256": status["execution_identity_sha256"],
        "contract_sha256": status["contract_sha256"],
        "source_receipt_sha256": SOURCE_RECEIPT_SHA256,
        "source_evidence_sha256": status["source_evidence_sha256"],
        "source_snapshot": status["source_snapshot"],
        "center_density": status["center_density"],
        "observation": observation,
    }


def _write_failed_terminal(
    claim: hessian_diagnostic._OutputRootClaim,
    status: dict[str, Any],
    exc: BaseException,
    *,
    dead_man: bool = False,
) -> dict[str, Any]:
    """Publish the failed terminal commit marker before converging status."""

    detail = f"{type(exc).__name__}: {exc}"
    finished = _now_utc()
    failed_status = {
        **status,
        "state": "failed",
        "accepted_campaign_result": False,
        "diagnostic_passed": False,
        "error": detail,
        "finished_utc": finished,
    }
    failed_status_raw = _json_bytes(failed_status)
    terminal = {
        "schema": SCHEMA,
        "state": "failed",
        "accepted_campaign_result": False,
        "diagnostic_passed": False,
        "route": ROUTE,
        "detail": detail,
        "dead_man_finalized": dead_man,
        "completed_point_count": len(status.get("completed_points", [])),
        "current_point": status.get("current_point"),
        "preflight_sha256": status.get("preflight_sha256"),
        "preflight_identity": status.get("preflight_identity"),
        "execution_identity_sha256": status.get("execution_identity_sha256"),
        "contract_sha256": status.get("contract_sha256"),
        "source_receipt_sha256": status.get("source_receipt_sha256"),
        "source_evidence_sha256": status.get("source_evidence_sha256"),
        "source_snapshot": status.get("source_snapshot"),
        "status": failed_status,
        "status_sha256": _sha(failed_status_raw),
        "started_utc": status.get("started_utc"),
        "finished_utc": finished,
    }
    if campaign._canonical_hash(_current_identity()) != status.get(
        "execution_identity_sha256"
    ):
        raise RuntimeError("execution identity changed before failed terminal commit")
    hessian_diagnostic._terminal_write_noreplace(
        claim, claim.path / "terminal.json", terminal
    )
    with suppress(BaseException):
        hessian_diagnostic._claimed_write(
            claim, claim.path / "status.json", failed_status_raw
        )
        # The terminal is the commit marker. finalize_if_running() repairs status.
    return terminal


def _run_locked(claim: hessian_diagnostic._OutputRootClaim) -> dict[str, Any]:
    if (
        _assert_resident_bindings is not _ATTESTED_RESIDENT_ASSERT
        or _assert_resident_bindings.__code__ is not _ATTESTED_RESIDENT_ASSERT_CODE
    ):
        raise RuntimeError("resident-boundary validator was replaced")
    _ATTESTED_RESIDENT_ASSERT()
    if claim.root_created:
        raise FileNotFoundError(
            "scalar diagnostic requires an immutable preflight root"
        )
    if (
        _evaluate_scalar_energy is not _ATTESTED_EVALUATOR
        or _evaluate_scalar_energy.__code__ is not _ATTESTED_EVALUATOR_CODE
    ):
        raise RuntimeError("attested scalar evaluator was replaced")
    preflight, preflight_raw = _validated_preflight(claim)
    assert claim.root_descriptor is not None
    if set(os.listdir(claim.root_descriptor)) != {"preflight.json"}:
        raise FileExistsError(
            "scalar diagnostic never resumes a used or partial root; "
            "create a fresh root"
        )
    identity = _current_identity()
    campaign._strict_json_equal(
        identity, preflight["execution_identity"], label="scalar execution identity"
    )
    evidence = _load_source_evidence(SOURCE_ROOT)
    campaign._strict_json_equal(
        _public_evidence(evidence),
        preflight["source_evidence"],
        label="scalar source evidence",
    )
    if _current_identity() != identity:
        raise RuntimeError("execution identity changed before scalar evaluation")
    plan = _point_plan(evidence)
    status = _initial_status(preflight, preflight_raw)
    hessian_diagnostic._claimed_write_noreplace(
        claim, claim.path / "status.json", _json_bytes(status)
    )
    points: list[dict[str, Any]] = []
    center_density: np.ndarray | None = None
    center_density_sha: str | None = None
    receipt_published = False
    try:
        template = evidence["transition_state"]
        for index, definition in enumerate(plan):
            key = definition["point_key"]
            if (index == 0) is not (key == "center"):
                raise RuntimeError(
                    "center must be evaluated before every displaced point"
                )
            status["current_point"] = key
            hessian_diagnostic._claimed_write(
                claim, claim.path / "status.json", _json_bytes(status)
            )
            cluster = replace(template, coords=np.asarray(definition["coords"]))
            geometry_sha = _sha(
                np.ascontiguousarray(cluster.coords, dtype="<f8").tobytes()
            )
            raw_observation, final_density = _evaluate_scalar_energy(
                cluster,
                evidence["settings"],
                dm0=None if key == "center" else center_density,
            )
            observation = _validated_observation(raw_observation)
            if key == "center" and not _center_source_state_is_accepted(observation):
                raise ScientificRejection(
                    "center SCF did not reproduce the pinned source electronic state"
                )
            verified_density, verified_sha, _ = _density_snapshot(
                final_density, label=f"{key} returned final density"
            )
            if (
                verified_sha != observation["density_final_sha256"]
                or list(verified_density.shape) != observation["density_shape"]
            ):
                raise RuntimeError("SCF final density disagrees with its observation")
            if key == "center":
                if observation["density_initial_guess_sha256"] is not None:
                    raise RuntimeError(
                        "center SCF unexpectedly used a supplied density"
                    )
                center_density = verified_density
                center_density.setflags(write=False)
                center_density_sha = verified_sha
                center_metadata = {
                    "path": "center-density.f64",
                    "sha256": center_density_sha,
                    "shape": list(center_density.shape),
                    "dtype": "little-endian float64",
                    "units": "electrons",
                }
                hessian_diagnostic._claimed_write_noreplace(
                    claim, claim.path / "center-density.f64", center_density.tobytes()
                )
                status["center_density"] = center_metadata
            elif (
                center_density is None
                or center_density_sha is None
                or observation["density_initial_guess_sha256"] != center_density_sha
                or observation["center_frobenius_norm"]
                != float(np.linalg.norm(center_density))
            ):
                raise RuntimeError(
                    "displaced SCF did not bind the exact converged center density"
                )
            point = _point_receipt(definition, observation, status, geometry_sha)
            point_raw = _json_bytes(point)
            point_sha = _sha(point_raw)
            hessian_diagnostic._claimed_write_noreplace(
                claim, claim.path / f"point-{key}.json", point_raw
            )
            points.append(point)
            status["completed_points"].append(
                {"point_key": key, "receipt_sha256": point_sha}
            )
            status["current_point"] = None
            hessian_diagnostic._claimed_write(
                claim, claim.path / "status.json", _json_bytes(status)
            )
        if center_density is None or center_density_sha is None or len(points) != 9:
            raise RuntimeError(
                "scalar evaluation did not produce the exact nine-point set"
            )
        analysis = _analyze(points)
        finished = _now_utc()
        receipt = {
            "schema": SCHEMA,
            "state": "completed",
            "accepted_campaign_result": False,
            "diagnostic_passed": analysis["diagnostic_passed"],
            "purpose": "mode-projected scalar-energy evidence only",
            "authority": (
                "no TS qualification, IRC, SCT, D2b gate flip, or D3b authorization"
            ),
            "route": ROUTE,
            "contract": CONTRACT,
            "contract_sha256": preflight["contract_sha256"],
            "preflight_sha256": status["preflight_sha256"],
            "preflight_identity": status["preflight_identity"],
            "execution_identity": identity,
            "execution_identity_sha256": status["execution_identity_sha256"],
            "source_snapshot": status["source_snapshot"],
            "source_evidence": preflight["source_evidence"],
            "source_evidence_sha256": status["source_evidence_sha256"],
            "center_density": status["center_density"],
            "point_receipts": status["completed_points"],
            "analysis": analysis,
            "started_utc": status["started_utc"],
            "finished_utc": finished,
        }
        receipt_raw = _json_bytes(receipt)
        receipt_sha = _sha(receipt_raw)
        if _current_identity() != identity:
            raise RuntimeError("execution identity changed before receipt commit")
        hessian_diagnostic._claimed_write_noreplace(
            claim, claim.path / "receipt.json", receipt_raw
        )
        receipt_published = True
        completed_status = {
            **status,
            "state": "completed",
            "diagnostic_passed": analysis["diagnostic_passed"],
            "receipt_sha256": receipt_sha,
            "finished_utc": finished,
        }
        completed_status_raw = _json_bytes(completed_status)
        terminal = {
            "schema": SCHEMA,
            "state": "completed",
            "accepted_campaign_result": False,
            "diagnostic_passed": analysis["diagnostic_passed"],
            "route": ROUTE,
            "detail": None,
            "dead_man_finalized": False,
            "completed_point_count": 9,
            "current_point": None,
            "preflight_sha256": status["preflight_sha256"],
            "preflight_identity": status["preflight_identity"],
            "execution_identity_sha256": status["execution_identity_sha256"],
            "contract_sha256": status["contract_sha256"],
            "source_receipt_sha256": SOURCE_RECEIPT_SHA256,
            "source_evidence_sha256": status["source_evidence_sha256"],
            "source_snapshot": status["source_snapshot"],
            "receipt_sha256": receipt_sha,
            "status": completed_status,
            "status_sha256": _sha(completed_status_raw),
            "started_utc": status["started_utc"],
            "finished_utc": finished,
        }
        if _current_identity() != identity:
            raise RuntimeError("execution identity changed before terminal commit")
        hessian_diagnostic._terminal_write_noreplace(
            claim, claim.path / "terminal.json", terminal
        )
        hessian_diagnostic._claimed_write(
            claim, claim.path / "status.json", completed_status_raw
        )
        return receipt
    except BaseException as exc:
        if not receipt_published and not hessian_diagnostic._claimed_root_entry_exists(
            claim, claim.path / "terminal.json"
        ):
            _write_failed_terminal(claim, status, exc)
        raise


def run(output_root: Path) -> dict[str, Any]:
    output = campaign._safe_absolute_root(output_root)
    with hessian_diagnostic._exclusive_output_claim(output) as claim:
        return _run_locked(claim)


def _validate_source_snapshot(snapshot: Any) -> dict[str, Any]:
    if type(snapshot) is not dict or set(snapshot) != {
        "git_commit",
        "git_tree",
        "source_files",
    }:
        raise ValueError("scalar source snapshot is invalid")
    for field in ("git_commit", "git_tree"):
        campaign._require_sha(
            campaign._require_json_string(snapshot[field], label=field),
            length=40,
            label=field,
        )
    source_files = snapshot["source_files"]
    if type(source_files) is not dict or not source_files:
        raise ValueError("scalar source snapshot file inventory is invalid")
    for relative, digest in source_files.items():
        if (
            type(relative) is not str
            or not relative
            or Path(relative).is_absolute()
            or ".." in Path(relative).parts
        ):
            raise ValueError("scalar source snapshot path is invalid")
        campaign._require_sha(
            campaign._require_json_string(digest, label="source file SHA-256"),
            length=64,
            label="source file SHA-256",
        )
    return snapshot


def _status_common_keys() -> set[str]:
    return {
        "schema",
        "state",
        "accepted_campaign_result",
        "diagnostic_passed",
        "route",
        "started_utc",
        "preflight_sha256",
        "preflight_identity",
        "execution_identity_sha256",
        "contract_sha256",
        "source_receipt_sha256",
        "source_evidence_sha256",
        "source_snapshot",
        "center_density",
        "completed_points",
        "current_point",
    }


def _validate_center_density_metadata(value: Any) -> dict[str, Any] | None:
    if value is None:
        return None
    if type(value) is not dict or set(value) != {
        "path",
        "sha256",
        "shape",
        "dtype",
        "units",
    }:
        raise ValueError("center density artifact metadata is invalid")
    if (
        value["path"] != "center-density.f64"
        or value["dtype"] != "little-endian float64"
        or value["units"] != "electrons"
        or type(value["shape"]) is not list
        or not value["shape"]
        or any(type(item) is not int or item <= 0 for item in value["shape"])
    ):
        raise ValueError("center density artifact contract is invalid")
    campaign._require_sha(
        campaign._require_json_string(value["sha256"], label="center density SHA-256"),
        length=64,
        label="center density SHA-256",
    )
    return value


def _validate_status(
    status: Any,
    raw: bytes,
    preflight: dict[str, Any],
    preflight_raw: bytes,
) -> dict[str, Any]:
    if type(status) is not dict or _json_bytes(status) != raw:
        raise ValueError("scalar diagnostic status is not strict canonical JSON")
    state = status.get("state")
    expected = _status_common_keys()
    if state == "completed":
        expected |= {"receipt_sha256", "finished_utc"}
    elif state == "failed":
        expected |= {"error", "finished_utc"}
    elif state != "running":
        raise ValueError("scalar diagnostic status state is invalid")
    completed = status.get("completed_points")
    if (
        set(status) != expected
        or status.get("schema") != SCHEMA
        or status.get("accepted_campaign_result") is not False
        or type(status.get("diagnostic_passed")) is not bool
        or status.get("route") != ROUTE
        or type(completed) is not list
        or len(completed) > 9
        or status.get("current_point") not in {None, *POINT_ORDER}
    ):
        raise ValueError("scalar diagnostic status contract is invalid")
    if state in {"running", "failed"} and status["diagnostic_passed"] is not False:
        raise ValueError("non-completed scalar status cannot pass")
    expected_bindings = {
        "preflight_sha256": _sha(preflight_raw),
        "preflight_identity": preflight["preflight_identity"],
        "execution_identity_sha256": preflight["execution_identity_sha256"],
        "contract_sha256": campaign._canonical_hash(CONTRACT),
        "source_receipt_sha256": SOURCE_RECEIPT_SHA256,
        "source_evidence_sha256": campaign._canonical_hash(
            preflight["source_evidence"]
        ),
        "source_snapshot": preflight["execution_identity"]["source_snapshot"],
    }
    for field, expected_value in expected_bindings.items():
        campaign._strict_json_equal(
            status.get(field), expected_value, label=f"status {field}"
        )
    _validate_source_snapshot(status["source_snapshot"])
    started = _utc_timestamp(status.get("started_utc"), label="status started_utc")
    created = _utc_timestamp(
        preflight.get("created_utc"), label="preflight created_utc"
    )
    if started < created:
        raise ValueError("scalar status predates its preflight")
    if state != "running":
        finished = _utc_timestamp(
            status.get("finished_utc"), label="status finished_utc"
        )
        if finished < started:
            raise ValueError("scalar status finished before it started")
    if state == "failed" and type(status.get("error")) is not str:
        raise ValueError("failed scalar status error is invalid")
    if state == "completed":
        campaign._require_sha(
            campaign._require_json_string(
                status.get("receipt_sha256"), label="status receipt SHA-256"
            ),
            length=64,
            label="status receipt SHA-256",
        )
        if len(completed) != 9 or status["current_point"] is not None:
            raise ValueError("completed scalar status ancestry is incomplete")
    for index, binding in enumerate(completed):
        if type(binding) is not dict or set(binding) != {
            "point_key",
            "receipt_sha256",
        }:
            raise ValueError("scalar status point ancestry is malformed")
        if binding.get("point_key") != POINT_ORDER[index]:
            raise ValueError("scalar status point ancestry is out-of-order")
        campaign._require_sha(
            campaign._require_json_string(
                binding.get("receipt_sha256"), label="point receipt SHA-256"
            ),
            length=64,
            label="point receipt SHA-256",
        )
    if status["current_point"] is not None:
        next_index = len(completed)
        if next_index >= 9 or status["current_point"] != POINT_ORDER[next_index]:
            raise ValueError("scalar status current point is inconsistent")
    center_density = _validate_center_density_metadata(status["center_density"])
    if center_density is None and completed:
        raise ValueError("scalar status center density ancestry is inconsistent")
    if (
        center_density is not None
        and not completed
        and status["current_point"] != "center"
    ):
        raise ValueError("scalar status center density ancestry is inconsistent")
    return status


def _read_validated_points(
    claim: hessian_diagnostic._OutputRootClaim,
    status: dict[str, Any],
    evidence: dict[str, Any],
) -> list[dict[str, Any]]:
    bindings = status["completed_points"]
    plan = _point_plan(evidence)
    center_metadata = status["center_density"]
    if center_metadata is not None:
        center_raw = campaign._read_bounded_regular_snapshot(
            claim.path / center_metadata["path"],
            label="persisted scalar center density",
            maximum_bytes=8 * math.prod(center_metadata["shape"]),
        )
        if (
            len(center_raw) != 8 * math.prod(center_metadata["shape"])
            or _sha(center_raw) != center_metadata["sha256"]
        ):
            raise ValueError("persisted scalar center density fingerprint drifted")
    points: list[dict[str, Any]] = []
    for index, binding in enumerate(bindings):
        definition = plan[index]
        key = definition["point_key"]
        point, point_raw = hessian_diagnostic._claimed_json_object(
            claim,
            claim.path / f"point-{key}.json",
            label=f"scalar point {key}",
            max_bytes=256 * 1024,
        )
        if (
            _json_bytes(point) != point_raw
            or _sha(point_raw) != binding["receipt_sha256"]
        ):
            raise ValueError("scalar point receipt is noncanonical or hash-mismatched")
        observation = _validated_observation(point.get("observation"))
        geometry_sha = _sha(
            np.ascontiguousarray(definition["coords"], dtype="<f8").tobytes()
        )
        expected_point = _point_receipt(definition, observation, status, geometry_sha)
        if point != expected_point:
            raise ValueError("scalar point receipt content is invalid")
        assert center_metadata is not None
        if index == 0:
            if (
                observation["density_initial_guess_sha256"] is not None
                or observation["density_final_sha256"] != center_metadata["sha256"]
                or observation["density_shape"] != center_metadata["shape"]
                or observation["robust_root_continuity_verified"] is not True
            ):
                raise ValueError("scalar center density/root binding is invalid")
        elif (
            observation["density_initial_guess_sha256"] != center_metadata["sha256"]
            or observation["density_shape"] != center_metadata["shape"]
            or observation["robust_root_continuity_metric"]
            != "unavailable-beyond-s2-and-density-proxies"
            or observation["robust_root_continuity_verified"] is not False
        ):
            raise ValueError("scalar displaced density/root binding is invalid")
        points.append(point)
    return points


def _validate_receipt(
    claim: hessian_diagnostic._OutputRootClaim,
    status: dict[str, Any],
    preflight: dict[str, Any],
    preflight_raw: bytes,
    evidence: dict[str, Any],
) -> tuple[dict[str, Any], bytes, list[dict[str, Any]]]:
    receipt, receipt_raw = hessian_diagnostic._claimed_json_object(
        claim,
        claim.path / "receipt.json",
        label="scalar diagnostic receipt",
        max_bytes=4 * 1024 * 1024,
    )
    expected_keys = {
        "schema",
        "state",
        "accepted_campaign_result",
        "diagnostic_passed",
        "purpose",
        "authority",
        "route",
        "contract",
        "contract_sha256",
        "preflight_sha256",
        "preflight_identity",
        "execution_identity",
        "execution_identity_sha256",
        "source_snapshot",
        "source_evidence",
        "source_evidence_sha256",
        "center_density",
        "point_receipts",
        "analysis",
        "started_utc",
        "finished_utc",
    }
    if (
        type(receipt) is not dict
        or set(receipt) != expected_keys
        or _json_bytes(receipt) != receipt_raw
        or receipt.get("schema") != SCHEMA
        or receipt.get("state") != "completed"
        or receipt.get("accepted_campaign_result") is not False
        or type(receipt.get("diagnostic_passed")) is not bool
        or receipt.get("route") != ROUTE
        or receipt.get("contract") != CONTRACT
        or receipt.get("contract_sha256") != campaign._canonical_hash(CONTRACT)
        or receipt.get("preflight_sha256") != _sha(preflight_raw)
        or receipt.get("preflight_identity") != preflight["preflight_identity"]
        or receipt.get("execution_identity") != preflight["execution_identity"]
        or receipt.get("execution_identity_sha256")
        != preflight["execution_identity_sha256"]
        or receipt.get("source_snapshot")
        != preflight["execution_identity"]["source_snapshot"]
        or receipt.get("source_evidence") != preflight["source_evidence"]
        or receipt.get("source_evidence_sha256")
        != campaign._canonical_hash(preflight["source_evidence"])
        or receipt.get("center_density") != status["center_density"]
        or receipt.get("point_receipts") != status["completed_points"]
        or receipt.get("started_utc") != status["started_utc"]
    ):
        raise ValueError("completed scalar receipt contract or ancestry is invalid")
    if len(status["completed_points"]) != 9:
        raise ValueError("completed scalar receipt lacks all nine point ancestors")
    started = _utc_timestamp(receipt["started_utc"], label="receipt started_utc")
    finished = _utc_timestamp(receipt["finished_utc"], label="receipt finished_utc")
    if finished < started:
        raise ValueError("completed scalar receipt finished before it started")
    points = _read_validated_points(claim, status, evidence)
    analysis = _analyze(points)
    if (
        receipt.get("analysis") != analysis
        or receipt["diagnostic_passed"] is not analysis["diagnostic_passed"]
    ):
        raise ValueError("completed scalar receipt analysis is not reproducible")
    return receipt, receipt_raw, points


def _validate_terminal(payload: dict[str, Any]) -> dict[str, Any]:
    common_keys = {
        "schema",
        "state",
        "accepted_campaign_result",
        "diagnostic_passed",
        "route",
        "detail",
        "dead_man_finalized",
        "completed_point_count",
        "current_point",
        "preflight_sha256",
        "preflight_identity",
        "execution_identity_sha256",
        "contract_sha256",
        "source_receipt_sha256",
        "source_evidence_sha256",
        "source_snapshot",
        "status",
        "status_sha256",
        "started_utc",
        "finished_utc",
    }
    state = payload.get("state")
    expected_keys = common_keys | (
        {"receipt_sha256"} if state == "completed" else set()
    )
    if (
        type(payload) is not dict
        or set(payload) != expected_keys
        or payload.get("schema") != SCHEMA
        or state not in {"completed", "failed"}
        or payload.get("accepted_campaign_result") is not False
        or type(payload.get("diagnostic_passed")) is not bool
        or payload.get("route") != ROUTE
        or type(payload.get("dead_man_finalized")) is not bool
        or type(payload.get("completed_point_count")) is not int
        or not 0 <= payload["completed_point_count"] <= 9
    ):
        raise ValueError("scalar diagnostic terminal is invalid")
    for field in (
        "preflight_sha256",
        "preflight_identity",
        "execution_identity_sha256",
        "contract_sha256",
        "source_receipt_sha256",
        "source_evidence_sha256",
        "status_sha256",
    ):
        campaign._require_sha(
            campaign._require_json_string(
                payload.get(field), label=f"terminal {field}"
            ),
            length=64,
            label=f"terminal {field}",
        )
    _validate_source_snapshot(payload.get("source_snapshot"))
    started = _utc_timestamp(payload.get("started_utc"), label="terminal started_utc")
    finished = _utc_timestamp(
        payload.get("finished_utc"), label="terminal finished_utc"
    )
    if finished < started:
        raise ValueError("scalar terminal finished before it started")
    embedded_status = payload.get("status")
    if (
        type(embedded_status) is not dict
        or _sha(_json_bytes(embedded_status)) != payload["status_sha256"]
        or embedded_status.get("state") != state
        or embedded_status.get("diagnostic_passed") is not payload["diagnostic_passed"]
        or embedded_status.get("started_utc") != payload.get("started_utc")
        or embedded_status.get("finished_utc") != payload.get("finished_utc")
    ):
        raise ValueError("scalar terminal embedded status binding is invalid")
    if (
        payload["contract_sha256"] != campaign._canonical_hash(CONTRACT)
        or payload["source_receipt_sha256"] != SOURCE_RECEIPT_SHA256
    ):
        raise ValueError("scalar terminal contract/source hash drifted")
    if state == "completed":
        campaign._require_sha(
            campaign._require_json_string(
                payload.get("receipt_sha256"), label="terminal receipt SHA-256"
            ),
            length=64,
            label="terminal receipt SHA-256",
        )
        if (
            payload["completed_point_count"] != 9
            or payload["current_point"] is not None
            or payload["detail"] is not None
            or payload["dead_man_finalized"] is not False
            or embedded_status.get("receipt_sha256") != payload.get("receipt_sha256")
        ):
            raise ValueError("completed scalar diagnostic terminal is inconsistent")
    elif (
        payload["diagnostic_passed"] is not False
        or type(payload["detail"]) is not str
        or embedded_status.get("error") != payload["detail"]
    ):
        raise ValueError("failed scalar diagnostic terminal is inconsistent")
    return payload


def _running_projection(status: dict[str, Any]) -> dict[str, Any]:
    return {key: status[key] for key in _status_common_keys()}


def _completed_status(
    status: dict[str, Any], receipt: dict[str, Any], receipt_sha: str
) -> dict[str, Any]:
    return {
        **_running_projection(status),
        "state": "completed",
        "diagnostic_passed": receipt["diagnostic_passed"],
        "receipt_sha256": receipt_sha,
        "finished_utc": receipt["finished_utc"],
    }


def _completed_terminal(
    status: dict[str, Any], receipt: dict[str, Any], receipt_sha: str
) -> tuple[dict[str, Any], bytes]:
    desired_status = _completed_status(status, receipt, receipt_sha)
    desired_raw = _json_bytes(desired_status)
    terminal = {
        "schema": SCHEMA,
        "state": "completed",
        "accepted_campaign_result": False,
        "diagnostic_passed": receipt["diagnostic_passed"],
        "route": ROUTE,
        "detail": None,
        "dead_man_finalized": False,
        "completed_point_count": 9,
        "current_point": None,
        "preflight_sha256": status["preflight_sha256"],
        "preflight_identity": status["preflight_identity"],
        "execution_identity_sha256": status["execution_identity_sha256"],
        "contract_sha256": status["contract_sha256"],
        "source_receipt_sha256": status["source_receipt_sha256"],
        "source_evidence_sha256": status["source_evidence_sha256"],
        "source_snapshot": status["source_snapshot"],
        "receipt_sha256": receipt_sha,
        "status": desired_status,
        "status_sha256": _sha(desired_raw),
        "started_utc": status["started_utc"],
        "finished_utc": receipt["finished_utc"],
    }
    return terminal, desired_raw


def _failed_status_from_terminal(
    status: dict[str, Any], terminal: dict[str, Any]
) -> dict[str, Any]:
    return {
        **_running_projection(status),
        "state": "failed",
        "diagnostic_passed": False,
        "error": terminal["detail"],
        "finished_utc": terminal["finished_utc"],
    }


def _validate_terminal_bindings(
    terminal: dict[str, Any],
    status: dict[str, Any],
    preflight: dict[str, Any],
    preflight_raw: bytes,
) -> None:
    expected = {
        "preflight_sha256": _sha(preflight_raw),
        "preflight_identity": preflight["preflight_identity"],
        "execution_identity_sha256": preflight["execution_identity_sha256"],
        "contract_sha256": campaign._canonical_hash(CONTRACT),
        "source_receipt_sha256": SOURCE_RECEIPT_SHA256,
        "source_evidence_sha256": campaign._canonical_hash(
            preflight["source_evidence"]
        ),
        "source_snapshot": preflight["execution_identity"]["source_snapshot"],
        "started_utc": status["started_utc"],
    }
    for field, expected_value in expected.items():
        campaign._strict_json_equal(
            terminal.get(field), expected_value, label=f"terminal {field}"
        )
    if (
        terminal["completed_point_count"] != len(status["completed_points"])
        or terminal["current_point"] != status["current_point"]
    ):
        raise ValueError("scalar terminal point status binding is invalid")


def _status_is_allowed_terminal_ancestor(
    observed: dict[str, Any], desired: dict[str, Any]
) -> bool:
    if observed == desired:
        return True
    if observed.get("state") != "running" or desired.get("state") not in {
        "completed",
        "failed",
    }:
        return False
    changing = {
        "state",
        "diagnostic_passed",
        "center_density",
        "completed_points",
        "current_point",
        "receipt_sha256",
        "error",
        "finished_utc",
    }
    common = _status_common_keys() - changing
    if any(observed.get(key) != desired.get(key) for key in common):
        return False
    observed_points = observed["completed_points"]
    desired_points = desired["completed_points"]
    same_progress = (
        observed_points == desired_points
        and observed["current_point"] == desired["current_point"]
        and observed["center_density"] == desired["center_density"]
    )
    current_advance = (
        len(observed_points) < len(POINT_ORDER)
        and observed_points == desired_points
        and observed["current_point"] is None
        and desired["current_point"] == POINT_ORDER[len(observed_points)]
        and observed["center_density"] == desired["center_density"]
    )
    center_artifact_advance = (
        observed_points == desired_points == []
        and observed["current_point"] == desired["current_point"] == "center"
        and observed["center_density"] is None
        and desired["center_density"] is not None
    )
    point_commit_advance = (
        len(desired_points) == len(observed_points) + 1
        and desired_points[:-1] == observed_points
        and observed["current_point"] == desired_points[-1]["point_key"]
        and desired["current_point"] is None
        and (
            observed["center_density"] == desired["center_density"]
            or (
                observed["center_density"] is None
                and desired_points[-1]["point_key"] == "center"
                and desired["center_density"] is not None
            )
        )
    )
    return (
        same_progress
        or current_advance
        or center_artifact_advance
        or point_commit_advance
    )


def _expected_root_entries(
    status: dict[str, Any], *, terminal: bool, receipt: bool
) -> set[str]:
    entries = {"preflight.json", "status.json"}
    if status["center_density"] is not None:
        entries.add("center-density.f64")
    entries.update(
        f"point-{binding['point_key']}.json" for binding in status["completed_points"]
    )
    if terminal:
        entries.add("terminal.json")
    if receipt:
        entries.add("receipt.json")
    return entries


def _assert_root_entries(
    claim: hessian_diagnostic._OutputRootClaim,
    status: dict[str, Any],
    *,
    terminal: bool,
    receipt: bool,
) -> None:
    assert claim.root_descriptor is not None
    observed = set(os.listdir(claim.root_descriptor))
    expected = _expected_root_entries(status, terminal=terminal, receipt=receipt)
    if observed != expected:
        raise ValueError("scalar diagnostic root contains malformed or orphan evidence")


def _recover_running_status_ancestry(
    claim: hessian_diagnostic._OutputRootClaim,
    status: dict[str, Any],
    preflight: dict[str, Any],
    preflight_raw: bytes,
    evidence: dict[str, Any],
) -> dict[str, Any]:
    """Recover only one exact publication-boundary advance from a running status."""

    if status["state"] != "running":
        return status
    assert claim.root_descriptor is not None
    expected = _expected_root_entries(status, terminal=False, receipt=False)
    observed = set(os.listdir(claim.root_descriptor))
    extras = observed - expected
    if not extras:
        return status
    key = status["current_point"]
    if key is None:
        raise ValueError("scalar diagnostic root contains malformed or orphan evidence")
    point_name = f"point-{key}.json"
    allowed = {point_name}
    if key == "center" and status["center_density"] is None:
        allowed.add("center-density.f64")
    if not extras <= allowed or observed - extras != expected:
        raise ValueError("scalar diagnostic root contains malformed or orphan evidence")
    candidate = {**status}
    if "center-density.f64" in extras:
        density_raw = campaign._read_bounded_regular_snapshot(
            claim.path / "center-density.f64",
            label="recoverable scalar center density",
            maximum_bytes=128 * 1024 * 1024,
        )
        if len(density_raw) % 16:
            raise ValueError("recoverable scalar center density shape is invalid")
        orbital_count = math.isqrt(len(density_raw) // 16)
        if orbital_count < 1 or 16 * orbital_count**2 != len(density_raw):
            raise ValueError("recoverable scalar center density shape is invalid")
        candidate["center_density"] = {
            "path": "center-density.f64",
            "sha256": _sha(density_raw),
            "shape": [2, orbital_count, orbital_count],
            "dtype": "little-endian float64",
            "units": "electrons",
        }
    if point_name in extras:
        point, point_raw = hessian_diagnostic._claimed_json_object(
            claim,
            claim.path / point_name,
            label=f"recoverable scalar point {key}",
            max_bytes=256 * 1024,
        )
        if key == "center" and candidate["center_density"] is None:
            candidate["center_density"] = _validate_center_density_metadata(
                point.get("center_density")
            )
        candidate["completed_points"] = [
            *status["completed_points"],
            {"point_key": key, "receipt_sha256": _sha(point_raw)},
        ]
        candidate["current_point"] = None
    candidate_raw = _json_bytes(candidate)
    candidate = _validate_status(candidate, candidate_raw, preflight, preflight_raw)
    _read_validated_points(claim, candidate, evidence)
    _assert_root_entries(claim, candidate, terminal=False, receipt=False)
    return candidate


def finalize_if_running(output_root: Path) -> dict[str, Any]:
    if (
        _assert_resident_bindings is not _ATTESTED_RESIDENT_ASSERT
        or _assert_resident_bindings.__code__ is not _ATTESTED_RESIDENT_ASSERT_CODE
    ):
        raise RuntimeError("resident-boundary validator was replaced")
    _ATTESTED_RESIDENT_ASSERT()
    output = campaign._safe_absolute_root(output_root)
    with hessian_diagnostic._exclusive_output_claim(output) as claim:
        if claim.root_created:
            raise FileNotFoundError("scalar diagnostic output root did not exist")
        preflight, preflight_raw = _validated_preflight(claim)
        _utc_timestamp(preflight["created_utc"], label="preflight created_utc")
        current_identity = _current_identity()
        campaign._strict_json_equal(
            current_identity,
            preflight["execution_identity"],
            label="finalizer resident execution identity",
        )
        evidence = _load_source_evidence(SOURCE_ROOT)
        campaign._strict_json_equal(
            _public_evidence(evidence),
            preflight["source_evidence"],
            label="finalizer source evidence",
        )
        status, status_raw = hessian_diagnostic._claimed_json_object(
            claim,
            output / "status.json",
            label="scalar diagnostic status",
            max_bytes=4 * 1024 * 1024,
        )
        status = _validate_status(status, status_raw, preflight, preflight_raw)
        terminal_exists = hessian_diagnostic._claimed_root_entry_exists(
            claim, output / "terminal.json"
        )
        receipt_exists = hessian_diagnostic._claimed_root_entry_exists(
            claim, output / "receipt.json"
        )
        if terminal_exists:
            terminal, terminal_raw = hessian_diagnostic._claimed_json_object(
                claim,
                output / "terminal.json",
                label="scalar diagnostic terminal",
                max_bytes=1024 * 1024,
            )
            if _json_bytes(terminal) != terminal_raw:
                raise ValueError(
                    "scalar diagnostic terminal is not strict canonical JSON"
                )
            terminal = _validate_terminal(terminal)
            terminal_status_raw = _json_bytes(terminal["status"])
            terminal_status = _validate_status(
                terminal["status"], terminal_status_raw, preflight, preflight_raw
            )
            if not _status_is_allowed_terminal_ancestor(status, terminal_status):
                raise ValueError("scalar terminal is not a valid status commit")
            status = terminal_status
            _validate_terminal_bindings(terminal, status, preflight, preflight_raw)
            if terminal["state"] == "completed":
                if not receipt_exists:
                    raise ValueError("completed terminal is orphaned from its receipt")
                if len(status["completed_points"]) != 9:
                    raise ValueError("completed terminal is orphaned from nine points")
                receipt_payload, receipt_raw, _ = _validate_receipt(
                    claim, status, preflight, preflight_raw, evidence
                )
                receipt_sha = _sha(receipt_raw)
                desired_terminal, desired_status_raw = _completed_terminal(
                    status, receipt_payload, receipt_sha
                )
                if terminal != desired_terminal:
                    raise ValueError("completed scalar terminal is not reproducible")
            else:
                if receipt_exists:
                    raise ValueError("failed terminal has an orphan completed receipt")
                _read_validated_points(claim, status, evidence)
                desired_status_raw = terminal_status_raw
                if terminal["status_sha256"] != _sha(desired_status_raw):
                    raise ValueError("failed scalar terminal status binding is invalid")
            if _sha(desired_status_raw) != terminal["status_sha256"]:
                raise ValueError("scalar terminal status hash is invalid")
            if status_raw != desired_status_raw:
                hessian_diagnostic._claimed_write(
                    claim, output / "status.json", desired_status_raw
                )
                status = json.loads(desired_status_raw)
            _assert_root_entries(
                claim,
                status,
                terminal=True,
                receipt=terminal["state"] == "completed",
            )
            return terminal

        if receipt_exists:
            if (
                status["state"] != "running"
                or len(status["completed_points"]) != 9
                or status["current_point"] is not None
            ):
                raise ValueError("receipt is orphaned from complete running ancestry")
            receipt_payload, receipt_raw, _ = _validate_receipt(
                claim, status, preflight, preflight_raw, evidence
            )
            receipt_sha = _sha(receipt_raw)
            terminal, desired_status_raw = _completed_terminal(
                status, receipt_payload, receipt_sha
            )
            if _current_identity() != current_identity:
                raise RuntimeError("identity changed before recovered terminal commit")
            hessian_diagnostic._terminal_write_noreplace(
                claim, output / "terminal.json", terminal
            )
            hessian_diagnostic._claimed_write(
                claim, output / "status.json", desired_status_raw
            )
            repaired_status = json.loads(desired_status_raw)
            _assert_root_entries(claim, repaired_status, terminal=True, receipt=True)
            return _validate_terminal(terminal)

        if status["state"] == "completed":
            raise ValueError(
                "completed status is orphaned from its receipt and terminal"
            )
        status = _recover_running_status_ancestry(
            claim, status, preflight, preflight_raw, evidence
        )
        status_raw = _json_bytes(status)
        _read_validated_points(claim, status, evidence)
        _assert_root_entries(claim, status, terminal=False, receipt=False)
        if status["state"] == "failed":
            failed_terminal = {
                "schema": SCHEMA,
                "state": "failed",
                "accepted_campaign_result": False,
                "diagnostic_passed": False,
                "route": ROUTE,
                "detail": status["error"],
                "dead_man_finalized": True,
                "completed_point_count": len(status["completed_points"]),
                "current_point": status["current_point"],
                "preflight_sha256": status["preflight_sha256"],
                "preflight_identity": status["preflight_identity"],
                "execution_identity_sha256": status["execution_identity_sha256"],
                "contract_sha256": status["contract_sha256"],
                "source_receipt_sha256": status["source_receipt_sha256"],
                "source_evidence_sha256": status["source_evidence_sha256"],
                "source_snapshot": status["source_snapshot"],
                "status": status,
                "status_sha256": _sha(status_raw),
                "started_utc": status["started_utc"],
                "finished_utc": status["finished_utc"],
            }
            if _current_identity() != current_identity:
                raise RuntimeError(
                    "identity changed before recovered failed terminal commit"
                )
            hessian_diagnostic._terminal_write_noreplace(
                claim, output / "terminal.json", failed_terminal
            )
            return _validate_terminal(failed_terminal)
        return _write_failed_terminal(
            claim,
            status,
            RuntimeError("dead-man finalized an interrupted scalar diagnostic"),
            dead_man=True,
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--create-preflight", action="store_true")
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--finalize-if-running", action="store_true")
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--nice", type=int, default=10)
    parser.add_argument("--log")
    args = parser.parse_args()
    if not 1 <= args.threads <= 16:
        parser.error("--threads must be between 1 and 16")
    if args.nice < 10:
        parser.error("--nice must be at least 10")
    if args.create_preflight:
        if args.source_root is None:
            parser.error("--source-root is required with --create-preflight")
        payload = create_preflight(args.source_root, args.output_root)
        result = {
            "state": payload["state"],
            "preflight": str(args.output_root / "preflight.json"),
            "accepted_campaign_result": False,
        }
    elif args.run:
        if args.source_root is not None:
            parser.error("--source-root is valid only with --create-preflight")
        payload = run(args.output_root)
        result = {
            "state": payload["state"],
            "receipt": str(args.output_root / "receipt.json"),
            "diagnostic_passed": payload["diagnostic_passed"],
            "accepted_campaign_result": False,
        }
    else:
        if args.source_root is not None:
            parser.error("--source-root is valid only with --create-preflight")
        result = finalize_if_running(args.output_root)
    print(json.dumps(result, sort_keys=True, allow_nan=False))
    return 0


_RESIDENT_LOCAL_CALLABLES = {
    name: value
    for name, value in globals().items()
    if inspect.isfunction(value) and value.__module__ == __name__
}
_RESIDENT_LOCAL_CODE = {
    name: value.__code__ for name, value in _RESIDENT_LOCAL_CALLABLES.items()
}
_ATTESTED_RESIDENT_ASSERT = _assert_resident_bindings
_ATTESTED_RESIDENT_ASSERT_CODE = _assert_resident_bindings.__code__


if __name__ == "__main__":
    raise SystemExit(main())
