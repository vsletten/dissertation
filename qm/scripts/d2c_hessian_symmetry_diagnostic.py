#!/usr/bin/env python
"""Reproducible Hessian numerical diagnostics for D2c's frozen TS.

This command is evidence-only: it never publishes a qualification receipt or an
accepted campaign result.  It preserves raw component matrices, evaluates the
strict pre-symmetrization gate, and reports the post-symmetrization vibrational
verdict for the original four analytic cases.  The explicit finite-difference
confirmation mode instead evaluates 72 dense-grid displaced gradients against
the hash-pinned analytic reference while preserving the same evidence-only rule.
"""

from __future__ import annotations

# Import order is part of the provenance witness bootstrap below.
# ruff: noqa: E402, I001

import argparse
import ast
import fcntl
import hashlib
import json
import math
import os
import stat
import subprocess
import sys
import time
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import asdict, replace
from pathlib import Path
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
    if not 1 <= _guard_args.threads <= 16:
        _guard.error("--threads must be <= 16 and at least 1")
    if _guard_args.nice < 10:
        _guard.error("--nice must be >= 10")
    _ETIQUETTE = bootstrap_cli(
        "d2c_hessian_symmetry_diagnostic",
        default_run_root=Path("/mnt/data/vsletten/dissertation-data"),
    )

# Import the campaign first: it installs the execution witness while loading every
# module covered by the production-boundary manifest.
from scripts import d2c_sct_campaign as campaign  # noqa: E402

import numpy as np  # noqa: E402

from quarry import native_hessian  # noqa: E402
from quarry.native_hessian import HARTREE_TO_EV, NativeHessianResult  # noqa: E402
from quarry.pipeline import (  # noqa: E402
    BOHR_TO_ANGSTROM,
    _gradient_method,
    _hessian_method,
    _make_scf,
    build_mol,
    frequency_geometry_fingerprint,
)
from quarry.reaction_path import (  # noqa: E402
    hessian_eigenvalues_to_wavenumbers_cm,
    project_vibrational_hessian,
)
from scripts.surface_rate_protocol import reactions  # noqa: E402

campaign._OBSERVED_MODULE_CODE.setdefault(__name__, _THIS_MODULE_IMPORT_CODE)

SCHEMA = "d2c-hessian-symmetry-diagnostic-v1"
FINITE_DIFFERENCE_SCHEMA = "d2c-hessian-fd-confirmation-v2"
HALF_STEP_SCHEMA = "d2c-hessian-fd-half-step-extension-v1"
TERMINAL = "terminal.json"
ROUTE = "h-co-1w-oside"
PRIOR_FD_ROOT = (
    "/mnt/data/vsletten/dissertation-data/task300-d2c-fd-0172406-confirmation"
)
PRIOR_FD_RECEIPT_SHA256 = (
    "1bde3b52f516dbdb8df2ef4ea5d06233a783edef91c961c2ba943837abc15554"
)
PRIOR_FD_TERMINAL_SHA256 = (
    "d11d73af132ca35909752be19c41105761638b0284f33e00cc47ef33165a5a7e"
)
PRIOR_FD_STATUS_SHA256 = (
    "151015b0eb118b0e592a246b2b19602b1e4c4110697f387c91e3e3486258d5c8"
)
PRIOR_FD_GIT_SHA = "0172406b7df089d2e2da35955e4a2a65a3f81065"
PRIOR_FD_SCRIPT_SHA256 = (
    "f5bc13acff0eed4b708ad4405ab8f2fb78f8bb4763470ca9509644a7b2428d8e"
)
PRIOR_FD_PREFLIGHT_ROOT = (
    "/mnt/data/vsletten/dissertation-data/task300-d2c-fd-0172406-preflight"
)
PRIOR_FD_PREFLIGHT_SHA256 = (
    "48383c18f8005fe949d611fcf34d01d40a2b80bcae43bae29c730883613077ea"
)
PRIOR_FD_REFERENCE_ROOT = (
    "/mnt/data/vsletten/dissertation-data/task300-d2c-hessian-5f74662-diagnostic"
)
REFERENCE_RECEIPT_SHA256 = (
    "b32f0a34939571679e653e4772c1a73172126acd7f58fdbae3ea425ef463e172"
)
REFERENCE_MATRIX_SHA256 = (
    "7994a5916b6b0232042d0911911ca1c68cf2072a16af7f2e5482eaf9a6834d8c"
)
FINITE_DIFFERENCE_CONTRACT: dict[str, Any] = {
    "route": ROUTE,
    "method": {
        "backend": "pyscf-cpu",
        "xc": "pwb6k",
        "dispersion": "d3bj",
        "basis": "def2-svp",
        "density_fit": True,
        "grid_level": 5,
        "scf_tolerance": 1.0e-12,
        "scf_max_cycle": 150,
        "gradient_grid_response": True,
        "density_initial_guess": "one identical center converged density matrix",
    },
    "stencil": {
        "coordinate_order": "atom-major Cartesian x,y,z",
        "orientation": "rows are gradient components; columns are displacements",
        "inner_step_bohr": 0.01,
        "outer_step_bohr": 0.02,
        "displacements_per_coordinate": [-2, -1, 1, 2],
        "displaced_gradient_count": 72,
        "central_gradient_count": 1,
        "maximum_total_scf_gradient_evaluations": 73,
        "scf_retry_count": 0,
        "h_h": "(g(+h)-g(-h))/(2h)",
        "h_2h": "(g(+2h)-g(-2h))/(4h)",
        "richardson": "(4*H_h-H_2h)/3",
    },
    "reference": {
        "receipt_sha256": REFERENCE_RECEIPT_SHA256,
        "dense_strict_symmetric_matrix_sha256": REFERENCE_MATRIX_SHA256,
        "case": "D-dense-strict-reference",
    },
    "thresholds": {
        "required_grid_point_count": 199560,
        "spin_square_inclusive_minimum": 0.74,
        "spin_square_inclusive_maximum": 0.80,
        "displaced_vs_center_spin_square_exclusive_maximum_delta": 0.01,
        (
            "displaced_final_density_vs_center_exclusive_maximum_"
            "normalized_frobenius_delta"
        ): 0.05,
        "center_energy_vs_reference_exclusive_maximum_delta_hartree": 1.0e-8,
        "center_s2_vs_reference_exclusive_maximum_delta": 1.0e-4,
        "center_fmax_ev_per_angstrom_exclusive_maximum": 0.02,
        "raw_fd_maximum_absolute_asymmetry_exclusive_maximum": 2.0e-5,
        "raw_fd_spectral_relative_asymmetry_exclusive_maximum": 1.0e-5,
        "symmetric_richardson_vs_h_h_maximum_absolute_delta_exclusive_maximum": 5.0e-5,
        "symmetric_richardson_vs_h_h_spectral_relative_delta_exclusive_maximum": 2.0e-5,
        "fd_vs_analytic_symmetric_maximum_absolute_delta_exclusive_maximum": 2.0e-4,
        "fd_vs_analytic_symmetric_spectral_relative_delta_exclusive_maximum": 1.0e-4,
        "negative_eigenvalue_threshold": -1.0e-8,
        "near_zero_eigenvalue_inclusive_maximum_absolute": 1.0e-8,
        "lowest_positive_eigenvalue_exclusive_minimum": 5.0e-6,
        "imaginary_wavenumber_cm_exclusive_minimum": 200.0,
        "per_mode_reference_eigenvalue_delta_absolute_floor": 1.0e-6,
        "per_mode_reference_eigenvalue_delta_relative_fraction": 0.005,
        "per_mode_reference_eigenvalue_delta_comparison": (
            "exclusive maximum of floor and fraction*abs(reference)"
        ),
        "fd_vs_reference_imaginary_frequency_delta_cm_exclusive_maximum": 5.0,
        "fd_vs_reference_low_positive_frequency_delta_cm_exclusive_maximum": 1.0,
        "maximum_overlap_assignment_must_equal_eigenvalue_order": True,
        "fd_vs_reference_unstable_overlap_exclusive_minimum": 0.995,
        "fd_vs_reference_low_positive_overlap_exclusive_minimum": 0.98,
        "fd_vs_reference_all_ordered_overlaps_exclusive_minimum": 0.90,
        "fd_mapped_reaction_overlap_exclusive_minimum": 0.5,
        "fd_mapped_reaction_overlap_vs_reference_exclusive_maximum_delta": 0.02,
        "h_h_vs_fd_unstable_overlap_exclusive_minimum": 0.999,
        "h_h_vs_fd_low_positive_overlap_exclusive_minimum": 0.995,
        "h_h_vs_fd_imaginary_frequency_delta_cm_exclusive_maximum": 2.0,
        "h_h_vs_fd_low_positive_frequency_delta_cm_exclusive_maximum": 0.5,
    },
    "policy": {
        "accepted_campaign_result": False,
        "confirmation_only": True,
        "passed_authority": (
            "narrow route-specific PySCF/PWB6K-D3BJ/def2-SVP/grid-5 "
            "production-policy proposal only; never qualification"
        ),
        "rejected_authority": (
            "diagnostic rejection only; no tolerance change, TS qualification, "
            "IRC launch, SCT result, or production-policy change"
        ),
    },
}

HALF_STEP_CONTRACT: dict[str, Any] = {
    "route": ROUTE,
    "method": {
        "backend": "pyscf-cpu",
        "xc": "pwb6k",
        "dispersion": "d3bj",
        "basis": "def2-svp",
        "density_fit": True,
        "grid_level": 5,
        "scf_tolerance": 1.0e-12,
        "scf_max_cycle": 150,
        "gradient_grid_response": True,
        "density_initial_guess": "prior receipt's identical center density matrix",
    },
    "stencil": {
        "coordinate_order": "atom-major Cartesian x,y,z",
        "orientation": "rows are gradient components; columns are displacements",
        "inner_step_bohr": 0.005,
        "outer_step_bohr": 0.01,
        "new_displacements_per_coordinate": [-1, 1],
        "reused_prior_displacements_per_coordinate": [-1, 1],
        "new_displaced_gradient_count": 36,
        "reused_displaced_gradient_count": 36,
        "maximum_new_scf_gradient_evaluations": 36,
        "scf_retry_count": 0,
        "h_inner": "(g(+0.005)-g(-0.005))/0.01",
        "h_outer": "(g(+0.01)-g(-0.01))/0.02",
        "richardson": "(4*H_0.005-H_0.01)/3",
    },
    "prior_fd": {
        "root": PRIOR_FD_ROOT,
        "schema": FINITE_DIFFERENCE_SCHEMA,
        "receipt_sha256": PRIOR_FD_RECEIPT_SHA256,
        "terminal_sha256": PRIOR_FD_TERMINAL_SHA256,
        "status_sha256": PRIOR_FD_STATUS_SHA256,
        "git_sha": PRIOR_FD_GIT_SHA,
        "script_sha256": PRIOR_FD_SCRIPT_SHA256,
        "preflight_root": PRIOR_FD_PREFLIGHT_ROOT,
        "preflight_receipt_sha256": PRIOR_FD_PREFLIGHT_SHA256,
        "reference_root": PRIOR_FD_REFERENCE_ROOT,
    },
    "reference": {
        "receipt_sha256": REFERENCE_RECEIPT_SHA256,
        "dense_strict_symmetric_matrix_sha256": REFERENCE_MATRIX_SHA256,
        "case": "D-dense-strict-reference",
    },
    "thresholds": {
        "required_grid_point_count": 199560,
        "spin_square_inclusive_minimum": 0.74,
        "spin_square_inclusive_maximum": 0.80,
        "displaced_vs_center_spin_square_exclusive_maximum_delta": 0.01,
        (
            "displaced_final_density_vs_center_exclusive_maximum_"
            "normalized_frobenius_delta"
        ): 0.05,
        "center_energy_vs_reference_exclusive_maximum_delta_hartree": 1.0e-8,
        "center_s2_vs_reference_exclusive_maximum_delta": 1.0e-4,
        "center_fmax_ev_per_angstrom_exclusive_maximum": 0.02,
        "raw_fd_maximum_absolute_asymmetry_exclusive_maximum": 2.0e-5,
        "raw_fd_spectral_relative_asymmetry_exclusive_maximum": 1.0e-5,
        "symmetric_richardson_vs_h_h_maximum_absolute_delta_exclusive_maximum": 5.0e-5,
        "symmetric_richardson_vs_h_h_spectral_relative_delta_exclusive_maximum": 2.0e-5,
        "fd_vs_analytic_symmetric_maximum_absolute_delta_exclusive_maximum": 2.0e-4,
        "fd_vs_analytic_symmetric_spectral_relative_delta_exclusive_maximum": 1.0e-4,
        "negative_eigenvalue_threshold": -1.0e-8,
        "near_zero_eigenvalue_inclusive_maximum_absolute": 1.0e-8,
        "lowest_positive_eigenvalue_exclusive_minimum": 5.0e-6,
        "imaginary_wavenumber_cm_exclusive_minimum": 200.0,
        "per_mode_reference_eigenvalue_delta_absolute_floor": 1.0e-6,
        "per_mode_reference_eigenvalue_delta_relative_fraction": 0.005,
        "per_mode_reference_eigenvalue_delta_comparison": (
            "exclusive maximum of floor and fraction*abs(reference)"
        ),
        "fd_vs_reference_imaginary_frequency_delta_cm_exclusive_maximum": 5.0,
        "fd_vs_reference_low_positive_frequency_delta_cm_exclusive_maximum": 1.0,
        "maximum_overlap_assignment_must_equal_eigenvalue_order": True,
        "fd_vs_reference_unstable_overlap_exclusive_minimum": 0.995,
        "fd_vs_reference_low_positive_overlap_exclusive_minimum": 0.98,
        "fd_vs_reference_all_ordered_overlaps_exclusive_minimum": 0.90,
        "fd_mapped_reaction_overlap_exclusive_minimum": 0.5,
        "fd_mapped_reaction_overlap_vs_reference_exclusive_maximum_delta": 0.02,
        "h_h_vs_fd_unstable_overlap_exclusive_minimum": 0.999,
        "h_h_vs_fd_low_positive_overlap_exclusive_minimum": 0.995,
        "h_h_vs_fd_imaginary_frequency_delta_cm_exclusive_maximum": 2.0,
        "h_h_vs_fd_low_positive_frequency_delta_cm_exclusive_maximum": 0.5,
    },
    "policy": {
        "accepted_campaign_result": False,
        "confirmation_only": True,
        "passed_authority": (
            "bounded half-step diagnostic evidence only; never qualification"
        ),
        "rejected_authority": (
            "diagnostic rejection only; no tolerance change, TS qualification, "
            "IRC launch, SCT result, D3b authorization, or production-policy change"
        ),
    },
}


class ScientificRejection(ValueError):
    """A predeclared scientific gate rejected otherwise valid evidence."""


class ConcurrentRunError(RuntimeError):
    """Another process owns a compatible run claim; do not touch its output."""


CASES = (
    {
        "name": "A-baseline",
        "grid_level": 3,
        "scf_tolerance": 1.0e-9,
        "cpscf_tolerance": 1.0e-8,
        "hessian_max_cycle": 50,
    },
    {
        "name": "B-strict-response",
        "grid_level": 3,
        "scf_tolerance": 1.0e-12,
        "cpscf_tolerance": 1.0e-12,
        "hessian_max_cycle": 150,
    },
    {
        "name": "C-dense-grid",
        "grid_level": 5,
        "scf_tolerance": 1.0e-9,
        "cpscf_tolerance": 1.0e-8,
        "hessian_max_cycle": 50,
    },
    {
        "name": "D-dense-strict-reference",
        "grid_level": 5,
        "scf_tolerance": 1.0e-12,
        "cpscf_tolerance": 1.0e-12,
        "hessian_max_cycle": 150,
    },
)


def _json_bytes(payload: dict[str, Any]) -> bytes:
    return (
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode()


def _atomic_write(
    path: Path, raw: bytes, *, claim: _OutputRootClaim | None = None
) -> None:
    temporary = path.parent / f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp"
    with temporary.open("xb") as handle:
        handle.write(raw)
        handle.flush()
        os.fsync(handle.fileno())
    if claim is not None:
        claim.verify()
    os.replace(temporary, path)
    directory_fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


class _OutputRootClaim:
    """External basename claim plus a continuously checked output inode binding."""

    def __init__(
        self, path: Path, parent_fd: int, descriptor: int, claim_name: str
    ) -> None:
        self.path = path
        self.parent = path.parent
        self.name = path.name
        self.parent_fd = parent_fd
        self.descriptor = descriptor
        self.claim_name = claim_name
        self.root_created = False
        claim_stat = os.fstat(descriptor)
        if not stat.S_ISREG(claim_stat.st_mode):
            raise ValueError("diagnostic external claim must be a regular file")
        self.claim_identity = (claim_stat.st_dev, claim_stat.st_ino)
        self.root_identity: tuple[int, int] | None = None

    def bind(self) -> None:
        try:
            os.mkdir(self.name, mode=0o700, dir_fd=self.parent_fd)
        except FileExistsError:
            self.root_created = False
        else:
            self.root_created = True
        root = os.stat(self.name, dir_fd=self.parent_fd, follow_symlinks=False)
        if not stat.S_ISDIR(root.st_mode):
            raise ValueError("diagnostic output root must be a regular directory")
        self.root_identity = (root.st_dev, root.st_ino)
        self.verify()

    def verify(self) -> None:
        parent_path = os.stat(self.parent, follow_symlinks=False)
        parent_fd = os.fstat(self.parent_fd)
        if not stat.S_ISDIR(parent_path.st_mode) or (
            parent_path.st_dev,
            parent_path.st_ino,
        ) != (parent_fd.st_dev, parent_fd.st_ino):
            raise RuntimeError("diagnostic output parent directory was replaced")
        try:
            claim_path = os.stat(
                self.claim_name, dir_fd=self.parent_fd, follow_symlinks=False
            )
        except FileNotFoundError as exc:
            raise RuntimeError("diagnostic external claim was removed") from exc
        claim_descriptor = os.fstat(self.descriptor)
        if (
            not stat.S_ISREG(claim_path.st_mode)
            or (claim_path.st_dev, claim_path.st_ino) != self.claim_identity
            or (claim_descriptor.st_dev, claim_descriptor.st_ino) != self.claim_identity
        ):
            raise RuntimeError("diagnostic external claim was replaced")
        if self.root_identity is None:
            return
        try:
            root = os.stat(self.name, dir_fd=self.parent_fd, follow_symlinks=False)
        except FileNotFoundError as exc:
            raise RuntimeError(
                "diagnostic output root was removed while claimed"
            ) from exc
        if (
            not stat.S_ISDIR(root.st_mode)
            or (root.st_dev, root.st_ino) != self.root_identity
        ):
            raise RuntimeError("diagnostic output root was replaced while claimed")
        path_root = os.stat(self.path, follow_symlinks=False)
        if (path_root.st_dev, path_root.st_ino) != self.root_identity:
            raise RuntimeError(
                "diagnostic output path no longer names the claimed root"
            )


@contextmanager
def _exclusive_output_claim(output_root: Path) -> Iterator[_OutputRootClaim]:
    """Claim a path outside its replaceable output directory and bind its inode."""

    output = campaign._safe_absolute_root(output_root)
    output.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    output = campaign._safe_absolute_root(output)
    parent_flags = (
        os.O_RDONLY
        | os.O_DIRECTORY
        | getattr(os, "O_NOFOLLOW", 0)
        | getattr(os, "O_CLOEXEC", 0)
    )
    parent_fd = os.open(output.parent, parent_flags)
    parent_locked = False
    try:
        try:
            fcntl.flock(parent_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise ConcurrentRunError(
                f"diagnostic stable output parent is already claimed: {output.parent}"
            ) from exc
        parent_locked = True
    except BaseException:
        os.close(parent_fd)
        raise
    claim_name = f".{output.name}.d2c-hessian-diagnostic.claim"
    claim_flags = (
        os.O_RDWR
        | os.O_CREAT
        | getattr(os, "O_NOFOLLOW", 0)
        | getattr(os, "O_CLOEXEC", 0)
    )
    try:
        descriptor = os.open(claim_name, claim_flags, 0o600, dir_fd=parent_fd)
    except BaseException:
        fcntl.flock(parent_fd, fcntl.LOCK_UN)
        os.close(parent_fd)
        raise
    locked = False
    try:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(
                f"diagnostic output path is already claimed: {output}"
            ) from exc
        locked = True
        claim = _OutputRootClaim(output, parent_fd, descriptor, claim_name)
        claim.bind()
        if claim.root_identity is None:
            raise RuntimeError("diagnostic output root was not bound")
        output_device, output_inode = claim.root_identity
        payload = _json_bytes(
            {
                "claimed_unix_ns": time.time_ns(),
                "output_basename": output.name,
                "output_st_dev": output_device,
                "output_st_ino": output_inode,
                "pid": os.getpid(),
            }
        )
        os.ftruncate(descriptor, 0)
        os.write(descriptor, payload)
        os.fsync(descriptor)
        claim.verify()
        yield claim
        claim.verify()
    finally:
        if locked:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        os.close(descriptor)
        if parent_locked:
            fcntl.flock(parent_fd, fcntl.LOCK_UN)
        os.close(parent_fd)


def _claimed_write(claim: _OutputRootClaim, path: Path, raw: bytes) -> None:
    claim.verify()
    _atomic_write(path, raw, claim=claim)
    claim.verify()


def _terminal_write_noreplace(
    claim: _OutputRootClaim, path: Path, payload: dict[str, Any]
) -> None:
    """Publish the authoritative terminal once; never replace prior evidence."""

    claim.verify()
    stage = path.parent / f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp"
    try:
        with stage.open("xb") as handle:
            handle.write(_json_bytes(payload))
            handle.flush()
            os.fsync(handle.fileno())
        claim.verify()
        campaign._renameat2_noreplace(stage, path)
        directory_fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        stage.unlink(missing_ok=True)
    claim.verify()


def _terminal_payload(
    status: dict[str, Any], *, state: str, detail: str | None
) -> dict[str, Any]:
    payload = {
        "schema": status.get("schema", SCHEMA),
        "state": state,
        "route": ROUTE,
        "git_sha": status.get("git_sha"),
        "script_sha256": status.get("script_sha256"),
        "preflight_receipt_sha256": status.get("preflight_receipt_sha256"),
        "completed_cases": status.get("completed_cases", []),
        "current_case": status.get("current_case"),
        "detail": detail,
        "finished_utc": status.get("finished_utc")
        or time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    if status.get("schema") in {FINITE_DIFFERENCE_SCHEMA, HALF_STEP_SCHEMA}:
        payload.update(
            {
                "accepted_campaign_result": False,
                "confirmation_passed": bool(
                    state == "completed" and status.get("confirmation_passed") is True
                ),
                "contract_sha256": status.get("contract_sha256"),
                "reference_receipt_sha256": status.get("reference_receipt_sha256"),
                "reference_matrix_sha256": status.get("reference_matrix_sha256"),
                "completed_point_count": len(status.get("completed_points", [])),
                "current_point": status.get("current_point"),
            }
        )
    if status.get("schema") == HALF_STEP_SCHEMA:
        payload.update(
            {
                "diagnostic_execution_identity_sha256": status.get(
                    "diagnostic_execution_identity_sha256"
                ),
                "prior_fd_receipt_sha256": status.get("prior_fd_receipt_sha256"),
                "prior_fd_terminal_sha256": status.get("prior_fd_terminal_sha256"),
                "prior_fd_status_sha256": status.get("prior_fd_status_sha256"),
            }
        )
    return payload


def _validated_terminal(path: Path) -> dict[str, Any]:
    terminal, _raw = campaign._read_json_object(
        path, label="diagnostic terminal receipt", max_bytes=1024 * 1024
    )
    schema = terminal.get("schema")
    base_keys = {
        "schema",
        "state",
        "route",
        "git_sha",
        "script_sha256",
        "preflight_receipt_sha256",
        "completed_cases",
        "current_case",
        "detail",
        "finished_utc",
    }
    if schema in {FINITE_DIFFERENCE_SCHEMA, HALF_STEP_SCHEMA}:
        base_keys |= {
            "accepted_campaign_result",
            "confirmation_passed",
            "contract_sha256",
            "reference_receipt_sha256",
            "reference_matrix_sha256",
            "completed_point_count",
            "current_point",
        }
        if schema == HALF_STEP_SCHEMA:
            base_keys |= {
                "diagnostic_execution_identity_sha256",
                "prior_fd_receipt_sha256",
                "prior_fd_terminal_sha256",
                "prior_fd_status_sha256",
            }
    elif schema != SCHEMA:
        raise ValueError("diagnostic terminal schema is invalid")
    if terminal.get("state") == "completed":
        base_keys |= {"receipt", "receipt_sha256"}
    if set(terminal) != base_keys:
        raise ValueError("diagnostic terminal fields are unexpected or incomplete")
    if terminal.get("route") != ROUTE or terminal.get("state") not in {
        "completed",
        "failed",
    }:
        raise ValueError("diagnostic terminal identity or state is invalid")
    if schema in {FINITE_DIFFERENCE_SCHEMA, HALF_STEP_SCHEMA} and (
        terminal.get("accepted_campaign_result") is not False
        or type(terminal.get("confirmation_passed")) is not bool
        or (
            terminal["state"] != "completed"
            and terminal["confirmation_passed"] is not False
        )
    ):
        raise ValueError("FD terminal scientific verdict is invalid")
    if schema == HALF_STEP_SCHEMA:
        for field, length in (
            ("git_sha", 40),
            ("script_sha256", 64),
            ("diagnostic_execution_identity_sha256", 64),
            ("contract_sha256", 64),
            ("preflight_receipt_sha256", 64),
            ("reference_receipt_sha256", 64),
            ("reference_matrix_sha256", 64),
            ("prior_fd_receipt_sha256", 64),
            ("prior_fd_terminal_sha256", 64),
            ("prior_fd_status_sha256", 64),
        ):
            campaign._require_sha(
                campaign._require_json_string(
                    terminal.get(field), label=f"half-step terminal {field}"
                ),
                length=length,
                label=f"half-step terminal {field}",
            )
        completed_count = terminal.get("completed_point_count")
        expected_point_keys = [item["key"] for item in _half_step_plan(6)]
        if (
            type(completed_count) is not int
            or not 0 <= completed_count <= len(expected_point_keys)
            or (
                terminal.get("current_point") is not None
                and terminal.get("current_point") not in expected_point_keys
            )
        ):
            raise ValueError("half-step terminal point progress is invalid")
    return terminal


def _validated_completed_fd_receipt(
    output_root: Path, *, verify_resident_identity: bool = True
) -> tuple[dict[str, Any], str]:
    receipt, raw = campaign._read_json_object(
        output_root / "receipt.json",
        label="completed FD confirmation receipt",
        max_bytes=4 * 1024 * 1024,
    )
    expected_keys = {
        "schema",
        "state",
        "accepted_campaign_result",
        "confirmation_passed",
        "purpose",
        "authority",
        "route",
        "git_sha",
        "script_sha256",
        "diagnostic_execution_identity",
        "diagnostic_execution_identity_sha256",
        "preflight_root",
        "preflight_receipt_sha256",
        "campaign_identity",
        "reference_root",
        "reference_receipt_sha256",
        "reference_matrix_sha256",
        "reference_binding",
        "reference_binding_sha256",
        "reference_git_sha",
        "reference_script_sha256",
        "contract",
        "contract_sha256",
        "point_receipts",
        "aggregate_receipt",
        "analysis",
        "finished_utc",
    }
    if set(receipt) != expected_keys:
        raise ValueError("completed FD receipt schema is not exact")
    execution_identity = receipt.get("diagnostic_execution_identity")
    reference_binding = receipt.get("reference_binding")
    if (
        type(execution_identity) is not dict
        or set(execution_identity)
        != {
            "kind",
            "loaded_module_code_sha256",
            "loaded_code_sha256",
            "loaded_state_sha256",
            "source_files",
        }
        or execution_identity.get("kind") != "python-source"
        or type(execution_identity.get("source_files")) is not dict
        or not execution_identity["source_files"]
    ):
        raise ValueError("completed FD loaded-code identity schema is invalid")
    for field in (
        "loaded_module_code_sha256",
        "loaded_code_sha256",
        "loaded_state_sha256",
    ):
        campaign._require_sha(
            campaign._require_json_string(
                execution_identity.get(field), label=f"FD {field}"
            ),
            length=64,
            label=f"FD {field}",
        )
    for path, digest in execution_identity["source_files"].items():
        if type(path) is not str or not path:
            raise ValueError("completed FD loaded-code source path is invalid")
        campaign._require_sha(
            campaign._require_json_string(digest, label="FD loaded source SHA-256"),
            length=64,
            label="FD loaded source SHA-256",
        )
    campaign._require_sha(
        campaign._require_json_string(receipt.get("git_sha"), label="FD Git SHA"),
        length=40,
        label="FD Git SHA",
    )
    campaign._require_sha(
        campaign._require_json_string(
            receipt.get("script_sha256"), label="FD script SHA-256"
        ),
        length=64,
        label="FD script SHA-256",
    )
    for field, length in (
        ("campaign_identity", 64),
        ("preflight_receipt_sha256", 64),
        ("reference_receipt_sha256", 64),
        ("reference_matrix_sha256", 64),
        ("reference_binding_sha256", 64),
        ("reference_git_sha", 40),
        ("reference_script_sha256", 64),
        ("contract_sha256", 64),
    ):
        campaign._require_sha(
            campaign._require_json_string(receipt.get(field), label=f"FD {field}"),
            length=length,
            label=f"FD {field}",
        )
    if (
        receipt.get("schema") != FINITE_DIFFERENCE_SCHEMA
        or receipt.get("state") != "completed"
        or receipt.get("accepted_campaign_result") is not False
        or type(receipt.get("confirmation_passed")) is not bool
        or receipt.get("route") != ROUTE
        or receipt.get("contract") != FINITE_DIFFERENCE_CONTRACT
        or receipt.get("contract_sha256")
        != campaign._canonical_hash(FINITE_DIFFERENCE_CONTRACT)
        or type(execution_identity) is not dict
        or receipt.get("diagnostic_execution_identity_sha256")
        != campaign._canonical_hash(execution_identity)
        or type(reference_binding) is not dict
        or receipt.get("reference_binding_sha256")
        != campaign._canonical_hash(reference_binding)
        or reference_binding.get("campaign_identity")
        != receipt.get("campaign_identity")
        or reference_binding.get("route") != ROUTE
        or reference_binding.get("fd_method") != FINITE_DIFFERENCE_CONTRACT["method"]
    ):
        raise ValueError("completed FD receipt identity is invalid")
    if verify_resident_identity:
        current_identity = _fd_execution_identity()
        campaign._strict_json_equal(
            {
                "git_sha": receipt.get("git_sha"),
                "script_sha256": receipt.get("script_sha256"),
                "diagnostic_execution_identity": execution_identity,
                "diagnostic_execution_identity_sha256": receipt.get(
                    "diagnostic_execution_identity_sha256"
                ),
            },
            current_identity,
            label="completed FD resident execution identity",
        )
    preflight_root_text = campaign._require_json_string(
        receipt.get("preflight_root"), label="completed FD preflight root"
    )
    reference_root_text = campaign._require_json_string(
        receipt.get("reference_root"), label="completed FD reference root"
    )
    preflight_root = campaign._safe_absolute_root(Path(preflight_root_text))
    reference_root = campaign._safe_absolute_root(Path(reference_root_text))
    if (
        str(preflight_root) != preflight_root_text
        or str(reference_root) != reference_root_text
    ):
        raise ValueError("completed FD provenance roots are not canonical")
    if verify_resident_identity:
        preflight, preflight_sha = campaign._validate_production_boundary(
            preflight_root, ROUTE
        )
    else:
        preflight, preflight_sha = campaign._validated_preflight(preflight_root, ROUTE)
    campaign._canonical_dft_settings(preflight)
    template = reactions(gpu=True, basis="def2-svp")[ROUTE].cluster
    fingerprints = campaign._trusted_input_fingerprints_from_route_record(
        preflight["routes"][ROUTE]
    )
    inputs = campaign._trusted_route_input_snapshots(
        campaign.DEFAULT_BUNDLE_ROOT, ROUTE, template, fingerprints
    )
    transition_state = inputs["transition_state"]
    if len(transition_state.symbols) != 6:
        raise ValueError("completed FD receipt does not bind the six-atom frozen TS")
    masses = np.asarray(preflight["routes"][ROUTE]["masses_amu"], dtype=float)
    (
        reference_receipt,
        analytic_matrix,
        reference_case,
        expected_reference_binding,
    ) = _load_reference(
        reference_root,
        6,
        preflight=preflight,
        transition_state=transition_state,
        reactant=inputs["reactant"],
        product=inputs["product"],
        masses=masses,
    )
    if (
        preflight_sha != receipt.get("preflight_receipt_sha256")
        or preflight.get("identity") != receipt.get("campaign_identity")
        or receipt.get("reference_receipt_sha256") != REFERENCE_RECEIPT_SHA256
        or receipt.get("reference_matrix_sha256") != REFERENCE_MATRIX_SHA256
        or receipt.get("reference_git_sha") != reference_receipt.get("git_sha")
        or receipt.get("reference_script_sha256")
        != reference_receipt.get("script_sha256")
    ):
        raise ValueError("completed FD external provenance binding is invalid")
    campaign._strict_json_equal(
        reference_binding,
        expected_reference_binding,
        label="completed FD reference binding",
    )
    points = receipt.get("point_receipts")
    if type(points) is not list or len(points) != 73:
        raise ValueError("completed FD receipt point ancestry is incomplete")
    expected_point_keys = [
        "center",
        *(item["key"] for item in _displacement_plan(6)),
    ]
    definitions = {item["key"]: item for item in _displacement_plan(6)}
    center_density: np.ndarray | None = None
    center_density_sha: str | None = None
    center_spin: float | None = None
    center_record: dict[str, Any] | None = None
    displaced_points: list[tuple[dict[str, Any], np.ndarray]] = []
    gradients: dict[tuple[int, int], np.ndarray] = {}
    for record, expected_key in zip(points, expected_point_keys, strict=True):
        if (
            type(record) is not dict
            or set(record) != {"point_key", "receipt_sha256"}
            or record.get("point_key") != expected_key
            or type(record.get("receipt_sha256")) is not str
            or len(record["receipt_sha256"]) != 64
            or any(
                character not in "0123456789abcdef"
                for character in record["receipt_sha256"]
            )
        ):
            raise ValueError("completed FD receipt point ancestry is invalid")
        receipt_keys = (
            _POINT_OBSERVATION_KEYS
            | _POINT_BINDING_KEYS
            | {"schema", "point_key", "kind", "artifacts"}
        )
        if expected_key == "center":
            receipt_keys |= _CENTER_REFERENCE_KEYS | {"physical_fmax_ev_per_angstrom"}
        else:
            receipt_keys |= set(definitions[expected_key]) | _DENSITY_CONTINUITY_KEYS
        point, point_sha, arrays = _read_record_bundle(
            output_root,
            Path("points") / expected_key,
            expected_receipt_keys=receipt_keys,
        )
        if point_sha != record["receipt_sha256"]:
            raise ValueError("completed FD point receipt SHA-256 hash mismatch")
        if (
            point.get("schema") != FINITE_DIFFERENCE_SCHEMA
            or point.get("point_key") != expected_key
            or any(
                point.get(field) != receipt.get(field) for field in _POINT_BINDING_KEYS
            )
            or set(arrays) != {"gradient.f64", "density.f64"}
            or arrays["gradient.f64"].shape != (6, 3)
            or point["artifacts"]["gradient.f64"]["units"] != "hartree / bohr"
            or point["artifacts"]["density.f64"]["units"] != "electrons"
        ):
            raise ValueError("completed FD point content is invalid")
        density = arrays["density.f64"]
        density_sha = hashlib.sha256(
            np.ascontiguousarray(density, dtype="<f8").tobytes()
        ).hexdigest()
        if point.get("density_final_sha256") != density_sha:
            raise ValueError("completed FD point final density hash mismatch")
        if expected_key == "center":
            if (
                point.get("kind") != "center"
                or density.ndim != 3
                or density.shape[0] != 2
                or density.shape[1] != density.shape[2]
                or point.get("density_initial_guess_sha256") != density_sha
                or point.get("geometry_fingerprint")
                != reference_binding.get("transition_state_geometry_fingerprint")
            ):
                raise ValueError("completed FD center point content is invalid")
            center_density = density
            center_density_sha = density_sha
            center_spin = _validate_point_record(
                point,
                center_spin_square=None,
                center=True,
                enforce_scientific_gates=False,
            )
            center_record = point
            gradient = arrays["gradient.f64"]
            observed_fmax = (
                float(np.max(np.linalg.norm(gradient, axis=1)))
                * HARTREE_TO_EV
                / BOHR_TO_ANGSTROM
            )
            if point.get("physical_fmax_ev_per_angstrom") != observed_fmax:
                raise ValueError("completed FD center fmax is not reproducible")
            campaign._strict_json_equal(
                {key: point.get(key) for key in _CENTER_REFERENCE_KEYS},
                _center_reference_metrics(point, reference_case),
                label="completed FD center reference metrics",
            )
        else:
            definition = definitions[expected_key]
            if (
                center_density is None
                or center_density_sha is None
                or point.get("kind") != "displacement"
                or any(point.get(field) != value for field, value in definition.items())
                or density.shape != center_density.shape
                or point.get("density_initial_guess_sha256") != center_density_sha
            ):
                raise ValueError("completed FD displaced point content is invalid")
            expected_coords = np.asarray(transition_state.coords, dtype=float).copy()
            expected_coords.reshape(-1)[definition["coordinate_index"]] += (
                definition["displacement_bohr"] * BOHR_TO_ANGSTROM
            )
            if point.get("geometry_fingerprint") != frequency_geometry_fingerprint(
                replace(transition_state, coords=expected_coords)
            ):
                raise ValueError("completed FD displaced geometry is invalid")
            _validate_point_record(
                point,
                center_spin_square=center_spin,
                center=False,
                enforce_scientific_gates=False,
            )
            expected_density_metrics = _density_continuity_metrics(
                density, center_density
            )
            campaign._strict_json_equal(
                {key: point.get(key) for key in _DENSITY_CONTINUITY_KEYS},
                expected_density_metrics,
                label="completed FD density continuity receipt",
            )
            displaced_points.append((point, density))
            gradients[
                (definition["coordinate_index"], definition["step_multiplier"])
            ] = arrays["gradient.f64"].reshape(-1)
    aggregate = receipt.get("aggregate_receipt")
    if (
        type(aggregate) is not dict
        or set(aggregate) != {"path", "sha256"}
        or aggregate.get("path") != "matrices/receipt.json"
    ):
        raise ValueError("completed FD aggregate ancestry is invalid")
    matrix_contract = {
        filename: ((18, 18), "hartree / bohr^2")
        for filename in (
            "H_h-raw.f64",
            "H_2h-raw.f64",
            "richardson-raw.f64",
            "H_h-symmetric.f64",
            "H_2h-symmetric.f64",
            "richardson-symmetric.f64",
        )
    }
    aggregate_payload, aggregate_sha, matrix_arrays = _read_record_bundle(
        output_root,
        Path("matrices"),
        expected_receipt_keys={
            "schema",
            "kind",
            "confirmation_passed",
            *_POINT_BINDING_KEYS,
            "point_receipts",
            "analysis",
            "artifacts",
        },
        expected_artifacts=matrix_contract,
    )
    if aggregate_sha != aggregate.get("sha256"):
        raise ValueError("completed FD aggregate receipt hash mismatch")
    if (
        aggregate_payload.get("schema") != FINITE_DIFFERENCE_SCHEMA
        or aggregate_payload.get("kind") != "finite-difference-Hessian-analysis"
        or aggregate_payload.get("confirmation_passed")
        is not receipt["confirmation_passed"]
        or aggregate_payload.get("point_receipts") != points
        or aggregate_payload.get("analysis") != receipt.get("analysis")
        or any(
            aggregate_payload.get(field) != receipt.get(field)
            for field in _POINT_BINDING_KEYS
        )
    ):
        raise ValueError("completed FD aggregate receipt content is invalid")
    if center_record is None or center_density is None:
        raise ValueError("completed FD center evidence is missing")
    h_h, h_2h, raw_fd = _finite_difference_matrices(gradients, 18)
    reconstructed = {
        "H_h-raw.f64": h_h,
        "H_2h-raw.f64": h_2h,
        "richardson-raw.f64": raw_fd,
        "H_h-symmetric.f64": 0.5 * (h_h + h_h.T),
        "H_2h-symmetric.f64": 0.5 * (h_2h + h_2h.T),
        "richardson-symmetric.f64": 0.5 * (raw_fd + raw_fd.T),
    }
    for filename, expected_matrix in reconstructed.items():
        if not np.array_equal(matrix_arrays[filename], expected_matrix):
            raise ValueError(
                "completed FD matrix diverges from 73-gradient reconstruction: "
                f"{filename}"
            )
    route_vector = campaign._mapped_route_vector(
        transition_state,
        masses,
        inputs["reactant"].coords,
        inputs["product"].coords,
    )
    confirmation_passed, expected_analysis = _scientific_analysis_verdict(
        transition_state,
        masses,
        route_vector,
        center_record,
        center_density,
        displaced_points,
        h_h,
        h_2h,
        raw_fd,
        analytic_matrix,
    )
    if confirmation_passed is not receipt["confirmation_passed"]:
        raise ValueError("completed FD scientific verdict is not reproducible")
    campaign._strict_json_equal(
        receipt.get("analysis"),
        expected_analysis,
        label="completed FD independently recomputed analysis",
    )
    return receipt, hashlib.sha256(raw).hexdigest()


def _completed_fd_status(
    receipt: dict[str, Any], receipt_sha256: str
) -> dict[str, Any]:
    return {
        "schema": FINITE_DIFFERENCE_SCHEMA,
        "state": "completed",
        "route": ROUTE,
        "git_sha": receipt["git_sha"],
        "script_sha256": receipt["script_sha256"],
        "diagnostic_execution_identity": receipt["diagnostic_execution_identity"],
        "diagnostic_execution_identity_sha256": receipt[
            "diagnostic_execution_identity_sha256"
        ],
        "preflight_root": receipt["preflight_root"],
        "reference_root": receipt["reference_root"],
        "contract": receipt["contract"],
        "contract_sha256": receipt["contract_sha256"],
        "preflight_receipt_sha256": receipt["preflight_receipt_sha256"],
        "reference_receipt_sha256": receipt["reference_receipt_sha256"],
        "reference_matrix_sha256": receipt["reference_matrix_sha256"],
        "reference_binding": receipt["reference_binding"],
        "reference_binding_sha256": receipt["reference_binding_sha256"],
        "completed_points": receipt["point_receipts"],
        "current_point": None,
        "confirmation_passed": receipt["confirmation_passed"],
        "receipt": "receipt.json",
        "receipt_sha256": receipt_sha256,
        "finished_utc": receipt["finished_utc"],
    }


def _validated_completed_half_step_receipt(
    output_root: Path, *, verify_resident_identity: bool = True
) -> tuple[dict[str, Any], str]:
    """Validate and independently reconstruct a completed half-step receipt."""

    receipt, raw = campaign._read_json_object(
        output_root / "receipt.json",
        label="completed half-step receipt",
        max_bytes=4 * 1024 * 1024,
    )
    expected_keys = {
        "schema",
        "state",
        "accepted_campaign_result",
        "confirmation_passed",
        "purpose",
        "authority",
        "route",
        "git_sha",
        "script_sha256",
        "diagnostic_execution_identity",
        "diagnostic_execution_identity_sha256",
        "prior_fd_root",
        "prior_fd_receipt_sha256",
        "prior_fd_terminal_sha256",
        "prior_fd_status_sha256",
        "prior_fd_git_sha",
        "prior_fd_script_sha256",
        "preflight_root",
        "preflight_receipt_sha256",
        "campaign_identity",
        "reference_root",
        "reference_receipt_sha256",
        "reference_matrix_sha256",
        "reference_git_sha",
        "reference_script_sha256",
        "reference_binding",
        "reference_binding_sha256",
        "center_receipt",
        "reused_prior_point_receipts",
        "new_point_receipts",
        "contract",
        "contract_sha256",
        "aggregate_receipt",
        "analysis",
        "finished_utc",
    }
    if set(receipt) != expected_keys:
        raise ValueError("completed half-step receipt schema is not exact")
    identity = {
        "git_sha": receipt.get("git_sha"),
        "script_sha256": receipt.get("script_sha256"),
        "diagnostic_execution_identity": receipt.get("diagnostic_execution_identity"),
        "diagnostic_execution_identity_sha256": receipt.get(
            "diagnostic_execution_identity_sha256"
        ),
    }
    execution_identity = identity["diagnostic_execution_identity"]
    reference_binding = receipt.get("reference_binding")
    if (
        receipt.get("schema") != HALF_STEP_SCHEMA
        or receipt.get("state") != "completed"
        or receipt.get("accepted_campaign_result") is not False
        or type(receipt.get("confirmation_passed")) is not bool
        or receipt.get("purpose")
        != "immutable half-step FD extension only; no production authority"
        or receipt.get("authority")
        != HALF_STEP_CONTRACT["policy"][
            (
                "passed_authority"
                if receipt.get("confirmation_passed") is True
                else "rejected_authority"
            )
        ]
        or receipt.get("route") != ROUTE
        or receipt.get("contract") != HALF_STEP_CONTRACT
        or receipt.get("contract_sha256")
        != campaign._canonical_hash(HALF_STEP_CONTRACT)
        or type(execution_identity) is not dict
        or receipt.get("diagnostic_execution_identity_sha256")
        != campaign._canonical_hash(execution_identity)
        or type(reference_binding) is not dict
        or receipt.get("reference_binding_sha256")
        != campaign._canonical_hash(reference_binding)
        or type(receipt.get("finished_utc")) is not str
        or not receipt["finished_utc"].endswith("Z")
    ):
        raise ValueError("completed half-step receipt identity is invalid")
    if (
        set(execution_identity)
        != {
            "kind",
            "loaded_module_code_sha256",
            "loaded_code_sha256",
            "loaded_state_sha256",
            "source_files",
        }
        or execution_identity.get("kind") != "python-source"
        or type(execution_identity.get("source_files")) is not dict
        or not execution_identity["source_files"]
    ):
        raise ValueError("completed half-step loaded-code identity schema is invalid")
    for field in (
        "loaded_module_code_sha256",
        "loaded_code_sha256",
        "loaded_state_sha256",
    ):
        campaign._require_sha(
            campaign._require_json_string(
                execution_identity.get(field), label=f"half-step {field}"
            ),
            length=64,
            label=f"half-step {field}",
        )
    for path, digest in execution_identity["source_files"].items():
        if type(path) is not str or not path:
            raise ValueError("completed half-step loaded-code source path is invalid")
        campaign._require_sha(
            campaign._require_json_string(
                digest, label="half-step loaded source SHA-256"
            ),
            length=64,
            label="half-step loaded source SHA-256",
        )
    for field, length in (
        ("git_sha", 40),
        ("script_sha256", 64),
        ("diagnostic_execution_identity_sha256", 64),
        ("prior_fd_receipt_sha256", 64),
        ("prior_fd_terminal_sha256", 64),
        ("prior_fd_status_sha256", 64),
        ("prior_fd_git_sha", 40),
        ("prior_fd_script_sha256", 64),
        ("preflight_receipt_sha256", 64),
        ("campaign_identity", 64),
        ("reference_receipt_sha256", 64),
        ("reference_matrix_sha256", 64),
        ("reference_git_sha", 40),
        ("reference_script_sha256", 64),
        ("reference_binding_sha256", 64),
        ("contract_sha256", 64),
    ):
        campaign._require_sha(
            campaign._require_json_string(
                receipt.get(field), label=f"half-step {field}"
            ),
            length=length,
            label=f"half-step {field}",
        )
    if verify_resident_identity:
        campaign._strict_json_equal(
            identity,
            _fd_execution_identity(),
            label="completed half-step resident execution identity",
        )

    prior_root_text = campaign._require_json_string(
        receipt.get("prior_fd_root"), label="completed half-step prior FD root"
    )
    prior_root = campaign._safe_absolute_root(Path(prior_root_text))
    if str(prior_root) != prior_root_text:
        raise ValueError("completed half-step prior FD root is not canonical")
    source = _load_half_step_source(prior_root)
    expected_status = _initial_half_step_status(identity, prior_root, source)
    prior_receipt = source["prior_receipt"]
    reference_receipt = source["reference_receipt"]
    expected_center_receipt = expected_status["center_receipt"]
    expected_provenance = {
        "prior_fd_root": str(prior_root),
        "prior_fd_receipt_sha256": source["prior_receipt_sha256"],
        "prior_fd_terminal_sha256": source["prior_terminal_sha256"],
        "prior_fd_status_sha256": PRIOR_FD_STATUS_SHA256,
        "prior_fd_git_sha": prior_receipt["git_sha"],
        "prior_fd_script_sha256": prior_receipt["script_sha256"],
        "preflight_root": str(source["preflight_root"]),
        "preflight_receipt_sha256": source["preflight_receipt_sha256"],
        "campaign_identity": prior_receipt["campaign_identity"],
        "reference_root": str(source["reference_root"]),
        "reference_receipt_sha256": REFERENCE_RECEIPT_SHA256,
        "reference_matrix_sha256": REFERENCE_MATRIX_SHA256,
        "reference_git_sha": reference_receipt["git_sha"],
        "reference_script_sha256": reference_receipt["script_sha256"],
        "reference_binding": source["reference_binding"],
        "reference_binding_sha256": expected_status["reference_binding_sha256"],
        "center_receipt": expected_center_receipt,
        "reused_prior_point_receipts": source["reused_points"],
    }
    campaign._strict_json_equal(
        {key: receipt.get(key) for key in expected_provenance},
        expected_provenance,
        label="completed half-step provenance",
    )

    plan = _half_step_plan(6)
    points = receipt.get("new_point_receipts")
    if type(points) is not list or len(points) != len(plan):
        raise ValueError("completed half-step receipt point ancestry is incomplete")
    center_density = np.asarray(source["center_density"])
    center_density_sha = hashlib.sha256(
        np.ascontiguousarray(center_density, dtype="<f8").tobytes()
    ).hexdigest()
    center_spin = _validate_point_record(
        source["center_record"],
        center_spin_square=None,
        center=True,
        enforce_scientific_gates=False,
    )
    gradients: dict[tuple[int, int], np.ndarray] = {}
    new_point_data: list[tuple[dict[str, Any], np.ndarray]] = []
    for ancestry, definition in zip(points, plan, strict=True):
        key = definition["key"]
        if (
            type(ancestry) is not dict
            or set(ancestry) != {"point_key", "receipt_sha256"}
            or ancestry.get("point_key") != key
        ):
            raise ValueError("completed half-step point ancestry is invalid")
        campaign._require_sha(
            campaign._require_json_string(
                ancestry.get("receipt_sha256"),
                label="completed half-step point receipt SHA-256",
            ),
            length=64,
            label="completed half-step point receipt SHA-256",
        )
        point, point_sha, arrays = _read_record_bundle(
            output_root,
            Path("points") / key,
            expected_receipt_keys=(
                _POINT_OBSERVATION_KEYS
                | _HALF_STEP_POINT_BINDING_KEYS
                | set(definition)
                | _DENSITY_CONTINUITY_KEYS
                | {"schema", "point_key", "kind", "artifacts"}
            ),
            expected_artifacts={
                "gradient.f64": ((6, 3), "hartree / bohr"),
                "density.f64": (center_density.shape, "electrons"),
            },
        )
        density = arrays["density.f64"]
        density_sha = hashlib.sha256(
            np.ascontiguousarray(density, dtype="<f8").tobytes()
        ).hexdigest()
        expected_coords = np.asarray(
            source["transition_state"].coords, dtype=float
        ).copy()
        expected_coords.reshape(-1)[definition["coordinate_index"]] += (
            definition["displacement_bohr"] * BOHR_TO_ANGSTROM
        )
        if (
            point_sha != ancestry["receipt_sha256"]
            or point.get("schema") != HALF_STEP_SCHEMA
            or point.get("point_key") != key
            or point.get("kind") != "new-half-step-displacement"
            or any(point.get(field) != value for field, value in definition.items())
            or any(
                point.get(field) != value
                for field, value in _half_step_point_bindings(expected_status).items()
            )
            or point.get("density_initial_guess_sha256") != center_density_sha
            or point.get("density_final_sha256") != density_sha
            or point.get("geometry_fingerprint")
            != frequency_geometry_fingerprint(
                replace(source["transition_state"], coords=expected_coords)
            )
        ):
            raise ValueError(f"completed half-step point content is invalid: {key}")
        _validate_point_record(
            point,
            center_spin_square=center_spin,
            center=False,
            enforce_scientific_gates=False,
        )
        campaign._strict_json_equal(
            {name: point.get(name) for name in _DENSITY_CONTINUITY_KEYS},
            _density_continuity_metrics(density, center_density),
            label=f"completed half-step density continuity {key}",
        )
        gradients[(definition["coordinate_index"], definition["step_multiplier"])] = (
            arrays["gradient.f64"].reshape(-1)
        )
        new_point_data.append((point, density))

    h_inner, h_outer, richardson = _half_step_matrices(
        gradients, source["reused_gradients"], 18
    )
    reconstructed = {
        "H_0.005-raw.f64": h_inner,
        "H_0.01-raw.f64": h_outer,
        "richardson-raw.f64": richardson,
        "H_0.005-symmetric.f64": 0.5 * (h_inner + h_inner.T),
        "H_0.01-symmetric.f64": 0.5 * (h_outer + h_outer.T),
        "richardson-symmetric.f64": 0.5 * (richardson + richardson.T),
    }
    matrix_contract = {
        filename: ((18, 18), "hartree / bohr^2") for filename in reconstructed
    }
    aggregate = receipt.get("aggregate_receipt")
    if (
        type(aggregate) is not dict
        or set(aggregate) != {"path", "sha256"}
        or aggregate.get("path") != "matrices/receipt.json"
    ):
        raise ValueError("completed half-step aggregate ancestry is invalid")
    campaign._require_sha(
        campaign._require_json_string(
            aggregate.get("sha256"), label="completed half-step aggregate SHA-256"
        ),
        length=64,
        label="completed half-step aggregate SHA-256",
    )
    aggregate_payload, aggregate_sha, matrix_arrays = _read_record_bundle(
        output_root,
        Path("matrices"),
        expected_receipt_keys={
            "schema",
            "kind",
            "confirmation_passed",
            *_HALF_STEP_POINT_BINDING_KEYS,
            "new_point_receipts",
            "reused_prior_point_receipts",
            "analysis",
            "artifacts",
        },
        expected_artifacts=matrix_contract,
    )
    expected_aggregate = {
        "schema": HALF_STEP_SCHEMA,
        "kind": "immutable-half-step-finite-difference-analysis",
        "confirmation_passed": receipt["confirmation_passed"],
        **_half_step_point_bindings(expected_status),
        "new_point_receipts": points,
        "reused_prior_point_receipts": source["reused_points"],
        "analysis": receipt["analysis"],
        "artifacts": aggregate_payload["artifacts"],
    }
    if aggregate_sha != aggregate["sha256"]:
        raise ValueError("completed half-step aggregate receipt hash mismatch")
    campaign._strict_json_equal(
        aggregate_payload,
        expected_aggregate,
        label="completed half-step aggregate receipt",
    )
    for filename, expected_matrix in reconstructed.items():
        if not np.array_equal(matrix_arrays[filename], expected_matrix):
            raise ValueError(
                f"completed half-step matrix diverges from reconstruction: {filename}"
            )
    route_vector = campaign._mapped_route_vector(
        source["transition_state"],
        source["masses"],
        source["reactant"].coords,
        source["product"].coords,
    )
    expected_analysis = _complete_half_step_analysis(
        source["transition_state"],
        source["masses"],
        route_vector,
        h_inner,
        h_outer,
        richardson,
        source["reference_matrix"],
        point_gate_failures=_half_step_point_gates(source, new_point_data),
    )
    if (
        type(expected_analysis) is not dict
        or type(expected_analysis.get("confirmation_passed")) is not bool
        or expected_analysis["confirmation_passed"]
        is not receipt["confirmation_passed"]
    ):
        raise ValueError("completed half-step scientific verdict is not reproducible")
    campaign._strict_json_equal(
        receipt.get("analysis"),
        expected_analysis,
        label="completed half-step independently recomputed analysis",
    )
    return receipt, hashlib.sha256(raw).hexdigest()


def _completed_half_step_status(
    receipt: dict[str, Any], receipt_sha256: str
) -> dict[str, Any]:
    return {
        "schema": HALF_STEP_SCHEMA,
        "state": "completed",
        "route": ROUTE,
        "git_sha": receipt["git_sha"],
        "script_sha256": receipt["script_sha256"],
        "diagnostic_execution_identity": receipt["diagnostic_execution_identity"],
        "diagnostic_execution_identity_sha256": receipt[
            "diagnostic_execution_identity_sha256"
        ],
        "prior_fd_root": receipt["prior_fd_root"],
        "prior_fd_receipt_sha256": receipt["prior_fd_receipt_sha256"],
        "prior_fd_terminal_sha256": receipt["prior_fd_terminal_sha256"],
        "prior_fd_status_sha256": receipt["prior_fd_status_sha256"],
        "preflight_root": receipt["preflight_root"],
        "preflight_receipt_sha256": receipt["preflight_receipt_sha256"],
        "reference_root": receipt["reference_root"],
        "reference_receipt_sha256": receipt["reference_receipt_sha256"],
        "reference_matrix_sha256": receipt["reference_matrix_sha256"],
        "reference_binding": receipt["reference_binding"],
        "reference_binding_sha256": receipt["reference_binding_sha256"],
        "center_receipt": receipt["center_receipt"],
        "reused_prior_point_receipts": receipt["reused_prior_point_receipts"],
        "contract": receipt["contract"],
        "contract_sha256": receipt["contract_sha256"],
        "completed_points": receipt["new_point_receipts"],
        "current_point": None,
        "accepted_campaign_result": False,
        "confirmation_passed": receipt["confirmation_passed"],
        "receipt": "receipt.json",
        "receipt_sha256": receipt_sha256,
        "finished_utc": receipt["finished_utc"],
    }


def _validate_failed_half_step_terminal(
    output_root: Path,
    terminal: dict[str, Any],
    status: dict[str, Any],
) -> None:
    """Bind a failed half-step terminal to its identity and durable progress."""

    identity = _half_step_failure_execution_identity()
    expected_identity = {
        "git_sha": identity["git_sha"],
        "script_sha256": identity["script_sha256"],
        "diagnostic_execution_identity_sha256": identity[
            "diagnostic_execution_identity_sha256"
        ],
        "contract_sha256": campaign._canonical_hash(HALF_STEP_CONTRACT),
        "preflight_receipt_sha256": PRIOR_FD_PREFLIGHT_SHA256,
        "prior_fd_receipt_sha256": PRIOR_FD_RECEIPT_SHA256,
        "prior_fd_terminal_sha256": PRIOR_FD_TERMINAL_SHA256,
        "prior_fd_status_sha256": PRIOR_FD_STATUS_SHA256,
        "reference_receipt_sha256": REFERENCE_RECEIPT_SHA256,
        "reference_matrix_sha256": REFERENCE_MATRIX_SHA256,
    }
    if any(terminal.get(key) != value for key, value in expected_identity.items()):
        raise ValueError("failed half-step terminal identity hash mismatch")

    required_status_fields = {
        "schema",
        "state",
        "accepted_campaign_result",
        "confirmation_passed",
        "git_sha",
        "script_sha256",
        "diagnostic_execution_identity_sha256",
        "contract_sha256",
        "preflight_receipt_sha256",
        "reference_receipt_sha256",
        "reference_matrix_sha256",
        "prior_fd_receipt_sha256",
        "prior_fd_terminal_sha256",
        "prior_fd_status_sha256",
        "completed_points",
        "current_point",
        "error",
        "finished_utc",
    }
    if (
        not required_status_fields <= set(status)
        or status.get("schema") != HALF_STEP_SCHEMA
        or status.get("state") != "failed"
        or status.get("accepted_campaign_result") is not False
        or status.get("confirmation_passed") is not False
        or terminal.get("detail") != status.get("error")
        or terminal.get("finished_utc") != status.get("finished_utc")
    ):
        raise ValueError("failed half-step terminal status contract is invalid")
    binding_fields = {
        "git_sha",
        "script_sha256",
        "diagnostic_execution_identity_sha256",
        "contract_sha256",
        "preflight_receipt_sha256",
        "reference_receipt_sha256",
        "reference_matrix_sha256",
        "prior_fd_receipt_sha256",
        "prior_fd_terminal_sha256",
        "prior_fd_status_sha256",
    }
    if any(terminal.get(key) != status.get(key) for key in binding_fields):
        raise ValueError("failed half-step terminal/status identity hash mismatch")

    completed = status.get("completed_points")
    expected_keys = [item["key"] for item in _half_step_plan(6)]
    if (
        type(completed) is not list
        or len(completed) != terminal["completed_point_count"]
    ):
        raise ValueError("failed half-step terminal point-count mismatch")
    for ancestry, expected_key in zip(completed, expected_keys, strict=False):
        if (
            type(ancestry) is not dict
            or set(ancestry) != {"point_key", "receipt_sha256"}
            or ancestry.get("point_key") != expected_key
        ):
            raise ValueError("failed half-step terminal progress is invalid")
        digest = campaign._require_sha(
            campaign._require_json_string(
                ancestry.get("receipt_sha256"),
                label="failed half-step point receipt SHA-256",
            ),
            length=64,
            label="failed half-step point receipt SHA-256",
        )
        raw = campaign._read_bounded_regular_snapshot(
            output_root / "points" / expected_key / "receipt.json",
            label="failed half-step durable point receipt",
            maximum_bytes=1024 * 1024,
        )
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError("failed half-step terminal progress receipt hash mismatch")
    if terminal.get("current_point") != status.get("current_point"):
        raise ValueError("failed half-step terminal current-point mismatch")


def _prepare_failed_half_step_status(status: dict[str, Any]) -> dict[str, Any]:
    """Fill the immutable half-step identity before persisting failure state."""

    identity = _half_step_failure_execution_identity()
    status.update(
        {
            "schema": HALF_STEP_SCHEMA,
            "route": ROUTE,
            **identity,
            "prior_fd_root": status.get("prior_fd_root", PRIOR_FD_ROOT),
            "prior_fd_receipt_sha256": PRIOR_FD_RECEIPT_SHA256,
            "prior_fd_terminal_sha256": PRIOR_FD_TERMINAL_SHA256,
            "prior_fd_status_sha256": PRIOR_FD_STATUS_SHA256,
            "preflight_root": status.get("preflight_root", PRIOR_FD_PREFLIGHT_ROOT),
            "preflight_receipt_sha256": PRIOR_FD_PREFLIGHT_SHA256,
            "reference_root": status.get("reference_root", PRIOR_FD_REFERENCE_ROOT),
            "reference_receipt_sha256": REFERENCE_RECEIPT_SHA256,
            "reference_matrix_sha256": REFERENCE_MATRIX_SHA256,
            "contract": HALF_STEP_CONTRACT,
            "contract_sha256": campaign._canonical_hash(HALF_STEP_CONTRACT),
            "completed_points": status.get("completed_points", []),
            "current_point": status.get("current_point"),
            "accepted_campaign_result": False,
            "confirmation_passed": False,
        }
    )
    return status


def _reconcile_status_from_terminal(
    claim: _OutputRootClaim,
    status_path: Path,
    status: dict[str, Any],
    terminal: dict[str, Any],
) -> None:
    status.update(
        {
            "schema": terminal["schema"],
            "state": terminal["state"],
            "route": ROUTE,
            "error": terminal["detail"] if terminal["state"] == "failed" else None,
            "finished_utc": terminal["finished_utc"],
        }
    )
    if terminal["schema"] in {FINITE_DIFFERENCE_SCHEMA, HALF_STEP_SCHEMA}:
        status["accepted_campaign_result"] = False
        status["confirmation_passed"] = terminal["confirmation_passed"]
    if terminal["state"] == "completed":
        status["receipt"] = terminal["receipt"]
        status["receipt_sha256"] = terminal["receipt_sha256"]
    _claimed_write(claim, status_path, _json_bytes(status))


def finalize_if_running(
    output_root: Path,
    *,
    finite_difference: bool = False,
    half_step: bool = False,
) -> dict[str, Any]:
    """Publish or reconcile one no-clobber terminal under the run's claim."""

    if finite_difference and half_step:
        raise ValueError("dead-man finalization mode must be unambiguous")
    output = campaign._safe_absolute_root(output_root)
    with _exclusive_output_claim(output) as claim:
        terminal_path = output / TERMINAL
        status_path = output / "status.json"
        status: dict[str, Any] = {}
        if status_path.exists() or status_path.is_symlink():
            try:
                status, _ = campaign._read_json_object(
                    status_path,
                    label="diagnostic status",
                    max_bytes=1024 * 1024,
                )
            except ValueError:
                status = {}
        if terminal_path.exists() or terminal_path.is_symlink():
            terminal = _validated_terminal(terminal_path)
            if half_step and terminal["schema"] != HALF_STEP_SCHEMA:
                raise ValueError(
                    "dead-man half-step terminal schema does not match requested mode"
                )
            if terminal["schema"] == HALF_STEP_SCHEMA and terminal["state"] == "failed":
                _validate_failed_half_step_terminal(output, terminal, status)
            if (
                terminal["state"] == "completed"
                and terminal["schema"] == FINITE_DIFFERENCE_SCHEMA
            ):
                completed_receipt, receipt_sha = _validated_completed_fd_receipt(output)
                completed_status = _completed_fd_status(completed_receipt, receipt_sha)
                expected_terminal = _terminal_payload(
                    completed_status, state="completed", detail=None
                )
                expected_terminal["receipt"] = "receipt.json"
                expected_terminal["receipt_sha256"] = receipt_sha
                campaign._strict_json_equal(
                    terminal,
                    expected_terminal,
                    label="authoritative completed FD terminal",
                )
                _claimed_write(claim, status_path, _json_bytes(completed_status))
                return terminal
            if (
                terminal["state"] == "completed"
                and terminal["schema"] == HALF_STEP_SCHEMA
            ):
                completed_receipt, receipt_sha = _validated_completed_half_step_receipt(
                    output
                )
                completed_status = _completed_half_step_status(
                    completed_receipt, receipt_sha
                )
                expected_terminal = _terminal_payload(
                    completed_status, state="completed", detail=None
                )
                expected_terminal["receipt"] = "receipt.json"
                expected_terminal["receipt_sha256"] = receipt_sha
                campaign._strict_json_equal(
                    terminal,
                    expected_terminal,
                    label="authoritative completed half-step terminal",
                )
                _claimed_write(claim, status_path, _json_bytes(completed_status))
                return terminal
            if terminal["state"] == "completed":
                receipt_raw = campaign._read_bounded_regular_snapshot(
                    output / terminal["receipt"],
                    label="completed diagnostic receipt",
                    maximum_bytes=4 * 1024 * 1024,
                )
                receipt_sha = hashlib.sha256(receipt_raw).hexdigest()
                if terminal["receipt_sha256"] != receipt_sha:
                    raise ValueError("terminal completed receipt hash mismatch")
            _reconcile_status_from_terminal(claim, status_path, status, terminal)
            return terminal

        completed_receipt: dict[str, Any] | None = None
        completed_receipt_sha: str | None = None
        completed_schema: str | None = None
        if (
            (output / "receipt.json").exists() or (output / "receipt.json").is_symlink()
        ) and (
            finite_difference
            or half_step
            or status.get("schema") in {FINITE_DIFFERENCE_SCHEMA, HALF_STEP_SCHEMA}
        ):
            requested_schema = (
                FINITE_DIFFERENCE_SCHEMA
                if finite_difference
                else HALF_STEP_SCHEMA
                if half_step
                else status.get("schema")
            )
            try:
                if requested_schema == HALF_STEP_SCHEMA:
                    completed_receipt, completed_receipt_sha = (
                        _validated_completed_half_step_receipt(output)
                    )
                else:
                    completed_receipt, completed_receipt_sha = (
                        _validated_completed_fd_receipt(output)
                    )
                completed_schema = requested_schema
            except (ValueError, RuntimeError) as exc:
                status.pop("receipt", None)
                status.pop("receipt_sha256", None)
                failure_schema = requested_schema
                if failure_schema == HALF_STEP_SCHEMA:
                    _prepare_failed_half_step_status(status)
                status.update(
                    {
                        "schema": failure_schema,
                        "state": "failed",
                        "route": ROUTE,
                        "accepted_campaign_result": False,
                        "confirmation_passed": False,
                        "current_point": None,
                        "finished_utc": time.strftime(
                            "%Y-%m-%dT%H:%M:%SZ", time.gmtime()
                        ),
                    }
                )
                detail = f"InvalidCompletedReceipt: {type(exc).__name__}: {exc}"
                terminal = _terminal_payload(status, state="failed", detail=detail)
                _terminal_write_noreplace(claim, terminal_path, terminal)
                _reconcile_status_from_terminal(claim, status_path, status, terminal)
                return terminal
        if completed_receipt is not None:
            if completed_receipt_sha is None:
                raise RuntimeError("completed receipt hash was not computed")
            if completed_schema == HALF_STEP_SCHEMA:
                status = _completed_half_step_status(
                    completed_receipt, completed_receipt_sha
                )
            else:
                status = _completed_fd_status(completed_receipt, completed_receipt_sha)
            terminal = _terminal_payload(status, state="completed", detail=None)
            terminal["receipt"] = "receipt.json"
            terminal["receipt_sha256"] = completed_receipt_sha
        else:
            if finite_difference and status.get("schema") != FINITE_DIFFERENCE_SCHEMA:
                status = {
                    "schema": FINITE_DIFFERENCE_SCHEMA,
                    "completed_points": [],
                    "current_point": None,
                    "confirmation_passed": False,
                }
            elif half_step or status.get("schema") == HALF_STEP_SCHEMA:
                status = _prepare_failed_half_step_status(status)
            terminal = _terminal_payload(
                status,
                state="failed",
                detail=(
                    "DiagnosticInterrupted: bounded systemd unit ended without a "
                    "terminal receipt"
                ),
            )
        _terminal_write_noreplace(claim, terminal_path, terminal)
        if completed_receipt is not None:
            _claimed_write(claim, status_path, _json_bytes(status))
        else:
            _reconcile_status_from_terminal(claim, status_path, status, terminal)
        return terminal


def _matrix_artifact(
    claim: _OutputRootClaim,
    root: Path,
    case: str,
    label: str,
    matrix: np.ndarray,
) -> dict[str, Any]:
    array = np.ascontiguousarray(matrix, dtype="<f8")
    raw = array.tobytes()
    relative = Path(case) / f"{label}.f64"
    destination = root / relative
    claim.verify()
    destination.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    claim.verify()
    _claimed_write(claim, destination, raw)
    return {
        "path": str(relative),
        "dtype": "little-endian float64",
        "shape": list(array.shape),
        "units": "hartree / bohr^2",
        "sha256": hashlib.sha256(raw).hexdigest(),
    }


def _canonical_native(native: Any, atom_count: int) -> np.ndarray:
    array = np.asarray(native.get() if hasattr(native, "get") else native, dtype=float)
    expected = (atom_count, atom_count, 3, 3)
    if array.shape != expected or not np.all(np.isfinite(array)):
        raise ValueError(
            f"diagnostic native Hessian must be finite {expected}, got {array.shape}"
        )
    return array.transpose(0, 2, 1, 3).reshape(3 * atom_count, 3 * atom_count)


def _metric_payload(matrix: np.ndarray) -> dict[str, Any]:
    payload = asdict(native_hessian.cartesian_hessian_symmetry_metrics(matrix))
    for key, value in tuple(payload.items()):
        if isinstance(value, float) and not math.isfinite(value):
            payload[key] = "infinity"
    return payload


def _diagnostic_execution_identity() -> dict[str, Any]:
    """Prove the loaded module code and literal state match its source bytes."""

    source = Path(__file__)
    if source.is_symlink():
        raise RuntimeError("diagnostic executable source cannot be a symlink")
    origin = source.resolve(strict=True)
    raw = campaign._read_bounded_regular_snapshot(
        origin,
        label="Hessian diagnostic executable module",
        maximum_bytes=2 * 1024 * 1024,
    )
    source_code = compile(
        raw.decode("utf-8"),
        str(origin),
        "exec",
        dont_inherit=True,
        optimize=sys.flags.optimize,
    )
    source_digests = campaign._nested_code_digests(source_code)

    def source_literal(node: ast.expr, known: dict[str, Any]) -> Any:
        if isinstance(node, ast.Constant):
            return node.value
        if isinstance(node, ast.Name) and node.id in known:
            return known[node.id]
        if isinstance(node, ast.Dict):
            if any(key is None for key in node.keys):
                raise ValueError(
                    "scientific literal dictionary unpacking is unsupported"
                )
            return {
                source_literal(key, known): source_literal(value, known)
                for key, value in zip(node.keys, node.values, strict=True)
                if key is not None
            }
        if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
            values = [source_literal(value, known) for value in node.elts]
            return {ast.List: list, ast.Tuple: tuple, ast.Set: set}[type(node)](values)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
            value = source_literal(node.operand, known)
            if type(value) not in {int, float, complex}:
                raise ValueError("non-numeric unary scientific literal")
            return value if isinstance(node.op, ast.UAdd) else -value
        raise ValueError("unsupported scientific literal expression")

    tree = ast.parse(raw.decode("utf-8"), filename=str(origin))
    known_literals: dict[str, Any] = {}
    for node in tree.body:
        if not isinstance(node, (ast.Assign, ast.AnnAssign)):
            continue
        targets = node.targets if isinstance(node, ast.Assign) else [node.target]
        if (
            len(targets) != 1
            or not isinstance(targets[0], ast.Name)
            or node.value is None
        ):
            continue
        try:
            known_literals[targets[0].id] = source_literal(node.value, known_literals)
        except ValueError:
            continue
    scientific_names = {
        "SCHEMA",
        "FINITE_DIFFERENCE_SCHEMA",
        "HALF_STEP_SCHEMA",
        "TERMINAL",
        "ROUTE",
        "PRIOR_FD_ROOT",
        "PRIOR_FD_RECEIPT_SHA256",
        "PRIOR_FD_TERMINAL_SHA256",
        "PRIOR_FD_STATUS_SHA256",
        "PRIOR_FD_GIT_SHA",
        "PRIOR_FD_SCRIPT_SHA256",
        "PRIOR_FD_PREFLIGHT_ROOT",
        "PRIOR_FD_PREFLIGHT_SHA256",
        "PRIOR_FD_REFERENCE_ROOT",
        "REFERENCE_RECEIPT_SHA256",
        "REFERENCE_MATRIX_SHA256",
        "FINITE_DIFFERENCE_CONTRACT",
        "HALF_STEP_CONTRACT",
        "CASES",
    }
    if scientific_names - known_literals.keys():
        raise RuntimeError("diagnostic scientific source literals are incomplete")
    scientific_state = {}
    for name in sorted(scientific_names):
        expected = campaign._literal_state_payload(known_literals[name])
        try:
            observed = campaign._literal_state_payload(
                vars(sys.modules[__name__])[name]
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise RuntimeError(
                f"loaded diagnostic scientific state is invalid: {__name__}.{name}"
            ) from exc
        if observed != expected:
            raise RuntimeError(
                "loaded diagnostic scientific state disagrees with source: "
                f"{__name__}.{name}"
            )
        scientific_state[name] = observed
    scientific_state["BOHR_TO_ANGSTROM"] = campaign._literal_state_payload(
        BOHR_TO_ANGSTROM
    )
    scientific_state["HARTREE_TO_EV"] = campaign._literal_state_payload(HARTREE_TO_EV)
    scientific_state["campaign.DFT_SETTINGS"] = campaign._literal_state_payload(
        campaign.DFT_SETTINGS
    )

    expected_slots: dict[str, str] = {}
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            expected_slots[node.name] = node.name
        elif isinstance(node, ast.ClassDef):
            for member in node.body:
                if isinstance(member, (ast.FunctionDef, ast.AsyncFunctionDef)):
                    expected_slots[f"{node.name}.{member.name}"] = (
                        f"{node.name}.{member.name}"
                    )

    def require_local(function: Any, *, label: str, qualname: str) -> None:
        if not campaign.inspect.isfunction(function):
            raise RuntimeError(
                f"loaded diagnostic executable slot is missing: {__name__}.{label}"
            )
        function = campaign.inspect.unwrap(function)
        if (
            not campaign.inspect.isfunction(function)
            or getattr(function, "__module__", None) != __name__
        ):
            raise RuntimeError(
                "loaded diagnostic executable slot has a foreign owner: "
                f"{__name__}.{label}"
            )
        code = function.__code__
        try:
            code_origin = Path(code.co_filename).resolve(strict=True)
        except (OSError, RuntimeError) as exc:
            raise RuntimeError(
                f"loaded diagnostic code has no source origin: {label}"
            ) from exc
        if (
            code_origin != origin
            or code.co_qualname != qualname
            or campaign._code_digest(code) not in source_digests.get(qualname, set())
        ):
            raise RuntimeError(
                f"loaded Python code disagrees with source: {__name__}.{label}"
            )

    module_values = vars(sys.modules[__name__])
    for label, qualname in expected_slots.items():
        if "." not in label:
            value = module_values.get(label)
        else:
            class_name, member_name = label.split(".", 1)
            owner = module_values.get(class_name)
            if (
                not campaign.inspect.isclass(owner)
                or getattr(owner, "__module__", None) != __name__
            ):
                raise RuntimeError(
                    "loaded diagnostic executable class is missing: "
                    f"{__name__}.{class_name}"
                )
            value = vars(owner).get(member_name)
            if isinstance(value, (staticmethod, classmethod)):
                value = value.__func__
            elif isinstance(value, property):
                value = value.fget
        require_local(value, label=label, qualname=qualname)
    identity = campaign._python_source_execution_identity(
        sys.modules[__name__],
        __name__,
        origin,
        raw,
        trusted_source_roots=(
            origin.parents[1],
            Path(sys.prefix).resolve(),
            Path(sys.base_prefix).resolve(),
        ),
    )
    identity["loaded_state_sha256"] = campaign._canonical_hash(
        {
            "generic_loaded_literal_state_sha256": identity["loaded_state_sha256"],
            "scientific_runtime_state": scientific_state,
        }
    )
    return identity


def _git_identity() -> dict[str, Any]:
    repository = Path(__file__).resolve().parents[2]
    status = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=all"],
        cwd=repository,
        check=True,
        capture_output=True,
        text=True,
    )
    if status.stdout:
        raise RuntimeError("diagnostic requires a clean Git worktree")
    sha = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repository,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    script_raw = campaign._read_bounded_regular_snapshot(
        Path(__file__),
        label="Hessian diagnostic committed source",
        maximum_bytes=2 * 1024 * 1024,
    )
    return {
        "git_sha": sha,
        "script_sha256": hashlib.sha256(script_raw).hexdigest(),
    }


def _fd_execution_identity() -> dict[str, Any]:
    """Bind loaded FD code to the exact clean committed diagnostic source."""

    identity = _git_identity()
    repository = Path(__file__).resolve().parents[2]
    source = Path(__file__).resolve(strict=True)
    relative_source = source.relative_to(repository)
    committed_source = subprocess.run(
        ["git", "show", f"{identity['git_sha']}:{relative_source.as_posix()}"],
        cwd=repository,
        check=True,
        capture_output=True,
    ).stdout
    committed_sha = hashlib.sha256(committed_source).hexdigest()
    if committed_sha != identity["script_sha256"]:
        raise RuntimeError("loaded diagnostic source is not the committed Git source")
    execution_identity = _diagnostic_execution_identity()
    if execution_identity.get("source_files", {}).get(str(source)) != committed_sha:
        raise RuntimeError(
            "loaded diagnostic identity does not bind its committed source"
        )
    if _git_identity() != identity:
        raise RuntimeError(
            "diagnostic Git identity changed during loaded-code attestation"
        )
    return {
        **identity,
        "diagnostic_execution_identity": execution_identity,
        "diagnostic_execution_identity_sha256": campaign._canonical_hash(
            execution_identity
        ),
    }


def _half_step_failure_execution_identity() -> dict[str, Any]:
    """Attest the loaded failure path even when clean-Git preflight rejected it."""

    try:
        return _fd_execution_identity()
    except RuntimeError:
        repository = Path(__file__).resolve().parents[2]
        git_sha = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repository,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        source = Path(__file__).resolve(strict=True)
        script_sha256 = hashlib.sha256(
            campaign._read_bounded_regular_snapshot(
                source,
                label="failed half-step loaded source",
                maximum_bytes=2 * 1024 * 1024,
            )
        ).hexdigest()
        execution_identity = {
            "kind": "python-source-failure",
            "loaded_module_code_sha256": script_sha256,
            "loaded_code_sha256": script_sha256,
            "loaded_state_sha256": campaign._canonical_hash(
                {
                    "schema": HALF_STEP_SCHEMA,
                    "route": ROUTE,
                    "prior_fd_root": PRIOR_FD_ROOT,
                    "contract_sha256": campaign._canonical_hash(HALF_STEP_CONTRACT),
                }
            ),
            "source_files": {str(source): script_sha256},
        }
        return {
            "git_sha": git_sha,
            "script_sha256": script_sha256,
            "diagnostic_execution_identity": execution_identity,
            "diagnostic_execution_identity_sha256": campaign._canonical_hash(
                execution_identity
            ),
        }


def _evaluate_case(
    claim: _OutputRootClaim,
    root: Path,
    definition: dict[str, Any],
    transition_state: Any,
    mapped_reactant: np.ndarray,
    mapped_product: np.ndarray,
    masses: np.ndarray,
    base_settings: Any,
) -> tuple[dict[str, Any], np.ndarray, Any]:
    started = time.monotonic()
    settings = replace(
        base_settings,
        use_gpu=False,
        grid_level=int(definition["grid_level"]),
    )
    mf = _make_scf(build_mol(transition_state, settings), settings)
    mf.conv_tol = float(definition["scf_tolerance"])
    mf.conv_tol_cpscf = float(definition["cpscf_tolerance"])
    electronic = float(mf.kernel())
    if not mf.converged:
        raise RuntimeError(f"SCF did not converge for {definition['name']}")

    hessian_method = _hessian_method(mf, settings)
    hessian_method.max_cycle = int(definition["hessian_max_cycle"])
    electronic_component = hessian_method.hess_elec(
        mf.mo_energy,
        mf.mo_coeff,
        mf.mo_occ,
        atmlst=hessian_method.atmlst,
    )
    nuclear_component = hessian_method.hess_nuc(mf.mol, atmlst=hessian_method.atmlst)
    dispersion_component = (
        hessian_method.get_dispersion()
        if mf.do_disp()
        else np.zeros_like(electronic_component)
    )
    native_total = electronic_component + nuclear_component + dispersion_component
    atom_count = len(transition_state.symbols)
    components = {
        "electronic": _canonical_native(electronic_component, atom_count),
        "nuclear": _canonical_native(nuclear_component, atom_count),
        "dispersion": _canonical_native(dispersion_component, atom_count),
        "total": _canonical_native(native_total, atom_count),
    }
    symmetric = 0.5 * (components["total"] + components["total"].T)

    gradient_value = _gradient_method(mf, settings).kernel()
    if hasattr(gradient_value, "get"):
        gradient_value = gradient_value.get()
    gradient = np.asarray(gradient_value, dtype=float)
    physical_fmax = float(np.max(np.linalg.norm(gradient, axis=1)))
    physical_fmax *= HARTREE_TO_EV / BOHR_TO_ANGSTROM

    diagnostic_fingerprint = campaign._canonical_hash(
        {
            "settings": campaign.DFT_SETTINGS,
            "actual_grid_level": int(definition["grid_level"]),
            "scf_tolerance": float(definition["scf_tolerance"]),
            "cpscf_tolerance": float(definition["cpscf_tolerance"]),
            "hessian_max_cycle": int(definition["hessian_max_cycle"]),
            "backend": "pyscf",
        }
    )
    result = NativeHessianResult(
        electronic_hartree=electronic,
        gradient_hartree_per_bohr=gradient,
        physical_fmax_ev_per_angstrom=physical_fmax,
        cartesian_hessian_hartree_per_bohr2=symmetric,
        requested_backend="pyscf",
        actual_backend="pyscf",
        gpu_fallback_used=False,
        geometry_fingerprint=frequency_geometry_fingerprint(transition_state),
        settings_fingerprint=diagnostic_fingerprint,
    )
    modes = project_vibrational_hessian(transition_state.coords, masses, symmetric)
    wavenumbers = hessian_eigenvalues_to_wavenumbers_cm(modes.eigenvalues)
    route_vector = campaign._mapped_route_vector(
        transition_state,
        masses,
        mapped_reactant,
        mapped_product,
    )
    try:
        gate = campaign.validate_transition_state_gate(
            transition_state,
            masses,
            result,
            expected_settings_fingerprint=diagnostic_fingerprint,
            reaction_vector_mass_scaled=route_vector,
            reaction_vector_source="diagnostic:trusted-mapped-reactant-product-displacement",
            mapped_reactant_coordinates_angstrom=mapped_reactant,
            mapped_product_coordinates_angstrom=mapped_product,
        )
        gate_verdict: dict[str, Any] = {"accepted": True, "evidence": gate}
    except ValueError as exc:
        gate_verdict = {"accepted": False, "detail": f"{type(exc).__name__}: {exc}"}

    component_records = {}
    for label, matrix in components.items():
        component_records[label] = {
            "artifact": _matrix_artifact(
                claim, root, definition["name"], label, matrix
            ),
            "symmetry": _metric_payload(matrix),
        }
    symmetric_artifact = _matrix_artifact(
        claim, root, definition["name"], "total-symmetric", symmetric
    )
    grid_points = int(np.asarray(mf.grids.coords).shape[0])
    spin_square = None
    if hasattr(mf, "spin_square"):
        observed = mf.spin_square()
        spin_square = [float(value) for value in observed]
    record = {
        **definition,
        "elapsed_seconds": time.monotonic() - started,
        "electronic_hartree": electronic,
        "scf_converged": bool(mf.converged),
        "grid_point_count": grid_points,
        "spin_square": spin_square,
        "physical_fmax_ev_per_angstrom": physical_fmax,
        "components": component_records,
        "symmetric_total_artifact": symmetric_artifact,
        "projected_eigenvalues_hartree_per_bohr2_amu": modes.eigenvalues.tolist(),
        "signed_wavenumbers_cm": wavenumbers.tolist(),
        "negative_mode_count_below_1e-8": int(
            np.count_nonzero(modes.eigenvalues < -1.0e-8)
        ),
        "near_zero_mode_count_at_1e-8": int(
            np.count_nonzero(np.abs(modes.eigenvalues) <= 1.0e-8)
        ),
        "post_symmetrization_gate": gate_verdict,
    }
    return record, symmetric, modes


def _finite_number(value: Any, *, label: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} must be a finite number")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be a finite number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} must be a finite number")
    return result


def _strict_less(value: Any, limit: float, *, label: str) -> float:
    result = _finite_number(value, label=label)
    if not result < limit:
        raise ScientificRejection(
            f"{label} must be strictly less than {limit}; got {result}"
        )
    return result


def _strict_greater(value: Any, limit: float, *, label: str) -> float:
    result = _finite_number(value, label=label)
    if not result > limit:
        raise ScientificRejection(
            f"{label} must be strictly greater than {limit}; got {result}"
        )
    return result


def _displacement_plan(atom_count: int) -> tuple[dict[str, Any], ...]:
    if type(atom_count) is not int or atom_count < 1:
        raise ValueError("finite-difference atom count must be a positive integer")
    plan = []
    for coordinate in range(3 * atom_count):
        atom, axis = divmod(coordinate, 3)
        for multiplier in (-2, -1, 1, 2):
            token = f"m{abs(multiplier)}" if multiplier < 0 else f"p{multiplier}"
            plan.append(
                {
                    "key": f"coordinate-{coordinate:03d}-{token}",
                    "coordinate_index": coordinate,
                    "atom_index": atom,
                    "axis_index": axis,
                    "step_multiplier": multiplier,
                    "displacement_bohr": 0.01 * multiplier,
                }
            )
    return tuple(plan)


def _finite_difference_matrices(
    gradients: dict[tuple[int, int], np.ndarray], dimension: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return H_h, H_2h, and fourth-order Richardson in gradient-row order."""

    expected = {
        (coordinate, step) for coordinate in range(dimension) for step in (-2, -1, 1, 2)
    }
    if set(gradients) != expected:
        raise ValueError(
            "finite-difference gradient inventory is incomplete or unexpected"
        )
    h_h = np.empty((dimension, dimension), dtype=float)
    h_2h = np.empty_like(h_h)
    for coordinate in range(dimension):
        values = {}
        for multiplier in (-2, -1, 1, 2):
            gradient = np.asarray(gradients[(coordinate, multiplier)], dtype=float)
            if gradient.shape != (dimension,) or not np.all(np.isfinite(gradient)):
                raise ValueError(
                    "finite-difference gradient has invalid shape or values"
                )
            values[multiplier] = gradient
        h_h[:, coordinate] = (values[1] - values[-1]) / 0.02
        h_2h[:, coordinate] = (values[2] - values[-2]) / 0.04
    richardson = (4.0 * h_h - h_2h) / 3.0
    return h_h, h_2h, richardson


def _half_step_plan(atom_count: int) -> tuple[dict[str, Any], ...]:
    """Return only the missing +/-0.005 Bohr displacement pairs."""

    if type(atom_count) is not int or atom_count < 1:
        raise ValueError("half-step atom count must be a positive integer")
    plan = []
    for coordinate in range(3 * atom_count):
        atom, axis = divmod(coordinate, 3)
        for multiplier in (-1, 1):
            token = "m1" if multiplier < 0 else "p1"
            plan.append(
                {
                    "key": f"coordinate-{coordinate:03d}-{token}",
                    "coordinate_index": coordinate,
                    "atom_index": atom,
                    "axis_index": axis,
                    "step_multiplier": multiplier,
                    "displacement_bohr": 0.005 * multiplier,
                }
            )
    return tuple(plan)


def _half_step_matrices(
    new_gradients: dict[tuple[int, int], np.ndarray],
    reused_gradients: dict[tuple[int, int], np.ndarray],
    dimension: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    expected = {
        (coordinate, sign) for coordinate in range(dimension) for sign in (-1, 1)
    }
    if set(new_gradients) != expected or set(reused_gradients) != expected:
        raise ValueError("half-step gradient inventory is incomplete or unexpected")
    inner = np.empty((dimension, dimension), dtype=float)
    outer = np.empty_like(inner)
    for coordinate in range(dimension):
        new_minus = np.asarray(new_gradients[(coordinate, -1)], dtype=float)
        new_plus = np.asarray(new_gradients[(coordinate, 1)], dtype=float)
        old_minus = np.asarray(reused_gradients[(coordinate, -1)], dtype=float)
        old_plus = np.asarray(reused_gradients[(coordinate, 1)], dtype=float)
        if any(
            value.shape != (dimension,) or not np.all(np.isfinite(value))
            for value in (new_minus, new_plus, old_minus, old_plus)
        ):
            raise ValueError("half-step gradient has invalid shape or values")
        inner[:, coordinate] = (new_plus - new_minus) / 0.01
        outer[:, coordinate] = (old_plus - old_minus) / 0.02
    richardson = (4.0 * inner - outer) / 3.0
    return inner, outer, richardson


def _array_record(array: Any, *, path: str, units: str) -> tuple[dict[str, Any], bytes]:
    value = np.ascontiguousarray(array, dtype="<f8")
    if value.size < 1 or not np.all(np.isfinite(value)):
        raise ValueError(f"artifact {path} must be a non-empty finite array")
    raw = value.tobytes()
    return (
        {
            "path": path,
            "dtype": "little-endian float64",
            "shape": list(value.shape),
            "units": units,
            "sha256": hashlib.sha256(raw).hexdigest(),
        },
        raw,
    )


def _publish_record_bundle(
    output_root: Path,
    relative_directory: Path,
    record: dict[str, Any],
    arrays: dict[str, tuple[Any, str]],
    *,
    claim: _OutputRootClaim | None = None,
) -> tuple[dict[str, Any], str]:
    """Atomically publish one receipt and all of its raw arrays as a directory."""

    if claim is not None:
        claim.verify()
    parent = output_root / relative_directory.parent
    parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    destination = output_root / relative_directory
    if destination.exists() or destination.is_symlink():
        raise FileExistsError(f"durable artifact bundle already exists: {destination}")
    stage = parent / f".{relative_directory.name}.{os.getpid()}.{time.time_ns()}.tmp"
    stage.mkdir(mode=0o700)
    try:
        artifacts = {}
        for filename, (array, units) in arrays.items():
            relative_path = str(relative_directory / filename)
            artifact, raw = _array_record(array, path=relative_path, units=units)
            with (stage / filename).open("xb") as handle:
                handle.write(raw)
                handle.flush()
                os.fsync(handle.fileno())
            artifacts[filename] = artifact
        payload = {**record, "artifacts": artifacts}
        receipt_raw = _json_bytes(payload)
        with (stage / "receipt.json").open("xb") as handle:
            handle.write(receipt_raw)
            handle.flush()
            os.fsync(handle.fileno())
        stage_fd = os.open(stage, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(stage_fd)
        finally:
            os.close(stage_fd)
        if claim is not None:
            claim.verify()
        campaign._renameat2_noreplace(stage, destination)
        parent_fd = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(parent_fd)
        finally:
            os.close(parent_fd)
        if claim is not None:
            claim.verify()
        return payload, hashlib.sha256(receipt_raw).hexdigest()
    except BaseException:
        if stage.is_dir():
            for child in stage.iterdir():
                child.unlink()
            stage.rmdir()
        raise


def _read_record_bundle(
    output_root: Path,
    relative_directory: Path,
    *,
    expected_receipt_keys: set[str] | None = None,
    expected_artifacts: dict[str, tuple[tuple[int, ...], str]] | None = None,
) -> tuple[dict[str, Any], str, dict[str, np.ndarray]]:
    directory = output_root / relative_directory
    if directory.is_symlink() or not directory.is_dir():
        raise ValueError(f"artifact bundle must be a regular directory: {directory}")
    receipt, receipt_raw = campaign._read_json_object(
        directory / "receipt.json",
        label=f"{relative_directory} receipt",
        max_bytes=1024 * 1024,
    )
    artifacts = receipt.get("artifacts")
    if type(artifacts) is not dict or not artifacts:
        raise ValueError(f"{relative_directory} artifact inventory is invalid")
    if expected_receipt_keys is not None and set(receipt) != expected_receipt_keys:
        raise ValueError(f"{relative_directory} receipt schema is not exact")
    if expected_artifacts is not None and set(artifacts) != set(expected_artifacts):
        raise ValueError(f"{relative_directory} artifact inventory is not exact")
    arrays = {}
    expected_children = {"receipt.json"}
    for filename, artifact in artifacts.items():
        if type(filename) is not str or Path(filename).name != filename:
            raise ValueError(f"{relative_directory} artifact filename is invalid")
        if type(artifact) is not dict or set(artifact) != {
            "path",
            "dtype",
            "shape",
            "units",
            "sha256",
        }:
            raise ValueError(f"{relative_directory}/{filename} metadata is invalid")
        expected_path = str(relative_directory / filename)
        if (
            artifact.get("path") != expected_path
            or artifact.get("dtype") != "little-endian float64"
            or type(artifact.get("shape")) is not list
            or type(artifact.get("units")) is not str
            or type(artifact.get("sha256")) is not str
        ):
            raise ValueError(f"{relative_directory}/{filename} metadata is invalid")
        shape = artifact["shape"]
        if not shape or any(type(size) is not int or size < 1 for size in shape):
            raise ValueError(f"{relative_directory}/{filename} shape is invalid")
        if expected_artifacts is not None:
            expected_shape, expected_units = expected_artifacts[filename]
            if shape != list(expected_shape) or artifact["units"] != expected_units:
                raise ValueError(
                    f"{relative_directory}/{filename} shape or units mismatch"
                )
        path = directory / filename
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"{relative_directory}/{filename} must be a regular file")
        raw = campaign._read_bounded_regular_snapshot(
            path,
            label=f"{relative_directory}/{filename}",
            maximum_bytes=16 * 1024 * 1024,
        )
        if hashlib.sha256(raw).hexdigest() != artifact["sha256"]:
            raise ValueError(f"{relative_directory}/{filename} hash mismatch")
        if len(raw) != math.prod(shape) * 8:
            raise ValueError(f"{relative_directory}/{filename} byte length mismatch")
        array = np.frombuffer(raw, dtype="<f8").reshape(shape)
        if not np.all(np.isfinite(array)):
            raise ValueError(f"{relative_directory}/{filename} is non-finite")
        arrays[filename] = array
        expected_children.add(filename)
    if {child.name for child in directory.iterdir()} != expected_children:
        raise ValueError(f"{relative_directory} contains unexpected artifacts")
    return receipt, hashlib.sha256(receipt_raw).hexdigest(), arrays


def _fd_settings(base_settings: Any) -> Any:
    settings = replace(base_settings, use_gpu=False, grid_level=5)
    observed = {
        "backend": "pyscf-cpu",
        "xc": str(settings.xc).lower(),
        "dispersion": str(settings.dispersion).lower(),
        "basis": str(settings.basis).lower(),
        "density_fit": settings.density_fit,
        "grid_level": settings.grid_level,
        "scf_tolerance": 1.0e-12,
        "scf_max_cycle": 150,
        "gradient_grid_response": True,
        "density_initial_guess": "one identical center converged density matrix",
    }
    if observed != FINITE_DIFFERENCE_CONTRACT["method"]:
        raise ValueError("finite-difference electronic-structure settings drifted")
    if settings.solvent is not None or settings.composite is not None:
        raise ValueError("finite-difference settings include an unexpected method term")
    return settings


def _evaluate_fd_gradient(
    cluster: Any, settings: Any, *, central_density: np.ndarray | None
) -> tuple[dict[str, Any], np.ndarray, np.ndarray]:
    """Run one bounded CPU SCF/gradient; tests replace PySCF below this seam."""

    started = time.monotonic()
    mf = _make_scf(build_mol(cluster, settings), settings)
    mf.conv_tol = 1.0e-12
    mf.max_cycle = 150
    if getattr(mf.grids, "level", None) != 5:
        raise RuntimeError("PySCF grid level 5 was not established")
    density_sha = None
    if central_density is None:
        electronic = float(mf.kernel())
    else:
        density_raw = np.ascontiguousarray(central_density, dtype="<f8").tobytes()
        density_sha = hashlib.sha256(density_raw).hexdigest()
        electronic = float(mf.kernel(dm0=central_density))
        if (
            hashlib.sha256(
                np.ascontiguousarray(central_density, dtype="<f8").tobytes()
            ).hexdigest()
            != density_sha
        ):
            raise RuntimeError("SCF mutated the shared central density-matrix guess")
    if not math.isfinite(electronic):
        raise ValueError("SCF energy is non-finite")
    if bool(getattr(mf, "converged", False)) is not True:
        raise RuntimeError("SCF did not converge for finite-difference gradient")
    gradient_method = _gradient_method(mf, settings)
    gradient_method.grid_response = True
    if gradient_method.grid_response is not True:
        raise RuntimeError("gradient grid_response=true could not be established")
    gradient_value = gradient_method.kernel()
    if hasattr(gradient_value, "get"):
        gradient_value = gradient_value.get()
    gradient = np.asarray(gradient_value, dtype=float)
    expected_shape = np.asarray(cluster.coords).shape
    if gradient.shape != expected_shape or not np.all(np.isfinite(gradient)):
        raise ValueError("finite-difference gradient shape or values are invalid")
    grid_count = int(np.asarray(mf.grids.coords).shape[0])
    if not hasattr(mf, "spin_square"):
        raise ValueError("SCF spin_square result is missing")
    spin = mf.spin_square()
    if not isinstance(spin, (tuple, list)) or len(spin) != 2:
        raise ValueError("SCF spin_square result is incomplete")
    spin_square = [_finite_number(value, label="SCF spin result") for value in spin]
    density = np.asarray(mf.make_rdm1(), dtype=float)
    if density.size < 1 or not np.all(np.isfinite(density)):
        raise ValueError("final density matrix is empty or non-finite")
    density = np.frombuffer(
        np.ascontiguousarray(density, dtype="<f8").tobytes(), dtype="<f8"
    ).reshape(density.shape)
    final_density_sha = hashlib.sha256(density.tobytes()).hexdigest()
    if central_density is None:
        density_sha = final_density_sha
    record = {
        "backend": "pyscf-cpu",
        "electronic_hartree": electronic,
        "scf_converged": True,
        "grid_point_count": grid_count,
        "spin_2s": cluster.spin,
        "spin_square": spin_square,
        "gradient_grid_response": True,
        "density_initial_guess_sha256": density_sha,
        "density_final_sha256": final_density_sha,
        "geometry_fingerprint": frequency_geometry_fingerprint(cluster),
        "elapsed_seconds": time.monotonic() - started,
    }
    return record, gradient, density


def _validate_point_record(
    record: dict[str, Any],
    *,
    center_spin_square: float | None,
    center: bool,
    enforce_scientific_gates: bool = True,
) -> float:
    thresholds = FINITE_DIFFERENCE_CONTRACT["thresholds"]
    if record.get("backend") != "pyscf-cpu":
        raise ValueError("finite-difference backend must be CPU PySCF")
    if type(record.get("scf_converged")) is not bool:
        raise ValueError("finite-difference SCF convergence evidence is malformed")
    if record.get("scf_converged") is not True and enforce_scientific_gates:
        raise ScientificRejection("finite-difference SCF convergence gate failed")
    if record.get("gradient_grid_response") is not True:
        raise ValueError("finite-difference gradient grid-response gate failed")
    if type(record.get("grid_point_count")) is not int:
        raise ValueError("finite-difference grid point count evidence is malformed")
    if enforce_scientific_gates and (
        record["grid_point_count"] != thresholds["required_grid_point_count"]
    ):
        raise ScientificRejection(
            "finite-difference grid point count must equal 199560"
        )
    _finite_number(record.get("electronic_hartree"), label="finite-difference energy")
    if type(record.get("spin_2s")) is not int or record["spin_2s"] != 1:
        raise ValueError(
            "finite-difference electronic spin must be the exact doublet 2S=1"
        )
    spin = record.get("spin_square")
    if type(spin) is not list or len(spin) != 2:
        raise ValueError("finite-difference spin receipt is missing or incomplete")
    s2 = _finite_number(spin[0], label="finite-difference S2")
    _finite_number(spin[1], label="finite-difference spin multiplicity")
    if enforce_scientific_gates and (
        not thresholds["spin_square_inclusive_minimum"]
        <= s2
        <= thresholds["spin_square_inclusive_maximum"]
    ):
        raise ScientificRejection(
            "finite-difference S2 is outside the inclusive [0.74, 0.80] range"
        )
    if center and enforce_scientific_gates:
        _strict_less(
            record.get("physical_fmax_ev_per_angstrom"),
            thresholds["center_fmax_ev_per_angstrom_exclusive_maximum"],
            label="center physical fmax",
        )
    elif not center and enforce_scientific_gates:
        if center_spin_square is None:
            raise ValueError("center S2 is required before displaced-point validation")
        _strict_less(
            abs(s2 - center_spin_square),
            thresholds["displaced_vs_center_spin_square_exclusive_maximum_delta"],
            label="displaced-vs-center S2 delta",
        )
    return s2


def _density_continuity_metrics(density: Any, center_density: Any) -> dict[str, float]:
    observed = np.asarray(density, dtype=float)
    center = np.asarray(center_density, dtype=float)
    if (
        observed.shape != center.shape
        or observed.size < 1
        or not np.all(np.isfinite(observed))
        or not np.all(np.isfinite(center))
    ):
        raise ValueError("final and center density matrices must have one finite shape")
    observed_norm = float(np.linalg.norm(observed))
    center_norm = float(np.linalg.norm(center))
    denominator = max(observed_norm, center_norm)
    if not math.isfinite(denominator) or denominator <= 0.0:
        raise ValueError("density continuity scale must be finite and positive")
    delta = float(np.linalg.norm(observed - center) / denominator)
    return {
        "final_vs_center_normalized_frobenius_delta": delta,
        "final_frobenius_norm": observed_norm,
        "center_frobenius_norm": center_norm,
    }


def _enforce_density_continuity(metrics: dict[str, Any]) -> None:
    _strict_less(
        metrics.get("final_vs_center_normalized_frobenius_delta"),
        FINITE_DIFFERENCE_CONTRACT["thresholds"][
            "displaced_final_density_vs_center_exclusive_maximum_normalized_frobenius_delta"
        ],
        label="displaced-final-vs-center normalized Frobenius density delta",
    )


def _center_reference_metrics(
    record: dict[str, Any], reference_case: dict[str, Any]
) -> dict[str, float]:
    energy = _finite_number(
        record.get("electronic_hartree"), label="center electronic energy"
    )
    reference_energy = _finite_number(
        reference_case.get("electronic_hartree"),
        label="reference center electronic energy",
    )
    spin = record.get("spin_square")
    reference_spin = reference_case.get("spin_square")
    if type(spin) is not list or len(spin) != 2:
        raise ValueError("center spin receipt is incomplete")
    if type(reference_spin) is not list or len(reference_spin) != 2:
        raise ValueError("reference center spin receipt is incomplete")
    s2 = _finite_number(spin[0], label="center S2")
    reference_s2 = _finite_number(reference_spin[0], label="reference center S2")
    thresholds = FINITE_DIFFERENCE_CONTRACT["thresholds"]
    if not (
        thresholds["spin_square_inclusive_minimum"]
        <= reference_s2
        <= thresholds["spin_square_inclusive_maximum"]
    ):
        raise ValueError("pinned reference center S2 is outside its valid range")
    energy_delta = abs(energy - reference_energy)
    s2_delta = abs(s2 - reference_s2)
    return {
        "reference_electronic_hartree": reference_energy,
        "electronic_energy_absolute_delta_hartree": energy_delta,
        "reference_s2": reference_s2,
        "s2_absolute_delta": s2_delta,
    }


def _enforce_center_reference_metrics(metrics: dict[str, Any]) -> None:
    thresholds = FINITE_DIFFERENCE_CONTRACT["thresholds"]
    _strict_less(
        metrics.get("electronic_energy_absolute_delta_hartree"),
        thresholds["center_energy_vs_reference_exclusive_maximum_delta_hartree"],
        label="center energy delta from pinned reference",
    )
    _strict_less(
        metrics.get("s2_absolute_delta"),
        thresholds["center_s2_vs_reference_exclusive_maximum_delta"],
        label="center S2 delta from pinned reference",
    )


def _reference_settings_fingerprint() -> str:
    return campaign._canonical_hash(
        {
            "settings": campaign.DFT_SETTINGS,
            "actual_grid_level": 5,
            "scf_tolerance": 1.0e-12,
            "cpscf_tolerance": 1.0e-12,
            "hessian_max_cycle": 150,
            "backend": "pyscf",
        }
    )


def _current_reference_binding(
    preflight: dict[str, Any],
    transition_state: Any,
    reactant: Any,
    product: Any,
    masses: np.ndarray,
) -> dict[str, Any]:
    """Materialize every current input whose identity the reference must share."""

    if type(preflight) is not dict:
        raise ValueError("current preflight campaign identity is invalid")
    campaign_identity = campaign._require_json_string(
        preflight.get("identity"), label="current preflight campaign identity"
    )
    campaign._require_sha(
        campaign_identity, length=64, label="current preflight campaign identity"
    )
    route_record = preflight.get("routes", {}).get(ROUTE)
    if type(route_record) is not dict:
        raise ValueError("current preflight route identity is missing")
    fingerprints = campaign._trusted_input_fingerprints_from_route_record(route_record)
    transition_fingerprint = frequency_geometry_fingerprint(transition_state)
    reactant_fingerprint = frequency_geometry_fingerprint(reactant)
    product_fingerprint = frequency_geometry_fingerprint(product)
    reactant_coordinate_sha256 = hashlib.sha256(
        np.ascontiguousarray(reactant.coords, dtype="<f8").tobytes()
    ).hexdigest()
    product_coordinate_sha256 = hashlib.sha256(
        np.ascontiguousarray(product.coords, dtype="<f8").tobytes()
    ).hexdigest()
    if (
        transition_fingerprint
        != fingerprints["transition_state"]["geometry_fingerprint"]
    ):
        raise ValueError("current frozen TS geometry fingerprint drifted")
    if reactant_fingerprint != fingerprints["reactant"]["geometry_fingerprint"]:
        raise ValueError("current mapped reactant geometry fingerprint drifted")
    if product_fingerprint != fingerprints["product"]["geometry_fingerprint"]:
        raise ValueError("current mapped product geometry fingerprint drifted")
    symbols = list(transition_state.symbols)
    atom_identity_labels = route_record.get("atom_identity_labels")
    if route_record.get("symbols") != symbols:
        raise ValueError("current route atom symbols drifted")
    if (
        type(atom_identity_labels) is not list
        or len(atom_identity_labels) != len(symbols)
        or any(type(label) is not str or not label for label in atom_identity_labels)
    ):
        raise ValueError("current route semantic atom identity is invalid")
    atom_mapping_sha256 = route_record.get("atom_mapping_sha256")
    if (
        type(atom_mapping_sha256) is not str
        or len(atom_mapping_sha256) != 64
        or any(character not in "0123456789abcdef" for character in atom_mapping_sha256)
    ):
        raise ValueError("current route atom-mapping hash is invalid")
    mass_values = np.asarray(masses, dtype=float)
    if (
        mass_values.shape != (len(symbols),)
        or not np.all(np.isfinite(mass_values))
        or np.any(mass_values <= 0.0)
    ):
        raise ValueError("current isotopic masses are invalid")
    campaign._strict_json_equal(
        mass_values.tolist(),
        route_record.get("masses_amu"),
        label="current route isotopic masses",
    )
    method = dict(FINITE_DIFFERENCE_CONTRACT["method"])
    if (
        method.get("xc") != str(campaign.DFT_SETTINGS["xc"]).lower()
        or method.get("basis") != str(campaign.DFT_SETTINGS["basis"]).lower()
        or method.get("dispersion") != str(campaign.DFT_SETTINGS["dispersion"]).lower()
        or method.get("density_fit") is not campaign.DFT_SETTINGS["density_fit"]
    ):
        raise ValueError("current FD method disagrees with the frozen campaign method")
    return {
        "campaign_identity": campaign_identity,
        "route": ROUTE,
        "symbols": symbols,
        "atom_identity_labels": atom_identity_labels,
        "atom_mapping_sha256": atom_mapping_sha256,
        "transition_state_file_sha256": fingerprints["transition_state"]["file_sha256"],
        "transition_state_geometry_sha256": fingerprints["transition_state"][
            "geometry_sha256"
        ],
        "transition_state_geometry_fingerprint": transition_fingerprint,
        "mapped_reactant_file_sha256": fingerprints["reactant"]["file_sha256"],
        "mapped_reactant_geometry_sha256": fingerprints["reactant"]["geometry_sha256"],
        "mapped_reactant_geometry_fingerprint": reactant_fingerprint,
        "mapped_reactant_coordinate_sha256": reactant_coordinate_sha256,
        "mapped_product_file_sha256": fingerprints["product"]["file_sha256"],
        "mapped_product_geometry_sha256": fingerprints["product"]["geometry_sha256"],
        "mapped_product_geometry_fingerprint": product_fingerprint,
        "mapped_product_coordinate_sha256": product_coordinate_sha256,
        "isotopic_masses_amu": mass_values.tolist(),
        "fd_method": method,
        "analytic_reference_settings_fingerprint": _reference_settings_fingerprint(),
    }


def _load_reference(
    reference_root: Path,
    atom_count: int,
    *,
    preflight: dict[str, Any],
    transition_state: Any,
    reactant: Any,
    product: Any,
    masses: np.ndarray,
) -> tuple[dict[str, Any], np.ndarray, dict[str, Any], dict[str, Any]]:
    current_binding = _current_reference_binding(
        preflight, transition_state, reactant, product, masses
    )
    root = campaign._safe_absolute_root(reference_root)
    receipt, raw = campaign._read_json_object(
        root / "receipt.json",
        label="analytic Hessian reference receipt",
        max_bytes=1024 * 1024,
    )
    receipt_sha = hashlib.sha256(raw).hexdigest()
    if receipt_sha != REFERENCE_RECEIPT_SHA256:
        raise ValueError("analytic Hessian reference receipt SHA-256 mismatch")
    expected_receipt_keys = {
        "accepted_campaign_result",
        "campaign_identity",
        "cases",
        "comparisons_to_reference",
        "finished_utc",
        "git_sha",
        "preflight_receipt_sha256",
        "purpose",
        "route",
        "schema",
        "script_sha256",
        "state",
    }
    if (
        set(receipt) != expected_receipt_keys
        or receipt.get("schema") != SCHEMA
        or receipt.get("state") != "completed"
        or receipt.get("accepted_campaign_result") is not False
        or receipt.get("route") != ROUTE
    ):
        raise ValueError("analytic Hessian reference receipt identity is invalid")
    for field, length in (
        ("campaign_identity", 64),
        ("git_sha", 40),
        ("preflight_receipt_sha256", 64),
        ("script_sha256", 64),
    ):
        campaign._require_sha(
            campaign._require_json_string(
                receipt.get(field), label=f"analytic reference {field}"
            ),
            length=length,
            label=f"analytic reference {field}",
        )
    cases = receipt.get("cases")
    if type(cases) is not list:
        raise ValueError("analytic Hessian reference cases are missing")
    matches = [
        case
        for case in cases
        if type(case) is dict and case.get("name") == "D-dense-strict-reference"
    ]
    if len(matches) != 1:
        raise ValueError("analytic dense-strict reference case must occur exactly once")
    case = matches[0]
    expected_case_keys = {
        "components",
        "cpscf_tolerance",
        "elapsed_seconds",
        "electronic_hartree",
        "grid_level",
        "grid_point_count",
        "hessian_max_cycle",
        "name",
        "near_zero_mode_count_at_1e-8",
        "negative_mode_count_below_1e-8",
        "physical_fmax_ev_per_angstrom",
        "post_symmetrization_gate",
        "projected_eigenvalues_hartree_per_bohr2_amu",
        "scf_converged",
        "scf_tolerance",
        "signed_wavenumbers_cm",
        "spin_square",
        "symmetric_total_artifact",
    }
    gate = case.get("post_symmetrization_gate")
    evidence = gate.get("evidence") if type(gate) is dict else None
    expected_evidence_keys = {
        "accepted",
        "actual_backend",
        "canonical_hessian_sha256",
        "eigenvalues_hartree_per_bohr2_amu",
        "electronic_hartree",
        "geometry_fingerprint",
        "gpu_fallback_used",
        "gradient_sha256",
        "imaginary_mode_count",
        "imaginary_wavenumber_cm",
        "mapped_product_geometry_sha256",
        "mapped_reactant_geometry_sha256",
        "mapped_reaction_vector_overlap",
        "masses_amu",
        "minimum_imaginary_wavenumber_cm",
        "minimum_mapped_reaction_vector_overlap",
        "physical_fmax_ev_per_angstrom",
        "physical_fmax_exclusive_limit_ev_per_angstrom",
        "reaction_vector_route_binding_overlap",
        "reaction_vector_sha256",
        "reaction_vector_source",
        "requested_backend",
        "settings_fingerprint",
        "unstable_mode_mass_scaled",
        "vibrational_mode_count",
    }
    if (
        set(case) != expected_case_keys
        or type(gate) is not dict
        or set(gate) != {"accepted", "evidence"}
        or gate.get("accepted") is not True
        or type(evidence) is not dict
        or set(evidence) != expected_evidence_keys
        or evidence.get("accepted") is not True
        or evidence.get("actual_backend") != "pyscf"
        or evidence.get("requested_backend") != "pyscf"
        or evidence.get("gpu_fallback_used") is not False
        or case.get("grid_level") != 5
        or case.get("scf_tolerance") != 1.0e-12
        or case.get("cpscf_tolerance") != 1.0e-12
        or case.get("hessian_max_cycle") != 150
        or case.get("grid_point_count") != 199560
        or case.get("scf_converged") is not True
        or evidence.get("geometry_fingerprint")
        != current_binding.get("transition_state_geometry_fingerprint")
        or evidence.get("settings_fingerprint")
        != current_binding.get("analytic_reference_settings_fingerprint")
        or evidence.get("canonical_hessian_sha256") != REFERENCE_MATRIX_SHA256
        or evidence.get("masses_amu") != current_binding.get("isotopic_masses_amu")
        or evidence.get("mapped_reactant_geometry_sha256")
        != current_binding.get("mapped_reactant_coordinate_sha256")
        or evidence.get("mapped_product_geometry_sha256")
        != current_binding.get("mapped_product_coordinate_sha256")
        or evidence.get("electronic_hartree") != case.get("electronic_hartree")
        or evidence.get("physical_fmax_ev_per_angstrom")
        != case.get("physical_fmax_ev_per_angstrom")
    ):
        raise ValueError(
            "analytic dense-strict reference geometry or settings identity is invalid"
        )
    _center_reference_metrics(case, case)
    artifact = case.get("symmetric_total_artifact")
    expected_path = "D-dense-strict-reference/total-symmetric.f64"
    expected = {
        "path": expected_path,
        "dtype": "little-endian float64",
        "shape": [3 * atom_count, 3 * atom_count],
        "units": "hartree / bohr^2",
        "sha256": REFERENCE_MATRIX_SHA256,
    }
    if artifact != expected:
        raise ValueError("analytic dense-strict matrix metadata or SHA-256 mismatch")
    path = root / expected_path
    if path.is_symlink() or not path.is_file():
        raise ValueError("analytic dense-strict matrix must be a regular file")
    matrix_raw = campaign._read_bounded_regular_snapshot(
        path,
        label="analytic dense-strict matrix",
        maximum_bytes=1024 * 1024,
    )
    if hashlib.sha256(matrix_raw).hexdigest() != REFERENCE_MATRIX_SHA256:
        raise ValueError("analytic dense-strict matrix content SHA-256 mismatch")
    if len(matrix_raw) != (3 * atom_count) ** 2 * 8:
        raise ValueError("analytic dense-strict matrix byte length mismatch")
    matrix = np.frombuffer(matrix_raw, dtype="<f8").reshape(
        3 * atom_count, 3 * atom_count
    )
    if not np.all(np.isfinite(matrix)) or not np.array_equal(matrix, matrix.T):
        raise ValueError(
            "analytic dense-strict matrix must be finite and exactly symmetric"
        )
    reference_identity = {
        "receipt_sha256": receipt_sha,
        "matrix_sha256": REFERENCE_MATRIX_SHA256,
        "campaign_identity": receipt["campaign_identity"],
        "preflight_receipt_sha256": receipt["preflight_receipt_sha256"],
        "git_sha": receipt["git_sha"],
        "script_sha256": receipt["script_sha256"],
    }
    return (
        receipt,
        matrix,
        case,
        {**current_binding, "analytic_reference_identity": reference_identity},
    )


def _matrix_delta_metrics(
    matrix: np.ndarray, reference: np.ndarray
) -> dict[str, float]:
    value = np.asarray(matrix, dtype=float)
    baseline = np.asarray(reference, dtype=float)
    delta = value - baseline
    if not np.all(np.isfinite(delta)):
        raise ValueError("matrix delta is non-finite")
    denominator = max(
        float(np.linalg.norm(value, ord=2)),
        float(np.linalg.norm(baseline, ord=2)),
    )
    if not math.isfinite(denominator) or denominator <= 0.0:
        raise ValueError("matrix comparison reference spectral norm is invalid")
    return {
        "maximum_absolute_delta": float(np.max(np.abs(delta))),
        "spectral_norm_delta": float(np.linalg.norm(delta, ord=2)),
        "spectral_relative_delta": float(np.linalg.norm(delta, ord=2) / denominator),
    }


def _fd_symmetry_metrics(matrix: np.ndarray) -> dict[str, Any]:
    value = np.asarray(matrix, dtype=float)
    if (
        value.ndim != 2
        or value.shape[0] != value.shape[1]
        or not np.all(np.isfinite(value))
    ):
        raise ValueError("raw finite-difference Hessian must be a finite square matrix")
    symmetric = 0.5 * (value + value.T)
    asymmetry = value - value.T
    symmetric_norm = float(np.linalg.norm(symmetric, ord=2))
    asymmetry_norm = float(np.linalg.norm(asymmetry, ord=2))
    if symmetric_norm <= 0.0 or not math.isfinite(symmetric_norm):
        raise ValueError("raw finite-difference Hessian has invalid symmetric scale")
    maximum = float(np.max(np.abs(asymmetry)))
    relative = asymmetry_norm / symmetric_norm
    thresholds = FINITE_DIFFERENCE_CONTRACT["thresholds"]
    return {
        "maximum_absolute_asymmetry": maximum,
        "asymmetry_spectral_norm": asymmetry_norm,
        "symmetric_spectral_norm": symmetric_norm,
        "spectral_relative_asymmetry": relative,
        "maximum_absolute_asymmetry_exclusive_limit": thresholds[
            "raw_fd_maximum_absolute_asymmetry_exclusive_maximum"
        ],
        "spectral_relative_asymmetry_exclusive_limit": thresholds[
            "raw_fd_spectral_relative_asymmetry_exclusive_maximum"
        ],
    }


def _maximum_overlap_assignment(overlaps: np.ndarray) -> list[int]:
    """Exact maximum-weight one-to-one assignment without a new runtime dependency."""

    matrix = np.asarray(overlaps, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or matrix.shape[0] < 1:
        raise ValueError("mode-overlap matrix must be non-empty and square")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("mode-overlap matrix must be finite")
    size = matrix.shape[0]
    scores: dict[int, tuple[float, tuple[int, ...]]] = {0: (0.0, ())}
    for row in range(size):
        next_scores = {}
        for mask, (score, assignment) in scores.items():
            for column in range(size):
                bit = 1 << column
                if mask & bit:
                    continue
                candidate = (score + float(matrix[row, column]), assignment + (column,))
                previous = next_scores.get(mask | bit)
                if previous is None or candidate[0] > previous[0]:
                    next_scores[mask | bit] = candidate
        scores = next_scores
    return list(scores[(1 << size) - 1][1])


def _mode_comparison(left: Any, right: Any) -> dict[str, Any]:
    left_vectors = np.asarray(left.mass_weighted_eigenvectors, dtype=float).reshape(
        len(left.eigenvalues), -1
    )
    right_vectors = np.asarray(right.mass_weighted_eigenvectors, dtype=float).reshape(
        len(right.eigenvalues), -1
    )
    if left_vectors.shape != right_vectors.shape:
        raise ValueError("mode sets have incompatible dimensions")
    overlaps = np.abs(left_vectors @ right_vectors.T)
    assignment = _maximum_overlap_assignment(overlaps)
    ordered = np.diag(overlaps)
    return {
        "maximum_overlap_assignment": assignment,
        "assignment_equals_eigenvalue_order": assignment
        == list(range(len(assignment))),
        "ordered_overlaps": ordered.tolist(),
        "assigned_overlaps": [
            float(overlaps[row, column]) for row, column in enumerate(assignment)
        ],
    }


def _spectrum_payload(name: str, modes: Any) -> dict[str, Any]:
    thresholds = FINITE_DIFFERENCE_CONTRACT["thresholds"]
    eigenvalues = np.asarray(modes.eigenvalues, dtype=float)
    if (
        eigenvalues.ndim != 1
        or eigenvalues.size < 2
        or not np.all(np.isfinite(eigenvalues))
    ):
        raise ValueError(f"{name} projected spectrum is invalid")
    negative = eigenvalues < thresholds["negative_eigenvalue_threshold"]
    near_zero = (
        np.abs(eigenvalues)
        <= thresholds["near_zero_eigenvalue_inclusive_maximum_absolute"]
    )
    positive = (
        eigenvalues > thresholds["near_zero_eigenvalue_inclusive_maximum_absolute"]
    )
    if np.count_nonzero(negative) != 1:
        raise ScientificRejection(
            f"{name} must have exactly one eigenvalue below -1e-8"
        )
    if np.any(near_zero):
        raise ScientificRejection(
            f"{name} contains an eigenvalue at the inclusive noise boundary"
        )
    if not np.any(positive):
        raise ScientificRejection(f"{name} has no positive vibrational eigenvalue")
    lowest_positive = float(np.min(eigenvalues[positive]))
    _strict_greater(
        lowest_positive,
        thresholds["lowest_positive_eigenvalue_exclusive_minimum"],
        label=f"{name} lowest positive eigenvalue",
    )
    frequencies = hessian_eigenvalues_to_wavenumbers_cm(eigenvalues)
    imaginary = float(abs(frequencies[np.flatnonzero(negative)[0]]))
    _strict_greater(
        imaginary,
        thresholds["imaginary_wavenumber_cm_exclusive_minimum"],
        label=f"{name} imaginary frequency cm^-1",
    )
    return {
        "eigenvalues_hartree_per_bohr2_amu": eigenvalues.tolist(),
        "signed_wavenumbers_cm": frequencies.tolist(),
        "negative_mode_count_below_negative_1e_minus_8": 1,
        "near_zero_mode_count_at_1e_minus_8": 0,
        "lowest_positive_eigenvalue": lowest_positive,
        "imaginary_wavenumber_cm": imaginary,
    }


def _confirmation_analysis(
    transition_state: Any,
    masses: np.ndarray,
    route_vector: np.ndarray,
    h_h: np.ndarray,
    h_2h: np.ndarray,
    raw_fd: np.ndarray,
    analytic: np.ndarray,
) -> dict[str, Any]:
    thresholds = FINITE_DIFFERENCE_CONTRACT["thresholds"]
    sym_h_h = 0.5 * (h_h + h_h.T)
    sym_h_2h = 0.5 * (h_2h + h_2h.T)
    sym_fd = 0.5 * (raw_fd + raw_fd.T)
    symmetry = _fd_symmetry_metrics(raw_fd)
    _strict_less(
        symmetry["maximum_absolute_asymmetry"],
        thresholds["raw_fd_maximum_absolute_asymmetry_exclusive_maximum"],
        label="raw FD maximum asymmetry",
    )
    _strict_less(
        symmetry["spectral_relative_asymmetry"],
        thresholds["raw_fd_spectral_relative_asymmetry_exclusive_maximum"],
        label="raw FD spectral-relative asymmetry",
    )
    richardson_delta = _matrix_delta_metrics(sym_fd, sym_h_h)
    _strict_less(
        richardson_delta["maximum_absolute_delta"],
        thresholds[
            "symmetric_richardson_vs_h_h_maximum_absolute_delta_exclusive_maximum"
        ],
        label="symmetric Richardson-vs-H_h maximum delta",
    )
    _strict_less(
        richardson_delta["spectral_relative_delta"],
        thresholds[
            "symmetric_richardson_vs_h_h_spectral_relative_delta_exclusive_maximum"
        ],
        label="symmetric Richardson-vs-H_h spectral-relative delta",
    )
    analytic_delta = _matrix_delta_metrics(sym_fd, analytic)
    _strict_less(
        analytic_delta["maximum_absolute_delta"],
        thresholds["fd_vs_analytic_symmetric_maximum_absolute_delta_exclusive_maximum"],
        label="FD-vs-analytic maximum delta",
    )
    _strict_less(
        analytic_delta["spectral_relative_delta"],
        thresholds[
            "fd_vs_analytic_symmetric_spectral_relative_delta_exclusive_maximum"
        ],
        label="FD-vs-analytic spectral-relative delta",
    )

    modes = {
        "H_h": project_vibrational_hessian(transition_state.coords, masses, sym_h_h),
        "H_2h": project_vibrational_hessian(transition_state.coords, masses, sym_h_2h),
        "FD": project_vibrational_hessian(transition_state.coords, masses, sym_fd),
        "analytic_reference": project_vibrational_hessian(
            transition_state.coords, masses, analytic
        ),
    }
    spectra = {
        name: _spectrum_payload(name, modes[name]) for name in ("H_h", "H_2h", "FD")
    }
    reference_eigenvalues = np.asarray(modes["analytic_reference"].eigenvalues)
    eigenvalue_limits = np.maximum(
        thresholds["per_mode_reference_eigenvalue_delta_absolute_floor"],
        thresholds["per_mode_reference_eigenvalue_delta_relative_fraction"]
        * np.abs(reference_eigenvalues),
    )
    comparisons_to_reference = {}
    for name in ("H_h", "H_2h", "FD"):
        comparison = _mode_comparison(modes[name], modes["analytic_reference"])
        if comparison["assignment_equals_eigenvalue_order"] is not True:
            raise ScientificRejection(
                f"{name}-vs-reference maximum-overlap assignment differs "
                "from eigenvalue order"
            )
        eigenvalue_deltas = np.abs(
            np.asarray(modes[name].eigenvalues) - reference_eigenvalues
        )
        if np.any(eigenvalue_deltas >= eigenvalue_limits):
            raise ScientificRejection(
                f"{name} per-mode eigenvalue delta reached its exclusive "
                "reference limit"
            )
        comparisons_to_reference[name] = {
            **comparison,
            "eigenvalue_absolute_deltas": eigenvalue_deltas.tolist(),
            "eigenvalue_exclusive_limits": eigenvalue_limits.tolist(),
        }
    fd_reference = comparisons_to_reference["FD"]
    ordered = fd_reference["ordered_overlaps"]
    if not all(
        value > thresholds["fd_vs_reference_all_ordered_overlaps_exclusive_minimum"]
        for value in ordered
    ):
        raise ScientificRejection(
            "an FD-vs-reference ordered mode overlap is at or below 0.90"
        )
    _strict_greater(
        ordered[0],
        thresholds["fd_vs_reference_unstable_overlap_exclusive_minimum"],
        label="FD-vs-reference unstable overlap",
    )
    _strict_greater(
        ordered[1],
        thresholds["fd_vs_reference_low_positive_overlap_exclusive_minimum"],
        label="FD-vs-reference low-positive overlap",
    )
    fd_frequencies = np.asarray(spectra["FD"]["signed_wavenumbers_cm"])
    reference_frequencies = hessian_eigenvalues_to_wavenumbers_cm(reference_eigenvalues)
    fd_reference_frequency_deltas = np.abs(fd_frequencies - reference_frequencies)
    _strict_less(
        fd_reference_frequency_deltas[0],
        thresholds["fd_vs_reference_imaginary_frequency_delta_cm_exclusive_maximum"],
        label="FD-vs-reference imaginary frequency delta",
    )
    _strict_less(
        fd_reference_frequency_deltas[1],
        thresholds["fd_vs_reference_low_positive_frequency_delta_cm_exclusive_maximum"],
        label="FD-vs-reference low-positive frequency delta",
    )

    h_h_fd = _mode_comparison(modes["H_h"], modes["FD"])
    if h_h_fd["assignment_equals_eigenvalue_order"] is not True:
        raise ScientificRejection(
            "H_h-vs-FD maximum-overlap assignment differs from eigenvalue order"
        )
    _strict_greater(
        h_h_fd["ordered_overlaps"][0],
        thresholds["h_h_vs_fd_unstable_overlap_exclusive_minimum"],
        label="H_h-vs-FD unstable overlap",
    )
    _strict_greater(
        h_h_fd["ordered_overlaps"][1],
        thresholds["h_h_vs_fd_low_positive_overlap_exclusive_minimum"],
        label="H_h-vs-FD low-positive overlap",
    )
    h_2h_fd = _mode_comparison(modes["H_2h"], modes["FD"])
    if h_2h_fd["assignment_equals_eigenvalue_order"] is not True:
        raise ScientificRejection(
            "H_2h-vs-FD maximum-overlap assignment differs from eigenvalue order"
        )
    h_h_frequencies = np.asarray(spectra["H_h"]["signed_wavenumbers_cm"])
    h_h_fd_frequency_deltas = np.abs(h_h_frequencies - fd_frequencies)
    _strict_less(
        h_h_fd_frequency_deltas[0],
        thresholds["h_h_vs_fd_imaginary_frequency_delta_cm_exclusive_maximum"],
        label="H_h-vs-FD imaginary frequency delta",
    )
    _strict_less(
        h_h_fd_frequency_deltas[1],
        thresholds["h_h_vs_fd_low_positive_frequency_delta_cm_exclusive_maximum"],
        label="H_h-vs-FD low-positive frequency delta",
    )

    normalized_route = np.asarray(route_vector, dtype=float).reshape(-1)
    route_norm = float(np.linalg.norm(normalized_route))
    if not math.isfinite(route_norm) or route_norm <= 1.0e-12:
        raise ValueError("mapped reaction vector is invalid")
    normalized_route /= route_norm
    fd_reaction_overlap = float(
        abs(
            np.dot(
                modes["FD"].mass_weighted_eigenvectors[0].reshape(-1), normalized_route
            )
        )
    )
    reference_reaction_overlap = float(
        abs(
            np.dot(
                modes["analytic_reference"].mass_weighted_eigenvectors[0].reshape(-1),
                normalized_route,
            )
        )
    )
    _strict_greater(
        fd_reaction_overlap,
        thresholds["fd_mapped_reaction_overlap_exclusive_minimum"],
        label="FD mapped reaction overlap",
    )
    _strict_less(
        abs(fd_reaction_overlap - reference_reaction_overlap),
        thresholds["fd_mapped_reaction_overlap_vs_reference_exclusive_maximum_delta"],
        label="FD mapped-reaction-overlap delta from reference",
    )
    return {
        "confirmation_passed": True,
        "raw_fd_symmetry": symmetry,
        "symmetric_richardson_vs_h_h": richardson_delta,
        "fd_vs_analytic_symmetric": analytic_delta,
        "spectra": spectra,
        "comparisons_to_analytic_reference": {
            **comparisons_to_reference,
            "FD": {
                **fd_reference,
                "frequency_absolute_deltas_cm": (
                    fd_reference_frequency_deltas.tolist()
                ),
            },
        },
        "h_h_vs_fd_modes": {
            **h_h_fd,
            "frequency_absolute_deltas_cm": h_h_fd_frequency_deltas.tolist(),
        },
        "h_2h_vs_fd_modes": h_2h_fd,
        "mapped_reaction_overlap": {
            "fd": fd_reaction_overlap,
            "analytic_reference": reference_reaction_overlap,
            "absolute_delta": abs(fd_reaction_overlap - reference_reaction_overlap),
        },
        "symmetric_matrices": {
            "H_h": sym_h_h,
            "H_2h": sym_h_2h,
            "FD": sym_fd,
        },
    }


def _scientific_analysis_verdict(
    transition_state: Any,
    masses: np.ndarray,
    route_vector: np.ndarray,
    center_record: dict[str, Any],
    center_density: np.ndarray,
    displaced_points: list[tuple[dict[str, Any], np.ndarray]],
    h_h: np.ndarray,
    h_2h: np.ndarray,
    raw_fd: np.ndarray,
    analytic_matrix: np.ndarray,
) -> tuple[bool, dict[str, Any]]:
    """Evaluate every declared gate after all durable evidence exists."""

    symmetric = {
        "H_h": 0.5 * (h_h + h_h.T),
        "H_2h": 0.5 * (h_2h + h_2h.T),
        "FD": 0.5 * (raw_fd + raw_fd.T),
    }
    try:
        center_spin = _validate_point_record(
            center_record, center_spin_square=None, center=True
        )
        _enforce_center_reference_metrics(
            {key: center_record.get(key) for key in _CENTER_REFERENCE_KEYS}
        )
        for point, density in displaced_points:
            _validate_point_record(point, center_spin_square=center_spin, center=False)
            expected_density = _density_continuity_metrics(density, center_density)
            campaign._strict_json_equal(
                {key: point.get(key) for key in _DENSITY_CONTINUITY_KEYS},
                expected_density,
                label=f"scientific density continuity {point.get('point_key')}",
            )
            _enforce_density_continuity(expected_density)
        analysis = _confirmation_analysis(
            transition_state,
            masses,
            route_vector,
            h_h,
            h_2h,
            raw_fd,
            analytic_matrix,
        )
    except ScientificRejection as exc:
        return False, {
            "confirmation_passed": False,
            "rejection": f"{type(exc).__name__}: {exc}",
        }
    if type(analysis) is not dict or analysis.get("confirmation_passed") is not True:
        raise RuntimeError("confirmation analysis returned an invalid passing verdict")
    analysis = dict(analysis)
    observed_symmetric = analysis.pop("symmetric_matrices", None)
    if type(observed_symmetric) is not dict or any(
        name not in observed_symmetric
        or not np.array_equal(observed_symmetric[name], symmetric[name])
        for name in symmetric
    ):
        raise RuntimeError(
            "confirmation analysis returned divergent symmetric matrices"
        )
    return True, analysis


def _analysis_gate(
    gates: list[dict[str, Any]],
    name: str,
    value: Any,
    limit: Any,
    comparison: str,
) -> None:
    """Record, rather than raise from, one strict scientific boundary."""

    if comparison == "strictly_less_than":
        if value is None:
            observed, passed = None, False
        else:
            observed = _finite_number(value, label=name)
            passed = bool(observed < float(limit))
    elif comparison == "strictly_greater_than":
        if value is None:
            observed, passed = None, False
        else:
            observed = _finite_number(value, label=name)
            passed = bool(observed > float(limit))
    elif comparison == "greater_than_or_equal":
        observed = _finite_number(value, label=name)
        passed = bool(observed >= float(limit))
    elif comparison == "less_than_or_equal":
        observed = _finite_number(value, label=name)
        passed = bool(observed <= float(limit))
    elif comparison == "equals":
        observed = value
        passed = bool(value == limit)
    elif comparison == "all_strictly_less_than":
        values = np.asarray(value, dtype=float)
        limits = np.asarray(limit, dtype=float)
        if (
            values.shape != limits.shape
            or not np.all(np.isfinite(values))
            or not np.all(np.isfinite(limits))
        ):
            raise ValueError(f"{name} gate values are malformed")
        observed = values.tolist()
        limit = limits.tolist()
        passed = bool(np.all(values < limits))
    elif comparison == "all_strictly_greater_than":
        values = np.asarray(value, dtype=float)
        if values.size < 1 or not np.all(np.isfinite(values)):
            raise ValueError(f"{name} gate values are malformed")
        observed = values.tolist()
        passed = bool(np.all(values > float(limit)))
    else:
        raise ValueError(f"unsupported scientific comparison: {comparison}")
    gates.append(
        {
            "name": name,
            "value": observed,
            "limit": limit,
            "comparison": comparison,
            "passed": passed,
        }
    )


def _complete_spectrum_metrics(name: str, modes: Any) -> dict[str, Any]:
    thresholds = HALF_STEP_CONTRACT["thresholds"]
    eigenvalues = np.asarray(modes.eigenvalues, dtype=float)
    if (
        eigenvalues.ndim != 1
        or eigenvalues.size < 2
        or not np.all(np.isfinite(eigenvalues))
    ):
        raise ValueError(f"{name} projected spectrum is invalid")
    frequencies = hessian_eigenvalues_to_wavenumbers_cm(eigenvalues)
    negative = eigenvalues < thresholds["negative_eigenvalue_threshold"]
    near_zero = (
        np.abs(eigenvalues)
        <= thresholds["near_zero_eigenvalue_inclusive_maximum_absolute"]
    )
    positive = (
        eigenvalues > thresholds["near_zero_eigenvalue_inclusive_maximum_absolute"]
    )
    negative_indices = np.flatnonzero(negative)
    return {
        "eigenvalues_hartree_per_bohr2_amu": eigenvalues.tolist(),
        "signed_wavenumbers_cm": frequencies.tolist(),
        "negative_mode_indices": negative_indices.tolist(),
        "negative_mode_count_below_negative_1e_minus_8": int(
            np.count_nonzero(negative)
        ),
        "near_zero_mode_indices": np.flatnonzero(near_zero).tolist(),
        "near_zero_mode_count_at_1e_minus_8": int(np.count_nonzero(near_zero)),
        "lowest_positive_eigenvalue": (
            float(np.min(eigenvalues[positive])) if np.any(positive) else None
        ),
        "imaginary_wavenumber_cm": (
            float(abs(frequencies[negative_indices[0]]))
            if negative_indices.size
            else None
        ),
    }


def _complete_half_step_analysis(
    transition_state: Any,
    masses: np.ndarray,
    route_vector: np.ndarray,
    h_inner: np.ndarray,
    h_outer: np.ndarray,
    raw_fd: np.ndarray,
    analytic: np.ndarray,
    *,
    point_gate_failures: list[dict[str, Any]],
) -> dict[str, Any]:
    """Persist every computable metric and reduce all unchanged gates at the end."""

    thresholds = HALF_STEP_CONTRACT["thresholds"]
    symmetric = {
        "H_0.005": 0.5 * (h_inner + h_inner.T),
        "H_0.01": 0.5 * (h_outer + h_outer.T),
        "FD": 0.5 * (raw_fd + raw_fd.T),
    }
    symmetry = _fd_symmetry_metrics(raw_fd)
    richardson_delta = _matrix_delta_metrics(symmetric["FD"], symmetric["H_0.005"])
    analytic_delta = _matrix_delta_metrics(symmetric["FD"], analytic)
    modes = {
        name: project_vibrational_hessian(transition_state.coords, masses, matrix)
        for name, matrix in {**symmetric, "analytic_reference": analytic}.items()
    }
    spectra = {
        name: _complete_spectrum_metrics(name, value) for name, value in modes.items()
    }
    reference_eigenvalues = np.asarray(modes["analytic_reference"].eigenvalues)
    eigenvalue_limits = np.maximum(
        thresholds["per_mode_reference_eigenvalue_delta_absolute_floor"],
        thresholds["per_mode_reference_eigenvalue_delta_relative_fraction"]
        * np.abs(reference_eigenvalues),
    )
    comparisons_to_reference = {}
    for name in ("H_0.005", "H_0.01", "FD"):
        comparison = _mode_comparison(modes[name], modes["analytic_reference"])
        comparisons_to_reference[name] = {
            **comparison,
            "eigenvalue_absolute_deltas": np.abs(
                np.asarray(modes[name].eigenvalues) - reference_eigenvalues
            ).tolist(),
            "eigenvalue_exclusive_limits": eigenvalue_limits.tolist(),
            "frequency_absolute_deltas_cm": np.abs(
                np.asarray(spectra[name]["signed_wavenumbers_cm"])
                - np.asarray(spectra["analytic_reference"]["signed_wavenumbers_cm"])
            ).tolist(),
        }
    inner_fd = {
        **_mode_comparison(modes["H_0.005"], modes["FD"]),
        "frequency_absolute_deltas_cm": np.abs(
            np.asarray(spectra["H_0.005"]["signed_wavenumbers_cm"])
            - np.asarray(spectra["FD"]["signed_wavenumbers_cm"])
        ).tolist(),
    }
    outer_fd = {
        **_mode_comparison(modes["H_0.01"], modes["FD"]),
        "frequency_absolute_deltas_cm": np.abs(
            np.asarray(spectra["H_0.01"]["signed_wavenumbers_cm"])
            - np.asarray(spectra["FD"]["signed_wavenumbers_cm"])
        ).tolist(),
    }
    normalized_route = np.asarray(route_vector, dtype=float).reshape(-1)
    route_norm = float(np.linalg.norm(normalized_route))
    if not math.isfinite(route_norm) or route_norm <= 1.0e-12:
        raise ValueError("mapped reaction vector is invalid")
    normalized_route /= route_norm
    reaction_overlaps = {
        name: float(
            abs(
                np.dot(
                    modes[name].mass_weighted_eigenvectors[0].reshape(-1),
                    normalized_route,
                )
            )
        )
        for name in ("H_0.005", "H_0.01", "FD", "analytic_reference")
    }
    mapped_reaction_overlap = {
        **reaction_overlaps,
        "fd_vs_analytic_reference_absolute_delta": abs(
            reaction_overlaps["FD"] - reaction_overlaps["analytic_reference"]
        ),
    }

    gates = [dict(item) for item in point_gate_failures]
    scalar_gates = (
        (
            "raw_fd_maximum_asymmetry",
            symmetry["maximum_absolute_asymmetry"],
            "raw_fd_maximum_absolute_asymmetry_exclusive_maximum",
        ),
        (
            "raw_fd_spectral_relative_asymmetry",
            symmetry["spectral_relative_asymmetry"],
            "raw_fd_spectral_relative_asymmetry_exclusive_maximum",
        ),
        (
            "symmetric_richardson_vs_h_inner_maximum_absolute_delta",
            richardson_delta["maximum_absolute_delta"],
            "symmetric_richardson_vs_h_h_maximum_absolute_delta_exclusive_maximum",
        ),
        (
            "symmetric_richardson_vs_h_inner_spectral_relative_delta",
            richardson_delta["spectral_relative_delta"],
            "symmetric_richardson_vs_h_h_spectral_relative_delta_exclusive_maximum",
        ),
        (
            "fd_vs_analytic_maximum_absolute_delta",
            analytic_delta["maximum_absolute_delta"],
            "fd_vs_analytic_symmetric_maximum_absolute_delta_exclusive_maximum",
        ),
        (
            "fd_vs_analytic_spectral_relative_delta",
            analytic_delta["spectral_relative_delta"],
            "fd_vs_analytic_symmetric_spectral_relative_delta_exclusive_maximum",
        ),
    )
    for gate_name, value, threshold_name in scalar_gates:
        _analysis_gate(
            gates,
            gate_name,
            value,
            thresholds[threshold_name],
            "strictly_less_than",
        )
    for name in ("H_0.005", "H_0.01", "FD"):
        spectrum = spectra[name]
        _analysis_gate(
            gates,
            f"{name}_negative_mode_count",
            spectrum["negative_mode_count_below_negative_1e_minus_8"],
            1,
            "equals",
        )
        _analysis_gate(
            gates,
            f"{name}_near_zero_mode_count",
            spectrum["near_zero_mode_count_at_1e_minus_8"],
            0,
            "equals",
        )
        _analysis_gate(
            gates,
            f"{name}_lowest_positive_eigenvalue",
            spectrum["lowest_positive_eigenvalue"],
            thresholds["lowest_positive_eigenvalue_exclusive_minimum"],
            "strictly_greater_than",
        )
        _analysis_gate(
            gates,
            f"{name}_imaginary_wavenumber_cm",
            spectrum["imaginary_wavenumber_cm"],
            thresholds["imaginary_wavenumber_cm_exclusive_minimum"],
            "strictly_greater_than",
        )
        comparison = comparisons_to_reference[name]
        _analysis_gate(
            gates,
            f"{name}_reference_assignment_order",
            comparison["assignment_equals_eigenvalue_order"],
            True,
            "equals",
        )
        _analysis_gate(
            gates,
            f"{name}_reference_eigenvalue_deltas",
            comparison["eigenvalue_absolute_deltas"],
            comparison["eigenvalue_exclusive_limits"],
            "all_strictly_less_than",
        )
    fd_reference = comparisons_to_reference["FD"]
    _analysis_gate(
        gates,
        "fd_reference_all_ordered_overlaps",
        fd_reference["ordered_overlaps"],
        thresholds["fd_vs_reference_all_ordered_overlaps_exclusive_minimum"],
        "all_strictly_greater_than",
    )
    for index, label, overlap_threshold, frequency_threshold in (
        (
            0,
            "unstable",
            "fd_vs_reference_unstable_overlap_exclusive_minimum",
            "fd_vs_reference_imaginary_frequency_delta_cm_exclusive_maximum",
        ),
        (
            1,
            "low_positive",
            "fd_vs_reference_low_positive_overlap_exclusive_minimum",
            "fd_vs_reference_low_positive_frequency_delta_cm_exclusive_maximum",
        ),
    ):
        _analysis_gate(
            gates,
            f"fd_reference_{label}_overlap",
            fd_reference["ordered_overlaps"][index],
            thresholds[overlap_threshold],
            "strictly_greater_than",
        )
        _analysis_gate(
            gates,
            f"fd_reference_{label}_frequency_delta",
            fd_reference["frequency_absolute_deltas_cm"][index],
            thresholds[frequency_threshold],
            "strictly_less_than",
        )
    _analysis_gate(
        gates,
        "h_inner_fd_assignment_order",
        inner_fd["assignment_equals_eigenvalue_order"],
        True,
        "equals",
    )
    _analysis_gate(
        gates,
        "h_outer_fd_assignment_order",
        outer_fd["assignment_equals_eigenvalue_order"],
        True,
        "equals",
    )
    for index, label, overlap_threshold, frequency_threshold in (
        (
            0,
            "unstable",
            "h_h_vs_fd_unstable_overlap_exclusive_minimum",
            "h_h_vs_fd_imaginary_frequency_delta_cm_exclusive_maximum",
        ),
        (
            1,
            "low_positive",
            "h_h_vs_fd_low_positive_overlap_exclusive_minimum",
            "h_h_vs_fd_low_positive_frequency_delta_cm_exclusive_maximum",
        ),
    ):
        _analysis_gate(
            gates,
            f"h_inner_fd_{label}_overlap",
            inner_fd["ordered_overlaps"][index],
            thresholds[overlap_threshold],
            "strictly_greater_than",
        )
        _analysis_gate(
            gates,
            f"h_inner_fd_{label}_frequency_delta",
            inner_fd["frequency_absolute_deltas_cm"][index],
            thresholds[frequency_threshold],
            "strictly_less_than",
        )
    _analysis_gate(
        gates,
        "fd_mapped_reaction_overlap",
        reaction_overlaps["FD"],
        thresholds["fd_mapped_reaction_overlap_exclusive_minimum"],
        "strictly_greater_than",
    )
    _analysis_gate(
        gates,
        "fd_mapped_reaction_overlap_vs_reference_delta",
        mapped_reaction_overlap["fd_vs_analytic_reference_absolute_delta"],
        thresholds["fd_mapped_reaction_overlap_vs_reference_exclusive_maximum_delta"],
        "strictly_less_than",
    )
    return {
        "confirmation_passed": bool(all(gate.get("passed") is True for gate in gates)),
        "gates": gates,
        "raw_fd_symmetry": symmetry,
        "symmetric_richardson_vs_h_inner": richardson_delta,
        "fd_vs_analytic_symmetric": analytic_delta,
        "spectra": spectra,
        "comparisons_to_analytic_reference": comparisons_to_reference,
        "h_inner_vs_fd_modes": inner_fd,
        "h_outer_vs_fd_modes": outer_fd,
        "mapped_reaction_overlap": mapped_reaction_overlap,
    }


def _initial_fd_status(
    identity: dict[str, Any], preflight_root: Path, reference_root: Path
) -> dict[str, Any]:
    contract_sha = campaign._canonical_hash(FINITE_DIFFERENCE_CONTRACT)
    return {
        "schema": FINITE_DIFFERENCE_SCHEMA,
        "state": "running",
        "route": ROUTE,
        **identity,
        "preflight_root": str(preflight_root),
        "reference_root": str(reference_root),
        "contract": FINITE_DIFFERENCE_CONTRACT,
        "contract_sha256": contract_sha,
        "preflight_receipt_sha256": None,
        "reference_receipt_sha256": None,
        "reference_matrix_sha256": None,
        "reference_binding": None,
        "reference_binding_sha256": None,
        "completed_points": [],
        "current_point": None,
        "confirmation_passed": False,
    }


def _point_bindings(status: dict[str, Any]) -> dict[str, Any]:
    return {
        "git_sha": status["git_sha"],
        "script_sha256": status["script_sha256"],
        "diagnostic_execution_identity_sha256": status[
            "diagnostic_execution_identity_sha256"
        ],
        "contract_sha256": status["contract_sha256"],
        "preflight_receipt_sha256": status["preflight_receipt_sha256"],
        "reference_receipt_sha256": status["reference_receipt_sha256"],
        "reference_matrix_sha256": status["reference_matrix_sha256"],
        "reference_binding_sha256": status["reference_binding_sha256"],
    }


_POINT_OBSERVATION_KEYS = {
    "backend",
    "electronic_hartree",
    "scf_converged",
    "grid_point_count",
    "spin_2s",
    "spin_square",
    "gradient_grid_response",
    "density_initial_guess_sha256",
    "density_final_sha256",
    "geometry_fingerprint",
    "elapsed_seconds",
}
_POINT_BINDING_KEYS = {
    "git_sha",
    "script_sha256",
    "diagnostic_execution_identity_sha256",
    "contract_sha256",
    "preflight_receipt_sha256",
    "reference_receipt_sha256",
    "reference_matrix_sha256",
    "reference_binding_sha256",
}
_CENTER_REFERENCE_KEYS = {
    "reference_electronic_hartree",
    "electronic_energy_absolute_delta_hartree",
    "reference_s2",
    "s2_absolute_delta",
}
_DENSITY_CONTINUITY_KEYS = {
    "final_vs_center_normalized_frobenius_delta",
    "final_frobenius_norm",
    "center_frobenius_norm",
}


def _load_half_step_source(prior_root: Path) -> dict[str, Any]:
    """Validate the immutable 73-point source chain before any calculator call."""

    root = campaign._safe_absolute_root(prior_root)
    if str(root) != PRIOR_FD_ROOT:
        raise ValueError("half-step extension requires the exact prior FD receipt root")
    receipt, receipt_raw = campaign._read_json_object(
        root / "receipt.json", label="prior FD receipt", max_bytes=4 * 1024 * 1024
    )
    receipt_sha = hashlib.sha256(receipt_raw).hexdigest()
    if receipt_sha != PRIOR_FD_RECEIPT_SHA256:
        raise ValueError("prior FD receipt SHA-256 mismatch")
    if (
        receipt.get("git_sha") != PRIOR_FD_GIT_SHA
        or receipt.get("script_sha256") != PRIOR_FD_SCRIPT_SHA256
        or receipt.get("preflight_root") != PRIOR_FD_PREFLIGHT_ROOT
        or receipt.get("preflight_receipt_sha256") != PRIOR_FD_PREFLIGHT_SHA256
        or receipt.get("reference_root") != PRIOR_FD_REFERENCE_ROOT
        or receipt.get("reference_receipt_sha256") != REFERENCE_RECEIPT_SHA256
        or receipt.get("reference_matrix_sha256") != REFERENCE_MATRIX_SHA256
    ):
        raise ValueError(
            "prior FD route, method, preflight, or reference identity is invalid"
        )
    repository = Path(__file__).resolve().parents[2]
    relative_source = Path(__file__).resolve().relative_to(repository)
    committed_source = subprocess.run(
        ["git", "show", f"{PRIOR_FD_GIT_SHA}:{relative_source.as_posix()}"],
        cwd=repository,
        check=True,
        capture_output=True,
    ).stdout
    if hashlib.sha256(committed_source).hexdigest() != PRIOR_FD_SCRIPT_SHA256:
        raise RuntimeError("prior FD Git source does not match its pinned SHA-256")
    source_files = receipt.get("diagnostic_execution_identity", {}).get("source_files")
    if (
        type(source_files) is not dict
        or source_files.get(str(Path(__file__).resolve())) != PRIOR_FD_SCRIPT_SHA256
    ):
        raise ValueError("prior FD execution identity does not bind its script source")

    validated_receipt, validated_sha = _validated_completed_fd_receipt(
        root, verify_resident_identity=False
    )
    if validated_sha != receipt_sha or validated_receipt != receipt:
        raise ValueError("prior FD complete receipt chain is inconsistent")
    terminal, terminal_raw = campaign._read_json_object(
        root / TERMINAL, label="prior FD terminal", max_bytes=1024 * 1024
    )
    if hashlib.sha256(terminal_raw).hexdigest() != PRIOR_FD_TERMINAL_SHA256:
        raise ValueError("prior FD terminal SHA-256 mismatch")
    if terminal != _validated_terminal(root / TERMINAL) or (
        terminal.get("schema") != FINITE_DIFFERENCE_SCHEMA
        or terminal.get("state") != "completed"
        or terminal.get("accepted_campaign_result") is not False
        or terminal.get("confirmation_passed") is not False
        or terminal.get("receipt_sha256") != receipt_sha
        or terminal.get("completed_point_count") != 73
        or terminal.get("git_sha") != PRIOR_FD_GIT_SHA
        or terminal.get("script_sha256") != PRIOR_FD_SCRIPT_SHA256
        or terminal.get("preflight_receipt_sha256") != PRIOR_FD_PREFLIGHT_SHA256
        or terminal.get("reference_receipt_sha256") != REFERENCE_RECEIPT_SHA256
        or terminal.get("reference_matrix_sha256") != REFERENCE_MATRIX_SHA256
        or terminal.get("contract_sha256")
        != campaign._canonical_hash(FINITE_DIFFERENCE_CONTRACT)
    ):
        raise ValueError("prior FD terminal does not bind the completed receipt")
    _status, status_raw = campaign._read_json_object(
        root / "status.json", label="prior FD status", max_bytes=4 * 1024 * 1024
    )
    if hashlib.sha256(status_raw).hexdigest() != PRIOR_FD_STATUS_SHA256:
        raise ValueError("prior FD status SHA-256 mismatch")

    preflight_root = campaign._safe_absolute_root(Path(PRIOR_FD_PREFLIGHT_ROOT))
    preflight, preflight_sha = campaign._validated_preflight(preflight_root, ROUTE)
    if preflight_sha != PRIOR_FD_PREFLIGHT_SHA256 or preflight.get(
        "identity"
    ) != receipt.get("campaign_identity"):
        raise ValueError("prior FD preflight identity mismatch")
    base_settings, _settings_receipt = campaign._canonical_dft_settings(preflight)
    settings = _fd_settings(base_settings)
    template = reactions(gpu=True, basis="def2-svp")[ROUTE].cluster
    fingerprints = campaign._trusted_input_fingerprints_from_route_record(
        preflight["routes"][ROUTE]
    )
    inputs = campaign._trusted_route_input_snapshots(
        campaign.DEFAULT_BUNDLE_ROOT, ROUTE, template, fingerprints
    )
    transition_state = inputs["transition_state"]
    atom_count = len(transition_state.symbols)
    if atom_count != 6:
        raise ValueError("prior FD source is not the exact six-atom frozen TS")
    masses = np.asarray(preflight["routes"][ROUTE]["masses_amu"], dtype=float)
    (
        reference_receipt,
        reference_matrix,
        reference_case,
        reference_binding,
    ) = _load_reference(
        Path(PRIOR_FD_REFERENCE_ROOT),
        atom_count,
        preflight=preflight,
        transition_state=transition_state,
        reactant=inputs["reactant"],
        product=inputs["product"],
        masses=masses,
    )
    campaign._strict_json_equal(
        receipt.get("reference_binding"),
        reference_binding,
        label="prior FD reference binding",
    )

    ancestry = receipt["point_receipts"]
    ancestry_by_key = {item["point_key"]: item for item in ancestry}
    center_record, center_sha, center_arrays = _read_record_bundle(
        root,
        Path("points/center"),
        expected_receipt_keys=(
            _POINT_OBSERVATION_KEYS
            | _POINT_BINDING_KEYS
            | _CENTER_REFERENCE_KEYS
            | {
                "schema",
                "point_key",
                "kind",
                "artifacts",
                "physical_fmax_ev_per_angstrom",
            }
        ),
    )
    center_density = center_arrays["density.f64"]
    center_density_sha = hashlib.sha256(
        np.ascontiguousarray(center_density, dtype="<f8").tobytes()
    ).hexdigest()
    if center_sha != ancestry_by_key["center"]["receipt_sha256"]:
        raise ValueError("prior FD center receipt chain mismatch")

    reused_gradients: dict[tuple[int, int], np.ndarray] = {}
    reused_points = []
    reused_point_data = []
    for definition in _displacement_plan(atom_count):
        multiplier = definition["step_multiplier"]
        if multiplier not in {-1, 1}:
            continue
        key = definition["key"]
        point, point_sha, arrays = _read_record_bundle(
            root,
            Path("points") / key,
            expected_receipt_keys=(
                _POINT_OBSERVATION_KEYS
                | _POINT_BINDING_KEYS
                | _DENSITY_CONTINUITY_KEYS
                | set(definition)
                | {"schema", "point_key", "kind", "artifacts"}
            ),
            expected_artifacts={
                "gradient.f64": ((atom_count, 3), "hartree / bohr"),
                "density.f64": (center_density.shape, "electrons"),
            },
        )
        if point_sha != ancestry_by_key[key]["receipt_sha256"]:
            raise ValueError(f"prior FD reused point receipt chain mismatch: {key}")
        reused_gradients[(definition["coordinate_index"], multiplier)] = arrays[
            "gradient.f64"
        ].reshape(-1)
        reused_points.append(
            {
                "point_key": key,
                "receipt_sha256": point_sha,
                "gradient_sha256": point["artifacts"]["gradient.f64"]["sha256"],
            }
        )
        reused_point_data.append((point, arrays["density.f64"]))
    if len(reused_points) != 36:
        raise ValueError("prior FD reused +/-0.01 Bohr payload chain is incomplete")
    return {
        "prior_receipt": receipt,
        "prior_receipt_sha256": receipt_sha,
        "prior_terminal_sha256": PRIOR_FD_TERMINAL_SHA256,
        "preflight": preflight,
        "preflight_root": preflight_root,
        "preflight_receipt_sha256": preflight_sha,
        "reference_root": Path(PRIOR_FD_REFERENCE_ROOT),
        "reference_receipt": reference_receipt,
        "reference_binding": reference_binding,
        "reference_matrix": reference_matrix,
        "reference_case": reference_case,
        "transition_state": transition_state,
        "reactant": inputs["reactant"],
        "product": inputs["product"],
        "masses": masses,
        "settings": settings,
        "center_record": center_record,
        "center_density": center_density,
        "center_receipt": {
            "point_key": "center",
            "receipt_sha256": center_sha,
            "density_sha256": center_density_sha,
        },
        "reused_points": reused_points,
        "reused_point_data": reused_point_data,
        "reused_gradients": reused_gradients,
    }


def _resume_fd_points(
    output_root: Path,
    status: dict[str, Any],
    plan: tuple[dict[str, Any], ...],
    *,
    claim: _OutputRootClaim,
) -> dict[str, tuple[dict[str, Any], dict[str, np.ndarray], str]]:
    expected_keys = ["center", *(item["key"] for item in plan)]
    points_root = output_root / "points"
    existing_names = set()
    if points_root.exists() or points_root.is_symlink():
        if points_root.is_symlink() or not points_root.is_dir():
            raise ValueError(
                "finite-difference points root must be a regular directory"
            )
        for child in points_root.iterdir():
            if (
                child.is_symlink()
                or not child.is_dir()
                or child.name not in expected_keys
            ):
                raise ValueError(
                    "finite-difference points root contains an unexpected artifact"
                )
            existing_names.add(child.name)
    prefix = []
    for key in expected_keys:
        if key in existing_names:
            prefix.append(key)
        else:
            break
    if existing_names != set(prefix):
        raise ValueError(
            "finite-difference point artifacts are not a contiguous prefix"
        )
    current_point = status.get("current_point")
    if current_point is not None and current_point not in expected_keys:
        raise ValueError("finite-difference status names an unexpected current point")
    if current_point is not None and current_point not in existing_names:
        raise RuntimeError(
            "finite-difference current point has an uncertain unreceipted attempt; "
            "refusing a duplicate evaluation"
        )
    loaded = {}
    disk_progress = []
    definitions = {item["key"]: item for item in plan}
    atom_count = len(plan) // 12
    center_density_shape: tuple[int, ...] | None = None
    for key in prefix:
        receipt_keys = (
            _POINT_OBSERVATION_KEYS
            | _POINT_BINDING_KEYS
            | {"schema", "point_key", "kind", "artifacts"}
        )
        artifact_contract: dict[str, tuple[tuple[int, ...], str]] = {
            "gradient.f64": ((atom_count, 3), "hartree / bohr"),
        }
        if key == "center":
            receipt_keys |= _CENTER_REFERENCE_KEYS | {"physical_fmax_ev_per_angstrom"}
        else:
            receipt_keys |= set(definitions[key]) | _DENSITY_CONTINUITY_KEYS
        receipt, receipt_sha, arrays = _read_record_bundle(
            output_root,
            Path("points") / key,
            expected_receipt_keys=receipt_keys,
        )
        if (
            receipt.get("schema") != FINITE_DIFFERENCE_SCHEMA
            or receipt.get("point_key") != key
        ):
            raise ValueError(
                f"finite-difference point receipt identity mismatch: {key}"
            )
        if any(
            receipt.get(field) != value
            for field, value in _point_bindings(status).items()
        ):
            raise ValueError(f"finite-difference point ancestry mismatch: {key}")
        if key == "center":
            if receipt.get("kind") != "center" or set(arrays) != {
                "gradient.f64",
                "density.f64",
            }:
                raise ValueError("finite-difference center artifacts are incomplete")
            density = arrays["density.f64"]
            if (
                density.ndim != 3
                or density.shape[0] != 2
                or density.shape[1] != density.shape[2]
            ):
                raise ValueError("finite-difference center density shape is invalid")
            center_density_shape = density.shape
            artifact_contract["density.f64"] = (density.shape, "electrons")
        else:
            definition = definitions[key]
            if receipt.get("kind") != "displacement" or any(
                receipt.get(field) != value for field, value in definition.items()
            ):
                raise ValueError(
                    f"finite-difference displacement receipt mismatch: {key}"
                )
            if set(arrays) != {"gradient.f64", "density.f64"}:
                raise ValueError(
                    f"finite-difference point artifacts are incomplete: {key}"
                )
            if (
                center_density_shape is None
                or arrays["density.f64"].shape != center_density_shape
            ):
                raise ValueError(
                    f"finite-difference final density shape mismatch: {key}"
                )
            artifact_contract["density.f64"] = (center_density_shape, "electrons")
        if set(receipt["artifacts"]) != set(artifact_contract):
            raise ValueError(
                f"finite-difference point artifact inventory mismatch: {key}"
            )
        for filename, (shape, units) in artifact_contract.items():
            artifact = receipt["artifacts"][filename]
            if artifact["shape"] != list(shape) or artifact["units"] != units:
                raise ValueError(
                    "finite-difference point artifact schema mismatch: "
                    f"{key}/{filename}"
                )
        loaded[key] = (receipt, arrays, receipt_sha)
        disk_progress.append({"point_key": key, "receipt_sha256": receipt_sha})
    declared = status.get("completed_points")
    if type(declared) is not list or declared != disk_progress[: len(declared)]:
        raise ValueError(
            "finite-difference status progress diverges from durable artifacts"
        )
    if len(declared) > len(disk_progress):
        raise ValueError("finite-difference status claims a missing durable point")
    if declared != disk_progress:
        status["completed_points"] = disk_progress
        _claimed_write(claim, output_root / "status.json", _json_bytes(status))
    return loaded


def _run_finite_difference_locked(
    preflight_root: Path,
    reference_root: Path,
    output_root: Path,
    claim: _OutputRootClaim,
) -> dict[str, Any]:
    status_path = output_root / "status.json"
    claim.verify()
    legacy_lock = output_root / ".preflight.lock"
    if legacy_lock.exists() or legacy_lock.is_symlink():
        claim.verify()
        legacy_descriptor = os.open(
            legacy_lock,
            os.O_RDWR | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_CLOEXEC", 0),
        )
        legacy_locked = False
        try:
            legacy_status = os.fstat(legacy_descriptor)
            if not stat.S_ISREG(legacy_status.st_mode):
                raise ValueError("legacy output-local claim is not a regular file")
            try:
                fcntl.flock(legacy_descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as exc:
                raise ConcurrentRunError(
                    "legacy output-local claim is still held; refusing a duplicate run"
                ) from exc
            legacy_locked = True
            path_status = os.stat(legacy_lock, follow_symlinks=False)
            if (path_status.st_dev, path_status.st_ino) != (
                legacy_status.st_dev,
                legacy_status.st_ino,
            ):
                raise RuntimeError("legacy output-local claim was replaced")
            legacy_lock.unlink()
        finally:
            if legacy_locked:
                fcntl.flock(legacy_descriptor, fcntl.LOCK_UN)
            os.close(legacy_descriptor)
        claim.verify()
    identity = _fd_execution_identity()
    if status_path.exists() or status_path.is_symlink():
        status, _ = campaign._read_json_object(
            status_path, label="FD confirmation status", max_bytes=1024 * 1024
        )
        expected = _initial_fd_status(identity, preflight_root, reference_root)
        if set(status) != set(expected):
            raise ValueError("finite-difference resume status schema is not exact")
        for key in (
            "schema",
            "route",
            "git_sha",
            "script_sha256",
            "diagnostic_execution_identity",
            "diagnostic_execution_identity_sha256",
            "preflight_root",
            "reference_root",
            "contract",
            "contract_sha256",
        ):
            if status.get(key) != expected[key]:
                raise ValueError(f"finite-difference resume status {key} mismatch")
        if (
            status.get("state") != "running"
            or (output_root / TERMINAL).exists()
            or (output_root / "receipt.json").exists()
        ):
            raise ValueError(
                "finite-difference output root is not resumable running state"
            )
    else:
        unexpected = {child.name for child in output_root.iterdir()}
        if unexpected:
            raise FileExistsError(
                "finite-difference output root contains unowned artifacts"
            )
        status = _initial_fd_status(identity, preflight_root, reference_root)
        _claimed_write(claim, status_path, _json_bytes(status))
    try:
        preflight, preflight_sha = campaign._validate_production_boundary(
            preflight_root, ROUTE
        )
        base_settings, _ = campaign._canonical_dft_settings(preflight)
        settings = _fd_settings(base_settings)
        template = reactions(gpu=True, basis="def2-svp")[ROUTE].cluster
        fingerprints = campaign._trusted_input_fingerprints_from_route_record(
            preflight["routes"][ROUTE]
        )
        inputs = campaign._trusted_route_input_snapshots(
            campaign.DEFAULT_BUNDLE_ROOT, ROUTE, template, fingerprints
        )
        transition_state = inputs["transition_state"]
        atom_count = len(transition_state.symbols)
        if atom_count != 6:
            raise ValueError(
                "finite-difference confirmation requires the exact six-atom frozen TS"
            )
        plan = _displacement_plan(atom_count)
        if (
            len(plan)
            != FINITE_DIFFERENCE_CONTRACT["stencil"]["displaced_gradient_count"]
        ):
            raise RuntimeError("finite-difference displacement count contract drifted")
        masses = np.asarray(preflight["routes"][ROUTE]["masses_amu"], dtype=float)
        if (
            masses.shape != (atom_count,)
            or not np.all(np.isfinite(masses))
            or np.any(masses <= 0.0)
        ):
            raise ValueError(
                "finite-difference masses are missing, non-finite, or non-positive"
            )
        (
            reference_receipt,
            analytic_matrix,
            reference_case,
            reference_binding,
        ) = _load_reference(
            reference_root,
            atom_count,
            preflight=preflight,
            transition_state=transition_state,
            reactant=inputs["reactant"],
            product=inputs["product"],
            masses=masses,
        )
        status.update(
            {
                "preflight_receipt_sha256": preflight_sha,
                "reference_receipt_sha256": REFERENCE_RECEIPT_SHA256,
                "reference_matrix_sha256": REFERENCE_MATRIX_SHA256,
                "reference_binding": reference_binding,
                "reference_binding_sha256": campaign._canonical_hash(reference_binding),
            }
        )
        _claimed_write(claim, status_path, _json_bytes(status))
        loaded = _resume_fd_points(output_root, status, plan, claim=claim)
        if "center" not in loaded:
            status["current_point"] = "center"
            _claimed_write(claim, status_path, _json_bytes(status))
            claim.verify()
            point, gradient, density = _evaluate_fd_gradient(
                transition_state, settings, central_density=None
            )
            claim.verify()
            point.update(
                {
                    "schema": FINITE_DIFFERENCE_SCHEMA,
                    "point_key": "center",
                    "kind": "center",
                    **_point_bindings(status),
                    "physical_fmax_ev_per_angstrom": float(
                        np.max(np.linalg.norm(gradient, axis=1))
                    )
                    * HARTREE_TO_EV
                    / BOHR_TO_ANGSTROM,
                    **_center_reference_metrics(point, reference_case),
                }
            )
            point, point_sha = _publish_record_bundle(
                output_root,
                Path("points/center"),
                point,
                {
                    "gradient.f64": (gradient, "hartree / bohr"),
                    "density.f64": (density, "electrons"),
                },
                claim=claim,
            )
            _validate_point_record(
                point,
                center_spin_square=None,
                center=True,
                enforce_scientific_gates=False,
            )
            status["completed_points"].append(
                {"point_key": "center", "receipt_sha256": point_sha}
            )
            _claimed_write(claim, status_path, _json_bytes(status))
            loaded["center"] = (
                point,
                {"gradient.f64": gradient, "density.f64": density},
                point_sha,
            )
        center_record, center_arrays, _ = loaded["center"]
        if center_record.get("geometry_fingerprint") != frequency_geometry_fingerprint(
            transition_state
        ):
            raise ValueError("center gradient geometry differs from the frozen TS")
        center_gradient = center_arrays["gradient.f64"]
        if center_gradient.shape != (atom_count, 3):
            raise ValueError("center gradient artifact shape is invalid")
        observed_center_fmax = (
            float(np.max(np.linalg.norm(center_gradient, axis=1)))
            * HARTREE_TO_EV
            / BOHR_TO_ANGSTROM
        )
        if center_record.get("physical_fmax_ev_per_angstrom") != observed_center_fmax:
            raise ValueError("center fmax does not reproduce from its raw gradient")
        center_spin = _validate_point_record(
            center_record,
            center_spin_square=None,
            center=True,
            enforce_scientific_gates=False,
        )
        campaign._strict_json_equal(
            {key: center_record.get(key) for key in _CENTER_REFERENCE_KEYS},
            _center_reference_metrics(center_record, reference_case),
            label="center-to-reference continuity receipt",
        )
        central_density = center_arrays["density.f64"]
        density_sha = hashlib.sha256(
            np.ascontiguousarray(central_density, dtype="<f8").tobytes()
        ).hexdigest()
        if (
            center_record.get("density_initial_guess_sha256") != density_sha
            or center_record.get("density_final_sha256") != density_sha
        ):
            raise ValueError(
                "central density artifact does not match its bound SHA-256"
            )
        gradients = {}
        displaced_points: list[tuple[dict[str, Any], np.ndarray]] = []
        for definition in plan:
            key = definition["key"]
            if key not in loaded:
                status["current_point"] = key
                _claimed_write(claim, status_path, _json_bytes(status))
                displaced_coords = np.asarray(
                    transition_state.coords, dtype=float
                ).copy()
                displaced_coords.reshape(-1)[definition["coordinate_index"]] += (
                    definition["displacement_bohr"] * BOHR_TO_ANGSTROM
                )
                displaced = replace(transition_state, coords=displaced_coords)
                claim.verify()
                point, gradient, returned_density = _evaluate_fd_gradient(
                    displaced, settings, central_density=central_density
                )
                claim.verify()
                point.update(
                    {
                        "schema": FINITE_DIFFERENCE_SCHEMA,
                        "point_key": key,
                        "kind": "displacement",
                        **_point_bindings(status),
                        **definition,
                        **_density_continuity_metrics(
                            returned_density, central_density
                        ),
                    }
                )
                point, point_sha = _publish_record_bundle(
                    output_root,
                    Path("points") / key,
                    point,
                    {
                        "gradient.f64": (gradient, "hartree / bohr"),
                        "density.f64": (returned_density, "electrons"),
                    },
                    claim=claim,
                )
                _validate_point_record(
                    point,
                    center_spin_square=center_spin,
                    center=False,
                    enforce_scientific_gates=False,
                )
                status["completed_points"].append(
                    {"point_key": key, "receipt_sha256": point_sha}
                )
                _claimed_write(claim, status_path, _json_bytes(status))
                loaded[key] = (
                    point,
                    {"gradient.f64": gradient, "density.f64": returned_density},
                    point_sha,
                )
            point, arrays, _ = loaded[key]
            _validate_point_record(
                point,
                center_spin_square=center_spin,
                center=False,
                enforce_scientific_gates=False,
            )
            expected_coords = np.asarray(transition_state.coords, dtype=float).copy()
            expected_coords.reshape(-1)[definition["coordinate_index"]] += (
                definition["displacement_bohr"] * BOHR_TO_ANGSTROM
            )
            expected_displaced = replace(transition_state, coords=expected_coords)
            if point.get("geometry_fingerprint") != frequency_geometry_fingerprint(
                expected_displaced
            ):
                raise ValueError(f"displacement geometry fingerprint mismatch: {key}")
            if arrays["gradient.f64"].shape != (atom_count, 3):
                raise ValueError(
                    f"displacement gradient artifact shape is invalid: {key}"
                )
            if point.get("density_initial_guess_sha256") != density_sha:
                raise ValueError(
                    f"displacement did not use the bound central density: {key}"
                )
            final_density = arrays["density.f64"]
            final_density_sha = hashlib.sha256(
                np.ascontiguousarray(final_density, dtype="<f8").tobytes()
            ).hexdigest()
            if point.get("density_final_sha256") != final_density_sha:
                raise ValueError(f"displacement final density hash mismatch: {key}")
            campaign._strict_json_equal(
                {name: point.get(name) for name in _DENSITY_CONTINUITY_KEYS},
                _density_continuity_metrics(final_density, central_density),
                label=f"displacement density continuity receipt {key}",
            )
            displaced_points.append((point, final_density))
            gradients[
                (definition["coordinate_index"], definition["step_multiplier"])
            ] = arrays["gradient.f64"].reshape(-1)
        status["current_point"] = None
        _claimed_write(claim, status_path, _json_bytes(status))
        h_h, h_2h, raw_fd = _finite_difference_matrices(gradients, 3 * atom_count)
        route_vector = campaign._mapped_route_vector(
            transition_state,
            masses,
            inputs["reactant"].coords,
            inputs["product"].coords,
        )
        symmetric = {
            "H_h": 0.5 * (h_h + h_h.T),
            "H_2h": 0.5 * (h_2h + h_2h.T),
            "FD": 0.5 * (raw_fd + raw_fd.T),
        }
        confirmation_passed, analysis = _scientific_analysis_verdict(
            transition_state,
            masses,
            route_vector,
            center_record,
            central_density,
            displaced_points,
            h_h,
            h_2h,
            raw_fd,
            analytic_matrix,
        )
        aggregate_record = {
            "schema": FINITE_DIFFERENCE_SCHEMA,
            "kind": "finite-difference-Hessian-analysis",
            "confirmation_passed": confirmation_passed,
            **_point_bindings(status),
            "point_receipts": status["completed_points"],
            "analysis": analysis,
        }
        matrix_arrays = {
            "H_h-raw.f64": (h_h, "hartree / bohr^2"),
            "H_2h-raw.f64": (h_2h, "hartree / bohr^2"),
            "richardson-raw.f64": (raw_fd, "hartree / bohr^2"),
            "H_h-symmetric.f64": (symmetric["H_h"], "hartree / bohr^2"),
            "H_2h-symmetric.f64": (symmetric["H_2h"], "hartree / bohr^2"),
            "richardson-symmetric.f64": (symmetric["FD"], "hartree / bohr^2"),
        }
        matrices_path = output_root / "matrices"
        matrix_artifact_contract = {
            filename: ((3 * atom_count, 3 * atom_count), units)
            for filename, (_array, units) in matrix_arrays.items()
        }
        if matrices_path.exists() or matrices_path.is_symlink():
            aggregate, aggregate_sha, observed_arrays = _read_record_bundle(
                output_root,
                Path("matrices"),
                expected_receipt_keys=set(aggregate_record) | {"artifacts"},
                expected_artifacts=matrix_artifact_contract,
            )
            campaign._strict_json_equal(
                {field: aggregate[field] for field in aggregate_record},
                aggregate_record,
                label="resumed matrix-analysis receipt",
            )
            for filename, (expected_array, _units) in matrix_arrays.items():
                if not np.array_equal(observed_arrays[filename], expected_array):
                    raise ValueError(
                        "resumed matrix artifact diverges from recomputation: "
                        f"{filename}"
                    )
        else:
            aggregate, aggregate_sha = _publish_record_bundle(
                output_root,
                Path("matrices"),
                aggregate_record,
                matrix_arrays,
                claim=claim,
            )
        finished = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        receipt = {
            "schema": FINITE_DIFFERENCE_SCHEMA,
            "state": "completed",
            "accepted_campaign_result": False,
            "confirmation_passed": confirmation_passed,
            "purpose": (
                "route-specific numerical confirmation only; no TS or SCT qualification"
            ),
            "authority": FINITE_DIFFERENCE_CONTRACT["policy"][
                "passed_authority" if confirmation_passed else "rejected_authority"
            ],
            "route": ROUTE,
            "git_sha": identity["git_sha"],
            "script_sha256": identity["script_sha256"],
            "diagnostic_execution_identity": identity["diagnostic_execution_identity"],
            "diagnostic_execution_identity_sha256": identity[
                "diagnostic_execution_identity_sha256"
            ],
            "preflight_root": str(preflight_root),
            "preflight_receipt_sha256": preflight_sha,
            "campaign_identity": preflight["identity"],
            "reference_root": str(reference_root),
            "reference_receipt_sha256": REFERENCE_RECEIPT_SHA256,
            "reference_matrix_sha256": REFERENCE_MATRIX_SHA256,
            "reference_binding": reference_binding,
            "reference_binding_sha256": status["reference_binding_sha256"],
            "reference_git_sha": reference_receipt["git_sha"],
            "reference_script_sha256": reference_receipt["script_sha256"],
            "contract": FINITE_DIFFERENCE_CONTRACT,
            "contract_sha256": status["contract_sha256"],
            "point_receipts": status["completed_points"],
            "aggregate_receipt": {
                "path": "matrices/receipt.json",
                "sha256": aggregate_sha,
            },
            "analysis": aggregate["analysis"],
            "finished_utc": finished,
        }
        receipt_raw = _json_bytes(receipt)
        receipt_sha = hashlib.sha256(receipt_raw).hexdigest()
        _claimed_write(claim, output_root / "receipt.json", receipt_raw)
        status.update(
            {
                "state": "completed",
                "confirmation_passed": confirmation_passed,
                "receipt": "receipt.json",
                "receipt_sha256": receipt_sha,
                "finished_utc": finished,
            }
        )
        terminal = _terminal_payload(status, state="completed", detail=None)
        terminal["receipt"] = "receipt.json"
        terminal["receipt_sha256"] = receipt_sha
        _terminal_write_noreplace(claim, output_root / TERMINAL, terminal)
        _claimed_write(claim, status_path, _json_bytes(status))
        return receipt
    except BaseException as exc:
        terminal_path = output_root / TERMINAL
        if terminal_path.exists() or terminal_path.is_symlink():
            campaign._read_json_object(
                terminal_path,
                label="authoritative FD terminal receipt",
                max_bytes=1024 * 1024,
            )
            raise
        detail = f"{type(exc).__name__}: {exc}"
        status.update(
            {
                "state": "failed",
                "confirmation_passed": False,
                "error": detail,
                "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            }
        )
        _claimed_write(claim, status_path, _json_bytes(status))
        if not terminal_path.exists() and not terminal_path.is_symlink():
            _terminal_write_noreplace(
                claim,
                terminal_path,
                _terminal_payload(status, state="failed", detail=detail),
            )
        raise


def run_finite_difference(
    preflight_root: Path, reference_root: Path, output_root: Path
) -> dict[str, Any]:
    """Run or resume the receipt-bound 72-displacement confirmation."""

    preflight = campaign._safe_absolute_root(preflight_root)
    reference = campaign._safe_absolute_root(reference_root)
    output = campaign._safe_absolute_root(output_root)
    with _exclusive_output_claim(output) as claim:
        try:
            return _run_finite_difference_locked(preflight, reference, output, claim)
        except BaseException as exc:
            if isinstance(exc, ConcurrentRunError):
                raise
            terminal_path = output / TERMINAL
            if terminal_path.exists() or terminal_path.is_symlink():
                raise
            status_path = output / "status.json"
            status: dict[str, Any] = {
                "schema": FINITE_DIFFERENCE_SCHEMA,
                "state": "failed",
                "route": ROUTE,
                "preflight_root": str(preflight),
                "reference_root": str(reference),
                "contract": FINITE_DIFFERENCE_CONTRACT,
                "contract_sha256": campaign._canonical_hash(FINITE_DIFFERENCE_CONTRACT),
                "completed_points": [],
                "current_point": None,
                "confirmation_passed": False,
            }
            if status_path.is_file() and not status_path.is_symlink():
                try:
                    loaded, _ = campaign._read_json_object(
                        status_path,
                        label="failed FD confirmation status",
                        max_bytes=1024 * 1024,
                    )
                except ValueError:
                    pass
                else:
                    status = loaded
            detail = f"{type(exc).__name__}: {exc}"
            status.update(
                {
                    "state": "failed",
                    "confirmation_passed": False,
                    "error": detail,
                    "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                }
            )
            _claimed_write(claim, status_path, _json_bytes(status))
            _terminal_write_noreplace(
                claim,
                terminal_path,
                _terminal_payload(status, state="failed", detail=detail),
            )
            raise


def _initial_half_step_status(
    identity: dict[str, Any], prior_root: Path, source: dict[str, Any]
) -> dict[str, Any]:
    center_receipt = {
        **source["center_receipt"],
        "density_shape": list(np.asarray(source["center_density"]).shape),
    }
    return {
        "schema": HALF_STEP_SCHEMA,
        "state": "running",
        "route": ROUTE,
        **identity,
        "prior_fd_root": str(prior_root),
        "prior_fd_receipt_sha256": source["prior_receipt_sha256"],
        "prior_fd_terminal_sha256": source["prior_terminal_sha256"],
        "prior_fd_status_sha256": PRIOR_FD_STATUS_SHA256,
        "preflight_root": str(source["preflight_root"]),
        "preflight_receipt_sha256": source["preflight_receipt_sha256"],
        "reference_root": str(source["reference_root"]),
        "reference_receipt_sha256": REFERENCE_RECEIPT_SHA256,
        "reference_matrix_sha256": REFERENCE_MATRIX_SHA256,
        "reference_binding": source["reference_binding"],
        "reference_binding_sha256": campaign._canonical_hash(
            source["reference_binding"]
        ),
        "center_receipt": center_receipt,
        "reused_prior_point_receipts": source["reused_points"],
        "contract": HALF_STEP_CONTRACT,
        "contract_sha256": campaign._canonical_hash(HALF_STEP_CONTRACT),
        "completed_points": [],
        "current_point": None,
        "accepted_campaign_result": False,
        "confirmation_passed": False,
    }


def _validate_half_step_resume_status(
    status: dict[str, Any], expected: dict[str, Any]
) -> None:
    """Require the exact running half-step status before reusing any field."""

    if set(status) != set(expected):
        raise ValueError("half-step resume status schema is not exact")
    mutable = {"state", "completed_points", "current_point", "confirmation_passed"}
    for key, value in expected.items():
        if key not in mutable and status.get(key) != value:
            raise ValueError(f"half-step resume status {key} mismatch")
    if status.get("state") != "running":
        raise ValueError("half-step resume status state is not running")
    if status.get("confirmation_passed") is not False:
        raise ValueError("half-step resume status scientific verdict is invalid")
    completed = status.get("completed_points")
    if type(completed) is not list or any(
        type(point) is not dict
        or set(point) != {"point_key", "receipt_sha256"}
        or type(point.get("point_key")) is not str
        or type(point.get("receipt_sha256")) is not str
        for point in completed
    ):
        raise ValueError("half-step resume status progress is invalid")
    expected_keys = [item["key"] for item in _half_step_plan(6)]
    if len(completed) > len(expected_keys):
        raise ValueError("half-step resume status progress is invalid")
    for point, expected_key in zip(completed, expected_keys, strict=False):
        if point["point_key"] != expected_key:
            raise ValueError("half-step resume status progress is not an exact prefix")
        campaign._require_sha(
            point["receipt_sha256"],
            length=64,
            label="half-step resume point receipt SHA-256",
        )
    current = status.get("current_point")
    if current is not None and current not in expected_keys:
        raise ValueError("half-step resume status current point is invalid")


def _half_step_point_bindings(status: dict[str, Any]) -> dict[str, Any]:
    return {
        "git_sha": status["git_sha"],
        "script_sha256": status["script_sha256"],
        "diagnostic_execution_identity_sha256": status[
            "diagnostic_execution_identity_sha256"
        ],
        "contract_sha256": status["contract_sha256"],
        "prior_fd_receipt_sha256": status["prior_fd_receipt_sha256"],
        "prior_fd_terminal_sha256": status["prior_fd_terminal_sha256"],
        "preflight_receipt_sha256": status["preflight_receipt_sha256"],
        "reference_receipt_sha256": status["reference_receipt_sha256"],
        "reference_matrix_sha256": status["reference_matrix_sha256"],
        "reference_binding_sha256": status["reference_binding_sha256"],
        "center_receipt_sha256": status["center_receipt"]["receipt_sha256"],
        "center_density_sha256": status["center_receipt"]["density_sha256"],
    }


_HALF_STEP_POINT_BINDING_KEYS = {
    "git_sha",
    "script_sha256",
    "diagnostic_execution_identity_sha256",
    "contract_sha256",
    "prior_fd_receipt_sha256",
    "prior_fd_terminal_sha256",
    "preflight_receipt_sha256",
    "reference_receipt_sha256",
    "reference_matrix_sha256",
    "reference_binding_sha256",
    "center_receipt_sha256",
    "center_density_sha256",
}


def _resume_half_step_points(
    output_root: Path,
    status: dict[str, Any],
    plan: tuple[dict[str, Any], ...],
    *,
    claim: _OutputRootClaim,
) -> dict[str, tuple[dict[str, Any], dict[str, np.ndarray], str]]:
    expected_keys = [item["key"] for item in plan]
    definitions = {item["key"]: item for item in plan}
    points_root = output_root / "points"
    existing_names: set[str] = set()
    if points_root.exists() or points_root.is_symlink():
        if points_root.is_symlink() or not points_root.is_dir():
            raise ValueError("half-step points root must be a regular directory")
        for child in points_root.iterdir():
            if (
                child.is_symlink()
                or not child.is_dir()
                or child.name not in expected_keys
            ):
                raise ValueError(
                    "half-step points root contains an unexpected artifact"
                )
            existing_names.add(child.name)
    prefix = []
    for key in expected_keys:
        if key not in existing_names:
            break
        prefix.append(key)
    if existing_names != set(prefix):
        raise ValueError("half-step point artifacts are not a contiguous prefix")
    current = status.get("current_point")
    if current is not None and current not in expected_keys:
        raise ValueError("half-step status names an unexpected current point")
    if current is not None and current not in existing_names:
        raise RuntimeError(
            "half-step current point has an uncertain unreceipted attempt; "
            "refusing a duplicate evaluation"
        )
    loaded = {}
    progress = []
    density_shape = tuple(status["center_receipt"]["density_shape"])
    if (
        len(density_shape) != 3
        or density_shape[0] != 2
        or density_shape[1] != density_shape[2]
    ):
        raise ValueError("half-step bound center density shape is invalid")
    for key in prefix:
        definition = definitions[key]
        point, point_sha, arrays = _read_record_bundle(
            output_root,
            Path("points") / key,
            expected_receipt_keys=(
                _POINT_OBSERVATION_KEYS
                | _HALF_STEP_POINT_BINDING_KEYS
                | set(definition)
                | _DENSITY_CONTINUITY_KEYS
                | {"schema", "point_key", "kind", "artifacts"}
            ),
            expected_artifacts={
                "gradient.f64": ((6, 3), "hartree / bohr"),
                "density.f64": (density_shape, "electrons"),
            },
        )
        if (
            point.get("schema") != HALF_STEP_SCHEMA
            or point.get("point_key") != key
            or point.get("kind") != "new-half-step-displacement"
            or any(point.get(field) != value for field, value in definition.items())
            or any(
                point.get(field) != value
                for field, value in _half_step_point_bindings(status).items()
            )
        ):
            raise ValueError(f"half-step durable point identity mismatch: {key}")
        loaded[key] = (point, arrays, point_sha)
        progress.append({"point_key": key, "receipt_sha256": point_sha})
    declared = status.get("completed_points")
    if type(declared) is not list or declared != progress[: len(declared)]:
        raise ValueError("half-step status progress diverges from durable artifacts")
    if len(declared) > len(progress):
        raise ValueError("half-step status claims a missing durable point")
    if declared != progress:
        status["completed_points"] = progress
        _claimed_write(claim, output_root / "status.json", _json_bytes(status))
    return loaded


def _half_step_point_gates(
    source: dict[str, Any],
    new_points: list[tuple[dict[str, Any], np.ndarray]],
) -> list[dict[str, Any]]:
    """Apply every unchanged point gate while retaining each point verdict."""

    gates = []
    center = source["center_record"]
    center_density = source["center_density"]
    center_spin = _validate_point_record(
        center,
        center_spin_square=None,
        center=True,
        enforce_scientific_gates=False,
    )
    center_checks = (
        (
            "prior_center_scf_converged",
            center["scf_converged"],
            True,
            "equals",
        ),
        (
            "prior_center_grid_point_count",
            center["grid_point_count"],
            HALF_STEP_CONTRACT["thresholds"]["required_grid_point_count"],
            "equals",
        ),
        (
            "prior_center_spin_square_minimum",
            center_spin,
            HALF_STEP_CONTRACT["thresholds"]["spin_square_inclusive_minimum"],
            "greater_than_or_equal",
        ),
        (
            "prior_center_spin_square_maximum",
            center_spin,
            HALF_STEP_CONTRACT["thresholds"]["spin_square_inclusive_maximum"],
            "less_than_or_equal",
        ),
        (
            "prior_center_fmax",
            center["physical_fmax_ev_per_angstrom"],
            HALF_STEP_CONTRACT["thresholds"][
                "center_fmax_ev_per_angstrom_exclusive_maximum"
            ],
            "strictly_less_than",
        ),
        (
            "prior_center_energy_delta",
            center["electronic_energy_absolute_delta_hartree"],
            HALF_STEP_CONTRACT["thresholds"][
                "center_energy_vs_reference_exclusive_maximum_delta_hartree"
            ],
            "strictly_less_than",
        ),
        (
            "prior_center_s2_delta",
            center["s2_absolute_delta"],
            HALF_STEP_CONTRACT["thresholds"][
                "center_s2_vs_reference_exclusive_maximum_delta"
            ],
            "strictly_less_than",
        ),
    )
    for name, value, limit, comparison in center_checks:
        _analysis_gate(gates, name, value, limit, comparison)
    for point, density in [*source.get("reused_point_data", []), *new_points]:
        spin = _validate_point_record(
            point,
            center_spin_square=center_spin,
            center=False,
            enforce_scientific_gates=False,
        )
        metrics = _density_continuity_metrics(density, center_density)
        campaign._strict_json_equal(
            {name: point.get(name) for name in _DENSITY_CONTINUITY_KEYS},
            metrics,
            label=f"half-step density continuity {point.get('point_key')}",
        )
        prefix = f"point:{point.get('point_key')}"
        checks = (
            (
                f"{prefix}:scf_converged",
                point["scf_converged"],
                True,
                "equals",
            ),
            (
                f"{prefix}:grid_point_count",
                point["grid_point_count"],
                HALF_STEP_CONTRACT["thresholds"]["required_grid_point_count"],
                "equals",
            ),
            (
                f"{prefix}:spin_square_minimum",
                spin,
                HALF_STEP_CONTRACT["thresholds"]["spin_square_inclusive_minimum"],
                "greater_than_or_equal",
            ),
            (
                f"{prefix}:spin_square_maximum",
                spin,
                HALF_STEP_CONTRACT["thresholds"]["spin_square_inclusive_maximum"],
                "less_than_or_equal",
            ),
            (
                f"{prefix}:spin_square_delta",
                abs(spin - center_spin),
                HALF_STEP_CONTRACT["thresholds"][
                    "displaced_vs_center_spin_square_exclusive_maximum_delta"
                ],
                "strictly_less_than",
            ),
            (
                f"{prefix}:density_delta",
                metrics["final_vs_center_normalized_frobenius_delta"],
                HALF_STEP_CONTRACT["thresholds"][
                    "displaced_final_density_vs_center_exclusive_maximum_"
                    "normalized_frobenius_delta"
                ],
                "strictly_less_than",
            ),
        )
        for name, value, limit, comparison in checks:
            _analysis_gate(gates, name, value, limit, comparison)
    return gates


def _run_half_step_locked(
    prior_root: Path,
    output_root: Path,
    claim: _OutputRootClaim,
) -> dict[str, Any]:
    """Run the immutable 36-call continuation after complete source validation."""

    source = _load_half_step_source(prior_root)
    identity = _fd_execution_identity()
    status_path = output_root / "status.json"
    expected_status = _initial_half_step_status(identity, prior_root, source)
    if status_path.exists() or status_path.is_symlink():
        status, _raw = campaign._read_json_object(
            status_path, label="half-step status", max_bytes=4 * 1024 * 1024
        )
        _validate_half_step_resume_status(status, expected_status)
        if (output_root / TERMINAL).exists() or (output_root / "receipt.json").exists():
            raise ValueError("half-step output root is not resumable running state")
    else:
        if not claim.root_created:
            raise FileExistsError(
                f"half-step extension requires a fresh output root: {output_root}"
            )
        if any(output_root.iterdir()):
            raise FileExistsError("half-step output root contains unowned artifacts")
        status = expected_status
        _claimed_write(claim, status_path, _json_bytes(status))
    try:
        transition_state = source["transition_state"]
        atom_count = len(transition_state.symbols)
        if atom_count != 6:
            raise ValueError(
                "half-step extension requires the exact six-atom frozen TS"
            )
        plan = _half_step_plan(atom_count)
        maximum_calls = HALF_STEP_CONTRACT["stencil"][
            "maximum_new_scf_gradient_evaluations"
        ]
        if len(plan) != maximum_calls:
            raise RuntimeError("half-step displacement-count contract drifted")
        loaded = _resume_half_step_points(output_root, status, plan, claim=claim)
        center_density = source["center_density"]
        center_density_sha = hashlib.sha256(
            np.ascontiguousarray(center_density, dtype="<f8").tobytes()
        ).hexdigest()
        if center_density_sha != source["center_receipt"]["density_sha256"]:
            raise ValueError("half-step bound center density changed after validation")
        center_spin = _validate_point_record(
            source["center_record"],
            center_spin_square=None,
            center=True,
            enforce_scientific_gates=False,
        )
        new_gradients = {}
        new_point_data = []
        new_calls = 0
        for definition in plan:
            key = definition["key"]
            if key not in loaded:
                status["current_point"] = key
                _claimed_write(claim, status_path, _json_bytes(status))
                coords = np.asarray(transition_state.coords, dtype=float).copy()
                coords.reshape(-1)[definition["coordinate_index"]] += (
                    definition["displacement_bohr"] * BOHR_TO_ANGSTROM
                )
                displaced = replace(transition_state, coords=coords)
                claim.verify()
                point, gradient, density = _evaluate_fd_gradient(
                    displaced, source["settings"], central_density=center_density
                )
                new_calls += 1
                if new_calls > maximum_calls:
                    raise RuntimeError("half-step calculator-call bound exceeded")
                claim.verify()
                point.update(
                    {
                        "schema": HALF_STEP_SCHEMA,
                        "point_key": key,
                        "kind": "new-half-step-displacement",
                        **_half_step_point_bindings(status),
                        **definition,
                        **_density_continuity_metrics(density, center_density),
                    }
                )
                point, point_sha = _publish_record_bundle(
                    output_root,
                    Path("points") / key,
                    point,
                    {
                        "gradient.f64": (gradient, "hartree / bohr"),
                        "density.f64": (density, "electrons"),
                    },
                    claim=claim,
                )
                status["completed_points"].append(
                    {"point_key": key, "receipt_sha256": point_sha}
                )
                _claimed_write(claim, status_path, _json_bytes(status))
                loaded[key] = (
                    point,
                    {"gradient.f64": gradient, "density.f64": density},
                    point_sha,
                )
            point, arrays, _point_sha = loaded[key]
            expected_coords = np.asarray(transition_state.coords, dtype=float).copy()
            expected_coords.reshape(-1)[definition["coordinate_index"]] += (
                definition["displacement_bohr"] * BOHR_TO_ANGSTROM
            )
            density = arrays["density.f64"]
            density_sha = hashlib.sha256(
                np.ascontiguousarray(density, dtype="<f8").tobytes()
            ).hexdigest()
            if (
                point.get("geometry_fingerprint")
                != frequency_geometry_fingerprint(
                    replace(transition_state, coords=expected_coords)
                )
                or point.get("density_initial_guess_sha256") != center_density_sha
                or point.get("density_final_sha256") != density_sha
            ):
                raise ValueError(f"half-step point payload binding mismatch: {key}")
            _validate_point_record(
                point,
                center_spin_square=center_spin,
                center=False,
                enforce_scientific_gates=False,
            )
            campaign._strict_json_equal(
                {name: point.get(name) for name in _DENSITY_CONTINUITY_KEYS},
                _density_continuity_metrics(density, center_density),
                label=f"half-step point density continuity {key}",
            )
            new_gradients[
                (definition["coordinate_index"], definition["step_multiplier"])
            ] = arrays["gradient.f64"].reshape(-1)
            new_point_data.append((point, density))
        if new_calls > maximum_calls:
            raise RuntimeError("half-step calculator-call bound exceeded")
        status["current_point"] = None
        _claimed_write(claim, status_path, _json_bytes(status))
        h_inner, h_outer, richardson = _half_step_matrices(
            new_gradients, source["reused_gradients"], 3 * atom_count
        )
        route_vector = campaign._mapped_route_vector(
            transition_state,
            source["masses"],
            source["reactant"].coords,
            source["product"].coords,
        )
        analysis = _complete_half_step_analysis(
            transition_state,
            source["masses"],
            route_vector,
            h_inner,
            h_outer,
            richardson,
            source["reference_matrix"],
            point_gate_failures=_half_step_point_gates(source, new_point_data),
        )
        if (
            type(analysis) is not dict
            or type(analysis.get("confirmation_passed")) is not bool
        ):
            raise RuntimeError("half-step analysis returned a malformed verdict")
        confirmation_passed = analysis["confirmation_passed"]
        symmetric = {
            "H_0.005": 0.5 * (h_inner + h_inner.T),
            "H_0.01": 0.5 * (h_outer + h_outer.T),
            "FD": 0.5 * (richardson + richardson.T),
        }
        matrix_arrays = {
            "H_0.005-raw.f64": (h_inner, "hartree / bohr^2"),
            "H_0.01-raw.f64": (h_outer, "hartree / bohr^2"),
            "richardson-raw.f64": (richardson, "hartree / bohr^2"),
            "H_0.005-symmetric.f64": (
                symmetric["H_0.005"],
                "hartree / bohr^2",
            ),
            "H_0.01-symmetric.f64": (
                symmetric["H_0.01"],
                "hartree / bohr^2",
            ),
            "richardson-symmetric.f64": (symmetric["FD"], "hartree / bohr^2"),
        }
        aggregate_record = {
            "schema": HALF_STEP_SCHEMA,
            "kind": "immutable-half-step-finite-difference-analysis",
            "confirmation_passed": confirmation_passed,
            **_half_step_point_bindings(status),
            "new_point_receipts": status["completed_points"],
            "reused_prior_point_receipts": source["reused_points"],
            "analysis": analysis,
        }
        matrix_contract = {
            name: ((3 * atom_count, 3 * atom_count), units)
            for name, (_value, units) in matrix_arrays.items()
        }
        matrix_path = output_root / "matrices"
        if matrix_path.exists() or matrix_path.is_symlink():
            aggregate, aggregate_sha, observed = _read_record_bundle(
                output_root,
                Path("matrices"),
                expected_receipt_keys=set(aggregate_record) | {"artifacts"},
                expected_artifacts=matrix_contract,
            )
            campaign._strict_json_equal(
                {key: aggregate[key] for key in aggregate_record},
                aggregate_record,
                label="resumed half-step matrix receipt",
            )
            if any(
                not np.array_equal(observed[name], value)
                for name, (value, _units) in matrix_arrays.items()
            ):
                raise ValueError(
                    "resumed half-step matrices diverge from reconstruction"
                )
        else:
            aggregate, aggregate_sha = _publish_record_bundle(
                output_root,
                Path("matrices"),
                aggregate_record,
                matrix_arrays,
                claim=claim,
            )
        finished = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        prior_receipt = source["prior_receipt"]
        reference_receipt = source["reference_receipt"]
        receipt = {
            "schema": HALF_STEP_SCHEMA,
            "state": "completed",
            "accepted_campaign_result": False,
            "confirmation_passed": confirmation_passed,
            "purpose": "immutable half-step FD extension only; no production authority",
            "authority": HALF_STEP_CONTRACT["policy"][
                "passed_authority" if confirmation_passed else "rejected_authority"
            ],
            "route": ROUTE,
            "git_sha": identity["git_sha"],
            "script_sha256": identity["script_sha256"],
            "diagnostic_execution_identity": identity["diagnostic_execution_identity"],
            "diagnostic_execution_identity_sha256": identity[
                "diagnostic_execution_identity_sha256"
            ],
            "prior_fd_root": str(prior_root),
            "prior_fd_receipt_sha256": source["prior_receipt_sha256"],
            "prior_fd_terminal_sha256": source["prior_terminal_sha256"],
            "prior_fd_status_sha256": PRIOR_FD_STATUS_SHA256,
            "prior_fd_git_sha": prior_receipt["git_sha"],
            "prior_fd_script_sha256": prior_receipt["script_sha256"],
            "preflight_root": str(source["preflight_root"]),
            "preflight_receipt_sha256": source["preflight_receipt_sha256"],
            "campaign_identity": prior_receipt["campaign_identity"],
            "reference_root": str(source["reference_root"]),
            "reference_receipt_sha256": REFERENCE_RECEIPT_SHA256,
            "reference_matrix_sha256": REFERENCE_MATRIX_SHA256,
            "reference_git_sha": reference_receipt["git_sha"],
            "reference_script_sha256": reference_receipt["script_sha256"],
            "reference_binding": source["reference_binding"],
            "reference_binding_sha256": status["reference_binding_sha256"],
            "center_receipt": status["center_receipt"],
            "reused_prior_point_receipts": source["reused_points"],
            "new_point_receipts": status["completed_points"],
            "contract": HALF_STEP_CONTRACT,
            "contract_sha256": status["contract_sha256"],
            "aggregate_receipt": {
                "path": "matrices/receipt.json",
                "sha256": aggregate_sha,
            },
            "analysis": aggregate["analysis"],
            "finished_utc": finished,
        }
        receipt_raw = _json_bytes(receipt)
        receipt_sha = hashlib.sha256(receipt_raw).hexdigest()
        _claimed_write(claim, output_root / "receipt.json", receipt_raw)
        status.update(
            {
                "state": "completed",
                "accepted_campaign_result": False,
                "confirmation_passed": confirmation_passed,
                "receipt": "receipt.json",
                "receipt_sha256": receipt_sha,
                "finished_utc": finished,
            }
        )
        terminal = _terminal_payload(status, state="completed", detail=None)
        terminal.update({"receipt": "receipt.json", "receipt_sha256": receipt_sha})
        _terminal_write_noreplace(claim, output_root / TERMINAL, terminal)
        _claimed_write(claim, status_path, _json_bytes(status))
        return receipt
    except BaseException as exc:
        terminal_path = output_root / TERMINAL
        if terminal_path.exists() or terminal_path.is_symlink():
            raise
        detail = f"{type(exc).__name__}: {exc}"
        status.update(
            {
                "state": "failed",
                "accepted_campaign_result": False,
                "confirmation_passed": False,
                "error": detail,
                "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            }
        )
        _claimed_write(claim, status_path, _json_bytes(status))
        _terminal_write_noreplace(
            claim,
            terminal_path,
            _terminal_payload(status, state="failed", detail=detail),
        )
        raise


def run_half_step_extension(prior_root: Path, output_root: Path) -> dict[str, Any]:
    """Run or resume the exact immutable +/-0.005 Bohr extension."""

    prior = campaign._safe_absolute_root(prior_root)
    output = campaign._safe_absolute_root(output_root)
    if str(prior) != PRIOR_FD_ROOT:
        raise ValueError("half-step extension requires the exact prior FD receipt root")
    if output == prior:
        raise ValueError("half-step extension requires a separate fresh output root")
    with _exclusive_output_claim(output) as claim:
        status_path = output / "status.json"
        if not claim.root_created and not (
            status_path.exists() or status_path.is_symlink()
        ):
            raise FileExistsError(
                f"half-step extension requires a fresh output root: {output}"
            )
        try:
            return _run_half_step_locked(prior, output, claim)
        except BaseException as exc:
            if isinstance(exc, ConcurrentRunError):
                raise
            terminal_path = output / TERMINAL
            if terminal_path.exists() or terminal_path.is_symlink():
                raise
            status_path = output / "status.json"
            failure_identity = _half_step_failure_execution_identity()
            status: dict[str, Any] = {
                "schema": HALF_STEP_SCHEMA,
                "state": "failed",
                "route": ROUTE,
                **failure_identity,
                "prior_fd_root": str(prior),
                "prior_fd_receipt_sha256": PRIOR_FD_RECEIPT_SHA256,
                "prior_fd_terminal_sha256": PRIOR_FD_TERMINAL_SHA256,
                "prior_fd_status_sha256": PRIOR_FD_STATUS_SHA256,
                "preflight_root": PRIOR_FD_PREFLIGHT_ROOT,
                "preflight_receipt_sha256": PRIOR_FD_PREFLIGHT_SHA256,
                "reference_root": PRIOR_FD_REFERENCE_ROOT,
                "reference_receipt_sha256": REFERENCE_RECEIPT_SHA256,
                "reference_matrix_sha256": REFERENCE_MATRIX_SHA256,
                "contract": HALF_STEP_CONTRACT,
                "contract_sha256": campaign._canonical_hash(HALF_STEP_CONTRACT),
                "completed_points": [],
                "current_point": None,
                "confirmation_passed": False,
                "accepted_campaign_result": False,
            }
            if status_path.is_file() and not status_path.is_symlink():
                try:
                    loaded, _raw = campaign._read_json_object(
                        status_path,
                        label="failed half-step status",
                        max_bytes=4 * 1024 * 1024,
                    )
                    source = _load_half_step_source(prior)
                    expected = _initial_half_step_status(
                        _fd_execution_identity(), prior, source
                    )
                    _validate_half_step_resume_status(loaded, expected)
                except Exception:
                    pass
                else:
                    status = loaded
            detail = f"{type(exc).__name__}: {exc}"
            status.update(
                {
                    "state": "failed",
                    "accepted_campaign_result": False,
                    "confirmation_passed": False,
                    "error": detail,
                    "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                }
            )
            _claimed_write(claim, status_path, _json_bytes(status))
            _terminal_write_noreplace(
                claim,
                terminal_path,
                _terminal_payload(status, state="failed", detail=detail),
            )
            raise


def _run_analytic_locked(
    preflight_root: Path,
    output_root: Path,
    claim: _OutputRootClaim,
) -> dict[str, Any]:
    claim.verify()
    if not claim.root_created:
        raise FileExistsError(f"diagnostic output root already exists: {output_root}")
    identity = _git_identity()
    status = {
        "schema": SCHEMA,
        "state": "running",
        "route": ROUTE,
        **identity,
        "completed_cases": [],
        "current_case": None,
    }
    _claimed_write(claim, output_root / "status.json", _json_bytes(status))
    try:
        preflight, preflight_sha = campaign._validate_production_boundary(
            preflight_root, ROUTE
        )
        status["preflight_receipt_sha256"] = preflight_sha
        _claimed_write(claim, output_root / "status.json", _json_bytes(status))
        base_settings, _ = campaign._canonical_dft_settings(preflight)
        template = reactions(gpu=True, basis="def2-svp")[ROUTE].cluster
        fingerprints = campaign._trusted_input_fingerprints_from_route_record(
            preflight["routes"][ROUTE]
        )
        inputs = campaign._trusted_route_input_snapshots(
            campaign.DEFAULT_BUNDLE_ROOT, ROUTE, template, fingerprints
        )
        transition_state = inputs["transition_state"]
        masses = np.asarray(preflight["routes"][ROUTE]["masses_amu"], dtype=float)

        cases: list[dict[str, Any]] = []
        symmetric_matrices: dict[str, np.ndarray] = {}
        mode_sets: dict[str, Any] = {}
        for definition in CASES:
            status["current_case"] = definition["name"]
            _claimed_write(claim, output_root / "status.json", _json_bytes(status))
            record, symmetric, modes = _evaluate_case(
                claim,
                output_root,
                definition,
                transition_state,
                inputs["reactant"].coords,
                inputs["product"].coords,
                masses,
                base_settings,
            )
            cases.append(record)
            symmetric_matrices[definition["name"]] = symmetric
            mode_sets[definition["name"]] = modes
            status["completed_cases"] = [case["name"] for case in cases]
            _claimed_write(claim, output_root / "status.json", _json_bytes(status))

        status["current_case"] = None
        _claimed_write(claim, output_root / "status.json", _json_bytes(status))
        reference_name = CASES[-1]["name"]
        reference_matrix = symmetric_matrices[reference_name]
        reference_modes = mode_sets[reference_name]
        comparisons = {}
        for definition in CASES[:-1]:
            name = definition["name"]
            matrix_delta = symmetric_matrices[name] - reference_matrix
            eigenvalue_delta = mode_sets[name].eigenvalues - reference_modes.eigenvalues
            comparisons[name] = {
                "reference": reference_name,
                "symmetric_matrix_delta_spectral_norm": float(
                    np.linalg.norm(matrix_delta, ord=2)
                ),
                "projected_eigenvalue_maximum_absolute_delta": float(
                    np.max(np.abs(eigenvalue_delta))
                ),
                "same_negative_mode_count": bool(
                    np.count_nonzero(mode_sets[name].eigenvalues < -1.0e-8)
                    == np.count_nonzero(reference_modes.eigenvalues < -1.0e-8)
                ),
            }
        receipt = {
            "schema": SCHEMA,
            "state": "completed",
            "accepted_campaign_result": False,
            "purpose": "numerical diagnostic only; no TS qualification or SCT result",
            "route": ROUTE,
            "preflight_receipt_sha256": preflight_sha,
            "campaign_identity": preflight["identity"],
            **identity,
            "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "cases": cases,
            "comparisons_to_reference": comparisons,
        }
        receipt_raw = _json_bytes(receipt)
        receipt_sha = hashlib.sha256(receipt_raw).hexdigest()
        _claimed_write(claim, output_root / "receipt.json", receipt_raw)
        terminal = _terminal_payload(status, state="completed", detail=None)
        terminal["receipt"] = "receipt.json"
        terminal["receipt_sha256"] = receipt_sha
        _terminal_write_noreplace(claim, output_root / TERMINAL, terminal)
        status.update(
            {
                "state": "completed",
                "receipt": "receipt.json",
                "receipt_sha256": receipt_sha,
                "finished_utc": receipt["finished_utc"],
            }
        )
        _claimed_write(claim, output_root / "status.json", _json_bytes(status))
        return receipt
    except BaseException as exc:
        terminal_path = output_root / TERMINAL
        if terminal_path.exists() or terminal_path.is_symlink():
            _validated_terminal(terminal_path)
            raise
        detail = f"{type(exc).__name__}: {exc}"
        status.update(
            {
                "state": "failed",
                "failed_case": status.get("current_case"),
                "error": detail,
                "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            }
        )
        _claimed_write(claim, output_root / "status.json", _json_bytes(status))
        terminal = _terminal_payload(status, state="failed", detail=detail)
        _terminal_write_noreplace(claim, terminal_path, terminal)
        raise


def run(preflight_root: Path, output_root: Path) -> dict[str, Any]:
    """Run the original analytic diagnostic under the external output claim."""

    preflight = campaign._safe_absolute_root(preflight_root)
    output = campaign._safe_absolute_root(output_root)
    with _exclusive_output_claim(output) as claim:
        return _run_analytic_locked(preflight, output, claim)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--create-preflight", action="store_true")
    mode.add_argument("--finalize-if-running", action="store_true")
    mode.add_argument("--finite-difference-confirmation", action="store_true")
    mode.add_argument("--half-step-extension", action="store_true")
    parser.add_argument("--preflight-root", type=Path)
    parser.add_argument("--reference-root", type=Path)
    parser.add_argument("--prior-fd-root", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--nice", type=int, default=10)
    parser.add_argument("--log")
    parser.add_argument("--finite-difference-finalization", action="store_true")
    parser.add_argument("--half-step-finalization", action="store_true")
    args = parser.parse_args()
    if args.half_step_extension:
        if args.prior_fd_root is None:
            parser.error("--prior-fd-root is required with --half-step-extension")
        if args.output_root is None:
            parser.error("--output-root is required with --half-step-extension")
        if args.preflight_root is not None or args.reference_root is not None:
            parser.error(
                "--preflight-root/--reference-root are invalid with "
                "--half-step-extension"
            )
        if args.finite_difference_finalization or args.half_step_finalization:
            parser.error("finalization flags require --finalize-if-running")
        receipt = run_half_step_extension(args.prior_fd_root, args.output_root)
        print(
            json.dumps(
                {
                    "state": receipt["state"],
                    "receipt": str(args.output_root / "receipt.json"),
                    "confirmation_passed": receipt["confirmation_passed"],
                    "accepted_campaign_result": False,
                },
                sort_keys=True,
            )
        )
        return 0
    if args.create_preflight:
        if args.reference_root is not None or args.prior_fd_root is not None:
            parser.error(
                "--reference-root/--prior-fd-root are invalid with --create-preflight"
            )
        if args.preflight_root is None:
            parser.error("--preflight-root is required with --create-preflight")
        if args.finite_difference_finalization or args.half_step_finalization:
            parser.error("finalization flags require --finalize-if-running")
        preflight_root = campaign._safe_absolute_root(args.preflight_root)
        receipt = campaign.create_preflight_receipt(
            campaign.DEFAULT_BUNDLE_ROOT, preflight_root
        )
        print(
            json.dumps(
                {
                    "state": "pending",
                    "identity": receipt["identity"],
                    "preflight": str(preflight_root / "preflight.json"),
                },
                sort_keys=True,
            )
        )
        return 0
    if args.finalize_if_running:
        if args.reference_root is not None or args.prior_fd_root is not None:
            parser.error("provenance roots are invalid with --finalize-if-running")
        if args.finite_difference_finalization and args.half_step_finalization:
            parser.error("finalization mode must be unambiguous")
        if args.output_root is None:
            parser.error("--output-root is required with --finalize-if-running")
        terminal = finalize_if_running(
            args.output_root,
            finite_difference=args.finite_difference_finalization,
            half_step=args.half_step_finalization,
        )
        print(json.dumps(terminal, sort_keys=True))
        return 0
    if args.preflight_root is None:
        parser.error("--preflight-root is required unless --finalize-if-running is set")
    if args.output_root is None:
        parser.error("--output-root is required for a diagnostic run")
    if args.prior_fd_root is not None:
        parser.error("--prior-fd-root requires --half-step-extension")
    if args.finite_difference_finalization or args.half_step_finalization:
        parser.error("finalization flags require --finalize-if-running")
    if args.finite_difference_confirmation:
        if args.reference_root is None:
            parser.error(
                "--reference-root is required with --finite-difference-confirmation"
            )
        receipt = run_finite_difference(
            args.preflight_root,
            args.reference_root,
            args.output_root,
        )
    else:
        if args.reference_root is not None:
            parser.error(
                "--reference-root is valid only with --finite-difference-confirmation"
            )
        receipt = run(args.preflight_root, args.output_root)
    output = {
        "state": receipt["state"],
        "receipt": str(args.output_root / "receipt.json"),
    }
    if receipt.get("schema") in {FINITE_DIFFERENCE_SCHEMA, HALF_STEP_SCHEMA}:
        output.update(
            {
                "confirmation_passed": receipt["confirmation_passed"],
                "accepted_campaign_result": receipt["accepted_campaign_result"],
            }
        )
    print(
        json.dumps(
            output,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
