#!/usr/bin/env python
"""Reproducible raw-Hessian convergence diagnostic for D2c's frozen TS.

This command is evidence-only: it never publishes a qualification receipt or an
accepted campaign result.  It preserves raw component matrices, evaluates the
strict pre-symmetrization gate, and reports the post-symmetrization vibrational
verdict for a fixed baseline/response/grid/reference matrix.
"""

from __future__ import annotations

# Import order is part of the provenance witness bootstrap below.
# ruff: noqa: E402, I001

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any

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

SCHEMA = "d2c-hessian-symmetry-diagnostic-v1"
TERMINAL = "terminal.json"
ROUTE = "h-co-1w-oside"
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


def _atomic_write(path: Path, raw: bytes) -> None:
    temporary = path.parent / f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp"
    with temporary.open("xb") as handle:
        handle.write(raw)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    directory_fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def _terminal_payload(
    status: dict[str, Any], *, state: str, detail: str | None
) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "state": state,
        "route": ROUTE,
        "git_sha": status.get("git_sha"),
        "script_sha256": status.get("script_sha256"),
        "preflight_receipt_sha256": status.get("preflight_receipt_sha256"),
        "completed_cases": status.get("completed_cases", []),
        "current_case": status.get("current_case"),
        "detail": detail,
        "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }


def finalize_if_running(output_root: Path) -> dict[str, Any]:
    """Emit the dead-man terminal receipt if the bounded unit left none."""

    terminal_path = output_root / TERMINAL
    if terminal_path.is_file() and not terminal_path.is_symlink():
        return json.loads(terminal_path.read_text())
    status_path = output_root / "status.json"
    status: dict[str, Any] = {}
    if status_path.is_file() and not status_path.is_symlink():
        loaded = json.loads(status_path.read_text())
        if isinstance(loaded, dict):
            status = loaded
    terminal = _terminal_payload(
        status,
        state="failed",
        detail=(
            "DiagnosticInterrupted: bounded systemd unit ended without a "
            "terminal receipt"
        ),
    )
    output_root.mkdir(mode=0o700, parents=True, exist_ok=True)
    _atomic_write(terminal_path, _json_bytes(terminal))
    status.update(
        {
            "schema": SCHEMA,
            "state": "failed",
            "route": ROUTE,
            "error": terminal["detail"],
            "finished_utc": terminal["finished_utc"],
        }
    )
    _atomic_write(status_path, _json_bytes(status))
    return terminal


def _matrix_artifact(
    root: Path, case: str, label: str, matrix: np.ndarray
) -> dict[str, Any]:
    array = np.ascontiguousarray(matrix, dtype="<f8")
    raw = array.tobytes()
    relative = Path(case) / f"{label}.f64"
    destination = root / relative
    destination.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    _atomic_write(destination, raw)
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


def _git_identity() -> dict[str, str]:
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
    script_raw = Path(__file__).read_bytes()
    return {"git_sha": sha, "script_sha256": hashlib.sha256(script_raw).hexdigest()}


def _evaluate_case(
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
            "artifact": _matrix_artifact(root, definition["name"], label, matrix),
            "symmetry": _metric_payload(matrix),
        }
    symmetric_artifact = _matrix_artifact(
        root, definition["name"], "total-symmetric", symmetric
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


def run(preflight_root: Path, output_root: Path) -> dict[str, Any]:
    if output_root.exists() or output_root.is_symlink():
        raise FileExistsError(f"diagnostic output root already exists: {output_root}")
    output_root.mkdir(mode=0o700, parents=True)
    identity = _git_identity()
    status = {
        "schema": SCHEMA,
        "state": "running",
        "route": ROUTE,
        **identity,
        "completed_cases": [],
        "current_case": None,
    }
    _atomic_write(output_root / "status.json", _json_bytes(status))
    try:
        preflight, preflight_sha = campaign._validate_production_boundary(
            preflight_root, ROUTE
        )
        status["preflight_receipt_sha256"] = preflight_sha
        _atomic_write(output_root / "status.json", _json_bytes(status))
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
            _atomic_write(output_root / "status.json", _json_bytes(status))
            record, symmetric, modes = _evaluate_case(
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
            _atomic_write(output_root / "status.json", _json_bytes(status))

        status["current_case"] = None
        _atomic_write(output_root / "status.json", _json_bytes(status))
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
        _atomic_write(output_root / "receipt.json", _json_bytes(receipt))
        terminal = _terminal_payload(status, state="completed", detail=None)
        terminal["receipt"] = "receipt.json"
        _atomic_write(output_root / TERMINAL, _json_bytes(terminal))
        status.update(
            {
                "state": "completed",
                "receipt": "receipt.json",
                "finished_utc": receipt["finished_utc"],
            }
        )
        _atomic_write(output_root / "status.json", _json_bytes(status))
        return receipt
    except BaseException as exc:
        detail = f"{type(exc).__name__}: {exc}"
        status.update(
            {
                "state": "failed",
                "failed_case": status.get("current_case"),
                "error": detail,
                "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            }
        )
        _atomic_write(output_root / "status.json", _json_bytes(status))
        terminal = _terminal_payload(status, state="failed", detail=detail)
        _atomic_write(output_root / TERMINAL, _json_bytes(terminal))
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--create-preflight", action="store_true")
    mode.add_argument("--finalize-if-running", action="store_true")
    parser.add_argument("--preflight-root", type=Path)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--nice", type=int, default=10)
    parser.add_argument("--log")
    args = parser.parse_args()
    if args.create_preflight:
        if args.preflight_root is None:
            parser.error("--preflight-root is required with --create-preflight")
        receipt = campaign.create_preflight_receipt(
            campaign.DEFAULT_BUNDLE_ROOT, args.preflight_root.resolve()
        )
        print(
            json.dumps(
                {
                    "state": "pending",
                    "identity": receipt["identity"],
                    "preflight": str(args.preflight_root.resolve() / "preflight.json"),
                },
                sort_keys=True,
            )
        )
        return 0
    if args.finalize_if_running:
        if args.output_root is None:
            parser.error("--output-root is required with --finalize-if-running")
        terminal = finalize_if_running(args.output_root.resolve())
        print(json.dumps(terminal, sort_keys=True))
        return 0
    if args.preflight_root is None:
        parser.error("--preflight-root is required unless --finalize-if-running is set")
    if args.output_root is None:
        parser.error("--output-root is required for a diagnostic run")
    receipt = run(args.preflight_root.resolve(), args.output_root.resolve())
    print(
        json.dumps(
            {
                "state": receipt["state"],
                "receipt": str(args.output_root / "receipt.json"),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
