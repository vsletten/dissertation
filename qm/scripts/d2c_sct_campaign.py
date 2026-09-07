#!/usr/bin/env python
"""D2c SCT campaign preflight and immutable resume contract.

Only ``--dry-run`` is implemented.  The command verifies the frozen D2b bundle,
validates all four transition-state atom mappings against the existing D2b reaction
templates, and writes one non-overwriting pending receipt.  It deliberately performs
no IRC, Hessian, high-level, VAG, or tunnelling calculation.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import subprocess
import sys
import time
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

from scripts import d2c_input_bundle
from scripts.production_energetics import load_xyz_like
from scripts.surface_rate_protocol import reactions

SCHEMA = "d2c-sct-campaign-preflight-v1"
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
# Ground-state isotopic masses.  Carbon-12 is exact by definition; the others are
# the neutral-atom masses used by this bounded H/C/O campaign.
ISOTOPIC_MASSES_AMU = {
    "H": 1.00782503223,
    "C": 12.0,
    "O": 15.99491461957,
}
REFERENCE_MASS_AMU = 1.0
BOUNDS: dict[str, Any] = {
    "irc": {
        "directions": ["forward", "reverse"],
        "algorithm": "sella-gonzalez-schlegel",
        "step_size_angstrom": 0.05,
        "maximum_steps_per_direction": 400,
        "outer_fmax_ev_per_angstrom": 0.05,
        "inner_fmax_ev_per_angstrom": 0.01,
    },
    "hessian": {
        "coverage": "every retained IRC point including TS and endpoints",
        "maximum_retained_points_per_direction": 201,
        "transverse_negative_eigenvalue_tolerance": 1.0e-8,
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
    },
}
ROUTE_STAGE_CONTRACT: dict[str, dict[str, Any]] = {
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
        "requirement": "all four routes compared from one explicitly common reference",
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


def _dependency_versions() -> dict[str, str]:
    versions: dict[str, str] = {}
    for distribution in DEPENDENCY_DISTRIBUTIONS:
        try:
            versions[distribution] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError as exc:
            raise RuntimeError(
                f"required D2c dependency is not installed: {distribution}"
            ) from exc
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
        inventory[route] = {
            "symbols": list(transition_state.symbols),
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
    if set(dependencies) != set(DEPENDENCY_DISTRIBUTIONS):
        raise ValueError("dependency version inventory must be exact and complete")
    if any(not isinstance(value, str) or not value for value in dependencies.values()):
        raise ValueError("dependency versions must be non-empty strings")
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

    run_root.mkdir(parents=True, exist_ok=True)
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
