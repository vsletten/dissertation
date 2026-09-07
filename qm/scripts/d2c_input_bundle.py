#!/usr/bin/env python
"""Freeze the D2b inputs consumed by the D2c SCT campaign.

D2b's production outputs live in a gitignored campaign worktree.  D2c must not
silently depend on that mutable directory, so this tool copies the four
predeclared direct-CC one-water benchmark routes into a small tracked bundle and
binds every byte, geometry, method, and generating revision in a manifest.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import sys
import time
from pathlib import Path
from typing import Any

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    bootstrap_cli(
        "d2c_input_bundle",
        default_run_root=Path(__file__).resolve().parent.parent / "runs",
    )

ROUTES = (
    "h-co-1w-oside",
    "h-co-1w-cside",
    "h-h2co-ch3o-1w",
    "h-h2co-h2-hco-1w",
)
FILES = (
    "ts.xyz",
    "irc_back.xyz",
    "irc_fwd.xyz",
    "results.json",
    "cc-reactant-tz.json",
    "cc-ts-tz.json",
    "cc-product-tz.json",
)
ROLES = ("reactant", "ts", "product")
SCHEMA = "d2c-input-bundle-v1"
SOURCE_REPOSITORY = "vsletten/dissertation"
SOURCE_BRANCH = "agents/D2b-explicit-surface-rates"


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require_digest(value: str, *, label: str) -> str:
    if len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise ValueError(f"{label} must be 64 lowercase hexadecimal digits")
    return value


def _require_revision(value: str) -> str:
    if len(value) != 40 or any(c not in "0123456789abcdef" for c in value):
        raise ValueError("source revision must be a full 40-character Git SHA")
    return value


def _regular_below(root: Path, relative: Path) -> Path:
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"unsafe relative path: {relative}")
    candidate = root / relative
    if candidate.is_symlink() or not candidate.is_file():
        raise ValueError(f"required regular file is missing: {relative}")
    resolved = candidate.resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise ValueError(f"file escapes source root: {relative}")
    return candidate


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object: {path.name}")
    return payload


def _geometry_hash_xyz(path: Path, *, fixed_comment: str | None) -> str:
    lines = path.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"invalid XYZ file: {path.name}")
    try:
        atom_count = int(lines[0])
    except ValueError as exc:
        raise ValueError(f"invalid XYZ atom count: {path.name}") from exc
    atom_lines = lines[2:]
    if len(atom_lines) != atom_count:
        raise ValueError(f"XYZ atom count mismatch: {path.name}")
    for line in atom_lines:
        fields = line.split()
        if len(fields) != 4:
            raise ValueError(f"invalid XYZ atom row: {path.name}")
        coordinates = [float(value) for value in fields[1:]]
        if not all(math.isfinite(value) for value in coordinates):
            raise ValueError(f"non-finite XYZ coordinate: {path.name}")
    comment = lines[1] if fixed_comment is None else fixed_comment
    canonical = f"{atom_count}\n{comment}\n" + "\n".join(atom_lines)
    return hashlib.sha256(canonical.encode()).hexdigest()


def geometry_hash_xyz(path: Path) -> str:
    """Reproduce current D2b's name-insensitive geometry hash."""

    return _geometry_hash_xyz(path, fixed_comment="geometry")


def checkpoint_geometry_hash_xyz(path: Path) -> str:
    """Reproduce pre-fix D2b hashes that included the XYZ comment/name."""

    return _geometry_hash_xyz(path, fixed_comment=None)


def _matches_d2b_geometry_hash(path: Path, expected: str) -> bool:
    # Early D2b outputs predate the name-insensitive hash fix at f60c810 and
    # legitimately bind the checkpoint comment.  Later outputs use "geometry".
    return expected in {geometry_hash_xyz(path), checkpoint_geometry_hash_xyz(path)}


def _route_receipt(source_root: Path, route: str) -> dict[str, Any]:
    paths = {name: _regular_below(source_root, Path(route) / name) for name in FILES}
    result = _load_json(paths["results.json"])
    if result.get("key") != route:
        raise ValueError(f"route identity mismatch: {route}")
    if result.get("classification") != "first-order-saddle":
        raise ValueError(f"route lacks a first-order saddle: {route}")
    if result.get("n_water") != 1 or result.get("cc_delta_source") != "direct-tz":
        raise ValueError(f"route is not a direct-CC one-water result: {route}")

    provenance = result.get("provenance")
    method = result.get("method")
    if not isinstance(provenance, dict) or not isinstance(method, dict):
        raise ValueError(f"route lacks method/provenance: {route}")
    geometry_hashes = provenance.get("geometry_sha256")
    if not isinstance(geometry_hashes, dict) or set(geometry_hashes) != set(ROLES):
        raise ValueError(f"route lacks complete geometry identity: {route}")
    for role in ROLES:
        _require_digest(
            str(geometry_hashes[role]), label=f"{route} {role} geometry hash"
        )

    if not _matches_d2b_geometry_hash(paths["ts.xyz"], geometry_hashes["ts"]):
        raise ValueError(f"transition-state geometry hash mismatch: {route}")
    endpoint_paths = (paths["irc_back.xyz"], paths["irc_fwd.xyz"])
    if not all(
        any(
            _matches_d2b_geometry_hash(path, geometry_hashes[role])
            for path in endpoint_paths
        )
        for role in ("reactant", "product")
    ):
        raise ValueError(f"IRC endpoint geometry hashes do not match receipt: {route}")

    for role in ROLES:
        cc = _load_json(paths[f"cc-{role}-tz.json"])
        identity = cc.get("identity")
        if (
            not isinstance(identity, dict)
            or identity.get("geometry_sha256") != geometry_hashes[role]
        ):
            raise ValueError(
                f"coupled-cluster geometry identity mismatch: {route}/{role}"
            )
        if (
            identity.get("method") != "uhf-uccsd(t)"
            or identity.get("basis") != "cc-pvtz"
        ):
            raise ValueError(f"unexpected coupled-cluster method: {route}/{role}")

    return {
        "family": result.get("family"),
        "site": result.get("site"),
        "n_water": result["n_water"],
        "method": method,
        "generated_utc": provenance.get("generated_utc"),
        "generating_git_sha": _require_revision(str(provenance.get("git_sha"))),
        "geometry_sha256": geometry_hashes,
        "canonical_checkpoint_geometry_sha256": {
            name: geometry_hash_xyz(paths[name])
            for name in ("irc_back.xyz", "ts.xyz", "irc_fwd.xyz")
        },
        "files": {
            name: {
                "bytes": paths[name].stat().st_size,
                "sha256": sha256_path(paths[name]),
            }
            for name in FILES
        },
    }


def create_bundle(
    source_root: Path,
    output_root: Path,
    *,
    source_revision: str,
    expected_aggregate_sha256: str,
) -> dict[str, Any]:
    source_root = source_root.resolve()
    output_root = output_root.resolve()
    source_revision = _require_revision(source_revision)
    expected_aggregate_sha256 = _require_digest(
        expected_aggregate_sha256, label="expected aggregate SHA-256"
    )
    aggregate = _regular_below(source_root, Path("results.json"))
    observed_aggregate_sha256 = sha256_path(aggregate)
    if observed_aggregate_sha256 != expected_aggregate_sha256:
        raise ValueError("D2b aggregate results SHA-256 mismatch")
    if output_root == source_root or output_root.is_relative_to(source_root):
        raise ValueError("output root must be outside the mutable source archive")
    if output_root.exists():
        raise FileExistsError(f"output bundle already exists: {output_root}")

    routes = {route: _route_receipt(source_root, route) for route in ROUTES}
    temporary = output_root.with_name(f".{output_root.name}.{time.time_ns()}.tmp")
    temporary.mkdir(parents=True)
    try:
        for route in ROUTES:
            destination = temporary / route
            destination.mkdir()
            for name in FILES:
                shutil.copyfile(source_root / route / name, destination / name)
        manifest = {
            "schema": SCHEMA,
            "purpose": "Immutable direct-CC one-water D2b inputs for D2c SCT",
            "source": {
                "repository": SOURCE_REPOSITORY,
                "branch": SOURCE_BRANCH,
                "revision": source_revision,
                "aggregate_results_sha256": observed_aggregate_sha256,
            },
            "selection": {
                "rule": (
                    "predeclared direct-CC one-water routes; both H+CO orientations "
                    "and matched-label H2CO channels"
                ),
                "routes": list(ROUTES),
            },
            "routes": routes,
        }
        (temporary / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
        temporary.rename(output_root)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    verify_bundle(output_root)
    return manifest


def verify_bundle(bundle_root: Path) -> dict[str, Any]:
    bundle_root = bundle_root.resolve()
    manifest_path = _regular_below(bundle_root, Path("manifest.json"))
    manifest = _load_json(manifest_path)
    if manifest.get("schema") != SCHEMA:
        raise ValueError("unexpected D2c input-bundle schema")
    if manifest.get("selection", {}).get("routes") != list(ROUTES):
        raise ValueError("D2c route selection drifted")
    routes = manifest.get("routes")
    if not isinstance(routes, dict) or set(routes) != set(ROUTES):
        raise ValueError("D2c manifest route inventory drifted")
    for route in ROUTES:
        files = routes[route].get("files")
        if not isinstance(files, dict) or set(files) != set(FILES):
            raise ValueError(f"D2c file inventory drifted: {route}")
        for name in FILES:
            path = _regular_below(bundle_root, Path(route) / name)
            receipt = files[name]
            if path.stat().st_size != receipt.get("bytes") or sha256_path(
                path
            ) != receipt.get("sha256"):
                raise ValueError(f"D2c bundled input drifted: {route}/{name}")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    # Parsed first by ``bootstrap_cli`` so limits apply before campaign imports;
    # retained here so the full parser accepts the same invocation.
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--nice", type=int, default=10)
    parser.add_argument("--log")
    subparsers = parser.add_subparsers(dest="command", required=True)
    create = subparsers.add_parser("create")
    create.add_argument("--source-root", type=Path, required=True)
    create.add_argument("--output-root", type=Path, required=True)
    create.add_argument("--source-revision", required=True)
    create.add_argument("--expected-aggregate-sha256", required=True)
    verify = subparsers.add_parser("verify")
    verify.add_argument("--bundle-root", type=Path, required=True)
    args = parser.parse_args(argv)

    if args.command == "create":
        payload = create_bundle(
            args.source_root,
            args.output_root,
            source_revision=args.source_revision,
            expected_aggregate_sha256=args.expected_aggregate_sha256,
        )
    else:
        payload = verify_bundle(args.bundle_root)
    print(json.dumps({"status": "ok", "schema": payload["schema"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
