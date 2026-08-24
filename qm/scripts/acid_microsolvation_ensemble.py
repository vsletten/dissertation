#!/usr/bin/env python
"""Screen the finite A1c matched proton-relay conformer ensemble.

The campaign is intentionally finite: four deterministic topologies at each of
3--6 explicit waters.  Si--O--Si is screened first.  Only deduplicated,
zero-imaginary protonated-bridge Si minima are screened on matched Si--O--Al
seeds.  Every seed is journaled atomically and all console output is teed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

A1D_SOURCE_MANIFEST_PATH = Path(
    "/mnt/data/vsletten/dissertation-data/"
    "task197-a1c-acid-conformers-20260824/"
    "acid-microsolvation-ensemble-v2-g1-b3lyp-def2-svp/manifest.json"
)
A1D_SOURCE_MANIFEST_SHA256 = (
    "471a099c01c3fd8107409fd087e033bdd5da167cc6d41f09475a6a661cae5ac6"
)
A1D_SOURCE_LOG_SHA256 = (
    "80729b9ed1dd61728f4959f29e3adb8151341a5764a779559ded1610005af312"
)
A1D_TARGET_SEED = "si-acid:3w:bridge-donor-chain"
A1D_ATTEMPT_PARENT_NAME = "task199-a1d-acid-3w-bridge-20260824"
A1D_ATTEMPT_NAME = "acid-3w-bridge-refinement-b3lyp-def2-svp"
A1D_LOG_PATH = (
    A1D_SOURCE_MANIFEST_PATH.parent.parent.parent
    / A1D_ATTEMPT_PARENT_NAME
    / "a1d-refinement.log"
)

_ETIQUETTE_SESSION = None
if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    refinement_requested = any(
        argument == "--refine-failed-endpoint"
        or argument.startswith("--refine-failed-endpoint=")
        for argument in sys.argv[1:]
    )
    if refinement_requested and any(
        argument == "--log" or argument.startswith("--log=")
        for argument in sys.argv[1:]
    ):
        raise SystemExit("refinement log path is fixed; --log is forbidden")
    bootstrap_argv = None
    if refinement_requested:
        bootstrap_argv = [*sys.argv[1:], "--log", str(A1D_LOG_PATH)]
    _ETIQUETTE_SESSION = bootstrap_cli(
        "acid_microsolvation_ensemble",
        default_run_root=Path(__file__).resolve().parent.parent / "runs",
        argv=bootstrap_argv,
    )

import numpy as np  # noqa: E402

from quarry.clusters import (  # noqa: E402
    ACID_MICROSOLVATION_FAMILIES,
    Cluster,
    aluminosilicate_dimer,
    disilicate,
    protonated_bridge_complex,
)
from quarry.pipeline import (  # noqa: E402
    DftSettings,
    frequencies,
    optimize_bounded,
)
from scripts import phase1_xiao_lasaga as phase1  # noqa: E402

SCHEMA_VERSION = 1
WATER_COUNTS = (3, 4, 5, 6)
TERMINAL_STATUSES = frozenset({"accepted", "rejected", "failed", "blocked"})
MIN_PAIR_DISTANCE_A = 0.85
HF_PREOPT_MAX_STEPS = 40
REFINEMENT_SCHEMA_VERSION = 1
REFINEMENT_MAX_STEPS = 160
REFINEMENT_ATTEMPT = "attempt-001"
OPTIMIZER_EXHAUSTION_REASON = (
    "RuntimeError: geometry optimization did not converge within 160 steps"
)
FINALIZATION_RESERVE_BYTES = 64 * 1024


@dataclass(frozen=True)
class BasinCandidate:
    label: str
    n_water: int
    cluster: Cluster
    electronic_hartree: float


@dataclass(frozen=True)
class RefinementSource:
    """Hash-pinned immutable A1c evidence for one A1d continuation."""

    manifest_path: Path
    manifest_sha256: str
    log_path: Path
    log_sha256: str
    manifest: dict[str, Any]
    settings: DftSettings
    target_key: str
    target_record: dict[str, Any]
    artifact_hashes: dict[str, str]


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{time.time_ns()}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    temporary.replace(path)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_xyz_strict(path: Path, template: Cluster) -> Cluster:
    """Load XYZ only when atom count, order, symbols, and coordinates are exact."""
    try:
        lines = path.read_text().splitlines()
        n_atoms = int(lines[0])
        atom_lines = lines[2 : 2 + n_atoms]
        fields = [line.split() for line in atom_lines]
        symbols = [field[0] for field in fields]
        coords = np.asarray(
            [[float(value) for value in field[1:4]] for field in fields],
            dtype=float,
        )
    except (OSError, UnicodeDecodeError, ValueError, IndexError) as exc:
        raise ValueError(f"XYZ artifact is unreadable: {path}: {exc}") from exc
    if n_atoms != len(template.symbols) or len(fields) != n_atoms:
        raise ValueError(f"XYZ atom-count mismatch: {path}")
    if any(len(field) < 4 for field in fields) or symbols != template.symbols:
        raise ValueError(f"XYZ atom-order or symbol mismatch: {path}")
    if coords.shape != template.coords.shape or not np.all(np.isfinite(coords)):
        raise ValueError(f"XYZ coordinates are malformed or non-finite: {path}")
    return phase1.load_xyz(path, template)


def terminal_record_is_reusable(root: Path, record: object) -> bool:
    if not isinstance(record, dict) or record.get("status") not in TERMINAL_STATUSES:
        return False
    artifacts = record.get("artifacts")
    if not isinstance(artifacts, dict) or not artifacts:
        return False
    if record.get("status") == "accepted":
        required = {"seed.xyz", "preoptimized.xyz", "optimized.xyz", "minimum.json"}
        if not required.issubset(artifacts):
            return False
        energy = record.get("electronic_hartree")
        if (
            not isinstance(energy, (int, float))
            or isinstance(energy, bool)
            or not np.isfinite(energy)
            or type(record.get("n_imaginary")) is not int
            or record.get("n_imaginary") != 0
        ):
            return False
    for relative, expected_hash in artifacts.items():
        if not isinstance(relative, str) or not isinstance(expected_hash, str):
            return False
        path = root / relative
        if not path.is_file() or sha256_path(path) != expected_hash:
            return False
    return True


def quarantine_seed_artifacts(seed_dir: Path, *, reason: str) -> Path | None:
    """Move stale/malformed checkpoints aside before recomputing one seed."""
    names = ("seed.xyz", "preoptimized.xyz", "optimized.xyz", "minimum.json")
    existing = [seed_dir / name for name in names if (seed_dir / name).exists()]
    if not existing:
        return None
    stamp = f"{time.strftime('%Y%m%dT%H%M%S')}-{time.time_ns()}-{reason}"
    quarantine = seed_dir / "quarantine" / stamp
    quarantine.mkdir(parents=True)
    for path in existing:
        path.replace(quarantine / path.name)
    return quarantine


def minimum_pair_distance(cluster: Cluster) -> float:
    delta = cluster.coords[:, None, :] - cluster.coords[None, :, :]
    distances = np.linalg.norm(delta, axis=2)
    distances[np.diag_indices_from(distances)] = np.inf
    return float(np.min(distances))


def deduplicate_basins(
    candidates: list[BasinCandidate],
) -> dict[str, str]:
    """Map each basin label to the lowest-energy equivalent heavy-atom basin."""
    assignments: dict[str, str] = {}
    canonicals: list[BasinCandidate] = []
    for candidate in sorted(
        candidates,
        key=lambda item: (item.n_water, item.electronic_hartree, item.label),
    ):
        if not np.isfinite(candidate.electronic_hartree):
            raise ValueError(f"basin {candidate.label} has non-finite energy")
        duplicate_of = None
        for canonical in canonicals:
            if candidate.n_water != canonical.n_water:
                continue
            heavy_atoms = [
                index
                for index, symbol in enumerate(candidate.cluster.symbols)
                if symbol != "H"
            ]
            reason = phase1.basin_equivalence_reason(
                candidate.cluster,
                canonical.cluster,
                heavy_atoms,
                candidate_energy_hartree=candidate.electronic_hartree,
                canonical_energy_hartree=canonical.electronic_hartree,
            )
            if reason is None:
                duplicate_of = canonical.label
                break
        if duplicate_of is None:
            duplicate_of = candidate.label
            canonicals.append(candidate)
        assignments[candidate.label] = duplicate_of
    return assignments


def _seed_key(reaction: str, n_water: int, family: str) -> str:
    return f"{reaction}:{n_water}w:{family}"


def _seed_dir(run_dir: Path, reaction: str, n_water: int, family: str) -> Path:
    return run_dir / reaction / f"{n_water}w" / family


def _template(reaction: str, n_water: int, family: str) -> Cluster:
    dimer = disilicate() if reaction == "si-acid" else aluminosilicate_dimer()
    return protonated_bridge_complex(
        dimer,
        n_water=n_water,
        conformer_family=family,
    )


def accepted_record_receipt_is_valid(
    seed_dir: Path,
    record: dict,
    template: Cluster,
    settings: DftSettings,
    *,
    reaction: str,
    n_water: int,
    family: str,
) -> bool:
    """Verify accepted state against identity and its bound minimum receipt."""
    if (
        record.get("reaction") != reaction
        or record.get("water_count") != n_water
        or record.get("family") != family
    ):
        return False
    try:
        cluster = load_xyz_strict(seed_dir / "optimized.xyz", template)
    except (OSError, ValueError, IndexError):
        return False
    receipt = phase1.load_reactant_minimum_receipt(
        seed_dir / "minimum.json", cluster, settings
    )
    if receipt is None or receipt.get("passed") is not True:
        return False
    return receipt.get("n_imaginary") == record.get("n_imaginary") == 0 and receipt.get(
        "electronic_hartree"
    ) == record.get("electronic_hartree")


def _artifact_hashes(seed_dir: Path, names: tuple[str, ...]) -> dict[str, str]:
    return {
        name: sha256_path(seed_dir / name)
        for name in names
        if (seed_dir / name).is_file()
    }


def _expected_si_keys() -> set[str]:
    return {
        _seed_key("si-acid", n_water, family)
        for n_water in WATER_COUNTS
        for family in ACID_MICROSOLVATION_FAMILIES
    }


def _require_refinement_bound(max_steps: object) -> int:
    if type(max_steps) is not int or not 1 <= max_steps <= REFINEMENT_MAX_STEPS:
        raise ValueError(
            f"refinement max_steps must be an integer in 1..{REFINEMENT_MAX_STEPS}"
        )
    return max_steps


def _artifact_path_below(root: Path, relative: object) -> Path:
    if not isinstance(relative, str):
        raise ValueError("artifact path must be a string")
    relative_path = Path(relative)
    if relative_path.is_absolute() or relative_path.name != relative:
        raise ValueError(f"artifact path must be one filename: {relative!r}")
    path = root / relative_path
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"artifact is missing or symlinked: {relative}")
    if not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"artifact escapes its seed directory: {relative}")
    return path


def validate_refinement_source(
    manifest_path: Path,
    *,
    expected_manifest_sha256: str,
    expected_log_sha256: str,
    target_key: str,
    attempt_dir: Path,
    max_steps: int,
) -> RefinementSource:
    """Load exactly one immutable optimizer-exhausted A1c source record."""
    _require_refinement_bound(max_steps)
    manifest_path = manifest_path.expanduser().resolve()
    attempt_dir = attempt_dir.expanduser().resolve()
    if not manifest_path.is_file() or manifest_path.is_symlink():
        raise ValueError("source manifest is missing or symlinked")
    source_root = manifest_path.parent
    if attempt_dir == source_root or attempt_dir.is_relative_to(source_root):
        raise ValueError("refinement attempt directory must be outside source archive")
    if (
        not isinstance(expected_manifest_sha256, str)
        or len(expected_manifest_sha256) != 64
        or any(
            character not in "0123456789abcdef"
            for character in expected_manifest_sha256
        )
    ):
        raise ValueError(
            "expected source manifest SHA-256 must be 64 lowercase hex digits"
        )
    manifest_sha256 = sha256_path(manifest_path)
    if manifest_sha256 != expected_manifest_sha256:
        raise ValueError("source manifest SHA-256 mismatch")
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"source manifest is unreadable: {exc}") from exc
    if not isinstance(manifest, dict):
        raise ValueError("source manifest must be a JSON object")

    required_identity = {
        "schema_version": SCHEMA_VERSION,
        "conformer_version": phase1.ACID_CONFORMER_VERSION,
        "gate_version": phase1.ACID_CONFORMER_GATE_VERSION,
        "max_steps": REFINEMENT_MAX_STEPS,
        "water_counts": list(WATER_COUNTS),
        "families": list(ACID_MICROSOLVATION_FAMILIES),
    }
    for field, expected in required_identity.items():
        if manifest.get(field) != expected:
            raise ValueError(f"source manifest {field} mismatch")
    settings = DftSettings(
        xc="b3lyp",
        basis="def2-svp",
        density_fit=True,
        use_gpu=True,
    )
    if manifest.get("settings") != asdict(settings):
        raise ValueError("source manifest settings mismatch")
    if manifest.get("preoptimization") != {
        "xc": "hf",
        "basis": "sto-3g",
        "max_steps": HF_PREOPT_MAX_STEPS,
        "convergence_is_advisory": True,
    }:
        raise ValueError("source manifest preoptimization mismatch")

    seeds = manifest.get("seeds")
    if not isinstance(seeds, dict) or set(seeds) != _expected_si_keys():
        raise ValueError("source manifest must contain the exact 16-seed Si population")
    if target_key not in seeds:
        raise ValueError("target refinement seed is absent from source manifest")
    failed_keys = {
        key for key, record in seeds.items() if record.get("status") == "failed"
    }
    if failed_keys != {target_key}:
        raise ValueError("target must be the sole failed source seed")

    artifact_hashes: dict[str, str] = {}
    for key, record in sorted(seeds.items()):
        if (
            not isinstance(record, dict)
            or record.get("status") not in TERMINAL_STATUSES
        ):
            raise ValueError(f"source seed {key} is not terminal")
        reaction, water_token, family = key.split(":", maxsplit=2)
        n_water = int(water_token.removesuffix("w"))
        if (
            record.get("reaction") != reaction
            or record.get("water_count") != n_water
            or record.get("family") != family
        ):
            raise ValueError(f"source seed identity mismatch: {key}")
        if key != target_key and record.get("status") != "rejected":
            raise ValueError(
                f"non-target source seed is not conclusively rejected: {key}"
            )
        seed_dir = _seed_dir(source_root, reaction, n_water, family)
        artifacts = record.get("artifacts")
        if not isinstance(artifacts, dict) or not artifacts:
            raise ValueError(f"source seed lacks artifacts: {key}")
        for name, expected_hash in sorted(artifacts.items()):
            path = _artifact_path_below(seed_dir, name)
            if not isinstance(expected_hash, str) or sha256_path(path) != expected_hash:
                raise ValueError(f"source artifact SHA-256 mismatch: {key}/{name}")
            artifact_hashes[str(path.relative_to(source_root))] = expected_hash

    target_record = seeds[target_key]
    if (
        target_record.get("production_converged") is not False
        or target_record.get("reason") != OPTIMIZER_EXHAUSTION_REASON
    ):
        raise ValueError("target source record is not the exact optimizer exhaustion")
    target_artifacts = target_record.get("artifacts", {})
    if not {"seed.xyz", "preoptimized.xyz", "optimized.xyz"}.issubset(target_artifacts):
        raise ValueError("target source record lacks immutable endpoint provenance")
    reaction, water_token, family = target_key.split(":", maxsplit=2)
    n_water = int(water_token.removesuffix("w"))
    target_dir = _seed_dir(source_root, reaction, n_water, family)
    endpoint = load_xyz_strict(
        target_dir / "optimized.xyz", _template(reaction, n_water, family)
    )
    if not np.all(np.isfinite(endpoint.coords)):
        raise ValueError("target source endpoint has non-finite coordinates")
    if minimum_pair_distance(endpoint) <= MIN_PAIR_DISTANCE_A:
        raise ValueError("target source endpoint fails the minimum-pair-distance gate")

    log_path = source_root / "ensemble.log"
    if log_path.is_symlink() or not log_path.is_file():
        raise ValueError("archived source log is missing or symlinked")
    log_sha256 = sha256_path(log_path)
    if log_sha256 != expected_log_sha256:
        raise ValueError("source log SHA-256 mismatch")
    return RefinementSource(
        manifest_path=manifest_path,
        manifest_sha256=manifest_sha256,
        log_path=log_path,
        log_sha256=log_sha256,
        manifest=manifest,
        settings=settings,
        target_key=target_key,
        target_record=target_record,
        artifact_hashes=artifact_hashes,
    )


def assert_refinement_source_unchanged(source: RefinementSource) -> None:
    """Recheck every frozen source byte after expensive work and before publish."""
    if sha256_path(source.manifest_path) != source.manifest_sha256:
        raise RuntimeError("source manifest drifted during refinement")
    if sha256_path(source.log_path) != source.log_sha256:
        raise RuntimeError("source log drifted during refinement")
    root = source.manifest_path.parent
    for relative, expected_hash in source.artifact_hashes.items():
        path = root / relative
        if (
            path.is_symlink()
            or not path.is_file()
            or sha256_path(path) != expected_hash
        ):
            raise RuntimeError(f"source artifact drifted during refinement: {relative}")


def refinement_attempt_dir(manifest_path: Path) -> Path:
    """Return the one card-owned attempt namespace for the pinned source."""
    source_root = manifest_path.expanduser().resolve().parent
    return source_root.parent.parent / A1D_ATTEMPT_PARENT_NAME / A1D_ATTEMPT_NAME


def load_a1d_refinement_source(*, max_steps: int) -> RefinementSource:
    """Load only the card-pinned A1c source; no caller chooses authority."""
    attempt_dir = refinement_attempt_dir(A1D_SOURCE_MANIFEST_PATH)
    return validate_refinement_source(
        A1D_SOURCE_MANIFEST_PATH,
        expected_manifest_sha256=A1D_SOURCE_MANIFEST_SHA256,
        expected_log_sha256=A1D_SOURCE_LOG_SHA256,
        target_key=A1D_TARGET_SEED,
        attempt_dir=attempt_dir,
        max_steps=max_steps,
    )


def _record_running(
    manifest: dict,
    manifest_path: Path,
    *,
    key: str,
    reaction: str,
    n_water: int,
    family: str,
    seed_dir: Path,
) -> None:
    manifest["seeds"][key] = {
        "status": "running",
        "reaction": reaction,
        "water_count": n_water,
        "family": family,
        "seed_dir": str(seed_dir.relative_to(manifest_path.parent)),
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    atomic_json(manifest_path, manifest)


def evaluate_optimized_endpoint(
    seed_dir: Path,
    optimized: Cluster,
    settings: DftSettings,
    *,
    reaction: str,
    n_water: int,
    family: str,
) -> dict[str, Any]:
    """Apply the shared structural, occupancy, and minimum gates."""
    template = _template(reaction, n_water, family)
    if (
        optimized.symbols != template.symbols
        or optimized.coords.shape != template.coords.shape
    ):
        raise RuntimeError("optimized endpoint changed atom count or ordering")
    if not np.all(np.isfinite(optimized.coords)):
        raise RuntimeError("optimized endpoint has non-finite coordinates")
    pair_distance = minimum_pair_distance(optimized)
    if pair_distance <= MIN_PAIR_DISTANCE_A:
        return {
            "status": "blocked",
            "reason": (
                f"optimized minimum pair distance {pair_distance:.6f} A is not above "
                f"{MIN_PAIR_DISTANCE_A:.2f} A"
            ),
            "basin_signature": None,
            "n_imaginary": None,
            "electronic_hartree": None,
            "minimum_pair_distance_a": pair_distance,
        }

    dimer = disilicate() if reaction == "si-acid" else aluminosilicate_dimer()
    ow_index = len(dimer.symbols)
    proton_indices, solvent_oxygen_indices = phase1.acid_mobile_indices(
        ow_index, n_water
    )
    signature = phase1.acid_basin_signature(
        optimized,
        ow_index,
        proton_indices,
        solvent_oxygen_indices=solvent_oxygen_indices,
    )
    reason = phase1.protonated_bridge_reason(
        optimized,
        ow_index,
        proton_indices,
        solvent_oxygen_indices=solvent_oxygen_indices,
    )
    n_imaginary = None
    electronic_hartree = None
    if reason is None:
        frequency = frequencies(optimized, settings)
        minimum_path = seed_dir / "minimum.json"
        phase1.write_reactant_minimum_receipt(
            minimum_path, optimized, settings, frequency
        )
        minimum_receipt = phase1.load_reactant_minimum_receipt(
            minimum_path, optimized, settings
        )
        if minimum_receipt is None:
            raise RuntimeError(
                "frequency result failed geometry/settings/Hessian receipt validation"
            )
        n_imaginary_value = minimum_receipt.get("n_imaginary")
        energy_value = minimum_receipt.get("electronic_hartree")
        if (
            type(n_imaginary_value) is not int
            or not isinstance(energy_value, (int, float))
            or isinstance(energy_value, bool)
            or not np.isfinite(energy_value)
        ):
            raise RuntimeError("validated minimum receipt has invalid numeric values")
        n_imaginary = n_imaginary_value
        electronic_hartree = float(energy_value)
        if n_imaginary:
            reason = f"optimized protonated bridge has {n_imaginary} imaginary modes"
    return {
        "status": "accepted" if reason is None else "rejected",
        "reason": reason,
        "basin_signature": list(signature),
        "n_imaginary": n_imaginary,
        "electronic_hartree": electronic_hartree,
        "minimum_pair_distance_a": pair_distance,
    }


def screen_seed(
    manifest: dict,
    manifest_path: Path,
    *,
    reaction: str,
    n_water: int,
    family: str,
    settings: DftSettings,
    max_steps: int,
) -> dict:
    run_dir = manifest_path.parent
    seed_dir = _seed_dir(run_dir, reaction, n_water, family)
    seed_dir.mkdir(parents=True, exist_ok=True)
    key = _seed_key(reaction, n_water, family)
    existing = manifest["seeds"].get(key)
    reusable = terminal_record_is_reusable(seed_dir, existing)
    if reusable and existing["status"] == "accepted":
        reusable = accepted_record_receipt_is_valid(
            seed_dir,
            existing,
            _template(reaction, n_water, family),
            settings,
            reaction=reaction,
            n_water=n_water,
            family=family,
        )
    if reusable:
        print(
            f"A1C_RESUME seed={key} status={existing['status']}",
            flush=True,
        )
        return existing
    if existing is not None:
        quarantine = quarantine_seed_artifacts(seed_dir, reason="invalid-checkpoint")
        if quarantine is not None:
            print(f"A1C_QUARANTINE seed={key} path={quarantine}", flush=True)

    seed = _template(reaction, n_water, family)
    pair_distance = minimum_pair_distance(seed)
    seed_path = seed_dir / "seed.xyz"
    phase1.save_xyz(seed, seed_path)
    _record_running(
        manifest,
        manifest_path,
        key=key,
        reaction=reaction,
        n_water=n_water,
        family=family,
        seed_dir=seed_dir,
    )
    print(
        f"A1C_START seed={key} atoms={len(seed.symbols)} min_pair={pair_distance:.6f}",
        flush=True,
    )

    if pair_distance <= MIN_PAIR_DISTANCE_A:
        record = {
            "status": "blocked",
            "reaction": reaction,
            "water_count": n_water,
            "family": family,
            "reason": (
                f"seed minimum pair distance {pair_distance:.6f} A is not above "
                f"{MIN_PAIR_DISTANCE_A:.2f} A"
            ),
            "artifacts": _artifact_hashes(seed_dir, ("seed.xyz",)),
        }
    else:
        preoptimization_converged = None
        production_converged = None
        try:
            preopt_result = optimize_bounded(
                seed,
                DftSettings(xc="hf", basis="sto-3g"),
                max_steps=HF_PREOPT_MAX_STEPS,
            )
            preoptimization_converged = preopt_result.converged
            preopt = preopt_result.cluster
            phase1.save_xyz(preopt, seed_dir / "preoptimized.xyz")
            production_result = optimize_bounded(preopt, settings, max_steps=max_steps)
            production_converged = production_result.converged
            optimized = production_result.cluster
            optimized_path = seed_dir / "optimized.xyz"
            phase1.save_xyz(optimized, optimized_path)
            optimized = load_xyz_strict(optimized_path, preopt)
            if not production_converged:
                raise RuntimeError(
                    f"geometry optimization did not converge within {max_steps} steps"
                )
            record = evaluate_optimized_endpoint(
                seed_dir,
                optimized,
                settings,
                reaction=reaction,
                n_water=n_water,
                family=family,
            )
            record.update(
                {
                    "reaction": reaction,
                    "water_count": n_water,
                    "family": family,
                    "preoptimization_converged": preoptimization_converged,
                    "production_converged": production_converged,
                    "artifacts": _artifact_hashes(
                        seed_dir,
                        (
                            "seed.xyz",
                            "preoptimized.xyz",
                            "optimized.xyz",
                            "minimum.json",
                        ),
                    ),
                }
            )
        except Exception as exc:  # every bounded seed keeps its own failure receipt
            record = {
                "status": "failed",
                "reaction": reaction,
                "water_count": n_water,
                "family": family,
                "reason": f"{type(exc).__name__}: {exc}",
                "preoptimization_converged": preoptimization_converged,
                "production_converged": production_converged,
                "artifacts": _artifact_hashes(
                    seed_dir,
                    ("seed.xyz", "preoptimized.xyz", "optimized.xyz", "minimum.json"),
                ),
            }
    record["finished_at"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    manifest["seeds"][key] = record
    atomic_json(manifest_path, manifest)
    print(
        f"A1C_DONE seed={key} status={record['status']} reason={record.get('reason')}",
        flush=True,
    )
    phase1.trim_gpu_pool()
    return record


def _load_candidate(
    run_dir: Path,
    *,
    reaction: str,
    n_water: int,
    family: str,
    record: dict,
    settings: DftSettings,
) -> BasinCandidate:
    template = _template(reaction, n_water, family)
    seed_dir = _seed_dir(run_dir, reaction, n_water, family)
    cluster = load_xyz_strict(seed_dir / "optimized.xyz", template)
    if not accepted_record_receipt_is_valid(
        seed_dir,
        record,
        template,
        settings,
        reaction=reaction,
        n_water=n_water,
        family=family,
    ):
        raise RuntimeError(
            f"accepted seed {_seed_key(reaction, n_water, family)} has invalid receipt"
        )
    energy = record.get("electronic_hartree")
    if (
        not isinstance(energy, (int, float))
        or isinstance(energy, bool)
        or not np.isfinite(energy)
    ):
        raise RuntimeError(
            f"accepted seed {_seed_key(reaction, n_water, family)} lacks finite energy"
        )
    return BasinCandidate(
        label=_seed_key(reaction, n_water, family),
        n_water=n_water,
        cluster=cluster,
        electronic_hartree=float(energy),
    )


def _expected_manifest(settings: DftSettings, max_steps: int, log_path: str) -> dict:
    return {
        "schema_version": SCHEMA_VERSION,
        "conformer_version": phase1.ACID_CONFORMER_VERSION,
        "gate_version": phase1.ACID_CONFORMER_GATE_VERSION,
        "settings": asdict(settings),
        "preoptimization": {
            "xc": "hf",
            "basis": "sto-3g",
            "max_steps": HF_PREOPT_MAX_STEPS,
            "convergence_is_advisory": True,
        },
        "max_steps": max_steps,
        "water_counts": list(WATER_COUNTS),
        "families": list(ACID_MICROSOLVATION_FAMILIES),
        "log": log_path,
        "seeds": {},
    }


def write_matched_candidate_receipt(
    run_dir: Path,
    manifest_path: Path,
    manifest: dict,
    settings: DftSettings,
) -> Path:
    """Bind matched Si/Al minima to one tamper-evident barrier handoff."""
    candidates = []
    for matched in manifest["summary"]["matched_candidates"]:
        n_water = matched["water_count"]
        family = matched["family"]
        models = {}
        for reaction in ("si-acid", "al-acid"):
            key = _seed_key(reaction, n_water, family)
            record = manifest["seeds"][key]
            seed_dir = _seed_dir(run_dir, reaction, n_water, family)
            template = _template(reaction, n_water, family)
            if not accepted_record_receipt_is_valid(
                seed_dir,
                record,
                template,
                settings,
                reaction=reaction,
                n_water=n_water,
                family=family,
            ):
                raise RuntimeError(f"matched candidate {key} lost its accepted receipt")
            models[reaction] = {
                "seed_key": key,
                "optimized_path": str(
                    (seed_dir / "optimized.xyz").relative_to(run_dir)
                ),
                "optimized_sha256": record["artifacts"]["optimized.xyz"],
                "minimum_path": str((seed_dir / "minimum.json").relative_to(run_dir)),
                "minimum_sha256": record["artifacts"]["minimum.json"],
                "electronic_hartree": record["electronic_hartree"],
            }
        candidates.append(
            {
                "water_count": n_water,
                "family": family,
                "models": models,
            }
        )
    payload = {
        "schema_version": 1,
        "conformer_version": phase1.ACID_CONFORMER_VERSION,
        "gate_version": phase1.ACID_CONFORMER_GATE_VERSION,
        "settings": asdict(settings),
        "manifest_path": manifest_path.name,
        "manifest_sha256": sha256_path(manifest_path),
        "candidates": candidates,
    }
    path = run_dir / "matched-candidates.json"
    atomic_json(path, payload)
    return path


def _refinement_summary(
    records: dict[str, dict[str, Any]],
    *,
    source: RefinementSource,
    al_record: dict[str, Any] | None,
) -> dict[str, Any]:
    if set(records) != _expected_si_keys():
        raise RuntimeError("effective refinement population is not exactly 16 Si seeds")
    si_records = [records[key] for key in sorted(records)]
    counts = {
        status: sum(record.get("status") == status for record in si_records)
        for status in sorted(TERMINAL_STATUSES)
    }
    accepted = [record for record in si_records if record.get("status") == "accepted"]
    matched = []
    if accepted:
        if len(accepted) != 1 or accepted[0].get("seed_key") != source.target_key:
            raise RuntimeError(
                "refinement produced an unexpected accepted Si population"
            )
        if al_record is None:
            verdict = "incomplete-matched-screen"
        elif al_record.get("status") == "accepted":
            _reaction, water_token, family = source.target_key.split(":", maxsplit=2)
            n_water = int(water_token.removesuffix("w"))
            matched = [
                {
                    "water_count": n_water,
                    "family": family,
                    "si_seed": source.target_key,
                    "al_seed": _seed_key("al-acid", n_water, family),
                }
            ]
            verdict = "matched-minimum-ready-for-barrier-ladder"
        elif al_record.get("status") in {"failed", "blocked", "rejected"}:
            verdict = "incomplete-matched-screen"
        else:
            raise RuntimeError("matched Al screen produced an invalid status")
    elif counts["rejected"] == len(si_records):
        verdict = "model-valid-no-go-pre-equilibrated-bridge"
    else:
        verdict = "incomplete-si-screen"
    return {
        "verdict": verdict,
        "si_seed_count": len(si_records),
        "si_status_counts": counts,
        "si_accepted_count": len(accepted),
        "si_unique_basin_count": len(accepted),
        "matched_candidates": matched,
        "matched_receipt": "matched-candidates.json" if matched else None,
        "source_summary": source.manifest.get("summary"),
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }


def _write_refinement_matched_receipt(
    attempt_dir: Path,
    resolution_path: Path,
    *,
    source: RefinementSource,
    si_record: dict[str, Any],
    al_record: dict[str, Any],
) -> Path:
    _reaction, water_token, family = source.target_key.split(":", maxsplit=2)
    n_water = int(water_token.removesuffix("w"))
    al_dir = _seed_dir(attempt_dir / "al-screen", "al-acid", n_water, family)
    models = {
        "si-acid": {
            "seed_key": source.target_key,
            "optimized_path": f"{REFINEMENT_ATTEMPT}/optimized.xyz",
            "optimized_sha256": si_record["artifacts"]["optimized.xyz"],
            "minimum_path": f"{REFINEMENT_ATTEMPT}/minimum.json",
            "minimum_sha256": si_record["artifacts"]["minimum.json"],
            "electronic_hartree": si_record["electronic_hartree"],
        },
        "al-acid": {
            "seed_key": _seed_key("al-acid", n_water, family),
            "optimized_path": str((al_dir / "optimized.xyz").relative_to(attempt_dir)),
            "optimized_sha256": al_record["artifacts"]["optimized.xyz"],
            "minimum_path": str((al_dir / "minimum.json").relative_to(attempt_dir)),
            "minimum_sha256": al_record["artifacts"]["minimum.json"],
            "electronic_hartree": al_record["electronic_hartree"],
        },
    }
    payload = {
        "schema_version": 1,
        "conformer_version": phase1.ACID_CONFORMER_VERSION,
        "gate_version": phase1.ACID_CONFORMER_GATE_VERSION,
        "settings": asdict(source.settings),
        "manifest_path": resolution_path.name,
        "manifest_sha256": sha256_path(resolution_path),
        "candidates": [{"water_count": n_water, "family": family, "models": models}],
    }
    path = attempt_dir / "matched-candidates.json"
    atomic_json(path, payload)
    phase1.validate_matched_ensemble_receipt(
        path,
        reaction="si-acid",
        n_water=n_water,
        family=family,
        settings=source.settings,
    )
    phase1.validate_matched_ensemble_receipt(
        path,
        reaction="al-acid",
        n_water=n_water,
        family=family,
        settings=source.settings,
    )
    return path


def _emergency_atomic_json(path: Path, payload: dict[str, Any]) -> None:
    """Independent last-ditch writer used after releasing reserved disk space."""
    temporary = path.with_name(f".{path.name}.emergency-{time.time_ns()}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    temporary.replace(path)


def _finalize_incomplete_resolution(
    resolution_path: Path,
    reserve_path: Path,
    resolution: dict[str, Any],
    exc: BaseException,
) -> dict[str, Any]:
    """Persist a spent-budget, non-success verdict after recoverable I/O failure."""
    reserve_path.unlink(missing_ok=True)
    finished_at = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    incomplete = {
        **resolution,
        "status": "incomplete-artifact-persistence",
        "attempt": {
            **resolution["attempt"],
            "status": "incomplete-artifact-persistence",
            "budget_spent": True,
            "finished_at": finished_at,
        },
        "summary": {
            "verdict": "incomplete-artifact-persistence",
            "completed_at": finished_at,
        },
        "finalization_error": f"{type(exc).__name__}: {exc}",
    }
    try:
        atomic_json(resolution_path, incomplete)
    except Exception:
        _emergency_atomic_json(resolution_path, incomplete)
    return incomplete


def run_failed_endpoint_refinement(
    source: RefinementSource,
    *,
    attempt_dir: Path,
    max_steps: int,
    log_path: str,
) -> dict[str, Any]:
    """Spend one bounded fresh optimizer attempt on the pinned failed endpoint."""
    max_steps = _require_refinement_bound(max_steps)
    attempt_dir = attempt_dir.expanduser().resolve()
    expected_attempt_dir = refinement_attempt_dir(source.manifest_path)
    if attempt_dir != expected_attempt_dir:
        raise RuntimeError(
            f"refinement attempt directory is fixed at {expected_attempt_dir}"
        )
    if attempt_dir.exists():
        raise RuntimeError(
            "refinement attempt directory already exists; refusing a second step budget"
        )
    attempt_dir.mkdir(parents=True)
    attempt_root = attempt_dir / REFINEMENT_ATTEMPT
    attempt_root.mkdir()
    resolution_path = attempt_dir / "refinement-manifest.json"
    receipt_path = attempt_root / "receipt.json"

    reaction, water_token, family = source.target_key.split(":", maxsplit=2)
    n_water = int(water_token.removesuffix("w"))
    source_seed_dir = _seed_dir(source.manifest_path.parent, reaction, n_water, family)
    source_endpoint = source_seed_dir / "optimized.xyz"
    input_path = attempt_root / "input-optimized.xyz"
    input_path.write_bytes(source_endpoint.read_bytes())
    source_endpoint_sha256 = source.target_record["artifacts"]["optimized.xyz"]
    if sha256_path(input_path) != source_endpoint_sha256:
        raise RuntimeError("copied refinement input is not byte-identical to source")
    input_cluster = load_xyz_strict(input_path, _template(reaction, n_water, family))

    resolution: dict[str, Any] = {
        "schema_version": REFINEMENT_SCHEMA_VERSION,
        "kind": "failed-endpoint-refinement",
        "status": "running",
        "source": {
            "manifest_path": str(source.manifest_path),
            "manifest_sha256": source.manifest_sha256,
            "log_path": str(source.log_path),
            "log_sha256": source.log_sha256,
        },
        "target_seed_key": source.target_key,
        "attempt": {
            "id": REFINEMENT_ATTEMPT,
            "status": "running",
            "input_path": str(input_path.relative_to(attempt_dir)),
            "input_sha256": source_endpoint_sha256,
            "optimizer": {
                "name": "geometric",
                "state_initialization": "fresh",
                "hessian_initialization": "fresh-default",
                "max_steps": max_steps,
            },
            "settings": asdict(source.settings),
            "log_path": log_path,
            "started_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        },
    }
    reserve_path = attempt_dir / ".finalization-reserve"
    reserve_path.write_bytes(b"\0" * FINALIZATION_RESERVE_BYTES)
    atomic_json(resolution_path, resolution)
    print(
        f"A1D_START seed={source.target_key} attempt={REFINEMENT_ATTEMPT} "
        f"max_steps={max_steps}",
        flush=True,
    )

    production_converged = False
    failure_kind = None
    try:
        result = optimize_bounded(input_cluster, source.settings, max_steps=max_steps)
        production_converged = result.converged
        output_path = attempt_root / "optimized.xyz"
        phase1.save_xyz(result.cluster, output_path)
        optimized = load_xyz_strict(output_path, input_cluster)
        if result.converged:
            try:
                record = evaluate_optimized_endpoint(
                    attempt_root,
                    optimized,
                    source.settings,
                    reaction=reaction,
                    n_water=n_water,
                    family=family,
                )
            except Exception as exc:
                failure_kind = "endpoint-evaluation-error"
                record = {
                    "status": "failed",
                    "reason": f"{type(exc).__name__}: {exc}",
                    "basin_signature": None,
                    "n_imaginary": None,
                    "electronic_hartree": None,
                    "minimum_pair_distance_a": minimum_pair_distance(optimized),
                }
        else:
            failure_kind = "optimizer-exhaustion"
            record = {
                "status": "failed",
                "reason": (
                    f"RuntimeError: geometry optimization did not converge within "
                    f"{max_steps} refinement steps"
                ),
                "basin_signature": None,
                "n_imaginary": None,
                "electronic_hartree": None,
                "minimum_pair_distance_a": minimum_pair_distance(optimized),
            }
    except Exception as exc:
        failure_kind = "optimizer-error"
        record = {
            "status": "failed",
            "reason": f"{type(exc).__name__}: {exc}",
            "basin_signature": None,
            "n_imaginary": None,
            "electronic_hartree": None,
            "minimum_pair_distance_a": None,
        }

    record.update(
        {
            "seed_key": source.target_key,
            "reaction": reaction,
            "water_count": n_water,
            "family": family,
            "preoptimization_converged": None,
            "production_converged": production_converged,
            "failure_kind": failure_kind,
            "refinement_attempt": REFINEMENT_ATTEMPT,
            "source_optimized_sha256": source_endpoint_sha256,
            "artifacts": _artifact_hashes(
                attempt_root,
                ("input-optimized.xyz", "optimized.xyz", "minimum.json"),
            ),
            "finished_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        }
    )
    try:
        assert_refinement_source_unchanged(source)
    except RuntimeError as exc:
        record["status"] = "blocked"
        record["reason"] = str(exc)

    receipt = {
        "schema_version": REFINEMENT_SCHEMA_VERSION,
        "kind": "failed-endpoint-refinement-attempt",
        "source_manifest_sha256": source.manifest_sha256,
        "source_seed_key": source.target_key,
        "input": {
            "path": "input-optimized.xyz",
            "sha256": source_endpoint_sha256,
        },
        "optimizer": resolution["attempt"]["optimizer"],
        "settings": asdict(source.settings),
        "record": record,
    }
    try:
        atomic_json(receipt_path, receipt)
        receipt_sha256 = sha256_path(receipt_path)

        effective_records = {
            key: dict(source_record)
            for key, source_record in source.manifest["seeds"].items()
        }
        effective_records[source.target_key] = record
        al_record = None
        al_manifest_path = None
        if record["status"] == "accepted":
            al_manifest_path = attempt_dir / "al-screen" / "manifest.json"
            al_manifest_path.parent.mkdir()
            al_manifest: dict[str, Any] = {"seeds": {}}
            try:
                al_record = screen_seed(
                    al_manifest,
                    al_manifest_path,
                    reaction="al-acid",
                    n_water=n_water,
                    family=family,
                    settings=source.settings,
                    max_steps=source.manifest["max_steps"],
                )
            except Exception as exc:
                al_manifest_path = None
                al_record = {
                    "status": "failed",
                    "reaction": "al-acid",
                    "water_count": n_water,
                    "family": family,
                    "reason": f"{type(exc).__name__}: {exc}",
                    "failure_kind": "matched-al-screen-error",
                }
            assert_refinement_source_unchanged(source)

        summary = _refinement_summary(
            effective_records, source=source, al_record=al_record
        )
        al_screen = None
        if al_record is not None:
            al_screen = {
                "manifest_path": (
                    None
                    if al_manifest_path is None
                    else str(al_manifest_path.relative_to(attempt_dir))
                ),
                "manifest_sha256": (
                    None if al_manifest_path is None else sha256_path(al_manifest_path)
                ),
                "record": al_record,
            }
        resolution.update(
            {
                "status": "evidence-complete",
                "attempt": {
                    **resolution["attempt"],
                    "status": record["status"],
                    "receipt_path": str(receipt_path.relative_to(attempt_dir)),
                    "receipt_sha256": receipt_sha256,
                    "finished_at": record["finished_at"],
                },
                "effective_seed_refs": {
                    key: (
                        {
                            "source": "a1d-refinement",
                            "receipt_path": str(receipt_path.relative_to(attempt_dir)),
                            "receipt_sha256": receipt_sha256,
                        }
                        if key == source.target_key
                        else {
                            "source": "a1c-manifest",
                            "manifest_sha256": source.manifest_sha256,
                            "record_key": key,
                        }
                    )
                    for key in sorted(effective_records)
                },
                "al_screen": al_screen,
                "summary": summary,
            }
        )
        evidence_path = attempt_dir / "evidence-manifest.json"
        atomic_json(evidence_path, resolution)
        matched_receipt = None
        if summary["matched_candidates"]:
            if al_record is None:
                raise RuntimeError("matched summary lacks an Al record")
            matched_path = _write_refinement_matched_receipt(
                attempt_dir,
                evidence_path,
                source=source,
                si_record=record,
                al_record=al_record,
            )
            matched_receipt = {
                "path": str(matched_path.relative_to(attempt_dir)),
                "sha256": sha256_path(matched_path),
            }
        assert_refinement_source_unchanged(source)
        final_resolution = {
            **resolution,
            "status": "completed",
            "evidence_manifest": {
                "path": evidence_path.name,
                "sha256": sha256_path(evidence_path),
            },
            "matched_receipt": matched_receipt,
        }
        atomic_json(resolution_path, final_resolution)
    except Exception as exc:
        incomplete = _finalize_incomplete_resolution(
            resolution_path, reserve_path, resolution, exc
        )
        print(
            "A1D_SUMMARY " + json.dumps(incomplete["summary"], sort_keys=True),
            flush=True,
        )
        return incomplete
    reserve_path.unlink(missing_ok=True)
    print("A1D_SUMMARY " + json.dumps(summary, sort_keys=True), flush=True)
    return final_resolution


def run_ensemble(
    run_dir: Path,
    settings: DftSettings,
    *,
    max_steps: int,
    log_path: str,
) -> dict:
    run_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = run_dir / "manifest.json"
    expected = _expected_manifest(settings, max_steps, log_path)
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text())
        for field in (
            "schema_version",
            "conformer_version",
            "gate_version",
            "settings",
            "preoptimization",
            "max_steps",
            "water_counts",
            "families",
        ):
            if manifest.get(field) != expected[field]:
                raise RuntimeError(f"existing ensemble manifest has different {field}")
        manifest["log"] = log_path
    else:
        manifest = expected
        atomic_json(manifest_path, manifest)

    for n_water in WATER_COUNTS:
        for family in ACID_MICROSOLVATION_FAMILIES:
            screen_seed(
                manifest,
                manifest_path,
                reaction="si-acid",
                n_water=n_water,
                family=family,
                settings=settings,
                max_steps=max_steps,
            )

    si_candidates = []
    for n_water in WATER_COUNTS:
        for family in ACID_MICROSOLVATION_FAMILIES:
            key = _seed_key("si-acid", n_water, family)
            record = manifest["seeds"][key]
            if record["status"] == "accepted":
                si_candidates.append(
                    _load_candidate(
                        run_dir,
                        reaction="si-acid",
                        n_water=n_water,
                        family=family,
                        record=record,
                        settings=settings,
                    )
                )
    assignments = deduplicate_basins(si_candidates)
    for label, canonical in assignments.items():
        manifest["seeds"][label]["basin_canonical"] = canonical
    atomic_json(manifest_path, manifest)

    canonical_si = [
        candidate
        for candidate in si_candidates
        if assignments[candidate.label] == candidate.label
    ]
    matched = []
    al_records = []
    for candidate in canonical_si:
        _reaction, water_token, family = candidate.label.split(":", maxsplit=2)
        n_water = int(water_token.removesuffix("w"))
        al_record = screen_seed(
            manifest,
            manifest_path,
            reaction="al-acid",
            n_water=n_water,
            family=family,
            settings=settings,
            max_steps=max_steps,
        )
        al_records.append(al_record)
        if al_record["status"] == "accepted":
            matched.append(
                {
                    "water_count": n_water,
                    "family": family,
                    "si_seed": candidate.label,
                    "al_seed": _seed_key("al-acid", n_water, family),
                    "barrier_commands": [
                        (
                            "python scripts/phase1_xiao_lasaga.py "
                            f"--reaction {reaction} --gpu "
                            f"--microsolvation-waters {n_water} "
                            f"--microsolvation-family {family} "
                            "--matched-ensemble-receipt "
                            f"{run_dir / 'matched-candidates.json'} "
                            "--threads 16 --nice 10"
                        )
                        for reaction in ("si-acid", "al-acid")
                    ],
                }
            )

    si_records = [
        manifest["seeds"][_seed_key("si-acid", n_water, family)]
        for n_water in WATER_COUNTS
        for family in ACID_MICROSOLVATION_FAMILIES
    ]
    si_failed = [record for record in si_records if record["status"] != "rejected"]
    if matched:
        verdict = "matched-minimum-ready-for-barrier-ladder"
    elif canonical_si and any(
        record["status"] in {"failed", "blocked"} for record in al_records
    ):
        verdict = "incomplete-matched-screen"
    elif canonical_si:
        verdict = "si-minimum-without-matched-al-minimum"
    elif not si_failed:
        verdict = "model-valid-no-go-pre-equilibrated-bridge"
    else:
        verdict = "incomplete-si-screen"
    manifest["summary"] = {
        "verdict": verdict,
        "si_seed_count": len(si_records),
        "si_status_counts": {
            status: sum(record["status"] == status for record in si_records)
            for status in sorted(TERMINAL_STATUSES)
        },
        "si_accepted_count": len(si_candidates),
        "si_unique_basin_count": len(canonical_si),
        "matched_candidates": matched,
        "matched_receipt": "matched-candidates.json" if matched else None,
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    atomic_json(manifest_path, manifest)
    matched_receipt = run_dir / "matched-candidates.json"
    if matched:
        write_matched_candidate_receipt(run_dir, manifest_path, manifest, settings)
    elif matched_receipt.exists():
        stale = matched_receipt.with_name(
            f"matched-candidates.stale-{time.time_ns()}.json"
        )
        matched_receipt.replace(stale)
    print("A1C_SUMMARY " + json.dumps(manifest["summary"], sort_keys=True), flush=True)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--nice", type=int, default=10)
    parser.add_argument("--log")
    parser.add_argument("--xc", default="b3lyp")
    parser.add_argument("--basis", default="def2-svp")
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--max-steps", type=int, default=100)
    parser.add_argument("--run-dir", type=Path)
    parser.add_argument("--refine-failed-endpoint", choices=(A1D_TARGET_SEED,))
    parser.add_argument(
        "--refinement-max-steps", type=int, default=REFINEMENT_MAX_STEPS
    )
    args = parser.parse_args()
    if args.refine_failed_endpoint:
        if (
            args.run_dir is not None
            or args.xc != "b3lyp"
            or args.basis != "def2-svp"
            or args.gpu
            or args.max_steps != 100
        ):
            parser.error(
                "refinement derives settings from the pinned source; do not pass "
                "--run-dir, --xc, --basis, --gpu, or --max-steps"
            )
        if _ETIQUETTE_SESSION is None:
            raise RuntimeError("campaign etiquette was not initialized")
        try:
            source = load_a1d_refinement_source(max_steps=args.refinement_max_steps)
            attempt_dir = refinement_attempt_dir(source.manifest_path)
            resolution = run_failed_endpoint_refinement(
                source,
                attempt_dir=attempt_dir,
                max_steps=args.refinement_max_steps,
                log_path=str(_ETIQUETTE_SESSION.log_path),
            )
        except (RuntimeError, ValueError) as exc:
            parser.error(str(exc))
        verdict = resolution["summary"]["verdict"]
        return (
            0
            if verdict
            in {
                "matched-minimum-ready-for-barrier-ladder",
                "model-valid-no-go-pre-equilibrated-bridge",
            }
            else 1
        )
    if args.max_steps <= 0 or args.max_steps > 240:
        parser.error("--max-steps must be in 1..240")
    settings = DftSettings(
        xc=args.xc,
        basis=args.basis,
        density_fit=True,
        use_gpu=args.gpu,
    )
    method_slug = f"{args.xc}-{args.basis}".lower().replace("/", "-")
    run_dir = args.run_dir or (
        Path(__file__).resolve().parent.parent
        / "runs"
        / "phase1"
        / (
            f"acid-microsolvation-ensemble-v{phase1.ACID_CONFORMER_VERSION}-"
            f"g{phase1.ACID_CONFORMER_GATE_VERSION}-{method_slug}"
        )
    )
    if _ETIQUETTE_SESSION is None:
        raise RuntimeError("campaign etiquette was not initialized")
    manifest = run_ensemble(
        run_dir,
        settings,
        max_steps=args.max_steps,
        log_path=str(_ETIQUETTE_SESSION.log_path),
    )
    verdict = manifest["summary"]["verdict"]
    return 1 if verdict.startswith("incomplete-") else 0


if __name__ == "__main__":
    raise SystemExit(main())
