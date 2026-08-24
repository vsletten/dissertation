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

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

_ETIQUETTE_SESSION = None
if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    _ETIQUETTE_SESSION = bootstrap_cli(
        "acid_microsolvation_ensemble",
        default_run_root=Path(__file__).resolve().parent.parent / "runs",
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
    optimize,
    optimize_bounded,
)
from scripts import phase1_xiao_lasaga as phase1  # noqa: E402

SCHEMA_VERSION = 1
WATER_COUNTS = (3, 4, 5, 6)
TERMINAL_STATUSES = frozenset({"accepted", "rejected", "failed", "blocked"})
MIN_PAIR_DISTANCE_A = 0.85
HF_PREOPT_MAX_STEPS = 40


@dataclass(frozen=True)
class BasinCandidate:
    label: str
    n_water: int
    cluster: Cluster
    electronic_hartree: float


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
        cluster = phase1.load_xyz(seed_dir / "optimized.xyz", template)
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
        try:
            preopt_result = optimize_bounded(
                seed,
                DftSettings(xc="hf", basis="sto-3g"),
                max_steps=HF_PREOPT_MAX_STEPS,
            )
            preoptimization_converged = preopt_result.converged
            preopt = preopt_result.cluster
            phase1.save_xyz(preopt, seed_dir / "preoptimized.xyz")
            optimized = optimize(preopt, settings, max_steps=max_steps)
            optimized_path = seed_dir / "optimized.xyz"
            phase1.save_xyz(optimized, optimized_path)
            optimized = phase1.load_xyz(optimized_path, preopt)
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
                        "frequency result failed geometry/settings/Hessian "
                        "receipt validation"
                    )
                n_imaginary_value = minimum_receipt.get("n_imaginary")
                energy_value = minimum_receipt.get("electronic_hartree")
                if type(n_imaginary_value) is not int or not isinstance(
                    energy_value, (int, float)
                ):
                    raise RuntimeError(
                        "validated minimum receipt has invalid numeric types"
                    )
                n_imaginary = n_imaginary_value
                electronic_hartree = float(energy_value)
                if n_imaginary:
                    reason = (
                        f"optimized protonated bridge has {n_imaginary} imaginary modes"
                    )
            status = "accepted" if reason is None else "rejected"
            record = {
                "status": status,
                "reaction": reaction,
                "water_count": n_water,
                "family": family,
                "reason": reason,
                "preoptimization_converged": preoptimization_converged,
                "basin_signature": list(signature),
                "n_imaginary": n_imaginary,
                "electronic_hartree": electronic_hartree,
                "minimum_pair_distance_a": minimum_pair_distance(optimized),
                "artifacts": _artifact_hashes(
                    seed_dir,
                    ("seed.xyz", "preoptimized.xyz", "optimized.xyz", "minimum.json"),
                ),
            }
        except Exception as exc:  # every bounded seed keeps its own failure receipt
            record = {
                "status": "failed",
                "reaction": reaction,
                "water_count": n_water,
                "family": family,
                "reason": f"{type(exc).__name__}: {exc}",
                "preoptimization_converged": preoptimization_converged,
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
    cluster = phase1.load_xyz(seed_dir / "optimized.xyz", template)
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
    args = parser.parse_args()
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
