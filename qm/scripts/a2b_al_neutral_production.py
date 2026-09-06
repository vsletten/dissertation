#!/usr/bin/env python3
"""A2b production re-tiering for the banked one-water Si-O-Al route.

This route-specific front end deliberately reuses the resume-safe A2 checkpoint
primitives.  It accepts only the frozen TASK-168 sequential mechanism, rejects
its two-imaginary-mode addition candidate, preserves physical hydrogen/topology
identity, and publishes terminal outputs atomically.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sqlite3
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    _ETIQUETTE = bootstrap_cli(
        "a2b-al-neutral-production",
        default_run_root=(
            "/mnt/data/vsletten/dissertation-data/"
            "task290-a2b-al-neutral-production-20260906"
        ),
        gpu_owner="a2b-al-neutral-production",
    )

import numpy as np

from quarry.clusters import Cluster
from quarry.pipeline import DftSettings, FrequencyResult
from quarry.rates import rate_from_thermo
from quarry.store import Store
from quarry.ts import find_ts, full_irc, reaction_path_vector
from scripts import production_energetics as a2

A2B_VERSION = "a2b-al-neutral-production-v1"
SVP_METHOD = "b3lyp/def2-svp/df"
ROLE_IDS = {"reactant": 1, "intermediate": 2, "transition-state": 3}
ROLE_NAMES = {
    "reactant": "al-neutral-complex",
    "intermediate": "al-neutral-intermediate",
    "transition-state": "al-neutral-cleavage-ts",
}
EXPECTED_BASINS = {
    "reactant": (False, True, False),
    "intermediate": (True, True, True),
    "product": (True, False, True),
}
SOURCE_ROOT = Path("/mnt/data/vsletten/dissertation-data/task168-al-neutral-20260823")


@dataclass(frozen=True)
class SourceContract:
    artifact_sha256: dict[str, str]
    role_ids: dict[str, int]
    role_names: dict[str, str]
    expected_basins: dict[str, tuple[bool, bool, bool]]


TS_XYZ_SHA256 = "40357de814b0d0113b1bc23e6194a083f15f5f93dfcfe55aa5391848eb095064"


ACCEPTED_SOURCE = SourceContract(
    artifact_sha256={
        "store.sqlite": (
            "d93aa48f4633b0f48d8e8a814246fca3adbc4aee179fbf3b835e91180f50dc88"
        ),
        "results.json": (
            "94ecd8ea45dbea17743d7ce9ad354beb4ed150ff97d988714e6f0dd895ccb148"
        ),
        "complex.xyz": (
            "b8b64d59d07cb0d27063de257bee912d2cdb1edbe4f1cec458215964892ede3e"
        ),
        "intermediate.xyz": (
            "db06060ef7b06cae63a26f7f8556af9b2549ac34a42903987971f69ff278b9a8"
        ),
        "ts.xyz": TS_XYZ_SHA256,
        "hydrolyzed_product.xyz": (
            "5510547d03c84095cd86533c6f82fdc84c865f9a7f95aeeff7da486cd46288a3"
        ),
    },
    role_ids=ROLE_IDS,
    role_names=ROLE_NAMES,
    expected_basins=EXPECTED_BASINS,
)


@dataclass(frozen=True)
class SourceEvidence:
    root: Path
    clusters: dict[str, Cluster]
    svp_energies: dict[str, float]
    receipt: dict[str, Any]
    contract: SourceContract


def method_slug(method: str) -> str:
    return method.split("/")[0].replace("(", "-").replace(")", "")


def method_contract(*, use_gpu: bool) -> dict[str, Any]:
    r2scan3c, production, b3lyp_d4 = a2.settings(use_gpu=use_gpu)
    methods = {
        a2.R2SCAN3C_METHOD: asdict(r2scan3c),
        a2.PRODUCTION_METHOD: asdict(production),
        a2.B3LYP_D4_METHOD: asdict(b3lyp_d4),
    }
    return {
        "a2_method_contract": "settled-a2-production-tier-v1",
        "methods": methods,
        "settings_fingerprints": {
            method: a2.frequency_settings_fingerprint(settings)
            for method, settings in zip(
                methods,
                (r2scan3c, production, b3lyp_d4),
                strict=True,
            )
        },
    }


def _read_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"{path}: expected a JSON object")
    return payload


def _require_artifact_hashes(root: Path, contract: SourceContract) -> None:
    for relative, expected in contract.artifact_sha256.items():
        path = root / relative
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"accepted source artifact missing or symlinked: {path}")
        actual = a2.sha256_path(path)
        if actual != expected:
            raise ValueError(
                f"accepted source artifact drift for {relative}: "
                f"expected {expected}, got {actual}"
            )


def _source_frequency_energies(
    store_path: Path, contract: SourceContract
) -> dict[str, float]:
    connection = sqlite3.connect(f"file:{store_path}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    try:
        integrity = connection.execute("PRAGMA integrity_check").fetchone()[0]
        if integrity != "ok":
            raise ValueError(f"accepted source store integrity failure: {integrity}")
        rows = connection.execute(
            """
            SELECT s.id AS structure_id, s.name, j.kind, j.method, j.status,
                   r.key, r.value, r.units
            FROM structures AS s
            JOIN jobs AS j ON j.structure_id = s.id
            JOIN results AS r ON r.job_id = j.id
            ORDER BY s.id, j.id, r.key
            """
        ).fetchall()
    finally:
        connection.close()
    if len(rows) != len(contract.role_ids):
        raise ValueError("accepted source store must contain exactly three energy rows")
    by_id = {int(row["structure_id"]): row for row in rows}
    energies: dict[str, float] = {}
    for role, structure_id in contract.role_ids.items():
        row = by_id.get(structure_id)
        if row is None:
            raise ValueError(f"accepted source store missing role {role}")
        expected = (
            contract.role_names[role],
            "freq",
            SVP_METHOD,
            "done",
            "electronic",
            "hartree",
        )
        actual = tuple(
            row[key] for key in ("name", "kind", "method", "status", "key", "units")
        )
        if actual != expected:
            raise ValueError(f"accepted source job drift for {role}: {actual!r}")
        value = float(row["value"])
        if not np.isfinite(value):
            raise ValueError(f"accepted source energy is non-finite for {role}")
        energies[role] = value
    return energies


def _require_source_mechanism(results: dict[str, Any]) -> None:
    required = {
        "reaction": "al-neutral",
        "method": SVP_METHOD,
        "mechanism": "sequential-associative-uphill-slide-no-resolved-addition-saddle",
    }
    for key, expected in required.items():
        if results.get(key) != expected:
            raise ValueError(f"accepted TASK-168 {key} drift")
    addition = results.get("addition_attempt", {})
    step = results.get("steps", {}).get("addition", {})
    if addition.get("status") != "no-saddle":
        raise ValueError("rejected addition candidate was reclassified as a saddle")
    if step.get("accepted_first_order_saddle") is not False:
        raise ValueError(
            "rejected two-imaginary-mode addition candidate cannot be used"
        )
    if step.get("classification") != "quasi-barrierless-uphill-slide":
        raise ValueError("accepted quasi-barrierless addition classification drift")
    if results.get("profile", {}).get("highest_profile_ts") != "cleavage":
        raise ValueError("accepted TASK-168 profile maximum is not cleavage")


def require_role_identities(
    clusters: dict[str, Cluster],
    reference: dict[str, Cluster],
    *,
    attacker_index: int,
    expected_basins: dict[str, tuple[bool, bool, bool]],
) -> dict[str, Any]:
    receipt: dict[str, Any] = {}
    for role in ("reactant", "intermediate", "product"):
        cluster = clusters[role]
        expected_cluster = reference[role]
        basin = a2.si_neutral_signature(cluster, attacker_index)
        if basin != expected_basins[role]:
            raise RuntimeError(
                f"{role}: expected basin {expected_basins[role]}, got {basin}"
            )
        actual_identity = a2.endpoint_identity(cluster, attacker_index)
        expected_identity = a2.endpoint_identity(expected_cluster, attacker_index)
        if actual_identity != expected_identity:
            raise RuntimeError(
                f"{role}: physical hydrogen ownership or heavy topology drift"
            )
        receipt[role] = {
            "basin": list(basin),
            "hydrogen_owners": [list(value) for value in actual_identity[1]],
            "heavy_atom_bonds": [list(value) for value in actual_identity[2]],
        }
    return receipt


def load_source_evidence(
    root: Path,
    *,
    attacker_index: int,
    contract: SourceContract = ACCEPTED_SOURCE,
) -> SourceEvidence:
    root = root.resolve()
    _require_artifact_hashes(root, contract)
    results = _read_json(root / "results.json")
    _require_source_mechanism(results)
    svp_energies = _source_frequency_energies(root / "store.sqlite", contract)
    clusters: dict[str, Cluster] = {}
    structure_receipts: dict[str, Any] = {}
    for role, structure_id in contract.role_ids.items():
        cluster, receipt = a2.load_store_structure(root / "store.sqlite", structure_id)
        if receipt["name"] != contract.role_names[role]:
            raise ValueError(f"accepted source structure-name drift for {role}")
        clusters[role] = cluster
        structure_receipts[role] = receipt
    reactant = clusters["reactant"]
    clusters["product"] = a2.load_xyz_like(
        root / "hydrolyzed_product.xyz",
        reactant,
        name="al-neutral-product",
    )
    if not all(
        cluster.symbols == reactant.symbols
        and cluster.charge == reactant.charge
        and cluster.spin == reactant.spin
        for cluster in clusters.values()
    ):
        raise ValueError("accepted source atom order, charge, or spin drift")
    identity = require_role_identities(
        clusters,
        clusters,
        attacker_index=attacker_index,
        expected_basins=contract.expected_basins,
    )
    receipt = {
        "source_contract": "task168-al-neutral-accepted-sequential-v1",
        "root": str(root),
        "artifact_sha256": dict(contract.artifact_sha256),
        "structures": structure_receipts,
        "source_method": SVP_METHOD,
        "source_electronic_hartree": svp_energies,
        "mechanism": results["mechanism"],
        "rejected_addition_saddle_used": False,
        "accepted_cleavage_transition_state": True,
        "typed_identity": identity,
    }
    return SourceEvidence(root, clusters, svp_energies, receipt, contract)


def assert_source_unchanged(source: SourceEvidence) -> None:
    _require_artifact_hashes(source.root, source.contract)


def require_sequential_irc(
    endpoints: tuple[Cluster, Cluster],
    intermediate: Cluster,
    product: Cluster,
    *,
    attacker_index: int,
) -> dict[str, Any]:
    expected = {
        a2.endpoint_identity(intermediate, attacker_index),
        a2.endpoint_identity(product, attacker_index),
    }
    actual = {a2.endpoint_identity(cluster, attacker_index) for cluster in endpoints}
    if actual != expected:
        raise RuntimeError(
            "cleavage IRC endpoints do not connect the exact accepted intermediate "
            "and product basins with physical-H ownership and heavy topology"
        )
    basins = {a2.si_neutral_signature(cluster, attacker_index) for cluster in endpoints}
    if basins != {EXPECTED_BASINS["intermediate"], EXPECTED_BASINS["product"]}:
        raise RuntimeError(f"cleavage IRC basin mismatch: {sorted(basins)}")
    return {
        "endpoint_basins": [list(value) for value in sorted(basins)],
        "typed_identity": "basin+physical-hydrogen-owners+heavy-atom-bonds",
    }


def significant_imaginary_count(freq: FrequencyResult, floor_cm: float) -> int:
    return int(np.count_nonzero(np.asarray(freq.imaginary_cm) > floor_cm))


def record_store(
    path: Path,
    clusters: dict[str, Cluster],
    source_clusters: dict[str, Cluster],
    frequencies: dict[str, FrequencyResult],
    energies: dict[str, dict[str, float]],
    summary: dict[str, Any],
    *,
    use_gpu: bool,
    source_store_sha256: str,
) -> None:
    temporary = path.with_name(f".{path.name}.{time.time_ns()}.tmp")
    temporary.unlink(missing_ok=True)
    with Store(temporary) as store:
        structure_ids = {
            role: store.add_structure(
                f"al-neutral-{role}",
                cluster.formula,
                a2.exact_xyz(cluster),
                charge=cluster.charge,
                spin=cluster.spin,
            )
            for role, cluster in clusters.items()
        }
        source_structure_ids = {
            role: store.add_structure(
                f"al-neutral-source-svp-{role}",
                source_clusters[role].formula,
                a2.exact_xyz(source_clusters[role]),
                charge=source_clusters[role].charge,
                spin=source_clusters[role].spin,
            )
            for role in ("reactant", "intermediate", "transition-state")
        }
        for role, freq in frequencies.items():
            job = store.add_job(
                structure_ids[role],
                "freq",
                a2.R2SCAN3C_METHOD,
                "gpu4pyscf" if use_gpu else "pyscf",
                detail=json.dumps({"a2b_version": A2B_VERSION}, sort_keys=True),
            )
            store.set_job_status(job, "done")
            store.add_result(job, "electronic", freq.electronic_hartree, "hartree")
            store.add_result(
                job,
                "imaginary_count_above_floor",
                int(summary["stationary_points"][role]["significant_imaginary_count"]),
                "count",
            )
        for method, values in energies.items():
            # r2SCAN-3c electronic energies are already persisted on the
            # frequency jobs above; do not duplicate them as synthetic SP jobs.
            if method == a2.R2SCAN3C_METHOD:
                continue
            for role, value in values.items():
                engine = (
                    "source-task168"
                    if method == SVP_METHOD
                    else ("gpu4pyscf" if use_gpu else "pyscf")
                )
                detail = (
                    {
                        "a2b_version": A2B_VERSION,
                        "source_store_sha256": source_store_sha256,
                    }
                    if method == SVP_METHOD
                    else {"a2b_version": A2B_VERSION}
                )
                job = store.add_job(
                    (
                        source_structure_ids[role]
                        if method == SVP_METHOD
                        else structure_ids[role]
                    ),
                    "sp",
                    method,
                    engine,
                    detail=json.dumps(detail, sort_keys=True),
                )
                store.set_job_status(job, "done")
                store.add_result(job, "electronic", value, "hartree")
        analysis = store.add_job(
            structure_ids["transition-state"],
            "analysis",
            a2.PRODUCTION_METHOD,
            "quarry",
            detail=json.dumps(
                {"a2b_version": A2B_VERSION, "temperature_k": a2.TEMPERATURE_K},
                sort_keys=True,
            ),
        )
        store.set_job_status(analysis, "done")
        for method, barrier in summary["electronic_barriers_kj_mol"].items():
            store.add_result(
                analysis,
                f"electronic_barrier::{method}",
                float(barrier),
                "kJ/mol",
            )
        store.add_result(
            analysis,
            "production_dg_kj_mol",
            float(summary["production_thermochemistry"]["overall_dg_dagger_kj_mol"]),
            "kJ/mol",
        )
        store.add_result(
            analysis,
            "production_cleavage_dg_kj_mol",
            float(summary["production_thermochemistry"]["cleavage_dg_dagger_kj_mol"]),
            "kJ/mol",
        )
    connection = sqlite3.connect(temporary)
    try:
        integrity = connection.execute("PRAGMA integrity_check").fetchone()[0]
    finally:
        connection.close()
    if integrity != "ok":
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"published store integrity failure: {integrity}")
    temporary.replace(path)


def quarantine_terminal_outputs(run_dir: Path) -> Path | None:
    names = ("results.json", "store.sqlite", "terminal-receipt.json")
    existing = [run_dir / name for name in names if (run_dir / name).exists()]
    if not existing:
        return None
    quarantine = run_dir / "quarantine" / f"stale-{time.time_ns()}"
    quarantine.mkdir(parents=True, exist_ok=False)
    for path in existing:
        shutil.move(str(path), quarantine / path.name)
    return quarantine


def run(args: argparse.Namespace) -> int:
    run_dir = args.run_dir.resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    source_root = args.source_root.resolve()
    if run_dir == source_root or source_root in run_dir.parents:
        raise ValueError("run directory cannot be inside the immutable source archive")
    source = load_source_evidence(
        source_root,
        attacker_index=args.attacker_index,
    )
    r2scan3c, production, b3lyp_d4 = a2.settings(use_gpu=args.gpu)
    contract = method_contract(use_gpu=args.gpu)
    source_receipt = {
        **source.receipt,
        "driver": str(Path(__file__).resolve()),
        "driver_sha256": a2.sha256_path(Path(__file__).resolve()),
        "a2_method_contract": contract,
    }
    a2.atomic_json(run_dir / "source-receipt.json", source_receipt)
    a2.atomic_json(
        run_dir / "settings.json",
        {
            "a2b_version": A2B_VERSION,
            **contract,
            "temperature_k": a2.TEMPERATURE_K,
            "imaginary_floor_cm": args.imaginary_floor,
            "bounds": {
                "minimum_steps": args.minimum_steps,
                "saddle_steps": args.saddle_steps,
                "irc_steps": args.irc_steps,
                "threads": args.threads,
                "gpu_memory_gb": args.gpu_mem_gb,
            },
        },
    )

    r2scan_identity = a2.frequency_settings_fingerprint(r2scan3c)
    clusters: dict[str, Cluster] = {}
    for role in ("reactant", "intermediate", "product"):
        template = source.clusters[role]
        clusters[role] = a2.checkpoint_cluster(
            run_dir / f"{role}.r2scan3c.xyz",
            template,
            lambda template=template, role=role: a2.optimize_minimum(
                template,
                r2scan3c,
                max_steps=args.minimum_steps,
                trajectory=run_dir / f"{role}.r2scan3c.traj",
            ),
            identity={
                "a2b_version": A2B_VERSION,
                "stage": "minimum",
                "role": role,
                "algorithm": "ase-bfgs-v1",
                "settings": r2scan_identity,
                "source_sha256": source.contract.artifact_sha256,
                "max_steps": args.minimum_steps,
            },
        )
    clusters["transition-state"] = a2.checkpoint_cluster(
        run_dir / "transition-state.r2scan3c.xyz",
        source.clusters["transition-state"],
        lambda: find_ts(
            source.clusters["transition-state"],
            r2scan3c,
            max_steps=args.saddle_steps,
            trajectory=str(run_dir / "transition-state.r2scan3c.traj"),
            initial_mode=reaction_path_vector(
                clusters["intermediate"], clusters["product"]
            ),
            internal=False,
        ),
        identity={
            "a2b_version": A2B_VERSION,
            "stage": "cleavage-transition-state",
            "algorithm": "sella-directed-cartesian-v1",
            "settings": r2scan_identity,
            "source_sha256": source.contract.artifact_sha256,
            "max_steps": args.saddle_steps,
            "intermediate_geometry": a2.frequency_geometry_fingerprint(
                clusters["intermediate"]
            ),
            "product_geometry": a2.frequency_geometry_fingerprint(clusters["product"]),
            "rejected_addition_saddle_used": False,
        },
    )
    typed_identity = require_role_identities(
        clusters,
        source.clusters,
        attacker_index=args.attacker_index,
        expected_basins=source.contract.expected_basins,
    )

    frequencies: dict[str, FrequencyResult] = {}
    stationary: dict[str, Any] = {}
    expected_indices = {
        "reactant": 0,
        "intermediate": 0,
        "product": 0,
        "transition-state": 1,
    }
    for role, cluster in clusters.items():
        freq = a2.checkpoint_frequency(
            run_dir / f"{role}.r2scan3c.frequency.json",
            cluster,
            r2scan3c,
            finite_difference=True,
        )
        count = significant_imaginary_count(freq, args.imaginary_floor)
        if count != expected_indices[role]:
            raise RuntimeError(
                f"r2SCAN-3c {role} has index {count} above "
                f"{args.imaginary_floor} cm^-1, expected {expected_indices[role]}"
            )
        frequencies[role] = freq
        stationary[role] = {
            "significant_imaginary_count": count,
            "all_imaginary_cm": [float(value) for value in freq.imaginary_cm],
            "geometry_fingerprint": a2.frequency_geometry_fingerprint(cluster),
        }
    significant = np.asarray(frequencies["transition-state"].imaginary_cm)
    significant = significant[significant > args.imaginary_floor]

    endpoints = a2.checkpoint_irc_endpoints(
        run_dir,
        clusters["transition-state"],
        lambda: full_irc(
            clusters["transition-state"],
            r2scan3c,
            max_steps=args.irc_steps,
            trajectory=run_dir / "full-irc.r2scan3c.traj",
            logfile=run_dir / "full-irc.r2scan3c.log",
        ),
        identity={
            "a2b_version": A2B_VERSION,
            "stage": "full-cleavage-irc",
            "algorithm": "sella-gonzalez-schlegel-full-irc-v1",
            "settings": r2scan_identity,
            "max_steps": args.irc_steps,
            "rejected_addition_saddle_used": False,
        },
    )
    irc_receipt = require_sequential_irc(
        endpoints,
        clusters["intermediate"],
        clusters["product"],
        attacker_index=args.attacker_index,
    )

    settings_by_method: dict[str, DftSettings] = {
        a2.PRODUCTION_METHOD: production,
        a2.B3LYP_D4_METHOD: b3lyp_d4,
    }
    energies: dict[str, dict[str, float]] = {SVP_METHOD: dict(source.svp_energies)}
    for method, method_settings in settings_by_method.items():
        slug = method_slug(method)
        energies[method] = {}
        for role, cluster in clusters.items():
            energies[method][role] = a2.checkpoint_converged_energy(
                run_dir / f"{role}.{slug}.energy.json",
                cluster,
                method_settings,
                method,
            )
    r2scan_energies = {
        role: float(freq.electronic_hartree) for role, freq in frequencies.items()
    }
    energies[a2.R2SCAN3C_METHOD] = r2scan_energies

    barriers = {
        method: (values["transition-state"] - values["reactant"]) * a2.HARTREE_TO_KJ
        for method, values in energies.items()
    }
    production_thermo = {
        role: a2.thermo(freq, energies[a2.PRODUCTION_METHOD][role])
        for role, freq in frequencies.items()
    }
    rate = rate_from_thermo(
        production_thermo["reactant"],
        production_thermo["transition-state"],
        imag_nu_cm=float(significant[0]),
        tunneling="wigner",
    )
    profile = {
        role: (value - energies[a2.PRODUCTION_METHOD]["reactant"]) * a2.HARTREE_TO_KJ
        for role, value in energies[a2.PRODUCTION_METHOD].items()
    }
    profile_maximum_role = max(profile, key=profile.__getitem__)
    if profile_maximum_role != "transition-state":
        raise RuntimeError(
            f"production profile maximum moved to {profile_maximum_role}; "
            "accepted cleavage saddle no longer rate-limiting"
        )
    summary: dict[str, Any] = {
        "a2b_version": A2B_VERSION,
        "reaction": "al-neutral",
        "mechanism": "sequential-associative-uphill-slide-cleavage-saddle",
        "rejected_addition_saddle_used": False,
        "temperature_k": a2.TEMPERATURE_K,
        "imaginary_floor_cm": args.imaginary_floor,
        "gpu": args.gpu,
        "stationary_points": stationary,
        "typed_identity": typed_identity,
        "irc": irc_receipt,
        "electronic_hartree": energies,
        "electronic_barriers_kj_mol": barriers,
        "barrier_shifts_vs_svp_kj_mol": {
            method: barrier - barriers[SVP_METHOD]
            for method, barrier in barriers.items()
            if method != SVP_METHOD
        },
        "production_profile_relative_kj_mol": profile,
        "production_profile_maximum_role": profile_maximum_role,
        "production_thermochemistry": {
            "overall_dg_dagger_kj_mol": rate.dg_kj,
            "overall_dh_dagger_kj_mol": rate.dh_kj,
            "overall_ds_dagger_kj_mol_k": rate.ds_kj_per_k,
            "cleavage_dg_dagger_kj_mol": (
                production_thermo["transition-state"].gibbs
                - production_thermo["intermediate"].gibbs
            ),
            "rate_s^-1": rate.k,
            "wigner_kappa": rate.kappa,
            "transition_state_imaginary_cm": float(significant[0]),
        },
        "source_receipt_sha256": a2.sha256_path(run_dir / "source-receipt.json"),
        "settings_sha256": a2.sha256_path(run_dir / "settings.json"),
        "source_unchanged_after_compute": True,
    }
    assert_source_unchanged(source)
    a2.atomic_json(run_dir / "results.json", summary)
    record_store(
        run_dir / "store.sqlite",
        clusters,
        source.clusters,
        frequencies,
        energies,
        summary,
        use_gpu=args.gpu,
        source_store_sha256=source.contract.artifact_sha256["store.sqlite"],
    )
    return 0


def _terminal_payload(
    run_dir: Path,
    *,
    outcome: str,
    error: BaseException | None = None,
) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "a2b_version": A2B_VERSION,
        "outcome": outcome,
        "finished_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "driver_sha256": a2.sha256_path(Path(__file__).resolve()),
        "run_status_sha256": a2.sha256_path(run_dir / "run_status.json"),
        "running_record_present": False,
    }
    if error is None:
        payload["artifacts"] = {
            name: a2.sha256_path(run_dir / name)
            for name in (
                "results.json",
                "store.sqlite",
                "source-receipt.json",
                "settings.json",
            )
        }
    else:
        payload.update({"error_type": type(error).__name__, "error": str(error)})
    return payload


def execute_with_status(args: argparse.Namespace) -> int:
    run_dir = args.run_dir.resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    quarantine = quarantine_terminal_outputs(run_dir)
    running = {
        "a2b_version": A2B_VERSION,
        "status": "running",
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "quarantine": str(quarantine) if quarantine else None,
    }
    a2.atomic_json(run_dir / "run_status.json", running)
    try:
        code = run(args)
    except BaseException as exc:
        quarantine_terminal_outputs(run_dir)
        failed = {
            **running,
            "status": "failed",
            "outcome": "incomplete-computational-failure",
            "finished_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "error_type": type(exc).__name__,
            "error": str(exc),
        }
        a2.atomic_json(run_dir / "run_status.json", failed)
        a2.atomic_json(
            run_dir / "terminal-receipt.json",
            _terminal_payload(
                run_dir,
                outcome="incomplete-computational-failure",
                error=exc,
            ),
        )
        raise
    completed = {
        **running,
        "status": "completed",
        "outcome": "success",
        "finished_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "results_sha256": a2.sha256_path(run_dir / "results.json"),
        "store_sha256": a2.sha256_path(run_dir / "store.sqlite"),
    }
    a2.atomic_json(run_dir / "run_status.json", completed)
    a2.atomic_json(
        run_dir / "terminal-receipt.json",
        _terminal_payload(run_dir, outcome="success"),
    )
    return code


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--source-root", type=Path, default=SOURCE_ROOT)
    result.add_argument("--run-dir", type=Path, required=True)
    result.add_argument("--attacker-index", type=int, default=15)
    result.add_argument("--gpu", action="store_true")
    result.add_argument("--gpu-mem-gb", type=float, default=18.0)
    result.add_argument("--minimum-steps", type=int, default=200)
    result.add_argument("--saddle-steps", type=int, default=400)
    result.add_argument("--irc-steps", type=int, default=400)
    result.add_argument("--imaginary-floor", type=float, default=30.0)
    result.add_argument("--threads", type=int, default=16)
    result.add_argument("--nice", type=int, default=10)
    result.add_argument("--log", type=Path)
    return result


def main() -> int:
    args = parser().parse_args()
    if args.minimum_steps <= 0 or args.saddle_steps <= 0 or args.irc_steps <= 0:
        raise ValueError("optimization and IRC bounds must be positive")
    if args.attacker_index < 0:
        raise ValueError("attacker-index must be nonnegative")
    if args.imaginary_floor <= 0:
        raise ValueError("imaginary-floor must be positive")
    if args.gpu_mem_gb <= 0:
        raise ValueError("gpu-mem-gb must be positive")
    return execute_with_status(args)


if __name__ == "__main__":
    raise SystemExit(main())
