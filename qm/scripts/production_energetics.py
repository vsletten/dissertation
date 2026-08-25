#!/usr/bin/env python3
"""Resume-safe A2 production re-tiering for one banked reaction pair.

The driver consumes immutable stationary-point coordinates from a quarry SQLite
store, optimizes and frequency-verifies them at exact r2SCAN-3c, traces the full
Gonzalez--Schlegel IRC, and applies the production single-point tier. Heavy
artifacts live below ``--run-dir`` (normally dissertation-data), never in Git.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sqlite3
import time
from collections.abc import Callable
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    _ETIQUETTE = bootstrap_cli(
        "production-energetics",
        default_run_root="/mnt/data/vsletten/dissertation-data/task207-a2-production",
    )

import numpy as np

from quarry.clusters import Cluster
from quarry.pipeline import (
    HARTREE_TO_KJ,
    DftSettings,
    FrequencyResult,
    energy,
    frequencies,
    frequency_geometry_fingerprint,
    frequency_settings_fingerprint,
    optimize,
)
from quarry.rates import rate_from_thermo, thermo_from_frequencies
from quarry.store import Store, geometry_hash
from quarry.ts import find_ts, full_irc

R2SCAN3C_METHOD = "r2scan-3c/def2-mtzvpp/d4/gcp"
PRODUCTION_METHOD = "wb97m-v/def2-tzvpd/smd(water)"
B3LYP_D4_METHOD = "b3lyp-d4/def2-tzvpd/smd(water)"
TEMPERATURE_K = 298.15


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_xyz(text: str, *, name: str, charge: int, spin: int) -> Cluster:
    lines = text.splitlines()
    if len(lines) < 2:
        raise ValueError(f"{name}: invalid XYZ header")
    count = int(lines[0])
    if len(lines) != count + 2:
        raise ValueError(f"{name}: XYZ line count drift")
    symbols: list[str] = []
    coords: list[list[float]] = []
    for line in lines[2:]:
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"{name}: malformed XYZ row")
        symbols.append(fields[0])
        coords.append([float(value) for value in fields[1:4]])
    array = np.asarray(coords, dtype=float)
    if array.shape != (count, 3) or not np.all(np.isfinite(array)):
        raise ValueError(f"{name}: non-finite or malformed coordinates")
    return Cluster(name=name, symbols=symbols, coords=array, charge=charge, spin=spin)


def load_store_structure(
    path: Path, structure_id: int
) -> tuple[Cluster, dict[str, Any]]:
    connection = sqlite3.connect(path)
    connection.row_factory = sqlite3.Row
    try:
        row = connection.execute(
            "SELECT id,name,formula,charge,spin,xyz,geometry_hash "
            "FROM structures WHERE id = ?",
            (structure_id,),
        ).fetchone()
    finally:
        connection.close()
    if row is None:
        raise KeyError(f"{path}: no structure id {structure_id}")
    payload = dict(row)
    xyz = str(payload.pop("xyz"))
    if geometry_hash(xyz) != payload["geometry_hash"]:
        raise ValueError(f"{path}: structure {structure_id} hash mismatch")
    cluster = parse_xyz(
        xyz,
        name=str(payload["name"]),
        charge=int(payload["charge"]),
        spin=int(payload["spin"]),
    )
    if cluster.formula != payload["formula"]:
        raise ValueError(f"{path}: structure {structure_id} formula mismatch")
    return cluster, payload


def load_xyz_like(path: Path, template: Cluster, *, name: str) -> Cluster:
    candidate = parse_xyz(
        path.read_text(),
        name=name,
        charge=template.charge,
        spin=template.spin,
    )
    if candidate.symbols != template.symbols:
        raise ValueError(f"{path}: atom order differs from reactant")
    return replace(candidate, frozen_indices=list(template.frozen_indices))


def checkpoint_cluster(
    path: Path,
    template: Cluster,
    compute: Callable[[], Cluster],
) -> Cluster:
    if path.exists():
        return load_xyz_like(path, template, name=path.stem)
    result = compute()
    if result.symbols != template.symbols:
        raise ValueError(f"{path.name}: computed atom order drift")
    temporary = path.with_name(f".{path.name}.{time.time_ns()}.tmp")
    temporary.write_text(result.to_xyz())
    temporary.replace(path)
    return result


def frequency_payload(result: FrequencyResult) -> dict[str, Any]:
    return {
        "frequencies_cm": result.frequencies_cm.tolist(),
        "imaginary_cm": result.imaginary_cm.tolist(),
        "electronic_hartree": result.electronic_hartree,
        "molar_mass_kg": result.molar_mass_kg,
        "rotational_temperatures_k": (
            list(result.rotational_temperatures_k)
            if result.rotational_temperatures_k is not None
            else None
        ),
        "linear": result.linear,
        "geometry_fingerprint": result.geometry_fingerprint,
        "settings_fingerprint": result.settings_fingerprint,
    }


def frequency_from_payload(payload: dict[str, Any]) -> FrequencyResult:
    rotational = payload["rotational_temperatures_k"]
    return FrequencyResult(
        frequencies_cm=np.asarray(payload["frequencies_cm"], dtype=float),
        imaginary_cm=np.asarray(payload["imaginary_cm"], dtype=float),
        electronic_hartree=float(payload["electronic_hartree"]),
        molar_mass_kg=float(payload["molar_mass_kg"]),
        rotational_temperatures_k=(
            tuple(float(value) for value in rotational)
            if rotational is not None
            else None
        ),
        linear=bool(payload["linear"]),
        geometry_fingerprint=str(payload["geometry_fingerprint"]),
        settings_fingerprint=str(payload["settings_fingerprint"]),
    )


def checkpoint_frequency(
    path: Path,
    cluster: Cluster,
    settings: DftSettings,
) -> FrequencyResult:
    expected_geometry = frequency_geometry_fingerprint(cluster)
    expected_settings = frequency_settings_fingerprint(settings)
    if path.exists():
        payload = json.loads(path.read_text())
        if payload.get("geometry_fingerprint") != expected_geometry:
            raise ValueError(f"{path.name}: cached geometry fingerprint drift")
        if payload.get("settings_fingerprint") != expected_settings:
            raise ValueError(f"{path.name}: cached settings fingerprint drift")
        return frequency_from_payload(payload)
    result = frequencies(cluster, settings)
    atomic_json(path, frequency_payload(result))
    return result


def checkpoint_energy(
    path: Path,
    cluster: Cluster,
    settings: DftSettings,
    method: str,
) -> float:
    expected_geometry = frequency_geometry_fingerprint(cluster)
    expected_settings = frequency_settings_fingerprint(settings)
    if path.exists():
        payload = json.loads(path.read_text())
        if payload.get("geometry_fingerprint") != expected_geometry:
            raise ValueError(f"{path.name}: cached geometry fingerprint drift")
        if payload.get("settings_fingerprint") != expected_settings:
            raise ValueError(f"{path.name}: cached settings fingerprint drift")
        return float(payload["electronic_hartree"])
    value = energy(cluster, settings)
    atomic_json(
        path,
        {
            "method": method,
            "electronic_hartree": value,
            "geometry_fingerprint": expected_geometry,
            "settings_fingerprint": expected_settings,
        },
    )
    return value


def thermo(freq: FrequencyResult, electronic_hartree: float):
    return thermo_from_frequencies(
        electronic_hartree * HARTREE_TO_KJ,
        freq.frequencies_cm,
        TEMPERATURE_K,
        molar_mass_kg=freq.molar_mass_kg,
        rotational_temperatures_k=(
            list(freq.rotational_temperatures_k)
            if freq.rotational_temperatures_k is not None
            else None
        ),
        linear=freq.linear,
    )


def settings(*, use_gpu: bool) -> tuple[DftSettings, DftSettings, DftSettings]:
    r2scan3c = DftSettings(
        xc="r2scan",
        basis="def2-mtzvpp",
        composite="r2scan3c",
        density_fit=True,
        use_gpu=use_gpu,
    )
    production = DftSettings(
        xc="wb97m-v",
        basis="def2-tzvpd",
        solvent="smd",
        density_fit=True,
        use_gpu=use_gpu,
    )
    b3lyp_d4 = DftSettings(
        xc="b3lyp",
        basis="def2-tzvpd",
        solvent="smd",
        dispersion="d4",
        density_fit=True,
        use_gpu=use_gpu,
    )
    return r2scan3c, production, b3lyp_d4


def si_neutral_signature(
    cluster: Cluster, attacker_index: int
) -> tuple[bool, bool, bool]:
    # Kept local so the production driver is not coupled to an argparse-heavy
    # historical campaign module.
    si_index = 1
    bridge_index = 0
    proton_indices = (attacker_index + 1, attacker_index + 2)
    return (
        float(np.linalg.norm(cluster.coords[si_index] - cluster.coords[attacker_index]))
        < 2.3,
        float(np.linalg.norm(cluster.coords[si_index] - cluster.coords[bridge_index]))
        < 2.3,
        min(
            float(np.linalg.norm(cluster.coords[bridge_index] - cluster.coords[proton]))
            for proton in proton_indices
        )
        < 1.25,
    )


def require_si_neutral_irc(
    endpoints: tuple[Cluster, Cluster],
    reactant: Cluster,
    product: Cluster,
    attacker_index: int,
) -> dict[str, Any]:
    expected = {
        si_neutral_signature(reactant, attacker_index),
        si_neutral_signature(product, attacker_index),
    }
    actual = {si_neutral_signature(endpoint, attacker_index) for endpoint in endpoints}
    if len(expected) != 2 or actual != expected:
        raise RuntimeError(
            f"full IRC endpoints {sorted(actual)} do not match reactant/product "
            f"basins {sorted(expected)}"
        )
    return {
        "expected": [list(value) for value in sorted(expected)],
        "actual": [list(value) for value in sorted(actual)],
    }


def record_store(
    path: Path,
    reaction: str,
    clusters: dict[str, Cluster],
    frequencies_by_role: dict[str, FrequencyResult],
    energies: dict[str, dict[str, float]],
    summary: dict[str, Any],
) -> None:
    temporary = path.with_name(f".{path.name}.{time.time_ns()}.tmp")
    if temporary.exists():
        temporary.unlink()
    with Store(temporary) as store:
        structure_ids: dict[str, int] = {}
        for role, cluster in clusters.items():
            structure_ids[role] = store.add_structure(
                f"{reaction}-{role}",
                cluster.formula,
                cluster.to_xyz(),
                charge=cluster.charge,
                spin=cluster.spin,
            )
        for role, freq in frequencies_by_role.items():
            job = store.add_job(
                structure_ids[role],
                "freq",
                R2SCAN3C_METHOD,
                "gpu4pyscf" if summary["gpu"] else "pyscf",
                detail=json.dumps({"composite": "r2scan3c"}, sort_keys=True),
            )
            store.set_job_status(job, "done")
            store.add_result(job, "electronic", freq.electronic_hartree, "hartree")
            store.add_result(job, "imaginary_count", freq.n_imaginary, "count")
        for method, values in energies.items():
            for role, value in values.items():
                job = store.add_job(
                    structure_ids[role],
                    "sp",
                    method,
                    "gpu4pyscf" if summary["gpu"] else "pyscf",
                )
                store.set_job_status(job, "done")
                store.add_result(job, "electronic", value, "hartree")
        analysis = store.add_job(
            structure_ids["transition-state"],
            "analysis",
            PRODUCTION_METHOD,
            "quarry",
            detail=json.dumps({"temperature_k": TEMPERATURE_K}, sort_keys=True),
        )
        store.set_job_status(analysis, "done")
        for key in ("production_dg_kj", "production_dh_kj", "production_ds_kj_per_k"):
            units = "kJ/mol/K" if key.endswith("per_k") else "kJ/mol"
            store.add_result(analysis, key, float(summary[key]), units)
    temporary.replace(path)


def run(args: argparse.Namespace) -> int:
    run_dir = args.run_dir.resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    source_store = args.source_store.resolve()
    reactant_source, reactant_receipt = load_store_structure(
        source_store, args.reactant_id
    )
    ts_source, ts_receipt = load_store_structure(source_store, args.ts_id)
    if reactant_source.symbols != ts_source.symbols:
        raise ValueError("reactant/TS atom order differs")
    product_source = load_xyz_like(
        args.product_reference.resolve(),
        reactant_source,
        name=f"{args.reaction}-product",
    )
    source = {
        "store": str(source_store),
        "store_sha256": sha256_path(source_store),
        "reactant": reactant_receipt,
        "transition_state": ts_receipt,
        "product_reference": str(args.product_reference.resolve()),
        "product_reference_sha256": sha256_path(args.product_reference.resolve()),
    }
    atomic_json(run_dir / "source-receipt.json", source)

    r2scan3c, production, b3lyp_d4 = settings(use_gpu=args.gpu)
    atomic_json(
        run_dir / "settings.json",
        {
            "r2scan3c": asdict(r2scan3c),
            "production": asdict(production),
            "b3lyp_d4": asdict(b3lyp_d4),
            "temperature_k": TEMPERATURE_K,
            "bounds": {
                "minimum_steps": args.minimum_steps,
                "saddle_steps": args.saddle_steps,
                "irc_steps": args.irc_steps,
            },
        },
    )

    reactant = checkpoint_cluster(
        run_dir / "reactant.r2scan3c.xyz",
        reactant_source,
        lambda: optimize(reactant_source, r2scan3c, max_steps=args.minimum_steps),
    )
    product = checkpoint_cluster(
        run_dir / "product.r2scan3c.xyz",
        product_source,
        lambda: optimize(product_source, r2scan3c, max_steps=args.minimum_steps),
    )
    transition_state = checkpoint_cluster(
        run_dir / "transition-state.r2scan3c.xyz",
        ts_source,
        lambda: find_ts(
            ts_source,
            r2scan3c,
            max_steps=args.saddle_steps,
            trajectory=str(run_dir / "transition-state.r2scan3c.traj"),
        ),
    )

    reactant_freq = checkpoint_frequency(
        run_dir / "reactant.r2scan3c.frequency.json", reactant, r2scan3c
    )
    product_freq = checkpoint_frequency(
        run_dir / "product.r2scan3c.frequency.json", product, r2scan3c
    )
    ts_freq = checkpoint_frequency(
        run_dir / "transition-state.r2scan3c.frequency.json",
        transition_state,
        r2scan3c,
    )
    if reactant_freq.n_imaginary != 0 or product_freq.n_imaginary != 0:
        raise RuntimeError("r2SCAN-3c reactant/product is not a true minimum")
    significant = ts_freq.imaginary_cm[ts_freq.imaginary_cm > args.imaginary_floor]
    if significant.size != 1:
        raise RuntimeError(
            "r2SCAN-3c TS does not have exactly one significant imaginary mode: "
            f"{ts_freq.imaginary_cm.tolist()}"
        )

    irc_forward_path = run_dir / "irc-forward.r2scan3c.xyz"
    irc_reverse_path = run_dir / "irc-reverse.r2scan3c.xyz"
    if irc_forward_path.exists() and irc_reverse_path.exists():
        irc_endpoints = (
            load_xyz_like(irc_forward_path, transition_state, name="irc-forward"),
            load_xyz_like(irc_reverse_path, transition_state, name="irc-reverse"),
        )
    else:
        irc_endpoints = full_irc(
            transition_state,
            r2scan3c,
            max_steps=args.irc_steps,
            trajectory=run_dir / "full-irc.r2scan3c.traj",
            logfile=run_dir / "full-irc.r2scan3c.log",
        )
        irc_forward_path.write_text(irc_endpoints[0].to_xyz())
        irc_reverse_path.write_text(irc_endpoints[1].to_xyz())
    irc_receipt = require_si_neutral_irc(
        irc_endpoints, reactant, product, args.attacker_index
    )

    clusters = {
        "reactant": reactant,
        "product": product,
        "transition-state": transition_state,
    }
    method_settings = {
        PRODUCTION_METHOD: production,
        B3LYP_D4_METHOD: b3lyp_d4,
    }
    energies: dict[str, dict[str, float]] = {}
    for method, method_settings_value in method_settings.items():
        method_slug = method.split("/")[0].replace("(", "-").replace(")", "")
        energies[method] = {}
        for role, cluster in clusters.items():
            energies[method][role] = checkpoint_energy(
                run_dir / f"{role}.{method_slug}.energy.json",
                cluster,
                method_settings_value,
                method,
            )

    production_rate = rate_from_thermo(
        thermo(reactant_freq, energies[PRODUCTION_METHOD]["reactant"]),
        thermo(ts_freq, energies[PRODUCTION_METHOD]["transition-state"]),
        imag_nu_cm=float(significant[0]),
        tunneling="wigner",
    )
    barriers: dict[str, float] = {
        R2SCAN3C_METHOD: (ts_freq.electronic_hartree - reactant_freq.electronic_hartree)
        * HARTREE_TO_KJ,
        PRODUCTION_METHOD: (
            energies[PRODUCTION_METHOD]["transition-state"]
            - energies[PRODUCTION_METHOD]["reactant"]
        )
        * HARTREE_TO_KJ,
        B3LYP_D4_METHOD: (
            energies[B3LYP_D4_METHOD]["transition-state"]
            - energies[B3LYP_D4_METHOD]["reactant"]
        )
        * HARTREE_TO_KJ,
    }
    summary = {
        "reaction": args.reaction,
        "gpu": args.gpu,
        "temperature_k": TEMPERATURE_K,
        "source": source,
        "geometry_method": R2SCAN3C_METHOD,
        "production_method": PRODUCTION_METHOD,
        "barrier_electronic_kj": barriers,
        "production_dg_kj": production_rate.dg_kj,
        "production_dh_kj": production_rate.dh_kj,
        "production_ds_kj_per_k": production_rate.ds_kj_per_k,
        "wigner_kappa": production_rate.kappa,
        "k_per_s": production_rate.k,
        "ts_imaginary_cm": float(significant[0]),
        "full_irc": irc_receipt,
        "geometry_hashes": {
            role: geometry_hash(cluster.to_xyz()) for role, cluster in clusters.items()
        },
    }
    atomic_json(run_dir / "results.json", summary)
    record_store(
        run_dir / "store.sqlite",
        args.reaction,
        clusters,
        {
            "reactant": reactant_freq,
            "product": product_freq,
            "transition-state": ts_freq,
        },
        energies,
        summary,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--reaction", required=True)
    result.add_argument("--source-store", type=Path, required=True)
    result.add_argument("--reactant-id", type=int, required=True)
    result.add_argument("--ts-id", type=int, required=True)
    result.add_argument("--product-reference", type=Path, required=True)
    result.add_argument("--attacker-index", type=int, required=True)
    result.add_argument("--run-dir", type=Path, required=True)
    result.add_argument("--gpu", action="store_true")
    result.add_argument("--minimum-steps", type=int, default=200)
    result.add_argument("--saddle-steps", type=int, default=400)
    result.add_argument("--irc-steps", type=int, default=400)
    result.add_argument("--imaginary-floor", type=float, default=30.0)
    # Parsed first by ``bootstrap_cli`` so limits apply before NumPy/PySCF import;
    # retained here so the full parser accepts the same invocation.
    result.add_argument("--threads", type=int, default=16)
    result.add_argument("--nice", type=int, default=10)
    result.add_argument("--log")
    return result


def main() -> int:
    args = parser().parse_args()
    for name in ("minimum_steps", "saddle_steps", "irc_steps"):
        if getattr(args, name) <= 0:
            raise ValueError(f"{name.replace('_', '-')} must be positive")
    if args.attacker_index < 0:
        raise ValueError("attacker-index must be non-negative")
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
