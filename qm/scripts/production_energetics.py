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
    frequencies_finite_difference,
    frequency_geometry_fingerprint,
    frequency_settings_fingerprint,
)
from quarry.rates import rate_from_thermo, thermo_from_frequencies
from quarry.store import Store, geometry_hash
from quarry.ts import (
    find_ts,
    full_irc,
    make_ase_calculator,
    reaction_path_vector,
)

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


def exact_xyz(cluster: Cluster) -> str:
    """Round-trip-safe XYZ for scientific checkpoints and cache identities."""
    lines = [str(len(cluster.symbols)), cluster.name]
    lines.extend(
        f"{symbol} {x:.17g} {y:.17g} {z:.17g}"
        for symbol, (x, y, z) in zip(
            cluster.symbols,
            cluster.coords,
            strict=True,
        )
    )
    return "\n".join(lines)


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
    *,
    identity: dict[str, Any],
) -> Cluster:
    receipt_path = path.with_name(f"{path.name}.checkpoint.json")
    expected = {
        "source_geometry_fingerprint": frequency_geometry_fingerprint(template),
        "identity": identity,
    }
    if path.exists() or receipt_path.exists():
        if not path.exists() or not receipt_path.exists():
            raise ValueError(f"{path.name}: incomplete checkpoint pair")
        receipt = json.loads(receipt_path.read_text())
        if (
            receipt.get("source_geometry_fingerprint")
            != expected["source_geometry_fingerprint"]
        ):
            raise ValueError(f"{path.name}: checkpoint source geometry drift")
        if receipt.get("identity") != identity:
            raise ValueError(f"{path.name}: checkpoint stage identity drift")
        loaded = load_xyz_like(path, template, name=path.stem)
        if receipt.get("output_geometry_fingerprint") != frequency_geometry_fingerprint(
            loaded
        ):
            raise ValueError(f"{path.name}: checkpoint output geometry drift")
        return loaded
    result = compute()
    if result.symbols != template.symbols:
        raise ValueError(f"{path.name}: computed atom order drift")
    temporary = path.with_name(f".{path.name}.{time.time_ns()}.tmp")
    temporary.write_text(exact_xyz(result))
    temporary.replace(path)
    loaded = load_xyz_like(path, template, name=path.stem)
    atomic_json(
        receipt_path,
        {
            **expected,
            "output_geometry_fingerprint": frequency_geometry_fingerprint(loaded),
            "output_sha256": sha256_path(path),
        },
    )
    return loaded


def checkpoint_irc_endpoints(
    run_dir: Path,
    transition_state: Cluster,
    compute: Callable[[], tuple[Cluster, Cluster]],
    *,
    identity: dict[str, Any],
) -> tuple[Cluster, Cluster]:
    """Reuse IRC endpoints only when bound to this TS, settings, and code."""
    forward_path = run_dir / "irc-forward.r2scan3c.xyz"
    reverse_path = run_dir / "irc-reverse.r2scan3c.xyz"
    receipt_path = run_dir / "irc.r2scan3c.receipt.json"
    expected = {
        "transition_state_geometry_fingerprint": frequency_geometry_fingerprint(
            transition_state
        ),
        "identity": identity,
    }
    paths = (forward_path, reverse_path, receipt_path)
    if any(path.exists() for path in paths):
        if not all(path.exists() for path in paths):
            raise ValueError("incomplete IRC checkpoint set")
        receipt = json.loads(receipt_path.read_text())
        if (
            receipt.get("transition_state_geometry_fingerprint")
            != expected["transition_state_geometry_fingerprint"]
        ):
            raise ValueError("IRC checkpoint transition-state geometry drift")
        if receipt.get("identity") != identity:
            raise ValueError("IRC checkpoint identity drift")
        loaded = (
            load_xyz_like(forward_path, transition_state, name="irc-forward"),
            load_xyz_like(reverse_path, transition_state, name="irc-reverse"),
        )
        directions = receipt.get("directions")
        if not isinstance(directions, dict):
            raise ValueError("IRC checkpoint missing direction receipts")
        for cluster, direction, path in (
            (loaded[0], "forward", forward_path),
            (loaded[1], "reverse", reverse_path),
        ):
            payload = directions.get(direction)
            if not isinstance(payload, dict):
                raise ValueError(f"IRC checkpoint missing {direction} receipt")
            if payload.get("path") != path.name:
                raise ValueError(f"IRC checkpoint {direction} path drift")
            if payload.get("geometry_fingerprint") != frequency_geometry_fingerprint(
                cluster
            ):
                raise ValueError(f"IRC checkpoint {direction} geometry drift")
            if payload.get("sha256") != sha256_path(path):
                raise ValueError(f"IRC checkpoint {direction} hash drift")
        return loaded
    endpoints = compute()
    if len(endpoints) != 2:
        raise ValueError("full IRC must return forward and reverse endpoints")
    for cluster, path in (
        (endpoints[0], forward_path),
        (endpoints[1], reverse_path),
    ):
        if cluster.symbols != transition_state.symbols:
            raise ValueError(f"{path.name}: computed atom order drift")
    staged: list[Path] = []
    published: list[Path] = []
    try:
        for cluster, path in (
            (endpoints[0], forward_path),
            (endpoints[1], reverse_path),
        ):
            temporary = path.with_name(f".{path.name}.{time.time_ns()}.tmp")
            temporary.write_text(exact_xyz(cluster))
            staged.append(temporary)
        for temporary, path in zip(
            staged,
            (forward_path, reverse_path),
            strict=True,
        ):
            temporary.replace(path)
            published.append(path)
        loaded = (
            load_xyz_like(forward_path, transition_state, name="irc-forward"),
            load_xyz_like(reverse_path, transition_state, name="irc-reverse"),
        )
        atomic_json(
            receipt_path,
            {
                **expected,
                "directions": {
                    "forward": {
                        "path": forward_path.name,
                        "geometry_fingerprint": frequency_geometry_fingerprint(
                            loaded[0]
                        ),
                        "sha256": sha256_path(forward_path),
                    },
                    "reverse": {
                        "path": reverse_path.name,
                        "geometry_fingerprint": frequency_geometry_fingerprint(
                            loaded[1]
                        ),
                        "sha256": sha256_path(reverse_path),
                    },
                },
            },
        )
    except Exception:
        for temporary in staged:
            temporary.unlink(missing_ok=True)
        if not receipt_path.exists():
            for path in published:
                path.unlink(missing_ok=True)
        raise
    return loaded


def optimize_minimum(
    cluster: Cluster,
    settings: DftSettings,
    *,
    max_steps: int,
    trajectory: Path,
    fmax_ev_a: float = 0.02,
) -> Cluster:
    """Converge a minimum on the exact composite surface with ASE BFGS."""
    from ase import Atoms
    from ase.constraints import FixAtoms
    from ase.optimize import BFGS

    atoms = Atoms(symbols=cluster.symbols, positions=cluster.coords)
    atoms.calc = make_ase_calculator(settings, cluster.charge, cluster.spin)
    if cluster.frozen_indices:
        atoms.set_constraint(FixAtoms(indices=cluster.frozen_indices))
    optimizer = BFGS(
        atoms,
        maxstep=0.05,
        trajectory=str(trajectory),
        logfile="-",
    )
    converged = optimizer.run(fmax=fmax_ev_a, steps=max_steps)
    projected_fmax = float(np.linalg.norm(atoms.get_forces(), axis=1).max())
    if not converged or projected_fmax >= fmax_ev_a:
        raise RuntimeError(
            f"ASE minimum did not converge within {max_steps} steps; "
            f"projected fmax={projected_fmax:.6f} eV/A"
        )
    return replace(
        cluster, coords=atoms.positions.copy(), name=f"{cluster.name}-minimum"
    )


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
    *,
    finite_difference: bool = False,
) -> FrequencyResult:
    expected_geometry = frequency_geometry_fingerprint(cluster)
    expected_settings = frequency_settings_fingerprint(settings)
    expected_hessian = "finite-difference-gradient" if finite_difference else "analytic"
    if path.exists():
        payload = json.loads(path.read_text())
        if payload.get("geometry_fingerprint") != expected_geometry:
            raise ValueError(f"{path.name}: cached geometry fingerprint drift")
        if payload.get("settings_fingerprint") != expected_settings:
            raise ValueError(f"{path.name}: cached settings fingerprint drift")
        if payload.get("hessian_method") != expected_hessian:
            raise ValueError(f"{path.name}: cached Hessian method drift")
        return frequency_from_payload(payload)
    result = (
        frequencies_finite_difference(cluster, settings)
        if finite_difference
        else frequencies(cluster, settings)
    )
    payload = frequency_payload(result)
    payload["hessian_method"] = expected_hessian
    atomic_json(path, payload)
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
        if payload.get("method") != method:
            raise ValueError(f"{path.name}: cached method drift")
        value = float(payload["electronic_hartree"])
        if not np.isfinite(value):
            raise ValueError(f"{path.name}: cached electronic energy is non-finite")
        return value
    value = float(energy(cluster, settings))
    if not np.isfinite(value):
        raise RuntimeError(f"{path.name}: computed electronic energy is non-finite")
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


def endpoint_identity(
    cluster: Cluster,
    attacker_index: int,
) -> tuple[
    tuple[bool, bool, bool],
    tuple[tuple[int, int], ...],
    tuple[tuple[int, int], ...],
]:
    """Typed neutral endpoint identity: basin, physical-H owners, heavy graph."""
    from ase.data import atomic_numbers, covalent_radii

    deltas = cluster.coords[:, None, :] - cluster.coords[None, :, :]
    distances = np.linalg.norm(deltas, axis=2)
    nonzero = distances[np.triu_indices(len(cluster.symbols), k=1)]
    if not np.all(np.isfinite(nonzero)) or float(np.min(nonzero)) < 0.55:
        raise ValueError(
            "endpoint has non-finite coordinates or a sub-0.55 A collision"
        )

    oxygen_indices = [
        index for index, symbol in enumerate(cluster.symbols) if symbol == "O"
    ]
    hydrogen_owners: list[tuple[int, int]] = []
    for hydrogen in (
        index for index, symbol in enumerate(cluster.symbols) if symbol == "H"
    ):
        candidates = sorted(
            (float(distances[hydrogen, oxygen]), oxygen) for oxygen in oxygen_indices
        )
        if candidates[0][0] >= 1.25:
            raise ValueError(f"hydrogen {hydrogen} is unassigned")
        if len(candidates) > 1 and candidates[1][0] < 1.05:
            raise ValueError(f"hydrogen {hydrogen} has ambiguous oxygen ownership")
        hydrogen_owners.append((hydrogen, candidates[0][1]))

    heavy = [index for index, symbol in enumerate(cluster.symbols) if symbol != "H"]
    heavy_bonds = []
    for position, left in enumerate(heavy):
        for right in heavy[position + 1 :]:
            cutoff = 1.25 * (
                covalent_radii[atomic_numbers[cluster.symbols[left]]]
                + covalent_radii[atomic_numbers[cluster.symbols[right]]]
            )
            if distances[left, right] < cutoff:
                heavy_bonds.append((left, right))
    return (
        si_neutral_signature(cluster, attacker_index),
        tuple(hydrogen_owners),
        tuple(heavy_bonds),
    )


def require_si_neutral_irc(
    endpoints: tuple[Cluster, Cluster],
    reactant: Cluster,
    product: Cluster,
    attacker_index: int,
) -> dict[str, Any]:
    for endpoint in (*endpoints, reactant, product):
        if endpoint.symbols != reactant.symbols:
            raise ValueError("IRC endpoint atom order differs from the reactant")
        if endpoint.charge != reactant.charge or endpoint.spin != reactant.spin:
            raise ValueError("IRC endpoint electronic state differs from the reactant")
    expected = {
        endpoint_identity(reactant, attacker_index),
        endpoint_identity(product, attacker_index),
    }
    actual = {endpoint_identity(endpoint, attacker_index) for endpoint in endpoints}
    if len(expected) != 2 or actual != expected:
        expected_basins = sorted(identity[0] for identity in expected)
        actual_basins = sorted(identity[0] for identity in actual)
        raise RuntimeError(
            f"full IRC endpoints {actual_basins} do not match reactant/product "
            f"basins {expected_basins} with exact H ownership/heavy topology"
        )
    return {
        "expected": [
            list(value) for value in sorted(identity[0] for identity in expected)
        ],
        "actual": [list(value) for value in sorted(identity[0] for identity in actual)],
        "typed_identity": "basin+physical-hydrogen-owners+heavy-atom-bonds",
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
                exact_xyz(cluster),
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

    r2scan_identity = frequency_settings_fingerprint(r2scan3c)
    reactant = checkpoint_cluster(
        run_dir / "reactant.r2scan3c.xyz",
        reactant_source,
        lambda: optimize_minimum(
            reactant_source,
            r2scan3c,
            max_steps=args.minimum_steps,
            trajectory=run_dir / "reactant.r2scan3c.traj",
        ),
        identity={
            "stage": "minimum",
            "role": "reactant",
            "algorithm": "ase-bfgs-v1",
            "settings": r2scan_identity,
            "max_steps": args.minimum_steps,
        },
    )
    product = checkpoint_cluster(
        run_dir / "product.r2scan3c.xyz",
        product_source,
        lambda: optimize_minimum(
            product_source,
            r2scan3c,
            max_steps=args.minimum_steps,
            trajectory=run_dir / "product.r2scan3c.traj",
        ),
        identity={
            "stage": "minimum",
            "role": "product",
            "algorithm": "ase-bfgs-v1",
            "settings": r2scan_identity,
            "max_steps": args.minimum_steps,
        },
    )
    transition_state = checkpoint_cluster(
        run_dir / "transition-state.r2scan3c.xyz",
        ts_source,
        lambda: find_ts(
            ts_source,
            r2scan3c,
            max_steps=args.saddle_steps,
            trajectory=str(run_dir / "transition-state.r2scan3c.traj"),
            initial_mode=reaction_path_vector(reactant, product),
            internal=False,
        ),
        identity={
            "stage": "transition-state",
            "algorithm": "sella-directed-cartesian-v1",
            "settings": r2scan_identity,
            "max_steps": args.saddle_steps,
            "reactant_geometry": frequency_geometry_fingerprint(reactant),
            "product_geometry": frequency_geometry_fingerprint(product),
        },
    )

    reactant_freq = checkpoint_frequency(
        run_dir / "reactant.r2scan3c.frequency.json",
        reactant,
        r2scan3c,
        finite_difference=True,
    )
    product_freq = checkpoint_frequency(
        run_dir / "product.r2scan3c.frequency.json",
        product,
        r2scan3c,
        finite_difference=True,
    )
    ts_freq = checkpoint_frequency(
        run_dir / "transition-state.r2scan3c.frequency.json",
        transition_state,
        r2scan3c,
        finite_difference=True,
    )
    if reactant_freq.n_imaginary != 0 or product_freq.n_imaginary != 0:
        raise RuntimeError("r2SCAN-3c reactant/product is not a true minimum")
    significant = ts_freq.imaginary_cm[ts_freq.imaginary_cm > args.imaginary_floor]
    if significant.size != 1:
        raise RuntimeError(
            "r2SCAN-3c TS does not have exactly one significant imaginary mode: "
            f"{ts_freq.imaginary_cm.tolist()}"
        )

    irc_endpoints = checkpoint_irc_endpoints(
        run_dir,
        transition_state,
        lambda: full_irc(
            transition_state,
            r2scan3c,
            max_steps=args.irc_steps,
            trajectory=run_dir / "full-irc.r2scan3c.traj",
            logfile=run_dir / "full-irc.r2scan3c.log",
        ),
        identity={
            "stage": "full-irc",
            "algorithm": "sella-gonzalez-schlegel-full-irc-v1",
            "settings": r2scan_identity,
            "max_steps": args.irc_steps,
        },
    )
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
            role: geometry_hash(exact_xyz(cluster))
            for role, cluster in clusters.items()
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
    result.add_argument("--reaction", choices=("si-neutral",), required=True)
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


def execute_with_status(args: argparse.Namespace) -> int:
    status_path = args.run_dir.resolve() / "run_status.json"
    atomic_json(
        status_path,
        {
            "status": "running",
            "reaction": args.reaction,
            "started_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        },
    )
    try:
        result = run(args)
    except Exception as exc:
        atomic_json(
            status_path,
            {
                "status": "failed",
                "reaction": args.reaction,
                "error_type": type(exc).__name__,
                "error": str(exc),
                "failed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            },
        )
        raise
    atomic_json(
        status_path,
        {
            "status": "completed",
            "reaction": args.reaction,
            "completed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "results_sha256": sha256_path(args.run_dir.resolve() / "results.json"),
            "store_sha256": sha256_path(args.run_dir.resolve() / "store.sqlite"),
        },
    )
    return result


def main() -> int:
    args = parser().parse_args()
    for name in ("minimum_steps", "saddle_steps", "irc_steps"):
        if getattr(args, name) <= 0:
            raise ValueError(f"{name.replace('_', '-')} must be positive")
    if args.attacker_index < 0:
        raise ValueError("attacker-index must be non-negative")
    return execute_with_status(args)


if __name__ == "__main__":
    raise SystemExit(main())
