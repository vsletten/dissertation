#!/usr/bin/env python3
"""Dedicated CALC-005 state probe, bounded production pilot, and validator."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import time
from contextlib import contextmanager
from dataclasses import asdict, replace
from pathlib import Path
from typing import Any, Protocol

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

DEFAULT_OUTPUT_ROOT = Path(
    "/mnt/data/vsletten/dissertation-data/a3h-calc005-si-n1-pilot"
)
MAX_STEPS = 150
GRADIENT_RMS_MAX = 3.0e-4
GRADIENT_MAX_MAX = 4.5e-4
IMAGINARY_NOISE_FLOOR_CM = 30.0
MIN_PAIR_DISTANCE_A = 0.75
MAX_OWNER_DISTANCE_A = 1.25
MIN_OWNER_MARGIN_A = 0.15
R2SCAN3C_METHOD = "r2scan-3c/def2-mtzvpp/d4/gcp"
PRODUCTION_METHOD = "wb97m-v/def2-tzvpd/smd(water)"
CLASSIFICATION = "thermodynamic/non-kinetic/non-emittable"
GPU_TTL_HOURS = 13.0
EXPECTED_MEMORY_MAX_BYTES = 32 * 1024**3
EXPECTED_MEMORY_SWAP_MAX_BYTES = 4 * 1024**3
EXPECTED_CPU_QUOTA_PERCENT = 1600
EXPECTED_NICE = 10
EXPECTED_THREADS = 16
EXPECTED_RUNTIME_MAX_SECONDS = 43200
EXPECTED_GPU_OWNER = "calc005_si_attachment"

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    bootstrap_cli(
        "calc005_si_attachment",
        default_run_root=DEFAULT_OUTPUT_ROOT / "logs",
        gpu_owner="calc005_si_attachment",
        gpu_ttl_hours=GPU_TTL_HOURS,
    )

import numpy as np  # noqa: E402

from quarry.calc005 import (  # noqa: E402
    STOICHIOMETRY,
    TEMPERATURE_K,
    Calc005Pair,
    atom_map_sha256,
    build_calc005_pair,
    compose_thermochemistry,
    probe_calc005_pairs,
)
from quarry.calc005_store import (  # noqa: E402
    calculation_receipt_sha256,
    receipt_payload_sha256,
    validate_calc005_store,
    validate_staged_calc005_store,
    write_calc005_store,
)
from quarry.clusters import Cluster  # noqa: E402
from quarry.pipeline import (  # noqa: E402
    HARTREE_TO_KJ,
    DftSettings,
    FrequencyResult,
    OptimizationResult,
    frequency_geometry_fingerprint,
    frequency_settings_fingerprint,
)
from quarry.rates import (  # noqa: E402
    surface_thermo_from_frequencies,
    thermo_from_frequencies,
)
from quarry.store import geometry_hash  # noqa: E402


class CalculatorBackend(Protocol):
    """Mockable electronic-structure boundary; production is the default CLI backend."""

    def optimize(
        self,
        role: str,
        cluster: Cluster,
        settings: DftSettings,
        *,
        max_steps: int,
    ) -> OptimizationResult: ...

    def gradient(
        self, role: str, cluster: Cluster, settings: DftSettings
    ) -> np.ndarray: ...

    def frequencies(
        self, role: str, cluster: Cluster, settings: DftSettings
    ) -> FrequencyResult: ...

    def energy(self, role: str, cluster: Cluster, settings: DftSettings) -> float: ...


class PipelineBackend:
    """Production PySCF/GPU4PySCF implementation of the calculator boundary."""

    name = "gpu4pyscf/pyscf"

    def optimize(self, role, cluster, settings, *, max_steps):
        from quarry.pipeline import optimize_bounded

        return optimize_bounded(cluster, settings, max_steps=max_steps)

    def gradient(self, role, cluster, settings):
        from quarry.pipeline import gradient

        return gradient(cluster, settings)

    def frequencies(self, role, cluster, settings):
        from quarry.pipeline import frequencies

        return frequencies(cluster, settings)

    def energy(self, role, cluster, settings):
        from quarry.pipeline import energy

        return energy(cluster, settings)


class AnalyticalBackend:
    """Deterministic no-QM fixture backend; never selected by the CLI."""

    name = "analytical-test-backend"

    def optimize(self, role, cluster, settings, *, max_steps):
        return OptimizationResult(cluster, True, max_steps)

    def gradient(self, role, cluster, settings):
        return np.zeros_like(cluster.coords)

    def frequencies(self, role, cluster, settings):
        mode_count = expected_mode_count(role, cluster)
        return FrequencyResult(
            frequencies_cm=np.linspace(100.0, 1700.0, mode_count),
            imaginary_cm=np.array([], dtype=float),
            electronic_hartree=-100.0,
            molar_mass_kg=0.1,
            rotational_temperatures_k=(1.0, 2.0, 3.0) if role == "SiOH4" else None,
            linear=False,
            geometry_fingerprint=frequency_geometry_fingerprint(cluster),
            settings_fingerprint=frequency_settings_fingerprint(settings),
        )

    def energy(self, role, cluster, settings):
        return {"C": -300.0, "V": -200.0, "SiOH4": -100.02}[role]


class ComputationalFailure(RuntimeError):
    """A bounded calculator call failed or exhausted its declared budget."""


class PhysicalStateFailure(RuntimeError):
    """An endpoint failed an independent physical or numerical acceptance gate."""


def now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        Path(temporary_name).replace(path)
    finally:
        Path(temporary_name).unlink(missing_ok=True)


def exact_xyz(cluster: Cluster) -> str:
    lines = [str(len(cluster.symbols)), cluster.name]
    lines.extend(
        f"{symbol} {x:.17g} {y:.17g} {z:.17g}"
        for symbol, (x, y, z) in zip(cluster.symbols, cluster.coords, strict=True)
    )
    return "\n".join(lines) + "\n"


def atomic_xyz(path: Path, cluster: Cluster) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "w") as stream:
            stream.write(exact_xyz(cluster))
            stream.flush()
            os.fsync(stream.fileno())
        Path(temporary_name).replace(path)
    finally:
        Path(temporary_name).unlink(missing_ok=True)


def read_xyz(path: Path, template: Cluster) -> Cluster:
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        raise RuntimeError(f"malformed XYZ: {path}")
    try:
        count = int(lines[0])
        rows = [line.split() for line in lines[2:]]
        symbols = [row[0] for row in rows]
        coords = np.asarray([[float(value) for value in row[1:]] for row in rows])
    except (IndexError, ValueError) as exc:
        raise RuntimeError(f"malformed XYZ: {path}") from exc
    if (
        count != len(template.symbols)
        or len(rows) != count
        or any(len(row) != 4 for row in rows)
        or symbols != template.symbols
        or coords.shape != template.coords.shape
        or not np.all(np.isfinite(coords))
    ):
        raise RuntimeError(f"XYZ identity/order/coordinate mismatch: {path}")
    return replace(template, name=lines[1], coords=coords)


def calc005_settings(*, use_gpu: bool) -> tuple[DftSettings, DftSettings]:
    return (
        DftSettings(
            xc="r2scan",
            basis="def2-mtzvpp",
            composite="r2scan3c",
            density_fit=True,
            use_gpu=use_gpu,
        ),
        DftSettings(
            xc="wb97m-v",
            basis="def2-tzvpd",
            solvent="smd",
            density_fit=True,
            use_gpu=use_gpu,
        ),
    )


def _settings_hash(settings: DftSettings) -> str:
    return hashlib.sha256(
        json.dumps(asdict(settings), separators=(",", ":"), sort_keys=True).encode()
    ).hexdigest()


def _artifact(path: Path, cluster: Cluster) -> dict[str, Any]:
    return {
        "path": str(path),
        "sha256": sha256_path(path),
        "geometry_hash": geometry_hash(exact_xyz(cluster)),
        "geometry_fingerprint": frequency_geometry_fingerprint(cluster),
        "formula": cluster.formula,
        "charge": cluster.charge,
        "spin": cluster.spin,
        "frozen_indices": list(cluster.frozen_indices),
    }


def _minimum_pair(coords: np.ndarray) -> float:
    delta = coords[:, None, :] - coords[None, :, :]
    distances = np.linalg.norm(delta, axis=2)
    distances[np.diag_indices_from(distances)] = np.inf
    return float(np.min(distances))


def _owner_map(cluster: Cluster) -> tuple[dict[int, int], float, float]:
    oxygens = [i for i, symbol in enumerate(cluster.symbols) if symbol == "O"]
    owners: dict[int, int] = {}
    max_distance = 0.0
    min_margin = math.inf
    for h, symbol in enumerate(cluster.symbols):
        if symbol != "H":
            continue
        distances = sorted(
            (float(np.linalg.norm(cluster.coords[h] - cluster.coords[o])), o)
            for o in oxygens
        )
        if len(distances) < 2 or not all(
            math.isfinite(value) for value, _ in distances
        ):
            raise PhysicalStateFailure(f"H{h} lacks finite oxygen-owner distances")
        max_distance = max(max_distance, distances[0][0])
        min_margin = min(min_margin, distances[1][0] - distances[0][0])
        if distances[0][0] > MAX_OWNER_DISTANCE_A:
            raise PhysicalStateFailure(
                f"H{h} owner distance exceeds {MAX_OWNER_DISTANCE_A} A"
            )
        if distances[1][0] - distances[0][0] < MIN_OWNER_MARGIN_A:
            raise PhysicalStateFailure(f"H{h} oxygen ownership is ambiguous")
        owners[h] = distances[0][1]
    return owners, max_distance, min_margin


def _origin_token(
    origin: Any,
) -> tuple[str, tuple[int, tuple[int, int, int]], int | None]:
    return origin.kind, origin.node, origin.ordinal


def _heavy_topology(
    cluster: Cluster, origins: tuple[Any, ...]
) -> set[tuple[tuple[str, tuple[int, tuple[int, int, int]], int | None], ...]]:
    """Infer cation-oxygen connectivity and bind every edge to atom provenance."""

    if len(origins) != len(cluster.symbols):
        raise PhysicalStateFailure("heavy topology lacks complete atom provenance")
    edges = set()
    for cation, symbol in enumerate(cluster.symbols):
        if symbol not in {"Al", "Si"}:
            continue
        for oxygen, oxygen_symbol in enumerate(cluster.symbols):
            if oxygen_symbol != "O":
                continue
            distance = float(
                np.linalg.norm(cluster.coords[cation] - cluster.coords[oxygen])
            )
            if distance <= 2.35:
                edges.add(
                    (_origin_token(origins[cation]), _origin_token(origins[oxygen]))
                )
    return edges


def structural_gate(
    role: str,
    endpoint: Cluster,
    reference: Cluster,
    origins: tuple[Any, ...],
) -> dict[str, Any]:
    if endpoint.symbols != reference.symbols:
        raise PhysicalStateFailure(f"{role} changed atom identity/order")
    if endpoint.charge != reference.charge or endpoint.spin != reference.spin:
        raise PhysicalStateFailure(f"{role} changed charge/spin")
    if endpoint.frozen_indices != reference.frozen_indices:
        raise PhysicalStateFailure(f"{role} changed frozen atom indices")
    if endpoint.coords.shape != reference.coords.shape or not np.all(
        np.isfinite(endpoint.coords)
    ):
        raise PhysicalStateFailure(f"{role} has malformed or non-finite coordinates")
    frozen = sorted(reference.frozen_indices)
    if frozen and not np.array_equal(endpoint.coords[frozen], reference.coords[frozen]):
        raise PhysicalStateFailure(f"{role} changed frozen coordinates")
    minimum_pair = _minimum_pair(endpoint.coords)
    if not math.isfinite(minimum_pair) or minimum_pair < MIN_PAIR_DISTANCE_A:
        raise PhysicalStateFailure(f"{role} has collision at {minimum_pair} A")
    reference_owners, _, _ = _owner_map(reference)
    observed_owners, max_owner, min_margin = _owner_map(endpoint)
    if observed_owners != reference_owners:
        raise PhysicalStateFailure(f"{role} changed proton owner")
    expected_topology = _heavy_topology(reference, origins)
    observed_topology = _heavy_topology(endpoint, origins)
    if observed_topology != expected_topology:
        raise PhysicalStateFailure(
            f"{role} changed provenance-bound heavy connectivity"
        )
    center_token = ("deck", (4, (0, 0, 0)), None)
    expected_center_oxygen = {
        edge[1] for edge in expected_topology if edge[0] == center_token
    }
    observed_center_oxygen = {
        edge[1] for edge in observed_topology if edge[0] == center_token
    }
    if role in {"C", "SiOH4"} and (
        len(expected_center_oxygen) != 4
        or observed_center_oxygen != expected_center_oxygen
    ):
        raise PhysicalStateFailure(
            f"{role} does not retain the exact four center Si-O bonds"
        )
    if role == "V" and (expected_center_oxygen or observed_center_oxygen):
        raise PhysicalStateFailure("V unexpectedly retains the removed center Si")
    topology_payload = sorted(
        (repr(cation), repr(oxygen)) for cation, oxygen in observed_topology
    )
    return {
        "minimum_pair_distance_a": minimum_pair,
        "max_owner_distance_a": max_owner,
        "min_owner_margin_a": min_margin,
        "owner_map": [f"H{h}:O{o}" for h, o in sorted(observed_owners.items())],
        "heavy_edge_count": len(observed_topology),
        "heavy_topology_sha256": hashlib.sha256(
            json.dumps(topology_payload, separators=(",", ":")).encode()
        ).hexdigest(),
        "center_si_o_bond_count": len(observed_center_oxygen),
        "frozen_coordinates_bitwise_unchanged": True,
    }


def gradient_gate(
    role: str,
    cluster: Cluster,
    gradient: np.ndarray,
) -> dict[str, Any]:
    projected = np.asarray(gradient, dtype=float).copy()
    if projected.shape != cluster.coords.shape or not np.all(np.isfinite(projected)):
        raise PhysicalStateFailure(
            f"{role} independent gradient is malformed or non-finite"
        )
    projected[cluster.frozen_indices] = 0.0
    free = sorted(set(range(len(cluster.symbols))) - set(cluster.frozen_indices))
    values = projected[free]
    rms = float(np.sqrt(np.mean(values**2))) if values.size else 0.0
    maximum = float(np.max(np.abs(values))) if values.size else 0.0
    if not math.isfinite(rms) or not math.isfinite(maximum):
        raise PhysicalStateFailure(f"{role} projected gradient transform is non-finite")
    if rms > GRADIENT_RMS_MAX or maximum > GRADIENT_MAX_MAX:
        raise PhysicalStateFailure(
            f"{role} projected gradient exceeds RMS/max gates: {rms}/{maximum}"
        )
    return {
        "rms_hartree_per_bohr": rms,
        "max_hartree_per_bohr": maximum,
        "rms_threshold_hartree_per_bohr": GRADIENT_RMS_MAX,
        "max_threshold_hartree_per_bohr": GRADIENT_MAX_MAX,
        "free_atom_indices": free,
        "passed": True,
    }


def expected_mode_count(role: str, cluster: Cluster) -> int:
    """Return the exact pipeline mode cardinality for the fixed Hessian contract."""

    if role == "SiOH4":
        return 3 * len(cluster.symbols) - 6
    if role in {"C", "V"}:
        return 3 * (len(cluster.symbols) - len(cluster.frozen_indices))
    raise PhysicalStateFailure(f"unknown CALC-005 role: {role}")


def _frequency_payload(
    role: str,
    result: FrequencyResult,
    cluster: Cluster,
    settings: DftSettings,
) -> dict[str, Any]:
    frequencies = np.asarray(result.frequencies_cm, dtype=float)
    imaginary = np.asarray(result.imaginary_cm, dtype=float)
    rotational = result.rotational_temperatures_k
    rotational_values = (
        np.asarray(rotational, dtype=float)
        if rotational is not None
        else np.asarray([], dtype=float)
    )
    expected_geometry = frequency_geometry_fingerprint(cluster)
    expected_settings = frequency_settings_fingerprint(settings)
    if (
        frequencies.ndim != 1
        or imaginary.ndim != 1
        or not np.all(np.isfinite(frequencies))
        or not np.all(np.isfinite(imaginary))
        or np.any(frequencies <= 0.0)
        or np.any(imaginary < 0.0)
        or not math.isfinite(float(result.electronic_hartree))
        or not math.isfinite(float(result.molar_mass_kg))
        or float(result.molar_mass_kg) <= 0.0
        or not np.all(np.isfinite(rotational_values))
        or np.any(rotational_values <= 0.0)
        or result.geometry_fingerprint != expected_geometry
        or result.settings_fingerprint != expected_settings
    ):
        raise PhysicalStateFailure(f"{role} frequency evidence is malformed or unbound")
    if np.any(imaginary > IMAGINARY_NOISE_FLOOR_CM):
        raise PhysicalStateFailure(
            f"{role} is not a minimum: {imaginary.tolist()} cm^-1"
        )
    hessian = "full" if role == "SiOH4" else "PHVA"
    if role != "SiOH4" and not cluster.frozen_indices:
        raise PhysicalStateFailure(f"{role} PHVA requires a frozen shell")
    if role == "SiOH4" and cluster.frozen_indices:
        raise PhysicalStateFailure("SiOH4 full Hessian cannot have frozen atoms")
    if role == "SiOH4" and (rotational is None or len(rotational) != 3):
        raise PhysicalStateFailure(
            "SiOH4 full thermochemistry requires three rotations"
        )
    if role != "SiOH4" and rotational is not None:
        raise PhysicalStateFailure(f"{role} PHVA cannot carry rigid-body rotations")
    expected_modes = expected_mode_count(role, cluster)
    observed_modes = len(frequencies) + len(imaginary)
    if observed_modes != expected_modes:
        raise PhysicalStateFailure(
            f"{role} frequency mode count is {observed_modes}; "
            f"expected {expected_modes}"
        )
    if role == "SiOH4" and result.linear:
        raise PhysicalStateFailure("SiOH4 is nonlinear and requires 3N-6 modes")
    return {
        "hessian": hessian,
        "frequencies_cm": frequencies.tolist(),
        "imaginary_cm": imaginary.tolist(),
        "real_mode_count": len(frequencies),
        "imaginary_mode_count": len(imaginary),
        "expected_mode_count": expected_modes,
        "observed_mode_count": observed_modes,
        "noise_floor_cm": IMAGINARY_NOISE_FLOOR_CM,
        "electronic_hartree": float(result.electronic_hartree),
        "molar_mass_kg": float(result.molar_mass_kg),
        "rotational_temperatures_k": (
            list(result.rotational_temperatures_k)
            if result.rotational_temperatures_k is not None
            else None
        ),
        "linear": bool(result.linear),
        "geometry_fingerprint": expected_geometry,
        "settings_fingerprint": expected_settings,
        "passed": True,
    }


def _thermochemistry(
    role: str,
    frequency: FrequencyResult,
    production_electronic_hartree: float,
) -> dict[str, float]:
    electronic_kj = production_electronic_hartree * HARTREE_TO_KJ
    if role == "SiOH4":
        thermo = thermo_from_frequencies(
            electronic_kj,
            frequency.frequencies_cm,
            TEMPERATURE_K,
            molar_mass_kg=frequency.molar_mass_kg,
            rotational_temperatures_k=(
                list(frequency.rotational_temperatures_k)
                if frequency.rotational_temperatures_k is not None
                else None
            ),
            linear=frequency.linear,
        )
    else:
        thermo = surface_thermo_from_frequencies(
            electronic_kj,
            frequency.frequencies_cm,
            TEMPERATURE_K,
        )
    values = {
        "electronic_kj_mol": electronic_kj,
        "zpe_kj_mol": thermo.zpe_kj,
        "thermal_kj_mol": thermo.thermal_kj,
        "entropy_kj_mol_k": thermo.entropy_kj_per_k,
    }
    if not all(math.isfinite(value) for value in values.values()):
        raise PhysicalStateFailure(f"{role} thermochemistry is non-finite")
    return values


def _component_signature(
    role: str,
    seed: Cluster,
    source: dict[str, Any],
    geometry_settings: DftSettings,
    production_settings: DftSettings,
) -> dict[str, Any]:
    return {
        "schema": "calc005-component-signature-v1",
        "role": role,
        "source": source,
        "seed_geometry_fingerprint": frequency_geometry_fingerprint(seed),
        "formula": seed.formula,
        "charge": seed.charge,
        "spin": seed.spin,
        "frozen_indices": list(seed.frozen_indices),
        "settings": {
            "geometry": asdict(geometry_settings),
            "geometry_sha256": _settings_hash(geometry_settings),
            "production": asdict(production_settings),
            "production_sha256": _settings_hash(production_settings),
        },
        "optimizer": {"max_steps": MAX_STEPS, "retry_allowed": False},
        "gates": {
            "gradient_rms_max_hartree_per_bohr": GRADIENT_RMS_MAX,
            "gradient_max_max_hartree_per_bohr": GRADIENT_MAX_MAX,
            "imaginary_noise_floor_cm": IMAGINARY_NOISE_FLOOR_CM,
            "minimum_pair_distance_a": MIN_PAIR_DISTANCE_A,
            "maximum_owner_distance_a": MAX_OWNER_DISTANCE_A,
            "minimum_owner_margin_a": MIN_OWNER_MARGIN_A,
        },
    }


def _checkpoint_sha256(record: dict[str, Any]) -> str:
    payload = {
        key: value
        for key, value in record.items()
        if key not in {"checkpoint_sha256", "accepted_at", "cluster", "exact_xyz"}
    }
    return hashlib.sha256(
        json.dumps(payload, separators=(",", ":"), sort_keys=True).encode()
    ).hexdigest()


def _load_accepted_component(
    directory: Path,
    role: str,
    seed: Cluster,
    origins: tuple[Any, ...],
    signature: dict[str, Any],
) -> dict[str, Any] | None:
    reservation_path = directory / "reservation.json"
    entered_path = directory / "optimizer-entered.json"
    receipt_path = directory / "receipt.json"
    existing = (
        [path for path in directory.glob("*") if path.is_file()]
        if directory.exists()
        else []
    )
    if not existing:
        return None
    if not reservation_path.is_file():
        raise ComputationalFailure(f"{role} has an incomplete checkpoint")
    try:
        reservation = json.loads(reservation_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ComputationalFailure(
            f"{role} has a corrupt reservation checkpoint: {exc}"
        ) from exc
    if (
        reservation.get("schema") != "calc005-component-reservation-v1"
        or reservation.get("status") != "optimizer-budget-reserved"
        or reservation.get("signature") != signature
    ):
        raise ComputationalFailure(f"{role} checkpoint identity drifted")
    if not receipt_path.is_file():
        if entered_path.exists():
            raise ComputationalFailure(f"{role} optimizer budget is already spent")
        return None
    if not entered_path.is_file():
        raise ComputationalFailure(f"{role} receipt lacks optimizer-entered evidence")
    try:
        entered = json.loads(entered_path.read_text())
        receipt = json.loads(receipt_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ComputationalFailure(f"{role} has corrupt spent evidence: {exc}") from exc
    if (
        entered.get("schema") != "calc005-optimizer-entered-v1"
        or entered.get("status") != "optimizer-entered"
        or entered.get("signature") != signature
        or receipt.get("signature") != signature
    ):
        raise ComputationalFailure(f"{role} spent checkpoint identity drifted")
    if receipt.get("status") != "accepted":
        raise ComputationalFailure(f"{role} optimizer budget is already spent")
    if receipt.get("checkpoint_sha256") != _checkpoint_sha256(receipt):
        raise ComputationalFailure(f"{role} accepted checkpoint payload drifted")
    raw_path = directory / "raw-endpoint.xyz"
    endpoint_path = directory / "endpoint.xyz"
    for key, path in (("raw_endpoint", raw_path), ("endpoint", endpoint_path)):
        if not path.is_file() or receipt.get(key, {}).get("sha256") != sha256_path(
            path
        ):
            raise ComputationalFailure(f"{role} accepted checkpoint artifact drifted")
    cluster = read_xyz(endpoint_path, seed)
    structure = structural_gate(role, cluster, seed, origins)
    if receipt.get("structure") != structure:
        raise ComputationalFailure(f"{role} accepted structural receipt drifted")
    optimizer = receipt.get("optimizer")
    if optimizer != {
        "converged": True,
        "observed_calls": 1,
        "observed_retries": 0,
        "observed_max_steps": MAX_STEPS,
    }:
        raise ComputationalFailure(f"{role} accepted optimizer receipt drifted")
    try:
        frequency = receipt["frequency"]
        production_energy = float(receipt["production_electronic_hartree"])
        frequency_result = FrequencyResult(
            frequencies_cm=np.asarray(frequency["frequencies_cm"], dtype=float),
            imaginary_cm=np.asarray(frequency["imaginary_cm"], dtype=float),
            electronic_hartree=float(frequency["electronic_hartree"]),
            molar_mass_kg=float(frequency["molar_mass_kg"]),
            rotational_temperatures_k=(
                tuple(float(value) for value in frequency["rotational_temperatures_k"])
                if frequency["rotational_temperatures_k"] is not None
                else None
            ),
            linear=bool(frequency["linear"]),
            geometry_fingerprint=frequency["geometry_fingerprint"],
            settings_fingerprint=frequency["settings_fingerprint"],
        )
        expected_frequency = _frequency_payload(
            role,
            frequency_result,
            cluster,
            DftSettings(**signature["settings"]["geometry"]),
        )
        if frequency != expected_frequency or receipt.get(
            "thermochemistry"
        ) != _thermochemistry(role, frequency_result, production_energy):
            raise ComputationalFailure(f"{role} accepted numerical receipt drifted")
    except (KeyError, TypeError, ValueError, PhysicalStateFailure) as exc:
        raise ComputationalFailure(
            f"{role} accepted numerical receipt is corrupt: {exc}"
        ) from exc
    record = dict(receipt)
    record["cluster"] = cluster
    record["exact_xyz"] = exact_xyz(cluster)
    return record


def _run_component(
    output_root: Path,
    role: str,
    seed: Cluster,
    origins: tuple[Any, ...],
    source: dict[str, Any],
    geometry_settings: DftSettings,
    production_settings: DftSettings,
    backend: CalculatorBackend,
) -> dict[str, Any]:
    directory = output_root / "components" / role
    signature = _component_signature(
        role, seed, source, geometry_settings, production_settings
    )
    resumed = _load_accepted_component(directory, role, seed, origins, signature)
    if resumed is not None:
        return resumed
    directory.mkdir(parents=True, exist_ok=True)
    reservation = {
        "schema": "calc005-component-reservation-v1",
        "status": "optimizer-budget-reserved",
        "reserved_at": now(),
        "signature": signature,
    }
    atomic_json(directory / "reservation.json", reservation)
    atomic_json(
        directory / "optimizer-entered.json",
        {
            "schema": "calc005-optimizer-entered-v1",
            "status": "optimizer-entered",
            "entered_at": now(),
            "signature": signature,
        },
    )
    try:
        optimized = backend.optimize(role, seed, geometry_settings, max_steps=MAX_STEPS)
    except Exception as exc:
        receipt = {
            "schema": "calc005-component-receipt-v1",
            "status": "optimizer-failed",
            "signature": signature,
            "optimizer": {
                "converged": False,
                "observed_calls": 1,
                "observed_retries": 0,
                "observed_max_steps": MAX_STEPS,
            },
            "detail": f"{type(exc).__name__}: {exc}",
        }
        atomic_json(directory / "receipt.json", receipt)
        raise ComputationalFailure(f"{role} optimizer failed: {exc}") from exc

    raw_path = directory / "raw-endpoint.xyz"
    atomic_xyz(raw_path, optimized.cluster)
    receipt: dict[str, Any] = {
        "schema": "calc005-component-receipt-v1",
        "status": "pending-gates",
        "signature": signature,
        "optimizer": {
            "converged": bool(optimized.converged),
            "observed_calls": 1,
            "observed_retries": 0,
            "observed_max_steps": MAX_STEPS,
        },
        "raw_endpoint": _artifact(raw_path, optimized.cluster),
    }
    atomic_json(directory / "receipt.json", receipt)
    if optimized.max_steps != MAX_STEPS:
        raise ComputationalFailure(f"{role} backend changed optimizer max_steps")
    if not optimized.converged:
        raise ComputationalFailure(f"{role} optimizer exhausted its 150-step budget")

    structure = structural_gate(role, optimized.cluster, seed, origins)
    endpoint_path = directory / "endpoint.xyz"
    atomic_xyz(endpoint_path, optimized.cluster)
    endpoint = read_xyz(endpoint_path, seed)
    structure = structural_gate(role, endpoint, seed, origins)
    endpoint_artifact = _artifact(endpoint_path, endpoint)
    gradient = gradient_gate(
        role, endpoint, backend.gradient(role, endpoint, geometry_settings)
    )
    frequency_result = backend.frequencies(role, endpoint, geometry_settings)
    frequency = _frequency_payload(role, frequency_result, endpoint, geometry_settings)
    production_energy = float(backend.energy(role, endpoint, production_settings))
    if not math.isfinite(production_energy):
        raise PhysicalStateFailure(f"{role} production single point is non-finite")
    thermochemistry = _thermochemistry(role, frequency_result, production_energy)
    receipt.update(
        {
            "status": "accepted",
            "accepted_at": now(),
            "structure": structure,
            "endpoint": endpoint_artifact,
            "gradient": gradient,
            "frequency": frequency,
            "production_single_point": {
                "method": PRODUCTION_METHOD,
                "electronic_hartree": production_energy,
                "geometry_fingerprint": frequency_geometry_fingerprint(endpoint),
                "settings_fingerprint": frequency_settings_fingerprint(
                    production_settings
                ),
                "converged": True,
            },
            "production_electronic_hartree": production_energy,
            "thermochemistry": thermochemistry,
        }
    )
    receipt["checkpoint_sha256"] = _checkpoint_sha256(receipt)
    atomic_json(directory / "receipt.json", receipt)
    receipt["cluster"] = endpoint
    receipt["exact_xyz"] = exact_xyz(endpoint)
    return receipt


def _source_map(
    deck: Path,
    pair: Calc005Pair,
    source_commit: str,
    source_provenance: dict[str, Any] | None = None,
) -> dict[str, Any]:
    repo = Path(__file__).resolve().parents[2]
    verifier = Path(__file__).with_name("calc005_si_attachment_verify.py")
    calculations = repo / "qm/CALCULATIONS.md"

    def node_value(node: tuple[int, tuple[int, int, int]]) -> list[Any]:
        return [node[0], list(node[1])]

    def origin_value(origin: Any) -> dict[str, Any]:
        return {
            "kind": origin.kind,
            "node": node_value(origin.node),
            "ordinal": origin.ordinal,
        }

    hydrolysis_origins = [
        origin_value(origin)
        for origin in pair.occupied_origins
        if origin.kind.startswith("hydrolysis-water-")
    ]
    result = {
        "source_commit": source_commit,
        "deck": str(deck),
        "deck_sha256": sha256_path(deck),
        "driver_sha256": sha256_path(Path(__file__)),
        "calc005_sha256": sha256_path(repo / "qm/quarry/calc005.py"),
        "calc005_store_sha256": sha256_path(repo / "qm/quarry/calc005_store.py"),
        "verifier_sha256": sha256_path(verifier) if verifier.is_file() else None,
        "calculations_sha256": sha256_path(calculations),
        "atom_map_sha256": atom_map_sha256(pair),
        "condensed_geometry_hash": geometry_hash(exact_xyz(pair.condensed.cluster)),
        "pilot_identity": {
            "center_node": [4, [0, 0, 0]],
            "metal_shells": 2,
            "x": pair.x,
            "y": pair.y,
            "environment_index": pair.environment_index,
            "occupied_states": pair.occupied_states,
            "vacancy_states": pair.vacancy_states,
            "occupied_frozen_origins": [
                origin_value(pair.occupied_origins[index])
                for index in pair.occupied.frozen_indices
            ],
            "vacancy_frozen_origins": [
                origin_value(pair.vacancy_origins[index])
                for index in pair.vacancy.frozen_indices
            ],
            "condensed_topology_mask": {
                "center_bridges": [
                    node_value(node) for node in pair.condensed.center_bridges
                ],
                "kept_center_bridges": [
                    node_value(node) for node in pair.condensed.kept_center_bridges
                ],
                "termination_log": pair.condensed.termination_log,
            },
            "hydrolysis_water_origins": hydrolysis_origins,
        },
    }
    if source_provenance is not None:
        result["git"] = source_provenance
    return result


def _public_component(record: dict[str, Any]) -> dict[str, Any]:
    return {
        key: value
        for key, value in record.items()
        if key not in {"cluster", "exact_xyz"}
    }


def _quarantine_canonical(output_root: Path) -> str | None:
    paths = [
        output_root / name
        for name in ("calc005-result.json", "store.sqlite", "terminal-receipt.json")
        if (output_root / name).exists()
    ]
    if not paths:
        return None
    target = output_root / "quarantine" / f"stale-{time.time_ns()}"
    target.mkdir(parents=True, exist_ok=False)
    for path in paths:
        shutil.move(path, target / path.name)
    return str(target)


def _budget(records: dict[str, dict[str, Any]]) -> dict[str, dict[str, int]]:
    result = {}
    for role in ("C", "V", "SiOH4"):
        optimizer = records.get(role, {}).get("optimizer", {})
        result[role] = {
            "calls": int(optimizer.get("observed_calls", 0)),
            "max_steps": int(optimizer.get("observed_max_steps", MAX_STEPS)),
            "retries": int(optimizer.get("observed_retries", 0)),
        }
    return result


def _failure_receipt(
    *,
    verdict: str,
    failed_component: str,
    detail: str,
    executor_identity: str,
    source_commit: str,
    source: dict[str, Any],
    records: dict[str, dict[str, Any]],
    quarantine: str | None,
) -> dict[str, Any]:
    return {
        "schema": "calc005-terminal-receipt-v1",
        "written_at": now(),
        "verdict": verdict,
        "failed_component": failed_component,
        "detail": detail,
        "executor_identity": executor_identity,
        "source_commit": source_commit,
        "source": source,
        "components": {
            role: _public_component(value) for role, value in records.items()
        },
        "optimizer_budget": _budget(records),
        "quarantine": quarantine,
        "canonical_value_exposed": False,
        "forbidden_outputs_emitted": False,
        "independent_verification_required": True,
    }


def _publish_success(
    output_root: Path,
    pair: Calc005Pair,
    records: dict[str, dict[str, Any]],
    thermochemistry: dict[str, float],
    *,
    executor_identity: str,
    source_commit: str,
    source: dict[str, Any],
    settings_receipt: dict[str, Any],
    engine: str,
    quarantine: str | None,
    execution_envelope: dict[str, Any],
) -> dict[str, Any]:
    public_components = {
        role: _public_component(record) for role, record in records.items()
    }
    result = {
        "schema": "calc005-result-v1",
        "verdict": "passed-protocol-pilot",
        "classification": CLASSIFICATION,
        "temperature_k": TEMPERATURE_K,
        "stoichiometry": STOICHIOMETRY,
        "thermochemistry": thermochemistry,
        "water_thermochemical_coefficient": 0,
    }
    generations = output_root / "generations"
    generations.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".staging-", dir=generations))
    staged_result = staging / "calc005-result.json"
    staged_calculation = staging / "calculation-receipt.json"
    staged_store = staging / "store.sqlite"
    atomic_json(staged_result, result)
    calculation = {
        "schema": "calc005-calculation-receipt-v1",
        "written_at": now(),
        "verdict": "passed-protocol-pilot",
        "classification": CLASSIFICATION,
        "temperature_k": TEMPERATURE_K,
        "executor_identity": executor_identity,
        "source_commit": source_commit,
        "source": source,
        "settings": settings_receipt,
        "engine": engine,
        "physical_state": {
            "center_node": [4, [0, 0, 0]],
            "metal_shells": 2,
            "x": 1,
            "y": 0,
            "environment_index": 1,
            "occupied_states": pair.occupied_states,
            "vacancy_states": pair.vacancy_states,
            "cycle": "Al6H38O30Si -> Al6H34O26 + H4O4Si",
        },
        "stoichiometry": STOICHIOMETRY,
        "components": public_components,
        "thermochemistry": thermochemistry,
        "optimizer_budget": _budget(records),
        "execution_envelope": execution_envelope,
        "result_sha256": sha256_path(staged_result),
        "water_thermochemical_coefficient": 0,
        "quarantine": quarantine,
        "forbidden_outputs_emitted": False,
        "independent_verification_required": True,
    }
    atomic_json(staged_calculation, calculation)
    calculation_hash = calculation_receipt_sha256(calculation)
    calculation_file_hash = sha256_path(staged_calculation)
    generation_name = f"generation-{calculation_hash[:20]}"
    generation = generations / generation_name
    try:
        write_calc005_store(staged_store, pair, records, calculation)
        validate_staged_calc005_store(staged_store, calculation, pair)
        store_hash = sha256_path(staged_store)
        if generation.exists():
            raise RuntimeError(f"immutable generation already exists: {generation}")
        staging.replace(generation)
        result_path = generation / "calc005-result.json"
        calculation_path = generation / "calculation-receipt.json"
        store_path = generation / "store.sqlite"
        receipt = {
            "schema": "calc005-terminal-receipt-v2",
            "written_at": now(),
            "verdict": "passed-protocol-pilot",
            "classification": CLASSIFICATION,
            "executor_identity": executor_identity,
            "source_commit": source_commit,
            "generation": generation_name,
            "artifacts": {
                "result": {
                    "path": str(result_path),
                    "sha256": sha256_path(result_path),
                },
                "calculation_receipt": {
                    "path": str(calculation_path),
                    "sha256": calculation_file_hash,
                    "canonical_sha256": calculation_hash,
                },
                "store": {"path": str(store_path), "sha256": store_hash},
            },
            "canonical_value_exposed": True,
            "forbidden_outputs_emitted": False,
            "independent_verification_required": True,
        }
        receipt["receipt_payload_sha256"] = receipt_payload_sha256(receipt)
        validate_calc005_store(store_path, receipt, pair)
        return receipt
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def _run_pilot_locked(
    output_root: Path,
    deck: Path,
    backend: CalculatorBackend,
    *,
    executor_identity: str = "hermes-custom-build-001",
    source_commit: str = "0000000000000000000000000000000000000000",
    environment_index: int = 1,
    use_gpu: bool = False,
    execution_envelope: dict[str, Any] | None = None,
    source_provenance: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Run only the exact C_10/V_10/SiOH4 pilot, once per component."""

    if type(environment_index) is not int or environment_index != 1:
        raise ValueError(
            "CALC-005 production authority covers only the fixed i=1 pilot"
        )
    if not isinstance(source_commit, str) or len(source_commit) != 40:
        raise ValueError("source_commit must be a 40-character Git revision")
    output_root = output_root.resolve()
    deck = deck.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    existing_terminal = output_root / "terminal-receipt.json"
    if existing_terminal.is_file():
        terminal: dict[str, Any]
        try:
            terminal = json.loads(existing_terminal.read_text())
            if not isinstance(terminal, dict):
                raise RuntimeError("terminal payload is not an object")
            verdict = terminal.get("verdict")
            if verdict == "passed-protocol-pilot":
                if terminal.get("schema") != "calc005-terminal-receipt-v2":
                    raise RuntimeError("successful terminal schema mismatch")
                if terminal.get("source_commit") != source_commit:
                    raise RuntimeError("terminal source commit mismatch")
                validate_run(output_root, deck)
            elif verdict not in {
                "incomplete-computational-failure",
                "rejected-physical-state",
            }:
                raise RuntimeError("terminal verdict is unknown")
            elif terminal.get("schema") != "calc005-terminal-receipt-v1":
                raise RuntimeError("failure terminal schema mismatch")
            return terminal
        except (
            OSError,
            json.JSONDecodeError,
            KeyError,
            TypeError,
            RuntimeError,
        ) as exc:
            quarantine = _quarantine_canonical(output_root)
            refusal_detail = (
                f"existing terminal receipt refused: {type(exc).__name__}: {exc}"
            )
            terminal = {
                "schema": "calc005-terminal-receipt-v1",
                "written_at": now(),
                "verdict": "incomplete-computational-failure",
                "failed_component": "preflight",
                "detail": refusal_detail,
                "executor_identity": executor_identity,
                "source_commit": source_commit,
                "optimizer_budget": _budget({}),
                "quarantine": quarantine,
                "canonical_value_exposed": False,
                "forbidden_outputs_emitted": False,
                "independent_verification_required": True,
            }
            atomic_json(existing_terminal, terminal)
            return terminal

    quarantine = _quarantine_canonical(output_root)
    pair = build_calc005_pair(deck, environment_index=1, metal_shells=2)
    source = _source_map(deck, pair, source_commit, source_provenance)
    geometry_settings, production_settings = calc005_settings(use_gpu=use_gpu)
    settings_receipt = {
        "geometry": asdict(geometry_settings),
        "geometry_method": R2SCAN3C_METHOD,
        "geometry_sha256": _settings_hash(geometry_settings),
        "geometry_fingerprint": frequency_settings_fingerprint(geometry_settings),
        "production": asdict(production_settings),
        "production_method": PRODUCTION_METHOD,
        "production_sha256": _settings_hash(production_settings),
        "production_fingerprint": frequency_settings_fingerprint(production_settings),
        "temperature_k": TEMPERATURE_K,
    }
    atomic_json(
        output_root / "source-settings.json",
        {"source": source, "settings": settings_receipt},
    )
    seeds = {
        "C": (pair.occupied, pair.occupied_origins),
        "V": (pair.vacancy, pair.vacancy_origins),
        "SiOH4": (pair.silicic_acid, pair.silicic_acid_origins),
    }
    records: dict[str, dict[str, Any]] = {}
    failed_component = "preflight"

    def retain_failed_component_receipt() -> None:
        path = output_root / "components" / failed_component / "receipt.json"
        if not path.is_file():
            return
        try:
            records[failed_component] = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            # The terminal still records the spent component and refuses replay.
            records[failed_component] = {
                "status": "invalid-stage-evidence",
                "optimizer": {
                    "observed_calls": 1,
                    "observed_retries": 0,
                    "observed_max_steps": MAX_STEPS,
                },
            }

    try:
        for role in ("C", "V", "SiOH4"):
            failed_component = role
            records[role] = _run_component(
                output_root,
                role,
                *seeds[role],
                source,
                geometry_settings,
                production_settings,
                backend,
            )
        thermochemistry = compose_thermochemistry(
            {role: records[role]["thermochemistry"] for role in STOICHIOMETRY}
        )
    except ComputationalFailure as exc:
        retain_failed_component_receipt()
        receipt = _failure_receipt(
            verdict="incomplete-computational-failure",
            failed_component=failed_component,
            detail=str(exc),
            executor_identity=executor_identity,
            source_commit=source_commit,
            source=source,
            records=records,
            quarantine=quarantine or "spent-checkpoint-retained-in-place",
        )
        atomic_json(existing_terminal, receipt)
        return receipt
    except (PhysicalStateFailure, ValueError, FloatingPointError) as exc:
        component_receipt = (
            output_root / "components" / failed_component / "receipt.json"
        )
        if component_receipt.is_file():
            partial = json.loads(component_receipt.read_text())
            partial.update(
                {
                    "status": "gate-rejected",
                    "rejected_at": now(),
                    "detail": f"{type(exc).__name__}: {exc}",
                }
            )
            atomic_json(component_receipt, partial)
        retain_failed_component_receipt()
        receipt = _failure_receipt(
            verdict="rejected-physical-state",
            failed_component=failed_component,
            detail=str(exc),
            executor_identity=executor_identity,
            source_commit=source_commit,
            source=source,
            records=records,
            quarantine=quarantine,
        )
        atomic_json(existing_terminal, receipt)
        return receipt
    except Exception as exc:
        component_receipt = (
            output_root / "components" / failed_component / "receipt.json"
        )
        if component_receipt.is_file():
            try:
                partial = json.loads(component_receipt.read_text())
            except (OSError, json.JSONDecodeError):
                partial = {}
            partial.update(
                {
                    "status": "calculator-failed-after-optimizer",
                    "failed_at": now(),
                    "detail": f"{type(exc).__name__}: {exc}",
                }
            )
            atomic_json(component_receipt, partial)
        retain_failed_component_receipt()
        receipt = _failure_receipt(
            verdict="incomplete-computational-failure",
            failed_component=failed_component,
            detail=f"{type(exc).__name__}: {exc}",
            executor_identity=executor_identity,
            source_commit=source_commit,
            source=source,
            records=records,
            quarantine=quarantine or "spent-checkpoint-retained-in-place",
        )
        atomic_json(existing_terminal, receipt)
        return receipt

    try:
        receipt = _publish_success(
            output_root,
            pair,
            records,
            thermochemistry,
            executor_identity=executor_identity,
            source_commit=source_commit,
            source=source,
            settings_receipt=settings_receipt,
            engine=getattr(backend, "name", type(backend).__name__),
            quarantine=quarantine,
            execution_envelope=execution_envelope or {"mode": "in-process-test"},
        )
    except Exception as exc:
        publication_quarantine = _quarantine_canonical(output_root)
        receipt = _failure_receipt(
            verdict="incomplete-computational-failure",
            failed_component="evidence-graph-publication",
            detail=f"{type(exc).__name__}: {exc}",
            executor_identity=executor_identity,
            source_commit=source_commit,
            source=source,
            records=records,
            quarantine=publication_quarantine or quarantine,
        )
    atomic_json(existing_terminal, receipt)
    return receipt


@contextmanager
def exclusive_run(output_root: Path):
    """Prevent concurrent processes from double-spending a component budget."""

    lock_path = output_root.resolve() / "run.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError("CALC-005 pilot is already active") from exc
        try:
            yield
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def run_pilot(
    output_root: Path,
    deck: Path,
    backend: CalculatorBackend,
    *,
    executor_identity: str = "hermes-custom-build-001",
    source_commit: str = "0000000000000000000000000000000000000000",
    environment_index: int = 1,
    use_gpu: bool = False,
    execution_envelope: dict[str, Any] | None = None,
    source_provenance: dict[str, Any] | None = None,
) -> dict[str, Any]:
    if type(environment_index) is not int or environment_index != 1:
        raise ValueError(
            "CALC-005 production authority covers only the fixed i=1 pilot"
        )
    if (
        not isinstance(source_commit, str)
        or len(source_commit) != 40
        or any(character not in "0123456789abcdef" for character in source_commit)
    ):
        raise ValueError("source_commit must be a 40-character lowercase Git revision")
    with exclusive_run(output_root):
        return _run_pilot_locked(
            output_root,
            deck,
            backend,
            executor_identity=executor_identity,
            source_commit=source_commit,
            environment_index=environment_index,
            use_gpu=use_gpu,
            execution_envelope=execution_envelope,
            source_provenance=source_provenance,
        )


def validate_run(output_root: Path, deck: Path) -> dict[str, Any]:
    """Validate the atomic terminal signal, immutable receipt, Store, and source."""

    output_root = output_root.resolve()
    terminal_path = output_root / "terminal-receipt.json"
    receipt = json.loads(terminal_path.read_text())
    verdict = receipt.get("verdict")
    if verdict != "passed-protocol-pilot":
        if receipt.get("schema") == "calc005-terminal-receipt-v1" and verdict in {
            "incomplete-computational-failure",
            "rejected-physical-state",
        }:
            return {
                "schema": "calc005-validator-v1",
                "status": "invalid",
                "verdict": verdict,
                "detail": receipt.get("detail"),
                "optimizer_calls": 0,
            }
        raise RuntimeError(
            "CALC-005 terminal receipt is malformed or has unknown verdict"
        )
    if receipt.get("schema") != "calc005-terminal-receipt-v2":
        raise RuntimeError("CALC-005 successful terminal schema mismatch")
    pair = build_calc005_pair(deck.resolve(), environment_index=1, metal_shells=2)
    calculation_artifact = receipt.get("artifacts", {}).get("calculation_receipt", {})
    calculation_path = Path(str(calculation_artifact.get("path", "")))
    calculation = json.loads(calculation_path.read_text())
    calculation_source = dict(calculation.get("source", {}))
    git_provenance = calculation_source.pop("git", None)
    if calculation_source != _source_map(
        deck.resolve(), pair, receipt.get("source_commit")
    ):
        raise RuntimeError("CALC-005 source/deck/atom-map hash drifted")
    if git_provenance is not None and (
        git_provenance.get("clean") is not True
        or git_provenance.get("head") != receipt.get("source_commit")
        or git_provenance.get("origin_head") != receipt.get("source_commit")
        or not git_provenance.get("branch")
    ):
        raise RuntimeError("CALC-005 committed source provenance drifted")
    store_path = Path(receipt["artifacts"]["store"]["path"])
    store = validate_calc005_store(store_path, receipt, pair)
    return {
        "schema": "calc005-validator-v1",
        "status": "valid",
        "verdict": receipt["verdict"],
        "receipt_sha256": sha256_path(terminal_path),
        "calculation_receipt_sha256": calculation_artifact["canonical_sha256"],
        "store": store,
        "optimizer_calls": 0,
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    probe = subparsers.add_parser("probe")
    probe.add_argument("--deck", type=Path, required=True)
    probe.add_argument("--repeat", type=int, default=2)
    probe.add_argument("--output", type=Path)

    run = subparsers.add_parser("run")
    run.add_argument("--deck", type=Path, required=True)
    run.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    run.add_argument("--executor-identity", default="hermes-custom-build-001")
    run.add_argument("--gpu", action="store_true")

    validate = subparsers.add_parser("validate")
    validate.add_argument("--deck", type=Path, required=True)
    validate.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    for subparser in (probe, run, validate):
        subparser.add_argument("--threads", type=int, default=16)
        subparser.add_argument("--nice", type=int, default=10)
        subparser.add_argument("--log")
    return parser


def _git(repo: Path, *arguments: str) -> str:
    return subprocess.run(
        ["git", *arguments],
        cwd=repo,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def verify_production_source(repo: Path, deck: Path) -> dict[str, Any]:
    """Require a clean, committed, origin-bound branch and tracked source set."""

    repo = repo.resolve()
    expected_repo = Path(__file__).resolve().parents[2]
    if repo != expected_repo:
        raise RuntimeError("production source must be this driver's repository")
    deck = deck.resolve()
    try:
        deck_relative = deck.relative_to(repo)
    except ValueError as exc:
        raise RuntimeError(
            "production deck must be inside the source repository"
        ) from exc
    dirty = _git(repo, "status", "--porcelain=v1", "--untracked-files=all")
    if dirty:
        raise RuntimeError("production CALC-005 refuses a dirty source worktree")
    branch = _git(repo, "branch", "--show-current")
    if not branch:
        raise RuntimeError("production CALC-005 refuses detached HEAD")
    head = _git(repo, "rev-parse", "HEAD")
    origin_head = _git(repo, "rev-parse", f"origin/{branch}")
    if head != origin_head:
        raise RuntimeError("production HEAD must exactly equal origin/<current-branch>")
    required = [
        deck_relative.as_posix(),
        "qm/CALCULATIONS.md",
        "qm/quarry/calc005.py",
        "qm/quarry/calc005_store.py",
        "qm/scripts/calc005_si_attachment.py",
        "qm/scripts/calc005_si_attachment_verify.py",
    ]
    _git(repo, "ls-files", "--error-unmatch", "--", *required)
    for relative in required:
        _git(repo, "cat-file", "-e", f"HEAD:{relative}")
    return {
        "repository": str(repo),
        "branch": branch,
        "head": head,
        "origin_head": origin_head,
        "clean": True,
        "tracked_source_files": required,
    }


def _timespan_seconds(value: str) -> int:
    value = value.strip()
    if value.isdigit():
        return int(value) // 1_000_000
    units = {"d": 86400, "h": 3600, "min": 60, "s": 1}
    total = 0.0
    remaining = value
    while remaining:
        for suffix in ("min", "d", "h", "s"):
            marker = remaining.find(suffix)
            if marker > 0:
                total += float(remaining[:marker]) * units[suffix]
                remaining = remaining[marker + len(suffix) :].strip()
                break
        else:
            raise RuntimeError(f"unrecognized systemd timespan: {value}")
    return int(total)


def measure_production_envelope() -> dict[str, Any]:
    """Read and fail closed on the exact live cgroup/process/QI2 lease envelope."""

    cgroup_lines = Path("/proc/self/cgroup").read_text().splitlines()
    unified = [line.split("::", 1)[1] for line in cgroup_lines if "::" in line]
    if len(unified) != 1:
        raise RuntimeError("CALC-005 requires one cgroup-v2 systemd unit")
    relative = unified[0].lstrip("/")
    cgroup = Path("/sys/fs/cgroup") / relative
    units = [
        part for part in Path(relative).parts if part.endswith((".service", ".scope"))
    ]
    if not units:
        raise RuntimeError("CALC-005 must run inside a systemd service/scope")
    unit = units[-1]

    def integer_limit(name: str) -> int:
        value = (cgroup / name).read_text().strip()
        if value == "max":
            raise RuntimeError(f"CALC-005 requires a finite {name}")
        return int(value)

    memory_max = integer_limit("memory.max")
    swap_max = integer_limit("memory.swap.max")
    cpu_fields = (cgroup / "cpu.max").read_text().split()
    if len(cpu_fields) != 2 or cpu_fields[0] == "max":
        raise RuntimeError("CALC-005 requires a finite cpu.max quota")
    cpu_percent = int(round(int(cpu_fields[0]) / int(cpu_fields[1]) * 100))
    runtime_text = subprocess.run(
        ["systemctl", "show", unit, "--property=RuntimeMaxUSec", "--value"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    runtime_seconds = _timespan_seconds(runtime_text)
    niceness = os.getpriority(os.PRIO_PROCESS, 0)
    thread_values = {
        name: int(os.environ.get(name, "0"))
        for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")
    }
    lease_path = Path(
        os.environ.get(
            "GPU_LEASE_PATH", str(Path.home() / ".local/state/gpu-lease/lease.json")
        )
    ).expanduser()
    lease = json.loads(lease_path.read_text())
    expected = {
        "runtime": runtime_seconds == EXPECTED_RUNTIME_MAX_SECONDS,
        "memory": memory_max == EXPECTED_MEMORY_MAX_BYTES,
        "swap": swap_max == EXPECTED_MEMORY_SWAP_MAX_BYTES,
        "cpu": cpu_percent == EXPECTED_CPU_QUOTA_PERCENT,
        "nice": niceness == EXPECTED_NICE,
        "threads": all(value == EXPECTED_THREADS for value in thread_values.values()),
        "lease_owner": lease.get("owner") == EXPECTED_GPU_OWNER,
        "lease_pid": lease.get("pid") == os.getpid(),
        "lease_ttl": lease.get("ttl") == GPU_TTL_HOURS,
        "lease_memory": lease.get("expected_gb") == 16.0,
    }
    failed = [name for name, passed in expected.items() if not passed]
    if failed:
        raise RuntimeError(
            "production execution envelope is not exact: " + ", ".join(failed)
        )
    return {
        "measured": True,
        "systemd_unit": unit,
        "cgroup_path": str(cgroup),
        "runtime_max_seconds": runtime_seconds,
        "memory_max_bytes": memory_max,
        "memory_swap_max_bytes": swap_max,
        "cpu_quota_percent": cpu_percent,
        "nice": niceness,
        "thread_environment": thread_values,
        "qi2_lease": {
            "path": str(lease_path),
            "owner": lease["owner"],
            "pid": lease["pid"],
            "ttl_hours": lease["ttl"],
            "expected_gb": lease["expected_gb"],
            "maximum_gb": 18.0,
        },
        "shared_service_mutation": False,
        "restoration": "pending-process-exit",
    }


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.threads < 1 or args.threads > 16:
        raise SystemExit("--threads must be within 1..16")
    if args.nice < 10:
        raise SystemExit("--nice must be >=10")
    if args.command == "probe":
        if args.repeat != 2:
            raise SystemExit("--repeat must be exactly 2 for the fixed proof")
        payload = probe_calc005_pairs(args.deck)
        payload["driver_sha256"] = sha256_path(Path(__file__))
        content = json.dumps(payload, indent=2, sort_keys=True) + "\n"
        if args.output is not None:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(content)
        print(content, end="")
        return 0
    if args.command == "validate":
        try:
            payload = validate_run(args.output_root, args.deck)
        except Exception as exc:
            payload = {
                "schema": "calc005-validator-v1",
                "status": "invalid",
                "detail": f"{type(exc).__name__}: {exc}",
                "optimizer_calls": 0,
            }
        print(json.dumps(payload, indent=2, sort_keys=True))
        return 0 if payload["status"] == "valid" else 1
    if not args.gpu:
        raise SystemExit("production CALC-005 run requires --gpu")
    if args.threads != EXPECTED_THREADS or args.nice != EXPECTED_NICE:
        raise SystemExit("production CALC-005 requires exactly --threads 16 --nice 10")
    repo = Path(__file__).resolve().parents[2]
    source_provenance = verify_production_source(repo, args.deck)
    execution_envelope = measure_production_envelope()
    result = run_pilot(
        args.output_root,
        args.deck,
        PipelineBackend(),
        executor_identity=args.executor_identity,
        source_commit=source_provenance["head"],
        use_gpu=True,
        execution_envelope=execution_envelope,
        source_provenance=source_provenance,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["verdict"] == "passed-protocol-pilot" else 2


if __name__ == "__main__":
    raise SystemExit(main())
