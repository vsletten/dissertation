"""Immutable evidence graph and validators for the CALC-005 protocol pilot."""

from __future__ import annotations

import hashlib
import json
import math
import os
import sqlite3
import tempfile
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np

from quarry.calc005 import STOICHIOMETRY, Calc005Pair, atom_map_sha256
from quarry.pipeline import (
    HARTREE_TO_KJ,
    DftSettings,
    frequency_geometry_fingerprint,
    frequency_settings_fingerprint,
)
from quarry.rates import surface_thermo_from_frequencies, thermo_from_frequencies
from quarry.store import Store, geometry_hash

R2SCAN3C_METHOD = "r2scan-3c/def2-mtzvpp/d4/gcp"
PRODUCTION_METHOD = "wb97m-v/def2-tzvpd/smd(water)"
CLASSIFICATION = "thermodynamic/non-kinetic/non-emittable"
LABEL = (
    "S_10 relative environment stabilization; thermodynamic/non-kinetic/non-emittable"
)
MAX_STEPS = 150
GRADIENT_RMS_MAX = 3.0e-4
GRADIENT_MAX_MAX = 4.5e-4
IMAGINARY_NOISE_FLOOR_CM = 30.0
TEMPERATURE_K = 298.15
STANDARD_STATE_CORRECTION_KJ_MOL = 7.958500693927389
R_KJ_MOL_K = 0.00831446261815324
MIN_PAIR_DISTANCE_A = 0.75
MAX_OWNER_DISTANCE_A = 1.25
MIN_OWNER_MARGIN_A = 0.15
GPU_TTL_HOURS = 13.0
EXPECTED_MEMORY_MAX_BYTES = 32 * 1024**3
EXPECTED_MEMORY_SWAP_MAX_BYTES = 4 * 1024**3
EXPECTED_CPU_QUOTA_PERCENT = 1600
EXPECTED_NICE = 10
EXPECTED_THREADS = 16
EXPECTED_RUNTIME_MAX_SECONDS = 43200
EXPECTED_GPU_OWNER = "calc005_si_attachment"

_EXTENSION_SCHEMA = """
CREATE TABLE calc005_cycles (
    id INTEGER PRIMARY KEY CHECK (id = 1),
    schema TEXT NOT NULL,
    x INTEGER NOT NULL,
    y INTEGER NOT NULL,
    environment_index INTEGER NOT NULL,
    center_state_from INTEGER NOT NULL,
    center_state_to INTEGER NOT NULL,
    osa_state_from INTEGER NOT NULL,
    osa_state_to INTEGER NOT NULL,
    temperature_k REAL NOT NULL,
    surface_standard_state TEXT NOT NULL,
    solute_standard_state TEXT NOT NULL,
    status TEXT NOT NULL,
    label TEXT NOT NULL,
    atom_map_sha256 TEXT NOT NULL,
    calculation_receipt_sha256 TEXT NOT NULL,
    condensed_structure_id INTEGER NOT NULL REFERENCES structures(id)
);
CREATE TABLE calc005_components (
    cycle_id INTEGER NOT NULL REFERENCES calc005_cycles(id),
    role TEXT NOT NULL UNIQUE,
    coefficient INTEGER NOT NULL,
    structure_id INTEGER NOT NULL REFERENCES structures(id),
    opt_job_id INTEGER NOT NULL UNIQUE REFERENCES jobs(id),
    freq_job_id INTEGER NOT NULL UNIQUE REFERENCES jobs(id),
    sp_job_id INTEGER NOT NULL UNIQUE REFERENCES jobs(id),
    raw_endpoint_sha256 TEXT NOT NULL,
    endpoint_sha256 TEXT NOT NULL,
    checkpoint_sha256 TEXT NOT NULL,
    PRIMARY KEY (cycle_id, role)
);
CREATE TABLE calc005_construction (
    cycle_id INTEGER NOT NULL REFERENCES calc005_cycles(id),
    c_atom_index INTEGER NOT NULL,
    element TEXT NOT NULL,
    origin_kind TEXT NOT NULL,
    origin_node TEXT NOT NULL,
    origin_ordinal INTEGER,
    product_role TEXT NOT NULL,
    product_atom_index INTEGER NOT NULL,
    thermochemical INTEGER NOT NULL CHECK (thermochemical = 0),
    PRIMARY KEY (cycle_id, c_atom_index)
);
"""


def _canonical_hash(payload: dict[str, Any]) -> str:
    return hashlib.sha256(
        json.dumps(
            payload, separators=(",", ":"), sort_keys=True, allow_nan=False
        ).encode()
    ).hexdigest()


def receipt_payload_sha256(receipt: dict[str, Any]) -> str:
    """Hash the complete terminal payload except its unavoidable self-hash field."""

    payload = json.loads(json.dumps(receipt, allow_nan=False))
    payload.pop("receipt_payload_sha256", None)
    return _canonical_hash(payload)


def calculation_receipt_sha256(receipt: dict[str, Any]) -> str:
    """Hash an immutable calculation receipt with no Store/terminal back-edge."""

    return _canonical_hash(receipt)


def _checkpoint_sha256(record: dict[str, Any]) -> str:
    payload = {
        key: value
        for key, value in record.items()
        if key not in {"checkpoint_sha256", "accepted_at", "cluster", "exact_xyz"}
    }
    return _canonical_hash(payload)


def _finite(value: Any) -> bool:
    return (
        not isinstance(value, bool)
        and isinstance(value, (int, float))
        and math.isfinite(float(value))
    )


def _close(left: Any, right: Any, *, tolerance: float = 1.0e-10) -> bool:
    return (
        _finite(left)
        and _finite(right)
        and math.isclose(float(left), float(right), rel_tol=1.0e-12, abs_tol=tolerance)
    )


def _sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _node_text(origin: Any) -> str:
    return json.dumps(origin.node, separators=(",", ":"))


def _expected_modes(role: str, pair: Calc005Pair) -> int:
    cluster = {"C": pair.occupied, "V": pair.vacancy, "SiOH4": pair.silicic_acid}[role]
    if role == "SiOH4":
        return 3 * len(cluster.symbols) - 6
    return 3 * (len(cluster.symbols) - len(cluster.frozen_indices))


def _heavy_edges(coords: np.ndarray, symbols: list[str]) -> set[tuple[int, int]]:
    return {
        (cation, oxygen)
        for cation, symbol in enumerate(symbols)
        if symbol in {"Al", "Si"}
        for oxygen, oxygen_symbol in enumerate(symbols)
        if oxygen_symbol == "O"
        and float(np.linalg.norm(coords[cation] - coords[oxygen])) <= 2.35
    }


def _owner_map(coords: np.ndarray, symbols: list[str]) -> dict[int, int]:
    oxygens = [index for index, symbol in enumerate(symbols) if symbol == "O"]
    owners = {}
    for hydrogen, symbol in enumerate(symbols):
        if symbol != "H":
            continue
        distances = sorted(
            (float(np.linalg.norm(coords[hydrogen] - coords[oxygen])), oxygen)
            for oxygen in oxygens
        )
        if (
            len(distances) < 2
            or distances[0][0] > MAX_OWNER_DISTANCE_A
            or distances[1][0] - distances[0][0] < MIN_OWNER_MARGIN_A
        ):
            raise RuntimeError("proton owner distance/margin gate failed")
        owners[hydrogen] = distances[0][1]
    return owners


def _validate_endpoint_structure(
    role: str,
    path: Path,
    pair: Calc005Pair,
) -> np.ndarray:
    template = {"C": pair.occupied, "V": pair.vacancy, "SiOH4": pair.silicic_acid}[role]
    origins = {
        "C": pair.occupied_origins,
        "V": pair.vacancy_origins,
        "SiOH4": pair.silicic_acid_origins,
    }[role]
    coords = _read_xyz_coords(path, template.symbols)
    if template.frozen_indices and not np.array_equal(
        coords[template.frozen_indices], template.coords[template.frozen_indices]
    ):
        raise RuntimeError(f"{role} frozen coordinates drifted")
    delta = coords[:, None, :] - coords[None, :, :]
    distances = np.linalg.norm(delta, axis=2)
    distances[np.diag_indices_from(distances)] = np.inf
    if float(np.min(distances)) < MIN_PAIR_DISTANCE_A:
        raise RuntimeError(f"{role} collision gate failed")
    if _owner_map(coords, template.symbols) != _owner_map(
        template.coords, template.symbols
    ):
        raise RuntimeError(f"{role} proton owner changed")
    expected_edges = _heavy_edges(template.coords, template.symbols)
    observed_edges = _heavy_edges(coords, template.symbols)
    if observed_edges != expected_edges:
        raise RuntimeError(f"{role} heavy connectivity drifted")
    center = [
        index
        for index, origin in enumerate(origins)
        if origin.kind == "deck" and origin.node == (4, (0, 0, 0))
    ]
    center_bonds = {
        oxygen for cation, oxygen in observed_edges if center and cation == center[0]
    }
    if role in {"C", "SiOH4"} and (len(center) != 1 or len(center_bonds) != 4):
        raise RuntimeError(f"{role} exact four center Si-O bonds failed")
    if role == "V" and center:
        raise RuntimeError("V retains the removed center Si")
    return coords


def _compose_thermochemistry(
    components: dict[str, dict[str, float]],
) -> dict[str, float]:
    def delta(key: str) -> float:
        return math.fsum(
            STOICHIOMETRY[role] * float(components[role][key])
            for role in ("C", "V", "SiOH4")
        )

    electronic = delta("electronic_kj_mol")
    zpe = delta("zpe_kj_mol")
    thermal = delta("thermal_kj_mol")
    minus_t_delta_s = -TEMPERATURE_K * delta("entropy_kj_mol_k")
    total = math.fsum(
        (electronic, zpe, thermal, minus_t_delta_s, STANDARD_STATE_CORRECTION_KJ_MOL)
    )
    independent = (
        math.fsum(
            STOICHIOMETRY[role]
            * math.fsum(
                (
                    float(components[role]["electronic_kj_mol"]),
                    float(components[role]["zpe_kj_mol"]),
                    float(components[role]["thermal_kj_mol"]),
                    -TEMPERATURE_K * float(components[role]["entropy_kj_mol_k"]),
                )
            )
            for role in ("C", "V", "SiOH4")
        )
        + STANDARD_STATE_CORRECTION_KJ_MOL
    )
    return {
        "electronic_kj_mol": electronic,
        "zpe_kj_mol": zpe,
        "thermal_kj_mol": thermal,
        "minus_t_delta_s_kj_mol": minus_t_delta_s,
        "standard_state_kj_mol": STANDARD_STATE_CORRECTION_KJ_MOL,
        "s_10_kj_mol": total,
        "independent_s_10_kj_mol": independent,
        "recomposition_error_kj_mol": abs(total - independent),
    }


def _validate_component_record(
    role: str,
    record: dict[str, Any],
    pair: Calc005Pair,
    settings: dict[str, Any],
) -> None:
    """Enforce gates from raw fields, never from an aggregate executor verdict."""

    if record.get("schema") != "calc005-component-receipt-v1":
        raise RuntimeError(f"{role} component schema drifted")
    if record.get("status") != "accepted":
        raise RuntimeError(f"{role} component status is not accepted")
    optimizer = record.get("optimizer")
    if optimizer != {
        "converged": True,
        "observed_calls": 1,
        "observed_retries": 0,
        "observed_max_steps": MAX_STEPS,
    }:
        raise RuntimeError(f"{role} optimizer budget/convergence drifted")
    if record.get("checkpoint_sha256") != _checkpoint_sha256(record):
        raise RuntimeError(f"{role} checkpoint hash drifted")

    gradient = record.get("gradient", {})
    rms = gradient.get("rms_hartree_per_bohr")
    maximum = gradient.get("max_hartree_per_bohr")
    if (
        gradient.get("passed") is not True
        or gradient.get("rms_threshold_hartree_per_bohr") != GRADIENT_RMS_MAX
        or gradient.get("max_threshold_hartree_per_bohr") != GRADIENT_MAX_MAX
        or not _finite(rms)
        or not _finite(maximum)
        or float(rms) > GRADIENT_RMS_MAX
        or float(maximum) > GRADIENT_MAX_MAX
    ):
        raise RuntimeError(f"{role} gradient gate failed")

    frequency = record.get("frequency", {})
    real = frequency.get("frequencies_cm")
    imaginary = frequency.get("imaginary_cm")
    expected = _expected_modes(role, pair)
    if not isinstance(real, list) or not isinstance(imaginary, list):
        raise RuntimeError(f"{role} frequency arrays are malformed")
    if (
        frequency.get("passed") is not True
        or frequency.get("hessian") != ("full" if role == "SiOH4" else "PHVA")
        or frequency.get("noise_floor_cm") != IMAGINARY_NOISE_FLOOR_CM
        or frequency.get("expected_mode_count") != expected
        or frequency.get("real_mode_count") != len(real)
        or frequency.get("imaginary_mode_count") != len(imaginary)
        or frequency.get("observed_mode_count") != len(real) + len(imaginary)
        or len(real) + len(imaginary) != expected
        or any(not _finite(value) or float(value) <= 0.0 for value in real)
        or any(
            not _finite(value)
            or float(value) < 0.0
            or float(value) > IMAGINARY_NOISE_FLOOR_CM
            for value in imaginary
        )
    ):
        raise RuntimeError(f"{role} frequency gate/count failed")
    rotational = frequency.get("rotational_temperatures_k")
    if role == "SiOH4":
        if (
            frequency.get("linear") is not False
            or not isinstance(rotational, list)
            or len(rotational) != 3
            or any(not _finite(value) or float(value) <= 0.0 for value in rotational)
        ):
            raise RuntimeError("SiOH4 nonlinear rotational evidence drifted")
    elif rotational is not None:
        raise RuntimeError(f"{role} PHVA carries forbidden rotations")
    if (
        not _finite(frequency.get("molar_mass_kg"))
        or float(frequency["molar_mass_kg"]) <= 0.0
    ):
        raise RuntimeError(f"{role} frequency molar mass is malformed")

    template = {"C": pair.occupied, "V": pair.vacancy, "SiOH4": pair.silicic_acid}[role]
    endpoint = record.get("endpoint", {})
    raw = record.get("raw_endpoint", {})
    for label, artifact in (("endpoint", endpoint), ("raw endpoint", raw)):
        path = Path(str(artifact.get("path", "")))
        if not path.is_file() or artifact.get("sha256") != _sha256_path(path):
            raise RuntimeError(f"{role} {label} artifact hash drifted")
        _validate_endpoint_structure(role, path, pair)
    endpoint_coords = _validate_endpoint_structure(role, Path(endpoint["path"]), pair)
    if endpoint.get("geometry_fingerprint") != frequency_geometry_fingerprint(
        template.__class__(
            name=template.name,
            symbols=template.symbols,
            coords=endpoint_coords,
            charge=template.charge,
            spin=template.spin,
            frozen_indices=template.frozen_indices,
            site_family=template.site_family,
            note=template.note,
        )
    ):
        raise RuntimeError(f"{role} endpoint geometry fingerprint drifted")
    if frequency.get("geometry_fingerprint") != endpoint.get("geometry_fingerprint"):
        raise RuntimeError(f"{role} frequency is not bound to endpoint")
    if frequency.get("settings_fingerprint") != settings.get("geometry_fingerprint"):
        raise RuntimeError(f"{role} frequency settings fingerprint drifted")

    single_point = record.get("production_single_point", {})
    energy = record.get("production_electronic_hartree")
    if (
        single_point.get("method") != PRODUCTION_METHOD
        or single_point.get("converged") is not True
        or not _close(single_point.get("electronic_hartree"), energy)
        or not _finite(energy)
        or single_point.get("geometry_fingerprint")
        != endpoint.get("geometry_fingerprint")
        or single_point.get("settings_fingerprint")
        != settings.get("production_fingerprint")
    ):
        raise RuntimeError(f"{role} production single point failed convergence/binding")
    thermochemistry = record.get("thermochemistry", {})
    required_terms = {
        "electronic_kj_mol",
        "zpe_kj_mol",
        "thermal_kj_mol",
        "entropy_kj_mol_k",
    }
    if set(thermochemistry) != required_terms or any(
        not _finite(value) for value in thermochemistry.values()
    ):
        raise RuntimeError(f"{role} thermochemistry is malformed")
    if not _close(
        thermochemistry["electronic_kj_mol"],
        float(energy) * HARTREE_TO_KJ,
        tolerance=1.0e-8,
    ):
        raise RuntimeError(f"{role} production energy conversion drifted")
    electronic_kj = float(energy) * HARTREE_TO_KJ
    if role == "SiOH4":
        independent = thermo_from_frequencies(
            electronic_kj,
            [float(value) for value in real],
            TEMPERATURE_K,
            molar_mass_kg=float(frequency["molar_mass_kg"]),
            rotational_temperatures_k=[float(value) for value in rotational],
            linear=False,
        )
    else:
        independent = surface_thermo_from_frequencies(
            electronic_kj,
            [float(value) for value in real],
            TEMPERATURE_K,
        )
    independent_terms = {
        "electronic_kj_mol": electronic_kj,
        "zpe_kj_mol": independent.zpe_kj,
        "thermal_kj_mol": independent.thermal_kj,
        "entropy_kj_mol_k": independent.entropy_kj_per_k,
    }
    if any(
        not _close(thermochemistry[key], value, tolerance=1.0e-8)
        for key, value in independent_terms.items()
    ):
        raise RuntimeError(f"{role} thermochemistry does not recompute from raw modes")


def _read_xyz_coords(path: Path, symbols: list[str]) -> np.ndarray:
    lines = path.read_text().splitlines()
    if len(lines) != len(symbols) + 2 or int(lines[0]) != len(symbols):
        raise RuntimeError(f"malformed XYZ artifact: {path}")
    rows = [line.split() for line in lines[2:]]
    if any(len(row) != 4 for row in rows) or [row[0] for row in rows] != symbols:
        raise RuntimeError(f"XYZ atom identity/order drifted: {path}")
    coords = np.asarray([[float(value) for value in row[1:]] for row in rows])
    if not np.all(np.isfinite(coords)):
        raise RuntimeError(f"XYZ coordinates are non-finite: {path}")
    return coords


def _validate_execution_envelope(envelope: Any, backend_kind: Any) -> None:
    if envelope == {"mode": "in-process-test"}:
        if backend_kind != "test":
            raise RuntimeError(
                "production backend cannot use an in-process-test execution envelope"
            )
        return
    if backend_kind != "production":
        raise RuntimeError("test backend cannot claim a production execution envelope")
    if not isinstance(envelope, dict):
        raise RuntimeError("production execution envelope is missing")
    threads = envelope.get("thread_environment")
    lease = envelope.get("qi2_lease")
    if (
        envelope.get("measured") is not True
        or not str(envelope.get("systemd_unit", "")).endswith((".service", ".scope"))
        or not str(envelope.get("cgroup_path", "")).startswith("/sys/fs/cgroup/")
        or envelope.get("runtime_max_seconds") != EXPECTED_RUNTIME_MAX_SECONDS
        or envelope.get("memory_max_bytes") != EXPECTED_MEMORY_MAX_BYTES
        or envelope.get("memory_swap_max_bytes") != EXPECTED_MEMORY_SWAP_MAX_BYTES
        or envelope.get("cpu_quota_percent") != EXPECTED_CPU_QUOTA_PERCENT
        or envelope.get("nice") != EXPECTED_NICE
        or threads
        != {
            "OMP_NUM_THREADS": EXPECTED_THREADS,
            "MKL_NUM_THREADS": EXPECTED_THREADS,
            "OPENBLAS_NUM_THREADS": EXPECTED_THREADS,
        }
        or not isinstance(lease, dict)
        or lease.get("owner") != EXPECTED_GPU_OWNER
        or type(lease.get("pid")) is not int
        or lease["pid"] <= 0
        or lease.get("ttl_hours") != GPU_TTL_HOURS
        or lease.get("expected_gb") != 16.0
        or lease.get("maximum_gb") != 18.0
        or not str(lease.get("path", ""))
        or envelope.get("shared_service_mutation") is not False
        or envelope.get("restoration") != "pending-process-exit"
    ):
        raise RuntimeError("production execution envelope evidence drifted")


def validate_calculation_receipt(receipt: dict[str, Any], pair: Calc005Pair) -> None:
    """Independently validate immutable pre-Store calculation content."""

    if (
        receipt.get("schema") != "calc005-calculation-receipt-v1"
        or receipt.get("verdict") != "passed-protocol-pilot"
        or receipt.get("classification") != CLASSIFICATION
        or receipt.get("temperature_k") != TEMPERATURE_K
        or receipt.get("stoichiometry") != STOICHIOMETRY
        or receipt.get("water_thermochemical_coefficient") != 0
        or receipt.get("forbidden_outputs_emitted") is not False
        or receipt.get("independent_verification_required") is not True
    ):
        raise RuntimeError("calculation receipt identity/classification drifted")
    _validate_execution_envelope(
        receipt.get("execution_envelope"), receipt.get("backend_kind")
    )
    components = receipt.get("components")
    if not isinstance(components, dict) or set(components) != set(STOICHIOMETRY):
        raise RuntimeError("calculation receipt component set drifted")
    settings = receipt.get("settings", {})
    geometry = settings.get("geometry")
    production = settings.get("production")
    if not isinstance(geometry, dict) or not isinstance(production, dict):
        raise RuntimeError("calculation settings payload is malformed")
    use_gpu = geometry.get("use_gpu")
    if type(use_gpu) is not bool or production.get("use_gpu") is not use_gpu:
        raise RuntimeError("calculation CPU/GPU settings disagree")
    expected_geometry = DftSettings(
        xc="r2scan",
        basis="def2-mtzvpp",
        composite="r2scan3c",
        density_fit=True,
        use_gpu=use_gpu,
    )
    expected_production = DftSettings(
        xc="wb97m-v",
        basis="def2-tzvpd",
        solvent="smd",
        density_fit=True,
        use_gpu=use_gpu,
    )
    if (
        geometry != asdict(expected_geometry)
        or production != asdict(expected_production)
        or settings.get("geometry_method") != R2SCAN3C_METHOD
        or settings.get("production_method") != PRODUCTION_METHOD
        or settings.get("temperature_k") != TEMPERATURE_K
        or settings.get("geometry_sha256") != _canonical_hash(geometry)
        or settings.get("production_sha256") != _canonical_hash(production)
        or settings.get("geometry_fingerprint")
        != frequency_settings_fingerprint(expected_geometry)
        or settings.get("production_fingerprint")
        != frequency_settings_fingerprint(expected_production)
    ):
        raise RuntimeError("calculation settings/method drifted")
    for role in STOICHIOMETRY:
        _validate_component_record(role, components[role], pair, settings)
    if receipt.get("optimizer_budget") != {
        role: {"calls": 1, "max_steps": MAX_STEPS, "retries": 0}
        for role in ("C", "V", "SiOH4")
    }:
        raise RuntimeError("calculation optimizer aggregate drifted")
    recomposed = _compose_thermochemistry(
        {role: components[role]["thermochemistry"] for role in STOICHIOMETRY}
    )
    claimed_thermochemistry = receipt.get("thermochemistry", {})
    if set(claimed_thermochemistry) != set(recomposed) or any(
        not _close(claimed_thermochemistry[key], value, tolerance=1.0e-8)
        for key, value in recomposed.items()
    ):
        raise RuntimeError("calculation thermochemistry arithmetic drifted")
    if recomposed["recomposition_error_kj_mol"] > 1.0e-8:
        raise RuntimeError("calculation independent recomposition exceeds tolerance")


def _add_done_job(
    store: Store,
    structure_id: int,
    kind: str,
    method: str,
    engine: str,
    results: dict[str, tuple[float, str]],
    *,
    detail: dict[str, Any],
) -> int:
    job = store.add_job(
        structure_id,
        kind,
        method,
        engine,
        detail=json.dumps(detail, separators=(",", ":"), sort_keys=True),
    )
    store.set_job_status(job, "done")
    for key, (value, units) in results.items():
        if not _finite(value):
            raise RuntimeError(f"refusing non-finite Store result {key}")
        store.add_result(job, key, float(value), units)
    return job


def write_calc005_store(
    path: str | Path,
    pair: Calc005Pair,
    components: dict[str, dict[str, Any]],
    calculation_receipt: dict[str, Any],
) -> None:
    """Atomically write a Store bound to an immutable calculation receipt."""

    path = Path(path)
    if calculation_receipt.get("verdict") != "passed-protocol-pilot":
        raise RuntimeError("Store only accepts passed-protocol-pilot calculations")
    if calculation_receipt.get("classification") != CLASSIFICATION:
        raise RuntimeError("Store refuses kinetic/emittable calculation payloads")
    if calculation_receipt.get("components") != {
        role: {
            key: value
            for key, value in components[role].items()
            if key not in {"cluster", "exact_xyz"}
        }
        for role in ("C", "V", "SiOH4")
    }:
        raise RuntimeError("Store components differ from immutable calculation receipt")
    validate_calculation_receipt(calculation_receipt, pair)
    if not all(
        bool(value) for value in pair.validate().values() if isinstance(value, bool)
    ):
        raise RuntimeError("Store refuses an invalid CALC-005 pair")

    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".sqlite", dir=path.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    temporary.unlink()
    try:
        with Store(temporary) as store:
            store.conn.executescript(_EXTENSION_SCHEMA)
            condensed_id = store.add_structure(
                pair.condensed.cluster.name,
                pair.condensed.cluster.formula,
                pair.condensed.cluster.to_xyz(),
                charge=pair.condensed.cluster.charge,
                spin=pair.condensed.cluster.spin,
            )
            structure_ids: dict[str, int] = {}
            job_ids: dict[str, tuple[int, int, int]] = {}
            engine = str(calculation_receipt["engine"])
            for role in ("C", "V", "SiOH4"):
                record = components[role]
                cluster = record["cluster"]
                structure_id = store.add_structure(
                    cluster.name,
                    cluster.formula,
                    record["exact_xyz"],
                    charge=cluster.charge,
                    spin=cluster.spin,
                )
                structure_ids[role] = structure_id
                optimizer = record["optimizer"]
                gradient = record["gradient"]
                frequency = record["frequency"]
                thermo = record["thermochemistry"]
                opt_job = _add_done_job(
                    store,
                    structure_id,
                    "opt",
                    R2SCAN3C_METHOD,
                    engine,
                    {
                        "optimizer_converged": (1.0, "boolean"),
                        "projected_gradient_rms": (
                            gradient["rms_hartree_per_bohr"],
                            "hartree/bohr",
                        ),
                        "projected_gradient_max": (
                            gradient["max_hartree_per_bohr"],
                            "hartree/bohr",
                        ),
                    },
                    detail={
                        "calls": optimizer["observed_calls"],
                        "max_steps": optimizer["observed_max_steps"],
                        "retries": optimizer["observed_retries"],
                    },
                )
                freq_job = _add_done_job(
                    store,
                    structure_id,
                    "freq",
                    R2SCAN3C_METHOD,
                    engine,
                    {
                        "electronic_hartree": (
                            frequency["electronic_hartree"],
                            "hartree",
                        ),
                        "frequency_passed": (1.0, "boolean"),
                        "imaginary_count": (frequency["imaginary_mode_count"], "count"),
                        "real_mode_count": (frequency["real_mode_count"], "count"),
                        "expected_mode_count": (
                            frequency["expected_mode_count"],
                            "count",
                        ),
                        "zpe": (thermo["zpe_kj_mol"], "kJ/mol"),
                        "thermal": (thermo["thermal_kj_mol"], "kJ/mol"),
                        "entropy": (thermo["entropy_kj_mol_k"], "kJ/mol/K"),
                    },
                    detail={
                        "hessian": frequency["hessian"],
                        "noise_floor_cm": IMAGINARY_NOISE_FLOOR_CM,
                    },
                )
                sp_job = _add_done_job(
                    store,
                    structure_id,
                    "sp",
                    PRODUCTION_METHOD,
                    engine,
                    {
                        "electronic_hartree": (
                            record["production_electronic_hartree"],
                            "hartree",
                        ),
                        "single_point_converged": (1.0, "boolean"),
                    },
                    detail={"solvent": "SMD(water)"},
                )
                job_ids[role] = (opt_job, freq_job, sp_job)

            calculation_hash = calculation_receipt_sha256(calculation_receipt)
            store.conn.execute(
                "INSERT INTO calc005_cycles VALUES (1,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (
                    "calc005-cycle-v2",
                    1,
                    0,
                    1,
                    205,
                    200,
                    402,
                    404,
                    TEMPERATURE_K,
                    "fixed-site unit activity; vibration-only",
                    "Si(OH)4(aq), 1 M from 1 bar",
                    "passed-protocol-pilot",
                    LABEL,
                    atom_map_sha256(pair),
                    calculation_hash,
                    condensed_id,
                ),
            )
            for role in ("C", "V", "SiOH4"):
                record = components[role]
                store.conn.execute(
                    "INSERT INTO calc005_components VALUES (1,?,?,?,?,?,?,?,?,?)",
                    (
                        role,
                        STOICHIOMETRY[role],
                        structure_ids[role],
                        *job_ids[role],
                        record["raw_endpoint"]["sha256"],
                        record["endpoint"]["sha256"],
                        record["checkpoint_sha256"],
                    ),
                )
            destinations = {
                (origin.kind, origin.node, origin.ordinal): ("V", index)
                for index, origin in enumerate(pair.vacancy_origins)
            }
            destinations.update(
                {
                    (origin.kind, origin.node, origin.ordinal): ("SiOH4", index)
                    for index, origin in enumerate(pair.silicic_acid_origins)
                }
            )
            for index, (element, origin) in enumerate(
                zip(pair.occupied.symbols, pair.occupied_origins, strict=True)
            ):
                product_role, product_index = destinations[
                    (origin.kind, origin.node, origin.ordinal)
                ]
                store.conn.execute(
                    "INSERT INTO calc005_construction VALUES (1,?,?,?,?,?,?,?,0)",
                    (
                        index,
                        element,
                        origin.kind,
                        _node_text(origin),
                        origin.ordinal,
                        product_role,
                        product_index,
                    ),
                )
            store.conn.commit()
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def _load_calculation_from_terminal(
    terminal: dict[str, Any], store_path: Path
) -> tuple[dict[str, Any], Path]:
    if (
        terminal.get("schema") != "calc005-terminal-receipt-v2"
        or terminal.get("verdict") != "passed-protocol-pilot"
        or terminal.get("classification") != CLASSIFICATION
        or terminal.get("canonical_value_exposed") is not True
        or terminal.get("forbidden_outputs_emitted") is not False
        or terminal.get("independent_verification_required") is not True
        or terminal.get("receipt_payload_sha256") != receipt_payload_sha256(terminal)
    ):
        raise RuntimeError("CALC-005 terminal identity/classification/hash drifted")
    artifacts = terminal.get("artifacts", {})
    if set(artifacts) != {"result", "calculation_receipt", "store"}:
        raise RuntimeError("CALC-005 terminal artifact manifest drifted")
    artifact_paths = {
        name: Path(str(artifact.get("path", ""))).resolve()
        for name, artifact in artifacts.items()
    }
    generation = artifact_paths["store"].parent
    if (
        {path.parent for path in artifact_paths.values()} != {generation}
        or generation.name != terminal.get("generation")
        or generation.parent.name != "generations"
        or artifact_paths["result"].name != "calc005-result.json"
        or artifact_paths["calculation_receipt"].name != "calculation-receipt.json"
        or artifact_paths["store"].name != "store.sqlite"
    ):
        raise RuntimeError("CALC-005 terminal does not point to one atomic generation")
    calculation_artifact = artifacts["calculation_receipt"]
    calculation_path = artifact_paths["calculation_receipt"]
    if not calculation_path.is_file() or calculation_artifact.get(
        "sha256"
    ) != _sha256_path(calculation_path):
        raise RuntimeError("CALC-005 immutable calculation receipt hash drifted")
    calculation = json.loads(calculation_path.read_text())
    if calculation_receipt_sha256(calculation) != calculation_artifact.get(
        "canonical_sha256"
    ):
        raise RuntimeError("CALC-005 calculation receipt canonical content drifted")
    store_artifact = artifacts["store"]
    if artifact_paths["store"] != store_path or store_artifact.get(
        "sha256"
    ) != _sha256_path(store_path):
        raise RuntimeError("CALC-005 closed Store terminal hash drifted")
    result_artifact = artifacts["result"]
    result_path = artifact_paths["result"]
    if not result_path.is_file() or result_artifact.get("sha256") != _sha256_path(
        result_path
    ):
        raise RuntimeError("CALC-005 result artifact hash drifted")
    result = json.loads(result_path.read_text())
    if (
        result.get("schema") != "calc005-result-v1"
        or result.get("verdict") != "passed-protocol-pilot"
        or result.get("classification") != CLASSIFICATION
        or result.get("temperature_k") != TEMPERATURE_K
        or result.get("thermochemistry") != calculation.get("thermochemistry")
        or result.get("stoichiometry") != STOICHIOMETRY
        or result.get("water_thermochemical_coefficient") != 0
        or calculation.get("result_sha256") != result_artifact.get("sha256")
    ):
        raise RuntimeError("CALC-005 result/calculation content drifted")
    return calculation, calculation_path


def _validate_store_graph(
    connection: sqlite3.Connection,
    calculation: dict[str, Any],
    pair: Calc005Pair,
) -> dict[str, Any]:
    integrity = connection.execute("PRAGMA integrity_check").fetchone()[0]
    if integrity != "ok" or connection.execute("PRAGMA foreign_key_check").fetchall():
        raise RuntimeError("CALC-005 Store integrity/foreign-key check failed")
    tables = {
        row[0]
        for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type='table' "
            "AND name NOT LIKE 'sqlite_%'"
        )
    }
    expected_tables = {
        "structures",
        "jobs",
        "results",
        "calc005_cycles",
        "calc005_components",
        "calc005_construction",
    }
    if tables != expected_tables:
        raise RuntimeError("CALC-005 Store table set drifted")
    counts = {
        table: connection.execute(f"SELECT count(*) FROM {table}").fetchone()[0]
        for table in ("structures", "jobs", "results")
    }
    if counts != {"structures": 4, "jobs": 9, "results": 39}:
        raise RuntimeError(f"CALC-005 base Store cardinality drifted: {counts}")
    cycles = connection.execute("SELECT * FROM calc005_cycles").fetchall()
    if len(cycles) != 1:
        raise RuntimeError("CALC-005 Store requires exactly one cycle")
    cycle = cycles[0]
    expected_cycle = {
        "schema": "calc005-cycle-v2",
        "x": 1,
        "y": 0,
        "environment_index": 1,
        "center_state_from": 205,
        "center_state_to": 200,
        "osa_state_from": 402,
        "osa_state_to": 404,
        "temperature_k": TEMPERATURE_K,
        "surface_standard_state": "fixed-site unit activity; vibration-only",
        "solute_standard_state": "Si(OH)4(aq), 1 M from 1 bar",
        "status": "passed-protocol-pilot",
        "label": LABEL,
        "atom_map_sha256": atom_map_sha256(pair),
        "calculation_receipt_sha256": calculation_receipt_sha256(calculation),
    }
    if any(cycle[key] != value for key, value in expected_cycle.items()):
        raise RuntimeError("CALC-005 cycle identity/calculation-hash/label drifted")

    condensed = connection.execute(
        "SELECT formula,charge,spin,xyz,geometry_hash FROM structures WHERE id=?",
        (cycle["condensed_structure_id"],),
    ).fetchone()
    if condensed is None or (
        condensed["formula"] != pair.condensed.cluster.formula
        or condensed["charge"] != pair.condensed.cluster.charge
        or condensed["spin"] != pair.condensed.cluster.spin
        or geometry_hash(condensed["xyz"]) != condensed["geometry_hash"]
        or condensed["geometry_hash"] != geometry_hash(pair.condensed.cluster.to_xyz())
    ):
        raise RuntimeError("CALC-005 condensed structure edge drifted")

    components = connection.execute(
        "SELECT c.*,s.formula,s.charge,s.spin,s.xyz,s.geometry_hash "
        "FROM calc005_components c JOIN structures s ON s.id=c.structure_id "
        "ORDER BY c.role"
    ).fetchall()
    if [row["role"] for row in components] != ["C", "SiOH4", "V"]:
        raise RuntimeError("CALC-005 component roles are missing or duplicated")
    if {row["role"]: row["coefficient"] for row in components} != STOICHIOMETRY:
        raise RuntimeError("CALC-005 component coefficients drifted")
    templates = {"C": pair.occupied, "V": pair.vacancy, "SiOH4": pair.silicic_acid}
    expected_job_ids: dict[int, tuple[str, str]] = {}
    for row in components:
        role = row["role"]
        claimed = calculation["components"][role]
        if (
            row["formula"] != templates[role].formula
            or row["charge"] != templates[role].charge
            or row["spin"] != templates[role].spin
            or geometry_hash(row["xyz"]) != row["geometry_hash"]
            or row["geometry_hash"] != claimed["endpoint"]["geometry_hash"]
            or row["endpoint_sha256"] != claimed["endpoint"]["sha256"]
            or row["raw_endpoint_sha256"] != claimed["raw_endpoint"]["sha256"]
            or row["checkpoint_sha256"] != claimed["checkpoint_sha256"]
        ):
            raise RuntimeError(
                f"CALC-005 {role} structure/artifact/checkpoint edge drifted"
            )
        expected_job_ids[row["opt_job_id"]] = (role, "opt")
        expected_job_ids[row["freq_job_id"]] = (role, "freq")
        expected_job_ids[row["sp_job_id"]] = (role, "sp")

    jobs = connection.execute(
        "SELECT j.*,c.role FROM jobs j JOIN calc005_components c ON "
        "j.id IN (c.opt_job_id,c.freq_job_id,c.sp_job_id) ORDER BY j.id"
    ).fetchall()
    if len(jobs) != 9 or len(expected_job_ids) != 9:
        raise RuntimeError("CALC-005 Store requires exactly nine distinct jobs")
    structure_by_role = {row["role"]: row["structure_id"] for row in components}
    method_by_kind = {
        "opt": R2SCAN3C_METHOD,
        "freq": R2SCAN3C_METHOD,
        "sp": PRODUCTION_METHOD,
    }
    detail_by_kind = {
        "opt": {"calls": 1, "max_steps": MAX_STEPS, "retries": 0},
        "freq": None,
        "sp": {"solvent": "SMD(water)"},
    }
    for job in jobs:
        role, kind = expected_job_ids.get(job["id"], (None, None))
        detail = json.loads(job["detail"])
        expected_detail = detail_by_kind[kind]
        if kind == "freq":
            expected_detail = {
                "hessian": "full" if role == "SiOH4" else "PHVA",
                "noise_floor_cm": IMAGINARY_NOISE_FLOOR_CM,
            }
        if (
            job["role"] != role
            or job["kind"] != kind
            or job["status"] != "done"
            or job["method"] != method_by_kind[kind]
            or job["engine"] != calculation["engine"]
            or job["structure_id"] != structure_by_role[role]
            or detail != expected_detail
        ):
            raise RuntimeError("CALC-005 job status/method/budget edge drifted")

    results = connection.execute(
        "SELECT r.job_id,r.key,r.value,r.units,j.kind FROM results r "
        "JOIN jobs j ON j.id=r.job_id"
    ).fetchall()
    required_units = {
        "optimizer_converged": "boolean",
        "projected_gradient_rms": "hartree/bohr",
        "projected_gradient_max": "hartree/bohr",
        "electronic_hartree": "hartree",
        "frequency_passed": "boolean",
        "imaginary_count": "count",
        "real_mode_count": "count",
        "expected_mode_count": "count",
        "zpe": "kJ/mol",
        "thermal": "kJ/mol",
        "entropy": "kJ/mol/K",
        "single_point_converged": "boolean",
    }
    if any(
        row["key"] not in required_units
        or row["units"] != required_units[row["key"]]
        or not _finite(row["value"])
        for row in results
    ):
        raise RuntimeError("CALC-005 result key/unit/value drifted")
    expected_keys = {
        "opt": {
            "optimizer_converged",
            "projected_gradient_rms",
            "projected_gradient_max",
        },
        "freq": {
            "electronic_hartree",
            "frequency_passed",
            "imaginary_count",
            "real_mode_count",
            "expected_mode_count",
            "zpe",
            "thermal",
            "entropy",
        },
        "sp": {"electronic_hartree", "single_point_converged"},
    }
    for job in jobs:
        observed = {
            row["key"]: row["value"] for row in results if row["job_id"] == job["id"]
        }
        role = job["role"]
        claimed = calculation["components"][role]
        if set(observed) != expected_keys[job["kind"]]:
            raise RuntimeError("CALC-005 per-job result key graph drifted")
        if job["kind"] == "opt":
            expected = {
                "optimizer_converged": 1.0,
                "projected_gradient_rms": claimed["gradient"]["rms_hartree_per_bohr"],
                "projected_gradient_max": claimed["gradient"]["max_hartree_per_bohr"],
            }
        elif job["kind"] == "freq":
            expected = {
                "electronic_hartree": claimed["frequency"]["electronic_hartree"],
                "frequency_passed": 1.0,
                "imaginary_count": claimed["frequency"]["imaginary_mode_count"],
                "real_mode_count": claimed["frequency"]["real_mode_count"],
                "expected_mode_count": claimed["frequency"]["expected_mode_count"],
                "zpe": claimed["thermochemistry"]["zpe_kj_mol"],
                "thermal": claimed["thermochemistry"]["thermal_kj_mol"],
                "entropy": claimed["thermochemistry"]["entropy_kj_mol_k"],
            }
        else:
            expected = {
                "electronic_hartree": claimed["production_electronic_hartree"],
                "single_point_converged": 1.0,
            }
        if any(
            not _close(observed[key], value, tolerance=1.0e-8)
            for key, value in expected.items()
        ):
            raise RuntimeError("CALC-005 per-job result values drifted")

    construction = connection.execute(
        "SELECT * FROM calc005_construction ORDER BY c_atom_index"
    ).fetchall()
    destinations = {
        (origin.kind, origin.node, origin.ordinal): ("V", index)
        for index, origin in enumerate(pair.vacancy_origins)
    }
    destinations.update(
        {
            (origin.kind, origin.node, origin.ordinal): ("SiOH4", index)
            for index, origin in enumerate(pair.silicic_acid_origins)
        }
    )
    expected_construction = []
    for index, (element, origin) in enumerate(
        zip(pair.occupied.symbols, pair.occupied_origins, strict=True)
    ):
        role, product_index = destinations[(origin.kind, origin.node, origin.ordinal)]
        expected_construction.append(
            (
                1,
                index,
                element,
                origin.kind,
                _node_text(origin),
                origin.ordinal,
                role,
                product_index,
                0,
            )
        )
    if [tuple(row) for row in construction] != expected_construction:
        raise RuntimeError("CALC-005 construction atom edges drifted")
    water_rows = [
        row
        for row in construction
        if row["origin_kind"].startswith("hydrolysis-water-")
    ]
    if len(water_rows) != 3 or any(row["thermochemical"] != 0 for row in construction):
        raise RuntimeError("CALC-005 construction-water provenance drifted")
    return {
        "status": cycle["status"],
        "roles": [row["role"] for row in components],
        "coefficients": {row["role"]: row["coefficient"] for row in components},
        "job_count": len(jobs),
        "construction_water_atoms": len(water_rows),
        "production_energy_scale_kj_mol": HARTREE_TO_KJ,
        "calculation_receipt_sha256": cycle["calculation_receipt_sha256"],
    }


def validate_calc005_store(
    path: str | Path,
    terminal: dict[str, Any],
    pair: Calc005Pair,
) -> dict[str, Any]:
    """Validate terminal, immutable calculation content, and Store independently."""

    path = Path(path).resolve()
    if not path.is_file():
        raise RuntimeError("CALC-005 Store is missing")
    try:
        calculation, _ = _load_calculation_from_terminal(terminal, path)
        validate_calculation_receipt(calculation, pair)
        restoration = terminal.get("restoration")
        if calculation.get("backend_kind") == "production":
            lease_path = str(
                calculation.get("execution_envelope", {})
                .get("qi2_lease", {})
                .get("path", "")
            )
            if (
                not isinstance(restoration, dict)
                or restoration.get("status") != "verified-restored"
                or not restoration.get("checked_at")
                or restoration.get("bootstrap_session_closed") is not True
                or restoration.get("qi2_lease_released") is not True
                or restoration.get("lease_path") != lease_path
            ):
                raise RuntimeError("CALC-005 successful restoration receipt drifted")
        elif restoration != {"status": "not-applicable-in-process-test"}:
            raise RuntimeError("CALC-005 in-process restoration receipt drifted")
        uri = f"{path.as_uri()}?mode=ro&immutable=1"
        with sqlite3.connect(uri, uri=True) as connection:
            connection.row_factory = sqlite3.Row
            return _validate_store_graph(connection, calculation, pair)
    except (
        sqlite3.DatabaseError,
        OSError,
        json.JSONDecodeError,
        KeyError,
        TypeError,
        ValueError,
    ) as exc:
        raise RuntimeError(f"CALC-005 Store validation failed: {exc}") from exc


def validate_staged_calc005_store(
    path: str | Path,
    calculation: dict[str, Any],
    pair: Calc005Pair,
) -> dict[str, Any]:
    """Validate a closed staged Store before its generation is promoted."""

    path = Path(path).resolve()
    validate_calculation_receipt(calculation, pair)
    uri = f"{path.as_uri()}?mode=ro&immutable=1"
    with sqlite3.connect(uri, uri=True) as connection:
        connection.row_factory = sqlite3.Row
        return _validate_store_graph(connection, calculation, pair)
