#!/usr/bin/env python3
"""Independent optimizer-free verifier for the CALC-005 Si pilot."""

from __future__ import annotations

import argparse
import fcntl
import hashlib
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import time
import tomllib
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import asdict, replace
from datetime import datetime
from pathlib import Path
from typing import Any, Protocol

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

GPU_TTL_HOURS = 13.0
EXPECTED_MEMORY_MAX_BYTES = 32 * 1024**3
EXPECTED_MEMORY_SWAP_MAX_BYTES = 4 * 1024**3
EXPECTED_CPU_QUOTA_PERCENT = 1600
EXPECTED_NICE = 10
EXPECTED_THREADS = 16
EXPECTED_RUNTIME_MAX_SECONDS = 43200
EXPECTED_GPU_OWNER = "calc005_si_attachment_verify"

import numpy as np  # noqa: E402

from quarry.calc005 import (  # noqa: E402
    R_KJ_MOL_K,
    STOICHIOMETRY,
    TEMPERATURE_K,
    Calc005Pair,
    build_calc005_pair,
)
from quarry.calc005_store import validate_calc005_store  # noqa: E402
from quarry.clusters import Cluster  # noqa: E402
from quarry.pipeline import (  # noqa: E402
    HARTREE_TO_KJ,
    DftSettings,
    FrequencyResult,
    frequency_geometry_fingerprint,
    frequency_settings_fingerprint,
)
from quarry.rates import (  # noqa: E402
    surface_thermo_from_frequencies,
    thermo_from_frequencies,
)
from quarry.store import geometry_hash  # noqa: E402

GRADIENT_RMS_MAX = 3.0e-4
GRADIENT_MAX_MAX = 4.5e-4
IMAGINARY_NOISE_FLOOR_CM = 30.0
MIN_PAIR_DISTANCE_A = 0.75
MAX_OWNER_DISTANCE_A = 1.25
MIN_OWNER_MARGIN_A = 0.15
STANDARD_STATE_CORRECTION_KJ_MOL = 7.958500693927389
CLASSIFICATION = "thermodynamic/non-kinetic/non-emittable"
R2SCAN3C_METHOD = "r2scan-3c/def2-mtzvpp/d4/gcp"
PRODUCTION_METHOD = "wb97m-v/def2-tzvpd/smd(water)"


class CalculatorBackend(Protocol):
    production_backend: bool

    def gradient(
        self, role: str, cluster: Cluster, settings: DftSettings
    ) -> np.ndarray: ...

    def frequencies(
        self, role: str, cluster: Cluster, settings: DftSettings
    ) -> FrequencyResult: ...

    def energy(self, role: str, cluster: Cluster, settings: DftSettings) -> float: ...


class PipelineBackend:
    """Verifier calculator boundary; deliberately has no optimizer method."""

    production_backend = True

    def gradient(self, role, cluster, settings):
        from quarry.pipeline import gradient

        return gradient(cluster, settings)

    def frequencies(self, role, cluster, settings):
        from quarry.pipeline import frequencies

        return frequencies(cluster, settings)

    def energy(self, role, cluster, settings):
        from quarry.pipeline import energy

        return energy(cluster, settings)


def _is_production_backend(backend: CalculatorBackend) -> bool:
    pipeline_type = PipelineBackend
    return backend.production_backend is True or (
        isinstance(pipeline_type, type) and isinstance(backend, pipeline_type)
    )


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
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        temporary.replace(path)
        directory = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        temporary.unlink(missing_ok=True)


def _canonical_hash(payload: dict[str, Any]) -> str:
    return hashlib.sha256(
        json.dumps(
            payload, separators=(",", ":"), sort_keys=True, allow_nan=False
        ).encode()
    ).hexdigest()


def _terminal_hash(receipt: dict[str, Any]) -> str:
    payload = json.loads(json.dumps(receipt, allow_nan=False))
    payload.pop("receipt_payload_sha256", None)
    return _canonical_hash(payload)


def _checkpoint_hash(record: dict[str, Any]) -> str:
    return _canonical_hash(
        {
            key: value
            for key, value in record.items()
            if key not in {"checkpoint_sha256", "accepted_at", "cluster", "exact_xyz"}
        }
    )


def _settings(use_gpu: bool) -> tuple[DftSettings, DftSettings]:
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


def _read_xyz(path: Path, template: Cluster) -> Cluster:
    lines = path.read_text().splitlines()
    if len(lines) != len(template.symbols) + 2:
        raise RuntimeError(f"malformed XYZ: {path}")
    try:
        count = int(lines[0])
        rows = [line.split() for line in lines[2:]]
        coords = np.asarray([[float(value) for value in row[1:]] for row in rows])
    except (IndexError, ValueError) as exc:
        raise RuntimeError(f"malformed XYZ: {path}") from exc
    if (
        count != len(template.symbols)
        or any(len(row) != 4 for row in rows)
        or [row[0] for row in rows] != template.symbols
        or coords.shape != template.coords.shape
        or not np.all(np.isfinite(coords))
    ):
        raise RuntimeError(f"XYZ identity/order/coordinate mismatch: {path}")
    return replace(template, name=lines[1], coords=coords)


def _owners(cluster: Cluster) -> dict[int, int]:
    oxygens = [index for index, symbol in enumerate(cluster.symbols) if symbol == "O"]
    owners = {}
    for hydrogen, symbol in enumerate(cluster.symbols):
        if symbol != "H":
            continue
        distances = sorted(
            (
                float(
                    np.linalg.norm(cluster.coords[hydrogen] - cluster.coords[oxygen])
                ),
                oxygen,
            )
            for oxygen in oxygens
        )
        if (
            len(distances) < 2
            or distances[0][0] > MAX_OWNER_DISTANCE_A
            or distances[1][0] - distances[0][0] < MIN_OWNER_MARGIN_A
        ):
            raise RuntimeError(f"H{hydrogen} fails owner distance/margin")
        owners[hydrogen] = distances[0][1]
    return owners


def _heavy_edges(cluster: Cluster) -> set[tuple[int, int]]:
    return {
        (cation, oxygen)
        for cation, symbol in enumerate(cluster.symbols)
        if symbol in {"Al", "Si"}
        for oxygen, oxygen_symbol in enumerate(cluster.symbols)
        if oxygen_symbol == "O"
        and float(np.linalg.norm(cluster.coords[cation] - cluster.coords[oxygen]))
        <= 2.35
    }


def _structure_gate(
    role: str,
    endpoint: Cluster,
    reference: Cluster,
    origins: tuple[Any, ...],
) -> dict[str, Any]:
    if (
        endpoint.symbols != reference.symbols
        or endpoint.charge != reference.charge
        or endpoint.spin != reference.spin
        or endpoint.frozen_indices != reference.frozen_indices
        or len(origins) != len(endpoint.symbols)
    ):
        raise RuntimeError(f"{role} identity/charge/spin/frozen/provenance drifted")
    frozen = sorted(reference.frozen_indices)
    if frozen and not np.array_equal(endpoint.coords[frozen], reference.coords[frozen]):
        raise RuntimeError(f"{role} frozen coordinates changed")
    delta = endpoint.coords[:, None, :] - endpoint.coords[None, :, :]
    distances = np.linalg.norm(delta, axis=2)
    distances[np.diag_indices_from(distances)] = np.inf
    minimum = float(np.min(distances))
    if minimum < MIN_PAIR_DISTANCE_A:
        raise RuntimeError(f"{role} collision at {minimum} A")
    owners = _owners(endpoint)
    if owners != _owners(reference):
        raise RuntimeError(f"{role} proton owner changed")
    edges = _heavy_edges(endpoint)
    if edges != _heavy_edges(reference):
        raise RuntimeError(f"{role} provenance-bound heavy connectivity changed")
    center = [
        index
        for index, origin in enumerate(origins)
        if origin.kind == "deck" and origin.node == (4, (0, 0, 0))
    ]
    center_bonds = {
        oxygen for cation, oxygen in edges if center and cation == center[0]
    }
    if role in {"C", "SiOH4"} and (len(center) != 1 or len(center_bonds) != 4):
        raise RuntimeError(f"{role} lacks exact four center Si-O bonds")
    if role == "V" and center:
        raise RuntimeError("V retains center Si")
    return {
        "minimum_pair_distance_a": minimum,
        "owner_map": owners,
        "heavy_edges": edges,
        "center_si_o_bond_count": len(center_bonds),
    }


def _gradient_gate(
    role: str, cluster: Cluster, gradient: np.ndarray
) -> dict[str, float]:
    projected = np.asarray(gradient, dtype=float).copy()
    if projected.shape != cluster.coords.shape or not np.all(np.isfinite(projected)):
        raise RuntimeError(f"{role} gradient is malformed/non-finite")
    projected[cluster.frozen_indices] = 0.0
    free = sorted(set(range(len(cluster.symbols))) - set(cluster.frozen_indices))
    values = projected[free]
    rms = float(np.sqrt(np.mean(values**2))) if values.size else 0.0
    maximum = float(np.max(np.abs(values))) if values.size else 0.0
    if rms > GRADIENT_RMS_MAX or maximum > GRADIENT_MAX_MAX:
        raise RuntimeError(f"{role} gradient exceeds RMS/max thresholds")
    return {"rms_hartree_per_bohr": rms, "max_hartree_per_bohr": maximum}


def _expected_modes(role: str, cluster: Cluster) -> int:
    if role == "SiOH4":
        return 3 * len(cluster.symbols) - 6
    return 3 * (len(cluster.symbols) - len(cluster.frozen_indices))


def _frequency_gate(
    role: str,
    result: FrequencyResult,
    cluster: Cluster,
    settings: DftSettings,
) -> dict[str, Any]:
    real = np.asarray(result.frequencies_cm, dtype=float)
    imaginary = np.asarray(result.imaginary_cm, dtype=float)
    rotational = result.rotational_temperatures_k
    rotations = np.asarray(rotational if rotational is not None else [], dtype=float)
    expected = _expected_modes(role, cluster)
    if (
        real.ndim != 1
        or imaginary.ndim != 1
        or len(real) + len(imaginary) != expected
        or not np.all(np.isfinite(real))
        or not np.all(np.isfinite(imaginary))
        or np.any(real <= 0.0)
        or np.any(imaginary < 0.0)
        or np.any(imaginary > IMAGINARY_NOISE_FLOOR_CM)
        or not math.isfinite(float(result.electronic_hartree))
        or not math.isfinite(float(result.molar_mass_kg))
        or float(result.molar_mass_kg) <= 0.0
        or result.geometry_fingerprint != frequency_geometry_fingerprint(cluster)
        or result.settings_fingerprint != frequency_settings_fingerprint(settings)
    ):
        raise RuntimeError(f"{role} frequency gate/count/binding failed")
    if role == "SiOH4":
        if (
            result.linear
            or rotational is None
            or len(rotational) != 3
            or np.any(rotations <= 0)
        ):
            raise RuntimeError("SiOH4 nonlinear rotational evidence failed")
    elif rotational is not None or not cluster.frozen_indices:
        raise RuntimeError(f"{role} PHVA frozen-shell/rotation evidence failed")
    return {
        "frequencies_cm": real.tolist(),
        "imaginary_cm": imaginary.tolist(),
        "real_mode_count": len(real),
        "imaginary_mode_count": len(imaginary),
        "expected_mode_count": expected,
        "observed_mode_count": len(real) + len(imaginary),
        "electronic_hartree": float(result.electronic_hartree),
        "molar_mass_kg": float(result.molar_mass_kg),
        "rotational_temperatures_k": list(rotational)
        if rotational is not None
        else None,
        "linear": bool(result.linear),
    }


def _derive_thermal_terms(
    role: str,
    frequency: FrequencyResult,
    production_electronic_hartree: float,
) -> dict[str, float]:
    electronic = production_electronic_hartree * HARTREE_TO_KJ
    if role == "SiOH4":
        thermo = thermo_from_frequencies(
            electronic,
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
            electronic, frequency.frequencies_cm, TEMPERATURE_K
        )
    values = {
        "electronic_kj_mol": electronic,
        "zpe_kj_mol": thermo.zpe_kj,
        "thermal_kj_mol": thermo.thermal_kj,
        "entropy_kj_mol_k": thermo.entropy_kj_per_k,
    }
    if not all(math.isfinite(value) for value in values.values()):
        raise RuntimeError(f"{role} thermochemistry is non-finite")
    return values


def _compose(components: dict[str, dict[str, float]]) -> dict[str, float]:
    if set(components) != set(STOICHIOMETRY):
        raise RuntimeError("CALC-005 requires exactly C, V, SiOH4")

    def delta(key: str) -> float:
        return math.fsum(
            STOICHIOMETRY[role] * components[role][key] for role in ("C", "V", "SiOH4")
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
                    components[role]["electronic_kj_mol"],
                    components[role]["zpe_kj_mol"],
                    components[role]["thermal_kj_mol"],
                    -TEMPERATURE_K * components[role]["entropy_kj_mol_k"],
                )
            )
            for role in ("C", "V", "SiOH4")
        )
        + STANDARD_STATE_CORRECTION_KJ_MOL
    )
    error = abs(total - independent)
    if error > 1.0e-8:
        raise RuntimeError("independent CALC-005 recomposition exceeds tolerance")
    return {
        "electronic_kj_mol": electronic,
        "zpe_kj_mol": zpe,
        "thermal_kj_mol": thermal,
        "minus_t_delta_s_kj_mol": minus_t_delta_s,
        "standard_state_kj_mol": STANDARD_STATE_CORRECTION_KJ_MOL,
        "s_10_kj_mol": total,
        "independent_s_10_kj_mol": independent,
        "recomposition_error_kj_mol": error,
    }


def _semantic_equal(left: Any, right: Any, *, tolerance: float = 1.0e-8) -> bool:
    if isinstance(left, bool) or isinstance(right, bool):
        return left is right
    if isinstance(left, (int, float)) and isinstance(right, (int, float)):
        return (
            math.isfinite(float(left))
            and math.isfinite(float(right))
            and math.isclose(
                float(left), float(right), rel_tol=1.0e-12, abs_tol=tolerance
            )
        )
    if isinstance(left, dict) and isinstance(right, dict):
        return set(left) == set(right) and all(
            _semantic_equal(left[key], right[key], tolerance=tolerance) for key in left
        )
    if isinstance(left, list) and isinstance(right, list):
        return len(left) == len(right) and all(
            _semantic_equal(a, b, tolerance=tolerance)
            for a, b in zip(left, right, strict=True)
        )
    return left == right


def _selector(value: dict[str, Any]) -> dict[str, Any]:
    return {
        "distance": value.get("distance", 1),
        "kind": value.get("kind"),
        "label": value.get("label"),
        "exclude_label": value.get("exclude_label"),
        "state": value.get("state"),
        "frozen": value.get("frozen"),
        "min": value.get("min"),
        "max": value.get("max"),
    }


def _effect(value: dict[str, Any]) -> dict[str, Any]:
    return {
        "target": value.get("target"),
        "select": _selector(value["select"]) if "select" in value else None,
        "select_mode": value.get("select_mode"),
        "set": value.get("set"),
        "shift": value.get("shift"),
        "map": value.get("map"),
        "missing": value.get("missing"),
    }


def _reaction(document: dict[str, Any], name: str) -> dict[str, Any]:
    rows = [row for row in document.get("reactions", []) if row.get("name") == name]
    if len(rows) != 1:
        raise RuntimeError(f"Petra deck requires exactly one {name} reaction")
    return rows[0]


def _cpp_function(source: str, signature: str) -> str:
    start = source.find(signature)
    if start < 0:
        raise RuntimeError(f"legacy source lacks {signature}")
    opening = source.find("{", start + len(signature))
    if opening < 0:
        raise RuntimeError(f"legacy source has malformed {signature}")
    depth = 0
    for index in range(opening, len(source)):
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
            if depth == 0:
                return source[opening + 1 : index]
    raise RuntimeError(f"legacy source has unterminated {signature}")


def _legacy_case_transitions(body: str) -> dict[int, int]:
    pairs = re.findall(
        r"case\s+(\d+)\s*:\s*"
        r"this->lattice->sites\[nbr\]\.state\s*=\s*(\d+)\s*;",
        body,
    )
    transitions = {int(source): int(target) for source, target in pairs}
    if len(transitions) != len(pairs):
        raise RuntimeError("legacy Si transitions contain duplicate switch cases")
    return transitions


def _validate_legacy_si_boundary(actions_source: str, envrn_source: str) -> None:
    adsorb = _cpp_function(actions_source, "bool Actions::AdsorbSi(int site)")
    desorb = _cpp_function(actions_source, "void Actions::DesorbSi(int site)")
    adsorb_center_updates = re.findall(r"sites\[site\]\.state\s*=\s*([^;]+);", adsorb)
    desorb_center_updates = re.findall(r"sites\[site\]\.state\s*=\s*([^;]+);", desorb)
    if (
        adsorb_center_updates != ["205", "100 + WRONG"]
        or "if (ISSI(this->lattice->sites[site]))" not in adsorb
        or "for (i = 0; i < 4; i++)" not in adsorb
        or [int(value) for value in re.findall(r"case\s+(\d+)\s*:", adsorb)]
        != [300, 303, 404, 405, 409, 400]
        or _legacy_case_transitions(adsorb)
        != {300: 303, 303: 302, 404: 402, 405: 403, 409: 407, 400: 408}
        or desorb_center_updates != ["200", "100"]
        or "if (ISSI(this->lattice->sites[site]))" not in desorb
        or "for (i = 0; i < 4; i++)" not in desorb
        or [int(value) for value in re.findall(r"case\s+(\d+)\s*:", desorb)]
        != [303, 302, 402, 403, 407, 408]
        or _legacy_case_transitions(desorb)
        != {303: 300, 302: 303, 402: 404, 403: 405, 407: 409, 408: 400}
    ):
        raise RuntimeError("legacy Si transitions drifted")

    check200 = _cpp_function(envrn_source, "int Environment::Check200(int site)")
    compact = re.sub(r"\s+", " ", check200).strip()
    required = (
        "y = 0; x = 1;",
        "for (i = 0; i < 4; i++)",
        "if (this->lattice->sites[nbr].state == EDGE)",
        "case 408: x = 0; break;",
        "case 302: y++; break;",
        "return (x + y);",
    )
    assignments = re.findall(r"\b([xy])\s*(=|\+\+|--|\+=|-=)\s*(\d*)\s*;", compact)
    cases = [int(value) for value in re.findall(r"case\s+(\d+)\s*:", compact)]
    returns = re.findall(r"return\s+([^;]+);", compact)
    if (
        any(fragment not in compact for fragment in required)
        or assignments
        != [("y", "=", "0"), ("x", "=", "1"), ("x", "=", "0"), ("y", "++", "")]
        or cases != [408, 302]
        or returns != ["-1", "(x + y)"]
    ):
        raise RuntimeError("legacy Check200 ladder drifted")


def validate_petra_boundary(repo_root: Path) -> dict[str, Any]:
    """Bind exact source selectors/effects/modifiers and prove non-publication."""

    repo_root = repo_root.resolve()
    deck_path = repo_root / "petra/examples/kaolinite.toml"
    calculations_path = repo_root / "qm/CALCULATIONS.md"
    actions_path = repo_root / "legacy/cpp-model/actions.cpp"
    envrn_path = repo_root / "legacy/cpp-model/envrn.cpp"
    actions_source = actions_path.read_text()
    envrn_source = envrn_path.read_text()
    _validate_legacy_si_boundary(actions_source, envrn_source)
    document = tomllib.loads(deck_path.read_text())
    adsorb = _reaction(document, "adsorb-si")
    desorb = _reaction(document, "desorb-si")
    expected_adsorb = {
        "center": {"kind": "Si", "state": ["empty"]},
        "guards": [
            {
                "distance": 1,
                "kind": None,
                "label": None,
                "exclude_label": None,
                "state": ["@o_occ"],
                "frozen": False,
                "min": 1,
                "max": None,
            }
        ],
        "consumes": ["Si"],
        "produces": [],
        "rate": {"arrhenius": {"prefactor": 100.0, "ea": 6.4}},
        "modifiers": [],
        "effects": [
            {
                "target": "center",
                "select": None,
                "select_mode": None,
                "set": "oh4",
                "shift": None,
                "map": None,
                "missing": None,
            },
            {
                "target": "neighbors",
                "select": {
                    "distance": 1,
                    "kind": "Oss",
                    "label": None,
                    "exclude_label": None,
                    "state": ["Oss.empty", "Oss.si1"],
                    "frozen": False,
                    "min": None,
                    "max": None,
                },
                "select_mode": None,
                "set": None,
                "shift": None,
                "map": {"empty": "si1", "si1": "hy"},
                "missing": None,
            },
            {
                "target": "neighbors",
                "select": {
                    "distance": 1,
                    "kind": "Osa",
                    "label": None,
                    "exclude_label": None,
                    "state": ["Osa.empty", "Osa.al1", "Osa.albr", "Osa.alhy"],
                    "frozen": False,
                    "min": None,
                    "max": None,
                },
                "select_mode": None,
                "set": None,
                "shift": None,
                "map": {"empty": "si1", "al1": "sialh", "albr": "sih", "alhy": "full"},
                "missing": None,
            },
        ],
    }
    expected_desorb = {
        "center": {"kind": "Si", "state": ["oh4"]},
        "guards": [],
        "consumes": [],
        "produces": ["Si"],
        "rate": {"arrhenius": {"prefactor": 1998.0, "ea": 6.0}},
        "modifiers": [
            {
                "select": {
                    "distance": 1,
                    "kind": None,
                    "label": None,
                    "exclude_label": None,
                    "state": ["Oss.hy"],
                    "frozen": False,
                    "min": None,
                    "max": None,
                },
                "by_count": {"dea": [0.0, 6.0, 12.0, 18.0, 24.0]},
                "when": None,
                "per_match": None,
            },
            {
                "select": {
                    "distance": 1,
                    "kind": None,
                    "label": None,
                    "exclude_label": None,
                    "state": ["Osa.si1"],
                    "frozen": False,
                    "min": None,
                    "max": None,
                },
                "by_count": None,
                "when": {"min": 1, "dea": -6.0},
                "per_match": None,
            },
        ],
        "effects": [
            {
                "target": "center",
                "select": None,
                "select_mode": None,
                "set": "empty",
                "shift": None,
                "map": None,
                "missing": None,
            },
            {
                "target": "neighbors",
                "select": {
                    "distance": 1,
                    "kind": "Oss",
                    "label": None,
                    "exclude_label": None,
                    "state": ["*"],
                    "frozen": False,
                    "min": None,
                    "max": None,
                },
                "select_mode": None,
                "set": None,
                "shift": None,
                "map": {"si1": "empty", "hy": "si1"},
                "missing": "skip",
            },
            {
                "target": "neighbors",
                "select": {
                    "distance": 1,
                    "kind": "Osa",
                    "label": None,
                    "exclude_label": None,
                    "state": ["*"],
                    "frozen": False,
                    "min": None,
                    "max": None,
                },
                "select_mode": None,
                "set": None,
                "shift": None,
                "map": {"sih": "albr", "full": "alhy", "sialh": "al1", "si1": "empty"},
                "missing": "skip",
            },
        ],
    }

    def normalized(reaction: dict[str, Any]) -> dict[str, Any]:
        return {
            "center": reaction.get("center"),
            "guards": [_selector(value) for value in reaction.get("guards", [])],
            "consumes": reaction.get("consumes", []),
            "produces": reaction.get("produces", []),
            "rate": reaction.get("rate"),
            "modifiers": [
                {
                    "select": _selector(value.get("select", {})),
                    "by_count": value.get("by_count"),
                    "when": value.get("when"),
                    "per_match": value.get("per_match"),
                }
                for value in reaction.get("modifiers", [])
            ],
            "effects": [_effect(value) for value in reaction.get("effects", [])],
        }

    if normalized(adsorb) != expected_adsorb:
        raise RuntimeError("Petra adsorb-si exact selector/effect semantics drifted")
    if normalized(desorb) != expected_desorb:
        raise RuntimeError(
            "Petra desorb-si exact selector/modifier/effect semantics drifted"
        )
    calc_lines = [
        line
        for line in calculations_path.read_text().splitlines()
        if "| CALC-005 |" in line
    ]
    if len(calc_lines) != 1:
        raise RuntimeError("CALCULATIONS.md lacks exactly one CALC-005 row")
    fields = [field.strip() for field in calc_lines[0].strip().strip("|").split("|")]
    if len(fields) < 5 or fields[3] != "needed":
        raise RuntimeError("CALC-005 was published to CALCULATIONS.md")
    return {
        "deck_sha256": sha256_path(deck_path),
        "calculations_sha256": sha256_path(calculations_path),
        "legacy_actions_cpp_sha256": sha256_path(actions_path),
        "legacy_envrn_cpp_sha256": sha256_path(envrn_path),
        "adsorb_si": {"center_transition": "200->205", "osa_transition": "404->402"},
        "desorb_si": {"center_transition": "205->200", "osa_transition": "402->404"},
        "legacy_si_ladder_kcal_mol": [0.0, 6.0, 12.0, 18.0, 24.0],
        "calc005_status": fields[3],
        "production_value_emitted": False,
    }


def _atom_map_hash(pair: Calc005Pair) -> str:
    def key(origin: Any) -> tuple[str, Any, int | None]:
        return origin.kind, origin.node, origin.ordinal

    payload = {
        "occupied": [key(origin) for origin in pair.occupied_origins],
        "vacancy": [key(origin) for origin in pair.vacancy_origins],
        "silicic_acid": [key(origin) for origin in pair.silicic_acid_origins],
    }
    return _canonical_hash(payload)


def _verify_source(source: dict[str, Any], deck: Path, pair: Calc005Pair) -> None:
    repo = Path(__file__).resolve().parents[2]
    expected = {
        "source_commit": source.get("source_commit"),
        "deck": str(deck),
        "deck_sha256": sha256_path(deck),
        "driver_sha256": sha256_path(
            Path(__file__).with_name("calc005_si_attachment.py")
        ),
        "calc005_sha256": sha256_path(repo / "qm/quarry/calc005.py"),
        "calc005_store_sha256": sha256_path(repo / "qm/quarry/calc005_store.py"),
        "verifier_sha256": sha256_path(Path(__file__)),
        "calculations_sha256": sha256_path(repo / "qm/CALCULATIONS.md"),
        "legacy_actions_cpp_sha256": sha256_path(repo / "legacy/cpp-model/actions.cpp"),
        "legacy_envrn_cpp_sha256": sha256_path(repo / "legacy/cpp-model/envrn.cpp"),
        "atom_map_sha256": _atom_map_hash(pair),
    }
    for key, value in expected.items():
        if source.get(key) != value:
            raise RuntimeError(f"source provenance drifted: {key}")
    git = source.get("git")
    if git is not None and (
        git.get("clean") is not True
        or git.get("head") != source.get("source_commit")
        or git.get("origin_head") != source.get("source_commit")
        or not git.get("branch")
    ):
        raise RuntimeError("committed Git provenance is internally inconsistent")


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


def _validate_production_envelope(envelope: Any) -> None:
    if not isinstance(envelope, dict):
        raise RuntimeError("production verifier execution envelope is missing")
    threads = envelope.get("thread_environment")
    lease = envelope.get("qi2_lease")
    checks = {
        "measured": envelope.get("measured") is True,
        "unit": str(envelope.get("systemd_unit", "")).endswith((".service", ".scope")),
        "cgroup": str(envelope.get("cgroup_path", "")).startswith("/sys/fs/cgroup/"),
        "runtime": envelope.get("runtime_max_seconds") == EXPECTED_RUNTIME_MAX_SECONDS,
        "memory": envelope.get("memory_max_bytes") == EXPECTED_MEMORY_MAX_BYTES,
        "swap": envelope.get("memory_swap_max_bytes") == EXPECTED_MEMORY_SWAP_MAX_BYTES,
        "cpu": envelope.get("cpu_quota_percent") == EXPECTED_CPU_QUOTA_PERCENT,
        "nice": envelope.get("nice") == EXPECTED_NICE,
        "threads": threads
        == {
            "OMP_NUM_THREADS": EXPECTED_THREADS,
            "MKL_NUM_THREADS": EXPECTED_THREADS,
            "OPENBLAS_NUM_THREADS": EXPECTED_THREADS,
        },
        "lease_owner": isinstance(lease, dict)
        and lease.get("owner") == EXPECTED_GPU_OWNER,
        "lease_pid": isinstance(lease, dict) and lease.get("pid") == os.getpid(),
        "lease_ttl": isinstance(lease, dict)
        and lease.get("ttl_hours") == GPU_TTL_HOURS,
        "lease_memory": isinstance(lease, dict)
        and lease.get("expected_gb") == 16.0
        and lease.get("maximum_gb") == 18.0,
        "lease_path": isinstance(lease, dict) and bool(lease.get("path")),
        "shared_service_mutation": envelope.get("shared_service_mutation") is False,
    }
    failed = [name for name, passed in checks.items() if not passed]
    if failed:
        raise RuntimeError(
            "production verifier execution envelope is not exact: " + ", ".join(failed)
        )


def _json_object(path: Path) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise RuntimeError(f"required regular JSON artifact is absent: {path}")
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise RuntimeError(f"JSON artifact is not an object: {path}")
    return payload


def _exact_xyz(cluster: Cluster) -> str:
    lines = [str(len(cluster.symbols)), cluster.name]
    lines.extend(
        f"{symbol} {x:.17g} {y:.17g} {z:.17g}"
        for symbol, (x, y, z) in zip(cluster.symbols, cluster.coords, strict=True)
    )
    return "\n".join(lines) + "\n"


def _origin_value(origin: Any) -> dict[str, Any]:
    return {
        "kind": origin.kind,
        "node": [origin.node[0], list(origin.node[1])],
        "ordinal": origin.ordinal,
    }


def _node_value(node: tuple[int, tuple[int, int, int]]) -> list[Any]:
    return [node[0], list(node[1])]


def _pilot_identity(pair: Calc005Pair) -> dict[str, Any]:
    return {
        "center_node": [4, [0, 0, 0]],
        "metal_shells": 2,
        "x": pair.x,
        "y": pair.y,
        "environment_index": pair.environment_index,
        "occupied_states": pair.occupied_states,
        "vacancy_states": pair.vacancy_states,
        "occupied_frozen_origins": [
            _origin_value(pair.occupied_origins[index])
            for index in pair.occupied.frozen_indices
        ],
        "vacancy_frozen_origins": [
            _origin_value(pair.vacancy_origins[index])
            for index in pair.vacancy.frozen_indices
        ],
        "condensed_topology_mask": {
            "center_bridges": [
                _node_value(node) for node in pair.condensed.center_bridges
            ],
            "kept_center_bridges": [
                _node_value(node) for node in pair.condensed.kept_center_bridges
            ],
            "termination_log": pair.condensed.termination_log,
        },
        "hydrolysis_water_origins": [
            _origin_value(origin)
            for origin in pair.occupied_origins
            if origin.kind.startswith("hydrolysis-water-")
        ],
    }


def _git_blob(repo: Path, commit: str, relative: str) -> bytes:
    return subprocess.run(
        ["git", "-C", str(repo), "show", f"{commit}:{relative}"],
        check=True,
        capture_output=True,
    ).stdout


def _default_restoration_probe(restoration: dict[str, Any]) -> dict[str, Any]:
    lease_path = Path(str(restoration.get("lease_path", ""))).expanduser()
    unit = str(restoration.get("unit", ""))
    output = subprocess.run(
        [
            "systemctl",
            "--user",
            "show",
            unit,
            "--property=ActiveState",
            "--property=SubState",
            "--property=Result",
            "--property=MainPID",
            "--property=ExecMainStatus",
            "--no-pager",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    unit_state = dict(line.split("=", 1) for line in output.splitlines() if "=" in line)
    return {"lease_present": lease_path.exists(), "unit": unit_state}


@contextmanager
def _locked_evidence(output_root: Path) -> Iterator[None]:
    lock_path = output_root / "run.lock"
    if lock_path.is_symlink() or not lock_path.is_file():
        raise RuntimeError("executor run.lock is absent or not a regular file")
    with lock_path.open("r") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_SH)
        try:
            yield
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def _verify_failure_source(
    source: dict[str, Any],
    source_commit: str,
    deck: Path,
    pair: Calc005Pair,
    repo: Path,
    blob_reader: Callable[[Path, str, str], bytes],
) -> dict[str, str]:
    if not re.fullmatch(r"[0-9a-f]{40}", source_commit):
        raise RuntimeError("executor source commit is not a full lowercase Git SHA")
    if source.get("source_commit") != source_commit:
        raise RuntimeError("source map commit differs from terminal commit")
    paths = {
        "deck_sha256": "petra/examples/kaolinite.toml",
        "driver_sha256": "qm/scripts/calc005_si_attachment.py",
        "calc005_sha256": "qm/quarry/calc005.py",
        "calc005_store_sha256": "qm/quarry/calc005_store.py",
        "verifier_sha256": "qm/scripts/calc005_si_attachment_verify.py",
        "calculations_sha256": "qm/CALCULATIONS.md",
        "legacy_actions_cpp_sha256": "legacy/cpp-model/actions.cpp",
        "legacy_envrn_cpp_sha256": "legacy/cpp-model/envrn.cpp",
    }
    pinned = {
        field: hashlib.sha256(blob_reader(repo, source_commit, relative)).hexdigest()
        for field, relative in paths.items()
    }
    independently_pinned = {
        "legacy_actions_cpp_sha256",
        "legacy_envrn_cpp_sha256",
    }
    for field, digest in pinned.items():
        if field not in independently_pinned and source.get(field) != digest:
            raise RuntimeError(f"pinned source hash differs: {field}")
        if field in source and source[field] != digest:
            raise RuntimeError(f"recorded source hash differs: {field}")
    current_paths = {
        "deck_sha256": deck,
        "calculations_sha256": repo / "qm/CALCULATIONS.md",
        "legacy_actions_cpp_sha256": repo / "legacy/cpp-model/actions.cpp",
        "legacy_envrn_cpp_sha256": repo / "legacy/cpp-model/envrn.cpp",
    }
    for field, path in current_paths.items():
        if sha256_path(path) != pinned[field]:
            raise RuntimeError(f"live source drifted from execution source: {field}")
    git = source.get("git")
    tracked = git.get("tracked_source_files") if isinstance(git, dict) else None
    core_tracked = {
        "petra/examples/kaolinite.toml",
        "qm/CALCULATIONS.md",
        "qm/quarry/calc005.py",
        "qm/quarry/calc005_store.py",
        "qm/scripts/calc005_si_attachment.py",
        "qm/scripts/calc005_si_attachment_verify.py",
    }
    allowed_tracked = core_tracked | {
        "legacy/cpp-model/actions.cpp",
        "legacy/cpp-model/envrn.cpp",
    }
    if (
        not isinstance(git, dict)
        or git.get("clean") is not True
        or git.get("head") != source_commit
        or git.get("origin_head") != source_commit
        or not git.get("branch")
        or not isinstance(tracked, list)
        or len(tracked) != len(set(tracked))
        or not core_tracked.issubset(tracked)
        or not set(tracked).issubset(allowed_tracked)
    ):
        raise RuntimeError("executor Git source provenance is not exact")
    if source.get("atom_map_sha256") != _atom_map_hash(pair):
        raise RuntimeError("reconstructed atom map differs from executor source")
    if source.get("condensed_geometry_hash") != geometry_hash(
        _exact_xyz(pair.condensed.cluster)
    ):
        raise RuntimeError("reconstructed condensed geometry differs")
    if source.get("pilot_identity") != _pilot_identity(pair):
        raise RuntimeError("reconstructed physical-state identity differs")
    return pinned


def _verify_failure_namespace(output_root: Path) -> dict[str, str]:
    component_root = output_root / "components"
    c_root = component_root / "C"
    required = {"reservation.json", "optimizer-entered.json", "receipt.json"}
    observed = {path.name for path in c_root.iterdir()}
    if observed != required or any(
        path.is_symlink() or not path.is_file() for path in c_root.iterdir()
    ):
        raise RuntimeError("failed C checkpoint namespace is not exact")
    for role in ("V", "SiOH4"):
        if (component_root / role).exists():
            raise RuntimeError(f"unspent component unexpectedly exists: {role}")
    forbidden_names = {
        "raw-endpoint.xyz",
        "endpoint.xyz",
        "calc005-result.json",
        "calculation-receipt.json",
        "store.sqlite",
    }
    for path in output_root.rglob("*"):
        if path.is_symlink():
            raise RuntimeError(f"evidence namespace contains a symlink: {path}")
        if path.name in forbidden_names or path.name.startswith(".staging-"):
            raise RuntimeError(f"forbidden executor artifact exists: {path}")
    generations = output_root / "generations"
    if generations.exists():
        raise RuntimeError("executor generation namespace exists after failed pilot")
    return {name: sha256_path(c_root / name) for name in sorted(required)}


def _verify_failure_restoration(
    terminal: dict[str, Any],
    probe: Callable[[dict[str, Any]], dict[str, Any]],
) -> dict[str, Any]:
    restoration = terminal.get("restoration")
    if not isinstance(restoration, dict):
        raise RuntimeError("failure terminal lacks restoration evidence")
    required = {
        "status": "verified-restored",
        "bootstrap_session_closed": True,
        "qi2_lease_released": True,
        "unit_active_state": "failed",
        "unit_result": "exit-code",
        "unit_exec_main_status": 2,
        "unit_main_pid": 0,
    }
    if any(restoration.get(key) != value for key, value in required.items()):
        raise RuntimeError("recorded executor restoration state is not exact")
    if not restoration.get("unit") or not restoration.get("lease_path"):
        raise RuntimeError("restoration unit or lease path is missing")
    prior = dict(terminal)
    prior.pop("restoration")
    prior_bytes = (json.dumps(prior, indent=2, sort_keys=True) + "\n").encode()
    prior_hash = hashlib.sha256(prior_bytes).hexdigest()
    if restoration.get("prior_terminal_sha256") != prior_hash:
        raise RuntimeError("pre-restoration terminal hash does not reproduce")
    if datetime.fromisoformat(restoration["checked_at"]) < datetime.fromisoformat(
        terminal["written_at"]
    ):
        raise RuntimeError("restoration predates the terminal failure")
    live = probe(restoration)
    expected_unit = {
        "MainPID": "0",
        "Result": "exit-code",
        "ExecMainStatus": "2",
        "ActiveState": "failed",
        "SubState": "failed",
    }
    if live != {"lease_present": False, "unit": expected_unit}:
        raise RuntimeError("live executor restoration state differs from receipt")
    return {"prior_terminal_sha256": prior_hash, **live}


def verify_failed_pilot(
    output_root: Path,
    deck: Path,
    *,
    verifier_identity: str,
    repo: Path | None = None,
    blob_reader: Callable[[Path, str, str], bytes] = _git_blob,
    restoration_probe: Callable[[dict[str, Any]], dict[str, Any]] = (
        _default_restoration_probe
    ),
) -> dict[str, Any]:
    """Adjudicate a terminal executor failure without calculator or optimizer calls."""

    output_root = output_root.resolve()
    deck = deck.resolve()
    repo = (repo or Path(__file__).resolve().parents[2]).resolve()
    result: dict[str, Any]
    with _locked_evidence(output_root):
        try:
            terminal_path = output_root / "terminal-receipt.json"
            terminal = _json_object(terminal_path)
            terminal_keys = {
                "schema",
                "written_at",
                "verdict",
                "failed_component",
                "detail",
                "executor_identity",
                "source_commit",
                "source",
                "components",
                "optimizer_budget",
                "quarantine",
                "canonical_value_exposed",
                "forbidden_outputs_emitted",
                "independent_verification_required",
                "restoration",
            }
            if set(terminal) != terminal_keys:
                raise RuntimeError("failure terminal field set is not exact")
            executor = terminal.get("executor_identity")
            if not verifier_identity or verifier_identity == executor:
                raise RuntimeError(
                    "verifier identity must differ from executor identity"
                )
            if (
                terminal.get("schema") != "calc005-terminal-receipt-v1"
                or terminal.get("verdict") != "incomplete-computational-failure"
                or terminal.get("failed_component") != "C"
                or terminal.get("quarantine") != "spent-checkpoint-retained-in-place"
                or terminal.get("canonical_value_exposed") is not False
                or terminal.get("forbidden_outputs_emitted") is not False
                or terminal.get("independent_verification_required") is not True
            ):
                raise RuntimeError("executor failure terminal identity is invalid")
            expected_budget = {
                "C": {"calls": 1, "max_steps": 150, "retries": 0},
                "V": {"calls": 0, "max_steps": 150, "retries": 0},
                "SiOH4": {"calls": 0, "max_steps": 150, "retries": 0},
            }
            if terminal.get("optimizer_budget") != expected_budget:
                raise RuntimeError("executor optimizer budget is not exact")

            checkpoint_hashes = _verify_failure_namespace(output_root)
            c_root = output_root / "components/C"
            reservation = _json_object(c_root / "reservation.json")
            entered = _json_object(c_root / "optimizer-entered.json")
            receipt = _json_object(c_root / "receipt.json")
            signature = receipt.get("signature")
            if (
                set(reservation) != {"schema", "status", "reserved_at", "signature"}
                or set(entered) != {"schema", "status", "entered_at", "signature"}
                or set(receipt)
                != {"schema", "status", "signature", "optimizer", "detail"}
                or reservation.get("schema") != "calc005-component-reservation-v1"
                or reservation.get("status") != "optimizer-budget-reserved"
                or entered.get("schema") != "calc005-optimizer-entered-v1"
                or entered.get("status") != "optimizer-entered"
                or receipt.get("schema") != "calc005-component-receipt-v1"
                or receipt.get("status") != "optimizer-failed"
                or not isinstance(signature, dict)
                or reservation.get("signature") != signature
                or entered.get("signature") != signature
                or terminal.get("components") != {"C": receipt}
            ):
                raise RuntimeError("failed component checkpoint chain is invalid")
            if receipt.get("optimizer") != {
                "converged": False,
                "observed_calls": 1,
                "observed_retries": 0,
                "observed_max_steps": 150,
            } or signature.get("optimizer") != {
                "max_steps": 150,
                "retry_allowed": False,
            }:
                raise RuntimeError("failed component optimizer ledger is invalid")
            if "Nuclear gradients" not in receipt.get("detail", "") or (
                "not converged" not in receipt.get("detail", "")
            ):
                raise RuntimeError(
                    "failed component does not preserve the SCF-gradient failure"
                )

            source_settings = _json_object(output_root / "source-settings.json")
            source = terminal.get("source")
            if (
                not isinstance(source, dict)
                or source_settings.get("source") != source
                or signature.get("source") != source
            ):
                raise RuntimeError("source provenance copies differ")
            geometry_settings, production_settings = _settings(True)
            expected_settings = {
                "geometry": asdict(geometry_settings),
                "geometry_method": R2SCAN3C_METHOD,
                "geometry_sha256": _canonical_hash(asdict(geometry_settings)),
                "geometry_fingerprint": frequency_settings_fingerprint(
                    geometry_settings
                ),
                "production": asdict(production_settings),
                "production_method": PRODUCTION_METHOD,
                "production_sha256": _canonical_hash(asdict(production_settings)),
                "production_fingerprint": frequency_settings_fingerprint(
                    production_settings
                ),
                "temperature_k": TEMPERATURE_K,
            }
            if source_settings.get("settings") != expected_settings:
                raise RuntimeError("executor settings receipt drifted")
            if signature.get("settings") != {
                "geometry": expected_settings["geometry"],
                "geometry_sha256": expected_settings["geometry_sha256"],
                "production": expected_settings["production"],
                "production_sha256": expected_settings["production_sha256"],
            }:
                raise RuntimeError("component settings differ from source settings")

            pair = build_calc005_pair(deck, environment_index=1, metal_shells=2)
            validation = pair.validate()
            if any(
                value is not True
                for value in validation.values()
                if isinstance(value, bool)
            ):
                raise RuntimeError(
                    "reconstructed CALC-005 pair failed a physical-state gate"
                )
            if (
                pair.condensed.cluster.formula != "Al6H36O29Si"
                or pair.occupied.formula != "Al6H38O30Si"
                or pair.vacancy.formula != "Al6H34O26"
                or pair.silicic_acid.formula != "H4O4Si"
                or signature.get("formula") != pair.occupied.formula
                or signature.get("charge") != 0
                or signature.get("spin") != 0
                or signature.get("frozen_indices") != pair.occupied.frozen_indices
                or signature.get("seed_geometry_fingerprint")
                != frequency_geometry_fingerprint(pair.occupied)
            ):
                raise RuntimeError("reconstructed component identity differs")
            pinned_hashes = _verify_failure_source(
                source,
                str(terminal["source_commit"]),
                deck,
                pair,
                repo,
                blob_reader,
            )
            petra = validate_petra_boundary(repo)

            log_path = output_root / "logs/a3i-executor.log"
            if log_path.is_symlink() or not log_path.is_file():
                raise RuntimeError("canonical executor log is absent")
            log_text = log_path.read_text()
            step_lines = re.findall(r"^Step\s+\d+\s*:", log_text, flags=re.MULTILINE)
            if (
                log_text.count("geomeTRIC started.") != 1
                or log_text.count("maxiter                   150") != 1
                or step_lines != ["Step    0 :"]
            ):
                raise RuntimeError(
                    "executor log does not prove one step-0 optimizer entry"
                )
            json_start = log_text.rfind("\n{")
            if json_start < 0 or json.loads(log_text[json_start + 1 :]) != {
                key: value for key, value in terminal.items() if key != "restoration"
            }:
                raise RuntimeError("executor log terminal payload differs")

            restoration = _verify_failure_restoration(terminal, restoration_probe)
            result = {
                "schema": "calc005-verified-terminal-v2",
                "written_at": now(),
                "status": "verified-pass",
                "pilot_disposition": "terminally-rejected",
                "executor_verdict": terminal["verdict"],
                "executor_identity": executor,
                "verifier_identity": verifier_identity,
                "executor_terminal_sha256": sha256_path(terminal_path),
                "checkpoint_sha256": checkpoint_hashes,
                "executor_log_sha256": sha256_path(log_path),
                "source_commit": terminal["source_commit"],
                "pinned_source_sha256": pinned_hashes,
                "physical_state": {
                    "condensed_state": 204,
                    "live_states": pair.occupied_states,
                    "vacancy_states": pair.vacancy_states,
                    "cycle": "Al6H38O30Si -> Al6H34O26 + H4O4Si",
                    "water_thermochemical_coefficient": 0,
                    "pair_validation": validation,
                },
                "optimizer_budget": expected_budget,
                "optimizer_calls": 0,
                "calculator_calls": 0,
                "petra_boundary": petra,
                "restoration": restoration,
                "canonical_value_exposed": False,
                "forbidden_outputs_emitted": False,
                "failed_edges": [],
                "artifact_writes": ["verified-terminal.json"],
            }
        except Exception as exc:
            result = {
                "schema": "calc005-verified-terminal-v2",
                "written_at": now(),
                "status": "verified-fail",
                "pilot_disposition": "evidence-rejected",
                "verifier_identity": verifier_identity,
                "detail": f"{type(exc).__name__}: {exc}",
                "failed_edges": [f"{type(exc).__name__}: {exc}"],
                "optimizer_calls": 0,
                "calculator_calls": 0,
                "artifact_writes": ["verified-terminal.json"],
            }
        atomic_json(output_root / "verified-terminal.json", result)
    return result


def measure_production_envelope() -> dict[str, Any]:
    """Measure the verifier's own process, cgroup, unit, threads, and QI2 lease."""

    cgroup_lines = Path("/proc/self/cgroup").read_text().splitlines()
    unified = [line.split("::", 1)[1] for line in cgroup_lines if "::" in line]
    if len(unified) != 1:
        raise RuntimeError("CALC-005 verifier requires one cgroup-v2 systemd unit")
    relative = unified[0].lstrip("/")
    cgroup = Path("/sys/fs/cgroup") / relative
    units = [
        part for part in Path(relative).parts if part.endswith((".service", ".scope"))
    ]
    if not units:
        raise RuntimeError("CALC-005 verifier must run inside a systemd service/scope")
    unit = units[-1]

    def integer_limit(name: str) -> int:
        value = (cgroup / name).read_text().strip()
        if value == "max":
            raise RuntimeError(f"CALC-005 verifier requires a finite {name}")
        return int(value)

    cpu_fields = (cgroup / "cpu.max").read_text().split()
    if len(cpu_fields) != 2 or cpu_fields[0] == "max":
        raise RuntimeError("CALC-005 verifier requires a finite cpu.max quota")
    runtime_text = subprocess.run(
        [
            "systemctl",
            "--user",
            "show",
            unit,
            "--property=RuntimeMaxUSec",
            "--value",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    lease_path = Path(
        os.environ.get(
            "GPU_LEASE_PATH", str(Path.home() / ".local/state/gpu-lease/lease.json")
        )
    ).expanduser()
    lease = json.loads(lease_path.read_text())
    envelope = {
        "measured": True,
        "systemd_unit": unit,
        "cgroup_path": str(cgroup),
        "runtime_max_seconds": _timespan_seconds(runtime_text),
        "memory_max_bytes": integer_limit("memory.max"),
        "memory_swap_max_bytes": integer_limit("memory.swap.max"),
        "cpu_quota_percent": int(round(int(cpu_fields[0]) / int(cpu_fields[1]) * 100)),
        "nice": os.getpriority(os.PRIO_PROCESS, 0),
        "thread_environment": {
            name: int(os.environ.get(name, "0"))
            for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")
        },
        "qi2_lease": {
            "path": str(lease_path),
            "owner": lease["owner"],
            "pid": lease["pid"],
            "ttl_hours": lease["ttl"],
            "expected_gb": lease["expected_gb"],
            "maximum_gb": 18.0,
        },
        "shared_service_mutation": False,
    }
    _validate_production_envelope(envelope)
    return envelope


def verify_pilot(
    output_root: Path,
    deck: Path,
    backend: CalculatorBackend,
    *,
    verifier_identity: str,
    execution_envelope: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Recompute all scientific gates from raw artifacts without executor helpers."""

    if _is_production_backend(backend):
        _validate_production_envelope(execution_envelope)
        verified_envelope = execution_envelope
    else:
        if execution_envelope != {"mode": "in-process-test"}:
            raise RuntimeError("fake verifier backend must run explicitly in-process")
        verified_envelope = execution_envelope
    output_root = output_root.resolve()
    try:
        terminal_path = output_root / "terminal-receipt.json"
        terminal = json.loads(terminal_path.read_text())
        if (
            terminal.get("schema") != "calc005-terminal-receipt-v2"
            or terminal.get("verdict") != "passed-protocol-pilot"
            or terminal.get("classification") != CLASSIFICATION
            or terminal.get("receipt_payload_sha256") != _terminal_hash(terminal)
        ):
            raise RuntimeError("executor terminal identity/hash is invalid")
        executor = terminal.get("executor_identity")
        if not verifier_identity or verifier_identity == executor:
            raise RuntimeError("verifier identity must differ from executor identity")
        artifacts = terminal.get("artifacts", {})
        for name in ("result", "calculation_receipt", "store"):
            artifact = artifacts.get(name, {})
            path = Path(str(artifact.get("path", "")))
            if not path.is_file() or artifact.get("sha256") != sha256_path(path):
                raise RuntimeError(f"terminal {name} artifact hash mismatch")
        calculation = json.loads(
            Path(artifacts["calculation_receipt"]["path"]).read_text()
        )
        if _canonical_hash(calculation) != artifacts["calculation_receipt"].get(
            "canonical_sha256"
        ):
            raise RuntimeError("calculation receipt canonical hash mismatch")
        if (
            calculation.get("schema") != "calc005-calculation-receipt-v1"
            or calculation.get("verdict") != "passed-protocol-pilot"
            or calculation.get("classification") != CLASSIFICATION
            or calculation.get("temperature_k") != TEMPERATURE_K
            or calculation.get("stoichiometry") != STOICHIOMETRY
            or calculation.get("water_thermochemical_coefficient") != 0
            or terminal.get("source_commit") != calculation.get("source_commit")
        ):
            raise RuntimeError("immutable calculation identity drifted")
        pair = build_calc005_pair(deck.resolve(), environment_index=1, metal_shells=2)
        _verify_source(calculation["source"], deck.resolve(), pair)
        geometry_settings, production_settings = _settings(
            bool(calculation["settings"]["geometry"]["use_gpu"])
        )
        expected_settings = {
            "geometry": asdict(geometry_settings),
            "geometry_method": R2SCAN3C_METHOD,
            "geometry_fingerprint": frequency_settings_fingerprint(geometry_settings),
            "production": asdict(production_settings),
            "production_method": PRODUCTION_METHOD,
            "production_fingerprint": frequency_settings_fingerprint(
                production_settings
            ),
            "temperature_k": TEMPERATURE_K,
        }
        for key, value in expected_settings.items():
            if calculation["settings"].get(key) != value:
                raise RuntimeError(f"calculation settings drifted: {key}")
        templates = {"C": pair.occupied, "V": pair.vacancy, "SiOH4": pair.silicic_acid}
        origins = {
            "C": pair.occupied_origins,
            "V": pair.vacancy_origins,
            "SiOH4": pair.silicic_acid_origins,
        }
        components = calculation.get("components", {})
        if set(components) != set(STOICHIOMETRY):
            raise RuntimeError("calculation component set drifted")
        component_terms: dict[str, dict[str, float]] = {}
        component_hashes: dict[str, str] = {}
        for role in ("C", "V", "SiOH4"):
            claimed = components[role]
            if claimed.get("status") != "accepted":
                raise RuntimeError(f"{role} status is not accepted")
            if claimed.get("optimizer") != {
                "converged": True,
                "observed_calls": 1,
                "observed_retries": 0,
                "observed_max_steps": 150,
            }:
                raise RuntimeError(f"{role} optimizer evidence drifted")
            if claimed.get("checkpoint_sha256") != _checkpoint_hash(claimed):
                raise RuntimeError(f"{role} checkpoint payload hash drifted")
            endpoint_path = Path(claimed["endpoint"]["path"])
            raw_path = Path(claimed["raw_endpoint"]["path"])
            for label, path in (("endpoint", endpoint_path), ("raw", raw_path)):
                if (
                    not path.is_file()
                    or sha256_path(path)
                    != claimed[f"{label}_endpoint" if label == "raw" else label][
                        "sha256"
                    ]
                ):
                    raise RuntimeError(f"{role} {label} artifact hash mismatch")
            endpoint = _read_xyz(endpoint_path, templates[role])
            raw = _read_xyz(raw_path, templates[role])
            _structure_gate(role, raw, templates[role], origins[role])
            structure = _structure_gate(role, endpoint, templates[role], origins[role])
            if (
                claimed.get("structure", {}).get("center_si_o_bond_count")
                != structure["center_si_o_bond_count"]
            ):
                raise RuntimeError(f"{role} structural receipt center-bond drifted")
            gradient = _gradient_gate(
                role, endpoint, backend.gradient(role, endpoint, geometry_settings)
            )
            if not _semantic_equal(
                {
                    "rms_hartree_per_bohr": claimed["gradient"]["rms_hartree_per_bohr"],
                    "max_hartree_per_bohr": claimed["gradient"]["max_hartree_per_bohr"],
                },
                gradient,
                tolerance=1.0e-12,
            ):
                raise RuntimeError(f"{role} independent gradient differs")
            frequency_result = backend.frequencies(role, endpoint, geometry_settings)
            frequency = _frequency_gate(
                role, frequency_result, endpoint, geometry_settings
            )
            frequency_claim = claimed["frequency"]
            comparable = {key: frequency_claim[key] for key in frequency}
            if not _semantic_equal(comparable, frequency, tolerance=1.0e-7):
                raise RuntimeError(f"{role} independent frequency evidence differs")
            production_energy = float(
                backend.energy(role, endpoint, production_settings)
            )
            single_point = claimed.get("production_single_point", {})
            if (
                not math.isfinite(production_energy)
                or single_point.get("converged") is not True
                or single_point.get("method") != PRODUCTION_METHOD
                or not math.isclose(
                    production_energy,
                    float(claimed["production_electronic_hartree"]),
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-10,
                )
            ):
                raise RuntimeError(f"{role} independent production energy/SP differs")
            terms = _derive_thermal_terms(role, frequency_result, production_energy)
            if not _semantic_equal(terms, claimed["thermochemistry"], tolerance=1.0e-8):
                raise RuntimeError(f"{role} independent thermochemistry differs")
            component_terms[role] = terms
            component_hashes[role] = claimed["endpoint"]["sha256"]
        recomposed = _compose(component_terms)
        if not _semantic_equal(
            recomposed, calculation["thermochemistry"], tolerance=1.0e-8
        ):
            raise RuntimeError("independent CALC-005 recomposition differs")
        log_relative_rate_ratio = recomposed["s_10_kj_mol"] / (
            R_KJ_MOL_K * TEMPERATURE_K
        )
        relative_rate_ratio = math.exp(log_relative_rate_ratio)
        if not all(
            math.isfinite(value)
            for value in (log_relative_rate_ratio, relative_rate_ratio)
        ):
            raise RuntimeError("CALC-005 detailed-balance transform is non-finite")
        store_path = Path(artifacts["store"]["path"])
        store = validate_calc005_store(store_path, terminal, pair)
        petra = validate_petra_boundary(Path(__file__).resolve().parents[2])
        if terminal.get("forbidden_outputs_emitted") is not False:
            raise RuntimeError("executor claimed a forbidden production emission")
        result = {
            "schema": "calc005-verified-terminal-v2",
            "written_at": now(),
            "status": "verified-pass",
            "classification": CLASSIFICATION,
            "executor_identity": executor,
            "verifier_identity": verifier_identity,
            "executor_terminal_sha256": sha256_path(terminal_path),
            "calculation_receipt_sha256": artifacts["calculation_receipt"][
                "canonical_sha256"
            ],
            "store_sha256": artifacts["store"]["sha256"],
            "component_endpoint_sha256": component_hashes,
            "thermochemistry": recomposed,
            "detailed_balance": {
                "equation": "(katt_10/kdet_10)/(katt_00/kdet_00)=exp(S_10/RT)",
                "log_relative_rate_ratio": log_relative_rate_ratio,
                "relative_rate_ratio": relative_rate_ratio,
                "temperature_k": TEMPERATURE_K,
            },
            "store": store,
            "petra_boundary": petra,
            "optimizer_calls": 0,
            "artifact_writes": ["verified-terminal.json"],
            "execution_envelope": verified_envelope,
        }
    except Exception as exc:
        result = {
            "schema": "calc005-verified-terminal-v2",
            "written_at": now(),
            "status": "verified-fail",
            "verifier_identity": verifier_identity,
            "detail": f"{type(exc).__name__}: {exc}",
            "optimizer_calls": 0,
            "artifact_writes": ["verified-terminal.json"],
            "execution_envelope": verified_envelope,
        }
    atomic_json(output_root / "verified-terminal.json", result)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("/mnt/data/vsletten/dissertation-data/a3h-calc005-si-n1-pilot"),
    )
    parser.add_argument(
        "--deck",
        type=Path,
        default=Path(__file__).resolve().parents[2] / "petra/examples/kaolinite.toml",
    )
    parser.add_argument("--verifier-identity", required=True)
    parser.add_argument("--gpu", action="store_true")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--nice", type=int, default=10)
    parser.add_argument("--log")
    args = parser.parse_args(argv)
    terminal_path = args.output_root.resolve() / "terminal-receipt.json"
    try:
        terminal = _json_object(terminal_path)
    except (OSError, json.JSONDecodeError, RuntimeError) as exc:
        parser.error(f"cannot inspect executor terminal receipt: {exc}")
    if terminal.get("schema") == "calc005-terminal-receipt-v1" and terminal.get(
        "verdict"
    ) in {"incomplete-computational-failure", "rejected-physical-state"}:
        result = verify_failed_pilot(
            args.output_root,
            args.deck,
            verifier_identity=args.verifier_identity,
        )
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0 if result["status"] == "verified-pass" else 1
    if not args.gpu:
        parser.error("production evidence recomputation requires --gpu")
    if args.threads != EXPECTED_THREADS or args.nice != 0:
        parser.error(
            "production verification requires --threads 16 --nice 0 "
            "inside a Nice=10 unit"
        )
    if argv is None:
        from quarry.etiquette import bootstrap_cli

        bootstrap_cli(
            "calc005_si_attachment_verify",
            default_run_root=Path(
                "/mnt/data/vsletten/dissertation-data/a3h-calc005-si-n1-pilot/logs"
            ),
            gpu_owner="calc005_si_attachment_verify",
            gpu_ttl_hours=GPU_TTL_HOURS,
        )
    execution_envelope = measure_production_envelope()
    result = verify_pilot(
        args.output_root,
        args.deck,
        PipelineBackend(),
        verifier_identity=args.verifier_identity,
        execution_envelope=execution_envelope,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "verified-pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
