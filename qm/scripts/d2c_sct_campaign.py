#!/usr/bin/env python
"""D2c SCT campaign preflight and immutable resume contract.

Only ``--dry-run`` is implemented.  The command verifies the frozen D2b bundle,
validates all four transition-state atom mappings against the existing D2b reaction
templates, and writes one non-overwriting pending receipt.  It deliberately performs
no IRC, Hessian, high-level, VAG, or tunnelling calculation.
"""

from __future__ import annotations

import argparse
import ctypes
import errno
import fcntl
import hashlib
import importlib.metadata
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

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

import numpy as np  # noqa: E402

from quarry.clusters import Cluster  # noqa: E402
from quarry.native_hessian import HARTREE_TO_EV, NativeHessianResult  # noqa: E402
from quarry.pipeline import (  # noqa: E402
    BOHR_TO_ANGSTROM,
    frequency_geometry_fingerprint,
)
from quarry.reaction_path import (  # noqa: E402
    MassScaledPath,
    build_mass_scaled_path,
    hessian_eigenvalues_to_wavenumbers_cm,
    project_vibrational_hessian,
)
from quarry.ts import IrcPoint, SellaIrcTrace  # noqa: E402
from scripts import d2c_input_bundle
from scripts.production_energetics import load_xyz_like
from scripts.surface_rate_protocol import reactions

SCHEMA = "d2c-sct-campaign-preflight-v3"
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
GPU4PYSCF_DISTRIBUTIONS = (
    "gpu4pyscf-cuda12x",
    "gpu4pyscf-cuda11x",
    "gpu4pyscf",
)
CUPY_DISTRIBUTIONS = ("cupy-cuda13x", "cupy-cuda12x", "cupy-cuda11x", "cupy")
DEPENDENCY_VERSION_KEYS = frozenset(
    (
        *DEPENDENCY_DISTRIBUTIONS,
        "gpu_backend",
        "gpu4pyscf_distribution",
        "gpu4pyscf",
        "cupy_distribution",
        "cupy",
        "cuda_runtime",
        "cuda_driver",
        "cuda_device_count",
    )
)
# These geometry hashes are trusted source constants, deliberately independent of
# the mutable manifest being verified.  They bind atom row identity even if a copied
# bundle and all of its ordinary byte receipts are coherently rewritten.
TRUSTED_CANONICAL_TS_GEOMETRY_SHA256 = {
    "h-co-1w-oside": "7c121ddaddfb47932b95301593524498fdb7ed98faaf741ee383e04c49880106",
    "h-co-1w-cside": "4265cd1ee5dfcc082027bafa8322587aa505f200589adcaf9cbb97e5b7eb7f2b",
    "h-h2co-ch3o-1w": (
        "c5602e65b972f2703d467781cca1c8e1299eb06d68d13fd3889cbf75c294d6bf"
    ),
    "h-h2co-h2-hco-1w": (
        "a1b804833d74cd77b7382cd516b84af7981fd37f18417ae7a858a7898a8225bf"
    ),
}
# Independent, source-reviewed endpoint constants. Neither these values nor their
# chemistry labels are inferred from the mutable bundle manifest or Sella direction.
TRUSTED_FROZEN_ENDPOINT_EVIDENCE: dict[str, dict[str, dict[str, str]]] = {
    "h-co-1w-oside": {
        "irc_fwd.xyz": {
            "basin": "reactant",
            "geometry_sha256": (
                "b86e012fdbf6e28e83eed17b802629358817a8215294a50f947a30a371bf5ce3"
            ),
            "file_sha256": (
                "4ac3c871b4177c921926e73ce454c01d2af4cb6c1f835e22ac487d7f64bcd448"
            ),
        },
        "irc_back.xyz": {
            "basin": "product",
            "geometry_sha256": (
                "4d5e450139117b75c0d0161331caf78a49e639b8c9fa7f6e09c30977f02f1b0b"
            ),
            "file_sha256": (
                "8e8bbd3a34a5e98fd95cdf3afe873112a7dc5a149b2bcb634a125a20a5dd2854"
            ),
        },
    },
    "h-co-1w-cside": {
        "irc_fwd.xyz": {
            "basin": "product",
            "geometry_sha256": (
                "9a678f1ed4551143aca887d7d66bcaee5977a88392d8624289013f934e3fea9a"
            ),
            "file_sha256": (
                "b980041316556cc2fc150e3887d073805bbc6c6f8e6b74b2f65722a5815c339c"
            ),
        },
        "irc_back.xyz": {
            "basin": "reactant",
            "geometry_sha256": (
                "2a1018e2f00bfb5dad4171fa66b7a4d224fc6d6111e769f5f8e67cf10be290a4"
            ),
            "file_sha256": (
                "d4025f53cbfa5723cb0aaff25dcd8539134ed462a82098f82a205dc09ddc8eec"
            ),
        },
    },
    "h-h2co-ch3o-1w": {
        "irc_fwd.xyz": {
            "basin": "reactant",
            "geometry_sha256": (
                "28243da36a389ff9f7ea0548887de4e7a77c6c1da0503b551c394a3e170b116a"
            ),
            "file_sha256": (
                "63c01a51689b5bce687786733dc699601512195e66cbc4d5f0b9d3313e11ebf2"
            ),
        },
        "irc_back.xyz": {
            "basin": "product",
            "geometry_sha256": (
                "dc27197ec50d88424a9f35500a51404fd2df0a77786a3c20dbfad13332f8301f"
            ),
            "file_sha256": (
                "21340d103a01f7b5206d326dde7dcc8fa5aadeb93b1626768e8aec2cb45c4130"
            ),
        },
    },
    "h-h2co-h2-hco-1w": {
        "irc_fwd.xyz": {
            "basin": "reactant",
            "geometry_sha256": (
                "8e038fc5c5b5c63460f9753e6be33b37b6eea5864d1e16455beb9b603582d6a5"
            ),
            "file_sha256": (
                "7a7a4efb4f1dc1ce53ada54fcf467a932e0a1ba82dfbc4730493bea47bebd7ef"
            ),
        },
        "irc_back.xyz": {
            "basin": "product",
            "geometry_sha256": (
                "3e6f426a379fc62f34402f37232d6d858bb5fd2b80687d0eaefac0417a8231e4"
            ),
            "file_sha256": (
                "0ca6a15447617548ca58ce12975d53e43e64d718af7579592b57f37fe0a0a796"
            ),
        },
    },
}

ENDPOINT_ROUTE_STATES: dict[str, dict[str, Any]] = {
    "h-co-1w-oside": {
        "symbols": ("C", "O", "H", "O", "H", "H"),
        "charge": 0,
        "spin": 1,
        "frozen_indices": (),
    },
    "h-co-1w-cside": {
        "symbols": ("C", "O", "H", "O", "H", "H"),
        "charge": 0,
        "spin": 1,
        "frozen_indices": (),
    },
    "h-h2co-ch3o-1w": {
        "symbols": ("C", "O", "H", "H", "H", "O", "H", "H"),
        "charge": 0,
        "spin": 1,
        "frozen_indices": (),
    },
    "h-h2co-h2-hco-1w": {
        "symbols": ("C", "O", "H", "H", "H", "O", "H", "H"),
        "charge": 0,
        "spin": 1,
        "frozen_indices": (),
    },
}

# Distances below the first threshold are covalent and distances above the second
# are non-covalent. The closed interval is deliberately unclassifiable.
ENDPOINT_PAIR_DISTANCE_POLICY: dict[tuple[str, str], tuple[float, float]] = {
    ("C", "C"): (1.70, 2.00),
    ("C", "H"): (1.25, 1.55),
    ("C", "O"): (1.50, 1.75),
    ("H", "H"): (0.90, 1.05),
    ("H", "O"): (1.15, 1.35),
    ("O", "O"): (1.60, 1.85),
}
ENDPOINT_COLLISION_FLOOR_ANGSTROM = 0.60

_CO_REACTANT_GRAPH = ((0, 1), (3, 4), (3, 5))
_CO_PRODUCT_GRAPH = ((0, 1), (0, 2), (3, 4), (3, 5))
_H2CO_REACTANT_GRAPH = ((0, 1), (0, 2), (0, 3), (5, 6), (5, 7))
ENDPOINT_GRAPH_CONTRACTS: dict[
    str, dict[Literal["reactant", "product"], tuple[tuple[int, int], ...]]
] = {
    "h-co-1w-oside": {
        "reactant": _CO_REACTANT_GRAPH,
        "product": _CO_PRODUCT_GRAPH,
    },
    "h-co-1w-cside": {
        "reactant": _CO_REACTANT_GRAPH,
        "product": _CO_PRODUCT_GRAPH,
    },
    "h-h2co-ch3o-1w": {
        "reactant": _H2CO_REACTANT_GRAPH,
        "product": ((0, 1), (0, 2), (0, 3), (0, 4), (5, 6), (5, 7)),
    },
    "h-h2co-h2-hco-1w": {
        "reactant": _H2CO_REACTANT_GRAPH,
        "product": ((0, 1), (0, 3), (2, 4), (5, 6), (5, 7)),
    },
}
# Ground-state isotopic masses.  Carbon-12 is exact by definition; the others are
# the neutral-atom masses used by this bounded H/C/O campaign.
ISOTOPIC_MASSES_AMU = {
    "H": 1.00782503223,
    "C": 12.0,
    "O": 15.99491461957,
}
REFERENCE_MASS_AMU = 1.0
BOUNDS: dict[str, Any] = {
    "transition_state_qualification": {
        "physical_fmax_ev_per_angstrom_exclusive_maximum": 0.02,
        "negative_eigenvalue_tolerance_hartree_per_bohr2_amu": 1.0e-8,
        "required_significant_imaginary_mode_count": 1,
        "minimum_reaction_imaginary_wavenumber_cm": 200.0,
        "minimum_mapped_reaction_vector_overlap": 0.20,
        "minimum_irc_tangent_overlap": 0.80,
        "spectator_imaginary_modes_forbidden": True,
        "full_index_gate_precedes_transverse_projection": True,
    },
    "irc": {
        "directions": ["forward", "reverse"],
        "algorithm": "sella-gonzalez-schlegel",
        "step_size_angstrom": 0.05,
        "maximum_steps_per_direction": 200,
        "outer_fmax_ev_per_angstrom": 0.05,
        "inner_fmax_ev_per_angstrom": 0.01,
    },
    "hessian": {
        "coverage": "every retained IRC point including TS and endpoints",
        "maximum_retained_points_per_direction": 201,
        "cartesian_hessian_units": "hartree / bohr^2",
        "transverse_negative_eigenvalue_tolerance": 1.0e-8,
        "transverse_negative_eigenvalue_tolerance_units": "hartree / bohr^2 / amu",
        "transverse_eigenvalue_units": "hartree / bohr^2 / amu",
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
        "straightness_tolerance_per_angstrom": 1.0e-12,
    },
}
ROUTE_STAGE_CONTRACT: dict[str, dict[str, Any]] = {
    "transition_state_qualification": {
        "required": True,
        "receipt": "ts-qualification/receipt.json",
        "requirement": (
            "fresh physical gradient and strict full 3N-6 first-order-saddle gate"
        ),
    },
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
    "typed_irc_path": {
        "required": True,
        "receipt": "path/receipt.json",
        "requirement": "typed endpoints and one oriented reactant-to-product path",
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
        "requirement": (
            "both competing H2CO channels use one exact common-reactant receipt"
        ),
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
_PATH_TEMPORARY_NAME = re.compile(r"\.path\.[0-9]+\.[0-9]+\.tmp\Z")
_HESSIAN_POINT_TEMPORARY_NAME = re.compile(
    r"\.(?P<index>[0-9]{6})\.[0-9]+\.[0-9]+\.tmp\Z"
)
_HESSIAN_AGGREGATE_TEMPORARY_NAME = re.compile(
    r"\.receipt\.json\.[0-9]+\.[0-9]+\.tmp\Z"
)


@dataclass(frozen=True)
class EndpointBasinClassification:
    """Typed endpoint basin established from a complete covalent graph."""

    route: str
    basin: Literal["reactant", "product"]
    covalent_edges: tuple[tuple[int, int], ...]
    minimum_distance_angstrom: float


def _validate_cluster_state_types(cluster: Cluster, *, label: str) -> None:
    """Reject bool-as-int and malformed index state before route comparisons."""

    if not isinstance(cluster, Cluster):
        raise TypeError(f"{label} must be a Cluster")
    if type(cluster.charge) is not int or type(cluster.spin) is not int:
        raise ValueError(f"{label} electronic state values must be integers")
    if type(cluster.frozen_indices) is not list or any(
        type(index) is not int for index in cluster.frozen_indices
    ):
        raise ValueError(f"{label} frozen indices must be integers")


def endpoint_classification_policy_payload() -> dict[str, Any]:
    """Return the complete JSON-safe endpoint-typing policy bound by preflight."""

    return {
        "collision_floor_angstrom_exclusive_minimum": (
            ENDPOINT_COLLISION_FLOOR_ANGSTROM
        ),
        "boundary_policy": (
            "bonded below bonded_exclusive_maximum; nonbonded above "
            "nonbonded_exclusive_minimum; closed interval rejected"
        ),
        "pair_distances_angstrom": {
            "-".join(pair): {
                "bonded_exclusive_maximum": thresholds[0],
                "nonbonded_exclusive_minimum": thresholds[1],
            }
            for pair, thresholds in sorted(ENDPOINT_PAIR_DISTANCE_POLICY.items())
        },
        "route_states": {
            route: {
                "symbols": list(state["symbols"]),
                "charge": state["charge"],
                "spin_2s": state["spin"],
                "frozen_indices": list(state["frozen_indices"]),
            }
            for route, state in ENDPOINT_ROUTE_STATES.items()
        },
        "exact_covalent_graphs": {
            route: {
                basin: [list(edge) for edge in graphs[basin]]
                for basin in ("reactant", "product")
            }
            for route, graphs in ENDPOINT_GRAPH_CONTRACTS.items()
        },
        "spectator_water_is_part_of_exact_graph": True,
        "sella_direction_is_not_a_chemistry_label": True,
    }


def classify_endpoint_basin(
    route: str, cluster: Cluster
) -> EndpointBasinClassification:
    """Classify one route endpoint from every pair distance and an exact graph.

    Ordered atom identity, electronic state, and frozen indices are route contracts.
    Every pair is typed as covalent or non-covalent; any gray-zone pair blocks the
    classification. This intentionally cannot infer chemistry from Sella direction or
    from whichever single reactive distance happens to be shortest.
    """

    if type(route) is not str or route not in ENDPOINT_ROUTE_STATES:
        raise ValueError(f"unsupported typed endpoint route: {route}")
    _validate_cluster_state_types(cluster, label="endpoint")
    state = ENDPOINT_ROUTE_STATES[route]
    contracts = ENDPOINT_GRAPH_CONTRACTS.get(route)
    if contracts is None:
        raise ValueError(f"unsupported typed endpoint route: {route}")
    if tuple(cluster.symbols) != state["symbols"]:
        raise ValueError(
            f"endpoint ordered symbols do not match route contract: {route}"
        )
    if cluster.charge != state["charge"] or cluster.spin != state["spin"]:
        raise ValueError(
            f"endpoint electronic state does not match route contract: {route}"
        )
    if tuple(cluster.frozen_indices) != state["frozen_indices"]:
        raise ValueError(
            f"endpoint frozen indices do not match route contract: {route}"
        )

    coordinates = np.asarray(cluster.coords, dtype=float)
    expected_shape = (len(state["symbols"]), 3)
    if coordinates.shape != expected_shape or not np.all(np.isfinite(coordinates)):
        raise ValueError(
            f"endpoint coordinates must be a finite {expected_shape} array: {route}"
        )
    covalent_edges: list[tuple[int, int]] = []
    minimum_distance = math.inf
    with np.errstate(over="ignore", invalid="ignore"):
        for first in range(len(cluster.symbols) - 1):
            for second in range(first + 1, len(cluster.symbols)):
                distance = float(
                    np.linalg.norm(coordinates[first] - coordinates[second])
                )
                if not math.isfinite(distance):
                    raise ValueError(f"endpoint pair distances must be finite: {route}")
                minimum_distance = min(minimum_distance, distance)
                if distance <= ENDPOINT_COLLISION_FLOOR_ANGSTROM:
                    raise ValueError(
                        "endpoint violates collision floor for atoms "
                        f"{first}-{second}: {distance:.12g} A"
                    )
                symbols = (cluster.symbols[first], cluster.symbols[second])
                pair = symbols if symbols[0] <= symbols[1] else (symbols[1], symbols[0])
                thresholds = ENDPOINT_PAIR_DISTANCE_POLICY.get(pair)
                if thresholds is None:
                    raise ValueError(
                        f"endpoint pair {pair} has no covalent-distance assignment"
                    )
                bonded_maximum, nonbonded_minimum = thresholds
                if distance < bonded_maximum:
                    covalent_edges.append((first, second))
                elif distance <= nonbonded_minimum:
                    raise ValueError(
                        f"endpoint gray-zone pair {first}-{second} "
                        f"({pair[0]}-{pair[1]}) at {distance:.12g} A"
                    )

    graph = tuple(covalent_edges)
    matches: list[Literal["reactant", "product"]] = [
        basin for basin, expected in contracts.items() if graph == expected
    ]
    if len(matches) > 1:
        raise ValueError(
            f"ambiguous endpoint covalent graph for route {route}: {graph}"
        )
    if not matches:
        raise ValueError(
            f"unassigned endpoint covalent graph for route {route}: {graph}"
        )
    return EndpointBasinClassification(
        route=route,
        basin=matches[0],
        covalent_edges=graph,
        minimum_distance_angstrom=minimum_distance,
    )


@dataclass(frozen=True)
class OrientedIrcPath:
    """One validated reactant-to-TS-to-product path with source provenance."""

    route: str
    coordinates_angstrom: np.ndarray
    electronic_energy_ev: np.ndarray
    projected_fmax_ev_per_angstrom: np.ndarray
    masses_amu: np.ndarray
    transition_state_index: int
    point_provenance: tuple[dict[str, Any], ...]
    endpoint_classifications: tuple[
        EndpointBasinClassification, EndpointBasinClassification
    ]
    mass_scaled_path: MassScaledPath


@dataclass(frozen=True)
class PublishedTypedIrcPath:
    """A validated canonical typed-path publication."""

    receipt_path: Path
    receipt_sha256: str
    receipt: dict[str, Any]
    coordinates_angstrom: np.ndarray
    masses_amu: np.ndarray
    transition_state_index: int
    point_provenance: tuple[dict[str, Any], ...]
    mass_scaled_path: MassScaledPath


def _immutable_little_f64(
    values: Any, shape: tuple[int, ...] | None = None
) -> np.ndarray:
    array = np.ascontiguousarray(values, dtype="<f8")
    if shape is not None and array.shape != shape:
        raise ValueError(f"float64 payload must have shape {shape}, got {array.shape}")
    if not np.all(np.isfinite(array)):
        raise ValueError("float64 payload must be finite")
    return np.frombuffer(array.tobytes(), dtype="<f8").reshape(array.shape)


def _route_identity_cluster(route: str, cluster: Cluster, *, label: str) -> None:
    if type(route) is not str or route not in ENDPOINT_ROUTE_STATES:
        raise ValueError(f"unsupported typed endpoint route: {route}")
    _validate_cluster_state_types(cluster, label=label)
    state = ENDPOINT_ROUTE_STATES[route]
    if tuple(cluster.symbols) != state["symbols"]:
        raise ValueError(f"{label} ordered symbols do not match route contract")
    if cluster.charge != state["charge"] or cluster.spin != state["spin"]:
        raise ValueError(f"{label} electronic state does not match route contract")
    if tuple(cluster.frozen_indices) != state["frozen_indices"]:
        raise ValueError(f"{label} frozen indices do not match route contract")
    coordinates = np.asarray(cluster.coords, dtype=float)
    if coordinates.shape != (len(state["symbols"]), 3) or not np.all(
        np.isfinite(coordinates)
    ):
        raise ValueError(f"{label} coordinates must be finite with route atom shape")


def _cluster_at_coordinates(
    template: Cluster, coordinates: Any, *, name: str
) -> Cluster:
    return Cluster(
        name=name,
        symbols=list(template.symbols),
        coords=np.asarray(coordinates, dtype=float),
        charge=template.charge,
        spin=template.spin,
        frozen_indices=list(template.frozen_indices),
        site_family=template.site_family,
        note=template.note,
    )


def orient_sella_trace(
    route: str,
    qualified_transition_state: Cluster,
    trace: SellaIrcTrace,
) -> OrientedIrcPath:
    """Orient a two-direction Sella trace using typed terminal chemistry only."""

    _route_identity_cluster(route, qualified_transition_state, label="qualified TS")
    if not isinstance(trace, SellaIrcTrace):
        raise TypeError("trace must be a SellaIrcTrace")
    if len(trace.directions) != 2:
        raise ValueError("Sella trace must contain exactly two directions")
    expected_shape = np.asarray(qualified_transition_state.coords).shape
    for direction_index, direction in enumerate(trace.directions):
        expected_direction, expected_sign = (
            ("forward", 1) if direction_index == 0 else ("reverse", -1)
        )
        if type(direction.sella_direction) is not str:
            raise ValueError("source Sella direction must be a string")
        if direction.sella_direction != expected_direction:
            raise ValueError("source Sella directions must be forward then reverse")
        if type(direction.algebraic_direction) is not int:
            raise ValueError("source Sella algebraic direction must be an integer")
        if direction.algebraic_direction != expected_sign:
            raise ValueError("source Sella direction and algebraic direction disagree")
        for expected_outer_step, point in enumerate(direction.points):
            if not isinstance(point, IrcPoint):
                raise TypeError("source IRC path entries must be IrcPoint values")
            if type(point.outer_step) is not int:
                raise ValueError("source IRC outer_step must be an integer")
            if point.outer_step != expected_outer_step:
                raise ValueError("source IRC outer_step sequence is invalid")
            point_coordinates = np.asarray(point.coordinates_angstrom)
            if point_coordinates.shape != expected_shape or not np.all(
                np.isfinite(point_coordinates)
            ):
                raise ValueError("source IRC point coordinates are invalid")
            if type(point.electronic_energy_ev) is not float or not math.isfinite(
                point.electronic_energy_ev
            ):
                raise ValueError("source IRC electronic energy must be a finite float")
            if (
                type(point.projected_fmax_ev_per_angstrom) is not float
                or not math.isfinite(point.projected_fmax_ev_per_angstrom)
                or point.projected_fmax_ev_per_angstrom < 0.0
            ):
                raise ValueError("source IRC projected fmax must be a finite float")
    expected_masses = np.asarray(
        [ISOTOPIC_MASSES_AMU[symbol] for symbol in qualified_transition_state.symbols]
    )
    if trace.masses_amu.shape != expected_masses.shape or not np.array_equal(
        trace.masses_amu, expected_masses
    ):
        raise ValueError("Sella trace masses do not match the route isotopic standard")

    starts = [path.points[0] for path in trace.directions]
    expected_ts = np.asarray(qualified_transition_state.coords, dtype=float)
    if any(
        point.coordinates_angstrom.shape != expected_ts.shape
        or not np.array_equal(point.coordinates_angstrom, expected_ts)
        for point in starts
    ):
        raise ValueError(
            "both Sella TS copies must exactly agree with the qualified TS"
        )
    if (
        starts[0].electronic_energy_ev != starts[1].electronic_energy_ev
        or starts[0].projected_fmax_ev_per_angstrom
        != starts[1].projected_fmax_ev_per_angstrom
    ):
        raise ValueError("both Sella TS copies must have identical point evidence")

    classified: list[tuple[Any, EndpointBasinClassification]] = []
    for direction in trace.directions:
        terminal = _cluster_at_coordinates(
            qualified_transition_state,
            direction.points[-1].coordinates_angstrom,
            name=f"{route}-{direction.sella_direction}-terminal",
        )
        classified.append((direction, classify_endpoint_basin(route, terminal)))
    basins = [classification.basin for _, classification in classified]
    if sorted(basins) != ["product", "reactant"]:
        raise ValueError(
            "Sella terminals must classify as exactly one reactant and one product"
        )
    by_basin = {
        classification.basin: direction for direction, classification in classified
    }

    reactant_direction = by_basin["reactant"]
    product_direction = by_basin["product"]
    reactant_points = list(reversed(reactant_direction.points))
    product_points = list(product_direction.points[1:])
    ordered_points = [*reactant_points, *product_points]
    transition_state_index = len(reactant_points) - 1
    coordinates = _immutable_little_f64(
        [point.coordinates_angstrom for point in ordered_points]
    )
    energies = _immutable_little_f64(
        [point.electronic_energy_ev for point in ordered_points]
    )
    fmax = _immutable_little_f64(
        [point.projected_fmax_ev_per_angstrom for point in ordered_points]
    )
    provenance: list[dict[str, Any]] = []
    for point in reactant_points[:-1]:
        provenance.append(
            {
                "source_sella_direction": reactant_direction.sella_direction,
                "source_outer_step": point.outer_step,
            }
        )
    provenance.append(
        {
            "source_sella_directions": ["forward", "reverse"],
            "source_outer_steps": [0, 0],
            "transition_state": True,
        }
    )
    for point in product_points:
        provenance.append(
            {
                "source_sella_direction": product_direction.sella_direction,
                "source_outer_step": point.outer_step,
            }
        )
    mass_scaled_path = build_mass_scaled_path(
        coordinates,
        expected_masses,
        transition_state_index=transition_state_index,
        reference_mass_amu=REFERENCE_MASS_AMU,
    )
    return OrientedIrcPath(
        route=route,
        coordinates_angstrom=coordinates,
        electronic_energy_ev=energies,
        projected_fmax_ev_per_angstrom=fmax,
        masses_amu=_immutable_little_f64(expected_masses),
        transition_state_index=transition_state_index,
        point_provenance=tuple(provenance),
        endpoint_classifications=(
            next(c for _, c in classified if c.basin == "reactant"),
            next(c for _, c in classified if c.basin == "product"),
        ),
        mass_scaled_path=mass_scaled_path,
    )


def _safe_absolute_root(path: Path) -> Path:
    if ".." in path.parts:
        raise ValueError(f"run-root path escape is forbidden: {path}")
    absolute = Path(os.path.abspath(path))
    for candidate in reversed((absolute, *absolute.parents)):
        if candidate.is_symlink():
            raise ValueError(f"run-root path contains a symlink: {candidate}")
    return absolute


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _write_fsync(path: Path, data: bytes) -> None:
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags, 0o600)
    with os.fdopen(descriptor, "wb") as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())


def _strict_json_equal(actual: Any, expected: Any, *, label: str) -> None:
    """Require exact JSON container, primitive type, and value equality."""

    if type(actual) is not type(expected):
        raise ValueError(f"{label} has invalid JSON primitive type")
    if isinstance(expected, dict):
        if set(actual) != set(expected):
            raise ValueError(f"{label} fields are unexpected or incomplete")
        for key, value in expected.items():
            _strict_json_equal(actual[key], value, label=f"{label}.{key}")
        return
    if isinstance(expected, list):
        if len(actual) != len(expected):
            raise ValueError(f"{label} list length mismatch")
        for index, value in enumerate(expected):
            _strict_json_equal(actual[index], value, label=f"{label}[{index}]")
        return
    if isinstance(expected, float) and not math.isfinite(actual):
        raise ValueError(f"{label} must be finite")
    if actual != expected:
        raise ValueError(f"{label} mismatch")


def _require_json_string(value: Any, *, label: str) -> str:
    if type(value) is not str:
        raise ValueError(f"{label} must be a JSON string")
    return value


def _require_json_float(value: Any, *, label: str) -> float:
    if type(value) is not float or not math.isfinite(value):
        raise ValueError(f"{label} must be a finite JSON float")
    return value


def _renameat2_noreplace(source: Path, destination: Path) -> None:
    """Atomically rename without replacement, or report unsupported platforms."""

    if os.name != "posix" or not sys.platform.startswith("linux"):
        raise NotImplementedError("renameat2(RENAME_NOREPLACE) requires Linux")
    library = ctypes.CDLL(None, use_errno=True)
    try:
        renameat2 = library.renameat2
    except AttributeError as exc:
        raise NotImplementedError("libc does not expose renameat2") from exc
    renameat2.argtypes = [
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    ]
    renameat2.restype = ctypes.c_int
    at_fdcwd = -100
    rename_noreplace = 1
    result = renameat2(
        at_fdcwd,
        os.fsencode(source),
        at_fdcwd,
        os.fsencode(destination),
        rename_noreplace,
    )
    if result == 0:
        return
    error_number = ctypes.get_errno()
    if error_number in {errno.ENOSYS, errno.EINVAL, errno.EOPNOTSUPP}:
        raise NotImplementedError(
            "renameat2(RENAME_NOREPLACE) is unavailable on this filesystem"
        )
    if error_number in {errno.EEXIST, errno.ENOTEMPTY}:
        raise FileExistsError(
            error_number,
            os.strerror(error_number),
            str(destination),
        )
    raise OSError(error_number, os.strerror(error_number), str(destination))


def _publish_noreplace(
    source: Path, destination: Path, *, source_is_directory: bool
) -> None:
    """Publish at one no-clobber boundary, failing closed for directories."""

    try:
        _renameat2_noreplace(source, destination)
        return
    except NotImplementedError:
        if source_is_directory:
            raise RuntimeError(
                "atomic directory publication requires renameat2(RENAME_NOREPLACE)"
            ) from None
    if source.is_symlink() or not source.is_file():
        raise ValueError("file publication source must be a regular file")
    try:
        os.link(source, destination, follow_symlinks=False)
    except FileExistsError:
        raise
    except OSError as exc:
        raise RuntimeError(
            "atomic file publication requires renameat2 or same-filesystem hard links"
        ) from exc
    source.unlink()


def _remove_owned_temporary_directories(
    parent: Path,
    *,
    name_pattern: re.Pattern[str],
    allowed_files: set[str],
    maximum_point_index: int | None = None,
) -> None:
    removed = False
    for temporary in list(parent.iterdir()):
        match = name_pattern.fullmatch(temporary.name)
        if match is None:
            continue
        if maximum_point_index is not None:
            point_index = int(match.group("index"))
            if point_index >= maximum_point_index:
                continue
        if (
            temporary.parent != parent
            or temporary.is_symlink()
            or not temporary.is_dir()
        ):
            raise ValueError(
                f"owned temporary must be a real confined directory: {temporary}"
            )
        children = list(temporary.iterdir())
        unexpected = {child.name for child in children} - allowed_files
        if unexpected or any(
            child.is_symlink() or not child.is_file() for child in children
        ):
            raise ValueError(f"owned temporary has invalid contents: {temporary}")
        shutil.rmtree(temporary)
        removed = True
    if removed:
        _fsync_directory(parent)


def _remove_owned_temporary_files(
    parent: Path, *, name_pattern: re.Pattern[str]
) -> None:
    removed = False
    for temporary in list(parent.iterdir()):
        if name_pattern.fullmatch(temporary.name) is None:
            continue
        if (
            temporary.parent != parent
            or temporary.is_symlink()
            or not temporary.is_file()
        ):
            raise ValueError(
                f"owned temporary must be a real confined file: {temporary}"
            )
        temporary.unlink()
        removed = True
    if removed:
        _fsync_directory(parent)


def _json_bytes(payload: dict[str, Any]) -> bytes:
    return (
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode()


@contextmanager
def _exclusive_route_claim(route_root: Path) -> Iterator[None]:
    route_existed = route_root.exists()
    route_root.mkdir(parents=True, exist_ok=True)
    if route_root.is_symlink() or not route_root.is_dir():
        raise ValueError(f"route root must be a real directory: {route_root}")
    if not route_existed:
        _fsync_directory(route_root.parent)
    lock_path = route_root / ".route.lock"
    if lock_path.is_symlink():
        raise ValueError(f"route lock must not be a symlink: {lock_path}")
    descriptor = os.open(
        lock_path,
        os.O_RDWR | os.O_CREAT | getattr(os, "O_NOFOLLOW", 0),
        0o600,
    )
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        yield
    finally:
        fcntl.flock(descriptor, fcntl.LOCK_UN)
        os.close(descriptor)


def _path_receipt_payload(
    campaign_identity: str,
    route: str,
    atom_mapping_sha256: str,
    qualified_transition_state: Cluster,
    oriented: OrientedIrcPath,
) -> tuple[dict[str, Any], bytes]:
    coordinates_bytes = oriented.coordinates_angstrom.tobytes()
    points = []
    for index, provenance in enumerate(oriented.point_provenance):
        point_bytes = oriented.coordinates_angstrom[index].tobytes()
        points.append(
            {
                "index": index,
                "geometry_sha256": hashlib.sha256(point_bytes).hexdigest(),
                "electronic_energy_ev": float(oriented.electronic_energy_ev[index]),
                "projected_fmax_ev_per_angstrom": float(
                    oriented.projected_fmax_ev_per_angstrom[index]
                ),
                "source": provenance,
            }
        )
    receipt = {
        "schema": "d2c-typed-irc-path-v1",
        "stage": "typed_irc_path",
        "state": "accepted",
        "accepted": True,
        "campaign_identity": campaign_identity,
        "route": route,
        "atom_mapping_sha256": atom_mapping_sha256,
        "qualified_transition_state_geometry_fingerprint": (
            frequency_geometry_fingerprint(qualified_transition_state)
        ),
        "endpoint_classification_policy_sha256": _canonical_hash(
            endpoint_classification_policy_payload()
        ),
        "symbols": list(qualified_transition_state.symbols),
        "charge": qualified_transition_state.charge,
        "spin_2s": qualified_transition_state.spin,
        "frozen_indices": list(qualified_transition_state.frozen_indices),
        "masses_amu": oriented.masses_amu.tolist(),
        "point_count": len(points),
        "transition_state_index": oriented.transition_state_index,
        "coordinates": {
            "path": "coordinates.f64",
            "dtype": "little-endian float64",
            "shape": list(oriented.coordinates_angstrom.shape),
            "units": "angstrom",
            "sha256": hashlib.sha256(coordinates_bytes).hexdigest(),
        },
        "endpoints": {
            classification.basin: {
                "point_index": (
                    0 if classification.basin == "reactant" else len(points) - 1
                ),
                "covalent_edges": [
                    list(edge) for edge in classification.covalent_edges
                ],
                "minimum_distance_angstrom": (classification.minimum_distance_angstrom),
            }
            for classification in oriented.endpoint_classifications
        },
        "points": points,
    }
    return receipt, coordinates_bytes


def _read_json_object(path: Path, *, label: str) -> tuple[dict[str, Any], bytes]:
    if path.is_symlink() or not path.is_file():
        raise ValueError(f"{label} must be a regular file: {path}")
    raw = path.read_bytes()
    try:
        payload = json.loads(raw)
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise ValueError(f"malformed {label}: {path}") from exc
    if not isinstance(payload, dict) or _json_bytes(payload) != raw:
        raise ValueError(f"{label} is not canonical JSON: {path}")
    return payload, raw


def _validate_published_path(
    path_dir: Path,
    *,
    campaign_identity: str,
    route: str,
    atom_mapping_sha256: str,
    qualified_transition_state: Cluster,
) -> PublishedTypedIrcPath:
    if path_dir.is_symlink() or not path_dir.is_dir():
        raise ValueError(
            f"canonical path publication must be a real directory: {path_dir}"
        )
    allowed = {"coordinates.f64", "receipt.json"}
    observed = {child.name for child in path_dir.iterdir()}
    unexpected = observed - allowed
    if unexpected:
        raise ValueError(f"unexpected path artifact(s): {sorted(unexpected)}")
    if observed != allowed:
        raise ValueError("canonical path publication is incomplete")
    receipt, receipt_bytes = _read_json_object(
        path_dir / "receipt.json", label="typed path receipt"
    )
    expected_receipt_keys = {
        "schema",
        "stage",
        "state",
        "accepted",
        "campaign_identity",
        "route",
        "atom_mapping_sha256",
        "qualified_transition_state_geometry_fingerprint",
        "endpoint_classification_policy_sha256",
        "symbols",
        "charge",
        "spin_2s",
        "frozen_indices",
        "masses_amu",
        "point_count",
        "transition_state_index",
        "coordinates",
        "endpoints",
        "points",
    }
    if set(receipt) != expected_receipt_keys:
        raise ValueError("typed path receipt fields are unexpected or incomplete")
    expected_identity = {
        "schema": "d2c-typed-irc-path-v1",
        "stage": "typed_irc_path",
        "state": "accepted",
        "accepted": True,
        "campaign_identity": campaign_identity,
        "route": route,
        "atom_mapping_sha256": atom_mapping_sha256,
        "qualified_transition_state_geometry_fingerprint": (
            frequency_geometry_fingerprint(qualified_transition_state)
        ),
        "endpoint_classification_policy_sha256": _canonical_hash(
            endpoint_classification_policy_payload()
        ),
        "symbols": list(qualified_transition_state.symbols),
        "charge": qualified_transition_state.charge,
        "spin_2s": qualified_transition_state.spin,
        "frozen_indices": list(qualified_transition_state.frozen_indices),
    }
    for key, value in expected_identity.items():
        _strict_json_equal(receipt.get(key), value, label=f"typed path receipt {key}")
    point_count = receipt.get("point_count")
    ts_index = receipt.get("transition_state_index")
    coordinate_record = receipt.get("coordinates")
    if not isinstance(coordinate_record, dict):
        raise ValueError("typed path coordinate metadata is invalid")
    if (
        type(point_count) is not int
        or point_count < 3
        or type(ts_index) is not int
        or not 0 < ts_index < point_count - 1
        or type(receipt.get("points")) is not list
        or len(receipt["points"]) != point_count
    ):
        raise ValueError("typed path receipt point layout is invalid")
    expected_shape = [point_count, len(qualified_transition_state.symbols), 3]
    coordinate_sha = _require_json_string(
        coordinate_record.get("sha256"), label="typed path coordinate SHA-256"
    )
    _require_sha(coordinate_sha, length=64, label="typed path coordinate SHA-256")
    _strict_json_equal(
        coordinate_record,
        {
            "path": "coordinates.f64",
            "dtype": "little-endian float64",
            "shape": expected_shape,
            "units": "angstrom",
            "sha256": coordinate_sha,
        },
        label="typed path coordinate metadata",
    )
    coordinates_path = path_dir / "coordinates.f64"
    if coordinates_path.is_symlink() or not coordinates_path.is_file():
        raise ValueError("typed path coordinates must be a regular file")
    raw_coordinates = coordinates_path.read_bytes()
    if hashlib.sha256(raw_coordinates).hexdigest() != coordinate_record["sha256"]:
        raise ValueError("typed path coordinates hash mismatch")
    expected_bytes = point_count * len(qualified_transition_state.symbols) * 3 * 8
    if len(raw_coordinates) != expected_bytes:
        raise ValueError("typed path coordinate byte length mismatch")
    coordinates = np.frombuffer(raw_coordinates, dtype="<f8").reshape(expected_shape)
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("typed path coordinates contain non-finite values")
    if not np.array_equal(coordinates[ts_index], qualified_transition_state.coords):
        raise ValueError("typed path TS does not match qualified TS")
    if (
        sum(
            np.array_equal(point, qualified_transition_state.coords)
            for point in coordinates
        )
        != 1
    ):
        raise ValueError("typed path must contain exactly one qualified TS copy")
    expected_masses = np.asarray(
        [ISOTOPIC_MASSES_AMU[symbol] for symbol in qualified_transition_state.symbols]
    )
    _strict_json_equal(
        receipt.get("masses_amu"),
        expected_masses.tolist(),
        label="typed path masses",
    )
    masses = _immutable_little_f64(receipt["masses_amu"])
    for index, point in enumerate(receipt["points"]):
        if type(point) is not dict or set(point) != {
            "index",
            "geometry_sha256",
            "electronic_energy_ev",
            "projected_fmax_ev_per_angstrom",
            "source",
        }:
            raise ValueError("typed path point fields are unexpected or incomplete")
        _strict_json_equal(
            point.get("index"), index, label=f"typed path point {index} index"
        )
        geometry_sha = _require_json_string(
            point.get("geometry_sha256"),
            label=f"typed path point {index} geometry SHA-256",
        )
        _require_sha(
            geometry_sha,
            length=64,
            label=f"typed path point {index} geometry SHA-256",
        )
        if geometry_sha != hashlib.sha256(coordinates[index].tobytes()).hexdigest():
            raise ValueError("typed path point geometry hash mismatch")
        _require_json_float(
            point.get("electronic_energy_ev"),
            label=f"typed path point {index} electronic energy",
        )
        fmax = _require_json_float(
            point.get("projected_fmax_ev_per_angstrom"),
            label=f"typed path point {index} projected fmax",
        )
        if fmax < 0.0 or type(point.get("source")) is not dict:
            raise ValueError("typed path point evidence is invalid")
    expected_ts_source = {
        "source_sella_directions": ["forward", "reverse"],
        "source_outer_steps": [0, 0],
        "transition_state": True,
    }
    _strict_json_equal(
        receipt["points"][ts_index]["source"],
        expected_ts_source,
        label="typed path TS source provenance",
    )
    left_direction = _require_json_string(
        receipt["points"][0]["source"].get("source_sella_direction"),
        label="typed path reactant endpoint source direction",
    )
    right_direction = _require_json_string(
        receipt["points"][-1]["source"].get("source_sella_direction"),
        label="typed path product endpoint source direction",
    )
    if {left_direction, right_direction} != {"forward", "reverse"}:
        raise ValueError("typed path endpoint source directions are invalid")
    for index in range(ts_index):
        _strict_json_equal(
            receipt["points"][index]["source"],
            {
                "source_sella_direction": left_direction,
                "source_outer_step": ts_index - index,
            },
            label="typed path reactant-side provenance",
        )
    for index in range(ts_index + 1, point_count):
        _strict_json_equal(
            receipt["points"][index]["source"],
            {
                "source_sella_direction": right_direction,
                "source_outer_step": index - ts_index,
            },
            label="typed path product-side provenance",
        )
    reactant = classify_endpoint_basin(
        route,
        _cluster_at_coordinates(
            qualified_transition_state,
            coordinates[0],
            name=f"{route}-published-reactant",
        ),
    )
    product = classify_endpoint_basin(
        route,
        _cluster_at_coordinates(
            qualified_transition_state,
            coordinates[-1],
            name=f"{route}-published-product",
        ),
    )
    if reactant.basin != "reactant" or product.basin != "product":
        raise ValueError("typed path endpoint orientation is invalid")
    expected_endpoints = {
        "reactant": {
            "point_index": 0,
            "covalent_edges": [list(edge) for edge in reactant.covalent_edges],
            "minimum_distance_angstrom": reactant.minimum_distance_angstrom,
        },
        "product": {
            "point_index": point_count - 1,
            "covalent_edges": [list(edge) for edge in product.covalent_edges],
            "minimum_distance_angstrom": product.minimum_distance_angstrom,
        },
    }
    _strict_json_equal(
        receipt.get("endpoints"),
        expected_endpoints,
        label="typed path endpoint evidence",
    )
    mass_scaled_path = build_mass_scaled_path(
        coordinates,
        masses,
        transition_state_index=ts_index,
        reference_mass_amu=REFERENCE_MASS_AMU,
    )
    return PublishedTypedIrcPath(
        receipt_path=path_dir / "receipt.json",
        receipt_sha256=hashlib.sha256(receipt_bytes).hexdigest(),
        receipt=receipt,
        coordinates_angstrom=coordinates,
        masses_amu=masses,
        transition_state_index=ts_index,
        point_provenance=tuple(point["source"] for point in receipt["points"]),
        mass_scaled_path=mass_scaled_path,
    )


def publish_typed_irc_path(
    run_root: Path,
    *,
    campaign_identity: str,
    route: str,
    atom_mapping_sha256: str,
    qualified_transition_state: Cluster,
    trace: SellaIrcTrace | None,
    _failure_injector: Callable[[str], None] | None = None,
) -> PublishedTypedIrcPath:
    """Atomically publish or validate/resume one typed path stage."""

    _require_sha(campaign_identity, length=64, label="campaign identity")
    _require_sha(atom_mapping_sha256, length=64, label="atom mapping SHA-256")
    _route_identity_cluster(route, qualified_transition_state, label="qualified TS")
    root = _safe_absolute_root(run_root)
    root.mkdir(parents=True, exist_ok=True)
    if root.is_symlink() or not root.is_dir():
        raise ValueError(f"run root must be a real directory: {root}")
    route_root = root / route
    with _exclusive_route_claim(route_root):
        _remove_owned_temporary_directories(
            route_root,
            name_pattern=_PATH_TEMPORARY_NAME,
            allowed_files={"coordinates.f64", "receipt.json"},
        )
        path_dir = route_root / "path"
        if path_dir.exists() or path_dir.is_symlink():
            return _validate_published_path(
                path_dir,
                campaign_identity=campaign_identity,
                route=route,
                atom_mapping_sha256=atom_mapping_sha256,
                qualified_transition_state=qualified_transition_state,
            )
        if trace is None:
            raise ValueError(
                "a Sella trace is required for a new typed path publication"
            )
        oriented = orient_sella_trace(route, qualified_transition_state, trace)
        receipt, coordinate_bytes = _path_receipt_payload(
            campaign_identity,
            route,
            atom_mapping_sha256,
            qualified_transition_state,
            oriented,
        )
        temporary = route_root / f".path.{os.getpid()}.{time.time_ns()}.tmp"
        temporary.mkdir(mode=0o700)
        try:
            _write_fsync(temporary / "coordinates.f64", coordinate_bytes)
            _write_fsync(temporary / "receipt.json", _json_bytes(receipt))
            _fsync_directory(temporary)
            _validate_published_path(
                temporary,
                campaign_identity=campaign_identity,
                route=route,
                atom_mapping_sha256=atom_mapping_sha256,
                qualified_transition_state=qualified_transition_state,
            )
            if _failure_injector is not None:
                _failure_injector("before_path_commit")
            if path_dir.exists() or path_dir.is_symlink():
                raise FileExistsError(
                    f"typed path publication already exists: {path_dir}"
                )
            _publish_noreplace(temporary, path_dir, source_is_directory=True)
            _fsync_directory(route_root)
        finally:
            if temporary.exists():
                shutil.rmtree(temporary)
        return _validate_published_path(
            path_dir,
            campaign_identity=campaign_identity,
            route=route,
            atom_mapping_sha256=atom_mapping_sha256,
            qualified_transition_state=qualified_transition_state,
        )


@dataclass(frozen=True)
class PublishedPathHessians:
    """A fully validated aggregate of immutable per-path-point Hessians."""

    receipt_path: Path
    receipt_sha256: str
    receipt: dict[str, Any]
    point_results: tuple[NativeHessianResult, ...]


def _validated_backend_policy(policy: dict[str, Any]) -> dict[str, Any]:
    if set(policy) != {"requested_backend", "allow_cpu_fallback"}:
        raise ValueError("native Hessian backend policy must have exact fields")
    requested = policy.get("requested_backend")
    allow_fallback = policy.get("allow_cpu_fallback")
    if requested not in {"pyscf", "gpu4pyscf"} or not isinstance(allow_fallback, bool):
        raise ValueError("native Hessian backend policy is invalid")
    if requested == "pyscf" and allow_fallback:
        raise ValueError("CPU native Hessian policy cannot allow GPU fallback")
    return {
        "requested_backend": requested,
        "allow_cpu_fallback": allow_fallback,
    }


def _validate_evaluated_native_hessian(
    result: NativeHessianResult,
    cluster: Cluster,
    *,
    settings_fingerprint: str,
    backend_policy: dict[str, Any],
) -> NativeHessianResult:
    if not isinstance(result, NativeHessianResult):
        raise TypeError("native Hessian evaluator must return NativeHessianResult")
    if type(result.gpu_fallback_used) is not bool:
        raise ValueError("native Hessian gpu_fallback_used must be a boolean")
    # Reconstruct to re-run every invariant even if a fake or deserializer bypassed
    # the frozen dataclass's constructor.
    checked = NativeHessianResult(
        electronic_hartree=result.electronic_hartree,
        gradient_hartree_per_bohr=result.gradient_hartree_per_bohr,
        physical_fmax_ev_per_angstrom=result.physical_fmax_ev_per_angstrom,
        cartesian_hessian_hartree_per_bohr2=(
            result.cartesian_hessian_hartree_per_bohr2
        ),
        requested_backend=result.requested_backend,
        actual_backend=result.actual_backend,
        gpu_fallback_used=result.gpu_fallback_used,
        geometry_fingerprint=result.geometry_fingerprint,
        settings_fingerprint=result.settings_fingerprint,
    )
    expected_geometry = frequency_geometry_fingerprint(cluster)
    if checked.geometry_fingerprint != expected_geometry:
        raise ValueError("native Hessian geometry fingerprint mismatch")
    if checked.settings_fingerprint != settings_fingerprint:
        raise ValueError("native Hessian settings fingerprint mismatch")
    requested = backend_policy["requested_backend"]
    if checked.requested_backend != requested:
        raise ValueError("native Hessian result violates backend policy")
    if checked.gpu_fallback_used and not backend_policy["allow_cpu_fallback"]:
        raise ValueError("native Hessian result violates backend policy fallback")
    expected_actual = "pyscf" if checked.gpu_fallback_used else requested
    if checked.actual_backend != expected_actual:
        raise ValueError("native Hessian result violates backend policy provenance")
    return checked


def _point_hessian_receipt(
    *,
    campaign_identity: str,
    route: str,
    atom_mapping_sha256: str,
    qualified_transition_state: Cluster,
    path: PublishedTypedIrcPath,
    point_index: int,
    cluster: Cluster,
    settings_fingerprint: str,
    backend_policy: dict[str, Any],
    result: NativeHessianResult,
) -> tuple[dict[str, Any], bytes, bytes]:
    gradient_bytes = result.gradient_hartree_per_bohr.tobytes()
    hessian_bytes = result.cartesian_hessian_hartree_per_bohr2.tobytes()
    point_coordinates = path.coordinates_angstrom[point_index]
    receipt = {
        "schema": "d2c-native-path-hessian-point-v1",
        "stage": "hessian_every_path_point",
        "state": "accepted",
        "accepted": True,
        "campaign_identity": campaign_identity,
        "route": route,
        "path_receipt_sha256": path.receipt_sha256,
        "atom_mapping_sha256": atom_mapping_sha256,
        "point_index": point_index,
        "point_count": len(path.coordinates_angstrom),
        "transition_state_index": path.transition_state_index,
        "is_transition_state": point_index == path.transition_state_index,
        "path_point_provenance": path.point_provenance[point_index],
        "symbols": list(qualified_transition_state.symbols),
        "charge": qualified_transition_state.charge,
        "spin_2s": qualified_transition_state.spin,
        "frozen_indices": list(qualified_transition_state.frozen_indices),
        "masses_amu": path.masses_amu.tolist(),
        "geometry_sha256": hashlib.sha256(point_coordinates.tobytes()).hexdigest(),
        "geometry_fingerprint": result.geometry_fingerprint,
        "settings_fingerprint": settings_fingerprint,
        "backend_policy": backend_policy,
        "electronic_hartree": result.electronic_hartree,
        "physical_fmax_ev_per_angstrom": result.physical_fmax_ev_per_angstrom,
        "requested_backend": result.requested_backend,
        "actual_backend": result.actual_backend,
        "gpu_fallback_used": result.gpu_fallback_used,
        "gradient": {
            "path": "gradient.f64",
            "dtype": "little-endian float64",
            "shape": list(result.gradient_hartree_per_bohr.shape),
            "units": "hartree / bohr",
            "sha256": hashlib.sha256(gradient_bytes).hexdigest(),
        },
        "cartesian_hessian": {
            "path": "hessian.f64",
            "dtype": "little-endian float64",
            "shape": list(result.cartesian_hessian_hartree_per_bohr2.shape),
            "units": "hartree / bohr^2",
            "sha256": hashlib.sha256(hessian_bytes).hexdigest(),
        },
    }
    return receipt, gradient_bytes, hessian_bytes


def _validate_hessian_point(
    point_dir: Path,
    *,
    campaign_identity: str,
    route: str,
    atom_mapping_sha256: str,
    qualified_transition_state: Cluster,
    path: PublishedTypedIrcPath,
    point_index: int,
    settings_fingerprint: str,
    backend_policy: dict[str, Any],
) -> tuple[NativeHessianResult, str]:
    if point_dir.is_symlink() or not point_dir.is_dir():
        raise ValueError(f"Hessian point must be a real directory: {point_dir}")
    allowed = {"gradient.f64", "hessian.f64", "receipt.json"}
    observed = {child.name for child in point_dir.iterdir()}
    if observed - allowed:
        raise ValueError(
            f"unexpected Hessian point artifact(s): {sorted(observed - allowed)}"
        )
    if observed != allowed:
        raise ValueError(f"incomplete Hessian point publication: {point_index}")
    receipt, receipt_bytes = _read_json_object(
        point_dir / "receipt.json", label=f"Hessian point {point_index} receipt"
    )
    expected_receipt_keys = {
        "schema",
        "stage",
        "state",
        "accepted",
        "campaign_identity",
        "route",
        "path_receipt_sha256",
        "atom_mapping_sha256",
        "point_index",
        "point_count",
        "transition_state_index",
        "is_transition_state",
        "path_point_provenance",
        "symbols",
        "charge",
        "spin_2s",
        "frozen_indices",
        "masses_amu",
        "geometry_sha256",
        "geometry_fingerprint",
        "settings_fingerprint",
        "backend_policy",
        "electronic_hartree",
        "physical_fmax_ev_per_angstrom",
        "requested_backend",
        "actual_backend",
        "gpu_fallback_used",
        "gradient",
        "cartesian_hessian",
    }
    if set(receipt) != expected_receipt_keys:
        raise ValueError(
            f"Hessian point {point_index} receipt fields are unexpected or incomplete"
        )
    cluster = _cluster_at_coordinates(
        qualified_transition_state,
        path.coordinates_angstrom[point_index],
        name=f"{route}-path-point-{point_index}",
    )
    expected_scalars = {
        "schema": "d2c-native-path-hessian-point-v1",
        "stage": "hessian_every_path_point",
        "state": "accepted",
        "accepted": True,
        "campaign_identity": campaign_identity,
        "route": route,
        "path_receipt_sha256": path.receipt_sha256,
        "atom_mapping_sha256": atom_mapping_sha256,
        "point_index": point_index,
        "point_count": len(path.coordinates_angstrom),
        "transition_state_index": path.transition_state_index,
        "is_transition_state": point_index == path.transition_state_index,
        "path_point_provenance": path.point_provenance[point_index],
        "symbols": list(qualified_transition_state.symbols),
        "charge": qualified_transition_state.charge,
        "spin_2s": qualified_transition_state.spin,
        "frozen_indices": list(qualified_transition_state.frozen_indices),
        "masses_amu": path.masses_amu.tolist(),
        "geometry_sha256": hashlib.sha256(
            path.coordinates_angstrom[point_index].tobytes()
        ).hexdigest(),
        "geometry_fingerprint": frequency_geometry_fingerprint(cluster),
        "settings_fingerprint": settings_fingerprint,
        "backend_policy": backend_policy,
    }
    for key, value in expected_scalars.items():
        _strict_json_equal(
            receipt.get(key), value, label=f"Hessian point {point_index} {key}"
        )
    atom_count = len(qualified_transition_state.symbols)
    records = {
        "gradient": (
            point_dir / "gradient.f64",
            (atom_count, 3),
            "hartree / bohr",
        ),
        "cartesian_hessian": (
            point_dir / "hessian.f64",
            (3 * atom_count, 3 * atom_count),
            "hartree / bohr^2",
        ),
    }
    arrays: dict[str, np.ndarray] = {}
    for label, (data_path, shape, units) in records.items():
        record = receipt.get(label)
        expected_filename = data_path.name
        if type(record) is not dict:
            raise ValueError(f"Hessian point {point_index} {label} metadata mismatch")
        data_sha = _require_json_string(
            record.get("sha256"),
            label=f"Hessian point {point_index} {label} SHA-256",
        )
        _require_sha(
            data_sha,
            length=64,
            label=f"Hessian point {point_index} {label} SHA-256",
        )
        _strict_json_equal(
            record,
            {
                "path": expected_filename,
                "dtype": "little-endian float64",
                "shape": list(shape),
                "units": units,
                "sha256": data_sha,
            },
            label=f"Hessian point {point_index} {label} metadata",
        )
        if data_path.is_symlink() or not data_path.is_file():
            raise ValueError(
                f"Hessian point {point_index} {label} must be a regular file"
            )
        data = data_path.read_bytes()
        if hashlib.sha256(data).hexdigest() != record["sha256"]:
            short_label = "gradient" if label == "gradient" else "Hessian"
            raise ValueError(f"Hessian point {point_index} {short_label} hash mismatch")
        if len(data) != math.prod(shape) * 8:
            raise ValueError(
                f"Hessian point {point_index} {label} byte length mismatch"
            )
        arrays[label] = np.frombuffer(data, dtype="<f8").reshape(shape)
    required_result_fields = {
        "electronic_hartree",
        "physical_fmax_ev_per_angstrom",
        "requested_backend",
        "actual_backend",
        "gpu_fallback_used",
    }
    if any(field not in receipt for field in required_result_fields):
        raise ValueError(f"Hessian point {point_index} result provenance is incomplete")
    electronic_hartree = _require_json_float(
        receipt["electronic_hartree"],
        label=f"Hessian point {point_index} electronic_hartree",
    )
    physical_fmax = _require_json_float(
        receipt["physical_fmax_ev_per_angstrom"],
        label=f"Hessian point {point_index} physical_fmax_ev_per_angstrom",
    )
    requested_backend = _require_json_string(
        receipt["requested_backend"],
        label=f"Hessian point {point_index} requested_backend",
    )
    actual_backend = _require_json_string(
        receipt["actual_backend"],
        label=f"Hessian point {point_index} actual_backend",
    )
    gpu_fallback_used = receipt["gpu_fallback_used"]
    if type(gpu_fallback_used) is not bool:
        raise ValueError(
            f"Hessian point {point_index} gpu_fallback_used must be a JSON boolean"
        )
    result = NativeHessianResult(
        electronic_hartree=electronic_hartree,
        gradient_hartree_per_bohr=arrays["gradient"],
        physical_fmax_ev_per_angstrom=physical_fmax,
        cartesian_hessian_hartree_per_bohr2=arrays["cartesian_hessian"],
        requested_backend=requested_backend,
        actual_backend=actual_backend,
        gpu_fallback_used=gpu_fallback_used,
        geometry_fingerprint=receipt["geometry_fingerprint"],
        settings_fingerprint=receipt["settings_fingerprint"],
    )
    result = _validate_evaluated_native_hessian(
        result,
        cluster,
        settings_fingerprint=settings_fingerprint,
        backend_policy=backend_policy,
    )
    return result, hashlib.sha256(receipt_bytes).hexdigest()


def _child_hash_chain(child_hashes: list[str]) -> list[dict[str, Any]]:
    previous = bytes(32)
    chain: list[dict[str, Any]] = []
    for index, child_hash in enumerate(child_hashes):
        _require_sha(child_hash, length=64, label=f"child receipt {index} SHA-256")
        chained = hashlib.sha256(previous + bytes.fromhex(child_hash)).digest()
        chain.append(
            {
                "point_index": index,
                "receipt_path": f"points/{index:06d}/receipt.json",
                "receipt_sha256": child_hash,
                "chain_sha256": chained.hex(),
            }
        )
        previous = chained
    return chain


def _aggregate_hessian_receipt(
    *,
    campaign_identity: str,
    route: str,
    atom_mapping_sha256: str,
    path: PublishedTypedIrcPath,
    settings_fingerprint: str,
    backend_policy: dict[str, Any],
    child_hashes: list[str],
) -> dict[str, Any]:
    chain = _child_hash_chain(child_hashes)
    return {
        "schema": "d2c-native-path-hessians-v1",
        "stage": "hessian_every_path_point",
        "state": "accepted",
        "accepted": True,
        "campaign_identity": campaign_identity,
        "route": route,
        "path_receipt_sha256": path.receipt_sha256,
        "atom_mapping_sha256": atom_mapping_sha256,
        "settings_fingerprint": settings_fingerprint,
        "backend_policy": backend_policy,
        "point_count": len(path.coordinates_angstrom),
        "transition_state_index": path.transition_state_index,
        "ordered_child_hash_chain": chain,
        "final_child_chain_sha256": chain[-1]["chain_sha256"],
    }


def _validate_hessian_tree_shape(
    hessian_root: Path,
    points_root: Path,
    point_count: int,
    *,
    staged_receipt_name: str | None = None,
) -> None:
    if hessian_root.is_symlink() or not hessian_root.is_dir():
        raise ValueError("Hessian stage root must be a real directory")
    root_names = {child.name for child in hessian_root.iterdir()}
    allowed_root_names = {"points", "receipt.json"}
    if staged_receipt_name is not None:
        allowed_root_names.add(staged_receipt_name)
    unexpected_root = root_names - allowed_root_names
    if unexpected_root:
        raise ValueError(
            f"unexpected Hessian stage artifact(s): {sorted(unexpected_root)}"
        )
    if points_root.is_symlink() or not points_root.is_dir():
        raise ValueError("Hessian points root must be a real directory")
    expected_names = {f"{index:06d}" for index in range(point_count)}
    point_names = {child.name for child in points_root.iterdir()}
    unexpected_points = point_names - expected_names
    if unexpected_points:
        raise ValueError(
            f"unexpected Hessian point artifact(s): {sorted(unexpected_points)}"
        )


def _atomic_non_overwriting_file(
    path: Path,
    data: bytes,
    *,
    validate_before_commit: Callable[[Path], None] | None = None,
) -> None:
    if path.exists() or path.is_symlink():
        raise FileExistsError(f"canonical receipt already exists: {path}")
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    try:
        _write_fsync(temporary, data)
        if validate_before_commit is not None:
            validate_before_commit(temporary)
        if path.exists() or path.is_symlink():
            raise FileExistsError(f"canonical receipt already exists: {path}")
        _publish_noreplace(temporary, path, source_is_directory=False)
        _fsync_directory(path.parent)
    finally:
        temporary.unlink(missing_ok=True)


def _validate_hessian_aggregate(
    hessian_root: Path,
    *,
    campaign_identity: str,
    route: str,
    atom_mapping_sha256: str,
    qualified_transition_state: Cluster,
    path: PublishedTypedIrcPath,
    settings_fingerprint: str,
    backend_policy: dict[str, Any],
    receipt_path: Path | None = None,
) -> PublishedPathHessians:
    canonical_receipt_path = hessian_root / "receipt.json"
    validated_receipt_path = receipt_path or canonical_receipt_path
    staged_receipt_name = None
    if validated_receipt_path != canonical_receipt_path:
        if (
            validated_receipt_path.parent != hessian_root
            or _HESSIAN_AGGREGATE_TEMPORARY_NAME.fullmatch(validated_receipt_path.name)
            is None
        ):
            raise ValueError("Hessian aggregate staged receipt path is invalid")
        staged_receipt_name = validated_receipt_path.name
    points_root = hessian_root / "points"
    _validate_hessian_tree_shape(
        hessian_root,
        points_root,
        len(path.coordinates_angstrom),
        staged_receipt_name=staged_receipt_name,
    )
    results: list[NativeHessianResult] = []
    child_hashes: list[str] = []
    for point_index in range(len(path.coordinates_angstrom)):
        result, child_hash = _validate_hessian_point(
            points_root / f"{point_index:06d}",
            campaign_identity=campaign_identity,
            route=route,
            atom_mapping_sha256=atom_mapping_sha256,
            qualified_transition_state=qualified_transition_state,
            path=path,
            point_index=point_index,
            settings_fingerprint=settings_fingerprint,
            backend_policy=backend_policy,
        )
        results.append(result)
        child_hashes.append(child_hash)
    receipt, receipt_bytes = _read_json_object(
        validated_receipt_path, label="Hessian aggregate receipt"
    )
    expected = _aggregate_hessian_receipt(
        campaign_identity=campaign_identity,
        route=route,
        atom_mapping_sha256=atom_mapping_sha256,
        path=path,
        settings_fingerprint=settings_fingerprint,
        backend_policy=backend_policy,
        child_hashes=child_hashes,
    )
    try:
        _strict_json_equal(receipt, expected, label="Hessian aggregate receipt")
    except ValueError as exc:
        try:
            _strict_json_equal(
                receipt.get("ordered_child_hash_chain"),
                expected["ordered_child_hash_chain"],
                label="Hessian aggregate child hash chain",
            )
        except ValueError:
            raise ValueError("Hessian aggregate child hash chain mismatch") from exc
        raise ValueError("Hessian aggregate receipt identity mismatch") from exc
    return PublishedPathHessians(
        receipt_path=validated_receipt_path,
        receipt_sha256=hashlib.sha256(receipt_bytes).hexdigest(),
        receipt=receipt,
        point_results=tuple(results),
    )


def publish_path_hessians(
    run_root: Path,
    *,
    campaign_identity: str,
    route: str,
    atom_mapping_sha256: str,
    qualified_transition_state: Cluster,
    path: PublishedTypedIrcPath,
    settings_fingerprint: str,
    backend_policy: dict[str, Any],
    evaluator: Callable[[Cluster], NativeHessianResult],
    _failure_injector: Callable[[str], None] | None = None,
) -> PublishedPathHessians:
    """Publish crash-resumable immutable native Hessians for every typed path point."""

    _require_sha(campaign_identity, length=64, label="campaign identity")
    _require_sha(atom_mapping_sha256, length=64, label="atom mapping SHA-256")
    if not isinstance(settings_fingerprint, str) or not settings_fingerprint:
        raise ValueError("native Hessian settings fingerprint must be non-empty")
    policy = _validated_backend_policy(backend_policy)
    _route_identity_cluster(route, qualified_transition_state, label="qualified TS")
    root = _safe_absolute_root(run_root)
    route_root = root / route
    with _exclusive_route_claim(route_root):
        current_path = _validate_published_path(
            route_root / "path",
            campaign_identity=campaign_identity,
            route=route,
            atom_mapping_sha256=atom_mapping_sha256,
            qualified_transition_state=qualified_transition_state,
        )
        if current_path.receipt_sha256 != path.receipt_sha256 or not np.array_equal(
            current_path.coordinates_angstrom, path.coordinates_angstrom
        ):
            raise ValueError(
                "supplied typed path does not match canonical path receipt"
            )
        path = current_path
        hessian_root = route_root / "hessians"
        if hessian_root.is_symlink():
            raise ValueError("Hessian stage root must not be a symlink")
        hessian_root_existed = hessian_root.exists()
        hessian_root.mkdir(mode=0o700, exist_ok=True)
        if not hessian_root_existed:
            _fsync_directory(route_root)
        points_root = hessian_root / "points"
        if points_root.is_symlink():
            raise ValueError("Hessian points root must not be a symlink")
        points_root_existed = points_root.exists()
        points_root.mkdir(mode=0o700, exist_ok=True)
        if not points_root_existed:
            _fsync_directory(hessian_root)
        _remove_owned_temporary_files(
            hessian_root, name_pattern=_HESSIAN_AGGREGATE_TEMPORARY_NAME
        )
        _remove_owned_temporary_directories(
            points_root,
            name_pattern=_HESSIAN_POINT_TEMPORARY_NAME,
            allowed_files={"gradient.f64", "hessian.f64", "receipt.json"},
            maximum_point_index=len(path.coordinates_angstrom),
        )
        _validate_hessian_tree_shape(
            hessian_root, points_root, len(path.coordinates_angstrom)
        )
        aggregate_path = hessian_root / "receipt.json"
        if aggregate_path.exists() or aggregate_path.is_symlink():
            return _validate_hessian_aggregate(
                hessian_root,
                campaign_identity=campaign_identity,
                route=route,
                atom_mapping_sha256=atom_mapping_sha256,
                qualified_transition_state=qualified_transition_state,
                path=path,
                settings_fingerprint=settings_fingerprint,
                backend_policy=policy,
            )

        results: list[NativeHessianResult] = []
        child_hashes: list[str] = []
        for point_index, coordinates in enumerate(path.coordinates_angstrom):
            point_dir = points_root / f"{point_index:06d}"
            if point_dir.exists() or point_dir.is_symlink():
                result, child_hash = _validate_hessian_point(
                    point_dir,
                    campaign_identity=campaign_identity,
                    route=route,
                    atom_mapping_sha256=atom_mapping_sha256,
                    qualified_transition_state=qualified_transition_state,
                    path=path,
                    point_index=point_index,
                    settings_fingerprint=settings_fingerprint,
                    backend_policy=policy,
                )
            else:
                cluster = _cluster_at_coordinates(
                    qualified_transition_state,
                    coordinates,
                    name=f"{route}-path-point-{point_index}",
                )
                result = _validate_evaluated_native_hessian(
                    evaluator(cluster),
                    cluster,
                    settings_fingerprint=settings_fingerprint,
                    backend_policy=policy,
                )
                receipt, gradient_bytes, hessian_bytes = _point_hessian_receipt(
                    campaign_identity=campaign_identity,
                    route=route,
                    atom_mapping_sha256=atom_mapping_sha256,
                    qualified_transition_state=qualified_transition_state,
                    path=path,
                    point_index=point_index,
                    cluster=cluster,
                    settings_fingerprint=settings_fingerprint,
                    backend_policy=policy,
                    result=result,
                )
                temporary = points_root / (
                    f".{point_index:06d}.{os.getpid()}.{time.time_ns()}.tmp"
                )
                temporary.mkdir(mode=0o700)
                try:
                    _write_fsync(temporary / "gradient.f64", gradient_bytes)
                    _write_fsync(temporary / "hessian.f64", hessian_bytes)
                    _write_fsync(temporary / "receipt.json", _json_bytes(receipt))
                    _fsync_directory(temporary)
                    _validate_hessian_point(
                        temporary,
                        campaign_identity=campaign_identity,
                        route=route,
                        atom_mapping_sha256=atom_mapping_sha256,
                        qualified_transition_state=qualified_transition_state,
                        path=path,
                        point_index=point_index,
                        settings_fingerprint=settings_fingerprint,
                        backend_policy=policy,
                    )
                    if _failure_injector is not None:
                        _failure_injector(f"before_hessian_point_commit:{point_index}")
                    if point_dir.exists() or point_dir.is_symlink():
                        raise FileExistsError(
                            f"Hessian point publication already exists: {point_dir}"
                        )
                    _publish_noreplace(temporary, point_dir, source_is_directory=True)
                    _fsync_directory(points_root)
                finally:
                    if temporary.exists():
                        shutil.rmtree(temporary)
                result, child_hash = _validate_hessian_point(
                    point_dir,
                    campaign_identity=campaign_identity,
                    route=route,
                    atom_mapping_sha256=atom_mapping_sha256,
                    qualified_transition_state=qualified_transition_state,
                    path=path,
                    point_index=point_index,
                    settings_fingerprint=settings_fingerprint,
                    backend_policy=policy,
                )
            results.append(result)
            child_hashes.append(child_hash)

        aggregate = _aggregate_hessian_receipt(
            campaign_identity=campaign_identity,
            route=route,
            atom_mapping_sha256=atom_mapping_sha256,
            path=path,
            settings_fingerprint=settings_fingerprint,
            backend_policy=policy,
            child_hashes=child_hashes,
        )
        if _failure_injector is not None:
            _failure_injector("before_hessian_aggregate_commit")

        def validate_staged_aggregate(receipt_path: Path) -> None:
            _validate_hessian_aggregate(
                hessian_root,
                campaign_identity=campaign_identity,
                route=route,
                atom_mapping_sha256=atom_mapping_sha256,
                qualified_transition_state=qualified_transition_state,
                path=path,
                settings_fingerprint=settings_fingerprint,
                backend_policy=policy,
                receipt_path=receipt_path,
            )

        _atomic_non_overwriting_file(
            aggregate_path,
            _json_bytes(aggregate),
            validate_before_commit=validate_staged_aggregate,
        )
        return _validate_hessian_aggregate(
            hessian_root,
            campaign_identity=campaign_identity,
            route=route,
            atom_mapping_sha256=atom_mapping_sha256,
            qualified_transition_state=qualified_transition_state,
            path=path,
            settings_fingerprint=settings_fingerprint,
            backend_policy=policy,
        )


def validate_transition_state_gate(
    cluster: Cluster,
    masses_amu: Any,
    native_hessian: NativeHessianResult,
    *,
    expected_settings_fingerprint: str,
    reaction_vector_mass_scaled: Any,
    reaction_vector_source: str,
    mapped_reactant_coordinates_angstrom: Any,
    mapped_product_coordinates_angstrom: Any,
) -> dict[str, Any]:
    """Apply the strict fresh first-order-saddle gate before any IRC call."""

    coordinates = np.asarray(cluster.coords, dtype=float)
    if native_hessian.geometry_fingerprint != frequency_geometry_fingerprint(cluster):
        raise ValueError(
            "native Hessian geometry fingerprint does not match TS geometry"
        )
    if (
        not expected_settings_fingerprint
        or native_hessian.settings_fingerprint != expected_settings_fingerprint
    ):
        raise ValueError("native Hessian settings fingerprint does not match campaign")
    masses = np.asarray(masses_amu, dtype=float)
    expected_masses = np.asarray(
        [ISOTOPIC_MASSES_AMU[symbol] for symbol in cluster.symbols], dtype=float
    )
    if masses.shape != expected_masses.shape or not np.array_equal(
        masses, expected_masses
    ):
        raise ValueError(
            "mass vector does not match the receipt-bound isotopic standard"
        )
    gradient_fmax = float(
        np.max(np.linalg.norm(native_hessian.gradient_hartree_per_bohr, axis=1))
    )
    gradient_fmax *= HARTREE_TO_EV / BOHR_TO_ANGSTROM
    if not math.isclose(
        gradient_fmax,
        native_hessian.physical_fmax_ev_per_angstrom,
        rel_tol=1.0e-12,
        abs_tol=1.0e-12,
    ):
        raise ValueError(
            "native Hessian physical fmax is stale relative to its gradient"
        )
    bounds = BOUNDS["transition_state_qualification"]
    fmax_limit = bounds["physical_fmax_ev_per_angstrom_exclusive_maximum"]
    if gradient_fmax >= fmax_limit:
        raise ValueError(
            "transition state is not stationary: physical fmax "
            f"{gradient_fmax:.12g} >= "
            f"{fmax_limit:.12g} eV/A"
        )
    modes = project_vibrational_hessian(
        coordinates,
        masses,
        native_hessian.cartesian_hessian_hartree_per_bohr2,
    )
    tolerance = bounds["negative_eigenvalue_tolerance_hartree_per_bohr2_amu"]
    negative = modes.eigenvalues < -tolerance
    near_zero = np.abs(modes.eigenvalues) <= tolerance
    if int(np.count_nonzero(negative)) != 1 or bool(np.any(near_zero)):
        raise ValueError(
            "transition state must have exactly one significant negative "
            "full-vibrational mode and no zero/noise-floor modes"
        )
    if bool(np.any(modes.eigenvalues[~negative] <= 0.0)):
        raise ValueError(
            "every non-reaction full-vibrational mode must be strictly positive"
        )
    wavenumbers = hessian_eigenvalues_to_wavenumbers_cm(modes.eigenvalues)
    imaginary = float(abs(wavenumbers[negative][0]))
    minimum_imaginary = bounds["minimum_reaction_imaginary_wavenumber_cm"]
    if imaginary < minimum_imaginary:
        raise ValueError(
            f"reaction imaginary mode {imaginary:.12g} cm^-1 is below "
            f"{minimum_imaginary:.12g} cm^-1"
        )
    reactant = np.asarray(mapped_reactant_coordinates_angstrom, dtype=float)
    product = np.asarray(mapped_product_coordinates_angstrom, dtype=float)
    if (
        reactant.shape != coordinates.shape
        or product.shape != coordinates.shape
        or not np.all(np.isfinite(reactant))
        or not np.all(np.isfinite(product))
    ):
        raise ValueError("mapped route basin geometries must be finite (N,3) arrays")
    mapped_path = build_mass_scaled_path(
        np.stack((reactant, coordinates, product)),
        masses,
        transition_state_index=1,
        reference_mass_amu=REFERENCE_MASS_AMU,
    )
    route_vector = (
        mapped_path.mass_scaled_coordinates[2] - mapped_path.mass_scaled_coordinates[0]
    )
    with np.errstate(over="ignore", invalid="ignore"):
        route_vector = modes.vibrational_basis @ (
            modes.vibrational_basis.T @ route_vector
        )
    if not np.all(np.isfinite(route_vector)):
        raise ValueError("mapped route vector projection must be finite")
    route_norm = float(np.linalg.norm(route_vector))
    if not math.isfinite(route_norm) or route_norm <= 1.0e-12:
        raise ValueError(
            "mapped route has no non-rigid reactant-to-product displacement"
        )
    route_vector /= route_norm

    reaction_vector = np.asarray(reaction_vector_mass_scaled, dtype=float)
    if reaction_vector.shape == coordinates.shape:
        reaction_vector = reaction_vector.reshape(-1)
    if reaction_vector.shape != (3 * len(cluster.symbols),) or not np.all(
        np.isfinite(reaction_vector)
    ):
        raise ValueError(
            "reaction vector must be a finite mass-scaled (N,3) or (3N,) vector"
        )
    with np.errstate(over="ignore", invalid="ignore"):
        reaction_vector = modes.vibrational_basis @ (
            modes.vibrational_basis.T @ reaction_vector
        )
    if not np.all(np.isfinite(reaction_vector)):
        raise ValueError("reaction vector projection must be finite")
    reaction_norm = float(np.linalg.norm(reaction_vector))
    if not math.isfinite(reaction_norm) or reaction_norm <= 1.0e-12:
        raise ValueError("reaction vector has no non-rigid vibrational component")
    reaction_vector /= reaction_norm
    if (
        not isinstance(reaction_vector_source, str)
        or not reaction_vector_source.strip()
    ):
        raise ValueError(
            "reaction vector source must identify its mapped route construction"
        )
    route_binding_overlap = float(abs(np.dot(reaction_vector, route_vector)))
    if (
        not math.isfinite(route_binding_overlap)
        or route_binding_overlap < 1.0 - 1.0e-10
    ):
        raise ValueError(
            "reaction vector does not match the receipt-bound mapped route geometries"
        )
    unstable_mode = modes.mass_weighted_eigenvectors[
        np.flatnonzero(negative)[0]
    ].reshape(-1)
    reaction_overlap = float(abs(np.dot(unstable_mode, reaction_vector)))
    minimum_overlap = bounds["minimum_mapped_reaction_vector_overlap"]
    if not math.isfinite(reaction_overlap) or reaction_overlap < minimum_overlap:
        raise ValueError(
            f"imaginary mode overlap {reaction_overlap:.12g} is below mapped-reaction "
            f"minimum {minimum_overlap:.12g}"
        )
    gradient_bytes = np.ascontiguousarray(
        native_hessian.gradient_hartree_per_bohr, dtype="<f8"
    ).tobytes()
    hessian_bytes = np.ascontiguousarray(
        native_hessian.cartesian_hessian_hartree_per_bohr2, dtype="<f8"
    ).tobytes()
    reaction_vector_bytes = np.ascontiguousarray(reaction_vector, dtype="<f8").tobytes()
    reactant_bytes = np.ascontiguousarray(reactant, dtype="<f8").tobytes()
    product_bytes = np.ascontiguousarray(product, dtype="<f8").tobytes()
    return {
        "accepted": True,
        "physical_fmax_ev_per_angstrom": gradient_fmax,
        "physical_fmax_exclusive_limit_ev_per_angstrom": fmax_limit,
        "vibrational_mode_count": int(modes.eigenvalues.size),
        "imaginary_mode_count": 1,
        "imaginary_wavenumber_cm": imaginary,
        "minimum_imaginary_wavenumber_cm": minimum_imaginary,
        "mapped_reaction_vector_overlap": reaction_overlap,
        "minimum_mapped_reaction_vector_overlap": minimum_overlap,
        "reaction_vector_source": reaction_vector_source,
        "reaction_vector_sha256": hashlib.sha256(reaction_vector_bytes).hexdigest(),
        "reaction_vector_route_binding_overlap": route_binding_overlap,
        "mapped_reactant_geometry_sha256": hashlib.sha256(reactant_bytes).hexdigest(),
        "mapped_product_geometry_sha256": hashlib.sha256(product_bytes).hexdigest(),
        "eigenvalues_hartree_per_bohr2_amu": modes.eigenvalues.tolist(),
        "masses_amu": masses.tolist(),
        "gradient_sha256": hashlib.sha256(gradient_bytes).hexdigest(),
        "canonical_hessian_sha256": hashlib.sha256(hessian_bytes).hexdigest(),
        "requested_backend": native_hessian.requested_backend,
        "actual_backend": native_hessian.actual_backend,
        "gpu_fallback_used": native_hessian.gpu_fallback_used,
        "geometry_fingerprint": native_hessian.geometry_fingerprint,
        "settings_fingerprint": native_hessian.settings_fingerprint,
    }


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


def _installed_distribution(label: str, candidates: tuple[str, ...]) -> tuple[str, str]:
    installed: list[tuple[str, str]] = []
    for distribution in candidates:
        try:
            installed.append((distribution, importlib.metadata.version(distribution)))
        except importlib.metadata.PackageNotFoundError:
            continue
    if not installed:
        raise RuntimeError(
            f"required D2c {label} distribution is not installed; expected one of "
            + ", ".join(candidates)
        )
    if len(installed) != 1:
        names = ", ".join(distribution for distribution, _ in installed)
        raise RuntimeError(
            f"ambiguous D2c {label} distributions are installed: {names}"
        )
    return installed[0]


def _dependency_versions() -> dict[str, str]:
    """Inventory the exact CPU packages, GPU backend, and live CUDA runtime."""

    versions: dict[str, str] = {}
    for distribution in DEPENDENCY_DISTRIBUTIONS:
        try:
            versions[distribution] = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError as exc:
            raise RuntimeError(
                f"required D2c dependency is not installed: {distribution}"
            ) from exc

    gpu_distribution, gpu_version = _installed_distribution(
        "GPU4PySCF", GPU4PYSCF_DISTRIBUTIONS
    )
    cupy_distribution, cupy_version = _installed_distribution(
        "CuPy", CUPY_DISTRIBUTIONS
    )
    try:
        __import__("gpu4pyscf")
        cupy = __import__("cupy")
        runtime = cupy.cuda.runtime
        runtime_version = runtime.runtimeGetVersion()
        driver_version = runtime.driverGetVersion()
        device_count = runtime.getDeviceCount()
    except Exception as exc:
        raise RuntimeError(
            "required D2c GPU backend or CUDA runtime is unavailable"
        ) from exc
    if not isinstance(device_count, int) or device_count < 1:
        raise RuntimeError("required D2c CUDA runtime exposes no GPU devices")

    versions.update(
        {
            "gpu_backend": "gpu4pyscf",
            "gpu4pyscf_distribution": gpu_distribution,
            "gpu4pyscf": gpu_version,
            "cupy_distribution": cupy_distribution,
            "cupy": cupy_version,
            "cuda_runtime": str(runtime_version),
            "cuda_driver": str(driver_version),
            "cuda_device_count": str(device_count),
        }
    )
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


def _semantic_atom_identity_labels(route: str, template: Any) -> tuple[str, ...]:
    """Bind the index semantics defined by the D2b reaction-template builders."""

    if template.n_water != 1:
        raise ValueError(f"D2c route must use exactly one template water: {route}")
    if template.family == "h-co":
        labels = (
            "carbon_monoxide_carbon",
            "carbon_monoxide_oxygen",
            "incoming_hydrogen",
            "water_1_oxygen",
            "water_1_hydrogen_a",
            "water_1_hydrogen_b",
        )
        expected_symbols = ("C", "O", "H", "O", "H", "H")
        expected_scan_pair = (0, 2)
    elif template.family == "h-h2co-ch3o":
        labels = (
            "formaldehyde_carbon",
            "formaldehyde_oxygen",
            "formaldehyde_hydrogen_a",
            "formaldehyde_hydrogen_b",
            "incoming_hydrogen",
            "water_1_oxygen",
            "water_1_hydrogen_a",
            "water_1_hydrogen_b",
        )
        expected_symbols = ("C", "O", "H", "H", "H", "O", "H", "H")
        expected_scan_pair = (0, 4)
    elif template.family == "h-h2co-h2-hco":
        labels = (
            "formaldehyde_carbon",
            "formaldehyde_oxygen",
            "abstracted_formaldehyde_hydrogen",
            "retained_formaldehyde_hydrogen",
            "incoming_hydrogen",
            "water_1_oxygen",
            "water_1_hydrogen_a",
            "water_1_hydrogen_b",
        )
        expected_symbols = ("C", "O", "H", "H", "H", "O", "H", "H")
        expected_scan_pair = (2, 4)
    else:
        raise ValueError(f"unsupported D2c reaction-template family: {route}")
    if (
        tuple(template.cluster.symbols) != expected_symbols
        or (template.scan_i, template.scan_j) != expected_scan_pair
    ):
        raise ValueError(f"semantic atom template order drifted: {route}")
    return labels


def _frozen_endpoint_evidence(
    bundle_root: Path,
    manifest_route: dict[str, Any],
    route: str,
    template: Any,
) -> dict[str, dict[str, Any]]:
    """Load and independently type both frozen endpoints for one route."""

    trusted = TRUSTED_FROZEN_ENDPOINT_EVIDENCE.get(route)
    if trusted is None or set(trusted) != {"irc_fwd.xyz", "irc_back.xyz"}:
        raise RuntimeError(f"trusted frozen endpoint inventory is incomplete: {route}")
    declared_geometry = manifest_route.get("canonical_checkpoint_geometry_sha256")
    declared_files = manifest_route.get("files")
    if not isinstance(declared_geometry, dict) or not isinstance(declared_files, dict):
        raise ValueError(f"frozen endpoint manifest evidence is missing: {route}")

    evidence: dict[str, dict[str, Any]] = {}
    for filename in ("irc_fwd.xyz", "irc_back.xyz"):
        path = bundle_root / route / filename
        expected = trusted[filename]
        declared_file = declared_files.get(filename)
        observed_geometry = d2c_input_bundle.geometry_hash_xyz(path)
        observed_file = d2c_input_bundle.sha256_path(path)
        if (
            declared_geometry.get(filename) != expected["geometry_sha256"]
            or not isinstance(declared_file, dict)
            or declared_file.get("sha256") != expected["file_sha256"]
            or observed_geometry != expected["geometry_sha256"]
            or observed_file != expected["file_sha256"]
        ):
            raise ValueError(
                f"trusted frozen endpoint evidence drifted: {route}/{filename}"
            )
        endpoint = load_xyz_like(path, template.cluster, name=f"{route}-{filename}")
        classification = classify_endpoint_basin(route, endpoint)
        if classification.basin != expected["basin"]:
            raise ValueError(
                f"trusted frozen endpoint chemistry label drifted: {route}/{filename}"
            )
        evidence[filename] = {
            "basin": classification.basin,
            "trusted_geometry_sha256": expected["geometry_sha256"],
            "trusted_file_sha256": expected["file_sha256"],
            "geometry_fingerprint": frequency_geometry_fingerprint(endpoint),
            "covalent_edges": [list(edge) for edge in classification.covalent_edges],
            "minimum_distance_angstrom": classification.minimum_distance_angstrom,
        }
    basins = [record["basin"] for record in evidence.values()]
    if sorted(basins) != ["product", "reactant"]:
        raise ValueError(
            f"frozen endpoints must prove exactly one reactant and one product: {route}"
        )
    return evidence


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
        expected_geometry_hash = TRUSTED_CANONICAL_TS_GEOMETRY_SHA256[route]
        declared_canonical_hashes = manifest_route.get(
            "canonical_checkpoint_geometry_sha256"
        )
        observed_geometry_hash = d2c_input_bundle.geometry_hash_xyz(
            bundle_root / route / "ts.xyz"
        )
        if (
            not isinstance(declared_canonical_hashes, dict)
            or declared_canonical_hashes.get("ts.xyz") != expected_geometry_hash
            or observed_geometry_hash != expected_geometry_hash
        ):
            raise ValueError(
                f"canonical transition-state geometry/mapping identity drifted: {route}"
            )
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
        atom_identity_labels = _semantic_atom_identity_labels(route, template)
        frozen_endpoint_evidence = _frozen_endpoint_evidence(
            bundle_root, manifest_route, route, template
        )
        if tuple(transition_state.symbols) != tuple(template.cluster.symbols):
            raise ValueError(f"transition-state atom symbols drifted: {route}")
        atom_mapping_sha256 = _canonical_hash(
            {
                "route": route,
                "canonical_transition_state_geometry_sha256": expected_geometry_hash,
                "atoms": [
                    {"index": index, "symbol": symbol, "semantic_identity": label}
                    for index, (symbol, label) in enumerate(
                        zip(
                            transition_state.symbols,
                            atom_identity_labels,
                            strict=True,
                        )
                    )
                ],
            }
        )
        inventory[route] = {
            "symbols": list(transition_state.symbols),
            "atom_identity_labels": list(atom_identity_labels),
            "atom_mapping_sha256": atom_mapping_sha256,
            "canonical_transition_state_geometry_sha256": expected_geometry_hash,
            "masses_amu": masses,
            "transition_state_sha256": d2c_input_bundle.sha256_path(
                bundle_root / route / "ts.xyz"
            ),
            "frozen_endpoint_evidence": frozen_endpoint_evidence,
            "stages": _pending_stages(ROUTE_STAGE_CONTRACT),
        }
    if tuple(inventory) != d2c_input_bundle.ROUTES or len(inventory) != 4:
        raise RuntimeError(
            "D2c preflight must enumerate exactly the four frozen routes"
        )
    return inventory


@contextmanager
def _exclusive_run_claim(run_root: Path) -> Iterator[None]:
    """Hold a crash-recoverable exclusive claim on one campaign run root."""

    run_root.mkdir(parents=True, exist_ok=True)
    claim_path = run_root / ".preflight.lock"
    flags = os.O_RDWR | os.O_CREAT | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(claim_path, flags, 0o600)
    locked = False
    try:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(
                f"D2c run root is already claimed by another preflight: {run_root}"
            ) from exc
        locked = True
        claim = json.dumps(
            {"pid": os.getpid(), "claimed_unix_ns": time.time_ns()},
            sort_keys=True,
        ).encode()
        os.ftruncate(descriptor, 0)
        os.write(descriptor, claim + b"\n")
        os.fsync(descriptor)
        yield
    finally:
        if locked:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        os.close(descriptor)


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
    if set(dependencies) != DEPENDENCY_VERSION_KEYS:
        raise ValueError("dependency version inventory must be exact and complete")
    if any(not isinstance(value, str) or not value for value in dependencies.values()):
        raise ValueError("dependency versions must be non-empty strings")
    if (
        dependencies["gpu_backend"] != "gpu4pyscf"
        or dependencies["gpu4pyscf_distribution"] not in GPU4PYSCF_DISTRIBUTIONS
        or dependencies["cupy_distribution"] not in CUPY_DISTRIBUTIONS
        or not dependencies["cuda_runtime"].isdigit()
        or not dependencies["cuda_driver"].isdigit()
    ):
        raise ValueError("dependency inventory has an invalid GPU backend identity")
    try:
        cuda_device_count = int(dependencies["cuda_device_count"])
    except ValueError as exc:
        raise ValueError("CUDA device count must be a positive integer") from exc
    if cuda_device_count < 1:
        raise ValueError("CUDA device count must be a positive integer")
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
                "atom_identity_labels": record["atom_identity_labels"],
                "atom_mapping_sha256": record["atom_mapping_sha256"],
                "canonical_transition_state_geometry_sha256": record[
                    "canonical_transition_state_geometry_sha256"
                ],
                "masses_amu": record["masses_amu"],
                "transition_state_sha256": record["transition_state_sha256"],
                "frozen_endpoint_evidence": record["frozen_endpoint_evidence"],
            }
            for route, record in routes.items()
        },
        "endpoint_classification_policy": endpoint_classification_policy_payload(),
        "trusted_frozen_endpoint_evidence": TRUSTED_FROZEN_ENDPOINT_EVIDENCE,
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

    with _exclusive_run_claim(run_root):
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
