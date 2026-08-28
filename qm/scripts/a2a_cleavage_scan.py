#!/usr/bin/env python3
"""Resume-safe 2-D r2SCAN-3c scan for proton-coupled A2a cleavage.

The scan binds the actual attacking-water proton H16, not the already
bridge-owned H17, and relaxes every spectator degree of freedom while holding
``d(Si-Obr)`` and ``d(H16-Obr)``.  A deterministic serpentine traversal makes
all seeds adjacent.  Every cell is independently checkpointed and fingerprinted
before a complete-grid minimax classifier is allowed to emit a diagnostic.
A coarse scan crest is only a localization seed; it is never a verified saddle.
"""

from __future__ import annotations

import argparse
import hashlib
import heapq
import json
import math
import time
from collections.abc import Callable
from dataclasses import asdict, dataclass
from itertools import pairwise
from pathlib import Path
from typing import Any

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    _ETIQUETTE = bootstrap_cli(
        "a2a-cleavage-scan",
        default_run_root=(
            "/mnt/data/vsletten/dissertation-data/task208-a2a-path-rebuild-20260825"
        ),
        gpu_owner="a2a-cleavage-scan",
    )

import numpy as np

from quarry.clusters import Cluster
from quarry.pipeline import (
    HARTREE_TO_KJ,
    DftSettings,
    frequency_geometry_fingerprint,
    frequency_settings_fingerprint,
)
from quarry.ts import relax_at_fixed_distances
from scripts import a2a_path_rebuild as a2a
from scripts import production_energetics as a2

SCAN_VERSION = "a2a-coupled-cleavage-scan-v1"
SCAN_RELATIVE_PATH = Path("cleavage") / "coupled-scan-v1"
DOWNHILL_RELEASE_VERSION = "a2a-barrierless-downhill-release-v2"
DOWNHILL_RELEASE_RELATIVE_PATH = Path("cleavage") / "barrierless-downhill-release-v2"
DEFAULT_AXIS_POINTS = 9
DEFAULT_BARRIER_THRESHOLD_KJ_MOL = 2.0
DEFAULT_FMAX_EV_A = 0.02
DEFAULT_MAX_STEPS = 120
DEFAULT_DOWNHILL_MAX_STEPS = 200
DEFAULT_OPTIMIZER_MAXSTEP_A = 0.03
DEFAULT_DISTANCE_TOLERANCE_A = 1.0e-4
DOWNHILL_MINIMUM_DISPLACEMENT_A = 1.0e-3
DOWNHILL_ENERGY_TOLERANCE_HARTREE = 1.0e-8
PROTON_SELECTION_MIN_CONTRACTION_A = 0.05
PROTON_SELECTION_MARGIN_A = 0.05

SI_INDEX = 1
BRIDGE_INDEX = 0
OW_INDEX = 15
HW_INDEX = 16
DEGENERATE_HW_INDEX = 17

PROTON_FIRST_MINIMUM = "proton-first-minimum"
INTERIOR_CREST = "interior-crest-above-threshold"
BARRIERLESS_SHELF = "barrierless-shelf"

RelaxFunction = Callable[..., Cluster]
OptimizeFunction = Callable[..., Cluster]
EnergyCheckpointFunction = Callable[[Path, Cluster, DftSettings, str], float]
FrequencyCheckpointFunction = Callable[..., Any]


@dataclass(frozen=True)
class CoordinateMapping:
    si_index: int
    bridge_index: int
    ow_index: int
    hw_index: int
    rejected_degenerate_hw_index: int
    si_obr_pair: tuple[int, int]
    hw_obr_pair: tuple[int, int]
    intermediate_si_obr_a: float
    product_si_obr_a: float
    intermediate_hw_obr_a: float
    product_hw_obr_a: float
    intermediate_h17_obr_a: float
    product_h17_obr_a: float
    hw_contraction_a: float
    h17_change_a: float


@dataclass(frozen=True)
class GridSpec:
    si_obr_targets_a: tuple[float, ...]
    hw_obr_targets_a: tuple[float, ...]
    product_index: int
    traversal: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class CellRecord:
    row: int
    column: int
    si_obr_target_a: float
    hw_obr_target_a: float
    electronic_hartree: float
    seed_geometry_fingerprint: str
    output_geometry_fingerprint: str
    previous_cell: tuple[int, int] | None
    geometry_path: str
    receipt_path: str
    topology: dict[str, Any]
    cluster: Cluster


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    """Reuse the campaign's atomic JSON primitive."""
    a2.atomic_json(path, payload)


def _json_stable(payload: Any) -> Any:
    """Return the exact JSON value used for on-disk equality checks."""
    return json.loads(json.dumps(payload, sort_keys=True, allow_nan=False))


def _payload_fingerprint(payload: Any) -> str:
    serialized = json.dumps(
        _json_stable(payload),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    return hashlib.sha256(serialized.encode()).hexdigest()


def _distance(cluster: Cluster, left: int, right: int) -> float:
    return float(np.linalg.norm(cluster.coords[left] - cluster.coords[right]))


def _owners(cluster: Cluster) -> dict[int, int]:
    """Obtain fail-closed physical-H ownership from the typed endpoint helper."""
    _basin, owners, _heavy_bonds = a2.endpoint_identity(cluster, OW_INDEX)
    result = dict(owners)
    if len(result) != len(owners):
        raise ValueError(f"{cluster.name}: duplicate physical-H ownership records")
    return result


def validate_coordinate_mapping(
    intermediate: Cluster,
    product: Cluster,
    *,
    hw_index: int = HW_INDEX,
) -> CoordinateMapping:
    """Validate and return the explicit H16 coupled-coordinate mapping.

    H17 is already Obr-owned at both endpoints and barely changes distance.  It
    is therefore a degenerate coordinate for the commissioned scan.  No nearest-
    proton heuristic or fallback is allowed: callers must use physical H16.
    """
    if hw_index == DEGENERATE_HW_INDEX:
        raise ValueError("H17 is the degenerate already-Obr-owned proton; bind H16")
    if hw_index != HW_INDEX:
        raise ValueError(
            f"coupled cleavage scan must explicitly bind H16, not H{hw_index}"
        )
    if intermediate.symbols != product.symbols:
        raise ValueError("I/P atom order differs")
    if intermediate.charge != product.charge or intermediate.spin != product.spin:
        raise ValueError("I/P electronic states differ")
    required = {
        BRIDGE_INDEX: "O",
        SI_INDEX: "Si",
        OW_INDEX: "O",
        HW_INDEX: "H",
        DEGENERATE_HW_INDEX: "H",
    }
    if max(required) >= len(intermediate.symbols):
        raise ValueError("I/P atom list is too short for the fixed cleavage mapping")
    for index, symbol in required.items():
        if intermediate.symbols[index] != symbol:
            raise ValueError(
                f"fixed cleavage mapping expected {symbol} at index {index}, "
                f"found {intermediate.symbols[index]}"
            )

    intermediate_owners = _owners(intermediate)
    product_owners = _owners(product)
    for role, owners in (("I", intermediate_owners), ("P", product_owners)):
        if owners.get(HW_INDEX) != OW_INDEX:
            raise ValueError(f"{role}: H16 is not uniquely water-owned by O15")
        if owners.get(DEGENERATE_HW_INDEX) != BRIDGE_INDEX:
            raise ValueError(f"{role}: H17 is not uniquely bridge-owned by O0")

    i_si = _distance(intermediate, SI_INDEX, BRIDGE_INDEX)
    p_si = _distance(product, SI_INDEX, BRIDGE_INDEX)
    i_hw = _distance(intermediate, HW_INDEX, BRIDGE_INDEX)
    p_hw = _distance(product, HW_INDEX, BRIDGE_INDEX)
    i_h17 = _distance(intermediate, DEGENERATE_HW_INDEX, BRIDGE_INDEX)
    p_h17 = _distance(product, DEGENERATE_HW_INDEX, BRIDGE_INDEX)
    values = (i_si, p_si, i_hw, p_hw, i_h17, p_h17)
    if any(not math.isfinite(value) or value <= 0.0 for value in values):
        raise ValueError(
            "I/P cleavage-coordinate distances must be finite and positive"
        )
    if p_si - i_si <= PROTON_SELECTION_MIN_CONTRACTION_A:
        raise ValueError("I/P does not span Si-Obr cleavage")
    hw_contraction = i_hw - p_hw
    h17_change = abs(i_h17 - p_h17)
    if hw_contraction <= max(
        PROTON_SELECTION_MIN_CONTRACTION_A,
        h17_change + PROTON_SELECTION_MARGIN_A,
    ):
        raise ValueError(
            "degenerate H16/H17 endpoint motion: H16 is not the unique "
            "contracting water-owned proton"
        )
    return CoordinateMapping(
        si_index=SI_INDEX,
        bridge_index=BRIDGE_INDEX,
        ow_index=OW_INDEX,
        hw_index=HW_INDEX,
        rejected_degenerate_hw_index=DEGENERATE_HW_INDEX,
        si_obr_pair=(SI_INDEX, BRIDGE_INDEX),
        hw_obr_pair=(HW_INDEX, BRIDGE_INDEX),
        intermediate_si_obr_a=i_si,
        product_si_obr_a=p_si,
        intermediate_hw_obr_a=i_hw,
        product_hw_obr_a=p_hw,
        intermediate_h17_obr_a=i_h17,
        product_h17_obr_a=p_h17,
        hw_contraction_a=hw_contraction,
        h17_change_a=h17_change,
    )


def _axis_targets(initial: float, product: float, points: int) -> tuple[float, ...]:
    if points < 3:
        raise ValueError("axis-points must be at least 3")
    if not math.isfinite(initial) or not math.isfinite(product):
        raise ValueError("grid endpoints must be finite")
    if initial <= 0.0 or product <= 0.0:
        raise ValueError("grid distances must be positive")
    span = product - initial
    if abs(span) < 1.0e-8:
        raise ValueError("grid endpoint targets are degenerate")
    product_index = points - 2
    step = span / product_index
    values = [initial + step * index for index in range(points)]
    values[0] = initial
    values[product_index] = product
    values[-1] = product + step
    if any(not math.isfinite(value) or value <= 0.0 for value in values):
        raise ValueError("grid extension produced a non-positive distance")
    return tuple(values)


def build_grid(
    intermediate_si_obr_a: float,
    product_si_obr_a: float,
    intermediate_hw_obr_a: float,
    product_hw_obr_a: float,
    *,
    axis_points: int = DEFAULT_AXIS_POINTS,
) -> GridSpec:
    """Build a square I-to-P grid with one exact step of extension past P."""
    si_targets = _axis_targets(intermediate_si_obr_a, product_si_obr_a, axis_points)
    hw_targets = _axis_targets(intermediate_hw_obr_a, product_hw_obr_a, axis_points)
    traversal: list[tuple[int, int]] = []
    for row in range(axis_points):
        columns = range(axis_points) if row % 2 == 0 else range(axis_points - 1, -1, -1)
        traversal.extend((row, column) for column in columns)
    if any(
        abs(left[0] - right[0]) + abs(left[1] - right[1]) != 1
        for left, right in pairwise(traversal)
    ):
        raise AssertionError("internal error: scan traversal is not adjacent")
    return GridSpec(
        si_obr_targets_a=si_targets,
        hw_obr_targets_a=hw_targets,
        product_index=axis_points - 2,
        traversal=tuple(traversal),
    )


def scan_root(run_dir: Path) -> Path:
    return run_dir / SCAN_RELATIVE_PATH


def _endpoint_payload(cluster: Cluster) -> dict[str, Any]:
    basin, owners, heavy_bonds = a2.endpoint_identity(cluster, OW_INDEX)
    return {
        "geometry_fingerprint": frequency_geometry_fingerprint(cluster),
        "basin": list(basin),
        "hydrogen_owners": [list(pair) for pair in owners],
        "heavy_bonds": [list(pair) for pair in heavy_bonds],
    }


def _cell_topology(cluster: Cluster) -> dict[str, Any]:
    """Record typed cell topology without promoting it to an endpoint gate."""
    try:
        basin, owners, heavy_bonds = a2.endpoint_identity(cluster, OW_INDEX)
    except ValueError as exc:
        return {
            "valid_typed_identity": False,
            "identity_error": str(exc),
            "h16_owner": None,
            "si_obr_bonded": None,
        }
    owner_map = dict(owners)
    normalized_heavy_bonds = {tuple(sorted(pair)) for pair in heavy_bonds}
    return {
        "valid_typed_identity": True,
        "basin": list(basin),
        "hydrogen_owners": [list(pair) for pair in owners],
        "heavy_bonds": [list(pair) for pair in heavy_bonds],
        "h16_owner": owner_map.get(HW_INDEX),
        "si_obr_bonded": tuple(sorted((SI_INDEX, BRIDGE_INDEX)))
        in normalized_heavy_bonds,
    }


def _require_exact_json(path: Path, expected: dict[str, Any], label: str) -> None:
    try:
        actual = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"{label} is unreadable: {exc}") from exc
    if actual != expected:
        raise ValueError(f"{label} drift")


def _manifest_payload(
    intermediate: Cluster,
    product: Cluster,
    settings: DftSettings,
    mapping: CoordinateMapping,
    grid: GridSpec,
    *,
    fmax_ev_a: float,
    max_steps: int,
    optimizer_maxstep_a: float,
    distance_tolerance_a: float,
) -> dict[str, Any]:
    return _json_stable(
        {
            "scan_version": SCAN_VERSION,
            "status": "configured",
            "coordinates": asdict(mapping),
            "endpoints": {
                "intermediate": _endpoint_payload(intermediate),
                "product": _endpoint_payload(product),
            },
            "settings": asdict(settings),
            "settings_fingerprint": frequency_settings_fingerprint(settings),
            "grid": {
                "axis_points": len(grid.si_obr_targets_a),
                "product_index": grid.product_index,
                "si_obr_targets_a": list(grid.si_obr_targets_a),
                "hw_obr_targets_a": list(grid.hw_obr_targets_a),
                "traversal": [list(cell) for cell in grid.traversal],
                "extension": "one-grid-step-beyond-exact-product-on-both-axes",
            },
            "relaxation": {
                "algorithm": "relax_at_fixed_distances",
                "constraint_method": "sella-internal",
                "fmax_ev_a": fmax_ev_a,
                "max_steps": max_steps,
                "optimizer_maxstep_a": optimizer_maxstep_a,
                "distance_tolerance_a": distance_tolerance_a,
            },
        }
    )


def run_scan_grid(
    scan_dir: Path,
    intermediate: Cluster,
    product: Cluster,
    settings: DftSettings,
    *,
    axis_points: int = DEFAULT_AXIS_POINTS,
    fmax_ev_a: float = DEFAULT_FMAX_EV_A,
    max_steps: int = DEFAULT_MAX_STEPS,
    optimizer_maxstep_a: float = DEFAULT_OPTIMIZER_MAXSTEP_A,
    distance_tolerance_a: float = DEFAULT_DISTANCE_TOLERANCE_A,
    relax_fn: RelaxFunction = relax_at_fixed_distances,
    checkpoint_energy_fn: EnergyCheckpointFunction = a2.checkpoint_energy,
) -> list[CellRecord]:
    """Run or resume every cell, failing closed on any identity drift."""
    if not math.isfinite(fmax_ev_a) or fmax_ev_a <= 0.0:
        raise ValueError("fmax must be finite and positive")
    if max_steps <= 0:
        raise ValueError("max-steps must be positive")
    if not math.isfinite(optimizer_maxstep_a) or optimizer_maxstep_a <= 0.0:
        raise ValueError("optimizer-maxstep must be finite and positive")
    if not math.isfinite(distance_tolerance_a) or distance_tolerance_a <= 0.0:
        raise ValueError("distance-tolerance must be finite and positive")

    mapping = validate_coordinate_mapping(intermediate, product)
    grid = build_grid(
        mapping.intermediate_si_obr_a,
        mapping.product_si_obr_a,
        mapping.intermediate_hw_obr_a,
        mapping.product_hw_obr_a,
        axis_points=axis_points,
    )
    scan_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = scan_dir / "scan-manifest.json"
    manifest = _manifest_payload(
        intermediate,
        product,
        settings,
        mapping,
        grid,
        fmax_ev_a=fmax_ev_a,
        max_steps=max_steps,
        optimizer_maxstep_a=optimizer_maxstep_a,
        distance_tolerance_a=distance_tolerance_a,
    )
    if manifest_path.exists():
        _require_exact_json(manifest_path, manifest, "scan manifest")
    else:
        atomic_json(manifest_path, manifest)

    settings_fp = frequency_settings_fingerprint(settings)
    endpoint_fps = {
        "intermediate": frequency_geometry_fingerprint(intermediate),
        "product": frequency_geometry_fingerprint(product),
    }
    records: list[CellRecord] = []
    seed = intermediate
    previous_cell: tuple[int, int] | None = None
    for row, column in grid.traversal:
        cell_dir = scan_dir / f"r{row:02d}-c{column:02d}"
        cell_dir.mkdir(parents=True, exist_ok=True)
        geometry_path = cell_dir / "optimized.xyz"
        energy_path = cell_dir / "electronic-energy.json"
        receipt_path = cell_dir / "cell-receipt.json"
        si_target = grid.si_obr_targets_a[row]
        hw_target = grid.hw_obr_targets_a[column]
        target_payload = {
            "si_obr_a": si_target,
            "hw_obr_a": hw_target,
            "pairs": [[SI_INDEX, BRIDGE_INDEX], [HW_INDEX, BRIDGE_INDEX]],
        }
        seed_fp = frequency_geometry_fingerprint(seed)
        identity = _json_stable(
            {
                "scan_version": SCAN_VERSION,
                "stage": "coupled-coordinate-relaxed-cell",
                "cell": [row, column],
                "previous_cell": (
                    list(previous_cell) if previous_cell is not None else None
                ),
                "target": target_payload,
                "target_fingerprint": _payload_fingerprint(target_payload),
                "settings_fingerprint": settings_fp,
                "seed_geometry_fingerprint": seed_fp,
                "endpoint_geometry_fingerprints": endpoint_fps,
                "relaxation": manifest["relaxation"],
            }
        )

        def compute(
            current_seed: Cluster = seed,
            current_si_target: float = si_target,
            current_hw_target: float = hw_target,
            current_cell_dir: Path = cell_dir,
        ) -> Cluster:
            return relax_fn(
                current_seed,
                settings,
                fixed_distances=[
                    (SI_INDEX, BRIDGE_INDEX, current_si_target),
                    (HW_INDEX, BRIDGE_INDEX, current_hw_target),
                ],
                fmax_ev_a=fmax_ev_a,
                max_steps=max_steps,
                optimizer_maxstep=optimizer_maxstep_a,
                distance_tolerance_a=distance_tolerance_a,
                constraint_method="sella-internal",
                trajectory=str(current_cell_dir / "relax.traj"),
                logfile=str(current_cell_dir / "relax.log"),
            )

        optimized = a2.checkpoint_cluster(
            geometry_path,
            seed,
            compute,
            identity=identity,
        )
        residuals = {
            "si_obr_a": abs(_distance(optimized, SI_INDEX, BRIDGE_INDEX) - si_target),
            "hw_obr_a": abs(_distance(optimized, HW_INDEX, BRIDGE_INDEX) - hw_target),
        }
        if max(residuals.values()) > distance_tolerance_a:
            raise RuntimeError(
                f"cell r{row:02d}-c{column:02d} target residual exceeds "
                f"{distance_tolerance_a:.3e} A: {residuals}"
            )
        electronic_hartree = checkpoint_energy_fn(
            energy_path,
            optimized,
            settings,
            a2.R2SCAN3C_METHOD,
        )
        if not math.isfinite(electronic_hartree):
            raise RuntimeError(f"cell r{row:02d}-c{column:02d} has non-finite energy")
        output_fp = frequency_geometry_fingerprint(optimized)
        topology = _cell_topology(optimized)
        receipt = _json_stable(
            {
                "status": "completed",
                "scan_version": SCAN_VERSION,
                "cell": [row, column],
                "previous_cell": (
                    list(previous_cell) if previous_cell is not None else None
                ),
                "target": target_payload,
                "target_fingerprint": identity["target_fingerprint"],
                "settings_fingerprint": settings_fp,
                "seed_geometry_fingerprint": seed_fp,
                "endpoint_geometry_fingerprints": endpoint_fps,
                "output_geometry_fingerprint": output_fp,
                "output_geometry_sha256": a2.sha256_path(geometry_path),
                "energy_receipt_sha256": a2.sha256_path(energy_path),
                "electronic_hartree": electronic_hartree,
                "target_residuals_a": residuals,
                "topology": topology,
                "coarse_grid_point_is_verified_saddle": False,
            }
        )
        if receipt_path.exists():
            _require_exact_json(
                receipt_path,
                receipt,
                f"cell r{row:02d}-c{column:02d} receipt",
            )
        else:
            atomic_json(receipt_path, receipt)
        records.append(
            CellRecord(
                row=row,
                column=column,
                si_obr_target_a=si_target,
                hw_obr_target_a=hw_target,
                electronic_hartree=electronic_hartree,
                seed_geometry_fingerprint=seed_fp,
                output_geometry_fingerprint=output_fp,
                previous_cell=previous_cell,
                geometry_path=str(geometry_path),
                receipt_path=str(receipt_path),
                topology=topology,
                cluster=optimized,
            )
        )
        seed = optimized
        previous_cell = (row, column)
    return records


def _neighbors(cell: tuple[int, int], size: int) -> list[tuple[int, int]]:
    row, column = cell
    candidates = (
        (row - 1, column),
        (row, column - 1),
        (row, column + 1),
        (row + 1, column),
    )
    return sorted(
        (left, right)
        for left, right in candidates
        if 0 <= left < size and 0 <= right < size
    )


def _minimax_path(
    energies: np.ndarray,
    start: tuple[int, int],
    goal: tuple[int, int],
) -> list[tuple[int, int]]:
    """Return a deterministic shortest lexicographic minimum-bottleneck path."""
    initial_path = (start,)
    initial_key = (float(energies[start]), 0, initial_path)
    best: dict[tuple[int, int], tuple[float, int, tuple[tuple[int, int], ...]]] = {
        start: initial_key
    }
    queue: list[tuple[float, int, tuple[tuple[int, int], ...], tuple[int, int]]] = [
        (initial_key[0], initial_key[1], initial_key[2], start)
    ]
    while queue:
        cost, steps, path, cell = heapq.heappop(queue)
        if best.get(cell) != (cost, steps, path):
            continue
        if cell == goal:
            return list(path)
        for neighbor in _neighbors(cell, len(energies)):
            candidate_path = (*path, neighbor)
            candidate = (
                max(cost, float(energies[neighbor])),
                steps + 1,
                candidate_path,
            )
            if neighbor not in best or candidate < best[neighbor]:
                best[neighbor] = candidate
                heapq.heappush(
                    queue,
                    (candidate[0], candidate[1], candidate[2], neighbor),
                )
    raise ValueError("complete grid does not contain an I-to-P path")


def _proton_first_minima(
    relative_energies: np.ndarray,
    product_index: int,
    cell_topologies: dict[tuple[int, int], dict[str, Any]],
) -> list[tuple[float, tuple[int, int]]]:
    candidates: list[tuple[float, tuple[int, int]]] = []
    size = len(relative_energies)
    for row in range(1, product_index):
        for column in range(1, product_index):
            cell = (row, column)
            topology = cell_topologies.get(cell, {})
            if not (
                topology.get("valid_typed_identity") is True
                and topology.get("h16_owner") == BRIDGE_INDEX
                and topology.get("si_obr_bonded") is True
            ):
                continue
            energy = float(relative_energies[cell])
            neighbor_energies = [
                float(relative_energies[neighbor])
                for neighbor in _neighbors(cell, size)
            ]
            if all(energy < neighbor for neighbor in neighbor_energies):
                candidates.append((energy, cell))
    return sorted(candidates, key=lambda item: (item[0], item[1]))


def classify_complete_grid(
    energies_kj_mol: np.ndarray,
    *,
    product_index: int,
    barrier_threshold_kj_mol: float = DEFAULT_BARRIER_THRESHOLD_KJ_MOL,
    cell_topologies: dict[tuple[int, int], dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Classify the complete I-to-P surface in the commissioned outcome order."""
    energies = np.asarray(energies_kj_mol, dtype=float)
    if energies.ndim != 2 or energies.shape[0] != energies.shape[1]:
        raise ValueError("complete-grid classifier requires a square energy matrix")
    size = energies.shape[0]
    if size < 3 or product_index <= 0 or product_index >= size - 1:
        raise ValueError("product index must be interior to an extended grid")
    if not np.all(np.isfinite(energies)):
        raise ValueError("complete-grid classifier rejects missing/non-finite cells")
    if not math.isfinite(barrier_threshold_kj_mol) or barrier_threshold_kj_mol <= 0.0:
        raise ValueError("barrier threshold must be finite and positive")

    start = (0, 0)
    goal = (product_index, product_index)
    relative = energies - float(energies[start])
    # Extension points diagnose behavior beyond P; they must not provide an
    # artificial low-energy detour around the physical I-to-P rectangle.
    physical_surface = relative[: product_index + 1, : product_index + 1]
    path = _minimax_path(physical_surface, start, goal)
    path_interior = path[1:-1]
    if path_interior:
        bottleneck_height = max(float(relative[cell]) for cell in path_interior)
        bottleneck = next(
            cell for cell in path_interior if float(relative[cell]) == bottleneck_height
        )
    else:
        bottleneck = start
        bottleneck_height = 0.0
    result: dict[str, Any] = {
        "classifier": "complete-grid-deterministic-minimax-v1",
        "barrier_threshold_kj_mol": barrier_threshold_kj_mol,
        "start_cell": list(start),
        "product_cell": list(goal),
        "minimax_path": [list(cell) for cell in path],
        "minimax_bottleneck": {
            "cell": list(bottleneck),
            "height_kj_mol": bottleneck_height,
            "is_grid_boundary": (
                bottleneck[0] in (0, product_index)
                or bottleneck[1] in (0, product_index)
            ),
        },
        "verified_saddle": False,
    }

    proton_minima = _proton_first_minima(
        physical_surface,
        product_index,
        cell_topologies or {},
    )
    if proton_minima:
        energy, cell = proton_minima[0]
        result.update(
            {
                "outcome": PROTON_FIRST_MINIMUM,
                "candidate_kind": "coarse-grid-proton-first-minimum",
                "proton_first_minimum": {
                    "cell": list(cell),
                    "relative_energy_kj_mol": energy,
                    "criterion": "H16-is-Obr-owned-while-Si-Obr-remains-bonded",
                },
            }
        )
        return result

    bottleneck_is_grid_interior = (
        0 < bottleneck[0] < product_index and 0 < bottleneck[1] < product_index
    )
    if (
        bottleneck not in (start, goal)
        and bottleneck_is_grid_interior
        and bottleneck_height > barrier_threshold_kj_mol
    ):
        result.update(
            {
                "outcome": INTERIOR_CREST,
                "candidate_kind": "coarse-grid-crest-seed-not-verified-saddle",
                "crest": {
                    "cell": list(bottleneck),
                    "height_kj_mol": bottleneck_height,
                    "follow_through": (
                        "seed one full-dimensional Sella localization; retain all "
                        "Hessian/IRC/basin gates"
                    ),
                },
            }
        )
        return result

    result.update(
        {
            "outcome": BARRIERLESS_SHELF,
            "candidate_kind": "complete-grid-valley-no-interior-crest-above-threshold",
            "shelf": {
                "maximum_interior_valley_height_kj_mol": bottleneck_height,
                "boundary_bottleneck_not_promoted": not bottleneck_is_grid_interior,
            },
        }
    )
    return result


def select_barrierless_release_seed(
    records: list[CellRecord],
    classification: dict[str, Any],
    product: Cluster,
    *,
    attacker_index: int = OW_INDEX,
) -> CellRecord:
    """Select the exact product cell at the end of the classified I-to-P valley."""
    if classification.get("outcome") != BARRIERLESS_SHELF:
        raise ValueError("downhill release requires a barrierless-shelf classification")
    if classification.get("verified_saddle") is not False:
        raise ValueError("barrierless classification must explicitly reject a saddle")

    def parse_cell(value: Any) -> tuple[int, int]:
        if not isinstance(value, (list, tuple)) or len(value) != 2:
            raise ValueError("cell must contain exactly two indices")
        return int(value[0]), int(value[1])

    try:
        path = [parse_cell(cell) for cell in classification["minimax_path"]]
        start = parse_cell(classification["start_cell"])
        goal = parse_cell(classification["product_cell"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(
            "barrierless classification has a malformed minimax path"
        ) from exc
    if not path or path[0] != start or path[-1] != goal:
        raise ValueError(
            "barrierless classification path does not span exact I-to-P cells"
        )

    by_cell: dict[tuple[int, int], CellRecord] = {}
    for record in records:
        cell = (record.row, record.column)
        if cell in by_cell:
            raise ValueError(f"duplicate scan record for cell {cell}")
        by_cell[cell] = record
    record = by_cell.get(goal)
    if record is None:
        raise ValueError(f"barrierless classification references missing cell {goal}")
    topology = record.topology
    if not (
        topology.get("valid_typed_identity") is True
        and tuple(topology.get("basin", ())) == a2a.PRODUCT_BASIN
        and a2.endpoint_identity(record.cluster, attacker_index)
        == a2.endpoint_identity(product, attacker_index)
    ):
        raise RuntimeError(
            "classified product cell is not the exact typed hydrolyzed-product basin"
        )
    return record


def run_barrierless_downhill_release(
    run_dir: Path,
    records: list[CellRecord],
    classification: dict[str, Any],
    product: Cluster,
    settings: DftSettings,
    *,
    attacker_index: int = OW_INDEX,
    max_steps: int = DEFAULT_DOWNHILL_MAX_STEPS,
    fmax_ev_a: float = DEFAULT_FMAX_EV_A,
    optimize_fn: OptimizeFunction = a2.optimize_minimum,
    checkpoint_energy_fn: EnergyCheckpointFunction = a2.checkpoint_energy,
    checkpoint_frequency_fn: FrequencyCheckpointFunction = a2.checkpoint_frequency,
) -> dict[str, Any]:
    """Release both scan constraints and prove downhill convergence into exact P."""
    if max_steps <= 0:
        raise ValueError("downhill max-steps must be positive")
    if not math.isfinite(fmax_ev_a) or fmax_ev_a <= 0.0:
        raise ValueError("downhill fmax must be finite and positive")
    seed = select_barrierless_release_seed(
        records,
        classification,
        product,
        attacker_index=attacker_index,
    )
    release_dir = run_dir / DOWNHILL_RELEASE_RELATIVE_PATH
    release_dir.mkdir(parents=True, exist_ok=True)
    geometry_path = release_dir / "released-product.xyz"
    trajectory_path = release_dir / "released-product.traj"
    energy_path = release_dir / "released-product.electronic-energy.json"
    frequency_path = release_dir / "released-product.fd.frequency.json"
    receipt_path = release_dir / "release-receipt.json"
    identity = _json_stable(
        {
            "release_version": DOWNHILL_RELEASE_VERSION,
            "stage": "fresh-unconstrained-downhill-release",
            "algorithm": "ase-bfgs-v1",
            "classification_fingerprint": _payload_fingerprint(classification),
            "seed_cell": [seed.row, seed.column],
            "seed_geometry_fingerprint": seed.output_geometry_fingerprint,
            "product_geometry_fingerprint": frequency_geometry_fingerprint(product),
            "settings_fingerprint": frequency_settings_fingerprint(settings),
            "constraints_released": [
                [SI_INDEX, BRIDGE_INDEX],
                [HW_INDEX, BRIDGE_INDEX],
            ],
            "fmax_ev_a": fmax_ev_a,
            "max_steps": max_steps,
        }
    )
    released = a2.checkpoint_cluster(
        geometry_path,
        seed.cluster,
        lambda: optimize_fn(
            seed.cluster,
            settings,
            max_steps=max_steps,
            trajectory=trajectory_path,
            fmax_ev_a=fmax_ev_a,
        ),
        identity=identity,
    )
    expected_identity = a2.endpoint_identity(product, attacker_index)
    actual_identity = a2.endpoint_identity(released, attacker_index)
    if actual_identity != expected_identity:
        raise RuntimeError(
            "unconstrained downhill release did not reach "
            "exact typed hydrolyzed product"
        )
    frequency = checkpoint_frequency_fn(
        frequency_path,
        released,
        settings,
        finite_difference=True,
    )
    imaginary_cm = np.asarray(frequency.imaginary_cm, dtype=float)
    if imaginary_cm.size:
        raise RuntimeError(
            "unconstrained downhill release did not reach an index-zero minimum: "
            f"{imaginary_cm.size} imaginary mode(s)"
        )
    final_energy = checkpoint_energy_fn(
        energy_path,
        released,
        settings,
        a2.R2SCAN3C_METHOD,
    )
    if not math.isfinite(final_energy):
        raise RuntimeError("unconstrained downhill release has non-finite energy")
    delta_hartree = final_energy - seed.electronic_hartree
    if delta_hartree > DOWNHILL_ENERGY_TOLERANCE_HARTREE:
        raise RuntimeError(
            "unconstrained release is not downhill: "
            f"delta={delta_hartree * HARTREE_TO_KJ:.6f} kJ/mol"
        )
    displacement = float(np.linalg.norm(released.coords - seed.cluster.coords))
    if (
        not math.isfinite(displacement)
        or displacement <= DOWNHILL_MINIMUM_DISPLACEMENT_A
    ):
        raise RuntimeError(
            "unconstrained downhill release did not make a nonzero Cartesian step"
        )
    receipt = _json_stable(
        {
            "status": "completed",
            **identity,
            "seed_electronic_hartree": seed.electronic_hartree,
            "released_electronic_hartree": final_energy,
            "electronic_delta_hartree": delta_hartree,
            "electronic_delta_kj_mol": delta_hartree * HARTREE_TO_KJ,
            "cartesian_displacement_a": displacement,
            "minimum_cartesian_displacement_a": DOWNHILL_MINIMUM_DISPLACEMENT_A,
            "typed_product_identity_matches": True,
            "zero_index_minimum": True,
            "imaginary_mode_count": int(imaginary_cm.size),
            "typed_product_identity": a2a.typed_identity_payload(
                released, attacker_index
            ),
            "released_geometry_fingerprint": frequency_geometry_fingerprint(released),
            "released_geometry_sha256": a2.sha256_path(geometry_path),
            "released_energy_sha256": a2.sha256_path(energy_path),
            "released_frequency_sha256": a2.sha256_path(frequency_path),
        }
    )
    if receipt_path.exists():
        _require_exact_json(receipt_path, receipt, "downhill release receipt")
    else:
        atomic_json(receipt_path, receipt)
    return receipt


def _energy_matrix(records: list[CellRecord], axis_points: int) -> np.ndarray:
    if len(records) != axis_points * axis_points:
        raise ValueError("classification requires every grid cell")
    matrix = np.full((axis_points, axis_points), np.nan)
    for record in records:
        if math.isfinite(matrix[record.row, record.column]):
            raise ValueError("duplicate grid cell record")
        matrix[record.row, record.column] = record.electronic_hartree
    if not np.all(np.isfinite(matrix)):
        raise ValueError("classification requires a finite complete grid")
    return (matrix - matrix[0, 0]) * HARTREE_TO_KJ


def run(args: argparse.Namespace) -> int:
    output_dir = scan_root(args.run_dir.resolve())
    output_dir.mkdir(parents=True, exist_ok=True)
    a2a.preflight_gpu_contraction_engine(args.gpu)
    intermediate = a2a.load_reference(
        args.intermediate_reference.resolve(), role="intermediate"
    )
    product = a2a.load_reference(
        args.product_reference.resolve(), role="product", template=intermediate
    )
    a2a.require_basin(intermediate, OW_INDEX, a2a.ASSOCIATIVE_BASIN)
    a2a.require_basin(product, OW_INDEX, a2a.PRODUCT_BASIN)
    settings, _, _ = a2.settings(use_gpu=args.gpu)
    records = run_scan_grid(
        output_dir,
        intermediate,
        product,
        settings,
        axis_points=args.axis_points,
        fmax_ev_a=args.fmax,
        max_steps=args.max_steps,
        optimizer_maxstep_a=args.optimizer_maxstep,
        distance_tolerance_a=args.distance_tolerance,
    )
    relative_energies = _energy_matrix(records, args.axis_points)
    classification = classify_complete_grid(
        relative_energies,
        product_index=args.axis_points - 2,
        barrier_threshold_kj_mol=args.barrier_threshold_kj_mol,
        cell_topologies={
            (record.row, record.column): record.topology for record in records
        },
    )
    payload_data: dict[str, Any] = {
        "status": "complete",
        "scan_version": SCAN_VERSION,
        "cell_count": len(records),
        "relative_energies_kj_mol": relative_energies.tolist(),
        "classification": classification,
        "coarse_grid_point_is_verified_saddle": False,
    }
    if classification["outcome"] == BARRIERLESS_SHELF:
        payload_data["downhill_release"] = run_barrierless_downhill_release(
            args.run_dir.resolve(),
            records,
            classification,
            product,
            settings,
            max_steps=args.downhill_max_steps,
            fmax_ev_a=args.fmax,
        )
    payload = _json_stable(payload_data)
    atomic_json(output_dir / "complete-grid-classification.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--intermediate-reference", type=Path, required=True)
    result.add_argument("--product-reference", type=Path, required=True)
    result.add_argument("--run-dir", type=Path, required=True)
    result.add_argument("--axis-points", type=int, default=DEFAULT_AXIS_POINTS)
    result.add_argument("--fmax", type=float, default=DEFAULT_FMAX_EV_A)
    result.add_argument("--max-steps", type=int, default=DEFAULT_MAX_STEPS)
    result.add_argument(
        "--downhill-max-steps", type=int, default=DEFAULT_DOWNHILL_MAX_STEPS
    )
    result.add_argument(
        "--optimizer-maxstep", type=float, default=DEFAULT_OPTIMIZER_MAXSTEP_A
    )
    result.add_argument(
        "--distance-tolerance", type=float, default=DEFAULT_DISTANCE_TOLERANCE_A
    )
    result.add_argument(
        "--barrier-threshold-kj-mol",
        type=float,
        default=DEFAULT_BARRIER_THRESHOLD_KJ_MOL,
    )
    result.add_argument("--gpu", action="store_true")
    result.add_argument("--gpu-mem-gb", type=float, default=16.0)
    result.add_argument("--threads", type=int, default=16)
    result.add_argument("--nice", type=int, default=10)
    result.add_argument("--log")
    return result


def execute_with_status(args: argparse.Namespace) -> int:
    output_dir = scan_root(args.run_dir.resolve())
    output_dir.mkdir(parents=True, exist_ok=True)
    status_path = output_dir / "scan-status.json"
    result_path = output_dir / "complete-grid-classification.json"
    result_path.unlink(missing_ok=True)
    atomic_json(
        status_path,
        {
            "status": "running",
            "scan_version": SCAN_VERSION,
            "started_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        },
    )
    try:
        code = run(args)
    except Exception as exc:
        atomic_json(
            status_path,
            {
                "status": "failed",
                "scan_version": SCAN_VERSION,
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
            "scan_version": SCAN_VERSION,
            "completed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "classification_sha256": a2.sha256_path(result_path),
        },
    )
    return code


def main() -> int:
    args = parser().parse_args()
    if args.axis_points < 3:
        raise ValueError("axis-points must be at least 3")
    if args.max_steps <= 0:
        raise ValueError("max-steps must be positive")
    if args.downhill_max_steps <= 0:
        raise ValueError("downhill-max-steps must be positive")
    return execute_with_status(args)


if __name__ == "__main__":
    raise SystemExit(main())
