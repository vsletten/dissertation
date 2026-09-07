"""CPU-only state construction and validation for CALC-005.

This module never imports a calculator.  It expands the condensed
``from_deck_cell`` Si seeds into the fully hydrolyzed state-205 reactant used by
the live adsorption/desorption reaction, then constructs the matched vacancy by
atom identity rather than distance.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from quarry.clusters import R_O_H, Cluster, silicic_acid
from quarry.crystal import (
    MIN_INTERATOMIC_A,
    AtomOrigin,
    CrystalCluster,
    DeckCell,
    Node,
    from_deck_cell,
    min_interatomic_distance,
)
from quarry.store import geometry_hash

CENTER_NODE: Node = (4, (0, 0, 0))
MAX_OWNER_DISTANCE_A = 1.25
MIN_OWNER_MARGIN_A = 0.15
ACID_OFFSET_A = 2.8
R_J_MOL_K = 8.31446261815324
R_KJ_MOL_K = R_J_MOL_K / 1000.0
EXPECTED_FORMULAS = {
    1: ("Al6H36O29Si", "Al6H38O30Si", "Al6H34O26"),
    2: ("Al6H38O36Si4", "Al6H42O38Si4", "Al6H38O34Si3"),
    3: ("Al7H43O46Si7", "Al7H49O49Si7", "Al7H45O45Si6"),
    4: ("Al7H43O52Si10", "Al7H51O56Si10", "Al7H47O52Si9"),
}


def _finite_number(name: str, value: float, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a real non-Boolean number")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    if positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def standard_state_1m_from_1bar_kj_mol(
    *,
    temperature_k: float,
    concentration_mol_m3: float = 1000.0,
    pressure_pa: float = 100000.0,
) -> float:
    """Ideal-solute ``1 bar -> 1 M`` free-energy correction."""

    temperature_k = _finite_number("temperature_k", temperature_k, positive=True)
    concentration_mol_m3 = _finite_number(
        "concentration_mol_m3", concentration_mol_m3, positive=True
    )
    pressure_pa = _finite_number("pressure_pa", pressure_pa, positive=True)
    ratio = concentration_mol_m3 * R_J_MOL_K * temperature_k / pressure_pa
    if not math.isfinite(ratio) or ratio <= 0:
        raise ValueError("standard-state concentration/pressure ratio is nonfinite")
    correction = R_KJ_MOL_K * temperature_k * math.log(ratio)
    if not math.isfinite(correction):
        raise ValueError("standard-state correction is nonfinite")
    return correction


def detailed_balance_rates(
    *,
    attachment_rate_s: float,
    detachment_free_energy_kj_mol: float,
    temperature_k: float,
    activity: float,
) -> dict[str, float]:
    """Close a reversible pair and expose the independent ratio equality."""

    attachment_rate_s = _finite_number(
        "attachment_rate_s", attachment_rate_s, positive=True
    )
    detachment_free_energy_kj_mol = _finite_number(
        "detachment_free_energy_kj_mol", detachment_free_energy_kj_mol
    )
    temperature_k = _finite_number("temperature_k", temperature_k, positive=True)
    activity = _finite_number("activity", activity, positive=True)
    exponent = detachment_free_energy_kj_mol / (R_KJ_MOL_K * temperature_k)
    if abs(exponent) > 700:
        raise ValueError("detachment_free_energy_kj_mol produces nonfinite rates")
    detachment_rate = attachment_rate_s * math.exp(-exponent)
    result = {
        "r_kj_mol_k": R_KJ_MOL_K,
        "attachment_rate_s": attachment_rate_s,
        "detachment_rate_s": detachment_rate,
        "kinetic_ratio": activity * attachment_rate_s / detachment_rate,
        "thermodynamic_ratio": activity * math.exp(exponent),
    }
    if not all(math.isfinite(value) for value in result.values()):
        raise ValueError("detailed-balance calculation produced nonfinite output")
    return result


def _origin_key(origin: AtomOrigin) -> tuple[str, Node, int | None]:
    return origin.kind, origin.node, origin.ordinal


def _nearest_oxygen_metrics(cluster: Cluster) -> tuple[float, float]:
    oxygen = [i for i, symbol in enumerate(cluster.symbols) if symbol == "O"]
    hydrogen = [i for i, symbol in enumerate(cluster.symbols) if symbol == "H"]
    max_owner = 0.0
    min_margin = float("inf")
    for h_index in hydrogen:
        distances = sorted(
            float(np.linalg.norm(cluster.coords[h_index] - cluster.coords[o_index]))
            for o_index in oxygen
        )
        if len(distances) < 2:
            raise ValueError("proton ownership requires at least two oxygens")
        max_owner = max(max_owner, distances[0])
        min_margin = min(min_margin, distances[1] - distances[0])
    return max_owner, min_margin


def _sha256_path(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _formula_counts(symbols: list[str]) -> Counter[str]:
    return Counter(symbols)


def _intended_owner_match(cluster: Cluster, origins: tuple[AtomOrigin, ...]) -> bool:
    if len(cluster.symbols) != len(origins):
        return False
    oxygen = [i for i, symbol in enumerate(cluster.symbols) if symbol == "O"]
    if not oxygen or not np.all(np.isfinite(cluster.coords)):
        return False
    for h_index, (symbol, origin) in enumerate(
        zip(cluster.symbols, origins, strict=True)
    ):
        if symbol != "H":
            continue
        if origin.kind in {"termination", "hydrolysis-water-h-support"}:
            expected = ("deck", origin.node)
        elif origin.kind == "hydrolysis-water-h-ligand":
            expected = ("hydrolysis-water-o", origin.node)
        else:
            return False
        owner = min(
            oxygen,
            key=lambda index: float(
                np.linalg.norm(cluster.coords[h_index] - cluster.coords[index])
            ),
        )
        actual_origin = origins[owner]
        if (actual_origin.kind, actual_origin.node) != expected:
            return False
    return True


def _cluster_exact(left: Cluster, right: Cluster) -> bool:
    return (
        left.name == right.name
        and left.symbols == right.symbols
        and np.array_equal(left.coords, right.coords)
        and left.charge == right.charge
        and left.spin == right.spin
        and left.frozen_indices == right.frozen_indices
        and left.site_family == right.site_family
        and left.note == right.note
    )


def _crystal_cluster_exact(left: CrystalCluster, right: CrystalCluster) -> bool:
    return (
        _cluster_exact(left.cluster, right.cluster)
        and left.attacked_index == right.attacked_index
        and left.bridge_index == right.bridge_index
        and left.site_kind == right.site_kind
        and left.center_site == right.center_site
        and left.metal_shells == right.metal_shells
        and left.n_intact_requested == right.n_intact_requested
        and left.n_intact == right.n_intact
        and left.atom_shell == right.atom_shell
        and left.atom_origins == right.atom_origins
        and left.center_bridges == right.center_bridges
        and left.kept_center_bridges == right.kept_center_bridges
        and left.termination_log == right.termination_log
        and left.metadata() == right.metadata()
    )


@dataclass(frozen=True)
class Calc005Pair:
    """One graph-bound, fully hydrolyzed 200s-Si environment pair."""

    condensed: CrystalCluster
    occupied: Cluster
    vacancy: Cluster
    silicic_acid: Cluster
    deck_path: Path
    occupied_origins: tuple[AtomOrigin, ...]
    vacancy_origins: tuple[AtomOrigin, ...]
    silicic_acid_origins: tuple[AtomOrigin, ...]
    x: int
    y: int
    environment_index: int
    hydrolysis_water_count: int
    occupied_states: dict[str, object]
    vacancy_states: dict[str, object]

    def validate(self) -> dict[str, bool | int | float]:
        provenance_lengths = (
            len(self.condensed.cluster.symbols) == len(self.condensed.atom_origins)
            and len(self.occupied.symbols) == len(self.occupied_origins)
            and len(self.vacancy.symbols) == len(self.vacancy_origins)
            and len(self.silicic_acid.symbols) == len(self.silicic_acid_origins)
        )
        occupied_keys = Counter(_origin_key(origin) for origin in self.occupied_origins)
        product_keys = Counter(
            _origin_key(origin)
            for origin in (*self.vacancy_origins, *self.silicic_acid_origins)
        )
        frozen_v = {
            _origin_key(self.vacancy_origins[index]): tuple(self.vacancy.coords[index])
            for index in self.vacancy.frozen_indices
        }
        frozen_c = {
            _origin_key(self.occupied_origins[index]): tuple(
                self.occupied.coords[index]
            )
            for index in self.occupied.frozen_indices
        }
        acid_keys = Counter(_origin_key(origin) for origin in self.silicic_acid_origins)
        occupied_support = occupied_keys - acid_keys
        vacancy_keys = Counter(_origin_key(origin) for origin in self.vacancy_origins)
        formula_balanced = _formula_counts(self.occupied.symbols) == (
            _formula_counts(self.vacancy.symbols)
            + _formula_counts(self.silicic_acid.symbols)
        )
        max_owner, min_margin = _nearest_oxygen_metrics(self.occupied)
        frozen = set(self.occupied.frozen_indices)
        center = self.occupied.coords[
            next(
                i
                for i, origin in enumerate(self.occupied_origins)
                if origin.kind == "deck" and origin.node == CENTER_NODE
            )
        ]
        frozen_radius = [
            float(np.linalg.norm(self.occupied.coords[i] - center)) for i in frozen
        ]
        free_metal_radius = [
            float(np.linalg.norm(self.occupied.coords[i] - center))
            for i, symbol in enumerate(self.occupied.symbols)
            if i not in frozen and symbol in {"Al", "Si"}
        ]
        hydrolysis_origins = [
            origin
            for origin in self.occupied_origins
            if origin.kind.startswith("hydrolysis-water")
        ]
        acid_oxygen_count = sum(
            1 for symbol in self.silicic_acid.symbols if symbol == "O"
        )
        acid_hydrogen_count = sum(
            1 for symbol in self.silicic_acid.symbols if symbol == "H"
        )
        expected_hydrolysis = Counter()
        for node in self.condensed.kept_center_bridges:
            expected_hydrolysis.update(
                {
                    ("hydrolysis-water-o", node, 0): 1,
                    ("hydrolysis-water-h-ligand", node, 1): 1,
                    ("hydrolysis-water-h-support", node, 2): 1,
                }
            )
        actual_hydrolysis = Counter(
            _origin_key(origin) for origin in hydrolysis_origins
        )
        expected_occupied_states = {
            "center": 205,
            "osa": 402,
            "oss": [302] * self.y + [303] * (3 - self.y),
        }
        expected_vacancy_states = {
            "center": 200,
            "osa": 404,
            "oss": [303] * self.y + [300] * (3 - self.y),
        }
        seed_formula, occupied_formula, vacancy_formula = EXPECTED_FORMULAS[
            self.environment_index
        ]
        finite_coordinates = all(
            np.all(np.isfinite(cluster.coords))
            for cluster in (
                self.condensed.cluster,
                self.occupied,
                self.vacancy,
                self.silicic_acid,
            )
        )
        canonical = _materialize_calc005_pair(
            self.deck_path, environment_index=self.environment_index
        )
        canonical_rebuild_match = (
            self.deck_path == canonical.deck_path
            and self.x == canonical.x
            and self.y == canonical.y
            and self.hydrolysis_water_count == canonical.hydrolysis_water_count
            and self.occupied_states == canonical.occupied_states
            and self.vacancy_states == canonical.vacancy_states
            and _crystal_cluster_exact(self.condensed, canonical.condensed)
            and _cluster_exact(self.occupied, canonical.occupied)
            and _cluster_exact(self.vacancy, canonical.vacancy)
            and _cluster_exact(self.silicic_acid, canonical.silicic_acid)
            and self.occupied_origins == canonical.occupied_origins
            and self.vacancy_origins == canonical.vacancy_origins
            and self.silicic_acid_origins == canonical.silicic_acid_origins
        )
        return {
            "canonical_rebuild_match": canonical_rebuild_match,
            "provenance_lengths": provenance_lengths,
            "atom_conservation": occupied_keys == product_keys and formula_balanced,
            "support_origin_match": occupied_support == vacancy_keys,
            "frozen_origin_match": set(frozen_c) == set(frozen_v),
            "frozen_coordinate_match": frozen_c == frozen_v,
            "four_center_oh_groups": acid_oxygen_count == acid_hydrogen_count == 4,
            "hydrolysis_origin_count": len(hydrolysis_origins),
            "hydrolysis_origin_match": actual_hydrolysis == expected_hydrolysis,
            "intended_owner_match": all(
                (
                    _intended_owner_match(self.occupied, self.occupied_origins),
                    _intended_owner_match(self.vacancy, self.vacancy_origins),
                    _intended_owner_match(self.silicic_acid, self.silicic_acid_origins),
                )
            ),
            "exact_state_match": self.x == 1
            and self.y == self.environment_index - 1
            and self.hydrolysis_water_count == self.environment_index
            and self.occupied_states == expected_occupied_states
            and self.vacancy_states == expected_vacancy_states,
            "exact_formula_match": self.condensed.cluster.formula == seed_formula
            and self.occupied.formula == occupied_formula
            and self.vacancy.formula == vacancy_formula
            and self.silicic_acid.formula == "H4O4Si",
            "charge_spin_match": all(
                cluster.charge == 0 and cluster.spin == 0
                for cluster in (
                    self.condensed.cluster,
                    self.occupied,
                    self.vacancy,
                    self.silicic_acid,
                )
            ),
            "finite_coordinates": bool(finite_coordinates),
            "min_interatomic_distance_a": min_interatomic_distance(
                self.occupied.coords
            ),
            "max_nearest_oh_a": max_owner,
            "min_owner_margin_a": min_margin,
            "frozen_heavy_only": all(self.occupied.symbols[i] != "H" for i in frozen),
            "frozen_peripheral": bool(frozen_radius)
            and bool(free_metal_radius)
            and float(np.mean(frozen_radius)) > float(np.mean(free_metal_radius)),
        }


def _support_hydrogen(
    oxygen: np.ndarray,
    center: np.ndarray,
    existing: np.ndarray,
    oxygen_index: int,
) -> np.ndarray:
    """Place the hydrolysis proton toward the vacated center, collision-aware."""

    primary = center - oxygen
    primary /= np.linalg.norm(primary)
    golden = np.pi * (3.0 - np.sqrt(5.0))
    directions = [primary]
    for sample in range(96):
        z = 1.0 - 2.0 * (sample + 0.5) / 96
        radius = np.sqrt(max(0.0, 1.0 - z * z))
        phi = golden * sample
        directions.append(np.array([radius * np.cos(phi), radius * np.sin(phi), z]))
    candidates = [oxygen + R_O_H * direction for direction in directions]
    other = np.delete(existing, oxygen_index, axis=0)
    return max(
        candidates,
        key=lambda point: (
            float(np.min(np.linalg.norm(other - point, axis=1))),
            float(np.dot((point - oxygen) / R_O_H, primary)),
        ),
    )


def _acid_origins(
    condensed: CrystalCluster,
) -> tuple[AtomOrigin, ...]:
    terminal = sorted(
        set(condensed.center_bridges) - set(condensed.kept_center_bridges)
    )
    retained = sorted(condensed.kept_center_bridges)
    ligand_origins: list[tuple[AtomOrigin, AtomOrigin]] = []
    for node in terminal:
        h_origin = next(
            origin
            for symbol, origin in zip(
                condensed.cluster.symbols, condensed.atom_origins, strict=True
            )
            if symbol == "H" and origin.kind == "termination" and origin.node == node
        )
        ligand_origins.append((AtomOrigin("deck", node), h_origin))
    for node in retained:
        ligand_origins.append(
            (
                AtomOrigin("hydrolysis-water-o", node, 0),
                AtomOrigin("hydrolysis-water-h-ligand", node, 1),
            )
        )
    if len(ligand_origins) != 4:
        raise RuntimeError("CALC-005 center must resolve to four Si-OH ligands")
    origins = [AtomOrigin("deck", CENTER_NODE)]
    for oxygen, hydrogen in ligand_origins:
        origins.extend([oxygen, hydrogen])
    return tuple(origins)


def _materialize_calc005_pair(
    deck_path: str | Path,
    *,
    environment_index: int,
    metal_shells: int = 2,
) -> Calc005Pair:
    """Build the deterministic ``x=1, y=i-1`` live-state pair for ``i=1..4``."""

    if type(environment_index) is not int:
        raise TypeError("environment_index must be an exact int")
    if environment_index not in {1, 2, 3, 4}:
        raise ValueError("environment_index must be one of 1, 2, 3, 4")
    condensed = from_deck_cell(
        deck_path,
        "Si",
        center_index=CENTER_NODE[0],
        metal_shells=metal_shells,
        n_intact=environment_index,
        target_charge=0,
        name=f"calc005-condensed-i{environment_index}",
    )
    if condensed.n_intact != environment_index:
        raise ValueError("environment_index did not resolve exactly")
    if condensed.kept_center_bridges[0][0] != 8:
        raise ValueError("x=1 series requires retained Osa node 8")

    center = np.asarray(condensed.cluster.coords[condensed.attacked_index])
    terminal = set(condensed.center_bridges) - set(condensed.kept_center_bridges)
    remove = {
        index
        for index, origin in enumerate(condensed.atom_origins)
        if index == condensed.attacked_index or origin.node in terminal
    }
    keep = [i for i in range(len(condensed.cluster.symbols)) if i not in remove]
    old_to_new = {old: new for new, old in enumerate(keep)}
    vacancy_symbols = [condensed.cluster.symbols[i] for i in keep]
    vacancy_coords = [np.asarray(condensed.cluster.coords[i]) for i in keep]
    vacancy_origins = [condensed.atom_origins[i] for i in keep]
    vacancy_frozen = [
        old_to_new[i] for i in condensed.cluster.frozen_indices if i in old_to_new
    ]

    for node in sorted(condensed.kept_center_bridges):
        oxygen_index = next(
            i
            for i, origin in enumerate(vacancy_origins)
            if origin.kind == "deck" and origin.node == node
        )
        oxygen = np.asarray(vacancy_coords[oxygen_index])
        proton = _support_hydrogen(
            oxygen,
            center,
            np.asarray(vacancy_coords),
            oxygen_index,
        )
        vacancy_symbols.append("H")
        vacancy_coords.append(proton)
        vacancy_origins.append(AtomOrigin("hydrolysis-water-h-support", node, 2))

    vacancy = Cluster(
        name=f"calc005-v-i{environment_index}-x1-y{environment_index - 1}",
        symbols=vacancy_symbols,
        coords=np.asarray(vacancy_coords),
        charge=0,
        spin=0,
        frozen_indices=vacancy_frozen,
        site_family=condensed.cluster.site_family,
        note="post-desorb-si matched support state",
    )

    acid = silicic_acid()
    cell = DeckCell.from_deck(deck_path)
    osa = cell.cart(condensed.kept_center_bridges[0])
    outward = center - osa
    outward /= np.linalg.norm(outward)
    y = environment_index - 1
    acid_offset = ACID_OFFSET_A + 0.6 * y
    acid = Cluster(
        name="calc005-silicic-acid",
        symbols=list(acid.symbols),
        coords=np.asarray(acid.coords) + center + acid_offset * outward,
        charge=0,
        spin=0,
        site_family=acid.site_family,
        note="neutral 1 M Si(OH)4 reference / occupied center fragment",
    )
    acid_origins = _acid_origins(condensed)

    occupied = Cluster(
        name=f"calc005-c-i{environment_index}-x1-y{environment_index - 1}",
        symbols=[*vacancy.symbols, *acid.symbols],
        coords=np.vstack((vacancy.coords, acid.coords)),
        charge=0,
        spin=0,
        frozen_indices=list(vacancy.frozen_indices),
        site_family=condensed.cluster.site_family,
        note="live fully hydrolyzed center state 205",
    )
    occupied_origins = tuple([*vacancy_origins, *acid_origins])
    pair = Calc005Pair(
        condensed=condensed,
        occupied=occupied,
        vacancy=vacancy,
        silicic_acid=acid,
        deck_path=Path(deck_path).resolve(),
        occupied_origins=occupied_origins,
        vacancy_origins=tuple(vacancy_origins),
        silicic_acid_origins=acid_origins,
        x=1,
        y=y,
        environment_index=environment_index,
        hydrolysis_water_count=environment_index,
        occupied_states={
            "center": 205,
            "osa": 402,
            "oss": [302] * y + [303] * (3 - y),
        },
        vacancy_states={
            "center": 200,
            "osa": 404,
            "oss": [303] * y + [300] * (3 - y),
        },
    )
    return pair


def build_calc005_pair(
    deck_path: str | Path,
    *,
    environment_index: int,
    metal_shells: int = 2,
) -> Calc005Pair:
    """Build and fail-closed validate one ``x=1, y=i-1`` live-state pair."""

    pair = _materialize_calc005_pair(
        deck_path,
        environment_index=environment_index,
        metal_shells=metal_shells,
    )
    report = pair.validate()
    failed = [
        key
        for key in (
            "canonical_rebuild_match",
            "provenance_lengths",
            "atom_conservation",
            "support_origin_match",
            "frozen_origin_match",
            "frozen_coordinate_match",
            "four_center_oh_groups",
            "hydrolysis_origin_match",
            "intended_owner_match",
            "exact_state_match",
            "exact_formula_match",
            "charge_spin_match",
            "finite_coordinates",
            "frozen_heavy_only",
            "frozen_peripheral",
        )
        if not report[key]
    ]
    if float(report["min_interatomic_distance_a"]) < MIN_INTERATOMIC_A:
        failed.append("min_interatomic_distance_a")
    if float(report["max_nearest_oh_a"]) > MAX_OWNER_DISTANCE_A:
        failed.append("max_nearest_oh_a")
    if float(report["min_owner_margin_a"]) < MIN_OWNER_MARGIN_A:
        failed.append("min_owner_margin_a")
    if failed:
        raise RuntimeError(f"CALC-005 pair failed CPU gates: {', '.join(failed)}")
    return pair


def probe_calc005_pairs(deck_path: str | Path) -> dict[str, object]:
    """Return deterministic JSON-ready evidence for all x=1 analytical rungs."""

    deck = Path(deck_path)

    def one_pass() -> list[dict[str, object]]:
        rows = []
        for index in range(1, 5):
            pair = build_calc005_pair(deck, environment_index=index)
            rows.append(
                {
                    "environment_index": index,
                    "x": pair.x,
                    "y": pair.y,
                    "condensed_formula": pair.condensed.cluster.formula,
                    "live_c_formula": pair.occupied.formula,
                    "matched_v_formula": pair.vacancy.formula,
                    "silicic_acid_formula": pair.silicic_acid.formula,
                    "occupied_geometry_hash": geometry_hash(pair.occupied.to_xyz()),
                    "vacancy_geometry_hash": geometry_hash(pair.vacancy.to_xyz()),
                    "atom_map_sha256": hashlib.sha256(
                        json.dumps(
                            [_origin_key(origin) for origin in pair.occupied_origins],
                            separators=(",", ":"),
                        ).encode()
                    ).hexdigest(),
                    "occupied_states": pair.occupied_states,
                    "vacancy_states": pair.vacancy_states,
                    **pair.validate(),
                }
            )
        return rows

    first = one_pass()
    second = one_pass()
    if first != second:
        raise RuntimeError("CALC-005 probe is not deterministic across two passes")
    occupied_hashes = {str(row["occupied_geometry_hash"]) for row in first}
    vacancy_hashes = {str(row["vacancy_geometry_hash"]) for row in first}
    if len(occupied_hashes) != 4 or len(vacancy_hashes) != 4:
        raise RuntimeError("CALC-005 live C/V geometries are not distinct by rung")
    return {
        "schema": "calc005-live-pair-probe-v1",
        "calc005_sha256": _sha256_path(Path(__file__)),
        "crystal_sha256": _sha256_path(Path(__file__).with_name("crystal.py")),
        "deck": str(deck),
        "deck_sha256": _sha256_path(deck),
        "center_site": 4,
        "metal_shells": 2,
        "topology": "x=1,y=i-1; Osa.sih(402)->Osa.albr(404)",
        "two_pass_deterministic": True,
        "pairs": first,
    }
