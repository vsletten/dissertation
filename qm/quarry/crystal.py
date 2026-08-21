"""Crystallographic cluster builder — Phase 2 (HANDOFF.md §2a).

Cuts terminated, peripherally-frozen clusters from the mineral unit cell
carried by a petra deck. The deck's ``[cell]`` section is the single
source of crystallographic truth: fractional sites, the frac→Cartesian
matrix (lattice vectors as *columns*, petra-core crystal.rs convention),
and the KMC bond graph with periodic ``dcell`` image offsets. The graph
carries full coordination (kaolinite: Si degree 4, Al degree 6, Osa is
the real 3-coordinate apical oxygen, Oaa the structural hydroxyl), so
termination and charge bookkeeping follow from Pauling bond valence
alone — no chemistry tables beyond Si–O = 1.0, Al–O = 0.5, O wants 2.

Termination rules (SURVEY.md §6.2, Pelmenschikov lattice resistance):
every bridge oxygen adjacent to an included metal is kept at its
crystallographic position; oxygens missing metal partners are capped
with H per bond-valence deficit, then protons are added/removed at the
most under/over-saturated peripheral oxygens until the target charge is
met. The outermost metal shell (and everything exclusively attached to
it) is frozen — lattice resistance is physics, not a detail.
"""

from __future__ import annotations

import math
import tomllib
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from quarry.clusters import R_O_H, Cluster, SiteFamily, _unit, merge

KIND_TO_FAMILY = {
    "Al": SiteFamily.AL,
    "Si": SiteFamily.SI,
    "Oss": SiteFamily.SI_O_SI,
    "Osa": SiteFamily.SI_O_AL2,
    "Oaa": SiteFamily.AL_OH_AL,
}
KIND_ELEMENT = {"Al": "Al", "Si": "Si", "Oss": "O", "Osa": "O", "Oaa": "O"}
METAL_KINDS = frozenset({"Al", "Si"})
# Pauling bond strength per M-O bond: cation valence / coordination.
PAULING_STRENGTH = {"Si": 1.0, "Al": 0.5}
FORMAL_CHARGE = {"Si": 4, "Al": 3}
# Structural protons the mineral itself carries (Al-OH-Al is a hydroxyl).
STRUCTURAL_H = {"Oaa": 1}
O_VALENCE = 2.0
MIN_INTERATOMIC_A = 0.75
# Determinant magnitude below which lattice vectors are treated as degenerate
# (zero cell volume). Absolute tolerance on |det|; scale conventions vary, so
# keep it a named constant for easy tuning.
DEGENERATE_DET_TOL = 1e-12

# A node is one site in one periodic image.
Node = tuple[int, tuple[int, int, int]]


@dataclass(frozen=True)
class CellSite:
    kind: str
    frac: tuple[float, float, float]


@dataclass(frozen=True)
class CellBond:
    i: int
    j: int
    dcell: tuple[int, int, int]


@dataclass
class DeckCell:
    """The ``[cell]`` section of a petra deck, with a lazy periodic graph."""

    matrix: np.ndarray  # 3x3, columns are the lattice vectors
    sites: list[CellSite]
    bonds: list[CellBond]
    _adjacency: dict[int, list[tuple[int, tuple[int, int, int]]]] = field(
        default_factory=dict, repr=False
    )

    def __post_init__(self) -> None:
        matrix = np.asarray(self.matrix, dtype=float)
        if matrix.shape != (3, 3):
            raise ValueError(f"deck cell matrix must be 3x3, got {matrix.shape}")
        if not np.all(np.isfinite(matrix)):
            raise ValueError("deck cell matrix must be finite")
        if abs(float(np.linalg.det(matrix))) < DEGENERATE_DET_TOL:
            raise ValueError(
                "deck cell lattice vectors are degenerate (zero cell volume)"
            )
        self.matrix = matrix
        unknown = {s.kind for s in self.sites} - set(KIND_TO_FAMILY)
        if unknown:
            raise ValueError(f"unknown site kinds in deck cell: {sorted(unknown)}")
        adj: dict[int, list[tuple[int, tuple[int, int, int]]]] = {
            i: [] for i in range(len(self.sites))
        }
        for b in self.bonds:
            adj[b.i].append((b.j, b.dcell))
            adj[b.j].append((b.i, tuple(-d for d in b.dcell)))
        self._adjacency = adj

    @classmethod
    def from_deck(cls, path: str | Path) -> DeckCell:
        with open(path, "rb") as fh:
            deck = tomllib.load(fh)
        cell = deck.get("cell")
        if cell is None:
            raise ValueError(f"{path}: deck has no [cell] section")
        sites = [
            CellSite(kind=s["kind"], frac=tuple(float(x) for x in s["frac"]))
            for s in cell["sites"]
        ]
        bonds = [
            CellBond(
                i=int(b["i"]),
                j=int(b["j"]),
                dcell=tuple(int(d) for d in b.get("dcell", [0, 0, 0])),
            )
            for b in cell.get("bonds", [])
        ]
        matrix = np.asarray(cell["matrix"], dtype=float)
        return cls(matrix=matrix, sites=sites, bonds=bonds)

    def kind(self, node: Node) -> str:
        return self.sites[node[0]].kind

    def cart(self, node: Node) -> np.ndarray:
        site, image = node
        frac = np.asarray(self.sites[site].frac) + np.asarray(image, dtype=float)
        return self.matrix @ frac

    def neighbors(self, node: Node) -> list[Node]:
        site, image = node
        out = []
        for j, dcell in self._adjacency[site]:
            out.append((j, tuple(a + d for a, d in zip(image, dcell, strict=True))))
        return sorted(out)

    def degree(self, site: int) -> int:
        return len(self._adjacency[site])

    def first_site_of_kind(self, kind: str) -> int:
        for i, s in enumerate(self.sites):
            if s.kind == kind:
                return i
        raise ValueError(f"deck cell has no site of kind '{kind}'")


@dataclass
class CrystalCluster:
    """A built cluster plus the bookkeeping the campaign driver needs."""

    cluster: Cluster
    site_kind: str
    center_site: int
    metal_shells: int
    # Atom index of the center bridge oxygen (None for metal-family centers).
    bridge_index: int | None
    # Atom index of the metal the attacker approaches.
    attacked_index: int
    n_intact_requested: int | None
    n_intact: int
    # Per-atom metal-shell provenance (H atoms inherit their oxygen's shell).
    atom_shell: list[int]
    termination_log: list[str] = field(default_factory=list)

    def metadata(self) -> dict:
        c = self.cluster
        return {
            "name": c.name,
            "site_kind": self.site_kind,
            "site_family": c.site_family.name if c.site_family else None,
            "center_site": self.center_site,
            "metal_shells": self.metal_shells,
            "bridge_index": self.bridge_index,
            "attacked_index": self.attacked_index,
            "n_intact_requested": self.n_intact_requested,
            "n_intact": self.n_intact,
            "n_atoms": len(c.symbols),
            "formula": c.formula,
            "charge": c.charge,
            "n_frozen": len(c.frozen_indices),
            "termination_log": self.termination_log,
        }


def _formal_charge(symbols: list[str], n_h: int) -> int:
    q = n_h  # each proton is +1 against the O2- lattice
    for s in symbols:
        if s in FORMAL_CHARGE:
            q += FORMAL_CHARGE[s]
        elif s == "O":
            q -= 2
        elif s != "H":
            raise ValueError(f"unexpected element '{s}' in formal-charge bookkeeping")
    return q


def _h_directions(away: np.ndarray, n_h: int) -> list[np.ndarray]:
    """Unit vectors for ``n_h`` protons on an O whose open side is ``away``.

    Tilted off the away-axis (25 deg for OH, 52.25 deg each for OH2 so the
    H-O-H angle lands near 104.5 deg). Only topology matters — the
    pipeline optimizes every guess.
    """
    if not 1 <= n_h <= 2:
        raise ValueError("an oxygen carries at most two protons")
    base = _unit(away)
    ref = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(ref, base))) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])
    e1 = _unit(np.cross(base, ref))
    e2 = np.cross(base, e1)
    theta = math.radians(25.0 if n_h == 1 else 52.25)
    phis = [0.0] if n_h == 1 else [0.0, math.pi]
    return [
        math.cos(theta) * base + math.sin(theta) * (math.cos(p) * e1 + math.sin(p) * e2)
        for p in phis
    ]


def min_interatomic_distance(coords: np.ndarray) -> float:
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    iu = np.triu_indices(len(coords), k=1)
    return float(d[iu].min()) if iu[0].size else math.inf


def from_deck_cell(
    deck_path: str | Path,
    site_kind: str,
    *,
    center_index: int | None = None,
    metal_shells: int = 2,
    n_intact: int | None = None,
    target_charge: int = 0,
    name: str | None = None,
) -> CrystalCluster:
    """Cut a terminated, peripherally-frozen cluster around one site.

    ``site_kind`` picks the KMC family (Al/Si/Oss/Osa/Oaa);
    ``center_index`` a specific cell site (default: first of that kind).
    ``metal_shells`` is the cutoff in metal shells from the center;
    the outermost shell is frozen (none frozen when ``metal_shells < 2``).
    ``n_intact`` is the number of intact bridges on the attacked metal,
    center bridge included — the connectivity rung of the ladder. Pruned
    bridges have their partner metals removed ("dissolved") and become
    terminal oxygens; bridges shared with a center-bridge partner
    (edge-sharing octahedra) are unprunable and stay intact.
    """
    cell = DeckCell.from_deck(deck_path)
    if site_kind not in KIND_TO_FAMILY:
        raise ValueError(f"unknown site kind '{site_kind}'")
    center_site = (
        cell.first_site_of_kind(site_kind) if center_index is None else center_index
    )
    if cell.kind((center_site, (0, 0, 0))) != site_kind:
        raise ValueError(
            f"center site {center_site} is kind "
            f"'{cell.sites[center_site].kind}', not '{site_kind}'"
        )
    center: Node = (center_site, (0, 0, 0))
    is_bridge_center = site_kind not in METAL_KINDS
    log: list[str] = []

    # --- pick the attacked metal and the sacrosanct center partners ----
    if is_bridge_center:
        center_partners = cell.neighbors(center)
        si_partners = [m for m in center_partners if cell.kind(m) == "Si"]
        attacked_node = si_partners[0] if si_partners else center_partners[0]
    else:
        center_partners = []
        attacked_node = center

    # --- connectivity pruning on the attacked metal ---------------------
    attacked_bridges = cell.neighbors(attacked_node)
    max_n = len(attacked_bridges)
    if n_intact is not None and not 1 <= n_intact <= max_n:
        raise ValueError(
            f"n_intact must be in 1..{max_n} for a "
            f"{cell.kind(attacked_node)} center (got {n_intact})"
        )
    # A pruned bridge's partner metals are "dissolved" (excluded from the
    # cluster); the bridge oxygen itself survives as a terminal OH on the
    # attacked metal — that is Q-notation chemistry. Bridges touching a
    # center-bridge partner are unprunable, and dissolving a partner
    # cascades to every bridge that shares it (edge-sharing octahedra).
    partners_of = {
        b: frozenset(
            m
            for m in cell.neighbors(b)
            if m != attacked_node and cell.kind(m) in METAL_KINDS
        )
        for b in attacked_bridges
    }
    unprunable = {
        b
        for b in attacked_bridges
        if b == center or partners_of[b] & set(center_partners)
    }
    excluded_metals: set[Node] = set()
    if n_intact is None:
        kept = list(attacked_bridges)
    else:
        kept = sorted(unprunable)
        for b in sorted(b for b in attacked_bridges if b not in unprunable):
            if len(kept) < n_intact and not partners_of[b] & excluded_metals:
                kept.append(b)
            else:
                excluded_metals |= partners_of[b]
                log.append(f"pruned bridge {b}: partner metals dissolved")
        # Cascade: a kept bridge whose partner dissolved is pruned too.
        changed = True
        while changed:
            changed = False
            for b in list(kept):
                if b != center and partners_of[b] & excluded_metals:
                    kept.remove(b)
                    excluded_metals |= partners_of[b]
                    log.append(f"pruned bridge {b}: cascade (shared partner)")
                    changed = True
        if len(kept) != n_intact:
            log.append(
                f"DEVIATION: requested n_intact={n_intact}, delivered "
                f"{len(kept)} (unprunable/edge-sharing constraints)"
            )
    resolved_n_intact = len(kept)

    # --- metal BFS to the shell cutoff ----------------------------------
    metal_shell: dict[Node, int] = {}
    queue: deque[Node] = deque()
    if is_bridge_center:
        for m in center_partners:
            metal_shell[m] = 1
            queue.append(m)
    else:
        metal_shell[attacked_node] = 0
        queue.append(attacked_node)
    while queue:
        m = queue.popleft()
        if metal_shell[m] >= metal_shells:
            continue
        for b in sorted(cell.neighbors(m)):
            for m2 in cell.neighbors(b):
                if (
                    m2 != m
                    and cell.kind(m2) in METAL_KINDS
                    and m2 not in excluded_metals
                    and m2 not in metal_shell
                ):
                    metal_shell[m2] = metal_shell[m] + 1
                    queue.append(m2)

    # --- collect bridge oxygens adjacent to included metals -------------
    # (a pruned bridge re-enters here as a terminal oxygen: its partner
    # metals are excluded, so only the attacked metal claims it)
    bridge_metals: dict[Node, list[Node]] = {}
    for m in metal_shell:
        for b in cell.neighbors(m):
            bridge_metals.setdefault(b, []).append(m)
    if is_bridge_center and center not in bridge_metals:
        raise RuntimeError("center bridge lost during cluster construction")

    def bridge_shell(b: Node) -> int:
        if is_bridge_center and b == center:
            return 0
        return min(metal_shell[m] for m in bridge_metals[b])

    # --- assemble heavy atoms (deterministic order) ---------------------
    symbols: list[str] = []
    coords: list[np.ndarray] = []
    shells: list[int] = []
    atom_of: dict[Node, int] = {}

    def add(node: Node, shell: int) -> int:
        atom_of[node] = len(symbols)
        symbols.append(KIND_ELEMENT[cell.kind(node)])
        coords.append(cell.cart(node))
        shells.append(shell)
        return atom_of[node]

    bridge_index: int | None = None
    if is_bridge_center:
        bridge_index = add(center, 0)
    for m in sorted(metal_shell, key=lambda n: (metal_shell[n], n)):
        add(m, metal_shell[m])
    for b in sorted(bridge_metals, key=lambda n: (bridge_shell(n), n)):
        if b not in atom_of:
            add(b, bridge_shell(b))
    attacked_index = atom_of[attacked_node]

    # --- proton bookkeeping by bond valence ------------------------------
    # deficit = 2 - sum(strengths of included metals) - structural H.
    o_records = []  # (atom_idx, node, n_h, residual)
    for b, metals in sorted(bridge_metals.items()):
        strengths = sum(PAULING_STRENGTH[cell.kind(m)] for m in metals)
        structural = STRUCTURAL_H.get(cell.kind(b), 0)
        deficit = O_VALENCE - strengths - structural
        n_base = max(0, math.floor(deficit + 1e-9))
        n_h = structural + n_base
        residual = deficit - n_base
        if len(metals) < cell.degree(b[0]):
            log.append(
                f"terminal O {b}: {len(metals)}/{cell.degree(b[0])} metals, "
                f"{n_h} H (residual {residual:+.1f})"
            )
        o_records.append([atom_of[b], b, n_h, residual])

    n_h_total = sum(r[2] for r in o_records)
    q0 = _formal_charge(symbols, n_h_total)
    need = target_charge - q0
    center_pos = cell.cart(center)

    def far_first(rec):
        return -float(np.linalg.norm(np.asarray(coords[rec[0]]) - center_pos))

    reactive = {attacked_index, bridge_index}
    if need > 0:
        candidates = sorted(
            (r for r in o_records if r[3] >= 0.5 - 1e-9 and r[2] < 2), key=far_first
        )
        candidates = [r for r in candidates if r[0] not in reactive]
        if len(candidates) < need:
            raise ValueError(
                f"cannot reach charge {target_charge}: need {need} extra "
                f"protons, only {len(candidates)} protonatable oxygens"
            )
        for r in candidates[:need]:
            r[2] += 1
            r[3] -= 1.0
            log.append(f"protonated O {r[1]} for charge balance")
    elif need < 0:
        candidates = sorted(
            (r for r in o_records if r[2] >= 1 and r[0] not in reactive),
            key=lambda r: (r[3], far_first(r)),
        )
        if len(candidates) < -need:
            raise ValueError(
                f"cannot reach charge {target_charge}: need {-need} "
                f"deprotonations, only {len(candidates)} candidates"
            )
        for r in candidates[:-need]:
            r[2] -= 1
            r[3] += 1.0
            log.append(f"deprotonated O {r[1]} for charge balance")

    # --- place protons ----------------------------------------------------
    for o_idx, b, n_h, _residual in o_records:
        if n_h == 0:
            continue
        o_pos = np.asarray(coords[o_idx])
        metal_dirs = [
            _unit(np.asarray(coords[atom_of[m]]) - o_pos) for m in bridge_metals[b]
        ]
        away = -np.sum(metal_dirs, axis=0)
        if float(np.linalg.norm(away)) < 1e-3:
            away = np.cross(metal_dirs[0], [1.0, 0.0, 0.0])
            if float(np.linalg.norm(away)) < 1e-3:
                away = np.cross(metal_dirs[0], [0.0, 1.0, 0.0])
        for d in _h_directions(away, n_h):
            symbols.append("H")
            coords.append(o_pos + R_O_H * d)
            shells.append(shells[o_idx])

    # --- frozen shell -----------------------------------------------------
    # Heavy atoms of the outermost metal shell freeze at crystallographic
    # positions (lattice resistance); protons never freeze — their
    # positions are constructed guesses, not crystallography.
    frozen = (
        [
            i
            for i, (s, sym) in enumerate(zip(shells, symbols, strict=True))
            if s >= metal_shells and sym != "H"
        ]
        if metal_shells >= 2
        else []
    )
    if not frozen:
        log.append("no frozen shell (metal_shells < 2) — free-cluster model")

    n_h_final = sum(1 for s in symbols if s == "H")
    q_final = _formal_charge(symbols, n_h_final)
    if q_final != target_charge:
        raise RuntimeError(
            f"charge bookkeeping failed: built {q_final}, target {target_charge}"
        )

    coords_arr = np.asarray(coords)
    dmin = min_interatomic_distance(coords_arr)
    if dmin < MIN_INTERATOMIC_A:
        raise RuntimeError(
            f"atom collision in built cluster: min distance {dmin:.2f} A"
        )
    reactive_frozen = attacked_index in frozen or (
        bridge_index is not None and bridge_index in frozen
    )
    if reactive_frozen:
        raise RuntimeError("reactive center landed in the frozen shell")

    tag = f"{site_kind.lower()}-c{center_site}-s{metal_shells}-n{resolved_n_intact}"
    cluster = Cluster(
        name=name or f"crystal-{tag}",
        symbols=symbols,
        coords=coords_arr,
        charge=target_charge,
        frozen_indices=frozen,
        site_family=KIND_TO_FAMILY[site_kind],
        note=f"cut from {Path(deck_path).name} cell",
    )
    return CrystalCluster(
        cluster=cluster,
        site_kind=site_kind,
        center_site=center_site,
        metal_shells=metal_shells,
        bridge_index=bridge_index,
        attacked_index=attacked_index,
        n_intact_requested=n_intact,
        n_intact=resolved_n_intact,
        atom_shell=shells,
        termination_log=log,
    )


def attack_complex(
    cc: CrystalCluster,
    attacker: Cluster,
    *,
    approach_a: float = 3.2,
    theta_deg: float = 80.0,
    n_azimuths: int = 24,
) -> tuple[Cluster, int]:
    """Place the attacker for hydrolysis of the center bridge.

    Flank geometry as in Phase 1 (attacker O ~80 deg off the M->Obr axis,
    first O-H aimed at the bridging O), but the azimuth is searched for
    maximum clearance from the rest of the cluster — a crystallographic
    cluster is crowded where a hand-built dimer was not. Returns the
    merged cluster and the attacker-oxygen atom index.
    """
    from quarry.clusters import _rotation_between

    c = cc.cluster
    m_pos = c.coords[cc.attacked_index]
    if cc.bridge_index is not None:
        axis = _unit(c.coords[cc.bridge_index] - m_pos)
    else:
        # Metal-family center: approach along the least-coordinated side.
        o_dirs = [
            _unit(c.coords[i] - m_pos)
            for i, s in enumerate(c.symbols)
            if s == "O" and np.linalg.norm(c.coords[i] - m_pos) < 2.3
        ]
        axis = -_unit(np.sum(o_dirs, axis=0))
    ref = np.array([1.0, 0.0, 0.0])
    if abs(float(np.dot(ref, axis))) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])
    e1 = _unit(np.cross(axis, ref))
    e2 = np.cross(axis, e1)

    # Crystallographic clusters are crowded: search azimuth first at the
    # requested polar angle, then widen the polar angle (nearest-to-
    # requested first) until a clear approach exists. The proton stays
    # deliverable to the bridge anywhere in this cone.
    if cc.bridge_index is not None:
        thetas = [theta_deg, 70.0, 90.0, 60.0, 100.0, 110.0]
    else:
        thetas = [0.0]
    best_target, best_score = None, -math.inf
    for th in thetas:
        theta = math.radians(th)
        for k in range(n_azimuths if theta > 0 else 1):
            phi = 2.0 * math.pi * k / n_azimuths
            direction = math.cos(theta) * axis + math.sin(theta) * (
                math.cos(phi) * e1 + math.sin(phi) * e2
            )
            target = m_pos + approach_a * direction
            score = float(np.min(np.linalg.norm(c.coords - target, axis=1)))
            if score > best_score:
                best_target, best_score = target, score
        if best_score >= 2.0:
            break
    assert best_target is not None
    if best_score < 1.6:
        raise RuntimeError(
            f"no attack direction with >=1.6 A clearance (best {best_score:.2f} A)"
        )

    shifted = attacker.coords - attacker.coords[0]
    if cc.bridge_index is not None and len(attacker.symbols) > 1:
        rot = _rotation_between(shifted[1], c.coords[cc.bridge_index] - best_target)
        shifted = shifted @ rot.T
    moved = Cluster(
        name=attacker.name,
        symbols=list(attacker.symbols),
        coords=shifted + best_target,
        charge=attacker.charge,
        spin=attacker.spin,
    )
    ow_index = len(c.symbols)
    return merge(c, moved, name=f"{c.name}+{attacker.name}"), ow_index
