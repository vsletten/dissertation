"""Cluster models per KMC site family and protonation state.

Phase 0/1 scope: the benchmark molecules (Xiao & Lasaga's dimers, their
monomer fragments, and the attacking species) built as clean initial
geometries — the pipeline pre-optimizes everything, so these only need
to be topologically right with sane bond lengths. The crystallographic
builder (terminated, constrained edge clusters cut from kaolinite
coordinates, SURVEY.md §6.2) is Phase 2; the ``frozen_indices``
constraint contract is already part of the dataclass so the pipeline
layer never changes.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum

import numpy as np

# Initial-guess bond lengths, Angstrom.
R_SI_O = 1.63
R_AL_O = 1.75
R_O_H = 0.96


class SiteFamily(Enum):
    """The KMC state taxonomy (legacy 100s-500s; petra kind names)."""

    AL = 100  # petra kind "Al"
    SI = 200  # petra kind "Si"
    SI_O_SI = 300  # petra kind "Oss"
    SI_O_AL2 = 400  # petra kind "Osa"
    AL_OH_AL = 500  # petra kind "Oaa"


@dataclass
class Cluster:
    """A molecular model: geometry + electronic state + constraint set."""

    name: str
    symbols: list[str]
    coords: np.ndarray  # (n, 3) Angstrom
    charge: int = 0
    spin: int = 0  # 2S, pyscf convention
    frozen_indices: list[int] = field(default_factory=list)
    site_family: SiteFamily | None = None
    note: str = ""

    def __post_init__(self) -> None:
        self.coords = np.asarray(self.coords, dtype=float).reshape(-1, 3)
        if len(self.symbols) != len(self.coords):
            raise ValueError("symbols and coords length mismatch")
        if any(not 0 <= i < len(self.symbols) for i in self.frozen_indices):
            raise ValueError("frozen index out of range")

    @property
    def formula(self) -> str:
        counts: dict[str, int] = {}
        for s in self.symbols:
            counts[s] = counts.get(s, 0) + 1
        return "".join(f"{el}{n if n > 1 else ''}" for el, n in sorted(counts.items()))

    def to_xyz(self, comment: str = "") -> str:
        lines = [str(len(self.symbols)), comment or self.name]
        lines += [
            f"{s} {x:.8f} {y:.8f} {z:.8f}"
            for s, (x, y, z) in zip(self.symbols, self.coords, strict=True)
        ]
        return "\n".join(lines)

    def to_pyscf_atom(self) -> str:
        """Geometry block for ``pyscf.gto.M(atom=...)`` (Angstrom)."""
        return "; ".join(
            f"{s} {x:.8f} {y:.8f} {z:.8f}"
            for s, (x, y, z) in zip(self.symbols, self.coords, strict=True)
        )


def _unit(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v)


def _tetrahedral_completion(axis: np.ndarray) -> list[np.ndarray]:
    """Three unit vectors completing a tetrahedron whose fourth arm is ``axis``."""
    a = _unit(axis)
    # Any vector not parallel to a.
    ref = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(ref, a)) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])
    e1 = _unit(np.cross(a, ref))
    e2 = np.cross(a, e1)
    cos_t = math.cos(math.radians(109.471))  # tetrahedral angle from `a`
    sin_t = math.sin(math.radians(109.471))
    dirs = []
    for k in range(3):
        phi = 2.0 * math.pi * k / 3.0
        dirs.append(cos_t * a + sin_t * (math.cos(phi) * e1 + math.sin(phi) * e2))
    return dirs


def _hydroxylate(
    center: np.ndarray, directions: list[np.ndarray], r_mo: float
) -> tuple[list[str], list[np.ndarray]]:
    """Attach OH groups along ``directions`` from ``center``."""
    symbols: list[str] = []
    coords: list[np.ndarray] = []
    for d in directions:
        o = center + r_mo * d
        # Tilt the O-H off the M-O axis so the guess isn't linear.
        tilt = _unit(np.cross(d, _tetrahedral_completion(d)[0]) + 0.4 * d)
        h = o + R_O_H * tilt
        symbols += ["O", "H"]
        coords += [o, h]
    return symbols, coords


def water() -> Cluster:
    """H2O at the experimental-ish geometry."""
    ang = math.radians(104.5)
    return Cluster(
        name="water",
        symbols=["O", "H", "H"],
        coords=np.array(
            [
                [0.0, 0.0, 0.0],
                [R_O_H, 0.0, 0.0],
                [R_O_H * math.cos(ang), R_O_H * math.sin(ang), 0.0],
            ]
        ),
    )


def hydronium() -> Cluster:
    """H3O+ (C3v pyramid) — the acid-pathway attacking species."""
    r = 0.98
    z = -0.25
    coords = [[0.0, 0.0, 0.0]]
    for k in range(3):
        phi = 2.0 * math.pi * k / 3.0
        coords.append([r * math.cos(phi), r * math.sin(phi), z])
    return Cluster(
        name="hydronium",
        symbols=["O", "H", "H", "H"],
        coords=np.array(coords),
        charge=1,
    )


ACID_MICROSOLVATION_WATER_COUNTS = frozenset({1, 3, 4, 5, 6})


def _validate_acid_water_count(n_water: int) -> None:
    if n_water not in ACID_MICROSOLVATION_WATER_COUNTS:
        allowed = ", ".join(
            str(value) for value in sorted(ACID_MICROSOLVATION_WATER_COUNTS)
        )
        raise ValueError(f"n_water must be one of {allowed}")


def hydronium_hydrate(n_water: int) -> Cluster:
    """H3O+ · (n-1)H2O fragment used by matched acid calculations.

    ``n_water=1`` preserves the historical one-water model. Production
    microsolvation uses 3--6 total explicit water oxygens. The extra waters
    form a deterministic first shell around hydronium; geometry optimization,
    not this seed, establishes the stationary point.
    """
    _validate_acid_water_count(n_water)
    base = hydronium()
    symbols = list(base.symbols)
    coords = list(base.coords)
    for direction in _shell_directions(n_water - 1):
        symbols += ["O", "H", "H"]
        coords += _water_at(direction, 2.75)
    return Cluster(
        name=f"hydronium-{n_water}w",
        symbols=symbols,
        coords=np.array(coords),
        charge=1,
        note=f"H3O+ with {n_water - 1} first-shell waters",
    )


def hydroxide() -> Cluster:
    """OH- — the base-pathway attacking species."""
    return Cluster(
        name="hydroxide",
        symbols=["O", "H"],
        coords=np.array([[0.0, 0.0, 0.0], [R_O_H, 0.0, 0.0]]),
        charge=-1,
    )


def silicic_acid() -> Cluster:
    """Si(OH)4 — the dissolved monomer and smoke-test workhorse."""
    si = np.zeros(3)
    corners = ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
    dirs = [_unit(np.array(d, dtype=float)) for d in corners]
    symbols, coords = _hydroxylate(si, dirs, R_SI_O)
    return Cluster(
        name="silicic-acid",
        symbols=["Si", *symbols],
        coords=np.array([si, *coords]),
        site_family=SiteFamily.SI,
    )


def _shell_directions(n: int) -> list[np.ndarray]:
    """n roughly-uniform unit vectors (golden-spiral points on a sphere)."""
    golden = math.pi * (3.0 - math.sqrt(5.0))
    dirs = []
    for k in range(n):
        z = 1.0 - 2.0 * (k + 0.5) / n
        r = math.sqrt(max(0.0, 1.0 - z * z))
        phi = golden * k
        dirs.append(np.array([r * math.cos(phi), r * math.sin(phi), z]))
    return dirs


MAX_HYDRATE_WATERS = 24


def _water_at(direction: np.ndarray, radius: float) -> list[np.ndarray]:
    """[O, H, H] for a water at ``radius`` along ``direction``, H2 pointing
    outward (dipole bisector along the shell normal) so shell waters never
    reach back toward the core."""
    d = _unit(direction)
    ref = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(ref, d)) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])
    e1 = _unit(np.cross(d, ref))
    half = math.radians(104.5 / 2.0)
    o = radius * d
    h1 = o + R_O_H * (math.cos(half) * d + math.sin(half) * e1)
    h2 = o + R_O_H * (math.cos(half) * d - math.sin(half) * e1)
    return [o, h1, h2]


def silicic_acid_hydrate(n_water: int) -> Cluster:
    """Si(OH)4 · n H2O — the Phase 0 smoke-test series (SURVEY.md §9).

    Waters are placed on a loose golden-spiral shell around the cluster;
    only the topology matters (the smoke test measures E+grad+Hessian
    cost, and real work pre-optimizes). Larger n exists purely to push
    the benchmark toward the 100+-atom regime where the GPU pays off.
    """
    if not 0 <= n_water <= MAX_HYDRATE_WATERS:
        raise ValueError(f"n_water must be in 0..{MAX_HYDRATE_WATERS}")
    base = silicic_acid()
    symbols = list(base.symbols)
    coords = list(base.coords)
    # Second shell sits further out so waters never collide.
    for d in _shell_directions(min(n_water, 12)):
        symbols += ["O", "H", "H"]
        coords += _water_at(d, 3.4)
    for d in _shell_directions(max(0, n_water - 12)):
        symbols += ["O", "H", "H"]
        coords += _water_at(d, 6.2)
    return Cluster(
        name=f"silicic-acid-{n_water}w",
        symbols=symbols,
        coords=np.array(coords),
        site_family=SiteFamily.SI,
    )


def _bridged_dimer(m2_symbol: str, r_m2_o: float, charge: int, name: str) -> Cluster:
    """(HO)3Si-O-M(OH)3 with a bridging oxygen and ~140 deg Si-O-M angle."""
    bridge = np.zeros(3)
    half = math.radians(140.0) / 2.0
    # Directions from the bridging O to the two centers, 140 deg apart.
    d_si = np.array([math.cos(half), math.sin(half), 0.0])
    d_m2 = np.array([math.cos(half), -math.sin(half), 0.0])
    si = bridge + R_SI_O * d_si
    m2 = bridge + r_m2_o * d_m2
    symbols: list[str] = ["O", "Si", m2_symbol]
    coords: list[np.ndarray] = [bridge, si, m2]
    s, c = _hydroxylate(si, _tetrahedral_completion(si - bridge), R_SI_O)
    symbols += s
    coords += c
    s, c = _hydroxylate(m2, _tetrahedral_completion(m2 - bridge), r_m2_o)
    symbols += s
    coords += c
    return Cluster(
        name=name,
        symbols=symbols,
        coords=np.array(coords),
        charge=charge,
    )


def disilicate() -> Cluster:
    """(HO)3Si-O-Si(OH)3 — the Xiao & Lasaga Si-O-Si benchmark cluster."""
    c = _bridged_dimer("Si", R_SI_O, 0, "disilicate")
    c.site_family = SiteFamily.SI_O_SI
    return c


def aluminosilicate_dimer() -> Cluster:
    """(HO)3Si-O-Al(OH)3^- — the Xiao & Lasaga Si-O-Al benchmark cluster."""
    c = _bridged_dimer("Al", R_AL_O, -1, "aluminosilicate-dimer")
    c.site_family = SiteFamily.SI_O_AL2
    return c


def merge(a: Cluster, b: Cluster, name: str | None = None) -> Cluster:
    """One cluster from two (frozen sets carried over, charges added)."""
    return Cluster(
        name=name or f"{a.name}+{b.name}",
        symbols=[*a.symbols, *b.symbols],
        coords=np.vstack([a.coords, b.coords]),
        charge=a.charge + b.charge,
        spin=a.spin + b.spin,
        frozen_indices=[
            *a.frozen_indices,
            *(i + len(a.symbols) for i in b.frozen_indices),
        ],
        site_family=a.site_family,
    )


def _rotation_between(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Rotation matrix taking unit direction of ``a`` onto that of ``b``."""
    a, b = _unit(a), _unit(b)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    s = np.linalg.norm(v)
    if s < 1e-12:
        if c > 0.0:
            return np.eye(3)
        # Antiparallel: rotate pi about any axis perpendicular to a.
        p = _unit(np.cross(a, [1.0, 0.0, 0.0] if abs(a[0]) < 0.9 else [0.0, 1.0, 0.0]))
        return 2.0 * np.outer(p, p) - np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1.0 - c) / s**2)


def hydrolysis_complex(
    dimer: Cluster,
    attacker: Cluster,
    *,
    approach_a: float = 3.2,
    mode: str = "flank",
) -> Cluster:
    """Attacker positioned to attack the dimer's Si (atom 1; atom 0 is
    the bridging O). The geometry is a starting guess; the pipeline
    optimizes it.

    ``mode="flank"`` (default): the 4-center X&L geometry — attacker O
    ~80 deg off the Si->Obr axis (same side as the bridge) with its
    first H aimed at the bridging O, so the concerted proton transfer
    is geometrically available. ``mode="backside"``: SN2-like approach
    opposite the bridge (proton transfer NOT reachable from here —
    measured live: nearest water H sat 3.9 A from the bridging O and
    driving it across never produced an interior maximum).
    """
    bridge, si = dimer.coords[0], dimer.coords[1]
    to_bridge = _unit(bridge - si)
    if mode == "backside":
        target = si + approach_a * -to_bridge
        rotated = attacker.coords - attacker.coords[0]
    elif mode == "flank":
        ref = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(ref, to_bridge)) > 0.9:
            ref = np.array([0.0, 1.0, 0.0])
        perp = _unit(np.cross(np.cross(to_bridge, ref), to_bridge))
        theta = math.radians(80.0)
        target = si + approach_a * (
            math.cos(theta) * to_bridge + math.sin(theta) * perp
        )
        shifted = attacker.coords - attacker.coords[0]
        if len(attacker.symbols) > 1:
            # Aim the first O-H bond at the bridging O.
            rot = _rotation_between(shifted[1], bridge - target)
            rotated = shifted @ rot.T
        else:
            rotated = shifted
    else:
        raise ValueError(f"unknown attack mode '{mode}'")
    moved = Cluster(
        name=attacker.name,
        symbols=list(attacker.symbols),
        coords=rotated + target,
        charge=attacker.charge,
        spin=attacker.spin,
    )
    return merge(dimer, moved, name=f"{dimer.name}+{attacker.name}")


def protonated_bridge_complex(
    dimer: Cluster,
    *,
    approach_a: float = 3.2,
    mode: str = "flank",
    n_water: int = 1,
) -> Cluster:
    """Pre-equilibrated acid complex with one or 3--6 explicit waters.

    Acid hydrolysis is not neutral hydrolysis with hydronium used as the
    nucleophile. H3O+ first transfers one proton to the bridging oxygen; the
    water left behind then attacks the weakened Si-O(bridge) bond. Preserve the
    hydronium atom order (attacker O followed by three H atoms), then append
    each additional solvent water as O/H/H so the Phase-1 driver can derive a
    stable mobile-proton/solvent index map.
    """
    _validate_acid_water_count(n_water)
    complex_ = hydrolysis_complex(
        dimer,
        hydronium(),
        approach_a=approach_a,
        mode=mode,
    )
    ow_index = len(dimer.symbols)
    bridge_proton = ow_index + 1
    coords = complex_.coords.copy()
    bridge_to_si = _unit(coords[1] - coords[0])
    toward_attacker = coords[ow_index] - coords[0]
    proton_direction = (
        toward_attacker - np.dot(toward_attacker, bridge_to_si) * bridge_to_si
    )
    if np.linalg.norm(proton_direction) < 1.0e-8:
        reference = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(reference, bridge_to_si)) > 0.9:
            reference = np.array([0.0, 1.0, 0.0])
        proton_direction = np.cross(bridge_to_si, reference)
    coords[bridge_proton] = coords[0] + R_O_H * _unit(proton_direction)
    for water_h in (ow_index + 2, ow_index + 3):
        coords[water_h] = coords[ow_index] + R_O_H * _unit(
            coords[water_h] - coords[ow_index]
        )

    symbols = list(complex_.symbols)
    if n_water > 1:
        # Put the transferred proton into a deterministic solvent-facing
        # orientation, then seed a compact network above and below the reactive
        # plane. The first added water accepts the Obr-H hydrogen bond; the
        # others solvate the residual attacker. Optimization establishes the
        # actual microsolvated minimum.
        attacker_vector = coords[ow_index] - coords[0]
        normal_raw = np.cross(bridge_to_si, attacker_vector)
        if np.linalg.norm(normal_raw) < 1.0e-8:
            reference = np.array([0.0, 0.0, 1.0])
            if abs(np.dot(reference, bridge_to_si)) > 0.9:
                reference = np.array([0.0, 1.0, 0.0])
            normal_raw = np.cross(bridge_to_si, reference)
        normal = _unit(normal_raw)
        coords[bridge_proton] = coords[0] + R_O_H * normal
        midpoint = 0.5 * (coords[0] + coords[ow_index])
        outward = _unit(attacker_vector)
        in_plane = _unit(np.cross(normal, outward))
        placements = [
            (coords[0], normal, 2.75),
            (coords[ow_index], normal, 2.75),
            (coords[ow_index], -normal, 2.75),
            (coords[0], -normal, 2.75),
            (midpoint, in_plane + 0.75 * normal, 3.10),
        ]
        expanded = list(coords)
        for origin, direction, radius in placements[: n_water - 1]:
            symbols += ["O", "H", "H"]
            expanded += [origin + point for point in _water_at(direction, radius)]
        coords = np.array(expanded)

    return Cluster(
        name=f"{dimer.name}+protonated-bridge+{n_water}water",
        symbols=symbols,
        coords=coords,
        charge=complex_.charge,
        spin=complex_.spin,
        frozen_indices=list(complex_.frozen_indices),
        site_family=complex_.site_family,
        note=(
            "H3O+ pre-equilibrated to a protonated bridging oxygen; "
            f"the residual H2O is the attacker in {n_water} explicit waters"
        ),
    )


BENCHMARKS = {
    "water": water,
    "hydronium": hydronium,
    "hydroxide": hydroxide,
    "silicic-acid": silicic_acid,
    "disilicate": disilicate,
    "aluminosilicate-dimer": aluminosilicate_dimer,
}
