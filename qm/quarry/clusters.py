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
ACID_MICROSOLVATION_FAMILIES = (
    "bridge-donor-chain",
    "attacker-centered-ring",
    "split-bridge-attacker",
    "compact-cyclic-relay",
)
CONCERTED_ACID_WATER_COUNTS = (3, 4)
CONCERTED_ACID_FAMILIES = ("bridge-donor-chain", "compact-cyclic-relay")


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
    return _water_centered(radius * d, d)


def _water_centered(oxygen: np.ndarray, direction: np.ndarray) -> list[np.ndarray]:
    """Return O/H/H coordinates with the HOH bisector along ``direction``."""
    d = _unit(direction)
    ref = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(ref, d)) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])
    e1 = _unit(np.cross(d, ref))
    half = math.radians(104.5 / 2.0)
    o = np.asarray(oxygen, dtype=float)
    h1 = o + R_O_H * (math.cos(half) * d + math.sin(half) * e1)
    h2 = o + R_O_H * (math.cos(half) * d - math.sin(half) * e1)
    return [o, h1, h2]


def _water_donating_to(
    oxygen: np.ndarray,
    acceptor: np.ndarray,
    plane_hint: np.ndarray,
) -> list[np.ndarray]:
    """Return O/H/H with one O-H bond aimed at an H-bond acceptor."""
    donor = np.asarray(oxygen, dtype=float)
    direction = _unit(np.asarray(acceptor, dtype=float) - donor)
    perpendicular = np.asarray(plane_hint, dtype=float)
    perpendicular -= np.dot(perpendicular, direction) * direction
    if np.linalg.norm(perpendicular) < 1.0e-8:
        reference = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(reference, direction)) > 0.9:
            reference = np.array([0.0, 1.0, 0.0])
        perpendicular = np.cross(direction, reference)
    perpendicular = _unit(perpendicular)
    angle = math.radians(104.5)
    h1 = donor + R_O_H * direction
    h2 = donor + R_O_H * (math.cos(angle) * direction + math.sin(angle) * perpendicular)
    return [donor, h1, h2]


def _rotate_vector(vector: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
    unit_axis = _unit(axis)
    return (
        vector * math.cos(angle)
        + np.cross(unit_axis, vector) * math.sin(angle)
        + unit_axis * np.dot(unit_axis, vector) * (1.0 - math.cos(angle))
    )


def _relay_arc_oxygens(
    start: np.ndarray,
    end: np.ndarray,
    count: int,
    *,
    segment_a: float,
    plane: np.ndarray,
) -> list[np.ndarray]:
    """Internal points on a major circular arc with fixed adjacent spacing."""
    chord = np.asarray(end, dtype=float) - np.asarray(start, dtype=float)
    chord_length = float(np.linalg.norm(chord))
    segments = count + 1
    if segment_a <= chord_length / segments:
        raise ValueError("relay segment is too short for the endpoint chord")

    def spacing(theta: float) -> float:
        return chord_length * math.sin(theta / (2.0 * segments)) / math.sin(theta / 2.0)

    low, high = 1.0e-8, 2.0 * math.pi - 1.0e-8
    for _ in range(100):
        midpoint_angle = 0.5 * (low + high)
        if spacing(midpoint_angle) < segment_a:
            low = midpoint_angle
        else:
            high = midpoint_angle
    theta = 0.5 * (low + high)
    radius = chord_length / (2.0 * math.sin(theta / 2.0))
    edge = chord / chord_length
    radial = np.asarray(plane, dtype=float)
    radial -= np.dot(radial, edge) * edge
    radial = _unit(radial)
    chord_midpoint = 0.5 * (
        np.asarray(start, dtype=float) + np.asarray(end, dtype=float)
    )
    center = chord_midpoint - radius * math.cos(theta / 2.0) * radial
    return [
        center
        + radius
        * (
            math.cos(-theta / 2.0 + index * theta / segments) * radial
            + math.sin(-theta / 2.0 + index * theta / segments) * edge
        )
        for index in range(1, segments)
    ]


def _attacker_ring_oxygens(
    bridge: np.ndarray,
    attacker: np.ndarray,
    count: int,
    *,
    normal: np.ndarray,
    in_plane: np.ndarray,
) -> list[np.ndarray]:
    """Complete a regular solvent ring whose first water also accepts Obr-H."""
    total_waters = count + 1
    side = 2.85
    chord = np.asarray(attacker, dtype=float) - np.asarray(bridge, dtype=float)
    chord_length = float(np.linalg.norm(chord))
    if chord_length >= 2.0 * side:
        raise ValueError("bridge and attacker are too far apart for the relay ring")
    ring_radial = _unit(normal + in_plane)
    first = (
        0.5 * (np.asarray(bridge, dtype=float) + np.asarray(attacker, dtype=float))
        + math.sqrt(side**2 - (0.5 * chord_length) ** 2) * ring_radial
    )
    edge = first - attacker
    edge_midpoint = 0.5 * (first + attacker)
    perpendicular = np.asarray(in_plane, dtype=float)
    perpendicular -= np.dot(perpendicular, edge) / np.dot(edge, edge) * edge
    perpendicular = _unit(perpendicular)
    apothem = side / (2.0 * math.tan(math.pi / total_waters))
    center = edge_midpoint + apothem * perpendicular
    first_vector = np.asarray(attacker, dtype=float) - center
    second_vector = first - center
    axis = _unit(np.cross(first_vector, second_vector))
    step = 2.0 * math.pi / total_waters
    vertices = [
        center + _rotate_vector(first_vector, axis, index * step)
        for index in range(total_waters)
    ]
    if not np.allclose(vertices[1], first, atol=1.0e-8):
        raise RuntimeError("regular relay ring construction lost its first edge")
    return vertices[1:]


def acid_microsolvation_hbond_edges(
    ow_index: int,
    n_water: int,
    family: str,
) -> tuple[tuple[int, int], ...]:
    """Directed donor-O -> acceptor-O graph encoded by one seed family."""
    if family not in ACID_MICROSOLVATION_FAMILIES:
        allowed = ", ".join(ACID_MICROSOLVATION_FAMILIES)
        raise ValueError(f"conformer_family must be one of {allowed}")
    extra = tuple(ow_index + 4 + 3 * index for index in range(n_water - 1))
    chain = tuple(zip(extra, extra[1:], strict=False))
    if family == "attacker-centered-ring":
        return ((0, extra[0]), (ow_index, extra[0]), *chain, (extra[-1], ow_index))
    if family == "split-bridge-attacker":
        return ((0, extra[0]), *chain, (ow_index, extra[-1]))
    if family == "compact-cyclic-relay":
        return ((0, extra[0]), *chain, (extra[-1], ow_index), (ow_index, 0))
    return ((0, extra[0]), *chain, (extra[-1], ow_index))


def _acid_conformer_waters(
    coords: np.ndarray,
    *,
    ow_index: int,
    n_water: int,
    family: str,
) -> tuple[list[np.ndarray], np.ndarray, np.ndarray | None]:
    """Return O/H/H blocks plus bridge/attacker H-bond targets."""
    if family not in ACID_MICROSOLVATION_FAMILIES:
        allowed = ", ".join(ACID_MICROSOLVATION_FAMILIES)
        raise ValueError(f"conformer_family must be one of {allowed}")
    count = n_water - 1
    bridge = coords[0]
    attacker = coords[ow_index]
    bridge_to_si = _unit(coords[1] - bridge)
    outward = _unit(attacker - bridge)
    normal_raw = np.cross(bridge_to_si, outward)
    if np.linalg.norm(normal_raw) < 1.0e-8:
        reference = np.array([0.0, 0.0, 1.0])
        if abs(np.dot(reference, bridge_to_si)) > 0.9:
            reference = np.array([0.0, 1.0, 0.0])
        normal_raw = np.cross(bridge_to_si, reference)
    normal = _unit(normal_raw)
    in_plane = _unit(np.cross(normal, outward))

    if family == "attacker-centered-ring":
        oxygens = _attacker_ring_oxygens(
            bridge,
            attacker,
            count,
            normal=normal,
            in_plane=in_plane,
        )
        attacker_target: np.ndarray | None = oxygens[0]
        bridge_target = oxygens[0] + in_plane
        targets = [*oxygens[1:], attacker]
    else:
        bridge_target = None
        plane = normal
        segment = 2.75
        if family == "split-bridge-attacker":
            plane = normal + 0.40 * in_plane
        elif family == "compact-cyclic-relay":
            plane = normal - 2.0 * in_plane
            segment = 2.60
        oxygens = _relay_arc_oxygens(
            bridge,
            attacker,
            count,
            segment_a=segment,
            plane=plane,
        )
        bridge_target = oxygens[0]
        if family == "split-bridge-attacker":
            attacker_target = oxygens[-1]
            targets = [*oxygens[1:], oxygens[-1] + normal]
        else:
            attacker_target = bridge if family == "compact-cyclic-relay" else None
            targets = [*oxygens[1:], attacker]

    waters = [
        point
        for oxygen, target in zip(oxygens, targets, strict=True)
        for point in _water_donating_to(oxygen, target, normal)
    ]
    return waters, bridge_target, attacker_target


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
    conformer_family: str | None = None,
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
    if n_water == 1 and conformer_family is not None:
        raise ValueError("conformer_family is not valid for n_water=1")
    if (
        conformer_family is not None
        and conformer_family not in ACID_MICROSOLVATION_FAMILIES
    ):
        allowed = ", ".join(ACID_MICROSOLVATION_FAMILIES)
        raise ValueError(f"conformer_family must be one of {allowed}")
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
        expanded = list(coords)
        if conformer_family is None:
            placements = [
                (coords[0], normal, 2.75),
                (coords[ow_index], normal, 2.75),
                (coords[ow_index], -normal, 2.75),
                (coords[0], -normal, 2.75),
                (midpoint, in_plane + 0.75 * normal, 3.10),
            ]
            for origin, direction, radius in placements[: n_water - 1]:
                symbols += ["O", "H", "H"]
                expanded += [origin + point for point in _water_at(direction, radius)]
        else:
            symbols += ["O", "H", "H"] * (n_water - 1)
            waters, bridge_target, attacker_target = _acid_conformer_waters(
                coords,
                ow_index=ow_index,
                n_water=n_water,
                family=conformer_family,
            )
            coords[bridge_proton] = coords[0] + R_O_H * _unit(bridge_target - coords[0])
            if attacker_target is not None:
                attacker_water = _water_donating_to(
                    coords[ow_index], attacker_target, normal
                )
                coords[ow_index + 2] = attacker_water[1]
                coords[ow_index + 3] = attacker_water[2]
            expanded = [*coords, *waters]
        coords = np.array(expanded)

    family_suffix = f"+{conformer_family}" if conformer_family else ""
    family_note = (
        f"; microsolvation family {conformer_family}" if conformer_family else ""
    )
    return Cluster(
        name=f"{dimer.name}+protonated-bridge+{n_water}water{family_suffix}",
        symbols=symbols,
        coords=coords,
        charge=complex_.charge,
        spin=complex_.spin,
        frozen_indices=list(complex_.frozen_indices),
        site_family=complex_.site_family,
        note=(
            "H3O+ pre-equilibrated to a protonated bridging oxygen; "
            f"the residual H2O is the attacker in {n_water} explicit waters"
            f"{family_note}"
        ),
    )


@dataclass(frozen=True)
class ConcertedRelayEndpoints:
    """Atom-matched endpoints for one concerted hydronium relay path."""

    reactant: Cluster
    product: Cluster
    ow_index: int
    transferred_h_index: int
    relay_h_indices: tuple[int, ...]
    solvent_oxygen_indices: tuple[int, ...]


def concerted_acid_relay_endpoints(
    dimer: Cluster,
    *,
    n_water: int,
    family: str,
) -> ConcertedRelayEndpoints:
    """Build unconstrained R/P seeds for concerted acid hydrolysis.

    Atom order matches the A1c contract: dimer atoms, hydronium ``O/H/H/H``,
    then ``O/H/H`` shell waters. The reactant has an intact, unprotonated
    bridge and H3O+. The product transfers exactly H1 to Obr, forms Si--Ow,
    and cleaves Si--Obr. These are seeds; no production coordinate is frozen.
    """
    if n_water not in CONCERTED_ACID_WATER_COUNTS:
        raise ValueError("concerted acid relay requires exactly 3 or 4 waters")
    if family not in CONCERTED_ACID_FAMILIES:
        allowed = ", ".join(CONCERTED_ACID_FAMILIES)
        raise ValueError(f"concerted acid relay family must be one of {allowed}")

    # Reuse A1c's non-colliding oxygen topology and atom order, then reverse its
    # proton wire: hydronium -> shell -> Obr. Returning H1 from Obr to Ow is the
    # load-bearing distinction from the retired pre-equilibrium mechanism.
    template = protonated_bridge_complex(
        dimer,
        n_water=n_water,
        conformer_family=family,
    )
    ow_index = len(dimer.symbols)
    transferred_h = ow_index + 1
    shell_oxygens = tuple(ow_index + 4 + 3 * index for index in range(n_water - 1))
    solvent_oxygens = (ow_index, *shell_oxygens)
    coords = template.coords.copy()

    attacker = coords[ow_index]
    donor_target = 0 if family == "compact-cyclic-relay" else shell_oxygens[-1]
    hydronium_seed = hydronium()
    relative = hydronium_seed.coords - hydronium_seed.coords[0]
    rotation = _rotation_between(relative[1], coords[donor_target] - attacker)
    oriented_hydronium = relative @ rotation.T + attacker
    coords[ow_index + 1 : ow_index + 4] = oriented_hydronium[1:]
    if family == "bridge-donor-chain":
        for position, oxygen_index in enumerate(shell_oxygens):
            target_index = 0 if position == 0 else shell_oxygens[position - 1]
            water_block = _water_donating_to(
                coords[oxygen_index], coords[target_index], np.array([0.0, 0.0, 1.0])
            )
            coords[oxygen_index + 1] = water_block[1]
            coords[oxygen_index + 2] = water_block[2]

    reactant = Cluster(
        name=f"{dimer.name}+concerted-acid-{n_water}w-{family}-reactant",
        symbols=list(template.symbols),
        coords=coords,
        charge=template.charge,
        spin=template.spin,
        frozen_indices=list(template.frozen_indices),
        site_family=template.site_family,
        note="intact unprotonated bridge; H3O+ -> solvent relay -> Obr",
    )

    product_coords = (
        coords.copy() if family == "bridge-donor-chain" else template.coords.copy()
    )
    # Dimer atom order is fixed: 1,3..8 are the attacked Si(OH)3 fragment.
    old_si = product_coords[1].copy()
    old_attacker = product_coords[ow_index].copy()
    separation = 1.15 * _unit(old_si - product_coords[0])
    for index in (1, 3, 4, 5, 6, 7, 8):
        product_coords[index] += separation
    attack_direction = _unit(old_attacker - old_si)
    new_attacker = product_coords[1] + 1.70 * attack_direction
    solvent_shift = new_attacker - old_attacker
    product_coords[ow_index:] += solvent_shift
    if family == "bridge-donor-chain":
        # Grotthuss-style physical-H relay: H1 moves from H3O+ to the last
        # shell water, and each shell donor proton moves one owner toward Obr.
        bridge_proton = shell_oxygens[0] + 1
        product_coords[bridge_proton] = product_coords[0] + R_O_H * _unit(
            product_coords[shell_oxygens[0]] - product_coords[0]
        )
        for position, oxygen_index in enumerate(shell_oxygens):
            previous_owner = 0 if position == 0 else shell_oxygens[position - 1]
            incoming_proton = (
                shell_oxygens[position + 1] + 1
                if position + 1 < len(shell_oxygens)
                else transferred_h
            )
            water_block = _water_centered(
                product_coords[oxygen_index],
                product_coords[oxygen_index] - product_coords[previous_owner],
            )
            product_coords[incoming_proton] = water_block[1]
            product_coords[oxygen_index + 2] = water_block[2]
        relay_h_indices = (transferred_h, *(oxygen + 1 for oxygen in shell_oxygens))
    else:
        product_coords[transferred_h] = product_coords[0] + R_O_H * _unit(
            product_coords[shell_oxygens[0]] - product_coords[0]
        )
        bridge_proton = transferred_h
        relay_h_indices = (transferred_h,)

    # Reorient the product Ow hydrogens away from both the attacked framework
    # and the translated shell. Pick the deterministic candidate with the
    # largest nearest-neighbor clearance; this is a geometry seed, not a gate.
    normal = np.cross(product_coords[1] - product_coords[0], attack_direction)
    if np.linalg.norm(normal) < 1.0e-8:
        normal = np.array([0.0, 0.0, 1.0])
    shell_center = np.mean(product_coords[list(shell_oxygens)], axis=0)
    candidate_directions = (
        new_attacker - product_coords[7],
        new_attacker - shell_center,
        new_attacker - product_coords[1],
        normal,
        -normal,
    )
    h2, h3 = ow_index + 2, ow_index + 3
    fixed = np.delete(product_coords, [ow_index, h2, h3], axis=0)
    other_oxygen_indices = [
        index
        for index, symbol in enumerate(template.symbols)
        if symbol == "O" and index != ow_index
    ]
    best_hydrogens = None
    best_score = (-1, -math.inf)
    for direction in candidate_directions:
        water_block = _water_centered(new_attacker, direction)
        hydrogens = np.asarray(water_block[1:])
        clearance = min(
            float(np.linalg.norm(hydrogens[0] - hydrogens[1])),
            float(
                np.min(
                    np.linalg.norm(hydrogens[:, None, :] - fixed[None, :, :], axis=2)
                )
            ),
        )
        oxygen_clearance = float(
            np.min(
                np.linalg.norm(
                    hydrogens[:, None, :]
                    - product_coords[other_oxygen_indices][None, :, :],
                    axis=2,
                )
            )
        )
        score = (int(oxygen_clearance >= 1.05), clearance)
        if score > best_score:
            best_score = score
            best_hydrogens = hydrogens
    assert best_hydrogens is not None
    product_coords[[h2, h3]] = best_hydrogens

    product = Cluster(
        name=f"{dimer.name}+concerted-acid-{n_water}w-{family}-product",
        symbols=list(template.symbols),
        coords=product_coords,
        charge=template.charge,
        spin=template.spin,
        frozen_indices=list(template.frozen_indices),
        site_family=template.site_family,
        note="hydrolyzed product; hydronium proton delivered to Obr during attack",
    )
    return ConcertedRelayEndpoints(
        reactant=reactant,
        product=product,
        ow_index=ow_index,
        transferred_h_index=bridge_proton,
        relay_h_indices=relay_h_indices,
        solvent_oxygen_indices=solvent_oxygens,
    )


BENCHMARKS = {
    "water": water,
    "hydronium": hydronium,
    "hydroxide": hydroxide,
    "silicic-acid": silicic_acid,
    "disilicate": disilicate,
    "aluminosilicate-dimer": aluminosilicate_dimer,
}
