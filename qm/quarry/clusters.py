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


def silicic_acid_hydrate(n_water: int) -> Cluster:
    """Si(OH)4 · n H2O — the Phase 0 smoke-test series (SURVEY.md §9).

    Waters are placed on a loose shell around the cluster; only the
    topology matters (the smoke test measures E+grad+Hessian cost, and
    real work pre-optimizes).
    """
    if not 0 <= n_water <= 6:
        raise ValueError("n_water must be in 0..6")
    base = silicic_acid()
    symbols = list(base.symbols)
    coords = list(base.coords)
    w = water()
    shell_dirs = [
        (1, 0, 0),
        (-1, 0, 0),
        (0, 1, 0),
        (0, -1, 0),
        (0, 0, 1),
        (0, 0, -1),
    ]
    for k in range(n_water):
        offset = 3.4 * _unit(np.array(shell_dirs[k], dtype=float))
        symbols += w.symbols
        coords += list(w.coords + offset)
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


BENCHMARKS = {
    "water": water,
    "hydronium": hydronium,
    "hydroxide": hydroxide,
    "silicic-acid": silicic_acid,
    "disilicate": disilicate,
    "aluminosilicate-dimer": aluminosilicate_dimer,
}
