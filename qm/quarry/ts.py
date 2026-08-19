"""Transition-state machinery: Sella saddle search over a PySCF calculator.

The stack per SURVEY.md §4.1: Sella does partitioned-RFO saddle-point
optimization over any ASE calculator; PySCF supplies energies/forces
(CPU or GPU via ``DftSettings.use_gpu``). Verification is quarry's own:
a TS must have exactly one imaginary mode, and the quick-IRC check
(displace ± along that mode, relax to minima) must connect reactant and
product basins. Full Gonzalez-Schlegel IRC can land later if a reaction
needs it; for hydrolysis barriers the displace-and-relax check is the
standard first line.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from quarry.clusters import Cluster
from quarry.pipeline import (
    DftSettings,
    FrequencyResult,
    _make_scf,
    build_mol,
    frequencies,
    optimize,
)

HARTREE_TO_EV = 27.211386245988
BOHR_TO_ANGSTROM = 0.529177210903


def make_ase_calculator(settings: DftSettings, charge: int, spin: int):
    """An ASE calculator that evaluates PySCF energy/forces for quarry.

    Constructed per-cluster (charge/spin are electronic state, not
    geometry), so Sella and other ASE drivers can move atoms freely.
    """
    from ase.calculators.calculator import Calculator, all_changes

    class PyscfCalculator(Calculator):
        implemented_properties = ["energy", "forces"]

        def calculate(
            self, atoms=None, properties=("energy",), system_changes=all_changes
        ):
            super().calculate(atoms, properties, system_changes)
            cluster = Cluster(
                name="ase-eval",
                symbols=list(atoms.symbols),
                coords=atoms.positions.copy(),
                charge=charge,
                spin=spin,
            )
            mf = _make_scf(build_mol(cluster, settings), settings)
            e = mf.kernel()
            if not mf.converged:
                raise RuntimeError("SCF did not converge during ASE evaluation")
            grad = mf.nuc_grad_method().kernel()
            grad = np.asarray(grad)
            if hasattr(grad, "get"):  # cupy -> numpy when running on GPU
                grad = grad.get()
            self.results = {
                "energy": float(e) * HARTREE_TO_EV,
                "forces": -grad * (HARTREE_TO_EV / BOHR_TO_ANGSTROM),
            }

    return PyscfCalculator()


def find_ts(
    cluster: Cluster,
    settings: DftSettings,
    *,
    fmax_ev_a: float = 0.02,
    max_steps: int = 300,
    trajectory: str | None = None,
) -> Cluster:
    """First-order saddle search (Sella, partitioned RFO) from a TS guess.

    ``fmax_ev_a`` is the ASE force-convergence threshold in eV/Angstrom
    (0.02 ~ 4e-4 Hartree/Bohr). Frozen atoms in the cluster are held
    fixed (the lattice-resistance contract). Raises if Sella does not
    converge within ``max_steps``.
    """
    from ase import Atoms
    from ase.constraints import FixAtoms
    from sella import Sella

    atoms = Atoms(symbols=cluster.symbols, positions=cluster.coords)
    atoms.calc = make_ase_calculator(settings, cluster.charge, cluster.spin)
    if cluster.frozen_indices:
        atoms.set_constraint(FixAtoms(indices=cluster.frozen_indices))
    # internal=True uses Sella's internal coordinates — the right choice
    # for molecular clusters; order=1 requests a first-order saddle.
    dyn = Sella(atoms, order=1, internal=True, trajectory=trajectory)
    converged = dyn.run(fmax=fmax_ev_a, steps=max_steps)
    if converged is None:
        # Older ASE returned None from run(); fall back to the force test
        # (forces are cached from the last step, no recompute).
        converged = float(np.abs(atoms.get_forces()).max()) < fmax_ev_a
    if not converged:
        raise RuntimeError(
            f"Sella did not converge a saddle for {cluster.name} "
            f"within {max_steps} steps"
        )
    return replace(cluster, coords=atoms.positions.copy(), name=f"{cluster.name}-ts")


def verify_ts(cluster: Cluster, settings: DftSettings) -> FrequencyResult:
    """Frequency-verify a saddle: exactly one imaginary mode or raise."""
    freq = frequencies(cluster, settings)
    if freq.n_imaginary != 1:
        raise RuntimeError(
            f"{cluster.name}: expected exactly 1 imaginary mode, "
            f"found {freq.n_imaginary} "
            f"(imaginary: {np.round(freq.imaginary_cm, 1).tolist()} cm^-1)"
        )
    return freq


def quick_irc(
    ts: Cluster,
    settings: DftSettings,
    *,
    displacement_a: float = 0.15,
    max_steps: int = 200,
) -> tuple[Cluster, Cluster]:
    """Displace ± along the imaginary mode and relax to the two minima.

    The cheap stand-in for a full IRC: confirms which basins the saddle
    connects. Returns (backward, forward) relaxed clusters.
    """
    freq = verify_ts(ts, settings)
    mode = freq.imaginary_mode
    if mode is None:
        raise RuntimeError(f"{ts.name}: no imaginary-mode vector available")
    step = displacement_a * mode / np.linalg.norm(mode)
    ends = []
    for sign, tag in ((-1.0, "back"), (+1.0, "fwd")):
        displaced = replace(ts, coords=ts.coords + sign * step, name=f"{ts.name}-{tag}")
        ends.append(optimize(displaced, settings, max_steps=max_steps))
    return ends[0], ends[1]


def constrained_scan(
    cluster: Cluster,
    settings: DftSettings,
    *,
    atom_i: int,
    atom_j: int,
    distances_a: list[float],
    fixed_distances: list[tuple[int, int, float]] = (),
    max_steps: int = 150,
) -> list[tuple[float, float, Cluster]]:
    """Relaxed scan of one interatomic distance (TS-guess generator).

    Optimizes the cluster with r(i,j) frozen at each value; returns
    (distance, energy_hartree, relaxed_cluster) tuples in scan order.
    ``fixed_distances`` holds additional (i, j, r) pairs constrained at
    every point — e.g. pinning the nucleophile in place while a proton
    is driven. The maximum-energy point is the natural Sella guess.
    """
    from pyscf.geomopt.geometric_solver import optimize as geometric_optimize

    from quarry.pipeline import constraints_file

    results: list[tuple[float, float, Cluster]] = []
    current = cluster
    for r in distances_a:
        constraint = f"$set\ndistance {atom_i + 1} {atom_j + 1} {r:.4f}\n"
        for fi, fj, fr in fixed_distances:
            constraint += f"distance {fi + 1} {fj + 1} {fr:.4f}\n"
        if current.frozen_indices:
            frozen = ",".join(str(k + 1) for k in sorted(current.frozen_indices))
            constraint += f"$freeze\nxyz {frozen}\n"
        mf = _make_scf(build_mol(current, settings), settings)
        with constraints_file(constraint) as path:
            mol_opt = geometric_optimize(mf, maxsteps=max_steps, constraints=path)
        coords = mol_opt.atom_coords() * BOHR_TO_ANGSTROM
        current = replace(current, coords=coords, name=f"{cluster.name}-r{r:.2f}")
        mf2 = _make_scf(build_mol(current, settings), settings)
        e = float(mf2.kernel())
        results.append((r, e, current))
    return results


def scan_maximum(scan: list[tuple[float, float, Cluster]]) -> Cluster:
    """The highest-energy scan point — the TS guess."""
    if not scan:
        raise ValueError("empty scan")
    return max(scan, key=lambda p: p[1])[2]


# Peaks must clear both neighbors by this much (Hartree) to count —
# well above converged-optimization noise, well below any real barrier.
_PEAK_TOL_HARTREE = 5e-5


def first_interior_maximum(
    scan: list[tuple[float, float, Cluster]], *, tol: float = _PEAK_TOL_HARTREE
) -> int | None:
    """Index of the first interior local energy maximum, else None.

    The *first* crest along the driven coordinate is the saddle guess
    for this elementary step. The global maximum is the wrong criterion:
    once the driven bond has actually formed, pushing further just
    compresses it and the energy rises without bound, so the global max
    lands on the compression endpoint (observed live: proton-transfer
    crest at r(Obr-H)=1.96 A, then bond formation, then a monotonic
    compression wall that out-climbed the crest).
    """
    energies = [e for _, e, _ in scan]
    for i in range(1, len(energies) - 1):
        if energies[i] > energies[i - 1] + tol and energies[i] > energies[i + 1] + tol:
            return i
    return None


def has_interior_maximum(scan: list[tuple[float, float, Cluster]]) -> bool:
    """True when the scan contains a genuine interior crest.

    A maximum at the tight end means the ridge has not been crossed and
    the "TS guess" is just the last constrained point — a saddle search
    from it slides back into the reactant basin (observed live: Sella
    relaxed a monotonic scan's endpoint into a 208i cm^-1 water-wag
    saddle with a 2.5 kJ/mol "barrier"). Extend the scan instead.
    """
    if len(scan) < 3:
        return False
    return first_interior_maximum(scan) is not None


def scan_ts_guess(scan: list[tuple[float, float, Cluster]]) -> Cluster:
    """The structure at the first interior crest — the Sella guess."""
    idx = first_interior_maximum(scan)
    if idx is None:
        raise ValueError("scan has no interior energy maximum")
    return scan[idx][2]


def neb_ts_guess(
    reactant: Cluster,
    product: Cluster,
    settings: DftSettings,
    *,
    n_images: int = 7,
    fmax_ev_a: float = 0.10,
    max_steps: int = 150,
) -> Cluster:
    """Climbing-image NEB between two basins; returns the peak image.

    The remedy for saddle searches that escape the reaction channel
    (observed live: Sella walked a crest guess from r(Si-Ow)=1.90 back
    out to 3.22 A): the band is pinned to both basins so it cannot
    escape, and the climbing image rides the path to the col, where the
    lowest curvature direction *is* the reaction mode. The peak image
    is then a Sella start it can refine rather than flee.

    ``fmax_ev_a`` is deliberately loose — this produces a guess, and
    Sella does the tight convergence.
    """
    from ase import Atoms

    try:
        from ase.mep import NEB
    except ImportError:  # older ASE layout
        from ase.neb import NEB
    from ase.optimize import FIRE

    if reactant.charge != product.charge or reactant.spin != product.spin:
        raise ValueError("reactant/product electronic states differ")
    images = [Atoms(symbols=reactant.symbols, positions=reactant.coords)]
    for _ in range(n_images - 2):
        images.append(images[0].copy())
    images.append(Atoms(symbols=product.symbols, positions=product.coords))
    for image in images:
        image.calc = make_ase_calculator(settings, reactant.charge, reactant.spin)
        if reactant.frozen_indices:
            from ase.constraints import FixAtoms

            image.set_constraint(FixAtoms(indices=reactant.frozen_indices))
    neb = NEB(images, climb=True)
    neb.interpolate(method="idpp")
    FIRE(neb).run(fmax=fmax_ev_a, steps=max_steps)
    energies = [float(image.get_potential_energy()) for image in images]
    peak = int(np.argmax(energies))
    if peak in (0, len(images) - 1):
        raise RuntimeError(
            "NEB peak landed on an endpoint — the two basins do not "
            "bracket a barrier at this method/geometry"
        )
    return replace(
        reactant,
        coords=images[peak].positions.copy(),
        name=f"{reactant.name}-neb-peak",
    )


class ScanNoMaximumError(RuntimeError):
    """A scan hit its distance floor with the energy still rising.

    Carries the completed scan so callers can pivot to a second driven
    coordinate from the best structure already computed.
    """

    def __init__(self, message: str, scan: list[tuple[float, float, Cluster]]):
        super().__init__(message)
        self.scan = scan


def scan_to_maximum(
    cluster: Cluster,
    settings: DftSettings,
    *,
    atom_i: int,
    atom_j: int,
    distances_a: list[float],
    fixed_distances: list[tuple[int, int, float]] = (),
    extend_step_a: float = 0.08,
    min_distance_a: float = 1.60,
    max_steps: int = 150,
    progress=None,
) -> list[tuple[float, float, Cluster]]:
    """Relaxed scan that keeps tightening r(i,j) until the maximum is
    interior (or ``min_distance_a`` is hit — then ``ScanNoMaximumError``
    carrying the scan, since a guess from an uncrossed ridge is worse
    than no guess).

    ``progress`` is an optional callable taking (r, energy) per point.
    """
    scan = constrained_scan(
        cluster,
        settings,
        atom_i=atom_i,
        atom_j=atom_j,
        distances_a=distances_a,
        fixed_distances=fixed_distances,
        max_steps=max_steps,
    )
    if progress:
        for r, e, _ in scan:
            progress(r, e)
    while not has_interior_maximum(scan):
        next_r = round(scan[-1][0] - extend_step_a, 3)
        if next_r < min_distance_a:
            raise ScanNoMaximumError(
                f"scan reached r={scan[-1][0]:.2f} A without an interior "
                "energy maximum — the driven coordinate alone does not "
                "cross the ridge; a different or combined reaction "
                "coordinate is needed",
                scan,
            )
        ext = constrained_scan(
            scan[-1][2],
            settings,
            atom_i=atom_i,
            atom_j=atom_j,
            distances_a=[next_r],
            fixed_distances=fixed_distances,
            max_steps=max_steps,
        )
        if progress:
            for r, e, _ in ext:
                progress(r, e)
        scan.extend(ext)
    return scan
