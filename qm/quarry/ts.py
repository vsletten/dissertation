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
    max_steps: int = 150,
) -> list[tuple[float, float, Cluster]]:
    """Relaxed scan of one interatomic distance (TS-guess generator).

    Optimizes the cluster with r(i,j) frozen at each value; returns
    (distance, energy_hartree, relaxed_cluster) tuples in scan order.
    The maximum-energy point is the natural Sella starting guess.
    """
    from pyscf.geomopt.geometric_solver import optimize as geometric_optimize

    results: list[tuple[float, float, Cluster]] = []
    current = cluster
    for r in distances_a:
        constraint = f"$set\ndistance {atom_i + 1} {atom_j + 1} {r:.4f}\n"
        if current.frozen_indices:
            frozen = ",".join(str(k + 1) for k in sorted(current.frozen_indices))
            constraint += f"$freeze\nxyz {frozen}\n"
        mf = _make_scf(build_mol(current, settings), settings)
        mol_opt = geometric_optimize(mf, maxsteps=max_steps, constraints=constraint)
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
