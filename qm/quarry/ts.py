"""Transition-state machinery: Sella saddle search over a PySCF calculator.

The stack per SURVEY.md §4.1: Sella does partitioned-RFO saddle-point
optimization over any ASE calculator; PySCF supplies energies/forces
(CPU or GPU via ``DftSettings.use_gpu``). Verification is quarry's own:
a TS must have exactly one imaginary mode.
Most campaigns use the cheap quick-IRC check (displace ± along that mode,
then relax to minima). Mechanisms whose contract requires the actual
mass-weighted path use :func:`full_irc`, backed by Sella's
Gonzalez--Schlegel integrator.
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict, replace
from pathlib import Path
from types import MethodType
from typing import Any

import numpy as np

from quarry.clusters import Cluster
from quarry.pipeline import (
    HARTREE_TO_KJ,
    DftSettings,
    FrequencyResult,
    _gradient_method,
    _make_scf,
    build_mol,
    frequencies,
    frequency_geometry_fingerprint,
    frequency_settings_fingerprint,
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
            grad = _gradient_method(mf, settings).kernel()
            if hasattr(grad, "get"):  # cupy -> numpy when running on GPU
                grad = grad.get()
            grad = np.asarray(grad)
            self.results = {
                "energy": float(e) * HARTREE_TO_EV,
                "forces": -grad * (HARTREE_TO_EV / BOHR_TO_ANGSTROM),
            }

    return PyscfCalculator()


def reaction_path_vector(
    reactant: Cluster,
    product: Cluster,
    *,
    active_indices: list[int] | None = None,
) -> np.ndarray:
    """Normalized Cartesian endpoint displacement with rigid motion removed.

    The vector is a physically motivated initial mode for a directed Cartesian
    Sella search.  Endpoint optimizations may independently rotate or translate
    a cluster, so that meaningless motion is removed before differencing.  Any
    frozen shell is excluded from the mode. ``active_indices`` can additionally
    restrict the displacement to a reactive core. This prevents unrelated
    endpoint conformer changes (for example, terminal-hydroxyl rotations) from
    choosing the initial unstable direction for a directed saddle search.
    """
    from ase import Atoms
    from ase.build.rotate import minimize_rotation_and_translation

    if reactant.symbols != product.symbols:
        raise ValueError("reactant/product atom order differs")
    if reactant.charge != product.charge or reactant.spin != product.spin:
        raise ValueError("reactant/product electronic states differ")
    if reactant.frozen_indices != product.frozen_indices:
        raise ValueError("reactant/product frozen atom sets differ")

    start = Atoms(symbols=reactant.symbols, positions=reactant.coords)
    end = Atoms(symbols=product.symbols, positions=product.coords)
    minimize_rotation_and_translation(start, end)
    mode = end.positions - start.positions
    if active_indices is not None:
        if not active_indices:
            raise ValueError("active_indices must not be empty")
        if len(set(active_indices)) != len(active_indices):
            raise ValueError("active_indices contains duplicates")
        if any(i < 0 or i >= len(reactant.symbols) for i in active_indices):
            raise ValueError("active_indices contains an out-of-range atom")
        masked = np.zeros_like(mode)
        masked[np.asarray(active_indices, dtype=int)] = mode[
            np.asarray(active_indices, dtype=int)
        ]
        mode = masked
    if reactant.frozen_indices:
        mode[np.asarray(reactant.frozen_indices, dtype=int)] = 0.0
    norm = float(np.linalg.norm(mode))
    if not np.isfinite(norm) or norm < 1e-12:
        raise ValueError("reactant/product displacement has no non-rigid component")
    return mode / norm


def saddle_search_frozen_indices(
    cluster: Cluster, active_indices: list[int] | None
) -> list[int]:
    """Return the temporary frozen set for an active-subspace saddle seed."""
    original = set(cluster.frozen_indices or [])
    if active_indices is None:
        return sorted(original)
    if not active_indices or len(set(active_indices)) != len(active_indices):
        raise ValueError("active_indices must be non-empty and unique")
    atom_count = len(cluster.symbols)
    active = set(active_indices)
    if any(index < 0 or index >= atom_count for index in active):
        raise ValueError("active_indices contains an out-of-range atom")
    overlap = active & original
    if overlap:
        raise ValueError("active_indices contains an originally frozen atom")
    return sorted(original | (set(range(atom_count)) - active))


def reaction_aligned_imaginary_mode(
    frequency: FrequencyResult,
    reaction_vector: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Imaginary mode with the largest absolute overlap to a reaction vector."""
    modes = frequency.imaginary_modes
    if modes is None:
        if frequency.n_imaginary == 1 and frequency.imaginary_mode is not None:
            modes = frequency.imaginary_mode[None, ...]
        else:
            raise ValueError("all imaginary modes are required for reaction alignment")
    vector = np.asarray(reaction_vector, dtype=float)
    if modes.ndim != 3 or vector.shape != modes.shape[1:]:
        raise ValueError("reaction vector and imaginary modes have incompatible shapes")
    vector_norm = float(np.linalg.norm(vector))
    if not np.isfinite(vector_norm) or vector_norm < 1e-12:
        raise ValueError("reaction vector is non-finite or zero")
    unit_vector = vector.reshape(-1) / vector_norm
    flattened = modes.reshape(modes.shape[0], -1)
    mode_norms = np.linalg.norm(flattened, axis=1)
    if np.any(~np.isfinite(mode_norms)) or np.any(mode_norms < 1e-12):
        raise ValueError("imaginary modes contain a non-finite or zero vector")
    overlaps = np.abs((flattened / mode_norms[:, None]) @ unit_vector)
    selected = int(np.argmax(overlaps))
    return modes[selected], float(overlaps[selected])


def find_ts(
    cluster: Cluster,
    settings: DftSettings,
    *,
    fmax_ev_a: float = 0.02,
    max_steps: int = 300,
    trajectory: str | None = None,
    initial_mode: np.ndarray | None = None,
    internal: bool = True,
    active_indices: list[int] | None = None,
) -> Cluster:
    """First-order saddle search (Sella, partitioned RFO) from a TS guess.

    ``fmax_ev_a`` is the ASE force-convergence threshold in eV/Angstrom
    (0.02 ~ 4e-4 Hartree/Bohr). Frozen atoms in the cluster are held
    fixed (the lattice-resistance contract). ``initial_mode`` supplies a
    normalized Cartesian reaction direction, normally from
    :func:`reaction_path_vector`; directed searches must use
    ``internal=False`` so the mode and optimizer basis agree.
    ``active_indices`` temporarily freezes every other atom during the search;
    the returned :class:`Cluster` retains only its original frozen shell. This
    is intended for a bounded active-subspace seed, not final TS acceptance:
    callers must subsequently relax the spectator space, release the temporary
    freezes, and verify a full-system stationary point and Hessian. Raises if
    Sella does not converge within ``max_steps``.
    """
    from ase import Atoms
    from ase.constraints import FixAtoms
    from sella import Sella

    mode = None
    if initial_mode is not None:
        mode = np.asarray(initial_mode, dtype=float)
        if mode.shape != cluster.coords.shape:
            expected = cluster.coords.shape
            detail = f"got {mode.shape}, expected {expected}"
            raise ValueError(f"initial_mode shape mismatch: {detail}")
        if internal:
            raise ValueError("a Cartesian initial_mode requires internal=False")
        if not np.all(np.isfinite(mode)):
            raise ValueError("initial_mode contains non-finite values")
        mode_norm = float(np.linalg.norm(mode))
        if mode_norm < 1e-12:
            raise ValueError("initial_mode has zero norm")
        mode = mode / mode_norm

    search_frozen = saddle_search_frozen_indices(cluster, active_indices)
    atoms = Atoms(symbols=cluster.symbols, positions=cluster.coords)
    atoms.calc = make_ase_calculator(settings, cluster.charge, cluster.spin)
    if search_frozen:
        atoms.set_constraint(FixAtoms(indices=search_frozen))
    # Internal coordinates are the default for an undirected molecular search.
    # A directed endpoint mode is Cartesian and therefore opts into the
    # Cartesian PES explicitly.
    dyn = Sella(atoms, order=1, internal=internal, trajectory=trajectory)
    if mode is not None:
        free_basis = dyn.pes.get_Ufree()
        if free_basis is None:
            raise RuntimeError("Sella did not expose a free-coordinate basis")
        projected = mode.reshape(-1) @ np.asarray(free_basis)
        projected_norm = float(np.linalg.norm(projected))
        if projected_norm < 1e-12:
            raise ValueError("initial_mode has no component in the free coordinates")
        dyn.pes.v0 = projected / projected_norm
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


def relax_at_fixed_distances(
    cluster: Cluster,
    settings: DftSettings,
    *,
    fixed_distances: list[tuple[int, int, float]],
    fmax_ev_a: float = 0.02,
    max_steps: int = 120,
    optimizer_maxstep: float = 0.05,
    distance_tolerance_a: float = 1e-6,
    constraint_method: str = "ase-projector",
    trajectory: str | None = None,
    logfile: str | None = "-",
) -> Cluster:
    """Relax transverse modes while holding exact Cartesian bond distances.

    This is a topology-agnostic alternative to a geomeTRIC ``$set``
    optimization when a transition-state guess sits on a known multidimensional
    reaction-coordinate crest. ``ase-projector`` uses ASE's
    ``FixBondLengths`` for independent pairs. ``sella-internal`` solves coupled
    bond equations simultaneously inside a Cartesian order-zero Sella search;
    use it when multiple pinned distances share atoms.

    The optimizer must converge and the final projected force and distance
    residuals are checked independently. Frozen atoms may be retained only when
    they do not participate in a pinned pair; overlapping exact constraints have
    order-dependent behavior in ASE and are rejected rather than guessed at.
    """
    from ase import Atoms
    from ase.constraints import FixAtoms, FixBondLengths
    from ase.optimize import BFGS
    from sella import Sella
    from sella.internal import Constraints

    if not fixed_distances:
        raise ValueError("fixed_distances must not be empty")
    if not np.isfinite(fmax_ev_a) or fmax_ev_a <= 0.0:
        raise ValueError("fmax_ev_a must be finite and positive")
    if max_steps <= 0:
        raise ValueError("max_steps must be positive")
    if not np.isfinite(optimizer_maxstep) or optimizer_maxstep <= 0.0:
        raise ValueError("optimizer_maxstep must be finite and positive")
    if not np.isfinite(distance_tolerance_a) or distance_tolerance_a <= 0.0:
        raise ValueError("distance_tolerance_a must be finite and positive")
    if constraint_method not in {"ase-projector", "sella-internal"}:
        raise ValueError(
            "constraint_method must be 'ase-projector' or 'sella-internal'"
        )

    pairs: list[tuple[int, int]] = []
    targets: list[float] = []
    seen: set[tuple[int, int]] = set()
    frozen = set(cluster.frozen_indices or [])
    for atom_i, atom_j, target in fixed_distances:
        if atom_i == atom_j:
            raise ValueError("a fixed distance requires two distinct atoms")
        if (
            atom_i < 0
            or atom_j < 0
            or atom_i >= len(cluster.symbols)
            or atom_j >= len(cluster.symbols)
        ):
            raise ValueError("fixed_distances contains an out-of-range atom")
        pair_key = (atom_i, atom_j) if atom_i < atom_j else (atom_j, atom_i)
        if pair_key in seen:
            raise ValueError("fixed_distances contains a duplicate atom pair")
        if atom_i in frozen or atom_j in frozen:
            raise ValueError("a fixed-distance atom may not also be frozen")
        if not np.isfinite(target) or target <= 0.0:
            raise ValueError("fixed-distance targets must be finite and positive")
        seen.add(pair_key)
        pairs.append((atom_i, atom_j))
        targets.append(float(target))

    atoms = Atoms(symbols=cluster.symbols, positions=cluster.coords)
    atoms.calc = make_ase_calculator(settings, cluster.charge, cluster.spin)
    targets_array = np.asarray(targets)
    if constraint_method == "ase-projector":
        constraints = []
        if cluster.frozen_indices:
            constraints.append(FixAtoms(indices=cluster.frozen_indices))
        constraints.append(
            FixBondLengths(
                pairs,
                bondlengths=targets,
                tolerance=distance_tolerance_a,
            )
        )
        atoms.set_constraint(constraints)
        # Explicit targets are not applied until ASE adjusts a position update.
        # Project before BFGS's initial force-only convergence check, otherwise
        # an off-target zero-force geometry can falsely return converged.
        atoms.set_positions(atoms.positions.copy())
        optimizer = BFGS(
            atoms,
            maxstep=optimizer_maxstep,
            trajectory=trajectory,
            logfile=logfile,
        )
        optimizer_reported_convergence = optimizer.run(
            fmax=fmax_ev_a,
            steps=max_steps,
        )
        # Reapply the holonomic projection at the release boundary. This only
        # closes ordinary sequential-projector drift; shared-atom pairs should
        # use simultaneous Sella constraint equations instead.
        for _ in range(20):
            atoms.set_positions(atoms.positions.copy())
            actual = np.asarray([atoms.get_distance(i, j) for i, j in pairs])
            if float(np.max(np.abs(actual - targets_array))) <= distance_tolerance_a:
                break
        projected_forces = np.asarray(atoms.get_forces(), dtype=float)
    else:
        if cluster.frozen_indices:
            atoms.set_constraint(FixAtoms(indices=cluster.frozen_indices))
        internal_constraints = Constraints(atoms)
        for pair, target in zip(pairs, targets, strict=True):
            internal_constraints.fix_bond(pair, target=target)
        optimizer = Sella(
            atoms,
            order=0,
            internal=False,
            constraints=internal_constraints,
            constraints_tol=distance_tolerance_a,
            trajectory=trajectory,
            logfile=logfile if logfile is not None else "/dev/null",
        )
        optimizer_reported_convergence = optimizer.run(
            fmax=fmax_ev_a,
            steps=max_steps,
        )
        projected_forces = np.asarray(optimizer.pes.get_projected_forces(), dtype=float)
    projected_fmax = float(np.linalg.norm(projected_forces, axis=1).max())
    if not optimizer_reported_convergence:
        raise RuntimeError(
            f"fixed-distance relaxation did not converge within {max_steps} steps"
        )

    actual = np.asarray([atoms.get_distance(i, j) for i, j in pairs])
    residual = float(np.max(np.abs(actual - np.asarray(targets))))
    if residual > distance_tolerance_a:
        raise RuntimeError(
            "fixed-distance relaxation violated a target by "
            f"{residual:.3e} A (limit {distance_tolerance_a:.3e} A)"
        )
    if projected_fmax >= fmax_ev_a:
        raise RuntimeError(
            "fixed-distance relaxation reported convergence with projected "
            f"fmax {projected_fmax:.6f} eV/A"
        )
    return replace(
        cluster,
        coords=atoms.positions.copy(),
        name=f"{cluster.name}-fixed-distance-relaxed",
    )


def _require_single_imaginary_mode(
    freq: FrequencyResult, *, name: str, noise_floor_cm: float
) -> None:
    """Raise unless exactly one imaginary mode clears the noise floor.

    Shared by :func:`verify_ts` and the precomputed-``FrequencyResult`` path
    in :func:`quick_irc` so the reuse path cannot drift from the main TS
    verification behavior.
    """
    significant = freq.imaginary_cm[freq.imaginary_cm > noise_floor_cm]
    if significant.size != 1:
        raise RuntimeError(
            f"{name}: expected exactly 1 imaginary mode "
            f"above {noise_floor_cm:.0f} cm^-1, found {significant.size} "
            f"(imaginary: {np.round(freq.imaginary_cm, 1).tolist()} cm^-1)"
        )


def verify_ts(
    cluster: Cluster, settings: DftSettings, *, noise_floor_cm: float = 0.0
) -> FrequencyResult:
    """Frequency-verify a saddle: exactly one imaginary mode or raise.

    ``noise_floor_cm`` ignores imaginary modes below the threshold —
    partial-Hessian analyses on peripherally-frozen clusters can carry a
    few numerically-imaginary soft modes that are not reaction modes.
    The default (0) keeps the strict free-cluster behavior.
    """
    freq = frequencies(cluster, settings)
    _require_single_imaginary_mode(
        freq, name=cluster.name, noise_floor_cm=noise_floor_cm
    )
    return freq


def quick_irc(
    ts: Cluster,
    settings: DftSettings,
    *,
    displacement_a: float = 0.15,
    max_steps: int = 200,
    frequency: FrequencyResult | None = None,
    noise_floor_cm: float = 0.0,
) -> tuple[Cluster, Cluster]:
    """Displace ± along the imaginary mode and relax to the two minima.

    The cheap stand-in for a full IRC: confirms which basins the saddle
    connects. Returns (backward, forward) relaxed clusters.
    """
    freq = (
        frequency
        if frequency is not None
        else verify_ts(ts, settings, noise_floor_cm=noise_floor_cm)
    )
    if frequency is not None:
        if freq.geometry_fingerprint != frequency_geometry_fingerprint(ts):
            raise ValueError("precomputed frequency belongs to a different geometry")
        if freq.settings_fingerprint != frequency_settings_fingerprint(settings):
            raise ValueError("precomputed frequency used different DFT settings")
        _require_single_imaginary_mode(
            freq, name=ts.name, noise_floor_cm=noise_floor_cm
        )
    mode = freq.imaginary_mode
    if mode is None:
        raise RuntimeError(f"{ts.name}: no imaginary-mode vector available")
    if mode.shape != ts.coords.shape:
        raise ValueError("imaginary-mode shape does not match TS coordinates")
    if not np.all(np.isfinite(mode)) or float(np.linalg.norm(mode)) < 1e-12:
        raise ValueError("imaginary-mode vector is non-finite or zero")
    step = displacement_a * mode / np.linalg.norm(mode)
    ends = []
    for sign, tag in ((-1.0, "back"), (+1.0, "fwd")):
        displaced = replace(ts, coords=ts.coords + sign * step, name=f"{ts.name}-{tag}")
        ends.append(optimize(displaced, settings, max_steps=max_steps))
    return ends[0], ends[1]


def full_irc(
    ts: Cluster,
    settings: DftSettings,
    *,
    fmax_ev_a: float = 0.05,
    fmax_inner_ev_a: float = 0.01,
    max_steps: int = 400,
    step_size_a: float = 0.10,
    trajectory: str | Path | None = None,
    logfile: str | Path | None = None,
) -> tuple[Cluster, Cluster]:
    """Trace both Gonzalez--Schlegel IRC directions to true minima.

    One :class:`sella.IRC` instance is deliberately reused for ``forward`` and
    ``reverse``. Sella caches the transition-state Hessian and restores it when
    the second direction starts; constructing two unrelated optimizers would
    lose that shared path identity. Both directions are bounded and must report
    convergence. The returned clusters preserve the input atom order,
    electronic state, and frozen-atom contract.
    """

    from ase import Atoms
    from ase.constraints import FixAtoms
    from sella import IRC

    if max_steps <= 0:
        raise ValueError("full IRC max_steps must be positive")
    if fmax_ev_a <= 0.0 or fmax_inner_ev_a <= 0.0:
        raise ValueError("full IRC force tolerances must be positive")
    if step_size_a <= 0.0:
        raise ValueError("full IRC step size must be positive")

    atoms = Atoms(symbols=ts.symbols, positions=ts.coords)
    atoms.calc = make_ase_calculator(settings, ts.charge, ts.spin)
    if ts.frozen_indices:
        atoms.set_constraint(FixAtoms(indices=ts.frozen_indices))

    irc = IRC(
        atoms,
        trajectory=str(trajectory) if trajectory is not None else None,
        logfile=str(logfile) if logfile is not None else "-",
        dx=step_size_a,
    )
    # Sella 2.5 implements its IRC stop rule in ``converged()``, but ASE 3.29
    # changed Optimizer.irun() to call ``gradient_converged()`` directly.  Left
    # unbridged, a force-converged transition state returns success at step zero
    # in both directions without ever taking the mass-weighted IRC kick.
    if hasattr(irc, "gradient_converged") and hasattr(irc, "converged"):

        def gradient_converged(method: Any, gradient: np.ndarray) -> bool:
            return bool(method.converged())

        irc.gradient_converged = MethodType(gradient_converged, irc)
    endpoints: list[Cluster] = []
    for direction, suffix in (("forward", "irc-forward"), ("reverse", "irc-reverse")):
        converged = irc.run(
            fmax=fmax_ev_a,
            fmax_inner=fmax_inner_ev_a,
            steps=max_steps,
            direction=direction,
        )
        if not converged:
            raise RuntimeError(
                f"full IRC {direction} direction did not converge within "
                f"{max_steps} steps"
            )
        endpoints.append(
            replace(ts, coords=atoms.positions.copy(), name=f"{ts.name}-{suffix}")
        )
    return endpoints[0], endpoints[1]


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
# A single constrained optimization can hop into a bad peripheral-hydroxyl
# conformer while its neighbors remain on the reaction path.  The live
# al-neutral scan produced a 117 kJ/mol spike followed by an 84 kJ/mol drop.
# A genuine crest that both neighboring seeds reproduce is retained; only the
# one-direction path-dependent result is replaced.
_SCAN_OUTLIER_THRESHOLD_KJ = 40.0


def _scan_outlier_indices(
    scan: list[tuple[float, float, Cluster]],
    *,
    threshold_kj: float = _SCAN_OUTLIER_THRESHOLD_KJ,
) -> list[int]:
    """Return interior points far above linear neighbor interpolation."""
    threshold_hartree = threshold_kj / HARTREE_TO_KJ
    outliers = []
    for i in range(1, len(scan) - 1):
        r_left, e_left, _ = scan[i - 1]
        r, energy, _ = scan[i]
        r_right, e_right, _ = scan[i + 1]
        if r_right == r_left:
            continue
        fraction = (r - r_left) / (r_right - r_left)
        interpolated = e_left + fraction * (e_right - e_left)
        if energy > interpolated + threshold_hartree:
            outliers.append(i)
    return outliers


def _repair_scan_outliers(
    scan: list[tuple[float, float, Cluster]],
    settings: DftSettings,
    *,
    atom_i: int,
    atom_j: int,
    fixed_distances: list[tuple[int, int, float]],
    max_steps: int,
    attempted: set[int],
    progress=None,
) -> None:
    """Reoptimize each new spike from both neighboring scan structures.

    ``attempted`` persists while a scan is extended.  A point reproduced as
    high energy from both directions may be a real sharp crest, so it is
    checked once and then left for crest detection rather than retried forever.
    """
    while True:
        indices = [i for i in _scan_outlier_indices(scan) if i not in attempted]
        if not indices:
            return
        for i in indices:
            attempted.add(i)
            r = scan[i][0]
            retries = []
            for seed in (scan[i - 1][2], scan[i + 1][2]):
                result = constrained_scan(
                    seed,
                    settings,
                    atom_i=atom_i,
                    atom_j=atom_j,
                    distances_a=[r],
                    fixed_distances=fixed_distances,
                    max_steps=max_steps,
                )[0]
                retries.append(result)
                if progress:
                    progress(result[0], result[1])
            scan[i] = min(retries, key=lambda point: point[1])


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


_NEB_CHECKPOINT_SCHEMA = 1


def _write_neb_checkpoint(
    images,
    checkpoint_root: Path,
    *,
    stage: str,
    converged: bool | None,
    settings: DftSettings,
    charge: int,
    spin: int,
    frozen_indices: list[int] | None,
    optimizer: str,
    fmax_ev_a: float,
    step_bound: int,
) -> Path:
    """Atomically persist every image plus the exact NEB stage contract.

    Checkpoints are immutable timestamped directories. ``latest.json`` is only
    an atomic pointer, so a later failed attempt cannot destroy an earlier
    pre-relaxed band. Energies are deliberately omitted: asking for them here
    could add expensive calculator calls merely to write a crash receipt.
    """
    checkpoint_root.mkdir(parents=True, exist_ok=True)
    checkpoint_name = f"{time.time_ns()}-{stage}"
    temporary = checkpoint_root / f".{checkpoint_name}.tmp"
    final = checkpoint_root / checkpoint_name
    temporary.mkdir()
    symbols = list(images[0].get_chemical_symbols())
    image_files = []
    for index, image in enumerate(images):
        if list(image.get_chemical_symbols()) != symbols:
            raise ValueError("NEB checkpoint images have inconsistent atom mappings")
        filename = f"image-{index:02d}.xyz"
        lines = [str(len(symbols)), f"NEB {stage} image {index}"]
        lines.extend(
            f"{symbol} {x:.16f} {y:.16f} {z:.16f}"
            for symbol, (x, y, z) in zip(symbols, image.positions, strict=True)
        )
        (temporary / filename).write_text("\n".join(lines) + "\n")
        image_files.append(filename)
    manifest = {
        "schema_version": _NEB_CHECKPOINT_SCHEMA,
        "stage": stage,
        "converged": converged,
        "created_at_ns": time.time_ns(),
        "settings": asdict(settings),
        "charge": charge,
        "spin": spin,
        "frozen_indices": list(frozen_indices or []),
        "symbols": symbols,
        "n_images": len(images),
        "optimizer": optimizer,
        "fmax_ev_a": fmax_ev_a,
        "step_bound": step_bound,
        "images": image_files,
    }
    (temporary / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    temporary.replace(final)
    latest_tmp = checkpoint_root / ".latest.json.tmp"
    latest_tmp.write_text(json.dumps({"checkpoint": checkpoint_name}, indent=2) + "\n")
    latest_tmp.replace(checkpoint_root / "latest.json")
    return final


def _load_neb_checkpoint(
    checkpoint: Path,
    *,
    settings: DftSettings,
    reactant: Cluster,
    n_images: int,
) -> tuple[list[np.ndarray], dict]:
    """Load and validate a complete pre-relaxed band for climb-only resume."""
    if (checkpoint / "latest.json").is_file():
        latest = json.loads((checkpoint / "latest.json").read_text())
        checkpoint = checkpoint / latest["checkpoint"]
    manifest = json.loads((checkpoint / "manifest.json").read_text())
    expected = {
        "schema_version": _NEB_CHECKPOINT_SCHEMA,
        "settings": asdict(settings),
        "charge": reactant.charge,
        "spin": reactant.spin,
        "frozen_indices": list(reactant.frozen_indices or []),
        "symbols": list(reactant.symbols),
        "n_images": n_images,
    }
    for key, value in expected.items():
        if manifest.get(key) != value:
            raise ValueError(f"NEB checkpoint {key} does not match this run")
    stage = str(manifest.get("stage", ""))
    terminal_stage = stage in {"pre-relax-final", "climb-final"} and (
        manifest.get("converged") is True
    )
    atomic_climb_step = stage.startswith("climb-step-") and (
        manifest.get("converged") is None
    )
    if not terminal_stage and not atomic_climb_step:
        raise ValueError(
            "NEB resume requires a converged stage boundary or atomic climb-step "
            "checkpoint"
        )
    coords = []
    for filename in manifest["images"]:
        lines = (checkpoint / filename).read_text().splitlines()
        if int(lines[0]) != len(reactant.symbols):
            raise ValueError("NEB checkpoint image atom count does not match")
        image_symbols = [line.split()[0] for line in lines[2:]]
        if image_symbols != reactant.symbols:
            raise ValueError("NEB checkpoint image atom order does not match")
        xyz = np.asarray(
            [[float(value) for value in line.split()[1:4]] for line in lines[2:]],
            dtype=float,
        )
        if xyz.shape != reactant.coords.shape or not np.all(np.isfinite(xyz)):
            raise ValueError("NEB checkpoint image coordinates are invalid")
        coords.append(xyz)
    return coords, manifest


def neb_ts_guess(
    reactant: Cluster,
    product: Cluster,
    settings: DftSettings,
    *,
    n_images: int = 7,
    fmax_ev_a: float = 0.10,
    max_steps: int = 150,
    pre_relax_fmax_ev_a: float = 0.30,
    pre_relax_steps: int = 80,
    optimizer_dt: float = 0.05,
    optimizer_maxstep: float = 0.05,
    climb_optimizer: str = "mdmin",
    checkpoint_dir: str | Path | None = None,
    checkpoint_interval: int = 5,
    resume_from: str | Path | None = None,
    initial_image_coords: list[np.ndarray] | None = None,
) -> Cluster:
    """Climbing-image NEB between two basins; returns the peak image.

    The remedy for saddle searches that escape the reaction channel
    (observed live: Sella walked a crest guess from r(Si-Ow)=1.90 back
    out to 3.22 A): the band is pinned to both basins so it cannot
    escape, and the climbing image rides the path to the col, where the
    lowest curvature direction *is* the reaction mode. The peak image
    is then a Sella start it can refine rather than flee.

    ``fmax_ev_a`` is deliberately loose — this produces a guess, and
    Sella does the tight convergence.  The band is first relaxed without a
    climbing image, then the climb is enabled.  Pre-relaxation uses bounded
    MDMin. The climbing stage uses MDMin by default, ASE's bounded ODE
    optimizer, or BFGS. ODE and BFGS are the measured fallbacks for a live
    al-neutral associative band whose MDMin climb stalled. Both stages must
    converge; returning an arbitrary peak from an exhausted optimizer is
    forbidden. When ``checkpoint_dir`` is set, every image and the optimizer
    contract are persisted periodically and at each stage boundary. A
    ``pre-relax-final`` checkpoint with ``converged=true`` can be supplied via
    ``resume_from`` to skip interpolation and pre-relaxation safely.
    ``initial_image_coords`` supplies a complete physically seeded band; it is
    mutually exclusive with resume and still undergoes bounded pre-relaxation.
    """
    from ase import Atoms
    from ase.build.rotate import minimize_rotation_and_translation

    try:
        from ase.mep import NEB
    except ImportError:  # older ASE layout
        from ase.neb import NEB
    from ase.optimize import BFGS, MDMin

    if reactant.charge != product.charge or reactant.spin != product.spin:
        raise ValueError("reactant/product electronic states differ")
    if reactant.symbols != product.symbols:
        raise ValueError("reactant/product atom order differs")
    if reactant.frozen_indices != product.frozen_indices:
        raise ValueError("reactant/product frozen atom sets differ")
    if climb_optimizer not in {"mdmin", "ode", "bfgs"}:
        raise ValueError(f"unknown NEB climb optimizer '{climb_optimizer}'")
    if n_images < 3:
        raise ValueError("NEB requires at least three images")
    if checkpoint_interval <= 0:
        raise ValueError("checkpoint_interval must be positive")
    if resume_from is not None and initial_image_coords is not None:
        raise ValueError("resume_from and initial_image_coords are mutually exclusive")
    checkpoint_root = Path(checkpoint_dir) if checkpoint_dir is not None else None

    images = [Atoms(symbols=reactant.symbols, positions=reactant.coords)]
    for _ in range(n_images - 2):
        images.append(images[0].copy())
    images.append(Atoms(symbols=product.symbols, positions=product.coords))
    # Independent optimizations can rigidly rotate/translate the product.
    # Interpolating those raw Cartesian frames sent the transferring proton
    # several Angstroms into vacuum and made the first DFT SCF fail. Remove
    # only that physically meaningless global motion before IDPP.
    # A frozen shell already pins both endpoints to one common frame —
    # rigid re-fitting there would *move* the frozen atoms, so skip it.
    frozen = bool(reactant.frozen_indices)
    resumed_manifest = None
    if resume_from is not None:
        resumed, resumed_manifest = _load_neb_checkpoint(
            Path(resume_from),
            settings=settings,
            reactant=reactant,
            n_images=n_images,
        )
        for image, coords in zip(images, resumed, strict=True):
            image.positions[:] = coords
    elif initial_image_coords is not None:
        if len(initial_image_coords) != n_images:
            raise ValueError("initial_image_coords length does not match n_images")
        seeded = [np.asarray(coords, dtype=float) for coords in initial_image_coords]
        if any(coords.shape != reactant.coords.shape for coords in seeded):
            raise ValueError("initial_image_coords contains an invalid shape")
        if any(not np.all(np.isfinite(coords)) for coords in seeded):
            raise ValueError("initial_image_coords contains non-finite coordinates")
        if not np.allclose(seeded[0], reactant.coords, atol=1e-8):
            raise ValueError("initial_image_coords first image is not the reactant")
        if not np.allclose(seeded[-1], product.coords, atol=1e-8):
            raise ValueError("initial_image_coords last image is not the product")
        for image, coords in zip(images, seeded, strict=True):
            image.positions[:] = coords
    elif not frozen:
        minimize_rotation_and_translation(images[0], images[-1])
    for image in images:
        image.calc = make_ase_calculator(settings, reactant.charge, reactant.spin)
        if frozen:
            from ase.constraints import FixAtoms

            image.set_constraint(FixAtoms(indices=reactant.frozen_indices))
    neb = NEB(
        images,
        climb=False,
        method="improvedtangent",
        remove_rotation_and_translation=not frozen,
    )
    if resumed_manifest is not None or initial_image_coords is not None:
        pass
    elif frozen:
        # ASE refuses to interpolate constrained images without an explicit
        # choice; applying is exact here (identical frozen coordinates).
        neb.interpolate(method="idpp", apply_constraint=True)
    else:
        neb.interpolate(method="idpp")

    if checkpoint_root is not None:
        _write_neb_checkpoint(
            images,
            checkpoint_root,
            stage="initialized",
            converged=None,
            settings=settings,
            charge=reactant.charge,
            spin=reactant.spin,
            frozen_indices=reactant.frozen_indices,
            optimizer="none",
            fmax_ev_a=pre_relax_fmax_ev_a,
            step_bound=pre_relax_steps,
        )

    if resumed_manifest is None:
        pre_relax = MDMin(
            neb,  # type: ignore[arg-type] -- ASE optimizers accept NEB objects
            dt=optimizer_dt,
            maxstep=optimizer_maxstep,
        )
        if checkpoint_root is not None and hasattr(pre_relax, "attach"):
            pre_relax.attach(
                lambda: _write_neb_checkpoint(
                    images,
                    checkpoint_root,
                    stage=f"pre-relax-step-{pre_relax.nsteps:06d}",
                    converged=None,
                    settings=settings,
                    charge=reactant.charge,
                    spin=reactant.spin,
                    frozen_indices=reactant.frozen_indices,
                    optimizer="mdmin",
                    fmax_ev_a=pre_relax_fmax_ev_a,
                    step_bound=pre_relax_steps,
                ),
                interval=checkpoint_interval,
            )
        relaxed = pre_relax.run(fmax=pre_relax_fmax_ev_a, steps=pre_relax_steps)
        if checkpoint_root is not None:
            _write_neb_checkpoint(
                images,
                checkpoint_root,
                stage="pre-relax-final",
                converged=bool(relaxed),
                settings=settings,
                charge=reactant.charge,
                spin=reactant.spin,
                frozen_indices=reactant.frozen_indices,
                optimizer="mdmin",
                fmax_ev_a=pre_relax_fmax_ev_a,
                step_bound=pre_relax_steps,
            )
        if not relaxed:
            raise RuntimeError(
                "NEB pre-relaxation did not converge "
                f"to {pre_relax_fmax_ev_a:.3f} eV/A within {pre_relax_steps} steps"
            )

    neb.climb = True
    if climb_optimizer == "mdmin":
        climb = MDMin(
            neb,  # type: ignore[arg-type] -- ASE optimizers accept NEB objects
            dt=optimizer_dt,
            maxstep=optimizer_maxstep,
        )
    elif climb_optimizer == "ode":
        from ase.mep.neb import NEBOptimizer

        climb = NEBOptimizer(
            neb,  # type: ignore[arg-type] -- ASE accepts NEB objects
            method="ODE",
        )
    else:
        climb = BFGS(
            neb,  # type: ignore[arg-type] -- ASE optimizers accept NEB objects
            maxstep=optimizer_maxstep,
        )
    if checkpoint_root is not None and hasattr(climb, "attach"):
        climb.attach(
            lambda: _write_neb_checkpoint(
                images,
                checkpoint_root,
                stage=f"climb-step-{climb.nsteps:06d}",
                converged=None,
                settings=settings,
                charge=reactant.charge,
                spin=reactant.spin,
                frozen_indices=reactant.frozen_indices,
                optimizer=climb_optimizer,
                fmax_ev_a=fmax_ev_a,
                step_bound=max_steps,
            ),
            interval=checkpoint_interval,
        )
    climbed = climb.run(fmax=fmax_ev_a, steps=max_steps)
    if checkpoint_root is not None:
        _write_neb_checkpoint(
            images,
            checkpoint_root,
            stage="climb-final",
            converged=bool(climbed),
            settings=settings,
            charge=reactant.charge,
            spin=reactant.spin,
            frozen_indices=reactant.frozen_indices,
            optimizer=climb_optimizer,
            fmax_ev_a=fmax_ev_a,
            step_bound=max_steps,
        )
    if not climbed:
        raise RuntimeError(
            "climbing-image NEB did not converge "
            f"to {fmax_ev_a:.3f} eV/A within {max_steps} steps"
        )
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
    repaired_indices: set[int] = set()
    _repair_scan_outliers(
        scan,
        settings,
        atom_i=atom_i,
        atom_j=atom_j,
        fixed_distances=fixed_distances,
        max_steps=max_steps,
        attempted=repaired_indices,
        progress=progress,
    )
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
        _repair_scan_outliers(
            scan,
            settings,
            atom_i=atom_i,
            atom_j=atom_j,
            fixed_distances=fixed_distances,
            max_steps=max_steps,
            attempted=repaired_indices,
            progress=progress,
        )
    return scan
