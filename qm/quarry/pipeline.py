"""The tiered barrier protocol as composable steps (SURVEY.md §6.1, §8).

Full ladder: xtb screen -> composite/hybrid DFT opt+freq -> geomeTRIC/
Sella TS + IRC verification -> high-level single points (+ solvation).
Phase 0 implements the PySCF-backed rungs (energy, gradient, opt, freq)
against the ``Cluster`` contract, engine-selectable so GPU4PySCF slots
in on the workstation without touching call sites. TS search and IRC
land with Phase 1; xtb screening is optional-detected.

Everything here is CPU-testable on tiny molecules at small basis — the
default test gate exercises each rung on H2O/Si(OH)4.
"""

from __future__ import annotations

import hashlib
import os
import tempfile
import warnings
from contextlib import contextmanager
from dataclasses import dataclass, replace
from typing import Any

import numpy as np

from quarry.clusters import Cluster

HARTREE_TO_KJ = 2625.4996394798254  # CODATA 2018
BOHR_TO_ANGSTROM = 0.529177210903


@dataclass(frozen=True)
class DftSettings:
    """One rung's electronic-structure settings.

    ``xc = "hf"`` runs Hartree-Fock (cheap tests); anything else is a
    libxc functional name. ``solvent`` of ``None`` is gas phase; "smd" /
    "pcm" attach the corresponding implicit model (water parameters).
    """

    xc: str = "b3lyp"
    basis: str = "def2-svp"
    solvent: str | None = None
    dispersion: str | None = None  # e.g. "d3bj", "d4" (needs pyscf[dispersion])
    grid_level: int | None = None
    # RI-JK. Off here so cheap tests and comparisons are explicit;
    # production call sites and the smoke benchmark turn it on
    # (SURVEY §2.2, §3.2).
    density_fit: bool = False
    use_gpu: bool = False


def build_mol(cluster: Cluster, settings: DftSettings) -> Any:
    """A pyscf Mole for the cluster (Angstrom in, pyscf defaults out)."""
    from pyscf import gto

    return gto.M(
        atom=cluster.to_pyscf_atom(),
        charge=cluster.charge,
        spin=cluster.spin,
        basis=settings.basis,
        verbose=0,
    )


def _make_scf(mol: Any, settings: DftSettings) -> Any:
    from pyscf import dft, scf

    restricted = mol.spin == 0
    if settings.xc.lower() == "hf":
        mf = scf.RHF(mol) if restricted else scf.UHF(mol)
    else:
        mf = dft.RKS(mol) if restricted else dft.UKS(mol)
        mf.xc = settings.xc
        if settings.grid_level is not None:
            mf.grids.level = settings.grid_level
    # Hand-built cluster guesses and stretched scan geometries need more
    # than the default 50 cycles; DIIS usually gets there given room.
    mf.max_cycle = 150
    if settings.dispersion is not None:
        mf.disp = settings.dispersion
    if settings.density_fit:
        mf = mf.density_fit()
    if settings.solvent is not None:
        mf = mf.SMD() if settings.solvent.lower() == "smd" else mf.PCM()
    if settings.use_gpu:
        mf = mf.to_gpu()
    return mf


def energy(cluster: Cluster, settings: DftSettings) -> float:
    """Single-point energy, Hartree."""
    mf = _make_scf(build_mol(cluster, settings), settings)
    e = mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"SCF did not converge for {cluster.name}")
    return float(e)


def gradient(cluster: Cluster, settings: DftSettings) -> np.ndarray:
    """Nuclear gradient, Hartree/Bohr, shape (n, 3)."""
    mf = _make_scf(build_mol(cluster, settings), settings)
    mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"SCF did not converge for {cluster.name}")
    return np.asarray(mf.nuc_grad_method().kernel())


@contextmanager
def constraints_file(text: str):
    """geomeTRIC takes constraints as a *file path*, not inline text
    (passing text raises FileNotFoundError on the constraint string —
    found live in the phase-1 scan). Yields a temp-file path."""
    fd, path = tempfile.mkstemp(suffix=".constraints", text=True)
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(text)
        yield path
    finally:
        os.unlink(path)


def optimize(
    cluster: Cluster, settings: DftSettings, *, max_steps: int = 100
) -> Cluster:
    """Minimum-energy geometry via geomeTRIC (frozen atoms respected)."""
    from pyscf.geomopt.geometric_solver import optimize as geometric_optimize

    mf = _make_scf(build_mol(cluster, settings), settings)
    if cluster.frozen_indices:
        # geomeTRIC constraint block: freeze xyz of the peripheral atoms
        # (the lattice-resistance contract, SURVEY.md §6.2).
        atoms = ",".join(str(i + 1) for i in sorted(cluster.frozen_indices))
        with constraints_file(f"$freeze\nxyz {atoms}\n") as path:
            mol_opt = geometric_optimize(mf, maxsteps=max_steps, constraints=path)
    else:
        mol_opt = geometric_optimize(mf, maxsteps=max_steps)
    return replace(
        cluster,
        symbols=[mol_opt.atom_symbol(i) for i in range(mol_opt.natm)],
        coords=mol_opt.atom_coords() * BOHR_TO_ANGSTROM,
    )


@dataclass(frozen=True)
class FrequencyResult:
    """Projected harmonic frequencies (cm^-1) and the pieces rates.py needs."""

    frequencies_cm: np.ndarray  # real (positive) modes only
    imaginary_cm: np.ndarray  # magnitudes of imaginary modes
    electronic_hartree: float
    molar_mass_kg: float
    rotational_temperatures_k: tuple[float, ...] | None
    linear: bool
    # Cartesian displacement of the largest imaginary mode, shape (n, 3);
    # None when there is no imaginary mode. This is the reaction mode a
    # TS hands to quick-IRC and to the tunneling correction.
    imaginary_mode: np.ndarray | None = None
    imaginary_modes: np.ndarray | None = None
    geometry_fingerprint: str | None = None
    settings_fingerprint: str | None = None

    @property
    def n_imaginary(self) -> int:
        return int(self.imaginary_cm.size)


def frequency_geometry_fingerprint(cluster: Cluster) -> str:
    """Exact geometry/electronic-state identity for reusable Hessian results.

    Includes the frozen-atom set: PHVA (partial Hessian) and TS verification
    depend on which atoms are frozen, so a Hessian computed with a different
    frozen set must never be reused for this geometry.
    """
    digest = hashlib.sha256()
    digest.update("\0".join(cluster.symbols).encode())
    digest.update(f"\0{cluster.charge}\0{cluster.spin}\0".encode())
    coords = np.ascontiguousarray(cluster.coords, dtype="<f8")
    digest.update(coords.tobytes())
    frozen = ",".join(str(i) for i in sorted(cluster.frozen_indices))
    digest.update(f"\0frozen:{frozen}\0".encode())
    return digest.hexdigest()


def frequency_settings_fingerprint(settings: DftSettings) -> str:
    """Scientific settings identity, independent of CPU/GPU execution backend."""
    return repr(replace(settings, use_gpu=False))


def frequencies(cluster: Cluster, settings: DftSettings) -> FrequencyResult:
    """Analytic Hessian -> harmonic analysis.

    Free clusters get the standard projected analysis (pyscf thermo).
    Clusters with ``frozen_indices`` get a partial-Hessian vibrational
    analysis (PHVA) over the free-atom block instead: no trans/rot
    projection (the frozen lattice removes those motions) and no
    rotational temperatures — the embedded cluster does not rotate.
    Translational terms downstream cancel in any barrier between states
    of identical composition, which is the only use PHVA results have.
    """
    from pyscf.hessian import thermo as pyscf_thermo

    mf, e, hess = _scf_hessian(cluster, settings)
    hess = np.asarray(hess.get() if hasattr(hess, "get") else hess)
    if cluster.frozen_indices:
        return _partial_hessian_analysis(cluster, settings, mf, hess, float(e))
    freq_info = pyscf_thermo.harmonic_analysis(mf.mol, hess)
    nu = np.asarray(freq_info["freq_wavenumber"])
    modes = np.asarray(freq_info["norm_mode"])  # (nmodes, natm, 3)
    imag_mode = None
    imag_modes = None
    if np.iscomplexobj(nu):
        is_imag = np.abs(np.imag(nu)) > 0
        imag = np.abs(np.imag(nu[is_imag]))
        real = np.real(nu[~is_imag])
        if imag.size:
            order = np.argsort(imag)
            imag = imag[order]
            imag_modes = np.real(modes[is_imag])[order]
            # Backward-compatible default: largest-magnitude imaginary mode.
            imag_mode = imag_modes[-1]
    else:
        imag = np.zeros(0)
        real = nu
    real = real[real > 0]

    mol = mf.mol
    mass_kg = float(np.sum(mol.atom_mass_list())) * 1e-3  # amu/mol -> kg/mol
    rot = _rotational_temperatures(cluster)
    return FrequencyResult(
        frequencies_cm=np.sort(real),
        imaginary_cm=imag,
        electronic_hartree=float(e),
        molar_mass_kg=mass_kg,
        rotational_temperatures_k=rot[0],
        linear=rot[1],
        imaginary_mode=imag_mode,
        imaginary_modes=imag_modes,
        geometry_fingerprint=frequency_geometry_fingerprint(cluster),
        settings_fingerprint=frequency_settings_fingerprint(settings),
    )


def _partial_hessian_analysis(
    cluster: Cluster,
    settings: DftSettings,
    mf: Any,
    hess: np.ndarray,
    e: float,
) -> FrequencyResult:
    """PHVA: diagonalize the mass-weighted free-atom Hessian block."""
    from pyscf.data import elements, nist

    n = len(cluster.symbols)
    free = sorted(set(range(n)) - set(cluster.frozen_indices))
    masses = np.array(
        [elements.MASSES[elements.charge(cluster.symbols[i])] for i in free]
    )
    sub = hess[np.ix_(free, free)]  # (nf, nf, 3, 3)
    nf = len(free)
    h = sub.transpose(0, 2, 1, 3).reshape(3 * nf, 3 * nf)
    h = 0.5 * (h + h.T)
    inv_sqrt_m = 1.0 / np.sqrt(np.repeat(masses, 3))
    h_mw = h * np.outer(inv_sqrt_m, inv_sqrt_m)  # Hartree/(Bohr^2 amu)
    eigvals, eigvecs = np.linalg.eigh(h_mw)
    # Hartree/(Bohr^2 amu) -> s^-2 -> cm^-1.
    to_si = nist.HARTREE2J / (nist.ATOMIC_MASS * (nist.BOHR_SI**2))
    nu_cm = np.sqrt(np.abs(eigvals) * to_si) / (2.0 * np.pi * nist.LIGHT_SPEED_SI * 100)
    imag_mask = eigvals < 0
    imag = np.asarray(nu_cm[imag_mask])
    real = np.asarray(nu_cm[~imag_mask])
    imag_mode = None
    imag_modes = None
    if imag.size:
        indices = np.flatnonzero(imag_mask)
        order = np.argsort(imag)
        imag = imag[order]
        full_modes = []
        for idx in indices[order]:
            vec = (eigvecs[:, idx] * inv_sqrt_m).reshape(nf, 3)
            full = np.zeros((n, 3))
            full[free] = vec
            full_modes.append(full)
        imag_modes = np.asarray(full_modes)
        imag_mode = imag_modes[-1]

    mass_kg = float(np.sum(mf.mol.atom_mass_list())) * 1e-3
    return FrequencyResult(
        frequencies_cm=np.sort(real[real > 0]),
        imaginary_cm=imag,
        electronic_hartree=e,
        molar_mass_kg=mass_kg,
        rotational_temperatures_k=None,
        linear=False,
        imaginary_mode=imag_mode,
        imaginary_modes=imag_modes,
        geometry_fingerprint=frequency_geometry_fingerprint(cluster),
        settings_fingerprint=frequency_settings_fingerprint(settings),
    )


def _is_gpu_hessian_contiguity_assertion(err: AssertionError) -> bool:
    """True only for the known GPU4PySCF DF-UKS Hessian contiguity defect.

    GPU4PySCF 1.8.1 aborts inside ``gpu4pyscf.hessian.uks`` on
    ``assert wv.flags.c_contiguous`` for open-shell molecules.  We match on the
    assertion's message or its origin (a gpu4pyscf frame) and re-raise anything
    else, so genuine programming errors are never masked by the CPU fallback.
    """
    if "contiguous" in str(err).lower():
        return True
    tb = err.__traceback__
    while tb is not None:
        if "gpu4pyscf" in tb.tb_frame.f_code.co_filename:
            return True
        tb = tb.tb_next
    return False


def _scf_hessian(cluster: Cluster, settings: DftSettings) -> tuple[Any, float, Any]:
    """Run SCF + Hessian, with a bounded CPU fallback for GPU backend defects.

    GPU4PySCF 1.8.1's density-fitted UKS Hessian can raise an internal
    C-contiguity assertion for open-shell molecules.  Geometry and saddle
    work remains GPU-first; only the tiny-molecule analytic Hessian retries on
    CPU when that exact backend assertion occurs.
    """
    mf = _make_scf(build_mol(cluster, settings), settings)
    e = mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"SCF did not converge for {cluster.name}")
    try:
        return mf, float(e), mf.Hessian().kernel()
    except AssertionError as gpu_error:
        if not settings.use_gpu or not _is_gpu_hessian_contiguity_assertion(gpu_error):
            raise
        warnings.warn(
            "GPU4PySCF Hessian assertion; retrying this Hessian on CPU",
            RuntimeWarning,
            stacklevel=2,
        )
        cpu_settings = replace(settings, use_gpu=False)
        cpu_mf = _make_scf(build_mol(cluster, cpu_settings), cpu_settings)
        cpu_e = cpu_mf.kernel()
        if not cpu_mf.converged:
            raise RuntimeError(
                f"CPU fallback SCF did not converge for {cluster.name}"
            ) from gpu_error
        return cpu_mf, float(cpu_e), cpu_mf.Hessian().kernel()


def _rotational_temperatures(
    cluster: Cluster,
) -> tuple[tuple[float, ...] | None, bool]:
    """(rotational temperatures in K, is_linear) from the geometry."""
    from pyscf.data import elements

    from quarry.rates import KB, H

    masses = np.array(
        [elements.MASSES[elements.charge(s)] for s in cluster.symbols]
    )  # amu
    if len(masses) == 1:
        return None, False
    amu = 1.66053906660e-27
    coords_m = cluster.coords * 1e-10
    com = np.average(coords_m, axis=0, weights=masses)
    r = coords_m - com
    m = masses * amu
    inertia = np.zeros((3, 3))
    for mi, ri in zip(m, r, strict=True):
        inertia += mi * (np.dot(ri, ri) * np.eye(3) - np.outer(ri, ri))
    moments = np.sort(np.linalg.eigvalsh(inertia))  # kg m^2
    hbar = H / (2.0 * np.pi)
    linear = moments[0] / moments[2] < 1e-8
    if linear:
        i_val = moments[2]
        return (float(hbar**2 / (2.0 * i_val * KB)),), True
    thetas = tuple(float(hbar**2 / (2.0 * i_val * KB)) for i_val in moments)
    return thetas, False
