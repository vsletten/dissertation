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


def optimize(
    cluster: Cluster, settings: DftSettings, *, max_steps: int = 100
) -> Cluster:
    """Minimum-energy geometry via geomeTRIC (frozen atoms respected)."""
    from pyscf.geomopt.geometric_solver import optimize as geometric_optimize

    mf = _make_scf(build_mol(cluster, settings), settings)
    params = {}
    if cluster.frozen_indices:
        # geomeTRIC constraint block: freeze xyz of the peripheral atoms
        # (the lattice-resistance contract, SURVEY.md §6.2).
        atoms = ",".join(str(i + 1) for i in sorted(cluster.frozen_indices))
        params["constraints"] = f"$freeze\nxyz {atoms}\n"
    mol_opt = geometric_optimize(mf, maxsteps=max_steps, **params)
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

    @property
    def n_imaginary(self) -> int:
        return int(self.imaginary_cm.size)


def frequencies(cluster: Cluster, settings: DftSettings) -> FrequencyResult:
    """Analytic Hessian -> projected harmonic analysis (pyscf thermo)."""
    from pyscf.hessian import thermo as pyscf_thermo

    mf = _make_scf(build_mol(cluster, settings), settings)
    e = mf.kernel()
    if not mf.converged:
        raise RuntimeError(f"SCF did not converge for {cluster.name}")
    hess = mf.Hessian().kernel()
    freq_info = pyscf_thermo.harmonic_analysis(mf.mol, hess)
    nu = np.asarray(freq_info["freq_wavenumber"])
    imag = np.abs(np.imag(nu)) if np.iscomplexobj(nu) else np.zeros(0)
    imag = imag[imag > 0] if imag.size else np.zeros(0)
    real = np.real(nu[np.isreal(nu)]) if np.iscomplexobj(nu) else nu
    real = real[real > 0]

    mol = mf.mol
    mass_kg = float(np.sum(mol.atom_mass_list())) * 1e-3  # amu/mol -> kg/mol
    rot = _rotational_temperatures(cluster)
    return FrequencyResult(
        frequencies_cm=np.sort(real),
        imaginary_cm=np.sort(imag),
        electronic_hartree=float(e),
        molar_mass_kg=mass_kg,
        rotational_temperatures_k=rot[0],
        linear=rot[1],
    )


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
