"""Thermochemistry -> k(T): the in-house rate module.

Quasi-RRHO Gibbs energies from harmonic frequencies, Eyring k(T) from
ΔG‡, and Wigner / asymmetric-Eckart tunneling corrections. This is the
~50 lines SURVEY.md §4.2 says are worth owning for transparency; the
audited reference implementations are GoodVibes (quasi-RRHO) and Arkane
(TST + tunneling), and the tests gate the closed-form pieces to 1e-12.

Unit conventions (chosen to match the emit targets, HANDOFF.md §2):
energies in kJ/mol unless a name says otherwise, temperatures in K,
frequencies in cm^-1, rate constants in s^-1 (unimolecular convention;
standard-state corrections for bimolecular steps live with the caller).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# CODATA 2018 exact/derived values.
KB = 1.380649e-23  # J/K
H = 6.62607015e-34  # J s
C_CM = 2.99792458e10  # speed of light, cm/s
NA = 6.02214076e23  # 1/mol
R = KB * NA  # J/(mol K)
R_KJ = R / 1000.0  # kJ/(mol K)
KCAL_TO_KJ = 4.184
EV_TO_KJ = 96.48533212331002  # NA * e / 1000

# Grimme's quasi-RRHO free-rotor moment-of-inertia average (JCC 2012).
_BAV = 1.0e-44  # kg m^2


def wavenumber_to_joule(nu_cm: float) -> float:
    """Photon energy h*c*nu of a wavenumber, in J."""
    return H * C_CM * nu_cm


def eyring_k(dg_kj: float, temperature: float) -> float:
    """Eyring TST rate constant from ΔG‡ (kJ/mol): (kB T / h) exp(-ΔG‡/RT)."""
    return (KB * temperature / H) * math.exp(-dg_kj / (R_KJ * temperature))


def eyring_k_dh_ds(dh_kj: float, ds_kj: float, temperature: float) -> float:
    """Eyring k(T) from ΔH‡ (kJ/mol) and ΔS‡ (kJ mol^-1 K^-1).

    The (dh, ds) pair is exactly what a petra `eyring` rate spec carries
    (in deck units), so this is the closed form the emitted decks will
    reproduce at runtime.
    """
    return eyring_k(dh_kj - temperature * ds_kj, temperature)


def wigner_kappa(imag_nu_cm: float, temperature: float) -> float:
    """Wigner tunneling correction: 1 + (1/24) (h nu‡ / kB T)^2."""
    u = wavenumber_to_joule(abs(imag_nu_cm)) / (KB * temperature)
    return 1.0 + u * u / 24.0


def eckart_kappa(
    imag_nu_cm: float,
    dv_forward_kj: float,
    dv_reverse_kj: float,
    temperature: float,
    *,
    n_grid: int = 20_000,
    e_max_factor: float = 40.0,
) -> float:
    """Asymmetric-Eckart tunneling correction (Johnston & Heicklen 1962).

    ``dv_forward_kj`` / ``dv_reverse_kj`` are the zero-point-corrected
    barrier heights (TS minus reactant / product) in kJ/mol; the
    imaginary frequency sets the barrier width. Transmission
    probabilities are Boltzmann-averaged on a fixed energy grid
    (trapezoid; converged well past the default tolerances for the
    frequencies and barriers this project meets).
    """
    if dv_forward_kj <= 0.0 or dv_reverse_kj <= 0.0:
        raise ValueError(
            "Eckart correction needs positive forward and reverse barriers"
        )
    kt = KB * temperature
    # Per-molecule quantities, in J.
    v1 = dv_forward_kj * 1000.0 / NA
    v2 = dv_reverse_kj * 1000.0 / NA
    hnu = wavenumber_to_joule(abs(imag_nu_cm))

    # Dimensionless barrier parameters (Johnston & Heicklen conventions,
    # as in Arkane's Eckart implementation).
    alpha1 = 2.0 * math.pi * v1 / hnu
    alpha2 = 2.0 * math.pi * v2 / hnu
    if alpha1 * alpha2 <= math.pi**2 / 4.0:
        # Barrier too thin/low for the Eckart form; Wigner is the safer bet.
        return wigner_kappa(imag_nu_cm, temperature)
    two_pi_d = 2.0 * math.sqrt(alpha1 * alpha2 - math.pi**2 / 4.0)
    inv_root = 1.0 / math.sqrt(alpha1) + 1.0 / math.sqrt(alpha2)

    # E measured from the reactant asymptote; transmission is nonzero
    # only above the higher of the two asymptotes.
    e_lo = max(0.0, v1 - v2)
    e_hi = v1 + e_max_factor * kt
    energies = np.linspace(e_lo, e_hi, n_grid)

    two_pi_a = 2.0 * np.sqrt(alpha1 * energies / v1) / inv_root
    two_pi_b = 2.0 * np.sqrt(alpha1 * (energies / v1 - 1.0) + alpha2) / inv_root
    # Use the overflow-safe form of (cosh(a+b)-cosh(a-b))/(cosh(a+b)+cosh(d)):
    # divide through by exp(a+b).
    ap_b = two_pi_a + two_pi_b
    numer = (
        1.0 + np.exp(-2.0 * ap_b) - np.exp(-2.0 * two_pi_b) - np.exp(-2.0 * two_pi_a)
    )
    denom = (
        1.0 + np.exp(-2.0 * ap_b) + np.exp(two_pi_d - ap_b) + np.exp(-two_pi_d - ap_b)
    )
    transmission = numer / denom
    boltz = np.exp(-(energies - v1) / kt)  # normalized to the classical TS flux
    kappa = np.trapezoid(transmission * boltz, energies) / kt
    return float(kappa)


@dataclass(frozen=True)
class Thermo:
    """Quasi-RRHO thermochemical summary of one stationary point.

    All energies kJ/mol at the given temperature; entropy included in
    ``gibbs``. ``zpe`` and the thermal terms come from the projected
    harmonic frequencies (imaginary modes must be removed by the caller
    — for a TS, drop the reaction mode and pass it to the tunneling
    correction instead).
    """

    electronic_kj: float
    zpe_kj: float
    thermal_kj: float  # E_trans + E_rot + E_vib(thermal) + RT (H = U + RT)
    entropy_kj_per_k: float  # total S, kJ/(mol K)
    temperature: float

    @property
    def enthalpy(self) -> float:
        return self.electronic_kj + self.zpe_kj + self.thermal_kj

    @property
    def gibbs(self) -> float:
        return self.enthalpy - self.temperature * self.entropy_kj_per_k


def _vib_thermal_and_entropy(
    frequencies_cm: np.ndarray, temperature: float, qrrho_cutoff_cm: float
) -> tuple[float, float]:
    """RRHO vibrational thermal energy (kJ/mol) and quasi-RRHO entropy.

    Entropy uses Grimme's damped rotor interpolation (Chem. Eur. J. 2012):
    below the cutoff region each mode's harmonic entropy is blended with a
    free-rotor entropy via the w(nu) = 1/(1+(nu0/nu)^4) switching function.
    Thermal energy stays harmonic (standard practice; the entropy is the
    quantity the low modes corrupt).
    """
    nu = np.asarray(frequencies_cm, dtype=float)
    if np.any(nu <= 0):
        raise ValueError("frequencies must be positive (project out imaginary modes)")
    theta = wavenumber_to_joule(nu) / KB  # vibrational temperatures, K
    x = theta / temperature
    # Harmonic thermal energy per mode (beyond ZPE): R*theta/(exp(x)-1)
    e_vib = R_KJ * np.sum(theta / np.expm1(x))
    # Harmonic entropy per mode.
    s_ho = R * (x / np.expm1(x) - np.log1p(-np.exp(-x)))  # J/(mol K)
    # Free-rotor entropy per mode with averaged moment Bav.
    mu = H / (8.0 * math.pi**2 * C_CM * nu)  # kg m^2
    mu_p = mu * _BAV / (mu + _BAV)
    s_fr = R * (0.5 + np.log(np.sqrt(8.0 * math.pi**3 * mu_p * KB * temperature) / H))
    w = 1.0 / (1.0 + (qrrho_cutoff_cm / nu) ** 4)
    s_vib = np.sum(w * s_ho + (1.0 - w) * s_fr) / 1000.0  # kJ/(mol K)
    return float(e_vib), float(s_vib)


def thermo_from_frequencies(
    electronic_kj: float,
    frequencies_cm: list[float] | np.ndarray,
    temperature: float,
    *,
    molar_mass_kg: float,
    rotational_temperatures_k: list[float] | None = None,
    symmetry_number: int = 1,
    linear: bool = False,
    pressure_pa: float = 1.0e5,
    qrrho_cutoff_cm: float = 100.0,
) -> Thermo:
    """Quasi-RRHO thermochemistry from projected harmonic frequencies.

    ``rotational_temperatures_k`` are the rotational constants expressed
    as temperatures (h c B / kB); ``None`` means an atom (no rotation).
    Standard state is the ideal gas at ``pressure_pa``; solution-phase
    standard-state corrections are applied downstream where the reaction
    stoichiometry is known.
    """
    t = temperature
    zpe = 0.0
    e_vib = 0.0
    s_vib = 0.0
    freqs = np.asarray(frequencies_cm, dtype=float)
    if freqs.size:
        zpe = 0.5 * NA * float(np.sum(wavenumber_to_joule(freqs))) / 1000.0
        e_vib, s_vib = _vib_thermal_and_entropy(freqs, t, qrrho_cutoff_cm)

    # Translation (Sackur-Tetrode) at the chosen standard state.
    q_trans = (
        (2.0 * math.pi * molar_mass_kg / NA * KB * t) ** 1.5
        / H**3
        * KB
        * t
        / pressure_pa
    )
    s_trans = R * (math.log(q_trans) + 2.5) / 1000.0
    e_trans = 1.5 * R_KJ * t

    # Rotation.
    if rotational_temperatures_k is None:
        s_rot, e_rot = 0.0, 0.0
    elif linear:
        (theta_r,) = rotational_temperatures_k
        q_rot = t / (symmetry_number * theta_r)
        s_rot = R * (math.log(q_rot) + 1.0) / 1000.0
        e_rot = R_KJ * t
    else:
        ta, tb, tc = rotational_temperatures_k
        q_rot = math.sqrt(math.pi * t**3 / (ta * tb * tc)) / symmetry_number
        s_rot = R * (math.log(q_rot) + 1.5) / 1000.0
        e_rot = 1.5 * R_KJ * t

    thermal = e_trans + e_rot + e_vib + R_KJ * t  # H = U + RT (ideal gas)
    entropy = s_trans + s_rot + s_vib
    return Thermo(
        electronic_kj=electronic_kj,
        zpe_kj=zpe,
        thermal_kj=thermal,
        entropy_kj_per_k=entropy,
        temperature=t,
    )


@dataclass(frozen=True)
class RateResult:
    """k(T) for one elementary step, with the pieces kept visible."""

    dg_kj: float
    dh_kj: float
    ds_kj_per_k: float
    kappa: float
    temperature: float

    @property
    def k(self) -> float:
        return self.kappa * eyring_k(self.dg_kj, self.temperature)


def rate_from_thermo(
    reactant: Thermo,
    ts: Thermo,
    *,
    imag_nu_cm: float | None = None,
    tunneling: str = "none",
    dv_reverse_kj: float | None = None,
) -> RateResult:
    """Eyring rate from reactant/TS thermochemistry (+ optional tunneling).

    ``tunneling``: "none", "wigner", or "eckart" (eckart additionally
    needs ``dv_reverse_kj``, the ZPE-corrected reverse barrier).
    """
    if reactant.temperature != ts.temperature:
        raise ValueError("reactant and TS thermochemistry at different temperatures")
    t = ts.temperature
    dh = ts.enthalpy - reactant.enthalpy
    ds = ts.entropy_kj_per_k - reactant.entropy_kj_per_k
    dg = dh - t * ds
    if tunneling == "none":
        kappa = 1.0
    elif tunneling == "wigner":
        if imag_nu_cm is None:
            raise ValueError("wigner correction needs the imaginary frequency")
        kappa = wigner_kappa(imag_nu_cm, t)
    elif tunneling == "eckart":
        if imag_nu_cm is None or dv_reverse_kj is None:
            raise ValueError("eckart correction needs imag_nu_cm and dv_reverse_kj")
        dv_forward = (ts.electronic_kj + ts.zpe_kj) - (
            reactant.electronic_kj + reactant.zpe_kj
        )
        kappa = eckart_kappa(imag_nu_cm, dv_forward, dv_reverse_kj, t)
    else:
        raise ValueError(f"unknown tunneling scheme '{tunneling}'")
    return RateResult(dg_kj=dg, dh_kj=dh, ds_kj_per_k=ds, kappa=kappa, temperature=t)


def arrhenius_fit(
    temperatures: list[float] | np.ndarray, ks: list[float] | np.ndarray
) -> tuple[float, float]:
    """Least-squares (A, Ea_kJ) over ln k = ln A - Ea/(R T).

    The bridge from computed k(T) points to a petra `arrhenius` rate
    spec when an explicit Eyring entry is not wanted.
    """
    t = np.asarray(temperatures, dtype=float)
    k = np.asarray(ks, dtype=float)
    if t.size < 2:
        raise ValueError("need at least two (T, k) points")
    if not (np.all(np.isfinite(k)) and np.all(k > 0.0)):
        raise ValueError("rate constants must be finite and positive")
    x = 1.0 / (R_KJ * t)
    slope, intercept = np.polyfit(x, np.log(k), 1)
    return float(np.exp(intercept)), float(-slope)
