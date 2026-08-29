#!/usr/bin/env python
"""D2b: surface-valid Langmuir-Hinshelwood rates for grain-surface channels.

Resolves the explicit NO-GO from D2a (docs/program/results/D2a-rate-reproduction.md):
gas-phase supermolecule rates compared 4-5 orders low against the Song & Kaestner
explicit-ASW instanton fits.  The protocol decomposes and closes that gap in
three auditable rungs, each with its own receipts:

  A. **Surface partition functions** — the LH convention of the field
     (Meisner, Lamberts & Kaestner 2017; Song & Kaestner 2017): unimolecular
     Eyring rates from the pre-reactive adsorbed complex with vibration-only
     partition functions (``quarry.rates.surface_thermo_from_frequencies``).
     Measured effect on D2a artifacts: factor ~2 — NOT the dominant term.
  B. **Coupled-cluster barrier correction** — UHF-UCCSD(T) single points on the
     DFT stationary points; CBS(cc-pVTZ, cc-pVQZ) for the gas-phase systems,
     additive gas-phase Delta-Delta transfer for water clusters (direct TZ on
     1-water clusters as the transferability receipt).  In deep tunneling,
     +-3 kJ/mol of barrier is orders of magnitude of rate; Simons, Lamberts &
     Cuppen 2020 (A&A 634, A52, sect. 3.2) document exactly this correction
     against the as-published Song & Kaestner DFT-based rates.
  C. **Explicit water clusters + site sampling** — each channel on (H2O)n,
     n = 1, 2, in distinct binding orientations; substrate atom order is
     invariant so scan indices carry across cluster sizes.

Literature anchors (the "documented spread" of the card's gate):
  * Song & Kaestner 2017 (ApJ 850, 118; arXiv:2009.05442): explicit-ASW QM/MM
    instanton LH fits, valid >= 59-75 K — the UPPER anchor (documented by
    Simons+20 as overestimated, most strongly for abstraction).
  * Simons, Lamberts & Cuppen 2020 (A&A 634, A52, Table 2): the plateau rates
    a benchmarked KMC actually adopts (2e5 / 1.5e5 s^-1) — the LOWER anchor.
  * Fuchs et al. 2009 (A&A 505, 629, Table 2): experimental *effective*
    Arrhenius barriers at 12.0-16.5 K with nu = kT/h ~ 2e11 s^-1 — compared in
    their own convention via E_eff(T) = -T ln(k / 2e11).
  * Andersson, Goumans & Arnaldsson 2011 (CPL 513, 31): H+CO crossover
    temperature 79 K — compared against our T_c = h nu_imag c / (2 pi kB).

Run on the workstation only::

    OMP_NUM_THREADS=16 uv run python scripts/surface_rate_protocol.py --gpu

Every expensive stage is checkpointed under
``qm/runs/D2b-explicit-surface-rates/<reaction>/``; reruns resume.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import time
from dataclasses import dataclass, replace
from pathlib import Path

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

if __name__ == "__main__":
    from quarry.etiquette import bootstrap_cli

    bootstrap_cli(
        "surface_rate_protocol",
        default_run_root=Path(__file__).resolve().parent.parent / "runs",
        gpu_owner="surface_rate_protocol",
    )

import numpy as np  # noqa: E402

from quarry.clusters import Cluster, _water_centered, _water_donating_to  # noqa: E402
from quarry.pipeline import (  # noqa: E402
    HARTREE_TO_KJ,
    DftSettings,
    FrequencyResult,
    frequencies,
    optimize,
)
from quarry.rates import (  # noqa: E402
    KB,
    R_KJ,
    H,
    eckart_kappa,
    surface_thermo_from_frequencies,
    thermo_from_frequencies,
)
from quarry.ts import (  # noqa: E402
    ScanNoMaximumError,
    find_ts,
    quick_irc,
    scan_to_maximum,
    scan_ts_guess,
)
from scripts.cc_calibration import (  # noqa: E402
    extrapolate_correlation,
    extrapolate_hf,
    frozen_core_orbitals,
)

TEMPERATURES = (50.0, 60.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0)
# Deep-tunneling grid: below every Eckart validity floor in this campaign.
# Reported separately, flagged, ONLY for the Fuchs/Simons low-T comparisons.
DEEP_TEMPERATURES = (12.0, 13.5, 15.0, 16.5, 20.0, 30.0, 40.0)
FUCHS_PREFACTOR_S = 2.0e11  # nu ~ kT/h convention of Fuchs et al. 2009
# Minima may carry numerically imaginary soft librations on clusters; a TS
# must have exactly one imaginary mode above the noise floor and that mode
# must be chemical (>= 200 cm^-1, the D2a gate).  Wet sites get the standard
# cluster-model spectator tolerance (observed live: a healthy 678i saddle on
# h-co-1w-cside carrying a 52.8i water libration); gas stays strict.
MINIMUM_NOISE_FLOOR_GAS_CM = 30.0
MINIMUM_NOISE_FLOOR_WET_CM = 100.0
TS_NOISE_FLOOR_GAS_CM = 50.0
TS_NOISE_FLOOR_WET_CM = 100.0
TS_CHEMICAL_MODE_CM = 200.0
RUN_ROOT = (
    Path(__file__).resolve().parent.parent / "runs" / "D2b-explicit-surface-rates"
)

SONG = "Song & Kaestner 2017, ApJ 850 118 (arXiv:2009.05442), Table 2/Table 4"
SIMONS = "Simons, Lamberts & Cuppen 2020, A&A 634 A52, Table 2"
FUCHS = "Fuchs et al. 2009, A&A 505 629, Table 2 (effective, nu=2e11 s^-1)"
ANDERSSON = "Andersson, Goumans & Arnaldsson 2011, CPL 513 31"


@dataclass(frozen=True)
class LiteratureFit:
    """LH fit k(T) = alpha (T/300)^beta exp(-gamma (T+T0) / (T^2+T0^2))."""

    alpha_s: float
    beta: float
    gamma_k: float
    t0_k: float
    source: str
    valid_floor_k: float

    def rate(self, temperature_k: float) -> float:
        t = temperature_k
        return (
            self.alpha_s
            * (t / 300.0) ** self.beta
            * np.exp(-self.gamma_k * (t + self.t0_k) / (t * t + self.t0_k * self.t0_k))
        )


@dataclass(frozen=True)
class FuchsEffective:
    """One experimental effective (barrier_K, rate_s) row at one temperature."""

    temperature_k: float
    barrier_k: float
    barrier_err_k: float
    rate_s: float


# Fuchs et al. 2009 Table 2.  These are *model-effective* values (classical
# Arrhenius with nu = 2e11 s^-1 inside a lattice-gas KMC that competes reaction
# against hopping); they are NOT intrinsic unimolecular rates and are compared
# only in their own convention.
FUCHS_H_CO = (
    FuchsEffective(12.0, 390.0, 40.0, 2e-3),
    FuchsEffective(13.5, 435.0, 50.0, 2e-3),
    FuchsEffective(15.0, 480.0, 60.0, 5e-3),
    FuchsEffective(16.5, 520.0, 70.0, 4e-3),
)
FUCHS_H_H2CO = (
    FuchsEffective(12.0, 415.0, 40.0, 2e-4),
    FuchsEffective(13.5, 435.0, 50.0, 1e-3),
    FuchsEffective(15.0, 470.0, 60.0, 1e-3),
    FuchsEffective(16.5, 500.0, 70.0, 1e-3),
)


@dataclass(frozen=True)
class Reaction:
    key: str
    label: str
    family: str  # gas-phase family key: 'h-co', 'h-h2co-ch3o', 'h-h2co-h2-hco'
    cluster: Cluster
    n_water: int
    site: str  # orientation/site label for the sampling receipt
    scan_i: int
    scan_j: int
    scan_distances_a: tuple[float, ...]
    scan_floor_a: float
    method: DftSettings
    literature_fit: LiteratureFit | None
    simons_adopted_s: float | None
    fuchs: tuple[FuchsEffective, ...]
    andersson_crossover_k: float | None = None


def log(message: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {message}", flush=True)


def _carbon_monoxide(incoming: np.ndarray, *, key: str) -> Cluster:
    # C, O, incoming H.  The H seeds a bent attack on carbon (HCO ~ 124 deg).
    return Cluster(
        name=key,
        symbols=["C", "O", "H"],
        coords=np.array([[0.0, 0.0, 0.0], [1.128, 0.0, 0.0], incoming]),
        spin=1,
    )


def _formaldehyde(incoming: np.ndarray, *, key: str) -> Cluster:
    # Planar H2CO seed identical to D2a: C, O, Ha, Hb, incoming H.
    return Cluster(
        name=key,
        symbols=["C", "O", "H", "H", "H"],
        coords=np.array(
            [
                [0.0, 0.0, 0.0],
                [1.205, 0.0, 0.0],
                [-0.575, 0.935, 0.0],
                [-0.575, -0.935, 0.0],
                incoming,
            ]
        ),
        spin=1,
    )


def _with_waters(base: Cluster, waters: list[list[np.ndarray]], *, key: str) -> Cluster:
    """Append O/H/H water blocks after the substrate so scan indices survive."""
    symbols = list(base.symbols)
    coords = list(base.coords)
    for block in waters:
        symbols += ["O", "H", "H"]
        coords += list(block)
    return replace(base, name=key, symbols=symbols, coords=np.array(coords))


def _co_site_waters(site: str) -> list[list[np.ndarray]]:
    """Water blocks for one CO binding orientation (site sampling)."""
    o_co = np.array([1.128, 0.0, 0.0])
    c_co = np.zeros(3)
    z = np.array([0.0, 0.0, 1.0])
    if site == "gas":
        return []
    # Donor water H-bonded to the CO oxygen, placed off-axis and on the
    # opposite flank from the incoming H so the attack channel stays open.
    w1_o = o_co + 2.85 * np.array([math.cos(-0.5), math.sin(-0.5), 0.0])
    donor = _water_donating_to(w1_o, o_co, z)
    if site == "1w-oside":
        return [donor]
    if site == "1w-cside":
        # Water oxygen presents a lone pair to the carbon from below.
        w_o = c_co + np.array([0.0, -3.0, 0.3])
        return [_water_centered(w_o, np.array([0.0, -1.0, 0.1]))]
    if site == "2w":
        w2_o = w1_o + 2.80 * np.array([0.45, -0.89, 0.0])
        return [donor, _water_donating_to(w2_o, w1_o, z)]
    raise ValueError(f"unknown CO site '{site}'")


def _h2co_site_waters(site: str) -> list[list[np.ndarray]]:
    """Water blocks for one H2CO binding orientation (site sampling)."""
    o_c = np.array([1.205, 0.0, 0.0])
    z = np.array([0.0, 0.0, 1.0])
    if site == "gas":
        return []
    # Classic in-plane H2CO...HOH complex: water donates to a carbonyl lone
    # pair on the Hb (-y) side, keeping Ha (the abstraction target) and the
    # out-of-plane addition channel clear.
    w1_o = o_c + 2.85 * np.array([math.cos(-1.05), math.sin(-1.05), 0.0])
    donor = _water_donating_to(w1_o, o_c, z)
    if site == "1w":
        return [donor]
    if site == "2w":
        w2_o = w1_o + 2.80 * np.array([0.60, -0.80, 0.0])
        return [donor, _water_donating_to(w2_o, w1_o, z)]
    raise ValueError(f"unknown H2CO site '{site}'")


def reactions(*, gpu: bool, basis: str) -> dict[str, Reaction]:
    pwb6k = DftSettings(
        xc="pwb6k", basis=basis, dispersion="d3bj", density_fit=True, use_gpu=gpu
    )
    song_ch3o = LiteratureFit(3146e10, 1.0, 830.0, 119.6, SONG, 59.0)
    song_abstraction = LiteratureFit(4.13e10, 1.0, 1222.0, 147.7, SONG, 59.0)

    # H+CO is an early TS on a very flat approach: the crest near 1.9 A
    # cleared its neighbors by only ~2e-5 Ha on the first campaign pass, so
    # the grid is fine near the top (see interior_global_maximum fallback).
    co_scan = (2.8, 2.5, 2.2, 2.05, 1.95, 1.9, 1.85, 1.8, 1.7, 1.6, 1.5, 1.4)
    ch3o_scan = (2.8, 2.5, 2.2, 2.0, 1.85, 1.7, 1.55, 1.4, 1.25)
    abstraction_scan = (2.6, 2.3, 2.0, 1.75, 1.5, 1.3, 1.15, 1.02, 0.9)

    result: dict[str, Reaction] = {}

    co_attack = np.array([-1.5, 2.598, 0.0])  # 3.0 A from C, 120 deg off C-O
    for site in ("gas", "1w-oside", "1w-cside", "2w"):
        n_water = {"gas": 0, "1w-oside": 1, "1w-cside": 1, "2w": 2}[site]
        key = "h-co" if site == "gas" else f"h-co-{site}"
        result[key] = Reaction(
            key=key,
            label=f"H + CO -> HCO [{site}]",
            family="h-co",
            cluster=_with_waters(
                _carbon_monoxide(co_attack, key=key), _co_site_waters(site), key=key
            ),
            n_water=n_water,
            site=site,
            scan_i=0,
            scan_j=2,
            scan_distances_a=co_scan,
            scan_floor_a=1.15,
            method=pwb6k,
            literature_fit=None,
            simons_adopted_s=2.0e5,
            fuchs=FUCHS_H_CO,
            andersson_crossover_k=79.0,
        )

    for site in ("gas", "1w", "2w"):
        n_water = {"gas": 0, "1w": 1, "2w": 2}[site]
        key = "h-h2co-ch3o" if site == "gas" else f"h-h2co-ch3o-{site}"
        result[key] = Reaction(
            key=key,
            label=f"H + H2CO -> CH3O [{site}]",
            family="h-h2co-ch3o",
            cluster=_with_waters(
                _formaldehyde(np.array([0.0, 0.0, 3.0]), key=key),
                _h2co_site_waters(site),
                key=key,
            ),
            n_water=n_water,
            site=site,
            scan_i=0,
            scan_j=4,
            scan_distances_a=ch3o_scan,
            scan_floor_a=1.1,
            method=pwb6k,
            literature_fit=song_ch3o,
            simons_adopted_s=1.5e5,
            fuchs=FUCHS_H_H2CO,
        )
        key = "h-h2co-h2-hco" if site == "gas" else f"h-h2co-h2-hco-{site}"
        result[key] = Reaction(
            key=key,
            label=f"H + H2CO -> H2 + HCO [{site}]",
            family="h-h2co-h2-hco",
            cluster=_with_waters(
                _formaldehyde(np.array([-0.575, 0.935, 2.8]), key=key),
                _h2co_site_waters(site),
                key=key,
            ),
            n_water=n_water,
            site=site,
            scan_i=2,
            scan_j=4,
            scan_distances_a=abstraction_scan,
            scan_floor_a=0.78,
            method=pwb6k,
            literature_fit=song_abstraction,
            simons_adopted_s=1.5e5,
            fuchs=FUCHS_H_H2CO,
        )
    return result


def save_xyz(cluster: Cluster, path: Path) -> None:
    path.write_text(cluster.to_xyz())


def load_xyz(path: Path, template: Cluster) -> Cluster:
    lines = path.read_text().splitlines()
    n = int(lines[0])
    coords = np.array(
        [[float(v) for v in line.split()[1:4]] for line in lines[2 : 2 + n]]
    )
    return replace(template, coords=coords)


def geometry_hash(cluster: Cluster) -> str:
    # A fixed comment keeps the hash name-insensitive: a checkpoint-resumed
    # cluster inherits its load template's name, and CC receipt identity must
    # survive that (observed live: cc-reactant-tz.json drift on resume).
    return hashlib.sha256(cluster.to_xyz(comment="geometry").encode()).hexdigest()


def interior_global_maximum(
    scan: list[tuple[float, float, Cluster]],
    *,
    min_rise_hartree: float = 1.0e-4,
) -> Cluster:
    """Shallow-crest fallback TS guess: the interior global energy maximum.

    An early transition state on a flat approach (H + CO: the crest near
    1.9 A clears its neighbors by ~2e-5 Ha) fails ``scan_to_maximum``'s
    strict per-neighbor tolerance while still being a genuine interior
    global maximum well above both scan endpoints.  Sella plus the
    chemical-imaginary-mode gate remain the actual saddle acceptance;
    this only hands over a physically sensible seed.
    """
    if len(scan) < 3:
        raise ValueError("scan too short for an interior maximum")
    energies = [e for _, e, _ in scan]
    peak = int(np.argmax(energies))
    if peak in (0, len(scan) - 1):
        raise ValueError("scan global maximum is an endpoint — ridge not crossed")
    if (
        energies[peak] - energies[0] < min_rise_hartree
        or energies[peak] - energies[-1] < min_rise_hartree
    ):
        raise ValueError(
            "scan global maximum does not clear the endpoints by "
            f"{min_rise_hartree:.1e} Ha — noise, not a crest"
        )
    return scan[peak][2]


def crossover_temperature_k(imag_cm: float) -> float:
    """T_c = h c nu / (2 pi kB): where deep tunneling takes over."""
    # h*c*100/kB = 1.438777 cm K
    return 1.4387769 * abs(imag_cm) / (2.0 * math.pi)


def fuchs_effective_barrier_k(k_s: float, temperature_k: float) -> float:
    """Our k(T) restated in the Fuchs et al. 2009 effective-barrier convention."""
    if k_s <= 0.0:
        raise ValueError("rate must be positive")
    return -temperature_k * math.log(k_s / FUCHS_PREFACTOR_S)


def minimum_noise_floor_cm(n_water: int) -> float:
    # Physisorbed clusters keep their soft water librations numerically
    # imaginary below ~100 cm^-1 at production convergence (observed live:
    # 76.5i on the h-co-2w pre-reactive complex).  Same spectator tolerance
    # as the wet TS gate, disclosed via the stored imaginary-mode lists.
    return MINIMUM_NOISE_FLOOR_GAS_CM if n_water == 0 else MINIMUM_NOISE_FLOOR_WET_CM


def require_minimum(
    freq: FrequencyResult,
    *,
    name: str,
    noise_floor_cm: float = MINIMUM_NOISE_FLOOR_GAS_CM,
) -> None:
    significant = freq.imaginary_cm[freq.imaginary_cm > noise_floor_cm]
    if significant.size:
        raise RuntimeError(
            f"{name}: endpoint is not a minimum — imaginary modes "
            f"{np.round(freq.imaginary_cm, 1).tolist()} cm^-1"
        )


def ts_noise_floor_cm(n_water: int) -> float:
    return TS_NOISE_FLOOR_GAS_CM if n_water == 0 else TS_NOISE_FLOOR_WET_CM


def require_chemical_saddle(
    freq: FrequencyResult, *, name: str, noise_floor_cm: float = TS_NOISE_FLOOR_GAS_CM
) -> float:
    significant = freq.imaginary_cm[freq.imaginary_cm > noise_floor_cm]
    if significant.size != 1:
        raise RuntimeError(
            f"{name}: expected exactly 1 imaginary mode above "
            f"{noise_floor_cm:.0f} cm^-1, got "
            f"{np.round(freq.imaginary_cm, 1).tolist()}"
        )
    mode = float(significant[0])
    if mode < TS_CHEMICAL_MODE_CM:
        raise RuntimeError(
            f"{name}: imaginary mode {mode:.1f} cm^-1 is below the chemical "
            f"gate ({TS_CHEMICAL_MODE_CM:.0f} cm^-1) — trivial rearrangement saddle"
        )
    return mode


def frequency_payload(freq: FrequencyResult) -> dict:
    return {
        "frequencies_cm": freq.frequencies_cm.tolist(),
        "imaginary_cm": freq.imaginary_cm.tolist(),
        "electronic_hartree": freq.electronic_hartree,
        "molar_mass_kg": freq.molar_mass_kg,
        "rotational_temperatures_k": list(freq.rotational_temperatures_k or []),
        "linear": freq.linear,
    }


def gas_thermo(freq: FrequencyResult, temperature: float):
    return thermo_from_frequencies(
        freq.electronic_hartree * HARTREE_TO_KJ,
        freq.frequencies_cm,
        temperature,
        molar_mass_kg=freq.molar_mass_kg,
        rotational_temperatures_k=(
            list(freq.rotational_temperatures_k)
            if freq.rotational_temperatures_k
            else None
        ),
        linear=freq.linear,
    )


def surface_thermo(freq: FrequencyResult, temperature: float):
    return surface_thermo_from_frequencies(
        freq.electronic_hartree * HARTREE_TO_KJ,
        freq.frequencies_cm,
        temperature,
    )


# ---------------------------------------------------------------------------
# Open-shell coupled-cluster single points (rung B).
# The A2 cc_calibration engines are closed-shell only; every D2b species is a
# doublet, so this path runs PySCF UHF -> UCCSD(T) in-process with the same
# receipt discipline (atomic JSON, identity-checked cache).


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def uccsdt_receipt(
    cluster: Cluster, basis: str, out_path: Path, *, threads: int = 16
) -> dict:
    """Frozen-core UHF-UCCSD(T) single point with a cached JSON receipt."""
    identity = {
        "method": "uhf-uccsd(t)",
        "basis": basis.lower(),
        "frozen_core_orbitals": frozen_core_orbitals(cluster.symbols),
        "charge": cluster.charge,
        "spin_2s": cluster.spin,
        "geometry_sha256": geometry_hash(cluster),
    }
    if out_path.exists():
        payload = json.loads(out_path.read_text())
        if payload.get("identity") == identity:
            return payload
        raise ValueError(f"{out_path}: existing CC receipt identity drift")

    from pyscf import cc as pyscf_cc
    from pyscf import gto, lib, scf

    lib.num_threads(threads)
    mol = gto.M(
        atom=cluster.to_pyscf_atom(),
        basis=basis,
        charge=cluster.charge,
        spin=cluster.spin,
        verbose=0,
    )
    started = time.time()
    mf = scf.UHF(mol)
    mf.max_cycle = 200
    mf.conv_tol = 1e-10
    scf_energy = float(mf.kernel())
    if not mf.converged:
        raise RuntimeError(f"UHF did not converge for {cluster.name}/{basis}")
    spin_square = float(mf.spin_square()[0])
    mycc = pyscf_cc.UCCSD(mf, frozen=identity["frozen_core_orbitals"] or None)
    mycc.max_cycle = 120
    mycc.conv_tol = 1e-8
    ccsd_correlation = float(mycc.kernel()[0])
    if not mycc.converged:
        raise RuntimeError(f"UCCSD did not converge for {cluster.name}/{basis}")
    triples = float(mycc.ccsd_t())
    payload = {
        "engine": "pyscf-uhf-uccsd(t)",
        "identity": identity,
        "scf_hartree": scf_energy,
        "ccsd_correlation_hartree": ccsd_correlation,
        "triples_hartree": triples,
        "correlation_hartree": ccsd_correlation + triples,
        "total_hartree": scf_energy + ccsd_correlation + triples,
        "uhf_spin_square": spin_square,
        "elapsed_seconds": time.time() - started,
        "completed_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    }
    atomic_json(out_path, payload)
    return payload


def cbs_total_hartree(tz: dict, qz: dict) -> float:
    return extrapolate_hf(
        float(tz["scf_hartree"]), float(qz["scf_hartree"])
    ) + extrapolate_correlation(
        float(tz["correlation_hartree"]), float(qz["correlation_hartree"])
    )


def cc_electronic_kj(
    stationary: dict[str, Cluster],
    run_dir: Path,
    *,
    mode: str,
    threads: int,
) -> dict | None:
    """CC electronic energies (kJ/mol) for reactant/ts/product, or None.

    ``mode``: 'cbs' (TZ+QZ extrapolation, gas-phase systems) or 'tz'
    (single-basis receipt, 1-water clusters).
    """
    if mode == "skip":
        return None
    bases = ("cc-pVTZ", "cc-pVQZ") if mode == "cbs" else ("cc-pVTZ",)
    receipts: dict[str, dict[str, dict]] = {}
    for role, cluster in stationary.items():
        receipts[role] = {}
        for basis in bases:
            tag = basis.lower().replace("cc-pv", "").replace("z", "") + "z"
            out = run_dir / f"cc-{role}-{tag}.json"
            log(f"{run_dir.name}: UCCSD(T)/{basis} single point on {role}")
            receipts[role][basis] = uccsdt_receipt(cluster, basis, out, threads=threads)
    energies = {}
    for role in stationary:
        if mode == "cbs":
            total = cbs_total_hartree(
                receipts[role]["cc-pVTZ"], receipts[role]["cc-pVQZ"]
            )
        else:
            total = float(receipts[role]["cc-pVTZ"]["total_hartree"])
        energies[role] = total * HARTREE_TO_KJ
    return {
        "mode": mode,
        "electronic_kj": energies,
        "uhf_spin_square": {
            role: {basis: receipts[role][basis]["uhf_spin_square"] for basis in bases}
            for role in stationary
        },
    }


# ---------------------------------------------------------------------------


def cc_electronic_barrier_deltas(
    cc_electronic_kj: dict[str, float],
    dft_electronic_kj: dict[str, float],
) -> tuple[float, float]:
    """Forward and reverse CC-minus-DFT electronic barrier corrections (kJ/mol).

    Each well (reactant / product) has its own CC vs DFT shift; reusing the
    TS-minus-reactant delta on the reverse barrier builds the wrong Eckart
    asymmetry whenever those corrections differ.
    """

    def delta(well: str) -> float:
        cc_barrier = cc_electronic_kj["ts"] - cc_electronic_kj[well]
        dft_barrier = dft_electronic_kj["ts"] - dft_electronic_kj[well]
        return cc_barrier - dft_barrier

    return delta("reactant"), delta("product")


def rate_table(
    reactant_freq: FrequencyResult,
    ts_freq: FrequencyResult,
    *,
    imag_cm: float,
    barrier_zpe_dft_kj: float,
    reverse_zpe_dft_kj: float,
    cc_delta_kj: float | None,
    literature_fit: LiteratureFit | None,
    temperatures: tuple[float, ...],
    extrapolated: bool,
    cc_reverse_delta_kj: float | None = None,
) -> list[dict]:
    """Per-temperature decomposition: gas vs surface PFs x DFT vs CC barrier."""
    rows = []
    barrier_cc = barrier_zpe_dft_kj + cc_delta_kj if cc_delta_kj is not None else None
    reverse_delta = (
        cc_reverse_delta_kj if cc_reverse_delta_kj is not None else cc_delta_kj
    )
    reverse_cc = (
        reverse_zpe_dft_kj + reverse_delta if reverse_delta is not None else None
    )
    for t in temperatures:
        classical = KB * t / H
        gas_dg = gas_thermo(ts_freq, t).gibbs - gas_thermo(reactant_freq, t).gibbs
        surf_dg = (
            surface_thermo(ts_freq, t).gibbs - surface_thermo(reactant_freq, t).gibbs
        )
        kappa_dft = eckart_kappa(imag_cm, barrier_zpe_dft_kj, reverse_zpe_dft_kj, t)
        row = {
            "temperature_k": t,
            "eckart_extrapolated": extrapolated,
            "delta_g_gas_kj": gas_dg,
            "delta_g_surface_kj": surf_dg,
            "eckart_kappa_dft": kappa_dft,
            "k_gas_dft_s": kappa_dft * classical * math.exp(-gas_dg / (R_KJ * t)),
            "k_surface_dft_s": kappa_dft * classical * math.exp(-surf_dg / (R_KJ * t)),
        }
        if barrier_cc is not None and reverse_cc is not None:
            if barrier_cc <= 0.0 or reverse_cc <= 0.0:
                raise RuntimeError(
                    "CC-corrected Eckart barrier is non-positive "
                    f"(forward={barrier_cc}, reverse={reverse_cc})"
                )
            kappa_cc = eckart_kappa(imag_cm, barrier_cc, reverse_cc, t)
            surf_dg_cc = surf_dg + (barrier_cc - barrier_zpe_dft_kj)
            row["eckart_kappa_cc"] = kappa_cc
            row["k_surface_cc_s"] = (
                kappa_cc * classical * math.exp(-surf_dg_cc / (R_KJ * t))
            )
        headline = row.get("k_surface_cc_s", row["k_surface_dft_s"])
        row["k_headline_s"] = headline
        row["fuchs_effective_barrier_k"] = fuchs_effective_barrier_k(headline, t)
        if literature_fit is not None:
            lit = literature_fit.rate(t)
            row["literature_fit_s"] = lit
            row["headline_over_fit"] = headline / lit
            row["fit_extrapolated"] = bool(t < literature_fit.valid_floor_k)
        rows.append(row)
    return rows


def provenance(reaction: Reaction, structures: dict[str, Cluster]) -> dict:
    try:
        git_sha = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parent, text=True
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        git_sha = "unknown"
    return {
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "git_sha": git_sha,
        "python": platform.python_version(),
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
        "gpu_requested": reaction.method.use_gpu,
        "site": reaction.site,
        "n_water": reaction.n_water,
        "literature_sources": {
            "fit": reaction.literature_fit.source if reaction.literature_fit else None,
            "simons_adopted": SIMONS if reaction.simons_adopted_s else None,
            "fuchs_effective": FUCHS if reaction.fuchs else None,
            "andersson_crossover": (
                ANDERSSON if reaction.andersson_crossover_k else None
            ),
        },
        "geometry_sha256": {
            name: geometry_hash(cluster) for name, cluster in structures.items()
        },
    }


def run_reaction(
    reaction: Reaction,
    *,
    cc_mode: str,
    gas_cc_delta_kj: dict[str, float],
    threads: int,
    force: bool = False,
) -> dict:
    run_dir = RUN_ROOT / reaction.key
    run_dir.mkdir(parents=True, exist_ok=True)
    result_path = run_dir / "results.json"
    if result_path.exists() and not force:
        log(f"{reaction.key}: results.json exists; resume complete")
        return json.loads(result_path.read_text())

    ts_guess_path = run_dir / "ts_guess.xyz"
    scan_path = run_dir / "scan.json"
    if ts_guess_path.exists() and not force:
        log(f"{reaction.key}: ts_guess.xyz exists; skipping completed scan")
        ts_guess = load_xyz(ts_guess_path, reaction.cluster)
    else:
        seed = reaction.cluster
        if reaction.n_water:
            # Hand-built water networks settle a long way on their first
            # relaxation; doing that inside a constrained scan point strained
            # the SCF to divergence (observed live: h-h2co-h2-hco-2w died 35
            # geomeTRIC steps into its first pinned point).  Relax the seed
            # unconstrained first, checkpointed, then drive the coordinate.
            preopt_path = run_dir / "preopt.xyz"
            if preopt_path.exists():
                seed = load_xyz(preopt_path, reaction.cluster)
            else:
                log(f"{reaction.key}: pre-optimizing wet seed")
                seed = optimize(reaction.cluster, reaction.method, max_steps=200)
                save_xyz(seed, preopt_path)
        log(f"{reaction.key}: relaxed scan")
        try:
            scan = scan_to_maximum(
                seed,
                reaction.method,
                atom_i=reaction.scan_i,
                atom_j=reaction.scan_j,
                distances_a=list(reaction.scan_distances_a),
                min_distance_a=reaction.scan_floor_a,
                progress=lambda r, e: log(f"{reaction.key}: r={r:.3f} A E={e:.10f} Ha"),
            )
        except ScanNoMaximumError as exc:
            scan = exc.scan
        scan_path.write_text(
            json.dumps(
                [{"distance_a": r, "energy_hartree": e} for r, e, _ in scan],
                indent=2,
            )
        )
        try:
            ts_guess = scan_ts_guess(scan)
        except ValueError:
            ts_guess = interior_global_maximum(scan)
            log(
                f"{reaction.key}: no tolerance-cleared crest; seeding Sella "
                "from the interior global maximum at "
                f"r={max(scan, key=lambda p: p[1])[0]:.2f} A"
            )
        save_xyz(ts_guess, ts_guess_path)

    ts_path = run_dir / "ts.xyz"
    ts = (
        load_xyz(ts_path, ts_guess)
        if ts_path.exists()
        else find_ts(ts_guess, reaction.method, trajectory=str(run_dir / "sella.traj"))
    )
    save_xyz(ts, ts_path)
    ts_freq = frequencies(ts, reaction.method)
    noise_floor = ts_noise_floor_cm(reaction.n_water)
    imag = require_chemical_saddle(
        ts_freq, name=reaction.key, noise_floor_cm=noise_floor
    )

    back_path, fwd_path = run_dir / "irc_back.xyz", run_dir / "irc_fwd.xyz"
    if back_path.exists() and fwd_path.exists():
        back, fwd = load_xyz(back_path, ts), load_xyz(fwd_path, ts)
    else:
        back, fwd = quick_irc(ts, reaction.method, noise_floor_cm=noise_floor)
        save_xyz(back, back_path)
        save_xyz(fwd, fwd_path)

    def distance(cluster: Cluster) -> float:
        return float(
            np.linalg.norm(
                cluster.coords[reaction.scan_i] - cluster.coords[reaction.scan_j]
            )
        )

    reactant, product = sorted((back, fwd), key=distance, reverse=True)
    reactant_freq = frequencies(reactant, reaction.method)
    product_freq = frequencies(product, reaction.method)
    minimum_floor = minimum_noise_floor_cm(reaction.n_water)
    require_minimum(
        reactant_freq,
        name=f"{reaction.key} reactant complex",
        noise_floor_cm=minimum_floor,
    )
    require_minimum(
        product_freq, name=f"{reaction.key} product", noise_floor_cm=minimum_floor
    )

    def zpe(freq: FrequencyResult) -> float:
        return surface_thermo(freq, 1.0).zpe_kj

    barrier_zpe = (
        ts_freq.electronic_hartree * HARTREE_TO_KJ
        + zpe(ts_freq)
        - reactant_freq.electronic_hartree * HARTREE_TO_KJ
        - zpe(reactant_freq)
    )
    reverse_zpe = (
        ts_freq.electronic_hartree * HARTREE_TO_KJ
        + zpe(ts_freq)
        - product_freq.electronic_hartree * HARTREE_TO_KJ
        - zpe(product_freq)
    )
    if barrier_zpe <= 0 or reverse_zpe <= 0:
        raise RuntimeError(
            f"{reaction.key}: non-positive Eckart barrier "
            f"(forward={barrier_zpe}, reverse={reverse_zpe})"
        )

    # Rung B: coupled-cluster electronic correction to the forward barrier.
    stationary = {"reactant": reactant, "ts": ts, "product": product}
    reaction_cc_mode = "skip"
    if cc_mode != "skip":
        reaction_cc_mode = (
            "cbs"
            if reaction.n_water == 0
            else ("tz" if reaction.n_water == 1 else "skip")
        )
    cc = cc_electronic_kj(stationary, run_dir, mode=reaction_cc_mode, threads=threads)
    dft_electronic_kj = {
        "reactant": reactant_freq.electronic_hartree * HARTREE_TO_KJ,
        "ts": ts_freq.electronic_hartree * HARTREE_TO_KJ,
        "product": product_freq.electronic_hartree * HARTREE_TO_KJ,
    }
    dft_barrier_electronic = dft_electronic_kj["ts"] - dft_electronic_kj["reactant"]
    cc_delta = None
    cc_reverse_delta = None
    cc_delta_source = None
    if cc is not None:
        cc_delta, cc_reverse_delta = cc_electronic_barrier_deltas(
            cc["electronic_kj"], dft_electronic_kj
        )
        cc_delta_source = f"direct-{cc['mode']}"
    elif cc_mode != "skip" and reaction.family in gas_cc_delta_kj:
        # Additive gas-phase Delta-Delta transfer for clusters without direct CC.
        cc_delta = gas_cc_delta_kj[reaction.family]
        cc_reverse_delta = cc_delta
        cc_delta_source = "gas-phase-transfer"
    require_requested_cc_delta(
        cc_mode=cc_mode,
        key=reaction.key,
        n_water=reaction.n_water,
        family=reaction.family,
        cc_delta=cc_delta,
    )

    rates = rate_table(
        reactant_freq,
        ts_freq,
        imag_cm=imag,
        barrier_zpe_dft_kj=barrier_zpe,
        reverse_zpe_dft_kj=reverse_zpe,
        cc_delta_kj=cc_delta,
        cc_reverse_delta_kj=cc_reverse_delta,
        literature_fit=reaction.literature_fit,
        temperatures=TEMPERATURES,
        extrapolated=False,
    )
    deep_rates = rate_table(
        reactant_freq,
        ts_freq,
        imag_cm=imag,
        barrier_zpe_dft_kj=barrier_zpe,
        reverse_zpe_dft_kj=reverse_zpe,
        cc_delta_kj=cc_delta,
        cc_reverse_delta_kj=cc_reverse_delta,
        literature_fit=reaction.literature_fit,
        temperatures=DEEP_TEMPERATURES,
        extrapolated=True,
    )

    result = {
        "reaction": reaction.label,
        "key": reaction.key,
        "family": reaction.family,
        "site": reaction.site,
        "n_water": reaction.n_water,
        "classification": "first-order-saddle",
        "barrier_zpe_dft_kj_mol": barrier_zpe,
        "reverse_barrier_zpe_dft_kj_mol": reverse_zpe,
        "dft_barrier_electronic_kj_mol": dft_barrier_electronic,
        "cc_delta_kj_mol": cc_delta,
        "cc_reverse_delta_kj_mol": cc_reverse_delta,
        "cc_delta_source": cc_delta_source,
        "cc": cc,
        "imaginary_frequency_cm": imag,
        "crossover_temperature_k": crossover_temperature_k(imag),
        "andersson_crossover_k": reaction.andersson_crossover_k,
        "simons_adopted_s": reaction.simons_adopted_s,
        "fuchs_effective": [vars(row) for row in reaction.fuchs],
        "rates": rates,
        "deep_tunneling_rates": deep_rates,
        "method": vars(reaction.method),
        "frequency_data": {
            "reactant": frequency_payload(reactant_freq),
            "ts": frequency_payload(ts_freq),
            "product": frequency_payload(product_freq),
        },
        "provenance": provenance(reaction, stationary),
    }
    result_path.write_text(json.dumps(result, indent=2))
    return result


def gas_family_key(all_reactions: dict[str, Reaction], family: str) -> str:
    for key, reaction in all_reactions.items():
        if reaction.family == family and reaction.n_water == 0:
            return key
    raise KeyError(f"no gas-phase reaction for family '{family}'")


def ensure_cc_transfer_prerequisites(
    selected: list[str],
    all_reactions: dict[str, Reaction],
    *,
    cc_mode: str,
) -> None:
    """Fail fast when --cc full cannot produce a 2-water transfer correction."""
    if cc_mode == "skip":
        return
    selected_set = set(selected)
    missing: list[str] = []
    for key in selected:
        reaction = all_reactions[key]
        if reaction.n_water < 2:
            continue
        gas_key = gas_family_key(all_reactions, reaction.family)
        if gas_key not in selected_set:
            missing.append(f"{key} needs {gas_key}")
    if missing:
        details = "; ".join(missing)
        raise SystemExit(
            f"--cc {cc_mode} requires the gas-phase family reaction for each "
            f"2-water selection so the additive CC transfer exists. Missing: "
            f"{details}. Add those --reaction keys or pass --cc skip."
        )


def require_requested_cc_delta(
    *,
    cc_mode: str,
    key: str,
    n_water: int,
    family: str,
    cc_delta: float | None,
) -> None:
    """Refuse a successful DFT-only result when the caller asked for rung B."""
    if cc_mode == "skip" or cc_delta is not None:
        return
    raise RuntimeError(
        f"{key}: --cc {cc_mode} requested but no coupled-cluster correction "
        f"is available (n_water={n_water}, family={family}). Include the "
        f"gas-phase family reaction in the same invocation, or pass --cc skip."
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--threads", type=int, default=16, help="numerical thread cap")
    parser.add_argument("--nice", type=int, default=10, help="niceness increment")
    parser.add_argument("--log", help="tee output to this path (default: qm/runs)")
    parser.add_argument(
        "--reaction", action="append", help="reaction key (repeatable; default all)"
    )
    parser.add_argument("--gpu", action="store_true", help="use GPU4PySCF")
    parser.add_argument("--gpu-mem-gb", type=float, help="CuPy pool cap (etiquette)")
    parser.add_argument("--basis", default="def2-svp", help="DFT basis")
    parser.add_argument(
        "--cc",
        choices=("full", "skip"),
        default="full",
        help="rung-B coupled-cluster single points",
    )
    parser.add_argument(
        "--force", action="store_true", help="recompute completed results"
    )
    args = parser.parse_args()
    all_reactions = reactions(gpu=args.gpu, basis=args.basis)
    selected = args.reaction or list(all_reactions)
    unknown = sorted(set(selected) - set(all_reactions))
    if unknown:
        raise SystemExit(f"unknown reactions: {unknown}")
    ensure_cc_transfer_prerequisites(selected, all_reactions, cc_mode=args.cc)
    RUN_ROOT.mkdir(parents=True, exist_ok=True)

    # Gas-phase reactions run first so cluster reactions can inherit the
    # additive Delta-Delta correction of their family.
    selected = sorted(selected, key=lambda key: all_reactions[key].n_water)
    results, failures = [], []
    gas_cc_delta_kj: dict[str, float] = {}
    for key in selected:
        reaction = all_reactions[key]
        try:
            result = run_reaction(
                reaction,
                cc_mode=args.cc,
                gas_cc_delta_kj=gas_cc_delta_kj,
                threads=args.threads,
                force=args.force,
            )
            results.append(result)
            if (
                reaction.n_water == 0
                and result.get("cc_delta_kj_mol") is not None
                and result.get("cc_delta_source", "").startswith("direct")
            ):
                gas_cc_delta_kj[reaction.family] = result["cc_delta_kj_mol"]
        except Exception as exc:  # campaign isolation: preserve and continue
            failure = {
                "key": key,
                "error_type": type(exc).__name__,
                "error": str(exc),
                "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            }
            failures.append(failure)
            (RUN_ROOT / key).mkdir(parents=True, exist_ok=True)
            (RUN_ROOT / key / "failure.json").write_text(json.dumps(failure, indent=2))
            log(f"{key}: FAILED but campaign continues: {exc}")
    (RUN_ROOT / "results.json").write_text(json.dumps(results, indent=2))
    (RUN_ROOT / "failures.json").write_text(json.dumps(failures, indent=2))
    log(
        f"complete: {len(results)} results, {len(failures)} failures -> "
        f"{RUN_ROOT / 'results.json'}"
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
