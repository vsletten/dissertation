#!/usr/bin/env python
"""Phase 0 smoke test: environment sanity + CPU-vs-GPU timing table.

Run on the workstation (SURVEY.md §9 Phase 0):

    uv run python scripts/phase0_smoke.py            # full table
    uv run python scripts/phase0_smoke.py --quick    # seconds-long sanity
    uv run python scripts/phase0_smoke.py --waters 2 # bigger cluster

Prints a paste-back-able report: environment, library versions, GPU
detection (gpu4pyscf + cuEST), then B3LYP/def2-TZVP energy, gradient,
and analytic Hessian timings on Si(OH)4·nH2O for CPU and (if present)
GPU. The GPU numbers only mean something on the 4090 box — a cloud
session container reports what it is and says so.
"""

from __future__ import annotations

import argparse
import os
import platform
import sys
import time
from dataclasses import dataclass

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from quarry.clusters import silicic_acid_hydrate, water  # noqa: E402
from quarry.pipeline import DftSettings, _make_scf, build_mol  # noqa: E402


def _try_version(module: str) -> str:
    try:
        mod = __import__(module)
        return getattr(mod, "__version__", "present")
    except ImportError as exc:
        return f"MISSING ({exc.msg})"


def detect_gpu() -> tuple[bool, list[str]]:
    notes = []
    try:
        import cupy

        n = cupy.cuda.runtime.getDeviceCount()
        for i in range(n):
            props = cupy.cuda.runtime.getDeviceProperties(i)
            notes.append(
                f"cuda:{i} {props['name'].decode()} "
                f"({props['totalGlobalMem'] / 2**30:.0f} GiB)"
            )
    except Exception as exc:  # noqa: BLE001 - report, don't crash
        notes.append(f"no CUDA devices via cupy ({type(exc).__name__}: {exc})")
        return False, notes
    try:
        import gpu4pyscf  # noqa: F401

        notes.append(f"gpu4pyscf {_try_version('gpu4pyscf')}")
    except ImportError:
        notes.append("gpu4pyscf MISSING (install the [gpu] extra)")
        return False, notes
    try:
        import cuest  # noqa: F401

        notes.append("nvidia-cuest present (INT8 FP64-emulation path available)")
    except ImportError:
        notes.append("nvidia-cuest not installed (FP64-native path only)")
    return len(notes) > 0 and n > 0, notes


@dataclass
class Timing:
    label: str
    energy_s: float | None = None
    gradient_s: float | None = None
    hessian_s: float | None = None
    e_tot: float | None = None


def time_engine(cluster, settings: DftSettings, *, hessian: bool) -> Timing:
    label = "GPU" if settings.use_gpu else "CPU"
    t = Timing(label=label)
    mf = _make_scf(build_mol(cluster, settings), settings)
    t0 = time.perf_counter()
    t.e_tot = float(mf.kernel())
    t.energy_s = time.perf_counter() - t0
    if not mf.converged:
        raise RuntimeError(f"{label} SCF did not converge")
    t0 = time.perf_counter()
    mf.nuc_grad_method().kernel()
    t.gradient_s = time.perf_counter() - t0
    if hessian:
        t0 = time.perf_counter()
        mf.Hessian().kernel()
        t.hessian_s = time.perf_counter() - t0
    return t


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--quick",
        action="store_true",
        help="H2O at HF/STO-3G, no Hessian — env sanity only",
    )
    ap.add_argument(
        "--waters",
        type=int,
        default=2,
        help="n in Si(OH)4·nH2O for the timing run (default 2)",
    )
    ap.add_argument("--skip-hessian", action="store_true")
    args = ap.parse_args()

    print("=" * 72)
    print("quarry phase-0 smoke test")
    print("=" * 72)
    print(f"host       : {platform.node()} ({platform.machine()})")
    print(f"platform   : {platform.platform()}")
    print(f"python     : {sys.version.split()[0]}")
    print(f"cpu threads: {os.cpu_count()}")
    for lib in ("numpy", "pyscf", "geometric", "sella", "ase"):
        print(f"{lib:<11}: {_try_version(lib)}")
    gpu_ok, notes = detect_gpu()
    for n in notes:
        print(f"gpu        : {n}")
    if not gpu_ok:
        print("gpu        : >>> CPU-only run — GPU timings need the 4090 box <<<")

    if args.quick:
        cluster = water()
        settings = DftSettings(xc="hf", basis="sto-3g")
        hessian = False
        print(f"\nquick mode: {cluster.name} at HF/STO-3G")
    else:
        cluster = silicic_acid_hydrate(args.waters)
        settings = DftSettings(xc="b3lyp", basis="def2-tzvp")
        hessian = not args.skip_hessian
        print(
            f"\ntarget: {cluster.name} ({len(cluster.symbols)} atoms), "
            f"B3LYP/def2-TZVP, analytic Hessian={'yes' if hessian else 'no'}"
        )

    rows = [time_engine(cluster, settings, hessian=hessian)]
    if gpu_ok:
        from dataclasses import replace

        rows.append(
            time_engine(cluster, replace(settings, use_gpu=True), hessian=hessian)
        )

    print()
    print(
        f"{'engine':<8}{'E (Hartree)':>16}{'energy s':>12}"
        f"{'gradient s':>12}{'hessian s':>12}"
    )
    print("-" * 60)
    for r in rows:
        hess = f"{r.hessian_s:12.1f}" if r.hessian_s is not None else f"{'—':>12}"
        print(
            f"{r.label:<8}{r.e_tot:>16.8f}{r.energy_s:>12.1f}"
            f"{r.gradient_s:>12.1f}{hess}"
        )
    if len(rows) == 2:
        for attr in ("energy_s", "gradient_s", "hessian_s"):
            c, g = getattr(rows[0], attr), getattr(rows[1], attr)
            if c and g:
                print(f"speedup {attr.removesuffix('_s'):<9}: {c / g:5.2f}x")
    print("\npaste this whole output back into the session. done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
