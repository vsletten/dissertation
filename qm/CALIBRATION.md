# The calibration layer — open-source coupled cluster (no ORCA)

*Decision 2026-08-23 (Victor): ORCA is descoped from the platform on
licensing grounds (registration-gated, closed, no redistribution — the
one component that could never be containerized or committed). Its only
role was DLPNO-CCSD(T) calibration single points (SURVEY §6.4). Two
open engines replace it, and cross-validate each other:*

| Engine | Method | License | Role |
|---|---|---|---|
| **Psi4 1.11** | DLPNO-CCSD(T) (native since 1.11, Jiang et al. JCP 161, 082502) | LGPL-3.0 | The like-for-like DLPNO replacement; CPU |
| **ByteQC** | **Canonical** CCSD(T) on GPU (cucc; ~1380 orbitals on 80 GB) | Apache-2.0 | The trust anchor — no local approximation at all |

Running both on the same geometry gives the calibration number *and* an
explicit local-approximation error bar (DLPNO − canonical), something the
single-engine ORCA plan never provided. Psi4's SURVEY §2.4 "dealbreakers"
(finite-difference Hessians, no SMD) are irrelevant here: calibration is
single points on frozen geometries; frequencies and solvation stay in the
PySCF/GPU4PySCF pipeline.

## Installed environments (workstation, 2026-08-23)

Both live outside uv (like ORCA would have) but unlike ORCA are fully
rebuildable from public sources.

### Psi4 — conda env `calib-psi4` (verified)

```bash
conda create -n calib-psi4 -c conda-forge psi4 -y
conda activate calib-psi4
```

Verified install: psi4 1.11. Smoke (water dimer, cc-pVTZ, frozen core,
8 threads): DLPNO-CCSD(T) = −152.67402179 Eh vs canonical CCSD(T) =
−152.67382195 Eh — **agreement −0.53 kJ/mol**. Etiquette: Psi4 is
CPU-only; `psi4.set_num_threads(≤16)` per the OMP cap.

### ByteQC — venv `~/venvs/byteqc` + source at `~/opt/byteqc`

```bash
uv venv ~/venvs/byteqc --python 3.12
source ~/venvs/byteqc/bin/activate
uv pip install cupy-cuda12x nvmath-python h5py numpy scipy "pyscf==2.5.0" \
    cutensor-cu12 nvidia-nccl-cu12 mpi4py cmake ninja
git clone https://github.com/bytedance/byteqc ~/opt/byteqc
# Kernel build — traps on this box, all hit and solved 2026-08-23
# (.claude/learnings/gotchas.md):
#   1. snap cmake fails CUDA compiler detection — use the pip cmake (venv).
#   2. plain `cmake` finds the distro CUDA 12.0 nvcc (/usr/bin/nvcc),
#      which cannot parse gcc-13 glibc headers (`_Float32` errors) —
#      export CUDACXX=/usr/local/cuda-12.6/bin/nvcc.
#   3. a failed configure leaves a poisoned CMakeCache pinning the bad
#      compiler — rm -rf the */lib/build dirs before retrying.
export CUDA_HOME=/usr/local/cuda-12.6 CUDACXX=/usr/local/cuda-12.6/bin/nvcc
rm -rf ~/opt/byteqc/cuobc/lib/build
cmake ~/opt/byteqc/cuobc/lib -B ~/opt/byteqc/cuobc/lib/build -DCMAKE_CUDA_ARCHITECTURES=89
make -j12 -C ~/opt/byteqc/cuobc/lib/build  # capped -j: 32-way nvcc can swap-thrash this box
# (cupbc — periodic-boundary mean field — is NOT needed for cluster
# calibration and its CMake wants a cuTENSOR tarball/conda layout; skip
# it unless PBC work appears. cucc/cump2 are pure Python+cupy, no build.)

# Runtime: cupy's cutensor backend dlopens libcutensor AND libcutensorMg
# (which needs NCCL) at import — both wheel lib dirs must be on the
# loader path (the phase2_ladder.preload_cutensor gotcha, extended):
SP=~/venvs/byteqc/lib/python3.12/site-packages
export LD_LIBRARY_PATH="$SP/cutensor/lib:$SP/nvidia/nccl/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH=~/opt:$PYTHONPATH   # the *parent* of the clone
python -c "from byteqc import cucc"
```

Pinned deps per upstream README: python ≥3.11, `pyscf==2.5.0` (hard pin
— keep this venv separate from quarry's), cupy, cutensor ≥2.1. GPU
support is explicitly tested upstream on consumer Ada (4070 Ti), so the
4090 (sm_89) is in-family.

API note: `(T)` does not hang off the CCSD object — use the module
kernel: `from byteqc.cucc.ccsd_t import kernel; e_t = kernel(mycc,
mycc.ao2mo())`. `cucc.CCSD(mf, frozen=..., gpulim=16 << 30)` caps GPU
memory (QI2 lease etiquette).

**Verified 2026-08-23 (water dimer, cc-pVTZ, frozen core, 4090):**
ByteQC canonical CCSD(T) = −152.67382194 Eh vs Psi4 canonical CCSD(T)
= −152.67382195 Eh — **two independent engines, CPU vs GPU, agree to
1e-8 Eh**. The cross-engine triangle (Psi4-DLPNO / Psi4-canonical /
ByteQC-canonical) is closed: DLPNO −0.53 kJ/mol off canonical, the two
canonicals identical.

## The calibration protocol (si-neutral first)

Target: SURVEY §6.4 — a CCSD(T)/CBS-quality electronic barrier for the
si-neutral reaction (18 atoms: H₂O + (HO)₃Si–O–Si(OH)₃), to rank
B3LYP-D4 / r²SCAN-3c / ωB97M-V for Si–O chemistry. Geometries: the
banked stationary points (complex, TS, product) from
`runs/phase1/si-neutral-b3lyp-def2-svp-flank/results.json` — single
points only, no re-optimization.

1. **Canonical anchor (ByteQC, GPU)**: CCSD(T)/cc-pVTZ and cc-pVQZ on
   all three structures. At ~440 basis functions (TZ) this is far inside
   cucc's demonstrated envelope; the 4090's 24 GB suffices at TZ — if QZ
   (~800 bf) OOMs, it is a few dollars on a rented H200/B200 (NOT a
   B300: Blackwell Ultra's FP64 is cut to ~1.2 TFLOPS, 4090-class;
   rent Hopper/B200 for FP64 work).
2. **Focal-point CBS**: HF and correlation extrapolated separately from
   the TZ/QZ pair (Helgaker two-point, Psi4's `cbs()` conventions).
3. **DLPNO cross-check (Psi4, CPU)**: `energy('dlpno-ccsd(t)')` with
   `PNO_CONVERGENCE TIGHT`, `freeze_core true`, same basis pair.
   Report ΔE‡(DLPNO) − ΔE‡(canonical); acceptance: |Δ| ≤ 2 kJ/mol on
   the barrier (the water-dimer smoke suggests well under that).
4. **Verdict table**: ΔE‡(elec) at CCSD(T)/CBS vs B3LYP/def2-SVP (the
   Phase-1/2 tier) vs the A2 production tier (ωB97M-V/def2-TZVPD) —
   the "which functional do we trust" number, per SURVEY §6.4.
5. **Provenance**: every energy into the run's `store.sqlite` with
   engine + version + basis + settings, like every quarry number.

### Executable adapter and fail-closed gate (2026-08-24)

`qm/scripts/cc_calibration.py` now owns both one-geometry jobs plus the
TZ/QZ CBS summarizer. Each job is atomic and idempotent on
engine+basis+input-SHA; the summarizer requires the complete
reactant/TS × TZ/QZ matrix and reports the 2 kJ/mol DLPNO gate. The adapters
were exercised independently on water/cc-pVTZ: Psi4 1.11 TightPNO
DLPNO-CCSD(T) `-76.3323355604 Eh`; ByteQC canonical CCSD(T)
`-76.3325745245 Eh` (PySCF 2.5.0, frozen O 1s).

### A2a directive-adjusted production route (2026-09-05)

The canonical reactant/cc-pVQZ slow-mode calculation eventually completed, but
the matching TS/QZ calculation was explicitly cancelled rather than spending a
second multi-day canonical triples campaign. `summarize-focal` therefore accepts
exactly three canonical receipts (reactant TZ/QZ and TS TZ), requires the full
four-receipt TightPNO reactant/TS × TZ/QZ matrix, and refuses a canonical TS/QZ
selector. The production focal point is

`canonical barrier(TZ) + [DLPNO barrier(CBS) - DLPNO barrier(TZ)]`.

The scientific agreement gate remains independent and like-for-like at TZ:
`DLPNO barrier(TZ) - canonical barrier(TZ)`. For A2a this is
`+0.245700 kJ/mol`, passing the `2 kJ/mol` gate; the focal barrier is
`132.960133 kJ/mol`. The large reactant-only QZ absolute method offset is retained
as anchor provenance and is never substituted for a barrier delta. Receipt
arithmetic, exact role geometry hashes, electronic state, engine, and basis are
validated before publication. Durable receipt hashes are in the A2a card.

Do not launch the si-neutral matrix until the candidate transition state passes
A2's repaired nonzero-step full IRC. The first production attempt proved the
banked TS reaches the associative basin in both directions, so it is ineligible
for a functional ranking or coupled-cluster barrier. READY card
`A2a-si-neutral-production-path-rebuild` owns the route reconstruction; this is
a scientific gate, not a missing-engine blocker.

Larger calibration points later (30–80-atom ladder clusters) follow the
same recipe; those are where DLPNO (near-linear scaling) carries the
load and ByteQC spot-checks at TZ where memory allows (rental GPUs
extend the envelope — see the H200/B200 note above).
