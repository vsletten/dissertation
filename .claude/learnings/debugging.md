# Debugging

### Saddle searches slide to trivial saddles from uncrossed-ridge guesses (2026-08-18)
First live hydrolysis TS hunt: the r(Si-Ow) scan ended (1.9 A) with energy
still rising, its endpoint became the Sella guess, and the search relaxed
back to a 208i cm^-1 water-wag saddle — mechanically valid (1 imaginary
mode, IRC converges) but chemically wrong (dG‡ 2.5 vs ~120 kJ/mol lit.).
Diagnosis signature: the saddle's driven coordinate far past the guess and
both quick-IRC ends identical to the reactant complex. Fixes in quarry:
scan_to_maximum (extend until interior max) + escaped-channel guard in
the phase-1 driver. A verified saddle is not necessarily the right saddle.

### CuPy pool fragmentation OOMs multi-hour GPU campaigns (2026-08-19)
First Phase-2 pilot: stage 1 (2 h of geomeTRIC steps on a 66-atom cluster)
ran fine at ~7 GB, then the first stage-2 DF-K gradient died with
cudaErrorMemoryAllocation on the 24 GB card — the CuPy memory pool had
fragmented across thousands of allocations. Fix: free_all_blocks() on the
default + pinned pools at every stage boundary and scan point
(phase2_ladder.trim_gpu_pool). Also: `cmd | tee log` reports tee's exit
code — a crashed driver looks like exit 0 unless `set -o pipefail`.

### GPU4PySCF 1.8.1 DF-UKS Hessian can assert on contiguity (2026-08-19)
Open-shell D2a saddles optimized normally on the GPU, but the analytic
density-fitted UKS Hessian aborted inside `gpu4pyscf.hessian.uks` at
`assert wv.flags.c_contiguous`. Quarry now retries only that failed Hessian on
CPU; geometry, scans, and saddle optimization remain GPU-first. The affected
astro systems contain 3–5 atoms, so this is a bounded tiny-molecule fallback,
not a license to move campaigns onto CPU.

### Composite gradients and meta-GGA Hessians must share one surface (2026-08-24)
The first A2 r2SCAN-3c run attached D4+gCP to SCF energies, but PySCF's
geomeTRIC bridge called `mf.nuc_grad_method()` directly and bypassed quarry's
patched derivative helper. It therefore optimized plain-r2SCAN forces on
composite energies. Patch the derivative factory itself, then independently
check the final composite force. GPU4PySCF 1.8.1's analytic r2SCAN Hessian also
returned 7--11 large imaginary modes for force-converged free-cluster minima;
a central-difference Hessian over the exact full composite gradient returned
zero, as physically required. Never weaken the index gate to accept the fast
Hessian; preserve the failed receipt and use the same-surface gradient fallback.

### Sella 2.5 IRC silently no-ops under ASE 3.29 without an API bridge (2026-08-24)
Sella's `IRC` implements its load-bearing stop rule in `converged()`, including
`first => False` so the mass-weighted kick must happen. ASE 3.29 changed
`Optimizer.irun()` to call `gradient_converged()` directly; the inherited
force-only rule sees the transition state already below `fmax` and returns at
step zero in both directions. Signature: a four-line IRC log containing only
`IRC: 0`, identical forward/reverse coordinates, and a falsely successful
return. Bridge `gradient_converged()` to Sella's `converged()` and regression
assert both directions reject the initial force-only decision.

### Shared-atom FixBondLengths residuals need an independent gate (2026-08-25)
ASE 3.29's `FixBondLengths` can report a converged projected optimization while
three coupled distances sharing atoms still miss their requested values by
`1.2e-4--1.7e-4 A`. Repeated `set_positions()` projection did not monotonically
close the residual in the A2a crest conditioner. Keep an independent final
distance check and fail closed; use a residual-controlled internal-coordinate
conditioner (or a hash-bound previously accepted crest) rather than loosening
the scientific gate.

### Active-subspace Sella can escape a converged local CI-NEB crest (2026-08-25)
On A2a I↔P, a converged seven-image pre-relax/climb and simultaneous internal-
distance conditioner reached `3.95e-6 A`, but unconstrained active-core Sella
then climbed roughly six Hartree in 27 steps before the SCF failed. A finite
single-point failure is not a retry license and the optimizer trajectory is not
a candidate checkpoint. Resume from the hash-bound conditioned crest, constrain
the trust/displacement envelope around the local band tangent, and persist only
a converged saddle; never restart the expensive band or weaken the SCF/index
gates to rescue the runaway search.

### A local trust wall can stop escape without creating a saddle (2026-08-26)
A2a I↔P bounded Cartesian Sella to the exact adjacent climb-image radii: the
physical PES stayed untouched inside `0.344314 A`, a flat-bottom restoring wall
acted beyond it, a hard guard refused pre-PES evaluations past `0.545962 A`,
and Sella's step radius was capped at `0.03 A`. The 300-step run never escaped
(`0.375573 A` maximum; `0.333471 A` final) but settled into a limit cycle near
`0.3--0.6 eV/A`. A trust wall is localization, not evidence of a stationary
point: reject wall-supported/nonconverged geometries and pivot to a bounded
full-system local eigenvector-following or dimer search rather than increasing
the radius, weakening `fmax`, or replaying the same active-subspace attempt.

### Early-TS crests can sit under the scan peak tolerance (2026-08-28)
H+CO's crest at r≈1.9 Å cleared its neighbor by only 2e-5 Ha — under
`scan_to_maximum`'s 5e-5 anti-noise tolerance — so the scan extended into
the product well and raised ScanNoMaximumError despite a genuine interior
maximum 9.6 kJ/mol above the entrance. For flat-approach early TSs, fall
back to the interior *global* maximum (guarded against endpoint maxima and
noise-level rises) and let Sella + the chemical-mode gate be the acceptance.

### Wet-cluster stationary points carry spectator imaginaries <100 cm⁻¹ (2026-08-28)
At production convergence, physisorbed water librations sit numerically
imaginary (observed live: 52.8i beside a healthy 678i TS mode; 76.5i and
32.8i on pre-reactive complexes). Flat noise floors sized for gas-phase
work (30–50 cm⁻¹) fail correct wet structures. Tier the gates: gas strict,
wet sites 100 cm⁻¹ spectator tolerance with the full imaginary list
disclosed in results; keep the ≥200 cm⁻¹ single-chemical-mode gate.

### Pre-optimize hand-built water networks before constrained scans (2026-08-28)
A first *pinned* scan point that also has to settle a seeded water network
walks through strained geometries and can diverge the SCF (DFUKS scanner
death 35 geomeTRIC steps into h-h2co-h2-hco-2w's first point). Relax wet
seeds unconstrained first (checkpointed), then drive the coordinate.

### Pre-optimize terminated ladder complexes before production DFT (2026-09-05)
The first Osa-neutral n=1 ladder attempt passed geometry/collision gates but its
raw crystallographic cluster plus attacker left the B3LYP/def2-SVP DFRKS
convergence basin after one large geomeTRIC move (0.340 A maximum displacement).
The initial gradient improved and constraints stayed within 0.0011 A, then the
next electronic gradient exhausted the existing 150-cycle SCF bound. Preserve
that failure and settle only the reactant complex with checkpointed constrained
HF/STO-3G before production optimization; do not mask it by raising the SCF
bound, weakening convergence, or treating the failed trajectory as a checkpoint.
