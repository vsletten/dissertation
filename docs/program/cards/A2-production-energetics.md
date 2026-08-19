# A2-production-energetics — production-tier energies + ORCA calibration

- status: blocked (ORCA install by Victor; A1b for the full reaction set)
- track: A (geochemistry)
- priority: P1
- machine: workstation
- depends: A1b
- claimed-by:

## Objective
Upgrade Phase-1 barriers to the survey protocol: r2SCAN-3c geometries
(needs pyscf-dispersion — already in the venv), wB97M-V/def2-TZVPD + SMD
single points on all stationary points, D4 dispersion; ORCA
DLPNO-CCSD(T)/CBS spot calibration on si-neutral (the "which functional
do we trust" half-day, SURVEY §6.4). Full Gonzalez-Schlegel IRC on the
soft-mode saddles replaces quick-IRC. Report the barrier shifts
SVP→production per reaction.

## Context
qm/SURVEY.md §6 (the tiered protocol), §2.3 (ORCA); quarry pipeline
supports solvent="smd" already. ORCA: registration-gated download,
NEVER committed/containerized (license — Appendix A).

## Acceptance
- Production ΔG‡ table for all Phase-1 reactions with method provenance;
  DLPNO calibration number for si-neutral; IRC-verified saddles.
- Unblocks: mark done only when Victor confirms ORCA installed.

## Progress
- 2026-08-18 — card created (fable).
