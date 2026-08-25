# A2a-si-neutral-production-path-rebuild — replace the invalid banked saddle

- status: active
- track: A (geochemistry)
- priority: P1
- machine: workstation
- depends: — (A2 infrastructure and exact r2SCAN-3c minima are already pushed)
- claimed-by: hermes-workstation

## Objective

Rebuild the si-neutral production-tier hydrolysis path after A2 proved that the
banked B3LYP saddle is not the reaction saddle on the exact r2SCAN-3c surface.
Start from A2's hash-pinned zero-index r2SCAN-3c reactant/product, construct an
endpoint-spanning path (relaxed scan and/or persisted pre-climb + climb CI-NEB),
then refine the interior crest with a reaction-aligned Sella search. If the
production surface is sequential associative, split it explicitly into
reactant↔intermediate and intermediate↔product elementary steps and verify both;
do not force one concerted R↔P saddle. Do not retry either the undirected or
endpoint-mode search from the rejected banked TS: both reconverged the same
associative-basin saddle.

Once a correct saddle exists, complete the si-neutral slice of A2: exact
finite-difference composite Hessian, repaired full Gonzalez--Schlegel IRC,
wB97M-V/def2-TZVPD+SMD and B3LYP-D4 comparison single points, and the open-source
CCSD(T)/CBS calibration triangle (ByteQC canonical + Psi4 TightPNO DLPNO).

## Context

- `qm/scripts/production_energetics.py` — resume-safe production driver
- `qm/scripts/cc_calibration.py` — verified Psi4/ByteQC job + CBS summarizer
- `qm/quarry/pipeline.py` — exact r2SCAN-3c D4+gCP and FD-Hessian fallback
- `qm/quarry/ts.py` — ASE-3.29-safe Sella full IRC
- `/mnt/data/vsletten/dissertation-data/task207-a2-production-20260824/`
  (workstation evidence root; manifest SHA-256
  `efdcfc48c533dbadc1a30023f3a539f607711df0cd5cd354ff50fca3ba3556d1`)

## Constraints

- One heavy GPU campaign at a time; driver-enforced threads <=16, niceness, and
  durable log. Use the established pipeline/Honcho isolation + bounded dead-man
  restoration protocol.
- Persist every NEB stage before climb/refinement. Terminal receipts bind code,
  settings, geometry, and artifact hashes. Never weaken a basin, index, or
  convergence gate to rescue a number.
- Large run artifacts remain outside Git. Commit/push implementation and card
  checkpoints incrementally.

## Acceptance

- r2SCAN-3c reactant and product remain exact zero-index minima; every
  replacement TS has exactly one significant imaginary mode from a
  central-difference Hessian over the full composite gradient. The index and
  reaction-mode assignment are stable under a bounded finite-difference step
  check. A sequential route reports and gates R↔I and I↔P separately.
- Full Gonzalez--Schlegel IRC performs nonzero steps in both directions for
  every accepted TS and terminates in the exact typed neighboring basins,
  including physical-H ownership, heavy-atom topology, and minimum-pair gates.
  The overall route must connect neutral reactant `(False,True,False)` to
  hydrolyzed product `(True,False,True)`; a zero-step or same-basin path fails.
- Production ΔG‡ and SVP→r2SCAN-3c→wB97M-V/B3LYP-D4 shifts are stored with
  geometry/method provenance; no number derives from the rejected banked TS.
- ByteQC canonical and Psi4 TightPNO cc-pVTZ/cc-pVQZ receipts yield a focal-point
  CBS barrier; `|DLPNO - canonical| <= 2 kJ/mol` or the gate fails honestly.
- Targeted + complete QM tests, whole-tree Ruff check/format, CLI smokes, and
  `git diff --check` pass. Result updates A2, CALCULATIONS.md, PLAN, and STATUS.

## Progress

- 2026-08-25 — atomic claim is remotely durable on
  `agents/A2a-si-neutral-production-path-rebuild`. Commit `13311ff` adds the
  sequential production-path driver and CPU-only orchestration contracts. It
  never accepts the rejected banked saddle as input: the exact archived
  r2SCAN-3c R/I/P references are independently reoptimized and minimum-gated,
  then R↔I and I↔P each receive persisted pre-climb/climb CI-NEB,
  reaction-aligned Sella, two-step central-difference Hessian stability, exact
  typed full-IRC endpoint gates, production single points, and one provenance
  store. The predeclared numerical contract and 59 focused A2/A2a/TS/CC tests,
  Ruff, CLI smoke, and diff checks are green. The bounded workstation campaign
  is running under `/mnt/data/vsletten/dissertation-data/task208-a2a-path-rebuild-20260825/`
  with extraction/Honcho GPU consumers isolated and a ten-hour dead-man restore.
- 2026-08-24 — created by hermes-workstation from A2 live production evidence.
  Exact r2SCAN-3c FD-Hessian indices are reactant/product/TS `0/0/1`, but the
  repaired 260-step full IRC closes in the associative basin on both sides.
  Endpoint-mode directed Sella changed the saddle by only 0.0124 Å RMS and
  reconverged the same basin. Route reconstruction, not another saddle retry,
  is the next bounded scientific move.
