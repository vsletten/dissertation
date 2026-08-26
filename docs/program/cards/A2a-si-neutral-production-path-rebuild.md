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

- 2026-08-26 12:30 PDT (hermes-custom-build-001; profile=workstation) — TASK-227 full-system localization implementation is crash-durable at pushed commit `7440a7df996a4ba06d8f278ac769dd2e836e5e61`. It adds a dedicated ASE 3.29 dimer path in which every non-frozen atom remains movable while the exact six-atom climb-coordinate envelope retains the `0.03 A` step cap and refuses before physical PES evaluation outside the hash-bound adjacent-image guard. Candidate persistence additionally requires negative curvature, physical full-system `fmax < 0.02 eV/A`, and final return to the flat restraint-free region. Focused A2a/TS verification passed 54 tests; whole-QM Ruff check and format gates pass. The single bounded r2SCAN-3c attempt is staged under `cleavage/full-system-local-dimer-v1`; runner SHA-256 is `8877101a4386e3c2892f26ce79a9329e94c25540ed845d526a514d4b8a1ec51b`. No candidate, Hessian, IRC, production single point, CC job, or barrier exists yet.
- 2026-08-26 05:08 PDT — **fourth bounded continuation prevents the cleavage
  runaway but does not localize a saddle; acceptance remains red.** Pushed
  commits `321f870`, `267dff8`, `826b99a`, and `9b9edb2` add a flat-bottom
  Cartesian trust envelope anchored to the exact adjacent climb-band images,
  a pre-PES hard guard, a `0.03 A` Sella step cap, legacy addition-receipt
  compatibility, and a mandatory cuTENSOR preload for GPU runs. Core focused
  verification passed 51 A2a/TS tests; the final cuTENSOR/receipt slice passed
  10 focused tests, Ruff, format, diff checks, and a live preload probe with
  gpu4pyscf's native cuTENSOR engine selected. The exact cleavage neighbor
  radii were `0.344314/0.545962 A`. One full 300-step r2SCAN-3c attempt stayed
  bounded (`0.375573 A` maximum, `0.333471 A` final; never crossed the
  `0.545962 A` pre-PES guard) instead of climbing six Hartree, but cycled at
  roughly `0.3--0.6 eV/A` and failed the convergence gate. No candidate,
  Hessian, IRC, production single point, CC job, or barrier was emitted.
  Durable evidence is **461 files / 2,645,283 bytes**, manifest SHA-256
  `8542944ade6f2c5aa339cf9dab487fd72a12d3fe712a9437b0ea87831a4c106e`.
  Both email pipelines and Honcho containers are restored, API health is `ok`,
  the dead-man is inactive, Ollama is empty, the GPU is idle, and no QM process
  remains. **Do not rerun this active-subspace trust search.** The next bounded
  scientific move is a separately tested full-system local eigenvector-following
  or dimer localization inside the same hash-bound climb neighborhood, followed
  by the unchanged Hessian/IRC gates.
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
- 2026-08-25 11:06 PDT — **bounded scientific continuation checkpoint; acceptance
  remains red.** Exact r2SCAN-3c R/I/P references independently pass zero-index
  finite-difference Hessians, and the associative intermediate is now a proved
  minimum rather than an inferred IRC endpoint. A converged seven-image R↔I
  CI-NEB has the exact electronic profile `0.0, 2.043, 10.744, 17.177, 51.021,
  125.472, 108.622 kJ/mol`: a genuine interior crest at image 5 before the
  associative image 6. Two bounded Sella refinements were rejected without
  weakening the gate: endpoint-vector mode produced a stable `147.72 cm^-1`
  spectator-OH saddle with overlap `0.006/0.006`; three-coordinate-conditioned
  local-band refinement produced a stable `101.89 cm^-1` spectator-H14 saddle
  with overlap `0.037/0.037`. Tight-climb continuations also failed honestly:
  ODE worsened `fmax 0.0780→0.1041 eV/A` over 12 steps; BFGS worsened
  `0.1700→0.3032` over 17. No full IRC, production single point, or CC job ran
  after that first elementary-step gate stayed red. Commit `ab9bf00` and every
  earlier checkpoint are pushed. Durable external evidence is **259 files /
  869,079 bytes**, manifest SHA-256
  `3852d6b5feea82e57f8b15eb87fc93150590fe711edba36da2ddd3f0b03a7efe`.
  Shared extraction/Honcho services are restored, Ollama is empty, and no QM
  process remains. **Continuation order:** reuse the proved minima and persisted
  bands; implement a spectator-torsion-aware active-subspace or internal-
  coordinate reaction-mode refinement. Do not rerun the rejected banked TS,
  recompute the loose band, weaken the mode gate, or start I↔P/high-level/CC
  work before R↔I closes.
- 2026-08-25 12:13 PDT — **second bounded continuation remains fail-closed.**
  Commits `b61efb6`, `630de91`, and `6778824` add and test an explicit
  active-subspace Sella seed, spectator-space reconditioning, and final
  full-system Sella release; 49 focused A2a/TS tests, Ruff, format, and diff
  checks passed. Three isolated GPU launches reused the exact minima and the
  converged loose band, but all stopped before active-subspace Sella: the
  three shared-atom `FixBondLengths` constraints reported final distance
  residuals `1.186e-4`, `1.216e-4`, and `1.749e-4 A`, above the predeclared
  `1.0e-4 A` gate. Reapplying ASE's coupled projection did not close the
  residual, so the gate was not loosened and no new TS/Hessian/IRC,
  I↔P, production single point, or CC job ran. Shared extraction/Honcho
  services are restored, Ollama is empty, and no QM process remains. Evidence
  is **266 files / 1,074,740 bytes**; manifest SHA-256
  `c8d143a849e6e225ab79844dadedeabb28ce7607dd3f866446363d22dd00042c`.
  **Continuation order:** replace the shared-atom Cartesian bond projector
  with independently residual-gated internal-coordinate conditioning, or
  hash-bind and reuse the previously accepted conditioned crest; then exercise
  the active-subspace→spectator-relaxation→full-release chain. Do not loosen
  the distance/mode gates or restart the loose NEB.
- 2026-08-25 14:29 PDT — **third bounded continuation closes R↔I and preserves
  a new cleavage-search failure honestly.** Commit `8ced110` replaces ASE's
  sequential shared-atom projector with simultaneous order-zero Sella bond
  equations; the shared-pin analytical regression plus 50 focused A2a/TS tests,
  Ruff, format, and diff gates pass. The addition conditioner now closes its
  three target residuals to `2.15e-6 A`, and the full-system R↔I saddle has one
  stable reaction mode at `714.831/714.821 cm^-1` across the two Hessian steps
  (mode cosine `0.99999999998`, local-tangent overlaps `0.854/0.854`). Full IRC
  displacements are `2.075/0.677 A` and terminate in the exact typed neutral
  reactant `(False,True,False)` and associative intermediate `(True,True,True)`
  basins. The independently persisted I↔P seven-image pre-relax and ODE climb
  also converge (`51` and `25` steps), and its internal conditioner reaches
  `3.95e-6 A`; however, the subsequent active-subspace Sella search escaped
  uphill and failed SCF at step 27 before any cleavage candidate was persisted.
  No cleavage Hessian/IRC, production single point, CC job, or barrier ran.
  Evidence is **457 files / 1,954,856 bytes**; manifest SHA-256
  `82af53d69809e8b688db183af5a914886178c8f82f928883f15bb3ffa7d01fb6`.
  Both email pipelines and Honcho containers are restored, API health is `ok`,
  Ollama is empty, the dead-man is inactive, the GPU is idle, and no QM process
  remains. **Continuation order:** reuse the exact cleavage `climb-final` band
  and `conditioned-crest-v3`; bound the active-subspace Sella release so it
  cannot outrun the local crest (or recover fail-closed from the last finite
  local geometry), then rerun only cleavage Sella→Hessian→IRC. Do not recompute
  R↔I, either loose band, or weaken any gate.
- 2026-08-24 — created by hermes-workstation from A2 live production evidence.
  Exact r2SCAN-3c FD-Hessian indices are reactant/product/TS `0/0/1`, but the
  repaired 260-step full IRC closes in the associative basin on both sides.
  Endpoint-mode directed Sella changed the saddle by only 0.0124 Å RMS and
  reconverged the same basin. Route reconstruction, not another saddle retry,
  is the next bounded scientific move.
