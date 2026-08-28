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

- 2026-08-28 00:58 PDT (hermes-custom-build-001; profile=workstation) — **the production 9x9 scan is durably in flight at 3/81 completed cells; this bounded worker continuation did not classify the surface or emit a saddle/barrier.** The first launch failed before any scientific work because a prior test-only `uv sync` had removed the workstation GPU extra; `uv sync --extra dev --extra gpu` restored CuPy/gpu4pyscf/cuest/cutensor and a live GPU bootstrap probe passed. The second launch acquired the QI2 lease and completed exact adjacent cells `(0,0)`, `(0,1)`, and `(0,2)` with finite r2SCAN-3c energies, typed associative topology, H16 still water-owned, Si--Obr still bonded, and maximum target residual `4.35e-6 A` against the unchanged `1e-4 A` gate. The scan manifest SHA-256 is `3b2e7b8fd7ab112e28ee1c3c7506cae842164a4a31350fd2c2a6b85535337251`; the three cell-receipt hashes are `98d3c6721797c3dfde6340ba24cc430494d29a7b87d2341db619378427acc397`, `14724c1761f0be7d8713444d6c737e73aa377e60b857564491dc1b790c375518`, and `07fb24eff76fb8692a00e8cc044948dc84ad5aece1aaa547ae22936fb3817fc1`. A watchdog timer restarted the intentionally stopped email pipeline during the run, although `ollama ps` remained empty and nvidia-smi showed only the QM Python process; the worker stopped at the durable cell boundary instead of accepting ambiguous isolation. External launcher `run-task220-coupled-scan.sh` is now watchdog-aware, hash `afeeca819bf8d0c74ea112cc03d708347dee654c8ec8c6ec1d20948d946eabd5`, and the next continuation resumes at `(0,3)` from the exact existing branch and evidence root. After bounded stop, both email pipelines, the watchdog timer, Honcho API/deriver, and API health are green; Ollama and the GPU are idle. No Hessian, IRC, high-level single point, CC calibration, final barrier, or PR exists yet.

- 2026-08-27 23:57 PDT (hermes-custom-build-001; profile=workstation) — **the commissioned coupled-coordinate redesign is implemented and mechanically green; the production 9x9 scan has not run yet.** The branch was first reconciled with current `origin/main` at `602a3ad`, preserving the merged PR #90/QI2 histories. New `qm/scripts/a2a_cleavage_scan.py` provides an 81-cell adjacent serpentine r2SCAN-3c scan over `d(Si1-Obr0)` and the physically relevant water proton `d(H16-Obr0)`, exact endpoint/settings/seed/target fingerprinting, atomic per-cell geometry/energy/topology receipts, strict resume-drift refusal, and a complete-grid minimax classifier whose precedence is topology-proved proton-first minimum → interior crest above 2 kJ/mol → barrierless shelf. Self-review fixed two false-green risks before checkpoint: extension cells cannot become a low-energy I→P detour, and proton-first classification requires actual H16 ownership plus intact Si-Obr rather than grid-index inference. `qm/A2A_CLEAVAGE_SCAN.md` is the decision matrix. Real archived endpoints pass the fixed mapping (`Si1/Obr0/Ow15/H16`; H16 contracts `1.127276 A`, H17 changes only `0.013367 A`); the archive-only old crest diagnostic is `+0.972660 kJ/mol` relative to I, below the predeclared threshold but not a substitute for the scan. Verification passed 11 new focused scan tests, 38 focused A2a/production/CC tests, the complete QM suite (**451 passed, 1 skipped**), whole-tree Ruff check/format (including a one-line wrap for an overlong QI2 lease message imported from current main), compile, CLI help, real-endpoint mapping, and diff checks. No new saddle, Hessian, IRC, production single point, CC result, or barrier was emitted. Next continuation runs the bounded production grid under the QI2 GPU lease and service-isolation/dead-man contract, then implements only the classifier-selected follow-through.

- 2026-08-26 14:58 PDT (hermes-custom-build-001; profile=workstation) — **the exact-envelope full-system local dimer family is exhausted; no cleavage saddle exists from this localization.** The corrected v2 run consumed all 120 bounded `0.03 A` translations with every one of 484 rotation-curvature estimates positive (minimum `+0.0274`, final `+0.0372 eV/A²`). It remained inside the exact `0.545962 A` pre-PES guard (`0.379622 A` maximum/final local radius) but outside the `0.344314 A` flat restraint-free radius, so no candidate was persisted and no Hessian, IRC, production single point, CC job, or barrier ran. Exact implementation head `95ea8df7f1fb7939dc2ecf0864f2e5a5a79c3543` also separates translation/rotation log files after self-review caught ASE's two file handles overwriting each other. Final verification is 425 complete QM tests plus the 54 focused A2a/TS slice; whole-QM Ruff, format, CLI, and diff gates pass. Durable evidence is **478 files / 3,088,395 bytes**, manifest SHA-256 `9e163528bb43b503b1f2f45d33cd8cc064d3effe66490ab4f13bc869705b6698`; v2 failure receipt SHA-256 `940d49067cc680a2233561da06d8fe07e7b3453069622efe99bd2a98d639987a`. Both email pipelines and Honcho are restored healthy, the dead-man is inactive, Ollama is empty, the GPU is idle, and no QM process remains. A2a remains scientifically blocked: do not run another local active-core/Sella/dimer retry from this conditioned crest without a new mechanism/path decision.
- 2026-08-26 12:58 PDT (hermes-custom-build-001; profile=workstation) — self-review stopped and quarantined the first live dimer attempt after 22 steps: v1 had incorrectly recomputed the envelope from the conditioned crest (`1.533078/1.642800 A`) instead of retaining the exact persisted climb-peak neighbor radii (`0.344314/0.545962 A`). No candidate was persisted; maximum observed translation was `0.030000 A`, maximum local displacement `0.561047 A`, final curvature remained positive `+0.055502 eV/A²`, and projected `fmax` remained red at `0.0387 eV/A`. Fix `375075d0e40a80ee011c6bed9a00eb340fe41a33` is pushed, 54 focused tests plus whole-QM Ruff/format/diff gates pass, and live preflight now reproduces the exact `0.344314/0.545962 A` contract. The corrected v2 runner SHA-256 is `a3bebe7196b06e7363ce2f21019e5e607a25fdac5d017d754d719da476afee69`; v1 failure evidence is preserved under `cleavage/full-system-local-dimer-v1/`.
- 2026-08-26 12:30 PDT (hermes-custom-build-001; profile=workstation) — TASK-227 full-system localization implementation is crash-durable at pushed commit `7440a7df996a4ba06d8f278ac769dd2e836e5e61`. It adds a dedicated ASE 3.29 dimer path in which every non-frozen atom remains movable while the exact six-atom climb-coordinate envelope retains the `0.03 A` step cap and refuses before physical PES evaluation outside the hash-bound adjacent-image guard. Candidate persistence additionally requires negative curvature, physical full-system `fmax < 0.02 eV/A`, and final return to the flat restraint-free region. Focused A2a/TS verification passed 54 tests; whole-QM Ruff check and format gates pass. The single bounded r2SCAN-3c attempt is staged under `cleavage/full-system-local-dimer-v1`; runner SHA-256 is `9f515fada0b8a21dfeffba568bb22f9bb29832ec1de3d9cc7b3d3e240614377e` and binds the exact implementation tree at `7440a7df996a4ba06d8f278ac769dd2e836e5e61`. No candidate, Hessian, IRC, production single point, CC job, or barrier exists yet.
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
