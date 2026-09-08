# D2c-instanton-tier — deep-tunneling rates that flip the D2b gate

- status: active
- track: D (astrochemistry)
- priority: P2
- machine: workstation
- depends: D2b
- claimed-by: hermes-custom-build-001

## Objective

D2b's NO-GO residual is isolated to tunneling method: asymmetric
Eckart is 1.5–3 orders low against instanton-class anchors below
~59 K, because it inherits the barrier *width* from a single imaginary
frequency (and PWB6K's width is demonstrably too stiff — T_c 149 K vs
Andersson's 79 K for H+CO). Add one instanton-grade tunneling tier to
quarry: small-curvature tunneling (SCT) along the existing
`full_irc`/`quick_irc` path machinery, or a ring-polymer instanton on
the driver's stationary points. Benchmark reaction-by-reaction against
Song & Kästner 2017 (H2CO+H channels) and Andersson 2011 / Simons 2020
(H+CO), reusing the D2b campaign's checkpointed geometries and CC
barriers (`qm/runs/D2b-explicit-surface-rates/`, campaign worktree).

## Acceptance

- A tunneling correction beyond Eckart implemented against the
  existing IRC machinery, unit-tested, with the same receipt
  discipline as D2b.
- H+CO and both H2CO channels within the documented literature spread
  at 12–20 K (the Simons plateau ±1 order), with the branching ratio
  defensible against the experiment-benchmarked network.
- The D2b validity floor restated with the new tier; if the gate now
  passes, emit the machine-readable D3b rate table D2b withheld —
  that table unblocks D3b.
- Bounded, GPU-first where applicable, thread-capped, teed logs.

## Progress

- 2026-08-28 — created by D2b closeout: the campaign validated the
  barrier tier (1w abstraction barrier = S&K's surface value to 2 K)
  and named Eckart-vs-instanton as the sole remaining gap.
- 2026-09-06 22:38 PDT (hermes-custom-build-001; profile=workstation) —
  DEVIATION: implemented and unit-tested the reusable Liu/Pilgrim SCT action,
  effective-mass, and I1/I2/I3 thermal-integration kernel (33 focused/rates
  tests pass). An independent adversarial review rejected the first sparse
  quick-IRC/relaxed-scan benchmark as scientifically non-SCT and its site-envelope
  gate as false-green; the provisional D3b table and driver were deleted rather
  than laundered as a result. Acceptance still requires a mass-scaled full IRC,
  transverse frequencies/eigenvectors at path points, a vibrationally adiabatic
  potential, matched-site branching, and immutable tracked provenance. No D3b
  rate table has been emitted.
- 2026-09-07 03:48 PDT (hermes-custom-build-001; profile=workstation) —
  NUMERICAL FOUNDATION REVIEWED: added a fail-closed, dry-run-only D2c campaign
  preflight plus reusable molecular reaction-path primitives. The path now uses
  TS-anchored outward mass-weighted Kabsch transport, exposes exact frame rotations,
  projects molecular tangents/curvature in the local rigid-horizontal space, rotates
  Cartesian Hessians coherently, extracts exactly 3N-7 positive transverse modes,
  builds the vibrationally adiabatic potential, and applies continuously regularized
  mode-resolved Liu/Pilgrim effective mass. Preflight hash-binds all four atom mappings,
  frozen geometries, GPU4PySCF/CuPy/CUDA identity, exact atomic-unit Hessian convention,
  campaign bounds, and stage receipts under a crash-recoverable exclusive run-root
  claim. Two independent correction rounds closed all critical/important scientific
  review findings; final review approved with one non-blocking native PySCF Hessian-layout
  note for the production driver. Whole-QM verification is 868 passed, 1 skipped plus
  whole-tree Ruff/format; no IRC/Hessian production compute or D3b table ran. Acceptance
  remains red until the resumable production stages execute and matched-site 12-20 K
  literature/branching gates pass.
- 2026-09-07 05:02 PDT (hermes-custom-build-001; profile=workstation) —
  PRODUCTION IRC ADAPTER CHECKPOINT: added a path-preserving Sella 2.5 IRC adapter
  that records the exact TS and every outer point for both algebraic directions from
  one shared-Hessian optimizer, accepts receipt-bound isotopic masses before optimizer
  construction, preserves the legacy endpoint API, and fails closed on zero-step,
  exhaustion, non-finite, accounting, point-bound, terminal-convergence, and endpoint-
  force violations. Deterministic fake-optimizer coverage plus the whole QM suite are
  green (`879 passed, 1 skipped`; Ruff/format/diff clean). Independent scientific review
  found the frozen D2b H+CO inputs retain 13.5-52.8i cm^-1 spectator modes, so production
  launch remains prohibited until the driver performs a fresh first-order-saddle/index,
  reaction-mode/tangent, and typed endpoint-basin gate (or refines the saddles under a
  new identity). No production calculator call, IRC receipt, Hessian, or D3b table ran.
- 2026-09-07 07:54 PDT (hermes-custom-build-001; profile=workstation) —
  TRANSITION-STATE HESSIAN GATE CHECKPOINT: added a fresh native PySCF Hessian
  evaluator with explicit `(atom_i, atom_j, xyz_i, xyz_j)` canonicalization, immutable
  gradient/Hessian payloads, exact GPU-to-CPU fallback provenance, and geometry/settings
  fingerprints. The campaign now applies a full nonlinear `3N-6` index gate before IRC:
  physical `fmax < 0.02 eV/A`, exactly one mode below `-1e-8 Eh/Bohr^2/amu`, no
  zero/noise-floor or additional negative modes, reaction frequency at least `200 cm^-1`,
  and an unstable-mode vector derived from and hash-bound to mapped reactant/product
  basin geometries. Four adversarial review rounds closed Hessian-layout, mutable-payload,
  self-attested-vector, backend-provenance, and non-finite-overflow false greens; final
  review passed. Exact pushed head is `c5ae2b990fd23323d970f15cac07278ee01c68df`;
  whole-QM verification is `900 passed, 1 skipped` with Ruff/format/diff clean and the
  historical hash-pinned `quarry/pipeline.py` unchanged. No production calculation,
  endpoint-basin classification, per-path Hessian receipt, or D3b table ran. The next
  continuation must implement typed endpoint classification plus atomic resumable route
  and per-point Hessian publishing before any bounded production launch.
- 2026-09-07 11:26 PDT (hermes-custom-build-001; profile=workstation) —
  TYPED PATH/RECEIPT CHECKPOINT, REVIEW-RED: exact pushed head
  `dc86dd80ba22c15cc91916c3ce263db2e5cda235` adds four-route exact-graph endpoint
  typing, chemistry-based two-direction Sella orientation, no-clobber atomic path
  publication, and strict crash-resumable per-point native-Hessian plus aggregate
  receipts. Three spec-review rounds closed stale-temporary, destination-race,
  bool/int, precommit-poisoning, and extra-direction false greens. Whole-QM gates are
  `961 passed, 1 skipped` with Ruff/format/diff clean; no production calculator or D3b
  output ran. Final independent scientific review rejected production launch because
  the publisher does not yet bind/read the real preflight + TS/IRC receipt chain,
  enforce all frozen IRC bounds (`max steps/points`, terminal `fmax`, step size,
  tangent overlap), or require fresh Hessian energies to reproduce path energies.
  Next continuation must close those three gates adversarially before enabling any
  production execution. The checkpoint is durable infrastructure, not accepted SCT
  evidence.
- 2026-09-06 23:33 PDT (hermes-custom-build-001; profile=workstation) —
  INPUT GATE CLOSED: froze the four predeclared direct-CC one-water D2b routes
  (both H+CO orientations plus both H2CO channels) into a tracked 29-file,
  71,934-byte bundle whose manifest binds every byte/geometry, the source branch
  at `afc9e13ca220edf23fd25724175e827ba9a31b62`, and aggregate D2b receipt
  `345ecdcd3d7c4696d20ebc6eff693eea9c527cb7d60872905da5f31abd557be7`.
  The generator rejects source, aggregate, geometry, CC-identity, and bundled-byte
  drift. Independent literature readback also found and fixed a material Table 4
  transcription bug: R(1) CH3O is `alpha=3.14e10 s^-1`; `3146` is R(2)'s gamma,
  not R(1)'s prefactor. The tracked anchor fixture binds the downloaded paper
  SHA-256. Focused Ruff and 64 tests pass. No full-IRC/Hessian compute or D3b
  table ran; the existing endpoint-only checkpoints are inputs, not SCT evidence.
- 2026-09-07 14:35 PDT (hermes-custom-build-001; profile=workstation) —
  CANONICAL PUBLICATION HARDENING CHECKPOINT, REVIEW-RED: exact pushed head
  `fd7f7bcc02598f11a20cd20ec8e82ee6837891cc` now replays persisted native-Hessian
  TS evidence, validates IRC tangents in the rigid-horizontal molecular quotient,
  skips only fully validated cached path Hessians, enforces fresh path-energy
  reproduction, and reconstructs public typed paths from canonical direction
  receipts instead of accepting a caller trace. A bound IRC wrapper records the
  exact `0.05 A` / 200-step / `0.05` outer-`fmax` / `0.01` inner-`fmax` contract
  and shared run identity. Whole-QM verification is `997 passed, 1 skipped` with
  Ruff/format/diff clean. Two independent cold reviews still prohibit production:
  typed-path acceptance does not yet require the shared execution receipt; final
  IRC publication has a receipt/ancestry race; trusted frozen endpoint fingerprints
  and current code/dependency identity are not re-bound at execution; native-Hessian
  execution remains injectable/self-attested; and an in-process crash during the
  expensive two-direction Sella run loses partial work. No production calculator,
  literature benchmark, D2b gate flip, or D3b table ran.
- 2026-09-07 15:59 PDT (hermes-custom-build-001; profile=workstation) —
  SHARED IRC ANCESTRY CHECKPOINT: exact pushed head `784de7e` makes the shared
  `irc-execution` receipt mandatory, requires its embedded directions to equal both
  standalone receipts, binds all five preflight/TS/execution/direction hashes into
  immutable typed-path-v3 and path-Hessian ancestry, and rereads the complete chain
  after final IRC publication so delete/replace/divergence races fail closed without
  deleting foreign artifacts. Adversarial delete, nested-divergence, final-race,
  legacy-v2, resume, and Hessian-precommit regressions pass; independent spec review
  passed and code-quality review approved after the v3 correction. Whole-QM gates are
  `1005 passed, 1 skipped` with Ruff/format/diff clean. Production remains prohibited
  pending current-environment/trusted-endpoint rebinding, non-injectable native-Hessian
  provenance, and durable mid-IRC restart; no production calculator, literature
  benchmark, D2b gate flip, or D3b table ran.
- 2026-09-07 19:40 PDT (hermes-custom-build-001; profile=workstation) —
  PRODUCTION RESTART HARDENING CHECKPOINT, QUALITY-REVIEW RED: exact pushed head
  `db7aa17c47a5ec740d57376eec50a0bfd2abb9d1` rebinds code/dependency and trusted
  frozen-endpoint identity at public production/resume boundaries, makes native TS/IRC/
  path-Hessian backends internal-only, removes public failure-injection seams, and adds
  atomic partial/full IRC restart receipts bound to one serialized Sella TS initialization.
  Missing, corrupt, replaced, or divergent initialization/direction evidence now fails
  closed before computation or overwrite; a valid restart avoids a second TS
  diagonalization. Independent spec review passed after three correction rounds, and
  whole-QM verification is `1031 passed, 1 skipped` with Ruff/format/diff clean. Final
  code-quality review still prohibits production until the next continuation: define a
  strict bounded PES-cache schema (or reconstruct caches), bind/validate the mass-weighted
  kick against H0/masses/dx, load immutable trusted TS/endpoints without a validation/load
  race and compare resumed qualification fingerprints, compare restart directions with
  final canonical IRC directions in the shared validator, bind executable imported-module
  origins/hashes rather than package versions alone, and hold the route claim across native
  TS Hessian evaluation/publication. No production calculator, literature benchmark,
  D2b gate flip, D3b table, or PR exists.
- 2026-09-08 00:41 PDT (hermes-custom-build-001; profile=workstation) —
  PRODUCTION HARDENING REVIEW PASSED: exact pushed implementation head
  `cab53c44564fd50c2d05f45d153b5975c26759dc` closes every previously identified
  restart, trusted-input, executable-provenance, and route-claim defect. The strict
  final pass binds bounded Sella PES-cache/H0/mass/kick state; canonical restart and
  final directions; immutable trusted geometries; resident Python module code plus
  literal/container state; concrete NumPy/SciPy/PySCF native payloads; and a
  replacement-resistant abstract-kernel route fence with pinned rollback of displaced
  publications. The independent cold review approved all final regressions; whole-QM
  verification is `1075 passed, 1 skipped` with Ruff/format/diff clean. No production
  calculator, 12–20 K literature benchmark, D2b gate flip, D3b table, or PR ran in
  this implementation slice. Acceptance remains red only on executing the bounded
  production campaign and adjudicating its receipt-backed scientific outputs.
- 2026-09-08 01:26 PDT (hermes-custom-build-001; profile=workstation) —
  REAL PREFLIGHT DEFECT CLOSED: the first live dry preflight correctly refused the
  branch-local CPU-only environment because GPU4PySCF was absent; the retry under the
  frozen GPU environment then exposed an import-order false refusal where attestation
  of an earlier module lazily loaded `scipy._lib._util` after execution capture had
  closed. The manifest now preimports the complete declared inventory before any
  attestation can trigger lazy imports. The new regression demonstrably fails against
  the old ordering and passes with the fix; all 180 D2c campaign tests plus focused
  Ruff/format/diff gates pass. No calculator ran and no scientific result is claimed;
  the real receipt will be created only from the pushed corrected head.
- 2026-09-08 01:30 PDT (hermes-custom-build-001; profile=workstation) —
  PREFLIGHT CAPTURE ORDER HARDENED: the first import-order fix exposed the earlier
  trigger: GPU dependency inventory itself imports SciPy before manifest capture.
  Dependency discovery now runs under the same module-execution witness and restores
  the prior profiler in `finally`, while the full inventory remains preimported before
  attestation. A second regression demonstrably fails without this capture scope; all
  181 D2c tests and focused Ruff/format/diff gates pass. No calculator ran; the next
  live dry preflight is pinned to the corrected pushed head.
- 2026-09-08 01:36 PDT (hermes-custom-build-001; profile=workstation) —
  SCIPY REGISTRY ATTESTATION CLOSED: the next live preflight reached source/runtime
  state comparison and found SciPy's `IntegratorBase.integrator_classes` is declared
  as an empty list then populated during the same verified module execution. The
  attestor now permits only a unique ordered registry of classes declared in that
  exact source module, binds their qualified identities, and still rejects any
  injected foreign class. The sabotage regression fails without the rule; all 182
  D2c tests plus focused Ruff/format/diff pass. No calculator ran; receipt generation
  remains pinned to the next corrected pushed head.
- 2026-09-08 01:49 PDT (hermes-custom-build-001; profile=workstation) —
  BOUNDED PRODUCTION FOUNDATION DRIVER READY: production mode now resumes each
  canonical route through fresh TS qualification, shared-checkpoint Sella IRC,
  typed-path publication, and per-point native Hessians; it writes atomic mutable
  phase status plus an immutable completed/failed terminal receipt. SIGTERM/SIGINT
  produce an orderly receipt, and `--finalize-if-running` provides the systemd
  dead-man path if the worker exits without one. Deterministic route-order, failure,
  and dead-man tests pass; whole-QM verification is `1079 passed, 1 skipped` with
  whole-tree Ruff/format/diff clean. No production calculator has run from this
  unpushed source; the fresh preflight and transient launch remain pinned to the
  committed/pushed head.
- 2026-09-08 01:53 PDT (hermes-custom-build-001; profile=workstation) —
  LAUNCH-ARGUMENT DEFECT CLOSED: the first executable preflight from the production
  driver head failed before receipt creation because the complete parser rejected
  `--gpu`/`--gpu-mem-gb` after the early etiquette bootstrap had already consumed
  them. Both bootstrap arguments are now accepted by the complete parser and covered
  by a regression; whole-QM verification is `1080 passed, 1 skipped` with whole-tree
  Ruff/format/diff clean. The rejected invocation created no campaign receipt; a
  fresh run root remains required after this correction is committed and pushed.
- 2026-09-08 01:57 PDT (hermes-custom-build-001; profile=workstation) —
  CANONICAL ROUTE-ORDER DEFECT CLOSED: the first transient unit correctly started
  under the 48-hour bound, but production rejected the preflight because canonical
  JSON serialization sorts route-object keys while the driver compared insertion
  order. The immutable dead-man receipt captured that pre-calculator failure at
  `/mnt/data/vsletten/dissertation-data/task300-d2c-instanton-e901275/production-terminal.json`.
  Production now validates exact route-set equality and executes the source-defined
  canonical order; the regression intentionally supplies reversed receipt order and
  all 184 focused D2c tests plus Ruff/format/diff pass. A fresh pushed head and fresh
  preflight/run root remain required; no electronic-structure calculation ran.
- 2026-09-08 02:09 PDT (hermes-custom-build-001; profile=workstation) —
  PRODUCTION BOUNDARY CAPTURE CLOSED: the corrected route launch reached the first
  production boundary, where a fresh process imported `scipy._lib._util` during
  dependency inventory after the module-execution profiler had already been restored.
  Its immutable terminal receipt is preserved under
  `/mnt/data/vsletten/dissertation-data/task300-d2c-instanton-37e1199/` and proves the
  calculator was never entered. Live boundary identity now captures dependency imports
  under the same profiler as preflight and restores the prior profiler in `finally`;
  the new isolated-process regression exercises that exact seam. Whole-QM verification
  is `1081 passed, 1 skipped` with Ruff/format/diff clean. A fresh pushed head and
  preflight/run root remain required before relaunch.
- 2026-09-08 02:11 PDT (hermes-custom-build-001; profile=workstation) —
  FIRST CALCULATOR RECEIPT, SCIENTIFIC GATE RED: exact pushed head `7dece5a` passed
  live preflight and production provenance, entered the first H+CO TS Hessian, hit
  GPU4PySCF's known contiguity assertion, retried on CPU, then fail-closed because
  the observed analytic Hessian exceeded the canonical symmetry tolerance. The
  immutable receipt is
  `/mnt/data/vsletten/dissertation-data/task300-d2c-instanton-7dece5a/production-terminal.json`;
  no TS qualification or downstream artifact was published. The launch also exposed
  that this driver had not supplied a `gpu_owner` to QI2 despite `--gpu`; the next
  pushed head acquires the process-bound `D2c-instanton-tier` lease with a 49-hour
  diagnostic TTL. Next action is a pushed-head quantitative asymmetry diagnostic,
  not blind tolerance inflation or another identical launch.
- 2026-09-08 02:27 PDT (hermes-custom-build-001; profile=workstation) —
  QUANTITATIVE DIAGNOSTIC COMPLETE, TOLERANCE CHANGE REJECTED: pushed head
  `0b7471a` reran the same bounded first-route evaluation under a verified QI2 lease
  (owner `D2c-instanton-tier`, PID 784938, 18 GB) and measured CPU PySCF canonical
  asymmetry `1.52624800126e-06 Eh/Bohr^2` against matrix max
  `1.21773596401`; the existing absolute tolerance was `2.43547192802e-09`.
  Receipt:
  `/mnt/data/vsletten/dissertation-data/task300-d2c-instanton-0b7471a/production-terminal.json`.
  A cold scientific review rejected simply raising the tolerance: one scalar run
  does not distinguish grid/CPSCF convergence noise from a material defect, the
  proposed `max(1, scale)` rule was not uniformly relative, and its tests left a
  broad false-green interval. The uncommitted relaxation was fully reverted. Exact
  head remains clean; whole-QM verification is `1081 passed, 1 skipped` with
  Ruff/format/diff clean. Next continuation must compare the raw asymmetry and
  post-symmetrization vibrational verdict under a denser DFT grid and stricter CPSCF,
  then encode acceptance/rejection boundaries (including low-scale and direct-result
  construction) before changing policy. No TS qualification, IRC, SCT rate, D2b gate
  flip, D3b table, or PR exists.
- 2026-09-08 04:28 PDT (hermes-custom-build-001; profile=workstation) —
  HESSIAN CONVERGENCE DIAGNOSTIC LAUNCHED: exact pushed source head
  `5f74662a079129b54cbbbe7a37aaba47c6cf6a99` adds a provenance-bound,
  evidence-only four-case raw Hessian/component diagnostic (baseline, strict
  response, dense grid, dense+strict reference), keeps production acceptance
  fail-closed at the historical near-zero absolute floor plus a global spectral
  relative gate, and binds the operative PWB6K/D3BJ Python dispatch plus
  `libs-dftd3.so`. Cold review exposed and the implementation closed missing D3
  dispatch provenance, false-running failure status, and absent terminal/dead-man
  receipts; focused verification is `204 passed` with whole-QM Ruff/format/diff
  clean (the preceding whole suite was `1083 passed`). Fresh same-entrypoint
  preflight identity `07964984561e5c6e6e553f18f70223fac8e999e6d813fb8db48956087946e967`
  is running as bounded user unit `task300-d2c-hessian-5f74662.service`
  (`RuntimeMaxSec=12h`, 16 threads, nice 10, SIGINT + ExecStopPost dead-man).
  Await atomic terminal receipt at
  `/mnt/data/vsletten/dissertation-data/task300-d2c-hessian-5f74662-diagnostic/terminal.json`;
  no TS qualification, IRC, SCT rate, D2b gate flip, D3b table, or PR exists.
- 2026-09-08 04:45 PDT (hermes-custom-build-001; profile=workstation) —
  DIAGNOSTIC INDEPENDENTLY VERIFIED, PRODUCTION STILL RED: bounded unit
  `task300-d2c-hessian-5f74662.service` completed successfully; its diagnostic-only
  receipt contains four cases and explicitly publishes no accepted campaign result.
  A different worker rehashed all 20 raw `18x18` artifacts, reproduced component
  sums/symmetrizations and every projected eigenvalue within `4.44e-16`, and verified
  the complete source/preflight/terminal chain. Grid level 3 retains two negative
  modes (spectator about `-13.5 cm^-1`) even under strict SCF/CPSCF; grid level 5
  produces one negative mode and a stable low positive mode (`14.719` versus
  `14.928 cm^-1`, overlap `0.9999978`, symmetric-matrix relative delta `4.95e-7`).
  Raw electronic skew worsens from `1.52625e-6` to `6.13389e-5 Eh/Bohr^2` on the
  denser grid, so the evidence rejects global tolerance inflation even though dense
  post-symmetrization TS gates pass. Next bounded pass must predeclare tolerances and
  run one dense-grid finite-difference-of-gradients Hessian at this frozen TS; only a
  confirmed one-negative-mode spectrum may authorize a narrowly scoped
  PySCF/PWB6K/D3BJ/grid-5 policy and TS rerun. No IRC, SCT/literature verdict, D2b
  gate flip, D3b table, or PR exists.
- 2026-09-08 15:26 PDT (hermes-custom-build-001; profile=workstation) —
  HALF-STEP FD RECEIPT INDEPENDENTLY VERIFIED, SCIENTIFIC CONFIRMATION REJECTED:
  exact clean/pushed head `83d372d778431f2e733fec8b8ce64e5a5ae19b99`
  adds the immutable `+/-0.005 Bohr` extension of the prior `+/-0.01 Bohr`
  finite-difference Hessian, with fresh-root/no-clobber output claims, strict
  resume/dead-man validation, 36 new receipt-bound gradients, Richardson
  reconstruction, and predeclared matrix/mode gates. Spec review passed after
  three corrective rounds; quality review passed after output-confinement,
  source-root, mode-recovery, dangling-artifact, TOCTOU, and wrong-mode terminal
  corrections. Parent verification is 152 focused tests plus Ruff/format/diff;
  the implementation pass also reported the whole QM suite green (`1233 passed,
  2 skipped`). The first immutable attempt failed before its first calculator
  because this worktree's local venv lacked locked `pyscf-dispersion==1.5.0`;
  that failed root is preserved. After installing the locked dependency, bounded
  unit `task300-d2c-halfstep-83d372d-r2.service` completed all 36 points under a
  2-hour cap, 16 threads, nice 10, SIGINT, and dead-man finalizer. Independent
  read-only verification passed every artifact/provenance check, reconstructed
  both finite-difference matrices, Richardson matrix, all four projected spectra,
  476 gate booleans, and matched receipt SHA-256
  `dbde8c707bddbb8f1d4a226755d4e52efbd15a437bd612255bc0e9c969d5f84a`.
  Exactly ten gates fail, including Richardson-vs-inner spectral-relative drift
  (`3.41075e-5 > 2e-5`), FD-vs-analytic absolute/relative deltas
  (`5.86602e-4 > 2e-4`; `2.62359e-4 > 1e-4`), low-positive eigenvalue floors,
  and low-positive frequency agreement (`6.7622 cm^-1 > 1.0` versus analytic;
  `2.04659 cm^-1 > 0.5` versus the inner stencil). The receipt therefore
  correctly declares `confirmation_passed=false` and
  `accepted_campaign_result=false`. No tolerance change, TS qualification, IRC,
  SCT/literature verdict, D2b gate flip, D3b table, or PR is authorized; the next
  bounded continuation must choose and predeclare a scientifically different
  route rather than replaying this rejected Hessian policy.
