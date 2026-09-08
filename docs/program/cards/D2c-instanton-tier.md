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
