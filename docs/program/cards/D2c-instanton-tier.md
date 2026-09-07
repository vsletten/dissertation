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
  note for the production driver. Whole-QM verification is 865 passed, 1 skipped plus
  whole-tree Ruff/format; no IRC/Hessian production compute or D3b table ran. Acceptance
  remains red until the resumable production stages execute and matched-site 12-20 K
  literature/branching gates pass.
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
