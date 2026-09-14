# E3b1-periodic-dft-evidence-contract — close calibration evidence gates

- status: done
- track: E (muscovite / the perfect circle)
- priority: P2
- machine: workstation
- depends: E3a-classical-neb-barriers
- claimed-by:

## Objective

Make E3b's valid 2×2×1 CP2K/LAMMPS calibration route executable without a
false-green result. Close the evidence defects found by E3b's one allowed cold
review before any DFT barrier or correction can be banked.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`.
- Preserve the valid periodic 2×2×1 route crop (wrapped window `(5,2)`) and its
  exact identity/neutrality/topology gates.
- No QM unit may exceed four hours. Any expected >30-minute run uses a bounded
  transient unit with an atomic receipt.
- This card produces trustworthy machinery and prerequisite evidence; it does
  not claim an E3b DFT barrier unless the raw outputs independently prove one.

## Acceptance

- Prepared inputs are immutable and runtime outputs live outside their manifest;
  a smoke run cannot invalidate later `smoke`/`analyze` commands.
- CP2K endpoint/BAND and matched-classical results are parsed from hash-bound raw
  outputs. Hand-entered convergence booleans or barriers cannot emit a
  calibration number.
- The three E3a source barriers used by E3b are either independently converged
  and hash-bound or remain typed incomplete; incomplete references cannot yield
  a complete correction/endorsement.
- The dehydroxylate check defines and tests what “involves five-coordinate Al”
  means chemically; proximity to the route midpoint alone is not sufficient.
- One completed 2×2 reconstructed-Ar `ENERGY_FORCE` timing probe runs under a
  bounded transient unit and records method, input/output hashes, wall time,
  resources, and cleanup. It is planning evidence, not an invented BAND lower
  bound.
- Tests/lint green; result doc, card, PLAN, and STATUS updated in the PR.

## Progress

- 2026-09-14 11:17 PDT — `(hermes-custom-build-001; profile=workstation)` —
  DONE after the single cold review's six must-fix findings were regression-
  closed and final verification passed: 30 focused tests, Ruff check/format,
  compile, and diff checks. A real 334-atom reconstructed-Ar PBE-D3
  `ENERGY_FORCE` probe converged in 780.622 s under a bounded transient unit
  (4 MPI × 4 OpenMP, 48 GiB, 3700 s ceiling); its atomic receipt binds the
  method, exact input/dependencies, raw output, elapsed time, resources, and
  cleanup. Final evidence SHA-256 is
  `cc199ab545fa29e9d6ceb5c24db84f9229f0d0208f12008631dda561f1d721e3`;
  details: `docs/program/results/E3b1-periodic-dft-evidence-contract.md`.
  This is planning/evidence-contract proof, not a DFT BAND barrier.

- 2026-09-14 10:42 PDT — `(hermes-custom-build-001; profile=workstation)` —
  Implemented and locally verified the E3b1 evidence machinery. Prepared inputs
  now execute only from an isolated copied runtime tree; timing receipts bind
  the prepared manifest, all referenced inputs, executable/method identity,
  raw output, elapsed time, declared resources, and cleanup. New parsers copy,
  hash, and reparse raw CP2K endpoint/BAND plus matched LAMMPS endpoint/NEB
  evidence before analysis; hand-entered fields cannot emit a calibration.
  The three E3a anchors independently reparse to `incomplete-convergence` and
  remain non-promotable. The dehydroxylate criterion now proves Al 10/11 changed
  6→5 coordination by removal of hydroxyl O 76 and intersects the explicit
  eight-image noble-gas path, rather than relying on midpoint proximity. Real
  input regeneration preserved wrapped window `(5,2)`, 334/331/334 atoms, and
  all preparation gates; 28 focused tests plus Ruff/format/compile/diff checks
  pass. A bounded reconstructed-Ar `ENERGY_FORCE` timing unit remains to run.

- 2026-09-14 09:01 PDT — `(hermes-custom-build-001; profile=workstation)` —
  Filed from E3b's cold review. E3b already corrected the cell to valid 2×2,
  disabled manual-observation promotion, and rejects incomplete E3a references;
  this card owns the remaining executable evidence path and review follow-ups.

## Result

Done. Native, hash-bound CP2K/LAMMPS parsing; immutable prepared/runtime
separation; incomplete E3a-anchor rejection; the chemical five-coordinate-Al
criterion; and the completed bounded 2×2 reconstructed-Ar timing probe are
verified. E3b is unblocked but remains at 0/3 DFT barriers; no barrier,
correction, endorsement, or BAND lower bound is claimed by this card.
