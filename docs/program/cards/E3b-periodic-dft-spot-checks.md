# E3b-periodic-dft-spot-checks — route 2 calibration of the classical set

- status: done
- track: E (muscovite / the perfect circle)
- priority: P2
- machine: workstation (periodic DFT; CP2K)
- depends: E3a-classical-neb-barriers
- claimed-by: hermes-custom-build-001

## Objective

Phase 3 route 2 of `docs/scoping/ar-muscovite.md` §6: periodic DFT
(CP2K, PBE-D3, ~2×2×1 supercells, climbing-image NEB) spot checks on the
subset of E3a barriers where classical potentials are least trustworthy —
above all the post-dehydroxylation 5-coordinate Al environments in the
dehydroxylate lattice. The deliverable is a calibration statement: for each
checked barrier, DFT vs classical with a correction or a documented
endorsement, so E4 runs on numbers whose error bars mean something.

**POLICY §12 scope amendment (2026-09-14):** the three fully relaxed endpoint
and CI-NEB calculations exceed the one-workstation platform-test envelope. This
card therefore banks three PBE-D3 **frozen-classical-path energy profiles** on
the immutable E3a images, paired with fresh matched-cell classical profiles.
The reported quantity is the energy rise along each tested path, not a
DFT-relaxed activation barrier or converged DFT minimum-energy path. Full
periodic CI-NEB is outside this phase's acceptance and is not claimed.

## Constraints

- Follow PROTOCOL.md and POLICY.md; QI2 GPU-lease etiquette; heavy
  artifacts archived with hashes.
- Spot checks only — this card does not re-run the full E3a matrix.
- No silent overrides: where DFT disagrees with the classical number, the
  deck fragment records both and states which one E4 uses and why.

## Acceptance

- Three eight-image PBE-D3 frozen-path profiles are complete, including the
  dehydroxylate-lattice Ar hop whose prepared identity gate proves interaction
  with the five-coordinate-Al environment.
- Every image has converged SCF, normal termination, immutable coordinate/input/
  native-output hashes, an exact bounded-controller receipt, and no surviving
  container; incomplete profiles emit no numeric result.
- Fresh matched-cell classical profiles are complete for the same three paths.
  The results doc compares like-for-like profile rises and gives a typed
  `close-enough`, `disagrees`, or `incomplete` survey verdict per route without
  relabeling the frozen-path rise as a DFT barrier or applying it as a correction.
- Tests/lint green; card/PLAN/STATUS bookkeeping in the PR.

## Progress

- 2026-09-14 16:40 PDT — `(hermes-custom-build-001; profile=workstation)` —
  `REVIEW FIX:` the one POLICY §13 cold review found two safety defects and one
  completion defect. The LAMMPS parser now accepts only the row bound to its
  exact declared thermo header, the runner recomputes and verifies the canonical
  atom-map identity, and `e3b-calibration-fragment.toml` records every E4 value
  selection. All three classical profiles were replayed through the corrected
  runner with unchanged numbers; comparison receipt v2 independently verifies
  the fixes. Non-blocking source-attestation hardening is captured by ready card
  `E3b2-classical-receipt-code-binding` per the review cost ceiling.

- 2026-09-14 16:23 PDT — `(hermes-custom-build-001; profile=workstation)` —
  DONE under POLICY §12's platform-test amendment. All three PBE-D3 frozen-path
  profiles and all three fresh matched-cell LAMMPS profiles completed 8/8
  images with hash-bound receipts. Coordinate-matched rises are respectively
  297.913659 vs 9877.633000, 305.946204 vs 9877.629000, and 425.854862 vs
  39166.519000 kcal/mol; all three typed verdicts are `disagrees` at the
  repository's 5.0 kcal/mol transfer tolerance. No correction, endorsement,
  relaxed DFT barrier, CI-NEB, or production result is claimed. Independent
  receipt replay verified all 48 image records and container cleanup; the
  E3-focused suite is 60/60 green after the review fixes.

- 2026-09-14 15:07 PDT — `(hermes-custom-build-001; profile=workstation)` —
  The reconstructed-replication profile completed 8/8 converged PBE-D3 single
  points in 2,248.029 s with a hash-bound frozen-path rise of 305.946204
  kcal/mol; independent replay verified every coordinate/input/output hash, all
  normal-termination/SCF gates, and the receipt arithmetic. Launched the third
  profile, `xenon-divacancy`, from exact clean pushed head
  `d910d5617eece3b12da64c8c3b71df42bfcae39f` in bounded transient unit
  `task320-e3b-xenon-frozen-profile.service` (invocation
  `e0960c9d933b4c1e918d7afc3ad5d621`, `RuntimeMaxSec=14,400`, 4 MPI × 4
  OpenMP, `Nice=10`, `CPUQuota=1600%`, 60 GiB controller / 48 GiB container
  caps). Prelaunch load15 was 2.64 with 40,752,040 KiB available and no swap
  growth; the pinned CP2K container and a growing image-00 native output were
  observed live. Planned atomic receipt:
  `/mnt/data/vsletten/dissertation-data/e3b-periodic-dft-calibration-20260914/receipts/xenon-frozen-profile-01.json`.
  Scientific acceptance is 2/3 pending this receipt, then the three fresh
  matched-cell classical profiles, synthesis, one allowed cold review, and PR.

- 2026-09-14 13:58 PDT — `(hermes-custom-build-001; profile=workstation)` —
  The required dehydroxylate profile completed 8/8 converged PBE-D3 single
  points in 2,412.217 s with a frozen-path rise of 297.913659 kcal/mol; the
  receipt hashes every immutable coordinate, rendered input, and native output
  and proves all containers absent. Launched the second profile,
  `reconstructed-replication`, from verified pushed head
  `e60fe5b427d89c11a300f1d062e746a7349bd720` in bounded transient unit
  `task320-e3b-reconstructed-frozen-profile.service` (invocation
  `08d86444aa3e4ace8248ab86e0970dd2`, `RuntimeMaxSec=14,400`, 4 MPI × 4
  OpenMP, `Nice=10`, `CPUQuota=1600%`, 60 GiB controller / 48 GiB container
  caps). Prelaunch load15 was 2.97 with 40,172,700 KiB available; the pinned
  CP2K container and growing native image-00 output were observed live. Atomic
  receipt:
  `/mnt/data/vsletten/dissertation-data/e3b-periodic-dft-calibration-20260914/receipts/reconstructed-frozen-profile-01.json`.
  Scientific acceptance is 1/3 pending this receipt, the Xe profile, and all
  three fresh matched-cell classical profiles.

- 2026-09-14 12:49 PDT — `(hermes-custom-build-001; profile=workstation)` —
  Launched the first narrowed survey profile, the required
  `dehydroxylate-lattice` five-coordinate-Al route, from pushed harness head
  `b56b1fcd6e381489f39ddf04e9751e5262d8e43d`. Bounded transient unit
  `task320-e3b-dehydroxylate-frozen-profile.service` has invocation
  `14444d8b45e44030844484bc7692290b`, `RuntimeMaxSec=14,400`, 4 MPI × 4
  OpenMP, `Nice=10`, `CPUQuota=1600%`, and a 48 GiB container limit. Runner
  SHA-256 is `6c99aa6a4143f720100300622eeef32bc7768f747b3fcdb26133417fd0345fc3`;
  prelaunch load15 was 5.82 with 39,895,036 KiB available. Native image-00
  stdout and the pinned CP2K container were observed live. Atomic receipt:
  `/mnt/data/vsletten/dissertation-data/e3b-periodic-dft-calibration-20260914/receipts/dehydroxylate-frozen-profile-01.json`.
  Scientific acceptance remains 0/3 until all eight images complete and the
  typed receipt emits a frozen-path profile rise.

- 2026-09-14 12:47 PDT — `(hermes-custom-build-001; profile=workstation)` —
  `DEVIATION:` POLICY §12 superseded the full three-route CI-NEB plan after the
  first dehydroxylate endpoint GEO_OPT used 1,830.0 s for only ten converged SCF/
  geometry cycles and timed out without endpoint convergence. The atomic
  receipt is `dehydroxylate-cp2k-initial-02.json`; it emitted no barrier.
  Narrowed the card to three frozen-classical-path PBE-D3 profiles plus fresh
  matched-cell classical comparisons. Added the fail-closed bounded profile
  runner and six tests; the existing + new focused suite is 36/36 green.

- 2026-09-14 12:04 PDT — `(hermes-custom-build-001; profile=workstation)` —
  Launched the first scientific stage for the required five-coordinate-Al route:
  dehydroxylate-lattice `cp2k_initial` in bounded transient unit
  `task320-e3b-dehydroxylate-initial.service` (invocation
  `6c9b8c4dd3d14d4893e29b65535a0d8f`). The envelope is 4 MPI × 4 OpenMP,
  `Nice=10`, `CPUQuota=1600%`, 48 GiB container memory, a 1,800 s internal
  timeout, and `RuntimeMaxSec=1,920`; atomic completion receipt:
  `/mnt/data/vsletten/dissertation-data/e3b-periodic-dft-calibration-20260914/receipts/dehydroxylate-cp2k-initial-02.json`.
  Prelaunch load15 was 0.87 with 40,608,508 KiB available; the pinned CP2K
  container and growing native stdout were observed live. Runner SHA-256 is
  `0277a43756af7d6f493a94a14726f57b7dc9e9b280f874e6027a223f269fd1e5`.
  A preceding launch with the obsolete model alias `dehydroxylate-screen`
  failed before container construction or compute and emitted no receipt; the
  corrected launch uses the exact prepared model identity.

- 2026-09-14 11:17 PDT — `(hermes-custom-build-001; profile=workstation)` —
  UNBLOCKED by E3b1's completed evidence contract and a converged 334-atom
  reconstructed-Ar timing probe (780.622 s). The corrected native-output gates
  are ready for bounded endpoint/BAND work. Scientific acceptance is unchanged:
  0/3 DFT barriers, no correction/endorsement, and no BAND lower bound yet.

- 2026-09-14 09:01 PDT — `(hermes-custom-build-001; profile=workstation)` —
  Cold review rejected the 3×2 compute-envelope verdict: a neutral wrapped
  2×2×1 route cell exists, and an unfinished force evaluation cannot establish
  a BAND lower bound. The harness now prepares valid 334/331/334-atom models,
  refuses manual observation numbers, and rejects all three incomplete E3a
  classical references. A real 2×2 PBE-D3 force probe remained unfinished at
  the execution channel's 420-second ceiling; no DFT number was emitted and no
  campaign-feasibility conclusion is claimed. `E3b1-periodic-dft-evidence-contract`
  owns the remaining raw-parser, immutable-artifact, classical-reference,
  interaction-criterion, and bounded timing work. Results:
  `docs/program/results/E3b-periodic-dft-spot-checks.md`.

- 2026-08-27 — filed by Fable alongside E3a.

## Result

Completed the six-profile coordinate-matched survey and its fail-closed
runner. `scripts/e3b_classical_frozen_path.py` binds every immutable XYZ row to
the prepared source atom identity, evaluates the unchanged LAMMPS topology at
those coordinates, and emits no numeric profile after a failed image or
manifest drift. The three PBE-D3 profile rises are 297.913659, 305.946204, and
425.854862 kcal/mol; the same-path classical rises are 9877.633000,
9877.629000, and 39166.519000 kcal/mol. All are typed `disagrees`; none is
promoted to a barrier or correction.

Evidence and claim boundaries are in
`docs/program/results/E3b-periodic-dft-spot-checks.md`. Durable receipts live
under `/mnt/data/vsletten/dissertation-data/e3b-periodic-dft-calibration-20260914/`;
the independent comparison receipt is
`receipts/comparison-verification-02.json` at SHA-256
`b0f7b58bd3bf7a79e6d86365f04cea16a57eaeb07f5efc0f0b92e2351e5a8cc9`.

## Prior blocked result (superseded by the 2026-09-14 11:17 PDT unblock)

Blocked on `E3b1-periodic-dft-evidence-contract`, not completed. The original
acceptance remains **0/3 DFT barriers checked**. A valid 2×2 preparation exists,
but the current E3a references are all typed incomplete and manual observation
JSON is no longer accepted as evidence. No compute-envelope verdict, DFT
barrier, correction, endorsement, or deck replacement is claimed.
