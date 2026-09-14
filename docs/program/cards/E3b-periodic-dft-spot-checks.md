# E3b-periodic-dft-spot-checks — route 2 calibration of the classical set

- status: ready
- track: E (muscovite / the perfect circle)
- priority: P2
- machine: workstation (periodic DFT; CP2K)
- depends: E3a-classical-neb-barriers
- claimed-by:

## Objective

Phase 3 route 2 of `docs/scoping/ar-muscovite.md` §6: periodic DFT
(CP2K, PBE-D3, ~2×2×1 supercells, climbing-image NEB) spot checks on the
subset of E3a barriers where classical potentials are least trustworthy —
above all the post-dehydroxylation 5-coordinate Al environments in the
dehydroxylate lattice. The deliverable is a calibration statement: for each
checked barrier, DFT vs classical with a correction or a documented
endorsement, so E4 runs on numbers whose error bars mean something.

## Constraints

- Follow PROTOCOL.md and POLICY.md; QI2 GPU-lease etiquette; heavy
  artifacts archived with hashes.
- Spot checks only — this card does not re-run the full E3a matrix.
- No silent overrides: where DFT disagrees with the classical number, the
  deck fragment records both and states which one E4 uses and why.

## Acceptance

- ≥ 3 E3a barriers checked, including at least one dehydroxylate-lattice
  Ar hop involving 5-coordinate Al.
- Results doc extends the E3a table with a DFT column and per-barrier
  calibration verdicts.
- Tests/lint green; card/PLAN/STATUS bookkeeping in the PR.

## Progress

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

## Prior blocked result (superseded by the 2026-09-14 11:17 PDT unblock)

Blocked on `E3b1-periodic-dft-evidence-contract`, not completed. The original
acceptance remains **0/3 DFT barriers checked**. A valid 2×2 preparation exists,
but the current E3a references are all typed incomplete and manual observation
JSON is no longer accepted as evidence. No compute-envelope verdict, DFT
barrier, correction, endorsement, or deck replacement is claimed.
