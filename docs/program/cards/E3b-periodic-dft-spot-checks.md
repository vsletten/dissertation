# E3b-periodic-dft-spot-checks — route 2 calibration of the classical set

- status: blocked
- track: E (muscovite / the perfect circle)
- priority: P2
- machine: workstation (periodic DFT; CP2K)
- depends: E3a-classical-neb-barriers
- blocked-on: E3a-classical-neb-barriers
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

- 2026-08-27 — filed by Fable alongside E3a.
