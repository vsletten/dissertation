# E2b-grain-size-sweep — statistics for the spectra

- status: done
- track: E (muscovite / the perfect circle)
- priority: P1
- machine: any (CPU ensembles)
- depends: E2-muscovite-full-mechanism, B4-ensembles-observables
- claimed-by: hermes-8d9e30fee91548

## Objective

Close E2's recorded DEVIATION: the first full-mechanism comparison ran a
compact 768-site deterministic demonstration whose small isotope counts are
not statistics. Sweep lattice size toward grain-relevant volumes and run
B4-style deterministic replica ensembles at each size, so synthetic age
spectra and release curves carry real distributions and bootstrap CIs.

1. Size ladder from the 4×4×6 demonstration cell up to the largest volume
   the engine sustains comfortably (document the practical ceiling and
   events/sec scaling as an engine datum).
2. At each size: replica ensemble over recoil/defect seeds; emit per-step
   release distributions, age-spectrum bands, and ³⁶Ar/⁴⁰Ar ratio evolution.
3. Test the grain-size claims from the scoping doc §5 — does the staircase
   shape sharpen/persist with volume, and at what size do the E2 features
   stabilize?

## Constraints

- Proxy barriers stay labelled as such (E3a upgrades them; do not wait on
  it — this card is about statistics, not kinetics).
- Deterministic replay must hold at every size; exact-seed regression per
  E2's gates.
- machine: any — laptop-class workers welcome; big sweeps may hand the top
  sizes to the workstation as a follow-up Progress note rather than
  blocking.

## Acceptance

- Results doc with the size ladder, ensemble bands on the synthetic age
  spectrum, and an explicit statement of the volume at which E2's
  qualitative features stabilize (or that they don't, and how).
- Sweep driver + regression landed; Petra tests, replay gates, lint green.
- Card/PLAN/STATUS bookkeeping in the PR.

## Progress

- 2026-08-27 — filed by Fable to close E2's deviation and feed the paper-1
  figure set.
- 2026-08-28 — (hermes-8d9e30fee91548; profile=laptop) — added a
  deterministic size-deck renderer and bounded ensemble driver; ran seven
  32-replica volumes from 768 through 559,872 sites with byte-identical
  exact-seed replay at every size. The `48×48×72` next rung was terminated by
  `SIGKILL`, establishing the evidence-backed host ceiling at 559,872 sites.
  The 5%-release/10%-age adjacent-size gate, strengthened to require ≤5%
  change in defined-age replica support on material (≥0.1% ³⁹Ar release)
  steps, first passes at 49,152 sites and remains green through 165,888 and
  559,872 sites.

## Result

- Verdict: the E2 old-first-step staircase is stable from 49,152 sites upward;
  at 559,872 sites the synthetic apparent-age means are 38.279, 21.315, 0.402,
  0, and 0 Ma over the five heating steps.
- Statistics: 32 replicas per size, deterministic 2,000-resample bootstrap
  bands, explicit defined/missing replica counts for conditional age and ratio
  distributions, and 9,689,761 primary-ensemble events.
- Practical ceiling: 559,872 sites completed 6,740,495 events in 138.420 s;
  1,327,104 sites did not produce output before `SIGKILL`.
- Products: `docs/program/results/E2b-grain-size-sweep.md` and
  `docs/program/results/E2b-grain-size-sweep/{spectrum-bands.csv,size-scaling.csv,stability.json}`.
- Verification: 116 Rust tests and 26 Python tests pass; Cargo fmt, Clippy
  `-D warnings`, Ruff check, and Ruff format are clean. Exact-seed replay is
  byte-identical at all seven completed sizes.
