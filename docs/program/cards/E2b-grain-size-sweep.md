# E2b-grain-size-sweep — statistics for the spectra

- status: ready
- track: E (muscovite / the perfect circle)
- priority: P1
- machine: any (CPU ensembles)
- depends: E2-muscovite-full-mechanism, B4-ensembles-observables
- claimed-by:

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
