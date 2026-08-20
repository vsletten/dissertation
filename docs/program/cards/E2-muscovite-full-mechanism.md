# E2-muscovite-full-mechanism — species, defects, delamination, spectra

- status: blocked
- track: E (muscovite / the perfect circle)
- priority: P1
- machine: any (CPU-light KMC; no barrier campaign in this card)
- depends: E1-muscovite-deck, B5-execution-schedule
- claimed-by:

## Objective
Extend the verified E1 phase-1 gallery into the full mechanism from
`docs/scoping/ar-muscovite.md` §6 phase 2: extended-defect zones,
delamination interfaces, species-resolved ⁴⁰Ar/³⁹Ar/³⁶Ar reservoirs,
recoil-initialized octahedral ³⁹Ar traps, and piecewise-isothermal execution.
Emit cumulative species release, synthetic ⁴⁰Ar/³⁹Ar age spectra, and the
first mechanism comparison against the 1998 staircase/contamination claims.
Use sensitivity brackets for barriers still awaiting the separate phase-3
campaign; do not present proxy values as computed kinetics.

## Acceptance
- Starts from the E1 deck and preserves its 500/700 °C qualitative gates.
- Uses merged B5 schedule support for a documented incremental-heating
  program with deterministic replay across segment boundaries.
- Initial populations distinguish gallery vacancies, extended defects, and
  octahedral ³⁹Ar traps; event output distinguishes ⁴⁰Ar, ³⁹Ar, and ³⁶Ar.
- Delamination is an explicit interface-state transition driven by local
  dehydroxylation, not an unlabelled global rate multiplier.
- Results include synthetic age spectrum, ³⁶Ar/⁴⁰Ar release-ratio evolution,
  sensitivity bands for uncomputed barriers, and direct pointers to the 1998
  target figures/data.
- Petra tests, deterministic replay, and the E1 regression gates are green.

## Progress
- 2026-08-20 — created from the verified E1 gap ledger. Blocked on
  B5-execution-schedule; E1 lands the validated fixed-temperature base in the
  same PR that adds this card.
