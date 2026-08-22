# E2-muscovite-full-mechanism — species, defects, delamination, spectra

- status: done
- track: E (muscovite / the perfect circle)
- priority: P1
- machine: any (CPU-light KMC; no barrier campaign in this card)
- depends: E1-muscovite-deck, B5-execution-schedule
- claimed-by: hermes-laptop

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
- 2026-08-22 — hermes-laptop — delivered the schema-v2 full mechanism and
  exact five-segment replay. Extended zones, locally dehydroxylation-driven
  interface delamination, three Ar isotopes, recoil ³⁹Ar traps, cumulative
  release, age spectrum, contamination ratio, and explicitly labelled proxy
  sensitivity bands are all exercised. DEVIATION: the first comparison uses a
  compact 768-site deterministic mechanism demonstration rather than a
  grain-size sweep; its small isotope counts are not presented as statistics
  or a fit, and all phase-3 barriers remain proxy-bracketed.
- 2026-08-20 — created from the verified E1 gap ledger. Blocked on
  B5-execution-schedule; E1 lands the validated fixed-temperature base in the
  same PR that adds this card.

## Result

- Deck: `petra/decks/muscovite-full-mechanism.toml`, deterministically rendered
  by `petra/scripts/build_muscovite_full_deck.py`.
- Analysis: `petra/scripts/muscovite_full_analysis.py`; tracked CSV/JSON/SVG and
  the candid 1998 comparison are under
  `docs/program/results/E2-muscovite-full-mechanism*`.
- Seed-1998 schedule replay is byte-identical across all segment boundaries.
  The 500/550 °C steps release 12 ⁴⁰Ar and 8 ³⁹Ar at a synthetic 26.86 Ma;
  600/650 °C then release 29 trap-delayed ³⁹Ar. The only ³⁶Ar leaves in the
  first step. This is the predicted differential-reservoir direction, not a
  quantitative fit.
- E1's independent 500/700 °C full trajectories still pass unchanged
  (rise 2.13×/13.55×; fall 275.5×/32.3×). Focused Python tests, Petra workspace
  tests, format, Clippy, and paranoid schedule replay are the merge gates.
