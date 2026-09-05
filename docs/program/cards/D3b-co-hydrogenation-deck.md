# D3b-co-hydrogenation-deck — surface-valid methanol ladder

- status: blocked (D2c)
- track: D (astrochemistry)
- priority: P2
- machine: any
- depends: D3, D2c
- claimed-by:

## Objective

Extend D3's verified open-system ice-mantle machinery from H/H2 to the
CO -> HCO -> H2CO -> CH3O/CH2OH -> CH3OH ladder, including validated
addition/abstraction competition, only after D2c emits the D2b-valid
surface-rate table after the D2c gate passes.
Compare CO/H2CO/CH3OH yields and temperature trends against the Simons 2020
microscopic-KMC figures and Watanabe/Fuchs experiments.

## Acceptance

- Deck consumes the D2b-valid provenance-carrying surface-rate table that D2c
  emits after the D2c gate passes; no gas-phase rate is silently relabeled as
  Langmuir–Hinshelwood.
- Deterministic ensemble outputs cover 8, 12, 16, and 20 K with documented
  comparisons to Simons 2020 and experimental yield trends.
- CO, H2CO, and CH3OH population/fluence products and uncertainty intervals are
  committed with bounded reproducibility commands.

## Progress

- 2026-08-29 — fable — re-blocked on D2c: D2b completed but its gate
  verdict was NO-GO, so the rate table this card consumes was withheld
  (docs/program/results/D2b-explicit-surface-rates.md). D2c-instanton-tier
  emits the table when its gate passes. This frontmatter had gone stale
  against PLAN.md and the board feeder dispatched a premature pointer
  (mission-control TASK-243, since parked); one source of truth is the
  card file — keep it current in the same PR that changes the board.
- 2026-08-22 — split from D3 because D2a's merged verdict is an explicit NO-GO
  for its gas-phase rate protocol on surfaces; blocked until D2b passes.
