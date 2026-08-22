# D3b-co-hydrogenation-deck — surface-valid methanol ladder

- status: blocked (D2b)
- track: D (astrochemistry)
- priority: P2
- machine: any
- depends: D3, D2b
- claimed-by:

## Objective

Extend D3's verified open-system ice-mantle machinery from H/H2 to the
CO -> HCO -> H2CO -> CH3O/CH2OH -> CH3OH ladder, including validated
addition/abstraction competition, only after D2b produces surface-valid rates.
Compare CO/H2CO/CH3OH yields and temperature trends against the Simons 2020
microscopic-KMC figures and Watanabe/Fuchs experiments.

## Acceptance

- Deck consumes D2b's provenance-carrying surface rate table; no gas-phase rate
  is silently relabeled as Langmuir–Hinshelwood.
- Deterministic ensemble outputs cover 8, 12, 16, and 20 K with documented
  comparisons to Simons 2020 and experimental yield trends.
- CO, H2CO, and CH3OH population/fluence products and uncertainty intervals are
  committed with bounded reproducibility commands.

## Progress

- 2026-08-22 — split from D3 because D2a's merged verdict is an explicit NO-GO
  for its gas-phase rate protocol on surfaces; blocked until D2b passes.
