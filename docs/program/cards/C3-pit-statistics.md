# C3-pit-statistics — quantitative pitting statistics from the C2 deck

- status: ready
- track: C (corrosion)
- priority: P2
- machine: any (CPU ensembles; macbot-friendly)
- depends: C2-corrosion-deck
- claimed-by:

## Objective

Turn C2's qualitative shape study into quantitative pitting statistics on
the same `petra/decks/passive-metal.toml` deck family:

1. Scale replica counts well past C2's 7×32 (target ≥ 256 replicas per
   potential point) and extend the potential sweep resolution around the
   metastable→stable transition.
2. Fit induction-time distributions properly: Weibull parameter estimates
   with bootstrap CIs per potential, plus the survival/hazard diagnostic
   family Shibata-style analyses use. Report goodness-of-fit honestly —
   rate-compressed model, chemical-potential proxy for volts, so the claim
   is distribution-family agreement, not calibrated prediction (keep C2's
   honesty language).
3. Pit-size/depth distributions and spatial clustering statistics vs
   potential; nucleation-rate scaling; a finite-size-corrected estimate of
   the critical-potential crossover extending C2's three-size scan.

## Constraints

- CPU-only ensembles (B4 machinery); deterministic replays and exact-seed
  regression gates per C2. Laptop-class workers welcome; if a sweep rung
  exceeds host memory, record the ceiling (E2b pattern) rather than
  escalating silently.
- No new mechanism states in this card — deck physics changes are a
  separate card. Rate compression stays; every plot labels it.

## Acceptance

- Results doc with Weibull fits + CIs, induction/survival diagnostics,
  size/clustering distributions, nucleation-rate scaling, and the
  finite-size-corrected crossover, each with the shape-level honesty
  caveat stated.
- Sweep driver + regression tests landed; Petra tests, replay gates,
  Clippy/lint green; card/PLAN/STATUS bookkeeping in the PR.

## Progress

- 2026-08-29 — filed by Fable per Victor's fleet-routing directive: first
  machine:any C-track science card so laptop/cloud workers have corrosion
  work while the workstation runs QM campaigns.
