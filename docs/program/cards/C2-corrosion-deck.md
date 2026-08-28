# C2-corrosion-deck — minimal viable pitting-corrosion model

- status: done
- track: C (corrosion)
- priority: P2
- machine: any (ensembles are CPU; big sweeps prefer workstation)
- depends: B3
- claimed-by: hermes-8d9e30fee91548

## Objective
Implement the minimal viable model from
docs/scoping/pitting-corrosion.md: film-patch lattice (intact →
Cl-adsorbed → thinned → bare → dissolving/repassivated), chloride
activity + electrode potential via [thermo] activity/mu, autocatalysis
via by_count, defect seeding via [[defects]]. Ensemble target:
pit-initiation-time distributions with Weibull shape vs Shibata's data;
metastable transients as event-rate series. CTMC strategy (this
showcase proves the ORIGINAL strategy in a new domain).

## Acceptance
- Deck + ensemble study reproducing at least the qualitative
  statistics targets in the scoping doc §"phased plan" phase 2;
  results doc in docs/program/results/.

## Progress
- 2026-08-27 — (hermes-8d9e30fee91548; profile=laptop) — implemented the
  24×24 passive-film CTMC deck, deterministic 7×32-replica potential study,
  three-size finite-size scan, current-transient exporter, and four acceptance
  regressions. The study spans isolated metastable patches (0/32 initiated at
  μ=-0.70 kcal/mol) through cooperative clusters (32/32 at +0.80; mean largest
  cluster 7.28), with explicitly descriptive uncensored Weibull-plot slopes,
  pit-number/spatial statistics, exact deadlines, and repair/dissolution pulses.
- 2026-08-22 — omnibus supervisor — unblocked after B3 merged in PR #60;
  PLAN.md already listed C2 as ready, so this reconciles the card metadata.
- 2026-08-18 — card created (fable).

## Result

Delivered `petra/decks/passive-metal.toml`, the bounded
`corrosion_study` executable, generated per-replica/statistical CSV artifacts,
and `docs/program/results/C2-corrosion-deck.md`. The model uses separate chloride
activity and electrochemical μ controls, nonlinear `by_count` autocatalysis,
strain-coupled defect nucleation, CTMC first-passage times, and event-rate current
series. The result document states the coarse-grained/rate-compressed scope and
makes no numerical experimental-calibration claim. `cargo test -p petra-deck
--test corrosion -- --nocapture` passes 4/4, including deterministic replay,
defect/thermo wiring, the low/high clustering transition, and metastable
repassivation drops.
