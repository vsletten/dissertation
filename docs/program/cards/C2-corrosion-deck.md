# C2-corrosion-deck — minimal viable pitting-corrosion model

- status: blocked (B3)
- track: C (corrosion)
- priority: P2
- machine: any (ensembles are CPU; big sweeps prefer workstation)
- depends: B3
- claimed-by:

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
- 2026-08-18 — card created (fable).
