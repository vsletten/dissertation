# B4-ensembles-observables — replicas + declarative observables

- status: ready
- track: B (platform)
- priority: P1
- machine: any
- depends: B2
- claimed-by:

## Objective
[observables] and ensemble execution per RFC-001: multi-seed parallel
replicas (rayon), aggregation (means/distributions/bootstrap CIs), and
the first observable implementations chosen to serve A5-Phase-0:
event-rate time series, RATE SPECTRA (distribution of instantaneous
site rates — the Fischer/Luttge observable), state counts, cluster-size
distribution, and surface-area accounting (geometric + BET-like
site-count proxy). Design doc note in the PR: how A5's exposure-age
tracking will attach.

## Acceptance
- Ising Binder-cumulant gate (B3) reimplemented ON this framework
  (dogfood); rate-spectra output validated on the kossel etch-pit deck
  (spectrum broadens as pits nucleate — qualitative gate + golden file).

## Progress
- 2026-08-21 — omnibus supervisor — unblocked after B2 merged in PR #33.
- 2026-08-18 — card created (fable).
