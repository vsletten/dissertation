# B4-ensembles-observables — replicas + declarative observables

- status: active
- track: B (platform)
- priority: P1
- machine: any
- depends: B2
- claimed-by: hermes-laptop

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
- 2026-08-21T23:57:01+00:00 — hermes-laptop — Claimed atomic branch `agents/B4-ensembles-observables`. Implemented the `petra-observables` framework, Rayon replica parallelism with ordered deterministic collection, full distributions/means/bootstrap CIs, declarative state-count/event-rate/rate-spectrum/cluster/surface-area sampling, CLI long-form outputs, shipped Ising/Kossel declarations, and A5 exposure-age attachment design. Focused schema/compatibility gates are green; full workspace verification remains.
- 2026-08-21 — omnibus supervisor — unblocked after B2 merged in PR #33.
- 2026-08-18 — card created (fable).
