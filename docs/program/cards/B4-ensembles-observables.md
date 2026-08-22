# B4-ensembles-observables — replicas + declarative observables

- status: done
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
- 2026-08-22T00:03:06+00:00 — hermes-laptop — Acceptance verified: the full Petra workspace has 88 passing tests; Clippy is warning-free with `-D warnings`; the parallel-framework Binder crossing, deterministic Kossel rate-spectrum broadening golden, single/ensemble CLI output, v1 compatibility, and all bitwise parity gates are green.
- 2026-08-21T23:57:01+00:00 — hermes-laptop — Claimed atomic branch `agents/B4-ensembles-observables`. Implemented the `petra-observables` framework, Rayon replica parallelism with ordered deterministic collection, full distributions/means/bootstrap CIs, declarative state-count/event-rate/rate-spectrum/cluster/surface-area sampling, CLI long-form outputs, shipped Ising/Kossel declarations, and A5 exposure-age attachment design. Focused schema/compatibility gates are green; full workspace verification remains.
- 2026-08-21 — omnibus supervisor — unblocked after B2 merged in PR #33.
- 2026-08-18 — card created (fable).

## Result

Added the `petra-observables` workspace crate and compiled observable declarations.
Replica trajectories run serially and independently while Rayon parallelizes only
across replicas; indexed collection preserves replica order. Aggregation retains
the complete distribution and reports means plus deterministic bootstrap 95% CIs.

The CLI now writes `observables.csv` for both single and ensemble runs, plus
`ensemble.csv` and `ensemble-summary.csv` for ensembles. Implemented declarations
are `state_counts`, `event_rates`, `rate_spectra`, filtered `cluster_sizes`, and
surface-area accounting with separate projected-geometric and BET-like exposed-site
measures. The v1 compatibility shim accepts explicit observables without changing
undeclared v1 decks.

Acceptance evidence:
- Ising Binder crossing runs on the Rayon ensemble + declarative state-count path.
- Kossel etch-pit rate-spectrum log-width broadens from `0.847370502` to
  `1.141527947` after 200 seeded steps and matches the committed golden file.
- Full `cargo test --workspace`: 88 passed, 0 failed.
- `cargo clippy --workspace --all-targets -- -D warnings`: clean.
- Existing CTMC bitwise parity and deterministic draw-order gates remain green.

The A5 exposure-age attachment contract is documented in
`petra/docs/OBSERVABLES.md`: core-owned deterministic exposure timestamps plus a
read-only distribution observable using the same ensemble aggregation machinery.
