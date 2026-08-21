# B3-conformance-decks — Conway, Ising, SIR with analytic gates

- status: active
- track: B (platform)
- priority: P1
- machine: any
- depends: B2
- claimed-by: hermes-laptop

## Objective
Implement SynchronousCA, AsyncMetropolis, and DiscreteTimePCA
strategies per RFC-001, each proven by a tiny conformance deck with an
analytic answer, kept as CI gates: Conway glider (period 4,
displacement (1,1)); 2-D Ising critical temperature via Binder-cumulant
crossing vs Onsager 2.269 (within stated tolerance); SIR epidemic
threshold vs mean-field. This is the "petra is general now" milestone.

## Acceptance
- Three decks in petra/examples/conformance/, three gates in cargo
  test, CTMC parity gates still green.

## Progress
- 2026-08-21 — hermes-laptop — atomically claimed `agents/B3-conformance-decks` from current `origin/main`; began grid compiler and three-strategy conformance implementation under RFC-001 with strict seeded-determinism and CTMC parity gates.
- 2026-08-21 — omnibus supervisor — unblocked after B2 merged in PR #33.
- 2026-08-18 — card created (fable).
