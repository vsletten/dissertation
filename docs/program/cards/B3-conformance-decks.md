# B3-conformance-decks — Conway, Ising, SIR with analytic gates

- status: blocked (B2)
- track: B (platform)
- priority: P1
- machine: any
- depends: B2
- claimed-by:

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
- 2026-08-18 — card created (fable).
