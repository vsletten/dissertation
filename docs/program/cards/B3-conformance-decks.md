# B3-conformance-decks — Conway, Ising, SIR with analytic gates

- status: done
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
- 2026-08-21 — hermes-laptop — completed B3. Independent pre-commit review first found batch-trajectory, simultaneous-write, tiny-periodic-grid, ensemble-dispatch, and SIR-evidence defects; each was corrected and the second independent review passed with no blocking findings.
- 2026-08-21 — hermes-laptop — atomically claimed `agents/B3-conformance-decks` from current `origin/main`; began grid compiler and three-strategy conformance implementation under RFC-001 with strict seeded-determinism and CTMC parity gates.
- 2026-08-21 — omnibus supervisor — unblocked after B2 merged in PR #33.
- 2026-08-18 — card created (fable).

## Result

Petra now compiles square/hex/cubic schema-v2 grids into deterministic simple
CSR graphs and executes deck-selected `SynchronousCA`, `AsyncMetropolis`, and
`DiscreteTimePCA` strategies alongside bitwise-frozen `ExactCtmc`. Batch
strategies read one pre-step state and commit atomically; conflicting writes
fail before mutation. Strategy RNG order has pinned trajectory fingerprints,
and tiny periodic grids are loop/duplicate-free.

The three shipped decks under `petra/examples/conformance/` are permanent CI
gates. Conway proves glider period 4 with displacement (1,1), blinker period 2,
and block still-life. The Ising Binder crossing is `2.273230` versus Onsager
`2.269185` (absolute error `0.004045`, tolerance `0.25`). SIR's interpolated
growth threshold is `R0=1.017094` (tolerance `0.15` around 1); its finite-lattice
macroscopic-outbreak probabilities are `0.000000` at `R0=0.5` and `0.476562`
at `R0=2.0` across the deck-declared 128 hash-seeded replicas.

Verification is green: `cargo test -j 16 --workspace` (including all 12 B3
conformance tests and all three 20k/full-lifetime CTMC parity gates),
`cargo clippy --workspace --all-targets -- -D warnings`,
`cargo fmt --all -- --check`, `cargo build -j 16 --workspace`, and
`git diff --check`. CLI execution honors the deck's strategy, replica count,
and seed policy. CTMC trajectory output remains supported; discrete batch
trajectory requests fail closed instead of emitting ambiguous replay rows.
