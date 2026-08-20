# B2-engine-refactor — petra core behind the UpdateStrategy trait

- status: in-progress
- track: B (platform)
- priority: P0
- machine: any
- depends: B1
- claimed-by: hermes workstation worker (agents/B2-engine-refactor)

## Objective
Implement RFC-001: refactor petra-core/petra-deck so the CTMC engine is
one UpdateStrategy implementation behind the trait, deck schema v2
parses (with the v1 compat shim), [execution] selects strategy. NO new
strategies in this card (that's B3) — this card is pure refactor.

## Acceptance
- Bitwise-parity gates: kaolinite.toml and kossel.toml produce
  IDENTICAL trajectories (same seeds) before/after, long runs — same
  discipline as the C++→Rust port. cargo test green incl. new
  schema-v2 round-trip tests; v1 decks load unchanged.

## Progress

- 2026-08-20 — hermes worker 1 (sol) — claimed; goldens captured BEFORE
  refactor; ExactCtmc behind the UpdateStrategy trait; v2 shim;
  splitmix64 vectors; 20k byte-identical parity; then adversarial
  self-review surfaced 7 real findings; ceiling hit with precise
  handoff (incremental-bookkeeping amendment held).
- 2026-08-20 — hermes worker 2 (sol) — continuation from dirty diff with
  supervision-dispositioned findings: all 7 addressed. DEVIATION
  (dispositioned): shipped Kossel 8^3 deck exhausts at 1,512 events, so
  the 20k parity gate uses enlarged dims AND a second gate pins the
  shipped deck's complete natural lifetime byte-identical
  (sha256 687fd5bf...) — full-lifetime parity is the stronger claim.
  Ceiling hit at final mechanical closeout.
- 2026-08-20 — fable — gates rerun on final state (workspace tests,
  clippy -D warnings, both --paranoid decks), final edit committed,
  closeout executed.

## Result

RFC-001 implemented: engine core refactored behind UpdateStrategy with
ExactCtmc as the sole executable strategy (others compile-rejected at
selection); schema v2 parses fully (all RFC surfaces incl. select_mode,
source rules, rekind, per-strategy blocks); v1 decks load via pure
load-time shim with full compiled-table equality tests; determinism
contract enforced (draw-order preserved — parity byte-identical on
kaolinite 20k, kossel 20k enlarged, kossel full lifetime; splitmix64
test vectors green). B3 (conformance strategies) is unblocked.
- 2026-08-20 — verification stage green: `cargo test -j 16 --workspace`
  (54 tests), `cargo build -j 16 --workspace`, and clippy across all targets
  with warnings denied. CLI `--paranoid` completed 5,000 kaolinite steps and
  correctly tracked Kossel through full dissolution (1,512 fired events).
  The two 20k parity logs and all analytic gates remain green.
- 2026-08-20 — deck schema v2 + v1 shim complete. A custom one-pass
  deserializer normalizes every schema-less v1 deck to v2/`ctmc`; explicit
  cell-form v2 decks compile through the same dense tables, mixed/unknown
  schemas fail closed, and future strategies parse but cannot silently run as
  CTMC. All four shipped v1 decks shim+compile. Added the RFC's exact
  splitmix64 replica-seed function/vectors and a v1↔v2 trajectory-equivalence
  gate. Grid-form v2 parses and is deliberately compile-rejected until B3,
  which owns the first grid conformance decks and non-CTMC strategies.
- 2026-08-20 — `UpdateStrategy`, `StepCtx`, `StepOutcome`, and stateless
  `ExactCtmc` seam landed in petra-core. The existing selection/wait/apply/
  dirty-propagation order moved intact behind the trait; `Engine::step` remains
  the compatibility wrapper. Core tests pass and both 20k Fired goldens remain
  byte-identical after the move.
- 2026-08-20 — baseline parity harness complete before engine changes:
  `parity_golden.rs` records full `Fired` values as fixed-width bytes
  (including exact `f64::to_bits()` time) for seed 42 and 20,000 steps.
  Kaolinite uses its shipped dimensions; Kossel uses the shipped deck physics
  on a 20³ lattice because the tutorial's 8³ crystal fully dissolves at step
  1,513. Both checked-in 520,031-byte goldens regenerate and compare cleanly.
- 2026-08-20 — claimed by hermes workstation worker. Read the contract,
  corrected RFC-001 in full, DESIGN.md, and the current petra-core/
  petra-deck implementation. Recovered the accepted RFC review commit
  (`28589de`, including pinned splitmix64 vectors) from the surviving local
  B1 worktree because main only carried the pre-review squash; cherry-picked
  it onto this claim branch unchanged. Baseline parity harness is next,
  before any engine mutation.
- 2026-08-18 — card created (fable).
