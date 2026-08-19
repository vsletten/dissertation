# B2-engine-refactor — petra core behind the UpdateStrategy trait

- status: blocked (B1 RFC must merge + get Victor's review ack)
- track: B (platform)
- priority: P0
- machine: any
- depends: B1
- claimed-by:

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
- 2026-08-18 — card created (fable).
