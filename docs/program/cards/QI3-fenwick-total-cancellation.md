# QI3-fenwick-total-cancellation — preserve slow events beside extreme rates

- status: done
- track: B (platform quality infrastructure)
- priority: P1
- machine: any
- depends: —
- claimed-by: hermes-cloud-default

## Objective
Harden Petra's incremental Fenwick event total against catastrophic
cancellation when a deck mixes a short-lived lumped rate near `1e12 s^-1`
with surviving events near `1e-6 s^-1`. E1-muscovite-deck reproduced a state
where `Engine::step()` returned "events exist but all rates are zero" even
though slow dehydroxylation events remained. Preserve the slow events without
weakening deterministic replay or ordinary-run performance.

Start from `.claude/learnings/gotchas.md` (2026-08-20 entry) and the E1 result
doc. Prefer a narrow tree rebuild/consistency repair at the impossible-total
boundary before redesigning the summation structure.

## Acceptance
- A deterministic regression mixes an extreme fast event with surviving slow
  events and proves the engine reaches the true no-event state rather than
  the false zero-total stop.
- The repair re-establishes a positive total from current leaves, or returns a
  truthful terminal reason if every leaf is actually zero.
- Existing Petra tests and `--paranoid` gates remain green.
- Card Progress/Result and `docs/program/STATUS.md` land in the implementation
  PR per `docs/program/PROTOCOL.md`.

## Progress
- 2026-08-20 — created from the E1 muscovite phase-1 numerical-range finding;
  the E1 deck uses a scientifically sufficient 1 s⁻¹ surface rate and is not
  blocked on this hardening.
- 2026-08-20 — hermes-cloud-default reproduced the defect with a deterministic
  two-site deck: removing a one-shot `1e12 s⁻¹` event left a live `1e-6 s⁻¹`
  leaf while the incrementally maintained Fenwick total became zero. Added the
  regression first (observed `Stop::ZeroRate`), then repaired the impossible
  total by rebuilding from authoritative leaves before terminal classification.

## Result
Petra now preserves slow events after catastrophic incremental-total
cancellation without changing ordinary selection or RNG draw order. The repair
is an O(N) leaf rebuild only at the impossible `total <= 0` boundary when cached
positive events still exist; a genuinely empty system still returns
`Stop::NoEvents`, and a still-nonpositive rebuilt total remains the truthful
`Stop::ZeroRate` fallback.

Verification on Rust 1.97.0:
- Focused deterministic regression: RED as `Stop::ZeroRate`, then GREEN.
- `cargo test --workspace`: 43 tests passed, including the kaolinite dynamics
  tests that call the paranoid differential oracle.
- Explicit `petra examples/kossel.toml --steps 1000 --paranoid`: completed 1000
  steps and wrote populations successfully.
- `git diff --check`: clean.
- Strict Clippy has no new findings from this change; both `origin/main` and this
  branch report the same three pre-existing `unnecessary_unwrap` findings in
  `petra-core/src/reaction.rs` under Rust 1.97.0.
