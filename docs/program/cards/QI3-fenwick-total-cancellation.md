# QI3-fenwick-total-cancellation — preserve slow events beside extreme rates

- status: ready
- track: B (platform quality infrastructure)
- priority: P1
- machine: any
- depends: —
- claimed-by:

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
