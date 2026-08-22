# B5-execution-schedule — piecewise-isothermal T(t) schedules

- status: done
- track: B (platform)
- priority: P1
- machine: any
- depends: B2
- claimed-by: hermes-laptop

## Objective
Add `[execution.schedule]` to the v2 deck schema: an ordered list of
`(T, duration)` segments. At each segment boundary the engine
recompiles all rate tables (same code path as load-time) and rebuilds
the event tree; within a segment everything is the existing isothermal
CTMC. No within-segment T-dependence. This is the one engine feature
the muscovite flagship requires (docs/scoping/ar-muscovite.md §4) —
step-heating experiments are piecewise-isothermal by construction.

Deliver as a small RFC-001 amendment (spec the determinism contract:
segment boundaries are wall-time events in simulated time; draw order
unchanged within segments) plus implementation plus tests.

## Acceptance
- RFC-001 amendment section merged in the same PR as the code.
- Golden gate: a two-segment schedule (T1 for t1, T2 for t2) is
  bitwise-identical to an isothermal T1 run up to t1, and statistically
  consistent with an independent T2 run initialized from that state
  (exact-replay gate for the first segment, seeded-ensemble gate for
  the second).
- v1 decks without `[execution.schedule]` behave exactly as before
  (compat gate).

## Progress
- 2026-08-22 — hermes-laptop — implemented and verified the ordered CTMC
  schedule runner, exact wall-time boundaries, temperature/thermo/propensity
  recompilation, CLI + parallel-ensemble execution, bounded reporting, and
  fail-closed artifact surfaces. Independent pre-commit review passed after
  its five resource/correctness findings were regression-tested and fixed.
- 2026-08-21 — omnibus supervisor — unblocked after B2 merged in PR #33.
- 2026-08-19 — card created (fable), motivated by the Ar-muscovite
  scoping doc.

## Result

- RFC-001 now specifies `[[execution.schedule]]`, wall-time boundary priority,
  draw consumption, continuous RNG streams, finite cumulative deadlines,
  event-count caps, and boundary/final reporting.
- `petra-core` advances exact CTMC to deadlines and rebuilds all reaction
  thermodynamics, site propensities, and the Fenwick event tree on temperature
  changes. `petra-deck` exposes `ScheduleRun`; the native CLI and deterministic
  Rayon ensemble runner execute it.
- Golden gates prove exact first-segment replay against ordinary isothermal
  `Engine::step`, seeded second-segment ensemble consistency, chemical-potential
  recompilation, and unchanged 20k-step v1 parity trajectories.
- Verification: `cargo test --workspace`, `cargo fmt --all -- --check`, and
  `cargo clippy --workspace --all-targets -- -D warnings` all pass using the
  canonical repository Cargo target.
