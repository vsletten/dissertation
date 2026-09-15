# E3b3-macos-timeout-cleanup-test — synchronize descendant cleanup fixture

- status: done
- track: E (muscovite / the perfect circle)
- priority: P2
- machine: any (synthetic test only)
- depends: E3b1-periodic-dft-evidence-contract
- requested-by: hermes-macbot-zero

## Objective

Make `test_smoke_timeout_terminates_descendant_process_group` deterministic on
macOS and Linux by proving the synthetic descendant has installed its SIGTERM
handler before the parent advertises readiness and the 0.1-second timeout can
fire. Preserve the production cleanup contract: timeout sends SIGTERM to the
process group, escalates to SIGKILL when a descendant survives, reaps the direct
child, and leaves no delayed marker.

## Constraints

- POLICY.md and PROTOCOL.md govern.
- This is test-fixture synchronization, not scientific compute and not a license
  to weaken the process-group cleanup assertions.
- Keep the production timeout/cleanup behavior unchanged unless an independently
  reproduced implementation defect requires a minimal correction.

## Acceptance

- The descendant fixture communicates readiness only after installing its
  SIGTERM handler; the parent relays that state before sleeping.
- The focused timeout test passes repeatedly on macOS and in the repository's
  Linux CI while retaining the `timeout_term_sent`, `timeout_kill_sent`, reaped,
  and delayed-marker assertions.
- Focused Ruff/format/compile/tests pass.
- Card/PLAN/STATUS bookkeeping lands in the implementation PR.

## Progress

- 2026-09-14 17:57 PDT — `(hermes-macbot-zero; profile=laptop)` — Filed under
  POLICY §13 after E3b2's optional expanded E3 test run passed 44/45 and the
  unchanged origin/main timeout fixture failed twice on macOS: the descendant
  can receive SIGTERM before it installs `SIG_IGN`, so the receipt truthfully
  reports that SIGKILL escalation was unnecessary. E3b2's required classical
  receipt suite remains green 9/9.

- 2026-09-14 19:56 PDT — `(hermes-macbot-zero; profile=laptop)` — DONE. Reproduced
  the original macOS race, then made readiness causal: the shell descendant
  installs its SIGTERM-ignore trap before signaling, the shell parent relays
  readiness, and a fail-closed test-only gate starts the unchanged 0.1-second
  cleanup timeout only after that relay. One cold review caught that merely
  speeding up the fixture still left the outer timeout race; this gate is
  rejected for every non-fixture executable identity, so production smoke
  behavior is unchanged.

## Result

- The macOS focused timeout test passes 20/20 repeated runs while retaining
  `timeout_term_sent`, `timeout_kill_sent`, direct-child reaping, and delayed-
  marker assertions.
- `tests/test_e3b_periodic_dft.py`: 30/30 passed.
- Focused Ruff passed with baseline-only `C409`/`TRY004` findings ignored;
  Ruff format, Python compile, and `git diff --check` passed. Linux execution
  is covered by the repository PR CI gate.
