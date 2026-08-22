# Program board — the Omnibus, tracked

*The live board for [`docs/OMNIBUS.md`](../OMNIBUS.md). Workers: read
[`PROTOCOL.md`](PROTOCOL.md) before touching anything. This file is
updated through the same PRs that do the work.*

## Milestone state

| Milestone | State | Notes |
|---|---|---|
| M0 — scoping + in-flight landed | ✅ done 2026-08-18 | Omnibus + 3 scoping docs merged; Phase-1 si-neutral barrier landed (ΔG‡ 27.0 kcal/mol vs X&L ~29) |
| M1 — B1 schema RFC reviewed | ✅ done 2026-08-19 | RFC-001 merged (PR #25); B2 unblocked |
| M2 — B2 refactor, parity green | ✅ done 2026-08-20 | B2 merged in PR #33; schema v2 + bitwise-parity gates green |
| M3 — B3 conformance decks green | ✅ done 2026-08-21 | Conway / Ising / SIR analytic gates green; all four update strategies executable |
| M4 — B4 ensembles + observables | ✅ done 2026-08-22 | Parallel replicas, distributions/bootstrap CIs, first A5 observables |
| M5 — first science on new platform | ready | M3/M4 done; A5 aging study ∥ C2 corrosion deck ∥ D2 rates campaign |
| M6 — first external deliverables | blocked (M5) | tool paper, corrosion statistics, KIDA submissions, field-lab mechanism paper |

## Cards

Status: `ready` (claimable now) · `blocked` (named blocker) ·
`active` (live `agents/<id>` branch exists — check
`git ls-remote --heads origin 'agents/*'`) · `done` (acceptance passes
on main). Priority P0 > P1 > P2 within READY.

| Card | Track | Priority | Machine | Status | Depends on |
|---|---|---|---|---|---|
| [B1-schema-rfc](cards/B1-schema-rfc.md) | B | P0 | any | done | — |
| [B2-engine-refactor](cards/B2-engine-refactor.md) | B | P0 | any | done | B1 ✅ |
| [QI1-driver-enforced-etiquette](cards/QI1-driver-enforced-etiquette.md) | A | P2 | any | done | — |
| [QI2-gpu-lease](cards/QI2-gpu-lease.md) | A | P1 | any | ready | — |
| [QI3-fenwick-total-cancellation](cards/QI3-fenwick-total-cancellation.md) | B | P1 | any | done | — |
| [E1-muscovite-deck](cards/E1-muscovite-deck.md) | E | P1 | any | done | — (phase 1 isothermal) |
| [E2-muscovite-full-mechanism](cards/E2-muscovite-full-mechanism.md) | E | P1 | any | done | B5 ✅ |
| [A8-thesis-archive-intake](cards/A8-thesis-archive-intake.md) | A | P1 | workstation | ready | — |
| [A1b-acid-mechanisms](cards/A1b-acid-mechanisms.md) | A | P1 | workstation | ready | — |
| [A3-barrier-ladder](cards/A3-barrier-ladder.md) | A | P1 | workstation | ready | — |
| [A7-kinetics-database](cards/A7-kinetics-database.md) | A | P2 | any | done | — |
| [A2-production-energetics](cards/A2-production-energetics.md) | A | P1 | workstation | blocked: ORCA install (Victor) + A1b | A1b |
| [B3-conformance-decks](cards/B3-conformance-decks.md) | B | P1 | any | done | B2 ✅ |
| [B4-ensembles-observables](cards/B4-ensembles-observables.md) | B | P1 | any | done | B2 ✅ |
| [B5-execution-schedule](cards/B5-execution-schedule.md) | B | P1 | any | done | B2 ✅ |
| [C2-corrosion-deck](cards/C2-corrosion-deck.md) | C | P2 | any | ready | B3 ✅ |
| [D2a-astro-rate-reproduction](cards/D2a-astro-rate-reproduction.md) | D | P1 | workstation | done | — |
| [D2b-explicit-surface-rates](cards/D2b-explicit-surface-rates.md) | D | P1 | workstation | ready | D2a ✅ |
| [D3-ice-mantle-deck](cards/D3-ice-mantle-deck.md) | D | P2 | any | done | B3 ✅, D2a ✅ |
| [D3b-co-hydrogenation-deck](cards/D3b-co-hydrogenation-deck.md) | D | P2 | any | blocked: D2b | D3 ✅, D2b |
| [A5p0-aging-observables](cards/A5p0-aging-observables.md) | A | P1 | any | done | B4 ✅ |

**Done** (acceptance verified on main): D3-ice-mantle-deck,
E2-muscovite-full-mechanism,
B5-execution-schedule,
A5p0-aging-observables,
B4-ensembles-observables,
B3-conformance-decks,
QI1-driver-enforced-etiquette,
QI3-fenwick-total-cancellation (PR #52), B1-schema-rfc (PR #25),
B2-engine-refactor (PR #33),
D2a-astro-rate-reproduction (PR #31 — verdict: GO gas-phase /
NO-GO surface-LH), A7-kinetics-database (PR #27 — 36 minerals, 74
mechanisms, validator green). **Active claim**: A1b-acid-mechanisms.
A3-barrier-ladder is ready for its remaining family campaigns after its
pilot result merged in PR #67.

**Tracked elsewhere**: TASK-164 (scan-smoothness repair + al-neutral
barrier) predates this board and lives in the mission-control queue —
do not duplicate; its completion should be logged in STATUS.md here by
whoever lands it.

## Standing context for cold workers

- Program rationale + track detail: [`docs/OMNIBUS.md`](../OMNIBUS.md)
- Research grounding: [`docs/scoping/`](../scoping/) (field-lab
  discrepancy, pitting corrosion, astrochemistry, Ar-in-muscovite)
- Engine: [`petra/docs/DESIGN.md`](../../petra/docs/DESIGN.md);
  QM lab: [`qm/SURVEY.md`](../../qm/SURVEY.md),
  [`qm/README.md`](../../qm/README.md)
- Institutional memory: [`.claude/learnings/`](../../.claude/learnings/)
  — read INDEX.md; it earns its keep (petra deck gotchas, PR-merge
  blockers, saddle-search traps, GPU benchmarking artifacts)
