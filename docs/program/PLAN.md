# Program board — the Omnibus, tracked

*The live board for [`docs/OMNIBUS.md`](../OMNIBUS.md). Workers: read
[`PROTOCOL.md`](PROTOCOL.md) before touching anything. This file is
updated through the same PRs that do the work.*

## Milestone state

| Milestone | State | Notes |
|---|---|---|
| M0 — scoping + in-flight landed | ✅ done 2026-08-18 | Omnibus + 3 scoping docs merged; Phase-1 si-neutral barrier landed (ΔG‡ 27.0 kcal/mol vs X&L ~29) |
| M1 — B1 schema RFC reviewed | ✅ done 2026-08-19 | RFC-001 merged (PR #25); B2 unblocked |
| M2 — B2 refactor, parity green | in progress | B2 READY as of B1 merge |
| M3 — B3 conformance decks green | blocked (B2) | Conway / Ising / SIR with analytic gates |
| M4 — B4 ensembles + observables | blocked (B2) | design against A5-Phase-0 needs |
| M5 — first science on new platform | blocked (M3/M4) | A5 aging study ∥ C2 corrosion deck ∥ D2 rates campaign |
| M6 — first external deliverables | blocked (M5) | tool paper, corrosion statistics, KIDA submissions, field-lab mechanism paper |

## Cards

Status: `ready` (claimable now) · `blocked` (named blocker) ·
`active` (live `agents/<id>` branch exists — check
`git ls-remote --heads origin 'agents/*'`) · `done` (acceptance passes
on main). Priority P0 > P1 > P2 within READY.

| Card | Track | Priority | Machine | Status | Depends on |
|---|---|---|---|---|---|
| [B1-schema-rfc](cards/B1-schema-rfc.md) | B | P0 | any | done | — |
| [B2-engine-refactor](cards/B2-engine-refactor.md) | B | P0 | any | ready | B1 ✅ |
| [A1b-acid-mechanisms](cards/A1b-acid-mechanisms.md) | A | P1 | workstation | ready | — |
| [A3-barrier-ladder](cards/A3-barrier-ladder.md) | A | P1 | workstation | ready | — |
| [A7-kinetics-database](cards/A7-kinetics-database.md) | A | P2 | any | done | — |
| [A2-production-energetics](cards/A2-production-energetics.md) | A | P1 | workstation | blocked: ORCA install (Victor) + A1b | A1b |
| [B3-conformance-decks](cards/B3-conformance-decks.md) | B | P1 | any | blocked: B2 | B2 |
| [B4-ensembles-observables](cards/B4-ensembles-observables.md) | B | P1 | any | blocked: B2 | B2 |
| [C2-corrosion-deck](cards/C2-corrosion-deck.md) | C | P2 | any | blocked: B3 | B3 |
| [D3-ice-mantle-deck](cards/D3-ice-mantle-deck.md) | D | P2 | any | blocked: B3 + D2a | B3, D2a |
| [A5p0-aging-observables](cards/A5p0-aging-observables.md) | A | P1 | any | blocked: B4 | B4 |

**Done** (acceptance verified on main): B1-schema-rfc (PR #25),
D2a-astro-rate-reproduction (PR #31 — verdict: GO gas-phase /
NO-GO surface-LH), A7-kinetics-database (PR #27 — 36 minerals, 74
mechanisms, validator green). **Active**: A3-barrier-ladder
(claim branch live).

**Tracked elsewhere**: TASK-164 (scan-smoothness repair + al-neutral
barrier) predates this board and lives in the mission-control queue —
do not duplicate; its completion should be logged in STATUS.md here by
whoever lands it.

## Standing context for cold workers

- Program rationale + track detail: [`docs/OMNIBUS.md`](../OMNIBUS.md)
- Research grounding: [`docs/scoping/`](../scoping/) (field-lab
  discrepancy, pitting corrosion, astrochemistry)
- Engine: [`petra/docs/DESIGN.md`](../../petra/docs/DESIGN.md);
  QM lab: [`qm/SURVEY.md`](../../qm/SURVEY.md),
  [`qm/README.md`](../../qm/README.md)
- Institutional memory: [`.claude/learnings/`](../../.claude/learnings/)
  — read INDEX.md; it earns its keep (petra deck gotchas, PR-merge
  blockers, saddle-search traps, GPU benchmarking artifacts)
