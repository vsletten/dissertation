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
| [QI2-gpu-lease](cards/QI2-gpu-lease.md) | A | P1 | any | done | — |
| [QI3-fenwick-total-cancellation](cards/QI3-fenwick-total-cancellation.md) | B | P1 | any | done | — |
| [E1-muscovite-deck](cards/E1-muscovite-deck.md) | E | P1 | any | done | — (phase 1 isothermal) |
| [E2-muscovite-full-mechanism](cards/E2-muscovite-full-mechanism.md) | E | P1 | any | done | B5 ✅ |
| [A8-thesis-archive-intake](cards/A8-thesis-archive-intake.md) | A | P1 | workstation | ready | — |
| [A1b-acid-mechanisms](cards/A1b-acid-mechanisms.md) | A | P1 | workstation | blocked: A1i | production-tier revisit after A2 |
| [A1c-acid-microsolvation-conformers](cards/A1c-acid-microsolvation-conformers.md) | A | P1 | workstation | done | 16/16 conclusive Si rejections via A1d |
| [A1d-acid-3w-bridge-convergence](cards/A1d-acid-3w-bridge-convergence.md) | A | P1 | workstation | done | A1c infrastructure + archived receipts |
| [A1e-acid-concerted-hydronium-relay](cards/A1e-acid-concerted-hydronium-relay.md) | A | P1 | workstation | done | 4/4 conclusive hydronium-as-nucleophile rejections |
| [A1f-acid-neutral-water-attacker-relay](cards/A1f-acid-neutral-water-attacker-relay.md) | A | P1 | workstation | done | exact donor migrated to relay |
| [A1g-acid-bridge-side-hydronium-neutral-attacker](cards/A1g-acid-bridge-side-hydronium-neutral-attacker.md) | A | P1 | workstation | done | conclusive product-family rejection |
| [A1i-acid-production-tier-bridge-side-revisit](cards/A1i-acid-production-tier-bridge-side-revisit.md) | A | P1 | workstation | blocked: A2 | A1g ✅, A2 |
| [A3-barrier-ladder](cards/A3-barrier-ladder.md) | A | P1 | workstation | ready | — |
| [A7-kinetics-database](cards/A7-kinetics-database.md) | A | P2 | any | done | — |
| [A2-production-energetics](cards/A2-production-energetics.md) | A | P1 | workstation | blocked: A2a | banked si-neutral TS fails production full-IRC |
| [A2a-si-neutral-production-path-rebuild](cards/A2a-si-neutral-production-path-rebuild.md) | A | P1 | workstation | active | exact r2SCAN-3c minima + A2 infrastructure ✅ |
| [B3-conformance-decks](cards/B3-conformance-decks.md) | B | P1 | any | done | B2 ✅ |
| [B4-ensembles-observables](cards/B4-ensembles-observables.md) | B | P1 | any | done | B2 ✅ |
| [B5-execution-schedule](cards/B5-execution-schedule.md) | B | P1 | any | done | B2 ✅ |
| [C2-corrosion-deck](cards/C2-corrosion-deck.md) | C | P2 | any | done | B3 ✅ |
| [C3-pit-statistics](cards/C3-pit-statistics.md) | C | P2 | any | ready | C2 ✅ |
| [D2a-astro-rate-reproduction](cards/D2a-astro-rate-reproduction.md) | D | P1 | workstation | done | — |
| [D2b-explicit-surface-rates](cards/D2b-explicit-surface-rates.md) | D | P1 | workstation | done | gate NO-GO with receipts; barrier tier validated |
| [D2c-instanton-tier](cards/D2c-instanton-tier.md) | D | P2 | workstation | ready | D2b ✅ |
| [D3-ice-mantle-deck](cards/D3-ice-mantle-deck.md) | D | P2 | any | done | B3 ✅, D2a ✅ |
| [D3b-co-hydrogenation-deck](cards/D3b-co-hydrogenation-deck.md) | D | P2 | any | blocked: D2c | D3 ✅, D2b ✅ (table withheld — D2c gate) |
| [A5p0-aging-observables](cards/A5p0-aging-observables.md) | A | P1 | any | done | B4 ✅ |
| [A5p1-aging-study](cards/A5p1-aging-study.md) | A | P1 | any | done | A5p0 ✅ |
| [E2b-grain-size-sweep](cards/E2b-grain-size-sweep.md) | E | P1 | any | done | E2 ✅, B4 ✅ |
| [E3a-classical-neb-barriers](cards/E3a-classical-neb-barriers.md) | E | P1 | any | ready | — |
| [E3b-periodic-dft-spot-checks](cards/E3b-periodic-dft-spot-checks.md) | E | P2 | workstation | blocked: E3a | E3a |
| [E4-1998-comparison](cards/E4-1998-comparison.md) | E | P1 | any | blocked: E3a, E2b, A8 | E3a, E2b, A8 |

**Done** (acceptance verified on main): E2b-grain-size-sweep,
A5p1-aging-study, C2-corrosion-deck,
QI2-gpu-lease,
A1g-acid-bridge-side-hydronium-neutral-attacker,
A1f-acid-neutral-water-attacker-relay,
A1e-acid-concerted-hydronium-relay,
A1d-acid-3w-bridge-convergence,
A1c-acid-microsolvation-conformers, D3-ice-mantle-deck,
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
mechanisms, validator green). **Blocked closeout branch**:
`agents/A1b-acid-microsolvation` for A1b-acid-mechanisms.
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
