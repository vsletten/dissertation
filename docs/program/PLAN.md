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
| [A8-thesis-archive-intake](cards/A8-thesis-archive-intake.md) | A | P1 | workstation | done | — |
| [A8a-legacy-kmc-conformance](cards/A8a-legacy-kmc-conformance.md) | A | P2 | any | done | A8 ✅ |
| [A8b-legacy-dft-transition-state-gap](cards/A8b-legacy-dft-transition-state-gap.md) | A | P2 | workstation | blocked: A3 | archival ledger → modern barriers |
| [A8c-thesis-table-4-10c-erratum](cards/A8c-thesis-table-4-10c-erratum.md) | A | P2 | any | done | A8 ✅; 185.6 source-supported, 183.6 transcription error |
| [A8d-ligand-promotion-track-scoping](cards/A8d-ligand-promotion-track-scoping.md) | A | P2 | any | done | A8 ✅; CONDITIONAL GO for bounded dual-site oxalate pilot |
| [A8e-dual-site-oxalate-pilot](cards/A8e-dual-site-oxalate-pilot.md) | A | P2 | workstation | blocked: A3 | A8d ✅; six-system falsification pilot only |
| [A1b-acid-mechanisms](cards/A1b-acid-mechanisms.md) | A | P1 | workstation | blocked: A1i | production-tier revisit after A2 |
| [A1c-acid-microsolvation-conformers](cards/A1c-acid-microsolvation-conformers.md) | A | P1 | workstation | done | 16/16 conclusive Si rejections via A1d |
| [A1d-acid-3w-bridge-convergence](cards/A1d-acid-3w-bridge-convergence.md) | A | P1 | workstation | done | A1c infrastructure + archived receipts |
| [A1e-acid-concerted-hydronium-relay](cards/A1e-acid-concerted-hydronium-relay.md) | A | P1 | workstation | done | 4/4 conclusive hydronium-as-nucleophile rejections |
| [A1f-acid-neutral-water-attacker-relay](cards/A1f-acid-neutral-water-attacker-relay.md) | A | P1 | workstation | done | exact donor migrated to relay |
| [A1g-acid-bridge-side-hydronium-neutral-attacker](cards/A1g-acid-bridge-side-hydronium-neutral-attacker.md) | A | P1 | workstation | done | conclusive product-family rejection |
| [A1i-acid-production-tier-bridge-side-revisit](cards/A1i-acid-production-tier-bridge-side-revisit.md) | A | P1 | workstation | blocked: A2 | A1g ✅, A2 |
| [A9-approximate-rate-closure](cards/A9-approximate-rate-closure.md) | A | P0 | any | done | honest partial: zero observed lattice release; all seven families censored/unrankable |
| [A9b-mechanism-reachability](cards/A9b-mechanism-reachability.md) | A | P0 | any | done | 260 live release paths + 60 boundary topology NO-GOs; finite biased sampler ✅ |
| [A9b-reservoir-origin-contract](cards/A9b-reservoir-origin-contract.md) | A | P0 | any | done | A9 ✅; pH 3–5 open-flow contract + origin-safe release accounting |
| [A9b-sensitivity-ranking](cards/A9b-sensitivity-ranking.md) | A | P0 | workstation | done | verified 29 × 8 campaign; all responses censored/unrankable; no irrelevant family |
| [A3-barrier-ladder](cards/A3-barrier-ladder.md) | A | P1 | workstation | blocked: victor | choose an observable-release regime or revise the model after A9b's all-censored result |
| [A3a-reactant-minimum-recovery](cards/A3a-reactant-minimum-recovery.md) | A | P1 | workstation | done | terminal microstate rejection; no minimum promoted |
| [A3b-osa-neutral-n1-proton-microstate-stability](cards/A3b-osa-neutral-n1-proton-microstate-stability.md) | A | P1 | workstation | done | verified inconclusive; H52 became third mobile owner |
| [A3c-osa-neutral-n1-mobile-proton-triad-conditioning](cards/A3c-osa-neutral-n1-mobile-proton-triad-conditioning.md) | A | P1 | workstation | done | verified triad-conditioning failure; owners retained, projected gradients red |
| [A3d-osa-neutral-n1-triad-failure-adjudication](cards/A3d-osa-neutral-n1-triad-failure-adjudication.md) | A | P1 | workstation | done | verified terminal rejection; A3 returns to remaining campaigns |
| [A3e-oss-neutral-n1-proton-microstate-stability](cards/A3e-oss-neutral-n1-proton-microstate-stability.md) | A | P1 | workstation | done | verified inconclusive; constrained max-gradient gate red, release not run |
| [A3f-oss-neutral-n1-terminal-adjudication](cards/A3f-oss-neutral-n1-terminal-adjudication.md) | A | P1 | any | done | verified terminal rejection; no global no-basin claim |
| [A3g-oaa-neutral-n2-proton-microstate-stability](cards/A3g-oaa-neutral-n2-proton-microstate-stability.md) | A | P1 | workstation | done | verified inconclusive; conditioning nonstationary/owner-changing, production not run |
| [A3h-calc005-si-attachment-protocol](cards/A3h-calc005-si-attachment-protocol.md) | A | P1 | any | done | live state-205 cycle + state-expanded C/V proof; no energy emitted |
| [A3i-calc005-si-n1-pilot](cards/A3i-calc005-si-n1-pilot.md) | A | P1 | workstation | done | terminal C-optimizer failure; no value emitted |
| [A3j-calc005-si-n1-verification](cards/A3j-calc005-si-n1-verification.md) | A | P1 | workstation | done | closed by ruling; A3i receipt stands |
| [A7-kinetics-database](cards/A7-kinetics-database.md) | A | P2 | any | done | — |
| [A2-production-energetics](cards/A2-production-energetics.md) | A | P1 | workstation | blocked: A9b-sensitivity-ranking | final banked survey/literature table after corrected ranking |
| [A2a-si-neutral-production-path-rebuild](cards/A2a-si-neutral-production-path-rebuild.md) | A | P1 | workstation | done | exact r2SCAN-3c minima + A2 infrastructure ✅ |
| [A2b-al-neutral-production-energetics](cards/A2b-al-neutral-production-energetics.md) | A | P1 | workstation | blocked: victor | banked 32.2 kcal/mol survey value; no supported replication target after A9b |
| [A2b1-wb97mv-reactant-scf-recovery](cards/A2b1-wb97mv-reactant-scf-recovery.md) | A | P1 | workstation | done | independently verified finite SCF failure; no retry authorized |
| [A2c-al-acid-production-energetics](cards/A2c-al-acid-production-energetics.md) | A | P1 | workstation | blocked: A2b | A1b banked one-water Al-acid route |
| [A2d-oss-neutral-n4-production-energetics](cards/A2d-oss-neutral-n4-production-energetics.md) | A | P1 | workstation | blocked: A2c | documented embedded pilot; deterministic rebuild |
| [B3-conformance-decks](cards/B3-conformance-decks.md) | B | P1 | any | done | B2 ✅ |
| [B4-ensembles-observables](cards/B4-ensembles-observables.md) | B | P1 | any | done | B2 ✅ |
| [B5-execution-schedule](cards/B5-execution-schedule.md) | B | P1 | any | done | B2 ✅ |
| [C2-corrosion-deck](cards/C2-corrosion-deck.md) | C | P2 | any | done | B3 ✅ |
| [C3-pit-statistics](cards/C3-pit-statistics.md) | C | P2 | any | done | C2 ✅ |
| [D2a-astro-rate-reproduction](cards/D2a-astro-rate-reproduction.md) | D | P1 | workstation | done | — |
| [D2b-explicit-surface-rates](cards/D2b-explicit-surface-rates.md) | D | P1 | workstation | done | gate NO-GO with receipts; barrier tier validated |
| [D2c-instanton-tier](cards/D2c-instanton-tier.md) | D | P2 | workstation | ready | D2b ✅ |
| [D3-ice-mantle-deck](cards/D3-ice-mantle-deck.md) | D | P2 | any | done | B3 ✅, D2a ✅ |
| [D3b-co-hydrogenation-deck](cards/D3b-co-hydrogenation-deck.md) | D | P2 | any | blocked: D2c | D3 ✅, D2b ✅ (table withheld — D2c gate) |
| [A5p0-aging-observables](cards/A5p0-aging-observables.md) | A | P1 | any | done | B4 ✅ |
| [A5p1-aging-study](cards/A5p1-aging-study.md) | A | P1 | any | done | A5p0 ✅ |
| [E2b-grain-size-sweep](cards/E2b-grain-size-sweep.md) | E | P1 | any | done | E2 ✅, B4 ✅ |
| [E3a-classical-neb-barriers](cards/E3a-classical-neb-barriers.md) | E | P1 | any | done | — |
| [E3b-periodic-dft-spot-checks](cards/E3b-periodic-dft-spot-checks.md) | E | P2 | workstation | done | E3b1 ✅; 3/3 PBE-D3 + matched-cell frozen-path profiles; all `disagrees`, no correction |
| [E3b1-periodic-dft-evidence-contract](cards/E3b1-periodic-dft-evidence-contract.md) | E | P2 | workstation | done | E3a ✅; bounded force timing + native evidence gates verified |
| [E3b2-classical-receipt-code-binding](cards/E3b2-classical-receipt-code-binding.md) | E | P2 | any | done | E3b ✅; receipt v2 binds exact runner bytes/revision plus explicit invocation/operator |
| [E3b3-macos-timeout-cleanup-test](cards/E3b3-macos-timeout-cleanup-test.md) | E | P2 | any | done | E3b1 ✅; readiness-gated descendant cleanup fixture |
| [E4-1998-comparison](cards/E4-1998-comparison.md) | E | P1 | any | done | E3a ✅, E2b ✅, A8 ✅ |
| [E4a-surface-gated-release](cards/E4a-surface-gated-release.md) | E | P1 | any | done | E4 ✅, B5 ✅; NO-GO — local basal gating did not recover the grain-size crossover |
| [E4a2-surface-connected-lateral-release](cards/E4a2-surface-connected-lateral-release.md) | E | P1 | any | done | E4a ✅; PRs #131 + #132; NO-GO under inherited gate timescale + proxy-size ladder |
| [E4b-isothermal-reservoir-discriminants](cards/E4b-isothermal-reservoir-discriminants.md) | E | P1 | any | done | E4a ✅; §5.1/.4/.7 partial, §5.2 not reproduced |

**Done** (acceptance verified on main): E3b2-classical-receipt-code-binding,
A9-approximate-rate-closure,
E4a-surface-gated-release (PR #129;
NO-GO — local basal gating did not recover the grain-size crossover),
E4-1998-comparison (PR #127),
A3j-calc005-si-n1-verification,
E3a-classical-neb-barriers
(bounded reconstructed route-1 gate 68.414811 vs 67.644151 kcal/mol, +1.1393%;
all six NEBs honestly `incomplete-convergence`), A8a-legacy-kmc-conformance (five replays complete;
DEVIATION: all five behaviorally mismatch the archive, so canonical conformance
remains red), A3g-oaa-neutral-n2-proton-microstate-stability,
A3f-oss-neutral-n1-terminal-adjudication,
A3d-osa-neutral-n1-triad-failure-adjudication,
A3c-osa-neutral-n1-mobile-proton-triad-conditioning,
A3b-osa-neutral-n1-proton-microstate-stability,
A2a-si-neutral-production-path-rebuild,
C3-pit-statistics, E2b-grain-size-sweep,
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
A9-approximate-rate-closure is DONE by the 2026-09-13 scope ruling with an honest partial/negative survey result: 46.4 million events observed zero original-lattice Si/Al release, the corrected gate types the mechanism unsampled, and all seven sensitivity families remain censored/unrankable. A9b mechanism reachability, reservoir origin, and the corrected 29 × 8 sensitivity campaign are now DONE. The 232/232 stationary A9b trajectories still observed zero biased original-lattice release, so A9b names no top-k and no irrelevant family. A3 and A2b are blocked on Victor's program-level choice of an observable-release regime, mechanism/deck revision, or accelerated rare-event method; neither card authorizes QM from this non-observation. The recommended cheapest discriminator is the approximate deck at 350–450 K and/or pH 3, followed by sensitivity ranking in a regime where release is observable. If a later decision and result identifies a sensitive family, A3 rebuilds the closest published cluster and method and reproduces its observable within ±3 kcal/mol; the survey-tier or literature-replication value is then banked. One family per PR carries the paper/replication/banked-value table and Petra fragment. The prior home-grown hydrated-cell
Osa/OSS/Oaa/CALC-005 attempts remain historical null evidence and are not to be
replayed. A3j is DONE by ruling; A3i's hash-bound terminal optimizer-failure
receipt stands and no second-worker null-result ceremony remains.

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
