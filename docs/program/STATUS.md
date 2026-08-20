# Program status log

*Newest first. One line per merged unit of work (or notable event),
appended by the PR that did it: `date — actor — what + pointer`.
Deviations get a `DEVIATION:` note; details live in the card.*

- 2026-08-20 — hermes-cloud-default — QI3 Fenwick cancellation hardening:
  impossible zero totals now rebuild from authoritative positive leaves; a
  deterministic `1e12`/`1e-6 s⁻¹` regression, all 43 Petra tests, and an
  explicit 1000-step Kossel `--paranoid` run are green.

- 2026-08-20 — hermes-workstation — E1 exact-main lint closeout: restored the
  one-line `strain_gates.rs` identity-op fix that the prior squash merge omitted;
  all Rust/Python/lint gates plus both full trajectories and generated-product
  equivalence are green.
- 2026-08-20 — hermes-workstation — E1 Ruff-format closeout: reformatted the
  cloud review-fix commit's long SVG `sy()` expression so merged-main Ruff
  0.15.1 check/format gates are clean; no simulation semantics changed.
- 2026-08-20 — hermes-workstation — E1 merged-main closeout: removed four
  default-Ruff E731 violations from the muscovite SVG post-processor with
  typed local helpers, cleared current Clippy findings, and reverified both
  full trajectories plus byte-identical CSV/SVG/JSON products.
- 2026-08-20 — hermes-workstation — E1-muscovite-deck done: schema-v1
  500/700 °C paired-OH + divacancy decks, Hames–Bowring cylinder `D/a²`
  inversion, reproducible CSV/SVG/JSON products, and qualitative gates green
  (rise 2.13×/13.55×; fall 275×/32.3×). Result:
  `docs/program/results/E1-muscovite-phase1.md`. Added executable E2 and QI3
  cards for the full mechanism and Fenwick dynamic-range hardening.
- 2026-08-20 — fable — A3-barrier-ladder claimed + platform half landed:
  crystallographic cluster builder (`qm/quarry/crystal.py`, deck-cell
  parser + bond-valence termination/charge + n_intact pruning + frozen
  shell), PHVA frequencies for frozen clusters, and the per-cell campaign
  driver (`qm/scripts/phase2_ladder.py`). Pilot cell oss-neutral running
  on the GPU; barriers land in follow-up per-family PRs.
- 2026-08-20 — fable — Card QI1: campaign drivers must enforce their own
  thread caps + niceness (a live worker probe ran 22 threads past the
  prompt-level cap — prompts don't bind subprocesses).
- 2026-08-20 — fable — Board refresh: B1+D2a+A7 done, B2 READY (keystone
  merged), A3 active. PROTOCOL amendment: incremental bookkeeping
  (the D2a iteration-ceiling lesson). Hermes worker-profile max_turns
  raised 90 -> 1000 fleet-side (infra, not repo; the raise gives workers
  more headroom, but the PROTOCOL closeout rule still applies — budget
  the last ~15% of turns so a capped worker leaves legible state).
- 2026-08-20 — local-hermes + fable — D2a DONE: 6-reaction astro rate
  campaign on the 4090. Verdict: GO gas-phase / NO-GO surface-LH rates
  (4–5 orders low without explicit-surface treatment). First non-Claude
  runtime to work the board end-to-end; worker hit iteration ceiling at
  closeout, fable finished bookkeeping. results/D2a-rate-reproduction.md.
- 2026-08-18 — hermes-workstation — A7-kinetics-database done: schema v0,
  36 P&K silicate records / 74 mechanisms, stdlib validator, and a
  49-source 2004–present expansion inventory under `kinetics-db/`.

- 2026-08-19 — hermes cloud worker — B1-schema-rfc done:
  `petra/docs/RFC-001-DECK-V2.md` written (deck schema v2 + UpdateStrategy
  trait + v1 compat shim + determinism/RNG contract + B2/B3 test plan;
  corrosion + ice-mantle showcase fragments). B2 unblocked pending Victor's
  review ack of the RFC.
- 2026-08-18 — fable — Card A7-kinetics-database added (machine:any,
  macbot-shaped): schema + extraction groundwork for a machine-readable
  successor to Palandri & Kharaka 2004. Pointer task queued in
  mission-control for laptop workers.
- 2026-08-18 — fable — Program tracking system created (PLAN, PROTOCOL,
  STATUS, 10 seed cards, inbox). Board opens with B1-schema-rfc as the
  keystone READY card.
- 2026-08-18 — fable — Omnibus program plan + three scoping docs merged
  (PR #20); root README refreshed (PR #21).
- 2026-08-18 — fable — Phase 1 first barrier: si-neutral ΔG‡(298) =
  113.0 kJ/mol (27.0 kcal/mol) vs X&L ~29, full pipeline (scan → product
  construction → CI-NEB → Sella → gates → quasi-RRHO), results +
  provenance archived. PRs #16–#18. TASK-164 queued for al-neutral.
- 2026-08-18 — fleet worker — qm/tests/test_phase1_driver.py landed
  (driver test coverage, via mission-control drain).
- 2026-08-18 — fable — Phase 0 closed: quarry package (PR #12), smoke
  v2 (PR #14); measured 4090: 74× analytic Hessians vs 32-thread CPU,
  identical energies. GPU owns all frequency work.
- 2026-08-19 — fable — qm/HANDOFF.md rewritten as the Phase-2 edition
  (barrier-ladder build plan: crystallographic cluster builder from the
  kaolinite deck cell, per-cell campaign driver, by_count emission,
  CALC-002..005 closure); new READY card A3-barrier-ladder.
- 2026-08-19 — fable — Track E opened: Ar-in-muscovite scoping doc
  landed (docs/scoping/ar-muscovite.md — the 1998 Sletten & Onstott
  circle, Nteme 2022–24 barriers, the uncomputed middle, Villa
  arbitration). New cards: E1-muscovite-deck (READY, published numbers
  only) and B5-execution-schedule (blocked: B2 — piecewise-isothermal
  T(t), the one engine feature the flagship needs).
- 2026-08-20 — Fable — Autopilot engaged: board_feeder.py (mission-
  control, hourly) auto-enqueues READY/unblocked cards as queue pointer
  tasks; omnibus-supervisor Hermes cron (daily 07:50) does bookkeeping
  + stale claims + Telegram digest. PROTOCOL.md amended. Humans and
  Fable are now exception handlers, not schedulers.
- 2026-08-20 — fable — card QI2-gpu-lease (READY, P1): explicit
  single-lane GPU lease + 16 GB cupy cap + shared worker preflight,
  after reviewing overnight A3-pilot/ollama cohabitation (22 ollama
  loads, all-GPU, zero errors — but pilot peaked 21.7 GB; luck, not
  design) and TASK-168's per-tick hand-rolled preflights.
- 2026-08-20 — fable — PROTOCOL.md merge-doctrine paragraph aligned to
  POLICY v7 (open PR and STOP; watcher owns lifecycle; human-merge
  label) — the stale auto-merge/Sourcery wording was flagged by the
  cloud Hermes worker during its enrollment (it correctly refused to
  follow it). Companion factory fixes landed in mission-control.
