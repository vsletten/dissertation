# Program status log

*Newest first. One line per merged unit of work (or notable event),
appended by the PR that did it: `date — actor — what + pointer`.
Deviations get a `DEVIATION:` note; details live in the card.*

- 2026-08-19 — fable — A3-barrier-ladder claimed + platform half landed:
  crystallographic cluster builder (`qm/quarry/crystal.py`, deck-cell
  parser + bond-valence termination/charge + n_intact pruning + frozen
  shell), PHVA frequencies for frozen clusters, and the per-cell campaign
  driver (`qm/scripts/phase2_ladder.py`). Pilot cell oss-neutral running
  on the GPU; barriers land in follow-up per-family PRs.

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
