# A3-barrier-ladder — connectivity × protonation barrier ladder (Phase 2)

- status: ready
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaigns; cluster-builder code is machine-any)
- depends: — (runs at the validated b3lyp/def2-svp/df tier; A2 re-tiers
  the finished ladder later; acid-state variants may reuse A1b's
  protonated-bridge machinery when it lands but are not blocked by it)
- claimed-by:

## Objective
Execute qm/HANDOFF.md (Phase 2 edition — the full build plan lives
there): build the crystallographic cluster builder (terminated,
peripherally-frozen edge clusters cut from the kaolinite cell in
petra/examples/kaolinite.toml, per KMC site family and protonation
state), generalize the Phase-1 driver into a per-ladder-cell campaign
script, compute the connectivity- and protonation-resolved barriers for
the five site families, and emit them as petra `by_count`/`when` tables
with provenance. Closes CALCULATIONS.md rows CALC-002..005 from
`needed` to `computed`.

## Context
- qm/HANDOFF.md — the complete Phase-2 plan (read first)
- qm/SURVEY.md §6.2 (cluster norms), §7.3 (pKa table, kaolinite targets)
- qm/CALCULATIONS.md CALC-002..005 (the shopping list)
- Repo skills: quarry-campaign, petra-deck, pr-merge-gauntlet
- Reference result: si-neutral 113.05 kJ/mol (27.0 kcal/mol) in the
  claude/qm-phase1-scan-extend worktree runs/ dir
- .claude/learnings/ — the measured TS-hunt failure modes and gates

## Acceptance
- Cluster-builder unit gates green (CPU-only): taxonomy match against
  the deck cell, termination/charge bookkeeping, frozen shell
  peripheral, no collisions.
- Pilot cell (300s Si–O–Si, neutral, reference connectivity) reproduces
  the si-neutral barrier modulo a logged lattice-resistance shift.
- Per completed family: gated barriers in runs/ + store, an emitted
  petra fragment that compiles via petra-cli round-trip, and the
  CALCULATIONS.md row updated — one PR per family.
- Full fast test suite + ruff green on every PR.

## Progress
- 2026-08-19 — card created with the Phase-2 handoff rewrite (fable).
