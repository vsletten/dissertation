# A1b-acid-mechanisms — si-acid / al-acid hydrolysis barriers

- status: active
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaign)
- depends: — (TASK-164's scan-smoothness repair helps but is not required)
- claimed-by: hermes-workstation (see agents/A1b-acid-mechanisms branch)

## Objective
Extend scripts/phase1_xiao_lasaga.py with the acid mechanism: H3O+ first
protonates the bridging O (barrierless or low-barrier proton transfer —
model as a pre-equilibrated protonated-bridge complex), THEN water
attacks the weakened Si-Obr bond. This is a different stage-2 design
from the neutral path: build the protonated complex directly
(clusters.py: protonated_bridge_complex()), then approach-scan + the
existing product-construction/CI-NEB/Sella ladder. Compute si-acid and
al-acid barriers at b3lyp/def2-svp/df; compare vs X&L acid ~24 kcal/mol
and the si-neutral 27.0 kcal/mol (acid must come out LOWER).

## Context
- docs/scoping + qm/HANDOFF.md Phase 1; X&L 1994 mechanism section
- runs/phase1/si-neutral-* in the qm-phase1-scan-extend worktree
  (reference results + checkpoints pattern)
- .claude/learnings/debugging.md (saddle traps), gotchas.md

## Acceptance
- Both barriers pass the structural gates (1 imaginary mode, in-channel
  TS), results.json + store.sqlite archived, logs teed.
- Ordering reported: acid < neutral, Si-O-Al vs Si-O-Si stated.
- New clusters/driver code unit-tested; suite green; PR merged.

## Progress
- 2026-08-21 00:07 PDT — Claimed atomically by `hermes-workstation`; pushed branch `agents/A1b-acid-mechanisms` from current `origin/main` `bd8af3b8b60d5b0df1e80498443090aa208316a9` and created the congruent dedicated worktree. GPU/resource preflight and code archaeology are next.
- 2026-08-18 — card created (fable).
