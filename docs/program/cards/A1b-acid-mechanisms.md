# A1b-acid-mechanisms — si-acid / al-acid hydrolysis barriers

- status: blocked
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
- 2026-08-22 12:29 PDT — Live GPU campaign closed the one-water `al-acid` model as a fully gated **sequential associative** mechanism: addition ΔG‡(298) = **82.193 kJ/mol = 19.645 kcal/mol** (205.95i cm⁻¹), cleavage local ΔG‡ = **8.047 kJ/mol = 1.923 kcal/mol** (97.46i cm⁻¹), and both quick-IRC endpoint pairs match the exact acid connectivity/proton-count basins. Mechanism-v2 `results.json` + integrity-checked `store.sqlite` are complete under `runs/phase1/al-acid-preprotonated-v2-b3lyp-def2-svp-flank/`; durable hashes and mechanism receipts are in `qm/ACID_MECHANISMS.md`. The matching one-water `si-acid` pre-equilibrium is **not an unconstrained minimum**: gas phase, a 4.5 Å water-placement probe, and aqueous PCM all return the bridge proton to H3O+ (`(False, True, True, 0, 3)`), while SMD/GPU optimization is numerically pathological. Card is blocked on a 3–6-water microsolvated Si-acid reactant; no constrained fake barrier or unsupported acid ordering is reported.
- 2026-08-21 18:20 PDT — Mandatory shared-GPU preflight still blocks the acid campaigns. Dissertation card A3's live bounded `ts2_directed_oss.py` process owns the RTX 4090 at 100% (PID 2891394, 4,244 seconds elapsed, 9,960 MiB process / 11,930 MiB total GPU memory); its durable `oss-neutral-n4-s2-pilot.log` is 3,032,960 bytes. No competing si-acid/al-acid DFT run launched. The implementation worktree remains clean before this bookkeeping checkpoint; acceptance is still red on both archived barrier results, so the mission-control pointer is parked until 2026-08-23 while this atomic board claim stays active.
- 2026-08-21 00:31 PDT — CPU implementation checkpoint complete. Added `protonated_bridge_complex()` (H3O+ pre-equilibrium → protonated intact bridge + residual H2O), a mechanism-versioned acid checkpoint namespace, acid-only approach→guarded product→CI-NEB routing, full Si/Al bridge/speciation gates for reactant/product/IRC, fail-closed cached-product validation, and crash-resumable product reuse. Cold review found four false-green/resume classes; all four were fixed with full 19-atom Si/Al regressions. Fresh exact-branch gates: **177 passed**, Ruff check, targeted Ruff format, CLI-help smoke, and diff checks green. The live A3 ladder pilot still exclusively occupies the RTX 4090 (PID 124433, 11,693 s elapsed, 100% utilization, 17,686 MiB process allocation at preflight); no competing acid DFT campaign launched. Code is ready and durable; si-acid/al-acid barriers remain pending the shared GPU lane.
- 2026-08-21 00:07 PDT — Claimed atomically by `hermes-workstation`; pushed branch `agents/A1b-acid-mechanisms` from current `origin/main` `bd8af3b8b60d5b0df1e80498443090aa208316a9` and created the congruent dedicated worktree. GPU/resource preflight and code archaeology are next.
- 2026-08-18 — card created (fable).
