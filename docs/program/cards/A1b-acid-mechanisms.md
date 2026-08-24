# A1b-acid-mechanisms — si-acid / al-acid hydrolysis barriers

- status: blocked
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaign)
- depends: — (TASK-164's scan-smoothness repair helps but is not required)
- claimed-by: hermes-workstation (continuation: agents/A1b-acid-microsolvation)
- blocked-on: A1c-acid-microsolvation-conformers

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
- 2026-08-23 20:21 PDT — **DEVIATION / SCIENTIFIC BLOCKER:** a bounded
  matched microsolvation count survey did not rescue the Si protonated bridge.
  Production B3LYP/def2-SVP/DF endpoints converged for 3, 4 (continued beyond
  its initial 100-step bound), 5, and 6 total waters. Every endpoint retained
  intact Si--O--Si but moved the bridge proton before any Hessian or TS search:
  signatures were respectively `(False, True, True, 0, 7, 0)`, `(False, True,
  True, 0, 9, 0)`, `(False, True, True, 0, 10, 1)`, and `(False, True, True, 0,
  12, 1)`. The exact required signatures are `(False, True, True, 1, 2n, 0)`.
  `qm/ACID_MECHANISMS.md` carries all geometry/status/log hashes and the bounded
  claim: these four deterministic seeds fail; this is not a proof over every
  microsolvation conformer. No Al TS run can make a reportable matched pair
  while all Si counts are red, so no Al barrier or ordering was fabricated.
  Two cold reviews found four false-green/state-corruption classes: aggregate
  shell counts hid OH-/H3O+ swaps, framework H atoms were untracked, shell O
  atoms were omitted from basin equivalence, and reactant-only runs could
  quarantine completed outputs/trust malformed receipts. All were fixed with
  physical-H occupancy gates, full solvent heavy-core equivalence, strict
  geometry/settings-bound receipt validation, exact archived-coordinate reuse,
  and non-destructive `reactant_status.json`. Fresh whole-suite proof is **226
  passed**, Ruff check, changed-file Ruff format, CLI refusal, and diff check.
  Email pipelines and Honcho inference were isolated under dead-man restarts for
  every GPU run and are restored healthy. Next move is a finite matched
  proton-relay conformer ensemble (`A1c-acid-microsolvation-conformers`); do
  not rerun water counts or constrain Obr-H.
- 2026-08-23 19:04 PDT — Victor selected matched microsolvation for **both**
  models. Continuation branch `agents/A1b-acid-microsolvation` was created from
  current `origin/main`; the historical PR #71 worktree was copied into this
  active campaign tree and removed with the mandatory review-ready teardown
  (`ok=true`, exact PR/local/remote head `6a103e4`, 3.24 GB reclaimed). The
  driver now accepts 3--6 total explicit waters with checkpoint-isolated run
  namespaces, deterministic non-colliding shell seeds, matched hydrated-H3O+
  fragment thermochemistry, and solvent-aware proton-relay basin gates. TDD
  first failed at import; the implementation then passed 63 focused and 220
  full no-cache tests, whole-tree Ruff check, changed-file Ruff format, CLI
  help, neutral-mode refusal, and diff check. The pre-existing whole-tree
  formatter debt is only `phase2_ladder.py` / `ts2_directed_oss.py` and is
  already owned by mission-control TASK-196. No production DFT campaign has
  launched yet; first live gate is a matched four-water Si/Al unconstrained
  reactant-minimum proof.
- 2026-08-22 12:29 PDT — Live GPU campaign closed the one-water `al-acid` model as a fully gated **sequential associative** mechanism: addition ΔG‡(298) = **82.193 kJ/mol = 19.645 kcal/mol** (205.95i cm⁻¹), cleavage local ΔG‡ = **8.047 kJ/mol = 1.923 kcal/mol** (97.46i cm⁻¹), and both quick-IRC endpoint pairs match the exact acid connectivity/proton-count basins. Mechanism-v2 `results.json` + integrity-checked `store.sqlite` are complete under `runs/phase1/al-acid-preprotonated-v2-b3lyp-def2-svp-flank/`; durable hashes and mechanism receipts are in `qm/ACID_MECHANISMS.md`. The matching one-water `si-acid` pre-equilibrium is **not an unconstrained minimum**: gas phase, a 4.5 Å water-placement probe, and aqueous PCM all return the bridge proton to H3O+ (`(False, True, True, 0, 3)`), while SMD/GPU optimization is numerically pathological. Card is blocked on a 3–6-water microsolvated Si-acid reactant; no constrained fake barrier or unsupported acid ordering is reported.
- 2026-08-21 18:20 PDT — Mandatory shared-GPU preflight still blocks the acid campaigns. Dissertation card A3's live bounded `ts2_directed_oss.py` process owns the RTX 4090 at 100% (PID 2891394, 4,244 seconds elapsed, 9,960 MiB process / 11,930 MiB total GPU memory); its durable `oss-neutral-n4-s2-pilot.log` is 3,032,960 bytes. No competing si-acid/al-acid DFT run launched. The implementation worktree remains clean before this bookkeeping checkpoint; acceptance is still red on both archived barrier results, so the mission-control pointer is parked until 2026-08-23 while this atomic board claim stays active.
- 2026-08-21 00:31 PDT — CPU implementation checkpoint complete. Added `protonated_bridge_complex()` (H3O+ pre-equilibrium → protonated intact bridge + residual H2O), a mechanism-versioned acid checkpoint namespace, acid-only approach→guarded product→CI-NEB routing, full Si/Al bridge/speciation gates for reactant/product/IRC, fail-closed cached-product validation, and crash-resumable product reuse. Cold review found four false-green/resume classes; all four were fixed with full 19-atom Si/Al regressions. Fresh exact-branch gates: **177 passed**, Ruff check, targeted Ruff format, CLI-help smoke, and diff checks green. The live A3 ladder pilot still exclusively occupies the RTX 4090 (PID 124433, 11,693 s elapsed, 100% utilization, 17,686 MiB process allocation at preflight); no competing acid DFT campaign launched. Code is ready and durable; si-acid/al-acid barriers remain pending the shared GPU lane.
- 2026-08-21 00:07 PDT — Claimed atomically by `hermes-workstation`; pushed branch `agents/A1b-acid-mechanisms` from current `origin/main` `bd8af3b8b60d5b0df1e80498443090aa208316a9` and created the congruent dedicated worktree. GPU/resource preflight and code archaeology are next.
- 2026-08-18 — card created (fable).
