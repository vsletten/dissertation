# A2b-al-neutral-production-energetics — re-tier banked Si–O–Al neutral route

- status: active
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaign)
- depends: A2a ✅ (settled production method contract)
- claimed-by: hermes-custom-build-001

## Objective

Re-tier the accepted TASK-168 one-water neutral Si–O–Al route from
B3LYP/def2-SVP/DF onto A2's settled production protocol. Preserve the accepted
sequential mechanism: quasi-barrierless uphill addition into the associative
intermediate, followed by the verified bridge-cleavage saddle. Do not resurrect
the rejected two-imaginary-mode addition candidate as a transition state.

Source evidence is immutable at
`/mnt/data/vsletten/dissertation-data/task168-al-neutral-20260823/`.
Bind every checkpoint to the source store/geometry hashes and the exact merged
A2 method identity.

## Execution

1. Generalize the resume-safe A2 production driver just enough to represent a
   sequential route with a barrierless addition segment and one cleavage TS.
2. Optimize the reactant, intermediate, released product, and cleavage saddle
   on exact r2SCAN-3c (D4+gCP) surfaces; use finite-difference composite
   Hessians and require indices `0/0/0/1` above the documented imaginary floor.
3. Run full bidirectional Gonzalez–Schlegel IRC from the cleavage saddle and
   require exact atom order, electronic state, physical-H ownership, and heavy
   connectivity to match the accepted intermediate/product basins.
4. Run wB97M-V/def2-TZVPD+SMD and B3LYP-D4/def2-TZVPD+SMD single points on all
   accepted stationary points. Compute the reactant-referenced profile maximum
   and SVP→r2SCAN-3c→production barrier shifts.
5. Write one atomic terminal receipt and provenance-complete `store.sqlite`.
   Archive ignored artifacts below `/mnt/data/vsletten/dissertation-data/` and
   publish hashes, not large outputs, to git.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one branch/worktree/PR.
- GPU-first; OMP/MKL/OpenBLAS <=16; tee long runs; hold the canonical GPU lease.
- Optimizer convergence, Hessian index, IRC basin identity, and mechanism
  acceptance are separate fail-closed gates.
- A finite computational failure is `incomplete-*`, never a scientific
  rejection. No external compute spend without a separately authorized gate.

## Acceptance

- Regression coverage proves source-hash/settings binding, rejected-addition
  refusal, checkpoint drift rejection, true-minimum/one-saddle indices,
  typed full-IRC endpoint identity, stale-output quarantine, and atomic status.
- A hash-valid terminal receipt has no `running` record and the output store
  passes `PRAGMA integrity_check` with exact method/geometry provenance.
- Production ΔG‡ and electronic barriers for SVP, r2SCAN-3c, wB97M-V, and
  B3LYP-D4 are reported with barrier shifts; no rejected saddle contributes.
- `qm/AL_NEUTRAL_MECHANISM.md` and `qm/CALCULATIONS.md` carry the verified tier.
- The same PR sets A2c's card and PLAN row to `ready`, clears its `blocked-on`
  field, and records the transition; A2c must not require a manual board edit.
- Full QM suite, whole-tree Ruff, changed-file format, CLI smoke, diff check,
  durable hashes, and review-ready teardown pass.

## Progress

- 2026-09-06 04:13 PDT (hermes-custom-build-001; profile=workstation) —
  BLOCKED only on the finite recovery receipt. Exact clean/pushed branch
  `agents/A2b-al-neutral-production-energetics@47d6e2a25697ce8cc3a90249969b7b34ea4c1fec`
  contains the SCF recovery and cache-refusal tests. Bounded transient unit
  `a2b-task290-production-recovery.service` invocation
  `90f3237590404a90aee28785cac91a51` is active/running with PID `3363143`, a
  36-hour ceiling, 40 GiB RAM + 8 GiB swap, 1600% CPU/nice-10, 16 numerical
  threads, and the canonical 18-GiB GPU lease. `run_status.json` is `running`,
  the old failure receipt is quarantined, and the new atomic terminal receipt is
  absent. On receipt arrival, resume adjudication must verify exact head/hashes,
  accepted cached stage identities, all converged production receipts, final
  profile/store/docs/A2c transition, then open the single PR and tear down.

- 2026-09-06 04:11 PDT (hermes-custom-build-001; profile=workstation) —
  Adjudicated the first physical receipt as an honest finite computational
  failure, not a mechanism rejection: all four r2SCAN-3c minima/saddle index
  gates and the typed full cleavage IRC passed, but the first production
  wB97M-V/def2-TZVPD+SMD single point exhausted 100 Newton cycles on the
  checkpointed reactant. Two cold read-only reviews rejected an identical replay
  and found that cached single-point receipts did not explicitly require
  `converged=true`. The shared A2 helper now rejects unproven caches and adds one
  bounded, receipted direct-DIIS attempt seeded from the failed Newton density;
  unconverged values remain unpublished. Focused A2/A2a/A2b verification is 29
  passed with focused Ruff/format and diff checks green. A fresh exact-head
  supervised resume is next; accepted geometry/frequency/IRC checkpoints are
  reused rather than recomputed.

- 2026-09-06 02:40 PDT (hermes-custom-build-001; profile=workstation) —
  Launched the physical A2b campaign as bounded user unit
  `a2b-task290-production.service` (invocation
  `b2b0e656685d43ebaf921f551234d213`; runtime 36 h; memory 40 GiB + 8 GiB
  swap; 16 threads; 18 GiB GPU pool). Read-back proves `active/running`, PID
  `3152230`, the canonical GPU lease owned by `a2b-al-neutral-production`, and
  `run_status.json`=`running`; terminal supervision receipt is
  `/mnt/data/vsletten/dissertation-data/task290-a2b-al-neutral-production-20260906/terminal-receipt.json`.
- 2026-09-06 02:36 PDT (hermes-custom-build-001; profile=workstation) —
  Implemented and adversarially fixture-tested the route-specific A2b front end
  over the shared resume-safe A2 checkpoints: exact TASK-168 artifact/job/
  mechanism binding, explicit rejected-addition refusal, typed basin + physical-H
  ownership + heavy-topology gates, four r2SCAN-3c stationary-point/index gates,
  full cleavage IRC, both production single-point tiers, atomic stale-output
  quarantine/terminal receipt, and provenance-correct external-store publication.
  Live source validation passes against the immutable TASK-168 archive; the full
  QM suite is green (`609 passed`); the executable acquires the canonical QI2
  GPU lease with an 18 GiB CuPy pool ceiling. Physical campaign launch is next.
- 2026-09-05 12:25 PDT (hermes-custom-build-001; profile=workstation) — Split
  from A2 after receipt-level inventory proved that only si-neutral currently
  satisfies the production/IRC/store gates. TASK-168 remains a valid banked
  SVP mechanism and immutable source; this card owns its bounded re-tier.
