# A2b-al-neutral-production-energetics — re-tier banked Si–O–Al neutral route

- status: ready
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaign)
- depends: A2a ✅ (settled production method contract)
- claimed-by:

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

- 2026-09-05 12:25 PDT (hermes-custom-build-001; profile=workstation) — Split
  from A2 after receipt-level inventory proved that only si-neutral currently
  satisfies the production/IRC/store gates. TASK-168 remains a valid banked
  SVP mechanism and immutable source; this card owns its bounded re-tier.
