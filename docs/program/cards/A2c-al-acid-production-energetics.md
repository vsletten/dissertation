# A2c-al-acid-production-energetics — re-tier banked one-water acid route

- status: blocked
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaign)
- depends: A2b-al-neutral-production-energetics
- blocked-on: A2b
- claimed-by:

## Objective

Re-tier A1b's accepted one-water `al-acid` sequential associative route from
B3LYP/def2-SVP/DF onto A2's settled production protocol. This card covers the
banked Al route only. It does not authorize another Si-acid topology search;
A1i owns the one exact production-tier Si revisit after A2 closes.

Source evidence is immutable at
`/mnt/data/vsletten/dissertation-data/task177-a1b-acid-microsolvation-20260823/phase1/al-acid-preprotonated-v2-b3lyp-def2-svp-flank/`.
Bind checkpoints to its source store/geometry hashes and exact merged A2 method
identity.

## Execution

1. Extend the resume-safe production route contract to the accepted
   reactant→addition-TS→intermediate→cleavage-TS→product sequence.
2. Optimize all five roles on exact r2SCAN-3c (D4+gCP); use finite-difference
   composite Hessians and require minima/TS indices `0/1/0/1/0` above the
   documented imaginary floor.
3. Run full bidirectional Gonzalez–Schlegel IRC for both saddles. Require exact
   atom order, electronic state, physical-H ownership, and heavy connectivity
   against their adjacent accepted basins.
4. Run wB97M-V/def2-TZVPD+SMD and B3LYP-D4/def2-TZVPD+SMD single points on every
   accepted stationary point. Report the reactant-referenced profile maximum,
   rate-limiting step, and SVP→r2SCAN-3c→production shifts.
5. Emit one atomic terminal receipt and provenance-complete `store.sqlite`;
   archive ignored artifacts outside git and publish their hashes.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one branch/worktree/PR.
- GPU-first; OMP/MKL/OpenBLAS <=16; tee long runs; hold the canonical GPU lease.
- Preserve acid mechanism-v2 physical-proton ownership and route identity.
- Computational failure is `incomplete-*`, never a scientific rejection. Do
  not silently substitute an acid conformer or Si topology.
- No external compute spend without a separately authorized gate.

## Acceptance

- Regressions prove source/settings binding, both-saddle full-IRC identity,
  proton conservation/ownership, checkpoint drift rejection, stale-output
  quarantine, and atomic terminal status.
- The terminal receipt is hash-valid with no `running` record; `store.sqlite`
  passes `PRAGMA integrity_check` and records all method/geometry provenance.
- Production ΔG‡ plus SVP, r2SCAN-3c, wB97M-V, and B3LYP-D4 electronic barriers,
  rate-limiting step, and shifts are reported without extrapolating to Si-acid.
- `qm/ACID_MECHANISMS.md` and `qm/CALCULATIONS.md` carry the verified tier.
- The same PR sets A2d's card and PLAN row to `ready`, clears its `blocked-on`
  field, and records the transition; A2d must not require a manual board edit.
- Full QM suite, Ruff, changed-file format, CLI smoke, diff check, durable
  hashes, and review-ready teardown pass.

## Progress

- 2026-09-05 12:25 PDT (hermes-custom-build-001; profile=workstation) — Split
  from A2 after receipt-level inventory found a valid four-stationary-point SVP
  store and product geometry, but no r2SCAN-3c, production single-point, or full
  IRC receipts. This card follows A2b so it extends one merged route framework
  instead of creating a conflicting second implementation.
