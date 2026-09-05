# A2d-oss-neutral-n4-production-energetics — re-tier embedded pilot route

- status: blocked
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaign)
- depends: A2c-al-acid-production-energetics
- blocked-on: A2c
- claimed-by:

## Objective

Reconstruct and re-tier the documented `oss-neutral-n4-s2` embedded Si–O–Si
pilot from B3LYP/def2-SVP/DF onto A2's settled production protocol. Preserve
the accepted sequential mechanism and PHVA/frozen-shell contract. This is the
historical pilot only; it does not consume or relabel current A3/Osa recovery
outputs.

The historical ignored stationary-point directory and a hash-pinned external
archive are not currently available. Reconstruct the route from the committed
pilot builder and the mechanism/PHVA contract in
`docs/program/cards/A3-barrier-ladder.md`; do not pretend that a missing ignored
artifact was banked. Any later-discovered source may be consumed only after its
identity and provenance are independently hash-bound.

## Execution

1. Rebuild and hash-bind the pilot reactant, associative intermediate,
   hydrolyzed product, frozen shell, and PHVA identity from committed inputs.
   If an archived source appears, compare it to the reconstruction and refuse
   unexplained identity drift.
2. Extend the resume-safe A2 route contract to periodic embedded/frozen-shell
   geometry and PHVA semantics without weakening full-system contracts.
3. Optimize and verify the required stationary points at exact r2SCAN-3c
   (D4+gCP), preserving the frozen shell. Require accepted PHVA indices and
   post-optimizer connectivity/proton/frozen-coordinate gates.
4. Locate and verify both addition TS1 and bridge-cleavage TS2 at production
   geometry. A bounded scan/minimax proof may classify addition as barrierless,
   but the SVP assertion that TS1 lies below TS2 is not inherited. Replace
   quick-IRC evidence for every accepted soft-mode saddle with full
   bidirectional Gonzalez–Schlegel IRC and exact adjacent-basin identity.
5. Run wB97M-V/def2-TZVPD+SMD and B3LYP-D4/def2-TZVPD+SMD single points on
   accepted roles; emit the profile maximum and method shifts to an atomic
   receipt and provenance store.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one branch/worktree/PR.
- GPU-first; OMP/MKL/OpenBLAS <=16; tee long runs; hold the canonical GPU lease.
- This route is materially larger than free dimers. Use checkpointed bounded
  stages and fail closed on 4090 OOM. Do not rent H200/B200 capacity without a
  separately authorized spend gate; B300 remains excluded.
- Never consume A3's currently invalidated or incomplete Osa outputs as pilot
  truth. Computational failure remains `incomplete-*`.

## Acceptance

- Deterministic reconstruction is hash-valid and identifies every role,
  frozen index, method setting, PHVA projection, and route receipt used.
- Regressions cover frozen-shell drift, PHVA/full-Hessian distinction, source
  drift, endpoint identity, checkpoint drift, stale-output quarantine, and
  atomic terminal status.
- A terminal receipt has no `running` state; `store.sqlite` passes
  `PRAGMA integrity_check` with complete method/geometry/frozen-shell provenance.
- Production ΔG‡ plus SVP, r2SCAN-3c, wB97M-V, and B3LYP-D4 electronic barriers
  and shifts are reported only after TS1 and TS2 are independently verified or
  the addition segment has a bounded barrierless/minimax proof. Otherwise the
  terminal verdict is `incomplete-*` and no production ΔG‡ is promoted.
- The kaolinite Petra fragment and `qm/CALCULATIONS.md` are updated only from
  verified re-tiered values; full QM/Ruff/format/CLI/diff/teardown gates pass.
- The same PR sets parent A2's card and PLAN row to `ready`, clears its
  `blocked-on` field, and records the transition; parent closeout must not
  require a manual board edit.

## Progress

- 2026-09-05 12:25 PDT (hermes-custom-build-001; profile=workstation) — Split
  from A2 because the 63-atom embedded pilot has independent frozen-shell/PHVA,
  deterministic-reconstruction, runtime, and 4090-memory hazards. No current
  A3/Osa output is accepted as a substitute for this documented pilot. This
  card follows A2c so the embedded extension lands on one reviewed multi-step
  route contract.
