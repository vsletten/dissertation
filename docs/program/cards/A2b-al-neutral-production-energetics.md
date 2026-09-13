# A2b-al-neutral-production-energetics — banked Si–O–Al neutral survey route

- status: blocked
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaign)
- depends: A9-approximate-rate-closure
- claimed-by: hermes-custom-build-001
- blocked-on: A9

## Objective

Bank the accepted neutral Si–O–Al survey result for the approximate-rate loop:
32.2 kcal/mol at r2SCAN-3c, with the TASK-168 one-water sequential mechanism and
its provenance. Preserve the quasi-barrierless uphill addition into the
associative intermediate followed by the verified bridge-cleavage saddle. Do
not resurrect the rejected two-imaginary-mode addition candidate, and do not
retry the exhausted wB97M-V reactant SCF routes.

Source evidence is immutable at
`/mnt/data/vsletten/dissertation-data/task168-al-neutral-20260823/`.
Use the immutable source store/geometry hashes to support the banked survey
value. A9's sensitivity verdict decides whether a later, separately scoped
published-family replication is warranted.

## Execution

1. Verify the existing r2SCAN-3c value, units, mechanism identity, and source
   hashes from the accepted evidence.
2. Record 32.2 kcal/mol as the banked A9 input with survey-tier provenance.
3. Make no new QM call on this card. If A9 ranks this family sensitive, a
   separate card may replicate the closest published family within POLICY v16's
   four-hour-per-unit envelope.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one branch/worktree/PR.
- No QM run, GPU lease, external compute spend, or SCF retry on this card.
- Historical branch, service, and artifact names retain their literal IDs; they
  are provenance, not authority for a higher-tier restart.

## Acceptance

- The 32.2 kcal/mol r2SCAN-3c value is recorded with exact source provenance and
  identified as survey tier.
- The accepted sequential mechanism remains unchanged and rejected saddles do
  not contribute.
- A9 consumes the banked value; this card remains blocked until A9's sensitivity
  verdict determines whether any published-family replication is worth running.
- No new higher-tier energy, SCF retry, or workstation quality claim is made.

## Progress

- 2026-09-13 12:04 PDT (hermes-custom-build-001; profile=workstation) — POLICY v16/A9 AMENDMENT: the accepted 32.2 kcal/mol r2SCAN-3c result is the banked survey-tier value for A9. The exhausted wB97M-V routes are historical evidence rather than a blocker to another workstation tier; A2b now waits only on A9's sensitivity verdict and authorizes no QM retry.

- 2026-09-06 12:41 PDT (hermes-custom-build-001; profile=workstation) — A2b1
  exhausted its independently reviewed, materially different Hückel → damped
  level-shifted Roothaan → fresh CDIIS route with both actual convergence flags
  false. The hash-valid terminal outcome is
  `incomplete-computational-failure`; exact density provenance and input
  identity rehash, the unit/process/running record are gone, and no canonical
  energy/result/store exists. A2b therefore remains terminally BLOCKED rather
  than becoming ready when A2b1 closes. This is not a mechanism rejection, and
  no additional solver, method, geometry, threshold, or cycle budget is
  authorized.
- 2026-09-06 05:02 PDT (hermes-custom-build-001; profile=workstation) —
  BLOCKED on READY `A2b1-wb97mv-reactant-scf-recovery` after the sole bounded
  Newton-to-direct-DIIS recovery also failed at the first wB97M-V reactant
  single point. Exact terminal receipt SHA-256 is
  `c4233567cc3a88a23683514531e89af43c1fe6743a020638bc4fb20281d2fecd`;
  its driver and run-status hashes re-match, the unit is terminal, and no
  canonical energy, result, barrier, store, or documentation value exists.
  DEVIATION: cold review found that the aggregate terminal identified the
  `newton-then-direct-diis` route but did not separately persist per-attempt
  diagnostics. A2b1 owns that evidence gap plus one materially different,
  reactant-only SCF hypothesis. Fresh closeout verification is `611 passed`
  with whole-QM Ruff, Ruff format, CLI smoke, and `git diff --check` green. No
  identical replay is authorized.

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
  bounded direct-DIIS attempt seeded from the failed Newton density; the
  aggregate terminal identifies that route, but no separate per-attempt receipt
  exists, and unconverged values remain unpublished. Focused A2/A2a/A2b verification is 29
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

## Result

- The route infrastructure and regression suite are review-ready, but A2b's
  scientific acceptance is incomplete. Both bounded wB97M-V reactant SCF routes
  failed before publishing a production energy: the original Newton route and
  its one direct-DIIS fallback seeded from the failed density. The current
  terminal receipt is hash-valid, has no running record, and reports
  `incomplete-computational-failure`; this is not a mechanism rejection.
- All four r2SCAN-3c stationary-point/index checkpoints and the typed full IRC
  remain reusable only after exact validator re-execution. They are not a
  production profile. `results.json`, `store.sqlite`, all production energy
  receipts, production barriers, thermochemistry, `qm/AL_NEUTRAL_MECHANISM.md`,
  and `qm/CALCULATIONS.md` updates remain absent.
- READY card `A2b1-wb97mv-reactant-scf-recovery` owns durable per-attempt failure
  evidence and one predeclared reactant-only recovery hypothesis. A2c remains
  blocked. A2b may resume only after A2b1 produces a hash-valid converged
  reactant receipt; another unchanged TASK-290 launch is prohibited.
- A2b1 has now closed with an independently verified finite failure: its two
  additional attempt receipts are hash-valid but unconverged, and no canonical
  production artifact exists. A2b remains blocked on terminal reactant SCF
  failure; this PR authorizes no further computational retry or A2c transition.
