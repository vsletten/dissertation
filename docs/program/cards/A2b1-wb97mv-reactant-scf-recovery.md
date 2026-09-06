# A2b1-wb97mv-reactant-scf-recovery — bounded reactant SCF adjudication

- status: active
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaign)
- depends: A2a-si-neutral-production-path-rebuild
- claimed-by: hermes-custom-build-001

## Objective

Recover or terminally adjudicate the exact A2b Al-neutral reactant
wB97M-V/def2-TZVPD+SMD(water) single point after both the 100-cycle Newton route
and its one 150-cycle density-seeded direct-DIIS fallback failed. This is a
reactant-only gate. Do not run the intermediate, product, saddle, barrier,
thermochemistry, or store publication stages until this exact reactant SCF has a
hash-valid converged receipt.

Reuse the accepted A2b r2SCAN-3c geometry/frequency/IRC checkpoints only after
exact source, settings, geometry, and validator revalidation. Preserve the A2
production method identity: functional, basis, SMD solvent, charge, spin,
geometry, density fitting, grids, and convergence acceptance do not change.

## Execution

1. Persist one atomic receipt for every SCF attempt, including failures. Bind
   exact geometry/settings/code fingerprints, solver configuration, cycle bound,
   actual convergence flag, available energy/residual/cycle diagnostics, and
   initial/final density or MO provenance. Unconverged values are diagnostic-only.
2. Add an adversarial both-attempts-fail regression proving no canonical energy,
   `results.json`, or `store.sqlite` can survive or be emitted.
3. Predeclare and independently review one scientifically motivated recovery
   route materially different from the exhausted Newton→direct-DIIS sequence.
   It must have a hard cycle/time bound and an explicit terminal failure outcome.
4. Run only the exact reactant gate under the canonical GPU lease and bounded
   systemd supervision. Rehash the terminal receipt and independently verify the
   convergence flag, cache identity, diagnostics, and canonical-output boundary.
5. On success, update A2b and its PLAN row to `ready`; A2b then owns the remaining
   seven production single points, profile/barrier/thermochemistry/store/docs, and
   A2c transition. On terminal failure, record `incomplete-computational-failure`
   and authorize no ad hoc solver, method, geometry, or cycle-bound expansion.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one branch/worktree/PR.
- GPU-first; OMP/MKL/OpenBLAS <=16; tee the long run; hold the canonical QI2 GPU
  lease; run compute longer than 30 minutes only as a bounded transient unit with
  an atomic terminal receipt.
- Preserve exact electronic state and A2 method identity. No weakened convergence,
  accepted unconverged energy, silent geometry change, or external compute spend.
- The two failed TASK-290 routes are immutable evidence. An identical replay is
  prohibited.

## Progress

- 2026-09-06 10:42 PDT — Claimed atomically from `origin/main@dc194e7ae9eecdc0d5c39e904f6f8511e3a82c4a`; independent scientific review approved the predeclared `a2b1-huckel-damped-level-shifted-roothaan-to-cdiis-v1` contract only. Implemented hash/settings/geometry/driver-bound atomic per-attempt receipts, strict cache revalidation, a fail-closed timeout backstop, and adversarial both-fail/tamper/timeout tests. A cold implementation review rejected launch for semantic receipt gaps and pre-start timeout behavior; those findings are fixed and focused verification is 17 passed. (hermes-custom-build-001; profile=workstation)

## Acceptance

### Verify

- Regression coverage proves atomic per-attempt success/failure receipts,
  geometry/settings/code binding, rejection of unconverged caches, stale-output
  quarantine, and zero canonical publication when all attempts fail.
- A cold independent verifier rehashes the selected recovery contract and the
  terminal evidence, checks the actual convergence flag and diagnostics, and
  confirms the canonical-output boundary from the original A2b complaint.
- Success requires a hash-valid converged reactant wB97M-V receipt with no running
  record and an exact A2b `ready` transition. Failure requires a hash-valid
  `incomplete-computational-failure` receipt, A2b still blocked, and no further
  retry authorization hidden in prose.
- Full QM suite, whole-tree Ruff, changed-file format, CLI smoke, diff check,
  durable hashes, and review-ready teardown pass.

## Progress

- 2026-09-06 05:02 PDT (hermes-custom-build-001; profile=workstation) — Split
  from A2b after a cold independent review verified that its one bounded
  Newton-to-direct-DIIS recovery failed before the first production energy and
  did not separately persist attempt diagnostics. This card owns one bounded,
  materially different reactant-only hypothesis; A2b and A2c remain blocked.
