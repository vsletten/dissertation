# A2b1-wb97mv-reactant-scf-recovery — bounded reactant SCF adjudication

- status: done
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

- 2026-09-06 12:41 PDT — DONE with the predeclared terminal-failure outcome.
  The Hückel → 50-cycle damped/0.5-Eh-level-shifted Roothaan stage and its
  density-seeded 150-cycle unshifted CDIIS finalization both ended with
  `converged=false`; final orbital-gradient norms were respectively
  `81.57863386820702` and `4.388151693060163`, and every energy remains
  diagnostic-only. A cold independent verifier recomputed the terminal,
  run-status, attempt, driver, legacy-source/settings, geometry, and frequency
  hashes; confirmed the exact density handoff; reran all nine focused
  publication-boundary tests; and observed no unit, process, running record,
  canonical energy, `results.json`, or SQLite store. Terminal receipt SHA-256 is
  `49d1fd95e82a8a178b3a33c8c946516d3a3655e0a2d108cb3afccf4f3f2bca6d`.
  A2b remains blocked on terminal reactant SCF failure; no solver, method,
  geometry, threshold, or cycle-bound expansion is authorized. Final fresh gates
  are `642 passed`, whole-QM Ruff, changed-file Ruff format, worktree-PYTHONPATH
  CLI smoke, board/PLAN cross-check, and `git diff --check` green.
  (hermes-custom-build-001; profile=workstation)
- 2026-09-06 11:02 PDT — Launched the exact reactant gate as bounded transient unit `task294-a2b1-reactant-scf-recovery.service` (invocation `120d66393dcc47729ea73c9fe17ba079`): driver wall 7,200 s, `RuntimeMaxSec=7,800`, `MemoryMax=48G`, 16 threads, 18 GiB GPU pool, canonical GPU lease acquired by PID 33585, and `ExecStopPost` atomic failure backstop installed. Prelaunch gates: Ruff clean, format clean, full suite 641 passed / 1 skipped, and both independent review rejection classes reproduced then fixed with adversarial tests. Terminal blocker: `/mnt/data/vsletten/dissertation-data/task294-a2b1-wb97mv-reactant-scf-recovery-20260906/terminal-receipt.json`. (hermes-custom-build-001; profile=workstation)
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

## Result

- `incomplete-computational-failure`, independently verified. Attempt 1
  exhausted 50 damped/level-shifted Roothaan cycles; attempt 2 exhausted 150
  fresh unshifted CDIIS cycles from the exact attempt-1 density. Both actual
  convergence flags are false. Attempt receipt hashes are
  `c6a980a031020e4c02c638221463260139319e30858cc93015df3ab5306afe7b`
  and `b0020febb2a653eafd9a27dfe8b046e1386687729e8ff0e69b88d7069c68cf1b`.
- The exact wB97M-V/def2-TZVPD+SMD(water), charge `-1`, spin `0`, geometry,
  settings, source, and driver bindings rehash. No running record or process
  remains, and neither the TASK-290 nor TASK-294 evidence root contains a
  canonical reactant energy, `results.json`, `store.sqlite`, WAL, or SHM file.
- This is a finite computational failure, not a mechanism rejection. A2b and
  A2c remain blocked. The authorized recovery budget is exhausted; no further
  retry, ad hoc solver, method/geometry change, weakened convergence, or longer
  cycle bound is authorized by this card. (hermes-custom-build-001; profile=workstation)
