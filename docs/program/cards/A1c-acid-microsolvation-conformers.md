# A1c-acid-microsolvation-conformers — finite matched proton-relay ensemble

- status: blocked
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaign)
- depends: A1b microsolvation driver/receipt hardening merged
- claimed-by: hermes-workstation
- blocked-on: A1d-acid-3w-bridge-convergence

## Objective

Resolve A1b's remaining scientific ambiguity with a **finite, matched conformer
ensemble**, not another single-shell or water-count retry.

1. Add deterministic proton-relay seed families for 3--6 total explicit waters.
   At minimum cover: bridge-donor chain, attacker-centered ring, split
   bridge/attacker solvation, and one compact cyclic relay. Atom order and the
   physical-H ownership contract must remain stable across every family.
2. Screen Si--O--Si first at B3LYP/def2-SVP/DF. For every seed, run bounded
   unconstrained optimization, exact connectivity/physical-occupancy gates, and
   a frequency calculation only after the protonated bridge survives geometry
   optimization. Deduplicate converged basins by heavy-atom RMSD + energy.
3. If at least one zero-imaginary protonated Si bridge survives, run the exact
   same topology/water count for Si--O--Al. Only a topology that yields valid
   minima for **both** models may enter A1b's existing addition/cleavage ladder.
4. If every bounded Si seed deprotonates, publish a model-valid NO-GO for the
   pre-equilibrated-bridge assumption and file the executable redesign card for
   a concerted microsolvated hydronium/proton-relay acid mechanism. Do not hold
   Obr--H by constraint and do not emit an unmatched Al ordering.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet POLICY.md; one branch/worktree/PR.
- GPU-first, OMP/MKL/OpenBLAS <=16, every long run teed, email extraction and
  Honcho inference isolated under a bounded dead-man restart.
- Checkpoint namespace must include water count, conformer family, method, and
  gate version. Rejected/stale outputs are quarantined; canonical success is
  emitted only after all route/minimum/IRC gates.
- Use the exact physical-H occupancy and non-destructive reactant receipt
  contracts landed by A1b. A count-aggregate gate is not acceptable.

## Acceptance

- At least four documented deterministic topologies at each water count 3--6,
  with non-collision/unit tests and exact atom-index contracts.
- Every Si seed has a durable converged/failed/blocked receipt and teed log;
  valid minima have zero imaginary modes and deduplicated basin identity.
- Either: one matched Si/Al topology closes both gated barriers and reports the
  acid ordering; or: every bounded Si seed is a receipted NO-GO and a concrete
  concerted-relay redesign card is filed.
- Full QM suite, Ruff check, changed-file format, CLI smoke, and diff check pass.
- Card Progress/Result plus STATUS.md are updated in the same PR.

## Progress

- 2026-08-24 05:59 PDT — **DEVIATION / SCIENTIFIC BLOCKER:** the finite
  production campaign completed all 16 Si seeds, but its exact final verdict is
  `incomplete-si-screen`, not a matched minimum and not a model-valid NO-GO.
  Fifteen B3LYP/def2-SVP/DF optimizations converged and were rejected by the
  exact protonated-bridge/physical-H occupancy gates. The sole exception,
  `si-acid:3w:bridge-donor-chain`, exhausted the 160-step production bound and
  remains a hash-receipted `failed` endpoint. Summary: accepted 0, rejected 15,
  failed 1, blocked 0; no Si basin, Al screen, Hessian, matched barrier, or acid
  ordering was emitted. The ignored run tree was copied byte-for-byte before
  teardown risk to
  `/mnt/data/vsletten/dissertation-data/task197-a1c-acid-conformers-20260824/acid-microsolvation-ensemble-v2-g1-b3lyp-def2-svp`;
  manifest SHA-256 is
  `471a099c01c3fd8107409fd087e033bdd5da167cc6d41f09475a6a661cae5ac6`
  and log SHA-256 is
  `80729b9ed1dd61728f4959f29e3adb8151341a5764a779559ded1610005af312`.
  Email pipelines and both Honcho inference containers were restored immediately
  after process exit; the 13:37 dead-man was cancelled. Executable card
  `A1d-acid-3w-bridge-convergence` owns the one-seed fresh-optimizer refinement
  and strict final-verdict closeout. A1c remains blocked rather than laundering
  optimizer exhaustion into the redesign-triggering all-seed NO-GO.
- 2026-08-24 04:21 PDT — bounded worker handoff from the live production
  campaign. PID `4160112` remains healthy on exact pushed code; the atomic
  manifest now has five terminal Si records (four connectivity/occupancy
  rejections and one 160-step optimizer failure) while
  `si-acid:4w:attacker-centered-ring` is running. No Al screen or scientific
  verdict has started. The durable run root/log and 13:37 PDT dead-man remain
  unchanged; do not launch a competing campaign. The current partial population
  cannot support the final matched-minimum or model-valid NO-GO claim.
- 2026-08-24 — production campaign launched from pushed head `6e4cd646` with
  exact run root `qm/runs/phase1/acid-microsolvation-ensemble-v2-g1-b3lyp-
  def2-svp`, durable `ensemble.log`, advisory HF bound 40, load-bearing B3LYP
  bound 160, OMP/MKL/OpenBLAS 16, nice 10, and process-scoped 16 GiB memlock.
  OS PID `4160112` is live under the GPU-isolated lane; both email extraction
  pipelines and Honcho inference containers are stopped, with bounded systemd
  dead-man `a1c-qm-deadman-restart.timer` due 2026-08-24 13:37 PDT. The
  manifest is atomic and every completed seed is hash-receipted, so the next
  worker must inspect this PID/manifest/log and never launch a competing copy.
- 2026-08-24 — campaign preflight correction: the first GPU-isolated launch was
  deliberately stopped during seed 1's HF pre-relaxation at step 76/100, before
  any production B3LYP result or accepted receipt. geomeTRIC convergence remains
  mandatory for every production endpoint, while the cheap HF stage is capped
  at 40 advisory steps, exposes and records its convergence bit, and may hand
  its bounded last geometry to the production optimizer instead of burning the
  entire ensemble on a non-load-bearing precondition. Fresh
  293-test/Ruff/format/diff gates pass; the stopped
  seed's running checkpoint was quarantined before relaunch. A second launch
  proved the first production endpoint was close but unconverged at step 100
  (RMS gradient `6.018e-4`, max `1.663e-3`) and was stopped during seed 2;
  the runner now preserves every bounded production endpoint plus its explicit
  convergence bit before failing closed, so a larger bounded retry has durable
  evidence instead of terminal scrollback.
- 2026-08-24 — implementation checkpoint: four deterministic 3--6-water
  families now carry explicit connected donor/acceptor graphs, stable O/H/H
  indices, collision gates, and exact physical-H ownership. The finite Si-first
  ensemble journals every bounded optimization, requires real geomeTRIC
  convergence, validates Hessian receipts, deduplicates by heavy-atom RMSD plus
  finite energy, and only emits an Al/barrier handoff when both exact minima are
  hash-bound in one matched receipt. A cold adversarial review found and the
  branch fixed optimizer-exhaustion, NaN/weak-receipt, ambiguous-proton,
  disconnected-family, Al-failure-verdict, lowest-energy-canonical, and
  unmatched-barrier false greens. Pre-campaign verification: 293 QM tests,
  whole-tree Ruff, changed-file format, both CLI help smokes, and diff check
  pass. Implementation campaign checkpoint is committed/pushed next; production
  GPU screening has not started yet.
- 2026-08-23 — created from A1b's 3--6-water deterministic-shell NO-GO. The
  count sweep proved four seeds fail but explicitly did not claim the complete
  microsolvated potential-energy surface was exhausted.

## Result

- Implementation is complete and remotely durable on branch
  `agents/A1c-acid-microsolvation-conformers`: four deterministic connected
  proton-relay families for each water count 3--6; exact atom/H ownership;
  collision/connectivity/occupancy gates; bounded optimization and frequency
  sequencing; atomic hash-bound receipts; heavy-atom-RMSD plus energy basin
  deduplication; and fail-closed matched Si/Al barrier handoff.
- Production ran the exact finite 16-seed Si ensemble at
  B3LYP/def2-SVP/DF. The machine result is `15 rejected / 1 failed / 0
  accepted`; the single failure is optimizer exhaustion for the three-water
  bridge-donor chain. Therefore A1c is blocked on
  `A1d-acid-3w-bridge-convergence`. No Al minimum, barrier, acid ordering, or
  concerted-mechanism redesign claim was fabricated.
- Durable external evidence:
  `/mnt/data/vsletten/dissertation-data/task197-a1c-acid-conformers-20260824/acid-microsolvation-ensemble-v2-g1-b3lyp-def2-svp`
  (1.3 MiB; manifest SHA-256
  `471a099c01c3fd8107409fd087e033bdd5da167cc6d41f09475a6a661cae5ac6`;
  log SHA-256
  `80729b9ed1dd61728f4959f29e3adb8151341a5764a779559ded1610005af312`).
- Final verification and PR/teardown receipts are recorded in the associated
  mission-control TASK-197 Result.
