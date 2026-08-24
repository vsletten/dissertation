# A1d-acid-3w-bridge-convergence — close the sole inconclusive A1c seed

- status: done
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaign)
- depends: A1c finite-ensemble implementation and archived production receipts
- claimed-by: hermes-workstation

## Objective

Resolve the only inconclusive member of A1c's finite 3--6-water Si-first
ensemble: `si-acid:3w:bridge-donor-chain`. The exact B3LYP/def2-SVP/DF run
reached its 160-step production bound without geomeTRIC convergence; the other
15 deterministic seeds converged and were rejected by the existing exact
connectivity/physical-H occupancy gates. This card must preserve those 15
terminal receipts, refine only the failed seed through a reproducible bounded
path, and let A1c emit a scientifically valid final verdict.

1. Recover and independently hash-verify the archived A1c run at
   `/mnt/data/vsletten/dissertation-data/task197-a1c-acid-conformers-20260824/acid-microsolvation-ensemble-v2-g1-b3lyp-def2-svp`.
   Its manifest SHA-256 is
   `471a099c01c3fd8107409fd087e033bdd5da167cc6d41f09475a6a661cae5ac6`.
2. Add a tested, explicit failed-endpoint refinement mode to
   `acid_microsolvation_ensemble.py`. It may operate only on a hash-valid failed
   record whose reason is optimizer exhaustion. Preserve the original seed,
   preoptimized endpoint, production endpoint, and manifest as immutable
   provenance; use a new attempt namespace and receipt. Restart from the exact
   failed production endpoint with a fresh optimizer/Hessian, the same
   B3LYP/def2-SVP/DF settings, and one additional predeclared bound of at most
   160 steps. Do not rerun the 15 conclusive seeds and do not silently relabel
   their original 160-step receipts.
3. If the refinement converges, apply the same connectivity, exact physical-H
   ownership, solvent occupancy, minimum-pair-distance, finite-energy, and
   zero-imaginary-mode gates as every other A1c seed. Only an accepted Si
   minimum may trigger the exact matched Si--O--Al topology; no unmatched Al
   result may be published.
4. Regenerate the atomic summary from the original plus refinement receipts.
   A model-valid pre-equilibrated-bridge NO-GO requires all 16 Si seeds to have
   conclusive rejected receipts. A matched handoff requires hash-bound accepted
   Si and Al minima. Any remaining optimizer failure stays `incomplete-*`; file
   the next executable mechanism/convergence card rather than laundering it
   into a NO-GO.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one branch/worktree/PR.
- GPU-first; OMP/MKL/OpenBLAS <=16; tee every long run. Isolate email extraction
  and Honcho inference under a bounded dead-man and restore them immediately
  when the run exits.
- Do not constrain Obr--H, weaken convergence, infer a minimum from endpoint
  geometry, alter the 15 conclusive receipts, or launch the full ensemble again.
- The implementation worktree is disposable at PR creation. Archive all ignored
  campaign artifacts under `/mnt/data/vsletten/dissertation-data/` and verify
  hashes before mandatory review-ready teardown.

## Acceptance

- Regression tests force the exact failed-record continuation route and prove:
  settings/artifact drift refusal; fresh optimizer state; immutable original
  receipts; strict one-seed scope; bounded retry; and the same downstream
  minimum/occupancy verdict gates as an ordinary seed.
- The archived A1c manifest contains 16 terminal Si records plus a hash-bound
  refinement receipt for `si-acid:3w:bridge-donor-chain`; no `running` record
  remains and the final summary agrees mechanically with the record population.
- Either a matched Si/Al minimum is ready for A1b's barrier ladder, or all 16 Si
  seeds are conclusive rejections and A1c publishes the model-valid NO-GO plus
  an executable concerted hydronium/proton-relay redesign card. If the bounded
  refinement remains inconclusive, this card reports that honestly and files
  the exact next executable card.
- Full QM suite, whole-tree Ruff, changed-file format, CLI smoke, and diff check
  pass. Card/PLAN/STATUS bookkeeping and durable archive hashes land in the same
  PR; worktree teardown receipt is recorded.

## Progress

- 2026-08-24 10:15 PDT — **DONE / model-valid NO-GO.** The one card-pinned
  B3LYP/def2-SVP/DF refinement restarted from the byte-identical failed
  production endpoint with a fresh geomeTRIC state and converged after 54
  additional steps. The endpoint is structurally finite and non-colliding but
  fails the exact protonated-bridge basin: signature `(False, True, True, 0, 7,
  0)` versus required `(False, True, True, 1, 6, 0)`. The bridge proton again
  moved into the solvent network. The effective finite population is now
  **16 rejected / 0 failed / 0 accepted / 0 blocked**, so A1c's
  pre-equilibrated protonated-bridge model is conclusively NO-GO. No Al screen,
  unmatched barrier, frequency, or acid ordering was emitted. Executable card
  `A1e-acid-concerted-hydronium-relay` owns the concerted mechanism redesign.
- 2026-08-24 09:41 PDT — implementation preflight complete on atomic branch
  `agents/A1d-acid-3w-bridge-convergence`. Independently re-hashed the archived
  A1c manifest (`471a099c01c3fd8107409fd087e033bdd5da167cc6d41f09475a6a661cae5ac6`),
  ensemble log (`80729b9ed1dd61728f4959f29e3adb8151341a5764a779559ded1610005af312`),
  and all 48 referenced seed artifacts: 16 terminal Si records, exactly 15
  rejected plus the named optimizer-exhausted failure, zero drift and zero
  running records. **DEVIATION:** the source manifest cannot both retain its
  pinned SHA-256 and be rewritten to contain the refinement. A1d will therefore
  emit a separate hash-bound resolution manifest with an effective 16-record
  overlay (15 immutable A1c references plus one refinement receipt), leaving
  every original byte untouched. A cold independent design audit agreed this is
  the only fail-closed interpretation; implementation and production refinement
  have not started yet.
- 2026-08-24 — created from A1c's exact production result: 15/16 Si seeds
  converged and failed the physical protonated-bridge basin, while only the
  three-water bridge-donor chain exhausted the 160-step optimizer bound. No Si
  minimum, Al screen, matched barrier, or model-valid NO-GO exists yet.

## Result

- Added an explicit failed-endpoint refinement mode to
  `qm/scripts/acid_microsolvation_ensemble.py`. Production authority is pinned
  in code to the exact A1c manifest/log/seed; the attempt namespace and log are
  deterministic, alternate step budgets and log overrides are refused, source
  artifacts are re-hashed before publication, and ordinary/refinement routes
  share the same exact geometry/connectivity/physical-H/solvent/frequency gates.
  The immutable A1c archive remains byte-for-byte fixed.
- Durable resolution root:
  `/mnt/data/vsletten/dissertation-data/task199-a1d-acid-3w-bridge-20260824/acid-3w-bridge-refinement-b3lyp-def2-svp`.
  Final manifest SHA-256 `6b02432f4920b2ff8406de5cd22d379ebc54ceed66800b1db8149efff7df62cf`;
  evidence manifest `e6806c9c59eabcd6abb470f268f1a6424f51475fdaf51de34e9b81b7751574ac`;
  attempt receipt `52c52ba7024047f6b00b282f132d43086c21d4be5e50644af783b83c0cd4a463`;
  teed log `15bc3a64abb8bafac34b96eb57a9ea8a21d0e791489e51cc4c7d854f5a746ac9`.
  The copied input hash is the original endpoint hash
  `ba58f9b9cf9c3b44c977c93dcfccd9bef0c526cc46287315e1c38c5c519b3421`;
  refined endpoint hash is
  `b06774630adf1252013fac45f94df58dad15def05f146c876685b0235e68d094`.
- Verification after the final review fix: **308 full QM tests passed** with
  cache disabled; whole-tree Ruff check, changed-file Ruff format, pinned CLI
  smoke, and diff check passed. A cold adversarial review found six authority,
  duplicate-budget, provenance-log, finalization, premature-success, and
  unmatched-Al verdict defects; all six were fixed with regressions before the
  production attempt.
- GPU isolation closed cleanly: both email pipelines are active; Honcho deriver
  is running; Honcho API is healthy and `/health` returns `{"status":"ok"}`;
  the four-hour dead-man was cancelled after immediate restoration. Source
  manifest/log hashes remain exactly `471a099c…` / `80729b9e…`.
