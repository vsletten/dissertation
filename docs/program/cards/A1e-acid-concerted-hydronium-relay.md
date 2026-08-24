# A1e-acid-concerted-hydronium-relay — retire the failed pre-equilibrium

- status: done
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaign)
- depends: A1c/A1d model-valid NO-GO for the pre-equilibrated protonated bridge
- claimed-by: hermes-workstation

## Objective

Replace A1b's retired "protonate Obr, then attack" assumption with a finite,
matched **concerted microsolvated hydronium/proton-relay** mechanism. The
reactant must be a genuine unconstrained minimum with hydronium in the solvent
network and an intact unprotonated bridge; proton transfer toward Obr must occur
along the same elementary path as water attack and Si--Obr weakening/cleavage.

1. Add a mechanism-versioned concerted builder/driver that preserves the A1c
   physical-H ownership and atom-order contracts. Construct exact reactant and
   hydrolyzed-product endpoints without holding Obr--H, Si--Ow, or Si--Obr by a
   production constraint.
2. Keep the first campaign finite: Si--O--Si at total waters 3 and 4, for exactly
   `bridge-donor-chain` and `compact-cyclic-relay` (four deterministic paths).
   Persist and hash every endpoint, pre-relaxed band image, climbing band,
   saddle, Hessian, IRC endpoint, log, and terminal record. Bounds and checkpoint
   identities include topology, water count, method, and mechanism/gate version.
3. Use double-ended CI-NEB (persist every pre-relaxed image before climb), refine
   the in-channel crest with a fresh saddle search, and require exactly one
   imaginary mode whose displacement couples proton relay, Si--Ow approach, and
   Si--Obr cleavage. Full bidirectional IRC must land in the exact typed
   reactant/product basins; a water wag, proton shuffle without attack, or
   escaped cleavage channel is a rejection.
4. Only after a Si path passes every gate may the exact same water count/family
   run on Si--O--Al. A barrier/order handoff requires accepted, hash-bound Si and
   Al saddles for one matched topology. Never publish an unmatched Al result.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one
  branch/worktree/PR, incremental pushed checkpoints, review-ready teardown.
- B3LYP/def2-SVP/DF, GPU-first; OMP/MKL/OpenBLAS <=16; tee every long run.
  Isolate email extraction and Honcho inference under a bounded dead-man and
  restore them immediately on exit.
- Do not constrain Obr--H, weaken convergence, infer a saddle from a scan crest,
  reuse an incompatible optimizer/Hessian after stage changes, or rerun A1c's
  retired pre-equilibrium ensemble.
- Archive ignored campaign artifacts under `/mnt/data/vsletten/dissertation-data/`
  with independent hashes before teardown.

## Acceptance

- Regression tests prove the exact four-path Si scope, physical-H conservation
  (no unassigned or double-owned H), collision/finite-coordinate gates,
  persisted pre-climb band checkpoints, mechanism-versioned resume refusal,
  coupled-mode classification, and full IRC endpoint identity.
- Every Si path has one hash-valid terminal accepted/rejected/failed/blocked
  record and no `running` record. Only an accepted Si path can trigger its exact
  matched Al path.
- Either one matched Si/Al topology produces two fully gated barriers and an
  honest acid ordering, or all four Si paths are conclusive rejections and a
  narrower executable next-mechanism card is filed. Any optimizer/NEB/Hessian/
  IRC failure remains `incomplete-*`, never a scientific NO-GO.
- Full QM suite, whole-tree Ruff, changed-file format, CLI smoke, diff check,
  card/PLAN/STATUS bookkeeping, durable hashes, and review-ready teardown pass.

## Progress

- 2026-08-24 14:28 PDT — **POST-MERGE ACTIVATION CLOSEOUT:** independently
  verified squash merge `f5a3accf9668f6d3eb96fdff319e9530ebd97414`, the
  exact four terminal receipts, and all 24 archived evidence files (286,860
  bytes; every declared size and SHA-256 exact). The customer outcome remains
  4/4 typed endpoint rejections, zero failed/blocked/running paths, no admitted
  Al run, and READY narrower card A1f. Fresh merged-tree verification exposed
  one source gate missed by PR CI: Sourcery's review-fix regression additions
  left `test_concerted_acid_relay.py` inconsistent with whole-tree Ruff lint and
  formatting. This follow-up normalizes that test file without changing the
  scientific result. Fresh post-fix gates: **341 QM tests passed**, whole-tree
  Ruff check/format, CLI help, and diff check all green.
- 2026-08-24 13:08 PDT — **DONE / four-path model-valid rejection.** The exact
  mechanism-v1/gate-v1 B3LYP/def2-SVP/DF campaign finished with four terminal
  Si records: **4 rejected / 0 accepted / 0 failed / 0 blocked**, and no
  `running` record. Both three-water reactants were genuine zero-imaginary
  minima, but each exact hydrolyzed product relaxed to complete Si--Ow
  dissociation while Si--Obr remained cleaved. The four-water bridge-chain
  reactant moved one proton from H3O+ to the last shell water; the compact-cycle
  reactant moved one solvent proton onto a terminal framework O. These are
  conclusive typed endpoint rejections, not optimizer/Hessian/NEB/IRC failures.
  No Si path was eligible for CI-NEB or matched Al, so no downstream saddle,
  barrier, or ordering was fabricated. **DEVIATION:** the card listed every
  downstream artifact class, but its fail-closed endpoint gates correctly
  terminated all four paths before pre-climb; manufacturing NEB/saddle/IRC
  artifacts from rejected endpoints would violate the same contract. Narrower
  executable card `A1f-acid-neutral-water-attacker-relay` is filed and separates
  the H3O+ proton donor from the neutral-water nucleophile.
- 2026-08-24 12:41 PDT — implementation checkpoint: added mechanism-v1/gate-v1
  unconstrained concerted endpoint builders for the exact 3/4-water ×
  bridge-chain/compact-cycle scope, including stable A1c atom order and an
  explicit physical-H Grotthuss relay. The finite driver now fail-closes every
  bounded stage, persists CI-NEB pre-climb/climb checkpoints, requires a fresh
  directed saddle, a coupled relay/attack/cleavage one-mode Hessian, and Sella's
  full bidirectional Gonzalez–Schlegel IRC before barrier admission. Twenty-seven
  new regressions plus the cluster/TS suites are green (**132 passed**); targeted
  Ruff/format/diff and CLI-help gates pass. No production QM path has started.
- 2026-08-24 — created from A1c/A1d's model-valid finite-ensemble NO-GO. All 16
  unconstrained pre-equilibrated Si bridge seeds are conclusive rejections; the
  concerted hydronium/proton-relay path is now the only authorized acid-mechanism
  continuation. No A1e computation has started.

## Result

- Post-merge live acceptance reconfirmed the archived campaign rather than
  rerunning expensive DFT: evidence manifest SHA-256
  `3162e39f2f0f263444741b03e0a023b1209bf515a9f8291c330ca9995811e622`,
  run-manifest SHA-256
  `dfea131556fddf907ccabfc52945d390695f99d372cfdeac6859a4f2c0549982`,
  and production-log SHA-256
  `627914bbb22ee655932508ea974f006a503389198f656b173fe225607a8ca784`
  all match. The follow-up branch only closes the merged review regression's
  Ruff format/lint gate; terminal scientific evidence and verdict are unchanged.
- Added `concerted_acid_relay_endpoints()` and the finite
  `scripts/concerted_acid_relay.py` campaign. The builder preserves A1c atom
  order and exact physical-H ownership; bridge chains shift every relay proton
  owner explicitly. The driver binds topology/water/method/mechanism/gates and
  every stage bound into resume identity; refuses drift; persists/hash-journals
  terminal artifacts; and admits Al only after its exact Si path is fully
  accepted. `quarry.ts.full_irc()` adds bounded two-direction Sella
  Gonzalez--Schlegel integration using one shared TS Hessian.
- Production verdict:
  `all-four-si-paths-conclusive-rejection`. Three-water bridge-chain and compact
  cycle both rejected at optimized product signature `(False, False, True, 1,
  5, 0)` versus required `(True, False, True, 1, 6, 0)`. Four-water bridge-chain
  rejected exact solvent occupancies `(2,2,2,3)` versus `(3,2,2,2)`; four-water
  compact cycle rejected reactant signature `(False, True, True, 0, 8, 1)`
  versus `(False, True, True, 0, 9, 0)`. No computational stage failed.
- Durable external root:
  `/mnt/data/vsletten/dissertation-data/task200-a1e-acid-concerted-20260824/`.
  Evidence manifest: 24 files / 286,860 bytes, SHA-256
  `3162e39f2f0f263444741b03e0a023b1209bf515a9f8291c330ca9995811e622`.
  Campaign log SHA-256
  `627914bbb22ee655932508ea974f006a503389198f656b173fe225607a8ca784`;
  run manifest SHA-256
  `dfea131556fddf907ccabfc52945d390695f99d372cfdeac6859a4f2c0549982`.
  Every terminal and referenced artifact re-hashed exactly after completion.
- GPU isolation used a process-scoped 16 GiB memlock cap, cuTENSOR assertion,
  OMP/MKL/OpenBLAS 16, nice 10, stopped email pipelines/Honcho inference, and a
  ten-hour dead-man. On exit both email pipelines restarted on fresh PIDs,
  Honcho deriver/API are running, API `/health` is `{"status":"ok"}`, the
  dead-man is inactive, Ollama is empty, and the GPU returned idle.
- Fresh implementation verification before closeout: **338 full QM tests
  passed** with cache disabled; whole-tree Ruff check/format, CLI help, and diff
  check are green. Card/PLAN/STATUS and executable A1f follow-through are in the
  same PR.
