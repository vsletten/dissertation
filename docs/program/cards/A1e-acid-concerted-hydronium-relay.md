# A1e-acid-concerted-hydronium-relay — retire the failed pre-equilibrium

- status: ready
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaign)
- depends: A1c/A1d model-valid NO-GO for the pre-equilibrated protonated bridge
- claimed-by:

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
