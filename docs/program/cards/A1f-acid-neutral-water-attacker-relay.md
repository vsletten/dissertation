# A1f-acid-neutral-water-attacker-relay — separate donor from nucleophile

- status: ready
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaign)
- depends: A1e's four-path model-valid rejection
- claimed-by: hermes-workstation

## Objective

Test the narrower acid mechanism exposed by A1e: **hydronium donates the
proton, but a distinct neutral water attacks silicon**. A1e used the hydronium
oxygen as both donor and eventual Si--Ow nucleophile. Both three-water product
endpoints relaxed to complete attacker dissociation, while both four-water
reactants redistributed the excess proton before path finding. Do not retry
those four topologies.

1. Build exactly one four-water Si--O--Si topology:
   `separated-donor-neutral-attacker`. Preserve the A1e dimer/solvent atom-order
   and physical-H contracts, but partition the solvent explicitly into H3O+
   donor, H2O attacker, one relay water, and one spectator/solvator. The
   unconstrained reactant has an intact unprotonated bridge, hydronium on the
   donor side, and a separate neutral attacker in the measured flank geometry.
2. Build and independently optimize the exact hydrolyzed product in which the
   neutral attacker remains Si-bound while the relay delivers one physical H to
   Obr. No Obr--H, Si--Ow, or Si--Obr production constraint may survive endpoint
   construction. Both endpoints require exact ownership, finite/collision,
   connectivity, zero-imaginary-mode, and family-identity gates.
3. Run one double-ended CI-NEB with every pre-climb image persisted, refine its
   crest with a fresh directed saddle search, require exactly one imaginary
   mode coupling donor/relay proton transfer, neutral-water attack, and
   Si--Obr cleavage, then run full bidirectional Gonzalez--Schlegel IRC into the
   exact typed endpoints.
4. Only an accepted Si path may trigger the atom-matched Si--O--Al path. Publish
   an ordering only if both saddles pass every gate. Any optimizer/NEB/Hessian/
   IRC failure is `incomplete-*`, never a rejection or barrier.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one
  branch/worktree/PR, incremental pushed checkpoints, review-ready teardown.
- Reuse A1e's mechanism-versioned journal, full-IRC implementation, hash-bound
  terminal records, and GPU isolation. Increment mechanism/gate identity; a
  resumed A1e checkpoint must be refused.
- B3LYP/def2-SVP/DF, GPU-first; OMP/MKL/OpenBLAS <=16; tee every long run.
  Isolate email extraction and Honcho inference under a bounded dead-man and
  restore immediately on exit.
- Scope is one Si topology plus its conditional exact Al match. No water-count
  sweep, no pre-equilibrated protonated bridge, no hydronium-as-nucleophile
  retry, and no unmatched Al result.
- Archive ignored campaign artifacts under `/mnt/data/vsletten/dissertation-data/`
  with independent hashes before teardown.

## Acceptance

- Regressions prove the donor/attacker atoms are distinct, exact four-water
  scope, physical-H conservation, no unassigned/double-owned H, finite and
  collision-safe endpoints, neutral-attacker product identity, persisted
  pre-climb checkpoints, resume/version refusal, coupled-mode classification,
  and exact full-IRC endpoint identity.
- The one Si path has one hash-valid terminal accepted/rejected/failed/blocked
  record and no `running` record. Al runs only after Si acceptance.
- Either the exact matched Si/Al pair yields two fully gated barriers and an
  honest ordering, or the one Si path is a conclusive rejection and the next
  mechanism decision is based on that receipt. Computational failure remains
  incomplete and is not laundered into science.
- Full QM suite, whole-tree Ruff, changed-file format, CLI smoke, diff check,
  card/PLAN/STATUS bookkeeping, durable hashes, and review-ready teardown pass.

## Progress

- 2026-08-24 16:00 PDT — implementation checkpoint: atomically claimed branch
  `agents/A1f-acid-neutral-water-attacker-relay`; added the exact four-water
  `separated-donor-neutral-attacker` R/P builder with A1e atom order, distinct
  H3O+ donor and neutral H2O attacker, one physical-H relay, spectator water,
  and unconstrained typed endpoints. Mechanism-v2/gate-v2 reuses A1e's bounded,
  hash-journaled CI-NEB/saddle/full-IRC driver while refusing v1 checkpoints and
  requiring all six donor/relay/attack/cleavage mode components. New and A1e/TS
  regressions are green (**86 passed**); targeted Ruff/format, CLI help, and
  diff check pass. No production DFT has started.
- 2026-08-24 — filed from A1e's exact four-path result. The narrower hypothesis
  is not “more waters”: it separates acid donor from neutral nucleophile so the
  attacking oxygen need not be the oxygen that began as H3O+.
