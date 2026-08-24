# A1f-acid-neutral-water-attacker-relay — separate donor from nucleophile

- status: done
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

- 2026-08-24 16:09 PDT — **DONE / one-path conclusive rejection.** The exact
  mechanism-v2/gate-v2 B3LYP/def2-SVP/DF Si reactant converged in 68 steps with
  valid aggregate basin `(False, True, True, 0, 9, 0)`, unique ownership, and a
  neutral attacker, but physical H16 moved from the declared outer H3O+ donor
  onto the bridge-side relay water. Exact solvent occupancies became
  `(2,2,3,2)` instead of required `(3,2,2,2)`. The terminal is therefore
  `rejected` at `optimize-reactant`, not a computational failure; no CI-NEB,
  saddle, IRC, Al screen, barrier, or ordering was fabricated. Executable A1g
  now tests the selected bridge-side-H3O+ / distinct-neutral-attacker topology
  directly. All six archived files re-hash exactly; services and GPU isolation
  are restored healthy.
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

## Result

- Added `separated_acid_relay_endpoints()` and
  `scripts/separated_acid_relay.py`. The exact four-water topology preserves
  A1e atom order but assigns distinct outer H3O+ donor, neutral measured-flank
  attacker, bridge-side relay, and spectator waters. Its product moves donor H
  through the relay onto Obr while the attacker's own two H remain on the
  Si-bound nucleophile. Mechanism-v2/gate-v2 reuses A1e's bounded,
  hash-journaled endpoint/CI-NEB/fresh-saddle/full-IRC transaction; resumes
  bound to A1e v1 are refused.
- Nineteen new regressions pin exact scope, role-distinct atoms, every physical-H
  owner, finite/collision-safe endpoints, neutral-attacker product identity,
  six-component donor/relay/attack/cleavage mode coupling, exact IRC basins,
  conditional atom-matched Al, terminal identity, and incomplete failure
  semantics. Fresh final gates: **360 QM tests passed** with cache disabled;
  whole-tree Ruff check/format, CLI help, and `git diff --check` are green.
- Production verdict: `si-path-conclusive-rejection`. The optimized reactant
  retained intact Si--O--Si, neutral attacker, nine solvent H, zero Obr H, zero
  framework migration, and unique ownership, but changed exact solvent roles
  from `(3,2,2,2)` to `(2,2,3,2)`: physical H16 moved from outer donor O15 to
  bridge-side relay O22. Minimum pair distance is `0.964661 A`; no computational
  stage failed and no `running` record remains. No Al path was eligible.
- Durable evidence root:
  `/mnt/data/vsletten/dissertation-data/task201-a1f-acid-separated-20260824/`.
  Evidence manifest covers 6 files / 36,212 bytes and has SHA-256
  `940cd4b754f8b67244fd9e51812d4504a7872bb4cf1bb797bdc81776c7777d9c`;
  production log `bcaa68f2c415653b235c61c9b537ea0fe8f8cae6f35911c3f500a9d40243073b`,
  run manifest `3a72c6fdaef2ac8fe0730a91b44604e8b83dded4f9a5dfa7275818709dee9efb`,
  and terminal `0b51c68336cf3dff739f51135664f6ac7918bd2374500f642c8f112bc1b5f095`
  all re-hash exactly.
- GPU isolation used a process-scoped 16 GiB memlock cap, OMP/MKL/OpenBLAS 16,
  nice 10, stopped email pipelines/Honcho inference, empty Ollama, and a bounded
  ten-hour dead-man. The first wrapper mistakenly `exec`'d the campaign, which
  discards the shell's EXIT trap; restoration was immediately completed and
  verified manually. Both email pipelines are active on fresh PIDs, Honcho API
  is healthy at `/health`, the dead-man is inactive, Ollama is empty, and GPU
  utilization returned idle. The trap/`exec` pitfall is captured in learnings.
- Follow-through is executable card
  `A1g-acid-bridge-side-hydronium-neutral-attacker`; it promotes A1f's selected
  relay-H3O+ topology under a new identity instead of weakening this card's
  exact physical-H/family acceptance.
