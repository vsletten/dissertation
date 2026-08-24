# A1g-acid-bridge-side-hydronium-neutral-attacker — follow the migrated proton

- status: ready
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaign)
- depends: A1f's exact donor-to-relay migration receipt
- claimed-by:

## Objective

Test the single acid topology selected by A1f's unconstrained reactant rather
than fighting its proton transfer: **bridge-side H3O+ donates directly to Obr
while a distinct neutral water attacks silicon**. A1f began with H3O+ on the
outer donor oxygen and a neutral relay water next to the bridge. Its exact
B3LYP/def2-SVP/DF reactant converged cleanly only after physical H16 moved from
the outer donor to the relay oxygen, yielding occupancies `(2,2,3,2)` for
outer donor / attacker / bridge-side relay / spectator. The neutral attacker
remained H2O. Do not relabel that A1f endpoint as accepted; it selects this new,
mechanism-versioned topology.

1. Build exactly one four-water Si--O--Si topology:
   `bridge-side-hydronium-neutral-attacker`. Preserve A1e/A1f atom order and
   physical-H identity, but seed H3O+ at the bridge-side solvent oxygen, keep the
   measured-flank attacker as a distinct neutral H2O, and keep the other two
   waters neutral spectators/solvators. The unconstrained reactant requires
   exact occupancies `(2,2,3,2)`, an intact unprotonated bridge, and a neutral
   attacker.
2. Build and independently optimize the atom-matched hydrolyzed product: one
   physical H moves from bridge-side H3O+ to Obr while the distinct attacker
   remains Si-bound. No Obr--H, Si--Ow, or Si--Obr production constraint may
   survive endpoint construction. Both endpoints require exact ownership,
   finite/collision, connectivity, zero-imaginary-mode, and family-identity
   gates.
3. Run one double-ended CI-NEB with every pre-climb image persisted; refine a
   fresh directed saddle; require exactly one mode coupling direct hydronium
   proton transfer, neutral-water attack, and Si--Obr cleavage; then run full
   bidirectional Gonzalez--Schlegel IRC into the exact typed endpoints.
4. Only accepted Si may trigger the atom-matched Si--O--Al path. Publish an
   ordering only if both saddles pass every gate. Computational failure remains
   `incomplete-*`, never rejection.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one
  branch/worktree/PR, incremental pushed checkpoints, review-ready teardown.
- Reuse A1f's mechanism-versioned journal, full-IRC implementation, hash-bound
  terminals, and GPU isolation. Increment mechanism/gate identity again; refuse
  A1e/A1f checkpoints.
- B3LYP/def2-SVP/DF, GPU-first; OMP/MKL/OpenBLAS <=16; tee every long run.
  Isolate email extraction and Honcho inference under a bounded dead-man and
  restore immediately on exit. A trap-owning shell must not `exec` the campaign.
- Scope is one Si topology plus its conditional exact Al match. No water-count
  sweep, no A1e hydronium-as-nucleophile retry, no A1f outer-donor retry, and no
  unmatched Al result.
- Archive ignored artifacts under `/mnt/data/vsletten/dissertation-data/` with
  independent hashes before teardown.

## Acceptance

- Regressions prove distinct donor/attacker atoms, exact four-water scope,
  physical-H conservation, no unassigned/double-owned H, finite/collision-safe
  endpoints, neutral-attacker product identity, persisted pre-climb
  checkpoints, A1e/A1f resume refusal, direct-transfer coupled-mode
  classification, and exact full-IRC endpoint identity.
- The one Si path has one hash-valid terminal accepted/rejected/failed/blocked
  record and no `running` record. Al runs only after Si acceptance.
- Either the matched Si/Al pair yields two fully gated barriers and an honest
  ordering, or the one Si path is a conclusive rejection and the next mechanism
  decision cites that receipt. Computational failure stays incomplete.
- Full QM suite, whole-tree Ruff, changed-file format, CLI smoke, diff check,
  card/PLAN/STATUS bookkeeping, durable hashes, and review-ready teardown pass.

## Progress

- 2026-08-24 — filed from A1f's exact unconstrained reactant: H16 migrated from
  the outer donor to the bridge-side relay oxygen, while the neutral attacker
  remained neutral. This card tests that selected topology directly instead of
  weakening A1f's physical-H/family gate.
