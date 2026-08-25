# A1g-acid-bridge-side-hydronium-neutral-attacker — follow the migrated proton

- status: done
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaign)
- depends: A1f's exact donor-to-relay migration receipt
- claimed-by: hermes-workstation

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

- 2026-08-24 18:20 PDT — **DONE / one-path conclusive product rejection.**
  Mechanism-v3/gate-v3 started from A1f's exact hash-pinned Si reactant. It
  reconverged in 6 steps, retained exact `(2,2,3,2)` solvent roles and a neutral
  attacker, and had zero imaginary modes. The independent product converged in
  89 steps but did not retain the hydrolyzed family: Si--Ow remained nonbonded
  at 2.643 A, Si--Obr cleaved to 3.535 A, and attacker H21 became unassigned.
  Exact basin `(False,False,True,1,7,0)` therefore failed required
  `(True,False,True,1,8,0)`. This is a typed `rejected` terminal at
  `optimize-product`, not a computational failure; no CI-NEB, saddle, IRC, Al
  screen, barrier, or ordering was fabricated. Follow-up A1i permits exactly
  one calibrated production-tier revisit after A2, rather than another SVP
  topology chase. All archived files re-hash; GPU-isolated services are restored.
- 2026-08-24 18:08 PDT — implementation checkpoint: atomically claimed branch
  `agents/A1g-acid-bridge-side-hydronium-neutral-attacker`; added the exact
  four-water bridge-side-H3O+ / distinct-neutral-attacker endpoints and
  mechanism-v3/gate-v3 campaign. The Si path is hash-bound to A1f's archived
  optimized reactant `3a510764...`, preserves migrated H16 on bridge-side O22,
  transfers bridge-facing H23 directly to Obr, and refuses A1e/A1f terminals
  plus source-hash drift. New and inherited acid-path regressions are green
  (**70 passed**); targeted Ruff/format, source receipt, CLI help, and diff
  checks pass. No production DFT has started.
- 2026-08-24 — filed from A1f's exact unconstrained reactant: H16 migrated from
  the outer donor to the bridge-side relay oxygen, while the neutral attacker
  remained neutral. This card tests that selected topology directly instead of
  weakening A1f's physical-H/family gate.

## Result

- Added `bridge_side_hydronium_endpoints()` and
  `scripts/bridge_side_acid_relay.py`. The four-water atom map keeps migrated
  physical H16 on bridge-side O22, keeps measured-flank O19/H20/H21 as the
  distinct neutral attacker, and transfers bridge-facing H23 directly to Obr.
  Mechanism-v3/gate-v3 binds the Si terminal to A1f source SHA-256
  `3a5107648bd44cbba15dd7482b4a9e8b3c6434f049442453744737e2c32aba75`
  and refuses A1e/A1f terminal identities or source-hash drift.
- Production verdict: `si-path-conclusive-rejection`. The exact source reactant
  reconverged in 6 steps with zero imaginary modes, intact Si--O--Si, exact
  solvent occupancies `(2,2,3,2)`, unique H ownership, and neutral attacker.
  The independent product converged in 89 steps but reached basin
  `(False,False,True,1,7,0)`: Si--Ow `2.642649 A`, Si--Obr `3.535115 A`, and
  attacker H21 unassigned. Minimum pair distance stayed safe at `0.964884 A`.
  The failure is chemical family identity, not optimizer/Hessian/path failure.
- The terminal is hash-valid `rejected` at `optimize-product`; the manifest has
  one Si terminal and no `running` record. Since Si did not pass, the guarded
  campaign emitted no pre-climb band, saddle, IRC, Al path, barrier, or ordering.
- Durable evidence root:
  `/mnt/data/vsletten/dissertation-data/task202-a1g-acid-bridge-side-20260824/`.
  Evidence manifest covers 10 files / 62,705 bytes and has SHA-256
  `76d1d0275a773d7941a855250696cedc78a169ad68d58305195aa9df05d4564f`;
  production log `7c6ab526c88b86a9f7de4031b2beedbd3a2f4b8c3e8276d8a05927b0ae133a46`,
  run manifest `27f7a7040cb5f087d9aee7126676e25bed309eed243fd024a1579aa1f946fba2`,
  and terminal `cf74851972575f2a734d03b38c95fea85bc49acf67598739d3f01aea191c7834`
  all re-hash exactly.
- Twenty-one new A1g regressions plus inherited A1e/A1f coverage pin exact
  water scope, source identity, physical-H ownership, collision safety,
  direct-transfer mode classification, full-IRC endpoint identity, pre-climb
  checkpoint plumbing, conditional Al, and incomplete-failure semantics. Fresh
  final gates: **381 QM tests passed** with cache cleared; whole-tree Ruff check,
  all 45 Python files' format check, CLI help, and `git diff --check` are green.
  Email pipelines and Honcho inference were isolated under a ten-hour dead-man;
  both pipelines and both containers are restored healthy, the dead-man is
  inactive, Ollama is empty, and the GPU returned idle.
- The next acid decision is executable card
  `A1i-acid-production-tier-bridge-side-revisit`, blocked on A2's calibrated
  production method. No further B3LYP/def2-SVP water-count, donor-role, or
  topology retries are authorized by this result.
