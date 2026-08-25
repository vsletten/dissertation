# A1i-acid-production-tier-bridge-side-revisit — one calibrated acid retest

- status: blocked
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaign)
- depends: A1g-acid-bridge-side-hydronium-neutral-attacker, A2-production-energetics
- blocked-on: A2-production-energetics
- claimed-by:

## Objective

Run exactly one production-tier revisit of A1g's hash-pinned
`bridge-side-hydronium-neutral-attacker` Si path after A2 has settled and
verified the production geometry/solvation protocol. The scientific question is
narrow: does the calibrated tier stabilize both the exact A1g reactant and its
atom-matched hydrolyzed product, where B3LYP/def2-SVP/DF retained the reactant
minimum but left Si--Ow at 2.643 A and attacker H21 unassigned in the product?

1. Start from A1g's exact atom order and archived endpoints. Bind the new run to
   A1g terminal SHA-256
   `cf74851972575f2a734d03b38c95fea85bc49acf67598739d3f01aea191c7834`
   plus the A2 method/settings identity. Increment mechanism/gate identity and
   refuse all A1e/A1f/A1g checkpoints.
2. Independently optimize and frequency-check the unconstrained Si reactant and
   product at A2's settled production geometry/solvation tier. Preserve exact
   physical-H ownership, distinct neutral attacker, finite/collision,
   connectivity, zero-imaginary-mode, and family gates.
3. If either endpoint fails, emit one conclusive production-tier rejection and
   stop the acid topology program. Do not add another water count, conformer,
   donor role, restraint, or B3LYP/def2-SVP retry.
4. Only if both endpoints pass may the existing persisted pre-climb CI-NEB,
   fresh directed saddle, direct-transfer coupled-mode, full bidirectional IRC,
   and conditional atom-matched Al path run. Publish ordering only for two fully
   accepted saddles.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one
  branch/worktree/PR, incremental pushed checkpoints, and review-ready teardown.
- A2's merged method contract is binding. Do not invent a provisional
  production tier or reinterpret a single-point protocol as a geometry method.
- GPU-first; OMP/MKL/OpenBLAS <=16; tee every long run; isolate email/Honcho
  inference under a bounded dead-man and restore immediately on exit.
- Scope is one exact Si topology plus conditional exact Al. Computational
  failure remains `incomplete-*`, never rejection.
- Archive ignored artifacts under `/mnt/data/vsletten/dissertation-data/` with
  independent hashes.

## Acceptance

- A2 is done and the run records its exact merged settings/method provenance.
- Regressions prove A1g source/terminal binding, prior-version refusal, endpoint
  ownership/connectivity/minimum gates, direct-transfer mode classification,
  persisted pre-climb checkpoints, exact full-IRC identity, conditional Al, and
  incomplete-failure semantics.
- The one Si path has one hash-valid terminal and no `running` record. Either it
  conclusively rejects at the calibrated endpoint gate, or a fully gated matched
  Si/Al pair yields two barriers and an honest ordering.
- Full QM suite, whole-tree Ruff, changed-file format, CLI smoke, diff check,
  card/PLAN/STATUS bookkeeping, durable hashes, and review-ready teardown pass.

## Progress

- 2026-08-24 — filed from A1g terminal `cf748519...`. B3LYP/def2-SVP/DF
  stabilizes the exact bridge-side-H3O+ reactant but not its independent
  hydrolyzed product. This card is the only authorized acid-topology revisit and
  remains blocked until A2 establishes the calibrated production tier.
