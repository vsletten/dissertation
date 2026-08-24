# A1c-acid-microsolvation-conformers — finite matched proton-relay ensemble

- status: active
- track: A (geochemistry)
- priority: P1
- machine: workstation (GPU campaign)
- depends: A1b microsolvation driver/receipt hardening merged
- claimed-by: hermes-workstation

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

- 2026-08-24 — campaign preflight correction: the first GPU-isolated launch was
  deliberately stopped during seed 1's HF pre-relaxation at step 76/100, before
  any production B3LYP result or accepted receipt. geomeTRIC convergence remains
  mandatory for every production endpoint, while the cheap HF stage now exposes
  and records its convergence bit but may hand its bounded last geometry to the
  production optimizer instead of burning the entire ensemble on a non-load-
  bearing precondition. Fresh 293-test/Ruff/format/diff gates pass; the stopped
  seed's running checkpoint will be quarantined automatically on relaunch.
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
