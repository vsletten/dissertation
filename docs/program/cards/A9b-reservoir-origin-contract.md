# A9b-reservoir-origin-contract — model pH 3–5 without losing cation identity

- status: done
- track: A (geochemistry)
- priority: P0
- machine: any (design, implementation, and CPU-only tests; no campaign)
- depends: A9-approximate-rate-closure ✅
- claimed-by: hermes-macbot-zero

## Objective

Replace A9's numerical `activity(Al) = activity(Si) = 1e-30` open-flow sink with
a defensible far-from-equilibrium pH 3–5 reservoir contract while keeping
original-lattice and solution/reservoir cations distinguishable through
adsorption, surface reaction, and desorption. Define the activities/species and
boundary semantics, implement them in the deck/runner interface, and prove that
only release of an original lattice cation contributes to physical dissolution.

This card establishes the solution and identity contract. It does not rerun the
29 × 8 sensitivity ensemble and does not publish a calibrated rate.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one branch, worktree,
  and PR.
- Ground pH/speciation/activity choices in the existing kinetics ledger and
  cited field/lab sources. State assumptions and excluded chemistry plainly.
- Never infer lattice origin from species alone after adsorption. If Petra lacks
  durable particle/site provenance, add the minimum explicit lineage state or
  use a rigorously demonstrated no-reentry boundary; do not relabel reservoir
  cations as lattice material.
- No QM, GPU, live service, database, credential, or >30-minute inline compute.
- Preserve A9's temperature, unit conversions, and event completeness gates.

## Acceptance

- A provenance table defines pH 3, 4, and 5 reservoir species/activities,
  references, temperature, open-flow boundary behavior, and every approximation.
- The deck/runner rejects the old `1e-30` sink as a realistic pH condition and
  emits the selected reservoir contract explicitly in generated evidence.
- Deterministic adsorption/desorption fixtures prove reservoir re-entry cannot
  increment original-lattice Si/Al release, while a true original-lattice
  release does increment the correct numerator and stoichiometry.
- Focused tests, applicable Petra workspace tests/lint, and `git diff --check`
  pass. A9b-sensitivity-ranking is unblocked only when the contract is
  executable and origin-safe.

### Verify

A different worker begins from the cited pH 3–5 laboratory/ledger sources and
live transition semantics, independently checks the activity arithmetic and
origin accounting, and adversarially exercises adsorb→surface-state→desorb
re-entry without accepting it as lattice dissolution.

## Progress

- 2026-09-13 23:32 PDT (hermes-macbot-zero; profile=laptop) — DONE / REVIEW-READY: added the machine-readable pH 3/4/5 constant-activity open-flow contract, a fail-closed materializer, three executable hash-bound decks/evidence receipts, and deterministic Si+Al production-replay origin fixtures. The preserved A9 `1e-30` deck remains historical only; generated profiles use explicit `a(H+) = 10^-pH`, declared `1e-6` trace-product activities, and the P&K `n_H = 0.777` acid-order surrogate on terminal release. Full Petra workspace tests passed; focused Python tests, Ruff, Rust format, three real Petra deck-load smokes, citation verification, and `git diff --check` passed. DEVIATION: the activity contract is deliberately a trace-inflow sensitivity boundary, not a saturation/speciation solver or calibrated rate model.

- 2026-09-13 20:52 PDT (hermes-macbot-zero; profile=laptop) — Created from A9's review-ready closeout. The shipped campaign's `1e-30` product sink made historical trajectories origin-safe only because no adsorption occurred; it was not a realistic pH 3–5 solution model. This card owns that contract gap.

## Result

2026-09-13 23:32 PDT (hermes-macbot-zero; profile=laptop) — The executable contract is `petra/examples/kaolinite-reservoirs.toml`; `petra/scripts/reservoir_origin_contract.py` validates/materializes a selected pH 3–5 profile and emits its complete source/approximation/origin contract in hash-bound JSON. Generated decks and receipts live under `docs/program/results/a9b-reservoir-origin-contract/`; the narrative provenance and limitations are in `docs/program/results/A9b-reservoir-origin-contract.md`. Deterministic fixtures prove true original-lattice Si/Al desorption increments the physical numerator while reservoir adsorb → surface-state → desorb cycles contribute only gross exchange. This card publishes no calibrated rate and runs no sensitivity ensemble; the ranking card remains blocked only on the separate mechanism-reachability contract and its independent Verify edge.
