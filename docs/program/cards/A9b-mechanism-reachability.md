# A9b-mechanism-reachability — close the lattice-desorption sampling gap

- status: ready
- track: A (geochemistry)
- priority: P0
- machine: any (analysis, deterministic fixtures, and CPU-only tests; no long campaign)
- depends: A9-approximate-rate-closure ✅
- claimed-by:

## Objective

Turn A9's `mechanism-unsampled` outcome into a mechanically decidable statement.
Starting from the exact 298 K `kaolinite-approx.toml` initial state, determine
whether each nonfrozen original-lattice Si and Al center has a legal reaction
path to its desorption-eligible state (`Si.oh4` / `Al.l6`) and then to the
corresponding release transition. Implement the smallest deterministic
reachability or state-seeding proof needed to exercise those paths while
preserving lattice/reservoir lineage.

This card selects and validates the rare-event/stiffness method needed before a
new sensitivity campaign. It does not rerun the 29 × 8 ensemble and does not
emit an absolute dissolution rate.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one branch, worktree,
  and PR.
- Work from the production deck and live Petra transition semantics, not from
  the prior result document's conclusions.
- Preserve original-lattice versus reservoir identity throughout every proof.
  Adsorbed reservoir cations may never count as lattice dissolution.
- Keep the method survey-tier and bounded. No QM, GPU, live service, database,
  credential, or >30-minute inline compute.
- Do not weaken A9's complete-event, PGIF, timestamp, or state-distribution
  validation gates.

## Acceptance

- A machine-readable graph/path report covers every initially nonfrozen Si and
  Al center and identifies either a concrete legal path to desorption or the
  exact missing transition/environment predicate.
- A deterministic integration fixture reaches each reachable Si/Al release
  transition and proves event-level lattice lineage, or the card records a
  verified topology-level NO-GO for any unreachable family.
- The selected finite rare-event/stiffness method is specified with its biasing
  or seeding semantics, unbiased estimator/reweighting rule if applicable,
  stopping rule, and failure/censor outcome. It must be executable on the
  existing deck without manufacturing an ordinal sensitivity rank.
- Focused tests, applicable Petra workspace tests/lint, and `git diff --check`
  pass. A9b-sensitivity-ranking is unblocked only if this card supplies an
  executable, independently verified sampling route.

### Verify

A different worker starts from the deck and Petra reaction graph, independently
recomputes the path census, adversarially checks lineage preservation and at
least one unreachable/tampered fixture, and confirms the proposed estimator
cannot turn non-observation into a positive rate or ordinal rank.

## Progress

- 2026-09-13 20:52 PDT (hermes-macbot-zero; profile=laptop) — Created from A9's review-ready closeout. The shipped A9 campaign observed zero original-lattice release and the corrected runner types a real trajectory `mechanism-unsampled`; this card owns the reachability and finite rare-event-method question rather than hiding it behind another A9 continuation.

## Result
