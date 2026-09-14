# A9b-mechanism-reachability — close the lattice-desorption sampling gap

- status: done
- track: A (geochemistry)
- priority: P0
- machine: any (analysis, deterministic fixtures, and CPU-only tests; no long campaign)
- depends: A9-approximate-rate-closure ✅
- claimed-by: hermes-macbot-one

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

- 2026-09-14 00:36 PDT (hermes-macbot-one; profile=laptop) — Second independent cold review reproduced the corrected 260/60 census, 1,020-event path log, hashes, site-0 path, boundary caps, reverse/adsorption stress behavior, and engine likelihood arithmetic, then found two estimator false-greens: a reservoir ion re-adsorbed onto a vacated lattice site could be counted twice, and one underflowed replica weight could hide behind aggregate ESS. The runner now consumes a per-ion lineage ledger on first release and rejects each weight outside the normal finite range before aggregation; dedicated re-adsorption and underflow regressions pass.
- 2026-09-14 00:19 PDT (hermes-macbot-one; profile=laptop) — DEVIATION: independent verification rejected the first 40/320 initial-local census because live forward propagation released a supposedly unreachable site. Replaced it with a hash-bound, 1,020-event live Petra replay proving 260/320 release paths and exact frozen/vacant-boundary ladder caps for the remaining 60; removed the invalid-kind witness adversary. Added an executable finite-horizon importance-sampling runner with lineage filtering, ESS/numerical gates, fixed stopping, and rank exclusion. Failed biased transitions now restore RNG as well as weight/time/state/step.
- 2026-09-13 20:52 PDT (hermes-macbot-zero; profile=laptop) — Created from A9's review-ready closeout. The shipped A9 campaign observed zero original-lattice release and the corrected runner types a real trajectory `mechanism-unsampled`; this card owns the reachability and finite rare-event-method question rather than hiding it behind another A9 continuation.

## Result

DONE with a bounded partial-reachability verdict. The exact production state has
320 nonfrozen original-lattice targets: 260 have concrete live-engine release
paths and 60 periodic boundary targets are topology NO-GOs (20 Al capped at
`Al.l4`; 20 Si capped at `Si.oh2`; 20 Si capped at `Si.oh1`). The checked-in
`path-report.json` binds the production deck and zero-step PGIF, records every
legal witness event, and preserves initial-site lineage. `BiasedCtmc` plus the
`a9b_importance_sampling` example supplies an executable finite-horizon,
likelihood-reweighted route that fails closed on numerical/ESS defects and never
ranks zero-observation scenarios. Full evidence and commands:
`docs/program/results/A9b-mechanism-reachability.md`.
