# A9b-sensitivity-ranking — rerun corrected approximate-rate ensemble

- status: blocked
- track: A (geochemistry)
- priority: P0
- machine: workstation
- depends: A9b-mechanism-reachability; A9b-reservoir-origin-contract
- claimed-by:

## Objective

After the reachability/sampling method and realistic pH 3–5 origin contract are
merged, rerun the A9 nominal plus one-at-a-time sensitivity design at 298 K:
±3 kcal/mol and ×/÷10 prefactors for all seven barrier families, at least eight
fixed-seed replicas per scenario. Produce the corrected absolute Si/Al rates,
stoichiometry, state-population evolution, laboratory comparison, and a
censor-aware family ranking.

This card owns the final ranking that A9 could not produce. It may return an
honest, independently verified all-censored or topology-level NO-GO if the new
method still cannot estimate a response; it must never force ordinal ranks from
zero or undefined values.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one branch, worktree,
  and PR.
- Launch any >30-minute campaign as one bounded workstation unit with an atomic
  `file:` receipt, thread caps ≤16, niceness ≥10, and durable tee logs. Obey the
  load/memory gate before launch.
- Use the merged reachability/sampling and reservoir/origin contracts unchanged.
  No new gate infrastructure, QM, GPU, live database, service activation, or
  credential work belongs here.
- Preserve complete event streams, PGIF/timestamp cross-checks, full per-kind
  state-distribution stationarity, lattice/reservoir lineage, and typed censoring.

## Acceptance

- A reproducible 29-scenario × ≥8-replica bundle passes the merged verifier and
  reports estimator, stopping, stationarity, lineage, and censor outcomes for
  every trajectory and family.
- The results document reports original-lattice Si and Al rates in mol m^-2 s^-1,
  ensemble bands or valid bounds, Si:Al stoichiometry, pH 3–5 laboratory gaps,
  and a censor-aware sensitivity table.
- The verdict explicitly names families the data can rank, families it cannot,
  the top-k survey targets if any, and families shown irrelevant at this tier.
  Undefined responses remain unranked.
- Focused tests, Petra workspace tests/lint, artifact verification, and
  `git diff --check` pass; one PR carries results, card, PLAN, and STATUS.

### Verify

A different worker starts from the raw campaign and merged contracts, rehashes
and independently recomputes conversions, rates/bounds, stationarity, lineage,
sensitivity deltas, censoring, and the top-k/irrelevant verdict. Until this pass,
all numerical conclusions are labeled `unverified`.

## Progress

- 2026-09-13 20:52 PDT (hermes-macbot-zero; profile=laptop) — Created from A9's review-ready closeout and blocked on the two upstream scientific contracts. A9's historical 29 × 8 bundle is preserved as negative evidence but predates the corrected mechanism-unsampled gate and cannot supply the requested ordinal ranking.

## Result
