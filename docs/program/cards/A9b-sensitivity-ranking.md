# A9b-sensitivity-ranking — rerun corrected approximate-rate ensemble

- status: done
- track: A (geochemistry)
- priority: P0
- machine: workstation
- depends: A9b-mechanism-reachability; A9b-reservoir-origin-contract
- claimed-by: hermes-custom-build-001

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

- 2026-09-14 06:01 PDT (hermes-custom-build-001; profile=workstation) — DONE: bounded systemd campaign completed all 232 trajectories in 19.101 s at niceness 10 with 106.0 MiB peak memory and no swap. The independent verifier regenerated 29 decks, replayed all raw replicas, and byte-matched the finalized analysis. All scenarios completed and passed stationarity, but all Si/Al responses remain typed `censored_zero_biased_observations`; no family is rankable or demonstrated irrelevant.

- 2026-09-14 05:53 PDT (hermes-custom-build-001; profile=workstation) — CHECKPOINT: implemented and independently reviewed the survey-tier 29 × 8 pH-4 importance-sampling campaign, complete raw event/likelihood/PGIF/lineage evidence, log-domain censor-aware analyzer, independent byte-regenerating verifier, and self-enforcing bounded systemd launcher. Focused Python tests, the full Petra script suite, Rust workspace tests/Clippy, formatting, and diff checks pass. A real transient-unit admission probe proved exact cgroup/InvocationID/RuntimeMaxUSec membership and niceness 10; forged out-of-unit admission failed closed. The full campaign has not started yet; next step is the bounded autonomous launch and receipt-backed independent result verification.

- 2026-09-13 20:52 PDT (hermes-macbot-zero; profile=laptop) — Created from A9's review-ready closeout and blocked on the two upstream scientific contracts. A9's historical 29 × 8 bundle is preserved as negative evidence but predates the corrected mechanism-unsampled gate and cannot supply the requested ordinal ranking.

## Result

Implemented the reproducible corrected sensitivity harness, independent verifier,
bounded launcher, contract regressions, and the required engine support. The real
29-scenario × 8-seed campaign completed with all 232 raw evidence sets and passed
independent artifact replay (`verification.json` SHA-256
`66fdbde89522a00e612007d924ff3bed489bc15c9bd84267fd1f14d2a264b04d`).

Scientific verdict: all 29 scenarios completed and passed full per-kind
stationarity, but every Si and Al response is censored because no biased
original-lattice release was observed. Si/Al rates, stoichiometry, laboratory
gaps, ordinal ranks, and top-k are therefore unavailable; none of the seven
families is demonstrated irrelevant. Evidence and exact reproduction pointers:
[`docs/program/results/A9b-sensitivity-ranking.md`](../results/A9b-sensitivity-ranking.md).
