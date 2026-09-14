# A9b sensitivity campaign operator guide

A9b is a **survey-tier platform test**, not production kinetics. It simulates
exactly 29 one-at-a-time scenarios at the pH 4 reservoir anchor with the fixed
seeds `90401,90403,90407,90409,90413,90419,90421,90427` (232 runs). pH 3 and
pH 5 are laboratory-comparison bounds only; no pH 3/5 trajectory is run.

## Bundle contract

Each single-replica directory is atomically published and contains:

- `initial.pgif.json` and `final.pgif.json`;
- complete `events.jsonl` state transitions and physical timestamps;
- complete `likelihood.jsonl` physical/biased hazards and event factors;
- fixed-time `checkpoints.jsonl` with full state counts, geometric area, and
  lineage-valid Si/Al releases;
- hash-bound `metadata.json` and non-overwrite `receipt.json`.

The campaign keeps both the 320-center whole finite-deck census (260 reachable
plus 60 topology NO-GOs in its denominator/area) and the 260-center reachable
census. The former is the primary estimand. Analysis uses unnormalized,
log-domain importance sampling; it never self-normalizes. The merged trajectory
sampling contract requires weight ESS fraction >=0.1. Separately, positive
species estimates require contribution ESS >=2 as an analysis-precision/ranking
gate; that gate is not a sampling contract and does not alter trajectories.
Zero biased observations remain censored with no Poisson bound.

Estimator means, confidence intervals, and state-population evolution retain
their authoritative natural-log values in `log_mean`, `log_ci95`, and
`log_state_counts`. If a positive linear value cannot be represented as a
finite float, its linear field is `null` and the result is fail-closed with a
typed `censored_numerical_overflow` or `censored_numerical_underflow` status;
the ESS gates and unnormalized estimator are unchanged.

`analysis.json` is always a preliminary analyzer product: it sets
`artifact_verified=false`, labels numerical conclusions unverified, and leaves
acceptance/ranking pending. Only the independently recomputing verifier may
publish acceptance. Its atomically written `verification.json` is the commit
receipt and embeds the hash-bound `finalized_analysis`; no other artifact is an
accepted conclusion. Only rankable families enter the finalized preregistered
top-3. A family is called irrelevant only when every complete uncertainty
interval lies inside the documented ±0.1 log10 equivalence band. Laboratory
comparison reports explicit pH 3 and pH 5 bound gaps as well as the pH 4 gap;
only pH 4 is simulated. Palandri-Kharaka rates are mineral formula-unit rates,
not summed elemental-cation fluxes. The report therefore retains separate Si and
Al fluxes and converts them to kaolinite `Al2Si2O5(OH)4` formula units as
`min(Si/2, Al/2)` when both species estimates are positive and valid. If either
species is censored, the formula-unit rate and every laboratory gap are censored;
the Si+Al sum is never compared with the laboratory formula-unit rate.

## Short smoke

Build the runner with the repository-wide canonical Cargo target, then run:

```bash
export CARGO_TARGET_DIR=/home/vsletten/.hermes/cache/cargo-target/vsletten-dissertation-da68ebaefa3a
cargo build --release -p petra-deck --example a9b_importance_sampling \
  --manifest-path petra/Cargo.toml
python3 petra/scripts/a9b_sensitivity_ranking.py smoke /absolute/outside-git/a9b-smoke \
  --runner "$CARGO_TARGET_DIR/release/examples/a9b_importance_sampling" --workers 4
```

Smoke creates all 29 scenario decks and one short replica per scenario. Its
`scientific_acceptance` is always `not_evaluated_smoke`; it cannot satisfy the
232-run dissertation acceptance gate.

## Full bounded launch (parent/operator only)

Do not run the full campaign in an interactive worker. From the merged source,
invoke the public launcher directly:

```bash
/absolute/repo/petra/scripts/launch_a9b_sensitivity_campaign.sh \
  /absolute/outside-git/a9b-sensitivity-campaign 4
```

The public launcher re-execs itself exactly once inside
`a9b-sensitivity-campaign.service` via `systemd-run --user`, with the fixed
`RuntimeMaxSec=86400` ceiling. Entry to the workload requires a live D-Bus proof
that the current process cgroup belongs to that exact active user unit and that
its numeric `RuntimeMaxUSec` is positive and no greater than 86400 seconds.
Caller environment markers and invocation IDs cannot bypass this proof.
Inside that bound it rechecks 15-minute load, available memory, and swap-out
growth; sets OMP/OpenBLAS/MKL/Rayon/NumExpr aggregate caps <=16; re-execs at
niceness 10; uses at most four workers; writes a durable tee log; and atomically
publishes a v2 success/failure launch receipt under `<output>-operator/`. The
receipt records the live unit invocation ID and runtime ceiling in seconds and
microseconds. The campaign is
resumable only from hash-valid per-job receipts and refuses unreceipted output.

Initial preparation is restart-safe: the complete input/deck/manifest tree is
materialized in a temporary sibling and atomically published as the campaign
root. A preparation failure therefore leaves no nonempty manifest-less root.

## Lint-enabling scope note

`petra-core` now exposes fired-event site IDs as `usize`. The matching removal
of the obsolete `as usize` cast in `corrosion_study.rs` is an unrelated-behavior,
minimal workspace-lint compatibility fix: restoring the cast fails canonical
workspace Clippy with `clippy::unnecessary_cast` under `-D warnings`.
