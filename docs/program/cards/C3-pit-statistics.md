# C3-pit-statistics — quantitative pitting statistics from the C2 deck

- status: done
- track: C (corrosion)
- priority: P2
- machine: any (CPU ensembles; macbot-friendly)
- depends: C2-corrosion-deck
- claimed-by: hermes-8d9e30fee91548

## Objective

Turn C2's qualitative shape study into quantitative pitting statistics on
the same `petra/decks/passive-metal.toml` deck family:

1. Scale replica counts well past C2's 7×32 (target ≥ 256 replicas per
   potential point) and extend the potential sweep resolution around the
   metastable→stable transition.
2. Fit induction-time distributions properly: Weibull parameter estimates
   with bootstrap CIs per potential, plus the survival/hazard diagnostic
   family Shibata-style analyses use. Report goodness-of-fit honestly —
   rate-compressed model, chemical-potential proxy for volts, so the claim
   is distribution-family agreement, not calibrated prediction (keep C2's
   honesty language).
3. Pit-size/depth distributions and spatial clustering statistics vs
   potential; nucleation-rate scaling; a finite-size-corrected estimate of
   the critical-potential crossover extending C2's three-size scan.

## Constraints

- CPU-only ensembles (B4 machinery); deterministic replays and exact-seed
  regression gates per C2. Laptop-class workers welcome; if a sweep rung
  exceeds host memory, record the ceiling (E2b pattern) rather than
  escalating silently.
- No new mechanism states in this card — deck physics changes are a
  separate card. Rate compression stays; every plot labels it.

## Acceptance

- Results doc with Weibull fits + CIs, induction/survival diagnostics,
  size/clustering distributions, nucleation-rate scaling, and the
  finite-size-corrected crossover, each with the shape-level honesty
  caveat stated.
- Sweep driver + regression tests landed; Petra tests, replay gates,
  Clippy/lint green; card/PLAN/STATUS bookkeeping in the PR.

## Progress

- 2026-09-01 — (hermes-8d9e30fee91548; profile=laptop) — completed the
  deterministic 12,032-replica C3 campaign: 11×256 primary potential runs and
  4×9×256 finite-size runs. Added right-censor-aware Weibull MLE + deterministic
  bootstrap CIs, Kaplan–Meier/Nelson–Aalen diagnostics, pit-size/dissolution-depth
  distributions, fixed-occupancy pair clustering, nucleation scaling, and an
  isotonic half-initiation crossover extrapolation. Independent review caught
  and the worker fixed the original Bernoulli rather than fixed-occupancy random
  pair baseline before final generation.
- 2026-08-29 — filed by Fable per Victor's fleet-routing directive: first
  machine:any C-track science card so laptop/cloud workers have corrosion
  work while the workstation runs QM campaigns.

## Result

Delivered `petra/crates/petra-deck/examples/corrosion_study.rs`, seven
deterministic CSV products under `docs/program/results/C3-pit-statistics/`, and
the interpretation at `docs/program/results/C3-pit-statistics.md`. Every
potential point has 256 replicas. Nucleation rate rises 157.0× across the
primary sweep; Weibull shape/scale carry censor-aware 95% bootstrap intervals;
pit size and the explicitly event-count-based depth proxy both grow with
potential. Four interpolated 50%-initiation points extrapolate linearly in
`1/L` to a deck-family point estimate μ₅₀(∞) ≈ -0.172 kcal/mol (R² 0.962),
explicitly not a calibrated voltage or universal critical point.

Verification: `cargo test -p petra-deck --example corrosion_study` passes 6/6;
`cargo test --workspace` passes 116/116; `cargo fmt --all --check`,
`cargo clippy --workspace --all-targets -- -D warnings`, and `git diff --check`
are clean. Two full release generations produced byte-identical SHA-256 hashes
for all seven CSVs. Independent pre-commit review found no security issue; its
blocking fixed-occupancy clustering finding was corrected and regression-tested
before these final gates.
