# C3 pit statistics — quantitative ensemble study of the C2 deck

## Verdict

**PASS for the card's quantitative, shape-level statistics gate.** The unchanged
`petra/decks/passive-metal.toml` mechanism now has a deterministic 12,032-replica
study: 11 driving-potential points × 256 replicas on the 24² deck, plus a
4-size × 9-potential × 256-replica finite-size sweep. Induction times are fit by
a right-censor-aware Weibull maximum-likelihood estimator with deterministic
2,000-resample bootstrap confidence intervals. The exported diagnostics include
Kaplan–Meier survival, Nelson–Aalen cumulative hazard, per-pit final size and a
coarse-grained dissolution-event depth proxy, nucleation-rate scaling, spatial
nearest-neighbor clustering, and a finite-size extrapolation of the 50%
initiation crossover.

This remains a **rate-compressed model-family study, not a calibrated alloy
prediction**. The driving coordinate is `mu(e_drive)` in kcal/mol, not electrode
voltage. No experimental curve was digitized or fit. The defensible claim is
that this declared CTMC mechanism produces a Weibull-compatible induction-time
diagnostic family and systematic potential/size trends; it is not numerical
agreement with Shibata, Valor, Punckt, a particular stainless steel, or a
critical voltage.

## Reproduction and scope

Run from `petra/` using the repository's mechanically selected canonical Cargo
target:

```bash
CARGO_TARGET_DIR=$(python3 /path/to/mission-control/scripts/review_ready_teardown.py \
  --repo vsletten/dissertation --print-cargo-target)
export CARGO_TARGET_DIR
OMP_NUM_THREADS=4 RAYON_NUM_THREADS=4 nice -n 10 \
  cargo run --release -p petra-deck --example corrosion_study -- \
  decks/passive-metal.toml ../docs/program/results/C3-pit-statistics
```

Every trajectory advances to the exact 600 s CTMC deadline and uses Petra's
hash seed policy. A 120,000-event fail-closed ceiling bounds each replica. The
primary sweep uses μ = -0.20, -0.05, +0.05, +0.15, +0.25, +0.35, +0.45,
+0.55, +0.65, +0.75, and +0.85 kcal/mol. The finite-size sweep uses 16²,
24², 32², and 40² grids at μ = -0.30 through +0.50 in 0.10 increments.
All potential points use 256 replicas.

A stable initiation remains C2's first nearest-neighbor `bare|dissolving`
cluster of at least six patches. Right-censored replicas contribute their full
600 s exposure to the likelihood and nucleation-rate denominator.

## Induction-time distributions and nucleation scaling

The Weibull likelihood includes both events and right-censored observations.
For shape `k`, scale `λ` is profiled analytically and `k` is solved from the
profile score. Confidence bands are deterministic nonparametric percentile
bootstrap intervals over whole replica observations. `weibull_ks_distance` is
a machine-visible, KS-style maximum distance between the fitted survival curve
and post-event Kaplan–Meier survival; it is a fit diagnostic, not a formal
model-selection p-value. The study does not claim Weibull is uniquely preferred
over exponential, lognormal, or nonhomogeneous-Poisson alternatives.

| μ (kcal/mol) | Initiated / 256 | Weibull k (95% CI) | Weibull λ, s (95% CI) | KS-style distance | Nucleation rate, patch⁻¹ s⁻¹ |
|---:|---:|---:|---:|---:|---:|
| -0.20 | 12 | 3.330 (2.343–5.647) | 1492.8 (1021.7–2421.3) | 0.010 | 1.372×10⁻⁷ |
| -0.05 | 33 | 2.631 (1.902–3.818) | 1274.0 (997.2–1833.4) | 0.014 | 3.871×10⁻⁷ |
| +0.05 | 75 | 2.336 (1.952–2.862) | 942.5 (835.0–1097.5) | 0.018 | 9.380×10⁻⁷ |
| +0.15 | 111 | 1.864 (1.658–2.124) | 799.6 (719.6–909.9) | 0.055 | 1.535×10⁻⁶ |
| +0.25 | 169 | 1.897 (1.744–2.083) | 553.9 (512.4–607.5) | 0.090 | 2.755×10⁻⁶ |
| +0.35 | 239 | 1.979 (1.835–2.150) | 344.8 (321.5–367.2) | 0.069 | 5.438×10⁻⁶ |
| +0.45 | 252 | 2.155 (1.966–2.406) | 252.3 (237.5–269.0) | 0.080 | 7.700×10⁻⁶ |
| +0.55 | 255 | 2.079 (1.885–2.340) | 197.5 (184.5–210.4) | 0.111 | 9.938×10⁻⁶ |
| +0.65 | 256 | 2.760 (2.450–3.265) | 137.8 (131.5–144.4) | 0.116 | 1.410×10⁻⁵ |
| +0.75 | 256 | 3.763 (3.488–4.119) | 103.6 (100.0–107.0) | 0.060 | 1.853×10⁻⁵ |
| +0.85 | 256 | 3.934 (3.601–4.380) | 88.7 (85.8–91.6) | 0.080 | 2.154×10⁻⁵ |

Across the sampled range, the exposure-normalized nucleation rate rises
157.0-fold. Scale `λ` falls correspondingly. The fitted shape is non-monotonic:
it is near 1.9–2.2 through the central transition and rises above 3.7 only after
all replicas initiate. That is reported as model behavior, not forced into one
potential-independent Weibull family.

`survival-hazard.csv` carries every unique event/censor time with at-risk,
event, and censor counts; Kaplan–Meier survival; Nelson–Aalen cumulative
hazard; and the matched fitted Weibull survival/hazard. This preserves the
Shibata-style diagnostic family without discarding censored runs.

## Pit size, depth proxy, and clustering

`pit-size-depth.csv` contains one row per final connected pit cluster. Size is
the number of final `bare|dissolving` patches. “Depth” is explicitly the count
of `metal_dissolution` events accumulated at sites that belong to that final
cluster; it is a coarse-grained event-depth proxy, not nanometres and not a
surface reconstruction. Dissolution at a site that later repassivates outside a
final pit is not assigned retrospectively to that pit.

Selected distribution summaries are:

| μ | Final pits | Size p50 / p95 / max (patches) | Max-depth p50 / p95 / max (events) |
|---:|---:|---:|---:|
| -0.20 | 2,564 | 1 / 2 / 5 | 0 / 2 / 6 |
| +0.05 | 3,738 | 1 / 2 / 7 | 1 / 4 / 10 |
| +0.25 | 5,040 | 1 / 3 / 7 | 1 / 6 / 15 |
| +0.45 | 6,344 | 1 / 3 / 11 | 2 / 9 / 24 |
| +0.65 | 7,987 | 1 / 4 / 15 | 5 / 15 / 37 |
| +0.85 | 9,599 | 1 / 5 / 18 | 9 / 24 / 56 |

Mean pit density rises from 0.01739 to 0.06510 patches⁻¹ (3.74×), while mean
largest-cluster size rises from 2.00 to 7.05 patches (3.53×) and its p95 rises
from 4 to 12 patches. The fixed-occupancy nearest-neighbor pair ratio stays above random mixing
but declines from 3.205 (defined in 255/256 replicas) to 1.986 (256/256) as
occupancy increases. Replicas with fewer than two active sites have no defined
pair baseline and are excluded rather than entered as zero. Thus absolute clusters
and depth grow with driving potential while relative pair enrichment weakens;
the study does not hide the distinction or claim percolation.

The fixed-occupancy baseline conditions on occupancy but not on location. In
particular, it does not distinguish interaction-driven clustering from spatial
enrichment near the deck's single declared defect; the ratio is a deck-family
spatial diagnostic, not an interaction-only estimator.

## Finite-size-corrected crossover

A susceptibility maximum proved grid-boundary-sensitive during development, so
the final driver fail-closes on a better-defined crossover: the potential at
which 50% of replicas have initiated by 600 s. For each size, equal-weight
pool-adjacent-violators regression removes only finite-replica reversals in the
initiation fractions; linear interpolation then locates μ₅₀. The sampled
brackets and crossings are:

| Grid | Lower bracket | Upper bracket | Interpolated μ₅₀ (kcal/mol) |
|---:|---:|---:|---:|
| 16² | +0.20: 88/256 | +0.30: 139/256 | +0.278431373 |
| 24² | +0.10: 92/256 | +0.20: 141/256 | +0.173469388 |
| 32² | 0.00: 96/256 | +0.10: 151/256 | +0.058181818 |
| 40² | -0.10: 54/256 | 0.00: 130/256 | -0.002631579 |

A linear `μ₅₀(L) = μ₅₀(∞) + a/L` point-estimate fit gives:

- **μ₅₀(∞) ≈ -0.172 kcal/mol**
- **a ≈ 7.44 kcal/mol·patches**
- **R² = 0.962**

This is a finite-size correction for this exact deck family, fixed 600 s
horizon, six-patch initiation definition, and one central declared defect per
grid. Because the defect density changes with grid area, it is not a universal
thermodynamic critical point and must not be translated to volts.
The CSV retains full deterministic point-estimate precision for replay, but
binomial uncertainty in each 256-replica crossing is not propagated through
the four-point `1/L` fit; the last digits are not an uncertainty claim.

## Artifacts

- Driver and statistical regressions:
  `petra/crates/petra-deck/examples/corrosion_study.rs`
- C2 physics/replay regressions:
  `petra/crates/petra-deck/tests/corrosion.rs`
- Primary aggregate:
  `docs/program/results/C3-pit-statistics/potential-summary.csv`
- Replica observations:
  `docs/program/results/C3-pit-statistics/induction-times.csv`
- Survival/hazard diagnostics:
  `docs/program/results/C3-pit-statistics/survival-hazard.csv`
- Pit size/depth distributions:
  `docs/program/results/C3-pit-statistics/pit-size-depth.csv`
- Finite-size grid and extrapolation:
  `docs/program/results/C3-pit-statistics/critical-scaling.csv` and
  `crossover-fit.csv`
- Representative metastable current trace:
  `docs/program/results/C3-pit-statistics/metastable-transient.csv`

The generated CSVs are deterministic products of the checked-in driver, deck,
and seeds. Regression coverage includes censor-aware Weibull recovery,
deterministic bootstrap intervals, survival/censor tie handling, isotonic
half-initiation interpolation, inverse-size extrapolation, the original
same-seed bitwise trajectory gate, and the C2 low/high corrosion transition.

Replica indices intentionally reuse the same derived seeds across potentials
and lattice sizes, providing common random numbers for condition contrasts.
Consequently, conditions are paired rather than independent. The 24×24,
μ=-0.20 scaling rung also repeats the corresponding 256 primary trajectories;
it is not an additional independent sample.
