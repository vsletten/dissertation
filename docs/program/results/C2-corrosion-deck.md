# C2 corrosion deck — passive-film breakdown ensemble

## Verdict

**PASS for the card's qualitative Phase-2 statistics gate.** Petra now ships a
coarse-grained CTMC passive-film deck whose local patches progress through
`intact → adsorbed → thinned → bare ↔ dissolving` while competing
repassivation returns active sites through `repassivated → intact`. A bounded
potential sweep produces the requested pit-number rise, spatial clustering,
induction-time distribution diagnostics, and dissolution/repassivation current
pulses. A separate three-size scan brackets the finite-size crossover regime.

This is a **shape/conformance study**, not a calibrated alloy prediction. The
rates are compressed for bounded execution and the electrochemical driving
term is expressed as a chemical-potential proxy in kcal/mol, not volts. The
comparison with Shibata is limited to reporting the same induction-time /
Weibull-plot diagnostic family; no experimental goodness-of-fit or numerical
agreement is claimed.

## Artifacts

- Deck: `petra/decks/passive-metal.toml`
- Deterministic study: `petra/crates/petra-deck/examples/corrosion_study.rs`
- Acceptance regressions: `petra/crates/petra-deck/tests/corrosion.rs`
- Potential summary: `docs/program/results/C2-corrosion-deck/potential-summary.csv`
- Replica induction times: `docs/program/results/C2-corrosion-deck/induction-times.csv`
- Finite-size scan: `docs/program/results/C2-corrosion-deck/critical-scaling.csv`
- Representative current series:
  `docs/program/results/C2-corrosion-deck/metastable-transient.csv`

Reproduce from `petra/` using the repository's mechanically selected canonical
Cargo target:

```bash
CARGO_TARGET_DIR=$(python3 /path/to/mission-control/scripts/review_ready_teardown.py \
  --repo vsletten/dissertation --print-cargo-target)
export CARGO_TARGET_DIR
cargo test -p petra-deck --test corrosion -- --nocapture
cargo run -p petra-deck --example corrosion_study -- \
  decks/passive-metal.toml \
  ../docs/program/results/C2-corrosion-deck
```

The primary study is bounded at 7 driving conditions × 32 deterministic
replicas. The finite-size scan adds 16 replicas at each combination of 16²,
24², and 32² patches with μ = +0.30, +0.55, and +0.80 kcal/mol. Every trajectory
uses an exact 600 s CTMC deadline via `advance_ctmc_until`; a fail-closed
40,000-event cap prevents runaway execution. No reaction at or beyond the
deadline is applied or counted.

A stable initiation is the first observed nearest-neighbor
`bare|dissolving` cluster of at least six patches. The Weibull columns are the
median-rank slope and R² of **uncensored initiations only**. They are descriptive
plots, not censor-aware parameter estimates or distribution-selection tests.
The conditional mean column is named `mean_induction_uncensored_s` for the same
reason; censored replicas remain explicit in `induction-times.csv`.

## Potential-driven statistics

| μ(e_drive) | Initiated / 32 | Descriptive Weibull slope / R² | Conditional mean induction (s) | Mean pit count | Pit density | Mean largest cluster | NN pair ratio | Susceptibility |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -0.70 | 0 | — | — | 3.47 | 0.0060 | 1.25 | 4.55 | 0.0088 |
| -0.45 | 1 | — | 569.45 | 6.91 | 0.0120 | 1.50 | 3.13 | 0.0172 |
| -0.20 | 2 | — | 423.64 | 9.66 | 0.0168 | 1.78 | 2.19 | 0.0162 |
| +0.05 | 16 | 2.82 / 0.959 | 395.75 | 14.72 | 0.0256 | 2.88 | 3.09 | 0.0412 |
| +0.30 | 26 | 2.43 / 0.969 | 269.04 | 20.28 | 0.0352 | 3.50 | 2.56 | 0.0393 |
| +0.55 | 32 | 2.42 / 0.845 | 173.08 | 29.44 | 0.0511 | 5.28 | 2.24 | 0.0923 |
| +0.80 | 32 | 4.23 / 0.939 | 81.58 | 36.22 | 0.0629 | 7.28 | 2.03 | 0.0699 |

The low-driving regime produces sparse active patches and no stable six-patch
pit within 600 s. The midpoint has 50% initiation at +0.05 kcal/mol, and all
replicas initiate by +0.55. Mean pit count increases 10.4-fold from 3.47 to
36.22, while mean largest cluster increases 5.8-fold from 1.25 to 7.28. The
nearest-neighbor pair ratio remains above the random-mixing value of one at
every condition; its decrease at high occupancy is reported rather than hidden,
while the absolute cluster-size and pit-density measures increase strongly.

The positive Weibull-plot slopes and reported R² values describe the uncensored
samples in the same coordinates commonly used for Shibata-style induction-time
analysis. They do **not** establish that Weibull is uniquely preferred over an
exponential or another survival model. A calibrated corrosion study would use a
censor-aware fit and compare alternative distributions against digitized data.

## Pit-number, clustering, and finite-size crossover

The nonlinear `by_count` tables lower thinning and rupture barriers when active
neighbors are present, while active neighbors raise repassivation barriers. The
primary sweep therefore records all three Phase-2 qualitative signatures:

- pit-number density rises from 0.0060 to 0.0629 patches⁻¹;
- largest connected pit clusters grow from 1.25 to 7.28 patches;
- the active-fraction susceptibility peaks at +0.55 kcal/mol on the 24² deck,
  identifying the bounded crossover region rather than selecting it by eye.

The finite-size scan is:

| Size | μ | Active fraction | Susceptibility | Pit density | Largest cluster | Largest / area | NN pair ratio |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 16² | +0.30 | 0.0425 | 0.0447 | 0.0359 | 2.13 | 0.00830 | 1.64 |
| 16² | +0.55 | 0.0723 | 0.1010 | 0.0476 | 3.69 | 0.01440 | 2.19 |
| 16² | +0.80 | 0.1108 | 0.1432 | 0.0632 | 5.31 | 0.02075 | 1.95 |
| 24² | +0.30 | 0.0475 | 0.0284 | 0.0350 | 3.94 | 0.00684 | 2.77 |
| 24² | +0.55 | 0.0770 | 0.0687 | 0.0499 | 5.69 | 0.00987 | 2.28 |
| 24² | +0.80 | 0.1108 | 0.0255 | 0.0655 | 7.31 | 0.01270 | 1.94 |
| 32² | +0.30 | 0.0460 | 0.0617 | 0.0366 | 3.31 | 0.00323 | 2.18 |
| 32² | +0.55 | 0.0745 | 0.1412 | 0.0496 | 5.13 | 0.00500 | 2.28 |
| 32² | +0.80 | 0.1044 | 0.1156 | 0.0608 | 7.25 | 0.00708 | 2.06 |

For 24² and 32², susceptibility is largest at +0.55; the 16² system peaks at
+0.80, showing finite-size drift and bracketing the crossover at roughly
+0.55–+0.80 kcal/mol. Largest-cluster fraction decreases with increasing area
at each sampled potential, so this bounded scan does not support a
system-spanning/percolating-cluster claim. It is a qualitative finite-size
scaling check, not a critical-exponent estimate.

The independent regression gate uses eight fixed seeds at a 400 s exact
horizon. It observes low-driving largest clusters
`[1, 1, 1, 1, 1, 2, 1, 2]` with 0/8 initiated and high-driving clusters
`[8, 7, 4, 5, 7, 4, 5, 6]` with 8/8 initiated. The test also requires the mean
pit count to rise by more than 3×, directly guarding the pit-number target.

## Metastable current transients

`metastable-transient.csv` records, every 20 events and at every rupture:
physical CTMC time, active-patch population, largest cluster, the instantaneous
`metal_dissolution` propensity, and cumulative dissolution/repassivation event
counts. A final exact-deadline row is always written. For the fixed +0.30
kcal/mol representative trajectory it contains 627 samples ending at exactly
600.000000000 s. The active population ranges from 0 to 36 and decreases
between successive samples 151 times; the current propensity ranges from 0 to
1.8003 s⁻¹ and decreases 195 times. Early examples include an active population
of six patches at 47.2 s collapsing to one by 53.2 s, and a current propensity
of 0.1825 s⁻¹ at 39.5 s falling to 0.0913 s⁻¹ at 40.8 s. Those rise-and-repair
episodes are the event-rate analogue of metastable pitting transients.

## Deck physics and scope

- **Chloride and potential:** `thermo.activity.Cl` is chloride activity;
  `thermo.mu.e_drive` is the additive electrochemical driving proxy. The
  adsorption rule consumes both reservoir species, and anodic rupture/
  dissolution consume `e_drive`.
- **Autocatalysis:** nonlinear `by_count` barrier tables encode local
  acidification/chloride-retention feedback and repassivation suppression.
- **Defect nucleation:** `[[structure.defects]]` declares a screw-dislocation
  strain field; `film_rupture` has negative strain coupling, so the defect core
  is a mechanically derived weak point rather than a hard-coded special site.
- **Exact timing:** the deck uses Petra's CTMC strategy and exact-deadline API,
  so induction times and current-rate samples are in simulation seconds rather
  than CA generations.
- **Determinism:** regression coverage proves same deck + seed gives an identical
  trajectory while a different seed diverges.
- **Idealization:** sites represent nanometre-scale passive-film patches. There
  is no explicit solution transport, ohmic drop, salt film, alloy chemistry, or
  physical prefactor calibration. The neighbor feedback is a declared proxy,
  and the outputs support qualitative mechanism/conformance claims only.

## Literature anchors

The interpretation follows the source map already vetted in
`docs/scoping/pitting-corrosion.md`: Shibata's stochastic induction-time theory,
Punckt et al.'s cooperative critical onset, and the passive-film
breakdown-versus-repair framing. Primary links carried by that scoping record:

- Shibata stochastic theory: https://content.ampp.org/corrosion/article-abstract/33/7/243/3898/Stochastic-Theory-of-Pitting-Corrosion
- Punckt et al. cooperative onset: https://science.sciencemag.org/content/305/5687/1133
- Valor et al. Weibull/nonhomogeneous-Poisson statistics: https://www.sciencedirect.com/science/article/abs/pii/S0010938X06002150
