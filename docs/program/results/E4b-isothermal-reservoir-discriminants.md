# E4b — Isothermal reservoir discriminants

**Status:** complete — source-backed discriminant campaign, not a parameter fit
**Board card:** [E4b](../../cards/E4b-isothermal-reservoir-discriminants.md)
**Source:** [Sletten & Onstott (1998)](https://doi.org/10.1016/S0016-7037(97)00323-2), especially Figures 6–11

## Bottom line

The E4a surface-gated mechanism now has the missing isothermal observables and
three explicit, provenance-labelled discriminator hypotheses. A 12-cell CPU
campaign (500/700 °C × vacuum/H₂O screen × three E2b volumes), eight replicas
per cell, completed with byte-identical two-replica replay pairs.

The result is useful because it is mixed rather than cosmetically positive:

- **§5.1 two-stage non-Fickian release: partially reproduced.** All three 500 °C
  volumes show the requested rise-then-fall apparent cylinder D/a²; none of the
  700 °C cells meets that shape test.
- **§5.2 distinct ³⁶Ar/⁴⁰Ar reservoirs: not reproduced.** After requiring at
  least four finite replica ratios at a boundary, no temperature/volume cell
  produces both the measured early 10–10⁵ ratio and a supported late-time
  approach to unity. Single-replica apparent matches are explicitly excluded.
- **§5.4 Ar/Xe discriminator: partially reproduced as a screen only.** The E3a
  Xe barrier produces a paired release prediction and clear Ar/Xe decoupling,
  but it cannot be called a calibrated reproduction without the Villa-series
  muscovite Xe time series.
- **§5.7 hydrothermal contrast: partially reproduced.** At 700 °C in the largest
  volume, final ⁴⁰Ar release falls from **0.4894** in vacuum to **0.4394** in the
  ideal 2-kbar H₂O-vacancy screen. The direction is right; the activity mapping
  is deliberately not fitted.

The machine-readable verdict and every confidence interval are in
[`verdict.json`](E4b-isothermal-reservoir-discriminants/verdict.json) and
[`isothermal-series.csv`](E4b-isothermal-reservoir-discriminants/isothermal-series.csv).

## Source-backed targets

[`figures6-11-time-series.csv`](../../../petra/data/muscovite-1998/figures6-11-time-series.csv)
contains sparse manual digitizations with explicit plot-reading uncertainty:

- Figures 6 and 7: 3000 µm, 500 and 700 °C cumulative release and the paper's
  infinite-cylinder D/a² branch;
- Figure 8: dehydroxylation K/a² and Ar-loss D/a² comparison;
- Figures 9 and 10: 150–180 and 52–64 µm, 500 °C cumulative release and
  cylinder D/a²; and
- Figure 11a: log₁₀[D/a²(³⁶Ar)/D/a²(⁴⁰Ar*)] by run and grain size.

The source README records direct publisher raster URLs and SHA-256 receipts for
Figures 6–11. The raster images are not redistributed. These are image-derived
point estimates, not recovered GD150 instrument files, and cannot support a
formal residual fit.

The generated overlay is
[`isothermal-overlays.svg`](E4b-isothermal-reservoir-discriminants/isothermal-overlays.svg).
Black points are the 1998 digitization; colored means and translucent 95% bands
are the three model volumes. The 4×4×6, 8×8×6, and 12×12×6 volumes are compared
with the three experimental size fractions as a resolution ladder; no physical
lattice-to-micrometre conversion is claimed.

## Mechanisms tested

| Hypothesis | Implementation | Provenance / limit |
|---|---|---|
| Radiogenic-Ar local access | pristine ⁴⁰Ar/³⁹Ar/³⁶Ar hops use 64.095991 kcal mol⁻¹ | E3a local-dehydroxylated Ar screen; incomplete-convergence, not calibrated |
| Radiogenic-Ar extended-zone access | extended ⁴⁰Ar/³⁹Ar hops use 68.410690 kcal mol⁻¹ | E3a interlayer-vacancy Ar screen; incomplete-convergence |
| Atmospheric-Ar fast reservoir | extended-zone ³⁶Ar retains E2's 40 kcal mol⁻¹ proxy | explicit kinetic access; tests the 1998 reservoir interpretation instead of relying on initialization |
| Xe paired release | Xe states use 90.744700 kcal mol⁻¹ hops plus the same surface gate | E3a Xe screen; incomplete-convergence; designed to be rejected/prioritized by external Villa-series data |
| H₂O control | dehydroxylation consumes an explicit `H2O_vacancy` reservoir; activity is 1 in vacuum and 1/2000 for the ideal 2-kbar comparison | directionality screen only; not a water fugacity equation or 2-kbar calibration |

No parameter above was tuned against Figures 6–11. E3a values retain their
screen labels, and the E2 40/58 kcal mol⁻¹ terms remain proxies.

## D/a² inversion

For every replica and isotope, cumulative loss `F` is inverted through the same
infinite-cylinder solution already implemented in `muscovite_release.py`:

`F = 1 - 4 Σ exp(-αₙ² τ)/αₙ²`, where `τ = Dt/a²`.

The per-interval apparent diffusivity is
`(τᵢ - τᵢ₋₁)/(tᵢ - tᵢ₋₁)`. Values at complete release are left undefined rather
than extrapolated through the singular end member. The test suite includes the
hand-checkable short-time vector `τ=10⁻⁴`,
`F=0.022467583341910253`, which inverts back to `τ` and gives
`D/a²=10⁻⁶ s⁻¹` over 100 s.

## Campaign receipt

- **Matrix:** 12 primary cells, eight replicas each; seeds 19981–19988.
- **Schedule:** eleven cumulative boundaries at √time = 2, 4, 6, 8, 10, 15,
  20, 30, 40, 50, and 60 min½ (ending at 216000 s).
- **Replay:** every cell rerun twice at two replicas; `ensemble.csv`,
  `ensemble-summary.csv`, and `observables.csv` were byte-identical.
- **Primary event count:** 121250 events across the 96 primary trajectories.
- **Runtime:** individual cells completed in 0.14–0.28 s on the recorded machine.
- **Receipt:** [`campaign-receipts.json`](E4b-isothermal-reservoir-discriminants/campaign-receipts.json).
- **Raw outputs:** `/private/tmp/e4b-campaign-20260913` (ephemeral, intentionally
  not committed); all analysis inputs needed for review are represented by the
  committed decks, seeds, receipt, digitization, and driver.

## Interpretation and negative evidence

The 500 °C rise-then-fall result demonstrates that coupled dehydroxylation,
interface opening, and finite release access can create the qualitative
non-Fickian signature. Failure at 700 °C matters: the current gates open too
quickly there, collapsing the two-stage sequence into an early release pulse.
That is a mechanism falsification signal, not a request to tune a barrier.

The reservoir-ratio result is negative. The largest 500 °C cell's superficially
high first ratio (21.21 at √time=2) has effective `n=1` and is rejected. With
the fixed `n≥4` gate, that cell begins at 0.75 and has no supported boundary
later than √time=15; the other cells likewise lack the required high-early plus
near-unity-late sequence. The small cells often contain only a handful of ³⁶Ar
atoms, producing undefined or noisy interval ratios. A future quantitative
claim requires larger species inventories and source-data recovery, not merely
more replicas of the same tiny state space.

The Xe and H₂O screens now make the intended external tests executable, but not
settled. The next scientifically defensible step is to obtain the Villa-series
Ar/Xe time series and use a proper H₂O fugacity/activity model. Until then both
claims remain **partial** by construction.

## Reproduction

```bash
python3.14 petra/scripts/muscovite_isothermal_discriminants.py generate \
  petra/decks/muscovite-isothermal

python3.14 petra/scripts/muscovite_isothermal_discriminants.py run \
  petra /path/to/petra petra/decks/muscovite-isothermal \
  /private/tmp/e4b-campaign

python3.14 petra/scripts/muscovite_isothermal_discriminants.py analyze \
  /private/tmp/e4b-campaign petra/decks/muscovite-isothermal \
  petra/data/muscovite-1998/figures6-11-time-series.csv \
  docs/program/results/E4b-isothermal-reservoir-discriminants
```

## Verification

Executed locally:

- `python3.14 -m unittest discover -s petra/scripts/tests -p 'test_*.py'`
  — 43/43 pass;
- focused `ruff check` and `ruff format --check` on the new driver/tests — pass;
- full 12-cell `--paranoid` campaign, byte-replay verification, and analysis
  regeneration — pass.

Python 3.13 and 3.14 are both supported. The reproduction commands and local
record above use 3.14 (the Mac that ran the campaign). PR verification used
Python 3.13. The regenerated overlay is byte-identical across those versions.
The 132-row series has 18 Python-version float-rendering differences with
maximum absolute delta `3.552713678800501e-15`; semantic JSON/CSV comparison
passes at `1e-14` relative tolerance.

This Mac does not have a Cargo toolchain. No Rust source changed, and the cached
release engine compiled every generated deck and executed the full campaign,
but the Rust workspace suite is honestly deferred to mandatory PR CI. The
repository-wide Ruff command also reports pre-existing findings in
`aging_study.py` and `muscovite_release.py`; the changed files are clean.
