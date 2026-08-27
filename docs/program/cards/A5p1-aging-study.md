# A5p1-aging-study — first science on the aging observables

- status: ready
- track: A (geochemistry)
- priority: P1
- machine: any (CPU ensembles)
- depends: A5p0-aging-observables
- claimed-by:

## Objective

Run the aging study M5 names: use A5p0's long-time rate spectra,
projected-geometric vs BET-like surface accounting, and per-site exposure
ages to characterize how kaolinite dissolution kinetics evolve as fast
sites exhaust — the field-lab discrepancy mechanism the track exists to
test (lab rates measure fresh surfaces; field rates measure aged ones).

1. Ensemble campaigns on the finite-defect kaolinite deck across defect
   densities and durations well past the fast-site exhaustion knee
   (A5p0 smoke: 33 fast sites gone by step 40, bulk propensity −287×).
2. Emit rate-vs-exposure-age curves with ensemble bands; quantify the
   apparent-rate drop between fresh-surface and aged-surface regimes and
   its dependence on defect density.
3. State the comparison against the published field-lab discrepancy range
   (orders of magnitude) — does surface aging of this mechanism class
   plausibly account for it, partially, or not?

## Constraints

- Current barrier set (with its provenance) is fine — this is a
  mechanism-class study, not a calibrated prediction; label accordingly.
- Deterministic replay and A5p0 regression gates must stay green.

## Acceptance

- Results doc with rate-spectrum evolution, exposure-age distributions,
  and the explicit fresh-vs-aged apparent-rate ratio with bands.
- An honest verdict paragraph on the field-lab question at this model
  tier.
- Petra tests, lint green; card/PLAN/STATUS bookkeeping in the PR.

## Progress

- 2026-08-27 — filed by Fable: A5p0's observables landed 2026-08-22 but the
  study itself was never carded; M5 lists it as first-science.
