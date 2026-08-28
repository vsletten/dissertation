# A5p1-aging-study — first science on the aging observables

- status: done
- track: A (geochemistry)
- priority: P1
- machine: any (CPU ensembles)
- depends: A5p0-aging-observables
- claimed-by: hermes-8d9e30fee91548

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

- 2026-08-28T04:02:39Z — (hermes-8d9e30fee91548; profile=laptop) —
  Acceptance complete. Four defect-density campaigns (0.05/0.15/0.25/0.50),
  each with 64 deterministic replicas, reached the finite slab's honest
  no-events stop. Step 140 is defect-free in every replica. Mean fresh/aged
  apparent-rate drops rise monotonically from 36.36× to 363.06×; rate spectra
  transiently broaden and then collapse onto the slow channel, while the
  BET-like area proxy rises despite the rate decline. Full results, bands,
  provenance, and the model-tier verdict are in
  `docs/program/results/A5p1-aging-study.md`.
- 2026-08-27 — filed by Fable: A5p0's observables landed 2026-08-22 but the
  study itself was never carded; M5 lists it as first-science.

## Result

- `petra/scripts/aging_study.py` provides a bounded, deterministic campaign
  runner and bootstrap analysis; `petra/scripts/tests/test_aging_study.py`
  covers deck rendering, aggregation/area normalization, and the defect-free
  aged-regime gate.
- Results: `docs/program/results/A5p1-aging-study.md`, with committed summary
  CSVs and SVG under `docs/program/results/a5p1-aging-study/`.
- Across 256 trajectories, the mean unnormalized drop spans 36.36×
  `[33.42, 39.28]` to 363.06× `[352.55, 373.41]`; the moderate/high-density
  runs cover 1.98–2.56 log units, supporting a partial-mechanism verdict at
  this explicitly uncalibrated tier.
- External evidence: 20 files / 340,515,688 bytes at
  `/opt/data/run-outputs/dissertation/A5p1-aging-study-20260828/`; manifest
  SHA-256 `e32d658b1857754aab06b5796c8a97a25221fc510621c24d0855876a5a674d58`.
- Verification: 23 Python script tests and 116 Petra workspace tests pass;
  Clippy with `-D warnings`, Rustfmt, deterministic parity, and SVG XML gates
  are green.
