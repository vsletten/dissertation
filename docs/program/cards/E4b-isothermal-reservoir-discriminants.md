# E4b-isothermal-reservoir-discriminants — finish §5 time-domain tests

- status: done
- track: E (muscovite / the perfect circle)
- priority: P1
- machine: any (CPU ensembles)
- depends: E4a-surface-gated-release
- claimed-by: hermes-macbot-zero

## Objective

Close the four E4 §5 verdicts that failed because the current deck lacks the
observable or mechanism rather than because a calibrated hypothesis survived:

1. run 500/700 °C isothermal schedules on E2b volumes and invert cumulative
   species release into the same D/a²(t) observables used in 1998 Figures 6–11;
2. test the ³⁶Ar/⁴⁰Ar diffusivity-ratio decay with explicit reservoir-specific
   kinetic access;
3. add a provenance-labelled Ar/Xe paired release hypothesis using E3a's Xe
   screen only as a screen (not a production calibration), so the Villa
   discriminator can reject or prioritize the mechanism; and
4. add an H₂O chemical-potential control and paired 700 °C vacuum/2 kbar
   hydrothermal executions.

## Constraints

- Depends on E4a so unconditional surface loss cannot trivially dominate every
  time-domain test.
- Recover or digitize the 1998 Figures 6–11 time-series support with uncertainty;
  do not invent a furnace-controller log or raw laboratory table.
- Every D/a² inversion states geometry and branch formula and includes a
  forward/inverse unit test against a hand-checkable vector.
- Xe and H₂O parameters remain provenance-labelled screens/proxies until a
  separate accepted barrier/calibration card upgrades them. No fitting knobs.

## Acceptance

- Source-backed cumulative-release and D/a² overlays at 500/700 °C for all three
  E2b volumes, with ensemble bands and deterministic replay.
- Explicit reproduced / partially / not verdicts for §5 claims 1, 2, 4, and 7,
  including the Villa Ar/Xe and hydrothermal-vacuum discriminators.
- Driver, generated decks/schedules, inversion tests, Petra suite, Ruff, card /
  PLAN / STATUS green.

## Progress

- 2026-09-13 12:45 PDT — (hermes-macbot-zero; profile=laptop) — Filed from E4's
  honest NO-GO verdicts. E4 proves the current stepped deck releases 91.6–100%
  of ⁴⁰Ar by 600 °C, has direct ³⁶Ar/⁴⁰Ar release ratios only 0–0.333, and
  contains neither a complete Xe mechanism nor H₂O chemical-potential control;
  those missing capabilities require an explicit follow-up rather than a
  paper-figure claim.
- 2026-09-13 19:43 PDT — (hermes-macbot-zero; profile=laptop) — Completed the
  12-cell isothermal/H₂O campaign with eight replicas per cell and byte-identical
  replay checks. Digitized Figures 6–11 with uncertainty and source hashes;
  added infinite-cylinder D/a² inversion, explicit ³⁶Ar reservoir access,
  provenance-labelled Xe and H₂O-vacancy screens, ensemble overlays, and honest
  mixed verdicts for §5 claims 1, 2, 4, and 7. Result:
  `docs/program/results/E4b-isothermal-reservoir-discriminants.md`.
- 2026-09-13 21:06 PDT — (hermes-custom-build-001; profile=workstation) —
  Independently closed the delivery gates on exact pushed head
  `94be8c6f7b375368efa17179b086870c86bd5343`: 43 Python tests, focused Ruff
  check/format, the complete Petra Cargo suite, and diff checks pass. A fresh
  12-cell/96-primary-trajectory campaign reproduced all four verdicts, all replay
  checks, 121250 primary events, and a byte-identical overlay; the 132-row series
  differs only in 18 cross-Python-version float renderings (maximum absolute delta
  `3.552713678800501e-15`). Corrected the card's non-protocol `complete` status to
  canonical `done` before PR creation.

## Result

- **Outcome:** complete. §5.1, §5.4, and §5.7 are partially reproduced; §5.2
  is not reproduced once diffusivity-ratio boundaries require at least four
  finite replica values.
- **Evidence:** 12 primary cells × 8 replicas and 24 two-replica replay runs;
  all replay file triplets are byte-identical. Source rasters for Figures 6–11
  are hash-receipted and sparse digitizations carry explicit plot-reading
  uncertainty.
- **Verification:** 43 Python tests pass; focused Ruff check and format pass;
  all 12 generated decks complete under the cached Petra release engine with
  `--paranoid`. This laptop has no Cargo toolchain, so the unchanged Rust
  workspace suite is explicitly deferred to required PR CI rather than claimed.
- **Workstation re-verification:** the complete Petra Cargo suite passes against
  the canonical repository target; an independent campaign rerun reproduced the
  verdicts and all replay checks with only sub-femtoscale cross-version float
  serialization differences described in Progress.
- **Artifacts:** driver, 12 generated decks, source CSV, bootstrap-mean 95%
  intervals with effective sample counts, overlay SVG, campaign receipts,
  machine-readable verdict, and full result memo are committed in this PR.
- **Deviation:** the explicit fast-access ³⁶Ar mechanism does not reproduce the
  observed early-high/late-near-unity D/a²-ratio sequence. Xe and H₂O remain
  provenance-labelled screens, not calibrations.
