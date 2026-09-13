# E4b-isothermal-reservoir-discriminants — finish §5 time-domain tests

- status: blocked: E4a
- track: E (muscovite / the perfect circle)
- priority: P1
- machine: any (CPU ensembles)
- depends: E4a-surface-gated-release
- claimed-by: -

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
