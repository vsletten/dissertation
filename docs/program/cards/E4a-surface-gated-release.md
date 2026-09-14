# E4a-surface-gated-release — repair the grain-size crossover

- status: done
- track: E (muscovite / the perfect circle)
- priority: P1
- machine: any (CPU ensembles)
- depends: E4-1998-comparison, B5-execution-schedule
- claimed-by: hermes-macbot-zero

## Objective

Close E4's falsification of the current full-mechanism boundary condition. The
always-open basal `Surface_gate` releases most ⁴⁰Ar by 500–600 °C for every E2b
volume, while the Sletten–Onstott (1998) Figure 4 digitization
moves the release-rate maximum from about 700 °C (`52–64 µm`) through 800 °C
(`150–180 µm`) to about 1025 °C (`3000 µm`).

Replace unconditional surface accessibility with a structural gate tied to the
existing dehydroxylation/delamination state machine, then rerun the exact E4
E3a-hop × retained-delamination-sensitivity × volume campaign without adding a grain-size-specific fit
parameter.

## Constraints

- Preserve E3a's 64.095991/68.410690 kcal mol⁻¹ incomplete-convergence Ar-hop
  sensitivities, E2's retained 58 ± 5 delamination proxy, and E4's tracked
  schedule, seeds, digitization, uncertainty, and qualitative gates.
- Surface accessibility must emerge from local state/topology; no per-volume,
  per-grain-size, or per-temperature release multiplier.
- Add a fail-closed regression proving intact basal surfaces cannot fire the
  ⁴⁰Ar/³⁹Ar/³⁶Ar release rules and delaminated accessible surfaces can.
- Keep residual E2 proxy barriers labelled as proxies. No fit claim.

## Acceptance

- A tested deck/runtime representation gates isotope release on explicit local
  accessibility without weakening existing schedule/replay guarantees.
- The six E4 ensembles rerun with byte-identical same-seed replay and committed
  comparison products.
- Explicit verdict on the 700→800→1025 °C grain-size trend. If still absent,
  report NO-GO and identify the next discriminating mechanism; do not tune around
  it.
- Petra Rust and Python suites, Ruff, formatting, card/PLAN/STATUS green.

## Progress

- 2026-09-13 17:07 PDT (hermes-custom-build-001; profile=workstation) — Closed
  the board state after independently verifying the tracked result's **NO-GO**
  bottom line, merged [PR #129](https://github.com/vsletten/dissertation/pull/129),
  and deletion of remote branch `agents/E4a-surface-gated-release`. Follow-ups
  `E4a2-surface-connected-lateral-release` and
  `E4b-isothermal-reservoir-discriminants` are now READY.

- 2026-09-13 15:17 PDT (hermes-macbot-zero; profile=laptop) — CONTINUATION:
  implementation, campaign evidence, tests, and bookkeeping are clean and pushed
  at `agents/E4a-surface-gated-release@8dfb6086e5f88c39ae38812c1813083c63e6199c`;
  opening the PR is the only remaining lifecycle step. Existing `gh` auth is
  invalid, and the authenticated Chrome profile is on an AX-unreachable Space;
  native background and approved foreground escalation could not focus the exact
  window, while typed-browser inspection correctly refused because
  `computer_use.grant_existing_profile` is not enabled. No API credential or
  browser-profile grant was configured.

- 2026-09-13 15:05 PDT (hermes-macbot-zero; profile=laptop) — Replaced the
  unconditional basal boundary with inert/closed/open `Surface_gate` states:
  only basal gates initialize closed, a bonded local delaminated interface opens
  them, and every isotope release rule now requires that explicit open state.
  The exact six-deck campaign completed 72 replicas and 98,264 primary events;
  all six same-seed replay pairs are byte-identical. The grain-size verdict is
  NO-GO: 53 kcal mol⁻¹ peaks are 550/500/500 °C and 63 kcal mol⁻¹ peaks are
  600/600/600 °C across the volume ladder. Follow-up
  `E4a2-surface-connected-lateral-release` owns the next parameter-free test.

- 2026-09-13 12:45 PDT — (hermes-macbot-zero; profile=laptop) — Filed from E4's
  direct experimental comparison. The retained 53 kcal mol⁻¹ delamination
  sensitivity peaks at 500 °C and the 63 kcal mol⁻¹ sensitivity at 600 °C in
  all three E2b volumes, isolating unconditional `Surface_gate` access as the
  next falsifiable mechanism rather than an Ar-hop or delamination choice.

## Result

- Structural gate and regression: `petra/scripts/build_muscovite_full_deck.py`
  and `petra/scripts/tests/test_muscovite_full.py`.
- Campaign evidence: `docs/program/results/E4a-surface-gated-release/` and
  `docs/program/results/E4a-surface-gated-release.md`.
- Verdict: **NO-GO** on the 700→800→1025 °C grain-size trend; local basal
  delamination gating changes release extent but produces no monotonic peak
  ordering. The next discriminator is a surface-connected lateral release front
  through the varying `a/b` dimensions, filed as
  `E4a2-surface-connected-lateral-release`.
