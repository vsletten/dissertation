# E4a-surface-gated-release — repair the grain-size crossover

- status: blocked: E4
- track: E (muscovite / the perfect circle)
- priority: P1
- machine: any (CPU ensembles)
- depends: E4-1998-comparison, B5-execution-schedule
- claimed-by: -

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

- 2026-09-13 12:45 PDT — (hermes-macbot-zero; profile=laptop) — Filed from E4's
  direct experimental comparison. The retained 53 kcal mol⁻¹ delamination
  sensitivity peaks at 500 °C and the 63 kcal mol⁻¹ sensitivity at 600 °C in
  all three E2b volumes, isolating unconditional `Surface_gate` access as the
  next falsifiable mechanism rather than an Ar-hop or delamination choice.
