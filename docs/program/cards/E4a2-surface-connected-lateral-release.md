# E4a2-surface-connected-lateral-release — test the missing in-plane length scale

- status: ready
- track: E (muscovite / the perfect circle)
- priority: P1
- machine: any (CPU ensembles)
- depends: E4a-surface-gated-release
- claimed-by: -

## Objective

Test the mechanism isolated by E4a's NO-GO: replace independently opening basal
cells with an explicit surface-connected release front that originates at the
open `a/b` edges and advances only through a connected delaminated-interface
path. Rerun the exact E4 campaign to determine whether the existing in-plane
`4×4`, `8×8`, and `12×12` ladder then produces the observed monotonic
700→800→1025 °C release-peak ordering.

## Constraints

- No per-volume, per-grain-size, or per-temperature release parameter. The only
  new control is graph connectivity to a declared lateral edge through existing
  `Interface.delaminated` state.
- Preserve E3a's 64.095991/68.410690 kcal mol⁻¹ incomplete-convergence Ar-hop
  sensitivities, E2's retained 58 ± 5 kcal mol⁻¹ delamination proxy, and E4's
  schedules, seeds, digitization, uncertainty, and qualitative gates.
- Keep residual E2 barriers labelled as proxies. No fit claim.

## Acceptance

- A fail-closed regression proves an interior delaminated pocket remains
  inaccessible until it connects to a lateral edge, while a connected path
  enables isotope release.
- The six E4 ensembles rerun with byte-identical same-seed replay and committed
  comparison products.
- Explicit verdict on the 700→800→1025 °C trend. If still absent, report NO-GO
  and identify the next discriminating mechanism rather than tuning around it.
- Petra Rust and Python suites, Ruff, formatting, card/PLAN/STATUS green.

## Progress

- 2026-09-13 15:05 PDT (hermes-macbot-zero; profile=laptop) — Filed from E4a's
  independently opening basal-gate NO-GO. The current comparison ladder varies
  the `a/b` dimensions but holds the open basal-axis depth at six layers, so a
  surface-connected in-plane front is the next parameter-free mechanism that
  can create the missing size-dependent traversal length.
