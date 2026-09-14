# E4a2-surface-connected-lateral-release — test the missing in-plane length scale

- status: done
- track: E (muscovite / the perfect circle)
- priority: P1
- machine: any (CPU ensembles)
- depends: E4a-surface-gated-release
- claimed-by: hermes-macbot-zero

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

- 2026-09-13 20:17 PDT (hermes-macbot-one; profile=laptop) — Cold review found
  that the replay receipt exposed only a boolean assertion. Added byte-count and
  SHA-256 evidence for the three primary/replay summary CSVs, cross-checked the
  eight-replica primary prefixes against both same-seed replays, bound analysis
  inputs to committed products, and added divergence/drift tests. An independent
  72-replica run on this host reproduced the exact 84,585-event total and
  unchanged NO-GO comparison. The compiled gate regression now proves a
  two-interior-hop path from both opposite `a` edges; the deck regression pins
  both axes and all four lateral faces. Also
  documented that the universal 1.0 s⁻¹ gate timescale is inherited unchanged
  from E4a rather than introduced or fitted by E4a2; this result excludes the
  connectivity intervention under that inherited timescale, not every possible
  propagation law.

- 2026-09-13 19:57 PDT (hermes-macbot-one; profile=laptop) — Opened
  [PR #131](https://github.com/vsletten/dissertation/pull/131) from the prepared
  pushed branch after confirming the committed user-visible outcome: the
  connected lateral front still yields 500/500/500 °C and 600/600/600 °C
  volume-ladder peaks rather than the observed 700→800→1025 °C ordering.
  Marked card and board DONE in the PR; cloud Hermes now owns review, CI, and
  merge.

- 2026-09-13 18:39 PDT (hermes-macbot-zero; profile=laptop) — Implementation,
  campaign evidence, final verification, and independent review are complete and
  pushed, but PR creation could not finish. Existing `gh` auth is invalid;
  typed-browser access correctly refused without the user-controlled
  `computer_use.grant_existing_profile` opt-in; and native Chrome control failed
  because the authenticated window is on another Space and AX-unreachable. No
  API credential or browser-profile grant was configured. Continue from branch
  `agents/E4a2-surface-connected-lateral-release`: open the PR, record it, and
  perform POLICY §2 teardown. Scientific verdict remains **NO-GO**.

- 2026-09-13 18:14 PDT (hermes-macbot-zero; profile=laptop) — Implemented the
  parameter-free lateral front: all `a/b` boundaries are open, edge gates seed
  only after local delamination, and interior gates require both local
  delamination and an already-open `lateral_front` neighbor. The exact six-deck
  campaign completed 48 primary plus 24 replay replicas with all six replay
  pairs byte-identical. The result is **NO-GO**: the 53 kcal mol⁻¹ proxy peaks
  at 500/500/500 °C and the 63 kcal mol⁻¹ proxy at 600/600/600 °C, not the
  observed 700→800→1025 °C sequence. E4b's isothermal reservoir discriminants
  are the next test; no front parameter was tuned.

- 2026-09-13 15:05 PDT (hermes-macbot-zero; profile=laptop) — Filed from E4a's
  independently opening basal-gate NO-GO. The current comparison ladder varies
  the `a/b` dimensions but holds the open basal-axis depth at six layers, so a
  surface-connected in-plane front is the next parameter-free mechanism that
  can create the missing size-dependent traversal length.

## Result

- 2026-09-13 19:57 PDT (hermes-macbot-one; profile=laptop) — Published the
  completed slice as [PR #131](https://github.com/vsletten/dissertation/pull/131).
- Implementation: `petra/scripts/build_muscovite_full_deck.py`, the tracked full
  deck, and all six tracked E4 campaign decks.
- Regressions: `petra/crates/petra-deck/tests/surface_connected_release.rs`
  executes controlled isolated and connected lattices through the compiled
  engine; `test_surface_release_requires_edge_connected_delamination_front`
  pins the generated campaign-deck structure.
- Evidence: `docs/program/results/E4a2-surface-connected-lateral-release.md` and
  its committed `campaign-receipts.json`, `analysis-receipt.json`,
  `comparison.json`, `comparison.svg`, and `synthetic-spectra.csv` products.
- Verification: 42 Python script tests passed; the Petra Rust workspace tests and
  doc-tests passed; changed-file Ruff check/format and `git diff --check` passed.
- Verdict: **NO-GO** on the 700→800→1025 °C trend for the inherited E4a gate
  timescale, tracked schedule, and proxy-size ladder. The next discriminator is
  the already-READY `E4b-isothermal-reservoir-discriminants` card.
