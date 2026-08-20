# E1-muscovite-deck — Ar release from muscovite, phase 1 (published numbers)

- status: done
- track: E (muscovite / the perfect circle)
- priority: P1
- machine: any (KMC is CPU-light; no QM in this phase)
- depends: — (phase 1 is isothermal; staircase spectra later need B5)
- claimed-by: hermes-workstation

## Objective
Phase 1 of docs/scoping/ar-muscovite.md §6: the minimal muscovite deck
using published barriers only. Single-gallery 2M₁ reduction — K/vacancy/Ar40
interlayer gallery with divacancy hops (Em ≈ 66 kcal/mol, Nteme
et al. 2022), paired OH sites with bulk dehydroxylation kinetics
(by_count steric table raising Ea with neighboring O_residual —
Guggenheim mechanism), open c-axis boundary with fast Ar desorption at
the surface. Isothermal runs at 500 °C and 700 °C (the Sletten &
Onstott 1998 in-vacuo conditions). Post-process the event log into
cumulative-release curves and D/a²(t) via the same Fickian inversion
the 1998 paper used.

Read the scoping doc first — it is the design authority for site
kinds, rules, and barrier provenance (§4 table).

## Acceptance
- Deck under petra/decks/ validates and runs on current main (v1
  schema is fine if v2 is not yet merged; note which).
- Release curves at both temperatures show the two-stage non-Fickian
  shape; D/a²(t) inversion shows the rise-then-fall signature
  (scoping §5 claim 1 — qualitative gate, not a fit).
- Results doc in docs/program/results/E1-muscovite-phase1.md with
  plots/data pointers and an honest gap list for phase 2.

## Progress
- 2026-08-19 — card created (fable); scoping doc landed same PR.
- 2026-08-20 — hermes-workstation — claimed atomically on
  `agents/E1-muscovite-deck`; pushed schema-v1 500/700 °C deck checkpoint
  `c977bfa` after both decks compiled and initialized identically.
- 2026-08-20 — hermes-workstation — completed both `--viz --paranoid`
  isothermal trajectories, Hames–Bowring cylinder inversion, machine-readable
  rise-then-fall gates, plots, gap ledger, full Petra suite, and focused
  post-processing tests. DEVIATION: phase 1 uses one tracked deck per fixed
  temperature because B5 schedule/schema-v2 support has not landed.
- 2026-08-20 — hermes-workstation — merged-main re-verification found four
  default-Ruff E731 violations in the SVG helper. Replaced the local lambdas
  with typed nested functions without changing output; the full Petra suite,
  focused post-processing tests, Ruff, and regenerated-result equivalence are
  the follow-up gate.
- 2026-08-20 — hermes-workstation — follow-up gate PASS: 42 Rust tests and
  10 focused Python tests passed; current Rust Clippy (`-D warnings`) and
  Python Ruff check/format passed after clearing four additional toolchain
  findings. Fresh 500/700 °C `--viz --paranoid` runs reproduced 630/831 and
  631/831 release with 2.13×/13.55× rise and 275×/32.3× fall. All five
  regenerated CSV/SVG/JSON products are byte-identical to the tracked
  artifacts, all nine relative report links resolve, and both rasterized SVGs
  were visually inspected with readable, unclipped two-panel PASS plots.
- 2026-08-20 — hermes-workstation — a fresh merged-main Ruff 0.15.1 format
  check caught the cloud review-fix commit's long `sy()` expression. Reformatted
  that expression only; no model, gate, or artifact semantics changed.
- 2026-08-20 — hermes-workstation — exact-main closeout after PR #48 found
  that the earlier pushed `strain_gates.rs` identity-op correction was absent
  from the squash-merged tree. Restored the one-line Clippy fix on a fresh
  follow-up branch; the five strain gates and workspace all-target Clippy pass.

## Result

- status: done
- Decks: `petra/decks/muscovite-phase1-{500c,700c}.toml` (Petra schema v1).
- Analysis: `petra/scripts/muscovite_release.py`; Hames–Bowring
  infinite-cylinder cumulative-loss inversion at 2%-release crossings.
- Durable report/artifacts:
  `docs/program/results/E1-muscovite-phase1.md` and
  `docs/program/results/E1-muscovite-phase1/` (CSV, SVG, JSON verdict).
- 500 °C: 630/831 Ar released (75.8%); apparent `D/a²` rises 2.13× then
  falls 275× from its early peak — PASS.
- 700 °C: 631/831 Ar released (75.9%); apparent `D/a²` rises 13.55× then
  falls 32.3× from its early peak — PASS.
- Verification: both complete runs exited at the truthful no-possible-event
  state under `--paranoid`; `cargo test` passed 42 tests; focused Python
  suite passed 10 tests; Clippy and Ruff passed; generated SVGs were rasterized
  and visually checked. No runtime activation applies; dissertation has no
  installed service for these paths.
- No quantitative fit is claimed. Every phase-2 limitation and proxy is
  enumerated in the result doc. Executable continuation is captured by
  `E2-muscovite-full-mechanism`; the numerical-range engine defect found
  during prototyping is captured by `QI3-fenwick-total-cancellation`; B5 is
  the existing schedule dependency.
