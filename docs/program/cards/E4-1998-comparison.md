# E4-1998-comparison — synthetic spectra over the founding data

- status: active
- track: E (muscovite / the perfect circle)
- priority: P1
- machine: any (CPU ensembles; needs the E3a barrier set)
- depends: E3a-classical-neb-barriers, E2b-grain-size-sweep, A8-thesis-archive-intake
- claimed-by: hermes-macbot-zero

## Objective

The point of the track: integrate E3a's computed barriers (with E3b
calibration verdicts where available) into the full-mechanism deck, run the
laboratory T(t) schedules from Sletten & Onstott (1998) at E2b-validated
volumes, and lay the synthetic ⁴⁰Ar/³⁹Ar age spectra, cumulative release
curves, and ³⁶Ar/⁴⁰Ar correlations directly over the 1998 target figures
(Figs. 3–11) and Villa's isotope-correlation test. Digitized 1998 release
data comes from the thesis archive intake (A8).

Testable claims on the line (scoping doc §5): staircase spectra from a
homogeneous-age lattice with defect + alteration reservoirs (no fossil
⁴⁰Ar* gradient); the rise-then-fall D/a²(t) signature; species-resolved
reservoir release ordering.

## Constraints

- Every plotted synthetic curve states its barrier provenance (computed vs
  calibrated vs residual proxy) — the comparison is only as honest as its
  labels.
- No fitting knobs beyond the physically carded parameters; if the spectra
  don't match, that result is reported as-is (NO-GO is a finding).
- This card produces the paper-1 figure set; manuscript work is a separate
  future card.

## Acceptance

- Synthetic-vs-1998 overlay figures for the untreated and hydrothermally
  treated splits, with ensemble bands, landed under
  `docs/program/results/` with the generating decks/schedules committed.
- An explicit verdict per §5 claim: reproduced / partially / not — with the
  mechanism attribution for each.
- Petra tests, replay gates, lint green; card/PLAN/STATUS bookkeeping.

## Progress

- 2026-09-13 13:43 PDT (hermes-macbot-zero; profile=laptop) — Implementation,
  ensemble evidence, docs, and follow-up cards are verified and pushed at
  `agents/E4-1998-comparison@13128786e3371c531d9d60596a4454378928fd7f`.
  PR creation is the only remaining lifecycle step. Existing `gh` auth and the
  macOS Keychain GitHub credential both return HTTP 401, and the authenticated
  Chrome window is on an AX-unreachable Space; no API credential was configured
  or exposed. Returning the queue task as a continuation for a browser/API-capable
  fleet worker to open the compare URL and complete POLICY §2 teardown.
- 2026-09-13 12:45 PDT (hermes-macbot-zero; profile=laptop) — Completed the
  E3a-hop × retained-delamination-sensitivity × E2b-volume campaign: six
  eight-replica ensembles plus six byte-identical two-replica replay pairs (72
  total replicas, 62,709 primary events). The source-backed Figure 3–5
  digitizations and overlays show a mixed verdict: recoil distortion and
  high-temperature merge reproduce; only one volume per delamination sensitivity
  has a materially old initial step; the observed 700→800→1025 °C grain-size
  crossover fails because the low-sensitivity family peaks at 500 °C and the
  high-sensitivity family at 600 °C regardless of volume. Full Petra Rust
  workspace tests and all 35 Python script tests pass; Ruff and `cargo fmt` are
  green; `cargo clippy --workspace --all-targets` exits 0 with three pre-existing
  Rust 1.98 lints in `petra-core/src/engine.rs` (lines 42, 58, 100), so the
  stricter `-D warnings` variant remains red for unrelated baseline code.
  Artifact and replay gates are green.
- 2026-09-13 12:29 PDT (hermes-macbot-zero; profile=laptop) — Claimed by atomic push of `agents/E4-1998-comparison` at exact `origin/main@1a15b3141a60df331ad41a21864186317c63fdf2`; all three declared prerequisites are done on main. Beginning the bounded CPU implementation and 1998-data overlay verification in the branch-matched worktree.

- 2026-09-13 11:09 PDT — (hermes-macbot-one; profile=laptop) — Unblocked after
  E3a's bounded result and the completed E2b/A8 prerequisites. E3b remains an
  optional calibration input "where available," not a hard dependency; this
  card's provenance contract still distinguishes computed, calibrated, and
  residual-proxy barriers explicitly.

- 2026-08-27 — filed by Fable. Blocked on the computed-barrier campaign
  (E3a), the statistics ladder (E2b), and thesis-archive intake (A8),
  which supplies the digitized 1998 release data.

## Result

- Verdict: partial mechanism validation, with the grain-size crossover a hard
  falsification of the current always-open surface-release boundary.
- Products: `docs/program/results/E4-1998-comparison.md`, its four committed
  result artifacts, six generated comparison decks, the source-backed
  digitization data, and the tested campaign/analysis driver.
