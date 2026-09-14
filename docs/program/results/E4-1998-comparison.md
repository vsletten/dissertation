# E4 — Sletten & Onstott (1998) comparison

**Card:** `E4-1998-comparison`
**Scope:** qualitative model discrimination, not parameter fitting

## Bottom line

The E2 mechanism with E3a's bounded local-dehydroxylate/extended-zone Ar-hop
estimates and E2's retained 58 ± 5 kcal mol⁻¹ delamination proxy reproduces
recoil-driven spectrum distortion and the eventual high-temperature merge. It
does **not** reproduce the grain-size crossover, and it produces a materially
old initial step in only one volume per delamination sensitivity, without a
monotonic size trend.

The failure is diagnostic, not cosmetic. In both retained delamination-proxy
variants and all three volumes, the synthetic ⁴⁰Ar release-rate maximum occurs
at 500 °C (53 kcal mol⁻¹) or 600 °C (63 kcal mol⁻¹), with no size ordering. The
1998 Figure 4 digitization instead moves the maximum from
about 700 °C (`52–64 µm`) through 800 °C (`150–180 µm`) to about 1025 °C
(`3000 µm`). Petra's always-open basal `Surface_gate` empties most ⁴⁰Ar before
the structural sequence can establish that crossover. The next mechanism
revision must make surface accessibility structural/temperature-dependent;
retuning the already-fast E3a delamination barrier cannot repair it. Follow-up
card `E4a-surface-gated-release` owns that tested mechanism change and exact
campaign rerun. `E4b-isothermal-reservoir-discriminants` then owns the explicit
isothermal D/a², ³⁶Ar/⁴⁰Ar, Ar/Xe, and hydrothermal-vacuum discriminators whose
current §5 verdicts are NO-GO because the required mechanisms/observables are
absent.

## Inputs and provenance

- E2 full mechanism: `petra/scripts/build_muscovite_full_deck.py`.
- E3a bounded classical estimates: 64.095991 kcal mol⁻¹ for the local
  dehydroxylate Ar hop and 68.410690 kcal mol⁻¹ for the lowest extended-zone
  geometry. Both remain `incomplete-convergence`; they are sensitivity values,
  not calibrated production barriers.
- Delamination: E2's 58 ± 5 kcal mol⁻¹ proxy (53/63 kcal mol⁻¹ variants) is
  retained because E3a ended `incomplete-input-contract` and emitted no
  delamination barrier. The 64 kcal mol⁻¹ octahedral-trap proxy is likewise
  retained after E3a's `incomplete-index-gate`.
- E2b comparison volumes: `4×4×6`, `8×8×6`, and `12×12×6` cells
  (768, 3,072, and 6,912 Petra sites).
- Experimental comparison points: manual image digitizations of Sletten &
  Onstott (1998) Figures 3–5, with plot-reading uncertainty and exact source
  image hashes recorded in `petra/data/muscovite-1998/README.md`.
- Comparison schedule: 500–1200 °C in 50 °C, 600 s increments. This matches the
  Figure 4 temperature support and preserves E2's dwell convention. It is **not**
  represented as a recovered furnace-controller log; A8 contains no such table.

The sparse digitization is adequate for sign/shape tests only. No interpolation,
least-squares score, or calibrated age fit is claimed.

## Ensemble and replay evidence

Six tracked decks cover the Cartesian product of two retained delamination
sensitivities and three E2b volumes. Each ran eight ordered seeds
(`19980`–`19987`) with `--paranoid`, plus
two independent two-replica same-seed replays. All six replay comparisons are
byte-identical for `ensemble.csv`, `ensemble-summary.csv`, and
`observables.csv`.

| Evidence | Value |
|---|---:|
| Primary ensemble replicas | 48 |
| Replay replicas | 24 |
| Total executed replicas | 72 |
| Primary KMC events | 62,709 |
| Primary elapsed time | 1.245 s |
| Replay checks | 6 / 6 byte-identical |

Raw logs and per-replica trajectories remain in the bounded host scratch root
`/private/tmp/e4-campaign-20260913-corrected`; the compact receipts and derived products
are committed here.

## §5 claim-by-claim verdict

| §5 claim | Experimental discriminator | E3a-hop / delamination-sensitivity result | Verdict |
|---|---|---|---|
| 1. Dehydroxylation-limited, two-stage non-Fickian loss and rise-then-fall D/a²(t) | 1998 Figs. 6–10 | The low sensitivity releases 100% of ⁴⁰Ar by 600 °C; the high sensitivity releases 91.6–94.4% by then. Neither produces a resolvable late stage or rise-then-fall inversion | **Not reproduced** — unconditional surface release outruns the structural sequence |
| 2. Distinct reservoirs produce ³⁶Ar/⁴⁰Ar diffusivity-ratio decay from 10–10⁵ toward 1 | 1998 Fig. 11a | The direct synthetic ³⁶Ar/⁴⁰Ar release ratio spans only 0–0.333 and the present analysis has no species-specific D/a² inversion | **Not reproduced** — distinct initialization alone does not create the measured kinetic separation |
| 3. Structure generates a staircase from homogeneous-age material | 1998 Fig. 3d hydrothermal staircase | All six ensembles generate ordered age steps without a fossil gradient, but the maximum synthetic age is only 56.7 Ma against the roughly 400–910 Ma digitized staircase | **Partial** — shape mechanism exists, magnitude/late structure do not |
| 4. Ar/Xe decoupling arbitrates delamination vs in-layer siting | Villa 2021–2026 | E2 has no Xe state/rules; E3a's 90.7447 kcal mol⁻¹ Xe screen is not a complete paired release mechanism | **Not reproduced** — discriminator is absent rather than silently proxied |
| 5. Recoil siting explains anomalously old first steps | 1998 Fig. 5 | Low-sensitivity excesses are 55.5, 44.1, and 40.2 Ma across the size ladder; high-sensitivity excesses are 32.2, 53.8, and 42.1 Ma. Every run shows early-age decline with later ³⁹Ar-rich release | **Partial** — recoil distortion reproduces, but only one volume per sensitivity crosses the 50 Ma old-step gate |
| 6. Grain size controls the delaminated/unstable release fraction | 1998 Fig. 4/11b peak shifts 700 → 800 → 1025 °C | All low-sensitivity profiles peak at 500 °C; all high-sensitivity profiles at 600 °C | **Not reproduced** — always-open `Surface_gate` erases the size dependence |
| 7. Hydrothermal H₂O suppresses breakdown-controlled release | 1998 700 °C/2 kbar split and Harrison 2009 | No H₂O chemical-potential term or paired hydrothermal execution exists; the synthetic curves sit far below the Fig. 3d age staircase | **Not reproduced** — required thermodynamic control is absent |

A secondary expected behavior, the high-temperature merge, **is reproduced**:
final cumulative ⁴⁰Ar release is 1.0 in every volume and barrier bracket. The
53 and 63 kcal mol⁻¹ delamination sensitivities move the entire synthetic family
together from a 500 to a 600 °C peak; neither creates size ordering because
direct surface release still dominates. The missing crossover is therefore not
resolved by choosing one retained delamination sensitivity.

## Products

- `petra/decks/muscovite-1998/*.toml` — six deterministic generated decks.
- `petra/scripts/muscovite_1998_comparison.py` — generation, execution, replay,
  analysis, verdict, and dependency-free SVG writer.
- `petra/data/muscovite-1998/` — source-backed sparse digitizations and limits.
- `docs/program/results/E4-1998-comparison/synthetic-spectra.csv` — 90
  barrier/volume/temperature rows.
- `docs/program/results/E4-1998-comparison/comparison.json` — machine-readable
  seven-claim verdict plus ensemble-derived diagnostics.
- `docs/program/results/E4-1998-comparison/comparison.svg` — normalized release
  profile and age-spectrum overlays.
- `docs/program/results/E4-1998-comparison/campaign-receipts.json` — seed,
  event-count, timing, and replay receipts.

## Reproduction

```bash
# If Rust is not installed, use an isolated toolchain; no API credentials needed.
PATH=/private/tmp/e4-cargo/bin:$PATH \
CARGO_HOME=/private/tmp/e4-cargo RUSTUP_HOME=/private/tmp/e4-rustup \
  cargo build --release -p petra-cli

python3.14 petra/scripts/muscovite_1998_comparison.py generate \
  petra/decks/muscovite-1998
python3.14 petra/scripts/muscovite_1998_comparison.py run \
  petra petra/target/release/petra petra/decks/muscovite-1998 \
  /private/tmp/e4-campaign-20260913-corrected --replicas 8 --base-seed 19980
python3.14 petra/scripts/muscovite_1998_comparison.py analyze \
  /private/tmp/e4-campaign-20260913-corrected petra/decks/muscovite-1998 \
  petra/data/muscovite-1998 docs/program/results/E4-1998-comparison
```

## Verification

- `cargo test --workspace`: pass.
- `cargo fmt --all -- --check`: pass.
- `cargo clippy --workspace --all-targets`: exits 0 with three pre-existing Rust
  1.98 suggestions in `petra-core/src/engine.rs` (lines 42, 58, and 100). The
  strict `-- -D warnings` form fails on those baseline warnings; E4 changes no
  Rust.
- `python3.14 -m unittest discover -s petra/scripts/tests -p 'test_*.py'`: 35
  tests pass.
- Ruff check and format-check: pass.
- Artifact gate: six decks, six receipts, 90 synthetic rows, valid SVG, and all
  replay comparisons verified.
