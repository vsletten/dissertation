# E4a2 — Surface-connected lateral release

**Card:** `E4a2-surface-connected-lateral-release`
**Scope:** parameter-free connectivity test; qualitative discrimination, not fitting

## Bottom line

**NO-GO.** Requiring release access to propagate inward from the open `a/b`
edges through a connected chain of `Interface.delaminated` cells does not
recover the observed 700→800→1025 °C grain-size ordering. The retained 53 kcal
mol⁻¹ delamination proxy peaks at 500 °C for all `4×4×6`, `8×8×6`, and
`12×12×6` volumes; the 63 kcal mol⁻¹ proxy peaks at 600 °C for all three.

The explicit lateral front therefore does not supply the missing experimental
length-scale response under the tracked E4 schedule. The next discriminating
mechanism is E4b's isothermal reservoir-kinetics test: separate reaction-limited
delamination from diffusion-limited Ar transport using time-domain observables,
rather than tuning this front or adding a grain-size-specific release rate.

## Surface-connected front

The comparison decks now use open `a`, `b`, and `c` boundaries. Every
`Surface_gate` starts `closed`; the four lateral faces are initialized as
`edge`. A lateral edge gate opens only after its bonded local interface is
`delaminated`. An interior gate opens only when both its local interface is
`delaminated` and a nearest-neighbor gate connected by a `lateral_front` bond is
already `open`. Isotope release remains possible only beside an explicitly
`open` gate.

The engine-level regression
`petra/crates/petra-deck/tests/surface_connected_release.rs` constructs two
controlled three-cell lattices and executes the compiled rules. In the first,
an isolated interior delaminated pocket enables no event and retains its ⁴⁰Ar.
In the second, delaminating the adjacent lateral-edge interface fires exactly
edge-open → interior-front-advance → isotope-release. The Python regression
`test_surface_release_requires_edge_connected_delamination_front` separately
checks that the generated campaign decks carry those exact boundaries, seeds,
bonds, guards, and release rule. No per-volume, per-grain-size, or
per-temperature parameter was added.

## Preserved inputs

- E3a Ar-hop sensitivities: 64.095991 and 68.410690 kcal mol⁻¹, retained as
  `incomplete-convergence` values.
- E2 delamination sensitivity: 58 ± 5 kcal mol⁻¹, retained and labelled as a
  proxy rather than computed kinetics.
- E4 comparison ladder: `4×4×6`, `8×8×6`, `12×12×6`.
- E4 schedule: 500–1200 °C in 50 °C increments, 600 s per increment.
- E4 seeds and statistics: eight primary replicas per deck, seeds
  `19980`–`19987`, plus two independent two-replica same-seed replays.
- Sletten–Onstott (1998) digitization, uncertainty, and qualitative gates are
  unchanged.

## Campaign and replay evidence

The exact six-deck Cartesian campaign was rerun: two retained delamination
sensitivities × three E2b volumes.

| Evidence | Value |
|---|---:|
| Primary ensemble replicas | 48 |
| Replay replicas | 24 |
| Total executed replicas | 72 |
| Primary KMC events | 84,585 |
| Primary elapsed time | 1.462 s |
| Same-seed replay checks | 6 / 6 byte-identical |

Byte identity holds for `ensemble.csv`, `ensemble-summary.csv`, and
`observables.csv` in every replay pair. Raw logs and trajectories are retained
at `/private/tmp/e4a2-campaign-20260913-1815-r2`; compact receipts and comparison
products are committed under
`docs/program/results/E4a2-surface-connected-lateral-release/`.

## Grain-size verdict

| Delamination sensitivity | `4×4×6` peak | `8×8×6` peak | `12×12×6` peak | 700→800→1025 °C trend |
|---|---:|---:|---:|---|
| 53 kcal mol⁻¹ proxy | 500 °C | 500 °C | 500 °C | **NO-GO** — no size ordering |
| 63 kcal mol⁻¹ proxy | 600 °C | 600 °C | 600 °C | **NO-GO** — no size ordering |

All six families reach a final cumulative ⁴⁰Ar release fraction of 1.0 under
the tracked schedule. The lateral front changes the explicit path to release,
but its traversal remains too fast relative to the retained delamination and
hop kinetics to shift the temperature of maximum release across this size
ladder. That is the model result; it was not tuned away.

## Products

- `petra/scripts/build_muscovite_full_deck.py` and the tracked full deck;
- six regenerated `petra/decks/muscovite-1998/*.toml` campaign decks;
- `petra/scripts/muscovite_1998_comparison.py` with E4a2 verdict wording;
- `synthetic-spectra.csv`, `comparison.json`, `comparison.svg`, and
  `campaign-receipts.json` in this result directory.

## Verification

- Python script suite: 38 tests passed.
- Petra Rust workspace: all tests and doc-tests passed, including the compiled
  isolated-pocket/connected-path behavioral regression.
- Changed Python files: Ruff check and format check passed.
- `git diff --check`: passed.

## Reproduction

```bash
python3.14 petra/scripts/build_muscovite_full_deck.py \
  --output petra/decks/muscovite-full-mechanism.toml
python3.14 petra/scripts/muscovite_1998_comparison.py generate \
  petra/decks/muscovite-1998
python3.14 petra/scripts/muscovite_1998_comparison.py run \
  . "$CARGO_TARGET_DIR/release/petra" petra/decks/muscovite-1998 \
  /private/tmp/e4a2-campaign-20260913-1815-r2 --replicas 8 --base-seed 19980
python3.14 petra/scripts/muscovite_1998_comparison.py analyze \
  /private/tmp/e4a2-campaign-20260913-1815-r2 \
  petra/decks/muscovite-1998 petra/data/muscovite-1998 \
  docs/program/results/E4a2-surface-connected-lateral-release
```
