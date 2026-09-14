# E4a — Delamination-gated surface release

**Card:** `E4a-surface-gated-release`  
**Scope:** structural boundary-condition test; qualitative discrimination, not fitting

## Bottom line

**NO-GO.** Replacing the always-open basal boundary with explicit local
accessibility changes the release inventories but does not recover the observed
700→800→1025 °C grain-size crossover. With the retained 53 kcal mol⁻¹
(delamination-low) proxy, the synthetic ⁴⁰Ar release-rate peaks are 550, 500,
and 500 °C for the `4×4×6`, `8×8×6`, and `12×12×6` volumes. With the retained
63 kcal mol⁻¹ (delamination-high) proxy, all three peak at 600 °C.

The result rules out unconditional basal accessibility as the sole cause of E4's
failure. It also shows that a gate driven independently at every basal cell does
not create a length scale across the existing comparison ladder: the ladder
changes the in-plane dimensions while retaining six layers along the open basal
axis. The next discriminating mechanism is a **surface-connected lateral
release front**—edge accessibility that advances through connected delamination
in the `a/b` plane, without a grain-size-specific parameter. Follow-up card
`E4a2-surface-connected-lateral-release` owns that test.

## Structural gate

`Surface_gate` now has three explicit states:

- `inert` for non-surface cells;
- `closed` for the two basal faces before local delamination; and
- `open` only after the gate's bonded local `Interface` reaches `delaminated`.

The former direct `release_*_delamination_*` bypass is removed. All six isotope ×
gallery release rules require a neighboring `Surface_gate.open`. The retained
58 ± 5 kcal mol⁻¹ interface barrier remains labelled as a proxy; E3a's
64.095991 and 68.410690 kcal mol⁻¹ Ar-hop sensitivities remain labelled
`incomplete-convergence`. No per-volume, per-grain-size, or per-temperature
multiplier was introduced.

The fail-closed regression
`test_surface_release_requires_local_delamination_accessibility` verifies that
bulk gates remain inert, basal gates initialize closed, local delamination is
the only opening path, and every ⁴⁰Ar/³⁹Ar/³⁶Ar release rule requires the open
gate.

## Campaign and replay evidence

The exact E4 Cartesian campaign was rerun: two retained delamination
sensitivities × three E2b volumes. Each primary ensemble used eight ordered seeds
(`19980`–`19987`) with `--paranoid`; each deck also ran two independent
two-replica same-seed replays.

| Evidence | Value |
|---|---:|
| Primary ensemble replicas | 48 |
| Replay replicas | 24 |
| Total executed replicas | 72 |
| Primary KMC events | 98,264 |
| Primary elapsed time | 1.368 s |
| Replay checks | 6 / 6 byte-identical |

Byte identity holds for `ensemble.csv`, `ensemble-summary.csv`, and
`observables.csv` for every replay pair. Raw logs and trajectories are retained
at `/private/tmp/e4a-campaign-20260913-1458`; compact receipts and comparison
products are committed under
`docs/program/results/E4a-surface-gated-release/`.

## Grain-size verdict

| Delamination sensitivity | `4×4×6` peak | `8×8×6` peak | `12×12×6` peak | 700→800→1025 °C trend |
|---|---:|---:|---:|---|
| 53 kcal mol⁻¹ proxy | 550 °C | 500 °C | 500 °C | **NO-GO** — unordered |
| 63 kcal mol⁻¹ proxy | 600 °C | 600 °C | 600 °C | **NO-GO** — no size ordering |

Final cumulative ⁴⁰Ar release is 0.359/0.548/0.544 for the low-sensitivity
volume ladder and 0.386/0.500/0.503 for the high-sensitivity ladder. The
previous high-temperature full-release merge is therefore not preserved by the
local gate within the tracked schedule. That is reported as a model result, not
tuned away.

## Products

- `petra/scripts/build_muscovite_full_deck.py` and the tracked full deck;
- six regenerated `petra/decks/muscovite-1998/*.toml` campaign decks;
- `petra/scripts/muscovite_1998_comparison.py`, with data-driven crossover verdicts;
- `synthetic-spectra.csv`, `comparison.json`, `comparison.svg`, and
  `campaign-receipts.json` in this result directory.

## Reproduction

```bash
python3.14 petra/scripts/build_muscovite_full_deck.py \
  --output petra/decks/muscovite-full-mechanism.toml
python3.14 petra/scripts/muscovite_1998_comparison.py generate \
  petra/decks/muscovite-1998
python3.14 petra/scripts/muscovite_1998_comparison.py run \
  . "$CARGO_TARGET_DIR/release/petra" petra/decks/muscovite-1998 \
  /private/tmp/e4a-campaign-20260913-1458 --replicas 8 --base-seed 19980
python3.14 petra/scripts/muscovite_1998_comparison.py analyze \
  /private/tmp/e4a-campaign-20260913-1458 petra/decks/muscovite-1998 \
  petra/data/muscovite-1998 docs/program/results/E4a-surface-gated-release
```
