# E2b grain-size sweep — ensemble statistics for the muscovite spectra

**Card:** `E2b-grain-size-sweep`
**Mechanism tier:** E2 proxy-barrier mechanism; this is a finite-volume/statistics study, not a kinetic calibration

## Verdict

**The E2 staircase is a stable mechanism-level ensemble feature from 49,152 sites upward.** Seven lattice volumes from `4×4×6` through `36×36×54` were each run as 32 deterministic replicas. The stabilization gate was declared in the analysis code as:

- maximum adjacent-size change in any mean cumulative isotope-release fraction ≤ 0.05; and
- maximum adjacent-size relative change in every defined mean apparent-age step ≤ 0.10; and
- maximum adjacent-size change in the fraction of replicas with a defined apparent age ≤ 0.05 for each material step (a step where either size releases at least 0.1% of its initial ³⁹Ar inventory).

The `20,736 → 49,152`, `49,152 → 165,888`, and `165,888 → 559,872` transitions all pass. At the largest transition the maxima are 0.01634 in released fraction, 0.01008 in relative step age, and 0.0 in defined-age fraction. The first stable volume is therefore **49,152 Petra sites (`16×16×24` cells)**, and that verdict survives two further volume increases totaling an 11.39× site-count expansion.

The qualitative spectrum persists rather than washing out: an old first step, a lower second step, a near-zero third step, then zero-age tail steps dominated by delayed ³⁹Ar release. At the largest completed volume the bootstrap means are 38.28 Ma, 21.31 Ma, 0.402 Ma, 0 Ma, and 0 Ma. These are synthetic `J = 0.01` observables from E2's labelled proxy mechanism, not geological age estimates.

## Size ladder and practical ceiling

Every row is 32 replicas under the same five-step 500–700 °C schedule. `events/s` includes the complete ensemble wall time on the laptop host. Each size also passed an independent two-replica exact-seed replay: `ensemble.csv`, `ensemble-summary.csv`, and `observables.csv` were byte-identical between repeats.

| cells | Petra sites | replicas | events | wall s | events/s | replay |
|---:|---:|---:|---:|---:|---:|:---:|
| 4×4×6 | 768 | 32 | 9,174 | 0.366 | 25,088 | exact |
| 6×6×9 | 2,592 | 32 | 30,865 | 0.367 | 84,144 | exact |
| 8×8×12 | 6,144 | 32 | 73,536 | 0.468 | 157,063 | exact |
| 12×12×18 | 20,736 | 32 | 248,955 | 0.867 | 287,054 | exact |
| 16×16×24 | 49,152 | 32 | 591,474 | 1.730 | 341,965 | exact |
| 24×24×36 | 165,888 | 32 | 1,995,262 | 7.206 | 276,874 | exact |
| 36×36×54 | 559,872 | 32 | 6,740,495 | 138.420 | 48,696 | exact |

The **largest comfortably completed ensemble is 559,872 sites**. The next rung, `48×48×72` (1,327,104 sites × 32 replicas), was terminated by `SIGKILL` before producing output; its log is empty. The run therefore treats 559,872 as this host/configuration's evidence-backed practical ceiling rather than extrapolating through the failed rung. The sharp throughput drop at 559,872 and the failed next rung make that boundary operationally meaningful.

## Ensemble spectrum and release bands

The stabilized volumes agree on the staircase and the release partitioning. Values below are means with deterministic 2,000-resample bootstrap 95% intervals for apparent age; release columns are mean cumulative inventory fractions. A ratio or age is undefined when its denominator is zero, so ratio/age bands are conditional on the replicas where the metric exists rather than padded with zeros. `spectrum-bands.csv` records both the defined and missing replica counts beside every distribution and confidence band.

### Largest completed volume: 36×36×54 cells (559,872 sites)

| °C step | apparent age, Ma (95% CI; defined n/32) | cumulative ⁴⁰Ar | cumulative ³⁹Ar | cumulative ³⁶Ar | step ³⁶Ar/⁴⁰Ar |
|---:|---:|---:|---:|---:|---:|
| 500 | 38.279 (37.729–38.836; 32/32) | 0.1751 | 0.0225 | 0.1714 | 0.1977 |
| 550 | 21.315 (21.168–21.453; 32/32) | 0.9586 | 0.2044 | 0.9584 | 0.2027 |
| 600 | 0.402 (0.394–0.411; 32/32) | 1.0000 | 0.7165 | 1.0000 | 0.2028 |
| 650 | 0.000 (0.000–0.000; 32/32) | 1.0000 | 0.9920 | 1.0000 | — |
| 700 | 0.000 (0.000–0.000; 25/32) | 1.0000 | 0.9921 | 1.0000 | — |

At 500 °C only 2.25% of ³⁹Ar has escaped while 17.51% of ⁴⁰Ar and 17.14% of ³⁶Ar have escaped. At 550 °C the two gallery species are ≈95.8% released while ³⁹Ar remains only 20.44% released. The delayed ³⁹Ar then dominates the 600–650 °C steps. That reservoir separation is the mechanism producing the old-first-step staircase; increasing volume narrows its uncertainty rather than removing it.

The E2 parameterization does **not** generate a large evolving ³⁶Ar/⁴⁰Ar step ratio: the largest-volume means are 0.1977, 0.2027, and 0.2028 for the three steps with ⁴⁰Ar release. This is an honest negative for scoping claim 2 at the current mechanism tier, while the staircase/³⁹Ar-retardation claim is positive. E3a owns replacement of the proxy barriers; no fitted interpretation is introduced here.

## Stabilization details

| previous sites → current sites | max release-fraction Δ | max relative age Δ | max defined-age-fraction Δ¹ | verdict |
|---:|---:|---:|---:|:---:|
| 768 → 2,592 | 0.10040 | 0.20388 | 0.09375 | changing |
| 2,592 → 6,144 | 0.04707 | 0.22512 | 0.00000 | changing |
| 6,144 → 20,736 | 0.08353 | 0.09594 | 0.00000 | changing |
| 20,736 → 49,152 | 0.02155 | 0.03419 | 0.00000 | stable |
| 49,152 → 165,888 | 0.02762 | 0.02724 | 0.00000 | stable |
| 165,888 → 559,872 | 0.01634 | 0.01008 | 0.00000 | stable |

¹ Evaluated only on material steps (either adjacent size releases ≥0.1% of initial ³⁹Ar). This excludes zero-denominator tail noise while preventing a materially changing set of contributing replicas from passing on similar conditional means alone.

## Artifacts and reproduction

Tracked products:

- `docs/program/results/E2b-grain-size-sweep/spectrum-bands.csv` — per-step means, bootstrap CIs, contributing values, and explicit defined/missing counts out of 32 replicas.
- `docs/program/results/E2b-grain-size-sweep/size-scaling.csv` — measured volumes, events, throughput, and replay receipts.
- `docs/program/results/E2b-grain-size-sweep/stability.json` — machine-readable adjacent-size gates and the 49,152-site stabilization verdict.

Campaign receipts and long-run logs are intentionally outside Git at:

- `/opt/data/logs/e2b-grain-size-sweep/campaign-20260828-v2/`
- failed-ceiling receipt: `logs/size-48x48x72.log` (empty after `SIGKILL`)

Reproduce with the repository's mechanically obtained canonical Cargo target:

```bash
export CARGO_TARGET_DIR="$(python3 /opt/data/workspaces/vsletten/mission-control/main/scripts/review_ready_teardown.py --repo vsletten/dissertation --print-cargo-target)"
export OMP_NUM_THREADS=16 MKL_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 RAYON_NUM_THREADS=16
cargo build --release -p petra-cli
python3 petra/scripts/muscovite_grain_size_sweep.py run \
  petra "$CARGO_TARGET_DIR/release/petra" /path/to/new/campaign
python3 petra/scripts/muscovite_grain_size_sweep.py analyze \
  /path/to/new/campaign docs/program/results/E2b-grain-size-sweep
```

## Verification

- Seven 32-replica scheduled ensembles completed through 559,872 sites.
- Same-seed replay is byte-identical at every completed size.
- The complete sweep executed 9,689,761 events across the primary ensembles.
- Focused Python tests cover deterministic size-specific deck rendering, ensemble isotope accounting (including zero-inventory sparse replicas), confidence/distribution products, stabilization missingness, ladder validation, event receipts, and incremental campaign receipts.
- Petra workspace tests, Rust lint, Python tests, and Python lint/format gates are recorded in the card closeout.
