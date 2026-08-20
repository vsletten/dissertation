# E1 muscovite phase 1 — published-number deck and isothermal release

**Card:** `E1-muscovite-deck`

**Scope:** qualitative mechanism gate, not a fit to the 1998 observations

**Engine contract:** Petra schema v1 (fixed-temperature CTMC)

## Verdict

**PASS at both 500 °C and 700 °C.** The phase-1 deck produces an early
release burst followed by a long structural tail. When each synthetic curve
is deliberately inverted through the Hames–Bowring cylindrical muscovite
geometry used by Sletten & Onstott (1998), the apparent `D/a²` rises early
and then falls by more than an order of magnitude:

| Run | Released | Early `D/a²` | Peak `D/a²` | Peak fraction | Tail median `D/a²` | Rise | Fall |
|---|---:|---:|---:|---:|---:|---:|---:|
| 500 °C | 630 / 831 (75.8%) | 7.104×10⁻⁷ s⁻¹ | 1.515×10⁻⁶ s⁻¹ | 4% | 5.501×10⁻⁹ s⁻¹ | 2.13× | 275× |
| 700 °C | 631 / 831 (75.9%) | 1.394×10⁻⁴ s⁻¹ | 1.889×10⁻³ s⁻¹ | 12% | 5.842×10⁻⁵ s⁻¹ | 13.55× | 32.3× |

The absolute clock accelerates sharply at 700 °C, but the same two-stage
shape survives. That is the phase-1 claim from
[`docs/scoping/ar-muscovite.md`](../../scoping/ar-muscovite.md) §5 claim 1:
decomposition-gated transport can look diffusion-like briefly while its
inferred diffusivity evolves strongly through the run.

## Artifacts

- Decks:
  - [`petra/decks/muscovite-phase1-500c.toml`](../../../petra/decks/muscovite-phase1-500c.toml)
  - [`petra/decks/muscovite-phase1-700c.toml`](../../../petra/decks/muscovite-phase1-700c.toml)
- Post-processor and qualitative gate:
  [`petra/scripts/muscovite_release.py`](../../../petra/scripts/muscovite_release.py)
- Machine-readable verdict:
  [`E1-muscovite-phase1/summary.json`](E1-muscovite-phase1/summary.json)
- 2%-release data:
  [`500-c.csv`](E1-muscovite-phase1/500-c.csv) and
  [`700-c.csv`](E1-muscovite-phase1/700-c.csv)
- Plots:
  - [cumulative release vs √time](E1-muscovite-phase1/release-curves.svg)
  - [apparent `D/a²` vs released fraction](E1-muscovite-phase1/apparent-da2.svg)

The SVGs were rasterized and visually inspected after generation; labels,
axes, the early rise, and the later tail are all visible without clipping.
The raw replay logs are intentionally not committed: each is reproducible
from the decks, while the compact CSV/SVG/JSON products are the durable
result.

## Model implemented

The deck is a deliberately small single-gallery reduction, not an atomistic
muscovite structure:

- `32 × 32` independent c-axis columns, eight 2M₁ repeats deep;
  32,768 total Petra sites.
- Alternating `K_gallery_A/B` sites form a symmetric one-dimensional gallery
  chain. Separate intra/inter edge labels give four directionally symmetric
  hop reactions and avoid Petra v1's first-matching-neighbor bias.
- The seeded gallery contains 94.5% vacancies, 5% ⁴⁰Ar, and 0.5% immobile K
  in expectation. Seed 1998 produced 831 Ar and 79 K sites in both runs.
- A dehydroxylation actor is an explicitly paired `OH_donor` +
  `OH_acceptor`: one event removes the donor as water and leaves one
  `O_residual`.
- A gallery hop requires (a) two adjacent vacancies and (b) local
  `O_residual`. Ar therefore cannot migrate through a repeat before that
  repeat dehydroxylates.
- The c axis is open. Surface Ar desorbs with a 1 s⁻¹ lumped rate, which is
  fast relative to migration at both temperatures without introducing an
  unnecessary 10¹²-to-10⁻⁸ rate span into the Fenwick total.

## Number provenance and assumptions

The phase-1 rule was **no fitted barrier campaign**. Barrier-valued inputs are
therefore published anchors already tabled in the scoping document:

1. **Ar divacancy hop:** 66 kcal/mol, Nteme et al. (2022).
2. **Bulk dehydroxylation:** 57 kcal/mol, Kodama & Brydon (1968), the value
   used by Sletten & Onstott (1998).
3. **Steric slowdown:** `by_count = [0, 3, 12]` kcal/mol on top of 57,
   producing 57/60/69 kcal/mol for zero/one/two residual-O neighbors. These
   are a phase-1 bracketing table spanning the published bulk-to-modern
   in-situ envelope (251–290 kJ/mol), **not** a claim that those three
   environment-resolved barriers were measured. Computing that table remains
   phase 3.

Only activation enthalpies are available, so both activated processes use
Petra's Eyring expression with `ΔS‡ = 0`. That is an explicit prefactor proxy,
not a measured activation entropy. The surface rate and initial population
fractions are modeling choices, not literature measurements. No unlabelled
calibration was performed against the target curves.

## Fickian inversion and mechanical gate

For every `release_*` event, the post-processor records cumulative fraction
released. At each 2% crossing it inverts the infinite-cylinder series from
Crank (1975), using the radial geometry of
[Hames & Bowring (1994)](https://doi.org/10.1016/0012-821X(94)00079-4),
cited by the 1998 paper:

`F(τ) = 1 − 4 Σ exp(−αₙ²τ)/αₙ²`, where `J₀(αₙ) = 0` and `τ = Dt/a²`.

The first 64 positive `J₀` roots evaluate the series for `τ ≥ 10⁻³`; the
cylindrical short-time expansion `F = 4√(τ/π) − τ` handles smaller `τ`.
Tests round-trip the inversion from 1% through 99% release.

The post-processor then reports the interval slope `Δτ/Δt` as apparent
`D/a²`. This intentionally forces the structurally controlled trajectory
through the continuum lens used to interpret the 1998 release experiments.

The qualitative gate is fixed in code:

- at least 100 initial Ar atoms and at least 50% release;
- an early peak by 12% release, later than the first 2% point, at least 1.5×
  the first inferred `D/a²`;
- the early peak is at least 10× the median of the final five 2% intervals.

Both runs pass. Reaction-direction counts are also balanced (all four hop
classes are within about 5%), so the shape is not a one-way topology artifact.

## Reproduction

From the repository root:

```bash
cd petra
export OMP_NUM_THREADS=16

set -o pipefail
nice -n 10 cargo run --release -p petra-cli -- \
  decks/muscovite-phase1-500c.toml \
  --out /tmp/muscovite-phase1-500c --viz --paranoid \
  2>&1 | tee /tmp/muscovite-phase1-500c.log

nice -n 10 cargo run --release -p petra-cli -- \
  decks/muscovite-phase1-700c.toml \
  --out /tmp/muscovite-phase1-700c --viz --paranoid \
  2>&1 | tee /tmp/muscovite-phase1-700c.log

cd ..
python3 petra/scripts/muscovite_release.py \
  --run '500 C:500:/tmp/muscovite-phase1-500c' \
  --run '700 C:700:/tmp/muscovite-phase1-700c' \
  --out-dir docs/program/results/E1-muscovite-phase1
```

Focused post-processing tests:

```bash
cd petra
python3 -m unittest discover -s scripts/tests -p 'test_*.py' -v
```

## Honest gap list for phase 2

1. **Qualitative, not fitted.** The 1998 curves/tables have not yet been
   digitized into a residual or likelihood; this phase proves shape only.
2. **One-dimensional gallery.** There is no in-plane percolation, K-vacancy
   migration, or explicit divacancy-complex state. About 24% of Ar remains
   trapped when no legal event remains; that residual is not asserted to be
   an experimental retained fraction.
3. **Steric table is bracketed, not computed.** Phase 3 must replace the
   published-envelope mapping with environment-resolved kinetics and
   sensitivity bands.
4. **No delamination or extended defects.** Those are deliberately phase 2,
   along with ³⁶Ar, ³⁹Ar, octahedral recoil traps, and Xe.
5. **No staircase schedule.** Petra schema v1 is fixed-temperature; B5 must
   land piecewise-isothermal execution before synthetic age spectra.
6. **No grain-scale claim.** This is a representative gallery stack, not a
   52–64 µm flake or 3000 µm disk. The Hames–Bowring cylinder inversion is a
   deliberately imposed laboratory interpretation, not a statement that the
   lattice dimensions equal an experimental radial diffusion length.
7. **No activation entropy or measured surface prefactor.** `ΔS‡ = 0` and
   1 s⁻¹ surface desorption are explicit phase-1 proxies.

Phase 1 does what it was supposed to do and no more: published barriers plus
structural gating are sufficient to generate the target non-Fickian
rise-then-fall signature. The quantitative fight belongs to phase 2/3.
