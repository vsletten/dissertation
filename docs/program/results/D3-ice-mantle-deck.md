# D3 ice-mantle deck — H₂ formation benchmark

## Verdict

**PASS for the card's H₂ acceptance gate.** Petra now ships a microscopic,
open-system grain deck with deposition, thermal and tunneling diffusion,
Langmuir–Hinshelwood recombination, H/H₂ desorption, and seeded quenched
binding-site disorder. A 16-replica temperature sweep produces the expected
finite efficiency window: mobility limits the cold edge, efficiency is near
unity from 10–16 K, and desorption suppresses formation above 20 K.

This is a **shape/conformance benchmark**, not an absolute astrochemical rate
prediction. The deck rate-compresses the activated-process scale to avoid the
many-hop stiffness that Petra cannot yet accelerate; it preserves the ordering
of deposition, diffusion, reaction, and desorption but must not be used to
publish physical H₂ rates.

## Artifacts

- Deck: `petra/decks/ice-mantle-h2.toml`
- Deterministic sweep: `petra/crates/petra-deck/examples/h2_efficiency.rs`
- Acceptance regressions: `petra/crates/petra-deck/tests/ice_mantle.rs`
- Generated sweep: `docs/program/results/D3-ice-mantle-deck/efficiency.csv`

Reproduce from `petra/` (with the repository's canonical Cargo target):

```bash
cargo test -p petra-deck --test ice_mantle -- --nocapture
cargo run -p petra-deck --example h2_efficiency -- \
  decks/ice-mantle-h2.toml \
  ../docs/program/results/D3-ice-mantle-deck/efficiency.csv
```

The sweep is bounded at 16 deterministic replicas × 300 arriving H atoms per
temperature, with a 250,000-event fail-closed cap per replica. Efficiency is
`2 × formed H₂ / arriving H`.

## Temperature sweep

| T (K) | Mean efficiency | 95% interval |
|---:|---:|---:|
| 6 | 0.310 | 0.291–0.330 |
| 8 | 0.859 | 0.850–0.868 |
| 10 | 0.941 | 0.935–0.947 |
| 12 | 0.962 | 0.957–0.966 |
| 14 | 0.960 | 0.955–0.965 |
| 16 | 0.916 | 0.910–0.922 |
| 18 | 0.806 | 0.792–0.819 |
| 20 | 0.559 | 0.535–0.583 |
| 22 | 0.308 | 0.291–0.325 |
| 24 | 0.148 | 0.134–0.163 |
| 26 | 0.055 | 0.043–0.066 |

The regression gate independently checks 6/12/24 K and observed efficiencies
0.273/0.960/0.127 for seed 7001.

## Comparison with published microscopic KMC

Cuppen, Morata, and Herbst describe continuous-time random-walk simulations in
which individual H atoms deposit, hop, evaporate, and form H₂; their olivine
parameters use multiple local binding classes, and they report that disorder /
roughness broadens efficient H₂ formation from roughly 6–9 K on a smooth
surface to about 20 K on a very rough surface.[2] Petra reproduces that
**qualitative topology**: a narrow cold mobility edge, a broad efficient middle
window, and a high-temperature desorption edge; its seeded 20% deep-site class
is the deck-level analogue of quenched local binding heterogeneity.

Simons, Lamberts, and Cuppen use a 50×50 periodic microscopic KMC surface with
the same four process classes—deposition, reaction, hopping, and desorption—and
show CO-hydrogenation mantle evolution from 8 to 20 K (their Figures 6–7), with
temperature effects becoming pronounced above 16 K.[3] The Petra deck is 24×24
periodic and its efficiency likewise turns downward after 16 K. This is an
architectural/process and temperature-trend comparison only: Simons 2020 does
**not** publish an H₂-efficiency-vs-temperature curve, so claiming numerical
agreement to those figures would be false.

## Model details and limitations

- **Quenched disorder:** 20% of sites are deterministically re-kinded to a deep
  binding class by the init RNG. Same seed gives the identical disorder map;
  another seed changes it. The class never mutates during the trajectory.
- **Open bath:** each vacant site carries an independent constant H deposition
  event. This exercises Petra's source-event semantics directly.
- **Diffusion:** thermal Arrhenius channels and an additive constant tunneling
  floor are separate events. A landed H becomes `H_mobile` after a hop;
  recombination requires that encounter state, preventing adjacent random
  deposition from masquerading as Langmuir–Hinshelwood diffusion.
- **Recombination:** nearest-neighbor H + H produces one adsorbed H₂ and one
  vacancy; H₂ subsequently desorbs.
- **Rate compression:** physical trial frequencies and diffuse-cloud fluxes
  create an extreme separation of timescales (many hops per arrival). This deck
  compresses that ratio for bounded CI. A physical-scale rerun waits on Petra's
  acceleration/superbasin work.
- **CO hydrogenation:** not included. D2a's merged verdict is NO-GO for using its
  cheap gas-phase rates as surface Langmuir–Hinshelwood rates (4–5 orders of
  magnitude low). `D3b-co-hydrogenation-deck` records the explicit-surface rate
  dependency rather than laundering invalid numbers into this deck.

## Sources

[2] https://arxiv.org/pdf/astro-ph/0601554 — Cuppen et al. (2006) H2 formation on stochastically heated grains
[3] https://www.aanda.org/articles/aa/full_html/2020/02/aa36522-19/aa36522-19.html — Formation of COMs through CO hydrogenation
