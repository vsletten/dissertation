# D3-ice-mantle-deck — grain-surface chemistry KMC

- status: done
- track: D (astrochemistry)
- priority: P2
- machine: any
- depends: B3, D2a
- claimed-by: hermes-laptop

## Objective
Ice-mantle deck per docs/scoping/astrochemistry.md: adsorption
(deposition/source events), thermal + tunneling diffusion,
Langmuir-Hinshelwood reaction, desorption; quenched per-site
binding-energy disorder. First target: H2 formation efficiency vs
temperature curve; then CO hydrogenation with D2a's rates. Headline
observable: stochastic small-grain effects (where rate equations fail).
Needs the RFC's non-Arrhenius rate laws + source events (flag any gap
back to a B-track card).

## Acceptance
- Deck runs; H2-efficiency-vs-T reproduces the literature shape;
  comparison vs published KMC figures (Simons 2020) documented.

## Progress
- 2026-08-22 — hermes-laptop — claimed atomically as
  `agents/D3-ice-mantle-deck`; added the deck and strict RED→GREEN acceptance
  gates for its process inventory, seeded disorder, and finite efficiency
  window.
- 2026-08-22 — hermes-laptop — generated a bounded 16-replica temperature
  sweep. Efficiency rises from 0.310 at 6 K to 0.962 at 12 K, stays >0.9
  through 16 K, then falls to 0.559/0.148/0.055 at 20/24/26 K.
- 2026-08-22 — hermes-laptop — DEVIATION: did not build CO hydrogenation from
  D2a's rates because D2a's merged verdict explicitly says those gas-phase
  rates are 4–5 orders low for surface LH chemistry. Created
  `D2b-explicit-surface-rates` and blocked `D3b-co-hydrogenation-deck` on that
  gate rather than laundering invalid rates into a production deck.
- 2026-08-18 — card created (fable).

## Result

- Shipped `petra/decks/ice-mantle-h2.toml`: 24×24 periodic grain, open H bath,
  thermal + tunneling hopping, LH H2 formation, H/H2 desorption, and a seeded
  20% deep-site quenched class.
- Shipped deterministic sweep tool
  `petra/crates/petra-deck/examples/h2_efficiency.rs`, regression suite
  `petra/crates/petra-deck/tests/ice_mantle.rs`, and generated data
  `docs/program/results/D3-ice-mantle-deck/efficiency.csv`.
- Comparison, model limits, and honest Simons-2020 disposition:
  `docs/program/results/D3-ice-mantle-deck.md`.
- Verification: targeted acceptance tests, full Petra workspace tests, format,
  Clippy, and exact generated-product replay are required green in this PR.
