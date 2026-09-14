# A9b reservoir/origin contract

## Bottom line

A9's `1e-30` Al/Si activities remain preserved only in the historical
`kaolinite-approx.toml` campaign deck. New work materializes one of three
explicit pH 3–5 open-flow profiles through
`petra/scripts/reservoir_origin_contract.py`; generated decks use a declared
`1e-6` trace-product activity rather than presenting a numerical underflow sink
as realistic chemistry.

The contract is executable and origin-safe for the physical dissolution
numerator. Initial occupied Si/Al sites are `original_lattice`; adsorption marks
the cation `reservoir`; intervening surface-state transitions preserve that
lineage; desorption increments physical dissolution only for
`original_lattice`. Reservoir desorption remains visible as gross exchange but
never enters the lattice-release numerator.

## Provenance table

The P&K ledger supplies the project’s far-from-equilibrium kaolinite
acid/neutral/base law and H+ reaction order.[1] Cama, Metz & Ganor provide the
acidic pH/temperature laboratory anchor.[2] Yang & Steefel provide the directly
relevant continuously stirred flow-through pH 4 experiment near ambient
temperature.[3]

| Profile | T (K) | a(H+) | Al reservoir species | a(Al) | Si reservoir species | a(Si) | Boundary |
|---|---:|---:|---|---:|---|---:|---|
| pH 3 | 298.0 | `1e-3` | `Al3+` (total-Al proxy) | `1e-6` | `H4SiO4(aq)` (total-Si proxy) | `1e-6` | constant-activity open flow |
| pH 4 | 298.0 | `1e-4` | `Al3+` (total-Al proxy) | `1e-6` | `H4SiO4(aq)` (total-Si proxy) | `1e-6` | constant-activity open flow |
| pH 5 | 298.0 | `1e-5` | `Al3+` (total-Al proxy) | `1e-6` | `H4SiO4(aq)` (total-Si proxy) | `1e-6` | constant-activity open flow |

`a(H+) = 10^(-pH)` is exact by contract. The product activity is a declared
1-micromolal, unit-activity-coefficient trace-inflow condition. It is not a
measured universal kaolinite composition, not a solubility result, and not a
claim that all three profiles have the same saturation index.

## Boundary semantics and approximations

- The reservoir is continuously refreshed: inflow activities remain constant,
  and released products do not feed back into the declared bath.
- Activity coefficients are fixed at one, so numerical activities are dilute
  molality proxies.
- Total Al is represented as `Al3+`; hydrolyzed and polynuclear Al are excluded.
- Total Si is represented as neutral `H4SiO4(aq)`; deprotonated and polymerized
  silica are excluded over this bounded pH 3–5 contract.
- H+ activity is emitted explicitly. The P&K acid order `n_H = 0.777` scales
  both terminal Si/Al release barriers relative to the pH 4 anchor. This is a
  bounded survey surrogate for the macroscopic acid-rate law, not a claim that
  Petra now has proton-explicit elementary chemistry.
- The historical `mu(Al) = mu(Si) = -1 kcal/mol` offset is retained to avoid
  silently changing two boundary knobs at once. The selected activity and mu
  therefore remain a survey sensitivity condition, not calibrated chemical
  potentials.

The authoritative machine-readable form is
`petra/examples/kaolinite-reservoirs.toml`. It lists all sources and every
approximation; the runner fails closed unless the three profiles cover exactly
pH 3, 4, and 5, the H+ arithmetic is exact, Al/Si use the declared trace value,
and the old `1e-30` sink is absent from generated decks.

## Executable evidence

Generated decks and hash-bound receipts:

| pH | Deck | SHA-256 | Evidence receipt SHA-256 |
|---:|---|---|---|
| 3 | `a9b-reservoir-origin-contract/kaolinite-ph3.toml` | `17e40e0d3bc3b388bc126a81f5385c81b953ea69b97b6e78b1295130c2269f9b` | `71bff310f056bfedddc1b57f9350299bfb747c10e9626553256a56c856fe87a1` |
| 4 | `a9b-reservoir-origin-contract/kaolinite-ph4.toml` | `9316eb842ebaf96c8f8cbf45796dff163a13cb9e6243dc67c8f5e3b84152767c` | `4cd16e348f84367139a357a0f3d4fadbb3f5a04616e80dfd8f0693396fc54b00` |
| 5 | `a9b-reservoir-origin-contract/kaolinite-ph5.toml` | `b4f928877bf14ef5840cea194ba6575c93310c8a7d7b716e4968d7ca08017ac9` | `cc9b470e7aeebb1d48fbea7879af075f9df2cae060e2d3dab5e2380f2a349174` |

The production A9 event-replay regression performs the deterministic
adversarial fixture for both Si and Al: true lattice desorption increments the
numerator; adsorption → repeated live surface-state transitions → desorption
does not. Contract tests additionally verify that all generated decks pass the
existing campaign validator and that their terminal release barriers order pH
3 faster than pH 4 faster than pH 5.

## Limits and next consumer

This card does not rerun the 29 × 8 sensitivity ensemble and does not publish a
calibrated rate. `A9b-sensitivity-ranking` may consume these generated profiles
only after the independent Verify stage has checked the source arithmetic and
adversarial origin semantics together with the separate reachability result.

## Sources

[1] https://pubs.usgs.gov/of/2004/1068/pdf/OFR_2004_1068.pdf — Palandri and Kharaka 2004
[2] https://doi.org/10.1016/S0016-7037%2802%2900966-3 — Cama, Metz and Ganor 2002
[3] https://doi.org/10.1016/j.gca.2007.10.011 — Yang and Steefel 2008
