# A8b-legacy-dft-transition-state-gap — connect the 1999 reaction ledger to modern barrier work

- status: blocked
- track: A (geochemistry)
- priority: P2
- machine: workstation
- depends: A3-barrier-ladder
- claimed-by: —

## Objective

A8 proves the 1999 archive contains total reaction energies, not transition-state calculations. Convert that negative finding into a bounded modern campaign rather than mislabeling negative-frequency structures as legacy barriers.

1. Rank the 45 Chapter 4 reaction-energy mappings in `legacy/thesis-archive/DFT-LEDGER.csv` by relevance to the current dissolution mechanism program.
2. Map overlapping reactant/product chemistry onto A3's production protocol and existing stationary-point bank.
3. Select a small first tranche with defensible endpoint atom maps; locate TSs and verify each with exactly one intended imaginary mode plus bilateral endpoint evidence.
4. Add only verified reactant/TS/product triplets to a new modern ledger layer. Never mutate A8's archival evidence rows.

## Acceptance

- Explicit crosswalk from each selected legacy reaction ID to modern structures/method tier.
- At least one independently verified barrier or a well-evidenced no-go verdict.
- No archival `unmatched` row promoted solely because its frequency job contains a negative mode.
- STATUS.md and this card updated in the delivering PR.

## Progress

- 2026-09-02 — card filed from A8 after exhaustive parsing found zero defensible TS logs among 83 Gaussian result logs. (hermes-custom-build-001; profile=workstation)
