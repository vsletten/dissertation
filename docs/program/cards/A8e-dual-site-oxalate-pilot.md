# A8e-dual-site-oxalate-pilot — test the adjacent-edge promotion mechanism

- status: blocked
- track: A (geochemistry / ligand-promotion)
- priority: P2
- machine: workstation
- depends: A3-barrier-ladder, A8d-ligand-promotion-track-scoping
- claimed-by: —

## Objective

Execute only Gate 1 from the A8d scoping memo: a six-system, two-pH-state pilot that tests whether paired oxalate occupancy on neighboring kaolinite-like edge Al sites lowers the connected Al-O-Si cleavage/Al-detachment barrier relative to both unliganded and singly occupied controls.

This card is not authorization for the malonate/succinate extension. That extension requires every promotion gate below to pass.

## Inputs

- [`docs/scoping/ligand-promoted-aluminosilicate-dissolution.md`](../../scoping/ligand-promoted-aluminosilicate-dissolution.md)
- A3 model/protonation/TS/IRC infrastructure and receipts
- ranked A8d legacy starting structures in `legacy/thesis-archive/DFT-LEDGER.csv`

## Required matrix

For each of two documented pH/speciation conventions (acidic near pH 3; near-neutral near pH 6-6.5):

1. hydrated adjacent-Al edge with no ligand;
2. the matched edge with one oxalate;
3. the matched edge with two oxalates on neighboring Al sites.

Before generating structures, document aqueous oxalate/bioxalate, relevant Al-oxalate, edge protonation, cluster charge, counterion, and explicit-water choices. Compare thermodynamic cycles with identical atom/proton counts.

## Acceptance

- Six-system matrix is complete or each bounded failure has a receipt.
- Every claimed barrier is a first-order saddle with one relevant imaginary mode and full IRC connectivity to the intended endpoint families.
- Adsorption `ΔG`, cleavage/detachment `ΔG‡`, reaction `ΔG`, key coordination distances, proton-transfer topology, and uncertainty/sensitivity envelopes are tabulated.
- Matched no-ligand and single-oxalate controls use identical edge, water-shell, charge, and method conventions.
- A larger explicit-water-shell check and the A2 production single-point tier are applied to the unliganded/dual-oxalate comparison.
- Promotion decision is explicit against all A8d gates: at least 10 kJ/mol robust barrier lowering, paired occupancy distinguishable from single occupancy, same-sign water/method sensitivity, and closed thermodynamic bookkeeping.
- If any gate fails, record NO-GO and do not launch homologous ligand work.
- If every gate passes, file a separate malonate/succinate extension card before doing that work.
- Compute receipts, artifact hashes, and the final memo are committed; `STATUS.md` and this card are updated.

## Bounds

- Planning envelope: six reaction systems; benchmark before bulk launch.
- All execution loops obey the bounded-worker and GPU-lease policy.
- No periodic slab escalation without a separate card justified by a sign-changing edge-boundary check.
- No broad ligand screen, ERW reactor model, or production campaign belongs in this card.

## Progress

- 2026-09-03 — filed from A8d's CONDITIONAL GO decision; blocked on A3 barrier infrastructure and deliberately limited to the adjacent-site oxalate falsification pilot. (hermes-macbot-zero; profile=laptop)
