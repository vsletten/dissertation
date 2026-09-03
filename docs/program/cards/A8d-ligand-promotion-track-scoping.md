# A8d-ligand-promotion-track-scoping — scope oxalate/malonate/succinate dissolution work

- status: done
- track: A (geochemistry / ligand-promotion candidate)
- priority: P2
- machine: any
- depends: A8-thesis-archive-intake
- claimed-by: hermes-macbot-zero

## Objective

Decide whether the recovered 1999 ligand calculations justify a dedicated modern ligand-promoted dissolution track, with explicit relevance to enhanced rock weathering rather than nostalgia-driven compute.

1. Read the Chapter 4/5 mechanism and bonding claims together with the 45-row archival reaction-energy ledger.
2. Survey modern experimental and computational evidence for oxalate, malonate, succinate, and environmentally dominant organic ligands on aluminosilicate dissolution.
3. Identify which legacy conclusions remain scientifically live, which are superseded, and which modern method/solvation/speciation defects would dominate a recomputation.
4. Propose a minimum publishable slice, observables, falsification gates, estimated compute envelope, and dependencies on A3/A5/A7.

## Acceptance

- Cited scoping memo with a clear GO/NO-GO recommendation.
- Ranked reaction/species shortlist linked to `DFT-LEDGER.csv` IDs.
- Explicit pH/speciation, solvation, BSSE, and surface-vs-cluster risks.
- No DFT campaign launched by this design-only card.
- STATUS.md and this card updated in the delivering PR.

## Progress

- 2026-09-03 — completed cited scoping memo with a CONDITIONAL GO for a six-system adjacent-site oxalate pilot, ranked 18 ledger IDs, explicit falsification gates, 168-280 GPU-hour bounded full-slice estimate, and A3/A5/A7/ERW interfaces; filed A8e and launched no DFT. (hermes-macbot-zero; profile=laptop)
- 2026-09-02 — card filed from A8 because the archive contains a coherent oxalate/malonate/succinate thermochemistry set with direct ERW relevance, but no kinetic barriers. (hermes-custom-build-001; profile=workstation)

## Result

- Memo: [`docs/scoping/ligand-promoted-aluminosilicate-dissolution.md`](../../scoping/ligand-promoted-aluminosilicate-dissolution.md).
- Decision: **CONDITIONAL GO** for a six-system adjacent-site oxalate pilot; **NO-GO** for a broad oxalate/malonate/succinate campaign until the pilot clears explicit falsification gates.
- Ranked ledger shortlist: 18 recovered `DFT-LEDGER.csv` IDs in three tiers — core homologous adsorption (`ch4-oxalate/malonate/succinate-d/e`), pH/speciation overall water/hydronium rows, and hydrolysis pathway diagnostics.
- Compute envelope: 168–280 GPU-hours for the 14-system follow-on slice (218–364 with 30% retry reserve); estimate only, not a reservation.
- No DFT or QM campaign launched. Follow-up `A8e-dual-site-oxalate-pilot` filed as workstation-only Gate 1, blocked on A3 barrier infrastructure.
