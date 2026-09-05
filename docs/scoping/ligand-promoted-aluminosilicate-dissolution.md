# Ligand-promoted aluminosilicate dissolution: publication-track scoping

**Date:** 2026-09-03

**Decision:** **CONDITIONAL GO** for one bounded, dual-site oxalate mechanism pilot; **NO-GO** for a broad oxalate/malonate/succinate campaign until the pilot clears explicit falsification gates.

**Scope:** A8d scoping only. No production DFT was launched.

## Bottom line

The recovered 1999 calculations contain a real and unusually timely mechanistic hypothesis: short dicarboxylates bind an Al-bearing aluminosilicate edge more strongly than longer homologues, and oxalate can rupture an Al-O bridge in geometries where succinate does not. The numerical energies are not publication evidence. They are gas-phase/limited-microsolvation B3LYP thermochemistry on a single Si-Al dimer, with no transition states, incomplete Chapter 5 analysis, and acid/speciation choices that do not map cleanly onto any one aqueous pH.

The post-1999 literature strengthens the *question* while invalidating the old model as a direct answer.

Flow-through kaolinite experiments found up to a 30-fold oxalate enhancement at matched pH and argued for a quadratic, edge-site mechanism involving two oxalates on neighboring aluminol sites—not the single-Al topology used in the archive.[1]

Later basaltic-glass work at pH 6.4 found oxalate and desferrioxamine raised initial dissolution rates by at least fourfold and some Al/Fe release rates by as much as sevenfold, while also exposing precipitation and non-stoichiometry hazards.[6]

Contemporary kaolinite and feldspar studies keep oxalate, citrate, malate, and hydrated surface adsorption active, but the direct modern evidence for a homologous oxalate/malonate/succinate rate series remains thin.[2][3][8]

That creates a defensible publication wedge: **test whether adjacent-site ligand occupancy changes the Al-edge detachment barrier, and whether ring size survives when protonation and explicit water are treated correctly.** The result is publishable whether it confirms or falsifies the legacy ordering, provided the pilot is benchmarked to matched-pH, no-ligand controls and the experimentally proposed dual-site mechanism.

## 1. Evidence separation

### 1.1 Recovered 1999 evidence

Primary archive sources:

- [`legacy/thesis-archive/DFT-LEDGER.csv`](../../legacy/thesis-archive/DFT-LEDGER.csv) is the machine-readable authority for the recovered energies and identifiers.
- [`legacy/thesis-archive/DFT-LEDGER.md`](../../legacy/thesis-archive/DFT-LEDGER.md) documents scope, provenance, and recovery limitations.
- [`Playing-Dice-with-the-Universe.rebuilt-2026.pdf`](../../legacy/thesis-archive/text/Playing-Dice-with-the-Universe.rebuilt-2026.pdf), Chapter 4, is the human-readable scientific narrative and table source.
- [`legacy/thesis-archive/CATALOG.md`](../../legacy/thesis-archive/CATALOG.md) records the rebuilt-PDF evidence chain and the incomplete Chapter 5 state.

What the archive actually shows:

1. **A homologous adsorption ordering on one protonated dimer.** The archived Table 4.11 entries give increasingly less favorable ligand-complex formation from oxalate to malonate to succinate: `ch4-oxalate-d` (-539.897 kJ/mol), `ch4-malonate-d` (-457.325 kJ/mol), and `ch4-succinate-d` (-304.550 kJ/mol). Adding one explicit water preserves the qualitative oxalate lead (`*-e`: -531.396, -428.027, and -359.567 kJ/mol, respectively).
2. **A qualitative mechanistic contrast.** The Chapter 4 narrative reports Al-O bridge rupture for optimized oxalate and malonate structures, but not for succinate. That is a structure-level observation from the recovered minima, not a demonstrated kinetic path.
3. **A water-only thermodynamic crossing.** The reconstructed overall hydrolysis rows are favorable for oxalate (`ch4-overall-oxalate-water`, -14.157 kJ/mol) but unfavorable for malonate and succinate (`ch4-overall-malonate-water`, +66.772; `ch4-overall-succinate-water`, +108.402 kJ/mol).
4. **Strong proton sensitivity.** The hydronium bookkeeping makes all three overall rows strongly favorable (`*-hydronium`: -772.774, -691.845, -650.215 kJ/mol). This is a warning that proton chemical potential/speciation dominates the comparison, not a rate prediction.
5. **No barrier evidence.** None of the Chapter 4 ledger rows is a transition state. Chapter 5 contains a chapter shell, fragment calculations, and basis-set-superposition/counterpoise intent, but it does not contain a completed topological analysis or publishable barrier ladder.

### 1.2 Later experimental and computational evidence

The robust, directly relevant findings are:

- **Kaolinite, matched-pH rate evidence:** Cama and Ganor measured steady-state Al and Si release in far-from-equilibrium flow-through experiments at pH 2.5-3.5 and 25-70 °C. Oxalate increased rates by up to 30 times. They rejected Al-inhibition/saturation changes as the sole explanation and fit a quadratic dependence on adsorbed oxalate, proposing two ligands on neighboring edge aluminol sites.[1]
- **Kaolinite, ligand identity:** experiments comparing citric, oxalic, and malic acids found ligand-specific dissolution behavior rather than a generic "organic acidity" effect; more recent work reported complementary citrate/oxalate effects on Si and Al release.[2][3]
- **Feldspar surface relevance:** plagioclase work tied oxalate adsorption to ligand-promoted dissolution, and later water-feldspar-interface work directly characterized oxalic-acid adsorption behavior.[4][8] The older review-level caution still stands: millimolar simple-ligand concentrations often are required for large laboratory effects, and field extrapolation is not automatic.[13]
- **Basalt/ERW relevance:** at pH 6.4 and 25 °C, oxalate and a trivalent-metal siderophore accelerated synthetic basaltic-glass dissolution, with preferential Al/Fe interactions and up to sevenfold release-rate increases.[6] Oxalate-promoted forsterite dissolution and a 2025 dunite-weathering study extend the question to Mg-silicate and enhanced-weathering settings, but they do not establish the aluminosilicate edge mechanism.[7][10]
- **Hydrated-model requirement:** modern DFT adsorption work on kaolinite explicitly treats hydrated surface adsorption, underscoring that the archive's gas-phase/one-water energetics are hypothesis-generating rather than transferable free energies.[9]
- **Near-neutral and near-equilibrium behavior is a separate regime:** modern kaolinite dissolution work demonstrates that saturation state and neutral-pH conditions must be controlled rather than inferred from acidic, far-from-equilibrium experiments.[11]

## 2. What remains live, and what is superseded

| 1999 conclusion or premise | Status | Reason |
|---|---|---|
| Organic ligands can promote aluminosilicate dissolution beyond proton-only behavior. | **Live** | Supported by matched-pH kaolinite and basaltic-glass experiments.[1][6] |
| Al-carboxylate coordination at edge sites is mechanistically central. | **Live, sharpened** | The surface-coordination framework is longstanding, and later work points specifically to edge aluminol adsorption and potentially adjacent-site occupancy.[14][1] |
| Oxalate is the strongest of oxalate/malonate/succinate because five-membered chelation is least strained. | **Plausible, unvalidated** | The archive supports the ordering thermodynamically; the surveyed post-1999 dissolution literature is oxalate-heavy and does not provide a clean matched-pH three-ligand rate series. |
| A single ligand on a single Al-bearing dimer is the publication-grade mechanism. | **Superseded** | The best direct kinetic evidence proposes two neighboring edge sites and a quadratic low-coverage rate law.[1] |
| Reaction energies can stand in for dissolution kinetics. | **Superseded** | A publication claim requires located first-order saddles, IRC connectivity, and activation free energies. |
| Fully deprotonated ligand + isolated cluster is a sufficient aqueous model. | **Superseded** | Protonation, surface charge, Al-ligand aqueous complexation, explicit hydration, counterions, and standard states can reorder charged complexes. |
| The Chapter 5 bond-topology argument was completed. | **Not supported** | Recovered Chapter 5 is incomplete; no finished topological results or barrier evidence were recovered. |

## 3. Ranked legacy ledger shortlist

Ranking is by *information value for a modern recomputation*, not by the magnitude of the old electronic energy.

### Tier 1 — core homologous adsorption set

1. `ch4-oxalate-d` / `ch4-oxalate-e` — strongest archive signal; one-water sensitivity and reported Al-O rupture.
2. `ch4-malonate-d` / `ch4-malonate-e` — nearest homolog; tests whether the six-membered chelate remains mechanistically competent.
3. `ch4-succinate-d` / `ch4-succinate-e` — negative-control homologue; no archived bridge rupture.

**Why first:** these six rows provide the cleanest common reaction definition across all three ligands and map directly onto a modern single-site versus adjacent-site adsorption screen.

### Tier 2 — pH/speciation stress test

4. `ch4-overall-oxalate-water`
5. `ch4-overall-malonate-water`
6. `ch4-overall-succinate-water`
7. `ch4-overall-oxalate-hydronium`
8. `ch4-overall-malonate-hydronium`
9. `ch4-overall-succinate-hydronium`

**Why second:** the water rows create a qualitative sign change across the homologues, while the hydronium rows expose how strongly proton bookkeeping drives the legacy numbers. They should be reinterpreted through pH-dependent free-energy cycles, not recomputed verbatim as isolated charged reactions.

### Tier 3 — pathway and interaction diagnostics

10. `ch4-oxalate-hydrolysis-b` / `ch4-oxalate-hydrolysis-c`
11. `ch4-malonate-hydrolysis-b` / `ch4-malonate-hydrolysis-c`
12. `ch4-succinate-hydrolysis-b` / `ch4-succinate-hydrolysis-c`

**Why third:** these rows isolate water attack and proton-assisted hydrolysis bookkeeping, but the archive contains minima and reaction energies rather than validated reaction paths. Use their structures only as starting guesses.

The standalone `strain` and `counterpoise_*` systems are useful diagnostics after a modern pathway exists; they are not first-wave reaction targets because the recovered Chapter 5 analysis is incomplete.

## 4. Falsifiable mechanistic hypothesis

> At a hydrated kaolinite-like Al edge, oxalate promotion arises primarily when two neighboring aluminol sites are simultaneously occupied; that occupancy lowers the activation free energy for Al-O-Si bridge cleavage/Al detachment relative to both an unliganded edge and a singly occupied edge. Malonate should produce a smaller barrier reduction, and succinate should provide a weak or null homologous control after matched protonation and hydration.

This hypothesis deliberately separates four effects that the archive mixed together:

1. proton-promoted hydrolysis;
2. ligand adsorption/coverage;
3. ligand-assisted bond rupture at the surface;
4. aqueous Al-ligand stabilization after detachment.

A negative result is informative. If dual occupancy does not lower the connected detachment barrier, the most direct modern interpretation of the archive is falsified and the broad homologous ligand track stops.

## 5. Minimum publishable recomputation slice

### 5.1 Surface and chemistry

- **Surface motif:** a hydrated kaolinite-like edge cluster with **two adjacent edge Al sites** and at least one bridging Al-O-Si bond; derive protonation and freeze/relax conventions from A3 rather than reusing the isolated 1999 dimer unchanged.
- **Matched control:** the identical hydrated edge, charge, pH-state convention, and water count without ligand.
- **Ligands:** oxalate is the mechanistic lead; malonate and succinate are homologous controls. Citrate and malate belong in the experimental-comparison table but not in the minimum QM pilot unless the oxalate gate passes.
- **Concentration/coverage bridge:** represent no ligand, one bound ligand, and two neighboring bound ligands. Map structures to a surface-coverage model rather than treating cluster stoichiometry as bulk concentration.
- **pH states:** two experimentally anchored regimes: acidic kaolinite conditions near pH 3 and a near-neutral weathering regime near pH 6-6.5. Before structure generation, calculate aqueous ligand, Al-ligand, and edge-site protonation populations; carry only species with material population or a defined mechanistic role.
- **Water:** an explicit first/second hydration shell sufficient to support proton relay and Al coordination, plus a consistent continuum correction. Require convergence against one larger-water-shell spot check.

### 5.2 Reaction coordinate and observables

For every retained model, locate the reactant complex, first-order saddle, and product complex for the same Al-O-Si cleavage/Al-detachment coordinate. Verify saddle identity by one imaginary frequency and full IRC connectivity. If water substitution is a distinct prerequisite, report it as a separate step rather than hiding it in a net reaction energy.

Primary observables:

- ligand adsorption free energy on the hydrated edge;
- activation free energy `ΔG‡` for the connected cleavage/detachment step;
- reaction free energy `ΔG` for that step;
- Al-O-Si, Al-O(ligand), and adjacent-site geometric changes along the IRC;
- proton-transfer topology and coordination number;
- model sensitivity to water-shell size, edge protonation, and production single-point method;
- predicted coverage/rate-law shape for qualitative comparison with the experimental quadratic dependence.[1]

Experimental/field observables to preserve for A5/A7 handoff:

- steady-state Si and Al release rates, not Al alone;
- Al/Si stoichiometry and dissolved Al-ligand speciation;
- adsorbed ligand or surface coverage where measurable;
- pH, ionic strength, temperature, saturation index, reactive surface area, and ligand concentration;
- Ca/Fe/Mg release and secondary precipitates, especially Ca-oxalate in basalt-bearing systems.[6]

### 5.3 Staged matrix and compute envelope

**Gate 1 — six-system oxalate pilot:** no ligand, one oxalate, and two adjacent oxalates at each of the two pH-state conventions. This is the minimum direct test of single-site versus dual-site promotion.

**Gate 2 — eight-system homolog extension:** if Gate 1 passes, add single and paired malonate and succinate occupancy at both pH-state conventions.

The full bounded slice is therefore 14 reaction systems. A3 budgets about 10 GPU-hours per reaction cell for its smaller B3LYP barrier tier. Allowing 12-20 GPU-hours per larger ligand/dual-site system gives a planning envelope of **168-280 GPU-hours**, or **218-364 GPU-hours including a 30% retry reserve**. This is an estimate, not a reservation; benchmark one unliganded and one dual-oxalate path before committing the extension.

Recommended method ladder, aligned with A2/A3:

1. conformer and protonation enumeration with cheap prescreening;
2. r2SCAN-3c optimizations and frequency checks where feasible;
3. production single points using the A2 protocol on all connected stationary points;
4. BSSE/basis/solvation sensitivity checks on the matched unliganded and dual-oxalate pair;
5. no periodic escalation unless the edge-cluster result changes sign under boundary-condition checks.

## 6. Promotion criteria and stop conditions

### Pilot promotes to the homolog extension only if all pass

1. **Connected kinetics:** the dual-oxalate saddle has one relevant imaginary mode and IRC endpoints matching the intended reactant/product families.
2. **Material barrier effect:** dual oxalate lowers `ΔG‡` by at least **10 kJ/mol** relative to the matched unliganded control, and the uncertainty/sensitivity envelope does not cross zero.
3. **Occupancy discrimination:** the paired-occupancy result is measurably more promoting than single occupancy, consistent in sign with the experimentally inferred nonlinear coverage dependence.[1]
4. **Model robustness:** the sign survives the larger explicit-water shell and production-tier single point.
5. **No hidden bookkeeping win:** the claimed lowering is not solely an aqueous Al-complexation or proton-count difference; all compared pathways close the same thermodynamic cycle.

### The dedicated track stops or narrows if any occurs

- no stable adjacent-site oxalate reactant exists under the retained pH/speciation model;
- no connected saddle can be located after the bounded attempt set;
- `|ΔΔG‡| < 10 kJ/mol` or method/water/protonation sensitivity exceeds the effect;
- single and paired occupancy are indistinguishable, defeating the proposed quadratic mechanism test;
- malonate/succinate ordering changes with routine hydration/protonation choices, making a ring-size series non-identifiable;
- precipitation or ligand depletion dominates the ERW-facing mass balance, in which case the work becomes a geochemical-speciation study rather than a surface-barrier track.

A pilot failure is a **NO-GO for a broad DFT campaign**, not permission to keep adding ligands or surface models.

## 7. Key model risks

1. **Wrong surface topology:** a single Al dimer cannot test the best later mechanism; the model must expose adjacent edge aluminol sites.
2. **pH/speciation aliasing:** oxalic, malonic, and succinic acids have different protonation distributions, and Al binding shifts them further. Comparing one arbitrarily chosen charge state can manufacture a trend.
3. **Charged-cluster artifacts:** diffuse basis functions, counterions, cavity choice, and standard-state corrections matter strongly for multiply charged complexes.
4. **Water-network non-uniqueness:** one explicit water is insufficient for proton relays and six-coordinate Al; conformers and water-shell size can move barriers.
5. **Adsorption versus detachment:** a favorable adsorption energy does not imply a lower dissolution barrier and may instead poison/reactively passivate a site.
6. **Surface coverage:** the experimental quadratic law is a collective/coverage observation. One cluster at fixed stoichiometry cannot recover it without an explicit occupancy model.
7. **Secondary phases:** Ca-oxalate and Al-bearing precipitates can alter solution release, particularly in basalt/ERW contexts.[6]
8. **Transfer across solids:** kaolinite edge Al is octahedral, feldspar Al is tetrahedral, and basaltic glass is disordered; a common ligand effect does not imply an identical elementary path.[1]
9. **Field relevance:** millimolar ligand effects in clean reactors do not automatically survive soils, competing cations, microbial turnover, or transport.[13]

## 8. Dependencies and interfaces

- **A3 barrier ladder — hard computational dependency.** Reuse its model-provenance, protonation enumeration, endpoint-family, TS/IRC, and barrier-receipt discipline. A8e needs an adjacent-Al edge derivative rather than silently substituting the old dimer. A3 should remain the owner of generic barrier infrastructure.
- **A5 field-lab discrepancy — translation dependency, not a pilot blocker.** Export `ΔΔG‡`, ligand coverage assumptions, Al/Si stoichiometry, and pH/speciation envelopes as mechanism modifiers that can be compared against aging/transport observables. Do not convert cluster adsorption energy directly into a field rate multiplier.
- **A7 kinetics database — evidence/provenance dependency.** Add ligand identity, total/free ligand concentration, adsorbed coverage, pH, metal-ligand speciation method, and no-ligand matched-control fields before using literature rates for calibration. Existing acid/neutral/base bins are not sufficient for ligand-promoted terms.
- **Enhanced rock weathering — downstream application filter.** Treat oxalate/citrate/siderophore chemistry as a possible interface-local accelerator, while preserving cation complexation, alkalinity yield, ligand persistence, and secondary-precipitate penalties. The dunite literature establishes application relevance, not transferability of the Al-edge mechanism.[10]

## 9. Decision

**CONDITIONAL GO:** create one workstation-bound implementation card for the six-system adjacent-site oxalate pilot, blocked on A3's barrier infrastructure. Do not begin the eight-system malonate/succinate extension until the oxalate pilot passes every promotion criterion above.

This decision is stronger than "recompute the thesis." The legacy archive contributes starting geometries, a homologous negative control, and an old ring-size hypothesis. The publishable question comes from confronting those artifacts with the later two-site kinetic mechanism, explicit hydration, pH-resolved speciation, and connected activation barriers.

## Sources

[1] https://doi.org/10.1016/j.gca.2006.01.028 — Cama & Ganor (2006): oxalate catalysis of kaolinite dissolution
[2] https://doi.org/10.1016/j.jcis.2005.04.066 — Wang et al. (2005): kaolinite dissolution by citric, oxalic, and malic acids
[3] https://doi.org/10.1016/j.clay.2020.105756 — Lin et al. (2020): citrate/oxalate synergy in kaolinite dissolution
[4] https://doi.org/10.1021/es980258d — Stillings et al. (1998): oxalate adsorption on plagioclase
[6] https://doi.org/10.1016/j.gca.2015.04.025 — Perez et al. (2015): iron chelators and basaltic-glass dissolution
[7] https://doi.org/10.1016/j.gca.2007.12.026 — Olsen & Rimstidt (2008): oxalate-promoted forsterite dissolution
[8] https://doi.org/10.1007/s10450-019-00111-8 — Xue et al. (2019): oxalic acid at the water-feldspar interface
[9] https://doi.org/10.1016/j.clay.2021.106075 — Heimann et al. (2021): hydrated adsorption on kaolinite by DFT
[10] https://doi.org/10.1016/j.ceja.2025.100902 — Niron et al. (2025): organic ligands in dunite weathering and carbon capture
[11] https://doi.org/10.1016/j.clay.2019.105284 — Gong et al. (2019): kaolinite dissolution near equilibrium and neutral pH
[13] https://www.sciencedirect.com/science/article/pii/S092777579603720X — Drever & Stillings (1997): organic acids in mineral weathering
[14] https://www.sciencedirect.com/science/article/pii/0016703786902437 — Furrer & Stumm (1986): coordination chemistry of weathering I
