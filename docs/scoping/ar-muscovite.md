# Closing the 1998 circle: a defect-resolved lattice KMC of Ar release from muscovite

*Scoping document, 2026-08-19. Founding document: Sletten & Onstott, GCA 62, 123–141 (1998). Toolchain: petra (declarative lattice-KMC engine, this repo) + quarry (QM barrier pipeline). Prepared for Victor Sletten.*

---

## 1. The question and the 27-year arc

**1998 — the empirical claim.** Sletten & Onstott (1998) heated a single, exceptionally well-characterized Benson Mines muscovite (3000 µm disks plus 150–180 µm and 52–64 µm size fractions) isothermally at 500 °C and 700 °C in vacuo, monitored dehydroxylation by XRD and delamination by SEM, and compared Ar loss against volume-diffusion models. The findings, stated as the paper stated them:

- In-vacuo Ar release below ~900 °C is **rate-limited by dehydroxylation and delamination**, not volume diffusion. Cumulative-loss curves are non-Fickian (two-stage, polymer-sorption-like); D/a² values track the *dehydroxylation* K/a² values through time (Fig. 8 of the paper).
- **Atmospheric Ar and radiogenic ⁴⁰Ar\* occupy different structural reservoirs**: D/a²(³⁶Ar)/D/a²(⁴⁰Ar\*) ranges 10–10⁵ early in the runs — atmospheric Ar sits in interlayer *extended defects* (Guinier–Preston-like zones), ⁴⁰Ar\* in interlayer *vacancies*, and recoil-implanted ³⁹Ar is uniform (some fraction driven into the octahedral layer, retarding it relative to ⁴⁰Ar\* and producing anomalously old first steps).
- **Staircase age spectra can be artifacts of differential release from alteration/intergrown phases**, not fossil ⁴⁰Ar\* gradients — the hydrothermally treated split's staircase did not match its laser-probe map (Fig. 13).
- Bottom line: *recovery of fossil ⁴⁰Ar\* gradients from muscovite by incremental heating is not possible*; the mechanism in the lab is not the mechanism in nature.

The paper's model was verbal-plus-continuum: a two-stage story (defect purge → dehydroxylate slab diffusion) fit with cylindrical/slab/2-D reaction approximations. It could not *compute* a spectrum from the mechanism. That is the open question: **can a mechanistic, structurally-evolving model quantitatively generate the observed release systematics — and therefore tell us what step-heating spectra of white mica actually mean?**

**2022–2024 — the ideal-lattice theory.** The Orléans group closed a different flank. [Nteme, Scaillet, Brault & Tassan-Got (2022, GCA 331, 123–142)](https://www.sciencedirect.com/science/article/abs/pii/S0016703722002198) ran classical MD of ⁴⁰Ar recoil plus NEB+TST migration barriers in 2M₁ muscovite: recoiled ⁴⁰Ar resides in the interlayer near its production site, and the cheapest migration pathway is a **divacancy mechanism in the interlayer, with a migration enthalpy of ~66 kcal/mol**. Crucially, ideal-lattice diffusivities computed this way are **too slow** to explain measured Ar mobility — lab and nature are both faster than the perfect crystal allows. Their companion paper ([Nteme et al. 2023, GCA — defect-controlled ⁴⁰Ar diffusion-domain structure from ⁴⁰Ar/³⁹Ar crystal-mapping](https://www.sciencedirect.com/science/article/abs/pii/S0016703722006512)) mapped a slowly-cooled Harney Peak muscovite into ~250–300 µm retentive sub-grain cores separated by high-diffusivity zones where diffusivity is enhanced **up to six orders of magnitude**, tied by TEM to lenticular voids and basal partings. The [2024 phlogopite PVT calibration (Nteme et al., GCA 382, 103–122)](https://www.sciencedirect.com/science/article/abs/pii/S0016703724003454) extended the machinery: pristine-lattice Ar diffusivity is orders of magnitude slower than hydrothermal estimates, with strong pressure dependence — again pointing at defects as the real carrier.

**The uncomputed middle.** So we now have (a) 1998's empirical demonstration that *structural breakdown* controls in-vacuo release and that distinct defect reservoirs carry distinct Ar species, and (b) 2022–24's first-principles barriers proving the *ideal lattice cannot do the job*, with defects invoked but not simulated. Nobody has built the middle layer: a **lattice model in which the crystal's structural state evolves (dehydroxylation fronts, delamination, defect populations) while Ar migrates through it with siting-specific barriers**, integrated over a laboratory T(t) schedule to emit a synthetic ⁴⁰Ar/³⁹Ar age spectrum. That model is exactly what petra was built to express, and the 1998 paper is its ready-made validation dataset.

There is also live controversy to arbitrate. [Villa's "In vacuo release of Ar from minerals" series (part 1: Chemical Geology, 2021](https://www.sciencedirect.com/science/article/pii/S0009254121000206); [part 3: apatite rare gases, 2023](https://www.researchgate.net/publication/376050172_The_in_vacuo_release_of_Ar_and_rare_gases_from_minerals_3_The_degassing_of_He_Ne_Ar_Kr_and_Xe_from_irradiated_apatite); [part 4: polymineralic samples, Geosciences 16(8):298, 2026](https://doi.org/10.3390/geosciences16080298)) argues that **Ar/Xe decoupling rules out delamination as the dominant degassing mechanism** in mica, that recoiled rare gases mostly sit *inside the T–O–T layers*, that white micas release Ar *above* the literature dehydroxylation temperature, and that a common apparent Ea of ~300 kJ/mol across tecto-, phyllo- and cyclosilicates demands a shared kinetics. This is a direct challenge to the 1998 delamination interpretation — and it is exactly the kind of dispute a mechanism-resolved simulation can settle, because the competing hypotheses predict *different* isotope- and species-resolved release curves.

## 2. Literature map

**Diffusion calibrations (continuum).**
- [Harrison et al. 2009, GCA 73, 1039–1051 — "Diffusion of ⁴⁰Ar in muscovite"](http://sims.ess.ucla.edu/PDF/Harrison_et_al_2009.pdf): hydrothermal piston-cylinder at 10 kbar, 600–730 °C; Ea = 63 ± 7 kcal/mol; activation volume ~14 cm³/mol (20 kbar runs ~10× slower); the field's standard muscovite calibration, Tc ≈ 405–445 °C. Note the near-coincidence with Nteme's 66 kcal/mol divacancy barrier — suspicious and interesting.
- [Hames & Bowring 1994, EPSL 124, 161–167](https://www.sciencedirect.com/science/article/abs/pii/0012821X94000794): empirical case for cylindrical geometry with the effective dimension being the in-plane grain radius.
- Robbins 1972 (Brown MS thesis): the older hydrothermal dataset the 1998 paper leaned on; never formally published — a standing embarrassment of the calibration chain worth noting in the paper.

**Defect/multipath line.**
- [Lee 1995, Contrib. Mineral. Petrol. 120, 60–82 — "Multipath diffusion in geochronology"](https://link.springer.com/article/10.1007/BF00311008): the canonical coupled lattice+short-circuit continuum model; successors applied it to biotite/hornblende. Continuum, two fixed pathways, no structural evolution.
- [Kramar, Cosca, Buffat & Baumgartner 2003, GSL SP 220 — stacking-fault-enhanced Ar diffusion in deformed muscovite](https://www.lyellcollection.org/doi/abs/10.1144/gsl.sp.2003.220.01.15): TEM-documented stacking faults as fast pathways; dilatation argument.
- [Nteme et al. 2022](https://www.sciencedirect.com/science/article/abs/pii/S0016703722002198), [2023](https://www.sciencedirect.com/science/article/abs/pii/S0016703722006512), [2024](https://www.sciencedirect.com/science/article/abs/pii/S0016703724003454) as above — the modern atomistic anchor (open-access versions on [HAL](https://insu.hal.science/insu-03668854)).
- Analogous multi-scale noble-gas work outside mica: [nanopore trapping of He in quartz (ACS Earth & Space Chem. 2020)](https://pubs.acs.org/doi/abs/10.1021/acsearthspacechem.0c00187), multi-method He/Ne in zircon — the "defects control noble gas retention" theme is now general.

**Release-mechanism / in-vacuo line (the 1998 paper's direct descendants).**
- [Villa series parts 1–4 (2021–2026)](https://doi.org/10.3390/geosciences16080298) as above — the delamination challenge; three-isotope correlation diagrams for polymineralic mixtures.
- Earlier: Gaber et al. 1988, Lo 1990, Lee 1993, Wartho 1995 (amphibole/biotite in-vacuo breakdown), all cited in and consistent with the 1998 mechanism.

**Age-spectrum modeling practice.**
- [Lister & Forster 2024, Aust. J. Earth Sci. 71, 950–964 — quantitative modelling of muscovite age-plateau morphology](https://www.tandfonline.com/doi/full/10.1080/08120099.2024.2421569): MacArgon/DiffArg-class forward modeling; the "method of asymptotes and limits"; ANU keeps step-heating interpretation alive as *the* deformation-dating tool.
- [MDD parameter-optimization tooling (Gchron 2024 preprint)](https://gchron.copernicus.org/preprints/gchron-2024-11/gchron-2024-11.pdf): the K-feldspar MDD world is still Crank-Nicolson + Monte Carlo over thermal paths — all continuum, domains fixed a priori.
- Practitioner state of the art for white mica interpretation: deformation/recrystallization control rather than pure thermal closure — [Barnes et al. 2023, J. Metam. Geol.](https://dx.doi.org/10.1111/jmg.12739) (in-situ vs single-grain fusion), [Beaudoin et al. 2020, Tectonics](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020TC006246), [combined in-situ Rb/Sr + Ar/Ar (Lithos 2024)](https://www.sciencedirect.com/science/article/pii/S0024493724002007), [Kellett/Schneider-school "approaches and best practices" chapter](https://www.sciencedirect.com/science/article/abs/pii/B9780443188039000080), [excess-Ar tests vs Rb/Sr in Himalayan mica (EPSL 2023)](https://www.sciencedirect.com/science/article/abs/pii/S0012821X23000717), [Erzgebirge microtectonic control of white-mica dates (GCA 2021)](https://www.sciencedirect.com/science/article/abs/pii/S0016703721005329).

**Who is active, 2020–2026.** Orléans/ISTO (Scaillet, Nteme — atomistic + crystal-mapping; the closest group to this project), Bern/Milano (Villa — in-vacuo mechanism revisionism), ANU (Forster & Lister — spectra morphology, eArgon), UCLA (Harrison lineage — MDD, calibrations), Durham/Open U. (Warren — excess Ar, subduction micas), Ottawa (Schneider) and GSC (Kellett) — deformation dating practice, USGS (Cosca), WiscAr (Jicha), Melbourne/Curtin in-situ labs. Nobody in this list runs lattice KMC.

**KMC adjacent.** Lattice/off-lattice KMC exists for mineral *dissolution* (Kurganskaya & Luttge; [clay-nanoparticle dissolution KMC, Minerals 14:900, 2024](https://doi.org/10.3390/min14090900)) and coupled chemo-mechanics ([MASKE](https://arxiv.org/pdf/2311.06533)), and for gas–surface problems in catalysis and astrochemistry. **No published KMC couples structural breakdown of a mineral with noble-gas migration, and no one has emitted synthetic ⁴⁰Ar/³⁹Ar step-heating spectra from any atomistic or lattice model.** The 2024 ANU modeling and all MDD work are continuum diffusion with static domain structure.

## 3. The precise gap

Every existing tool holds one of these fixed: either the *structure* (continuum diffusion, MDD, multipath — domains and pathways are static parameters) or the *time integration* (atomistic MD/NEB — barriers, but no mesoscale kinetics, no spectra). The gap, precisely:

> **A mesoscale, defect-resolved lattice model of muscovite in which (i) dehydroxylation propagates as an autocatalytic/steric-limited reaction, (ii) delamination opens fast-release interfaces along the stacking axis as a *consequence* of dehydroxylation strain, (iii) Ar species (⁴⁰Ar\*, ³⁹Ar_K, ³⁶Ar_atm — and Xe as discriminant) occupy siting-specific reservoirs (interlayer vacancies, extended defects, T–O–T interior) with mechanism-specific barriers, and (iv) a laboratory T(t) staircase drives the whole system — emitting cumulative release curves and synthetic age spectra that can be laid directly over 1998 Figs. 3–11 and Villa's isotope-correlation data.**

The payoff question it answers is the one the field still argues about: *what does a muscovite staircase spectrum mean?* The 1998 answer (mixing + differential structural release), the Villa answer (in-layer siting, universal ~300 kJ/mol kinetics, no delamination), and the ANU answer (recoverable geological signal, with care) make different predictions once you can compute spectra from mechanisms. That's a falsification engine, not a parameter fit.

## 4. The muscovite petra deck sketch

Petra today (DESIGN.md; RFC-001): declarative TOML deck → site kinds with (occupant, config) states, guarded local rewrite rules with Arrhenius/Eyring rates, `by_count`/`when` environment modifiers, labeled bonds, multilayer cells (kaolinite-multilayer deck exists), `[[defects]]` analytic strain fields with per-reaction strain coupling, `[[init]]` probability passes (seeded, deterministic), ensembles, event-log outputs; schema v2 adds `[execution]`/`[observables]`.

**Structure.** 2M₁ muscovite reduced to the KMC-relevant graph: per layer, a T–O–T slab represented by (a) octahedral OH sites (the dehydroxylation actors), (b) interlayer K sites (the Ar/K/vacancy gallery), (c) interlayer-interface sites on the c-axis (the delamination actors). A `[[cell]]` with these three kinds, tiled `na × nb × n_layers`, open boundary on c (basal surfaces), periodic or open in-plane depending on whether we model disk interiors or flake edges.

**Site kinds and states.**
- `K_gallery`: occupant ∈ {K, vacant, Ar40, Ar39, Ar36, Xe}; config tracks siting class: `normal`, `extended_defect` (assigned by `[[init]]` patches — lenticular-void zones à la Nteme 2023), `adjacent_dehydroxylate`.
- `OH_site`: config ∈ {`OH`, `O_residual` (post-dehydroxylation, Al 5-coordinate)} — paired via labeled bonds so the OH+OH → H₂O + O rewrite is a two-site effect.
- `interface`: config ∈ {`bonded`, `delaminated`}.
- Octahedral-interior trap sites for the recoil-implanted ³⁹Ar fraction (the 1998 explanation for ³⁹Ar retardation) — occupant Ar39, high escape barrier.

**Rules (with barrier provenance).**

| Rule | Form | Barrier | Status |
|---|---|---|---|
| Ar interlayer hop (divacancy) | guard: ≥2 vacant K-neighbors; hop effect | Em ≈ 66 kcal/mol | **published** (Nteme 2022) |
| Ar hop within extended-defect zone | same, config-gated `when` | ~26 kcal/mol lower at 700 °C if the 10⁶× enhancement (Nteme 2023) is barrier-equivalent | **anchored, needs computing** |
| K vacancy migration | hop rule | in Nteme 2022's supplement (vacancy mechanisms compared) | **published/partial** |
| Dehydroxylation | OH-pair → H₂O(removed) + O_residual; `by_count` on neighboring `O_residual` *raising* Ea (Guggenheim steric mechanism → diminishing rate, exactly the 1998 Fig. 8 behavior) | apparent Ea 57 kcal/mol (Kodama & Brydon 1968, used in 1998) up to 251–290 kJ/mol modern in-situ XRPD ([Gualtieri-school kinetics](https://link.springer.com/article/10.1007/s002690050197)); environment-resolved table | **bulk published; per-environment table needs computing** |
| Delamination propagation | `interface` bonded→delaminated; guard: local dehydroxylation fraction ≥ threshold via distance-2 counts; strain coupling via existing `[[defects]]`/strain machinery | none published | **needs computing or calibration** |
| Ar release | Ar adjacent to `delaminated` interface or open c-boundary → removed (fast desorption) | small; effectively instant | modeling choice |
| Ar in dehydroxylate lattice | hop with collapsed-interlayer barrier | none published (1998 used Robbins 1972 as proxy) | **needs computing — the key missing number** |
| ³⁹Ar recoil siting | `[[init]]` probability pass: uniform over gallery + fraction into octahedral traps | geometry from Nteme 2022 recoil MD | **published** |
| Xe analogs of all hops | same rules, Xe occupant | none for muscovite | **needs computing (the Villa discriminant)** |

**Computing the missing barriers.** The interlayer problems are periodic-slab problems, not cluster problems — quarry's current cluster workflow (kaolinite edge clusters) is the wrong shape for them. Two viable routes: (1) replicate Nteme's classical-potential NEB setup (cheap, directly comparable, runs on the workstation GPU/CPU today) for dehydroxylate-interlayer Ar, defect-zone Ar, and Xe; (2) periodic DFT (CP2K, PBE-D3, ~2×2×1 supercells, climbing-image NEB) for the subset where classical potentials are dubious (post-dehydroxylation 5-coordinate Al environments). Route 1 first; route 2 as calibration. This is a bounded campaign: ~6–10 distinct barrier types × 2 species.

**The T(t) engine feature (required work, flagged).** Petra runs fixed T per run. Step-heating is *piecewise-isothermal by construction*, so the minimal feature is exact: `[execution.schedule]` as a list of `(T, duration)` segments; at each boundary, recompile the rate tables (same code path as load-time) and rebuild the Fenwick tree. No within-segment T-dependence needed. This slots cleanly into RFC-001's `[execution]` section and should be specified as a small RFC amendment; estimated effort is days, not weeks. A second engine need is already scheduled elsewhere: cumulative per-species release counting is just the event log + a post-processor (exists today); a first-class `release_curve` observable belongs in the B4 observables work.

**Observables → synthetic spectra.** Cumulative ³⁹Ar and ⁴⁰Ar release per schedule segment → apparent age per step → staircase plot; D/a² per segment via the same Fickian inversion formulas the 1998 paper used (deliberately inverting the synthetic data with the *wrong* continuum model, to reproduce what a lab would infer). ³⁶Ar/⁴⁰Ar contamination-vs-step and Ar/Xe release decoupling come free from species-resolved counting.

**Scale honesty.** A 3000 µm disk is ~10¹⁹ atoms; petra's practical range is 10⁶–10⁸ sites (~0.1–1 µm slabs — which is, usefully, exactly the 0.5 µm dehydroxylate-lamella scale the 1998 slab model inferred). Strategy: simulate representative volumes (single gallery stacks, n_layers ~ 10–100), demonstrate the *mechanism-level* kinetics (two-stage curves, isotope decoupling, staircase generation), and bridge to grain scale either by the fine-grained size fractions (52–64 µm data in 1998 — closest to simulable) or by extracting effective time-dependent release kinetics per structural state and feeding a 1-D grain-scale integrator. State this limitation in the paper before a reviewer does.

## 5. Testable claims, ranked

1. **Dehydroxylation-limited release reproduces the two-stage non-Fickian loss curves.** Observable: cumulative Ar fraction vs √t at 500/700 °C, and D/a²(t) inverted as in 1998. Validation: 1998 Figs. 6, 7, 9, 10 + the paper's own XRD K/a² kinetics (Fig. 8), Kodama & Brydon / Vedder & Wilkins dehydroxylation data.
2. **Distinct reservoirs produce the observed ³⁶Ar/⁴⁰Ar\* diffusivity ratio decay** (10–10⁵ → ~1). Observable: D/a²(³⁶)/D/a²(⁴⁰) vs time by grain size. Validation: 1998 Fig. 11a; contamination-vs-step from Table/Fig. 3 systematics.
3. **A structurally-evolving lattice generates staircase age spectra without any fossil ⁴⁰Ar\* gradient** — the 1998 central claim, computed. Observable: synthetic spectrum from a homogeneous-age lattice with defect + alteration-phase reservoirs, mirrored against the hydrothermally treated split (1998 Fig. 3d) and Wijbrans & McDougall-type staircases. This is the headline figure of paper 1.
4. **Delamination vs in-layer siting are distinguishable by Ar/Xe decoupling** — the Villa arbitration. Observable: simulated Ar and Xe release curves under both siting hypotheses; whichever matches the [Villa 2021–2026 series data](https://doi.org/10.3390/geosciences16080298) wins. Honest risk: this test may *refute* the delamination-dominant reading of 1998 — which is fine; either outcome is a result, and the founding author testing his own 28-year-old interpretation is a good look.
5. **Recoil siting explains the anomalously old first steps** (³⁹Ar octahedral trapping). Observable: first-4%-of-³⁹Ar age anomaly vs trap fraction. Validation: 1998 Fig. 5; Nteme 2022 recoil distributions.
6. **Grain-size dependence of the delamination fraction** (1998 Fig. 11b, Brandt & Voronovskiy). Observable: released-fraction-in-unstable-region vs lattice extent. Validation: the three 1998 size fractions.
7. (Stretch) **Hydrothermal vs in-vacuo contrast**: with dehydroxylation suppressed (H₂O chemical-potential term), release collapses to slow lattice diffusion — consistency check against Harrison 2009 and 1998's 2-kbar run, and the bridge from lab mechanism to nature.

## 6. Phased plan (one person + agents + one GPU)

**Phase 0 — engine (petra work, ~1–2 weeks of agent time).** `[execution.schedule]` piecewise-isothermal T(t) (RFC amendment + implementation + golden gate: two-segment schedule reproduces two independent isothermal runs seeded identically). Release-count post-processor over the existing JSONL event log. *Exists already:* multilayer cells, labeled-bond pair rewrites, `by_count` environment tables, `[[init]]` probability passes, defect strain fields, ensembles, deterministic replay.

**Phase 1 — minimal deck, published numbers only (~2–4 weeks).** Single-gallery muscovite deck: K/vacancy/Ar40 gallery with divacancy hops (66 kcal/mol), OH pairs with bulk dehydroxylation kinetics, c-axis open boundary. Deliverable: isothermal 500/700 °C release curves; qualitative match to two-stage behavior. Milestone gate: D/a²(t) inversion shows the 1998 rise-then-fall signature.

**Phase 2 — the full mechanism (~1–2 months).** Add extended-defect zones (init patches), delamination interface rules, ³⁹Ar/³⁶Ar species and recoil init, T(t) staircase schedules, synthetic age spectra. Sweep lattice size for the grain-size claims. This phase produces the paper-1 figure set (claims 1–3, 5, 6).

**Phase 3 — barrier campaign (parallel with 2; GPU).** Classical-NEB replication of Nteme's setup, then the missing barriers: Ar-in-dehydroxylate, Ar-near-defect, Xe everything, environment-resolved dehydroxylation table. CP2K spot checks. Everything provenance-stamped per house rules. Feeds claim 4 and upgrades phases 1–2 numbers.

**Phase 4 — paper 1.** *"Synthetic ⁴⁰Ar/³⁹Ar step-heating spectra from a structurally-evolving lattice model of muscovite"* — the minimal, self-contained first paper: mechanism → spectra, validated against Sletten & Onstott (1998) release data and confronted with the Villa Ar/Xe test. Not a new calibration, not a thermochronology overhaul — a quantitative test of the 1998 interpretation, 28 years later, by its author. Paper 2 (later): geological implications — what staircase spectra can and cannot record, with the ANU forward-modeling community as the audience.

## 7. Risks and reviewer skepticism

- **"Your barriers are classical-potential numbers."** Mitigate: inherit Nteme's published values where possible (their credibility, not ours), DFT spot-checks, sensitivity analysis showing which conclusions survive ±5 kcal/mol.
- **"KMC lattice ≪ grain; you can't claim grain-scale spectra."** Mitigate: lead with the fine-grained fractions; frame large-grain results as mechanism demonstrations plus continuum bridging; never claim a simulated 3 mm disk.
- **"Delamination as a lattice rule is a cartoon."** True — it's an interface-state idealization of a mechanical process. Mitigate: calibrate the threshold/rate against the 1998 SEM + XRD conversion data and present it as a parameterized hypothesis being tested, not derived physics.
- **"Villa has already shown delamination is wrong."** The model doesn't presuppose the 1998 mechanism — it can run *both* siting hypotheses and report which reproduces the joint Ar+Xe+dehydroxylation data. Design the paper so this is a feature.
- **Is Orléans already heading here?** They have the barriers (2022, 2024), the defect maps (2023), and continuum domain modeling — the natural next step for them is exactly this KMC. Nothing published yet (checked 2026-08). Speed matters; so does the fact that petra already exists and they'd be building from scratch. Collaboration is also a live option — they hold the barrier expertise, Victor holds the engine and the 1998 dataset.
- **Timescale stiffness**: fast defect-zone hops vs slow dehydroxylation could waste steps; petra's stiffness detection exists, superbasin acceleration (Omnibus B7) may be needed for the 45-day hydrothermal comparison. In-vacuo lab schedules (minutes–days) are likely tractable as-is.
- **The 1998 raw data**: figures exist in the paper; step tables may need digitizing (or recovering from Victor's own archives — worth checking for the original Varian-MAT GD150 data sheets).

## 8. Venues

Paper 1: **GCA** (the 1998 paper's home, the 2022/23/24 papers' home — the conversation is already there) or **Chemical Geology** (Villa series home; mechanism-focused). **Geochronology (Copernicus)** as the modern, open, methods-friendly alternative. EPSL if the headline sharpens to "step-heating spectra of white mica cannot record fossil gradients — quantitatively." Tool/engine exposition stays with the petra JOSS/Computers & Geosciences paper (Omnibus B6). AGU/Goldschmidt abstract early (Thermo2026-adjacent sessions) to flag the territory.

---

*Sources consulted are linked inline; primary anchors: [Sletten & Onstott 1998](https://www.sciencedirect.com/science/article/pii/S0016703797003232), [Nteme et al. 2022](https://www.sciencedirect.com/science/article/abs/pii/S0016703722002198), [Nteme et al. 2023](https://www.sciencedirect.com/science/article/abs/pii/S0016703722006512), [Nteme et al. 2024](https://www.sciencedirect.com/science/article/abs/pii/S0016703724003454), [Harrison et al. 2009](http://sims.ess.ucla.edu/PDF/Harrison_et_al_2009.pdf), [Lee 1995](https://link.springer.com/article/10.1007/BF00311008), [Kramar et al. 2003](https://www.lyellcollection.org/doi/abs/10.1144/gsl.sp.2003.220.01.15), [Villa series part 4 2026](https://doi.org/10.3390/geosciences16080298), [Lister & Forster 2024](https://www.tandfonline.com/doi/full/10.1080/08120099.2024.2421569).*
