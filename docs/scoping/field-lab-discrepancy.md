# Attacking the Field–Lab Discrepancy with Defect-Resolved, Surface-Aging KMC

**Scoping document — 2026-08-18**
Toolchain: `petra` (declarative lattice KMC, Rust; defects + strain fields + multilayer landed) + `quarry` (GPU QM barrier pipeline, validated against Xiao & Lasaga 1994 within 2 kcal/mol). One researcher, one agent, one RTX 4090.

---

## 1. The question and why it matters

Silicate minerals dissolve **2–5 orders of magnitude slower in the field than in the laboratory** when rates are normalized the same way. This is one of the oldest unsolved quantitative problems in geochemistry — canonically framed by [White & Brantley (2003, *Chemical Geology* 202:479–506)](https://www.sciencedirect.com/science/article/abs/pii/S0009254103002560), who showed the gap is not a constant offset but a *decline of rate with the age of the reacting surface*: fresh Panola granite in laboratory columns declined parabolically over 6 years to 7.0×10⁻¹⁴ mol m⁻² s⁻¹, while naturally pre-weathered granite from the same outcrop ran at 2.1×10⁻¹⁵ — and extrapolating the fresh-material decline implied *thousands of years* to reach the field rate. Something about the surface itself keeps changing on timescales no flask experiment reaches.

The working hypothesis of this project: **a substantial, quantifiable fraction of the discrepancy is a lattice-state phenomenon** — depletion of the reactive-site inventory (kinks, steps, strained material at dislocation outcrops, high-energy sites from grinding damage), annealing of surface topography toward low-energy configurations, and the consequent drift of the site-type distribution toward less-reactive populations. Continuum rate laws (`rate = k·A·f(ΔG)`) structurally cannot represent this: they carry one `k` and one `A` where the real system carries a *distribution* over sites that evolves under its own dynamics. A kinetic Monte Carlo model tracks exactly this state, natively.

**Why now, and why the stakes are larger than a legacy puzzle.** Enhanced rock weathering (ERW) — spreading crushed silicate on cropland as CO₂ removal — is a gigatonne-scale industry being built directly on top of this unsolved problem. Published CDR estimates from ERW field trials [span four orders of magnitude](https://www.frontiersin.org/journals/climate/articles/10.3389/fclim.2025.1657058/full); a [three-year multi-site trial (ES&T 2025)](https://pubs.acs.org/doi/10.1021/acs.est.5c09820) found a sustained weathering signal but *limited CO₂ removal*; basalt in the field dissolves far slower than lab kinetics predict, exactly the classic discrepancy replayed with carbon credits attached. Freshly *ground* feedstock is the maximally lab-like surface — maximum grinding damage, maximum reactive-site inventory — so the aging trajectory from "fresh lab surface" to "field surface" is precisely the trajectory an ERW deployment lives through in its first years. [Brantley's 2025 *Reviews of Geophysics* review](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2025RG000881) ("Understanding the Lab-Field Discrepancy in Mineral Dissolution From Flasks to Enhanced Rock Weathering") makes this connection the framing of the field's current agenda, and lists as its *first* unresolved category "changing mineral surface attributes difficult to model." That is a standing invitation for exactly this project.

---

## 2. Literature map

### 2.1 Established

**The phenomenon.** White & Brantley 2003 (above; the PDF circulates freely, e.g. [here](https://remineralize.org/wp-content/uploads/2022/03/1-s2.0-S0009254103002560-main.pdf)): lab rates up to 5 orders above field; rate declines with surface age in *both* lab columns and soil chronosequences; BET surface area *increases* while reactivity *decreases* (their "time-dependent roughness factor" — an internal contradiction for any model where rate ∝ area). Extended by [White et al.'s 13.8-year flow-through columns on four granitoids](https://www.sciencedirect.com/science/article/abs/pii/S0016703716306950) (Merced, Panola, Loch Vale, Rio Icacos) — the longest weathering experiment ever run, and still declining.

**The extrinsic (transport/thermodynamic) explanations.** [Maher (2010, EPSL)](https://www.sciencedirect.com/science/article/abs/pii/S0012821X10001810): weathering rates depend most strongly on fluid residence time; [Maher & Chamberlain (2014, *Science*)](https://www.science.org/doi/10.1126/science.1250770): hydrologic regulation with a thermodynamic ceiling on flux. In this view the field is slow because pore fluids sit near equilibrium — the mineral surface is largely irrelevant. This is the strongest competing hypothesis family and is *partially* correct by consensus; the open question is the attribution split.

**The intrinsic (surface-state) explanations.** Three strands:

1. *Defect-driven dissolution theory*: Lasaga & Blum (1986) on dissolution at dislocations; Cabrera–Levine critical undersaturation for etch-pit opening; [Lasaga & Lüttge (2001, *Science* 291:2400)](https://www.ncbi.nlm.nih.gov/pubmed/11264534) dissolution **stepwave model** — steps emanating from etch pits at dislocation outcrops carry the dissolution flux, giving a mechanistic ΔG dependence far richer than TST rate laws.
2. *Rate spectra*: [Fischer, Arvidson & Lüttge (2012, GCA)](https://www.sciencedirect.com/science/article/abs/pii/S001670371830108X) — the measured "rate" of a surface is a **distribution** (rate spectrum) over surface sites, resolvable by vertical scanning interferometry (VSI); surface-area normalization *obscures* the variance that carries the mechanism. Follow-ons: [holographic-microscopy calcite rate spectra (2017, open access)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5590660/), [statistical analysis frameworks (Minerals 2019)](https://doi.org/10.3390/min9080458), ["pulsating dissolution" (GCA 2024)](https://www.sciencedirect.com/science/article/pii/S0016703724000796) showing intrinsic temporal rate oscillation even at fixed conditions.
3. *Passivating altered layers*: Hellmann's [interfacial dissolution–reprecipitation](https://www.sciencedirect.com/science/article/abs/pii/S0009254111004712) and Daval's [amorphous silica layer passivation](https://www.sciencedirect.com/science/article/abs/pii/S0009254111000982) — nm-scale coatings that restrict fluid access. Established as real for olivine/glass; contested in generality ([Daval: passivation is "not universal"](https://www.sciencedirect.com/science/article/abs/pii/S0009254116304508)).

**The synthesis experiment.** [Gruber, Zhu, Georg, Zakon & Ganor (2014, GCA 147:90–106)](https://www.sciencedirect.com/science/article/abs/pii/S0016703714006218): using Si-isotope-doped solutions on *naturally aged* albite near equilibrium, they argue the gap closes when you combine reduced surface reactivity of aged material with near-equilibrium ΔG — i.e., **both camps are right and the discrepancy is a product of factors**. The Zhu lab's [Navajo Sandstone aquifer work](https://hydrogeochem.earth.indiana.edu/publications/index.html) gives in-situ albite rates ~10⁻¹⁵ mol m⁻² s⁻¹ and matching lab near-equilibrium unidirectional rates via isotope doping.

**The key negative constraint.** [Blum, Yund & Lasaga (1990, GCA)](https://www.osti.gov/biblio/6935823): raising quartz dislocation density from <10⁵ to 5×10¹⁰ cm⁻² left *far-from-equilibrium* dissolution rates **indistinguishable**, despite dramatic etch pitting. Any defect-based story must reproduce this: defect inventory controls rates *near equilibrium* (where only strained material clears the Cabrera–Levine threshold) and controls the *rate spectrum tail*, not the far-from-equilibrium mean. This constraint is a feature for us — it is a sharp, falsifiable two-regime prediction a KMC with strain fields either reproduces or doesn't.

### 2.2 Who is active now (2015–2026)

- **Susan Brantley** (Penn State) — the 2025 RoG review; empirical/field synthesis; ERW framing.
- **Inna Kurganskaya** (Bremen, formerly with Lüttge) — *the* active KMC person in this space: [silicate KMC parametrization (2013)](https://www.researchgate.net/publication/257924556_Kinetic_Monte_Carlo_Simulations_of_Silicate_Dissolution_Model_Complexity_and_Parametrization), [muscovite stochastic model](https://www.researchgate.net/publication/257292502), [Kurganskaya & Lüttge 2021 "Pathways to Equilibrium" (ACS ESC)](https://pubs.acs.org/doi/10.1021/acsearthspacechem.1c00017) — generic rate laws from statistical mechanics with explicit defect-type/amount dependence, [pH/protonation-resolved clay-nanoparticle KMC (Minerals 2024)](https://doi.org/10.3390/min14090900).
- **Cornelius Fischer** (HZDR) & **Andreas Lüttge** (Bremen) — rate spectra, VSI, [upscaling frameworks](https://www.sciencedirect.com/science/article/abs/pii/S000925411830576X), pulsating dissolution.
- **Chen Zhu** (Indiana) — near-equilibrium isotope-doping kinetics, field rates.
- **Damien Daval** (CNRS), **Roland Hellmann** (Grenoble) — altered surface layers.
- **KMC/DFT pipelines**: [Martin et al. 2021 quartz KMC](https://pubs.acs.org/doi/10.1021/acsearthspacechem.0c00303) jointly reproducing rate *and* activation energy from DFT barriers, with the open-source [KIMERA](https://www.researchgate.net/publication/345192590_KIMERA_A_Kinetic_Montecarlo_Code_for_Mineral_Dissolution) code; Izadifar (periodic DFT barriers); Schliemann & Churakov, Liu & Ruiz Pestana (AIMD barriers) — see `qm/SURVEY.md` §7.
- **ERW MRV community** — [Frontiers in Climate](https://www.frontiersin.org/journals/climate/articles/10.3389/fclim.2025.1657058/full), [ES&T field trials](https://pubs.acs.org/doi/10.1021/acs.est.5c09820), integrated frameworks like [T&C-SMEW](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12710202/).

**Current consensus**: multi-causal; the 2025 review's five categories are (1) changing surface attributes, (2) hydrology/equilibration/armoring, (3) spatial heterogeneity, (4) biology, (5) system feedbacks. **No quantitative attribution exists** — nobody can say "surface aging buys N of the 2–4 log units." That attribution is the prize.

---

## 3. The gap: what hasn't been done that this toolchain can do

Precisely stated: **no one has run a long-duration, defect-inventory-resolved KMC in which the surface site-type distribution evolves freely from a "fresh, damaged" initial state toward an "aged" state, with ab-initio site-resolved barriers, emitting the observables the experimental record is written in (rate-vs-time, rate spectra, BET-vs-geometric area), to put a number on the intrinsic-aging contribution to the lab–field discrepancy.**

What exists, and why each stops short:

| Prior work | What it did | Why it isn't this project |
|---|---|---|
| Kurganskaya & Lüttge 2021 | Analytic rate laws + KMC with defect types/amounts | Steady-state kinetic pathways, not multi-year transient aging; no depletion of a finite defect inventory; no comparison to the White/Brantley decline curves |
| Lasaga & Lüttge 2001 stepwave | ΔG-dependence from etch-pit step trains | Continuum-ized immediately; static dislocation source; never run to source exhaustion |
| Martin et al. 2021 / KIMERA | DFT→KMC, reproduces quartz rate + Ea | Far-from-equilibrium steady state; no defects, no aging, no rate spectra |
| Fischer/Lüttge rate spectra | Measured the site-rate distribution | Experimental concept; the KMC that *predicts* spectrum evolution under aging doesn't exist |
| White/Brantley columns | The 13.8-year decline data | Purely empirical; "deep etching of dislocations and grain boundaries" invoked, never simulated |
| Gruber/Ganor 2014 | Showed aged-surface × ΔG closes the gap | A two-factor bookkeeping argument; no mechanistic model of *why* aged surfaces are 10–100× less reactive |

Why *this* toolchain is unusually positioned:

1. **petra already has the hard parts nobody else's KMC has together**: declarative decks (any mineral, no recompile), crystallographic defects with analytic strain-energy fields coupled to barriers via BEP transfer coefficients (STRAIN.md §2 — the Blum–Lasaga β=1 treatment is a one-line deck choice), finite defect populations via `[[init]]` (probability or explicit sites), etch-pit gates already green, multilayer structure, Δμ/activity solution coupling for near-equilibrium runs, trajectory record/replay, and incremental O(1)-per-step updates that make big-lattice long runs feasible.
2. **quarry closes the parameter-provenance hole** that killed credibility for 1990s-era KMC: site-resolved, connectivity- and protonation-resolved barriers at ~1 GPU-hour each, already validated against Xiao & Lasaga within 2 kcal/mol, emitted directly as petra deck fragments with provenance.
3. The combination — **QM-grounded barriers + defect strain + finite inventories + years of simulated time + rate-spectrum output** — is the experiment the 2025 review implicitly asks for and no group has assembled.

---

## 4. Testable hypotheses (ranked)

**H1 — Inventory depletion reproduces the parabolic decline.** Starting from a "ground fresh" surface (high step/kink density + damage layer + dislocation outcrops), gross dissolution rate declines quasi-parabolically over simulated years, by 1–2 orders of magnitude, as the high-energy site inventory is consumed faster than lattice-plane retreat regenerates it.
*Observable*: rate vs. time master curve → White & Brantley 2003 Fig-format comparison; White et al. 2017 column data.
*Needs*: long-time runs (P6 performance + stiffness handling), damaged-surface init recipes (exists: `[[init]]` passes), rate time-series output (exists).

**H2 — The rate spectrum narrows and down-shifts with age.** Fresh surfaces show a broad spectrum with a fat fast tail (strained/kink sites); aged surfaces show a narrow, slow spectrum. The *mean* moves 1–2 orders while the *mode* barely moves — explaining why "the same mineral" gives irreproducible lab rates (the 2–3 order lab-to-lab variance noted by Kurganskaya & Lüttge 2021).
*Observable*: per-site instantaneous rate histograms vs. VSI/holography rate spectra (Fischer's calcite data on PANGAEA; feldspar VSI in the literature).
*Needs*: **rate-spectra exporter** (new, small: the engine already holds per-site total rates in the Fenwick leaves — dump + histogram).

**H3 — Two-regime defect dependence (the Blum–Yund–Lasaga test).** Far from equilibrium, total rate is insensitive to dislocation density (10⁵→10¹⁰ cm⁻² changes rate <2×) because step retreat everywhere dominates; near equilibrium (|ΔG| below the Cabrera–Levine threshold), rate is *proportional* to surviving defect outcrop density, and **depletes as pits consume their cores or age anneals the strained inventory** — the field regime. This single mechanism yields both the null 1990 result and multi-order field slowdown.
*Observable*: rate vs. ΔG curves at several dislocation densities and ages; hysteresis between fresh and aged surfaces at fixed ΔG → Gruber/Ganor 2014 and Zhu near-equilibrium data.
*Needs*: strain fields (done), Δμ sweeps (done), pit-source exhaustion over long runs (needs H1 machinery).

**H4 — The roughness paradox falls out.** BET-proxy area (accessible site count, including pit walls and micropores) *rises* during aging while reactivity falls, reproducing White & Brantley's observation that BET-normalized field rates look even slower — because new area created by etching is predominantly *low-energy* area.
*Observable*: geometric vs. accessible surface-area time series; roughness factor vs. surface age → White & Brantley's time-dependent roughness compilation.
*Needs*: **surface-area observables** (new: face-counting for geometric, probe-accessible site census for BET-proxy — both computable from the lattice graph).

**H5 — Annealing moves accelerate aging beyond pure depletion.** Adding reattachment/surface-diffusion reactions (Ostwald-ripening-like moves: detach from high-strain, reattach at low-energy sites — expressible today as paired rewrites with detailed-balance rates) speeds the decline and makes aged morphology smoother, matching the coarsening seen in long-term experiments; without them, decline plateaus too early.
*Observable*: decline exponent and terminal morphology with/without annealing channel.
*Needs*: deck authoring only (mechanism exists); detailed-balance linter (P2 remainder) becomes important here.

**H6 (stretch, ERW-facing) — Grinding resets the clock, with a computable half-life.** A comminution event (re-damaged init state) restores the fast spectrum; the model yields an "aging half-life" per mineral/condition — the quantity ERW MRV actually needs (how long does fresh-feedstock kinetics last?).
*Observable*: rate recovery factor and decay time vs. particle damage parameterization; qualitative comparison to ERW trial time series.

---

## 5. Calibration / validation datasets (all accessible to an independent researcher)

1. **White & Brantley 2003 compilation** — rate-vs-duration for lab + field silicates (their Fig. 1 master plot), Panola column series, roughness factors. Paper PDF is [openly mirrored](https://remineralize.org/wp-content/uploads/2022/03/1-s2.0-S0009254103002560-main.pdf); data are in tables.
2. **White et al. 2017 GCA 13.8-year columns** ([link](https://www.sciencedirect.com/science/article/abs/pii/S0016703716306950)) — four granitoids, fresh vs. weathered pairs; supplement has effluent time series. The single best transient target.
3. **Merced chronosequence** ([White et al. 1996, GCA](https://www.sciencedirect.com/science/article/abs/pii/0016703796001068)) — 10 ka–3 Ma soils: mineralogy, BET area, and primary-silicate rates vs. age. The long-time anchor.
4. **Fischer et al. calcite rate-spectra data on PANGAEA** ([doi:10.1594/PANGAEA.829556](https://doi.pangaea.de/10.1594/PANGAEA.829556)) — open download; plus [open-access holographic rate spectra](https://pmc.ncbi.nlm.nih.gov/articles/PMC5590660/). Calcite, not silicate — but the *spectrum-shape* validation target for the H2 observable.
5. **Kaolinite lab kinetics** — Ganor et al. 1995, Cama et al. 2002 (acid, Ea ≈ 29 kJ/mol), Bauer & Berger 1998 (alkaline): the absolute-rate targets for the quarry-parameterized kaolinite deck (already catalogued in `qm/SURVEY.md` §7.3).
6. **Gruber et al. 2014 aged-albite** ([GCA 147:90](https://www.sciencedirect.com/science/article/abs/pii/S0016703714006218)) and **Zhu-lab near-equilibrium Si-isotope series** ([publications](https://hydrogeochem.earth.indiana.edu/publications/index.html)) — the H3 near-equilibrium targets.
7. **Blum, Yund & Lasaga 1990 quartz dislocation experiments** ([OSTI](https://www.osti.gov/biblio/6935823)) — the falsification constraint.
8. **ERW field-trial supplements** ([ES&T 2025](https://pubs.acs.org/doi/10.1021/acs.est.5c09820); [Frontiers 2025](https://www.frontiersin.org/journals/climate/articles/10.3389/fclim.2025.1657058/full)) — qualitative H6 context.

Access note: most primary papers are paywalled at Elsevier but obtainable (author sites, preprint mirrors, interlibrary channels); the *data* needed are largely in-paper tables/figures at digitizable fidelity, plus the open PANGAEA sets.

---

## 6. Phased project plan (one person + agent + one GPU)

**Phase 0 — Observables (2–4 weeks; petra work only).**
Build the three missing outputs on existing decks (Kossel + kaolinite golden): (a) rate-spectrum exporter (per-site rate histogram time series — data already in the Fenwick leaves), (b) geometric + BET-proxy surface-area observables from the lattice graph, (c) a "surface age" bookkeeping column (time since site exposure). Exit test: spectra/area on the existing `kossel-etchpit` demo behave sanely (pit growth fattens the fast tail; area rises).
*Exists*: engine, defects/strain, etch-pit gate, replay/WASM viz (which makes the aging movies nearly free). *New*: the three observables.

**Phase 1 — Generic aging study (1–2 months; CPU-bound ensembles).**
Kossel-family decks with finite defect inventories (`[[defects]]` lines + `[[init]]` damage layers), with and without annealing moves (H5 decks), far from and near equilibrium. Produce: decline master curves (H1), spectrum evolution (H2), the two-regime defect-density map (H3), roughness paradox (H4). This phase is publishable alone as a mechanism paper: *"Surface aging and the intrinsic component of the lab–field discrepancy: a defect-resolved KMC study."* Needs P6-style performance attention (10⁶-site lattices, stiffness detection for fast reversible pairs) — the one real engineering risk (see §7).

**Phase 2 — QM-grounded kaolinite (2–4 months; GPU barrier campaign in parallel).**
quarry campaign: connectivity- and protonation-resolved barriers for the five site families (~30–60 barriers × ~1 GPU-hr — weeks of unattended 4090 time), emitted as deck `by_count` tables replacing the flat legacy numbers. Validate absolute rate + Ea against Ganor/Cama/Bauer-Berger. Then rerun Phase 1 protocols on the real mineral. The multilayer deck + stacking-fault areal energies (STRAIN.md §3.2, open question 1) let kaolinite's actual defect style — stacking disorder — be tested against line-defect fields.

**Phase 3 — The discrepancy budget (1–2 months).**
Near-equilibrium, aged-surface, defect-depleted runs → a stacked attribution: log-units from (i) site-inventory depletion, (ii) defect exhaustion at field ΔG, (iii) annealing, (iv) normalization artifact (H4) — versus the log-units the extrinsic literature assigns to transport/equilibrium. Deliverable: *"How much of the field–lab discrepancy is the lattice? A quantitative ceiling from defect-resolved KMC."* Compare against Gruber/Ganor and Zhu targets.

**Phase 4 — ERW translation (optional, 1 month).**
H6: aging half-life of ground feedstock, comminution reset, implications for rate-based MRV. Short, applied, high-visibility.

Feature ledger: **exists** — engine, decks, strain/defects, init substitution, Δμ, multilayer, replay/viz, quarry pipeline + emitters. **To build** — rate-spectra output, area observables, exposure-age tracking, long-run performance/superbasin handling (P6, the only substantial engine work), annealing decks, stacking-fault re-registry (only if Phase 2 demands it), quarry barrier campaign itself.

---

## 7. Risks and honest novelty assessment

**Timescale bridging is the hard technical risk.** Years of simulated time with μs-scale fast events is the classic KMC stiffness wall. Mitigations: near-equilibrium conditions suppress raw event rates; superbasin/τ-leaping is a known (roadmap) technique; and the *claim* can be formulated on decline *exponents and mechanisms* extrapolated over decades rather than brute-forced. If 13.8 simulated years is unreachable, 13.8 simulated *days* with correct scaling behavior plus the analytic envelope may still carry the argument — but a reviewer will press here, and the scoping should budget real P6 work.

**A skeptical reviewer will say:**
- *"Maher showed the field is transport/equilibrium-controlled; your intrinsic effect is second-order."* Answer: we do not claim exclusivity; we compute the intrinsic **ceiling** and note Gruber/Ganor showed the factors multiply. The near-equilibrium runs put both effects in one model.
- *"Blum & Lasaga 1990 killed defect-controlled dissolution."* Answer: only far from equilibrium — and reproducing that null result is our H3 exit test, turning the objection into a validation gate.
- *"Aging in real systems is dominated by altered/passivating layers and secondary coatings (Hellmann, Daval), which a lattice KMC cannot represent."* This is the **honest structural limit**: petra has no amorphous interphase, no coupled reprecipitation film. Response: state the scope — we quantify the *lattice-state* channel; where our budget falls short of the observed decline, the residual is an *estimate of the coating/transport channels' share*. That reframing makes the limitation a result.
- *"Your barriers are cluster-QM with known systematics (lattice resistance +20–70 kJ/mol)."* Mitigation: quarry's calibration ladder (DLPNO-CCSD(T) spot checks, comparison to Izadifar's periodic values, and the Ganor/Cama absolute-rate anchor); sensitivity analysis over barrier uncertainty is cheap in petra (deck sweeps).
- *"Kaolinite is a poor test mineral — basal-inert, slow, and the classic discrepancy datasets are feldspathic."* Fair. Hence the two-track design: the generic (Kossel/feldspar-like) mechanism papers carry the discrepancy argument; kaolinite is the QM-grounded absolute-rate showcase. A feldspar deck is a natural follow-on — petra is chemistry-free by construction.

**Novelty verdict**: the mechanism-level ingredients are all published separately (stepwave, rate spectra, DFT→KMC, aging experiments); the *assembly* — transient aging KMC with finite defect inventories, QM-provenance barriers, and experiment-native observables, aimed at a quantitative attribution — is not in the literature as of this scoping (checked against the 2025 RoG review and 2020s KMC output of Kurganskaya, Fischer, Martin/KIMERA groups). It would survive review as a modeling contribution if the timescale and validation gates hold; the attribution claim will be contested, which is what makes it citable.

**Independent-researcher risk**: no lab, no VSI instrument — all validation is against published data; this is workable (the datasets exist, §5) but means no bespoke experiment can rescue a gap. Collaboration outreach (Kurganskaya, Fischer, Zhu groups) after the Phase 1 preprint is the mitigation and likely the route to a stronger Phase 3 paper.

---

## 8. Publication and venue targets

| Deliverable | Venue (1st / alt) | Timing |
|---|---|---|
| Phase 1 mechanism paper (generic aging KMC) | *ACS Earth & Space Chemistry* / *Geochimica et Cosmochimica Acta* | after Phase 1 |
| petra tool paper (declarative KMC engine) | *Computers & Geosciences* / *JOSS* | anytime; boosts credibility |
| Phase 2 kaolinite QM→KMC (quarry + deck) | *GCA* / *Chemical Geology* | after Phase 2 |
| Phase 3 attribution ("how much is the lattice") | *Chemical Geology* (the White & Brantley home turf) / *EPSL* | capstone |
| Phase 4 ERW aging half-life note | *Environmental Science & Technology* / *Frontiers in Climate* | opportunistic |
| Conference presence | Goldschmidt, AGU (Brantley/ERW sessions) | Phase 1 onward |

The in-browser WASM visualization is an underrated asset for all of these: reviewable, shareable, live simulations of etch pits aging a surface are a communication weapon no competing group offers.
