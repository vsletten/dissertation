# QM Platform Survey

**Rebuilding the dissertation's quantum chemistry on a 2026 workstation**

*Surveyed 2026-08-12. Target machine: AMD Ryzen 9 7950X (16C/32T), 64 GB DDR5,
NVIDIA RTX 4090 24 GB, Linux. Constraint: free software only — open source
preferred, free-to-use closed source acceptable with flags, paid options
excluded.*

The original goal, restated: bridge ab initio quantum mechanics to the kinetic
Monte Carlo model — activation energies and rate constants for elementary
hydrolysis reactions (Si–O–Si, Si–O–Al, Al–OH–Al bridges) from transition
state theory via electronic structure calculations on small aluminosilicate
clusters. In the late 1990s that meant Gaussian DFT on a DEC Alpha. This
survey answers: what do we run today, on this box, for free?

---

## TL;DR — the recommended stack

**We do not need to roll our own QM engine.** The physics is well known, but a
production DFT code (integrals, SCF, analytic gradients, analytic Hessians,
solvation) is a decade-scale effort, and the open-source ecosystem has already
converged on excellent answers. The roll-our-own part is the thin
*orchestration layer* — the code that maps the KMC model's reaction table to
barrier calculations. That's the platform we build.

| Layer | Choice | Why |
|---|---|---|
| Workhorse DFT engine | **PySCF + GPU4PySCF** | Apache-2.0, Python-native, analytic gradients *and* Hessians on the 4090, SMD/PCM solvation, D3/D4 dispersion. The only free stack that makes the GPU a real asset. |
| Cross-check + benchmark engine | **ORCA 6.1.1** | Free for academic/personal use (closed source). Turnkey OptTS / NEB-TS / IRC / analytic frequencies, and the only free access to DLPNO-CCSD(T) for calibrating DFT barriers. |
| TS search / IRC | **geomeTRIC** (drives PySCF natively) + **Sella** (over ASE) | Both maintained, both do saddle points and IRC, backend-agnostic. |
| Pre-screening | **xtb (GFN2-xTB) + CREST**; optionally an OMol25-class ML potential (**UMA** or **OrbMol**) | Seconds-per-structure geometry pre-optimization, conformer/water-placement sampling, TS guesses. |
| Thermochemistry → k(T) | **Shermo or GoodVibes** (quasi-RRHO ΔG) + **Arkane** (TST with Wigner/Eckart tunneling) or a small in-house Eyring module | All free, all scriptable. |
| Glue | **ASE + plain Python** (optionally quacc) | Everything above speaks ASE or Python; heavier workflow engines (AiiDA etc.) are not worth it on one box. |
| Future periodic/slab work | **CP2K** | GPL, NEB + dimer TS search, the natural bridge from clusters to kaolinite slabs. |

Realistic hardware expectation up front: the RTX 4090's double-precision units
are deliberately capped at 1/64 of FP32 (~1.3 TFLOPS FP64 — roughly equal to
the 7950X's peak). GPU4PySCF runs in FP64, so expect **~2–5× over the 16-core
CPU for SCF/gradients and ~5–10× for Hessians/frequency jobs** — a big
practical win for the frequency-heavy TST workflow, but not the 20–50× seen on
datacenter cards. NVIDIA's new cuEST library (FP64 emulation on INT8 tensor
cores, already integrated into GPU4PySCF) may claw back part of the gap; see §3.

---

## 1. Then vs now

The late-1990s reality: HF or early DFT with 3-21G/6-31G* on 10–30-atom
clusters, MP2 single points if you were lucky, days of wall time per job on a
~1 GFLOPS DEC Alpha with a few hundred MB of RAM. Xiao & Lasaga's benchmark
disilicate clusters were ~15–25 atoms.

This workstation is roughly **a thousand times** that machine in raw FLOPS
(CPU alone), with 64 GB of RAM, plus a GPU adding another multiple for the
right workloads. Concretely, on one node today:

- Every calculation in the original dissertation could be re-run, at a much
  higher level of theory, in an afternoon.
- Hybrid-DFT (e.g. ωB97X-D4 or B3LYP-D4/def2-TZVP) geometry optimization +
  analytic frequencies on **50–150-atom clusters**: overnight on the CPU,
  hours with the GPU. Composite-method (r²SCAN-3c) opt+freq on 100–300 atoms:
  hours.
- **DLPNO-CCSD(T)/CBS single points on 30–80-atom clusters** (ORCA): hours to
  a day — the "gold standard" calibration layer that was unaffordable in 1999.
- Whole new modalities fit on one box: ab initio MD with metadynamics on
  ~200-atom explicit-water systems (CP2K), and ML-potential MD on thousands of
  atoms at near-DFT quality on the 4090 — the biggest qualitative capability
  that did not exist then.

Cluster size at fixed method grew ~10–30×; method quality moved from
HF/MP2-in-small-basis to dispersion-corrected range-separated hybrids with
coupled-cluster calibration; and the field now computes *site-resolved barrier
ladders* (by connectivity and protonation state) rather than single numbers —
which is exactly what a KMC reaction table wants.

---

## 2. The DFT engines

### 2.1 At a glance

| Package | License | TS opt | Analytic DFT Hessian | IRC | Solvation | Python | GPU | Verdict |
|---|---|---|---|---|---|---|---|---|
| **PySCF** 2.14 | Apache-2.0 | via geomeTRIC/Sella | **yes** | via Sella/geomeTRIC | PCM, **SMD**, ddCOSMO | native | **GPU4PySCF** | **Primary** |
| **ORCA** 6.1.1 | free academic/personal, closed | **native OptTS, NEB-TS** | **yes** | **native** | CPCM, **SMD** | no (wrappers) | no | **Cross-check + DLPNO-CCSD(T)** |
| **Psi4** 1.11 | LGPL-3.0 | native (OptKing) | no (finite-diff) | native | PCM, ddCOSMO; **no SMD** | native | no | Easiest install; freq-heavy work pays 2×3N gradients per Hessian |
| **NWChem** 7.3.1 | ECL-2.0 | native (`saddle`) | yes (in-core, small systems) | yes | COSMO, SMD | weak | not for DFT | Open-source fallback; unique Polyrate/VTST interface; archaic ergonomics |
| **GAMESS-US** 2024R2 | free, no-redistribution | native | HF only | native | PCM, SMD | no | partial | Dominated by the above; license blocks containers |
| **xtb** 6.7.1 | LGPL-3.0 | (screening level) | numerical | — | ALPB/GBSA | via tblite/ASE | no | Pre-screening layer, not for barriers |

Skipped after review: **OpenMolcas** (multireference specialist — keep in the
back pocket for CASPT2 sanity checks on a suspicious TS), **Dalton** (dormant
since 2020), **deMon2k** (gated, stale), **DFTB+** (redundant with xtb here),
**BAGEL/CFOUR/MRCC** (wrong specialty). Commercial, excluded: Gaussian,
Q-Chem, Molpro, TURBOMOLE, ADF/AMS, Jaguar, TeraChem.

### 2.2 PySCF — the primary engine

- Apache-2.0, `pip install pyscf`; v2.14.0 (July 2026), very active
  (10-year retrospective paper in 2026). https://github.com/pyscf/pyscf
- Full Kohn–Sham DFT via libxc (LDA → range-separated hybrids incl.
  ωB97M-V); D3/D4 dispersion via `pyscf[dispersion]`; def2 basis family
  shipped; density fitting (RI-J/RI-JK) throughout.
- **Analytic Hessians** for restricted/unrestricted HF and DFT, with a
  built-in thermochemistry module (ZPE, H, S, G) — the core of the TST
  workflow. Solvation: IEF-PCM/C-PCM **and SMD**, with analytic gradients and
  semi-analytic Hessians. https://pyscf.org/user/solvent.html
- Geometry optimization drives **geomeTRIC** natively; TS optimization comes
  from geomeTRIC's `transition=True` mode, the `qsdopt` extension, or Sella
  over ASE — assembled, not turnkey (this is the one ergonomic gap vs ORCA).
- OpenMP-threaded, single-node; 64 GB RAM is comfortable for DF hybrid DFT at
  def2-TZVP up to ~150 atoms.
- Being a library rather than a program is the point: the KMC-facing
  orchestration layer can call it directly, in-process, with no input-file
  generation or output parsing.

### 2.3 ORCA 6.1.1 — the turnkey cross-check

- Closed source, **free for academic and personal use**; commercial use
  requires a paid FACCTs license; download via registration. No
  redistribution (so no containers/conda). v6.1.1 December 2025.
  https://www.faccts.de/orca/
- This is the closest modern equivalent of the dissertation-era Gaussian
  workflow, and it is *built* for exactly our pipeline: `Opt` → `OptTS`
  (eigenvector following) / `NEB-TS` → analytic frequencies → `IRC` →
  quasi-RRHO thermochemistry, with CPCM/SMD, RIJCOSX-accelerated hybrids, and
  native GFN2-xTB pre-optimization. MPI-parallel; 16 cores gives ~14–15× on
  typical DFT jobs.
- Uniquely valuable here: **DLPNO-CCSD(T)** — near-CCSD(T) accuracy on
  30–80-atom clusters on this hardware. No other free code offers this. It is
  the calibration standard for DFT barrier heights in the modern literature
  (see §6.4).
- Use it for: validating PySCF barriers, NEB-TS on awkward reactions, and all
  coupled-cluster single points.

### 2.4 Psi4 — honorable mention

LGPL, one-line conda install, native TS+IRC via OptKing, clean Python API.
Two dealbreakers for this project: DFT frequencies are finite-difference
(a ~100-atom Hessian costs ~600 gradient evaluations), and there is no SMD
solvation (a long-standing gap). Keep as a third opinion; don't build on it.

### 2.5 NWChem — the all-open fallback

Fully open (ECL-2.0), analytic DFT Hessians (in-core — small systems only),
COSMO+SMD, `task dft saddle`, and one unique asset: the **DIRDYVTST module
interfacing Truhlar's Polyrate** for variational TST with multidimensional
tunneling — the most rigorous free rate-constant machinery anywhere. Costs:
1990s input-deck ergonomics, manual memory tuning, and slower hybrid DFT (no
COSX-class acceleration). Worth keeping installed for the VTST option alone.

### 2.6 The semi-empirical layer: xtb / tblite / CREST

GFN2-xTB (LGPL, Grimme lab) handles Al/Si/O/H fine at the *screening* level:
pre-optimize clusters in seconds, sample conformers and explicit-water
placements with CREST, generate initial Hessians and TS guesses — then refine
everything at DFT. Barrier heights at this level are **not** trustworthy;
that's not its job. `tblite` (v0.7.0, 2026) is the maintained library
implementation with Python/ASE APIs. Watch **g-xTB** — Grimme's next-gen
method approximating ωB97M-V quality, landing in tblite — as a likely future
upgrade of this layer.

### 2.7 Periodic codes, for later

When cluster work graduates to kaolinite slabs/edges: **CP2K** (GPL,
Gaussian-and-plane-waves, NEB + dimer TS search, hybrid functionals via ADMM,
annual releases, v2026.1) is the natural home — same Gaussian-basis mindset,
huge overlap with the mineral-surface literature (§7). Quantum ESPRESSO,
ABINIT, GPAW, SIESTA are all healthy but are plane-wave/PAW-shaped and less
convenient for this chemistry; FHI-aims is licensed (€2k academic) and
excluded.

---

## 3. The GPU question: what the RTX 4090 is actually worth

### 3.1 The FP64 reality

The 4090 (AD102): ~82.6 TFLOPS FP32, but FP64 capped at 1/64 → **~1.29 TFLOPS
FP64**, which is *roughly the same as the 7950X's peak* (16 Zen4 cores ≈
1.1–1.2 TFLOPS FP64). Quantum chemistry traditionally demands FP64. What the
card still brings: ~12× the memory bandwidth of dual-channel DDR5 (1008 vs
~83 GB/s), massive latency hiding, and 660+ INT8 tensor TOPS if a code can use
them. Every conclusion below follows from this asymmetry.

### 3.2 GPU4PySCF — the only serious free option, and it covers everything

- Apache-2.0, part of the PySCF org; v1.8.1 (Aug 2026), roughly monthly
  releases; `pip install gpu4pyscf-cuda12x`. Also the engine behind Rowan's
  commercial cloud DFT. https://github.com/pyscf/gpu4pyscf
- Production coverage is exactly our workflow: direct and density-fitted SCF,
  analytic gradients, **analytic Hessians** (HF + DFT incl. hybrids and
  range-separated, restricted and unrestricted), TS search via
  geomeTRIC/pysisyphus interfaces, **PCM and SMD with analytic gradients and
  Hessians**, D3/D4, ECPs, TDDFT. Known gaps: no double hybrids, basis up to
  g functions (fine for Al/Si/O/H), DF module tuned for &lt;~200 atoms.
- Published benchmarks (A100-80G vs Q-Chem on 32 vCPUs, def2-TZVPP): SCF ~20×,
  gradients 20–40×, **Hessians ~50×** (higher with PCM). But GPU4PySCF runs
  in **FP64 throughout**, and the same authors measured the **RTX 4090 at
  5–10× slower than an A100** for these workloads — the FP64 wall.
  (arXiv:2404.09452; arXiv:2605.06489)
- **Chained estimate for this box: ~2–5× over the 16 cores for SCF/gradients,
  ~5–10× for Hessian/frequency jobs**, best case at 100+ atoms with density
  fitting + solvation, worst case near parity on small molecules. Frequency
  jobs are our dominant cost, and that's where the GPU helps most.
- Memory: the 24 GB VRAM spills DF tensors to host RAM asynchronously; in
  practice the **64 GB of system RAM** is the harder ceiling for ~150-atom
  def2-TZVP analytic Hessians.

### 3.3 The mixed-precision escape hatch (watch this space)

The literature is unambiguous that *dynamic* mixed precision preserves
sub-kcal/mol accuracy at 2–7× speedups (Luehr/Martínez JCTC 2011 — TeraChem's
foundation; Kussmann & Ochsenfeld JCP 2021: max error 1.8 µEh at 7× on a
consumer card). Two 2026 developments finally bring this to the open stack:

- **NVIDIA cuEST** (beta, 2026): a CUDA-X Gaussian-basis DFT primitives
  library with **FP64 emulation via INT8 tensor cores** (Ozaki scheme),
  explicitly targeting GPUs with FP32:FP64 ≥ 32 — i.e., exactly this card.
  Already integrated into GPU4PySCF (v1.7.0+) for DF-J/K builds.
  `pip install nvidia-cuest-cu12`. No published 4090 QC benchmark yet —
  worth testing empirically on our own workload.
- Adaptive-precision INT8 density fitting in GPU-PySCF (arXiv:2601.08077,
  2026): **up to ~2–3× over FP64 on an RTX 4090** with unchanged converged
  energies.

(TeraChem — commercial, excluded — remains the proof of what full mixed
precision achieves on a 4090.)

### 3.4 Everything else, briefly

**QUICK** (open source, ships with AmberTools, GPU HF/DFT in FP64): no
Hessians, no range-separated hybrids or meta-GGAs, Cartesian-only basis —
dominated by GPU4PySCF here. **NWChem** GPU = coupled-cluster triples only.
**CP2K/QE/GPAW** GPU paths are FP64-dependent (their own docs warn against
consumer cards) and periodic-shaped. **Psi4** has no native GPU SCF. And the
**ML-potential ecosystem** (§5) runs on the 4090's FP32/tensor cores at full
speed — the card is unambiguously first-class there, whatever the FP64 story.

---

## 4. The tooling layer

### 4.1 Optimizers, TS finders, reaction paths

| Tool | License | Status (mid-2026) | Role here |
|---|---|---|---|
| **geomeTRIC** 1.1 | BSD-3 (+no-AI-training clause) | active | Minimization, constrained scans, **TS optimization, NEB, IRC**; TRIC internal coordinates; PySCF drives it natively. The core optimizer. |
| **Sella** 2.5 | LGPL-3.0 | active (perf work ongoing) | Partitioned-RFO **saddle-point search + IRC** over any ASE calculator — the standard TS refiner for xtb/MLIP-level work; accepts externally supplied Hessians. |
| **autodE** 1.4.5 | MIT | active | **Automated reaction profiles** (conformers → NEB/CI-NEB → TS location → imaginary-mode verification) driving ORCA/xtb. No PySCF backend. |
| **pysisyphus** 1.0 | GPL-3.0 | **unmaintained since Nov 2024** | Broadest single toolbox (growing string, NEB variants, RS-I-RFO, multi-integrator IRC). Use pinned, don't build on it. |
| ASE built-ins | LGPL | active (3.29, Jun 2026) | NEB/CI-NEB/AutoNEB + dimer method; serviceable but bare-bones. |

ML transition-state guessers (React-OT, NewtonNet-based workflows) are
impressive but trained on organic C/H/N/O chemistry — out of domain for
Al/Si; noted as inspiration only.

### 4.2 Frequencies → ΔG‡ → k(T)

- **Shermo 2.6.1** (free, Tian Lu) and **GoodVibes 4.3** (MIT, Paton lab):
  quasi-RRHO Gibbs energies (both Grimme entropy-interpolation and Truhlar
  frequency-raising schemes), frequency scale factors, standard-state
  corrections. Shermo parses ORCA/NWChem/GAMESS output; GoodVibes is pure
  Python.
- **Arkane** (MIT, inside RMG-Py 4.0, actively maintained): the most complete
  open TST engine — canonical TST with **Wigner and asymmetric Eckart
  tunneling**, 1D/hindered-rotor corrections, master-equation pressure
  dependence; parses ORCA and Psi4 output; fully scriptable.
- **Pilgrim** (Fernández-Ramos/Truhlar): free **variational TST with
  small-curvature tunneling** driving Gaussian/ORCA — the free answer if SCT
  becomes scientifically necessary (plausible for proton-shuttle steps in
  hydrolysis). NWChem's DIRDYVTST→Polyrate is the alternative route.
- For the platform's core loop, the Eyring equation from quasi-RRHO ΔG‡ with
  a Wigner/Eckart correction is a ~50-line in-house module — worth owning for
  transparency, with Arkane as the audited reference implementation.
- Dormant/avoid: TAMkin (2019), KiSThelP (GUI-era, no ORCA parser),
  overreact (stale).

### 4.3 Workflow glue

**ASE + plain Python + files/SQLite** is the right substrate for one user on
one box. **quacc** (BSD-3, very active, v1.5.2 Aug 2026) is the one framework
worth considering: ASE-native recipes for exactly our job shapes (TS, IRC,
freq with Sella/NewtonNet/ORCA), zero required services. AiiDA, pyiron, and
QCFractal are healthy projects whose provenance/scheduling machinery pays off
on clusters and multi-thousand-job campaigns, not here — skip.

---

## 5. ML potentials: the layer that didn't exist in 1999

Universal ML interatomic potentials (MLIPs) now provide near-DFT
energies/forces at ~10⁵–10⁶ steps/day on the 4090 — nanosecond reactive
sampling of a solvated cluster or slab as a workstation job. For this
chemistry (Al/Si/O/H, bond breaking, charged species):

| Model | Weights license | Al/Si coverage | Notes |
|---|---|---|---|
| **UMA 1.2** (Meta FairChem) | gated "FAIR Chemistry" (commercial OK, geo-restricted) | yes — `omol` head, ωB97M-V-level OMol25 training, **charge/spin aware** | Arguably the best reactive pre-screener available in 2026. |
| **OrbMol** (Orbital Materials) | **Apache-2.0** | yes — OMol25-trained, charge/spin inputs | The license-clean OMol25-class option; use `-conservative-` variants. |
| **MACE-MP / MPA** | MIT | yes (materials PES) | Good for periodic/mineral checks; soft on water H-bonds; MACE-OMOL is stronger but ASL non-commercial. |
| **MatterSim** (Microsoft) | MIT | yes (materials) | Bulk/mineral phases across T/P; not reaction-specialized. |
| **AIMNet2** | MIT | **no Al — disqualified** | Si-only test cases at most. |
| **NequIP/Allegro** | MIT | train-your-own | The later-phase option: fit a bespoke potential to our accumulated DFT data. |

Two hard-won caveats from the 2025 benchmark literature (ReactBench):
pretrained universal MLIPs **frequently fail to localize transition states
without task-specific fine-tuning**, and non-conservative (direct-force)
variants produce **broken Hessians** — use energy-conserving checkpoints for
any saddle-point work. So the sane division of labor: MLIPs generate paths,
guesses, and conformers; **DFT owns every number that enters the KMC** (TS
refinement, frequencies, barriers). NVIDIA's ALCHEMI toolkit (batched MLIP
MD/relaxation) is relevant if this layer grows.

---

## 6. Methods: what the 2020s literature actually runs

### 6.1 The standard protocol (the Grimme-school consensus)

The de facto best-practice paper (Bursch, Mewes, Hansen, Grimme, *Angew.
Chem.* 2022) prescribes a tiered workflow that maps directly onto our stack:

1. **Screening**: GFN2-xTB / CREST for conformers and water placement.
2. **Geometries + frequencies**: a "3c" composite — **r²SCAN-3c** is the
   workhorse ("Swiss army knife", near hybrid/QZ quality at mGGA cost) — or a
   dispersion-corrected hybrid (PBE0-D4, B3LYP-D4/def2-TZVP). ωB97X-3c is the
   newer hybrid-level composite.
3. **Final energies**: range-separated hybrid — **ωB97M-V** is the
   GMTKN55-best hybrid — at def2-TZVPD or better. Diffuse functions matter
   for anionic pathways (deprotonated silanols, Al(OH)₄⁻).
4. **Always**: dispersion correction (D3(BJ)/D4/VV10); quasi-RRHO treatment
   of low frequencies; implicit solvation on top (§6.3). BSSE handled by
   ≥TZ basis + the gCP built into the 3c methods.

For barrier heights specifically, M06-2X, MN15, and ωB97X-D keep appearing at
the top of main-group benchmarks — and M06-2X is what the field's silicate
papers actually use (§6.2).

### 6.2 What silicate dissolution papers use

The direct 2025 descendant of the dissertation approach — Watts & Kubicki
(*Sci. Rep.* 2025) — runs **M06-2X/6-311G(d,p) with SMD single points on an
Si₈ cluster with 5 explicit waters**, constrained-scan TS location, and gets
Ea = 56 kJ/mol for SN2-type attack at a Q1 silicon vs the experimental quartz
70–75 kJ/mol, explicitly defending small-cluster proxies. The reference level
chosen to train the 2024 Nature Comms reactive silica–water ML potential was
**ωB97X-D3/def2-TZVP** on 400k cluster configurations — a good statement of
what that community considers the target quality. Current cluster-practice
norms: 8–20 tetrahedral/octahedral centers (~50–150 atoms), OH/H₂O
termination, peripheral-atom constraints to mimic lattice resistance, 3–6
explicit waters at the reactive site, implicit solvent beyond that.

### 6.3 Solvation

Standard practice is **cluster-continuum**: implicit solvent (SMD in
Gaussian-world, CPCM in ORCA) *plus* explicit first-shell waters that shuttle
protons — pure continuum can be qualitatively wrong for water-mediated
reactions. Both PySCF/GPU4PySCF and ORCA provide SMD/CPCM with analytic
gradients (and Hessians in GPU4PySCF), so this is fully covered. COSMO-RS
(openCOSMO-RS in ORCA) is the top tier for ΔG_solv of dissolved products
(Si(OH)₄, Al(OH)₄⁻); charged transition states remain the known
several-kcal/mol error source across the field.

### 6.4 The calibration layer

Modern barrier papers calibrate DFT against **DLPNO-CCSD(T)/CBS (TightPNO)**
single points — MAD ≤ 0.4 kcal/mol vs canonical CCSD(T). On this hardware
that's hours-to-a-day for 30–80-atom clusters in ORCA. The 1990s dissertation
could never afford this layer; we get it for free, and it converts "which
functional do we trust for Si–O–Al?" from a literature argument into a
half-day calculation.

---

## 7. The domain literature since the dissertation

### 7.1 The classics still cited

Xiao & Lasaga (GCA 1994, 1996): HF/MP2 disilicate and Si–O–Al dimer clusters;
H₃O⁺-catalyzed ~24 kcal/mol, OH⁻-catalyzed ~19 kcal/mol, pentacoordinate Si
intermediate; Si–O–Al easier than Si–O–Si. Pelmenschikov et al. (2000, 2001):
"lattice resistance" (surface-embedded barriers 5–16 kcal/mol above free
dimers) and "self-healing" (first-bond rupture ~205 kJ/mol because the bond
re-forms) — the reason cluster termination/constraint strategy matters.
Criscenti, Kubicki & Brantley (2006): the modernized cluster template (Q3
site, H₃O⁺ + 4 explicit waters, B3LYP). Nangia & Garrison (2008): pH-resolved
quartz pathways + TST, and a stochastic aluminosilicate dissolution model —
the closest published analog of the dissertation's architecture.

### 7.2 What replaced them (2010–2026)

- **Periodic DFT + dimer-method TS search** (Izadifar et al. 2022–2025,
  VASP): metakaolinite Si–O–Si 169 kJ/mol vs **Si–O–Al 86 kJ/mol**; alkaline
  aluminate detachment 27–165 kJ/mol depending on cation; Al release cheaper
  than Si — now accelerated by on-the-fly ML force fields.
- **AIMD + metadynamics in explicit water** — the current center of gravity:
  Schliemann & Churakov (GCA 2021) on pyrophyllite edges (multi-step
  detachment with intermediates); Guo et al. (2023/24) on gibbsite step edges
  (two Al-detachment pathways, ~50% barrier reduction at high pH — directly
  relevant since kaolinite's octahedral sheet is gibbsite-like); Liu & Ruiz
  Pestana (2024): the connectivity ladder **Q1 54 / Q2 71 / Q3 81 kJ/mol** for
  Si–O–T hydrolysis, SN2-type pentacoordinate TS.
- **DFT→KMC pipelines** — the dissertation's architecture, validated: Martin
  et al. (ACS Earth Space Chem. 2021) build a quartz KMC from
  connectivity-resolved DFT cluster barriers and jointly reproduce rate *and*
  activation energy (open-source code: **KIMERA**); Kurganskaya (2024) adds
  explicit protonation-state tracking (site pKas → pH-dependent rates) for
  phyllosilicate KMC.
- **ReaxFF** exists for these systems but the MLIP camp documents it can't
  reproduce aqueous condensation reactions quantitatively; **reactive MLIPs**
  (Nature Comms 2024 silica–water, trained on ωB97X-D3 clusters) are the
  emerging bridge.

**Cluster models remain respectable** — actively published and defended by
Kubicki's group in 2025 — with known systematics: free clusters underestimate
lattice resistance (+20–70 kJ/mol possible), connectivity shifts barriers by
30–50 kJ/mol, protonation state by up to ~50%.

### 7.3 Kaolinite targets for the modern study

- Structure: 1:1 dioctahedral clay; reactivity dominated by **pH-dependent
  edge sites** (silanol, aluminol, bridging Si–O–Al / Al–OH–Al); basal faces
  near-inert (KMC + experiment agree).
- Site pKas from FPMD (Liu, Sprik et al. 2013; Kurganskaya 2024 compilation):
  Si–OH ≈ 6.9, Al–OH₂OH ≈ 5.7, bridging Si–O ≈ 2, octahedral Al–OH₂ ≈ 10 —
  these set the cluster charge/protonation states worth computing.
- Experimental calibration targets: acidic kaolinite dissolution Ea ≈ 29
  kJ/mol (Ganor et al. 1995; Cama et al. 2002), alkaline 33–51 kJ/mol (Bauer
  & Berger 1998) — well below quartz's 66–90 kJ/mol, consistent with easier
  Si–O–Al hydrolysis and Al-first leaching.
- The mechanistic template maps 1:1 onto the KMC's site taxonomy (100s Al,
  200s Si, 300s Si–O–Si, 400s Si–O–Al₂, 500s Al–OH–Al): the modern upgrade is
  computing that reaction table's barriers **per connectivity and protonation
  state** at ωB97M-V//r²SCAN-3c + SMD with DLPNO-CCSD(T) calibration, rather
  than HF/MP2 gas-phase dimers.

---

## 8. Verdict: assemble, don't build

Rolling our own DFT engine would be rediscovering libcint, libxc, and twenty
years of SCF convergence tricks. The physics being public doesn't make the
engineering small: analytic Hessians with solvation — the feature this
project lives on — exist in exactly three free packages (PySCF/GPU4PySCF,
ORCA, NWChem), each representing decades of work.

What *doesn't* exist anywhere, and is genuinely ours to build, is the
platform layer:

1. **Cluster builder**: generate terminated, constrained aluminosilicate edge
   clusters from kaolinite crystallography, per site type in the KMC state
   encoding (100s–500s), per protonation state.
2. **Barrier pipeline**: xtb/CREST screening → r²SCAN-3c or B3LYP-D4
   opt+freq (GPU4PySCF) → geomeTRIC/Sella TS + IRC verification → ωB97M-V +
   SMD single points → quasi-RRHO ΔG‡ → Eyring k(T) with tunneling
   correction; ORCA DLPNO-CCSD(T) spot checks.
3. **The bridge**: emit the site-resolved rate table in exactly the format
   `kmc-rs` consumes (the modern `data.rxn`), with provenance for every
   number.

That's a modest Python package — the QM analog of what `kmc-rs` is to the
C++ model — standing on PySCF + ASE + geomeTRIC/Sella, with every heavy
component free and (except ORCA) open source.

## 9. Pilot plan

**Phase 0 — install & smoke test** (an evening): `pip install pyscf
gpu4pyscf-cuda12x geometric sella ase xtb goodvibes` (+ ORCA download after
registration; + `nvidia-cuest-cu12` to test the emulation path). Verify the
4090: B3LYP/def2-TZVP energy+gradient+Hessian on Si(OH)₄·nH₂O, CPU vs GPU
timings.

**Phase 1 — reproduce the dissertation's benchmark** (a weekend):
(HO)₃Si–O–Si(OH)₃ and (HO)₃Si–O–Al(OH)₃⁻ + H₂O/H₃O⁺ hydrolysis — the Xiao &
Lasaga clusters — at r²SCAN-3c geometries, ωB97M-V/def2-TZVPD + SMD energies,
DLPNO-CCSD(T)/CBS calibration. Compare against the 1990s numbers and the
modern literature values (§7). This single exercise validates every layer of
the stack.

**Phase 2 — the barrier ladder**: connectivity- and protonation-resolved
barriers for the KMC's five site families; benchmark one or two against
CP2K metadynamics or a fine-tuned MLIP as an independent check.

**Phase 3 — close the loop**: regenerate `data.rxn` from computed k(T) and
run `kmc-rs` against the dissertation's original parameterization.

---

## Appendix A — license flags

| Item | Flag |
|---|---|
| ORCA | Free academic/personal only; closed; no redistribution; commercial use needs paid license |
| GAMESS-US | Free but no redistribution (blocks containers); registration |
| geomeTRIC | BSD-3 + added "no AI training datasets" clause |
| pysisyphus, SevenNet | GPL-3.0; pysisyphus unmaintained |
| MACE-OMOL / OMAT / MH / OFF | ASL — non-commercial research only (MACE-MP/MPA are MIT) |
| UMA weights | Gated FAIR Chemistry License — commercial OK, geo-restrictions, access request |
| ReactBench | CC BY-NC-SA (non-commercial) |
| AIMNet2 | MIT, but **no aluminum** — unusable here |

## Appendix B — key sources

Engines & GPU: PySCF (github.com/pyscf/pyscf; pyscf.org/user/gpu.html);
GPU4PySCF papers arXiv:2404.09452, arXiv:2407.09700; RTX 4090 vs A100 data
arXiv:2605.06489; INT8-DF on 4090 arXiv:2601.08077; NVIDIA cuEST
(docs.nvidia.com/cuda/cuest); Rowan benchmarks (rowansci.com/blog/gpu4pyscf);
ORCA (faccts.de/orca); Psi4 (psicode.org); NWChem (nwchemgit.github.io);
mixed-precision accuracy: Luehr et al. JCTC 2011 (10.1021/ct100701w),
Kussmann & Ochsenfeld JCP 2021.

Tooling: geomeTRIC (github.com/leeping/geomeTRIC); Sella
(github.com/zadorlab/sella); autodE (github.com/duartegroup/autodE); Arkane
(RMG-Py docs); Shermo (sobereva.com/soft/shermo); GoodVibes
(github.com/patonlab/GoodVibes); Pilgrim (github.com/cathedralpkg/Pilgrim);
quacc (github.com/Quantum-Accelerators/quacc); ReactBench findings
(PMC12442704).

Methods & domain: Bursch et al., Angew. Chem. 61, e202205735 (2022); Grimme
r²SCAN-3c JCP 154, 064103 (2021); Watts & Kubicki, Sci. Rep. 15 (2025); Xiao
& Lasaga GCA 58, 5379 (1994) & GCA 60, 2283 (1996); Pelmenschikov JPC B 104,
5779 (2000) & JPC A 105 (2001); Criscenti et al. JPC A 110, 198 (2006);
Nangia & Garrison JPC A 112, 2027 (2008); Izadifar et al. Nanomaterials 13,
1196 (2023) & Nanoscale Adv. 7 (2025); Schliemann & Churakov GCA 293 (2021);
Guo et al. (gibbsite, 2023/24); Liu & Ruiz Pestana npj Mater. Degrad. 8
(2024); Martin et al. ACS Earth Space Chem. 5, 516 (2021) + KIMERA
(github.com/hegoimanzano/KIMERA); Kurganskaya Minerals 14, 900 (2024);
silica–water MLIP Nat. Commun. 15 (2024); kaolinite edge pKas: Liu et al. GCA
117, 180 (2013); Pouvreau et al. JPC C 123, 11628 (2019); experimental Ea:
Ganor et al. GCA 59 (1995), Cama et al. GCA 66 (2002), Bauer & Berger Appl.
Geochem. 13 (1998).
