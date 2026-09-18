# MSR AI-for-Science scan — what fits petra / quarry (2026-09-17)

*Scoping note. Victor asked (2026-09-17) for a pass over the Microsoft
Research AI for Science group's recent publications to see what is
useful to the petra/quarry platform, or could seed another platform
variant. Verdict first, evidence after. Written under fleet POLICY §12
(platform-test compute envelope): nothing here proposes a production
tier; the one executable item is a bounded survey-tier functional test.*

## TL;DR

| Item | What it is | Fit | Action |
|---|---|---|---|
| **Skala** (neural exchange–correlation functional) | ML XC functional at meta-GGA cost; GMTKN55 WTMAD-2 2.8 kcal/mol incl. barrier heights; MIT; **native PySCF / GPU4PySCF 1.8.1** | Direct hit on quarry. Same engine, same CUDA-12 stack; `skala-cuda12x==2026.9` resolves against quarry's exact pins with torch as the only new dependency. | **Card A2e** — single points on A2a's four hash-bound si-neutral structures vs the banked CCSD(T) focal barrier. Minutes of GPU. |
| Anisotropic message-passing MLIP for enzymatic *reactions* (Riniker group, JACS 2026) | An MLIP architecture built to get reaction barriers right, not just equilibrium PES | Method-relevant to quarry's TS-hunt (the CI-NEB pre-climb that burned A2a/A3 GPU hours), not to petra directly | Note only. Becomes a card if the A-track resumes and MLIP pre-climb is worth a bounded unit. |
| MatterSim (universal MLIP, MSR) | Bulk/mineral phases 0–5000 K, to 10 Mbar; MIT; LAMMPS multi-GPU | Already in SURVEY §5 as "not reaction-specialized". Universal MLIPs are weakest on surfaces and barriers — exactly petra's domain. | None. Reaffirms SURVEY §5's division of labor (MLIPs propose, DFT owns numbers). |
| MatterGen (generative crystal structures under property constraints) | Diffusion-style inverse design of materials | Not our chemistry, but a real **platform-variant seed**: generative structure proposer → automatic petra deck → kinetic screen | Design question for after M5, not a card. §4 below. |
| Skala in CP2K via GauXC (Aug 2026), FHI-aims roadmap (Apr 2026) | Skala reaching periodic codes | Relevant only if Track E's periodic DFT (E3b, CP2K) ever wants a Skala tier | None now; E3b closed `disagrees` on PBE-D3 frozen paths, and the question there was the path, not the functional. |
| Everything else on the feed | Aurora weather model, protein fitness RL, synthesis-procedure reasoning, SimPoly polymer force fields, "mindless molecules" chemical space, diffeomorphic optimization, optical-property scaling laws | Not relevant | — |

The one thing worth spending compute on is Skala, and only because the
comparison it needs already exists: A2a banked a canonical-TZ + TightPNO
CBS focal-point barrier for si-neutral with receipts. A Skala single
point on those geometries is a survey-tier data point that either earns
Skala a place in the §12 survey-tier list or rejects it for this
chemistry — in either case a banked, honest result for under an hour of
workstation time.

## 1. Posture: POLICY §12 and what "interesting" is allowed to mean

POLICY v16 §12 says the program is *testing a platform* on one consumer
GPU: survey-tier methods are the banked values, no CC calibration, no
re-tiering, 4 h wall per QM unit, nothing labelled production. Read
against that, the MSR feed sorts itself:

- A cheaper functional with barrier-height accuracy is *exactly* a
  survey-tier question, and it can be answered against an anchor that
  already exists (A2a), so no new CC is implied.
- Anything that would only matter for a production-grade barrier (the
  MLIP pre-climb, Skala-in-periodic-codes) is noted, not filed.
- A platform-variant idea (MatterGen → petra) is a design conversation
  for after M5 (first science on the platform), which is where the
  Omnibus already puts new-track questions.

## 2. Skala — the fit for quarry

### 2.1 What it is

Skala is Microsoft Research's deep-learning exchange–correlation
functional (Luise, Huang, Vogels et al., arXiv 2506.14665, first posted
Jun 2025, revised Apr 2026; PyPI/GitHub `microsoft/skala`; weights
`microsoft/skala-1.0` and `-1.1` on Hugging Face, MIT). A ~385k-parameter
network takes meta-GGA-level grid features of the density and learns the
non-local representation directly, so it runs at semi-local cost (the
model card says ~3× r²SCAN) rather than hybrid or double-hybrid cost.
Reported accuracy: ~1 kcal/mol MAE on W4-17 atomization energies and
**2.8 kcal/mol WTMAD-2 on GMTKN55**, the 55-subset benchmark whose
categories include reaction barrier heights and thermochemistry. Skala
1.1 is the recommended version (2.5× the training data of 1.0, per-atom
packed grids). Training data: ~78k reactions from MSR-ACC plus
conformers, IP/EA/PA, noncovalent, distorted geometries and elementary
reactions at CCSD(T)/CBS quality, elements H–Xe.

### 2.2 Why it fits our stack (verified tonight, not assumed)

The primary reference implementation is **PySCF and GPU4PySCF** — the
engines quarry already runs (SURVEY §2.2, §3.2). Checked on the
workstation on 2026-09-17:

- quarry's `uv.lock` pins `pyscf==2.14.0`, `gpu4pyscf-cuda12x==1.8.1`,
  `cupy-cuda12x==14.1.1`, Python 3.13.1. Skala's documented support
  matrix is PySCF 2.14, GPU4PySCF 1.8.1, CUDA 12/13, Python 3.11–3.13.
- A `uv pip install --dry-run` of `skala-cuda12x` against those exact
  pins resolves with **no conflicts**: it adds `skala==2026.9`,
  `skala-cuda12x==2026.9`, `torch==2.14.0`, and pulls
  `pyscf-dispersion==1.5.0` (already in the venv). torch is the only
  material new dependency (~3 GB; `/mnt/data` has 1.4 TB free). The qm
  venv currently has no torch.
- API is a drop-in KS object: `from skala.gpu4pyscf import SkalaKS;
  ks = SkalaKS(mol, xc="skala-1.1"); ks.kernel()` (CPU:
  `from skala.pyscf import SkalaKS`). Gradients exist (an ASE calculator
  with geometry optimization is documented); Hessians are not documented.

So the integration cost is an optional `skala` extra in `qm/pyproject.toml`
and a single-point driver; no ORCA, no CP2K, no new engine.

### 2.3 Known gaps to respect (from the model card and README)

These bound what the spike may claim:

- **Molecular, gas-phase training.** GMTKN55 is molecular. Nothing is
  known about Si–O–Si hydrolysis transfer — that is the question the
  spike asks, so the answer must be reported as one data point, not a
  validation.
- **No solvation model documented.** PySCF's SMD/PCM are separate
  objects that wrap a KS instance, so composition may simply work, but
  it is untested by the authors. The comparison must be gas-phase
  against a gas-phase anchor; the CC focal barrier *is* gas-phase, so
  that is natural. Any SMD-on-Skala number is a curiosity, not a result.
- **No dispersion treatment documented.** Skala's training includes
  noncovalent sets, so it may already carry medium-range correlation; a
  post-hoc D4 could double-count. The spike runs bare and +D4 as two
  labelled rows and does not choose between them.
- **Closed-shell RKS only documented**; no PBC; largest tested system
  ~180 atoms at def2-TZVP. The si-neutral dimer is 18 atoms (H8O8Si2),
  so none of this binds.
- The model card says "this is not a production model" and "test the
  functional further before applying it" — which under §12 is the same
  standing every method on this box has.

### 2.4 The anchor that makes the spike cheap

A2a (done 2026-09-05) left, hash-bound, in
`/mnt/data/vsletten/dissertation-data/task208-a2a-path-rebuild-20260825/production-closeout/`:

| Method (all on the same four r²SCAN-3c geometries) | Solvent | ΔE‡ reactant → addition-TS, kJ/mol |
|---|---|---|
| canonical CCSD(T)/cc-pVTZ + TightPNO TZ/QZ CBS correction (**focal point**) | gas | **132.960133** |
| canonical CCSD(T)/cc-pVTZ | gas | 128.085164 |
| TightPNO DLPNO-CCSD(T)/cc-pVTZ | gas | 128.330864 (DLPNO − canonical = +0.2457, gate pass) |
| r²SCAN-3c / def2-mTZVPP / D4 / gCP | gas | 125.742278 |
| B3LYP / def2-SVP / DF | gas | 114.173360 |
| ωB97M-V / def2-TZVPD | SMD(water) | 134.504248 *(not gas-comparable)* |
| B3LYP-D4 / def2-TZVPD | SMD(water) | 148.293616 *(not gas-comparable)* |

Structures in `store.sqlite` (SHA-256
`480cc244cf06bf7fd8edce86e4d57c3a8466e24e2793daf39be3fe27207dcd89`):
ids 1 reactant, 2 intermediate, 3 addition-transition-state, 4
released-product; geometry hashes listed in card A2e. Receipts
`dft-summary.json` (`b6f074c7…`) and
`cc-calibration/focal-point-summary.json` (`2df64b9f…`) rehash to the A2a
card's values (checked 2026-09-17).

Two things fall out. First, the gas-phase functional spread on this one
barrier is already 19 kJ/mol (B3LYP/SVP to the focal point), so a
survey-tier functional that lands within ~4 kJ/mol would be a real
improvement over what the A-track banks today. Second, the existing
ωB97M-V number is SMD and therefore *not* comparable to the CC anchor;
the spike adds a gas-phase ωB97M-V/def2-TZVPD single point (two SCFs) so
the ranking is solvent-consistent — a cheap repair of an inconsistency
the ledger already has.

### 2.5 The spike (card A2e), in one paragraph

Sixteen single points — the four A2a structures × {Skala-1.1/def2-TZVP,
Skala-1.1/def2-TZVPD, Skala-1.1/def2-TZVP + D4, ωB97M-V/def2-TZVPD
gas} — into a fresh quarry store, receipts hashed, the A2a store
untouched. Pre-declared verdict bands on |ΔE‡(Skala) − 132.960| for the
gas-phase rows: ≤ 4.2 kJ/mol (1 kcal/mol) *adopt as a survey-tier
single-point functional candidate*; ≤ 8.4 *usable with stated
uncertainty*; > 8.4 *reject for this chemistry*. One reaction, one data
point; no geometry changes, no frequencies, no CC, no SMD claims, no
production language. Expected wall time is minutes; the §12 cap of 4 h
is stated for form.

### 2.6 What a pass would change

- POLICY §12's survey-tier list would gain "Skala-1.1 single points on
  r²SCAN-3c geometries" — that is a POLICY edit and Victor's call, not
  the card's.
- A2's functional ranking gains a column; A2 stays blocked on Victor per
  the A9b ruling — A2e does not close or unblock it.
- Future A-track units (if the observable-release regime question is
  ever ruled) get a barrier tier that is cheaper than a hybrid and
  closer to CC than r²SCAN-3c on this reaction, still inside §12.

A fail is equally useful: it bounds how far a molecular-trained neural
functional transfers to silicate hydrolysis, which is worth one line in
SURVEY §6.

## 3. MLIPs — MatterSim and the reaction-capable MLIP

SURVEY §5 already lists MatterSim as "bulk/mineral phases across T/P;
not reaction-specialized" and carries the ReactBench caveat: pretrained
universal MLIPs frequently fail to localize transition states without
fine-tuning, and direct-force variants give broken Hessians. Nothing in
MSR's 2026 MatterSim updates (3–5× inference, LAMMPS multi-GPU,
MatterSim-MT for stress/magnetic/dielectric properties) changes that for
Al/Si/O/H bond breaking. The companion literature is explicit that
surfaces are where universal MLIPs are weakest — petra's whole domain.

The genuinely new method item is the **anisotropic message-passing
multiscale MLIP for enzymatic reactions** (Thürlemann, Pultar, Gordiy,
Ruijsenaars, Riniker; JACS, July 2026). It is a reaction-first MLIP
design. Where it would plug into this platform is quarry's saddle
search, not petra: A2a spent its GPU budget on CI-NEB pre-climb and
directed Sella refinement at r²SCAN-3c; an MLIP that reproduces barrier
topology could supply the path and the TS guess, with DFT owning every
number that reaches a deck (SURVEY §5's rule). That is a bounded card
*if* the A-track resumes; today the track is blocked on the
observable-release regime, so it is a note.

## 4. MatterGen → petra: the platform-variant seed

MatterGen proposes crystal structures under property constraints. On its
own it is materials discovery. The reason it belongs in this note is
petra's shape: a chemistry-free lattice KMC engine whose mineral is a
declarative TOML deck, plus quarry, which already *emits* deck fragments
(`by_count`/`when` tables) from calculations. That is most of a pipeline
in which a generator proposes defected or doped structures, a compiler
turns each into a deck, and petra screens their alteration kinetics —
inverse design of *reactivity* rather than forward simulation of one
known lattice. It also lines up with the defect/strain machinery petra
already has.

What is missing, honestly: a general structure → deck compiler (the
kaolinite deck was hand-derived from the legacy cell; `quarry/crystal.py`
cuts clusters *from* a deck, not the reverse), a site-taxonomy inference
step, and the rates themselves (which is the whole A-track problem
again, multiplied by every generated structure). None of that is a
weekend, and none of it is M5. The right time to raise it is when the
Omnibus asks what Track B does after B5 — file it there, not on the
board now.

## 5. Sources

- MSR AI for Science publications feed:
  https://www.microsoft.com/en-us/research/lab/microsoft-research-ai-for-science/publications/
- Skala paper: https://arxiv.org/abs/2506.14665 ; project page
  https://www.microsoft.com/en-us/research/publication/accurate-and-scalable-exchange-correlation-with-deep-learning/
- Skala code + README (install, `SkalaKS`, support matrix):
  https://github.com/microsoft/skala ; weights
  https://huggingface.co/microsoft/skala-1.1
- Skala in CP2K via GauXC (Poschel, Pototschnig, Stein et al., Aug 2026) — on the MSR feed.
- Thürlemann et al., "Multiscale Neural Network Potential with
  Anisotropic Message Passing for the Fast and Accurate Simulation of
  Protein Dynamics and Enzymatic Reactions", JACS, July 2026 — on the MSR feed.
- MatterSim: https://www.microsoft.com/en-us/research/blog/mattersim-a-deep-learning-model-for-materials-under-real-world-conditions/ ;
  MatterSim-MT (May 2026) on the MSR feed.
- Universal-MLIP surface caveat: https://arxiv.org/abs/2403.04217
- Local anchors: A2a card and evidence root (§2.4); `qm/SURVEY.md` §5, §6.4;
  mission-control `POLICY.md` §12 (v16).
