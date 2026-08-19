# The Omnibus — one engine, many worlds

*Drafted 2026-08-18. The program plan unifying the in-flight KMC/QM work,
the generalization of petra into a multi-strategy lattice simulation
platform, and the first two out-of-domain showcases: pitting corrosion
and astrochemistry. Owner: Victor + Fable. This document is the map;
each track keeps its own detailed scoping doc under `docs/scoping/`.*

**The thesis, in one paragraph.** Lasaga's career move was importing
physical chemistry's standard tools (QM, TST, Monte Carlo) into
geochemistry. The 2026 version of that move is: *importing
physics-grade stochastic simulation rigor — declarative, reproducible,
exactly-sampled, massively-ensembled — into every field still running
bespoke single-shot lattice scripts.* Petra becomes the platform;
quarry supplies first-principles rates where chemistry is involved; and
each new field costs one deck plus one validation campaign, not a new
codebase. The dissertation's architecture (QM barriers → stochastic
lattice dynamics), built for kaolinite in 1999, turns out to be a
general-purpose research instrument.

---

## Program structure

Four tracks. A and B are load-bearing; C and D are the showcases that
prove generality. Tracks run concurrently but B1 (the schema RFC) gates
the showcases.

```
A. Geochemistry science        (in-flight; uses petra+quarry as they are today)
B. Platform generalization     (petra → multi-strategy lattice engine)
C. Showcase: pitting corrosion (new field, same math)
D. Showcase: astrochemistry    (new sky, same architecture — QM rates + KMC)
```

---

## Track A — Geochemistry (in-flight and flagship)

The scientific home base. Everything here runs on current petra/quarry;
nothing waits on Track B.

**A1. Phase 1 completion — the Xiao & Lasaga ladder.** ✅ si-neutral
DONE (ΔG‡ = 113.0 kJ/mol = 27.0 kcal/mol vs X&L ~29; merged, archived
with provenance). Remaining:
- a. **TASK-164** (queued in mission-control): scan-smoothness repair
  (outlier re-seeding from both neighbors) + al-neutral barrier. The
  KMC-relevant deliverable: confirm Si–O–Al < Si–O–Si ordering.
- b. **Acid mechanisms** (si-acid, al-acid): needs the
  protonated-bridge pathway (H₃O⁺ protonates O_br first, then water
  attacks) — a different scan design, not just a different deck.
- c. **Base pathway** (OH⁻ attack) to complete the X&L trilogy (~19
  kcal/mol anchor).

**A2. Production-tier energetics.** Upgrade from B3LYP/def2-SVP first
pass to the survey's protocol: r²SCAN-3c geometries, ωB97M-V/def2-TZVPD
+ SMD single points, D4 dispersion; ORCA DLPNO-CCSD(T) spot calibration.
*Victor action: ORCA registration + install (never committed).* Also:
full Gonzalez-Schlegel IRC to replace quick-IRC on soft-mode saddles.

**A3. Phase 2 — the barrier ladder.** The crystallographic cluster
builder (terminated, constrained edge clusters cut from kaolinite
coordinates per KMC site family; the frozen-atom contract is already
plumbed end-to-end), cluster-continuum microsolvation (3–6 explicit
waters), barriers by connectivity × protonation state.

**A4. Phase 3 — close the loop.** Emit the kaolinite deck's rate table
from computed k(T) (+ regenerate legacy `data.rxn`), run petra and
kmc-rs, compare against experimental Ea (kaolinite ~29 kJ/mol acidic,
33–51 alkaline) and the dissertation's original parameterization.

**A5. FLAGSHIP — the field–lab discrepancy.** Scoped in
[`docs/scoping/field-lab-discrepancy.md`](scoping/field-lab-discrepancy.md).
The short version: nobody has quantitatively attributed how many of the
2–4 missing orders of magnitude in field weathering rates come from
intrinsic lattice-state evolution (defect depletion, surface aging,
rate-spectrum narrowing). Petra's defect+strain machinery plus quarry's
provenance-carrying barriers make the first defensible attribution
possible. Phases: (0) three new observables — rate-spectra exporter,
BET-vs-geometric area accounting, exposure-age tracking (these are
Track-B observables-framework features; build once, share); (1) generic
aging study incl. the Blum–Yund–Lasaga two-regime validation gate;
(2) QM-grounded kaolinite deck (consumes A3/A4); (3) the "discrepancy
budget" capstone paper; (4) optional ERW feedstock-aging note.
Known hard problem: simulated-years timescales need superbasin/
acceleration work (see B7).

---

## Track B — Platform: petra → general lattice simulation engine

**The design insight.** A petra deck is already a declarative spatial
stochastic state machine: state vocabulary + neighborhood predicates +
guarded transition rules + rates. Continuous-time KMC is just one
**update strategy**. Generalizing the strategy while keeping the deck
formalism turns petra into an engine for every lattice-model dialect in
science — and the moat is not the engine but the *practice* it enables:
models as portable documents, deterministic replay, and ensembles with
real statistics in fields that publish n=1 animations.

**B1. Architecture RFC (gates C and D).** One design doc, reviewed
before code. Contents:

- **Deck schema v2**, split into orthogonal sections:
  - `[structure]` — lattice/topology: unit cell + dims (today's case),
    plus simple-grid shorthands (square/hex/cubic with Moore or von
    Neumann neighborhoods) and, later, arbitrary graphs. Site kinds,
    state vocabularies unchanged in spirit.
  - `[dynamics]` — the transition system: today's reactions
    (center/guards/effects/modifiers) generalized to *rules*, with the
    rate field interpreted per strategy (propensity for CTMC,
    probability for PCA, energy delta for Metropolis, truth-table for
    deterministic CA).
  - `[execution]` — strategy selection + parameters:
    `strategy = "ctmc" | "synchronous" | "metropolis" | "pca"`,
    temperature/thermo where meaningful, stopping conditions, ensemble
    spec (n_replicas, seed policy).
  - `[observables]` — declarative outputs: state counts, event rates,
    rate spectra, spatial statistics (cluster-size distributions,
    correlation lengths, interface roughness), time series cadence.
  - v1 decks must load unchanged (compat shim: today's layout implies
    `strategy = "ctmc"`).
- **UpdateStrategy trait** (Rust): the engine core owns lattice state,
  neighborhood queries, rule matching, RNG streams, and event/trajectory
  recording; strategies own *selection and time*: ExactCTMC (today's
  Gillespie), SynchronousCA (double-buffered), AsyncMetropolis
  (Glauber/Metropolis acceptance), DiscreteTimePCA (per-site independent
  probabilities). Tau-leaping later.
- **Determinism contract**: same deck + seed + strategy ⇒ bitwise
  trajectory, per strategy, forever (the golden-gate discipline that
  made kmc-rs trustworthy, promoted to platform law).

**B2. Engine refactor behind the trait.** Mechanical but delicate;
gated by bitwise parity of the kaolinite and kossel decks against
current petra over long trajectories (the same milestone discipline as
the C++ → Rust port).

**B3. Conformance decks — the "it's general now" proof.** Three tiny
decks with analytic answers, kept as CI gates:
- Conway's Life (synchronous): glider period/displacement, known
  oscillators.
- 2-D Ising (Metropolis): critical temperature via Binder cumulant
  crossing vs Onsager's exact result.
- SIR epidemic (PCA): threshold behavior vs mean-field prediction.
The Ising gate doubles as the ensemble-statistics smoke test.

**B4. Ensembles + observables as first-class.** Multi-seed parallel
runner (rayon across replicas), aggregation (means, distributions,
bootstrap CIs), and the declarative `[observables]` implementations.
*Design with A5-Phase-0 in the room:* the rate-spectra exporter and
area observables the field–lab project needs are the first real
customers of this framework.

**B5. Visualization generalization.** The WASM viewer learns to render
arbitrary decks (state-color mapping declared in the deck), ensemble
summary views (small-multiple replicas, distribution panels), and keeps
trajectory record/replay working across strategies.

**B6. Tool paper + possible repo extraction.** JOSS or Computers &
Geosciences. Decision point (deliberately deferred until B3 is green):
extract petra to its own repository for standalone citation, or keep
the monorepo. Extraction also forces the licensing/naming decisions.

**B7. Acceleration (serves A5 and C).** Superbasin escape /
first-passage-time methods for stiff systems (years of simulated time).
This is real algorithmic work — scheduled after B3, informed by A5's
concrete stiffness profile.

---

## Track C — Showcase 1: pitting corrosion

*Scoping in flight → `docs/scoping/pitting-corrosion.md`.*

The bet: pitting is mathematically isomorphic to mineral dissolution —
localized attack at defects, film-breakdown-vs-repair competition,
stochastic initiation, pit morphology evolution — and the corrosion
modeling literature (ad-hoc CA rules, mean-field point-defect models,
phase-field) lacks exactly what petra brings: exact stochastic
dynamics, environment-dependent barriers with provenance, and matched
ensemble statistics against experimentally-Weibull pit-initiation
distributions.

- **C1. Scoping** (subagent, in flight): modeling literature map, the
  gap, minimal-viable-model sketch, validation datasets, reviewer-
  skepticism assessment.
- **C2. Minimal viable deck**: passive-film lattice (states: intact
  film / thinned / bare / dissolving / repassivated), chloride attack
  as guarded transitions, electrochemical potential as the Δμ-like
  knob, defect sites seeding initiation. Runs on B2 engine, CTMC
  strategy — corrosion is genuinely kinetic, so the showcase also
  proves the *original* strategy in a new domain.
- **C3. The statistics paper**: ensemble pit-initiation-time
  distributions vs published stainless-steel data; metastable pitting
  transients as event-rate time series; morphology comparisons.
- **C4. Stretch**: quarry on oxide clusters (Cl⁻ attack on
  Fe/Cr-oxide) for first-principles barrier inputs — the full
  dissertation architecture, transplanted.
- Venue targets: Corrosion Science, npj Materials Degradation.

## Track D — Showcase 2: astrochemistry

*Scoping in flight → `docs/scoping/astrochemistry.md`.*

The bet: the astrochemical rate databases (KIDA, UMIST) are full of
order-of-magnitude guesses and openly solicit computed rates; grain-
surface chemistry at 10 K is *tunneling-dominated* (quarry's Eckart
machinery is the whole story, not a correction) and is modeled by
lattice KMC on ice mantles — the same QM-rates-feed-lattice-KMC
architecture, in a field where both halves are undersupplied and every
rate is a citable community contribution. JWST's ice observations are
creating demand for exactly this kinetics *right now*.

- **D1. Scoping** (subagent, in flight): database contribution
  pathways, the honest niche map (where an outsider's GPU rates are
  welcomed vs crowded), Cuppen-group landscape, validation anchors.
- **D2. Rates campaign**: quarry pointed at flagged-uncertain
  reactions — H-abstractions on ice analogs, the CO→CH₃OH
  hydrogenation ladder — with tunneling; systems are tiny (the GPU
  yawns), throughput is bounded by TS-hunting judgment, not compute.
- **D3. Ice-mantle deck**: adsorption / thermal+tunneling diffusion /
  Langmuir-Hinshelwood reaction / desorption on an ice lattice;
  stochastic small-grain effects (the reason rate equations fail and
  Monte Carlo is required) as the headline observable.
- **D4. Validation**: Watanabe & Kouchi CO-hydrogenation experiments,
  TPD binding energies, comparison against published grain-chemistry
  models; JWST ice abundances as the observational anchor.
- Venue targets: A&A / ApJ (rates + model), KIDA database submissions.

---

## Sequencing and milestones

```
M0  (now)        A1a TASK-164 + scoping docs land (C1, D1)
M1               B1 RFC written and reviewed          ← gates everything B
M2               B2 refactor, bitwise parity green
M3               B3 conformance decks green (Conway, Ising, SIR)  ← "general" proven
M4               B4 ensembles/observables (built against A5-Phase-0 needs)
M5   parallel:   A5 Phase 1 (aging study)  |  C2 corrosion deck  |  D2 rates campaign
M6               First external deliverables: petra tool paper, corrosion statistics
                 paper, first KIDA submissions, field-lab mechanism paper
M7               A2-A4 production geochem + A5 capstone; C3/C4; D3/D4
```

Concurrency rule of thumb: one heavy compute campaign on the box at a
time (compute etiquette below); writing/scoping/refactoring tracks run
freely alongside.

## House rules (standing, learned the hard way)

- **Compute etiquette**: GPU-first; CPU stages capped
  (`OMP_NUM_THREADS≤16`); every long job teed to a log file — *nothing*
  depends on terminal scrollback; no multi-hour full-CPU jobs without
  Victor's explicit go.
- **Repo discipline**: worktree-per-branch, PR-per-unit-of-work,
  auto-merge, Sourcery threads resolved (conversation resolution is a
  required check), learnings captured in `.claude/learnings/`.
- **Numbers discipline**: barriers cross tool boundaries, never
  pre-exponentiated rates; every emitted number carries provenance
  (method, geometry hash, date); golden/bitwise gates on every engine
  refactor; a verified saddle is not necessarily the right saddle —
  structural gates stay on.
- **Big outputs** (checkpoints, trajectories, cubes) stay out of git.

## Risks, honestly

| Risk | Track | Mitigation |
|---|---|---|
| Timescale stiffness (simulated years) | A5, C | B7 acceleration; frame claims as intrinsic-ceiling attribution |
| Schema v2 second-system effect | B | RFC discipline; v1 compat gate; three tiny conformance decks, not grand abstraction |
| Outsider skepticism in corrosion/astro journals | C, D | Lead with validation against *their* canonical datasets; contribute to *their* databases first |
| Cuppen group overlap (grain KMC) | D | Differentiate on openness + declarative decks + QM-rate provenance; collaborate, don't compete |
| Amorphous passive films / ice ≠ crystal lattice | C, D | State the lattice idealization up front; it's the standard first approximation in both fields' own literatures |
| One human's attention across four tracks | all | The queue + agent workers; milestone gates prevent half-finished everything |

---

*The unifying joke being that the dissertation's architecture was right
all along — it just needed twenty-five years, a thousand-fold hardware
speedup, and an agent to carry the water. Let's go.*
