# The Whole Thing, Explained

*What's actually going on in this repository — from a 1999 dissertation
to a working research program, written for the intelligent reader who
knows what kaolinite is but hasn't been in the computational weeds.
(Companion documents: [OMNIBUS.md](OMNIBUS.md) is the program plan;
[scoping/](scoping/) holds the per-project research briefs.)*

---

## The short version

Twenty-seven years ago a PhD dissertation proposed a simple, ambitious
architecture: compute the energy barriers of individual chemical
reactions on a mineral surface with quantum mechanics, feed those
barriers into a Monte Carlo simulation that plays out millions of
reactions on a virtual crystal lattice, and watch macroscopic
geochemistry — dissolution rates, surface textures, activation
energies — *emerge* from first principles.

The architecture was right. The hardware wasn't. The quantum half of
the plan needed weeks of supercomputer time per number, so the
dissertation shipped with placeholder rates, and the idea waited.

In 2026, a consumer graphics card computes those numbers in about a
minute each. This repository is what happens when the original
architecture finally meets adequate hardware — plus a staff of AI
agents to do the labor. In one week it: reproduced the field's
classic 1994 benchmark calculation to within 2 kcal/mol, generalized
the simulation engine into a platform that can model systems far
beyond geochemistry, and identified two flagship scientific questions
— one of global relevance (why lab and field weathering rates disagree
by orders of magnitude), and one intensely personal (what ⁴⁰Ar/³⁹Ar
step-heating of muscovite actually measures — the subject of a 1998
paper that started as a Princeton senior thesis).

Everything below is the long version.

---

## 1. Where it started: kaolinite, 1999

The dissertation modeled the dissolution of kaolinite — the humble
clay at the end of nearly every silicate weathering sequence — as a
**kinetic Monte Carlo (KMC)** simulation. The idea: represent the
mineral surface as a lattice of sites (aluminum sites, silicon sites,
the bridging oxygens between them), give every possible elementary
reaction (protonate this oxygen, hydrolyze that bond, detach this
silicon) a rate, and let the simulation choose events at random in
proportion to their rates, one at a time, exactly as nature does. Run
long enough, and etch pits form, surfaces retreat, and a dissolution
*rate* emerges that no one put in by hand.

The catch is those input rates. Transition-state theory says each one
follows from an **activation energy barrier** — the energy cost of the
molecular contortion at the moment a bond breaks. Those barriers come
from quantum chemistry, and in the late 1990s a single barrier
calculation on a realistic mineral fragment consumed weeks of a DEC
Alpha workstation. (A labmate of the era, Yitian Xiao, monopolized one
of the group's two Alphas for roughly a year producing the founding
numbers of this field — his 1994-96 papers with Antonio Lasaga are
still cited as the benchmarks.) So the dissertation's KMC ran on
estimated rates, the QM-to-KMC bridge was described rather than built,
and the code was preserved like a time capsule: C, then C++, on the
assumption that someday the other half would become possible.

## 2. What changed

Three things, compounding:

**Hardware.** The workstation this program runs on — a 16-core Ryzen
and an RTX 4090 — measured at **74× faster** than its own 16 CPU cores
for the dominant quantum-chemistry operation (the "Hessian," the
calculation frequency analysis needs). A barrier that took Xiao's Alpha
weeks takes this machine about **42 minutes end to end**, including all
the failed attempts along the way. The entire 1994 landmark study could
be rerun, at higher accuracy, in a day.

**Software.** The open-source quantum chemistry stack matured into
something assembled rather than written: PySCF for the quantum
mechanics, GPU acceleration bolted on, saddle-point searches and
reaction-path methods as libraries. Nobody has to write an integrals
engine anymore; the craft moved up a level, to orchestration.

**Agents.** The labor of research computing — writing the glue code,
running the campaigns, chasing the failures, keeping the books — is
now done by AI agents working from a shared task board in this
repository. More on that at the end, because how this program *runs*
turned out to be a story of its own.

## 3. The two platforms

Everything scientific here stands on two in-house pieces of software,
both descended directly from the dissertation.

**Petra** (Greek: rock) is the dissertation's KMC engine, reborn with
one deep change: *the mineral is data, not code.* A "deck" — a single
declarative text file — describes the crystal structure, the states
each lattice site can be in, and the reactions with their barriers.
The same engine runs any deck. The original kaolinite model exists as
one such deck, verified to reproduce the old C++ code's dynamics
*bit-for-bit* — the past is a regression test. Petra has since grown
crystallographic defects, strain energy, multilayer structure, live
in-browser visualization, and (this week) a generalization that makes
the KMC algorithm itself just one plug-in "update strategy" among
several — the significance of which appears in §6.

**Quarry** (the place stone comes from) is the half of the
dissertation that couldn't be built in 1999: the quantum-chemistry
pipeline. Given a reaction on a mineral-fragment cluster, it finds the
transition state (the molecular configuration at the top of the energy
hill), verifies it honestly (a real transition state has exactly one
unstable vibration — and the pipeline includes hard-won safeguards
against the many ways a search finds a *plausible wrong answer*),
computes the barrier and the quantum tunneling correction, and emits
the resulting rate directly into petra's deck format, every number
stamped with its provenance. The bridge the dissertation described is
now a tested piece of software: a deck of computed rates round-trips
into the engine as an automated test.

## 4. The dissertation, finally

With both halves in hand, the original program is simply... running.

The first result landed this week: the neutral hydrolysis barrier of
the Si–O–Si bridge — water attacking the bond that holds silicate
minerals together, via a five-coordinate silicon transition state with
a proton in mid-flight — computed at **27.0 kcal/mol**, against Xiao &
Lasaga's 1994 benchmark of ~29. Getting there took eight attempts in
three hours, each failure diagnostic and now encoded as a safeguard;
the same failures would each have cost a week of Alpha time in 1996.

The roadmap from here is the dissertation's chapter that never got
written: the full **barrier ladder** — every reaction in the kaolinite
model's repertoire, resolved by how many neighbors each site has and
how protonated it is (both matter enormously: 30–50 kJ/mol shifts) —
then closing the loop by running the kaolinite deck with computed
rates against the original 1999 parameterization and against
experimental dissolution data. A campaign for the Si–O–Si family is
running on the GPU as this document is written.

Beyond kaolinite, the same machinery points at any mineral: olivine
and wollastonite are queued interests because of **enhanced rock
weathering** — the carbon-removal strategy of spreading crushed
silicates on farmland, a gigatonne-scale industry whose central
uncertainty is exactly the dissolution kinetics this program computes.
A side product already exists: an agent extracted the field's standard
rate compilation (Palandri & Kharaka 2004, a PDF-bound USGS report)
into a machine-readable, provenance-tracked database — 36 minerals, 74
rate laws — as groundwork for a modern successor.

## 5. Flagship one: why the field is 1000× slower than the lab

The great embarrassment of weathering kinetics: minerals dissolve
**two to four orders of magnitude slower in nature** than in laboratory
experiments, and after four decades nobody has quantitatively
attributed the gap. (A 2025 review by Susan Brantley re-opened the
problem; its first-listed unresolved cause — "changing mineral surface
attributes difficult to model" — is this project's target painted on a
barn.) The leading intrinsic explanation is that fresh, ground,
defect-rich lab surfaces are simply more reactive than ancient field
surfaces whose reactive sites were spent long ago — a hypothesis about
the *statistical state of the lattice*, which continuum rate equations
structurally cannot represent and a defect-resolved KMC tracks
natively.

The plan (fully scoped in [scoping/field-lab-discrepancy.md](scoping/field-lab-discrepancy.md)):
run petra from "freshly ground" to "geologically aged," watch the rate
decline as the defect inventory depletes, and produce the first
quantitative budget of how many of those missing orders of magnitude
are intrinsic surface aging versus everything else. Validation data is
all published — including a famous null result from Lasaga's own lab
(dislocation density doesn't matter far from equilibrium) that the
model must reproduce to be believed. If it works, it says something
about every soil profile, river flux, and ERW carbon credit on Earth.

## 6. The double Lasaga: one engine, other worlds

Antonio Lasaga's signature move was importing physical chemistry's
tools into geochemistry. The generalized petra enables the same move
*outward*: a deck is really a declarative spatial stochastic state
machine, and the KMC update rule is just one strategy. Swap it for
synchronous updating and a deck is Conway's Game of Life; for
Metropolis dynamics and it's the Ising model of magnetism; for
per-step probabilities and it's an epidemic. Those three — each with
an exact analytic answer — are the engine's new conformance tests:
when they pass, "general" is proven, not claimed.

The wedge into other fields is a quiet scandal: most published lattice
models in ecology, social science, and materials degradation are
single-run bespoke scripts — no reproducibility, no statistics, no
shared format. Petra's offer is models-as-documents, deterministic
replay, and ten-thousand-run ensembles with real error bars. Two
showcases are scoped:

**Pitting corrosion** ([scoping/pitting-corrosion.md](scoping/pitting-corrosion.md)).
The localized corrosion of stainless steel — a trillion-dollar-a-year
problem — is mathematically the dissertation's problem wearing a hard
hat: passive-film breakdown at defects, repair racing attack,
stochastic pit initiation. The field models it with mean-field theory
or ad hoc cellular automata; nobody does exact stochastic dynamics
with ensemble statistics matched to the experimental pit-initiation
distributions. The headline target: reproduce a famous 2004 *Science*
result showing pitting onset is a cooperative critical phenomenon.

**Astrochemistry** ([scoping/astrochemistry.md](scoping/astrochemistry.md)).
The chemistry of interstellar clouds — how molecules like methanol
form on ice-coated dust grains at 10 K — runs on reaction-rate
databases in which roughly **60% of the rates are educated guesses**,
and on grain-surface KMC codes that are all closed and single-group.
At 10 K chemistry proceeds by quantum tunneling, which quarry treats
properly; and JWST's new ice observations have made the models the
bottleneck. First blood was drawn this week: a de-risking campaign
computed five benchmark reactions, matched the gas-phase literature
within honest error bars, exposed exactly why cheap methods can't yet
claim *surface* rates (a 10⁴ gap traceable to missing surface
physics), and rejected two of its own transition states with receipts
rather than publish them. That last sentence is the scientific
culture of this program in miniature.

## 7. The perfect circle: argon in muscovite

Before the dissertation, there was a senior thesis. Sletten & Onstott
(*GCA*, 1998) asked what ⁴⁰Ar/³⁹Ar step-heating of muscovite — a
cornerstone technique of geochronology — actually measures, and
answered with experiments: when you heat muscovite in a vacuum, it
doesn't politely diffuse its argon out; it **falls apart** —
dehydroxylates and delaminates — and the argon release tracks the
structural breakdown, not lattice diffusion. The radiogenic argon
lives in interlayer vacancies; atmospheric argon hides in extended
defects; the elegant "diffusion profiles" read from age spectra may
be artifacts of decomposition. Conclusion: recovering true thermal
histories from muscovite step-heating, as classically interpreted,
is not possible.

Twenty-four years later (2022), the first atomistic simulations of Ar
in muscovite computed the ideal-lattice migration barrier — ~66
kcal/mol, divacancy mechanism — and concluded that pristine-lattice
diffusion is *too slow to matter*. The 1998 empirical argument and
the 2022 ideal-lattice theory agree from opposite directions, and
between them sits a gap nobody has computed: **a model of the
structurally-evolving, defect-bearing lattice itself** — dehydroxylation
propagating, layers delaminating, vacancies migrating, argon riding
the chaos out of the crystal.

That model is, feature for feature, what petra now is: interlayer
sites and vacancies, extended defects, multilayer structure,
competing reactions, and (one small engine feature away) programmable
temperature schedules. The experiment that closes the loop: generate
**synthetic ⁴⁰Ar/³⁹Ar age spectra** from a structurally-evolving
muscovite deck, seeded with the 2022 published barriers, and test
them against the 1998 paper's own release curves — the rare study
whose hypothesis and validation data share an author, twenty-seven
years apart. The dissertation became the engine; the senior thesis
becomes its flagship application.

## 8. How it actually runs

The program's day-to-day labor is done by a small fleet of AI agents
— cloud workers, local workers, and a pair of otherwise-idle laptops —
coordinating through nothing but this git repository. Work lives on a
card board (`docs/program/`); an agent claims a card by pushing a
branch (atomic, race-free, by git's own rules); the work, the status
update, and the record all land in one reviewed pull request; a cloud
agent handles the merge ceremony. There is no server, no scheduler, no
meeting. In the program's first thirty hours the board processed the
quantum pipeline's construction, the engine generalization, the astro
campaign, and the kinetics database — across three different agent
runtimes, two of which found and honestly reported their own failures
before a human ever looked.

Every number carries provenance. Every engine change must reproduce
the past bit-for-bit. Every claimed transition state must survive
structural gates that exist because a plausible wrong answer fooled
us once. And every lesson — scientific or operational — is written
down in the repo where the next worker, human or otherwise, cold-starts
from it.

## 9. The point

Strip away the machinery and the program is one idea, held for
twenty-seven years: that the macroscopic behavior of minerals —
how they dissolve, corrode, and leak their argon — is the statistics
of countable atomic events, and that if you can price each event with
quantum mechanics and let a lattice play them out, the geology
computes itself. In 1999 that was a dissertation's reach exceeding
its grasp. In 2026 it's a barrier per hour on a machine whose coolant
runs at 45 °C, an engine that treats minerals, steel, and stardust as
the same kind of object, and a to-do list that ends with a 1998
paper finally getting its sequel.

*— written 2026-08-20, while the Si–O–Si barrier ladder ran on the GPU
in the background, which is the whole point.*
