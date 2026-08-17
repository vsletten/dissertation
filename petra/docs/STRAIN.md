# Defects and strain energy in a repeated lattice

The puzzle, as posed: kinetic Monte Carlo on a lattice gets its tractability
from the *repeated* lattice — every site is a copy of a template, topology is
fixed, energetics are local. Crystallographic defects violate all three at
once: a dislocation is a place where the tiling fails to close, and it carries
a long-range elastic strain field that modifies reaction energetics far from
the core. How can a model built on perfect repetition host the very
imperfections that dominate real mineral reactivity (etch pits nucleate at
dislocation outcrops; growth spirals wind around screw cores)?

This document works the problem out. The resolution is that "defects" are
not one problem but **three, of sharply different difficulty**, and the
repeated lattice blocks only part of one of them:

1. **Energetics** — the strain field changes reaction rates. This does NOT
   require breaking the lattice at all: strain energy enters KMC as a
   *per-site scalar field* modifying barriers. Linear elasticity gives the
   field analytically. **Implemented** (§2, §5).
2. **Topology** — the defect changes connectivity itself. For screw
   dislocations this is *exactly representable* on an unbroken lattice by
   re-wiring bonds across a cut surface — the trick behind every
   Burton–Cabrera–Frank spiral-growth KMC. **Designed, not yet implemented**
   (§3).
3. **Self-consistency** — strain that relaxes as material dissolves, moves
   with the dislocation, responds to the evolving surface. This genuinely
   needs a mechanics solver coupled to KMC. **Deferred, with the seam
   designed so it can land without touching the reaction language** (§4).

The classic geochemistry literature (Lasaga & Blum 1986 on dissolution at
dislocations; Frank's hollow core; Cabrera–Levine; BCF spiral growth;
Lasaga & Lüttge's stepwave model) lives almost entirely in category 1, which
is why category 1 is where we start: it is both tractable *and* the part
with established, testable physics.

---

## 1. Decomposing the puzzle

A defect affects a KMC model through exactly two channels:

- **the rate of every event near it** (strain energy stored in the solid
  raises the free energy of strained material → detachment is easier,
  attachment is harder), and
- **the set of events that exist at all** (a spiral step emanating from a
  screw core never grows away — new step is always exposed — because the
  *connectivity* winds around the core).

The first channel needs no geometry: rates in Petra are already
`base × exp(−(Ea + ΣΔEa)/RT)`, and a strain contribution is just one more
ΔEa term — one that depends on *where the site is* rather than on its
neighbors' states. The lattice stays perfectly repeated; only the rate
field over it is not.

The second channel needs graph surgery, not geometry: KMC never uses atomic
positions for dynamics (positions are output decoration in Petra, exactly as
in the legacy model). A dislocation's topological content is which sites
bond to which — and that is data Petra already owns per-bond.

Positions only lie to the *viewer* (atoms drawn at ideal sites are slightly
wrong near a core). That is a rendering refinement, not a physics problem.

## 2. The energetic route: static analytic strain fields

### 2.1 Physics

Isotropic linear elasticity gives the strain energy density around a
dislocation line. Integrated over one site's atomic volume Ω, the energy
stored *in the material at that site*, at perpendicular distance r from the
line:

- **screw**: `u(r) = μ b² Ω / (8 π² r²)`
- **edge** (θ-averaged, the Blum–Lasaga form): `u(r) = μ b² Ω / (8 π² (1−ν) r²)`

with μ the shear modulus, b the Burgers vector magnitude, ν Poisson's ratio.
Multiple defects superpose (linear elasticity). Two singular-core corrections,
both physical:

- **core clamp**: continuum elasticity fails inside r ≲ b; we clamp
  `r → max(r, r_core)` with `r_core` a deck parameter (~b by default);
- **hollow core** (Frank 1951): where u would exceed the cost of creating
  new surface, the core is energetically hollow; the clamp doubles as the
  cap, and an explicit `cap` (kcal/mol) is available.

This is precisely the energetics behind etch-pit formation: dissolution
accelerates as exp(u/RT) toward the core, and whether a macroscopic pit
opens (vs. mere point-depression) is the Cabrera–Levine criterion —
strain energy vs. surface energy at the critical undersaturation. Our KMC
does not *assume* that criterion; it should *reproduce* pit opening from
the site-level rates, which is the point of doing this with KMC at all.

### 2.2 Thermodynamic bookkeeping (the part that's easy to get wrong)

Strain energy is a **state function of the solid**: material sitting at
site i carries extra free energy u_i. So for any reaction that removes
material from site i (dissolution, hydrolysis of its bonds) vs. its reverse:

```
ΔG_i = ΔG_0 − u_i          (the reaction is more downhill by u_i)
```

How that splits between forward and reverse barriers is a transfer
coefficient β ∈ [0,1] (Brønsted–Evans–Polanyi):

```
Ea_forward(i) = Ea_f0 − β u_i
Ea_reverse(i) = Ea_r0 + (1−β) u_i
```

Any β gives correct thermodynamics (the *difference* is u_i); β sets how
much of the destabilization is felt in the transition state. Blum & Lasaga's
classic treatment is β = 1 (all of it in the dissolution rate,
`k = k_0 exp(u/RT)`). Petra makes β a per-reaction deck choice: each
reaction declares `strain_scale`, the multiplier on u_i added to its Ea —
so a forward/reverse pair writes `strain_scale = -β` and `+(1−β)`
respectively, and the detailed-balance linter (P2 roadmap) can eventually
check that pairs sum to −1.

**Convention**: `u_i ≥ 0` is stored elastic energy (kcal/mol, in the deck's
energy units) at site i; `Ea_eff = Ea + strain_scale · u_i`.

### 2.3 What "static" honestly buys and costs

The field is computed once at build time from the deck's defect list and
never changes. That is the standard approximation of the dissolution
literature, and it is *appropriate* when: dislocations are immobile at the
temperature of interest (true for silicates far below annealing
temperatures), dissolution does not appreciably relax the remaining
lattice's long-range field (true until the pit has consumed a substantial
neighborhood of the core), and defect density is low enough that fields
superpose linearly. It is wrong for: dislocation motion/multiplication,
strain relaxation at the opening pit walls, misfit strain from growth.
Those are §4's problem. The engine cost of a static field is zero — it is
a per-site constant folded into rate resolution, fully compatible with the
incremental event updates.

## 3. The topological route: defects as graph surgery

### 3.1 Screw dislocations — exact on a repeated lattice

A screw dislocation with Burgers vector **b** equal to a full lattice
translation along its line is *pure topology*: take the perfect lattice,
choose a half-plane H bounded by the line, and for every bond crossing H,
re-target it to the equivalent site displaced by **b**. Every site keeps its
template, coordination, and states; only the adjacency winds. Going around
the core once now advances you by one cell along the line — the lattice has
become a helicoid *in its connectivity*, which is the only sense dynamics
ever sees. This is exactly how BCF spiral growth is realized in
solid-on-solid KMC, and it composes with §2's strain field (the topology
supplies the perpetual step; the field supplies the core energetics).

Petra mechanics (designed): a `[[defects]] topology = true` entry triggers a
post-build **rewiring pass**: for each adjacency entry (i→j) whose bond
crosses the cut half-plane (a geometric test on the two cell coordinates
against the declared line + cut direction), re-target j → j ± b_cells along
the line axis (min-image on periodic axes; drop at open boundaries). Bond
labels ride along untouched; CSR is rebuilt once. The `seam`-style
convention extends to mark rewired bonds for the viewer. Constraints to
enforce at compile: **b** must be an integer number of cells along the line
axis (partial Burgers vectors ⇒ stacking fault, below), and the line axis
must be periodic or the line must terminate on open faces.

Caveat for the current kaolinite deck: it is a single (001) sheet (the
legacy model never tiled c), so a screw threading the basal surface has no
lattice to wind through — screw topology becomes interesting exactly when a
deck tiles the third dimension. The engine supports 3D tiling already; the
kaolinite deck's 2D-ness is inherited from the legacy model, not a Petra
limit.

### 3.2 Edge dislocations, stacking faults, twins — harder, in order

- **Stacking fault / partial-Burgers re-registry** (includes twins): two
  half-crystals with a shifted registry across a plane. On-lattice: bonds
  crossing the fault plane use *different cell offsets* than the bulk
  template. Representable today, clumsily, by hand-authoring; cleanly, by a
  region-scoped bond-template override in the deck (a bounded extension of
  the rewiring pass). The fault also carries a constant *areal* energy —
  a per-site field stamp on the fault plane (γ_SF × area/site), which is
  §2 machinery.
- **Edge dislocations**: an inserted/removed half-plane changes the *site
  count*, not just wiring — the honest on-lattice representation removes a
  half-plane of sites and reconnects across the gap, leaving a row of
  under-coordinated core sites. Feasible as graph surgery (delete sites,
  stitch adjacency), but there is no unique canonical stitching; it needs a
  per-crystal-family recipe. The *energetics* of an edge (which drive etch
  pits) need none of this — §2's edge field works on the perfect lattice.
  Recommended posture: edge cores enter first through §2 + optionally a
  hand-authored core state via init passes; true edge surgery waits for a
  motivating question.
- **Grain boundaries**: out of scope for a single-lattice engine; two
  lattices glued by an interface graph is a different (interesting) project.

## 4. The self-consistent frontier — and the seam that keeps it open

Real strain relaxes: as the etch pit opens, the field that drove it decays.
Three escalation levels, all published in some KMC literature:

1. **Quasi-static re-evaluation**: recompute the analytic field every N
   events with defect *geometry* updated (e.g. the dissolved region excluded
   by image constructions). Cheap, crude.
2. **Harmonic lattice relaxation**: give the lattice ball-and-spring
   (Keating-type) elasticity; after events near defects, minimize locally
   (sparse linear solve over a relaxation radius); u_i falls out of the
   spring energies. This is the natural on-lattice mechanics — no FEM mesh,
   and it reuses the site graph. The realistic next step if relaxation ever
   matters.
3. **Coupled KMC–FEM**: full continuum solve on the evolving geometry
   (the strained-epitaxy growth literature). Heavy; only if quantitative
   long-range relaxation becomes the science question.

**The architectural decision that matters today**: rates couple to strain
only through the per-site field `u` and the per-reaction `strain_scale`.
The field's *provider* — analytic formulas now, a relaxation solver later —
is invisible to the reaction language, the deck reactions, and the engine's
incremental machinery (a solver that updates u_i for a region just marks
those sites dirty, which the engine already knows how to absorb). Nothing
authored today is thrown away when the provider upgrades.

## 5. As implemented (this change)

- `Lattice.strain: Vec<f64>` — per-site stored elastic energy, kcal/mol
  internal (deck units converted at compile), zeros unless the deck declares
  defects.
- Deck section:

```toml
[[defects]]
type = "screw"              # or "edge"
line_axis = 2               # dislocation line direction: 0=a, 1=b, 2=c
at = [10.0, 1.5, 0.0]       # a point the line passes through, in CELL coords
                            # (the component along line_axis is irrelevant)
burgers = 7.4               # |b|, Å
shear_modulus = 20.0        # μ, GPa
poisson = 0.25              # ν, edge only
core_radius = 5.0           # Å; r is clamped to this (continuum cutoff /
                            # Frank hollow core), default = burgers
# strain_prefactor = 30.0   # alternative: set A directly (deck-energy·Å²),
                            # u = A/r², overriding the physical inputs
```

  `u_i = A / max(r_i, r_core)²`, `A = μ b² Ω / (8π²)` (screw; edge divides
  by (1−ν)), Ω = cell volume / sites-per-cell, GPa·Å³ → energy via
  0.1439 kcal/mol. `r_i` is the perpendicular distance from site i's
  Cartesian position to the line, **minimum-image on periodic axes**.
  Multiple `[[defects]]` superpose.
- Rate coupling: per-reaction `strain = { scale = s }` →
  `Ea_eff = Ea + s·u_center` (§2.2 conventions; dissolution-forward
  reactions use negative s).
- Export: snapshots gain a `strain` f32 node column when any defect exists —
  the viewer's inspector shows it per site, and it is available for future
  color-by-strain rendering.
- Gates: (a) field values against hand-computed formulas, including clamp
  and min-image; (b) *exactness* of the rate coupling
  (`k(r)/k(∞) = exp(β u(r)/RT)` to f64 precision); (c) an ensemble
  etch-pit gate — Kossel crystal with a through-screw dissolves
  preferentially near the core (near-core dissolved fraction exceeds
  far-field at matched depth).
- Demo: `examples/kossel-etchpit.toml` — watchable in graph-viz (recorded
  or live): a pit drilling down the dislocation line.

## 6. Open questions (Victor)

1. **Which defects matter in kaolinite?** Kaolinite plates are famously
   defect-rich in *stacking* (the b/3 interlayer shifts) rather than
   classic screw/edge lines. If the target science is kaolinite etch
   pitting, the stacking-fault re-registry (§3.2) plus its areal energy
   stamp may matter more than line-defect fields — and it needs a c-tiled
   kaolinite deck (the current one is the legacy single sheet).
2. **Elastic inputs**: μ ≈ 20 GPa-ish for kaolinite books, but anisotropy
   is large in layered silicates. Is isotropic-μ good enough for the
   questions you'd ask, or does the field need the anisotropic form?
3. **β (transfer coefficient)**: default the decks to β = 1
   (Blum–Lasaga: all strain in the forward rate) or β = ½? It changes
   kinetics near the core, not thermodynamics.
4. **Hollow-core handling**: clamp at r_core (implemented) vs. an explicit
   removed-core (init pass emptying sites within r_h) — the latter is more
   faithful to Frank for large-b dislocations and is already expressible
   with existing init passes; worth a demo?
