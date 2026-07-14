# 02 — Faithful Model Spec (behavior as implemented)

The physics and the KMC algorithm **exactly as the C++ code executes them**, so the
Rust port can reproduce the same numbers before anyone "improves" anything. Two
strictly separated sections:

- **Part A — Intended behavior**: the model as designed, faithfully described.
- **Part B — Warts, bugs, dead code**: places where the code does something other
  than what the docs/comments/README claim, or something clearly accidental. These
  are flagged, not fixed. The port should decide case-by-case whether to preserve
  bug-for-bug (default for reproducibility) or repair (only with Victor's sign-off).

File references are to `dissertation/main/model/`.

---

# PART A — Intended behavior

## A1. The system being modeled

A **kaolinite** surface — a 1:1 Al-Si clay: one sheet of Si-O tetrahedra bonded to
one sheet of Al-O(H) octahedra. The simulation models the atomistic surface reactions
of dissolution and (re)growth in water: hydrolysis (water breaking metal-oxygen-metal
bridges), and adsorption/desorption of dissolved Si(OH)₄ and Al(OH,H₂O)₆ species.
The scientific question the model was built for: **how does the chemical state of a
site's nth-nearest neighbors change the rate at which that site reacts?** — i.e.
neighbor-dependent kinetics, the thing that distinguishes a real surface from a
mean-field rate law.

## A2. Lattice as a graph

The crystal is represented as a fixed graph:

- A **unit cell** (motif) of `Npos` positions (26 in the sample `data.cell`), each
  with a class (Al / Si / one of three oxygen types) and up to 6 neighbors given as
  (target-position, Δa, Δb, Δc) offsets.
- The motif is tiled `Num_aCells × Num_bCells` times into a flat array of
  `Num_aCells · Num_bCells · Npos` **lattice sites**. Site linear index
  `i = a·(Num_bCells·Npos) + b·Npos + n`.
- Each site's 6 neighbors are resolved once at build time into integer indices
  (`nbr[j]`), or `-1` where a neighbor is absent (open boundary or a template slot
  with target `n < 0`). Degree is fixed per class but ≤ 6 in general.

So the lattice **is** an adjacency-list graph with fixed topology and a mutable
integer `state` per node. Chemistry never changes topology (except the dead cluster
removal); it only rewrites `state` labels on nodes.

### A2.1 Boundary conditions (`Lattice::GetNeighbor`)
A **slab**: one direction is the surface normal (open), the two in-plane directions
wrap (periodic), governed by `SurfacePlane`:

- `SurfacePlane == 0` (ac surface): the **b** direction is open (a neighbor with
  `b == Num_bCells` or `b < 0` returns `-1`); the **a** direction is periodic
  (wraps modulo `Num_aCells`).
- `SurfacePlane == 1` (bc surface): the **a** direction is open; **b** periodic.

Concretely: compute neighbor cell `(a+Δa, b+Δb)`. If the template slot's target
`n < 0`, no neighbor. Else if the open direction went out of range, `-1`. Else wrap
the periodic direction into range and compute the linear index. The third crystal
axis (`Δc`) is carried in the data but not used in resolution — this is a
single-sheet 2D tiling.

### A2.2 The `EDGE` sentinel and frozen boundary
`TerminateLattice` marks the topmost and bottommost rows of the open direction as
`EDGE` (state code 9) when a site there doesn't have its full expected neighbor count.
`EDGE` sites are frozen: they never react, are skipped in output, and every
environment check bails out if it "runs into" an EDGE. They are the model's fixed
boundary — the bulk crystal below the simulated surface.

## A3. State encoding (`common.hpp`)

`state` is a 3-digit integer. First digit = class; last two = occupancy /
coordination / protonation. Empty = `state % 100 == 0`.

| Range | Class | Meaning |
|---|---|---|
| 100 | Al | empty Al site |
| 101–107 | Al | Al with 0–6 OH/H₂O ligands |
| 199 | Al-class | (intended) Si adsorbed on an Al site — **never produced, see B**|
| 200 | Si | empty Si site |
| 201–205 | Si | Si with 0–4 OH groups |
| 299 | Si-class | (intended) Al adsorbed on a Si site — **never produced, see B** |
| 300–303 | O | Si-O-Si bridging oxygen (empty / bridge / both-OH / one-OH) |
| 400–410 | O | Si-O-Al₂ bridging oxygen — the richest class, 10 live states |
| 500–503 | O | Al-OH-Al bridging oxygen |
| 9 (`EDGE`) | — | frozen boundary sentinel |
| x99 (`WRONG`) | — | wrong-cation marker (dead) |

Helper macros (`lattice.hpp`): `TYPE(A)=state/100`, `ISOCC(A)=state%100>0`,
`ISSI`, `ISAL`, `ISOX(state>299)`, `ISEDGE(state==9)`.

The full semantic gloss of every 400-series state (which of the two Al is present,
which O-H is protonated) is in the `common.hpp` header comment and should be copied
verbatim into the Rust port as documentation — it is the authoritative key.

## A4. Building the initial configuration

Order (from `main`): `FindPairs → PopulateSolid → TerminateSurface → TerminateLattice`.

### A4.1 `FindPairs` (`lattice.cpp`)
Each 400/500 bridging oxygen belongs to a **double bridge** (two oxygens bridging the
same pair of Al). This links each such O to its partner via `pair`. Algorithm: for
each O with `state ≥ 400` and no partner yet, take its first two neighbors (`al1`,
`al2`); search `al1`'s neighbors for another O (`o2`) whose first two neighbors are
also `{al1, al2}`; link them. If a boundary is hit first (`nbr < 0`), it gives up on
that O (leaves `pair = -1`). `pair` is later essential to several environment checks.

### A4.2 `PopulateSolid(dmSi, dmAl)` (`lattice.cpp`)
Decides the **fill fraction** (how deep the initial solid slab is) from the sum of
chemical-potential deltas:

- `dmSi + dmAl > 0.5` → supersaturated → fill fraction **0.3**
- `dmSi + dmAl < -0.5` → undersaturated → fill fraction **0.7**
- otherwise → near equilibrium → **0.5**

Then it fills the first `frac ·` (depth) unit-cell layers by **incrementing** each
site's state by 1 (100→101, 200→201, 300→301, 400→401, 500→501), i.e. promoting
empty→first-occupied. Special guard: an empty Si-O-Si oxygen (state 300) is only
promoted to 301 if its full connectivity chain exists (both Si neighbors present, and
each of those Si has its first neighbor present); otherwise it is left empty. (The
"depth" axis is whichever direction is open; the in-plane axis is filled fully.)

### A4.3 `TerminateSurface` (`lattice.cpp`)
Caps dangling bonds where the solid meets water, in three passes:
1. **Dangle the O's missing a cation**: any occupied bridging O adjacent to an empty
   cation site is demoted along a fixed state map (e.g. 401→404 or 406 depending on
   whether the missing cation was Si or Al; 301→303; 501→503; etc.).
2. **Turn on dangling oxygens to terminate bonds**: for occupied cations whose
   coordination is minimal (`state % 100 == 1`), promote adjacent empty/partial
   oxygens to their OH-terminated forms (300→303, 400→408/409, 500→503, …).
3. **Update cation states to reflect terminal OH's**: increment Al/Si state for each
   adjacent terminal-OH oxygen (503/409 for Al; 303/408 for Si).

The net effect is a chemically plausible hydroxylated surface as the starting point.

### A4.4 `TerminateLattice` (`lattice.cpp`)
Two passes: (1) set boundary rows (first/last along the open axis) to `EDGE` where the
actual neighbor count ≠ the expected count for that class; (2) revert any 301 Si-O-Si
bridge whose connectivity chain now touches an EDGE back to empty (300). This freezes
the crystal's cut faces.

## A5. Reactions

28 reactions, indices 0–27, loaded from `data.rxn`.

### A5.1 Hydrolysis (0–15): 8 forward/reverse pairs
Even index = forward (bond-breaking, +H₂O), odd index = its reverse (bond-forming).
The pairs, with the reactant→product state change on the oxygen site (full list in
`common.hpp`):

| Fwd | Rev | O transition (fwd) | Chemistry |
|---|---|---|---|
| R0 | R1 | 301→302 | Si-O-Si + H₂O ⇌ 2 Si-OH |
| R2 | R3 | 401→402 | Si-O-Al₂ + H₂O ⇌ Si-OH + Al-OH-Al |
| R4 | R5 | 401→410 | Si-O-Al₂ + H₂O ⇌ Si-OH-Al + Al-OH |
| R6 | R7 | 402→403 | second-stage Al-O-Al hydrolysis |
| R8 | R9 | 410→403 | |
| R10 | R11 | 406→407 | Si-OH-Al + H₂O ⇌ Si-OH + Al-H₂O |
| R12 | R13 | 404→405 | Al-OH-Al (400-side) + H₂O ⇌ Al-OH + Al-H₂O |
| R14 | R15 | 501→502 | Al-OH-Al (500-side) + H₂O ⇌ Al-OH + Al-H₂O |

Each forward reaction mutates the oxygen site's state and increments/decrements the
coordination of its participating cation neighbor(s); the reverse does the inverse.
The exact neighbor bookkeeping — including the stochastic choice of *which* of two
symmetric Al takes the proton (R4, R9 use a `ran2() < 0.5` coin, recording the choice
in `lostal` so the reverse can undo it) — is implemented one function per reaction in
`actions.cpp` (`DoReaction0`…`DoReaction15`). These must be ported literally; they
encode the model's mechanism.

### A5.2 Adsorption (16–19) and desorption (20–23)
- **AdsorbAl** (16/17): an empty Al site becomes Al(OH,H₂O)₆ (107); its 6 neighboring
  oxygens each advance one protonation step along a fixed map (500→503→502, 400→409,
  408→407, 406→410 with `lostal` recorded, etc.). **AdsorbSi** (18/19): empty Si→
  Si(OH)₄ (205); its 4 oxygens advance (300→303→302, 404→402, 405→403, 409→407,
  400→408).
- **DesorbAl** (20/23) and **DesorbSi** (21/22) are the exact inverses, restoring the
  cation to empty (100/200) and walking the neighboring oxygens back down.
- The reactant states are: adsorb Al reactant 100 (rxn 16), adsorb Si reactant 200
  (rxn 19); desorb Al reactant 107 (rxn 20), desorb Si reactant 205 (rxn 22). The
  "cross" variants (17, 18, 21, 23: wrong cation on wrong site, states 199/299) are
  present in the tables but gated off — see Part B.

### A5.3 Diffusion (24–27)
Intended as desorb-here + adsorb-there (a surface hop), with the event's `rxn` field
carrying `OFFSET + target`. **Never fires** in practice — see Part B.

### A5.4 Rate construction (`rxnlist.cpp`)
Input for each reaction is a set of `(pre-exponential k, activation energy ΔE)` pairs,
one per environment bucket. Constants: `R = 1.987e-3 kcal/mol/K`, and `rt = R·T`.

- **Hydrolysis forward**: `rate = k` (the pre-exponential *as given* — no Arrhenius
  factor applied; see B for why this is notable).
- **Hydrolysis reverse**: `rate = k · exp(ΔE / rt)`. So forward/reverse ratio encodes
  the equilibrium constant `exp(−ΔE/RT)`.
- **Adsorption**: `rate = a · exp(ΔE/rt) · exp(Δμ/rt)`, where Δμ is `dmsi` for the Si
  adsorption reactions and `dmal` for Al — coupling adsorption rate to solution
  saturation. (Index test in code: `i < 18` uses `dmsi`, else `dmal`.)
- **Desorption**: `rate = a · exp(−ΔE/rt)` per bucket.
- **Diffusion**: copied from the corresponding desorption reaction's rates.

In the shipped `data.rxn`, every bucket within a reaction has **identical** (k, ΔE),
so the per-environment rate variation is currently flat (see B4 — this neutralizes the
nth-neighbor effect in the sample data, even though the machinery to express it is
fully present).

## A6. The nth-nearest-neighbor environment (`envrn.cpp`)

This is the scientific heart. For each candidate (site, reaction), two questions:

### A6.1 `IsActive(site, rxn)` — is the reaction allowed here?
- **Forward hydrolysis** (even, <16): allowed only if the site is "at the surface" —
  defined operationally as having a **2nd-neighbor** (neighbor-of-neighbor) in a
  hydrolyzed/OH-bearing state (303, 404, 405, 406, 408, or 409). This double loop over
  `nbr[i].nbr[j]` is the concrete "nth-nearest-neighbor" reach: a reaction's
  legality depends on the chemical state two hops out.
- **Reverse hydrolysis** (odd, <16): always allowed.
- **Adsorption** (rxn 16 or 19 only): allowed if the site has ≥1 occupied,
  non-EDGE neighbor (can't nucleate on nothing).
- **Desorption** (rxn 20 or 22 only): always allowed.
- **Everything else** (17, 18, 21, 23, and all diffusion 24–27): returns FALSE →
  disabled. This is how the cross-cation and diffusion reactions are switched off.

### A6.2 `CheckEnv(site)` — which environment bucket → which rate?
Dispatches by class to `Check100/200/300/400/500`. Each counts specific
neighbor/2nd-neighbor states and folds them into a single integer index into the
reaction's `rate[]` array. Intended index formulas are documented in `common.hpp`
(e.g. 300s: `index = 5·x + y` with x = # broken 400s among 2nd-neighbors, y = # broken
300s). The 300/400/500 checkers reach through `pair` and through Si→O second neighbors,
so the bucket genuinely reflects nth-neighbor chemical state. (Note: the 100/200
checkers' returned index does **not** match the documented formula — see B3.)

If a checker "runs into" an EDGE or a missing neighbor, it prints a diagnostic and
returns `-1`, which the event builder treats as a fatal invalid-environment error and
aborts the run. So a correctly terminated lattice is a correctness precondition.

## A7. The KMC step (`actions.cpp DoEvent`)

Standard rejection-free (n-fold-way / Gillespie direct) selection:
1. `ratesum = Σ event.rate` over the whole event list.
2. If `ratesum == 0` → error, stop.
3. **Time increment**: `dt = −ln(ran2()) / ratesum` (exponential waiting time for the
   aggregate Poisson process).
4. **Selection**: draw `eps = ran2()`; walk the list accumulating `rate/ratesum`;
   the first event whose cumulative sum ≥ eps is chosen.
5. Execute it (`DoReaction`), which rewrites the site and its affected neighbors.
6. Return `dt`; `main` accumulates `time += dt`.

The event list is regenerated from scratch on the next iteration — every rate is
recomputed, so state changes propagate implicitly. Simulation ends after `nsteps`
iterations (fixed count), or early if no events are possible or a reaction fails.

## A8. Time-step handling

- `dt` is a `float`, and `time` (in `main`) is a `float` accumulator over up to
  millions of steps. Waiting time is the KMC-canonical `−ln(u)/R_total`.
- Output snapshots are triggered on **step count**, not simulated time: `writeData`
  every `wsteps` steps, `writeMSI` every `msteps` steps (with `msteps == -1`
  reinterpreted as "final only" via `msteps = nsteps + 1`).
- There is no wall-clock, no real-time coupling; simulated `time` is only reported in
  the output rows.

## A9. Output formats

### A9.1 Cerius2 MSI (`writeMSI`) — the viewer's input
Text, parenthesized object tree, one `Model` containing `Atom` and (optionally)
`Bond` objects:
```
# MSI CERIUS2 DataModel File Version 3 5
(1 Model
(A I Id 1)
 (A C Label "start")
 (2 Atom
  (A C ACL "13 Al")            ; element as "Z Symbol"
  (A D XYZ (1.17181 4.38134 2.84548))
  (A I Id 2)
  (A C Label "Al0")            ; class label + site index
  (A I LabelType 0)
 )
 …
 (9002 Bond
  (A O Atom1 10)               ; references atom object IDs
  (A O Atom2 2)
 )
)
```
Rules the writer follows:
- Object IDs start at 2 and increment across *all* atoms then *all* bonds
  (`nthing`), so bonds and atoms share one ID space; atom IDs are remembered in an
  `id[]` map for bond emission.
- Only occupied, non-EDGE sites are emitted.
- Element mapping by class: 1→"13 Al", 2→"14 Si"; oxygen classes 3/4 → "1 H" if
  `state%100 > 1` (protonated) else "8 O"; class 5 → "1 H" or **"9 F"** (fluorine is
  used as a visual stand-in for the Al-OH-Al bridge oxygen — a rendering hack, not
  chemistry). Class 0 → "10 Ne" (edge, shouldn't occur since EDGE is skipped).
- Bonds are emitted only when `drawbonds` is set; each undirected neighbor pair once
  (`i2 < i`), and **periodic wrap bonds are suppressed** (an explicit check against
  `a==0/aCells-1` etc. with the template offset sign) so the rendered structure
  doesn't draw bonds across the periodic seam.
- Coordinates are Cartesian, computed from the unit-cell fractional coords and a
  **hard-coded** cell matrix `cd` in `output.cpp` (see B5).

### A9.2 XYZ (`writeXYZ`)
Comma-separated `Z,x,y,z` per occupied site (not whitespace-delimited standard XYZ,
and no header/count line — a nonstandard dialect). Al→13, Si→14, O→8.

### A9.3 Population time series (`writeData`)
One CSV row: `time, Si_total, Si(OH)0, …, Si(OH)4, Al_total, Al(OH,H2O)0, …,
Al(OH,H2O)6`. Counts derived from cation states. Written to a **separate file per
snapshot** (`step{i}.dat`, `end.dat`) — see B1.

### A9.4 Surface projections (`writeSurf`)
`surfSi.out` / `surfAl.out`: CSV `x,y` for cations that have ≥1 surface-exposed
oxygen neighbor (303, 404–409, or 503). A 2D map of the exposed surface.

---

# PART B — Warts, bugs, and dead code

Flagged, not fixed. Each notes the likely disposition for the port.

### B1. `results.dat` is documented but never written
`output::initDatafile()` does `remove("results.dat")` at startup, and **nothing ever
writes it** — `writeData` only ever targets `step{i}.dat` / `end.dat`, each of which
it *truncates* and writes a single row to. The README and `common`/CLAUDE docs
describe `results.dat` as the appended population time series; the committed
`results.dat` in the tree is a real multi-row series left over from an **older version
that appended to one file**. So today the "time series" is scattered across thousands
of one-row files (`step0.dat … step4999000.dat` for a 5M-step, wsteps=1000 run).
→ **Port disposition**: this is almost certainly not intended. Recommend the Rust
port write one appended time-series file (restore the original intent), but call it
out to Victor rather than silently matching the current per-file behavior.

### B2. The random seed in `data.sim` is silently ignored — runs are deterministic
`Simulation::CreateSimulation` reads fields in this order: `nsteps`, `wsteps`,
`msteps`, then **`drawbonds` twice** (two consecutive `f >> s->drawbonds`). The file
layout is `nsteps / wsteps / msteps / seed / drawbonds`. So the first `drawbonds`
read consumes the **seed line** (−2) and discards it; the second read consumes the
real `drawbonds`. `s->ranseed` is **never read from the file**; because `Simulation`
is heap-allocated with `new Simulation()` (value-initialization, no user
constructor), `ranseed` is zero-initialized to 0. `ran2` then treats seed ≤ 0 as
"seed me", mapping 0→1. **Net effect: every run uses the same PRNG seed regardless of
the file; results are bit-identical run to run** (there is exactly one intended
`Simulation` field read wrong). The `sim.hpp` comment even documents `ranseed` as
"seed for initran2", confirming the intent was to read it.
→ **Port disposition**: a real bug defeating the entire point of a seed field. The
Rust port should read the seed properly (parse all 5 fields), and can offer a
"legacy" flag to reproduce the fixed-seed behavior if bitwise comparison is ever
needed. Flag to Victor.

### B3. `Check100`/`Check200` environment indices don't match the documented scheme
`common.hpp` documents 100s as `index = x·3 + y` and 200s as `index = x·4 + y`.
The code returns `(x + y) / 2` for `Check100` and `(x + y)` for `Check200` — neither
matches. The counted quantities also differ subtly (e.g. `Check100` bumps `x` on
state 502 and `y` on 403/405/407/409/410, then integer-divides the sum by 2). This is
inconsistent with the header's own spec.
→ **Port disposition**: preserve as-is for bitwise reproduction, but annotate loudly.
It matters only if a `data.rxn` ever provides non-flat adsorption rate tables (it
doesn't today — see B4), so it's currently latent. Get Victor's read on what the
*intended* index was before "fixing".

### B4. The nth-neighbor rate machinery is currently inert (data, not code)
Every reaction's rate table in the shipped `data.rxn` has **identical (k, ΔE) across
all buckets** (all `1 2.6`, or all `100 4.8`, etc.). So although `CheckEnv` computes a
per-site environment bucket and indexes `rate[env]`, the selected rate is the same
regardless of environment. **The central scientific feature — neighbor-state-dependent
rates — has no effect on the dynamics in the sample input**; only `IsActive` gating
(the surface-reachability test) currently differentiates sites. This is not a code
bug — it's the state of the tuning data — but it is essential context: any port
"validated" only against this data cannot claim the environment machinery is correct,
because the data never exercises it. → **Port disposition**: build the machinery
faithfully; add a test `data.rxn` with distinct per-bucket rates to actually exercise
`CheckEnv`, and treat matching *those* outputs as the real validation.

### B5. Output geometry uses a hard-coded cell matrix, not `data.cell`
`writeMSI`/`writeXYZ`/`writeSurf` each define a local `float cd[3][3]` (≈
`{4.9725,-0.026,-1.336},{0,8.929,-0.301},{0,0,7.384}`) and derive `al,bl,cl` from it,
ignoring the `A/B/C/alpha/beta/gamma` actually read from `data.cell`
(`5.140 8.930 7.370`). The fractional coords come from `data.cell`, but the
Cartesian transform is baked in. So changing the cell dimensions in `data.cell`
changes the physics indexing but **not** the output coordinates. Three copies of the
same magic matrix. → **Port disposition**: obvious latent inconsistency; the port
should derive the Cartesian basis from the cell parameters in one place. Flag; likely
safe to fix, but confirm the magic matrix isn't a deliberate "known-good render cell".

### B6. Diffusion is entirely dead code (and the cluster/BFS machinery with it)
`IsActive` returns FALSE for every diffusion reaction (24–27), so diffusion events are
never added to the event list, so `DoReaction`'s default case (`Diffuse` +
`RemoveUnattachedClusters`) never runs, so `Lattice::RemoveUnattachedClusters`,
`BreadthFirstSearch::ColorNodes`, and the entire `bfsearch.cpp` queue implementation
are **unreachable at runtime**. The diffusion rate tables (rxn 24–27) are built but
never read. → **Port disposition**: do **not** port diffusion/BFS in the first
milestone — it's dead. Preserve the design note (diffusion = desorb+adsorb hop, floaters
removed by BFS) for a possible future milestone, but don't carry the dead weight.

### B7. Out-of-bounds read building diffusion rates (`rxnlist.cpp`)
In the diffusion loop, `numRates` still holds the value from the **last desorption**
read (5, from the 299 table), while `reactions[i-4].nrates` for rxn 24/25 is 4. The
loop allocates `new float[numRates]` (5) and copies `j < numRates` (5) entries from a
source array of length 4 → reads `reactions[20].rate[4]` / `[21].rate[4]` out of
bounds. Harmless in practice only because diffusion rates are never used (B6), but
it's a genuine buffer over-read. → **Port disposition**: won't exist in Rust (no
diffusion port); note it so nobody faithfully re-introduces it.

### B8. Everything is `float`; `-ffast-math` is on
`state` is `int`, but all rates, `dt`, and the `time` accumulator are `float`
(32-bit). Over millions of steps, `time += dt` loses precision (small `dt` added to a
large `time` can vanish); rate sums accumulate rounding. Production builds add
`-ffast-math`, which relaxes IEEE semantics (reassociation, no strict NaN handling),
so results can even depend on optimization level. → **Port disposition**: use `f64`
in Rust for `dt`/`time`/rates and Kahan-or-just-`f64` accumulation. This *changes the
numbers*, so it can't be validated by bitwise diff against the C++ output — validate
statistically (population trajectories, rates) instead. Flag as an intended
improvement, not a silent one.

### B9. `EventList : public Event` and `int`-as-bool / `-1`-as-null idioms
Public inheritance of a data class purely to bolt on a `next` pointer; `int` return
values used as booleans; `-1` sentinel neighbor indices everywhere; raw
`new`/`delete`/`malloc` mix (e.g. `writeMSI` uses `malloc`/`free` for `id[]` amid C++
`new` elsewhere). Not bugs — just C-in-C++ idioms that the Rust port replaces with
`Option<usize>`, real `bool`, `enum` states, and `Vec`. Listed so the port doesn't
mistake them for meaningful design.

### B10. Minor: doubled `return result;` and dead locals
`Environment::IsActive` ends with `return result; return result;` (second is dead).
`rxnlist.cpp` has commented-out `kt`/`Kb` machinery (a Boltzmann prefactor that was
considered and abandoned — consistent with A5.4's "forward rate = k with no Arrhenius
factor"). `Check400` computes `p` twice. Cosmetic; noted for completeness.

### B11. README/spec drift
The README claims "28 reaction types including … diffusion", "Gillespie first-reaction
method", and a `results.dat` time series — all three are imprecise: diffusion is dead
(B6), the selection is the **direct** method not first-reaction (A7), and `results.dat`
isn't written (B1). The README also lists a `data.cell` sample with 5 site types and a
neighbor block format that matches the code, but its `## Algorithm` step 5 ("Cleanup:
Remove unattached clusters") describes the dead BFS path as if it runs each step.
Trust the code (Part A), not the README, wherever they disagree.

---

## C. Validation strategy for the port (summary)

1. **Structural build** (`FindPairs`/`Populate`/`Terminate*`): deterministic, no RNG.
   The Rust port should reproduce the initial `state[]` array **exactly** for a given
   `data.cell`/`data.lattice`/`data.rxn`. This is the strongest, cleanest checkpoint —
   diff the start-of-run state array and the `start.msi` atom set.
2. **Dynamics**: because of B2 (fixed seed) *and* B8 (float), exact trajectory
   reproduction is fragile. Options: (a) port `ran2` faithfully and use `f32` to chase
   bitwise parity as a one-time confidence check; then (b) switch to `f64` + a proper
   seed and validate **statistically** (population curves, surface roughness) against
   ensembles. Recommend (a) once, then live in (b).
3. **The environment machinery** must be validated against a **non-flat** `data.rxn`
   (B4), which does not exist yet — create one.
4. **Output**: MSI atom/bond sets and population counts are deterministic given the
   state array; diff them directly.
