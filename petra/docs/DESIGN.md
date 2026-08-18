# Petra — a general lattice KMC engine for mineral elementary reactions

**Status: v0 design + compiling skeleton.** Greenfield by decision (2026-08-12):
not an extension of `kmc-rs`, though it borrows every lesson the port surfaced.
The name is a placeholder (πέτρα, "rock"); rename freely.

## 1. What this is, and where it comes from

The dissertation model (`legacy/cpp-model`, ported in `kmc-rs`) is a kinetic
Monte Carlo simulation of one mineral — kaolinite — with one process family —
dissolution/growth — where the crystal structure, the site-state vocabulary
(the 100–599 codes), the 28 reactions, and the environment logic that selects
among rate variants are all **hard-coded**. Adding a reaction means editing
five source files; modeling a different mineral means writing a new program.

Petra inverts that: the program is a fixed, chemistry-free engine, and *the
mineral is data*. A user supplies one input deck describing:

1. **Structure** — an arbitrary crystalline unit cell at atomic resolution:
   cell parameters, sites, occupying species, and bond topology (including
   bonds that cross cell boundaries).
2. **State space** — the local states each kind of site can be in
   (protonation, hydration, hydrolysis progress, occupancy…), as *named
   declarations*, not integer codes baked into code.
3. **Reactions** — a set of elementary steps written as local graph-rewrite
   rules, each with transition-state-theory kinetic parameters (activation
   energy, prefactor, reaction energy) and optional environment-dependent
   modifiers.

The engine compiles the deck into dense runtime tables and runs rejection-free
KMC on it. Dissolution is one deck. Growth, ion exchange, surface diffusion,
Ostwald ripening of a surface phase, isotope/tracer exchange — other decks,
same binary.

### Lessons imported from mckaol / kmc-rs

| Lesson | Consequence in Petra |
| --- | --- |
| The engine/model seam works (`kmc-engine::Model` proved it) | Same seam, but the "model" side is *interpreted from the deck*, not hand-written per mineral |
| Hard-coded state integers (100–599) made every change invasive | States are deck-declared names, compiled to dense ids at load |
| Environment rate variants (Check100–500, up to 40 variants) are the real complexity | A declarative guard/modifier language designed to express exactly that class of logic |
| f32 accumulation order decided trajectories (spec §C2a) | f64 everywhere; no ordering-sensitive summation contracts |
| The seed swallow (spec B2) cost the model its ensembles | Seedable RNG, seed in the deck, CLI override, ensemble-first validation |
| O(N) event-list rebuild every step (`CreateEventList`) | Incremental updates: only sites inside the fired event's influence radius are re-enumerated; O(log N) selection |
| The `nbr[6]` fixed array + out-of-bounds phantom (reform R1) | Variable-coordination CSR adjacency; bounds are structural |
| Golden-gate testing was the port's backbone | Analytic-solution gates + detailed-balance linting + a kaolinite deck cross-checked against `kmc-rs` corrected mode |

## 2. Scope

**In scope (v1).** On-lattice KMC on arbitrary periodic crystals: any
triclinic cell, any number of site kinds and species, arbitrary bond
topology and coordination, slab or bulk geometry, per-axis periodic/open
boundaries. Elementary events are state transitions on that fixed site graph
with TST rates. This covers dissolution, growth, adsorption/desorption,
hydrolysis chains, surface diffusion (as paired vacate/occupy rewrites), and
exchange reactions.

**Designed-for but deferred.** Defects and substitutions: the data model
keeps *occupant species* separate from *site kind* precisely so a Fe³⁺-for-Al
substitution, a vacancy population, or a doped deck is initialization data
plus (optionally) extra reactions — no schema change. See §6.

**Out of scope (v1), acknowledged.** Off-lattice/relaxed positions, lattice
strain, long-range electrostatics, explicit solvent, non-axis-aligned Miller
surfaces (workaround: supply a re-oriented cell), and on-the-fly barrier
calculation (barriers are inputs; Petra is not a DFT driver).

## 3. Conceptual model

### 3.1 Crystal template → lattice instance

A **unit cell** declares:

- cell parameters `a b c α β γ` (the fractional→Cartesian matrix is derived);
- **sites**: index, *site kind* (a name like `Al_oct`, `Si_tet`, `O_bridge`),
  fractional coordinates;
- **bonds**: `(site_i, site_j, Δcell)` where `Δcell = (da, db, dc)` says which
  periodic image of `j` the bond reaches. Bonds may carry an optional *label*
  (e.g. `pair` — the second bridging oxygen linking the same two cations,
  generalizing the legacy `pair` field into ordinary labeled topology).

A **lattice** is the cell tiled `na × nb × nc` with a per-axis boundary
condition:

- `periodic` — wrap (the legacy in-plane behavior);
- `open` — bonds crossing the face are dropped (the legacy free surface);
- `fixed` — like `open`, but sites in the terminal layer are frozen: they
  participate in neighbors' environments but never react (a clean way to
  anchor a slab, and the anchor set for cluster removal, §5.5).

Adjacency is stored CSR-style (offset table + flat neighbor array), so a site
has however many neighbors its chemistry gives it — 2, 4, 6, 12 — with no
fixed-width array and no sentinel `-1`s. Distance-2 neighborhoods (needed by
kaolinite's surface-reachability guards) are iterated on demand through the
CSR rows, not stored.

### 3.2 Site state = kind + occupant + configuration

The dynamic state of a site is the pair **(occupant, configuration)**:

- **occupant**: a species from the deck's species list, or `vacant`;
- **configuration**: a deck-declared name for the local arrangement —
  ligands, protonation, hydrolysis progress (`"3OH"`, `"bridging"`,
  `"hydrolyzed_Si_side"`, …).

Each site kind declares its allowed states. The legacy encoding maps over
directly — `103` becomes `Al_oct: occupant=Al, config="2OH"`; `301`/`302`
become `O_bridge_SiSi: occupant=O, config="intact"/"hydrolyzed"`; `100`
becomes `Al_oct: vacant`. At load time every declared `(kind, occupant,
config)` triple is interned to a dense `StateId(u16)`, and every state *set*
mentioned by a rule becomes a bitset — the hot loop never touches strings.

Keeping occupant separate from configuration (rather than one flat state
name) is the defect/substitution design-ahead: rules may match "any trivalent
occupant of an octahedral site," and a substituted deck reuses every rule
that doesn't distinguish occupants.

### 3.3 Reactions as guarded local rewrites

An elementary reaction is:

```
center selector  +  guards over the neighborhood  →  effects  @ rate(T, environment)
```

- **Center selector**: site kind + a set of states the center may be in.
- **Guards**: conjunctions of *neighbor count predicates*: "at least/exactly/
  at most `n` neighbors (optionally: of kind K, via bonds labeled L, at graph
  distance 1 or 2) whose state is in set S." This is the declarative form of
  `IsActive`/`CheckEnv`: e.g. the corrected R1 surface rule is "≥1 distance-2
  neighbor in the hydrolyzed set."
- **Effects**: state assignments to the center and to matched neighbors
  (selected by the same selector machinery), including occupant changes
  (adsorption writes an occupant into a vacant site; desorption removes it).
  An effect may have **weighted random branches** — generalizing the legacy
  R4/R9 "which side keeps the proton" coin flip into deck syntax.
- **Rate**: §4.

Every reaction gets a compiled **influence radius**: the maximum graph
distance at which its effects can change some other site's guard evaluation
(effect reach + max guard distance). The engine uses it for incremental
event-list maintenance (§5.2).

### 3.4 What a deck looks like

The full schema lives in `petra-deck`; an annotated Kossel-crystal
dissolution deck ships in `examples/kossel.toml` and doubles as the tutorial.
The shape:

```toml
[deck]                    # name, comment, schema version
[cell]                    # a b c alpha beta gamma
[[species]]               # name, optional mass/charge (bookkeeping only in v1)
[[kinds]]                 # site-kind names + allowed states
[[cell.sites]]            # kind, fractional xyz
[[cell.bonds]]            # i, j, dcell, optional label
[lattice]                 # dims, per-axis boundary, initial fill rules
[thermo]                  # temperature, per-species chemical potentials/activities
[[reactions]]             # the rewrite rules (see below)
[simulation]              # steps/t_max, seed, output cadence
```

One reaction, concretely (Kossel dissolution with a bond-counting barrier):

```toml
[[reactions]]
name   = "detach"
center = { kind = "M_site", state = ["occupied"] }
rate   = { arrhenius = { prefactor = 1.0e13, ea = 10.0 } }   # kcal/mol
# each occupied neighbor raises the barrier — BEP-style bond counting:
[[reactions.modifiers]]
per_match = { dea = 5.0 }
select    = { distance = 1, state = ["occupied"] }
[[reactions.effects]]
target = "center"
set    = { occupant = "vacant", config = "empty" }
```

And the legacy Si–O–Si hydrolysis pair, as a fragment of a future kaolinite
deck (states renamed from 301/302):

```toml
[[reactions]]
name    = "SiOSi_hydrolysis"
center  = { kind = "O_SiSi", state = ["intact"] }
guards  = [ { distance = 2, min = 1, state = ["@hydrolyzed_any"] } ]  # surface-reachable
rate    = { arrhenius = { prefactor = 1.0, ea = 0.0 } }
reverse = { name = "SiOSi_condensation", de = 2.6 }   # k⁻ = k⁺·exp(ΔE/RT), detailed balance
[[reactions.effects]]
target = "center"
set    = { config = "hydrolyzed" }
```

(`@hydrolyzed_any` is a deck-defined state *alias* — named unions keep guard
lists readable where the legacy code had literal lists like
`303|404|405|406|408|409`.)

## 4. Rate physics

All rates are computed at load/update time in f64. The runtime is canonical
**kcal/mol and Kelvin** (the dissertation's convention;
`R = 1.987×10⁻³ kcal·mol⁻¹·K⁻¹`); each deck declares its own energy unit —
`units = "kcal/mol" | "kJ/mol" | "eV"` in `[deck]` — and every energy-valued
field (Ea, ΔEa modifiers, Δμ, ΔH‡, ΔS‡ per K) is converted once at compile
time. Temperatures are always Kelvin, prefactors always s⁻¹. (Decided
2026-08-12: selectable from day 1.)

**Base rate expressions** (per reaction, one of):

- `constant`: `k = k₀` — for lumped/experimental rate constants;
- `arrhenius`: `k = A·exp(−Ea/RT)`;
- `eyring`: `k = (k_B T/h)·exp(ΔS‡/R)·exp(−ΔH‡/RT)` — full TST when the
  activation entropy/enthalpy are known (e.g. from DFT + vibrational
  analysis).

**Solution coupling.** A reaction may declare `consumes`/`produces` species
from solution. Attachment rates get the standard supersaturation factor:
`k_att = k_base · a_X · exp(Δμ_X/RT)` with activity and chemical potential
from `[thermo]` — the generalization of the legacy adsorption law
`a·exp(dE/RT)·exp(dm/RT)`. v1 treats the solution as an infinite reservoir
(fixed activities); a finite, depleting reservoir is roadmap.

**Environment modifiers.** After the base rate, an ordered list of modifiers
adjusts the barrier per local environment:

- `by_count { dea = [e₀, e₁, e₂, …] }` with a neighbor selector — a
  **tabulated** ΔEa indexed by match count, last entry extending to all
  higher counts. This is the general form: barriers are rarely linear in
  coordination (Victor, 2026-08-12: "I would bet that it is usually
  nonlinear"), so a measured or computed per-coordination table goes in
  verbatim rather than being forced through a slope;
- `per_match { dea }` — the linear convenience case (`by_count` with a
  constant increment): bond-counting / BEP-style, used by the Kossel
  example;
- `when { guard } → { dea | factor }` — a discrete override when a guard
  holds (for cases keyed on something other than a count, e.g. the legacy
  Check400 variants that key on the *pair* partner's exact state via the
  `pair` bond label).

This language was chosen by test: it can express every variant family in
`envrn.cpp` (which are all either counts of neighbors in a state set,
membership tests on a specific labeled neighbor, or small combinations),
while staying analyzable — the compiler can enumerate a reaction's complete
variant table for review, which is exactly the table `data.rxn` used to hold
as opaque numbers.

**Reversibility and detailed balance.** A reaction may declare a `reverse`
with reaction energy `ΔE`; the compiler generates the reverse rule
(effects inverted) with `k⁻ = k⁺·exp(ΔE/RT)` — the legacy convention — or
accepts an explicitly parameterized reverse and then **lints**: for every
reversible pair, and for every closed thermodynamic cycle it can find in the
state graph of each site kind, it checks that reaction energies sum
consistently (ΣΔE ≈ 0 around a cycle). Violations are warnings with the
cycle printed. This makes the classic silent KMC sin — a deck that pumps
free energy out of nothing — loud.

## 5. Engine

### 5.1 Selection: rejection-free, two-level

Standard n-fold-way/direct-method KMC, structured for scale:

- Per site: a small cached list of currently-possible events (reaction id +
  resolved rate), rebuilt only when the site is dirtied.
- Global: a **Fenwick (binary-indexed) tree** over per-site total rates.
  Selection draws `u·R_total`, descends the tree to a site in O(log N), then
  scans that site's short event list. Time advances `Δt = −ln(u′)/R_total`.
- f64 throughout. The tree maintains sums incrementally; to bound float
  drift from many incremental updates, the tree rebuilds from leaves on a
  cadence (and the rebuild is also the cheap correctness oracle in tests).

### 5.2 Incremental updates via influence radius

When an event fires, its effects change states at known sites. The dirty set
is those sites plus every site within `influence_radius(rxn)` of them
(precomputed per reaction, §3.3). Only dirty sites re-enumerate events and
push new totals into the tree. With local chemistry this is O(z·r) per step
— independent of lattice size — versus the legacy O(N·z) full rebuild. A
`--paranoid` mode re-enumerates everything every step and asserts agreement;
that mode *is* the migration-grade test that the locality analysis is right.

### 5.3 RNG and reproducibility

PCG64, seed from the deck (CLI-overridable), one explicitly-threaded stream;
branch draws (§3.3) come from the same stream in a documented order. Runs are
reproducible given (deck, seed, binary version); ensembles are seed sweeps.
No global RNG, no hidden draws — the seed-swallow lesson, inverted.

### 5.4 Termination and stiffness

Runs end on step count, simulated-time horizon, or a deck-declared
predicate (e.g. "fewer than n sites of state S remain"). Known limitation,
stated up front: fast reversible pairs (stiffness) waste steps flip-flopping.
v1 ships detection (report the top rate-sum consumers); mitigation
(superbasin/τ-leaping methods) is explicitly roadmap, not silently absent.

### 5.5 Cluster removal as policy, not hard-code

The legacy BFS removal of detached clusters becomes an optional post-event
policy: `detach_clusters = { anchor = <site selector>, action = "remove" }`,
run only when an event actually severed a bond (cheap incremental check
first, BFS only on suspicion). Off by default — detachment is a modeling
*choice* (those clusters are colloids, arguably), and Petra makes it
visible. (Confirmed 2026-08-12: default doesn't matter as long as it's a
one-line deck toggle — which it is.)

## 6. Defects and substitutions (the design-ahead)

Requested for "down the road"; the commitments that make it cheap later:

1. **Occupant ≠ kind** (§3.2). A substitutional defect is a site whose
   occupant differs from the pristine deck — same kind, same topology, so
   every occupant-agnostic rule still applies, and occupant-specific rules
   (an Fe–O–Si bridge hydrolyzing faster) are just additional reactions with
   a narrower selector.
2. **Initial-condition rules** (implemented): any `[[init]]` pass takes
   `probability = p` (random substitution by fraction — `Fe3+ replaces Al
   on Al_oct at 0.02` is `probability = 0.02` + `set` to an Fe-occupant
   state) and/or `sites = [[a, b, c, t], …]` (explicit site lists; seeded
   vacancy populations are the same with a vacant-state `set`). Both
   compose with the existing region/state/guard filters. Deterministic
   under the run seed: one init RNG stream, salted apart from the dynamics
   stream, drawn in pass order then site-index order.
3. **Vacancies are ordinary states** (`occupant = vacant`), so vacancy
   diffusion is a rewrite rule pair, not a mechanism.
4. **Charge is declared** on species now (inert bookkeeping), so a future
   charge-balance linter or coupled-substitution constraint (Al³⁺ + Na⁺ ↔
   Si⁴⁺-style) has the data it needs without a schema break.

Non-stoichiometric *structural* disorder (stacking faults, dislocations) is
harder — it breaks the "tiled unit cell" premise — and is honestly deferred;
the CSR graph doesn't care, but deck *authoring* of such structures needs a
builder tool (roadmap: accept a pre-built graph JSON as an alternative to
cell tiling).

## 7. Architecture

```
petra/
  crates/
    petra-core/   chemistry-free runtime: lattice graph, interned states,
                  rate eval, guards, compiled reactions, KMC engine
    petra-deck/   the TOML schema (serde), validation, and the compiler
                  deck → petra-core runtime tables
    petra-cli/    binary: load deck → build lattice → run → write outputs
  examples/       decks (kossel.toml first; kaolinite deck is milestone P4)
  docs/           this file
```

The core is engine-only and knows no chemistry names; the deck crate owns
all vocabulary and produces only dense ids — the same seam `kmc-engine`
proved, moved one level up (model = interpreter + tables instead of
hand-written trait impl). Everything is a library; the CLI is thin, so a
future Python binding (PyO3) or a `graph-viz` live socket wraps the same
core.

**Outputs:** population time series per (kind, state) as CSV; snapshots as
extended-XYZ (positions from the cell matrix); an append-only JSONL event log
(step, time, site, reaction) that is both the replay record and the
`graph-viz` feed — PGIF export from a snapshot + event log is a planned
`petra-io` module, wired to this repo's existing viewer rather than a new one.

## 8. Validation plan

Petra cannot be validated by bitwise parity (different algorithm, different
float width — deliberately). The replacement ladder:

1. **Unit gates** (in-tree now): Fenwick tree vs. linear scan on random
   workloads; rate formulas vs. hand-computed values; lattice tiling vs.
   hand-enumerated neighbor tables for a known cell (periodic and open).
2. **Analytic gates**: two-state site ensemble → occupancy ratio
   `exp(−ΔE/RT)` within statistical error; Langmuir adsorption isotherm from
   an attachment/detachment deck; mean event waiting times vs. `1/R_total`.
3. **Detailed-balance linter** (§4) as a load-time gate on all shipped decks.
4. **The kaolinite cross-check** (milestone P4): author the kaolinite deck
   from the state map in `kmc-rs`, run ensembles, and compare *observables*
   (population trajectories, surface roughness, dissolution rate vs. Δμ)
   against `kmc-rs` **corrected mode** (post-REFORM_PLAN — legacy mode's
   phantom-activation wart R1 is precisely what Petra's guards won't
   reproduce, and shouldn't). Agreement bands, not bitwise.
5. **`--paranoid` differential mode** (§5.2) on every example deck in CI.

## 9. Milestones

| # | Lands | Exit test |
| --- | --- | --- |
| P0 | ✅ This doc + workspace skeleton: core types, Fenwick engine, deck schema, Kossel example compiles & runs | `cargo test` green; smoke test runs 10⁴ steps |
| P1 | Mostly landed with P4: distance-2 guards, labeled bonds, branches, frozen filters, set/shift/map effect ops. Remaining: variant-table dump, auto-reverse generation | variant-table dump matches hand analysis |
| P2 | ✅ Analytic gates (two-state, Langmuir, waiting times); ensemble runner. Remaining: detailed-balance linter | gates green in CI |
| P3 | ✅ `--viz` writes a PGIF snapshot + JSONL event log (`petra-io`); graph-viz grew a trajectory playback mode (play/scrub both directions, speed dial, in-place instance updates) with the format spec'd in graph-viz/docs/PGIF.md §6a; kaolinite trajectory bundled as a viewer demo, verified end-to-end in a real browser. Plain-XYZ via `--xyz`. Stage 2: `petra-wasm` (wasm-bindgen) runs the engine in-browser; graph-viz live mode generates trajectory events on demand (play/scrub = simulate), deck editable in the page. Remaining: binary event log for very long runs | kaolinite run plays in graph-viz ✅ |
| P4 | ✅ The kaolinite deck (`examples/kaolinite.toml`) with the legacy `data.rxn` numbers, **corrected physics** (bounded surface test, R13 restored, Δμ un-swapped; as-implemented termination maps). Structural golden gate: Petra's built lattice matches the kmc-rs (bitwise-parity-tested) builder **site-by-site** — states, EDGE⟺frozen, degrees — plus a 5k-step dynamics gate with paranoid differential checks. The build machinery this required is general: deck-declared init passes (region clears, per-neighbor termination maps, hydroxyl counting), the pr1/pr2 state split replacing the hidden `lostal` field, and labeled-bond addressing | `kaolinite_golden.rs` green |
| P5 | Partial: **defect strain fields landed** — `[[defects]]` analytic screw/edge elastic fields, per-reaction `strain` coupling with detailed-balance conventions, etch-pit gates and demo; full treatment (topological screw rewiring, self-consistent relaxation) designed in **docs/STRAIN.md**. **Substitution init rules landed** — `[[init]]` `probability`/`sites`, seeded + deterministic, gated in `init_substitution.rs` (§6). Remaining: Fe-substituted kaolinite deck | etch-pit gate green ✅ |
| P6 | Performance: profiling, superbasin detection report, big-lattice runs | 10⁶-site lattice at interactive step rates |

**P4 cross-check caveat, stated honestly:** trajectory-level agreement bands
against kmc-rs are *not possible yet* — kmc-rs today implements only legacy
mode, whose forward-hydrolysis surface test carries the out-of-bounds
phantom (reform R1) that Petra deliberately does not reproduce (~4× the
forward-hydrolysis event availability). The site-by-site structural gate is
the strongest cross-check available until kmc-rs's corrected mode lands;
at that point ensemble observables (population trajectories, dissolution
rate vs Δμ) become comparable.

## 10. Decisions (resolved 2026-08-12)

The v0 open questions, answered by Victor:

1. **Units: selectable per deck from day 1.** Implemented — `[deck] units`
   accepts `kcal/mol` (default), `kJ/mol`, `eV`; conversion happens once at
   compile time, runtime stays canonical kcal/mol (§4).
2. **Kaolinite deck numbers: the old model's numbers are all there is.**
   Context, since the original phrasing confused: the legacy `data.rxn`
   carries an environment-variant table per reaction (15–40 variants), but
   every variant within a reaction holds *identical* numbers — e.g. all 15
   Si-O-Si hydrolysis variants read `(1, 2.6)` — so the environment
   machinery never actually differentiates rates in the shipped inputs
   ("flat tables", kmc-rs spec B4). There are no other numbers. Resolution:
   the P4 kaolinite deck uses the legacy values as-is — which is *better*
   for the cross-check anyway, since both codes then run identical physics;
   flat tables also collapse cleanly to plain base rates in deck form.
   Differentiated per-environment barriers become a science question on top
   of a validated deck (`by_count` tables are the natural slot for them).
3. **Modifiers: nonlinear in general.** ΔEa is usually *not* additive per
   matching neighbor. Implemented — `by_count` tabulated ΔEa is the general
   primitive; `per_match` remains as the linear convenience (§4).
4. **Cluster removal: default is a don't-care; changeability is the
   requirement.** Stays off by default as a one-line deck toggle (§5.5).
