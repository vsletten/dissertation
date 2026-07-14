# 03 — Rust Workspace Design

How to rebuild the KMC model in Rust as a **shareable, generic KMC engine** cleanly
separated from the **specific kaolinite model** (expressed as config/data, not code)
and from **IO**. Written for co-building with Victor, who is learning Rust *through*
this port — so it is organized as a milestone ladder, each rung a small, runnable,
testable step that teaches one Rust idea and lands one piece of the model.

Design intent (from `projects/kmc-port.md`): **shareable tool, not bespoke model;
faithful-serial first (`04-parallelism-memo.md`); don't gold-plate a nostalgia
project.** So: minimal crates, boring dependencies, clear seams, no premature traits.

---

## 1. The core separation

The old code fused three things that should be pulled apart:

| Concern | Old code | Rust home |
|---|---|---|
| **Generic KMC** — a graph of sites, an event list, rejection-free select + advance time | tangled through `evtlist.cpp`, `actions.cpp`, `mckaol.cpp` | crate `kmc-engine` |
| **The kaolinite model** — states 100–599, the 28 reactions, the environment buckets, the build/terminate rules | hard-coded in `envrn.cpp`, `actions.cpp`, `lattice.cpp` | crate `kaolinite` + **data files** |
| **IO** — MSI/XYZ/dat readers & writers | `output.cpp`, `futil.cpp`, `ucell.cpp` readers | crate `kmc-io` |

The engine must know **nothing** about aluminum, silicon, or state code 401. It knows:
"there is a graph of sites, each with an opaque state; a model can, given a site and
its neighborhood, enumerate possible events with rates; when an event is chosen, the
model mutates state." That contract is one trait (§4).

"Model as config/data, not code" is a spectrum. Full data-driven (a reaction DSL the
engine interprets) is over-engineering for a nostalgia project. The realistic target:
**rates, state tables, the unit cell, and lattice dims live in data files; the reaction
*mechanisms* (the 16 hand-written transitions) live in the `kaolinite` crate as Rust,
behind the engine's trait.** The engine stays reusable for a *different* mineral by
writing a new model crate, not by editing the engine. §7 revisits how far to push this.

---

## 2. Workspace layout

A Cargo **workspace** (one repo, several crates, shared `Cargo.lock`). Proposed
location: its own repo (e.g. `vsletten/kmc`), *not* inside the dissertation repo —
the dissertation stays read-only reference. Extract the viewer to pair with it later
(that's the graph-viz project's call).

```
kmc/                          # workspace root
├── Cargo.toml                # [workspace] members = [...]
├── crates/
│   ├── kmc-engine/           # generic, model-agnostic KMC core
│   │   ├── src/lib.rs
│   │   ├── src/graph.rs      # SiteGraph: nodes + fixed-degree adjacency
│   │   ├── src/event.rs      # Event, EventList/Vec, rate accumulation
│   │   ├── src/select.rs     # rejection-free selection + dt = -ln(u)/R
│   │   ├── src/rng.rs        # Rng trait + a faithful ran2 impl for parity
│   │   └── src/model.rs      # the Model trait (the seam)
│   ├── kaolinite/            # THE specific model: states, reactions, build rules
│   │   ├── src/lib.rs
│   │   ├── src/state.rs      # SiteState newtype + class helpers (100–599)
│   │   ├── src/build.rs      # find_pairs, populate_solid, terminate_* 
│   │   ├── src/reactions.rs  # the 16 hydrolysis + ads/desorb transitions
│   │   ├── src/environment.rs# is_active + check_env (nth-neighbor buckets)
│   │   └── src/model_impl.rs # impl Model for Kaolinite
│   ├── kmc-io/               # readers (data.*) + writers (msi/xyz/dat/pgif)
│   │   ├── src/lib.rs
│   │   ├── src/cell.rs       # read data.cell → UnitCell
│   │   ├── src/inputs.rs     # read data.sim / data.rxn / data.lattice
│   │   ├── src/msi.rs        # write Cerius2 MSI (bug-compatible)
│   │   ├── src/xyz.rs, dat.rs, surf.rs
│   │   └── src/pgif.rs       # the new format (see 05-pgif-draft.md)
│   └── mckaol-cli/           # the binary: wire it all together (old mckaol.cpp)
│       └── src/main.rs
├── data/                     # sample inputs ported from model/ (read-only fixtures)
└── tests/                    # integration + golden-file tests
```

Dependency direction (no cycles):
```
mckaol-cli → kaolinite → kmc-engine
mckaol-cli → kmc-io    → kmc-engine   (io depends on engine types only)
kaolinite  → kmc-engine
```
`kmc-engine` depends on nothing project-specific (maybe `rand`-free; we bring our own
`ran2`). `kmc-io` knows engine + kaolinite state types to serialize them. Keep
`kaolinite` unaware of IO (it takes already-parsed structs).

---

## 3. Data model in Rust (replacing the C structs)

Teach the newtype + enum idioms here; they carry most of the safety win.

```rust
// kmc-engine: opaque to the engine
pub type SiteId = usize;                 // index into the site array
pub type SiteState = i32;                // engine treats state as an opaque tag...
// ...or better, make the engine generic over `S: Copy + Eq` so kaolinite can pass a newtype.

// kaolinite: the meaning
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct State(pub i32);               // newtype over the 3-digit code
impl State {
    pub fn class(self) -> u8 { (self.0 / 100) as u8 } // 1 Al, 2 Si, 3/4/5 O
    pub fn occupancy(self) -> i32 { self.0 % 100 }
    pub fn is_occupied(self) -> bool { self.occupancy() > 0 }
    pub const EDGE: State = State(9);
}
```
The site record becomes:
```rust
// kmc-engine
pub struct Site<S> {
    pub state: S,
    pub nbr: [Option<SiteId>; 6],   // -1 sentinel → Option; None = no neighbor
    // model-specific extras (pair, lostal) do NOT live here — see below
}
```
`pair` and `lostal` are **kaolinite-specific**, so they don't belong on the engine's
`Site`. Two clean options:
- (a) Engine stores `state: S` only; kaolinite keeps parallel `Vec<Pair>` /
  `Vec<LostAl>` arrays indexed by `SiteId`.
- (b) Engine is generic over a per-site `Aux` payload the model defines.
Start with **(a)** — simplest, and it forces the engine to stay generic. Revisit (b)
only if the parallel-array bookkeeping gets annoying.

The graph itself:
```rust
pub struct SiteGraph<S> { pub sites: Vec<Site<S>> }   // flat array, index = old linear index
```
Same flat layout and indexing as the C++ (`a*bCells*npos + b*npos + n`) so build-time
output matches for validation.

---

## 4. The engine↔model seam (the one trait that matters)

The whole "generic engine" claim rests on this trait. Keep it small.

```rust
// kmc-engine/src/model.rs
pub struct ProposedEvent {
    pub site: SiteId,
    pub rxn: u16,       // model-defined reaction id
    pub rate: f64,      // f64: see spec B8
}

pub trait Model {
    type State: Copy + Eq;

    /// Enumerate every possible event at `site` given the current graph.
    /// (Engine calls this for every site each step — the event-list rebuild.)
    fn events_at(&self, graph: &SiteGraph<Self::State>, site: SiteId,
                 out: &mut Vec<ProposedEvent>);

    /// Apply a chosen event: mutate `site` and affected neighbors.
    fn apply(&mut self, graph: &mut SiteGraph<Self::State>, ev: &ProposedEvent);
}
```
The engine's loop then reads (this is `mckaol.cpp`'s loop, generic):
```rust
pub fn step<M: Model>(graph: &mut SiteGraph<M::State>, model: &mut M,
                      rng: &mut impl Rng, scratch: &mut Vec<ProposedEvent>) -> Option<f64> {
    scratch.clear();
    for s in 0..graph.sites.len() { model.events_at(graph, s, scratch); }
    if scratch.is_empty() { return None; }
    let total: f64 = scratch.iter().map(|e| e.rate).sum();
    if total == 0.0 { return None; }
    let dt = -rng.uniform().ln() / total;
    let eps = rng.uniform();
    let mut acc = 0.0;
    let chosen = scratch.iter().find(|e| { acc += e.rate / total; eps <= acc })
                        .unwrap_or(scratch.last().unwrap());
    model.apply(graph, chosen);
    Some(dt)
}
```
Everything kaolinite — the states, `IsActive`, `CheckEnv`, the 16 transitions —
lives behind `events_at`/`apply`. Swap `Kaolinite` for another `Model` and the engine
is unchanged. This is the payoff of the separation and worth getting right early
(milestone 4), even though the first version of `Kaolinite` will be a stub.

Two deliberate simplifications vs. the old code: (1) the engine owns the event vector
and reuses it (no per-step `new`/`delete` linked list — a `Vec` cleared each step is
both faster and idiomatic); (2) `f64` throughout (spec B8).

---

## 5. RNG parity

To validate against the C++ (spec §C), port `ran2` faithfully behind a trait:
```rust
pub trait Rng { fn uniform(&mut self) -> f64; }   // (0,1)
pub struct Ran2 { /* idum, idum2, iv[32], iy */ }  // NR ran2, ported literally
```
This lets a one-time bitwise-ish parity run use `Ran2`, then everyday runs swap in a
modern PRNG (`rand`'s `SmallRng`) without touching the engine. Note spec B2: the C++
ignores the file seed and always effectively seeds from 0→1; the Rust port should read
the seed *correctly* but offer a `--legacy-seed` toggle to reproduce the old fixed
behavior when chasing parity.

---

## 6. Milestone ladder (co-build order)

Each milestone is a self-contained session: compiles, runs, has a test, teaches one or
two Rust concepts. Ordered so the **deterministic, RNG-free structural build comes
first** — it's the cleanest validation checkpoint (spec §C1) and the gentlest Rust
on-ramp.

**M0 — Workspace skeleton & the CLI that reads `data.sim`.**
Create the workspace, the five crates, a `main.rs` that parses `data.sim` and prints
it. *Teaches*: cargo workspaces, modules, `Result`/`?`, reading files, structs.
*Lands*: `kmc-io::inputs` for `data.sim` (read all 5 fields — fixes spec B2 from day
one). *Test*: parse the sample `data.sim`, assert fields.

**M1 — The unit cell.** Port `data.cell` reader → `UnitCell` with positions +
neighbor template. *Teaches*: `Vec`, fixed arrays `[T; 6]`, nested parsing. *Test*:
round-trip the sample cell, assert `Npos == 26` and a couple neighbor offsets.

**M2 — Build the lattice graph (deterministic).** Tile the cell into `SiteGraph`,
resolve neighbors with the exact boundary-condition logic (periodic in-plane, open
along `SurfacePlane`). *Teaches*: `Option<SiteId>` replacing `-1`, ownership of a big
`Vec`, indexing math. *Test (golden)*: for the sample inputs, the built `state[]` and
`nbr[]` arrays match a fixture dumped from the C++ (add a tiny C++ debug dump once, or
compare against `start.msi` atom positions).

**M3 — Structural setup: `find_pairs`, `populate_solid`, `terminate_surface`,
`terminate_lattice`.** The gnarliest deterministic code; port literally, one function
per session if needed. *Teaches*: mutable iteration, split borrows (you'll fight the
borrow checker here — good; §8), match on state codes. *Test (golden)*: full initial
`state[]` array matches the C++ start state, and `start.msi` atom set matches
byte-for-relevant-fields. **This is the milestone that proves faithfulness** — after
M3 the initial configuration is provably identical.

**M4 — The engine seam + a trivial model.** Define `Model`, `ProposedEvent`, the
generic `step`, and `Ran2`. Implement a **stub** `Kaolinite` whose `events_at`
proposes nothing, just to compile the loop. *Teaches*: traits, associated types,
generics, closures. *Test*: `step` on an empty model returns `None`.

**M5 — Reactions, part 1: adsorption/desorption + R0/R1.** Implement `is_active` and
`check_env` for the classes these touch, plus `apply` for the simplest transitions.
Run the loop; watch populations change. *Teaches*: bigger `match`, mutating neighbors
through the graph, the parallel `pair`/`lostal` arrays. *Test*: a hand-checked
1-cell scenario; population counts move the expected direction.

**M6 — Reactions, part 2: the full hydrolysis set R2–R15.** The `lostal`/`pair`
bookkeeping and the `ran2()<0.5` proton coin (R4, R9). *Teaches*: nothing new, lots of
care. *Test*: with `Ran2` + `--legacy-seed`, a short run reproduces the C++ trajectory
for N steps (parity check, spec §C2a).

**M7 — Output writers.** `writeData` (as **one appended file** — spec B1, with a note),
`writeMSI` (bug-compatible: shared atom/bond ID space, suppressed periodic bonds,
"9 F" hack, but derive Cartesian basis from `data.cell` — spec B5, flagged), `writeXYZ`,
`writeSurf`. *Teaches*: `Write`/`fmt`, formatting floats. *Test*: golden MSI/dat diffs
against fixtures.

**M8 — The environment actually matters.** Author a **non-flat** `data.rxn` (spec B4)
and add tests that exercise `check_env` buckets. *Teaches*: test design. *Lands*: the
first real validation of the nth-neighbor feature — the scientific point of the model.

**M9 (optional) — PGIF output.** Emit the new property-graph format (`05-pgif-draft.md`)
alongside MSI, so the viewer can migrate. *Teaches*: serde, columnar layout.

**M10+ (later, appetite-permitting) — parallelism.** See `04-parallelism-memo.md`.
Not before M8; the serial version must be trusted first.

Diffusion / BFS (spec B6) are **not** on this ladder — dead code, don't port.

---

## 7. How far to push "model as data"

Three tiers; pick per appetite, don't over-invest up front:
1. **Tier 1 (target for the port):** dims, unit cell, states-as-initial-classes, and
   *all rate tables* are data (`data.*` files, parsed by `kmc-io`). Reaction
   *mechanisms* and environment buckets are Rust in `kaolinite`. The engine is fully
   generic. → This already delivers "shareable tool": a new mineral = a new model
   crate + new data, engine untouched.
2. **Tier 2 (nice, optional):** describe the state-transition maps
   (the big `match` tables in `terminate_*`, adsorb/desorb) as data tables the
   `kaolinite` crate loads, shrinking hand-written match arms. Reduces the port's most
   error-prone code to a checked table.
3. **Tier 3 (over-engineering — resist):** a general reaction DSL the engine
   interprets. Not worth it for one mineral and a nostalgia budget.

Recommend Tier 1 firmly, Tier 2 only if the transition matches prove painful to
maintain, never Tier 3.

---

## 8. Rust-learning notes specific to this port

- **The borrow checker will bite in M3/M5**: mutating `graph.sites[nbr]` while looking
  at `graph.sites[site]`. Idiomatic fixes: read what you need into locals first, then
  mutate; or use `split_at_mut`/indices carefully; or a small helper that takes
  `&mut Vec` + two indices. This friction is *teaching* the aliasing rules the C++
  ignored — lean into it rather than reaching for `RefCell`.
- **`Option<SiteId>` everywhere** replaces the `-1` sentinel; `if let Some(nbr) = ...`
  replaces `if (nbr >= 0)`. This makes the many "ran into edge" guards explicit.
- **Newtypes** (`State`, maybe `SiteId`) prevent the class of bug where an int meaning
  "state" gets used as an index.
- **`enum` for site class** (`Al, Si, SiOSi, SiOAl2, AlOHAl, Edge`) can replace
  `state/100` dispatch with an exhaustive `match` the compiler checks — a real upgrade
  over the C++ `switch` with silent `default`.
- **Keep `unsafe` out.** Nothing here needs it.
- **Tests as the safety net for a faithful port**: golden fixtures at M2/M3/M7 are how
  you *know* the port matches without re-reading C++ each session. Write them as you go.

---

## 9. Dependencies (keep it boring)

- `kmc-engine`: std only (bring our own `ran2`); optional `rand` behind a feature for
  the non-parity RNG.
- `kmc-io`: std for the legacy text formats; `serde` + `serde_json` (or `rmp-serde`
  for MessagePack) only when PGIF lands (M9).
- No async, no web framework, no proc-macro-heavy deps. This is a batch simulator.

Right-sized: four small crates + a CLI, boring deps, a milestone ladder that keeps
every session runnable and every faithful piece tested against the old output.
