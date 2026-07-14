# 01 — Codebase Map

Archaeology of the C++ kinetic Monte Carlo (KMC) model at
`/mnt/data/vsletten/src/vsletten/dissertation/main/model`. Read-only reference;
nothing in the dissertation repo was modified. This document maps *what is there*:
files, control flow, data structures, and build. Physics and algorithm detail go
in `02-model-spec.md`; warts are cataloged there too (this map only points at them).

## 1. What the program is

`mckaol` ("Monte Carlo kaolinite") is a single-binary, single-threaded C++11
command-line program. It simulates dissolution/growth of a kaolinite (Al-Si
clay mineral) surface as a lattice of reaction sites, using a rejection-free
kinetic Monte Carlo loop (Gillespie / n-fold-way flavor — see spec). It reads
four fixed-name input files from the current working directory, runs a fixed
number of MC steps, and writes snapshot/geometry/population files, including
**Cerius2 MSI** structure files that the `viewer/` React app consumes.

There are **no command-line arguments** and **no config path**: input file names
(`data.sim`, `data.rxn`, `data.cell`, `data.lattice`) and all output file names
are hard-coded. You configure a run by editing files in the run directory.

Total source: ~2,700 lines across 13 `.cpp` + 13 `.hpp` translation units.

## 2. Directory contents

```
model/
├── Makefile                 # g++ -std=c++11 -O3 -ffast-math, static list of .cpp
├── README.md                # user-facing docs (partly aspirational — see notes)
├── *.cpp / *.hpp            # the 13 translation units (below)
├── data.sim/.rxn/.cell      # active input files (fixed names)
├── data.lattice[.NxM]       # active lattice dims + saved variants (50x10, 100x10, 500x10)
├── input/prot.*             # "prototype"/template copies of the four inputs
├── step*.dat / *.msi        # committed sample OUTPUTS from past runs (not inputs)
├── results.dat, 8000.out,   # more stale sample outputs
│   surfSi.out, surfAl.out
├── *.o, mckaol              # committed build artifacts (should be gitignored)
├── RCS/                     # 1990s RCS history (*.c,v / *.h,v) — the C ancestor,
│                            #   incl. reactions.c/reactionlist.h (renamed since)
└── .vscode/, .cursorignore  # editor cruft
```

The `RCS/` directory is a genuine time capsule: it holds the original **C** (not
C++) sources under RCS version control (`reactions.c,v`, `reactionlist.h,v`,
`mckaol.c,v`, …). The current tree is a C++ re-skin of that C code — which explains
much of the style (raw `malloc`, `int` booleans, `-1`-as-null-index, class methods
that read like free functions with a `lattice->` prefix). Useful as provenance if a
behavior in the current code looks accidental; the RCS originals can disambiguate
"always did this" from "broke during the C→C++ port". Not required for the Rust port.

## 3. Translation units (by role)

Grouped by concern. Line counts approximate; "→" = "depends on / calls".

### Entry point / orchestration
| File | Lines | Role |
|---|---|---|
| `mckaol.cpp` | 115 | `main()`. Wires everything: build sim/reactions/cell/lattice, terminate surface, then the MC loop, then final output. |

### Domain data structures
| File | Lines | Role |
|---|---|---|
| `ucell.hpp/.cpp` | 43/86 | `UnitCell`, `CellSite`, `NeighborSite`, `CellDimensions`. Reads `data.cell`: fractional atom positions + per-position neighbor template (6 slots, with unit-cell offsets). The *motif* that tiles into the lattice. |
| `lattice.hpp/.cpp` | 80/381 | `LatticeSite` (the per-site record) + `Lattice` (flat array of sites). Builds the lattice by tiling the unit cell, resolves neighbor indices with boundary conditions, and does all the one-time structural setup: `FindPairs`, `PopulateSolid`, `TerminateSurface`, `TerminateLattice`. Also `RemoveUnattachedClusters`. |
| `sim.hpp/.cpp` | 18/45 | `Simulation`: run parameters from `data.sim` (nsteps, wsteps, msteps, drawbonds, ranseed). |

### Reaction model
| File | Lines | Role |
|---|---|---|
| `rxnlist.hpp/.cpp` | 41/123 | `Reaction` + `ReactionList`. Reads `data.rxn`: temperature, chemical potentials, and the 28 reaction rate tables. Converts (pre-exponential, ΔE) pairs into concrete float rates (Arrhenius / equilibrium math). |
| `envrn.hpp/.cpp` | 25/318 | `Environment`. Two jobs: `IsActive(site,rxn)` (is this reaction allowed here at all?) and `CheckEnv(site)` (which nth-neighbor environment bucket → which rate variant?). This is where the "nth-nearest-neighbor chemical state" physics lives. |
| `evtlist.hpp/.cpp` | 30/71 | `Event` + `EventList` (singly-linked list). `CreateEventList` scans every site, and for each allowed (site, reaction) pair appends an event carrying its specific rate. Rebuilt from scratch every step. |
| `actions.hpp/.cpp` | 87/480 | `Actions`. `DoEvent(list)`: pick one event weighted by rate, execute it, return Δt. `DoReaction` dispatches to 16 hand-written `DoReactionN` state-transition routines plus adsorb/desorb/diffuse. This is the largest and most intricate file. |

### Support / IO / utilities
| File | Lines | Role |
|---|---|---|
| `output.hpp/.cpp` | 18/256 | `writeMSI` (Cerius2), `writeXYZ`, `writeData` (population counts), `writeSurf` (2D surface projections), `initDatafile`. |
| `bfsearch.hpp/.cpp` | 25/81 | `BreadthFirstSearch::ColorNodes`: BFS from a seed site to mark the connected occupied cluster; used by `RemoveUnattachedClusters` to delete floaters. **Reachable only via diffusion, which never fires — effectively dead (see spec).** |
| `futil.hpp/.cpp` | 18/64 | File open/close helpers + `EatComment` (skip `#`-to-EOL comments in inputs). |
| `ran2.hpp/.cpp` | 14/77 | Numerical Recipes `ran2` PRNG (two combined LCGs + shuffle). Global `static long *idum`. |
| `myerr.hpp/.cpp` | 11/13 | `Myerr::die(msg)` → stderr + `exit(1)`. |
| `common.hpp` | 175 | No code — the canonical **state-encoding documentation** (site classes 100–599), the reaction list R0–R15, and the intended environment-indexing scheme, all as a comment block. Plus the `Color` enum for BFS. The single most useful file for understanding the model. |

## 4. Control flow

### 4.1 Startup (in `main`, top to bottom)
1. `Simulation::CreateSimulation()` — read `data.sim`.
2. `initran2(&sc->ranseed)` — hand the PRNG its seed slot. *(See spec: the file's
   seed value is not actually what ends up here — a read-ordering wart.)*
3. `output::initDatafile()` — `remove("results.dat")`. *(Nothing ever writes
   `results.dat` afterward — see spec.)*
4. `ReactionList::CreateReactionList()` — read `data.rxn`, build 28 rate tables.
5. `UnitCell::CreateUnitCell()` — read `data.cell`, build the motif.
6. `Lattice::CreateLattice(uc)` — read `data.lattice`, tile the cell into a flat
   `LatticeSite[]`, resolve each site's 6 neighbor indices (with periodic / open BCs).
7. Structural setup, in this exact order (order matters):
   - `FindPairs()` — link each bridging-O (400/500) to its "double-bridge" partner.
   - `PopulateSolid(dmSi, dmAl)` — fill a fraction (0.3/0.5/0.7 of the slab depth,
     chosen by chemical-potential sign) by promoting empty→first-occupied states.
   - `TerminateSurface()` — cap dangling bonds at the exposed surface with OH/H2O,
     adjusting cation and oxygen states accordingly.
   - `TerminateLattice()` — mark the top/bottom rows as `EDGE` (frozen boundary),
     and revert now-invalid bridges.
8. Construct `Actions(lattice)` and `Environment(lattice)`.
9. `output::writeMSI(lattice, "start", drawbonds)` — initial structure.

### 4.2 The MC loop (`for i in 0..nsteps`)
```
eventList = EventList::CreateEventList(lattice, environment, reactions)
    → for each site s:
        pick the reaction-index window [lo,hi) from s's state class
        for each reaction r in [lo,hi):
            if state == reactions[r].reactant AND environment.IsActive(s,r):
                env = environment.CheckEnv(s)          # nth-neighbor bucket
                append Event{site:s, rxn:r, rate:reactions[r].rate[env]}
if eventList == null: break                            # no moves possible → stop

dt = actions.DoEvent(eventList)
    → ratesum = Σ rate over list
    → dt = -ln(ran2()) / ratesum                       # KMC time step
    → pick event: eps=ran2(); walk list accumulating rate/ratesum until eps ≤ partsum
    → DoReaction(site, rxn): mutate this site + affected neighbors' states
if dt < 0: break                                       # reaction failed
time += dt

if wsteps and (i % wsteps == 0): writeData("step{i}.dat", …, time)   # one-row file
if msteps and i and (i % msteps == 0): writeMSI("step{i}", …)         # snapshot

DisposeEventList(eventList)                             # free the whole list
```
Key shape: **the event list is fully regenerated every step** — no incremental /
locality-aware update. That O(N_sites) rebuild per step is the dominant cost and the
central thing the Rust port's trait seams should isolate (see `04-parallelism-memo.md`).

### 4.3 Shutdown
`writeData("end.dat")`, `writeSurf`, `writeXYZ`, `writeMSI("end")`, then explicit
`delete`/`Dispose` of every object.

## 5. Core data structures

### `LatticeSite` (the hot record — `lattice.hpp`)
```cpp
int a, b;      // unit-cell coordinates in the 2D tiling
int n;         // which position within the unit cell (0..Npos-1)
int state;     // the 3-digit state code (100..599, or EDGE=9) — the whole model
int pair;      // for 400/500 bridging O: index of the partner O in the double bridge (-1 if none)
int lostal;    // for 400s: which Al was lost when the bridge half-broke (-1 if none)
int nbr[6];    // resolved neighbor site indices (-1 = no neighbor / open boundary)
Color color;   // BFS marker (UNREACHABLE/ENQUEUED/DISCOVERED)
```
The lattice is `LatticeSite sites[Num_Sites]`, a flat array; index
`i = a*(Num_bCells*Npos) + b*Npos + n`. Everything is integer-index-based; there are
no pointers between sites, only `nbr[]` indices. This is already a graph in
adjacency-list form (fixed degree ≤ 6). `state` is the single mutable field during
simulation; `nbr`, `pair`, `a/b/n` are set once at build time. `lostal`/`pair` are
mutated by a few reactions.

**State encoding** (full table in `common.hpp`; summarized in spec): first digit =
site class (1 Al, 2 Si, 3 Si-O-Si oxygen, 4 Si-O-Al₂ oxygen, 5 Al-OH-Al oxygen); the
last two digits = occupancy/coordination/protonation. `state % 100 == 0` means empty;
`EDGE = 9` is a sentinel for frozen boundary sites; `WRONG = 99` marks a
wrong-cation adsorption (never actually produced — see spec).

### `UnitCell` / `CellSite` (`ucell.hpp`)
The motif read from `data.cell`: `Npos` positions, each with fractional `x,y,z`, a
`state` (initial class), and a 6-entry neighbor template `nbr[j] = {n, a, b, c}`
(target position + integer unit-cell offset). `GetNeighbor` in the lattice consumes
these offsets to resolve concrete neighbor indices per tiled site. `c` (third-axis
offset) is read but the surface is a single 2D sheet, so `c` is effectively unused
in neighbor resolution (only `a`,`b` participate).

### `Reaction` / `ReactionList` (`rxnlist.hpp`)
```cpp
int reactant;   // the state code this reaction consumes
int info;       // product state (used only for adsorption bookkeeping)
int nrates;     // length of rate[]
float *rate;    // one rate per nth-neighbor environment bucket
```
28 reactions: 16 hydrolysis (8 forward/reverse pairs, indices 0–15), 4 adsorption
(16–19), 4 desorption (20–23), 4 diffusion (24–27, copies of desorption, never used).
`rate[]` is indexed by the environment bucket that `CheckEnv` returns.

### `Event` / `EventList` (`evtlist.hpp`)
`Event{site, rxn, rate}`; `EventList : public Event { EventList* next; }`. A raw
singly-linked list, heap-allocated node per possible move, rebuilt each step. `rxn`
for diffusion would carry `OFFSET(100)+target` but that path is dead.

## 6. Input files (formats live in the spec; here just the map)

All four are whitespace-delimited text with `#`-to-EOL comments, read positionally
via `Futil::EatComment` between fields. Fixed names in CWD.

| File | Consumed by | Contents (one-liner) |
|---|---|---|
| `data.sim` | `sim.cpp` | nsteps, wsteps, msteps, seed, drawbonds |
| `data.rxn` | `rxnlist.cpp` | T, Δμ_Si, Δμ_Al, then 28 reaction rate tables |
| `data.cell` | `ucell.cpp` | cell dims/angles, Npos, per-position coords + neighbor template |
| `data.lattice` | `lattice.cpp` | aCells, bCells, surfacePlane |

`input/prot.*` are template/prototype copies (a "clean" starting set). The
`data.lattice.NxM` files are saved size presets you copy over `data.lattice`.

## 7. Output files (formats in the spec)

| File | Writer | Format |
|---|---|---|
| `start.msi`, `end.msi`, `step{i}.msi` | `writeMSI` | Cerius2 MSI (atoms + optional bonds) — the viewer's input |
| `start.xyz` | `writeXYZ` | comma-separated `Z,x,y,z` (note: comma, not standard XYZ) |
| `step{i}.dat`, `end.dat` | `writeData` | one CSV row: `time, Si_tot, Si(OH)0..4, Al_tot, Al..0..6` |
| `surfSi.out`, `surfAl.out` | `writeSurf` | CSV `x,y` of surface-exposed Si / Al |
| `results.dat` | — | **removed at startup, never written** (stale committed copy is from an older appending version) |

## 8. Build system

`Makefile`, GNU make, one binary:
- `CXX = g++`, `CXXFLAGS = -std=c++11 -Wall -Wextra -O3 -ffast-math`, `LDFLAGS = -lm`.
- **Static** `SOURCES` list (13 `.cpp`), `%.o: %.cpp $(HEADERS)` — every object
  depends on *all* headers, so any header edit rebuilds everything.
- Targets: `all` (default), `clean`. No install, no tests, no dependency
  auto-generation. `-ffast-math` is on in production (relevant to reproducibility —
  it relaxes IEEE float semantics; combined with the all-`float` state of the
  simulation, see spec §warts).
- Debug build = manually swap the commented `CXXFLAGS` line to `-ggdb`, `make clean && make`.
- Committed `.o` files and the `mckaol` binary are in the tree (build cruft; a Rust
  port's repo hygiene should ignore equivalents).

No external dependencies beyond libstdc++/libm. No package manager, no CMake.

## 9. The viewer (output-consumer, recon only)

`viewer/` is a separate React 19 + TypeScript + Vite + three.js app (its own
`package.json`). It parses MSI text via `src/utils/cerius2Parser.ts` into
`{atoms:[{id,element,position,label,charge}], bonds:[{from,to}]}` and renders atoms
as spheres, bonds as cylinders. `viewer/msi-schema.md` documents the MSI dialect.
It is the downstream consumer that `05-pgif-draft.md` proposes to feed with a new
property-graph format instead of MSI. Only skimmed here; full renderer recon already
exists at `projects/graph-viz-recon.md`.

## 10. Reading order for a newcomer

1. `common.hpp` — the state encoding and reaction list (the vocabulary).
2. `mckaol.cpp` — the whole program shape in 115 lines.
3. `lattice.hpp` (`LatticeSite`) + `ucell.hpp` — the data model.
4. `evtlist.cpp` + `actions.cpp` `DoEvent` — the KMC core loop.
5. `envrn.cpp` — how neighbor state selects rates (the science).
6. `rxnlist.cpp` — how input numbers become rates.
7. `output.cpp` `writeMSI` — the format the viewer eats.

Then read `02-model-spec.md` for the precise semantics and the wart list.
