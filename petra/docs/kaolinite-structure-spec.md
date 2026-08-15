# Kaolinite KMC Model — Structural Specification

Extracted 2026-08-12 from `/home/user/dissertation`. Authoritative sources:

- `legacy/cpp-model/data.cell`, `data.lattice` (inputs; byte-identical copy at `kmc-rs/data/golden/inputs/`)
- `legacy/cpp-model/ucell.cpp/.hpp`, `lattice.cpp/.hpp`, `bfsearch.cpp`, `envrn.cpp`, `evtlist.cpp`, `actions.cpp`, `output.cpp`, `common.hpp`
- `kmc-rs/crates/kaolinite/src/{cell.rs,build.rs,model_impl.rs,state.rs,reactions.rs}`, `kmc-rs/crates/kmc-engine/src/graph.rs`, `kmc-rs/crates/kmc-io/src/view.rs`, `kmc-rs/docs/REFORM_PLAN.md`
- Original dissertation-era C (pre-C++ refactor): git commit `8004955` (`mckaol.c`, `lattice.c`, `actions.c`) — cited where the C++ diverges from it.

Conventions in this document: "class" = `state / 100` (1=Al, 2=Si, 3=Si-O-Si O, 4=Si-O-Al2 O, 5=Al-OH-Al O); "occupancy" = `state % 100`. `EDGE = 9`, `WRONG = 99` (lattice.hpp:10-11).

---

## 1. Cell geometry

### 1.1 What data.cell declares

`data.cell` lines 1–3:

```
5.140  8.930  7.370        # a, b, c dimensions (angstroms)
-0.03141  -0.25038  0.0    # alpha, beta, gamma (radians)
26                         # Number of positions in unit cell
```

Parsed into `UnitCell::{A,B,C}` and `{Alpha,Beta,Gamma}` (ucell.cpp:65-67).

### 1.2 CRITICAL: the declared cell parameters are DEAD INPUT

Neither the lengths (5.140, 8.930, 7.370) nor the "angles" (-0.03141, -0.25038, 0.0) are ever used by any code path:

- The only accessor, `UnitCell::GetCellDimensions()` (ucell.cpp:18-26), has **no callers** anywhere in the codebase (verified by grep over all .cpp/.hpp).
- The three output writers (`writeMSI`, `writeXYZ`, `writeSurf`, output.cpp:54-56, 95-97, 144-146) each carry their own **hard-coded** Cartesian cell matrix and derive lengths from it, ignoring the file values.
- The Rust port documents this as wart "spec B5" (cell.rs:67-80; view.rs:67-72: "this matrix is presumably a once-correct render cell that outlived its inputs").

So the second line's three small numbers are *not* shear, *not* conventional angles, *not* anything: they are parsed, stored, and never read. Do **not** try to reconstruct geometry from them. The angles do not even reproduce the tilts of the hard-coded matrix (b-tilt of the matrix is atan(0.0262374/8.92893) ≈ 0.0029 rad, not 0.03141; c-tilt ≈ 0.179 rad, not 0.25038).

Also note: the simulation dynamics never touch coordinates at all — the lattice is purely topological. Coordinates exist only for output.

### 1.3 The effective fractional→Cartesian matrix (hard-coded in output writers)

output.cpp:144-146 (identical copies at :54-56 and :95-97); ported once as `CD` in kmc-io/src/view.rs:73-77:

```
CD = | 4.9725   -0.0262374   -1.3362  |
     | 0.0       8.92893     -0.30084 |
     | 0.0       0.0          7.384   |
```

Derived row-norm lengths (output.cpp:147-149; view.rs:82-88):

```
al = sqrt(4.9725² + 0.0262374² + 1.3362²) ≈ 5.148983
bl = sqrt(8.92893² + 0.30084²)            ≈ 8.933996
cl = 7.384
```

Note al ≈ 5.149 ≠ data.cell's a = 5.140 (etc.) — the matrix and the file disagree; the matrix wins because the file values are unused.

### 1.4 Exact per-axis coordinate formulas

For lattice site `i` at cell coords `(a, b)` stamped from unit-cell position `n` with stored position `(x, y, z)` (output.cpp:195-200, all three writers identical; view.rs:107-123):

```
fx = x/al + a          (pseudo-fractional a-coordinate)
fy = y/bl + b
fz = z/cl              (no c offset: single sheet, c index always 0)

X = fx*CD[0][0] + fy*CD[0][1] + fz*CD[0][2]
  = (x/al + a)*4.9725 + (y/bl + b)*(-0.0262374) + (z/cl)*(-1.3362)
Y = fy*CD[1][1] + fz*CD[1][2]
  = (y/bl + b)*8.92893 + (z/cl)*(-0.30084)
Z = fz*CD[2][2]
  = (z/cl)*7.384
```

i.e. cart = M·f with M = CD (upper triangular). The cell **vectors** are the columns of CD:

```
a_vec = ( 4.9725,     0.0,      0.0   )
b_vec = (-0.0262374,  8.92893,  0.0   )
c_vec = (-1.3362,    -0.30084,  7.384 )
```

**Are stored positions Cartesian or fractional?** They are Cartesian angstroms *within one cell* in an orthonormal-ish frame, converted to pseudo-fractional by dividing per-axis by the row norms `al, bl, cl` (which only approximate the true per-axis scale, since the rows are not the cell vectors — the columns are). This is a lossy/incoherent convention inherited from the original code; the golden outputs are defined by exactly the formulas above, so transcribe them literally.

**f32 subtlety for bitwise parity** (view.rs:90-122): the golden binary was compiled with `-ffast-math`; g++ replaced `x / al` with `x * (1.0f/al)` (reciprocal computed once). For unit-cell position 20 this differs by 1 ulp and changes one printed digit in `start.msi`. The Rust port deliberately uses reciprocal-multiply. Irrelevant unless you chase bitwise output parity.

---

## 2. The 26 unit-cell positions

Parsing (ucell.cpp:72-77): each position line is `id state x y z`, followed by exactly **6** neighbor entries `nbr_id da db dc`; `nbr_id = -1` means "no j-th neighbor" (slot unused). The `dc` column is parsed and stored but **never used** in neighbor resolution — the lattice is a single 2D sheet (cell.rs:10-12,31; every `dc` in the file is 0 anyway).

Expected neighbor count by class (`UnitCell::GetNumNeighbors`, ucell.cpp:36-48): **Al (1xx): 6; Si (2xx): 4; Si-O-Si (3xx): 2; Si-O-Al2 (4xx): 3; Al-OH-Al (5xx): 2** (any other class, incl. EDGE: -1).

Complete table, verbatim from `legacy/cpp-model/data.cell` (lines 5-160). Neighbor entries as `(nbr_id, da, db, dc)` in slot order j=0..5.

| pos | class | x | y | z | nbr[0] | nbr[1] | nbr[2] | nbr[3] | nbr[4] | nbr[5] | #real |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 100 (Al) | 1.76027 | 4.49982 | 2.84548 | (8,0,0,0) | (9,0,0,0) | (20,0,0,0) | (22,0,0,0) | (23,0,0,0) | (25,0,0,0) | 6 |
| 1 | 100 (Al) | 1.85474 | 7.3923 | 2.58231 | (9,0,0,0) | (13,0,0,0) | (18,0,0,0) | (21,0,0,0) | (23,0,0,0) | (24,0,0,0) | 6 |
| 2 | 100 (Al) | 4.35591 | 8.87916 | 1.9999 | (13,0,0,0) | (19,1,1,0) | (14,0,1,0) | (21,1,0,0) | (18,1,0,0) | (24,0,0,0) | 6 |
| 3 | 100 (Al) | 4.21668 | 2.93323 | 2.16386 | (19,1,0,0) | (14,0,0,0) | (20,1,0,0) | (8,1,0,0) | (22,0,0,0) | (25,0,0,0) | 6 |
| 4 | 200 (Si) | 0.308295 | 3.11457 | 0.225645 | (8,0,0,0) | (10,0,0,0) | (11,0,0,0) | (17,-1,0,0) | (-1,0,0,0) | (-1,0,0,0) | 4 |
| 5 | 200 (Si) | 0.412718 | 5.80163 | 0.328536 | (9,0,0,0) | (10,0,0,0) | (12,0,0,0) | (16,-1,0,0) | (-1,0,0,0) | (-1,0,0,0) | 4 |
| 6 | 200 (Si) | 2.70007 | 7.49498 | -0.417473 | (13,0,0,0) | (12,0,0,0) | (15,0,1,0) | (16,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 4 |
| 7 | 200 (Si) | 2.78958 | 1.60142 | -0.309399 | (14,0,0,0) | (11,0,0,0) | (15,0,0,0) | (17,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 4 |
| 8 | 400 (Si-O-Al2) | 0.546975 | 3.05081 | 2.01699 | (0,0,0,0) | (3,-1,0,0) | (4,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 3 |
| 9 | 400 (Si-O-Al2) | 0.735931 | 5.89814 | 1.7447 | (0,0,0,0) | (1,0,0,0) | (5,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 3 |
| 10 | 300 (Si-O-Si) | 0.407745 | 4.35517 | -0.45608 | (4,0,0,0) | (5,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 2 |
| 11 | 300 (Si-O-Si) | 1.28788 | 1.95757 | -0.55998 | (4,0,0,0) | (7,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 2 |
| 12 | 300 (Si-O-Si) | 1.16357 | 6.85128 | -0.80221 | (5,0,0,0) | (6,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 2 |
| 13 | 400 (Si-O-Al2) | 3.11279 | 7.51066 | 1.08072 | (1,0,0,0) | (2,0,0,0) | (6,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 3 |
| 14 | 400 (Si-O-Al2) | 3.40617 | 1.38387 | 1.30428 | (2,0,-1,0) | (3,0,0,0) | (7,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 3 |
| 15 | 300 (Si-O-Si) | 2.72493 | 0.0124087 | -0.91041 | (6,0,-1,0) | (7,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 2 |
| 16 | 300 (Si-O-Si) | 3.8338 | 6.45325 | -1.02691 | (5,1,0,0) | (6,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 2 |
| 17 | 300 (Si-O-Si) | 3.80894 | 2.31928 | -1.23535 | (4,1,0,0) | (7,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 2 |
| 18 | 500 (Al-OH-Al) | 0.581783 | 8.63121 | 1.85651 | (1,0,0,0) | (2,-1,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 2 |
| 19 | 500 (Al-OH-Al) | 0.11934 | 1.53515 | 4.22843 | (3,-1,0,0) | (2,-1,-1,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 2 |
| 20 | 500 (Al-OH-Al) | 0.636481 | 4.14859 | 4.00868 | (0,0,0,0) | (3,-1,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 2 |
| 21 | 500 (Al-OH-Al) | 0.527085 | 7.52431 | 4.02772 | (1,0,0,0) | (2,-1,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 2 |
| 22 | 500 (Al-OH-Al) | 3.01831 | 4.16281 | 1.29279 | (0,0,0,0) | (3,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 2 |
| 23 | 500 (Al-OH-Al) | 2.60559 | 5.87935 | 3.46514 | (0,0,0,0) | (1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 2 |
| 24 | 500 (Al-OH-Al) | 2.9487 | 8.68122 | 3.17508 | (1,0,0,0) | (2,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 2 |
| 25 | 500 (Al-OH-Al) | 2.894 | 3.27951 | 3.43088 | (0,0,0,0) | (3,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | (-1,0,0,0) | 2 |

Composition per cell: 4 Al (0-3), 4 Si (4-7), 4 Si-O-Al2 O (8, 9, 13, 14), 6 Si-O-Si O (10-12, 15-17), 8 Al-OH-Al O (18-25). Total 26.

### 2.1 Load-bearing slot ordering conventions

These slot positions are relied upon by the physics code — an input deck must preserve them:

- **400 sites:** `nbr[0]`, `nbr[1]` are the two Al; `nbr[2]` is the Si. FindPairs reads nbr[0]/nbr[1] as the Al pair (lattice.cpp:102-103, 111-112); environment code reads them likewise.
- **500 sites:** `nbr[0]`, `nbr[1]` are the two Al (FindPairs applies to state ≥ 400, same slots).
- **300 sites:** `nbr[0]`, `nbr[1]` are the two Si.
- **Si (200) sites:** `nbr[0]` is that Si's **400** (Si-O-Al2) oxygen; nbr[1..3] are its three 300 oxygens. envrn.cpp:124 comment: "Si-O-Al is Si nbr[0]"; used by Check300 (envrn.cpp:124-125), PopulateSolid's chain guard (lattice.cpp:198-199), and TerminateLattice pass 2 (lattice.cpp:245-246).
- **Al (100) sites:** 6 oxygens, two 400-class + four 500-class, in **no particular slot order** (site 2 interleaves them: 13(400), 19(500), 14(400), 21, 18, 24).
- Real neighbors are packed **first**; `-1` entries pad the tail. `CountNbrs` counts the leading run of non-negative entries and stops at the first -1 (lattice.cpp:87-93; graph.rs:105-111 preserves the quirk) — a real neighbor after a -1 would be silently dropped from the count.

---

## 3. Pairs (the double bridge)

`Lattice::FindPairs()` (lattice.cpp:96-121; Rust `find_pairs`, build.rs:196-230).

**Rule:** every adjacent Al-Al pair in kaolinite's octahedral sheet is bridged by exactly *two* oxygen sites (edge-sharing octahedra). `pair` links those two bridging oxygens to each other.

**Algorithm as implemented:**
1. Scan all sites `o1` with `state >= 400` (i.e. 400 *and* 500 classes) not yet paired.
2. Let `al1 = o1.nbr[0]`, `al2 = o1.nbr[1]` (the two Al). If either is -1 (open boundary), skip — `o1` stays unpaired.
3. Scan `al1`'s neighbor slots j = 0..5 in order:
   - if `al1.nbr[j] < 0`: stop searching ("no neighbor, must be edge") — `o1` stays unpaired;
   - else let `o2 = al1.nbr[j]`; if `{o2.nbr[0], o2.nbr[1]} == {al1, al2}` (set comparison, lattice.cpp:113), then `pair[o1] = o2; pair[o2] = o1` and stop.

**The pair partner is NOT in the site's own nbr[] list.** It is a separate relation on top of the graph: `o1` and `o2` are both neighbors *of the two Al*, not of each other. (E.g. pos 8's nbr list is {0, 3, 4}; its pair is pos 20, absent from that list.) `pair` is stored in its own field (`LatticeSite::pair`, lattice.hpp:36; `Structure::pair: Vec<Option<SiteId>>`, build.rs:72).

**Known quirk (ported deliberately, build.rs:178-183):** the scan does not exclude `o1` itself; `o1` trivially matches its own criterion. Whether it self-pairs transiently depends on whether the Al's slot order presents `o1` before the true partner; in the golden cell the final result is correct either way (a later scan of the partner overwrites), and the golden lattice has **zero** self-paired sites (verified computationally, below).

**Derived pair topology** (computed by replaying the exact algorithm on data.cell over a 20x3 lattice; interior cells):

| oxygen | partner | partner cell offset | bridge type |
|---|---|---|---|
| pos 8 (400) | pos 20 (500) | (0,0) | Al0-Al3(a-1) |
| pos 9 (400) | pos 23 (500) | (0,0) | Al0-Al1 |
| pos 13 (400) | pos 24 (500) | (0,0) | Al1-Al2 |
| pos 14 (400) | pos 19 (500) | (+1,0) | Al2(b-1)-Al3 |
| pos 18 (500) | pos 21 (500) | (0,0) | Al1-Al2(a-1) |
| pos 22 (500) | pos 25 (500) | (0,0) | Al0-Al3 |

So each 400 pairs with a 500, and the remaining 500s pair among themselves (18↔21, 22↔25). Every Al participates in 6 bridges = 3 double bridges... (4 Al x 6 O)/2 = 12 bridging O per cell = 6 double bridges per cell, matching the 6 rows above.

**Boundary effect:** in a 20x3 surface_plane=0 lattice, positions 14 and 19 in the b=0 row stay **unpaired** (their Al templates reach db=-1 off the open boundary): 40 unpaired 400/500 sites total (20 a-cells x 2 positions).

---

## 4. Lattice build

### 4.1 data.lattice

```
20	3	0
acells  bcells  ac/bc surface plane
```

`aCells=20, bCells=3, SurfacePlane=0` (lattice.cpp:27). Alternates in repo: `data.lattice.{50x10,100x10,500x10}`. kmc-rs golden uses the same 20 3 0.

**SurfacePlane semantics** (lattice.cpp:138-139; build.rs:25-30):
- `0` = **ac surface plane**: **b is the open (surface-normal) axis**, a wraps periodically.
- `1` (or any nonzero — the C++ tests truthiness) = **bc surface plane**: a open, b periodic.

### 4.2 c direction

There is **no c tiling — exactly 1 cell in c**. `LatticeSite` has only `a, b` cell coords (lattice.hpp:33); the template's `dc` offsets are parsed but never used in `GetNeighbor` (lattice.cpp:130-153 reads only `.a` and `.b`); all `dc` in data.cell are 0. The simulated system is a single (001) sheet.

### 4.3 Site indexing

`i = a * (bCells * npos) + b * npos + n` with npos = 26 (lattice.cpp:35, build.rs:85). `Num_Sites = aCells * bCells * npos` = 20*3*26 = **1560** for the golden lattice.

### 4.4 Neighbor resolution (GetNeighbor, lattice.cpp:130-153; resolve_neighbor, build.rs:127-165)

For site (a,b,n), slot j with template (tn, ta, tb):
1. If `tn < 0` → no neighbor (C++ returns the negative value itself; Rust `None`).
2. Let `na = a + ta`, `nb = b + tb`.
3. **Open-boundary check:** if `(na == aCells || na < 0) && SurfacePlane`, or `(nb == bCells || nb < 0) && !SurfacePlane` → neighbor = -1/None. (Tests `== aCells`, not `>=` — safe because template offsets are at most ±1.)
4. **Periodic wrap on both axes otherwise:** `na >= aCells → 0`; `na < 0 → aCells-1`; same for b. (The wrap on the open axis is unreachable after step 3.)
5. `nbr = na * bCells * npos + nb * npos + tn`.

Neighbors are resolved once at build time and never recomputed.

---

## 5. Initial state (PopulateSolid)

lattice.cpp:159-214; build.rs:265-329.

1. **Fill fraction from chemical potentials** (from data.rxn: dmSi, dmAl; golden = -1.0, -1.0):
   - `dmSi + dmAl > 0.5` → frac = 0.3 ("supersaturated")
   - `dmSi + dmAl < -0.5` → frac = 0.7 ("undersaturated") ← golden case
   - otherwise → frac = 0.5 ("near equilibrium")
2. **Filled depth** `top = (int)(cells_along_open_axis * frac)` where the open axis is a if SurfacePlane else b. **f32 truncation is load-bearing** (build.rs:250-260): `3 * 0.7f32 = 2.0999999 → 2` filled layers in the golden run, but `20 * 0.7f32 = exactly 14.0 → 14`. Reproduce `int = int * float` truncation exactly.
3. **Fill:** for each layer z in `0..top` (from the row-0 side of the open axis) and every cell x along the periodic axis, every one of the 26 positions gets `state++`:
   - 100 → **101** (Al, fully coordinated Al(OH,H2O)0)
   - 200 → **201** (Si(OH)0)
   - 300 → **301** (Si-O-Si bridge)
   - 400 → **401** (Si-O-Al2)
   - 500 → **501** (Al-OH-Al)

   So the solid starts as fully-occupied bulk (occupancy 1 in every class); everything above `top` stays empty (100/200/300/400/500).
4. **Guard (modern C++/Rust only — NOT in the original dissertation C, git 8004955 `lattice.c` increments unconditionally):** a 300 site is only promoted to 301 if its full chain exists: `nbr[0]` and `nbr[1]` (both Si) present, and each Si's `nbr[0]` (its 400 O) present (lattice.cpp:189-205; build.rs:302-316). With data.lattice=20 3 0 filling starts at b=0, where several sites reach off the open boundary, so this guard fires.

---

## 6. Surface termination (TerminateSurface)

lattice.cpp:262-381; Rust terminate_surface, build.rs:335-454. Runs after FindPairs+PopulateSolid, **before** TerminateLattice (order fixed by mckaol.cpp:45-49). Three passes over all sites, in order; each pass's writes are visible to later iterations of the same pass (sequential, single sweep in index order).

### Pass 1 — demote oxygens missing a cation (lattice.cpp:266-304)

For each site `i` with occupied state (`%100 != 0`), class ≥ 3, and not EDGE; for each existing neighbor `i2` (slots 0..5 in order) whose occupancy is 0 (`sites[i2].state % 100 == 0` — an empty cation site), step `i`'s state down per the map (applied once per empty neighbor found, cumulatively within the slot loop):

**As implemented (the dead-ternary wart):** the code assigns `type = sites[i2].state % 100` inside the `== 0` test (lattice.cpp:272), so `type` is always 0, making both `type == 2` ternaries dead:

| current | → new (as implemented) |
|---|---|
| 401 | **406** (always; the `type==2 ? 404` arm is dead) |
| 404 | 409 |
| 406 | **409** (always; the `type==2 ? 408` arm is dead) |
| 409 | 400 |
| 408 | 400 |
| 301 | 303 |
| 303 | 300 |
| 501 | 503 |
| 503 | 500 |
| anything else | unchanged |

Rust reproduces this verbatim (build.rs:354-380, "ported AS-WRITTEN, quirk included").

**Intended/corrected reading (kmc-rs docs/REFORM_PLAN.md R3, lines 102-129):** the author evidently meant to branch on the **class of the missing cation** (`sites[i2].state / 100`: 2 = empty Si site, 1 = empty Al site):

| current | missing cation is Si (empty 200) | missing cation is Al (empty 100) |
|---|---|---|
| 401 | **404** (drop the Si → HO<Al2) | **406** (drop one Al → Si-OH-Al) |
| 406 | **408** (drop the remaining Al's counterpart → Si-OH) | **409** (→ HO-Al) |

(remaining rows unchanged — i.e. the ternaries as literally written, with `type` computed as the empty neighbor's class digit instead of its occupancy.) REFORM_PLAN's confidence is explicitly **low on what the correction should be** ("Do not guess the intended map ... Blocked on Victor['s] domain read", REFORM_PLAN.md:122-129); high confidence only that current behavior is accidental. The corrected mode changes the *initial* hydroxylated surface and breaks the start.msi bitwise gate. For a faithful input deck, transcribe the as-implemented map; flag the corrected map as an unresolved variant.

### Pass 2 — hydroxylate oxygens next to minimally-coordinated cations (lattice.cpp:307-343)

For each site `i` that is an **occupied cation at occupancy exactly 1** (state 101 or 201; condition: class ≤ 2 and `%100 == 1`), with `type = class of i` (1=Al, 2=Si — here the ternaries are **live**); for each existing neighbor `i2` with occupancy ≠ 1, map `i2`:

| current i2 | → if cation i is Si (type==2) | → if cation i is Al (type==1) |
|---|---|---|
| 300 | 303 | 303 |
| 400 | **408** | **409** |
| 408 | 406 | 406 |
| 409 | **406** | **404** |
| 404 | 401 | 401 |
| 406 | 401 | 401 |
| 500 | 503 | 503 |
| 503 | 501 | 501 |
| else | unchanged | unchanged |

### Pass 3 — count terminal OH into cation states (lattice.cpp:346-380)

For each occupied, non-EDGE cation (class 1 or 2, occupancy ≥ 1): for each existing neighbor,
- **Al** sites: `state++` for each neighbor in state **503 or 409**;
- **Si** sites: `state++` for each neighbor in state **303 or 408**.

(E.g. an Al with three 503/409 neighbors becomes 104 = Al(OH,H2O)3.)

---

## 7. EDGE / TerminateLattice

lattice.cpp:218-258; Rust terminate_lattice, build.rs:464-528. Runs **after** TerminateSurface (order matters: the neighbor-count test reads post-termination states).

### Pass 1 — mark boundary rows EDGE

The two extreme rows along the **open** axis — row index `top = (openAxisCells - 1)` and row 0 (both rows, all cells along the periodic axis, all 26 positions) — get `state = EDGE (9)` **iff** `GetNumNeighbors(state) != CountNbrs(site)`, i.e. the expected count for the site's current class (Al 6 / Si 4 / 3xx,5xx 2 / 4xx 3) differs from the actual leading-run count of resolved neighbors. Sites in the boundary rows whose neighbor template happens to be fully satisfied are *not* marked.

### Pass 2 — revert orphaned Si-O-Si bridges (modern C++/Rust only; NOT in original C `terminateLattice`, git 8004955)

Any site in state **301** whose chain — `si1 = nbr[0]`, `si2 = nbr[1]`, `sio1 = si1.nbr[0]`, `sio2 = si2.nbr[0]` — contains an EDGE site reverts to **300** (lattice.cpp:239-257).

### What EDGE means for dynamics

- **EDGE sites never react:** state 9 falls through the event-window selection into the adsorption/desorption window (evtlist.cpp:26-28) but matches no reaction's reactant, so proposes no events. (state.rs:78-81: "EDGE sites never react ... they are the bulk crystal below the simulated surface".)
- **EDGE counts as "occupied"** by the `ISOCC` macro (`9 % 100 = 9 > 0`, lattice.hpp:20; state.rs:95-100) — relevant to BFS (below) — but:
  - adsorption's is-active check explicitly requires an occupied neighbor that is **not** EDGE (envrn.cpp:295-307);
  - every environment-classification routine (`Check100/200/300/400/500`) **returns -1 ("ran into lattice edge") if it touches an EDGE site** (envrn.cpp:56-58, 89-91, 128-134, 172-176, 225-227, ...), which `CreateEventList` treats as fatal (evtlist.cpp:45-52 returns nullptr → run aborts; Rust: `KaolError::InvalidEnvironment`, model_impl.rs:168-174). So EDGE is a hard wall: any reactive site whose environment reaches an EDGE kills the run — the buffer of empty rows between the solid and the frozen boundary is what keeps this from happening.
- EDGE sites are skipped by all output writers (output.cpp:165, 231).

---

## 8. Cluster removal (RemoveUnattachedClusters / BFS)

lattice.cpp:70-83; bfsearch.cpp:13-38.

- **Anchor/root:** BFS starts from **site index 0** — hard-coded (`ColorNodes(this, 0)`, lattice.cpp:73), i.e. cell (0,0), position 0 (an Al in the filled bottom layer). The root is enqueued unconditionally, without an occupancy check on the root itself.
- **Connectivity:** from a dequeued site, each of the 6 `nbr[]` entries is followed iff it exists (`>= 0`) **and** is occupied (`ISOCC`: `state % 100 > 0`) (bfsearch.cpp:27-29). Note **EDGE (9) passes ISOCC**, so EDGE sites are traversable/conductive in the BFS. Bond state is irrelevant — any occupied adjacency conducts (a hydrolyzed bridge, e.g. 302, still connects).
- **After BFS:** every site left UNREACHABLE whose state is not 9 is reset to its empty class value: `state = (state / 100) * 100` (→ 100/200/300/400/500) (lattice.cpp:74-76). EDGE sites are preserved. A site found in ENQUEUED color is an internal-error condition (prints "bad color", returns false → run aborts).
- Note the reset does **not** touch `pair`/`lostal` bookkeeping — only the state.

### When it runs — three different answers across generations

1. **Original dissertation C (git 8004955, actions.c:69-86):** `clusters()` ran **after every desorption event (reactions 20-23) and after every diffusion event** — i.e. effectively after every event that can detach material.
2. **Current C++ (legacy/cpp-model/actions.cpp:106-117):** the desorption cases no longer call it; only the `default` branch (Diffuse, reactions ≥ 24) calls `RemoveUnattachedClusters()`. But diffusion is **never active** (`IsActive` returns FALSE for all diffusion reactions, envrn.cpp:312-315), so diffusion events are never proposed and **cluster removal never runs at all** in the current C++. (The top-level CLAUDE.md's "3. RemoveUnattachedClusters() - BFS cleanup [every step]" describes the original design, not the current code.) The `LatticeSite::color` field is dead in practice.
3. **kmc-rs:** neither diffusion nor BFS is ported (RUST_TOUR.md:341: "diffusion + BFS dead code — not ported"; graph.rs:9 calls color "the dead BFS `color`"). Reactions array holds only 24 entries (reactions.rs:11-17, N_DES = 24) — reactions 24-27 (diffusion) are omitted along with a genuine out-of-bounds read in the C++ diffusion table copy (rxnlist.cpp:110-116 copies with a stale `numRates`; documented as spec B7).

**For a new engine:** decide explicitly. Bitwise parity with the current C++/Rust golden = no cluster removal. Physical fidelity to the dissertation-era model = run BFS cleanup after every desorption.

---

## Appendix: reaction/window facts needed to interpret the above

- 28 reactions in C++ (24 live in Rust): 0-15 hydrolysis (8 forward/reverse pairs), 16-19 adsorption (16 Al→107, 19 Si→205; 17/18 are never-active cross variants), 20-23 desorption (20/23 Al 107→100, 21/22 Si 205→200), 24-27 diffusion (dead) (reactions.rs:5-17, common.hpp:57-112).
- Event windows by site state (evtlist.cpp:14-29; model_impl.rs:127-162): empty O (`%100==0 && >200`) skipped; `>500` → R14-15; `>400` → R2-13; `>300` → R0-1; everything else (incl. empty 100/200) → adsorb/desorb window.
- 404/405 use `Check500` for environment despite being class 4 (envrn.cpp:27-33).
- data.rxn wart R4 (REFORM_PLAN.md:133-157): R13's reactant token is `40100` (typo, likely `405`), so R13 never fires in legacy.
- State catalog verbatim in common.hpp:1-56 (mirrored in state.rs): 100-107+199, 200-205+299, 300-303, 400-410, 500-503, EDGE=9, WRONG=99.
