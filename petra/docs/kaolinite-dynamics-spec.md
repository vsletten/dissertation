# Kaolinite KMC Dynamics Specification

Extracted from the dissertation repo for authoring a declarative input deck (petra).
Authoritative sources: the Rust port `kmc-rs/crates/kaolinite/` (bitwise-parity-tested
against the C++ over 20,000 steps) and the legacy C++ in `legacy/cpp-model/`.
All paths below are absolute repo paths; line numbers are cited for every nontrivial claim.

Abbreviations:
- `state.rs` = /home/user/dissertation/kmc-rs/crates/kaolinite/src/state.rs
- `environment.rs` = /home/user/dissertation/kmc-rs/crates/kaolinite/src/environment.rs
- `mechanisms.rs` = /home/user/dissertation/kmc-rs/crates/kaolinite/src/mechanisms.rs
- `reactions.rs` = /home/user/dissertation/kmc-rs/crates/kaolinite/src/reactions.rs
- `model_impl.rs` = /home/user/dissertation/kmc-rs/crates/kaolinite/src/model_impl.rs
- `envrn.cpp` / `actions.cpp` / `rxnlist.cpp` / `evtlist.cpp` = /home/user/dissertation/legacy/cpp-model/…
- `data.rxn` = /home/user/dissertation/legacy/cpp-model/data.rxn
- `REFORM_PLAN` = /home/user/dissertation/kmc-rs/docs/REFORM_PLAN.md

---

## 1. State vocabulary

A site state is a 3-digit integer: hundreds digit = site class, last two digits =
occupancy/coordination/protonation. `state % 100 == 0` means the site is **vacant**.
The table is verbatim from `common.hpp` via `state.rs:13-41`.

### Class 1xx — Al cation sites (vacant = 100; expects 6 neighbors, state.rs:119-127)

| code | meaning |
|---|---|
| 100 | empty Al site |
| 101 | Al with 0 OH/H2O ligands — Al(OH,H2O)0 |
| 102 | Al(OH,H2O)1 |
| 103 | Al(OH,H2O)2 |
| 104 | Al(OH,H2O)3 |
| 105 | Al(OH,H2O)4 |
| 106 | Al(OH,H2O)5 |
| 107 | Al(OH,H2O)6 — fully solvated, the adsorbed/desorbable form |
| 199 | Si(OH)4 sitting on an Al site (wrong-cation; never produced — cross adsorption rxn 18 is gated off) |

The last digit of an occupied 1xx literally counts OH/H2O ligands; mechanisms do
`state++`/`state--` on the Al neighbor to add/remove one ligand (e.g. actions.cpp:298,
mechanisms.rs:65-67).

### Class 2xx — Si cation sites (vacant = 200; expects 4 neighbors)

| code | meaning |
|---|---|
| 200 | empty Si site |
| 201 | Si(OH)0 — fully lattice-bonded Si |
| 202 | Si(OH)1 |
| 203 | Si(OH)2 |
| 204 | Si(OH)3 |
| 205 | Si(OH)4 — fully hydrolyzed, the adsorbed/desorbable form |
| 299 | Al(OH,H2O)6 on a Si site (wrong-cation; never produced) |

### Class 3xx — Si-O-Si bridging oxygen (vacant = 300; expects 2 neighbors)

| code | meaning |
|---|---|
| 300 | empty (both Si absent) |
| 301 | Si-O-Si intact bridge |
| 302 | Si-OH HO-Si — hydrolyzed bridge (both Si present, bond broken) |
| 303 | Si-OH — one Si present (dangling hydroxyl) |

### Class 4xx — Si-O<Al2 bridging oxygen (vacant = 400; expects 3 neighbors: nbr[0]=Al, nbr[1]=Al, nbr[2]=Si)

| code | meaning |
|---|---|
| 400 | empty (no cations attached) |
| 401 | Si-O<Al2 — intact triple bridge (Si + 2 Al) |
| 402 | Si-OH HO<Al2 — Si-side hydrolyzed (Si-O bond broken, O still bridges both Al) |
| 403 | Si-OH HO-Al H2O-Al — fully hydrolyzed (all three cation-O bonds broken) |
| 404 | HO<Al2 — Si absent, O bridges both Al (behaves as an Al-OH-Al; routed to Check500) |
| 405 | HO-Al H2O-Al — Si absent, Al-O-Al bridge hydrolyzed (routed to Check500) |
| 406 | Si-OH-Al — one Al absent, Si-O-Al bridge intact |
| 407 | Si-OH HO-Al — one Al absent, Si-O bond broken |
| 408 | Si-OH — only the Si present |
| 409 | HO-Al — only one Al present |
| 410 | Si-OH-Al HO-Al — Si-O-Al intact on one Al, other Al carries the proton (produced by R4/R9 and by AdsorbAl on a 406) |

### Class 5xx — Al-OH-Al bridging oxygen (vacant = 500; expects 2 neighbors: nbr[0]=Al, nbr[1]=Al)

| code | meaning |
|---|---|
| 500 | empty |
| 501 | Al-OH-Al intact bridge |
| 502 | Al-OH H2O-Al — hydrolyzed bridge |
| 503 | Al-OH — one Al present |

### Sentinels

| code | meaning |
|---|---|
| 9 (EDGE) | frozen boundary site: never reacts, skipped in output; represents bulk crystal below the surface (state.rs:78-81). **Trap:** `EDGE.class_code() == 0` and `ISOCC(EDGE)` is true (9 % 100 = 9 > 0), so every EDGE check is explicit (state.rs:95-105, 143-148). |
| x99 (WRONG) | wrong-cation marker (199/299), reachable only via disabled cross reactions |

### Hidden per-site fields carrying dynamical information

Beyond `state`, the C++ `LatticeSite` carries (legacy/cpp-model/CLAUDE.md struct;
Rust: model_impl.rs:88-95):

1. **`lostal: Option<SiteId>`** (C++ `int lostal`, -1 = none) — for a 4xx oxygen in
   state 410, *which of its two Al neighbors carries the proton* (i.e. which Al was
   "lost" from the bridge). Set by R4's coin flip (mechanisms.rs:154-170), by R9's
   coin flip (mechanisms.rs:226-244), and by `AdsorbAl` when it converts a 406→410
   (mechanisms.rs:336-339, actions.cpp:149-152). Cleared (to None/-1) by R5, R8, and
   `DesorbAl`'s 410→406 arm.
   **Legality/rate impact:** none — no rate and no `is_active`/`check_env` test reads
   `lostal`. It affects *only which neighbor site* R5 and R8 mutate (mechanisms.rs:179,
   215-219). It is state the rewrite rules read and write, so a declarative engine must
   either carry it as an auxiliary site attribute or encode it in the state (e.g. split
   410 into 410a/410b by which Al holds the proton).

2. **`pair: Option<SiteId>`** — the double-bridge partner oxygen of a 4xx/5xx site
   (kaolinite has pairs of bridging O between the same two Al). Built once by
   `FindPairs` at setup, **never mutated by dynamics** (mechanisms.rs:14-16). Read only
   by `check_env` (Check400/Check500 use the partner's state to bucket the rate,
   environment.rs:279, 327). With the shipped flat rate tables, its only dynamical
   effect is the fatal `-1` when the partner is missing/EDGE.

3. **`color`** (C++ BFS marking field) — used only by `RemoveUnattachedClusters`,
   which the live main loop **never calls** (it is reachable only from the dead
   diffusion arm, actions.cpp:114-117; mckaol.cpp:61-99 has no call). Irrelevant to
   dynamics — except that the `is_active` out-of-bounds read `nbr[6]` aliases the
   `color` field (REFORM_PLAN:41-43), which is what makes the phantom deterministic.

---

## 2. Event enumeration

Per step the entire event list is rebuilt from scratch over all sites
(evtlist.cpp:14-57; per-site body in model_impl.rs:139-184).

### Window selection by site state (evtlist.cpp:15-29, model_impl.rs:147-162)

Reaction indices: 0–15 hydrolysis (even = forward), 16–19 adsorption, 20–23
desorption, 24–27 diffusion (dead). Constants `N300 = 2`, `N400 = 14`, `NHYD = 16`,
`NADS = 20`, `NDES = 24`, `NRXN = 28` (reactions.rs:30-40).

| site state | window | candidate reactions |
|---|---|---|
| empty oxygen: `state % 100 == 0 && state > 200` (300/400/500) | — | none, site skipped |
| `> 500` (501–503) | 14..16 | R14, R15 |
| `> 400` (401–410) | 2..14 | R2–R13 |
| `> 300` (301–303) | 0..2 | R0, R1 |
| else (all 1xx and 2xx, incl. **empty** 100 and 200, and EDGE=9) | 16..28 (Rust caps at 24) | adsorption/desorption (16–23); diffusion 24–27 always inactive |

Within the window, an event `(site, rxn)` is proposed iff **both**:
1. `state == reactions[rxn].reactant` (exact match, evtlist.cpp:32-33), and
2. `IsActive(site, rxn)` (evtlist.cpp:34).

Then `env = CheckEnv(site)` picks the rate: `rate = reactions[rxn].rate[env]`; if
`env < 0 || env >= nrates` the run **aborts** ("invalid environment",
evtlist.cpp:44-52, model_impl.rs:167-175).

Note EDGE sites fall in the `else` window but match no reactant (9 ≠ 100/107/200/205…),
so they propose nothing.

### `is_active` — exact role (envrn.cpp:259-318, environment.rs:51-119)

| rxn | rule |
|---|---|
| even hydrolysis (0,2,4,…,14) — forward | "at the surface" test (below) |
| odd hydrolysis (1,3,…,15) — reverse | always active |
| 16 (adsorb Al→Al site), 19 (adsorb Si→Si site) | active iff ≥1 neighbor is occupied (`ISOCC`) and not EDGE (envrn.cpp:295-307) |
| 20 (desorb Al), 22 (desorb Si) | always active |
| 17, 18, 21, 23 (cross-cation) and 24–27 (diffusion) | **never active** — this is how those reactions are switched off (environment.rs:115-118) |

#### Forward-hydrolysis surface test — LEGACY (as-implemented, parity-bearing)

C++ loop (envrn.cpp:263-289):
```c
for (i = 0; (nbr = sites[site].nbr[i]) >= 0 && !result; i++)
  for (j = 0; (nbr2 = sites[nbr].nbr[j]) >= 0 && j < 6 && !result; j++)
    if (sites[nbr2].state in {303,404,405,406,408,409}) result = TRUE;
```
The inner loop reads `nbr[j]` **before** checking `j < 6`. When a neighbor `nbr` has
all six slots filled (every Al does), the loop reaches `j == 6` and reads
`sites[nbr].nbr[6]` — out of bounds, aliasing `color`. Under the golden build
(`g++ -O3 -ffast-math`) this deterministically makes `result = TRUE`
(environment.rs:56-96; REFORM_PLAN:35-47).

**Legacy effective rule** (verified across all 20,000 parity steps,
environment.rs:68-74): forward hydrolysis is active at `site` iff for some neighbor
`nbr` (scanned in order, stopping at the first absent slot):
- some *present* neighbor of `nbr` (scanned in slot order up to the first absent
  slot) has state in the hydrolyzed set **{303, 404, 405, 406, 408, 409}**, OR
- **all six** of `nbr`'s neighbor slots are present (the phantom: the inner loop
  reaches j == 6). Since topology is fixed, this phantom condition is static per site.

Consequence: every site adjacent to a fully-6-coordinated neighbor (any Al) is always
"active" for forward hydrolysis. At step 0 this inflates the event list from a
bounds-clean 180 to the golden 660 (REFORM_PLAN:46-47); it roughly quadruples the
legal forward-hydrolysis events and partly defeats the surface-only premise
(REFORM_PLAN:53-60).

#### Forward-hydrolysis surface test — CORRECTED (intended, spec A6.1)

Active iff at least one **second-neighbor** (a neighbor of a neighbor) has state in
the hydrolyzed set **{303, 404, 405, 406, 408, 409}** (REFORM_PLAN:48-51). Note the
set is "dangling/partially-hydrolyzed oxygens": 303 (Si-OH), 404 (HO<Al2),
405 (HO-Al H2O-Al), 406 (Si-OH-Al), 408 (Si-OH), 409 (HO-Al). Notably 302, 402, 403,
407, 410, 502, 503 are NOT in the set.

For the new engine: this is a predicate over the 2-hop neighborhood of the center site.

### `check_env` — what it computes (envrn.cpp:8-255, environment.rs:154-352)

Returns an integer "environment bucket" indexing the reaction's rate table. Dispatch
(envrn.cpp:11-42):
- class 1: 0 if state==100, else `Check100`
- class 2: 0 if state==200, else `Check200`
- class 3: `Check300`
- class 4: `Check500` if state ∈ {404, 405} (they are chemically Al-OH-Al bridges),
  else `Check400` — a deliberate quirk (environment.rs:172-178)
- class 5: `Check500`

Bucket formulas (all preserved warts included):
- **Check100** (envrn.cpp:47-76): over all 6 neighbors, x = #(502), y = #(states in
  {403,405,407,409,410}); returns `(x + y) / 2` (integer division). WART: the header
  documents `x*3 + y` (REFORM_PLAN:161-167). Returns -1 if any neighbor is EDGE.
- **Check200** (envrn.cpp:79-105): over 4 neighbors, x starts at 1 and is set to 0 by
  any 408 neighbor; y = #(302); returns `x + y`. WART: header says `x*4 + y`.
  -1 on any EDGE neighbor.
- **Check300** (envrn.cpp:108-156): reaches through both Si (`nbr[0]`, `nbr[1]`) to
  each Si's `nbr[0]` (its 4xx oxygen). x = #(those two oxygens in {402,403,407,408});
  y = (si1.state-201)+(si2.state-201) − x, minus 2 more if the center is > 301;
  returns `5*x + y`. −1 if any of the chain (si1, si2, their nbr[0]s) is missing or EDGE.
- **Check400** (envrn.cpp:159-212): needs al1=nbr[0], al2=nbr[1], si=nbr[2], and the
  `pair` partner p, all present and non-EDGE, else −1; also needs si.nbr[1], si.nbr[2]
  present. x = #(si.nbr[1], si.nbr[2] with state > 301). y = (al1−101)+(al2−101),
  minus 2 if center ∈ {403,405}; minus 1 if center ∈ {407,410} or p == 502; minus 1
  more if p == 502 (i.e. p==502 subtracts 2 total, and stacks with the 407/410 term).
  Returns `x*5 + y` if center ∈ {406,407}, else `x*10 + y`.
- **Check500** (envrn.cpp:215-255): needs al1=nbr[0], al2=nbr[1], and pair p, present,
  non-EDGE, else −1. x = 1 if p ∈ {502,403,405,410} else 0. y = (al1−101)+(al2−101),
  minus 1 if p ∈ {502,403,405}, minus 2 if center ∈ {502,405}. Returns `x*9 + y`.

**When −1 (fatal):** any required neighbor / second-neighbor / pair partner missing or
EDGE (per the lists above). In the event-list builder, −1 (or an index ≥ nrates)
aborts the whole run (evtlist.cpp:44-52). In practice the initial `TerminateLattice`
freezes boundary sites so live sites never trip this.

**Does check_env affect dynamics beyond legality?** With the shipped `data.rxn`,
**no**: every rate table is flat (all buckets within a reaction identical —
see §5), so the bucket value never changes a rate (reactions.rs:19-22,
REFORM_PLAN:219-227). Its only live effects are (a) the fatal −1/out-of-range abort
gate, and (b) bookkeeping. The moment a non-flat deck is authored, the (warty)
formulas above select rates — and the Check100/Check200 formulas disagree with their
own documentation (REFORM_PLAN R5), so a new engine should treat the *intended*
formulas as open questions ("blocked on Victor", REFORM_PLAN:180-184).

---

## 3. The 16 hydrolysis mechanisms (R0–R15)

Ground truth: mechanisms.rs:114-312 ≡ actions.cpp:291-480. Conventions:
- "center" = the reacting oxygen site (the one whose state matched the reactant).
- Neighbor slot conventions: 3xx center: nbr[0], nbr[1] = the two Si. 4xx center:
  nbr[0], nbr[1] = the two Al, nbr[2] = the Si. 5xx center: nbr[0], nbr[1] = the two Al.
- "+1"/"−1" on a cation = one protonation step (`state++`/`state--`); the mechanism
  never verifies the cation's current state (except R10/R11's ==100 test).
- No mechanism checks any guard at application time beyond what's listed; legality
  was fully decided at enumeration time (reactant match + is_active).
- Forward = even index; Reverse(k) = exact inverse of reaction k.

| R | name / meaning | center | other sites touched | random draw |
|---|---|---|---|---|
| **R0** | Si-O-Si hydrolysis, forward: =Si-O-Si= + H2O → 2(=Si-OH) | 301 → 302 | both Si (nbr[0], nbr[1]): each +1 | — |
| **R1** | Si-O-Si condensation (reverse of R0) | 302 → 301 | both Si: each −1 | — |
| **R2** | Si-O<Al2 hydrolysis, Si side, forward: =Si-O<Al2 + H2O → =Si-OH + Al-O(H)-Al | 401 → 402 | Si (nbr[2]): +1. **Al untouched** | — |
| **R3** | reverse of R2 | 402 → 401 | Si (nbr[2]): −1 | — |
| **R4** | Si-O<Al2 hydrolysis, Al side, forward: =Si-O<Al2 + H2O → =Si-OH-Al + Al-OH | 401 → 410 | ONE Al gains +1 — chosen by coin (see draw); `lostal[center] := that Al` | one `ran2()`: r < 0.5 → nbr[0] (Al1) gains, else nbr[1] (Al2) gains (actions.cpp:336-343) |
| **R5** | reverse of R4 | 410 → 401 | the Al recorded in `lostal[center]`: −1; `lostal[center] := none` | — (deterministic undo via lostal) |
| **R6** | second-stage: =Si-OH + Al-O(H)-Al + H2O → =Si-OH + Al-OH + Al-H2O | 402 → 403 | both Al (nbr[0], nbr[1]): each +1. Si untouched | — |
| **R7** | reverse of R6 | 403 → 402 | both Al: each −1 | — |
| **R8** | =Si-OH-Al + Al-OH + H2O → =Si-OH + Al-OH + Al-H2O | 410 → 403 | Si (nbr[2]): +1; the **other** Al — the one that is NOT `lostal[center]` (`lostal==nbr[0] ? nbr[1] : nbr[0]`): +1; `lostal[center] := none` (actions.cpp:379-390) | — |
| **R9** | reverse of R8: =Si-OH + Al-OH + Al-H2O → =Si-OH-Al + Al-OH + H2O | 403 → 410 | Si (nbr[2]): −1; ONE Al loses −1 — chosen by coin; the **other** Al is recorded as `lostal[center]` | one `ran2()`: r < 0.5 → nbr[0] loses (−1) and lostal := nbr[1]; else nbr[1] loses and lostal := nbr[0] (actions.cpp:392-410) |
| **R10** | Si-OH-Al hydrolysis (one Al vacant): =Si-OH-Al + H2O → =Si-OH + Al-H2O | 406 → 407 | Si (nbr[2]): +1; the **occupied** Al: `(al1.state == 100) ? al2 : al1` gets +1 (actions.cpp:412-423) | — |
| **R11** | reverse of R10 | 407 → 406 | Si (nbr[2]): −1; the occupied Al (same ==100 test): −1 | — |
| **R12** | Al-OH-Al hydrolysis at a 4xx site (Si vacant): HO<Al2 + H2O → HO-Al + H2O-Al | 404 → 405 | both Al (nbr[0], nbr[1]): each +1 | — |
| **R13** | reverse of R12 | 405 → 404 | both Al: each −1. **NEVER FIRES from shipped data**: data.rxn:48 mis-types the product as `40100`, so `reactions[13].reactant = 40100` matches no site (REFORM_PLAN:135-141; confirmed 0 events over 20,000 oracle steps). Mechanism is correct and ported (mechanisms.rs:281-292). | — |
| **R14** | Al-OH-Al hydrolysis, normal 5xx site: Al-OH-Al + H2O → Al-OH + Al-H2O | 501 → 502 | both Al (nbr[0], nbr[1]): each +1 | — |
| **R15** | reverse of R14 | 502 → 501 | both Al: each −1 | — |

**Proton-coin semantics (R4/R9), exact:** a single uniform draw `r = ran2()` from the
same stream that generated the step's `dt` and event-selection `eps`, drawn *after*
them, mid-mutation. `r < 0.5` selects the first Al slot (strictly less-than). R4:
selected Al gains the proton, `lostal := selected`. R9: selected Al *loses* (−1),
`lostal := the other one`. Draw-position in the RNG stream is the most parity-fragile
point in the model (mechanisms.rs:18-27) — irrelevant for a new engine's physics, but
it means R4/R9 consume one extra RNG draw per firing.

Sanity: every reverse is an exact structural inverse of its forward given the `lostal`
record (round-trip pinned by tests, mechanisms.rs:449-493).

---

## 4. Adsorption / desorption (reactions 16–23)

Rate-table reactants/products from data.rxn:58-86 and rxnlist.cpp:79-107; mechanisms
from actions.cpp:122-278 (mechanisms.rs:314-430).

| rxn | name | center transition | active? |
|---|---|---|---|
| 16 | adsorb Al onto empty Al site | 100 → 107 | ≥1 occupied non-EDGE neighbor |
| 17 | adsorb Al onto Si site (cross) | 200 → 299 | **never** |
| 18 | adsorb Si onto Al site (cross) | 100 → 199 | **never** |
| 19 | adsorb Si onto empty Si site | 200 → 205 | ≥1 occupied non-EDGE neighbor |
| 20 | desorb Al | 107 → 100 | always |
| 21 | desorb Si from Al site (cross) | reactant 199 | never proposed (no 199 exists); would run DesorbSi |
| 22 | desorb Si | 205 → 200 | always |
| 23 | desorb Al from Si site (cross) | reactant 299 | never proposed; would run DesorbAl |

Guards: adsorption requires the site empty (reactant match 100/200) AND `is_active`'s
"at least one occupied, non-EDGE neighbor" (can't nucleate on nothing,
envrn.cpp:295-307). Desorption requires only reactant match (107 or 205 exactly —
partially bonded cations 101–106/201–204 can NOT desorb).

### Multi-site side effects — the oxygen-shell rewrites (critical)

Adsorbing/desorbing a cation rewrites **every** neighboring oxygen by a fixed
transition map. An oxygen in a state not in the map is a **fatal error** for
adsorption (C++ prints "invalid state" and the event fails, actions.cpp:153-156) but
is silently **left unchanged** for desorption (`default: break`, actions.cpp:235-236,
271-272).

**AdsorbAl (rxn 16): center 100 → 107; each of the 6 oxygen neighbors:**
(actions.cpp:126-158, mechanisms.rs:320-354)

| O before | O after | note |
|---|---|---|
| 500 | 503 | empty Al-OH-Al gains one Al |
| 503 | 502 | Al-OH + new Al → Al-OH H2O-Al |
| 400 | 409 | empty 4xx gains one Al → HO-Al |
| 409 | 405 | HO-Al + new Al → HO-Al H2O-Al |
| 408 | 407 | Si-OH + new Al → Si-OH HO-Al |
| 407 | 403 | Si-OH HO-Al + new Al → fully hydrolyzed 3-cation config |
| 406 | 410 | Si-OH-Al + new Al → Si-OH-Al HO-Al; **also sets `lostal[O] := adsorbing site`** (the new Al holds the proton) |
| anything else | FATAL "invalid state in adsorbAl" |

Note: the adsorbing Al arrives fully protonated (107 = 6 ligands) and every O it
touches lands in a *hydrolyzed* (not bridging-bonded) configuration — adsorption
never directly forms 501/401 bridges; those form via reverse hydrolysis afterwards.

**AdsorbSi (rxn 19): center 200 → 205; each of the 4 oxygen neighbors:**
(actions.cpp:166-202, mechanisms.rs:358-382)

| O before | O after |
|---|---|
| 300 | 303 |
| 303 | 302 |
| 400 | 408 |
| 409 | 407 |
| 404 | 402 |
| 405 | 403 |
| anything else | FATAL "invalid state in adsorbSi" |

**DesorbAl (rxn 20): center 107 → 100; each of the 6 oxygen neighbors (exact inverse
of AdsorbAl):** (actions.cpp:205-242, mechanisms.rs:387-409)

| O before | O after |
|---|---|
| 502 | 503 |
| 503 | 500 |
| 403 | 407 |
| 405 | 409 |
| 407 | 408 |
| 409 | 400 |
| 410 | 406, and **`lostal[O] := none`** |
| anything else | unchanged (silent) |

**DesorbSi (rxn 22): center 205 → 200; each of the 4 oxygen neighbors:**
(actions.cpp:245-278, mechanisms.rs:412-430)

| O before | O after |
|---|---|
| 303 | 300 |
| 302 | 303 |
| 402 | 404 |
| 403 | 405 |
| 407 | 409 |
| 408 | 400 |
| anything else | unchanged (silent) |

Implicit invariant: desorption is only legal at 107/205, i.e. all the cation's bonds
are already hydrolyzed, which is why the O maps only need hydrolyzed input states.
The desorb maps are *not* total inverses at the map level (e.g. DesorbAl leaves an
unexpected O silently alone) — but on states reachable in the dynamics they invert
exactly (round-trip test mechanisms.rs:496-517).

Cross variants (17/18/21/23): mechanistically `AdsorbAl`/`AdsorbSi` on a wrong-class
site just set state to 299/199 with **no neighbor updates** (actions.cpp:159-161,
198-200); `DesorbAl`/`DesorbSi` on a wrong-class site set 200/100. All dead code —
`is_active` never allows 17/18/21/23 (environment.rs:115-118), and no 199/299 ever
exists to match 21/23's reactant.

---

## 5. Rates (from data.rxn, formulas from rxnlist.cpp)

### Global conditions (data.rxn:1-3)

| parameter | value |
|---|---|
| Temperature T | 8000.0 K (a relative-rate knob, not a physical claim — reactions.rs:76-77) |
| Δμ_Si | −1.0 kcal/mol |
| Δμ_Al | −1.0 kcal/mol |
| R (gas constant) | 1.987e-3 kcal/mol/K (`#define R` — reactions.rs:46); RT = R·T ≈ 15.896 kcal/mol |

All rate math is `float` (f32), computed as `a * exp(de/rt)` etc. in C float
precision (rxnlist.cpp:39, 72-73; f32 kept in the port for parity, reactions.rs:63-68).

### Formulas (rxnlist.cpp:60-107)

- **Forward hydrolysis** (even i): `rate = k+` (the raw table value).
- **Reverse hydrolysis** (odd i = even+1): `rate = k+ · exp(dE / RT)` — computed from
  the SAME (k+, dE) pairs read for the forward reaction (one table serves the pair,
  rxnlist.cpp:70-74). With dE > 0, reverse is *faster* than forward.
- **Adsorption** (i = 16..19): single rate, `rate = a · exp(dE/RT) · exp(Δμ/RT)`,
  where `a` = pre-exponential ("ν·exp(μ_solid/RT)", rxnlist.cpp:77) and — **WART** —
  the Δμ is selected by `i < 18 ? Δμ_Si : Δμ_Al` (rxnlist.cpp:89-92), which routes
  reaction 16 (**Al** adsorption) through **Δμ_Si**. Inert today because
  Δμ_Si == Δμ_Al == −1.0, but a deck with distinct potentials hits the swap
  (ported verbatim: kmc-rs/crates/kmc-io/src/rxn.rs:91-94; flagged in RUST_TOUR.md:343).
- **Desorption** (i = 20..23): `rate = a · exp(−dE/RT)` per bucket (rxnlist.cpp:103-106).
- **Diffusion** (24–27): tables copied from 20–23 with a stale-`numRates` buffer
  over-read (rxnlist.cpp:110-116); dead code, never active; not ported
  (reactions.rs:10-17).

### Per-reaction parameters (data.rxn:7-86 — every variant table is FLAT: all buckets
within a reaction carry the identical (k, dE) pair, so one value per reaction suffices)

| rxn pair | reactant → product | n buckets | k+ | dE (kcal/mol) | forward rate | reverse rate = k+·e^(dE/RT) |
|---|---|---|---|---|---|---|
| R0/R1 | 301 ↔ 302 (Si-O-Si) | 15 | 1 | 2.6 | 1.0 | ≈ 1.1777 |
| R2/R3 | 401 ↔ 402 (Si-O-Al2, Si side) | 40 | 1 | 2.6 | 1.0 | ≈ 1.1777 |
| R4/R5 | 401 ↔ 410 (Si-O-Al2, Al side) | 40 | 100 | 4.8 | 100.0 | ≈ 135.25 |
| R6/R7 | 402 ↔ 403 (Si-OH Al-OH-Al) | 40 | 100 | 4.8 | 100.0 | ≈ 135.25 |
| R8/R9 | 410 ↔ 403 (Si-OH-Al Al-OH) | 40 | 1 | 2.6 | 1.0 | ≈ 1.1777 |
| R10/R11 | 406 ↔ 407 (Si-OH-Al) | 20 | 1 | 2.6 | 1.0 | ≈ 1.1777 |
| R12/R13 | 404 ↔ **40100** (Al-OH-Al on 4xx; product token is the typo — should be 405) | 20 | 100 | 4.8 | 100.0 | ≈ 135.25 (assigned to reactant 40100 → never used) |
| R14/R15 | 501 ↔ 502 (Al-OH-Al normal) | 20 | 100 | 4.8 | 100.0 | ≈ 135.25 |

(Reverse-rate numbers use RT = 15.896; exact stored values are f32 of `k·exp(dE/RT)`.)

| rxn | reaction | a | dE | Δμ used | rate formula | value |
|---|---|---|---|---|---|---|
| 16 | adsorb Al, 100→107 | 100.0 | −14.5 | Δμ_Si (**swap wart**) = −1.0 | a·e^(dE/RT)·e^(Δμ/RT) | ≈ 100·e^(−14.5/15.896)·e^(−1/15.896) ≈ 37.7 |
| 17 | adsorb Al→Si site, 200→299 | 0.0 | −14.5 | Δμ_Si | same | 0 (and never active) |
| 18 | adsorb Si→Al site, 100→199 | 0.0 | −6.4 | Δμ_Al | same | 0 (and never active) |
| 19 | adsorb Si, 200→205 | 100.0 | −6.4 | Δμ_Al = −1.0 | same | ≈ 100·e^(−6.4/15.896)·e^(−1/15.896) ≈ 62.8 |

| rxn | reaction | buckets (a, dE) | rate per bucket a·e^(−dE/RT) |
|---|---|---|---|
| 20 | desorb Al, 107→100 | 4: (1998, 0), (1998, 12), (1998, 24), (1998, 36) | 1998.0, ≈939.2, ≈441.5, ≈207.5 |
| 21 | desorb Si from Al site (199) | 4: same as rxn 20 | same (dead) |
| 22 | desorb Si, 205→200 | 5: (1998, 0), (1998, 6), (1998, 12), (1998, 18), (1998, 24) | 1998.0, ≈1369.8, ≈939.2, ≈643.9, ≈441.5 |
| 23 | desorb Al from Si site (299) | 5: same as rxn 22 | same (dead) |

**Desorption is the ONE place the environment bucket is not flat**: the 4/5 desorption
buckets carry distinct dE (0/12/24/36 for Al via Check100's `(x+y)/2`; 0/6/12/18/24
for Si via Check200's `x + y`) — so desorption rate genuinely depends on how many
intact bridges the environment counts. All hydrolysis and adsorption tables are flat.

**Forward/reverse partners:** (R0,R1), (R2,R3), (R4,R5), (R6,R7), (R8,R9), (R10,R11),
(R12,R13 — broken by the 40100 typo), (R14,R15). Adsorption/desorption pairs by
species: (16, 20) for Al, (19, 22) for Si. Detailed balance for the 404↔405 channel
is broken in the shipped deck (R13 inert), making R12 effectively irreversible
(REFORM_PLAN:142-149) — except that AdsorbSi/DesorbSi also map 404↔402/405↔403, so
405 is not a total sink.

### Event selection & clock (context, actions.cpp:10-45)

Standard Gillespie direct method per step: `ratesum = Σ rates`;
`dt = −ln(ran2())/ratesum`; second draw `eps` selects the event by cumulative
normalized rate (`eps <= partsum`, first crossing). Fatal if ratesum == 0. Time and
rates are f32 under `-ffast-math` (REFORM_PLAN R6).

---

## 6. Special cases — what breaks a simple "guarded local rewrite" model

Loudly flagged, in rough order of pain for a declarative engine:

1. **The `lostal` memory (R4/R5, R8/R9, AdsorbAl/DesorbAl).** A 410 site carries
   hidden state: WHICH Al holds the proton. R5 and R8 rewrite *the neighbor selected
   by that memory*, not by a state test. A pure state-rewrite engine must either
   (a) support per-site auxiliary tags written/read by rules, or (b) split state 410
   into two states (proton-on-nbr[0] vs proton-on-nbr[1]) — option (b) makes lostal
   pure state but bakes a slot-ordering convention into the state space. Note the
   information is *almost* recoverable from Al states (the lostal Al is the more
   protonated one) but not reliably: both Al evolve independently while the center
   sits in 410.

2. **Stochastic rewrite outcomes (R4, R9).** One rule, two outcomes chosen by a fair
   coin at application time (which Al gains/loses). The engine needs
   "choose-one-of-N-symmetric-targets uniformly" in its rule language, plus the
   lostal write tied to the choice.

3. **Adsorption/desorption are 1+k-site rewrites with per-neighbor case maps.** One
   event rewrites the center AND all 6 (Al) / 4 (Si) neighbors, each through its own
   state-transition map (§4), including a nested `lostal` write on the 406→410 arm of
   AdsorbAl. Also asymmetric error semantics: unknown neighbor state is fatal on
   adsorb, silently ignored on desorb.

4. **`is_active` reaches TWO hops out** (forward hydrolysis: any 2nd-neighbor in
   {303,404,405,406,408,409}). Guards are not nearest-neighbor-only. Worse, the
   *legacy* rule adds the phantom "OR any neighbor has all 6 neighbor slots filled" —
   a static topological predicate born from an out-of-bounds read (§2). A new engine
   should implement the corrected 2nd-neighbor rule and treat legacy-phantom as an
   optional compatibility quirk, not a modeling primitive.

5. **`check_env` reaches through specific neighbor slots and the `pair` partner.**
   Check300 reads each Si neighbor's `nbr[0]` (a hard-coded "the 4xx oxygen is slot 0"
   convention); Check400 reads the Si's `nbr[1]`/`nbr[2]` and the pair partner;
   Check500 reads the pair partner. So rates (once non-flat) depend on: 1st neighbors,
   selected 2nd neighbors via fixed slot positions, and a non-adjacency "pair" relation
   computed at setup. The engine needs (a) named neighbor roles (slot semantics differ
   by class: 4xx = [Al, Al, Si]; 3xx = [Si, Si]; 5xx = [Al, Al]) and (b) a static
   binary site relation (`pair`) usable in rate predicates.

6. **Multi-neighbor simultaneous conditions.** Check400's bucket is a joint function
   of BOTH Al states, the center, the pair partner, and two of the Si's other oxygens
   (a ~7-site predicate). Check300 likewise joins 6 sites. As long as tables stay flat
   this is only a legality/abort gate; with real rates it is full n-body rate
   modulation.

7. **The 404/405 class-4→Check500 reroute** (envrn.cpp:27-32): the bucket function is
   chosen by *state*, not class — 404/405 are "chemically 500s living at a 400 site".
   A declarative deck keyed purely on site class will get this wrong.

8. **Reaction asymmetries in the shipped deck**: R13 dead via the `40100` typo
   (detailed balance broken for 404↔405); the adsorption Δμ swap (`i<18` → Al uses
   Δμ_Si); Check100/Check200 formulas contradict their documentation. All inert with
   the shipped flat/equal parameters; all live the moment you author real numbers.
   Decide explicitly which behavior the new deck encodes.

9. **Slot-order sensitivity everywhere.** Mechanisms address neighbors positionally
   (nbr[0], nbr[1], nbr[2]) with class-specific meanings; the R4/R9 coin maps r<0.5 to
   slot 0. Any engine that canonicalizes/permutes neighbor lists changes the
   trajectory (though not the statistics, for the coin).

10. **Non-issues, for the record:** the BFS `RemoveUnattachedClusters` is NEVER called
    in live dynamics (only from the dead diffusion arm, actions.cpp:114-117) — despite
    top-level docs suggesting otherwise, no cluster-removal / action-at-a-distance
    step exists in the working model. Diffusion (24–27) and all cross-cation reactions
    (17/18/21/23) are permanently disabled by `is_active`. EDGE sites are inert but
    must be *visible* to guards (occupied-neighbor tests exclude them explicitly;
    check_env treats them as fatal).
