//! Building the lattice: tile the unit cell, resolve neighbors, then (M3)
//! the deterministic structural setup. Port of `lattice.cpp`.
//!
//! Everything in this module is **deterministic and RNG-free** — the reason
//! the milestone ladder front-loads it: its output can be diffed bitwise
//! against the C++ golden capture (spec §C1), which is the strongest
//! validation checkpoint the whole port has.
//!
//! Order matters and is fixed by `mckaol.cpp`:
//! [`create_lattice`] → `find_pairs` → `populate_solid` →
//! `terminate_surface` → `terminate_lattice` (the last four land at M3).

use kmc_engine::{Site, SiteGraph, SiteId};

use crate::cell::UnitCell;
use crate::state::State;

/// Lattice dimensions from `data.lattice` (read by `kmc-io`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LatticeParams {
    /// Number of unit cells along a.
    pub a_cells: i32,
    /// Number of unit cells along b.
    pub b_cells: i32,
    /// Which face is the exposed surface: 0 = ac plane (b is the open,
    /// surface-normal direction; a wraps periodically), 1 = bc plane (a
    /// open, b periodic). Kept as `i32` because the C++ tests its
    /// truthiness (`SurfacePlane` as a flag) rather than comparing to 1,
    /// and the port follows suit.
    pub surface_plane: i32,
}

/// Where a site sits in the tiling: which unit cell, which position.
///
/// In the C++ these are the `a`, `b`, `n` members of `LatticeSite`. They
/// are kaolinite bookkeeping (output writers and boundary logic use them),
/// not generic-KMC data, so they live here in a parallel array rather than
/// on the engine's `Site` (design doc §3, option (a)).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CellCoord {
    /// Unit-cell index along a.
    pub a: i32,
    /// Unit-cell index along b.
    pub b: i32,
    /// Position within the unit cell (0..npos).
    pub n: i32,
}

/// The built lattice: the engine-visible graph plus kaolinite's per-site
/// side data.
///
/// \[IDIOM\] **Parallel arrays instead of a fat site struct.** The C++
/// `LatticeSite` carries everything every subsystem might want (`state`,
/// `nbr`, `pair`, `lostal`, `a/b/n`, a dead BFS `color`). We split it:
/// the engine's `Site` holds what generic KMC needs, and each
/// kaolinite-only attribute is a `Vec` indexed by the same `SiteId`.
/// What this buys besides layering: the borrow checker can prove that
/// mutating `pair` doesn't alias `graph`, because they are separate
/// values — one `&mut` each. A fat struct would make every mutation a
/// borrow of *everything*. (This is also the seam the parallelism memo
/// §6 wants: a future domain decomposition reasons about `pair` crossing
/// strip boundaries; a flat `Vec<Option<SiteId>>` is easy to reason about.)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Structure {
    /// States + topology (the engine's view).
    pub graph: SiteGraph<State>,
    /// Tiling coordinates per site (C++ `a`, `b`, `n` members).
    pub coord: Vec<CellCoord>,
    /// Double-bridge partner per site (C++ `pair`, -1 → `None`): for
    /// 400/500-class bridging oxygens, the other oxygen bridging the same
    /// two Al. Filled by `find_pairs` (M3).
    pub pair: Vec<Option<SiteId>>,
    /// For 400s that lost one Al: which Al it was (C++ `lostal`, -1 →
    /// `None`). Mutated only by the reactions (M5+); built as all-`None`.
    pub lostal: Vec<Option<SiteId>>,
    /// The dimensions this structure was built with.
    pub params: LatticeParams,
}

/// Tile the unit cell into the flat site array and resolve every neighbor.
/// Port of `Lattice::CreateLattice` (lattice.cpp).
///
/// Linear index scheme — identical to the C++ so states can be compared
/// index-for-index against the golden capture:
/// `i = a * (b_cells * npos) + b * npos + n`.
pub fn create_lattice(uc: &UnitCell, params: LatticeParams) -> Structure {
    let npos = uc.npos();
    let num_sites = params.a_cells as usize * params.b_cells as usize * npos;

    let mut sites = Vec::with_capacity(num_sites);
    let mut coord = Vec::with_capacity(num_sites);
    // Same triple loop, same order, as the C++ — push order IS the index
    // scheme, so getting this loop nest wrong would scramble every site id.
    for a in 0..params.a_cells {
        for b in 0..params.b_cells {
            for n in 0..npos {
                let mut nbr = [None; 6];
                for (j, slot) in nbr.iter_mut().enumerate() {
                    *slot = resolve_neighbor(uc, params, a, b, n, j);
                }
                sites.push(Site {
                    state: uc.sites[n].state,
                    nbr,
                });
                coord.push(CellCoord { a, b, n: n as i32 });
            }
        }
    }

    Structure {
        graph: SiteGraph { sites },
        coord,
        pair: vec![None; num_sites],
        lostal: vec![None; num_sites],
        params,
    }
}

/// Resolve one neighbor-template slot into a concrete site index, applying
/// the slab boundary conditions. Port of `Lattice::GetNeighbor`
/// (lattice.cpp) — the periodic/open logic below is transliterated, `==`
/// comparisons and all.
///
/// The slab: the open (surface-normal) direction returns `None` past its
/// ends; the in-plane direction wraps periodically. Which is which is
/// chosen by `surface_plane` (see [`LatticeParams`]).
fn resolve_neighbor(
    uc: &UnitCell,
    p: LatticeParams,
    a: i32,
    b: i32,
    n: usize,
    j: usize,
) -> Option<SiteId> {
    let npos = uc.npos() as i32;
    let t = uc.sites[n].nbr[j];
    let mut na = a + t.a;
    let mut nb = b + t.b;
    if t.n < 0 {
        // Template says "no j-th neighbor". (C++ returns the negative n
        // itself; any negative collapses to None here.)
        return None;
    }
    // Open-boundary check on whichever axis is the surface normal. The C++
    // tests `a == Num_aCells` (not >=) because template offsets are ±1 at
    // most; ported as-is.
    if ((na == p.a_cells || na < 0) && p.surface_plane != 0)
        || ((nb == p.b_cells || nb < 0) && p.surface_plane == 0)
    {
        return None;
    }
    // Periodic wrap on the in-plane axis (and vacuously on the open axis,
    // which the check above already fenced off).
    if na >= p.a_cells {
        na = 0;
    } else if na < 0 {
        na = p.a_cells - 1;
    }
    if nb >= p.b_cells {
        nb = 0;
    } else if nb < 0 {
        nb = p.b_cells - 1;
    }
    Some((na * p.b_cells * npos + nb * npos + t.n) as SiteId)
}

/// Link each 400/500-class bridging oxygen to its double-bridge partner —
/// the other oxygen bridging the *same* two Al. Port of `Lattice::FindPairs`
/// (lattice.cpp).
///
/// Algorithm, as implemented: for each unpaired site with state ≥ 400, take
/// its first two neighbors (`al1`, `al2`); scan `al1`'s neighbors for an
/// oxygen whose own first two neighbors are also `{al1, al2}`; link the two.
/// Hitting a missing neighbor (`None`) while scanning **gives up on this
/// oxygen** — the C++ sets `found` on a boundary and moves on, leaving
/// `pair = None`.
///
/// Ported quirk (kept, it's load-bearing behavior): the scan does not
/// exclude the oxygen itself. `o1` is one of `al1`'s neighbors, and `o1`'s
/// first two neighbors are trivially `{al1, al2}` — so if the template
/// order presents `o1` before its true partner, a site pairs with
/// *itself*. Whether that ever happens depends on the unit cell's neighbor
/// ordering; the C++ semantics, not our opinion of them, are the contract.
///
/// \[IDIOM\] The **read-then-write borrow dance**. Inside the loop we read
/// from `s.graph` (immutable borrows) and write to `s.pair` (mutable) —
/// fine, they are different fields, and the borrow checker tracks field
/// borrows separately. What it would *reject* is holding `&s.graph.sites[o1]`
/// while mutating `s.graph.sites[o2]` — the aliasing the C++ does freely
/// (and gets away with here, but exactly such patterns hide the
/// iterator-invalidation class of C++ bugs). The idiomatic fix used
/// throughout this module: copy the small values you need into locals
/// (`al1`, `al2` are just `usize`s), drop the borrow, then mutate by index.
/// When the borrow checker fights you, it is usually pointing at a real
/// aliasing question — answer it, don't `RefCell` it away (design doc §8).
pub fn find_pairs(s: &mut Structure) {
    for o1 in 0..s.graph.len() {
        if s.graph.sites[o1].state.0 < 400 || s.pair[o1].is_some() {
            continue;
        }
        // C++: al1/al2 = nbr[0]/nbr[1]; if either is -1, skip this oxygen.
        let (Some(al1), Some(al2)) = (s.graph.sites[o1].nbr[0], s.graph.sites[o1].nbr[1]) else {
            continue;
        };
        // [IDIOM] `let ... else`: destructure or bail. The happy path stays
        // unindented; the "pattern didn't match" arm must diverge (here:
        // `continue`), so after this line al1/al2 are plain usize — no
        // Option left to unwrap, no -1 to forget to test.
        let mut found = false;
        let mut j = 0;
        while j < 6 && !found {
            match s.graph.sites[al1].nbr[j] {
                // No neighbor: C++ comments "must be edge" and gives up.
                None => found = true,
                Some(o2) => {
                    let al3 = s.graph.sites[o2].nbr[0];
                    let al4 = s.graph.sites[o2].nbr[1];
                    if (al3 == Some(al1) || al4 == Some(al1))
                        && (al3 == Some(al2) || al4 == Some(al2))
                    {
                        found = true;
                        s.pair[o1] = Some(o2);
                        s.pair[o2] = Some(o1);
                    }
                }
            }
            j += 1;
        }
    }
}

/// What [`populate_solid`] decided, for the caller to report (the C++
/// prints this from inside the function; we return it so the library stays
/// silent and the CLI owns stdout).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PopulateReport {
    /// "supersaturated" / "undersaturated" / "near equilibrium".
    pub condition: &'static str,
    /// How many unit-cell layers were filled (the C++ `top`).
    pub filled_cells: i32,
}

/// Fill the initial solid slab by promoting empty sites to first-occupied
/// (100→101, 200→201, ...). Port of `Lattice::PopulateSolid` (lattice.cpp).
///
/// The fill *fraction* comes from the chemical-potential sum (spec A4.2):
/// supersaturated (> 0.5) fills 0.3 of the depth, undersaturated (< -0.5)
/// fills 0.7, otherwise 0.5. The golden inputs (Δμ = -1 + -1 = -2) are
/// undersaturated → 0.7 → with 3 b-cells, `(3 × 0.7f) as int`.
///
/// **Numeric care (M3 gate):** that layer count is `int = int × float` in
/// C++ — the product is computed in f32 and truncated toward zero, and
/// whether `cells × 0.7f` lands just above or just below the integer is at
/// the mercy of f32 rounding (3 × 0.7f32 = 2.0999999… → 2, while
/// 20 × 0.7f32 rounds *up* to exactly 14.0 → 14). Any "cleaner" arithmetic
/// — f64, or `(cells * 7) / 10` — could fill a different slab on some
/// lattice size. Rust's `as i32` truncates toward zero exactly like the
/// C++ float→int conversion, so the transliteration below is bit-faithful.
/// This is the first of the port's "f32 or bust" spots — see
/// docs/RUST_TOUR.md §"The bitwise gate".
///
/// One guard (spec A4.2): an empty Si-O-Si oxygen (300) is only promoted if
/// its full connectivity chain exists — both Si neighbors present, and each
/// Si's *first* neighbor present.
pub fn populate_solid(s: &mut Structure, dm_si: f32, dm_al: f32) -> PopulateReport {
    let p = s.params;
    // graph.len() = a_cells * b_cells * npos by construction, so npos is
    // recoverable — the C++ asks the unit cell, we ask the arithmetic.
    let npos = s.graph.len() as i32 / (p.a_cells * p.b_cells);

    // C++ compares the f32 sum against double literals (usual arithmetic
    // conversions promote the sum); f64 comparison ported as-is.
    let sum = (dm_si + dm_al) as f64;
    let (frac, condition): (f32, &'static str) = if sum > 0.5 {
        (0.3, "supersaturated")
    } else if sum < -0.5 {
        (0.7, "undersaturated")
    } else {
        (0.5, "near equilibrium")
    };

    // int = int * float, truncating. See the doc comment.
    let top = if p.surface_plane != 0 {
        (p.a_cells as f32 * frac) as i32
    } else {
        (p.b_cells as f32 * frac) as i32
    };
    let len = if p.surface_plane != 0 {
        p.b_cells
    } else {
        p.a_cells
    };

    for z in 0..top {
        for x in 0..len {
            // if bc surface, x runs over b; else x runs over a (C++ comment).
            let mut i = if p.surface_plane != 0 {
                (z * p.b_cells * npos + x * npos) as usize
            } else {
                (x * p.b_cells * npos + z * npos) as usize
            };
            for _ in 0..npos {
                let mut can_increment = true;
                if s.graph.sites[i].state == State(300) {
                    // 301 (Si-O-Si) needs its complete connectivity chain.
                    match (s.graph.sites[i].nbr[0], s.graph.sites[i].nbr[1]) {
                        (Some(si1), Some(si2)) => {
                            let sio1 = s.graph.sites[si1].nbr[0];
                            let sio2 = s.graph.sites[si2].nbr[0];
                            if sio1.is_none() || sio2.is_none() {
                                can_increment = false;
                            }
                        }
                        _ => can_increment = false,
                    }
                }
                if can_increment {
                    s.graph.sites[i].state.0 += 1;
                }
                i += 1;
            }
        }
    }

    PopulateReport {
        condition,
        filled_cells: top,
    }
}

/// Cap dangling bonds where the solid meets water: demote oxygens that lost
/// a cation, hydroxylate oxygens next to minimally-coordinated cations, and
/// bump cation states for each terminal OH. Port of
/// `Lattice::TerminateSurface` (lattice.cpp) — three passes, transliterated.
pub fn terminate_surface(s: &mut Structure) {
    let n_sites = s.graph.len();

    // Pass 1 — "dangle" the O's which are missing cations: every occupied
    // oxygen adjacent to an empty cation site steps down a fixed state map.
    for i in 0..n_sites {
        let st = s.graph.sites[i].state;
        if st.0 % 100 == 0 || st.0 / 100 < 3 || st == State::EDGE {
            continue;
        }
        for j in 0..6 {
            // [IDIOM] `if let Some(i2) = ...` is the C++ `if (i2 >= 0)`
            // guard, except unskippable: there is no way to reach the body
            // with a missing neighbor. First of ~a dozen in this module;
            // every -1 check the C++ author remembered to write becomes one
            // the compiler would have demanded anyway.
            if let Some(i2) = s.graph.sites[i].nbr[j] {
                let t = s.graph.sites[i2].state.0 % 100;
                if t == 0 {
                    // Ported AS-WRITTEN, quirk included: the C++ reads
                    //   (type = sites[i2].state % 100) == 0
                    // and then branches on `type == 2` — but `type` is 0 in
                    // this block by construction, so the `? 404 :`/`? 408 :`
                    // arms are unreachable and 401 ALWAYS demotes to 406,
                    // 406 always to 409. The spec's Part A (A4.3) describes
                    // the seemingly *intended* Si/Al distinction; the code
                    // never had it. Not in the spec's Part B wart list —
                    // discovered during this port, reported in the task
                    // Result. We reproduce the code, not the intent.
                    let cur = s.graph.sites[i].state.0;
                    let new = match cur {
                        401 => {
                            if t == 2 { 404 } else { 406 } // t==0: always 406
                        }
                        404 => 409,
                        406 => {
                            if t == 2 { 408 } else { 409 } // t==0: always 409
                        }
                        409 | 408 => 400,
                        301 => 303,
                        303 => 300,
                        501 => 503,
                        503 => 500,
                        other => other, // C++ default: leave unchanged
                    };
                    s.graph.sites[i].state = State(new);
                }
            }
        }
    }

    // Pass 2 — turn on dangling oxygens to terminate bonds: for occupied
    // cations at minimal coordination (occupancy == 1), promote adjacent
    // not-minimally-occupied oxygens to their OH-terminated forms.
    for i in 0..n_sites {
        let st = s.graph.sites[i].state;
        let cation_class = st.0 / 100;
        if cation_class > 2 || st.0 % 100 != 1 {
            continue;
        }
        for j in 0..6 {
            // [IDIOM] A let-chain: `if let Some(x) = ... && condition` —
            // pattern match and guard in one breath (stabilized in edition
            // 2024). Ports the C++'s `if (i2 >= 0 && sites[i2].state % 100
            // != 1)` as a single condition instead of nested ifs.
            if let Some(i2) = s.graph.sites[i].nbr[j]
                && s.graph.sites[i2].state.0 % 100 != 1
            {
                let cur = s.graph.sites[i2].state.0;
                let new = match cur {
                    300 => 303,
                    400 => {
                        if cation_class == 2 {
                            408
                        } else {
                            409
                        }
                    }
                    408 => 406,
                    409 => {
                        if cation_class == 2 {
                            406
                        } else {
                            404
                        }
                    }
                    404 => 401,
                    406 => 401,
                    500 => 503,
                    503 => 501,
                    other => other,
                };
                s.graph.sites[i2].state = State(new);
            }
        }
    }

    // Pass 3 — update Si and Al states to count their terminal OH's:
    // +1 per adjacent 503/409 for Al, +1 per adjacent 303/408 for Si.
    for i in 0..n_sites {
        let st = s.graph.sites[i].state;
        if st.0 % 100 == 0 || st.0 / 100 > 2 || st == State::EDGE {
            continue;
        }
        let is_al = st.0 / 100 == 1;
        for j in 0..6 {
            if let Some(i2) = s.graph.sites[i].nbr[j] {
                let nb = s.graph.sites[i2].state.0;
                let bump = if is_al {
                    nb == 503 || nb == 409
                } else {
                    nb == 303 || nb == 408
                };
                if bump {
                    s.graph.sites[i].state.0 += 1;
                }
            }
        }
    }
}

/// Freeze the slab's cut faces: mark first/last rows along the open axis as
/// [`State::EDGE`] where a site is missing neighbors, then revert any
/// Si-O-Si bridge (301) whose connectivity chain touches an EDGE back to
/// empty (300). Port of `Lattice::TerminateLattice` (lattice.cpp).
///
/// Must run **after** [`terminate_surface`] — `mckaol.cpp` fixes the order,
/// and the state-vs-expected-neighbor-count test below depends on the
/// states `terminate_surface` left behind.
pub fn terminate_lattice(s: &mut Structure) {
    let p = s.params;
    let npos = s.graph.len() as i32 / (p.a_cells * p.b_cells);
    let top = if p.surface_plane != 0 {
        p.a_cells - 1
    } else {
        p.b_cells - 1
    };
    let len = if p.surface_plane != 0 {
        p.b_cells
    } else {
        p.a_cells
    };

    // Pass 1: boundary rows → EDGE wherever the actual neighbor count
    // differs from what the site's class requires.
    for x in 0..len {
        let mut i = if p.surface_plane != 0 {
            (top * p.b_cells * npos + x * npos) as usize
        } else {
            (x * p.b_cells * npos + top * npos) as usize
        };
        for _ in 0..npos {
            if s.graph.sites[i].state.expected_neighbor_count() != s.graph.count_nbrs(i) {
                s.graph.sites[i].state = State::EDGE;
            }
            i += 1;
        }
        let mut i = if p.surface_plane != 0 {
            (x * npos) as usize
        } else {
            (x * p.b_cells * npos) as usize
        };
        for _ in 0..npos {
            if s.graph.sites[i].state.expected_neighbor_count() != s.graph.count_nbrs(i) {
                s.graph.sites[i].state = State::EDGE;
            }
            i += 1;
        }
    }

    // Pass 2: any 301 whose Si→O connectivity chain now touches an EDGE
    // reverts to empty (300).
    for i in 0..s.graph.len() {
        if s.graph.sites[i].state != State(301) {
            continue;
        }
        if let (Some(si1), Some(si2)) = (s.graph.sites[i].nbr[0], s.graph.sites[i].nbr[1]) {
            let sio1 = s.graph.sites[si1].nbr[0];
            let sio2 = s.graph.sites[si2].nbr[0];
            // [IDIOM] `Option::is_some_and` — "present AND passes this
            // test" in one breath, porting `(sio1 >= 0 && sites[sio1].state
            // == EDGE)` without a nested if. Combinators like this are
            // just `match` in a dress; prefer them exactly as long as the
            // sentence stays readable.
            let touches_edge = s.graph.sites[si1].state == State::EDGE
                || s.graph.sites[si2].state == State::EDGE
                || sio1.is_some_and(|k| s.graph.sites[k].state == State::EDGE)
                || sio2.is_some_and(|k| s.graph.sites[k].state == State::EDGE);
            if touches_edge {
                s.graph.sites[i].state = State(300);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::{CellSite, NeighborTemplate};

    /// A minimal 2-position cell for BC tests: position 0 links to
    /// position 1 one cell over in +a, and to position 1 one cell over in
    /// +b; position 1 links back. Degenerate chemistry (all Al) — these
    /// tests are about *topology* only.
    fn tiny_cell() -> UnitCell {
        let t = |n, a, b| NeighborTemplate { n, a, b, c: 0 };
        let none = t(-1, 0, 0);
        UnitCell {
            a: 1.0,
            b: 1.0,
            c: 1.0,
            alpha: 0.0,
            beta: 0.0,
            gamma: 0.0,
            sites: vec![
                CellSite {
                    x: 0.0,
                    y: 0.0,
                    z: 0.0,
                    n: 0,
                    state: State(100),
                    nbr: [t(1, 1, 0), t(1, 0, 1), t(1, 0, -1), none, none, none],
                },
                CellSite {
                    x: 0.5,
                    y: 0.5,
                    z: 0.0,
                    n: 1,
                    state: State(100),
                    nbr: [t(0, -1, 0), t(0, 0, -1), t(0, 0, 1), none, none, none],
                },
            ],
        }
    }

    #[test]
    fn index_scheme_matches_the_cpp_formula() {
        let uc = tiny_cell();
        let s = create_lattice(
            &uc,
            LatticeParams {
                a_cells: 4,
                b_cells: 3,
                surface_plane: 0,
            },
        );
        assert_eq!(s.graph.len(), 4 * 3 * 2);
        // coord of index a*(b_cells*npos) + b*npos + n
        let (a, b, n) = (2usize, 1usize, 1usize);
        let i = a * (3 * 2) + b * 2 + n;
        assert_eq!(s.coord[i], CellCoord { a: 2, b: 1, n: 1 });
    }

    #[test]
    fn surface_plane_0_is_periodic_in_a_open_in_b() {
        let uc = tiny_cell();
        let p = LatticeParams {
            a_cells: 4,
            b_cells: 3,
            surface_plane: 0,
        };
        let s = create_lattice(&uc, p);
        let idx = |a: i32, b: i32, n: i32| (a * 6 + b * 2 + n) as usize;
        // +a from the last a-cell wraps to a=0 (periodic).
        assert_eq!(s.graph.sites[idx(3, 0, 0)].nbr[0], Some(idx(0, 0, 1)));
        // +b from the last b-cell falls off the slab (open).
        assert_eq!(s.graph.sites[idx(0, 2, 0)].nbr[1], None);
        // -b from b=0 falls off too.
        assert_eq!(s.graph.sites[idx(0, 0, 0)].nbr[2], None);
        // Interior +b resolves normally.
        assert_eq!(s.graph.sites[idx(0, 1, 0)].nbr[1], Some(idx(0, 2, 1)));
        // Template slot with n < 0 is None everywhere.
        assert_eq!(s.graph.sites[idx(1, 1, 0)].nbr[3], None);
    }

    #[test]
    fn surface_plane_1_swaps_the_axes() {
        let uc = tiny_cell();
        let p = LatticeParams {
            a_cells: 4,
            b_cells: 3,
            surface_plane: 1,
        };
        let s = create_lattice(&uc, p);
        let idx = |a: i32, b: i32, n: i32| (a * 6 + b * 2 + n) as usize;
        // Now +a from the last a-cell is the open direction...
        assert_eq!(s.graph.sites[idx(3, 0, 0)].nbr[0], None);
        // ...and b wraps.
        assert_eq!(s.graph.sites[idx(0, 0, 0)].nbr[2], Some(idx(0, 2, 1)));
    }

    #[test]
    fn build_is_deterministic() {
        // The premise of the whole M3 gate, pinned cheaply here.
        let uc = tiny_cell();
        let p = LatticeParams {
            a_cells: 4,
            b_cells: 3,
            surface_plane: 0,
        };
        assert_eq!(create_lattice(&uc, p), create_lattice(&uc, p));
    }

    #[test]
    fn populate_fill_fraction_truncates_in_f32_like_the_cpp() {
        let uc = tiny_cell();
        let p = LatticeParams {
            a_cells: 4,
            b_cells: 3,
            surface_plane: 0,
        };

        // Undersaturated: frac 0.7, (3 * 0.7f32) as i32 == 2 layers.
        let mut s = create_lattice(&uc, p);
        let r = populate_solid(&mut s, -1.0, -1.0);
        assert_eq!(
            r,
            PopulateReport {
                condition: "undersaturated",
                filled_cells: 2
            }
        );
        let occupied = s
            .graph
            .sites
            .iter()
            .filter(|x| x.state.is_occupied())
            .count();
        assert_eq!(occupied, 4 * 2 * 2); // all a-cells, first 2 b-layers, 2 sites each

        // Supersaturated: frac 0.3, (3 * 0.3f32) as i32 == 0 layers.
        let mut s = create_lattice(&uc, p);
        let r = populate_solid(&mut s, 1.0, 1.0);
        assert_eq!(
            r,
            PopulateReport {
                condition: "supersaturated",
                filled_cells: 0
            }
        );

        // Near equilibrium: frac 0.5 → 1 layer.
        let mut s = create_lattice(&uc, p);
        let r = populate_solid(&mut s, 0.0, 0.0);
        assert_eq!(
            r,
            PopulateReport {
                condition: "near equilibrium",
                filled_cells: 1
            }
        );

        // The truncation lives at the mercy of f32 rounding, either way:
        // 3 * 0.7f32 lands BELOW the integer (2.0999999 → 2), while
        // 20 * 0.7f32 rounds UP to exactly 14.0 → 14. Pinned so nobody
        // "simplifies" the arithmetic and shifts a slab boundary.
        assert_eq!((3f32 * 0.7f32) as i32, 2);
        assert_eq!((20f32 * 0.7f32) as i32, 14);
    }

    /// Hand-built 4-site double bridge: two Al (0, 1) bridged by two
    /// 400-class oxygens (2, 3).
    fn bridge_structure(al1_nbrs: [Option<SiteId>; 6]) -> Structure {
        let o = |a, b| Site {
            state: State(400),
            nbr: [Some(a), Some(b), None, None, None, None],
        };
        let sites = vec![
            Site {
                state: State(100),
                nbr: al1_nbrs,
            },
            Site {
                state: State(100),
                nbr: [Some(2), Some(3), None, None, None, None],
            },
            o(0, 1),
            o(0, 1),
        ];
        let n = sites.len();
        Structure {
            graph: SiteGraph { sites },
            coord: vec![CellCoord { a: 0, b: 0, n: 0 }; n],
            pair: vec![None; n],
            lostal: vec![None; n],
            params: LatticeParams {
                a_cells: 1,
                b_cells: 1,
                surface_plane: 0,
            },
        }
    }

    #[test]
    fn find_pairs_links_the_double_bridge() {
        // al1's neighbor list presents the PARTNER oxygen (3) before o1
        // itself (2), so the scan finds the real partner.
        let mut s = bridge_structure([Some(3), Some(2), None, None, None, None]);
        find_pairs(&mut s);
        assert_eq!(s.pair[2], Some(3));
        assert_eq!(s.pair[3], Some(2));
    }

    #[test]
    fn find_pairs_survives_meeting_o1_first() {
        // The ported C++ quirk in action: al1's list presents o1 (2) before
        // its partner (3), so o1 transiently pairs with ITSELF
        // (pair[2] = 2). Site 3 is still unpaired, scans, hits 2 first,
        // matches, and overwrites: pair[3] = 2, pair[2] = 3. Net effect in
        // this topology: same link, ugly route — identical to the C++,
        // which is the point.
        let mut s = bridge_structure([Some(2), Some(3), None, None, None, None]);
        find_pairs(&mut s);
        assert_eq!(s.pair[2], Some(3));
        assert_eq!(s.pair[3], Some(2));
    }
}
