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
                sites.push(Site { state: uc.sites[n].state, nbr });
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
                    x: 0.0, y: 0.0, z: 0.0, n: 0, state: State(100),
                    nbr: [t(1, 1, 0), t(1, 0, 1), t(1, 0, -1), none, none, none],
                },
                CellSite {
                    x: 0.5, y: 0.5, z: 0.0, n: 1, state: State(100),
                    nbr: [t(0, -1, 0), t(0, 0, -1), t(0, 0, 1), none, none, none],
                },
            ],
        }
    }

    #[test]
    fn index_scheme_matches_the_cpp_formula() {
        let uc = tiny_cell();
        let s = create_lattice(&uc, LatticeParams { a_cells: 4, b_cells: 3, surface_plane: 0 });
        assert_eq!(s.graph.len(), 4 * 3 * 2);
        // coord of index a*(b_cells*npos) + b*npos + n
        let (a, b, n) = (2usize, 1usize, 1usize);
        let i = a * (3 * 2) + b * 2 + n;
        assert_eq!(s.coord[i], CellCoord { a: 2, b: 1, n: 1 });
    }

    #[test]
    fn surface_plane_0_is_periodic_in_a_open_in_b() {
        let uc = tiny_cell();
        let p = LatticeParams { a_cells: 4, b_cells: 3, surface_plane: 0 };
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
        let p = LatticeParams { a_cells: 4, b_cells: 3, surface_plane: 1 };
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
        let p = LatticeParams { a_cells: 4, b_cells: 3, surface_plane: 0 };
        assert_eq!(create_lattice(&uc, p), create_lattice(&uc, p));
    }
}
