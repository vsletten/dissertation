//! The site-state encoding — the model's entire vocabulary.
//!
//! A site's `state` is a 3-digit integer: the hundreds digit is the site
//! *class* (what kind of atom position this is), the last two digits are
//! occupancy / coordination / protonation. `state % 100 == 0` means the
//! site is empty. The special value 9 ([`State::EDGE`]) marks frozen
//! boundary sites.
//!
//! The table below is copied **verbatim** from `common.hpp` (the spec calls
//! it "the authoritative key", A3) so this crate carries its own ground
//! truth:
//!
//! ```text
//! EDGE   used to identify site at edge of lattice
//!
//! Al sites            Si sites           Si-O-Si type O sites
//! --------            --------           --------------------
//! 100  empty          200  empty         300  empty
//! 101  Al(OH,H2O)0    201  Si(OH)0       301  Si-O-Si
//! 102  Al(OH,H2O)1    202  Si(OH)1       302  Si-OH HO-Si
//! 103  Al(OH,H2O)2    203  Si(OH)2       303  Si-OH
//! 104  Al(OH,H2O)3    204  Si(OH)3
//! 105  Al(OH,H2O)4    205  Si(OH)4
//! 106  Al(OH,H2O)5    299  Al(OH,H2O)6   (wrong-cation, never produced)
//! 107  Al(OH,H20)6
//! 199  Si(OH)4        (wrong-cation, never produced)
//!
//! Si-O<Al2 type O sites                  Al-OH-Al type O(H) sites
//! ---------------------                  ------------------------
//! 400  empty                             500  empty
//! 401  Si-O<Al2                          501  Al-OH-Al
//! 402  Si-OH HO<Al2                      502  Al-OH H2O-Al
//! 403  Si-OH HO-Al H2O-Al                503  Al-OH
//! 404        HO<Al2
//! 405        HO-Al H2O-Al
//! 406  Si-OH-Al
//! 407  Si-OH HO-Al
//! 408  Si-OH
//! 409        HO-Al
//! 410  Si-OH-Al HO-Al
//! ```
//!
//! # Why not a Rust `enum`?
//!
//! The design doc (§8) floats an `enum SiteClass` so `match` can be
//! exhaustive. We deliberately stop at a **newtype over the raw code**
//! instead, because the legacy model does *arithmetic* on states — it
//! increments them to occupy a site (`PopulateSolid`: 100→101, 400→401),
//! compares them with `<` (`FindPairs`: `state < 400`), and bumps cation
//! coordination with `state++`. An enum would force a lossy translation
//! layer through the very code the M3 bitwise gate must prove faithful.
//! The newtype keeps the arithmetic honest while still giving us a place
//! to hang meaning (`class_code`, `is_occupied`) and making it a *type
//! error* to confuse a state with a site index — the exact bug class the
//! all-`int` C++ left open. Enum-per-class remains the right move for a
//! *reformed* model (M8+), not a faithful port.

/// A site state code (newtype over the raw 3-digit integer).
///
/// \[IDIOM\] The **newtype pattern**: `struct State(pub i32)` is a distinct
/// type that costs nothing at runtime (same size, same registers as `i32`)
/// but partitions the world at compile time. The C++ site record holds
/// `int state`, `int pair`, `int nbr[6]` — three integers with three
/// meanings and nothing stopping `sites[state]` or `state = pair`. Here,
/// using a `State` where a site index belongs (or vice versa) refuses to
/// compile. This is the cheapest safety Rust sells, and for a reviewer the
/// tell is inverted: bare `i32`/`usize` fields whose *names* carry the
/// meaning are the smell; newtypes are the fix.
///
/// The inner value is `pub` on purpose: the faithful port needs literal
/// access (`state.0 < 400`, `state.0 += 1`) to stay line-for-line with the
/// C++, and hiding it behind accessors would only add distance between
/// this code and the spec it must mirror.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct State(pub i32);

impl State {
    /// Frozen-boundary sentinel (`#define EDGE 9` in `lattice.hpp`).
    /// EDGE sites never react and are skipped in output; they are the bulk
    /// crystal below the simulated surface.
    pub const EDGE: State = State(9);

    /// The class digit: 1 Al, 2 Si, 3 Si-O-Si, 4 Si-O<Al2, 5 Al-OH-Al.
    /// (Port of the `TYPE(A)` macro, `state / 100`.) Note `EDGE` (9) maps
    /// to class 0 — the C++ leans on that in output switches.
    pub fn class_code(self) -> i32 {
        self.0 / 100
    }

    /// The occupancy/coordination digits (`state % 100`).
    pub fn occupancy(self) -> i32 {
        self.0 % 100
    }

    /// Is the site occupied? (Port of the `ISOCC` macro.) Careful at the
    /// call sites: `EDGE` (9) counts as "occupied" by this test — the C++
    /// always pairs it with an explicit EDGE check, and so do we.
    pub fn is_occupied(self) -> bool {
        self.occupancy() > 0
    }

    /// Is this the EDGE sentinel? (Port of the `ISEDGE` macro.)
    pub fn is_edge(self) -> bool {
        self == State::EDGE
    }

    /// How many neighbors a site of this class is *supposed* to have —
    /// port of `UnitCell::GetNumNeighbors` (ucell.cpp), used by
    /// `TerminateLattice` to detect boundary sites (actual count ≠
    /// expected count → freeze as EDGE). Returns -1 for classes outside
    /// 1–5, exactly like the C++ `default:` arm (EDGE lands here: 9/100=0).
    ///
    /// \[IDIOM\] `match` vs `switch`: no fallthrough to forget a `break` on,
    /// arms are expressions (the whole match *returns* the value), multiple
    /// patterns share an arm with `|`, and — had we matched on an enum —
    /// the compiler would reject a missing case. Matching on `i32` needs
    /// the `_` catch-all, which is precisely the C++ `default:` we are
    /// porting.
    pub fn expected_neighbor_count(self) -> i32 {
        match self.class_code() {
            1 => 6,     // Al
            2 => 4,     // Si
            3 | 5 => 2, // Si-O-Si or Al-OH-Al type O
            4 => 3,     // Si-O<Al2 type O
            _ => -1,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn class_and_occupancy_split_the_code() {
        assert_eq!(State(401).class_code(), 4);
        assert_eq!(State(401).occupancy(), 1);
        assert!(State(401).is_occupied());
        assert!(!State(400).is_occupied());
    }

    #[test]
    fn edge_is_class_zero_and_counts_as_occupied() {
        // Two traps the C++ code lives with; pin them so nobody "fixes"
        // them mid-port.
        assert_eq!(State::EDGE.class_code(), 0);
        assert!(State::EDGE.is_occupied());
        assert!(State::EDGE.is_edge());
    }

    #[test]
    fn expected_neighbor_counts_match_get_num_neighbors() {
        assert_eq!(State(100).expected_neighbor_count(), 6);
        assert_eq!(State(205).expected_neighbor_count(), 4);
        assert_eq!(State(301).expected_neighbor_count(), 2);
        assert_eq!(State(404).expected_neighbor_count(), 3);
        assert_eq!(State(503).expected_neighbor_count(), 2);
        assert_eq!(State::EDGE.expected_neighbor_count(), -1); // default arm
    }
}
