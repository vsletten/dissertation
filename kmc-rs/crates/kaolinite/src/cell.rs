//! The unit cell — the crystallographic *motif* that tiles into the lattice.
//! Port of the data structures in `ucell.hpp`.
//!
//! `data.cell` describes `Npos` positions (26 in the golden cell), each with
//! coordinates, an initial state class, and a 6-slot **neighbor template**:
//! "my j-th neighbor is position `n` of the unit cell offset `(a, b, c)`
//! cells away". The lattice build (M2) tiles this motif `aCells × bCells`
//! times and resolves every template entry into a concrete site index.
//!
//! The third-axis offset `c` is read and stored but never participates in
//! neighbor resolution — the simulated surface is a single 2D sheet (spec
//! A2). Kept because the file has it and faithfulness starts at the parse.

use crate::state::State;

/// One slot of a position's neighbor template (C++ `struct NeighborSite`).
///
/// `n < 0` means "no j-th neighbor" — the C++ `-1`-as-null convention. We
/// keep the raw `i32` *here*, at the data boundary, because that is what
/// the file contains; the lattice build converts it to `Option<SiteId>`
/// where the meaning actually matters (see M2).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct NeighborTemplate {
    /// Target position within the neighboring unit cell (`< 0` = none).
    pub n: i32,
    /// Unit-cell offset along the a axis.
    pub a: i32,
    /// Unit-cell offset along the b axis.
    pub b: i32,
    /// Unit-cell offset along the c axis (read, never used — 2D sheet).
    pub c: i32,
}

/// One position of the unit-cell motif (C++ `struct CellSite`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CellSite {
    /// Coordinates of the site within the unit cell, as given in
    /// `data.cell`. (The MSI writer divides these by hard-coded cell-vector
    /// lengths before use — see spec B5; they are *not* clean fractional
    /// coordinates.)
    pub x: f32,
    /// See `x`.
    pub y: f32,
    /// See `x`.
    pub z: f32,
    /// Site number as recorded in the file (redundant with the array
    /// position, but the C++ reads and stores it, so we do too).
    pub n: i32,
    /// Initial state class for sites stamped from this position (100, 200,
    /// 300, 400, or 500 in the golden cell).
    pub state: State,
    /// The 6-slot neighbor template.
    ///
    /// \[IDIOM\] `[NeighborTemplate; 6]` — a fixed-size array whose length is
    /// part of the **type**. The C++ `struct NeighborSite nbr[6]` gives the
    /// same layout but decays to a bare pointer the moment it's passed
    /// anywhere; the Rust array never loses its length, and every index is
    /// bounds-checked (`nbr[7]` panics deterministically instead of reading
    /// a neighbor struct that isn't there). Guarantee worth naming: Rust
    /// has no silent buffer over-read — the exact bug class of spec B7.
    pub nbr: [NeighborTemplate; 6],
}

/// The unit cell: dimensions plus the motif positions (C++ `class UnitCell`).
#[derive(Debug, Clone, PartialEq)]
pub struct UnitCell {
    /// Cell edge a (Å). Read from `data.cell` and — WART (spec B5) — never
    /// used by the output writers, which carry their own hard-coded cell
    /// matrix. Stored because parsing it is part of the file contract.
    pub a: f32,
    /// Cell edge b (Å). Same B5 caveat.
    pub b: f32,
    /// Cell edge c (Å). Same B5 caveat.
    pub c: f32,
    /// Cell angle alpha (radians, per the file's own comment). Same B5 caveat.
    pub alpha: f32,
    /// Cell angle beta (radians). Same B5 caveat.
    pub beta: f32,
    /// Cell angle gamma (radians). Same B5 caveat.
    pub gamma: f32,
    /// The motif positions, index = position number `n`.
    ///
    /// \[IDIOM\] `Vec<CellSite>` replaces `new CellSite[Npos + 1]` + manual
    /// `delete[]`. Ownership is the point: the Vec's memory lives exactly
    /// as long as the `UnitCell` that owns it, freed automatically when it
    /// drops — there is no `DisposeUnitCell` to forget to call and no
    /// double-free to write. The C++ even allocates one *extra* slot to
    /// plant an `n = -1` end sentinel; a Vec knows its own length, so the
    /// sentinel (and the off-by-one bugs it invites) has nothing to do.
    pub sites: Vec<CellSite>,
}

impl UnitCell {
    /// Number of positions in the motif (C++ `GetNumPositions`).
    pub fn npos(&self) -> usize {
        self.sites.len()
    }
}
