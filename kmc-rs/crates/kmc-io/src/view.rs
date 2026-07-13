//! The writers' shared view of the lattice, plus the (bug-compatible)
//! render geometry all three C++ output writers duplicate.
//!
//! `output.cpp` has no such module: `writeMSI`, `writeXYZ`, and `writeSurf`
//! each take a `Lattice*` — a pointer to *everything* — and each carries its
//! own local copy of the hard-coded cell matrix (spec B5). The port keeps
//! the *values* faithfully but hosts them once, here, and narrows the
//! writers' input to exactly what they read.
//!
//! # Why a view struct exists at all
//!
//! Through M3 the writers took `&Structure` — fine while the CLI still owned
//! one. But at the dynamics handoff (`Kaolinite::from_structure`) the
//! `Structure` is *dismembered*: the graph moves to the engine's loop, and
//! `pair`/`lostal` move into the model. By the time `end.msi` is written
//! there is no `Structure` left to point at — only the live graph plus the
//! immutable tiling data the CLI kept aside. [`LatticeView`] is that
//! regrouping: borrowed pieces, reassembled without copying anything.

use kaolinite::build::{CellCoord, LatticeParams, Structure};
use kaolinite::cell::UnitCell;
use kaolinite::state::State;
use kmc_engine::SiteGraph;

/// Everything an output writer may look at, by reference.
///
/// \[IDIOM\] **A borrowed view struct** — a struct whose fields are `&'a`
/// references, the Rust counterpart of a C++ "parameter object" of
/// pointers, with one decisive difference: the lifetime `'a` ties the view
/// to its sources *in the type*. A `LatticeView` cannot outlive the graph
/// it points into, cannot be stashed in a global, and freezes its sources
/// for as long as it lives (shared borrows: nobody mutates the lattice
/// mid-write). The C++ equivalent — `Lattice*` — promises none of that.
///
/// Reviewer's lens: view structs are the honest fix when a function's
/// parameter list grows past ~4 borrowed things that always travel
/// together. The dishonest fixes to watch for in agent code are (a)
/// cloning the world into an owned struct "to make lifetimes easy", and
/// (b) passing the original fat object so the signature hides what is
/// actually read. A view names the real dependency set and costs nothing.
#[derive(Clone, Copy)]
pub struct LatticeView<'a> {
    /// States + topology — *live*: after the run this is the engine's
    /// mutated graph, not the structural build.
    pub graph: &'a SiteGraph<State>,
    /// Tiling coordinates per site (C++ `a`/`b`/`n` members), fixed at
    /// build time.
    pub coord: &'a [CellCoord],
    /// Lattice dimensions (for the periodic-seam bond suppression).
    pub params: LatticeParams,
    /// The unit-cell motif (fractional-ish coordinates per position).
    pub cell: &'a UnitCell,
}

impl<'a> LatticeView<'a> {
    /// View a still-assembled [`Structure`] (pre-dynamics: `start.msi`).
    pub fn of_structure(st: &'a Structure, cell: &'a UnitCell) -> Self {
        LatticeView {
            graph: &st.graph,
            coord: &st.coord,
            params: st.params,
            cell,
        }
    }
}

/// WART (spec B5): the hard-coded Cartesian cell matrix, duplicated in
/// three writers in the C++ (`writeMSI`, `writeXYZ`, `writeSurf` each carry
/// their own copy). Rust gets one copy with a name; the *values* stay
/// magic, faithfully. Note row 0 disagrees with `data.cell`'s a = 5.140 —
/// this matrix is presumably a once-correct render cell that outlived its
/// inputs.
pub(crate) const CD: [[f32; 3]; 3] = [
    [4.9725, -0.0262374, -1.3362],
    [0.0, 8.92893, -0.30084],
    [0.0, 0.0, 7.384],
];

/// The derived basis lengths, exactly as the C++ computes them per call
/// (f32 sums of f32 products; `sqrt` on the promoted double, narrowed back
/// — bit-identical to `sqrtf` for these values, see docs/RUST_TOUR.md).
pub(crate) fn basis_lengths() -> (f32, f32, f32) {
    let al =
        (((CD[0][0] * CD[0][0] + CD[0][1] * CD[0][1] + CD[0][2] * CD[0][2]) as f64).sqrt()) as f32;
    let bl = (((CD[1][1] * CD[1][1] + CD[1][2] * CD[1][2]) as f64).sqrt()) as f32;
    let cl = CD[2][2];
    (al, bl, cl)
}

/// Site coordinates in the render basis — the f32 expression the C++
/// repeats in all three writers (output.cpp), term order preserved, with
/// ONE deliberate deviation from the written C++:
///
/// **Numeric care (the hard-won ulp):** the source says `cs.x / al`, but
/// the golden binary was compiled with `-ffast-math`, which includes
/// `-freciprocal-math` — g++ hoisted the loop-invariant divisors and
/// emitted `cs.x * (1.0f/al)` instead. For most sites the two round
/// identically; for cell position 20 (and only through the 6-digit `%g`
/// window at atom "OH20") they differ by one ulp: division gives
/// -0.12292254 → "-0.122923", reciprocal-multiply gives -0.12292248 →
/// "-0.122922", and the golden file says -0.122922. So the port
/// reproduces the *compiled* semantics: reciprocals computed once
/// (correctly rounded f32 constants, exactly what g++ folds), then
/// multiplied. Writing `/ al` here would be truer to the C++ text and
/// wrong by one byte — faithfulness follows the binary that produced the
/// reference, not the source that suggested it.
pub(crate) fn site_xyz(
    view: &LatticeView<'_>,
    i: usize,
    al: f32,
    bl: f32,
    cl: f32,
) -> (f32, f32, f32) {
    let c = view.coord[i];
    let cs = &view.cell.sites[c.n as usize];
    let (a, b) = (c.a as f32, c.b as f32);
    // -freciprocal-math, reproduced literally.
    let (ral, rbl, rcl) = (1.0f32 / al, 1.0f32 / bl, 1.0f32 / cl);
    let x = (cs.x * ral + a) * CD[0][0] + (cs.y * rbl + b) * CD[0][1] + (cs.z * rcl) * CD[0][2];
    let y = (cs.y * rbl + b) * CD[1][1] + (cs.z * rcl) * CD[1][2];
    let z = (cs.z * rcl) * CD[2][2];
    (x, y, z)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basis_lengths_match_hand_computation() {
        // The same constants the golden capture was rendered with.
        let (al, bl, cl) = basis_lengths();
        assert_eq!(cl, 7.384);
        // al = |row 0|, bl = |(row 1).yz| — sanity envelope, exact values
        // are certified by the golden gates.
        assert!((al - 5.149).abs() < 1e-3);
        assert!((bl - 8.934).abs() < 1e-3);
    }
}
