//! Surface-projection writer — port of `output::writeSurf` (output.cpp),
//! bug-compatible.
//!
//! Two CSV files of `x,y` rows (spec A9.4): `surfSi.out` for silicon,
//! `surfAl.out` for aluminum. A cation makes the cut when it is occupied
//! **beyond its minimal state** (Si 202–205, Al 102–107 — the `> 201` /
//! `> 101` bounds exclude the bare 201/101, faithfully) *and* at least one
//! neighbor oxygen is in a surface-exposed state (303, 404–409, or 503).
//! The result is a 2D map of where the reactive surface pokes through.
//!
//! Faithful subtleties:
//! * The neighbor scan is `for (j = 0; j < 6 && nbr[j] >= 0; j++)` — the
//!   guard *precedes* the read, and the loop stops at the first missing
//!   slot. This is the same loop shape whose guard-after-read twin in
//!   `IsActive` produced the `nbr[6]` phantom (REFORM_PLAN R1) — here the
//!   author wrote it correctly, and the port mirrors the correct one just
//!   as literally (`iter().map_while(...)` — yields neighbors while the
//!   slot is `Some`, stops at the first `None`).
//! * Only `x,y` are written; the z expression is commented out in the C++
//!   ("for 3D") and stays unwritten here.
//! * Both files are truncated fresh each call; row order is site order,
//!   routed to one file or the other by class.

use std::io::{self, Write};

use crate::fmt::format_g6;
use crate::view::{LatticeView, basis_lengths, site_xyz};

/// Oxygen states that mark a neighboring cation as surface-exposed:
/// 303, 404..=409, 503 (the C++'s `== 303 || (> 403 && < 410) || == 503`).
fn exposes_surface(s: i32) -> bool {
    matches!(s, 303 | 404..=409 | 503)
}

/// Write the Si and Al surface maps to two sinks.
///
/// \[IDIOM\] Two independent generic sinks, two type parameters. `si` and
/// `al` may be different concrete types (a `File` and a `Vec<u8>` in a
/// test) precisely because each gets its own parameter; a single `W` would
/// force them identical. The C++ signature — `void writeSurf(Lattice*)` —
/// names neither sink at all: the filenames are buried in the body, which
/// is why the C++ version is untestable without a filesystem.
pub fn write_surf<WS: Write, WA: Write>(
    si: &mut WS,
    al: &mut WA,
    view: &LatticeView<'_>,
) -> io::Result<()> {
    let (bal, bbl, bcl) = basis_lengths();

    for i in 0..view.graph.len() {
        let s = view.graph.sites[i].state.0;
        // Occupied beyond minimal: Si 202..=205 or Al 102..=107 (the C++'s
        // `> 201 && < 206` / `> 101 && < 108`, exclusive bounds preserved).
        if !matches!(s, 202..=205 | 102..=107) {
            continue;
        }
        // Stop at the first absent neighbor slot, like the C++'s
        // `j < 6 && nbr[j] >= 0` (see module docs).
        let exposed = view.graph.sites[i]
            .nbr
            .iter()
            .map_while(|n| *n)
            .any(|nbr| exposes_surface(view.graph.sites[nbr].state.0));
        if !exposed {
            continue;
        }
        let (x, y, _z) = site_xyz(view, i, bal, bbl, bcl);
        let row_target: &mut dyn Write = if s > 200 { si } else { al };
        writeln!(row_target, "{},{}", format_g6(x), format_g6(y))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exposure_states_match_the_cpp_ranges() {
        for s in [303, 404, 405, 406, 407, 408, 409, 503] {
            assert!(exposes_surface(s), "{s} should expose");
        }
        // Boundary fenceposts: 403 and 410 are OUTSIDE `> 403 && < 410`.
        for s in [302, 403, 410, 501, 502] {
            assert!(!exposes_surface(s), "{s} should not expose");
        }
    }
}
