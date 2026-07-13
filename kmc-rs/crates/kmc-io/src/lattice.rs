//! `data.lattice` reader — lattice dimensions. Port of the read at the top
//! of `Lattice::CreateLattice` (lattice.cpp).
//!
//! The file is three integers and then free text the C++ never reads:
//!
//! ```text
//! 20  3   0
//! acells  bcells  ac/bc surface plane
//! ```
//!
//! Yes — the "column headers" are on the line *after* the values, and they
//! are not a `#` comment; the C++ simply closes the file after the third
//! integer. We read exactly three values and stop, same as the original.

use std::path::Path;

use kaolinite::build::LatticeParams;

use crate::error::ReadError;
use crate::scan::Scanner;

/// Read `data.lattice` into [`LatticeParams`].
pub fn read_lattice(path: &Path) -> Result<LatticeParams, ReadError> {
    let mut s = Scanner::open(path)?;
    let a_cells = s.next_i32("invalid number of a cells")?;
    let b_cells = s.next_i32("invalid number of b cells")?;
    let surface_plane = s.next_i32("invalid surface plane")?;
    Ok(LatticeParams {
        a_cells,
        b_cells,
        surface_plane,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn golden_lattice_is_20_by_3_plane_0() {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../data/golden/inputs/data.lattice"
        );
        let p = read_lattice(Path::new(path)).unwrap();
        assert_eq!(
            p,
            LatticeParams {
                a_cells: 20,
                b_cells: 3,
                surface_plane: 0
            }
        );
    }
}
