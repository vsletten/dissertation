//! `data.cell` reader — the unit-cell motif. Port of `UnitCell::CreateUnitCell`
//! (ucell.cpp).
//!
//! File shape (whitespace-delimited, `#` comments only where the C++ calls
//! `EatComment`):
//!
//! ```text
//! 5.140  8.930  7.370        # a, b, c dimensions (angstroms)
//! -0.03141  -0.25038  0.0    # alpha, beta, gamma (radians)
//! 26                         # Number of positions in unit cell
//! #id type     x       y      z      nbr:  id  da db dc
//!  0  100   1.76027 4.49982 2.84548        8   0  0  0
//!                                          9   0  0  0
//!                                          ...(6 neighbor rows)...
//! ```
//!
//! Inside the position block there are **no** comments and the C++ performs
//! **no** `EatComment` calls — one `#` in there would kill the legacy parser.
//! Our scanner is faithful to that (see [`crate::scan`]).

use std::path::Path;

use kaolinite::{CellSite, NeighborTemplate, State, UnitCell};

use crate::error::ReadError;
use crate::scan::Scanner;

/// Read `data.cell` into a [`UnitCell`].
///
/// Same read order as the C++: `A B C`, comment, `alpha beta gamma`,
/// comment, `Npos`, comment, then for each position `n state x y z`
/// followed by six `n a b c` neighbor rows, with no comment handling until
/// the block ends.
pub fn read_cell(path: &Path) -> Result<UnitCell, ReadError> {
    let mut s = Scanner::open(path)?;
    parse_cell(&mut s)
}

/// The parse, on any scanner (tests feed strings).
fn parse_cell(s: &mut Scanner) -> Result<UnitCell, ReadError> {
    let a = s.next_f32("invalid cell dimension a")?;
    let b = s.next_f32("invalid cell dimension b")?;
    let c = s.next_f32("invalid cell dimension c")?;
    s.eat_comment();
    let alpha = s.next_f32("invalid cell angle alpha")?;
    let beta = s.next_f32("invalid cell angle beta")?;
    let gamma = s.next_f32("invalid cell angle gamma")?;
    s.eat_comment();
    let npos = s.next_i32("invalid number of unit cell positions")?;
    s.eat_comment();

    // [IDIOM] Building a Vec by pushing vs. `new CellSite[Npos + 1]`.
    // `with_capacity` pre-sizes the allocation (the one performance fact
    // worth knowing about Vec: growth is amortized doubling, and telling it
    // the size up front skips the copies). Note what we DON'T port: the
    // C++ allocates Npos+1 and writes an n = -1 sentinel into the extra
    // slot; a Vec carries its length, so the sentinel is structurally
    // meaningless here.
    let mut sites = Vec::with_capacity(npos.max(0) as usize);
    for _ in 0..npos {
        let n = s.next_i32("invalid unit cell site number")?;
        let state = State(s.next_i32("invalid unit cell site state")?);
        let x = s.next_f32("invalid unit cell site x")?;
        let y = s.next_f32("invalid unit cell site y")?;
        let z = s.next_f32("invalid unit cell site z")?;
        // [IDIOM] `Default` + array init, then fill in place. C++ default-
        // constructs the whole CellSite array and assigns fields one
        // stream-read at a time; Rust makes the "not yet meaningful" window
        // explicit and local — the array exists only inside this loop body
        // and is fully overwritten before anyone else can see it.
        let mut nbr = [NeighborTemplate::default(); 6];
        for slot in &mut nbr {
            slot.n = s.next_i32("invalid neighbor site number")?;
            slot.a = s.next_i32("invalid neighbor offset a")?;
            slot.b = s.next_i32("invalid neighbor offset b")?;
            slot.c = s.next_i32("invalid neighbor offset c")?;
        }
        sites.push(CellSite {
            x,
            y,
            z,
            n,
            state,
            nbr,
        });
    }

    Ok(UnitCell {
        a,
        b,
        c,
        alpha,
        beta,
        gamma,
        sites,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn golden() -> UnitCell {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../data/golden/inputs/data.cell"
        );
        read_cell(Path::new(path)).unwrap()
    }

    #[test]
    fn golden_cell_has_26_positions_and_the_documented_dims() {
        let uc = golden();
        assert_eq!(uc.npos(), 26);
        assert_eq!(uc.a, 5.140);
        assert_eq!(uc.b, 8.930);
        assert_eq!(uc.c, 7.370);
        assert_eq!(uc.alpha, -0.03141);
        assert_eq!(uc.beta, -0.25038);
        assert_eq!(uc.gamma, 0.0);
    }

    #[test]
    fn golden_cell_spot_checks() {
        let uc = golden();
        // Position 0: Al site at the documented coordinates.
        let s0 = &uc.sites[0];
        assert_eq!(s0.n, 0);
        assert_eq!(s0.state, State(100));
        assert_eq!((s0.x, s0.y, s0.z), (1.76027, 4.49982, 2.84548));
        assert_eq!(
            s0.nbr[0],
            NeighborTemplate {
                n: 8,
                a: 0,
                b: 0,
                c: 0
            }
        );
        assert_eq!(
            s0.nbr[5],
            NeighborTemplate {
                n: 25,
                a: 0,
                b: 0,
                c: 0
            }
        );
        // Position 2: neighbor with a nonzero cell offset (crosses cells).
        assert_eq!(
            uc.sites[2].nbr[1],
            NeighborTemplate {
                n: 19,
                a: 1,
                b: 1,
                c: 0
            }
        );
        // Position 4: first Si site.
        assert_eq!(uc.sites[4].state, State(200));
        assert_eq!(uc.sites[4].x, 0.308295);
    }

    #[test]
    fn class_mix_matches_the_kaolinite_stoichiometry() {
        // 26 positions: 4 Al + 4 Si + 18 oxygens of three flavors.
        let uc = golden();
        let count = |class: i32| {
            uc.sites
                .iter()
                .filter(|s| s.state.class_code() == class)
                .count()
        };
        // [IDIOM] Iterator chains (`iter().filter().count()`) — the Rust
        // for-loop-with-accumulator. Zero-cost: compiles to the same loop,
        // but the *intent* (count matching) is the code, not a pattern the
        // reviewer must reverse-engineer from mutation.
        assert_eq!(count(1) + count(2) + count(3) + count(4) + count(5), 26);
        assert_eq!(count(1), 4); // Al
        assert_eq!(count(2), 4); // Si
    }
}
