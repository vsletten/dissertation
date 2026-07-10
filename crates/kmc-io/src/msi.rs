//! Cerius2 MSI writer — the viewer's input format and the M3 golden
//! artifact. Port of `output::writeMSI` (output.cpp), bug-compatible.
//!
//! Format sketch (spec A9.1): a parenthesized object tree; one `Model`
//! containing `Atom` objects then `Bond` objects, all sharing a single
//! object-ID counter that starts at 2 (`nthing` in the C++). Only occupied,
//! non-EDGE sites are emitted. Bonds appear only when `drawbonds` is set,
//! once per undirected pair (`i2 < i`), with bonds across the periodic
//! seam suppressed so the rendered structure doesn't draw springs across
//! the box.
//!
//! Preserved warts, each marked at its line:
//! * **WART (spec B5)** — coordinates come from a hard-coded 3×3 cell
//!   matrix, NOT from the `data.cell` dimensions the reader dutifully
//!   parsed. Changing `data.cell` changes the physics but not the picture.
//! * **"9 F"** — Al-OH-Al bridge oxygens render as fluorine so the viewer
//!   colors them distinctly. A rendering hack, not chemistry (spec A9.1).
//! * The `ACL` prefix is written *before* the element switch; a site class
//!   outside 0–5 would emit a malformed half-line, exactly like the C++
//!   `default: break` does. Unreachable with real data; ported anyway.
//!
//! # Numeric care (the bitwise gate)
//!
//! Every arithmetic step below is `f32`, in the C++'s exact evaluation
//! order, because the golden file was produced by `g++ -O3 -ffast-math`
//! from float expressions:
//!
//! * The basis lengths `al`/`bl`/`cl` are compile-time constants in both
//!   implementations (g++ constant-folds them with correct per-operation
//!   rounding; rustc evaluates our `f32` ops the same way at runtime).
//! * The x/y/z expressions are sums of three/two/one `f32` products,
//!   evaluated left to right. `-ffast-math` *licenses* g++ to reassociate,
//!   but for these short scalar chains it emits the written order — which
//!   the byte-identical golden diff then certifies. The one fast-math
//!   transform g++ DID apply is `-freciprocal-math` (divide → multiply by
//!   reciprocal); see `site_xyz` for the ulp story. If a future compiler
//!   re-capture changes association, the gate (not silent drift) will say
//!   so.
//! * No FMA: the golden binary targeted baseline x86-64 (SSE2, no
//!   `-march`), so multiply-add pairs round twice; Rust's `a * b + c`
//!   does too (it never contracts implicitly).

use std::io::{self, Write};

use kaolinite::build::Structure;
use kaolinite::cell::UnitCell;
use kaolinite::state::State;

use crate::fmt::format_g6;

/// WART (spec B5): the hard-coded Cartesian cell matrix, duplicated in
/// three writers in the C++ (`writeMSI`, `writeXYZ`, `writeSurf` each carry
/// their own copy). Rust gets one copy with a name; the *values* stay
/// magic, faithfully. Note row 0 disagrees with `data.cell`'s a = 5.140 —
/// this matrix is presumably a once-correct render cell that outlived its
/// inputs.
const CD: [[f32; 3]; 3] = [
    [4.9725, -0.0262374, -1.3362],
    [0.0, 8.92893, -0.30084],
    [0.0, 0.0, 7.384],
];

/// The derived basis lengths, exactly as the C++ computes them per call
/// (f32 sums of f32 products; `sqrt` on the promoted double, narrowed back
/// — bit-identical to `sqrtf` for these values, see docs/RUST_TOUR.md).
fn basis_lengths() -> (f32, f32, f32) {
    let al =
        (((CD[0][0] * CD[0][0] + CD[0][1] * CD[0][1] + CD[0][2] * CD[0][2]) as f64).sqrt()) as f32;
    let bl = (((CD[1][1] * CD[1][1] + CD[1][2] * CD[1][2]) as f64).sqrt()) as f32;
    let cl = CD[2][2];
    (al, bl, cl)
}

/// Site coordinates in the render basis — the f32 expression from
/// output.cpp lines 195–200, term order preserved, with ONE deliberate
/// deviation from the written C++:
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
fn site_xyz(st: &Structure, uc: &UnitCell, i: usize, al: f32, bl: f32, cl: f32) -> (f32, f32, f32) {
    let c = st.coord[i];
    let cs = &uc.sites[c.n as usize];
    let (a, b) = (c.a as f32, c.b as f32);
    // -freciprocal-math, reproduced literally.
    let (ral, rbl, rcl) = (1.0f32 / al, 1.0f32 / bl, 1.0f32 / cl);
    let x = (cs.x * ral + a) * CD[0][0] + (cs.y * rbl + b) * CD[0][1] + (cs.z * rcl) * CD[0][2];
    let y = (cs.y * rbl + b) * CD[1][1] + (cs.z * rcl) * CD[1][2];
    let z = (cs.z * rcl) * CD[2][2];
    (x, y, z)
}

/// Write the structure as Cerius2 MSI to any sink.
///
/// `name` is the model label ("start", "end", "step42000"); the caller owns
/// file naming (`{name}.msi`) and creation. `drawbonds` is the C++ int flag
/// (nonzero = emit bonds).
///
/// \[IDIOM\] `W: Write` — a *generic* sink, Rust's answer to "takes an
/// ostream&". The golden test hands this a `Vec<u8>` and diffs bytes in
/// memory; the CLI hands it a `File`. Traits make the seam: the function
/// body can only call what `Write` declares, so swapping the sink cannot
/// change behavior. (Reviewer's counterpoint: a trait with exactly one
/// implementor in the codebase is usually premature abstraction — `Write`
/// earns its keep because std provides the many implementors.)
///
/// \[IDIOM\] `io::Result<()>` + `writeln!(...)?` on every line: each write
/// can fail (disk full is real), and the `?` forest makes that visible
/// without drowning the format logic. Compare the C++, which checks no
/// stream write anywhere and would silently truncate on a full disk.
pub fn write_msi<W: Write>(
    w: &mut W,
    st: &Structure,
    uc: &UnitCell,
    name: &str,
    drawbonds: i32,
) -> io::Result<()> {
    let (al, bl, cl) = basis_lengths();
    let natom = st.graph.len();
    let (acells, bcells) = (st.params.a_cells, st.params.b_cells);

    writeln!(w, "# MSI CERIUS2 DataModel File Version 3 5")?;
    writeln!(w, "(1 Model")?;
    writeln!(w, "(A I Id 1)")?;
    writeln!(w, " (A C Label \"{name}\")")?;

    // Object IDs: atoms first, then bonds, one shared counter starting
    // after the Model's id 1. `id[i]` remembers each atom's object id for
    // the bond pass — the C++ mallocs this uninitialized; only occupied
    // sites' entries are ever written OR read, so zero-init is
    // observationally identical (and Rust wouldn't sell us uninitialized
    // memory at this price anyway).
    let mut nthing: i32 = 1;
    let mut id = vec![0i32; natom];

    // [IDIOM] `.iter_mut().enumerate()` where the C++ writes `for (i = 0;
    // i < natom; i++) ... id[i] = ...`. Clippy (needless_range_loop) pushes
    // this rewrite because iterating the slice directly proves every access
    // in-bounds once, instead of bounds-checking id[i] per hit — and it
    // makes the loop's write target part of its declaration.
    for (i, id_slot) in id.iter_mut().enumerate() {
        let s = st.graph.sites[i].state;
        // Occupied and not EDGE. (EDGE is 9: 9 % 100 > 0, so the explicit
        // EDGE test is load-bearing — see State::is_occupied's warning.)
        if s.0 % 100 > 0 && s != State::EDGE {
            nthing += 1;
            writeln!(w, " ({nthing} Atom")?;
            *id_slot = nthing;
            // Element line. The prefix goes out before the class switch —
            // C++ structure preserved (see module docs on the malformed-
            // half-line behavior for impossible classes).
            write!(w, "  (A C ACL ")?;
            match s.0 / 100 {
                0 => writeln!(w, "\"10 Ne\")")?, // "Edge" class — unreachable, ported
                1 => writeln!(w, "\"13 Al\")")?,
                2 => writeln!(w, "\"14 Si\")")?,
                3 | 4 => {
                    if s.0 % 100 > 1 {
                        writeln!(w, "\"1 H\")")?; // protonated O renders as H
                    } else {
                        writeln!(w, "\"8 O\")")?;
                    }
                }
                5 => {
                    if s.0 % 100 > 1 {
                        writeln!(w, "\"1 H\")")?;
                    } else {
                        writeln!(w, "\"9 F\")")?; // the fluorine rendering hack
                    }
                }
                _ => {} // C++ default: break (half-written line, faithfully)
            }
            let (x, y, z) = site_xyz(st, uc, i, al, bl, cl);
            writeln!(
                w,
                "  (A D XYZ ({} {} {}))",
                format_g6(x),
                format_g6(y),
                format_g6(z)
            )?;
            writeln!(w, "  (A I Id {nthing})")?;
            write!(w, "  (A C Label \"")?;
            match s.0 / 100 {
                0 => write!(w, "Edge")?,
                1 => write!(w, "Al")?,
                2 => write!(w, "Si")?,
                3 | 4 => write!(w, "O")?,
                5 => write!(w, "OH")?,
                _ => {}
            }
            writeln!(w, "{i}\")")?;
            writeln!(w, "  (A I LabelType 0)")?;
            writeln!(w, " )")?;
        }
    }

    // Bond pass: every undirected neighbor pair once (i2 < i), both ends
    // occupied and non-EDGE, minus the periodic-seam suppressions: a bond
    // whose template offset wrapped around the box (site at a == 0 linking
    // to a == acells-1 via a Δa = -1 template, and the three mirror cases)
    // is real topology but would render as a spring across the image.
    for i in 0..natom {
        let s = st.graph.sites[i].state;
        if drawbonds != 0 && s.0 % 100 > 0 && s != State::EDGE {
            for j in 0..6 {
                let Some(i2) = st.graph.sites[i].nbr[j] else {
                    continue;
                };
                if i2 >= i {
                    continue;
                }
                let s2 = st.graph.sites[i2].state;
                if s2.0 % 100 == 0 || s2 == State::EDGE {
                    continue;
                }
                let t = uc.sites[st.coord[i].n as usize].nbr[j];
                let (ca, cb) = (st.coord[i], st.coord[i2]);
                let wraps = (ca.a == 0 && cb.a == acells - 1 && t.a == -1)
                    || (ca.a == acells - 1 && cb.a == 0 && t.a == 1)
                    || (ca.b == 0 && cb.b == bcells - 1 && t.b == -1)
                    || (ca.b == bcells - 1 && cb.b == 0 && t.b == 1);
                if wraps {
                    continue;
                }
                nthing += 1;
                writeln!(w, " ({nthing} Bond")?;
                writeln!(w, "  (A O Atom1 {})", id[i])?;
                writeln!(w, "  (A O Atom2 {})", id[i2])?;
                writeln!(w, " )")?;
            }
        }
    }

    writeln!(w, ")")?;
    Ok(())
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
        // are certified by the golden gate.
        assert!((al - 5.149).abs() < 1e-3);
        assert!((bl - 8.934).abs() < 1e-3);
    }
}
