//! Cerius2 MSI writer — the viewer's input format, the M3 golden artifact,
//! and (M7) the movie-frame and end-state snapshots. Port of
//! `output::writeMSI` (output.cpp), bug-compatible.
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
//!   matrix (shared with the other writers via [`crate::view`]), NOT from
//!   the `data.cell` dimensions the reader dutifully parsed. Changing
//!   `data.cell` changes the physics but not the picture.
//! * **"9 F"** — Al-OH-Al bridge oxygens render as fluorine so the viewer
//!   colors them distinctly. A rendering hack, not chemistry (spec A9.1).
//! * The `ACL` prefix is written *before* the element switch; a site class
//!   outside 0–5 would emit a malformed half-line, exactly like the C++
//!   `default: break` does. Unreachable with real data; ported anyway.
//!
//! # Numeric care (the bitwise gates)
//!
//! Every arithmetic step below is `f32`, in the C++'s exact evaluation
//! order, because the golden files were produced by `g++ -O3 -ffast-math`
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
//!   reciprocal); see `view::site_xyz` for the ulp story. If a
//!   future compiler re-capture changes association, the gate (not silent
//!   drift) will say so.
//! * No FMA: the golden binary targeted baseline x86-64 (SSE2, no
//!   `-march`), so multiply-add pairs round twice; Rust's `a * b + c`
//!   does too (it never contracts implicitly).

use std::io::{self, Write};

use kaolinite::state::State;

use crate::fmt::format_g6;
use crate::view::{LatticeView, basis_lengths, site_xyz};

/// Write the lattice as Cerius2 MSI to any sink.
///
/// `name` is the model label ("start", "end", "step5000"); the caller owns
/// file naming (`{name}.msi`) and creation. `drawbonds` is the C++ int flag
/// (nonzero = emit bonds).
///
/// \[IDIOM\] `W: Write` — a *generic* sink, Rust's answer to "takes an
/// ostream&". The golden tests hand this a `Vec<u8>` and diff bytes in
/// memory; the CLI hands it a buffered `File`. Traits make the seam: the
/// function body can only call what `Write` declares, so swapping the sink
/// cannot change behavior. (Reviewer's counterpoint: a trait with exactly
/// one implementor in the codebase is usually premature abstraction —
/// `Write` earns its keep because std provides the many implementors.)
///
/// \[IDIOM\] `io::Result<()>` + `writeln!(...)?` on every line: each write
/// can fail (disk full is real), and the `?` forest makes that visible
/// without drowning the format logic. Compare the C++, which checks no
/// stream write anywhere and would silently truncate on a full disk.
pub fn write_msi<W: Write>(
    w: &mut W,
    view: &LatticeView<'_>,
    name: &str,
    drawbonds: i32,
) -> io::Result<()> {
    let (al, bl, cl) = basis_lengths();
    let natom = view.graph.len();
    let (acells, bcells) = (view.params.a_cells, view.params.b_cells);

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
        let s = view.graph.sites[i].state;
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
            let (x, y, z) = site_xyz(view, i, al, bl, cl);
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
        let s = view.graph.sites[i].state;
        if drawbonds != 0 && s.0 % 100 > 0 && s != State::EDGE {
            for j in 0..6 {
                let Some(i2) = view.graph.sites[i].nbr[j] else {
                    continue;
                };
                if i2 >= i {
                    continue;
                }
                let s2 = view.graph.sites[i2].state;
                if s2.0 % 100 == 0 || s2 == State::EDGE {
                    continue;
                }
                let t = view.cell.sites[view.coord[i].n as usize].nbr[j];
                let (ca, cb) = (view.coord[i], view.coord[i2]);
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
