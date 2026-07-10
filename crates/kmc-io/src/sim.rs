//! `data.sim` reader — run parameters. Port of `sim.cpp` / `sim.hpp`.
//!
//! The file, as documented in its own comments, is five lines:
//!
//! ```text
//! 20000     # Number of steps in simulation
//! 1000      # Number of steps between data writes (0 for no data)
//! 1000000   # Number of steps between movie frames (0 for no movie)
//! -2        # Seed for random number generator
//! 1         # Draw bonds? (0 No, 1 Yes)
//! ```
//!
//! The C++ does **not** read it that way — see the WART note in [`read_sim`].

use std::path::Path;

use crate::error::ReadError;
use crate::scan::Scanner;

/// Run parameters from `data.sim` (C++ `class Simulation`).
///
/// \[IDIOM\] Plain data as a `struct` with public fields — no getters, no
/// constructor boilerplate. Rust structs are not classes: there is no
/// inheritance, and "just data" is idiomatic where C++ style pressure would
/// demand accessors. Behavior, when it exists, goes in `impl` blocks, and
/// *visibility* (`pub` or not) is the encapsulation tool, not method
/// indirection.
///
/// Field types are `i32` on purpose: the C++ members are `int`, and the port
/// keeps integer widths explicit so overflow behavior can't silently differ.
/// (Rust makes you say `i32`/`i64`; there is no width-varies-by-platform
/// `int` to inherit surprises from.)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SimParams {
    /// Total number of MC steps to run.
    pub nsteps: i32,
    /// Steps between population-data writes (0 = never).
    pub wsteps: i32,
    /// Steps between MSI movie frames (0 = never; `-1` in the file is
    /// rewritten to `nsteps + 1`, i.e. "final frame only").
    pub msteps: i32,
    /// Draw bonds in MSI output? Kept as `i32`, not `bool`, because the C++
    /// field is an `int` that transiently holds the *seed* value mid-parse
    /// (see the WART below) — faithfulness beats prettiness here.
    pub drawbonds: i32,
    /// PRNG seed. **Always 0** — see the WART note in [`read_sim`]. The C++
    /// member is `long`; on the 64-bit Linux the golden run used, that is
    /// 64 bits, hence `i64`.
    pub ranseed: i64,
}

/// Read `data.sim`, reproducing the C++ read sequence **exactly** —
/// including its most consequential bug.
///
/// # WART (spec B2): the doubled `drawbonds` read swallows the seed
///
/// `Simulation::CreateSimulation` (sim.cpp lines 14–29) reads, in order:
/// `nsteps`, `wsteps`, `msteps`, then **`drawbonds` twice**. The first
/// `drawbonds` read lands on the file's *seed* line and discards it (the
/// golden `data.sim` has `-2` there); the second lands on the real
/// drawbonds value. `ranseed` is never read from the file at all: the C++
/// heap-allocates with `new Simulation()` (value-initialization), so
/// `ranseed` is zero — and `ran2` treats a seed of 0 as "seed yourself",
/// mapping it to 1. **Net effect: every run of the legacy model ever made
/// used the same PRNG stream**, which is precisely why the golden reference
/// capture is deterministic and why this port can chase bitwise parity.
///
/// Port disposition: preserved AS-IS by task decree (TASK-016 overrides the
/// design doc's day-one fix). A correct 5-field reader is a later, flagged
/// change; until then this function is bug-compatible and `ranseed` is
/// always 0. Full analysis: `projects/kmc/02-model-spec.md` Part B, §B2.
pub fn read_sim(path: &Path) -> Result<SimParams, ReadError> {
    let mut s = Scanner::open(path)?;
    parse_sim(&mut s)
}

/// The parse itself, split from file-opening so tests can feed it strings.
///
/// \[IDIOM\] `&mut Scanner` — a mutable *borrow*. The caller keeps ownership
/// of the scanner; we get temporary, exclusive permission to advance it.
/// Exclusivity is the guarantee C++ references never gave: while this
/// function holds the `&mut`, the compiler proves nobody else can even
/// *look* at the scanner, so "who else is moving my cursor?" is not a
/// question that can arise.
fn parse_sim(s: &mut Scanner) -> Result<SimParams, ReadError> {
    // Line-by-line transliteration of sim.cpp — same reads, same order,
    // same die-messages.
    let nsteps = s.next_i32("invalid number of steps in input file")?;
    s.eat_comment();
    let wsteps = s.next_i32("invalid number of data write steps in input file")?;
    s.eat_comment();
    let mut msteps = s.next_i32("invalid number of movie write steps in input file")?;
    if msteps == -1 {
        msteps = nsteps + 1; // "-1" means: only the final movie frame
    }
    s.eat_comment();
    // WART (spec B2), preserved: this read consumes the SEED line into
    // drawbonds and throws it away...
    //
    // [IDIOM] The `#[allow]` below is load-bearing pedagogy: rustc flags
    // this line with `unused_assignments` — the compiler FOUND the
    // 25-year-old bug (a value read and immediately overwritten) that
    // g++ -Wall -Wextra never mentioned. We must silence the lint because
    // we are preserving the bug on purpose; the annotation is the receipt.
    // This is the "strict compiler as a defense line" story in one line of
    // code: agent- or human-written, dead stores don't pass unnoticed.
    #[allow(unused_assignments)]
    let mut drawbonds = s.next_i32("invalid number of steps in input file")?;
    s.eat_comment();
    // ...and this one reads the actual drawbonds value over it.
    drawbonds = s.next_i32("invalid draw bonds? parameter in input file")?;

    Ok(SimParams {
        nsteps,
        wsteps,
        msteps,
        drawbonds,
        // WART (spec B2), preserved: never read from the file. The C++
        // `new Simulation()` zero-initializes the member; we write the zero
        // explicitly, which is the only honest way to say it in Rust —
        // there is no "forgot to initialize" state to reproduce, because
        // the language has no uninitialized reads to offer.
        ranseed: 0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The golden fixture, byte-for-byte (data/golden/inputs/data.sim).
    const GOLDEN: &str = "\
20000     # Number of steps in simulation
1000        # Number of steps between data writes (0 for no data)
1000000     # Number of steps between movie frames (0 for no movie)
-2\t    # Seed for random number generator
1           # Draw bonds? (0 No - only occupied, 1 Yes - entire lattice)
";

    #[test]
    fn golden_data_sim_parses_like_the_cpp() {
        let mut s = Scanner::from_text(GOLDEN, "data.sim");
        let p = parse_sim(&mut s).unwrap();
        assert_eq!(p.nsteps, 20000);
        assert_eq!(p.wsteps, 1000);
        assert_eq!(p.msteps, 1_000_000);
        assert_eq!(p.drawbonds, 1);
        // The wart, pinned as a test: the file says seed = -2, the parse
        // says 0. If someone "fixes" B2 without meaning to, this fails.
        assert_eq!(p.ranseed, 0);
    }

    #[test]
    fn golden_fixture_file_matches_the_inline_copy() {
        // Belt and suspenders: the checked-in fixture parses identically.
        let root = concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/golden/inputs/data.sim");
        let p = read_sim(Path::new(root)).unwrap();
        assert_eq!(
            p,
            SimParams { nsteps: 20000, wsteps: 1000, msteps: 1_000_000, drawbonds: 1, ranseed: 0 }
        );
    }

    #[test]
    fn seed_line_is_swallowed_even_with_distinct_values() {
        // Synthetic file where seed and drawbonds differ, to prove which
        // token wins: drawbonds must be the 5th value, seed discarded.
        let text = "10 # a\n2 # b\n3 # c\n777 # seed\n0 # drawbonds\n";
        let mut s = Scanner::from_text(text, "data.sim");
        let p = parse_sim(&mut s).unwrap();
        assert_eq!(p.drawbonds, 0); // not 777
        assert_eq!(p.ranseed, 0); // never 777 either
    }

    #[test]
    fn msteps_minus_one_becomes_nsteps_plus_one() {
        let text = "500 #\n0 #\n-1 #\n-2 #\n1 #\n";
        let mut s = Scanner::from_text(text, "data.sim");
        let p = parse_sim(&mut s).unwrap();
        assert_eq!(p.msteps, 501);
    }

    #[test]
    fn truncated_file_reports_the_cpp_die_message() {
        let mut s = Scanner::from_text("10 # only nsteps\n", "data.sim");
        let err = parse_sim(&mut s).unwrap_err();
        assert!(err.to_string().contains("invalid number of data write steps"));
    }
}
