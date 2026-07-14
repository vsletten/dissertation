//! Population time-series row — port of `output::writeData` (output.cpp),
//! bug-compatible.
//!
//! One CSV row per call (spec A9.3):
//! `time, Si_total, Si(OH)0..4, Al_total, Al(OH,H2O)0..6` — fourteen
//! comma-separated fields, `time` in the C++'s `%g` float rendering, the
//! thirteen counts as plain integers.
//!
//! **WART (spec B1), preserved by decree**: the C++ never appends. Each
//! snapshot goes to its own file (`step{i}.dat`, `end.dat`), truncated and
//! written with this single row — a 20k-step golden run leaves 20 one-row
//! files, a 5M-step production run leaves *thousands*. The documented
//! `results.dat` appended series is never written (only `remove`d at
//! startup — see the CLI). This module writes **one row to one sink** and
//! the caller owns the scattering, exactly like the C++; the one-appended-
//! file repair is designed in docs/REFORM_PLAN.md (R8) and lands with the
//! reform era, not here — matching the golden capture byte-for-byte is
//! M7's whole job.
//!
//! Counting quirks kept faithfully:
//! * Only cation classes are counted; oxygens are invisible to this output.
//! * The histogram index is `state - 101` / `state - 201` — the C++ would
//!   index out of bounds for the never-produced `WRONG` markers (199/299);
//!   Rust would panic instead. Both are "can't happen" per spec A3; the
//!   golden gate certifies they don't.

use std::io::{self, Write};

use kaolinite::state::State;
use kmc_engine::SiteGraph;

use crate::fmt::format_g6;

/// Write one population row for the current lattice state.
///
/// `t` is the simulated-time accumulator — an `f32` because the C++'s is
/// (spec B8; the f64 clock is reform R6). Takes the bare graph, not a
/// [`crate::view::LatticeView`]: populations are pure state counting, and
/// taking only the graph *proves* geometry can't leak in (contrast
/// `writeData(const char*, Lattice*, int, float)`, which hands the writer
/// the whole world and lets you wonder).
pub fn write_data<W: Write>(w: &mut W, graph: &SiteGraph<State>, t: f32) -> io::Result<()> {
    // nsi[oh] = count of Si with `oh` OH groups; nal likewise for Al.
    let mut nsi = [0i32; 5];
    let mut nal = [0i32; 7];
    let (mut sitot, mut altot) = (0i32, 0i32);

    for site in &graph.sites {
        let s = site.state.0;
        if s % 100 == 0 {
            continue; // skip empty sites (EDGE=9 survives this, but 9/100=0 matches no class below)
        }
        match s / 100 {
            1 => {
                nal[(s - 101) as usize] += 1;
                altot += 1;
            }
            2 => {
                nsi[(s - 201) as usize] += 1;
                sitot += 1;
            }
            _ => {}
        }
    }

    // Field order is the C++'s exactly: t, Si total, the five Si bins,
    // Al total, the seven Al bins.
    write!(w, "{}", format_g6(t))?;
    write!(w, ",{sitot}")?;
    for n in nsi {
        write!(w, ",{n}")?;
    }
    write!(w, ",{altot}")?;
    for n in nal {
        write!(w, ",{n}")?;
    }
    writeln!(w)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use kmc_engine::Site;

    fn graph_of(states: &[i32]) -> SiteGraph<State> {
        SiteGraph {
            sites: states
                .iter()
                .map(|&s| Site {
                    state: State(s),
                    nbr: [None; 6],
                })
                .collect(),
        }
    }

    #[test]
    fn counts_cations_into_the_right_bins() {
        // 2× Si(OH)0 (201), 1× Si(OH)4 (205), 1× Al(·)0 (101), 1× Al(·)6
        // (107); oxygens and empties invisible; EDGE ignored.
        let g = graph_of(&[201, 201, 205, 101, 107, 301, 404, 100, 200, 9]);
        let mut out = Vec::new();
        write_data(&mut out, &g, 0.5).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "0.5,3,2,0,0,0,1,2,1,0,0,0,0,0,1\n"
        );
    }

    #[test]
    fn time_uses_the_cpp_g_rendering() {
        // The golden step0.dat leads with "2.71473e-05" — %g scientific.
        let g = graph_of(&[]);
        let mut out = Vec::new();
        write_data(&mut out, &g, 2.71473e-05).unwrap();
        assert!(String::from_utf8(out).unwrap().starts_with("2.71473e-05,"));
    }
}
