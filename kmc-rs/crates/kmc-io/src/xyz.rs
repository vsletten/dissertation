//! XYZ writer — port of `output::writeXYZ` (output.cpp), bug-compatible.
//!
//! "XYZ" generously: the C++ emits a **nonstandard dialect** (spec A9.2) —
//! comma-separated `Z,x,y,z` rows, no atom-count header, no comment line,
//! whitespace nowhere. Standard XYZ tooling will not read it; the golden
//! gate compares it byte-for-byte, so the dialect is preserved verbatim.
//! Al→13, Si→14, every oxygen class→8 (no H/F rendering games here,
//! unlike the MSI writer).
//!
//! Two warts live at the *call site* (the CLI), not here, but matter for
//! anyone hunting for this file's output:
//! * The C++ hard-codes the output name **`start.xyz`** and calls the
//!   writer once, at shutdown — so `start.xyz` is an END-state artifact
//!   under a misleading name, and no `end.xyz` ever exists (TASK-004
//!   confirmed; the golden manifest says so explicitly).
//! * Site selection here is by explicit state *ranges* rather than the MSI
//!   writer's `%100 > 0 && != EDGE` — the two encode the same "occupied,
//!   non-EDGE" set for every state the model can produce, but the C++
//!   wrote the test twice in two styles, and the port keeps each writer's
//!   own style so line-by-line review against output.cpp stays 1:1.

use std::io::{self, Write};

use crate::fmt::format_g6;
use crate::view::{LatticeView, basis_lengths, site_xyz};

/// Write every occupied site as a `Z,x,y,z` row (the legacy dialect).
pub fn write_xyz<W: Write>(w: &mut W, view: &LatticeView<'_>) -> io::Result<()> {
    let (al, bl, cl) = basis_lengths();

    for i in 0..view.graph.len() {
        let s = view.graph.sites[i].state.0;
        // The C++'s five range tests, verbatim: occupied Si (201..=205),
        // occupied Al (101..=107), occupied O of each class (x01..x99).
        // [IDIOM] `matches!` — one expression, five patterns, no repeated
        // `lattice->sites[idxSite].state` spelling. Range patterns are
        // *inclusive* (`..=`), so the C++'s `> 200 && < 206` becomes
        // `201..=205` — same set, bounds visible at a glance.
        let selected = matches!(s, 201..=205 | 101..=107 | 301..=399 | 401..=499 | 501..=599);
        if !selected {
            continue;
        }
        let (x, y, z) = site_xyz(view, i, al, bl, cl);
        let zed = match s / 100 {
            1 => 13, // Al
            2 => 14, // Si
            _ => 8,  // all oxygen classes
        };
        writeln!(
            w,
            "{},{},{},{}",
            zed,
            format_g6(x),
            format_g6(y),
            format_g6(z)
        )?;
    }
    Ok(())
}
