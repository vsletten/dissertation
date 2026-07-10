//! **The M3 golden gate**: the deterministic structural build must
//! reproduce the C++ initial state **bitwise**.
//!
//! Reference: `data/golden/` — TASK-004's capture of the legacy binary
//! (g++ 13.3.0, `-std=c++11 -O3 -ffast-math`), SHA-256-manifested in
//! mission-control `projects/kmc/golden/manifest.md`, reproducibility
//! confirmed by duplicate run. `start.msi` is written by the C++ *before*
//! its first MC step (mckaol.cpp line 58), so it is a pure rendering of the
//! structural build: lattice tiling → find_pairs → populate_solid →
//! terminate_surface → terminate_lattice, plus the (bug-compatible,
//! spec B5) coordinate transform and `%g` formatting. No RNG is involved,
//! which is why *exact* equality is a fair demand (spec §C1).
//!
//! One caution from the M0 session plan (§4): the golden `start.xyz` is NOT
//! an initial-state artifact — the C++ writes it at shutdown (after 20,000
//! steps) under a misleading name. The initial state has exactly one
//! rendering: `start.msi`. That is what we diff.
//!
//! \[IDIOM\] Integration tests live in `tests/` beside `src/`, compiled as
//! *external* crates: they can only touch the public API, exactly like a
//! downstream user. Unit tests (in `#[cfg(test)]` modules) see internals;
//! integration tests certify the contract. This gate belongs out here.

use std::path::Path;

use kaolinite::{create_lattice, find_pairs, populate_solid, terminate_lattice, terminate_surface};

fn golden_dir() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/golden"))
}

/// Run the full M3 pipeline on the golden inputs and return the MSI bytes.
fn build_start_msi() -> Vec<u8> {
    let inputs = golden_dir().join("inputs");
    let sim = kmc_io::read_sim(&inputs.join("data.sim")).expect("golden data.sim reads");
    let rxns = kmc_io::read_rxn(&inputs.join("data.rxn")).expect("golden data.rxn reads");
    let cell = kmc_io::read_cell(&inputs.join("data.cell")).expect("golden data.cell reads");
    let params = kmc_io::read_lattice(&inputs.join("data.lattice")).expect("golden data.lattice");

    let mut structure = create_lattice(&cell, params);
    find_pairs(&mut structure);
    populate_solid(&mut structure, rxns.dm_si, rxns.dm_al);
    terminate_surface(&mut structure);
    terminate_lattice(&mut structure);

    let mut buf = Vec::new();
    kmc_io::write_msi(&mut buf, &structure, &cell, "start", sim.drawbonds)
        .expect("writing to a Vec cannot fail");
    buf
}

/// THE gate. Byte-for-byte equality with the C++ output — floats, IDs,
/// bond order, line endings, everything.
#[test]
fn start_msi_is_bitwise_identical_to_the_cpp_golden() {
    let got = build_start_msi();
    let want =
        std::fs::read(golden_dir().join("outputs/start.msi")).expect("golden start.msi exists");

    if got != want {
        // Byte equality failed: report the first differing LINE so the
        // failure is debuggable without a hex editor.
        let got_s = String::from_utf8_lossy(&got);
        let want_s = String::from_utf8_lossy(&want);
        for (n, (g, w)) in got_s.lines().zip(want_s.lines()).enumerate() {
            assert_eq!(g, w, "first divergence at line {} (1-based {})", n, n + 1);
        }
        // Same common prefix but different lengths:
        panic!(
            "outputs share {} lines but differ in length: got {} bytes, want {} bytes",
            got_s.lines().count().min(want_s.lines().count()),
            got.len(),
            want.len()
        );
    }
}

/// Cross-checks that make gate failures diagnosable and pin the headline
/// numbers independently of the byte diff: 1,000 atoms, 1,426 bonds, and
/// the C++'s final object id (1 model + atoms + bonds).
#[test]
fn golden_msi_shape_counts() {
    let got = build_start_msi();
    let text = String::from_utf8(got).expect("MSI output is ASCII");
    let atoms = text.lines().filter(|l| l.ends_with(" Atom")).count();
    let bonds = text.lines().filter(|l| l.ends_with(" Bond")).count();
    assert_eq!(atoms, 1000);
    assert_eq!(bonds, 1426);
    assert!(text.contains("(2427 Bond"), "last object id should be 2427");
}

/// The build is deterministic end to end (the premise of bitwise gating):
/// run the whole pipeline twice, byte-identical both times.
#[test]
fn pipeline_is_deterministic() {
    assert_eq!(build_start_msi(), build_start_msi());
}
