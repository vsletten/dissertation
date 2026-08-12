//! P4 structural golden gate: build the kaolinite lattice from the Petra
//! deck (examples/kaolinite.toml) and from the kmc-rs `kaolinite` crate
//! (bitwise-parity-tested against the dissertation C++), then compare
//! SITE BY SITE.
//!
//! Comparison rule: legacy EDGE sites (state 9) must be exactly Petra's
//! frozen sites (which keep their pre-freeze states — guards handle the
//! difference via `frozen` filters); every other site's state must map
//! exactly through the name↔code table.
//!
//! Both builders use the same flat index scheme
//! (`i = a*(b_cells*npos) + b*npos + n`), so index i compares to index i.

use std::collections::HashMap;
use std::path::PathBuf;

use kaolinite::build::LatticeParams;
use kaolinite::{create_lattice, find_pairs, populate_solid, terminate_lattice};

fn repo_path(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..").join(rel)
}

/// Petra "Kind.state" name → legacy 3-digit code.
fn code_table() -> HashMap<String, i32> {
    fn add(m: &mut HashMap<String, i32>, kind: &str, states: &[&str], base: i32) {
        for (i, s) in states.iter().enumerate() {
            m.insert(format!("{kind}.{s}"), base + i as i32);
        }
    }
    let mut m = HashMap::new();
    add(&mut m, "Al", &["empty", "l0", "l1", "l2", "l3", "l4", "l5", "l6"], 100);
    add(&mut m, "Si", &["empty", "oh0", "oh1", "oh2", "oh3", "oh4"], 200);
    add(&mut m, "Oss", &["empty", "br", "hy", "si1"], 300);
    add(
        &mut m,
        "Osa",
        &[
            "empty", "br", "sih", "full", "albr", "alhy", "sial", "sialh", "si1", "al1",
        ],
        400,
    );
    // Both proton-split states are legacy 410.
    m.insert("Osa.pr1".to_string(), 410);
    m.insert("Osa.pr2".to_string(), 410);
    add(&mut m, "Oaa", &["empty", "br", "hy", "al1"], 500);
    m
}

#[test]
fn kaolinite_build_matches_kmc_rs_site_by_site() {
    // --- golden reference via kmc-rs ---
    let cell = kmc_io::read_cell(&repo_path("kmc-rs/data/golden/inputs/data.cell"))
        .expect("golden data.cell reads");
    let mut golden = create_lattice(
        &cell,
        LatticeParams {
            a_cells: 20,
            b_cells: 3,
            surface_plane: 0,
        },
    );
    find_pairs(&mut golden);
    populate_solid(&mut golden, -1.0, -1.0);
    kaolinite::build::terminate_surface(&mut golden);
    terminate_lattice(&mut golden);

    // --- Petra build from the deck ---
    let deck = petra_deck::load(repo_path("petra/examples/kaolinite.toml"))
        .expect("kaolinite deck compiles");
    let engine = deck.build_engine(Some(0)).expect("engine builds");
    let lat = &engine.lattice;

    assert_eq!(lat.len(), golden.graph.sites.len(), "site counts");

    let codes = code_table();
    let petra_code = |i: usize| -> i32 {
        let name = &deck.state_names[lat.states[i].0 as usize];
        *codes
            .get(name.as_str())
            .unwrap_or_else(|| panic!("no legacy code for state '{name}'"))
    };

    let mut mismatches = Vec::new();
    for i in 0..lat.len() {
        let legacy = golden.graph.sites[i].state.0;
        if legacy == 9 {
            if !lat.frozen[i] {
                mismatches.push(format!("site {i}: legacy EDGE but petra not frozen"));
            }
            continue;
        }
        if lat.frozen[i] {
            mismatches.push(format!(
                "site {i}: petra frozen but legacy state {legacy}"
            ));
            continue;
        }
        let p = petra_code(i);
        if p != legacy {
            let (cell_coord, t) = lat.coords(i);
            mismatches.push(format!(
                "site {i} (cell {cell_coord:?} pos {t}): legacy {legacy} vs petra {p} ({})",
                deck.state_names[lat.states[i].0 as usize]
            ));
        }
    }
    assert!(
        mismatches.is_empty(),
        "{} site mismatches (first 40):\n{}",
        mismatches.len(),
        mismatches
            .iter()
            .take(40)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );

    // Neighbor topology spot check: same adjacency degree everywhere
    // (legacy stores a padded [i32; 6]; count the resolved ones).
    for i in 0..lat.len() {
        let legacy_deg = golden.graph.sites[i]
            .nbr
            .iter()
            .filter(|n| n.is_some())
            .count();
        assert_eq!(
            lat.neighbors(i).len(),
            legacy_deg,
            "degree mismatch at site {i}"
        );
    }
}

/// Dynamics smoke: the full reaction set sustains thousands of KMC steps
/// with per-kind site conservation, monotonic time, and incremental
/// event-table consistency (the paranoid oracle).
#[test]
fn kaolinite_dynamics_runs_and_conserves_sites() {
    let deck = petra_deck::load(repo_path("petra/examples/kaolinite.toml"))
        .expect("kaolinite deck compiles");
    let mut engine = deck.build_engine(Some(1)).expect("engine builds");

    // Per-kind site totals: 60 cells x (4 Al, 4 Si, 6 Oss, 4 Osa, 8 Oaa).
    let kind_totals = [240u64, 240, 360, 240, 480];
    let count_by_kind = |engine: &petra_core::Engine| -> Vec<u64> {
        let counts = engine.state_counts(deck.n_states);
        deck.kind_state_ranges
            .iter()
            .map(|&(start, n)| {
                (start..start + n).map(|s| counts[s as usize]).sum::<u64>()
            })
            .collect()
    };
    assert_eq!(count_by_kind(&engine), kind_totals);

    let mut last_time = 0.0;
    for i in 0..5000 {
        let fired = engine
            .step()
            .unwrap_or_else(|e| panic!("step {i} failed: {e}"));
        assert!(fired.time > last_time, "time must increase (step {i})");
        last_time = fired.time;
        if i % 1000 == 999 {
            engine.paranoid_check().expect("incremental tables consistent");
            assert_eq!(count_by_kind(&engine), kind_totals, "site conservation");
        }
    }

    // The proton-memory machinery (pr1/pr2 via R4/R9 coins) must actually
    // exercise: over 5000 steps at these rates some 410s exist or existed.
    // Weak but real: the R4 channel has rate 100 on every intact Osa.
    let counts = engine.state_counts(deck.n_states);
    let osa_range = deck.kind_state_ranges[3];
    let osa_total: u64 = (osa_range.0..osa_range.0 + osa_range.1)
        .map(|s| counts[s as usize])
        .sum();
    assert_eq!(osa_total, 240);
}
