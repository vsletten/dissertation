//! **The M6 dynamics parity gate**: the Rust simulation reproduces the C++
//! trajectory **step-by-step, bitwise**, under the legacy fixed seed.
//!
//! # The oracle
//!
//! `data/golden/trajectory/traj_seed0_full.txt` is a per-step trace captured
//! from the read-only C++ model (`dissertation/main/model`, rebuilt verbatim
//! in scratch — no source edits to the model logic, only an added trace
//! `fprintf`; the same build reproduces the golden `start.msi`/`end.msi`
//! hashes, proving the instrumentation is inert). Each line is:
//!
//! ```text
//! step  site  rxn  0xDT_BITS  occupied  0xSTATE_HASH
//! ```
//!
//! * `site`, `rxn` — the event the C++ selected and applied that step;
//! * `DT_BITS` — the raw 32-bit pattern of the C++ `float dt` (bitwise, so an
//!   ulp of drift fails loudly rather than hiding under an epsilon);
//! * `occupied` — count of occupied non-EDGE sites after the step;
//! * `STATE_HASH` — FNV-1a-64 over every site's 4-byte state code, in site
//!   order — a whole-lattice fingerprint. Two runs agreeing on this every
//!   step **are** in the same state (spec §C2a "step-by-step state
//!   agreement").
//!
//! # Why bitwise is attainable here (and the honest caveats)
//!
//! Three legacy facts conspire to make exact reproduction possible, not
//! merely statistical: the seed is fixed (spec B2 — every legacy run shares
//! one `ran2` stream), the arithmetic is `f32` throughout (kept `f32` in the
//! port for exactly this, spec §C2a), and the event-list summation order is
//! deterministic (reproduced by the reversed fold in `kmc_engine::step`).
//! The one thing that had to be earned rather than assumed: the *order* of
//! the `f32` `ratesum`/`partsum` adds. Get it wrong and a single ulp flips
//! which event crosses `eps`, and the trajectories diverge irrecoverably a
//! few steps later. This test is the proof that the order is right — for all
//! 20,000 steps of the reference run.

use std::path::Path;

use kaolinite::{
    Kaolinite, State, create_lattice, find_pairs, populate_solid, terminate_lattice,
    terminate_surface,
};
use kmc_engine::{Ran2, step};

fn repo_root() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../.."))
}

/// FNV-1a-64 over each site's state code (4 little-endian bytes), in site
/// index order — byte-identical to the C++ capture harness.
///
/// Note the offset basis: the capture harness used `1469598103934665603`, a
/// digit-truncated form of the canonical FNV-1a-64 basis
/// (`14695981039346656037`). The hash is only a whole-lattice *fingerprint* —
/// any fixed basis fingerprints equally well — so what matters is that this
/// matches the oracle bit-for-bit, not that it is the textbook constant. We
/// deliberately reproduce the harness's basis rather than "correct" it.
fn state_hash(graph: &kmc_engine::SiteGraph<State>) -> u64 {
    let mut h: u64 = 1469598103934665603; // the harness's (non-canonical) basis
    for site in &graph.sites {
        for b in site.state.0.to_le_bytes() {
            h ^= b as u64;
            h = h.wrapping_mul(0x100000001b3); // FNV prime
        }
    }
    h
}

/// Count occupied, non-EDGE sites (C++ trace's `occ`: `state%100>0 && !EDGE`).
fn occupied(graph: &kmc_engine::SiteGraph<State>) -> i64 {
    graph
        .sites
        .iter()
        .filter(|s| s.state.is_occupied() && !s.state.is_edge())
        .count() as i64
}

/// One expected trajectory row.
struct Row {
    step: usize,
    site: usize,
    rxn: u16,
    dt_bits: u32,
    occ: i64,
    hash: u64,
}

fn parse_oracle(text: &str) -> Vec<Row> {
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split_whitespace().collect();
            Row {
                step: f[0].parse().unwrap(),
                site: f[1].parse().unwrap(),
                rxn: f[2].parse().unwrap(),
                dt_bits: u32::from_str_radix(f[3].trim_start_matches("0x"), 16).unwrap(),
                occ: f[4].parse().unwrap(),
                hash: u64::from_str_radix(f[5].trim_start_matches("0x"), 16).unwrap(),
            }
        })
        .collect()
}

/// Build the initial state from the golden inputs (the M3 structural build),
/// then hand off to the dynamics model.
fn build_simulation() -> (kmc_engine::SiteGraph<State>, Kaolinite) {
    let inputs = repo_root().join("data/golden/inputs");
    let rxns = kmc_io::read_rxn(&inputs.join("data.rxn")).expect("data.rxn");
    let cell = kmc_io::read_cell(&inputs.join("data.cell")).expect("data.cell");
    let params = kmc_io::read_lattice(&inputs.join("data.lattice")).expect("data.lattice");

    let mut structure = create_lattice(&cell, params);
    find_pairs(&mut structure);
    populate_solid(&mut structure, rxns.dm_si, rxns.dm_al);
    terminate_surface(&mut structure);
    terminate_lattice(&mut structure);

    Kaolinite::from_structure(structure, rxns)
}

/// THE gate. Run the legacy-seeded simulation and check every step against
/// the C++ oracle: chosen (site, rxn), the `dt` bit pattern, the occupied
/// count, and the whole-lattice state hash.
#[test]
fn trajectory_is_bitwise_identical_to_the_cpp_under_the_legacy_seed() {
    let oracle_text =
        std::fs::read_to_string(repo_root().join("data/golden/trajectory/traj_seed0_full.txt"))
            .expect("trajectory oracle present");
    let oracle = parse_oracle(&oracle_text);
    assert_eq!(
        oracle.len(),
        20_000,
        "oracle should hold the full 20k-step run"
    );

    let (mut graph, mut model) = build_simulation();
    let mut rng = Ran2::legacy(); // spec B2 fixed seed
    let mut scratch = Vec::new();

    for (i, row) in oracle.iter().enumerate() {
        let adv = step(&mut graph, &mut model, &mut rng, &mut scratch)
            .unwrap_or_else(|e| panic!("step {i} failed: {e:?}"));

        assert_eq!(row.step, i, "oracle step index out of order at {i}");
        assert_eq!(
            (adv.event.site, adv.event.rxn),
            (row.site, row.rxn),
            "step {i}: chose (site {}, rxn {}), C++ chose (site {}, rxn {})",
            adv.event.site,
            adv.event.rxn,
            row.site,
            row.rxn
        );
        assert_eq!(
            adv.dt.to_bits(),
            row.dt_bits,
            "step {i}: dt bits diverge — got 0x{:08x} ({}), want 0x{:08x}",
            adv.dt.to_bits(),
            adv.dt,
            row.dt_bits
        );
        assert_eq!(
            occupied(&graph),
            row.occ,
            "step {i}: occupied count diverges"
        );
        assert_eq!(
            state_hash(&graph),
            row.hash,
            "step {i}: whole-lattice state hash diverges"
        );
    }
}

/// A cheaper smoke test that the first steps exercise several reaction
/// classes (so the gate above is not trivially satisfied by one reaction) —
/// documents the coverage the parity claim rests on.
#[test]
fn early_trajectory_exercises_multiple_reaction_classes() {
    let (mut graph, mut model) = build_simulation();
    let mut rng = Ran2::legacy();
    let mut scratch = Vec::new();

    let mut seen = std::collections::BTreeSet::new();
    for _ in 0..2000 {
        let adv = step(&mut graph, &mut model, &mut rng, &mut scratch).expect("step");
        seen.insert(adv.event.rxn);
    }
    // The first 2000 steps hit hydrolysis (both directions), the R4 proton
    // coin, adsorption and desorption — a broad spread, not one channel.
    for rxn in [0u16, 2, 3, 4, 5, 6, 7, 14, 15, 16, 19, 20, 22] {
        assert!(
            seen.contains(&rxn),
            "expected reaction {rxn} within 2000 steps"
        );
    }
}
