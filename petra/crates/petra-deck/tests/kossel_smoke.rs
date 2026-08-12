//! End-to-end smoke test on the tutorial deck (design doc P0 exit test):
//! compile examples/kossel.toml, run real KMC steps, and check physical
//! sanity — time advances, populations are conserved, undersaturated
//! solution dissolves the crystal, and the incremental event tables agree
//! with a from-scratch re-derivation (the paranoid oracle).

use std::path::PathBuf;

fn deck_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/kossel.toml")
}

#[test]
fn kossel_deck_compiles_with_expected_shape() {
    let deck = petra_deck::load(deck_path()).expect("deck should load");
    assert_eq!(deck.kind_names, ["M_site"]);
    assert_eq!(deck.state_names, ["M_site.occupied", "M_site.empty"]);
    assert_eq!(deck.reactions.len(), 2);
    assert_eq!(deck.unit_cell.sites[0].bonds.len(), 6, "expanded to ±a ±b ±c");
    // Attach folds mu = -2 into ln_thermo: ln(1) + (-2)/RT.
    let attach = deck
        .reactions
        .iter()
        .find(|r| r.name == "attach")
        .expect("attach reaction");
    let rt = petra_core::rate::R_KCAL * deck.temperature;
    assert!((attach.ln_thermo - (-2.0 / rt)).abs() < 1e-12);
}

#[test]
fn kossel_run_dissolves_under_undersaturation() {
    let deck = petra_deck::load(deck_path()).expect("deck should load");
    let mut engine = deck.build_engine(Some(7));
    let n_sites = engine.lattice.len() as u64;
    assert_eq!(n_sites, 8 * 8 * 8);

    let occupied_start = engine.state_counts(deck.n_states)[0];
    assert_eq!(occupied_start, n_sites, "perfect-crystal fill");

    let mut last_time = 0.0;
    let steps = 5000;
    let mut ran = 0;
    for i in 0..steps {
        match engine.step() {
            Ok(fired) => {
                assert!(
                    fired.time > last_time,
                    "time must strictly increase (step {i})"
                );
                last_time = fired.time;
                ran += 1;
            }
            Err(petra_core::Stop::NoEvents) => break, // fully dissolved
            Err(e) => panic!("unexpected stop at step {i}: {e}"),
        }
        if i % 1000 == 999 {
            engine
                .paranoid_check()
                .expect("incremental tables must match from-scratch derivation");
        }
    }
    assert!(ran > 1000, "should sustain KMC steps, ran only {ran}");

    let counts = engine.state_counts(deck.n_states);
    assert_eq!(
        counts.iter().sum::<u64>(),
        n_sites,
        "sites are conserved: every site is exactly one state"
    );
    let occupied_end = counts[0];
    assert!(
        occupied_end < occupied_start,
        "mu = -2 kcal/mol must dissolve: {occupied_start} -> {occupied_end}"
    );
}

#[test]
fn same_seed_reproduces_same_seed_diverges_different() {
    let deck = petra_deck::load(deck_path()).expect("deck should load");

    let run = |seed: u64| {
        let mut e = deck.build_engine(Some(seed));
        for _ in 0..500 {
            if e.step().is_err() {
                break;
            }
        }
        (e.time, e.state_counts(deck.n_states))
    };

    let (t1, c1) = run(11);
    let (t2, c2) = run(11);
    let (t3, c3) = run(12);
    assert_eq!(t1.to_bits(), t2.to_bits(), "same seed → identical trajectory");
    assert_eq!(c1, c2);
    assert!(
        t1.to_bits() != t3.to_bits() || c1 != c3,
        "different seed → different trajectory"
    );
}
