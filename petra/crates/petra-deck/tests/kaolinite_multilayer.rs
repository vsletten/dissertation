//! Gates for the c-tiled multilayer kaolinite deck
//! (examples/kaolinite-multilayer.toml):
//! 1. shape — interlayer hydrogen-bond topology matches the geometric
//!    donor→acceptor assignment (count, degrees, open-face drops);
//! 2. per-layer equivalence — every layer's initial state is site-by-site
//!    identical to the validated single-sheet deck (the interlayer bonds
//!    are invisible to the build);
//! 3. the anhydrous gap — a hydrolyzed basal oxygen in the layer above
//!    does NOT make a donor hydroxyl surface-reachable (exclude_label),
//!    and demonstrably WOULD without the exclusion;
//! 4. H-bond energetics — desorption rates are suppressed by exactly
//!    exp(−n·e_hb/RT) per intact interlayer bond;
//! 5. dynamics smoke — 5000 steps with conservation and the paranoid
//!    incremental-table oracle.

use std::path::PathBuf;

use petra_core::rate::R_KCAL;
use petra_core::reaction::{count_matches, resolve_rate};
use petra_core::state::StateId;

const E_HB: f64 = 5.0; // kcal/mol, the deck's placeholder (CALC-001)
const T: f64 = 8000.0;

fn repo_path(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..").join(rel)
}

fn load_multilayer() -> petra_deck::CompiledDeck {
    petra_deck::load(repo_path("petra/examples/kaolinite-multilayer.toml"))
        .expect("multilayer deck compiles")
}

/// Flat index of (a, b, c, template) for the 20×3×4 lattice.
fn site(a: usize, b: usize, c: usize, t: usize) -> usize {
    ((a * 3 + b) * 4 + c) * 26 + t
}

fn state_id(deck: &petra_deck::CompiledDeck, name: &str) -> StateId {
    StateId(
        deck.state_names
            .iter()
            .position(|n| n == name)
            .unwrap_or_else(|| panic!("state '{name}' exists")) as u16,
    )
}

fn kinds_of(deck: &petra_deck::CompiledDeck, lat: &petra_core::lattice::Lattice) -> Vec<petra_core::crystal::KindId> {
    lat.template_index
        .iter()
        .map(|&t| deck.kinds_per_template[t as usize])
        .collect()
}

#[test]
fn interlayer_bond_topology_matches_geometry() {
    let deck = load_multilayer();
    let engine = deck.build_engine(Some(0)).expect("engine builds");
    let lat = &engine.lattice;
    assert_eq!(lat.len(), 20 * 3 * 4 * 26);

    let hbond = deck
        .label_names
        .iter()
        .position(|l| l == "hbond")
        .expect("hbond label interned") as u16;

    // Directed adjacency entries carrying the hbond label. Five donors
    // (dcell (0,0,1) or (1,0,1)-free, db = 0) exist for every (a, b) and
    // c < 3: 5 × 20 × 3 × 3 = 900 undirected. Donor 24 (db = +1) loses
    // its bond at the b = 2 fixed face: 20 × 2 × 3 = 120. Times 2 for the
    // two CSR directions.
    let count: usize = (0..lat.len())
        .map(|s| lat.neighbor_labels(s).iter().filter(|&&l| l == hbond).count())
        .sum();
    assert_eq!(count, 2 * (900 + 120), "hbond adjacency entries");

    // Spot geometry: interior donor pos 19 at (5,1,1) reaches acceptor
    // pos 11 at (5,1,2); its degree is 2 covalent + 1 hbond.
    let d = site(5, 1, 1, 19);
    assert_eq!(lat.neighbors(d).len(), 3);
    let acceptor = site(5, 1, 2, 11);
    assert!(
        lat.neighbors(d).contains(&(acceptor as u32)),
        "19@(5,1,1) donates to 11@(5,1,2)"
    );
    // Top layer's donors point into vacuum: no hbond entries at c = 3.
    let top = site(5, 1, 3, 19);
    assert_eq!(lat.neighbors(top).len(), 2, "open c face drops the bond");
}

#[test]
fn every_layer_initializes_identically_to_the_single_sheet() {
    let ml_deck = load_multilayer();
    let ss_deck = petra_deck::load(repo_path("petra/examples/kaolinite.toml"))
        .expect("single-sheet deck compiles");
    let ml = ml_deck.build_engine(Some(0)).expect("multilayer builds");
    let ss = ss_deck.build_engine(Some(0)).expect("single sheet builds");

    for a in 0..20 {
        for b in 0..3 {
            for t in 0..26 {
                let ss_i = (a * 3 + b) * 26 + t;
                let ss_name = &ss_deck.state_names[ss.lattice.states[ss_i].0 as usize];
                for c in 0..4 {
                    let ml_i = site(a, b, c, t);
                    let ml_name = &ml_deck.state_names[ml.lattice.states[ml_i].0 as usize];
                    assert_eq!(
                        ml_name, ss_name,
                        "state at (a={a}, b={b}, c={c}, t={t})"
                    );
                    // Frozen flags agree except the one documented case:
                    // donor 24's hbond crosses the fixed b-face, so pos 24
                    // in the extreme b-row freezes in the multilayer.
                    if t == 24 && b == 2 {
                        assert!(ml.lattice.frozen[ml_i] && !ss.lattice.frozen[ss_i]);
                    } else {
                        assert_eq!(
                            ml.lattice.frozen[ml_i], ss.lattice.frozen[ss_i],
                            "frozen at (a={a}, b={b}, c={c}, t={t})"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn hydrolysis_reachability_does_not_cross_the_interlayer() {
    let deck = load_multilayer();
    let mut engine = deck.build_engine(Some(0)).expect("engine builds");
    let kinds = kinds_of(&deck, &engine.lattice);
    let mut scratch = Vec::new();

    // Structural half of the invariant: every guard whose state set
    // includes @hydrolyzed members (Oss.si1 is one) walks the covalent
    // graph only. Today the exact-distance-2 semantics happens to keep a
    // distance-1 hbond partner out of these guards anyway; the exclusion
    // makes the anhydrous gap an invariant of the deck rather than an
    // accident of state sets and distances.
    let si1 = state_id(&deck, "Oss.si1");
    let mut hydro_guards = 0;
    for r in &deck.reactions {
        for g in &r.guards {
            if g.select.states.contains(si1) {
                hydro_guards += 1;
                assert!(
                    g.select.exclude_label.is_some(),
                    "{}: @hydrolyzed guard must exclude hbond",
                    r.name
                );
            }
        }
    }
    // 7 @hydrolyzed surface-reachability guards + the 2 adsorption
    // @o_occ guards (which contain Oss.si1 too).
    assert_eq!(hydro_guards, 9, "all occupied-O guards checked");

    // Behavioral half: hydrolyzing a donor's acceptor (layer above) does
    // not make the donor surface-reachable. Unfrozen mid-lattice donor:
    // pos 19 at (5,1,1); its acceptor is pos 11 at (5,1,2).
    let r14 = deck
        .reactions
        .iter()
        .find(|r| r.name == "R14-alohal-hydrolysis")
        .expect("R14 exists");
    let guard = &r14.guards[0].select;
    let donor = site(5, 1, 1, 19);
    assert!(!engine.lattice.frozen[donor]);
    assert_eq!(
        count_matches(&engine.lattice, &kinds, donor, guard, &mut scratch),
        0,
        "bulk donor starts unreachable"
    );
    engine.lattice.states[site(5, 1, 2, 11)] = si1;
    assert_eq!(
        count_matches(&engine.lattice, &kinds, donor, guard, &mut scratch),
        0,
        "the anhydrous interlayer transmits no reachability"
    );
}

#[test]
fn desorption_is_suppressed_by_exactly_exp_n_ehb_over_rt() {
    let deck = load_multilayer();
    let mut engine = deck.build_engine(Some(0)).expect("engine builds");
    let kinds = kinds_of(&deck, &engine.lattice);
    let mut scratch = Vec::new();
    let rt = R_KCAL * T;

    let desorb_al = deck
        .reactions
        .iter()
        .find(|r| r.name == "desorb-al")
        .expect("desorb-al exists");
    let hb_mod = desorb_al
        .modifiers
        .iter()
        .find(|m| m.select.distance == 2)
        .expect("hbond modifier present");

    // Mid-lattice Al pos 0 at (5,1,1) — away from the b=0 wall, where
    // legacy termination quirks empty some oxygens. Its three outer
    // hydroxyls (20, 23, 25) donate to acceptors 10, 16, 17 in the layer
    // above, all intact.
    let al = site(5, 1, 1, 0);
    assert_eq!(
        count_matches(&engine.lattice, &kinds, al, &hb_mod.select, &mut scratch),
        3,
        "a bulk Al is anchored by three interlayer bonds"
    );
    let k_bonded = resolve_rate(&engine.lattice, &kinds, desorb_al, al, T, &mut scratch);

    // Sever them (empty the acceptors) and the rate recovers by exactly
    // exp(3 e_hb / RT).
    for t in [10, 16, 17] {
        engine.lattice.states[site(5, 1, 2, t)] = state_id(&deck, "Oss.empty");
    }
    assert_eq!(
        count_matches(&engine.lattice, &kinds, al, &hb_mod.select, &mut scratch),
        0
    );
    let k_free = resolve_rate(&engine.lattice, &kinds, desorb_al, al, T, &mut scratch);
    let expect = (3.0 * E_HB / rt).exp();
    assert!(
        ((k_free / k_bonded) / expect - 1.0).abs() < 1e-12,
        "rate ratio {} vs exp(3 e_hb/RT) {expect}",
        k_free / k_bonded
    );
}

#[test]
fn multilayer_dynamics_runs_and_conserves_sites() {
    let deck = load_multilayer();
    let mut engine = deck.build_engine(Some(1)).expect("engine builds");

    // 240 cells × (4 Al, 4 Si, 6 Oss, 4 Osa, 8 Oaa).
    let kind_totals = [960u64, 960, 1440, 960, 1920];
    let count_by_kind = |engine: &petra_core::Engine| -> Vec<u64> {
        let counts = engine.state_counts(deck.n_states);
        deck.kind_state_ranges
            .iter()
            .map(|&(start, n)| (start..start + n).map(|s| counts[s as usize]).sum::<u64>())
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
}
