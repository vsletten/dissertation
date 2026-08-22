use std::path::PathBuf;

fn deck_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../decks/ice-mantle-h2.toml")
}

fn load_deck() -> petra_deck::CompiledDeck {
    petra_deck::load(deck_path()).expect("the shipped H2 ice-mantle deck must load")
}

fn efficiency_at(deck: &petra_deck::CompiledDeck, temperature: f64, seed: u64) -> f64 {
    let mut engine = deck.build_engine(Some(seed)).expect("engine builds");
    engine
        .set_temperature(temperature)
        .expect("temperature is valid");

    let mut arrivals = 0_u64;
    let mut molecules = 0_u64;
    for _ in 0..250_000 {
        let fired = engine.step().expect("the open surface keeps advancing");
        let name = &deck.reactions[fired.reaction as usize].name;
        if name.starts_with("deposit_") {
            arrivals += 1;
        } else if name.starts_with("form_h2_") {
            molecules += 1;
        }
        if arrivals == 300 {
            return 2.0 * molecules as f64 / arrivals as f64;
        }
    }
    panic!("did not reach 300 H arrivals within the bounded event budget");
}

#[test]
fn deck_exercises_every_required_grain_surface_process() {
    let deck = load_deck();
    let names: Vec<&str> = deck
        .reactions
        .iter()
        .map(|reaction| reaction.name.as_str())
        .collect();

    for process in [
        "deposit_shallow",
        "deposit_deep",
        "thermal_hop_shallow",
        "thermal_hop_deep",
        "tunnel_hop_shallow",
        "tunnel_hop_deep",
        "form_h2_shallow",
        "form_h2_deep",
        "desorb_h_shallow",
        "desorb_h_deep",
        "desorb_h2_shallow",
        "desorb_h2_deep",
    ] {
        assert!(
            names.contains(&process),
            "missing required process {process}"
        );
    }
}

#[test]
fn quenched_disorder_is_seeded_and_contains_both_site_classes() {
    let deck = load_deck();
    let a = deck.build_engine(Some(41)).expect("first engine");
    let b = deck.build_engine(Some(41)).expect("second engine");
    let c = deck.build_engine(Some(42)).expect("different-seed engine");

    assert_eq!(
        a.lattice.states, b.lattice.states,
        "same seed is bitwise stable"
    );
    assert_ne!(
        a.lattice.states, c.lattice.states,
        "different seeds change disorder"
    );

    let deep_empty = deck
        .state_names
        .iter()
        .position(|name| name == "deep.empty")
        .expect("deep vacancy state") as u16;
    let deep = a
        .lattice
        .states
        .iter()
        .filter(|state| state.0 == deep_empty)
        .count();
    assert!(
        deep > 20 && deep < a.lattice.len() / 2,
        "plausible mixed surface"
    );
}

#[test]
fn h2_efficiency_has_the_canonical_finite_temperature_window() {
    let deck = load_deck();
    let cold = efficiency_at(&deck, 6.0, 7001);
    let window = efficiency_at(&deck, 12.0, 7001);
    let hot = efficiency_at(&deck, 24.0, 7001);
    eprintln!("efficiency diagnostics: cold={cold:.3}, window={window:.3}, hot={hot:.3}");

    assert!(window > 0.55, "efficient middle window: {window:.3}");
    assert!(
        window > cold + 0.12,
        "mobility-limited cold edge: cold={cold:.3}, window={window:.3}"
    );
    assert!(
        window > hot + 0.25,
        "desorption-limited hot edge: hot={hot:.3}, window={window:.3}"
    );
}
