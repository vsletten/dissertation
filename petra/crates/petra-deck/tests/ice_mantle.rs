use std::{collections::BTreeSet, path::PathBuf};

use petra_core::{Fired, StateId, StepCtx, StepOutcome, Stop, UpdateStrategy};

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

struct FireRule(u16);

impl UpdateStrategy for FireRule {
    fn step(&mut self, ctx: &mut StepCtx<'_>) -> Result<StepOutcome, Stop> {
        let site = ctx
            .apply
            .live_sites()
            .iter()
            .copied()
            .find(|&site| ctx.apply.enabled_rules(site).contains(&self.0))
            .ok_or(Stop::NoEvents)?;
        ctx.apply.apply_transition(site, self.0, ctx.rng)?;
        Ok(StepOutcome {
            fired: vec![Fired {
                step: 0,
                time: 0.0,
                site,
                reaction: self.0,
            }],
            dt: 0.0,
        })
    }
}

fn state_id(deck: &petra_deck::CompiledDeck, name: &str) -> StateId {
    StateId(
        deck.state_names
            .iter()
            .position(|candidate| candidate == name)
            .unwrap_or_else(|| panic!("missing state {name}")) as u16,
    )
}

#[allow(clippy::too_many_arguments)]
fn exercise_rule(
    deck: &petra_deck::CompiledDeck,
    rule_name: &str,
    center_empty: &str,
    neighbor_empty: &str,
    center_before: &str,
    neighbor_before: &str,
    center_after: &str,
    neighbor_after: &str,
) {
    let mut engine = deck.build_engine(Some(41)).expect("engine builds");
    let center_empty = state_id(deck, center_empty);
    let neighbor_empty = state_id(deck, neighbor_empty);
    let (center, neighbor) = engine
        .lattice
        .states
        .iter()
        .enumerate()
        .find_map(|(center, &state)| {
            (state == center_empty).then(|| {
                engine
                    .lattice
                    .neighbors(center)
                    .iter()
                    .map(|&site| site as usize)
                    .find(|&site| engine.lattice.states[site] == neighbor_empty)
                    .map(|neighbor| (center, neighbor))
            })?
        })
        .unwrap_or_else(|| {
            panic!("seeded surface has no {center_empty:?}/{neighbor_empty:?} pair")
        });

    let center_before = state_id(deck, center_before);
    let neighbor_before = state_id(deck, neighbor_before);
    let center_after = state_id(deck, center_after);
    let neighbor_after = state_id(deck, neighbor_after);
    engine.lattice.states[center] = center_before;
    engine.lattice.states[neighbor] = neighbor_before;
    engine
        .set_temperature(engine.temperature())
        .expect("refreshes propensities after arranging the transition fixture");

    let reaction = deck
        .reactions
        .iter()
        .position(|candidate| candidate.name == rule_name)
        .unwrap_or_else(|| panic!("missing reaction {rule_name}")) as u16;
    let outcome = engine
        .step_with(&mut FireRule(reaction))
        .unwrap_or_else(|stop| panic!("{rule_name} must be reachable: {stop}"));
    assert_eq!(outcome.fired.len(), 1, "{rule_name} fires once");
    assert_eq!(outcome.fired[0].reaction, reaction, "{rule_name} fired");

    let changes = engine.last_changes();
    assert!(
        changes.contains(&(center, center_before, center_after)),
        "{rule_name} applies its center rewrite"
    );
    assert!(
        changes.iter().any(|&(site, old, new)| {
            engine.lattice.neighbors(center).contains(&(site as u32))
                && old == neighbor_before
                && new == neighbor_after
        }),
        "{rule_name} applies its neighbor rewrite"
    );
}

#[test]
fn deck_declares_the_complete_process_inventory() {
    let deck = load_deck();
    let names: BTreeSet<&str> = deck
        .reactions
        .iter()
        .map(|reaction| reaction.name.as_str())
        .collect();
    let expected = BTreeSet::from([
        "deposit_shallow",
        "deposit_deep",
        "thermal_hop_shallow",
        "thermal_hop_shallow_to_deep",
        "thermal_hop_deep",
        "thermal_hop_deep_to_shallow",
        "tunnel_hop_shallow",
        "tunnel_hop_shallow_to_deep",
        "tunnel_hop_deep",
        "tunnel_hop_deep_to_shallow",
        "form_h2_shallow",
        "form_h2_shallow_with_deep",
        "form_h2_deep",
        "form_h2_deep_with_shallow",
        "desorb_h_shallow",
        "desorb_h_deep",
        "desorb_h2_shallow",
        "desorb_h2_deep",
    ]);
    assert_eq!(names, expected, "deck process inventory changed");
}

#[test]
fn every_diffusion_and_recombination_rule_is_reachable_and_rewrites_both_sites() {
    let deck = load_deck();
    for (rule, source_kind, target_kind) in [
        ("thermal_hop_shallow", "shallow", "shallow"),
        ("thermal_hop_shallow_to_deep", "shallow", "deep"),
        ("thermal_hop_deep", "deep", "deep"),
        ("thermal_hop_deep_to_shallow", "deep", "shallow"),
        ("tunnel_hop_shallow", "shallow", "shallow"),
        ("tunnel_hop_shallow_to_deep", "shallow", "deep"),
        ("tunnel_hop_deep", "deep", "deep"),
        ("tunnel_hop_deep_to_shallow", "deep", "shallow"),
    ] {
        exercise_rule(
            &deck,
            rule,
            &format!("{source_kind}.empty"),
            &format!("{target_kind}.empty"),
            &format!("{source_kind}.H"),
            &format!("{target_kind}.empty"),
            &format!("{source_kind}.empty"),
            &format!("{target_kind}.H_mobile"),
        );
    }

    for (rule, center_kind, neighbor_kind) in [
        ("form_h2_shallow", "shallow", "shallow"),
        ("form_h2_shallow_with_deep", "shallow", "deep"),
        ("form_h2_deep", "deep", "deep"),
        ("form_h2_deep_with_shallow", "deep", "shallow"),
    ] {
        exercise_rule(
            &deck,
            rule,
            &format!("{center_kind}.empty"),
            &format!("{neighbor_kind}.empty"),
            &format!("{center_kind}.H_mobile"),
            &format!("{neighbor_kind}.H"),
            &format!("{center_kind}.H2"),
            &format!("{neighbor_kind}.empty"),
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
