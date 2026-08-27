use petra_core::{CtmcAdvance, StateId};
use std::collections::BTreeSet;
use std::path::PathBuf;

fn deck_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../decks/passive-metal.toml")
}

fn deck_text() -> String {
    std::fs::read_to_string(deck_path()).expect("passive-metal deck fixture")
}

fn compiled_at(potential_mu: f64) -> petra_deck::CompiledDeck {
    let mut parsed: petra_deck::DeckFile =
        toml::from_str(&deck_text()).expect("passive-metal deck parses");
    parsed.thermo.mu.insert("e_drive".to_string(), potential_mu);
    petra_deck::compile(&parsed).expect("passive-metal deck compiles")
}

fn state_ids(deck: &petra_deck::CompiledDeck, names: &[&str]) -> Vec<StateId> {
    names
        .iter()
        .map(|name| {
            StateId(
                deck.state_names
                    .iter()
                    .position(|candidate| candidate == name)
                    .unwrap_or_else(|| panic!("missing state {name}")) as u16,
            )
        })
        .collect()
}

fn cluster_stats(engine: &petra_core::Engine, selected: &[StateId]) -> (u64, u64) {
    let mut seen = vec![false; engine.lattice.len()];
    let mut stack = Vec::new();
    let mut count = 0_u64;
    let mut largest = 0_u64;
    for start in 0..engine.lattice.len() {
        if seen[start] || !selected.contains(&engine.lattice.states[start]) {
            continue;
        }
        seen[start] = true;
        stack.push(start);
        count += 1;
        let mut size = 0_u64;
        while let Some(site) = stack.pop() {
            size += 1;
            for &neighbor in engine.lattice.neighbors(site) {
                let neighbor = neighbor as usize;
                if !seen[neighbor] && selected.contains(&engine.lattice.states[neighbor]) {
                    seen[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
        largest = largest.max(size);
    }
    (count, largest)
}

fn largest_cluster(engine: &petra_core::Engine, selected: &[StateId]) -> u64 {
    cluster_stats(engine, selected).1
}

fn reaction_rate(engine: &petra_core::Engine, reaction: u16) -> f64 {
    engine
        .site_event_rates()
        .iter()
        .flat_map(|events| events.iter())
        .filter_map(|&(rule, rate)| (rule == reaction).then_some(rate))
        .sum()
}

#[test]
fn passive_film_deck_declares_thermo_autocatalysis_defect_and_full_state_machine() {
    let deck = compiled_at(0.0);
    assert_eq!(deck.strategy.as_str(), "ctmc");
    assert_eq!(deck.n_replicas, 32);

    let states: BTreeSet<_> = deck.state_names.iter().map(String::as_str).collect();
    for state in [
        "film_patch.intact",
        "film_patch.adsorbed",
        "film_patch.thinned",
        "film_patch.bare",
        "film_patch.dissolving",
        "film_patch.repassivated",
    ] {
        assert!(states.contains(state), "missing {state}");
    }
    let reactions: BTreeSet<_> = deck
        .reactions
        .iter()
        .map(|reaction| reaction.name.as_str())
        .collect();
    for reaction in [
        "chloride_adsorption",
        "chloride_desorption",
        "film_thinning",
        "film_rupture",
        "dissolution_activate",
        "metal_dissolution",
        "repassivate_bare",
        "repassivate_dissolving",
        "film_healing",
    ] {
        assert!(reactions.contains(reaction), "missing {reaction}");
    }

    let engine = deck.build_engine(Some(17)).expect("engine builds");
    let max_strain = engine
        .lattice
        .strain
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let min_strain = engine
        .lattice
        .strain
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    assert!(max_strain > 0.25, "declared defect has a strong core");
    assert!(
        max_strain > min_strain * 10.0,
        "defect field must be spatially resolved: min={min_strain}, max={max_strain}"
    );

    let low = compiled_at(-0.7);
    let high = compiled_at(0.8);
    let adsorption = low
        .reactions
        .iter()
        .position(|reaction| reaction.name == "chloride_adsorption")
        .expect("adsorption rule") as u16;
    let low_engine = low.build_engine(Some(17)).expect("low-potential engine");
    let high_engine = high.build_engine(Some(17)).expect("high-potential engine");
    let ratio = reaction_rate(&high_engine, adsorption) / reaction_rate(&low_engine, adsorption);
    assert!(
        ratio > 10.0,
        "thermo.mu(e_drive) must accelerate chloride attack: ratio={ratio}"
    );
}

fn cluster_at_horizon(deck: &petra_deck::CompiledDeck, seed: u64) -> (u64, u64, bool) {
    let pit = state_ids(deck, &["film_patch.bare", "film_patch.dissolving"]);
    let mut engine = deck.build_engine(Some(seed)).expect("engine builds");
    let mut initiated = false;
    let mut event = 0_u64;
    while event < 15_000 {
        match engine
            .advance_ctmc_until(400.0)
            .expect("open-film trajectory advances")
        {
            CtmcAdvance::Fired(_) => event += 1,
            CtmcAdvance::Deadline { time } => {
                assert_eq!(time.to_bits(), 400.0_f64.to_bits());
                break;
            }
        }
        if event.is_multiple_of(20) || !initiated {
            initiated |= largest_cluster(&engine, &pit) >= 6;
        }
    }
    assert_eq!(engine.time.to_bits(), 400.0_f64.to_bits());
    let (count, largest) = cluster_stats(&engine, &pit);
    (largest, count, initiated)
}

#[test]
fn potential_sweep_has_a_cooperative_pit_clustering_transition() {
    let low = compiled_at(-0.7);
    let high = compiled_at(0.8);
    let mut low_clusters = Vec::new();
    let mut high_clusters = Vec::new();
    let mut low_pit_counts = Vec::new();
    let mut high_pit_counts = Vec::new();
    let mut low_initiated = 0;
    let mut high_initiated = 0;
    for replica in 0..8 {
        let seed = petra_deck::replica_seed(low.seed, replica, low.seed_policy);
        let (cluster, pit_count, initiated) = cluster_at_horizon(&low, seed);
        low_clusters.push(cluster);
        low_pit_counts.push(pit_count);
        low_initiated += usize::from(initiated);

        let (cluster, pit_count, initiated) = cluster_at_horizon(&high, seed);
        high_clusters.push(cluster);
        high_pit_counts.push(pit_count);
        high_initiated += usize::from(initiated);
    }
    let low_mean = low_clusters.iter().sum::<u64>() as f64 / low_clusters.len() as f64;
    let high_mean = high_clusters.iter().sum::<u64>() as f64 / high_clusters.len() as f64;
    let low_pit_mean = low_pit_counts.iter().sum::<u64>() as f64 / low_pit_counts.len() as f64;
    let high_pit_mean = high_pit_counts.iter().sum::<u64>() as f64 / high_pit_counts.len() as f64;
    eprintln!(
        "corrosion transition: low={low_clusters:?} ({low_initiated}/8), high={high_clusters:?} ({high_initiated}/8)"
    );
    assert_eq!(low_initiated, 0, "subcritical film remains metastable");
    assert!(
        high_initiated >= 6,
        "supercritical film initiates in most replicas"
    );
    assert!(
        high_mean > low_mean + 2.0,
        "clusters grow cooperatively: low={low_mean}, high={high_mean}"
    );
    assert!(
        high_pit_mean > low_pit_mean * 3.0,
        "pit-number density rises sharply: low={low_pit_mean}, high={high_pit_mean}"
    );
}

#[test]
fn same_seed_reproduces_the_same_corrosion_trajectory_and_another_seed_diverges() {
    let deck = compiled_at(0.3);
    let run = |seed| {
        let mut engine = deck.build_engine(Some(seed)).expect("engine builds");
        for _ in 0..1_000 {
            engine.step().expect("trajectory advances");
        }
        (engine.time.to_bits(), engine.lattice.states)
    };
    let first = run(91);
    let repeat = run(91);
    let different = run(92);
    assert_eq!(first, repeat, "same seed must be bitwise reproducible");
    assert_ne!(
        first, different,
        "different seeds must select different trajectories"
    );
}

#[test]
fn mid_transition_trajectory_contains_current_pulses_and_repassivation_drops() {
    let deck = compiled_at(0.3);
    let seed = petra_deck::replica_seed(deck.seed, 3, deck.seed_policy);
    let pit = state_ids(&deck, &["film_patch.bare", "film_patch.dissolving"]);
    let dissolution = deck
        .reactions
        .iter()
        .position(|reaction| reaction.name == "metal_dissolution")
        .expect("metal dissolution") as u16;
    let mut engine = deck.build_engine(Some(seed)).expect("engine builds");
    let mut previous_active = 0_u64;
    let mut previous_rate = 0.0;
    let mut active_drop = false;
    let mut current_drop = false;
    let mut dissolution_events = 0;
    let mut repassivation_events = 0;

    let mut event = 0_u64;
    while event < 15_000 {
        let fired = match engine
            .advance_ctmc_until(300.0)
            .expect("trajectory advances")
        {
            CtmcAdvance::Fired(fired) => fired,
            CtmcAdvance::Deadline { time } => {
                assert_eq!(time.to_bits(), 300.0_f64.to_bits());
                break;
            }
        };
        event += 1;
        let name = &deck.reactions[fired.reaction as usize].name;
        dissolution_events += usize::from(name == "metal_dissolution");
        repassivation_events += usize::from(name.starts_with("repassivate_"));
        if event.is_multiple_of(20) {
            let active = engine
                .lattice
                .states
                .iter()
                .filter(|&&state| pit.contains(&state))
                .count() as u64;
            let rate = reaction_rate(&engine, dissolution);
            active_drop |= active < previous_active;
            current_drop |= rate < previous_rate;
            previous_active = active;
            previous_rate = rate;
        }
    }
    assert_eq!(engine.time.to_bits(), 300.0_f64.to_bits());
    assert!(dissolution_events > 25, "anodic current events occur");
    assert!(repassivation_events > 25, "repair competes with attack");
    assert!(active_drop, "a metastable active population must collapse");
    assert!(
        current_drop,
        "event-rate current proxy must contain downward pulses"
    );
}
