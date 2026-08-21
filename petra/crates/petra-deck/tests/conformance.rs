use std::path::PathBuf;

fn conformance_text(name: &str) -> String {
    std::fs::read_to_string(
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../..")
            .join("petra/examples/conformance")
            .join(name),
    )
    .unwrap_or_else(|error| panic!("read conformance deck {name}: {error}"))
}

#[test]
fn square_grid_compiles_to_periodic_csr_adjacency() {
    let text = r#"
[deck]
name = "grid-compile"
schema = 2

[structure]
kind = "grid"
grid = "square"
neighborhood = "von_neumann"
dims = [3, 3]
boundary = ["periodic", "periodic"]
default_kind = "site"

[[structure.kinds]]
name = "site"
initial = "empty"
[[structure.kinds.states]]
name = "empty"
occupant = "vacant"

[dynamics.thermo]
temperature = 300.0

[execution]
strategy = "ctmc"
[execution.stop]
steps = 1
"#;
    let parsed: petra_deck::DeckFile = toml::from_str(text).expect("grid parses");
    let deck = petra_deck::compile(&parsed).expect("grid compiles");
    let engine = deck.build_engine(Some(7)).expect("grid engine builds");

    assert_eq!(engine.lattice.len(), 9);
    assert_eq!(engine.lattice.dims, [3, 3, 1]);
    for site in 0..engine.lattice.len() {
        assert_eq!(engine.lattice.neighbors(site).len(), 4, "site {site}");
    }
}

#[test]
fn tiny_periodic_grids_have_simple_adjacency_without_loops_or_duplicates() {
    for (dims, neighborhood, expected_degree) in [
        ("[1, 1]", "von_neumann", 0),
        ("[2, 2]", "von_neumann", 2),
        ("[2, 2]", "moore", 3),
    ] {
        let text = format!(
            r#"
[deck]
name = "tiny-grid"
schema = 2
[structure]
kind = "grid"
grid = "square"
neighborhood = "{neighborhood}"
dims = {dims}
boundary = ["periodic", "periodic"]
default_kind = "site"
[[structure.kinds]]
name = "site"
initial = "empty"
[[structure.kinds.states]]
name = "empty"
occupant = "vacant"
[dynamics.thermo]
temperature = 1.0
[execution]
strategy = "ctmc"
"#
        );
        let parsed: petra_deck::DeckFile = toml::from_str(&text).expect("tiny grid parses");
        let deck = petra_deck::compile(&parsed).expect("tiny grid compiles");
        let engine = deck.build_engine(Some(1)).expect("tiny grid builds");
        for site in 0..engine.lattice.len() {
            let neighbors = engine.lattice.neighbors(site);
            assert_eq!(neighbors.len(), expected_degree, "{dims} {neighborhood}");
            assert!(!neighbors.contains(&(site as u32)), "self-loop at {site}");
            let mut unique = neighbors.to_vec();
            unique.sort_unstable();
            unique.dedup();
            assert_eq!(unique.len(), neighbors.len(), "duplicate at {site}");
        }
    }
}

#[test]
fn every_rfc_grid_family_compiles_with_expected_coordination() {
    for (family, neighborhood, dims, expected_degree) in [
        ("square", Some("moore"), "[4, 4]", 8),
        ("hex", None, "[4, 4]", 6),
        ("cubic", Some("von_neumann"), "[4, 4, 4]", 6),
        ("cubic", Some("moore"), "[4, 4, 4]", 26),
    ] {
        let neighborhood = neighborhood
            .map(|name| format!("neighborhood = \"{name}\""))
            .unwrap_or_default();
        let boundary = if family == "cubic" {
            "[\"periodic\", \"periodic\", \"periodic\"]"
        } else {
            "[\"periodic\", \"periodic\"]"
        };
        let text = format!(
            r#"
[deck]
name = "{family}-grid"
schema = 2
[structure]
kind = "grid"
grid = "{family}"
{neighborhood}
dims = {dims}
boundary = {boundary}
default_kind = "site"
[[structure.kinds]]
name = "site"
initial = "empty"
[[structure.kinds.states]]
name = "empty"
occupant = "vacant"
[dynamics.thermo]
temperature = 1.0
[execution]
strategy = "ctmc"
"#
        );
        let parsed: petra_deck::DeckFile = toml::from_str(&text).expect("grid parses");
        let deck = petra_deck::compile(&parsed).expect("grid compiles");
        let engine = deck.build_engine(Some(1)).expect("grid builds");
        assert!(
            (0..engine.lattice.len())
                .all(|site| { engine.lattice.neighbors(site).len() == expected_degree }),
            "{family} {neighborhood}"
        );
    }
}

const CONWAY: &str = r#"
[deck]
name = "conway-glider"
schema = 2

[structure]
kind = "grid"
grid = "square"
neighborhood = "moore"
dims = [8, 8]
boundary = ["periodic", "periodic"]
default_kind = "cell"

[[structure.species]]
name = "life"
[[structure.kinds]]
name = "cell"
initial = "dead"
[[structure.kinds.states]]
name = "dead"
occupant = "vacant"
[[structure.kinds.states]]
name = "alive"
occupant = "life"

[[structure.init]]
name = "glider"
center = { kind = "cell", state = ["dead"] }
sites = [[1, 0, 0, 0], [2, 1, 0, 0], [0, 2, 0, 0], [1, 2, 0, 0], [2, 2, 0, 0]]
set = "alive"

[dynamics.thermo]
temperature = 1.0

[[dynamics.rules]]
name = "underpopulation"
center = { kind = "cell", state = ["alive"] }
guards = [{ distance = 1, kind = "cell", state = ["alive"], max = 1 }]
[[dynamics.rules.effects]]
target = "center"
set = "dead"

[[dynamics.rules]]
name = "overpopulation"
center = { kind = "cell", state = ["alive"] }
guards = [{ distance = 1, kind = "cell", state = ["alive"], min = 4 }]
[[dynamics.rules.effects]]
target = "center"
set = "dead"

[[dynamics.rules]]
name = "birth"
center = { kind = "cell", state = ["dead"] }
guards = [{ distance = 1, kind = "cell", state = ["alive"], min = 3, max = 3 }]
[[dynamics.rules.effects]]
target = "center"
set = "alive"

[execution]
strategy = "synchronous"
[execution.stop]
steps = 8
"#;

fn alive_sites(
    deck: &petra_deck::CompiledDeck,
    engine: &petra_core::Engine,
) -> Vec<(usize, usize)> {
    let alive = deck
        .state_names
        .iter()
        .position(|name| name == "cell.alive")
        .expect("alive state") as u16;
    engine
        .lattice
        .states
        .iter()
        .enumerate()
        .filter(|(_, state)| state.0 == alive)
        .map(|(site, _)| {
            let ([x, y, _], _) = engine.lattice.coords(site);
            (x, y)
        })
        .collect()
}

#[test]
fn conway_glider_has_period_four_and_displacement_one_one() {
    let text = conformance_text("conway.toml");
    let parsed: petra_deck::DeckFile = toml::from_str(&text).expect("Conway deck parses");
    let deck = petra_deck::compile(&parsed).expect("Conway deck compiles");
    let mut engine = deck.build_engine(Some(42)).expect("engine builds");
    let initial = alive_sites(&deck, &engine);
    let mut strategy = deck.strategy();

    for _ in 0..4 {
        engine
            .step_with(&mut strategy)
            .expect("synchronous step succeeds");
    }

    let shifted: Vec<_> = initial.iter().map(|&(x, y)| (x + 1, y + 1)).collect();
    assert_eq!(alive_sites(&deck, &engine), shifted);
    assert_eq!(engine.time.to_bits(), 4.0f64.to_bits());
}

#[test]
fn conway_blinker_has_period_two_and_block_is_still() {
    let text = conformance_text("conway.toml");
    let build_pattern = |cells: &[(usize, usize)], seed| {
        let mut parsed: petra_deck::DeckFile = toml::from_str(&text).expect("Conway deck parses");
        parsed.init[0].sites = Some(cells.iter().map(|&(x, y)| [x, y, 0, 0]).collect());
        let deck = petra_deck::compile(&parsed).expect("Conway pattern compiles");
        let engine = deck.build_engine(Some(seed)).expect("pattern engine");
        (deck, engine)
    };

    let horizontal = vec![(1, 2), (2, 2), (3, 2)];
    let vertical = vec![(2, 1), (2, 2), (2, 3)];
    let (blinker_deck, mut blinker) = build_pattern(&horizontal, 1);
    let mut strategy = blinker_deck.strategy();
    blinker
        .step_with(&mut strategy)
        .expect("first blinker step");
    assert_eq!(alive_sites(&blinker_deck, &blinker), vertical);
    blinker
        .step_with(&mut strategy)
        .expect("second blinker step");
    assert_eq!(alive_sites(&blinker_deck, &blinker), horizontal);

    let block_pattern = vec![(1, 1), (1, 2), (2, 1), (2, 2)];
    let (block_deck, mut block) = build_pattern(&block_pattern, 2);
    let mut strategy = block_deck.strategy();
    block.step_with(&mut strategy).expect("block step");
    assert_eq!(alive_sites(&block_deck, &block), block_pattern);
}

#[test]
fn synchronous_rules_reject_hidden_random_draws() {
    let random_branches = CONWAY.replace(
        "[[dynamics.rules.effects]]\ntarget = \"center\"\nset = \"dead\"",
        "[[dynamics.rules.branches]]\nweight = 1.0\n[[dynamics.rules.branches.effects]]\ntarget = \"center\"\nset = \"dead\"\n[[dynamics.rules.branches]]\nweight = 1.0\n[[dynamics.rules.branches.effects]]\ntarget = \"center\"\nset = \"dead\"",
    );
    let error = toml::from_str::<petra_deck::DeckFile>(&random_branches)
        .expect_err("synchronous CA must remain draw-free");
    assert!(error.to_string().contains("draw-free"), "{error}");
}

const ISING: &str = r#"
[deck]
name = "ising-2d"
schema = 2
units = "kcal/mol"

[structure]
kind = "grid"
grid = "square"
neighborhood = "von_neumann"
dims = [8, 8]
boundary = ["periodic", "periodic"]
default_kind = "spin"

[[structure.species]]
name = "plus"
[[structure.species]]
name = "minus"
[[structure.kinds]]
name = "spin"
initial = "up"
[[structure.kinds.states]]
name = "up"
occupant = "plus"
[[structure.kinds.states]]
name = "down"
occupant = "minus"

[[structure.init]]
name = "random-spins"
center = { kind = "spin", state = ["up"] }
probability = 0.5
set = "down"

[dynamics.thermo]
temperature = 2.269185314

[[dynamics.rules]]
name = "up-to-down"
center = { kind = "spin", state = ["up"] }
rate = { energy = { delta = -0.015896 } }
[[dynamics.rules.modifiers]]
select = { distance = 1, kind = "spin", state = ["up"] }
per_match = { dea = 0.007948 }
[[dynamics.rules.effects]]
target = "center"
set = "down"

[[dynamics.rules]]
name = "down-to-up"
center = { kind = "spin", state = ["down"] }
rate = { energy = { delta = -0.015896 } }
[[dynamics.rules.modifiers]]
select = { distance = 1, kind = "spin", state = ["down"] }
per_match = { dea = 0.007948 }
[[dynamics.rules.effects]]
target = "center"
set = "up"

[execution]
strategy = "metropolis"
[execution.metropolis]
temperature = 2.269185314
[execution.stop]
steps = 100000
[execution.ensemble]
seed = 42
n_replicas = 10
seed_policy = "hash"
"#;

fn binder_cumulant(size: usize, temperature: f64) -> f64 {
    let text = conformance_text("ising.toml");
    let mut parsed: petra_deck::DeckFile = toml::from_str(&text).expect("Ising deck parses");
    parsed.grid.as_mut().expect("grid").dims = vec![size, size];
    parsed
        .execution
        .metropolis
        .as_mut()
        .expect("Metropolis parameters")
        .temperature = Some(temperature);
    let deck = petra_deck::compile(&parsed).expect("Ising deck compiles");
    let up = deck
        .state_names
        .iter()
        .position(|name| name == "spin.up")
        .expect("up state") as u16;
    let sites = size * size;
    let mut m2 = 0.0;
    let mut m4 = 0.0;
    let mut samples = 0.0;

    for replica in 0..deck.n_replicas {
        let seed = petra_deck::replica_seed(deck.seed, replica, deck.seed_policy);
        let mut engine = deck.build_engine(Some(seed)).expect("engine builds");
        let mut strategy = deck.strategy();
        for _ in 0..300 * sites {
            engine.step_with(&mut strategy).expect("burn-in step");
        }
        for _ in 0..500 {
            for _ in 0..sites {
                engine.step_with(&mut strategy).expect("sample sweep");
            }
            let n_up = engine
                .lattice
                .states
                .iter()
                .filter(|state| state.0 == up)
                .count() as f64;
            let magnetization = (2.0 * n_up - sites as f64) / sites as f64;
            let square = magnetization * magnetization;
            m2 += square;
            m4 += square * square;
            samples += 1.0;
        }
    }
    let mean2 = m2 / samples;
    let mean4 = m4 / samples;
    1.0 - mean4 / (3.0 * mean2 * mean2)
}

#[test]
fn ising_binder_crossing_is_near_onsager_critical_temperature() {
    const ONSAGER_TC: f64 = 2.269_185_314_213_022;
    const TOLERANCE: f64 = 0.25;
    let temperatures = [2.0, ONSAGER_TC, 2.6];
    let differences: Vec<_> = temperatures
        .iter()
        .map(|&temperature| binder_cumulant(8, temperature) - binder_cumulant(12, temperature))
        .collect();
    let crossing = temperatures
        .windows(2)
        .zip(differences.windows(2))
        .find_map(|(temperature, difference)| {
            if difference[0] * difference[1] <= 0.0 {
                Some(
                    temperature[0]
                        - difference[0] * (temperature[1] - temperature[0])
                            / (difference[1] - difference[0]),
                )
            } else {
                None
            }
        })
        .unwrap_or_else(|| panic!("Binder curves did not cross: differences={differences:?}"));
    eprintln!(
        "Binder crossing={crossing:.6}, Onsager={ONSAGER_TC:.6}, differences={differences:?}"
    );

    assert!(
        (crossing - ONSAGER_TC).abs() <= TOLERANCE,
        "crossing {crossing:.6}, Onsager {ONSAGER_TC:.6}, tolerance {TOLERANCE}, differences={differences:?}"
    );
}

const SIR: &str = r#"
[deck]
name = "sir-pca"
schema = 2

[structure]
kind = "grid"
grid = "square"
neighborhood = "von_neumann"
dims = [32, 32]
boundary = ["periodic", "periodic"]
default_kind = "person"

[[structure.species]]
name = "host"
[[structure.kinds]]
name = "person"
initial = "susceptible"
[[structure.kinds.states]]
name = "susceptible"
occupant = "host"
[[structure.kinds.states]]
name = "infected"
occupant = "host"
[[structure.kinds.states]]
name = "recovered"
occupant = "host"

[[structure.init]]
name = "seed-infections"
center = { kind = "person", state = ["susceptible"] }
probability = 0.01
set = "infected"

[dynamics.thermo]
temperature = 1.0

[[dynamics.rules]]
name = "infection"
center = { kind = "person", state = ["susceptible"] }
guards = [{ distance = 1, kind = "person", state = ["infected"], min = 1 }]
rate = { probability = 0.05 }
[[dynamics.rules.effects]]
target = "center"
set = "infected"

[[dynamics.rules]]
name = "recovery"
center = { kind = "person", state = ["infected"] }
rate = { probability = 0.2 }
[[dynamics.rules.effects]]
target = "center"
set = "recovered"

[execution]
strategy = "pca"
[execution.pca]
[execution.stop]
steps = 100
[execution.ensemble]
seed = 42
n_replicas = 128
seed_policy = "hash"
"#;

fn sir_one_step_growth(r0: f64) -> f64 {
    const GAMMA: f64 = 0.2;
    const COORDINATION: f64 = 4.0;
    let text = conformance_text("sir.toml");
    let mut parsed: petra_deck::DeckFile = toml::from_str(&text).expect("SIR deck parses");
    parsed.reactions[0]
        .rate
        .as_mut()
        .expect("infection rate")
        .probability = Some(r0 * GAMMA / COORDINATION);
    let deck = petra_deck::compile(&parsed).expect("SIR deck compiles");
    let infected = deck
        .state_names
        .iter()
        .position(|name| name == "person.infected")
        .expect("infected state") as u16;
    let mut before = 0u64;
    let mut after = 0u64;
    for replica in 0..deck.n_replicas {
        let seed = petra_deck::replica_seed(deck.seed, replica, deck.seed_policy);
        let mut engine = deck.build_engine(Some(seed)).expect("engine builds");
        before += engine
            .lattice
            .states
            .iter()
            .filter(|state| state.0 == infected)
            .count() as u64;
        let mut strategy = deck.strategy();
        engine.step_with(&mut strategy).expect("PCA step");
        after += engine
            .lattice
            .states
            .iter()
            .filter(|state| state.0 == infected)
            .count() as u64;
    }
    after as f64 / before as f64
}

#[test]
fn sir_growth_threshold_matches_mean_field_r0_one() {
    const TOLERANCE: f64 = 0.15;
    let low_r0 = 0.5;
    let high_r0 = 1.5;
    let low_growth = sir_one_step_growth(low_r0);
    let high_growth = sir_one_step_growth(high_r0);
    assert!(low_growth < 1.0, "R0={low_r0}: growth={low_growth}");
    assert!(high_growth > 1.0, "R0={high_r0}: growth={high_growth}");
    let threshold = low_r0 + (1.0 - low_growth) * (high_r0 - low_r0) / (high_growth - low_growth);
    eprintln!(
        "SIR threshold={threshold:.6}, growth({low_r0})={low_growth:.6}, growth({high_r0})={high_growth:.6}"
    );
    assert!(
        (threshold - 1.0).abs() <= TOLERANCE,
        "estimated R0 threshold {threshold:.6}, tolerance {TOLERANCE}; growth({low_r0})={low_growth:.6}, growth({high_r0})={high_growth:.6}"
    );
}

fn sir_outbreak_probability(r0: f64) -> f64 {
    const GAMMA: f64 = 0.2;
    const COORDINATION: f64 = 4.0;

    let text = conformance_text("sir.toml");
    let mut parsed: petra_deck::DeckFile = toml::from_str(&text).expect("SIR deck parses");
    parsed.reactions[0]
        .rate
        .as_mut()
        .expect("infection rate")
        .probability = Some(r0 * GAMMA / COORDINATION);
    let steps = parsed.execution.stop.steps.expect("SIR stop steps");
    let deck = petra_deck::compile(&parsed).expect("SIR deck compiles");
    let infected = deck
        .state_names
        .iter()
        .position(|name| name == "person.infected")
        .expect("infected state") as u16;
    let recovered = deck
        .state_names
        .iter()
        .position(|name| name == "person.recovered")
        .expect("recovered state") as u16;
    let outbreak_size = engine_outbreak_size(deck.dims[0] * deck.dims[1]);
    let outbreaks = (0..deck.n_replicas)
        .filter(|&replica| {
            let seed = petra_deck::replica_seed(deck.seed, replica, deck.seed_policy);
            let mut engine = deck.build_engine(Some(seed)).expect("engine builds");
            let mut strategy = deck.strategy();
            for _ in 0..steps {
                engine.step_with(&mut strategy).expect("PCA step");
            }
            engine
                .lattice
                .states
                .iter()
                .filter(|state| state.0 == infected || state.0 == recovered)
                .count()
                >= outbreak_size
        })
        .count();
    outbreaks as f64 / deck.n_replicas as f64
}

fn engine_outbreak_size(site_count: usize) -> usize {
    // A macroscopic outbreak must reach at least 10% of the finite lattice;
    // the seeded 1% initial prevalence cannot satisfy this by itself.
    site_count.div_ceil(10)
}

#[test]
fn sir_outbreak_probability_changes_across_mean_field_threshold() {
    let below = sir_outbreak_probability(0.5);
    let above = sir_outbreak_probability(2.0);
    eprintln!("SIR outbreak P(R0=0.5)={below:.6}, P(R0=2.0)={above:.6}");
    assert!(below <= 0.05, "subcritical outbreak probability {below}");
    assert!(above >= 0.40, "supercritical outbreak probability {above}");
}

#[test]
fn shipped_decks_compile_to_the_models_exercised_by_the_analytic_gates() {
    for (name, expected) in [
        ("conway.toml", CONWAY),
        ("ising.toml", ISING),
        ("sir.toml", SIR),
    ] {
        let actual_text = conformance_text(name);
        let actual: petra_deck::DeckFile =
            toml::from_str(&actual_text).unwrap_or_else(|error| panic!("{name} parses: {error}"));
        let expected: petra_deck::DeckFile =
            toml::from_str(expected).unwrap_or_else(|error| panic!("{name} gate parses: {error}"));
        let actual =
            petra_deck::compile(&actual).unwrap_or_else(|error| panic!("{name} compiles: {error}"));
        let expected = petra_deck::compile(&expected)
            .unwrap_or_else(|error| panic!("{name} gate compiles: {error}"));
        assert_eq!(actual.state_names, expected.state_names, "{name} states");
        assert_eq!(actual.dims, expected.dims, "{name} dimensions");
        assert_eq!(actual.boundary, expected.boundary, "{name} boundary");
        assert_eq!(actual.init_passes, expected.init_passes, "{name} init");
        assert_eq!(actual.reactions, expected.reactions, "{name} rules");
        assert_eq!(actual.strategy, expected.strategy, "{name} strategy");
    }
}

fn trajectory_fingerprint(name: &str, seed: u64, steps: usize) -> u64 {
    fn mix(hash: &mut u64, value: u64) {
        for byte in value.to_le_bytes() {
            *hash ^= u64::from(byte);
            *hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }

    let text = conformance_text(name);
    let parsed: petra_deck::DeckFile = toml::from_str(&text).expect("deck parses");
    let deck = petra_deck::compile(&parsed).expect("deck compiles");
    let mut engine = deck.build_engine(Some(seed)).expect("engine builds");
    let mut strategy = deck.strategy();
    let mut hash = 0xcbf2_9ce4_8422_2325;
    for _ in 0..steps {
        let outcome = engine.step_with(&mut strategy).expect("strategy step");
        mix(&mut hash, outcome.dt.to_bits());
        mix(&mut hash, outcome.fired.len() as u64);
        for fired in outcome.fired {
            mix(&mut hash, fired.step);
            mix(&mut hash, fired.time.to_bits());
            mix(&mut hash, fired.site as u64);
            mix(&mut hash, u64::from(fired.reaction));
        }
    }
    for state in &engine.lattice.states {
        mix(&mut hash, u64::from(state.0));
    }
    hash
}

#[test]
fn pca_rejects_conflicting_simultaneous_writes() {
    let text = r#"
[deck]
name = "pca-conflict"
schema = 2
[structure]
kind = "grid"
grid = "square"
neighborhood = "von_neumann"
dims = [1, 1]
boundary = ["periodic", "periodic"]
default_kind = "site"
[[structure.species]]
name = "x"
[[structure.kinds]]
name = "site"
initial = "a"
[[structure.kinds.states]]
name = "a"
occupant = "x"
[[structure.kinds.states]]
name = "b"
occupant = "x"
[[structure.kinds.states]]
name = "c"
occupant = "x"
[dynamics.thermo]
temperature = 1.0
[[dynamics.rules]]
name = "a-to-b"
center = { kind = "site", state = ["a"] }
rate = { probability = 1.0 }
[[dynamics.rules.effects]]
target = "center"
set = "b"
[[dynamics.rules]]
name = "a-to-c"
center = { kind = "site", state = ["a"] }
rate = { probability = 1.0 }
[[dynamics.rules.effects]]
target = "center"
set = "c"
[execution]
strategy = "pca"
"#;
    let parsed: petra_deck::DeckFile = toml::from_str(text).expect("conflict deck parses");
    let deck = petra_deck::compile(&parsed).expect("conflict deck compiles");
    let mut engine = deck.build_engine(Some(4)).expect("conflict engine");
    let mut strategy = deck.strategy();
    let error = engine
        .step_with(&mut strategy)
        .expect_err("different simultaneous writes must fail closed");
    assert!(
        matches!(error, petra_core::Stop::EffectFailed { reason, .. } if reason == "conflicting simultaneous writes")
    );
    assert_eq!(engine.lattice.states[0].0, 0, "failure is atomic");
}

#[test]
fn stochastic_strategy_draw_order_is_pinned() {
    assert_eq!(
        trajectory_fingerprint("ising.toml", 42, 64),
        2_331_448_806_470_377_247
    );
    assert_eq!(
        trajectory_fingerprint("sir.toml", 42, 8),
        7_104_091_214_398_920_014
    );
}
