use std::path::PathBuf;

use petra_deck::{replica_seed, SeedPolicy, StructureKind};
use rand::Rng as _;

const V1: &str = r#"
[deck]
name = "shim-equivalence"
units = "kcal/mol"

[cell]
a = 5.0
b = 5.0
c = 5.0
alpha = 90.0
beta = 90.0
gamma = 90.0
[[cell.sites]]
kind = "S"
frac = [0.0, 0.0, 0.0]

[[species]]
name = "A"
[[kinds]]
name = "S"
initial = "a"
[[kinds.states]]
name = "a"
occupant = "A"
[[kinds.states]]
name = "b"
occupant = "A"

[lattice]
dims = [2, 2, 2]
boundary = ["periodic", "periodic", "periodic"]
[thermo]
temperature = 300.0

[[reactions]]
name = "flip"
center = { kind = "S", state = ["a"] }
rate = { constant = 2.0 }
[[reactions.effects]]
target = "center"
set = "b"

[simulation]
steps = 20
seed = 42
report_every = 5
"#;

const V2: &str = r#"
[deck]
name = "shim-equivalence"
schema = 2
units = "kcal/mol"

[structure]
kind = "cell"
[[structure.species]]
name = "A"
[[structure.kinds]]
name = "S"
initial = "a"
[[structure.kinds.states]]
name = "a"
occupant = "A"
[[structure.kinds.states]]
name = "b"
occupant = "A"

[structure.cell]
a = 5.0
b = 5.0
c = 5.0
alpha = 90.0
beta = 90.0
gamma = 90.0
[[structure.cell.sites]]
kind = "S"
frac = [0.0, 0.0, 0.0]

[structure.lattice]
dims = [2, 2, 2]
boundary = ["periodic", "periodic", "periodic"]

[dynamics.thermo]
temperature = 300.0
[[dynamics.rules]]
name = "flip"
center = { kind = "S", state = ["a"] }
rate = { constant = 2.0 }
[[dynamics.rules.effects]]
target = "center"
set = "b"

[execution]
strategy = "ctmc"
[execution.stop]
steps = 20
[execution.ensemble]
seed = 42
n_replicas = 1
seed_policy = "increment"

[observables]
report_every = 5
"#;

fn repo_path(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

fn explicit_v2_from_v1(text: &str) -> String {
    let mut value: toml::Value = text.parse().expect("v1 TOML value");
    let root = value.as_table_mut().expect("deck root table");
    root.get_mut("deck")
        .and_then(toml::Value::as_table_mut)
        .expect("deck metadata")
        .insert("schema".to_string(), toml::Value::Integer(2));

    let mut structure = toml::map::Map::new();
    structure.insert("kind".to_string(), toml::Value::String("cell".to_string()));
    for key in ["cell", "lattice", "species", "kinds", "init", "defects"] {
        if let Some(section) = root.remove(key) {
            structure.insert(key.to_string(), section);
        }
    }

    let mut dynamics = toml::map::Map::new();
    for (v1, v2) in [
        ("thermo", "thermo"),
        ("aliases", "aliases"),
        ("reactions", "rules"),
    ] {
        if let Some(section) = root.remove(v1) {
            dynamics.insert(v2.to_string(), section);
        }
    }

    let simulation = root
        .remove("simulation")
        .and_then(|value| value.as_table().cloned())
        .expect("v1 simulation");
    let mut stop = toml::map::Map::new();
    stop.insert("steps".to_string(), simulation["steps"].clone());
    let mut ensemble = toml::map::Map::new();
    ensemble.insert("seed".to_string(), simulation["seed"].clone());
    ensemble.insert("n_replicas".to_string(), toml::Value::Integer(1));
    ensemble.insert(
        "seed_policy".to_string(),
        toml::Value::String("increment".to_string()),
    );
    let mut execution = toml::map::Map::new();
    execution.insert(
        "strategy".to_string(),
        toml::Value::String("ctmc".to_string()),
    );
    execution.insert("stop".to_string(), toml::Value::Table(stop));
    execution.insert("ensemble".to_string(), toml::Value::Table(ensemble));
    let mut observables = toml::map::Map::new();
    observables.insert(
        "report_every".to_string(),
        simulation
            .get("report_every")
            .cloned()
            .unwrap_or(toml::Value::Integer(0)),
    );

    root.insert("structure".to_string(), toml::Value::Table(structure));
    root.insert("dynamics".to_string(), toml::Value::Table(dynamics));
    root.insert("execution".to_string(), toml::Value::Table(execution));
    root.insert("observables".to_string(), toml::Value::Table(observables));
    toml::to_string(&value).expect("serialize explicit v2")
}

#[test]
fn v1_without_schema_normalizes_to_ctmc_v2() {
    let deck: petra_deck::DeckFile = toml::from_str(V1).expect("v1 parses through shim");
    assert_eq!(deck.deck.schema, Some(2));
    assert_eq!(deck.structure_kind, StructureKind::Cell);
    assert_eq!(deck.execution.strategy, "ctmc");
    assert_eq!(deck.execution.stop.steps, Some(20));
    assert_eq!(deck.execution.ensemble.seed, 42);
    assert_eq!(deck.observables.report_every, 5);
}

#[test]
fn explicit_cell_v2_compiles() {
    let parsed: petra_deck::DeckFile = toml::from_str(V2).expect("v2 parses");
    let deck = petra_deck::compile(&parsed).expect("v2 compiles");
    assert_eq!(deck.steps, 20);
    assert_eq!(deck.seed, 42);
    assert_eq!(deck.report_every, 5);
    assert_eq!(deck.reactions.len(), 1);
}

#[test]
fn v1_and_equivalent_v2_compile_to_identical_dense_tables_and_trajectory() {
    let v1: petra_deck::DeckFile = toml::from_str(V1).expect("v1 parses");
    let v2: petra_deck::DeckFile = toml::from_str(V2).expect("v2 parses");
    let a = petra_deck::compile(&v1).expect("v1 compiles");
    let b = petra_deck::compile(&v2).expect("v2 compiles");

    assert_eq!(a.kind_names, b.kind_names);
    assert_eq!(a.state_names, b.state_names);
    assert_eq!(a.kind_state_ranges, b.kind_state_ranges);
    assert_eq!(a.initial_per_template, b.initial_per_template);
    assert_eq!(a.dims, b.dims);
    assert_eq!(a.boundary, b.boundary);
    assert_eq!(a.steps, b.steps);
    assert_eq!(a.seed, b.seed);
    assert_eq!(a.report_every, b.report_every);

    let mut ea = a.build_engine(None).expect("v1 engine");
    let mut eb = b.build_engine(None).expect("v2 engine");
    for _ in 0..8 {
        let fa = ea.step().expect("v1 fires");
        let fb = eb.step().expect("v2 fires");
        assert_eq!(fa.step, fb.step);
        assert_eq!(fa.time.to_bits(), fb.time.to_bits());
        assert_eq!(fa.site, fb.site);
        assert_eq!(fa.reaction, fb.reaction);
    }
}

#[test]
fn every_shipped_v1_deck_matches_explicit_v2_in_every_compiled_table() {
    for name in [
        "kaolinite.toml",
        "kaolinite-multilayer.toml",
        "kossel.toml",
        "kossel-etchpit.toml",
    ] {
        let path = repo_path(&format!("petra/examples/{name}"));
        let text = std::fs::read_to_string(&path).expect("read shipped deck");
        let shimmed: petra_deck::DeckFile =
            toml::from_str(&text).unwrap_or_else(|error| panic!("{name} shim parse: {error}"));
        let explicit_text = explicit_v2_from_v1(&text);
        let explicit: petra_deck::DeckFile = toml::from_str(&explicit_text)
            .unwrap_or_else(|error| panic!("{name} explicit v2 parse: {error}"));
        let shimmed = petra_deck::compile(&shimmed)
            .unwrap_or_else(|error| panic!("{name} shim compiles: {error}"));
        let explicit = petra_deck::compile(&explicit)
            .unwrap_or_else(|error| panic!("{name} explicit v2 compiles: {error}"));
        assert_eq!(shimmed, explicit, "all compiled tables for {name}");
    }
}

#[test]
fn schema_2_rejects_flat_v1_keys() {
    let mixed = V1.replace("units = \"kcal/mol\"", "schema = 2\nunits = \"kcal/mol\"");
    let error = toml::from_str::<petra_deck::DeckFile>(&mixed).expect_err("mixed schema rejected");
    assert!(error.to_string().contains("unknown field"), "{error}");
}

#[test]
fn unknown_schema_version_is_rejected() {
    let bad = V2.replace("schema = 2", "schema = 99");
    let error = toml::from_str::<petra_deck::DeckFile>(&bad).expect_err("unknown schema rejected");
    assert!(
        error.to_string().contains("unsupported deck schema 99"),
        "{error}"
    );
}

#[test]
fn metropolis_strategy_compiles_in_b3() {
    let text = V2
        .replace("strategy = \"ctmc\"", "strategy = \"metropolis\"")
        .replace(
            "rate = { constant = 2.0 }",
            "rate = { energy = { delta = 2.0 } }",
        );
    let parsed: petra_deck::DeckFile = toml::from_str(&text).expect("future strategy parses");
    let deck = petra_deck::compile(&parsed).expect("Metropolis compiles");
    assert!(matches!(
        deck.strategy,
        petra_deck::ExecutionStrategy::Metropolis { temperature: 300.0 }
    ));
}

#[test]
fn every_rfc_strategy_parameter_block_parses() {
    let text = V2.replace(
        "[execution.stop]",
        "[execution.ctmc]\n[execution.synchronous]\nconflict_resolution = \"first_match\"\n[execution.metropolis]\ntemperature = 310.0\n[execution.pca]\n[execution.stop]",
    );
    let parsed: petra_deck::DeckFile = toml::from_str(&text).expect("all blocks parse");
    assert!(parsed.execution.ctmc.is_some());
    let synchronous = parsed
        .execution
        .synchronous
        .as_ref()
        .expect("synchronous block present");
    assert_eq!(
        synchronous.conflict_resolution,
        petra_deck::schema::SynchronousConflictResolution::FirstMatch
    );
    assert_eq!(
        parsed
            .execution
            .metropolis
            .as_ref()
            .and_then(|spec| spec.temperature),
        Some(310.0)
    );
    assert!(parsed.execution.pca.is_some());
}

#[test]
fn execution_parameter_blocks_validate_at_parse_time() {
    let invalid_conflict_resolution = V2.replace(
        "[execution.stop]",
        "[execution.synchronous]\nconflict_resolution = \"last_match\"\n[execution.stop]",
    );
    let error = toml::from_str::<petra_deck::DeckFile>(&invalid_conflict_resolution)
        .expect_err("unknown conflict resolution rejected");
    assert!(error.to_string().contains("first_match"), "{error}");

    let bad_temperature = V2.replace(
        "[execution.stop]",
        "[execution.metropolis]\ntemperature = -1.0\n[execution.stop]",
    );
    let error = toml::from_str::<petra_deck::DeckFile>(&bad_temperature)
        .expect_err("negative temperature rejected");
    assert!(error.to_string().contains("finite and positive"), "{error}");

    let no_replicas = V2.replace("n_replicas = 1", "n_replicas = 0");
    let error =
        toml::from_str::<petra_deck::DeckFile>(&no_replicas).expect_err("zero replicas rejected");
    assert!(error.to_string().contains("at least 1"), "{error}");
}

#[test]
fn b3_strategy_rate_surfaces_parse_and_compile() {
    let synchronous = V2
        .replace("strategy = \"ctmc\"", "strategy = \"synchronous\"")
        .replace("rate = { constant = 2.0 }\n", "truth = \"and\"\n")
        .replace(
            "[execution.stop]",
            "[execution.synchronous]\nconflict_resolution = \"first_match\"\n[execution.stop]",
        );
    let parsed: petra_deck::DeckFile =
        toml::from_str(&synchronous).expect("synchronous truth surface parses");
    assert!(petra_deck::compile(&parsed).is_ok());

    let pca = V2
        .replace("strategy = \"ctmc\"", "strategy = \"pca\"")
        .replace("rate = { constant = 2.0 }", "rate = { probability = 0.25 }");
    let parsed: petra_deck::DeckFile = toml::from_str(&pca).expect("PCA probability parses");
    assert!(petra_deck::compile(&parsed).is_ok());
}

#[test]
fn rfc_surfaces_fail_closed_at_parse_time() {
    let bad_mode = V2.replace(
        "target = \"center\"",
        "target = \"center\"\nselect_mode = \"one\"",
    );
    assert!(toml::from_str::<petra_deck::DeckFile>(&bad_mode).is_err());

    let bad_probability = V2
        .replace("strategy = \"ctmc\"", "strategy = \"pca\"")
        .replace("rate = { constant = 2.0 }", "rate = { probability = 1.5 }");
    assert!(toml::from_str::<petra_deck::DeckFile>(&bad_probability).is_err());

    let bad_observable = V2.replace(
        "report_every = 5",
        "report_every = 5\n[[observables.series]]\nkind = \"interface_roughness\"",
    );
    assert!(toml::from_str::<petra_deck::DeckFile>(&bad_observable).is_err());
}

#[test]
fn grid_v2_compiles_to_the_shared_lattice_representation() {
    let text = r#"
[deck]
name = "grid-parse"
schema = 2
[structure]
kind = "grid"
grid = "square"
neighborhood = "moore"
dims = [8, 8]
boundary = ["periodic", "periodic"]
default_kind = "S"
[[structure.kinds]]
name = "S"
initial = "a"
[[structure.kinds.states]]
name = "a"
occupant = "vacant"
[dynamics.thermo]
temperature = 300.0
[execution]
strategy = "ctmc"
[execution.stop]
steps = 10
[observables]
"#;
    let parsed: petra_deck::DeckFile = toml::from_str(text).expect("grid schema parses");
    assert_eq!(parsed.structure_kind, StructureKind::Grid);
    let deck = petra_deck::compile(&parsed).expect("grid compiles in B3");
    let engine = deck.build_engine(Some(1)).expect("grid lattice builds");
    assert_eq!(engine.lattice.dims, [8, 8, 1]);
    assert!(engine
        .lattice
        .adj_off
        .windows(2)
        .all(|offset| offset[1] - offset[0] == 8));
}

struct PublicSeamCtmc;

impl petra_core::UpdateStrategy for PublicSeamCtmc {
    fn step(
        &mut self,
        ctx: &mut petra_core::StepCtx<'_>,
    ) -> Result<petra_core::StepOutcome, petra_core::Stop> {
        assert_eq!(ctx.lattice.len(), 8);
        assert_eq!(ctx.rules.len(), 1);
        let total = ctx.apply.total_rate();
        let draw = ctx.rng.gen::<f64>() * total;
        let (site, reaction) = ctx
            .apply
            .select_event(draw)
            .ok_or(petra_core::Stop::NoEvents)?;
        let wait_draw = ctx.rng.gen::<f64>();
        let dt = -(1.0 - wait_draw).ln() / total;
        ctx.apply.apply_transition(site, reaction, ctx.rng)?;
        Ok(petra_core::StepOutcome {
            fired: vec![petra_core::Fired {
                step: 0,
                time: 0.0,
                site,
                reaction,
            }],
            dt,
        })
    }
}

#[test]
fn step_context_public_fields_are_a_complete_strategy_seam() {
    let parsed: petra_deck::DeckFile = toml::from_str(V2).expect("v2 parses");
    let deck = petra_deck::compile(&parsed).expect("v2 compiles");
    let mut engine = deck.build_engine(None).expect("engine builds");
    let outcome = engine
        .step_with(&mut PublicSeamCtmc)
        .expect("public-only strategy advances");
    assert_eq!(outcome.fired[0].step, 1);
    assert_eq!(outcome.fired[0].time.to_bits(), outcome.dt.to_bits());
}

#[test]
fn explicit_exact_ctmc_strategy_matches_compatibility_wrapper() {
    let parsed: petra_deck::DeckFile = toml::from_str(V1).expect("v1 parses");
    let deck = petra_deck::compile(&parsed).expect("compiles");
    let mut explicit = deck.build_engine(None).expect("explicit engine");
    let mut compatibility = deck.build_engine(None).expect("compat engine");

    let mut strategy = petra_core::ExactCtmc;
    let outcome = explicit.step_with(&mut strategy).expect("strategy step");
    let fired = compatibility.step().expect("compatibility step");
    assert_eq!(outcome.fired.len(), 1);
    assert_eq!(outcome.dt.to_bits(), outcome.fired[0].time.to_bits());
    assert_eq!(outcome.fired[0].step, fired.step);
    assert_eq!(outcome.fired[0].time.to_bits(), fired.time.to_bits());
    assert_eq!(outcome.fired[0].site, fired.site);
    assert_eq!(outcome.fired[0].reaction, fired.reaction);
}

#[test]
fn select_mode_one_changes_exactly_one_matching_neighbor() {
    let text = V2
        .replace("dims = [2, 2, 2]", "dims = [3, 1, 1]")
        .replace(
            "frac = [0.0, 0.0, 0.0]",
            "frac = [0.0, 0.0, 0.0]\n[[structure.cell.bonds]]\ni = 0\nj = 0\ndcell = [1, 0, 0]",
        )
        .replace(
            "center = { kind = \"S\", state = [\"a\"] }",
            "center = { kind = \"S\", state = [\"b\"] }",
        )
        .replace(
            "initial = \"a\"",
            "initial = \"a\"\n[[structure.init]]\nname = \"seed-center\"\ncenter = { kind = \"S\", state = [\"a\"] }\nsites = [[0, 0, 0, 0]]\nset = \"b\"",
        )
        .replace(
            "target = \"center\"\nset = \"b\"",
            "target = \"neighbors\"\nselect = { distance = 1, kind = \"S\", state = [\"a\"] }\nselect_mode = \"one\"\nset = \"b\"",
        );
    let parsed: petra_deck::DeckFile = toml::from_str(&text).expect("select_mode deck parses");
    let deck = petra_deck::compile(&parsed).expect("select_mode deck compiles");
    let mut engine = deck.build_engine(Some(17)).expect("engine builds");
    engine.step().expect("one-neighbor event fires");
    let b = deck
        .state_names
        .iter()
        .position(|name| name == "S.b")
        .expect("state exists") as u16;
    assert_eq!(
        engine
            .lattice
            .states
            .iter()
            .filter(|state| state.0 == b)
            .count(),
        2,
        "the center plus exactly one of its two neighbors changes"
    );
}

#[test]
fn source_rule_is_one_independent_event_per_vacant_site() {
    let text = V2
        .replace("dims = [2, 2, 2]", "dims = [4, 1, 1]")
        .replace(
            "name = \"a\"\noccupant = \"A\"",
            "name = \"a\"\noccupant = \"vacant\"",
        )
        .replace("center = { kind = \"S\", state = [\"a\"] }\n", "")
        .replace("target = \"center\"", "target = \"source\"");
    let parsed: petra_deck::DeckFile = toml::from_str(&text).expect("source deck parses");
    let deck = petra_deck::compile(&parsed).expect("source deck compiles");
    let mut engine = deck.build_engine(Some(9)).expect("engine builds");
    let mut fired = 0;
    while engine.step().is_ok() {
        fired += 1;
    }
    assert_eq!(fired, 4, "one source event for each initially vacant site");
}

#[test]
fn init_rekind_changes_runtime_site_kind_and_state() {
    let text = V2
        .replace(
            "[[structure.kinds.states]]\nname = \"b\"\noccupant = \"A\"",
            "[[structure.kinds.states]]\nname = \"b\"\noccupant = \"A\"\n\n[[structure.kinds]]\nname = \"T\"\ninitial = \"a\"\n[[structure.kinds.states]]\nname = \"a\"\noccupant = \"A\"\n[[structure.kinds.states]]\nname = \"b\"\noccupant = \"A\"",
        )
        .replace(
            "center = { kind = \"S\", state = [\"a\"] }",
            "center = { kind = \"T\", state = [\"a\"] }",
        )
        .replace(
            "[structure.cell]",
            "[[structure.init]]\nname = \"rekind-all\"\ncenter = { kind = \"S\", state = [\"a\"] }\nmap = { a = \"a\" }\nrekind = \"T\"\n\n[structure.cell]",
        );
    let parsed: petra_deck::DeckFile = toml::from_str(&text).expect("rekind deck parses");
    let deck = petra_deck::compile(&parsed).expect("rekind deck compiles");
    let mut engine = deck.build_engine(Some(3)).expect("rekind init runs");
    engine.step().expect("T-kind reaction sees rekindled sites");
    let t_b = deck
        .state_names
        .iter()
        .position(|name| name == "T.b")
        .expect("T.b exists") as u16;
    assert!(engine.lattice.states.iter().any(|state| state.0 == t_b));
}

#[test]
fn splitmix64_replica_seed_matches_rfc_vectors() {
    assert_eq!(replica_seed(42, 0, SeedPolicy::Hash), 0x57E1_FABA_6510_7204);
    assert_eq!(replica_seed(42, 1, SeedPolicy::Hash), 0xF34F_E924_8C93_42E5);
    assert_eq!(replica_seed(42, 2, SeedPolicy::Hash), 0x7253_9538_8690_AE46);
    assert_eq!(replica_seed(42, 3, SeedPolicy::Increment), 45);
}
