use std::path::PathBuf;

use petra_deck::{replica_seed, SeedPolicy, StructureKind};

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
fn every_shipped_v1_deck_shims_and_compiles() {
    for name in [
        "kaolinite.toml",
        "kaolinite-multilayer.toml",
        "kossel.toml",
        "kossel-etchpit.toml",
    ] {
        let path = repo_path(&format!("petra/examples/{name}"));
        let text = std::fs::read_to_string(&path).expect("read shipped deck");
        let parsed: petra_deck::DeckFile =
            toml::from_str(&text).unwrap_or_else(|e| panic!("{name} shim parse: {e}"));
        assert_eq!(parsed.deck.schema, Some(2), "{name}");
        petra_deck::compile(&parsed).unwrap_or_else(|e| panic!("{name} compiles: {e}"));
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
fn non_ctmc_strategy_is_rejected_in_b2() {
    let text = V2.replace("strategy = \"ctmc\"", "strategy = \"metropolis\"");
    let parsed: petra_deck::DeckFile = toml::from_str(&text).expect("future strategy parses");
    let error = petra_deck::compile(&parsed).expect_err("B2 implements only CTMC");
    assert!(
        error.to_string().contains("not implemented in B2"),
        "{error}"
    );
}

#[test]
fn grid_v2_parses_but_waits_for_b3_compilation() {
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
    let error = petra_deck::compile(&parsed).expect_err("grid compile is B3");
    assert!(error.to_string().contains("B3"), "{error}");
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
fn splitmix64_replica_seed_matches_rfc_vectors() {
    assert_eq!(replica_seed(42, 0, SeedPolicy::Hash), 0x57E1_FABA_6510_7204);
    assert_eq!(replica_seed(42, 1, SeedPolicy::Hash), 0xF34F_E924_8C93_42E5);
    assert_eq!(replica_seed(42, 2, SeedPolicy::Hash), 0x7253_9538_8690_AE46);
    assert_eq!(replica_seed(42, 3, SeedPolicy::Increment), 45);
}
