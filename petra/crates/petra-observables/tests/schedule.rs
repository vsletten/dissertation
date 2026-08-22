use petra_observables::{run_ensemble, EnsembleConfig};

const SCHEDULE: &str = r#"
[deck]
name = "ensemble-schedule"
schema = 2
[structure]
kind = "grid"
grid = "square"
neighborhood = "von_neumann"
dims = [4, 4]
boundary = ["periodic", "periodic"]
default_kind = "S"
[[structure.kinds]]
name = "S"
initial = "a"
[[structure.kinds.states]]
name = "a"
occupant = "vacant"
[[structure.kinds.states]]
name = "b"
occupant = "vacant"
[dynamics.thermo]
temperature = 1.0
[[dynamics.rules]]
name = "flip"
center = { kind = "S", state = ["a"] }
rate = { arrhenius = { prefactor = 100.0, ea = 1000.0 } }
[[dynamics.rules.effects]]
target = "center"
set = "b"
[execution]
strategy = "ctmc"
[[execution.schedule]]
temperature = 1.0
duration = 0.75
[[execution.schedule]]
temperature = 100000.0
duration = 1.25
[execution.stop]
steps = 1000
[observables]
report_every = 1
[[observables.series]]
kind = "state_counts"
"#;

#[test]
fn ensemble_runner_executes_every_schedule_to_the_final_wall_time() {
    let parsed: petra_deck::DeckFile = toml::from_str(SCHEDULE).expect("schedule parses");
    let deck = petra_deck::compile(&parsed).expect("schedule compiles");
    let result = run_ensemble(
        &deck,
        &EnsembleConfig {
            replicas: 8,
            base_seed: 50,
            steps: 1000,
            burn_in: 0,
            sample_every: 5,
            bootstrap_resamples: 10,
            bootstrap_seed: 99,
        },
    )
    .expect("scheduled ensemble runs");

    assert_eq!(result.replicas.len(), 8);
    for replica in &result.replicas {
        let final_sample = replica.samples.last().expect("final sample");
        assert_eq!(final_sample.time.to_bits(), 2.0f64.to_bits());
        assert!(
            replica.stop.is_none(),
            "normal schedule completion is not a stop"
        );
        assert!(
            replica.final_state_counts[1] > 0,
            "T2 table rebuild must make the reaction live"
        );
    }
}
