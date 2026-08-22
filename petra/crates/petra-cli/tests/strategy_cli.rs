use std::path::PathBuf;
use std::process::Command;

fn write_schedule_deck(
    root: &std::path::Path,
    steps: Option<u64>,
    report_every: Option<u64>,
) -> PathBuf {
    std::fs::create_dir_all(root).expect("temp dir");
    let deck = root.join("schedule.toml");
    let stop = steps.map_or_else(String::new, |steps| format!("steps = {steps}"));
    let cadence = report_every.map_or_else(String::new, |value| format!("report_every = {value}"));
    std::fs::write(
        &deck,
        format!(
            r#"
[deck]
name = "cli-schedule-gates"
schema = 2
[structure]
kind = "grid"
grid = "square"
neighborhood = "von_neumann"
dims = [2, 2]
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
temperature = 100000.0
[[dynamics.rules]]
name = "flip"
center = {{ kind = "S", state = ["a"] }}
rate = {{ arrhenius = {{ prefactor = 100.0, ea = 1000.0 }} }}
[[dynamics.rules.effects]]
target = "center"
set = "b"
[execution]
strategy = "ctmc"
[[execution.schedule]]
temperature = 100000.0
duration = 10.0
[[execution.schedule]]
temperature = 100000.0
duration = 1.0
[execution.stop]
{stop}
[observables]
{cadence}
"#
        ),
    )
    .expect("write temp deck");
    deck
}

#[test]
fn cli_dispatches_the_deck_selected_strategy() {
    let repo = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..");
    let deck = repo.join("petra/examples/conformance/conway.toml");
    let out = std::env::temp_dir().join(format!("petra-cli-strategy-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&out);

    let output = Command::new(env!("CARGO_BIN_EXE_petra"))
        .arg(deck)
        .arg("--steps")
        .arg("4")
        .arg("--out")
        .arg(&out)
        .output()
        .expect("petra CLI runs");
    let _ = std::fs::remove_dir_all(&out);
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("strategy = synchronous"), "{stdout}");
    assert!(stdout.contains("completed 4 steps"), "{stdout}");
}

#[test]
fn cli_rejects_batch_trajectory_output_for_discrete_strategies() {
    let repo = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..");
    let deck = repo.join("petra/examples/conformance/conway.toml");
    let out = std::env::temp_dir().join(format!("petra-cli-viz-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&out);

    let output = Command::new(env!("CARGO_BIN_EXE_petra"))
        .arg(deck)
        .arg("--steps")
        .arg("1")
        .arg("--viz")
        .arg("--out")
        .arg(&out)
        .output()
        .expect("petra CLI runs");
    let _ = std::fs::remove_dir_all(&out);
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("batch trajectory format"), "{stderr}");
}

#[test]
fn cli_parallel_ensemble_writes_distribution_ci_and_declared_observables() {
    let repo = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..");
    let source = repo.join("petra/examples/conformance/ising.toml");
    let root = std::env::temp_dir().join(format!("petra-cli-ensemble-{}", std::process::id()));
    let deck = root.join("ising.toml");
    let out = root.join("out");
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(&root).expect("temp dir");
    let text = std::fs::read_to_string(source).expect("fixture");
    std::fs::write(&deck, text).expect("temp deck");
    let output = Command::new(env!("CARGO_BIN_EXE_petra"))
        .arg(&deck)
        .arg("--steps")
        .arg("20")
        .arg("--ensemble")
        .arg("4")
        .arg("--out")
        .arg(&out)
        .output()
        .expect("petra CLI runs");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let summary = std::fs::read_to_string(out.join("ensemble-summary.csv"))
        .expect("distribution summary exists");
    assert!(
        summary.starts_with("state,mean,ci95_low,ci95_high,distribution\n"),
        "{summary}"
    );
    let observables = std::fs::read_to_string(out.join("observables.csv"))
        .expect("declarative observable output exists");
    assert!(
        observables.starts_with("replica,seed,step,time,kind,index,value\n"),
        "{observables}"
    );
    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn cli_single_run_writes_all_declared_kossel_observables() {
    let repo = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..");
    let deck = repo.join("petra/examples/kossel-etchpit.toml");
    let root = std::env::temp_dir().join(format!("petra-cli-observables-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&root);
    let output = Command::new(env!("CARGO_BIN_EXE_petra"))
        .arg(deck)
        .arg("--steps")
        .arg("200")
        .arg("--out")
        .arg(&root)
        .output()
        .expect("petra CLI runs");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let observables = std::fs::read_to_string(root.join("observables.csv"))
        .expect("single-run declarative observable output exists");
    for kind in [
        "state_counts",
        "event_rates",
        "rate_spectra",
        "cluster_sizes",
        "surface_area",
    ] {
        assert!(
            observables.contains(&format!(",{kind},")),
            "missing {kind}: {observables}"
        );
    }
    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn cli_replica_zero_is_invariant_when_hash_ensemble_cardinality_changes() {
    let repo = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..");
    let deck = repo.join("petra/examples/conformance/ising.toml");
    let root = std::env::temp_dir().join(format!("petra-cli-seed-policy-{}", std::process::id()));
    let single_out = root.join("single");
    let ensemble_out = root.join("ensemble");
    let _ = std::fs::remove_dir_all(&root);

    let run = |replicas: &str, out: &std::path::Path| {
        Command::new(env!("CARGO_BIN_EXE_petra"))
            .arg(&deck)
            .arg("--seed")
            .arg("99")
            .arg("--steps")
            .arg("20")
            .arg("--ensemble")
            .arg(replicas)
            .arg("--out")
            .arg(out)
            .output()
            .expect("petra CLI runs")
    };
    let single = run("1", &single_out);
    let ensemble = run("2", &ensemble_out);
    assert!(
        single.status.success(),
        "single stderr: {}",
        String::from_utf8_lossy(&single.stderr)
    );
    assert!(
        ensemble.status.success(),
        "ensemble stderr: {}",
        String::from_utf8_lossy(&ensemble.stderr)
    );

    let populations = std::fs::read_to_string(single_out.join("populations.csv"))
        .expect("single populations exist");
    let single_fields: Vec<_> = populations
        .lines()
        .last()
        .expect("final row")
        .split(',')
        .collect();
    let ensemble_csv = std::fs::read_to_string(ensemble_out.join("ensemble.csv"))
        .expect("ensemble populations exist");
    let replica_zero_fields: Vec<_> = ensemble_csv
        .lines()
        .nth(1)
        .expect("replica zero row")
        .split(',')
        .collect();
    assert_eq!(
        &single_fields[2..],
        &replica_zero_fields[3..],
        "replica zero trajectory must not depend on ensemble cardinality"
    );

    let expected_seed = petra_deck::replica_seed(99, 0, petra_deck::SeedPolicy::Hash).to_string();
    assert_eq!(replica_zero_fields[0], expected_seed);
    let observables = std::fs::read_to_string(single_out.join("observables.csv"))
        .expect("single observables exist");
    let reported_seed = observables
        .lines()
        .nth(1)
        .expect("observable data row")
        .split(',')
        .nth(1)
        .expect("seed field");
    assert_eq!(reported_seed, expected_seed);

    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn cli_runs_piecewise_isothermal_schedule_to_its_wall_time() {
    let root = std::env::temp_dir().join(format!("petra-cli-schedule-{}", std::process::id()));
    let deck = root.join("schedule.toml");
    let out = root.join("out");
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(&root).expect("temp dir");
    std::fs::write(
        &deck,
        r#"
[deck]
name = "cli-schedule"
schema = 2
[structure]
kind = "grid"
grid = "square"
neighborhood = "von_neumann"
dims = [2, 2]
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
rate = { arrhenius = { prefactor = 1.0, ea = 1000.0 } }
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
"#,
    )
    .expect("temp deck");

    let output = Command::new(env!("CARGO_BIN_EXE_petra"))
        .arg(&deck)
        .arg("--out")
        .arg(&out)
        .output()
        .expect("petra CLI runs");
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert!(stdout.contains("completed schedule"), "{stdout}");
    let populations =
        std::fs::read_to_string(out.join("populations.csv")).expect("populations exist");
    let final_time: f64 = populations
        .lines()
        .last()
        .expect("final row")
        .split(',')
        .nth(1)
        .expect("time column")
        .parse()
        .expect("numeric time");
    assert_eq!(final_time.to_bits(), 2.0f64.to_bits());
    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn scheduled_cli_honors_an_explicit_zero_event_cap() {
    let root = std::env::temp_dir().join(format!("petra-cli-zero-cap-{}", std::process::id()));
    let deck = write_schedule_deck(&root, Some(0), Some(100));
    let out = root.join("out");
    let output = Command::new(env!("CARGO_BIN_EXE_petra"))
        .arg(&deck)
        .arg("--out")
        .arg(&out)
        .output()
        .expect("petra CLI runs");
    assert!(output.status.success());
    let populations = std::fs::read_to_string(out.join("populations.csv")).expect("populations");
    let fields: Vec<_> = populations
        .lines()
        .last()
        .expect("final row")
        .split(',')
        .collect();
    assert_eq!(fields[0], "0");
    assert_eq!(
        fields[1].parse::<f64>().expect("time").to_bits(),
        0.0f64.to_bits()
    );
    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn scheduled_cli_reports_the_final_state_when_cap_falls_between_cadences() {
    let root = std::env::temp_dir().join(format!("petra-cli-final-cap-{}", std::process::id()));
    let deck = write_schedule_deck(&root, Some(1), Some(100));
    let out = root.join("out");
    let output = Command::new(env!("CARGO_BIN_EXE_petra"))
        .arg(&deck)
        .arg("--seed")
        .arg("7")
        .arg("--out")
        .arg(&out)
        .output()
        .expect("petra CLI runs");
    assert!(output.status.success());
    let populations = std::fs::read_to_string(out.join("populations.csv")).expect("populations");
    let fields: Vec<_> = populations
        .lines()
        .last()
        .expect("final row")
        .split(',')
        .collect();
    assert_eq!(fields[0], "1", "final event-cap state must be reported");
    assert!(fields[1].parse::<f64>().expect("time") > 0.0);
    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn scheduled_cli_without_cadence_reports_boundaries_not_every_event() {
    let root = std::env::temp_dir().join(format!("petra-cli-no-cadence-{}", std::process::id()));
    let deck = write_schedule_deck(&root, None, None);
    let out = root.join("out");
    let output = Command::new(env!("CARGO_BIN_EXE_petra"))
        .arg(&deck)
        .arg("--out")
        .arg(&out)
        .output()
        .expect("petra CLI runs");
    assert!(output.status.success());
    let populations = std::fs::read_to_string(out.join("populations.csv")).expect("populations");
    assert!(
        populations.lines().count() <= 5,
        "default cadence must not write every event: {populations}"
    );
    let final_time: f64 = populations
        .lines()
        .last()
        .expect("final row")
        .split(',')
        .nth(1)
        .expect("time")
        .parse()
        .expect("numeric time");
    assert_eq!(final_time.to_bits(), 11.0f64.to_bits());
    let _ = std::fs::remove_dir_all(&root);
}

#[test]
fn scheduled_cli_rejects_viz_until_boundaries_have_an_artifact_encoding() {
    let root = std::env::temp_dir().join(format!("petra-cli-schedule-viz-{}", std::process::id()));
    let deck = write_schedule_deck(&root, None, Some(1));
    let output = Command::new(env!("CARGO_BIN_EXE_petra"))
        .arg(&deck)
        .arg("--viz")
        .arg("--out")
        .arg(root.join("out"))
        .output()
        .expect("petra CLI runs");
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("schedule"), "{stderr}");
    let _ = std::fs::remove_dir_all(&root);
}
