use std::path::PathBuf;
use std::process::Command;

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
