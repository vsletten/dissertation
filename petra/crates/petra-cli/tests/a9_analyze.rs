use std::path::PathBuf;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn repo_path(relative: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(relative)
}

#[test]
fn canonical_real_cli_run_reaches_typed_zero_analysis() {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock after epoch")
        .as_nanos();
    let root = std::env::temp_dir().join(format!("petra-a9-real-{nonce}"));
    std::fs::create_dir(&root).expect("temporary root");
    let deck = repo_path("petra/examples/kaolinite-approx.toml");
    let run = root.join("run");
    let petra = Command::new(env!("CARGO_BIN_EXE_petra"))
        .arg(&deck)
        .args(["--seed", "90401", "--ensemble", "1", "--out"])
        .arg(&run)
        .args(["--viz", "--paranoid"])
        .output()
        .expect("run Petra CLI");
    assert!(
        petra.status.success(),
        "petra stdout={} stderr={}",
        String::from_utf8_lossy(&petra.stdout),
        String::from_utf8_lossy(&petra.stderr)
    );

    let analysis = root.join("single-run-analysis.json");
    let python = Command::new(repo_path("petra/scripts/approximate_rate_closure.py"))
        .args(["analyze-run"])
        .arg(&deck)
        .arg(&run)
        .arg(&analysis)
        .args(["--seed", "90401"])
        .output()
        .expect("run A9 analyzer");
    assert!(
        python.status.success(),
        "analyzer stdout={} stderr={}",
        String::from_utf8_lossy(&python.stdout),
        String::from_utf8_lossy(&python.stderr)
    );
    let payload = std::fs::read_to_string(&analysis).expect("analysis payload");
    assert!(payload.contains("\"steady_state_status\": \"steady-zero\""));
    assert!(payload.contains("\"outcome\": \"no-dissolution\""));
    assert!(payload.contains("\"acceptance_passed\": true"));
    assert!(payload.contains("\"upper_95_mol_m2_s\":"));
    assert!(payload.contains("\"propensity_estimator_basis\": \"integrated_ctmc_hazard\""));
    assert!(payload.contains("\"expected_lattice_origin_si_flux_from_propensity_mol_m2_s\":"));
    assert!(payload.contains("not observed event release"));
    std::fs::remove_dir_all(root).expect("temporary cleanup");
}
