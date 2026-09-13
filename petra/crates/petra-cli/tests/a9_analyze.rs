use std::path::PathBuf;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn repo_path(relative: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(relative)
}

#[test]
fn tiny_real_cli_run_reaches_typed_zero_analysis() {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("clock after epoch")
        .as_nanos();
    let root = std::env::temp_dir().join(format!("petra-a9-real-{nonce}"));
    std::fs::create_dir(&root).expect("temporary root");
    let source = std::fs::read_to_string(repo_path("petra/examples/kaolinite-approx.toml"))
        .expect("canonical A9 deck");
    let tiny = source
        .replace("steps = 20000", "steps = 8")
        .replace("report_every = 2000", "report_every = 1");
    assert_ne!(source, tiny, "tiny execution bounds changed");
    let deck = root.join("tiny.toml");
    std::fs::write(&deck, tiny).expect("tiny deck");
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
    let python = Command::new("python3")
        .arg(repo_path("petra/scripts/approximate_rate_closure.py"))
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
    assert!(payload.contains("\"acceptance_passed\": false"));
    assert!(payload.contains("\"upper_95_mol_m2_s\":"));
    std::fs::remove_dir_all(root).expect("temporary cleanup");
}
