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
