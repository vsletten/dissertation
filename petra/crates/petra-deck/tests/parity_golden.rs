//! B2's migration-grade parity gate.
//!
//! The checked-in byte logs were captured from the pre-refactor engine on
//! main (4522ecc). Each record is the complete public `Fired` value encoded
//! without textual float conversion: step/time bits/site/reaction. The
//! ExactCtmc refactor must reproduce every byte for both shipped decks.

use std::path::{Path, PathBuf};

const SEED: u64 = 42;
const STEPS: u64 = 20_000;
const MAGIC: &[u8] = b"PETRA-FIRED-V1\0";
const RECORD_BYTES: usize = 8 + 8 + 8 + 2;

fn repo_path(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

fn fired_log(deck_path: &Path, dims_override: Option<[usize; 3]>) -> Vec<u8> {
    let text = std::fs::read_to_string(deck_path).expect("read parity deck");
    let mut parsed: petra_deck::DeckFile = toml::from_str(&text).expect("parity deck parses");
    if let Some(dims) = dims_override {
        parsed.lattice.dims = dims;
    }
    let deck = petra_deck::compile(&parsed).expect("parity deck compiles");
    let mut engine = deck.build_engine(Some(SEED)).expect("engine builds");
    let mut bytes = Vec::with_capacity(MAGIC.len() + 16 + STEPS as usize * RECORD_BYTES);
    bytes.extend_from_slice(MAGIC);
    bytes.extend_from_slice(&SEED.to_le_bytes());
    bytes.extend_from_slice(&STEPS.to_le_bytes());

    for expected_step in 1..=STEPS {
        let fired = engine
            .step()
            .unwrap_or_else(|e| panic!("deck stopped at step {expected_step}: {e}"));
        bytes.extend_from_slice(&fired.step.to_le_bytes());
        bytes.extend_from_slice(&fired.time.to_bits().to_le_bytes());
        bytes.extend_from_slice(&(fired.site as u64).to_le_bytes());
        bytes.extend_from_slice(&fired.reaction.to_le_bytes());
    }
    bytes
}

fn assert_or_update(deck: &str, golden: &str, dims_override: Option<[usize; 3]>) {
    let actual = fired_log(&repo_path(deck), dims_override);
    let golden_path = repo_path(golden);
    if std::env::var_os("UPDATE_PETRA_PARITY_GOLDENS").is_some() {
        std::fs::create_dir_all(golden_path.parent().expect("golden has parent"))
            .expect("create golden directory");
        std::fs::write(&golden_path, &actual).expect("write parity golden");
        return;
    }

    let expected = std::fs::read(&golden_path).unwrap_or_else(|e| {
        panic!(
            "cannot read {}: {e}; regenerate only from the pinned pre-refactor baseline",
            golden_path.display()
        )
    });
    if actual != expected {
        let first = actual
            .iter()
            .zip(&expected)
            .position(|(a, b)| a != b)
            .unwrap_or(actual.len().min(expected.len()));
        let header = MAGIC.len() + 16;
        let record = first.saturating_sub(header) / RECORD_BYTES;
        panic!(
            "Fired log drift for {deck}: first differing byte {first} (event index {record}); \
             actual {} bytes, golden {} bytes",
            actual.len(),
            expected.len()
        );
    }
}

#[test]
fn kaolinite_fired_log_is_bitwise_stable_for_20k_steps() {
    assert_or_update(
        "petra/examples/kaolinite.toml",
        "petra/crates/petra-deck/tests/golden/kaolinite-seed42-20000.fired.bin",
        None,
    );
}

#[test]
fn kossel_fired_log_is_bitwise_stable_for_20k_steps() {
    // The tutorial's 8^3 undersaturated crystal fully dissolves before 20k
    // events. Scale only its lattice extent so the exact same deck physics
    // sustains the RFC-mandated long trajectory.
    assert_or_update(
        "petra/examples/kossel.toml",
        "petra/crates/petra-deck/tests/golden/kossel-20cubed-seed42-20000.fired.bin",
        Some([20, 20, 20]),
    );
}
