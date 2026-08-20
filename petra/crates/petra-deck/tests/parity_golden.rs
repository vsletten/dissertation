//! B2's migration-grade parity gates.
//!
//! Checked-in byte logs were captured from the pre-refactor engine on main
//! (4522ecc). Tests are deliberately read-only: oracle regeneration lives in
//! the explicit `generate-parity-goldens` maintenance binary.

use std::path::{Path, PathBuf};

const SEED: u64 = 42;
const STEPS: u64 = 20_000;
const COMPLETE_KOSSEL_STEPS: u64 = 1_512;
const MAGIC: &[u8] = b"PETRA-FIRED-V1\0";
const RECORD_BYTES: usize = 8 + 8 + 8 + 2;

fn repo_path(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

fn fired_log(deck_path: &Path, dims_override: Option<[usize; 3]>, limit: Option<u64>) -> Vec<u8> {
    let text = std::fs::read_to_string(deck_path).expect("read parity deck");
    let mut parsed: petra_deck::DeckFile = toml::from_str(&text).expect("parity deck parses");
    if let Some(dims) = dims_override {
        parsed.lattice.dims = dims;
    }
    let deck = petra_deck::compile(&parsed).expect("parity deck compiles");
    let mut engine = deck.build_engine(Some(SEED)).expect("engine builds");
    let mut records = Vec::new();
    while limit.is_none_or(|limit| engine.step_count < limit) {
        let fired = match engine.step() {
            Ok(fired) => fired,
            Err(petra_core::Stop::NoEvents | petra_core::Stop::ZeroRate) if limit.is_none() => {
                break
            }
            Err(error) => panic!("deck stopped at step {}: {error}", engine.step_count + 1),
        };
        records.extend_from_slice(&fired.step.to_le_bytes());
        records.extend_from_slice(&fired.time.to_bits().to_le_bytes());
        records.extend_from_slice(&(fired.site as u64).to_le_bytes());
        records.extend_from_slice(&fired.reaction.to_le_bytes());
    }
    if let Some(limit) = limit {
        assert_eq!(
            engine.step_count, limit,
            "deck stopped before requested limit"
        );
    }
    let mut bytes = Vec::with_capacity(MAGIC.len() + 16 + records.len());
    bytes.extend_from_slice(MAGIC);
    bytes.extend_from_slice(&SEED.to_le_bytes());
    bytes.extend_from_slice(&engine.step_count.to_le_bytes());
    bytes.extend_from_slice(&records);
    bytes
}

fn assert_golden(deck: &str, golden: &str, dims_override: Option<[usize; 3]>, limit: Option<u64>) {
    let actual = fired_log(&repo_path(deck), dims_override, limit);
    let golden_path = repo_path(golden);
    let expected = std::fs::read(&golden_path).unwrap_or_else(|error| {
        panic!(
            "cannot read {}: {error}; use only the explicit generator on the pinned baseline",
            golden_path.display()
        )
    });
    if actual != expected {
        let first = actual
            .iter()
            .zip(&expected)
            .position(|(actual, expected)| actual != expected)
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
    assert_golden(
        "petra/examples/kaolinite.toml",
        "petra/crates/petra-deck/tests/golden/kaolinite-seed42-20000.fired.bin",
        None,
        Some(STEPS),
    );
}

#[test]
fn enlarged_kossel_fired_log_is_bitwise_stable_for_20k_steps() {
    assert_golden(
        "petra/examples/kossel.toml",
        "petra/crates/petra-deck/tests/golden/kossel-20cubed-seed42-20000.fired.bin",
        Some([20, 20, 20]),
        Some(STEPS),
    );
}

#[test]
fn shipped_kossel_complete_natural_trajectory_is_bitwise_stable() {
    let actual = fired_log(&repo_path("petra/examples/kossel.toml"), None, None);
    let encoded_steps = u64::from_le_bytes(
        actual[MAGIC.len() + 8..MAGIC.len() + 16]
            .try_into()
            .expect("step header"),
    );
    assert_eq!(encoded_steps, COMPLETE_KOSSEL_STEPS);
    assert_golden(
        "petra/examples/kossel.toml",
        "petra/crates/petra-deck/tests/golden/kossel-8cubed-seed42-complete.fired.bin",
        None,
        None,
    );
}
