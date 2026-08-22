//! Explicit maintenance tool for Petra's pinned parity oracles.
//!
//! Run only from the reviewed pre-refactor baseline (or a commit already
//! proven byte-identical to it):
//! `cargo run -p petra-deck --bin generate-parity-goldens -- --regenerate-pinned-oracles`

use std::path::{Path, PathBuf};

const SEED: u64 = 42;
const MAGIC: &[u8] = b"PETRA-FIRED-V1\0";

fn repo_path(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

fn fired_records(
    deck_path: &Path,
    dims_override: Option<[usize; 3]>,
    limit: Option<u64>,
) -> Vec<u8> {
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
    let steps = engine.step_count;
    if let Some(limit) = limit {
        assert_eq!(steps, limit, "long-run deck stopped before requested limit");
    }
    let mut bytes = Vec::with_capacity(MAGIC.len() + 16 + records.len());
    bytes.extend_from_slice(MAGIC);
    bytes.extend_from_slice(&SEED.to_le_bytes());
    bytes.extend_from_slice(&steps.to_le_bytes());
    bytes.extend_from_slice(&records);
    bytes
}

fn write_golden(deck: &str, golden: &str, dims_override: Option<[usize; 3]>, limit: Option<u64>) {
    let bytes = fired_records(&repo_path(deck), dims_override, limit);
    let path = repo_path(golden);
    std::fs::create_dir_all(path.parent().expect("golden has parent"))
        .expect("create golden directory");
    std::fs::write(&path, &bytes).expect("write parity golden");
    println!("wrote {} ({} bytes)", path.display(), bytes.len());
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    assert_eq!(
        args,
        ["--regenerate-pinned-oracles"],
        "refusing to overwrite parity oracles; pass the explicit maintenance flag"
    );
    write_golden(
        "petra/examples/kaolinite.toml",
        "petra/crates/petra-deck/tests/golden/kaolinite-seed42-20000.fired.bin",
        None,
        Some(20_000),
    );
    write_golden(
        "petra/examples/kossel.toml",
        "petra/crates/petra-deck/tests/golden/kossel-20cubed-seed42-20000.fired.bin",
        Some([20, 20, 20]),
        Some(20_000),
    );
    write_golden(
        "petra/examples/kossel.toml",
        "petra/crates/petra-deck/tests/golden/kossel-8cubed-seed42-complete.fired.bin",
        None,
        None,
    );
}
