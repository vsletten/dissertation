//! petra-deck — the human-facing side of Petra: the TOML deck schema and
//! the compiler that lowers a deck to `petra-core`'s dense runtime tables.
//! Design doc: `petra/docs/DESIGN.md` §3–§4.

pub mod compile;
pub mod schema;

use std::path::Path;

pub use compile::{compile, replica_seed, CompileError, CompiledDeck, ExecutionStrategy};
pub use schema::{DeckFile, SeedPolicy, StructureKind};

#[derive(Debug, thiserror::Error)]
pub enum LoadError {
    #[error("cannot read deck file: {0}")]
    Io(#[from] std::io::Error),
    #[error("deck is not valid TOML/schema: {0}")]
    Parse(#[from] toml::de::Error),
    #[error(transparent)]
    Compile(#[from] CompileError),
}

/// Read, parse, and compile a deck file in one call.
pub fn load(path: impl AsRef<Path>) -> Result<CompiledDeck, LoadError> {
    let text = std::fs::read_to_string(path)?;
    let deck: DeckFile = toml::from_str(&text)?;
    Ok(compile(&deck)?)
}
