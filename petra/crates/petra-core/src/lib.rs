//! petra-core — the chemistry-free KMC runtime.
//!
//! This crate knows nothing about minerals: no species names, no state
//! names, no reaction names. It sees a site graph, dense `StateId`s, and
//! compiled [`reaction::Reaction`] tables, all produced by `petra-deck`
//! from a human-readable input deck. Design doc: `petra/docs/DESIGN.md`.

pub mod crystal;
pub mod engine;
pub mod lattice;
pub mod rate;
pub mod reaction;
pub mod state;

pub use engine::{Engine, Fired, Stop};
pub use lattice::{Boundary, Lattice, SiteId};
pub use state::{StateId, StateSet};
