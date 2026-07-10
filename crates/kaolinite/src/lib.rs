//! The kaolinite model — the *specific* chemistry the generic KMC engine
//! knows nothing about.
//!
//! Kaolinite is a 1:1 Al-Si clay: a sheet of Si-O tetrahedra bonded to a
//! sheet of Al-O(H) octahedra. The model simulates its surface dissolving
//! and regrowing in water, one reaction at a time, on a fixed lattice of
//! sites whose only mutable property is a 3-digit integer *state code*
//! (see [`state`]). The scientific point of the model is that a site's
//! reaction rates depend on the chemical state of its *nth-nearest
//! neighbors* — real surfaces are not mean-field.
//!
//! This crate owns everything kaolinite-specific:
//!
//! | module | C++ counterpart | job |
//! |---|---|---|
//! | [`state`] | `common.hpp` + `lattice.hpp` macros | the state encoding |
//! | [`cell`] | `ucell.hpp` | the unit-cell motif that tiles into the lattice |
//! | [`build`] | `lattice.cpp` | tile the motif into a site graph; structural setup |
//!
//! (The ladder adds `reactions`/`environment`/the engine `Model` impl at M4+.)
//!
//! Behavior contract: `mission-control/projects/kmc/02-model-spec.md`.
//! Warts preserved on purpose are marked `WART (spec BN)` in place.

#![warn(missing_docs)]

pub mod build;
pub mod cell;
pub mod state;

pub use build::{create_lattice, LatticeParams, Structure};
pub use cell::{CellSite, NeighborTemplate, UnitCell};
pub use state::State;
