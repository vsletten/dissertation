//! Readers and writers for the legacy KMC text formats.
//!
//! The C++ model reads four fixed-name, whitespace-delimited text files from
//! the current working directory (`data.sim`, `data.cell`, `data.lattice`,
//! `data.rxn`) and writes Cerius2 MSI structure files plus a few CSV-ish
//! outputs. This crate ports those formats **behavior-for-behavior**,
//! including the input-reading warts cataloged in the faithful spec
//! (`projects/kmc/02-model-spec.md`, Part B). Where a wart is preserved, the
//! code carries a `WART (spec BN)` comment at the exact spot.
//!
//! Layering: this crate depends on `kaolinite` (the domain types it parses
//! into / serializes from) but never on the simulation loop. Nothing here
//! mutates model state; it is the IO shell around the hexagon, same role the
//! `db`/`extractor` adapters play in Victor's other codebases.
//!
//! # Module map
//!
//! | module | C++ counterpart | job |
//! |---|---|---|
//! | [`scan`] | `futil.cpp` (`EatComment` + `operator>>`) | tokenizer for the input dialect |
//! | [`sim`] | `sim.cpp` | `data.sim` → [`sim::SimParams`] |
//! | [`cell`] | `ucell.cpp` | `data.cell` → [`kaolinite::UnitCell`] |
//! | [`lattice`] | `lattice.cpp` (top) | `data.lattice` → [`kaolinite::LatticeParams`] |
//! | [`rxn`] | `rxnlist.cpp` | `data.rxn` → [`kaolinite::reactions::ReactionSet`] |
//! | [`view`] | (none — the C++ passes `Lattice*`) | the writers' borrowed view + shared render geometry |
//! | [`msi`] | `output.cpp` (`writeMSI`) | lattice → Cerius2 MSI (bug-compatible) |
//! | [`xyz`] | `output.cpp` (`writeXYZ`) | lattice → legacy comma-XYZ dialect |
//! | [`dat`] | `output.cpp` (`writeData`) | population row (one per snapshot file, spec B1) |
//! | [`surf`] | `output.cpp` (`writeSurf`) | surface x,y projections (`surfSi.out`/`surfAl.out`) |
//! | [`fmt`] | `ostream <<` defaults | C++ `%g`-equivalent float formatting |
//! | [`error`] | `myerr.cpp` (`Myerr::die`) | error type all readers return |
//!
//! (A later milestone adds PGIF.)

// [IDIOM] `#![warn(missing_docs)]` — a crate-wide lint gate. The compiler
// itself becomes the documentation-coverage reviewer: any public item
// without a doc comment is flagged at build time. There is no C++ analogue
// short of bolting Doxygen coverage checks onto CI. For agent-written code
// this class of lint is a cheap, machine-enforced review pass — turn it on
// early and the codebase can never silently grow an undocumented API.
#![warn(missing_docs)]

pub mod cell;
pub mod dat;
pub mod error;
pub mod fmt;
pub mod lattice;
pub mod msi;
pub mod rxn;
pub mod scan;
pub mod sim;
pub mod surf;
pub mod view;
pub mod xyz;

// [IDIOM] `pub use` re-exports. C++ headers leak every name they include;
// Rust modules hide everything by default and you *choose* the public face.
// Re-exporting the two most-used names at the crate root gives callers
// `kmc_io::SimParams` without knowing the module layout — the module tree is
// an implementation detail, the re-exports are the API.
pub use cell::read_cell;
pub use dat::write_data;
pub use error::ReadError;
pub use lattice::read_lattice;
pub use msi::write_msi;
pub use rxn::read_rxn;
pub use sim::{SimParams, read_sim};
pub use surf::write_surf;
pub use view::LatticeView;
pub use xyz::write_xyz;
