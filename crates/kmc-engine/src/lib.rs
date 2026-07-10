//! Generic, model-agnostic kinetic Monte Carlo core.
//!
//! What "generic" means here (design doc §1): the engine knows there is a
//! graph of sites, each carrying an *opaque* state; a model can enumerate
//! possible events with rates and apply a chosen one. It knows nothing about
//! kaolinite. At M2 the crate holds only the graph; the event loop, `Model`
//! trait, and the faithful `ran2` PRNG arrive at M4 per the milestone ladder
//! (design doc §6).
//!
//! # The concurrency road not taken (yet)
//!
//! Everything here is single-threaded on purpose. KMC is intrinsically
//! serial — one global clock, one event at a time, each event's
//! probabilities conditioned on the last (parallelism memo §2). If
//! parallelism ever pays (memo §4: synchronous sublattice decomposition,
//! only for lattices far larger than any committed sample), *this* is the
//! type it would partition: strips of the flat site array with 2-hop halos,
//! one worker per strip, `rayon`-or-`std::thread` over disjoint
//! `&mut [Site<S>]` chunks. Rust's ownership model is what makes that future
//! tractable: `split_at_mut` hands each worker a provably disjoint slice,
//! and the compiler rejects any halo bookkeeping that aliases — the data
//! races a C++ decomposition would debug at runtime are compile errors
//! here. Nothing in the serial design forecloses it; that is the seam
//! discipline the memo asks for, and it costs nothing today.

#![warn(missing_docs)]

pub mod graph;

pub use graph::{Site, SiteGraph, SiteId};
