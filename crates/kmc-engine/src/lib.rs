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
//!
//! Where `rayon` would first enter (memo 04 tie-in): the event-list rebuild
//! in [`step`] — `for s in 0..graph.len() { model.events_at(graph, s, ..) }`
//! — is a read-only fan-out over an immutable `&graph`, the textbook
//! `par_iter().flat_map(..).collect()` shape. It is the cheapest parallelism
//! win and the one to reach for *first* if a large lattice ever makes the
//! per-step rebuild the bottleneck, precisely because `events_at` takes
//! `&self`/`&graph`: the borrow checker already proves the workers share no
//! mutable state. **But** it must stay bit-parity-gated — a parallel collect
//! reorders the event list, and §the selection-order note in [`model::step`]
//! shows why reordering breaks the f32 summation the legacy oracle pins. So
//! parallel rebuild is a *reform-era* (post-parity) move: correct-by-default
//! reorders freely, `--legacy` keeps the serial order. The `apply` step is
//! intrinsically serial (it mutates and draws from the shared RNG); memo
//! 04's sublattice decomposition is the only route to parallelizing *that*,
//! and it is deliberately out of scope until the serial port is trusted.

#![warn(missing_docs)]

pub mod graph;
pub mod model;
pub mod rng;

pub use graph::{Site, SiteGraph, SiteId};
pub use model::{Model, ProposedEvent, StepStop, step};
pub use rng::{Ran2, Rng};
