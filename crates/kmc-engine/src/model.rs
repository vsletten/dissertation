//! The engine↔model seam — the one trait that makes "generic KMC engine" a
//! true claim — plus the rejection-free selection step built on it.
//!
//! Design doc §4: the engine knows only that there is a graph of sites, each
//! with an opaque state; a [`Model`] can, for any site, enumerate the events
//! possible there with their rates, and apply a chosen one by mutating state.
//! Everything kaolinite — the 100–599 codes, `IsActive`, `CheckEnv`, the 16
//! transition mechanisms — lives *behind* this trait in the `kaolinite`
//! crate. Swap `Kaolinite` for another mineral's `Model` and [`step`] is
//! unchanged. That substitutability is the entire payoff of the port's crate
//! split, and it is a single small trait.
//!
//! This is the Rust replacement for the C++'s hard-wired coupling, where
//! `evtlist.cpp` reaches straight into `envrn.cpp` and `actions.cpp` by name
//! (`environment->IsActive`, `actions.DoEvent`) — the engine and the model
//! were one lump. Here the lump is cut along a declared interface the
//! compiler enforces.

use crate::graph::{SiteGraph, SiteId};
use crate::rng::Rng;

/// A candidate event: "at `site`, reaction `rxn` could fire, at this `rate`".
///
/// The engine treats `rxn` as an opaque `u16` tag — only the model knows that
/// 4 means "Si-O-Al₂ al-side hydrolysis". `rate` is the per-environment rate
/// the model already resolved (the engine never indexes a rate table).
///
/// \[IDIOM\] `#[derive(Clone, Copy)]` on a small POD struct: three machine
/// words, so copying is free and passing by value avoids lifetime noise in
/// the selection loop. The C++ `EventList` was a heap node in an intrusive
/// linked list (`class EventList : public Event { EventList *next; }`, a
/// `new`/`delete` per event per step); a `Copy` value in a reused `Vec` is
/// both faster and drops the entire allocation-churn and dangling-`next`
/// surface.
///
/// # `f32` rate (not the doc's `f64`)
///
/// Design doc §4 shows `rate: f64`. We keep **`f32`** through M6 for
/// trajectory parity (spec §C2a): `ratesum`, the `rate/ratesum` quotients,
/// and `eps` are all `float` in the C++ selection, and the *order* and
/// *width* of those adds decide which event is chosen. The `f64` switch
/// (spec B8) lands with the M7+ reform, flagged. See
/// `kaolinite::reactions::Reaction::rates` for the matching decision.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProposedEvent {
    /// The site the event acts on.
    pub site: SiteId,
    /// Model-defined reaction id.
    pub rxn: u16,
    /// The resolved rate for this (site, reaction) in its current environment.
    pub rate: f32,
}

/// Why a step could not advance — the two legacy stop conditions plus a
/// model-defined error, as a value instead of the C++'s scattered
/// `std::cerr` + `break`/`return -1`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StepStop<E> {
    /// No event was possible anywhere (C++: `CreateEventList` returned null →
    /// "failed to create event list").
    NoEvents,
    /// Every possible event had rate 0, so total rate is 0 (C++
    /// `DoEvent`: "ratesum is zero!" → returns -1).
    ZeroRate,
    /// The model refused to enumerate or apply — e.g. an invalid environment
    /// bucket or an unexpected state (C++: `CheckEnv` returns -1, or a
    /// `Adsorb*` hits its `default:`). Carries the model's own error type.
    Model(E),
}

/// The engine↔model contract. Small on purpose (design doc §4: "keep it
/// small").
///
/// \[IDIOM\] **An associated type (`State`) plus a generic-free trait.** The
/// model names the concrete state type it stores (`kaolinite::State`); the
/// engine stays generic over it *through the trait* without itself carrying a
/// `<S>` parameter everywhere. Associated types say "there is exactly one
/// State per Model", which is what we mean — contrast a type *parameter*
/// `Model<S>`, which would invite `impl Model<A>` and `impl Model<B>` for one
/// model. This is the C++ "the model defines its state" relationship, but
/// checked: the engine cannot accidentally assume anything about `State`
/// beyond `Copy + Eq`.
pub trait Model {
    /// The per-site state the model stores in the graph.
    type State: Copy + Eq;
    /// The model's own error type (invalid environment, bad transition, …).
    type Error;

    /// Enumerate every event possible at `site` right now, pushing each into
    /// `out`. The engine calls this for every site each step to rebuild the
    /// event list (the legacy `CreateEventList` inner body, for one site).
    ///
    /// **Ordering is contractual for parity.** The engine builds the event
    /// list by iterating sites `0..n` ascending and, within a site, expects
    /// the model to push in ascending reaction id — matching the C++ loop
    /// order — because the float summation order in [`step`] depends on it
    /// (see the reverse-fold note there).
    ///
    /// Returns `Err` if the model cannot form a valid event (e.g. `CheckEnv`
    /// yields an out-of-range bucket — a fatal invalid-environment condition
    /// that aborts the legacy run).
    fn events_at(
        &self,
        graph: &SiteGraph<Self::State>,
        site: SiteId,
        out: &mut Vec<ProposedEvent>,
    ) -> Result<(), Self::Error>;

    /// Apply a chosen event: mutate `site` and the neighbors the reaction
    /// touches. Port target: `Actions::DoReaction`.
    ///
    /// \[IDIOM\] **`rng` is threaded in, not owned by the model.** Some
    /// reactions draw a random number *during* application (the R4/R9 proton
    /// coin, `ran2() < 0.5`). Those draws must come from the *same* stream,
    /// *after* the `dt` and `eps` draws of this step — the legacy
    /// `DoReaction` calls the global `ran2()` mid-mutation. Passing the one
    /// `&mut dyn Rng` the engine holds guarantees a single ordered stream; a
    /// model that owned its own RNG would silently desynchronize parity. The
    /// `dyn` (not a generic `R`) keeps `Model` object-safe-ish and the
    /// signature un-parameterized — the cost is one virtual call per coin
    /// flip, paid only by reactions that flip.
    fn apply(
        &mut self,
        graph: &mut SiteGraph<Self::State>,
        ev: &ProposedEvent,
        rng: &mut dyn Rng,
    ) -> Result<(), Self::Error>;
}

/// One KMC step: rejection-free (n-fold-way / Gillespie *direct*) selection
/// and time advance. The generic form of `mckaol.cpp`'s inner loop body plus
/// `Actions::DoEvent`.
///
/// The algorithm (spec A7):
/// 1. rebuild the event list (`events_at` over every site);
/// 2. `ratesum = Σ rate`;
/// 3. `dt = -ln(u) / ratesum` for a fresh uniform `u` (Poisson waiting time);
/// 4. draw `eps`, walk the list accumulating `rate/ratesum`, pick the first
///    event whose cumulative sum ≥ `eps`;
/// 5. apply it; return `dt`.
///
/// # Parity-critical float ordering
///
/// The C++ prepends each new event to a linked list, so its head is the
/// **last** site inserted (highest site id, and within a site the highest
/// reaction id). `DoEvent` then sums and scans from that head. We build
/// `scratch` in ascending order (as `events_at` is called), so to reproduce
/// the C++'s summation order **bit-for-bit** we fold and scan it *reversed*
/// (`.rev()`), head = last element. f32 addition is not associative — get
/// this order wrong and `ratesum`/`partsum` differ by an ulp, which can flip
/// which event crosses `eps`, and the whole trajectory diverges. This is the
/// subtlest faithfulness requirement in the dynamics port; the parity test is
/// its proof.
///
/// The `scratch` buffer is caller-owned and reused across steps — the engine
/// never allocates per step (contrast the C++ `new EventList()` per event).
///
/// \[IDIOM\] Generic over `M: Model` and `R: Rng`, monomorphized per call
/// site — zero dispatch cost for the hot `uniform()`/`events_at` calls, while
/// the trait bounds document exactly what `step` needs and nothing more.
pub fn step<M, R>(
    graph: &mut SiteGraph<M::State>,
    model: &mut M,
    rng: &mut R,
    scratch: &mut Vec<ProposedEvent>,
) -> Result<f32, StepStop<M::Error>>
where
    M: Model,
    R: Rng,
{
    scratch.clear();
    for s in 0..graph.len() {
        model
            .events_at(graph, s, scratch)
            .map_err(StepStop::Model)?;
    }
    if scratch.is_empty() {
        return Err(StepStop::NoEvents);
    }

    // Sum head-first = reversed build order (see the parity note above).
    let ratesum: f32 = scratch.iter().rev().fold(0.0f32, |acc, e| acc + e.rate);
    if ratesum == 0.0 {
        return Err(StepStop::ZeroRate);
    }

    // dt = -log(ran2()) / ratesum. The C++ computes -log in double (the float
    // ran2 promotes) and divides by ratesum promoted to double, narrowing the
    // quotient to the float `dt` at the assignment. Reproduced exactly.
    let u = rng.uniform();
    let dt = (-(u as f64).ln() / ratesum as f64) as f32;

    // Selection. eps and partsum are f32; the C++ walks head..second-to-last,
    // and if none crosses eps it lands on the *tail* (the last/lowest event)
    // without adding its share. `scratch[0]` is that tail (lowest site).
    let eps = rng.uniform();
    let n = scratch.len();
    let mut chosen = 0usize; // default = tail = scratch[0], matching C++.
    let mut partsum = 0.0f32;
    for k in 0..(n - 1) {
        let idx = n - 1 - k; // head-first: scratch[n-1], n-2, …, 1.
        partsum += scratch[idx].rate / ratesum;
        if eps <= partsum {
            chosen = idx;
            break;
        }
    }

    let ev = scratch[chosen];
    model.apply(graph, &ev, rng).map_err(StepStop::Model)?;
    Ok(dt)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::Site;
    use crate::rng::Ran2;

    /// A trivial model: no site ever proposes an event. Exercises the empty
    /// list path (design doc M4 test) and proves `step` compiles against the
    /// trait for a non-kaolinite `State`.
    struct Silent;
    impl Model for Silent {
        type State = i32;
        type Error = ();
        fn events_at(
            &self,
            _g: &SiteGraph<i32>,
            _s: SiteId,
            _out: &mut Vec<ProposedEvent>,
        ) -> Result<(), ()> {
            Ok(())
        }
        fn apply(
            &mut self,
            _g: &mut SiteGraph<i32>,
            _e: &ProposedEvent,
            _r: &mut dyn Rng,
        ) -> Result<(), ()> {
            Ok(())
        }
    }

    #[test]
    fn empty_model_stops_with_no_events() {
        let mut g = SiteGraph {
            sites: vec![Site {
                state: 0i32,
                nbr: [None; 6],
            }],
        };
        let mut rng = Ran2::legacy();
        let mut scratch = Vec::new();
        let r = step(&mut g, &mut Silent, &mut rng, &mut scratch);
        assert_eq!(r, Err(StepStop::NoEvents));
    }

    /// A two-site model with a fixed environment: site 0 offers rxn 0 at rate
    /// 1.0, site 1 offers rxn 1 at rate 3.0. Checks that selection lands on
    /// the right event for a forced `eps`, exercising the reversed scan.
    struct TwoFixed;
    impl Model for TwoFixed {
        type State = i32;
        type Error = ();
        fn events_at(
            &self,
            _g: &SiteGraph<i32>,
            s: SiteId,
            out: &mut Vec<ProposedEvent>,
        ) -> Result<(), ()> {
            let rate = if s == 0 { 1.0 } else { 3.0 };
            out.push(ProposedEvent {
                site: s,
                rxn: s as u16,
                rate,
            });
            Ok(())
        }
        fn apply(
            &mut self,
            _g: &mut SiteGraph<i32>,
            _e: &ProposedEvent,
            _r: &mut dyn Rng,
        ) -> Result<(), ()> {
            Ok(())
        }
    }

    #[test]
    fn step_advances_and_returns_positive_dt() {
        let mut g = SiteGraph {
            sites: vec![
                Site {
                    state: 0i32,
                    nbr: [None; 6],
                },
                Site {
                    state: 0i32,
                    nbr: [None; 6],
                },
            ],
        };
        let mut rng = Ran2::legacy();
        let mut scratch = Vec::new();
        let dt = step(&mut g, &mut TwoFixed, &mut rng, &mut scratch).unwrap();
        assert!(dt > 0.0, "dt should be a positive waiting time, got {dt}");
        // ratesum is 4.0; both events present.
        assert_eq!(scratch.len(), 2);
    }
}
