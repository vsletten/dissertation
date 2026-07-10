//! Reaction tables: the 24 live reactions and their per-environment rates.
//! Port of the data structures in `rxnlist.hpp` (the *mechanisms* — what a
//! reaction does to the lattice — arrive at M5 per the ladder).
//!
//! Index layout, inherited from the C++ `#define`s and load order:
//!
//! * 0–15  — hydrolysis, 8 forward/reverse pairs (even = forward)
//! * 16–19 — adsorption (16 Al→Al-site, 17/18 the never-active "wrong
//!   cation" cross variants, 19 Si→Si-site)
//! * 20–23 — desorption (20 Al, 21/23 cross variants, 22 Si)
//! * 24–27 — diffusion: **not ported.** Dead code (spec B6: `IsActive`
//!   never allows it) built with a genuine out-of-bounds read (spec B7:
//!   the C++ copies 5 rates out of a 4-rate table because `numRates` is
//!   stale). The task decree is "don't carry the dead weight, don't
//!   reintroduce the over-read" — so [`ReactionSet::reactions`] holds 24
//!   entries and any future diffusion port must re-derive its tables on
//!   purpose rather than inherit a buffer over-read by accident.
//!
//! Note on the shipped `data.rxn` (spec B4): every bucket within each
//! reaction carries identical (k, ΔE), so the per-environment machinery is
//! *inert* in the sample data. Faithful ≠ validated: exercising the
//! environment buckets needs an authored non-flat `data.rxn` (M8).

use crate::state::State;

/// Reaction index where the "400-class" hydrolysis block begins (C++
/// `N300`): reactions 0..2 are the Si-O-Si (300) pair R0/R1, 2..14 the
/// 400-class R2..R13. Used by the event-list builder to pick the reaction
/// window for a given site class (spec A7 / `evtlist.cpp`).
pub const N_300: usize = 2;
/// Reaction index where the "500-class" hydrolysis block begins (C++ `N400`):
/// reactions 14..16 are the Al-OH-Al (500) pair R14/R15.
pub const N_400: usize = 14;
/// Number of hydrolysis reactions (C++ `NHYD`).
pub const N_HYD: usize = 16;
/// End of the adsorption block (C++ `NADS`).
pub const N_ADS: usize = 20;
/// End of the desorption block (C++ `NDES`) — and, with diffusion unported,
/// the length of [`ReactionSet::reactions`].
pub const N_DES: usize = 24;

/// Gas constant in kcal/mol/K (C++ `#define R 1.987e-3`). An `f64` because
/// the C++ macro is a double literal — `rt = R * t` multiplies in double
/// before narrowing to float, and the reader must reproduce that exact
/// rounding path (see `kmc-io::rxn`).
pub const R_KCAL: f64 = 1.987e-3;

/// One reaction: what it consumes and its rate per environment bucket
/// (C++ `class Reaction`).
#[derive(Debug, Clone, PartialEq)]
pub struct Reaction {
    /// The state code this reaction consumes.
    pub reactant: State,
    /// Product state code — set only for adsorption reactions, where the
    /// event executor needs it for bookkeeping (the C++ leaves this member
    /// *uninitialized* for hydrolysis/desorption and never reads it there;
    /// we store 0, which is the honest Rust spelling of "never meaningful").
    pub info: i32,
    /// One rate per nth-neighbor environment bucket; `CheckEnv` (M5)
    /// produces the index. Length varies by reaction (15/40/20 for
    /// hydrolysis families, 1 for adsorption, 4/5 for desorption).
    ///
    /// `f32` — NOT the design doc's recommended `f64` (§4). The doc's
    /// upgrade is right for the *reformed* model, but rates feed `ratesum`
    /// and `dt`, and the dynamics-parity milestone (M6) chases the C++
    /// trajectory bit-for-bit first (spec §C2a: f32 once, then switch).
    /// Faithful now, improved later, never silently.
    pub rates: Vec<f32>,
}

/// Everything `data.rxn` provides: conditions plus the reaction tables
/// (C++ `class ReactionList`, minus the dead diffusion tail).
#[derive(Debug, Clone, PartialEq)]
pub struct ReactionSet {
    /// Temperature in Kelvin (8000.0 in the golden file — yes, really;
    /// it's a relative-rate knob, not a physical claim).
    pub temperature: f32,
    /// Δμ for Si — solution supersaturation driving Si adsorption
    /// (C++ `GetSiPotential`, consumed by `populate_solid`).
    pub dm_si: f32,
    /// Δμ for Al (C++ `GetAlPotential`).
    pub dm_al: f32,
    /// Reactions 0–23. Diffusion (24–27) intentionally absent — see the
    /// module docs.
    pub reactions: Vec<Reaction>,
}
