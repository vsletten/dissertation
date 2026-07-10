//! The 16 hydrolysis mechanisms plus adsorption/desorption — *what a chosen
//! reaction does to the lattice*. Port of `actions.cpp`
//! (`DoReaction0`…`DoReaction15`, `AdsorbAl`/`AdsorbSi`,
//! `DesorbAl`/`DesorbSi`).
//!
//! Each function rewrites the reacting oxygen's `state` and adjusts the
//! coordination of the participating cation neighbors, exactly as the C++.
//! The bookkeeping that makes the reverse reactions exact inverses lives in
//! two parallel arrays the model owns:
//!
//! * `lostal[o]` — for a 400 oxygen that lost one of its two Al, *which* Al
//!   (C++ `LatticeSite::lostal`). Set by R4/R9 and `AdsorbAl`, read by
//!   R5/R8 and `DesorbAl`.
//! * `pair` — the double-bridge partner (built once by `find_pairs`, never
//!   mutated here) — used only by the environment checks, not these
//!   mutations, so it is absent from this module.
//!
//! # The proton coin (R4, R9) and RNG order
//!
//! R4 and R9 each draw **one** `ran2()` mid-mutation to decide which of two
//! symmetric Al takes the proton, recording the choice in `lostal` so the
//! reverse can undo it. That draw comes from the *same* stream as the step's
//! `dt` and `eps` draws, *after* them — so these functions take the shared
//! `&mut dyn Rng`. Getting this draw's position in the stream wrong
//! desynchronizes every subsequent step (this is the single most parity-
//! fragile spot in the dynamics; the trajectory oracle pins it).
//!
//! # Neighbor access
//!
//! The C++ reads `sites[site].nbr[i]` as a raw index and dereferences it
//! unconditionally — a valid reaction always has the neighbors it needs. We
//! mirror that with [`nb`], which `expect`s the neighbor present: the
//! invariant is "the event-list builder only proposed this reaction because
//! its reactant state matched, which by construction means the required
//! neighbors exist". A panic here would mean the *builder* is wrong, not the
//! input — exactly the class of bug an `expect`-with-invariant should surface
//! loudly rather than paper over.

use kmc_engine::{Rng, SiteGraph, SiteId};

use crate::model_impl::KaolError;
use crate::state::State;

/// The `nbr[i]` of `site`, asserted present (see the module note on why
/// `expect` is correct here).
#[inline]
fn nb(graph: &SiteGraph<State>, site: SiteId, i: usize) -> SiteId {
    graph.sites[site].nbr[i].expect("reaction neighbor must exist for a proposed event")
}

/// Read a site's raw state code.
#[inline]
fn st(graph: &SiteGraph<State>, site: SiteId) -> i32 {
    graph.sites[site].state.0
}

/// Set a site's raw state code.
#[inline]
fn set(graph: &mut SiteGraph<State>, site: SiteId, s: i32) {
    graph.sites[site].state.0 = s;
}

/// Add `d` to a site's raw state code (the C++ `state++` / `state--`).
#[inline]
fn bump(graph: &mut SiteGraph<State>, site: SiteId, d: i32) {
    graph.sites[site].state.0 += d;
}

/// Dispatch a chosen reaction to its mechanism. Port of
/// `Actions::DoReaction`'s `switch`. Reactions 16/17 → adsorb Al, 18/19 →
/// adsorb Si, 20/23 → desorb Al, 21/22 → desorb Si (the cross variants
/// 17/18/21/23 are never proposed — `is_active` gates them off — but the
/// mapping is ported so the `switch` reads identically).
pub fn do_reaction(
    graph: &mut SiteGraph<State>,
    lostal: &mut [Option<SiteId>],
    site: SiteId,
    rxn: u16,
    rng: &mut dyn Rng,
) -> Result<(), KaolError> {
    match rxn {
        0 => do_reaction0(graph, site),
        1 => do_reaction1(graph, site),
        2 => do_reaction2(graph, site),
        3 => do_reaction3(graph, site),
        4 => do_reaction4(graph, lostal, site, rng),
        5 => do_reaction5(graph, lostal, site),
        6 => do_reaction6(graph, site),
        7 => do_reaction7(graph, site),
        8 => do_reaction8(graph, lostal, site),
        9 => do_reaction9(graph, lostal, site, rng),
        10 => do_reaction10(graph, site),
        11 => do_reaction11(graph, site),
        12 => do_reaction12(graph, site),
        13 => do_reaction13(graph, site),
        14 => do_reaction14(graph, site),
        15 => do_reaction15(graph, site),
        16 | 17 => adsorb_al(graph, lostal, site),
        18 | 19 => adsorb_si(graph, site),
        20 | 23 => {
            desorb_al(graph, lostal, site);
            Ok(())
        }
        21 | 22 => {
            desorb_si(graph, site);
            Ok(())
        }
        // Diffusion (24–27): dead code (spec B6), never proposed. A faithful
        // "reached the impossible default" signal rather than silent success.
        other => Err(KaolError::UnsupportedReaction(other)),
    }
}

// --- Hydrolysis R0–R15 -----------------------------------------------------

/// R0: Si-O-Si + H₂O → 2 Si-OH. 301 → 302; each Si `nbr[0..2]` gains an OH.
fn do_reaction0(graph: &mut SiteGraph<State>, site: SiteId) -> Result<(), KaolError> {
    set(graph, site, 302);
    for i in 0..2 {
        let si = nb(graph, site, i);
        bump(graph, si, 1);
    }
    Ok(())
}

/// R1: reverse of R0. 302 → 301; each Si loses an OH.
fn do_reaction1(graph: &mut SiteGraph<State>, site: SiteId) -> Result<(), KaolError> {
    set(graph, site, 301);
    for i in 0..2 {
        let si = nb(graph, site, i);
        bump(graph, si, -1);
    }
    Ok(())
}

/// R2: Si-O-Al₂ + H₂O → Si-OH + Al-OH-Al. 401 → 402; the Si (`nbr[2]`) gains.
fn do_reaction2(graph: &mut SiteGraph<State>, site: SiteId) -> Result<(), KaolError> {
    set(graph, site, 402);
    let si = nb(graph, site, 2);
    bump(graph, si, 1);
    Ok(())
}

/// R3: reverse of R2. 402 → 401; the Si loses.
fn do_reaction3(graph: &mut SiteGraph<State>, site: SiteId) -> Result<(), KaolError> {
    set(graph, site, 401);
    let si = nb(graph, site, 2);
    bump(graph, si, -1);
    Ok(())
}

/// R4: Si-O-Al₂ + H₂O → Si-OH-Al + Al-OH. 401 → 410; **the proton coin**
/// picks which Al (`nbr[0]` or `nbr[1]`) gains, recorded in `lostal`.
fn do_reaction4(
    graph: &mut SiteGraph<State>,
    lostal: &mut [Option<SiteId>],
    site: SiteId,
    rng: &mut dyn Rng,
) -> Result<(), KaolError> {
    set(graph, site, 410);
    let r = rng.uniform();
    let nbr = if r < 0.5 {
        nb(graph, site, 0)
    } else {
        nb(graph, site, 1)
    };
    bump(graph, nbr, 1);
    lostal[site] = Some(nbr);
    Ok(())
}

/// R5: reverse of R4. 410 → 401; the `lostal` Al loses, `lostal` cleared.
fn do_reaction5(
    graph: &mut SiteGraph<State>,
    lostal: &mut [Option<SiteId>],
    site: SiteId,
) -> Result<(), KaolError> {
    set(graph, site, 401);
    let nbr = lostal[site].expect("R5 requires the R4 lostal record");
    bump(graph, nbr, -1);
    lostal[site] = None;
    Ok(())
}

/// R6: second-stage Al-O-Al hydrolysis. 402 → 403; both Al (`nbr[0..2]`) gain.
fn do_reaction6(graph: &mut SiteGraph<State>, site: SiteId) -> Result<(), KaolError> {
    set(graph, site, 403);
    for i in 0..2 {
        let nbr = nb(graph, site, i);
        bump(graph, nbr, 1);
    }
    Ok(())
}

/// R7: reverse of R6. 403 → 402; both Al lose.
fn do_reaction7(graph: &mut SiteGraph<State>, site: SiteId) -> Result<(), KaolError> {
    set(graph, site, 402);
    for i in 0..2 {
        let nbr = nb(graph, site, i);
        bump(graph, nbr, -1);
    }
    Ok(())
}

/// R8: Si-OH-Al + Al-OH → Si-OH + Al-OH + Al-H₂O. 410 → 403; the Si gains,
/// and the *other* Al (the one that isn't `lostal`) gains; `lostal` cleared.
fn do_reaction8(
    graph: &mut SiteGraph<State>,
    lostal: &mut [Option<SiteId>],
    site: SiteId,
) -> Result<(), KaolError> {
    set(graph, site, 403);
    let si = nb(graph, site, 2);
    bump(graph, si, 1);
    // C++: nbr = (lostal == nbr[0]) ? nbr[1] : nbr[0].
    let al0 = nb(graph, site, 0);
    let al1 = nb(graph, site, 1);
    let other = if lostal[site] == Some(al0) { al1 } else { al0 };
    bump(graph, other, 1);
    lostal[site] = None;
    Ok(())
}

/// R9: reverse of R8. 403 → 410; the Si loses, the **coin**-chosen Al loses,
/// the other is recorded as `lostal`.
fn do_reaction9(
    graph: &mut SiteGraph<State>,
    lostal: &mut [Option<SiteId>],
    site: SiteId,
    rng: &mut dyn Rng,
) -> Result<(), KaolError> {
    set(graph, site, 410);
    let r = rng.uniform();
    let si = nb(graph, site, 2);
    bump(graph, si, -1);
    let (nbr, lost) = if r < 0.5 {
        (nb(graph, site, 0), nb(graph, site, 1))
    } else {
        (nb(graph, site, 1), nb(graph, site, 0))
    };
    bump(graph, nbr, -1);
    lostal[site] = Some(lost);
    Ok(())
}

/// R10: Si-OH-Al + H₂O → Si-OH + Al-H₂O. 406 → 407; the Si gains and the
/// *occupied* Al (the one that isn't empty `100`) gains.
fn do_reaction10(graph: &mut SiteGraph<State>, site: SiteId) -> Result<(), KaolError> {
    set(graph, site, 407);
    let si = nb(graph, site, 2);
    bump(graph, si, 1);
    let al1 = nb(graph, site, 0);
    let al2 = nb(graph, site, 1);
    let nbr = if st(graph, al1) == 100 { al2 } else { al1 };
    bump(graph, nbr, 1);
    Ok(())
}

/// R11: reverse of R10. 407 → 406; the Si loses and the occupied Al loses.
fn do_reaction11(graph: &mut SiteGraph<State>, site: SiteId) -> Result<(), KaolError> {
    set(graph, site, 406);
    let si = nb(graph, site, 2);
    bump(graph, si, -1);
    let al1 = nb(graph, site, 0);
    let al2 = nb(graph, site, 1);
    let nbr = if st(graph, al1) == 100 { al2 } else { al1 };
    bump(graph, nbr, -1);
    Ok(())
}

/// R12: Al-OH-Al + H₂O → Al-OH + Al-H₂O (400 side). 404 → 405; both Al gain.
fn do_reaction12(graph: &mut SiteGraph<State>, site: SiteId) -> Result<(), KaolError> {
    set(graph, site, 405);
    for i in 0..2 {
        let nbr = nb(graph, site, i);
        bump(graph, nbr, 1);
    }
    Ok(())
}

/// R13: reverse of R12. 405 → 404; both Al lose. (Note: the shipped
/// `data.rxn` mis-types R13's reactant as 40100 — spec-noted — so this
/// mechanism is never *selected* from the sample data, but the transition is
/// ported faithfully for a corrected input.)
fn do_reaction13(graph: &mut SiteGraph<State>, site: SiteId) -> Result<(), KaolError> {
    set(graph, site, 404);
    for i in 0..2 {
        let nbr = nb(graph, site, i);
        bump(graph, nbr, -1);
    }
    Ok(())
}

/// R14: Al-OH-Al + H₂O → Al-OH + Al-H₂O (500 side). 501 → 502; both Al gain.
fn do_reaction14(graph: &mut SiteGraph<State>, site: SiteId) -> Result<(), KaolError> {
    set(graph, site, 502);
    for i in 0..2 {
        let nbr = nb(graph, site, i);
        bump(graph, nbr, 1);
    }
    Ok(())
}

/// R15: reverse of R14. 502 → 501; both Al lose.
fn do_reaction15(graph: &mut SiteGraph<State>, site: SiteId) -> Result<(), KaolError> {
    set(graph, site, 501);
    for i in 0..2 {
        let nbr = nb(graph, site, i);
        bump(graph, nbr, -1);
    }
    Ok(())
}

// --- Adsorption / desorption ----------------------------------------------

/// `AdsorbAl`: an empty Al site becomes Al(OH,H₂O)₆ (107); each of its 6
/// oxygen neighbors advances one protonation step along a fixed map. The
/// 406→410 case records the losing Al in `lostal`. An oxygen in an
/// unexpected state is a fatal error (C++ prints and returns false).
fn adsorb_al(
    graph: &mut SiteGraph<State>,
    lostal: &mut [Option<SiteId>],
    site: SiteId,
) -> Result<(), KaolError> {
    if graph.sites[site].state.class_code() == 1 {
        set(graph, site, 107);
        for i in 0..6 {
            let nbr = nb(graph, site, i);
            match st(graph, nbr) {
                503 => set(graph, nbr, 502),
                500 => set(graph, nbr, 503),
                407 => set(graph, nbr, 403),
                409 => set(graph, nbr, 405),
                408 => set(graph, nbr, 407),
                400 => set(graph, nbr, 409),
                406 => {
                    set(graph, nbr, 410);
                    lostal[nbr] = Some(site);
                }
                other => return Err(KaolError::InvalidAdsorbState { site: nbr, state: other }),
            }
        }
    } else {
        // Cross variant (wrong cation): 200 + WRONG. Never reached in
        // practice (is_active gates rxn 17 off); ported for completeness.
        set(graph, site, 299);
    }
    Ok(())
}

/// `AdsorbSi`: an empty Si site becomes Si(OH)₄ (205); each of its 4 oxygen
/// neighbors advances one step.
fn adsorb_si(graph: &mut SiteGraph<State>, site: SiteId) -> Result<(), KaolError> {
    if graph.sites[site].state.class_code() == 2 {
        set(graph, site, 205);
        for i in 0..4 {
            let nbr = nb(graph, site, i);
            match st(graph, nbr) {
                300 => set(graph, nbr, 303),
                303 => set(graph, nbr, 302),
                404 => set(graph, nbr, 402),
                405 => set(graph, nbr, 403),
                409 => set(graph, nbr, 407),
                400 => set(graph, nbr, 408),
                other => return Err(KaolError::InvalidAdsorbState { site: nbr, state: other }),
            }
        }
    } else {
        set(graph, site, 199);
    }
    Ok(())
}

/// `DesorbAl`: the exact inverse of [`adsorb_al`], restoring the Al to empty
/// (100) and walking the oxygens back down. Unlisted oxygen states are left
/// unchanged (C++ `default: break;` — no error).
fn desorb_al(graph: &mut SiteGraph<State>, lostal: &mut [Option<SiteId>], site: SiteId) {
    if graph.sites[site].state.class_code() == 1 {
        set(graph, site, 100);
        for i in 0..6 {
            let nbr = nb(graph, site, i);
            match st(graph, nbr) {
                502 => set(graph, nbr, 503),
                503 => set(graph, nbr, 500),
                403 => set(graph, nbr, 407),
                405 => set(graph, nbr, 409),
                407 => set(graph, nbr, 408),
                409 => set(graph, nbr, 400),
                410 => {
                    set(graph, nbr, 406);
                    lostal[nbr] = None;
                }
                _ => {}
            }
        }
    } else {
        set(graph, site, 200);
    }
}

/// `DesorbSi`: the exact inverse of [`adsorb_si`].
fn desorb_si(graph: &mut SiteGraph<State>, site: SiteId) {
    if graph.sites[site].state.class_code() == 2 {
        set(graph, site, 200);
        for i in 0..4 {
            let nbr = nb(graph, site, i);
            match st(graph, nbr) {
                303 => set(graph, nbr, 300),
                302 => set(graph, nbr, 303),
                402 => set(graph, nbr, 404),
                403 => set(graph, nbr, 405),
                407 => set(graph, nbr, 409),
                408 => set(graph, nbr, 400),
                _ => {}
            }
        }
    } else {
        set(graph, site, 100);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use kmc_engine::{Ran2, Site};

    fn graph_of(states_nbrs: Vec<(i32, [Option<SiteId>; 6])>) -> SiteGraph<State> {
        SiteGraph {
            sites: states_nbrs
                .into_iter()
                .map(|(s, nbr)| Site {
                    state: State(s),
                    nbr,
                })
                .collect(),
        }
    }

    #[test]
    fn r0_and_r1_are_inverses() {
        // Si(0)-O(2)-Si(1): oxygen 301, two Si at 201.
        let make = || {
            graph_of(vec![
                (201, [Some(2), None, None, None, None, None]),
                (201, [Some(2), None, None, None, None, None]),
                (301, [Some(0), Some(1), None, None, None, None]),
            ])
        };
        let mut g = make();
        let mut lostal = vec![None; 3];
        let mut rng = Ran2::legacy();
        do_reaction(&mut g, &mut lostal, 2, 0, &mut rng).unwrap();
        assert_eq!(st(&g, 2), 302);
        assert_eq!(st(&g, 0), 202);
        assert_eq!(st(&g, 1), 202);
        do_reaction(&mut g, &mut lostal, 2, 1, &mut rng).unwrap();
        assert_eq!(g, make()); // full round-trip
    }

    #[test]
    fn r4_coin_records_lostal_and_r5_undoes_it() {
        // 401 oxygen bridging Al(0), Al(1); nbr[2] = Si(3).
        let make = || {
            graph_of(vec![
                (101, [None; 6]),
                (101, [None; 6]),
                (401, [Some(0), Some(1), Some(3), None, None, None]),
                (201, [None; 6]),
            ])
        };
        let mut g = make();
        let mut lostal = vec![None; 4];
        let mut rng = Ran2::legacy();
        do_reaction(&mut g, &mut lostal, 2, 4, &mut rng).unwrap();
        assert_eq!(st(&g, 2), 410);
        // One of the two Al gained; lostal records which.
        let chosen = lostal[2].expect("R4 sets lostal");
        assert!(chosen == 0 || chosen == 1);
        assert_eq!(st(&g, chosen), 102);
        do_reaction(&mut g, &mut lostal, 2, 5, &mut rng).unwrap();
        assert_eq!(g, make());
        assert_eq!(lostal[2], None);
    }

    #[test]
    fn adsorb_then_desorb_al_round_trips() {
        // Empty Al(0) with six 400 oxygens around it.
        let make = || {
            graph_of(vec![
                (100, [Some(1), Some(2), Some(3), Some(4), Some(5), Some(6)]),
                (400, [None; 6]),
                (400, [None; 6]),
                (400, [None; 6]),
                (400, [None; 6]),
                (400, [None; 6]),
                (400, [None; 6]),
            ])
        };
        let mut g = make();
        let mut lostal = vec![None; 7];
        let mut rng = Ran2::legacy();
        do_reaction(&mut g, &mut lostal, 0, 16, &mut rng).unwrap();
        assert_eq!(st(&g, 0), 107);
        assert_eq!(st(&g, 1), 409); // 400 -> 409
        do_reaction(&mut g, &mut lostal, 0, 20, &mut rng).unwrap();
        assert_eq!(g, make()); // 409 -> 400, Al -> 100
    }

    #[test]
    fn adsorb_al_rejects_an_unexpected_oxygen_state() {
        let mut g = graph_of(vec![
            (100, [Some(1), Some(1), Some(1), Some(1), Some(1), Some(1)]),
            (301, [None; 6]), // not a valid Al-neighbor oxygen
        ]);
        let mut lostal = vec![None; 2];
        let mut rng = Ran2::legacy();
        let e = do_reaction(&mut g, &mut lostal, 0, 16, &mut rng);
        assert!(matches!(e, Err(KaolError::InvalidAdsorbState { .. })));
    }
}
