//! The nth-nearest-neighbor environment — the scientific heart of the model.
//! Port of `envrn.cpp` (`Environment::IsActive` and `CheckEnv` + its five
//! `Check{100..500}` helpers).
//!
//! Two questions per candidate (site, reaction):
//!
//! * [`is_active`] — is the reaction *allowed* here? For forward hydrolysis
//!   this is the "at the surface?" test: does a **second**-neighbor bear a
//!   hydrolyzed/OH state? That double loop over `nbr[i].nbr[j]` is the
//!   concrete "reach two hops out" that distinguishes this from a mean-field
//!   rate law (spec A6.1).
//! * [`check_env`] — *which* environment bucket (→ which rate)? Each class
//!   folds specific neighbor / second-neighbor / `pair`-partner states into a
//!   single integer index into the reaction's `rates[]` (spec A6.2).
//!
//! # Faithfulness over cleanliness
//!
//! These functions are transliterated with their integer arithmetic and
//! short-circuits intact, **including** the documented warts (B3: the
//! `Check100`/`Check200` indices don't match the header's own formula) — the
//! rates are flat in the sample data (B4), so the exact bucket value is inert
//! today, but the *in-range-ness* of the returned index gates whether the run
//! aborts, and the whole point of parity is to reproduce that gate. A "-1"
//! return is the C++ "ran into a lattice edge / invalid environment"
//! signal, which aborts the legacy run; here it surfaces as
//! [`EnvError::InvalidEnvironment`](crate::model_impl::KaolError) at the call
//! site (the library stays silent; no `std::cerr`).
//!
//! Where the C++ reads `sites[nbr]` without a `nbr >= 0` guard (`Check100`,
//! `Check200`), a missing neighbor would be undefined behavior — it never
//! happens on a correctly terminated lattice (an interior occupied cation has
//! all its neighbors; boundary ones are frozen `EDGE`), so we return -1 on
//! the impossible `None` rather than reproduce the UB.

use kmc_engine::{SiteGraph, SiteId};

use crate::reactions::N_HYD;
use crate::state::State;

/// `Environment::IsActive` — may reaction `rxn` fire at `site`?
///
/// The five cases (spec A6.1), in the C++'s order:
/// * forward hydrolysis (even, `< NHYD`): only "at the surface" — some
///   second-neighbor is in {303, 404, 405, 406, 408, 409};
/// * reverse hydrolysis (odd, `< NHYD`): always;
/// * adsorption (`rxn == 16` or `19` only): needs ≥1 occupied non-EDGE
///   neighbor;
/// * desorption (`rxn == 20` or `22` only): always;
/// * everything else (the cross-cation 17/18/21/23 and diffusion 24–27):
///   never — this is how those reactions are switched off (spec B6).
pub fn is_active(graph: &SiteGraph<State>, site: SiteId, rxn: u16) -> bool {
    let n_hyd = N_HYD as u16;
    // C++ `rxn % 2 == 0` (even index = forward hydrolysis); `is_multiple_of`
    // is the clippy-preferred spelling of the same test.
    if rxn < n_hyd && rxn.is_multiple_of(2) {
        // Forward hydrolysis: the surface-reachability double loop — with a
        // WART that is *load-bearing for parity* (spec §C2a).
        //
        // # WART (new — the `nbr[6]` out-of-bounds phantom)
        //
        // The C++ inner loop is
        //   `for (j=0; (nbr2 = sites[nbr].nbr[j]) >= 0 && j < 6 && !result; j++)`
        // — the `nbr[j]` read happens *before* the `j < 6` guard in the
        // `&&` chain. When a neighbor `nbr` has all six neighbor slots
        // filled, the loop reaches `j == 6` and evaluates `sites[nbr].nbr[6]`,
        // one element past the fixed `nbr[6]` array. That read is undefined
        // behavior; under the golden build (`g++ -O3 -ffast-math`) it
        // manifests, deterministically, as **`result` becoming TRUE** whenever
        // the inner loop reaches `j == 6`. (Reverse-engineered against the
        // golden trajectory: the rule "a real 2nd-neighbor is hydrolyzed OR a
        // neighbor is fully six-coordinated" reproduces every `IsActive`
        // decision the golden binary makes — verified across all 20,000 steps
        // by the parity gate. A bounds-clean loop drops step-0's event count
        // from 660 to 180 and desynchronizes immediately.)
        //
        // So faithful behavior is: for each neighbor `nbr` of the site, the
        // reaction is active if any of `nbr`'s *present* neighbors bears a
        // hydrolyzed state {303,404,405,406,408,409}, **or** `nbr` has all six
        // neighbor slots present (Al sites always do) — the latter is the
        // phantom. Since topology is fixed for the whole run, the phantom
        // condition is static per neighbor; only the real-state test is
        // dynamic. This is not in spec Part B — discovered during the M6
        // parity chase; a candidate for a spec 02 addendum and a REFORM_PLAN
        // entry (bounds-clean is the corrected behavior).
        //
        // The site's own outer loop also lacks an `i < 6` bound, but every
        // forward-hydrolysis reactant is an O-class state (≤3 neighbors), so
        // the outer loop always meets a `None` within six slots and never
        // reaches its own phantom — reproduced by stopping at the first
        // `None` below.
        for &nbr in graph.sites[site].nbr.iter() {
            let Some(nbr) = nbr else { break };
            if inner_surface_active(graph, nbr) {
                return true;
            }
        }
        false
    } else if rxn < n_hyd {
        // Reverse hydrolysis: always allowed.
        true
    } else if rxn == 16 || rxn == 19 {
        // Adsorption to the correct site: at least one occupied, non-EDGE
        // neighbor (can't nucleate on nothing).
        for &nbr in graph.sites[site].nbr.iter() {
            let Some(nbr) = nbr else { break };
            let st = graph.sites[nbr].state;
            if st.is_occupied() && !st.is_edge() {
                return true;
            }
        }
        false
    } else if rxn == 20 || rxn == 22 {
        // Desorption: always allowed.
        true
    } else {
        // Cross-cation (17/18/21/23) and diffusion (24–27): disabled.
        false
    }
}

/// The inner loop of the forward-hydrolysis surface test, over one neighbor
/// `nbr` of the reacting site — **including the `nbr[6]` OOB phantom** (see
/// the big WART note in [`is_active`]).
///
/// Returns TRUE if any of `nbr`'s present neighbors is in the hydrolyzed
/// matchset, OR if all six of `nbr`'s neighbor slots are present (the loop
/// would reach `j == 6` and the golden binary's out-of-bounds read forces the
/// result true). If a `None` slot is met before the sixth, the loop stops
/// with no phantom.
fn inner_surface_active(graph: &SiteGraph<State>, nbr: SiteId) -> bool {
    for &nbr2 in graph.sites[nbr].nbr.iter() {
        match nbr2 {
            Some(nbr2) => {
                if matches!(graph.sites[nbr2].state.0, 303 | 404 | 405 | 406 | 408 | 409) {
                    return true;
                }
            }
            // A missing slot within the six: the C++ loop stops here, before
            // the phantom `j == 6` read. No trigger from this neighbor.
            None => return false,
        }
    }
    // All six slots were present → the loop reached `j == 6` → the OOB
    // phantom fires TRUE in the golden build.
    true
}

/// `Environment::CheckEnv` — the environment bucket (rate index) for `site`.
///
/// Dispatches by class to the `check_*` helpers. Note the class-4 split: 404
/// and 405 (the two Al-OH-Al-like 400 states) route to `check500`, every
/// other 400 to `check400` — a quirk of the original worth preserving
/// exactly. Returns `-1` for an invalid environment (the C++ abort signal).
pub fn check_env(graph: &SiteGraph<State>, pair: &[Option<SiteId>], site: SiteId) -> i32 {
    let st = graph.sites[site].state.0;
    match st / 100 {
        1 => {
            if st == 100 {
                0
            } else {
                check100(graph, site)
            }
        }
        2 => {
            if st == 200 {
                0
            } else {
                check200(graph, site)
            }
        }
        3 => check300(graph, site),
        4 => {
            if st == 404 || st == 405 {
                check500(graph, pair, site)
            } else {
                check400(graph, pair, site)
            }
        }
        5 => check500(graph, pair, site),
        _ => -1,
    }
}

/// `Check100` — environment index for an occupied Al site.
///
/// WART (spec B3): the header documents `index = x*3 + y`, but the code
/// returns `(x + y) / 2`. Preserved verbatim; inert while rates are flat
/// (B4), but pinned so a "fix" can't sneak in unnoticed.
fn check100(graph: &SiteGraph<State>, site: SiteId) -> i32 {
    let (mut x, mut y) = (0, 0);
    for j in 0..6 {
        let Some(nbr) = graph.sites[site].nbr[j] else {
            return -1;
        };
        let ns = graph.sites[nbr].state;
        if ns.is_edge() {
            return -1;
        }
        match ns.0 {
            502 => x += 1,
            403 | 405 | 407 | 409 | 410 => y += 1,
            _ => {}
        }
    }
    (x + y) / 2
}

/// `Check200` — environment index for an occupied Si site.
///
/// WART (spec B3): header says `index = x*4 + y`; code returns `x + y` with
/// `x` starting at 1 and cleared to 0 by any adjacent 408. Preserved.
fn check200(graph: &SiteGraph<State>, site: SiteId) -> i32 {
    let mut y = 0;
    let mut x = 1;
    for j in 0..4 {
        let Some(nbr) = graph.sites[site].nbr[j] else {
            return -1;
        };
        let ns = graph.sites[nbr].state;
        if ns.is_edge() {
            return -1;
        }
        match ns.0 {
            408 => x = 0,
            302 => y += 1,
            _ => {}
        }
    }
    x + y
}

/// `Check300` — environment index for a Si-O-Si oxygen. Reaches through both
/// Si to their first oxygen neighbor (`Si-O-Al` is Si's `nbr[0]`), so the
/// bucket reflects second-neighbor chemistry (spec A6.2).
fn check300(graph: &SiteGraph<State>, site: SiteId) -> i32 {
    let (Some(si1), Some(si2)) = (graph.sites[site].nbr[0], graph.sites[site].nbr[1]) else {
        return -1;
    };
    let (sio1, sio2) = (graph.sites[si1].nbr[0], graph.sites[si2].nbr[0]);
    // Any missing second-neighbor oxygen, or any EDGE along the chain, aborts.
    let (Some(sio1), Some(sio2)) = (sio1, sio2) else {
        return -1;
    };
    if graph.sites[si1].state.is_edge()
        || graph.sites[si2].state.is_edge()
        || graph.sites[sio1].state.is_edge()
        || graph.sites[sio2].state.is_edge()
    {
        return -1;
    }

    let mut x = 0;
    let broke = |s: i32| matches!(s, 402 | 403 | 407 | 408);
    if broke(graph.sites[sio1].state.0) {
        x += 1;
    }
    if broke(graph.sites[sio2].state.0) {
        x += 1;
    }
    let mut y = (graph.sites[si1].state.0 - 201) + (graph.sites[si2].state.0 - 201) - x;
    if graph.sites[site].state.0 > 301 {
        y -= 2;
    }
    5 * x + y
}

/// `Check400` — environment index for a Si-O-Al₂ oxygen (the richest class).
/// Reaches through the Si's other oxygens (its `nbr[1]`, `nbr[2]`) for the
/// broken-300 count and through both Al plus the `pair` partner for the
/// broken-400/500 count.
fn check400(graph: &SiteGraph<State>, pair: &[Option<SiteId>], site: SiteId) -> i32 {
    let (Some(al1), Some(al2), Some(si)) = (
        graph.sites[site].nbr[0],
        graph.sites[site].nbr[1],
        graph.sites[site].nbr[2],
    ) else {
        return -1;
    };
    let Some(p) = pair[site] else {
        return -1;
    };
    if graph.sites[al1].state.is_edge()
        || graph.sites[al2].state.is_edge()
        || graph.sites[si].state.is_edge()
        || graph.sites[p].state.is_edge()
    {
        return -1;
    }

    let mut x = 0;
    for i in 1..3 {
        let Some(o1) = graph.sites[si].nbr[i] else {
            return -1;
        };
        if graph.sites[o1].state.0 > 301 {
            x += 1;
        }
    }

    let mut y = (graph.sites[al1].state.0 - 101) + (graph.sites[al2].state.0 - 101);
    let ss = graph.sites[site].state.0;
    let ps = graph.sites[p].state.0;
    if ss == 403 || ss == 405 {
        y -= 2;
    }
    if ss == 407 || ss == 410 || ps == 502 {
        y -= 1;
    }
    if ps == 502 {
        y -= 1;
    }

    if ss == 406 || ss == 407 {
        x * 5 + y
    } else {
        x * 10 + y
    }
}

/// `Check500` — environment index for an Al-OH-Al oxygen (and the 404/405
/// 400-states routed here by [`check_env`]). Uses the `pair` partner's state
/// to detect a hydrolyzed sibling bridge.
fn check500(graph: &SiteGraph<State>, pair: &[Option<SiteId>], site: SiteId) -> i32 {
    let (Some(al1), Some(al2)) = (graph.sites[site].nbr[0], graph.sites[site].nbr[1]) else {
        return -1;
    };
    let Some(p) = pair[site] else {
        return -1;
    };
    if graph.sites[al1].state.is_edge()
        || graph.sites[al2].state.is_edge()
        || graph.sites[p].state.is_edge()
    {
        return -1;
    }

    let ps = graph.sites[p].state.0;
    let x = if matches!(ps, 502 | 403 | 405 | 410) {
        1
    } else {
        0
    };
    let mut y = (graph.sites[al1].state.0 - 101) + (graph.sites[al2].state.0 - 101);
    if matches!(ps, 502 | 403 | 405) {
        y -= 1;
    }
    let ss = graph.sites[site].state.0;
    if ss == 502 || ss == 405 {
        y -= 2;
    }
    x * 9 + y
}

#[cfg(test)]
mod tests {
    use super::*;
    use kmc_engine::Site;

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
    fn reverse_hydrolysis_and_desorption_are_always_active() {
        let g = graph_of(vec![(302, [None; 6])]);
        assert!(is_active(&g, 0, 1)); // odd hydrolysis
        assert!(is_active(&g, 0, 20)); // desorb Al
        assert!(is_active(&g, 0, 22)); // desorb Si
    }

    #[test]
    fn cross_and_diffusion_reactions_are_never_active() {
        let g = graph_of(vec![(401, [None; 6])]);
        for rxn in [17u16, 18, 21, 23, 24, 25, 26, 27] {
            assert!(!is_active(&g, 0, rxn), "rxn {rxn} should be disabled");
        }
    }

    #[test]
    fn forward_hydrolysis_needs_a_hydrolyzed_second_neighbor() {
        // site 0 -> nbr 1 -> nbr 2. Only when site 2 bears an OH-state
        // (e.g. 303) is forward hydrolysis (rxn 0) active at site 0.
        let inactive = graph_of(vec![
            (301, [Some(1), None, None, None, None, None]),
            (200, [Some(0), Some(2), None, None, None, None]),
            (200, [Some(1), None, None, None, None, None]),
        ]);
        assert!(!is_active(&inactive, 0, 0));

        let active = graph_of(vec![
            (301, [Some(1), None, None, None, None, None]),
            (200, [Some(0), Some(2), None, None, None, None]),
            (303, [Some(1), None, None, None, None, None]),
        ]);
        assert!(is_active(&active, 0, 0));
    }

    #[test]
    fn adsorption_needs_an_occupied_nonedge_neighbor() {
        let none = graph_of(vec![
            (100, [Some(1), None, None, None, None, None]),
            (400, [Some(0), None, None, None, None, None]), // empty
        ]);
        assert!(!is_active(&none, 0, 16));
        let edge_only = graph_of(vec![
            (100, [Some(1), None, None, None, None, None]),
            (9, [Some(0), None, None, None, None, None]), // EDGE, excluded
        ]);
        assert!(!is_active(&edge_only, 0, 16));
        let ok = graph_of(vec![
            (100, [Some(1), None, None, None, None, None]),
            (401, [Some(0), None, None, None, None, None]),
        ]);
        assert!(is_active(&ok, 0, 16));
    }

    #[test]
    fn check_env_edge_neighbor_returns_minus_one() {
        // An occupied Al (101) whose neighbor is EDGE → Check100 bails to -1.
        let g = graph_of(vec![
            (101, [Some(1), Some(1), Some(1), Some(1), Some(1), Some(1)]),
            (9, [None; 6]),
        ]);
        assert_eq!(check_env(&g, &[None, None], 0), -1);
    }

    #[test]
    fn empty_cation_buckets_are_zero() {
        let g = graph_of(vec![(100, [None; 6]), (200, [None; 6])]);
        assert_eq!(check_env(&g, &[None, None], 0), 0); // empty Al
        assert_eq!(check_env(&g, &[None, None], 1), 0); // empty Si
    }
}
