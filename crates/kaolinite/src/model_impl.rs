//! `impl Model for Kaolinite` — the kaolinite chemistry bolted onto the
//! engine's [`Model`] seam. This is where the C++'s
//! `evtlist.cpp` (event enumeration), `envrn.cpp` (allow/bucket), and
//! `actions.cpp` (mutation) meet, behind the one trait the engine speaks.
//!
//! The C++ fused these three files by direct call: `CreateEventList` calls
//! `environment->IsActive` and `CheckEnv`; `DoEvent` calls
//! `actions.DoReaction`. Here [`events_at`](Kaolinite::events_at) is the
//! per-site body of `CreateEventList` and [`apply`](Kaolinite::apply) is
//! `DoReaction` — same logic, cut along the engine interface so the engine
//! never names a state code.

use kmc_engine::{Model, ProposedEvent, Rng, SiteGraph, SiteId};

use crate::build::Structure;
use crate::environment::{check_env, is_active};
use crate::mechanisms::do_reaction;
use crate::reactions::{N_300, N_400, N_DES, N_HYD, ReactionSet};
use crate::state::State;

/// Errors the kaolinite model can raise mid-simulation — the value form of
/// the C++'s scattered `std::cerr` + abort. All three are "the run cannot
/// continue" conditions the legacy code also stops on.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum KaolError {
    /// `CheckEnv` produced an index that is negative or ≥ the reaction's rate
    /// count (C++ `evtlist.cpp`: prints "invalid environment" and returns a
    /// null event list, aborting the run).
    InvalidEnvironment {
        /// The site whose environment could not be classified.
        site: SiteId,
        /// Its state code (for diagnosis).
        state: i32,
        /// The out-of-range bucket index `CheckEnv` returned.
        env: i32,
        /// How many buckets the reaction actually has.
        nrates: usize,
    },
    /// An adsorption mechanism met an oxygen in a state its transition map
    /// doesn't cover (C++ `AdsorbAl`/`AdsorbSi`: prints "invalid state" and
    /// returns false).
    InvalidAdsorbState {
        /// The oxygen neighbor with the unexpected state.
        site: SiteId,
        /// Its state code.
        state: i32,
    },
    /// `do_reaction` was handed a reaction id it does not implement — only
    /// the dead diffusion range (24–27) can reach this, and only via a bug
    /// (they are never proposed). Faithful "impossible default" signal.
    UnsupportedReaction(u16),
}

impl std::fmt::Display for KaolError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            KaolError::InvalidEnvironment {
                site,
                state,
                env,
                nrates,
            } => write!(
                f,
                "invalid environment: site {site} state {state} env {env} nrates {nrates}"
            ),
            KaolError::InvalidAdsorbState { site, state } => {
                write!(f, "invalid state in adsorb: site {site} state {state}")
            }
            KaolError::UnsupportedReaction(r) => write!(f, "unsupported reaction id {r}"),
        }
    }
}

impl std::error::Error for KaolError {}

/// The kaolinite model: the reaction tables plus the two kaolinite-specific
/// per-site arrays the engine's `Site` deliberately does not carry (design
/// doc §3, option (a)). The `graph` itself lives with the engine and is
/// passed to every trait method; `pair`/`lostal` ride here because they are
/// the *model's* bookkeeping, not generic-KMC data.
///
/// \[IDIOM\] Splitting mutable state across the `&mut self` model and the
/// `&mut graph` argument is what lets the borrow checker bless the whole
/// step: `apply` mutates `self.lostal` and `graph` as two independent `&mut`
/// borrows, never aliasing. A fat C++-style site struct (state + pair +
/// lostal in one array) would force every mutation to borrow the lot.
#[derive(Debug, Clone)]
pub struct Kaolinite {
    /// Temperature, chemical potentials, and the 24 reaction rate tables.
    pub rxn: ReactionSet,
    /// Double-bridge partner per site (built by `find_pairs`, read-only here).
    pub pair: Vec<Option<SiteId>>,
    /// The lost-Al record per site (mutated by R4/R5/R8/R9 and adsorb/desorb).
    pub lostal: Vec<Option<SiteId>>,
}

impl Kaolinite {
    /// Consume a built [`Structure`] and a [`ReactionSet`] into the graph the
    /// engine will own plus the model that drives it. This is the handoff
    /// from the deterministic structural build (M3) to the dynamics (M4–M6).
    ///
    /// Returns the graph separately (not inside the model) precisely because
    /// [`kmc_engine::step`] wants `&mut graph` and `&mut model` as distinct
    /// borrows — see the split-state idiom on [`Kaolinite`].
    pub fn from_structure(structure: Structure, rxn: ReactionSet) -> (SiteGraph<State>, Self) {
        let model = Kaolinite {
            rxn,
            pair: structure.pair,
            lostal: structure.lostal,
        };
        (structure.graph, model)
    }
}

impl Model for Kaolinite {
    type State = State;
    type Error = KaolError;

    /// Enumerate the events possible at `site` — the per-site body of
    /// `evtlist.cpp::CreateEventList`.
    ///
    /// The C++ picks a *window* of reaction ids from the site's class, then
    /// keeps those whose reactant state matches and that `IsActive` allows,
    /// resolving each one's rate through `CheckEnv`:
    ///
    /// | site state | window (C++ `lo..hi`) | reactions |
    /// |---|---|---|
    /// | empty O (`>200`, `%100==0`) | — | skipped |
    /// | `>500` | `N400..NHYD` (14..16) | R14/R15 |
    /// | `>400` | `N300..N400` (2..14) | R2..R13 |
    /// | `>300` | `0..N300` (0..2) | R0/R1 |
    /// | else (cations, empty Al/Si) | `NHYD..NRXN` (16..28) | ads/desorb |
    ///
    /// We cap the last window at [`N_DES`] (24) instead of 28: reactions
    /// 24–27 are diffusion, which `is_active` always rejects, so iterating
    /// them is a guaranteed no-op — omitting the dead tail changes nothing
    /// observable (spec B6) and keeps the port from indexing tables it does
    /// not hold.
    fn events_at(
        &self,
        graph: &SiteGraph<State>,
        site: SiteId,
        out: &mut Vec<ProposedEvent>,
    ) -> Result<(), KaolError> {
        let s = graph.sites[site].state.0;

        // Skip empty non-Al oxygens/Si (C++: `state%100==0 && state>200`).
        // Empty Al (100) and empty Si (200) are NOT skipped — 200 is not
        // `> 200` — so they still reach the adsorption window below.
        if s % 100 == 0 && s > 200 {
            return Ok(());
        }

        let (lo, hi) = if s > 500 {
            (N_400, N_HYD)
        } else if s > 400 {
            (N_300, N_400)
        } else if s > 300 {
            (0, N_300)
        } else {
            (N_HYD, N_DES)
        };

        for i in lo..hi {
            let rxn = &self.rxn.reactions[i];
            if graph.sites[site].state == rxn.reactant && is_active(graph, site, i as u16) {
                let env = check_env(graph, &self.pair, site);
                if env < 0 || env as usize >= rxn.rates.len() {
                    return Err(KaolError::InvalidEnvironment {
                        site,
                        state: s,
                        env,
                        nrates: rxn.rates.len(),
                    });
                }
                out.push(ProposedEvent {
                    site,
                    rxn: i as u16,
                    rate: rxn.rates[env as usize],
                });
            }
        }
        Ok(())
    }

    /// Apply a chosen event — `actions.cpp::DoReaction`, threaded with the
    /// engine's shared RNG for the R4/R9 proton coin.
    fn apply(
        &mut self,
        graph: &mut SiteGraph<State>,
        ev: &ProposedEvent,
        rng: &mut dyn Rng,
    ) -> Result<(), KaolError> {
        do_reaction(graph, &mut self.lostal, ev.site, ev.rxn, rng)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::build::{
        LatticeParams, create_lattice, find_pairs, populate_solid, terminate_lattice,
        terminate_surface,
    };
    use crate::cell::{CellSite, NeighborTemplate, UnitCell};
    use crate::reactions::Reaction;

    // A tiny hand-made reaction set is awkward; instead reuse the crate's own
    // types with minimal tables where a specific event is exercised. Most of
    // the model is covered end-to-end by the integration parity test; these
    // unit tests pin the window-selection and error-surfacing seams.

    fn empty_rxn_set() -> ReactionSet {
        // 24 reactions with the golden reactants but 1-bucket flat rates, so
        // window selection can be tested without a full data.rxn.
        let reactants = [
            301, 302, 401, 402, 401, 410, 402, 403, 410, 403, 406, 407, 404, 40100, 501, 502, 100,
            200, 100, 200, 107, 199, 205, 299,
        ];
        let reactions = reactants
            .iter()
            .map(|&r| Reaction {
                reactant: State(r),
                info: 0,
                rates: vec![1.0],
            })
            .collect();
        ReactionSet {
            temperature: 8000.0,
            dm_si: -1.0,
            dm_al: -1.0,
            reactions,
        }
    }

    #[test]
    fn empty_oxygen_proposes_nothing_but_empty_al_proposes_adsorption() {
        // Two isolated sites: empty Al (100, has an occupied neighbor) and
        // empty 400 oxygen (skipped entirely).
        let graph = SiteGraph {
            sites: vec![
                kmc_engine::Site {
                    state: State(100),
                    nbr: [Some(1), None, None, None, None, None],
                },
                kmc_engine::Site {
                    state: State(205),
                    nbr: [Some(0), None, None, None, None, None],
                },
                kmc_engine::Site {
                    state: State(400),
                    nbr: [None; 6],
                },
            ],
        };
        let model = Kaolinite {
            rxn: empty_rxn_set(),
            pair: vec![None; 3],
            lostal: vec![None; 3],
        };
        let mut out = Vec::new();
        model.events_at(&graph, 2, &mut out).unwrap(); // empty 400 O: skipped
        assert!(out.is_empty());
        model.events_at(&graph, 0, &mut out).unwrap(); // empty Al: adsorb (rxn 16)
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].rxn, 16);
    }

    #[test]
    fn from_structure_splits_graph_and_model() {
        // A minimal all-Al cell just to exercise the handoff plumbing.
        let t = |n, a, b| NeighborTemplate { n, a, b, c: 0 };
        let none = t(-1, 0, 0);
        let uc = UnitCell {
            a: 1.0,
            b: 1.0,
            c: 1.0,
            alpha: 0.0,
            beta: 0.0,
            gamma: 0.0,
            sites: vec![CellSite {
                x: 0.0,
                y: 0.0,
                z: 0.0,
                n: 0,
                state: State(100),
                nbr: [none; 6],
            }],
        };
        let params = LatticeParams {
            a_cells: 2,
            b_cells: 2,
            surface_plane: 0,
        };
        let mut structure = create_lattice(&uc, params);
        find_pairs(&mut structure);
        populate_solid(&mut structure, -1.0, -1.0);
        terminate_surface(&mut structure);
        terminate_lattice(&mut structure);
        let n = structure.graph.len();
        let (graph, model) = Kaolinite::from_structure(structure, empty_rxn_set());
        assert_eq!(graph.len(), n);
        assert_eq!(model.pair.len(), n);
        assert_eq!(model.lostal.len(), n);
    }
}
