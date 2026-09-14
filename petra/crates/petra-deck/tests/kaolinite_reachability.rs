use std::collections::HashSet;
use std::path::PathBuf;

use petra_core::{
    BiasedCtmc, CtmcAdvance, Fired, StateId, StepCtx, StepOutcome, Stop, UpdateStrategy,
};

fn deck_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/kaolinite-approx.toml")
}

struct FireExact {
    site: usize,
    reaction: u16,
}

impl UpdateStrategy for FireExact {
    fn step(&mut self, ctx: &mut StepCtx<'_>) -> Result<StepOutcome, Stop> {
        if !ctx.apply.enabled_rules(self.site).contains(&self.reaction) {
            return Err(Stop::NoEvents);
        }
        ctx.apply
            .apply_transition(self.site, self.reaction, ctx.rng)?;
        Ok(StepOutcome {
            fired: vec![Fired {
                step: 0,
                time: 0.0,
                site: self.site,
                reaction: self.reaction,
            }],
            dt: 1.0,
        })
    }
}

struct ForwardClosure<'a> {
    forward: &'a [u16],
    releases: [u16; 2],
    original: &'a HashSet<usize>,
}

impl UpdateStrategy for ForwardClosure<'_> {
    fn step(&mut self, ctx: &mut StepCtx<'_>) -> Result<StepOutcome, Stop> {
        for &reaction in self.releases.iter().chain(self.forward.iter()) {
            for site in 0..ctx.lattice.len() {
                if !ctx.apply.enabled_rules(site).contains(&reaction) {
                    continue;
                }
                if self.releases.contains(&reaction) && !self.original.contains(&site) {
                    continue;
                }
                if ctx.apply.apply_transition(site, reaction, ctx.rng).is_ok() {
                    return Ok(StepOutcome {
                        fired: vec![Fired {
                            step: 0,
                            time: 0.0,
                            site,
                            reaction,
                        }],
                        dt: 1.0,
                    });
                }
            }
        }
        Err(Stop::NoEvents)
    }
}

fn state_id(deck: &petra_deck::CompiledDeck, name: &str) -> StateId {
    StateId(
        deck.state_names
            .iter()
            .position(|candidate| candidate == name)
            .expect("state exists") as u16,
    )
}

fn exercise_original_lattice_release(kind: &str, site: usize, eligible: &str, reaction: &str) {
    let deck = petra_deck::load(deck_path()).expect("production deck compiles");
    let mut engine = deck.build_engine(Some(42)).expect("initial state builds");
    let original = engine.lattice.states[site];
    assert_eq!(
        deck.state_occupants[original.0 as usize].as_deref(),
        Some(kind)
    );
    assert!(!engine.lattice.frozen[site]);

    let eligible = state_id(&deck, eligible);
    engine.lattice.states[site] = eligible;
    engine
        .set_temperature(deck.temperature)
        .expect("refresh seeded event table");
    let reaction = deck
        .reactions
        .iter()
        .position(|candidate| candidate.name == reaction)
        .expect("release reaction exists") as u16;
    let outcome = engine
        .step_with(&mut FireExact { site, reaction })
        .expect("seeded release fires");

    assert_eq!(outcome.fired.len(), 1);
    assert_eq!(outcome.fired[0].site, site);
    assert_eq!(outcome.fired[0].reaction, reaction);
    assert!(engine.last_changes().contains(&(
        site,
        eligible,
        state_id(&deck, &format!("{kind}.empty"))
    )));
}

#[test]
fn reachable_si_and_al_release_transitions_preserve_original_site_lineage() {
    // These are the first reachable Al/Si site ids in the independently
    // generated A9b production path census. Each was occupied at step zero.
    exercise_original_lattice_release("Al", 28, "Al.l6", "desorb-al");
    exercise_original_lattice_release("Si", 32, "Si.oh4", "desorb-si");
}

#[test]
fn live_forward_closure_proves_260_paths_and_60_boundary_no_gos() {
    let deck = petra_deck::load(deck_path()).expect("production deck compiles");
    let mut engine = deck.build_engine(Some(42)).expect("initial state builds");
    let original: HashSet<_> = engine
        .lattice
        .states
        .iter()
        .enumerate()
        .filter(|(site, state)| {
            !engine.lattice.frozen[*site]
                && matches!(
                    deck.state_occupants[state.0 as usize].as_deref(),
                    Some("Al" | "Si")
                )
        })
        .map(|(site, _)| site)
        .collect();
    let reaction_id = |name: &str| {
        deck.reactions
            .iter()
            .position(|reaction| reaction.name == name)
            .expect("reaction exists") as u16
    };
    let forward: Vec<_> = [
        "R0-sio-si-hydrolysis",
        "R2-sioal2-si-hydrolysis",
        "R4-sioal2-al-hydrolysis",
        "R6-second-stage-hydrolysis",
        "R8a-sioha-hydrolysis",
        "R8b-sioha-hydrolysis",
        "R10-sial-hydrolysis",
        "R12-albr-hydrolysis",
        "R14-alohal-hydrolysis",
    ]
    .iter()
    .map(|name| reaction_id(name))
    .collect();
    let releases = [reaction_id("desorb-al"), reaction_id("desorb-si")];
    let mut released = HashSet::new();

    for _ in 0..2_000 {
        let mut closure = ForwardClosure {
            forward: &forward,
            releases,
            original: &original,
        };
        match engine.step_with(&mut closure) {
            Ok(outcome) => {
                let fired = outcome.fired[0];
                if releases.contains(&fired.reaction) {
                    released.insert(fired.site);
                }
            }
            Err(Stop::NoEvents) => break,
            Err(error) => panic!("closure failed: {error}"),
        }
    }

    assert_eq!(original.len(), 320);
    assert_eq!(released.len(), 260);
    assert!(
        released.contains(&0),
        "initial-local analysis missed this path"
    );
    let expected_no_go: HashSet<_> = (0..20)
        .flat_map(|cell| [78 * cell + 3, 78 * cell + 4, 78 * cell + 7])
        .collect();
    assert_eq!(
        original
            .difference(&released)
            .copied()
            .collect::<HashSet<_>>(),
        expected_no_go
    );

    assert_eq!(engine.lattice.states[3], state_id(&deck, "Al.l4"));
    assert_eq!(engine.lattice.states[4], state_id(&deck, "Si.oh2"));
    assert_eq!(engine.lattice.states[7], state_id(&deck, "Si.oh1"));
}

#[test]
fn reservoir_readsorption_on_a_vacated_lattice_site_counts_only_once() {
    let deck = petra_deck::load(deck_path()).expect("production deck compiles");
    let mut engine = deck.build_engine(Some(42)).expect("initial state builds");
    let site = 28;
    let reaction_id = |name: &str| {
        deck.reactions
            .iter()
            .position(|reaction| reaction.name == name)
            .expect("reaction exists") as u16
    };
    let mut lattice_ion_present = HashSet::from([site]);
    let mut counted_releases = 0;

    engine.lattice.states[site] = state_id(&deck, "Al.l6");
    engine.set_temperature(deck.temperature).unwrap();
    engine
        .step_with(&mut FireExact {
            site,
            reaction: reaction_id("desorb-al"),
        })
        .expect("original lattice ion desorbs");
    counted_releases += usize::from(lattice_ion_present.remove(&site));

    engine
        .step_with(&mut FireExact {
            site,
            reaction: reaction_id("adsorb-al"),
        })
        .expect("reservoir ion readsorbs onto the same site");
    engine
        .step_with(&mut FireExact {
            site,
            reaction: reaction_id("desorb-al"),
        })
        .expect("readsorbed reservoir ion desorbs");
    counted_releases += usize::from(lattice_ion_present.remove(&site));

    assert_eq!(counted_releases, 1);
}

#[test]
fn identical_seed_on_an_initially_empty_site_is_not_lattice_dissolution() {
    let deck = petra_deck::load(deck_path()).expect("production deck compiles");
    let mut engine = deck.build_engine(Some(42)).expect("initial state builds");
    let empty_si = state_id(&deck, "Si.empty");
    let reservoir_site = engine
        .lattice
        .states
        .iter()
        .enumerate()
        .find_map(|(site, &state)| {
            (state == empty_si && !engine.lattice.frozen[site]).then_some(site)
        })
        .expect("nonfrozen reservoir site exists");
    let initially_lattice = deck.state_occupants[empty_si.0 as usize].is_some();

    engine.lattice.states[reservoir_site] = state_id(&deck, "Si.oh4");
    engine
        .set_temperature(deck.temperature)
        .expect("refresh seeded event table");
    let reaction = deck
        .reactions
        .iter()
        .position(|candidate| candidate.name == "desorb-si")
        .expect("release reaction exists") as u16;
    engine
        .step_with(&mut FireExact {
            site: reservoir_site,
            reaction,
        })
        .expect("reservoir release fires physically");

    assert!(!initially_lattice, "event must remain reservoir-origin");
}

#[test]
fn finite_horizon_biased_ctmc_executes_on_the_production_deck() {
    let deck = petra_deck::load(deck_path()).expect("production deck compiles");
    let mut engine = deck
        .build_engine(Some(90401))
        .expect("initial state builds");
    let overrides = deck
        .reactions
        .iter()
        .enumerate()
        .filter(|(_, reaction)| {
            reaction.name.contains("hydrolysis") || reaction.name.starts_with("desorb-")
        })
        .map(|(reaction, _)| (reaction as u16, 1.0e6));
    let mut sampler =
        BiasedCtmc::new(deck.reactions.len(), overrides).expect("finite positive bias");
    let horizon = 1.0e-6;
    let mut releases = 0_u64;

    for _ in 0..256 {
        match engine
            .advance_biased_ctmc_until(&mut sampler, horizon)
            .expect("biased production advance")
        {
            CtmcAdvance::Fired(event) => {
                let name = &deck.reactions[event.reaction as usize].name;
                releases += u64::from(name == "desorb-al" || name == "desorb-si");
            }
            CtmcAdvance::Deadline { time } => {
                assert_eq!(time, horizon);
                break;
            }
        }
    }

    assert!(sampler.log_likelihood_ratio().is_finite());
    let estimate = sampler.likelihood_ratio() * releases as f64 / horizon;
    if releases == 0 {
        assert_eq!(
            estimate, 0.0,
            "non-observation cannot become a positive rate"
        );
    }
}
