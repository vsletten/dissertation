use petra_core::crystal::{Cell, KindId, TemplateSite, UnitCell};
use petra_core::rate::RateExpr;
use petra_core::reaction::{
    Branch, Effect, EffectOp, EffectTarget, NeighborSelect, Reaction, RuleValue,
};
use petra_core::{BiasedCtmc, Boundary, CtmcAdvance, Engine, Lattice, StateId, StateSet};

fn one_shot_engine(seed: u64) -> Engine {
    let cell = UnitCell {
        cell: Cell::from_params(1.0, 1.0, 1.0, 90.0, 90.0, 90.0),
        sites: vec![TemplateSite {
            kind: KindId(0),
            frac: [0.0; 3],
            bonds: Vec::new(),
        }],
    };
    let lattice = Lattice::build(&cell, [1, 1, 1], [Boundary::Periodic; 3], |_| StateId(0));
    let mut center_states = StateSet::new(2);
    center_states.insert(StateId(0));
    let reaction = Reaction {
        name: "release".into(),
        center_kind: KindId(0),
        center_states,
        guards: Vec::new(),
        rate: RateExpr::Constant { k: 2.0 },
        value: RuleValue::Ctmc,
        ln_activity: 0.0,
        mu_energy: 0.0,
        ln_thermo: 0.0,
        strain_scale: 0.0,
        modifiers: Vec::new(),
        branches: vec![Branch {
            weight: 1.0,
            effects: vec![Effect {
                target: EffectTarget::Center,
                op: EffectOp::Set(StateId(1)),
            }],
        }],
    };
    Engine::new(
        lattice,
        &[KindId(0)],
        1,
        vec![(0, 2)],
        vec![reaction],
        298.0,
        seed,
    )
}

#[test]
fn unit_bias_is_bitwise_identical_to_exact_ctmc() {
    let mut exact = one_shot_engine(17);
    let mut biased = one_shot_engine(17);
    let mut sampler = BiasedCtmc::new(1, []).expect("unit bias");

    let exact_event = exact.advance_ctmc_until(1.0e9).expect("exact event");
    let biased_event = biased
        .advance_biased_ctmc_until(&mut sampler, 1.0e9)
        .expect("biased event");

    assert_eq!(biased_event, exact_event);
    assert_eq!(biased.lattice.states, exact.lattice.states);
    assert_eq!(sampler.log_likelihood_ratio(), 0.0);
}

#[test]
fn biased_event_accumulates_exact_path_likelihood_ratio() {
    let mut engine = one_shot_engine(31);
    let mut sampler = BiasedCtmc::new(1, [(0, 4.0)]).expect("valid bias");

    let CtmcAdvance::Fired(event) = engine
        .advance_biased_ctmc_until(&mut sampler, 1.0e9)
        .expect("release fires")
    else {
        panic!("deadline unexpectedly reached");
    };
    let expected = (8.0 - 2.0) * event.time - 4.0_f64.ln();
    assert!((sampler.log_likelihood_ratio() - expected).abs() < 1.0e-14);
}

#[test]
fn biased_deadline_accumulates_survival_likelihood_without_an_event() {
    let mut engine = one_shot_engine(41);
    let mut sampler = BiasedCtmc::new(1, [(0, 4.0)]).expect("valid bias");
    let deadline = 1.0e-12;

    assert_eq!(
        engine
            .advance_biased_ctmc_until(&mut sampler, deadline)
            .expect("deadline"),
        CtmcAdvance::Deadline { time: deadline }
    );
    assert_eq!(engine.step_count, 0);
    assert_eq!(engine.lattice.states, vec![StateId(0)]);
    assert!((sampler.log_likelihood_ratio() - (8.0 - 2.0) * deadline).abs() < 1.0e-20);
}

#[test]
fn invalid_bias_contract_fails_closed() {
    assert!(BiasedCtmc::new(1, [(1, 2.0)]).is_err());
    assert!(BiasedCtmc::new(1, [(0, 0.0)]).is_err());
    assert!(BiasedCtmc::new(1, [(0, f64::NAN)]).is_err());
}

#[test]
fn failed_transition_is_atomic_across_weight_time_state_and_rng() {
    let mut engine = one_shot_engine(53);
    let mut clean = one_shot_engine(53);
    let original_target = engine.reactions[0].branches[0].effects[0].target.clone();
    let mut impossible_states = StateSet::new(2);
    impossible_states.insert(StateId(0));
    engine.reactions[0].branches[0].effects[0].target = EffectTarget::FirstMatch(NeighborSelect {
        distance: 1,
        kind: None,
        label: None,
        exclude_label: None,
        frozen: None,
        states: impossible_states,
    });
    let mut sampler = BiasedCtmc::new(1, [(0, 4.0)]).expect("valid bias");

    assert!(engine
        .advance_biased_ctmc_until(&mut sampler, 1.0e9)
        .is_err());
    assert_eq!(sampler.log_likelihood_ratio(), 0.0);
    assert_eq!(engine.time, 0.0);
    assert_eq!(engine.step_count, 0);
    assert_eq!(engine.lattice.states, vec![StateId(0)]);

    engine.reactions[0].branches[0].effects[0].target = original_target;
    let retried = engine
        .advance_biased_ctmc_until(&mut sampler, 1.0e9)
        .expect("retry succeeds");
    let mut clean_sampler = BiasedCtmc::new(1, [(0, 4.0)]).expect("valid bias");
    let baseline = clean
        .advance_biased_ctmc_until(&mut clean_sampler, 1.0e9)
        .expect("clean run succeeds");

    assert_eq!(retried, baseline, "failed transition must not consume RNG");
    assert_eq!(sampler, clean_sampler);
}
