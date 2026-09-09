use std::path::Path;

use petra_core::rate::RateExpr;
use petra_core::reaction::{EffectOp, EffectTarget, ModifierKind, NeighborSelect};
use petra_core::state::StateId;

fn state_names(deck: &petra_deck::CompiledDeck, select: &NeighborSelect) -> Vec<String> {
    deck.state_names
        .iter()
        .enumerate()
        .filter(|(index, _)| select.states.contains(StateId(*index as u16)))
        .map(|(_, name)| name.clone())
        .collect()
}

fn kind_name<'a>(deck: &'a petra_deck::CompiledDeck, select: &NeighborSelect) -> Option<&'a str> {
    select
        .kind
        .map(|kind| deck.kind_names[kind.0 as usize].as_str())
}

fn map_names(deck: &petra_deck::CompiledDeck, op: &EffectOp) -> (Vec<(String, String)>, bool) {
    match op {
        EffectOp::Map {
            entries,
            missing_is_error,
        } => (
            entries
                .iter()
                .map(|(left, right)| {
                    (
                        deck.state_names[left.0 as usize].clone(),
                        deck.state_names[right.0 as usize].clone(),
                    )
                })
                .collect(),
            *missing_is_error,
        ),
        other => panic!("expected compiled map, got {other:?}"),
    }
}

#[test]
fn calc005_live_si_pair_compiles_to_exact_runtime_semantics() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/kaolinite.toml");
    let deck = petra_deck::load(path).expect("kaolinite deck compiles");
    let adsorb = deck
        .reactions
        .iter()
        .find(|reaction| reaction.name == "adsorb-si")
        .expect("compiled adsorb-si");
    let desorb = deck
        .reactions
        .iter()
        .find(|reaction| reaction.name == "desorb-si")
        .expect("compiled desorb-si");

    assert_eq!(deck.kind_names[adsorb.center_kind.0 as usize], "Si");
    assert!(adsorb.center_states.contains(StateId(
        deck.state_names
            .iter()
            .position(|name| name == "Si.empty")
            .unwrap() as u16
    )));
    assert!(matches!(
        adsorb.rate,
        RateExpr::Arrhenius {
            prefactor: 100.0,
            ea: 6.4
        }
    ));
    assert!(adsorb.modifiers.is_empty());
    assert_eq!(adsorb.guards.len(), 1);
    let guard = &adsorb.guards[0];
    assert_eq!(guard.select.distance, 1);
    assert_eq!(kind_name(&deck, &guard.select), None);
    assert_eq!(guard.select.frozen, Some(false));
    assert_eq!(guard.min, 1);
    assert_eq!(guard.max, u32::MAX);
    assert_eq!(
        state_names(&deck, &guard.select),
        [
            "Oss.br",
            "Oss.hy",
            "Oss.si1",
            "Osa.br",
            "Osa.sih",
            "Osa.full",
            "Osa.albr",
            "Osa.alhy",
            "Osa.sial",
            "Osa.sialh",
            "Osa.si1",
            "Osa.al1",
            "Osa.pr1",
            "Osa.pr2",
            "Oaa.br",
            "Oaa.hy",
            "Oaa.al1",
        ]
    );

    assert_eq!(adsorb.branches.len(), 1);
    assert_eq!(adsorb.branches[0].effects.len(), 3);
    assert!(matches!(
        adsorb.branches[0].effects[0].target,
        EffectTarget::Center
    ));
    match adsorb.branches[0].effects[0].op {
        EffectOp::Set(state) => assert_eq!(deck.state_names[state.0 as usize], "Si.oh4"),
        ref other => panic!("expected center set, got {other:?}"),
    }
    let (oss_select, oss_op) = match &adsorb.branches[0].effects[1] {
        petra_core::reaction::Effect {
            target: EffectTarget::AllMatches(select),
            op,
        } => (select, op),
        other => panic!("expected Oss all-match map, got {other:?}"),
    };
    assert_eq!(oss_select.distance, 1);
    assert_eq!(kind_name(&deck, oss_select), Some("Oss"));
    assert_eq!(oss_select.frozen, Some(false));
    assert_eq!(state_names(&deck, oss_select), ["Oss.empty", "Oss.si1"]);
    assert_eq!(
        map_names(&deck, oss_op),
        (
            vec![
                ("Oss.empty".into(), "Oss.si1".into()),
                ("Oss.si1".into(), "Oss.hy".into()),
            ],
            true,
        )
    );
    let (osa_select, osa_op) = match &adsorb.branches[0].effects[2] {
        petra_core::reaction::Effect {
            target: EffectTarget::AllMatches(select),
            op,
        } => (select, op),
        other => panic!("expected Osa all-match map, got {other:?}"),
    };
    assert_eq!(osa_select.distance, 1);
    assert_eq!(kind_name(&deck, osa_select), Some("Osa"));
    assert_eq!(osa_select.frozen, Some(false));
    assert_eq!(
        state_names(&deck, osa_select),
        ["Osa.empty", "Osa.albr", "Osa.alhy", "Osa.al1"]
    );
    assert_eq!(
        map_names(&deck, osa_op),
        (
            vec![
                ("Osa.al1".into(), "Osa.sialh".into()),
                ("Osa.albr".into(), "Osa.sih".into()),
                ("Osa.alhy".into(), "Osa.full".into()),
                ("Osa.empty".into(), "Osa.si1".into()),
            ],
            true,
        )
    );

    assert_eq!(deck.kind_names[desorb.center_kind.0 as usize], "Si");
    assert!(desorb.center_states.contains(StateId(
        deck.state_names
            .iter()
            .position(|name| name == "Si.oh4")
            .unwrap() as u16
    )));
    assert!(matches!(
        desorb.rate,
        RateExpr::Arrhenius {
            prefactor: 1998.0,
            ea: 6.0
        }
    ));
    assert_eq!(desorb.modifiers.len(), 2);
    assert_eq!(desorb.modifiers[0].select.distance, 1);
    assert_eq!(desorb.modifiers[0].select.frozen, Some(false));
    assert_eq!(state_names(&deck, &desorb.modifiers[0].select), ["Oss.hy"]);
    assert!(matches!(
        &desorb.modifiers[0].kind,
        ModifierKind::ByCount { dea }
            if dea == &[0.0, 6.0, 12.0, 18.0, 24.0]
    ));
    assert_eq!(desorb.modifiers[1].select.distance, 1);
    assert_eq!(desorb.modifiers[1].select.frozen, Some(false));
    assert_eq!(state_names(&deck, &desorb.modifiers[1].select), ["Osa.si1"]);
    assert!(matches!(
        desorb.modifiers[1].kind,
        ModifierKind::When {
            min: 1,
            max: u32::MAX,
            dea: -6.0,
            factor: 1.0,
        }
    ));
    assert_eq!(desorb.branches.len(), 1);
    assert_eq!(desorb.branches[0].effects.len(), 3);
    assert!(matches!(
        desorb.branches[0].effects[0].target,
        EffectTarget::Center
    ));
    match desorb.branches[0].effects[0].op {
        EffectOp::Set(state) => assert_eq!(deck.state_names[state.0 as usize], "Si.empty"),
        ref other => panic!("expected center empty set, got {other:?}"),
    }
    let (desorb_oss_select, desorb_oss_op) = match &desorb.branches[0].effects[1] {
        petra_core::reaction::Effect {
            target: EffectTarget::AllMatches(select),
            op,
        } => (select, op),
        other => panic!("expected desorb Oss all-match map, got {other:?}"),
    };
    assert_eq!(desorb_oss_select.distance, 1);
    assert_eq!(kind_name(&deck, desorb_oss_select), Some("Oss"));
    assert_eq!(desorb_oss_select.frozen, Some(false));
    assert_eq!(
        state_names(&deck, desorb_oss_select),
        ["Oss.empty", "Oss.br", "Oss.hy", "Oss.si1"]
    );
    assert_eq!(
        map_names(&deck, desorb_oss_op),
        (
            vec![
                ("Oss.hy".into(), "Oss.si1".into()),
                ("Oss.si1".into(), "Oss.empty".into()),
            ],
            false,
        )
    );
    let (desorb_osa_select, desorb_osa_op) = match &desorb.branches[0].effects[2] {
        petra_core::reaction::Effect {
            target: EffectTarget::AllMatches(select),
            op,
        } => (select, op),
        other => panic!("expected desorb Osa all-match map, got {other:?}"),
    };
    assert_eq!(desorb_osa_select.distance, 1);
    assert_eq!(kind_name(&deck, desorb_osa_select), Some("Osa"));
    assert_eq!(desorb_osa_select.frozen, Some(false));
    assert_eq!(state_names(&deck, desorb_osa_select).len(), 12);
    assert_eq!(
        map_names(&deck, desorb_osa_op),
        (
            vec![
                ("Osa.full".into(), "Osa.alhy".into()),
                ("Osa.si1".into(), "Osa.empty".into()),
                ("Osa.sialh".into(), "Osa.al1".into()),
                ("Osa.sih".into(), "Osa.albr".into()),
            ],
            false,
        )
    );
}
