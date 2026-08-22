use petra_deck::{ScheduleAdvance, ScheduleSegment};

type EventKey = (u64, u64, usize, u16);
type SegmentEvents = Vec<Vec<EventKey>>;

fn deck_text(ea: f64, segments: &[(f64, f64)]) -> String {
    let schedule = segments
        .iter()
        .map(|(temperature, duration)| {
            format!("[[execution.schedule]]\ntemperature = {temperature}\nduration = {duration}\n")
        })
        .collect::<Vec<_>>()
        .join("\n");
    format!(
        r#"
[deck]
name = "scheduled-decay"
schema = 2
units = "kcal/mol"

[structure]
kind = "grid"
grid = "square"
neighborhood = "von_neumann"
dims = [8, 8]
boundary = ["periodic", "periodic"]
default_kind = "S"

[[structure.species]]
name = "A"

[[structure.kinds]]
name = "S"
initial = "a"
[[structure.kinds.states]]
name = "a"
occupant = "A"
[[structure.kinds.states]]
name = "b"
occupant = "vacant"

[dynamics.thermo]
temperature = {initial_temperature}

[[dynamics.rules]]
name = "decay"
center = {{ kind = "S", state = ["a"] }}
rate = {{ arrhenius = {{ prefactor = 100.0, ea = {ea} }} }}
[[dynamics.rules.effects]]
target = "center"
set = "b"

[execution]
strategy = "ctmc"
{schedule}
[execution.stop]
steps = 1000

[execution.ensemble]
seed = 17
n_replicas = 1
seed_policy = "increment"

[observables]
report_every = 10
"#,
        initial_temperature = segments[0].0,
    )
}

fn compile_schedule(ea: f64, segments: &[(f64, f64)]) -> petra_deck::CompiledDeck {
    let parsed: petra_deck::DeckFile =
        toml::from_str(&deck_text(ea, segments)).expect("schedule parses");
    petra_deck::compile(&parsed).expect("schedule compiles")
}

fn run_schedule(deck: &petra_deck::CompiledDeck, seed: u64) -> (SegmentEvents, Vec<u64>) {
    let mut run = deck.build_schedule(Some(seed)).expect("scheduled engine");
    let mut by_segment = vec![Vec::new(); deck.schedule.len()];
    for _ in 0..10_000 {
        match run.advance().expect("schedule advances") {
            ScheduleAdvance::Fired { segment, fired } => by_segment[segment].push((
                fired.step,
                fired.time.to_bits(),
                fired.site,
                fired.reaction,
            )),
            ScheduleAdvance::Boundary { .. } => {}
            ScheduleAdvance::Complete => {
                let counts = run.engine().state_counts(deck.n_states);
                return (by_segment, counts);
            }
        }
    }
    panic!("bounded schedule did not complete");
}

#[test]
fn schedule_schema_preserves_order_and_rejects_invalid_segments() {
    let parsed: petra_deck::DeckFile =
        toml::from_str(&deck_text(1.0, &[(300.0, 0.4), (600.0, 0.2)]))
            .expect("valid schedule parses");
    assert_eq!(
        parsed.execution.schedule,
        vec![
            ScheduleSegment {
                temperature: 300.0,
                duration: 0.4,
            },
            ScheduleSegment {
                temperature: 600.0,
                duration: 0.2,
            },
        ]
    );

    let zero_duration = deck_text(1.0, &[(300.0, 0.0)]);
    let error = toml::from_str::<petra_deck::DeckFile>(&zero_duration)
        .expect_err("zero-duration segment rejected");
    assert!(error.to_string().contains("duration"), "{error}");

    let non_ctmc =
        deck_text(1.0, &[(300.0, 0.4)]).replace("strategy = \"ctmc\"", "strategy = \"pca\"");
    let error = toml::from_str::<petra_deck::DeckFile>(&non_ctmc)
        .expect_err("schedule on a discrete strategy rejected");
    assert!(error.to_string().contains("CTMC"), "{error}");

    let overflow = deck_text(1.0, &[(300.0, 1.0), (600.0, 1.0)])
        .replace("duration = 1\n", "duration = 1.0e308\n");
    let error = toml::from_str::<petra_deck::DeckFile>(&overflow)
        .expect_err("cumulative duration overflow rejected");
    assert!(error.to_string().contains("cumulative"), "{error}");

    let rounded_away = deck_text(1.0, &[(300.0, 1.0), (600.0, 1.0)]).replacen(
        "duration = 1\n",
        "duration = 1.0e300\n",
        1,
    );
    let error = toml::from_str::<petra_deck::DeckFile>(&rounded_away)
        .expect_err("a segment that cannot advance wall time is rejected");
    assert!(error.to_string().contains("advance"), "{error}");
}

#[test]
fn schedule_preserves_an_explicit_zero_event_cap() {
    let text = deck_text(1.0, &[(300.0, 0.4)]).replace("steps = 1000", "steps = 0");
    let parsed: petra_deck::DeckFile = toml::from_str(&text).expect("zero cap parses");
    let deck = petra_deck::compile(&parsed).expect("zero cap compiles");
    assert_eq!(deck.step_limit, Some(0));
}

#[test]
fn first_schedule_temperature_recompiles_chemical_potential_coupling() {
    let text = deck_text(1.0, &[(300.0, 0.4)])
        .replacen(
            "[dynamics.thermo]\ntemperature = 300",
            "[dynamics.thermo]\ntemperature = 100\n[dynamics.thermo.mu]\nA = -10.0",
            1,
        )
        .replace(
            "rate = { arrhenius = { prefactor = 100.0, ea = 1 } }",
            "rate = { arrhenius = { prefactor = 100.0, ea = 1 } }\nconsumes = [\"A\"]",
        );
    let parsed: petra_deck::DeckFile = toml::from_str(&text).expect("thermo schedule parses");
    let deck = petra_deck::compile(&parsed).expect("thermo schedule compiles");
    let run = deck.build_schedule(Some(4)).expect("schedule builds");
    let expected = -10.0 / (petra_core::rate::R_KCAL * 300.0);
    assert_eq!(run.engine().temperature().to_bits(), 300.0f64.to_bits());
    assert!((run.engine().reactions[0].ln_thermo - expected).abs() < 1.0e-12);
}

#[test]
fn first_segment_is_bitwise_identical_to_an_independent_isothermal_run() {
    let scheduled = compile_schedule(1.0, &[(300.0, 0.4), (600.0, 0.2)]);
    let isothermal_text = deck_text(1.0, &[(300.0, 0.4)]).replace(
        "[[execution.schedule]]\ntemperature = 300\nduration = 0.4\n",
        "",
    );
    let parsed: petra_deck::DeckFile =
        toml::from_str(&isothermal_text).expect("independent isothermal deck parses");
    let isothermal = petra_deck::compile(&parsed).expect("independent isothermal deck compiles");

    let (scheduled_events, _) = run_schedule(&scheduled, 91);
    let mut engine = isothermal
        .build_engine(Some(91))
        .expect("independent isothermal engine");
    let mut isothermal_events = Vec::new();
    let boundary_counts = loop {
        let before = engine.state_counts(isothermal.n_states);
        match engine.step() {
            Ok(fired) if fired.time >= 0.4 => break before,
            Ok(fired) => isothermal_events.push((
                fired.step,
                fired.time.to_bits(),
                fired.site,
                fired.reaction,
            )),
            Err(petra_core::Stop::NoEvents) => break before,
            Err(error) => panic!("unexpected isothermal stop: {error}"),
        }
    };

    assert!(
        !scheduled_events[0].is_empty(),
        "gate must exercise real events"
    );
    assert_eq!(scheduled_events[0], isothermal_events);
    let mut first_segment = scheduled
        .build_schedule(Some(91))
        .expect("scheduled engine");
    loop {
        match first_segment.advance().expect("first segment advances") {
            ScheduleAdvance::Boundary {
                completed_segment: 0,
                ..
            } => break,
            ScheduleAdvance::Fired { .. } => {}
            other => panic!("unexpected first-segment outcome: {other:?}"),
        }
    }
    assert_eq!(
        first_segment.engine().state_counts(scheduled.n_states),
        boundary_counts,
        "state at the wall-time boundary is exact-replay identical"
    );
}

#[test]
fn second_segment_is_seeded_ensemble_consistent_with_independent_t2_runs() {
    // At 1 K this barrier underflows to a zero rate, so every scheduled run
    // reaches T2 from exactly the declared initial state. The boundary still
    // exercises the zero-rate-to-live-rate table rebuild path.
    let scheduled = compile_schedule(1_000.0, &[(1.0, 1.0), (100_000.0, 1.0)]);
    let independent_t2 = compile_schedule(1_000.0, &[(100_000.0, 1.0)]);
    let mut scheduled_b = Vec::new();
    let mut independent_b = Vec::new();
    for seed in 0..128 {
        let (events, counts) = run_schedule(&scheduled, seed);
        assert!(events[0].is_empty(), "T1 state must be unchanged");
        scheduled_b.push(counts[1] as f64);
        independent_b.push(run_schedule(&independent_t2, seed + 10_000).1[1] as f64);
    }
    let mean = |values: &[f64]| values.iter().sum::<f64>() / values.len() as f64;
    let scheduled_mean = mean(&scheduled_b);
    let independent_mean = mean(&independent_b);
    assert!(
        scheduled_mean > 20.0,
        "T2 rates were not rebuilt: {scheduled_mean}"
    );
    assert!(
        (scheduled_mean - independent_mean).abs() < 2.0,
        "scheduled T2 mean {scheduled_mean} vs independent T2 mean {independent_mean}"
    );
}
