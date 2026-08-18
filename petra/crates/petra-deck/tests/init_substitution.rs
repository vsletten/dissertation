//! Gates for substitution init rules (design doc §6, roadmap P5):
//! 1. probabilistic substitution hits the declared fraction (binomial bounds)
//!    and is deterministic given (deck, seed) — seed sweeps re-roll it;
//! 2. explicit `[a, b, c, t]` site lists substitute exactly those sites;
//! 3. composition with the other init filters (state sets, pass ordering);
//! 4. compile-time validation of probability range and site coordinates.

/// 40×40 single-site sheet: occupant A everywhere, an Fe-analogue state
/// `sub` on the same kind, and a 25% substitution pass. No reactions —
/// these gates only exercise the build.
const SUB_DECK: &str = r#"
[deck]
name = "substitution-check"

[cell]
a = 3.0
b = 3.0
c = 3.0
alpha = 90.0
beta = 90.0
gamma = 90.0

[[cell.sites]]
kind = "M"
frac = [0.0, 0.0, 0.0]

[[cell.bonds]]
i = 0
j = 0
dcell = [1, 0, 0]

[[cell.bonds]]
i = 0
j = 0
dcell = [0, 1, 0]

[[species]]
name = "A"

[[species]]
name = "F"

[[kinds]]
name = "M"
initial = "occupied"

[[kinds.states]]
name = "occupied"
occupant = "A"

[[kinds.states]]
name = "sub"
occupant = "F"

[[kinds.states]]
name = "empty"
occupant = "vacant"

[lattice]
dims = [40, 40, 1]
boundary = ["periodic", "periodic", "open"]

[thermo]
temperature = 300.0

[[init]]
name = "F substitution"
center = { kind = "M", state = ["occupied"] }
probability = 0.25
set = "sub"

[simulation]
steps = 1
seed = 7
"#;

fn parse(src: &str) -> petra_deck::DeckFile {
    toml::from_str(src).expect("deck parses")
}

fn compile(src: &str) -> petra_deck::CompiledDeck {
    petra_deck::compile(&parse(src)).expect("deck compiles")
}

/// Count of sites left in each state, keyed by qualified state name.
fn count(deck: &petra_deck::CompiledDeck, engine: &petra_core::Engine, state: &str) -> usize {
    engine
        .lattice
        .states
        .iter()
        .filter(|s| deck.state_names[s.0 as usize] == state)
        .count()
}

#[test]
fn substitution_fraction_lands_in_binomial_bounds() {
    let deck = compile(SUB_DECK);
    let engine = deck.build_engine(None).expect("engine builds");
    let n_sub = count(&deck, &engine, "M.sub");
    // N = 1600, p = 0.25: mean 400, σ ≈ 17.3; accept ±5σ.
    assert!(
        (313..=487).contains(&n_sub),
        "expected ≈400 of 1600 substituted, got {n_sub}"
    );
    assert_eq!(count(&deck, &engine, "M.occupied"), 1600 - n_sub);
}

#[test]
fn substitution_is_deterministic_given_seed_and_rerolls_across_seeds() {
    let deck = compile(SUB_DECK);
    let a = deck.build_engine(Some(11)).expect("engine builds");
    let b = deck.build_engine(Some(11)).expect("engine builds");
    let c = deck.build_engine(Some(12)).expect("engine builds");
    assert_eq!(
        a.lattice.states, b.lattice.states,
        "same seed must reproduce the same substitution pattern"
    );
    assert_ne!(
        a.lattice.states, c.lattice.states,
        "different seeds must re-roll the pattern"
    );
}

#[test]
fn explicit_site_list_substitutes_exactly_those_sites() {
    let mut deck = parse(SUB_DECK);
    let pass = &mut deck.init[0];
    pass.probability = None;
    pass.sites = Some(vec![[0, 0, 0, 0], [2, 1, 0, 0], [39, 39, 0, 0]]);
    let deck = petra_deck::compile(&deck).expect("deck compiles");
    let engine = deck.build_engine(None).expect("engine builds");
    // Flat index for this single-template sheet is a*40 + b.
    for (s, state) in engine.lattice.states.iter().enumerate() {
        let expect = matches!(s, 0 | 81 | 1599);
        assert_eq!(
            deck.state_names[state.0 as usize] == "M.sub",
            expect,
            "site {s}"
        );
    }
}

#[test]
fn substitution_respects_state_filter_and_pass_order() {
    // Clear a slab first; the p = 1 substitution that follows only sees
    // `occupied`, so the cleared slab must stay empty.
    let mut deck = parse(SUB_DECK);
    deck.init.insert(
        0,
        toml::from_str(
            r#"
name = "clear slab"
center = { kind = "M", state = ["occupied"] }
region = { axis = 0, max = 9 }
set = "empty"
"#,
        )
        .expect("pass parses"),
    );
    deck.init[1].probability = Some(1.0);
    let deck = petra_deck::compile(&deck).expect("deck compiles");
    let engine = deck.build_engine(None).expect("engine builds");
    for (s, state) in engine.lattice.states.iter().enumerate() {
        let name = &deck.state_names[state.0 as usize];
        let (cell, _) = engine.lattice.coords(s);
        let expect = if cell[0] <= 9 { "M.empty" } else { "M.sub" };
        assert_eq!(name, expect, "site {s} at a={}", cell[0]);
    }
}

#[test]
fn probability_and_site_list_validation() {
    let cases: &[(&str, &str)] = &[
        ("probability = 1.5", "probability must be in [0, 1]"),
        ("probability = -0.1", "probability must be in [0, 1]"),
        ("sites = []", "sites list must be non-empty"),
        ("sites = [[40, 0, 0, 0]]", "outside lattice dims"),
        ("sites = [[0, 0, 1, 0]]", "outside lattice dims"),
        ("sites = [[0, 0, 0, 1]]", "names template site"),
    ];
    for (line, expect) in cases {
        let src = SUB_DECK.replace("probability = 0.25", line);
        let e = petra_deck::compile(&parse(&src)).expect_err(line);
        assert!(
            e.to_string().contains(expect),
            "for `{line}`: got '{e}', wanted '{expect}'"
        );
    }
}
