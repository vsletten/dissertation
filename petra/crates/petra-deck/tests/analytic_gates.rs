//! P2 analytic gates (design doc §8.2): run real KMC on systems with
//! closed-form answers and require statistical agreement. These are the
//! replacement for bitwise parity — the engine is validated against
//! textbook statistical mechanics, not a reference trajectory.
//!
//! Gate 1 — two-state equilibrium: independent sites flipping a↔b with
//! constant rates k_ab, k_ba must spend a time-fraction
//! k_ba/(k_ab + k_ba) in state a.
//!
//! Gate 2 — Langmuir isotherm: adsorption (with activity and Δμ folded in)
//! against desorption on independent sites must give coverage
//! θ = k_ads/(k_ads + k_des).

use petra_core::rate::R_KCAL;

const TWO_STATE: &str = r#"
[deck]
name = "two-state"
comment = "independent flippers, closed-form occupancy"

[cell]
a = 5.0
b = 5.0
c = 5.0
alpha = 90.0
beta = 90.0
gamma = 90.0

[[cell.sites]]
kind = "S"
frac = [0.0, 0.0, 0.0]

[[species]]
name = "A"

[[kinds]]
name = "S"
initial = "a"

[[kinds.states]]
name = "a"
occupant = "A"

[[kinds.states]]
name = "b"
occupant = "A"

[lattice]
dims = [4, 4, 4]
boundary = ["periodic", "periodic", "periodic"]

[thermo]
temperature = 300.0

[[reactions]]
name = "flip_ab"
center = { kind = "S", state = ["a"] }
rate = { constant = 1.0 }

[[reactions.effects]]
target = "center"
set = "b"

[[reactions]]
name = "flip_ba"
center = { kind = "S", state = ["b"] }
rate = { constant = 3.0 }

[[reactions.effects]]
target = "center"
set = "a"

[simulation]
steps = 100000
seed = 5
"#;

const LANGMUIR: &str = r#"
[deck]
name = "langmuir"
comment = "independent adsorption sites, closed-form isotherm"

[cell]
a = 5.0
b = 5.0
c = 5.0
alpha = 90.0
beta = 90.0
gamma = 90.0

[[cell.sites]]
kind = "S"
frac = [0.0, 0.0, 0.0]

[[species]]
name = "M"

[[kinds]]
name = "S"
initial = "vacant_site"

[[kinds.states]]
name = "vacant_site"
occupant = "vacant"

[[kinds.states]]
name = "adsorbed"
occupant = "M"

[lattice]
dims = [4, 4, 4]
boundary = ["periodic", "periodic", "periodic"]

[thermo]
temperature = 320.0

[thermo.mu]
M = -0.5

[thermo.activity]
M = 0.7

[[reactions]]
name = "adsorb"
center = { kind = "S", state = ["vacant_site"] }
rate = { arrhenius = { prefactor = 1.0e6, ea = 1.0 } }
consumes = ["M"]

[[reactions.effects]]
target = "center"
set = "adsorbed"

[[reactions]]
name = "desorb"
center = { kind = "S", state = ["adsorbed"] }
rate = { arrhenius = { prefactor = 1.0e6, ea = 1.5 } }
produces = ["M"]

[[reactions.effects]]
target = "center"
set = "vacant_site"

[simulation]
steps = 100000
seed = 9
"#;

/// Run `steps` KMC steps and return the time-weighted average fraction of
/// sites in state 0, discarding the first `burn_in` steps as equilibration.
fn time_averaged_state0_fraction(
    deck: &petra_deck::CompiledDeck,
    seed: u64,
    steps: u64,
    burn_in: u64,
) -> f64 {
    let mut engine = deck.build_engine(Some(seed));
    let n = engine.lattice.len() as f64;
    let mut last_time = 0.0;
    let mut weighted = 0.0;
    let mut window = 0.0;
    for i in 0..steps {
        // State *before* the step holds for the waiting time the step draws.
        let count0 = engine.state_counts(deck.n_states)[0] as f64;
        let fired = engine.step().expect("both directions always available");
        let dt = fired.time - last_time;
        last_time = fired.time;
        if i >= burn_in {
            weighted += count0 * dt;
            window += dt;
        }
    }
    weighted / (window * n)
}

#[test]
fn two_state_equilibrium_matches_closed_form() {
    let parsed: petra_deck::DeckFile = toml::from_str(TWO_STATE).expect("valid deck");
    let deck = petra_deck::compile(&parsed).expect("compiles");

    // k_ab = 1, k_ba = 3 → equilibrium fraction in a = 3/4.
    let frac_a = time_averaged_state0_fraction(&deck, 5, 100_000, 20_000);
    let expect = 3.0 / 4.0;
    assert!(
        (frac_a - expect).abs() < 0.01,
        "fraction in state a: got {frac_a:.4}, closed form {expect}"
    );
}

#[test]
fn langmuir_coverage_matches_closed_form() {
    let parsed: petra_deck::DeckFile = toml::from_str(LANGMUIR).expect("valid deck");
    let deck = petra_deck::compile(&parsed).expect("compiles");

    // θ = k_ads/(k_ads + k_des), with the solution factors in k_ads.
    let rt = R_KCAL * 320.0;
    let k_ads = 1.0e6 * (-1.0 / rt).exp() * 0.7 * (-0.5 / rt).exp();
    let k_des = 1.0e6 * (-1.5 / rt).exp();
    let theta_expect = k_ads / (k_ads + k_des);

    // State 0 is vacant_site; coverage is the complement.
    let frac_vacant = time_averaged_state0_fraction(&deck, 9, 100_000, 20_000);
    let theta = 1.0 - frac_vacant;
    assert!(
        (theta - theta_expect).abs() < 0.01,
        "coverage: got {theta:.4}, closed form {theta_expect:.4}"
    );
}

#[test]
fn mean_waiting_time_tracks_total_rate() {
    // For the two-state system the total rate is exactly
    // N_a·k_ab + N_b·k_ba at every instant; over many steps the mean of
    // dt·R_total must be 1 (unit-mean exponential draws).
    let parsed: petra_deck::DeckFile = toml::from_str(TWO_STATE).expect("valid deck");
    let deck = petra_deck::compile(&parsed).expect("compiles");
    let mut engine = deck.build_engine(Some(13));

    let mut last_time = 0.0;
    let mut sum = 0.0;
    let steps = 50_000;
    for _ in 0..steps {
        let counts = engine.state_counts(deck.n_states);
        let total_rate = counts[0] as f64 * 1.0 + counts[1] as f64 * 3.0;
        let fired = engine.step().expect("system never empties");
        sum += (fired.time - last_time) * total_rate;
        last_time = fired.time;
    }
    let mean = sum / steps as f64;
    assert!(
        (mean - 1.0).abs() < 0.02,
        "mean dt·R_total: got {mean:.4}, expected 1.0"
    );
}
