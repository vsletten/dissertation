//! Gates for defect strain fields (docs/STRAIN.md §5):
//! 1. the analytic field — values, core clamp, minimum-image wrap;
//! 2. exactness of the rate coupling, k(u)/k(0) = exp(−scale·u/RT);
//! 3. the physics: an etch pit opens preferentially around the screw core
//!    in the Kossel demo deck (ensemble check).

use std::path::PathBuf;

use petra_core::rate::R_KCAL;
use petra_core::reaction::resolve_rate;

/// 9×9 single-site sheet, screw along c with a bare prefactor: u = 16/r²,
/// clamp at 2 Å. Cell is a 3 Å cube so cartesian = 3 × cell coords.
const FIELD_DECK: &str = r#"
[deck]
name = "field-check"

[cell]
a = 3.0
b = 3.0
c = 3.0
alpha = 90.0
beta = 90.0
gamma = 90.0

[[cell.sites]]
kind = "X"
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

[[kinds]]
name = "X"
initial = "occupied"

[[kinds.states]]
name = "occupied"
occupant = "A"

[[kinds.states]]
name = "empty"
occupant = "vacant"

[lattice]
dims = [9, 9, 1]
boundary = ["periodic", "periodic", "open"]

[[defects]]
type = "screw"
line_axis = 2
at = [1.0, 1.0, 0.0]
strain_prefactor = 16.0
core_radius = 2.0

[thermo]
temperature = 300.0

[[reactions]]
name = "detach"
center = { kind = "X", state = ["occupied"] }
rate = { arrhenius = { prefactor = 1.0e10, ea = 3.0 } }
strain = { scale = -1.0 }

[[reactions.effects]]
target = "center"
set = "empty"

[simulation]
steps = 10
seed = 1
"#;

fn compile(src: &str) -> petra_deck::CompiledDeck {
    let deck: petra_deck::DeckFile = toml::from_str(src).expect("deck parses");
    petra_deck::compile(&deck).expect("deck compiles")
}

/// Row-major site index for the 9×9 single-site sheet: index = a*9 + b.
const fn site(a: usize, b: usize) -> usize {
    a * 9 + b
}

#[test]
fn strain_field_matches_hand_formula_with_clamp_and_min_image() {
    let deck = compile(FIELD_DECK);
    let engine = deck.build_engine(Some(1)).expect("engine builds");
    let lat = &engine.lattice;
    let u = |a: usize, b: usize| lat.strain[site(a, b)];

    // Line pierces cell (1,1) → cartesian (3,3). Distances in Å:
    // site (1,1): r = 0 → clamped to 2 → u = 16/4 = 4.
    assert!((u(1, 1) - 4.0).abs() < 1e-12, "core clamp: {}", u(1, 1));
    // site (2,1): r = 3 → u = 16/9.
    assert!((u(2, 1) - 16.0 / 9.0).abs() < 1e-12);
    // site (3,5): dr = (6, 12) Å → r² = 180 → u = 16/180.
    assert!((u(3, 5) - 16.0 / 180.0).abs() < 1e-12);
    // Minimum image: site (8,1) is 7 cells = 21 Å away directly, but only
    // 2 cells = 6 Å through the periodic wrap → u = 16/36.
    assert!((u(8, 1) - 16.0 / 36.0).abs() < 1e-12, "min-image: {}", u(8, 1));

    // Superposition sanity: exactly one defect, all sites strained > 0.
    assert!(lat.strain.iter().all(|&v| v > 0.0));
}

#[test]
fn strain_cap_bounds_the_field() {
    // Large prefactor, small cap: without the cap the clamped-core value
    // would be 1e6/4; with it, no site may exceed 0.5.
    let mut deck: petra_deck::DeckFile = toml::from_str(FIELD_DECK).expect("deck parses");
    deck.defects[0].strain_prefactor = Some(1.0e6);
    deck.defects[0].cap = Some(0.5);
    let deck = petra_deck::compile(&deck).expect("deck compiles");
    let engine = deck.build_engine(Some(1)).expect("engine builds");
    let max_u = engine.lattice.strain.iter().cloned().fold(0.0f64, f64::max);
    assert!(
        max_u <= 0.5 + 1e-12 && max_u > 0.0,
        "cap not enforced: max u = {max_u}"
    );
}

#[test]
fn no_defects_means_zero_strain_everywhere() {
    let mut deck: petra_deck::DeckFile = toml::from_str(FIELD_DECK).expect("deck parses");
    deck.defects.clear();
    let deck = petra_deck::compile(&deck).expect("deck compiles");
    let engine = deck.build_engine(Some(1)).expect("engine builds");
    assert!(engine.lattice.strain.iter().all(|&u| u == 0.0));
}

#[test]
fn strain_rate_coupling_is_exact() {
    let deck = compile(FIELD_DECK);
    let engine = deck.build_engine(Some(1)).expect("engine builds");
    let lat = &engine.lattice;
    let kinds: Vec<_> = lat
        .template_index
        .iter()
        .map(|&t| deck.kinds_per_template[t as usize])
        .collect();
    let rxn = &deck.reactions[0];
    let rt = R_KCAL * 300.0;
    let mut scratch = Vec::new();

    // Both sites are interior with identical 4-neighbor environments; the
    // only rate difference is the strain term: k1/k2 = exp((u1 − u2)/RT)
    // for scale = −1.
    let s1 = site(1, 1); // at the core (u = 4)
    let s2 = site(5, 5); // far away
    let k1 = resolve_rate(lat, &kinds, rxn, s1, 300.0, &mut scratch);
    let k2 = resolve_rate(lat, &kinds, rxn, s2, 300.0, &mut scratch);
    let expect = ((lat.strain[s1] - lat.strain[s2]) / rt).exp();
    assert!(
        ((k1 / k2) / expect - 1.0).abs() < 1e-12,
        "rate ratio {} vs exp(Δu/RT) {expect}",
        k1 / k2
    );
}

#[test]
fn screw_dislocation_opens_an_etch_pit() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../examples/kossel-etchpit.toml");
    let deck = petra_deck::load(path).expect("etch-pit deck compiles");

    // In-plane distance (Å) from a site to the dislocation line at cell
    // (7,7): orthogonal 3 Å cell, periodic min-image in a and b.
    let dist = |a: usize, b: usize| -> f64 {
        let d = |x: usize| {
            let raw = (x as f64 - 7.0).abs();
            raw.min(14.0 - raw) * 3.0
        };
        (d(a).powi(2) + d(b).powi(2)).sqrt()
    };

    // Ensemble of 3 seeds: fraction dissolved near the core (r < 7 Å) vs
    // far field (r > 15 Å) after 1200 events (the far field is largely intact then; by ~4000 the whole slab is gone at this undersaturation).
    let mut near_frac = 0.0;
    let mut far_frac = 0.0;
    for seed in [3u64, 4, 5] {
        let mut engine = deck.build_engine(Some(seed)).expect("engine builds");
        for _ in 0..1200 {
            if engine.step().is_err() {
                break;
            }
        }
        let lat = &engine.lattice;
        let (mut near, mut near_empty, mut far, mut far_empty) = (0u32, 0u32, 0u32, 0u32);
        for s in 0..lat.len() {
            let (cell, _) = lat.coords(s);
            let r = dist(cell[0], cell[1]);
            let empty = deck.state_names[lat.states[s].0 as usize].ends_with(".empty");
            if r < 7.0 {
                near += 1;
                near_empty += u32::from(empty);
            } else if r > 15.0 {
                far += 1;
                far_empty += u32::from(empty);
            }
        }
        near_frac += near_empty as f64 / near as f64;
        far_frac += far_empty as f64 / far as f64;
    }
    near_frac /= 3.0;
    far_frac /= 3.0;
    assert!(
        near_frac > far_frac + 0.10,
        "expected preferential dissolution at the core: near {near_frac:.3} vs far {far_frac:.3}"
    );
}
