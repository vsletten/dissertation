//! A9 contract: the approximate-rate deck changes kinetics/solution context,
//! never the validated kaolinite topology, initialization, or reaction effects.

use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;

use petra_core::rate::RateExpr;

fn repo_path(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../..")
        .join(rel)
}

fn parsed(rel: &str) -> toml::Value {
    let text = std::fs::read_to_string(repo_path(rel)).expect("deck reads");
    toml::from_str(&text).expect("deck parses as TOML")
}

#[test]
fn approximate_deck_preserves_structure_init_and_reaction_semantics() {
    let legacy = parsed("petra/examples/kaolinite.toml");
    let approximate = parsed("petra/examples/kaolinite-approx.toml");

    for key in ["species", "kinds", "aliases", "cell", "lattice", "init"] {
        assert_eq!(approximate.get(key), legacy.get(key), "{key} drifted");
    }
    assert_eq!(
        approximate["simulation"], legacy["simulation"],
        "run length/cadence must remain the canonical 20k-step shape"
    );
    assert_eq!(approximate["observables"], legacy["observables"]);

    let legacy_reactions = legacy["reactions"].as_array().expect("legacy reactions");
    let approximate_reactions = approximate["reactions"]
        .as_array()
        .expect("approximate reactions");
    assert_eq!(approximate_reactions.len(), 22);
    assert_eq!(legacy_reactions.len(), approximate_reactions.len());
    for (old, new) in legacy_reactions.iter().zip(approximate_reactions) {
        let mut old = old.as_table().expect("legacy reaction table").clone();
        let mut new = new.as_table().expect("approximate reaction table").clone();
        old.remove("rate");
        new.remove("rate");
        assert_eq!(new, old, "non-rate reaction semantics drifted");
    }
}

#[test]
fn approximate_deck_has_exact_temperature_rate_contract_and_family_partition() {
    let path = repo_path("petra/examples/kaolinite-approx.toml");
    let source = std::fs::read_to_string(&path).expect("deck reads");
    let parsed: toml::Value = toml::from_str(&source).expect("deck parses");
    assert_eq!(parsed["deck"]["units"].as_str(), Some("kcal/mol"));
    assert_eq!(parsed["thermo"]["temperature"].as_float(), Some(298.0));
    assert!(
        parsed.get("execution").is_none(),
        "v1 deck must have no schedule"
    );

    let deck = petra_deck::load(path).expect("approximate deck compiles");
    assert_eq!(deck.temperature, 298.0);
    assert!(deck.schedule.is_empty());
    assert_eq!(deck.reactions.len(), 22);

    let names: BTreeSet<_> = deck
        .reactions
        .iter()
        .map(|reaction| reaction.name.as_str())
        .collect();
    assert_eq!(names.len(), 22, "reaction names must be unique");
    for reaction in &deck.reactions {
        match reaction.rate {
            RateExpr::Arrhenius { prefactor, ea } => {
                assert_eq!(prefactor, 1.0e13, "{} prefactor", reaction.name);
                assert!(ea.is_finite() && ea >= 0.0, "{} barrier", reaction.name);
            }
            ref other => panic!("{} is not Arrhenius: {other:?}", reaction.name),
        }
    }

    let mut families: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    let lines = source.lines().collect::<Vec<_>>();
    for window in lines.windows(3) {
        let Some(name) = window[0]
            .strip_prefix("name = \"")
            .and_then(|value| value.strip_suffix('\"'))
        else {
            continue;
        };
        let Some(family) = window[1]
            .strip_prefix("# a9-family = \"")
            .and_then(|value| value.strip_suffix('\"'))
        else {
            continue;
        };
        let provenance = window[2]
            .strip_prefix("# a9-provenance = \"")
            .and_then(|value| value.strip_suffix('\"'))
            .expect("family annotation must be adjacent to provenance");
        assert!(["computed", "literature", "heuristic"].contains(&provenance));
        families
            .entry(family.to_owned())
            .or_default()
            .insert(name.to_owned());
    }
    assert_eq!(families.values().map(BTreeSet::len).sum::<usize>(), 22);
    assert_eq!(families.len(), 7);
    assert_eq!(
        families.keys().map(String::as_str).collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "adsorption",
            "al-o-al-analogue",
            "cation-desorption",
            "connectivity-ladder",
            "siloxane-neutral",
            "sioal-al-neutral",
            "sioal-si-neutral",
        ])
    );
    assert!(source.contains("B3LYP/def2-SVP/DF"));
    assert!(source.contains("r2SCAN-3c Si-neutral electronic barrier is 30.053"));
    assert!(source.contains("pH 3-5 is contextual and is NOT"));
}
