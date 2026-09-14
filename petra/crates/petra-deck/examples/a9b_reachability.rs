//! Deterministic A9b mechanism-reachability proof for the production deck.
//!
//! Usage: cargo run -p petra-deck --example a9b_reachability -- DECK SNAPSHOT OUT

use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::env;
use std::fs;
use std::path::{Path, PathBuf};

use petra_core::{Fired, StepCtx, StepOutcome, Stop, UpdateStrategy};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

const FORWARD_REACTIONS: &[&str] = &[
    "R0-sio-si-hydrolysis",
    "R2-sioal2-si-hydrolysis",
    "R4-sioal2-al-hydrolysis",
    "R6-second-stage-hydrolysis",
    "R8a-sioha-hydrolysis",
    "R8b-sioha-hydrolysis",
    "R10-sial-hydrolysis",
    "R12-albr-hydrolysis",
    "R14-alohal-hydrolysis",
];
const MAX_EVENTS: usize = 2_000;

struct ForwardClosure<'a> {
    releases: [u16; 2],
    forward: &'a [u16],
    original: &'a HashSet<usize>,
}

impl UpdateStrategy for ForwardClosure<'_> {
    fn step(&mut self, ctx: &mut StepCtx<'_>) -> Result<StepOutcome, Stop> {
        // Release as soon as a lattice target is eligible. This prevents an
        // already-proved target from affecting later FirstMatch selection.
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

fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn reaction_id(deck: &petra_deck::CompiledDeck, name: &str) -> Result<u16, String> {
    deck.reactions
        .iter()
        .position(|reaction| reaction.name == name)
        .map(|index| index as u16)
        .ok_or_else(|| format!("required reaction is missing: {name}"))
}

fn categorical(column: &Value) -> Result<Vec<String>, String> {
    let data = column["data"]
        .as_array()
        .ok_or("categorical column has no data array")?;
    let dictionary = column["dict"]
        .as_array()
        .ok_or("categorical column has no dictionary")?;
    data.iter()
        .map(|index| {
            let index = index
                .as_u64()
                .ok_or("categorical index is not an integer")? as usize;
            dictionary
                .get(index)
                .and_then(Value::as_str)
                .map(str::to_owned)
                .ok_or_else(|| format!("categorical index {index} is outside its dictionary"))
        })
        .collect()
}

fn validate_snapshot(
    deck: &petra_deck::CompiledDeck,
    engine: &petra_core::Engine,
    bytes: &[u8],
) -> Result<(), String> {
    let snapshot: Value = serde_json::from_slice(bytes).map_err(|error| error.to_string())?;
    if snapshot["pgif"] != 1 || snapshot["meta"]["directed"] != false {
        return Err("snapshot must be undirected PGIF v1".into());
    }
    let nodes = &snapshot["nodes"];
    let count = nodes["count"]
        .as_u64()
        .ok_or("snapshot node count is missing")? as usize;
    if count != engine.lattice.len() {
        return Err(format!(
            "snapshot has {count} nodes; live deck initializes {}",
            engine.lattice.len()
        ));
    }
    let states = categorical(&nodes["columns"]["state"])?;
    let kinds = categorical(&nodes["columns"]["kind"])?;
    let frozen = nodes["columns"]["frozen"]["data"]
        .as_array()
        .ok_or("snapshot frozen column has no data array")?;
    for site in 0..count {
        let live_state = &deck.state_names[engine.lattice.states[site].0 as usize];
        let template = engine.lattice.template_index[site] as usize;
        let live_kind = &deck.kind_names[deck.kinds_per_template[template].0 as usize];
        if states.get(site) != Some(live_state) || kinds.get(site) != Some(live_kind) {
            return Err(format!("snapshot/live deck state mismatch at site {site}"));
        }
        if frozen.get(site).and_then(Value::as_bool) != Some(engine.lattice.frozen[site]) {
            return Err(format!("snapshot/live deck frozen mismatch at site {site}"));
        }
    }
    Ok(())
}

fn ladder(kind: &str) -> &'static [&'static str] {
    match kind {
        "Al" => &[
            "Al.l0", "Al.l1", "Al.l2", "Al.l3", "Al.l4", "Al.l5", "Al.l6",
        ],
        "Si" => &["Si.oh0", "Si.oh1", "Si.oh2", "Si.oh3", "Si.oh4"],
        _ => unreachable!("only lattice cations are targets"),
    }
}

fn compatible_bridge(kind: &str, oxygen_kind: &str, state: &str) -> bool {
    state.ends_with(".br")
        && match kind {
            "Al" => matches!(oxygen_kind, "Oaa" | "Osa"),
            "Si" => matches!(oxygen_kind, "Oss" | "Osa"),
            _ => false,
        }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args_os().skip(1).collect();
    if args.len() != 3 {
        return Err("usage: a9b_reachability DECK SNAPSHOT OUT".into());
    }
    let deck_path = PathBuf::from(&args[0]);
    let snapshot_path = PathBuf::from(&args[1]);
    let out_path = PathBuf::from(&args[2]);
    let deck_bytes = fs::read(&deck_path)?;
    let snapshot_bytes = fs::read(&snapshot_path)?;
    let deck = petra_deck::load(&deck_path)?;
    if deck.temperature != 298.0 {
        return Err(format!("A9b requires exact 298 K deck, got {}", deck.temperature).into());
    }
    let mut engine = deck.build_engine(None)?;
    validate_snapshot(&deck, &engine, &snapshot_bytes)?;

    let initial_states = engine.lattice.states.clone();
    let mut kinds = Vec::with_capacity(engine.lattice.len());
    for &template in &engine.lattice.template_index {
        let kind = deck.kinds_per_template[template as usize];
        kinds.push(deck.kind_names[kind.0 as usize].clone());
    }
    let original: HashSet<usize> = initial_states
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

    let forward: Vec<u16> = FORWARD_REACTIONS
        .iter()
        .map(|name| reaction_id(&deck, name))
        .collect::<Result<_, _>>()?;
    let releases = [
        reaction_id(&deck, "desorb-al")?,
        reaction_id(&deck, "desorb-si")?,
    ];
    let mut release_event = BTreeMap::new();
    let mut events = Vec::new();

    for event_index in 0..MAX_EVENTS {
        let mut closure = ForwardClosure {
            releases,
            forward: &forward,
            original: &original,
        };
        match engine.step_with(&mut closure) {
            Ok(outcome) => {
                let fired = outcome.fired[0];
                let name = deck.reactions[fired.reaction as usize].name.clone();
                let changes: Vec<_> = engine
                    .last_changes()
                    .iter()
                    .map(|&(site, old, new)| {
                        json!({
                            "site": site,
                            "from": deck.state_names[old.0 as usize],
                            "to": deck.state_names[new.0 as usize],
                        })
                    })
                    .collect();
                if releases.contains(&fired.reaction) {
                    release_event.insert(fired.site, event_index);
                }
                events.push(json!({
                    "index": event_index,
                    "reaction": name,
                    "center_site": fired.site,
                    "changes": changes,
                }));
            }
            Err(Stop::NoEvents) => break,
            Err(error) => return Err(format!("closure failed: {error}").into()),
        }
    }
    if events.len() == MAX_EVENTS {
        return Err(format!("forward closure exceeded bounded event cap {MAX_EVENTS}").into());
    }

    let mut target_rows = Vec::new();
    let mut unreachable = BTreeSet::new();
    let mut by_kind: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for &site in &original {
        let kind = &kinds[site];
        by_kind.entry(kind.clone()).or_default().push(site);
        let initial = &deck.state_names[initial_states[site].0 as usize];
        let states = ladder(kind);
        let initial_index = states.iter().position(|state| *state == initial).unwrap();
        let release = release_event.get(&site).copied();
        let row = if let Some(event) = release {
            json!({
                "site": site,
                "kind": kind,
                "origin": "original-lattice",
                "initial_state": initial,
                "eligible_state": states.last().unwrap(),
                "reachable": true,
                "release_event_index": event,
                "legal_path": {"event_prefix_inclusive": [0, event]},
                "missing_predicates": [],
            })
        } else {
            let mut initial_bridges = 0_usize;
            let mut frozen_bridges = 0_usize;
            for &neighbor in engine.lattice.neighbors(site) {
                let neighbor = neighbor as usize;
                let oxygen_kind = &kinds[neighbor];
                let state = &deck.state_names[initial_states[neighbor].0 as usize];
                if compatible_bridge(kind, oxygen_kind, state) {
                    if engine.lattice.frozen[neighbor] {
                        frozen_bridges += 1;
                    } else {
                        initial_bridges += 1;
                    }
                }
            }
            let max_index = (initial_index + initial_bridges).min(states.len() - 1);
            if max_index == states.len() - 1 {
                return Err(format!(
                    "site {site} was not released but structural upper bound did not prove NO-GO"
                )
                .into());
            }
            unreachable.insert(site);
            json!({
                "site": site,
                "kind": kind,
                "origin": "original-lattice",
                "initial_state": initial,
                "eligible_state": states.last().unwrap(),
                "reachable": false,
                "release_event_index": null,
                "legal_path": null,
                "structural_upper_bound": {
                    "initial_ladder_index": initial_index,
                    "nonfrozen_initial_bridge_increments": initial_bridges,
                    "frozen_initial_bridges": frozen_bridges,
                    "maximum_state": states[max_index],
                },
                "missing_predicates": [format!(
                    "requires {}, but frozen/vacant boundary oxygen topology caps this center at {}; reverse condensation and adsorption can only restore previously removed bridge increments",
                    states.last().unwrap(), states[max_index]
                )],
            })
        };
        target_rows.push(row);
    }
    target_rows.sort_by_key(|row| row["site"].as_u64().unwrap());
    for sites in by_kind.values_mut() {
        sites.sort_unstable();
    }

    let reachable_count = release_event.len();
    let report = json!({
        "schema": "a9b-mechanism-reachability-v2",
        "deck": display_path(&deck_path),
        "deck_sha256": sha256(&deck_bytes),
        "snapshot": display_path(&snapshot_path),
        "snapshot_sha256": sha256(&snapshot_bytes),
        "temperature_kelvin": deck.temperature,
        "lineage_anchor": "initially nonfrozen occupied Si/Al site id in the exact zero-step snapshot",
        "proof_method": {
            "live_semantics": "every witness event is applied through petra_core::Engine::step_with using the compiled production deck",
            "closure": "desorb eligible original targets, then exhaust all live forward hydrolysis transitions in reaction/site order",
            "bound": MAX_EVENTS,
            "no_go_rule": "a forward reaction can remove each initially intact adjacent bridge at most once; reverse condensation restores the same increment and adsorption-created oxygen states have zero net additional bridge-removal capacity",
        },
        "target_count": original.len(),
        "reachable_count": reachable_count,
        "unreachable_count": unreachable.len(),
        "original_lattice_sites": by_kind,
        "legal_events": events,
        "rare_event_method": {
            "name": "finite-horizon reaction-tilted CTMC importance sampling",
            "implementation": "petra_core::BiasedCtmc + Engine::advance_biased_ctmc_until",
            "path_reweighting_rule": "log(P/Q) += (Lambda_biased-Lambda_physical)*dt - ln(factor_fired); deadline censor segments omit the fired-event term",
            "estimator_contract": "for preregistered T and N, validate finite log weights, compute ESS=(sum w)^2/sum(w^2), and estimate mean(w*N_original_lattice_release(T))/T; zero releases remain exactly zero and cannot create a rank",
            "stopping_rule": "fixed physical-time horizon T and fixed replica count N; each replica stops at T or a declared numerical failure",
            "failure_or_censor": "nonfinite weights or failed preregistered ESS invalidate the scenario; zero lineage-valid releases are censored/unrankable",
        },
        "targets": target_rows,
    });
    if let Some(parent) = out_path.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(&out_path, serde_json::to_vec_pretty(&report)?)?;
    println!(
        "wrote {}: {reachable_count}/{} reachable, {} topology NO-GO, {} legal events",
        out_path.display(),
        original.len(),
        unreachable.len(),
        report["legal_events"].as_array().unwrap().len()
    );
    Ok(())
}

fn display_path(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}
