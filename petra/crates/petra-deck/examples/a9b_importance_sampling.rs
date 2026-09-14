//! Durable single-replica A9b finite-horizon importance-sampling runner.
//!
//! Usage:
//! a9b_importance_sampling DECK PATH_REPORT SNAPSHOT OUT --seed N --horizon T
//!   [--checkpoints N] [--bias-factor F] [--max-events N]

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::env;
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use petra_core::{BiasedCtmc, CtmcAdvance, LikelihoodSegment};
use petra_observables::{observe, ObservableValue};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

const DEFAULT_HORIZON: f64 = 1.0e-6;
const DEFAULT_BIAS_FACTOR: f64 = 1.0e6;
const DEFAULT_CHECKPOINTS: usize = 21;
const DEFAULT_MAX_EVENTS: usize = 1_000_000;
const EXPECTED_PATH_SCHEMA: &str = "a9b-mechanism-reachability-v2";
const PACKAGE_SCHEMA: &str = "a9b-replica-package-v2";
const RECEIPT_SCHEMA: &str = "a9b-replica-receipt-v2";

#[derive(Debug)]
struct Config {
    deck: PathBuf,
    path_report: PathBuf,
    snapshot: PathBuf,
    out: PathBuf,
    seed: u64,
    horizon: f64,
    checkpoints: usize,
    bias_factor: f64,
    max_events: usize,
    argv: Vec<String>,
}

fn finite_positive(value: &str, label: &str) -> Result<f64, String> {
    let parsed: f64 = value
        .parse()
        .map_err(|_| format!("{label} must be a number"))?;
    if !parsed.is_finite() || parsed <= 0.0 {
        return Err(format!("{label} must be finite and positive"));
    }
    Ok(parsed)
}

fn parse_args(args: Vec<String>) -> Result<Config, String> {
    if args.len() < 8 {
        return Err(
            "usage: a9b_importance_sampling DECK PATH_REPORT SNAPSHOT OUT --seed N --horizon T [--checkpoints N] [--bias-factor F] [--max-events N]".into(),
        );
    }
    let argv = args.clone();
    let deck = PathBuf::from(&args[0]);
    let path_report = PathBuf::from(&args[1]);
    let snapshot = PathBuf::from(&args[2]);
    let out = PathBuf::from(&args[3]);
    let mut seed = None;
    let mut horizon = DEFAULT_HORIZON;
    let mut checkpoints = DEFAULT_CHECKPOINTS;
    let mut bias_factor = DEFAULT_BIAS_FACTOR;
    let mut max_events = DEFAULT_MAX_EVENTS;
    let mut index = 4;
    while index < args.len() {
        let flag = args[index].as_str();
        let value = args
            .get(index + 1)
            .ok_or_else(|| format!("{flag} requires a value"))?;
        match flag {
            "--seed" => {
                if seed.is_some() {
                    return Err("--seed may be specified only once".into());
                }
                seed = Some(value.parse().map_err(|_| "seed must be u64")?);
            }
            "--horizon" => horizon = finite_positive(value, "horizon")?,
            "--checkpoints" => {
                checkpoints = value
                    .parse()
                    .map_err(|_| "checkpoints must be an integer")?;
                if !(2..=1001).contains(&checkpoints) {
                    return Err("checkpoints must be in [2,1001]".into());
                }
            }
            "--bias-factor" => bias_factor = finite_positive(value, "bias factor")?,
            "--max-events" => {
                max_events = value.parse().map_err(|_| "max-events must be an integer")?;
                if max_events == 0 {
                    return Err("max-events must be positive".into());
                }
            }
            _ => return Err(format!("unknown argument {flag}")),
        }
        index += 2;
    }
    Ok(Config {
        deck,
        path_report,
        snapshot,
        out,
        seed: seed.ok_or("--seed is required")?,
        horizon,
        checkpoints,
        bias_factor,
        max_events,
        argv,
    })
}

fn sha256_bytes(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn sha256_file(path: &Path) -> Result<String, Box<dyn std::error::Error>> {
    Ok(sha256_bytes(&fs::read(path)?))
}

fn write_bytes(path: &Path, bytes: &[u8]) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = OpenOptions::new().write(true).create_new(true).open(path)?;
    file.write_all(bytes)?;
    file.sync_all()?;
    Ok(())
}

fn write_json(path: &Path, value: &Value) -> Result<(), Box<dyn std::error::Error>> {
    let mut bytes = serde_json::to_vec_pretty(value)?;
    bytes.push(b'\n');
    write_bytes(path, &bytes)
}

fn append_json_line(
    writer: &mut BufWriter<File>,
    value: &Value,
) -> Result<(), Box<dyn std::error::Error>> {
    serde_json::to_writer(&mut *writer, value)?;
    writer.write_all(b"\n")?;
    Ok(())
}

fn categorical<'a>(value: &'a Value, label: &str) -> Result<(&'a [Value], &'a [Value]), String> {
    let dict = value["dict"]
        .as_array()
        .ok_or_else(|| format!("{label}.dict is missing"))?;
    let data = value["data"]
        .as_array()
        .ok_or_else(|| format!("{label}.data is missing"))?;
    Ok((dict, data))
}

fn validate_snapshot(
    deck: &petra_deck::CompiledDeck,
    engine: &petra_core::Engine,
    snapshot: &Value,
) -> Result<(), String> {
    if snapshot["pgif"] != 1 || snapshot["meta"]["directed"] != false {
        return Err("snapshot must be undirected PGIF v1".into());
    }
    if snapshot["meta"]["petra"]["step"] != 0 || snapshot["meta"]["petra"]["time"] != 0.0 {
        return Err("lineage snapshot must be at step/time zero".into());
    }
    if snapshot["meta"]["petra"]["states"] != json!(deck.state_names) {
        return Err("snapshot state table does not match deck".into());
    }
    let nodes = &snapshot["nodes"];
    if nodes["count"].as_u64() != Some(engine.lattice.len() as u64) {
        return Err("snapshot node count does not match deck".into());
    }
    let (state_dict, states) = categorical(&nodes["columns"]["state"], "state")?;
    let (kind_dict, kinds) = categorical(&nodes["columns"]["kind"], "kind")?;
    let frozen = nodes["columns"]["frozen"]["data"]
        .as_array()
        .ok_or("snapshot frozen data is missing")?;
    if state_dict != json!(deck.state_names).as_array().unwrap()
        || kind_dict != json!(deck.kind_names).as_array().unwrap()
        || states.len() != engine.lattice.len()
        || kinds.len() != engine.lattice.len()
        || frozen.len() != engine.lattice.len()
    {
        return Err("snapshot categorical columns are malformed".into());
    }
    for site in 0..engine.lattice.len() {
        let template = engine.lattice.template_index[site] as usize;
        let expected_kind = deck.kinds_per_template[template].0 as u64;
        if states[site].as_u64() != Some(engine.lattice.states[site].0 as u64)
            || kinds[site].as_u64() != Some(expected_kind)
            || frozen[site].as_bool() != Some(engine.lattice.frozen[site])
        {
            return Err(format!("snapshot/live deck mismatch at site {site}"));
        }
    }
    let canonical: Value = serde_json::from_str(&petra_io::snapshot_json(deck, engine))
        .map_err(|error| format!("failed to construct canonical deck snapshot: {error}"))?;
    if snapshot != &canonical {
        return Err(
            "lineage snapshot differs from the full canonical deck PGIF (metadata, nodes, coordinates, or edges)"
                .into(),
        );
    }
    Ok(())
}

fn validate_path_report(
    report: &Value,
    report_hash: &str,
    snapshot_hash: &str,
    deck: &petra_deck::CompiledDeck,
    engine: &petra_core::Engine,
) -> Result<BTreeMap<usize, String>, String> {
    if report["schema"] != EXPECTED_PATH_SCHEMA
        || report["temperature_kelvin"] != 298.0
        || report["target_count"] != 320
        || report["reachable_count"] != 260
        || report["unreachable_count"] != 60
        || report["snapshot_sha256"] != snapshot_hash
    {
        return Err(format!(
            "path report census/hash contract failed ({report_hash})"
        ));
    }
    let targets = report["targets"]
        .as_array()
        .ok_or("path report targets are missing")?;
    if targets.len() != 320 {
        return Err("path report must contain exactly 320 targets".into());
    }
    let mut original = BTreeMap::new();
    let mut reachable = 0;
    let mut unreachable = 0;
    for target in targets {
        let site = target["site"].as_u64().ok_or("target site is missing")? as usize;
        let kind = target["kind"].as_str().ok_or("target kind is missing")?;
        if !matches!(kind, "Al" | "Si") || site >= engine.lattice.len() {
            return Err("target identity is invalid".into());
        }
        if engine.lattice.frozen[site]
            || deck.state_occupants[engine.lattice.states[site].0 as usize].as_deref() != Some(kind)
            || original.insert(site, kind.to_owned()).is_some()
        {
            return Err(format!(
                "target {site} is not a unique original lattice cation"
            ));
        }
        match target["reachable"].as_bool() {
            Some(true)
                if target["legal_path"].is_object()
                    && target["release_event_index"].is_u64()
                    && target["missing_predicates"]
                        .as_array()
                        .is_some_and(Vec::is_empty) =>
            {
                reachable += 1;
            }
            Some(false)
                if target["legal_path"].is_null()
                    && target["release_event_index"].is_null()
                    && target["structural_upper_bound"].is_object()
                    && target["missing_predicates"]
                        .as_array()
                        .is_some_and(|items| !items.is_empty()) =>
            {
                unreachable += 1;
            }
            Some(_) => return Err("target reachability evidence is incomplete".into()),
            None => return Err("target reachable value must be boolean".into()),
        }
    }
    if original.len() != 320 || reachable != 260 || unreachable != 60 {
        return Err("path report target census disagrees with summary".into());
    }
    Ok(original)
}

fn state_map(deck: &petra_deck::CompiledDeck, engine: &petra_core::Engine) -> Value {
    let counts = engine.state_counts(deck.n_states);
    Value::Object(
        deck.state_names
            .iter()
            .cloned()
            .zip(counts)
            .map(|(name, count)| (name, json!(count)))
            .collect(),
    )
}

fn kind_totals(deck: &petra_deck::CompiledDeck, engine: &petra_core::Engine) -> Value {
    let mut totals = BTreeMap::<String, u64>::new();
    for &template in &engine.lattice.template_index {
        let kind = deck.kind_names[deck.kinds_per_template[template as usize].0 as usize].clone();
        *totals.entry(kind).or_default() += 1;
    }
    json!(totals)
}

fn geometric_area(
    deck: &petra_deck::CompiledDeck,
    engine: &petra_core::Engine,
) -> Result<f64, String> {
    let sample = observe(engine, deck);
    let areas: Vec<_> = sample
        .values
        .iter()
        .filter_map(|value| match value {
            ObservableValue::SurfaceArea(area) => Some(area.geometric),
            _ => None,
        })
        .collect();
    if areas.len() != 1 || !areas[0].is_finite() || areas[0] <= 0.0 {
        return Err("deck must emit exactly one positive finite geometric surface area".into());
    }
    Ok(areas[0])
}

fn segment_json(index: usize, segment: &LikelihoodSegment) -> Result<Value, String> {
    let values = [
        segment.start_time,
        segment.end_time,
        segment.physical_total_rate,
        segment.biased_total_rate,
        segment.log_likelihood_increment,
        segment.cumulative_log_likelihood,
    ];
    if values.iter().any(|value| !value.is_finite()) {
        return Err("likelihood segment contains a nonfinite value".into());
    }
    Ok(json!({
        "segment": index,
        "start_time_s": segment.start_time,
        "end_time_s": segment.end_time,
        "physical_total_rate_s-1": segment.physical_total_rate,
        "biased_total_rate_s-1": segment.biased_total_rate,
        "fired_reaction_id": segment.fired_reaction,
        "fired_bias_factor": segment.fired_bias_factor,
        "log_likelihood_increment": segment.log_likelihood_increment,
        "cumulative_log_likelihood": segment.cumulative_log_likelihood,
    }))
}

fn sync_writer(mut writer: BufWriter<File>) -> Result<(), Box<dyn std::error::Error>> {
    writer.flush()?;
    writer.get_ref().sync_all()?;
    Ok(())
}

fn run(config: &Config, temp: &Path) -> Result<Value, Box<dyn std::error::Error>> {
    let deck_bytes = fs::read(&config.deck)?;
    let report_bytes = fs::read(&config.path_report)?;
    let snapshot_bytes = fs::read(&config.snapshot)?;
    let report_hash = sha256_bytes(&report_bytes);
    let snapshot_hash = sha256_bytes(&snapshot_bytes);
    let report: Value = serde_json::from_slice(&report_bytes)?;
    let reference_snapshot: Value = serde_json::from_slice(&snapshot_bytes)?;
    let deck_source = std::str::from_utf8(&deck_bytes)?;
    let deck_file: petra_deck::DeckFile = toml::from_str(deck_source)?;
    let deck = petra_deck::compile(&deck_file)?;
    if deck.temperature != 298.0 {
        return Err(format!("A9b requires exact 298 K, got {}", deck.temperature).into());
    }
    let mut engine = deck.build_engine(Some(config.seed))?;
    validate_snapshot(&deck, &engine, &reference_snapshot)?;
    let original_ions =
        validate_path_report(&report, &report_hash, &snapshot_hash, &deck, &engine)?;
    let mut ion_ledger = original_ions.clone();

    let initial_pgif = petra_io::snapshot_json(&deck, &engine);
    write_bytes(&temp.join("initial.pgif.json"), initial_pgif.as_bytes())?;

    let reaction_ids: HashMap<_, _> = deck
        .reactions
        .iter()
        .enumerate()
        .map(|(index, reaction)| (reaction.name.as_str(), index as u16))
        .collect();
    let desorb_al = *reaction_ids.get("desorb-al").ok_or("missing desorb-al")?;
    let desorb_si = *reaction_ids.get("desorb-si").ok_or("missing desorb-si")?;
    let overrides: Vec<_> = deck
        .reactions
        .iter()
        .enumerate()
        .filter(|(_, reaction)| {
            reaction.name.contains("hydrolysis") || reaction.name.starts_with("desorb-")
        })
        .map(|(reaction, _)| (reaction as u16, config.bias_factor))
        .collect();
    let override_ids: BTreeSet<_> = overrides.iter().map(|(id, _)| *id).collect();
    let mut sampler = BiasedCtmc::new(deck.reactions.len(), overrides.iter().copied())?;

    let events_file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(temp.join("events.jsonl"))?;
    let likelihood_file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(temp.join("likelihood.jsonl"))?;
    let checkpoints_file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(temp.join("checkpoints.jsonl"))?;
    let mut events = BufWriter::new(events_file);
    let mut likelihood = BufWriter::new(likelihood_file);
    let mut checkpoints = BufWriter::new(checkpoints_file);

    append_json_line(
        &mut events,
        &json!({
            "schema": "a9b-event-stream-v2",
            "seed": config.seed,
            "states": deck.state_names,
            "state_occupants": deck.state_occupants,
            "reactions": deck.reactions.iter().map(|reaction| reaction.name.as_str()).collect::<Vec<_>>(),
        }),
    )?;
    append_json_line(
        &mut likelihood,
        &json!({
            "schema": "a9b-likelihood-stream-v2",
            "seed": config.seed,
            "formula": "(biased_total_rate-physical_total_rate)*(end_time-start_time)-ln(fired_bias_factor); deadline segments omit the event term",
        }),
    )?;

    let mut releases = BTreeMap::from([("Al", 0_u64), ("Si", 0_u64)]);
    let mut event_count = 0usize;
    let mut segment_count = 0usize;
    let fixed_kind_totals = kind_totals(&deck, &engine);
    let initial_area = geometric_area(&deck, &engine)?;
    append_json_line(
        &mut checkpoints,
        &json!({
            "checkpoint": 0,
            "checkpoint_role": "fixed_physical_time",
            "target_time_s": 0.0,
            "actual_time_s": engine.time,
            "step": engine.step_count,
            "geometric_area_a2": initial_area,
            "state_counts": state_map(&deck, &engine),
            "kind_totals": fixed_kind_totals,
            "original_lattice_releases": releases,
            "cumulative_log_likelihood": sampler.log_likelihood_ratio(),
        }),
    )?;

    let mut stop_reason = "physical_horizon".to_owned();
    let mut complete = true;
    let mut recorded_checkpoint_count = 1_usize;
    let mut last_checkpoint_time = engine.time;
    'checkpoints: for checkpoint in 1..config.checkpoints {
        let deadline = config.horizon * checkpoint as f64 / (config.checkpoints - 1) as f64;
        loop {
            if event_count >= config.max_events {
                stop_reason = "max_events".to_owned();
                complete = false;
                break 'checkpoints;
            }
            let advance = match engine.advance_biased_ctmc_until(&mut sampler, deadline) {
                Ok(value) => value,
                Err(error) => {
                    stop_reason = format!("engine_error:{error}");
                    complete = false;
                    break 'checkpoints;
                }
            };
            let segment = *sampler
                .last_segment()
                .ok_or("successful advance did not retain a likelihood segment")?;
            append_json_line(&mut likelihood, &segment_json(segment_count, &segment)?)?;
            let current_segment = segment_count;
            segment_count += 1;
            match advance {
                CtmcAdvance::Fired(fired) => {
                    event_count += 1;
                    let reaction_name = &deck.reactions[fired.reaction as usize].name;
                    let expected_kind = if fired.reaction == desorb_al {
                        Some("Al")
                    } else if fired.reaction == desorb_si {
                        Some("Si")
                    } else {
                        None
                    };
                    let lattice_release_kind = expected_kind.and_then(|kind| {
                        (ion_ledger.get(&fired.site).map(String::as_str) == Some(kind)).then(|| {
                            ion_ledger.remove(&fired.site);
                            *releases.get_mut(kind).expect("known cation kind") += 1;
                            kind
                        })
                    });
                    let changes = engine
                        .last_changes()
                        .iter()
                        .map(|(site, old, new)| {
                            json!({
                                "site": site,
                                "old_state_id": old.0,
                                "new_state_id": new.0,
                                "old_state": deck.state_names[old.0 as usize],
                                "new_state": deck.state_names[new.0 as usize],
                            })
                        })
                        .collect::<Vec<_>>();
                    append_json_line(
                        &mut events,
                        &json!({
                            "event": event_count - 1,
                            "step": fired.step,
                            "time_s": fired.time,
                            "reaction_id": fired.reaction,
                            "reaction": reaction_name,
                            "center_site": fired.site,
                            "changes": changes,
                            "likelihood_segment": current_segment,
                            "lattice_release_kind": lattice_release_kind,
                        }),
                    )?;
                }
                CtmcAdvance::Deadline { time } => {
                    if time != deadline {
                        return Err("engine deadline timestamp drifted".into());
                    }
                    break;
                }
            }
        }
        append_json_line(
            &mut checkpoints,
            &json!({
                "checkpoint": checkpoint,
                "checkpoint_role": "fixed_physical_time",
                "target_time_s": deadline,
                "actual_time_s": engine.time,
                "step": engine.step_count,
                "geometric_area_a2": geometric_area(&deck, &engine)?,
                "state_counts": state_map(&deck, &engine),
                "kind_totals": fixed_kind_totals,
                "original_lattice_releases": releases,
                "cumulative_log_likelihood": sampler.log_likelihood_ratio(),
            }),
        )?;
        recorded_checkpoint_count += 1;
        last_checkpoint_time = engine.time;
    }
    if !complete && engine.time > last_checkpoint_time {
        append_json_line(
            &mut checkpoints,
            &json!({
                "checkpoint": recorded_checkpoint_count,
                "checkpoint_role": "terminal_censor_snapshot",
                "target_time_s": engine.time,
                "actual_time_s": engine.time,
                "step": engine.step_count,
                "geometric_area_a2": geometric_area(&deck, &engine)?,
                "state_counts": state_map(&deck, &engine),
                "kind_totals": fixed_kind_totals,
                "original_lattice_releases": releases,
                "cumulative_log_likelihood": sampler.log_likelihood_ratio(),
            }),
        )?;
        recorded_checkpoint_count += 1;
    }
    sync_writer(events)?;
    sync_writer(likelihood)?;
    sync_writer(checkpoints)?;

    let final_pgif = petra_io::snapshot_json(&deck, &engine);
    write_bytes(&temp.join("final.pgif.json"), final_pgif.as_bytes())?;

    let reaction_map = deck
        .reactions
        .iter()
        .enumerate()
        .map(|(id, reaction)| {
            json!({
                "id": id,
                "name": reaction.name,
                "bias_factor": if override_ids.contains(&(id as u16)) { config.bias_factor } else { 1.0 },
            })
        })
        .collect::<Vec<_>>();
    let executable = env::current_exe()?;
    let metadata = json!({
        "schema": PACKAGE_SCHEMA,
        "survey_tier_platform_test": true,
        "not_production": true,
        "seed": config.seed,
        "temperature_k": deck.temperature,
        "horizon_s": config.horizon,
        "fixed_checkpoint_count": config.checkpoints,
        "recorded_checkpoint_count": recorded_checkpoint_count,
        "bias_factor": config.bias_factor,
        "max_events": config.max_events,
        "complete": complete,
        "stop_reason": stop_reason,
        "event_count": event_count,
        "likelihood_segment_count": segment_count,
        "final_time_s": engine.time,
        "final_step": engine.step_count,
        "final_log_likelihood": sampler.log_likelihood_ratio(),
        "original_lattice_releases": releases,
        "remaining_original_ion_ledger": ion_ledger.len(),
        "reaction_mapping": reaction_map,
        "census": {
            "whole_finite_deck_original_centers": 320,
            "reachable_original_centers": 260,
            "topology_no_go_original_centers": 60,
            "whole_deck_estimand_includes_topology_no_go": true,
        },
        "inputs": {
            "deck": config.deck,
            "deck_sha256": sha256_bytes(&deck_bytes),
            "path_report": config.path_report,
            "path_report_sha256": report_hash,
            "lineage_snapshot": config.snapshot,
            "lineage_snapshot_sha256": snapshot_hash,
            "runner_executable": executable,
            "runner_executable_sha256": sha256_file(&executable)?,
        },
        "argv": config.argv,
    });
    write_json(&temp.join("metadata.json"), &metadata)?;
    Ok(metadata)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let config = parse_args(env::args().skip(1).collect())?;
    if config.out.exists() {
        return Err(format!(
            "refusing to overwrite existing output {}",
            config.out.display()
        )
        .into());
    }
    let parent = config.out.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(parent)?;
    let name = config
        .out
        .file_name()
        .and_then(|value| value.to_str())
        .ok_or("output requires a normal final path component")?;
    let temp = parent.join(format!(".{name}.tmp-{}", std::process::id()));
    fs::create_dir(&temp)?;

    let result = run(&config, &temp);
    if let Err(error) = result {
        let _ = fs::remove_dir_all(&temp);
        return Err(error);
    }
    let metadata = result.expect("checked above");
    let files = [
        "metadata.json",
        "initial.pgif.json",
        "events.jsonl",
        "likelihood.jsonl",
        "checkpoints.jsonl",
        "final.pgif.json",
    ];
    let hashes = files
        .iter()
        .map(|name| Ok((name.to_string(), sha256_file(&temp.join(name))?)))
        .collect::<Result<BTreeMap<_, _>, Box<dyn std::error::Error>>>()?;
    let receipt = json!({
        "schema": RECEIPT_SCHEMA,
        "status": if metadata["complete"] == true { "complete" } else { "censored" },
        "seed": config.seed,
        "files_sha256": hashes,
    });
    write_json(&temp.join("receipt.json"), &receipt)?;
    File::open(&temp)?.sync_all()?;
    if config.out.exists() {
        fs::remove_dir_all(&temp)?;
        return Err("output appeared during run; refusing to overwrite".into());
    }
    fs::rename(&temp, &config.out)?;
    File::open(parent)?.sync_all()?;
    println!(
        "wrote {}: seed={} events={} status={}",
        config.out.display(),
        config.seed,
        metadata["event_count"],
        receipt["status"]
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{finite_positive, parse_args};

    #[test]
    fn explicit_seed_is_required_and_numerics_fail_closed() {
        let base = vec!["d", "p", "s", "o", "--horizon", "1e-6", "--seed", "90401"]
            .into_iter()
            .map(str::to_owned)
            .collect();
        let parsed = parse_args(base).expect("valid arguments");
        assert_eq!(parsed.seed, 90401);
        assert_eq!(parsed.checkpoints, 21);
        assert!(finite_positive("nan", "x").is_err());
        assert!(finite_positive("0", "x").is_err());
        assert!(parse_args(
            vec!["d", "p", "s", "o", "--horizon", "1e-6", "--max-events", "0"]
                .into_iter()
                .map(str::to_owned)
                .collect()
        )
        .is_err());
    }
}
