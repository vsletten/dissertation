//! Bounded finite-horizon A9b importance-sampling runner.
//!
//! Usage: cargo run -p petra-deck --example a9b_importance_sampling -- DECK OUT [N] [T]

use std::collections::{BTreeMap, HashMap};
use std::env;
use std::fs;
use std::path::PathBuf;

use petra_core::{BiasedCtmc, CtmcAdvance};
use serde_json::json;

const DEFAULT_REPLICAS: usize = 8;
const DEFAULT_HORIZON: f64 = 1.0e-6;
const BIAS_FACTOR: f64 = 1.0e6;
const MAX_EVENTS_PER_REPLICA: usize = 1_000_000;
const MIN_ESS_FRACTION: f64 = 0.1;

fn representable_weight(log_weight: f64) -> Option<f64> {
    if !log_weight.is_finite() || log_weight < f64::MIN_POSITIVE.ln() || log_weight > f64::MAX.ln()
    {
        return None;
    }
    let weight = log_weight.exp();
    weight.is_normal().then_some(weight)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = env::args().skip(1).collect();
    if !(2..=4).contains(&args.len()) {
        return Err("usage: a9b_importance_sampling DECK OUT [N] [T]".into());
    }
    let deck_path = PathBuf::from(&args[0]);
    let out_path = PathBuf::from(&args[1]);
    let replicas = args
        .get(2)
        .map(|value| value.parse())
        .transpose()?
        .unwrap_or(DEFAULT_REPLICAS);
    let horizon = args
        .get(3)
        .map(|value| value.parse())
        .transpose()?
        .unwrap_or(DEFAULT_HORIZON);
    if replicas == 0 || !horizon.is_finite() || horizon <= 0.0 {
        return Err("N and T must be finite and positive".into());
    }

    let deck = petra_deck::load(&deck_path)?;
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
        .map(|(reaction, _)| (reaction as u16, BIAS_FACTOR))
        .collect();

    let mut rows = Vec::new();
    let mut failure = None;
    for replica in 0..replicas {
        let seed = deck.seed.wrapping_add(replica as u64);
        let mut engine = deck.build_engine(Some(seed))?;
        // This is an ion ledger, not a permanent site label. Consuming its
        // entry on release prevents reservoir re-adsorption onto a vacated
        // lattice site from becoming a second lattice-origin release.
        let mut lattice_ions: BTreeMap<_, _> = engine
            .lattice
            .states
            .iter()
            .enumerate()
            .filter_map(|(site, state)| {
                if engine.lattice.frozen[site] {
                    return None;
                }
                match deck.state_occupants[state.0 as usize].as_deref() {
                    Some(kind @ ("Al" | "Si")) => Some((site, kind.to_owned())),
                    _ => None,
                }
            })
            .collect();
        let mut sampler = BiasedCtmc::new(deck.reactions.len(), overrides.iter().copied())?;
        let mut release_count = 0_u64;
        let mut releases_by_kind = BTreeMap::from([("Al", 0_u64), ("Si", 0_u64)]);
        let mut event_count = 0_usize;

        loop {
            if event_count == MAX_EVENTS_PER_REPLICA {
                failure = Some(format!(
                    "replica {replica} exceeded {MAX_EVENTS_PER_REPLICA} events before T"
                ));
                break;
            }
            match engine.advance_biased_ctmc_until(&mut sampler, horizon) {
                Ok(CtmcAdvance::Fired(event)) => {
                    event_count += 1;
                    let expected_kind = if event.reaction == desorb_al {
                        Some("Al")
                    } else if event.reaction == desorb_si {
                        Some("Si")
                    } else {
                        None
                    };
                    if let Some(kind) = expected_kind {
                        if lattice_ions.get(&event.site).map(String::as_str) == Some(kind) {
                            release_count += 1;
                            *releases_by_kind.get_mut(kind).unwrap() += 1;
                            lattice_ions.remove(&event.site);
                        }
                    }
                }
                Ok(CtmcAdvance::Deadline { .. }) => break,
                Err(error) => {
                    failure = Some(format!("replica {replica} failed: {error}"));
                    break;
                }
            }
        }
        let log_weight = sampler.log_likelihood_ratio();
        if representable_weight(log_weight).is_none() {
            failure = Some(format!(
                "replica {replica} has a nonfinite representable weight"
            ));
        }
        rows.push(json!({
            "replica": replica,
            "seed": seed,
            "event_count": event_count,
            "original_lattice_release_count": release_count,
            "releases_by_kind": releases_by_kind,
            "log_likelihood_ratio": log_weight,
        }));
        if failure.is_some() {
            break;
        }
    }

    let mut estimate = None;
    let mut ess = None;
    let mut status = "invalid";
    if failure.is_none() && rows.len() == replicas {
        let weights: Vec<_> = rows
            .iter()
            .map(|row| {
                representable_weight(row["log_likelihood_ratio"].as_f64().unwrap())
                    .expect("weights were validated before aggregation")
            })
            .collect();
        let sum_w: f64 = weights.iter().sum();
        let sum_w2: f64 = weights.iter().map(|weight| weight * weight).sum();
        if !sum_w.is_finite() || !sum_w2.is_finite() || sum_w <= 0.0 || sum_w2 <= 0.0 {
            failure = Some("importance weights underflowed or overflowed".into());
        } else {
            let computed_ess = sum_w * sum_w / sum_w2;
            ess = Some(computed_ess);
            if computed_ess / replicas as f64 >= MIN_ESS_FRACTION {
                let weighted_releases: f64 = rows
                    .iter()
                    .zip(weights.iter())
                    .map(|(row, weight)| {
                        weight * row["original_lattice_release_count"].as_u64().unwrap() as f64
                    })
                    .sum();
                let computed_estimate = weighted_releases / replicas as f64 / horizon;
                if !weighted_releases.is_finite() || !computed_estimate.is_finite() {
                    failure = Some("weighted release estimator overflowed".into());
                } else {
                    estimate = Some(computed_estimate);
                    status = if weighted_releases == 0.0 {
                        "censored-zero-observation"
                    } else {
                        "estimated"
                    };
                }
            } else {
                failure = Some(format!(
                    "effective sample size {computed_ess} is below preregistered fraction {MIN_ESS_FRACTION}"
                ));
            }
        }
    }

    let result = json!({
        "schema": "a9b-importance-sampling-v1",
        "deck": deck_path,
        "bias_factor": BIAS_FACTOR,
        "biased_reactions": overrides.iter().map(|(id, _)| &deck.reactions[*id as usize].name).collect::<Vec<_>>(),
        "replicas_preregistered": replicas,
        "horizon_preregistered": horizon,
        "max_events_per_replica": MAX_EVENTS_PER_REPLICA,
        "minimum_ess_fraction": MIN_ESS_FRACTION,
        "status": status,
        "failure": failure,
        "effective_sample_size": ess,
        "physical_release_rate_estimate": estimate,
        "rank_eligible": false,
        "replica_results": rows,
    });
    if let Some(parent) = out_path.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(&out_path, serde_json::to_vec_pretty(&result)?)?;
    println!("wrote {}: status={status}", out_path.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::representable_weight;

    #[test]
    fn importance_weights_fail_closed_before_overflow_or_underflow() {
        assert_eq!(representable_weight(0.0), Some(1.0));
        assert!(representable_weight(f64::NAN).is_none());
        assert!(representable_weight(1_000.0).is_none());
        assert!(representable_weight(-1_000.0).is_none());
    }
}
