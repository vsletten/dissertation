use std::io::Write;
use std::path::{Path, PathBuf};

use petra_core::{CtmcAdvance, StateId};

const POTENTIALS: [f64; 7] = [-0.70, -0.45, -0.20, 0.05, 0.30, 0.55, 0.80];
const REPLICAS: u64 = 32;
const SCALING_SIZES: [usize; 3] = [16, 24, 32];
const SCALING_POTENTIALS: [f64; 3] = [0.30, 0.55, 0.80];
const SCALING_REPLICAS: u64 = 16;
const HORIZON: f64 = 600.0;
const MAX_EVENTS: u64 = 40_000;
const STABLE_CLUSTER: u64 = 6;

#[derive(Debug)]
struct ReplicaResult {
    seed: u64,
    induction_time: Option<f64>,
    final_active: u64,
    pit_count: u64,
    largest_cluster: u64,
    nn_pair_ratio: f64,
    dissolution_events: u64,
    repassivation_events: u64,
}

#[derive(Debug, Clone, Copy)]
struct ClusterStats {
    count: u64,
    largest: u64,
}

fn compiled_at(text: &str, potential_mu: f64) -> Result<petra_deck::CompiledDeck, String> {
    compiled_at_size(text, potential_mu, 24)
}

fn compiled_at_size(
    text: &str,
    potential_mu: f64,
    size: usize,
) -> Result<petra_deck::CompiledDeck, String> {
    let mut parsed: petra_deck::DeckFile = toml::from_str(text).map_err(|e| e.to_string())?;
    parsed.thermo.mu.insert("e_drive".to_string(), potential_mu);
    parsed.grid.as_mut().expect("passive-film grid").dims = vec![size, size];
    parsed.defects[0].at = [size as f64 / 2.0, size as f64 / 2.0, 0.0];
    petra_deck::compile(&parsed).map_err(|e| e.to_string())
}

fn state_ids(deck: &petra_deck::CompiledDeck, names: &[&str]) -> Vec<StateId> {
    names
        .iter()
        .map(|name| {
            StateId(
                deck.state_names
                    .iter()
                    .position(|candidate| candidate == name)
                    .unwrap_or_else(|| panic!("deck is missing state {name}"))
                    as u16,
            )
        })
        .collect()
}

fn is_selected(state: StateId, selected: &[StateId]) -> bool {
    selected.contains(&state)
}

fn active_count(engine: &petra_core::Engine, active: &[StateId]) -> u64 {
    engine
        .lattice
        .states
        .iter()
        .filter(|&&state| is_selected(state, active))
        .count() as u64
}

fn cluster_stats(engine: &petra_core::Engine, active: &[StateId]) -> ClusterStats {
    let mut seen = vec![false; engine.lattice.len()];
    let mut stack = Vec::new();
    let mut largest = 0_u64;
    let mut count = 0_u64;
    for start in 0..engine.lattice.len() {
        if seen[start] || !is_selected(engine.lattice.states[start], active) {
            continue;
        }
        seen[start] = true;
        stack.push(start);
        count += 1;
        let mut size = 0_u64;
        while let Some(site) = stack.pop() {
            size += 1;
            for &neighbor in engine.lattice.neighbors(site) {
                let neighbor = neighbor as usize;
                if !seen[neighbor] && is_selected(engine.lattice.states[neighbor], active) {
                    seen[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
        largest = largest.max(size);
    }
    ClusterStats { count, largest }
}

// Nearest-neighbor pair correlation relative to a random surface with the
// same active fraction. Values above one mean pit patches cluster spatially.
fn nn_pair_ratio(engine: &petra_core::Engine, active: &[StateId]) -> f64 {
    let active_sites = active_count(engine, active);
    if active_sites < 2 {
        return 0.0;
    }
    let mut total_edges = 0_u64;
    let mut active_edges = 0_u64;
    for site in 0..engine.lattice.len() {
        for &neighbor in engine.lattice.neighbors(site) {
            let neighbor = neighbor as usize;
            if neighbor <= site {
                continue;
            }
            total_edges += 1;
            active_edges += u64::from(
                is_selected(engine.lattice.states[site], active)
                    && is_selected(engine.lattice.states[neighbor], active),
            );
        }
    }
    let fraction = active_sites as f64 / engine.lattice.len() as f64;
    (active_edges as f64 / total_edges as f64) / (fraction * fraction)
}

fn reaction_rate(engine: &petra_core::Engine, reaction: u16) -> f64 {
    engine
        .site_event_rates()
        .iter()
        .flat_map(|events| events.iter())
        .filter_map(|&(rule, rate)| (rule == reaction).then_some(rate))
        .sum()
}

fn run_replica(
    deck: &petra_deck::CompiledDeck,
    seed: u64,
    mut transient: Option<&mut dyn Write>,
) -> Result<ReplicaResult, String> {
    let mut engine = deck
        .build_engine(Some(seed))
        .map_err(|error| error.to_string())?;
    let pit = state_ids(deck, &["film_patch.bare", "film_patch.dissolving"]);
    let names: Vec<_> = deck
        .reactions
        .iter()
        .map(|reaction| reaction.name.as_str())
        .collect();
    let dissolution_rule = names
        .iter()
        .position(|name| *name == "metal_dissolution")
        .expect("metal_dissolution rule") as u16;
    let mut induction_time = None;
    let mut dissolution_events = 0_u64;
    let mut repassivation_events = 0_u64;

    if let Some(writer) = transient.as_mut() {
        writeln!(
            writer,
            "event,time_s,active_patches,largest_cluster,current_propensity,dissolution_events,repassivation_events"
        )
        .map_err(|error| error.to_string())?;
    }

    let mut event = 0_u64;
    while event < MAX_EVENTS {
        let fired = match engine
            .advance_ctmc_until(HORIZON)
            .map_err(|stop| stop.to_string())?
        {
            CtmcAdvance::Fired(fired) => fired,
            CtmcAdvance::Deadline { time } => {
                debug_assert_eq!(time.to_bits(), HORIZON.to_bits());
                break;
            }
        };
        event += 1;
        let name = names[fired.reaction as usize];
        if name == "metal_dissolution" {
            dissolution_events += 1;
        } else if name.starts_with("repassivate_") {
            repassivation_events += 1;
        }

        let sample = event.is_multiple_of(20) || name == "film_rupture";
        if sample {
            let clusters = cluster_stats(&engine, &pit);
            if induction_time.is_none() && clusters.largest >= STABLE_CLUSTER {
                induction_time = Some(engine.time);
            }
            if let Some(writer) = transient.as_mut() {
                writeln!(
                    writer,
                    "{event},{:.9},{},{},{:.9},{dissolution_events},{repassivation_events}",
                    engine.time,
                    active_count(&engine, &pit),
                    clusters.largest,
                    reaction_rate(&engine, dissolution_rule),
                )
                .map_err(|error| error.to_string())?;
            }
        }
    }

    if engine.time < HORIZON {
        return Err(format!(
            "seed {seed}: bounded event cap reached at t={:.3} < {HORIZON}",
            engine.time
        ));
    }
    debug_assert_eq!(engine.time.to_bits(), HORIZON.to_bits());
    let clusters = cluster_stats(&engine, &pit);
    if let Some(writer) = transient.as_mut() {
        writeln!(
            writer,
            "{event},{:.9},{},{},{:.9},{dissolution_events},{repassivation_events}",
            engine.time,
            active_count(&engine, &pit),
            clusters.largest,
            reaction_rate(&engine, dissolution_rule),
        )
        .map_err(|error| error.to_string())?;
    }
    Ok(ReplicaResult {
        seed,
        induction_time,
        final_active: active_count(&engine, &pit),
        pit_count: clusters.count,
        largest_cluster: clusters.largest,
        nn_pair_ratio: nn_pair_ratio(&engine, &pit),
        dissolution_events,
        repassivation_events,
    })
}

// Descriptive median-rank Weibull plot: ln(-ln(1-F)) = k ln(t) - k ln(lambda).
// This ignores right-censored replicas and is not a distribution-selection
// test; report both its slope and R² so the limitation remains machine-visible.
fn weibull_plot(times: &[f64]) -> Option<(f64, f64)> {
    if times.len() < 3 {
        return None;
    }
    let mut times = times.to_vec();
    times.sort_by(f64::total_cmp);
    let n = times.len() as f64;
    let points: Vec<_> = times
        .iter()
        .enumerate()
        .map(|(index, &time)| {
            let rank = (index as f64 + 0.7) / (n + 0.4);
            (time.ln(), (-(1.0 - rank).ln()).ln())
        })
        .collect();
    let x_mean = points.iter().map(|point| point.0).sum::<f64>() / n;
    let y_mean = points.iter().map(|point| point.1).sum::<f64>() / n;
    let numerator = points
        .iter()
        .map(|(x, y)| (x - x_mean) * (y - y_mean))
        .sum::<f64>();
    let denominator = points
        .iter()
        .map(|(x, _)| (x - x_mean).powi(2))
        .sum::<f64>();
    if denominator <= 0.0 {
        return None;
    }
    let slope = numerator / denominator;
    let intercept = y_mean - slope * x_mean;
    let residual = points
        .iter()
        .map(|(x, y)| (y - (intercept + slope * x)).powi(2))
        .sum::<f64>();
    let total = points
        .iter()
        .map(|(_, y)| (y - y_mean).powi(2))
        .sum::<f64>();
    Some((slope, 1.0 - residual / total))
}

fn mean(values: impl Iterator<Item = f64>) -> f64 {
    let values: Vec<_> = values.collect();
    values.iter().sum::<f64>() / values.len() as f64
}

fn sample_variance(values: &[f64]) -> f64 {
    let average = mean(values.iter().copied());
    values
        .iter()
        .map(|value| (value - average).powi(2))
        .sum::<f64>()
        / (values.len() - 1) as f64
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let deck_path = args
        .next()
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("decks/passive-metal.toml"));
    let output_dir = args
        .next()
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("corrosion-study"));
    if args.next().is_some() {
        return Err("usage: corrosion_study [deck.toml] [output-dir]".into());
    }

    let text = std::fs::read_to_string(&deck_path)?;
    std::fs::create_dir_all(&output_dir)?;
    let mut summary = std::io::BufWriter::new(std::fs::File::create(
        output_dir.join("potential-summary.csv"),
    )?);
    let mut induction = std::io::BufWriter::new(std::fs::File::create(
        output_dir.join("induction-times.csv"),
    )?);
    let mut scaling = std::io::BufWriter::new(std::fs::File::create(
        output_dir.join("critical-scaling.csv"),
    )?);
    writeln!(summary, "potential_mu_kcal_mol,replicas,initiated,initiation_fraction,weibull_plot_slope_uncensored_descriptive,weibull_plot_r2_uncensored_descriptive,mean_induction_uncensored_s,mean_final_active,active_fraction_susceptibility,mean_pit_count,mean_pit_density,mean_largest_cluster,mean_nn_pair_ratio,mean_dissolution_events,mean_repassivation_events")?;
    writeln!(
        induction,
        "potential_mu_kcal_mol,replica,seed,induction_time_s,censored_at_s"
    )?;
    writeln!(scaling, "size,potential_mu_kcal_mol,replicas,mean_active_fraction,active_fraction_susceptibility,mean_pit_density,mean_largest_cluster,mean_largest_cluster_fraction,mean_nn_pair_ratio")?;

    for potential in POTENTIALS {
        let deck = compiled_at(&text, potential)?;
        let mut results = Vec::with_capacity(REPLICAS as usize);
        for replica in 0..REPLICAS {
            let seed = petra_deck::replica_seed(deck.seed, replica, deck.seed_policy);
            results.push(run_replica(&deck, seed, None)?);
        }
        let initiated: Vec<_> = results
            .iter()
            .filter_map(|result| result.induction_time)
            .collect();
        for (replica, result) in results.iter().enumerate() {
            match result.induction_time {
                Some(time) => writeln!(
                    induction,
                    "{potential:.2},{replica},{},{time:.9},",
                    result.seed
                )?,
                None => writeln!(
                    induction,
                    "{potential:.2},{replica},{},,{HORIZON:.1}",
                    result.seed
                )?,
            }
        }
        let weibull = weibull_plot(&initiated);
        let mean_induction = (!initiated.is_empty()).then(|| mean(initiated.iter().copied()));
        let active_fractions: Vec<_> = results
            .iter()
            .map(|result| result.final_active as f64 / 576.0)
            .collect();
        let susceptibility = 576.0 * sample_variance(&active_fractions);
        writeln!(
            summary,
            "{potential:.2},{REPLICAS},{},{:.6},{},{},{},{:.6},{susceptibility:.6},{:.6},{:.9},{:.6},{:.6},{:.6},{:.6}",
            initiated.len(),
            initiated.len() as f64 / REPLICAS as f64,
            weibull.map_or_else(String::new, |value| format!("{:.6}", value.0)),
            weibull.map_or_else(String::new, |value| format!("{:.6}", value.1)),
            mean_induction.map_or_else(String::new, |value| format!("{value:.6}")),
            mean(results.iter().map(|result| result.final_active as f64)),
            mean(results.iter().map(|result| result.pit_count as f64)),
            mean(results.iter().map(|result| result.pit_count as f64)) / 576.0,
            mean(results.iter().map(|result| result.largest_cluster as f64)),
            mean(results.iter().map(|result| result.nn_pair_ratio)),
            mean(
                results
                    .iter()
                    .map(|result| result.dissolution_events as f64)
            ),
            mean(
                results
                    .iter()
                    .map(|result| result.repassivation_events as f64)
            ),
        )?;
        println!(
            "mu={potential:+.2}: initiated={}/{REPLICAS}, k={}, largest={:.2}",
            initiated.len(),
            weibull.map_or_else(|| "n/a".to_string(), |value| format!("{:.2}", value.0)),
            mean(results.iter().map(|result| result.largest_cluster as f64)),
        );
    }

    // A bounded three-size scan distinguishes occupancy from finite-size
    // clustering and locates the susceptibility peak without claiming a
    // critical exponent from this coarse grid.
    for size in SCALING_SIZES {
        let sites = (size * size) as f64;
        for potential in SCALING_POTENTIALS {
            let deck = compiled_at_size(&text, potential, size)?;
            let mut results = Vec::with_capacity(SCALING_REPLICAS as usize);
            for replica in 0..SCALING_REPLICAS {
                let seed = petra_deck::replica_seed(deck.seed, replica, deck.seed_policy);
                results.push(run_replica(&deck, seed, None)?);
            }
            let active_fractions: Vec<_> = results
                .iter()
                .map(|result| result.final_active as f64 / sites)
                .collect();
            writeln!(
                scaling,
                "{size},{potential:.2},{SCALING_REPLICAS},{:.9},{:.9},{:.9},{:.6},{:.9},{:.6}",
                mean(active_fractions.iter().copied()),
                sites * sample_variance(&active_fractions),
                mean(results.iter().map(|result| result.pit_count as f64)) / sites,
                mean(results.iter().map(|result| result.largest_cluster as f64)),
                mean(results.iter().map(|result| result.largest_cluster as f64)) / sites,
                mean(results.iter().map(|result| result.nn_pair_ratio)),
            )?;
        }
    }

    // A mid-transition representative event-rate series records both the
    // anodic-current propensity and repassivation competition.
    let transient_deck = compiled_at(&text, 0.30)?;
    let transient_path = output_dir.join("metastable-transient.csv");
    let mut transient = std::io::BufWriter::new(std::fs::File::create(&transient_path)?);
    let seed = petra_deck::replica_seed(transient_deck.seed, 3, transient_deck.seed_policy);
    run_replica(&transient_deck, seed, Some(&mut transient))?;
    transient.flush()?;
    summary.flush()?;
    induction.flush()?;
    scaling.flush()?;
    println!("wrote {}", display(&output_dir));
    Ok(())
}

fn display(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}
