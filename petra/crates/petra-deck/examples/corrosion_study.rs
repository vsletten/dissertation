use std::io::Write;
use std::path::{Path, PathBuf};

use petra_core::{CtmcAdvance, StateId};
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64Mcg;

const POTENTIALS: [f64; 11] = [
    -0.20, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85,
];
const REPLICAS: u64 = 256;
const SCALING_SIZES: [usize; 4] = [16, 24, 32, 40];
const SCALING_POTENTIALS: [f64; 9] = [-0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50];
const SCALING_REPLICAS: u64 = 256;
const BOOTSTRAP_RESAMPLES: usize = 2_000;
const HORIZON: f64 = 600.0;
const MAX_EVENTS: u64 = 120_000;
const STABLE_CLUSTER: u64 = 6;

#[derive(Debug)]
struct ReplicaResult {
    seed: u64,
    induction_time: Option<f64>,
    final_active: u64,
    pit_count: u64,
    largest_cluster: u64,
    nn_pair_ratio: Option<f64>,
    dissolution_events: u64,
    repassivation_events: u64,
    pits: Vec<PitProfile>,
}

#[derive(Debug, Clone, PartialEq)]
struct PitProfile {
    size: u64,
    total_depth: u64,
    max_depth: u64,
    mean_depth: f64,
}

#[derive(Debug, Clone, Copy)]
struct InductionObservation {
    time: f64,
    initiated: bool,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct WeibullFit {
    shape: f64,
    scale: f64,
    log_likelihood: f64,
    aic: f64,
    ks_distance: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct WeibullBootstrap {
    fit: WeibullFit,
    shape_ci: (f64, f64),
    scale_ci: (f64, f64),
    successful_resamples: usize,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct SurvivalPoint {
    time: f64,
    at_risk: usize,
    events: usize,
    censored: usize,
    survival: f64,
    cumulative_hazard: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct CrossoverFit {
    infinite_size: f64,
    slope: f64,
    r_squared: f64,
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

fn pit_profiles(
    engine: &petra_core::Engine,
    active: &[StateId],
    dissolution_depth: &[u64],
) -> Vec<PitProfile> {
    let mut seen = vec![false; engine.lattice.len()];
    let mut stack = Vec::new();
    let mut profiles = Vec::new();
    for start in 0..engine.lattice.len() {
        if seen[start] || !is_selected(engine.lattice.states[start], active) {
            continue;
        }
        seen[start] = true;
        stack.push(start);
        let mut sites = Vec::new();
        while let Some(site) = stack.pop() {
            sites.push(site);
            for &neighbor in engine.lattice.neighbors(site) {
                let neighbor = neighbor as usize;
                if !seen[neighbor] && is_selected(engine.lattice.states[neighbor], active) {
                    seen[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
        let total_depth = sites.iter().map(|&site| dissolution_depth[site]).sum();
        let max_depth = sites
            .iter()
            .map(|&site| dissolution_depth[site])
            .max()
            .unwrap_or(0);
        profiles.push(PitProfile {
            size: sites.len() as u64,
            total_depth,
            max_depth,
            mean_depth: total_depth as f64 / sites.len() as f64,
        });
    }
    profiles
}

// Nearest-neighbor pair correlation relative to a random surface with the
// same active fraction. Values above one mean pit patches cluster spatially.
fn fixed_occupancy_pair_probability(active_sites: u64, total_sites: usize) -> Option<f64> {
    if active_sites < 2 || total_sites < 2 || active_sites as usize > total_sites {
        return None;
    }
    Some(
        active_sites as f64 * (active_sites - 1) as f64
            / (total_sites as f64 * (total_sites - 1) as f64),
    )
}

fn nn_pair_ratio(engine: &petra_core::Engine, active: &[StateId]) -> Option<f64> {
    let active_sites = active_count(engine, active);
    let random_pair_fraction =
        fixed_occupancy_pair_probability(active_sites, engine.lattice.len())?;
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
    (total_edges > 0).then(|| (active_edges as f64 / total_edges as f64) / random_pair_fraction)
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
    let mut dissolution_depth = vec![0_u64; engine.lattice.len()];

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
            dissolution_depth[fired.site as usize] += 1;
        } else if name.starts_with("repassivate_") {
            repassivation_events += 1;
        }

        // Record first formation immediately; keep sparse sampling only for CSV.
        let sample = event.is_multiple_of(20) || name == "film_rupture";
        if sample || induction_time.is_none() {
            let clusters = cluster_stats(&engine, &pit);
            if induction_time.is_none() && clusters.largest >= STABLE_CLUSTER {
                induction_time = Some(engine.time);
            }
            if sample {
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
    }

    if engine.time < HORIZON {
        return Err(format!(
            "seed {seed}: bounded event cap reached at t={:.3} < {HORIZON}",
            engine.time
        ));
    }
    debug_assert_eq!(engine.time.to_bits(), HORIZON.to_bits());
    let clusters = cluster_stats(&engine, &pit);
    let pits = pit_profiles(&engine, &pit, &dissolution_depth);
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
        pits,
    })
}

fn survival_diagnostics(observations: &[InductionObservation]) -> Vec<SurvivalPoint> {
    let mut ordered = observations.to_vec();
    ordered.sort_by(|left, right| left.time.total_cmp(&right.time));
    let mut at_risk = ordered.len();
    let mut survival = 1.0;
    let mut cumulative_hazard = 0.0;
    let mut points = Vec::new();
    let mut cursor = 0;
    while cursor < ordered.len() {
        let time = ordered[cursor].time;
        let mut events = 0;
        let mut censored = 0;
        while cursor < ordered.len() && ordered[cursor].time == time {
            if ordered[cursor].initiated {
                events += 1;
            } else {
                censored += 1;
            }
            cursor += 1;
        }
        if events > 0 {
            survival *= 1.0 - events as f64 / at_risk as f64;
            cumulative_hazard += events as f64 / at_risk as f64;
        }
        points.push(SurvivalPoint {
            time,
            at_risk,
            events,
            censored,
            survival,
            cumulative_hazard,
        });
        at_risk -= events + censored;
    }
    points
}

fn fit_weibull_censored(observations: &[InductionObservation]) -> Option<WeibullFit> {
    let event_logs: Vec<_> = observations
        .iter()
        .filter(|observation| observation.initiated)
        .map(|observation| observation.time.ln())
        .collect();
    if event_logs.len() < 2
        || observations.is_empty()
        || observations
            .iter()
            .any(|observation| !observation.time.is_finite() || observation.time <= 0.0)
    {
        return None;
    }
    let logs: Vec<_> = observations
        .iter()
        .map(|observation| observation.time.ln())
        .collect();
    let event_log_mean = event_logs.iter().sum::<f64>() / event_logs.len() as f64;
    let max_log = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let score = |shape: f64| {
        let weights: Vec<_> = logs
            .iter()
            .map(|log_time| (shape * (log_time - max_log)).exp())
            .collect();
        let total = weights.iter().sum::<f64>();
        let weighted_log = weights
            .iter()
            .zip(&logs)
            .map(|(weight, log_time)| weight * log_time)
            .sum::<f64>()
            / total;
        1.0 / shape + event_log_mean - weighted_log
    };

    let mut low = 0.05;
    let mut high = 20.0;
    if score(low) <= 0.0 {
        return None;
    }
    while score(high) > 0.0 && high < 1_000.0 {
        high *= 2.0;
    }
    if score(high) > 0.0 {
        return None;
    }
    for _ in 0..100 {
        let middle = (low + high) / 2.0;
        if score(middle) > 0.0 {
            low = middle;
        } else {
            high = middle;
        }
    }
    let shape = (low + high) / 2.0;
    let scaled_sum = logs
        .iter()
        .map(|log_time| (shape * (log_time - max_log)).exp())
        .sum::<f64>();
    let log_scale = max_log + (scaled_sum / event_logs.len() as f64).ln() / shape;
    let scale = log_scale.exp();
    let normalized_sum = logs
        .iter()
        .map(|log_time| (shape * (log_time - log_scale)).exp())
        .sum::<f64>();
    let log_likelihood = event_logs.len() as f64 * shape.ln()
        - event_logs.len() as f64 * shape * log_scale
        + (shape - 1.0) * event_logs.iter().sum::<f64>()
        - normalized_sum;
    let ks_distance = survival_diagnostics(observations)
        .iter()
        .filter(|point| point.events > 0)
        .map(|point| {
            let fitted_survival = (-(point.time / scale).powf(shape)).exp();
            (point.survival - fitted_survival).abs()
        })
        .fold(0.0_f64, f64::max);
    Some(WeibullFit {
        shape,
        scale,
        log_likelihood,
        aic: 4.0 - 2.0 * log_likelihood,
        ks_distance,
    })
}

fn percentile(sorted: &[f64], probability: f64) -> f64 {
    let index = ((sorted.len() - 1) as f64 * probability).round() as usize;
    sorted[index]
}

fn bootstrap_weibull(
    observations: &[InductionObservation],
    resamples: usize,
    seed: u64,
) -> Option<WeibullBootstrap> {
    let fit = fit_weibull_censored(observations)?;
    if observations.is_empty() || resamples == 0 {
        return None;
    }
    let mut rng = Pcg64Mcg::seed_from_u64(seed);
    let mut shapes = Vec::with_capacity(resamples);
    let mut scales = Vec::with_capacity(resamples);
    for _ in 0..resamples {
        let sample: Vec<_> = (0..observations.len())
            .map(|_| observations[rng.gen_range(0..observations.len())])
            .collect();
        if let Some(sample_fit) = fit_weibull_censored(&sample) {
            shapes.push(sample_fit.shape);
            scales.push(sample_fit.scale);
        }
    }
    if shapes.is_empty() {
        return None;
    }
    shapes.sort_by(f64::total_cmp);
    scales.sort_by(f64::total_cmp);
    Some(WeibullBootstrap {
        fit,
        shape_ci: (percentile(&shapes, 0.025), percentile(&shapes, 0.975)),
        scale_ci: (percentile(&scales, 0.025), percentile(&scales, 0.975)),
        successful_resamples: shapes.len(),
    })
}

fn interpolate_crossing(samples: &[(f64, f64)], target: f64) -> Result<f64, String> {
    if samples.len() < 2 || !(0.0..=1.0).contains(&target) {
        return Err("crossing requires two samples and a target in [0, 1]".to_string());
    }
    for pair in samples.windows(2) {
        if !pair[0].0.is_finite()
            || !pair[0].1.is_finite()
            || !(0.0..=1.0).contains(&pair[0].1)
            || pair[0].0 >= pair[1].0
        {
            return Err(
                "potentials must increase and fractions must be finite in [0, 1]".to_string(),
            );
        }
    }
    let &(last_potential, last_fraction) = samples.last().expect("two samples");
    if !last_potential.is_finite()
        || !last_fraction.is_finite()
        || !(0.0..=1.0).contains(&last_fraction)
    {
        return Err("potentials and fractions must be finite".to_string());
    }

    // Equal-weight pool-adjacent-violators fit keeps finite-replica noise from
    // creating a false reverse crossing while preserving the deterministic data.
    let mut blocks: Vec<(usize, usize, f64, usize)> = Vec::new();
    for (index, &(_, fraction)) in samples.iter().enumerate() {
        blocks.push((index, index + 1, fraction, 1));
        while blocks.len() >= 2 {
            let right = blocks[blocks.len() - 1];
            let left = blocks[blocks.len() - 2];
            if left.2 / left.3 as f64 <= right.2 / right.3 as f64 {
                break;
            }
            blocks.pop();
            blocks.pop();
            blocks.push((left.0, right.1, left.2 + right.2, left.3 + right.3));
        }
    }
    let mut fitted = vec![0.0; samples.len()];
    for (start, end, sum, count) in blocks {
        fitted[start..end].fill(sum / count as f64);
    }
    if fitted[0] > target || fitted[fitted.len() - 1] < target {
        return Err(format!(
            "target fraction {target:.3} is not bracketed by the potential grid"
        ));
    }
    for index in 0..fitted.len() - 1 {
        if fitted[index] <= target && target <= fitted[index + 1] {
            let low = samples[index].0;
            let high = samples[index + 1].0;
            let span = fitted[index + 1] - fitted[index];
            return Ok(if span == 0.0 {
                (low + high) / 2.0
            } else {
                low + (target - fitted[index]) / span * (high - low)
            });
        }
    }
    Err("no crossing found after isotonic fit".to_string())
}

fn extrapolate_crossover(peaks: &[(usize, f64)]) -> Option<CrossoverFit> {
    if peaks.len() < 3 {
        return None;
    }
    let points: Vec<_> = peaks
        .iter()
        .map(|&(size, potential)| (1.0 / size as f64, potential))
        .collect();
    let count = points.len() as f64;
    let x_mean = points.iter().map(|point| point.0).sum::<f64>() / count;
    let y_mean = points.iter().map(|point| point.1).sum::<f64>() / count;
    let denominator = points
        .iter()
        .map(|(x, _)| (x - x_mean).powi(2))
        .sum::<f64>();
    if denominator == 0.0 {
        return None;
    }
    let slope = points
        .iter()
        .map(|(x, y)| (x - x_mean) * (y - y_mean))
        .sum::<f64>()
        / denominator;
    let infinite_size = y_mean - slope * x_mean;
    let residual = points
        .iter()
        .map(|(x, y)| (y - (infinite_size + slope * x)).powi(2))
        .sum::<f64>();
    let total = points
        .iter()
        .map(|(_, y)| (y - y_mean).powi(2))
        .sum::<f64>();
    Some(CrossoverFit {
        infinite_size,
        slope,
        r_squared: if total == 0.0 {
            1.0
        } else {
            1.0 - residual / total
        },
    })
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
    let mut survival = std::io::BufWriter::new(std::fs::File::create(
        output_dir.join("survival-hazard.csv"),
    )?);
    let mut pits = std::io::BufWriter::new(std::fs::File::create(
        output_dir.join("pit-size-depth.csv"),
    )?);
    let mut scaling = std::io::BufWriter::new(std::fs::File::create(
        output_dir.join("critical-scaling.csv"),
    )?);
    writeln!(summary, "potential_mu_kcal_mol,replicas,initiated,initiation_fraction,weibull_shape_mle,weibull_shape_ci95_low,weibull_shape_ci95_high,weibull_scale_s_mle,weibull_scale_ci95_low_s,weibull_scale_ci95_high_s,weibull_ks_distance,weibull_log_likelihood,weibull_aic,bootstrap_successful_resamples,nucleation_rate_per_patch_s,mean_final_active,active_fraction_susceptibility,mean_pit_count,mean_pit_density,mean_largest_cluster,p95_largest_cluster,mean_nn_pair_ratio,nn_pair_ratio_defined_replicas,mean_dissolution_events,mean_repassivation_events")?;
    writeln!(
        induction,
        "potential_mu_kcal_mol,replica,seed,observed_time_s,initiated,censored_at_s"
    )?;
    writeln!(survival, "potential_mu_kcal_mol,time_s,at_risk,events,censored,kaplan_meier_survival,nelson_aalen_cumulative_hazard,fitted_weibull_survival,fitted_weibull_cumulative_hazard")?;
    writeln!(pits, "potential_mu_kcal_mol,replica,seed,pit_index,final_cluster_size_patches,total_dissolution_depth_events,max_dissolution_depth_events,mean_dissolution_depth_events,replica_nn_pair_ratio")?;
    writeln!(scaling, "size,potential_mu_kcal_mol,replicas,initiated,initiation_fraction,mean_active_fraction,active_fraction_susceptibility,mean_pit_density,mean_largest_cluster,mean_largest_cluster_fraction,mean_nn_pair_ratio,nn_pair_ratio_defined_replicas")?;

    for potential in POTENTIALS {
        let deck = compiled_at(&text, potential)?;
        let sites = engine_sites(&deck)? as f64;
        let mut results = Vec::with_capacity(REPLICAS as usize);
        for replica in 0..REPLICAS {
            let seed = petra_deck::replica_seed(deck.seed, replica, deck.seed_policy);
            results.push(run_replica(&deck, seed, None)?);
        }
        let observations: Vec<_> = results
            .iter()
            .map(|result| InductionObservation {
                time: result.induction_time.unwrap_or(HORIZON),
                initiated: result.induction_time.is_some(),
            })
            .collect();
        let initiated = observations
            .iter()
            .filter(|observation| observation.initiated)
            .count();
        let weibull = bootstrap_weibull(
            &observations,
            BOOTSTRAP_RESAMPLES,
            deck.seed ^ potential.to_bits() ^ 0xC3C3_2026,
        );
        for (replica, (result, observation)) in results.iter().zip(observations.iter()).enumerate()
        {
            writeln!(
                induction,
                "{potential:.2},{replica},{},{:.9},{},{}",
                result.seed,
                observation.time,
                observation.initiated,
                if observation.initiated {
                    String::new()
                } else {
                    format!("{HORIZON:.1}")
                }
            )?;
            for (pit_index, pit) in result.pits.iter().enumerate() {
                writeln!(
                    pits,
                    "{potential:.2},{replica},{},{pit_index},{},{},{},{:.9},{}",
                    result.seed,
                    pit.size,
                    pit.total_depth,
                    pit.max_depth,
                    pit.mean_depth,
                    result
                        .nn_pair_ratio
                        .map_or_else(String::new, |value| format!("{value:.9}")),
                )?;
            }
        }
        for point in survival_diagnostics(&observations) {
            let fitted_survival =
                weibull.map(|value| (-(point.time / value.fit.scale).powf(value.fit.shape)).exp());
            let fitted_hazard =
                weibull.map(|value| (point.time / value.fit.scale).powf(value.fit.shape));
            writeln!(
                survival,
                "{potential:.2},{:.9},{},{},{},{:.9},{:.9},{},{}",
                point.time,
                point.at_risk,
                point.events,
                point.censored,
                point.survival,
                point.cumulative_hazard,
                fitted_survival.map_or_else(String::new, |value| format!("{value:.9}")),
                fitted_hazard.map_or_else(String::new, |value| format!("{value:.9}")),
            )?;
        }
        let active_fractions: Vec<_> = results
            .iter()
            .map(|result| result.final_active as f64 / sites)
            .collect();
        let susceptibility = sites * sample_variance(&active_fractions);
        let mut largest: Vec<_> = results
            .iter()
            .map(|result| result.largest_cluster as f64)
            .collect();
        largest.sort_by(f64::total_cmp);
        let total_time_at_risk = observations
            .iter()
            .map(|observation| observation.time)
            .sum::<f64>();
        let nn_pair_ratios: Vec<_> = results
            .iter()
            .filter_map(|result| result.nn_pair_ratio)
            .collect();
        let fit_fields = weibull.map(|value| {
            vec![
                format!("{:.9}", value.fit.shape),
                format!("{:.9}", value.shape_ci.0),
                format!("{:.9}", value.shape_ci.1),
                format!("{:.9}", value.fit.scale),
                format!("{:.9}", value.scale_ci.0),
                format!("{:.9}", value.scale_ci.1),
                format!("{:.9}", value.fit.ks_distance),
                format!("{:.9}", value.fit.log_likelihood),
                format!("{:.9}", value.fit.aic),
                value.successful_resamples.to_string(),
            ]
        });
        let mut row = vec![
            format!("{potential:.2}"),
            REPLICAS.to_string(),
            initiated.to_string(),
            format!("{:.9}", initiated as f64 / REPLICAS as f64),
        ];
        row.extend(fit_fields.unwrap_or_else(|| vec![String::new(); 10]));
        row.extend([
            format!("{:.12}", initiated as f64 / total_time_at_risk / sites),
            format!(
                "{:.9}",
                mean(results.iter().map(|result| result.final_active as f64))
            ),
            format!("{susceptibility:.9}"),
            format!(
                "{:.9}",
                mean(results.iter().map(|result| result.pit_count as f64))
            ),
            format!(
                "{:.12}",
                mean(results.iter().map(|result| result.pit_count as f64)) / sites
            ),
            format!(
                "{:.9}",
                mean(results.iter().map(|result| result.largest_cluster as f64))
            ),
            format!("{:.9}", percentile(&largest, 0.95)),
            format!("{:.9}", mean(nn_pair_ratios.iter().copied())),
            nn_pair_ratios.len().to_string(),
            format!(
                "{:.9}",
                mean(
                    results
                        .iter()
                        .map(|result| result.dissolution_events as f64)
                )
            ),
            format!(
                "{:.9}",
                mean(
                    results
                        .iter()
                        .map(|result| result.repassivation_events as f64)
                )
            ),
        ]);
        writeln!(summary, "{}", row.join(","))?;
        println!(
            "mu={potential:+.2}: initiated={initiated}/{REPLICAS}, k={}, largest={:.2}",
            weibull.map_or_else(
                || "n/a".to_string(),
                |value| format!("{:.2}", value.fit.shape)
            ),
            mean(results.iter().map(|result| result.largest_cluster as f64)),
        );
    }

    let mut peaks = Vec::new();
    for size in SCALING_SIZES {
        let sites = (size * size) as f64;
        let mut size_samples = Vec::new();
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
            let susceptibility = sites * sample_variance(&active_fractions);
            let initiated = results
                .iter()
                .filter(|result| result.induction_time.is_some())
                .count();
            let initiation_fraction = initiated as f64 / SCALING_REPLICAS as f64;
            let nn_pair_ratios: Vec<_> = results
                .iter()
                .filter_map(|result| result.nn_pair_ratio)
                .collect();
            size_samples.push((potential, initiation_fraction));
            writeln!(
                scaling,
                "{size},{potential:.2},{SCALING_REPLICAS},{initiated},{initiation_fraction:.9},{:.9},{susceptibility:.9},{:.12},{:.9},{:.12},{:.9},{}",
                mean(active_fractions.iter().copied()),
                mean(results.iter().map(|result| result.pit_count as f64)) / sites,
                mean(results.iter().map(|result| result.largest_cluster as f64)),
                mean(results.iter().map(|result| result.largest_cluster as f64)) / sites,
                mean(nn_pair_ratios.iter().copied()),
                nn_pair_ratios.len(),
            )?;
        }
        let crossing_potential = interpolate_crossing(&size_samples, 0.5)?;
        peaks.push((size, crossing_potential));
    }
    let crossover = extrapolate_crossover(&peaks).expect("four finite sizes identify a fit");
    let mut crossover_file =
        std::io::BufWriter::new(std::fs::File::create(output_dir.join("crossover-fit.csv"))?);
    writeln!(
        crossover_file,
        "infinite_size_half_initiation_potential_mu_kcal_mol,inverse_size_slope,r_squared,sizes,half_initiation_potentials"
    )?;
    writeln!(
        crossover_file,
        "{:.9},{:.9},{:.9},\"{}\",\"{}\"",
        crossover.infinite_size,
        crossover.slope,
        crossover.r_squared,
        peaks
            .iter()
            .map(|(size, _)| size.to_string())
            .collect::<Vec<_>>()
            .join(";"),
        peaks
            .iter()
            .map(|(_, potential)| format!("{potential:.9}"))
            .collect::<Vec<_>>()
            .join(";"),
    )?;

    let transient_deck = compiled_at(&text, 0.30)?;
    let transient_path = output_dir.join("metastable-transient.csv");
    let mut transient = std::io::BufWriter::new(std::fs::File::create(&transient_path)?);
    let seed = petra_deck::replica_seed(transient_deck.seed, 3, transient_deck.seed_policy);
    run_replica(&transient_deck, seed, Some(&mut transient))?;
    transient.flush()?;
    summary.flush()?;
    induction.flush()?;
    survival.flush()?;
    pits.flush()?;
    scaling.flush()?;
    crossover_file.flush()?;
    println!("wrote {}", display(&output_dir));
    Ok(())
}

fn engine_sites(deck: &petra_deck::CompiledDeck) -> Result<usize, String> {
    deck.build_engine(Some(deck.seed))
        .map(|engine| engine.lattice.len())
        .map_err(|error| error.to_string())
}

fn display(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn censored_weibull_mle_recovers_a_known_quantile_population() {
        let shape = 2.0;
        let scale = 10.0;
        let horizon = 12.0;
        let observations: Vec<_> = (0..400)
            .map(|index| {
                let probability = (index as f64 + 0.5) / 400.0;
                let time = scale * (-(1.0 - probability).ln()).powf(1.0 / shape);
                InductionObservation {
                    time: time.min(horizon),
                    initiated: time <= horizon,
                }
            })
            .collect();
        let fit = fit_weibull_censored(&observations).expect("events identify a fit");
        assert!((fit.shape - shape).abs() < 0.03, "shape={}", fit.shape);
        assert!((fit.scale - scale).abs() < 0.05, "scale={}", fit.scale);
        assert!(fit.ks_distance < 0.01, "ks={}", fit.ks_distance);
    }

    #[test]
    fn survival_diagnostics_handle_events_and_censoring_at_the_same_time() {
        let observations = vec![
            InductionObservation {
                time: 1.0,
                initiated: true,
            },
            InductionObservation {
                time: 2.0,
                initiated: true,
            },
            InductionObservation {
                time: 2.0,
                initiated: false,
            },
        ];
        let points = survival_diagnostics(&observations);
        assert_eq!(points.len(), 2);
        assert_eq!(points[0].at_risk, 3);
        assert_eq!(points[0].events, 1);
        assert!((points[0].survival - 2.0 / 3.0).abs() < 1e-12);
        assert_eq!(points[1].at_risk, 2);
        assert_eq!(points[1].events, 1);
        assert_eq!(points[1].censored, 1);
        assert!((points[1].survival - 1.0 / 3.0).abs() < 1e-12);
        assert!((points[1].cumulative_hazard - 5.0 / 6.0).abs() < 1e-12);
    }

    #[test]
    fn bootstrap_weibull_interval_is_deterministic_and_contains_the_estimate() {
        let observations: Vec<_> = (1..=40)
            .map(|value| InductionObservation {
                time: value as f64,
                initiated: value <= 30,
            })
            .collect();
        let first = bootstrap_weibull(&observations, 200, 73).expect("bootstrap fit");
        let repeat = bootstrap_weibull(&observations, 200, 73).expect("bootstrap fit");
        assert_eq!(first, repeat);
        assert!(first.shape_ci.0 <= first.fit.shape && first.fit.shape <= first.shape_ci.1);
        assert!(first.scale_ci.0 <= first.fit.scale && first.fit.scale <= first.scale_ci.1);
    }

    #[test]
    fn inverse_size_extrapolation_recovers_infinite_size_crossover() {
        let peaks = [(16_usize, 0.50), (24, 0.40), (32, 0.35), (40, 0.32)];
        let fit = extrapolate_crossover(&peaks).expect("four sizes fit");
        assert!((fit.infinite_size - 0.20).abs() < 1e-12);
        assert!((fit.slope - 4.80).abs() < 1e-12);
        assert!((fit.r_squared - 1.0).abs() < 1e-12);
    }

    #[test]
    fn finite_size_crossover_interpolates_the_half_initiation_point() {
        let samples = [(-0.1, 0.2), (0.1, 0.4), (0.3, 0.8)];
        let crossing = interpolate_crossing(&samples, 0.5).expect("bracketed crossing");
        assert!((crossing - 0.15).abs() < 1e-12);
        assert!(interpolate_crossing(&samples, 0.9).is_err());
    }

    #[test]
    fn random_pair_baseline_conditions_on_the_observed_active_count() {
        assert!((fixed_occupancy_pair_probability(2, 4).unwrap() - 1.0 / 6.0).abs() < 1e-12);
        assert!(fixed_occupancy_pair_probability(1, 4).is_none());
        assert!(fixed_occupancy_pair_probability(4, 4).is_some_and(|value| value == 1.0));
    }
}
