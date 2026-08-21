//! Declarative observables and deterministic replica-level parallelism.
//!
//! Rayon is used only across replicas. Each trajectory remains serial, and
//! indexed parallel collection restores replica-index order independent of
//! scheduling.

use petra_core::{Engine, StateId, Stop};
use petra_deck::{replica_seed, CompiledDeck, CompiledObservable};
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64Mcg;
use rayon::prelude::*;

#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceArea {
    /// Selected solid sites exposed to a declared surface or a non-solid neighbor.
    pub exposed_sites: u64,
    /// Projected geometric area in the deck cell's squared length unit.
    pub geometric: f64,
    /// BET-like site-count proxy. Kept separately so later adsorbate weighting
    /// can evolve without changing the geometric measure.
    pub bet_site_proxy: u64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ObservableValue {
    StateCounts(Vec<u64>),
    EventRates(Vec<f64>),
    RateSpectrum(Vec<f64>),
    ClusterSizes(Vec<u64>),
    SurfaceArea(SurfaceArea),
}

#[derive(Debug, Clone, PartialEq)]
pub struct Sample {
    pub step: u64,
    pub time: f64,
    /// Same order as the deck's `[[observables.series]]` declarations.
    pub values: Vec<ObservableValue>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct DistributionSummary {
    pub values: Vec<f64>,
    pub mean: f64,
    pub ci95: (f64, f64),
}

/// Deterministic percentile-bootstrap confidence interval for the mean.
pub fn summarize(values: &[f64], bootstrap_resamples: usize, seed: u64) -> DistributionSummary {
    assert!(!values.is_empty(), "cannot summarize an empty distribution");
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    if bootstrap_resamples == 0 {
        return DistributionSummary {
            values: values.to_vec(),
            mean,
            ci95: (mean, mean),
        };
    }
    let mut rng = Pcg64Mcg::seed_from_u64(seed);
    let mut bootstrap = Vec::with_capacity(bootstrap_resamples);
    for _ in 0..bootstrap_resamples {
        let mut total = 0.0;
        for _ in 0..values.len() {
            total += values[rng.gen_range(0..values.len())];
        }
        bootstrap.push(total / values.len() as f64);
    }
    bootstrap.sort_by(f64::total_cmp);
    let lower = ((bootstrap_resamples - 1) as f64 * 0.025).round() as usize;
    let upper = ((bootstrap_resamples - 1) as f64 * 0.975).round() as usize;
    DistributionSummary {
        values: values.to_vec(),
        mean,
        ci95: (bootstrap[lower], bootstrap[upper]),
    }
}

fn selected(state: StateId, states: &[StateId]) -> bool {
    states.contains(&state)
}

fn cluster_sizes(engine: &Engine, states: &[StateId]) -> Vec<u64> {
    let lattice = &engine.lattice;
    let mut seen = vec![false; lattice.len()];
    let mut sizes = Vec::new();
    let mut stack = Vec::new();
    for start in 0..lattice.len() {
        if seen[start] || !selected(lattice.states[start], states) {
            continue;
        }
        seen[start] = true;
        stack.push(start);
        let mut size = 0u64;
        while let Some(site) = stack.pop() {
            size += 1;
            // Reverse insertion makes the lowest site id pop first. Cluster
            // sizes are order-independent, but deterministic traversal is a
            // maintained simulation/reporting invariant.
            for &neighbor in lattice.neighbors(site).iter().rev() {
                let neighbor = neighbor as usize;
                if !seen[neighbor] && selected(lattice.states[neighbor], states) {
                    seen[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
        sizes.push(size);
    }
    sizes.sort_unstable();
    sizes
}

fn face_area(deck: &CompiledDeck, axis: usize) -> f64 {
    let matrix = deck.unit_cell.cell.matrix();
    let other: Vec<_> = (0..3).filter(|&candidate| candidate != axis).collect();
    let a = [
        matrix[0][other[0]],
        matrix[1][other[0]],
        matrix[2][other[0]],
    ];
    let b = [
        matrix[0][other[1]],
        matrix[1][other[1]],
        matrix[2][other[1]],
    ];
    let cross = [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ];
    (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt()
        / deck.unit_cell.sites.len().max(1) as f64
}

fn surface_area(engine: &Engine, deck: &CompiledDeck, states: &[StateId], axis: u8) -> SurfaceArea {
    let lattice = &engine.lattice;
    let axis = axis as usize;
    let mut exposed_sites = 0u64;
    for site in 0..lattice.len() {
        if !selected(lattice.states[site], states) {
            continue;
        }
        let (cell, _) = lattice.coords(site);
        let declared_face = lattice.boundary[axis] == petra_core::Boundary::Open
            && (cell[axis] == 0 || cell[axis] + 1 == lattice.dims[axis]);
        let internal_face = lattice
            .neighbors(site)
            .iter()
            .any(|&neighbor| !selected(lattice.states[neighbor as usize], states));
        if declared_face || internal_face {
            exposed_sites += 1;
        }
    }
    SurfaceArea {
        exposed_sites,
        geometric: exposed_sites as f64 * face_area(deck, axis),
        bet_site_proxy: exposed_sites,
    }
}

/// Evaluate every declared observable without mutating trajectory state.
pub fn observe(engine: &Engine, deck: &CompiledDeck) -> Sample {
    let values = deck
        .observables
        .iter()
        .filter_map(|observable| match observable {
            CompiledObservable::StateCounts => Some(ObservableValue::StateCounts(
                engine.state_counts(deck.n_states),
            )),
            CompiledObservable::EventRates => {
                let mut rates = vec![0.0; engine.reactions.len()];
                for events in engine.site_event_rates() {
                    for &(rule, rate) in events {
                        rates[rule as usize] += rate;
                    }
                }
                Some(ObservableValue::EventRates(rates))
            }
            CompiledObservable::RateSpectra => Some(ObservableValue::RateSpectrum(
                engine
                    .site_event_rates()
                    .iter()
                    .map(|events| events.iter().map(|(_, rate)| rate).sum::<f64>())
                    .filter(|rate| *rate > 0.0)
                    .collect(),
            )),
            CompiledObservable::ClusterSizes { states } => {
                Some(ObservableValue::ClusterSizes(cluster_sizes(engine, states)))
            }
            CompiledObservable::SurfaceArea { states, axis } => Some(ObservableValue::SurfaceArea(
                surface_area(engine, deck, states, *axis),
            )),
            // Existing exporters own these two outputs. They remain accepted
            // schema declarations but do not produce scalar ensemble values.
            CompiledObservable::Snapshot | CompiledObservable::InterfaceRoughness { .. } => None,
        })
        .collect();
    Sample {
        step: engine.step_count,
        time: engine.time,
        values,
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EnsembleConfig {
    pub replicas: u64,
    pub base_seed: u64,
    pub steps: u64,
    pub burn_in: u64,
    pub sample_every: u64,
    pub bootstrap_resamples: usize,
    pub bootstrap_seed: u64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ReplicaRun {
    pub replica: u64,
    pub seed: u64,
    pub samples: Vec<Sample>,
    pub final_state_counts: Vec<u64>,
    pub stop: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct EnsembleResult {
    pub replicas: Vec<ReplicaRun>,
    pub final_state_counts: Vec<DistributionSummary>,
}

fn run_replica(
    deck: &CompiledDeck,
    config: &EnsembleConfig,
    replica: u64,
) -> Result<ReplicaRun, String> {
    let seed = replica_seed(config.base_seed, replica, deck.seed_policy);
    let mut engine = deck
        .build_engine(Some(seed))
        .map_err(|error| error.to_string())?;
    let mut strategy = deck.strategy();
    let mut stop: Option<Stop> = None;
    for _ in 0..config.burn_in {
        if let Err(reason) = engine.step_with(&mut strategy) {
            stop = Some(reason);
            break;
        }
    }
    let mut samples = vec![observe(&engine, deck)];
    if stop.is_none() {
        for step in 1..=config.steps {
            if let Err(reason) = engine.step_with(&mut strategy) {
                stop = Some(reason);
                break;
            }
            if step % config.sample_every == 0 {
                samples.push(observe(&engine, deck));
            }
        }
    }
    if samples
        .last()
        .is_none_or(|sample| sample.step != engine.step_count)
    {
        samples.push(observe(&engine, deck));
    }
    Ok(ReplicaRun {
        replica,
        seed,
        samples,
        final_state_counts: engine.state_counts(deck.n_states),
        stop: stop.map(|reason| reason.to_string()),
    })
}

/// Run independent trajectories using Rayon only across replicas.
pub fn run_ensemble(
    deck: &CompiledDeck,
    config: &EnsembleConfig,
) -> Result<EnsembleResult, String> {
    if config.replicas == 0 {
        return Err("ensemble requires at least one replica".to_string());
    }
    if config.sample_every == 0 {
        return Err("sample_every must be at least one".to_string());
    }
    let results: Vec<Result<ReplicaRun, String>> = (0..config.replicas)
        .into_par_iter()
        .map(|replica| run_replica(deck, config, replica))
        .collect();
    let replicas: Vec<ReplicaRun> = results.into_iter().collect::<Result<_, _>>()?;
    let final_state_counts = (0..deck.n_states)
        .map(|state| {
            let values: Vec<f64> = replicas
                .iter()
                .map(|replica| replica.final_state_counts[state] as f64)
                .collect();
            summarize(
                &values,
                config.bootstrap_resamples,
                config.bootstrap_seed.wrapping_add(state as u64),
            )
        })
        .collect();
    Ok(EnsembleResult {
        replicas,
        final_state_counts,
    })
}
