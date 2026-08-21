use petra_observables::{observe, run_ensemble, summarize, EnsembleConfig, ObservableValue};

fn ising_deck() -> petra_deck::CompiledDeck {
    let path = format!(
        "{}/../../examples/conformance/ising.toml",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(path).expect("Ising fixture");
    let parsed: petra_deck::DeckFile = toml::from_str(&text).expect("deck parses");
    petra_deck::compile(&parsed).expect("deck compiles")
}

#[test]
fn summary_keeps_distribution_and_has_reproducible_bootstrap_ci() {
    let first = summarize(&[1.0, 2.0, 3.0, 4.0], 2_000, 91);
    let second = summarize(&[1.0, 2.0, 3.0, 4.0], 2_000, 91);
    assert_eq!(first, second);
    assert_eq!(first.values, vec![1.0, 2.0, 3.0, 4.0]);
    assert_eq!(first.mean, 2.5);
    assert!(first.ci95.0 <= 2.5 && first.ci95.1 >= 2.5);
}

#[test]
fn ensemble_runs_replicas_in_parallel_but_returns_replica_order() {
    let deck = ising_deck();
    let config = EnsembleConfig {
        replicas: 4,
        base_seed: 99,
        steps: 20,
        burn_in: 0,
        sample_every: 5,
        bootstrap_resamples: 200,
        bootstrap_seed: 7,
    };
    let first = run_ensemble(&deck, &config).expect("ensemble runs");
    let second = run_ensemble(&deck, &config).expect("ensemble repeats");
    let expected: Vec<_> = (0..4)
        .map(|k| petra_deck::replica_seed(config.base_seed, k, deck.seed_policy))
        .collect();
    assert_eq!(
        first.replicas.iter().map(|r| r.seed).collect::<Vec<_>>(),
        expected
    );
    assert_eq!(first, second, "parallel scheduling must not affect results");
    assert!(first.replicas.iter().all(|r| r.samples.len() == 5));
    assert_eq!(first.final_state_counts.len(), deck.n_states);
}

#[test]
fn ensemble_always_records_the_final_state_between_cadence_boundaries() {
    let deck = ising_deck();
    let run = run_ensemble(
        &deck,
        &EnsembleConfig {
            replicas: 1,
            base_seed: deck.seed,
            steps: 7,
            burn_in: 0,
            sample_every: 5,
            bootstrap_resamples: 10,
            bootstrap_seed: 3,
        },
    )
    .expect("ensemble runs");
    assert_eq!(
        run.replicas[0]
            .samples
            .iter()
            .map(|sample| sample.step)
            .collect::<Vec<_>>(),
        vec![0, 5, 7]
    );
}

#[test]
fn declared_kossel_observables_cover_rates_clusters_and_area() {
    let path = format!(
        "{}/../../examples/kossel-etchpit.toml",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(path).expect("Kossel fixture");
    let parsed: petra_deck::DeckFile = toml::from_str(&text).expect("deck parses");
    let deck = petra_deck::compile(&parsed).expect("observable declarations compile");
    let engine = deck.build_engine(Some(42)).expect("engine builds");
    let sample = observe(&engine, &deck);
    assert!(matches!(sample.values[0], ObservableValue::StateCounts(_)));
    assert!(matches!(sample.values[1], ObservableValue::EventRates(_)));
    match &sample.values[2] {
        ObservableValue::RateSpectrum(rates) => assert!(!rates.is_empty()),
        other => panic!("expected rate spectrum, got {other:?}"),
    }
    assert_eq!(sample.values[3], ObservableValue::ClusterSizes(Vec::new()));
    match sample.values[4] {
        ObservableValue::SurfaceArea(ref area) => {
            assert!(area.exposed_sites > 0);
            assert!(area.geometric > 0.0);
            assert_eq!(area.bet_site_proxy, area.exposed_sites);
        }
        ref other => panic!("expected surface area, got {other:?}"),
    }
}

fn spectrum_log10_width(sample: &petra_observables::Sample) -> f64 {
    let rates = sample
        .values
        .iter()
        .find_map(|value| match value {
            ObservableValue::RateSpectrum(rates) => Some(rates),
            _ => None,
        })
        .expect("rate spectrum");
    let logs: Vec<_> = rates.iter().map(|rate| rate.log10()).collect();
    let mean = logs.iter().sum::<f64>() / logs.len() as f64;
    (logs.iter().map(|value| (value - mean).powi(2)).sum::<f64>() / logs.len() as f64).sqrt()
}

#[test]
fn kossel_rate_spectrum_broadens_as_the_etch_pit_nucleates_and_matches_golden() {
    let path = format!(
        "{}/../../examples/kossel-etchpit.toml",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(path).expect("Kossel fixture");
    let parsed: petra_deck::DeckFile = toml::from_str(&text).expect("deck parses");
    let deck = petra_deck::compile(&parsed).expect("deck compiles");
    let mut engine = deck.build_engine(Some(42)).expect("engine builds");
    let initial = spectrum_log10_width(&observe(&engine, &deck));
    let mut strategy = deck.strategy();
    for _ in 0..200 {
        engine
            .step_with(&mut strategy)
            .expect("etch-pit trajectory continues");
    }
    let nucleated = spectrum_log10_width(&observe(&engine, &deck));
    assert!(
        nucleated > initial,
        "initial={initial}, nucleated={nucleated}"
    );
    let actual = format!("initial_log10_std={initial:.9}\nnucleated_log10_std={nucleated:.9}\n");
    let golden_path = format!(
        "{}/tests/golden/kossel-rate-spectrum.txt",
        env!("CARGO_MANIFEST_DIR")
    );
    let expected = std::fs::read_to_string(golden_path)
        .unwrap_or_else(|error| panic!("rate-spectrum golden file ({error}); actual:\n{actual}"));
    assert_eq!(actual, expected);
}

fn binder_cumulant(size: usize, temperature: f64) -> f64 {
    let path = format!(
        "{}/../../examples/conformance/ising.toml",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(path).expect("Ising fixture");
    let mut parsed: petra_deck::DeckFile = toml::from_str(&text).expect("deck parses");
    parsed.grid.as_mut().expect("grid").dims = vec![size, size];
    parsed
        .execution
        .metropolis
        .as_mut()
        .expect("Metropolis")
        .temperature = Some(temperature);

    let deck = petra_deck::compile(&parsed).expect("deck compiles");
    let sites = size * size;
    let run = run_ensemble(
        &deck,
        &EnsembleConfig {
            replicas: deck.n_replicas,
            base_seed: deck.seed,
            steps: (500 * sites) as u64,
            burn_in: (300 * sites) as u64,
            sample_every: sites as u64,
            bootstrap_resamples: 200,
            bootstrap_seed: 17,
        },
    )
    .expect("Ising ensemble runs");
    let up = deck
        .state_names
        .iter()
        .position(|name| name == "spin.up")
        .expect("up");
    let mut m2 = 0.0;
    let mut m4 = 0.0;
    let mut samples = 0.0;
    for replica in &run.replicas {
        for sample in replica.samples.iter().skip(1) {
            let ObservableValue::StateCounts(counts) = &sample.values[0] else {
                panic!("state-count observable missing")
            };
            let magnetization = (2.0 * counts[up] as f64 - sites as f64) / sites as f64;
            let square = magnetization * magnetization;
            m2 += square;
            m4 += square * square;
            samples += 1.0;
        }
    }
    let mean2 = m2 / samples;
    1.0 - (m4 / samples) / (3.0 * mean2 * mean2)
}

#[test]
fn ising_binder_crossing_is_dogfooded_on_the_parallel_ensemble_framework() {
    const ONSAGER_TC: f64 = 2.269_185_314_213_022;
    let temperatures = [2.0, ONSAGER_TC, 2.6];
    let differences: Vec<_> = temperatures
        .iter()
        .map(|&temperature| binder_cumulant(8, temperature) - binder_cumulant(12, temperature))
        .collect();
    let crossing = temperatures
        .windows(2)
        .zip(differences.windows(2))
        .find_map(|(temperature, difference)| {
            (difference[0] * difference[1] <= 0.0).then(|| {
                temperature[0]
                    - difference[0] * (temperature[1] - temperature[0])
                        / (difference[1] - difference[0])
            })
        })
        .unwrap_or_else(|| panic!("Binder curves did not cross: {differences:?}"));
    assert!(
        (crossing - ONSAGER_TC).abs() <= 0.25,
        "crossing={crossing}, differences={differences:?}"
    );
}
