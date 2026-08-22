use std::io::Write as _;
use std::path::PathBuf;

const TEMPERATURES: [f64; 11] = [
    6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0,
];
const REPLICAS: u64 = 16;
const ARRIVALS_PER_REPLICA: u64 = 300;
const MAX_EVENTS: u64 = 250_000;

fn efficiency(deck: &petra_deck::CompiledDeck, temperature: f64, seed: u64) -> Result<f64, String> {
    let mut engine = deck
        .build_engine(Some(seed))
        .map_err(|error| error.to_string())?;
    engine.set_temperature(temperature)?;

    let mut arrivals = 0_u64;
    let mut molecules = 0_u64;
    for _ in 0..MAX_EVENTS {
        let fired = engine.step().map_err(|stop| {
            format!("T={temperature} K seed={seed}: stopped after {arrivals} arrivals: {stop}")
        })?;
        let name = &deck.reactions[fired.reaction as usize].name;
        if name.starts_with("deposit_") {
            arrivals += 1;
        } else if name.starts_with("form_h2_") {
            molecules += 1;
        }
        if arrivals == ARRIVALS_PER_REPLICA {
            return Ok(2.0 * molecules as f64 / arrivals as f64);
        }
    }
    Err(format!(
        "T={temperature} K seed={seed}: bounded event budget exhausted at {arrivals} arrivals"
    ))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let deck_path = args
        .next()
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("decks/ice-mantle-h2.toml"));
    let output_path = args
        .next()
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("ice-mantle-efficiency.csv"));
    if args.next().is_some() {
        return Err("usage: h2_efficiency [deck.toml] [output.csv]".into());
    }

    let deck = petra_deck::load(&deck_path)?;
    if let Some(parent) = output_path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let mut output = std::io::BufWriter::new(std::fs::File::create(&output_path)?);
    writeln!(
        output,
        "temperature_k,mean_efficiency,ci95_low,ci95_high,replicas,arrivals_per_replica"
    )?;

    for temperature in TEMPERATURES {
        let mut values = Vec::with_capacity(REPLICAS as usize);
        for replica in 0..REPLICAS {
            let seed = petra_deck::replica_seed(7001, replica, petra_deck::SeedPolicy::Hash);
            values.push(efficiency(&deck, temperature, seed)?);
        }
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let variance = values
            .iter()
            .map(|value| (value - mean).powi(2))
            .sum::<f64>()
            / (values.len() - 1) as f64;
        let half_width = 1.96 * (variance / values.len() as f64).sqrt();
        writeln!(
            output,
            "{temperature:.1},{mean:.6},{:.6},{:.6},{REPLICAS},{ARRIVALS_PER_REPLICA}",
            (mean - half_width).max(0.0),
            (mean + half_width).min(1.0),
        )?;
        println!("T={temperature:>4.1} K  efficiency={mean:.3} ± {half_width:.3}");
    }
    output.flush()?;
    println!("wrote {}", output_path.display());
    Ok(())
}
