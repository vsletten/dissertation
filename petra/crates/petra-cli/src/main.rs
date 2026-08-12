//! petra — run a deck.
//!
//! Usage: petra <deck.toml> [--steps N] [--seed S] [--out DIR] [--paranoid]
//!
//! Writes `populations.csv` (step, time, one column per state) to the
//! output directory and prints a run summary. `--paranoid` re-derives every
//! site's events from scratch each report interval and asserts agreement
//! with the incrementally maintained tables (design doc §5.2).

use std::io::Write as _;
use std::process::ExitCode;

use petra_core::Stop;

struct Args {
    deck: String,
    steps: Option<u64>,
    seed: Option<u64>,
    out: String,
    paranoid: bool,
}

fn parse_args() -> Result<Args, String> {
    let mut args = std::env::args().skip(1);
    let mut deck = None;
    let mut steps = None;
    let mut seed = None;
    let mut out = ".".to_string();
    let mut paranoid = false;
    while let Some(a) = args.next() {
        match a.as_str() {
            "--steps" => {
                steps = Some(
                    args.next()
                        .ok_or("--steps needs a value")?
                        .parse()
                        .map_err(|e| format!("--steps: {e}"))?,
                )
            }
            "--seed" => {
                seed = Some(
                    args.next()
                        .ok_or("--seed needs a value")?
                        .parse()
                        .map_err(|e| format!("--seed: {e}"))?,
                )
            }
            "--out" => out = args.next().ok_or("--out needs a value")?,
            "--paranoid" => paranoid = true,
            other if deck.is_none() && !other.starts_with('-') => deck = Some(other.to_string()),
            other => return Err(format!("unexpected argument '{other}'")),
        }
    }
    Ok(Args {
        deck: deck.ok_or("usage: petra <deck.toml> [--steps N] [--seed S] [--out DIR] [--paranoid]")?,
        steps,
        seed,
        out,
        paranoid,
    })
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("petra: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run() -> Result<(), String> {
    let args = parse_args()?;
    let deck = petra_deck::load(&args.deck).map_err(|e| e.to_string())?;
    let mut engine = deck.build_engine(args.seed);
    let steps = args.steps.unwrap_or(deck.steps);
    let report_every = if deck.report_every == 0 {
        steps.max(1)
    } else {
        deck.report_every
    };

    std::fs::create_dir_all(&args.out).map_err(|e| e.to_string())?;
    let csv_path = format!("{}/populations.csv", args.out);
    let mut csv = std::fs::File::create(&csv_path).map_err(|e| e.to_string())?;
    writeln!(csv, "step,time,{}", deck.state_names.join(",")).map_err(|e| e.to_string())?;

    let report = |engine: &petra_core::Engine, csv: &mut std::fs::File| -> Result<(), String> {
        let counts = engine.state_counts(deck.n_states);
        let row: Vec<String> = counts.iter().map(|c| c.to_string()).collect();
        writeln!(
            csv,
            "{},{:.6e},{}",
            engine.step_count,
            engine.time,
            row.join(",")
        )
        .map_err(|e| e.to_string())
    };

    println!(
        "deck '{}': {} sites, {} reactions, T = {} K, seed = {}",
        deck.name,
        engine.lattice.len(),
        engine.reactions.len(),
        deck.temperature,
        args.seed.unwrap_or(deck.seed),
    );
    report(&engine, &mut csv)?;

    let mut stopped: Option<Stop> = None;
    for i in 1..=steps {
        match engine.step() {
            Ok(_) => {}
            Err(stop) => {
                stopped = Some(stop);
                break;
            }
        }
        if i % report_every == 0 {
            report(&engine, &mut csv)?;
            if args.paranoid {
                engine.paranoid_check().map_err(|e| format!("paranoid: {e}"))?;
            }
        }
    }
    report(&engine, &mut csv)?;

    match stopped {
        Some(stop) => println!(
            "stopped early at step {}: {stop} (t = {:.6e})",
            engine.step_count, engine.time
        ),
        None => println!(
            "completed {} steps, simulated time {:.6e}",
            engine.step_count, engine.time
        ),
    }
    let counts = engine.state_counts(deck.n_states);
    for (name, count) in deck.state_names.iter().zip(&counts) {
        println!("  {name}: {count}");
    }
    println!("wrote {csv_path}");
    Ok(())
}
