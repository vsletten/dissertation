//! petra — run a deck.
//!
//! Usage: petra <deck.toml> [--steps N] [--seed S] [--out DIR] [--paranoid]
//!                          [--ensemble N]
//!
//! Writes `populations.csv` (step, time, one column per state) to the
//! output directory and prints a run summary. `--paranoid` re-derives every
//! site's events from scratch each report interval and asserts agreement
//! with the incrementally maintained tables (design doc §5.2).
//!
//! `--ensemble N` runs N independent trajectories with seeds S, S+1, …,
//! S+N-1, writes per-seed final states to `ensemble.csv`, and prints
//! mean ± std per state — the ensemble-first validation workflow of design
//! doc §8 (the thing the legacy seed swallow made impossible).

use std::io::Write as _;
use std::process::ExitCode;

use petra_core::Stop;

struct Args {
    deck: String,
    steps: Option<u64>,
    seed: Option<u64>,
    out: String,
    paranoid: bool,
    ensemble: Option<u64>,
    xyz: bool,
    /// Write trajectory artifacts (snapshot.pgif.json + events.jsonl) for
    /// graph-viz playback.
    viz: bool,
}

fn parse_args() -> Result<Args, String> {
    let mut args = std::env::args().skip(1);
    let mut deck = None;
    let mut steps = None;
    let mut seed = None;
    let mut out = ".".to_string();
    let mut paranoid = false;
    let mut ensemble = None;
    let mut xyz = false;
    let mut viz = false;
    while let Some(a) = args.next() {
        match a.as_str() {
            "--xyz" => xyz = true,
            "--viz" => viz = true,
            "--ensemble" => {
                let value = args
                    .next()
                    .ok_or("--ensemble needs a value")?
                    .parse()
                    .map_err(|e| format!("--ensemble: {e}"))?;
                if value == 0 {
                    return Err("--ensemble must be at least 1".into());
                }
                ensemble = Some(value);
            }
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
        deck: deck.ok_or(
            "usage: petra <deck.toml> [--steps N] [--seed S] [--out DIR] [--paranoid] [--ensemble N] [--xyz]",
        )?,
        steps,
        seed,
        out,
        paranoid,
        ensemble,
        xyz,
        viz,
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
    if args.viz && !matches!(deck.strategy, petra_deck::ExecutionStrategy::Ctmc) {
        return Err(
            "the current batch trajectory format supports only CTMC; discrete strategies fail closed"
                .to_string(),
        );
    }
    let ensemble = args.ensemble.unwrap_or(deck.n_replicas);
    if ensemble > 1 {
        return run_ensemble(&args, &deck, ensemble);
    }
    let mut engine = deck.build_engine(args.seed).map_err(|e| e.to_string())?;
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

    // Trajectory artifacts: initial snapshot now, events as they fire.
    let mut event_log = if args.viz {
        let snap_path = format!("{}/snapshot.pgif.json", args.out);
        std::fs::write(&snap_path, petra_io::snapshot_json(&deck, &engine))
            .map_err(|e| e.to_string())?;
        println!("wrote {snap_path}");
        let file = std::fs::File::create(format!("{}/events.jsonl", args.out))
            .map_err(|e| e.to_string())?;
        let writer = std::io::BufWriter::new(file);
        Some(
            petra_io::EventLogWriter::new(
                writer,
                &deck,
                args.seed.unwrap_or(deck.seed),
                engine.lattice.len(),
            )
            .map_err(|e| e.to_string())?,
        )
    } else {
        None
    };

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
        "deck '{}': {} sites, {} reactions, T = {} K, seed = {}, strategy = {}",
        deck.name,
        engine.lattice.len(),
        engine.reactions.len(),
        deck.temperature,
        args.seed.unwrap_or(deck.seed),
        deck.strategy.as_str(),
    );
    report(&engine, &mut csv)?;

    let mut stopped: Option<Stop> = None;
    let mut strategy = deck.strategy();
    for i in 1..=steps {
        match engine.step_with(&mut strategy) {
            Ok(outcome) => {
                if let Some(log) = &mut event_log {
                    for fired in &outcome.fired {
                        log.record(fired, &engine).map_err(|e| e.to_string())?;
                    }
                }
            }
            Err(stop) => {
                stopped = Some(stop);
                break;
            }
        }
        if i % report_every == 0 {
            report(&engine, &mut csv)?;
            if args.paranoid {
                engine
                    .paranoid_check()
                    .map_err(|e| format!("paranoid: {e}"))?;
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
    if let Some(log) = event_log {
        let events = log.events_written();
        let mut writer = log.into_inner();
        use std::io::Write as _;
        writer.flush().map_err(|e| e.to_string())?;
        println!("wrote {}/events.jsonl ({events} events)", args.out);
    }
    if args.xyz {
        let path = format!("{}/final.xyz", args.out);
        write_xyz(&engine, &deck, &path)?;
        println!("wrote {path}");
    }
    Ok(())
}

/// Snapshot the occupied sites as plain XYZ: element = occupant species,
/// position from the cell matrix, state name in the comment column.
fn write_xyz(
    engine: &petra_core::Engine,
    deck: &petra_deck::CompiledDeck,
    path: &str,
) -> Result<(), String> {
    let lat = &engine.lattice;
    let mut lines = Vec::new();
    for s in 0..lat.len() {
        let state = lat.states[s];
        let Some(species) = &deck.state_occupants[state.0 as usize] else {
            continue; // vacant
        };
        let (cell_coord, t) = lat.coords(s);
        let pos = deck.unit_cell.cell.to_cartesian(
            deck.unit_cell.sites[t].frac,
            [
                cell_coord[0] as i32,
                cell_coord[1] as i32,
                cell_coord[2] as i32,
            ],
        );
        lines.push(format!(
            "{species} {:.6} {:.6} {:.6} # {}",
            pos[0], pos[1], pos[2], deck.state_names[state.0 as usize]
        ));
    }
    let mut text = String::new();
    text.push_str(&format!("{}\n", lines.len()));
    text.push_str(&format!(
        "petra deck '{}' step {} time {:.6e}\n",
        deck.name, engine.step_count, engine.time
    ));
    for l in &lines {
        text.push_str(l);
        text.push('\n');
    }
    std::fs::write(path, text).map_err(|e| e.to_string())
}

fn run_ensemble(args: &Args, deck: &petra_deck::CompiledDeck, ensemble: u64) -> Result<(), String> {
    let steps = args.steps.unwrap_or(deck.steps);
    let base_seed = args.seed.unwrap_or(deck.seed);

    std::fs::create_dir_all(&args.out).map_err(|e| e.to_string())?;
    let csv_path = format!("{}/ensemble.csv", args.out);
    let mut csv = std::fs::File::create(&csv_path).map_err(|e| e.to_string())?;
    writeln!(csv, "seed,steps,time,{}", deck.state_names.join(",")).map_err(|e| e.to_string())?;

    println!(
        "deck '{}': ensemble of {} runs, seed policy {:?}, {} steps each",
        deck.name, ensemble, deck.seed_policy, steps
    );

    // Per-state running sums for mean/std over ensemble members.
    let mut sum = vec![0.0f64; deck.n_states];
    let mut sumsq = vec![0.0f64; deck.n_states];
    for k in 0..ensemble {
        let seed = petra_deck::replica_seed(base_seed, k, deck.seed_policy);
        let mut engine = deck.build_engine(Some(seed)).expect("engine builds");
        let mut strategy = deck.strategy();
        for _ in 0..steps {
            if engine.step_with(&mut strategy).is_err() {
                break;
            }
        }
        let counts = engine.state_counts(deck.n_states);
        let row: Vec<String> = counts.iter().map(|c| c.to_string()).collect();
        writeln!(
            csv,
            "{seed},{},{:.6e},{}",
            engine.step_count,
            engine.time,
            row.join(",")
        )
        .map_err(|e| e.to_string())?;
        for (i, &c) in counts.iter().enumerate() {
            sum[i] += c as f64;
            sumsq[i] += (c as f64) * (c as f64);
        }
    }

    let n = ensemble as f64;
    println!("final populations, mean ± std over {} members:", ensemble);
    for (i, name) in deck.state_names.iter().enumerate() {
        let mean = sum[i] / n;
        let var = (sumsq[i] / n - mean * mean).max(0.0);
        println!("  {name}: {mean:.2} ± {:.2}", var.sqrt());
    }
    println!("wrote {csv_path}");
    Ok(())
}
