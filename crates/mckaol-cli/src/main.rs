//! `mckaol` — the CLI that wires the port together. Port of `mckaol.cpp`.
//!
//! The C++ binary takes **no arguments**: it reads four fixed-name input
//! files from the current working directory and writes outputs beside them.
//! This port keeps the fixed file names (they are part of the model's
//! contract) but accepts one optional argument — the run directory — so you
//! don't have to `cd` into a fixture tree to use it. With no argument it
//! behaves exactly like the C++: everything relative to the CWD.
//!
//! Wired through M3: read all four inputs, run the deterministic structural
//! build in the C++'s exact order, and write `start.msi` — the file the
//! golden gate diffs bitwise. The MC loop itself (event list, `DoEvent`,
//! the step outputs) is M4–M6 territory and intentionally absent.

use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;
use std::process::ExitCode;

use kaolinite::{create_lattice, find_pairs, populate_solid, terminate_lattice, terminate_surface};

/// Entry point.
///
/// \[IDIOM\] `main() -> ExitCode` and errors printed at the *edge*. The C++
/// model exits from deep inside parsers via `Myerr::die`; here every layer
/// below `main` returns `Result` and only this function turns an error into
/// a process exit code. One place decides policy; everything else reports
/// facts. When reviewing agent-written binaries, look for `panic!`/`unwrap`
/// in library code — policy decisions hiding below the edge.
fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("mckaol: {e}");
            ExitCode::FAILURE
        }
    }
}

/// The fallible body of the program — `mckaol.cpp`'s `main`, stage by
/// stage, through the structural build.
///
/// \[IDIOM\] `Box<dyn Error>` — "some error, I only need to display it".
/// Library crates define precise error enums (see `kmc_io::ReadError`);
/// binaries that just report-and-exit can erase the type behind a trait
/// object. `dyn` is Rust's *opt-in* dynamic dispatch: unlike C++ virtual
/// functions, dynamism is visible in the type and paid for only where
/// declared. The `?`s below convert `ReadError` and `io::Error` into the
/// box automatically (any `Error` type coerces).
fn run() -> Result<(), Box<dyn std::error::Error>> {
    // [IDIOM] Iterator chains over index arithmetic. `args().nth(1)` is
    // "second element if present" — no argc bounds check to get wrong, no
    // argv[1] to dereference on faith. `unwrap_or_else` supplies the
    // default lazily.
    let run_dir: PathBuf = std::env::args_os()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    // Input stage — same file order as mckaol.cpp lines 23–44.
    let sim = kmc_io::read_sim(&run_dir.join("data.sim"))?;
    let rxns = kmc_io::read_rxn(&run_dir.join("data.rxn"))?;
    let cell = kmc_io::read_cell(&run_dir.join("data.cell"))?;
    let lattice_params = kmc_io::read_lattice(&run_dir.join("data.lattice"))?;

    println!(
        "data.sim: nsteps={} wsteps={} msteps={} drawbonds={} ranseed={} (seed swallowed, spec B2)",
        sim.nsteps, sim.wsteps, sim.msteps, sim.drawbonds, sim.ranseed
    );
    println!(
        "data.rxn: T={} K, dmu_si={}, dmu_al={}, {} reactions (diffusion tail unported, spec B6)",
        rxns.temperature,
        rxns.dm_si,
        rxns.dm_al,
        rxns.reactions.len()
    );
    println!("data.cell: {} positions", cell.npos());

    // Structural build — mckaol.cpp lines 45–49, order is law.
    let mut structure = create_lattice(&cell, lattice_params);
    find_pairs(&mut structure);
    let report = populate_solid(&mut structure, rxns.dm_si, rxns.dm_al);
    // The C++ prints this from inside PopulateSolid; same words, same data.
    println!(
        "Detected {} conditions -- filling {} unit cells",
        report.condition, report.filled_cells
    );
    terminate_surface(&mut structure);
    terminate_lattice(&mut structure);

    let occupied = structure
        .graph
        .sites
        .iter()
        .filter(|s| s.state.is_occupied() && !s.state.is_edge())
        .count();
    println!(
        "lattice: {} x {} cells, surface plane {} -> {} sites, {} occupied",
        lattice_params.a_cells,
        lattice_params.b_cells,
        lattice_params.surface_plane,
        structure.graph.len(),
        occupied
    );

    // Initial-state snapshot — mckaol.cpp line 58. This is the M3 artifact.
    let msi_path = run_dir.join("start.msi");
    // [IDIOM] BufWriter: a raw File is unbuffered — every write! would be
    // a syscall; wrapping it buffers the way ofstream does internally.
    // Classic reviewer catch in agent-written code: unbuffered writers
    // inside a loop.
    let mut out = BufWriter::new(File::create(&msi_path)?);
    kmc_io::write_msi(&mut out, &structure, &cell, "start", sim.drawbonds)?;
    println!("wrote {}", msi_path.display());
    println!("(MC dynamics not yet ported — M0-M3 stop at the initial state)");

    Ok(())
}
