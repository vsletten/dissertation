//! `mckaol` — the CLI that wires the port together. Port of `mckaol.cpp`.
//!
//! The C++ binary takes **no arguments**: it reads four fixed-name input
//! files from the current working directory and writes outputs beside them.
//! This port keeps the fixed file names (they are part of the model's
//! contract) but accepts one optional argument — the run directory — so you
//! don't have to `cd` into a fixture tree to use it. With no argument it
//! behaves exactly like the C++: everything relative to the CWD.
//!
//! At M0 the wiring stops after `data.sim`: read it, print it. Each later
//! milestone extends `main` one stage further down the C++ `main()`.

use std::path::PathBuf;
use std::process::ExitCode;

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

/// The fallible body of the program.
///
/// \[IDIOM\] `Box<dyn Error>` — "some error, I only need to display it".
/// Library crates define precise error enums (see `kmc_io::ReadError`);
/// binaries that just report-and-exit can erase the type behind a trait
/// object. `dyn` is Rust's *opt-in* dynamic dispatch: unlike C++ virtual
/// functions, dynamism is visible in the type and paid for only where
/// declared.
fn run() -> Result<(), Box<dyn std::error::Error>> {
    // [IDIOM] Iterator chains over index arithmetic. `args().nth(1)` is
    // "second element if present" — no argc bounds check to get wrong, no
    // argv[1] to dereference on faith. `unwrap_or_else` supplies the
    // default lazily.
    let run_dir: PathBuf = std::env::args_os()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));

    let sim = kmc_io::read_sim(&run_dir.join("data.sim"))?;
    let cell = kmc_io::read_cell(&run_dir.join("data.cell"))?;

    println!("data.sim ({}):", run_dir.join("data.sim").display());
    println!("  nsteps    = {}", sim.nsteps);
    println!("  wsteps    = {}", sim.wsteps);
    println!("  msteps    = {}", sim.msteps);
    println!("  drawbonds = {}", sim.drawbonds);
    println!(
        "  ranseed   = {}   (always 0: the legacy reader swallows the seed — spec B2)",
        sim.ranseed
    );

    println!("data.cell: {} positions, a/b/c = {} {} {}", cell.npos(), cell.a, cell.b, cell.c);

    Ok(())
}
