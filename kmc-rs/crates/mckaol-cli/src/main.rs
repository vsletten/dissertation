//! `mckaol` — the CLI that wires the port together. Port of `mckaol.cpp`.
//!
//! The C++ binary takes **no arguments**: it reads four fixed-name input
//! files from the current working directory and writes outputs beside them.
//! This port keeps the fixed file names (they are part of the model's
//! contract) but accepts one optional argument — the run directory — so you
//! don't have to `cd` into a fixture tree to use it. With no argument it
//! behaves exactly like the C++: everything relative to the CWD.
//!
//! Wired through M7: read all four inputs, run the deterministic structural
//! build in the C++'s exact order, write `start.msi` (the bitwise golden
//! gate), run the Monte Carlo loop — the generic `kmc_engine::step` driving
//! the `Kaolinite` model under the legacy `ran2` seed (spec B2) — with the
//! in-loop snapshot writers (`step{i}.dat` population rows, `step{i}.msi`
//! movie frames), then the shutdown quartet: `end.dat`, the surface maps,
//! the end-state XYZ (to `start.xyz` — the legacy misnomer, see below), and
//! `end.msi`. Every artifact is byte-gated against the C++ golden capture
//! (`tests/golden_m7.rs`).

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use kaolinite::{
    Kaolinite, State, create_lattice, find_pairs, populate_solid, terminate_lattice,
    terminate_surface,
};
use kmc_engine::{Ran2, SiteGraph, StepStop, step};
use kmc_io::LatticeView;

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

/// Create `path`, stream a writer body into it buffered, then **flush
/// explicitly** and surface every error.
///
/// \[IDIOM\] Buffered IO and the Drop/flush trap. A raw `File` is
/// unbuffered — every `write!` would be a syscall — so writers get a
/// `BufWriter`. But `BufWriter`'s buffer is flushed in its `Drop` impl,
/// and **`Drop` cannot return an error**: whatever the final flush hits
/// (disk full, quota, NFS hiccup) is silently swallowed if you let the
/// writer fall off scope. The C++ has the same hole (`ofstream` flushes in
/// its destructor and nobody checks `badbit`), and there it's invisible;
/// Rust makes the hole visible by *typing* flush as fallible. So: explicit
/// `flush()?` before drop, and the drop becomes a no-op. Reviewer's lens
/// for agent-written IO: a `BufWriter` that is never `flush()`ed (or
/// `into_inner()`ed) means write errors are being dropped on the floor —
/// clippy has no lint for it; only review catches it.
///
/// \[IDIOM\] `impl FnOnce(&mut ...)` — a closure parameter instead of C++'s
/// "open the stream at the top of every function" copy-paste. The
/// open/flush policy lives once, here; each call site contributes only the
/// body. `FnOnce` (not `Fn`) because the body runs exactly once and may
/// consume captures.
fn write_file(
    path: &Path,
    body: impl FnOnce(&mut BufWriter<File>) -> io::Result<()>,
) -> io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    body(&mut w)?;
    w.flush()
}

/// The fallible body of the program — `mckaol.cpp`'s `main`, stage by
/// stage, end to end.
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

    // WART (spec B1), preserved: `output::initDatafile()` deletes
    // `results.dat` at startup — and *nothing ever writes it*. The
    // documented appended time series is a ghost; the real data is
    // scattered across the one-row `step{i}.dat` files below. The delete
    // is kept because a directory carrying a stale `results.dat` from the
    // pre-bug C++ era should end a legacy run in the same state the C++
    // would leave it. Repair design: docs/REFORM_PLAN.md R8.
    let _ = std::fs::remove_file(run_dir.join("results.dat"));

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
    write_file(&msi_path, |w| {
        kmc_io::write_msi(
            w,
            &LatticeView::of_structure(&structure, &cell),
            "start",
            sim.drawbonds,
        )
    })?;
    println!("wrote {}", msi_path.display());

    // The writers need the tiling data after the Structure is dismembered
    // below (graph → engine, pair/lostal → model), so keep our own copy of
    // the immutable pieces. `coord` never changes after the build; cloning
    // it once here is the whole cost of the split-ownership design.
    let coord = structure.coord.clone();
    let params = structure.params;

    // Dynamics — mckaol.cpp's `for (i=0; i<nsteps; i++)` loop, generic.
    // Hand the built structure to the model, seed the legacy RNG, and step.
    let (mut graph, mut model) = Kaolinite::from_structure(structure, rxns);
    // [IDIOM] The legacy fixed seed (spec B2). A future `--legacy` flag will
    // toggle this against a correctly-seeded modern RNG (see docs/REFORM_PLAN);
    // today the port reproduces the 2001 behavior by default so parity holds.
    let mut rng = Ran2::legacy();
    let mut scratch = Vec::new();
    let mut time: f32 = 0.0;

    let mut done = 0i32;
    for i in 0..sim.nsteps {
        match step(&mut graph, &mut model, &mut rng, &mut scratch) {
            Ok(adv) => {
                time += adv.dt;
                done = i + 1;
            }
            Err(StepStop::NoEvents) => {
                eprintln!("no events possible; stopping after {i} steps");
                break;
            }
            Err(StepStop::ZeroRate) => {
                eprintln!("total rate is zero; stopping after {i} steps");
                break;
            }
            Err(StepStop::Model(e)) => {
                eprintln!("reaction failed at step {i}: {e}");
                break;
            }
        }

        let view = LatticeView {
            graph: &graph,
            coord: &coord,
            params,
            cell: &cell,
        };

        // Population snapshot — mckaol.cpp lines 84–89, AFTER the event
        // (step0.dat already reflects one applied event and time = dt0).
        // WART (spec B1): one row into its own truncated `step{i}.dat`; a
        // production 5M-step run scatters the series across thousands of
        // one-row files. Preserved for the byte gate; repair is REFORM_PLAN
        // R8. (`i == 0` is redundant beside `i % wsteps == 0` — kept, like
        // the C++, so the condition reads 1:1.)
        if sim.wsteps != 0 && (i % sim.wsteps == 0 || i == 0) {
            write_file(&run_dir.join(format!("step{i}.dat")), |w| {
                kmc_io::write_data(w, &graph, time)
            })?;
        }

        // Movie frame — mckaol.cpp lines 91–96: `step{i}.msi` every msteps,
        // skipping i == 0 (`msteps && i && i % msteps == 0`); the C++ label
        // and the file share the "step{i}" name.
        if sim.msteps != 0 && i != 0 && i % sim.msteps == 0 {
            let name = format!("step{i}");
            write_file(&run_dir.join(format!("{name}.msi")), |w| {
                kmc_io::write_msi(w, &view, &name, sim.drawbonds)
            })?;
        }
    }

    let (si_total, al_total) = cation_totals(&graph);
    println!(
        "ran {done} steps; simulated time {time:.6}; final populations: Si={si_total} Al={al_total}"
    );

    // Shutdown writes — mckaol.cpp lines 100–103, same order: end.dat,
    // the surface maps, the XYZ, end.msi. These run even after an early
    // break (the C++ falls through to them too).
    let view = LatticeView {
        graph: &graph,
        coord: &coord,
        params,
        cell: &cell,
    };
    write_file(&run_dir.join("end.dat"), |w| {
        kmc_io::write_data(w, &graph, time)
    })?;
    // writeSurf opens both files itself; the port opens them here so the
    // writer stays filesystem-free. Same names, same truncation.
    {
        let mut si = BufWriter::new(File::create(run_dir.join("surfSi.out"))?);
        let mut al = BufWriter::new(File::create(run_dir.join("surfAl.out"))?);
        kmc_io::write_surf(&mut si, &mut al, &view)?;
        si.flush()?;
        al.flush()?;
    }
    // WART: the C++ hard-codes the XYZ output name "start.xyz" and writes
    // it HERE, at shutdown — so the file named start.xyz holds the END
    // state, and no end.xyz exists (golden manifest: "no output/end.xyz
    // file is produced by this code path"). Preserved: the golden capture
    // has an end-state start.xyz and the byte gate demands one.
    write_file(&run_dir.join("start.xyz"), |w| kmc_io::write_xyz(w, &view))?;
    write_file(&run_dir.join("end.msi"), |w| {
        kmc_io::write_msi(w, &view, "end", sim.drawbonds)
    })?;
    println!("wrote end.dat, surfSi.out, surfAl.out, start.xyz (end-state; legacy name), end.msi");

    Ok(())
}

/// Count occupied Si and Al cation sites for the end-of-run summary line
/// (the full population breakdown lives in `end.dat`).
fn cation_totals(graph: &SiteGraph<State>) -> (usize, usize) {
    let mut si = 0;
    let mut al = 0;
    for s in &graph.sites {
        if s.state.is_edge() || !s.state.is_occupied() {
            continue;
        }
        match s.state.class_code() {
            1 => al += 1,
            2 => si += 1,
            _ => {}
        }
    }
    (si, al)
}
