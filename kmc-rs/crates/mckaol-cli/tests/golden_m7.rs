//! **The M7 output-writer gate**: run the real binary on the golden inputs
//! and byte-compare **every artifact class** against the C++ capture.
//!
//! The M3 gate certifies the structural build through one artifact
//! (`start.msi`); the M6 gate certifies the trajectory through an
//! instrumented trace. This gate certifies the *user-visible contract*: a
//! directory the legacy binary ran in and a directory this port ran in are
//! indistinguishable, file for file, byte for byte —
//!
//! * `start.msi` — initial state (M3's artifact, re-proven end-to-end);
//! * `step{0,1000,…,19000}.dat` — the 20 one-row population snapshots
//!   (WART spec B1: one file per row, faithfully), plus the exact count;
//! * `end.dat` — the shutdown population row;
//! * `surfSi.out` / `surfAl.out` — the surface x,y projections;
//! * `start.xyz` — the END-state XYZ under its legacy misnomer;
//! * `end.msi` — the end-state structure;
//! * `step{5000,10000,15000}.msi` — movie frames, gated in a second run
//!   against a supplementary C++ capture (the primary golden run's
//!   `msteps = 1e6` never fires the movie path; see
//!   `data/golden/outputs/movie-msteps5000/README.md`).
//!
//! Absence is part of the contract too: no `end.xyz`, no `results.dat`
//! (spec B1: it is deleted, never written), and no movie frames under the
//! primary config.
//!
//! \[IDIOM\] `env!("CARGO_BIN_EXE_<name>")` — cargo builds the crate's
//! binary for integration tests and hands its path over as a compile-time
//! env var, so the test drives the *actual shipped executable* (argument
//! parsing, exit code, file creation — things a library-level test can't
//! see) with no PATH games. `CARGO_TARGET_TMPDIR` is the matching
//! target-local scratch dir: unlike `std::env::temp_dir()`, it lives under
//! `target/` where a stray artifact can't outlive `cargo clean`.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

const MCKAOL: &str = env!("CARGO_BIN_EXE_mckaol-cli");

fn golden_dir() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/golden"))
}

/// Set up a run directory with the golden inputs (optionally swapping in a
/// different `data.sim`), run the binary, and return the directory.
fn run_in_fresh_dir(tag: &str, data_sim: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(tag);
    // Recreate from scratch so reruns can't inherit stale artifacts.
    if dir.exists() {
        fs::remove_dir_all(&dir).expect("clear stale run dir");
    }
    fs::create_dir_all(&dir).expect("create run dir");

    let inputs = golden_dir().join("inputs");
    for f in ["data.cell", "data.lattice", "data.rxn"] {
        fs::copy(inputs.join(f), dir.join(f)).expect(f);
    }
    fs::copy(inputs.join(data_sim), dir.join("data.sim")).expect("data.sim");

    let status = Command::new(MCKAOL)
        .arg(&dir)
        .status()
        .expect("mckaol binary runs");
    assert!(status.success(), "mckaol exited nonzero: {status:?}");
    dir
}

/// Byte-compare one artifact against its golden counterpart, with a
/// line-level report on mismatch (these are all text formats).
fn assert_matches_golden(run_dir: &Path, name: &str, golden_rel: &str) {
    let got = fs::read(run_dir.join(name))
        .unwrap_or_else(|e| panic!("output {name} missing from run dir: {e}"));
    let want = fs::read(golden_dir().join("outputs").join(golden_rel))
        .unwrap_or_else(|e| panic!("golden {golden_rel} missing: {e}"));
    if got != want {
        let got_s = String::from_utf8_lossy(&got);
        let want_s = String::from_utf8_lossy(&want);
        for (n, (g, w)) in got_s.lines().zip(want_s.lines()).enumerate() {
            assert_eq!(g, w, "{name}: first divergence at line {}", n + 1);
        }
        panic!(
            "{name}: same leading lines but different sizes (got {} bytes, want {})",
            got.len(),
            want.len()
        );
    }
}

/// THE gate: the primary golden config (20k steps, wsteps=1000,
/// msteps=1e6), every artifact class byte-identical to the C++.
#[test]
fn full_run_reproduces_every_golden_artifact_bitwise() {
    let dir = run_in_fresh_dir("m7-primary", "data.sim");

    // The structural artifact, now proven through the real binary.
    assert_matches_golden(&dir, "start.msi", "start.msi");

    // All 20 population snapshots (manifest policy "first 10 + last 10 +
    // every 100th" ≡ all of them at this scale) — and exactly 20 of them.
    for i in (0..20_000).step_by(1000) {
        let name = format!("step{i}.dat");
        assert_matches_golden(&dir, &name, &name);
    }
    let step_dats = fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| {
            let n = e.file_name().to_string_lossy().into_owned();
            n.starts_with("step") && n.ends_with(".dat")
        })
        .count();
    assert_eq!(step_dats, 20, "exactly the C++'s 20 step files, no extras");

    // Shutdown artifacts.
    assert_matches_golden(&dir, "end.dat", "end.dat");
    assert_matches_golden(&dir, "surfSi.out", "surfSi.out");
    assert_matches_golden(&dir, "surfAl.out", "surfAl.out");
    // The legacy misnomer: start.xyz holds the END state (module docs).
    assert_matches_golden(&dir, "start.xyz", "start.xyz");
    assert_matches_golden(&dir, "end.msi", "end.msi");

    // Absences the C++ also guarantees.
    assert!(
        !dir.join("end.xyz").exists(),
        "no end.xyz exists in the legacy contract"
    );
    assert!(
        !dir.join("results.dat").exists(),
        "results.dat is deleted, never written (spec B1)"
    );
    let movie_frames = fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| {
            let n = e.file_name().to_string_lossy().into_owned();
            n.starts_with("step") && n.ends_with(".msi")
        })
        .count();
    assert_eq!(
        movie_frames, 0,
        "msteps=1e6 > nsteps: the movie path must not fire"
    );
}

/// The movie-frame gate: same trajectory, `msteps = 5000`, three frames
/// byte-identical to the supplementary C++ capture — plus a spot check
/// that changing `msteps` perturbed nothing else.
#[test]
fn movie_frames_match_the_cpp() {
    let dir = run_in_fresh_dir("m7-movie", "data.sim.msteps5000");

    for i in [5000, 10000, 15000] {
        let name = format!("step{i}.msi");
        assert_matches_golden(&dir, &name, &format!("movie-msteps5000/{name}"));
    }
    // Output cadence must not feed back into the dynamics: the movie run's
    // trajectory-derived artifacts are the primary golden ones.
    assert_matches_golden(&dir, "end.msi", "end.msi");
    assert_matches_golden(&dir, "end.dat", "end.dat");

    // The C++ writes a frame only when `i % msteps == 0` for i in
    // 1..nsteps: 5000/10000/15000 — never 0 (`i &&` guard) and never
    // 20000 (loop ends first).
    let frames: Vec<String> = fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.starts_with("step") && n.ends_with(".msi"))
        .collect();
    assert_eq!(frames.len(), 3, "exactly three frames: {frames:?}");
}

/// A stale `results.dat` in the run directory is removed at startup —
/// the one observable effect of the C++'s `initDatafile` (spec B1).
#[test]
fn stale_results_dat_is_removed_like_initdatafile_does() {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("m7-initdat");
    if dir.exists() {
        fs::remove_dir_all(&dir).unwrap();
    }
    fs::create_dir_all(&dir).unwrap();
    let inputs = golden_dir().join("inputs");
    for f in ["data.cell", "data.lattice", "data.rxn", "data.sim"] {
        fs::copy(inputs.join(f), dir.join(f)).unwrap();
    }
    fs::write(dir.join("results.dat"), "ghost of the pre-bug era\n").unwrap();

    let status = Command::new(MCKAOL).arg(&dir).status().unwrap();
    assert!(status.success());
    assert!(
        !dir.join("results.dat").exists(),
        "initDatafile's remove() must be reproduced"
    );
}
