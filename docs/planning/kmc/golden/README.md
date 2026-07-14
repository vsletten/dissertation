# KMC Golden Reference Capture (TASK-004)

This directory contains the captured C++ reference artifacts for KMC model outputs.

## Source model and constraints
- Upstream source copied into scratch at
  `.scratch/task-004/dissertation-model` from:
  `/mnt/data/vsletten/src/vsletten/dissertation/main/model`
- Dissertation repo itself was not modified.
- Build done entirely from scratch copy.
- No source/toolchain changes were required to compile this code.

## Run used
- Truncated deterministic run for runtime control: `data.sim` set to `nsteps=20000`
- `wsteps=1000`, `msteps=1000000`, `drawbonds=1`, seed token retained as `-2` (legacy path)
- Program emitted `start.msi`, `start.xyz`, `end.msi`, `end.dat`, `surfSi.out`, `surfAl.out`, and `step{0..19000}.dat`
- `end.xyz` was not emitted by this code path.

## Files
- `inputs/` — exact input files used (with a preserved set of lattice presets)
- `outputs/` — model outputs from the run
- `manifest.md` — SHA-256 digests + configuration and reproducibility notes
- `build-patch.diff` — patch summary (no source patches required)
- `meta/` — compiler version and run/build logs
