# Kaolinite

Kinetic Monte Carlo simulation of kaolinite surface reactions — the
dissertation project, three generations of code in one repo.

## The lineage

| Generation | Where | Status |
|---|---|---|
| **C → C++** (`mckaol`) | [`legacy/cpp-model/`](legacy/cpp-model/) | Historical reference — the dissertation-era simulation, preserved whole. The C→C++ evolution lives in this repo's history: `git log --follow legacy/cpp-model/<file>` walks back to the original C conversion commits. |
| **Rust port** | [`kmc-rs/`](kmc-rs/) | The living model. Milestone-gated faithful port: dynamics bitwise-identical to the C++ over 20,000 steps (M6), full-artifact golden gate green (M7). Doubles as a Rust teaching codebase — see [`kmc-rs/docs/RUST_TOUR.md`](kmc-rs/docs/RUST_TOUR.md). |
| **Viewer, v1** | [`legacy/viewer/`](legacy/viewer/) | The original React/TypeScript structure visualizer. Superseded. |
| **Viewer, v2** | [`graph-viz/`](graph-viz/) | The overhaul: PGIF-native instanced property-graph viewer. |

`kmc-rs/` and `graph-viz/` were imported from their standalone repos as
subtree merges on 2026-07-13 — their full commit histories are in this
repo's log.

## Quick start

```bash
# The living model
cd kmc-rs && cargo build --release && cargo test

# The visualizer
cd graph-viz && npm install && npm run dev

# The historical C++ (for parity reference / nostalgia)
cd legacy/cpp-model && make && ./mckaol
```

See [`CLAUDE.md`](CLAUDE.md) for the full map and per-component detail.
