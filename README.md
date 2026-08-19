# Kaolinite

Kinetic Monte Carlo simulation of kaolinite surface reactions — the
dissertation project, now several generations of code in one repo:
from the 1990s C model to a declarative lattice-KMC engine and a
GPU-accelerated quantum chemistry pipeline that feeds it.

## The lineage

| Generation | Where | Status |
|---|---|---|
| **C → C++** (`mckaol`) | [`legacy/cpp-model/`](legacy/cpp-model/) | Historical reference — the dissertation-era simulation, preserved whole. The C→C++ evolution lives in this repo's history: `git log --follow legacy/cpp-model/<file>` walks back to the original C conversion commits. |
| **Rust port** | [`kmc-rs/`](kmc-rs/) | The faithful port and golden standard. Milestone-gated: dynamics bitwise-identical to the C++ over 20,000 steps (M6), full-artifact golden gate green (M7). Doubles as a Rust teaching codebase — see [`kmc-rs/docs/RUST_TOUR.md`](kmc-rs/docs/RUST_TOUR.md). |
| **Petra** | [`petra/`](petra/) | The living engine: a chemistry-free, declarative lattice-KMC platform where **the mineral is data** — crystal structure, site states, and TST reactions declared in a TOML deck. Carries the full kaolinite model as a deck (golden-gated against kmc-rs), crystallographic defects + strain energy, trajectory record/replay, and live in-browser (WASM) simulation. Start at [`petra/README.md`](petra/README.md) and [`petra/docs/DESIGN.md`](petra/docs/DESIGN.md). |
| **Quarry (QM lab)** | [`qm/`](qm/) | The quantum chemistry side: GPU-accelerated (PySCF/GPU4PySCF) pipeline computing site-resolved activation energies — cluster building, saddle search (Sella/CI-NEB), quasi-RRHO → Eyring rates with tunneling — emitted as petra deck rate tables with per-number provenance. Design authority: [`qm/SURVEY.md`](qm/SURVEY.md); build plan: [`qm/HANDOFF.md`](qm/HANDOFF.md). |
| **Viewer, v1** | [`legacy/viewer/`](legacy/viewer/) | The original React/TypeScript structure visualizer. Superseded. |
| **Viewer, v2** | [`graph-viz/`](graph-viz/) | The overhaul: PGIF-native instanced property-graph viewer; renders petra trajectories live. |

`kmc-rs/` and `graph-viz/` were imported from their standalone repos as
subtree merges on 2026-07-13 — their full commit histories are in this
repo's log.

## Quick start

```bash
# The engine — run the kaolinite deck (or any deck)
cd petra && cargo run -p petra-cli -- examples/kaolinite.toml --steps 20000

# The QM lab — CPU-only test gate, or the GPU smoke test on a CUDA box
cd qm && uv sync --extra dev && uv run pytest

# The faithful port (parity reference)
cd kmc-rs && cargo build --release && cargo test

# The visualizer
cd graph-viz && npm install && npm run dev

# The historical C++ (for nostalgia)
cd legacy/cpp-model && make && ./mckaol
```

See [`CLAUDE.md`](CLAUDE.md) for the full map and per-component detail.
