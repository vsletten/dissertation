# Petra — general lattice KMC for mineral elementary reactions

A greenfield generalization of the dissertation's kaolinite Monte Carlo
model (`legacy/cpp-model`, ported in `kmc-rs`): a chemistry-free
rejection-free KMC engine where **the mineral is data**. One input deck
declares an arbitrary crystal structure (unit cell, species, bond
topology), the local states its sites can be in, and a set of elementary
reactions with transition-state-theory kinetics; the same binary runs any
of them.

**Start with [`docs/DESIGN.md`](docs/DESIGN.md)** — the architecture,
the deck format, the rate physics, the engine algorithms, and the roadmap
(P0–P6). Landed so far: P0 (design + tested skeleton), per-deck energy
units (kcal/mol, kJ/mol, eV), tabulated nonlinear ΔEa modifiers
(`by_count`), the P2 analytic gates (two-state equilibrium, Langmuir
isotherm, waiting-time statistics), the `--ensemble` seed-sweep runner,
and **P4: the full kaolinite deck** (`examples/kaolinite.toml`) — the
dissertation model as data, structurally golden-gated site-by-site
against the kmc-rs builder and running 20k-step dynamics with paranoid
differential checks green. **P3 trajectory export** is in too: run any
deck with `--viz` and drop the resulting `events.jsonl` on graph-viz (or
use its bundled `kaolinite dissolution` demo) to watch the run with
play/scrub/speed controls. **Stage 2 is live too**: `petra-wasm` compiles
the engine to WebAssembly, and graph-viz's `kaolinite LIVE` demo runs the
simulation *in the browser* — edit the deck TOML in the panel, re-seed,
re-run; events generate on demand as you play or scrub. **Defects and
strain energy** are in too: declare `[[defects]]` (screw/edge dislocations)
in a deck and their analytic elastic fields modify reaction barriers —
the `kossel-etchpit` example drills an etch pit down a screw core, gated
by ensemble tests and watchable live. The design space (topology
rewiring, self-consistent relaxation) is worked out in docs/STRAIN.md.

## Layout

```
crates/
  petra-core/   chemistry-free runtime: site graph (CSR, variable
                coordination), interned states, TST rates, guarded
                rewrites, Fenwick-tree engine with incremental updates
  petra-deck/   the TOML deck schema + compiler to core runtime tables
  petra-cli/    the `petra` binary
examples/
  kossel.toml   tutorial deck: simple-cubic dissolution/growth
docs/
  DESIGN.md     the design document
```

## Quick start

```bash
cd petra
cargo test                                  # unit + end-to-end smoke gates
cargo run --release -p petra-cli -- examples/kossel.toml --out /tmp/kossel
```

The Kossel run dissolves an 8×8×8 simple-cubic crystal under an
undersaturated solution (Δμ = −2 kcal/mol) and writes per-state population
time series to `populations.csv`. `--paranoid` cross-checks the engine's
incremental event tables against from-scratch re-derivation every report
interval.
