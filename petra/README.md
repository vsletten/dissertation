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
(P0–P6). This tree is the P0 milestone: the design plus a compiling,
tested skeleton that proves the shape end to end.

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
