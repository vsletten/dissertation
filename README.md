# kmc-rs — Rust port of the kaolinite KMC model

A faithful Rust port of `mckaol`, the 25-year-old C++ kinetic Monte Carlo
model of kaolinite dissolution/growth from Victor's dissertation
(`/mnt/data/vsletten/src/vsletten/dissertation/main/model`, read-only
reference). Ported behavior **as implemented**, warts included, per the
faithful spec — and written to double as a **Rust teaching codebase** for a
C++ veteran who reviews more than he writes.

## Repo status: greenfield, direct commits to `main`

This repo is greenfield and single-author (agent-built, milestone-per-commit:
M0, M1, M2, M3, ...). Victor's house worktree/branch rules resume the moment
this stops being greenfield — i.e. once anyone else has work in flight here,
changes go back to one-PR-one-branch-one-worktree.

## Start here (reading order)

1. **`docs/RUST_TOUR.md`** — the guided C++ → Rust tour, written for the
   reviewing architect. Every Rust idiom used in this codebase, what it buys,
   and what to look for when reviewing agent-written Rust.
2. `crates/kaolinite/src/state.rs` — the state encoding (the model's
   vocabulary, ported from `common.hpp`).
3. `crates/mckaol-cli/src/main.rs` — the whole program shape.
4. The `[IDIOM]` comment blocks throughout — each Rust idiom is explained at
   its first use, contrasted with its C++ counterpart.

## Layout

```
crates/
  kmc-engine/   generic KMC core: site graph (model-agnostic; no chemistry)
  kaolinite/    THE model: state codes 100–599, build rules, reactions
  kmc-io/       readers/writers for the legacy text formats + MSI
  mckaol-cli/   the binary (port of mckaol.cpp's main)
data/golden/    exact golden inputs + reference outputs (SHA-256-verified
                copies from mission-control projects/kmc/golden/, TASK-004)
docs/           RUST_TOUR.md and porting notes
```

Design authority: `mission-control/projects/kmc/` docs 01–05. The behavior
contract is `02-model-spec.md` (Part A = intended behavior, Part B = warts).
Preserved warts are marked `WART (spec BN)` at the exact line they live on.

## Milestones

| Milestone | Lands | Status |
|---|---|---|
| M0 | workspace, `data.sim` reader (bug-compatible, spec B2), CLI skeleton | done |
| M1 | `data.cell` → `UnitCell` | pending |
| M2 | lattice tiling + neighbor resolution (periodic/open BCs) | pending |
| M3 | structural build (`find_pairs`, `populate_solid`, `terminate_*`), `data.rxn` reader, MSI writer, **bitwise golden gate** | pending |
| M4+ | engine `Model` trait seam, `ran2`, reactions, dynamics parity | next |

## The M3 golden check

The structural build (everything before the first MC step) is deterministic
and RNG-free, so the port must reproduce the C++ **exactly** — not
statistically. The gate: build the initial lattice from the golden inputs and
write `start.msi`; the file must be **byte-identical** to the golden
`start.msi` captured from the C++ binary (g++ 13.3, `-O3 -ffast-math`).

Run it:

```bash
cargo test -p mckaol-cli --test golden_m3   # the gate itself
cargo test                                  # everything
cargo run -p mckaol-cli -- data/golden/inputs   # writes start.msi there; diff it yourself:
diff data/golden/inputs/start.msi data/golden/outputs/start.msi
```

Numeric-care notes (why bitwise equality is even possible, and where it was
delicate) live in `docs/RUST_TOUR.md` §"The bitwise gate".

## Toolchain

Stable Rust (built with 1.92). No dependencies outside std — by design
(design doc §9: "keep it boring"). `cargo test`, `cargo clippy`, `cargo doc`
all run clean.
