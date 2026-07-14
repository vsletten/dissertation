# 06 — M0 Co-Build Session Plan (90 minutes)

The first weekend co-build session for the KMC Rust port. Victor drives the keyboard;
Claude navigates. Nothing is written ahead of the session — the port is Victor-paced
Rust learning (`projects/kmc-port.md` goal #2), not an autonomous build.

**Session goal (M0 from `03-rust-workspace-design.md` §6):** by the end of 90 minutes,
a Cargo workspace exists in a new `kmc` repo, and `cargo run -p mckaol-cli -- data/data.sim`
parses the golden `data.sim` and prints all **five** fields — which fixes legacy bug B2
(the swallowed seed) on day one. One passing test proves it.

---

## 1. Before the session (Victor's checklist, ~10 min solo)

The workstation already has the toolchain — this is a verify list, not an install list.

1. **Toolchain** (already present: rustc 1.92.0 via rustup):
   ```bash
   rustup update            # refresh stable
   rustc --version && cargo --version
   ```
   If working from another machine instead: install via https://rustup.rs
   (one `curl | sh`), then the same two version checks.
2. **Editor**: install **rust-analyzer** in VS Code (or the JetBrains Rust plugin).
   This matters more than anything else on this list — inline type hints and
   borrow-checker errors in the editor are half the learning value.
3. **Repo**: create the new port repo (design §2: its own repo, *not* inside the
   dissertation repo). Per house worktree rules, M0 work happens on a branch worktree:
   ```bash
   mkdir -p /mnt/data/vsletten/src/vsletten/kmc/main && cd $_ && git init -b main
   git commit --allow-empty -m "init"
   git -C /mnt/data/vsletten/src/vsletten/kmc/main worktree add \
     /mnt/data/vsletten/src/vsletten/kmc/m0-workspace-skeleton -b m0-workspace-skeleton
   ```
   One milestone = one branch = one PR — the ladder maps cleanly onto the house rules
   and gives each milestone a reviewable unit.
4. **Fixture**: copy the golden input into the new repo (verify the hash first):
   ```bash
   sha256sum projects/kmc/golden/inputs/data.sim   # expect a7db7cb729ec…
   mkdir -p /mnt/data/vsletten/src/vsletten/kmc/m0-workspace-skeleton/data
   cp projects/kmc/golden/inputs/data.sim /mnt/data/vsletten/src/vsletten/kmc/m0-workspace-skeleton/data/
   ```
5. **Skim (optional, 15 min)**: `03-rust-workspace-design.md` §1–2 (the three-crate
   separation) and `02-model-spec.md` B2 (the seed bug we're about to fix). Don't
   pre-study Rust — the session teaches what it uses.

No other tools. No cargo plugins, no clap, no serde — M0 is std-only by design
(design §9: keep it boring).

---

## 2. The 90 minutes

| Time | Block | What happens | Rust ideas taught |
|---|---|---|---|
| 0:00–0:10 | Orientation | Tour the five archaeology docs + the golden capture; the co-build contract (Victor types, Claude explains-then-watches); today's finish line | — |
| 0:10–0:25 | Workspace skeleton | Root `Cargo.toml` with `[workspace]`; `cargo new --lib crates/kmc-io`, `cargo new crates/mckaol-cli`; wire members; `cargo build` green | workspaces, crates vs modules, `Cargo.toml` anatomy |
| 0:25–0:55 | The `data.sim` reader | `SimParams` struct + `read_sim(path) -> Result<SimParams, …>` in `kmc-io::inputs` (see §3) | structs, `Result` + `?`, `BufRead` lines, `str::parse`, ownership of `String` vs `&str` |
| 0:55–1:05 | The test | `#[test]` parsing `data/data.sim`; assert all five fields (20000 / 1000 / 1000000 / **−2** / 1); `cargo test` | test modules, `assert_eq!`, fixture paths |
| 1:05–1:20 | The CLI | `mckaol-cli/src/main.rs`: take the path from `std::env::args`, call `read_sim`, pretty-print; run it against the fixture | binary vs library crate, `main() -> Result`, error display |
| 1:20–1:30 | Land it | `cargo fmt`, commit on the `m0-workspace-skeleton` branch, push, open the PR; preview M1 (`data.cell` → `UnitCell`); jot what was confusing | the milestone rhythm |

**Skeleton note:** M0 creates only the two crates it needs (`kmc-io`, `mckaol-cli`).
The design's other crates (`kmc-engine`, `kaolinite`) are added when their milestones
arrive (M4 and M2/M3) — empty crates now would just be clutter to explain.

**If running behind** (the reader block is where time goes): cut the CLI block —
the test already proves the reader; the CLI is a 10-minute solo follow-up.
**If running ahead**: start M1 — open `golden/inputs/data.cell` and sketch the
`UnitCell` struct together; or discuss where the `--legacy-seed` toggle (design §5)
will live.

---

## 3. The centerpiece: `read_sim`, and the story that motivates it

`data.sim` is five lines, each `value  # comment` (the seed line has a tab):

```
20000     # Number of steps in simulation
1000        # Number of steps between data writes (0 for no data)
1000000     # Number of steps between movie frames (0 for no movie)
-2	    # Seed for random number generator
1           # Draw bonds? (0 No - only occupied, 1 Yes - entire lattice)
```

Target shape (Victor writes it; this is the reference, not a handout to paste):

```rust
pub struct SimParams {
    pub nsteps: u64,
    pub wsteps: u64,     // 0 = no data writes
    pub msteps: u64,     // 0 = no movie frames
    pub seed: i64,       // negative allowed — ran2 convention
    pub drawbonds: bool,
}
```

Parse rule: first whitespace-separated token of each line, `parse()`d, `?` on error.

**Tell the B2 story while writing it** (spec `02-model-spec.md` B2): the C++ reads
`drawbonds` twice — the first read eats the seed line, `ranseed` stays 0, and every
run of the legacy model ever made used the same seed. The five-field Rust reader is
the port's first faithful *improvement*: read the seed correctly now, add a
`--legacy-seed` flag later (M6) for parity runs. Day one, the port already fixes a
25-year-old bug — that's the hook that makes the milestone memorable.

Definition of done for M0:
- [ ] `cargo build` and `cargo test` green from the workspace root
- [ ] test asserts all five fields from the golden fixture, including `seed == -2`
- [ ] `cargo run -p mckaol-cli -- data/data.sim` prints the params
- [ ] PR open from `m0-workspace-skeleton`

---

## 4. The M3 faithfulness checkpoint (where this ladder is headed)

M0–M3 climb to the port's first proof of faithfulness. Worth previewing in the
session's orientation block so every milestone has a visible destination.

**What M3 proves** (spec §C1): the structural setup — lattice build, `find_pairs`,
`populate_solid`, `terminate_surface`, `terminate_lattice` — is **deterministic and
RNG-free**. So the Rust port must reproduce the C++ initial configuration *exactly*,
not statistically. It's the strongest checkpoint in the whole validation strategy, and
it lands before any KMC dynamics exist.

**The reference**: `projects/kmc/golden/` (TASK-004, captured 2026-07-08) — exact
inputs with SHA-256 manifest (`data.sim` / `data.cell` / `data.lattice` 20×3 surface
plane 0 / `data.rxn`), outputs including `start.msi` (the full atom+bond model) and
`start.xyz` (the ~1,019-atom occupied set), compiler metadata, and a duplicate-run
byte-identical reproducibility check.

**Pass criteria at M3:**
1. **State-array diff (primary)**: the Rust `state[]` (and `nbr[]`) after structural
   setup matches the C++ bit-for-bit, same flat indexing (`a*bCells*npos + b*npos + n`).
2. **Atom-set cross-check (secondary)**: per-class occupied-site counts and the
   occupied set match `start.xyz` / `start.msi`.

**Prep task before M3 (bounded, queueable — Fable or Hermes):** the golden capture has
no raw `state[]` dump; `start.msi`/`start.xyz` are *renderings* of it. One-time task,
same scratch method as TASK-004 (copy the model source to scratch, dissertation repo
untouched): add a debug dump of `state[]` + `nbr[]` right after `terminate_lattice`,
run, archive as `golden/outputs/state0.dump` + `nbr0.dump` with hashes and a rerun
check, then discard the patched source. Queue it around M2 so it's waiting when M3 needs it.

**Two caveats, so M3 isn't over-claimed:**
- **Don't diff Cartesian coordinates against a `data.cell`-derived basis.** The C++
  writers use a hard-coded cell matrix and ignore `data.cell` (spec B5). Compare state
  arrays and fractional/site-index space at M3; bug-compatible MSI output is M7's job.
- **M3 parity says nothing about the nth-neighbor machinery.** The shipped `data.rxn`
  is flat (spec B4) — the environment feature is inert in-sample. Its real validation
  is M8, against an authored non-flat `data.rxn`.

---

## 5. Guardrails

- Dissertation repo (`dissertation/main`) stays read-only reference, always.
- The port repo follows house worktree rules: one milestone = one branch = one PR;
  never work on the `main` worktree.
- Co-build pace: Claude does not write milestone code outside sessions. Between
  sessions, Claude/Hermes may do *non-port* support tasks (like the state-dump
  capture above) via the queue.
- Don't gold-plate (kmc-port goal #3): std-only until PGIF (M9), Tier-1
  model-as-data (design §7), no parallelism talk before M8.
