# Petra session handoff — new decks and next features

*Written 2026-08-18 for a fresh session. Everything referenced is on `main`
(or in PR #11, auto-merging). Read this, then skim the docs in §1 order.*

## 0. One-paragraph state of the world

Petra is a chemistry-free lattice KMC engine where the mineral is a TOML
deck: arbitrary unit cell + bonds, named site states, guarded-rewrite
reactions with TST rates, build-time init passes, and (new) crystallographic
line defects with analytic strain fields. The dissertation's kaolinite model
runs as a deck, golden-gated **site-by-site** against the parity-tested
kmc-rs builder. Trajectories record to a PGIF snapshot + JSONL event log and
play in graph-viz (scrub both directions); the engine also compiles to WASM
and runs live in the browser with an editable deck panel. 30 Rust tests + 17
viewer tests green.

## 1. Read these, in order

1. `petra/docs/DESIGN.md` — architecture, deck format, roadmap table
   (P0–P6 with ✅ status) and recorded decisions (§10).
2. `petra/docs/STRAIN.md` — defects/strain: what's implemented (analytic
   fields), what's designed-not-built (screw topology rewiring §3.1,
   stacking-fault re-registry §3.2), open science questions (§6).
3. `petra/examples/kaolinite.toml` — the reference for deck authoring style;
   `kossel-etchpit.toml` for defects; `kossel.toml` for the minimal shape.
4. `petra/docs/kaolinite-{structure,dynamics}-spec.md` — the legacy-model
   extraction specs (needed only when touching kaolinite fidelity).

## 2. House workflow (Victor's standing instructions)

- Develop on `claude/<topic>-<suffix>` branches; **every push ends in a PR
  to `main`, merged** (enable auto-merge — branch protection requires 4
  CodeQL checks that report slowly and read "expected" for a while; don't
  hand-retry the merge, arm auto-merge and watch). **After merge, sync the
  local `main`** (`git fetch origin main:main`).
- A merged branch restarts from `origin/main` (same name is fine).
- GitHub's API intermittently 500s; retry pushes with backoff.
- Rust: every feature lands with gates (`cd petra && cargo test`). Viewer:
  `cd graph-viz && npm test && npm run typecheck` (needs `npm install`
  fresh; wasm rebuild via `npm run build:wasm` needs
  `rustup target add wasm32-unknown-unknown` + the repo-pinned wasm-pack —
  note the sandbox proxy blocks binaryen, which is why `wasm-opt = false`
  is set in petra-wasm's Cargo.toml).
- Rebuild + re-commit `graph-viz/public/petra-wasm/` whenever the deck
  schema or engine changes, or the live demos silently reject new decks.

## 3. The ask: new decks

Priority order balancing science value and machinery exercised:

1. **Gibbsite Al(OH)₃** — kaolinite's octahedral sheet alone. Reuses the
   Al/Al-OH-Al chemistry patterns; directly relevant literature (Guo et al.
   gibbsite step-edge barriers). Good first deck: medium complexity, real
   mineral.
2. **c-tiled multilayer kaolinite** — the current deck is the legacy single
   sheet (dims [20,3,1]). Tiling c (and deciding the interlayer bonds —
   the specs note the legacy model has none between layers; interlayer
   H-bonds would be new, deck-authored physics) unlocks: screw-dislocation
   *topology* (STRAIN.md §3.1), stacking faults, and honest basal vs edge
   comparisons. Discuss interlayer bonding with Victor before authoring.
3. **Quartz-like SiO₂ framework** — exercises the Q1/Q2/Q3 connectivity
   ladder; Liu & Ruiz Pestana 2024 published per-connectivity barriers
   (54/71/81 kJ/mol) that drop straight into `by_count` tables. This is the
   deck that meets the QM pipeline (qm/HANDOFF.md) halfway.
4. **Fe³⁺-substituted kaolinite** (roadmap P5) — needs the one missing
   feature: a substitution init rule (§4.1). Occupant/kind split means no
   schema change for states/reactions.
5. Fun/simple: two-species rock-salt Kossel (exercises multi-species
   adsorption competition), or an etch-pit *array* (multiple `[[defects]]`
   — superposition already works).

Deck-authoring gotchas learned so far: states of one kind are contiguous
(declaration order = `shift` ladder — put protonation counters in order);
qualified `Kind.state` names avoid ambiguity; forward reactions carry the
2-hop hydrolyzed-set guard, reverses don't; `frozen = false` filters on any
selector that could see the fixed boundary; reverse rates are Arrhenius with
**negative ea = −ΔE** to stay temperature-consistent; kaolinite runs at the
legacy T = 8000 K "relative-rate knob" — new decks at real T need real
barriers (see the qm handoff — that's where they'll come from).

## 4. Feature backlog (in rough value order)

1. **Substitution init rules** — `[[init]]` variant: replace occupant on a
   kind with probability p (seeded, deterministic) and/or explicit site
   lists. Small compiler feature; unlocks deck #4. Design sketch in
   DESIGN.md §6.
2. **Detailed-balance linter** (P2 leftover) — check forward/reverse pairs:
   rate consistency, `strain_scale` pairs differing by exactly −1
   (STRAIN.md §2.2), thermo-cycle sums per kind. Compile-time warnings.
3. **Screw-dislocation topology rewiring** — fully specified in STRAIN.md
   §3.1 (cut half-plane, retarget bonds ±b cells along the line, integer
   Burgers only, min-image). Needs a c-tiled deck to demo (BCF spiral!).
4. **Variant-table dump** (P1 leftover) — CLI mode printing each reaction's
   resolved rate over its reachable environments; review tool for decks.
5. **Auto-reverse generation** (`reverse = { de }` on a reaction) — the
   compiler writes the inverted reaction; only for effects it can invert.
6. Cosmetic viewer bug: legend counts are computed at bind time and go
   stale during playback (screenshot showed `vacant 0` mid-dissolution).
   Refresh counts on pause/scrub-end (not per frame).
7. Binary event log for very long runs; ensemble-runner parallelism;
   10⁵–10⁶-site perf pass (the engine's incremental updates should scale —
   unproven).

## 5. Open questions parked with Victor

- STRAIN.md §6: which defects matter in kaolinite (stacking vs line),
  isotropic vs anisotropic elasticity, β default, hollow-core handling.
- DESIGN.md §10.2: environment-differentiated barriers are a science
  question awaiting real numbers (the qm pipeline's job).
- Interlayer bonding for a c-tiled kaolinite (§3.2 above).

## 6. Quick map of the code

```
petra/crates/petra-core/     engine: lattice.rs (CSR + strain field),
                             reaction.rs (selectors/guards/modifiers/effects),
                             engine.rs (Fenwick tree, incremental updates,
                             paranoid_check), rate.rs, crystal.rs
petra/crates/petra-deck/     schema.rs (TOML shapes), compile.rs (interning,
                             init passes, defects, build_engine)
petra/crates/petra-deck/tests/  the gates — copy their patterns for new decks
petra/crates/petra-io/       PGIF snapshot + event log writers
petra/crates/petra-wasm/     browser bindings (WasmSim)
graph-viz/src/traj.ts,live.ts,render/scene.ts   playback + live sim
```

A new deck's definition of done: compiles; a shape test (counts, key
rates); a dynamics smoke with conservation + `paranoid_check`; if it has an
interesting visual, a bundled viewer demo (`graph-viz/public/<name>.deck.toml`
+ a `petra:live-<name>` option) — and a PR that merges.
