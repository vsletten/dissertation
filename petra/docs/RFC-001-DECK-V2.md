# RFC-001 — Deck schema v2 and the `UpdateStrategy` trait

- status: draft (open for Victor's review)
- author: hermes cloud worker (agents/B1-schema-rfc)
- date: 2026-08-19
- supersedes: `petra/docs/DESIGN.md` §3.4 (deck shape) and §5 (engine) — this
  RFC *extends* those sections; it does not fork them. Where this document and
  DESIGN.md disagree, this document wins after it is accepted.
- gates: B2 (engine refactor), B3 (conformance decks), B4 (ensembles +
  observables), and by extension C2/D3 (the showcase decks).

---

## 0. Summary

Petra is today a chemistry-free **continuous-time KMC** engine: a deck
declares a lattice, a state vocabulary, and guarded rewrite rules with TST
rates, and the engine runs rejection-free Gillespie dynamics on it. The
Omnibus program generalizes petra into a *multi-strategy* lattice simulation
engine — the same deck formalism, but the **update strategy** (how a state is
selected and how time advances) becomes a swappable component.

This RFC specifies, in implementable detail:

1. **Deck schema v2** — four orthogonal sections (`[structure]`,
   `[dynamics]`, `[execution]`, `[observables]`) that replace today's flat
   layout.
2. **The `UpdateStrategy` trait** — the Rust seam that owns selection-and-time,
   with the engine core owning everything else.
3. **v1 compatibility** — the shim that makes `kaolinite.toml` and
   `kossel.toml` load unchanged.
4. **The determinism contract** — how RNG streams are assigned per strategy so
   that `same deck + seed + strategy ⇒ bitwise trajectory` holds forever.
5. **Migration + test plan** for B2/B3, including the three conformance decks
   and their analytic acceptance numbers.

The two showcase models (pitting corrosion, ice-mantle astrochemistry) are
demonstrated as ~20-line deck fragments in §7, proving the schema serves decks
it has never met.

---

## 1. Design principles (non-negotiable)

1. **The deck is a document, not a program.** Everything a reviewer needs to
   understand a model — structure, rules, rates, strategy, observables — is
   declarative TOML. No per-model code.
2. **Orthogonality.** Structure, dynamics, execution, and observables are
   independent sections. A strategy change must not touch `[structure]` or
   `[dynamics]`; an observable change must not touch `[execution]`.
3. **v1 decks load unchanged.** The compat shim is a *load-time* transform,
   not a migration tool users run.
4. **Determinism is a law, not a feature.** Bitwise reproducibility per
   (deck, seed, strategy) is the platform's golden-gate discipline, promoted
   from kmc-rs's port discipline to a standing contract.
5. **The engine core owns state; strategies own selection-and-time.** The
   core owns the lattice, neighborhood queries, rule matching, RNG streams,
   and trajectory recording. A strategy is given a read-only view of the
   current state plus a handle to apply transitions; it decides *what fires*
   and *when*.

---

## 2. Deck schema v2

The top-level table gains a `schema = 2` marker and splits into four
sections. The v1 layout is accepted via the shim in §5.

```toml
[deck]
name = "…"
schema = 2            # NEW: absent ⇒ v1 ⇒ shim (§5)
units = "kcal/mol"    # unchanged

[structure]           # was: [cell] + [lattice] + [[kinds]] + [[species]] + [[init]] + [[defects]]
[dynamics]            # was: [[reactions]] + [[aliases]] + [thermo]
[execution]           # NEW: strategy + parameters + stopping + ensemble
[observables]         # NEW: declarative outputs
```

### 2.1 `[structure]` — lattice/topology

Carries over today's crystallographic machinery **unchanged in spirit**, and
adds simple-grid shorthands. Two mutually exclusive forms:

**Form A — explicit cell (today's form, moved under `[structure]`):**

```toml
[structure]
kind = "cell"                       # default; explicit for clarity

[structure.cell]                    # was [cell]
a = 3.0; b = 3.0; c = 3.0
alpha = 90.0; beta = 90.0; gamma = 90.0
# or: matrix = [[…],[…],[…]]

[[structure.cell.sites]]            # was [[cell.sites]]
kind = "M_site"
frac = [0.0, 0.0, 0.0]

[[structure.cell.bonds]]            # was [[cell.bonds]]
i = 0; j = 0; dcell = [1, 0, 0]
label = "pair"                      # optional, unchanged

[structure.lattice]                 # was [lattice]
dims = [8, 8, 8]
boundary = ["periodic", "periodic", "open"]
```

**Form B — simple-grid shorthand (NEW):**

```toml
[structure]
kind = "grid"
grid = "square" | "hex" | "cubic"   # lattice family
neighborhood = "moore" | "von_neumann"   # square/cubic only; hex is fixed 6-neighbor
dims = [64, 64]                     # 2D for square/hex, 3D for cubic
boundary = ["periodic", "periodic"] # per-axis, as today
```

The grid shorthand is **compiled into the same CSR adjacency** the explicit
cell produces — it is sugar, not a second code path. `square`+`moore` yields
8-neighbor adjacency, `square`+`von_neumann` 4-neighbor, `cubic`+`moore`
26-neighbor, `cubic`+`von_neumann` 6-neighbor, `hex` 6-neighbor. Each grid
site is a single site of a single kind (declared once, below); multi-kind
grids are expressed by declaring multiple site kinds and using `[[init]]`
passes to assign them (see §2.1.1). Grid sites have no fractional
coordinates; positions are synthesized from `dims` and the family for
snapshot export.

**Site kinds and state vocabularies carry over unchanged** — `[[kinds]]` with
`initial` and `[[kinds.states]]` (name + occupant) move under `[structure]`
verbatim. `[[species]]`, `[[init]]`, and `[[defects]]` also move under
`[structure]` unchanged. The occupant/kind split (§3.2 of DESIGN.md) is
preserved — it is what makes substitutional defects and the corrosion/ice
decks cheap.

#### 2.1.1 Grid site-kind assignment

A grid declares its default site kind:

```toml
[structure]
kind = "grid"
grid = "square"
neighborhood = "moore"
dims = [64, 64]
boundary = ["periodic", "periodic"]
default_kind = "surf"              # every grid site starts as this kind
```

Multi-kind grids (e.g. a corrosion deck with a `metal` substrate and a `film`
surface) are built by declaring the extra kinds and assigning them with
`[[init]]` region passes — the same machinery that already exists. No new
init semantics are introduced by the grid shorthand.

### 2.2 `[dynamics]` — the transition system

Today's `[[reactions]]` generalize to **rules**. The center/guards/effects/
modifiers vocabulary is **unchanged in spirit** — a rule is still a guarded
local rewrite. What changes is that the **rate field is interpreted per
strategy** (§2.2.1), and the section gains a `[dynamics.thermo]` block
(today's `[thermo]`, moved) plus the state aliases (`[dynamics.aliases]`,
today's top-level `[aliases]`).

```toml
[dynamics]

[dynamics.thermo]                  # was [thermo]
temperature = 298.15
[dynamics.thermo.mu]
M = -2.0
[dynamics.thermo.activity]
M = 1.0

[dynamics.aliases]                 # was [aliases]
hydrolyzed_any = ["O_SiSi.hy", "O_SiAl.alhy", "…"]

[[dynamics.rules]]                 # was [[reactions]]
name = "detach"
center = { kind = "M_site", state = ["occupied"] }
guards = [ { distance = 1, min = 1, state = ["occupied"] } ]
rate = { arrhenius = { prefactor = 1.0e13, ea = 10.0 } }
[[dynamics.rules.modifiers]]
select = { distance = 1, state = ["occupied"] }
per_match = { dea = 5.0 }
[[dynamics.rules.effects]]
target = "center"
set = { occupant = "vacant", config = "empty" }
```

**Effect selection modes.** An effect with `target = "neighbors"` applies
to **all** matching neighbors (`select_mode = "all"`, the v1 default and
v1-compatible behavior). Diffusion hops and pair reactions need
exactly-one semantics: `select_mode = "one"` applies the effect to a
single matching neighbor chosen by **one uniform draw** over the matches,
routed through the strategy's dynamics stream at the documented position
in the draw order (§4.1). A rule whose physics is "move/react with one
neighbor" and whose effects use `"all"` is almost certainly a bug; the
compiler warns when a rule both removes its center occupant and applies an
occupant-writing `"all"` effect to neighbors.

```toml
[[dynamics.rules.effects]]
target = "neighbors"
select = { distance = 1, state = ["empty"] }
select_mode = "one"      # exactly one match receives the effect
set = "H"
```

The `rate` field is **required on every rule** (it is the propensity for
CTMC, the probability for PCA, the energy delta for Metropolis, the truth
table for deterministic CA). A rule may additionally carry a
`rate_table`/`truth` field for strategies that need a different quantity —
see §2.2.1.

#### 2.2.1 Rate-field interpretation per strategy

| Strategy | `rate` field means | Extra fields a rule may carry |
|---|---|---|
| `ctmc` | propensity `k` (s⁻¹), as today | — |
| `synchronous` | deterministic CA: `rate` is ignored; the rule fires iff its guards hold | `truth = "and"` (default) — all guards must hold; no rate math |
| `metropolis` | energy delta `ΔE` (deck units) for the proposed transition | `rate = { energy = { delta = 1.0 } }` — see below |
| `pca` | probability `p ∈ [0,1]` per step | `rate = { probability = 0.3 }` |

Concretely, the `rate` field becomes a tagged union with **one variant per
strategy family**, and the compiler validates that the variant matches the
selected strategy:

```toml
# CTMC (unchanged from today)
rate = { arrhenius = { prefactor = 1.0e13, ea = 10.0 } }
rate = { eyring = { dh = 60.0, ds = -0.02 } }
rate = { constant = 1.0e-3 }

# Metropolis — the rule proposes a transition with an energy change
rate = { energy = { delta = 1.0 } }        # ΔE in deck units; acceptance = min(1, exp(-ΔE/kT))

# PCA — per-step probability
rate = { probability = 0.3 }

# Deterministic CA — no rate; guards alone decide
# (rate field omitted; `truth = "and"` is the default)
```

**Non-Arrhenius / fixed-rate escape hatch (Track D requirement).** The
astrochemistry deck needs a temperature-independent tunneling floor rate and
flux-driven source events. Two additions cover this without a schema break:

- `rate = { constant = 1.0e-3 }` already exists and is temperature-independent
  — the tunneling floor is a `constant` rate.
- **Source/bath events** (deposition, flux-driven adsorption) are a new rule
  *effect target*: `target = "source"` writes an occupant into a vacant site
  at a per-site rate `F·σ` with no center selector. See §7.2.

**Quenched per-site disorder (Track D requirement) — binned kinds.** The
ice deck needs site-dependent binding energies sampled from a distribution.
The state model has no per-site scalar (sites carry discrete kind/state/
occupant — deliberately), so disorder is represented as **binned discrete
kinds**: the deck declares N site kinds (`site_b1..site_bN`, N ≈ 8 energy
bins, each bin's nominal energy documented in the deck), a seeded
`[[structure.init]]` pass assigns each grid site a kind by the declared bin
probabilities (drawing from the init stream, §4.1), and rules carry
per-bin rates via the existing `when` machinery keyed on the center's
kind. This is expressible with today's machinery, deterministic, and
honest about energy resolution; continuous per-site scalars would be a
state-model change and are a *new RFC*, not a B2 detail. B2 finalizes only
the table *syntax*, not the representation.

### 2.3 `[execution]` — strategy + parameters + stopping + ensemble

```toml
[execution]
strategy = "ctmc" | "synchronous" | "metropolis" | "pca"

# strategy parameters (per-strategy; see §3)
[execution.ctmc]        # nothing required — defaults are the current engine
[execution.synchronous]
conflict_resolution = "first_match" # required: first rule in definition order wins per site
[execution.metropolis]
temperature = 298.15   # defaults to [dynamics.thermo].temperature
[execution.pca]

# stopping conditions (any combination; first hit wins)
[execution.stop]
steps = 20000          # max steps (all strategies)
time = 1.0e6           # simulated time (ctmc only; ignored elsewhere)
predicate = "…"        # deck-declared predicate, as today's roadmap (§5.4 DESIGN.md)

# ensemble spec
[execution.ensemble]
n_replicas = 1         # default 1
seed_policy = "increment" | "hash"   # see §4.2
```

`strategy` is the single source of truth for which `UpdateStrategy`
implementation runs. The `[execution.stop]` block replaces today's
`[simulation] steps`/`report_every` (see §5 for the shim).

### 2.4 `[observables]` — declarative outputs

```toml
[observables]
report_every = 2000            # steps between time-series samples (all strategies)

[[observables.series]]
kind = "state_counts"          # population per (kind, state) — today's CSV
[[observables.series]]
kind = "event_rates"           # per-rule firing rate time series
[[observables.series]]
kind = "rate_spectra"          # distribution of instantaneous per-site rates (Fischer/Luttge)
[[observables.series]]
kind = "cluster_sizes"         # cluster-size distribution (connected-component histogram)
[[observables.series]]
kind = "interface_roughness"   # surface width W(L,t) — needs a surface axis declared
[[observables.series]]
kind = "snapshot"              # extended-XYZ / PGIF snapshot at cadence
```

Each `series` may carry parameters (e.g. `interface_roughness` needs
`axis = 2`; `cluster_sizes` needs a `state` filter). The **first four** are
the B4 deliverables (see §6.2); `interface_roughness` and `snapshot` are
specified here so B4's implementer makes no design decisions, but their
implementation may land in B4 or later. The `[observables]` section is
**strategy-agnostic**: every observable is computed from the trajectory
record, which the engine core produces identically for every strategy.

---

## 3. The `UpdateStrategy` trait

The engine core (`petra-core`) is refactored so that the current CTMC loop is
**one implementation** of a trait. The core owns:

- lattice state (`Lattice`), site kinds, interned `StateId`s;
- neighborhood queries (`sites_at_distance`, `first_match`, `all_matches`);
- rule matching (`guards_pass`, `resolve_rate` for CTMC);
- RNG streams (§4);
- trajectory recording (the `Fired` event log + `last_changes`).

Strategies own **selection and time**. The trait:

```rust
/// A strategy decides, given the current lattice state, which transition
/// fires and how time advances. The core owns state, matching, RNG, and
/// recording; the strategy owns selection-and-time.
pub trait UpdateStrategy {
    /// Advance the simulation by one strategy step. Returns the fired
    /// transition(s) and the time increment, or a stop reason.
    ///
    /// `ctx` gives read access to the lattice and the compiled rules, and a
    /// handle to apply a transition (which the core records). The strategy
    /// must draw all randomness from `ctx.rng` in the documented order
    /// (§4) — never from its own RNG.
    fn step(&mut self, ctx: &mut StepCtx) -> Result<StepOutcome, Stop>;
}

/// What the core hands a strategy each step.
pub struct StepCtx<'a> {
    pub lattice: &'a Lattice,          // read-only current state
    pub rules: &'a [Rule],             // compiled rules (dense ids)
    pub rng: &'a mut dyn RngCore,      // the strategy's assigned stream (§4)
    pub apply: ApplyHandle<'a>,        // apply a transition; core records it
}

/// The result of one strategy step.
pub struct StepOutcome {
    pub fired: Vec<Fired>,             // one entry for CTMC/Metropolis/PCA; N for synchronous
    pub dt: f64,                       // time increment (0.0 for discrete-time strategies)
}
```

The four initial implementations:

| Strategy | `step` semantics | `dt` |
|---|---|---|
| `ExactCtmc` | today's two-level Fenwick selection + Poisson wait | `-ln(u)/R_total` |
| `SynchronousCA` | evaluate every rule at every site against a **double buffer**; apply all that fire simultaneously | `1.0` (or `0.0`; see below) |
| `AsyncMetropolis` | pick a site (uniform or weighted), propose a rule, accept with `min(1, exp(-ΔE/kT))` | `1/N` (per-site sweep convention) |
| `DiscreteTimePCA` | per-site independent Bernoulli draws from `rate.probability` | `1.0` |

**SynchronousCA double-buffer contract.** All guard evaluations read the
pre-step buffer; all writes land in the post-step buffer; the buffers swap at
the end of the step. This is what makes Conway's glider deterministic and
update-order-free. **Conflict resolution**: when multiple rules match the
same center site in the same step, the **first matching rule in rule-index
order wins** for that site — deterministic, draw-free, and documented as
part of the contract (§4.1); Conway's disjoint guards never exercise this,
which is exactly why it must be specified rather than discovered.
Mechanically, the strategy evaluates every site against the pre-step
buffer, collects `StepOutcome.fired`, and the core's `ApplyHandle`
**batch-applies** the whole set at end-of-step followed by one
dirty-propagation pass — the batch apply is what makes the double-buffer
contract real (and is the same capability tau-leaping reserves, §3).
`dt` is `1.0` per step (discrete time); the trajectory record stores
`step` as the time coordinate for discrete strategies.

**AsyncMetropolis selection.** Site selection is uniform over live sites
(one draw), then a rule is proposed by **one uniform draw over the site's
enabled rules** — a fixed proposal order was considered and rejected:
always proposing the first enabled rule silently breaks ergodicity for any
site with ≥2 enabled rules, and the single-rule Ising gate would never
catch it. The proposal draw is routed through the dynamics stream like
every other draw, so determinism is unaffected. Acceptance is
`min(1, exp(-ΔE/kT))` with `ΔE` from `rate.energy.delta`. `dt = 1/N`
(one Monte Carlo sweep = N site updates).

**What tau-leaping needs later (named, not built).** Tau-leaping is a CTMC
*acceleration*, not a new strategy: it needs (a) a bounded propensity-change
criterion (the "leap condition"), (b) Poisson sampling of per-rule event
counts over a leap `τ`, and (c) a rejection/rounding pass for negative
populations. It slots in as a *mode* on `ExactCtmc` (or a fifth strategy
`ctmc-tau`) and is **out of scope for B2/B3** — this RFC only reserves the
name and the requirement that the core's `apply` handle be batch-capable
(apply N events, then one dirty-propagation pass), which the trait already
supports via `StepOutcome.fired: Vec<Fired>`.

---

## 4. Determinism contract

**Law:** same deck + seed + strategy ⇒ bitwise-identical trajectory, per
strategy, forever. This is the golden-gate discipline that made kmc-rs
trustworthy, promoted to platform law.

### 4.1 RNG stream assignment

Every strategy gets its randomness from **one PCG64 stream per run**, seeded
from the deck seed, with **fixed per-stream salts** so that (a) init draws,
(b) dynamics draws, and (c) per-replica draws never collide, and (d) the
draw *order* is a documented contract.

```
seed_dynamics = seed ^ DYNAMICS_STREAM_SALT        # 0x… (fixed constant)
seed_init     = seed ^ INIT_STREAM_SALT            # already exists (0x9E37_79B9_7F4A_7C15)
seed_replica_k = seed ^ REPLICA_STREAM_SALT ^ k    # per-replica, §4.2
```

The **draw order within a step** is part of the contract and is specified per
strategy:

- **ExactCtmc** (unchanged from today): (1) site selection draw
  `u·R_total`, (2) within-site event draw, (3) Poisson wait draw, (4) branch
  draws in effect order, (5) `select_mode = "one"` draws in effect order
  (one uniform draw per such effect, over its matches in site-index order).
  Draws 1–4 are the current `Engine::step` order — the refactor must
  preserve them **exactly**, or the B2 parity gates fail (v1 decks contain
  no `select_mode = "one"` effects, so draw 5 cannot perturb parity).
- **SynchronousCA**: no draws (deterministic). The double-buffer swap is
  order-independent by construction; same-site multi-rule conflicts resolve
  first-match in rule-index order (§3), draw-free.
- **AsyncMetropolis**: (1) site draw, (2) proposal draw — one uniform draw
  over the site’s enabled rules in rule-index order (§3; a fixed proposal
  order was rejected as ergodicity-breaking for multi-rule sites),
  (3) acceptance draw, (4) `select_mode = "one"` draws as in ExactCtmc.
- **DiscreteTimePCA**: one Bernoulli draw per site per enabled rule, in
  site-index order then rule-index order; then `select_mode = "one"` draws
  for fired rules in the same order.

**No hidden draws, no global RNG.** The seed-swallow lesson (DESIGN.md §3.4)
is inverted and extended: every strategy's draw sequence is enumerable from
the deck.

### 4.2 Ensemble seed policy

`[execution.ensemble] seed_policy`:

- `"increment"` (default): replica `k` uses `seed + k` (today's CLI
  `--ensemble` behavior).
- `"hash"`: replica `k` uses `hash(seed, k)` — a fixed 64-bit mix (e.g.
  splitmix64 of `seed` and `k`), so replicas are decorrelated even for
  adjacent seeds. This is the recommended policy for B4's bootstrap CIs.

Both policies are deterministic functions of `(seed, k)`; the choice is a
deck field, not a runtime flag, so a deck's ensemble is reproducible.

### 4.3 What breaks the contract (and is therefore forbidden)

- Any `HashMap`/`HashSet` iteration in a hot path (use `BTreeMap`/`Vec`).
- Any `f32` in rate or time math (f64 everywhere — already the rule).
- Any parallelism *within* a single trajectory (rayon is for *across*
  replicas, never within one — see B4).
- Any draw not routed through the strategy's assigned stream.

---

## 5. v1 compatibility shim

v1 decks (`kaolinite.toml`, `kossel.toml`, `kossel-etchpit.toml`,
`kaolinite-multilayer.toml`) load **unchanged**. The shim is a load-time
transform in `petra-deck`, applied when `[deck] schema` is absent:

| v1 key | v2 location | transform |
|---|---|---|
| `[cell]`, `[[cell.sites]]`, `[[cell.bonds]]` | `[structure.cell]` etc. | move |
| `[lattice]` | `[structure.lattice]` | move |
| `[[kinds]]`, `[[species]]`, `[[init]]`, `[[defects]]` | `[structure.*]` | move |
| `[[reactions]]` | `[[dynamics.rules]]` | move |
| `[aliases]` | `[dynamics.aliases]` | move |
| `[thermo]` | `[dynamics.thermo]` | move |
| `[simulation] steps/seed/report_every` | `[execution.stop.steps]`, `[execution.ensemble]`, `[observables.report_every]` | split |
| *(absent)* | `[execution] strategy` | **imply `strategy = "ctmc"`** |

The shim is **pure and total**: it produces a v2 `DeckFile` that compiles
identically to the v1 path, so the B2 parity gates (§6.1) are meaningful.
The shim is unit-tested by round-tripping every shipped v1 deck through
`parse → shim → compile` and asserting the compiled tables are unchanged.

**`strategy = "ctmc"` is the implied default** for any deck that does not
declare one — this is the single most important compat rule, because it means
every existing deck is a CTMC deck and the CTMC strategy must reproduce
today's engine bit-for-bit.

---

## 6. Migration + test plan (B2/B3)

### 6.1 B2 — engine refactor behind the trait

Pure refactor, no new strategies. Acceptance:

- **Bitwise-parity gates**: `kaolinite.toml` and `kossel.toml` produce
  **identical trajectories** (same seeds) before and after the refactor, over
  long runs (≥20,000 steps, matching the existing golden-gate discipline).
  The gate is: run the pre-refactor engine and the post-refactor
  `ExactCtmc` strategy on the same deck+seed, assert the full `Fired` event
  log is byte-identical.
- **Schema-v2 round-trip tests**: every shipped v1 deck parses through the
  shim and compiles to the same tables as the v1 path.
- `cargo test` green, including the existing analytic gates (two-state,
  Langmuir, waiting-time) and `--paranoid` differential mode.

The refactor is **mechanical but delicate**: the `ExactCtmc` strategy must
preserve the exact draw order and float arithmetic of today's `Engine::step`
(§4.1). Any deviation is a parity-gate failure and is wrong until proven
otherwise.

### 6.2 B3 — conformance decks (the "it's general now" proof)

Three tiny decks with analytic answers, kept as CI gates. Each is a
`[structure] kind = "grid"` deck — the grid shorthand is exercised here for
the first time.

| Deck | Strategy | Analytic gate | Acceptance number |
|---|---|---|---|
| **Conway's Life** | `synchronous` | glider period + displacement | period 4, displacement (1,1) per period; known oscillators (blinker period 2, block still-life) |
| **2-D Ising** | `metropolis` | critical temperature via Binder-cumulant crossing | `T_c = 2/ln(1+√2) ≈ 2.269` (Onsager), within stated tolerance (Binder cumulant `U₄` curves for different `L` cross at `T_c`) |
| **SIR epidemic** | `pca` | epidemic threshold vs mean-field | threshold `R₀ = 1` (mean-field `β/γ = 1`); outbreak probability vs `R₀` matches the mean-field prediction within stated tolerance |

The **Ising gate doubles as the ensemble-statistics smoke test** (B4
dogfoods it — see B4 card): the Binder cumulant requires an ensemble of
replicas, so B3's Ising gate is the first real consumer of
`[execution.ensemble]`.

### 6.3 B4 — ensembles + observables

`[observables]` and `[execution.ensemble]` per this RFC: multi-seed parallel
replicas (rayon **across** replicas, never within one — §4.3), aggregation
(means, distributions, bootstrap CIs), and the first observable
implementations chosen to serve A5-Phase-0: `event_rates`, `rate_spectra`,
`state_counts`, `cluster_sizes`, and surface-area accounting. The Ising
Binder-cumulant gate is reimplemented **on** this framework (dogfood), and
`rate_spectra` is validated on the kossel etch-pit deck (spectrum broadens as
pits nucleate — qualitative gate + golden file).

---

## 7. Showcase deck fragments

Both showcase models must be expressible in the proposed schema. Each
fragment is ~20 lines and demonstrates the schema serving a deck it has
never met.

### 7.1 Pitting corrosion — passive-film patch (Track C)

From `docs/scoping/pitting-corrosion.md` §4. A 2D grid of film sites with a
breakdown-vs-repair competition; chloride enters via activity, potential via
the `exp(αFE/RT)` factor folded into `constant` rates.

```toml
[deck]
name = "passive-metal-v0"
schema = 2
units = "kcal/mol"

[structure]
kind = "grid"
grid = "square"
neighborhood = "moore"
dims = [128, 128]
boundary = ["periodic", "periodic"]
default_kind = "surf"

[[structure.kinds]]
name = "surf"
initial = "intact"
[[structure.kinds.states]]
name = "intact";      occupant = "film"
[[structure.kinds.states]]
name = "adsorbed";    occupant = "film"
[[structure.kinds.states]]
name = "thinned";     occupant = "film"
[[structure.kinds.states]]
name = "bare";        occupant = "metal"
[[structure.kinds.states]]
name = "repassivated"; occupant = "film"

[dynamics]
[dynamics.thermo]
temperature = 298.15
[dynamics.thermo.activity]
Cl = 0.5

[[dynamics.rules]]
name = "Cl_adsorption"
center = { kind = "surf", state = ["intact"] }
rate = { constant = 1.0e-3 }          # ∝ a_Cl · exp(αFE/RT), folded in
consumes = ["Cl"]
[[dynamics.rules.effects]]
target = "center"; set = "adsorbed"

[[dynamics.rules]]
name = "film_thinning"
center = { kind = "surf", state = ["adsorbed"] }
rate = { arrhenius = { prefactor = 1.0e12, ea = 8.0 } }
[[dynamics.rules.modifiers]]
select = { distance = 1, state = ["thinned", "bare"] }
by_count = { dea = [0.0, -2.0, -4.0] }   # autocatalysis: active neighbors lower Ea
[[dynamics.rules.effects]]
target = "center"; set = "thinned"

[[dynamics.rules]]
name = "repassivation"
center = { kind = "surf", state = ["bare"] }
rate = { arrhenius = { prefactor = 1.0e11, ea = 6.0 } }
[[dynamics.rules.modifiers]]
select = { distance = 1, state = ["thinned", "bare"] }
by_count = { dea = [0.0, 2.0, 4.0] }    # active neighborhood suppresses repair
[[dynamics.rules.effects]]
target = "center"; set = "intact"

[execution]
strategy = "ctmc"
[execution.stop]
steps = 100000
[execution.ensemble]
n_replicas = 1000
seed_policy = "hash"

[observables]
report_every = 1000
[[observables.series]]
kind = "state_counts"
[[observables.series]]
kind = "cluster_sizes"
```

### 7.2 Ice-mantle astrochemistry (Track D)

From `docs/scoping/astrochemistry.md` §5. A 2D grid of adsorption sites with
quenched binding-energy disorder, flux-driven source events, a
temperature-independent tunneling hop, and a Langmuir–Hinshelwood reaction.

```toml
[deck]
name = "ice-mantle-v0"
schema = 2
units = "kcal/mol"

[structure]
kind = "grid"
grid = "square"
neighborhood = "von_neumann"
dims = [100, 100]
boundary = ["periodic", "periodic"]
default_kind = "site"

[[structure.kinds]]
name = "site"
initial = "empty"
[[structure.kinds.states]]
name = "empty"; occupant = "vacant"
[[structure.kinds.states]]
name = "H";      occupant = "H"
[[structure.kinds.states]]
name = "CO";     occupant = "CO"
[[structure.kinds.states]]
name = "HCO";    occupant = "HCO"

# quenched binding-energy disorder via binned kinds (§2.2): two bins
# shown for brevity; production decks use ~8. A seeded init pass assigns
# each site a bin kind by probability, drawing from the init stream.
[[structure.init]]
name = "E_bind_disorder"
center = { kind = "site", state = ["empty"] }
probability = 0.5                      # half the sites become deep-bin sites
map = { empty = "empty" }
rekind = "site_deep"                   # bin kind; per-bin rates via `when`

[dynamics]
[dynamics.thermo]
temperature = 10.0

# flux-driven source event (Track D requirement): no center selector
[[dynamics.rules]]
name = "H_deposition"
rate = { constant = 1.0e-4 }            # F·σ per site
[[dynamics.rules.effects]]
target = "source"                      # NEW effect target: vacant site ← H
set = "H"

# temperature-independent tunneling hop (Track D requirement)
[[dynamics.rules]]
name = "H_tunnel_hop"
center = { kind = "site", state = ["H"] }
rate = { constant = 1.0e-2 }            # tunneling floor rate
[[dynamics.rules.effects]]
target = "center"; set = "empty"
[[dynamics.rules.effects]]
target = "neighbors"
select = { distance = 1, state = ["empty"] }
select_mode = "one"                    # hop to exactly one empty neighbor
set = "H"

# Langmuir–Hinshelwood: neighbor pair (H, CO) → (empty, HCO)
[[dynamics.rules]]
name = "LH_H_CO"
center = { kind = "site", state = ["H"] }
guards = [ { distance = 1, min = 1, state = ["CO"] } ]
rate = { constant = 1.0e-3 }            # computed k incl. Eckart tunneling
[[dynamics.rules.effects]]
target = "center"; set = "empty"
[[dynamics.rules.effects]]
target = "neighbors"
select = { distance = 1, state = ["CO"] }
select_mode = "one"                    # react with exactly one CO partner
set = "HCO"

[execution]
strategy = "ctmc"
[execution.stop]
time = 1.0e6
[execution.ensemble]
n_replicas = 100
seed_policy = "hash"

[observables]
report_every = 1000
[[observables.series]]
kind = "state_counts"
[[observables.series]]
kind = "event_rates"
```

Both fragments compile against the schema in §2. The corrosion deck exercises
the grid shorthand, `by_count` autocatalysis, and the ensemble spec; the ice
deck exercises `target = "source"`, `constant` (tunneling-floor) rates, and
quenched disorder via `[[init]]`. Neither requires a new rate law or a new
strategy — the two named Track D requirements (source events, non-Arrhenius
rates) are covered by `target = "source"` and `rate = { constant = … }`
respectively.

---

## 8. Resolved decisions (recommendations adopted in review; Victor may veto before merge)

1. **`target = "source"` semantics: per-site.** Each vacant site is an
   independent source event at rate `F·σ` — the physically standard flux
   model. It composes directly with the existing Fenwick machinery, keeps
   the determinism contract simple, and makes total deposition rate
   respond correctly to occupancy (as coverage grows, fewer vacant sites,
   lower total flux — for free).
2. **`interface_roughness` axis: explicit `axis` required.** Implicit
   defaults are how decks silently change meaning.
3. **Metropolis `dt` convention: `1/N`** (one MC sweep = N site updates,
   the standard convention).
4. **`seed_policy = "hash"` mix: splitmix64, pinned.** Definition:

   ```
   splitmix64(x): x += 0x9E3779B97F4A7C15  (wrapping)
                  z = x
                  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9  (wrapping)
                  z = (z ^ (z >> 27)) * 0x94D049BB133111EB  (wrapping)
                  return z ^ (z >> 31)

   replica_seed(seed, k) = splitmix64(splitmix64(seed) ^ k)
   ```

   **Test vectors** (any implementation must reproduce these exactly):

   | seed | k | replica_seed |
   |---|---|---|
   | 42 | 0 | `0x57E1FABA65107204` |
   | 42 | 1 | `0xF34FE9248C9342E5` |
   | 42 | 2 | `0x725395388690AE46` |

   The vectors go into the B2 test suite verbatim so implementations
   cannot drift.

---

## 9. Non-goals (explicitly out of scope for this RFC)

- Tau-leaping (named in §3, built later — B7).
- Off-lattice/relaxed positions, long-range electrostatics, explicit solvent
  (already out of scope in DESIGN.md §2).
- A fifth strategy beyond the four named. New strategies are a *new RFC*.
- The `[observables]` implementations themselves (B4), except the four named
  in §6.3 which B4 must deliver.
- The grid shorthand's snapshot-position synthesis (B2 may emit a placeholder;
  real positions are a B4/B5 concern).
