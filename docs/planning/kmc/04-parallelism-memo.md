# 04 — Parallelism Memo

**Standing recommendation: port faithful-serial first, behind clean trait seams.
Treat parallel KMC as an optional later milestone (M10+), never a prerequisite for a
working Rust port.** This memo justifies that from what the code actually does, and
sketches the one parallel approach worth considering later.

---

## 1. The recommendation, up front

1. **Serial, faithful, correct — first.** Reproduce the C++ model in Rust (milestones
   M0–M8 of `03-rust-workspace-design.md`) with a single-threaded loop. This is the
   whole near-term goal. A parallel simulator that gives subtly different physics is
   worse than a serial one that is provably faithful.
2. **Design the seams so parallelism is *addable*, not *required*.** The `Model` trait
   (engine §4) and the flat `SiteGraph` already isolate the two things a parallel
   scheme needs to reason about: which sites an event *reads* and which it *writes*.
   Nothing in the serial design forecloses a later domain-decomposition pass.
3. **Only pursue parallelism if appetite survives contact with reality** (per
   `kmc-port.md`) and only after the serial version is trusted and its performance is
   actually the bottleneck for a run Victor wants to do.

The rest of this memo explains why KMC resists parallelism in general, why *this*
model's structure makes one specific scheme plausible, and what it would cost.

---

## 2. Why KMC is intrinsically serial

Kinetic Monte Carlo advances a single global clock by selecting **one** event at a
time, weighted by all current rates, and advancing time by `dt = -ln(u)/R_total`
(spec A7). The next event's probabilities depend on the state left by the previous
event. This is a strict sequential dependency: event *k+1* cannot be selected until
event *k*'s state mutation is applied. `kmc-port.md` states this plainly — "KMC
serializes at time-step boundaries, baked into the foundational algorithm." That is
correct and it is not a limitation of *this* code; it is the definition of the method.

So you cannot parallelize KMC by "doing several events at once" without changing the
statistics — unless you can prove the events you do concurrently are **independent**
(their read/write footprints don't overlap). That proof is exactly what
domain-decomposition parallel KMC tries to provide.

## 3. What the code actually does — the read/write footprint of an event

This is the crux, and it's favorable to a halo scheme *if* one is ever wanted.

- **Selection is global** (`DoEvent` sums every event's rate). That's the sequential
  part.
- **Each event's *effect* is local.** Every `DoReactionN`, `AdsorbX`, `DesorbX`
  mutates the chosen site and a bounded set of its **direct neighbors** (≤6), plus in
  a couple of cases the `pair` partner of a bridging oxygen. Writes never reach beyond
  1 hop (plus the paired-O, which is itself a near neighbor found at build time).
- **Rate *inputs* reach 2 hops.** `IsActive` for forward hydrolysis inspects
  neighbors-of-neighbors (the surface-reachability test), and `CheckEnv` for the
  300/400/500 classes reads neighbors and, through `pair` and Si→O links, effectively
  2 hops out (spec A6). So an event at site *s* can change the *rate* of events up to
  ~2 hops away.

Net: **write radius ≈ 1 hop (+pair), rate-dependency radius ≈ 2 hops.** This bounded,
short-range coupling is the structural fact that makes spatial decomposition
*conceivable* — and it is precisely the "nth-neighbor halo structure" the
`kmc-port.md` brief points at. Two events whose 2-hop neighborhoods are disjoint are
independent and could, in principle, be executed concurrently without altering the
statistics.

## 4. The one parallel approach worth considering later: sublattice / domain decomposition

If parallelism is ever pursued, the appropriate family is **synchronous sublattice KMC
(the Shim–Amar / Lubachevsky style)**, not thread-level event-list parallelism.

Sketch, specialized to this model:
- Partition the lattice's periodic in-plane direction into strips (domains), one per
  worker. The open (surface-normal) direction is short and shouldn't be split.
- Give each domain a **halo** of width = the rate-dependency radius (~2 hops) mirroring
  neighboring domains' state.
- In a synchronous sub-step, each worker runs KMC **within its domain interior**
  (excluding a boundary buffer of halo width), so no worker's events touch another's
  interior. Interior events are provably independent → run them concurrently.
- Handle boundary-region events with a conservative synchronization rule (e.g. process
  boundary events serially / with rollback, or the standard "leave a dead zone and
  reconcile each sync"). Exchange halos, advance the shared clock, repeat.

Why this and not "parallelize the event list": the event list is rebuilt every step and
selection is inherently global; threading *that* buys almost nothing and risks the
statistics. Domain decomposition instead exploits the **spatial locality of effects**
(§3), which is the only real concurrency this problem has.

### Feasibility notes grounded in the code
- **Favorable:** effects are genuinely local (§3); the graph is static (topology never
  changes — diffusion/cluster-removal that *would* change it is dead, spec B6), so
  domain boundaries and halos are fixed once; the lattice is already a flat array with
  computable neighbor indices, so strip partitioning is arithmetic.
- **Unfavorable / costly:** (a) the halo must cover 2 hops because rates depend 2 hops
  out, not just the 1-hop write radius — a wider, more error-prone halo than a naive
  nearest-neighbor model; (b) the `pair`/`lostal` bridging-oxygen bookkeeping crosses
  sites and must be kept consistent across a boundary — a real correctness hazard; (c)
  synchronous sublattice KMC changes the time-advancement bookkeeping (per-domain vs
  global clock), so it is **not** bitwise-faithful to the serial model and needs its
  own statistical validation; (d) the sample lattices are tiny (e.g. 20×3, 100×10,
  500×10) — at these sizes the serial loop is already fast and parallel overhead would
  dominate. Parallelism only pays for lattices far larger than anything in the repo.

## 5. Cost/benefit and the trigger condition

- **Benefit** is real only for **large** lattices run for **many** steps — sizes beyond
  the committed samples. There is no evidence in the repo that Victor currently needs
  those sizes; the biggest committed dims are 500×10.
- **Cost** is a genuine research-grade effort: correct halo width (2 hops), cross-domain
  `pair` consistency, boundary-event reconciliation, and a *separate* validation regime
  because it's not bitwise-faithful. That is much larger than the serial port itself.
- **Trigger to reconsider:** the serial Rust port exists, is trusted, and a *specific*
  run Victor wants (large lattice, long time) is too slow serially. Absent that
  concrete need, parallelism is speculative and should stay a "someday" milestone.

Cheaper wins to reach *before* any parallel KMC, if speed ever bites:
- **Incremental event list.** The C++ rebuilds the entire event list every step
  (O(N_sites) allocation + full rescan), which the Rust port's reusable `Vec` already
  improves. The bigger win is only re-evaluating events in the ~2-hop neighborhood of
  the site that just changed, keeping a running `ratesum`. This is a *serial* algorithmic
  speedup that can be 10–100× on large lattices and preserves exact statistics — do this
  before parallelism.
- **Better selection structure.** Replace linear cumulative-sum selection with a Fenwick
  tree / grouped rates for O(log N) selection. Also serial, also exact.

Both of the above live behind the same `Model`/engine seam and don't require
decomposition. They are the right next performance step after faithful-serial; parallel
KMC is the step after *that*, if ever.

## 6. Seam requirements so the door stays open

To keep domain decomposition addable later without redesigning the engine, the serial
design should (and, per doc 03, does) ensure:
- The `Model` exposes an event's **footprint** implicitly via `events_at(site)` /
  `apply(site)` being expressed purely in terms of a site and its neighborhood — no
  hidden global mutable state. (Global state would break any decomposition.)
- `pair`/`lostal` live in engine-visible parallel arrays indexed by `SiteId` (doc 03
  §3 option a), so a future partitioner can reason about cross-boundary bridges.
- The RNG is a trait (doc 03 §5); parallel schemes need per-domain independent streams,
  which a trait makes swappable.

None of these cost anything in the serial build; they're just discipline.

---

## 7. Bottom line

The model's short-range, static-topology, nth-neighbor structure makes synchronous
sublattice domain decomposition the *right* parallel approach **if** one is ever
needed — but the 2-hop rate-dependency halo and cross-boundary bridge bookkeeping make
it costly and non-faithful, and the current problem sizes don't demand it. Port
faithful-serial first; add an incremental event list and better selection as the
serial speed milestone; hold parallel KMC as an explicit, optional, appetite-gated
milestone with its own validation plan. Do not redesign the core algorithm before it
runs, correctly, in Rust.
