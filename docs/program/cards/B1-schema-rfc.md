# B1-schema-rfc — Deck schema v2 + UpdateStrategy RFC

- status: ready
- track: B (platform)
- priority: P0
- machine: any
- depends: —
- claimed-by: (see agents/B1-schema-rfc branch)

## Objective

Write the architecture RFC that gates the entire platform
generalization: `petra/docs/RFC-001-DECK-V2.md`. It must specify, in
enough detail that B2's implementer makes no design decisions:

1. **Deck schema v2** as four orthogonal sections:
   - `[structure]`: today's unit-cell crystallography, PLUS simple-grid
     shorthands (square/hex/cubic; Moore / von Neumann neighborhoods).
     Site kinds and state vocabularies carry over.
   - `[dynamics]`: today's reactions generalized to *rules*
     (center/guards/effects/modifiers unchanged in spirit); specify how
     the rate field is interpreted per strategy — propensity (CTMC),
     probability-per-step (PCA), energy delta (Metropolis), truth table
     (deterministic CA).
   - `[execution]`: `strategy = "ctmc"|"synchronous"|"metropolis"|"pca"`,
     strategy parameters, stopping conditions, ensemble spec
     (n_replicas, seed policy).
   - `[observables]`: declarative outputs — state counts, event-rate
     time series, rate spectra, cluster-size distributions, interface
     roughness, snapshot cadence.
2. **The `UpdateStrategy` trait** (Rust signatures): engine core owns
   lattice state, neighborhood queries, rule matching, RNG streams,
   trajectory recording; strategies own selection-and-time. Sketch all
   four initial strategies; name what tau-leaping would need later.
3. **v1 compatibility**: existing decks (kaolinite.toml, kossel.toml)
   load unchanged, implying `strategy = "ctmc"` — spell out the shim.
4. **Determinism contract**: same deck + seed + strategy ⇒ bitwise
   trajectory; how RNG streams are assigned per strategy so this holds.
5. **Migration + test plan for B2/B3**: the parity gates and the three
   conformance decks (Conway glider, Ising/Onsager via Binder cumulant,
   SIR threshold) with their analytic acceptance numbers.

## Context

- docs/OMNIBUS.md Track B; petra/docs/DESIGN.md (read fully — the RFC
  extends it, doesn't fork it); petra/crates/petra-deck/src/schema.rs
  (current schema); docs/scoping/pitting-corrosion.md §"minimal viable
  model" and docs/scoping/astrochemistry.md §"petra angle" (the RFC
  must be sufficient for both showcase decks — check their needs:
  deposition/source events, non-Arrhenius fixed-rate laws, quenched
  per-site disorder are named requirements from Track D).

## Acceptance

- `petra/docs/RFC-001-DECK-V2.md` exists on main and covers all five
  numbered items above.
- Both showcase model sketches (corrosion film patch, ice mantle) are
  expressible in the proposed schema — the RFC demonstrates each in a
  ~20-line deck fragment.
- An inbox message announces the RFC is ready for Victor's review.

## Progress

- 2026-08-18 — card created (fable).
