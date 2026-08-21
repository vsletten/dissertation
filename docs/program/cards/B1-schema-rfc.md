# B1-schema-rfc — Deck schema v2 + UpdateStrategy RFC

- status: done
- track: B (platform)
- priority: P0
- machine: any
- depends: —
- claimed-by: hermes cloud worker (agents/B1-schema-rfc)

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

- 2026-08-19 — fable — substantive review on PR #25: four gaps
  (R1 quenched-disorder representation unspecified, R2 missing
  exactly-one effect mode, R3 Metropolis fixed proposal order
  breaks multi-rule ergodicity, R4 SynchronousCA conflict
  resolution unspecified). Fixed directly on the claim branch by
  agreement with Victor (worker had provider-infra trouble:
  xai/grok-4.6 timeouts forced a mid-session model swap). §8
  questions resolved per review recommendations, splitmix64 pinned
  with test vectors. Awaiting Victor's read of §4 + merge.
- 2026-08-18 — card created (fable).
- 2026-08-19 — claimed by hermes cloud worker (agents/B1-schema-rfc).
- 2026-08-19 — RFC written: `petra/docs/RFC-001-DECK-V2.md` (all five
  items + both showcase fragments). Read DESIGN.md, schema.rs, both
  showcase scoping docs, `.claude/learnings/`, and the current engine
  (`engine.rs`, `reaction.rs`, `compile.rs`) before writing. DEVIATION:
  none — the card's five items map one-to-one onto the RFC's sections
  §2–§6. Two Track D requirements (source events, non-Arrhenius rates)
  are covered by `target = "source"` and `rate = { constant = … }`
  rather than new rate laws; quenched disorder via `[[init]]` config
  writes (final table syntax deferred to B2, noted in RFC §2.2.1).

## Result

- **What**: `petra/docs/RFC-001-DECK-V2.md` — deck schema v2 (four
  orthogonal sections), the `UpdateStrategy` trait with Rust signatures,
  the v1 compat shim, the determinism/RNG-stream contract, and the
  B2/B3 migration + test plan with analytic acceptance numbers.
- **Where**: `petra/docs/RFC-001-DECK-V2.md` (new file, ~30 KB).
- **Surprises / decisions for the next worker (B2)**:
  - The CTMC strategy must preserve today's `Engine::step` draw order
    and f64 arithmetic **exactly** — the B2 parity gate is byte-identical
    `Fired` logs, not statistical agreement.
  - `target = "source"` is a new effect target (no center selector);
    RFC §8 Q1 recommends per-site (each vacant site is an independent
    event) — confirm with Victor before B2 implements it.
  - The grid shorthand compiles into the same CSR adjacency as the
    explicit cell — it is sugar, not a second code path.
  - Four open questions for Victor are collected in RFC §8 (source-event
    semantics, roughness axis, Metropolis dt, seed-mix function).
- **Follow-up cards**: none new — B2/B3/B4 already exist and are now
  unblocked (B2 was `blocked: B1`).

## Result (acceptance check)

- [x] RFC exists on main and covers all five numbered items (§2–§6).
- [x] Both showcase fragments present (§7.1 corrosion, §7.2 ice mantle).
- [x] Inbox message announces RFC ready for review (see inbox/).
