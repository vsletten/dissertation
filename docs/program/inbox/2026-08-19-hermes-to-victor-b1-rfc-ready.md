# B1 RFC ready for review

From: hermes cloud worker (agents/B1-schema-rfc)
To: victor (and fable, for the board)

`petra/docs/RFC-001-DECK-V2.md` is written and ready for your review.
It is the keystone that gates B2/B3/B4 and both showcase tracks.

What it specifies, in implementable detail:

1. **Deck schema v2** — four orthogonal sections (`[structure]`,
   `[dynamics]`, `[execution]`, `[observables]`), with a simple-grid
   shorthand (square/hex/cubic, Moore/von Neumann) compiled into the
   same CSR adjacency as the explicit cell.
2. **The `UpdateStrategy` trait** — Rust signatures; core owns state +
   matching + RNG + recording, strategies own selection-and-time. Four
   initial strategies sketched (ExactCtmc, SynchronousCA,
   AsyncMetropolis, DiscreteTimePCA); tau-leaping named, not built.
3. **v1 compat shim** — kaolinite.toml / kossel.toml load unchanged,
   implying `strategy = "ctmc"`.
4. **Determinism contract** — per-strategy RNG stream salts + documented
   draw order; same deck + seed + strategy ⇒ bitwise trajectory.
5. **B2/B3 test plan** — byte-identical parity gates for the CTMC
   refactor, and the three conformance decks with analytic numbers
   (Conway glider period 4 / (1,1); Ising T_c = 2.269 via Binder
   cumulant; SIR threshold R₀ = 1).

Both showcase models are demonstrated as ~20-line deck fragments (§7):
the corrosion film patch (grid + by_count autocatalysis + ensemble) and
the ice mantle (source events + tunneling-floor constant rates + quenched
disorder). The two Track D requirements — source events and non-Arrhenius
rates — are covered by `target = "source"` and `rate = { constant = … }`,
no new rate law needed.

Four open questions for you are collected in RFC §8 (source-event
semantics, roughness axis, Metropolis dt convention, seed-mix function).
I recommend: per-site source events, explicit roughness axis, dt = 1/N,
splitmix64 seed mix.

Once you ack, B2 (engine refactor) is unblocked.

- ack: (append yours here)
