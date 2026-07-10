# A C++ Veteran's Tour of kmc-rs

**Audience**: a 25-year C++/Python architect who will *read and review* Rust
(much of it agent-written) far more often than he writes it. Every section
explains one idiom used in this codebase: what the C++ counterpart was, why
Rust works this way, what the idiom buys in *guarantees* — and a
**Reviewer's lens**: what to look for when this pattern crosses your desk in
agent-written code. File:line pointers land on the `[IDIOM]` comment blocks,
which carry the full discussion in place.

Read order: this file top to bottom (≈30–40 minutes), then
`crates/kaolinite/src/state.rs` and `crates/mckaol-cli/src/main.rs` in full.
The `[IDIOM]` blocks are placed at each idiom's *first* use, so reading the
code cold also works.

---

## 1. Orientation: crates instead of translation units

The C++ model is 13 `.cpp` files, a hand-written Makefile, and headers that
anything may include. The Rust port is a **Cargo workspace** of four crates
with a declared, acyclic dependency graph:

```
mckaol-cli → kmc-io → kaolinite → kmc-engine
     └────────┴──────────┘  (cli sees all three)
```

- `Cargo.toml:3` — workspace vs. Makefile, and how the crate graph makes the
  layering rules *compiler-enforced*: `kmc-engine`'s dependency list is
  empty, so it physically cannot come to know about aluminum.
- `crates/mckaol-cli/Cargo.toml:7` — binary vs. library crates.
- `crates/kmc-io/src/lib.rs:48` — `pub use` re-exports: modules are private
  plumbing, the crate root is the API. C++ headers leak; Rust modules hide
  by default and you *choose* the public face.
- `crates/kmc-io/src/lib.rs:31` — `#![warn(missing_docs)]`: the compiler as
  documentation-coverage bot.

**Reviewer's lens.** In Rust the architecture *is* in the manifests: review
`Cargo.toml` diffs like you review interface changes. A new dependency edge
(especially model→IO or engine→anything) is a layering decision, not
plumbing. Agents love adding crates to make an error go away — every new
`[dependencies]` line should have a sentence of justification. Also scan for
`pub` creep: anything `pub` is API forever; agents default to `pub` when the
compiler complains about visibility, where `pub(crate)` was the right answer.

## 2. Ownership: the end of `Dispose*`

Every C++ structure here had a `CreateX`/`DisposeX` pair, and `Lattice` had
*two* paths to `delete[]` (destructor + `DisposeLattice`) coordinated by a
nulling convention. In Rust, a value has exactly one owner; when the owner
goes out of scope, the value is freed; borrowing (`&`, `&mut`) lends access
without transferring ownership, and the compiler proves no borrow outlives
the owner and no `&mut` coexists with any other borrow.

- `crates/kmc-engine/src/graph.rs:65` — `SiteGraph` owns its `Vec<Site>`:
  the whole Create/Dispose ceremony, and its double-free surface, deleted.
- `crates/kaolinite/src/cell.rs:83` — `Vec<CellSite>` vs
  `new CellSite[Npos+1]` + end sentinel: a Vec knows its length, so the C++
  sentinel slot (and its off-by-one invitations) has no job.
- `crates/kmc-io/src/sim.rs:79` — `&mut Scanner`: exclusivity as a
  *guarantee*, not a comment. While a function holds `&mut`, nobody else can
  even look.
- `crates/kaolinite/src/build.rs:185` — the **read-then-write borrow
  dance**, the one place the borrow checker genuinely pushes back in this
  port: you cannot hold a reference into `graph.sites` while mutating
  another element. The fix is always the same: copy the small values you
  need into locals, then mutate by index. This is the compiler pointing at
  a real aliasing question the C++ never had to answer.

Why the language works this way: aliasing XOR mutation is the invariant
that makes both memory safety *and* fearless concurrency (§10) fall out of
one rule. The C++ model is small enough to hold in one head; the point of
the rule is code too big for that — or written by something that doesn't
have a head.

**Reviewer's lens.** `.clone()` density is the smell. Cloning compiles away
borrow-checker friction at runtime cost and, worse, hides design confusion
about who owns what — ask "why does this need its own copy?" of every one.
`RefCell`/`Rc` in single-threaded batch code is the same friction laundered
through runtime checks; in this codebase's domain there is no excuse for
either. And `unsafe`: this workspace has **zero**, which is the correct
number for a batch simulator. If an agent ever introduces one, the house
rule is a `// SAFETY:` comment stating the invariant that makes it sound —
reject any unsafe block arriving without its proof.

## 3. `Option`: the `-1` sentinel, retired

The C++ encodes "no neighbor" as `-1` and guards with `if (nbr >= 0)` at
every use site — remembering is the programmer's job. `Option<SiteId>` moves
the flag into the type: you *cannot* use a maybe-missing neighbor without
writing the `None` case.

- `crates/kmc-engine/src/graph.rs:43` — the core swap, including the honest
  cost accounting (16 bytes vs 4; why nobody cares at 1,560 sites; what the
  fix is if a lattice ever has millions).
- `crates/kaolinite/src/build.rs:346` — `if let Some(i2)` as the compulsory
  `>= 0` guard.
- `crates/kaolinite/src/build.rs:205` — `let ... else`: destructure or bail;
  the happy path stays unindented and *typed* (no Option left downstream).
- `crates/kaolinite/src/build.rs:396` — let-chains (`if let Some(x) = e &&
  cond`), porting `if (i2 >= 0 && state % 100 != 1)` as one condition.
- `crates/kaolinite/src/build.rs:514` — `is_some_and`: "present AND passes
  this test" as a single readable clause.
- `crates/kmc-io/src/scan.rs:117` — `Option<&str>` for stream exhaustion vs
  the C++ stream fail-state you must remember to test.

Note where the sentinel *survives*: `crates/kaolinite/src/cell.rs:19` keeps
raw `n < 0` in `NeighborTemplate` because that is what the *file* contains —
the conversion to `Option` happens once, at the lattice build, and the raw
form never escapes. Sentinels at data boundaries, types everywhere else.

**Reviewer's lens.** `.unwrap()` on an Option is the agent shortcut to flag:
each one converts "the type made me think about absence" back into "crash if
I was wrong". In model/library code demand the `match`/`if let`/`?`; unwrap
belongs in tests. Also watch for sentinel smuggling: an agent that stores
`usize::MAX` as "no value" has re-invented `-1` and defeated the type
system on purpose.

## 4. Newtypes and integer honesty

- `crates/kaolinite/src/state.rs:60` — `struct State(pub i32)`: zero runtime
  cost, but "state used as index" is now a compile error. The C++
  `LatticeSite` is six `int`s with six meanings; the type system can carry
  that meaning instead of the variable names.
- `crates/kmc-engine/src/graph.rs:16` — the deliberate counter-example:
  `SiteId` is a plain type *alias*, because id-ness is exercised at every
  subscript and the newtype tax would exceed its safety dividend there. The
  rule of thumb is recorded at the site.
- `crates/kmc-io/src/sim.rs:22` — explicit widths (`i32`, `i64`): Rust has
  no platform-width `int` to inherit surprises from, and the port pins the
  C++ widths deliberately.

**Reviewer's lens.** Bare `i32`/`usize`/`String` fields whose *names* carry
the meaning ("user_id: String") are the smell; a newtype per identity is the
fix. Inverted check: a newtype with one use and no invariant is ceremony —
the `SiteId` alias above shows the honest middle. Ask which bugs the type
makes impossible; if the answer is "none", it's costume jewelry.

## 5. `match`, exhaustiveness — and why `State` is *not* an enum

- `crates/kaolinite/src/state.rs:113` — `match` vs `switch`: no fallthrough,
  arms are expressions, multiple patterns per arm, and (on enums) missing
  cases are compile errors.
- `crates/kaolinite/src/state.rs:43` — the design-judgment section worth
  reading even if you skip everything else: the "correct" Rust idiom
  (an exhaustive `SiteClass` enum) was **considered and rejected** because
  the legacy model does arithmetic on state codes (`state++` to occupy,
  `state < 400` to classify) and a faithful port must keep that arithmetic
  literal. Idioms serve the domain, not the other way around. The enum is
  the right move for the *reformed* model (M8+), and the newtype is the
  bridge that loses nothing meanwhile.

**Reviewer's lens.** In agent code, `match` catch-alls (`_ => {}`) on
*enums* are where bugs hide — they silently absorb newly added variants;
demand explicit arms unless the default has a stated meaning. Here the
`_` arms on raw `i32` codes are correct (they port C++ `default:`) and each
says so. The general test: every `_` arm should be able to explain itself.

## 6. Errors: values, not control flow

The C++ error strategy is `Myerr::die` — print and `exit(1)` from wherever
the problem was noticed. Rust separates detecting (return `Err`) from
deciding (the binary's edge).

- `crates/kmc-io/src/error.rs:15` — the error enum: structured context, and
  a `Result` return the caller *cannot* ignore.
- `crates/kmc-io/src/error.rs:31` — error chaining via `source()` instead of
  string-flattening.
- `crates/kmc-io/src/scan.rs:47` — `map_err` + `?`: catch-wrap-rethrow with
  no hidden control flow. `?` is sugar for an early `return Err(e.into())`,
  nothing more; there is no unwinding to reason about on the happy path.
- `crates/mckaol-cli/src/main.rs:24` and `:43` — policy at the edge:
  `main() -> ExitCode`, `Box<dyn Error>` for "some error, display it",
  precise enums below.
- `crates/kmc-io/src/fmt.rs:41` — `expect` with the invariant stated:
  production panics must carry their proof.

**Reviewer's lens.** The failure modes to hunt in agent code, in order:
(1) `unwrap()` where an error path exists (the #1 shortcut — grep first);
(2) library code that panics on bad *input* rather than returning `Err`
(policy below the edge); (3) `Box<dyn Error>` in a *library's* public API,
which erases the error taxonomy callers need (fine in binaries, lazy in
libraries); (4) errors flattened to `String` early, destroying the chain.

## 7. Generics and traits: templates with a contract

- `crates/kmc-engine/src/graph.rs:28` — `Site<S>`: same monomorphization
  (zero cost, code per instantiation) as a C++ template, but type-checked
  once against declared bounds, so errors are local instead of
  template-backtrace archaeology. Absence of bounds is documentation too:
  this struct promises it will never compare or print an `S`.
- `crates/kmc-io/src/msi.rs:109` — `W: Write`: "takes an ostream&", as a
  trait bound. The golden test writes to a `Vec<u8>`, the CLI to a `File`,
  and the writer cannot tell — the seam is the trait.
- `crates/kmc-io/src/error.rs:61` — traits vs base classes: `Error` is an
  interface a plain value opts into; there is no hierarchy.
- `crates/mckaol-cli/src/main.rs:43` — `dyn` = *opt-in* dynamic dispatch,
  visible in the type, paid only where declared (contrast: C++ virtual is
  a class-wide property).
- Coming at M4: the one trait that matters here, the engine↔model `Model`
  seam (design doc §4) — the Rust replacement for the C++'s hard-wired
  `envrn.cpp`/`actions.cpp` coupling.

**Reviewer's lens.** Over-abstraction is the agent disease on this axis: a
trait with exactly one implementor in the codebase is usually premature
(std traits like `Write` are exempt — std supplies the other impls). Ask
"what second implementation is planned, and when?" — for `Model`, the answer
exists (another mineral, M-later); for a hypothetical `InputFileReader`
trait over four readers, it wouldn't. Generic parameters follow the same
rule: `fn f(x: impl Into<String>)` on an internal helper is flexibility
nobody ordered.

## 8. Iterators: loops with the intent left in

- `crates/kmc-engine/src/graph.rs:99` — `take_while(...).count()` porting
  `CountNbrs`, quirk included; the chain *is* the sentence.
- `crates/kmc-io/src/cell.rs:127` — `iter().filter().count()` vs
  loop-with-accumulator.
- `crates/kmc-io/src/msi.rs:146` — `iter_mut().enumerate()` where C++
  writes `for (i...) id[i] =`: bounds proven once, write target part of the
  loop's declaration.
- `crates/mckaol-cli/src/main.rs:51` — `args().nth(1)`: no argc arithmetic.

Note what this codebase does *not* do: the structural-build loops in
`build.rs` stay as index loops, because they transliterate C++ whose index
arithmetic **is the specification** (the tiling order defines site ids).
Idiomatic iteration would be a rewrite risk with zero payoff there.

**Reviewer's lens.** Performance is never the question (chains compile to
the same loops); *readability* is. The test: can you still read the chain
aloud as one sentence? Two or three combinators, yes; five deep with a
`fold` and two closures, demand the for-loop back. Also check agent code
for `collect::<Vec<_>>()` used only to iterate again — an allocation as
punctuation.

## 9. Tests as a language feature

- `crates/kmc-io/src/scan.rs:174` — `#[cfg(test)] mod tests`: unit tests in
  the same file, zero setup, compiled out of release builds. The cultural
  consequence matters more than the mechanism: there is no excuse layer
  between "wrote a parser" and "tested it".
- `crates/kmc-io/tests/golden_m3.rs:19` — integration tests are *external
  crates*: they exercise only the public API, like a real downstream user.
  The bitwise golden gate lives there, because faithfulness is a public
  contract, not an implementation detail.
- `crates/kmc-io/src/sim.rs:140` — tests that pin *warts*: if someone
  "fixes" the seed-swallowing bug (spec B2) without meaning to, a test
  fails and says so by name.

**Reviewer's lens.** For agent-written Rust, read tests *before*
implementation: agents produce false-green tests that assert what the code
does rather than what it should do (this port's populate-solid test
originally asserted my own wrong arithmetic — the golden gate caught it,
teaching the lesson twice: derived-from-implementation assertions are
worthless without an independent oracle). Prize tests with an external
truth source — a golden file, a hand computation, a spec number. A module
with no `#[cfg(test)]` block is the first review comment.

## 10. Concurrency: the road deliberately not taken

There are no threads in M0–M3, and none coming soon — KMC is intrinsically
serial (one global clock, each event conditioned on the last; parallelism
memo §2). But the *seams* are cut so the door stays open, and they are
worth reading now because they are Rust's concurrency story in miniature:

- `crates/kmc-engine/src/lib.rs:10` — where parallelism WOULD enter: the
  synchronous-sublattice decomposition (memo §4) would partition
  `SiteGraph`'s flat `Vec` into strips with 2-hop halos, one worker per
  strip. The Rust-specific point: `split_at_mut` hands workers provably
  disjoint `&mut [Site<S>]` slices, and the `Send`/`Sync` marker traits
  make "is it safe to move/share this across threads?" a *compile-time*
  question answered per type. The data races a C++ decomposition debugs at
  runtime are type errors here — that is the actual reason to consider
  Rust for the parallel milestone, and for agentic code generally: the
  compiler holds the line that review capacity can't.
- `crates/kaolinite/src/build.rs:52` — the parallel-arrays layout keeps
  `pair`/`lostal` as flat `Vec`s a future partitioner can reason about
  (memo §6's seam requirement), instead of pointers woven through a fat
  struct.
- The RNG (M4) arrives behind a trait for the same reason: per-domain
  streams are a swap, not a rewrite.

**Reviewer's lens.** When threading eventually appears in agent code, the
smells are: `Arc<Mutex<EverythingImportant>>` (a global lock wearing a
safety vest — correct but serial; ask where the *partition* is);
`unsafe impl Send`/`Sync` (a hand-written promise the compiler couldn't
verify — reject without a SAFETY proof); and locks held across `.await` or
long loops. The absence of `unsafe` in a parallel Rust diff is the point of
using Rust; treat its presence as the review.

## 11. The compiler as a review layer

Two moments in this port where rustc/clippy did reviewer work:

- `crates/kmc-io/src/sim.rs:100` — porting the C++ seed bug *verbatim*,
  rustc flagged `unused_assignments`: the compiler **found the 25-year-old
  bug** (a dead store) that `g++ -Wall -Wextra` never mentioned. The
  `#[allow]` that silences it is load-bearing pedagogy: we preserve the bug
  on purpose, and the annotation is the receipt.
- `crates/kmc-io/src/fmt.rs:53` and `crates/kmc-io/src/msi.rs:146` —
  clippy's nudges (range-contains, needless_range_loop), accepted where
  they improved the sentence.

**Reviewer's lens.** `#[allow(...)]` is the audit trail: every one should
carry a reason on the same lines (like both of ours do). An agent diff that
adds allows to get to green without explanations is hiding findings, not
resolving them. `cargo clippy --all-targets` clean and `#![warn(missing_docs)]`
on are cheap machine reviewers — insist on both before human eyes spend
time.

## 12. The bitwise gate: numeric care ledger

M3's acceptance is byte equality of the generated `start.msi` against the
golden capture from `g++ -O3 -ffast-math`. Everything that had to be exactly
right, and where it lives:

| Care point | Where | Story |
|---|---|---|
| decimal→f32 parsing | `crates/kmc-io/src/scan.rs:154` | Rust `parse::<f32>()` and libstdc++ `num_get` are both correctly rounded — parsed constants match bit-for-bit |
| fill-layer truncation | `crates/kaolinite/src/build.rs:251` | `int = cells × 0.7f` truncates an f32 product; 3×0.7f → 2 but 20×0.7f → exactly 14.0 → 14; any "cleaner" arithmetic shifts slab boundaries |
| coordinate transform | `crates/kmc-io/src/msi.rs:78` | **the hard-won ulp**: the C++ *source* divides, but `-ffast-math` includes `-freciprocal-math`, so the golden *binary* multiplies by folded reciprocals; one atom ("OH20") rounds differently through the 6-digit window. The port follows the binary, not the source — with the receipt written at the site |
| sum association | `crates/kmc-io/src/msi.rs:22` | fast-math *licenses* reassociation but g++ kept written order for these short chains; the gate certifies it |
| no FMA | `crates/kmc-io/src/msi.rs:39` | golden binary is baseline SSE2; Rust never contracts implicitly — both round twice |
| float→text | `crates/kmc-io/src/fmt.rs:1` | `ostream << float` ≡ `printf("%g", (double)f)`; implemented per C17 7.21.6.1 from Rust's correctly-rounded primitives |
| rate math f32/f64 dance | `crates/kmc-io/src/rxn.rs:63` | `f32` divide → f64 `exp` (same glibc libm) → narrow at assignment, per C++ promotion rules — matters at M6, pinned now |

The meta-lesson for reviewing numerics: "faithful to the source code" and
"faithful to the binary that produced the reference" can differ by one ulp,
and only a bitwise oracle tells you which one you achieved.

## 13. Warts ledger (preserved bugs, with receipts)

Every legacy wart this port keeps, marked `WART (spec BN)` in code, full
analysis in mission-control `projects/kmc/02-model-spec.md` Part B:

| Wart | Where preserved | Spec |
|---|---|---|
| doubled `drawbonds` read swallows the seed; `ranseed` always 0 | `crates/kmc-io/src/sim.rs:55` (+ pinning test `:140`) | B2 |
| hard-coded render cell matrix ignores `data.cell` | `crates/kmc-io/src/msi.rs:51` | B5 |
| diffusion + BFS dead code — **not ported**; OOB read documented, not reproduced | `crates/kaolinite/src/reactions.rs:1` | B6, B7 |
| flat rate tables make the environment machinery inert in-sample | `crates/kmc-io/src/rxn.rs` tests | B4 |
| adsorption Δμ selected by `i < 18`, routing Al adsorption through dm_si | `crates/kmc-io/src/rxn.rs:91` | A5.4 |
| "9 F" fluorine rendering hack for Al-OH-Al oxygens | `crates/kmc-io/src/msi.rs` element match | A9.1 |
| `TerminateSurface` pass 1: `type` is always 0, so the 401→404 / 406→408 "Si branch" is unreachable (401 always →406, 406 always →409) |  `crates/kaolinite/src/build.rs:354` | **new** — found during this port, not in spec Part B; reported in TASK-016 Result |
| `FindPairs` can transiently self-pair an oxygen | `crates/kaolinite/src/build.rs:178` (+ test) | new, same report |
| `data.rxn`'s "40100" product token (likely a 405 typo) parsed as-is | `crates/kmc-io/src/rxn.rs` tests | new, same report |
| `IsActive` forward-hydrolysis reads `sites[nbr].nbr[6]` (one past the array); under `-O3` the OOB forces the surface test TRUE for any fully-6-coordinated neighbor | `crates/kaolinite/src/environment.rs` (`inner_surface_active` + the WART note) | **new** — found during the M6 parity chase; the biggest behavioral wart in the port; REFORM_PLAN R1 |

Each of these now has a corresponding entry in `docs/REFORM_PLAN.md` — the
fix, its expected physics effect, its test strategy, and the
corrected-by-default-plus-`--legacy` design. The plan is drafted; **no fix is
implemented** (the sequencing decree: parity first, reform after).

## 14. Dynamics (M4–M6): traits in anger, and the parity chase

M4–M6 add the KMC engine loop and the kaolinite reactions behind it. The new
teaching surface, section by section of the code:

### 14.1 The `Rng` trait and a bit-faithful `ran2`
`crates/kmc-engine/src/rng.rs`. The legacy `ran2` is a free function with
**file-static** state (`static long iy`, `iv[]`) — untestable, unresettable,
unsubstitutable. The port makes it a `struct Ran2` behind an `Rng` trait, and
the win is threefold: the generator's state is a *value* (constructible,
clonable for a side-experiment), the trait is a **seam** (parity `Ran2` today,
a modern PRNG tomorrow, engine unchanged), and the whole thing is pinned by a
20-value **bitwise** test against the C++ stream. Reviewer's lens: this is the
rare "trait with one implementor" that is *not* premature — the second
implementor is named and dated in the design doc, and the trait is why parity
and production won't fork the engine. Note the deliberate deviation receipt:
the trait returns `f32`, not the design doc's `f64`, because the C++
`float ran2()` is 32-bit and every downstream consumer (`eps`, the `r < 0.5`
coin) is a float — bit-parity demands it. Faithful now, `f64` at reform.

### 14.2 The `Model` seam — the one trait that matters
`crates/kmc-engine/src/model.rs`. This is the payoff of the whole crate split
(tour §1): the C++ fused engine and model by direct call (`evtlist.cpp` names
`environment->IsActive`, `actions.DoEvent`); here the lump is cut along a
declared interface. Teaching points: an **associated type** (`type State`)
says "exactly one state type per model" (contrast a type *parameter*
`Model<S>`, which would invite two impls for one model); the generic `step<M,
R>` monomorphizes to zero dispatch cost while the bounds document exactly what
it needs; and `ProposedEvent` is a `Copy` POD in a reused `Vec`, retiring the
C++'s `new EventList()`-per-event intrusive linked list (`EventList : public
Event`, tour would call it the B9 idiom) — both faster and free of the
dangling-`next` surface. Reviewer's lens on the `apply(&mut self, graph,
ev, rng: &mut dyn Rng)` signature: the RNG is **threaded in**, not owned by
the model, because the R4/R9 proton coin must draw from the *same* stream as
the step's `dt`/`eps`, *after* them. A model that owned its own RNG would
silently desync parity — the shared `&mut dyn Rng` makes the single ordered
stream a type-level fact, at the cost of one virtual call per coin flip.

### 14.3 The parity chase: three ways the trajectory tried to diverge
`crates/mckaol-cli/tests/parity_m6.rs` is the M6 gate — 20,000 steps, bitwise.
Getting there taught three lessons a reviewer of numeric/agentic ports should
carry:

1. **Float summation *order* is behavior.** The C++ builds its event list by
   prepending (head = last/highest site) and sums from the head; `step` builds
   a `Vec` ascending and must therefore fold **reversed** to match. `f32`
   addition is not associative — the wrong order shifts `ratesum` by an ulp,
   which flips which event crosses `eps`, and the trajectories part ways a few
   steps later. The fix is one `.rev()`, but *finding* it required the oracle.

2. **A bug you must reproduce, not fix.** The `is_active` `nbr[6]`
   out-of-bounds read (warts ledger, above) makes the golden binary allow far
   more forward-hydrolysis events than a correct loop would. A bounds-clean
   port compiles, passes every unit test, and is *wrong* — it desyncs at step
   0 (180 events vs the golden's 660). The only thing that caught it was
   diffing against an independent oracle, and the only way to reproduce it was
   to reverse-engineer the UB's *effect* ("TRUE whenever a neighbor is fully
   coordinated") from the trajectory itself. Lesson: when the reference is a
   compiled binary, "faithful to the source" and "faithful to the binary" can
   differ by a whole code path, not just an ulp — and the binary is the
   oracle.

3. **The oracle's own quirks are part of the contract.** The trajectory's
   state-hash used a digit-truncated FNV basis in the capture harness; the
   Rust checker reproduces *that* constant, not the textbook one, because a
   fingerprint only has to agree between capture and check. Don't "correct" a
   value the oracle depends on.

Reviewer's lens for agent-written faithful ports generally: **prize the
external oracle above the tests the code ships with.** Every green unit test in
this crate was consistent with a wrong `is_active`; the golden trajectory was
not. A port that only tests itself certifies its own misreadings.

### 14.4 Where `rayon` would enter (memo 04)
`crates/kmc-engine/src/lib.rs` marks it: the per-step event rebuild
(`for s in 0..len { events_at(&graph, s, ..) }`) is a read-only fan-out over
an immutable `&graph` — the textbook `par_iter().flat_map().collect()`, and
the borrow checker already proves the workers share no mutable state. But it
must stay parity-gated: a parallel collect reorders the event list, and §14.3
lesson 1 shows why reordering breaks the `f32` summation the oracle pins. So
parallel rebuild is a *reform-era* move (correct-by-default reorders freely,
`--legacy` keeps serial order); `apply` is intrinsically serial (it mutates
and draws the shared RNG), and only memo 04's sublattice decomposition
parallelizes *that* — deliberately out of scope until the serial port is
trusted.
