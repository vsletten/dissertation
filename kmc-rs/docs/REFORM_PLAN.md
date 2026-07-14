# REFORM_PLAN — correcting the warts, without losing the oracle

**Status: design only. No fix in this document is implemented.** Through M6 the
port is *bug-compatible by decree*: every legacy wart is preserved so the C++'s
accidental determinism (spec B2, fixed seed) can serve as a bitwise parity
oracle. Fixing anything before parity would have destroyed the only ground
truth we have. This plan is the M7+ contract for landing the fixes **without
throwing the oracle away.**

## The governing design: corrected-by-default + `--legacy`

Every reform below ships behind one switch:

- **Corrected-by-default.** A fresh run does the scientifically intended thing:
  reads the seed, uses `f64`, bounds its loops, distinguishes Si from Al where
  the code meant to.
- **`--legacy` mode** reproduces the 2001 behavior bit-for-bit, so the parity
  gate (`tests/parity_m6.rs`) stays green *forever* and the dissertation's
  published results remain reproducible. `--legacy` is not deprecated cruft;
  it is a first-class, permanently-tested mode — the reproducibility guarantee
  for a 25-year-old result.

Mechanically this is a small `RunMode { Legacy, Corrected }` (or a bitflags
`Quirks`) threaded to the few sites that branch: the RNG constructor, the
`is_active` inner loop, the `terminate_surface` state map, the `data.rxn`
parser, and the float width of `rate`/`dt`/`time`. The parity test pins
`Legacy`; a parallel `corrected_*` test suite pins the intended physics
against authored fixtures (a non-flat `data.rxn`, spec B4).

The table is ordered by physics impact — how much fixing each wart moves the
trajectory.

---

## R1. The `is_active` `nbr[6]` out-of-bounds phantom  (NEW — found at M6)

**The wart.** The forward-hydrolysis surface test's inner loop
(`envrn.cpp::IsActive`) is
`for (j=0; (nbr2 = sites[nbr].nbr[j]) >= 0 && j < 6 && !result; j++)`. The
`nbr[j]` read precedes the `j < 6` guard, so when a neighbor has all six
neighbor slots filled the loop evaluates `sites[nbr].nbr[6]` — one past the
array, aliasing the `color` field. Under the golden build (`-O3 -ffast-math`)
this UB manifests, deterministically, as **`result = TRUE` whenever the loop
reaches `j == 6`**. Net effect: forward hydrolysis is allowed at any site with
a fully-coordinated neighbor (all Al neighbors qualify), not only at genuinely
surface-reachable sites. (At step 0 this inflates the event list from a
bounds-clean 180 to the golden 660.)

**The fix.** Bound the loop (`for j in 0..6`, break on the first absent slot),
so the surface test means what A6.1 says: active iff a *real* 2nd-neighbor is
hydrolyzed.

**Expected physics effect.** Large. This is the single biggest behavioral wart
found in the port. It roughly quadruples the set of legal forward-hydrolysis
events, systematically over-activating dissolution away from true surface
sites — it partly defeats the model's central "reactions happen at the
surface" premise. Corrected, the surface-reachability gate becomes the real
differentiator between sites it was designed to be. The corrected trajectory
will diverge from legacy almost immediately and is expected to dissolve more
slowly / more surface-locally.

**Test strategy.** `--legacy` keeps the phantom (parity gate stays green).
Corrected mode gets a hand-built fixture where a 501 site has *no* hydrolyzed
2nd-neighbor but a fully-coordinated Al neighbor: legacy → active, corrected →
inactive. Add an ensemble check that corrected dissolution is surface-biased
(surface roughness / depth profile) vs legacy's bulk-leaking behavior.

**Confidence / open question.** The corrected reading is unambiguous (it is
what the comments and A6.1 describe). Flag for Victor: confirm the intent was
"real 2nd-neighbor hydrolyzed", not some third rule.

---

## R2. The seed swallow  (spec B2)

**The wart.** `Simulation::CreateSimulation` reads `drawbonds` twice,
consuming the seed line into the first read and discarding it; `ranseed` is
never read and defaults to 0, which `ran2` maps to a fixed stream. Every
legacy run shares one PRNG sequence.

**The fix.** Parse all five `data.sim` fields; seed `ran2` from the real value;
allow a CLI/env override for ensembles.

**Expected physics effect.** Total, and the *point* of the reform: a single
fixed trajectory becomes a sampleable ensemble. Individual runs change
entirely; ensemble-averaged observables (population curves, surface roughness)
should be *stable* across seeds — that stability is the validation, replacing
the single-trajectory bitwise diff.

**Test strategy.** `--legacy` forces seed 0 (parity gate). Corrected mode:
(a) different seeds give different trajectories; (b) `initran2(seed)` matches
Numerical Recipes reference values for several seeds (extend the existing
`ran2` bit test beyond seed 0); (c) an ensemble smoke test that mean
populations converge with N runs.

**Confidence.** High; this is a documented, uncontroversial bug. The only care
is keeping the `f32`/`f64` decision (R6) orthogonal, so a seeded legacy-float
run is still possible for A/B comparison.

---

## R3. `TerminateSurface` pass-1 dead Si/Al ternaries  (NEW — found at M3)

**The wart.** Pass 1 computes `type = state % 100`, which is 0 in that block by
construction, then branches on `type == 2` — so the `? 404 : 406` /
`? 408 : 409` arms are dead: 401 **always** demotes to 406, 406 always to 409.
Spec A4.3 describes an intended "depending on whether the missing cation was
Si" distinction the code never had.

**The fix.** Compute the missing-cation class from the *empty neighbor's*
class (the value the author presumably meant), and take the Si vs Al branch
accordingly — restoring the intended state map.

**Expected physics effect.** Medium, and structural (it changes the *initial*
hydroxylated surface, before any dynamics). More chemically faithful capping;
shifts which oxygens start as 404/408 vs 406/409, which then changes early
event availability. Because it moves the *starting configuration*, it also
breaks the M3 `start.msi` bitwise gate — so corrected mode needs its own M3
golden (regenerated from a corrected reference, or asserted against hand
computation), while `--legacy` keeps the current gate.

**Test strategy.** This is the wart most in need of **Victor's domain read**
before implementation (memory: "Victor's read on wart #1 before any reform").
Do not guess the intended map. Once confirmed: legacy M3 gate unchanged;
corrected M3 gets a new golden; a unit test pins the corrected state map arm
by arm.

**Confidence.** Low on *what* the correction should be (the intent is
inferred), high that the current behavior is accidental. **Blocked on Victor.**

---

## R4. The `data.rxn` `40100` product token  (NEW — found at M3)

**The wart.** R12/R13's product is written `40100` in `data.rxn` (reads like a
typo for `405`). Parsed as-is, `reactions[13].reactant = 40100`, which matches
no site's state, so R13 (the reverse of the 404→405 hydrolysis) can never fire.
(Confirmed inert in the oracle: R13 has 0 events across all 20,000 steps.)

**The fix.** This is a **data** fix, not code: correct the token to `405` in a
corrected `data.rxn` (keep the legacy file verbatim for `--legacy`). The
mechanism (`do_reaction13`) is already ported faithfully and correct; it simply
never gets selected today.

**Expected physics effect.** Small-to-medium and *asymmetric*: it restores the
reverse of one hydrolysis channel (404 Al-OH-Al reforming). Legacy, that back-
reaction is missing, so the 404→405 forward is effectively irreversible —
corrected, detailed balance for that channel is restored, nudging the
equilibrium of the Al-OH-Al population.

**Test strategy.** Ship two `data.rxn`: `data.rxn` (corrected, 405) and
`data.rxn.legacy` (40100). Parity gate reads the legacy file. A corrected-mode
test asserts R13 now fires on a 405 site and is the exact inverse of R12
(round-trip a hand-built configuration).

**Confidence.** High that `40100` is a typo; medium that `405` is the intended
value (vs `404`). Cheap to test both against detailed balance.

---

## R5. `Check100` / `Check200` index formulas  (spec B3)

**The wart.** The header documents 100s as `x*3 + y`, 200s as `x*4 + y`; the
code returns `(x+y)/2` and `x+y`. Inert today because rates are flat (B4), so
every bucket holds the same rate — but the moment a non-flat `data.rxn` is
authored (which is the *whole point* of the model), the wrong index selects the
wrong rate.

**The fix.** Implement the header's documented formulas in corrected mode.

**Expected physics effect.** Zero on the sample data (flat rates), potentially
large on any real parameterization — this is latent, not active. It only
"turns on" together with R7 (a non-flat table). Reforming it is a precondition
for the model doing the neighbor-dependent-kinetics science it was built for.

**Test strategy.** Requires the authored non-flat `data.rxn` (see R7).
Unit-test `check100`/`check200` against the header formula on hand-built
neighborhoods; integration-test that adsorption rate now varies with the
counted environment. `--legacy` keeps `(x+y)/2`.

**Confidence.** Medium: the header states an intent, but B3 notes the *counted
quantities* also differ subtly from the header — so "implement the header"
needs Victor to confirm both the formula and the count. **Partially blocked.**

---

## R6. Everything is `f32`; `-ffast-math`  (spec B8)

**The wart.** `rate`, `dt`, and the `time` accumulator are 32-bit; over
millions of steps `time += dt` loses small increments, and `-ffast-math`
relaxes IEEE semantics so results depend on optimization level. The port keeps
`f32` through M6 *specifically* to chase bitwise parity (spec §C2a).

**The fix.** Move `rate`/`dt`/`time` (and the rate-construction math in
`kmc-io::rxn`) to `f64` in corrected mode; keep Kahan-or-just-`f64` time
accumulation.

**Expected physics effect.** Small per step, cumulative over a long run: less
drift in the clock and rate sums, so long-time observables are trustworthy.
Crucially, this **cannot be validated by bitwise diff** against the C++ (it
changes the numbers by construction) — it must be validated statistically,
which is why it is bundled with R2 (seeding → ensembles). This is the reform
that formally retires the bitwise oracle in favor of statistical validation.

**Test strategy.** `--legacy` keeps `f32` (parity gate). Corrected mode:
`f64` throughout; validate that `f32` and `f64` ensembles agree *within
statistical error* on population curves (they should, if the model is not
pathologically ill-conditioned — if they don't, that itself is a finding about
the model's numerical sensitivity worth reporting).

**Confidence.** High on the mechanism, medium on appetite: `f64` is clearly
better, but the `f32`↔`f64` seam touches the most code (every rate/time site),
so it lands last, after the cheaper reforms prove the corrected-mode
machinery.

---

## R7. Author a non-flat `data.rxn`  (spec B4 — data authoring, not a bug fix)

**The gap.** Every rate table in the shipped `data.rxn` is flat (identical
(k, ΔE) across all environment buckets), so `CheckEnv`'s bucket selection
never changes a rate. The environment machinery — the model's scientific
point — is *unvalidated by construction*: no golden artifact can certify
code the reference data never exercises.

**The work.** Author a `data.rxn` with distinct per-bucket rates (physically
motivated values need **Victor's domain read** — this is the M8 objective),
plus a synthetic sharply-non-flat fixture for tests. Fixture runs become the
validation target for R5 and the corrected `CheckEnv` path.

**Expected physics effect.** None on legacy runs (data file, not code). It
*enables* the neighbor-dependent kinetics the model was built to study.

**Test strategy.** With a non-flat fixture: unit-test that `check_env`
bucket selection changes selected rates; integration-test that two sites in
different environments get different propensities. Until R5 lands, corrected
mode still runs the legacy index formulas — R7 is what makes R5 observable.

**Confidence.** High on mechanics; the *values* are Victor's. **Blocked on
Victor for real parameters; unblocked for the synthetic fixture.**

---

## R8. `results.dat`: one appended time series  (spec B1)

**The wart.** `output::initDatafile()` deletes `results.dat` at startup and
nothing ever writes it. `writeData` truncates a fresh one-row file per
snapshot (`step{i}.dat`, `end.dat`) — a 5M-step production run scatters the
population series across thousands of one-row files. The committed
`results.dat` in the legacy tree is a relic of an older version that
appended. The M7 port reproduces all of it faithfully — the scatter *and*
the startup delete — because the golden byte gate demands it
(`tests/golden_m7.rs`).

**The fix.** Corrected mode: open `results.dat` once, append one row per
`wsteps` snapshot (restoring the documented, clearly original intent);
keep `end.dat` or fold it into the series' final row. `--legacy` keeps the
per-file scatter so the golden directory shape stays reproducible.

**Expected physics effect.** None — output plumbing only. Big usability
effect: the population time series becomes one plottable file instead of a
directory-glob-and-concatenate exercise.

**Test strategy.** Corrected mode: run N steps, assert `results.dat` has
`N/wsteps + 1` rows and no `step*.dat` exist; row contents equal the legacy
per-file rows for the same trajectory (same counting code — pinned by
sharing `write_data`'s row serializer between both modes). Legacy mode: the
M7 gate, unchanged.

**Confidence.** High; the C++'s own README documents the appended series.
Low risk; can land any time after the `RunMode` switch exists.

---

## Sequencing for M7+

(M7 landed the output writers bug-compatible and added the directory-shape
gate `tests/golden_m7.rs`; still zero reforms implemented.)

1. **R2 (seed)** first — it is the cheapest, most clearly-correct fix and it
   *creates the ensemble machinery* the statistical validations (R1, R6) need.
2. **R1 (phantom)** next — biggest physics impact, unambiguous correction,
   isolated to one function.
3. **R4 (40100)** — trivial data fix, restores one back-reaction.
4. **R8 (results.dat)** — output plumbing, independent of physics; any time
   after the `RunMode` switch exists.
5. **R7 (non-flat `data.rxn`)** = the M8 objective — *needs Victor's
   physics for real values*; the synthetic fixture can land first. **R3
   (dead ternaries)** and **R5 (Check100/200)** remain *blocked on Victor's
   domain read*; R7's fixture is what makes R5 observable at all.
6. **R6 (f64)** last — the widest seam, and the formal handoff from bitwise to
   statistical validation.

Throughout: `--legacy` stays green against `tests/parity_m6.rs`. The day that
test can no longer pass is the day we have silently broken the reproducibility
promise — so it is load-bearing forever, not a milestone artifact.
