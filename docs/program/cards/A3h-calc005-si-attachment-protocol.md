# A3h-calc005-si-attachment-protocol — define the Si coordination thermodynamic ladder

- status: ready
- track: A (geochemistry)
- priority: P1
- machine: any (design, builder probes, and analytical tests only; no production calculator)
- depends: A3g ✅; parent A3 pre-launch review
- claimed-by:

## Objective

Turn CALC-005's under-specified "Si and Al attachment/detachment energetics vs
coordination" into one executable, independently reviewable **200s Si** pilot
contract before another workstation campaign spends GPU time.

The design must define the balanced atomistic cycle, the exact observable and
reference state, and the mapping from that observable to Petra's attachment /
detachment kinetics. It must then specify the smallest one-rung pilot that can
falsify the construction before authorizing n=2–4. This card produces a design
and deterministic CPU-only proof; it does **not** produce or publish an energy,
barrier, Store row, Petra fragment, or CALCULATIONS value.

## Why this gate exists

Three chemically different neutral bridge families (Osa, OSS, Oaa) spent seven
bounded experiments without reaching a barrier. Their independent verifiers
all found the failure upstream at reactant-basin construction or the fresh
stationarity gate. CALC-005 is the best remaining A3 lane because HANDOFF calls
it minima-only and cheaper, but the repository does not yet define:

- the matched occupied and detached structures;
- the aqueous Si reference species, hydration, protonation/pH, solvent model,
  or standard-state correction;
- whether the computed quantity is a reaction free energy, attachment barrier,
  detachment barrier, or a detailed-balance constraint;
- how the builder's `n_intact` coordinate maps onto Petra's split `Oss.hy`
  `by_count` plus `Osa.si1` override.

The current hydrolysis drivers accept only `oss`, `osa`, and `oaa` bridge
families and require an attacking molecule plus bridge oxygen. They are not a
CALC-005 protocol. Adding `si` to their choice table would be a semantic lie.

## Candidate contract to adjudicate

Use this as the concrete starting hypothesis, not as a pre-approved result:

1. Scope exactly neutral 200s Si at deck center site 4, `metal_shells=2`,
   connectivity n=1..4, charge 0, spin 0. Al and protonation variants remain
   separate later cards.
2. Compare a provenance-matched occupied cluster `C_n` and vacancy cluster
   `V_n` with a balanced cycle such as
   `C_n + n H2O -> V_n + Si(OH)4`. The design must prove atom, charge, spin,
   frozen-shell, and deck-node conservation and may replace this equation only
   with a better fully balanced one.
3. Define a coordination stabilization `D_n` from the accepted component
   minima, with n=0 as an explicit separated-state reference. Never call
   `D_n` a kinetic barrier unless the design supplies and verifies the kinetic
   closure and detailed-balance relation that makes it one.
4. Prefer the smallest representative n=1 pilot. Authorize n=2–4 only after
   that pilot's occupied and vacancy minima, component mapping, arithmetic,
   and independent verifier all pass.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`.
- **No electronic-structure production call.** CPU-only builder probes,
  analytical fake calculators, source inspection, and Petra compile checks are
  allowed.
- Do not replay or reinterpret neutral Osa n=1–4, OSS n=1–4, or Oaa n=2/4/6.
  Do not relax their owner, projected-gradient, PHVA, or zero-retry gates.
- Do not route Si through `phase2_ladder.py` or `a3_family_campaign.py`; define
  a dedicated minima-only driver contract.
- Do not convert a binding/reaction free energy directly into an activation
  barrier without a declared kinetic model and a detailed-balance check.
- Keep every unresolved scientific assumption explicit. The design worker is
  expected to make recommendations rather than create another human gate.

## Acceptance

1. **Physical cycle:** define exact occupied, vacancy, water, and dissolved-Si
   structures; atom/deck-node mapping; stoichiometric coefficients; charge,
   spin, protonation/pH, hydration, solvent, temperature, and standard states.
   State precisely which energy/free-energy terms are included and excluded.
2. **Kinetic semantics:** define whether the output is ΔE, ΔH, ΔG, an
   attachment barrier, a detachment barrier, or a detailed-balance constraint.
   Derive the sign convention and the Petra rate/modifier mapping. Prove the
   emitted forward/reverse pair reproduces the selected equilibrium relation;
   otherwise keep the result out of Petra.
3. **Connectivity mapping:** prove, for every n=0..4 table index, how the QM
   structure maps to the 200s KMC state and to the current `Oss.hy` plus
   `Osa.si1` environment expression. Do not conflate the builder's generic
   `n_intact` with a Petra selector count.
4. **Builder feasibility:** run deterministic CPU-only probes for one fixed Si
   center and proposed n=1..4 structures. Require requested n equals resolved
   n, distinct geometry hashes, expected formulas, finite coordinates,
   collision clearance, peripheral heavy-only frozen shells, and unambiguous
   proton ownership. Specify graph/node provenance needed for a safe matched
   vacancy builder; distance-only atom deletion is rejected.
5. **Numerical protocol:** choose and justify geometry, electronic-structure,
   solvent, and thermochemistry tiers. Define conditioning without blindly
   inheriting the failed HF/STO-3G bridge protocol. Predeclare optimizer budgets,
   projected RMS/max-gradient thresholds, frozen-shell/collision/ownership
   gates, PHVA minimum criteria, checkpoint fingerprints, failure semantics,
   and zero-retry behavior for the pilot.
6. **Evidence contract:** specify raw-endpoint-before-gate persistence, exact
   source/deck/settings/geometry hashes, Store evidence graph, canonical output
   quarantine, atomic terminal receipt, and a different-worker optimizer-free
   verification that recomputes geometry identity, gradients, frequencies,
   thermodynamic arithmetic, connectivity mapping, and typed verdict.
7. **Executable plan:** name exact implementation and test paths, including a
   dedicated CALC-005 driver, matched-pair builder, analytical arithmetic tests,
   orchestration failure tests, Store validator, and Petra round-trip. Include
   focused and full test/lint commands plus a bounded workstation runtime /
   memory / GPU-lease / service-restoration envelope for the later pilot.
8. **Board handoff:** update this card, parent A3, PLAN, and STATUS in one PR.
   A passing design unblocks parent A3 only for the one-rung Si pilot; the full
   n=2–4 ladder and Al family remain unauthorized until their named gates pass.

### Verify

A different worker performs a context-cold scientific-and-KMC review. Starting
from CALC-005 and the live Petra reactions—not from the design's conclusions—it
must independently check the balanced cycle, component identity, sign and units,
detailed-balance derivation, n-to-selector mapping, failure boundaries, and that
the planned pilot cannot emit a value from an incomplete or merely
optimizer-converged state.

## Result

Pending.
