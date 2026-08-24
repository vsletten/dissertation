# Al-neutral sequential hydrolysis acceptance contract

**Rule:** change this matrix first, then make the Phase-1 driver and its tests match.

The al-neutral calculation is not accepted as a concerted complex → hydrolyzed-product event. Live CI-NEB/Sella/quick-IRC evidence resolved a sequential associative mechanism. The normal publication lane requires both elementary saddles to connect through one common intermediate. A separate mechanism-finding lane may publish an unresolved addition only after the bounded gates below prove it is an uphill slide rather than laundering a false saddle.

## Basin vocabulary

`hydrolysis_basin_signature` derives three booleans in this order:

1. Si–Ow bonded (`r < 2.3 Å`)
2. Si–Obr bonded (`r < 2.3 Å`)
3. transferred proton on Obr (`r < 1.25 Å`)

| Basin | Signature | Meaning |
|---|---:|---|
| R | `(False, True, False)` | optimized reactant complex |
| I | `(True, True, True)` | associative pentacoordinate intermediate |
| P | `(True, False, True)` | hydrolyzed product |

## Allowed transition matrix

| Reaction path | Accepted for al-neutral? | Required test/gate |
|---|---:|---|
| R ↔ I | Yes, addition step | `test_sequential_quick_irc_gates_each_elementary_step`; `quick_irc_addition_reason` |
| I ↔ P | Yes, cleavage step | `test_sequential_quick_irc_gates_each_elementary_step`; `quick_irc_cleavage_reason` |
| R ↔ P | **No** | al-neutral routing rejects a concerted result even if the generic hydrolysis gate passes |
| R ↔ R, I ↔ I, P ↔ P | **No** | exact endpoint-set mismatch |
| Any other signature | **No** | exact endpoint-set mismatch |

### Mechanism-finding lane (no accepted R ↔ I saddle)

An al-neutral result may omit `TS_add` only when all of these are true:

1. a bounded proton/hydroxyl conformer sweep finds the lowest sampled I,
   records every seed/result, and proves the selected I has zero imaginary
   modes and the exact associative signature;
2. the endpoint-revised R→I band persists every image, converges its bounded
   pre-relaxation/climbing guess, and a bounded Sella refinement is attempted;
3. the candidate fails the exactly-one-imaginary-mode or R↔I quick-IRC gate and
   is explicitly rejected (never stored as an addition saddle);
4. there is no unbounded retry campaign (TASK-168's approved stop was "no
   attempt #5"); and
5. a separately rebuilt I→P saddle passes every ordinary saddle/frequency/IRC
   gate below.

That result is classified `sequential-associative-uphill-slide-no-resolved-
addition-saddle`. Its profile barrier is `G(TS_cleave) - G(R)`. Any sampled
addition crest is diagnostic/bracketing evidence only, never a TS or Gibbs
barrier.

A basin signature identifies bonding topology, not conformer identity. The addition IRC's R and I endpoints must also match the canonical optimized R and cleavage-derived I using both:

- aligned reactive-core RMSD ≤ **0.30 Å** over Obr, Si, Ow, and the two water hydrogens;
- electronic-energy difference ≤ **5.0 kJ/mol** at the same DFT settings.

`test_basin_equivalence_rejects_same_signature_different_conformer` prevents the three booleans from laundering a discontinuous R → I₂ / I₁ → P profile.

## Saddle gates

Each elementary saddle must satisfy all of the following:

1. Sella reports convergence within its bounded step count.
2. `r(Si–Ow) < 2.6 Å` and the refined saddle has not escaped more than 0.5 Å from its NEB guess.
3. The analytic Hessian has exactly one imaginary mode.
4. Quick-IRC at 0.50 Å or 1.00 Å reaches exactly the allowed endpoint pair.
5. Reused frequency/mode data carry exact geometry and scientific-settings fingerprints; a stale same-sized Hessian is rejected.

The addition search is initialized with the aligned R → I Cartesian reaction
vector masked to the reactive Si/Al/Obr/water core; endpoint terminal-hydroxyl
conformer changes may not choose the unstable direction. The loose CI-NEB peak
is first relaxed with exact Cartesian pins on `r(Si–Ow)` and
`r(Obr–Htransfer)`, allowing every transverse spectator coordinate to minimize.
That preparation must converge at projected `fmax ≤ 0.02 eV Å⁻¹`, preserve both
pins within `10⁻⁶ Å`, and fail closed after 120 BFGS steps. The pins are then
removed before a fresh Cartesian Sella search. CI-NEB remains a bounded guess
generator and fails closed on either unconverged relaxation stage.

## Barrier reporting

For quasi-RRHO Gibbs energies `G`:

- addition local barrier: `G(TS_add) - G(R)`;
- cleavage local barrier: `G(TS_cleave) - G(I)`;
- overall reactant-referenced profile barrier: `max(G(TS_add), G(TS_cleave)) - G(R)`.

The largest local barrier identifies the slow elementary step. The overall profile barrier is the like-for-like value compared with si-neutral **113.048262 kJ/mol (27.019 kcal/mol)**. Both values are emitted because a stabilized intermediate can make them differ.

## Publication transaction

Starting a sequential closeout quarantines any prior `results.json` and `store.sqlite` and writes `run_status.json = running`. A failed gate leaves no canonical result/store and marks the status `failed`. The normal lane publishes only after both saddle/IRC gates, minimum checks, conformer equivalence, thermochemistry, and SQLite closeout succeed. The mechanism-finding lane publishes only after its rejected-addition receipt, selected-I minimum gate, rebuilt cleavage saddle/IRC gate, thermochemistry, and SQLite closeout all succeed. Only then is `run_status.json` changed to `completed`.

Primary regression: `test_sequential_closeout_writes_gated_profile`.
