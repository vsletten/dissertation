# CALC-005 neutral-200s-Si attachment protocol

Status: **design complete; no electronic-structure calculation or kinetic value**  
Scope: live kaolinite 200s Si adsorption/desorption, deck center site `4`,
neutral singlet, one `(x=1,y=0)` construction/minimum pilot  
CPU evidence:
[`results/A3h-calc005-builder-probe.json`](results/A3h-calc005-builder-probe.json)

## 1. Decision and quantity

The live reversible reaction is exchange of a **fully hydrolyzed** center Si:

\[
V_{x,y} + \mathrm{Si(OH)_4(aq)} \rightleftharpoons C_{x,y}.
\]

`C` has center state `Si.oh4` (legacy `205`), not a center with `n` intact
Si–O–framework bridges. `V` is the exact post-`desorb-si` support state. Define

\[
\Delta G^\circ_{det,x,y}(T)=G(V_{x,y})+G^\circ[\mathrm{Si(OH)_4(aq)}]-G(C_{x,y}).
\]

Positive `ΔG°det` means the occupied site is stabilized against desorption. The
CALC-005 ladder is the relative environment stabilization

\[
S_{x,y}(T)=\Delta G^\circ_{det,x,y}(T)-\Delta G^\circ_{det,0,0}(T),
\qquad S_{0,0}=0.
\]

The KMC table index is `i=x+y`; it is an **environment index**, not the center
Si's covalent coordination or the generic builder's `n_intact`. The legacy
placeholder `[0,6,12,18,24] kcal/mol` is a relative desorption-barrier modifier,
so `S_i` is the thermodynamic quantity that could constrain it. The absolute
base adsorption/desorption pair remains a separate kinetic calibration.

For the finite-cluster protocol, `(0,0)` is the separated algebraic reference:
`C_00 = Si(OH)4`, `V_00 = ∅`, and `ΔG°det,00=0` by identical-component
cancellation. This deliberately assigns all empty-site/configurational standard-
state physics to the empirical base rate; CALC-005 measures only the incremental
framework stabilization. The one-rung pilot therefore measures `S_10`.

This is a thermodynamic minimum-to-minimum free-energy difference, not `ΔE‡` or
`ΔG‡`. It may not be labeled a barrier or emitted to Petra without §5's kinetic
closure.

## 2. Source-state semantics

The contract is fixed by the live and legacy transition tables:

- center `200 -> 205` on adsorption and `205 -> 200` on desorption;
- `205` means `Si(OH)4` (`legacy/cpp-model/common.hpp`);
- each neighboring state changes under `AdsorbSi`/`DesorbSi`
  (`legacy/cpp-model/actions.cpp`, mirrored by `adsorb-si`/`desorb-si` effects in
  `petra/examples/kaolinite.toml`).

For the three Oss neighbors:

| occupied `C` state | meaning | post-desorption `V` state |
|---|---|---|
| `Oss.hy` / `302` | `Si-OH HO-Si`: center OH plus neighboring-Si OH | `Oss.si1` / `303` |
| `Oss.si1` / `303` | center Si-OH only | `Oss.empty` / `300` |

For the unique Osa neighbor:

| selector class | occupied `C` state | post-desorption `V` state |
|---|---|---|
| `x=0` | `Osa.si1` / `408` | `Osa.empty` / `400` |
| `x=1`, pilot | `Osa.sih` / `402` (`Si-OH HO<Al2`) | `Osa.albr` / `404` |
| `x=1`, other live variants | `Osa.full` / `403` or `Osa.sialh` / `407` | `Osa.alhy` / `405` or `Osa.al1` / `409` |

The pilot fixes `Osa.sih -> Osa.albr`; the other `x=1` states are different
protonation/Al environments and are not aliases that the QM protocol may silently
average.

For an interior 200s site, the selector is exactly

- `x=0` if at least one distance-1 neighbor is `Osa.si1`; otherwise `x=1`;
- `y=count(distance=1, state=Oss.hy)` over the three Oss neighbors;
- `i=x+y`.

This is legacy `Environment::Check200` and the Petra `Oss.hy by_count` plus
`Osa.si1 when(min=1, dea=-6)` implementation.

| i | selector-equivalent environments | allowed interpretation |
|---:|---|---|
| 0 | `(0,0)` | no retained Osa-Al environment and no neighboring Si across a hydrolyzed Oss site |
| 1 | `(1,0)` or `(0,1)` | Al-side environment or one Si-side environment |
| 2 | `(1,1)` or `(0,2)` | Al side + one Si side, or two Si sides |
| 3 | `(1,2)` or `(0,3)` | Al side + two Si sides, or three Si sides |
| 4 | `(1,3)` | Al side + all three Si sides |

Index multiplicities are `i=1: 1+3`, `i=2: 3+3`, `i=3: 3+1`; each `x=1`
class further contains distinct Osa protonation/Al states. No single
`n_intact=i` energy can represent those collapsed KMC classes without an
explicit equivalence study or a split selector.

## 3. Exact pilot structures and atom provenance

### 3.1 Fixed physical state

The pilot uses deck `petra/examples/kaolinite.toml`, center
`Node(4,(0,0,0))`, `metal_shells=2`, `(x=1,y=0)`, Osa transition
`sih(402) -> albr(404)`, charge `0`, and spin `0` (`2S`).

- `C_10`: center `205`; Osa `402`; all three Oss `303`. The central species has
  four Si-OH groups. Formula: **`Al6H38O30Si`**.
- `V_10`: center `200`; Osa `404`; all three Oss `300`. Formula:
  **`Al6H34O26`**.
- dissolved reference: neutral singlet orthosilicic acid
  **`H4O4Si`** in implicit water, nominal pH 7.

The balanced live cycle is therefore exactly

`Al6H38O30Si -> Al6H34O26 + H4O4Si`.

No water is a reactant or product of **desorption**. H2O appears only in the
construction provenance because the condensed crystallographic seed must first
be expanded to the fully hydrolyzed state that the live reaction consumes.
Counting that preparation water again in `ΔGdet` would double-count hydrolysis.

### 3.2 State expansion from the current builder

`from_deck_cell(..., site_kind="Si", n_intact=i)` currently builds a condensed
mineral fragment: center Si has `i` intact Si-O-framework connections and only
`4-i` OH groups, corresponding to center state `205-i`, not live state `205`.
For every retained condensed bridge, the state-aware builder applies

`framework-M-O-Si + H2O -> framework-M-OH + HO-Si`.

The graph mapping is fixed, not distance-selected:

1. the original deck oxygen stays with the framework and receives one named
   water H;
2. the named water O and its other H become the new center Si-OH;
3. terminal deck O/H groups already belonging only to center Si remain center
   ligands;
4. after all `i` expansions, the center Si and its four ligand O/H pairs are a
   bijective `Si(OH)4` subgraph;
5. `V` is produced by moving that exact subgraph out of `C`, not by deleting the
   nearest Si or four nearest oxygens.

For the current deterministic `(x=1,y=i-1)` seed series:

| i | condensed seed (`205-i`) | expansion | required live `C` (`205`) | required matched `V` |
|---:|---|---|---|---|
| 0 | `H4O4Si` (algebraic) | none | `H4O4Si` | empty |
| 1 | `Al6H36O29Si` | `+ H2O` | `Al6H38O30Si` | `Al6H34O26` |
| 2 | `Al6H38O36Si4` | `+ 2 H2O` | `Al6H42O38Si4` | `Al6H38O34Si3` |
| 3 | `Al7H43O46Si7` | `+ 3 H2O` | `Al7H49O49Si7` | `Al7H45O45Si6` |
| 4 | `Al7H43O52Si10` | `+ 4 H2O` | `Al7H51O56Si10` | `Al7H47O52Si9` |

Every row is neutral singlet and conserves atoms in `live C -> V + SiOH4`.
Only the `i=1` row is pilot authority; the rest are analytical test vectors.

### 3.3 Required identity model

`qm/quarry/crystal.py` currently drops `DeckCell Node=(site_index,image)` after
building coordinates. CALC-005 must retain:

- `deck(node)` for every crystallographic heavy atom;
- `termination(parent_oxygen,ordinal)` for existing caps;
- `hydrolysis_water(bridge_node,atom_name)` for each expansion;
- `center_si(Node(4,(0,0,0)))`;
- explicit condensed-edge, hydrolyzed-state, and product-fragment membership.

The returned pair contains condensed seed, live `C`, `V`, dissolved `Si(OH)4`,
construction-water origins, and a complete atom-bijection manifest. The
validator independently checks exact formulas, charge/spin, all state labels,
retained support-node equality, frozen-heavy equality/coordinates, four center
OH groups, and every `C -> V+SiOH4` atom edge. Distance-only deletion,
nearest-O assignment, generic `n_intact` as a state label, or missing water
origin is a hard refusal.

## 4. Standard states and numerical protocol

Fixed choices:

- **Temperature:** `298.15 K`.
- **Dissolved species:** neutral `Si(OH)4(aq)` at 1 M standard state, SMD water.
  With the existing 1-bar ideal-gas thermochemistry convention, apply exactly
  `ΔG°(1M<-1bar) = RT ln(C°RT/p°) = +7.958500693927389 kJ/mol` at
  `T=298.15 K`, `C°=1000 mol m^-3`, `p°=100000 Pa`, once with the
  `+1 Si(OH)4` product coefficient in `ΔG°det`. No bare `Si4+`, `H3SiO4-`,
  generic “Si”, counterion, or pH variant.
- **Surface clusters:** fixed-site unit activity. Use vibration-only quasi-RRHO
  over the movable subspace; exclude whole-cluster translation, rigid rotation,
  and `pV`.
- **H2O:** construction provenance only, with no energy or standard-state term in
  `ΔGdet`. A future hydrolysis-preparation calculation would use pure-liquid
  activity one and must be a separate balanced cycle.
- **Geometry/frequency tier:** exact r2SCAN-3c
  (`r2scan-3c/def2-mtzvpp/d4/gcp`). Direct start from graph-built state-expanded
  seeds; do not inherit HF/STO-3G bridge-conditioning routes.
- **Production electronic tier:**
  `wb97m-v/def2-tzvpd/smd(water)` single point on each accepted minimum.
- **Included:** production electronic difference, r2SCAN-3c ZPE and movable-mode
  quasi-RRHO terms, and the `+7.958500693927389 kJ/mol` solute correction with
  the sign and coefficient above.
- **Excluded from the pilot:** kinetic barrier/tunneling, explicit
  microsolvation/electrolyte, configurational degeneracy, BSSE correction,
  cluster-size convergence, conformer populations, pH variants, and
  coupled-cluster calibration. These exclusions keep the value
  non-production/non-emittable.

One fresh optimizer per `C_10` and `V_10`, maximum 150 steps, zero continuation
and zero retry. `Si(OH)4` receives one fresh full optimization under the same
budget. A predeclared bounded SCF solver sequence may change solver mechanics,
but never method, basis, solvent, geometry, or optimizer count.

Persist exact raw XYZ and calculator/optimizer receipt atomically **before**
applying gates. A failed gate quarantines the endpoint and cannot expose a
canonical value.

Required acceptance:

- independent final projected gradient on `C` and `V`: RMS
  `<=3.0e-4 Eh/Bohr`, maximum `<=4.5e-4 Eh/Bohr`; optimizer-local convergence
  alone never passes;
- frozen coordinates bitwise unchanged after projection;
- exact atom identities/order/map/state/formula/charge/spin/frozen set;
- finite coordinates/energies; minimum pair distance `>=0.75 Å`;
- each intended O-H owner distance `<=1.25 Å` and second-owner margin
  `>=0.15 Å`; no owner change;
- PHVA minima for `C` and `V`, full Hessian minimum for `Si(OH)4`, with no
  significant imaginary mode above the existing `30 cm^-1` noise floor;
- frequency receipts bound to the accepted endpoint and settings.

Record separate `ΔE`, `ΔZPE`, `ΔHthermal`, `-TΔS`, solute standard-state term,
and `S_10` in kJ/mol from raw Hartree/frequency data. The stoichiometry is
`[-1 C10, +1 V10, +1 SiOH4]`; no H2O term exists. Independent recomposition
must match within `1e-8 kJ/mol`. Verdicts are
`passed-protocol-pilot`, `incomplete-computational-failure`, or
`rejected-physical-state`. Only the first may reveal `S_10`, still marked
thermodynamic/non-kinetic/non-emittable.

## 5. Detailed balance and Petra boundary

Let `a_Si` be dimensionless neutral-`Si(OH)4` activity relative to 1 M. For the
absolute reversible pair,

`P(C_xy)/P(V_xy) = a_Si * exp(ΔG°det,xy/RT)`

and

`a_Si*k°att,xy/k°det,xy = a_Si*exp(ΔG°det,xy/RT)`.

Hence `k°det,xy = k°att,xy*exp(-ΔG°det,xy/RT)`. For a common Eyring transition
state,

- `ΔH‡det,xy = ΔH‡att,xy + ΔHdet,xy`,
- `ΔS‡det,xy = ΔS‡att,xy + ΔSdet,xy`.

Relative to the calibrated `(0,0)` base pair,

`(katt,xy/kdet,xy)/(katt,00/kdet,00) = exp(S_xy/RT)`.

If attachment is demonstrably environment-independent, this reduces to
`kdet,xy = kdet,00*exp(-S_xy/RT)`, the sign implied by a positive desorption
barrier modifier.

No CALC-005 number may yet enter the live deck because:

1. `adsorb-si`/`desorb-si` have independently inherited prefactors/barriers and
   the deck species `Si` is not yet the declared neutral `Si(OH)4` standard;
2. `by_count.dea` is energy-only and cannot represent a coordination-dependent
   entropy over multiple temperatures;
3. the current table collapses `(x,y)` and Osa-state variants that are physically
   different;
4. the demo deck temperature is not the fixed 298.15 K pilot convention.

Emission requires an absolute base kinetic anchor, compatible activity/standard-
state declaration, topology resolution, and either separate enthalpy/entropy
support or an explicitly fixed physical temperature. A compile/round-trip test
must recompute the detailed-balance equality for every emitted class. Replacing
legacy `dea` with `S_i`, or calling `S_i` an activation energy, is prohibited.

## 6. CPU-only state-expanded proof

The committed `calc005-live-pair-probe-v1` JSON is generated by the runnable,
calculator-free `qm/scripts/calc005_si_attachment.py` CLI. Two complete passes
materialize both the condensed seed and graph-matched live `C/V` structures for
the deterministic `(x=1,y=i-1)` analytical series `i=1..4`.

It proves:

- every condensed formula expands by exactly `i H2O` to center state 205, while
  every live desorption formula satisfies `C_i = V_i + H4O4Si`;
- the pilot states are exactly
  `205/Osa.sih(402)/Oss.si1(303)^3 ->`
  `200/Osa.albr(404)/Oss.empty(300)^3`;
- the complete C-to-products atom-origin multiset, retained support identities,
  frozen-origin set and frozen coordinates match for every rung;
- validation rebuilds the canonical pair from the bound deck and requires exact
  element order, coordinates, atom origins and condensed-seed identity; it checks
  intended proton owners separately in C, V and Si(OH)4;
- each live C has four center OH groups, and each expansion water contributes
  exactly three provenance-tagged atoms (support H plus center-ligand O/H);
- all live C/V formulas and geometry hashes are distinct; coordinates are finite;
  minimum interatomic distance is `0.96 Å`, above the builder's `0.75 Å` floor;
- frozen atoms are heavy and peripheral; nearest O-H distances are `0.96 Å`;
  the worst first/second-owner margin is `0.38982445786760844 Å`, above `0.15 Å`;
- `qm/quarry/calc005.py`, `qm/quarry/crystal.py`, the probe driver, and deck
  SHA-256 values bind the proof, and both full passes are byte-identical.

The Si(OH)4 fragment is translated deterministically away from the retained
support (`2.8 + 0.6*y Å` along the Osa-to-center direction) solely to create a
collision-free optimizer seed. That offset is not an optimized adsorption
distance and carries no energetic meaning. The proof is a physical-state,
formula, geometry and provenance gate—not a stationary-minimum or energy claim.

Reproduction:

```bash
cd qm
uv run --frozen --extra dev python -m pytest \
  tests/test_crystal.py tests/test_calc005_pair.py \
  tests/test_calc005_arithmetic.py tests/test_calc005_probe_cli.py -q
uv run --frozen python scripts/calc005_si_attachment.py probe \
  --deck ../petra/examples/kaolinite.toml --repeat 2 \
  --output ../docs/program/results/A3h-calc005-builder-probe.json
```

## 7. Evidence graph and independent verification

Every checkpoint fingerprint includes source commit, driver/deck/settings hashes,
center node, exact `C/V` KMC states, `(x,y,i)`, condensed mask, hydrolysis-water
origins, atom-map hash, seed/endpoint hashes, formula, charge/spin, frozen-node
set, optimizer budget and gate thresholds. Resume refuses any mismatch.

Store extensions:

- `calc005_cycles`: schema, `(x,y,i)`, explicit state transition, temperature,
  standard states, status, atom-map and terminal-receipt hashes;
- `calc005_components`: stoichiometric coefficient, role (`C`,`V`,`SiOH4`),
  structure id, and opt/freq/SP job ids;
- `calc005_construction`: condensed seed, each hydrolysis-water atom origin, and
  live-state expansion edges (non-thermochemical);
- ordinary structures/jobs/results preserve geometries, methods and raw terms.

Canonical output is atomically promoted only after complete graph/arithmetic
validation. The read-only validator refuses orphan/missing/duplicate roles,
wrong coefficients/state/formulas, hash/unit drift, non-done jobs, absent raw
endpoint, or kinetic labeling.

A different worker starts from CALC-005 and live Petra/legacy source, not this
Result. It independently rehashes artifacts; reconstructs state transitions,
formulas and the atom bijection; re-evaluates accepted endpoint energies and
analytic gradients; recomputes PHVA/full frequencies and all terms; queries the
Store read-only; derives the selector and detailed-balance equation; and proves
that no Petra/CALCULATIONS output exists. It performs zero optimizer calls and
returns one typed verdict.

## 8. Exact implementation and tests

Implemented by A3h (CPU-only):

- `qm/quarry/crystal.py`: typed atom origins and exact topology-mask support;
- `qm/quarry/calc005.py`: hydrolyzed-state expansion, matched `C/V` transform,
  canonical-rebuild atom/state/formula validator and thermodynamic arithmetic;
- `qm/scripts/calc005_si_attachment.py`: dedicated `probe` CLI;
- `qm/tests/test_calc005_pair.py`, `test_calc005_arithmetic.py`, and
  `test_calc005_probe_cli.py`: exact i=1..4 builds, algebraic n=0 reference,
  balance/detailed-balance arithmetic, adversarial refusal, deterministic CLI;
- `docs/program/results/A3h-calc005-builder-probe.json`: generated proof.

A3i/A3j must add the production surfaces at these exact paths:

- extend `qm/scripts/calc005_si_attachment.py` with `run`/`validate`;
- `qm/scripts/calc005_si_attachment_verify.py`: independent optimizer-free
  verifier;
- `qm/quarry/calc005_store.py`: Store evidence-graph writer and read-only
  validator;
- `qm/tests/test_calc005_driver.py`: stale/incomplete/corrupt checkpoints,
  optimizer exhaustion, owner change, nonminimum/nonfinite/interruption,
  quarantine and zero retry with analytical fake calculators;
- `qm/tests/test_calc005_store.py`: evidence graph and sabotage cases;
- `qm/tests/test_calc005_petra.py`: temporary reversible-pair round trip and
  production-emission refusal under the current incomplete schema.

Gates:

```bash
cd qm
uv run --frozen --extra dev pytest \
  tests/test_crystal.py tests/test_calc005_pair.py \
  tests/test_calc005_arithmetic.py tests/test_calc005_driver.py \
  tests/test_calc005_store.py tests/test_calc005_petra.py -q
uv run --frozen --extra dev ruff check quarry scripts tests
uv run --frozen --extra dev ruff format --check quarry scripts tests
uv run --frozen python -m compileall -q quarry scripts tests
cd ../petra && cargo test -p petra-deck -p petra-core
cd .. && git diff --check
cd qm && uv run --frozen --extra dev pytest -q
cd ../petra && cargo test --workspace
```

## 9. Bounded workstation envelope

Run no DFT inline in a queue-drain session. Launch one bounded transient unit:
`RuntimeMaxSec=12h`, `MemoryMax=32G`, `MemorySwapMax=4G`, `CPUQuota=1600%`,
`Nice=10`, `OMP_NUM_THREADS=16` plus matching BLAS caps. Acquire the QI2
process-bound lease before GPU imports (`expected_gb=16`, maximum `18 GB`, TTL
`13h`, retaining 6 GB Ollama headroom), and run `scripts/gpu_preflight.sh`.

Durable root:
`/mnt/data/vsletten/dissertation-data/a3h-calc005-si-n1-pilot/`; atomic terminal
signal: `terminal-receipt.json`. Default is no shared-service mutation. If a
named service must be isolated, snapshot exact active/enabled state, restore in
an exit trap, arm an independent 13-hour restoration dead-man before the first
stop, and invalidate the receipt until service state and lease release are read
back. Fleet load gates apply. A preflight refusal spends no science budget; an
optimizer/calculator attempt does.

## 10. Authorization boundary

This design authorizes implementation plus one state-expanded
`(x=1,y=0), i=1`, `Osa.sih -> Osa.albr` pilot. It does not authorize the
alternate i=1 class, other Osa states, `i=2..4`, Al, pH/microsolvation variants,
kinetics, Petra emission, or CALCULATIONS publication. Independent pilot
verification is a separate board edge before any next rung.

## 11. Independent review outcome

- 2026-09-06 23:16 PDT (hermes-macbot-zero; profile=laptop) — A different,
  context-cold worker returned **PASS with no blockers** after independently
  reading the live Petra and legacy reactions. It verified the `205->200`,
  `402->404`, `302->303`, and `303->300` transitions; the no-water live cycle;
  exact pilot formula; x/y/index mapping; standard-state sign/value; and detailed
  balance. Its coordinated-tamper matrix covered C/V and condensed provenance,
  element swaps, owner changes, charge/spin/NaN, numeric overflow,
  nondeterminism, and nondistinct rungs and returned PASS. The probe regenerated
  byte-identically. Focused QA passed; the full laptop-compatible suite passed
  with the Linux `/proc/self/exe`-dependent A3g file excluded. Cargo is unavailable
  on this laptop, so the unchanged Petra workspace was not retested here; A3i
  retains the mandatory Petra round-trip gate before any emission.
