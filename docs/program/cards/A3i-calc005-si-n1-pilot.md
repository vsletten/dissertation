# A3i-calc005-si-n1-pilot — implement and run the balanced Si attachment pilot

- status: ready
- track: A (geochemistry)
- priority: P1
- machine: workstation
- depends: A3h ✅
- claimed-by:

## Objective

Implement the exact CALC-005 contract in
`docs/program/A3h-calc005-si-attachment-protocol.md` and execute its single
authorized neutral-200s-Si pilot: center site 4, `metal_shells=2`, environment
isomer `(x=1,y=0)`, table rung `i=1`, `Osa.sih -> Osa.albr`, charge 0, spin 0.
Build a state-expanded, graph-provenance matched
`C_10 -> V_10 + Si(OH)4` component set, run the bounded minima-only
thermodynamic protocol, and produce one typed, hash-bound terminal receipt.

This card tests whether the construction is physically and operationally viable.
It does not authorize the alternate `(x=0,y=1)` topology, `n=2..4`, Al,
protonation variants, a kinetic barrier, a Petra fragment, or a CALCULATIONS
value.

## Constraints

- Follow fleet `POLICY.md`, `docs/program/PROTOCOL.md`, and the A3h design.
- No route through `phase2_ladder.py` or `a3_family_campaign.py`; this is a
  dedicated matched-minima driver, not a hydrolysis/TS campaign.
- Direct r2SCAN-3c geometry/minimum tier and
  `wB97M-V/def2-TZVPD+SMD(water)` single points; do not inherit the failed
  HF/STO-3G bridge conditioning path.
- One fresh optimizer budget per component, at most 150 steps, zero retry.
  Persist raw endpoints before gates. Missing/failed components are terminal and
  cannot emit a value.
- Run only as a bounded transient unit using the A3h 12 h / 32 GiB / 4 GiB swap /
  16-core / nice-10 envelope, QI2 lease, durable log, atomic receipt, and verified
  restoration. Never hold a queue-drain lane inline for the calculation.
- Large evidence stays outside Git under the exact A3h evidence root.

## Acceptance

- `qm/quarry/crystal.py` retains typed atom origins/topology masks and the CALC-005
  builder expands the condensed `Al6H36O29Si` state-204 seed with one provenance-
  mapped construction water into live state-205 `C_10`. The matched builder proves
  the complete atom bijection and exact live cycle
  `Al6H38O30Si -> Al6H34O26 + H4O4Si`, neutral singlets, identical retained
  deck/frozen nodes, four center OH groups, and exact hydrolysis-water origins.
- Dedicated CALC-005 builder/arithmetic/driver/Store/verifier paths and the exact
  focused tests named by A3h exist. Distance-only deletion is rejected.
- Every accepted component passes exact identity/order, formula, charge/spin,
  finite/collision/frozen-shell/proton-owner gates; fresh projected RMS/max
  gradients are `<=3.0e-4/4.5e-4 Eh/Bohr`; PHVA/full Hessians have no significant
  imaginary frequency above `30 cm^-1`.
- The receipt exposes separate electronic, ZPE, thermal, entropy, and standard-
  state terms for stoichiometry `[-1 C10,+1 V10,+1 SiOH4]` with no water energy;
  independent recomposition matches within `1e-8 kJ/mol`; it labels `S_10` a
  non-kinetic/non-emittable relative environment stabilization.
- Store evidence graph validation passes with the three thermochemical roles plus
  construction provenance and exact job/hash edges, or the run terminates with
  no canonical value. Petra/CALCULATIONS output is absent.
- Focused and full QM tests, whole-QM Ruff check/format, Python compile, Petra
  round-trip refusal/compile gates, and `git diff --check` pass.
- The terminal receipt records source/settings/deck/atom-map/artifact hashes,
  cgroup/QI2 envelope, zero-retry counts, and successful restoration. Update this
  card, parent A3, PLAN and STATUS in the same PR.

### Verify

A3j performs the required different-worker optimizer-free scientific/KMC
verification from raw artifacts and the live Petra reactions. Until A3j passes,
report only the executor's typed pilot outcome; do not authorize or publish a
CALC-005 value.

## Result

Pending.
