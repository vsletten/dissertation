# A3i-calc005-si-n1-pilot — implement and run the balanced Si attachment pilot

- status: done
- track: A (geochemistry)
- priority: P1
- machine: workstation
- depends: A3h ✅
- claimed-by: hermes-custom-build-001

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

## Progress

- 2026-09-09 02:22 PDT — The bounded production unit reached its atomic terminal failure at exact execution source `8c6cad747f1dbc6c8fd9486b9bf96d636ad2df55`: the sole `C` r2SCAN-3c optimizer call failed because the DFRKS scanner's nuclear gradients did not converge. `C` spent one 150-step budget with zero retries; `V` and `SiOH4` spent zero calls. No accepted/raw endpoint, immutable generation, Store, CALC-005 value, Petra fragment, or CALCULATIONS value was emitted. The QI2 lease is absent and the terminal now carries the independently checked restoration receipt. (hermes-custom-build-001; profile=workstation)
- 2026-09-09 02:22 PDT — Two context-cold adversarial reviews were reconciled before PR: final implementation head hardens production-only API authority, optimizer-entry crash accounting, post-promotion quarantine, terminal-after-restoration publication, live legacy `Check200`/Si-transition binding, and the verifier's own systemd/QI2 envelope. Fresh closeout is 72 focused tests plus 927 passed/1 skipped full QM tests, whole-QM Ruff/format/compileall, full Petra workspace, and diff checks; no scientific budget was replayed. (hermes-custom-build-001; profile=workstation)
- 2026-09-08 23:59 PDT — Implemented the dedicated fixed-`i=1` CALC-005 minima driver, hash-bound immutable evidence generation, read-only Store validator, optimizer-free verifier, and Petra semantic/refusal gates. A cold adversarial pass found and regression-closed false greens in heavy-atom topology, mode cardinality, Store acceptance, verifier independence, Petra selector binding, source/envelope provenance, zero-retry crash timing, and atomic publication. Pre-production gates pass: 114 focused QM tests; 913 passed/1 skipped full QM tests; whole-QM Ruff/format/compile; focused and full Petra workspace; and `git diff --check`. No production calculator call or scientific value is claimed yet. (hermes-custom-build-001; profile=workstation)

## Result

- 2026-09-09 02:22 PDT — **Terminal outcome: `incomplete-computational-failure`; no scientific value.** The exact authorized `(x=1,y=0), i=1` run failed inside the first and only `C` optimizer call with `RuntimeError: Nuclear gradients of <pyscf.scf.hf.DFRKS_Scanner ...> not converged`. The zero-retry ledger is `C=1`, `V=0`, `SiOH4=0`; no endpoint survived the calculator exception and no generation, Store, Petra fragment, or CALCULATIONS value exists. The final terminal SHA-256 is `b8e9da71bc4a98dd2e07250af4750b4c6aee253a17ac37c8dba78eaa2521c169` (pre-restoration terminal `b9ced238a034f80f9333a7e58d44c8c0b51f2c133d793395a2747a9295c9e259`), executor log SHA-256 is `4f443d7b76bc3ef26cfe4d2f8944c3ea2e4053b07c1fd026435ae9261acc7c91`, QI2 is released, and the bounded unit has no live PID. A3j must independently adjudicate this terminal failure; replay and surrogate publication remain prohibited. (hermes-custom-build-001; profile=workstation)
