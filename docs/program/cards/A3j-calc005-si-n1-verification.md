# A3j-calc005-si-n1-verification — independently verify the Si pilot

- status: ready
- track: A (geochemistry)
- priority: P1
- machine: workstation
- depends: A3i-calc005-si-n1-pilot
- claimed-by:

## Objective

From a context-cold start at CALC-005 and the live Petra adsorption/desorption
reactions, independently verify A3i's raw neutral-200s-Si `i=1`, `(x=1,y=0)`, `Osa.sih -> Osa.albr`
pilot. Use the dedicated optimizer-free verifier to re-establish component and
atom identity, stationary minima, thermodynamic arithmetic, Store provenance,
and the no-emission boundary. Return a typed PASS/FAIL verdict that either
releases a narrowly scoped next decision or terminally rejects the pilot.

## Constraints

- Follow fleet `POLICY.md`, `docs/program/PROTOCOL.md`, and
  `docs/program/A3h-calc005-si-attachment-protocol.md`.
- Do not trust A3i's Result or terminal verdict as evidence. Start from the raw
  artifacts, source/deck/settings hashes, CALC-005, `Check200`, and live Petra
  reactions.
- Zero optimizer calls and zero artifact writes outside the verifier's separate
  atomic verification receipt. Recomputed energies/analytic gradients and
  PHVA/full frequencies may use the GPU under a bounded QI2 envelope; do not
  mutate the executor's Store.
- No retry, threshold relaxation, alternate isomer, `n=2..4`, Petra emission,
  Store republish, or CALCULATIONS update.

## Acceptance

### Verify

- Rehash source, deck, settings, component XYZs, atom map, executor receipt, and
  Store; reject any mismatch or incomplete edge.
- Prove the condensed seed is state 204, reconstruct its one-water state expansion
  to live center state 205, and verify the complete atom bijection and exact
  desorption cycle `Al6H38O30Si -> Al6H34O26 + H4O4Si` with no thermochemical
  water term; verify charge/spin, retained deck nodes, frozen shell, proton
  ownership, finite coordinates and collision gates.
- Independently evaluate accepted endpoint energies/analytic gradients and
  PHVA/full frequencies. Reproduce the projected-gradient and minimum gates
  without invoking an optimizer.
- Recompute every electronic, ZPE, thermal, entropy, 1-bar-to-1-M, and
  solute standard-state term and final `S_10` sum from raw receipts.
- Independently derive the `x+y` selector mapping and detailed-balance equation
  from the live Petra/legacy source. Confirm `(x=1,y=0)` is only one of two
  table-index-1 topologies and no Petra/CALCULATIONS value was emitted.
- Emit one typed atomic verification receipt: `verified-pass` only if every
  check passes; otherwise `verified-fail` with exact failed edges. Update A3i,
  this card, parent A3, PLAN and STATUS in one PR.
- A PASS may recommend the next topology/sensitivity gate, but cannot directly
  authorize the full ladder. Any deferred recommendation must be a real board
  card before A3j is marked done.

## Result

Pending independent adjudication of A3i's terminal optimizer failure. A3i emitted no accepted endpoint or value, so A3j must verify the exact failure receipt, zero-retry ledger, artifact absence, source/legacy/Petra boundary, and restoration without replaying any optimizer or inventing missing raw evidence.

## Progress

- 2026-09-09 04:03 PDT — CONTINUATION: implemented an optimizer-free terminal-
  failure adjudication path and four focused tests, but a context-cold adversarial
  review correctly rejected completion. A3i's executor and this queue worker are
  both `hermes-custom-build-001`; changing only a caller-supplied identity suffix
  does not satisfy the required different-worker edge. The review also found that
  the draft receipt was not bound to the verifier implementation, coordinated
  signature tampering could false-green, and current-source reconstruction was
  not isolated from the pinned executor commit. Card remains READY and parent A3
  remains BLOCKED. A genuinely separate worker must harden those edges and perform
  the final receipt-emitting run; no optimizer/calculator replay is authorized.
  (hermes-custom-build-001; profile=workstation)
