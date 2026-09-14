# E3b2-classical-receipt-code-binding — bind runner identity to survey receipts

- status: ready
- track: E (muscovite / the perfect circle)
- priority: P2
- machine: any (synthetic tests only; no scientific compute)
- depends: E3b-periodic-dft-spot-checks
- requested-by: hermes-custom-build-001

## Objective

Close the non-blocking provenance gap identified by E3b's one allowed cold
review. A matched-cell classical survey receipt currently binds the prepared
and runtime artifacts plus the LAMMPS executable, but not the exact committed
Python runner or repository revision that produced it. Add code/revision and
invocation identity to future receipts without replaying the completed E3b
scientific campaign.

## Constraints

- POLICY.md and PROTOCOL.md govern.
- This is receipt provenance/attestation only. It must not change E3b's stored
  scientific numbers or typed `disagrees` verdicts.
- No CP2K, LAMMPS production data, GPU, or live campaign run is required.
- Do not retrofit or rewrite immutable historical receipts.

## Acceptance

- New classical survey receipts bind the runner file SHA-256 and exact Git
  revision, and identify the invocation/operator from explicit runtime inputs
  rather than relying only on a hard-coded module string.
- A verifier rejects a receipt whose runner hash or revision does not match its
  claimed source identity.
- Synthetic tests prove successful binding and fail-closed tampering behavior;
  focused Ruff/format/compile/tests pass.
- Card/PLAN/STATUS bookkeeping lands in the implementation PR.

## Progress

- 2026-09-14 16:40 PDT — `(hermes-custom-build-001; profile=workstation)` —
  Filed from E3b's single cold-review round under POLICY §13. This is explicitly
  non-blocking provenance hardening; E3b's numbers, verdict, and safety fixes
  are complete without replaying scientific compute.
