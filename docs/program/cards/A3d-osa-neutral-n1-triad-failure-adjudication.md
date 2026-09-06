# A3d-osa-neutral-n1-triad-failure-adjudication — decide the terminal path

- status: active
- track: A (geochemistry)
- priority: P1
- machine: workstation
- depends: A3c-osa-neutral-n1-mobile-proton-triad-conditioning verified failure
- claimed-by: hermes-custom-build-001

## Objective

Adjudicate the independently verified A3c terminal evidence without another
calculator call. Explain why geomeTRIC reported convergence while the separately
recomputed projected RMS/max gradients remained far above the production gates,
and decide exactly one terminal path for the original Osa-neutral n=1 proton
assignment:

1. retire it as lacking a usable constrained stationary basin and state how the
   parent A3 ladder should proceed or terminate; or
2. authorize one materially different, finite recovery mechanism with an
   executable specification and a new follow-up card.

This is a decision/evidence slice, not another optimization budget.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one branch,
  branch-matched worktree, and PR.
- Read and rehash the external A3a/A3b/A3c receipts and A3c executor/verifier logs
  on the workstation. Do not edit external evidence.
- No electronic-structure, optimizer, Hessian, release, n=2-4, saddle, barrier,
  store, or Petra calculation may run.
- Do not treat geomeTRIC's convergence flag as scientific stationarity; reconcile
  the exact coordinate/force semantics and thresholds from source plus receipts.
- If the decision is another recovery, it must differ materially from the spent
  A3a/A3b/A3c budgets, declare one finite budget, preserve rejected endpoints,
  and remain fail-closed on ownership and independent projected gradients.

## Acceptance

- Rehash and bind the A3c stage receipt, raw/projected endpoints, candidate
  terminal, independently verified terminal, implementation revision, and the
  A3a/A3b source chain.
- Recompute the A3c projected RMS/max gradient metrics and constraint residuals
  from evidence; identify the exact semantic reason geomeTRIC convergence did not
  satisfy the card's stationary-seed gate.
- Record one explicit decision: terminal rejection with parent disposition, or a
  materially distinct finite recovery with an executable follow-up card. No vague
  "try more steps" or replay is acceptable.
- Update A3/PLAN/STATUS in the same PR. Full relevant CPU-only tests, Ruff
  check/format, and `git diff --check` pass.

### Verify

- A different worker performs a cold, read-only verification from the original
  A3c complaint and raw receipts, independently recomputes the numerical gates,
  and confirms that the chosen terminal path follows from the evidence. Until
  that pass lands, any causal diagnosis is labeled `unverified`.

## Progress

- 2026-09-06 01:12 PDT (hermes-custom-build-001; profile=workstation) — Rehashed the complete A3a/A3b/A3c source and terminal chain with every pinned digest matching, rebound implementation commit `065e845bd0d5174091561b03f0cdd188dbd2647d`, and recomputed the A3c metrics from persisted evidence without a calculator call. The 165 free Cartesian components give `1.619628443366007e-3/1.229258979843806e-2 Eh/Bohr` RMS/max (`5.3988x/27.3169x` over gate); direct endpoint geometry gives a maximum owner-bond residual of `5.134240676830e-6 A`, no owner changes, and exact H50:O26/H52:O27/H57:O29 ownership. The unverified candidate decision retires this exact Osa-neutral n=1 assignment and returns A3 to the remaining non-replay campaigns; cold independent verification is pending. Full evidence and convergence-semantics analysis: `docs/program/A3d-osa-neutral-n1-adjudication.md`.
- 2026-09-05 23:46 PDT (hermes-custom-build-001; profile=workstation) — Filed from A3c's independently verified `triad-conditioning failure`: all three original proton owners were retained with residuals below `1e-4 A`, but independently recomputed projected RMS/max gradients were `1.619628e-3/1.229259e-2 Eh/Bohr`, far above the `3.0e-4/4.5e-4` gates despite geomeTRIC reporting convergence. This card owns the evidence-only terminal-path decision; no A3a/A3b/A3c replay is authorized.

## Result
