# A3b-osa-neutral-n1-proton-microstate-stability — test the requested reactant basin

- status: done
- track: A (geochemistry)
- priority: P1
- machine: workstation
- depends: A3a-reactant-minimum-recovery terminal failure
- claimed-by: hermes-custom-build-001

## Objective

Determine whether the original Osa-neutral n=1 termination-proton assignment is
a metastable B3LYP/def2-SVP/DF reactant basin. A3a's production optimizer
converged geometrically but then transferred `H50:O26->O31` and
`H57:O29->O20`; this is a microstate rejection, not optimizer exhaustion. Run
one bounded, provenance-bound two-order release experiment. Do not launch n=2–4,
a saddle, a barrier, a store, or a Petra publication.

1. Start only from the hash-pinned conditioned geometry
   `279d168cab2f604c5bb9c3eef83e04adbd0d66bff15a2a35237f750c1d087a47`
   in A3a's external evidence root.
2. At B3LYP/def2-SVP/DF, obtain one constrained stationary structure holding
   only the two observed original owner bonds `H50-O26` and `H57-O29`, plus the
   unchanged frozen shell. Preserve every raw endpoint before chemistry gates.
3. Release the two owner constraints in both orders. Every release boundary
   gets a fresh optimizer/Hessian and one predeclared finite budget; neither
   route may reuse an incompatible Hessian or receive an automatic continuation.
4. Classify exactly one terminal outcome. **Recovery** requires a fully
   unconstrained minimum retaining every original proton owner. **No basin**
   requires both release orders to leave the original owner assignment with
   production-level energy/force evidence of the downhill transfer. Anything
   else is an honest inconclusive failure and must produce the next executable
   mechanism card rather than another replay.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one branch,
  branch-matched worktree, and PR.
- GPU-first; OMP/MKL/OpenBLAS threads <=16; nice long CPU stages; tee every long
  run. Use the QI2 GPU lease and bounded dead-man restoration.
- Rehash A3a's closeout manifest and conditioned geometry before use. The
  rejected A3a production geometry was not preserved and must not be
  reconstructed from rounded log text or represented as an input artifact.
- Do not weaken SCF, optimizer, microstate, frozen-shell, collision, finite
  energy/gradient, or PHVA gates to rescue a number.

## Acceptance

- CPU-only regressions prove both constraint-release orders, fresh optimizer
  state at each boundary, exact one-budget bounds, and raw endpoint persistence
  before every chemistry rejection.
- Every constrained/unconstrained endpoint is bound to source, settings, seed,
  geometry hash, optimizer status, finite energy, and exact terminal stage.
- **Recovery outcome:** a fully unconstrained endpoint preserves all original
  owners and passes independent finite projected-gradient thresholds, exact
  frozen-shell/collision gates, and zero significant PHVA imaginary modes.
- **No-basin outcome:** both release orders reproducibly transfer away from the
  original owners, with preserved endpoints and finite production-level
  energy/force evidence. Record that the requested microstate is not a
  B3LYP/def2-SVP/DF basin and keep A3 blocked; do not retry it.
- No n=2–4, saddle, barrier, results store, or Petra fragment is emitted. Full QM
  suite, whole-tree Ruff check/format, CLI dry-run, and `git diff --check` pass.
  A3/PLAN/STATUS are updated in the same PR.

### Verify

A worker other than the executor independently rehashes both release paths from
the conditioned checkpoint, recomputes proton ownership and all final gates,
and confirms the terminal classification from the original A3 user outcome.
Until that edge passes, report only `unverified: recovery candidate` or
`unverified: no-basin candidate`; do not publish a scientific finding.

## Progress

- 2026-09-05 20:47 PDT (hermes-custom-build-001; profile=workstation) — The exact one-budget production experiment terminated honestly at its prerequisite common stage. The dual-owner B3LYP/def2-SVP/DF optimizer exhausted all 100 steps, retained the pinned `H50-O26` and `H57-O29` distances within `3.43e-5 A`, but transferred the previously mobile third proton `H52:O27->O32`; its independently projected RMS/max gradients remained `9.514e-4/9.102e-3 Eh/Bohr`, so no constrained stationary seed existed and neither release route was launched. The independent verifier rehashed source/stage artifacts, recomputed calculator evidence, and confirmed `inconclusive`. Follow-up A3c adds the newly evidenced third owner constraint; no A3b stage may be replayed.
- 2026-09-05 19:00 PDT (hermes-custom-build-001; profile=workstation) — Claimed the hash-pinned A3a conditioned geometry and began the bounded five-stage B3LYP/def2-SVP/DF experiment: one dual-owner constrained stationary seed followed by both single-owner release orders, one fresh optimizer/Hessian and one finite budget per boundary. No n=2–4, saddle, barrier, store, or Petra output is authorized.
- 2026-09-05 17:37 PDT (hermes-custom-build-001; profile=workstation) — Filed from A3a's independently reviewed terminal evidence. The bounded next experiment tests whether the requested termination-proton assignment is a metastable production-level basin; it does not spend another unconstrained convergence continuation.

## Result

**Verified terminal outcome: inconclusive.** The source evidence was rehashed before
execution: A3a closeout manifest `88039a9c...`, conditioned receipt `e42635ff...`,
and conditioned XYZ / geometry `279d168c...`. Under the QI2 GPU lease, the fresh
B3LYP/def2-SVP/DF common-stage optimizer spent its sole 100-step budget. Its raw
endpoint was persisted before gates (`cb61c135...`), and its projected endpoint
was persisted as `d0ef2489...` with finite energy `-4022.867021975885 Eh`.

The two requested owner constraints remained within `1.88e-5 A` and `3.43e-5 A`
of target, frozen-shell projection drift was zero, and the minimum pair distance
was `0.959168 A`. The structure nevertheless was not stationary: projected
RMS/max gradients were `9.514471e-4/9.101941e-3 Eh/Bohr`, respectively `3.17x`
and `20.23x` over the unchanged gates. It also changed the original proton
assignment at `H52:O27->O32`. Because the common constrained stationary seed was
unavailable, both release orders correctly remained `not-run`; no PHVA, saddle,
barrier, result store, n=2-4 calculation, or Petra output was emitted.

External evidence root:
`/mnt/data/vsletten/dissertation-data/a3b-osa-neutral-n1-proton-microstate-stability/`.
Candidate terminal SHA-256: `a575400e36e6d6dc9addb5d176eb80520233e68e6500f7d86e7ad547c792df30`.
Independent verified-terminal SHA-256:
`7734ebee4fc03546ccd7fc5470934d91474d6f7a2308567e63e5188f5649e113`;
the optimizer-free verifier used distinct identity
`hermes-a3b-independent-verifier`, rehashed the source and common-stage receipt,
recomputed the finite energy and projected gradient, and confirmed the
`inconclusive` classification. Executable follow-up
[`A3c-osa-neutral-n1-mobile-proton-triad-conditioning`](A3c-osa-neutral-n1-mobile-proton-triad-conditioning.md)
owns the new three-owner constrained-stationarity mechanism. A3 remains blocked;
this bounded experiment must not be replayed.
