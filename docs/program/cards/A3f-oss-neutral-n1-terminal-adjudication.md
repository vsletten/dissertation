# A3f-oss-neutral-n1-terminal-adjudication — close the exhausted H35 owner route

- status: done
- track: A (geochemistry)
- priority: P1
- machine: any
- depends: A3e-oss-neutral-n1-proton-microstate-stability (done)
- claimed-by: hermes-custom-build-001

## Objective

Independently adjudicate the exact, exhausted A3e evidence chain and close the
`oss-neutral-n1-s2` parent disposition without another calculator call. Determine
whether the requested proton assignment satisfies the fixed reactant-minimum
acceptance contract. If it does not, record a terminal rejection of this exact
campaign route and return parent A3 to its other independently gated families.

This is a bounded evidence-only terminal card, not authority for a new scientific
hypothesis. It must not manufacture an absolute claim that no mathematical basin
exists anywhere on the production PES; it decides only whether the requested
assignment produced an acceptable, provenance-backed reactant minimum within the
predeclared finite campaign.

## Evidence contract

- A3e execution source: `333917de87c32b6be4ecc886662dd6666947740c`
- root: `/mnt/data/vsletten/dissertation-data/a3e-oss-neutral-n1-proton-microstate-stability/`
- candidate terminal SHA-256: `b9fa9c4b0940b4e5ff7c40b83ffb106e8b3a4f4f2cbac355030c125eda9ab378`
- independently verified terminal SHA-256: `756e1e885150556ac9b2da4a63ab4e4965246cd1072f83b09a9b5ca6629c86ff`
- independent recomputation SHA-256: `82a5e9d0530568fd5040ffa2ca43fdcd4c6ff799fb94bed110f3faef291316fa`
- owner-conditioning receipt SHA-256: `9d1d725d4261fe35aa8ae9e50d557f8ec0f3e3974aad1b7124817f9bcd5c3778`
- constrained-production receipt SHA-256: `09feaa861412a97ed21e5cc5a4a5ec4baca04cc6fa8581334fefaead09db814c`

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one branch/worktree/PR.
- Evidence-only and CPU-only: zero optimizer, energy, gradient, Hessian, PHVA, or
  other electronic-structure calls. Do not alter external evidence.
- No A3e replay, retry, threshold relaxation, alternate optimizer, surrogate
  barrier, n=2-4 launch, store, Petra fragment, or CALCULATIONS value.
- The fixed A3e release gates remain RMS/max projected gradients
  `<=3.0e-4/4.5e-4 Eh/Bohr` plus owner/residual/shell/collision/finite-energy gates.
- Terminal default: if the exact receipts do not prove every fixed release gate,
  reject this campaign route and release A3 to other families. Do not create an
  A3g retry card from this adjudication.

## Acceptance

1. Rehash every Evidence-contract artifact and prove the independent verifier
   identity differs from the executor identity.
2. From raw receipt arrays—not prose—independently recompute RMS/max projected
   gradients, threshold ratios, H35-O21 residual, owner changes, frozen-shell
   drift, minimum pair distance, and finite-energy predicates for both spent stages.
3. Prove the release stage is absent/not-run, no PHVA was required, each spent
   stage used exactly one 100-step budget, retry count is zero, and no forbidden
   downstream artifact exists.
4. Record the narrow scientific conclusion: the requested `oss-neutral-n1-s2`
   assignment did not produce an accepted reactant minimum under the finite fixed
   contract. Do not promote the earlier HF owner transfer into a production-level
   no-basin claim.
5. Update A3 and PLAN consistently: mark this exact OSS-neutral route terminally
   rejected, preserve the evidence hashes, and return A3 to `ready` for a different
   family without authorizing OSS n=2-4.

### Verify

A different worker performs an adversarial arithmetic/hash pass from the raw stage
receipts and original H35-owner complaint. Until that pass agrees, report only
`unverified: terminal-rejection candidate`.

## Progress

- 2026-09-06 11:46 PDT (hermes-custom-build-001; profile=workstation) — Final isolated CPU-only QA passed with `PYTHONPATH=<worktree>/qm` and `CUDA_VISIBLE_DEVICES=''`: `629 passed, 4 deselected`; whole-QM Ruff check passed, Ruff format reported all 75 files formatted, and `git diff --check` passed. The evidence adjudication itself made zero production electronic-structure, optimizer, gradient, Hessian, or PHVA calls and did not alter the external evidence root.

- 2026-09-06 11:44 PDT (hermes-custom-build-001; profile=workstation) — Independently adjudicated the immutable A3e chain with zero calculator calls. All five contract hashes match; both spent stages used one 100-step budget with zero retries; raw-array arithmetic reproduces owner-conditioning RMS/max `4.853548703717787e-4/2.6190074055318924e-3` and constrained-production `1.0358709556610716e-4/4.836372724428273e-4 Eh/Bohr`. The constrained endpoint passes RMS, H35-O21 residual, owner, shell, collision, and finite-energy gates but fails the fixed max-gradient limit by `1.07474949432x`. Independent worker `a3f-adversarial-verifier` reproduced the hashes, arithmetic, budgets, and downstream absence and returned PASS. Full ledger: `docs/program/A3f-oss-neutral-n1-adjudication.md`.

- 2026-09-06 09:02 PDT (hermes-custom-build-001; profile=workstation) — Created from A3e's independently verified inconclusive outcome to make the parent disposition executable without another compute/retry loop. The terminal default is explicit: failure of any fixed release gate rejects this exhausted campaign route; this card may not create A3g.

## Result

- 2026-09-06 11:44 PDT (hermes-custom-build-001; profile=workstation) — **DONE, independently verified terminal rejection.** The exact `oss-neutral-n1-s2` assignment did not produce an accepted reactant minimum under the finite fixed A3e contract: constrained B3LYP retained `H35:O21` and passed its RMS-gradient (`1.0358709556610716e-4`), residual (`4.201582921470326e-6 A`), shell, collision, and finite-energy gates, but its maximum projected-gradient component `4.836372724428273e-4 Eh/Bohr` exceeded the immutable `4.5e-4` limit. Release/PHVA/downstream artifacts are absent, each spent stage used one 100-step budget, and retries are zero. This is not a global no-basin claim and does not promote the earlier HF owner transfer. A3 returns to ready only for a different independently gated family; A3e replay, threshold relaxation, OSS n=2-4, surrogate barriers/stores/Petra/CALCULATIONS output, and an A3g retry are forbidden. Evidence and the independent PASS are in `docs/program/A3f-oss-neutral-n1-adjudication.md`.