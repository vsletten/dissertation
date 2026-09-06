# A3f OSS-neutral n=1 terminal adjudication

- observation time: 2026-09-06 11:44 PDT
- executor: `(hermes-custom-build-001; profile=workstation)`
- independent verifier: `a3f-adversarial-verifier`
- status: **independently verified terminal decision**
- production electronic-structure calls in this adjudication: **none**
- decision: **retire the exact `oss-neutral-n1-s2` campaign route**

## Scope

This decision is deliberately narrow. It asks whether the requested
`oss-neutral-n1-s2` proton assignment produced an accepted reactant minimum under
the predeclared finite A3e contract. It does **not** claim that no mathematical
basin exists anywhere on the B3LYP production potential-energy surface, and it
does not promote the earlier HF/STO-3G `H35:O21->O14` transfer into a
production-level no-basin conclusion.

## Evidence binding

All five artifacts named by the A3f Evidence contract were rehashed from the
immutable external root. Every digest matched:

| Evidence | SHA-256 |
|---|---|
| executor candidate terminal | `b9fa9c4b0940b4e5ff7c40b83ffb106e8b3a4f4f2cbac355030c125eda9ab378` |
| independently verified terminal | `756e1e885150556ac9b2da4a63ab4e4965246cd1072f83b09a9b5ca6629c86ff` |
| independent recomputation | `82a5e9d0530568fd5040ffa2ca43fdcd4c6ff799fb94bed110f3faef291316fa` |
| owner-conditioning receipt | `9d1d725d4261fe35aa8ae9e50d557f8ec0f3e3974aad1b7124817f9bcd5c3778` |
| constrained-production receipt | `09feaa861412a97ed21e5cc5a4a5ec4baca04cc6fa8581334fefaead09db814c` |

The receipts bind code revision
`333917de87c32b6be4ecc886662dd6666947740c`, which exists in this repository.
The A3e executor is `hermes-custom-build-001`; its original cold verifier is
`hermes-task292-cold-verifier`, a different identity. A3f's fresh adversarial
pass used the independent identity `a3f-adversarial-verifier`.

## Independent arithmetic from raw receipt arrays

The contract metric uses the 81 Cartesian components belonging to the 27
projected free atoms. Each value below was recomputed from the raw projected
arrays, not copied from status prose.

| Stage / gate | Recomputed | Limit | Ratio | Verdict |
|---|---:|---:|---:|---|
| owner conditioning RMS | `4.853548703717787e-4 Eh/Bohr` | `3.0e-4` | `1.61784956791x` | fail |
| owner conditioning max | `2.6190074055318924e-3 Eh/Bohr` | `4.5e-4` | `5.82001645674x` | fail |
| constrained production RMS | `1.0358709556610716e-4 Eh/Bohr` | `3.0e-4` | `0.345290318554x` | pass |
| constrained production max | `4.836372724428273e-4 Eh/Bohr` | `4.5e-4` | `1.07474949432x` | **fail** |

The literal RMS over all 120 serialized components, including zeroed frozen
components, is `3.987597163464285e-4` for owner conditioning and
`8.51054833620004e-5 Eh/Bohr` for constrained production. That alternate
normalization does not change the decisive constrained-stage result: RMS passes,
maximum component fails.

| Structural/integrity gate | Owner conditioning | Constrained production | Verdict |
|---|---:|---:|---|
| H35-O21 absolute residual | `2.173663363858047e-6 A` | `4.201582921470326e-6 A` | pass (`<=1.0e-4 A`) |
| observed H35 owner | `O21` | `O21` | pass |
| owner changes vs reference | none | none | pass |
| projected frozen-shell drift | `0.0 A` | `0.0 A` | pass |
| raw frozen-shell drift | `1.5331326016276847e-4 A` | `3.710563175651771e-5 A` | projected exactly before acceptance |
| minimum pair distance | `0.9600021762253066 A` | `0.9600042041448641 A` | pass; no collision |
| energy | `-2720.955132238384 Eh` | `-2764.200841973155 Eh` | finite |

Both stage reservations and receipts specify exactly one `max_steps: 100`
budget and `retry_allowed: false`; candidate-terminal records one spend for each
stage and zero retries. The released-production stage has a null receipt and
`status: not-run`; its directory is absent. PHVA is correctly `not-required`.
Searches under the evidence root found no result/store/Petra/PHVA artifact and no
OSS n=2-4 path.

## Decision and parent disposition

**Terminal rejection.** The exact `oss-neutral-n1-s2` route did not produce an
accepted reactant minimum under its finite fixed contract. The constrained B3LYP
endpoint retained the requested proton owner and passed its RMS-gradient,
residual, frozen-shell, collision, and finite-energy gates, but its maximum
projected-gradient component exceeded the immutable limit by
`1.07474949432x`. Release was therefore correctly withheld.

This closes the exhausted A3e route:

- do not replay A3e or relax either gradient threshold;
- do not launch OSS-neutral n=2-4 without an accepted n=1 baseline;
- do not emit a barrier, store row, Petra fragment, PHVA result, or CALCULATIONS
  value for this rejected route;
- return parent A3 to `ready` only for a **different**, independently gated
  site-family/protonation campaign already within the parent card's objective.

No A3g retry card is created. Any materially changed termination, proton
assignment, solvation model, or electronic tier is a new scientific design
slice, not continuation authority from this adjudication.

## Independent verification

`a3f-adversarial-verifier` started from the raw A3e receipts and original
H35-owner complaint, modified no file, and made no scientific calculator call.
It independently matched all five contract hashes; reproduced the 81-component
RMS values, maximum components, ratios, residuals, owner sets, frozen-shell
drifts, minimum pair distances, finite energies, stage budgets, and zero retry
count; and confirmed release/downstream absence. Its verdict was **PASS**: the
terminal rejection follows from the fixed max-gradient gate and remains narrower
than a global no-basin claim.
