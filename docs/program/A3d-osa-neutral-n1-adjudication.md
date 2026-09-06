# A3d Osa-neutral n=1 terminal-path adjudication

- observation time: 2026-09-06 01:12 PDT
- executor: `(hermes-custom-build-001; profile=workstation)`
- status: **unverified terminal-decision candidate**
- production-calculator calls in this adjudication: **none**
- decision candidate: **retire the original Osa-neutral n=1 proton assignment**

## Scope

This adjudication is deliberately narrow. It decides whether the original
Osa-neutral n=1 proton assignment is a usable constrained stationary reactant
seed at B3LYP/def2-SVP/DF under the current crystallographic frozen-shell
contract. It does **not** claim that no Si–O–Al neutral reactant basin exists
under every termination model, proton assignment, solvation model, or electronic
structure tier.

## Evidence binding

Every owner-writable external artifact was rehashed at the observation time.
The expected hashes embedded in the A3a/A3b/A3c source chain all matched.

| Evidence | SHA-256 |
|---|---|
| A3a closeout manifest | `88039a9cb5f1db921a4bf49d93b195a07df023275199c4d987592b327943ecec` |
| A3a conditioned receipt | `e42635ffd41a72623fe3c13aceb1dc0f5b51926e2da279cd7a577efe4127e8d0` |
| A3a conditioned XYZ | `279d168cab2f604c5bb9c3eef83e04adbd0d66bff15a2a35237f750c1d087a47` |
| A3b candidate terminal | `a575400e36e6d6dc9addb5d176eb80520233e68e6500f7d86e7ad547c792df30` |
| A3b independently verified terminal | `7734ebee4fc03546ccd7fc5470934d91474d6f7a2308567e63e5188f5649e113` |
| A3b common-dual stage receipt | `a0226c4b260a198934523a7e02b4ce7807b7d9b77137874517bf0c0c57d93fa2` |
| A3b raw endpoint | `cb61c135e5f0e4b2fcb877eee6d9020c67890937fb1eb7ac6f08cc0a93ae3105` |
| A3b projected endpoint | `d0ef248967812ec7461f5581769a0eddc54f17e5d1e1ef1b8d2c0be8dd899bdb` |
| A3c stage receipt | `2f4a95dd3511f7981f33475f443e75c2e9c761cfd321c6a232fcde4f2ff242c6` |
| A3c raw endpoint | `831198e8343fcfdb74e20a639caee459c6c57416f18b914b6f94d715180e26b2` |
| A3c projected endpoint | `9a171cd7657ef65a34baf382e03faf85254a90d5b7116255d385f09aefa2f9ec` |
| A3c candidate terminal | `db1e8cc3730730ecf6886c9054f8c8f75c7ede57d2ea78740256b298adf4bc47` |
| A3c independently verified terminal | `fee49de7528aad175b185da24498234ecd34e36f3d650623d45a8909d40521cb` |
| A3c launch terminal receipt | `52d9ec7f1c4664c058a501aa6b13cd9c643e89d2081125f1b49036799ec096eb` |
| A3c executor log | `6424b5baa3bb3d221ef70f50a8fe0b6d07483b59781f8a4463e407c1fa2eb86d` |
| A3c verifier log | `f8c6a0b090902894d28d6615f40bb7008cd0ac9245fb8f6f35d5c6a59d6f07ad` |

The A3c receipt binds implementation revision
`065e845bd0d5174091561b03f0cdd188dbd2647d`, which exists as a commit in this
repository. The run signature also binds one 100-step budget, no continuation,
no retry, the exact A3a/A3b source hashes above, B3LYP/def2-SVP/DF GPU settings,
and the three owner-bond constraints.

## Independent arithmetic from the persisted evidence

The persisted projected gradient contains 165 Cartesian components over 55 free
atoms. Recomputing the card's metric from those components gives:

| Gate | Recomputed | Limit | Ratio | Verdict |
|---|---:|---:|---:|---|
| projected RMS gradient | `1.619628443366007e-3 Eh/Bohr` | `3.0e-4` | `5.39876147789x` | fail |
| projected max gradient | `1.229258979843806e-2 Eh/Bohr` | `4.5e-4` | `27.3168662188x` | fail |

Recomputing the constrained O–H distances directly from the persisted projected
endpoint gives:

| Owner bond | Target (Å) | Observed (Å) | Absolute residual (Å) |
|---|---:|---:|---:|
| H50–O26 | `0.960015972955185` | `0.960018750663280` | `2.777708094648e-6` |
| H52–O27 | `0.959998360582948` | `0.959993226342271` | `5.134240676830e-6` |
| H57–O29 | `0.960005635592378` | `0.960006381377162` | `7.457847847059e-7` |

The maximum residual is `5.134240676830e-6 Å`, safely below the `1.0e-4 Å`
gate. A nearest-oxygen recomputation assigns the triad exactly as `H50:O26`,
`H52:O27`, and `H57:O29`; the receipt records no owner changes. Thus the
failure is not a broken distance constraint or proton-owner drift. It is the
independent stationarity gate.

## Why geomeTRIC convergence is not the stationary-seed gate

The A3c executor log shows geomeTRIC 1.1.1 declaring convergence at step 91. Its
last line reports `Grad_T = 1.122e-05/3.712e-05`, displacement
`1.632e-04/5.213e-04 Å`, and an energy change of `-6.117e-07 Eh`.

That flag and the card's scientific gate are different measurements:

1. geomeTRIC's `Grad_T` is its optimizer-local gradient after
   `InternalCoordinates.calcGradProj`, evaluated in the live constrained run;
   `calcGradNorm` reduces per-atom Cartesian vector magnitudes across its full
   optimizer coordinate set.
2. The A3 gate starts from the canonical persisted endpoint after exact frozen-
   shell projection, obtains a fresh production gradient, zeros frozen atoms,
   explicitly removes the three constrained bond-normal components, and computes
   a component-wise RMS/max over the 55 free atoms.
3. The card explicitly requires that second, independently recomputed quantity.
   Optimizer termination is necessary operational evidence, not permission to
   promote a scientific stationary seed.

The evidence therefore supports a precise semantic conclusion: geomeTRIC
converged according to its own in-run projected-gradient/displacement/energy
criteria, but the canonical endpoint failed the separately defined and
independently repeated production-gradient acceptance gate. This adjudication
makes no unsupported claim about whether the numerical discrepancy arose from
electronic-state restart behavior, endpoint canonicalization, or another
calculator-level cause; resolving that would require the additional calculator
work this card forbids and is unnecessary to the acceptance decision.

## Decision candidate and parent disposition

**Terminal rejection.** Retire the original Osa-neutral n=1 proton assignment as
lacking a usable constrained stationary basin under the current builder,
frozen-shell, and B3LYP/def2-SVP/DF contract.

The parent A3 ladder should return to `ready`, with this exact Osa-neutral series
closed rather than replayed:

- do not rerun A3a, A3b, or A3c;
- do not launch Osa-neutral n=2–4 from an unavailable n=1 baseline;
- do not emit a barrier, store row, Petra fragment, or surrogate value for this
  rejected series;
- continue A3's other declared site-family/protonation cells only when each has
  its own valid reactant minimum and normal scientific gates.

The existing A3 card already owns those remaining campaigns, so this decision
does not create a separate execution follow-up. A future materially changed
termination, solvation, proton-assignment, or electronic-tier model would be a
new scientific design slice, not an A3a/A3b/A3c continuation.

## Verification status

The arithmetic and bindings above are the executor's candidate analysis. Per the
A3d verification edge, the causal diagnosis and terminal decision remain
**unverified** until a different cold worker starts from the original A3c
complaint and raw receipts, reproduces the gates, and confirms the disposition.
