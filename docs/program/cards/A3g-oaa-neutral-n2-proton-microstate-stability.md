# A3g-oaa-neutral-n2-proton-microstate-stability — test the H52 owner basin

- status: done
- track: A (geochemistry)
- priority: P1
- machine: workstation
- depends: A3 Oaa-neutral n=2 advisory failure recorded and supervisor merged
- claimed-by: hermes-custom-build-001

## Objective

Determine whether the exact `oaa-neutral-n2-s2` reactant proton assignment has a
usable B3LYP/def2-SVP/DF stationary basin. The parent A3 campaign's bounded
supervisor stopped correctly after its HF/STO-3G advisory preoptimization spent
100 steps and the unpersisted endpoint made proton `H52` ambiguous between builder
owner `O15` and `O9` (`1.134/1.195 A`). No B3LYP production optimizer ran.

Execute exactly one materially different, finite `H52-O15` owner-constrained →
released production-PES experiment for this single cell. Produce a hash-bound
typed verdict: accepted reactant minimum, verified lower-energy alternative
microstate, or inconclusive terminal failure. This card is not authority to replay
the unconstrained advisory route.

## Evidence contract

Rehash before use and bind every derived receipt to exact branch/source, atom order,
charge/spin, frozen shell, method settings, seed geometry, and evidence-root identity:

- root: `/mnt/data/vsletten/dissertation-data/task274-a3-oaa-neutral-family-20260906/`
- family terminal SHA-256: `e6cfb7063d1dc6bdfe9ce88324183f54838560aa379673a3c1ca2ce20fb531e6`
- family progress SHA-256: `2d75c859be5fe2496432537bfb7e9892b5ccf228377c762caad98d050457a92b`
- child log SHA-256: `ff607898882e7f16c1fd6d52f1283c01646ab32e9dccda603bc80ef54af746a2`
- exact seed `complex_guess.xyz` SHA-256: `bcf597e78dede88afe0f2b7af360c71ebb90b1fbacbb476eeb02ba6365dfbedd`
- metadata SHA-256: `6fdd24d73fb1d2cd5dd6145345562c5803abcf13e7286d57ae5bf6865a2a6190`
- launch/restoration SHA-256: `48bb6e018372caf56828d37aeb18915b2a579a58eafd62a6e75bb2bd6d5821a2` / `ad5a78f88fc84be7ba794b274af431d35b1d97da8cab855fccc054d9755d9646`
- execution source before closeout merge: `97cea5b585c95f18e499f53317ca8864a4c542d2`

The rejected advisory endpoint does not exist. Do not invent its geometry, hash,
energy, gradient, or a completed `H52->O9` transfer from the log.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one branch,
  branch-matched worktree, and one PR carrying implementation, evidence summary,
  card, PLAN, and STATUS updates.
- Exactly one `oaa-neutral-n2-s2` experiment. No n=4/6, saddle, barrier,
  thermochemical barrier, provenance store, Petra emission, or family publication.
- No identical replay of the unconstrained HF/STO-3G advisory preoptimization.
- Use the canonical GPU lease; threads <=16, nice >=10, durable tee logs, finite
  cgroup/runtime ceiling, atomic terminal receipt, and restoration dead-man.
- Persist every raw endpoint and its finite energy/gradient/settings receipt before
  applying a structural, convergence, or ownership gate. Rejected evidence survives.
- Never promote optimizer-local convergence alone. Independently recompute projected
  gradients, constraint residuals, ownership, frozen shell, collisions, and PHVA.
- Zero retries. Any non-accepted terminal outcome closes this exact n=2 route and
  therefore the Oaa-neutral n=2/4/6 serial family. Do not create an A3h retry card.

## Acceptance

1. Rehash the full Evidence contract before launch; mismatch stops without a
   calculator call.
2. Starting from the exact seed, constrain `H52-O15` at its seed distance and run
   one owner-preserving conditioning stage, then one fresh constrained
   B3LYP/def2-SVP/DF optimization bounded at 100 steps.
3. Release is allowed only when the constrained production endpoint has orthogonal
   projected RMS/max gradients <=`3.0e-4`/`4.5e-4 Eh/Bohr`, `H52-O15` residual
   <=`1.0e-4 A`, every original proton assigned to exactly one oxygen within
   `1.25 A` with owner margin >=`0.15 A`, exact atom order/charge/spin/frozen shell,
   minimum pair distance >`0.60 A`, and finite energy/gradient receipts. Release
   uses one fresh B3LYP optimizer bounded at 100 steps with no inherited Hessian.
4. Record exactly one typed terminal outcome:
   - **accepted reactant minimum:** the released endpoint retains every original
     proton owner and passes fresh projected-gradient, frozen-shell, collision,
     finite-energy, and PHVA zero-significant-imaginary-mode gates;
   - **verified alternative microstate:** the released endpoint converges lower in
     finite B3LYP energy to unambiguous changed ownership, with independently
     reproduced endpoint forces/energy and every non-owner integrity gate intact;
   - **inconclusive terminal failure:** any other result, including exhaustion,
     ambiguity, or unavailable release. It authorizes no replay or surrogate output.
5. Raw/projected endpoints, settings, stage receipts, terminal verdict, and complete
   hashes are durable outside Git; source contains only payload-free provenance and
   the scientific conclusion.
6. Full relevant CPU-only tests, whole-QM Ruff check/format, and `git diff --check`
   pass. An accepted minimum may unblock parent A3 only from that exact promoted
   artifact. Either other outcome terminally closes Oaa-neutral n=2/4/6 and returns
   parent A3 only to a different independently gated family.

### Verify

A different worker performs a cold, optimizer-free verification from the original
H52 owner complaint and raw receipts. It independently recomputes ownership,
constraint residual, projected-gradient, finite-energy, frozen-shell/collision, and
PHVA gates and confirms the typed terminal outcome. Before that verification lands,
report only `unverified: candidate outcome`.

## Progress

- 2026-09-06 15:28 PDT (hermes-custom-build-001; profile=workstation) — Filed from the fail-closed Oaa-neutral family attempt after an independent evidence verifier confirmed the n=2 supervisor failure and a separate cold scientific review rejected both an identical replay and a terminal inference from the missing HF-only endpoint. This is the only authorized next Oaa-neutral n=2 calculation: one finite H52-O15-constrained → released B3LYP basin test with raw-endpoint persistence and independent verification.

### Implementation evidence contract (2026-09-06 repair)

Verification requires `--worktree` to identify the exact clean branch revision
matching its remote-tracking ref. The verifier independently pins the audited
executor and pipeline sources and compares recorded loaded optimizer/stage
instructions against code compiled from that revision. These paths implement one
bounded call per stage, a new SCF/kernel invocation, and no Hessian input or retry.
Changing either audited file requires review and an explicit verifier pin update.

The experiment reservation binds executor identity, executable/source/loaded-code
hashes, hostname, PID, boot ID and process start ticks. Per-stage profiling records
optimizer and kernel entries; requested step/Hessian settings are explicitly named
as requests. The API does not return an iteration count or Hessian object identity.
These are local process records and independent code checks, not cryptographic
attestation against a writer with control of the entire machine or receipt store.

`--verifier-identity` supplies only durable `worker_id` and `profile`, with optional
expected `hostname` and `implementation_sha256`. Process facts supplied by the
caller are rejected. The verifier internally observes PID, boot ID, process start
ticks, hostname, executable hash and verifier implementation hash, and embeds them
as `verifier_provenance`. Both the worker ID and the process identity must differ
from the executor reservation. A fresh process on the same host is allowed;
renaming the worker in the executor process is rejected. This does not establish
independent human operators: worker/profile labels still depend on operator
honesty, and a host or evidence-store owner can alter code or forge reservations.
Local process observation is not remote attestation.

Energy and gradient returns are persisted separately before subsequent operations;
raw PHVA includes the complete returned object before gates. Nonfinite numbers are
retained as explicit strings in raw evidence. Executor replay writes only a
separate `replay-attempts/` receipt. Every rejected verifier invocation writes a
unique, crash-durable `verification-attempts/` receipt. Canonical candidate and
verified terminals are never quarantined or overwritten by verifier rejection;
an existing verified terminal rejects replay before source checks or calculators.

Executor and verifier independently define the same recursive forbidden-name
contract: `results.json`, `store.sqlite` and its WAL/SHM/journal sidecars,
`store.task168.tmp.sqlite`, `store.sequential.tmp.sqlite`, `ts.xyz`, `barrier.json`,
`petra.toml`, `family-progress.json`, `terminal-receipt.json`, and `terminal.json`.
No evidence or attempt subtree is exempt. The executor inventories before
reservation and candidate publication; dirty roots yield a terminal error with
the observed inventory. The verifier inventories before recomputation and again
before publication and requires the candidate's explicit empty inventory. A missing first-stage attempt cannot become a verified
experiment; verified optimizer failure explicitly reports zero recomputed items.

At this implementation checkpoint, the repair had run only mocked CPU tests; no
real calculator or downstream scientific output had yet been produced.

- 2026-09-06 19:45 PDT (hermes-custom-build-001; profile=workstation) — The bounded experiment and a distinct delegated verifier process completed. The sole owner-conditioning optimizer spent its one 100-step HF/STO-3G budget without convergence, changed three original proton owners (`H41:O11->O15`, `H48:O26->O28`, `H54:O28->O11`), and failed independently recomputed projected-gradient gates (`1.3994921840965815e-3/4.996643654196831e-3 Eh/Bohr` RMS/max versus `3.0e-4/4.5e-4`). `H52` itself remained unambiguously owned by `O15` (0.960017 A; 0.686038 A margin) and the constraint residual, finite-energy, frozen-shell, and collision gates passed, but those narrower passes cannot authorize production. Constrained B3LYP, release, and PHVA therefore did not run. The independent verifier rehashed source and receipts, recomputed four calculator quantities, found no forbidden artifacts, and published verified terminal SHA-256 `6786d6796a4b233f6a73b1272cf6804783b8f6cb3634b4ec6cfe9bd96c57171c`.

## Result

**Independently verified outcome: inconclusive terminal failure.** The exact
`oaa-neutral-n2-s2` H52-O15-constrained route is exhausted with one
owner-conditioning call and zero retries. Its endpoint was persisted with finite
energy (`-3316.690344523382 Eh`) and passed the H52 residual, frozen-shell, and
collision gates, but the optimizer did not converge, three other proton owners
changed, and independently recomputed projected gradients were above both fixed
thresholds. No constrained-production or released-production optimizer, PHVA,
barrier, store, Petra fragment, family publication, or CALCULATIONS value ran or
was emitted.

The verified terminal receipt is
`/mnt/data/vsletten/dissertation-data/a3g-oaa-neutral-n2-proton-microstate-stability/verified-terminal.json`
(SHA-256 `6786d6796a4b233f6a73b1272cf6804783b8f6cb3634b4ec6cfe9bd96c57171c`);
the hash-bound stage receipt is `7ff5632ef99476284a8c1254a74a8a11b7f0144c36053aca7faa7639bb50866f`.
Per this card's zero-retry terminal contract, Oaa-neutral n=2/4/6 is closed and
parent A3 returns only to a different independently gated family. No A3h retry
card is authorized or filed.
