# A3e-oss-neutral-n1-proton-microstate-stability — test the H35 owner basin

- status: done
- track: A (geochemistry)
- priority: P1
- machine: workstation
- depends: A3 OSS-neutral n=1 advisory failure recorded and supervisor merged
- claimed-by: hermes-custom-build-001

## Objective

Determine whether the requested `oss-neutral-n1-s2` reactant proton assignment has
a usable B3LYP/def2-SVP/DF stationary basin. The parent A3 campaign's bounded
supervisor stopped correctly after a converged HF/STO-3G advisory preoptimization
moved termination proton `H35` from builder owner `O21` to `O14`, but the rejected
endpoint was discarded before persistence and no production-level calculation ran.

Execute exactly one materially different, finite owner-constrained → released
production-PES experiment for this single cell. Produce a hash-bound typed verdict:
accepted reactant minimum, production-level no-basin rejection, or terminally
inconclusive. This card is not authority to replay the unconstrained advisory route.

## Evidence contract

Rehash before use and bind every derived receipt to exact branch/source, atom order,
charge/spin, frozen shell, method settings, seed geometry, and evidence-root identity:

- root: `/mnt/data/vsletten/dissertation-data/task274-a3-oss-neutral-family-20260906/`
- family receipt SHA-256: `c392dafe3c60f75687954295ae5b54da4c97a24ae1e54c9ec5b7cb8fdf9d9e53`
- child log SHA-256: `b04a2991be5621925756222ae66bd41490c3152dc1ee442fcdfece590d346173`
- exact seed `complex_guess.xyz` SHA-256: `e67374b19a681ea6420d3b048a60907ef0848cd2d8f184581c67ca74ee4e570a`
- metadata SHA-256: `09416f72a809f8208a69163a4f0c9dd6be16883186d76b5012c1cbfcae384354`
- launch/restoration SHA-256: `59e43dc58ccb2835319713f12fb171d14d10449a19acb02294d9eabd279d1416` / `4dff92f549962c2e6c801f159d1f351a12cc00bc55e6e002057b660b1bae5b34`
- execution source before closeout merge: `ec5bdb059911e32e4b790f849ff97595587c8519`

The old rejected advisory endpoint does not exist. Do not invent its geometry, hash,
energy, or production meaning from the log.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one branch, branch-matched
  worktree, and one PR carrying implementation, evidence summary, card, PLAN, and
  STATUS updates.
- Exactly one `oss-neutral-n1-s2` experiment. No n=2–4, saddle, barrier,
  thermochemical barrier, provenance store, Petra emission, or family publication.
- No identical replay of the unconstrained HF/STO-3G advisory preoptimization.
- Use the canonical GPU lease; threads <=16, nice >=10, durable tee logs, finite
  cgroup/runtime ceiling, atomic terminal receipt, and restoration dead-man.
- Persist every raw endpoint and its finite energy/gradient/settings receipt before
  applying a structural or ownership gate. Rejected evidence must survive.
- Never promote optimizer-local convergence alone. Independently recompute projected
  gradients, constraint residuals, ownership, frozen shell, collisions, and PHVA.

## Acceptance

1. Rehash the full Evidence contract before launch; mismatch stops without a
   calculator call.
2. Starting from the exact seed, constrain `H35-O21` at its seed distance and run
   one owner-preserving conditioning stage, then one fresh constrained
   B3LYP/def2-SVP/DF optimization bounded at 100 steps.
3. Release is allowed only when the constrained production endpoint has orthogonal
   projected RMS/max gradients <=`3.0e-4`/`4.5e-4 Eh/Bohr`, `H35-O21` residual
   <=`1.0e-4 A`, exact original proton ownership, exact frozen shell, no collision,
   and finite energy/gradient receipts. Release uses one fresh B3LYP optimizer,
   bounded at 100 steps, with no inherited incompatible Hessian.
4. Record exactly one typed terminal outcome:
   - **accepted reactant minimum:** the unconstrained endpoint retains every original
     proton owner and passes fresh projected-gradient, frozen-shell, collision,
     finite-energy, and PHVA zero-significant-imaginary-mode gates;
   - **production no-basin:** the released endpoint converges lower in finite
     B3LYP energy to unambiguous `H35->O14` ownership, with independently reproduced
     endpoint forces/energy and all non-owner integrity gates intact;
   - **inconclusive terminal failure:** any other result, including exhaustion or an
     unavailable release. It authorizes no replay or surrogate output.
5. Raw/projected endpoints, settings, stage receipts, terminal verdict, and complete
   hashes are durable outside Git; source contains only payload-free provenance and
   the scientific conclusion.
6. Full relevant CPU-only tests, whole-QM Ruff check/format, and `git diff --check`
   pass. If an accepted minimum is produced, parent A3 may resume only from that
   exact promoted artifact; otherwise update A3's terminal disposition explicitly.

### Verify

A different worker performs a cold, optimizer-free verification from the original
H35 owner complaint and raw receipts. It independently recomputes ownership,
constraint residual, projected-gradient, finite-energy, frozen-shell/collision, and
PHVA gates and confirms the typed terminal outcome. Before that verification lands,
report only `unverified: candidate outcome`.

## Progress

- 2026-09-06 09:02 PDT (hermes-custom-build-001; profile=workstation) — DONE with an independently verified `inconclusive terminal failure`. Cold verifier `hermes-task292-cold-verifier` rehashed the complete source/receipt chain and recomputed B3LYP energy plus projected gradients with zero optimizer calls. The constrained production endpoint retained every original proton owner and passed RMS gradient, H35-O21 residual, frozen-shell, finite-energy, and collision gates, but its independently recomputed max projected gradient `4.8363991055655475e-4 Eh/Bohr` reproducibly exceeds the fixed `4.5e-4` gate; release therefore correctly remained unspent and no PHVA or no-basin claim exists. READY A3f owns one CPU-only terminal-path adjudication; A3e may not replay.
- 2026-09-06 08:10 PDT (hermes-custom-build-001; profile=workstation) — Atomically claimed `agents/A3e-oss-neutral-n1-proton-microstate-stability` from exact `origin/main@d460f1b78defb5477452fff7cbe5ae847dbe462d`. Implemented the single authorized three-stage route: H35-O21-constrained HF/STO-3G conditioning, a fresh constrained B3LYP/def2-SVP/DF optimizer, and a fresh released B3LYP optimizer only after independent projected-gradient/residual/owner/shell/collision gates pass. Every spent stage reserves its one 100-step budget before launch, persists a raw endpoint and pending receipt before structural gates, binds parent/source/settings/code hashes, and emits only an unverified typed candidate. A separate optimizer-free verifier requires a genuinely different identity and recomputes production energy, projected gradients, ownership, residuals, integrity, and accepted-minimum PHVA. Ten focused adversarial tests plus focused Ruff/format and `git diff --check` pass; no calculator call has run yet.

- 2026-09-06 07:06 PDT (hermes-custom-build-001; profile=workstation) — Filed from the fail-closed parent family attempt after two cold read-only reviews agreed that HF-level owner transfer is strong evidence against the builder assignment but insufficient for a production-PES no-basin verdict because the rejected endpoint was not persisted. This card is the only authorized next OSS-neutral n=1 calculation; it spends one materially different finite constrained-to-released production route and forbids all replay and downstream family work.

## Result

- 2026-09-06 09:02 PDT (hermes-custom-build-001; profile=workstation) — **Verified terminal outcome: `inconclusive terminal failure`.** Exact execution source `333917de87c32b6be4ecc886662dd6666947740c`; executor candidate SHA-256 `b9fa9c4b0940b4e5ff7c40b83ffb106e8b3a4f4f2cbac355030c125eda9ab378`; independent terminal SHA-256 `756e1e885150556ac9b2da4a63ab4e4965246cd1072f83b09a9b5ca6629c86ff`; recomputed-evidence SHA-256 `82a5e9d0530568fd5040ffa2ca43fdcd4c6ff799fb94bed110f3faef291316fa`. Durable evidence root: `/mnt/data/vsletten/dissertation-data/a3e-oss-neutral-n1-proton-microstate-stability/`.
- 2026-09-06 09:02 PDT (hermes-custom-build-001; profile=workstation) — Conditioning and constrained production each spent exactly one 100-step budget with zero retries. The constrained endpoint is owner-retaining at `H35:O21`, residual `4.201582921470326e-6 A`, RMS/max projected gradients `1.0358722090842382e-4/4.8363991055655475e-4 Eh/Bohr`; the fixed maximum-gradient gate is red, so release is `not-run`. No accepted minimum, production no-basin verdict, PHVA, n=2-4 run, saddle, barrier, store, Petra fragment, or CALCULATIONS value was emitted.
- 2026-09-06 09:02 PDT (hermes-custom-build-001; profile=workstation) — The finite compute unit is terminal, restoration succeeded, the canonical GPU lease is absent, both email pipelines and queue/board timers are active, and API health is HTTP 200 `{"status":"ok"}`. Fresh full fast QA is `626 passed, 4 deselected`; whole-QM Ruff check/format and `git diff --check` passed. This terminal outcome closes A3e without loosening thresholds or replaying a stage; READY `A3f-oss-neutral-n1-terminal-adjudication` owns the exact evidence-only parent disposition.
