# E3a-classical-neb-barriers — the uncomputed middle, route 1

- status: blocked
- track: E (muscovite / the perfect circle)
- priority: P1
- machine: any (classical-potential NEB is CPU work; large-cell sweeps or DFT escalation go to the workstation)
- depends: TASK-276-e3a-atomistic-input-contract (mission-control)
- claimed-by: hermes-custom-build-001

## Objective

Execute phase 3 route 1 of `docs/scoping/ar-muscovite.md` §6: replicate
Nteme et al. (2022)'s classical-potential NEB setup in 2M₁ muscovite, then
compute the barriers the E2 deck currently carries as labelled proxies. This
is a periodic-slab problem, not a cluster problem — quarry's kaolinite
cluster workflow is the wrong shape; use classical potentials (Nteme's
published setup as the reference) with NEB+TST.

Validation gate first, then the campaign:

1. **Replication gate**: reproduce the published interlayer divacancy
   migration barrier (~66 kcal/mol, Nteme 2022) in the pristine gallery
   before any new number is trusted. This is the go/no-go for the setup.
2. **The missing barriers** (~6–10 types × 2 species, provenance-stamped):
   - Ar hop in the **dehydroxylate lattice** (collapsed interlayer) — the
     key missing number; 1998 had only the Robbins 1972 proxy.
   - Ar hop in **extended-defect zones** — test whether Nteme 2023's 10⁶×
     diffusivity enhancement is barrier-equivalent (~26 kcal/mol lower).
   - **Octahedral trap escape** for recoil-implanted ³⁹Ar (E2 proxy: 64).
   - **Delamination-adjacent** gallery hops and surface release steps
     (E2 proxy: 58 ± 5).
   - **Xe analogues** for the discriminant set (Villa arbitration).
   - Environment-resolved **dehydroxylation table** refinement vs the
     57/60/69 published envelope.

## Constraints

- Follow `docs/program/PROTOCOL.md` and fleet `POLICY.md`; one
  branch/worktree/PR; QI2 GPU-lease etiquette applies to any GPU stage.
- Every barrier lands with method provenance (potential, cell, images,
  convergence) per house rules; no proxy value may be silently replaced —
  each replacement is an explicit deck diff with the computed source cited.
- Computational failure is `incomplete-*`, never a rejection; barriers that
  resist convergence stay proxies with the failure documented.
- Archive heavy artifacts under `/mnt/data/vsletten/dissertation-data/`
  with independent hashes.

## Acceptance

- The replication gate reproduces Nteme's divacancy barrier within a stated,
  justified tolerance, documented before campaign numbers are reported.
- ≥ the dehydroxylate-lattice, extended-zone, and trap-escape barriers are
  computed for Ar (Xe where feasible), provenance-stamped, and delivered as
  a deck-ready fragment replacing the corresponding E2 proxies.
- A results doc under `docs/program/results/` tabulates computed vs proxy
  vs published values with honest uncertainty statements.
- Tests/lint green per repo gates; card/PLAN/STATUS bookkeeping in the PR.

## Progress

- 2026-09-02 22:17 PDT (hermes-custom-build-001; profile=workstation) —
  Replication gate executed against the recovered Nteme 2022 supplement. An
  eight-image ClayFF/LAMMPS CI-NEB run reproduced exact divacancy route 1 at
  68.956 vs 67.644 kcal/mol (+1.312; +1.94%), inside the predeclared
  ±5 kcal/mol gate. The bounded 13,000-step run reached max atom force
  5.08e-5 but whole-replica norm 0.0411 vs requested 0.01, so tighter
  convergence remains explicitly `incomplete-convergence`. Campaign numbers
  were not fabricated: the sources contain no atomistic dehydroxylate,
  extended-zone, octahedral-trap, or Xe endpoints and omit route-specific
  charge-compensation sites. BLOCKED on mission-control TASK-276 to deliver
  that source-backed atomic input contract. Harness, exact 1,080-row source
  catalogue, hashes, commands, and run receipt are in
  `docs/program/results/E3a-classical-neb-barriers.md`.

- 2026-08-29 — Fable — re-scoped machine: any per Victor's fleet-routing
  directive (POLICY v13 §10): the route-1 replication gate and barrier matrix
  run on classical potentials (CPU, laptop-feasible at unit-cell scale). Record
  any host ceiling per the E2b pattern; DFT-tier work stays in E3b.

- 2026-08-27 — filed by Fable. Phase 2 (E2) is done on proxy barriers; this
  card is the route-1 campaign that upgrades them to computed kinetics.
  Route-2 CP2K spot checks are carded separately (E3b).
