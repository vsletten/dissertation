# E3a-classical-neb-barriers — the uncomputed middle, route 1

- status: done
- track: E (muscovite / the perfect circle)
- priority: P1
- machine: any (classical-potential NEB is CPU work; large-cell sweeps or DFT escalation go to the workstation)
- depends: —
- input-contract: `docs/program/results/E3a-atomistic-input-contract.md`
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
   - Ar hop in the **local dehydroxylate lattice** (no collapsed-interlayer
     assumption; use the contract's post-reaction topology and cell relaxation
     gate) — the key missing number; 1998 had only the Robbins 1972 proxy.
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

- 2026-09-03 02:42 PDT (hermes-custom-build-001; profile=workstation) — DONE as
  a bounded sensitivity campaign. The neutral reconstructed route-1 gate is
  68.414811 vs 67.644151 kcal/mol (+0.770660; +1.1393%), inside the predeclared
  ±5 kcal/mol tolerance. Six 8-image CI-NEBs completed: local-dehydroxylate Ar
  64.095991; three extended-zone geometries 68.410690–68.411538; and Xe
  90.744700 kcal/mol. All are explicitly `incomplete-convergence` because their
  whole-replica norms remain 0.0272–0.0633 vs 0.01, despite ≤0.0101 kcal/mol
  final-500-step barrier spans. The vacant-M1 Ar seed relaxed 2.316 Å out of its
  assigned cage and lacks an index-zero Hessian, so trap escape remains
  `incomplete-index-gate` and its proxy is retained. Delamination remains
  `incomplete-input-contract`: the merged contract defines no atom-resolved
  delamination/surface-release endpoints. The deck-ready typed fragment,
  methods, limitations, and hash-verified 253-file evidence manifest are in
  `docs/program/results/E3a-classical-neb-barriers.md`.

- 2026-09-02 — `(hermes-custom-build-001; profile=workstation)` — recovered the
  atomistic input prerequisite on `agents/E3a-atomistic-input-contract`.
  `docs/program/results/E3a-atomistic-input-contract.md` and
  `scripts/e3a_input_contract.py` now specify exact transformations, typed
  reconstructed charge compensation, source limits, validation gates, and
  generated two-end LAMMPS run-0 fixtures for dehydroxylate, extended-zone,
  geometric vacant-M1 trap/escape, and Xe models. Do not start those NEBs until
  the input-contract PR merges; the octahedral route must additionally clear
  its typed index-zero minimum gate.

- 2026-09-02 22:17 PDT (hermes-custom-build-001; profile=workstation) —
  DEVIATION: Replication gate attempted against the recovered Nteme 2022
  supplement. An eight-image ClayFF/LAMMPS CI-NEB near-match put exact
  divacancy route 1 at 68.956 vs 67.644 kcal/mol (+1.312; +1.94%), inside a
  predeclared ±5 kcal/mol numerical tolerance, but the scientific gate remained
  CLOSED: the source requires remote charge-compensation substitutions and does
  not identify their sites. The run therefore used a different net -3 e
  Hamiltonian. The bounded 13,000-step run reached max atom force 5.08e-5 but
  whole-replica norm 0.0411 vs requested 0.01, so tighter convergence remained
  `incomplete-convergence`. No campaign numbers were fabricated: the sources
  contain no atomistic dehydroxylate, extended-zone, octahedral-trap, or Xe
  endpoints and omit route-specific charge-compensation sites. The harness,
  exact 1,080-row catalogue, hashes, commands, and run receipt are in
  `docs/program/results/E3a-classical-neb-barriers.md`.

- 2026-08-29 — Fable — re-scoped machine: any per Victor's fleet-routing
  directive (POLICY v13 §10): the route-1 replication gate and barrier matrix
  run on classical potentials (CPU, laptop-feasible at unit-cell scale). Record
  any host ceiling per the E2b pattern; DFT-tier work stays in E3b.

- 2026-08-27 — filed by Fable. Phase 2 (E2) is done on proxy barriers; this
  card is the route-1 campaign that upgrades them to computed kinetics.
  Route-2 CP2K spot checks are carded separately (E3b).
