# E3b periodic-DFT spot checks — completed frozen-path survey

**State:** done; 3/3 PBE-D3 frozen-path profiles and 3/3 matched-cell classical profiles complete
**Run date:** 2026-09-14
**Operator:** `(hermes-custom-build-001; profile=workstation)`

## Verdict

E3b is complete under the POLICY §12 platform-test amendment. All three
immutable eight-image paths have hash-bound PBE-D3 single-point profiles and
fresh LAMMPS profiles on the exact same coordinates. The survey comparison is:

| Route | PBE-D3 frozen-path rise (kcal/mol) | Classical same-path rise (kcal/mol) | DFT − classical (kcal/mol) | Verdict |
|---|---:|---:|---:|---|
| dehydroxylate-lattice Ar | 297.913659 | 9877.633000 | -9579.719341 | `disagrees` |
| reconstructed-replication Ar | 305.946204 | 9877.629000 | -9571.682796 | `disagrees` |
| xenon-divacancy | 425.854862 | 39166.519000 | -38740.664138 | `disagrees` |

`close-enough` means an absolute difference no greater than the repository's
5.0 kcal/mol transfer tolerance. Every route exceeds it by thousands of
kcal/mol, so all three verdicts are `disagrees`. The classical potential is not
endorsed by this coordinate-matched survey, and these numbers are not applied
as corrections or deck replacements.

These are energy rises along immutable classical paths. They are **not**
DFT-relaxed activation barriers, converged DFT minimum-energy paths, CI-NEB
results, or production energetics. The unusually large rises are reported
rather than laundered into a more attractive quantity: the coordinate-matched
classical evaluator and PBE-D3 both see highly unfavorable interior images, but
their scales disagree catastrophically.

## Acceptance evidence

### PBE-D3 profiles

Each route completed 8/8 CP2K single points with converged SCF, normal
termination, coordinate/input/native-output hashes, bounded controller
identity, and confirmed container cleanup:

| Route | Elapsed (s) | Frozen-path rise (kcal/mol) | Receipt SHA-256 |
|---|---:|---:|---|
| dehydroxylate-lattice | 2412.217 | 297.913659 | `87a6775e8d78efc663ffdaa5abce71ba9906f49bd07d8fe8934e23e6506c85a8` |
| reconstructed-replication | 2248.029 | 305.946204 | `dcc8a92c2b66f212baefd962a1431d1e0b1b3f8c01abc810a9b86b5c1d55bbb3` |
| xenon-divacancy | 2618.104 | 425.854862 | `aa5d6860b19855ba91c9bcea8d32daf065f35bf72d05cf9cf6a9bb58bd311459` |

The dehydroxylate preparation's typed criterion identifies source Al atoms 10
and 11 as transformed from six- to five-coordinate by loss of source oxygen 76
with no gained neighbor. The moving Ar enters the declared 6.0 Å interaction
shell in images 0–5 for both sites. Thus the required five-coordinate-Al
interaction gate passes; this is a local post-dehydroxylation Ar hop, not the
dehydroxylation reaction itself.

### Matched-cell classical profiles

`scripts/e3b_classical_frozen_path.py` binds each ordered XYZ row to the
explicit source atom ID in `atom-map.json`, rewrites only those coordinates in
the prepared LAMMPS topology, and runs a bounded zero-step energy evaluation.
All three routes completed 8/8 images in about five seconds each. The runner
emits no profile number after any timeout, nonzero return, malformed identity,
non-finite energy, missing normal termination, or prepared-manifest drift.

| Route | Runtime-manifest SHA-256 | Receipt SHA-256 |
|---|---|---|
| dehydroxylate-lattice | `4ac1048830471a14836d953c9220f91c57c520c6080356f396b4024d300b97a4` | `700004e307ac47056c3c28d4d0c4b556bce25f48f3086fc0ab0d7ba6351296a5` |
| reconstructed-replication | `ed325856e21bd28236556e9a064376fb1e1af9c6996937800ca0970e892e985d` | `e25840dc27bc289098a89c526721a0c5a4c96d382a1a869aad010f4f1dfc0ba6` |
| xenon-divacancy | `a89dacff84fae60e6b79d9b737638288833652d4f29873a2649abf15d23b0d14` | `729189c9bba8e3298b6f16473bfeeaab8e7ce51beba6c84bbc251e78e7ded9b7` |

## Scope deviation

The original card requested three fully relaxed endpoint calculations and
periodic CI-NEBs. The first 331-atom dehydroxylate endpoint GEO_OPT ran for
1830.018 s, completed ten converged SCF/geometry cycles, and timed out without
endpoint convergence or normal termination. Its typed
`incomplete-timeout` receipt emitted no numeric barrier.

POLICY §12 supersedes that multi-day route on this one-workstation platform
test: survey-tier methods are the banked value, one QM unit is capped at four
hours, and a card requiring more than roughly one day of workstation compute
must be narrowed. E3b therefore banks the six complete coordinate-matched
profiles above and makes no full-CI-NEB feasibility claim.

The earlier 3×2 timing-based “outside the workstation envelope” conclusion was
also rejected by the card's one prior cold review. A wrapped, neutral 2×2 route
cell exists, and an unfinished force evaluation cannot establish a BAND
iteration floor. No compute-rental decision follows from this card.

## Preparation and method

The authoritative crop is the periodic wrapped 2×2×1 window beginning at
source tile `(5,2)`, containing tiles `(5,2)`, `(0,2)`, `(5,0)`, and `(0,0)`.
It retains Nteme route 1 (moving site 3; vacancies 4 and 88), exact four-unit-
cell stoichiometry, net charge `1.78e-15 e`, and minimum pair distance
`1.007 Å`. The three prepared systems contain 331, 334, and 334 atoms.

The PBE-D3 method is CP2K 2024.3 periodic GPW, DZVP-MOLOPT-SR-GTH/GTH-PBE,
Γ point, 600 Ry cutoff, 60 Ry relative cutoff, D3 zero damping, with image
`cp2k/cp2k@sha256:979e011c2ea15a2a56eae16b5f3a2d214142bb74fe644bfb302d366dddc2a443`.
The classical method is the prepared E3a matched-cell LAMMPS force field in
`real` units, evaluated without minimization on the same eight coordinates.

Prepared evidence root:

`/mnt/data/vsletten/dissertation-data/e3b1-periodic-dft-evidence-contract-20260914/prepared-reviewed/`

Prepared manifest SHA-256:
`d796a30ce3e49b0150b4dcd42e625177ae7f44c4ab0ce5969e3e6d770c4a7ab0`.

## Durable receipts

Calibration root:

`/mnt/data/vsletten/dissertation-data/e3b-periodic-dft-calibration-20260914/`

The independent replay of all 48 image records, both runtime manifests per
route, all receipt arithmetic, every immutable coordinate hash, all CP2K
termination/SCF/cleanup gates, and all LAMMPS normal-termination/energy gates is:

`receipts/comparison-verification-02.json`

SHA-256:
`b0f7b58bd3bf7a79e6d86365f04cea16a57eaeb07f5efc0f0b92e2351e5a8cc9`.

The single POLICY §13 cold review found and closed three delivery-changing
issues: arbitrary five-number diagnostics could supersede the declared LAMMPS
thermo row; the receipt copied rather than recomputed atom-map identity; and no
deck decision fragment recorded which value E4 used. The corrected runner was
replayed for all 24 classical images with unchanged results. Non-blocking
runner/revision attestation is deliberately deferred to
`E3b2-classical-receipt-code-binding`.

## Verification

- E3-focused PBE-D3/evidence-contract/classical-runner suite: **60 passed**.
- Ruff check and format, Python compilation, and `git diff --check`: pass.
- Receipt replay: **48/48 images** hash-bound and complete; profile arithmetic
  independently recomputed.
- Cleanup: every recorded CP2K container ID is absent; no CP2K container or
  solver process survives.
- No numeric correction, classical endorsement, relaxed DFT barrier, CI-NEB,
  or production-energy claim is emitted.
