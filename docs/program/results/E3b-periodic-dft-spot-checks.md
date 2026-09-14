# E3b periodic-DFT spot checks — active survey calibration

**State:** active; frozen-path survey profiles in progress; no DFT barrier or correction emitted
**Run date:** 2026-09-14  
**Operator:** `(hermes-custom-build-001; profile=workstation)`

## Verdict

E3b is **not complete**: **0/3 frozen-path survey profiles are complete** and
no DFT barrier or correction is emitted. The immutable 2×2×1 preparation and
evidence contract remain valid. POLICY §12 now governs this platform-test phase:
three full endpoint optimizations plus converged eight-image CI-NEBs would
exceed the card's roughly one-day compute envelope, so the active route is three
PBE-D3 single-point profiles on the immutable classical images, paired with
fresh matched-cell classical profiles.

The only supportable quantity from that route is the energy rise along each
tested frozen path. It is not a DFT-relaxed activation barrier, a converged DFT
minimum-energy path, or a production result. The result stays incomplete until
all eight images and the corresponding matched-cell classical path pass.

## Current execution evidence

The first dehydroxylate endpoint GEO_OPT ran in a bounded transient unit for
1,830.018 s. Ten SCF/geometry cycles converged, but the geometry did not; the
typed receipt is `incomplete-timeout`, return code `-9`, with no normal
termination or barrier. Cleanup proved the process group and container absent,
and the immutable prepared manifest remained valid. This timing is the measured
basis for narrowing the card rather than replaying a multi-day route.

The earlier 3×2 timing-based “outside the workstation envelope” conclusion was
rejected by the card's one allowed cold review. A wrapped, neutral 2×2 route
cell exists, and an unfinished single-image timing attempt cannot establish a
BAND iteration floor. This result therefore makes **no claim that E3b cannot fit
the workstation envelope** and does not request a compute-rental decision.

The E3a classical values remain unchanged and explicitly incomplete:
68.414811 kcal/mol for reconstructed divacancy Ar, 64.095991 kcal/mol for the
local-dehydroxylate Ar hop, and 90.744700 kcal/mol for the Xe screen. They are
planning references, not converged calibration anchors.

## Corrected 2×2 preparation

The authoritative crop is the periodic wrapped window beginning at source tile
`(5,2)`, containing tiles `(5,2)`, `(0,2)`, `(5,0)`, and `(0,0)`. It retains
Nteme route 1 (moving site 3; vacancies 4 and 88), exact four-unit-cell
stoichiometry, net charge `1.78e-15 e`, and minimum pair distance `1.007 Å`.

| Check | Atoms | E3a reference status | Preparation gate |
|---|---:|---|---|
| reconstructed divacancy Ar | 334 | incomplete-convergence | pass |
| local-dehydroxylate Ar | 331 | incomplete-convergence | pass |
| reconstructed divacancy Xe | 334 | incomplete-convergence | pass |

The cell is 10.3002083 × 17.8316813 × 19.549356 Å with the source triclinic
angles. The dehydroxylate model preserves five-coordinate Al source IDs 10 and
11. All three models pass composition, neutrality, route, atom-identity,
minimum-distance, electron-parity, and topology preparation gates.

The generated method is CP2K 2024.3, periodic GPW PBE-D3 zero damping,
DZVP-MOLOPT-SR-GTH/GTH-PBE, Γ point, 600 Ry cutoff, 60 Ry relative cutoff, and
eight CI-NEB replicas. The container is pinned to
`cp2k/cp2k@sha256:979e011c2ea15a2a56eae16b5f3a2d214142bb74fe644bfb302d366dddc2a443`.

## Live CP2K evidence

A 4-MPI × 4-OpenMP `ENERGY_FORCE` attempt on the 334-atom reconstructed-Ar
initial state remained unfinished when the execution channel reached its
420-second hard ceiling. No `ENERGY| Total FORCE_EVAL` line or normal program
termination was emitted. The still-running container was detected and stopped;
no CP2K container or process remained afterward.

This establishes only that the full-method single-image probe takes more than
the available inline channel window on this machine. It is **not** a completed
timing, a BAND lower bound, or a campaign-feasibility verdict.

## Fail-closed corrections

The cold review found two paths that could have emitted an unsupported number.
Both are closed in the branch:

1. the analyzer now refuses manually entered observation booleans/barriers and
   returns `incomplete-unverified-observation` until raw CP2K and LAMMPS parsers
   with artifact hashes exist;
2. an E3a source reference not typed `converged` returns
   `incomplete-source-classical`, so the three current incomplete E3a values
   cannot be promoted into a correction or endorsement.

`E3b1-periodic-dft-evidence-contract` owns the remaining evidence work: separate
immutable inputs from runtime outputs, parse/hash raw solver evidence, resolve
the incomplete classical anchors, define the five-coordinate interaction
criterion, and run a bounded completed 2×2 timing probe.

## Artifacts

Authoritative 2×2 evidence root:

`/mnt/data/vsletten/dissertation-data/e3b-periodic-dft-spot-checks-2x2-20260914/`

- root evidence manifest: 71 files / 57,403,263 bytes
- manifest SHA-256: `1c2e837b8c3dc4afeda8d31a564695a9eb8ec1675c958029ea77ea5ece7eb529`
- `prepared/preparation.json`: `f7b776a19b138b68228809fde75201042abc3208ff8956cb4f06667dfeb97093`
- initial prepared-manifest SHA-256:
  `46fe4d03c50693e6a681f9a6bda500455d7695289412c6f3e84debf75755c5fe`
- incomplete 2×2 smoke log SHA-256:
  `fa4d4a541d1e3409dc085b4db805ed73f5c8db22414c2117097a8e15bab2f5a8`

The initial prepared manifest is intentionally reported as **not reusable after
execution** because CP2K wrote runtime restart data inside that tree. E3b1 must
separate immutable inputs from outputs before further campaign work.

The superseded 3×2 exploratory evidence remains under
`/mnt/data/vsletten/dissertation-data/e3b-periodic-dft-spot-checks-20260914/` for
audit only; it does not support this verdict.

## Verification

- `uv run --with openpyxl --with pytest ...`: **33 passed** after the cold-review
  corrections.
- Ruff check/format, Python compilation, and `git diff --check`: pass.
- The 2×2 real-input preparation generated 334/331/334-atom models with all
  preparation gates true.
- Post-probe cleanup: zero CP2K containers and zero CP2K/`mpirun` processes.
