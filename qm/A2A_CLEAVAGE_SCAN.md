# A2a coupled-coordinate cleavage scan contract

> Change this contract first, then make `scripts/a2a_cleavage_scan.py` and its
> tests match. A coarse grid point is never a verified transition state.

## Fixed physical mapping

| Role | Atom index | Contract |
|---|---:|---|
| bridging oxygen `Obr` | 0 | shared atom in both constrained distances |
| rupturing silicon `Si` | 1 | `xi1 = d(Si1, Obr0)` |
| attacking water oxygen `Ow` | 15 | H16 must remain uniquely owned here at exact I and P |
| coupled water proton `Hw` | 16 | `xi2 = d(H16, Obr0)`; contracts by 1.127 Å from I to P |
| already bridge-owned proton | 17 | rejected as a degenerate scan coordinate |

The exact archived I/P endpoints are ordered identically. `Si1-Obr0` expands
from 1.953932 to 3.732998 Å. `H16-Obr0` contracts from 2.979697 to
1.852421 Å, while `H17-Obr0` changes by only 0.013367 Å. The scan must fail
closed if this mapping or physical-H ownership drifts.

## Grid and persistence invariants

- Default grid is 9x9. Exact I is cell `(0, 0)` and exact P is `(7, 7)`.
  Row/column 8 extends one grid step beyond P on each axis.
- Traversal is a deterministic adjacent-cell serpentine path. Every cell records
  its predecessor and exact seed geometry fingerprint.
- Every cell binds scan version, endpoint fingerprints, settings, targets,
  predecessor, seed, convergence contract, optimized geometry, energy receipt,
  target residuals, typed topology, and content hashes.
- Existing checkpoints are accepted only on exact identity equality. Partial or
  drifted checkpoints fail; they are never silently reseeded.
- Both constraints use simultaneous Sella internal coordinates. Optimizer
  convergence and both distance residuals are separate acceptance gates.
- The extension cells diagnose behavior beyond P but may not be used as a
  low-energy detour in the physical I-to-P minimax path.

## Decision matrix

Classification runs only after all cells have finite energies. The physical
minimax path is restricted to rows/columns `0..product_index`.

| Typed proton-first local minimum? | Physical-path bottleneck location | Height above I | Outcome | Required follow-through |
|---|---|---:|---|---|
| yes: H16 is Obr-owned **and** Si-Obr remains bonded | any | any | `proton-first-minimum` | Minimum-gate I-prime; split and gate both elementary steps |
| no | strictly interior to physical I-P rectangle | `> 2.0 kJ/mol` | `interior-crest-above-threshold` | Seed exactly one full-dimensional Sella; unchanged two-step Hessian, full IRC, ownership, topology, and basin gates |
| no | strictly interior | `<= 2.0 kJ/mol` | `barrierless-shelf` | Fresh downhill release along the valley into exact typed P; accepted R-I saddle is the route's production TS |
| no | physical-grid boundary, including P | any | `barrierless-shelf` | Never promote a boundary maximum to a saddle; require the same downhill typed-P proof |

### Precedence

1. **Proton-first minimum first:** it changes the elementary-step graph and must
   not be hidden by a higher path crest.
2. **Interior crest second:** only a physical-domain interior bottleneck above
   the predeclared threshold can justify one new localization.
3. **Barrierless shelf otherwise:** threshold equality is barrierless, boundary
   maxima are not saddles, and extension points cannot rescue or condemn the
   physical route.

## Barrierless downhill release proof

A `barrierless-shelf` classification is incomplete until one fresh
full-dimensional minimum optimization releases **both** scan distances. The seed
is deterministic: the first physical minimax-path cell whose recomputed exact
endpoint identity equals the archived hydrolyzed product (basin, every physical
hydrogen owner, and every heavy-atom bond). Index position or the coarse basin
triple alone is insufficient.

The release is checkpointed under
`cleavage/barrierless-downhill-release-v1/` and binds the classification,
settings, seed, exact product reference, optimizer, force gate, and step bound.
Acceptance requires all of the following:

- no `Si1-Obr0` or `H16-Obr0` distance constraint remains;
- a nonzero Cartesian displacement greater than `1e-3 A` occurs;
- the final electronic energy is no higher than the constrained seed beyond a
  `1e-8 Eh` numerical tolerance;
- exact product typed identity matches the archived P endpoint; and
- geometry, energy, identity, and content hashes are durably receipted.

An uphill, zero-step, identity-drifted, partial, or checkpoint-drifted release
fails closed. No cleavage Hessian or IRC is commissioned after this proof: the
accepted R-I saddle is the route production transition state.

## Named proof

`tests/test_a2a_cleavage_scan.py` covers:

- exact endpoint inclusion, one-step extension, and adjacent traversal;
- fixed H16 mapping and fail-closed H17/endpoint drift;
- checkpoint resume without recomputation and predecessor-seed chaining;
- topology-backed proton-first precedence with no index-only inference;
- interior-crest, threshold/barrierless, boundary-maximum, and extension-detour
  cases;
- deterministic first exact-product valley seed selection plus fresh
  unconstrained downhill, nonzero-step, exact-identity, and uphill-rejection
  gates; and
- GPU-memory CLI wiring and the fixed external evidence namespace.

The aggregate QM suite must also remain green because the driver reuses shared
checkpoint, constrained-relaxation, endpoint-identity, and etiquette machinery.
