# Architecture decision record — graph-viz rebuild

Decided 2026-07-12 (the diagnosis + rebuild night). All numbers measured that
evening on Victor's workstation: RTX 4090, Chrome/WebGL2 via ANGLE, 60 Hz
display (frame times at ~16.7 ms are vsync-capped, not GPU-limited).

## The diagnosis

**Legacy viewer** (`dissertation/main/viewer`: React 19 + react-three-fiber,
one React component + one Mesh per atom and per bond, drei `Text` labels on by
default):

| workload | fps | mount time |
|---|---|---|
| kaolinite golden `start.msi` (1,000 atoms) | 28.5 | ~1 s |
| synthetic 10k atoms | 3.2 | 13 s |
| synthetic 50k atoms | 0.5 | 71 s |

**Instanced benchmark** (`graph-viz-bench`: raw three.js, one InstancedMesh
for nodes + one for bonds — measured after fixing a compile error that proved
the scaffold had never run):

| workload | fps | draw calls |
|---|---|---|
| 10k / 100k / 500k nodes | 60 (vsync cap) | 3 |
| 1M nodes + 750k edges | 30 | 3 |

Same GPU, same night, same browser. At ~30 fps the legacy architecture
handles **1 thousand** nodes and instancing handles **1 million** — a 1000×
architecture gap. The 4090 was never the bottleneck; per-object scene graphs
and React reconciliation were.

## Decisions

1. **Instanced WebGL via three.js. WebGPU rejected** — the bench shows WebGL
   instancing exceeds the target (10⁶ nodes interactive) with naive geometry;
   WebGPU buys nothing this project needs and costs browser-support risk.
2. **No React anywhere near the render path.** The app is vanilla TypeScript;
   the UI is a thin DOM overlay. (React wasn't just slow per-atom — it added
   a reconciliation layer whose job the columnar data path does for free.)
3. **Columnar end-to-end.** PGIF columns → typed arrays → instance buffers,
   no per-node objects at any stage (`src/graph.ts` is the single in-memory
   shape). Picking is a ray/sphere sweep over the flat arrays (~ms at 1M) —
   no raycaster, no BVH.
4. **Rebuild instance buffers on visibility change** (O(N), a few ms at 1M)
   instead of per-instance visibility state — the columnar graph stays the
   single source of truth.
5. **Adaptive geometry detail** by node count (24×16 sphere segments ≤20k,
   8×6 ≥200k) — the 1M case renders 87.5M triangles at 30 fps; impostor
   billboards are the known next lever if 1M@60 ever matters. Deliberately
   not built: no current workload needs it.
6. **Labels are opt-in and singular**: the picked node gets a DOM label
   projected per frame. Never N in-scene text objects (a first-order cause of
   the legacy collapse).
7. **Layout is a Web Worker**, only for graphs that arrive without positions:
   seeded deterministic force-directed (springs + sampled repulsion), typed
   arrays in, transferable Float32Array out. Positioned graphs never pay it.
8. **Own repo** (`vsletten/graph-viz`), not `dissertation/viewer` — settles
   the mission-control open question. The generalization to arbitrary
   property graphs (PGIF, hashed type styling, layout for position-less
   graphs) has no business in the dissertation repo; the legacy viewer stays
   untouched there as a reference until PGIF output lands in the Rust port
   (kmc M9) and it can be archived.
9. **CERIUS2 is an import path, not a format.** The MSI loader
   (grammar-compatible with the legacy parser, quirks preserved: shared
   atom/bond id space, object-id bond references, dangling bonds dropped)
   exists so old files open; `npm run msi2pgif` converts them out. Nothing
   new ever writes MSI.

## Measured results (this app, same hardware)

| workload | fps | draw calls | notes |
|---|---|---|---|
| kaolinite `start.msi` (1,000 nodes, 1,426 edges) | 60 (vsync) | 2 | legacy did 28.5 |
| lattice 100k | 60 (vsync) | 2 | 26.4M tris |
| lattice 1M + 750k edges | 30 | 2 | 87.5M tris |
| social 10k, no positions | 60 after ~2 s layout | 2 | worker layout |

Verified live 2026-07-12 (unit suite 13 green, tsc strict clean, zero console
errors): MSI + PGIF-JSON + PGIF-binary loads (binary through the real
drag-drop path), pick-with-occlusion, legend visibility toggles
(10,000 → 4,045 visible and back), inspector + projected label, both PGIF
exports.
