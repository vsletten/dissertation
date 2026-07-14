# graph-viz — diagnosis + rebuild, completed 2026-07-12

Victor's Sunday-night directive: all graph-viz work done tonight, nothing left.
Done. This doc is the record; the code is the deliverable.

**Repo (answers the open repo-location question): `/mnt/data/vsletten/src/vsletten/graph-viz/main`**
(own repo, local like kmc-rs; `gh repo create vsletten/graph-viz --private --source . --push`
is one command if a remote is wanted). The dissertation repo was not touched;
its viewer stays as-is until kmc M9 emits PGIF, then it can be archived.

## 1. The diagnosis (the planned Saturday session, run Sunday)

All on the 4090, Chrome/WebGL2 via ANGLE, 60 Hz display:

| workload | legacy viewer | instanced bench |
|---|---|---|
| 1k atoms (real kaolinite golden) | 28.5 fps | — |
| 10k | 3.2 fps, 13 s mount | 60 fps (vsync) |
| 50k | 0.5 fps, 71 s mount | — |
| 500k | — | 60 fps (vsync cap) |
| 1M + 750k edges | — | 30 fps, 3 draw calls |

Same GPU, same night: at equal ~30 fps the architecture gap is **1,000 vs
1,000,000 nodes**. The 4090 was never the ceiling; per-atom React components
were. Verdict: instanced WebGL suffices, **WebGPU rejected**. Bench numbers
committed to `graph-viz-bench/main/RESULTS.md` (note recorded there: the
scaffold as delivered had a compile error and had never run; fixed).

## 2. The rebuild (all four brief goals closed)

1. **Performance**: kaolinite 60 fps (legacy 28.5), lattice 1M @ 30 fps /
   2 draw calls, 500k @ vsync. Columnar → InstancedMesh; typed-array picking;
   worker force layout for position-less graphs.
2. **UI**: dark glass panels, type legend with visibility toggles, node
   inspector with every column value, projected pick label, live HUD,
   drag-drop + demo menu. Screenshots in `graph-viz/main/docs/img/`.
3. **Generalized**: any PGIF property graph — chemistry gets the element
   palette, arbitrary types get stable hashed colors, graphs without
   positions get laid out. Demo: 10k-node social graph, 60 fps.
4. **CERIUS2 killed at the consumer**: MSI is import-only (grammar-compatible
   loader + `npm run msi2pgif` converter; golden kaolinite → binary PGIF is
   5× smaller). Producer side is kmc-port M9, which now has a finalized spec
   to target: **PGIF v0 is normative in `graph-viz/main/docs/PGIF.md`**
   (draft §7 questions decided: edges by node index, bespoke binary
   container, `meta.directed`, trajectories deferred to v1).

Quality: 13 vitest tests (codecs round-trip incl. binary, MSI grammar +
dangling-bond handling, layout determinism), tsc strict clean, production
build clean, zero console errors. Verified live end-to-end via chrome-devtools
(all three load paths incl. drag-drop of binary PGIF, picking with occlusion,
legend toggles, exports).

## 3. Left over for other projects (not graph-viz work)

- **kmc-port M9**: the Rust model emits PGIF (producer notes in the spec §6).
- Optional someday-items recorded in `graph-viz/main/docs/ARCHITECTURE.md`
  (impostor billboards for 1M@60, trajectory PGIF) — explicitly decided
  against building now; there is no current workload for them.

Dev server left running for morning inspection: `http://localhost:5197`
(`?demo=lattice:1m`, `?demo=social:10000`; restart: `npm run dev` in the repo).
