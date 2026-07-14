# graph-viz Recon — TASK-001

## Repo location
- Primary candidate: `/mnt/data/vsletten/src/vsletten/dissertation/main/viewer`
- Git remote: `https://github.com/vsletten/dissertation.git` (`origin`, confirmed in `git remote -v`)
- Confidence: high (README + CLAUDE docs explicitly describe Cerius2 molecular visualization and kaolinite context)

## Why this repo matches the objective
- The repository is a browser viewer for Cerius2 `.msi` molecular model files used by the KMC dissertation on kaolin mineral surfaces (the mineral-surface/graph visualization app profile in TASK-001).
- It is not the emergent-graphs `graph-viz` app (`/mnt/data/vsletten/src/vsletten/emergent-graphs/main/graph-viz`), which is a React + D3 knowledge-graph demo and does not handle CERIUS, kaolinite, or mineral surfaces.

## Tech stack
- Framework: React 19 + TypeScript
- Renderer: `three` via `@react-three/fiber` and `@react-three/drei`
- Build tool: Vite
- Data ingestion UI: drag-and-drop + file input (`FileUploader.tsx`) in browser
- Layout/state logic: compute center/camera distance from model bounds in `MoleculeViewer.tsx` (no custom force-layout; uses 3D scene transforms)

## CERIUS 2 ingestion path (code-read)
1. User uploads a CERIUS2 file in `App.tsx` via `FileUploader`.
2. `App.tsx` calls `parseModelFile(content)` from `src/utils/cerius2Parser.ts`.
3. `parseModelFile` scans MSI text line-by-line (`split('\n')` + trim/filter), tracks `Atom`/`Bond` object contexts (`(N Atom`, `(N Bond)`), and extracts:
   - Atoms: `id`, `element`, `position`, `label`, `charge`
   - Bonds: `Atom1`, `Atom2`
4. Parsed `{ atoms, bonds }` is stored in React state and forwarded to `MoleculeViewer.tsx`.
5. `MoleculeViewer` centers positions, builds a lookup map from atom IDs, renders:
   - Atoms as spheres (`AtomSphere`)
   - Bonds as cylinders / lines (`Bond`)
   - Optional labels (`Billboard + Text`)

## Rough LOC estimate (from source files)
- `viewer/src/main.tsx`: 8
- `viewer/src/App.tsx`: 282
- `viewer/src/components/FileUploader.tsx`: 201
- `viewer/src/components/MoleculeViewer.tsx`: 614
- `viewer/src/components/ModelInfo.tsx`: 358
- `viewer/src/utils/cerius2Parser.ts`: 195
- `viewer/src/utils/elementColors.ts`: 112
- `viewer/src/utils/testParser.ts`: 78
- `viewer/src/types/cerius2.ts`: 31
- Approximate source LOC: **1,879**

## Top 3 performance bottlenecks visible from code read
1. **Full in-memory conversion of entire MSI to React objects** in `parseModelFile` (`O(lines)` and pushes all atoms/bonds up front), with no paging or streaming for very large snapshots.
2. **No render instancing / batching strategy**: each atom and bond becomes an individual React/Three component/sphere/cylinder, which can saturate render-thread and GPU submission overhead at scale.
3. **Dense labels and text in scene by default** (`showLabels` defaults true; text labels use drei `Text` for every atom), plus per-frame object updates from React-rendered groups, which compounds cost as node count rises.

