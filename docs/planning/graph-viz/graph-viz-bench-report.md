# Graph-viz instanced rendering benchmark prototype report

## What I built

I created a throwaway benchmark repo at:

- `/mnt/data/vsletten/src/vsletten/graph-viz-bench/main`

The repo contains a Vite + TypeScript app using raw `three.js` with instanced geometry.
No React, no react-three-fiber.

## Repository contents

- `src/main.ts` — benchmark scene and render loop
  - Orbit controls
  - Instanced atoms (spheres)
  - Instanced bonds (cylinders)
  - Rolling 60-frame FPS + frame-time overlay
  - Node/edge/draw-call counts
- `src/pgif.ts` — synthetic PGIF-shaped graph generator
  - Lattice-like coordinates for predictable geometry
  - Categorical `type` column for per-instance atom color
  - 10k / 100k / 500k / 1M size presets
  - ~1.5 edges-per-node target (`~0.75 * N` undirected edges)
- `script/capture-headless.mjs` — optional headless capture path (Playwright)
- `RESULTS.md` — placeholder to capture outputs
- `README.md` — run instructions

## How to run

```bash
cd /mnt/data/vsletten/src/vsletten/graph-viz-bench/main
npm i
npm run dev
```

Open the served page and use the size dropdown (or append `?size=100k`, etc.).

## Headless capture status

I added an optional capture path:

```bash
npm i -D playwright
npm run dev
npm run capture
```

In this environment I did not execute a live benchmark capture run; if you run `npm run capture`,
it will append a timestamped capture section to `RESULTS.md`.

## Comparison hooks for Saturday's session

The benchmark surface is intentionally narrow to keep comparisons straightforward:

- Compare instanced draw-call behavior directly against legacy viewer architecture.
- Compare FPS and frame-time at 10k/100k/500k/1M across same hardware.
- Inspect node/edge scaling and confirm label rendering is off by default.
- Carry results forward into the Saturday diagnosis session before deciding renderer
  migration paths.

Headless metrics should be treated as indicative only if GPU differs from Saturday's
`4090` machine.
