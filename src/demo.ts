/**
 * Built-in demo graphs: a positioned lattice (the bench workload, sizes to 1M)
 * and a position-less scale-free graph that exercises the layout engine.
 * Both are deterministic (mulberry32) so numbers are comparable across runs.
 */

import type { Graph, GraphColumn } from './graph';

export const LATTICE_SIZES = { '10k': 10_000, '100k': 100_000, '500k': 500_000, '1m': 1_000_000 } as const;
export type LatticeSize = keyof typeof LATTICE_SIZES;

function mulberry32(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a |= 0;
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

const LATTICE_TYPES = ['Al', 'Si', 'O', 'OH', 'Fe', 'Mg'];

/** Lattice-ish point cloud with ~0.75 edges/node — the bench's synthetic workload. */
export function generateLattice(nodeCount: number, seed = 0x9e37): Graph {
  const rand = mulberry32(seed);
  const side = Math.max(1, Math.round(Math.cbrt(nodeCount)));
  const positions = new Float32Array(nodeCount * 3);
  const typeIndex = new Int32Array(nodeCount);
  const ids = new Int32Array(nodeCount);
  for (let i = 0; i < nodeCount; i++) {
    positions[i * 3] = (i % side) * 2.5 + (rand() - 0.5) * 0.6;
    positions[i * 3 + 1] = (Math.floor(i / side) % side) * 2.5 + (rand() - 0.5) * 0.6;
    positions[i * 3 + 2] = Math.floor(i / (side * side)) * 2.5 + (rand() - 0.5) * 0.6;
    typeIndex[i] = i % LATTICE_TYPES.length;
    ids[i] = i;
  }
  const edgeCount = Math.floor(nodeCount * 0.75);
  const edgeSrc = new Uint32Array(edgeCount);
  const edgeDst = new Uint32Array(edgeCount);
  for (let e = 0; e < edgeCount; e++) {
    // connect near neighbors in index space — index-adjacent is lattice-adjacent
    const a = Math.floor(rand() * nodeCount);
    let b = a + 1 + Math.floor(rand() * 3) * side;
    if (b >= nodeCount) b = a > 0 ? a - 1 : 0;
    edgeSrc[e] = a;
    edgeDst[e] = b === a ? (a + 1) % nodeCount : b;
  }
  return withStandardColumns({
    meta: { producer: 'graph-viz demo', kind: 'lattice' },
    count: nodeCount,
    edgeCount,
    ids,
    positions,
    typeIndex,
    typeDict: [...LATTICE_TYPES],
    edgeSrc,
    edgeDst,
    nodeColumns: {},
    edgeColumns: {},
    directed: false,
  });
}

const SOCIAL_TYPES = ['Person', 'Org', 'Place', 'Topic'];

/** Scale-free-ish graph with NO positions — the layout engine's demo input. */
export function generateSocial(nodeCount: number, seed = 0x51f5): Graph {
  const rand = mulberry32(seed);
  const typeIndex = new Int32Array(nodeCount);
  const ids = new Int32Array(nodeCount);
  const names: string[] = new Array(nodeCount);
  for (let i = 0; i < nodeCount; i++) {
    const t = rand() < 0.6 ? 0 : 1 + Math.floor(rand() * 3);
    typeIndex[i] = t;
    ids[i] = i;
    names[i] = `${SOCIAL_TYPES[t]} ${i}`;
  }
  // preferential attachment: each new node links to 1-3 earlier nodes, biased low
  const src: number[] = [], dst: number[] = [];
  for (let i = 1; i < nodeCount; i++) {
    const links = 1 + Math.floor(rand() * 3);
    for (let l = 0; l < links; l++) {
      const target = Math.floor(Math.pow(rand(), 2.2) * i);
      src.push(i);
      dst.push(target);
    }
  }
  return withStandardColumns({
    meta: { producer: 'graph-viz demo', kind: 'social' },
    count: nodeCount,
    edgeCount: src.length,
    ids,
    positions: null,
    typeIndex,
    typeDict: [...SOCIAL_TYPES],
    edgeSrc: Uint32Array.from(src),
    edgeDst: Uint32Array.from(dst),
    nodeColumns: { name: { type: 'str', data: names } },
    edgeColumns: {},
    directed: false,
  });
}

/** Mirror packed arrays into inspectable columns (x/y/z/type) like loaders do. */
function withStandardColumns(g: Graph): Graph {
  const cols: Record<string, GraphColumn> = { ...g.nodeColumns };
  if (g.positions) {
    const n = g.count;
    const xs = new Float32Array(n), ys = new Float32Array(n), zs = new Float32Array(n);
    for (let i = 0; i < n; i++) {
      xs[i] = g.positions[i * 3];
      ys[i] = g.positions[i * 3 + 1];
      zs[i] = g.positions[i * 3 + 2];
    }
    cols['x'] = { type: 'f32', data: xs };
    cols['y'] = { type: 'f32', data: ys };
    cols['z'] = { type: 'f32', data: zs };
  }
  cols['type'] = { type: 'categorical', data: g.typeIndex, dict: g.typeDict };
  g.nodeColumns = cols;
  return g;
}
