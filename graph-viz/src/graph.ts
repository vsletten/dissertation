/**
 * The internal columnar graph: what every loader produces and the renderer consumes.
 *
 * Structure-of-arrays throughout — positions and edge endpoints live in typed
 * arrays that go straight into GPU instance buffers. Non-numeric columns are
 * kept for inspection (tooltips/filters), never iterated per-frame.
 */

export interface GraphColumn {
  type: string;
  /** Typed array for numerics, plain array for str. Categorical: Int32Array + dict. */
  data: Float32Array | Float64Array | Int32Array | Uint32Array | Int8Array | Uint8Array | string[];
  dict?: string[];
}

export interface Graph {
  meta: Record<string, unknown>;
  count: number;
  edgeCount: number;
  /** External node ids; identity when the source omitted them. */
  ids: Int32Array;
  /** xyz interleaved [x0,y0,z0, x1,...]; null when the source had no positions. */
  positions: Float32Array | null;
  /** Categorical type per node (index into typeDict); all-zeros + ['node'] when untyped. */
  typeIndex: Int32Array;
  typeDict: string[];
  /** Edge endpoints as node indices. */
  edgeSrc: Uint32Array;
  edgeDst: Uint32Array;
  /** Every node column by name (including x/y/z/type when present), for inspection. */
  nodeColumns: Record<string, GraphColumn>;
  edgeColumns: Record<string, GraphColumn>;
  directed: boolean;
}

export interface Bounds {
  min: [number, number, number];
  max: [number, number, number];
  center: [number, number, number];
  extent: number;
}

export function computeBounds(positions: Float32Array, count: number): Bounds {
  let minX = Infinity, minY = Infinity, minZ = Infinity;
  let maxX = -Infinity, maxY = -Infinity, maxZ = -Infinity;
  for (let i = 0; i < count; i++) {
    const x = positions[i * 3], y = positions[i * 3 + 1], z = positions[i * 3 + 2];
    if (x < minX) minX = x;
    if (x > maxX) maxX = x;
    if (y < minY) minY = y;
    if (y > maxY) maxY = y;
    if (z < minZ) minZ = z;
    if (z > maxZ) maxZ = z;
  }
  if (count === 0) {
    minX = minY = minZ = -1;
    maxX = maxY = maxZ = 1;
  }
  const center: [number, number, number] = [(minX + maxX) / 2, (minY + maxY) / 2, (minZ + maxZ) / 2];
  const extent = Math.max(maxX - minX, maxY - minY, maxZ - minZ, 1e-6);
  return { min: [minX, minY, minZ], max: [maxX, maxY, maxZ], center, extent };
}

/** Median edge length — used to scale node radius so any graph reads well. */
export function medianEdgeLength(g: Graph): number {
  if (!g.positions || g.edgeCount === 0) return 1;
  const sample = Math.min(g.edgeCount, 5000);
  const step = Math.max(1, Math.floor(g.edgeCount / sample));
  const lengths: number[] = [];
  for (let e = 0; e < g.edgeCount; e += step) {
    const a = g.edgeSrc[e] * 3, b = g.edgeDst[e] * 3;
    const dx = g.positions[a] - g.positions[b];
    const dy = g.positions[a + 1] - g.positions[b + 1];
    const dz = g.positions[a + 2] - g.positions[b + 2];
    lengths.push(Math.sqrt(dx * dx + dy * dy + dz * dz));
  }
  lengths.sort((p, q) => p - q);
  return lengths[Math.floor(lengths.length / 2)] || 1;
}

/** All column values for one node, dicts resolved — the inspector payload. */
export function nodeProperties(g: Graph, index: number): Record<string, unknown> {
  const out: Record<string, unknown> = { '#index': index, id: g.ids[index] };
  for (const [name, col] of Object.entries(g.nodeColumns)) {
    if (col.type === 'categorical' && col.dict) {
      out[name] = col.dict[(col.data as Int32Array)[index]] ?? null;
    } else if (col.type === 'bool') {
      out[name] = !!(col.data as Uint8Array)[index];
    } else {
      out[name] = (col.data as ArrayLike<unknown>)[index];
    }
  }
  return out;
}

/** Node counts per type — legend data. */
export function typeCounts(g: Graph): Map<number, number> {
  const counts = new Map<number, number>();
  for (let i = 0; i < g.count; i++) {
    const t = g.typeIndex[i];
    counts.set(t, (counts.get(t) ?? 0) + 1);
  }
  return counts;
}
