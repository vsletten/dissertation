/**
 * Force layout for graphs that arrive without positions (pure property graphs).
 *
 * Approximate 3D force-directed: springs along edges, sampled repulsion
 * (k random partners per node per iteration — O(kN) instead of O(N²)), weak
 * center gravity, velocity damping. Deterministic: seeded init on a golden-
 * spiral sphere and a seeded sampler, so the same graph always lands the same.
 * Runs in a worker; posts progress and returns positions as a transferable.
 */

interface LayoutRequest {
  count: number;
  edgeSrc: Uint32Array;
  edgeDst: Uint32Array;
  seed: number;
}

interface LayoutProgress {
  kind: 'progress';
  done: number;
  total: number;
}

interface LayoutResult {
  kind: 'done';
  positions: Float32Array;
}

export type LayoutMessage = LayoutProgress | LayoutResult;

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

export function runLayout(
  req: LayoutRequest,
  onProgress?: (done: number, total: number) => void,
): Float32Array {
  const { count, edgeSrc, edgeDst, seed } = req;
  const rand = mulberry32(seed);
  const pos = new Float32Array(count * 3);
  const vel = new Float32Array(count * 3);

  // golden-spiral sphere init, radius grows with n so density stays sane
  const radius = Math.max(4, Math.cbrt(count) * 2.2);
  const golden = Math.PI * (3 - Math.sqrt(5));
  for (let i = 0; i < count; i++) {
    const t = count > 1 ? i / (count - 1) : 0.5;
    const inclination = Math.acos(1 - 2 * t);
    const azimuth = golden * i;
    const r = radius * (0.6 + 0.4 * rand());
    pos[i * 3] = r * Math.sin(inclination) * Math.cos(azimuth);
    pos[i * 3 + 1] = r * Math.sin(inclination) * Math.sin(azimuth);
    pos[i * 3 + 2] = r * Math.cos(inclination);
  }

  const iterations = count <= 2_000 ? 300 : count <= 10_000 ? 150 : count <= 50_000 ? 80 : 40;
  const springLength = 2.0;
  const springK = 0.06;
  const repulseK = 6.0;
  const gravity = 0.012;
  const damping = 0.85;
  const repulseSamples = 6;
  const maxStep = 1.5;

  for (let iter = 0; iter < iterations; iter++) {
    // springs
    for (let e = 0; e < edgeSrc.length; e++) {
      const a = edgeSrc[e] * 3, b = edgeDst[e] * 3;
      const dx = pos[b] - pos[a], dy = pos[b + 1] - pos[a + 1], dz = pos[b + 2] - pos[a + 2];
      const dist = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-6;
      const f = springK * (dist - springLength) / dist;
      const fx = f * dx, fy = f * dy, fz = f * dz;
      vel[a] += fx; vel[a + 1] += fy; vel[a + 2] += fz;
      vel[b] -= fx; vel[b + 1] -= fy; vel[b + 2] -= fz;
    }
    // sampled repulsion + gravity
    for (let i = 0; i < count; i++) {
      const ia = i * 3;
      for (let s = 0; s < repulseSamples; s++) {
        const j = Math.floor(rand() * count);
        if (j === i) continue;
        const jb = j * 3;
        const dx = pos[ia] - pos[jb], dy = pos[ia + 1] - pos[jb + 1], dz = pos[ia + 2] - pos[jb + 2];
        const d2 = dx * dx + dy * dy + dz * dz + 0.05;
        const f = repulseK / d2;
        const inv = 1 / Math.sqrt(d2);
        vel[ia] += f * dx * inv; vel[ia + 1] += f * dy * inv; vel[ia + 2] += f * dz * inv;
      }
      vel[ia] -= pos[ia] * gravity;
      vel[ia + 1] -= pos[ia + 1] * gravity;
      vel[ia + 2] -= pos[ia + 2] * gravity;
    }
    // integrate
    for (let i = 0; i < count * 3; i += 3) {
      let sx = vel[i], sy = vel[i + 1], sz = vel[i + 2];
      const mag = Math.sqrt(sx * sx + sy * sy + sz * sz);
      if (mag > maxStep) {
        const scale = maxStep / mag;
        sx *= scale; sy *= scale; sz *= scale;
      }
      pos[i] += sx; pos[i + 1] += sy; pos[i + 2] += sz;
      vel[i] *= damping; vel[i + 1] *= damping; vel[i + 2] *= damping;
    }
    onProgress?.(iter + 1, iterations);
  }
  return pos;
}

// worker entry (absent in unit tests, which import runLayout directly)
if (typeof self !== 'undefined' && typeof (self as unknown as Worker).postMessage === 'function' && typeof document === 'undefined') {
  self.onmessage = (ev: MessageEvent<LayoutRequest>) => {
    let lastReport = 0;
    const positions = runLayout(ev.data, (done, total) => {
      const now = Date.now();
      if (now - lastReport > 120 || done === total) {
        lastReport = now;
        (self as unknown as Worker).postMessage({ kind: 'progress', done, total } satisfies LayoutProgress);
      }
    });
    (self as unknown as Worker).postMessage({ kind: 'done', positions } satisfies LayoutResult, [
      positions.buffer,
    ]);
  };
}
