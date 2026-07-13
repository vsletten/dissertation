import { describe, expect, it } from 'vitest';

import { generateSocial } from '../src/demo';
import { computeBounds } from '../src/graph';
import { runLayout } from '../src/layout.worker';

describe('force layout', () => {
  it('produces finite, deterministic, non-degenerate positions', () => {
    const g = generateSocial(400);
    const req = { count: g.count, edgeSrc: g.edgeSrc, edgeDst: g.edgeDst, seed: 0x51f5 };
    const a = runLayout(req);
    const b = runLayout(req);

    expect(a.length).toBe(400 * 3);
    for (let i = 0; i < a.length; i++) {
      expect(Number.isFinite(a[i])).toBe(true);
      expect(a[i]).toBe(b[i]); // deterministic
    }
    // non-degenerate: the cloud has real extent and isn't collapsed to a point
    const bounds = computeBounds(a, 400);
    expect(bounds.extent).toBeGreaterThan(1);
  });

  it('pulls connected nodes closer than the global scale', () => {
    const g = generateSocial(400);
    const pos = runLayout({ count: g.count, edgeSrc: g.edgeSrc, edgeDst: g.edgeDst, seed: 1 });
    const bounds = computeBounds(pos, 400);
    let sum = 0;
    for (let e = 0; e < g.edgeCount; e++) {
      const a = g.edgeSrc[e] * 3, b = g.edgeDst[e] * 3;
      sum += Math.hypot(pos[a] - pos[b], pos[a + 1] - pos[b + 1], pos[a + 2] - pos[b + 2]);
    }
    const meanEdge = sum / g.edgeCount;
    expect(meanEdge).toBeLessThan(bounds.extent / 2); // edges are local structure
  });
});
