/**
 * Node picking on the columnar arrays directly — no three.js raycaster, no
 * per-instance objects. A ray/sphere test over N nodes in flat typed arrays
 * runs in a few ms at 1M nodes, and honoring the same visibility/radius rules
 * the renderer uses keeps hits WYSIWYG.
 */

import type { PerspectiveCamera } from 'three';
import { Vector3 } from 'three';

import type { Graph } from '../graph';
import type { TypeStyle } from './style';

const origin = new Vector3();
const dir = new Vector3();

export function pickNode(
  ndcX: number,
  ndcY: number,
  camera: PerspectiveCamera,
  graph: Graph,
  styles: TypeStyle[],
  radiusScale: number,
  isHidden: (type: number) => boolean,
): number | null {
  if (!graph.positions) return null;
  origin.setFromMatrixPosition(camera.matrixWorld);
  dir.set(ndcX, ndcY, 0.5).unproject(camera).sub(origin).normalize();

  const pos = graph.positions;
  const ox = origin.x, oy = origin.y, oz = origin.z;
  const dx = dir.x, dy = dir.y, dz = dir.z;

  let best = -1;
  let bestT = Infinity;
  for (let i = 0; i < graph.count; i++) {
    const t = graph.typeIndex[i];
    if (isHidden(t)) continue;
    const r = styles[t].radius * radiusScale;
    const px = pos[i * 3] - ox, py = pos[i * 3 + 1] - oy, pz = pos[i * 3 + 2] - oz;
    const tc = px * dx + py * dy + pz * dz; // distance along ray to closest point
    if (tc < 0 || tc > bestT) continue;
    const d2 = px * px + py * py + pz * pz - tc * tc;
    if (d2 <= r * r) {
      best = i;
      bestT = tc;
    }
  }
  return best >= 0 ? best : null;
}
