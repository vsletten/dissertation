/**
 * Legacy CERIUS2 MSI reader — the migration path off the old format.
 *
 * Grammar-compatible with the dissertation viewer's parser (atoms via `(N Atom`,
 * ACL element extraction, `(A D XYZ (...))`, bonds via Atom1/Atom2 object-id
 * references), but emits the columnar Graph directly. MSI quirks handled:
 * atom and bond ids share one integer space, and bonds reference atom OBJECT
 * ids, which we remap to node indices (PGIF's edge convention).
 */

import type { Graph, GraphColumn } from '../graph';

const ATOM_START = /^\((\d+)\s+Atom/;
const BOND_START = /^\((\d+)\s+Bond/;
const QUOTED = /"([^"]+)"/;
const ACL_ELEMENT = /\d+\s+([A-Za-z]+)/;
const XYZ = /\(([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\)/;
const INT_ATTR = /(\d+)\s*\)?\s*$/;
const CHARGE = /Charge\s+([-\d.eE+]+)/;

export function loadMsi(content: string): Graph {
  const xs: number[] = [], ys: number[] = [], zs: number[] = [];
  const objectIds: number[] = [];
  const elements: string[] = [];
  const labels: string[] = [];
  const charges: number[] = [];
  const bondA: number[] = [], bondB: number[] = [];

  let mode: 'atom' | 'bond' | null = null;
  let el = '', label = '';
  let x = 0, y = 0, z = 0, charge = 0, objectId = 0;
  let atom1 = 0, atom2 = 0;

  const flush = () => {
    if (mode === 'atom') {
      objectIds.push(objectId);
      elements.push(el || 'X');
      labels.push(label);
      charges.push(charge);
      xs.push(x); ys.push(y); zs.push(z);
    } else if (mode === 'bond' && atom1 > 0 && atom2 > 0) {
      bondA.push(atom1);
      bondB.push(atom2);
    }
    mode = null;
  };

  for (const raw of content.split('\n')) {
    const line = raw.trim();
    if (!line || line.startsWith('#')) continue;

    const atomStart = line.match(ATOM_START);
    const bondStart = atomStart ? null : line.match(BOND_START);
    if (atomStart || bondStart) {
      flush();
      mode = atomStart ? 'atom' : 'bond';
      objectId = parseInt((atomStart ?? bondStart)![1], 10);
      el = ''; label = ''; x = y = z = 0; charge = 0;
      atom1 = 0; atom2 = 0;
      continue;
    }
    if (line === ')') {
      flush();
      continue;
    }
    if (mode === 'atom') {
      if (line.includes('ACL')) {
        const q = line.match(QUOTED);
        const em = q?.[1].match(ACL_ELEMENT);
        if (em) el = em[1];
      } else if (line.includes('(A D XYZ')) {
        const m = line.match(XYZ);
        if (m) {
          x = parseFloat(m[1]); y = parseFloat(m[2]); z = parseFloat(m[3]);
        }
      } else if (line.includes('(A C Label')) {
        const m = line.match(QUOTED);
        if (m) label = m[1];
      } else if (line.includes('(A F Charge')) {
        const m = line.match(CHARGE);
        if (m) charge = parseFloat(m[1]);
      }
    } else if (mode === 'bond') {
      if (line.includes('(A O Atom1')) {
        const m = line.match(INT_ATTR);
        if (m) atom1 = parseInt(m[1], 10);
      } else if (line.includes('(A O Atom2')) {
        const m = line.match(INT_ATTR);
        if (m) atom2 = parseInt(m[1], 10);
      }
    }
  }
  flush();

  const n = objectIds.length;

  // element strings → categorical type column
  const typeDict: string[] = [];
  const dictIndex = new Map<string, number>();
  const typeIndex = new Int32Array(n);
  for (let i = 0; i < n; i++) {
    let t = dictIndex.get(elements[i]);
    if (t === undefined) {
      t = typeDict.length;
      typeDict.push(elements[i]);
      dictIndex.set(elements[i], t);
    }
    typeIndex[i] = t;
  }

  // bonds reference atom object ids → node indices; drop dangling references
  const byObjectId = new Map<number, number>();
  for (let i = 0; i < n; i++) byObjectId.set(objectIds[i], i);
  const src: number[] = [], dst: number[] = [];
  for (let e = 0; e < bondA.length; e++) {
    const a = byObjectId.get(bondA[e]);
    const b = byObjectId.get(bondB[e]);
    if (a !== undefined && b !== undefined) {
      src.push(a);
      dst.push(b);
    }
  }

  const positions = new Float32Array(n * 3);
  for (let i = 0; i < n; i++) {
    positions[i * 3] = xs[i];
    positions[i * 3 + 1] = ys[i];
    positions[i * 3 + 2] = zs[i];
  }

  const nodeColumns: Record<string, GraphColumn> = {
    x: { type: 'f32', data: Float32Array.from(xs) },
    y: { type: 'f32', data: Float32Array.from(ys) },
    z: { type: 'f32', data: Float32Array.from(zs) },
    type: { type: 'categorical', data: typeIndex, dict: typeDict },
    label: { type: 'str', data: labels },
  };
  if (charges.some(c => c !== 0)) nodeColumns['charge'] = { type: 'f32', data: Float32Array.from(charges) };

  return {
    meta: { producer: 'graph-viz msi loader', kind: 'cerius2-msi' },
    count: n,
    edgeCount: src.length,
    ids: Int32Array.from(objectIds),
    positions,
    typeIndex,
    typeDict,
    edgeSrc: Uint32Array.from(src),
    edgeDst: Uint32Array.from(dst),
    nodeColumns,
    edgeColumns: {},
    directed: false,
  };
}
