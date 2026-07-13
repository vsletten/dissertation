import { describe, expect, it } from 'vitest';

import { generateLattice, generateSocial } from '../src/demo';
import { loadPgif } from '../src/pgif/decode';
import { encodeBinary, encodeJson } from '../src/pgif/encode';
import { PGIF_MAGIC } from '../src/pgif/types';

const SMALL_JSON = JSON.stringify({
  pgif: 1,
  meta: { producer: 'test', kind: 'fixture' },
  nodes: {
    count: 3,
    ids: [2, 3, 6],
    columns: {
      x: { type: 'f32', data: [1.17181, 1.30217, 0.24775] },
      y: { type: 'f32', data: [4.38134, 7.2829, 3.10361] },
      z: { type: 'f32', data: [2.84548, 2.58231, 0.225645] },
      type: { type: 'categorical', dict: ['Al', 'Si', 'O', 'OH'], data: [0, 0, 1] },
      state: { type: 'i32', data: [101, 101, 201] },
      label: { type: 'str', data: ['Al0', 'Al1', 'Si4'] },
    },
  },
  edges: {
    count: 2,
    src: [0, 0],
    dst: [2, 1],
    columns: { seam: { type: 'bool', data: [false, true] } },
  },
});

describe('pgif json decode', () => {
  it('decodes the spec example shape', () => {
    const g = loadPgif(SMALL_JSON);
    expect(g.count).toBe(3);
    expect(g.edgeCount).toBe(2);
    expect(Array.from(g.ids)).toEqual([2, 3, 6]);
    expect(g.positions).not.toBeNull();
    expect(g.positions![0]).toBeCloseTo(1.17181, 4);
    expect(g.typeDict).toEqual(['Al', 'Si', 'O', 'OH']);
    expect(Array.from(g.typeIndex)).toEqual([0, 0, 1]);
    expect(Array.from(g.edgeSrc)).toEqual([0, 0]);
    expect(Array.from(g.edgeDst)).toEqual([2, 1]);
    expect(g.nodeColumns['state'].data[1]).toBe(101);
    expect((g.nodeColumns['label'].data as string[])[2]).toBe('Si4');
    expect(g.edgeColumns['seam'].data[1]).toBe(1); // bool → u8
    expect(g.directed).toBe(false);
  });

  it('defaults ids to index identity and type to a single bucket', () => {
    const g = loadPgif(
      JSON.stringify({
        pgif: 1,
        meta: {},
        nodes: { count: 2, columns: {} },
        edges: { count: 1, src: [0], dst: [1] },
      }),
    );
    expect(Array.from(g.ids)).toEqual([0, 1]);
    expect(g.typeDict).toEqual(['node']);
    expect(g.positions).toBeNull();
  });

  it('rejects wrong version, bad counts, and out-of-range edges', () => {
    expect(() => loadPgif(JSON.stringify({ pgif: 2, meta: {}, nodes: { count: 0, columns: {} }, edges: { count: 0, src: [], dst: [] } }))).toThrow(/version/);
    expect(() =>
      loadPgif(JSON.stringify({ pgif: 1, meta: {}, nodes: { count: 1, columns: {} }, edges: { count: 1, src: [0], dst: [5] } })),
    ).toThrow(/out of range/);
    expect(() =>
      loadPgif(JSON.stringify({ pgif: 1, meta: {}, nodes: { count: 1, columns: {} }, edges: { count: 2, src: [0], dst: [0] } })),
    ).toThrow(/edges.count|lengths/);
  });
});

describe('pgif binary roundtrip', () => {
  it('encodes with magic and decodes back bit-equal', () => {
    const g = loadPgif(SMALL_JSON);
    const bytes = encodeBinary(g);
    for (let i = 0; i < PGIF_MAGIC.length; i++) expect(bytes[i]).toBe(PGIF_MAGIC[i]);

    const g2 = loadPgif(bytes);
    expect(g2.count).toBe(g.count);
    expect(g2.edgeCount).toBe(g.edgeCount);
    expect(Array.from(g2.ids)).toEqual(Array.from(g.ids));
    expect(Array.from(g2.typeIndex)).toEqual(Array.from(g.typeIndex));
    expect(g2.typeDict).toEqual(g.typeDict);
    expect(Array.from(g2.edgeSrc)).toEqual(Array.from(g.edgeSrc));
    expect(Array.from(g2.edgeDst)).toEqual(Array.from(g.edgeDst));
    // f32 positions survive exactly (no re-quantization on a f32→f32 trip)
    expect(Array.from(g2.positions!)).toEqual(Array.from(g.positions!));
    // str column stays inline and intact
    expect(g2.nodeColumns['label'].data).toEqual(g.nodeColumns['label'].data);
    // i32 model column survives
    expect(Array.from(g2.nodeColumns['state'].data as Int32Array)).toEqual([101, 101, 201]);
  });

  it('json roundtrip preserves the graph too', () => {
    const g = loadPgif(SMALL_JSON);
    const g2 = loadPgif(encodeJson(g));
    expect(g2.count).toBe(3);
    expect(Array.from(g2.edgeDst)).toEqual([2, 1]);
    expect(g2.typeDict).toEqual(['Al', 'Si', 'O', 'OH']);
  });

  it('roundtrips a 10k demo graph through binary', () => {
    const g = generateLattice(10_000);
    const g2 = loadPgif(encodeBinary(g));
    expect(g2.count).toBe(10_000);
    expect(g2.edgeCount).toBe(g.edgeCount);
    expect(g2.positions![12345]).toBeCloseTo(g.positions![12345], 5);
    expect(Array.from(g2.typeIndex.slice(0, 50))).toEqual(Array.from(g.typeIndex.slice(0, 50)));
  });
});

describe('demo generators', () => {
  it('lattice is deterministic and edge-valid', () => {
    const a = generateLattice(1000);
    const b = generateLattice(1000);
    expect(Array.from(a.positions!.slice(0, 30))).toEqual(Array.from(b.positions!.slice(0, 30)));
    for (let e = 0; e < a.edgeCount; e++) {
      expect(a.edgeSrc[e]).toBeLessThan(1000);
      expect(a.edgeDst[e]).toBeLessThan(1000);
      expect(a.edgeSrc[e]).not.toBe(a.edgeDst[e]);
    }
  });

  it('social graph has no positions and valid preferential edges', () => {
    const g = generateSocial(500);
    expect(g.positions).toBeNull();
    expect(g.count).toBe(500);
    expect(g.edgeCount).toBeGreaterThanOrEqual(499);
    for (let e = 0; e < g.edgeCount; e++) {
      expect(g.edgeDst[e]).toBeLessThan(g.edgeSrc[e] + 1); // attaches to earlier nodes
    }
  });
});
