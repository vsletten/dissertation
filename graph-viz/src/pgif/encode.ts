/**
 * Internal Graph → PGIF, both serializations. The binary encoder writes the
 * bespoke container (magic + u32 header length + JSON header + aligned blob)
 * with every numeric column as a raw little-endian buffer.
 */

import type { Graph, GraphColumn } from '../graph';
import { BLOB_ALIGN, PGIF_MAGIC, type PgifColumn, type PgifDocument } from './types';

function columnToJson(col: GraphColumn): PgifColumn {
  if (col.type === 'str') return { type: 'str', data: col.data as string[] };
  const data = Array.from(col.data as ArrayLike<number>);
  return col.type === 'categorical'
    ? { type: 'categorical', data, dict: col.dict ?? [] }
    : { type: col.type as PgifColumn['type'], data };
}

function toDocument(g: Graph): PgifDocument {
  const nodeColumns: Record<string, PgifColumn> = {};
  for (const [name, col] of Object.entries(g.nodeColumns)) nodeColumns[name] = columnToJson(col);
  // Loaders keep x/y/z/type in nodeColumns; graphs built programmatically may
  // only carry the packed arrays — reconstitute their columns for serialization.
  if (g.positions && !nodeColumns['x']) {
    const xs = new Array<number>(g.count), ys = new Array<number>(g.count), zs = new Array<number>(g.count);
    for (let i = 0; i < g.count; i++) {
      xs[i] = g.positions[i * 3];
      ys[i] = g.positions[i * 3 + 1];
      zs[i] = g.positions[i * 3 + 2];
    }
    nodeColumns['x'] = { type: 'f32', data: xs };
    nodeColumns['y'] = { type: 'f32', data: ys };
    nodeColumns['z'] = { type: 'f32', data: zs };
  }
  if (!nodeColumns['type'] && g.typeDict.length > 1) {
    nodeColumns['type'] = { type: 'categorical', dict: g.typeDict, data: Array.from(g.typeIndex) };
  }
  const edgeColumns: Record<string, PgifColumn> = {};
  for (const [name, col] of Object.entries(g.edgeColumns)) edgeColumns[name] = columnToJson(col);

  return {
    pgif: 1,
    meta: { ...g.meta, directed: g.directed },
    nodes: { count: g.count, ids: Array.from(g.ids), columns: nodeColumns },
    edges: {
      count: g.edgeCount,
      src: Array.from(g.edgeSrc),
      dst: Array.from(g.edgeDst),
      columns: Object.keys(edgeColumns).length ? edgeColumns : undefined,
    },
  };
}

export function encodeJson(g: Graph): string {
  return JSON.stringify(toDocument(g));
}

type NumericCtor =
  | Float32ArrayConstructor
  | Float64ArrayConstructor
  | Int32ArrayConstructor
  | Uint32ArrayConstructor
  | Int8ArrayConstructor
  | Uint8ArrayConstructor;

function ctorFor(type: string): NumericCtor {
  switch (type) {
    case 'f64': return Float64Array;
    case 'f32': return Float32Array;
    case 'u32': return Uint32Array;
    case 'i8': return Int8Array;
    case 'bool': return Uint8Array;
    default: return Int32Array; // i32 + categorical
  }
}

interface BlobSlot {
  bytes: Uint8Array;
  elemCount: number;
  assign: (ref: { offset: number; len: number }) => void;
}

export function encodeBinary(g: Graph): Uint8Array {
  const doc = toDocument(g);
  const slots: BlobSlot[] = [];

  const pushArray = (arr: Float32Array | Float64Array | Int32Array | Uint32Array | Int8Array | Uint8Array, assign: BlobSlot['assign']) => {
    slots.push({
      bytes: new Uint8Array(arr.buffer, arr.byteOffset, arr.byteLength),
      elemCount: arr.length,
      assign,
    });
  };

  const pushColumn = (col: PgifColumn, len: number) => {
    if (col.type === 'str' || !Array.isArray(col.data)) return; // strings stay inline
    const ctor = ctorFor(col.type);
    const arr = new ctor(len);
    const src = col.data as (number | boolean)[];
    for (let i = 0; i < len; i++) arr[i] = Number(src[i]);
    pushArray(arr, ref => {
      col.data = ref;
    });
  };

  pushArray(Int32Array.from(doc.nodes.ids as number[]), ref => {
    doc.nodes.ids = ref;
  });
  for (const col of Object.values(doc.nodes.columns)) pushColumn(col, g.count);
  pushArray(Uint32Array.from(doc.edges.src as number[]), ref => {
    doc.edges.src = ref;
  });
  pushArray(Uint32Array.from(doc.edges.dst as number[]), ref => {
    doc.edges.dst = ref;
  });
  for (const col of Object.values(doc.edges.columns ?? {})) pushColumn(col, g.edgeCount);

  // lay out the blob, 8-byte aligned per slot
  const offsets: number[] = [];
  let cursor = 0;
  for (const slot of slots) {
    cursor = Math.ceil(cursor / BLOB_ALIGN) * BLOB_ALIGN;
    offsets.push(cursor);
    slot.assign({ offset: cursor, len: slot.elemCount });
    cursor += slot.bytes.byteLength;
  }
  const blobLen = cursor;

  const headerBytes = new TextEncoder().encode(JSON.stringify(doc));
  const headerStart = PGIF_MAGIC.length + 4;
  const blobStart = Math.ceil((headerStart + headerBytes.length) / BLOB_ALIGN) * BLOB_ALIGN;
  const out = new Uint8Array(blobStart + blobLen);
  out.set(PGIF_MAGIC, 0);
  new DataView(out.buffer).setUint32(PGIF_MAGIC.length, headerBytes.length, true);
  out.set(headerBytes, headerStart);
  slots.forEach((slot, i) => out.set(slot.bytes, blobStart + offsets[i]));
  return out;
}
