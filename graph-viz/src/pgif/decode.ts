/**
 * PGIF → internal Graph. One entry point (`loadPgif`) sniffs the serialization:
 * bytes starting with the magic are the binary container; anything else is JSON.
 */

import type { Graph, GraphColumn } from '../graph';
import {
  BLOB_ALIGN,
  isBufferRef,
  PGIF_MAGIC,
  TYPE_BYTES,
  TYPED_ARRAYS,
  type BufferRef,
  type ColumnType,
  type PgifColumn,
  type PgifDocument,
} from './types';

export function loadPgif(input: ArrayBuffer | Uint8Array | string): Graph {
  if (typeof input === 'string') return fromDocument(parseJsonDocument(input), null);
  const bytes = input instanceof Uint8Array ? input : new Uint8Array(input);
  if (hasMagic(bytes)) return decodeBinary(bytes);
  return fromDocument(parseJsonDocument(new TextDecoder().decode(bytes)), null);
}

function hasMagic(bytes: Uint8Array): boolean {
  if (bytes.length < PGIF_MAGIC.length) return false;
  for (let i = 0; i < PGIF_MAGIC.length; i++) if (bytes[i] !== PGIF_MAGIC[i]) return false;
  return true;
}

function parseJsonDocument(text: string): PgifDocument {
  const doc = JSON.parse(text) as PgifDocument;
  if (doc.pgif !== 1) throw new Error(`unsupported pgif version: ${String(doc.pgif)}`);
  if (!doc.nodes || typeof doc.nodes.count !== 'number') throw new Error('pgif: missing nodes.count');
  if (!doc.edges || typeof doc.edges.count !== 'number') throw new Error('pgif: missing edges.count');
  return doc;
}

function decodeBinary(bytes: Uint8Array): Graph {
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const headerLen = view.getUint32(PGIF_MAGIC.length, true);
  const headerStart = PGIF_MAGIC.length + 4;
  const headerBytes = bytes.subarray(headerStart, headerStart + headerLen);
  const doc = parseJsonDocument(new TextDecoder().decode(headerBytes));
  const blobStart = align(headerStart + headerLen, BLOB_ALIGN);
  // Copy the blob when the source offset breaks alignment for typed-array views.
  const blob =
    (bytes.byteOffset + blobStart) % BLOB_ALIGN === 0
      ? new Uint8Array(bytes.buffer, bytes.byteOffset + blobStart, bytes.byteLength - blobStart)
      : bytes.slice(blobStart);
  return fromDocument(doc, blob);
}

function align(n: number, a: number): number {
  return Math.ceil(n / a) * a;
}

type NumericArray = InstanceType<(typeof TYPED_ARRAYS)[keyof typeof TYPED_ARRAYS]>;

function numericData(col: { type: ColumnType; data: PgifColumn['data'] }, blob: Uint8Array | null): NumericArray {
  const ctor = TYPED_ARRAYS[col.type as Exclude<ColumnType, 'str'>];
  if (!ctor) throw new Error(`pgif: column type ${col.type} is not numeric`);
  if (isBufferRef(col.data)) {
    if (!blob) throw new Error('pgif: buffer ref in JSON document (no blob)');
    return viewBuffer(blob, col.data, col.type as Exclude<ColumnType, 'str'>);
  }
  const arr = col.data as (number | boolean)[];
  const out = new ctor(arr.length);
  for (let i = 0; i < arr.length; i++) out[i] = Number(arr[i]);
  return out;
}

function viewBuffer(blob: Uint8Array, ref: BufferRef, type: Exclude<ColumnType, 'str'>): NumericArray {
  const ctor = TYPED_ARRAYS[type];
  const byteOffset = blob.byteOffset + ref.offset;
  if (byteOffset % ctor.BYTES_PER_ELEMENT !== 0) {
    // misaligned ref — copy instead of view
    const copy = blob.slice(ref.offset, ref.offset + ref.len * TYPE_BYTES[type]);
    return new ctor(copy.buffer as ArrayBuffer, 0, ref.len);
  }
  return new ctor(blob.buffer as ArrayBuffer, byteOffset, ref.len);
}

function toGraphColumn(col: PgifColumn, blob: Uint8Array | null): GraphColumn {
  if (col.type === 'str') {
    if (isBufferRef(col.data)) throw new Error('pgif: str columns must stay inline in the header');
    return { type: 'str', data: col.data as string[] };
  }
  const data = numericData(col, blob);
  return col.type === 'categorical'
    ? { type: col.type, data: data as Int32Array, dict: col.dict ?? [] }
    : { type: col.type, data };
}

function fromDocument(doc: PgifDocument, blob: Uint8Array | null): Graph {
  const n = doc.nodes.count;
  const m = doc.edges.count;

  const nodeColumns: Record<string, GraphColumn> = {};
  for (const [name, col] of Object.entries(doc.nodes.columns ?? {})) {
    nodeColumns[name] = toGraphColumn(col, blob);
  }
  const edgeColumns: Record<string, GraphColumn> = {};
  for (const [name, col] of Object.entries(doc.edges.columns ?? {})) {
    edgeColumns[name] = toGraphColumn(col, blob);
  }

  // ids default to index identity
  let ids: Int32Array;
  if (doc.nodes.ids) {
    const raw = numericData({ type: 'i32', data: doc.nodes.ids }, blob);
    ids = raw instanceof Int32Array ? raw : Int32Array.from(raw as ArrayLike<number>);
  } else {
    ids = new Int32Array(n);
    for (let i = 0; i < n; i++) ids[i] = i;
  }
  if (ids.length !== n) throw new Error(`pgif: ids length ${ids.length} != nodes.count ${n}`);

  // positions: x/y/z columns → interleaved Float32Array (absent → null)
  const xc = nodeColumns['x'], yc = nodeColumns['y'], zc = nodeColumns['z'];
  let positions: Float32Array | null = null;
  if (xc && yc && zc) {
    positions = new Float32Array(n * 3);
    const xs = xc.data as ArrayLike<number>, ys = yc.data as ArrayLike<number>, zs = zc.data as ArrayLike<number>;
    for (let i = 0; i < n; i++) {
      positions[i * 3] = xs[i];
      positions[i * 3 + 1] = ys[i];
      positions[i * 3 + 2] = zs[i];
    }
  }

  // type: categorical column named `type` (absent → single implicit type)
  let typeIndex: Int32Array;
  let typeDict: string[];
  const tc = nodeColumns['type'];
  if (tc && tc.type === 'categorical' && tc.dict) {
    typeIndex = tc.data as Int32Array;
    typeDict = tc.dict;
  } else {
    typeIndex = new Int32Array(n);
    typeDict = ['node'];
  }

  const srcRaw = numericData({ type: 'u32', data: doc.edges.src }, blob);
  const dstRaw = numericData({ type: 'u32', data: doc.edges.dst }, blob);
  const edgeSrc = srcRaw instanceof Uint32Array ? srcRaw : Uint32Array.from(srcRaw as ArrayLike<number>);
  const edgeDst = dstRaw instanceof Uint32Array ? dstRaw : Uint32Array.from(dstRaw as ArrayLike<number>);
  if (edgeSrc.length !== m || edgeDst.length !== m) {
    throw new Error(`pgif: edge endpoint lengths (${edgeSrc.length}/${edgeDst.length}) != edges.count ${m}`);
  }
  for (let e = 0; e < m; e++) {
    if (edgeSrc[e] >= n || edgeDst[e] >= n) {
      throw new Error(`pgif: edge ${e} references node index out of range (src=${edgeSrc[e]}, dst=${edgeDst[e]}, count=${n})`);
    }
  }

  return {
    meta: doc.meta ?? {},
    count: n,
    edgeCount: m,
    ids,
    positions,
    typeIndex,
    typeDict,
    edgeSrc,
    edgeDst,
    nodeColumns,
    edgeColumns,
    directed: doc.meta?.directed === true,
  };
}
