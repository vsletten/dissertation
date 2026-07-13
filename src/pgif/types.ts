/**
 * PGIF v0 — Property-Graph Interchange Format, logical model.
 *
 * Finalized from the mission-control draft (projects/kmc/05-pgif-draft.md §8)
 * with the §7 open questions decided; see docs/PGIF.md for the normative spec.
 * Decisions: edges reference node INDICES (not ids) in both serializations;
 * binary container is bespoke magic + JSON header + aligned blob;
 * `meta.directed` defaults to false; trajectories are out of v0.
 */

export type ColumnType = 'f32' | 'f64' | 'i32' | 'u32' | 'i8' | 'bool' | 'str' | 'categorical';

/** A column in the JSON serialization (binary uses BufferRef for numeric data). */
export interface PgifColumn {
  type: ColumnType;
  /** JSON form: plain array. Binary form: replaced by a BufferRef into the blob. */
  data: number[] | string[] | boolean[] | BufferRef;
  /** Only for `categorical`: data holds i32 indices into this dict. */
  dict?: string[];
}

/** Pointer into the binary blob region (element count, not bytes). */
export interface BufferRef {
  offset: number;
  len: number;
}

export interface PgifNodes {
  count: number;
  /** Stable external ids; defaults to index identity when omitted. */
  ids?: number[] | BufferRef;
  columns: Record<string, PgifColumn>;
}

export interface PgifEdges {
  count: number;
  /** Node indices (0-based positions in the node arrays), NOT external ids. */
  src: number[] | BufferRef;
  dst: number[] | BufferRef;
  columns?: Record<string, PgifColumn>;
}

export interface PgifDocument {
  pgif: 1;
  meta: Record<string, unknown> & { directed?: boolean };
  nodes: PgifNodes;
  edges: PgifEdges;
}

/** Binary container framing (docs/PGIF.md §4). */
export const PGIF_MAGIC = new Uint8Array([0x50, 0x47, 0x49, 0x46, 0x00, 0x00, 0x00, 0x01]); // "PGIF\0\0\0\1"
export const BLOB_ALIGN = 8;

export function isBufferRef(data: unknown): data is BufferRef {
  return (
    typeof data === 'object' &&
    data !== null &&
    !Array.isArray(data) &&
    typeof (data as BufferRef).offset === 'number' &&
    typeof (data as BufferRef).len === 'number'
  );
}

/** Bytes per element for numeric column types (str/categorical dict stay in JSON). */
export const TYPE_BYTES: Record<Exclude<ColumnType, 'str'>, number> = {
  f32: 4,
  f64: 8,
  i32: 4,
  u32: 4,
  i8: 1,
  bool: 1,
  categorical: 4, // i32 dict indices
};

export const TYPED_ARRAYS = {
  f32: Float32Array,
  f64: Float64Array,
  i32: Int32Array,
  u32: Uint32Array,
  i8: Int8Array,
  bool: Uint8Array,
  categorical: Int32Array,
} as const;
