# PGIF v0 — Property-Graph Interchange Format

The shared format between producers (the Rust KMC port, any future exporter)
and this viewer. Finalized 2026-07-12 from the mission-control draft
(`projects/kmc/05-pgif-draft.md`); this file is the normative implementer spec.
The reference implementation is `src/pgif/` (decode + encode, both
serializations, round-trip tested).

## 1. Logical model

A PGIF document is a property graph in structure-of-arrays form:

```
Document = {
  pgif: 1,
  meta:  { free-form; `directed?: bool` (default false) is the only spec'd key },
  nodes: { count: N, ids?: [N] (default: index identity), columns: {name → Column} },
  edges: { count: M, src: [M], dst: [M], columns?: {name → Column} }
}
Column = { type: "f32"|"f64"|"i32"|"u32"|"i8"|"bool"|"str"|"categorical",
           data: [...],           // JSON form
           dict?: [str] }         // categorical only: data = i32 indices into dict
```

Conventions the viewer understands (all optional):

- **`x`,`y`,`z` node columns** → positions. Absent → the viewer runs its force
  layout (any pure topological property graph is renderable).
- **`type` categorical node column** → color/size via the style table
  (chemistry names get the classic element palette; anything else gets a
  stable hashed color).
- Every other column is carried and shown in the node inspector.

## 2. Decisions on the draft's §7 open questions

1. **Edges reference node *indices*** (0-based positions in the node arrays),
   in BOTH serializations. `ids` is a display/identity column, not an edge key.
   Rationale: O(1) consumer lookup, no id-map allocation on the hot path; the
   producer already has indices. (The MSI loader remaps legacy object-id bonds
   to indices at load time — dangling references are dropped, not errors.)
2. **Binary container is bespoke** magic + header + blob (§3). No MessagePack,
   no Arrow — dependency-free on both sides, replaceable later.
3. **`meta.directed: bool`, default false.** Undirected edges are emitted once.
4. **Trajectories are out of v0** — and now specified as a **sidecar**, not a
   container change (§7): a producer (petra) needed them, and the design that
   won is static topology + per-event state *deltas*, not per-frame columns —
   deltas are ~100× smaller for KMC (a few sites change per event) and give
   backward scrubbing for free.
5. **`meta` stays free-form** (`producer`, `kind` recommended). Fields get
   promoted to the spec only when two producers share them.

## 3. Binary container (`.pgif.bin`)

```
[ 8 bytes   magic: "PGIF\0\0\0\x01" (50 47 49 46 00 00 00 01) ]
[ u32 LE    header_len ]
[ bytes     header: UTF-8 JSON of the Document, where every numeric column's
            `data` (and nodes.ids / edges.src / edges.dst) is replaced by
            { "offset": B, "len": N } — element offsets into the blob region.
            `str` columns and categorical `dict`s stay inline in the header. ]
[ padding   to the next 8-byte boundary ]
[ blob      concatenated little-endian typed-array buffers, each slot
            8-byte aligned ]
```

Element widths: f32/i32/u32 = 4, f64 = 8, i8/bool = 1, categorical = 4 (i32
dict indices). `offset` is in **bytes** from the blob start; `len` is in
**elements**. Consumers construct `new Float32Array(blob, offset, len)` —
zero-copy, zero per-element parsing. A consumer whose blob lands misaligned
(e.g. an offset view over a larger buffer) must copy the slice instead.

Type mapping on load: `ids` → Int32Array, `src`/`dst` → Uint32Array,
`bool` → Uint8Array (0/1).

## 4. Serialization choice

- `.pgif.json` — human-readable, diff-friendly. Fine to ~10⁵ nodes.
- `.pgif.bin` — the large-graph path. The kaolinite golden file is 5× smaller
  as binary PGIF than as MSI, and a 1M-node graph loads as typed arrays with
  no parse loop.

One consumer entry point (`loadPgif`) sniffs the magic and branches only at
the byte-decoding layer.

## 5. Validation rules (what the reference decoder enforces)

- `pgif === 1`; `nodes.count`/`edges.count` present.
- Column/ids lengths must equal their `count`.
- Every edge endpoint `< nodes.count` (indices, see §2.1).
- `str` columns may not be buffer refs (they stay in the header).
- Unknown `meta` keys and unknown columns are preserved, never errors.

## 6a. Trajectory sidecar (petra_traj v1)

A PGIF snapshot plus a JSONL event log make a replayable KMC trajectory.
Reference producer: `petra/crates/petra-io`; reference consumer:
`src/traj.ts`. The snapshot must contain **every** lattice site (vacant
included — occupancy flips during playback; the `vacant` display type
renders at radius 0) and a categorical `state` node column whose dict is
indexed by the dense state ids the log uses.

**Log format** — line 1 is a header object, every further line one event:

```
{"petra_traj":1, "deck":str, "seed":int, "n_sites":int,
 "states":[str],       // state names, index = state id (must match the
                       // snapshot's `state` dict)
 "state_types":[str],  // display type per state id (element name | "vacant")
 "reactions":[str]}    // reaction names, index = event rxn id
[step, time, rxn, [[site, old, new], ...]]   // actual state changes only
```

`old` is mandatory: it is what makes backward scrubbing a pure replay
(apply deltas inverted, in reverse order) instead of a re-simulation.
The consumer applies deltas to the `state` column and derives node display
type via `state_types`; types unseen at t=0 are appended to the type dict
at attach time.

**Edge conventions:** topology is static across a trajectory. Bonds that
wrap a periodic boundary carry `seam: true` in a bool edge column; the
viewer allocates them instance slots but never draws them (a wrap bond is
topology, not geometry).

## 6. Producer notes (for the Rust KMC port — milestone M9)

Emitting PGIF from the engine's flat `SiteId`-indexed arrays is a direct
column write: filter occupied sites, write each field array as a column,
emit bonds as index pairs (apply the existing `i2 < i` once-per-edge rule),
mark periodic wrap edges with a `seam: bool` edge column instead of dropping
them. Keep the real `state` codes as an `i32` column — PGIF exists so the
chemistry is no longer flattened into MSI element strings.
