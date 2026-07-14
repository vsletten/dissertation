# 05 — PGIF Draft (Property-Graph Interchange Format)

A simple, general interchange format to replace **CERIUS 2 MSI at the source**: the
Rust KMC port emits PGIF, and the viewer app (`dissertation/main/viewer`, today an MSI
text parser) consumes it generically. Design goals, in priority order:

1. **General property graph** — nodes and edges with arbitrary typed properties. Not
   crystal-specific; a molecule, a KMC lattice, a social graph, or a knowledge graph
   all serialize the same way. (This is also why `graph-viz` wants it: one input format,
   many producers.)
2. **Large-graph friendly** — columnar / typed-array layout so millions of nodes load
   without per-node parsing, and so a browser can push columns straight into typed
   arrays / GPU buffers.
3. **Simple** — a producer (Rust) and a consumer (TypeScript) should each be ~a
   screen of code. Resist ontology/schema-registry ambition.

This is a **draft for co-design**, not a frozen spec. Open questions are called out in §7.

---

## 1. Why replace MSI

MSI (spec A9.1, `viewer/msi-schema.md`) is a 1990s parenthesized object-tree text
format: every atom is `(N Atom (A C ACL "13 Al") (A D XYZ (...)) ...)`, parsed line by
line with regexes (`viewer/src/utils/cerius2Parser.ts`). Problems for our use:
- **Row-oriented text** → O(atoms) string parsing; poor for large graphs.
- **Crystal/molecule-specific vocabulary** (`ACL`, `Atom`, `Bond`, element strings) —
  not a general graph, and the viewer hard-codes element→color/size from it.
- **Atom and bond IDs share one integer space** and bonds reference atom object IDs
  (spec A9.1) — a quirk the parser must special-case.
- **Lossy for our model**: the KMC site's real identity is its `state` code (401,
  502, …) and its class; MSI flattens that to an element symbol plus a text label, and
  even abuses "9 F" to render Al-OH-Al oxygens (spec B5/A9.1). The viewer can't recover
  the chemistry.

PGIF keeps what the viewer needs (positions, a type per node, edges, per-node
properties for coloring) while being a generic property graph and cheap to load.

---

## 2. Format at a glance

PGIF is **columnar**: a graph is a set of named columns for nodes and a set for edges,
plus small metadata. Two serializations, same logical model:

- **`.pgif.json`** — human-readable, diff-friendly, the default for small/medium
  graphs and for getting the viewer working. Columns are JSON arrays.
- **`.pgif.msgpack`** (or a tiny binary container) — the large-graph path: the same
  structure, but numeric columns are raw typed-array buffers (Float32/Int32) so the
  consumer maps them directly into `Float32Array` with zero per-element parsing.

Start with JSON (M9 of doc 03); add the binary encoding only when a graph is big enough
to need it. The logical schema is identical, so the viewer's loader branches only at the
byte-decoding layer.

### 2.1 Logical model
```
Graph = {
  meta:  { … free-form provenance … },
  nodes: {
    count: N,
    ids:   [len N]              // stable node ids (ints; index i has id ids[i])
    columns: { <name>: Column } // one entry per node property, each length N
  },
  edges: {
    count: M,
    src:   [len M],             // node id (or node index — see §7 Q1)
    dst:   [len M],
    columns: { <name>: Column } // optional per-edge properties, each length M
  }
}
Column = { type: "f32"|"f64"|"i32"|"u32"|"i8"|"bool"|"str"|"categorical",
           data: [...] | <buffer ref>,
           // for categorical: dict: [strings], data: [i32 indices]
}
```
Everything is **structure-of-arrays**: to read node *i*'s properties you index each
column at *i*. This is the columnar/typed-array-friendly core.

---

## 3. Concrete example — the KMC lattice as PGIF (JSON form)

A tiny lattice, showing how the kaolinite model maps on generically:

```jsonc
{
  "pgif": 1,
  "meta": {
    "producer": "mckaol-rs 0.1",
    "kind": "kmc-lattice",              // free-form; consumers may ignore
    "step": 200000,
    "time": 0.0622065,
    "cell": { "a": 5.140, "b": 8.930, "c": 7.370,
              "alpha": -0.03141, "beta": -0.25038, "gamma": 0.0 },
    "periodic": ["a"]                   // which axes wrap (for edge-seam handling)
  },
  "nodes": {
    "count": 3,
    "ids": [2, 3, 6],
    "columns": {
      "x":     { "type": "f32", "data": [1.17181, 1.30217, 0.24775] },
      "y":     { "type": "f32", "data": [4.38134, 7.28290, 3.10361] },
      "z":     { "type": "f32", "data": [2.84548, 2.58231, 0.225645] },
      "type":  { "type": "categorical", "dict": ["Al","Si","O","OH"],
                 "data": [0, 0, 1] },   // generic node type → drives color/size
      "state": { "type": "i32", "data": [101, 101, 201] },  // model-specific detail
      "label": { "type": "str", "data": ["Al0","Al1","Si4"] }
    }
  },
  "edges": {
    "count": 2,
    "src": [2, 2],
    "dst": [6, 3],
    "columns": {
      "seam": { "type": "bool", "data": [false, false] }  // true = periodic wrap edge
    }
  }
}
```

Notes on the mapping:
- **`type` (categorical)** is the generic coloring key — the viewer maps dict entries to
  colors/radii via a small style table, replacing today's hard-coded element→color map.
  For a non-crystal graph the dict is whatever the producer wants ("Person","Org",…).
- **`state`** (and any other model fields) ride along as ordinary columns — the viewer
  can ignore them or use them for tooltips/filters. Chemistry is **no longer lost**
  (fixes MSI's lossiness): the real state code is preserved, not flattened to "9 F".
- **`seam`** replaces the MSI writer's special periodic-bond suppression (spec A9.1): the
  producer *emits* wrap edges but marks them, so the consumer decides whether to draw
  them, rather than the producer silently dropping them. More faithful and more flexible.
- **positions are optional** — a graph with no `x/y/z` columns is a pure topological
  property graph; the viewer can lay it out itself. This is what makes PGIF general
  beyond crystals.

---

## 4. The binary path (large graphs)

Same logical schema; numeric columns become buffers instead of JSON arrays. Minimal
container:

```
[ 8 bytes  magic "PGIF\0\0\0\1" ]
[ u32      header_len ]
[ header   JSON (UTF-8): the structure above, but every column's "data"
           is replaced by { "type": ..., "offset": B, "len": N } pointing
           into the blob region; str/categorical dicts stay inline in JSON ]
[ blob     concatenated little-endian typed-array buffers, 8-byte aligned ]
```
The consumer reads the JSON header, then for each numeric column does
`new Float32Array(blob, offset, len)` — **zero-copy, zero per-element parsing**. This is
the property that lets millions of nodes load in a browser: positions and scalar
properties become GPU-ready typed arrays directly. (MessagePack is a fine off-the-shelf
alternative to a bespoke container; the point is columns-as-buffers, not the exact
framing — §7 Q2.)

Rule of thumb: JSON form under ~10⁴–10⁵ nodes; binary form above. The viewer supports
both behind one `loadPgif(bytes)` that sniffs the magic.

---

## 5. Why this shape is right for both producer and consumer

- **Rust producer**: the KMC engine already stores state as flat arrays indexed by
  `SiteId` (doc 03 §3). Emitting PGIF is "write each field array as a column" — nearly a
  direct memcpy for the binary form, and `serde` for JSON. No tree building, no shared
  atom/bond ID juggling (edges just reference node ids). One occupied-site filter pass
  (as `writeMSI` already does) selects which nodes to emit.
- **TS consumer**: columnar → the renderer wants `Float32Array`s of positions for
  instanced meshes anyway (the graph-viz recon flagged per-atom components + no
  instancing as the current bottleneck). PGIF hands the renderer exactly that. The
  categorical `type` column drives a style lookup; extra columns become hover/tooltip
  data. The generic loader is ~a screen of code and has no crystal knowledge.
- **Generality**: nothing above says "atom". Nodes, edges, typed columns, optional
  positions. A KMC lattice, a molecule, or an arbitrary property graph all fit.

---

## 6. Migration plan (low-risk, incremental)

1. **Producer emits PGIF alongside MSI** (doc 03 M9). Nothing breaks; MSI still works.
2. **Viewer gains `loadPgif`** next to the existing MSI parser; a file-type sniff picks
   the loader. Render path is shared (both produce `{nodes, edges}` in memory).
3. **Viewer migrates its render layer to columns/instancing** driven by PGIF (this is a
   graph-viz project deliverable, not a KMC one — coordinate there).
4. **Deprecate MSI output** once the viewer defaults to PGIF and old MSI files are no
   longer needed. Keep an MSI reader in the viewer for legacy files if convenient.

PGIF is designed once here, consumed by two projects (kmc-port produces, graph-viz
consumes) — the shared-format goal from the mission-control registry.

---

## 7. Open questions for co-design

1. **Edges reference node *ids* or node *indices*?** Ids are stable and human-readable;
   indices are cheaper for the consumer (direct array lookup) and natural for the binary
   form. Leaning **indices in the binary form, ids in JSON** — or standardize on
   `ids[i]==i` (dense) and make them equivalent. Decide with the viewer's needs.
2. **Binary container: bespoke vs MessagePack vs a slice of Arrow/Parquet?** Bespoke is
   simplest and dependency-free; MessagePack is a one-line dep both sides; Arrow is the
   "real" columnar standard but heavy for a nostalgia project. Recommend **bespoke
   magic+JSON-header+blob** (§4) to start — trivially replaceable later.
3. **Directed vs undirected edges?** Crystal bonds are undirected; KG edges are
   directed. Add `meta.directed: bool` (default false) and emit each undirected edge
   once, mirroring `writeMSI`'s `i2 < i` rule.
4. **Trajectories / time series** (many KMC steps): a sequence of PGIF frames? A single
   file with per-frame state columns (topology fixed, only `state` changes — true for
   this model since topology is static, spec B6)? For KMC specifically, a "static
   topology + per-frame state column" variant would be very compact — worth a follow-up
   once the single-frame format is settled.
5. **How much `meta` to standardize?** Keep it free-form for now (`producer`, `kind`,
   plus whatever the producer wants); promote fields to the spec only once two producers
   actually share them.

---

## 8. Minimal spec summary (v0, for implementers)

- A PGIF document has `pgif: 1`, `meta` (object), `nodes`, `edges`.
- `nodes`: `count`, `ids` (length count), `columns` (map name→Column, each length count).
- `edges`: `count`, `src`, `dst` (each length count), optional `columns` (length count).
- `Column`: `type` ∈ {f32,f64,i32,u32,i8,bool,str,categorical}; `data` array (JSON) or
  `{offset,len}` (binary). `categorical` adds `dict: [str]`, `data` = indices.
- Positions are the columns `x,y,z` (optional). Node type/coloring is the categorical
  column `type` (optional). Everything else is free-form columns.
- JSON serialization for small graphs; magic+header+blob binary for large ones; one
  consumer entry point sniffs which.

That is the whole format. Simple enough to implement in an afternoon on each side,
general enough to outlive the crystal model, and columnar enough to scale.
