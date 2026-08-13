//! PGIF v0 JSON snapshot of an engine's lattice (spec:
//! graph-viz/docs/PGIF.md). All sites are emitted, vacant included, so a
//! trajectory player can flip occupancy without adding/removing nodes.

use serde_json::{json, Map, Value};

use petra_core::Engine;
use petra_deck::CompiledDeck;

/// Display type per state: the occupant species name, or `"vacant"`.
/// graph-viz's style table colors chemistry names (Al/Si/O/...) with the
/// element palette and renders `"vacant"` at radius 0.
fn display_types(deck: &CompiledDeck) -> Vec<String> {
    deck.state_occupants
        .iter()
        .map(|o| o.clone().unwrap_or_else(|| "vacant".to_string()))
        .collect()
}

/// Serialize the engine's current lattice as a PGIF v0 JSON document.
pub fn snapshot_json(deck: &CompiledDeck, engine: &Engine) -> String {
    let lat = &engine.lattice;
    let n = lat.len();

    // Positions through the cell matrix.
    let mut xs = Vec::with_capacity(n);
    let mut ys = Vec::with_capacity(n);
    let mut zs = Vec::with_capacity(n);
    for s in 0..n {
        let (cell, t) = lat.coords(s);
        let p = deck.unit_cell.cell.to_cartesian(
            deck.unit_cell.sites[t].frac,
            [cell[0] as i32, cell[1] as i32, cell[2] as i32],
        );
        xs.push(round6(p[0]));
        ys.push(round6(p[1]));
        zs.push(round6(p[2]));
    }

    // Per-node categorical columns.
    let state_types = display_types(deck);
    let mut type_dict: Vec<String> = Vec::new();
    let mut type_of_state = Vec::with_capacity(deck.n_states);
    for t in &state_types {
        let idx = match type_dict.iter().position(|d| d == t) {
            Some(i) => i,
            None => {
                type_dict.push(t.clone());
                type_dict.len() - 1
            }
        };
        type_of_state.push(idx as i64);
    }
    let states: Vec<i64> = lat.states.iter().map(|s| s.0 as i64).collect();
    let types: Vec<i64> = states.iter().map(|&s| type_of_state[s as usize]).collect();
    let kinds: Vec<i64> = lat
        .template_index
        .iter()
        .map(|&t| deck.kinds_per_template[t as usize].0 as i64)
        .collect();
    let frozen: Vec<bool> = lat.frozen.clone();

    // Edges once per pair (adjacency stores both directions), with a
    // best-effort `seam` flag for bonds that wrap a periodic boundary
    // (cell distance > 1 along a periodic axis — viz-only heuristic; the
    // viewer hides seam edges rather than drawing box-length cylinders).
    let mut src = Vec::new();
    let mut dst = Vec::new();
    let mut seam = Vec::new();
    for i in 0..n {
        let (ci, _) = lat.coords(i);
        for &j in lat.neighbors(i) {
            let j = j as usize;
            if j <= i {
                continue;
            }
            let (cj, _) = lat.coords(j);
            let wraps = (0..3).any(|ax| {
                let d = (ci[ax] as i64 - cj[ax] as i64).abs();
                d > 1
            });
            src.push(i as i64);
            dst.push(j as i64);
            seam.push(wraps);
        }
    }

    let mut meta = Map::new();
    meta.insert("producer".into(), json!("petra"));
    meta.insert("kind".into(), json!("kmc-lattice"));
    meta.insert("directed".into(), json!(false));
    meta.insert(
        "petra".into(),
        json!({
            "deck": deck.name,
            "temperature": deck.temperature,
            "states": deck.state_names,
            "state_types": state_types,
            "step": engine.step_count,
            "time": engine.time,
        }),
    );

    let doc = json!({
        "pgif": 1,
        "meta": Value::Object(meta),
        "nodes": {
            "count": n,
            "columns": {
                "x": { "type": "f32", "data": xs },
                "y": { "type": "f32", "data": ys },
                "z": { "type": "f32", "data": zs },
                "type": { "type": "categorical", "dict": type_dict, "data": types },
                "state": { "type": "categorical", "dict": deck.state_names, "data": states },
                "kind": { "type": "categorical", "dict": deck.kind_names, "data": kinds },
                "frozen": { "type": "bool", "data": frozen },
            },
        },
        "edges": {
            "count": src.len(),
            "src": src,
            "dst": dst,
            "columns": {
                "seam": { "type": "bool", "data": seam },
            },
        },
    });
    doc.to_string()
}

fn round6(v: f64) -> f64 {
    (v * 1e6).round() / 1e6
}

#[cfg(test)]
mod tests {
    use super::*;

    const MINI: &str = r#"
[deck]
name = "mini"

[cell]
a = 2.0
b = 2.0
c = 2.0
alpha = 90.0
beta = 90.0
gamma = 90.0

[[cell.sites]]
kind = "X"
frac = [0.0, 0.0, 0.0]

[[cell.bonds]]
i = 0
j = 0
dcell = [1, 0, 0]

[[species]]
name = "M"

[[kinds]]
name = "X"
initial = "occ"

[[kinds.states]]
name = "occ"
occupant = "M"

[[kinds.states]]
name = "gone"
occupant = "vacant"

[lattice]
dims = [4, 1, 1]
boundary = ["periodic", "periodic", "periodic"]

[thermo]
temperature = 300.0

[[reactions]]
name = "leave"
center = { kind = "X", state = ["occ"] }
rate = { constant = 1.0 }

[[reactions.effects]]
target = "center"
set = "gone"

[simulation]
steps = 10
seed = 1
"#;

    #[test]
    fn snapshot_is_valid_pgif_with_expected_shape() {
        let parsed: petra_deck::DeckFile = toml::from_str(MINI).unwrap();
        let deck = petra_deck::compile(&parsed).unwrap();
        let engine = deck.build_engine(Some(1)).unwrap();
        let doc: serde_json::Value = serde_json::from_str(&snapshot_json(&deck, &engine)).unwrap();

        assert_eq!(doc["pgif"], 1);
        assert_eq!(doc["nodes"]["count"], 4);
        // 4-ring: 4 undirected edges, each emitted once; the 0-3 bond wraps.
        assert_eq!(doc["edges"]["count"], 4);
        let seam: Vec<bool> = serde_json::from_value(
            doc["edges"]["columns"]["seam"]["data"].clone(),
        )
        .unwrap();
        assert_eq!(seam.iter().filter(|&&s| s).count(), 1, "one wrap edge");

        let dict: Vec<String> =
            serde_json::from_value(doc["nodes"]["columns"]["state"]["dict"].clone()).unwrap();
        assert_eq!(dict, ["X.occ", "X.gone"]);
        let types: Vec<String> =
            serde_json::from_value(doc["meta"]["petra"]["state_types"].clone()).unwrap();
        assert_eq!(types, ["M", "vacant"]);
    }
}
