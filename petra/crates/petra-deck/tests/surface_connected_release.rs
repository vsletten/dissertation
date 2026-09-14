//! Behavioral gate for E4a2's surface-connected lateral release front.

fn fixture(edge_cell: Option<usize>) -> String {
    let connected_init = if let Some(cell_a) = edge_cell {
        format!(
            r#"
[[structure.init]]
name = "connect-left-edge"
center = {{ kind = "Interface", state = ["bonded"] }}
sites = [[{cell_a}, 0, 0, 1]]
set = "delaminated"
"#,
        )
    } else {
        String::new()
    };
    let isotope_cell = if edge_cell == Some(3) { 1 } else { 2 };
    format!(
        r#"
[deck]
name = "surface-connected-release"
schema = 2

[structure]
kind = "cell"
[[structure.species]]
name = "Ar40"

[[structure.kinds]]
name = "Gallery"
initial = "vacant"
[[structure.kinds.states]]
name = "vacant"
occupant = "vacant"
[[structure.kinds.states]]
name = "Ar40"
occupant = "Ar40"

[[structure.kinds]]
name = "Interface"
initial = "bonded"
[[structure.kinds.states]]
name = "bonded"
occupant = "vacant"
[[structure.kinds.states]]
name = "delaminated"
occupant = "vacant"

[[structure.kinds]]
name = "Surface_gate"
initial = "closed"
[[structure.kinds.states]]
name = "closed"
occupant = "vacant"
[[structure.kinds.states]]
name = "edge"
occupant = "vacant"
[[structure.kinds.states]]
name = "open"
occupant = "vacant"

[structure.cell]
a = 1.0
b = 1.0
c = 1.0
alpha = 90.0
beta = 90.0
gamma = 90.0
[[structure.cell.sites]]
kind = "Gallery"
frac = [0.0, 0.0, 0.0]
[[structure.cell.sites]]
kind = "Interface"
frac = [0.25, 0.0, 0.0]
[[structure.cell.sites]]
kind = "Surface_gate"
frac = [0.5, 0.0, 0.0]
[[structure.cell.bonds]]
i = 0
j = 2
dcell = [0, 0, 0]
label = "surface_gate"
[[structure.cell.bonds]]
i = 1
j = 2
dcell = [0, 0, 0]
label = "gate_interface"
[[structure.cell.bonds]]
i = 2
j = 2
dcell = [1, 0, 0]
label = "lateral_front"

[structure.lattice]
dims = [4, 1, 1]
boundary = ["open", "periodic", "periodic"]

[[structure.init]]
name = "lateral-edges"
center = {{ kind = "Surface_gate", state = ["closed"] }}
sites = [[0, 0, 0, 2], [3, 0, 0, 2]]
set = "edge"
[[structure.init]]
name = "interior-isotope"
center = {{ kind = "Gallery", state = ["vacant"] }}
sites = [[{isotope_cell}, 0, 0, 0]]
set = "Ar40"
[[structure.init]]
name = "isolated-interior-pocket"
center = {{ kind = "Interface", state = ["bonded"] }}
sites = [[1, 0, 0, 1], [2, 0, 0, 1]]
set = "delaminated"
{connected_init}

[dynamics.thermo]
temperature = 300.0
[[dynamics.rules]]
name = "open_lateral_edge_after_delamination"
center = {{ kind = "Surface_gate", state = ["edge"] }}
guards = [{{ kind = "Interface", label = "gate_interface", state = ["delaminated"], min = 1 }}]
rate = {{ constant = 1.0 }}
[[dynamics.rules.effects]]
target = "center"
set = "open"
[[dynamics.rules]]
name = "advance_surface_connected_front"
center = {{ kind = "Surface_gate", state = ["closed"] }}
guards = [
  {{ kind = "Interface", label = "gate_interface", state = ["delaminated"], min = 1 }},
  {{ kind = "Surface_gate", label = "lateral_front", state = ["open"], min = 1 }},
]
rate = {{ constant = 1.0 }}
[[dynamics.rules.effects]]
target = "center"
set = "open"
[[dynamics.rules]]
name = "release_Ar40_surface"
center = {{ kind = "Gallery", state = ["Ar40"] }}
guards = [{{ kind = "Surface_gate", label = "surface_gate", state = ["open"], min = 1 }}]
rate = {{ constant = 1.0 }}
[[dynamics.rules.effects]]
target = "center"
set = "vacant"

[execution]
strategy = "ctmc"
[execution.stop]
steps = 10
[execution.ensemble]
seed = 1998
n_replicas = 1
seed_policy = "increment"

[observables]
report_every = 1
"#
    )
}

fn compiled(edge_cell: Option<usize>) -> petra_deck::CompiledDeck {
    let parsed: petra_deck::DeckFile = toml::from_str(&fixture(edge_cell)).expect("fixture parses");
    petra_deck::compile(&parsed).expect("fixture compiles")
}

fn state_name<'a>(
    deck: &'a petra_deck::CompiledDeck,
    engine: &petra_core::Engine,
    cell_a: usize,
    template: usize,
) -> &'a str {
    let site = engine.lattice.index([cell_a, 0, 0], template);
    &deck.state_names[engine.lattice.states[site].0 as usize]
}

#[test]
fn isotope_release_requires_a_delaminated_path_connected_to_a_lateral_edge() {
    let isolated = compiled(None);
    let mut isolated_engine = isolated.build_engine(None).expect("isolated engine builds");
    assert!(
        matches!(isolated_engine.step(), Err(petra_core::Stop::NoEvents)),
        "an isolated interior delaminated pocket must stop specifically because no event exists"
    );
    assert_eq!(
        state_name(&isolated, &isolated_engine, 2, 2),
        "Surface_gate.closed"
    );
    assert_eq!(
        state_name(&isolated, &isolated_engine, 2, 0),
        "Gallery.Ar40"
    );

    for (edge_cell, first_inner, second_inner, isotope_cell) in [(0, 1, 2, 2), (3, 2, 1, 1)] {
        let connected = compiled(Some(edge_cell));
        let mut connected_engine = connected
            .build_engine(None)
            .expect("connected engine builds");
        let mut fired = Vec::new();
        while let Ok(event) = connected_engine.step() {
            fired.push(connected.reactions[event.reaction as usize].name.as_str());
            assert!(
                fired.len() <= 4,
                "fixture must terminate after the four intended events"
            );
        }
        assert_eq!(
            fired,
            [
                "open_lateral_edge_after_delamination",
                "advance_surface_connected_front",
                "advance_surface_connected_front",
                "release_Ar40_surface",
            ]
        );
        for cell_a in [edge_cell, first_inner, second_inner] {
            assert_eq!(
                state_name(&connected, &connected_engine, cell_a, 2),
                "Surface_gate.open"
            );
        }
        assert_eq!(
            state_name(&connected, &connected_engine, isotope_cell, 0),
            "Gallery.vacant"
        );
    }
}
