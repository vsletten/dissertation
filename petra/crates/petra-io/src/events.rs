//! The JSONL event log — the trajectory sidecar to a PGIF snapshot.
//!
//! Line 1 is a header object; every following line is one fired event:
//!
//! ```text
//! {"petra_traj":1,"deck":"kaolinite","seed":42,"n_sites":1560,
//!  "states":[...],"state_types":[...],"reactions":[...]}
//! [step, time, rxn, [[site, old, new], ...]]
//! ```
//!
//! `old` is recorded so a player can scrub *backward* by applying deltas
//! inverted. State values are dense ids — indices into the header's
//! `states` array, which matches the snapshot's `state` column dict.

use std::io::Write;

use petra_core::Engine;
use petra_deck::CompiledDeck;

pub struct EventLogWriter<W: Write> {
    out: W,
    events: u64,
}

impl<W: Write> EventLogWriter<W> {
    /// Write the header line and return the writer.
    pub fn new(mut out: W, deck: &CompiledDeck, seed: u64, n_sites: usize) -> std::io::Result<Self> {
        let state_types: Vec<String> = deck
            .state_occupants
            .iter()
            .map(|o| o.clone().unwrap_or_else(|| "vacant".to_string()))
            .collect();
        let header = serde_json::json!({
            "petra_traj": 1,
            "deck": deck.name,
            "seed": seed,
            "n_sites": n_sites,
            "states": deck.state_names,
            "state_types": state_types,
            "reactions": deck.reactions.iter().map(|r| r.name.clone()).collect::<Vec<_>>(),
        });
        writeln!(out, "{header}")?;
        Ok(EventLogWriter { out, events: 0 })
    }

    /// Append the event the engine just applied (call after each
    /// successful `step`, before the next one).
    pub fn record(
        &mut self,
        fired: &petra_core::Fired,
        engine: &Engine,
    ) -> std::io::Result<()> {
        // Hand-formatted compact row — [step, time, rxn, [[site,old,new],..]]
        // — to keep 10^5-event logs cheap to write and parse.
        let mut line = format!("[{},{:.9e},{}", fired.step, fired.time, fired.reaction);
        line.push_str(",[");
        for (i, (site, old, new)) in engine.last_changes().iter().enumerate() {
            if i > 0 {
                line.push(',');
            }
            line.push_str(&format!("[{site},{},{}]", old.0, new.0));
        }
        line.push_str("]]");
        writeln!(self.out, "{line}")?;
        self.events += 1;
        Ok(())
    }

    pub fn events_written(&self) -> u64 {
        self.events
    }

    pub fn into_inner(self) -> W {
        self.out
    }
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
dims = [3, 1, 1]
boundary = ["open", "periodic", "periodic"]

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

    /// Run the mini deck to exhaustion, log it, and verify the log replays
    /// forward to the final state and backward to the initial state.
    #[test]
    fn log_replays_forward_and_backward()  {
        let parsed: petra_deck::DeckFile = toml::from_str(MINI).unwrap();
        let deck = petra_deck::compile(&parsed).unwrap();
        let mut engine = deck.build_engine(Some(3)).unwrap();
        let initial: Vec<u16> = engine.lattice.states.iter().map(|s| s.0).collect();

        let mut log = EventLogWriter::new(Vec::new(), &deck, 3, engine.lattice.len()).unwrap();
        while let Ok(fired) = engine.step() {
            log.record(&fired, &engine).unwrap();
        }
        let final_states: Vec<u16> = engine.lattice.states.iter().map(|s| s.0).collect();
        assert_eq!(log.events_written(), 3, "3 sites dissolve in 3 events");

        let text = String::from_utf8(log.into_inner()).unwrap();
        let mut lines = text.lines();
        let header: serde_json::Value = serde_json::from_str(lines.next().unwrap()).unwrap();
        assert_eq!(header["petra_traj"], 1);
        assert_eq!(header["n_sites"], 3);
        assert_eq!(header["reactions"][0], "leave");

        // Replay forward over the initial state.
        let mut replay = initial.clone();
        let mut events = Vec::new();
        for line in lines {
            let row: (u64, f64, u16, Vec<(usize, u16, u16)>) =
                serde_json::from_str(line).unwrap();
            for &(site, old, new) in &row.3 {
                assert_eq!(replay[site], old, "forward replay consistency");
                replay[site] = new;
            }
            events.push(row);
        }
        assert_eq!(replay, final_states);

        // And backward to the start.
        for row in events.iter().rev() {
            for &(site, old, new) in row.3.iter().rev() {
                assert_eq!(replay[site], new, "backward replay consistency");
                replay[site] = old;
            }
        }
        assert_eq!(replay, initial);
    }
}
