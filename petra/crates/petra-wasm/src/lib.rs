//! WASM bindings for live in-browser simulation (viz stage 2).
//!
//! The browser gets the same artifacts the CLI's `--viz` writes — a PGIF
//! snapshot and trajectory-header/event rows in the docs/PGIF.md §6a
//! format — but generated on demand, so the viewer's existing trajectory
//! player consumes a *growing* event stream: pause the actual simulation,
//! scrub backward through its history, resume, re-seed, or rebuild from
//! an edited deck, all client-side.
//!
//! Determinism note: petra is IEEE f64 + PCG64 throughout, so a given
//! (deck, seed) reproduces the same trajectory here as in the native CLI.

use wasm_bindgen::prelude::*;

use petra_core::{Engine, Stop};
use petra_deck::CompiledDeck;

#[wasm_bindgen]
pub struct WasmSim {
    deck: CompiledDeck,
    engine: Engine,
    stopped: Option<String>,
}

#[wasm_bindgen]
impl WasmSim {
    /// Compile a deck from TOML text and build the engine. `seed`
    /// overrides the deck's seed when >= 0.
    #[wasm_bindgen(constructor)]
    pub fn new(deck_toml: &str, seed: f64) -> Result<WasmSim, JsError> {
        let parsed: petra_deck::DeckFile =
            toml::from_str(deck_toml).map_err(|e| JsError::new(&format!("deck parse: {e}")))?;
        let deck = petra_deck::compile(&parsed).map_err(|e| JsError::new(&e.to_string()))?;
        let seed_override = if seed >= 0.0 { Some(seed as u64) } else { None };
        let engine = deck
            .build_engine(seed_override)
            .map_err(|e| JsError::new(&e.to_string()))?;
        Ok(WasmSim {
            deck,
            engine,
            stopped: None,
        })
    }

    /// PGIF v0 snapshot of the CURRENT lattice state.
    pub fn snapshot_json(&self) -> String {
        petra_io::snapshot_json(&self.deck, &self.engine)
    }

    /// The trajectory header line (docs/PGIF.md §6a) for this sim.
    pub fn header_json(&self) -> String {
        let mut buf = Vec::new();
        // The header writer needs a seed for the record; report the deck's.
        petra_io::EventLogWriter::new(&mut buf, &self.deck, self.deck.seed, self.engine.lattice.len())
            .expect("writing to a Vec cannot fail");
        String::from_utf8(buf).expect("header is UTF-8").trim_end().to_string()
    }

    /// Run up to `n` events, returning them as a JSON array of §6a rows:
    /// `[[step, time, rxn, [[site, old, new], ...]], ...]`. Returns fewer
    /// than `n` (possibly zero) rows when the simulation stops; see
    /// `stop_reason()`.
    pub fn step_batch(&mut self, n: u32) -> String {
        let mut out = String::from("[");
        let mut first = true;
        for _ in 0..n {
            if self.stopped.is_some() {
                break;
            }
            match self.engine.step() {
                Ok(fired) => {
                    if !first {
                        out.push(',');
                    }
                    first = false;
                    out.push_str(&format!(
                        "[{},{:.9e},{}",
                        fired.step, fired.time, fired.reaction
                    ));
                    out.push_str(",[");
                    for (i, (site, old, new)) in self.engine.last_changes().iter().enumerate() {
                        if i > 0 {
                            out.push(',');
                        }
                        out.push_str(&format!("[{site},{},{}]", old.0, new.0));
                    }
                    out.push_str("]]");
                }
                Err(stop) => {
                    self.stopped = Some(match stop {
                        Stop::NoEvents => "no events possible — simulation complete".into(),
                        other => other.to_string(),
                    });
                    break;
                }
            }
        }
        out.push(']');
        out
    }

    /// Why the simulation stopped, or null while it can still advance.
    pub fn stop_reason(&self) -> Option<String> {
        self.stopped.clone()
    }

    pub fn step_count(&self) -> f64 {
        self.engine.step_count as f64
    }

    pub fn time(&self) -> f64 {
        self.engine.time
    }

    pub fn n_sites(&self) -> u32 {
        self.engine.lattice.len() as u32
    }
}
