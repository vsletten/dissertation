//! petra-io — trajectory export for visualization (design doc P3).
//!
//! Two artifacts per run, consumed by graph-viz's trajectory mode:
//!
//! - a **PGIF v0 snapshot** (`pgif.rs`): every lattice site as a node
//!   (vacant included — occupancy changes during playback), bonds as
//!   edges, with `state`/`type`/`kind`/`frozen` columns and the
//!   state→display-type map in `meta.petra`;
//! - a **JSONL event log** (`events.rs`): one header line, then one line
//!   per fired event carrying `(site, old, new)` state deltas — enough to
//!   replay forward *and backward*.
//!
//! Format spec: `graph-viz/docs/PGIF.md` (trajectory sidecar section).

pub mod events;
pub mod pgif;

pub use events::EventLogWriter;
pub use pgif::snapshot_json;
