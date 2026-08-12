//! The deck file schema, one-to-one with the TOML the user writes.
//! All names are strings here; `compile` resolves them to dense ids.
//! Annotated example: `petra/examples/kossel.toml`; design doc §3.4.

use std::collections::BTreeMap;

use serde::Deserialize;

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DeckFile {
    pub deck: Meta,
    pub cell: CellSpec,
    #[serde(default)]
    pub species: Vec<SpeciesSpec>,
    pub kinds: Vec<KindSpec>,
    /// Named state-set aliases, referenced as `"@name"` in selectors.
    #[serde(default)]
    pub aliases: BTreeMap<String, Vec<String>>,
    pub lattice: LatticeSpec,
    pub thermo: ThermoSpec,
    #[serde(default)]
    pub reactions: Vec<ReactionSpec>,
    pub simulation: SimSpec,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Meta {
    pub name: String,
    #[serde(default)]
    pub comment: Option<String>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CellSpec {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    pub alpha: f64,
    pub beta: f64,
    pub gamma: f64,
    pub sites: Vec<SiteSpec>,
    #[serde(default)]
    pub bonds: Vec<BondSpec>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SiteSpec {
    pub kind: String,
    pub frac: [f64; 3],
}

/// One declared bond; the compiler expands it onto both endpoints.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BondSpec {
    pub i: usize,
    pub j: usize,
    pub dcell: [i32; 3],
    #[serde(default)]
    pub label: Option<String>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SpeciesSpec {
    pub name: String,
    /// Bookkeeping only in v1 (design doc §6): stored so future
    /// charge-balance linting has the data without a schema break.
    #[serde(default)]
    pub charge: Option<f64>,
    #[serde(default)]
    pub mass: Option<f64>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct KindSpec {
    pub name: String,
    /// State every site of this kind starts in (perfect-crystal fill;
    /// defect fills layer on top, design doc §6).
    pub initial: String,
    pub states: Vec<StateSpec>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StateSpec {
    pub name: String,
    /// A species name, or `"vacant"`.
    pub occupant: String,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LatticeSpec {
    pub dims: [usize; 3],
    /// Per axis: `"periodic"`, `"open"`, or `"fixed"`.
    pub boundary: [String; 3],
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ThermoSpec {
    /// Kelvin.
    pub temperature: f64,
    /// Chemical potential offsets Δμ per species, kcal/mol.
    #[serde(default)]
    pub mu: BTreeMap<String, f64>,
    /// Solution activities per species (default 1).
    #[serde(default)]
    pub activity: BTreeMap<String, f64>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReactionSpec {
    pub name: String,
    pub center: CenterSpec,
    #[serde(default)]
    pub guards: Vec<SelectorSpec>,
    pub rate: RateSpec,
    /// Species drawn from solution: each contributes activity and Δμ
    /// factors to the forward rate (design doc §4).
    #[serde(default)]
    pub consumes: Vec<String>,
    /// Species released to solution (bookkeeping in v1; reverse rates come
    /// from explicit reverse reactions until the P1 auto-reverse lands).
    #[serde(default)]
    pub produces: Vec<String>,
    #[serde(default)]
    pub modifiers: Vec<ModifierSpec>,
    /// Deterministic outcome (exactly one of `effects` / `branches`).
    #[serde(default)]
    pub effects: Vec<EffectSpec>,
    /// Weighted alternative outcomes.
    #[serde(default)]
    pub branches: Vec<BranchSpec>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CenterSpec {
    pub kind: String,
    pub state: Vec<String>,
}

/// A neighbor selector; doubles as a guard when `min`/`max` are read.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SelectorSpec {
    /// Graph distance, 1 (default) or 2.
    #[serde(default)]
    pub distance: Option<u8>,
    #[serde(default)]
    pub kind: Option<String>,
    /// Bond label filter (distance-1 selectors only).
    #[serde(default)]
    pub label: Option<String>,
    /// State names, `"Kind.state"` qualified names, or `"@alias"` refs.
    pub state: Vec<String>,
    /// Guard bounds: default min=1, max=unbounded.
    #[serde(default)]
    pub min: Option<u32>,
    #[serde(default)]
    pub max: Option<u32>,
}

/// Exactly one variant must be present.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RateSpec {
    #[serde(default)]
    pub constant: Option<f64>,
    #[serde(default)]
    pub arrhenius: Option<ArrheniusSpec>,
    #[serde(default)]
    pub eyring: Option<EyringSpec>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ArrheniusSpec {
    pub prefactor: f64,
    /// kcal/mol.
    pub ea: f64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EyringSpec {
    /// ΔH‡, kcal/mol.
    pub dh: f64,
    /// ΔS‡, kcal·mol⁻¹·K⁻¹.
    pub ds: f64,
}

/// Exactly one of `per_match` / `when` must be present.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ModifierSpec {
    pub select: SelectorSpec,
    #[serde(default)]
    pub per_match: Option<PerMatchSpec>,
    #[serde(default)]
    pub when: Option<WhenSpec>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PerMatchSpec {
    /// ΔEa added per matching neighbor, kcal/mol.
    pub dea: f64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WhenSpec {
    #[serde(default)]
    pub min: Option<u32>,
    #[serde(default)]
    pub max: Option<u32>,
    #[serde(default)]
    pub dea: Option<f64>,
    #[serde(default)]
    pub factor: Option<f64>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectSpec {
    /// `"center"` or `"neighbor"`.
    pub target: String,
    /// Required when `target = "neighbor"`; must name a `kind` so `set`
    /// resolves unambiguously (v0 restriction).
    #[serde(default)]
    pub select: Option<SelectorSpec>,
    /// New state name for the target site.
    pub set: String,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BranchSpec {
    pub weight: f64,
    pub effects: Vec<EffectSpec>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SimSpec {
    pub steps: u64,
    pub seed: u64,
    #[serde(default)]
    pub report_every: Option<u64>,
}
