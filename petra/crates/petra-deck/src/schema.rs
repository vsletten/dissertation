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
    /// Ordered build-time passes applied after the uniform per-kind fill
    /// and before dynamics: surface termination, region clearing, defect
    /// seeding. Each pass sweeps all sites in index order with writes
    /// immediately visible (the legacy TerminateSurface convention).
    #[serde(default)]
    pub init: Vec<InitPassSpec>,
    #[serde(default)]
    pub reactions: Vec<ReactionSpec>,
    pub simulation: SimSpec,
}

/// One build-time pass. The operation applies to the *center* site — once,
/// or once per neighbor matching `foreach` (in adjacency order), so
/// "step the map per missing cation" and "increment per terminal OH" are
/// both expressible.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InitPassSpec {
    pub name: String,
    /// Which sites this pass rewrites: kind (optional) + state set.
    pub center: CenterInitSpec,
    /// Restrict to a cell-coordinate slab: axis 0/1/2 and an inclusive
    /// range (either bound optional).
    #[serde(default)]
    pub region: Option<RegionSpec>,
    /// Additional guards on the center's neighborhood.
    #[serde(default)]
    pub guards: Vec<SelectorSpec>,
    /// Apply the op once per neighbor matching this selector.
    #[serde(default)]
    pub foreach: Option<SelectorSpec>,
    #[serde(default)]
    pub set: Option<String>,
    #[serde(default)]
    pub shift: Option<i32>,
    #[serde(default)]
    pub map: Option<BTreeMap<String, String>>,
    /// `map` miss policy; init defaults to `"skip"` (the termination maps
    /// leave unlisted states alone).
    #[serde(default)]
    pub missing: Option<String>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CenterInitSpec {
    #[serde(default)]
    pub kind: Option<String>,
    pub state: Vec<String>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RegionSpec {
    /// 0 = a, 1 = b, 2 = c.
    pub axis: u8,
    #[serde(default)]
    pub min: Option<usize>,
    #[serde(default)]
    pub max: Option<usize>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Meta {
    pub name: String,
    #[serde(default)]
    pub comment: Option<String>,
    /// Energy unit for every energy-valued field in this deck (activation
    /// energies, ΔEa modifiers, chemical potentials, ΔH‡, and ΔS‡ per K):
    /// `"kcal/mol"` (default), `"kJ/mol"`, or `"eV"`. Temperatures are
    /// always Kelvin; prefactors always 1/s. Converted once at compile
    /// time — the runtime stays in kcal/mol internally.
    #[serde(default)]
    pub units: Option<String>,
}

/// Cell geometry: either conventional parameters (`a`..`gamma`, angles in
/// degrees) or an explicit fractional→Cartesian `matrix` (columns are the
/// cell vectors) for nonstandard conventions. Exactly one form.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CellSpec {
    #[serde(default)]
    pub a: Option<f64>,
    #[serde(default)]
    pub b: Option<f64>,
    #[serde(default)]
    pub c: Option<f64>,
    #[serde(default)]
    pub alpha: Option<f64>,
    #[serde(default)]
    pub beta: Option<f64>,
    #[serde(default)]
    pub gamma: Option<f64>,
    #[serde(default)]
    pub matrix: Option<[[f64; 3]; 3]>,
    pub sites: Vec<SiteSpec>,
    #[serde(default)]
    pub bonds: Vec<BondSpec>,
}

/// A site position: fractional (`frac`) or Cartesian (`cart`, converted
/// through the cell matrix at compile time). Exactly one.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SiteSpec {
    pub kind: String,
    #[serde(default)]
    pub frac: Option<[f64; 3]>,
    #[serde(default)]
    pub cart: Option<[f64; 3]>,
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
    /// State names, `"Kind.state"` qualified names, `"@alias"` refs, or the
    /// wildcard `"*"` (every state in scope — with `kind` set this makes
    /// the selector a degree/coordination counter).
    pub state: Vec<String>,
    /// Restrict to frozen boundary sites (`true`) or live sites (`false`) —
    /// e.g. "occupied AND not part of the frozen wall".
    #[serde(default)]
    pub frozen: Option<bool>,
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

/// Exactly one of `per_match` / `by_count` / `when` must be present.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ModifierSpec {
    pub select: SelectorSpec,
    #[serde(default)]
    pub per_match: Option<PerMatchSpec>,
    #[serde(default)]
    pub by_count: Option<ByCountSpec>,
    #[serde(default)]
    pub when: Option<WhenSpec>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PerMatchSpec {
    /// ΔEa added per matching neighbor (linear convenience form of
    /// `by_count`), in deck units.
    pub dea: f64,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ByCountSpec {
    /// Tabulated ΔEa by match count: `dea[n]` applies when `n` neighbors
    /// match; the last entry extends to all higher counts. The general
    /// nonlinear form — barriers are rarely linear in coordination.
    pub dea: Vec<f64>,
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

/// One state rewrite. Exactly one operation: `set` (fixed state), `shift`
/// (±n along the kind's declared state ladder — the protonation-counter
/// pattern), or `map` (per-state transition table — the adsorption/
/// desorption oxygen-shell pattern).
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectSpec {
    /// `"center"`, `"neighbor"` (first match; matching nothing at apply
    /// time is an error), or `"neighbors"` (all matches; zero is legal).
    pub target: String,
    /// Required for neighbor targets; must name a `kind` so state names in
    /// `set`/`map` resolve unambiguously (v0 restriction).
    #[serde(default)]
    pub select: Option<SelectorSpec>,
    #[serde(default)]
    pub set: Option<String>,
    #[serde(default)]
    pub shift: Option<i32>,
    #[serde(default)]
    pub map: Option<BTreeMap<String, String>>,
    /// `map` policy for a matched site whose state has no entry:
    /// `"error"` (default — the legacy adsorb fatal) or `"skip"` (leave
    /// unchanged — the legacy desorb silence).
    #[serde(default)]
    pub missing: Option<String>,
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
