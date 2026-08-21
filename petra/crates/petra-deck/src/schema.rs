//! The deck file schema, one-to-one with the TOML the user writes.
//! All names are strings here; `compile` resolves them to dense ids.
//! Annotated example: `petra/examples/kossel.toml`; design doc §3.4.

use std::collections::BTreeMap;

use serde::Deserialize;

/// Canonical in-memory deck. Both v1 and v2 TOML deserialize to this shape;
/// v1 is upgraded by the custom deserializer before compilation.
#[derive(Debug)]
pub struct DeckFile {
    pub deck: Meta,
    pub structure_kind: StructureKind,
    pub grid: Option<GridSpec>,
    pub cell: CellSpec,
    pub species: Vec<SpeciesSpec>,
    pub kinds: Vec<KindSpec>,
    pub aliases: BTreeMap<String, Vec<String>>,
    pub lattice: LatticeSpec,
    pub defects: Vec<DefectSpec>,
    pub thermo: ThermoSpec,
    pub init: Vec<InitPassSpec>,
    pub reactions: Vec<ReactionSpec>,
    pub execution: ExecutionSpec,
    pub observables: ObservablesSpec,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StructureKind {
    Cell,
    Grid,
}

#[derive(Debug)]
pub struct GridSpec {
    pub family: String,
    pub neighborhood: Option<String>,
    pub dims: Vec<usize>,
    pub boundary: Vec<String>,
    pub default_kind: String,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExecutionSpec {
    pub strategy: String,
    #[serde(default)]
    pub ctmc: Option<CtmcSpec>,
    #[serde(default)]
    pub synchronous: Option<SynchronousSpec>,
    #[serde(default)]
    pub metropolis: Option<MetropolisSpec>,
    #[serde(default)]
    pub pca: Option<PcaSpec>,
    #[serde(default)]
    pub stop: StopSpec,
    #[serde(default)]
    pub ensemble: EnsembleSpec,
}

#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CtmcSpec {}

#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
/// Synchronous CA parameters. When multiple rules are enabled at one site,
/// the first rule in deck declaration order wins (RFC-001 §3).
pub struct SynchronousSpec {
    /// Must be `first_match`: synchronous updates are deterministic and fire
    /// at most one rule per site, prioritized by deck declaration order.
    pub conflict_resolution: String,
}

#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MetropolisSpec {
    #[serde(default)]
    pub temperature: Option<f64>,
}

#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PcaSpec {}

#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StopSpec {
    #[serde(default)]
    pub steps: Option<u64>,
    #[serde(default)]
    pub time: Option<f64>,
    #[serde(default)]
    pub predicate: Option<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SeedPolicy {
    Increment,
    Hash,
}

fn one_replica() -> u64 {
    1
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EnsembleSpec {
    /// Base seed. v1 `[simulation] seed` is moved here by the shim.
    #[serde(default)]
    pub seed: u64,
    #[serde(default = "one_replica")]
    pub n_replicas: u64,
    #[serde(default = "default_seed_policy")]
    pub seed_policy: SeedPolicy,
}

fn default_seed_policy() -> SeedPolicy {
    SeedPolicy::Increment
}

impl Default for EnsembleSpec {
    fn default() -> Self {
        Self {
            seed: 0,
            n_replicas: 1,
            seed_policy: SeedPolicy::Increment,
        }
    }
}

#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ObservablesSpec {
    #[serde(default)]
    pub report_every: u64,
    #[serde(default)]
    pub series: Vec<ObservableSpec>,
}

#[derive(Debug, Deserialize)]
pub struct ObservableSpec {
    pub kind: String,
    #[serde(flatten)]
    pub parameters: BTreeMap<String, toml::Value>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct DeckV1 {
    deck: Meta,
    cell: CellSpec,
    #[serde(default)]
    species: Vec<SpeciesSpec>,
    kinds: Vec<KindSpec>,
    #[serde(default)]
    aliases: BTreeMap<String, Vec<String>>,
    lattice: LatticeSpec,
    #[serde(default)]
    defects: Vec<DefectSpec>,
    thermo: ThermoSpec,
    #[serde(default)]
    init: Vec<InitPassSpec>,
    #[serde(default)]
    reactions: Vec<ReactionSpec>,
    simulation: SimSpec,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct DeckV2 {
    deck: Meta,
    structure: StructureV2,
    dynamics: DynamicsV2,
    execution: ExecutionSpec,
    #[serde(default)]
    observables: ObservablesSpec,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct StructureV2 {
    #[serde(default)]
    kind: Option<String>,
    #[serde(default)]
    cell: Option<CellSpec>,
    #[serde(default)]
    lattice: Option<LatticeSpec>,
    #[serde(default)]
    grid: Option<String>,
    #[serde(default)]
    neighborhood: Option<String>,
    #[serde(default)]
    dims: Option<Vec<usize>>,
    #[serde(default)]
    boundary: Option<Vec<String>>,
    #[serde(default)]
    default_kind: Option<String>,
    #[serde(default)]
    species: Vec<SpeciesSpec>,
    kinds: Vec<KindSpec>,
    #[serde(default)]
    init: Vec<InitPassSpec>,
    #[serde(default)]
    defects: Vec<DefectSpec>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct DynamicsV2 {
    thermo: ThermoSpec,
    #[serde(default)]
    aliases: BTreeMap<String, Vec<String>>,
    #[serde(default)]
    rules: Vec<ReactionSpec>,
}

impl DeckFile {
    fn from_v1(mut v1: DeckV1) -> Self {
        v1.deck.schema = Some(2);
        let execution = ExecutionSpec {
            strategy: "ctmc".to_string(),
            ctmc: None,
            synchronous: None,
            metropolis: None,
            pca: None,
            stop: StopSpec {
                steps: Some(v1.simulation.steps),
                time: None,
                predicate: None,
            },
            ensemble: EnsembleSpec {
                seed: v1.simulation.seed,
                n_replicas: 1,
                seed_policy: SeedPolicy::Increment,
            },
        };
        let observables = ObservablesSpec {
            report_every: v1.simulation.report_every.unwrap_or(0),
            series: Vec::new(),
        };
        Self {
            deck: v1.deck,
            structure_kind: StructureKind::Cell,
            grid: None,
            cell: v1.cell,
            species: v1.species,
            kinds: v1.kinds,
            aliases: v1.aliases,
            lattice: v1.lattice,
            defects: v1.defects,
            thermo: v1.thermo,
            init: v1.init,
            reactions: v1.reactions,
            execution,
            observables,
        }
    }

    fn from_v2(v2: DeckV2) -> Result<Self, String> {
        if v2.deck.schema != Some(2) {
            return Err("v2 deck must declare [deck] schema = 2".to_string());
        }
        validate_v2_surfaces(&v2)?;
        let StructureV2 {
            kind,
            cell,
            lattice,
            grid,
            neighborhood,
            dims,
            boundary,
            default_kind,
            species,
            kinds,
            init,
            defects,
        } = v2.structure;
        let (structure_kind, grid_spec, cell, lattice) = match kind.as_deref().unwrap_or("cell") {
            "cell" => {
                if grid.is_some()
                    || neighborhood.is_some()
                    || dims.is_some()
                    || boundary.is_some()
                    || default_kind.is_some()
                {
                    return Err("cell structure cannot contain grid-form fields".to_string());
                }
                (
                    StructureKind::Cell,
                    None,
                    cell.ok_or("[structure] kind = 'cell' requires [structure.cell]")?,
                    lattice.ok_or("[structure] kind = 'cell' requires [structure.lattice]")?,
                )
            }
            "grid" => {
                if cell.is_some() || lattice.is_some() {
                    return Err("grid structure cannot contain cell-form fields".to_string());
                }
                let dims_vec = dims.ok_or("[structure] kind = 'grid' requires dims")?;
                let boundary_vec = boundary.ok_or("[structure] kind = 'grid' requires boundary")?;
                let grid_spec = GridSpec {
                    family: grid.ok_or("[structure] kind = 'grid' requires grid")?,
                    neighborhood,
                    dims: dims_vec,
                    boundary: boundary_vec,
                    default_kind: default_kind
                        .ok_or("[structure] kind = 'grid' requires default_kind")?,
                };
                // Grid compilation is B3. Preserve the parsed grid verbatim
                // and carry inert placeholders that compile() never reaches.
                (
                    StructureKind::Grid,
                    Some(grid_spec),
                    CellSpec {
                        a: None,
                        b: None,
                        c: None,
                        alpha: None,
                        beta: None,
                        gamma: None,
                        matrix: Some([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
                        sites: Vec::new(),
                        bonds: Vec::new(),
                    },
                    LatticeSpec {
                        dims: [1, 1, 1],
                        boundary: [
                            "periodic".to_string(),
                            "periodic".to_string(),
                            "periodic".to_string(),
                        ],
                    },
                )
            }
            other => return Err(format!("unknown structure kind '{other}'")),
        };
        Ok(Self {
            deck: v2.deck,
            structure_kind,
            grid: grid_spec,
            cell,
            species,
            kinds,
            aliases: v2.dynamics.aliases,
            lattice,
            defects,
            thermo: v2.dynamics.thermo,
            init,
            reactions: v2.dynamics.rules,
            execution: v2.execution,
            observables: v2.observables,
        })
    }
}

fn validate_v2_surfaces(deck: &DeckV2) -> Result<(), String> {
    validate_execution(&deck.execution)?;
    validate_structure(&deck.structure)?;

    for rule in &deck.dynamics.rules {
        let ctx = format!("rule '{}'", rule.name);
        if let Some(truth) = rule.truth.as_deref() {
            if truth != "and" {
                return Err(format!("{ctx}: truth must be 'and'"));
            }
        }
        let effects: Vec<&EffectSpec> = rule
            .effects
            .iter()
            .chain(rule.branches.iter().flat_map(|branch| &branch.effects))
            .collect();
        if deck.execution.strategy == "synchronous"
            && (rule.branches.len() > 1
                || effects
                    .iter()
                    .any(|effect| effect.select_mode.as_deref() == Some("one")))
        {
            return Err(format!(
                "{ctx}: synchronous rules are draw-free and cannot use weighted branches or select_mode = 'one'"
            ));
        }
        let source = rule.center.is_none();
        if source && effects.iter().any(|effect| effect.target != "source") {
            return Err(format!(
                "{ctx}: a rule without center may contain only source effects"
            ));
        }
        if !source && effects.iter().any(|effect| effect.target == "source") {
            return Err(format!(
                "{ctx}: a source effect requires the center to be absent"
            ));
        }
        if source && (!rule.guards.is_empty() || !rule.modifiers.is_empty()) {
            return Err(format!(
                "{ctx}: source rules cannot declare center-relative guards or modifiers"
            ));
        }
        for effect in effects {
            match effect.target.as_str() {
                "center" | "neighbor" | "neighbors" | "source" => {}
                other => return Err(format!("{ctx}: unknown effect target '{other}'")),
            }
            match effect.select_mode.as_deref() {
                None | Some("all") | Some("one") => {}
                Some(other) => {
                    return Err(format!(
                        "{ctx}: select_mode must be 'all' or 'one', not '{other}'"
                    ));
                }
            }
            if effect.select_mode.is_some() && effect.target != "neighbors" {
                return Err(format!(
                    "{ctx}: select_mode is valid only for target = 'neighbors'"
                ));
            }
            if effect.target == "source" && effect.select.is_some() {
                return Err(format!("{ctx}: a source effect takes no selector"));
            }
        }
        validate_rate_for_strategy(rule, &deck.execution.strategy, &ctx)?;
    }

    for observable in &deck.observables.series {
        match observable.kind.as_str() {
            "state_counts" | "event_rates" | "rate_spectra" | "cluster_sizes" | "snapshot" => {}
            "interface_roughness" => {
                let axis = observable
                    .parameters
                    .get("axis")
                    .and_then(toml::Value::as_integer)
                    .ok_or("interface_roughness observable requires integer axis")?;
                if !(0..=2).contains(&axis) {
                    return Err("interface_roughness axis must be 0, 1, or 2".to_string());
                }
            }
            other => return Err(format!("unknown observable kind '{other}'")),
        }
    }
    Ok(())
}

fn validate_structure(structure: &StructureV2) -> Result<(), String> {
    if structure.kind.as_deref().unwrap_or("cell") != "grid" {
        return Ok(());
    }
    let family = structure
        .grid
        .as_deref()
        .ok_or("grid structure requires grid")?;
    let dimensions = structure
        .dims
        .as_ref()
        .ok_or("grid structure requires dims")?;
    let boundary = structure
        .boundary
        .as_ref()
        .ok_or("grid structure requires boundary")?;
    let expected_dimensions = match family {
        "square" | "hex" => 2,
        "cubic" => 3,
        other => return Err(format!("unknown grid family '{other}'")),
    };
    if dimensions.len() != expected_dimensions || dimensions.contains(&0) {
        return Err(format!(
            "{family} grid requires {expected_dimensions} nonzero dimensions"
        ));
    }
    if boundary.len() != expected_dimensions {
        return Err(format!(
            "{family} grid boundary must have {expected_dimensions} values (one per axis), got {}",
            boundary.len()
        ));
    }
    if boundary
        .iter()
        .any(|value| !matches!(value.as_str(), "periodic" | "open" | "fixed"))
    {
        return Err(format!(
            "{family} grid boundary values must be one of 'periodic', 'open', or 'fixed'"
        ));
    }
    match (family, structure.neighborhood.as_deref()) {
        ("hex", None) => {}
        ("hex", Some(_)) => return Err("hex grid has a fixed neighborhood".to_string()),
        (_, Some("moore" | "von_neumann")) => {}
        (_, None) => return Err("square/cubic grid requires neighborhood".to_string()),
        (_, Some(other)) => return Err(format!("unknown grid neighborhood '{other}'")),
    }
    Ok(())
}

fn validate_rate_for_strategy(
    rule: &ReactionSpec,
    strategy: &str,
    ctx: &str,
) -> Result<(), String> {
    let rate = rule.rate.as_ref();
    let variants = rate.map_or(0, |rate| {
        usize::from(rate.constant.is_some())
            + usize::from(rate.arrhenius.is_some())
            + usize::from(rate.eyring.is_some())
            + usize::from(rate.energy.is_some())
            + usize::from(rate.probability.is_some())
    });
    match strategy {
        "ctmc"
            if variants == 1
                && rate.is_some_and(|rate| {
                    rate.constant.is_some() || rate.arrhenius.is_some() || rate.eyring.is_some()
                }) => {}
        "synchronous" if variants == 0 => {}
        "metropolis" if variants == 1 && rate.is_some_and(|rate| rate.energy.is_some()) => {}
        "pca" if variants == 1 && rate.is_some_and(|rate| rate.probability.is_some()) => {
            let probability = rate.and_then(|rate| rate.probability).unwrap();
            if !probability.is_finite() || !(0.0..=1.0).contains(&probability) {
                return Err(format!("{ctx}: PCA probability must be in [0, 1]"));
            }
        }
        _ => {
            return Err(format!(
                "{ctx}: rate variant does not match execution strategy '{strategy}'"
            ));
        }
    }
    Ok(())
}

fn validate_execution(execution: &ExecutionSpec) -> Result<(), String> {
    match execution.strategy.as_str() {
        "ctmc" | "synchronous" | "metropolis" | "pca" => {}
        other => return Err(format!("unknown execution strategy '{other}'")),
    }
    if execution.strategy == "synchronous"
        && execution
            .synchronous
            .as_ref()
            .map(|spec| spec.conflict_resolution.as_str())
            != Some("first_match")
    {
        return Err(
            "[execution.synchronous] conflict_resolution = 'first_match' is required; \
             when multiple rules match one site, the first rule in definition order wins"
                .to_string(),
        );
    }
    if let Some(temperature) = execution
        .metropolis
        .as_ref()
        .and_then(|spec| spec.temperature)
    {
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(
                "[execution.metropolis] temperature must be finite and positive".to_string(),
            );
        }
    }
    if let Some(time) = execution.stop.time {
        if !time.is_finite() || time < 0.0 {
            return Err("[execution.stop] time must be finite and non-negative".to_string());
        }
    }
    if execution.ensemble.n_replicas == 0 {
        return Err("[execution.ensemble] n_replicas must be at least 1".to_string());
    }
    Ok(())
}

impl<'de> Deserialize<'de> for DeckFile {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        use serde::de::Error as _;

        let value = toml::Value::deserialize(deserializer)?;
        let schema = value
            .get("deck")
            .and_then(|deck| deck.get("schema"))
            .and_then(toml::Value::as_integer);
        match schema {
            None => {
                let v1: DeckV1 = value.try_into().map_err(D::Error::custom)?;
                Ok(Self::from_v1(v1))
            }
            Some(2) => {
                let v2: DeckV2 = value.try_into().map_err(D::Error::custom)?;
                Self::from_v2(v2).map_err(D::Error::custom)
            }
            Some(other) => Err(D::Error::custom(format!(
                "unsupported deck schema {other}; expected absent (v1) or 2"
            ))),
        }
    }
}

/// One build-time pass. The operation applies to the *center* site — once,
/// or once per neighbor matching `foreach` (in adjacency order), so
/// "step the map per missing cation" and "increment per terminal OH" are
/// both expressible. `probability` and `sites` are the substitution/defect
/// fill rules of design doc §6: "Fe replaces Al on Al_oct at 0.02" is a
/// pass with `probability = 0.02` and `set = <an Fe-occupant state>`.
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
    /// Apply the op at each qualifying site with this probability, in
    /// [0, 1]. Draws are deterministic given (deck, run seed): one shared
    /// init RNG stream, consumed in pass order then site-index order, one
    /// draw per site that passes every other filter. With `foreach`, the
    /// draw gates the whole site (all repetitions), not each repetition.
    #[serde(default)]
    pub probability: Option<f64>,
    /// Restrict to an explicit site list; each entry is `[a, b, c, t]` —
    /// cell coordinates plus template-site index (the order sites appear
    /// under `[[cell.sites]]`). Composes with every other filter.
    #[serde(default)]
    pub sites: Option<Vec<[usize; 4]>>,
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
    /// Change the site's kind after applying the state operation. Used by
    /// seeded binned-disorder init passes.
    #[serde(default)]
    pub rekind: Option<String>,
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
    /// Schema marker. The v1 shim normalizes an absent marker to `Some(2)`.
    #[serde(default)]
    pub schema: Option<u32>,
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
pub struct DefectSpec {
    /// `"screw"` or `"edge"` (isotropic θ-averaged edge form).
    #[serde(rename = "type")]
    pub kind: String,
    /// Dislocation line direction: 0 = a, 1 = b, 2 = c.
    pub line_axis: u8,
    /// A point the line passes through, in CELL coordinates (the component
    /// along `line_axis` is irrelevant).
    pub at: [f64; 3],
    /// |Burgers vector|, Å. Required unless `strain_prefactor` is given.
    #[serde(default)]
    pub burgers: Option<f64>,
    /// Shear modulus μ, GPa. Required unless `strain_prefactor` is given.
    #[serde(default)]
    pub shear_modulus: Option<f64>,
    /// Poisson's ratio ν (edge only; default 0.25).
    #[serde(default)]
    pub poisson: Option<f64>,
    /// Continuum cutoff / hollow-core clamp, Å: r is clamped to this.
    /// Default = burgers.
    #[serde(default)]
    pub core_radius: Option<f64>,
    /// Optional hard cap on u per site, in deck energy units.
    #[serde(default)]
    pub cap: Option<f64>,
    /// Directly set A in u = A/r² (deck-energy·Å²), overriding the
    /// physical inputs — for illustrative decks and tests.
    #[serde(default)]
    pub strain_prefactor: Option<f64>,
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
    /// Absent only for a per-site `target = "source"` rule.
    #[serde(default)]
    pub center: Option<CenterSpec>,
    #[serde(default)]
    pub guards: Vec<SelectorSpec>,
    #[serde(default)]
    pub rate: Option<RateSpec>,
    /// Deterministic-CA truth-table mode. Parsed in B2; executed in B3.
    #[serde(default)]
    pub truth: Option<String>,
    /// Species drawn from solution: each contributes activity and Δμ
    /// factors to the forward rate (design doc §4).
    #[serde(default)]
    pub consumes: Vec<String>,
    /// Species released to solution (bookkeeping in v1; reverse rates come
    /// from explicit reverse reactions until the P1 auto-reverse lands).
    #[serde(default)]
    pub produces: Vec<String>,
    /// Strain coupling: `strain = { scale = s }` adds `s · u_center` to the
    /// activation energy (docs/STRAIN.md §2.2; dissolution-forward = −β,
    /// reverse = +(1−β)).
    #[serde(default)]
    pub strain: Option<StrainSpec>,
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
    /// Bonds with this label are invisible to the selector: skipped at
    /// distance 1, never traversed on either hop at distance 2 (distances
    /// are measured on the filtered graph). Keeps chemistry from
    /// propagating through non-reactive contacts — e.g. a surface-
    /// reachability guard that must not cross anhydrous interlayer
    /// hydrogen bonds. Mutually exclusive with `label`.
    #[serde(default)]
    pub exclude_label: Option<String>,
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

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StrainSpec {
    /// Dimensionless multiplier on the center's strain energy.
    pub scale: f64,
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
    #[serde(default)]
    pub energy: Option<EnergySpec>,
    #[serde(default)]
    pub probability: Option<f64>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EnergySpec {
    pub delta: f64,
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
    /// Neighbor selection: `all` (default) or exactly `one` uniform match.
    #[serde(default)]
    pub select_mode: Option<String>,
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
