//! Deck → runtime compilation: intern names to dense ids, expand bonds,
//! build state sets, fold thermodynamic factors, and validate everything
//! validatable before the engine ever runs (design doc §3, §4).

use std::collections::BTreeMap;

use petra_core::crystal::{Cell, KindId, TemplateBond, TemplateSite, UnitCell, NO_LABEL};
use petra_core::lattice::{Boundary, Lattice};
use petra_core::rate::{RateExpr, R_KCAL};
use petra_core::reaction::{
    count_matches, Branch, Effect, EffectOp, EffectTarget, Guard, Modifier, ModifierKind,
    NeighborSelect, Reaction,
};
use petra_core::state::{StateId, StateSet};
use petra_core::Engine;
use rand::{Rng as _, SeedableRng};
use rand_pcg::Pcg64Mcg;

use crate::schema::{DeckFile, EffectSpec, RateSpec, SeedPolicy, SelectorSpec, StructureKind};

#[derive(Debug, thiserror::Error)]
#[error("deck error: {0}")]
pub struct CompileError(pub String);

fn err<T>(msg: impl Into<String>) -> Result<T, CompileError> {
    Err(CompileError(msg.into()))
}

/// Everything the engine and the reporting layer need, with the name
/// tables kept for output (the runtime itself sees only dense ids).
#[derive(Debug)]
pub struct CompiledDeck {
    pub name: String,
    pub unit_cell: UnitCell,
    pub kinds_per_template: Vec<KindId>,
    pub kind_names: Vec<String>,
    /// `"Kind.state"`, indexed by `StateId`.
    pub state_names: Vec<String>,
    /// Occupant species name per state (`None` = vacant), indexed by
    /// `StateId` — for snapshots/exports that need chemical identity.
    pub state_occupants: Vec<Option<String>>,
    /// Bond label names, indexed by the interned label id.
    pub label_names: Vec<String>,
    pub n_states: usize,
    /// Contiguous StateId range (start, count) per kind, for Shift effects.
    pub kind_state_ranges: Vec<(u16, u16)>,
    pub initial_per_template: Vec<StateId>,
    pub init_passes: Vec<InitPass>,
    pub defects: Vec<CompiledDefect>,
    pub reactions: Vec<Reaction>,
    pub temperature: f64,
    pub dims: [usize; 3],
    pub boundary: [Boundary; 3],
    pub steps: u64,
    pub seed: u64,
    pub report_every: u64,
}

/// One compiled build-time pass (design doc §3.1 fill rules): sweep all
/// sites in index order, writes immediately visible — the legacy
/// TerminateSurface convention.
#[derive(Debug)]
pub struct InitPass {
    pub name: String,
    pub center_kind: Option<KindId>,
    pub center_states: StateSet,
    /// (axis, min, max) inclusive cell-coordinate slab.
    pub region: Option<(u8, usize, usize)>,
    /// Apply at each qualifying site with this probability (schema doc for
    /// the draw-order contract); `None` = always.
    pub probability: Option<f64>,
    /// Explicit `[a, b, c, t]` site list restriction; `None` = all sites.
    pub sites: Option<Vec<[usize; 4]>>,
    pub guards: Vec<Guard>,
    /// Apply the op once per matching neighbor instead of once.
    pub foreach: Option<NeighborSelect>,
    pub op: EffectOp,
}

/// One compiled line defect: u(r) = prefactor / max(r, core_radius)²,
/// capped, superposed over all defects (docs/STRAIN.md §5). All values in
/// internal units (kcal/mol, Å).
#[derive(Debug, Clone)]
pub struct CompiledDefect {
    pub line_axis: u8,
    /// A point on the line, cell coordinates.
    pub at: [f64; 3],
    /// A in kcal·Å²/mol.
    pub prefactor: f64,
    pub core_radius: f64,
    pub cap: Option<f64>,
}

impl CompiledDeck {
    /// Instantiate the lattice (uniform per-kind fill), run the init
    /// passes, compute defect strain fields, and build the engine.
    pub fn build_engine(&self, seed_override: Option<u64>) -> Result<Engine, CompileError> {
        let seed = seed_override.unwrap_or(self.seed);
        let initial = self.initial_per_template.clone();
        let mut lattice = Lattice::build(&self.unit_cell, self.dims, self.boundary, |t| initial[t]);
        let kinds: Vec<KindId> = lattice
            .template_index
            .iter()
            .map(|&t| self.kinds_per_template[t as usize])
            .collect();
        self.run_init(&mut lattice, &kinds, seed)?;
        self.compute_strain(&mut lattice)?;
        Ok(Engine::new(
            lattice,
            &self.kinds_per_template,
            self.kind_names.len(),
            self.kind_state_ranges.clone(),
            self.reactions.clone(),
            self.temperature,
            seed,
        ))
    }

    /// Superpose every defect's analytic field into `lattice.strain`
    /// (docs/STRAIN.md §5). Perpendicular distance to the line uses
    /// minimum-image displacement on periodic axes.
    fn compute_strain(&self, lattice: &mut Lattice) -> Result<(), CompileError> {
        if self.defects.is_empty() {
            return Ok(());
        }
        let m = self.unit_cell.cell.matrix();
        let mul = |f: [f64; 3]| -> [f64; 3] {
            [
                m[0][0] * f[0] + m[0][1] * f[1] + m[0][2] * f[2],
                m[1][0] * f[0] + m[1][1] * f[1] + m[1][2] * f[2],
                m[2][0] * f[0] + m[2][1] * f[1] + m[2][2] * f[2],
            ]
        };
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

        for d in &self.defects {
            let ax = d.line_axis as usize;
            let mut dir = [m[0][ax], m[1][ax], m[2][ax]];
            let len = dot(dir, dir).sqrt();
            for c in &mut dir {
                *c /= len;
            }
            let q = mul(d.at);
            for s in 0..lattice.len() {
                let (cell, t) = lattice.coords(s);
                let p = self.unit_cell.cell.to_cartesian(
                    self.unit_cell.sites[t].frac,
                    [cell[0] as i32, cell[1] as i32, cell[2] as i32],
                );
                let v = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
                let mut f = self.unit_cell.cell.to_fractional(v).map_err(CompileError)?;
                for (axis, fc) in f.iter_mut().enumerate() {
                    if self.boundary[axis] == Boundary::Periodic {
                        let period = self.dims[axis] as f64;
                        *fc -= period * (*fc / period).round();
                    }
                }
                let w = mul(f);
                let along = dot(w, dir);
                let r2 = (dot(w, w) - along * along).max(0.0);
                let r = r2.sqrt().max(d.core_radius);
                let mut u = d.prefactor / (r * r);
                if let Some(cap) = d.cap {
                    u = u.min(cap);
                }
                lattice.strain[s] += u;
            }
        }
        Ok(())
    }

    fn run_init(
        &self,
        lattice: &mut Lattice,
        kinds: &[KindId],
        seed: u64,
    ) -> Result<(), CompileError> {
        let mut scratch = Vec::new();
        // One RNG stream for every probabilistic pass, decorrelated from the
        // dynamics stream by a fixed salt; draw order is the documented
        // contract (schema.rs): pass order, then site-index order, one draw
        // per site that passed every other filter.
        const INIT_STREAM_SALT: u64 = 0x9E37_79B9_7F4A_7C15;
        let mut rng = Pcg64Mcg::seed_from_u64(seed ^ INIT_STREAM_SALT);
        for pass in &self.init_passes {
            // An explicit site list, sorted+deduped, IS the iteration order —
            // site-index order, so the RNG draw contract is unchanged.
            let site_list: Option<Vec<usize>> = pass.sites.as_ref().map(|list| {
                let mut v: Vec<usize> = list
                    .iter()
                    .map(|&[a, b, c, t]| lattice.index([a, b, c], t))
                    .collect();
                v.sort_unstable();
                v.dedup();
                v
            });
            let sweep: Box<dyn Iterator<Item = usize>> = match &site_list {
                Some(list) => Box::new(list.iter().copied()),
                None => Box::new(0..lattice.len()),
            };
            for s in sweep {
                if let Some((axis, min, max)) = pass.region {
                    let (cell, _) = lattice.coords(s);
                    let coord = cell[axis as usize];
                    if coord < min || coord > max {
                        continue;
                    }
                }
                if let Some(k) = pass.center_kind {
                    if kinds[s] != k {
                        continue;
                    }
                }
                if !pass.center_states.contains(lattice.states[s]) {
                    continue;
                }
                let guards_ok = pass.guards.iter().all(|g| {
                    let n = count_matches(lattice, kinds, s, &g.select, &mut scratch);
                    n >= g.min && n <= g.max
                });
                if !guards_ok {
                    continue;
                }
                if let Some(p) = pass.probability {
                    if rng.gen::<f64>() >= p {
                        continue;
                    }
                }
                let times = match &pass.foreach {
                    None => 1,
                    Some(sel) => count_matches(lattice, kinds, s, sel, &mut scratch),
                };
                let range = self.kind_state_ranges[kinds[s].0 as usize];
                for _ in 0..times {
                    match pass.op.resolve(lattice.states[s], range) {
                        Ok(Some(new_state)) => lattice.states[s] = new_state,
                        Ok(None) => {}
                        Err(reason) => {
                            return err(format!(
                                "init pass '{}' failed at site {s}: {reason}",
                                pass.name
                            ))
                        }
                    }
                }
            }
        }
        Ok(())
    }
}

/// Name-resolution tables built once, used everywhere.
struct Names {
    kind_ids: BTreeMap<String, u16>,
    /// state name → all (kind, id) declaring it; qualified "Kind.state" keys
    /// are also present, each with exactly one entry.
    state_refs: BTreeMap<String, Vec<(u16, StateId)>>,
    /// per kind: state name → id.
    kind_states: Vec<BTreeMap<String, StateId>>,
    aliases: BTreeMap<String, Vec<String>>,
    labels: BTreeMap<String, u16>,
    n_states: usize,
}

impl Names {
    /// Resolve one state reference (plain, qualified, or `@alias`) into ids,
    /// pushed into `out`.
    fn resolve_ref(
        &self,
        r: &str,
        out: &mut Vec<(u16, StateId)>,
        depth: u8,
    ) -> Result<(), CompileError> {
        if depth > 4 {
            return err(format!("alias nesting too deep at '{r}'"));
        }
        if r == "*" {
            // Wildcard: every declared state. The caller's kind restriction
            // (if any) still applies, turning this into a degree counter.
            for ids in self.state_refs.values() {
                out.extend(ids.iter().copied());
            }
            return Ok(());
        }
        if let Some(alias) = r.strip_prefix('@') {
            let members = self
                .aliases
                .get(alias)
                .ok_or_else(|| CompileError(format!("unknown alias '@{alias}'")))?;
            for m in members {
                self.resolve_ref(m, out, depth + 1)?;
            }
            return Ok(());
        }
        match self.state_refs.get(r) {
            Some(ids) => {
                out.extend(ids.iter().copied());
                Ok(())
            }
            None => err(format!("unknown state '{r}'")),
        }
    }

    /// Resolve a list of refs into a `StateSet`, optionally restricted to
    /// one kind (with an error if a plain ref is ambiguous *within* the
    /// restriction's scope — i.e. resolves to nothing there).
    fn state_set(
        &self,
        refs: &[String],
        restrict_kind: Option<u16>,
        ctx: &str,
    ) -> Result<StateSet, CompileError> {
        let mut ids = Vec::new();
        for r in refs {
            self.resolve_ref(r, &mut ids, 0)?;
        }
        let mut set = StateSet::new(self.n_states);
        let mut any = false;
        for (kind, id) in ids {
            if restrict_kind.is_none_or(|k| k == kind) {
                set.insert(id);
                any = true;
            }
        }
        if !any {
            return err(format!("{ctx}: state list resolves to nothing in scope"));
        }
        Ok(set)
    }
}

/// kcal/mol per deck energy unit (design doc §4: the runtime is canonical
/// kcal/mol; decks choose their unit and are converted once, here).
fn energy_factor(units: Option<&str>) -> Result<f64, CompileError> {
    match units.unwrap_or("kcal/mol") {
        "kcal/mol" => Ok(1.0),
        "kJ/mol" => Ok(1.0 / 4.184),
        "eV" => Ok(23.060_548),
        other => err(format!(
            "unknown units '{other}' (expected kcal/mol, kJ/mol, or eV)"
        )),
    }
}

pub fn replica_seed(seed: u64, replica: u64, policy: SeedPolicy) -> u64 {
    match policy {
        SeedPolicy::Increment => seed.wrapping_add(replica),
        SeedPolicy::Hash => splitmix64(splitmix64(seed) ^ replica),
    }
}

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

pub fn compile(deck: &DeckFile) -> Result<CompiledDeck, CompileError> {
    if deck.structure_kind != StructureKind::Cell {
        return err(
            "grid structures parse in schema v2 but compile with B3's conformance strategies",
        );
    }
    if deck.execution.strategy != "ctmc" {
        return err(format!(
            "strategy '{}' is not implemented in B2 (only 'ctmc' is available)",
            deck.execution.strategy
        ));
    }
    if deck.execution.ensemble.n_replicas == 0 {
        return err("execution.ensemble.n_replicas must be at least 1");
    }
    let eunit = energy_factor(deck.deck.units.as_deref())?;

    // --- species ---
    let mut species = BTreeMap::new();
    for s in &deck.species {
        if species.insert(s.name.clone(), ()).is_some() {
            return err(format!("duplicate species '{}'", s.name));
        }
        if s.name == "vacant" {
            return err("'vacant' is reserved and cannot be a species name");
        }
    }

    // --- kinds and states: intern in declaration order ---
    let mut names = Names {
        kind_ids: BTreeMap::new(),
        state_refs: BTreeMap::new(),
        kind_states: Vec::new(),
        aliases: deck.aliases.clone(),
        labels: BTreeMap::new(),
        n_states: 0,
    };
    let mut state_names = Vec::new();
    let mut state_occupants = Vec::new();
    let mut kind_names = Vec::new();
    let mut initial_by_kind = Vec::new();
    let mut kind_state_ranges = Vec::new();

    for (ki, k) in deck.kinds.iter().enumerate() {
        let ki = ki as u16;
        kind_state_ranges.push((state_names.len() as u16, k.states.len() as u16));
        if names.kind_ids.insert(k.name.clone(), ki).is_some() {
            return err(format!("duplicate kind '{}'", k.name));
        }
        kind_names.push(k.name.clone());
        let mut per_kind = BTreeMap::new();
        for st in &k.states {
            if st.occupant != "vacant" && !species.contains_key(&st.occupant) {
                return err(format!(
                    "state '{}.{}' has unknown occupant '{}'",
                    k.name, st.name, st.occupant
                ));
            }
            let id = StateId(
                u16::try_from(state_names.len())
                    .map_err(|_| CompileError("more than 65535 states".into()))?,
            );
            if per_kind.insert(st.name.clone(), id).is_some() {
                return err(format!(
                    "duplicate state '{}' in kind '{}'",
                    st.name, k.name
                ));
            }
            state_names.push(format!("{}.{}", k.name, st.name));
            state_occupants.push(if st.occupant == "vacant" {
                None
            } else {
                Some(st.occupant.clone())
            });
            names
                .state_refs
                .entry(st.name.clone())
                .or_default()
                .push((ki, id));
            names
                .state_refs
                .insert(format!("{}.{}", k.name, st.name), vec![(ki, id)]);
        }
        let initial = *per_kind.get(&k.initial).ok_or_else(|| {
            CompileError(format!(
                "kind '{}': initial state '{}' is not declared",
                k.name, k.initial
            ))
        })?;
        initial_by_kind.push(initial);
        names.kind_states.push(per_kind);
    }
    names.n_states = state_names.len();

    // --- unit cell: geometry, sites, then bonds expanded onto both endpoints ---
    let params = (
        deck.cell.a,
        deck.cell.b,
        deck.cell.c,
        deck.cell.alpha,
        deck.cell.beta,
        deck.cell.gamma,
    );
    let cell =
        match (deck.cell.matrix, params) {
            (Some(m), (None, None, None, None, None, None)) => Cell::from_matrix(m),
            (None, (Some(a), Some(b), Some(c), Some(al), Some(be), Some(ga))) => {
                Cell::from_params(a, b, c, al, be, ga)
            }
            _ => return err(
                "cell geometry must be either all of a,b,c,alpha,beta,gamma or a matrix, not a mix",
            ),
        };
    let mut tsites = Vec::new();
    let mut kinds_per_template = Vec::new();
    let mut initial_per_template = Vec::new();
    for (si, s) in deck.cell.sites.iter().enumerate() {
        let ki = *names
            .kind_ids
            .get(&s.kind)
            .ok_or_else(|| CompileError(format!("cell site has unknown kind '{}'", s.kind)))?;
        kinds_per_template.push(KindId(ki));
        initial_per_template.push(initial_by_kind[ki as usize]);
        let frac = match (s.frac, s.cart) {
            (Some(f), None) => f,
            (None, Some(c)) => cell
                .to_fractional(c)
                .map_err(|e| CompileError(format!("cell site {si}: {e}")))?,
            _ => return err(format!("cell site {si}: exactly one of frac/cart")),
        };
        tsites.push(TemplateSite {
            kind: KindId(ki),
            frac,
            bonds: Vec::new(),
        });
    }
    for b in &deck.cell.bonds {
        if b.i >= tsites.len() || b.j >= tsites.len() {
            return err(format!("bond ({}, {}) references a missing site", b.i, b.j));
        }
        if b.i == b.j && b.dcell == [0, 0, 0] {
            return err(format!("site {} bonded to itself in the same cell", b.i));
        }
        let label = match &b.label {
            None => NO_LABEL,
            Some(l) => {
                let next = names.labels.len() as u16;
                *names.labels.entry(l.clone()).or_insert(next)
            }
        };
        let fwd = TemplateBond {
            to: b.j,
            dcell: b.dcell,
            label,
        };
        let rev = TemplateBond {
            to: b.i,
            dcell: [-b.dcell[0], -b.dcell[1], -b.dcell[2]],
            label,
        };
        if !tsites[b.i].bonds.contains(&fwd) {
            tsites[b.i].bonds.push(fwd);
        }
        if !tsites[b.j].bonds.contains(&rev) {
            tsites[b.j].bonds.push(rev);
        }
    }
    let unit_cell = UnitCell {
        cell,
        sites: tsites,
    };
    unit_cell
        .check_reciprocity()
        .map_err(|e| CompileError(format!("internal bond expansion bug: {e}")))?;

    // --- thermo ---
    let temperature = deck.thermo.temperature;
    if temperature <= 0.0 {
        return err("temperature must be positive");
    }
    let rt = R_KCAL * temperature;

    // --- boundary ---
    let mut boundary = [Boundary::Periodic; 3];
    for (i, b) in deck.lattice.boundary.iter().enumerate() {
        boundary[i] = match b.as_str() {
            "periodic" => Boundary::Periodic,
            "open" => Boundary::Open,
            "fixed" => Boundary::Fixed,
            other => return err(format!("unknown boundary '{other}'")),
        };
    }

    // --- reactions ---
    let mut reactions = Vec::new();
    for r in &deck.reactions {
        let ctx = format!("reaction '{}'", r.name);
        let center_kind = *names.kind_ids.get(&r.center.kind).ok_or_else(|| {
            CompileError(format!("{ctx}: unknown center kind '{}'", r.center.kind))
        })?;
        let center_states = names.state_set(&r.center.state, Some(center_kind), &ctx)?;

        let mut guards = Vec::new();
        for g in &r.guards {
            let select = compile_selector(&names, g, &ctx)?;
            guards.push(Guard {
                select,
                min: g.min.unwrap_or(1),
                max: g.max.unwrap_or(u32::MAX),
            });
        }

        let rate = compile_rate(&r.rate, eunit, &ctx)?;

        // Fold solution coupling into ln_thermo (design doc §4):
        // each consumed species contributes ln(activity) + Δμ/RT.
        let mut ln_thermo = 0.0;
        for sp in &r.consumes {
            if !species.contains_key(sp) {
                return err(format!("{ctx}: consumes unknown species '{sp}'"));
            }
            let activity = deck.thermo.activity.get(sp).copied().unwrap_or(1.0);
            if activity <= 0.0 {
                return err(format!("{ctx}: activity of '{sp}' must be positive"));
            }
            ln_thermo +=
                activity.ln() + deck.thermo.mu.get(sp).copied().unwrap_or(0.0) * eunit / rt;
        }
        for sp in &r.produces {
            if !species.contains_key(sp) {
                return err(format!("{ctx}: produces unknown species '{sp}'"));
            }
        }

        let mut modifiers = Vec::new();
        for m in &r.modifiers {
            let select = compile_selector(&names, &m.select, &ctx)?;
            let kind = match (&m.per_match, &m.by_count, &m.when) {
                (Some(p), None, None) => ModifierKind::PerMatch { dea: p.dea * eunit },
                (None, Some(b), None) => {
                    if b.dea.is_empty() {
                        return err(format!("{ctx}: by_count.dea must be non-empty"));
                    }
                    ModifierKind::ByCount {
                        dea: b.dea.iter().map(|d| d * eunit).collect(),
                    }
                }
                (None, None, Some(w)) => ModifierKind::When {
                    min: w.min.unwrap_or(1),
                    max: w.max.unwrap_or(u32::MAX),
                    dea: w.dea.unwrap_or(0.0) * eunit,
                    factor: w.factor.unwrap_or(1.0),
                },
                _ => {
                    return err(format!(
                        "{ctx}: modifier needs exactly one of per_match/by_count/when"
                    ))
                }
            };
            modifiers.push(Modifier { select, kind });
        }

        let branches = match (r.effects.is_empty(), r.branches.is_empty()) {
            (false, true) => vec![Branch {
                weight: 1.0,
                effects: compile_effects(&names, &r.effects, center_kind, &ctx)?,
            }],
            (true, false) => {
                let mut out = Vec::new();
                for b in &r.branches {
                    if b.weight <= 0.0 {
                        return err(format!("{ctx}: branch weight must be positive"));
                    }
                    out.push(Branch {
                        weight: b.weight,
                        effects: compile_effects(&names, &b.effects, center_kind, &ctx)?,
                    });
                }
                out
            }
            _ => return err(format!("{ctx}: declare exactly one of effects/branches")),
        };

        reactions.push(Reaction {
            name: r.name.clone(),
            center_kind: KindId(center_kind),
            center_states,
            guards,
            rate,
            ln_thermo,
            strain_scale: r.strain.as_ref().map(|s| s.scale).unwrap_or(0.0),
            modifiers,
            branches,
        });
    }

    // --- defects (docs/STRAIN.md §5) ---
    // kcal/mol per GPa·Å³: 1 GPa·Å³ = 1e-21 J → × N_A / 4184.
    const KCAL_PER_GPA_A3: f64 = 0.143_932_6;
    let mm = unit_cell.cell.matrix();
    let cell_volume = (mm[0][0] * (mm[1][1] * mm[2][2] - mm[1][2] * mm[2][1])
        - mm[0][1] * (mm[1][0] * mm[2][2] - mm[1][2] * mm[2][0])
        + mm[0][2] * (mm[1][0] * mm[2][1] - mm[1][1] * mm[2][0]))
        .abs();
    let omega = cell_volume / unit_cell.sites.len().max(1) as f64;
    let mut defects = Vec::new();
    for (di, d) in deck.defects.iter().enumerate() {
        let ctx = format!("defect {di}");
        if d.line_axis > 2 {
            return err(format!("{ctx}: line_axis must be 0, 1, or 2"));
        }
        let edge = match d.kind.as_str() {
            "screw" => false,
            "edge" => true,
            other => return err(format!("{ctx}: type must be screw or edge, not '{other}'")),
        };
        let prefactor = match (d.strain_prefactor, d.burgers, d.shear_modulus) {
            (Some(a), _, _) => a * eunit,
            (None, Some(b), Some(mu)) => {
                let nu = d.poisson.unwrap_or(0.25);
                if edge && nu >= 1.0 {
                    return err(format!("{ctx}: poisson must be < 1"));
                }
                let mut a = KCAL_PER_GPA_A3 * mu * b * b * omega
                    / (8.0 * std::f64::consts::PI * std::f64::consts::PI);
                if edge {
                    a /= 1.0 - nu;
                }
                a
            }
            _ => {
                return err(format!(
                    "{ctx}: give burgers + shear_modulus, or strain_prefactor"
                ))
            }
        };
        let core_radius = match d.core_radius.or(d.burgers) {
            Some(r) if r > 0.0 => r,
            _ => return err(format!("{ctx}: core_radius (or burgers) must be positive")),
        };
        defects.push(CompiledDefect {
            line_axis: d.line_axis,
            at: d.at,
            prefactor,
            core_radius,
            cap: d.cap.map(|c| c * eunit),
        });
    }

    let dims = deck.lattice.dims;
    if dims.contains(&0) {
        return err("lattice dims must be nonzero");
    }

    // --- init passes ---
    let mut init_passes = Vec::new();
    for p in &deck.init {
        let ctx = format!("init pass '{}'", p.name);
        let center_kind = match &p.center.kind {
            None => None,
            Some(k) => Some(
                *names
                    .kind_ids
                    .get(k)
                    .ok_or_else(|| CompileError(format!("{ctx}: unknown center kind '{k}'")))?,
            ),
        };
        let center_states = names.state_set(&p.center.state, center_kind, &ctx)?;
        let region = match &p.region {
            None => None,
            Some(r) => {
                if r.axis > 2 {
                    return err(format!("{ctx}: region axis must be 0, 1, or 2"));
                }
                Some((r.axis, r.min.unwrap_or(0), r.max.unwrap_or(usize::MAX)))
            }
        };
        if let Some(prob) = p.probability {
            if !(0.0..=1.0).contains(&prob) {
                return err(format!("{ctx}: probability must be in [0, 1], got {prob}"));
            }
        }
        if let Some(sites) = &p.sites {
            if sites.is_empty() {
                return err(format!("{ctx}: sites list must be non-empty"));
            }
            for &[a, b, c, t] in sites {
                if a >= dims[0] || b >= dims[1] || c >= dims[2] {
                    return err(format!(
                        "{ctx}: site [{a}, {b}, {c}, {t}] outside lattice dims {dims:?}"
                    ));
                }
                if t >= unit_cell.sites.len() {
                    return err(format!(
                        "{ctx}: site [{a}, {b}, {c}, {t}] names template site {t}, \
                         but the cell has only {}",
                        unit_cell.sites.len()
                    ));
                }
            }
        }
        let mut guards = Vec::new();
        for g in &p.guards {
            let select = compile_selector(&names, g, &ctx)?;
            guards.push(Guard {
                select,
                min: g.min.unwrap_or(1),
                max: g.max.unwrap_or(u32::MAX),
            });
        }
        let foreach = match &p.foreach {
            None => None,
            Some(sel) => Some(compile_selector(&names, sel, &ctx)?),
        };
        // Init maps default to skip (termination maps leave unlisted
        // states alone).
        let op = compile_op(
            &names,
            center_kind,
            &p.set,
            &p.shift,
            &p.map,
            &p.missing,
            false,
            &ctx,
        )?;
        init_passes.push(InitPass {
            name: p.name.clone(),
            center_kind: center_kind.map(KindId),
            center_states,
            region,
            probability: p.probability,
            sites: p.sites.clone(),
            guards,
            foreach,
            op,
        });
    }

    let mut label_names = vec![String::new(); names.labels.len()];
    for (name, &id) in &names.labels {
        label_names[id as usize] = name.clone();
    }

    Ok(CompiledDeck {
        name: deck.deck.name.clone(),
        unit_cell,
        kinds_per_template,
        kind_names,
        state_names,
        state_occupants,
        label_names,
        n_states: names.n_states,
        kind_state_ranges,
        initial_per_template,
        init_passes,
        defects,
        reactions,
        temperature,
        dims,
        boundary,
        steps: deck.simulation.steps,
        seed: deck.simulation.seed,
        report_every: deck.simulation.report_every.unwrap_or(0),
    })
}

fn compile_selector(
    names: &Names,
    s: &SelectorSpec,
    ctx: &str,
) -> Result<NeighborSelect, CompileError> {
    let distance = s.distance.unwrap_or(1);
    if !(1..=2).contains(&distance) {
        return err(format!("{ctx}: selector distance must be 1 or 2"));
    }
    if s.label.is_some() && distance != 1 {
        return err(format!("{ctx}: bond labels only select at distance 1"));
    }
    let kind = match &s.kind {
        None => None,
        Some(k) => Some(KindId(*names.kind_ids.get(k).ok_or_else(|| {
            CompileError(format!("{ctx}: selector names unknown kind '{k}'"))
        })?)),
    };
    let label = match &s.label {
        None => None,
        Some(l) => Some(*names.labels.get(l).ok_or_else(|| {
            CompileError(format!("{ctx}: selector names unknown bond label '{l}'"))
        })?),
    };
    let exclude_label = match &s.exclude_label {
        None => None,
        Some(l) => Some(*names.labels.get(l).ok_or_else(|| {
            CompileError(format!("{ctx}: selector excludes unknown bond label '{l}'"))
        })?),
    };
    if label.is_some() && exclude_label.is_some() {
        return err(format!(
            "{ctx}: a selector cannot both require and exclude a bond label"
        ));
    }
    let states = names.state_set(&s.state, kind.map(|k| k.0), ctx)?;
    Ok(NeighborSelect {
        distance,
        kind,
        label,
        exclude_label,
        frozen: s.frozen,
        states,
    })
}

fn compile_rate(r: &RateSpec, eunit: f64, ctx: &str) -> Result<RateExpr, CompileError> {
    match (r.constant, &r.arrhenius, &r.eyring) {
        (Some(k), None, None) => Ok(RateExpr::Constant { k }),
        (None, Some(a), None) => Ok(RateExpr::Arrhenius {
            prefactor: a.prefactor,
            ea: a.ea * eunit,
        }),
        (None, None, Some(e)) => Ok(RateExpr::Eyring {
            dh: e.dh * eunit,
            ds: e.ds * eunit,
        }),
        _ => err(format!(
            "{ctx}: rate needs exactly one of constant/arrhenius/eyring"
        )),
    }
}

/// Resolve one state name within a specific kind (plain or "Kind."-prefixed).
fn state_in_kind(names: &Names, kind: u16, name: &str, ctx: &str) -> Result<StateId, CompileError> {
    names
        .kind_states
        .get(kind as usize)
        .and_then(|m| {
            m.get(name).copied().or_else(|| {
                name.strip_prefix(&format!("{}.", names_kind(names, kind)))
                    .and_then(|plain| m.get(plain).copied())
            })
        })
        .ok_or_else(|| {
            CompileError(format!(
                "{ctx}: unknown state '{name}' for kind '{}'",
                names_kind(names, kind)
            ))
        })
}

/// Compile the set/shift/map trio into an EffectOp. `kind` scopes state-name
/// resolution and is required for `set`/`map`.
#[allow(clippy::too_many_arguments)]
fn compile_op(
    names: &Names,
    kind: Option<u16>,
    set: &Option<String>,
    shift: &Option<i32>,
    map: &Option<std::collections::BTreeMap<String, String>>,
    missing: &Option<String>,
    missing_default_error: bool,
    ctx: &str,
) -> Result<EffectOp, CompileError> {
    let missing_is_error = match missing.as_deref() {
        None => missing_default_error,
        Some("error") => true,
        Some("skip") => false,
        Some(other) => {
            return err(format!(
                "{ctx}: missing must be 'error' or 'skip', not '{other}'"
            ))
        }
    };
    match (set, shift, map) {
        (Some(name), None, None) => {
            let k = kind.ok_or_else(|| {
                CompileError(format!(
                    "{ctx}: 'set' needs a kind in scope to resolve '{name}'"
                ))
            })?;
            Ok(EffectOp::Set(state_in_kind(names, k, name, ctx)?))
        }
        (None, Some(n), None) => Ok(EffectOp::Shift(*n)),
        (None, None, Some(m)) => {
            let k =
                kind.ok_or_else(|| CompileError(format!("{ctx}: 'map' needs a kind in scope")))?;
            let mut entries = Vec::new();
            for (from, to) in m {
                entries.push((
                    state_in_kind(names, k, from, ctx)?,
                    state_in_kind(names, k, to, ctx)?,
                ));
            }
            if entries.is_empty() {
                return err(format!("{ctx}: map must be non-empty"));
            }
            Ok(EffectOp::Map {
                entries,
                missing_is_error,
            })
        }
        _ => err(format!("{ctx}: exactly one of set/shift/map")),
    }
}

fn compile_effects(
    names: &Names,
    effects: &[EffectSpec],
    center_kind: u16,
    ctx: &str,
) -> Result<Vec<Effect>, CompileError> {
    if effects.is_empty() {
        return err(format!("{ctx}: at least one effect is required"));
    }
    let mut out = Vec::new();
    for e in effects {
        let (target, target_kind) = match e.target.as_str() {
            "center" => {
                if e.select.is_some() {
                    return err(format!("{ctx}: a center effect takes no selector"));
                }
                (EffectTarget::Center, Some(center_kind))
            }
            t @ ("neighbor" | "neighbors") => {
                let sel = e.select.as_ref().ok_or_else(|| {
                    CompileError(format!("{ctx}: a neighbor effect requires a selector"))
                })?;
                let compiled = compile_selector(names, sel, ctx)?;
                if compiled.kind.is_none() && (e.set.is_some() || e.map.is_some()) {
                    return err(format!(
                        "{ctx}: neighbor-effect selectors must name a kind for set/map to resolve"
                    ));
                }
                let k = compiled.kind.map(|k| k.0);
                let target = if t == "neighbor" {
                    EffectTarget::FirstMatch(compiled)
                } else {
                    EffectTarget::AllMatches(compiled)
                };
                (target, k)
            }
            other => return err(format!("{ctx}: unknown effect target '{other}'")),
        };
        let op = compile_op(
            names,
            target_kind,
            &e.set,
            &e.shift,
            &e.map,
            &e.missing,
            true,
            ctx,
        )?;
        out.push(Effect { target, op });
    }
    Ok(out)
}

fn names_kind(names: &Names, id: u16) -> String {
    names
        .kind_ids
        .iter()
        .find(|(_, &v)| v == id)
        .map(|(k, _)| k.clone())
        .unwrap_or_else(|| format!("<kind {id}>"))
}
