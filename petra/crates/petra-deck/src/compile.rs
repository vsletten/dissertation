//! Deck → runtime compilation: intern names to dense ids, expand bonds,
//! build state sets, fold thermodynamic factors, and validate everything
//! validatable before the engine ever runs (design doc §3, §4).

use std::collections::BTreeMap;

use petra_core::crystal::{Cell, KindId, TemplateBond, TemplateSite, UnitCell, NO_LABEL};
use petra_core::lattice::{Boundary, Lattice};
use petra_core::rate::{RateExpr, R_KCAL};
use petra_core::reaction::{
    Branch, Effect, EffectTarget, Guard, Modifier, ModifierKind, NeighborSelect, Reaction,
};
use petra_core::state::{StateId, StateSet};
use petra_core::Engine;

use crate::schema::{DeckFile, EffectSpec, RateSpec, SelectorSpec};

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
    pub n_states: usize,
    pub initial_per_template: Vec<StateId>,
    pub reactions: Vec<Reaction>,
    pub temperature: f64,
    pub dims: [usize; 3],
    pub boundary: [Boundary; 3],
    pub steps: u64,
    pub seed: u64,
    pub report_every: u64,
}

impl CompiledDeck {
    /// Instantiate the lattice and engine (perfect-crystal fill; defect
    /// fills mutate `engine.lattice.states` before stepping).
    pub fn build_engine(&self, seed_override: Option<u64>) -> Engine {
        let initial = self.initial_per_template.clone();
        let lattice = Lattice::build(&self.unit_cell, self.dims, self.boundary, |t| initial[t]);
        Engine::new(
            lattice,
            &self.kinds_per_template,
            self.kind_names.len(),
            self.reactions.clone(),
            self.temperature,
            seed_override.unwrap_or(self.seed),
        )
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

pub fn compile(deck: &DeckFile) -> Result<CompiledDeck, CompileError> {
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
    let mut kind_names = Vec::new();
    let mut initial_by_kind = Vec::new();

    for (ki, k) in deck.kinds.iter().enumerate() {
        let ki = ki as u16;
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
                return err(format!("duplicate state '{}' in kind '{}'", st.name, k.name));
            }
            state_names.push(format!("{}.{}", k.name, st.name));
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

    // --- unit cell: sites, then bonds expanded onto both endpoints ---
    let cell = Cell {
        a: deck.cell.a,
        b: deck.cell.b,
        c: deck.cell.c,
        alpha: deck.cell.alpha,
        beta: deck.cell.beta,
        gamma: deck.cell.gamma,
    };
    let mut tsites = Vec::new();
    let mut kinds_per_template = Vec::new();
    let mut initial_per_template = Vec::new();
    for s in &deck.cell.sites {
        let ki = *names
            .kind_ids
            .get(&s.kind)
            .ok_or_else(|| CompileError(format!("cell site has unknown kind '{}'", s.kind)))?;
        kinds_per_template.push(KindId(ki));
        initial_per_template.push(initial_by_kind[ki as usize]);
        tsites.push(TemplateSite {
            kind: KindId(ki),
            frac: s.frac,
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
        let center_kind = *names
            .kind_ids
            .get(&r.center.kind)
            .ok_or_else(|| CompileError(format!("{ctx}: unknown center kind '{}'", r.center.kind)))?;
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
            modifiers,
            branches,
        });
    }

    let dims = deck.lattice.dims;
    if dims.iter().any(|&d| d == 0) {
        return err("lattice dims must be nonzero");
    }

    Ok(CompiledDeck {
        name: deck.deck.name.clone(),
        unit_cell,
        kinds_per_template,
        kind_names,
        state_names,
        n_states: names.n_states,
        initial_per_template,
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
    let states = names.state_set(&s.state, kind.map(|k| k.0), ctx)?;
    Ok(NeighborSelect {
        distance,
        kind,
        label,
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
                (EffectTarget::Center, center_kind)
            }
            "neighbor" => {
                let sel = e.select.as_ref().ok_or_else(|| {
                    CompileError(format!("{ctx}: a neighbor effect requires a selector"))
                })?;
                let compiled = compile_selector(names, sel, ctx)?;
                let k = compiled.kind.ok_or_else(|| {
                    CompileError(format!(
                        "{ctx}: neighbor-effect selectors must name a kind so 'set' resolves"
                    ))
                })?;
                (EffectTarget::FirstMatch(compiled), k.0)
            }
            other => return err(format!("{ctx}: unknown effect target '{other}'")),
        };
        // Resolve `set` within the target's kind.
        let set = names
            .kind_states
            .get(target_kind as usize)
            .and_then(|m| {
                m.get(&e.set).copied().or_else(|| {
                    e.set
                        .strip_prefix(&format!("{}.", names_kind(names, target_kind)))
                        .and_then(|plain| m.get(plain).copied())
                })
            })
            .ok_or_else(|| {
                CompileError(format!(
                    "{ctx}: effect sets unknown state '{}' for kind '{}'",
                    e.set,
                    names_kind(names, target_kind)
                ))
            })?;
        out.push(Effect { target, set });
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
