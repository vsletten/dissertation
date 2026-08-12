//! Compiled reactions: guarded local rewrites with environment-modified TST
//! rates (design doc §3.3–§4). Everything here is dense ids and bitsets —
//! the deck compiler has already resolved all names.

use crate::crystal::KindId;
use crate::lattice::{Lattice, SiteId};
use crate::rate::{RateExpr, R_KCAL};
use crate::state::{StateId, StateSet};

/// Selects neighbors of a center site: sites at exact graph `distance`
/// (1 or 2), optionally restricted by kind, bond label (distance 1 only),
/// and a state set.
#[derive(Debug, Clone)]
pub struct NeighborSelect {
    pub distance: u8,
    pub kind: Option<KindId>,
    pub label: Option<u16>,
    pub states: StateSet,
}

/// "Between `min` and `max` neighbors match the selector."
#[derive(Debug, Clone)]
pub struct Guard {
    pub select: NeighborSelect,
    pub min: u32,
    pub max: u32,
}

/// Environment-dependent rate adjustment (design doc §4).
#[derive(Debug, Clone)]
pub enum ModifierKind {
    /// Each matching neighbor adds `dea` (kcal/mol) to the activation
    /// energy — bond-counting / BEP-style. The *linear* convenience case
    /// of [`ModifierKind::ByCount`].
    PerMatch { dea: f64 },
    /// Tabulated ΔEa by match count: `dea[n]` is added when `n` neighbors
    /// match, with the last entry extending to all higher counts. The
    /// general nonlinear form — real barriers are rarely linear in
    /// coordination, so a measured/computed table goes here verbatim.
    /// Non-empty by construction (deck validation).
    ByCount { dea: Vec<f64> },
    /// If the match count lies in `[min, max]`, add `dea` and multiply the
    /// rate by `factor` — discrete overrides for non-additive cases.
    When {
        min: u32,
        max: u32,
        dea: f64,
        factor: f64,
    },
}

#[derive(Debug, Clone)]
pub struct Modifier {
    pub select: NeighborSelect,
    pub kind: ModifierKind,
}

/// Where an effect lands.
#[derive(Debug, Clone)]
pub enum EffectTarget {
    Center,
    /// The first neighbor (ascending site id) matching the selector.
    /// v0 semantics; a random-match option is planned (design doc §3.3).
    FirstMatch(NeighborSelect),
}

/// One state assignment.
#[derive(Debug, Clone)]
pub struct Effect {
    pub target: EffectTarget,
    pub set: StateId,
}

/// A weighted alternative outcome (generalizes the legacy R4/R9 proton
/// coin flip). A deterministic reaction has exactly one branch.
#[derive(Debug, Clone)]
pub struct Branch {
    pub weight: f64,
    pub effects: Vec<Effect>,
}

/// A fully compiled elementary reaction.
#[derive(Debug, Clone)]
pub struct Reaction {
    pub name: String,
    pub center_kind: KindId,
    pub center_states: StateSet,
    pub guards: Vec<Guard>,
    pub rate: RateExpr,
    /// ln of the solution-coupling factor (activities, chemical potentials),
    /// folded once at compile time: k *= exp(ln_thermo).
    pub ln_thermo: f64,
    pub modifiers: Vec<Modifier>,
    pub branches: Vec<Branch>,
}

impl Reaction {
    /// Largest guard/modifier/effect-selector distance this reaction reads.
    pub fn max_read_distance(&self) -> u8 {
        let g = self.guards.iter().map(|g| g.select.distance).max();
        let m = self.modifiers.iter().map(|m| m.select.distance).max();
        let e = self
            .branches
            .iter()
            .flat_map(|b| &b.effects)
            .filter_map(|e| match &e.target {
                EffectTarget::Center => None,
                EffectTarget::FirstMatch(s) => Some(s.distance),
            })
            .max();
        g.max(m).max(e).unwrap_or(0)
    }
}

/// Collect the unique sites at exact graph distance 1 or 2 from `center`
/// (distance 2 excludes the center and all distance-1 sites). Coordination
/// numbers are small, so linear membership checks beat hashing here.
pub fn sites_at_distance(lat: &Lattice, center: SiteId, distance: u8, out: &mut Vec<SiteId>) {
    out.clear();
    match distance {
        1 => out.extend(lat.neighbors(center).iter().map(|&n| n as SiteId)),
        2 => {
            let first: Vec<SiteId> = lat.neighbors(center).iter().map(|&n| n as SiteId).collect();
            for &n in &first {
                for &nn in lat.neighbors(n) {
                    let nn = nn as SiteId;
                    if nn != center && !first.contains(&nn) && !out.contains(&nn) {
                        out.push(nn);
                    }
                }
            }
        }
        d => panic!("unsupported guard distance {d} (deck validation admits 1 or 2)"),
    }
}

/// Count neighbors of `center` matching `sel`. `scratch` is reused across
/// calls to avoid per-evaluation allocation.
pub fn count_matches(
    lat: &Lattice,
    kinds: &[KindId],
    center: SiteId,
    sel: &NeighborSelect,
    scratch: &mut Vec<SiteId>,
) -> u32 {
    if sel.distance == 1 && sel.label.is_some() {
        // Label filters only make sense on direct bonds; walk the CSR row
        // with its parallel label array.
        let label = sel.label.unwrap();
        let nbrs = lat.neighbors(center);
        let labels = lat.neighbor_labels(center);
        return nbrs
            .iter()
            .zip(labels)
            .filter(|&(&n, &l)| {
                l == label && site_matches(lat, kinds, n as SiteId, sel)
            })
            .count() as u32;
    }
    sites_at_distance(lat, center, sel.distance, scratch);
    scratch
        .iter()
        .filter(|&&s| site_matches(lat, kinds, s, sel))
        .count() as u32
}

/// First matching neighbor for an effect target, if any.
pub fn first_match(
    lat: &Lattice,
    kinds: &[KindId],
    center: SiteId,
    sel: &NeighborSelect,
    scratch: &mut Vec<SiteId>,
) -> Option<SiteId> {
    if sel.distance == 1 && sel.label.is_some() {
        let label = sel.label.unwrap();
        let nbrs = lat.neighbors(center);
        let labels = lat.neighbor_labels(center);
        return nbrs
            .iter()
            .zip(labels)
            .find(|&(&n, &l)| l == label && site_matches(lat, kinds, n as SiteId, sel))
            .map(|(&n, _)| n as SiteId);
    }
    sites_at_distance(lat, center, sel.distance, scratch);
    scratch
        .iter()
        .copied()
        .find(|&s| site_matches(lat, kinds, s, sel))
}

#[inline]
fn site_matches(lat: &Lattice, kinds: &[KindId], s: SiteId, sel: &NeighborSelect) -> bool {
    if let Some(k) = sel.kind {
        if kinds[s] != k {
            return false;
        }
    }
    sel.states.contains(lat.states[s])
}

/// Resolve the rate of `rxn` at `center` in its current environment.
pub fn resolve_rate(
    lat: &Lattice,
    kinds: &[KindId],
    rxn: &Reaction,
    center: SiteId,
    temperature: f64,
    scratch: &mut Vec<SiteId>,
) -> f64 {
    let rt = R_KCAL * temperature;
    let mut extra_ea = 0.0;
    let mut factor = 1.0;
    for m in &rxn.modifiers {
        let n = count_matches(lat, kinds, center, &m.select, scratch);
        match &m.kind {
            ModifierKind::PerMatch { dea } => extra_ea += dea * n as f64,
            ModifierKind::ByCount { dea } => {
                extra_ea += dea[(n as usize).min(dea.len() - 1)];
            }
            ModifierKind::When {
                min,
                max,
                dea,
                factor: f,
            } => {
                if n >= *min && n <= *max {
                    extra_ea += dea;
                    factor *= f;
                }
            }
        }
    }
    rxn.rate.base_rate(temperature) * rxn.ln_thermo.exp() * (-extra_ea / rt).exp() * factor
}

/// Do all guards of `rxn` pass at `center`? (Center kind/state are checked
/// by the caller against its per-kind tables.)
pub fn guards_pass(
    lat: &Lattice,
    kinds: &[KindId],
    rxn: &Reaction,
    center: SiteId,
    scratch: &mut Vec<SiteId>,
) -> bool {
    rxn.guards.iter().all(|g| {
        let n = count_matches(lat, kinds, center, &g.select, scratch);
        n >= g.min && n <= g.max
    })
}
