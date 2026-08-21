//! The rejection-free KMC engine (design doc §5): two-level selection over
//! a Fenwick tree of per-site rate sums, incremental event maintenance via
//! the deck's maximum read distance, f64 time, seedable PCG64 RNG.

use rand::{Rng as _, RngCore, SeedableRng};
use rand_pcg::Pcg64Mcg;

use crate::crystal::KindId;
use crate::lattice::{Lattice, SiteId};
use crate::rate::R_KCAL;
use crate::reaction::{
    all_matches, first_match, guards_pass, resolve_energy_delta, resolve_rate, sites_at_distance,
    EffectTarget, Reaction, RuleValue,
};

/// Fenwick (binary-indexed) tree over per-site total rates: O(log N) point
/// update, prefix search, and total. Leaf values are kept alongside so
/// updates are set-by-difference and periodic full rebuilds (float-drift
/// hygiene, design doc §5.1) are cheap.
#[derive(Debug, Clone)]
pub struct RateTree {
    tree: Vec<f64>, // 1-based
    leaf: Vec<f64>,
    n: usize,
}

impl RateTree {
    pub fn new(n: usize) -> Self {
        RateTree {
            tree: vec![0.0; n + 1],
            leaf: vec![0.0; n],
            n,
        }
    }

    pub fn set(&mut self, i: usize, value: f64) {
        let delta = value - self.leaf[i];
        self.leaf[i] = value;
        let mut j = i + 1;
        while j <= self.n {
            self.tree[j] += delta;
            j += j & j.wrapping_neg();
        }
    }

    pub fn get(&self, i: usize) -> f64 {
        self.leaf[i]
    }

    pub fn total(&self) -> f64 {
        self.prefix(self.n)
    }

    fn prefix(&self, mut i: usize) -> f64 {
        let mut s = 0.0;
        while i > 0 {
            s += self.tree[i];
            i -= i & i.wrapping_neg();
        }
        s
    }

    /// Smallest index whose prefix sum exceeds `target`; also returns the
    /// residual target within that leaf (for second-level selection).
    /// `target` must be in `[0, total())`.
    ///
    /// Float care: in exact arithmetic the descent can only land on a
    /// nonzero leaf, but the descent path and `prefix()` sum the same nodes
    /// in different orders, so at a boundary adjacent to zero-rate leaves an
    /// ulp of disagreement can strand the result on a zero leaf. Zero leaves
    /// carry no probability mass, so we snap to the nearest nonzero leaf
    /// (forward first — the direction exact arithmetic would have taken).
    /// Returns `None` when no leaf is positive at all — possible when
    /// internal-node drift leaves `total()` slightly above zero with every
    /// leaf actually zero; the caller must treat that as "no events".
    pub fn find(&self, mut target: f64) -> Option<(usize, f64)> {
        let mut pos = 0usize;
        let mut mask = self.n.next_power_of_two();
        while mask > 0 {
            let next = pos + mask;
            if next <= self.n && self.tree[next] <= target {
                target -= self.tree[next];
                pos = next;
            }
            mask >>= 1;
        }
        // pos = count of leaves fully below target; clamp for float edge.
        let mut pos = pos.min(self.n - 1);
        if self.leaf[pos] == 0.0 {
            pos = (pos + 1..self.n)
                .find(|&i| self.leaf[i] > 0.0)
                .or_else(|| (0..pos).rev().find(|&i| self.leaf[i] > 0.0))?;
            target = 0.0;
        }
        Some((pos, target.clamp(0.0, self.leaf[pos])))
    }

    /// Recompute internal sums from leaves (drift hygiene).
    pub fn rebuild(&mut self) {
        for v in &mut self.tree {
            *v = 0.0;
        }
        for i in 0..self.n {
            let val = self.leaf[i];
            self.leaf[i] = 0.0;
            self.set(i, val);
        }
    }
}

/// One fired event.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Fired {
    pub step: u64,
    pub time: f64,
    pub site: SiteId,
    pub reaction: u16,
}

/// Why the simulation cannot advance.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum Stop {
    #[error("no site has any possible event")]
    NoEvents,
    #[error("events exist but all rates are zero")]
    ZeroRate,
    #[error("reaction {reaction} at site {site}: effect selector matched no neighbor")]
    EffectTargetMissing { site: SiteId, reaction: u16 },
    #[error("reaction {reaction} at site {site}: {reason}")]
    EffectFailed {
        site: SiteId,
        reaction: u16,
        reason: &'static str,
    },
}

/// Steps between full tree rebuilds (float-drift hygiene).
const REBUILD_EVERY: u64 = 1 << 16;

/// The result of one strategy step. CTMC fires exactly one transition;
/// future synchronous strategies may return several.
#[derive(Debug, Clone, PartialEq)]
pub struct StepOutcome {
    pub fired: Vec<Fired>,
    pub dt: f64,
}

#[derive(Debug)]
struct PreparedTransition {
    site: SiteId,
    reaction: u16,
    writes: Vec<(SiteId, crate::state::StateId)>,
}

/// Core-owned transition capability exposed to update strategies. Strategies
/// may inspect the CTMC selection tables and schedule transitions, but only the
/// engine commits writes, records changes, and refreshes dirty propensities.
pub struct ApplyHandle<'a> {
    lattice: &'a Lattice,
    rules: &'a [Reaction],
    kinds: &'a [KindId],
    kind_state_ranges: &'a [(u16, u16)],
    site_events: &'a [Vec<(u16, f64)>],
    enabled_rules: &'a [Vec<u16>],
    live_sites: &'a [SiteId],
    tree: &'a RateTree,
    scratch: &'a mut Vec<SiteId>,
    pending: &'a mut Vec<PreparedTransition>,
}

impl ApplyHandle<'_> {
    /// Total enabled CTMC propensity.
    pub fn total_rate(&self) -> f64 {
        self.tree.total()
    }

    /// Whether any site has an enabled event, including zero-rate events.
    pub fn has_events(&self) -> bool {
        self.site_events.iter().any(|events| !events.is_empty())
    }

    /// Enabled rules at `site`, in rule-index order, evaluated against the
    /// immutable pre-step state.
    pub fn enabled_rules(&self, site: SiteId) -> &[u16] {
        &self.enabled_rules[site]
    }

    pub fn energy_delta(&mut self, site: SiteId, reaction: u16) -> f64 {
        resolve_energy_delta(
            self.lattice,
            self.kinds,
            &self.rules[reaction as usize],
            site,
            self.scratch,
        )
    }

    pub fn probability(&self, reaction: u16) -> f64 {
        let RuleValue::Probability(probability) = self.rules[reaction as usize].value else {
            panic!("probability called for a non-PCA rule")
        };
        probability
    }

    pub fn live_sites(&self) -> &[SiteId] {
        self.live_sites
    }

    /// Select one enabled event using the legacy coupled site/event draw.
    /// Returning `None` means the positive tree total was only float drift.
    pub fn select_event(&self, draw: f64) -> Option<(SiteId, u16)> {
        let (site, mut residual) = self.tree.find(draw)?;
        let events = &self.site_events[site];
        debug_assert!(!events.is_empty(), "tree selected an event-less site");
        let mut chosen = events.len() - 1;
        for (i, &(_, rate)) in events.iter().enumerate() {
            if residual < rate {
                chosen = i;
                break;
            }
            residual -= rate;
        }
        Some((site, events[chosen].0))
    }

    /// Resolve and queue a transition against the current pre-step state.
    /// Branch and effect-selection randomness comes only from `rng`.
    pub fn apply_transition(
        &mut self,
        site: SiteId,
        reaction: u16,
        rng: &mut dyn RngCore,
    ) -> Result<(), Stop> {
        let rxn = &self.rules[reaction as usize];
        let branch = if rxn.branches.len() == 1 {
            &rxn.branches[0]
        } else {
            let wsum: f64 = rxn.branches.iter().map(|branch| branch.weight).sum();
            let mut draw = rng.gen::<f64>() * wsum;
            let mut pick = rxn.branches.len() - 1;
            for (index, branch) in rxn.branches.iter().enumerate() {
                if draw < branch.weight {
                    pick = index;
                    break;
                }
                draw -= branch.weight;
            }
            &rxn.branches[pick]
        };

        let mut targets: Vec<(SiteId, &crate::reaction::EffectOp)> =
            Vec::with_capacity(branch.effects.len());
        let mut matched = Vec::new();
        for effect in &branch.effects {
            match &effect.target {
                EffectTarget::Center | EffectTarget::Source => targets.push((site, &effect.op)),
                EffectTarget::FirstMatch(selector) => {
                    let target =
                        first_match(self.lattice, self.kinds, site, selector, self.scratch)
                            .ok_or(Stop::EffectTargetMissing { site, reaction })?;
                    targets.push((target, &effect.op));
                }
                EffectTarget::RandomMatch(selector) => {
                    all_matches(
                        self.lattice,
                        self.kinds,
                        site,
                        selector,
                        self.scratch,
                        &mut matched,
                    );
                    if matched.is_empty() {
                        return Err(Stop::EffectTargetMissing { site, reaction });
                    }
                    matched.sort_unstable();
                    let draw = rng.gen::<f64>();
                    let index = ((draw * matched.len() as f64) as usize).min(matched.len() - 1);
                    targets.push((matched[index], &effect.op));
                }
                EffectTarget::AllMatches(selector) => {
                    all_matches(
                        self.lattice,
                        self.kinds,
                        site,
                        selector,
                        self.scratch,
                        &mut matched,
                    );
                    targets.extend(matched.iter().map(|&target| (target, &effect.op)));
                }
            }
        }

        let mut writes = Vec::with_capacity(targets.len());
        for (target, op) in targets {
            if self.lattice.frozen[target] {
                return Err(Stop::EffectFailed {
                    site,
                    reaction,
                    reason: "effect writes a frozen site",
                });
            }
            let range = self.kind_state_ranges[self.kinds[target].0 as usize];
            match op.resolve(self.lattice.states[target], range) {
                Ok(Some(new_state)) => writes.push((target, new_state)),
                Ok(None) => {}
                Err(reason) => {
                    return Err(Stop::EffectFailed {
                        site,
                        reaction,
                        reason,
                    });
                }
            }
        }
        self.pending.push(PreparedTransition {
            site,
            reaction,
            writes,
        });
        Ok(())
    }
}

/// What the core hands a strategy each step, matching RFC-001 §3. State and
/// rules are read-only; all randomness and mutation use the public seams.
pub struct StepCtx<'a> {
    pub lattice: &'a Lattice,
    pub rules: &'a [Reaction],
    pub rng: &'a mut dyn RngCore,
    pub apply: ApplyHandle<'a>,
}

/// A strategy decides which transition fires and how simulation time advances.
pub trait UpdateStrategy {
    fn step(&mut self, ctx: &mut StepCtx<'_>) -> Result<StepOutcome, Stop>;
}

/// Exact rejection-free CTMC: the pre-RFC Petra algorithm, moved intact
/// behind [`UpdateStrategy`].
#[derive(Debug, Default)]
pub struct ExactCtmc;

/// Deterministic, double-buffered cellular automaton. Every rule match reads
/// the pre-step lattice; at most one rule fires per site, with the first
/// enabled rule in rule-index (deck declaration) order taking priority. The
/// engine commits all selected writes as one batch.
#[derive(Debug, Default)]
pub struct SynchronousCA;

/// Asynchronous single-site Metropolis updates. Temperature is Kelvin and
/// rule energies are canonical kcal/mol, so acceptance uses exp(-ΔE/RT).
#[derive(Debug)]
pub struct AsyncMetropolis {
    temperature: f64,
}

impl AsyncMetropolis {
    pub fn new(temperature: f64) -> Self {
        assert!(temperature.is_finite() && temperature > 0.0);
        Self { temperature }
    }
}

/// Per-site probabilistic cellular automaton with Bernoulli decisions made
/// against one shared pre-step state and committed as one batch.
#[derive(Debug, Default)]
pub struct DiscreteTimePCA;

/// A deck-selected strategy value suitable for library and CLI callers.
#[derive(Debug)]
pub enum Strategy {
    ExactCtmc(ExactCtmc),
    Synchronous(SynchronousCA),
    Metropolis(AsyncMetropolis),
    Pca(DiscreteTimePCA),
}

pub struct Engine {
    pub lattice: Lattice,
    pub reactions: Vec<Reaction>,
    /// Per-site kind, resolved once from the template.
    kinds: Vec<KindId>,
    /// Reaction ids grouped by center kind.
    by_kind: Vec<Vec<u16>>,
    /// Contiguous StateId range (start, count) per kind, for Shift effects.
    kind_state_ranges: Vec<(u16, u16)>,
    temperature: f64,
    /// Cached candidate events per site: (reaction id, resolved rate).
    site_events: Vec<Vec<(u16, f64)>>,
    /// Center/guard-enabled rules in rule-index order, independent of rate.
    enabled_rules: Vec<Vec<u16>>,
    /// Live-site ids in ascending order, fixed for the lifetime of a lattice.
    live_sites: Vec<SiteId>,
    tree: RateTree,
    /// Global maximum read distance over all reactions: after an event
    /// changes some sites, every site within this graph distance of a
    /// changed site must re-enumerate (design doc §5.2; per-reaction radii
    /// are a planned refinement — global max is conservative and correct).
    max_read: u8,
    pub time: f64,
    pub step_count: u64,
    rng: Pcg64Mcg,
    scratch: Vec<SiteId>,
    /// (site, old state, new state) for every site the LAST applied event
    /// actually changed — the trajectory-export/streaming feed. Reused
    /// across steps; read it before calling [`Engine::step`] again.
    last_changes: Vec<(SiteId, crate::state::StateId, crate::state::StateId)>,
}

impl Engine {
    pub fn new(
        lattice: Lattice,
        kinds_per_template: &[KindId],
        n_kinds: usize,
        kind_state_ranges: Vec<(u16, u16)>,
        reactions: Vec<Reaction>,
        temperature: f64,
        seed: u64,
    ) -> Self {
        let kinds = lattice
            .template_index
            .iter()
            .map(|&template| kinds_per_template[template as usize])
            .collect();
        Self::new_with_site_kinds(
            lattice,
            kinds,
            n_kinds,
            kind_state_ranges,
            reactions,
            temperature,
            seed,
        )
    }

    /// Construct an engine with explicit per-site kinds after init passes.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_site_kinds(
        lattice: Lattice,
        kinds: Vec<KindId>,
        n_kinds: usize,
        kind_state_ranges: Vec<(u16, u16)>,
        reactions: Vec<Reaction>,
        temperature: f64,
        seed: u64,
    ) -> Self {
        assert_eq!(lattice.len(), kinds.len(), "one kind per lattice site");
        let mut by_kind = vec![Vec::new(); n_kinds];
        for (i, r) in reactions.iter().enumerate() {
            by_kind[r.center_kind.0 as usize].push(i as u16);
        }
        let max_read = reactions
            .iter()
            .map(|r| r.max_read_distance())
            .max()
            .unwrap_or(0);
        let n = lattice.len();
        let live_sites = lattice
            .frozen
            .iter()
            .enumerate()
            .filter_map(|(site, &frozen)| (!frozen).then_some(site))
            .collect();
        let mut engine = Engine {
            lattice,
            reactions,
            kinds,
            by_kind,
            kind_state_ranges,
            temperature,
            site_events: vec![Vec::new(); n],
            enabled_rules: vec![Vec::new(); n],
            live_sites,
            tree: RateTree::new(n),
            max_read,
            time: 0.0,
            step_count: 0,
            rng: Pcg64Mcg::seed_from_u64(seed),
            scratch: Vec::new(),
            last_changes: Vec::new(),
        };
        for s in 0..n {
            engine.refresh_site(s);
        }
        engine
    }

    /// Re-enumerate the candidate events at one site and push its new total
    /// into the tree.
    fn refresh_site(&mut self, s: SiteId) {
        let mut events = std::mem::take(&mut self.site_events[s]);
        events.clear();
        let mut enabled = std::mem::take(&mut self.enabled_rules[s]);
        enabled.clear();
        if !self.lattice.frozen[s] {
            let kind = self.kinds[s];
            let state = self.lattice.states[s];
            for &ri in &self.by_kind[kind.0 as usize] {
                let rxn = &self.reactions[ri as usize];
                if rxn.center_states.contains(state)
                    && guards_pass(&self.lattice, &self.kinds, rxn, s, &mut self.scratch)
                {
                    enabled.push(ri);
                    if rxn.value == RuleValue::Ctmc {
                        let rate = resolve_rate(
                            &self.lattice,
                            &self.kinds,
                            rxn,
                            s,
                            self.temperature,
                            &mut self.scratch,
                        );
                        if rate > 0.0 {
                            events.push((ri, rate));
                        }
                    }
                }
            }
        }
        let total: f64 = events.iter().map(|&(_, r)| r).sum();
        self.tree.set(s, total);
        self.site_events[s] = events;
        self.enabled_rules[s] = enabled;
    }

    /// Advance with an explicitly supplied strategy. The strategy can only
    /// schedule writes through [`ApplyHandle`]; the core commits and records
    /// them after the strategy returns.
    pub fn step_with(&mut self, strategy: &mut impl UpdateStrategy) -> Result<StepOutcome, Stop> {
        // Guard against impossible zero totals caused by Fenwick cancellation:
        // if there are positive-rate events but the tree reports zero, rebuild
        // from the authoritative leaves before handing control to the strategy.
        if self.tree.total() <= 0.0 && self.site_events.iter().any(|e| !e.is_empty()) {
            self.tree.rebuild();
        }
        let mut pending = Vec::new();
        let result = {
            let mut ctx = StepCtx {
                lattice: &self.lattice,
                rules: &self.reactions,
                rng: &mut self.rng,
                apply: ApplyHandle {
                    lattice: &self.lattice,
                    rules: &self.reactions,
                    kinds: &self.kinds,
                    kind_state_ranges: &self.kind_state_ranges,
                    site_events: &self.site_events,
                    enabled_rules: &self.enabled_rules,
                    live_sites: &self.live_sites,
                    tree: &self.tree,
                    scratch: &mut self.scratch,
                    pending: &mut pending,
                },
            };
            strategy.step(&mut ctx)
        };

        let mut outcome = match result {
            Ok(outcome) => outcome,
            Err(stop) => {
                debug_assert!(pending.is_empty(), "strategy queued writes before stopping");
                return Err(stop);
            }
        };
        debug_assert_eq!(outcome.fired.len(), pending.len());
        for (fired, transition) in outcome.fired.iter().zip(&pending) {
            debug_assert_eq!(fired.site, transition.site);
            debug_assert_eq!(fired.reaction, transition.reaction);
        }
        if pending.len() > 1 {
            // A batch is evaluated against one pre-step state. Conflicting
            // writes cannot be made simultaneous by choosing an arbitrary
            // iteration winner, so fail atomically before mutating anything.
            let mut writes = vec![None; self.lattice.len()];
            for transition in &pending {
                for &(target, state) in &transition.writes {
                    if writes[target].is_some_and(|prior| prior != state) {
                        return Err(Stop::EffectFailed {
                            site: transition.site,
                            reaction: transition.reaction,
                            reason: "conflicting simultaneous writes",
                        });
                    }
                    writes[target] = Some(state);
                }
            }
        }
        self.commit_transitions(pending);
        self.time += outcome.dt;
        self.step_count += 1;
        if self.step_count.is_multiple_of(REBUILD_EVERY) {
            self.tree.rebuild();
        }
        for fired in &mut outcome.fired {
            fired.step = self.step_count;
            fired.time = self.time;
        }
        Ok(outcome)
    }

    /// Compatibility wrapper for the default exact-CTMC strategy.
    pub fn step(&mut self) -> Result<Fired, Stop> {
        let mut strategy = ExactCtmc;
        let outcome = self.step_with(&mut strategy)?;
        debug_assert_eq!(outcome.fired.len(), 1);
        outcome.fired.into_iter().next().ok_or(Stop::NoEvents)
    }

    /// Commit a strategy step's prepared transitions, record actual state
    /// changes, and refresh the union of their dirty neighborhoods once.
    fn commit_transitions(&mut self, transitions: Vec<PreparedTransition>) {
        self.last_changes.clear();
        // Marking array gives O(1) dedup while preserving insertion order, so
        // refresh order (and thus RNG consumption) stays deterministic.
        let mut seen = vec![false; self.lattice.len()];
        let mut changed = Vec::new();
        for transition in transitions {
            let before = changed.len();
            for (target, new_state) in transition.writes {
                let old = self.lattice.states[target];
                if old != new_state {
                    self.lattice.states[target] = new_state;
                    if !seen[target] {
                        seen[target] = true;
                        changed.push(target);
                    }
                    self.last_changes.push((target, old, new_state));
                }
            }
            if changed.len() == before && !seen[transition.site] {
                seen[transition.site] = true;
                changed.push(transition.site);
            }
        }

        let mut dirty = changed.clone();
        let mut ring = Vec::new();
        for &site in &changed {
            for distance in 1..=self.max_read {
                sites_at_distance(&self.lattice, site, distance, None, &mut ring);
                for &neighbor in &ring {
                    if !seen[neighbor] {
                        seen[neighbor] = true;
                        dirty.push(neighbor);
                    }
                }
            }
        }
        for site in dirty {
            self.refresh_site(site);
        }
    }

    /// The per-site state changes of the most recently applied event:
    /// `(site, old, new)`, actual changes only. Valid until the next
    /// [`Engine::step`].
    pub fn last_changes(&self) -> &[(SiteId, crate::state::StateId, crate::state::StateId)] {
        &self.last_changes
    }

    /// Population count per state id (scan; used at report cadence only).
    pub fn state_counts(&self, n_states: usize) -> Vec<u64> {
        let mut counts = vec![0u64; n_states];
        for &s in &self.lattice.states {
            counts[s.0 as usize] += 1;
        }
        counts
    }

    /// Differential-testing oracle (design doc §5.2 `--paranoid`):
    /// re-enumerate every site from scratch and compare totals with the
    /// incrementally maintained tree.
    pub fn paranoid_check(&mut self) -> Result<(), String> {
        for s in 0..self.lattice.len() {
            let before = self.tree.get(s);
            let events_before = self.site_events[s].clone();
            self.refresh_site(s);
            let after = self.tree.get(s);
            if (before - after).abs() > 1e-9 * after.abs().max(1.0)
                || events_before != self.site_events[s]
            {
                return Err(format!(
                    "site {s}: incremental total {before} != fresh {after}"
                ));
            }
        }
        Ok(())
    }
}

impl UpdateStrategy for ExactCtmc {
    fn step(&mut self, ctx: &mut StepCtx<'_>) -> Result<StepOutcome, Stop> {
        let total = ctx.apply.total_rate();
        if total <= 0.0 {
            return Err(if ctx.apply.has_events() {
                Stop::ZeroRate
            } else {
                Stop::NoEvents
            });
        }

        // Preserve the pre-refactor coupled site/event draw exactly.
        let draw = ctx.rng.gen::<f64>() * total;
        let Some((site, reaction)) = ctx.apply.select_event(draw) else {
            return Err(if ctx.apply.has_events() {
                Stop::ZeroRate
            } else {
                Stop::NoEvents
            });
        };

        // Poisson waiting time; map u∈[0,1) to (0,1] so ln never sees 0.
        let u: f64 = ctx.rng.gen();
        let dt = -(1.0 - u).ln() / total;
        ctx.apply.apply_transition(site, reaction, &mut *ctx.rng)?;

        Ok(StepOutcome {
            fired: vec![Fired {
                step: 0,
                time: 0.0,
                site,
                reaction,
            }],
            dt,
        })
    }
}

impl UpdateStrategy for SynchronousCA {
    fn step(&mut self, ctx: &mut StepCtx<'_>) -> Result<StepOutcome, Stop> {
        let mut selected = Vec::new();
        for site in 0..ctx.lattice.len() {
            // Same-site conflicts are intentionally draw-free: the first
            // enabled rule in deck declaration order wins (RFC-001 §3).
            if let Some(&reaction) = ctx.apply.enabled_rules(site).first() {
                selected.push((site, reaction));
            }
        }

        let mut fired = Vec::with_capacity(selected.len());
        for (site, reaction) in selected {
            ctx.apply.apply_transition(site, reaction, &mut *ctx.rng)?;
            fired.push(Fired {
                step: 0,
                time: 0.0,
                site,
                reaction,
            });
        }
        Ok(StepOutcome { fired, dt: 1.0 })
    }
}

impl UpdateStrategy for AsyncMetropolis {
    fn step(&mut self, ctx: &mut StepCtx<'_>) -> Result<StepOutcome, Stop> {
        let live_sites = ctx.apply.live_sites();
        if live_sites.is_empty() {
            return Err(Stop::NoEvents);
        }
        let site_draw = ctx.rng.gen::<f64>();
        let site =
            live_sites[((site_draw * live_sites.len() as f64) as usize).min(live_sites.len() - 1)];
        let enabled = ctx.apply.enabled_rules(site);
        let dt = 1.0 / live_sites.len() as f64;
        if enabled.is_empty() {
            return Ok(StepOutcome {
                fired: Vec::new(),
                dt,
            });
        }

        let proposal_draw = ctx.rng.gen::<f64>();
        let reaction =
            enabled[((proposal_draw * enabled.len() as f64) as usize).min(enabled.len() - 1)];
        let delta = ctx.apply.energy_delta(site, reaction);
        let acceptance = (-delta / (R_KCAL * self.temperature)).exp().min(1.0);
        let acceptance_draw = ctx.rng.gen::<f64>();
        if acceptance_draw >= acceptance {
            return Ok(StepOutcome {
                fired: Vec::new(),
                dt,
            });
        }

        ctx.apply.apply_transition(site, reaction, &mut *ctx.rng)?;
        Ok(StepOutcome {
            fired: vec![Fired {
                step: 0,
                time: 0.0,
                site,
                reaction,
            }],
            dt,
        })
    }
}

impl UpdateStrategy for DiscreteTimePCA {
    fn step(&mut self, ctx: &mut StepCtx<'_>) -> Result<StepOutcome, Stop> {
        let mut selected = Vec::new();
        for site in 0..ctx.lattice.len() {
            for &reaction in ctx.apply.enabled_rules(site) {
                let draw = ctx.rng.gen::<f64>();
                if draw < ctx.apply.probability(reaction) {
                    selected.push((site, reaction));
                }
            }
        }

        let mut fired = Vec::with_capacity(selected.len());
        for (site, reaction) in selected {
            ctx.apply.apply_transition(site, reaction, &mut *ctx.rng)?;
            fired.push(Fired {
                step: 0,
                time: 0.0,
                site,
                reaction,
            });
        }
        Ok(StepOutcome { fired, dt: 1.0 })
    }
}

impl UpdateStrategy for Strategy {
    fn step(&mut self, ctx: &mut StepCtx<'_>) -> Result<StepOutcome, Stop> {
        match self {
            Strategy::ExactCtmc(strategy) => strategy.step(ctx),
            Strategy::Synchronous(strategy) => strategy.step(ctx),
            Strategy::Metropolis(strategy) => strategy.step(ctx),
            Strategy::Pca(strategy) => strategy.step(ctx),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rate_tree_matches_linear_scan() {
        let mut tree = RateTree::new(7);
        let rates = [0.5, 0.0, 2.0, 1.25, 0.0, 3.0, 0.25];
        for (i, &r) in rates.iter().enumerate() {
            tree.set(i, r);
        }
        let total: f64 = rates.iter().sum();
        assert!((tree.total() - total).abs() < 1e-12);

        // Every probe must land on a nonzero-rate leaf whose prefix window
        // contains the target, matching a linear scan.
        let probes = [
            0.0, 0.4999, 0.5, 1.0, 2.4999, 2.5, 3.7, 5.7499, 5.75, 6.9999,
        ];
        for &p in &probes {
            let (idx, residual) = tree.find(p).expect("nonzero leaves exist");
            // linear reference
            let mut acc = 0.0;
            let mut want = rates.len() - 1;
            for (i, &r) in rates.iter().enumerate() {
                if p < acc + r {
                    want = i;
                    break;
                }
                acc += r;
            }
            assert_eq!(idx, want, "probe {p}");
            assert!(
                residual >= 0.0 && residual <= rates[idx] + 1e-12,
                "probe {p}"
            );
        }
    }

    /// Hammer `find` with sparse trees (many zero leaves) and boundary
    /// probes — the failure mode is landing on a zero-rate leaf.
    #[test]
    fn rate_tree_find_never_lands_on_zero_leaf() {
        // Deterministic pseudo-random workload, no external RNG needed.
        let mut x: u64 = 0x9e3779b97f4a7c15;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        for n in [1usize, 2, 3, 7, 8, 64, 100] {
            let mut tree = RateTree::new(n);
            let mut rates = vec![0.0f64; n];
            for _ in 0..200 {
                // Random sparse update.
                let i = (next() % n as u64) as usize;
                let r = if next() % 3 == 0 {
                    0.0
                } else {
                    (next() % 1000) as f64 / 7.0
                };
                rates[i] = r;
                tree.set(i, r);
                let total = tree.total();
                if total <= 0.0 {
                    continue;
                }
                // Probe interior points and the near-total boundary.
                for k in 0..8u64 {
                    let target = match k {
                        7 => total * (1.0 - 1e-16),
                        _ => total * (next() % 10_000) as f64 / 10_000.0,
                    };
                    match tree.find(target.min(total * 0.999_999_999_999_999)) {
                        Some((idx, residual)) => {
                            assert!(rates[idx] > 0.0, "n={n}: landed on zero leaf {idx}");
                            assert!(residual >= 0.0 && residual <= rates[idx]);
                        }
                        // Legal only when the positive total is pure drift
                        // over an all-zero leaf array.
                        None => assert!(rates.iter().all(|&r| r == 0.0)),
                    }
                }
            }
        }
    }

    #[test]
    fn rate_tree_update_and_rebuild() {
        let mut tree = RateTree::new(5);
        for i in 0..5 {
            tree.set(i, i as f64);
        }
        tree.set(2, 10.0);
        assert!((tree.total() - (0.0 + 1.0 + 10.0 + 3.0 + 4.0)).abs() < 1e-12);
        tree.rebuild();
        assert!((tree.total() - 18.0).abs() < 1e-12);
        assert_eq!(tree.get(2), 10.0);
    }
}
