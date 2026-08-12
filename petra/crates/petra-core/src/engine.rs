//! The rejection-free KMC engine (design doc §5): two-level selection over
//! a Fenwick tree of per-site rate sums, incremental event maintenance via
//! the deck's maximum read distance, f64 time, seedable PCG64 RNG.

use rand::{Rng as _, SeedableRng};
use rand_pcg::Pcg64Mcg;

use crate::crystal::KindId;
use crate::lattice::{Lattice, SiteId};
use crate::reaction::{
    first_match, guards_pass, resolve_rate, sites_at_distance, EffectTarget, Reaction,
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
}

/// Steps between full tree rebuilds (float-drift hygiene).
const REBUILD_EVERY: u64 = 1 << 16;

pub struct Engine {
    pub lattice: Lattice,
    pub reactions: Vec<Reaction>,
    /// Per-site kind, resolved once from the template.
    kinds: Vec<KindId>,
    /// Reaction ids grouped by center kind.
    by_kind: Vec<Vec<u16>>,
    temperature: f64,
    /// Cached candidate events per site: (reaction id, resolved rate).
    site_events: Vec<Vec<(u16, f64)>>,
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
}

impl Engine {
    pub fn new(
        lattice: Lattice,
        kinds_per_template: &[KindId],
        n_kinds: usize,
        reactions: Vec<Reaction>,
        temperature: f64,
        seed: u64,
    ) -> Self {
        let kinds: Vec<KindId> = lattice
            .template_index
            .iter()
            .map(|&t| kinds_per_template[t as usize])
            .collect();
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
        let mut engine = Engine {
            lattice,
            reactions,
            kinds,
            by_kind,
            temperature,
            site_events: vec![Vec::new(); n],
            tree: RateTree::new(n),
            max_read,
            time: 0.0,
            step_count: 0,
            rng: Pcg64Mcg::seed_from_u64(seed),
            scratch: Vec::new(),
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
        if !self.lattice.frozen[s] {
            let kind = self.kinds[s];
            let state = self.lattice.states[s];
            for &ri in &self.by_kind[kind.0 as usize] {
                let rxn = &self.reactions[ri as usize];
                if rxn.center_states.contains(state)
                    && guards_pass(&self.lattice, &self.kinds, rxn, s, &mut self.scratch)
                {
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
        let total: f64 = events.iter().map(|&(_, r)| r).sum();
        self.tree.set(s, total);
        self.site_events[s] = events;
    }

    /// One KMC step: select, apply, propagate, advance time.
    pub fn step(&mut self) -> Result<Fired, Stop> {
        let total = self.tree.total();
        if total <= 0.0 {
            let any = self.site_events.iter().any(|e| !e.is_empty());
            return Err(if any { Stop::ZeroRate } else { Stop::NoEvents });
        }

        // Site, then event within site. A `None` here means the positive
        // `total` was pure accumulated drift over zero leaves — no events.
        let Some((site, mut residual)) = self.tree.find(self.rng.gen::<f64>() * total) else {
            let any = self.site_events.iter().any(|e| !e.is_empty());
            return Err(if any { Stop::ZeroRate } else { Stop::NoEvents });
        };
        let events = &self.site_events[site];
        debug_assert!(!events.is_empty(), "tree selected an event-less site");
        let mut chosen = events.len() - 1; // clamp to last on float edge
        for (i, &(_, rate)) in events.iter().enumerate() {
            if residual < rate {
                chosen = i;
                break;
            }
            residual -= rate;
        }
        let (ri, _) = events[chosen];

        // Poisson waiting time; map u∈[0,1) to (0,1] so ln never sees 0.
        let u: f64 = self.rng.gen();
        let dt = -(1.0 - u).ln() / total;

        let changed = self.apply(site, ri)?;

        // Dirty propagation: changed sites plus everything within max_read.
        let mut dirty: Vec<SiteId> = changed.clone();
        let mut ring = Vec::new();
        for &c in &changed {
            for d in 1..=self.max_read {
                sites_at_distance(&self.lattice, c, d, &mut ring);
                for &s in &ring {
                    if !dirty.contains(&s) {
                        dirty.push(s);
                    }
                }
            }
        }
        for s in dirty {
            self.refresh_site(s);
        }

        self.time += dt;
        self.step_count += 1;
        if self.step_count % REBUILD_EVERY == 0 {
            self.tree.rebuild();
        }
        Ok(Fired {
            step: self.step_count,
            time: self.time,
            site,
            reaction: ri,
        })
    }

    /// Apply the chosen reaction's effects (choosing a branch if several),
    /// returning the sites whose state changed.
    fn apply(&mut self, site: SiteId, ri: u16) -> Result<Vec<SiteId>, Stop> {
        let rxn = &self.reactions[ri as usize];
        let branch = if rxn.branches.len() == 1 {
            &rxn.branches[0]
        } else {
            let wsum: f64 = rxn.branches.iter().map(|b| b.weight).sum();
            let mut draw = self.rng.gen::<f64>() * wsum;
            let mut pick = rxn.branches.len() - 1;
            for (i, b) in rxn.branches.iter().enumerate() {
                if draw < b.weight {
                    pick = i;
                    break;
                }
                draw -= b.weight;
            }
            &rxn.branches[pick]
        };

        let mut changed = Vec::with_capacity(branch.effects.len());
        // Resolve all targets against the *pre-effect* state, then write —
        // an effect must not see a sibling effect's result.
        let mut writes: Vec<(SiteId, crate::state::StateId)> =
            Vec::with_capacity(branch.effects.len());
        for eff in &branch.effects {
            let target = match &eff.target {
                EffectTarget::Center => site,
                EffectTarget::FirstMatch(sel) => first_match(
                    &self.lattice,
                    &self.kinds,
                    site,
                    sel,
                    &mut self.scratch,
                )
                .ok_or(Stop::EffectTargetMissing { site, reaction: ri })?,
            };
            writes.push((target, eff.set));
        }
        for (target, new_state) in writes {
            if self.lattice.states[target] != new_state {
                self.lattice.states[target] = new_state;
                changed.push(target);
            }
        }
        if changed.is_empty() {
            // Self-transition: nothing to dirty beyond the center itself
            // (its own rate may depend on its state — refresh regardless).
            changed.push(site);
        }
        Ok(changed)
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
        let probes = [0.0, 0.4999, 0.5, 1.0, 2.4999, 2.5, 3.7, 5.7499, 5.75, 6.9999];
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
            assert!(residual >= 0.0 && residual <= rates[idx] + 1e-12, "probe {p}");
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
