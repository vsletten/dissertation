//! Lattice instantiation: tile the unit cell, resolve bonds across cell
//! boundaries per the boundary conditions, store adjacency CSR-style.
//!
//! Design doc §3.1. Variable coordination is structural here — a site has
//! exactly the neighbors its bonds resolve to, no fixed-width array, no
//! sentinel indices (the legacy `nbr[6]` out-of-bounds phantom, kmc-rs
//! reform R1, is unrepresentable).

use crate::crystal::UnitCell;
use crate::state::StateId;

pub type SiteId = usize;

/// Per-axis boundary condition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Boundary {
    /// Wrap around (bulk direction).
    Periodic,
    /// Bonds crossing the face are dropped (free surface on both faces).
    Open,
    /// Like `Open`, but sites in the layer of cells at coordinate 0 along
    /// this axis are frozen: they appear in neighbors' environments but
    /// never react. Anchors a slab from below.
    Fixed,
}

/// The instantiated site graph plus per-site metadata and dynamic state.
#[derive(Debug, Clone)]
pub struct Lattice {
    pub dims: [usize; 3],
    pub boundary: [Boundary; 3],
    /// Sites per unit cell (template length).
    pub n_template: usize,
    /// Which template site each lattice site instantiates.
    pub template_index: Vec<u16>,
    /// Frozen sites (from `Boundary::Fixed`) never host events.
    pub frozen: Vec<bool>,
    /// CSR adjacency: neighbors of site `s` are `adj[adj_off[s]..adj_off[s+1]]`.
    pub adj_off: Vec<u32>,
    pub adj: Vec<u32>,
    /// Bond label per adjacency entry (parallel to `adj`); `NO_LABEL` if none.
    pub adj_label: Vec<u16>,
    /// Dynamic per-site state, initialized by the caller.
    pub states: Vec<StateId>,
}

impl Lattice {
    /// Flat index of (cell a,b,c, template site t): fastest-varying is t.
    #[inline]
    pub fn index(&self, cell: [usize; 3], t: usize) -> SiteId {
        ((cell[0] * self.dims[1] + cell[1]) * self.dims[2] + cell[2]) * self.n_template + t
    }

    /// Inverse of [`Lattice::index`].
    #[inline]
    pub fn coords(&self, s: SiteId) -> ([usize; 3], usize) {
        let t = s % self.n_template;
        let cell = s / self.n_template;
        let c = cell % self.dims[2];
        let b = (cell / self.dims[2]) % self.dims[1];
        let a = cell / (self.dims[2] * self.dims[1]);
        ([a, b, c], t)
    }

    pub fn len(&self) -> usize {
        self.states.len()
    }

    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }

    #[inline]
    pub fn neighbors(&self, s: SiteId) -> &[u32] {
        &self.adj[self.adj_off[s] as usize..self.adj_off[s + 1] as usize]
    }

    #[inline]
    pub fn neighbor_labels(&self, s: SiteId) -> &[u16] {
        &self.adj_label[self.adj_off[s] as usize..self.adj_off[s + 1] as usize]
    }

    /// Tile `ucell` into `dims` cells with the given boundary conditions.
    /// `initial_state(template_index)` supplies each site's starting state;
    /// richer fills (substitution fractions, explicit defect lists — design
    /// doc §6) layer on top by mutating `states` afterwards.
    pub fn build(
        ucell: &UnitCell,
        dims: [usize; 3],
        boundary: [Boundary; 3],
        initial_state: impl Fn(usize) -> StateId,
    ) -> Self {
        let n_template = ucell.sites.len();
        let n = dims[0] * dims[1] * dims[2] * n_template;

        let mut lat = Lattice {
            dims,
            boundary,
            n_template,
            template_index: Vec::with_capacity(n),
            frozen: Vec::with_capacity(n),
            adj_off: Vec::with_capacity(n + 1),
            adj: Vec::new(),
            adj_label: Vec::new(),
            states: Vec::with_capacity(n),
        };

        lat.adj_off.push(0);
        for a in 0..dims[0] {
            for b in 0..dims[1] {
                for c in 0..dims[2] {
                    for (t, tsite) in ucell.sites.iter().enumerate() {
                        lat.template_index.push(t as u16);
                        lat.states.push(initial_state(t));
                        let frozen = (0..3).any(|ax| {
                            boundary[ax] == Boundary::Fixed && [a, b, c][ax] == 0
                        });
                        lat.frozen.push(frozen);

                        for bond in &tsite.bonds {
                            let mut target = [0usize; 3];
                            let mut in_range = true;
                            for ax in 0..3 {
                                let raw = [a, b, c][ax] as i64 + bond.dcell[ax] as i64;
                                let d = dims[ax] as i64;
                                match boundary[ax] {
                                    Boundary::Periodic => {
                                        target[ax] = raw.rem_euclid(d) as usize;
                                    }
                                    Boundary::Open | Boundary::Fixed => {
                                        if raw < 0 || raw >= d {
                                            in_range = false;
                                            break;
                                        }
                                        target[ax] = raw as usize;
                                    }
                                }
                            }
                            if in_range {
                                let to = ((target[0] * dims[1] + target[1]) * dims[2]
                                    + target[2])
                                    * n_template
                                    + bond.to;
                                lat.adj.push(to as u32);
                                lat.adj_label.push(bond.label);
                            }
                        }
                        lat.adj_off.push(lat.adj.len() as u32);
                    }
                }
            }
        }
        lat
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::crystal::{Cell, KindId, TemplateBond, TemplateSite, NO_LABEL};

    /// Simple-cubic single-site cell bonded to ±a, ±b, ±c.
    fn sc_cell() -> UnitCell {
        let mut bonds = Vec::new();
        for ax in 0..3usize {
            for sign in [1i32, -1] {
                let mut d = [0i32; 3];
                d[ax] = sign;
                bonds.push(TemplateBond {
                    to: 0,
                    dcell: d,
                    label: NO_LABEL,
                });
            }
        }
        UnitCell {
            cell: Cell {
                a: 3.0,
                b: 3.0,
                c: 3.0,
                alpha: 90.0,
                beta: 90.0,
                gamma: 90.0,
            },
            sites: vec![TemplateSite {
                kind: KindId(0),
                frac: [0.0; 3],
                bonds,
            }],
        }
    }

    #[test]
    fn periodic_cube_every_site_has_six_neighbors() {
        let lat = Lattice::build(
            &sc_cell(),
            [3, 3, 3],
            [Boundary::Periodic; 3],
            |_| StateId(0),
        );
        assert_eq!(lat.len(), 27);
        for s in 0..lat.len() {
            assert_eq!(lat.neighbors(s).len(), 6, "site {s}");
            assert!(!lat.frozen[s]);
        }
        // Reciprocity of the instantiated graph.
        for s in 0..lat.len() {
            for &n in lat.neighbors(s) {
                assert!(lat.neighbors(n as usize).contains(&(s as u32)));
            }
        }
    }

    #[test]
    fn open_axis_drops_face_bonds() {
        let lat = Lattice::build(
            &sc_cell(),
            [3, 3, 4],
            [Boundary::Periodic, Boundary::Periodic, Boundary::Open],
            |_| StateId(0),
        );
        for s in 0..lat.len() {
            let ([_, _, c], _) = lat.coords(s);
            let expect = if c == 0 || c == 3 { 5 } else { 6 };
            assert_eq!(lat.neighbors(s).len(), expect, "site {s} at c={c}");
        }
    }

    #[test]
    fn fixed_axis_freezes_bottom_layer_only() {
        let lat = Lattice::build(
            &sc_cell(),
            [2, 2, 3],
            [Boundary::Periodic, Boundary::Periodic, Boundary::Fixed],
            |_| StateId(0),
        );
        for s in 0..lat.len() {
            let ([_, _, c], _) = lat.coords(s);
            assert_eq!(lat.frozen[s], c == 0, "site {s} at c={c}");
        }
    }
}
