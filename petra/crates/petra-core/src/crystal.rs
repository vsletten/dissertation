//! The crystal template: a triclinic unit cell with sites and bonds.
//!
//! Design doc §3.1. The deck compiler builds a [`UnitCell`] with *expanded*
//! bonds (each declared bond appears on both endpoints, with mirrored cell
//! offsets); [`UnitCell::check_reciprocity`] verifies that invariant.

/// Dense id for a site kind (e.g. "Al_oct", "O_bridge_SiSi" on the deck side).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct KindId(pub u16);

/// Dense id for a bond label ("pair", …); labels are optional deck vocabulary
/// for selecting a specific bonded partner.
pub const NO_LABEL: u16 = u16::MAX;

/// Cell geometry, stored as the fractional → Cartesian matrix (columns are
/// the cell vectors). Construct from conventional parameters
/// ([`Cell::from_params`]) or directly from a matrix ([`Cell::from_matrix`])
/// — the latter admits nonstandard conventions like the legacy kaolinite
/// shear cell, which conventional (a b c α β γ) cannot express verbatim.
#[derive(Debug, Clone, Copy)]
pub struct Cell {
    m: [[f64; 3]; 3],
}

impl Cell {
    /// Standard crystallographic construction with **a** along x and **b**
    /// in the xy plane. Angles in degrees.
    pub fn from_params(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> Self {
        let (al, be, ga) = (alpha.to_radians(), beta.to_radians(), gamma.to_radians());
        let (ca, cb, cg, sg) = (al.cos(), be.cos(), ga.cos(), ga.sin());
        let cx = c * cb;
        let cy = c * (ca - cb * cg) / sg;
        let cz = (c * c - cx * cx - cy * cy).max(0.0).sqrt();
        Cell {
            m: [[a, b * cg, cx], [0.0, b * sg, cy], [0.0, 0.0, cz]],
        }
    }

    pub fn from_matrix(m: [[f64; 3]; 3]) -> Self {
        Cell { m }
    }

    pub fn matrix(&self) -> [[f64; 3]; 3] {
        self.m
    }

    /// Cartesian → fractional, via the matrix inverse. Errors on a singular
    /// (degenerate) cell.
    pub fn to_fractional(&self, cart: [f64; 3]) -> Result<[f64; 3], String> {
        let m = &self.m;
        let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        if det.abs() < 1e-12 {
            return Err("cell matrix is singular".into());
        }
        let inv = [
            [
                (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det,
                (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / det,
                (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det,
            ],
            [
                (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det,
                (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det,
                (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / det,
            ],
            [
                (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det,
                (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / det,
                (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det,
            ],
        ];
        Ok([
            inv[0][0] * cart[0] + inv[0][1] * cart[1] + inv[0][2] * cart[2],
            inv[1][0] * cart[0] + inv[1][1] * cart[1] + inv[1][2] * cart[2],
            inv[2][0] * cart[0] + inv[2][1] * cart[1] + inv[2][2] * cart[2],
        ])
    }

    /// Fractional coordinates (with integer cell offset) → Cartesian.
    pub fn to_cartesian(&self, frac: [f64; 3], cell_coord: [i32; 3]) -> [f64; 3] {
        let m = self.matrix();
        let f = [
            frac[0] + cell_coord[0] as f64,
            frac[1] + cell_coord[1] as f64,
            frac[2] + cell_coord[2] as f64,
        ];
        [
            m[0][0] * f[0] + m[0][1] * f[1] + m[0][2] * f[2],
            m[1][0] * f[0] + m[1][1] * f[1] + m[1][2] * f[2],
            m[2][0] * f[0] + m[2][1] * f[1] + m[2][2] * f[2],
        ]
    }
}

/// One directed half of a bond: from the owning site to template site `to`
/// in the cell displaced by `dcell`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TemplateBond {
    pub to: usize,
    pub dcell: [i32; 3],
    /// Interned label id, or [`NO_LABEL`].
    pub label: u16,
}

/// A site position in the unit cell.
#[derive(Debug, Clone)]
pub struct TemplateSite {
    pub kind: KindId,
    pub frac: [f64; 3],
    pub bonds: Vec<TemplateBond>,
}

/// The full template: cell geometry + sites with expanded bonds.
#[derive(Debug, Clone)]
pub struct UnitCell {
    pub cell: Cell,
    pub sites: Vec<TemplateSite>,
}

impl UnitCell {
    /// Every bond i→(j, d) must have a mirror j→(i, −d) with the same label.
    /// The deck compiler expands declared bonds to satisfy this; a violation
    /// here is a compiler bug, not a user error.
    pub fn check_reciprocity(&self) -> Result<(), String> {
        for (i, site) in self.sites.iter().enumerate() {
            for b in &site.bonds {
                let mirror = TemplateBond {
                    to: i,
                    dcell: [-b.dcell[0], -b.dcell[1], -b.dcell[2]],
                    label: b.label,
                };
                let target = self
                    .sites
                    .get(b.to)
                    .ok_or_else(|| format!("bond from site {i} to missing site {}", b.to))?;
                if !target.bonds.contains(&mirror) {
                    return Err(format!(
                        "bond {i}→({},{:?}) has no mirror on site {}",
                        b.to, b.dcell, b.to
                    ));
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orthorhombic_matrix_is_diagonal() {
        let cell = Cell::from_params(3.0, 4.0, 5.0, 90.0, 90.0, 90.0);
        let m = cell.matrix();
        assert!((m[0][0] - 3.0).abs() < 1e-12);
        assert!((m[1][1] - 4.0).abs() < 1e-12);
        assert!((m[2][2] - 5.0).abs() < 1e-12);
        assert!(m[0][1].abs() < 1e-12 && m[0][2].abs() < 1e-9 && m[1][2].abs() < 1e-9);
        let p = cell.to_cartesian([0.5, 0.5, 0.5], [1, 0, 0]);
        assert!((p[0] - 4.5).abs() < 1e-12);
        assert!((p[1] - 2.0).abs() < 1e-12);
        assert!((p[2] - 2.5).abs() < 1e-12);
    }

    #[test]
    fn reciprocity_detects_missing_mirror() {
        let good = UnitCell {
            cell: Cell::from_params(1.0, 1.0, 1.0, 90.0, 90.0, 90.0),
            sites: vec![TemplateSite {
                kind: KindId(0),
                frac: [0.0; 3],
                bonds: vec![
                    TemplateBond {
                        to: 0,
                        dcell: [1, 0, 0],
                        label: NO_LABEL,
                    },
                    TemplateBond {
                        to: 0,
                        dcell: [-1, 0, 0],
                        label: NO_LABEL,
                    },
                ],
            }],
        };
        assert!(good.check_reciprocity().is_ok());

        let mut bad = good.clone();
        bad.sites[0].bonds.pop();
        assert!(bad.check_reciprocity().is_err());
    }
}
