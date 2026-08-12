//! TST rate expressions (design doc §4). All energies in kcal/mol,
//! temperatures in Kelvin, f64 throughout — the f32 accumulation-order
//! hazards of the legacy code (kmc-rs spec §C2a) are out of scope by
//! construction.

/// Gas constant, kcal·mol⁻¹·K⁻¹ (the dissertation's 1.987e-3 convention).
pub const R_KCAL: f64 = 1.987e-3;

/// k_B / h in s⁻¹·K⁻¹, for the Eyring prefactor.
pub const KB_OVER_H: f64 = 2.083_661_2e10;

/// The base rate law of one reaction, before environment modifiers and
/// solution-coupling factors.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RateExpr {
    /// A lumped/experimental rate constant, temperature-independent.
    Constant { k: f64 },
    /// k = A·exp(−Ea/RT).
    Arrhenius { prefactor: f64, ea: f64 },
    /// Full TST: k = (k_B·T/h)·exp(ΔS‡/R)·exp(−ΔH‡/RT).
    /// `ds` in kcal·mol⁻¹·K⁻¹, `dh` in kcal/mol.
    Eyring { dh: f64, ds: f64 },
}

impl RateExpr {
    pub fn base_rate(&self, temperature: f64) -> f64 {
        let rt = R_KCAL * temperature;
        match *self {
            RateExpr::Constant { k } => k,
            RateExpr::Arrhenius { prefactor, ea } => prefactor * (-ea / rt).exp(),
            RateExpr::Eyring { dh, ds } => {
                KB_OVER_H * temperature * (ds / R_KCAL).exp() * (-dh / rt).exp()
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arrhenius_hand_value() {
        // A=1e13, Ea=10 kcal/mol, T=298.15 K → RT=0.59245...,
        // k = 1e13 * exp(-10/0.5924...) = 1e13 * exp(-16.879...)
        let k = RateExpr::Arrhenius {
            prefactor: 1e13,
            ea: 10.0,
        }
        .base_rate(298.15);
        let rt = R_KCAL * 298.15;
        let expect = 1e13 * (-10.0 / rt).exp();
        assert!((k - expect).abs() / expect < 1e-12);
        assert!(k > 4.0e5 && k < 5.0e5, "k = {k}"); // ballpark sanity
    }

    #[test]
    fn eyring_zero_barrier_is_kbt_over_h() {
        let k = RateExpr::Eyring { dh: 0.0, ds: 0.0 }.base_rate(300.0);
        assert!((k - KB_OVER_H * 300.0).abs() / k < 1e-12);
    }

    #[test]
    fn constant_ignores_temperature() {
        let e = RateExpr::Constant { k: 2.5 };
        assert_eq!(e.base_rate(100.0), 2.5);
        assert_eq!(e.base_rate(1000.0), 2.5);
    }
}
