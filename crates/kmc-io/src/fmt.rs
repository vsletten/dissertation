//! Reproducing C++ `ostream << float` formatting, digit for digit.
//!
//! This tiny module carries a disproportionate share of the M3 bitwise
//! gate. The C++ writers print coordinates with plain `fxyz << x` — i.e.
//! `std::ostream`'s **defaultfloat, precision 6**: up to six *significant*
//! digits, trailing zeros stripped, scientific notation when the exponent
//! leaves [-4, 6). That is exactly C `printf("%g", v)` — and libstdc++
//! implements it by promoting the `float` to `double` and handing it to the
//! `%g` machinery. So the contract to match is: **`printf("%g", (double)f)`
//! for every f32 the model prints.**
//!
//! [`format_g6`] implements the C17 `%g` algorithm verbatim (7.21.6.1):
//! with precision P = 6, take the `%e` exponent X; if −4 ≤ X < P, format
//! as `%f` with precision P−1−X, else as `%e` with precision P−1; strip
//! trailing zeros either way. Both glibc's printf and Rust's `{:.*}` /
//! `{:.*e}` produce correctly rounded digits (glibc via exact multi-
//! precision, Rust via Dragon4 in its precision paths), so building `%g`
//! from Rust's primitives yields byte-identical output — verified against
//! all 3,000 coordinates of the golden `start.msi` by the M3 gate, and
//! spot-pinned in the tests below.

/// Format an `f32` exactly as the C++ model's `ofstream << x` does.
///
/// Promotes to `f64` first — that is what `ostream::operator<<(float)`
/// does before formatting, and skipping the promotion would round twice.
pub fn format_g6(v: f32) -> String {
    let d = v as f64;
    if d == 0.0 {
        // %g prints a bare "0" (and "-0" for negative zero — never seen in
        // the golden outputs, but cheap to get right).
        return if d.is_sign_negative() {
            "-0".to_string()
        } else {
            "0".to_string()
        };
    }
    // Step 1: the %e exponent AFTER rounding to 6 significant digits.
    // Deriving X from the already-rounded string (not from log10) is what
    // makes the 9.999996 → 1e+01 carry cases come out right.
    let sci = format!("{:.5e}", d);
    // [IDIOM] `expect` with a message vs bare `unwrap`: both panic, but
    // `expect` states the invariant being relied on ("{:.5e} always
    // contains 'e'"), which is precisely what a reviewer needs to judge
    // whether the panic is truly unreachable. House rule worth adopting:
    // unwrap is for tests; production panics must carry their proof.
    let epos = sci
        .find('e')
        .expect("{:.5e} output always contains an exponent");
    let exp: i32 = sci[epos + 1..]
        .parse()
        .expect("exponent is a valid integer");

    // [IDIOM] `!(-4..6).contains(&exp)` — clippy nudges C-style compound
    // comparisons toward the Range type, whose bounds read unambiguously
    // (half-open, like all Rust ranges). Same machine code either way.
    if !(-4..6).contains(&exp) {
        // %e style: mantissa with 5 decimals, zeros stripped, exponent
        // re-dressed in C's clothing (e+NN, at least two digits — Rust
        // writes `1.23456e7`, C writes `1.23456e+07`).
        let mut mant = sci[..epos].to_string();
        if mant.contains('.') {
            while mant.ends_with('0') {
                mant.pop();
            }
            if mant.ends_with('.') {
                mant.pop();
            }
        }
        format!(
            "{}e{}{:02}",
            mant,
            if exp < 0 { '-' } else { '+' },
            exp.abs()
        )
    } else {
        // %f style with precision 5 - X, zeros stripped.
        let prec = (5 - exp) as usize;
        let mut s = format!("{:.*}", prec, d);
        if s.contains('.') {
            while s.ends_with('0') {
                s.pop();
            }
            if s.ends_with('.') {
                s.pop();
            }
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_golden_msi_samples() {
        // Values lifted from data/golden/outputs/start.msi — each already
        // round-tripped through the C++ f32 pipeline and %g.
        assert_eq!(format_g6(1.17181), "1.17181");
        assert_eq!(format_g6(4.38134), "4.38134");
        assert_eq!(format_g6(1.9999), "1.9999"); // trailing-zero strip
        assert_eq!(format_g6(0.225645), "0.225645");
        assert_eq!(format_g6(-0.417473), "-0.417473");
    }

    #[test]
    fn integers_lose_the_point() {
        assert_eq!(format_g6(100.0), "100");
        assert_eq!(format_g6(-3.0), "-3");
        assert_eq!(format_g6(0.0), "0");
        assert_eq!(format_g6(-0.0), "-0");
    }

    #[test]
    fn scientific_kicks_in_outside_the_g_window() {
        assert_eq!(format_g6(0.0000123), "1.23e-05");
        assert_eq!(format_g6(1234567.0), "1.23457e+06");
        assert_eq!(format_g6(-0.00001), "-1e-05");
        // Six digits exactly at the boundary stays fixed-notation.
        assert_eq!(format_g6(999999.0), "999999");
    }

    #[test]
    fn rounding_carry_promotes_the_exponent() {
        // 999999.9f32 rounds to 1.00000e6 at 6 sig figs → "1e+06", the
        // classic %g carry case the two-step algorithm exists for.
        assert_eq!(format_g6(999999.9), "1e+06");
    }
}
