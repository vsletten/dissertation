//! The random-number seam: an [`Rng`] trait and a **bit-faithful** port of
//! Numerical Recipes' `ran2` (the legacy model's only source of randomness,
//! `ran2.cpp`).
//!
//! # Why a trait at all (design doc §5)
//!
//! The whole point of dynamics parity (spec §C2a) is to chase the C++'s
//! *exact* PRNG stream once, prove the port faithful, then be free to swap in
//! a modern generator for production ensembles without touching the engine.
//! A trait is the seam that makes "swap the generator" a one-line change:
//! [`step`](crate::step) takes `&mut impl Rng`, so [`Ran2`] (for parity) and,
//! later, a `rand::SmallRng` wrapper (for speed) are interchangeable.
//!
//! \[IDIOM\] **A trait with two intended implementors, one shipped today.**
//! The RUST_TOUR reviewer's-lens rule is "a trait with one implementor is
//! premature." This one earns its keep: the second implementor (a modern
//! PRNG) is named and dated in the design doc — the trait exists precisely so
//! parity and production don't fork the engine. Contrast the *legacy* code,
//! where `ran2()` is a free function with file-static state that nothing can
//! substitute; testing or reseeding it means editing it.
//!
//! # Why `f32`, not the design doc's `f64`
//!
//! Design doc §5 sketches `fn uniform(&mut self) -> f64`. We return **`f32`**
//! instead, deliberately, for the same reason rates stay `f32` through M6
//! (spec §C2a): the C++ `float ran2(void)` returns a 32-bit float, and every
//! downstream consumer — `eps` in the selection loop, the `r < 0.5` proton
//! coin — is a `float`. To reproduce the trajectory *bit-for-bit* the draws
//! must be the same 32-bit values. The `f64` upgrade belongs with the whole
//! `f64` dynamics switch (spec B8), which is a *reform* (M7+), announced not
//! smuggled. Faithful now, improved later, never silently.

/// A uniform `(0, 1)` random source.
///
/// One method, by design: the legacy model only ever needs a uniform draw.
/// `&mut self` because a PRNG mutates its internal state on every call — the
/// signature makes that visible, where the C++ hides it in file-static
/// variables a caller cannot see or reset.
pub trait Rng {
    /// Next uniform variate in the open interval `(0, 1)`, as `f32`.
    fn uniform(&mut self) -> f32;
}

// ran2 constants (Numerical Recipes, verbatim from ran2.cpp's #defines).
// `long` in the C++ is 64-bit on the golden run's platform; `i64` here.
const IM1: i64 = 2_147_483_563;
const IM2: i64 = 2_147_483_399;
const IMM1: i64 = IM1 - 1;
const IA1: i64 = 40_014;
const IA2: i64 = 40_692;
const IQ1: i64 = 53_668;
const IQ2: i64 = 52_774;
const IR1: i64 = 12_211;
const IR2: i64 = 3_791;
const NTAB: usize = 32;
const NDIV: i64 = 1 + IMM1 / (NTAB as i64);
/// `AM = 1.0/IM1`, computed in `f64` exactly as the C++ macro — the final
/// `AM * iy` product is a double, narrowed to `f32` at the `temp` assignment.
const AM: f64 = 1.0 / (IM1 as f64);
const EPS: f64 = 1.2e-7;
const RNMX: f64 = 1.0 - EPS;

/// Numerical Recipes `ran2`, ported literally — the long-period (>2×10¹⁸)
/// combined LCG the legacy model runs on.
///
/// State that the C++ keeps in **file-static** variables (`idum2`, `iy`,
/// `iv[]`) lives here as struct fields, which is the whole readability win:
/// the generator's state is a value you can construct, clone for a
/// side-experiment, or reset — none of which the `static long iy` version
/// permits. `idum` (the C++ passes it by pointer and mutates through it) is
/// just another field.
///
/// # The legacy seed (spec B2)
///
/// The C++ seeds from `Simulation::ranseed`, which the input bug (B2) leaves
/// at **0** every run. `ran2`'s self-seed branch maps a non-positive seed to
/// 1, so *every legacy run shares one stream*. That accidental determinism is
/// exactly what makes trajectory parity possible — [`Ran2::legacy`] builds
/// the generator in that state.
#[derive(Debug, Clone)]
pub struct Ran2 {
    idum: i64,
    idum2: i64,
    iy: i64,
    iv: [i64; NTAB],
}

impl Ran2 {
    /// Construct with an explicit seed (the value the C++ would pass to
    /// `initran2`). A seed `<= 0` triggers `ran2`'s self-seed on first draw,
    /// exactly as the original.
    pub fn new(seed: i64) -> Self {
        Ran2 {
            idum: seed,
            // C++ file-static initializers: `idum2 = 123456789`, `iy = 0`.
            idum2: 123_456_789,
            iy: 0,
            iv: [0; NTAB],
        }
    }

    /// The legacy generator: seed 0, reproducing the B2 fixed stream every
    /// run of the 2001 model ever produced. This is the constructor the
    /// parity test and `--legacy` mode use.
    pub fn legacy() -> Self {
        Ran2::new(0)
    }
}

impl Rng for Ran2 {
    /// One `ran2` draw. Transliterated statement-for-statement from
    /// `ran2.cpp` so the bit pattern of every return value matches; see the
    /// per-line notes.
    fn uniform(&mut self) -> f32 {
        // Self-seed branch: taken on the first call when `idum <= 0`.
        if self.idum <= 0 {
            // `if (-(*idum) < 1) *idum = 1; else *idum = -(*idum);`
            self.idum = if -self.idum < 1 { 1 } else { -self.idum };
            self.idum2 = self.idum;
            // Warm-up: `for (j = NTAB+7; j >= 0; j--)`. The loop shuffles the
            // table `iv` and only writes entries once `j < NTAB`.
            for j in (0..=(NTAB + 7)).rev() {
                let k = self.idum / IQ1;
                self.idum = IA1 * (self.idum - k * IQ1) - k * IR1;
                if self.idum < 0 {
                    self.idum += IM1;
                }
                if j < NTAB {
                    self.iv[j] = self.idum;
                }
            }
            self.iy = self.iv[0];
        }
        // Main step of the first LCG.
        let k = self.idum / IQ1;
        self.idum = IA1 * (self.idum - k * IQ1) - k * IR1;
        if self.idum < 0 {
            self.idum += IM1;
        }
        // Second LCG.
        let k = self.idum2 / IQ2;
        self.idum2 = IA2 * (self.idum2 - k * IQ2) - k * IR2;
        if self.idum2 < 0 {
            self.idum2 += IM2;
        }
        // Combine via the shuffle table.
        let j = (self.iy / NDIV) as usize;
        self.iy = self.iv[j] - self.idum2;
        self.iv[j] = self.idum;
        if self.iy < 1 {
            self.iy += IMM1;
        }
        // `if ((temp = AM*iy) > RNMX) return RNMX; else return temp;`
        // AM*iy is an f64 product; the C++ narrows it to the `float temp` at
        // the assignment, then compares the (re-promoted) float against the
        // f64 RNMX. We reproduce both the narrowing and the clamp.
        let temp = (AM * self.iy as f64) as f32;
        if temp as f64 > RNMX {
            RNMX as f32
        } else {
            temp
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference values captured from the C++ `ran2` (seed 0, legacy stream),
    /// compiled `g++ -O3 -ffast-math` on the golden toolchain — the raw
    /// 32-bit hex bit patterns of the first 20 draws. If a single bit of the
    /// port drifts, one of these fails. This is the RNG's golden gate.
    ///
    /// Capture harness (scratch, not committed): a `main` that does
    /// `long seed = 0; initran2(&seed); for i in 0..20 print ran2()` against
    /// the read-only `ran2.cpp`.
    const CPP_RAN2_SEED0_BITS: [u32; 20] = [
        0x3e921d72, 0x3e81b82a, 0x3dbf6c6e, 0x3f1bc674, 0x3f67468d, 0x3e4892fb, 0x3eed083f,
        0x3f7063b4, 0x3e0244d9, 0x3ed4f4ec, 0x3f08a31e, 0x3ddc0cb0, 0x3ecca1af, 0x3f267eb7,
        0x3cddc5b0, 0x3d9da4f8, 0x3d8d4886, 0x3f5a191e, 0x3f232917, 0x3f12b677,
    ];

    #[test]
    fn ran2_is_bit_faithful_to_the_cpp_legacy_stream() {
        let mut rng = Ran2::legacy();
        for (i, &want_bits) in CPP_RAN2_SEED0_BITS.iter().enumerate() {
            let got = rng.uniform();
            let got_bits = got.to_bits();
            assert_eq!(
                got_bits, want_bits,
                "draw {i}: got {got:.9} (0x{got_bits:08x}), want 0x{want_bits:08x}"
            );
        }
    }

    #[test]
    fn ran2_stays_in_the_open_unit_interval() {
        let mut rng = Ran2::legacy();
        for _ in 0..10_000 {
            let v = rng.uniform();
            assert!(v > 0.0 && v < 1.0, "ran2 escaped (0,1): {v}");
        }
    }

    #[test]
    fn legacy_is_seed_zero() {
        // Pin the B2 tie: legacy() == new(0), the fixed-seed stream.
        let mut a = Ran2::legacy();
        let mut b = Ran2::new(0);
        for _ in 0..100 {
            assert_eq!(a.uniform().to_bits(), b.uniform().to_bits());
        }
    }
}
