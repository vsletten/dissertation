//! `data.rxn` reader — temperature, chemical potentials, and the reaction
//! rate tables. Port of `ReactionList::CreateReactionList` (rxnlist.cpp).
//!
//! Rate construction (spec A5.4), reproduced with the C++'s exact
//! float/double dance — every narrowing documented at the line it happens:
//!
//! * hydrolysis forward: `rate = k` (the pre-exponential as given; the
//!   Boltzmann prefactor was considered and abandoned upstream — the
//!   commented-out `kt` machinery, spec B10)
//! * hydrolysis reverse: `rate = k · exp(ΔE / rt)`
//! * adsorption: `rate = a · exp(ΔE/rt) · exp(Δμ/rt)`, where reactions
//!   below index 18 take `dm_si` and the rest `dm_al` — ported by index
//!   test, exactly as written, even though index 16 is *Al* adsorption
//!   (harmless today: the golden file has dm_si == dm_al)
//! * desorption: `rate = a · exp(−ΔE / rt)`
//! * diffusion: not ported (dead + out-of-bounds read; spec B6/B7 — see
//!   `kaolinite::reactions`)

use std::path::Path;

use kaolinite::reactions::{N_ADS, N_DES, N_HYD, R_KCAL, Reaction, ReactionSet};
use kaolinite::state::State;

use crate::error::ReadError;
use crate::scan::Scanner;

/// Read `data.rxn` into a [`ReactionSet`].
pub fn read_rxn(path: &Path) -> Result<ReactionSet, ReadError> {
    let mut s = Scanner::open(path)?;
    parse_rxn(&mut s)
}

fn parse_rxn(s: &mut Scanner) -> Result<ReactionSet, ReadError> {
    let temperature = s.next_f32("invalid temperature")?;
    s.eat_comment();
    let dm_si = s.next_f32("invalid si chemical potential")?;
    s.eat_comment();
    let dm_al = s.next_f32("invalid al chemical potential")?;
    s.eat_comment();

    // C++: `rt = R * t` with float rt, double R → the product is computed
    // in f64 and narrowed to f32. Every exp() below divides by THIS f32.
    let rt = (R_KCAL * temperature as f64) as f32;

    let mut reactions: Vec<Reaction> = Vec::with_capacity(N_DES);

    // Hydrolysis: 8 forward/reverse pairs sharing one (k, dE) table.
    // C++ fills reactions[i] and reactions[i+1] in the same file pass.
    for _pair in 0..(N_HYD / 2) {
        let reactant = s.next_i32("invalid hydrolysis reactant")?;
        let product = s.next_i32("invalid hydrolysis product")?;
        s.eat_comment();
        let num_rates = s.next_i32("invalid hydrolysis rate count")?;
        s.eat_comment();
        let mut fwd = Vec::with_capacity(num_rates.max(0) as usize);
        let mut rev = Vec::with_capacity(num_rates.max(0) as usize);
        for _ in 0..num_rates {
            let kp = s.next_f32("invalid hydrolysis pre-exponential")?;
            let de = s.next_f32("invalid hydrolysis activation energy")?;
            // Forward: the pre-exponential verbatim (no Arrhenius factor —
            // faithful to A5.4, surprising as it reads).
            fwd.push(kp);
            // Reverse: kp * exp(de / rt). C++ evaluation order: de/rt in
            // f32, promoted to f64 for exp (the C library's exp(double) —
            // <cmath> unqualified call on a float promotes), multiplied by
            // kp promoted to f64, narrowed once at the assignment. Rust's
            // f64::exp calls the same platform libm, so the bits match.
            rev.push(((kp as f64) * (((de / rt) as f64).exp())) as f32);
        }
        reactions.push(Reaction {
            reactant: State(reactant),
            info: 0,
            rates: fwd,
        });
        reactions.push(Reaction {
            reactant: State(product),
            info: 0,
            rates: rev,
        });
    }

    // Adsorption: one rate each, coupled to solution saturation.
    s.eat_comment();
    for i in N_HYD..N_ADS {
        let reactant = s.next_i32("invalid adsorption reactant")?;
        let product = s.next_i32("invalid adsorption product")?;
        s.eat_comment();
        let a = s.next_f32("invalid adsorption pre-exponential")?;
        let de = s.next_f32("invalid adsorption energy")?;
        s.eat_comment();
        // WART-adjacent, ported as-is: the C++ selects the potential by
        // `i < 18`, which routes reaction 16 (adsorb *Al*) through dm_si.
        // Same f32-divide-then-f64-exp chain as above, two factors.
        let dm = if i < 18 { dm_si } else { dm_al };
        let rate = ((a as f64) * (((de / rt) as f64).exp()) * (((dm / rt) as f64).exp())) as f32;
        reactions.push(Reaction {
            reactant: State(reactant),
            info: product,
            rates: vec![rate],
        });
    }

    // Desorption: per-bucket tables, rate = a * exp(-de / rt).
    s.eat_comment();
    for _i in N_ADS..N_DES {
        let reactant = s.next_i32("invalid desorption reactant")?;
        s.eat_comment();
        let num_rates = s.next_i32("invalid desorption rate count")?;
        let mut rates = Vec::with_capacity(num_rates.max(0) as usize);
        for _ in 0..num_rates {
            let a = s.next_f32("invalid desorption pre-exponential")?;
            let de = s.next_f32("invalid desorption energy")?;
            rates.push(((a as f64) * ((((-de) / rt) as f64).exp())) as f32);
        }
        reactions.push(Reaction {
            reactant: State(reactant),
            info: 0,
            rates,
        });
    }

    // Diffusion (24–27): NOT ported — dead code with a buffer over-read
    // (spec B6/B7). The C++ would copy each desorption table here, reading
    // one element past the end of the 4-entry tables because `numRates`
    // still holds 5 from the last desorption block.

    Ok(ReactionSet {
        temperature,
        dm_si,
        dm_al,
        reactions,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn golden() -> ReactionSet {
        let path = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../data/golden/inputs/data.rxn"
        );
        read_rxn(Path::new(path)).unwrap()
    }

    #[test]
    fn golden_conditions_parse() {
        let r = golden();
        assert_eq!(r.temperature, 8000.0);
        assert_eq!(r.dm_si, -1.0);
        assert_eq!(r.dm_al, -1.0);
        assert_eq!(r.reactions.len(), N_DES); // 24: no diffusion tail
    }

    #[test]
    fn golden_reaction_shapes_match_the_file() {
        let r = golden();
        // R0/R1: si-o-si hydrolysis, 15 buckets, reactant 301 fwd / 302 rev.
        assert_eq!(r.reactions[0].reactant, State(301));
        assert_eq!(r.reactions[1].reactant, State(302));
        assert_eq!(r.reactions[0].rates.len(), 15);
        // The file's quirky "40100" product token for R13 (reads like a
        // typo for 405; parsed as-is, it just makes R13's reactant match
        // no site ever — data wart, preserved by parsing not fixing).
        assert_eq!(r.reactions[13].reactant, State(40100));
        // Adsorption block: reactants and products per the file comments.
        assert_eq!(r.reactions[16].reactant, State(100)); // adsorb Al
        assert_eq!(r.reactions[16].info, 107);
        assert_eq!(r.reactions[19].reactant, State(200)); // adsorb Si
        assert_eq!(r.reactions[19].info, 205);
        // Desorption tables: 4,4,5,5 buckets.
        let lens: Vec<usize> = (20..24).map(|i| r.reactions[i].rates.len()).collect();
        assert_eq!(lens, [4, 4, 5, 5]);
    }

    #[test]
    fn golden_rates_follow_the_cpp_float_paths() {
        let r = golden();
        let rt = (R_KCAL * 8000.0f64) as f32;
        // Forward hydrolysis: pre-exponential verbatim.
        assert_eq!(r.reactions[0].rates[0], 1.0);
        // Reverse: kp * exp(de/rt) with de/rt narrowed to f32 first.
        let expect = ((1.0f64) * (((2.6f32 / rt) as f64).exp())) as f32;
        assert_eq!(r.reactions[1].rates[0], expect);
        // All buckets within a reaction are identical in the shipped file
        // (spec B4 — the environment machinery is inert in-sample).
        assert!(r.reactions[1].rates.iter().all(|&x| x == expect));
        // Desorption: a * exp(-de/rt).
        let expect = ((1998.0f64) * ((((-12.0f32) / rt) as f64).exp())) as f32;
        assert_eq!(r.reactions[20].rates[1], expect);
    }
}
