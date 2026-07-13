//! The site graph: a flat array of sites with fixed-degree adjacency.
//!
//! The C++ lattice (`lattice.hpp`) is `LatticeSite sites[Num_Sites]` where
//! each site holds `int nbr[6]` — neighbor *indices*, `-1` meaning "no
//! neighbor". That is already a graph in adjacency-list form with degree
//! ≤ 6; this module is the same shape with the C-isms translated:
//! sentinel indices become `Option<SiteId>`, the malloc'd array becomes an
//! owned `Vec`, and the *chemistry-specific* fields of `LatticeSite`
//! (`pair`, `lostal`, the `a/b/n` tiling coordinates, the dead BFS `color`)
//! do **not** live here — they are the model's business, kept in parallel
//! arrays on the kaolinite side (design doc §3, option (a)). What remains
//! is the part any KMC model needs: states and topology.

/// Index of a site in the flat array.
///
/// \[IDIOM\] A type alias, not a newtype — deliberately weaker than
/// `kaolinite::State`. Site ids index slices constantly (`graph.sites[id]`),
/// and a newtype would force `.0` at every subscript for little safety
/// gain *within this crate*; the alias still documents intent at API
/// boundaries. Rule of thumb worth stealing: newtype where confusion is
/// plausible and arithmetic is rare (states), alias where the value's whole
/// job is indexing (ids). Judgment call, recorded so the tour can disagree
/// with it later.
pub type SiteId = usize;

/// One site: an opaque state plus up to six neighbors.
///
/// \[IDIOM\] `Site<S>` is a *generic* struct — compile-time polymorphism, the
/// honest analogue of a C++ class template, with two differences that
/// matter to a reviewer. (1) Monomorphization is identical (zero runtime
/// cost, one copy of the code per concrete `S`), but the definition is
/// type-checked *once, here*, against the declared bounds — not re-checked
/// per instantiation at each call site, so misuse errors are small and
/// local instead of template-backtrace archaeology. (2) The absence of
/// bounds is itself information: this struct promises it will never
/// compare, copy, or print an `S` — it only stores one.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Site<S> {
    /// The model-defined state tag. The engine never interprets it.
    pub state: S,
    /// Resolved neighbor indices.
    ///
    /// \[IDIOM\] `Option<SiteId>` replaces the C++ `-1` sentinel, and the
    /// trade is worth spelling out. What it buys: "no neighbor" is now a
    /// distinct *type state* — you cannot use a missing neighbor as an
    /// index without writing the `None` case, so the whole family of
    /// "forgot the `>= 0` guard" bugs is unrepresentable. Every one of the
    /// C++ code's many `if (nbr >= 0)` checks becomes an `if let
    /// Some(nbr)` the compiler *demands*. What it costs: `Option<usize>`
    /// is 16 bytes to the C++ `int`'s 4 (usize has no spare bit-pattern
    /// for the niche optimization), so this array is 96 bytes/site instead
    /// of 24. At 1,560 sites nobody cares; if a future lattice is millions
    /// of sites, the fix is a `NonMaxUsize`-style index type, not a return
    /// to sentinels.
    pub nbr: [Option<SiteId>; 6],
}

/// The whole lattice: sites in a flat `Vec`, index = [`SiteId`].
///
/// The linear-index scheme is the model's to define (kaolinite tiles
/// `a * bCells * npos + b * npos + n`, same as the C++ so build output can
/// be compared index-for-index); the engine only requires that neighbors
/// refer to valid indices.
///
/// \[IDIOM\] Ownership as memory management: `SiteGraph` *owns* its `Vec`,
/// whoever owns the graph frees it by going out of scope, and any borrow
/// of a site pins the graph alive for exactly the borrow's duration —
/// checked at compile time. The C++ pairs `new LatticeSite[Num_Sites]`
/// with a hand-written `DisposeLattice` and a destructor that *also*
/// deletes, giving it two paths to the same free (guarded by a nulling
/// convention). None of that ceremony exists here, and none of its bug
/// surface does either.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SiteGraph<S> {
    /// The sites. Public: the model mutates states directly, and hiding a
    /// `Vec` behind delegating methods would be abstraction theater.
    pub sites: Vec<Site<S>>,
}

impl<S> SiteGraph<S> {
    /// Number of sites (C++ `GetNsites`).
    pub fn len(&self) -> usize {
        self.sites.len()
    }

    /// True if the graph has no sites. (Exists because clippy insists every
    /// `len` has an `is_empty` — a fair API-consistency lint.)
    pub fn is_empty(&self) -> bool {
        self.sites.is_empty()
    }

    /// Count a site's neighbors: the length of the leading run of `Some`
    /// entries. Port of `Lattice::CountNbrs`, **quirk included**: the C++
    /// stops at the first `-1`, so a hypothetical `[Some, None, Some, ...]`
    /// pattern would undercount. In practice templates pack real neighbors
    /// first, but the port must count the same way the original did —
    /// `TerminateLattice` (M3) freezes sites based on this number.
    ///
    /// \[IDIOM\] `take_while(...).count()` over a manual `for i in 0..6`
    /// with a break: the iterator chain *is* the sentence "count the
    /// leading Somes", and the optimizer emits the same loop. When
    /// reviewing iterator chains, the question is never performance — it's
    /// whether the chain still reads as one sentence. This one does; five
    /// combinators deep, it wouldn't.
    pub fn count_nbrs(&self, site: SiteId) -> i32 {
        self.sites[site]
            .nbr
            .iter()
            .take_while(|n| n.is_some())
            .count() as i32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(nbr: [Option<SiteId>; 6]) -> Site<i32> {
        Site { state: 0, nbr }
    }

    #[test]
    fn count_nbrs_counts_the_leading_run_only() {
        let g = SiteGraph {
            sites: vec![
                site([Some(1), Some(2), None, None, None, None]),
                // The quirk: a Some AFTER a None is not counted (C++ parity).
                site([Some(0), None, Some(2), None, None, None]),
                site([None; 6]),
            ],
        };
        assert_eq!(g.count_nbrs(0), 2);
        assert_eq!(g.count_nbrs(1), 1);
        assert_eq!(g.count_nbrs(2), 0);
    }
}
