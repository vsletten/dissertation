//! Dense site-state ids and state sets.
//!
//! A state is the pair (occupant species, configuration) a site is in.
//! The deck compiler interns every declared pair to a `StateId`; the
//! occupant/config split and all names live on the deck side. Keeping the
//! runtime id opaque is the defect/substitution design-ahead (design doc
//! §3.2, §6): the engine never needs to know whether two states share an
//! occupant.

/// Dense id for one (occupant, configuration) pair. u16: no real deck has
/// 65k distinct local states (kaolinite has ~40).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct StateId(pub u16);

/// A set of states, as a bitset over `StateId`. Guards and selectors are
/// membership tests on these; building them once at compile time keeps the
/// hot loop free of hashing and string work.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StateSet {
    words: Vec<u64>,
}

impl StateSet {
    pub fn new(n_states: usize) -> Self {
        StateSet {
            words: vec![0; n_states.div_ceil(64)],
        }
    }

    pub fn insert(&mut self, s: StateId) {
        let i = s.0 as usize;
        self.words[i / 64] |= 1 << (i % 64);
    }

    pub fn contains(&self, s: StateId) -> bool {
        let i = s.0 as usize;
        self.words
            .get(i / 64)
            .is_some_and(|w| w & (1 << (i % 64)) != 0)
    }

    pub fn is_empty(&self) -> bool {
        self.words.iter().all(|w| *w == 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn insert_and_contains() {
        let mut s = StateSet::new(130);
        s.insert(StateId(0));
        s.insert(StateId(64));
        s.insert(StateId(129));
        assert!(s.contains(StateId(0)));
        assert!(s.contains(StateId(64)));
        assert!(s.contains(StateId(129)));
        assert!(!s.contains(StateId(1)));
        assert!(!s.contains(StateId(500)));
    }
}
