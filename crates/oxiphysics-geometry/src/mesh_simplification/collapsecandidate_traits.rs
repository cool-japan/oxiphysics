//! # CollapseCandidate - Trait Implementations
//!
//! This module contains trait implementations for `CollapseCandidate`.
//!
//! ## Implemented Traits
//!
//! - `PartialEq`
//! - `Eq`
//! - `PartialOrd`
//! - `Ord`
//!
//! The order is that of a min-heap on the cost (so that `BinaryHeap` pops the
//! cheapest collapse first). It is a total order: NaN costs sort last, and equal
//! costs are ordered by the vertex indices, so the collapse sequence does not
//! depend on the order in which candidates were pushed.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::CollapseCandidate;

impl PartialEq for CollapseCandidate {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == std::cmp::Ordering::Equal
    }
}

impl Eq for CollapseCandidate {}

impl PartialOrd for CollapseCandidate {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for CollapseCandidate {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other
            .cost
            .total_cmp(&self.cost)
            .then_with(|| (other.v0, other.v1).cmp(&(self.v0, self.v1)))
    }
}
