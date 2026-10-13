//! # MeshEdge - Trait Implementations
//!
//! This module contains trait implementations for `MeshEdge`.
//!
//! ## Implemented Traits
//!
//! - `Eq`
//! - `PartialOrd`
//! - `Ord`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::MeshEdge;

impl Eq for MeshEdge {}

impl PartialOrd for MeshEdge {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for MeshEdge {
    /// Orders edges so that a `BinaryHeap` pops the cheapest first.  Equal costs
    /// are broken on the vertex pair (smaller `v0`, then smaller `v1` first), so the
    /// pop order is deterministic.
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other
            .cost
            .partial_cmp(&self.cost)
            .unwrap_or_else(|| other.cost.total_cmp(&self.cost))
            .then_with(|| other.v0.cmp(&self.v0))
            .then_with(|| other.v1.cmp(&self.v1))
    }
}
