//! # State - Trait Implementations
//!
//! This module contains trait implementations for `State`.
//!
//! ## Implemented Traits
//!
//! - `Eq`
//! - `PartialOrd`
//! - `Ord`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::cmp::Ordering;

use super::types::State;

impl Eq for State {}

impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for State {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .cost
            .partial_cmp(&self.cost)
            .unwrap_or(Ordering::Equal)
    }
}
