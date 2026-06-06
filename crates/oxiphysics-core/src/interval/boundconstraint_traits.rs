//! # BoundConstraint - Trait Implementations
//!
//! This module contains trait implementations for `BoundConstraint`.
//!
//! ## Implemented Traits
//!
//! - `IntervalConstraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::IntervalConstraint;
use super::types::{BoundConstraint, Interval};

impl IntervalConstraint for BoundConstraint {
    fn propagate(&self, box_: &mut [Interval]) -> bool {
        if let Some(narrowed) = Interval::intersection(box_[self.idx], self.bounds)
            && (narrowed.lo > box_[self.idx].lo + 1e-15 || narrowed.hi < box_[self.idx].hi - 1e-15)
        {
            box_[self.idx] = narrowed;
            return true;
        }
        false
    }
}
