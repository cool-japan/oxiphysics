//! # LinearConstraint - Trait Implementations
//!
//! This module contains trait implementations for `LinearConstraint`.
//!
//! ## Implemented Traits
//!
//! - `IntervalConstraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::IntervalConstraint;
use super::types::{Interval, LinearConstraint};

impl IntervalConstraint for LinearConstraint {
    fn propagate(&self, box_: &mut [Interval]) -> bool {
        let mut changed = false;
        for &(target_idx, target_coeff) in &self.coeffs {
            if target_coeff.abs() < 1e-30 {
                continue;
            }
            let mut other_sum = Interval::point(0.0);
            for &(idx, coeff) in &self.coeffs {
                if idx != target_idx {
                    other_sum = other_sum + box_[idx] * coeff;
                }
            }
            let needed = (self.rhs - other_sum) / Interval::point(target_coeff);
            if let Some(narrowed) = Interval::intersection(box_[target_idx], needed)
                && (narrowed.lo > box_[target_idx].lo + 1e-15
                    || narrowed.hi < box_[target_idx].hi - 1e-15)
            {
                box_[target_idx] = narrowed;
                changed = true;
            }
        }
        changed
    }
}
