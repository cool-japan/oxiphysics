//! # PiecewiseLinearOpacityMap - Trait Implementations
//!
//! This module contains trait implementations for `PiecewiseLinearOpacityMap`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::PiecewiseLinearOpacityMap;

impl Default for PiecewiseLinearOpacityMap {
    fn default() -> Self {
        let mut m = Self::new();
        m.add_knot(0.0, 0.0);
        m.add_knot(1.0, 1.0);
        m
    }
}
