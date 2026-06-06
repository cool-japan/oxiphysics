//! # TensorTopology - Trait Implementations
//!
//! This module contains trait implementations for `TensorTopology`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::TensorTopology;

impl Default for TensorTopology {
    fn default() -> Self {
        Self {
            degeneracy_tolerance: 0.01,
            analyze_2d_slices: false,
        }
    }
}
