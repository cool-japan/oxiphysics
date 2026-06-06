//! # GradientColormap - Trait Implementations
//!
//! This module contains trait implementations for `GradientColormap`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{ColorStop, GradientColormap};

impl Default for GradientColormap {
    fn default() -> Self {
        let mut g = Self::new();
        g.add_stop(ColorStop::new(0.0, 0.0, 0.0, 1.0));
        g.add_stop(ColorStop::new(1.0, 1.0, 0.0, 0.0));
        g
    }
}
