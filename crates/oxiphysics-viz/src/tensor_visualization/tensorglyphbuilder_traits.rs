//! # TensorGlyphBuilder - Trait Implementations
//!
//! This module contains trait implementations for `TensorGlyphBuilder`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::TensorGlyphBuilder;

impl Default for TensorGlyphBuilder {
    fn default() -> Self {
        Self {
            scale: 1.0,
            min_axis: 1e-4,
            max_axis: 10.0,
            abs_eigenvalues: true,
            positive_color: [0.8, 0.2, 0.2, 1.0],
            negative_color: [0.2, 0.2, 0.8, 1.0],
        }
    }
}
