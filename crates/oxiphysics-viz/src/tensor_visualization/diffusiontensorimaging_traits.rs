//! # DiffusionTensorImaging - Trait Implementations
//!
//! This module contains trait implementations for `DiffusionTensorImaging`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{DiffusionTensorImaging, FaColorMode};

impl Default for DiffusionTensorImaging {
    fn default() -> Self {
        Self {
            clamp_fa: true,
            tracking_fa_threshold: 0.15,
            color_mode: FaColorMode::Directional,
        }
    }
}
