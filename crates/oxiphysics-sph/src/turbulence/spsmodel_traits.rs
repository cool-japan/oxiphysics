//! # SpsModel - Trait Implementations
//!
//! This module contains trait implementations for `SpsModel`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::SpsModel;

impl Default for SpsModel {
    fn default() -> Self {
        Self {
            cs: 0.12,
            ci: 0.0066,
            particle_spacing: 0.1,
        }
    }
}
