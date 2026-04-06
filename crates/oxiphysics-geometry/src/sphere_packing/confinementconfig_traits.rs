//! # ConfinementConfig - Trait Implementations
//!
//! This module contains trait implementations for `ConfinementConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::{ConfinementConfig, ConfinementShape};

impl Default for ConfinementConfig {
    fn default() -> Self {
        Self {
            shape: ConfinementShape::Box,
            dimensions: [1.0, 1.0, 1.0],
            sphere_radius: 0.05,
            max_attempts: 5000,
        }
    }
}
