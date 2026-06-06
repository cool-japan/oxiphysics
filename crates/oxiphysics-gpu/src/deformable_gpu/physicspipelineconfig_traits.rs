//! # PhysicsPipelineConfig - Trait Implementations
//!
//! This module contains trait implementations for `PhysicsPipelineConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::types::PhysicsPipelineConfig;

impl Default for PhysicsPipelineConfig {
    fn default() -> Self {
        Self {
            gravity: [0.0, -9.81, 0.0],
            dt: 1.0 / 60.0,
            substeps: 4,
            collision_enabled: true,
            fem_enabled: true,
            xpbd_enabled: true,
        }
    }
}
