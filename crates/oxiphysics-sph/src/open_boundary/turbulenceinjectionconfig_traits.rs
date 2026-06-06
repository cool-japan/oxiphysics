//! # TurbulenceInjectionConfig - Trait Implementations
//!
//! This module contains trait implementations for `TurbulenceInjectionConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::TurbulenceInjectionConfig;

impl Default for TurbulenceInjectionConfig {
    fn default() -> Self {
        Self {
            intensity: 0.05,
            length_scale: 0.01,
            seed: 42,
        }
    }
}
