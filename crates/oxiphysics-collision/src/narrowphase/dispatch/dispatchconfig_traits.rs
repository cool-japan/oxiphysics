//! # DispatchConfig - Trait Implementations
//!
//! This module contains trait implementations for `DispatchConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::DispatchConfig;

impl Default for DispatchConfig {
    fn default() -> Self {
        DispatchConfig {
            use_gjk_fallback: true,
            max_gjk_iterations: 64,
            max_epa_iterations: 64,
            contact_tolerance: 1e-6,
        }
    }
}
