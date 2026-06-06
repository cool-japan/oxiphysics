//! # TrustRegionConfig - Trait Implementations
//!
//! This module contains trait implementations for `TrustRegionConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::TrustRegionConfig;

impl Default for TrustRegionConfig {
    fn default() -> Self {
        Self {
            initial_radius: 1.0,
            max_radius: 10.0,
            min_radius: 1e-8,
            accept_ratio: 0.25,
            expand_factor: 2.0,
            shrink_factor: 0.5,
        }
    }
}
