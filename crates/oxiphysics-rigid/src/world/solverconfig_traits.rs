//! # SolverConfig - Trait Implementations
//!
//! This module contains trait implementations for `SolverConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::types::SolverConfig;

impl Default for SolverConfig {
    fn default() -> Self {
        Self {
            velocity_iterations: 8,
            position_iterations: 3,
            baumgarte_factor: 0.2,
            restitution_threshold: 1.0,
            substeps: 1,
            restitution: 0.3,
            contact_hertz: 0.0,
            contact_damping_ratio: 1.0,
            use_soft_contacts: false,
            use_speculative: false,
            speculative_margin: 0.0,
        }
    }
}
