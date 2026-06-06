//! # PenaltyParams - Trait Implementations
//!
//! This module contains trait implementations for `PenaltyParams`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::PenaltyParams;

impl Default for PenaltyParams {
    fn default() -> Self {
        Self {
            stiffness: 1e4,
            damping: 100.0,
            friction: 0.3,
        }
    }
}
