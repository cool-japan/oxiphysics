//! # EscState - Trait Implementations
//!
//! This module contains trait implementations for `EscState`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::types::EscState;

impl Default for EscState {
    fn default() -> Self {
        Self {
            is_active: false,
            intervention_type: String::from("none"),
            yaw_rate_error: 0.0,
        }
    }
}
