//! # JointLimits - Trait Implementations
//!
//! This module contains trait implementations for `JointLimits`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::JointLimits;

impl Default for JointLimits {
    fn default() -> Self {
        Self {
            lower: -std::f64::consts::PI,
            upper: std::f64::consts::PI,
            effort: 100.0,
            velocity: 2.0,
        }
    }
}
