//! # AngularMotorData - Trait Implementations
//!
//! This module contains trait implementations for `AngularMotorData`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::AngularMotorData;

impl Default for AngularMotorData {
    fn default() -> Self {
        Self::new(f64::INFINITY)
    }
}
