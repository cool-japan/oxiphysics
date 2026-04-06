//! # UrdfInertial - Trait Implementations
//!
//! This module contains trait implementations for `UrdfInertial`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::UrdfInertial;

impl Default for UrdfInertial {
    fn default() -> Self {
        Self {
            mass: 1.0,
            com: [0.0; 3],
            inertia_diag: [0.001; 3],
            inertia_off: [0.0; 3],
        }
    }
}
