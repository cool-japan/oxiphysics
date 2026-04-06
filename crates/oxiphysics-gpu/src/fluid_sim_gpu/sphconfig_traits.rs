//! # SphConfig - Trait Implementations
//!
//! This module contains trait implementations for `SphConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::SphConfig;

impl Default for SphConfig {
    fn default() -> Self {
        Self {
            h: 0.1,
            rest_density: 1000.0,
            pressure_k: 100.0,
            viscosity: 0.1,
            surface_tension: 0.0728,
            gravity: [0.0, -9.81, 0.0],
            dt: 0.001,
        }
    }
}
