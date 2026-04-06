//! # FluidSimConfig - Trait Implementations
//!
//! This module contains trait implementations for `FluidSimConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::types::FluidSimConfig;

impl Default for FluidSimConfig {
    fn default() -> Self {
        Self {
            grid_size: [16, 16, 16],
            dx: 0.1,
            dt: 0.01,
            density: 1000.0,
            gravity: [0.0, -9.81, 0.0],
            vorticity_eps: 0.5,
            surface_tension: 0.0728,
            pressure_iters: 20,
            flip_ratio: 0.95,
        }
    }
}
