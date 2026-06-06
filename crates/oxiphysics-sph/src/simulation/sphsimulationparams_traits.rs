//! # SphSimulationParams - Trait Implementations
//!
//! This module contains trait implementations for `SphSimulationParams`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::math::Vec3;

use super::types::{SolverType, SphSimulationParams};

impl Default for SphSimulationParams {
    fn default() -> Self {
        Self {
            gravity: Vec3::new(0.0, -9.81, 0.0),
            smoothing_length: 0.1,
            sound_speed: 100.0,
            rest_density: 1000.0,
            viscosity: 0.01,
            solver_type: SolverType::Wcsph,
            boundary_stiffness: 50_000.0,
            boundary_damping: 500.0,
            min_dt: 1e-6,
            max_dt: 0.005,
        }
    }
}
