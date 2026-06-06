//! # PhysicsConfig - Trait Implementations
//!
//! This module contains trait implementations for `PhysicsConfig`.
//!
//! ## Implemented Traits
//!
//! - `Default`
//! - `Default`
//! - `Default`
//! - `Default`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

pub use crate::math::{Mat3, Quat, Real, Vec3};

use super::types_impl::{PhysicsConfig, PhysicsConfigBuilder, Transform, Transform3D};

impl Default for PhysicsConfig {
    fn default() -> Self {
        Self {
            gravity: Vec3::new(0.0, -9.81, 0.0),
            solver_iterations: 8,
            linear_sleep_threshold: 0.01,
            angular_sleep_threshold: 0.01,
            time_before_sleep: 0.5,
            ccd_enabled: true,
        }
    }
}

impl Default for PhysicsConfigBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl Default for Transform {
    fn default() -> Self {
        Self {
            position: Vec3::zeros(),
            rotation: Quat::identity(),
        }
    }
}

impl Default for Transform3D {
    fn default() -> Self {
        Self {
            position: [0.0; 3],
            rotation: [0.0, 0.0, 0.0, 1.0],
            scale: 1.0,
        }
    }
}
