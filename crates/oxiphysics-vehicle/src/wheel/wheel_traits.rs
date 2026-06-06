// # Wheel - Trait Implementations
//
// This module contains trait implementations for `Wheel`.
//
// ## Implemented Traits
//
// - `Default`
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::math::Vec3;

use super::types::Wheel;

impl Default for Wheel {
    fn default() -> Self {
        Self {
            radius: 0.3,
            width: 0.2,
            mass: 15.0,
            suspension_rest_length: 0.3,
            suspension_stiffness: 20000.0,
            damping_compression: 4000.0,
            damping_relaxation: 2500.0,
            max_suspension_travel: 0.2,
            friction_slip: 1.0,
            connection_point: Vec3::zeros(),
            suspension_direction: Vec3::new(0.0, -1.0, 0.0),
            axle_direction: Vec3::new(1.0, 0.0, 0.0),
        }
    }
}
