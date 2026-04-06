#[allow(unused_imports)]
use super::functions::*;
// # WheelState - Trait Implementations
//
// This module contains trait implementations for `WheelState`.
//
// ## Implemented Traits
//
// - `Default`
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::math::Vec3;

use super::types::WheelState;

impl Default for WheelState {
    fn default() -> Self {
        Self {
            suspension_length: 0.3,
            suspension_force: 0.0,
            contact_point: Vec3::zeros(),
            contact_normal: Vec3::new(0.0, 1.0, 0.0),
            is_in_contact: false,
            steering_angle: 0.0,
            rotation_angle: 0.0,
            spin_velocity: 0.0,
            slip_angle: 0.0,
            slip_ratio: 0.0,
            lateral_force: 0.0,
            longitudinal_force: 0.0,
            world_position: Vec3::zeros(),
            forward_direction: Vec3::new(0.0, 0.0, 1.0),
            right_direction: Vec3::new(1.0, 0.0, 0.0),
            angular_velocity: 0.0,
            load: 0.0,
            temperature: 25.0,
        }
    }
}
