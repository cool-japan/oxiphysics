// # Body - Trait Implementations
//
// This module contains trait implementations for `Body`.
//
// ## Implemented Traits
//
// - `Default`
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::types::Transform;

use super::types::{Body, BodyMode, BodyVelocity};

impl Default for Body {
    fn default() -> Self {
        Self {
            transform: Transform::default(),
            prev_transform: Transform::default(),
            velocity: BodyVelocity::default(),
            force: [0.0; 3],
            torque: [0.0; 3],
            mass: 1.0,
            inv_mass: 1.0,
            inertia: 1.0,
            linear_damping: 0.01,
            angular_damping: 0.01,
            restitution: 0.3,
            mode: BodyMode::Dynamic,
            active: true,
            sleep_timer: 0.0,
            generation: 0,
        }
    }
}
