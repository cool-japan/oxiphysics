//! # RigidBody - kinetic_energy_group Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use oxiphysics_core::math::Real;

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Kinetic energy: 0.5 * m * v^2 + 0.5 * omega^T * I * omega.
    #[allow(dead_code)]
    pub fn kinetic_energy(&self) -> Real {
        let lin = 0.5 * self.mass * self.velocity.norm_squared();
        let ang = 0.5
            * self
                .angular_velocity
                .dot(&(self.local_inertia * self.angular_velocity));
        lin + ang
    }
}
