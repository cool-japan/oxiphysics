//! # RigidBody - angular_momentum_group Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use oxiphysics_core::math::Vec3;

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Angular momentum about the center of mass: I * omega.
    pub fn angular_momentum(&self) -> Vec3 {
        self.local_inertia * self.angular_velocity
    }
}
