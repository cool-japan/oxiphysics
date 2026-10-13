//! # RigidBody - angular_momentum_group Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use oxiphysics_core::math::Vec3;

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Angular momentum about the center of mass in the world frame:
    /// `I_w * omega`.
    ///
    /// `omega` is the world-frame angular velocity, so the product has to use
    /// the world-space inertia tensor [`RigidBody::world_inertia`] and not the
    /// local one.
    pub fn angular_momentum(&self) -> Vec3 {
        self.world_inertia() * self.angular_velocity
    }
}
