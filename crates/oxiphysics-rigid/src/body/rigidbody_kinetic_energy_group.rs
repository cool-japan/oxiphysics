//! # RigidBody - kinetic_energy_group Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use oxiphysics_core::math::Real;

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Kinetic energy: `0.5 * m * v^2 + 0.5 * omega . (I_w * omega)`.
    ///
    /// `omega` is the world-frame angular velocity, so the rotational term uses
    /// the world-space inertia tensor [`RigidBody::world_inertia`] (the local
    /// tensor rotated into the world frame) and is therefore correct for a
    /// rotated, anisotropic body.
    pub fn kinetic_energy(&self) -> Real {
        let lin = 0.5 * self.mass * self.velocity.norm_squared();
        let ang = 0.5
            * self
                .angular_velocity
                .dot(&(self.world_inertia() * self.angular_velocity));
        lin + ang
    }
}
