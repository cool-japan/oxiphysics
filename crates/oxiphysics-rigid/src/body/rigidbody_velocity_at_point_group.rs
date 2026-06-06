//! # RigidBody - velocity_at_point_group Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use oxiphysics_core::math::Vec3;

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// World-space velocity at a world-space point on the body.
    pub fn velocity_at_point(&self, point: Vec3) -> Vec3 {
        let r = point - self.transform.position;
        self.velocity + self.angular_velocity.cross(&r)
    }
}
