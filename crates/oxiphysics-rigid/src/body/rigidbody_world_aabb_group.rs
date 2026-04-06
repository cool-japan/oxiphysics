//! # RigidBody - world_aabb_group Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use oxiphysics_core::math::Vec3;

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Compute the world-space AABB of this body given a half-extent size.
    #[allow(dead_code)]
    pub fn world_aabb(&self, half_extents: &Vec3) -> (Vec3, Vec3) {
        let p = self.transform.position;
        let min = p - half_extents;
        let max = p + half_extents;
        (min, max)
    }
}
