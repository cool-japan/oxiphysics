//! # RigidBody - accessors Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use oxiphysics_core::math::{Quat, Vec3};

use super::types::KinematicTarget;

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Set a kinematic motion target (position + rotation).
    #[allow(dead_code)]
    pub fn set_kinematic_target(&mut self, position: Vec3, rotation: Quat) {
        self.kinematic_target = Some(KinematicTarget {
            target_position: position,
            target_rotation: rotation,
        });
    }
}
