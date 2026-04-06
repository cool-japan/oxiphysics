//! # RigidBody - accessors Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use oxiphysics_core::math::Mat3;

use super::types::BodyType;

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Switch body type at runtime (e.g. dynamic <-> kinematic).
    #[allow(dead_code)]
    pub fn set_body_type(&mut self, new_type: BodyType) {
        self.body_type = new_type;
        match new_type {
            BodyType::Static | BodyType::Kinematic => {
                self.inverse_mass = 0.0;
                self.world_inverse_inertia = Mat3::zeros();
            }
            BodyType::Dynamic => {
                if self.mass > 0.0 {
                    self.inverse_mass = 1.0 / self.mass;
                }
                self.update_world_inertia();
            }
        }
    }
}
