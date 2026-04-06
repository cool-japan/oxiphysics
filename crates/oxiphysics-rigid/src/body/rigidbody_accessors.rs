//! # RigidBody - accessors Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use oxiphysics_core::MassProperties;
use oxiphysics_core::math::{Mat3, Real};
use oxiphysics_geometry::Shape;

use super::types::BodyType;

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Set mass properties from computed values.
    pub fn set_mass_properties(&mut self, props: &MassProperties) {
        self.mass = props.mass;
        self.inverse_mass = props.inverse_mass();
        self.local_inertia = props.local_inertia;
        self.update_world_inertia();
    }
    /// Compute and set mass and inertia from a shape and density.
    /// Overrides current mass and inertia.
    pub fn set_mass_from_shape(&mut self, shape: &dyn Shape, density: Real) {
        let props = shape.mass_properties(density);
        self.set_mass_properties(&props);
    }
    /// Update world-space inverse inertia from local inertia and current rotation.
    pub fn update_world_inertia(&mut self) {
        if self.body_type == BodyType::Static || self.body_type == BodyType::Kinematic {
            self.world_inverse_inertia = Mat3::zeros();
            return;
        }
        let rot = self.transform.rotation.to_rotation_matrix();
        let inv_local = self.local_inertia.try_inverse().unwrap_or_else(Mat3::zeros);
        self.world_inverse_inertia = rot.matrix() * inv_local * rot.matrix().transpose();
    }
}
