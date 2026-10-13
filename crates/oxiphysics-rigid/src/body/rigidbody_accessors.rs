//! # RigidBody - accessors Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
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
    /// World-space inertia tensor `R * I_local * R^T` for the current orientation.
    ///
    /// The angular velocity stored in [`RigidBody`] is expressed in the world
    /// frame, so every quantity that pairs it with the inertia tensor (angular
    /// momentum `L = I_w * omega`, rotational kinetic energy
    /// `0.5 * omega . (I_w * omega)`) has to use this tensor and not
    /// `local_inertia`. For a body that has never been rotated the two are
    /// identical.
    pub fn world_inertia(&self) -> Mat3 {
        let rot = self.transform.rotation.to_rotation_matrix();
        rot.matrix() * self.local_inertia * rot.matrix().transpose()
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
