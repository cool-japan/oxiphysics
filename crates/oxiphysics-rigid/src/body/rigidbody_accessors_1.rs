//! # RigidBody - accessors Methods
//!
//! This module contains method implementations for `RigidBody`.
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use oxiphysics_geometry::Shape;

use super::rigidbody_type::RigidBody;

impl RigidBody {
    /// Set inertia tensor from a shape, keeping the current mass.
    pub fn set_inertia_from_shape(&mut self, shape: &dyn Shape) {
        self.local_inertia = shape.inertia_tensor(self.mass);
        self.update_world_inertia();
    }
}
