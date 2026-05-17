// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Articulated model: a tree of rigid bodies connected by joints.

use crate::{body::RigidBody, joint::Joint, spatial::SpatialVec};

/// An articulated rigid-body model: a kinematic tree.
///
/// Bodies are stored in topological order (parent index < child index), which
/// guarantees that a forward pass (outward from root) can be done by iterating
/// indices 0..n in order, and a backward pass by iterating in reverse.
pub struct ArticulatedModel {
    /// Rigid bodies in topological order.
    pub bodies: Vec<RigidBody>,
    /// Joints connecting body\[i\] to its parent; joints\[i\] belongs to body\[i\].
    pub joints: Vec<Box<dyn Joint>>,
    /// Gravity expressed as a spatial acceleration (linear part only, angular=0).
    ///
    /// Convention: the gravity acceleration is the spatial acceleration of the
    /// world (fixed) frame, i.e. `a_world = -g_vec` (upward for bodies above).
    /// In RNEA we add this as a "fictitious" acceleration to the root's parent.
    pub gravity: SpatialVec,
    /// Cumulative DOF offset for each body: `dof_start[i]` is the index into
    /// `q`/`q_dot`/`tau` slices where body i's coordinates begin.
    dof_start: Vec<usize>,
    /// Total DOF across all bodies.
    total_dof: usize,
}

impl ArticulatedModel {
    /// Create an empty model with the given gravity vector [gx, gy, gz] m/s².
    ///
    /// Gravity is stored as a spatial acceleration applied to the root's
    /// fictitious parent: `a_base = [0; −g_vec]` in the world frame.
    pub fn new(gravity: [f64; 3]) -> Self {
        // Gravity as a spatial vector (no angular component)
        let g_spatial = SpatialVec::new([0.0; 3], gravity);
        Self {
            bodies: Vec::new(),
            joints: Vec::new(),
            gravity: g_spatial,
            dof_start: Vec::new(),
            total_dof: 0,
        }
    }

    /// Add a body with its joint to the model.
    ///
    /// The body must already have its `parent_id` set correctly.
    /// The joint must match the body's DOF intent.
    ///
    /// Returns the index of the newly added body.
    pub fn add_body(&mut self, body: RigidBody, joint: Box<dyn Joint>) -> usize {
        let idx = self.bodies.len();
        self.dof_start.push(self.total_dof);
        self.total_dof += joint.dof();
        self.bodies.push(body);
        self.joints.push(joint);
        idx
    }

    /// Number of bodies in the model.
    pub fn num_bodies(&self) -> usize {
        self.bodies.len()
    }

    /// Total number of degrees of freedom across all joints.
    pub fn total_dof(&self) -> usize {
        self.total_dof
    }

    /// Starting index of body `i`'s joint coordinates in the `q` / `q_dot` / `tau` slices.
    pub fn dof_start(&self, i: usize) -> usize {
        self.dof_start[i]
    }

    /// DOF count for body `i`'s joint.
    pub fn dof_count(&self, i: usize) -> usize {
        self.joints[i].dof()
    }

    /// Slice of `q` belonging to body `i`'s joint.
    pub fn q_slice<'a>(&self, q: &'a [f64], i: usize) -> &'a [f64] {
        let start = self.dof_start[i];
        &q[start..start + self.joints[i].dof()]
    }
}
