// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Articulated model: a tree of rigid bodies connected by joints.
//!
//! # Configuration vs. velocity coordinates
//!
//! A model has two coordinate spaces:
//!
//! - **configuration** `q` of length [`ArticulatedModel::total_nq`] (`nq`), the
//!   concatenation of every joint's [`Joint::nq`] coordinates — quaternion
//!   joints (spherical, free-floating) carry 4 rotation coordinates here;
//! - **velocity / tangent** `q̇`, `q̈`, `τ` of length
//!   [`ArticulatedModel::total_dof`] (`nv`), the concatenation of every joint's
//!   [`Joint::dof`] coordinates.
//!
//! Use [`ArticulatedModel::q_slice`] / [`ArticulatedModel::v_slice`] to address
//! a joint's part of each, and [`ArticulatedModel::integrate`] /
//! [`ArticulatedModel::difference`] to move between them.

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
    /// Fictitious base acceleration that models gravity (Featherstone §5.3).
    ///
    /// For a gravity field `g_vec` this stores `a_0 = [0; −g_vec]`: the world
    /// frame is given an upward acceleration `−g_vec`, which is dynamically
    /// equivalent to every body being pulled by `g_vec`. RNEA and ABA use it
    /// directly as the root's parent acceleration.
    pub gravity: SpatialVec,
    /// Cumulative DOF offset for each body: `dof_start[i]` is the index into
    /// `q_dot`/`q_ddot`/`tau` slices where body i's coordinates begin.
    dof_start: Vec<usize>,
    /// Total DOF (velocity dimension `nv`) across all bodies.
    total_dof: usize,
    /// Cumulative configuration offset for each body: `q_start[i]` is the
    /// index into `q` where body i's configuration coordinates begin.
    q_start: Vec<usize>,
    /// Total configuration dimension `nq` across all bodies.
    total_nq: usize,
}

impl ArticulatedModel {
    /// Create an empty model with the given gravity vector `[gx, gy, gz]` m/s²
    /// (e.g. `[0.0, 0.0, -9.81]` for gravity pointing along −z).
    ///
    /// Gravity is stored as the fictitious base acceleration
    /// `a_0 = [0; −g_vec]` (see [`ArticulatedModel::gravity`]).
    pub fn new(gravity: [f64; 3]) -> Self {
        let g_spatial = SpatialVec::new([0.0; 3], [-gravity[0], -gravity[1], -gravity[2]]);
        Self {
            bodies: Vec::new(),
            joints: Vec::new(),
            gravity: g_spatial,
            dof_start: Vec::new(),
            total_dof: 0,
            q_start: Vec::new(),
            total_nq: 0,
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
        self.q_start.push(self.total_nq);
        self.total_nq += joint.nq();
        self.bodies.push(body);
        self.joints.push(joint);
        idx
    }

    /// Number of bodies in the model.
    pub fn num_bodies(&self) -> usize {
        self.bodies.len()
    }

    /// Total number of degrees of freedom across all joints (velocity
    /// dimension `nv`: the length of `q_dot`, `q_ddot` and `tau`).
    pub fn total_dof(&self) -> usize {
        self.total_dof
    }

    /// Total number of configuration coordinates (`nq`: the length of `q`).
    ///
    /// Equals [`ArticulatedModel::total_dof`] unless the model contains
    /// quaternion-parameterised joints (spherical: +1, free-floating: +1).
    pub fn total_nq(&self) -> usize {
        self.total_nq
    }

    /// Starting index of body `i`'s joint coordinates in the `q_dot` /
    /// `q_ddot` / `tau` slices.
    pub fn dof_start(&self, i: usize) -> usize {
        self.dof_start[i]
    }

    /// DOF count (velocity dimension) for body `i`'s joint.
    pub fn dof_count(&self, i: usize) -> usize {
        self.joints[i].dof()
    }

    /// Starting index of body `i`'s joint coordinates in `q`.
    pub fn q_start(&self, i: usize) -> usize {
        self.q_start[i]
    }

    /// Configuration-coordinate count for body `i`'s joint.
    pub fn nq_count(&self, i: usize) -> usize {
        self.joints[i].nq()
    }

    /// Slice of the configuration `q` (length [`ArticulatedModel::total_nq`])
    /// belonging to body `i`'s joint (length [`Joint::nq`]).
    pub fn q_slice<'a>(&self, q: &'a [f64], i: usize) -> &'a [f64] {
        let start = self.q_start[i];
        &q[start..start + self.joints[i].nq()]
    }

    /// Slice of a velocity-space vector (`q_dot`, `q_ddot`, `tau`; length
    /// [`ArticulatedModel::total_dof`]) belonging to body `i`'s joint
    /// (length [`Joint::dof`]).
    pub fn v_slice<'a>(&self, v: &'a [f64], i: usize) -> &'a [f64] {
        let start = self.dof_start[i];
        &v[start..start + self.joints[i].dof()]
    }

    /// Neutral configuration: every joint at its [`Joint::neutral`] value
    /// (zero angles/displacements, identity quaternions).
    pub fn neutral_configuration(&self) -> Vec<f64> {
        let mut q = Vec::with_capacity(self.total_nq);
        for joint in &self.joints {
            q.extend(joint.neutral());
        }
        q
    }

    /// Integrate a constant velocity `v` over `dt` on the configuration
    /// manifold: returns `q ⊕ v·dt` (joint-wise [`Joint::integrate`]).
    ///
    /// # Panics
    /// Panics if `q.len() != total_nq()` or `v.len() != total_dof()`.
    pub fn integrate(&self, q: &[f64], v: &[f64], dt: f64) -> Vec<f64> {
        assert_eq!(q.len(), self.total_nq, "q length mismatch");
        assert_eq!(v.len(), self.total_dof, "v length mismatch");
        let mut q_out = q.to_vec();
        for (i, joint) in self.joints.iter().enumerate() {
            let qs = self.q_start[i];
            let nq = joint.nq();
            joint.integrate(
                self.q_slice(q, i),
                self.v_slice(v, i),
                dt,
                &mut q_out[qs..qs + nq],
            );
        }
        q_out
    }

    /// Tangent-space difference `q1 ⊖ q0` (length `total_dof()`), the exact
    /// inverse of [`ArticulatedModel::integrate`]:
    /// `integrate(q0, difference(q0, q1), 1.0) == q1`.
    ///
    /// # Panics
    /// Panics if `q0.len()` or `q1.len()` differs from `total_nq()`.
    pub fn difference(&self, q0: &[f64], q1: &[f64]) -> Vec<f64> {
        assert_eq!(q0.len(), self.total_nq, "q0 length mismatch");
        assert_eq!(q1.len(), self.total_nq, "q1 length mismatch");
        let mut v = vec![0.0; self.total_dof];
        for (i, joint) in self.joints.iter().enumerate() {
            let vs = self.dof_start[i];
            let nv = joint.dof();
            joint.difference(
                self.q_slice(q0, i),
                self.q_slice(q1, i),
                &mut v[vs..vs + nv],
            );
        }
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        joint::{FreeFloatingJoint, RevoluteJoint, SphericalJoint},
        spatial::{SpatialInertia, SpatialTransform},
    };

    fn mixed_model() -> ArticulatedModel {
        let mut model = ArticulatedModel::new([0.0, 0.0, -9.81]);
        let inertia = SpatialInertia::from_com(
            1.0,
            [0.0; 3],
            [[0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]],
        );
        model.add_body(
            RigidBody::new("base", inertia, None, SpatialTransform::IDENTITY),
            Box::new(FreeFloatingJoint),
        );
        model.add_body(
            RigidBody::new("ball", inertia, Some(0), SpatialTransform::IDENTITY),
            Box::new(SphericalJoint),
        );
        model.add_body(
            RigidBody::new("hinge", inertia, Some(1), SpatialTransform::IDENTITY),
            Box::new(RevoluteJoint::new([0.0, 1.0, 0.0])),
        );
        model
    }

    #[test]
    fn test_gravity_stored_as_negated_base_acceleration() {
        let model = ArticulatedModel::new([0.5, -1.0, -9.81]);
        assert_eq!(model.gravity.angular, [0.0; 3]);
        assert_eq!(model.gravity.linear, [-0.5, 1.0, 9.81]);
    }

    #[test]
    fn test_q_and_v_offsets_are_tracked_separately() {
        let model = mixed_model();
        assert_eq!(model.total_nq(), 7 + 4 + 1);
        assert_eq!(model.total_dof(), 6 + 3 + 1);
        assert_eq!(
            (model.q_start(0), model.q_start(1), model.q_start(2)),
            (0, 7, 11)
        );
        assert_eq!(
            (model.dof_start(0), model.dof_start(1), model.dof_start(2)),
            (0, 6, 9)
        );
        let q: Vec<f64> = (0..12).map(f64::from).collect();
        let v: Vec<f64> = (0..10).map(|k| f64::from(k) * 10.0).collect();
        assert_eq!(model.q_slice(&q, 1), &[7.0, 8.0, 9.0, 10.0]);
        assert_eq!(model.v_slice(&v, 1), &[60.0, 70.0, 80.0]);
        assert_eq!(model.q_slice(&q, 2), &[11.0]);
        assert_eq!(model.v_slice(&v, 2), &[90.0]);
        assert_eq!(
            model.neutral_configuration(),
            vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
        );
    }

    #[test]
    fn test_integrate_difference_are_inverse() {
        let model = mixed_model();
        let q0 = model.integrate(
            &model.neutral_configuration(),
            &[0.3, -0.2, 0.5, 0.4, 0.1, -0.7, 0.2, 0.9, -0.3, 1.1],
            1.0,
        );
        let v = [1.0, -0.5, 0.25, -0.3, 0.8, 0.6, -1.2, 0.4, 0.7, -0.9];
        let q1 = model.integrate(&q0, &v, 0.3);
        let back = model.difference(&q0, &q1);
        for (k, (b, vk)) in back.iter().zip(v.iter()).enumerate() {
            assert!(
                (b - 0.3 * vk).abs() < 1e-12,
                "component {k}: {b} vs {}",
                0.3 * vk
            );
        }
        // Quaternions stay unit length.
        for range in [3..7, 7..11] {
            let n: f64 = q1[range].iter().map(|x| x * x).sum::<f64>().sqrt();
            assert!((n - 1.0).abs() < 1e-14);
        }
        // Euclidean coordinates integrate linearly.
        assert!((q1[11] - (q0[11] + 0.3 * v[9])).abs() < 1e-15);
    }
}
