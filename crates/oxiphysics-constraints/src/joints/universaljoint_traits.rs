//! # UniversalJoint - Trait Implementations
//!
//! This module contains trait implementations for `UniversalJoint`.
//!
//! ## Implemented Traits
//!
//! - `Constraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::traits::Constraint;
use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::Vec3;
use oxiphysics_rigid::RigidBodySet;

use super::constraint_math::{
    apply_orientation_correction, solve_angular_rows, solve_point_velocity,
};
use super::functions::{JOINT_BAUMGARTE, read_body};
use super::types::UniversalJoint;

impl Constraint for UniversalJoint {
    fn prepare(&mut self, bodies: &RigidBodySet, _dt: f64) {
        let a = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return,
        };
        let b = match read_body(bodies, self.body_b) {
            Some(p) => p,
            None => return,
        };
        self.r_a = a.rotation * self.local_anchor_a;
        self.r_b = b.rotation * self.local_anchor_b;
        self.world_axis_a = (a.rotation * self.local_axis_a).normalize();
        self.world_axis_b = (b.rotation * self.local_axis_b).normalize();
        self.constrained_axis = self
            .world_axis_a
            .cross(&self.world_axis_b)
            .try_normalize(1e-12)
            .unwrap_or_else(Vec3::zeros);
        let _ = b;
    }
    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, _dt: f64) {
        solve_point_velocity(bodies, self.body_a, self.body_b, &self.r_a, &self.r_b);
        solve_angular_rows(bodies, self.body_a, self.body_b, &[self.constrained_axis]);
    }
    fn solve_position(&mut self, bodies: &mut RigidBodySet, _dt: f64) {
        let a = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return,
        };
        let b = match read_body(bodies, self.body_b) {
            Some(p) => p,
            None => return,
        };
        let world_a = a.position + a.rotation * self.local_anchor_a;
        let world_b = b.position + b.rotation * self.local_anchor_b;
        let error = world_a - world_b;
        let k = a.inv_mass + b.inv_mass;
        if k > 1e-12 {
            let correction = error * JOINT_BAUMGARTE;
            if let Some(body) = bodies.get_mut(self.body_a) {
                body.transform.position -= correction * (a.inv_mass / k);
            }
            if let Some(body) = bodies.get_mut(self.body_b) {
                body.transform.position += correction * (b.inv_mass / k);
            }
        }
        // The free axes must stay perpendicular.  Rotating B by +ε about
        // a × b opens the angle between the axes by ε, so a deviation
        // δ = asin(a·b) from the perpendicular is removed by rotating B by +δ
        // about that axis, i.e. a world-frame error vector of −δ·(a × b)/|a × b|.
        let axis_a = (a.rotation * self.local_axis_a).normalize();
        let axis_b = (b.rotation * self.local_axis_b).normalize();
        let swing_axis = axis_a.cross(&axis_b);
        let sine = swing_axis.norm();
        if sine > 1e-9 {
            let deviation = axis_a.dot(&axis_b).clamp(-1.0, 1.0).asin();
            if deviation.abs() > 1e-6 {
                apply_orientation_correction(
                    bodies,
                    self.body_a,
                    self.body_b,
                    swing_axis * (-deviation / sine),
                    JOINT_BAUMGARTE,
                );
            }
        }
    }
    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_a, self.body_b]
    }
}
