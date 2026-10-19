//! # ConeJoint - Trait Implementations
//!
//! This module contains trait implementations for `ConeJoint`.
//!
//! ## Implemented Traits
//!
//! - `Constraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::traits::Constraint;
use oxiphysics_core::BodyHandle;
use oxiphysics_rigid::RigidBodySet;

use super::constraint_math::{apply_orientation_correction, solve_point_velocity};
use super::functions::{JOINT_BAUMGARTE, read_body};
use super::types::ConeJoint;

impl Constraint for ConeJoint {
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
        self.world_axis = (a.rotation * self.local_axis).normalize();
        let _ = b;
    }
    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, _dt: f64) {
        solve_point_velocity(bodies, self.body_a, self.body_b, &self.r_a, &self.r_b);
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
        let axis_a = (a.rotation * self.local_axis).normalize();
        let axis_b = (b.rotation * self.local_axis).normalize();
        let swing_axis = axis_a.cross(&axis_b);
        let sine = swing_axis.norm();
        let swing_angle = sine.atan2(axis_a.dot(&axis_b));
        let excess = swing_angle - self.half_angle;
        if excess > 1e-6 && sine > 1e-10 {
            // The world rotation that takes axis A onto axis B has the axis
            // a × b and the swing angle; removing the part beyond the cone
            // turns B back (and A forward) by a share set by the inertias.
            apply_orientation_correction(
                bodies,
                self.body_a,
                self.body_b,
                swing_axis * (excess / sine),
                JOINT_BAUMGARTE,
            );
        }
    }
    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_a, self.body_b]
    }
}
