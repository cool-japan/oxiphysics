//! # FixedJoint - Trait Implementations
//!
//! This module contains trait implementations for `FixedJoint`.
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

use super::constraint_math::apply_orientation_correction;
use super::functions::{
    JOINT_BAUMGARTE, apply_angular_impulse, apply_pair_impulse, linear_effective_mass, read_body,
};
use super::types::FixedJoint;

impl Constraint for FixedJoint {
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
        let ex = Vec3::new(1.0, 0.0, 0.0);
        let ey = Vec3::new(0.0, 1.0, 0.0);
        let ez = Vec3::new(0.0, 0.0, 1.0);
        self.eff_mass_x = linear_effective_mass(
            a.inv_mass,
            b.inv_mass,
            &a.inv_inertia,
            &b.inv_inertia,
            &self.r_a,
            &self.r_b,
            &ex,
        );
        self.eff_mass_y = linear_effective_mass(
            a.inv_mass,
            b.inv_mass,
            &a.inv_inertia,
            &b.inv_inertia,
            &self.r_a,
            &self.r_b,
            &ey,
        );
        self.eff_mass_z = linear_effective_mass(
            a.inv_mass,
            b.inv_mass,
            &a.inv_inertia,
            &b.inv_inertia,
            &self.r_a,
            &self.r_b,
            &ez,
        );
    }
    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, _dt: f64) {
        let a = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return,
        };
        let b = match read_body(bodies, self.body_b) {
            Some(p) => p,
            None => return,
        };
        let va = a.vel + a.ang_vel.cross(&self.r_a);
        let vb = b.vel + b.ang_vel.cross(&self.r_b);
        let dv = va - vb;
        let lambda_x = -dv.x * self.eff_mass_x;
        let lambda_y = -dv.y * self.eff_mass_y;
        let lambda_z = -dv.z * self.eff_mass_z;
        let impulse = Vec3::new(lambda_x, lambda_y, lambda_z);
        apply_pair_impulse(
            bodies,
            self.body_a,
            self.body_b,
            impulse,
            &self.r_a,
            &self.r_b,
        );
        let a2 = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return,
        };
        let b2 = match read_body(bodies, self.body_b) {
            Some(p) => p,
            None => return,
        };
        let d_omega = a2.ang_vel - b2.ang_vel;
        let ang_k = a2.inv_inertia + b2.inv_inertia;
        if let Some(ang_k_inv) = ang_k.try_inverse() {
            let ang_impulse = -(ang_k_inv * d_omega);
            apply_angular_impulse(bodies, self.body_a, self.body_b, ang_impulse);
        }
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
        if k >= 1e-12 {
            let correction = error * JOINT_BAUMGARTE;
            if let Some(body) = bodies.get_mut(self.body_a) {
                body.transform.position -= correction * (a.inv_mass / k);
            }
            if let Some(body) = bodies.get_mut(self.body_b) {
                body.transform.position += correction * (b.inv_mass / k);
            }
        }
        // With B = A·target·e the orientation error e (axis · angle in B's
        // frame) is the rotation that must be removed; in the world frame it is
        // the same rotation conjugated by A·target.
        let target_frame = a.rotation * self.target_rotation;
        let rot_error = self.target_rotation.inverse() * (a.rotation.inverse() * b.rotation);
        let error_world = target_frame * rot_error.scaled_axis();
        apply_orientation_correction(
            bodies,
            self.body_a,
            self.body_b,
            error_world,
            JOINT_BAUMGARTE,
        );
    }
    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_a, self.body_b]
    }
}
