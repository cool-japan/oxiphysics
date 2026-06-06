//! # WeldJoint - Trait Implementations
//!
//! This module contains trait implementations for `WeldJoint`.
//!
//! ## Implemented Traits
//!
//! - `Constraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::traits::Constraint;
use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::{Unit, Vec3};
use oxiphysics_rigid::RigidBodySet;

use super::functions::{
    apply_angular_impulse, apply_pair_impulse, linear_effective_mass, read_body,
    small_rotation_quat,
};
use super::types::WeldJoint;

impl Constraint for WeldJoint {
    fn prepare(&mut self, bodies: &RigidBodySet, _dt: f64) {
        let a = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return,
        };
        let b = match read_body(bodies, self.body_b) {
            Some(p) => p,
            None => return,
        };
        self.r_a = Vec3::zeros();
        self.r_b = Vec3::zeros();
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
        let _ = b;
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
        let dv = a.vel - b.vel;
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
        let current_offset = b.position - a.position;
        let pos_error = current_offset - self.desired_offset;
        let k = a.inv_mass + b.inv_mass;
        if k > 1e-12 {
            let correction = pos_error * self.baumgarte;
            if let Some(body) = bodies.get_mut(self.body_a) {
                body.transform.position += correction * (a.inv_mass / k);
            }
            if let Some(body) = bodies.get_mut(self.body_b) {
                body.transform.position -= correction * (b.inv_mass / k);
            }
        }
        let current_rel = a.rotation.inverse() * b.rotation;
        let rot_error = self.desired_relative_rotation.inverse() * current_rel;
        let (axis, angle) = rot_error
            .axis_angle()
            .unwrap_or((Unit::new_normalize(Vec3::new(1.0, 0.0, 0.0)), 0.0));
        if angle.abs() > 1e-6 {
            let ang_correction = axis.into_inner() * angle * self.baumgarte;
            let ang_k = a.inv_inertia + b.inv_inertia;
            if let Some(ang_k_inv) = ang_k.try_inverse() {
                let delta = ang_k_inv * ang_correction;
                if let Some(body) = bodies.get_mut(self.body_a) {
                    let dq = small_rotation_quat(&(-delta * 0.5));
                    body.transform.rotation = dq * body.transform.rotation;
                }
                if let Some(body) = bodies.get_mut(self.body_b) {
                    let dq = small_rotation_quat(&(delta * 0.5));
                    body.transform.rotation = dq * body.transform.rotation;
                }
            }
        }
    }
    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_a, self.body_b]
    }
}
