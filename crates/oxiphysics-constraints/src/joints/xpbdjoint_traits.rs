//! # XpbdJoint - Trait Implementations
//!
//! This module contains trait implementations for `XpbdJoint`.
//!
//! ## Implemented Traits
//!
//! - `Constraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::traits::Constraint;
use oxiphysics_core::BodyHandle;
use oxiphysics_rigid::RigidBodySet;

use super::constraint_math::{SoftRow, spring_softness};
use super::functions::read_body;
use super::types::XpbdJoint;

impl Constraint for XpbdJoint {
    fn prepare(&mut self, bodies: &RigidBodySet, dt: f64) {
        self.lambda = 0.0;
        self.velocity_row = SoftRow::inactive();
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
        let diff = (b.position + self.r_b) - (a.position + self.r_a);
        let dist = diff.norm();
        if dist <= 1e-10 {
            return;
        }
        // Velocity pass: a compliant row regularised by γ = α̃ = α/dt² (with
        // the damping folded in), without a position bias because
        // `solve_position` already removes the positional error.  A rigid joint
        // (α = 0) has γ = 0 and drives the relative speed along the joint to 0.
        let softness = if self.compliance > 1e-20 {
            spring_softness(1.0 / self.compliance, self.damping, dt).map(|(gamma, _)| gamma)
        } else {
            Some(0.0)
        };
        if let Some(gamma) = softness {
            self.velocity_row = SoftRow::new(&a, &b, &self.r_a, &self.r_b, diff / dist, gamma, 0.0);
        }
    }
    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, _dt: f64) {
        self.velocity_row
            .iterate(bodies, self.body_a, self.body_b, &self.r_a, &self.r_b, 1.0);
    }
    fn solve_velocity_relaxed(&mut self, bodies: &mut RigidBodySet, _dt: f64, omega: f64) -> f64 {
        self.velocity_row.iterate(
            bodies,
            self.body_a,
            self.body_b,
            &self.r_a,
            &self.r_b,
            omega,
        )
    }
    fn residual(&self) -> f64 {
        self.velocity_row.last_update()
    }
    fn solve_position(&mut self, bodies: &mut RigidBodySet, dt: f64) {
        if dt < 1e-12 {
            return;
        }
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
        let diff = world_b - world_a;
        let dist = diff.norm();
        let dir = if dist > 1e-10 { diff / dist } else { return };
        let constraint_val = dist - self.rest_length;
        let w = a.inv_mass + b.inv_mass;
        let alpha_tilde = self.compliance / (dt * dt);
        let denom = w + alpha_tilde;
        if denom < 1e-20 {
            return;
        }
        let delta_lambda = (-constraint_val - alpha_tilde * self.lambda) / denom;
        self.lambda += delta_lambda;
        let correction = dir * delta_lambda;
        if let Some(body) = bodies.get_mut(self.body_a) {
            body.transform.position -= correction * a.inv_mass;
        }
        if let Some(body) = bodies.get_mut(self.body_b) {
            body.transform.position += correction * b.inv_mass;
        }
    }
    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_a, self.body_b]
    }
}
