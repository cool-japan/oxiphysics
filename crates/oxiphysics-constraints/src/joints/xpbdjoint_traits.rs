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

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{apply_pair_impulse, read_body};
use super::types::XpbdJoint;

impl Constraint for XpbdJoint {
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
        self.lambda = 0.0;
    }
    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, dt: f64) {
        let a = match read_body(bodies, self.body_a) {
            Some(p) => p,
            None => return,
        };
        let b = match read_body(bodies, self.body_b) {
            Some(p) => p,
            None => return,
        };
        let world_a = a.position + self.r_a;
        let world_b = b.position + self.r_b;
        let diff = world_b - world_a;
        let dist = diff.norm();
        let dir = if dist > 1e-10 { diff / dist } else { return };
        let extension = dist - self.rest_length;
        let rel_v = (b.vel - a.vel).dot(&dir);
        let k = a.inv_mass + b.inv_mass;
        if k < 1e-12 {
            return;
        }
        let k_spring = if self.compliance > 1e-20 {
            1.0 / (self.compliance * dt * dt)
        } else {
            1e9
        };
        let force = k_spring * extension + self.damping * rel_v;
        let lambda = force * dt;
        let impulse = dir * lambda;
        apply_pair_impulse(
            bodies,
            self.body_a,
            self.body_b,
            impulse,
            &self.r_a,
            &self.r_b,
        );
        self.lambda += lambda;
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
