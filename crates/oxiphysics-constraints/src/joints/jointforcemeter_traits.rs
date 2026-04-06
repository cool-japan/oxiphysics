//! # JointForceMeter - Trait Implementations
//!
//! This module contains trait implementations for `JointForceMeter`.
//!
//! ## Implemented Traits
//!
//! - `Constraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use crate::traits::Constraint;
use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::Vec3;
use oxiphysics_rigid::RigidBodySet;

use super::types::JointForceMeter;

impl<C: Constraint + Clone + std::fmt::Debug> Constraint for JointForceMeter<C> {
    fn prepare(&mut self, bodies: &RigidBodySet, dt: f64) {
        self.last_impulse_magnitude = 0.0;
        self.last_torque_magnitude = 0.0;
        self.inner.prepare(bodies, dt);
    }
    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, dt: f64) {
        let vel_before: Vec<Vec3> = self
            .inner
            .body_handles()
            .iter()
            .filter_map(|&h| bodies.get(h).map(|b| b.velocity))
            .collect();
        let ang_before: Vec<Vec3> = self
            .inner
            .body_handles()
            .iter()
            .filter_map(|&h| bodies.get(h).map(|b| b.angular_velocity))
            .collect();
        self.inner.solve_velocity(bodies, dt);
        let handles = self.inner.body_handles();
        let mut total_dv = 0.0;
        let mut total_dw = 0.0;
        for (i, &h) in handles.iter().enumerate() {
            if let Some(body) = bodies.get(h) {
                if i < vel_before.len() {
                    let m = if body.inverse_mass > 1e-12 {
                        1.0 / body.inverse_mass
                    } else {
                        0.0
                    };
                    total_dv += m * (body.velocity - vel_before[i]).norm();
                }
                if i < ang_before.len() {
                    total_dw += (body.angular_velocity - ang_before[i]).norm();
                }
            }
        }
        self.last_impulse_magnitude = total_dv;
        self.last_torque_magnitude = total_dw;
    }
    fn solve_position(&mut self, bodies: &mut RigidBodySet, dt: f64) {
        self.inner.solve_position(bodies, dt);
    }
    fn body_handles(&self) -> Vec<BodyHandle> {
        self.inner.body_handles()
    }
}
