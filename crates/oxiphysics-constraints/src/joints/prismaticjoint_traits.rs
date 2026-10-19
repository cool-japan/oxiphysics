//! # PrismaticJoint - Trait Implementations
//!
//! This module contains trait implementations for `PrismaticJoint`.
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

use super::constraint_math::{solve_angular_rows, solve_point_velocity_in_plane};
use super::functions::{JOINT_BAUMGARTE, read_body};
use super::types::PrismaticJoint;

impl Constraint for PrismaticJoint {
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
        let reference = if self.world_axis.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        self.perp1 = self.world_axis.cross(&reference).normalize();
        self.perp2 = self.world_axis.cross(&self.perp1);
        let _ = b;
    }
    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, _dt: f64) {
        // Two linear rows perpendicular to the sliding axis, solved with the
        // exact 2x2 effective mass at the anchors ...
        solve_point_velocity_in_plane(
            bodies,
            self.body_a,
            self.body_b,
            &self.r_a,
            &self.r_b,
            &self.perp1,
            &self.perp2,
        );
        // ... and all three relative rotations locked.
        solve_angular_rows(
            bodies,
            self.body_a,
            self.body_b,
            &[Vec3::x(), Vec3::y(), Vec3::z()],
        );
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
        let along = error.dot(&self.world_axis) * self.world_axis;
        let perp_error = error - along;
        let k = a.inv_mass + b.inv_mass;
        if k < 1e-12 {
            return;
        }
        let correction = perp_error * JOINT_BAUMGARTE;
        if let Some(body) = bodies.get_mut(self.body_a) {
            body.transform.position -= correction * (a.inv_mass / k);
        }
        if let Some(body) = bodies.get_mut(self.body_b) {
            body.transform.position += correction * (b.inv_mass / k);
        }
        if let (Some(lower), Some(upper)) = (self.lower_limit, self.upper_limit) {
            let dist = along.dot(&self.world_axis);
            if dist < lower {
                let limit_correction = self.world_axis * ((lower - dist) * JOINT_BAUMGARTE);
                if let Some(body) = bodies.get_mut(self.body_a) {
                    body.transform.position += limit_correction * (a.inv_mass / k);
                }
                if let Some(body) = bodies.get_mut(self.body_b) {
                    body.transform.position -= limit_correction * (b.inv_mass / k);
                }
            } else if dist > upper {
                let limit_correction = self.world_axis * ((upper - dist) * JOINT_BAUMGARTE);
                if let Some(body) = bodies.get_mut(self.body_a) {
                    body.transform.position += limit_correction * (a.inv_mass / k);
                }
                if let Some(body) = bodies.get_mut(self.body_b) {
                    body.transform.position -= limit_correction * (b.inv_mass / k);
                }
            }
        }
    }
    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_a, self.body_b]
    }
}
