//! # RevoluteJoint - Trait Implementations
//!
//! This module contains trait implementations for `RevoluteJoint`.
//!
//! ## Implemented Traits
//!
//! - `Constraint`
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use crate::traits::Constraint;
use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::{Real, Vec3};
use oxiphysics_rigid::RigidBodySet;

use super::constraint_math::{
    apply_orientation_correction, hinge_alignment_error, hinge_angle, perpendicular_basis,
    solve_angular_rows, solve_point_velocity,
};
use super::functions::{JOINT_BAUMGARTE, apply_angular_impulse, read_body};
use super::types::RevoluteJoint;

/// Relative angular speed `(ω_B − ω_A)·axis` of the hinge and the inverse
/// effective inertia `axisᵀ(I_A⁻¹ + I_B⁻¹)axis` along it.
fn hinge_rate(
    bodies: &RigidBodySet,
    handle_a: BodyHandle,
    handle_b: BodyHandle,
    axis: &Vec3,
) -> Option<(Real, Real)> {
    let a = read_body(bodies, handle_a)?;
    let b = read_body(bodies, handle_b)?;
    let rate = (b.ang_vel - a.ang_vel).dot(axis);
    let inverse_inertia = axis.dot(&((a.inv_inertia + b.inv_inertia) * axis));
    Some((rate, inverse_inertia))
}

impl RevoluteJoint {
    /// Lower and upper angle limits clamped into the representable range
    /// `[-π, π]` of the hinge angle.
    fn clamped_limits(&self) -> (Option<Real>, Option<Real>) {
        (
            self.lower_limit.map(|limit| limit.max(-PI)),
            self.upper_limit.map(|limit| limit.min(PI)),
        )
    }

    /// Unilateral angular rows that keep the hinge inside its limits.
    ///
    /// Each limit is a row `ċ + b ≥ 0` on the rate of approach `ċ` to the limit
    /// with an accumulated impulse clamped to `≥ 0` (a limit can push but never
    /// pull).  While the hinge is still `gap = θ − lower` (or `upper − θ`) away
    /// from a limit the row is speculative, `b = gap/dt`, so the hinge may
    /// approach at full speed but cannot cross the limit during the step; once
    /// the limit is reached the row stops any further approach.  Correcting a
    /// violation that is already there is the job of `solve_position`.
    fn solve_limit_rows(&mut self, bodies: &mut RigidBodySet, dt: Real) {
        if !dt.is_finite() || dt <= 0.0 {
            return;
        }
        let (lower, upper) = self.clamped_limits();
        let axis = self.world_axis;
        if let Some(lower) = lower
            && let Some((rate, inverse_inertia)) =
                hinge_rate(bodies, self.body_a, self.body_b, &axis)
            && inverse_inertia > 1e-12
        {
            let gap = self.cached_angle - lower;
            let speculative_speed = if gap > 0.0 { gap / dt } else { 0.0 };
            let wanted = -(rate + speculative_speed) / inverse_inertia;
            let accumulated = (self.lower_impulse + wanted).max(0.0);
            let applied = accumulated - self.lower_impulse;
            self.lower_impulse = accumulated;
            // B spins up about +axis, A about −axis: θ̇ grows by K·applied.
            apply_angular_impulse(bodies, self.body_a, self.body_b, axis * -applied);
        }
        if let Some(upper) = upper
            && let Some((rate, inverse_inertia)) =
                hinge_rate(bodies, self.body_a, self.body_b, &axis)
            && inverse_inertia > 1e-12
        {
            let gap = upper - self.cached_angle;
            let speculative_speed = if gap > 0.0 { gap / dt } else { 0.0 };
            let wanted = -(-rate + speculative_speed) / inverse_inertia;
            let accumulated = (self.upper_impulse + wanted).max(0.0);
            let applied = accumulated - self.upper_impulse;
            self.upper_impulse = accumulated;
            // B spins down about +axis, A up: θ̇ shrinks by K·applied.
            apply_angular_impulse(bodies, self.body_a, self.body_b, axis * applied);
        }
    }
}

impl Constraint for RevoluteJoint {
    fn prepare(&mut self, bodies: &RigidBodySet, _dt: f64) {
        self.lower_impulse = 0.0;
        self.upper_impulse = 0.0;
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
        self.cached_angle = hinge_angle(
            &a.rotation,
            &b.rotation,
            &self.reference_rotation,
            &self.local_axis,
        );
    }
    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, dt: f64) {
        // The anchors coincide: three linear rows with the exact 3x3 effective
        // mass (lever arms and inertia included).
        solve_point_velocity(bodies, self.body_a, self.body_b, &self.r_a, &self.r_b);
        // Only rotation about the hinge axis is free: lock the two rows
        // perpendicular to it.
        let (first, second) = perpendicular_basis(&self.world_axis);
        solve_angular_rows(bodies, self.body_a, self.body_b, &[first, second]);
        self.solve_limit_rows(bodies, dt);
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
        // Orientation: bring the two hinge axes back into line and push the
        // hinge angle back inside its limits.
        let mut error_world = hinge_alignment_error(
            &a.rotation,
            &b.rotation,
            &self.reference_rotation,
            &self.local_axis,
        );
        let (lower, upper) = self.clamped_limits();
        let angle = hinge_angle(
            &a.rotation,
            &b.rotation,
            &self.reference_rotation,
            &self.local_axis,
        );
        let violation = match (lower, upper) {
            (Some(lower), _) if angle < lower => angle - lower,
            (_, Some(upper)) if angle > upper => angle - upper,
            _ => 0.0,
        };
        error_world += (a.rotation * self.local_axis).normalize() * violation;
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
