// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Lateral (side-slip and yaw) dynamics of the single-track motorcycle model.
//!
//! The bike is a rigid body with the body-fixed frame `x` forward and `y` to the
//! right; yaw rate `r`, steer angle `δ` and lateral force all point to the
//! right.  The front axle is `a` ahead of the centre of mass and the rear axle
//! `b` behind it.  With the lateral velocity `v_y` of the centre of mass and the
//! forward speed `v` the two axles see the slip angles
//!
//! * `α_f = δ − (v_y + a·r) / v`
//! * `α_r = (b·r − v_y) / v`
//!
//! and the equations of motion are
//!
//! * `m (v̇_y + v·r) = F_y,f · cos δ + F_y,r`
//! * `I_z ṙ = a · F_y,f · cos δ − b · F_y,r`
//!
//! The tyre forces follow the slip angles through the tyre relaxation filter
//! ([`PacejkaMotoTire::update_relaxation`]), which is integrated exactly for a
//! slip angle that is constant over a step.
//!
//! # Stiffness
//!
//! The relaxation filter turns the slip-angle dynamics into the oscillator
//! `α̇ = −(v_y + a·r)/σ − (v/σ) α`, `v̇_y = (C_α / m) α`, whose natural frequency
//! `√(C_α / (m σ))` does not depend on the speed: the `1/v` of the slip angle
//! cancels against the `v/σ` relaxation rate.  The explicit scheme used here is
//! stable for steps below `2 / ω`; [`advance`] splits a step into as many equal
//! sub-steps as needed to stay well below that limit.  Explicit integration also
//! follows a physically unstable (oversteering) mode correctly, which a linearly
//! implicit step would not.
//!
//! # Standstill
//!
//! Below [`GRIP_SPEED`] the contact patches do not slide: the slip angles lose
//! their meaning and the bike follows the kinematic no-side-slip motion
//! `r = v δ / L`, `v_y = b r`, which is zero at rest.

use super::PacejkaMotoTire;

/// Forward speed (m/s) up to which the tyres are taken to grip the road.
///
/// At and below it the bike follows the kinematic single-track motion instead
/// of the slip-angle dynamics.
pub(super) const GRIP_SPEED: f64 = 0.1;

/// Safety factor of the sub-step count: a sub-step is at most
/// `1 / (STABILITY_MARGIN · ω)` long, which keeps `h·ω` a factor
/// `2 · STABILITY_MARGIN` below the explicit stability limit `h·ω < 2`.
const STABILITY_MARGIN: f64 = 2.0;

/// Upper bound on the number of sub-steps of one call to [`advance`].
const MAX_SUBSTEPS: usize = 20_000;

/// Smallest relaxation length (m) taken into account for the sub-step count.
const MIN_RELAXATION_LENGTH: f64 = 1.0e-3;

/// Mass properties and axle positions of the bike.
#[derive(Debug, Clone, Copy)]
pub(super) struct Chassis {
    /// Mass (kg).
    pub mass: f64,
    /// Yaw moment of inertia (kg·m²).
    pub yaw_inertia: f64,
    /// Distance from the centre of mass to the front axle (m).
    pub a: f64,
    /// Distance from the centre of mass to the rear axle (m).
    pub b: f64,
}

/// What acts on the lateral dynamics and is held constant during a step.
#[derive(Debug, Clone, Copy)]
pub(super) struct LateralInput {
    /// Forward speed (m/s), never negative.
    pub speed: f64,
    /// Steering angle (rad) at the end of the step, positive = right.
    pub steer: f64,
    /// Lean angle (rad), positive = right; it produces camber thrust.
    pub lean: f64,
    /// Front tyre normal force (N).
    pub front_normal: f64,
    /// Rear tyre normal force (N).
    pub rear_normal: f64,
    /// Front tyre longitudinal slip ratio, which shares the grip with the lateral force.
    pub front_slip_ratio: f64,
    /// Rear tyre longitudinal slip ratio.
    pub rear_slip_ratio: f64,
}

/// Lateral velocity and yaw rate of the bike.
#[derive(Debug, Clone, Copy)]
pub(super) struct LateralState {
    /// Lateral velocity of the centre of mass (m/s), positive = right.
    pub lateral_velocity: f64,
    /// Yaw rate (rad/s), positive = right.
    pub yaw_rate: f64,
}

/// Initial slope `∂F_y/∂α` (N/rad) of the lateral Magic Formula of `tire` at the
/// given normal force.
///
/// `d/dα [D sin(C atan(φ(Bα)))]` at `α = 0` is `B·C·D`, whatever the curvature
/// factor `E` is; combined slip and the peak only lower the slope.
fn cornering_stiffness(tire: &PacejkaMotoTire, normal_force: f64) -> f64 {
    (tire.lat.b * tire.lat.c * tire.lat.d * normal_force).max(0.0)
}

/// Number of equal sub-steps that keep the explicit scheme well inside its
/// stability limit, see the module documentation.
fn substep_count(
    dt: f64,
    chassis: Chassis,
    input: LateralInput,
    front: &PacejkaMotoTire,
    rear: &PacejkaMotoTire,
) -> usize {
    let k_front = cornering_stiffness(front, input.front_normal);
    let k_rear = cornering_stiffness(rear, input.rear_normal);
    let sigma = front
        .relaxation_length
        .min(rear.relaxation_length)
        .max(MIN_RELAXATION_LENGTH);
    let side_force = (k_front + k_rear) / chassis.mass;
    let yaw_moment =
        (chassis.a * chassis.a * k_front + chassis.b * chassis.b * k_rear) / chassis.yaw_inertia;
    let omega = ((side_force + yaw_moment) / sigma).max(0.0).sqrt();
    let count = (STABILITY_MARGIN * dt * omega).ceil();
    if count.is_finite() {
        (count as usize).clamp(1, MAX_SUBSTEPS)
    } else {
        MAX_SUBSTEPS
    }
}

/// Lateral forces `(F_y,f, F_y,r)` (N) of the two tyres at their current
/// filtered slip angles.
fn lateral_forces(
    input: LateralInput,
    front: &PacejkaMotoTire,
    rear: &PacejkaMotoTire,
) -> (f64, f64) {
    let front_force = front
        .forces(
            front.filtered_slip_angle,
            input.front_slip_ratio,
            input.front_normal,
            input.lean,
        )
        .0;
    let rear_force = rear
        .forces(
            rear.filtered_slip_angle,
            input.rear_slip_ratio,
            input.rear_normal,
            input.lean,
        )
        .0;
    (front_force, rear_force)
}

/// Advance the lateral velocity and the yaw rate by `dt` seconds.
///
/// `front` and `rear` carry the relaxation state of the tyres; their filtered
/// slip angles are advanced as well.  Returns the new state and the mean
/// lateral acceleration of the centre of mass over the step, `(F_y,f cos δ +
/// F_y,r) / m` (m/s²), which the roll dynamics need.
///
/// A step of zero (or negative, or non-finite) length changes nothing.
pub(super) fn advance(
    dt: f64,
    chassis: Chassis,
    input: LateralInput,
    state: LateralState,
    front: &mut PacejkaMotoTire,
    rear: &mut PacejkaMotoTire,
) -> (LateralState, f64) {
    let speed = input.speed;
    if dt <= 0.0 || !dt.is_finite() {
        let (front_force, rear_force) = lateral_forces(input, front, rear);
        let lateral_accel = (front_force * input.steer.cos() + rear_force) / chassis.mass;
        return (state, lateral_accel);
    }

    let wheelbase = chassis.a + chassis.b;
    if speed.abs() <= GRIP_SPEED {
        // The contact patches grip: no side slip, the bike pivots about the
        // rear contact point.  The tyres lose their slip with the speed.
        front.update_relaxation(0.0, speed, dt);
        rear.update_relaxation(0.0, speed, dt);
        let yaw_rate = if wheelbase > 0.0 {
            speed * input.steer / wheelbase
        } else {
            0.0
        };
        let kinematic = LateralState {
            lateral_velocity: chassis.b * yaw_rate,
            yaw_rate,
        };
        return (kinematic, speed * yaw_rate);
    }

    let steps = substep_count(dt, chassis, input, front, rear);
    let h = dt / steps as f64;
    let cos_steer = input.steer.cos();
    let mut v_y = state.lateral_velocity;
    let mut r = state.yaw_rate;
    let mut accel_sum = 0.0;
    for _ in 0..steps {
        let raw_front = input.steer - (v_y + chassis.a * r) / speed;
        let raw_rear = (chassis.b * r - v_y) / speed;
        front.update_relaxation(raw_front, speed, h);
        rear.update_relaxation(raw_rear, speed, h);
        let (front_force, rear_force) = lateral_forces(input, front, rear);
        let side_force = front_force * cos_steer + rear_force;
        let yaw_moment = chassis.a * front_force * cos_steer - chassis.b * rear_force;
        let lateral_accel = side_force / chassis.mass;
        v_y += h * (lateral_accel - speed * r);
        r += h * yaw_moment / chassis.yaw_inertia;
        accel_sum += lateral_accel;
    }
    (
        LateralState {
            lateral_velocity: v_y,
            yaw_rate: r,
        },
        accel_sum / steps as f64,
    )
}

#[cfg(test)]
mod tests;
