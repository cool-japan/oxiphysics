// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Planar vehicle dynamics behind [`PyVehicle::step`].
//!
//! # Frames and sign conventions
//!
//! The chassis moves on the horizontal `X`-`Z` plane of a `Y`-up world.
//!
//! * At the heading angle `ψ = orientation[2]` the chassis points along
//!   `(cos ψ, 0, sin ψ)`; `ψ = 0` faces `+X`, which is the direction
//!   [`PyVehicle::forward_speed`] measures. The lateral axis is the heading
//!   turned by `+90°` in the plan view, `(-sin ψ, 0, cos ψ)`.
//! * The yaw rate is `angular_velocity[2]`, the rate of change of `ψ`. A positive
//!   steering angle, yaw rate and slip angle all turn toward the lateral axis.
//! * `wheel_positions[i]` is `[forward, up, lateral]` relative to the centre of
//!   mass (wheels are ordered `[FL, FR, RL, RR]`; the left wheels have a positive
//!   lateral offset).
//!
//! # Tyre kinematics
//!
//! A wheel with steering angle `δ` sees the ground slide under its contact patch
//! with the velocity of the chassis point above it,
//! `(u - r·b, v + r·a)` for the body-frame velocity `(u, v)`, yaw rate `r` and
//! wheel offset `(a, b)`, rotated by `δ` into the wheel frame. With the
//! longitudinal and lateral components `(vx, vy)` of that velocity
//!
//! * the slip ratio is `κ = (ωR - vx) / max(|vx|, |ωR|)`,
//! * the slip angle is `α = atan(-vy / |vx|)`,
//!
//! the definitions of [`compute_slip_ratio`](super::compute_slip_ratio) and
//! [`compute_slip_angle`](super::compute_slip_angle). Both reference speeds are
//! floored at [`SLIP_REFERENCE_SPEED`] so that the tyre forces stay finite,
//! continuous and well conditioned when the vehicle stands still; above that
//! speed the definitions are exactly the classical ones. The tyre models return
//! `(Fx, Fy)` in the wheel frame with `Fy` having the sign of `α`, so the
//! lateral force always opposes the lateral sliding, also when reversing.
//!
//! # Integration
//!
//! Tyre forces are stiff: a wheel's slip changes by `R²/(I·v)` per unit tyre
//! stiffness, which makes an explicit update unstable at any practical step.
//! Each call of [`PyVehicle::step`] is therefore split into sub-steps of at most
//! [`MAX_SUBSTEP`] seconds, and within a sub-step
//!
//! * the spin of every wheel is advanced *implicitly*: the new wheel speed solves
//!   `ω' = ω + h/I · (T_drive - T_brake(ω') - R·Fx(κ(ω')))`, so the force that
//!   brakes the wheel is the force that accelerates the chassis (momentum is
//!   exchanged exactly between the two),
//! * the chassis velocity, yaw rate, heading and position are then advanced with
//!   symplectic Euler using those tyre forces.

use super::{PyTireModel, PyVehicle};
use std::f64::consts::TAU;

/// Standard gravity (m/s²).
const GRAVITY: f64 = 9.81;

/// Aerodynamic drag coefficient `F = c·|v|·v` (N·s²/m²).
const DRAG_COEFFICIENT: f64 = 0.3;

/// Brake torque of one wheel at full pedal (N·m).
pub(super) const MAX_BRAKE_TORQUE: f64 = 2000.0;

/// Longest internal integration step (s).
pub(super) const MAX_SUBSTEP: f64 = 1.0e-3;

/// Largest number of internal steps of one [`PyVehicle::step`] call. Longer
/// calls use proportionally longer sub-steps instead of blocking the
/// interpreter.
const MAX_SUBSTEPS: usize = 20_000;

/// Floor (m/s) of the reference speed of the slip ratio and the slip angle.
pub(super) const SLIP_REFERENCE_SPEED: f64 = 0.5;

/// Wheel speed (rad/s) below which the brake torque fades linearly to zero.
///
/// A brake torque has the sign of the wheel speed; the fade replaces the
/// discontinuity at standstill, so that a braked wheel comes to rest and stays
/// there instead of chattering around zero.
const BRAKE_FADE_SPEED: f64 = 0.5;

/// Number of doublings the root bracket may grow by.
const MAX_BRACKET_DOUBLINGS: usize = 80;

/// Iteration cap of the root finder.
const MAX_ROOT_ITERATIONS: usize = 100;

/// Relative width of the bracket at which the root finder stops.
const ROOT_TOLERANCE: f64 = 1.0e-13;

/// Slip ratio of a wheel with the given peripheral speed over a ground that
/// moves under it at `longitudinal_speed` (see the module documentation).
pub(super) fn regularised_slip_ratio(peripheral_speed: f64, longitudinal_speed: f64) -> f64 {
    let reference = longitudinal_speed
        .abs()
        .max(peripheral_speed.abs())
        .max(SLIP_REFERENCE_SPEED);
    (peripheral_speed - longitudinal_speed) / reference
}

/// Slip angle of a wheel whose contact patch slides sideways with
/// `lateral_speed` while rolling at `longitudinal_speed`.
pub(super) fn regularised_slip_angle(lateral_speed: f64, longitudinal_speed: f64) -> f64 {
    (-lateral_speed / longitudinal_speed.abs().max(SLIP_REFERENCE_SPEED)).atan()
}

/// Control input limited to `[lo, hi]`; NaN counts as zero.
fn limited_input(value: f64, lo: f64, hi: f64) -> f64 {
    if value.is_nan() {
        0.0
    } else {
        value.clamp(lo, hi)
    }
}

/// `value`, or zero if it is NaN or infinite.
fn finite_or_zero(value: f64) -> f64 {
    if value.is_finite() { value } else { 0.0 }
}

/// Angle wrapped into `(-π, π]`.
fn wrap_angle(angle: f64) -> f64 {
    if !angle.is_finite() {
        return angle;
    }
    let wrapped = angle - TAU * (angle / TAU).round();
    if wrapped <= -std::f64::consts::PI {
        wrapped + TAU
    } else {
        wrapped
    }
}

/// Root of `f` near `start`, for a continuous `f` that is negative far below
/// and positive far above `start` (for example `x - c - g(x)` with a bounded
/// `g`).
///
/// A bracket of half-width `half_width` around `start` is doubled until it
/// encloses a sign change; the Illinois variant of regula falsi then closes in
/// on the root with a guaranteed bracket. `start` is returned if no sign change
/// exists (non-finite function values).
fn bracketed_root(f: impl Fn(f64) -> f64, start: f64, half_width: f64) -> f64 {
    let mut width = if half_width.is_finite() && half_width > 0.0 {
        half_width
    } else {
        1.0
    };
    let (mut lo, mut hi) = (start - width, start + width);
    let (mut f_lo, mut f_hi) = (f(lo), f(hi));
    for _ in 0..MAX_BRACKET_DOUBLINGS {
        if f_lo.is_nan() || f_hi.is_nan() {
            return start;
        }
        if f_lo <= 0.0 && f_hi >= 0.0 {
            break;
        }
        width *= 2.0;
        if f_lo > 0.0 {
            lo = start - width;
            f_lo = f(lo);
        }
        if f_hi < 0.0 {
            hi = start + width;
            f_hi = f(hi);
        }
    }
    if !(f_lo <= 0.0 && f_hi >= 0.0) {
        return start;
    }

    // Which end of the bracket was retained by the previous iteration.
    let mut retained = 0i8;
    let mut x = 0.5 * (lo + hi);
    for _ in 0..MAX_ROOT_ITERATIONS {
        let slope = f_hi - f_lo;
        x = if slope > 0.0 {
            (lo * f_hi - hi * f_lo) / slope
        } else {
            0.5 * (lo + hi)
        };
        if !(x > lo && x < hi) {
            x = 0.5 * (lo + hi);
        }
        let f_x = f(x);
        if f_x == 0.0 || (hi - lo) <= ROOT_TOLERANCE * (1.0 + x.abs()) {
            return x;
        }
        if f_x > 0.0 {
            hi = x;
            f_hi = f_x;
            if retained == -1 {
                f_lo *= 0.5;
            }
            retained = -1;
        } else {
            lo = x;
            f_lo = f_x;
            if retained == 1 {
                f_hi *= 0.5;
            }
            retained = 1;
        }
    }
    x
}

/// Velocity of the ground under a wheel contact patch, in the wheel frame.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
struct PatchVelocity {
    /// Component along the wheel heading (m/s).
    longitudinal: f64,
    /// Component along the wheel axle, positive toward the lateral axis (m/s).
    lateral: f64,
}

/// The spin equation of one wheel over one sub-step.
struct WheelSpin<'a> {
    /// Tyre model producing the longitudinal force.
    tire: &'a PyTireModel,
    /// Wheel speed at the start of the sub-step (rad/s).
    omega: f64,
    /// Wheel moment of inertia (kg·m²).
    inertia: f64,
    /// Effective rolling radius (m).
    radius: f64,
    /// Sub-step length (s).
    step: f64,
    /// Drive torque at the wheel (N·m).
    drive_torque: f64,
    /// Brake torque at full lock (N·m).
    brake_torque: f64,
    /// Ground speed under the contact patch along the wheel heading (m/s).
    longitudinal_speed: f64,
    /// Slip angle of the contact patch (rad).
    slip_angle: f64,
    /// Normal load (N).
    normal_load: f64,
}

impl WheelSpin<'_> {
    /// Longitudinal tyre force when the wheel turns at `omega` (N).
    fn longitudinal_force(&self, omega: f64) -> f64 {
        let kappa = regularised_slip_ratio(omega * self.radius, self.longitudinal_speed);
        finite_or_zero(
            self.tire
                .compute_forces(self.slip_angle, kappa, self.normal_load)
                .0,
        )
    }

    /// Residual of the implicit spin equation at the candidate speed `omega`.
    fn residual(&self, omega: f64) -> f64 {
        let brake = self.brake_torque * (omega / BRAKE_FADE_SPEED).clamp(-1.0, 1.0);
        let net_torque = self.drive_torque - brake - self.radius * self.longitudinal_force(omega);
        omega - self.omega - self.step / self.inertia * net_torque
    }

    /// Wheel speed at the end of the sub-step (rad/s).
    fn solve(&self) -> f64 {
        let reach = self.step / self.inertia * (self.drive_torque.abs() + self.brake_torque.abs());
        bracketed_root(|omega| self.residual(omega), self.omega, reach.max(1.0e-3))
    }
}

impl PyVehicle {
    /// Body-frame planar velocity `(u, v, r)`: forward speed, lateral speed and
    /// yaw rate.
    pub(super) fn planar_velocity(&self) -> (f64, f64, f64) {
        let (sin_yaw, cos_yaw) = self.orientation[2].sin_cos();
        let forward = self.velocity[0] * cos_yaw + self.velocity[2] * sin_yaw;
        let lateral = -self.velocity[0] * sin_yaw + self.velocity[2] * cos_yaw;
        (forward, lateral, self.angular_velocity[2])
    }

    /// Normal load of every wheel at rest (N), `[FL, FR, RL, RR]`.
    ///
    /// The weight is shared between the axles by the longitudinal position of
    /// the axles relative to the centre of mass, and between the two wheels of
    /// an axle by their lateral positions, so that the loads balance the weight
    /// and its moments. A wheel is never loaded negatively.
    pub(super) fn static_load_array(&self) -> [f64; 4] {
        let weight = self.mass * GRAVITY;
        let front_axle = 0.5 * (self.wheel_positions[0][0] + self.wheel_positions[1][0]);
        let rear_axle = 0.5 * (self.wheel_positions[2][0] + self.wheel_positions[3][0]);
        let wheelbase = front_axle - rear_axle;
        let (front_weight, rear_weight) = if wheelbase.abs() > 1.0e-9 {
            (
                weight * (-rear_axle) / wheelbase,
                weight * front_axle / wheelbase,
            )
        } else {
            (0.5 * weight, 0.5 * weight)
        };
        let split = |axle_weight: f64, left: usize, right: usize| -> [f64; 2] {
            let left_offset = self.wheel_positions[left][2];
            let right_offset = self.wheel_positions[right][2];
            let track = left_offset - right_offset;
            if track.abs() > 1.0e-9 {
                [
                    (axle_weight * (-right_offset) / track).max(0.0),
                    (axle_weight * left_offset / track).max(0.0),
                ]
            } else {
                [0.5 * axle_weight, 0.5 * axle_weight]
            }
        };
        let front = split(front_weight, 0, 1);
        let rear = split(rear_weight, 2, 3);
        [front[0], front[1], rear[0], rear[1]]
    }

    /// Heading of every wheel relative to the chassis (rad): the Ackermann
    /// angles of the front axle, the rear steering angle if four-wheel steering
    /// is active, plus each wheel's toe.
    fn wheel_steer_angles(&self) -> [f64; 4] {
        let (left, right) = self.steering.ackermann_front_angles();
        let rear = if self.steering.four_ws {
            self.steering.rear_angle()
        } else {
            0.0
        };
        [
            left + self.wheels[0].toe,
            right + self.wheels[1].toe,
            rear + self.wheels[2].toe,
            rear + self.wheels[3].toe,
        ]
    }

    /// Ground velocity under every contact patch, in the wheel frames.
    fn patch_velocities(&self, steer: &[f64; 4]) -> [PatchVelocity; 4] {
        let (forward, lateral, yaw_rate) = self.planar_velocity();
        let mut patches = [PatchVelocity::default(); 4];
        for (patch, (position, &angle)) in patches
            .iter_mut()
            .zip(self.wheel_positions.iter().zip(steer))
        {
            let ahead = forward - yaw_rate * position[2];
            let sideways = lateral + yaw_rate * position[0];
            let (sin_steer, cos_steer) = angle.sin_cos();
            *patch = PatchVelocity {
                longitudinal: ahead * cos_steer + sideways * sin_steer,
                lateral: -ahead * sin_steer + sideways * cos_steer,
            };
        }
        patches
    }

    /// Advance the vehicle by `dt` seconds; see [`PyVehicle::step`].
    pub(super) fn advance(&mut self, dt: f64, throttle: f64, brake: f64, steer: f64) {
        if !(dt.is_finite() && dt > 0.0) {
            return;
        }
        let throttle = limited_input(throttle, 0.0, 1.0);
        let brake = limited_input(brake, 0.0, 1.0);
        self.steering.input = limited_input(steer, -1.0, 1.0);
        let steer_angles = self.wheel_steer_angles();

        // Stability control reacts to the state the wheels are in right now.
        let patches = self.patch_velocities(&steer_angles);
        let slip_ratios: Vec<f64> = self
            .wheels
            .iter()
            .zip(&patches)
            .map(|(wheel, patch)| {
                regularised_slip_ratio(wheel.peripheral_speed(), patch.longitudinal)
            })
            .collect();
        let (forward_speed, _, yaw_rate) = self.planar_velocity();
        let desired_yaw = if self.steering.wheelbase > 0.0 {
            forward_speed * self.steering.front_angle().tan() / self.steering.wheelbase
        } else {
            0.0
        };
        let (throttle_scale, brake_scale) =
            self.stability.update(slip_ratios, desired_yaw, yaw_rate);
        let throttle = throttle * throttle_scale;
        let brake_torque = brake * brake_scale * MAX_BRAKE_TORQUE;

        let loads = self.static_load_array();
        let substeps = ((dt / MAX_SUBSTEP).ceil() as usize).clamp(1, MAX_SUBSTEPS);
        let h = dt / substeps as f64;
        for _ in 0..substeps {
            self.integrate_substep(h, throttle, brake_torque, &steer_angles, &loads);
        }

        self.position[1] += self.velocity[1] * dt;
        self.orientation[2] = wrap_angle(self.orientation[2]);
        self.gear_display = (self.drivetrain.current_gear + 1) as u32;
    }

    /// One sub-step of length `h`: tyre forces, wheel spin and chassis motion.
    fn integrate_substep(
        &mut self,
        h: f64,
        throttle: f64,
        brake_torque: f64,
        steer: &[f64; 4],
        loads: &[f64; 4],
    ) {
        let patches = self.patch_velocities(steer);
        // Rear-wheel drive through an open differential: equal torque split.
        let drive_total = self.drivetrain.wheel_torque(throttle);
        let drive = [0.0, 0.0, 0.5 * drive_total, 0.5 * drive_total];

        let mut force_forward = 0.0;
        let mut force_lateral = 0.0;
        let mut yaw_moment = 0.0;
        for i in 0..4 {
            let patch = patches[i];
            let slip_angle = regularised_slip_angle(patch.lateral, patch.longitudinal);
            let radius = self.wheels[i].radius;
            let omega = WheelSpin {
                tire: &self.tires[i],
                omega: self.wheels[i].omega,
                inertia: self.wheels[i].inertia.max(1.0e-6),
                radius,
                step: h,
                drive_torque: drive[i],
                brake_torque,
                longitudinal_speed: patch.longitudinal,
                slip_angle,
                normal_load: loads[i],
            }
            .solve();
            let slip_ratio = regularised_slip_ratio(omega * radius, patch.longitudinal);
            let (fx, fy) = self.tires[i].compute_forces(slip_angle, slip_ratio, loads[i]);
            let (fx, fy) = (finite_or_zero(fx), finite_or_zero(fy));

            let wheel = &mut self.wheels[i];
            wheel.omega = omega;
            wheel.slip_ratio = slip_ratio;
            wheel.slip_angle = slip_angle;
            wheel.fx = fx;
            wheel.fy = fy;
            wheel.fz = loads[i];
            wheel.grounded = true;

            // Wheel frame to chassis frame, then the moment about the centre of mass.
            let (sin_steer, cos_steer) = steer[i].sin_cos();
            let wheel_forward = fx * cos_steer - fy * sin_steer;
            let wheel_lateral = fx * sin_steer + fy * cos_steer;
            force_forward += wheel_forward;
            force_lateral += wheel_lateral;
            let [ahead, _, aside] = self.wheel_positions[i];
            yaw_moment += ahead * wheel_lateral - aside * wheel_forward;
        }
        // A rolled chassis stands for a banked road: gravity then has a component
        // along the lateral axis.
        force_lateral -= self.mass * GRAVITY * self.orientation[0].sin();

        let (sin_yaw, cos_yaw) = self.orientation[2].sin_cos();
        let ground_speed = self.velocity[0].hypot(self.velocity[2]);
        let force_x = force_forward * cos_yaw
            - force_lateral * sin_yaw
            - DRAG_COEFFICIENT * ground_speed * self.velocity[0];
        let force_z = force_forward * sin_yaw + force_lateral * cos_yaw
            - DRAG_COEFFICIENT * ground_speed * self.velocity[2];

        let inverse_mass = if self.mass.is_finite() && self.mass > 0.0 {
            1.0 / self.mass
        } else {
            0.0
        };
        let yaw_inertia = self.inertia[2];
        let inverse_yaw_inertia = if yaw_inertia.is_finite() && yaw_inertia > 0.0 {
            1.0 / yaw_inertia
        } else {
            0.0
        };
        self.velocity[0] += force_x * inverse_mass * h;
        self.velocity[2] += force_z * inverse_mass * h;
        self.angular_velocity[2] += yaw_moment * inverse_yaw_inertia * h;
        self.orientation[2] += self.angular_velocity[2] * h;
        self.position[0] += self.velocity[0] * h;
        self.position[2] += self.velocity[2] * h;

        let driven_speed = 0.5 * (self.wheels[2].omega + self.wheels[3].omega);
        self.drivetrain.update_rpm(driven_speed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bracketed_root_finds_the_root_of_a_monotone_function() {
        // x + 0.5 sin(x) - 2 is increasing; its root is where x + 0.5 sin(x) = 2.
        let root = bracketed_root(|x| x + 0.5 * x.sin() - 2.0, 0.0, 0.1);
        assert!((root + 0.5 * root.sin() - 2.0).abs() < 1.0e-12, "{root}");
    }

    #[test]
    fn test_bracketed_root_handles_kinks_and_far_roots() {
        // The root of a capped piecewise-linear function with a kink, 4e6 away.
        let f = |x: f64| (x - 4.0e6).max(2.0 * (x - 4.0e6)).min(5.0);
        let root = bracketed_root(f, 0.0, 1.0e-3);
        assert!((root - 4.0e6).abs() < 1.0e-6, "{root}");
        // A bounded force term, as in the wheel spin equation.
        let spin = |w: f64| w - 7.0 - 3.0 * (w / 0.5).clamp(-1.0, 1.0);
        let root = bracketed_root(spin, 7.0, 1.0e-3);
        assert!((root - 10.0).abs() < 1.0e-10, "{root}");
    }

    #[test]
    fn test_bracketed_root_gives_up_on_non_finite_functions() {
        assert_eq!(bracketed_root(|_| f64::NAN, 2.5, 1.0), 2.5);
        // An unusable initial width falls back to a unit bracket.
        assert!((bracketed_root(|x| x - 1.0, 3.0, f64::NAN) - 1.0).abs() < 1.0e-9);
    }

    #[test]
    fn test_regularised_slip_matches_the_classical_definitions_above_the_floor() {
        use super::super::{compute_slip_angle, compute_slip_ratio};
        for (wheel_speed, speed) in [(12.0, 10.0), (8.0, 10.0), (3.0, 3.0), (-6.0, -5.0)] {
            assert!(
                (regularised_slip_ratio(wheel_speed, speed)
                    - compute_slip_ratio(wheel_speed, speed))
                .abs()
                    < 1.0e-15
            );
        }
        for (lateral, longitudinal) in [(1.0, 10.0), (-2.0, 7.5), (0.3, -4.0)] {
            assert!(
                (regularised_slip_angle(lateral, longitudinal)
                    - compute_slip_angle(lateral, longitudinal))
                .abs()
                    < 1.0e-15
            );
        }
    }

    #[test]
    fn test_regularised_slip_stays_finite_and_small_at_standstill() {
        let kappa = regularised_slip_ratio(1.0e-3, 0.0);
        assert!(kappa > 0.0 && kappa < 0.01, "{kappa}");
        let alpha = regularised_slip_angle(1.0e-3, 0.0);
        assert!(alpha < 0.0 && alpha > -0.01, "{alpha}");
        assert_eq!(regularised_slip_ratio(0.0, 0.0), 0.0);
        assert_eq!(regularised_slip_angle(0.0, 0.0), 0.0);
    }

    #[test]
    fn test_wrap_angle_stays_in_a_half_open_turn() {
        use std::f64::consts::PI;
        for angle in [
            0.0,
            1.0,
            -1.0,
            PI,
            -PI,
            3.0 * PI,
            -3.0 * PI + 0.1,
            100.0,
            -100.0,
        ] {
            let wrapped = wrap_angle(angle);
            assert!(
                wrapped > -PI && wrapped <= PI + 1.0e-12,
                "{angle} -> {wrapped}"
            );
            assert!(
                (wrapped.sin() - angle.sin()).abs() < 1.0e-12
                    && (wrapped.cos() - angle.cos()).abs() < 1.0e-12
            );
        }
        assert!(wrap_angle(f64::NAN).is_nan());
    }

    #[test]
    fn test_limited_input_clamps_and_rejects_nan() {
        assert_eq!(limited_input(2.0, 0.0, 1.0), 1.0);
        assert_eq!(limited_input(-2.0, -1.0, 1.0), -1.0);
        assert_eq!(limited_input(f64::INFINITY, 0.0, 1.0), 1.0);
        assert_eq!(limited_input(f64::NAN, 0.0, 1.0), 0.0);
    }
    #[test]
    fn test_wheel_steer_angles_follow_ackermann_rear_steering_and_toe() {
        let mut vehicle = PyVehicle::sedan();
        vehicle.steering.input = 1.0;
        let front = vehicle.steering.front_angle();
        let (inner, outer) = vehicle.steering.ackermann_front_angles();
        // The left wheel is the inner wheel of a left turn.
        let angles = vehicle.wheel_steer_angles();
        assert_eq!(angles[0], inner);
        assert_eq!(angles[1], outer);
        assert!(inner > outer && outer > 0.0);
        assert_eq!(angles[2], 0.0);
        assert_eq!(angles[3], 0.0);
        // A right turn mirrors them.
        vehicle.steering.input = -1.0;
        let angles = vehicle.wheel_steer_angles();
        assert_eq!(angles[0], -outer);
        assert_eq!(angles[1], -inner);
        // Four-wheel steering turns both rear wheels by a fraction of the front angle.
        vehicle.steering.input = 1.0;
        vehicle.steering.four_ws = true;
        vehicle.steering.rear_steer_ratio = -0.5;
        let angles = vehicle.wheel_steer_angles();
        assert!((angles[2] + 0.5 * front).abs() < 1e-15);
        assert!((angles[3] + 0.5 * front).abs() < 1e-15);
        // Toe adds to the heading of its wheel.
        vehicle.wheels[1].toe = 0.01;
        assert!((vehicle.wheel_steer_angles()[1] - (outer + 0.01)).abs() < 1e-15);
    }

    #[test]
    fn test_planar_velocity_rotates_the_world_velocity_into_the_chassis_frame() {
        let mut vehicle = PyVehicle::sedan();
        vehicle.orientation[2] = 0.3;
        vehicle.velocity = [3.0, 0.5, 4.0];
        vehicle.angular_velocity[2] = 0.2;
        let (forward, lateral, yaw_rate) = vehicle.planar_velocity();
        assert!((forward - (3.0 * 0.3_f64.cos() + 4.0 * 0.3_f64.sin())).abs() < 1e-14);
        assert!((lateral - (-3.0 * 0.3_f64.sin() + 4.0 * 0.3_f64.cos())).abs() < 1e-14);
        assert_eq!(yaw_rate, 0.2);
        assert!((forward.hypot(lateral) - 5.0).abs() < 1e-14);
        assert_eq!(vehicle.forward_speed(), forward);
        assert_eq!(vehicle.lateral_speed(), lateral);
    }

    #[test]
    fn test_patch_velocity_includes_the_rotation_of_the_chassis() {
        let mut vehicle = PyVehicle::sedan();
        // Pure yaw: the left front wheel (1.3 ahead, 0.75 to the left) moves
        // backwards and sideways, the rear left wheel backwards and to the other side.
        vehicle.angular_velocity[2] = 1.0;
        let patches = vehicle.patch_velocities(&[0.0; 4]);
        assert!((patches[0].longitudinal + 0.75).abs() < 1e-14);
        assert!((patches[0].lateral - 1.3).abs() < 1e-14);
        assert!((patches[1].longitudinal - 0.75).abs() < 1e-14);
        assert!((patches[2].lateral + 1.3).abs() < 1e-14);
        // A steered wheel sees the ground velocity in its own frame.
        vehicle.angular_velocity[2] = 0.0;
        vehicle.velocity = [10.0, 0.0, 0.0];
        let steer = 0.2_f64;
        let patch = vehicle.patch_velocities(&[steer, 0.0, 0.0, 0.0])[0];
        assert!((patch.longitudinal - 10.0 * steer.cos()).abs() < 1e-12);
        assert!((patch.lateral + 10.0 * steer.sin()).abs() < 1e-12);
        // ... which is a positive slip angle for a wheel steered toward the lateral axis.
        assert!(regularised_slip_angle(patch.lateral, patch.longitudinal) > 0.0);
    }

    #[test]
    fn test_wheel_spin_solution_satisfies_the_implicit_equation() {
        let vehicle = PyVehicle::sedan();
        let spin = WheelSpin {
            tire: &vehicle.tires[2],
            omega: 19.0,
            inertia: 1.5,
            radius: 0.32,
            step: 1.0e-3,
            drive_torque: 900.0,
            brake_torque: 0.0,
            longitudinal_speed: 6.0,
            slip_angle: 0.0,
            normal_load: 3678.75,
        };
        let omega = spin.solve();
        assert!(
            spin.residual(omega).abs() < 1.0e-9,
            "residual {}",
            spin.residual(omega)
        );
        // Driven faster than the ground passes, the wheel speeds up, by less than
        // the free acceleration because the road pushes back.
        assert!(
            omega > 19.0 && omega < 19.0 + 900.0 / 1.5 * 1.0e-3,
            "omega {omega}"
        );
        // A locked brake on a wheel that is turning pulls it to a stop quickly ...
        let braked = WheelSpin {
            brake_torque: 2000.0,
            drive_torque: 0.0,
            ..spin
        };
        let omega = braked.solve();
        assert!(braked.residual(omega).abs() < 1.0e-9);
        assert!(omega < 19.0);
        // ... and never past zero, however large the step.
        let long_step = WheelSpin {
            step: 10.0,
            ..braked
        };
        let omega = long_step.solve();
        assert!(omega > -1.0e-6 && omega < 1.0, "omega {omega}");
    }
}
