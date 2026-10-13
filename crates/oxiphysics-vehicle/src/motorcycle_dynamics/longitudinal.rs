// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Longitudinal wheel/ground coupling of the single-track motorcycle model.
//!
//! The forward speed `v` of the bike and the spin of its two wheels interact
//! through the longitudinal tyre forces `F_x` only (positive = forward on the
//! bike):
//!
//! * `m v̇ = F_x,r + F_x,f − F_aero`
//! * `I_r ω̇_r = T_drive − T_brake,r − r_r F_x,r`
//! * `I_f ω̇_f = −T_brake,f − r_f F_x,f`
//!
//! The tyre force is the Magic Formula force at the slip ratio
//! `κ = (ω r − v) / max(|v|, v_floor)`.  It is the *signed* force that reacts on
//! the wheel: a driving force (`F_x > 0`) retards the wheel, a braking force
//! (`F_x < 0`) spins it up.  The drive and brake torques therefore reach the bike
//! only through the slip they create; adding them to the chassis force as well
//! would count them twice.
//!
//! A friction brake opposes the rotation of its wheel.  It can stop a wheel and
//! hold it stopped, but it never reverses it: when the brake outweighs the drive
//! torque and the ground reaction of the locked wheel (the wheel is *held*), the
//! wheel stays at zero speed and the bike sees the force of a locked tyre.
//!
//! A tyre produces `C_κ ≈ 20 F_z` newtons per unit slip ratio and a unit slip
//! ratio is only `|v|` of slip speed, so the wheel and chassis rates are stiff,
//! increasingly so at low speed.  An explicit Euler step is unstable below roughly
//! 10 m/s at ordinary step sizes.  [`advance`] therefore takes a linearly implicit
//! (Rosenbrock-Euler) step, `(I − Δt J) Δx = Δt f(x)`, with `J` the Jacobian of the
//! rates `f` for `x = (v, ω_r, ω_f)`.  The tyre slopes enter `J` clamped to be
//! non-negative, which keeps the coefficient matrix well conditioned for every
//! `Δt`.

use super::PacejkaMotoTire;

/// Lowest speed (m/s) used as the denominator of the longitudinal slip ratio.
///
/// `(ω r − v) / |v|` is undefined at standstill, and its sensitivity to the wheel
/// speed grows without bound as the bike slows down.  With a hard "no slip below
/// some speed" cut-off the bike could never pull away, because no slip means no
/// tyre force; with the bare formula the wheel/ground coupling becomes infinitely
/// stiff.  Flooring the denominator keeps the coupling finite.  Above this speed
/// the slip ratio is the textbook `(ω r − v) / |v|`.
pub(super) const SLIP_SPEED_FLOOR: f64 = 1.0;

/// Speed (m/s) under which a bike that is still being slowed down comes to rest.
///
/// The force of a locked wheel fades out linearly with the speed below
/// [`SLIP_SPEED_FLOOR`], so a braking bike would otherwise only approach
/// standstill asymptotically.
pub(super) const REST_SPEED: f64 = 1.0e-3;

/// Longitudinal slip ratio `(ω r − v) / max(|v|, SLIP_SPEED_FLOOR)` of a wheel of
/// radius `radius` spinning at `omega` on ground moving at `speed`.
pub(super) fn slip_ratio(omega: f64, radius: f64, speed: f64) -> f64 {
    (omega * radius - speed) / speed.abs().max(SLIP_SPEED_FLOOR)
}

/// Partial derivatives `(∂κ/∂ω, ∂κ/∂v)` of [`slip_ratio`].
pub(super) fn slip_ratio_derivatives(omega: f64, radius: f64, speed: f64) -> (f64, f64) {
    let d_omega = radius / speed.abs().max(SLIP_SPEED_FLOOR);
    let d_speed = if speed.abs() > SLIP_SPEED_FLOOR {
        -omega * radius / (speed * speed.abs())
    } else {
        -1.0 / SLIP_SPEED_FLOOR
    };
    (d_omega, d_speed)
}

/// Slope `∂F_x/∂κ` (N per unit slip ratio) of the combined-slip longitudinal tyre
/// force, from a central difference of [`PacejkaMotoTire::forces`].
///
/// Beyond the force peak the slope is negative; that branch is physically
/// unstable and must not enter an implicit step, so it is reported as zero.
pub(super) fn longitudinal_slope(
    tire: &PacejkaMotoTire,
    slip_angle: f64,
    slip: f64,
    normal_force: f64,
    lean_angle: f64,
) -> f64 {
    const STEP: f64 = 1.0e-3;
    let force = |kappa: f64| tire.forces(slip_angle, kappa, normal_force, lean_angle).1;
    ((force(slip + STEP) - force(slip - STEP)) / (2.0 * STEP)).max(0.0)
}

/// Longitudinal tyre force at the start of a step, with its slope, for the wheel
/// as it rolls and for the wheel locked.
#[derive(Debug, Clone, Copy)]
pub(super) struct TyreLongitudinal {
    /// Force on the bike (N), positive = forward, at the current slip ratio.
    pub force: f64,
    /// `∂F_x/∂κ` at the current slip ratio (N), never negative.
    pub slope: f64,
    /// Force on the bike (N) when the wheel is locked (`ω = 0`).
    pub locked_force: f64,
    /// `∂F_x/∂κ` of the locked wheel (N), never negative.
    pub locked_slope: f64,
}

impl TyreLongitudinal {
    /// Evaluate `tire` for a wheel of radius `radius` spinning at `omega` on ground
    /// moving at `speed`.
    pub(super) fn evaluate(
        tire: &PacejkaMotoTire,
        slip_angle: f64,
        normal_force: f64,
        lean_angle: f64,
        (omega, radius): (f64, f64),
        speed: f64,
    ) -> Self {
        let slip = slip_ratio(omega, radius, speed);
        let locked_slip = slip_ratio(0.0, radius, speed);
        let force = |kappa: f64| tire.forces(slip_angle, kappa, normal_force, lean_angle).1;
        Self {
            force: force(slip),
            slope: longitudinal_slope(tire, slip_angle, slip, normal_force, lean_angle),
            locked_force: force(locked_slip),
            locked_slope: longitudinal_slope(
                tire,
                slip_angle,
                locked_slip,
                normal_force,
                lean_angle,
            ),
        }
    }
}

/// One wheel as seen by [`advance`].
#[derive(Debug, Clone, Copy)]
pub(super) struct WheelInput {
    /// Rolling radius (m).
    pub radius: f64,
    /// Moment of inertia about the axle (kg·m²).
    pub inertia: f64,
    /// Angular velocity (rad/s), never negative.
    pub omega: f64,
    /// Torque that drives the wheel forward (N·m).
    pub drive_torque: f64,
    /// Torque of the friction brake (N·m, ≥ 0).
    pub brake_torque: f64,
    /// Tyre force and slopes at the start of the step.
    pub tyre: TyreLongitudinal,
}

impl WheelInput {
    /// Whether the friction brake holds this wheel at standstill during a step of
    /// `dt` seconds.
    ///
    /// That is the case when the brake could stop the wheel within the step on its
    /// own (`ω ≤ Δt T_brake / I`) and its torque outweighs what turns the stopped
    /// wheel forward: the drive torque plus the ground reaction of the locked tyre.
    fn is_held(&self, dt: f64) -> bool {
        self.brake_torque > 0.0
            && dt > 0.0
            && self.omega <= dt * self.brake_torque / self.inertia
            && self.drive_torque - self.radius * self.tyre.locked_force <= self.brake_torque
    }
}

/// The bike as seen by [`advance`].
#[derive(Debug, Clone, Copy)]
pub(super) struct ChassisInput {
    /// Mass (kg).
    pub mass: f64,
    /// Forward speed (m/s), never negative.
    pub speed: f64,
    /// Aerodynamic drag at `speed` (N, ≥ 0); it grows with `speed²`.
    pub drag: f64,
}

/// Forward speed and wheel speeds of the bike.
#[derive(Debug, Clone, Copy)]
pub(super) struct LongitudinalState {
    /// Forward speed (m/s).
    pub speed: f64,
    /// Rear wheel angular velocity (rad/s).
    pub rear_omega: f64,
    /// Front wheel angular velocity (rad/s).
    pub front_omega: f64,
}

/// A wheel linearised about its state at the start of a step.
#[derive(Debug, Clone, Copy)]
struct WheelRates {
    /// Force on the bike (N).
    force: f64,
    /// `∂F/∂ω` (N per rad/s).
    d_force_d_omega: f64,
    /// `∂F/∂v` (N per m/s).
    d_force_d_speed: f64,
    /// Angular acceleration (rad/s²); zero for a held wheel.
    omega_rate: f64,
}

impl WheelRates {
    fn new(wheel: &WheelInput, held: bool, speed: f64) -> Self {
        if held {
            let (_, d_speed) = slip_ratio_derivatives(0.0, wheel.radius, speed);
            Self {
                force: wheel.tyre.locked_force,
                d_force_d_omega: 0.0,
                d_force_d_speed: wheel.tyre.locked_slope * d_speed,
                omega_rate: 0.0,
            }
        } else {
            let (d_omega, d_speed) = slip_ratio_derivatives(wheel.omega, wheel.radius, speed);
            let torque = wheel.drive_torque - wheel.brake_torque - wheel.radius * wheel.tyre.force;
            Self {
                force: wheel.tyre.force,
                d_force_d_omega: wheel.tyre.slope * d_omega,
                d_force_d_speed: wheel.tyre.slope * d_speed,
                omega_rate: torque / wheel.inertia,
            }
        }
    }
}

/// Advance forward speed and wheel speeds by `dt` seconds.
///
/// Returns the new state and the longitudinal acceleration of the bike (m/s²), see
/// the module documentation for the model and the integration scheme.  A bike that
/// is still being slowed down and ends the step slower than [`REST_SPEED`] comes to
/// rest.
pub(super) fn advance(
    dt: f64,
    chassis: ChassisInput,
    rear: WheelInput,
    front: WheelInput,
) -> (LongitudinalState, f64) {
    let v = chassis.speed;
    let rear_held = rear.is_held(dt);
    let front_held = front.is_held(dt);
    let r = WheelRates::new(&rear, rear_held, v);
    let f = WheelRates::new(&front, front_held, v);

    // Aerodynamic drag D = k v² has the slope dD/dv = 2 D / v.
    let drag_slope = if v.abs() > 1.0e-9 {
        2.0 * chassis.drag / v
    } else {
        0.0
    };

    // Rates f(x) and Jacobian J = ∂f/∂x of the speed (v) and wheel (r, f) rows.
    // The wheels couple to each other only through the speed.
    let f_v = (r.force + f.force - chassis.drag) / chassis.mass;
    let j_vv = (r.d_force_d_speed + f.d_force_d_speed - drag_slope) / chassis.mass;
    let j_vr = r.d_force_d_omega / chassis.mass;
    let j_vf = f.d_force_d_omega / chassis.mass;
    let j_rr = -rear.radius * r.d_force_d_omega / rear.inertia;
    let j_rv = -rear.radius * r.d_force_d_speed / rear.inertia;
    let j_ff = -front.radius * f.d_force_d_omega / front.inertia;
    let j_fv = -front.radius * f.d_force_d_speed / front.inertia;

    // Solve (I − Δt J) Δx = Δt f by eliminating the two wheel rows.  Both pivots
    // are at least 1, and so is the speed denominator (the wheel terms it collects
    // are non-negative), so no division can fail.
    let pivot_r = 1.0 - dt * j_rr;
    let pivot_f = 1.0 - dt * j_ff;
    let speed_denominator =
        1.0 - dt * j_vv - dt * dt * (j_vr * j_rv / pivot_r + j_vf * j_fv / pivot_f);
    let speed_numerator =
        dt * f_v + dt * dt * (j_vr * r.omega_rate / pivot_r + j_vf * f.omega_rate / pivot_f);
    let d_v = speed_numerator / speed_denominator;
    let d_omega_r = dt * (r.omega_rate + j_rv * d_v) / pivot_r;
    let d_omega_f = dt * (f.omega_rate + j_fv * d_v) / pivot_f;

    let mut speed = (v + d_v).max(0.0);
    if speed < REST_SPEED && d_v <= 0.0 {
        speed = 0.0;
    }
    let state = LongitudinalState {
        speed,
        rear_omega: if rear_held {
            0.0
        } else {
            (rear.omega + d_omega_r).max(0.0)
        },
        front_omega: if front_held {
            0.0
        } else {
            (front.omega + d_omega_f).max(0.0)
        },
    };
    let accel = if dt > 0.0 { (speed - v) / dt } else { f_v };
    (state, accel)
}

#[cfg(test)]
mod tests;
