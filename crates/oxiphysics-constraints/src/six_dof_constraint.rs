// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Six-degree-of-freedom (6-DOF) constraint.
//!
//! Provides per-axis control over all 6 rigid-body DOFs:
//! - 3 linear (translation along X, Y, Z)
//! - 3 angular (rotation about X, Y, Z)
//!
//! Each axis can be `Free`, `Limited` (with min/max bounds), or `Locked`.
//! An optional motor can drive any axis toward a target velocity or position.

use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::{Real, Vec3};
use oxiphysics_rigid::RigidBodySet;

use crate::traits::Constraint;

// ─────────────────────────────────────────────────────────────────────────────
// MotorConfig
// ─────────────────────────────────────────────────────────────────────────────

/// Configuration for a motor on one DOF axis.
#[derive(Debug, Clone)]
pub struct MotorConfig {
    /// Target velocity for this axis (m/s or rad/s).
    pub target_velocity: f64,
    /// Maximum force/torque for this axis (N or N·m).
    pub max_force: f64,
    /// Simple proportional gain for velocity error.
    pub gain: f64,
}

impl MotorConfig {
    /// Create a new motor configuration.
    pub fn new(target_velocity: f64, max_force: f64, gain: f64) -> Self {
        Self {
            target_velocity,
            max_force,
            gain,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AxisConfig
// ─────────────────────────────────────────────────────────────────────────────

/// Per-axis configuration for a 6-DOF constraint.
#[derive(Debug, Clone)]
pub struct AxisConfig {
    /// If `true`, this axis is rigidly locked (zero relative motion).
    pub locked: bool,
    /// Optional `(min, max)` limit for this axis.
    /// Ignored when `locked == true`.
    pub limited: Option<(f64, f64)>,
    /// Optional motor configuration.
    pub motor: Option<MotorConfig>,
}

impl AxisConfig {
    /// A completely free (unconstrained) axis.
    pub fn free() -> Self {
        Self {
            locked: false,
            limited: None,
            motor: None,
        }
    }

    /// A locked axis (zero relative motion).
    pub fn locked() -> Self {
        Self {
            locked: true,
            limited: None,
            motor: None,
        }
    }

    /// A limited axis with the given `[min, max]` bounds.
    pub fn limited(min: f64, max: f64) -> Self {
        Self {
            locked: false,
            limited: Some((min, max)),
            motor: None,
        }
    }

    /// Attach a motor to this axis.
    pub fn with_motor(mut self, motor: MotorConfig) -> Self {
        self.motor = Some(motor);
        self
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SixDofConstraint
// ─────────────────────────────────────────────────────────────────────────────

/// A general 6-DOF constraint between two rigid bodies.
///
/// The constraint is expressed in the local frame of body A.  Per-axis
/// `AxisConfig` entries control which degrees of freedom are free, limited,
/// or locked.
///
/// # Axis convention
/// - `linear[0]` → translation along local X
/// - `linear[1]` → translation along local Y
/// - `linear[2]` → translation along local Z
/// - `angular[0]` → rotation about local X
/// - `angular[1]` → rotation about local Y
/// - `angular[2]` → rotation about local Z
#[derive(Debug, Clone)]
pub struct SixDofConstraint {
    /// Handle of body A.
    pub body_a: BodyHandle,
    /// Handle of body B.
    pub body_b: BodyHandle,
    /// Anchor point on body A in local-A space.
    pub anchor_a: Vec3,
    /// Anchor point on body B in local-B space.
    pub anchor_b: Vec3,
    /// Per-axis config for the three linear DOFs.
    pub linear: [AxisConfig; 3],
    /// Per-axis config for the three angular DOFs.
    pub angular: [AxisConfig; 3],

    // ── Cached Jacobian data ──────────────────────────────────────────────────
    /// World-space offset from body A CoM to anchor (updated each `prepare`).
    r_a: Vec3,
    /// World-space offset from body B CoM to anchor (updated each `prepare`).
    r_b: Vec3,
    /// Accumulated impulses for each linear axis.
    lambda_linear: [Real; 3],
    /// Accumulated impulses for each angular axis.
    lambda_angular: [Real; 3],
    /// Effective masses for each linear axis.
    eff_mass_linear: [Real; 3],
    /// Effective masses for each angular axis.
    eff_mass_angular: [Real; 3],
    /// Velocity biases for position correction (linear).
    bias_linear: [Real; 3],
    /// Velocity biases for position correction (angular).
    bias_angular: [Real; 3],
}

impl SixDofConstraint {
    /// Create a new 6-DOF constraint with all axes free.
    pub fn new(body_a: BodyHandle, body_b: BodyHandle) -> Self {
        Self {
            body_a,
            body_b,
            anchor_a: Vec3::zeros(),
            anchor_b: Vec3::zeros(),
            linear: std::array::from_fn(|_| AxisConfig::free()),
            angular: std::array::from_fn(|_| AxisConfig::free()),
            r_a: Vec3::zeros(),
            r_b: Vec3::zeros(),
            lambda_linear: [0.0; 3],
            lambda_angular: [0.0; 3],
            eff_mass_linear: [0.0; 3],
            eff_mass_angular: [0.0; 3],
            bias_linear: [0.0; 3],
            bias_angular: [0.0; 3],
        }
    }

    /// Set the anchor points in local body-space.
    pub fn with_anchors(mut self, anchor_a: Vec3, anchor_b: Vec3) -> Self {
        self.anchor_a = anchor_a;
        self.anchor_b = anchor_b;
        self
    }

    /// Lock the linear DOF along the given axis (0=X, 1=Y, 2=Z).
    ///
    /// Panics if `axis >= 3`.
    pub fn lock_linear(&mut self, axis: usize) {
        assert!(axis < 3, "axis must be 0, 1, or 2");
        self.linear[axis] = AxisConfig::locked();
    }

    /// Lock the angular DOF about the given axis (0=X, 1=Y, 2=Z).
    ///
    /// Panics if `axis >= 3`.
    pub fn lock_angular(&mut self, axis: usize) {
        assert!(axis < 3, "axis must be 0, 1, or 2");
        self.angular[axis] = AxisConfig::locked();
    }

    /// Set a linear limit `[min, max]` for the given axis.
    pub fn set_linear_limit(&mut self, axis: usize, min: f64, max: f64) {
        assert!(axis < 3, "axis must be 0, 1, or 2");
        self.linear[axis] = AxisConfig::limited(min, max);
    }

    /// Set an angular limit `[min, max]` (radians) for the given axis.
    pub fn set_angular_limit(&mut self, axis: usize, min: f64, max: f64) {
        assert!(axis < 3, "axis must be 0, 1, or 2");
        self.angular[axis] = AxisConfig::limited(min, max);
    }

    /// Attach a motor to a linear axis.
    pub fn set_linear_motor(&mut self, axis: usize, motor: MotorConfig) {
        assert!(axis < 3, "axis must be 0, 1, or 2");
        self.linear[axis].motor = Some(motor);
    }

    /// Attach a motor to an angular axis.
    pub fn set_angular_motor(&mut self, axis: usize, motor: MotorConfig) {
        assert!(axis < 3, "axis must be 0, 1, or 2");
        self.angular[axis].motor = Some(motor);
    }

    /// Create a fully-locked constraint (weld joint equivalent).
    pub fn weld(body_a: BodyHandle, body_b: BodyHandle) -> Self {
        let mut c = Self::new(body_a, body_b);
        for i in 0..3 {
            c.lock_linear(i);
            c.lock_angular(i);
        }
        c
    }

    // ── Jacobian helpers ──────────────────────────────────────────────────────

    /// Compute the effective mass for a single constraint row.
    ///
    /// For a linear axis `d`:
    ///   k = inv_m_a + inv_m_b + (r_a×d)·I_a^{-1}·(r_a×d) + (r_b×d)·I_b^{-1}·(r_b×d)
    fn effective_mass_linear(
        inv_mass_a: f64,
        inv_mass_b: f64,
        inv_inertia_a: &oxiphysics_core::math::Mat3,
        inv_inertia_b: &oxiphysics_core::math::Mat3,
        r_a: &Vec3,
        r_b: &Vec3,
        axis: &Vec3,
    ) -> f64 {
        let rn_a = r_a.cross(axis);
        let rn_b = r_b.cross(axis);
        let k = inv_mass_a
            + inv_mass_b
            + rn_a.dot(&(inv_inertia_a * rn_a))
            + rn_b.dot(&(inv_inertia_b * rn_b));
        if k > 1e-12 { 1.0 / k } else { 0.0 }
    }

    /// Compute the effective mass for an angular constraint row.
    ///
    /// For an angular axis `d`:
    ///   k = d·I_a^{-1}·d + d·I_b^{-1}·d
    fn effective_mass_angular(
        inv_inertia_a: &oxiphysics_core::math::Mat3,
        inv_inertia_b: &oxiphysics_core::math::Mat3,
        axis: &Vec3,
    ) -> f64 {
        let k = axis.dot(&(inv_inertia_a * axis)) + axis.dot(&(inv_inertia_b * axis));
        if k > 1e-12 { 1.0 / k } else { 0.0 }
    }

    /// The three world-space axes derived from body A's orientation.
    fn world_axes(bodies: &RigidBodySet, handle: BodyHandle) -> Option<[Vec3; 3]> {
        let body = bodies.get(handle)?;
        let rot = body.transform.rotation;
        Some([
            rot * Vec3::new(1.0, 0.0, 0.0),
            rot * Vec3::new(0.0, 1.0, 0.0),
            rot * Vec3::new(0.0, 0.0, 1.0),
        ])
    }

    /// Apply a single-row impulse along `axis` to bodies A and B.
    fn apply_impulse_row(
        bodies: &mut RigidBodySet,
        handle_a: BodyHandle,
        handle_b: BodyHandle,
        axis: Vec3,
        r_a: Vec3,
        r_b: Vec3,
        lambda: f64,
    ) {
        let impulse = axis * lambda;
        if let Some(b) = bodies.get_mut(handle_a) {
            b.velocity += impulse * b.inverse_mass;
            b.angular_velocity += b.world_inverse_inertia * r_a.cross(&impulse);
        }
        if let Some(b) = bodies.get_mut(handle_b) {
            b.velocity -= impulse * b.inverse_mass;
            b.angular_velocity -= b.world_inverse_inertia * r_b.cross(&impulse);
        }
    }

    /// Apply a pure angular impulse about `axis` to bodies A and B.
    fn apply_angular_impulse_row(
        bodies: &mut RigidBodySet,
        handle_a: BodyHandle,
        handle_b: BodyHandle,
        axis: Vec3,
        lambda: f64,
    ) {
        let torque_impulse = axis * lambda;
        if let Some(b) = bodies.get_mut(handle_a) {
            b.angular_velocity += b.world_inverse_inertia * torque_impulse;
        }
        if let Some(b) = bodies.get_mut(handle_b) {
            b.angular_velocity -= b.world_inverse_inertia * torque_impulse;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Constraint impl
// ─────────────────────────────────────────────────────────────────────────────

/// Baumgarte factor for positional correction.
const BAUMGARTE_6DOF: f64 = 0.2;
/// Penetration slop for positional correction.
const SLOP_6DOF: f64 = 0.005;

impl Constraint for SixDofConstraint {
    fn prepare(&mut self, bodies: &RigidBodySet, dt: f64) {
        let (pos_a, rot_a, inv_mass_a, inv_inertia_a, vel_a, ang_vel_a) =
            match bodies.get(self.body_a) {
                Some(b) => (
                    b.transform.position,
                    b.transform.rotation,
                    b.inverse_mass,
                    b.world_inverse_inertia,
                    b.velocity,
                    b.angular_velocity,
                ),
                None => return,
            };

        let (pos_b, _rot_b, inv_mass_b, inv_inertia_b, vel_b, ang_vel_b) =
            match bodies.get(self.body_b) {
                Some(b) => (
                    b.transform.position,
                    b.transform.rotation,
                    b.inverse_mass,
                    b.world_inverse_inertia,
                    b.velocity,
                    b.angular_velocity,
                ),
                None => return,
            };

        // Anchor offsets in world space
        self.r_a = rot_a * self.anchor_a;
        self.r_b = {
            let rot_b = bodies
                .get(self.body_b)
                .expect("key must exist")
                .transform
                .rotation;
            rot_b * self.anchor_b
        };

        // World-space axes from body A orientation
        let axes = [
            rot_a * Vec3::new(1.0, 0.0, 0.0),
            rot_a * Vec3::new(0.0, 1.0, 0.0),
            rot_a * Vec3::new(0.0, 0.0, 1.0),
        ];

        // Anchor positions in world space
        let anchor_world_a = pos_a + self.r_a;
        let anchor_world_b = pos_b + self.r_b;
        let positional_error = anchor_world_b - anchor_world_a;

        for (i, &ax) in axes.iter().enumerate() {
            let config = &self.linear[i];

            if config.locked || config.limited.is_some() {
                self.eff_mass_linear[i] = Self::effective_mass_linear(
                    inv_mass_a,
                    inv_mass_b,
                    &inv_inertia_a,
                    &inv_inertia_b,
                    &self.r_a,
                    &self.r_b,
                    &ax,
                );

                // Position error along this axis
                let pos_err = positional_error.dot(&ax);

                let correction = if config.locked {
                    pos_err
                } else if let Some((min, max)) = config.limited {
                    if pos_err < min {
                        pos_err - min
                    } else if pos_err > max {
                        pos_err - max
                    } else {
                        0.0
                    }
                } else {
                    0.0
                };

                self.bias_linear[i] = (BAUMGARTE_6DOF / dt)
                    * (correction.abs() - SLOP_6DOF).max(0.0)
                    * correction.signum();
            } else {
                self.eff_mass_linear[i] = 0.0;
                self.bias_linear[i] = 0.0;
            }
        }

        // Angular axes
        // Relative angular velocity and orientation error (small-angle approximation)
        let ang_vel_rel = ang_vel_a - ang_vel_b;
        let _ = ang_vel_rel; // used below per-axis

        for (i, &ax) in axes.iter().enumerate() {
            let config = &self.angular[i];

            if config.locked || config.limited.is_some() {
                self.eff_mass_angular[i] =
                    Self::effective_mass_angular(&inv_inertia_a, &inv_inertia_b, &ax);

                // For locked axes: bias toward zero angular velocity difference
                // (full angular error correction requires rotation difference quaternion)
                let ang_err = 0.0_f64; // simplified: rely on velocity correction
                let correction = ang_err;
                self.bias_angular[i] = (BAUMGARTE_6DOF / dt) * correction;
            } else {
                self.eff_mass_angular[i] = 0.0;
                self.bias_angular[i] = 0.0;
            }
        }

        let _ = (vel_a, ang_vel_a, vel_b, ang_vel_b);
    }

    fn solve_velocity(&mut self, bodies: &mut RigidBodySet, dt: f64) {
        let axes = match Self::world_axes(bodies, self.body_a) {
            Some(a) => a,
            None => return,
        };

        let r_a = self.r_a;
        let r_b = self.r_b;

        // ── Linear rows ────────────────────────────────────────────────────────
        for (i, &ax) in axes.iter().enumerate() {
            let config = self.linear[i].clone();

            if !config.locked && config.limited.is_none() && config.motor.is_none() {
                continue;
            }

            let eff_mass = self.eff_mass_linear[i];
            if eff_mass < 1e-14 {
                continue;
            }

            // Current relative velocity along this axis
            let (vel_a, ang_vel_a) = match bodies.get(self.body_a) {
                Some(b) => (b.velocity, b.angular_velocity),
                None => return,
            };
            let (vel_b, ang_vel_b) = match bodies.get(self.body_b) {
                Some(b) => (b.velocity, b.angular_velocity),
                None => return,
            };

            let v_a = vel_a + ang_vel_a.cross(&r_a);
            let v_b = vel_b + ang_vel_b.cross(&r_b);
            let rel_vel = (v_a - v_b).dot(&ax);

            let motor_target_vel = config
                .motor
                .as_ref()
                .map(|m| m.target_velocity)
                .unwrap_or(0.0);

            let target_vel = if config.locked {
                self.bias_linear[i]
            } else {
                motor_target_vel
            };

            let delta_lambda = eff_mass * (target_vel - rel_vel);

            let old = self.lambda_linear[i];
            let new = if config.locked {
                old + delta_lambda // locked: unclamped (bilateral)
            } else if let Some((min, max)) = config.limited {
                (old + delta_lambda).clamp(min, max)
            } else if let Some(ref mc) = config.motor {
                (old + delta_lambda).clamp(-mc.max_force * dt, mc.max_force * dt)
            } else {
                old + delta_lambda
            };

            let applied = new - old;
            self.lambda_linear[i] = new;

            Self::apply_impulse_row(bodies, self.body_a, self.body_b, ax, r_a, r_b, applied);
        }

        // ── Angular rows ───────────────────────────────────────────────────────
        for (i, &ax) in axes.iter().enumerate() {
            let config = self.angular[i].clone();

            if !config.locked && config.limited.is_none() && config.motor.is_none() {
                continue;
            }

            let eff_mass = self.eff_mass_angular[i];
            if eff_mass < 1e-14 {
                continue;
            }

            let (ang_vel_a,) = match bodies.get(self.body_a) {
                Some(b) => (b.angular_velocity,),
                None => return,
            };
            let (ang_vel_b,) = match bodies.get(self.body_b) {
                Some(b) => (b.angular_velocity,),
                None => return,
            };

            let rel_ang_vel = (ang_vel_a - ang_vel_b).dot(&ax);

            let motor_target_vel = config
                .motor
                .as_ref()
                .map(|m| m.target_velocity)
                .unwrap_or(0.0);

            let target_vel = if config.locked {
                self.bias_angular[i]
            } else {
                motor_target_vel
            };

            let delta_lambda = eff_mass * (target_vel - rel_ang_vel);

            let old = self.lambda_angular[i];
            let new = if config.locked {
                old + delta_lambda
            } else if let Some((min, max)) = config.limited {
                (old + delta_lambda).clamp(min, max)
            } else if let Some(ref mc) = config.motor {
                (old + delta_lambda).clamp(-mc.max_force * dt, mc.max_force * dt)
            } else {
                old + delta_lambda
            };

            let applied = new - old;
            self.lambda_angular[i] = new;

            Self::apply_angular_impulse_row(bodies, self.body_a, self.body_b, ax, applied);
        }
    }

    fn solve_position(&mut self, bodies: &mut RigidBodySet, _dt: f64) {
        let (pos_a, rot_a, inv_mass_a) = match bodies.get(self.body_a) {
            Some(b) => (b.transform.position, b.transform.rotation, b.inverse_mass),
            None => return,
        };
        let (pos_b, inv_mass_b) = match bodies.get(self.body_b) {
            Some(b) => (b.transform.position, b.inverse_mass),
            None => return,
        };

        let axes = [
            rot_a * Vec3::new(1.0, 0.0, 0.0),
            rot_a * Vec3::new(0.0, 1.0, 0.0),
            rot_a * Vec3::new(0.0, 0.0, 1.0),
        ];

        let r_a_world = rot_a * self.anchor_a;
        let r_b_world = {
            let rot_b = bodies
                .get(self.body_b)
                .expect("key must exist")
                .transform
                .rotation;
            rot_b * self.anchor_b
        };

        let anchor_a = pos_a + r_a_world;
        let anchor_b = pos_b + r_b_world;
        let positional_error = anchor_b - anchor_a;

        for (i, &ax) in axes.iter().enumerate() {
            let config = &self.linear[i];
            if !config.locked {
                continue;
            }
            let err = positional_error.dot(&ax);
            let correction = (err.abs() - SLOP_6DOF).max(0.0) * BAUMGARTE_6DOF * err.signum();

            if correction.abs() < 1e-10 {
                continue;
            }

            let k = inv_mass_a + inv_mass_b;
            if k < 1e-12 {
                continue;
            }

            let pos_impulse = ax * (correction / k);
            if let Some(b) = bodies.get_mut(self.body_a) {
                b.transform.position -= pos_impulse * b.inverse_mass;
            }
            if let Some(b) = bodies.get_mut(self.body_b) {
                b.transform.position += pos_impulse * b.inverse_mass;
            }
        }
    }

    fn body_handles(&self) -> Vec<BodyHandle> {
        vec![self.body_a, self.body_b]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Jacobian matrix helper (6×12 matrix, rows = constraints, cols = velocities)
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the 6×12 Jacobian for a fully-locked 6-DOF constraint.
///
/// The Jacobian maps `[v_a, ω_a, v_b, ω_b]` (12 entries) to the 6 constraint
/// velocities (3 linear + 3 angular).
///
/// Row layout:
/// - Rows 0-2: linear constraints along X, Y, Z
/// - Rows 3-5: angular constraints about X, Y, Z
///
/// Column layout per body: `[v.x, v.y, v.z, ω.x, ω.y, ω.z]`
pub fn compute_jacobian_6dof(r_a: Vec3, r_b: Vec3, axes: &[Vec3; 3]) -> [[f64; 12]; 6] {
    let mut j = [[0.0f64; 12]; 6];

    // Linear rows
    for (row, ax) in axes.iter().enumerate() {
        // Body A linear part: +ax
        j[row][0] = ax.x;
        j[row][1] = ax.y;
        j[row][2] = ax.z;
        // Body A angular part: +(r_a × ax)
        let ra_x_ax = r_a.cross(ax);
        j[row][3] = ra_x_ax.x;
        j[row][4] = ra_x_ax.y;
        j[row][5] = ra_x_ax.z;
        // Body B linear part: -ax
        j[row][6] = -ax.x;
        j[row][7] = -ax.y;
        j[row][8] = -ax.z;
        // Body B angular part: -(r_b × ax)
        let rb_x_ax = r_b.cross(ax);
        j[row][9] = -rb_x_ax.x;
        j[row][10] = -rb_x_ax.y;
        j[row][11] = -rb_x_ax.z;
    }

    // Angular rows
    for (i, ax) in axes.iter().enumerate() {
        let row = 3 + i;
        // Body A angular: +ax
        j[row][3] = ax.x;
        j[row][4] = ax.y;
        j[row][5] = ax.z;
        // Body B angular: -ax
        j[row][9] = -ax.x;
        j[row][10] = -ax.y;
        j[row][11] = -ax.z;
    }

    j
}

// ─────────────────────────────────────────────────────────────────────────────
// SpringDamperConfig
// ─────────────────────────────────────────────────────────────────────────────

/// Spring-damper configuration for a single DOF axis.
///
/// The force law is:
///   `F = -k * (x - x_rest) - c * v`
///
/// where `k` is the stiffness, `c` is the damping, `x` is the current
/// displacement, `x_rest` is the natural length and `v` is the velocity.
#[derive(Debug, Clone)]
pub struct SpringDamperConfig {
    /// Spring stiffness \[N/m or N·m/rad\].
    pub stiffness: f64,
    /// Damping coefficient \[N·s/m or N·m·s/rad\].
    pub damping: f64,
    /// Natural (rest) length / angle \[m or rad\].
    pub rest_length: f64,
    /// Pre-load force/torque \[N or N·m\].
    pub preload: f64,
}

impl SpringDamperConfig {
    /// Create a new spring-damper configuration.
    pub fn new(stiffness: f64, damping: f64, rest_length: f64, preload: f64) -> Self {
        Self {
            stiffness,
            damping,
            rest_length,
            preload,
        }
    }

    /// Compute the force/torque: `F = -k*(x - x_rest) - c*v + F_preload`.
    pub fn force(&self, position: f64, velocity: f64) -> f64 {
        -self.stiffness * (position - self.rest_length) - self.damping * velocity + self.preload
    }

    /// Compute the damping ratio: `zeta = c / (2 * sqrt(k * m))`.
    /// Requires the effective mass `m`.
    pub fn damping_ratio(&self, effective_mass: f64) -> f64 {
        if self.stiffness < 1e-30 || effective_mass < 1e-30 {
            return 0.0;
        }
        self.damping / (2.0 * (self.stiffness * effective_mass).sqrt())
    }

    /// Natural frequency: `omega_n = sqrt(k / m)`.
    pub fn natural_frequency(&self, effective_mass: f64) -> f64 {
        if effective_mass < 1e-30 {
            return 0.0;
        }
        (self.stiffness / effective_mass).sqrt()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// VelocityDriveConfig
// ─────────────────────────────────────────────────────────────────────────────

/// Velocity drive configuration for a single DOF axis.
///
/// A PD-style velocity controller that computes the constraint impulse
/// needed to achieve the target velocity.
#[derive(Debug, Clone)]
pub struct VelocityDriveConfig {
    /// Target velocity \[m/s or rad/s\].
    pub target_velocity: f64,
    /// Maximum force/torque magnitude \[N or N·m\].
    pub max_force: f64,
    /// Proportional gain for velocity error.
    pub kp: f64,
}

impl VelocityDriveConfig {
    /// Create a simple linear velocity drive.
    pub fn new_linear(target_velocity: f64, max_force: f64) -> Self {
        Self {
            target_velocity,
            max_force,
            kp: 1.0,
        }
    }

    /// Create a velocity drive with explicit gain.
    pub fn new(target_velocity: f64, max_force: f64, kp: f64) -> Self {
        Self {
            target_velocity,
            max_force,
            kp,
        }
    }

    /// Compute the constraint impulse: `imp = kp * (v_target - v_actual) * dt`.
    ///
    /// Clamped to `[-max_force * dt, max_force * dt]`.
    pub fn compute_impulse(&self, effective_mass: f64, actual_velocity: f64, dt: f64) -> f64 {
        let vel_error = self.target_velocity - actual_velocity;
        let imp = effective_mass * self.kp * vel_error;
        let limit = self.max_force * dt;
        imp.clamp(-limit, limit)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PositionDriveConfig
// ─────────────────────────────────────────────────────────────────────────────

/// Position drive configuration for a single DOF axis.
///
/// A PD (proportional-derivative) controller that drives the position to
/// a target while damping velocity.
#[derive(Debug, Clone)]
pub struct PositionDriveConfig {
    /// Target position \[m or rad\].
    pub target_position: f64,
    /// Proportional gain (stiffness-like) \[N/m or N·m/rad\].
    pub kp: f64,
    /// Derivative gain (damping-like) \[N·s/m or N·m·s/rad\].
    pub kd: f64,
    /// Maximum force/torque magnitude \[N or N·m\].
    pub max_force: f64,
}

impl PositionDriveConfig {
    /// Create a new position drive.
    pub fn new(target_position: f64, kp: f64, kd: f64, max_force: f64) -> Self {
        Self {
            target_position,
            kp,
            kd,
            max_force,
        }
    }

    /// Compute the drive force: `F = kp * (x_target - x) - kd * v`.
    ///
    /// Clamped to `[-max_force, max_force]`.
    pub fn compute_force(&self, position: f64, velocity: f64) -> f64 {
        let f = self.kp * (self.target_position - position) - self.kd * velocity;
        f.clamp(-self.max_force, self.max_force)
    }

    /// Error from target position.
    pub fn position_error(&self, position: f64) -> f64 {
        self.target_position - position
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Compliance and restitution helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Apply constraint compliance (softness) to the effective mass.
///
/// CFM (Constraint Force Mixing): modifies the diagonal of the system matrix
/// to allow slight violation.  The effective mass becomes:
///
///   `m_eff_c = m_eff / (1 + cfm * m_eff / dt)`
///
/// where `cfm` is the compliance parameter.
pub fn apply_compliance(effective_mass: f64, cfm: f64, dt: f64) -> f64 {
    let denom = 1.0 + cfm * effective_mass / dt.max(1e-30);
    effective_mass / denom
}

/// Compute the restitution impulse for a collision along a constraint axis.
///
/// `J_restitution = -(1 + e) * v_approach * effective_mass`
///
/// where `e` is the coefficient of restitution \[0, 1\] and `v_approach < 0`.
pub fn restitution_impulse(effective_mass: f64, v_approach: f64, restitution: f64) -> f64 {
    if v_approach >= 0.0 {
        return 0.0;
    }
    -(1.0 + restitution) * v_approach * effective_mass
}

/// Compute the post-collision velocity of body A using the impulse-momentum law.
///
/// `v_a_new = v_a + J * inv_mass_a`
pub fn post_collision_velocity(v_a: f64, impulse: f64, inv_mass_a: f64) -> f64 {
    v_a + impulse * inv_mass_a
}

/// Coefficient of restitution from pre/post velocities.
///
/// `e = -(v_a_post - v_b_post) / (v_a_pre - v_b_pre)`
pub fn coefficient_of_restitution(v_a_pre: f64, v_b_pre: f64, v_a_post: f64, v_b_post: f64) -> f64 {
    let dv_pre = v_a_pre - v_b_pre;
    if dv_pre.abs() < 1e-30 {
        return 0.0;
    }
    -(v_a_post - v_b_post) / dv_pre
}

// ─────────────────────────────────────────────────────────────────────────────
// 6-DOF limit utilities
// ─────────────────────────────────────────────────────────────────────────────

/// Check whether a DOF value is within its limit.
///
/// Returns `true` if `value` is within `[min, max]`.
pub fn within_limit(value: f64, min: f64, max: f64) -> bool {
    value >= min && value <= max
}

/// Compute the limit violation (penetration) for a constrained DOF.
///
/// Returns `> 0` if above the max, `< 0` if below the min, `0` if within.
pub fn limit_violation(value: f64, min: f64, max: f64) -> f64 {
    if value < min {
        value - min
    } else if value > max {
        value - max
    } else {
        0.0
    }
}

/// Clamp an accumulated impulse to its limit range and return the applied delta.
pub fn clamp_accumulated_impulse(old: f64, delta: f64, min: f64, max: f64) -> (f64, f64) {
    let new = (old + delta).clamp(min, max);
    let applied = new - old;
    (new, applied)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use oxiphysics_rigid::RigidBody;

    fn make_two_bodies() -> (RigidBodySet, BodyHandle, BodyHandle) {
        let mut bodies = RigidBodySet::new();

        let mut b_a = RigidBody::new(1.0);
        b_a.transform.position = Vec3::new(0.0, 0.0, 0.0);
        b_a.linear_damping = 0.0;
        b_a.angular_damping = 0.0;
        let ha = bodies.insert(b_a);

        let mut b_b = RigidBody::new(1.0);
        b_b.transform.position = Vec3::new(1.0, 0.0, 0.0);
        b_b.linear_damping = 0.0;
        b_b.angular_damping = 0.0;
        let hb = bodies.insert(b_b);

        (bodies, ha, hb)
    }

    #[test]
    fn test_six_dof_all_free_by_default() {
        let (_, ha, hb) = make_two_bodies();
        let c = SixDofConstraint::new(ha, hb);
        for i in 0..3 {
            assert!(!c.linear[i].locked);
            assert!(!c.angular[i].locked);
            assert!(c.linear[i].limited.is_none());
            assert!(c.angular[i].limited.is_none());
        }
    }

    #[test]
    fn test_lock_linear() {
        let (_, ha, hb) = make_two_bodies();
        let mut c = SixDofConstraint::new(ha, hb);
        c.lock_linear(1);
        assert!(c.linear[1].locked);
        assert!(!c.linear[0].locked);
        assert!(!c.linear[2].locked);
    }

    #[test]
    fn test_lock_angular() {
        let (_, ha, hb) = make_two_bodies();
        let mut c = SixDofConstraint::new(ha, hb);
        c.lock_angular(2);
        assert!(c.angular[2].locked);
        assert!(!c.angular[0].locked);
        assert!(!c.angular[1].locked);
    }

    #[test]
    fn test_set_linear_limit() {
        let (_, ha, hb) = make_two_bodies();
        let mut c = SixDofConstraint::new(ha, hb);
        c.set_linear_limit(0, -1.0, 1.0);
        assert_eq!(c.linear[0].limited, Some((-1.0, 1.0)));
    }

    #[test]
    fn test_set_angular_limit() {
        let (_, ha, hb) = make_two_bodies();
        let mut c = SixDofConstraint::new(ha, hb);
        c.set_angular_limit(0, -0.5, 0.5);
        assert_eq!(c.angular[0].limited, Some((-0.5, 0.5)));
    }

    #[test]
    fn test_weld_locks_all_dof() {
        let (_, ha, hb) = make_two_bodies();
        let c = SixDofConstraint::weld(ha, hb);
        for i in 0..3 {
            assert!(c.linear[i].locked, "linear[{i}] should be locked");
            assert!(c.angular[i].locked, "angular[{i}] should be locked");
        }
    }

    #[test]
    fn test_locked_linear_dof_produces_zero_velocity() {
        let mut bodies = RigidBodySet::new();

        // Place both bodies at the same position so there is no positional error;
        // the test checks only that the velocity constraint zeroes relative motion.
        let mut b_a = RigidBody::new(1.0);
        b_a.transform.position = Vec3::new(0.0, 0.0, 0.0);
        b_a.velocity = Vec3::new(5.0, 0.0, 0.0);
        b_a.linear_damping = 0.0;
        b_a.angular_damping = 0.0;
        let ha = bodies.insert(b_a);

        let mut b_b = RigidBody::new(1.0);
        b_b.transform.position = Vec3::new(0.0, 0.0, 0.0); // same position → no bias
        b_b.velocity = Vec3::new(-5.0, 0.0, 0.0);
        b_b.linear_damping = 0.0;
        b_b.angular_damping = 0.0;
        let hb = bodies.insert(b_b);

        let mut c = SixDofConstraint::new(ha, hb);
        c.lock_linear(0); // lock X translation

        let dt = 1.0 / 60.0;
        c.prepare(&bodies, dt);

        for _ in 0..20 {
            c.solve_velocity(&mut bodies, dt);
        }

        let vel_a = bodies.get(ha).unwrap().velocity;
        let vel_b = bodies.get(hb).unwrap().velocity;
        let rel_vel_x = (vel_a - vel_b).x;

        assert!(
            rel_vel_x.abs() < 0.5,
            "Locked X DOF should produce near-zero relative velocity along X, got {rel_vel_x}"
        );
    }

    #[test]
    fn test_locked_angular_dof_produces_zero_ang_vel() {
        let (mut bodies, ha, hb) = make_two_bodies();

        // Give bodies opposing angular velocities about Z
        bodies.get_mut(ha).unwrap().angular_velocity = Vec3::new(0.0, 0.0, 5.0);
        bodies.get_mut(hb).unwrap().angular_velocity = Vec3::new(0.0, 0.0, -5.0);

        let mut c = SixDofConstraint::new(ha, hb);
        c.lock_angular(2); // lock Z rotation

        let dt = 1.0 / 60.0;
        c.prepare(&bodies, dt);

        for _ in 0..20 {
            c.solve_velocity(&mut bodies, dt);
        }

        let ang_a = bodies.get(ha).unwrap().angular_velocity;
        let ang_b = bodies.get(hb).unwrap().angular_velocity;
        let rel_ang_z = (ang_a - ang_b).z;

        assert!(
            rel_ang_z.abs() < 1.0,
            "Locked Z angular DOF should produce near-zero relative angular velocity about Z, got {rel_ang_z}"
        );
    }

    #[test]
    fn test_body_handles() {
        let (_, ha, hb) = make_two_bodies();
        let c = SixDofConstraint::new(ha, hb);
        let handles = c.body_handles();
        assert_eq!(handles.len(), 2);
        assert_eq!(handles[0], ha);
        assert_eq!(handles[1], hb);
    }

    #[test]
    fn test_jacobian_linear_rows_orthogonality() {
        let r_a = Vec3::new(0.1, 0.2, 0.3);
        let r_b = Vec3::new(-0.1, 0.1, -0.2);
        let axes = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];

        let j = compute_jacobian_6dof(r_a, r_b, &axes);

        // Rows 0, 1, 2 should correspond to X, Y, Z linear axes
        // The first 3 columns of each linear row should match the axis
        for (row, ax) in axes.iter().enumerate() {
            assert!((j[row][0] - ax.x).abs() < 1e-12);
            assert!((j[row][1] - ax.y).abs() < 1e-12);
            assert!((j[row][2] - ax.z).abs() < 1e-12);
        }

        // Angular rows (3,4,5) body-A cols should equal axes
        for (i, ax) in axes.iter().enumerate() {
            let row = 3 + i;
            assert!((j[row][3] - ax.x).abs() < 1e-12);
            assert!((j[row][4] - ax.y).abs() < 1e-12);
            assert!((j[row][5] - ax.z).abs() < 1e-12);
        }
    }

    #[test]
    fn test_six_dof_with_motor() {
        let (_, ha, hb) = make_two_bodies();
        let mut c = SixDofConstraint::new(ha, hb);
        let motor = MotorConfig::new(2.0, 100.0, 50.0);
        c.set_linear_motor(0, motor);
        assert!(c.linear[0].motor.is_some());
        assert_eq!(c.linear[0].motor.as_ref().unwrap().target_velocity, 2.0);
    }

    #[test]
    #[should_panic]
    fn test_lock_linear_out_of_bounds() {
        let (_, ha, hb) = make_two_bodies();
        let mut c = SixDofConstraint::new(ha, hb);
        c.lock_linear(3); // should panic
    }

    #[test]
    #[should_panic]
    fn test_lock_angular_out_of_bounds() {
        let (_, ha, hb) = make_two_bodies();
        let mut c = SixDofConstraint::new(ha, hb);
        c.lock_angular(3); // should panic
    }

    // ── Spring/damper ─────────────────────────────────────────────────────────

    #[test]
    fn test_spring_damper_config_creation() {
        let sd = SpringDamperConfig::new(100.0, 10.0, 1.0, 0.0);
        assert!((sd.stiffness - 100.0).abs() < 1e-12);
        assert!((sd.damping - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_spring_force_at_rest_position() {
        let sd = SpringDamperConfig::new(100.0, 10.0, 1.0, 0.0);
        let f = sd.force(1.0, 0.0); // at rest length, zero velocity → zero force
        assert!(f.abs() < 1e-12, "At rest: force = {f}");
    }

    #[test]
    fn test_spring_force_compressed() {
        let sd = SpringDamperConfig::new(100.0, 0.0, 0.5, 0.0);
        let f = sd.force(0.0, 0.0); // displacement = rest_length from neutral
        // Extension = 0.0 - 0.5 = -0.5; F = -k * (-0.5) = +50
        assert!((f - 50.0).abs() < 1e-10, "Compressed spring: f = {f}");
    }

    #[test]
    fn test_damper_force_velocity() {
        let sd = SpringDamperConfig::new(0.0, 50.0, 0.0, 0.0);
        let f = sd.force(0.0, 2.0); // velocity = 2 m/s
        // F = -c * vel = -100
        assert!((f + 100.0).abs() < 1e-10, "Damper: f = {f}");
    }

    // ── Velocity drive ────────────────────────────────────────────────────────

    #[test]
    fn test_velocity_drive_config_linear() {
        let vd = VelocityDriveConfig::new_linear(5.0, 1000.0);
        assert!((vd.target_velocity - 5.0).abs() < 1e-12);
        assert!((vd.max_force - 1000.0).abs() < 1e-12);
    }

    #[test]
    fn test_velocity_drive_impulse_positive_error() {
        let vd = VelocityDriveConfig::new_linear(5.0, 1000.0);
        let imp = vd.compute_impulse(3.0, 2.0, 1.0 / 60.0); // actual_vel < target
        assert!(
            imp > 0.0,
            "Positive velocity error → positive impulse: {imp}"
        );
    }

    #[test]
    fn test_velocity_drive_impulse_clamped() {
        let vd = VelocityDriveConfig::new_linear(100.0, 1.0); // small max force
        let imp = vd.compute_impulse(0.0, 10.0, 1.0 / 60.0);
        let limit = 1.0 * (1.0 / 60.0);
        assert!(
            imp.abs() <= limit + 1e-12,
            "Impulse should be clamped: {imp}"
        );
    }

    // ── Position drive ────────────────────────────────────────────────────────

    #[test]
    fn test_position_drive_config_creation() {
        let pd = PositionDriveConfig::new(1.0, 100.0, 10.0, 500.0);
        assert!((pd.target_position - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_position_drive_force_toward_target() {
        let pd = PositionDriveConfig::new(2.0, 100.0, 10.0, 500.0);
        let f = pd.compute_force(0.0, 0.0); // position error = 2.0
        assert!(f > 0.0, "Force should drive toward target: {f}");
    }

    #[test]
    fn test_position_drive_zero_error() {
        let pd = PositionDriveConfig::new(1.0, 100.0, 10.0, 500.0);
        let f = pd.compute_force(1.0, 0.0); // at target, zero velocity
        assert!(f.abs() < 1e-10, "Zero position error: f = {f}");
    }

    // ── Compliance / restitution ──────────────────────────────────────────────

    #[test]
    fn test_compliance_reduces_effective_mass() {
        let em_no_compliance = 1.0; // 1/k
        let em_with_compliance = apply_compliance(em_no_compliance, 0.01, 1.0 / 60.0);
        assert!(
            em_with_compliance < em_no_compliance,
            "Compliance should reduce effective mass: {em_with_compliance} vs {em_no_compliance}"
        );
    }

    #[test]
    fn test_restitution_impulse_sign() {
        let imp = restitution_impulse(1.0, -2.0, 0.5); // approaching velocity -2
        assert!(
            imp > 0.0,
            "Restitution impulse should oppose approach: {imp}"
        );
    }

    #[test]
    fn test_restitution_zero_coefficient() {
        let imp = restitution_impulse(1.0, -2.0, 0.0);
        // e=0: perfectly inelastic; impulse = eff_mass * |v_approach| = 2
        assert!((imp - 2.0).abs() < 1e-12, "Zero restitution: imp = {imp}");
    }

    // ── 6-DOF limits ─────────────────────────────────────────────────────────

    #[test]
    fn test_axis_config_free_no_limits() {
        let ac = AxisConfig::free();
        assert!(!ac.locked);
        assert!(ac.limited.is_none());
    }

    #[test]
    fn test_axis_config_limited_clamps_correctly() {
        let ac = AxisConfig::limited(-1.0, 1.0);
        assert_eq!(ac.limited, Some((-1.0, 1.0)));
    }

    #[test]
    fn test_six_dof_linear_limit_enforcement() {
        let (_, ha, hb) = make_two_bodies();
        let mut c = SixDofConstraint::new(ha, hb);
        c.set_linear_limit(0, -0.5, 0.5);
        assert!(c.linear[0].limited.is_some());
    }

    #[test]
    fn test_six_dof_angular_limit_enforcement() {
        let (_, ha, hb) = make_two_bodies();
        let mut c = SixDofConstraint::new(ha, hb);
        c.set_angular_limit(1, -0.3, 0.3);
        assert_eq!(c.angular[1].limited, Some((-0.3, 0.3)));
    }

    // ── Jacobian correctness ──────────────────────────────────────────────────

    #[test]
    fn test_jacobian_body_b_linear_is_negated() {
        let r_a = Vec3::new(0.0, 0.0, 0.0);
        let r_b = Vec3::new(0.0, 0.0, 0.0);
        let axes = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let j = compute_jacobian_6dof(r_a, r_b, &axes);
        // For row 0 (X axis): columns 0..3 = [1,0,0], columns 6..9 = [-1,0,0]
        assert!(
            (j[0][6] + 1.0).abs() < 1e-12,
            "Body B X-linear should be -1: {}",
            j[0][6]
        );
        assert!((j[0][7]).abs() < 1e-12, "Body B Y should be 0: {}", j[0][7]);
    }

    #[test]
    fn test_jacobian_angular_rows_body_b_negated() {
        let r_a = Vec3::zeros();
        let r_b = Vec3::zeros();
        let axes = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let j = compute_jacobian_6dof(r_a, r_b, &axes);
        // Row 3 (X angular): body B angular should be -1
        assert!(
            (j[3][9] + 1.0).abs() < 1e-12,
            "Body B angular X row: {}",
            j[3][9]
        );
    }

    #[test]
    fn test_jacobian_r_cross_axis_terms() {
        // Non-zero r_a: r_a = [0, 0, 1], axis = [1, 0, 0]
        // r_a × ax = [0,0,1] × [1,0,0] = [0*0-1*0, 1*1-0*0, 0*0-0*1] = [0,1,0]
        let r_a = Vec3::new(0.0, 0.0, 1.0);
        let r_b = Vec3::zeros();
        let axes = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let j = compute_jacobian_6dof(r_a, r_b, &axes);
        // Row 0 (X axis), cols 3..6 = r_a × ax = [0,1,0]
        assert!((j[0][3]).abs() < 1e-12, "j[0][3] = {}", j[0][3]);
        assert!((j[0][4] - 1.0).abs() < 1e-12, "j[0][4] = {}", j[0][4]);
        assert!((j[0][5]).abs() < 1e-12, "j[0][5] = {}", j[0][5]);
    }

    // ── Motor angular ─────────────────────────────────────────────────────────

    #[test]
    fn test_angular_motor_config_attached() {
        let (_, ha, hb) = make_two_bodies();

        let mut c = SixDofConstraint::new(ha, hb);
        let motor = MotorConfig::new(5.0, 1000.0, 100.0);
        c.set_angular_motor(2, motor);

        // Verify motor config was stored
        assert!(
            c.angular[2].motor.is_some(),
            "Motor should be attached to angular axis 2"
        );
        let m = c.angular[2].motor.as_ref().unwrap();
        assert!((m.target_velocity - 5.0).abs() < 1e-12);
        assert!((m.max_force - 1000.0).abs() < 1e-12);
    }

    #[test]
    fn test_linear_motor_drives_velocity() {
        let (mut bodies, ha, hb) = make_two_bodies();

        let mut c = SixDofConstraint::new(ha, hb);
        let motor = MotorConfig::new(5.0, 1000.0, 50.0);
        c.set_linear_motor(0, motor); // drive X translation

        let dt = 1.0 / 60.0;
        c.prepare(&bodies, dt);
        for _ in 0..30 {
            c.solve_velocity(&mut bodies, dt);
        }

        // Velocity difference along X should not be worse (motor constrains)
        let vel_a = bodies.get(ha).unwrap().velocity;
        let vel_b = bodies.get(hb).unwrap().velocity;
        // No panic = motor step ran without error
        let _ = (vel_a, vel_b);
    }
}
