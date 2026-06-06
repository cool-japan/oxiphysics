// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Constraint system API for Python interop.
//!
//! Provides Python-friendly types for rigid body constraint solving,
//! joints, contact constraints, motors, island management, CCD,
//! PBD/XPBD solvers, control systems, friction models, and warm starting.

use pyo3::prelude::*;
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

/// Perform one PGS (Projected Gauss-Seidel) iteration over a set of constraints.
///
/// `lambda` — accumulated impulse per constraint (modified in place).
/// `rhs`    — constraint violation / desired velocity change per constraint.
/// `diag`   — diagonal of the effective-mass matrix per constraint.
/// `lo`     — lower bound per constraint (clamping).
/// `hi`     — upper bound per constraint (clamping).
///
/// Returns the total residual after the iteration.
pub fn solve_pgs_iteration(
    lambda: &mut [f64],
    rhs: &[f64],
    diag: &[f64],
    lo: &[f64],
    hi: &[f64],
) -> f64 {
    let n = lambda.len();
    let mut residual = 0.0;
    for i in 0..n {
        if diag[i].abs() < 1e-15 {
            continue;
        }
        let delta = (rhs[i] - diag[i] * lambda[i]) / diag[i];
        let new_lambda = (lambda[i] + delta).clamp(lo[i], hi[i]);
        let actual_delta = new_lambda - lambda[i];
        lambda[i] = new_lambda;
        residual += actual_delta * actual_delta;
    }
    residual.sqrt()
}

/// Python-exposed wrapper for PGS iteration.
///
/// Takes lambda as a Vec (returns modified vec) and all other slices as Vec.
/// Returns `(updated_lambda, residual)`.
#[pyfunction]
pub fn py_solve_pgs_iteration(
    mut lambda: Vec<f64>,
    rhs: Vec<f64>,
    diag: Vec<f64>,
    lo: Vec<f64>,
    hi: Vec<f64>,
) -> (Vec<f64>, f64) {
    let residual = solve_pgs_iteration(&mut lambda, &rhs, &diag, &lo, &hi);
    (lambda, residual)
}

/// Compute the Jacobian row for a distance constraint between two bodies.
///
/// `r_a` — vector from body A's center to the contact point.
/// `r_b` — vector from body B's center to the contact point.
/// `n`   — constraint normal direction (unit vector).
///
/// Returns \[J_lin_a (3), J_ang_a (3), J_lin_b (3), J_ang_b (3)\] = 12 values.
#[pyfunction]
pub fn compute_jacobian(r_a: [f64; 3], r_b: [f64; 3], n: [f64; 3]) -> [f64; 12] {
    // J_lin_a = n, J_ang_a = r_a × n, J_lin_b = -n, J_ang_b = -r_b × n
    let ang_a = cross(r_a, n);
    let ang_b = cross(r_b, n);
    [
        n[0], n[1], n[2], ang_a[0], ang_a[1], ang_a[2], -n[0], -n[1], -n[2], -ang_b[0], -ang_b[1],
        -ang_b[2],
    ]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Compute the effective mass (inverse denominator) for a constraint.
///
/// `inv_mass_a`, `inv_mass_b` — inverse masses.
/// `inv_inertia_a`, `inv_inertia_b` — diagonal inverse inertia tensors (3 values each).
/// `j` — Jacobian row (12 values from `compute_jacobian`).
///
/// Returns the effective mass = 1 / (J M⁻¹ Jᵀ).
#[pyfunction]
pub fn compute_effective_mass(
    inv_mass_a: f64,
    inv_mass_b: f64,
    inv_inertia_a: [f64; 3],
    inv_inertia_b: [f64; 3],
    j: [f64; 12],
) -> f64 {
    // linear part
    let lin_a = j[0] * j[0] * inv_mass_a + j[1] * j[1] * inv_mass_a + j[2] * j[2] * inv_mass_a;
    let ang_a = j[3] * j[3] * inv_inertia_a[0]
        + j[4] * j[4] * inv_inertia_a[1]
        + j[5] * j[5] * inv_inertia_a[2];
    let lin_b = j[6] * j[6] * inv_mass_b + j[7] * j[7] * inv_mass_b + j[8] * j[8] * inv_mass_b;
    let ang_b = j[9] * j[9] * inv_inertia_b[0]
        + j[10] * j[10] * inv_inertia_b[1]
        + j[11] * j[11] * inv_inertia_b[2];
    let denom = lin_a + ang_a + lin_b + ang_b;
    if denom.abs() < 1e-15 {
        0.0
    } else {
        1.0 / denom
    }
}

/// Clamp an impulse to the friction cone.
///
/// `lambda_n` — normal impulse (must be ≥ 0).
/// `lambda_t` — current tangential impulse vector \[lt_x, lt_y\].
/// `mu`       — friction coefficient.
///
/// Returns the clamped tangential impulse.
#[pyfunction]
pub fn clamp_impulse(lambda_n: f64, lambda_t: [f64; 2], mu: f64) -> [f64; 2] {
    let max_t = mu * lambda_n.max(0.0);
    let mag = (lambda_t[0] * lambda_t[0] + lambda_t[1] * lambda_t[1]).sqrt();
    if mag > max_t && mag > 1e-15 {
        [lambda_t[0] / mag * max_t, lambda_t[1] / mag * max_t]
    } else {
        lambda_t
    }
}

// ---------------------------------------------------------------------------
// Solver type
// ---------------------------------------------------------------------------

/// Constraint solver algorithm selection.
#[pyclass(eq, eq_int, from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum SolverType {
    /// Projected Gauss-Seidel.
    Pgs,
    /// Temporal Gauss-Seidel (velocity + position level).
    Tgs,
    /// Sequential impulse (SI).
    Si,
}

/// Main constraint solver configuration.
#[pyclass(get_all, set_all, from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyConstraintSolver {
    /// Solver algorithm.
    pub solver_type: SolverType,
    /// Number of velocity-level iterations.
    pub velocity_iterations: u32,
    /// Number of position-level iterations (for TGS/non-penetration correction).
    pub position_iterations: u32,
    /// Whether to use warm starting from the previous frame.
    pub warm_start: bool,
    /// Successive over-relaxation factor (1.0 = standard PGS).
    pub sor_factor: f64,
    /// Convergence tolerance.
    pub tolerance: f64,
    /// Maximum allowed penetration before applying a correction impulse.
    pub slop: f64,
}

#[pymethods]
impl PyConstraintSolver {
    /// Create a default PGS solver.
    #[staticmethod]
    pub fn default_pgs() -> Self {
        Self {
            solver_type: SolverType::Pgs,
            velocity_iterations: 10,
            position_iterations: 5,
            warm_start: true,
            sor_factor: 1.0,
            tolerance: 1e-4,
            slop: 0.005,
        }
    }

    /// Create a TGS solver (fewer iterations needed).
    #[staticmethod]
    pub fn default_tgs() -> Self {
        Self {
            solver_type: SolverType::Tgs,
            velocity_iterations: 4,
            position_iterations: 2,
            warm_start: true,
            sor_factor: 1.3,
            tolerance: 1e-4,
            slop: 0.005,
        }
    }

    /// Solve a simple set of constraints given accumulated lambdas, RHS and diagonal effective masses.
    ///
    /// Returns the final accumulated impulse array.
    pub fn solve(&self, rhs: Vec<f64>, diag: Vec<f64>, lo: Vec<f64>, hi: Vec<f64>) -> Vec<f64> {
        let n = rhs.len();
        let mut lambda = vec![0.0f64; n];
        for _ in 0..self.velocity_iterations {
            let res = solve_pgs_iteration(&mut lambda, &rhs, &diag, &lo, &hi);
            if res < self.tolerance {
                break;
            }
        }
        lambda
    }
}

// ---------------------------------------------------------------------------
// PyJoint
// ---------------------------------------------------------------------------

/// Joint axis and limit specification.
#[pyclass(get_all, set_all, from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AxisLimits {
    /// Axis direction (unit vector).
    pub axis: [f64; 3],
    /// Lower limit (radians or meters).
    pub lower: f64,
    /// Upper limit.
    pub upper: f64,
    /// Whether limits are enabled.
    pub enabled: bool,
}

#[pymethods]
impl AxisLimits {
    /// Unlimited joint along an axis.
    #[staticmethod]
    pub fn unlimited(axis: [f64; 3]) -> Self {
        Self {
            axis,
            lower: -f64::INFINITY,
            upper: f64::INFINITY,
            enabled: false,
        }
    }

    /// Limited joint.
    #[staticmethod]
    pub fn limited(axis: [f64; 3], lower: f64, upper: f64) -> Self {
        Self {
            axis,
            lower,
            upper,
            enabled: true,
        }
    }
}

/// The kind of joint.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum JointKind {
    /// Fully fixed — all 6 DOF locked.
    Fixed,
    /// Revolute (hinge) joint — 1 rotational DOF.
    Revolute(AxisLimits),
    /// Prismatic (sliding) joint — 1 translational DOF.
    Prismatic(AxisLimits),
    /// Ball-and-socket joint with swing and twist limits.
    Ball { swing_limit: f64, twist_limit: f64 },
    /// Spring joint with stiffness and damping.
    Spring {
        stiffness: f64,
        damping: f64,
        rest_length: f64,
    },
}

/// A joint connecting two rigid bodies.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyJoint {
    /// Unique joint identifier.
    pub id: u32,
    /// Body A handle.
    pub body_a: u32,
    /// Body B handle.
    pub body_b: u32,
    /// Anchor point on body A in local space.
    pub anchor_a: [f64; 3],
    /// Anchor point on body B in local space.
    pub anchor_b: [f64; 3],
    /// Joint type.
    pub kind: JointKind,
    /// Whether this joint is enabled.
    pub enabled: bool,
    /// Joint breaking force (inf = unbreakable).
    pub break_force: f64,
}

#[pymethods]
impl PyJoint {
    /// Create a fixed joint between two bodies.
    #[staticmethod]
    pub fn fixed(
        id: u32,
        body_a: u32,
        body_b: u32,
        anchor_a: [f64; 3],
        anchor_b: [f64; 3],
    ) -> Self {
        Self {
            id,
            body_a,
            body_b,
            anchor_a,
            anchor_b,
            kind: JointKind::Fixed,
            enabled: true,
            break_force: f64::INFINITY,
        }
    }

    /// Create a revolute joint.
    #[staticmethod]
    pub fn revolute(id: u32, body_a: u32, body_b: u32, anchor: [f64; 3], axis: [f64; 3]) -> Self {
        Self {
            id,
            body_a,
            body_b,
            anchor_a: anchor,
            anchor_b: anchor,
            kind: JointKind::Revolute(AxisLimits::unlimited(axis)),
            enabled: true,
            break_force: f64::INFINITY,
        }
    }

    /// Create a spring joint.
    #[staticmethod]
    pub fn spring(
        id: u32,
        body_a: u32,
        body_b: u32,
        anchor_a: [f64; 3],
        anchor_b: [f64; 3],
        stiffness: f64,
        damping: f64,
        rest_length: f64,
    ) -> Self {
        Self {
            id,
            body_a,
            body_b,
            anchor_a,
            anchor_b,
            kind: JointKind::Spring {
                stiffness,
                damping,
                rest_length,
            },
            enabled: true,
            break_force: f64::INFINITY,
        }
    }

    /// Check if the joint is breakable.
    pub fn is_breakable(&self) -> bool {
        self.break_force.is_finite()
    }

    /// Get id.
    #[getter]
    pub fn get_id(&self) -> u32 {
        self.id
    }

    /// Get body_a.
    #[getter]
    pub fn get_body_a(&self) -> u32 {
        self.body_a
    }

    /// Get body_b.
    #[getter]
    pub fn get_body_b(&self) -> u32 {
        self.body_b
    }

    /// Get anchor_a.
    #[getter]
    pub fn get_anchor_a(&self) -> [f64; 3] {
        self.anchor_a
    }

    /// Get anchor_b.
    #[getter]
    pub fn get_anchor_b(&self) -> [f64; 3] {
        self.anchor_b
    }

    /// Get enabled.
    #[getter]
    pub fn get_enabled(&self) -> bool {
        self.enabled
    }

    /// Set enabled.
    #[setter]
    pub fn set_enabled(&mut self, v: bool) {
        self.enabled = v;
    }

    /// Get break_force.
    #[getter]
    pub fn get_break_force(&self) -> f64 {
        self.break_force
    }

    /// Set break_force.
    #[setter]
    pub fn set_break_force(&mut self, v: f64) {
        self.break_force = v;
    }
}

// ---------------------------------------------------------------------------
// PyContactConstraint
// ---------------------------------------------------------------------------

/// A single contact point constraint between two bodies.
#[pyclass(get_all, set_all, from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyContactConstraint {
    /// Body A handle.
    pub body_a: u32,
    /// Body B handle.
    pub body_b: u32,
    /// Contact normal (from B to A, unit vector).
    pub normal: [f64; 3],
    /// Penetration depth (positive = overlap).
    pub penetration: f64,
    /// Contact position in world space.
    pub position: [f64; 3],
    /// Accumulated normal impulse (warm start).
    pub lambda_n: f64,
    /// Accumulated tangent impulse in primary direction.
    pub lambda_t1: f64,
    /// Accumulated tangent impulse in secondary direction.
    pub lambda_t2: f64,
    /// Coefficient of restitution.
    pub restitution: f64,
    /// Coulomb friction coefficient.
    pub friction: f64,
}

#[pymethods]
impl PyContactConstraint {
    /// Create a new contact constraint.
    #[new]
    pub fn new(
        body_a: u32,
        body_b: u32,
        normal: [f64; 3],
        penetration: f64,
        position: [f64; 3],
        restitution: f64,
        friction: f64,
    ) -> Self {
        Self {
            body_a,
            body_b,
            normal,
            penetration,
            position,
            lambda_n: 0.0,
            lambda_t1: 0.0,
            lambda_t2: 0.0,
            restitution,
            friction,
        }
    }

    /// Clamp the tangential impulses to the friction cone.
    pub fn clamp_friction(&mut self) {
        let clamped = clamp_impulse(
            self.lambda_n,
            [self.lambda_t1, self.lambda_t2],
            self.friction,
        );
        self.lambda_t1 = clamped[0];
        self.lambda_t2 = clamped[1];
    }

    /// Compute the primary and secondary tangent directions orthogonal to the normal.
    pub fn tangent_basis(&self) -> ([f64; 3], [f64; 3]) {
        let n = self.normal;
        let t1 = if n[0].abs() < 0.9 {
            let raw = [0.0 - n[1] * n[1], n[0] * n[1] - 0.0, 0.0];
            // Just use a simple perpendicular
            let t = [1.0 - n[0] * n[0], -n[0] * n[1], -n[0] * n[2]];
            let len = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
            if len > 1e-12 {
                [t[0] / len, t[1] / len, t[2] / len]
            } else {
                raw
            }
        } else {
            let t = [-n[1] * n[0], 1.0 - n[1] * n[1], -n[1] * n[2]];
            let len = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
            if len > 1e-12 {
                [t[0] / len, t[1] / len, t[2] / len]
            } else {
                [0.0, 1.0, 0.0]
            }
        };
        let t2 = cross(n, t1);
        (t1, t2)
    }
}

// ---------------------------------------------------------------------------
// PyMotorConstraint
// ---------------------------------------------------------------------------

/// PID controller for motor target tracking.
#[pyclass(get_all, set_all, from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PidGains {
    /// Proportional gain.
    pub kp: f64,
    /// Integral gain.
    pub ki: f64,
    /// Derivative gain.
    pub kd: f64,
    /// Accumulated integral error.
    pub integral: f64,
    /// Previous error (for derivative).
    pub prev_error: f64,
}

#[pymethods]
impl PidGains {
    /// Create new PID gains.
    #[new]
    pub fn new(kp: f64, ki: f64, kd: f64) -> Self {
        Self {
            kp,
            ki,
            kd,
            integral: 0.0,
            prev_error: 0.0,
        }
    }

    /// Compute PID output given the current error and time step.
    pub fn compute(&mut self, error: f64, dt: f64) -> f64 {
        self.integral += error * dt;
        let derivative = if dt > 1e-15 {
            (error - self.prev_error) / dt
        } else {
            0.0
        };
        self.prev_error = error;
        self.kp * error + self.ki * self.integral + self.kd * derivative
    }

    /// Reset the integral and derivative state.
    pub fn reset(&mut self) {
        self.integral = 0.0;
        self.prev_error = 0.0;
    }
}

/// Motor thermal model to cap torque under sustained load.
#[pyclass(get_all, set_all, from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MotorThermal {
    /// Current motor temperature (°C).
    pub temperature: f64,
    /// Ambient temperature (°C).
    pub ambient: f64,
    /// Thermal resistance (°C/W).
    pub resistance: f64,
    /// Thermal capacitance (J/°C).
    pub capacitance: f64,
    /// Motor resistance for heat generation (Ω).
    pub motor_resistance: f64,
    /// Temperature derating threshold (°C).
    pub derate_threshold: f64,
    /// Maximum allowed temperature (°C).
    pub max_temperature: f64,
}

#[pymethods]
impl MotorThermal {
    /// Create a default thermal model.
    #[new]
    pub fn new() -> Self {
        Self {
            temperature: 25.0,
            ambient: 25.0,
            resistance: 0.5,
            capacitance: 100.0,
            motor_resistance: 0.1,
            derate_threshold: 80.0,
            max_temperature: 120.0,
        }
    }

    /// Update temperature given current (A) and time step (s).
    pub fn update(&mut self, current: f64, dt: f64) {
        let heat_in = current * current * self.motor_resistance;
        let heat_out = (self.temperature - self.ambient) / self.resistance;
        self.temperature += (heat_in - heat_out) / self.capacitance * dt;
    }

    /// Compute the torque scale factor (1 = full, < 1 = derated).
    pub fn torque_scale(&self) -> f64 {
        if self.temperature <= self.derate_threshold {
            1.0
        } else {
            let range = self.max_temperature - self.derate_threshold;
            if range < 1e-9 {
                0.0
            } else {
                (1.0 - (self.temperature - self.derate_threshold) / range).max(0.0)
            }
        }
    }
}

impl Default for MotorThermal {
    fn default() -> Self {
        Self::new()
    }
}

/// A motor/actuator constraint driving a joint.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyMotorConstraint {
    /// Target angular velocity (rad/s). Used in velocity mode.
    pub target_velocity: f64,
    /// Target position (rad or m). Used in position mode.
    pub target_position: f64,
    /// Maximum torque/force (N·m or N).
    pub max_torque: f64,
    /// PID controller.
    pub pid: PidGains,
    /// Whether in position mode (true) or velocity mode (false).
    pub position_mode: bool,
    /// Current applied torque.
    pub applied_torque: f64,
    /// Thermal model.
    pub thermal: MotorThermal,
}

#[pymethods]
impl PyMotorConstraint {
    /// Create a velocity-mode motor.
    #[staticmethod]
    pub fn velocity_mode(target_velocity: f64, max_torque: f64, kp: f64) -> Self {
        Self {
            target_velocity,
            target_position: 0.0,
            max_torque,
            pid: PidGains::new(kp, 0.0, 0.0),
            position_mode: false,
            applied_torque: 0.0,
            thermal: MotorThermal::new(),
        }
    }

    /// Create a position-mode motor with PID gains.
    #[staticmethod]
    pub fn position_mode(target_position: f64, max_torque: f64, pid: PidGains) -> Self {
        Self {
            target_velocity: 0.0,
            target_position,
            max_torque,
            pid,
            position_mode: true,
            applied_torque: 0.0,
            thermal: MotorThermal::new(),
        }
    }

    /// Step the motor given the current measured value and time step.
    pub fn step(&mut self, current_value: f64, dt: f64) -> f64 {
        let error = if self.position_mode {
            self.target_position - current_value
        } else {
            self.target_velocity - current_value
        };
        let raw_torque = self.pid.compute(error, dt);
        let scale = self.thermal.torque_scale();
        let torque = raw_torque.clamp(-self.max_torque, self.max_torque) * scale;
        let current = (torque / self.max_torque.max(1e-9)).abs() * 10.0;
        self.thermal.update(current, dt);
        self.applied_torque = torque;
        torque
    }

    /// Get target_velocity.
    #[getter]
    pub fn get_target_velocity(&self) -> f64 {
        self.target_velocity
    }

    /// Set target_velocity.
    #[setter]
    pub fn set_target_velocity(&mut self, v: f64) {
        self.target_velocity = v;
    }

    /// Get target_position.
    #[getter]
    pub fn get_target_position(&self) -> f64 {
        self.target_position
    }

    /// Set target_position.
    #[setter]
    pub fn set_target_position(&mut self, v: f64) {
        self.target_position = v;
    }

    /// Get max_torque.
    #[getter]
    pub fn get_max_torque(&self) -> f64 {
        self.max_torque
    }

    /// Get pid (cloned).
    #[getter]
    pub fn get_pid(&self) -> PidGains {
        self.pid.clone()
    }

    /// Set pid.
    #[setter]
    pub fn set_pid(&mut self, v: PidGains) {
        self.pid = v;
    }

    /// Get position_mode.
    #[getter]
    pub fn get_position_mode(&self) -> bool {
        self.position_mode
    }

    /// Get applied_torque.
    #[getter]
    pub fn get_applied_torque(&self) -> f64 {
        self.applied_torque
    }

    /// Get thermal model (cloned).
    #[getter]
    pub fn get_thermal(&self) -> MotorThermal {
        self.thermal.clone()
    }

    /// Set thermal model.
    #[setter]
    pub fn set_thermal(&mut self, v: MotorThermal) {
        self.thermal = v;
    }
}

// ---------------------------------------------------------------------------
// PyIsland
// ---------------------------------------------------------------------------

/// Sleeping state of an island.
#[pyclass(eq, eq_int, from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum IslandSleepState {
    /// Island is awake and being simulated.
    Awake,
    /// Island is dormant (below velocity threshold).
    Sleeping,
    /// Island is in the process of going to sleep.
    Drowsy,
}

/// An island groups bodies and constraints that can interact.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyIsland {
    /// Unique island id.
    pub id: u32,
    /// Handles of rigid bodies in this island.
    pub bodies: Vec<u32>,
    /// Handles of constraints in this island.
    pub constraints: Vec<u32>,
    /// Sleep state.
    pub sleep_state: IslandSleepState,
    /// Frames below the velocity threshold (used for sleep triggering).
    pub quiet_frames: u32,
    /// Velocity threshold for sleeping (m/s and rad/s combined).
    pub sleep_threshold: f64,
    /// Number of frames needed to trigger sleep.
    pub sleep_frames: u32,
}

#[pymethods]
impl PyIsland {
    /// Create a new island.
    #[new]
    pub fn new(id: u32) -> Self {
        Self {
            id,
            bodies: Vec::new(),
            constraints: Vec::new(),
            sleep_state: IslandSleepState::Awake,
            quiet_frames: 0,
            sleep_threshold: 0.05,
            sleep_frames: 60,
        }
    }

    /// Add a body to the island.
    pub fn add_body(&mut self, body: u32) {
        if !self.bodies.contains(&body) {
            self.bodies.push(body);
        }
    }

    /// Add a constraint to the island.
    pub fn add_constraint(&mut self, constraint: u32) {
        if !self.constraints.contains(&constraint) {
            self.constraints.push(constraint);
        }
    }

    /// Update sleep state given the maximum body speed this frame.
    pub fn update_sleep(&mut self, max_speed: f64) {
        match self.sleep_state {
            IslandSleepState::Sleeping => {
                if max_speed > self.sleep_threshold * 2.0 {
                    self.sleep_state = IslandSleepState::Awake;
                    self.quiet_frames = 0;
                }
            }
            _ => {
                if max_speed < self.sleep_threshold {
                    self.quiet_frames += 1;
                    if self.quiet_frames >= self.sleep_frames {
                        self.sleep_state = IslandSleepState::Sleeping;
                    } else {
                        self.sleep_state = IslandSleepState::Drowsy;
                    }
                } else {
                    self.quiet_frames = 0;
                    self.sleep_state = IslandSleepState::Awake;
                }
            }
        }
    }

    /// Merge another island into this one.
    pub fn merge(&mut self, other: &PyIsland) {
        for &b in &other.bodies {
            self.add_body(b);
        }
        for &c in &other.constraints {
            self.add_constraint(c);
        }
        if other.sleep_state == IslandSleepState::Awake {
            self.sleep_state = IslandSleepState::Awake;
        }
    }

    /// Whether this island is currently active.
    pub fn is_active(&self) -> bool {
        self.sleep_state != IslandSleepState::Sleeping
    }

    /// Get id.
    #[getter]
    pub fn get_id(&self) -> u32 {
        self.id
    }

    /// Get bodies.
    #[getter]
    pub fn get_bodies(&self) -> Vec<u32> {
        self.bodies.clone()
    }

    /// Set bodies.
    #[setter]
    pub fn set_bodies(&mut self, v: Vec<u32>) {
        self.bodies = v;
    }

    /// Get constraints.
    #[getter]
    pub fn get_constraints(&self) -> Vec<u32> {
        self.constraints.clone()
    }

    /// Set constraints.
    #[setter]
    pub fn set_constraints(&mut self, v: Vec<u32>) {
        self.constraints = v;
    }

    /// Get sleep_state.
    #[getter]
    pub fn get_sleep_state(&self) -> IslandSleepState {
        self.sleep_state.clone()
    }

    /// Set sleep_state.
    #[setter]
    pub fn set_sleep_state(&mut self, v: IslandSleepState) {
        self.sleep_state = v;
    }

    /// Get quiet_frames.
    #[getter]
    pub fn get_quiet_frames(&self) -> u32 {
        self.quiet_frames
    }

    /// Get sleep_threshold.
    #[getter]
    pub fn get_sleep_threshold(&self) -> f64 {
        self.sleep_threshold
    }

    /// Set sleep_threshold.
    #[setter]
    pub fn set_sleep_threshold(&mut self, v: f64) {
        self.sleep_threshold = v;
    }

    /// Get sleep_frames.
    #[getter]
    pub fn get_sleep_frames(&self) -> u32 {
        self.sleep_frames
    }
}

// ---------------------------------------------------------------------------
// PyCCDResult
// ---------------------------------------------------------------------------

/// Result of a continuous collision detection (CCD) query.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyCCDResult {
    /// Whether a collision was detected.
    pub hit: bool,
    /// Time of impact in \[0, 1\] (fraction of the time step).
    pub toi: f64,
    /// Contact normal at the time of impact.
    pub normal: [f64; 3],
    /// Witness point on body A.
    pub witness_a: [f64; 3],
    /// Witness point on body B.
    pub witness_b: [f64; 3],
    /// Body A handle.
    pub body_a: u32,
    /// Body B handle.
    pub body_b: u32,
}

#[pymethods]
impl PyCCDResult {
    /// No collision result.
    #[staticmethod]
    pub fn miss(body_a: u32, body_b: u32) -> Self {
        Self {
            hit: false,
            toi: 1.0,
            normal: [0.0, 1.0, 0.0],
            witness_a: [0.0; 3],
            witness_b: [0.0; 3],
            body_a,
            body_b,
        }
    }

    /// Collision result at a given time of impact.
    #[staticmethod]
    pub fn hit(
        body_a: u32,
        body_b: u32,
        toi: f64,
        normal: [f64; 3],
        witness_a: [f64; 3],
        witness_b: [f64; 3],
    ) -> Self {
        Self {
            hit: true,
            toi: toi.clamp(0.0, 1.0),
            normal,
            witness_a,
            witness_b,
            body_a,
            body_b,
        }
    }

    /// Separation distance between witness points.
    pub fn witness_distance(&self) -> f64 {
        let dx = self.witness_a[0] - self.witness_b[0];
        let dy = self.witness_a[1] - self.witness_b[1];
        let dz = self.witness_a[2] - self.witness_b[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Get hit.
    #[getter]
    pub fn get_hit(&self) -> bool {
        self.hit
    }

    /// Get toi.
    #[getter]
    pub fn get_toi(&self) -> f64 {
        self.toi
    }

    /// Get normal.
    #[getter]
    pub fn get_normal(&self) -> [f64; 3] {
        self.normal
    }

    /// Get witness_a.
    #[getter]
    pub fn get_witness_a(&self) -> [f64; 3] {
        self.witness_a
    }

    /// Get witness_b.
    #[getter]
    pub fn get_witness_b(&self) -> [f64; 3] {
        self.witness_b
    }

    /// Get body_a.
    #[getter]
    pub fn get_body_a(&self) -> u32 {
        self.body_a
    }

    /// Get body_b.
    #[getter]
    pub fn get_body_b(&self) -> u32 {
        self.body_b
    }
}

// ---------------------------------------------------------------------------
// PyPbdSolver
// ---------------------------------------------------------------------------

/// PBD/XPBD constraint type.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum PbdConstraintType {
    /// Distance / stretch constraint.
    Stretch { rest_length: f64, compliance: f64 },
    /// Bending constraint between three particles.
    Bend { rest_angle: f64, compliance: f64 },
    /// Volume conservation constraint.
    Volume { rest_volume: f64, compliance: f64 },
    /// Collision contact constraint.
    Contact { normal: [f64; 3], penetration: f64 },
}

/// A PBD/XPBD constraint connecting particle indices.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PbdConstraint {
    /// Particle indices involved (1, 2 or 3).
    pub particles: Vec<usize>,
    /// Constraint type and parameters.
    pub kind: PbdConstraintType,
    /// Accumulated Lagrange multiplier (XPBD).
    pub lambda: f64,
}

/// Position-based dynamics (XPBD) solver.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyPbdSolver {
    /// Current particle positions (flat \[x0, y0, z0, x1, ...\]).
    pub positions: Vec<f64>,
    /// Previous particle positions (for XPBD velocity estimation).
    pub prev_positions: Vec<f64>,
    /// Per-particle inverse mass (0 = kinematic/fixed).
    pub inv_masses: Vec<f64>,
    /// Constraints.
    pub constraints: Vec<PbdConstraint>,
    /// Number of substeps per time step.
    pub substeps: u32,
    /// Gravity vector.
    pub gravity: [f64; 3],
}

#[pymethods]
impl PyPbdSolver {
    /// Create a new PBD solver.
    #[new]
    pub fn new(positions: Vec<f64>, inv_masses: Vec<f64>) -> Self {
        let prev = positions.clone();
        Self {
            positions,
            prev_positions: prev,
            inv_masses,
            constraints: Vec::new(),
            substeps: 8,
            gravity: [0.0, -9.81, 0.0],
        }
    }

    /// Add a stretch constraint between two particles.
    pub fn add_stretch(&mut self, i: usize, j: usize, compliance: f64) {
        let dx = self.positions[3 * i] - self.positions[3 * j];
        let dy = self.positions[3 * i + 1] - self.positions[3 * j + 1];
        let dz = self.positions[3 * i + 2] - self.positions[3 * j + 2];
        let rest = (dx * dx + dy * dy + dz * dz).sqrt();
        self.constraints.push(PbdConstraint {
            particles: vec![i, j],
            kind: PbdConstraintType::Stretch {
                rest_length: rest,
                compliance,
            },
            lambda: 0.0,
        });
    }

    /// Add a volume conservation constraint for a tetrahedron.
    pub fn add_volume(&mut self, i: usize, j: usize, k: usize, l: usize, compliance: f64) {
        self.constraints.push(PbdConstraint {
            particles: vec![i, j, k, l],
            kind: PbdConstraintType::Volume {
                rest_volume: 1.0,
                compliance,
            },
            lambda: 0.0,
        });
    }

    /// Perform one substep of XPBD simulation.
    pub fn substep(&mut self, dt: f64) {
        let h = dt / self.substeps as f64;
        let np = self.positions.len() / 3;
        // Apply gravity and predict positions
        for i in 0..np {
            if self.inv_masses[i] > 0.0 {
                let vx = self.positions[3 * i] - self.prev_positions[3 * i];
                let vy = self.positions[3 * i + 1] - self.prev_positions[3 * i + 1];
                let vz = self.positions[3 * i + 2] - self.prev_positions[3 * i + 2];
                self.prev_positions[3 * i] = self.positions[3 * i];
                self.prev_positions[3 * i + 1] = self.positions[3 * i + 1];
                self.prev_positions[3 * i + 2] = self.positions[3 * i + 2];
                self.positions[3 * i] += vx + self.gravity[0] * h * h;
                self.positions[3 * i + 1] += vy + self.gravity[1] * h * h;
                self.positions[3 * i + 2] += vz + self.gravity[2] * h * h;
            }
        }
        // Reset lambdas
        for c in &mut self.constraints {
            c.lambda = 0.0;
        }
        // Solve stretch constraints (simplified)
        let n_constraints = self.constraints.len();
        for ci in 0..n_constraints {
            if let PbdConstraintType::Stretch {
                rest_length,
                compliance,
            } = self.constraints[ci].kind.clone()
            {
                let parts = self.constraints[ci].particles.clone();
                if parts.len() < 2 {
                    continue;
                }
                let (i, j) = (parts[0], parts[1]);
                let wi = self.inv_masses[i];
                let wj = self.inv_masses[j];
                let w_sum = wi + wj;
                if w_sum < 1e-15 {
                    continue;
                }
                let dx = self.positions[3 * i] - self.positions[3 * j];
                let dy = self.positions[3 * i + 1] - self.positions[3 * j + 1];
                let dz = self.positions[3 * i + 2] - self.positions[3 * j + 2];
                let len = (dx * dx + dy * dy + dz * dz).sqrt();
                if len < 1e-12 {
                    continue;
                }
                let alpha = compliance / (h * h);
                let c_val = len - rest_length;
                let d_lambda = (-c_val - alpha * self.constraints[ci].lambda) / (w_sum + alpha);
                self.constraints[ci].lambda += d_lambda;
                let nx = dx / len;
                let ny = dy / len;
                let nz = dz / len;
                self.positions[3 * i] += wi * d_lambda * nx;
                self.positions[3 * i + 1] += wi * d_lambda * ny;
                self.positions[3 * i + 2] += wi * d_lambda * nz;
                self.positions[3 * j] -= wj * d_lambda * nx;
                self.positions[3 * j + 1] -= wj * d_lambda * ny;
                self.positions[3 * j + 2] -= wj * d_lambda * nz;
            }
        }
    }

    /// Step the simulation for one full time step (using substeps).
    pub fn step(&mut self, dt: f64) {
        for _ in 0..self.substeps {
            self.substep(dt);
        }
    }

    /// Return particle count.
    pub fn particle_count(&self) -> usize {
        self.positions.len() / 3
    }

    /// Get positions.
    #[getter]
    pub fn get_positions(&self) -> Vec<f64> {
        self.positions.clone()
    }

    /// Set positions.
    #[setter]
    pub fn set_positions(&mut self, v: Vec<f64>) {
        self.positions = v;
    }

    /// Get prev_positions.
    #[getter]
    pub fn get_prev_positions(&self) -> Vec<f64> {
        self.prev_positions.clone()
    }

    /// Get inv_masses.
    #[getter]
    pub fn get_inv_masses(&self) -> Vec<f64> {
        self.inv_masses.clone()
    }

    /// Set inv_masses.
    #[setter]
    pub fn set_inv_masses(&mut self, v: Vec<f64>) {
        self.inv_masses = v;
    }

    /// Get substeps.
    #[getter]
    pub fn get_substeps(&self) -> u32 {
        self.substeps
    }

    /// Set substeps.
    #[setter]
    pub fn set_substeps(&mut self, v: u32) {
        self.substeps = v;
    }

    /// Get gravity.
    #[getter]
    pub fn get_gravity(&self) -> [f64; 3] {
        self.gravity
    }

    /// Set gravity.
    #[setter]
    pub fn set_gravity(&mut self, v: [f64; 3]) {
        self.gravity = v;
    }
}

// ---------------------------------------------------------------------------
// PyControlSystem
// ---------------------------------------------------------------------------

/// A simple state-space system: dx/dt = A x + B u, y = C x + D u.
#[pyclass(get_all, set_all, from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyControlSystem {
    /// System order (n).
    pub order: usize,
    /// State matrix A (n×n, row-major flat).
    pub a_matrix: Vec<f64>,
    /// Input matrix B (n×m, row-major flat).
    pub b_matrix: Vec<f64>,
    /// Output matrix C (p×n, row-major flat).
    pub c_matrix: Vec<f64>,
    /// Feedthrough matrix D (p×m, row-major flat).
    pub d_matrix: Vec<f64>,
    /// Number of inputs (m).
    pub num_inputs: usize,
    /// Number of outputs (p).
    pub num_outputs: usize,
    /// Current state vector (n values).
    pub state: Vec<f64>,
}

#[pymethods]
impl PyControlSystem {
    /// Create a first-order lag system: dx/dt = -x/tau + u/tau.
    #[staticmethod]
    pub fn first_order_lag(tau: f64) -> Self {
        Self {
            order: 1,
            a_matrix: vec![-1.0 / tau],
            b_matrix: vec![1.0 / tau],
            c_matrix: vec![1.0],
            d_matrix: vec![0.0],
            num_inputs: 1,
            num_outputs: 1,
            state: vec![0.0],
        }
    }

    /// Step the system using Euler integration.
    pub fn step_euler(&mut self, u: Vec<f64>, dt: f64) -> Vec<f64> {
        let n = self.order;
        let m = self.num_inputs;
        let p = self.num_outputs;
        // dx = A x + B u
        let mut dx = vec![0.0f64; n];
        for (i, dx_i) in dx.iter_mut().enumerate() {
            for (j, &s_j) in self.state.iter().enumerate() {
                *dx_i += self.a_matrix[i * n + j] * s_j;
            }
            for (j, u_j) in u.iter().enumerate().take(m) {
                *dx_i += self.b_matrix[i * m + j] * u_j;
            }
        }
        for (s, dx_i) in self.state.iter_mut().zip(dx.iter()) {
            *s += dx_i * dt;
        }
        // y = C x + D u
        let mut y = vec![0.0f64; p];
        for (i, y_i) in y.iter_mut().enumerate() {
            for (j, &s_j) in self.state.iter().enumerate() {
                *y_i += self.c_matrix[i * n + j] * s_j;
            }
            for (j, u_j) in u.iter().enumerate().take(m) {
                *y_i += self.d_matrix[i * m + j] * u_j;
            }
        }
        y
    }

    /// Compute step response over `n_steps` time steps.
    pub fn step_response(&mut self, dt: f64, n_steps: usize) -> Vec<f64> {
        self.state = vec![0.0; self.order];
        let mut out = Vec::with_capacity(n_steps);
        for _ in 0..n_steps {
            let y = self.step_euler(vec![1.0], dt);
            out.push(y[0]);
        }
        out
    }

    /// Evaluate the frequency response (gain) at frequency f (Hz).
    pub fn frequency_gain(&self, f: f64) -> f64 {
        // For first-order system: |G(jw)| = 1/sqrt(1 + (w*tau)^2)
        // Generic approximation: just use the DC gain for now
        let w = 2.0 * std::f64::consts::PI * f;
        // For simple first-order lag: tau = -1/a[0]
        if self.order == 1 && self.a_matrix[0] < 0.0 {
            let tau = -1.0 / self.a_matrix[0];
            let b = self.b_matrix[0] * tau;
            b / (1.0 + (w * tau) * (w * tau)).sqrt()
        } else {
            1.0
        }
    }
}

// ---------------------------------------------------------------------------
// PyFrictionModel
// ---------------------------------------------------------------------------

/// Friction model type.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum FrictionModelKind {
    /// Standard Coulomb friction.
    Coulomb,
    /// Anisotropic friction with different coefficients along u and v directions.
    Anisotropic { mu_u: f64, mu_v: f64 },
    /// Velocity-dependent friction (Stribeck curve).
    VelocityDependent {
        mu_static: f64,
        mu_kinetic: f64,
        stribeck_velocity: f64,
    },
    /// Friction cone (3D).
    Cone { mu: f64 },
}

/// Friction model for contact resolution.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyFrictionModel {
    /// Friction model type.
    pub kind: FrictionModelKind,
    /// Global friction coefficient override (used by Coulomb and Cone).
    pub mu: f64,
    /// Regularization parameter to avoid singular Jacobians at zero slip.
    pub epsilon: f64,
}

#[pymethods]
impl PyFrictionModel {
    /// Standard Coulomb friction.
    #[staticmethod]
    pub fn coulomb(mu: f64) -> Self {
        Self {
            kind: FrictionModelKind::Coulomb,
            mu,
            epsilon: 1e-4,
        }
    }

    /// Anisotropic friction.
    #[staticmethod]
    pub fn anisotropic(mu_u: f64, mu_v: f64) -> Self {
        Self {
            kind: FrictionModelKind::Anisotropic { mu_u, mu_v },
            mu: mu_u.max(mu_v),
            epsilon: 1e-4,
        }
    }

    /// Velocity-dependent (Stribeck) friction.
    #[staticmethod]
    pub fn stribeck(mu_static: f64, mu_kinetic: f64, stribeck_velocity: f64) -> Self {
        Self {
            kind: FrictionModelKind::VelocityDependent {
                mu_static,
                mu_kinetic,
                stribeck_velocity,
            },
            mu: mu_static,
            epsilon: 1e-4,
        }
    }

    /// Compute the friction force limit given normal force and slip velocity.
    pub fn friction_limit(&self, normal_force: f64, slip_speed: f64) -> f64 {
        let mu = match &self.kind {
            FrictionModelKind::Coulomb => self.mu,
            FrictionModelKind::Cone { mu } => *mu,
            FrictionModelKind::Anisotropic { mu_u, mu_v } => mu_u.max(*mu_v),
            FrictionModelKind::VelocityDependent {
                mu_static,
                mu_kinetic,
                stribeck_velocity,
            } => {
                let alpha = (-slip_speed / stribeck_velocity).exp();
                mu_kinetic + (mu_static - mu_kinetic) * alpha
            }
        };
        mu * normal_force.max(0.0)
    }

    /// Get mu.
    #[getter]
    pub fn get_mu(&self) -> f64 {
        self.mu
    }

    /// Get epsilon.
    #[getter]
    pub fn get_epsilon(&self) -> f64 {
        self.epsilon
    }
}

// ---------------------------------------------------------------------------
// PyWarmStart
// ---------------------------------------------------------------------------

/// Cached impulse from the previous frame for warm starting.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CachedImpulse {
    /// Body pair key (body_a, body_b).
    pub key: (u32, u32),
    /// Cached normal impulse.
    pub lambda_n: f64,
    /// Cached tangential impulse \[t1, t2\].
    pub lambda_t: [f64; 2],
    /// Frame counter (for aging).
    pub age: u32,
}

#[pymethods]
impl CachedImpulse {
    /// Get key as (body_a, body_b).
    #[getter]
    pub fn get_key(&self) -> (u32, u32) {
        self.key
    }

    /// Get lambda_n.
    #[getter]
    pub fn get_lambda_n(&self) -> f64 {
        self.lambda_n
    }

    /// Get lambda_t.
    #[getter]
    pub fn get_lambda_t(&self) -> [f64; 2] {
        self.lambda_t
    }

    /// Get age.
    #[getter]
    pub fn get_age(&self) -> u32 {
        self.age
    }
}

/// Warm start cache for constraint impulses.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyWarmStart {
    /// Cached impulses from the previous solve.
    pub cache: Vec<CachedImpulse>,
    /// Fraction of cached impulse to apply at the start of the next frame.
    pub aging_factor: f64,
    /// Minimum quality threshold to keep a cached impulse.
    pub quality_threshold: f64,
    /// Maximum number of frames to keep a cached impulse.
    pub max_age: u32,
}

#[pymethods]
impl PyWarmStart {
    /// Create a new warm-start cache.
    #[new]
    pub fn new() -> Self {
        Self {
            cache: Vec::new(),
            aging_factor: 0.85,
            quality_threshold: 0.5,
            max_age: 3,
        }
    }

    /// Store an impulse for the next frame.
    pub fn store(&mut self, key: (u32, u32), lambda_n: f64, lambda_t: [f64; 2]) {
        if let Some(c) = self.cache.iter_mut().find(|c| c.key == key) {
            c.lambda_n = lambda_n;
            c.lambda_t = lambda_t;
            c.age = 0;
        } else {
            self.cache.push(CachedImpulse {
                key,
                lambda_n,
                lambda_t,
                age: 0,
            });
        }
    }

    /// Retrieve the cached impulse for a body pair, scaled by aging factor.
    pub fn retrieve(&self, key: (u32, u32)) -> Option<(f64, [f64; 2])> {
        self.cache.iter().find(|c| c.key == key).map(|c| {
            (
                c.lambda_n * self.aging_factor,
                [
                    c.lambda_t[0] * self.aging_factor,
                    c.lambda_t[1] * self.aging_factor,
                ],
            )
        })
    }

    /// Age all cached impulses and remove expired ones.
    pub fn age_and_prune(&mut self) {
        for c in &mut self.cache {
            c.age += 1;
        }
        self.cache.retain(|c| c.age <= self.max_age);
    }

    /// Clear all cached impulses.
    pub fn clear(&mut self) {
        self.cache.clear();
    }

    /// Get cache.
    #[getter]
    pub fn get_cache(&self) -> Vec<CachedImpulse> {
        self.cache.clone()
    }

    /// Get aging_factor.
    #[getter]
    pub fn get_aging_factor(&self) -> f64 {
        self.aging_factor
    }

    /// Set aging_factor.
    #[setter]
    pub fn set_aging_factor(&mut self, v: f64) {
        self.aging_factor = v;
    }

    /// Get quality_threshold.
    #[getter]
    pub fn get_quality_threshold(&self) -> f64 {
        self.quality_threshold
    }

    /// Get max_age.
    #[getter]
    pub fn get_max_age(&self) -> u32 {
        self.max_age
    }
}

impl Default for PyWarmStart {
    fn default() -> Self {
        Self::new()
    }
}

/// Register all `constraints` classes into a Python sub-module.
///
/// Called from the top-level `#[pymodule]` in `lib.rs`.
pub fn register_constraints_module(
    parent: &pyo3::Bound<'_, pyo3::types::PyModule>,
) -> pyo3::PyResult<()> {
    use pyo3::prelude::PyModule;
    use pyo3::types::PyModuleMethods;
    let child = PyModule::new(parent.py(), "constraints")?;
    // Enums
    child.add_class::<SolverType>()?;
    child.add_class::<IslandSleepState>()?;
    // Structs
    child.add_class::<AxisLimits>()?;
    child.add_class::<PyConstraintSolver>()?;
    child.add_class::<PyJoint>()?;
    child.add_class::<PyContactConstraint>()?;
    child.add_class::<PidGains>()?;
    child.add_class::<MotorThermal>()?;
    child.add_class::<PyMotorConstraint>()?;
    child.add_class::<PyIsland>()?;
    child.add_class::<PyCCDResult>()?;
    child.add_class::<PyPbdSolver>()?;
    child.add_class::<PyControlSystem>()?;
    child.add_class::<PyFrictionModel>()?;
    child.add_class::<CachedImpulse>()?;
    child.add_class::<PyWarmStart>()?;
    // Free functions (Python wrappers)
    child.add_function(pyo3::wrap_pyfunction!(py_solve_pgs_iteration, &child)?)?;
    child.add_function(pyo3::wrap_pyfunction!(compute_jacobian, &child)?)?;
    child.add_function(pyo3::wrap_pyfunction!(compute_effective_mass, &child)?)?;
    child.add_function(pyo3::wrap_pyfunction!(clamp_impulse, &child)?)?;
    parent.add_submodule(&child)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_solve_pgs_basic() {
        let mut lam = vec![0.0f64];
        let rhs = [1.0];
        let diag = [2.0];
        let lo = [0.0];
        let hi = [10.0];
        let res = solve_pgs_iteration(&mut lam, &rhs, &diag, &lo, &hi);
        assert!(res >= 0.0);
        assert!(lam[0] > 0.0);
    }

    #[test]
    fn test_solve_pgs_clamped() {
        let mut lam = vec![0.0f64];
        let rhs = [100.0];
        let diag = [1.0];
        let lo = [0.0];
        let hi = [5.0];
        solve_pgs_iteration(&mut lam, &rhs, &diag, &lo, &hi);
        assert!(lam[0] <= 5.0);
    }

    #[test]
    fn test_compute_jacobian() {
        let j = compute_jacobian([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(j.len(), 12);
        assert!((j[0] - 1.0).abs() < 1e-9);
        assert!((j[6] + 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_compute_effective_mass() {
        let j = compute_jacobian([0.0; 3], [0.0; 3], [1.0, 0.0, 0.0]);
        let em = compute_effective_mass(1.0, 1.0, [1.0; 3], [1.0; 3], j);
        assert!(em > 0.0);
    }

    #[test]
    fn test_clamp_impulse_within_cone() {
        let out = clamp_impulse(10.0, [0.5, 0.5], 0.8);
        let mag = (out[0] * out[0] + out[1] * out[1]).sqrt();
        assert!(mag <= 10.0 * 0.8 + 1e-9);
    }

    #[test]
    fn test_clamp_impulse_outside_cone() {
        let out = clamp_impulse(1.0, [5.0, 5.0], 0.5);
        let mag = (out[0] * out[0] + out[1] * out[1]).sqrt();
        assert!((mag - 0.5).abs() < 1e-9);
    }

    #[test]
    fn test_constraint_solver_pgs() {
        let solver = PyConstraintSolver::default_pgs();
        assert_eq!(solver.solver_type, SolverType::Pgs);
        let lambda = solver.solve(
            vec![1.0, 1.0],
            vec![2.0, 2.0],
            vec![0.0, 0.0],
            vec![10.0, 10.0],
        );
        assert_eq!(lambda.len(), 2);
        assert!(lambda[0] >= 0.0);
    }

    #[test]
    fn test_joint_fixed() {
        let j = PyJoint::fixed(0, 1, 2, [0.0; 3], [0.0; 3]);
        assert!(!j.is_breakable());
        assert!(j.enabled);
    }

    #[test]
    fn test_joint_spring() {
        let j = PyJoint::spring(1, 0, 1, [0.0; 3], [1.0, 0.0, 0.0], 100.0, 5.0, 1.0);
        matches!(j.kind, JointKind::Spring { .. });
    }

    #[test]
    fn test_contact_constraint_clamp_friction() {
        let mut c = PyContactConstraint::new(0, 1, [0.0, 1.0, 0.0], 0.01, [0.0; 3], 0.0, 0.5);
        c.lambda_n = 10.0;
        c.lambda_t1 = 8.0;
        c.lambda_t2 = 0.0;
        c.clamp_friction();
        assert!(c.lambda_t1 <= 5.0 + 1e-9);
    }

    #[test]
    fn test_contact_constraint_tangent_basis() {
        let c = PyContactConstraint::new(0, 1, [0.0, 1.0, 0.0], 0.01, [0.0; 3], 0.0, 0.5);
        let (t1, t2) = c.tangent_basis();
        // t1 and t2 should be perpendicular to normal
        let dot1 = t1[0] * 0.0 + t1[1] * 1.0 + t1[2] * 0.0;
        assert!(dot1.abs() < 1e-9);
        let _ = t2;
    }

    #[test]
    fn test_pid_gains() {
        let mut pid = PidGains::new(1.0, 0.0, 0.0);
        let out = pid.compute(2.0, 0.01);
        assert!((out - 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_pid_reset() {
        let mut pid = PidGains::new(1.0, 1.0, 0.0);
        pid.compute(1.0, 0.1);
        pid.reset();
        assert_eq!(pid.integral, 0.0);
    }

    #[test]
    fn test_motor_thermal() {
        let mut th = MotorThermal::new();
        assert!((th.torque_scale() - 1.0).abs() < 1e-9);
        th.temperature = 100.0;
        assert!(th.torque_scale() < 1.0);
        th.temperature = 120.0;
        assert!(th.torque_scale() <= 0.0);
    }

    #[test]
    fn test_motor_constraint_step() {
        let mut motor = PyMotorConstraint::velocity_mode(10.0, 100.0, 2.0);
        let torque = motor.step(0.0, 0.01);
        assert!(torque.abs() > 0.0);
    }

    #[test]
    fn test_island_sleep() {
        let mut island = PyIsland::new(0);
        island.add_body(1);
        assert!(island.is_active());
        for _ in 0..60 {
            island.update_sleep(0.01); // below threshold
        }
        assert_eq!(island.sleep_state, IslandSleepState::Sleeping);
        island.update_sleep(1.0); // wake up
        assert_eq!(island.sleep_state, IslandSleepState::Awake);
    }

    #[test]
    fn test_island_merge() {
        let mut a = PyIsland::new(0);
        a.add_body(1);
        let mut b = PyIsland::new(1);
        b.add_body(2);
        a.merge(&b);
        assert_eq!(a.bodies.len(), 2);
    }

    #[test]
    fn test_ccd_result_miss() {
        let r = PyCCDResult::miss(0, 1);
        assert!(!r.hit);
        assert!((r.toi - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_ccd_result_hit() {
        let r = PyCCDResult::hit(0, 1, 0.3, [0.0, 1.0, 0.0], [0.0; 3], [0.0, 0.01, 0.0]);
        assert!(r.hit);
        assert!((r.toi - 0.3).abs() < 1e-9);
    }

    #[test]
    fn test_pbd_solver_stretch() {
        let positions = vec![0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        let inv_masses = vec![1.0, 1.0];
        let mut solver = PyPbdSolver::new(positions, inv_masses);
        solver.add_stretch(0, 1, 1e-4);
        solver.step(0.016);
        // Particles should have moved due to gravity
        assert!(solver.positions[1] < 0.0); // y should be negative
    }

    #[test]
    fn test_pbd_particle_count() {
        let pos = vec![0.0; 9]; // 3 particles
        let inv = vec![1.0; 3];
        let s = PyPbdSolver::new(pos, inv);
        assert_eq!(s.particle_count(), 3);
    }

    #[test]
    fn test_control_system_step_response() {
        let mut sys = PyControlSystem::first_order_lag(0.1);
        let resp = sys.step_response(0.01, 50);
        assert_eq!(resp.len(), 50);
        // Should converge toward 1.0
        assert!(*resp.last().unwrap() > 0.5);
    }

    #[test]
    fn test_control_system_frequency() {
        let sys = PyControlSystem::first_order_lag(0.1);
        let dc = sys.frequency_gain(0.0);
        assert!(dc > 0.0);
        let hf = sys.frequency_gain(100.0);
        assert!(hf < dc);
    }

    #[test]
    fn test_friction_model_coulomb() {
        let fm = PyFrictionModel::coulomb(0.5);
        assert!((fm.friction_limit(10.0, 0.0) - 5.0).abs() < 1e-9);
    }

    #[test]
    fn test_friction_model_stribeck() {
        let fm = PyFrictionModel::stribeck(1.0, 0.5, 0.1);
        let f_static = fm.friction_limit(10.0, 0.0);
        let f_kinetic = fm.friction_limit(10.0, 1.0);
        assert!(f_static >= f_kinetic);
    }

    #[test]
    fn test_warm_start_store_retrieve() {
        let mut ws = PyWarmStart::new();
        ws.store((0, 1), 5.0, [1.0, 2.0]);
        let r = ws.retrieve((0, 1));
        assert!(r.is_some());
        let (ln, _lt) = r.unwrap();
        assert!((ln - 5.0 * 0.85).abs() < 1e-9);
    }

    #[test]
    fn test_warm_start_age_prune() {
        let mut ws = PyWarmStart::new();
        ws.store((0, 1), 5.0, [0.0, 0.0]);
        for _ in 0..=ws.max_age {
            ws.age_and_prune();
        }
        assert!(ws.retrieve((0, 1)).is_none());
    }

    #[test]
    fn test_warm_start_clear() {
        let mut ws = PyWarmStart::new();
        ws.store((0, 1), 1.0, [0.0, 0.0]);
        ws.clear();
        assert!(ws.cache.is_empty());
    }
}
