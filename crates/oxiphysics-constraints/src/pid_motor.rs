// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! PID-controlled motors, advanced joint limits, spring joints, torsion springs,
//! gear constraints, cascade PID, feed-forward control, setpoint profiles,
//! motor system identification, and gain scheduling.

// -----------------------------------------------------------------------
// PidController
// -----------------------------------------------------------------------

/// Standard PID controller with optional integral clamping.
pub struct PidController {
    /// Proportional gain.
    pub kp: f64,
    /// Integral gain.
    pub ki: f64,
    /// Derivative gain.
    pub kd: f64,
    /// Accumulated integral term.
    pub integral: f64,
    /// Previous error value (for derivative).
    pub prev_error: f64,
    /// Maximum absolute value of the integral term (anti-windup).
    pub integral_clamp: f64,
}

impl PidController {
    /// Create a new PID controller with the given gains.
    /// The integral clamp defaults to `f64::INFINITY` (no clamping).
    pub fn new(kp: f64, ki: f64, kd: f64) -> Self {
        Self {
            kp,
            ki,
            kd,
            integral: 0.0,
            prev_error: 0.0,
            integral_clamp: f64::INFINITY,
        }
    }

    /// Set an integral clamp for anti-windup.
    pub fn with_integral_clamp(mut self, clamp: f64) -> Self {
        self.integral_clamp = clamp;
        self
    }

    /// Update the PID controller with the current error and time step.
    /// Returns the control output.
    pub fn update(&mut self, error: f64, dt: f64) -> f64 {
        if dt <= 0.0 {
            return 0.0;
        }
        self.integral += error * dt;
        self.integral = self
            .integral
            .clamp(-self.integral_clamp, self.integral_clamp);
        let derivative = (error - self.prev_error) / dt;
        self.prev_error = error;
        self.kp * error + self.ki * self.integral + self.kd * derivative
    }

    /// Reset the integral accumulator and previous error to zero.
    pub fn reset(&mut self) {
        self.integral = 0.0;
        self.prev_error = 0.0;
    }
}

// -----------------------------------------------------------------------
// PidMotorJoint
// -----------------------------------------------------------------------

/// A revolute motor joint driven by a PID controller.
pub struct PidMotorJoint {
    /// Embedded PID controller.
    pub pid: PidController,
    /// Maximum torque the motor can produce.
    pub max_torque: f64,
    /// Desired angle in radians.
    pub target_angle: f64,
    /// Current angle in radians.
    pub current_angle: f64,
    /// Current angular velocity in rad/s.
    pub current_velocity: f64,
}

impl PidMotorJoint {
    /// Create a new PID motor joint.
    pub fn new(kp: f64, ki: f64, kd: f64, max_torque: f64) -> Self {
        Self {
            pid: PidController::new(kp, ki, kd),
            max_torque,
            target_angle: 0.0,
            current_angle: 0.0,
            current_velocity: 0.0,
        }
    }

    /// Set the desired target angle in radians.
    pub fn set_target(&mut self, angle: f64) {
        self.target_angle = angle;
    }

    /// Compute the clamped torque output for this time step.
    pub fn update(&mut self, dt: f64) -> f64 {
        let error = self.target_angle - self.current_angle;
        let raw = self.pid.update(error, dt);
        raw.clamp(-self.max_torque, self.max_torque)
    }

    /// Advance the simulation by `dt` seconds given an external torque and
    /// moment of inertia. Updates `current_angle` and `current_velocity`.
    pub fn simulate_step(&mut self, external_torque: f64, inertia: f64, dt: f64) {
        if inertia <= 0.0 || dt <= 0.0 {
            return;
        }
        let motor_torque = self.update(dt);
        let total_torque = motor_torque + external_torque;
        let angular_acc = total_torque / inertia;
        self.current_velocity += angular_acc * dt;
        self.current_angle += self.current_velocity * dt;
    }
}

// -----------------------------------------------------------------------
// JointLimit
// -----------------------------------------------------------------------

/// Soft/hard joint limits using a penalty spring.
pub struct JointLimit {
    /// Lower bound (radians or metres).
    pub lower: f64,
    /// Upper bound (radians or metres).
    pub upper: f64,
    /// Stiffness of the penalty spring.
    pub stiffness: f64,
    /// Damping of the penalty spring.
    pub damping: f64,
}

impl JointLimit {
    /// Create a hard joint limit with no soft spring (stiffness = 1e6, damping = 0).
    pub fn new(lower: f64, upper: f64) -> Self {
        Self {
            lower,
            upper,
            stiffness: 1.0e6,
            damping: 0.0,
        }
    }

    /// Attach soft-stop parameters (spring stiffness and damping).
    pub fn with_soft(mut self, k: f64, d: f64) -> Self {
        self.stiffness = k;
        self.damping = d;
        self
    }

    /// Return a restoring force if `pos` is outside `[lower, upper]`.
    /// The force pushes back toward the bound using a penalty spring.
    pub fn limit_force(&self, pos: f64, vel: f64) -> f64 {
        if pos < self.lower {
            let penetration = self.lower - pos;
            self.stiffness * penetration - self.damping * vel
        } else if pos > self.upper {
            let penetration = pos - self.upper;
            -(self.stiffness * penetration + self.damping * vel)
        } else {
            0.0
        }
    }
}

// -----------------------------------------------------------------------
// SpringJoint
// -----------------------------------------------------------------------

/// A 1-D / 3-D spring-damper joint.
pub struct SpringJoint {
    /// Natural (rest) length of the spring.
    pub rest_length: f64,
    /// Spring stiffness coefficient.
    pub stiffness: f64,
    /// Damping coefficient.
    pub damping: f64,
}

impl SpringJoint {
    /// Create a new spring joint.
    pub fn new(rest_length: f64, stiffness: f64, damping: f64) -> Self {
        Self {
            rest_length,
            stiffness,
            damping,
        }
    }

    /// Compute the scalar spring-damper force for a 1-D spring.
    /// Positive force means compression (pushing apart).
    pub fn force_1d(&self, length: f64, velocity: f64) -> f64 {
        -self.stiffness * (length - self.rest_length) - self.damping * velocity
    }

    /// Compute the 3-D spring-damper forces on each endpoint.
    /// Returns `(force_on_a, force_on_b)` where `force_on_b = -force_on_a`.
    pub fn force_3d(
        &self,
        pa: [f64; 3],
        pb: [f64; 3],
        va: [f64; 3],
        vb: [f64; 3],
    ) -> ([f64; 3], [f64; 3]) {
        let dx = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
        let length = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
        if length < 1e-14 {
            return ([0.0; 3], [0.0; 3]);
        }
        let dir = [dx[0] / length, dx[1] / length, dx[2] / length];
        let rel_vel = [vb[0] - va[0], vb[1] - va[1], vb[2] - va[2]];
        let vel_along = rel_vel[0] * dir[0] + rel_vel[1] * dir[1] + rel_vel[2] * dir[2];
        let scalar = self.stiffness * (length - self.rest_length) + self.damping * vel_along;
        let fa = [scalar * dir[0], scalar * dir[1], scalar * dir[2]];
        let fb = [-fa[0], -fa[1], -fa[2]];
        (fa, fb)
    }
}

// -----------------------------------------------------------------------
// TorsionSpring
// -----------------------------------------------------------------------

/// A rotational spring-damper.
pub struct TorsionSpring {
    /// Natural (rest) angle in radians.
    pub rest_angle: f64,
    /// Torsional stiffness.
    pub stiffness: f64,
    /// Torsional damping.
    pub damping: f64,
}

impl TorsionSpring {
    /// Create a new torsion spring.
    pub fn new(rest_angle: f64, stiffness: f64, damping: f64) -> Self {
        Self {
            rest_angle,
            stiffness,
            damping,
        }
    }

    /// Compute the restoring torque for the given angle and angular velocity.
    pub fn torque(&self, angle: f64, angular_velocity: f64) -> f64 {
        -self.stiffness * (angle - self.rest_angle) - self.damping * angular_velocity
    }
}

// -----------------------------------------------------------------------
// GearConstraint
// -----------------------------------------------------------------------

/// A gear constraint coupling two rotational joints with a fixed ratio.
pub struct GearConstraint {
    /// Gear ratio: `v_b = ratio * v_a`.
    pub ratio: f64,
    /// Index of joint A.
    pub joint_a: usize,
    /// Index of joint B.
    pub joint_b: usize,
}

impl GearConstraint {
    /// Create a new gear constraint.
    pub fn new(joint_a: usize, joint_b: usize, ratio: f64) -> Self {
        Self {
            ratio,
            joint_a,
            joint_b,
        }
    }

    /// Return the constrained velocity of joint B given joint A's velocity.
    pub fn constrained_velocity_b(&self, vel_a: f64) -> f64 {
        self.ratio * vel_a
    }

    /// Residual of the gear constraint: `vel_b - ratio * vel_a`.
    /// Zero when the constraint is satisfied.
    pub fn residual(&self, vel_a: f64, vel_b: f64) -> f64 {
        vel_b - self.ratio * vel_a
    }
}

// -----------------------------------------------------------------------
// Tuning helpers
// -----------------------------------------------------------------------

/// Compute "classic" PID gains using the Ziegler-Nichols method.
///
/// * `ku` - ultimate gain
/// * `tu` - oscillation period at ultimate gain
///
/// Returns `(kp, ki, kd)`.
pub fn pid_tune_ziegler_nichols(ku: f64, tu: f64) -> (f64, f64, f64) {
    let kp = 0.6 * ku;
    let ti = 0.5 * tu;
    let td = 0.125 * tu;
    let ki = kp / ti;
    let kd = kp * td;
    (kp, ki, kd)
}

/// Compute IMC (Internal Model Control / Lambda) PID gains.
///
/// * `process_gain` - steady-state gain of the process
/// * `time_constant` - dominant time constant of the process
/// * `lambda`        - desired closed-loop time constant (tuning knob)
///
/// Returns `(kp, ki, kd)`.
pub fn pid_tune_lambda(process_gain: f64, time_constant: f64, lambda: f64) -> (f64, f64, f64) {
    let kp = time_constant / (process_gain * lambda);
    let ki = kp / time_constant;
    let kd = 0.0;
    (kp, ki, kd)
}

// -----------------------------------------------------------------------
// CascadePid
// -----------------------------------------------------------------------

/// Cascade PID controller with an outer (primary) and inner (secondary) loop.
///
/// The outer loop sets the setpoint for the inner loop:
///
/// ```text
/// outer_output = outer_pid(outer_error)
/// inner_setpoint = outer_output
/// inner_output = inner_pid(inner_error = inner_setpoint - inner_measurement)
/// ```
pub struct CascadePid {
    /// Outer (primary) PID controller.
    pub outer: PidController,
    /// Inner (secondary) PID controller.
    pub inner: PidController,
    /// Scale factor from outer output to inner setpoint.
    pub outer_to_inner_scale: f64,
}

impl CascadePid {
    /// Create a new cascade PID.
    pub fn new(
        outer_kp: f64,
        outer_ki: f64,
        outer_kd: f64,
        inner_kp: f64,
        inner_ki: f64,
        inner_kd: f64,
    ) -> Self {
        Self {
            outer: PidController::new(outer_kp, outer_ki, outer_kd),
            inner: PidController::new(inner_kp, inner_ki, inner_kd),
            outer_to_inner_scale: 1.0,
        }
    }

    /// Set scale factor from outer output to inner setpoint.
    pub fn with_scale(mut self, scale: f64) -> Self {
        self.outer_to_inner_scale = scale;
        self
    }

    /// Update the cascade controller.
    ///
    /// * `outer_error` - error in the primary (slow) loop
    /// * `inner_measurement` - measurement of the secondary (fast) variable
    /// * `dt` - time step
    ///
    /// Returns the final control output from the inner loop.
    pub fn update(&mut self, outer_error: f64, inner_measurement: f64, dt: f64) -> f64 {
        let outer_output = self.outer.update(outer_error, dt);
        let inner_setpoint = outer_output * self.outer_to_inner_scale;
        let inner_error = inner_setpoint - inner_measurement;
        self.inner.update(inner_error, dt)
    }

    /// Reset both controllers.
    pub fn reset(&mut self) {
        self.outer.reset();
        self.inner.reset();
    }
}

// -----------------------------------------------------------------------
// FeedForwardController
// -----------------------------------------------------------------------

/// Feed-forward + PID controller.
///
/// The total output is: u = ff_gain * setpoint + pid_output
///
/// Feed-forward provides immediate response to setpoint changes,
/// while PID handles disturbance rejection.
pub struct FeedForwardController {
    /// PID controller for feedback.
    pub pid: PidController,
    /// Feed-forward gain.
    pub ff_gain: f64,
    /// Optional feed-forward derivative gain (for velocity feed-forward).
    pub ff_derivative_gain: f64,
    /// Previous setpoint (for derivative computation).
    prev_setpoint: f64,
}

impl FeedForwardController {
    /// Create a new feed-forward + PID controller.
    pub fn new(kp: f64, ki: f64, kd: f64, ff_gain: f64) -> Self {
        Self {
            pid: PidController::new(kp, ki, kd),
            ff_gain,
            ff_derivative_gain: 0.0,
            prev_setpoint: 0.0,
        }
    }

    /// Set derivative feed-forward gain.
    pub fn with_ff_derivative(mut self, gain: f64) -> Self {
        self.ff_derivative_gain = gain;
        self
    }

    /// Compute the control output.
    ///
    /// * `setpoint` - desired value
    /// * `measurement` - current measured value
    /// * `dt` - time step
    pub fn update(&mut self, setpoint: f64, measurement: f64, dt: f64) -> f64 {
        let error = setpoint - measurement;
        let pid_out = self.pid.update(error, dt);
        let ff_out = self.ff_gain * setpoint;
        let ff_deriv = if dt > 0.0 {
            self.ff_derivative_gain * (setpoint - self.prev_setpoint) / dt
        } else {
            0.0
        };
        self.prev_setpoint = setpoint;
        pid_out + ff_out + ff_deriv
    }

    /// Reset the controller.
    pub fn reset(&mut self) {
        self.pid.reset();
        self.prev_setpoint = 0.0;
    }
}

// -----------------------------------------------------------------------
// SetpointProfile
// -----------------------------------------------------------------------

/// A setpoint profile that varies over time.
///
/// Defined by a sequence of (time, value) pairs with linear interpolation.
#[derive(Debug, Clone)]
pub struct SetpointProfile {
    /// Time-value pairs, sorted by time.
    pub points: Vec<(f64, f64)>,
}

impl SetpointProfile {
    /// Create a new profile from a list of (time, value) pairs.
    ///
    /// Points are sorted by time automatically.
    pub fn new(mut points: Vec<(f64, f64)>) -> Self {
        points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        Self { points }
    }

    /// Create a step profile: constant value before `step_time`, then `final_value`.
    pub fn step(initial_value: f64, final_value: f64, step_time: f64) -> Self {
        Self {
            points: vec![(0.0, initial_value), (step_time, final_value)],
        }
    }

    /// Create a ramp profile from (0, start) to (ramp_time, end).
    pub fn ramp(start: f64, end: f64, ramp_time: f64) -> Self {
        Self {
            points: vec![(0.0, start), (ramp_time, end)],
        }
    }

    /// Evaluate the profile at time `t` using linear interpolation.
    pub fn evaluate(&self, t: f64) -> f64 {
        if self.points.is_empty() {
            return 0.0;
        }
        if t <= self.points[0].0 {
            return self.points[0].1;
        }
        if t >= self.points[self.points.len() - 1].0 {
            return self.points[self.points.len() - 1].1;
        }
        // Find the interval
        for i in 0..(self.points.len() - 1) {
            let (t0, v0) = self.points[i];
            let (t1, v1) = self.points[i + 1];
            if t >= t0 && t <= t1 {
                let frac = (t - t0) / (t1 - t0);
                return v0 + frac * (v1 - v0);
            }
        }
        self.points[self.points.len() - 1].1
    }

    /// Duration of the profile (time of last point minus time of first point).
    pub fn duration(&self) -> f64 {
        if self.points.len() < 2 {
            return 0.0;
        }
        self.points[self.points.len() - 1].0 - self.points[0].0
    }
}

// -----------------------------------------------------------------------
// SystemIdentification
// -----------------------------------------------------------------------

/// Simple first-order system identification from step response data.
///
/// Identifies: K (steady-state gain), tau (time constant), L (dead time)
/// from the step response of a first-order plus dead-time (FOPDT) model:
///
/// ```text
/// G(s) = K * exp(-L*s) / (tau*s + 1)
/// ```
#[derive(Debug, Clone)]
pub struct SystemIdentification {
    /// Steady-state gain K.
    pub gain: f64,
    /// Time constant tau.
    pub time_constant: f64,
    /// Dead time (transport delay) L.
    pub dead_time: f64,
}

impl SystemIdentification {
    /// Create from known parameters.
    pub fn new(gain: f64, time_constant: f64, dead_time: f64) -> Self {
        Self {
            gain,
            time_constant,
            dead_time,
        }
    }

    /// Identify FOPDT model from step response data.
    ///
    /// * `times` - time points
    /// * `values` - measured response values
    /// * `step_magnitude` - size of the input step
    /// * `initial_value` - value before the step
    ///
    /// Uses the 63.2% method for time constant estimation.
    pub fn from_step_response(
        times: &[f64],
        values: &[f64],
        step_magnitude: f64,
        initial_value: f64,
    ) -> Self {
        assert_eq!(times.len(), values.len());
        if times.is_empty() || step_magnitude.abs() < 1e-20 {
            return Self::new(1.0, 1.0, 0.0);
        }

        // Steady-state gain
        let final_value = values[values.len() - 1];
        let gain = (final_value - initial_value) / step_magnitude;

        // Find dead time: first time value starts changing significantly
        let threshold = initial_value + 0.05 * (final_value - initial_value);
        let mut dead_time = 0.0;
        for i in 0..values.len() {
            if (values[i] - initial_value).abs() > (threshold - initial_value).abs() {
                dead_time = times[i];
                break;
            }
        }

        // Find time constant: time to reach 63.2% of final value
        let target_63 = initial_value + 0.632 * (final_value - initial_value);
        let mut time_constant = 1.0;
        for i in 0..values.len() {
            if (values[i] - initial_value).abs() >= (target_63 - initial_value).abs() {
                time_constant = (times[i] - dead_time).max(0.01);
                break;
            }
        }

        Self {
            gain,
            time_constant,
            dead_time,
        }
    }

    /// Compute recommended PID gains using the Cohen-Coon method.
    ///
    /// Returns `(kp, ki, kd)`.
    pub fn cohen_coon_pid(&self) -> (f64, f64, f64) {
        let r = self.dead_time / self.time_constant;
        if r < 1e-10 {
            // No dead time: just use lambda tuning
            return pid_tune_lambda(self.gain, self.time_constant, self.time_constant);
        }
        let kp = (1.0 / (self.gain * r)) * (1.33 + r / 4.0);
        let ti = self.dead_time * (32.0 + 6.0 * r) / (13.0 + 8.0 * r);
        let td = self.dead_time * 4.0 / (11.0 + 2.0 * r);
        let ki = kp / ti;
        let kd = kp * td;
        (kp, ki, kd)
    }

    /// Simulate the FOPDT step response at time `t`.
    ///
    /// ```text
    /// y(t) = K * (1 - exp(-(t - L) / tau))  for t > L
    /// y(t) = 0                                for t <= L
    /// ```
    pub fn step_response(&self, t: f64) -> f64 {
        if t <= self.dead_time {
            return 0.0;
        }
        self.gain * (1.0 - (-(t - self.dead_time) / self.time_constant).exp())
    }
}

// -----------------------------------------------------------------------
// GainScheduler
// -----------------------------------------------------------------------

/// Gain scheduler that selects PID parameters based on an operating condition.
///
/// Uses a lookup table of (condition, kp, ki, kd) with linear interpolation.
#[derive(Debug, Clone)]
pub struct GainScheduler {
    /// Table of (operating_point, kp, ki, kd), sorted by operating_point.
    pub schedule: Vec<(f64, f64, f64, f64)>,
}

impl GainScheduler {
    /// Create a new gain scheduler from a schedule table.
    pub fn new(mut schedule: Vec<(f64, f64, f64, f64)>) -> Self {
        schedule.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        Self { schedule }
    }

    /// Look up gains for the given operating condition using linear interpolation.
    ///
    /// Returns `(kp, ki, kd)`.
    pub fn lookup(&self, condition: f64) -> (f64, f64, f64) {
        if self.schedule.is_empty() {
            return (1.0, 0.0, 0.0);
        }
        if condition <= self.schedule[0].0 {
            let s = &self.schedule[0];
            return (s.1, s.2, s.3);
        }
        let last = self.schedule.len() - 1;
        if condition >= self.schedule[last].0 {
            let s = &self.schedule[last];
            return (s.1, s.2, s.3);
        }
        // Find interval and interpolate
        for i in 0..last {
            let (c0, kp0, ki0, kd0) = self.schedule[i];
            let (c1, kp1, ki1, kd1) = self.schedule[i + 1];
            if condition >= c0 && condition <= c1 {
                let frac = (condition - c0) / (c1 - c0);
                return (
                    kp0 + frac * (kp1 - kp0),
                    ki0 + frac * (ki1 - ki0),
                    kd0 + frac * (kd1 - kd0),
                );
            }
        }
        let s = &self.schedule[last];
        (s.1, s.2, s.3)
    }

    /// Apply the scheduled gains to a PID controller.
    pub fn apply_to(&self, pid: &mut PidController, condition: f64) {
        let (kp, ki, kd) = self.lookup(condition);
        pid.kp = kp;
        pid.ki = ki;
        pid.kd = kd;
    }

    /// Number of entries in the schedule.
    pub fn num_entries(&self) -> usize {
        self.schedule.len()
    }
}

// -----------------------------------------------------------------------
// MotorModel (for system identification)
// -----------------------------------------------------------------------

/// Simple DC motor model for simulation and identification.
///
/// ```text
/// J * d(omega)/dt = K_t * i - B * omega - T_load
/// L * di/dt = V - R*i - K_e * omega
/// ```
///
/// Simplified to first-order (ignoring electrical dynamics):
/// ```text
/// J * d(omega)/dt = (K_t/R) * V - (B + K_t*K_e/R) * omega - T_load
/// ```
#[derive(Debug, Clone)]
pub struct MotorModel {
    /// Moment of inertia (kg*m^2).
    pub inertia: f64,
    /// Torque constant (N*m/A).
    pub torque_constant: f64,
    /// Back-EMF constant (V*s/rad).
    pub back_emf_constant: f64,
    /// Armature resistance (Ohm).
    pub resistance: f64,
    /// Viscous friction coefficient (N*m*s/rad).
    pub friction: f64,
    /// Current angular velocity (rad/s).
    pub omega: f64,
    /// Current angle (rad).
    pub angle: f64,
}

impl MotorModel {
    /// Create a new motor model.
    pub fn new(
        inertia: f64,
        torque_constant: f64,
        back_emf_constant: f64,
        resistance: f64,
        friction: f64,
    ) -> Self {
        Self {
            inertia,
            torque_constant,
            back_emf_constant,
            resistance,
            friction,
            omega: 0.0,
            angle: 0.0,
        }
    }

    /// Compute the effective gain K_eff = K_t / R.
    pub fn effective_gain(&self) -> f64 {
        self.torque_constant / self.resistance
    }

    /// Compute the effective damping B_eff = B + K_t*K_e/R.
    pub fn effective_damping(&self) -> f64 {
        self.friction + self.torque_constant * self.back_emf_constant / self.resistance
    }

    /// Steady-state speed for a given voltage (no load).
    pub fn steady_state_speed(&self, voltage: f64) -> f64 {
        let b_eff = self.effective_damping();
        if b_eff < 1e-20 {
            return f64::INFINITY;
        }
        self.effective_gain() * voltage / b_eff
    }

    /// Simulate one time step with applied voltage and load torque.
    pub fn step(&mut self, voltage: f64, load_torque: f64, dt: f64) {
        if self.inertia <= 0.0 || dt <= 0.0 {
            return;
        }
        let torque_motor = self.effective_gain() * voltage;
        let torque_damping = self.effective_damping() * self.omega;
        let alpha = (torque_motor - torque_damping - load_torque) / self.inertia;
        self.omega += alpha * dt;
        self.angle += self.omega * dt;
    }

    /// Time constant of the motor (J / B_eff).
    pub fn time_constant(&self) -> f64 {
        let b_eff = self.effective_damping();
        if b_eff < 1e-20 {
            return f64::INFINITY;
        }
        self.inertia / b_eff
    }
}

// -----------------------------------------------------------------------
// Tests
// -----------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -- PidController --------------------------------------------------

    #[test]
    fn pid_new_gains() {
        let pid = PidController::new(1.0, 0.5, 0.1);
        assert_eq!(pid.kp, 1.0);
        assert_eq!(pid.ki, 0.5);
        assert_eq!(pid.kd, 0.1);
    }

    #[test]
    fn pid_default_integral_clamp_is_inf() {
        let pid = PidController::new(1.0, 0.0, 0.0);
        assert!(pid.integral_clamp.is_infinite());
    }

    #[test]
    fn pid_with_integral_clamp() {
        let pid = PidController::new(1.0, 1.0, 0.0).with_integral_clamp(10.0);
        assert_eq!(pid.integral_clamp, 10.0);
    }

    #[test]
    fn pid_pure_proportional() {
        let mut pid = PidController::new(2.0, 0.0, 0.0);
        let out = pid.update(3.0, 0.1);
        assert!((out - 6.0).abs() < 1e-10);
    }

    #[test]
    fn pid_integral_accumulates() {
        let mut pid = PidController::new(0.0, 1.0, 0.0);
        pid.update(1.0, 0.1);
        pid.update(1.0, 0.1);
        let out = pid.update(1.0, 0.1);
        assert!((pid.integral - 0.3).abs() < 1e-10);
        assert!((out - 0.3).abs() < 1e-10);
    }

    #[test]
    fn pid_integral_clamped() {
        let mut pid = PidController::new(0.0, 1.0, 0.0).with_integral_clamp(0.15);
        pid.update(1.0, 0.1);
        pid.update(1.0, 0.1);
        assert!((pid.integral - 0.15).abs() < 1e-10);
    }

    #[test]
    fn pid_derivative_term() {
        let mut pid = PidController::new(0.0, 0.0, 1.0);
        pid.update(0.0, 0.1);
        let out = pid.update(1.0, 0.1);
        assert!((out - 10.0).abs() < 1e-10);
    }

    #[test]
    fn pid_reset_clears_state() {
        let mut pid = PidController::new(1.0, 1.0, 1.0);
        pid.update(5.0, 0.1);
        pid.reset();
        assert_eq!(pid.integral, 0.0);
        assert_eq!(pid.prev_error, 0.0);
    }

    #[test]
    fn pid_zero_dt_returns_zero() {
        let mut pid = PidController::new(1.0, 1.0, 1.0);
        let out = pid.update(1.0, 0.0);
        assert_eq!(out, 0.0);
    }

    // -- PidMotorJoint --------------------------------------------------

    #[test]
    fn motor_joint_new_defaults() {
        let mj = PidMotorJoint::new(1.0, 0.0, 0.0, 10.0);
        assert_eq!(mj.target_angle, 0.0);
        assert_eq!(mj.current_angle, 0.0);
        assert_eq!(mj.max_torque, 10.0);
    }

    #[test]
    fn motor_joint_set_target() {
        let mut mj = PidMotorJoint::new(1.0, 0.0, 0.0, 10.0);
        mj.set_target(1.5);
        assert_eq!(mj.target_angle, 1.5);
    }

    #[test]
    fn motor_joint_torque_clamped() {
        let mut mj = PidMotorJoint::new(1000.0, 0.0, 0.0, 5.0);
        mj.set_target(1.0);
        let torque = mj.update(0.01);
        assert!((torque - 5.0).abs() < 1e-10);
    }

    #[test]
    fn motor_joint_negative_torque_clamped() {
        let mut mj = PidMotorJoint::new(1000.0, 0.0, 0.0, 5.0);
        mj.set_target(-1.0);
        let torque = mj.update(0.01);
        assert!((torque + 5.0).abs() < 1e-10);
    }

    #[test]
    fn motor_joint_simulate_step_moves_angle() {
        let mut mj = PidMotorJoint::new(10.0, 0.0, 0.0, 100.0);
        mj.set_target(1.0);
        mj.simulate_step(0.0, 1.0, 0.01);
        assert!(mj.current_angle > 0.0);
    }

    #[test]
    fn motor_joint_simulate_step_zero_inertia_no_crash() {
        let mut mj = PidMotorJoint::new(10.0, 0.0, 0.0, 100.0);
        mj.set_target(1.0);
        mj.simulate_step(0.0, 0.0, 0.01);
    }

    // -- JointLimit -----------------------------------------------------

    #[test]
    fn joint_limit_inside_bounds() {
        let lim = JointLimit::new(-1.0, 1.0);
        assert_eq!(lim.limit_force(0.0, 0.0), 0.0);
    }

    #[test]
    fn joint_limit_below_lower() {
        let lim = JointLimit::new(-1.0, 1.0).with_soft(100.0, 0.0);
        let f = lim.limit_force(-2.0, 0.0);
        assert!(f > 0.0, "force should be positive to push back");
    }

    #[test]
    fn joint_limit_above_upper() {
        let lim = JointLimit::new(-1.0, 1.0).with_soft(100.0, 0.0);
        let f = lim.limit_force(2.0, 0.0);
        assert!(f < 0.0, "force should be negative to push back");
    }

    #[test]
    fn joint_limit_soft_stop_magnitude() {
        let lim = JointLimit::new(0.0, 1.0).with_soft(50.0, 10.0);
        let f = lim.limit_force(-0.5, 1.0);
        assert!((f - (50.0 * 0.5 - 10.0 * 1.0)).abs() < 1e-10);
    }

    // -- SpringJoint ----------------------------------------------------

    #[test]
    fn spring_1d_at_rest() {
        let s = SpringJoint::new(1.0, 100.0, 0.0);
        assert_eq!(s.force_1d(1.0, 0.0), 0.0);
    }

    #[test]
    fn spring_1d_compression() {
        let s = SpringJoint::new(1.0, 100.0, 0.0);
        let f = s.force_1d(0.5, 0.0);
        assert!((f - 50.0).abs() < 1e-10);
    }

    #[test]
    fn spring_1d_extension() {
        let s = SpringJoint::new(1.0, 100.0, 0.0);
        let f = s.force_1d(1.5, 0.0);
        assert!((f + 50.0).abs() < 1e-10);
    }

    #[test]
    fn spring_3d_action_reaction() {
        let s = SpringJoint::new(1.0, 100.0, 10.0);
        let pa = [0.0, 0.0, 0.0];
        let pb = [2.0, 0.0, 0.0];
        let (fa, fb) = s.force_3d(pa, pb, [0.0; 3], [0.0; 3]);
        assert!((fa[0] + fb[0]).abs() < 1e-10);
        assert!((fa[1] + fb[1]).abs() < 1e-10);
        assert!((fa[2] + fb[2]).abs() < 1e-10);
    }

    #[test]
    fn spring_3d_at_rest_length_no_spring_force() {
        let s = SpringJoint::new(2.0, 100.0, 0.0);
        let pa = [0.0, 0.0, 0.0];
        let pb = [2.0, 0.0, 0.0];
        let (fa, _) = s.force_3d(pa, pb, [0.0; 3], [0.0; 3]);
        assert!(fa[0].abs() < 1e-10);
    }

    #[test]
    fn spring_3d_zero_length_no_crash() {
        let s = SpringJoint::new(1.0, 100.0, 10.0);
        let (fa, fb) = s.force_3d([0.0; 3], [0.0; 3], [0.0; 3], [0.0; 3]);
        assert_eq!(fa, [0.0; 3]);
        assert_eq!(fb, [0.0; 3]);
    }

    // -- TorsionSpring --------------------------------------------------

    #[test]
    fn torsion_at_rest() {
        let ts = TorsionSpring::new(0.0, 10.0, 1.0);
        assert_eq!(ts.torque(0.0, 0.0), 0.0);
    }

    #[test]
    fn torsion_restoring_torque() {
        let ts = TorsionSpring::new(0.0, 10.0, 0.0);
        let t = ts.torque(1.0, 0.0);
        assert!((t + 10.0).abs() < 1e-10);
    }

    #[test]
    fn torsion_damping_term() {
        let ts = TorsionSpring::new(0.0, 0.0, 5.0);
        let t = ts.torque(0.0, 2.0);
        assert!((t + 10.0).abs() < 1e-10);
    }

    // -- GearConstraint -------------------------------------------------

    #[test]
    fn gear_constrained_velocity() {
        let g = GearConstraint::new(0, 1, 2.0);
        assert!((g.constrained_velocity_b(3.0) - 6.0).abs() < 1e-10);
    }

    #[test]
    fn gear_residual_satisfied() {
        let g = GearConstraint::new(0, 1, 2.0);
        assert!((g.residual(3.0, 6.0)).abs() < 1e-10);
    }

    #[test]
    fn gear_residual_violation() {
        let g = GearConstraint::new(0, 1, 2.0);
        let r = g.residual(3.0, 5.0);
        assert!((r + 1.0).abs() < 1e-10);
    }

    #[test]
    fn gear_indices_stored() {
        let g = GearConstraint::new(3, 7, 0.5);
        assert_eq!(g.joint_a, 3);
        assert_eq!(g.joint_b, 7);
    }

    // -- Tuning helpers -------------------------------------------------

    #[test]
    fn ziegler_nichols_classic() {
        let (kp, ki, kd) = pid_tune_ziegler_nichols(10.0, 2.0);
        assert!((kp - 6.0).abs() < 1e-10);
        assert!((ki - 6.0).abs() < 1e-10);
        assert!((kd - 1.5).abs() < 1e-10);
    }

    #[test]
    fn lambda_tuning_basic() {
        let (kp, ki, kd) = pid_tune_lambda(2.0, 4.0, 1.0);
        assert!((kp - 2.0).abs() < 1e-10);
        assert!((ki - 0.5).abs() < 1e-10);
        assert_eq!(kd, 0.0);
    }

    #[test]
    fn lambda_tuning_kd_zero() {
        let (_, _, kd) = pid_tune_lambda(1.0, 1.0, 1.0);
        assert_eq!(kd, 0.0);
    }

    // -- CascadePid tests -----------------------------------------------

    #[test]
    fn cascade_pid_basic() {
        let mut cascade = CascadePid::new(1.0, 0.0, 0.0, 2.0, 0.0, 0.0);
        // Outer error = 1.0 -> outer output = 1.0
        // Inner setpoint = 1.0, inner measurement = 0.0 -> inner error = 1.0
        // Inner output = 2.0 * 1.0 = 2.0
        let out = cascade.update(1.0, 0.0, 0.1);
        assert!(
            (out - 2.0).abs() < 1e-10,
            "cascade output = {out}, expected 2.0"
        );
    }

    #[test]
    fn cascade_pid_reset() {
        let mut cascade = CascadePid::new(1.0, 1.0, 0.0, 1.0, 1.0, 0.0);
        cascade.update(1.0, 0.0, 0.1);
        cascade.reset();
        assert_eq!(cascade.outer.integral, 0.0);
        assert_eq!(cascade.inner.integral, 0.0);
    }

    #[test]
    fn cascade_pid_with_scale() {
        let mut cascade = CascadePid::new(1.0, 0.0, 0.0, 1.0, 0.0, 0.0).with_scale(0.5);
        // Outer output = 1.0, scaled by 0.5 -> inner setpoint = 0.5
        // Inner output = 1.0 * 0.5 = 0.5
        let out = cascade.update(1.0, 0.0, 0.1);
        assert!(
            (out - 0.5).abs() < 1e-10,
            "scaled cascade output = {out}, expected 0.5"
        );
    }

    // -- FeedForwardController tests ------------------------------------

    #[test]
    fn feedforward_basic() {
        let mut ff = FeedForwardController::new(1.0, 0.0, 0.0, 0.5);
        // PID output = 1.0 * error = 1.0 * (2.0 - 0.0) = 2.0
        // FF output = 0.5 * setpoint = 0.5 * 2.0 = 1.0
        // Total = 3.0
        let out = ff.update(2.0, 0.0, 0.1);
        assert!((out - 3.0).abs() < 1e-10, "FF output = {out}, expected 3.0");
    }

    #[test]
    fn feedforward_with_derivative() {
        let mut ff = FeedForwardController::new(0.0, 0.0, 0.0, 0.0).with_ff_derivative(1.0);
        ff.prev_setpoint = 0.0;
        // setpoint changes from 0 to 1 in dt=0.1
        // ff_deriv = 1.0 * (1.0 - 0.0) / 0.1 = 10.0
        let out = ff.update(1.0, 0.0, 0.1);
        assert!(
            (out - 10.0).abs() < 1e-10,
            "FF derivative output = {out}, expected 10.0"
        );
    }

    #[test]
    fn feedforward_reset() {
        let mut ff = FeedForwardController::new(1.0, 1.0, 0.0, 1.0);
        ff.update(1.0, 0.0, 0.1);
        ff.reset();
        assert_eq!(ff.pid.integral, 0.0);
        assert_eq!(ff.prev_setpoint, 0.0);
    }

    // -- SetpointProfile tests ------------------------------------------

    #[test]
    fn setpoint_profile_step() {
        let sp = SetpointProfile::step(0.0, 1.0, 5.0);
        assert!((sp.evaluate(0.0) - 0.0).abs() < 1e-14);
        assert!((sp.evaluate(10.0) - 1.0).abs() < 1e-14);
        // At step time: interpolated
        assert!((sp.evaluate(5.0) - 1.0).abs() < 1e-14);
    }

    #[test]
    fn setpoint_profile_ramp() {
        let sp = SetpointProfile::ramp(0.0, 10.0, 10.0);
        let mid = sp.evaluate(5.0);
        assert!(
            (mid - 5.0).abs() < 1e-10,
            "ramp midpoint = {mid}, expected 5.0"
        );
    }

    #[test]
    fn setpoint_profile_interpolation() {
        let sp = SetpointProfile::new(vec![(0.0, 0.0), (10.0, 100.0), (20.0, 50.0)]);
        assert!((sp.evaluate(5.0) - 50.0).abs() < 1e-10);
        assert!((sp.evaluate(15.0) - 75.0).abs() < 1e-10);
    }

    #[test]
    fn setpoint_profile_duration() {
        let sp = SetpointProfile::new(vec![(1.0, 0.0), (5.0, 10.0)]);
        assert!((sp.duration() - 4.0).abs() < 1e-14);
    }

    #[test]
    fn setpoint_profile_empty() {
        let sp = SetpointProfile::new(vec![]);
        assert!((sp.evaluate(1.0)).abs() < 1e-14);
    }

    // -- SystemIdentification tests -------------------------------------

    #[test]
    fn sysid_step_response() {
        let si = SystemIdentification::new(2.0, 1.0, 0.5);
        // At t <= 0.5 (dead time): response = 0
        assert!((si.step_response(0.3)).abs() < 1e-14);
        // At t >> dead_time + tau: response -> K = 2.0
        let resp = si.step_response(100.0);
        assert!(
            (resp - 2.0).abs() < 1e-6,
            "step response at t=100 should be ~2.0, got {resp}"
        );
    }

    #[test]
    fn sysid_from_step_response_data() {
        // Generate synthetic FOPDT data: K=2, tau=1, L=0
        let k = 2.0;
        let tau = 1.0;
        let times: Vec<f64> = (0..100).map(|i| i as f64 * 0.1).collect();
        let values: Vec<f64> = times
            .iter()
            .map(|&t| k * (1.0 - (-t / tau).exp()))
            .collect();
        let si = SystemIdentification::from_step_response(&times, &values, 1.0, 0.0);
        assert!(
            (si.gain - k).abs() < 0.5,
            "identified gain = {}, expected ~{k}",
            si.gain
        );
    }

    #[test]
    fn sysid_cohen_coon() {
        let si = SystemIdentification::new(2.0, 1.0, 0.5);
        let (kp, ki, kd) = si.cohen_coon_pid();
        assert!(kp > 0.0, "Cohen-Coon kp should be positive");
        assert!(ki > 0.0, "Cohen-Coon ki should be positive");
        assert!(kd > 0.0, "Cohen-Coon kd should be positive");
    }

    // -- GainScheduler tests --------------------------------------------

    #[test]
    fn gain_scheduler_lookup_exact() {
        let gs = GainScheduler::new(vec![(0.0, 1.0, 0.1, 0.01), (1.0, 2.0, 0.2, 0.02)]);
        let (kp, ki, kd) = gs.lookup(0.0);
        assert!((kp - 1.0).abs() < 1e-14);
        assert!((ki - 0.1).abs() < 1e-14);
        assert!((kd - 0.01).abs() < 1e-14);
    }

    #[test]
    fn gain_scheduler_lookup_interpolated() {
        let gs = GainScheduler::new(vec![(0.0, 1.0, 0.0, 0.0), (10.0, 3.0, 1.0, 0.5)]);
        let (kp, ki, kd) = gs.lookup(5.0);
        assert!((kp - 2.0).abs() < 1e-10, "kp = {kp}, expected 2.0");
        assert!((ki - 0.5).abs() < 1e-10, "ki = {ki}, expected 0.5");
        assert!((kd - 0.25).abs() < 1e-10, "kd = {kd}, expected 0.25");
    }

    #[test]
    fn gain_scheduler_lookup_below_range() {
        let gs = GainScheduler::new(vec![(5.0, 10.0, 1.0, 0.1)]);
        let (kp, _, _) = gs.lookup(0.0);
        assert!(
            (kp - 10.0).abs() < 1e-14,
            "below range should use first entry"
        );
    }

    #[test]
    fn gain_scheduler_apply_to() {
        let gs = GainScheduler::new(vec![(0.0, 5.0, 0.5, 0.05)]);
        let mut pid = PidController::new(1.0, 0.0, 0.0);
        gs.apply_to(&mut pid, 0.0);
        assert!((pid.kp - 5.0).abs() < 1e-14);
        assert!((pid.ki - 0.5).abs() < 1e-14);
        assert!((pid.kd - 0.05).abs() < 1e-14);
    }

    // -- MotorModel tests -----------------------------------------------

    #[test]
    fn motor_model_steady_state() {
        let motor = MotorModel::new(0.01, 0.1, 0.1, 1.0, 0.01);
        let ss = motor.steady_state_speed(12.0);
        // K_eff = 0.1/1.0 = 0.1, B_eff = 0.01 + 0.1*0.1/1.0 = 0.02
        // ss = 0.1 * 12.0 / 0.02 = 60.0
        assert!(
            (ss - 60.0).abs() < 1e-10,
            "steady state speed = {ss}, expected 60.0"
        );
    }

    #[test]
    fn motor_model_time_constant() {
        let motor = MotorModel::new(0.01, 0.1, 0.1, 1.0, 0.01);
        let tc = motor.time_constant();
        // J / B_eff = 0.01 / 0.02 = 0.5
        assert!(
            (tc - 0.5).abs() < 1e-10,
            "time constant = {tc}, expected 0.5"
        );
    }

    #[test]
    fn motor_model_step_simulation() {
        let mut motor = MotorModel::new(0.01, 0.1, 0.1, 1.0, 0.01);
        // Apply voltage and simulate for 5 seconds (10x time constant = 5s)
        for _ in 0..5000 {
            motor.step(12.0, 0.0, 0.001);
        }
        let ss = motor.steady_state_speed(12.0);
        assert!(
            (motor.omega - ss).abs() / ss < 0.01,
            "after 5s, omega={} should be near ss={}",
            motor.omega,
            ss
        );
    }

    #[test]
    fn motor_model_angle_increases() {
        let mut motor = MotorModel::new(0.01, 0.1, 0.1, 1.0, 0.01);
        motor.step(12.0, 0.0, 0.01);
        assert!(
            motor.angle > 0.0,
            "angle should increase with positive voltage"
        );
    }
}
