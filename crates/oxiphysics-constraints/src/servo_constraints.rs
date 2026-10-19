// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Servo motor constraints for physics simulation.
//!
//! Provides PID controllers, impedance control, cascade control, feedforward
//! compensation, adaptive auto-tuning, and full multi-axis servo actuators.

// ── Servo mode ────────────────────────────────────────────────────────────────

/// Operating mode of a servo axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ServoMode {
    /// Control joint position to a target angle / distance.
    Position,
    /// Control joint velocity to a target speed.
    Velocity,
    /// Output a constant torque / force.
    Torque,
    /// Control via a prescribed force trajectory.
    ForceControl,
    /// Impedance (mass–damper–spring) control.
    Impedance,
}

// ── Servo parameters ──────────────────────────────────────────────────────────

/// PID and mechanical parameters for one servo axis.
#[derive(Debug, Clone)]
pub struct ServoParams {
    /// Proportional gain.
    pub kp: f64,
    /// Integral gain.
    pub ki: f64,
    /// Derivative gain.
    pub kd: f64,
    /// Maximum output force / torque (N or N·m).
    pub max_force: f64,
    /// Maximum joint velocity (m/s or rad/s).
    pub max_velocity: f64,
    /// Gear ratio (output / input).
    pub gear_ratio: f64,
    /// Viscous damping coefficient.
    pub damping: f64,
    /// Spring stiffness (for impedance mode).
    pub stiffness: f64,
}

impl ServoParams {
    /// Create default servo parameters.
    pub fn new(kp: f64, ki: f64, kd: f64, max_force: f64, max_velocity: f64) -> Self {
        Self {
            kp,
            ki,
            kd,
            max_force,
            max_velocity,
            gear_ratio: 1.0,
            damping: 0.0,
            stiffness: 0.0,
        }
    }
}

impl Default for ServoParams {
    fn default() -> Self {
        Self::new(100.0, 10.0, 5.0, 1000.0, 10.0)
    }
}

// ── Servo state ───────────────────────────────────────────────────────────────

/// Runtime state of a servo axis.
#[derive(Debug, Clone, Default)]
pub struct ServoState {
    /// Desired target position / velocity / torque.
    pub target: f64,
    /// Current measured value.
    pub current: f64,
    /// Current tracking error (`target - current`).
    pub error: f64,
    /// Accumulated integral of the error.
    pub integral: f64,
    /// Derivative of the error (filtered).
    pub derivative: f64,
    /// Most recent output force / torque.
    pub output_force: f64,
    /// Whether the output is currently clamped by the force limit.
    pub saturated: bool,
}

// ── ServoConstraint ───────────────────────────────────────────────────────────

/// A single-axis servo constraint connecting body `body_a` to body `body_b`.
#[derive(Debug, Clone)]
pub struct ServoConstraint {
    /// Index of body A.
    pub body_a: usize,
    /// Index of body B.
    pub body_b: usize,
    /// Constraint axis in world space (unit vector).
    pub axis: [f64; 3],
    /// Servo operating mode.
    pub mode: ServoMode,
    /// Servo parameters.
    pub params: ServoParams,
    /// Current runtime state.
    pub state: ServoState,
    /// Internal PID controller.
    pid: PidController,
}

impl ServoConstraint {
    /// Construct a new `ServoConstraint`.
    pub fn new(
        body_a: usize,
        body_b: usize,
        axis: [f64; 3],
        mode: ServoMode,
        params: ServoParams,
    ) -> Self {
        let pid = PidController::new(params.kp, params.ki, params.kd);
        Self {
            body_a,
            body_b,
            axis,
            mode,
            params,
            state: ServoState::default(),
            pid,
        }
    }

    /// Set the servo target (position, velocity, or torque depending on mode).
    pub fn set_target(&mut self, target: f64) {
        self.state.target = target;
        self.pid.reset();
    }

    /// Compute the impulse to apply for this time step.
    ///
    /// * `dt` — time step (s)
    /// * `current` — current measured position / velocity
    /// * `current_vel` — current joint velocity (used in velocity-mode)
    pub fn compute_impulse(&mut self, dt: f64, current: f64, current_vel: f64) -> f64 {
        self.state.current = current;
        self.state.error = self.state.target - current;

        let raw = match self.mode {
            ServoMode::Position => {
                self.pid.step(self.state.error, dt) - self.params.damping * current_vel
            }
            ServoMode::Velocity => {
                let vel_error = self.state.target - current_vel;
                self.pid.step(vel_error, dt)
            }
            ServoMode::Torque => self.state.target,
            ServoMode::ForceControl => self.state.target,
            ServoMode::Impedance => {
                // Spring-damper: F = k*e + d*vel_e
                let vel_error = -current_vel; // target vel assumed 0 for now
                self.params.stiffness * self.state.error + self.params.damping * vel_error
            }
        };

        // Clamp to max force
        let clamped = raw.max(-self.params.max_force).min(self.params.max_force);
        self.state.saturated = raw.abs() > self.params.max_force;
        self.state.output_force = clamped;

        // Anti-windup when saturated
        if self.state.saturated {
            self.pid
                .anti_windup(raw, -self.params.max_force, self.params.max_force);
        }

        // Scale by gear ratio and dt for impulse
        clamped * self.params.gear_ratio * dt
    }

    /// Check whether the servo has reached its target within tolerance `tol`.
    pub fn is_at_target(&self, tol: f64) -> bool {
        self.state.error.abs() <= tol
    }
}

// ── PID Controller ────────────────────────────────────────────────────────────

/// A standard discrete-time PID controller with anti-windup.
#[derive(Debug, Clone)]
pub struct PidController {
    /// Proportional gain.
    pub kp: f64,
    /// Integral gain.
    pub ki: f64,
    /// Derivative gain.
    pub kd: f64,
    /// Accumulated integral term.
    integral: f64,
    /// Previous error sample.
    prev_error: f64,
    /// Output lower limit.
    min_output: f64,
    /// Output upper limit.
    max_output: f64,
    /// Whether limits have been set.
    limits_set: bool,
}

impl PidController {
    /// Create a new PID controller with given gains.
    pub fn new(kp: f64, ki: f64, kd: f64) -> Self {
        Self {
            kp,
            ki,
            kd,
            integral: 0.0,
            prev_error: 0.0,
            min_output: f64::NEG_INFINITY,
            max_output: f64::INFINITY,
            limits_set: false,
        }
    }

    /// Compute the controller output for error `error` over time step `dt`.
    pub fn step(&mut self, error: f64, dt: f64) -> f64 {
        self.integral += error * dt;
        let derivative = if dt > 1e-15 {
            (error - self.prev_error) / dt
        } else {
            0.0
        };
        self.prev_error = error;

        let output = self.kp * error + self.ki * self.integral + self.kd * derivative;

        if self.limits_set {
            output.max(self.min_output).min(self.max_output)
        } else {
            output
        }
    }

    /// Reset the controller state (integral and previous error).
    pub fn reset(&mut self) {
        self.integral = 0.0;
        self.prev_error = 0.0;
    }

    /// Set output saturation limits and return `&mut Self` for chaining.
    pub fn set_limits(&mut self, min: f64, max: f64) -> &mut Self {
        self.min_output = min;
        self.max_output = max;
        self.limits_set = true;
        self
    }

    /// Back-calculate anti-windup: if output is saturated, clamp the integral.
    pub fn anti_windup(&mut self, output: f64, min: f64, max: f64) {
        if output > max {
            self.integral -= (output - max) / (self.ki.abs() + 1e-15);
        } else if output < min {
            self.integral -= (output - min) / (self.ki.abs() + 1e-15);
        }
    }
}

// ── Impedance Controller ──────────────────────────────────────────────────────

/// Impedance controller with desired mass, damping, and stiffness.
#[derive(Debug, Clone)]
pub struct ImpedanceController {
    /// Desired virtual mass M_d.
    pub m_d: f64,
    /// Desired damping B_d.
    pub b_d: f64,
    /// Desired stiffness K_d.
    pub k_d: f64,
}

impl ImpedanceController {
    /// Construct a new `ImpedanceController`.
    pub fn new(m_d: f64, b_d: f64, k_d: f64) -> Self {
        Self { m_d, b_d, k_d }
    }

    /// Compute the desired interaction force from the impedance model.
    ///
    /// `F = M_d * acc + B_d * vel_error + K_d * pos_error`
    pub fn compute_force(&self, pos_error: f64, vel_error: f64, acc: f64) -> f64 {
        self.m_d * acc + self.b_d * vel_error + self.k_d * pos_error
    }

    /// Bandwidth (natural frequency) of the impedance filter.
    pub fn natural_frequency(&self) -> f64 {
        if self.m_d < 1e-15 {
            return 0.0;
        }
        (self.k_d / self.m_d).sqrt()
    }

    /// Damping ratio ζ.
    pub fn damping_ratio(&self) -> f64 {
        let omega_n = self.natural_frequency();
        if omega_n < 1e-15 || self.m_d < 1e-15 {
            return 0.0;
        }
        self.b_d / (2.0 * self.m_d * omega_n)
    }
}

// ── Adaptive servo (Ziegler–Nichols auto-tune) ────────────────────────────────

/// Number of relay half-cycles (period and peak amplitude) the auto-tuner
/// keeps; older measurements are discarded.
const MAX_RELAY_HISTORY: usize = 16;

/// Number of most recent half-cycles [`AdaptiveServo::auto_tune`] averages
/// (two full oscillation periods).
const TUNING_HALF_CYCLES: usize = 4;

/// Smallest oscillation amplitude from which an ultimate gain is estimated.
const MIN_OSCILLATION_AMPLITUDE: f64 = 1e-12;

/// Self-tuning PID controller using the Ziegler–Nichols relay method.
///
/// [`relay_step`](AdaptiveServo::relay_step) drives the plant with a relay
/// (`±relay_amplitude`, switching with the sign of the error) and records the
/// period and the peak error of every half cycle of the resulting oscillation.
/// [`auto_tune`](AdaptiveServo::auto_tune) turns the measurements into the
/// ultimate period `T_u` and, through the describing function of the relay,
/// the ultimate gain `K_u = 4·d / (π·a)` (`d` the relay amplitude, `a` the
/// amplitude of the oscillation of the error), and from them the classic
/// Ziegler–Nichols PID gains.
#[derive(Debug, Clone)]
pub struct AdaptiveServo {
    /// Current PID parameters.
    pub params: ServoParams,
    /// Estimated ultimate gain K_u.
    ultimate_gain: f64,
    /// Estimated ultimate period T_u (s).
    ultimate_period: f64,
    /// Durations of the most recent complete oscillation half-cycles.
    half_periods: Vec<f64>,
    /// Peak absolute error of each half cycle in `half_periods`.
    peaks: Vec<f64>,
    /// Peak absolute error seen so far in the current half cycle.
    half_cycle_peak: f64,
    /// Whether the relay has switched at least once; the half cycle before the
    /// first switch starts at `t = 0`, not at a switch, so it is not measured.
    switched: bool,
    /// Last relay output.
    relay_output: f64,
    /// Relay amplitude.
    relay_amplitude: f64,
    /// Last sign of the error.
    last_sign: f64,
    /// Time of the last sign change.
    last_switch_time: f64,
    /// Current simulation time accumulator.
    time: f64,
    /// Whether tuning is complete.
    tuned: bool,
}

impl AdaptiveServo {
    /// Construct a new `AdaptiveServo` with given initial params.
    pub fn new(params: ServoParams) -> Self {
        let relay_amplitude = params.max_force * 0.1;
        Self {
            params,
            ultimate_gain: 0.0,
            ultimate_period: 0.0,
            half_periods: Vec::new(),
            peaks: Vec::new(),
            half_cycle_peak: 0.0,
            switched: false,
            relay_output: relay_amplitude,
            relay_amplitude,
            last_sign: 1.0,
            last_switch_time: 0.0,
            time: 0.0,
            tuned: false,
        }
    }

    /// Feed one error sample at the current dt and return the relay output.
    ///
    /// Internally accumulates oscillation data for Ziegler–Nichols computation:
    /// at every sign change of the error the duration and the peak absolute
    /// error of the half cycle that just ended are recorded (the first half
    /// cycle, which starts at `t = 0` and not at a switch, is skipped).  Only
    /// the last 16 half cycles are kept.
    pub fn relay_step(&mut self, error: f64, dt: f64) -> f64 {
        self.time += dt;
        let sign = if error >= 0.0 { 1.0 } else { -1.0 };
        if sign != self.last_sign {
            if self.switched {
                self.record_half_cycle(self.time - self.last_switch_time, self.half_cycle_peak);
            }
            self.switched = true;
            self.half_cycle_peak = 0.0;
            self.last_switch_time = self.time;
            self.last_sign = sign;
            self.relay_output = sign * self.relay_amplitude;
        }
        // The sample that triggered a switch belongs to the new half cycle.
        self.half_cycle_peak = self.half_cycle_peak.max(error.abs());
        self.relay_output
    }

    /// Store one finished half cycle and drop the oldest ones beyond the
    /// history limit.
    fn record_half_cycle(&mut self, duration: f64, peak: f64) {
        self.half_periods.push(duration);
        self.peaks.push(peak);
        let excess = self.half_periods.len().saturating_sub(MAX_RELAY_HISTORY);
        if excess > 0 {
            self.half_periods.drain(..excess);
            self.peaks.drain(..excess);
        }
    }

    /// Perform Ziegler–Nichols PID auto-tune based on relay oscillation data.
    ///
    /// Requires at least 4 complete half-cycle measurements.  The ultimate
    /// period is twice the mean of the last four half periods, the ultimate
    /// gain is `4·d / (π·a)` with `d` the relay amplitude and `a` the mean peak
    /// error of the same four half cycles.  Returns `ServoParams` with the
    /// suggested PID gains (PI-D variant).
    ///
    /// When there is not enough data, or the measured oscillation is
    /// degenerate (amplitude below `1e-12`, no positive period, no relay
    /// amplitude), the current parameters are returned unchanged and
    /// [`is_tuned`](AdaptiveServo::is_tuned) is not set.
    pub fn auto_tune(&mut self) -> ServoParams {
        let measured = self.half_periods.len().min(self.peaks.len());
        if measured < TUNING_HALF_CYCLES {
            return self.params.clone();
        }

        let recent_mean = |values: &[f64]| {
            let tail = &values[values.len() - TUNING_HALF_CYCLES..];
            tail.iter().sum::<f64>() / TUNING_HALF_CYCLES as f64
        };
        // Average full period from the last half periods (2 full cycles).
        let tu = 2.0 * recent_mean(&self.half_periods);
        // Oscillation amplitude of the error over the same half cycles.
        let amplitude = recent_mean(&self.peaks);
        let usable = tu.is_finite()
            && tu > 0.0
            && amplitude.is_finite()
            && amplitude >= MIN_OSCILLATION_AMPLITUDE
            && self.relay_amplitude.is_finite()
            && self.relay_amplitude > 0.0;
        if !usable {
            return self.params.clone();
        }
        self.ultimate_period = tu;

        // Describing-function estimate of the ultimate gain: a relay of
        // amplitude d that produces an oscillation of amplitude a has the
        // equivalent gain 4·d / (π·a).
        let ku = 4.0 * self.relay_amplitude / (std::f64::consts::PI * amplitude);
        self.ultimate_gain = ku;

        // Ziegler–Nichols PID (classic):
        // Kp = 0.6 Ku, Ki = 2 Kp / Tu, Kd = Kp * Tu / 8
        let kp = 0.6 * ku;
        let ki = 2.0 * kp / tu;
        let kd = kp * tu / 8.0;

        self.tuned = true;
        let mut p = self.params.clone();
        p.kp = kp;
        p.ki = ki;
        p.kd = kd;
        p
    }

    /// Whether the auto-tune procedure has completed.
    pub fn is_tuned(&self) -> bool {
        self.tuned
    }

    /// Ultimate gain `K_u` estimated by the last successful
    /// [`auto_tune`](AdaptiveServo::auto_tune) (`0.0` before that).
    pub fn ultimate_gain(&self) -> f64 {
        self.ultimate_gain
    }

    /// Ultimate period `T_u` in seconds estimated by the last successful
    /// [`auto_tune`](AdaptiveServo::auto_tune) (`0.0` before that).
    pub fn ultimate_period(&self) -> f64 {
        self.ultimate_period
    }

    /// Number of complete relay half cycles currently held (at most 16).
    pub fn measured_half_cycles(&self) -> usize {
        self.half_periods.len()
    }
}

// ── Feedforward compensation ──────────────────────────────────────────────────

/// Gravity and friction feedforward compensator.
pub struct FeedforwardCompensation;

impl FeedforwardCompensation {
    /// Gravity compensation torque for a link of `mass` kg at angle `angle` (rad)
    /// with gravitational acceleration `g` (m/s²).
    ///
    /// `T_grav = mass * g * cos(angle)`
    pub fn gravity_comp(mass: f64, g: f64, angle: f64) -> f64 {
        mass * g * angle.cos()
    }

    /// Friction compensation force combining Coulomb and viscous friction.
    ///
    /// `F_fric = coulomb * sign(vel) + viscous * vel`
    pub fn friction_comp(vel: f64, coulomb: f64, viscous: f64) -> f64 {
        if vel.abs() < 1e-9 {
            // Static case: no motion
            0.0
        } else {
            coulomb * vel.signum() + viscous * vel
        }
    }

    /// Combined feedforward force.
    pub fn combined(mass: f64, g: f64, angle: f64, vel: f64, coulomb: f64, viscous: f64) -> f64 {
        Self::gravity_comp(mass, g, angle) + Self::friction_comp(vel, coulomb, viscous)
    }
}

// ── Cascade controller ────────────────────────────────────────────────────────

/// Two-loop cascade controller: outer position loop driving inner velocity loop.
#[derive(Debug, Clone)]
pub struct CascadeController {
    /// Outer (position) loop PID.
    pub outer_pid: PidController,
    /// Inner (velocity) loop PID.
    pub inner_pid: PidController,
    /// Velocity command limit (m/s or rad/s).
    pub vel_limit: f64,
}

impl CascadeController {
    /// Construct a `CascadeController`.
    pub fn new(outer: PidController, inner: PidController, vel_limit: f64) -> Self {
        Self {
            outer_pid: outer,
            inner_pid: inner,
            vel_limit,
        }
    }

    /// Compute the force/torque output.
    ///
    /// 1. Outer PID generates a velocity command from position error.
    /// 2. Inner PID generates force from velocity error.
    pub fn step(&mut self, pos_target: f64, pos_current: f64, vel_current: f64, dt: f64) -> f64 {
        let pos_error = pos_target - pos_current;
        let vel_cmd = self
            .outer_pid
            .step(pos_error, dt)
            .max(-self.vel_limit)
            .min(self.vel_limit);

        let vel_error = vel_cmd - vel_current;
        self.inner_pid.step(vel_error, dt)
    }

    /// Reset both loops.
    pub fn reset(&mut self) {
        self.outer_pid.reset();
        self.inner_pid.reset();
    }
}

// ── Servo actuator ────────────────────────────────────────────────────────────

/// A complete single-axis servo with position limits and soft stops.
#[derive(Debug, Clone)]
pub struct ServoActuator {
    /// Name / label for debugging.
    pub name: String,
    /// Current joint position.
    pub position: f64,
    /// Current joint velocity.
    pub velocity: f64,
    /// Lower hard stop (rad or m).
    pub pos_min: f64,
    /// Upper hard stop (rad or m).
    pub pos_max: f64,
    /// Soft stop margin (additional boundary).
    pub soft_margin: f64,
    /// Actuator mass / inertia (kg).
    pub mass: f64,
    /// Internal servo constraint.
    servo: ServoConstraint,
    /// Whether homing is complete.
    pub homed: bool,
    /// Accumulated homing travel.
    homing_travel: f64,
}

impl ServoActuator {
    /// Create a new `ServoActuator`.
    pub fn new(
        name: impl Into<String>,
        mass: f64,
        pos_min: f64,
        pos_max: f64,
        soft_margin: f64,
        params: ServoParams,
    ) -> Self {
        let axis = [1.0, 0.0, 0.0];
        let servo = ServoConstraint::new(0, 1, axis, ServoMode::Position, params);
        Self {
            name: name.into(),
            position: 0.0,
            velocity: 0.0,
            pos_min,
            pos_max,
            soft_margin,
            mass,
            servo,
            homed: false,
            homing_travel: 0.0,
        }
    }

    /// Set a target position.
    pub fn set_target(&mut self, target: f64) {
        let clamped = target
            .max(self.pos_min + self.soft_margin)
            .min(self.pos_max - self.soft_margin);
        self.servo.set_target(clamped);
    }

    /// Advance the actuator by one time step.
    ///
    /// * `dt` — time step (s)
    /// * `external_load` — external force / torque (N or N·m)
    pub fn step(&mut self, dt: f64, external_load: f64) {
        let impulse = self.servo.compute_impulse(dt, self.position, self.velocity);
        let force = impulse / dt + external_load;

        // Soft stops: exponential repulsion near limits
        let force = force + self.soft_stop_force();

        // Integrate velocity and position (forward Euler)
        let acc = force / (self.mass.max(1e-15));
        self.velocity += acc * dt;

        // Velocity limit
        let vmax = self.servo.params.max_velocity;
        self.velocity = self.velocity.max(-vmax).min(vmax);

        self.position += self.velocity * dt;

        // Hard stop clamp
        self.position = self.position.max(self.pos_min).min(self.pos_max);
    }

    /// Compute soft-stop repulsion force.
    fn soft_stop_force(&self) -> f64 {
        let k_soft = self.servo.params.stiffness.max(1.0);
        let mut f = 0.0;
        let lo = self.pos_min + self.soft_margin;
        let hi = self.pos_max - self.soft_margin;
        if self.position < lo {
            f += k_soft * (lo - self.position);
        } else if self.position > hi {
            f -= k_soft * (self.position - hi);
        }
        f
    }

    /// Execute a homing sequence: move toward `pos_min` until the limit is reached.
    ///
    /// Returns `true` when homing is complete.
    pub fn home_sequence(&mut self, dt: f64) -> bool {
        if self.homed {
            return true;
        }
        self.servo.set_target(self.pos_min);
        self.step(dt, 0.0);
        self.homing_travel += self.velocity.abs() * dt;

        // Done when close to the lower limit
        if (self.position - self.pos_min).abs() < 1e-4 {
            self.position = self.pos_min;
            self.velocity = 0.0;
            self.homed = true;
            return true;
        }
        false
    }

    /// Check if the actuator is at its target position.
    pub fn is_at_target(&self, tol: f64) -> bool {
        self.servo.is_at_target(tol)
    }
}

// ── Multi-axis servo ──────────────────────────────────────────────────────────

/// Coordinated multi-axis servo system.
pub struct MultiAxisServo {
    /// All servo axes.
    pub axes: Vec<ServoActuator>,
}

impl MultiAxisServo {
    /// Create a `MultiAxisServo` from a vector of actuators.
    pub fn new(axes: Vec<ServoActuator>) -> Self {
        Self { axes }
    }

    /// Linearly interpolate between `start` and `end` at parameter `t ∈ [0, 1]`.
    ///
    /// Panics if `start.len() != end.len()`.
    pub fn linear_interpolation(start: &[f64], end: &[f64], t: f64) -> Vec<f64> {
        assert_eq!(
            start.len(),
            end.len(),
            "start and end must have equal length"
        );
        let t = t.clamp(0.0, 1.0);
        start
            .iter()
            .zip(end.iter())
            .map(|(a, b)| a + (b - a) * t)
            .collect()
    }

    /// Compute the interpolated target positions for a coordinated move
    /// of `speed` (fraction of full travel per second) at current time `t`.
    ///
    /// Returns the position vector at parameter `t * speed`.
    pub fn coordinated_move(targets: &[f64], speed: f64) -> Vec<f64> {
        // Normalize speed to [0, 1]
        let t = speed.clamp(0.0, 1.0);
        // Return linearly interpolated positions from zero
        targets.iter().map(|&v| v * t).collect()
    }

    /// Step all axes simultaneously.
    pub fn step_all(&mut self, dt: f64, loads: &[f64]) {
        for (axis, load) in self.axes.iter_mut().zip(loads.iter()) {
            axis.step(dt, *load);
        }
    }

    /// Set targets for all axes.
    ///
    /// `targets` must have the same length as `self.axes`.
    pub fn set_targets(&mut self, targets: &[f64]) {
        for (axis, &t) in self.axes.iter_mut().zip(targets.iter()) {
            axis.set_target(t);
        }
    }

    /// Check if all axes are at their targets.
    pub fn all_at_target(&self, tol: f64) -> bool {
        self.axes.iter().all(|a| a.is_at_target(tol))
    }

    /// Number of axes.
    pub fn num_axes(&self) -> usize {
        self.axes.len()
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn default_params() -> ServoParams {
        ServoParams::new(100.0, 10.0, 5.0, 500.0, 5.0)
    }

    fn make_actuator(pos_min: f64, pos_max: f64) -> ServoActuator {
        ServoActuator::new("test", 1.0, pos_min, pos_max, 0.01, default_params())
    }

    // 1. PID step returns proportional term only when ki=kd=0.
    #[test]
    fn test_pid_proportional_only() {
        let mut pid = PidController::new(10.0, 0.0, 0.0);
        let output = pid.step(2.0, 0.01);
        assert!((output - 20.0).abs() < 1e-10, "P-only output={output}");
    }

    // 2. PID integral accumulates over multiple steps.
    #[test]
    fn test_pid_integral_accumulates() {
        let mut pid = PidController::new(0.0, 1.0, 0.0);
        let dt = 0.1;
        for _ in 0..10 {
            pid.step(1.0, dt);
        }
        // integral = 10 * 1.0 * 0.1 = 1.0; output = ki * integral = 1.0
        let output = pid.step(0.0, dt);
        assert!((output - 1.0).abs() < 1e-9, "integral output={output}");
    }

    // 3. PID derivative term works on error change.
    #[test]
    fn test_pid_derivative() {
        let mut pid = PidController::new(0.0, 0.0, 1.0);
        let dt = 0.1;
        pid.step(0.0, dt); // prev_error = 0
        let output = pid.step(1.0, dt); // de = 1.0 / dt = 10
        assert!((output - 10.0).abs() < 1e-9, "D output={output}");
    }

    // 4. PID reset clears integral and previous error.
    #[test]
    fn test_pid_reset() {
        let mut pid = PidController::new(1.0, 1.0, 1.0);
        pid.step(5.0, 0.1);
        pid.reset();
        let output = pid.step(0.0, 0.1);
        assert!(output.abs() < 1e-10, "output after reset should be 0");
    }

    // 5. PID set_limits clamps output.
    #[test]
    fn test_pid_limits() {
        let mut pid = PidController::new(100.0, 0.0, 0.0);
        pid.set_limits(-10.0, 10.0);
        let output = pid.step(5.0, 0.01); // would be 500 without limits
        assert!((output - 10.0).abs() < 1e-10, "clamped output={output}");
    }

    // 6. PID set_limits lower clamp.
    #[test]
    fn test_pid_limits_lower() {
        let mut pid = PidController::new(100.0, 0.0, 0.0);
        pid.set_limits(-10.0, 10.0);
        let output = pid.step(-5.0, 0.01);
        assert!((output + 10.0).abs() < 1e-10, "lower clamp output={output}");
    }

    // 7. ServoConstraint starts with zero output.
    #[test]
    fn test_servo_constraint_initial_zero() {
        let sc = ServoConstraint::new(0, 1, [1.0, 0.0, 0.0], ServoMode::Position, default_params());
        assert_eq!(sc.state.output_force, 0.0);
    }

    // 8. ServoConstraint compute_impulse is non-zero when there is an error.
    #[test]
    fn test_servo_impulse_nonzero_error() {
        let mut sc =
            ServoConstraint::new(0, 1, [1.0, 0.0, 0.0], ServoMode::Position, default_params());
        sc.set_target(1.0);
        let impulse = sc.compute_impulse(0.01, 0.0, 0.0);
        assert!(impulse.abs() > 0.0, "impulse should be nonzero");
    }

    // 9. is_at_target returns true when error ≤ tol.
    #[test]
    fn test_servo_at_target() {
        let mut sc =
            ServoConstraint::new(0, 1, [1.0, 0.0, 0.0], ServoMode::Position, default_params());
        sc.set_target(1.0);
        sc.compute_impulse(0.01, 1.0, 0.0); // zero error
        assert!(sc.is_at_target(1e-6), "should be at target");
    }

    // 10. is_at_target returns false when error > tol.
    #[test]
    fn test_servo_not_at_target() {
        let mut sc =
            ServoConstraint::new(0, 1, [1.0, 0.0, 0.0], ServoMode::Position, default_params());
        sc.set_target(1.0);
        sc.compute_impulse(0.01, 0.5, 0.0); // error = 0.5
        assert!(!sc.is_at_target(0.1), "should not be at target");
    }

    // 11. Impulse is clamped by max_force.
    #[test]
    fn test_servo_impulse_clamped() {
        let params = ServoParams::new(1e9, 0.0, 0.0, 100.0, 100.0);
        let mut sc = ServoConstraint::new(0, 1, [1.0, 0.0, 0.0], ServoMode::Position, params);
        sc.set_target(1.0e6); // huge target → huge P error
        let impulse = sc.compute_impulse(0.01, 0.0, 0.0);
        assert!(
            impulse.abs() <= 100.0 * 0.01 + 1e-9,
            "impulse should be clamped: {impulse}"
        );
    }

    // 12. Torque mode passes target directly.
    #[test]
    fn test_torque_mode() {
        let mut sc =
            ServoConstraint::new(0, 1, [1.0, 0.0, 0.0], ServoMode::Torque, default_params());
        sc.set_target(50.0);
        let impulse = sc.compute_impulse(0.01, 0.0, 0.0);
        assert!(
            (impulse - 50.0 * 0.01).abs() < 1e-9,
            "torque impulse={impulse}"
        );
    }

    // 13. ImpedanceController compute_force at rest.
    #[test]
    fn test_impedance_zero_force_at_rest() {
        let ic = ImpedanceController::new(1.0, 10.0, 100.0);
        let f = ic.compute_force(0.0, 0.0, 0.0);
        assert!(f.abs() < 1e-12, "force at rest should be 0");
    }

    // 14. ImpedanceController spring term.
    #[test]
    fn test_impedance_spring_term() {
        let ic = ImpedanceController::new(0.0, 0.0, 100.0);
        let f = ic.compute_force(0.1, 0.0, 0.0); // pos_error = 0.1
        assert!((f - 10.0).abs() < 1e-10, "spring force={f}");
    }

    // 15. ImpedanceController natural frequency.
    #[test]
    fn test_impedance_natural_frequency() {
        let ic = ImpedanceController::new(1.0, 0.0, 100.0);
        let wn = ic.natural_frequency();
        assert!((wn - 10.0).abs() < 1e-9, "omega_n={wn}");
    }

    // 16. ImpedanceController damping ratio.
    #[test]
    fn test_impedance_damping_ratio() {
        let ic = ImpedanceController::new(1.0, 20.0, 100.0);
        let zeta = ic.damping_ratio();
        // B / (2 * m * omega_n) = 20 / (2 * 1 * 10) = 1.0
        assert!((zeta - 1.0).abs() < 1e-9, "damping_ratio={zeta}");
    }

    // 17. FeedforwardCompensation gravity at 0° is max.
    #[test]
    fn test_gravity_comp_zero_angle() {
        let f = FeedforwardCompensation::gravity_comp(10.0, 9.81, 0.0);
        assert!((f - 10.0 * 9.81).abs() < 1e-10, "gravity at 0° f={f}");
    }

    // 18. FeedforwardCompensation gravity at 90° is zero.
    #[test]
    fn test_gravity_comp_ninety_deg() {
        let f = FeedforwardCompensation::gravity_comp(10.0, 9.81, std::f64::consts::FRAC_PI_2);
        assert!(f.abs() < 1e-9, "gravity at 90° f={f}");
    }

    // 19. FeedforwardCompensation friction Coulomb sign.
    #[test]
    fn test_friction_coulomb_sign() {
        let f_pos = FeedforwardCompensation::friction_comp(1.0, 5.0, 0.0);
        let f_neg = FeedforwardCompensation::friction_comp(-1.0, 5.0, 0.0);
        assert!(f_pos > 0.0, "positive velocity → positive Coulomb");
        assert!(f_neg < 0.0, "negative velocity → negative Coulomb");
    }

    // 20. FeedforwardCompensation friction zero at rest.
    #[test]
    fn test_friction_zero_at_rest() {
        let f = FeedforwardCompensation::friction_comp(0.0, 5.0, 2.0);
        assert!(f.abs() < 1e-12, "no friction at zero velocity");
    }

    // 21. CascadeController step drives position toward target over multiple steps.
    #[test]
    fn test_cascade_drives_to_target() {
        let outer = PidController::new(10.0, 0.0, 0.0);
        let inner = PidController::new(5.0, 0.0, 0.0);
        let mut cas = CascadeController::new(outer, inner, 2.0);

        let mut pos = 0.0_f64;
        let mut vel = 0.0_f64;
        let dt = 0.01;
        for _ in 0..500 {
            let force = cas.step(1.0, pos, vel, dt);
            vel += force * dt;
            pos += vel * dt;
        }
        assert!(
            (pos - 1.0).abs() < 0.5,
            "cascade should drive toward target, pos={pos:.4}"
        );
    }

    // 22. CascadeController reset clears both loops.
    #[test]
    fn test_cascade_reset() {
        let outer = PidController::new(10.0, 1.0, 0.0);
        let inner = PidController::new(5.0, 1.0, 0.0);
        let mut cas = CascadeController::new(outer, inner, 10.0);
        cas.step(1.0, 0.0, 0.0, 0.01);
        cas.reset();
        let out = cas.step(0.0, 0.0, 0.0, 0.01);
        assert!(out.abs() < 1e-9, "after reset output should be 0");
    }

    // 23. ServoActuator clamps position to hard stops.
    #[test]
    fn test_actuator_hard_stop() {
        let mut act = make_actuator(-1.0, 1.0);
        act.position = 1.5; // beyond limit
        act.step(0.01, 0.0);
        assert!(act.position <= 1.0, "position should be clamped to max");
    }

    // 24. ServoActuator soft_margin clamps target.
    #[test]
    fn test_actuator_target_clamped_by_soft_margin() {
        let mut act = make_actuator(-1.0, 1.0);
        act.set_target(0.995); // within soft margin of 0.01 from 1.0
        // Target should be clamped to 0.99
        assert!(act.servo.state.target <= 1.0 - 0.01 + 1e-9);
    }

    // 25. ServoActuator step moves toward target.
    #[test]
    fn test_actuator_moves_toward_target() {
        let mut act = make_actuator(-5.0, 5.0);
        act.set_target(1.0);
        let pos0 = act.position;
        for _ in 0..100 {
            act.step(0.01, 0.0);
        }
        let pos1 = act.position;
        assert!(
            (pos1 - 1.0).abs() < (pos0 - 1.0).abs() + 1e-9,
            "actuator should move toward target: pos={pos1:.4}"
        );
    }

    // 26. ServoActuator home_sequence returns true eventually.
    #[test]
    fn test_home_sequence_completes() {
        let params = ServoParams::new(1000.0, 100.0, 50.0, 2000.0, 20.0);
        let mut act = ServoActuator::new("home_test", 0.1, 0.0, 2.0, 0.01, params);
        act.position = 1.0;
        let mut homed = false;
        for _ in 0..5000 {
            if act.home_sequence(0.01) {
                homed = true;
                break;
            }
        }
        assert!(homed, "home sequence should complete within 5000 steps");
    }

    // 27. MultiAxisServo linear_interpolation at t=0 returns start.
    #[test]
    fn test_linear_interp_t0() {
        let start = vec![0.0, 1.0, 2.0];
        let end = vec![3.0, 4.0, 5.0];
        let result = MultiAxisServo::linear_interpolation(&start, &end, 0.0);
        for (r, s) in result.iter().zip(start.iter()) {
            assert!((r - s).abs() < 1e-12);
        }
    }

    // 28. MultiAxisServo linear_interpolation at t=1 returns end.
    #[test]
    fn test_linear_interp_t1() {
        let start = vec![0.0, 1.0, 2.0];
        let end = vec![3.0, 4.0, 5.0];
        let result = MultiAxisServo::linear_interpolation(&start, &end, 1.0);
        for (r, e) in result.iter().zip(end.iter()) {
            assert!((r - e).abs() < 1e-12);
        }
    }

    // 29. MultiAxisServo linear_interpolation at t=0.5 is midpoint.
    #[test]
    fn test_linear_interp_midpoint() {
        let start = vec![0.0, 0.0];
        let end = vec![2.0, 4.0];
        let result = MultiAxisServo::linear_interpolation(&start, &end, 0.5);
        assert!((result[0] - 1.0).abs() < 1e-12);
        assert!((result[1] - 2.0).abs() < 1e-12);
    }

    // 30. MultiAxisServo coordinated_move scales target by speed.
    #[test]
    fn test_coordinated_move_half_speed() {
        let targets = vec![2.0, 4.0];
        let result = MultiAxisServo::coordinated_move(&targets, 0.5);
        assert!((result[0] - 1.0).abs() < 1e-12);
        assert!((result[1] - 2.0).abs() < 1e-12);
    }

    // 31. MultiAxisServo coordinated_move at speed=1 returns targets.
    #[test]
    fn test_coordinated_move_full_speed() {
        let targets = vec![3.0, 6.0];
        let result = MultiAxisServo::coordinated_move(&targets, 1.0);
        for (r, t) in result.iter().zip(targets.iter()) {
            assert!((r - t).abs() < 1e-12);
        }
    }

    // 32. MultiAxisServo all_at_target initial check.
    #[test]
    fn test_multi_axis_all_at_target_initial() {
        let axes = vec![make_actuator(-1.0, 1.0), make_actuator(-1.0, 1.0)];
        let multi = MultiAxisServo::new(axes);
        // Targets default to 0; positions are 0 → should be at target
        assert!(multi.all_at_target(0.01));
    }

    // 33. MultiAxisServo num_axes.
    #[test]
    fn test_multi_axis_num_axes() {
        let axes = vec![
            make_actuator(-1.0, 1.0),
            make_actuator(-2.0, 2.0),
            make_actuator(-3.0, 3.0),
        ];
        let multi = MultiAxisServo::new(axes);
        assert_eq!(multi.num_axes(), 3);
    }

    // 34. AdaptiveServo relay_step flips sign.
    #[test]
    fn test_adaptive_servo_relay_flip() {
        let mut adap = AdaptiveServo::new(default_params());
        let out1 = adap.relay_step(1.0, 0.01); // positive error
        let out2 = adap.relay_step(-1.0, 0.01); // negative error → flip
        assert!(
            out1 > 0.0 && out2 < 0.0,
            "relay should flip sign: {out1}, {out2}"
        );
    }

    // 35. AdaptiveServo auto_tune returns valid params after sufficient oscillation.
    #[test]
    fn test_adaptive_servo_auto_tune() {
        let mut adap = AdaptiveServo::new(default_params());
        // Simulate several oscillation half-cycles by alternating sign every step
        let dt = 0.05;
        for i in 0..20 {
            let err = if i % 2 == 0 { 1.0 } else { -1.0 };
            adap.relay_step(err, dt);
        }
        let p = adap.auto_tune();
        assert!(
            p.kp >= 0.0 && p.ki >= 0.0 && p.kd >= 0.0,
            "tuned gains must be non-negative"
        );
        assert!(adap.is_tuned());
        // Every sample is its own half cycle: |error| = 1, half period 0.05 s,
        // relay amplitude 0.1 * 500 N.
        assert!((adap.ultimate_period() - 0.1).abs() < 1e-12);
        let expected_ku = 4.0 * 50.0 / std::f64::consts::PI;
        assert!((adap.ultimate_gain() - expected_ku).abs() < 1e-9);
        assert!((p.kp - 0.6 * expected_ku).abs() < 1e-9);
    }

    // ── F0468: ultimate gain from the measured oscillation amplitude ──────────

    /// Run a relay experiment on the plant `1 / (s + 1)^3` (three cascaded
    /// first-order lags, Euler-integrated) with the given relay amplitude.
    ///
    /// The plant has the exact ultimate gain `K_u = 8` and ultimate period
    /// `T_u = 2π / √3`.
    fn relay_experiment_on_third_order_lag(relay_amplitude: f64) -> AdaptiveServo {
        // The relay amplitude of `AdaptiveServo` is a tenth of `max_force`.
        let mut servo =
            AdaptiveServo::new(ServoParams::new(1.0, 0.0, 0.0, relay_amplitude * 10.0, 1.0));
        let dt = 1e-3;
        let (mut x1, mut x2, mut x3) = (0.0_f64, 0.0_f64, 0.0_f64);
        for _ in 0..60_000 {
            // Set-point 0: the error is the negated plant output.
            let u = servo.relay_step(-x3, dt);
            x1 += dt * (u - x1);
            x2 += dt * (x1 - x2);
            x3 += dt * (x2 - x3);
        }
        servo
    }

    #[test]
    fn test_adaptive_servo_estimates_ultimate_gain_of_known_plant() {
        let true_ku = 8.0;
        let true_tu = 2.0 * std::f64::consts::PI / 3.0_f64.sqrt();
        // For a linear plant the estimate must not depend on the relay size.
        for relay_amplitude in [0.5, 1.0, 4.0] {
            let mut servo = relay_experiment_on_third_order_lag(relay_amplitude);
            let tuned = servo.auto_tune();
            assert!(servo.is_tuned(), "relay amplitude {relay_amplitude}");
            let ku = servo.ultimate_gain();
            let tu = servo.ultimate_period();
            assert!(
                (ku - true_ku).abs() / true_ku < 0.10,
                "K_u = {ku}, expected {true_ku} (relay amplitude {relay_amplitude})"
            );
            assert!(
                (tu - true_tu).abs() / true_tu < 0.10,
                "T_u = {tu}, expected {true_tu} (relay amplitude {relay_amplitude})"
            );
            // Classic Ziegler-Nichols PID from the estimate.
            assert!((tuned.kp - 0.6 * ku).abs() < 1e-12);
            assert!((tuned.ki - 2.0 * tuned.kp / tu).abs() < 1e-9);
            assert!((tuned.kd - tuned.kp * tu / 8.0).abs() < 1e-9);
        }
    }

    /// A servo with the default relay (amplitude 50) tuned on a sine error of
    /// the given amplitude and a 2 s period.
    fn servo_tuned_on_sine_error(amplitude: f64) -> AdaptiveServo {
        let mut servo = AdaptiveServo::new(default_params());
        let dt = 0.01;
        for k in 0..1000 {
            // Sample between the grid points so no sample is exactly zero.
            let t = (k as f64 + 0.5) * dt;
            servo.relay_step(amplitude * (std::f64::consts::PI * t).sin(), dt);
        }
        servo.auto_tune();
        servo
    }

    #[test]
    fn test_adaptive_servo_gain_is_inversely_proportional_to_amplitude() {
        let small = servo_tuned_on_sine_error(0.25);
        let large = servo_tuned_on_sine_error(1.0);
        assert!(small.is_tuned() && large.is_tuned());
        let ratio = small.ultimate_gain() / large.ultimate_gain();
        assert!(
            (ratio - 4.0).abs() < 0.05,
            "a 4x smaller oscillation means a 4x larger K_u, got ratio {ratio}"
        );
        // K_u = 4 d / (pi a) with d = 50 and a = 1 (the sampled peak is
        // within 0.02 % of the true amplitude).
        let expected = 4.0 * 50.0 / std::f64::consts::PI;
        assert!((large.ultimate_gain() - expected).abs() / expected < 0.01);
        assert!((large.ultimate_period() - 2.0).abs() < 1e-9);
    }

    /// Feed `samples` errors of constant magnitude `amplitude` and sign `sign`.
    fn feed_half_cycle(servo: &mut AdaptiveServo, sign: f64, amplitude: f64, samples: usize) {
        for _ in 0..samples {
            servo.relay_step(sign * amplitude, 0.01);
        }
    }

    #[test]
    fn test_adaptive_servo_skips_the_half_cycle_measured_from_time_zero() {
        let mut servo = AdaptiveServo::new(default_params());
        // 3 s of large positive error before the first switch: this half cycle
        // starts at t = 0, not at a relay switch, so it must not be measured.
        feed_half_cycle(&mut servo, 1.0, 5.0, 300);
        assert_eq!(servo.measured_half_cycles(), 0);
        // Four more half cycles (1 s each, amplitude 0.5): the first switch
        // ends the unmeasured half cycle, so four switches record three.
        for half in 0..4 {
            let sign = if half % 2 == 0 { -1.0 } else { 1.0 };
            feed_half_cycle(&mut servo, sign, 0.5, 100);
        }
        assert_eq!(servo.measured_half_cycles(), 3);
        let before = servo.params.clone();
        let unchanged = servo.auto_tune();
        assert_eq!(unchanged.kp, before.kp);
        assert_eq!(unchanged.ki, before.ki);
        assert_eq!(unchanged.kd, before.kd);
        assert!(!servo.is_tuned(), "three half cycles are not enough");

        // The fifth switch completes the fourth measured half cycle.
        feed_half_cycle(&mut servo, -1.0, 0.5, 100);
        assert_eq!(servo.measured_half_cycles(), 4);
        servo.auto_tune();
        assert!(servo.is_tuned());
        assert!((servo.ultimate_period() - 2.0).abs() < 1e-9);
        let expected_ku = 4.0 * 50.0 / (std::f64::consts::PI * 0.5);
        assert!(
            (servo.ultimate_gain() - expected_ku).abs() < 1e-9,
            "K_u = {}, the transient 5.0 peak must not leak in",
            servo.ultimate_gain()
        );
    }

    #[test]
    fn test_adaptive_servo_history_is_bounded() {
        let mut servo = AdaptiveServo::new(default_params());
        for half in 0..200 {
            let sign = if half % 2 == 0 { 1.0 } else { -1.0 };
            // Early half cycles are long and large, the last ones short.
            let (amplitude, samples) = if half < 150 { (4.0, 50) } else { (1.0, 20) };
            feed_half_cycle(&mut servo, sign, amplitude, samples);
        }
        assert_eq!(servo.measured_half_cycles(), 16);
        servo.auto_tune();
        // Only the most recent half cycles (0.2 s, amplitude 1) are used.
        assert!((servo.ultimate_period() - 0.4).abs() < 1e-9);
        let expected_ku = 4.0 * 50.0 / std::f64::consts::PI;
        assert!((servo.ultimate_gain() - expected_ku).abs() < 1e-9);
    }

    #[test]
    fn test_adaptive_servo_negligible_oscillation_keeps_params() {
        let mut servo = AdaptiveServo::new(default_params());
        for half in 0..12 {
            let sign = if half % 2 == 0 { 1.0 } else { -1.0 };
            feed_half_cycle(&mut servo, sign, 1e-15, 10);
        }
        assert!(servo.measured_half_cycles() >= 4);
        let before = servo.params.clone();
        let tuned = servo.auto_tune();
        assert_eq!(
            (tuned.kp, tuned.ki, tuned.kd),
            (before.kp, before.ki, before.kd)
        );
        assert!(!servo.is_tuned(), "an unmeasurable oscillation cannot tune");
        assert_eq!(servo.ultimate_gain(), 0.0);
    }

    #[test]
    fn test_adaptive_servo_auto_tune_needs_four_half_cycles() {
        let mut servo = AdaptiveServo::new(default_params());
        let before = servo.params.clone();
        let tuned = servo.auto_tune();
        assert_eq!(
            (tuned.kp, tuned.ki, tuned.kd),
            (before.kp, before.ki, before.kd)
        );
        assert!(!servo.is_tuned());
    }
}
