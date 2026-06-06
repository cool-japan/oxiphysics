//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{discretize_zoh, mat2_mul, mat2_mul_transpose};

/// Smith predictor wrapper for a PID controller with process dead-time.
///
/// Compensates for dead-time L by subtracting the predicted delayed response
/// from the feedback signal.
pub struct SmithPredictor {
    /// Inner PID controller.
    pub pid: BumplessPid,
    /// Dead-time in number of discrete steps.
    pub delay_steps: usize,
    /// Model output buffer for dead-time simulation.
    pub(super) model_buffer: Vec<f64>,
    /// Model gain (approximate plant gain at DC).
    pub model_gain: f64,
    /// Model time constant τ.
    pub model_tau: f64,
    /// Model state (first-order approximation).
    pub(super) model_state: f64,
    /// Time step dt.
    pub dt: f64,
}
impl SmithPredictor {
    /// Create a new Smith predictor.
    ///
    /// * `pid` - Inner PID controller.
    /// * `delay_steps` - Number of steps of dead-time.
    /// * `model_gain`, `model_tau` - First-order plant model parameters.
    /// * `dt` - Simulation time step.
    pub fn new(
        pid: BumplessPid,
        delay_steps: usize,
        model_gain: f64,
        model_tau: f64,
        dt: f64,
    ) -> Self {
        let model_buffer = vec![0.0; delay_steps + 1];
        Self {
            pid,
            delay_steps,
            model_buffer,
            model_gain,
            model_tau,
            model_state: 0.0,
            dt,
        }
    }
    /// Compute one step of the Smith predictor.
    ///
    /// * `setpoint` - Desired output.
    /// * `actual_output` - Actual (delayed) plant output.
    ///
    /// Returns the control signal.
    pub fn update(&mut self, setpoint: f64, actual_output: f64) -> f64 {
        let u_prev = *self.model_buffer.last().unwrap_or(&0.0);
        let dx = (-self.model_state + self.model_gain * u_prev) / self.model_tau;
        self.model_state += dx * self.dt;
        let buf_len = self.model_buffer.len();
        if buf_len > 1 {
            for i in (1..buf_len).rev() {
                self.model_buffer[i] = self.model_buffer[i - 1];
            }
        }
        let delayed_model = if buf_len > 0 {
            self.model_buffer[buf_len.saturating_sub(1)]
        } else {
            0.0
        };
        let corrected_feedback = actual_output + self.model_state - delayed_model;
        let error = setpoint - corrected_feedback;
        let u = self.pid.update(error, self.dt);
        if !self.model_buffer.is_empty() {
            self.model_buffer[0] = u;
        }
        u
    }
}
/// A classic PID controller with anti-windup and output saturation.
pub struct PidController {
    /// Proportional gain.
    pub kp: f64,
    /// Integral gain.
    pub ki: f64,
    /// Derivative gain.
    pub kd: f64,
    /// Accumulated integral error.
    pub integral: f64,
    /// Error from previous update step.
    pub prev_error: f64,
    /// Anti-windup clamp on integral term.
    pub integral_limit: f64,
    /// Output saturation limit.
    pub output_limit: f64,
}
impl PidController {
    /// Create a new PID controller with the given gains.
    ///
    /// Defaults: `integral_limit = f64::MAX`, `output_limit = f64::MAX`.
    pub fn new(kp: f64, ki: f64, kd: f64) -> Self {
        Self {
            kp,
            ki,
            kd,
            integral: 0.0,
            prev_error: 0.0,
            integral_limit: f64::MAX,
            output_limit: f64::MAX,
        }
    }
    /// Compute the PID output for the given error and time-step.
    pub fn update(&mut self, error: f64, dt: f64) -> f64 {
        self.integral += error * dt;
        self.integral = self
            .integral
            .clamp(-self.integral_limit, self.integral_limit);
        let derivative = if dt > 0.0 {
            (error - self.prev_error) / dt
        } else {
            0.0
        };
        self.prev_error = error;
        let output = self.kp * error + self.ki * self.integral + self.kd * derivative;
        output.clamp(-self.output_limit, self.output_limit)
    }
    /// Reset accumulated state (integral and previous error).
    pub fn reset(&mut self) {
        self.integral = 0.0;
        self.prev_error = 0.0;
    }
    /// Set the anti-windup and output saturation limits.
    pub fn set_limits(&mut self, integral_limit: f64, output_limit: f64) {
        self.integral_limit = integral_limit;
        self.output_limit = output_limit;
    }
}
/// Relay feedback experiment result for Ziegler-Nichols auto-tuning.
///
/// In a relay feedback test a relay of amplitude `d` forces the plant into
/// sustained oscillation.  The critical gain and period are estimated from
/// the oscillation amplitude and period.
pub struct RelayFeedbackResult {
    /// Estimated critical gain Ku.
    pub ku: f64,
    /// Estimated critical period Tu (seconds).
    pub tu: f64,
}
/// Linear spring-damper for soft joint limits or elastic connections.
pub struct SpringDamper {
    /// Spring stiffness coefficient (N/m or N·m/rad).
    pub stiffness: f64,
    /// Damping coefficient (N·s/m or N·m·s/rad).
    pub damping: f64,
    /// Natural rest length / rest angle.
    pub rest_length: f64,
}
impl SpringDamper {
    /// Compute the restoring force/torque.
    ///
    /// `F = -k * (displacement - rest_length) - c * velocity`
    pub fn force(&self, displacement: f64, velocity: f64) -> f64 {
        -self.stiffness * (displacement - self.rest_length) - self.damping * velocity
    }
    /// Total mechanical energy stored (potential + kinetic proxy via damping is
    /// dissipative, so only spring potential is returned).
    ///
    /// `E = 0.5 * k * x^2 + 0.5 * c * v^2`  (where `x = displacement - rest_length`)
    pub fn energy(&self, displacement: f64, velocity: f64) -> f64 {
        let x = displacement - self.rest_length;
        0.5 * self.stiffness * x * x + 0.5 * self.damping * velocity * velocity
    }
}
/// Ziegler-Nichols PID tuning rules.
///
/// Given the ultimate gain `ku` and ultimate period `tu` (from relay or
/// frequency response experiments), computes PID gains.
pub struct ZieglerNichols;
impl ZieglerNichols {
    /// Compute PID gains using the classic Ziegler-Nichols closed-loop method.
    ///
    /// Returns `(kp, ki, kd)`.
    pub fn classic(ku: f64, tu: f64) -> (f64, f64, f64) {
        let kp = 0.6 * ku;
        let ti = 0.5 * tu;
        let td = 0.125 * tu;
        let ki = kp / ti;
        let kd = kp * td;
        (kp, ki, kd)
    }
    /// Compute PI gains using Ziegler-Nichols.
    ///
    /// Returns `(kp, ki, 0.0)`.
    pub fn pi_only(ku: f64, tu: f64) -> (f64, f64, f64) {
        let kp = 0.45 * ku;
        let ti = tu / 1.2;
        let ki = kp / ti;
        (kp, ki, 0.0)
    }
    /// Compute P-only gain using Ziegler-Nichols.
    ///
    /// Returns `(kp, 0.0, 0.0)`.
    pub fn p_only(ku: f64) -> (f64, f64, f64) {
        (0.5 * ku, 0.0, 0.0)
    }
    /// Some-overshoot tuning: reduced gains for less overshoot.
    ///
    /// Returns `(kp, ki, kd)`.
    pub fn some_overshoot(ku: f64, tu: f64) -> (f64, f64, f64) {
        let kp = ku / 3.0;
        let ti = tu / 2.0;
        let td = tu / 3.0;
        let ki = kp / ti;
        let kd = kp * td;
        (kp, ki, kd)
    }
}
/// PID controller with bumpless transfer between Auto and Manual modes.
///
/// When switching from Manual to Auto, the integral term is pre-loaded so
/// that the controller output matches the manual output without a bump.
pub struct BumplessPid {
    /// Proportional gain.
    pub kp: f64,
    /// Integral gain.
    pub ki: f64,
    /// Derivative gain.
    pub kd: f64,
    /// Anti-windup integral limit.
    pub integral_limit: f64,
    /// Output saturation limit.
    pub output_limit: f64,
    /// Accumulated integral.
    pub integral: f64,
    /// Previous error.
    pub prev_error: f64,
    /// Current operating mode.
    pub mode: PidMode,
    /// Last manual override value (used for bumpless transfer).
    pub manual_output: f64,
}
impl BumplessPid {
    /// Create a new bumpless PID controller in Auto mode.
    pub fn new(kp: f64, ki: f64, kd: f64, integral_limit: f64, output_limit: f64) -> Self {
        Self {
            kp,
            ki,
            kd,
            integral_limit,
            output_limit,
            integral: 0.0,
            prev_error: 0.0,
            mode: PidMode::Auto,
            manual_output: 0.0,
        }
    }
    /// Switch to Manual mode with the given override output.
    pub fn set_manual(&mut self, manual_output: f64) {
        self.mode = PidMode::Manual;
        self.manual_output = manual_output;
    }
    /// Switch to Auto mode with bumpless transfer.
    ///
    /// Pre-loads the integral so that the first Auto output equals the last
    /// manual output (assuming current error and derivative are zero).
    pub fn set_auto(&mut self, current_error: f64) {
        if self.mode == PidMode::Manual {
            if self.ki.abs() > 1e-15 {
                let desired_integral = (self.manual_output - self.kp * current_error) / self.ki;
                self.integral = desired_integral.clamp(-self.integral_limit, self.integral_limit);
            }
            self.prev_error = current_error;
            self.mode = PidMode::Auto;
        }
    }
    /// Compute the controller output.
    ///
    /// In Manual mode returns the manual override value.
    /// In Auto mode runs the standard PID algorithm.
    pub fn update(&mut self, error: f64, dt: f64) -> f64 {
        match self.mode {
            PidMode::Manual => self.manual_output,
            PidMode::Auto => {
                self.integral += error * dt;
                self.integral = self
                    .integral
                    .clamp(-self.integral_limit, self.integral_limit);
                let derivative = if dt > 1e-15 {
                    (error - self.prev_error) / dt
                } else {
                    0.0
                };
                self.prev_error = error;
                let raw = self.kp * error + self.ki * self.integral + self.kd * derivative;
                raw.clamp(-self.output_limit, self.output_limit)
            }
        }
    }
    /// Reset the controller to zero state.
    pub fn reset(&mut self) {
        self.integral = 0.0;
        self.prev_error = 0.0;
        self.mode = PidMode::Auto;
    }
}
/// First-order transfer function: G(s) = K / (tau*s + 1).
pub struct FirstOrderTf {
    /// Static gain.
    pub gain: f64,
    /// Time constant.
    pub tau: f64,
    /// Internal state (output value).
    pub y: f64,
}
impl FirstOrderTf {
    /// Create a new first-order transfer function.
    pub fn new(gain: f64, tau: f64) -> Self {
        Self { gain, tau, y: 0.0 }
    }
    /// Step the system with input `u` and time step `dt`.
    pub fn step(&mut self, u: f64, dt: f64) {
        self.y += (self.gain * u - self.y) / self.tau * dt;
    }
    /// Get the current output.
    pub fn output(&self) -> f64 {
        self.y
    }
    /// Reset the output to zero.
    pub fn reset(&mut self) {
        self.y = 0.0;
    }
    /// DC gain (gain at zero frequency).
    pub fn dc_gain(&self) -> f64 {
        self.gain
    }
    /// Bandwidth frequency: omega_bw = 1 / tau.
    pub fn bandwidth(&self) -> f64 {
        1.0 / self.tau
    }
}
/// Discrete-time Kalman filter for a 2-state, 1-output system.
///
/// State equation:  x\[k+1\] = A x\[k\] + B u\[k\] + w\[k\]
/// Measurement:     y\[k\]   = C x\[k\] + v\[k\]
///
/// w ~ N(0, Q), v ~ N(0, R)
pub struct KalmanFilter {
    /// Discrete-time system matrix A.
    pub a: [[f64; 2]; 2],
    /// Input matrix B.
    pub b: [f64; 2],
    /// Output matrix C.
    pub c: [f64; 2],
    /// Process noise covariance Q (2×2, stored as \[q00, q01, q10, q11\]).
    pub q: [[f64; 2]; 2],
    /// Measurement noise variance R.
    pub r: f64,
    /// Current state estimate.
    pub x_hat: [f64; 2],
    /// Current error covariance P (2×2).
    pub p: [[f64; 2]; 2],
}
impl KalmanFilter {
    /// Create a new Kalman filter.
    pub fn new(a: [[f64; 2]; 2], b: [f64; 2], c: [f64; 2], q: [[f64; 2]; 2], r: f64) -> Self {
        let p = [[1.0, 0.0], [0.0, 1.0]];
        Self {
            a,
            b,
            c,
            q,
            r,
            x_hat: [0.0; 2],
            p,
        }
    }
    /// Create a Kalman filter for a double-integrator (discretized via Euler).
    pub fn double_integrator(dt: f64, process_noise: f64, measurement_noise: f64) -> Self {
        let a = [[1.0, dt], [0.0, 1.0]];
        let b = [0.5 * dt * dt, dt];
        let c = [1.0, 0.0];
        let qv = process_noise;
        let q = [
            [qv * dt.powi(4) / 4.0, qv * dt.powi(3) / 2.0],
            [qv * dt.powi(3) / 2.0, qv * dt * dt],
        ];
        Self::new(a, b, c, q, measurement_noise * measurement_noise)
    }
    /// Predict step: propagate state and covariance.
    pub fn predict(&mut self, u: f64) {
        let x0 = self.x_hat[0];
        let x1 = self.x_hat[1];
        self.x_hat[0] = self.a[0][0] * x0 + self.a[0][1] * x1 + self.b[0] * u;
        self.x_hat[1] = self.a[1][0] * x0 + self.a[1][1] * x1 + self.b[1] * u;
        let ap = mat2_mul(&self.a, &self.p);
        let apat = mat2_mul_transpose(&ap, &self.a);
        for (p_row, (apat_row, q_row)) in self.p.iter_mut().zip(apat.iter().zip(self.q.iter())) {
            for (p_ij, (apat_ij, q_ij)) in p_row.iter_mut().zip(apat_row.iter().zip(q_row.iter())) {
                *p_ij = apat_ij + q_ij;
            }
        }
    }
    /// Update step: incorporate measurement.
    pub fn update(&mut self, measurement: f64) {
        let y = measurement - (self.c[0] * self.x_hat[0] + self.c[1] * self.x_hat[1]);
        let cp = [
            self.c[0] * self.p[0][0] + self.c[1] * self.p[1][0],
            self.c[0] * self.p[0][1] + self.c[1] * self.p[1][1],
        ];
        let s = cp[0] * self.c[0] + cp[1] * self.c[1] + self.r;
        if s.abs() < 1e-60 {
            return;
        }
        let k = [
            (self.p[0][0] * self.c[0] + self.p[0][1] * self.c[1]) / s,
            (self.p[1][0] * self.c[0] + self.p[1][1] * self.c[1]) / s,
        ];
        self.x_hat[0] += k[0] * y;
        self.x_hat[1] += k[1] * y;
        let kc = [
            [k[0] * self.c[0], k[0] * self.c[1]],
            [k[1] * self.c[0], k[1] * self.c[1]],
        ];
        let i_kc = [[1.0 - kc[0][0], -kc[0][1]], [-kc[1][0], 1.0 - kc[1][1]]];
        let p_new = mat2_mul(&i_kc, &self.p);
        self.p = p_new;
    }
    /// Get the current state estimate.
    pub fn state(&self) -> [f64; 2] {
        self.x_hat
    }
}
/// Follows a piecewise-linear position trajectory using a PID controller.
pub struct TrajectoryFollower {
    /// `(time, position)` waypoints, must be sorted by time.
    pub waypoints: Vec<[f64; 2]>,
    /// Index of the last waypoint segment entered.
    pub current_idx: usize,
    /// PID controller for tracking.
    pub pid: PidController,
}
impl TrajectoryFollower {
    /// Create a new trajectory follower.
    pub fn new(waypoints: Vec<[f64; 2]>, pid: PidController) -> Self {
        Self {
            waypoints,
            current_idx: 0,
            pid,
        }
    }
    /// Linearly interpolate the desired position at time `t`.
    pub fn desired_position(&self, t: f64) -> f64 {
        let n = self.waypoints.len();
        if n == 0 {
            return 0.0;
        }
        if t <= self.waypoints[0][0] {
            return self.waypoints[0][1];
        }
        if t >= self.waypoints[n - 1][0] {
            return self.waypoints[n - 1][1];
        }
        for i in 0..n - 1 {
            let t0 = self.waypoints[i][0];
            let t1 = self.waypoints[i + 1][0];
            if t >= t0 && t <= t1 {
                let alpha = (t - t0) / (t1 - t0);
                return self.waypoints[i][1]
                    + alpha * (self.waypoints[i + 1][1] - self.waypoints[i][1]);
            }
        }
        self.waypoints[n - 1][1]
    }
    /// Estimate the desired velocity at time `t` via finite difference.
    pub fn desired_velocity(&self, t: f64) -> f64 {
        let eps = 1e-6;
        (self.desired_position(t + eps) - self.desired_position(t - eps)) / (2.0 * eps)
    }
    /// Update the follower: compute the PID output to track the trajectory.
    pub fn update(&mut self, t: f64, current_pos: f64, _current_vel: f64, dt: f64) -> f64 {
        let desired = self.desired_position(t);
        let error = desired - current_pos;
        self.pid.update(error, dt)
    }
}
/// State-feedback (LQR-style) controller: u = -K * x.
pub struct StateFeedbackController {
    /// Gain matrix K (n_inputs × n_states).
    pub k: Vec<Vec<f64>>,
    /// Number of state variables.
    pub n_states: usize,
    /// Number of control inputs.
    pub n_inputs: usize,
}
impl StateFeedbackController {
    /// Create a new state-feedback controller from a gain matrix.
    pub fn new(k: Vec<Vec<f64>>) -> Self {
        let n_inputs = k.len();
        let n_states = if n_inputs > 0 { k[0].len() } else { 0 };
        Self {
            k,
            n_states,
            n_inputs,
        }
    }
    /// Compute the control input u = -K * x.
    pub fn control(&self, state: &[f64]) -> Vec<f64> {
        (0..self.n_inputs)
            .map(|i| {
                -self.k[i]
                    .iter()
                    .zip(state.iter())
                    .map(|(k, x)| k * x)
                    .sum::<f64>()
            })
            .collect()
    }
}
/// Operating mode of a [`JointMotor`].
#[derive(Clone, Debug, PartialEq)]
pub enum MotorMode {
    /// Track a target position using the internal PID.
    PositionControl,
    /// Track a target velocity using the internal PID.
    VelocityControl,
    /// Apply a constant torque directly.
    TorqueControl {
        /// Desired output torque (before gear-ratio and clamping).
        target_torque: f64,
    },
    /// Motor is off; outputs zero torque.
    Disabled,
}
/// Position-velocity cascade PID controller.
///
/// The outer position loop produces a desired velocity, which the inner
/// velocity loop tracks.
pub struct CascadedPid {
    /// Outer position PID.
    pub position_pid: PidController,
    /// Inner velocity PID.
    pub velocity_pid: PidController,
}
impl CascadedPid {
    /// Create a cascaded PID with separate gains for position and velocity loops.
    pub fn new(
        pos_kp: f64,
        pos_ki: f64,
        pos_kd: f64,
        vel_kp: f64,
        vel_ki: f64,
        vel_kd: f64,
    ) -> Self {
        Self {
            position_pid: PidController::new(pos_kp, pos_ki, pos_kd),
            velocity_pid: PidController::new(vel_kp, vel_ki, vel_kd),
        }
    }
    /// Update the cascade: outer loop converts position error to desired
    /// velocity; inner loop tracks the velocity error.
    pub fn update(&mut self, pos_error: f64, velocity: f64, dt: f64) -> f64 {
        let desired_velocity = self.position_pid.update(pos_error, dt);
        let vel_error = desired_velocity - velocity;
        self.velocity_pid.update(vel_error, dt)
    }
}
/// Second-order transfer function: G(s) = K*wn^2 / (s^2 + 2*zeta*wn*s + wn^2).
pub struct SecondOrderTf {
    /// Static gain K.
    pub gain: f64,
    /// Natural frequency wn.
    pub wn: f64,
    /// Damping ratio zeta.
    pub zeta: f64,
    /// Internal state-space representation.
    pub(super) ss: StateSpace2x1,
}
impl SecondOrderTf {
    /// Create a new second-order transfer function.
    pub fn new(gain: f64, wn: f64, zeta: f64) -> Self {
        let ss = StateSpace2x1::damped_oscillator(wn, zeta);
        Self { gain, wn, zeta, ss }
    }
    /// Step the system.
    pub fn step(&mut self, u: f64, dt: f64) {
        self.ss.step(self.gain * self.wn * self.wn * u, dt);
    }
    /// Get the current output.
    pub fn output(&self) -> f64 {
        self.ss.output(0.0)
    }
    /// Reset internal state.
    pub fn reset(&mut self) {
        self.ss.reset();
    }
}
/// Discrete-time state-space system simulated with ZOH input.
///
/// x\[k+1\] = Ad * x\[k\] + Bd * u\[k\]
/// y\[k\]   = C  * x\[k\]
pub struct DiscreteStateSpace {
    /// Discrete-time system matrix (n × n).
    pub ad: [[f64; 2]; 2],
    /// Discrete-time input matrix (n × 1).
    pub bd: [f64; 2],
    /// Output matrix (1 × n).
    pub c: [f64; 2],
    /// Current state.
    pub x: [f64; 2],
}
impl DiscreteStateSpace {
    /// Create a new discrete-time state-space system.
    pub fn new(ad: [[f64; 2]; 2], bd: [f64; 2], c: [f64; 2]) -> Self {
        Self {
            ad,
            bd,
            c,
            x: [0.0; 2],
        }
    }
    /// Create a discrete double-integrator via ZOH with given dt.
    pub fn double_integrator_zoh(dt: f64) -> Self {
        let a_c = [[0.0, 1.0], [0.0, 0.0]];
        let b_c = [0.0, 1.0];
        let (ad, bd) = discretize_zoh(&a_c, &b_c, dt);
        Self::new(ad, bd, [1.0, 0.0])
    }
    /// Advance one step with input u.
    pub fn step(&mut self, u: f64) {
        let x0 = self.x[0];
        let x1 = self.x[1];
        self.x[0] = self.ad[0][0] * x0 + self.ad[0][1] * x1 + self.bd[0] * u;
        self.x[1] = self.ad[1][0] * x0 + self.ad[1][1] * x1 + self.bd[1] * u;
    }
    /// Compute the output y = C * x.
    pub fn output(&self) -> f64 {
        self.c[0] * self.x[0] + self.c[1] * self.x[1]
    }
    /// Reset state to zero.
    pub fn reset(&mut self) {
        self.x = [0.0; 2];
    }
    /// Run n steps with constant input u, returning output at each step.
    pub fn simulate(&mut self, u: f64, steps: usize) -> Vec<f64> {
        (0..steps)
            .map(|_| {
                self.step(u);
                self.output()
            })
            .collect()
    }
}
/// PID controller backed by a `PidGains` struct with anti-windup.
pub struct PidControllerV2 {
    /// Controller gains.
    pub gains: PidGains,
    /// Accumulated integral term.
    pub integral: f64,
    /// Error from the previous update.
    pub prev_error: f64,
    /// Anti-windup clamp on the integral term.
    pub windup_limit: f64,
}
impl PidControllerV2 {
    /// Create a new controller with the given gains and windup limit.
    pub fn new(kp: f64, ki: f64, kd: f64, windup_limit: f64) -> Self {
        Self {
            gains: PidGains { kp, ki, kd },
            integral: 0.0,
            prev_error: 0.0,
            windup_limit,
        }
    }
    /// Compute the control output for the given error and timestep.
    pub fn update(&mut self, error: f64, dt: f64) -> f64 {
        self.integral += error * dt;
        self.integral = self.integral.clamp(-self.windup_limit, self.windup_limit);
        let derivative = if dt > 0.0 {
            (error - self.prev_error) / dt
        } else {
            0.0
        };
        self.prev_error = error;
        self.gains.kp * error + self.gains.ki * self.integral + self.gains.kd * derivative
    }
    /// Reset the controller state (integral and previous error).
    pub fn reset(&mut self) {
        self.integral = 0.0;
        self.prev_error = 0.0;
    }
}
/// Linear Quadratic Regulator for a 2-state system \[position, velocity\].
pub struct LqrController {
    /// LQR gain vector `K = [k0, k1]`.
    pub k: [f64; 2],
}
impl LqrController {
    /// Analytical LQR gains for the double-integrator plant `x'' = u`.
    ///
    /// Uses Bryson's rules:
    /// - `K[0] = sqrt(q1 / r)`
    /// - `K[1] = sqrt(2 * K[0] + q2 / r)`
    pub fn new_double_integrator(q1: f64, q2: f64, r: f64) -> Self {
        let k0 = (q1 / r).sqrt();
        let k1 = (2.0 * k0 + q2 / r).sqrt();
        Self { k: [k0, k1] }
    }
    /// Compute the control input: `u = -K · (x - x_ref)`.
    pub fn compute_input(&self, state: [f64; 2], reference: [f64; 2]) -> f64 {
        let e0 = state[0] - reference[0];
        let e1 = state[1] - reference[1];
        -(self.k[0] * e0 + self.k[1] * e1)
    }
}
/// Motorised joint with position, velocity, or torque control modes.
pub struct JointMotor {
    /// Desired joint position (used in [`MotorMode::PositionControl`]).
    pub target_position: f64,
    /// Desired joint velocity (used in [`MotorMode::VelocityControl`]).
    pub target_velocity: f64,
    /// Maximum output torque (after gear ratio).
    pub max_torque: f64,
    /// Maximum joint velocity.
    pub max_velocity: f64,
    /// Gear ratio: motor-side turns per joint-side turn.
    pub gear_ratio: f64,
    /// Internal PID controller.
    pub controller: PidController,
    /// Active control mode.
    pub mode: MotorMode,
}
impl JointMotor {
    /// Compute the torque to apply at the joint for this time-step.
    pub fn compute_torque(&mut self, current_pos: f64, current_vel: f64, dt: f64) -> f64 {
        let raw = match &self.mode.clone() {
            MotorMode::PositionControl => {
                let error = self.target_position - current_pos;
                self.controller.update(error, dt)
            }
            MotorMode::VelocityControl => {
                let error = self.target_velocity - current_vel;
                self.controller.update(error, dt)
            }
            MotorMode::TorqueControl { target_torque } => *target_torque,
            MotorMode::Disabled => 0.0,
        };
        let torque = raw * self.gear_ratio;
        torque.clamp(-self.max_torque, self.max_torque)
    }
}
/// Linear state-space system: x' = Ax + Bu, y = Cx + Du.
///
/// Fixed to 2-state, 1-input, 1-output for simplicity.
pub struct StateSpace2x1 {
    /// System matrix A (2x2).
    pub a: [[f64; 2]; 2],
    /// Input matrix B (2x1).
    pub b: [f64; 2],
    /// Output matrix C (1x2).
    pub c: [f64; 2],
    /// Feedthrough scalar D.
    pub d: f64,
    /// Current state vector.
    pub state: [f64; 2],
}
impl StateSpace2x1 {
    /// Create a new state-space system.
    pub fn new(a: [[f64; 2]; 2], b: [f64; 2], c: [f64; 2], d: f64) -> Self {
        Self {
            a,
            b,
            c,
            d,
            state: [0.0; 2],
        }
    }
    /// Create a double-integrator system (x'' = u).
    ///
    /// State: \[position, velocity\], input: acceleration.
    pub fn double_integrator() -> Self {
        Self::new([[0.0, 1.0], [0.0, 0.0]], [0.0, 1.0], [1.0, 0.0], 0.0)
    }
    /// Create a damped oscillator: x'' + 2*zeta*wn*x' + wn^2*x = u.
    pub fn damped_oscillator(wn: f64, zeta: f64) -> Self {
        Self::new(
            [[0.0, 1.0], [-wn * wn, -2.0 * zeta * wn]],
            [0.0, 1.0],
            [1.0, 0.0],
            0.0,
        )
    }
    /// Advance the system one step using forward Euler.
    pub fn step(&mut self, u: f64, dt: f64) {
        let x0 = self.state[0];
        let x1 = self.state[1];
        self.state[0] = x0 + (self.a[0][0] * x0 + self.a[0][1] * x1 + self.b[0] * u) * dt;
        self.state[1] = x1 + (self.a[1][0] * x0 + self.a[1][1] * x1 + self.b[1] * u) * dt;
    }
    /// Compute the output y = Cx + Du.
    pub fn output(&self, u: f64) -> f64 {
        self.c[0] * self.state[0] + self.c[1] * self.state[1] + self.d * u
    }
    /// Reset state to zero.
    pub fn reset(&mut self) {
        self.state = [0.0; 2];
    }
    /// Check if the system is stable (both eigenvalues have negative real part).
    ///
    /// For a 2x2 matrix, eigenvalues are roots of:
    /// lambda^2 - tr(A)*lambda + det(A) = 0
    /// Stable if tr(A) < 0 and det(A) > 0.
    pub fn is_stable(&self) -> bool {
        let tr = self.a[0][0] + self.a[1][1];
        let det = self.a[0][0] * self.a[1][1] - self.a[0][1] * self.a[1][0];
        tr < 0.0 && det > 0.0
    }
    /// Eigenvalues of the A matrix (may be complex; returns (real, imag) pairs).
    pub fn eigenvalues(&self) -> [(f64, f64); 2] {
        let tr = self.a[0][0] + self.a[1][1];
        let det = self.a[0][0] * self.a[1][1] - self.a[0][1] * self.a[1][0];
        let disc = tr * tr - 4.0 * det;
        if disc >= 0.0 {
            let sq = disc.sqrt();
            [((tr + sq) / 2.0, 0.0), ((tr - sq) / 2.0, 0.0)]
        } else {
            let sq = (-disc).sqrt();
            [(tr / 2.0, sq / 2.0), (tr / 2.0, -sq / 2.0)]
        }
    }
}
/// A single point on a Bode plot.
pub struct BodePoint {
    /// Frequency (rad/s).
    pub freq: f64,
    /// Magnitude in decibels.
    pub magnitude_db: f64,
    /// Phase in degrees.
    pub phase_deg: f64,
}
/// Mode of a PID controller for bumpless transfer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PidMode {
    /// Controller is actively computing outputs.
    Auto,
    /// Controller output is being overridden by an external signal.
    Manual,
}
/// Stability margins for an open-loop transfer function.
pub struct StabilityMargins {
    /// Gain margin in decibels (how much gain can increase before instability).
    pub gain_margin_db: f64,
    /// Phase margin in degrees (how much phase can decrease before instability).
    pub phase_margin_deg: f64,
    /// Phase-crossover frequency ωpc (rad/s): where phase = -180°.
    pub phase_crossover_freq: f64,
    /// Gain-crossover frequency ωgc (rad/s): where |G| = 1 (0 dB).
    pub gain_crossover_freq: f64,
}
/// MRAC-inspired adaptive PID that adjusts `kp` based on tracking error.
pub struct AdaptivePid {
    /// Underlying PID controller.
    pub base_pid: PidController,
    /// Rate at which `kp` is adapted.
    pub adaptation_rate: f64,
}
impl AdaptivePid {
    /// Adapt the proportional gain based on the current tracking error.
    ///
    /// If the error is growing (error times previous error > 0), `kp` increases.
    pub fn update_gains(&mut self, tracking_error: f64, dt: f64) {
        let prev = self.base_pid.prev_error;
        if tracking_error * prev > 0.0 {
            self.base_pid.kp += self.adaptation_rate * tracking_error.abs() * dt;
        }
    }
    /// Compute the PID output, then adapt gains.
    pub fn update(&mut self, error: f64, dt: f64) -> f64 {
        let output = self.base_pid.update(error, dt);
        self.update_gains(error, dt);
        output
    }
}
/// Proportional-Integral-Derivative gain triplet.
pub struct PidGains {
    /// Proportional gain.
    pub kp: f64,
    /// Integral gain.
    pub ki: f64,
    /// Derivative gain.
    pub kd: f64,
}
/// Stateless feed-forward compensation functions.
pub struct FeedforwardController;
impl FeedforwardController {
    /// Gravity compensation torque for a revolute arm/pendulum.
    ///
    /// `torque = mass * g * L * cos(angle)`
    ///
    /// The caller supplies the effective link length via `g` (i.e. pass `g * L`
    /// or use `g = 9.81` and treat the result as normalised by `L`).
    pub fn gravity_compensation(mass: f64, g: f64, angle: f64) -> f64 {
        mass * g * angle.cos()
    }
    /// Inertia feed-forward: `torque = I * alpha`.
    pub fn inertia_compensation(inertia: f64, desired_acc: f64) -> f64 {
        inertia * desired_acc
    }
    /// Friction compensation.
    ///
    /// Returns static friction scaled by sign(v) when `|v| < threshold`,
    /// otherwise kinetic friction scaled by sign(v).
    pub fn friction_compensation(
        static_friction: f64,
        kinetic_friction: f64,
        velocity: f64,
        threshold: f64,
    ) -> f64 {
        let sign = if velocity > 0.0 {
            1.0
        } else if velocity < 0.0 {
            -1.0
        } else {
            0.0
        };
        if velocity.abs() < threshold {
            static_friction * sign
        } else {
            kinetic_friction * sign
        }
    }
}
/// Continuous-time LTI state-space model: dx/dt = Ax + Bu, y = Cx.
pub struct StateSpace {
    /// System matrix (n × n).
    pub a: Vec<Vec<f64>>,
    /// Input matrix (n × m).
    pub b: Vec<Vec<f64>>,
    /// Output matrix (p × n).
    pub c: Vec<Vec<f64>>,
    /// Number of states.
    pub n: usize,
    /// Number of inputs.
    pub m: usize,
    /// Number of outputs.
    pub p: usize,
}
impl StateSpace {
    /// Create a new state-space model.
    pub fn new(a: Vec<Vec<f64>>, b: Vec<Vec<f64>>, c: Vec<Vec<f64>>) -> Self {
        let n = a.len();
        let m = if n > 0 { b[0].len() } else { 0 };
        let p = c.len();
        Self { a, b, c, n, m, p }
    }
    /// Euler integration: x += (a*x + b*u) * dt.
    pub fn step(&mut self, x: &mut [f64], u: &[f64], dt: f64) {
        let mut dx = vec![0.0; self.n];
        for (i, (dx_i, (a_row, b_row))) in dx
            .iter_mut()
            .zip(self.a.iter().zip(self.b.iter()))
            .enumerate()
        {
            let _ = i;
            for (a_ij, x_j) in a_row.iter().zip(x.iter()) {
                *dx_i += a_ij * x_j;
            }
            for (b_ij, u_j) in b_row.iter().zip(u.iter()) {
                *dx_i += b_ij * u_j;
            }
        }
        for (x_i, dx_i) in x.iter_mut().zip(dx.iter()) {
            *x_i += dx_i * dt;
        }
    }
    /// Compute output y = C*x (D=0).
    pub fn output(&self, x: &[f64], _u: &[f64]) -> Vec<f64> {
        (0..self.p)
            .map(|i| (0..self.n).map(|j| self.c[i][j] * x[j]).sum())
            .collect()
    }
}
/// Discrete-time Luenberger state observer.
///
/// x_hat\[k+1\] = (A - L C) x_hat\[k\] + B u\[k\] + L y\[k\]
///
/// where L is the observer gain vector.
pub struct LuenbergerObserver {
    /// System matrix A (2×2).
    pub a: [[f64; 2]; 2],
    /// Input matrix B (2×1).
    pub b: [f64; 2],
    /// Output matrix C (1×2).
    pub c: [f64; 2],
    /// Observer gain L (2×1).
    pub l: [f64; 2],
    /// State estimate.
    pub x_hat: [f64; 2],
}
impl LuenbergerObserver {
    /// Create a new Luenberger observer.
    pub fn new(a: [[f64; 2]; 2], b: [f64; 2], c: [f64; 2], l: [f64; 2]) -> Self {
        Self {
            a,
            b,
            c,
            l,
            x_hat: [0.0; 2],
        }
    }
    /// Update observer with latest input u and measurement y.
    pub fn update(&mut self, u: f64, y: f64) {
        let x0 = self.x_hat[0];
        let x1 = self.x_hat[1];
        let y_hat = self.c[0] * x0 + self.c[1] * x1;
        let innov = y - y_hat;
        self.x_hat[0] = self.a[0][0] * x0 + self.a[0][1] * x1 + self.b[0] * u + self.l[0] * innov;
        self.x_hat[1] = self.a[1][0] * x0 + self.a[1][1] * x1 + self.b[1] * u + self.l[1] * innov;
    }
    /// Return the current state estimate.
    pub fn estimate(&self) -> [f64; 2] {
        self.x_hat
    }
    /// Reset state estimate to zero.
    pub fn reset(&mut self) {
        self.x_hat = [0.0; 2];
    }
    /// Check observer stability: eigenvalues of (A - L C) should be inside unit circle.
    pub fn is_stable(&self) -> bool {
        let ao = [
            [
                self.a[0][0] - self.l[0] * self.c[0],
                self.a[0][1] - self.l[0] * self.c[1],
            ],
            [
                self.a[1][0] - self.l[1] * self.c[0],
                self.a[1][1] - self.l[1] * self.c[1],
            ],
        ];
        let tr = ao[0][0] + ao[1][1];
        let det = ao[0][0] * ao[1][1] - ao[0][1] * ao[1][0];
        let disc = tr * tr - 4.0 * det;
        if disc < 0.0 {
            det.abs() < 1.0
        } else {
            let sqrt_d = disc.sqrt();
            let l1 = ((tr + sqrt_d) / 2.0).abs();
            let l2 = ((tr - sqrt_d) / 2.0).abs();
            l1 < 1.0 && l2 < 1.0
        }
    }
}
