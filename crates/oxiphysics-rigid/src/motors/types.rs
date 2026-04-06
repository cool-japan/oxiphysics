//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
/// Enhanced thermal model for a motor with winding and housing nodes.
///
/// Two-node lumped model:
/// - Node 1: winding (temperature rises fastest)
/// - Node 2: housing / case (coupled to ambient)
#[allow(dead_code)]
pub struct MotorThermalTwoNode {
    /// Winding temperature (°C).
    pub t_winding: f64,
    /// Case temperature (°C).
    pub t_case: f64,
    /// Ambient temperature (°C).
    pub t_ambient: f64,
    /// Maximum winding temperature (°C).
    pub t_max_winding: f64,
    /// Thermal resistance winding→case (°C/W).
    pub r_wc: f64,
    /// Thermal resistance case→ambient (°C/W).
    pub r_ca: f64,
    /// Winding thermal capacitance (J/°C).
    pub c_winding: f64,
    /// Case thermal capacitance (J/°C).
    pub c_case: f64,
}
impl MotorThermalTwoNode {
    /// Create a two-node thermal model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        t_ambient: f64,
        t_max_winding: f64,
        r_wc: f64,
        r_ca: f64,
        c_winding: f64,
        c_case: f64,
    ) -> Self {
        Self {
            t_winding: t_ambient,
            t_case: t_ambient,
            t_ambient,
            t_max_winding,
            r_wc,
            r_ca,
            c_winding,
            c_case,
        }
    }
    /// Advance by one time step with copper-loss power `p_loss` (W).
    pub fn step(&mut self, p_loss: f64, dt: f64) {
        let q_wc = (self.t_winding - self.t_case) / self.r_wc;
        let q_ca = (self.t_case - self.t_ambient) / self.r_ca;
        let dt_winding = (p_loss - q_wc) / self.c_winding * dt;
        let dt_case = (q_wc - q_ca) / self.c_case * dt;
        self.t_winding += dt_winding;
        self.t_case += dt_case;
    }
    /// Winding temperature margin before trip (°C).
    pub fn thermal_margin(&self) -> f64 {
        (self.t_max_winding - self.t_winding).max(0.0)
    }
    /// True when winding exceeds maximum temperature.
    pub fn is_overheated(&self) -> bool {
        self.t_winding >= self.t_max_winding
    }
    /// Steady-state winding temperature rise for constant loss `p` (°C above ambient).
    pub fn steady_state_winding_rise(&self, p_loss: f64) -> f64 {
        p_loss * (self.r_wc + self.r_ca)
    }
    /// Steady-state case temperature rise for constant loss `p` (°C above ambient).
    pub fn steady_state_case_rise(&self, p_loss: f64) -> f64 {
        p_loss * self.r_ca
    }
    /// Reset both nodes to ambient.
    pub fn reset(&mut self) {
        self.t_winding = self.t_ambient;
        self.t_case = self.t_ambient;
    }
}
/// Piecewise-linear motor efficiency curve.
///
/// Stores (speed, efficiency) data points and interpolates between them.
#[allow(dead_code)]
pub struct MotorEfficiencyCurve {
    /// Speed sample points (rad/s), must be sorted ascending.
    pub speeds: Vec<f64>,
    /// Efficiency at each speed point (0..1).
    pub efficiencies: Vec<f64>,
}
impl MotorEfficiencyCurve {
    /// Create from parallel speed and efficiency vectors.
    ///
    /// Panics in debug if lengths differ.
    pub fn new(speeds: Vec<f64>, efficiencies: Vec<f64>) -> Self {
        debug_assert_eq!(
            speeds.len(),
            efficiencies.len(),
            "speeds and efficiencies must match"
        );
        Self {
            speeds,
            efficiencies,
        }
    }
    /// A flat 85% efficiency curve for quick testing.
    pub fn flat(eta: f64) -> Self {
        Self {
            speeds: vec![0.0, 1e6],
            efficiencies: vec![eta, eta],
        }
    }
    /// Interpolate efficiency at `omega` (rad/s).
    pub fn efficiency_at(&self, omega: f64) -> f64 {
        let n = self.speeds.len();
        if n == 0 {
            return 0.0;
        }
        if omega <= self.speeds[0] {
            return self.efficiencies[0];
        }
        if omega >= self.speeds[n - 1] {
            return self.efficiencies[n - 1];
        }
        let idx = self.speeds.partition_point(|&s| s <= omega);
        let i = idx.saturating_sub(1);
        let j = i + 1;
        let t = (omega - self.speeds[i]) / (self.speeds[j] - self.speeds[i]);
        self.efficiencies[i] + t * (self.efficiencies[j] - self.efficiencies[i])
    }
    /// Output power given input power and speed.
    pub fn output_power(&self, p_input: f64, omega: f64) -> f64 {
        p_input * self.efficiency_at(omega)
    }
    /// Speed at peak efficiency (O(n) scan).
    pub fn peak_efficiency_speed(&self) -> f64 {
        self.speeds
            .iter()
            .zip(self.efficiencies.iter())
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(s, _)| *s)
            .unwrap_or(0.0)
    }
}
/// Proportional velocity controller for joint motors.
pub struct VelocityMotor {
    /// Desired joint velocity.
    pub target_velocity: f64,
    /// Current joint velocity.
    pub current_velocity: f64,
    /// Proportional gain.
    pub k_p: f64,
    /// Maximum allowable torque magnitude.
    pub max_torque: f64,
}
impl VelocityMotor {
    /// Creates a new `VelocityMotor` with the given gain and torque limit.
    pub fn new(k_p: f64, max_torque: f64) -> Self {
        Self {
            target_velocity: 0.0,
            current_velocity: 0.0,
            k_p,
            max_torque,
        }
    }
    /// Computes the control torque.
    pub fn compute_torque(&self) -> f64 {
        let raw = self.k_p * (self.target_velocity - self.current_velocity);
        raw.clamp(-self.max_torque, self.max_torque)
    }
    /// Updates the current velocity.
    pub fn update_velocity(&mut self, velocity: f64) {
        self.current_velocity = velocity;
    }
    /// Velocity error.
    pub fn velocity_error(&self) -> f64 {
        self.target_velocity - self.current_velocity
    }
}
/// Full PID (Proportional-Integral-Derivative) controller.
///
/// The controller integrates error over time and uses a finite-difference
/// derivative.  Output is clamped to `+-max_output`.
pub struct PidController {
    /// Proportional gain.
    pub kp: f64,
    /// Integral gain.
    pub ki: f64,
    /// Derivative gain.
    pub kd: f64,
    /// Accumulated integral of the error.
    pub integral: f64,
    /// Error from the previous `update` call.
    pub prev_error: f64,
    /// Maximum absolute output value.
    pub max_output: f64,
    /// Maximum absolute integral value (anti-windup).
    pub max_integral: f64,
}
impl PidController {
    /// Create a new PID controller.
    pub fn new(kp: f64, ki: f64, kd: f64, max_output: f64) -> Self {
        Self {
            kp,
            ki,
            kd,
            integral: 0.0,
            prev_error: 0.0,
            max_output,
            max_integral: max_output * 10.0,
        }
    }
    /// Create a PID controller with anti-windup integral clamp.
    pub fn with_integral_clamp(mut self, max_integral: f64) -> Self {
        self.max_integral = max_integral;
        self
    }
    /// Compute the control output for the given `error` and time step `dt`.
    ///
    /// Returns the output clamped to `+-max_output`.
    pub fn update(&mut self, error: f64, dt: f64) -> f64 {
        self.integral += error * dt;
        self.integral = self.integral.clamp(-self.max_integral, self.max_integral);
        let derivative = if dt > 1e-15 {
            (error - self.prev_error) / dt
        } else {
            0.0
        };
        self.prev_error = error;
        let raw = self.kp * error + self.ki * self.integral + self.kd * derivative;
        raw.clamp(-self.max_output, self.max_output)
    }
    /// Reset the controller state (integral and previous error).
    pub fn reset(&mut self) {
        self.integral = 0.0;
        self.prev_error = 0.0;
    }
    /// Apply Ziegler-Nichols tuning from the ultimate gain `ku` and oscillation
    /// period `tu`.
    pub fn ziegler_nichols_tune(&mut self, ku: f64, tu: f64) {
        self.kp = 0.6 * ku;
        self.ki = 2.0 * self.kp / tu;
        self.kd = self.kp * tu / 8.0;
    }
    /// Apply Cohen-Coon tuning from process parameters.
    ///
    /// * `k_process` - process gain
    /// * `tau` - process time constant
    /// * `theta` - process dead time
    pub fn cohen_coon_tune(&mut self, k_process: f64, tau: f64, theta: f64) {
        let r = theta / tau;
        self.kp = (1.0 / (k_process * r)) * (1.35 + 0.25 * r);
        self.ki = self.kp * (2.5 * theta) / (1.0 + 0.6 * r);
        self.kd = self.kp * 0.37 * theta / (1.0 + 0.2 * r);
    }
    /// Compute individual P, I, D contributions for diagnostics.
    pub fn contributions(&self, error: f64, dt: f64) -> PidContributions {
        let derivative = if dt > 1e-15 {
            (error - self.prev_error) / dt
        } else {
            0.0
        };
        PidContributions {
            proportional: self.kp * error,
            integral: self.ki * self.integral,
            derivative: self.kd * derivative,
        }
    }
}
/// DC brushed motor electrical model.
///
/// The motor is described by:
/// - V_a = R_a * I_a + L_a * dI_a/dt + K_e * ω
/// - τ = K_t * I_a
/// - J * dω/dt = τ - B * ω - τ_load
#[derive(Debug, Clone)]
pub struct DcMotor {
    /// Armature resistance (Ω).
    pub resistance: f64,
    /// Armature inductance (H).
    pub inductance: f64,
    /// Back-EMF constant (V·s/rad).
    pub ke: f64,
    /// Torque constant (N·m/A).
    pub kt: f64,
    /// Rotor inertia (kg·m²).
    pub rotor_inertia: f64,
    /// Viscous friction coefficient (N·m·s/rad).
    pub viscous_friction: f64,
    /// Current armature current (A).
    pub current: f64,
    /// Current angular velocity (rad/s).
    pub omega: f64,
}
impl DcMotor {
    /// Create a new DC motor.
    #[allow(clippy::too_many_arguments)]
    #[allow(non_snake_case)]
    pub fn new(
        resistance: f64,
        inductance: f64,
        ke: f64,
        kt: f64,
        rotor_inertia: f64,
        viscous_friction: f64,
    ) -> Self {
        Self {
            resistance,
            inductance,
            ke,
            kt,
            rotor_inertia,
            viscous_friction,
            current: 0.0,
            omega: 0.0,
        }
    }
    /// Integrate one time step with applied voltage `v_a` and load torque `tau_load`.
    pub fn step(&mut self, v_a: f64, tau_load: f64, dt: f64) {
        let v_emf = self.ke * self.omega;
        let di_dt =
            (v_a - self.resistance * self.current - v_emf) / self.inductance.max(f64::EPSILON);
        self.current += di_dt * dt;
        let tau = self.kt * self.current;
        let domega_dt = (tau - self.viscous_friction * self.omega - tau_load)
            / self.rotor_inertia.max(f64::EPSILON);
        self.omega += domega_dt * dt;
    }
    /// Compute back-EMF at current speed.
    pub fn back_emf(&self) -> f64 {
        self.ke * self.omega
    }
    /// Compute current torque output.
    pub fn torque(&self) -> f64 {
        self.kt * self.current
    }
    /// No-load speed at voltage `v` (rad/s).
    pub fn no_load_speed(&self, v: f64) -> f64 {
        v / self.ke.max(f64::EPSILON)
    }
    /// Stall torque at voltage `v` (N·m).
    pub fn stall_torque(&self, v: f64) -> f64 {
        self.kt * v / self.resistance.max(f64::EPSILON)
    }
}
impl DcMotor {
    /// Mechanical output power (W) = torque × angular velocity.
    pub fn mechanical_power(&self) -> f64 {
        self.torque() * self.omega
    }
    /// Electrical input power (W) = V_a × I.
    pub fn electrical_power(&self, v_a: f64) -> f64 {
        v_a * self.current
    }
    /// Copper loss (Ohmic loss, W) = I² × R.
    pub fn copper_loss(&self) -> f64 {
        self.current * self.current * self.resistance
    }
    /// Motor efficiency η = P_mech / P_elec (clamped to \[0, 1\]).
    pub fn efficiency(&self, v_a: f64) -> f64 {
        let p_elec = self.electrical_power(v_a).abs();
        let p_mech = self.mechanical_power().abs();
        if p_elec < 1e-30 {
            return 0.0;
        }
        (p_mech / p_elec).min(1.0)
    }
    /// Motor operating point: (torque, speed) at steady state for voltage `v`.
    ///
    /// Returns (τ, ω) assuming zero load torque.
    pub fn steady_state_no_load(&self, v: f64) -> (f64, f64) {
        let omega_nl = self.no_load_speed(v);
        (0.0, omega_nl)
    }
    /// Motor operating point under load torque `tau_load` at voltage `v`.
    ///
    /// Assumes steady-state (dI/dt = 0, dω/dt = 0).
    pub fn steady_state(&self, v: f64, tau_load: f64) -> (f64, f64) {
        let i_ss = tau_load / self.kt.max(f64::EPSILON);
        let omega_ss = (v - self.resistance * i_ss) / self.ke.max(f64::EPSILON);
        (tau_load, omega_ss.max(0.0))
    }
    /// Speed-torque slope: dω/dτ = -R / (kt * ke).
    pub fn speed_torque_slope(&self) -> f64 {
        -(self.resistance) / (self.kt * self.ke).max(f64::EPSILON)
    }
    /// Reset motor state (current and speed to zero).
    pub fn reset_state(&mut self) {
        self.current = 0.0;
        self.omega = 0.0;
    }
    /// Electrical time constant τ_e = L/R.
    pub fn electrical_time_constant(&self) -> f64 {
        self.inductance / self.resistance.max(f64::EPSILON)
    }
    /// Mechanical time constant τ_m = J*R / (kt*ke).
    pub fn mechanical_time_constant(&self) -> f64 {
        self.rotor_inertia * self.resistance / (self.kt * self.ke).max(f64::EPSILON)
    }
}
/// Permanent Magnet Synchronous Motor in the d-q rotating reference frame.
///
/// The d-q model eliminates time-varying mutual inductances and provides a
/// convenient framework for field-oriented control (FOC).
///
/// # Electrical equations (continuous time)
///
/// ```text
/// v_d = R·i_d + L_d·(di_d/dt) − ω_e·L_q·i_q
/// v_q = R·i_q + L_q·(di_q/dt) + ω_e·(L_d·i_d + λ_pm)
/// ```
///
/// # Torque equation
///
/// ```text
/// τ = (3/2) · P/2 · [λ_pm·i_q + (L_d − L_q)·i_d·i_q]
/// ```
///
/// where `P` is the number of poles.
#[allow(dead_code)]
#[allow(non_snake_case)]
pub struct PmsmMotor {
    /// Number of poles (must be even).
    pub P: u32,
    /// Stator resistance per phase (Ω).
    pub R_s: f64,
    /// d-axis inductance (H).
    pub L_d: f64,
    /// q-axis inductance (H).
    pub L_q: f64,
    /// Permanent magnet flux linkage λ_pm (Wb).
    pub lambda_pm: f64,
    /// d-axis current (A).
    pub i_d: f64,
    /// q-axis current (A).
    pub i_q: f64,
    /// Electrical angular velocity ω_e (rad/s).
    pub omega_e: f64,
    /// Mechanical angular velocity ω_m (rad/s).
    pub omega_m: f64,
    /// Rotor inertia (kg·m²).
    pub J: f64,
    /// Viscous friction (N·m·s/rad).
    pub B: f64,
}
impl PmsmMotor {
    /// Create a new PMSM.
    #[allow(clippy::too_many_arguments)]
    #[allow(non_snake_case)]
    pub fn new(P: u32, R_s: f64, L_d: f64, L_q: f64, lambda_pm: f64, J: f64, B: f64) -> Self {
        Self {
            P,
            R_s,
            L_d,
            L_q,
            lambda_pm,
            i_d: 0.0,
            i_q: 0.0,
            omega_e: 0.0,
            omega_m: 0.0,
            J,
            B,
        }
    }
    /// Electrical-to-mechanical speed ratio.
    pub fn pole_pairs(&self) -> f64 {
        self.P as f64 / 2.0
    }
    /// Update electrical speed from mechanical speed.
    pub fn sync_electrical_speed(&mut self) {
        self.omega_e = self.omega_m * self.pole_pairs();
    }
    /// Electromagnetic torque from the d-q currents (N·m).
    pub fn torque(&self) -> f64 {
        let pp = self.pole_pairs();
        (3.0 / 2.0) * pp * (self.lambda_pm * self.i_q + (self.L_d - self.L_q) * self.i_d * self.i_q)
    }
    /// d-axis voltage required to maintain `i_d` at steady state.
    pub fn steady_state_vd(&self) -> f64 {
        self.R_s * self.i_d - self.omega_e * self.L_q * self.i_q
    }
    /// q-axis voltage required to maintain `i_q` at steady state.
    pub fn steady_state_vq(&self) -> f64 {
        self.R_s * self.i_q + self.omega_e * (self.L_d * self.i_d + self.lambda_pm)
    }
    /// Maximum torque per ampere (MTPA) current angle (rad).
    ///
    /// For a non-salient machine (L_d = L_q), MTPA sets i_d = 0.
    /// For salient machines this returns the optimal current angle `β`:
    ///
    /// ```text
    /// tan(2β) = −λ_pm / ((L_d − L_q) · |I|)
    /// ```
    ///
    /// `I_magnitude` is the total current amplitude in A.
    #[allow(non_snake_case)]
    pub fn mtpa_angle(&self, I_magnitude: f64) -> f64 {
        let saliency = self.L_d - self.L_q;
        if saliency.abs() < 1e-12 {
            return 0.0;
        }
        0.5 * (-self.lambda_pm / (saliency * I_magnitude)).atan()
    }
    /// Set i_d and i_q from a total current magnitude and MTPA angle.
    #[allow(non_snake_case)]
    pub fn set_mtpa_currents(&mut self, I_magnitude: f64) {
        let beta = self.mtpa_angle(I_magnitude);
        self.i_d = I_magnitude * beta.sin();
        self.i_q = I_magnitude * beta.cos();
    }
    /// Flux-weakening: reduce i_d to operate above base speed.
    ///
    /// Sets i_d to limit the stator voltage magnitude to `V_max`:
    ///
    /// ```text
    /// i_d = (V_max / ω_e − λ_pm) / L_d
    /// ```
    ///
    /// Clamps i_d to zero if it would become positive (no field strengthening).
    #[allow(non_snake_case)]
    pub fn flux_weaken(&mut self, V_max: f64) {
        if self.omega_e.abs() < 1e-12 {
            self.i_d = 0.0;
            return;
        }
        let i_d_fw = (V_max / self.omega_e.abs() - self.lambda_pm) / self.L_d;
        self.i_d = i_d_fw.min(0.0);
    }
    /// Integrate rotor dynamics for timestep `dt` (s) with applied load
    /// torque `tau_load` (N·m) opposing motion.
    pub fn integrate_mechanics(&mut self, tau_load: f64, dt: f64) {
        let tau_em = self.torque();
        let alpha = (tau_em - tau_load - self.B * self.omega_m) / self.J;
        self.omega_m += alpha * dt;
        self.sync_electrical_speed();
    }
    /// Integrate d-axis current using forward Euler.
    ///
    /// `v_d` is the commanded d-axis voltage.
    pub fn integrate_id(&mut self, v_d: f64, dt: f64) {
        let di_d = (v_d - self.R_s * self.i_d + self.omega_e * self.L_q * self.i_q) / self.L_d;
        self.i_d += di_d * dt;
    }
    /// Integrate q-axis current using forward Euler.
    ///
    /// `v_q` is the commanded q-axis voltage.
    pub fn integrate_iq(&mut self, v_q: f64, dt: f64) {
        let di_q =
            (v_q - self.R_s * self.i_q - self.omega_e * (self.L_d * self.i_d + self.lambda_pm))
                / self.L_q;
        self.i_q += di_q * dt;
    }
    /// Copper loss: P_cu = (3/2) · R_s · (i_d² + i_q²)  \[W\]
    pub fn copper_loss(&self) -> f64 {
        1.5 * self.R_s * (self.i_d * self.i_d + self.i_q * self.i_q)
    }
    /// Back-EMF magnitude |E| = ω_e · λ_pm  \[V\]
    pub fn back_emf(&self) -> f64 {
        self.omega_e.abs() * self.lambda_pm
    }
}
/// DC motor operating point with efficiency lookup table.
///
/// Stores a 2D grid of (torque, speed) → efficiency values.
/// Bilinear interpolation is used for points between grid nodes.
#[derive(Debug, Clone)]
pub struct DcMotorEfficiencyMap {
    /// Torque breakpoints (N·m), sorted ascending.
    pub torque_points: Vec<f64>,
    /// Speed breakpoints (rad/s), sorted ascending.
    pub speed_points: Vec<f64>,
    /// Efficiency grid: `grid[i][j]` = efficiency at (torque_points\[i\], speed_points\[j\]).
    pub grid: Vec<Vec<f64>>,
    /// No-load speed (rad/s).
    pub no_load_speed: f64,
    /// Stall torque (N·m).
    pub stall_torque: f64,
}
impl DcMotorEfficiencyMap {
    /// Create a flat efficiency map (uniform efficiency for all operating points).
    pub fn uniform(
        no_load_speed: f64,
        stall_torque: f64,
        efficiency: f64,
        n_torque: usize,
        n_speed: usize,
    ) -> Self {
        let torque_points: Vec<f64> = (0..n_torque)
            .map(|i| stall_torque * i as f64 / (n_torque - 1).max(1) as f64)
            .collect();
        let speed_points: Vec<f64> = (0..n_speed)
            .map(|i| no_load_speed * i as f64 / (n_speed - 1).max(1) as f64)
            .collect();
        let grid = vec![vec![efficiency.clamp(0.0, 1.0); n_speed]; n_torque];
        Self {
            torque_points,
            speed_points,
            grid,
            no_load_speed,
            stall_torque,
        }
    }
    /// Look up efficiency at the given (torque, speed) via bilinear interpolation.
    pub fn efficiency_at(&self, torque: f64, speed: f64) -> f64 {
        let t = torque.clamp(0.0, *self.torque_points.last().unwrap_or(&0.0));
        let s = speed.clamp(0.0, *self.speed_points.last().unwrap_or(&0.0));
        let ti = self.find_index(&self.torque_points, t);
        let si = self.find_index(&self.speed_points, s);
        if self.torque_points.len() < 2 || self.speed_points.len() < 2 {
            return self.grid[0][0];
        }
        let t0 = self.torque_points[ti];
        let t1 = self.torque_points[(ti + 1).min(self.torque_points.len() - 1)];
        let s0 = self.speed_points[si];
        let s1 = self.speed_points[(si + 1).min(self.speed_points.len() - 1)];
        let ft = if (t1 - t0).abs() > 1e-30 {
            (t - t0) / (t1 - t0)
        } else {
            0.0
        };
        let fs = if (s1 - s0).abs() > 1e-30 {
            (s - s0) / (s1 - s0)
        } else {
            0.0
        };
        let j1 = (si + 1).min(self.speed_points.len() - 1);
        let i1 = (ti + 1).min(self.torque_points.len() - 1);
        let e00 = self.grid[ti][si];
        let e10 = self.grid[i1][si];
        let e01 = self.grid[ti][j1];
        let e11 = self.grid[i1][j1];
        let e = (1.0 - ft) * (1.0 - fs) * e00
            + ft * (1.0 - fs) * e10
            + (1.0 - ft) * fs * e01
            + ft * fs * e11;
        e.clamp(0.0, 1.0)
    }
    fn find_index(&self, pts: &[f64], val: f64) -> usize {
        if pts.len() < 2 {
            return 0;
        }
        for i in 0..pts.len() - 1 {
            if val <= pts[i + 1] {
                return i;
            }
        }
        pts.len() - 2
    }
    /// Output power (W) at the given operating point.
    pub fn output_power(&self, torque: f64, speed: f64) -> f64 {
        torque * speed
    }
    /// Input power (W) including losses.
    pub fn input_power(&self, torque: f64, speed: f64) -> f64 {
        let eff = self.efficiency_at(torque, speed);
        if eff > 1e-12 {
            self.output_power(torque, speed) / eff
        } else {
            0.0
        }
    }
}
/// S-curve (jerk-limited) velocity profile generator.
///
/// Limits jerk to produce smooth acceleration transitions, reducing mechanical
/// vibration in precision positioning systems.
#[allow(dead_code)]
pub struct SCurveProfile {
    /// Maximum velocity (units/s).
    pub v_max: f64,
    /// Maximum acceleration (units/s²).
    pub a_max: f64,
    /// Maximum jerk (units/s³).
    pub j_max: f64,
    /// Start position.
    pub x_start: f64,
    /// Target position.
    pub x_target: f64,
    /// Direction sign (+1 / -1).
    pub direction: f64,
    /// Duration of jerk-up segment (s).
    pub t_j: f64,
    /// Duration of constant-acceleration segment (s).
    pub t_a: f64,
    /// Duration of cruise segment (s).
    pub t_c: f64,
    /// Total move duration (s).
    pub t_total: f64,
    /// Actual peak velocity.
    pub v_cruise: f64,
    /// Actual peak acceleration.
    pub a_peak: f64,
}
impl SCurveProfile {
    /// Plan an S-curve move from `x_start` to `x_target`.
    pub fn plan(x_start: f64, x_target: f64, v_max: f64, a_max: f64, j_max: f64) -> Self {
        let dist = x_target - x_start;
        let direction = if dist >= 0.0 { 1.0 } else { -1.0 };
        let s = dist.abs();
        let t_j = a_max / j_max;
        let v_jerk = 0.5 * j_max * t_j * t_j;
        let v_remaining = (v_max - 2.0 * v_jerk).max(0.0);
        let t_a = v_remaining / a_max;
        let a_peak = a_max;
        let v_cruise = 2.0 * v_jerk + v_remaining;
        let s_jerk_up = (1.0 / 6.0) * j_max * t_j * t_j * t_j;
        let s_const_a = v_jerk * t_a + 0.5 * a_max * t_a * t_a;
        let s_jerk_down = v_cruise * t_j - (1.0 / 6.0) * j_max * t_j * t_j * t_j;
        let s_accel = 2.0 * s_jerk_up + s_const_a + s_jerk_down;
        let s_cruise = (s - 2.0 * s_accel).max(0.0);
        let t_c = if v_cruise > 1e-15 {
            s_cruise / v_cruise
        } else {
            0.0
        };
        let t_total = 4.0 * t_j + 2.0 * t_a + t_c;
        Self {
            v_max,
            a_max,
            j_max,
            x_start,
            x_target,
            direction,
            t_j,
            t_a,
            t_c,
            t_total,
            v_cruise,
            a_peak,
        }
    }
    /// Query velocity at time `t` (simplified symmetric 7-segment model).
    pub fn velocity(&self, t: f64) -> f64 {
        if t <= 0.0 || t >= self.t_total {
            return 0.0;
        }
        let t_half = self.t_total / 2.0;
        let t_norm = if t <= t_half { t } else { self.t_total - t };
        let v = self.velocity_accel_half(t_norm);
        if t <= t_half {
            self.direction * v
        } else {
            -self.direction * -v
        }
    }
    fn velocity_accel_half(&self, t: f64) -> f64 {
        if t <= self.t_j {
            return 0.5 * self.j_max * t * t;
        }
        let t2 = self.t_j + self.t_a;
        if t <= t2 {
            let v1 = 0.5 * self.j_max * self.t_j * self.t_j;
            return v1 + self.a_max * (t - self.t_j);
        }
        let t3 = t2 + self.t_j;
        if t <= t3 {
            let dt = t - t2;
            let v2 = 0.5 * self.j_max * self.t_j * self.t_j + self.a_max * self.t_a;
            return v2 + self.a_max * dt - 0.5 * self.j_max * dt * dt;
        }
        self.v_cruise
    }
    /// Total profile duration (s).
    pub fn duration(&self) -> f64 {
        self.t_total
    }
}
/// PD position controller for joint motors.
pub struct PositionMotor {
    /// Desired joint position.
    pub target_position: f64,
    /// Current joint position.
    pub current_position: f64,
    /// Proportional gain.
    pub k_p: f64,
    /// Derivative gain.
    pub k_d: f64,
    /// Maximum allowable force magnitude.
    pub max_force: f64,
}
impl PositionMotor {
    /// Creates a new `PositionMotor` with the given gains and force limit.
    pub fn new(k_p: f64, k_d: f64, max_force: f64) -> Self {
        Self {
            target_position: 0.0,
            current_position: 0.0,
            k_p,
            k_d,
            max_force,
        }
    }
    /// Computes the PD control force given the current joint velocity.
    pub fn compute_force(&self, velocity: f64) -> f64 {
        let raw = self.k_p * (self.target_position - self.current_position) - self.k_d * velocity;
        raw.clamp(-self.max_force, self.max_force)
    }
    /// Updates the current position.
    pub fn update_position(&mut self, position: f64) {
        self.current_position = position;
    }
    /// Current position error.
    pub fn position_error(&self) -> f64 {
        self.target_position - self.current_position
    }
    /// Power delivered by the motor at given velocity.
    pub fn power(&self, velocity: f64) -> f64 {
        self.compute_force(velocity) * velocity
    }
}
/// Servo motor combining PD position control with an enable/disable switch.
pub struct ServoMotor {
    /// Underlying position motor.
    pub position_motor: PositionMotor,
    /// Maximum allowed velocity (informational / for external limiting).
    pub max_velocity: f64,
    /// Whether the servo is active.
    pub enabled: bool,
}
impl ServoMotor {
    /// Creates a new `ServoMotor`.
    pub fn new(k_p: f64, k_d: f64, max_force: f64, max_velocity: f64) -> Self {
        Self {
            position_motor: PositionMotor::new(k_p, k_d, max_force),
            max_velocity,
            enabled: true,
        }
    }
    /// Sets the target position.
    pub fn set_target(&mut self, target: f64) {
        self.position_motor.target_position = target;
    }
    /// Computes the PD control force; returns zero if the servo is disabled.
    pub fn compute_force(&self, velocity: f64) -> f64 {
        if self.enabled {
            self.position_motor.compute_force(velocity)
        } else {
            0.0
        }
    }
    /// Enables the servo.
    pub fn enable(&mut self) {
        self.enabled = true;
    }
    /// Disables the servo (output force becomes zero).
    pub fn disable(&mut self) {
        self.enabled = false;
    }
    /// Toggle enabled state.
    pub fn toggle(&mut self) {
        self.enabled = !self.enabled;
    }
}
/// PID auto-tuner that records step-response data and estimates ku, tu via
/// relay feedback.  Uses a simple bang-bang relay of amplitude `relay_amp`.
#[allow(dead_code)]
pub struct RelayAutoTuner {
    /// Relay output amplitude.
    pub relay_amp: f64,
    /// Current relay output sign (+1 or -1).
    pub(super) relay_sign: f64,
    /// Previous process output.
    pub(super) prev_output: f64,
    /// Timestamps of zero crossings (output oscillation).
    pub(super) crossings: Vec<f64>,
    /// Accumulated process output for ultimate gain estimate.
    pub(super) output_sum: f64,
    /// Count of output samples.
    pub(super) sample_count: usize,
    /// Elapsed time.
    pub(super) time: f64,
}
impl RelayAutoTuner {
    /// Create a new relay auto-tuner.
    pub fn new(relay_amp: f64) -> Self {
        Self {
            relay_amp,
            relay_sign: 1.0,
            prev_output: 0.0,
            crossings: Vec::new(),
            output_sum: 0.0,
            sample_count: 0,
            time: 0.0,
        }
    }
    /// Feed the current process output and advance time.
    ///
    /// Returns the relay control signal to apply.
    pub fn update(&mut self, process_output: f64, dt: f64) -> f64 {
        if process_output * self.prev_output < 0.0 {
            self.crossings.push(self.time);
        }
        self.prev_output = process_output;
        self.output_sum += process_output.abs();
        self.sample_count += 1;
        self.time += dt;
        self.relay_sign = if process_output >= 0.0 { -1.0 } else { 1.0 };
        self.relay_sign * self.relay_amp
    }
    /// Estimated oscillation period Tu from detected zero-crossings.
    ///
    /// Returns `None` when fewer than two pairs of crossings are available.
    pub fn estimated_tu(&self) -> Option<f64> {
        if self.crossings.len() < 2 {
            return None;
        }
        let n = self.crossings.len();
        let period = 2.0 * (self.crossings[n - 1] - self.crossings[0]) / (n as f64 - 1.0);
        Some(period)
    }
    /// Estimated ultimate gain Ku = 4*d / (π * a_out).
    ///
    /// `a_out` is the estimated amplitude of the output oscillation.
    pub fn estimated_ku(&self, a_out: f64) -> f64 {
        4.0 * self.relay_amp / (std::f64::consts::PI * a_out.max(1e-30))
    }
    /// Number of oscillation crossings detected.
    pub fn crossing_count(&self) -> usize {
        self.crossings.len()
    }
    /// Reset the tuner state.
    pub fn reset(&mut self) {
        self.relay_sign = 1.0;
        self.prev_output = 0.0;
        self.crossings.clear();
        self.output_sum = 0.0;
        self.sample_count = 0;
        self.time = 0.0;
    }
}
/// Two-phase hybrid stepper motor model.
///
/// Models microstepping, holding torque, detent torque, and full torque–speed
/// characteristic using the simplified 2-phase winding equations.
///
/// # Physical model
///
/// The electromagnetic torque is:
///
/// ```text
/// τ_em = K_t · (iA · sin(N_r · θ) - iB · cos(N_r · θ))
/// ```
///
/// where `N_r` is the number of rotor teeth, `θ` is mechanical angle, `iA` and
/// `iB` are the phase currents.
#[allow(dead_code)]
#[allow(non_snake_case)]
pub struct StepperMotor {
    /// Rotor tooth count (determines step resolution).
    pub N_r: u32,
    /// Phase resistance in Ω.
    pub R_phase: f64,
    /// Phase inductance in H.
    pub L_phase: f64,
    /// Torque constant K_t (N·m/A).
    pub K_t: f64,
    /// Detent torque amplitude (N·m).
    pub T_detent: f64,
    /// Current supply voltage (V) per phase.
    pub V_supply: f64,
    /// Phase A current (A).
    pub iA: f64,
    /// Phase B current (A).
    pub iB: f64,
    /// Mechanical angle (rad).
    pub theta: f64,
    /// Mechanical angular velocity (rad/s).
    pub omega: f64,
    /// Rotor + load inertia (kg·m²).
    pub inertia: f64,
    /// Viscous damping (N·m·s/rad).
    pub damping: f64,
}
impl StepperMotor {
    /// Create a new stepper motor.
    ///
    /// `N_r` is typically 50 for a standard 200-step/rev motor.
    #[allow(clippy::too_many_arguments)]
    #[allow(non_snake_case)]
    pub fn new(
        N_r: u32,
        R_phase: f64,
        L_phase: f64,
        K_t: f64,
        T_detent: f64,
        V_supply: f64,
        inertia: f64,
        damping: f64,
    ) -> Self {
        Self {
            N_r,
            R_phase,
            L_phase,
            K_t,
            T_detent,
            V_supply,
            iA: 0.0,
            iB: 0.0,
            theta: 0.0,
            omega: 0.0,
            inertia,
            damping,
        }
    }
    /// Commanded full-step electrical angle in radians for step `n`.
    ///
    /// For a 2-phase stepper, a full step advances the electrical angle by
    /// π/2.  One mechanical step = π / (2 · N_r).
    #[allow(non_snake_case)]
    pub fn electrical_angle(n: i64, N_r: u32) -> f64 {
        n as f64 * std::f64::consts::FRAC_PI_2 / N_r as f64
    }
    /// Set the phase currents for a target electrical angle `phi_e` with peak
    /// current `I_peak`.
    ///
    /// ```text
    /// iA = I_peak · cos(phi_e)
    /// iB = I_peak · sin(phi_e)
    /// ```
    #[allow(non_snake_case)]
    pub fn set_phase_currents(&mut self, phi_e: f64, I_peak: f64) {
        self.iA = I_peak * phi_e.cos();
        self.iB = I_peak * phi_e.sin();
    }
    /// Electromagnetic torque at the current state.
    ///
    /// ```text
    /// τ_em = K_t · (iA · sin(N_r θ) − iB · cos(N_r θ))
    /// ```
    pub fn electromagnetic_torque(&self) -> f64 {
        let nr_theta = self.N_r as f64 * self.theta;
        self.K_t * (self.iA * nr_theta.sin() - self.iB * nr_theta.cos())
    }
    /// Detent torque at the current angle.
    ///
    /// Detent torque has twice the spatial frequency of the main torque:
    ///
    /// ```text
    /// τ_det = T_detent · sin(2 · N_r · θ)
    /// ```
    pub fn detent_torque(&self) -> f64 {
        (2.0 * self.N_r as f64 * self.theta).sin() * self.T_detent
    }
    /// Net torque (electromagnetic + detent − damping).
    pub fn net_torque(&self) -> f64 {
        self.electromagnetic_torque() + self.detent_torque() - self.damping * self.omega
    }
    /// Holding torque for a commanded electrical angle `phi_e` and peak
    /// current `I_peak` when the rotor is displaced by `delta_theta` from
    /// the equilibrium.
    ///
    /// Linearised about the equilibrium:
    ///
    /// ```text
    /// τ_hold ≈ K_t · I_peak · N_r · sin(N_r · delta_theta)
    /// ```
    #[allow(non_snake_case)]
    pub fn holding_torque(&self, I_peak: f64, delta_theta: f64) -> f64 {
        self.K_t * I_peak * (self.N_r as f64 * delta_theta).sin()
    }
    /// Maximum static holding torque = K_t · I_peak.
    #[allow(non_snake_case)]
    pub fn max_holding_torque(&self, I_peak: f64) -> f64 {
        self.K_t * I_peak
    }
    /// Integrate the stepper dynamics for a time step `dt` using semi-implicit
    /// Euler.  External load torque `tau_load` opposes motion.
    pub fn step(&mut self, tau_load: f64, dt: f64) {
        let tau_net = self.net_torque() - tau_load;
        let alpha = tau_net / self.inertia;
        self.omega += alpha * dt;
        self.theta += self.omega * dt;
    }
    /// Number of full steps per revolution.
    pub fn steps_per_rev(&self) -> u32 {
        4 * self.N_r
    }
    /// Mechanical angle per full step (rad).
    pub fn step_angle_rad(&self) -> f64 {
        2.0 * std::f64::consts::PI / self.steps_per_rev() as f64
    }
}
/// Dynamic current limiter that reduces the current command to protect the
/// motor from thermal overload while tracking a thermal state.
///
/// Uses a simple exponential thermal model:
/// ```text
/// I_max(T) = I_rated * sqrt((T_max - T) / (T_max - T_ambient))
/// ```
#[derive(Debug, Clone)]
pub struct MotorCurrentLimiter {
    /// Rated continuous current (A).
    pub i_rated: f64,
    /// Peak current limit (A) for short bursts.
    pub i_peak: f64,
    /// Maximum winding temperature (°C).
    pub t_max: f64,
    /// Ambient temperature (°C).
    pub t_ambient: f64,
    /// Current winding temperature (°C).
    pub t_winding: f64,
    /// Thermal time constant (s).
    pub tau_thermal: f64,
}
impl MotorCurrentLimiter {
    /// Create a new current limiter.
    pub fn new(i_rated: f64, i_peak: f64, t_max: f64, t_ambient: f64, tau_thermal: f64) -> Self {
        Self {
            i_rated,
            i_peak,
            t_max,
            t_ambient,
            t_winding: t_ambient,
            tau_thermal,
        }
    }
    /// Maximum allowed current given the current thermal state.
    pub fn max_current(&self) -> f64 {
        let margin = self.t_max - self.t_winding;
        let full_margin = self.t_max - self.t_ambient;
        if full_margin < 1e-12 || margin < 0.0 {
            return 0.0;
        }
        self.i_rated * (margin / full_margin).sqrt()
    }
    /// Clamp the requested current `i_req` to the thermal limit.
    pub fn clamp_current(&self, i_req: f64) -> f64 {
        let limit = self.max_current().min(self.i_peak);
        i_req.clamp(-limit, limit)
    }
    /// Integrate winding temperature for one time step.
    ///
    /// Heating: P = I² * R.  Cooling: exponential decay with τ.
    pub fn step(&mut self, i_actual: f64, resistance: f64, dt: f64) {
        let p_heat = i_actual * i_actual * resistance;
        let t_ss = self.t_ambient + p_heat * self.tau_thermal;
        let alpha = 1.0 - (-dt / self.tau_thermal.max(f64::EPSILON)).exp();
        self.t_winding += alpha * (t_ss - self.t_winding);
    }
    /// Whether the motor is currently thermally derated (below rated current).
    pub fn is_derated(&self) -> bool {
        self.max_current() < self.i_rated
    }
    /// Thermal headroom fraction ∈ \[0, 1\].
    pub fn thermal_headroom(&self) -> f64 {
        let margin = self.t_max - self.t_winding;
        let full_margin = self.t_max - self.t_ambient;
        if full_margin < 1e-12 {
            return 0.0;
        }
        (margin / full_margin).clamp(0.0, 1.0)
    }
}
/// A single gear stage with a ratio and efficiency.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct GearStage {
    /// Gear ratio (output_speed / input_speed).  Values < 1 reduce speed.
    pub ratio: f64,
    /// Mechanical efficiency η ∈ (0, 1].
    pub efficiency: f64,
}
impl GearStage {
    /// Create a new gear stage.
    pub fn new(ratio: f64, efficiency: f64) -> Self {
        Self { ratio, efficiency }
    }
    /// Output torque for a given input torque.
    ///
    /// ```text
    /// τ_out = τ_in · η / ratio
    /// ```
    pub fn output_torque(&self, tau_in: f64) -> f64 {
        tau_in * self.efficiency / self.ratio
    }
    /// Output speed for a given input speed.
    pub fn output_speed(&self, omega_in: f64) -> f64 {
        omega_in * self.ratio
    }
    /// Power loss in the stage at input power `P_in`.
    #[allow(non_snake_case)]
    pub fn power_loss(&self, P_in: f64) -> f64 {
        P_in * (1.0 - self.efficiency)
    }
}
/// Chain of gear stages (planetary, spur, worm, etc.) with series-efficiency
/// computation.
///
/// The overall gear ratio is the product of all individual ratios.
/// The overall efficiency is the product of all individual efficiencies.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct GearboxChain {
    /// Ordered list of gear stages from input to output.
    pub stages: Vec<GearStage>,
}
impl GearboxChain {
    /// Create an empty gearbox chain.
    pub fn new() -> Self {
        Self { stages: Vec::new() }
    }
    /// Append a stage.
    pub fn add_stage(&mut self, ratio: f64, efficiency: f64) {
        self.stages.push(GearStage::new(ratio, efficiency));
    }
    /// Overall gear ratio (output/input speed).
    pub fn overall_ratio(&self) -> f64 {
        self.stages.iter().map(|s| s.ratio).product()
    }
    /// Overall mechanical efficiency (product of stage efficiencies).
    pub fn overall_efficiency(&self) -> f64 {
        self.stages.iter().map(|s| s.efficiency).product()
    }
    /// Output torque given input torque `tau_in`.
    pub fn output_torque(&self, tau_in: f64) -> f64 {
        let eta = self.overall_efficiency();
        let ratio = self.overall_ratio();
        tau_in * eta / ratio
    }
    /// Output speed given input speed `omega_in` (rad/s).
    pub fn output_speed(&self, omega_in: f64) -> f64 {
        omega_in * self.overall_ratio()
    }
    /// Total power loss across all stages at input power `P_in` (W).
    #[allow(non_snake_case)]
    pub fn total_power_loss(&self, P_in: f64) -> f64 {
        P_in * (1.0 - self.overall_efficiency())
    }
    /// Back-drive torque: torque required at the output to overcome static
    /// friction and back-drive the input.
    ///
    /// This is the minimum load torque required for back-driveability.
    /// A worm gear with efficiency < 0.5 is typically not back-driveable.
    ///
    /// Formula: `τ_bd = τ_in_static / (η · ratio)` where `η` and `ratio`
    /// are the *reverse* pass values.  For a simple approximation:
    ///
    /// ```text
    /// τ_bd = tau_out / (2η - 1)   when η > 0.5
    /// ```
    pub fn back_drive_torque(&self, tau_out: f64) -> Option<f64> {
        let eta = self.overall_efficiency();
        if eta <= 0.5 {
            None
        } else {
            Some(tau_out / (2.0 * eta - 1.0))
        }
    }
}
/// Phase of a trapezoidal velocity profile.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub enum TrapPhase {
    /// Accelerating toward cruise velocity.
    Accelerating,
    /// Cruising at constant velocity.
    Cruising,
    /// Decelerating toward target.
    Decelerating,
    /// Move complete.
    Done,
}
/// Hydraulic linear actuator driven by fluid pressure.
pub struct HydraulicActuator {
    /// Piston bore area (m2).
    pub bore_area: f64,
    /// Maximum stroke length (m).
    pub stroke: f64,
    /// Maximum operating pressure (Pa).
    pub max_pressure: f64,
    /// Current extension (0..stroke).
    pub position: f64,
}
impl HydraulicActuator {
    /// Create a new hydraulic actuator.
    pub fn new(bore_area: f64, stroke: f64, max_pressure: f64) -> Self {
        Self {
            bore_area,
            stroke,
            max_pressure,
            position: 0.0,
        }
    }
    /// Compute the linear force produced by the given pressure.
    ///
    /// The pressure is clamped to `[0, max_pressure]`.
    pub fn force(&self, pressure: f64) -> f64 {
        let p = pressure.clamp(0.0, self.max_pressure);
        p * self.bore_area
    }
    /// Extend or retract the actuator.  Position is clamped to `[0, stroke]`.
    pub fn set_position(&mut self, pos: f64) {
        self.position = pos.clamp(0.0, self.stroke);
    }
    /// True when fully extended.
    pub fn at_max_extension(&self) -> bool {
        self.position >= self.stroke
    }
    /// True when fully retracted.
    pub fn at_min_extension(&self) -> bool {
        self.position <= 0.0
    }
    /// Volume of fluid displaced at current position (m3).
    pub fn displaced_volume(&self) -> f64 {
        self.bore_area * self.position
    }
    /// Flow rate needed for a desired velocity (m3/s).
    pub fn flow_rate_for_velocity(&self, velocity: f64) -> f64 {
        self.bore_area * velocity.abs()
    }
    /// Compute the volumetric flow rate (m³/s) through the actuator from a
    /// pump/valve model.
    ///
    /// The flow rate is determined by:
    ///
    /// ```text
    /// Q = Cd * A_valve * sqrt(2 * |ΔP| / rho)  ×  sign(ΔP)
    /// ```
    ///
    /// where:
    /// * `Cd`            — discharge coefficient (dimensionless, typically 0.6–0.8).
    /// * `A_valve`       — effective valve opening area (m²).
    /// * `delta_pressure`— pressure differential across the valve (Pa).
    ///   A positive value drives flow into the actuator
    ///   (extension); negative drives retraction.
    /// * `fluid_density` — hydraulic fluid density (kg/m³), e.g. 850 for oil.
    ///
    /// The result is additionally clamped so that the actuator cannot exceed
    /// the maximum flow it can physically absorb (bore area × maximum speed of
    /// `stroke / 0.1 s` as a conservative limit).
    #[allow(clippy::too_many_arguments)]
    pub fn compute_flow_rate(
        &self,
        cd: f64,
        a_valve: f64,
        delta_pressure: f64,
        fluid_density: f64,
    ) -> f64 {
        if fluid_density < 1e-10 || a_valve < 0.0 {
            return 0.0;
        }
        let cd = cd.clamp(0.0, 1.0);
        let a_valve = a_valve.max(0.0);
        let abs_dp = delta_pressure.abs();
        let sign = if delta_pressure >= 0.0 {
            1.0_f64
        } else {
            -1.0_f64
        };
        let q_raw = cd * a_valve * (2.0 * abs_dp / fluid_density).sqrt() * sign;
        let max_piston_vel = self.stroke / 0.1;
        let q_max = self.bore_area * max_piston_vel;
        q_raw.clamp(-q_max, q_max)
    }
    /// Velocity of the piston (m/s) resulting from a given flow rate (m³/s).
    ///
    /// `v = Q / bore_area`.  Returns zero if bore area is negligible.
    pub fn piston_velocity_from_flow(&self, flow_rate: f64) -> f64 {
        if self.bore_area < 1e-30 {
            return 0.0;
        }
        flow_rate / self.bore_area
    }
}
/// Linear actuator that integrates position from a set velocity.
pub struct LinearActuator {
    /// Current position along the actuation axis.
    pub position: f64,
    /// Current velocity along the actuation axis.
    pub velocity: f64,
    /// Minimum allowed position.
    pub min_pos: f64,
    /// Maximum allowed position.
    pub max_pos: f64,
    /// Applied actuation force.
    pub force: f64,
}
impl LinearActuator {
    /// Creates a new `LinearActuator` with the given position limits.
    pub fn new(min_pos: f64, max_pos: f64) -> Self {
        Self {
            position: 0.0,
            velocity: 0.0,
            min_pos,
            max_pos,
            force: 0.0,
        }
    }
    /// Integrates position from velocity over time step `dt`.
    pub fn step(&mut self, dt: f64) {
        self.position += self.velocity * dt;
        self.position = self.position.clamp(self.min_pos, self.max_pos);
    }
    /// Sets the actuation force.
    pub fn set_force(&mut self, f: f64) {
        self.force = f;
    }
    /// Returns `true` if the actuator is at or beyond either position limit.
    pub fn at_limit(&self) -> bool {
        self.position <= self.min_pos || self.position >= self.max_pos
    }
    /// Remaining travel in the positive direction.
    pub fn travel_remaining_positive(&self) -> f64 {
        self.max_pos - self.position
    }
    /// Remaining travel in the negative direction.
    pub fn travel_remaining_negative(&self) -> f64 {
        self.position - self.min_pos
    }
    /// Fraction of total range covered (0..1).
    pub fn range_fraction(&self) -> f64 {
        let range = self.max_pos - self.min_pos;
        if range > 1e-15 {
            (self.position - self.min_pos) / range
        } else {
            0.0
        }
    }
}
/// Decomposed PID contributions for diagnostics.
pub struct PidContributions {
    /// Proportional term.
    pub proportional: f64,
    /// Integral term.
    pub integral: f64,
    /// Derivative term.
    pub derivative: f64,
}
impl PidContributions {
    /// Total output (sum of all contributions).
    pub fn total(&self) -> f64 {
        self.proportional + self.integral + self.derivative
    }
}
