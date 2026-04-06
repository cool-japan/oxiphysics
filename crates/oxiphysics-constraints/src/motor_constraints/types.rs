//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
/// Embedded PID controller used by servo motors.
#[derive(Debug, Clone)]
pub struct MotorPid {
    /// Proportional gain.
    pub kp: f64,
    /// Integral gain.
    pub ki: f64,
    /// Derivative gain.
    pub kd: f64,
    /// Accumulated integral error.
    pub(super) integral: f64,
    /// Previous error (for derivative term).
    pub(super) prev_error: f64,
    /// Anti-windup integral clamp.
    pub integral_clamp: f64,
}
impl MotorPid {
    /// Create a new PID controller with the given gains.
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
    /// Attach anti-windup integral clamping.
    pub fn with_integral_clamp(mut self, clamp: f64) -> Self {
        self.integral_clamp = clamp;
        self
    }
    /// Update the controller and return the control output.
    pub fn update(&mut self, error: f64, dt: f64) -> f64 {
        if dt <= 0.0 {
            return 0.0;
        }
        self.integral =
            (self.integral + error * dt).clamp(-self.integral_clamp, self.integral_clamp);
        let derivative = (error - self.prev_error) / dt;
        self.prev_error = error;
        self.kp * error + self.ki * self.integral + self.kd * derivative
    }
    /// Reset controller state.
    pub fn reset(&mut self) {
        self.integral = 0.0;
        self.prev_error = 0.0;
    }
    /// Get the current integral accumulator.
    pub fn integral_value(&self) -> f64 {
        self.integral
    }
    /// Get the previous error.
    pub fn prev_error_value(&self) -> f64 {
        self.prev_error
    }
}
/// Configurable limits for motor outputs.
#[derive(Debug, Clone)]
pub struct MotorLimits {
    /// Maximum force/torque output (N or N·m).
    pub max_force: f64,
    /// Maximum velocity (m/s or rad/s). `None` = unlimited.
    pub max_velocity: Option<f64>,
    /// Minimum position (m or rad). `None` = unlimited.
    pub min_position: Option<f64>,
    /// Maximum position (m or rad). `None` = unlimited.
    pub max_position: Option<f64>,
    /// Maximum rate of change of force/torque (N/s or N·m/s). `None` = unlimited.
    pub max_force_rate: Option<f64>,
}
impl MotorLimits {
    /// Create limits with only a force limit.
    pub fn force_only(max_force: f64) -> Self {
        Self {
            max_force,
            max_velocity: None,
            min_position: None,
            max_position: None,
            max_force_rate: None,
        }
    }
    /// Create limits with force and velocity limits.
    pub fn force_and_velocity(max_force: f64, max_velocity: f64) -> Self {
        Self {
            max_force,
            max_velocity: Some(max_velocity),
            min_position: None,
            max_position: None,
            max_force_rate: None,
        }
    }
    /// Create limits with force and position range.
    pub fn with_position_range(mut self, min_pos: f64, max_pos: f64) -> Self {
        self.min_position = Some(min_pos);
        self.max_position = Some(max_pos);
        self
    }
    /// Clamp a force value to the limits.
    pub fn clamp_force(&self, force: f64) -> f64 {
        force.clamp(-self.max_force, self.max_force)
    }
    /// Clamp a velocity to the velocity limit if set.
    pub fn clamp_velocity(&self, velocity: f64) -> f64 {
        if let Some(max_v) = self.max_velocity {
            velocity.clamp(-max_v, max_v)
        } else {
            velocity
        }
    }
    /// Check if a position is within the allowed range.
    pub fn position_in_range(&self, position: f64) -> bool {
        if let Some(min_p) = self.min_position
            && position < min_p
        {
            return false;
        }
        if let Some(max_p) = self.max_position
            && position > max_p
        {
            return false;
        }
        true
    }
    /// Apply rate limiting to force.
    pub fn rate_limit_force(&self, current: f64, previous: f64, dt: f64) -> f64 {
        if let Some(max_rate) = self.max_force_rate
            && dt > 1e-15
        {
            let rate = (current - previous) / dt;
            if rate.abs() > max_rate {
                return previous + max_rate * rate.signum() * dt;
            }
        }
        current
    }
}
/// Multi-axis motor coordination controller.
///
/// Coordinates N motors to track a joint trajectory while enforcing
/// per-axis force/torque limits and synchronization.
#[derive(Debug, Clone)]
pub struct MultiAxisCoordinator {
    /// Number of axes.
    pub n_axes: usize,
    /// Per-axis maximum force/torque limits.
    pub max_forces: Vec<f64>,
    /// Per-axis PID controllers.
    pub pids: Vec<MotorPid>,
    /// Per-axis target positions.
    pub targets: Vec<f64>,
    /// Per-axis coupling gains (cross-coupling compensation).
    pub coupling_gains: Vec<Vec<f64>>,
}
impl MultiAxisCoordinator {
    /// Create a new multi-axis coordinator with decoupled P controllers.
    pub fn new(n_axes: usize, kp: f64, max_force: f64) -> Self {
        Self {
            n_axes,
            max_forces: vec![max_force; n_axes],
            pids: (0..n_axes).map(|_| MotorPid::new(kp, 0.0, 0.0)).collect(),
            targets: vec![0.0; n_axes],
            coupling_gains: vec![vec![0.0; n_axes]; n_axes],
        }
    }
    /// Set the target positions for all axes.
    pub fn set_targets(&mut self, targets: &[f64]) {
        for (i, &t) in targets.iter().enumerate().take(self.n_axes) {
            self.targets[i] = t;
        }
    }
    /// Set a cross-coupling compensation gain between axes i and j.
    pub fn set_coupling(&mut self, i: usize, j: usize, gain: f64) {
        if i < self.n_axes && j < self.n_axes {
            self.coupling_gains[i][j] = gain;
        }
    }
    /// Compute output forces/torques for all axes given current positions.
    ///
    /// Returns a Vec of forces/torques, one per axis.
    #[allow(clippy::too_many_arguments)]
    pub fn step(&mut self, dt: f64, positions: &[f64]) -> Vec<f64> {
        let n = self.n_axes;
        let mut errors = vec![0.0f64; n];
        for i in 0..n {
            let pos = positions.get(i).copied().unwrap_or(0.0);
            errors[i] = self.targets[i] - pos;
        }
        let mut forces = vec![0.0f64; n];
        for i in 0..n {
            let base = self.pids[i].update(errors[i], dt);
            let mut coupling = 0.0;
            for j in 0..n {
                if j != i {
                    coupling += self.coupling_gains[i][j] * errors[j];
                }
            }
            forces[i] = (base + coupling).clamp(-self.max_forces[i], self.max_forces[i]);
        }
        forces
    }
    /// Synchronization error: max difference between axis position errors.
    pub fn sync_error(&self, positions: &[f64]) -> f64 {
        let errors: Vec<f64> = (0..self.n_axes)
            .map(|i| {
                let pos = positions.get(i).copied().unwrap_or(0.0);
                (self.targets[i] - pos).abs()
            })
            .collect();
        let max_e = errors.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let min_e = errors.iter().copied().fold(f64::INFINITY, f64::min);
        max_e - min_e
    }
    /// Reset all PID controllers.
    pub fn reset_all(&mut self) {
        for pid in &mut self.pids {
            pid.reset();
        }
    }
}
/// Active compliance controller.
///
/// Implements virtual spring-damper behavior in software: the controller
/// modifies the motor command to achieve a desired compliance profile.
///
/// F_cmd = K_v * (x_d - x) + D_v * (dx_d/dt - dx/dt) + F_ext_comp
#[derive(Debug, Clone)]
pub struct ActiveCompliance {
    /// Virtual stiffness (N/m or N·m/rad).
    pub virtual_stiffness: f64,
    /// Virtual damping (N·s/m or N·m·s/rad).
    pub virtual_damping: f64,
    /// Maximum compliance force/torque (N or N·m).
    pub max_compliance_force: f64,
    /// Desired equilibrium position.
    pub equilibrium_position: f64,
    /// Desired equilibrium velocity.
    pub equilibrium_velocity: f64,
}
impl ActiveCompliance {
    /// Create a new active compliance controller.
    pub fn new(virtual_stiffness: f64, virtual_damping: f64, max_force: f64) -> Self {
        Self {
            virtual_stiffness,
            virtual_damping,
            max_compliance_force: max_force,
            equilibrium_position: 0.0,
            equilibrium_velocity: 0.0,
        }
    }
    /// Set the equilibrium setpoint.
    pub fn set_equilibrium(&mut self, position: f64, velocity: f64) {
        self.equilibrium_position = position;
        self.equilibrium_velocity = velocity;
    }
    /// Compute the compliance force for current state.
    pub fn compliance_force(&self, position: f64, velocity: f64) -> f64 {
        let spring = self.virtual_stiffness * (self.equilibrium_position - position);
        let damper = self.virtual_damping * (self.equilibrium_velocity - velocity);
        (spring + damper).clamp(-self.max_compliance_force, self.max_compliance_force)
    }
    /// Effective stiffness at a given frequency ω (rad/s).
    ///
    /// Impedance: Z(jω) = K + jωD  → magnitude |Z| = sqrt(K² + (ωD)²)
    pub fn impedance_magnitude(&self, omega: f64) -> f64 {
        let k = self.virtual_stiffness;
        let d = self.virtual_damping * omega;
        (k * k + d * d).sqrt()
    }
    /// Damping ratio: ζ = D / (2 * sqrt(K * M)).
    pub fn damping_ratio(&self, mass: f64) -> f64 {
        if mass < 1e-20 || self.virtual_stiffness < 1e-20 {
            return 0.0;
        }
        self.virtual_damping / (2.0 * (self.virtual_stiffness * mass).sqrt())
    }
}
/// Worm gear model.
///
/// A worm gear converts high-speed rotary motion to low-speed high-torque
/// output with inherent self-locking at low lead angles.
#[derive(Debug, Clone)]
pub struct WormGear {
    /// Number of worm starts (threads).
    pub starts: u32,
    /// Number of worm wheel teeth.
    pub wheel_teeth: u32,
    /// Lead angle (radians) — determines self-locking behavior.
    pub lead_angle: f64,
    /// Thread friction coefficient (μ).
    pub friction_coeff: f64,
    /// Normal force on the tooth (N) — for efficiency calculation.
    pub normal_pressure_angle: f64,
}
impl WormGear {
    /// Create a new worm gear.
    pub fn new(starts: u32, wheel_teeth: u32, lead_angle_deg: f64, friction_coeff: f64) -> Self {
        Self {
            starts,
            wheel_teeth,
            lead_angle: lead_angle_deg.to_radians(),
            friction_coeff,
            normal_pressure_angle: 20.0_f64.to_radians(),
        }
    }
    /// Gear ratio: wheel_teeth / starts.
    pub fn gear_ratio(&self) -> f64 {
        self.wheel_teeth as f64 / self.starts as f64
    }
    /// Driving efficiency (worm → wheel): η_d = tan(λ) / tan(λ + φ)
    /// where λ = lead angle, φ = friction angle = atan(μ/cos(αn)).
    pub fn driving_efficiency(&self) -> f64 {
        let phi = (self.friction_coeff / self.normal_pressure_angle.cos()).atan();
        let num = self.lead_angle.tan();
        let den = (self.lead_angle + phi).tan();
        if den.abs() < 1e-20 {
            return 0.0;
        }
        (num / den).clamp(0.0, 1.0)
    }
    /// Back-driving efficiency (wheel → worm): η_b = tan(λ - φ) / tan(λ)
    pub fn back_driving_efficiency(&self) -> f64 {
        let phi = (self.friction_coeff / self.normal_pressure_angle.cos()).atan();
        if self.lead_angle <= phi {
            return 0.0;
        }
        let num = (self.lead_angle - phi).tan();
        let den = self.lead_angle.tan();
        if den.abs() < 1e-20 {
            return 0.0;
        }
        (num / den).clamp(0.0, 1.0)
    }
    /// Whether the gear is self-locking (back-drive efficiency ≤ 0).
    pub fn is_self_locking(&self) -> bool {
        let phi = (self.friction_coeff / self.normal_pressure_angle.cos()).atan();
        self.lead_angle <= phi
    }
    /// Output torque given input torque (worm → wheel), accounting for efficiency.
    pub fn output_torque(&self, input_torque: f64) -> f64 {
        input_torque * self.gear_ratio() * self.driving_efficiency()
    }
}
/// Motor thermal model for heat dissipation.
///
/// Uses a lumped-capacitance (RC) thermal model:
/// `C * dT/dt = P_loss - (T - T_ambient) / R_th`
///
/// where:
/// - `C` is the thermal capacitance (J/K)
/// - `R_th` is the thermal resistance (K/W)
/// - `P_loss` is the instantaneous power loss (W)
#[derive(Debug, Clone)]
pub struct MotorThermal {
    /// Thermal capacitance of the motor winding (J/K).
    pub thermal_capacitance: f64,
    /// Thermal resistance winding-to-ambient (K/W).
    pub thermal_resistance: f64,
    /// Ambient temperature (K or °C).
    pub ambient_temperature: f64,
    /// Current winding temperature (K or °C).
    pub temperature: f64,
    /// Maximum permissible temperature (K or °C).
    pub max_temperature: f64,
    /// Motor electrical resistance (Ω) for copper loss computation.
    pub resistance_ohm: f64,
}
impl MotorThermal {
    /// Create a new thermal model.
    pub fn new(
        thermal_capacitance: f64,
        thermal_resistance: f64,
        ambient_temperature: f64,
        resistance_ohm: f64,
        max_temperature: f64,
    ) -> Self {
        Self {
            thermal_capacitance,
            thermal_resistance,
            ambient_temperature,
            temperature: ambient_temperature,
            max_temperature,
            resistance_ohm,
        }
    }
    /// Compute copper (I²R) power loss.
    pub fn copper_loss(&self, current_amps: f64) -> f64 {
        current_amps * current_amps * self.resistance_ohm
    }
    /// Compute iron (core) loss approximation: `P_iron = k_e * omega^2 + k_h * omega`.
    pub fn iron_loss(&self, angular_velocity: f64, ke: f64, kh: f64) -> f64 {
        let omega = angular_velocity.abs();
        ke * omega * omega + kh * omega
    }
    /// Total power loss (copper + iron).
    pub fn total_loss(&self, current_amps: f64, angular_velocity: f64, ke: f64, kh: f64) -> f64 {
        self.copper_loss(current_amps) + self.iron_loss(angular_velocity, ke, kh)
    }
    /// Advance the thermal model by one step `dt` (s) with given power loss.
    pub fn step(&mut self, p_loss: f64, dt: f64) {
        if self.thermal_capacitance < 1e-30 {
            return;
        }
        let dt_cool = (self.temperature - self.ambient_temperature) / self.thermal_resistance;
        let delta_t = (p_loss - dt_cool) / self.thermal_capacitance * dt;
        self.temperature += delta_t;
    }
    /// Returns true if the motor is overheated.
    pub fn is_overheated(&self) -> bool {
        self.temperature > self.max_temperature
    }
    /// Steady-state temperature for a constant power loss: `T_ss = T_amb + P * R_th`.
    pub fn steady_state_temperature(&self, p_loss: f64) -> f64 {
        self.ambient_temperature + p_loss * self.thermal_resistance
    }
    /// Thermal time constant: `tau = C * R_th`.
    pub fn time_constant(&self) -> f64 {
        self.thermal_capacitance * self.thermal_resistance
    }
    /// Reset temperature to ambient.
    pub fn reset(&mut self) {
        self.temperature = self.ambient_temperature;
    }
}
/// Series Elastic Actuator (SEA) model.
///
/// An SEA places a compliant spring between the motor and load.  This enables
/// accurate force control and provides passive compliance for safety.
///
/// The spring deflection `x_s = x_motor - x_load` gives the output force:
/// `F = k_s * x_s`.
#[derive(Debug, Clone)]
pub struct SeriesElasticActuator {
    /// Series spring stiffness (N/m or N·m/rad).
    pub spring_stiffness: f64,
    /// Motor-side position (m or rad).
    pub motor_position: f64,
    /// Load-side position (m or rad).
    pub load_position: f64,
    /// Target output force/torque (N or N·m).
    pub target_force: f64,
    /// Force controller gain (motor position rate per force error, 1/N).
    pub force_controller_gain: f64,
}
impl SeriesElasticActuator {
    /// Create a new SEA.
    pub fn new(spring_stiffness: f64, force_controller_gain: f64) -> Self {
        Self {
            spring_stiffness,
            motor_position: 0.0,
            load_position: 0.0,
            target_force: 0.0,
            force_controller_gain,
        }
    }
    /// Spring deflection.
    pub fn deflection(&self) -> f64 {
        self.motor_position - self.load_position
    }
    /// Current output force/torque (N or N·m).
    pub fn output_force(&self) -> f64 {
        self.spring_stiffness * self.deflection()
    }
    /// Force error: target - actual.
    pub fn force_error(&self) -> f64 {
        self.target_force - self.output_force()
    }
    /// Desired motor velocity for force control (m/s or rad/s).
    pub fn desired_motor_velocity(&self) -> f64 {
        self.force_controller_gain * self.force_error()
    }
    /// Step the SEA: advance motor and load positions.
    ///
    /// * `motor_velocity` — motor velocity command (m/s).
    /// * `load_velocity`  — load velocity (from dynamics, m/s).
    pub fn step(&mut self, motor_velocity: f64, load_velocity: f64, dt: f64) {
        self.motor_position += motor_velocity * dt;
        self.load_position += load_velocity * dt;
    }
    /// Elastic energy stored in the series spring (J).
    pub fn elastic_energy(&self) -> f64 {
        0.5 * self.spring_stiffness * self.deflection().powi(2)
    }
}
/// A constraint that drives a rotational DOF toward a target angular velocity.
#[derive(Debug, Clone)]
pub struct AngularMotorConstraint {
    /// Desired angular velocity (rad/s).
    pub target_angular_velocity: f64,
    /// Maximum torque this motor can exert (N·m).
    pub max_torque: f64,
    /// Gain for angular velocity error.
    pub velocity_gain: f64,
}
impl AngularMotorConstraint {
    /// Create a new angular motor constraint.
    pub fn new(target_angular_velocity: f64, max_torque: f64) -> Self {
        Self {
            target_angular_velocity,
            max_torque,
            velocity_gain: 100.0,
        }
    }
    /// Set the velocity error gain.
    pub fn with_gain(mut self, gain: f64) -> Self {
        self.velocity_gain = gain;
        self
    }
    /// Compute the motor torque for the current time step.
    pub fn step(&mut self, dt: f64, actual_vel: f64, actual_pos: f64) -> f64 {
        let _ = (dt, actual_pos);
        let error = self.target_angular_velocity - actual_vel;
        let raw = self.velocity_gain * error;
        raw.clamp(-self.max_torque, self.max_torque)
    }
    /// Set a new target angular velocity.
    pub fn set_target(&mut self, target: f64) {
        self.target_angular_velocity = target;
    }
}
/// Cable/pulley transmission system.
///
/// A single fixed pulley of radius `r` converts drum rotation to linear
/// cable travel.  Cable stiffness introduces compliance.
#[derive(Debug, Clone)]
pub struct CablePulleySystem {
    /// Pulley radius (m).
    pub radius: f64,
    /// Cable stiffness (N/m).
    pub cable_stiffness: f64,
    /// Cable natural length (m).
    pub natural_length: f64,
    /// Pulley mechanical efficiency (0..1).
    pub efficiency: f64,
    /// Current cable elongation (m) — positive = stretched.
    pub elongation: f64,
}
impl CablePulleySystem {
    /// Create a cable/pulley system.
    pub fn new(radius: f64, cable_stiffness: f64, natural_length: f64) -> Self {
        Self {
            radius,
            cable_stiffness,
            natural_length,
            efficiency: 1.0,
            elongation: 0.0,
        }
    }
    /// Cable tension from elongation (N).
    pub fn tension(&self) -> f64 {
        (self.cable_stiffness * self.elongation).max(0.0)
    }
    /// Load torque on the drum for current cable tension.
    pub fn drum_torque(&self) -> f64 {
        self.tension() * self.radius / self.efficiency.max(1e-12)
    }
    /// Update elongation given drum angle and load displacement.
    ///
    /// * `drum_angle` — accumulated drum rotation (rad).
    /// * `load_displacement` — load end displacement (m, positive = pulling cable).
    pub fn update_elongation(&mut self, drum_angle: f64, load_displacement: f64) {
        let cable_paid_out = drum_angle * self.radius;
        self.elongation = (load_displacement - cable_paid_out).max(0.0);
    }
    /// Kinetic energy stored in cable elasticity (J).
    pub fn elastic_energy(&self) -> f64 {
        0.5 * self.cable_stiffness * self.elongation * self.elongation
    }
}
/// Batch solver for a collection of motor constraints.
#[derive(Debug, Clone, Default)]
pub struct MotorSolver {
    /// The motor constraints managed by this solver.
    pub motors: Vec<MotorConstraint>,
}
impl MotorSolver {
    /// Create an empty motor solver.
    pub fn new() -> Self {
        Self { motors: Vec::new() }
    }
    /// Add a motor constraint to the solver.
    pub fn add(&mut self, motor: MotorConstraint) -> usize {
        let idx = self.motors.len();
        self.motors.push(motor);
        idx
    }
    /// Remove a motor by index. Returns `None` if out of bounds.
    pub fn remove(&mut self, index: usize) -> Option<MotorConstraint> {
        if index < self.motors.len() {
            Some(self.motors.remove(index))
        } else {
            None
        }
    }
    /// Number of motors.
    pub fn len(&self) -> usize {
        self.motors.len()
    }
    /// Whether the solver has no motors.
    pub fn is_empty(&self) -> bool {
        self.motors.is_empty()
    }
    /// Advance all motors by `dt` and return a vector of output forces/torques.
    pub fn step(&mut self, dt: f64, actual_vel: &[f64], actual_pos: &[f64]) -> Vec<f64> {
        self.motors
            .iter_mut()
            .enumerate()
            .map(|(i, m)| {
                let vel = actual_vel.get(i).copied().unwrap_or(0.0);
                let pos = actual_pos.get(i).copied().unwrap_or(0.0);
                m.step(dt, vel, pos)
            })
            .collect()
    }
    /// Compute total energy expenditure (sum of |force * velocity| * dt).
    pub fn energy_expenditure(&self, forces: &[f64], velocities: &[f64], dt: f64) -> f64 {
        forces
            .iter()
            .zip(velocities.iter())
            .map(|(f, v)| (f * v).abs() * dt)
            .sum()
    }
}
/// Regenerative braking model.
///
/// During braking, kinetic energy is converted to electrical energy
/// (stored in a battery or capacitor) instead of being wasted as heat.
///
/// The regenerative torque is:
/// `T_regen = -min(eta_regen * T_brake_requested, T_regen_max)`
#[derive(Debug, Clone)]
pub struct RegenerativeBraking {
    /// Regenerative efficiency (0..1).
    pub efficiency: f64,
    /// Maximum regenerative torque (N·m).
    pub max_regen_torque: f64,
    /// Minimum speed below which regen is disabled (rad/s).
    pub min_speed: f64,
    /// Accumulated regenerated energy (J).
    pub energy_recovered: f64,
}
impl RegenerativeBraking {
    /// Create a new regenerative braking model.
    pub fn new(efficiency: f64, max_regen_torque: f64, min_speed: f64) -> Self {
        Self {
            efficiency,
            max_regen_torque,
            min_speed,
            energy_recovered: 0.0,
        }
    }
    /// Compute the regenerative braking torque and update recovered energy.
    ///
    /// Returns the effective braking torque (negative = opposing motion).
    pub fn compute_regen_torque(
        &mut self,
        requested_brake_torque: f64,
        angular_velocity: f64,
        dt: f64,
    ) -> f64 {
        let speed = angular_velocity.abs();
        if speed < self.min_speed {
            return 0.0;
        }
        let regen = (self.efficiency * requested_brake_torque.abs()).min(self.max_regen_torque);
        let torque = -regen * angular_velocity.signum();
        let power = regen * speed;
        self.energy_recovered += power * dt;
        torque
    }
    /// Reset recovered energy counter.
    pub fn reset_energy(&mut self) {
        self.energy_recovered = 0.0;
    }
    /// Instantaneous regenerative power at given speed and torque.
    pub fn regen_power(&self, torque: f64, angular_velocity: f64) -> f64 {
        self.efficiency * torque.abs() * angular_velocity.abs()
    }
}
/// A piecewise-linear drive cycle for motor testing.
///
/// Stores `(time, target_velocity)` setpoints. The velocity is linearly
/// interpolated between adjacent setpoints.
#[derive(Debug, Clone)]
pub struct DriveCycle {
    /// List of `(t, v_target)` setpoints in ascending time order.
    pub setpoints: Vec<(f64, f64)>,
}
impl DriveCycle {
    /// Create a new empty drive cycle.
    pub fn new() -> Self {
        Self {
            setpoints: Vec::new(),
        }
    }
    /// Add a setpoint `(time, velocity)`.
    pub fn add_setpoint(&mut self, time: f64, velocity: f64) {
        self.setpoints.push((time, velocity));
        self.setpoints
            .sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    }
    /// Get the target velocity at time `t` by linear interpolation.
    pub fn velocity_at(&self, t: f64) -> f64 {
        if self.setpoints.is_empty() {
            return 0.0;
        }
        if t <= self.setpoints[0].0 {
            return self.setpoints[0].1;
        }
        if t >= self
            .setpoints
            .last()
            .expect("collection should not be empty")
            .0
        {
            return self
                .setpoints
                .last()
                .expect("collection should not be empty")
                .1;
        }
        for i in 1..self.setpoints.len() {
            let (t0, v0) = self.setpoints[i - 1];
            let (t1, v1) = self.setpoints[i];
            if t <= t1 {
                let alpha = (t - t0) / (t1 - t0);
                return v0 + alpha * (v1 - v0);
            }
        }
        self.setpoints
            .last()
            .expect("collection should not be empty")
            .1
    }
    /// Duration of the drive cycle (time of last setpoint).
    pub fn duration(&self) -> f64 {
        self.setpoints.last().map(|s| s.0).unwrap_or(0.0)
    }
    /// Peak (maximum absolute) velocity in the cycle.
    pub fn peak_velocity(&self) -> f64 {
        self.setpoints
            .iter()
            .map(|s| s.1.abs())
            .fold(0.0_f64, f64::max)
    }
    /// Root-mean-square velocity (trapezoidal integration).
    pub fn rms_velocity(&self, dt: f64) -> f64 {
        if self.setpoints.is_empty() {
            return 0.0;
        }
        let duration = self.duration();
        if duration < 1e-30 {
            return 0.0;
        }
        let steps = (duration / dt).ceil() as usize + 1;
        let mut sum_sq = 0.0;
        let mut count = 0;
        for k in 0..steps {
            let t = k as f64 * dt;
            let v = self.velocity_at(t);
            sum_sq += v * v;
            count += 1;
        }
        (sum_sq / count as f64).sqrt()
    }
}
/// A constraint that drives one body along a linear axis toward a target velocity
/// and/or target position.
#[derive(Debug, Clone)]
pub struct LinearMotorConstraint {
    /// Desired velocity along the motor axis (m/s).  `None` disables velocity mode.
    pub target_velocity: Option<f64>,
    /// Desired position along the motor axis (m).  `None` disables position mode.
    pub target_position: Option<f64>,
    /// Maximum force this motor can exert (N).
    pub max_force: f64,
    /// Current position along the motor axis (m).
    pub current_position: f64,
    /// Gain for velocity error in velocity-mode.
    pub velocity_gain: f64,
    /// Gain for position error when `target_position` is set.
    pub position_gain: f64,
}
impl LinearMotorConstraint {
    /// Create a velocity-only linear motor.
    pub fn velocity_mode(target_velocity: f64, max_force: f64) -> Self {
        Self {
            target_velocity: Some(target_velocity),
            target_position: None,
            max_force,
            current_position: 0.0,
            velocity_gain: 100.0,
            position_gain: 0.0,
        }
    }
    /// Create a position-only linear motor.
    pub fn position_mode(target_position: f64, max_force: f64) -> Self {
        Self {
            target_velocity: None,
            target_position: Some(target_position),
            max_force,
            current_position: 0.0,
            velocity_gain: 0.0,
            position_gain: 200.0,
        }
    }
    /// Create a combined position+velocity motor.
    pub fn combined_mode(target_position: f64, target_velocity: f64, max_force: f64) -> Self {
        Self {
            target_velocity: Some(target_velocity),
            target_position: Some(target_position),
            max_force,
            current_position: 0.0,
            velocity_gain: 50.0,
            position_gain: 200.0,
        }
    }
    /// Compute the motor force for the current time step.
    pub fn step(&mut self, dt: f64, actual_vel: f64, actual_pos: f64) -> f64 {
        self.current_position = actual_pos;
        let _ = dt;
        let mut force = 0.0;
        if let Some(tv) = self.target_velocity {
            let vel_error = tv - actual_vel;
            force += self.velocity_gain * vel_error;
        }
        if let Some(tp) = self.target_position {
            let pos_error = tp - actual_pos;
            force += self.position_gain * pos_error;
        }
        force.clamp(-self.max_force, self.max_force)
    }
    /// Compute force with friction compensation.
    pub fn step_with_friction(
        &mut self,
        dt: f64,
        actual_vel: f64,
        actual_pos: f64,
        friction: &MotorFriction,
    ) -> f64 {
        let raw = self.step(dt, actual_vel, actual_pos);
        let friction_comp = friction.friction_force(actual_vel);
        (raw - friction_comp).clamp(-self.max_force, self.max_force)
    }
}
/// Impedance control motor that regulates the mechanical impedance
/// (stiffness + damping) between the motor and its environment.
///
/// F = K * (x_d - x) + D * (v_d - v)
#[derive(Debug, Clone)]
pub struct ImpedanceControlMotor {
    /// Desired position (m or rad).
    pub target_position: f64,
    /// Desired velocity (m/s or rad/s).
    pub target_velocity: f64,
    /// Virtual stiffness K (N/m or N·m/rad).
    pub stiffness: f64,
    /// Virtual damping D (N·s/m or N·m·s/rad).
    pub damping: f64,
    /// Maximum output force/torque.
    pub max_force: f64,
}
impl ImpedanceControlMotor {
    /// Create a new impedance control motor.
    pub fn new(stiffness: f64, damping: f64, max_force: f64) -> Self {
        Self {
            target_position: 0.0,
            target_velocity: 0.0,
            stiffness,
            damping,
            max_force,
        }
    }
    /// Set desired trajectory.
    pub fn set_desired(&mut self, position: f64, velocity: f64) {
        self.target_position = position;
        self.target_velocity = velocity;
    }
    /// Compute the impedance control output force.
    pub fn step(&self, actual_vel: f64, actual_pos: f64) -> f64 {
        let pos_error = self.target_position - actual_pos;
        let vel_error = self.target_velocity - actual_vel;
        let force = self.stiffness * pos_error + self.damping * vel_error;
        force.clamp(-self.max_force, self.max_force)
    }
    /// Compute the mechanical impedance at a given frequency.
    /// Z(jω) = K/(jω) + D   →  |Z| = sqrt(K²/ω² + D²)
    pub fn impedance_magnitude(&self, omega: f64) -> f64 {
        if omega.abs() < 1e-15 {
            return f64::INFINITY;
        }
        ((self.stiffness / omega).powi(2) + self.damping * self.damping).sqrt()
    }
}
/// A position-servo motor using an embedded PID controller.
#[derive(Debug, Clone)]
pub struct ServoMotorConstraint {
    /// Desired position (m or rad depending on DOF type).
    pub target_position: f64,
    /// Maximum output torque/force (N·m or N).
    pub max_torque: f64,
    /// Embedded PID controller.
    pub pid: MotorPid,
}
impl ServoMotorConstraint {
    /// Create a new servo motor with explicit PID gains.
    pub fn new(target_position: f64, max_torque: f64, kp: f64, ki: f64, kd: f64) -> Self {
        Self {
            target_position,
            max_torque,
            pid: MotorPid::new(kp, ki, kd),
        }
    }
    /// Create a servo with a typical stiff PD configuration.
    pub fn with_pd_gains(target_position: f64, max_torque: f64, kp: f64, kd: f64) -> Self {
        Self::new(target_position, max_torque, kp, 0.0, kd)
    }
    /// Set the target position.
    pub fn set_target(&mut self, target: f64) {
        self.target_position = target;
    }
    /// Compute the servo output torque/force for the current time step.
    pub fn step(&mut self, dt: f64, actual_vel: f64, actual_pos: f64) -> f64 {
        let _ = actual_vel;
        let error = self.target_position - actual_pos;
        let raw = self.pid.update(error, dt);
        raw.clamp(-self.max_torque, self.max_torque)
    }
    /// Reset the embedded PID state.
    pub fn reset(&mut self) {
        self.pid.reset();
    }
}
/// Ball screw mechanism model.
///
/// Converts rotary motion to linear motion with high efficiency and
/// low backlash.  Used in CNC machines and robotic actuators.
#[derive(Debug, Clone)]
pub struct BallScrew {
    /// Lead (linear travel per revolution, meters/rev).
    pub lead: f64,
    /// Mechanical efficiency (0..1, typically 0.9–0.95).
    pub efficiency: f64,
    /// Maximum axial force (N).
    pub max_axial_force: f64,
    /// Current position (m).
    pub position: f64,
    /// Current velocity (m/s).
    pub velocity: f64,
}
impl BallScrew {
    /// Create a new ball screw.
    pub fn new(lead: f64, efficiency: f64, max_axial_force: f64) -> Self {
        Self {
            lead,
            efficiency,
            max_axial_force,
            position: 0.0,
            velocity: 0.0,
        }
    }
    /// Axial force from applied torque (N·m).
    pub fn axial_force(&self, torque: f64) -> f64 {
        let force = torque * 2.0 * std::f64::consts::PI * self.efficiency / self.lead.max(1e-12);
        force.clamp(-self.max_axial_force, self.max_axial_force)
    }
    /// Required torque for a given axial force (N·m).
    pub fn required_torque(&self, axial_force: f64) -> f64 {
        axial_force * self.lead / (2.0 * std::f64::consts::PI * self.efficiency.max(1e-12))
    }
    /// Linear velocity (m/s) from angular velocity (rad/s).
    pub fn linear_velocity(&self, angular_velocity: f64) -> f64 {
        angular_velocity * self.lead / (2.0 * std::f64::consts::PI)
    }
    /// Angular velocity (rad/s) from linear velocity (m/s).
    pub fn angular_velocity_from_linear(&self, linear_velocity: f64) -> f64 {
        linear_velocity * 2.0 * std::f64::consts::PI / self.lead.max(1e-12)
    }
    /// Advance the screw position.
    pub fn step(&mut self, angular_velocity: f64, dt: f64) {
        self.velocity = self.linear_velocity(angular_velocity);
        self.position += self.velocity * dt;
    }
}
/// Selects which type of motor constraint to use.
#[derive(Debug, Clone)]
pub enum MotorConstraint {
    /// Linear velocity/position motor.
    Linear(LinearMotorConstraint),
    /// Angular velocity motor.
    Angular(AngularMotorConstraint),
    /// Position servo with PID controller.
    Servo(ServoMotorConstraint),
    /// Impedance control motor.
    Impedance(ImpedanceControlMotor),
}
impl MotorConstraint {
    /// Advance the motor by one time step and return the output force/torque.
    pub fn step(&mut self, dt: f64, actual_vel: f64, actual_pos: f64) -> f64 {
        match self {
            MotorConstraint::Linear(m) => m.step(dt, actual_vel, actual_pos),
            MotorConstraint::Angular(m) => m.step(dt, actual_vel, actual_pos),
            MotorConstraint::Servo(m) => m.step(dt, actual_vel, actual_pos),
            MotorConstraint::Impedance(m) => m.step(actual_vel, actual_pos),
        }
    }
}
/// Harmonic drive (strain wave gear) model.
///
/// A harmonic drive achieves very high gear ratios in a compact package by
/// using the elastic deformation of a flexspline.  It has near-zero backlash
/// and is used extensively in robotics.
///
/// The ratio is: `output_angle = input_angle / gear_ratio`.
#[derive(Debug, Clone)]
pub struct HarmonicDrive {
    /// Gear ratio (input / output rotations, typically 50–320).
    pub gear_ratio: f64,
    /// Torsional stiffness of the flexspline (N·m/rad).
    pub flexspline_stiffness: f64,
    /// Torsional damping (N·m·s/rad).
    pub flexspline_damping: f64,
    /// Coulomb friction torque (N·m), referred to output shaft.
    pub friction_torque: f64,
    /// Current output angle (rad).
    pub output_angle: f64,
    /// Current output angular velocity (rad/s).
    pub output_velocity: f64,
    /// Output shaft inertia (kg·m²).
    pub output_inertia: f64,
}
impl HarmonicDrive {
    /// Create a new harmonic drive.
    pub fn new(gear_ratio: f64, flexspline_stiffness: f64, output_inertia: f64) -> Self {
        Self {
            gear_ratio,
            flexspline_stiffness,
            flexspline_damping: 0.0,
            friction_torque: 0.0,
            output_angle: 0.0,
            output_velocity: 0.0,
            output_inertia,
        }
    }
    /// Set torsional damping.
    pub fn with_damping(mut self, damping: f64) -> Self {
        self.flexspline_damping = damping;
        self
    }
    /// Set friction torque.
    pub fn with_friction(mut self, friction: f64) -> Self {
        self.friction_torque = friction;
        self
    }
    /// Compute the output torque for a given input torque.
    ///
    /// Torque is amplified by the gear ratio minus losses.
    pub fn output_torque(&self, input_torque: f64) -> f64 {
        let amplified = input_torque * self.gear_ratio;
        let friction = if self.output_velocity.abs() < 1e-6 {
            0.0
        } else {
            self.friction_torque * self.output_velocity.signum()
        };
        amplified - friction
    }
    /// Advance dynamics by `dt`.
    ///
    /// * `input_torque` — motor input torque (N·m).
    pub fn step(&mut self, input_torque: f64, dt: f64) {
        let torque = self.output_torque(input_torque);
        let spring = -self.flexspline_stiffness * self.output_angle;
        let damp = -self.flexspline_damping * self.output_velocity;
        let net = torque + spring + damp;
        let alpha = net / self.output_inertia.max(1e-20);
        self.output_velocity += alpha * dt;
        self.output_angle += self.output_velocity * dt;
    }
    /// Input angle for a given output angle.
    pub fn input_angle(&self) -> f64 {
        self.output_angle * self.gear_ratio
    }
    /// Reflected inertia: output_inertia / ratio².
    pub fn reflected_inertia(&self) -> f64 {
        self.output_inertia / (self.gear_ratio * self.gear_ratio).max(1e-20)
    }
}
/// Rack-and-pinion linear actuator driven by a rotary motor.
///
/// The pinion of radius `r` converts motor torque to linear rack force:
/// F = T / r * efficiency.
#[derive(Debug, Clone)]
pub struct RackAndPinionMotor {
    /// Pinion radius (m).
    pub pinion_radius: f64,
    /// Mechanical efficiency (0..1).
    pub efficiency: f64,
    /// Maximum rack force (N).
    pub max_force: f64,
    /// Current rack position (m).
    pub position: f64,
    /// Current rack velocity (m/s).
    pub velocity: f64,
}
impl RackAndPinionMotor {
    /// Create a new rack-and-pinion motor.
    pub fn new(pinion_radius: f64, efficiency: f64, max_force: f64) -> Self {
        Self {
            pinion_radius,
            efficiency,
            max_force,
            position: 0.0,
            velocity: 0.0,
        }
    }
    /// Rack force from motor torque.
    pub fn rack_force(&self, motor_torque: f64) -> f64 {
        let raw = motor_torque * self.efficiency / self.pinion_radius.max(1e-12);
        raw.clamp(-self.max_force, self.max_force)
    }
    /// Motor torque needed for a given rack force.
    pub fn required_torque(&self, rack_force: f64) -> f64 {
        rack_force * self.pinion_radius / self.efficiency.max(1e-12)
    }
    /// Step the rack position given motor angular velocity.
    pub fn step(&mut self, angular_velocity: f64, dt: f64) {
        self.velocity = angular_velocity * self.pinion_radius;
        self.position += self.velocity * dt;
    }
}
/// Motor friction model combining Coulomb, viscous, and Stribeck effects.
#[derive(Debug, Clone)]
pub struct MotorFriction {
    /// Coulomb (dry) friction level (N or N·m).
    pub coulomb_friction: f64,
    /// Viscous friction coefficient (N·s/m or N·m·s/rad).
    pub viscous_coefficient: f64,
    /// Static friction level (N or N·m). Should be >= coulomb_friction.
    pub static_friction: f64,
    /// Stribeck velocity parameter (m/s or rad/s).
    pub stribeck_velocity: f64,
}
impl MotorFriction {
    /// Create a simple Coulomb + viscous friction model.
    pub fn coulomb_viscous(coulomb: f64, viscous: f64) -> Self {
        Self {
            coulomb_friction: coulomb,
            viscous_coefficient: viscous,
            static_friction: coulomb,
            stribeck_velocity: 0.01,
        }
    }
    /// Create a full Stribeck friction model.
    pub fn stribeck(coulomb: f64, viscous: f64, static_f: f64, stribeck_vel: f64) -> Self {
        Self {
            coulomb_friction: coulomb,
            viscous_coefficient: viscous,
            static_friction: static_f,
            stribeck_velocity: stribeck_vel,
        }
    }
    /// Compute friction force/torque for a given velocity.
    ///
    /// Uses the Stribeck friction curve:
    /// F_friction = (F_c + (F_s - F_c) * exp(-(v/v_s)^2)) * sign(v) + b * v
    pub fn friction_force(&self, velocity: f64) -> f64 {
        if velocity.abs() < 1e-15 {
            return 0.0;
        }
        let stribeck = self.coulomb_friction
            + (self.static_friction - self.coulomb_friction)
                * (-(velocity / self.stribeck_velocity).powi(2)).exp();
        stribeck * velocity.signum() + self.viscous_coefficient * velocity
    }
    /// Compute friction force with a velocity threshold for static regime.
    pub fn friction_force_with_threshold(&self, velocity: f64, threshold: f64) -> f64 {
        if velocity.abs() < threshold {
            return 0.0;
        }
        self.friction_force(velocity)
    }
}
/// Type of motor fault detected.
#[derive(Debug, Clone, PartialEq)]
pub enum MotorFault {
    /// No fault detected.
    None,
    /// Motor is overheated.
    OverTemperature,
    /// Motor current exceeds limit.
    OverCurrent,
    /// Position tracking error exceeds limit.
    PositionError,
    /// Velocity exceeds limit.
    OverSpeed,
    /// Stall detected (high current, near-zero speed).
    Stall,
}
/// Full impedance controller with virtual mass, stiffness, and damping.
///
/// The impedance relation is:
///   M_v * ẍ + D_v * ẋ + K_v * x = F_ext
///
/// The controller computes the motor force to realize this virtual dynamics.
#[derive(Debug, Clone)]
pub struct ImpedanceController {
    /// Virtual mass (kg or kg·m²).
    pub virtual_mass: f64,
    /// Virtual stiffness (N/m or N·m/rad).
    pub virtual_stiffness: f64,
    /// Virtual damping (N·s/m or N·m·s/rad).
    pub virtual_damping: f64,
    /// Maximum motor force/torque.
    pub max_force: f64,
    /// Desired end-effector position.
    pub desired_position: f64,
    /// Desired end-effector velocity.
    pub desired_velocity: f64,
    /// Desired end-effector acceleration.
    pub desired_accel: f64,
}
impl ImpedanceController {
    /// Create a new impedance controller.
    pub fn new(
        virtual_mass: f64,
        virtual_stiffness: f64,
        virtual_damping: f64,
        max_force: f64,
    ) -> Self {
        Self {
            virtual_mass,
            virtual_stiffness,
            virtual_damping,
            max_force,
            desired_position: 0.0,
            desired_velocity: 0.0,
            desired_accel: 0.0,
        }
    }
    /// Set the desired trajectory.
    pub fn set_desired(&mut self, pos: f64, vel: f64, accel: f64) {
        self.desired_position = pos;
        self.desired_velocity = vel;
        self.desired_accel = accel;
    }
    /// Compute the motor force given actual state and external force.
    ///
    /// F_cmd = M_v * ẍ_d + D_v * (ẋ_d - ẋ) + K_v * (x_d - x) - F_ext
    pub fn compute_force(&self, actual_pos: f64, actual_vel: f64, external_force: f64) -> f64 {
        let inertia_term = self.virtual_mass * self.desired_accel;
        let damping_term = self.virtual_damping * (self.desired_velocity - actual_vel);
        let spring_term = self.virtual_stiffness * (self.desired_position - actual_pos);
        let raw = inertia_term + damping_term + spring_term - external_force;
        raw.clamp(-self.max_force, self.max_force)
    }
    /// Natural frequency ω_n = sqrt(K_v / M_v).
    pub fn natural_frequency(&self) -> f64 {
        if self.virtual_mass < 1e-20 {
            return 0.0;
        }
        (self.virtual_stiffness / self.virtual_mass).sqrt()
    }
    /// Critical damping: D_c = 2 * sqrt(K_v * M_v).
    pub fn critical_damping(&self) -> f64 {
        2.0 * (self.virtual_stiffness * self.virtual_mass).sqrt()
    }
    /// Damping ratio.
    pub fn damping_ratio(&self) -> f64 {
        let dc = self.critical_damping();
        if dc < 1e-20 {
            return 0.0;
        }
        self.virtual_damping / dc
    }
}
/// Motor fault detection system.
#[derive(Debug, Clone)]
pub struct MotorFaultDetector {
    /// Maximum permissible winding temperature.
    pub max_temperature: f64,
    /// Maximum permissible current (A).
    pub max_current: f64,
    /// Maximum permissible position error (m or rad).
    pub max_position_error: f64,
    /// Maximum permissible speed (m/s or rad/s).
    pub max_speed: f64,
    /// Speed threshold below which stall check is triggered (m/s or rad/s).
    pub stall_speed_threshold: f64,
    /// Current threshold above which stall is triggered.
    pub stall_current_threshold: f64,
    /// Last detected fault.
    pub last_fault: MotorFault,
}
impl MotorFaultDetector {
    /// Create a new fault detector with given limits.
    pub fn new(
        max_temperature: f64,
        max_current: f64,
        max_position_error: f64,
        max_speed: f64,
    ) -> Self {
        Self {
            max_temperature,
            max_current,
            max_position_error,
            max_speed,
            stall_speed_threshold: 0.01,
            stall_current_threshold: 0.9,
            last_fault: MotorFault::None,
        }
    }
    /// Check all conditions and return the detected fault (or None).
    pub fn check(
        &mut self,
        temperature: f64,
        current: f64,
        position_error: f64,
        speed: f64,
    ) -> &MotorFault {
        self.last_fault = if temperature > self.max_temperature {
            MotorFault::OverTemperature
        } else if current.abs() > self.max_current {
            MotorFault::OverCurrent
        } else if speed.abs() > self.max_speed {
            MotorFault::OverSpeed
        } else if position_error.abs() > self.max_position_error {
            MotorFault::PositionError
        } else if speed.abs() < self.stall_speed_threshold
            && current.abs() > self.stall_current_threshold * self.max_current
        {
            MotorFault::Stall
        } else {
            MotorFault::None
        };
        &self.last_fault
    }
    /// Return true if a non-None fault is active.
    pub fn has_fault(&self) -> bool {
        self.last_fault != MotorFault::None
    }
    /// Reset the last detected fault to None.
    pub fn clear(&mut self) {
        self.last_fault = MotorFault::None;
    }
}
