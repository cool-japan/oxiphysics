//! Extended motor types: ElectricMotor, MotorBank, MotorWindingThermal,
//! CascadedPid, PidContributions, TrapezoidalProfile, PwmController,
//! GearTrain, TorqueSpring, MotorRundownIdentifier, BackEmfObserver,
//! RegenerativeBrakeController, HarmonicDrive, BldcMotor, SpringDamper,
//! MotorThermalModel, ZnTuningMethod.
use super::functions::{regenerated_power, regenerative_braking_torque, trapezoidal_wave};
use super::types::{PidController, TrapPhase};

/// Ziegler-Nichols tuning rule variant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ZnTuningMethod {
    /// Classic Z-N: highest loop bandwidth, quarter-decay ratio.
    Classic,
    /// Overshoot-reduced Z-N (less aggressive).
    OvershootReduced,
    /// No-overshoot Z-N (conservative).
    NoOvershoot,
    /// PI controller (Kd = 0).
    PiOnly,
}
/// Electric motor model with RPM-dependent efficiency and thermal dynamics.
pub struct ElectricMotor {
    /// Peak torque at stall (N*m).
    pub max_torque: f64,
    /// No-load speed (RPM).
    pub max_rpm: f64,
    /// Efficiency map as `(rpm, efficiency)` pairs sorted by rpm.
    pub efficiency_map: Vec<(f64, f64)>,
    /// Thermal state of the motor.
    pub thermal: MotorThermalModel,
    /// Nominal voltage (V).
    pub voltage: f64,
    /// Motor resistance (Ohms).
    pub resistance: f64,
}
impl ElectricMotor {
    /// Create a motor with uniform 90% efficiency if no map is provided.
    pub fn new(max_torque: f64, max_rpm: f64) -> Self {
        Self {
            max_torque,
            max_rpm,
            efficiency_map: vec![(0.0, 0.9), (max_rpm, 0.9)],
            thermal: MotorThermalModel::new(25.0, 150.0, 50.0, 0.5),
            voltage: 48.0,
            resistance: 0.1,
        }
    }
    /// Create a motor with a custom efficiency map.
    pub fn with_efficiency_map(mut self, map: Vec<(f64, f64)>) -> Self {
        self.efficiency_map = map;
        self
    }
    /// Set the thermal model parameters.
    pub fn with_thermal(mut self, thermal: MotorThermalModel) -> Self {
        self.thermal = thermal;
        self
    }
    /// Set voltage and resistance for electrical calculations.
    pub fn with_electrical(mut self, voltage: f64, resistance: f64) -> Self {
        self.voltage = voltage;
        self.resistance = resistance;
        self
    }
    /// Interpolate efficiency at the given `rpm`.
    fn efficiency_at(&self, rpm: f64) -> f64 {
        let map = &self.efficiency_map;
        if map.is_empty() {
            return 1.0;
        }
        if rpm <= map[0].0 {
            return map[0].1;
        }
        for i in 1..map.len() {
            if rpm <= map[i].0 {
                let t = (rpm - map[i - 1].0) / (map[i].0 - map[i - 1].0);
                return map[i - 1].1 + t * (map[i].1 - map[i - 1].1);
            }
        }
        map[map.len() - 1].1
    }
    /// Torque output at the given `rpm` and `throttle` (0..1).
    ///
    /// Models a linear speed-torque curve: torque decreases from `max_torque`
    /// at stall to 0 at `max_rpm`, scaled by `throttle` and efficiency.
    pub fn torque_at(&self, rpm: f64, throttle: f64) -> f64 {
        let throttle = throttle.clamp(0.0, 1.0);
        let rpm = rpm.max(0.0);
        let base = self.max_torque * (1.0 - (rpm / self.max_rpm).min(1.0));
        base * throttle * self.efficiency_at(rpm)
    }
    /// Mechanical power output (W) at the given `rpm` and `throttle`.
    pub fn power(&self, rpm: f64, throttle: f64) -> f64 {
        let torque = self.torque_at(rpm, throttle);
        let omega = rpm * (2.0 * std::f64::consts::PI / 60.0);
        torque * omega
    }
    /// Electrical power input (W) for the given `rpm` and `throttle`.
    ///
    /// P_elec = P_mech / efficiency.
    pub fn electrical_power(&self, rpm: f64, throttle: f64) -> f64 {
        let eff = self.efficiency_at(rpm);
        if eff > 1e-10 {
            self.power(rpm, throttle) / eff
        } else {
            0.0
        }
    }
    /// Current draw (A) at the given operating point.
    ///
    /// I = P_elec / V.
    pub fn current(&self, rpm: f64, throttle: f64) -> f64 {
        if self.voltage > 1e-10 {
            self.electrical_power(rpm, throttle) / self.voltage
        } else {
            0.0
        }
    }
    /// Resistive heat loss (W) = I^2 * R.
    pub fn heat_loss(&self, rpm: f64, throttle: f64) -> f64 {
        let i = self.current(rpm, throttle);
        i * i * self.resistance
    }
    /// Update thermal model with one time step.
    pub fn step_thermal(&mut self, rpm: f64, throttle: f64, dt: f64) {
        let heat = self.heat_loss(rpm, throttle);
        self.thermal.step(heat, dt);
    }
    /// Whether the motor is overheating (temperature above max).
    pub fn is_overheating(&self) -> bool {
        self.thermal.is_overheating()
    }
    /// Peak RPM at which the motor can still produce the given torque.
    pub fn rpm_at_torque(&self, target_torque: f64) -> f64 {
        let eff = self.efficiency_at(self.max_rpm * 0.5);
        let ratio = target_torque / (self.max_torque * eff);
        if ratio >= 1.0 {
            return 0.0;
        }
        self.max_rpm * (1.0 - ratio).max(0.0)
    }
    /// Compute regenerative braking torque (N·m) when the motor acts as a
    /// generator during deceleration.
    ///
    /// Models the generator braking torque as:
    ///
    /// ```text
    /// T_regen = throttle_braking * max_torque * min(rpm / rpm_peak_regen, 1)
    ///           * efficiency(rpm) * thermal_derating_factor
    /// ```
    ///
    /// * `rpm` — current shaft speed (RPM), must be positive.
    /// * `throttle_braking` — normalised braking demand in `[0, 1]`.
    /// * `rpm_peak_regen` — speed at which regen torque reaches its peak.
    ///
    /// Returns the magnitude of the braking torque (positive value that
    /// opposes rotation).
    pub fn compute_regenerative_braking(
        &self,
        rpm: f64,
        throttle_braking: f64,
        rpm_peak_regen: f64,
    ) -> f64 {
        let rpm = rpm.max(0.0);
        let throttle = throttle_braking.clamp(0.0, 1.0);
        if rpm < 1e-6 || rpm_peak_regen < 1e-6 {
            return 0.0;
        }
        let ramp = (rpm / rpm_peak_regen).min(1.0);
        let eff = self.efficiency_at(rpm);
        let thermal = self.thermal_derating_factor();
        self.max_torque * throttle * ramp * eff * thermal
    }
    /// Thermal derating factor in `[0, 1]`.
    ///
    /// Returns `1.0` when the motor is cool, linearly decreasing to `0.0`
    /// when the temperature reaches `max_temperature`.
    ///
    /// The derating curve activates only above a knee temperature defined as
    /// `ambient + 0.7 * (max - ambient)` (i.e. at 70 % of the thermal
    /// headroom).  Below the knee the factor is `1.0`.
    pub fn thermal_derating_factor(&self) -> f64 {
        let t = self.thermal.temperature;
        let t_ambient = self.thermal.ambient_temperature;
        let t_max = self.thermal.max_temperature;
        let headroom = t_max - t_ambient;
        if headroom < 1e-10 {
            return 0.0;
        }
        let knee = t_ambient + 0.7 * headroom;
        if t <= knee {
            return 1.0;
        }
        let derating_range = t_max - knee;
        if derating_range < 1e-10 {
            return 0.0;
        }
        (1.0 - (t - knee) / derating_range).max(0.0)
    }
    /// Compute the temperature derating curve at a given temperature.
    ///
    /// Stateless version of `compute_regenerative_braking` that returns the
    /// derating factor for an arbitrary `temperature` (C) given the motor's
    /// thermal parameters.
    ///
    /// * `temperature`       — temperature to evaluate (C).
    /// * `ambient`           — ambient temperature (C).
    /// * `max_temperature`   — maximum allowable temperature (C).
    ///
    /// Returns the derating factor in `[0.0, 1.0]`.
    pub fn compute_temperature_derating(
        temperature: f64,
        ambient: f64,
        max_temperature: f64,
    ) -> f64 {
        let headroom = max_temperature - ambient;
        if headroom < 1e-10 {
            return 0.0;
        }
        let knee = ambient + 0.7 * headroom;
        if temperature <= knee {
            return 1.0;
        }
        let derating_range = max_temperature - knee;
        if derating_range < 1e-10 {
            return 0.0;
        }
        (1.0 - (temperature - knee) / derating_range).max(0.0)
    }
}
/// A collection of PID controllers for multi-joint systems.
pub struct MotorBank {
    /// PID controllers, one per joint.
    pub controllers: Vec<PidController>,
}
impl MotorBank {
    /// Create a motor bank with `n` identical PID controllers.
    pub fn new(n: usize, kp: f64, ki: f64, kd: f64, max_output: f64) -> Self {
        let controllers = (0..n)
            .map(|_| PidController::new(kp, ki, kd, max_output))
            .collect();
        Self { controllers }
    }
    /// Update all controllers with corresponding errors.
    ///
    /// Returns the control outputs.
    pub fn update(&mut self, errors: &[f64], dt: f64) -> Vec<f64> {
        self.controllers
            .iter_mut()
            .zip(errors.iter())
            .map(|(ctrl, &err)| ctrl.update(err, dt))
            .collect()
    }
    /// Reset all controllers.
    pub fn reset_all(&mut self) {
        for ctrl in &mut self.controllers {
            ctrl.reset();
        }
    }
    /// Number of controllers in the bank.
    pub fn len(&self) -> usize {
        self.controllers.len()
    }
    /// Whether the bank is empty.
    pub fn is_empty(&self) -> bool {
        self.controllers.is_empty()
    }
}
/// Simple thermal model for a motor winding.
///
/// Thermal circuit: winding → housing → ambient.
/// dT_w/dt = P_loss / C_w - (T_w - T_h) / R_wh
/// dT_h/dt = (T_w - T_h) / (R_wh * C_h) - (T_h - T_amb) / (R_ha * C_h)
#[derive(Debug, Clone)]
pub struct MotorWindingThermal {
    /// Winding thermal resistance winding→housing (K/W).
    pub r_winding_housing: f64,
    /// Housing thermal resistance housing→ambient (K/W).
    pub r_housing_ambient: f64,
    /// Winding thermal capacitance (J/K).
    pub c_winding: f64,
    /// Housing thermal capacitance (J/K).
    pub c_housing: f64,
    /// Ambient temperature (°C).
    pub t_ambient: f64,
    /// Current winding temperature (°C).
    pub t_winding: f64,
    /// Current housing temperature (°C).
    pub t_housing: f64,
}
impl MotorWindingThermal {
    /// Create a new thermal model.
    pub fn new(
        r_winding_housing: f64,
        r_housing_ambient: f64,
        c_winding: f64,
        c_housing: f64,
        t_ambient: f64,
    ) -> Self {
        Self {
            r_winding_housing,
            r_housing_ambient,
            c_winding,
            c_housing,
            t_ambient,
            t_winding: t_ambient,
            t_housing: t_ambient,
        }
    }
    /// Integrate one time step with winding power loss `p_loss` (W).
    pub fn step(&mut self, p_loss: f64, dt: f64) {
        let q_wh = (self.t_winding - self.t_housing) / self.r_winding_housing.max(f64::EPSILON);
        let q_ha = (self.t_housing - self.t_ambient) / self.r_housing_ambient.max(f64::EPSILON);
        let dt_w = (p_loss - q_wh) / self.c_winding.max(f64::EPSILON);
        let dt_h = (q_wh - q_ha) / self.c_housing.max(f64::EPSILON);
        self.t_winding += dt_w * dt;
        self.t_housing += dt_h * dt;
    }
    /// Check whether the winding temperature is within limits.
    pub fn within_limit(&self, max_temp: f64) -> bool {
        self.t_winding <= max_temp
    }
    /// Steady-state winding temperature rise (K) for constant power `p`.
    pub fn steady_state_rise(&self, p: f64) -> f64 {
        p * (self.r_winding_housing + self.r_housing_ambient)
    }
}
/// Gain configuration for a [`CascadedPid`] controller.
///
/// Groups the 12 PID parameters (Kp, Ki, Kd, max for each of three loops) into
/// a single struct so that [`CascadedPid::new`] stays within the argument-count
/// limit.
#[derive(Debug, Clone, Copy)]
pub struct CascadedPidGains {
    /// Position loop proportional gain
    pub pos_kp: f64,
    /// Position loop integral gain
    pub pos_ki: f64,
    /// Position loop derivative gain
    pub pos_kd: f64,
    /// Position loop output saturation (rad/s)
    pub pos_max: f64,
    /// Velocity loop proportional gain
    pub vel_kp: f64,
    /// Velocity loop integral gain
    pub vel_ki: f64,
    /// Velocity loop derivative gain
    pub vel_kd: f64,
    /// Velocity loop output saturation (A)
    pub vel_max: f64,
    /// Current loop proportional gain
    pub cur_kp: f64,
    /// Current loop integral gain
    pub cur_ki: f64,
    /// Current loop derivative gain
    pub cur_kd: f64,
    /// Current loop output saturation (V)
    pub cur_max: f64,
}

/// Cascaded PID controller: outer position loop drives inner velocity loop,
/// which in turn drives a current (torque) command.
///
/// # Physical motivation
/// High-performance servo drives use a cascade architecture to separately tune
/// bandwidth for position tracking, velocity bandwidth, and current bandwidth.
pub struct CascadedPid {
    /// Outer position loop controller.
    pub position_pid: PidController,
    /// Middle velocity loop controller.
    pub velocity_pid: PidController,
    /// Inner current loop controller.
    pub current_pid: PidController,
    /// Saturation limit for velocity command (rad/s).
    pub vel_cmd_limit: f64,
    /// Saturation limit for current command (A).
    pub current_cmd_limit: f64,
}
impl CascadedPid {
    /// Create a cascaded PID from a [`CascadedPidGains`] configuration bundle.
    pub fn new(g: CascadedPidGains) -> Self {
        Self {
            position_pid: PidController::new(g.pos_kp, g.pos_ki, g.pos_kd, g.pos_max),
            velocity_pid: PidController::new(g.vel_kp, g.vel_ki, g.vel_kd, g.vel_max),
            current_pid: PidController::new(g.cur_kp, g.cur_ki, g.cur_kd, g.cur_max),
            vel_cmd_limit: g.vel_max,
            current_cmd_limit: g.cur_max,
        }
    }
    /// Run one update step.
    ///
    /// * `pos_target` – desired position (rad or m)
    /// * `pos_actual` – measured position
    /// * `vel_actual` – measured velocity (rad/s or m/s)
    /// * `cur_actual` – measured current (A)
    /// * `dt`         – timestep (s)
    ///
    /// Returns the voltage/torque command from the innermost loop.
    pub fn update(
        &mut self,
        pos_target: f64,
        pos_actual: f64,
        vel_actual: f64,
        cur_actual: f64,
        dt: f64,
    ) -> f64 {
        let pos_error = pos_target - pos_actual;
        let vel_cmd = self
            .position_pid
            .update(pos_error, dt)
            .clamp(-self.vel_cmd_limit, self.vel_cmd_limit);
        let vel_error = vel_cmd - vel_actual;
        let cur_cmd = self
            .velocity_pid
            .update(vel_error, dt)
            .clamp(-self.current_cmd_limit, self.current_cmd_limit);
        let cur_error = cur_cmd - cur_actual;
        self.current_pid.update(cur_error, dt)
    }
    /// Reset all three controllers.
    pub fn reset(&mut self) {
        self.position_pid.reset();
        self.velocity_pid.reset();
        self.current_pid.reset();
    }
}
/// Diagnostic breakdown of PID controller contributions.
#[derive(Debug, Clone)]
/// Trapezoidal velocity profile generator.
///
/// Generates position, velocity, and acceleration commands for smooth
/// point-to-point motion with bounded acceleration.
pub struct TrapezoidalProfile {
    /// Maximum cruise velocity (units/s).
    pub v_max: f64,
    /// Maximum acceleration / deceleration (units/s²).
    pub a_max: f64,
    /// Start position.
    pub x_start: f64,
    /// Target position.
    pub x_target: f64,
    /// Direction (+1 or -1).
    pub direction: f64,
    /// Time of acceleration phase end.
    pub t1: f64,
    /// Time of deceleration phase start.
    pub t2: f64,
    /// Total move duration.
    pub t_total: f64,
    /// Actual cruise velocity (may be less than v_max for short moves).
    pub v_cruise: f64,
}
impl TrapezoidalProfile {
    /// Plan a move from `x_start` to `x_target`.
    pub fn plan(x_start: f64, x_target: f64, v_max: f64, a_max: f64) -> Self {
        let dist = x_target - x_start;
        let direction = if dist >= 0.0 { 1.0 } else { -1.0 };
        let s = dist.abs();
        let s_triangle = v_max * v_max / a_max;
        let (v_cruise, t1, t2, t_total) = if s >= s_triangle {
            let t_acc = v_max / a_max;
            let s_cruise = s - s_triangle;
            let t_cru = s_cruise / v_max;
            (v_max, t_acc, t_acc + t_cru, 2.0 * t_acc + t_cru)
        } else {
            let v_peak = (a_max * s).sqrt();
            let t_acc = v_peak / a_max;
            (v_peak, t_acc, t_acc, 2.0 * t_acc)
        };
        Self {
            v_max,
            a_max,
            x_start,
            x_target,
            direction,
            t1,
            t2,
            t_total,
            v_cruise,
        }
    }
    /// Query position, velocity, and acceleration at time `t`.
    ///
    /// Returns `(position, velocity, acceleration)`.
    pub fn query(&self, t: f64) -> (f64, f64, f64) {
        if t <= 0.0 {
            return (self.x_start, 0.0, 0.0);
        }
        if t >= self.t_total {
            return (self.x_target, 0.0, 0.0);
        }
        let (pos_rel, vel, acc) = if t <= self.t1 {
            let vel = self.a_max * t;
            let pos = 0.5 * self.a_max * t * t;
            (pos, vel, self.a_max)
        } else if t <= self.t2 {
            let pos = 0.5 * self.a_max * self.t1 * self.t1 + self.v_cruise * (t - self.t1);
            (pos, self.v_cruise, 0.0)
        } else {
            let dt = t - self.t2;
            let vel = self.v_cruise - self.a_max * dt;
            let pos = 0.5 * self.a_max * self.t1 * self.t1
                + self.v_cruise * (self.t2 - self.t1)
                + self.v_cruise * dt
                - 0.5 * self.a_max * dt * dt;
            (pos, vel.max(0.0), -self.a_max)
        };
        (
            self.x_start + self.direction * pos_rel,
            self.direction * vel,
            self.direction * acc,
        )
    }
    /// Return the current phase for time `t`.
    pub fn phase(&self, t: f64) -> TrapPhase {
        if t >= self.t_total {
            TrapPhase::Done
        } else if t < self.t1 {
            TrapPhase::Accelerating
        } else if t < self.t2 {
            TrapPhase::Cruising
        } else {
            TrapPhase::Decelerating
        }
    }
}
/// PWM (Pulse-Width Modulation) controller that converts a duty cycle to an
/// average bus voltage, with optional dead-time compensation.
///
/// The average output voltage is:
/// ```text
/// V_avg = duty * V_bus  (0 ≤ duty ≤ 1)
/// ```
///
/// Dead-time compensation reduces the effective duty cycle to account for the
/// blanking time in the gate drive circuitry.
#[derive(Debug, Clone)]
pub struct PwmController {
    /// DC bus voltage (V).
    pub v_bus: f64,
    /// PWM switching frequency (Hz).
    pub switching_freq: f64,
    /// Dead-time (s) per half-period.
    pub dead_time: f64,
    /// Current duty cycle ∈ \[0, 1\].
    pub duty: f64,
}
impl PwmController {
    /// Create a new PWM controller.
    pub fn new(v_bus: f64, switching_freq: f64, dead_time: f64) -> Self {
        Self {
            v_bus,
            switching_freq,
            dead_time,
            duty: 0.0,
        }
    }
    /// Set the duty cycle.  Clamped to \[0, 1\].
    pub fn set_duty(&mut self, duty: f64) {
        self.duty = duty.clamp(0.0, 1.0);
    }
    /// Average output voltage without dead-time compensation.
    pub fn average_voltage(&self) -> f64 {
        self.duty * self.v_bus
    }
    /// Effective duty cycle after dead-time compensation.
    ///
    /// For a full-bridge inverter, the dead-time reduces the effective duty
    /// by `dead_time * switching_freq` per half-cycle:
    ///
    /// ```text
    /// duty_eff = duty - dead_time * switching_freq
    /// ```
    ///
    /// Clamped to \[0, 1\].
    pub fn effective_duty(&self) -> f64 {
        (self.duty - self.dead_time * self.switching_freq).clamp(0.0, 1.0)
    }
    /// Average output voltage with dead-time compensation.
    pub fn compensated_voltage(&self) -> f64 {
        self.effective_duty() * self.v_bus
    }
    /// Peak-to-peak current ripple for a given inductance `l` (H) and
    /// current `i_mean` (A) at the current operating point.
    ///
    /// ```text
    /// ΔI = V_bus * duty * (1 - duty) / (l * f_sw)
    /// ```
    pub fn current_ripple(&self, l: f64) -> f64 {
        if l < f64::EPSILON || self.switching_freq < f64::EPSILON {
            return 0.0;
        }
        self.v_bus * self.duty * (1.0 - self.duty) / (l * self.switching_freq)
    }
    /// Switching loss estimate (W) given conduction/switching energy per cycle
    /// `e_sw` (J).
    pub fn switching_loss(&self, e_sw: f64) -> f64 {
        e_sw * self.switching_freq
    }
    /// Conduction loss (W) for resistive load `r` at average duty cycle.
    pub fn conduction_loss(&self, r: f64) -> f64 {
        if r < f64::EPSILON {
            return 0.0;
        }
        let v = self.average_voltage();
        v * v / r
    }
    /// Modulate: given a target average voltage, compute and set the required
    /// duty cycle.  Returns the actual set duty.
    pub fn modulate(&mut self, v_target: f64) -> f64 {
        if self.v_bus.abs() < f64::EPSILON {
            self.duty = 0.0;
        } else {
            self.duty = (v_target / self.v_bus).clamp(0.0, 1.0);
        }
        self.duty
    }
}
/// Simple spur-gear (or belt-pulley) transmission model.
///
/// Relates input shaft to output shaft through a gear ratio and efficiency.
/// gear_ratio = ω_in / ω_out  (> 1 means speed reduction, torque amplification).
pub struct GearTrain {
    /// Gear ratio N = ω_in / ω_out.
    pub gear_ratio: f64,
    /// Transmission efficiency η ∈ (0, 1].
    pub efficiency: f64,
    /// Maximum continuous output torque (N·m).
    pub max_torque: f64,
    /// Reflected inertia factor: J_reflected = J_load / N².
    pub load_inertia: f64,
}
impl GearTrain {
    /// Create a new gear train.
    pub fn new(gear_ratio: f64, efficiency: f64, max_torque: f64, load_inertia: f64) -> Self {
        Self {
            gear_ratio,
            efficiency,
            max_torque,
            load_inertia,
        }
    }
    /// Output speed for a given input speed.
    pub fn output_speed(&self, omega_in: f64) -> f64 {
        omega_in / self.gear_ratio
    }
    /// Input speed for a given desired output speed.
    pub fn input_speed(&self, omega_out: f64) -> f64 {
        omega_out * self.gear_ratio
    }
    /// Output torque amplified and derated by efficiency.
    ///
    /// τ_out = τ_in × N × η
    pub fn output_torque(&self, tau_in: f64) -> f64 {
        tau_in * self.gear_ratio * self.efficiency
    }
    /// Input torque required to produce a given output torque.
    ///
    /// τ_in = τ_out / (N × η)
    pub fn input_torque_required(&self, tau_out: f64) -> f64 {
        tau_out / (self.gear_ratio * self.efficiency).max(f64::EPSILON)
    }
    /// Reflected load inertia seen at the motor shaft: J_refl = J_load / N².
    pub fn reflected_inertia(&self) -> f64 {
        self.load_inertia / (self.gear_ratio * self.gear_ratio).max(f64::EPSILON)
    }
    /// True when output torque would exceed the maximum rating.
    pub fn is_overloaded(&self, tau_in: f64) -> bool {
        self.output_torque(tau_in) > self.max_torque
    }
    /// Power loss in the gear train (W).
    pub fn power_loss(&self, tau_in: f64, omega_in: f64) -> f64 {
        let p_in = tau_in * omega_in;
        let p_out = self.output_torque(tau_in) * self.output_speed(omega_in);
        (p_in - p_out).max(0.0)
    }
    /// Combined motor+gear inertia at motor shaft.
    pub fn total_inertia(&self, motor_inertia: f64) -> f64 {
        motor_inertia + self.reflected_inertia()
    }
}
/// Rotational spring-damper element.
pub struct TorqueSpring {
    /// Torsional stiffness (N*m/rad).
    pub stiffness: f64,
    /// Torsional damping coefficient (N*m*s/rad).
    pub damping: f64,
    /// Rest angle (rad).
    pub rest_angle: f64,
}
impl TorqueSpring {
    /// Create a new torque spring.
    pub fn new(stiffness: f64, damping: f64, rest_angle: f64) -> Self {
        Self {
            stiffness,
            damping,
            rest_angle,
        }
    }
    /// Compute the restoring torque for the given angle and angular velocity.
    ///
    /// `torque = -stiffness * (angle - rest_angle) - damping * angular_vel`
    pub fn torque(&self, angle: f64, angular_vel: f64) -> f64 {
        -self.stiffness * (angle - self.rest_angle) - self.damping * angular_vel
    }
    /// Elastic potential energy stored in the spring.
    pub fn energy(&self, angle: f64) -> f64 {
        let da = angle - self.rest_angle;
        0.5 * self.stiffness * da * da
    }
    /// Natural frequency (rad/s) for a given inertia.
    pub fn natural_frequency(&self, inertia: f64) -> f64 {
        if inertia > 1e-15 {
            (self.stiffness / inertia).sqrt()
        } else {
            0.0
        }
    }
    /// Damping ratio for a given inertia.
    pub fn damping_ratio(&self, inertia: f64) -> f64 {
        let omega_n = self.natural_frequency(inertia);
        if omega_n > 1e-15 {
            self.damping / (2.0 * inertia * omega_n)
        } else {
            0.0
        }
    }
    /// Whether the system is overdamped for a given inertia.
    pub fn is_overdamped(&self, inertia: f64) -> bool {
        self.damping_ratio(inertia) > 1.0
    }
    /// Whether the system is critically damped for a given inertia.
    pub fn is_critically_damped(&self, inertia: f64) -> bool {
        (self.damping_ratio(inertia) - 1.0).abs() < 1e-6
    }
}
/// Identifies motor parameters (viscous friction `B`, inertia `J`, and
/// back-EMF constant `K_e`) from a free-coast-down speed profile.
///
/// During a coast-down test the motor shaft decelerates under its own friction
/// and back-EMF braking with open-circuit terminals:
///
/// ```text
/// J · dω/dt = −B · ω
/// ω(t) = ω₀ · exp(−(B/J) · t)
/// ```
///
/// A least-squares exponential fit over the sampled `(t, ω)` pairs estimates
/// the time constant `τ = J/B` and initial speed `ω₀`.
#[derive(Debug, Clone)]
pub struct MotorRundownIdentifier {
    /// Sampled (time, speed) pairs collected during coast-down.
    pub samples: Vec<(f64, f64)>,
}
impl MotorRundownIdentifier {
    /// Create a new identifier with no samples.
    pub fn new() -> Self {
        Self {
            samples: Vec::new(),
        }
    }
    /// Record a speed sample at time `t` (s), speed `omega` (rad/s).
    pub fn push(&mut self, t: f64, omega: f64) {
        self.samples.push((t, omega));
    }
    /// Estimate the exponential decay time constant `τ = J/B` (s) from the
    /// collected samples using a linear regression on `ln(ω)` vs `t`.
    ///
    /// Returns `None` if fewer than 2 samples are available or if all speeds
    /// are non-positive.
    pub fn time_constant(&self) -> Option<f64> {
        let n = self.samples.len();
        if n < 2 {
            return None;
        }
        let mut sum_t = 0.0f64;
        let mut sum_y = 0.0f64;
        let mut sum_tt = 0.0f64;
        let mut sum_ty = 0.0f64;
        let mut count = 0usize;
        for &(t, omega) in &self.samples {
            if omega > 1e-12 {
                let y = omega.ln();
                sum_t += t;
                sum_y += y;
                sum_tt += t * t;
                sum_ty += t * y;
                count += 1;
            }
        }
        if count < 2 {
            return None;
        }
        let n_f = count as f64;
        let denom = n_f * sum_tt - sum_t * sum_t;
        if denom.abs() < 1e-30 {
            return None;
        }
        let slope = (n_f * sum_ty - sum_t * sum_y) / denom;
        if slope >= 0.0 {
            return None;
        }
        Some(-1.0 / slope)
    }
    /// Estimate the initial speed `ω₀` (rad/s) from the regression intercept.
    pub fn initial_speed(&self) -> Option<f64> {
        let n = self.samples.len();
        if n < 2 {
            return None;
        }
        let mut sum_t = 0.0f64;
        let mut sum_y = 0.0f64;
        let mut sum_tt = 0.0f64;
        let mut sum_ty = 0.0f64;
        let mut count = 0usize;
        for &(t, omega) in &self.samples {
            if omega > 1e-12 {
                let y = omega.ln();
                sum_t += t;
                sum_y += y;
                sum_tt += t * t;
                sum_ty += t * y;
                count += 1;
            }
        }
        if count < 2 {
            return None;
        }
        let n_f = count as f64;
        let denom = n_f * sum_tt - sum_t * sum_t;
        if denom.abs() < 1e-30 {
            return None;
        }
        let intercept = (sum_y * sum_tt - sum_t * sum_ty) / denom;
        Some(intercept.exp())
    }
    /// Identify J and B given measured stall torque `tau_stall` (N·m) and
    /// voltage `v` (V) and the coast-down time constant.
    ///
    /// From `τ = J/B` and the stall torque relationship `B = Kt*I_stall / ω_ss`,
    /// this provides a rough `(J, B)` pair.
    ///
    /// Returns `(J, B)` or `None` if the time constant cannot be computed.
    pub fn identify_jb(&self, tau_stall: f64, no_load_speed: f64) -> Option<(f64, f64)> {
        let tau = self.time_constant()?;
        if no_load_speed < 1e-12 {
            return None;
        }
        let b = tau_stall / no_load_speed;
        let j = tau * b;
        Some((j, b))
    }
    /// Predicted speed at time `t` given the fitted parameters.
    pub fn predict_speed(&self, t: f64) -> Option<f64> {
        let tau = self.time_constant()?;
        let omega0 = self.initial_speed()?;
        Some(omega0 * (-t / tau).exp())
    }
    /// Clear all samples.
    pub fn clear(&mut self) {
        self.samples.clear();
    }
}
/// A single gear stage with a gear ratio and efficiency.
#[derive(Debug, Clone)]
/// Simple back-EMF observer for sensorless DC/BLDC motor speed estimation.
///
/// Estimates rotor speed `ω` from terminal voltage `V`, armature current `I`,
/// and motor constants using:
///
/// ```text
/// E_emf = V - R_a·I - L_a·(dI/dt)
/// ω_est = E_emf / K_e
/// ```
pub struct BackEmfObserver {
    /// Armature resistance (Ω).
    pub r_a: f64,
    /// Armature inductance (H).
    pub l_a: f64,
    /// Back-EMF constant K_e (V·s/rad).
    pub k_e: f64,
    /// Low-pass filter coefficient α ∈ (0, 1].
    pub alpha: f64,
    /// Previous current sample (A) for derivative estimation.
    prev_current: f64,
    /// Estimated speed (rad/s).
    pub omega_est: f64,
}
impl BackEmfObserver {
    /// Create a new back-EMF observer.
    ///
    /// `alpha` is the IIR filter coefficient.  Set to 1.0 for no filtering.
    pub fn new(r_a: f64, l_a: f64, k_e: f64, alpha: f64) -> Self {
        Self {
            r_a,
            l_a,
            k_e,
            alpha,
            prev_current: 0.0,
            omega_est: 0.0,
        }
    }
    /// Update the speed estimate given terminal voltage `V`, current `I`, and
    /// timestep `dt`.
    ///
    /// Returns the new estimated speed (rad/s).
    pub fn update(&mut self, v: f64, i: f64, dt: f64) -> f64 {
        let d_i_dt = if dt > 1e-15 {
            (i - self.prev_current) / dt
        } else {
            0.0
        };
        let e_emf = v - self.r_a * i - self.l_a * d_i_dt;
        let omega_raw = if self.k_e.abs() > 1e-15 {
            e_emf / self.k_e
        } else {
            0.0
        };
        self.omega_est = self.alpha * omega_raw + (1.0 - self.alpha) * self.omega_est;
        self.prev_current = i;
        self.omega_est
    }
    /// Reset observer state.
    pub fn reset(&mut self) {
        self.prev_current = 0.0;
        self.omega_est = 0.0;
    }
}
/// State machine for regenerative braking control.
///
/// Decides whether to apply regenerative braking or friction braking based
/// on speed and battery state of charge.
pub struct RegenerativeBrakeController {
    /// Minimum speed below which regeneration is disabled (rad/s).
    pub min_regen_speed: f64,
    /// Back-EMF constant K_e (V·s/rad).
    pub ke: f64,
    /// Torque constant K_t (N·m/A).
    pub kt: f64,
    /// Armature resistance (Ω).
    pub r_armature: f64,
    /// External braking resistance (Ω).
    pub r_brake: f64,
    /// Maximum regenerated current (A).
    pub max_regen_current: f64,
    /// Battery state of charge (0..1). If > `soc_cutoff`, regen is limited.
    pub soc: f64,
    /// SoC above which regen is cut off (e.g. 0.95).
    pub soc_cutoff: f64,
    /// Total energy recovered (J).
    pub energy_recovered: f64,
}
impl RegenerativeBrakeController {
    /// Create a new regenerative brake controller.
    pub fn new(
        min_regen_speed: f64,
        ke: f64,
        kt: f64,
        r_armature: f64,
        r_brake: f64,
        max_regen_current: f64,
        soc_cutoff: f64,
    ) -> Self {
        Self {
            min_regen_speed,
            ke,
            kt,
            r_armature,
            r_brake,
            max_regen_current,
            soc: 0.5,
            soc_cutoff,
            energy_recovered: 0.0,
        }
    }
    /// Compute regenerative braking torque and update energy counter.
    ///
    /// Returns `(braking_torque N·m, regen_power W)`.
    pub fn update(&mut self, omega: f64, dt: f64) -> (f64, f64) {
        if omega.abs() < self.min_regen_speed || self.soc >= self.soc_cutoff {
            return (0.0, 0.0);
        }
        let torque =
            regenerative_braking_torque(omega, self.ke, self.kt, self.r_armature, self.r_brake);
        let power = regenerated_power(omega, self.ke, self.r_armature, self.r_brake).max(0.0);
        self.energy_recovered += power * dt;
        (torque, power)
    }
    /// Reset energy counter.
    pub fn reset_energy(&mut self) {
        self.energy_recovered = 0.0;
    }
}
/// Harmonic drive (strain wave gear) compliance model.
///
/// Models the torsional spring behaviour of the flex spline:
/// τ = K_hd * (θ_in / ratio - θ_out)
#[derive(Debug, Clone)]
pub struct HarmonicDrive {
    /// Gear ratio (e.g. 100 for a 100:1 harmonic drive).
    pub ratio: f64,
    /// Torsional stiffness of the flex spline (N·m/rad).
    pub stiffness: f64,
    /// Torsional damping (N·m·s/rad).
    pub damping: f64,
    /// Input angle (rad).
    pub theta_in: f64,
    /// Output angle (rad).
    pub theta_out: f64,
    /// Output angular velocity (rad/s).
    pub omega_out: f64,
}
impl HarmonicDrive {
    /// Create a harmonic drive model.
    pub fn new(ratio: f64, stiffness: f64, damping: f64) -> Self {
        Self {
            ratio,
            stiffness,
            damping,
            theta_in: 0.0,
            theta_out: 0.0,
            omega_out: 0.0,
        }
    }
    /// Compute the output torque given current state.
    pub fn output_torque(&self) -> f64 {
        let deflection = self.theta_in / self.ratio - self.theta_out;
        self.stiffness * deflection - self.damping * self.omega_out
    }
    /// Integrate one step: update input angle, then output.
    pub fn step(&mut self, d_theta_in: f64, load_inertia: f64, dt: f64) {
        self.theta_in += d_theta_in;
        let tau = self.output_torque();
        let domega = tau / load_inertia.max(f64::EPSILON);
        self.omega_out += domega * dt;
        self.theta_out += self.omega_out * dt;
    }
}
/// A simplified 3-phase BLDC motor model.
///
/// Uses trapezoidal back-EMF and a six-step commutation table.
#[derive(Debug, Clone)]
pub struct BldcMotor {
    /// Phase resistance (Ω) per phase.
    pub phase_resistance: f64,
    /// Phase inductance (H) per phase.
    pub phase_inductance: f64,
    /// Back-EMF constant Ke (V·s/rad).
    pub ke: f64,
    /// Torque constant Kt (N·m/A peak).
    pub kt: f64,
    /// Number of pole pairs.
    pub pole_pairs: u32,
    /// Current electrical angle (rad).
    pub theta_e: f64,
    /// Current angular velocity (rad/s, mechanical).
    pub omega: f64,
    /// Phase currents \[I_a, I_b, I_c\].
    pub currents: [f64; 3],
}
impl BldcMotor {
    /// Create a BLDC motor.
    pub fn new(
        phase_resistance: f64,
        phase_inductance: f64,
        ke: f64,
        kt: f64,
        pole_pairs: u32,
    ) -> Self {
        Self {
            phase_resistance,
            phase_inductance,
            ke,
            kt,
            pole_pairs,
            theta_e: 0.0,
            omega: 0.0,
            currents: [0.0; 3],
        }
    }
    /// Advance electrical angle by mechanical step.
    pub fn update_angle(&mut self, dt: f64) {
        self.theta_e += self.omega * self.pole_pairs as f64 * dt;
        let two_pi = 2.0 * std::f64::consts::PI;
        self.theta_e = ((self.theta_e % two_pi) + two_pi) % two_pi;
    }
    /// Get trapezoidal back-EMF waveform values for each phase at `theta_e`.
    pub fn back_emf_phases(&self) -> [f64; 3] {
        let two_pi = 2.0 * std::f64::consts::PI;
        [
            self.ke * trapezoidal_wave(self.theta_e),
            self.ke * trapezoidal_wave(self.theta_e - two_pi / 3.0),
            self.ke * trapezoidal_wave(self.theta_e + two_pi / 3.0),
        ]
    }
    /// Compute total electromagnetic torque.
    pub fn torque(&self) -> f64 {
        let emf = self.back_emf_phases();
        let i_mag = (self.currents[0] * self.currents[0]
            + self.currents[1] * self.currents[1]
            + self.currents[2] * self.currents[2])
            .sqrt();
        let _ = emf;
        self.kt * i_mag
    }
}
/// Spring-damper actuator.
pub struct SpringDamper {
    /// Spring stiffness coefficient.
    pub k: f64,
    /// Damping coefficient.
    pub c: f64,
    /// Natural (rest) length of the spring.
    pub rest_length: f64,
}
impl SpringDamper {
    /// Creates a new `SpringDamper`.
    pub fn new(k: f64, c: f64, rest_length: f64) -> Self {
        Self { k, c, rest_length }
    }
    /// Computes the spring-damper force.
    pub fn force(&self, length: f64, velocity: f64) -> f64 {
        self.k * (self.rest_length - length) - self.c * velocity
    }
    /// Computes the elastic potential energy stored in the spring.
    pub fn energy(&self, length: f64) -> f64 {
        let dx = length - self.rest_length;
        0.5 * self.k * dx * dx
    }
    /// Natural frequency (rad/s) for a given attached mass.
    pub fn natural_frequency(&self, mass: f64) -> f64 {
        if mass > 1e-15 {
            (self.k / mass).sqrt()
        } else {
            0.0
        }
    }
    /// Damping ratio for a given attached mass.
    pub fn damping_ratio(&self, mass: f64) -> f64 {
        let omega_n = self.natural_frequency(mass);
        if omega_n > 1e-15 {
            self.c / (2.0 * mass * omega_n)
        } else {
            0.0
        }
    }
    /// Critical damping coefficient for a given mass.
    pub fn critical_damping(&self, mass: f64) -> f64 {
        2.0 * (self.k * mass).sqrt()
    }
    /// Power dissipated by the damper.
    pub fn damping_power(&self, velocity: f64) -> f64 {
        self.c * velocity * velocity
    }
}
/// Simple lumped-parameter thermal model for a motor.
///
/// Models the motor as a single thermal mass with convective cooling
/// to the ambient environment.
#[derive(Debug, Clone)]
pub struct MotorThermalModel {
    /// Current temperature (C).
    pub temperature: f64,
    /// Maximum allowed temperature (C).
    pub max_temperature: f64,
    /// Thermal mass (J/C) - energy required to raise temp by 1 degree.
    pub thermal_mass: f64,
    /// Thermal resistance to ambient (C/W) - resistance to heat dissipation.
    pub thermal_resistance: f64,
    /// Ambient temperature (C).
    pub ambient_temperature: f64,
}
impl MotorThermalModel {
    /// Create a new thermal model.
    pub fn new(
        ambient_temperature: f64,
        max_temperature: f64,
        thermal_mass: f64,
        thermal_resistance: f64,
    ) -> Self {
        Self {
            temperature: ambient_temperature,
            max_temperature,
            thermal_mass,
            thermal_resistance,
            ambient_temperature,
        }
    }
    /// Update temperature for one time step.
    ///
    /// `heat_input` is the power dissipated as heat (W).
    pub fn step(&mut self, heat_input: f64, dt: f64) {
        let cooling = (self.temperature - self.ambient_temperature) / self.thermal_resistance;
        let net_power = heat_input - cooling;
        let d_temp = net_power * dt / self.thermal_mass;
        self.temperature += d_temp;
    }
    /// Whether the temperature exceeds the maximum.
    pub fn is_overheating(&self) -> bool {
        self.temperature > self.max_temperature
    }
    /// Temperature margin before overheating.
    pub fn margin(&self) -> f64 {
        self.max_temperature - self.temperature
    }
    /// Thermal time constant (seconds).
    pub fn time_constant(&self) -> f64 {
        self.thermal_mass * self.thermal_resistance
    }
    /// Steady-state temperature for a given constant heat input.
    pub fn steady_state_temperature(&self, heat_input: f64) -> f64 {
        self.ambient_temperature + heat_input * self.thermal_resistance
    }
    /// Reset temperature to ambient.
    pub fn reset(&mut self) {
        self.temperature = self.ambient_temperature;
    }
}
