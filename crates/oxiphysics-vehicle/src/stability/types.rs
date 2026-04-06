//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Anti-lock Braking System controller (legacy high-level wrapper).
///
/// Prevents wheel lockup during hard braking by modulating brake torque.
pub struct AbsController {
    /// Slip ratio threshold for ABS activation (typical: 0.1-0.2).
    pub slip_threshold: f64,
    /// Brake torque reduction factor per ABS cycle (typical: 0.8).
    pub reduction_factor: f64,
}
impl AbsController {
    /// Create a new ABS controller with the given slip threshold and reduction factor.
    pub fn new(slip_threshold: f64, reduction_factor: f64) -> Self {
        Self {
            slip_threshold,
            reduction_factor,
        }
    }
    /// Given current slip ratio and requested brake torque, return actual brake torque.
    ///
    /// If `|slip| > threshold`: reduce torque by `reduction_factor`.
    /// Otherwise return the requested torque unchanged.
    pub fn apply(&self, slip_ratio: f64, requested_torque: f64) -> f64 {
        if slip_ratio.abs() > self.slip_threshold {
            requested_torque * self.reduction_factor
        } else {
            requested_torque
        }
    }
    /// Returns `true` when ABS is currently active (slip exceeds threshold).
    pub fn is_active(&self, slip_ratio: f64) -> bool {
        slip_ratio.abs() > self.slip_threshold
    }
}
/// Rollover warning system based on lateral acceleration and roll angle.
///
/// Triggers when either the lateral G-force exceeds `g_threshold` **or** the
/// roll angle exceeds `roll_threshold_rad`.
#[derive(Debug, Clone)]
pub struct RolloverDetection {
    /// Lateral acceleration threshold in m/s² (positive).
    pub g_threshold: f64,
    /// Roll angle threshold in radians.
    pub roll_threshold_rad: f64,
}
impl RolloverDetection {
    /// Create a detection module with typical passenger-car defaults.
    ///
    /// - g threshold: 8.5 m/s² (~0.87 g)
    /// - roll threshold: 0.35 rad (~20°)
    pub fn default_car() -> Self {
        Self {
            g_threshold: 8.5,
            roll_threshold_rad: 0.35,
        }
    }
    /// Returns `true` when the vehicle is at risk of rollover.
    pub fn is_at_risk(&self, lateral_accel: f64, roll_rad: f64) -> bool {
        lateral_accel.abs() > self.g_threshold || roll_rad.abs() > self.roll_threshold_rad
    }
    /// Risk level in `[0, 1]`; 0 = no risk, 1 = at or beyond rollover threshold.
    pub fn risk_level(&self, lateral_accel: f64, roll_rad: f64) -> f64 {
        let g_ratio = (lateral_accel.abs() / self.g_threshold).min(1.0);
        let r_ratio = (roll_rad.abs() / self.roll_threshold_rad).min(1.0);
        g_ratio.max(r_ratio)
    }
}
impl RolloverDetection {
    /// Compute the Time-to-Limit (TTL) rollover metric (seconds).
    ///
    /// TTL estimates the time remaining before the lateral acceleration reaches
    /// the rollover threshold, assuming a constant acceleration build-up rate.
    ///
    /// `TTL = (g_threshold - |lat_accel|) / |accel_rate|`
    ///
    /// Returns `f64::INFINITY` when the vehicle is not building lateral G, or
    /// when it is already within `headroom_fraction` of the threshold.
    ///
    /// # Arguments
    /// * `lat_accel`    – current lateral acceleration (m/s²)
    /// * `accel_rate`   – rate of change of lateral acceleration (m/s³)
    /// * `roll_rad`     – current roll angle (radians)
    pub fn compute_ttl(&self, lat_accel: f64, accel_rate: f64, roll_rad: f64) -> f64 {
        let g_margin = self.g_threshold - lat_accel.abs();
        let roll_margin = self.roll_threshold_rad - roll_rad.abs();
        let margin = g_margin.min(roll_margin);
        if margin <= 0.0 {
            return 0.0;
        }
        if accel_rate.abs() < 1e-9 {
            return f64::INFINITY;
        }
        if accel_rate * lat_accel.signum() <= 0.0 {
            return f64::INFINITY;
        }
        (margin / accel_rate.abs()).max(0.0)
    }
    /// Risk index in `[0, 1]` combining lateral G and roll, scaled to threshold.
    pub fn combined_risk(&self, lat_accel: f64, roll_rad: f64) -> f64 {
        let g_ratio = (lat_accel.abs() / self.g_threshold).min(1.0);
        let r_ratio = (roll_rad.abs() / self.roll_threshold_rad).min(1.0);
        (0.7 * g_ratio + 0.3 * r_ratio).min(1.0)
    }
}
/// Estimated vehicle attitude: pitch and roll angles.
#[derive(Debug, Clone, Default)]
pub struct PitchRollState {
    /// Roll angle in radians (positive = left side up).
    pub roll: f64,
    /// Pitch angle in radians (positive = nose up).
    pub pitch: f64,
    /// Roll rate in rad/s.
    pub roll_rate: f64,
    /// Pitch rate in rad/s.
    pub pitch_rate: f64,
}
impl PitchRollState {
    /// Integrate the state using Euler forward integration with `dt` seconds.
    ///
    /// `roll_dot` and `pitch_dot` are gyroscope readings (rad/s).
    pub fn integrate(&mut self, roll_dot: f64, pitch_dot: f64, dt: f64) {
        self.roll_rate = roll_dot;
        self.pitch_rate = pitch_dot;
        self.roll += roll_dot * dt;
        self.pitch += pitch_dot * dt;
    }
    /// Correct the attitude estimate using accelerometer readings.
    ///
    /// The accelerometer provides a gravity-referenced roll/pitch.  A
    /// complementary filter blends gyroscope integration with the
    /// accelerometer correction using `alpha` as the high-pass weight
    /// (`alpha` close to 1 → trust gyro more; close to 0 → trust accel).
    pub fn correct_with_accel(&mut self, accel_roll: f64, accel_pitch: f64, alpha: f64) {
        let alpha = alpha.clamp(0.0, 1.0);
        self.roll = alpha * self.roll + (1.0 - alpha) * accel_roll;
        self.pitch = alpha * self.pitch + (1.0 - alpha) * accel_pitch;
    }
}
/// Simplified bicycle model for lateral vehicle dynamics.
///
/// State variables: sideslip angle β (rad) and yaw rate r (rad/s).
#[derive(Debug, Clone)]
pub struct BicycleModel {
    /// Vehicle mass (kg).
    pub mass: f64,
    /// Front cornering stiffness (N/rad).
    pub cf: f64,
    /// Rear cornering stiffness (N/rad).
    pub cr: f64,
    /// Distance from CG to front axle (m).
    pub a: f64,
    /// Distance from CG to rear axle (m).
    pub b: f64,
    /// Yaw moment of inertia (kg·m²).
    pub iz: f64,
}
impl BicycleModel {
    /// Create a new bicycle model. Yaw inertia is estimated as `Iz = m·(a²+b²)/3`.
    pub fn new(mass: f64, cf: f64, cr: f64, a: f64, b: f64) -> Self {
        let iz = mass * (a * a + b * b) / 3.0;
        Self {
            mass,
            cf,
            cr,
            a,
            b,
            iz,
        }
    }
    /// Compute front and rear lateral tire forces.
    ///
    /// Slip angles (linear model):
    /// - `alpha_f = delta - beta - a·r/v`
    /// - `alpha_r = -beta + b·r/v`
    ///
    /// Returns `(Fyf, Fyr)`.
    pub fn lateral_forces(&self, speed: f64, beta: f64, r: f64, delta: f64) -> (f64, f64) {
        if speed.abs() < 1e-6 {
            return (0.0, 0.0);
        }
        let alpha_f = delta - beta - self.a * r / speed;
        let alpha_r = -beta + self.b * r / speed;
        let fyf = self.cf * alpha_f;
        let fyr = self.cr * alpha_r;
        (fyf, fyr)
    }
    /// Integrate the bicycle model by one time step `dt` using Euler integration.
    ///
    /// Equations of motion:
    /// - `β̇ = (Fyf + Fyr) / (m·v) - r`
    /// - `ṙ  = (a·Fyf - b·Fyr) / Iz`
    ///
    /// Returns `(beta_new, r_new)`.
    pub fn step(&self, speed: f64, beta: f64, r: f64, delta: f64, dt: f64) -> (f64, f64) {
        let (fyf, fyr) = self.lateral_forces(speed, beta, r, delta);
        let beta_dot = if speed.abs() > 1e-6 {
            (fyf + fyr) / (self.mass * speed) - r
        } else {
            0.0
        };
        let r_dot = (self.a * fyf - self.b * fyr) / self.iz;
        let beta_new = beta + beta_dot * dt;
        let r_new = r + r_dot * dt;
        (beta_new, r_new)
    }
}
/// Electronic Brake-force Distribution configuration.
///
/// Adjusts the front/rear brake bias based on vehicle load and deceleration
/// to maintain optimal brake efficiency and stability.
#[derive(Debug, Clone)]
pub struct EbdConfig {
    /// Base front brake bias \[0, 1\]. Typical: 0.65 for front-heavy car.
    pub front_bias: f64,
    /// Additional front bias per g of deceleration (increases with load transfer).
    pub decel_bias_gain: f64,
    /// Minimum rear brake fraction (to prevent complete cut-off).
    pub min_rear_fraction: f64,
    /// Maximum front brake fraction.
    pub max_front_fraction: f64,
}
/// Runtime state for the Electronic Stability Control.
#[derive(Debug, Clone)]
pub struct EscState {
    /// Whether ESC is currently active.
    pub is_active: bool,
    /// Human-readable label for the current intervention type.
    pub intervention_type: String,
    /// Signed yaw-rate error (actual − target) in rad/s.
    pub yaw_rate_error: f64,
}
/// Vehicle handling balance based on front and rear slip angle difference.
///
/// `delta_alpha = alpha_front - alpha_rear`
/// - Positive → understeer (front slides more).
/// - Negative → oversteer (rear slides more).
#[derive(Debug, Clone, PartialEq)]
pub enum HandlingBalance {
    /// |Δα| is within the neutral threshold.
    Neutral,
    /// Front axle slip angle exceeds rear (understeer).
    Understeer,
    /// Rear axle slip angle exceeds front (oversteer).
    Oversteer,
}
/// Enhanced ABS controller that also monitors wheel angular deceleration.
///
/// Supplements the basic slip ratio trigger with a wheel deceleration check:
/// if `|d(omega)/dt| > decel_threshold`, ABS activates even before slip builds.
#[derive(Debug, Clone)]
pub struct AbsEnhanced {
    /// Slip ratio threshold for activation.
    pub slip_threshold: f64,
    /// Wheel angular deceleration threshold (rad/s²).
    pub decel_threshold: f64,
    /// Brake pressure reduction factor per step.
    pub pressure_reduction: f64,
    /// Brake pressure recovery rate per second.
    pub pressure_recovery_rate: f64,
    /// Current brake pressures per wheel \[FL, FR, RL, RR\].
    pub pressures: [f64; 4],
    /// Previous wheel angular velocities (rad/s) \[FL, FR, RL, RR\].
    pub prev_wheel_speeds: [f64; 4],
    /// Whether ABS is active on each wheel.
    pub active: [bool; 4],
}
impl AbsEnhanced {
    /// Create an enhanced ABS controller.
    pub fn new(slip_threshold: f64, decel_threshold: f64) -> Self {
        Self {
            slip_threshold,
            decel_threshold,
            pressure_reduction: 0.85,
            pressure_recovery_rate: 0.5,
            pressures: [1.0; 4],
            prev_wheel_speeds: [0.0; 4],
            active: [false; 4],
        }
    }
    /// Default ABS for a passenger car.
    pub fn default_car() -> Self {
        Self::new(0.12, 50.0)
    }
    /// Update ABS state and return effective brake pressures \[FL, FR, RL, RR\].
    ///
    /// `slip_ratios` – per-wheel longitudinal slip ratios.
    /// `wheel_speeds_rad_s` – per-wheel angular velocities (rad/s).
    /// `dt` – timestep (s).
    pub fn update(
        &mut self,
        slip_ratios: &[f64; 4],
        wheel_speeds_rad_s: &[f64; 4],
        dt: f64,
    ) -> [f64; 4] {
        let mut result = [0.0f64; 4];
        for i in 0..4 {
            let decel = if dt > 1e-12 {
                (self.prev_wheel_speeds[i] - wheel_speeds_rad_s[i]) / dt
            } else {
                0.0
            };
            let slip_exceeded = slip_ratios[i].abs() > self.slip_threshold;
            let decel_exceeded = decel > self.decel_threshold;
            if slip_exceeded || decel_exceeded {
                self.active[i] = true;
                self.pressures[i] = (self.pressures[i] * self.pressure_reduction).max(0.0);
            } else {
                self.active[i] = false;
                self.pressures[i] = (self.pressures[i] + self.pressure_recovery_rate * dt).min(1.0);
            }
            self.prev_wheel_speeds[i] = wheel_speeds_rad_s[i];
            result[i] = self.pressures[i];
        }
        result
    }
    /// Whether any wheel has ABS active.
    pub fn any_active(&self) -> bool {
        self.active.iter().any(|&a| a)
    }
}
/// Configuration for Active Roll Control (ARC).
#[derive(Debug, Clone)]
pub struct ArcConfig {
    /// Proportional gain for roll angle correction (N·m/rad).
    pub roll_gain: f64,
    /// Derivative gain for roll rate correction (N·m·s/rad).
    pub roll_rate_gain: f64,
    /// Maximum anti-roll bar torque (N·m).
    pub max_torque: f64,
    /// Front axle fraction of total torque (0–1).
    pub front_fraction: f64,
}
/// Hill hold control: retains brake pressure on a gradient until the vehicle
/// begins to move forward (throttle applied and clutch engaged).
#[derive(Debug, Clone)]
pub struct HillHold {
    /// Gradient threshold (degrees) above which hill hold activates.
    pub slope_threshold_deg: f64,
    /// Retained brake pressure fraction when active \[0, 1\].
    pub hold_pressure: f64,
    /// Whether hill hold is currently engaged.
    pub is_engaged: bool,
}
impl HillHold {
    /// Create a hill hold system.
    pub fn new(slope_threshold_deg: f64, hold_pressure: f64) -> Self {
        Self {
            slope_threshold_deg,
            hold_pressure: hold_pressure.clamp(0.0, 1.0),
            is_engaged: false,
        }
    }
    /// Standard hill hold settings.
    pub fn default_car() -> Self {
        Self::new(3.0, 0.3)
    }
    /// Update hill hold state.
    ///
    /// * `slope_deg` – current road slope in degrees (positive = uphill).
    /// * `vehicle_speed` – current vehicle speed (m/s).
    /// * `throttle` – driver throttle input \[0, 1\].
    ///
    /// Returns the effective brake pressure to apply \[0, 1\].
    pub fn update(&mut self, slope_deg: f64, vehicle_speed: f64, throttle: f64) -> f64 {
        let on_slope = slope_deg.abs() > self.slope_threshold_deg;
        let stationary = vehicle_speed.abs() < 0.1;
        let applying_throttle = throttle > 0.05;
        if on_slope && stationary && !applying_throttle {
            self.is_engaged = true;
        } else if applying_throttle || vehicle_speed.abs() > 0.5 {
            self.is_engaged = false;
        }
        if self.is_engaged {
            self.hold_pressure
        } else {
            0.0
        }
    }
}
/// Luenberger observer state for estimating yaw rate and sideslip angle.
///
/// State: x = \[beta (sideslip, rad), r (yaw rate, rad/s)\]
/// Measurement: y = r (from yaw sensor)
#[derive(Debug, Clone, Default)]
pub struct YawRateObserver {
    /// Estimated sideslip angle β (rad).
    pub beta_est: f64,
    /// Estimated yaw rate r (rad/s).
    pub r_est: f64,
    /// Observer gain L1 for sideslip correction.
    pub l1: f64,
    /// Observer gain L2 for yaw rate correction.
    pub l2: f64,
}
impl YawRateObserver {
    /// Create a new observer with given Luenberger gains.
    #[allow(dead_code)]
    pub fn new(l1: f64, l2: f64) -> Self {
        Self {
            beta_est: 0.0,
            r_est: 0.0,
            l1,
            l2,
        }
    }
    /// Propagate the observer one step using bicycle model dynamics.
    ///
    /// `bm`        – bicycle model (provides A matrix)
    /// `speed`     – longitudinal velocity (m/s)
    /// `delta`     – front wheel steer angle (rad)
    /// `r_meas`    – measured yaw rate (rad/s)
    /// `dt`        – time step (s)
    #[allow(dead_code)]
    pub fn update(&mut self, bm: &BicycleModel, speed: f64, delta: f64, r_meas: f64, dt: f64) {
        if speed.abs() < 1e-6 {
            return;
        }
        let m = bm.mass;
        let cf = bm.cf;
        let cr = bm.cr;
        let a = bm.a;
        let b = bm.b;
        let iz = bm.iz;
        let v = speed;
        let a11 = -(cf + cr) / (m * v);
        let a12 = (cr * b - cf * a) / (m * v * v) - 1.0;
        let a21 = (cr * b - cf * a) / iz;
        let a22 = -(cf * a * a + cr * b * b) / (iz * v);
        let b1 = cf / (m * v);
        let b2 = cf * a / iz;
        let err = r_meas - self.r_est;
        let beta_dot = a11 * self.beta_est + a12 * self.r_est + b1 * delta + self.l1 * err;
        let r_dot = a21 * self.beta_est + a22 * self.r_est + b2 * delta + self.l2 * err;
        self.beta_est += beta_dot * dt;
        self.r_est += r_dot * dt;
    }
}
/// ESC yaw moment demand calculator.
///
/// Computes the corrective yaw moment `Mz` (N·m) required to bring the vehicle
/// yaw rate toward the target, considering sideslip angle.
#[derive(Debug, Clone)]
pub struct EscYawMoment {
    /// Proportional gain on yaw-rate error (N·m / (rad/s)).
    pub kp_yaw: f64,
    /// Proportional gain on sideslip angle error (N·m / deg).
    pub kp_sideslip: f64,
    /// Maximum allowable yaw moment magnitude (N·m).
    pub max_moment: f64,
    /// Vehicle yaw moment of inertia (kg·m²), used for normalisation.
    pub iz: f64,
}
impl EscYawMoment {
    /// Typical mid-size sedan ESC parameters.
    pub fn default_sedan() -> Self {
        Self {
            kp_yaw: 3500.0,
            kp_sideslip: 200.0,
            max_moment: 2000.0,
            iz: 2500.0,
        }
    }
    /// Compute the corrective yaw moment demand (N·m).
    ///
    /// `Mz = kp_yaw * (r_actual - r_target) + kp_sideslip * beta_deg`
    ///
    /// The sign convention: positive Mz = counter-clockwise (corrects oversteer
    /// to the right / understeer to the left).
    ///
    /// Result is clamped to `[-max_moment, +max_moment]`.
    pub fn compute_yaw_moment_demand(&self, r_actual: f64, r_target: f64, beta_deg: f64) -> f64 {
        let yaw_error = r_actual - r_target;
        let mz = self.kp_yaw * yaw_error + self.kp_sideslip * beta_deg;
        mz.clamp(-self.max_moment, self.max_moment)
    }
    /// Normalised yaw acceleration demand (rad/s²) from the moment.
    pub fn yaw_acceleration_demand(&self, r_actual: f64, r_target: f64, beta_deg: f64) -> f64 {
        let mz = self.compute_yaw_moment_demand(r_actual, r_target, beta_deg);
        mz / self.iz.max(1.0)
    }
}
/// State for an advanced ABS pressure modulation controller.
#[derive(Debug, Clone)]
pub struct AbsPressureModulator {
    /// Slip ratio threshold above which ABS applies pressure reduction.
    pub slip_threshold: f64,
    /// Re-apply threshold: slip ratio below which ABS allows pressure to increase.
    pub regen_threshold: f64,
    /// Pressure increase rate (fraction per second) when below regen threshold.
    pub pressure_build_rate: f64,
    /// Pressure decrease rate (fraction per second) when above slip threshold.
    pub pressure_dump_rate: f64,
    /// Current brake pressure fraction `[0, 1]` for this wheel.
    pub pressure: f64,
}
impl AbsPressureModulator {
    /// Construct a modulator with typical passenger-car ABS parameters.
    pub fn default_car() -> Self {
        Self {
            slip_threshold: 0.15,
            regen_threshold: 0.08,
            pressure_build_rate: 5.0,
            pressure_dump_rate: 20.0,
            pressure: 1.0,
        }
    }
    /// Modulate brake pressure based on current wheel slip and time step.
    ///
    /// ABS logic:
    /// - If `|slip| > slip_threshold` → dump pressure at `pressure_dump_rate`.
    /// - If `|slip| < regen_threshold` → build pressure at `pressure_build_rate`.
    /// - Otherwise → hold current pressure.
    ///
    /// Returns the updated pressure fraction `[0, 1]`.
    pub fn modulate_brake_pressure(&mut self, slip: f64, dt: f64) -> f64 {
        if slip.abs() > self.slip_threshold {
            self.pressure = (self.pressure - self.pressure_dump_rate * dt).max(0.0);
        } else if slip.abs() < self.regen_threshold {
            self.pressure = (self.pressure + self.pressure_build_rate * dt).min(1.0);
        }
        self.pressure
    }
    /// Returns `true` when ABS intervention is active (pressure below 1.0).
    pub fn is_active(&self) -> bool {
        self.pressure < 1.0 - 1e-9
    }
}
/// Torque vectoring output: additional torque on each rear wheel (N·m).
#[derive(Debug, Clone, Default)]
pub struct TorqueVectoringOutput {
    /// Additional torque on rear-left wheel (positive = drive).
    pub rear_left: f64,
    /// Additional torque on rear-right wheel (positive = drive).
    pub rear_right: f64,
}
/// Runtime state for the ABS controller.
#[derive(Debug, Clone)]
pub struct AbsState {
    /// Whether ABS is currently active on any wheel.
    pub is_active: bool,
    /// Hydraulic brake pressure fraction per wheel `[FL, FR, RL, RR]` in \[0, 1\].
    pub pressure: Vec<f64>,
    /// Current phase label: `"increase"` or `"decrease"`.
    pub cycle_phase: String,
}
/// Electronic Stability Control (legacy high-level wrapper).
///
/// Applies individual wheel braking to counteract oversteer/understeer.
pub struct ElectronicStabilityControl {
    /// Yaw-rate error threshold (rad/s) above which ESC intervenes.
    pub yaw_rate_threshold: f64,
    /// Embedded ABS controller used when applying corrective braking.
    pub abs: AbsController,
}
impl ElectronicStabilityControl {
    /// Create an ESC with the given yaw-rate threshold.
    /// Defaults to a standard ABS sub-controller (threshold 0.15, reduction 0.8).
    pub fn new(yaw_rate_threshold: f64) -> Self {
        Self {
            yaw_rate_threshold,
            abs: AbsController::new(0.15, 0.8),
        }
    }
    /// Given measured vs target yaw rate, return corrective brake torque for each wheel.
    ///
    /// Returns `(front_left, front_right, rear_left, rear_right)` brake torque corrections.
    ///
    /// * Oversteer  (`measured > target`): brake the front-left and rear-left wheels.
    /// * Understeer (`measured < target`): brake the front-right and rear-right wheels.
    /// * Within threshold: no correction (all zeros).
    pub fn correction_torques(
        &self,
        measured_yaw: f64,
        target_yaw: f64,
        base_torque: f64,
    ) -> (f64, f64, f64, f64) {
        let error = measured_yaw - target_yaw;
        if error.abs() <= self.yaw_rate_threshold {
            return (0.0, 0.0, 0.0, 0.0);
        }
        let correction = base_torque * (error.abs() / (error.abs() + 1.0));
        if error > 0.0 {
            (correction, 0.0, correction, 0.0)
        } else {
            (0.0, correction, 0.0, correction)
        }
    }
}
/// Per-wheel ESC brake and throttle intervention output.
#[derive(Debug, Clone, Default)]
pub struct EscIntervention {
    /// Additional brake pressure fraction applied to front-left wheel \[0, 1\].
    pub brake_front_left: f64,
    /// Additional brake pressure fraction applied to front-right wheel \[0, 1\].
    pub brake_front_right: f64,
    /// Additional brake pressure fraction applied to rear-left wheel \[0, 1\].
    pub brake_rear_left: f64,
    /// Additional brake pressure fraction applied to rear-right wheel \[0, 1\].
    pub brake_rear_right: f64,
    /// Throttle reduction fraction \[0, 1\]; 0 = no reduction, 1 = full cut.
    pub throttle_reduction: f64,
}
/// Configuration parameters for the Traction Control System.
#[derive(Debug, Clone)]
pub struct TcsConfig {
    /// Wheel slip ratio above which TCS intervenes (default 0.1).
    pub slip_threshold: f64,
    /// Multiplicative torque reduction applied each step while slipping (default 0.8).
    pub torque_reduction_rate: f64,
    /// Minimum allowed torque fraction (floor) when TCS is active (default 0.2).
    pub min_torque_fraction: f64,
}
/// Unified stability controller combining TCS, ABS, and ESC.
///
/// Acts as the top-level interface: given vehicle state, it runs all three
/// sub-controllers and returns combined intervention commands.
pub struct StabilityController {
    /// Traction control sub-controller.
    pub tcs_config: TcsConfig,
    /// ABS sub-controller configuration.
    pub abs_config: AbsConfig,
    /// ESC configuration.
    pub esc_config: EscConfig,
    /// ESC runtime state.
    pub esc_state: EscState,
    /// ABS runtime state.
    pub abs_state: AbsState,
    /// TCS runtime state.
    pub tcs_state: TcsState,
    /// Rollover detection module.
    pub rollover: RolloverDetection,
}
impl StabilityController {
    /// Create a controller with default sub-system configurations.
    pub fn new() -> Self {
        Self {
            tcs_config: TcsConfig::default(),
            abs_config: AbsConfig::default(),
            esc_config: EscConfig::default(),
            esc_state: EscState::default(),
            abs_state: AbsState::default(),
            tcs_state: TcsState::default(),
            rollover: RolloverDetection::default_car(),
        }
    }
    /// Run one update step.
    ///
    /// # Arguments
    /// * `driven_slip`   – slip ratio of the driven wheel(s)
    /// * `wheel_slips`   – per-wheel slip ratios `[FL, FR, RL, RR]`
    /// * `drive_torque`  – requested drive torque (N·m)
    /// * `brake_torque`  – requested brake torque (N·m, positive = braking)
    /// * `actual_yaw`    – measured yaw rate (rad/s)
    /// * `target_yaw`    – desired yaw rate (rad/s)
    /// * `sideslip_deg`  – vehicle sideslip angle (degrees)
    /// * `dt`            – time step (s)
    ///
    /// Returns `(actual_drive_torque, brake_pressures, esc_intervention)`.
    #[allow(clippy::too_many_arguments)]
    pub fn update(
        &mut self,
        driven_slip: f64,
        wheel_slips: &[f64; 4],
        drive_torque: f64,
        _brake_torque: f64,
        actual_yaw: f64,
        target_yaw: f64,
        sideslip_deg: f64,
        dt: f64,
    ) -> (f64, [f64; 4], EscIntervention) {
        let tcs_factor = tcs_update(&mut self.tcs_state, driven_slip, &self.tcs_config);
        let actual_drive = drive_torque * tcs_factor;
        let brake_pressures = abs_update(&mut self.abs_state, wheel_slips, dt, &self.abs_config);
        let esc = esc_update(
            &mut self.esc_state,
            actual_yaw,
            target_yaw,
            sideslip_deg,
            &self.esc_config,
        );
        (actual_drive, brake_pressures, esc)
    }
    /// Returns `true` if any sub-system is currently active.
    pub fn any_active(&self) -> bool {
        self.tcs_state.is_active || self.abs_state.is_active || self.esc_state.is_active
    }
}
/// Runtime state for the Traction Control System.
#[derive(Debug, Clone)]
pub struct TcsState {
    /// Whether TCS is currently intervening.
    pub is_active: bool,
    /// Current torque reduction factor in \[0, 1\]; 1.0 = no reduction.
    pub torque_reduction_factor: f64,
    /// Measured slip ratio at the last update.
    pub slip_ratio: f64,
}
/// Traction control system (legacy high-level wrapper).
///
/// Reduces drive torque when driven wheels spin (slip > threshold).
pub struct TractionControl {
    /// Slip ratio threshold above which traction control intervenes.
    pub slip_threshold: f64,
    /// Drive torque reduction factor applied when slip exceeds the threshold.
    pub reduction_factor: f64,
}
impl TractionControl {
    /// Create a new traction control system.
    pub fn new(slip_threshold: f64, reduction_factor: f64) -> Self {
        Self {
            slip_threshold,
            reduction_factor,
        }
    }
    /// Given current wheel slip ratio and requested drive torque, return actual drive torque.
    ///
    /// When `|slip| > threshold` the torque is scaled by `reduction_factor`.
    pub fn apply(&self, slip_ratio: f64, drive_torque: f64) -> f64 {
        if slip_ratio.abs() > self.slip_threshold {
            drive_torque * self.reduction_factor
        } else {
            drive_torque
        }
    }
}
impl TractionControl {
    /// Compute the longitudinal slip ratio for a driven wheel.
    ///
    /// `slip = (omega_wheel * r - v_vehicle) / v_vehicle`
    ///
    /// where `omega_wheel` is the wheel angular velocity (rad/s), `r` is the
    /// effective rolling radius (m), and `v_vehicle` is the vehicle speed (m/s).
    ///
    /// Returns 0.0 when vehicle speed is near zero to avoid division errors.
    /// Result is clamped to `[-1, 1]`.
    pub fn compute_slip_ratio(omega_wheel: f64, wheel_radius: f64, v_vehicle: f64) -> f64 {
        if v_vehicle.abs() < 1e-6 {
            return 0.0;
        }
        let v_wheel = omega_wheel * wheel_radius;
        ((v_wheel - v_vehicle) / v_vehicle).clamp(-1.0, 1.0)
    }
}
/// Rollover prevention system based on lateral acceleration.
///
/// When lateral G exceeds the threshold, throttle is cut and a warning
/// brake intervention is requested.
#[derive(Debug, Clone)]
pub struct RolloverPrevention {
    /// Lateral acceleration threshold (m/s²) above which to intervene.
    pub lat_accel_threshold: f64,
    /// Throttle cut fraction applied during intervention \[0, 1\].
    pub throttle_cut: f64,
    /// Brake fraction applied during intervention \[0, 1\].
    pub brake_fraction: f64,
}
impl RolloverPrevention {
    /// Create a rollover prevention system.
    pub fn new(lat_accel_threshold: f64, throttle_cut: f64, brake_fraction: f64) -> Self {
        Self {
            lat_accel_threshold: lat_accel_threshold.max(0.1),
            throttle_cut: throttle_cut.clamp(0.0, 1.0),
            brake_fraction: brake_fraction.clamp(0.0, 1.0),
        }
    }
    /// Default SUV rollover prevention (threshold ~7 m/s² ≈ 0.71 g).
    pub fn default_suv() -> Self {
        Self::new(7.0, 0.8, 0.15)
    }
    /// Evaluate the lateral acceleration and return `(throttle_reduction, brake_pressure)`.
    ///
    /// Returns `(0.0, 0.0)` when below threshold.
    pub fn evaluate(&self, lateral_accel: f64) -> (f64, f64) {
        if lateral_accel.abs() > self.lat_accel_threshold {
            (self.throttle_cut, self.brake_fraction)
        } else {
            (0.0, 0.0)
        }
    }
    /// Whether the rollover prevention is currently active.
    pub fn is_active(&self, lateral_accel: f64) -> bool {
        lateral_accel.abs() > self.lat_accel_threshold
    }
}
/// Represents the vehicle's position within its stability envelope.
#[derive(Debug, Clone, PartialEq)]
pub enum StabilityEnvelopeStatus {
    /// Normal: all margins are comfortable.
    Normal,
    /// Warning: approaching one or more limits.
    Warning,
    /// Critical: near or beyond stability limits, intervention required.
    Critical,
}
/// Vehicle stability envelope monitor.
///
/// Aggregates multiple stability metrics to give an overall status.
#[derive(Debug, Clone)]
pub struct StabilityEnvelopeMonitor {
    /// SSF threshold below which rollover risk is elevated.
    pub ssf_warning: f64,
    /// SSF threshold below which status is Critical.
    pub ssf_critical: f64,
    /// Slip ratio above which Warning.
    pub slip_warning: f64,
    /// Slip ratio above which Critical.
    pub slip_critical: f64,
    /// Sideslip angle (deg) above which Warning.
    pub sideslip_warning_deg: f64,
    /// Sideslip angle (deg) above which Critical.
    pub sideslip_critical_deg: f64,
}
impl StabilityEnvelopeMonitor {
    /// Evaluate the stability envelope given current measurements.
    #[allow(dead_code)]
    pub fn evaluate(&self, ssf: f64, max_slip: f64, sideslip_deg: f64) -> StabilityEnvelopeStatus {
        let slip = max_slip.abs();
        let ss = sideslip_deg.abs();
        if ssf < self.ssf_critical || slip > self.slip_critical || ss > self.sideslip_critical_deg {
            return StabilityEnvelopeStatus::Critical;
        }
        if ssf < self.ssf_warning || slip > self.slip_warning || ss > self.sideslip_warning_deg {
            return StabilityEnvelopeStatus::Warning;
        }
        StabilityEnvelopeStatus::Normal
    }
    /// Return a normalised stability score in \[0, 1\] where 1 = perfectly stable.
    #[allow(dead_code)]
    pub fn stability_score(&self, ssf: f64, max_slip: f64, sideslip_deg: f64) -> f64 {
        let slip = max_slip.abs();
        let ss = sideslip_deg.abs();
        let ssf_score =
            ((ssf - self.ssf_critical) / (self.ssf_warning - self.ssf_critical)).clamp(0.0, 1.0);
        let slip_score = 1.0 - (slip / self.slip_critical).clamp(0.0, 1.0);
        let ss_score = 1.0 - (ss / self.sideslip_critical_deg).clamp(0.0, 1.0);
        (ssf_score + slip_score + ss_score) / 3.0
    }
}
/// Anti-roll bar torque output for front and rear axles.
#[derive(Debug, Clone, Default)]
pub struct ArcOutput {
    /// Torque applied to front anti-roll bar (N·m).
    pub front_torque: f64,
    /// Torque applied to rear anti-roll bar (N·m).
    pub rear_torque: f64,
}
/// Emergency Brake Assist state.
///
/// Detects a fast pedal application (panic stop) and holds or increases
/// brake pressure to achieve maximum deceleration.
#[derive(Debug, Clone)]
pub struct EmergencyBrakeAssist {
    /// Pedal speed threshold (pressure/s) to trigger EBA.
    pub trigger_threshold: f64,
    /// Pressure amplification factor when active.
    pub amplification: f64,
    /// Whether EBA is currently active.
    pub is_active: bool,
    /// Previous brake pressure (used to compute speed of application).
    pub(super) prev_pressure: f64,
}
impl EmergencyBrakeAssist {
    /// Create a new EBA with given threshold and amplification.
    pub fn new(trigger_threshold: f64, amplification: f64) -> Self {
        Self {
            trigger_threshold: trigger_threshold.max(0.1),
            amplification: amplification.clamp(1.0, 4.0),
            is_active: false,
            prev_pressure: 0.0,
        }
    }
    /// Default EBA settings.
    pub fn default_car() -> Self {
        Self::new(5.0, 1.5)
    }
    /// Update EBA state and return the effective brake pressure.
    ///
    /// `requested_pressure` in \[0, 1\], `dt` in seconds.
    pub fn update(&mut self, requested_pressure: f64, dt: f64) -> f64 {
        let pressure = requested_pressure.clamp(0.0, 1.0);
        let dp_dt = if dt > 1e-12 {
            (pressure - self.prev_pressure) / dt
        } else {
            0.0
        };
        self.prev_pressure = pressure;
        if dp_dt > self.trigger_threshold {
            self.is_active = true;
        } else if pressure < 0.05 {
            self.is_active = false;
        }
        if self.is_active {
            (pressure * self.amplification).min(1.0)
        } else {
            pressure
        }
    }
}
/// Configuration for the ABS controller.
#[derive(Debug, Clone)]
pub struct AbsConfig {
    /// Target slip ratio (typically 0.1 – 0.15 for maximum braking).
    pub slip_target: f64,
    /// Hydraulic pressure increase rate (fraction per second) during recovery phase.
    pub pressure_increase_rate: f64,
    /// Hydraulic pressure decrease rate (fraction per second) during reduction phase.
    pub pressure_decrease_rate: f64,
}
/// Configuration parameters for the Electronic Stability Control.
#[derive(Debug, Clone)]
pub struct EscConfig {
    /// Yaw-rate error (rad/s) above which ESC intervenes.
    pub yaw_rate_threshold: f64,
    /// Vehicle sideslip angle (degrees) above which ESC intervenes.
    pub sideslip_threshold: f64,
    /// Gain applied to understeer correction (positive).
    pub understeer_gain: f64,
    /// Gain applied to oversteer correction (positive).
    pub oversteer_gain: f64,
}
/// ESP (Electronic Stability Program) intervention level.
#[derive(Debug, Clone, PartialEq)]
pub enum EspLevel {
    /// No intervention.
    None,
    /// Light intervention (gentle braking).
    Light,
    /// Heavy intervention (full braking + throttle cut).
    Heavy,
}
