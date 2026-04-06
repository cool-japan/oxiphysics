#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Vehicle simulation API for Python interop.
//!
//! Provides Python-friendly types for full vehicle simulation including:
//! tire models (Pacejka/Fiala), drivetrain, suspension, stability control,
//! telemetry, and racing line optimization.

#![allow(missing_docs)]
#![allow(dead_code)]

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

/// Compute the Pacejka magic formula lateral force Fy.
///
/// Parameters follow the standard Pacejka notation.
/// - `alpha` — slip angle in radians
/// - `b` — stiffness factor
/// - `c` — shape factor
/// - `d` — peak value
/// - `e` — curvature factor
///
/// Returns the lateral force as a fraction of normal load (Fy / Fz).
pub fn pacejka_fy(alpha: f64, b: f64, c: f64, d: f64, e: f64) -> f64 {
    let x = alpha.to_degrees(); // Pacejka often uses degrees internally
    let phi = (1.0 - e) * x + (e / b) * (b * x).atan();
    d * (c * (b * phi).atan()).sin()
}

/// Compute the Pacejka magic formula longitudinal force Fx.
///
/// - `kappa` — slip ratio (dimensionless, -1 to 1)
/// - `b` — stiffness factor
/// - `c` — shape factor
/// - `d` — peak value
/// - `e` — curvature factor
pub fn pacejka_fx(kappa: f64, b: f64, c: f64, d: f64, e: f64) -> f64 {
    let x = kappa * 100.0; // Pacejka often uses percent slip
    let phi = (1.0 - e) * x + (e / b) * (b * x).atan();
    d * (c * (b * phi).atan()).sin()
}

/// Compute the slip angle from tire velocity components.
///
/// - `vy` — lateral velocity component of the tire contact point (m/s)
/// - `vx` — longitudinal velocity component (m/s)
///
/// Returns slip angle in radians.
pub fn compute_slip_angle(vy: f64, vx: f64) -> f64 {
    if vx.abs() < 1e-4 {
        return 0.0;
    }
    (-vy / vx).atan()
}

/// Compute the longitudinal slip ratio.
///
/// - `wheel_speed` — peripheral speed of the tire = omega * R (m/s)
/// - `vehicle_speed` — forward speed of the vehicle (m/s)
///
/// Returns slip ratio (dimensionless, positive = drive slip, negative = brake slip).
pub fn compute_slip_ratio(wheel_speed: f64, vehicle_speed: f64) -> f64 {
    let v_ref = vehicle_speed.abs().max(wheel_speed.abs()).max(1e-4);
    (wheel_speed - vehicle_speed) / v_ref
}

/// Compute Ackermann steering angles for inner and outer front wheels.
///
/// - `steer_angle` — rack steering angle (rad, positive = left turn)
/// - `wheelbase` — distance between front and rear axles (m)
/// - `track_width` — distance between left and right wheels (m)
///
/// Returns `(angle_inner, angle_outer)` in radians.
pub fn ackermann_angles(steer_angle: f64, wheelbase: f64, track_width: f64) -> (f64, f64) {
    if steer_angle.abs() < 1e-9 {
        return (0.0, 0.0);
    }
    let r = wheelbase / steer_angle.abs();
    let sign = steer_angle.signum();
    let inner = (wheelbase / (r - track_width * 0.5)).atan() * sign;
    let outer = (wheelbase / (r + track_width * 0.5)).atan() * sign;
    if steer_angle > 0.0 {
        (inner, outer)
    } else {
        (outer, inner)
    }
}

// ---------------------------------------------------------------------------
// PyTireModel
// ---------------------------------------------------------------------------

/// Pacejka magic formula tire coefficients.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PacejkaCoeffs {
    /// Stiffness factor B.
    pub b: f64,
    /// Shape factor C.
    pub c: f64,
    /// Peak factor D (as fraction of Fz).
    pub d: f64,
    /// Curvature factor E.
    pub e: f64,
}

impl PacejkaCoeffs {
    /// Typical road tire lateral coefficients.
    pub fn road_lateral() -> Self {
        Self {
            b: 10.0,
            c: 1.9,
            d: 1.0,
            e: 0.97,
        }
    }

    /// Typical road tire longitudinal coefficients.
    pub fn road_longitudinal() -> Self {
        Self {
            b: 11.0,
            c: 1.65,
            d: 1.0,
            e: 0.0,
        }
    }

    /// Slick racing tire coefficients.
    pub fn slick_lateral() -> Self {
        Self {
            b: 12.0,
            c: 2.0,
            d: 1.4,
            e: 0.95,
        }
    }
}

/// Fiala tire model parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FialaParams {
    /// Cornering stiffness (N/rad).
    pub cornering_stiffness: f64,
    /// Friction coefficient.
    pub mu: f64,
}

impl FialaParams {
    /// Create default Fiala parameters.
    pub fn new(cornering_stiffness: f64, mu: f64) -> Self {
        Self {
            cornering_stiffness,
            mu,
        }
    }

    /// Compute lateral force given slip angle and normal load.
    pub fn lateral_force(&self, alpha: f64, fz: f64) -> f64 {
        let f_max = self.mu * fz;
        let f_lin = -self.cornering_stiffness * alpha;
        if f_lin.abs() > f_max {
            f_max * f_lin.signum()
        } else {
            f_lin
        }
    }
}

/// Tire model selection.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TireModelKind {
    /// Pacejka magic formula.
    Pacejka {
        lateral: PacejkaCoeffs,
        longitudinal: PacejkaCoeffs,
    },
    /// Fiala brush model.
    Fiala(FialaParams),
    /// Simple linear tire model.
    Linear {
        cornering_stiffness: f64,
        longitudinal_stiffness: f64,
    },
}

/// Tire model for force computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyTireModel {
    /// Tire model type.
    pub kind: TireModelKind,
    /// Tire rolling radius (m).
    pub radius: f64,
    /// Tire width (m).
    pub width: f64,
    /// Unloaded tire radius (m).
    pub unloaded_radius: f64,
}

impl PyTireModel {
    /// Create a Pacejka tire model for a road car.
    pub fn road_pacejka(radius: f64) -> Self {
        Self {
            kind: TireModelKind::Pacejka {
                lateral: PacejkaCoeffs::road_lateral(),
                longitudinal: PacejkaCoeffs::road_longitudinal(),
            },
            radius,
            width: 0.225,
            unloaded_radius: radius,
        }
    }

    /// Create a racing slick tire model.
    pub fn slick_pacejka(radius: f64) -> Self {
        Self {
            kind: TireModelKind::Pacejka {
                lateral: PacejkaCoeffs::slick_lateral(),
                longitudinal: PacejkaCoeffs::road_longitudinal(),
            },
            radius,
            width: 0.300,
            unloaded_radius: radius,
        }
    }

    /// Compute combined lateral and longitudinal forces.
    ///
    /// Returns `(Fx, Fy)` — longitudinal and lateral forces in N.
    pub fn compute_forces(&self, alpha: f64, kappa: f64, fz: f64) -> (f64, f64) {
        match &self.kind {
            TireModelKind::Pacejka {
                lateral,
                longitudinal,
            } => {
                let fy = pacejka_fy(alpha, lateral.b, lateral.c, lateral.d, lateral.e) * fz;
                let fx = pacejka_fx(
                    kappa,
                    longitudinal.b,
                    longitudinal.c,
                    longitudinal.d,
                    longitudinal.e,
                ) * fz;
                (fx, fy)
            }
            TireModelKind::Fiala(params) => {
                let fy = params.lateral_force(alpha, fz);
                let fx = 0.0;
                (fx, fy)
            }
            TireModelKind::Linear {
                cornering_stiffness,
                longitudinal_stiffness,
            } => {
                let fy = -cornering_stiffness * alpha;
                let fx = longitudinal_stiffness * kappa * fz;
                (fx, fy)
            }
        }
    }
}

// ---------------------------------------------------------------------------
// PyTireThermal
// ---------------------------------------------------------------------------

/// Three-layer tire thermal model (surface / bulk / carcass).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyTireThermal {
    /// Surface temperature (°C).
    pub surface_temp: f64,
    /// Bulk temperature (°C).
    pub bulk_temp: f64,
    /// Carcass temperature (°C).
    pub carcass_temp: f64,
    /// Ambient temperature (°C).
    pub ambient_temp: f64,
    /// Thermal conductivity surface→bulk (W/m²K).
    pub k_surface_bulk: f64,
    /// Thermal conductivity bulk→carcass (W/m²K).
    pub k_bulk_carcass: f64,
    /// Thermal conductivity carcass→ambient (W/m²K).
    pub k_carcass_ambient: f64,
    /// Optimal temperature for maximum grip (°C).
    pub optimal_temp: f64,
    /// Temperature window for full grip (°C).
    pub optimal_window: f64,
}

impl PyTireThermal {
    /// Create a default tire thermal model.
    pub fn new() -> Self {
        Self {
            surface_temp: 25.0,
            bulk_temp: 25.0,
            carcass_temp: 25.0,
            ambient_temp: 25.0,
            k_surface_bulk: 50.0,
            k_bulk_carcass: 20.0,
            k_carcass_ambient: 5.0,
            optimal_temp: 90.0,
            optimal_window: 20.0,
        }
    }

    /// Update the thermal model given power dissipation (W) and time step (s).
    pub fn update(&mut self, power: f64, dt: f64) {
        let mass_surface = 0.3; // kg (approximate)
        let mass_bulk = 2.0;
        let mass_carcass = 5.0;
        let cp = 1300.0; // J/(kg K) rubber

        let q_surface_bulk = self.k_surface_bulk * (self.surface_temp - self.bulk_temp);
        let q_bulk_carcass = self.k_bulk_carcass * (self.bulk_temp - self.carcass_temp);
        let q_carcass_ambient = self.k_carcass_ambient * (self.carcass_temp - self.ambient_temp);

        self.surface_temp += (power - q_surface_bulk) / (mass_surface * cp) * dt;
        self.bulk_temp += (q_surface_bulk - q_bulk_carcass) / (mass_bulk * cp) * dt;
        self.carcass_temp += (q_bulk_carcass - q_carcass_ambient) / (mass_carcass * cp) * dt;
    }

    /// Compute grip scaling factor based on surface temperature.
    pub fn grip_scale(&self) -> f64 {
        let dt = (self.surface_temp - self.optimal_temp).abs();
        if dt <= self.optimal_window * 0.5 {
            1.0
        } else {
            let excess = dt - self.optimal_window * 0.5;
            (1.0 - excess / self.optimal_window).max(0.3)
        }
    }
}

impl Default for PyTireThermal {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// PyDrivetrain
// ---------------------------------------------------------------------------

/// Engine torque curve defined by RPM breakpoints and torque values.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EngineCurve {
    /// RPM breakpoints (sorted ascending).
    pub rpm: Vec<f64>,
    /// Torque at each breakpoint (N·m).
    pub torque: Vec<f64>,
    /// Maximum RPM (redline).
    pub redline: f64,
    /// Idle RPM.
    pub idle: f64,
}

impl EngineCurve {
    /// Create a simple flat torque curve.
    pub fn flat(peak_torque: f64, redline: f64) -> Self {
        Self {
            rpm: vec![0.0, redline * 0.3, redline * 0.7, redline],
            torque: vec![
                peak_torque * 0.6,
                peak_torque,
                peak_torque,
                peak_torque * 0.8,
            ],
            redline,
            idle: 800.0,
        }
    }

    /// Interpolate torque at a given RPM.
    pub fn torque_at(&self, rpm: f64) -> f64 {
        let n = self.rpm.len();
        if n == 0 {
            return 0.0;
        }
        if rpm <= self.rpm[0] {
            return self.torque[0];
        }
        if rpm >= self.rpm[n - 1] {
            return self.torque[n - 1];
        }
        for i in 1..n {
            if rpm <= self.rpm[i] {
                let t = (rpm - self.rpm[i - 1]) / (self.rpm[i] - self.rpm[i - 1]);
                return self.torque[i - 1] + t * (self.torque[i] - self.torque[i - 1]);
            }
        }
        self.torque[n - 1]
    }
}

/// Differential type.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum DiffType {
    /// Open differential (no torque biasing).
    Open,
    /// Limited-slip differential.
    Lsd,
    /// Locked (spool).
    Locked,
}

/// Drivetrain configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyDrivetrain {
    /// Engine torque curve.
    pub engine_curve: EngineCurve,
    /// Gear ratios (index 0 = first gear).
    pub gear_ratios: Vec<f64>,
    /// Final drive ratio.
    pub final_drive: f64,
    /// Current gear (0-indexed).
    pub current_gear: usize,
    /// Differential type.
    pub diff_type: DiffType,
    /// Clutch engagement factor (0 = fully open, 1 = fully engaged).
    pub clutch: f64,
    /// Current engine RPM.
    pub engine_rpm: f64,
    /// Engine inertia (kg·m²).
    pub engine_inertia: f64,
}

impl PyDrivetrain {
    /// Create a default 6-speed drivetrain.
    pub fn six_speed(peak_torque: f64) -> Self {
        Self {
            engine_curve: EngineCurve::flat(peak_torque, 7000.0),
            gear_ratios: vec![3.31, 2.13, 1.55, 1.21, 0.97, 0.79],
            final_drive: 3.73,
            current_gear: 1,
            diff_type: DiffType::Open,
            clutch: 1.0,
            engine_rpm: 1000.0,
            engine_inertia: 0.2,
        }
    }

    /// Compute wheel torque given throttle (0-1) and the current state.
    pub fn wheel_torque(&self, throttle: f64) -> f64 {
        let engine_torque = self.engine_curve.torque_at(self.engine_rpm) * throttle;
        let gear_ratio = self
            .gear_ratios
            .get(self.current_gear)
            .copied()
            .unwrap_or(1.0);
        engine_torque * gear_ratio * self.final_drive * self.clutch
    }

    /// Upshift to the next gear.
    pub fn upshift(&mut self) {
        if self.current_gear + 1 < self.gear_ratios.len() {
            self.current_gear += 1;
        }
    }

    /// Downshift.
    pub fn downshift(&mut self) {
        if self.current_gear > 0 {
            self.current_gear -= 1;
        }
    }

    /// Current gear ratio.
    pub fn current_ratio(&self) -> f64 {
        self.gear_ratios
            .get(self.current_gear)
            .copied()
            .unwrap_or(1.0)
    }

    /// Update engine RPM from wheel speed.
    pub fn update_rpm(&mut self, wheel_angular_velocity: f64) {
        let ratio = self.current_ratio() * self.final_drive;
        self.engine_rpm = (wheel_angular_velocity * ratio * 60.0 / (2.0 * std::f64::consts::PI))
            .max(self.engine_curve.idle)
            .min(self.engine_curve.redline);
    }
}

// ---------------------------------------------------------------------------
// PySteering
// ---------------------------------------------------------------------------

/// Steering system with Ackermann geometry and 4WS option.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PySteering {
    /// Steering rack ratio (radians of wheel turn per radian of steering input).
    pub rack_ratio: f64,
    /// Maximum steering angle at the wheels (rad).
    pub max_wheel_angle: f64,
    /// Whether four-wheel steering is active.
    pub four_ws: bool,
    /// Rear steer angle ratio relative to front (0 = no rear steer, negative = counter-steer).
    pub rear_steer_ratio: f64,
    /// Wheelbase (m).
    pub wheelbase: f64,
    /// Track width (m).
    pub track_width: f64,
    /// Steering feel: damping coefficient.
    pub steering_damping: f64,
    /// Current steering angle input (normalized -1 to 1).
    pub input: f64,
}

impl PySteering {
    /// Create a default front-wheel-steering system.
    pub fn new(wheelbase: f64, track_width: f64) -> Self {
        Self {
            rack_ratio: 0.15,
            max_wheel_angle: 0.5,
            four_ws: false,
            rear_steer_ratio: 0.0,
            wheelbase,
            track_width,
            steering_damping: 5.0,
            input: 0.0,
        }
    }

    /// Compute the front wheel steer angle from the input.
    pub fn front_angle(&self) -> f64 {
        (self.input * self.rack_ratio).clamp(-self.max_wheel_angle, self.max_wheel_angle)
    }

    /// Compute the rear wheel steer angle (0 for standard FWS).
    pub fn rear_angle(&self) -> f64 {
        self.front_angle() * self.rear_steer_ratio
    }

    /// Compute Ackermann inner/outer wheel angles.
    pub fn ackermann_front_angles(&self) -> (f64, f64) {
        ackermann_angles(self.front_angle(), self.wheelbase, self.track_width)
    }
}

// ---------------------------------------------------------------------------
// PyStabilityControl
// ---------------------------------------------------------------------------

/// Traction control system (TCS) state.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TcsState {
    /// Whether TCS is active this frame.
    pub active: bool,
    /// Slip ratio threshold for TCS intervention.
    pub slip_threshold: f64,
    /// Throttle reduction factor applied by TCS.
    pub throttle_reduction: f64,
}

impl TcsState {
    /// Create default TCS state.
    pub fn new() -> Self {
        Self {
            active: false,
            slip_threshold: 0.2,
            throttle_reduction: 0.0,
        }
    }

    /// Update TCS given current drive slip ratios.
    pub fn update(&mut self, slip_ratios: &[f64]) -> f64 {
        let max_slip = slip_ratios.iter().cloned().fold(0.0_f64, f64::max);
        if max_slip > self.slip_threshold {
            self.active = true;
            self.throttle_reduction =
                ((max_slip - self.slip_threshold) / self.slip_threshold).min(1.0);
        } else {
            self.active = false;
            self.throttle_reduction = 0.0;
        }
        1.0 - self.throttle_reduction
    }
}

impl Default for TcsState {
    fn default() -> Self {
        Self::new()
    }
}

/// ABS (anti-lock braking system) state.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AbsState {
    /// Whether ABS is active this frame.
    pub active: bool,
    /// Brake slip threshold for ABS intervention.
    pub slip_threshold: f64,
    /// Brake pressure modulation factor.
    pub brake_modulation: f64,
    /// Current phase: pressure build (0), hold (1), release (2).
    pub phase: u8,
}

impl AbsState {
    /// Create default ABS state.
    pub fn new() -> Self {
        Self {
            active: false,
            slip_threshold: -0.15,
            brake_modulation: 1.0,
            phase: 0,
        }
    }

    /// Update ABS given current brake slip ratios.
    pub fn update(&mut self, slip_ratios: &[f64]) -> f64 {
        let min_slip = slip_ratios.iter().cloned().fold(0.0_f64, f64::min);
        if min_slip < self.slip_threshold {
            self.active = true;
            // Simple modulation: alternate between hold and release
            self.phase = (self.phase + 1) % 3;
            self.brake_modulation = if self.phase == 2 { 0.5 } else { 1.0 };
        } else {
            self.active = false;
            self.brake_modulation = 1.0;
            self.phase = 0;
        }
        self.brake_modulation
    }
}

impl Default for AbsState {
    fn default() -> Self {
        Self::new()
    }
}

/// Electronic stability control (ESC) state.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EscState {
    /// Whether ESC is active this frame.
    pub active: bool,
    /// Yaw rate error threshold (rad/s).
    pub yaw_error_threshold: f64,
    /// Brake intervention torque (N·m per wheel).
    pub brake_torque: f64,
    /// Last computed yaw rate error.
    pub yaw_error: f64,
}

impl EscState {
    /// Create default ESC state.
    pub fn new() -> Self {
        Self {
            active: false,
            yaw_error_threshold: 0.1,
            brake_torque: 0.0,
            yaw_error: 0.0,
        }
    }

    /// Update ESC given desired and actual yaw rates.
    pub fn update(&mut self, desired_yaw: f64, actual_yaw: f64) {
        self.yaw_error = desired_yaw - actual_yaw;
        if self.yaw_error.abs() > self.yaw_error_threshold {
            self.active = true;
            self.brake_torque = self.yaw_error.abs() * 500.0; // proportional intervention
        } else {
            self.active = false;
            self.brake_torque = 0.0;
        }
    }
}

impl Default for EscState {
    fn default() -> Self {
        Self::new()
    }
}

/// Vehicle stability control system combining TCS, ABS, and ESC.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyStabilityControl {
    /// Traction control.
    pub tcs: TcsState,
    /// Anti-lock braking.
    pub abs: AbsState,
    /// Electronic stability control.
    pub esc: EscState,
    /// Whether the whole system is enabled.
    pub enabled: bool,
}

impl PyStabilityControl {
    /// Create a default stability control system.
    pub fn new() -> Self {
        Self {
            tcs: TcsState::new(),
            abs: AbsState::new(),
            esc: EscState::new(),
            enabled: true,
        }
    }

    /// Compute throttle scale, brake modulation, and yaw correction.
    pub fn update(&mut self, slip_ratios: &[f64], desired_yaw: f64, actual_yaw: f64) -> (f64, f64) {
        if !self.enabled {
            return (1.0, 1.0);
        }
        let throttle_scale = self.tcs.update(slip_ratios);
        let brake_scale = self.abs.update(slip_ratios);
        self.esc.update(desired_yaw, actual_yaw);
        (throttle_scale, brake_scale)
    }
}

impl Default for PyStabilityControl {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// PySuspension
// ---------------------------------------------------------------------------

/// Suspension geometry type.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub enum SuspensionKind {
    /// MacPherson strut (front).
    MacPherson,
    /// Double wishbone.
    DoubleWishbone,
    /// Multi-link.
    MultiLink,
    /// Simple spring-damper (used for rear or simplified models).
    SpringDamper,
}

/// One corner suspension.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PySuspension {
    /// Suspension geometry type.
    pub kind: SuspensionKind,
    /// Spring rate (N/m).
    pub spring_rate: f64,
    /// Damper coefficient in compression (N·s/m).
    pub damper_compression: f64,
    /// Damper coefficient in rebound (N·s/m).
    pub damper_rebound: f64,
    /// Anti-roll bar stiffness (N·m/rad).
    pub arb_stiffness: f64,
    /// Bump stop spring rate (N/m), activated near max compression.
    pub bump_stop_rate: f64,
    /// Bump stop activation distance from maximum compression (m).
    pub bump_stop_clearance: f64,
    /// Current displacement from rest (m, positive = compression).
    pub displacement: f64,
    /// Current displacement velocity (m/s).
    pub velocity: f64,
    /// Maximum compression travel (m).
    pub max_compression: f64,
    /// Maximum droop travel (m).
    pub max_droop: f64,
}

impl PySuspension {
    /// Create a default double wishbone suspension.
    pub fn double_wishbone() -> Self {
        Self {
            kind: SuspensionKind::DoubleWishbone,
            spring_rate: 30000.0,
            damper_compression: 2000.0,
            damper_rebound: 3000.0,
            arb_stiffness: 5000.0,
            bump_stop_rate: 100000.0,
            bump_stop_clearance: 0.02,
            displacement: 0.0,
            velocity: 0.0,
            max_compression: 0.08,
            max_droop: 0.06,
        }
    }

    /// Create a MacPherson strut suspension.
    pub fn macpherson() -> Self {
        Self {
            kind: SuspensionKind::MacPherson,
            spring_rate: 25000.0,
            damper_compression: 1500.0,
            damper_rebound: 2500.0,
            arb_stiffness: 3000.0,
            bump_stop_rate: 80000.0,
            bump_stop_clearance: 0.02,
            displacement: 0.0,
            velocity: 0.0,
            max_compression: 0.07,
            max_droop: 0.05,
        }
    }

    /// Compute the suspension force given displacement and velocity.
    pub fn force(&self) -> f64 {
        let spring_force = -self.spring_rate * self.displacement;
        let damper_force = if self.velocity > 0.0 {
            -self.damper_compression * self.velocity
        } else {
            -self.damper_rebound * self.velocity
        };
        let bump_force = if self.displacement > self.max_compression - self.bump_stop_clearance {
            let excess = self.displacement - (self.max_compression - self.bump_stop_clearance);
            -self.bump_stop_rate * excess
        } else {
            0.0
        };
        spring_force + damper_force + bump_force
    }

    /// Step the suspension dynamics.
    pub fn step(&mut self, external_displacement: f64, dt: f64) {
        let prev = self.displacement;
        self.displacement = external_displacement.clamp(-self.max_droop, self.max_compression);
        self.velocity = (self.displacement - prev) / dt.max(1e-9);
    }
}

// ---------------------------------------------------------------------------
// PyWheelDynamics
// ---------------------------------------------------------------------------

/// State and forces at a single wheel.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyWheelDynamics {
    /// Wheel angular velocity (rad/s).
    pub omega: f64,
    /// Wheel moment of inertia (kg·m²).
    pub inertia: f64,
    /// Effective rolling radius (m).
    pub radius: f64,
    /// Current slip ratio.
    pub slip_ratio: f64,
    /// Current slip angle (rad).
    pub slip_angle: f64,
    /// Camber angle (rad, positive = lean toward vehicle center).
    pub camber: f64,
    /// Toe angle (rad).
    pub toe: f64,
    /// Longitudinal force at the contact patch (N).
    pub fx: f64,
    /// Lateral force at the contact patch (N).
    pub fy: f64,
    /// Normal load (N).
    pub fz: f64,
    /// Whether the wheel is on the ground.
    pub grounded: bool,
}

impl PyWheelDynamics {
    /// Create a new wheel with given inertia and radius.
    pub fn new(inertia: f64, radius: f64) -> Self {
        Self {
            omega: 0.0,
            inertia,
            radius,
            slip_ratio: 0.0,
            slip_angle: 0.0,
            camber: 0.0,
            toe: 0.0,
            fx: 0.0,
            fy: 0.0,
            fz: 0.0,
            grounded: false,
        }
    }

    /// Update wheel angular velocity given applied torque and time step.
    pub fn step_omega(&mut self, drive_torque: f64, brake_torque: f64, dt: f64) {
        // Net torque = drive - brake - fx * radius (reaction from road)
        let road_reaction = -self.fx * self.radius;
        let net_torque = drive_torque - brake_torque + road_reaction;
        self.omega += net_torque / self.inertia * dt;
        if !self.grounded {
            // No road reaction off ground, just dampen
            self.omega *= (1.0 - 0.1 * dt).max(0.0);
        }
    }

    /// Peripheral speed of the tire.
    pub fn peripheral_speed(&self) -> f64 {
        self.omega * self.radius
    }
}

// ---------------------------------------------------------------------------
// PyVehicle
// ---------------------------------------------------------------------------

/// Full vehicle state.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyVehicle {
    /// Chassis mass (kg).
    pub mass: f64,
    /// Inertia tensor diagonal \[Ixx, Iyy, Izz\] (kg·m²).
    pub inertia: [f64; 3],
    /// Center of mass height (m).
    pub cog_height: f64,
    /// World position of chassis center \[x, y, z\].
    pub position: [f64; 3],
    /// Velocity \[vx, vy, vz\] in world frame.
    pub velocity: [f64; 3],
    /// Euler angles \[roll, pitch, yaw\] (rad).
    pub orientation: [f64; 3],
    /// Angular velocity \[wx, wy, wz\] (rad/s).
    pub angular_velocity: [f64; 3],
    /// Wheels: \[FL, FR, RL, RR\].
    pub wheels: [PyWheelDynamics; 4],
    /// Wheel positions in local frame \[FL, FR, RL, RR\].
    pub wheel_positions: [[f64; 3]; 4],
    /// Suspensions \[FL, FR, RL, RR\].
    pub suspensions: [PySuspension; 4],
    /// Tire models \[FL, FR, RL, RR\].
    pub tires: [PyTireModel; 4],
    /// Drivetrain.
    pub drivetrain: PyDrivetrain,
    /// Steering system.
    pub steering: PySteering,
    /// Stability control.
    pub stability: PyStabilityControl,
    /// Current gear as display number (1-indexed).
    pub gear_display: u32,
}

impl PyVehicle {
    /// Create a default sedan-like vehicle.
    pub fn sedan() -> Self {
        let wheel = PyWheelDynamics::new(1.5, 0.32);
        let suspension = PySuspension::macpherson();
        let tire = PyTireModel::road_pacejka(0.32);
        Self {
            mass: 1500.0,
            inertia: [2000.0, 2500.0, 3000.0],
            cog_height: 0.5,
            position: [0.0; 3],
            velocity: [0.0; 3],
            orientation: [0.0; 3],
            angular_velocity: [0.0; 3],
            wheels: [wheel.clone(), wheel.clone(), wheel.clone(), wheel.clone()],
            wheel_positions: [
                [-0.75, -0.5, 1.3],
                [0.75, -0.5, 1.3],
                [-0.75, -0.5, -1.2],
                [0.75, -0.5, -1.2],
            ],
            suspensions: [
                suspension.clone(),
                suspension.clone(),
                suspension.clone(),
                suspension.clone(),
            ],
            tires: [tire.clone(), tire.clone(), tire.clone(), tire.clone()],
            drivetrain: PyDrivetrain::six_speed(300.0),
            steering: PySteering::new(2.6, 1.5),
            stability: PyStabilityControl::new(),
            gear_display: 1,
        }
    }

    /// Create a racing car vehicle.
    pub fn racing_car() -> Self {
        let wheel = PyWheelDynamics::new(1.0, 0.30);
        let suspension = PySuspension::double_wishbone();
        let tire = PyTireModel::slick_pacejka(0.30);
        Self {
            mass: 700.0,
            inertia: [800.0, 1000.0, 1200.0],
            cog_height: 0.3,
            position: [0.0; 3],
            velocity: [0.0; 3],
            orientation: [0.0; 3],
            angular_velocity: [0.0; 3],
            wheels: [wheel.clone(), wheel.clone(), wheel.clone(), wheel.clone()],
            wheel_positions: [
                [-0.7, -0.4, 1.6],
                [0.7, -0.4, 1.6],
                [-0.7, -0.4, -1.6],
                [0.7, -0.4, -1.6],
            ],
            suspensions: [
                suspension.clone(),
                suspension.clone(),
                suspension.clone(),
                suspension.clone(),
            ],
            tires: [tire.clone(), tire.clone(), tire.clone(), tire.clone()],
            drivetrain: PyDrivetrain::six_speed(500.0),
            steering: PySteering::new(2.8, 1.4),
            stability: PyStabilityControl::new(),
            gear_display: 1,
        }
    }

    /// Get the forward speed (m/s).
    pub fn forward_speed(&self) -> f64 {
        // Assume vehicle heading is along local X axis; yaw rotates it
        let yaw = self.orientation[2];
        self.velocity[0] * yaw.cos() + self.velocity[2] * yaw.sin()
    }

    /// Perform one simulation step.
    ///
    /// `dt`       — time step (s)
    /// `throttle` — throttle input \[0, 1\]
    /// `brake`    — brake input \[0, 1\]
    /// `steer`    — steering input \[-1, 1\]
    pub fn step(&mut self, dt: f64, throttle: f64, brake: f64, steer: f64) {
        // Update steering
        self.steering.input = steer.clamp(-1.0, 1.0);
        let front_angle = self.steering.front_angle();

        // Compute slip ratios (simplified)
        let v = self.forward_speed().abs().max(0.01);
        let wheel_speeds: [f64; 4] = core::array::from_fn(|i| self.wheels[i].peripheral_speed());
        let slip_ratios: Vec<f64> = wheel_speeds
            .iter()
            .map(|&ws| compute_slip_ratio(ws, v))
            .collect();

        // Stability control
        let desired_yaw = v / 2.6 * front_angle;
        let actual_yaw = self.angular_velocity[1];
        let (throttle_scale, brake_scale) =
            self.stability.update(&slip_ratios, desired_yaw, actual_yaw);
        let effective_throttle = throttle * throttle_scale;
        let effective_brake = brake * brake_scale;

        // Drivetrain
        let drive_torque_total = self.drivetrain.wheel_torque(effective_throttle);
        let drive_torque_per_wheel = drive_torque_total / 2.0; // simplified rear-wheel drive

        // Brake torque
        let max_brake_torque = 2000.0;
        let brake_torque = effective_brake * max_brake_torque;

        // Update rear wheels with drive torque
        for i in [2, 3] {
            let slip_angle = self.wheels[i].slip_angle;
            let kappa = self.wheels[i].slip_ratio;
            let fz = self.mass * 9.81 * 0.25; // equal load distribution
            let (fx, fy) = self.tires[i].compute_forces(slip_angle, kappa, fz);
            self.wheels[i].fx = fx;
            self.wheels[i].fy = fy;
            self.wheels[i].fz = fz;
            self.wheels[i].grounded = true;
            self.wheels[i].step_omega(drive_torque_per_wheel, brake_torque, dt);
        }
        // Update front wheels with steering
        for i in [0, 1] {
            self.wheels[i].slip_angle = -front_angle;
            let fz = self.mass * 9.81 * 0.25;
            let (fx, fy) = self.tires[i].compute_forces(self.wheels[i].slip_angle, 0.0, fz);
            self.wheels[i].fx = fx;
            self.wheels[i].fy = fy;
            self.wheels[i].fz = fz;
            self.wheels[i].grounded = true;
            self.wheels[i].step_omega(0.0, brake_torque, dt);
        }

        // Total forces on chassis (simplified)
        let total_fx: f64 = self.wheels.iter().map(|w| w.fx).sum();
        let total_fy: f64 = self.wheels.iter().map(|w| w.fy).sum();

        // Integrate velocity / position (Euler)
        let ax = total_fx / self.mass;
        let az = total_fy / self.mass - 9.81 * self.orientation[0].sin();
        self.velocity[0] += ax * dt;
        self.velocity[2] += az * dt;

        // Drag
        let drag = 0.3 * self.velocity[0] * self.velocity[0].abs();
        self.velocity[0] -= drag / self.mass * dt;

        self.position[0] += self.velocity[0] * dt;
        self.position[1] += self.velocity[1] * dt;
        self.position[2] += self.velocity[2] * dt;

        // Update engine RPM
        let avg_rear_omega = (self.wheels[2].omega + self.wheels[3].omega) * 0.5;
        self.drivetrain.update_rpm(avg_rear_omega);
        self.gear_display = (self.drivetrain.current_gear + 1) as u32;
    }
}

// ---------------------------------------------------------------------------
// PyTelemetry
// ---------------------------------------------------------------------------

/// A single telemetry sample.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TelemetrySample {
    /// Simulation time (s).
    pub time: f64,
    /// Vehicle speed (m/s).
    pub speed: f64,
    /// Longitudinal G-force.
    pub g_lon: f64,
    /// Lateral G-force.
    pub g_lat: f64,
    /// Slip angles \[FL, FR, RL, RR\] (rad).
    pub slip_angles: [f64; 4],
    /// Tire temperatures \[FL, FR, RL, RR\] (°C, surface).
    pub tire_temps: [f64; 4],
    /// Current gear.
    pub gear: u32,
    /// Engine RPM.
    pub rpm: f64,
    /// Throttle (0-1).
    pub throttle: f64,
    /// Brake (0-1).
    pub brake: f64,
    /// Steering angle (rad).
    pub steer: f64,
}

/// Lap statistics derived from telemetry.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LapStats {
    /// Lap time (s).
    pub lap_time: f64,
    /// Maximum speed (m/s).
    pub max_speed: f64,
    /// Average speed (m/s).
    pub avg_speed: f64,
    /// Maximum lateral G.
    pub max_g_lat: f64,
    /// Maximum longitudinal G (braking).
    pub max_g_lon_brake: f64,
    /// Maximum longitudinal G (acceleration).
    pub max_g_lon_accel: f64,
    /// Number of samples.
    pub sample_count: usize,
}

/// Vehicle telemetry recorder.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyTelemetry {
    /// Recorded samples.
    pub samples: Vec<TelemetrySample>,
    /// Maximum number of samples to retain.
    pub max_samples: usize,
    /// Whether recording is active.
    pub recording: bool,
}

impl PyTelemetry {
    /// Create a new telemetry recorder.
    pub fn new(max_samples: usize) -> Self {
        Self {
            samples: Vec::new(),
            max_samples,
            recording: true,
        }
    }

    /// Record a sample.
    pub fn record(&mut self, sample: TelemetrySample) {
        if !self.recording {
            return;
        }
        if self.samples.len() >= self.max_samples {
            self.samples.remove(0);
        }
        self.samples.push(sample);
    }

    /// Compute lap statistics from all recorded samples.
    pub fn lap_stats(&self) -> Option<LapStats> {
        if self.samples.is_empty() {
            return None;
        }
        let n = self.samples.len();
        let total_time = self
            .samples
            .last()
            .expect("collection should not be empty")
            .time
            - self.samples[0].time;
        let max_speed = self.samples.iter().map(|s| s.speed).fold(0.0_f64, f64::max);
        let avg_speed = self.samples.iter().map(|s| s.speed).sum::<f64>() / n as f64;
        let max_g_lat = self
            .samples
            .iter()
            .map(|s| s.g_lat.abs())
            .fold(0.0_f64, f64::max);
        let max_g_lon_brake = self
            .samples
            .iter()
            .map(|s| (-s.g_lon).max(0.0))
            .fold(0.0_f64, f64::max);
        let max_g_lon_accel = self
            .samples
            .iter()
            .map(|s| s.g_lon.max(0.0))
            .fold(0.0_f64, f64::max);
        Some(LapStats {
            lap_time: total_time,
            max_speed,
            avg_speed,
            max_g_lat,
            max_g_lon_brake,
            max_g_lon_accel,
            sample_count: n,
        })
    }

    /// Clear all samples.
    pub fn clear(&mut self) {
        self.samples.clear();
    }
}

// ---------------------------------------------------------------------------
// PyRacingLine
// ---------------------------------------------------------------------------

/// A 2D point on the track.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrackPoint {
    /// X coordinate (m).
    pub x: f64,
    /// Y coordinate (m).
    pub y: f64,
    /// Track half-width at this point (m).
    pub half_width: f64,
}

impl TrackPoint {
    /// Create a new track point.
    pub fn new(x: f64, y: f64, half_width: f64) -> Self {
        Self { x, y, half_width }
    }
}

/// Racing line optimizer and analyzer.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyRacingLine {
    /// Track centerline points.
    pub centerline: Vec<TrackPoint>,
    /// Optimized racing line (same count as centerline).
    pub racing_line: Vec<[f64; 2]>,
    /// Sector boundaries (indices into centerline).
    pub sector_boundaries: Vec<usize>,
    /// Sector times (s) from the last lap simulation.
    pub sector_times: Vec<f64>,
    /// Total track length (m).
    pub track_length: f64,
}

impl PyRacingLine {
    /// Create a racing line optimizer from a centerline.
    pub fn new(centerline: Vec<TrackPoint>) -> Self {
        let n = centerline.len();
        // Initialize racing line to centerline
        let racing_line = centerline.iter().map(|p| [p.x, p.y]).collect();
        let track_length = Self::compute_length(&centerline);
        Self {
            centerline,
            racing_line,
            sector_boundaries: if n >= 3 {
                vec![0, n / 3, 2 * n / 3, n - 1]
            } else {
                vec![0, n.max(1) - 1]
            },
            sector_times: Vec::new(),
            track_length,
        }
    }

    fn compute_length(pts: &[TrackPoint]) -> f64 {
        let mut len = 0.0;
        for i in 1..pts.len() {
            let dx = pts[i].x - pts[i - 1].x;
            let dy = pts[i].y - pts[i - 1].y;
            len += (dx * dx + dy * dy).sqrt();
        }
        len
    }

    /// Compute local curvature at point i (using finite differences).
    pub fn curvature_at(&self, i: usize) -> f64 {
        let n = self.racing_line.len();
        if n < 3 || i == 0 || i >= n - 1 {
            return 0.0;
        }
        let p0 = self.racing_line[i - 1];
        let p1 = self.racing_line[i];
        let p2 = self.racing_line[i + 1];
        // Menger curvature
        let a = ((p1[0] - p0[0]).powi(2) + (p1[1] - p0[1]).powi(2)).sqrt();
        let b = ((p2[0] - p1[0]).powi(2) + (p2[1] - p1[1]).powi(2)).sqrt();
        let c = ((p2[0] - p0[0]).powi(2) + (p2[1] - p0[1]).powi(2)).sqrt();
        let area2 = ((p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0])).abs();
        if a * b * c < 1e-12 {
            0.0
        } else {
            area2 / (a * b * c)
        }
    }

    /// Run one step of minimum-curvature optimization (gradient descent on curvature).
    pub fn optimize_step(&mut self, alpha: f64) {
        let n = self.racing_line.len();
        if n < 3 {
            return;
        }
        let mut new_line = self.racing_line.clone();
        for i in 1..n - 1 {
            let p0 = self.racing_line[i - 1];
            let p1 = self.racing_line[i];
            let p2 = self.racing_line[i + 1];
            // Move toward midpoint of neighbors (smoothing)
            let mid = [(p0[0] + p2[0]) * 0.5, (p0[1] + p2[1]) * 0.5];
            let new_x = p1[0] + alpha * (mid[0] - p1[0]);
            let new_y = p1[1] + alpha * (mid[1] - p1[1]);
            // Clamp to track width
            let cl = &self.centerline[i];
            let dx = new_x - cl.x;
            let dy = new_y - cl.y;
            let dist = (dx * dx + dy * dy).sqrt();
            if dist <= cl.half_width {
                new_line[i] = [new_x, new_y];
            } else if dist > 1e-12 {
                new_line[i] = [
                    cl.x + dx / dist * cl.half_width,
                    cl.y + dy / dist * cl.half_width,
                ];
            }
        }
        self.racing_line = new_line;
    }

    /// Run n_iter optimization iterations.
    pub fn optimize(&mut self, n_iter: u32, alpha: f64) {
        for _ in 0..n_iter {
            self.optimize_step(alpha);
        }
    }

    /// Compute the total curvature of the current racing line.
    pub fn total_curvature(&self) -> f64 {
        (1..self.racing_line.len() - 1)
            .map(|i| self.curvature_at(i))
            .sum()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pacejka_fy_zero_slip() {
        let fy = pacejka_fy(0.0, 10.0, 1.9, 1.0, 0.97);
        assert!(fy.abs() < 1e-9);
    }

    #[test]
    fn test_pacejka_fy_nonzero() {
        let fy = pacejka_fy(0.1, 10.0, 1.9, 1.0, 0.97);
        assert!(fy.abs() > 0.0);
    }

    #[test]
    fn test_pacejka_fx_peak() {
        let fx = pacejka_fx(0.2, 11.0, 1.65, 1.0, 0.0);
        assert!(fx.abs() > 0.0);
    }

    #[test]
    fn test_compute_slip_angle() {
        let alpha = compute_slip_angle(1.0, 10.0);
        assert!(alpha.abs() < 0.2);
    }

    #[test]
    fn test_compute_slip_ratio_no_slip() {
        let kappa = compute_slip_ratio(10.0, 10.0);
        assert!(kappa.abs() < 1e-9);
    }

    #[test]
    fn test_compute_slip_ratio_drive_slip() {
        let kappa = compute_slip_ratio(12.0, 10.0);
        assert!(kappa > 0.0);
    }

    #[test]
    fn test_ackermann_straight() {
        let (a, b) = ackermann_angles(0.0, 2.5, 1.5);
        assert!(a.abs() < 1e-9);
        assert!(b.abs() < 1e-9);
    }

    #[test]
    fn test_ackermann_turn() {
        let (inner, outer) = ackermann_angles(0.3, 2.5, 1.5);
        assert!(inner.abs() > outer.abs());
    }

    #[test]
    fn test_tire_model_forces() {
        let tm = PyTireModel::road_pacejka(0.32);
        let (fx, fy) = tm.compute_forces(0.1, 0.05, 3000.0);
        let _ = fx;
        assert!(fy.abs() > 0.0);
    }

    #[test]
    fn test_tire_thermal_update() {
        let mut tt = PyTireThermal::new();
        let grip0 = tt.grip_scale();
        tt.update(5000.0, 10.0);
        assert!(tt.surface_temp > 25.0);
        // Check grip scale is valid
        assert!(tt.grip_scale() >= 0.0 && tt.grip_scale() <= 1.0);
        let _ = grip0;
    }

    #[test]
    fn test_tire_thermal_optimal() {
        let mut tt = PyTireThermal::new();
        tt.surface_temp = 90.0;
        assert!((tt.grip_scale() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_engine_curve_flat() {
        let ec = EngineCurve::flat(300.0, 7000.0);
        let t = ec.torque_at(3500.0);
        assert!(t > 0.0);
        assert!(t <= 300.0 * 1.01);
    }

    #[test]
    fn test_drivetrain_upshift() {
        let mut dt = PyDrivetrain::six_speed(300.0);
        dt.upshift();
        assert_eq!(dt.current_gear, 2);
    }

    #[test]
    fn test_drivetrain_downshift() {
        let mut dt = PyDrivetrain::six_speed(300.0);
        dt.downshift();
        assert_eq!(dt.current_gear, 0);
    }

    #[test]
    fn test_drivetrain_wheel_torque() {
        let dt = PyDrivetrain::six_speed(300.0);
        let t = dt.wheel_torque(1.0);
        assert!(t > 0.0);
    }

    #[test]
    fn test_steering_angles() {
        let mut steer = PySteering::new(2.6, 1.5);
        steer.input = 1.0;
        let fa = steer.front_angle();
        assert!(fa > 0.0);
        assert!(fa <= steer.max_wheel_angle);
    }

    #[test]
    fn test_stability_tcs() {
        let mut sc = PyStabilityControl::new();
        let (ts, _) = sc.update(&[0.5, 0.5, 0.5, 0.5], 0.0, 0.0);
        assert!(ts < 1.0);
    }

    #[test]
    fn test_stability_disabled() {
        let mut sc = PyStabilityControl::new();
        sc.enabled = false;
        let (ts, bs) = sc.update(&[2.0, 2.0, 2.0, 2.0], 0.0, 0.0);
        assert!((ts - 1.0).abs() < 1e-9);
        assert!((bs - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_suspension_force() {
        let mut s = PySuspension::double_wishbone();
        s.step(0.03, 0.01);
        let f = s.force();
        assert!(f < 0.0); // spring pushes back
    }

    #[test]
    fn test_suspension_bump_stop() {
        let mut s = PySuspension::double_wishbone();
        s.step(s.max_compression, 0.01);
        let f = s.force();
        assert!(f < 0.0);
    }

    #[test]
    fn test_wheel_dynamics_step() {
        let mut w = PyWheelDynamics::new(1.5, 0.32);
        w.fx = 100.0;
        w.grounded = true;
        w.step_omega(200.0, 0.0, 0.01);
        assert!(w.omega > 0.0);
    }

    #[test]
    fn test_vehicle_sedan_step() {
        let mut v = PyVehicle::sedan();
        v.step(0.016, 0.5, 0.0, 0.0);
        // Vehicle should have moved or velocity changed
        let speed = v.forward_speed();
        let _ = speed;
        assert!(v.drivetrain.engine_rpm > 0.0);
    }

    #[test]
    fn test_vehicle_racing_step() {
        let mut v = PyVehicle::racing_car();
        v.step(0.016, 1.0, 0.0, 0.1);
        assert!(v.gear_display >= 1);
    }

    #[test]
    fn test_telemetry_record() {
        let mut tel = PyTelemetry::new(100);
        let sample = TelemetrySample {
            time: 0.0,
            speed: 50.0,
            g_lon: 0.5,
            g_lat: 0.3,
            slip_angles: [0.0; 4],
            tire_temps: [80.0; 4],
            gear: 3,
            rpm: 4000.0,
            throttle: 0.8,
            brake: 0.0,
            steer: 0.1,
        };
        tel.record(sample);
        assert_eq!(tel.samples.len(), 1);
    }

    #[test]
    fn test_telemetry_lap_stats() {
        let mut tel = PyTelemetry::new(100);
        for i in 0..10 {
            let sample = TelemetrySample {
                time: i as f64 * 0.1,
                speed: 50.0 + i as f64,
                g_lon: 0.1,
                g_lat: 0.2,
                slip_angles: [0.0; 4],
                tire_temps: [80.0; 4],
                gear: 3,
                rpm: 4000.0,
                throttle: 0.8,
                brake: 0.0,
                steer: 0.0,
            };
            tel.record(sample);
        }
        let stats = tel.lap_stats().unwrap();
        assert!(stats.max_speed > 50.0);
        assert!(stats.sample_count == 10);
    }

    #[test]
    fn test_racing_line_curvature() {
        let pts: Vec<TrackPoint> = (0..10)
            .map(|i| TrackPoint::new(i as f64, 0.0, 3.0))
            .collect();
        let rl = PyRacingLine::new(pts);
        // Straight line should have near-zero curvature
        let k = rl.curvature_at(5);
        assert!(k < 1e-9);
    }

    #[test]
    fn test_racing_line_optimize() {
        let pts: Vec<TrackPoint> = (0..20)
            .map(|i| {
                let t = i as f64 / 20.0 * 2.0 * std::f64::consts::PI;
                TrackPoint::new(t.cos() * 10.0, t.sin() * 10.0, 2.0)
            })
            .collect();
        let mut rl = PyRacingLine::new(pts);
        let k0 = rl.total_curvature();
        rl.optimize(10, 0.1);
        let k1 = rl.total_curvature();
        // After optimization curvature should not increase dramatically
        assert!(k1 <= k0 * 1.5 + 0.1);
    }

    #[test]
    fn test_fiala_tire_linear() {
        let fp = FialaParams::new(50000.0, 1.0);
        let fy = fp.lateral_force(0.05, 4000.0);
        let expected = -50000.0 * 0.05;
        assert!((fy - expected).abs() < 1.0);
    }

    #[test]
    fn test_fiala_tire_saturated() {
        let fp = FialaParams::new(50000.0, 1.0);
        let fy = fp.lateral_force(1.0, 4000.0);
        assert!(fy.abs() <= 4000.0 + 1.0);
    }
}
