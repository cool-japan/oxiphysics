// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Motorcycle dynamics simulation.
//!
//! This module implements a comprehensive single-track motorcycle model
//! including:
//!
//! - **Lean dynamics**: tilt angle, gyroscopic precession effects
//! - **Counter-steering model**: steering torque to lean coupling
//! - **Single-track (Pacejka) tire forces**: Magic Formula for both tyres
//! - **Tank-slapper / weave / wobble modes**: linearised eigenvalue analysis
//! - **Rider balance control**: PD lean-angle controller
//! - **Telescopic fork / monoshock suspension**: spring-damper models
//! - **Wheelie / stoppie detection**: threshold-based detection
//! - **Chain drive dynamics**: torque transmission with chain compliance
//! - **Aerodynamic fairing model**: drag & downforce for fairing geometries
//! - **Lean-angle limiter**: hard and soft lean limits with warning thresholds

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Gravitational acceleration (m/s²).
const G: f64 = 9.81;

/// Air density at sea level (kg/m³).
const RHO_AIR: f64 = 1.225;

// ---------------------------------------------------------------------------
// Utility helpers
// ---------------------------------------------------------------------------

/// Clamp `x` to `[lo, hi]`.
#[inline]
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}

/// Sign of `x` (returns +1.0, -1.0, or 0.0).
#[inline]
fn sign(x: f64) -> f64 {
    if x > 0.0 {
        1.0
    } else if x < 0.0 {
        -1.0
    } else {
        0.0
    }
}

/// Degrees to radians.
#[inline]
fn deg2rad(deg: f64) -> f64 {
    deg * PI / 180.0
}

/// Radians to degrees.
#[inline]
fn rad2deg(rad: f64) -> f64 {
    rad * 180.0 / PI
}

// ---------------------------------------------------------------------------
// PacejkaMotoTire
// ---------------------------------------------------------------------------

/// Pacejka Magic Formula coefficients for a single axis.
#[derive(Debug, Clone)]
pub struct PacejkaMotoCoeffs {
    /// Stiffness factor B.
    pub b: f64,
    /// Shape factor C.
    pub c: f64,
    /// Peak value factor D.
    pub d: f64,
    /// Curvature factor E.
    pub e: f64,
}

impl Default for PacejkaMotoCoeffs {
    fn default() -> Self {
        Self {
            b: 9.0,
            c: 1.65,
            d: 1.0,
            e: 0.95,
        }
    }
}

impl PacejkaMotoCoeffs {
    /// Evaluate `F = D * sin(C * atan(B*x - E*(B*x - atan(B*x))))`.
    pub fn evaluate(&self, x: f64) -> f64 {
        let bx = self.b * x;
        let phi = bx - self.e * (bx - bx.atan());
        self.d * (self.c * phi.atan()).sin()
    }

    /// Create a set of default motorcycle rear tyre lateral coefficients.
    pub fn default_lateral() -> Self {
        Self {
            b: 8.5,
            c: 1.55,
            d: 1.1,
            e: 0.9,
        }
    }

    /// Create default motorcycle tyre longitudinal coefficients.
    pub fn default_longitudinal() -> Self {
        Self {
            b: 11.0,
            c: 1.75,
            d: 1.05,
            e: 0.97,
        }
    }
}

/// Single-track Pacejka motorcycle tyre.
///
/// Computes lateral and longitudinal forces plus self-aligning torque,
/// scaling peak values by the current normal force.
#[derive(Debug, Clone)]
pub struct PacejkaMotoTire {
    /// Lateral Pacejka coefficients.
    pub lat: PacejkaMotoCoeffs,
    /// Longitudinal Pacejka coefficients.
    pub lon: PacejkaMotoCoeffs,
    /// Tyre relaxation length (m).
    pub relaxation_length: f64,
    /// Current filtered (relaxed) slip angle (rad).
    pub filtered_slip_angle: f64,
    /// Pneumatic trail constant (m) for self-aligning torque.
    pub pneumatic_trail: f64,
}

impl PacejkaMotoTire {
    /// Create a new tyre with default motorcycle coefficients.
    pub fn new() -> Self {
        Self {
            lat: PacejkaMotoCoeffs::default_lateral(),
            lon: PacejkaMotoCoeffs::default_longitudinal(),
            relaxation_length: 0.3,
            filtered_slip_angle: 0.0,
            pneumatic_trail: 0.03,
        }
    }

    /// Update the relaxation filter for the slip angle.
    pub fn update_relaxation(&mut self, raw_slip: f64, velocity: f64, dt: f64) {
        if velocity.abs() < 1e-3 {
            self.filtered_slip_angle = raw_slip;
            return;
        }
        // First-order relaxation: dα/ds = (α_raw - α) / σ
        let rate = velocity / self.relaxation_length;
        self.filtered_slip_angle += dt * rate * (raw_slip - self.filtered_slip_angle);
        self.filtered_slip_angle = clamp(self.filtered_slip_angle, -0.5, 0.5);
    }

    /// Compute lateral (Fy) and longitudinal (Fx) forces and self-aligning torque (Mz).
    ///
    /// # Arguments
    /// * `slip_angle` - lateral slip angle (rad)
    /// * `slip_ratio` - longitudinal slip ratio (dimensionless)
    /// * `normal_force` - vertical load on tyre (N)
    /// * `lean_angle` - camber/lean angle (rad), positive = right lean
    pub fn forces(
        &self,
        slip_angle: f64,
        slip_ratio: f64,
        normal_force: f64,
        lean_angle: f64,
    ) -> (f64, f64, f64) {
        // Camber thrust: lateral force from lean (linear model)
        let camber_stiffness = 0.2 * normal_force;
        let camber_thrust = camber_stiffness * lean_angle;

        let fy_pure = normal_force * self.lat.evaluate(slip_angle) + camber_thrust;
        let fx_pure = normal_force * self.lon.evaluate(slip_ratio);

        // Combined slip via Pythagorean friction ellipse
        let sigma = {
            let denom = (fy_pure * fy_pure + fx_pure * fx_pure).sqrt();
            if denom < 1.0 { 1.0 } else { denom }
        };
        let fy = fy_pure * (fy_pure.abs() / sigma).min(1.0);
        let fx = fx_pure * (fx_pure.abs() / sigma).min(1.0);

        // Self-aligning torque
        let mz = -self.pneumatic_trail * fy;

        (fy, fx, mz)
    }
}

impl Default for PacejkaMotoTire {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// TelescopicFork
// ---------------------------------------------------------------------------

/// Telescopic fork front suspension model.
///
/// Two co-axial spring-damper tubes providing both suspension travel and
/// anti-dive under braking.
#[derive(Debug, Clone)]
pub struct TelescopicFork {
    /// Spring stiffness (N/m).
    pub spring_rate: f64,
    /// Damping coefficient (N·s/m).
    pub damping: f64,
    /// Preload displacement (m).
    pub preload: f64,
    /// Maximum fork travel (m).
    pub max_travel: f64,
    /// Minimum fork travel — limits compression (m).
    pub min_travel: f64,
    /// Current compression displacement (m).
    pub compression: f64,
    /// Current velocity of compression (m/s).
    pub velocity: f64,
    /// Anti-dive coefficient — reduces dive under braking torque.
    pub anti_dive: f64,
}

impl TelescopicFork {
    /// Construct a new telescopic fork with typical sport-bike defaults.
    pub fn new() -> Self {
        Self {
            spring_rate: 9500.0,
            damping: 1200.0,
            preload: 0.01,
            max_travel: 0.12,
            min_travel: 0.0,
            compression: 0.02,
            velocity: 0.0,
            anti_dive: 0.3,
        }
    }

    /// Compute spring-damper force given compression displacement and velocity.
    ///
    /// Returns the upward (reaction) force in Newtons.
    pub fn force(&self, compression: f64, vel: f64) -> f64 {
        let spring_f = self.spring_rate * (compression + self.preload);
        let damper_f = self.damping * vel;
        (spring_f + damper_f).max(0.0)
    }

    /// Advance the fork by one time step.
    ///
    /// `input_vel` is the relative velocity of the axle with respect to the chassis.
    pub fn step(&mut self, input_vel: f64, dt: f64) {
        self.velocity = input_vel;
        self.compression += dt * input_vel;
        self.compression = clamp(self.compression, self.min_travel, self.max_travel);
    }

    /// Apply anti-dive correction based on braking deceleration.
    pub fn anti_dive_correction(&self, brake_torque: f64, wheel_radius: f64) -> f64 {
        // Anti-dive force proportional to brake torque
        self.anti_dive * brake_torque / wheel_radius.max(0.01)
    }
}

impl Default for TelescopicFork {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Monoshock
// ---------------------------------------------------------------------------

/// Rear monoshock (single shock absorber) suspension model.
///
/// Incorporates a motion-ratio linkage so wheel travel differs from shock travel.
#[derive(Debug, Clone)]
pub struct Monoshock {
    /// Spring stiffness at the shock (N/m).
    pub spring_rate: f64,
    /// Damping coefficient (N·s/m).
    pub damping: f64,
    /// Preload displacement (m).
    pub preload: f64,
    /// Motion ratio (wheel travel / shock travel), typically 0.6 – 0.75.
    pub motion_ratio: f64,
    /// Maximum shock travel (m).
    pub max_travel: f64,
    /// Current shock compression (m).
    pub compression: f64,
}

impl Monoshock {
    /// Construct a monoshock with typical street-bike defaults.
    pub fn new() -> Self {
        Self {
            spring_rate: 14000.0,
            damping: 1600.0,
            preload: 0.015,
            motion_ratio: 0.68,
            max_travel: 0.10,
            compression: 0.02,
        }
    }

    /// Compute the wheel-rate (effective spring at the axle).
    pub fn wheel_rate(&self) -> f64 {
        self.spring_rate * self.motion_ratio * self.motion_ratio
    }

    /// Compute the force at the wheel given wheel compression and velocity.
    pub fn wheel_force(&self, wheel_compression: f64, wheel_vel: f64) -> f64 {
        let shock_comp = wheel_compression * self.motion_ratio;
        let shock_vel = wheel_vel * self.motion_ratio;
        let spring_f = self.spring_rate * (shock_comp + self.preload);
        let damper_f = self.damping * shock_vel;
        ((spring_f + damper_f) * self.motion_ratio).max(0.0)
    }

    /// Advance the monoshock by one time step given wheel velocity.
    pub fn step(&mut self, wheel_vel: f64, dt: f64) {
        let shock_vel = wheel_vel * self.motion_ratio;
        self.compression += dt * shock_vel;
        self.compression = clamp(self.compression, 0.0, self.max_travel);
    }
}

impl Default for Monoshock {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// ChainDrive
// ---------------------------------------------------------------------------

/// Roller chain drive dynamics.
///
/// Models the torque multiplication through sprocket ratio and chain
/// compliance (modelled as a torsional spring-damper).
#[derive(Debug, Clone)]
pub struct ChainDrive {
    /// Number of teeth on the engine (front) sprocket.
    pub front_teeth: u32,
    /// Number of teeth on the rear (wheel) sprocket.
    pub rear_teeth: u32,
    /// Chain compliance stiffness (N·m/rad).
    pub stiffness: f64,
    /// Chain damping (N·m·s/rad).
    pub chain_damping: f64,
    /// Current angular wrap of chain (rad).
    pub wrap_angle: f64,
    /// Angular velocity of wrap (rad/s).
    pub wrap_vel: f64,
}

impl ChainDrive {
    /// Create a chain drive with typical sport-bike sprocket ratio and stiffness.
    pub fn new(front_teeth: u32, rear_teeth: u32) -> Self {
        Self {
            front_teeth,
            rear_teeth,
            stiffness: 4000.0,
            chain_damping: 50.0,
            wrap_angle: 0.0,
            wrap_vel: 0.0,
        }
    }

    /// Sprocket (overall) ratio: rear / front.
    pub fn sprocket_ratio(&self) -> f64 {
        self.rear_teeth as f64 / self.front_teeth as f64
    }

    /// Compute chain torque transmitted to the rear wheel given engine output
    /// torque and current gearbox ratio.
    pub fn transmitted_torque(&self, engine_torque: f64, gearbox_ratio: f64) -> f64 {
        engine_torque * gearbox_ratio * self.sprocket_ratio()
    }

    /// Compute chain tension force (N) given chain torque and rear sprocket radius.
    pub fn chain_tension(&self, chain_torque: f64, rear_sprocket_radius: f64) -> f64 {
        chain_torque / rear_sprocket_radius.max(0.001)
    }

    /// Advance chain compliance dynamics by one step.
    ///
    /// `engine_angle` is the cumulative engine shaft rotation (rad),
    /// `wheel_angle` is the rear wheel rotation (rad).
    pub fn step(&mut self, engine_angle: f64, wheel_angle: f64, dt: f64) {
        let target = engine_angle * self.sprocket_ratio();
        let error = target - wheel_angle - self.wrap_angle;
        let acc = (self.stiffness * error - self.chain_damping * self.wrap_vel)
            / (self.stiffness.max(1.0));
        self.wrap_vel += dt * acc;
        self.wrap_angle += dt * self.wrap_vel;
    }
}

// ---------------------------------------------------------------------------
// AerodynamicFairing
// ---------------------------------------------------------------------------

/// Aerodynamic fairing model for a motorcycle.
///
/// Computes drag force, lift (downforce), and pitching moment as functions of
/// speed and angle of attack.
#[derive(Debug, Clone)]
pub struct AerodynamicFairing {
    /// Frontal area (m²).
    pub frontal_area: f64,
    /// Drag coefficient at zero angle of attack.
    pub cd_zero: f64,
    /// Drag-area slope with angle of attack (1/rad).
    pub cd_alpha: f64,
    /// Lift coefficient slope (1/rad).
    pub cl_alpha: f64,
    /// Pitching moment coefficient (dimensionless, negative = nose-down).
    pub cm: f64,
    /// Wheelbase (m), used to compute pitching moment arm.
    pub wheelbase: f64,
}

impl AerodynamicFairing {
    /// Create a fairing model for a litre-class sport bike.
    pub fn new_sport_bike() -> Self {
        Self {
            frontal_area: 0.38,
            cd_zero: 0.55,
            cd_alpha: 0.10,
            cl_alpha: -0.30,
            cm: -0.05,
            wheelbase: 1.43,
        }
    }

    /// Compute drag force (N) at speed `v` (m/s) and angle of attack `aoa` (rad).
    pub fn drag_force(&self, v: f64, aoa: f64) -> f64 {
        let cd = self.cd_zero + self.cd_alpha * aoa.abs();
        0.5 * RHO_AIR * v * v * self.frontal_area * cd
    }

    /// Compute downforce (positive = down, N) at speed `v` and pitch angle `aoa`.
    pub fn downforce(&self, v: f64, aoa: f64) -> f64 {
        let cl = self.cl_alpha * aoa;
        -0.5 * RHO_AIR * v * v * self.frontal_area * cl
    }

    /// Compute pitching moment (N·m, positive = nose-up) at speed `v`.
    pub fn pitching_moment(&self, v: f64) -> f64 {
        0.5 * RHO_AIR * v * v * self.frontal_area * self.cm * self.wheelbase
    }
}

// ---------------------------------------------------------------------------
// LeanAngleLimiter
// ---------------------------------------------------------------------------

/// Lean-angle limiter: enforces a maximum safe lean angle and issues
/// graduated warnings as the bike approaches the limit.
#[derive(Debug, Clone)]
pub struct LeanAngleLimiter {
    /// Hard maximum lean angle (rad).  Beyond this, the lean is clamped.
    pub hard_limit: f64,
    /// Warning threshold (rad): below hard limit where warnings begin.
    pub warn_threshold: f64,
    /// Soft-limit zone width (rad): range over which a restoring torque is applied.
    pub soft_zone: f64,
    /// Soft-limit restoring torque gain (N·m/rad).
    pub soft_gain: f64,
}

impl LeanAngleLimiter {
    /// Create a limiter with typical sport-bike lean limits.
    ///
    /// Hard limit 58°, warning at 50°.
    pub fn new() -> Self {
        Self {
            hard_limit: deg2rad(58.0),
            warn_threshold: deg2rad(50.0),
            soft_zone: deg2rad(8.0),
            soft_gain: 80.0,
        }
    }

    /// Returns `true` if the lean angle exceeds the warning threshold.
    pub fn is_warning(&self, lean: f64) -> bool {
        lean.abs() > self.warn_threshold
    }

    /// Returns `true` if the lean angle exceeds the hard limit.
    pub fn is_over_limit(&self, lean: f64) -> bool {
        lean.abs() > self.hard_limit
    }

    /// Clamp lean angle to the hard limit.
    pub fn apply_hard_limit(&self, lean: f64) -> f64 {
        clamp(lean, -self.hard_limit, self.hard_limit)
    }

    /// Compute a soft-limit restoring torque that opposes lean in the danger zone.
    pub fn restoring_torque(&self, lean: f64) -> f64 {
        let abs_lean = lean.abs();
        if abs_lean < self.warn_threshold {
            return 0.0;
        }
        let excess = (abs_lean - self.warn_threshold).min(self.soft_zone);
        -sign(lean) * self.soft_gain * excess
    }
}

impl Default for LeanAngleLimiter {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// RiderBalanceController
// ---------------------------------------------------------------------------

/// PD lean-angle controller representing the rider's balance effort.
///
/// Applies a torque about the steering axis to maintain a target lean angle.
#[derive(Debug, Clone)]
pub struct RiderBalanceController {
    /// Proportional gain (N·m/rad).
    pub kp: f64,
    /// Derivative gain (N·m·s/rad).
    pub kd: f64,
    /// Target lean angle (rad).
    pub target_lean: f64,
    /// Previous lean error for derivative term.
    pub prev_error: f64,
    /// Maximum steering torque output (N·m).
    pub torque_limit: f64,
}

impl RiderBalanceController {
    /// Create a controller with typical sport-bike rider gains.
    pub fn new(kp: f64, kd: f64) -> Self {
        Self {
            kp,
            kd,
            target_lean: 0.0,
            prev_error: 0.0,
            torque_limit: 40.0,
        }
    }

    /// Compute and return the steering torque given current lean angle and `dt`.
    pub fn update(&mut self, current_lean: f64, dt: f64) -> f64 {
        let error = self.target_lean - current_lean;
        let d_error = if dt > 1e-9 {
            (error - self.prev_error) / dt
        } else {
            0.0
        };
        self.prev_error = error;
        let torque = self.kp * error + self.kd * d_error;
        clamp(torque, -self.torque_limit, self.torque_limit)
    }

    /// Set a new target lean angle in radians.
    pub fn set_target(&mut self, lean: f64) {
        self.target_lean = lean;
    }
}

// ---------------------------------------------------------------------------
// CounterSteeringModel
// ---------------------------------------------------------------------------

/// Counter-steering model: relates handlebar steering torque to the resulting
/// lean moment for a single-track vehicle.
///
/// When turning at speed, pushing the handlebar forward on one side (i.e.,
/// steering away from the intended lean direction) is required to initiate lean.
#[derive(Debug, Clone)]
pub struct CounterSteeringModel {
    /// Vehicle forward speed (m/s).
    pub speed: f64,
    /// Wheelbase (m).
    pub wheelbase: f64,
    /// Front wheel gyroscopic coefficient (kg·m²/s).
    pub gyro_coeff_front: f64,
    /// Rear wheel gyroscopic coefficient (kg·m²/s).
    pub gyro_coeff_rear: f64,
    /// Trail (m): distance from steering axis to contact patch along ground.
    pub trail: f64,
}

impl CounterSteeringModel {
    /// Create a counter-steering model for a typical sport bike.
    pub fn new() -> Self {
        Self {
            speed: 15.0,
            wheelbase: 1.43,
            gyro_coeff_front: 3.8,
            gyro_coeff_rear: 4.5,
            trail: 0.096,
        }
    }

    /// Compute the gyroscopic lean moment (N·m) for a given steering rate (rad/s).
    ///
    /// A positive steering rate (turning handlebars right) creates a rightward lean.
    pub fn gyroscopic_lean_moment(&self, steer_rate: f64) -> f64 {
        let gyro_total = self.gyro_coeff_front + self.gyro_coeff_rear;
        gyro_total * self.speed * steer_rate
    }

    /// Compute required steering torque to achieve a target lean rate (rad/s).
    pub fn required_steer_torque(&self, lean_rate: f64) -> f64 {
        let gyro_total = self.gyro_coeff_front + self.gyro_coeff_rear;
        if gyro_total * self.speed < 1e-6 {
            return 0.0;
        }
        lean_rate / (gyro_total * self.speed)
    }

    /// Compute equilibrium steering angle (rad) for steady circular turn of radius `r`.
    ///
    /// Positive steering angle = right turn.
    pub fn steady_state_steer(&self, radius: f64, lean_angle: f64) -> f64 {
        if radius.abs() < 0.01 {
            return 0.0;
        }
        // δ ≈ (1/R - (v² tan(φ))/(g L)) * L
        let v2 = self.speed * self.speed;
        let turn_angle = self.wheelbase / radius;
        let gravity_term = v2 * lean_angle.tan() / (G * self.wheelbase);
        turn_angle - gravity_term
    }
}

impl Default for CounterSteeringModel {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// WeelieStopiDetector
// ---------------------------------------------------------------------------

/// Wheelie and stoppie detection.
///
/// Monitors wheel normal forces to detect when a wheel lifts off the ground.
#[derive(Debug, Clone)]
pub struct WheelieStoppieDetector {
    /// Normal force threshold below which a wheel is considered lifted (N).
    pub lift_threshold: f64,
    /// Hysteresis: force must exceed this to confirm ground contact (N).
    pub contact_threshold: f64,
    /// `true` if the front wheel is currently lifted (stoppie state).
    pub front_lifted: bool,
    /// `true` if the rear wheel is currently lifted (wheelie state).
    pub rear_lifted: bool,
    /// Duration (seconds) for which the front has been lifted continuously.
    pub front_lift_duration: f64,
    /// Duration (seconds) for which the rear has been lifted continuously.
    pub rear_lift_duration: f64,
}

impl WheelieStoppieDetector {
    /// Create a new detector with default force thresholds.
    pub fn new() -> Self {
        Self {
            lift_threshold: 5.0,
            contact_threshold: 20.0,
            front_lifted: false,
            rear_lifted: false,
            front_lift_duration: 0.0,
            rear_lift_duration: 0.0,
        }
    }

    /// Update detector state given current front and rear normal forces.
    pub fn update(&mut self, front_normal: f64, rear_normal: f64, dt: f64) {
        // Front wheel
        if front_normal < self.lift_threshold {
            self.front_lifted = true;
            self.front_lift_duration += dt;
        } else if front_normal > self.contact_threshold {
            self.front_lifted = false;
            self.front_lift_duration = 0.0;
        }

        // Rear wheel
        if rear_normal < self.lift_threshold {
            self.rear_lifted = true;
            self.rear_lift_duration += dt;
        } else if rear_normal > self.contact_threshold {
            self.rear_lifted = false;
            self.rear_lift_duration = 0.0;
        }
    }

    /// Returns `true` if a wheelie (rear wheel on ground, front lifted) is active.
    pub fn is_wheelie(&self) -> bool {
        self.front_lifted && !self.rear_lifted
    }

    /// Returns `true` if a stoppie (front wheel on ground, rear lifted) is active.
    pub fn is_stoppie(&self) -> bool {
        self.rear_lifted && !self.front_lifted
    }
}

impl Default for WheelieStoppieDetector {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// TankSlapperAnalysis (linearised oscillation modes)
// ---------------------------------------------------------------------------

/// Linearised single-track motorcycle weave and wobble mode analysis.
///
/// The weave mode is a coupled lean-yaw-steer oscillation at ~2–4 Hz.
/// The wobble (tank-slapper) mode is a high-frequency steering oscillation
/// at ~7–12 Hz.
#[derive(Debug, Clone)]
pub struct TankSlapperAnalysis {
    /// Steering system natural frequency (rad/s).
    pub steer_natural_freq: f64,
    /// Steering damping ratio.
    pub steer_damping_ratio: f64,
    /// Weave mode natural frequency (rad/s).
    pub weave_natural_freq: f64,
    /// Weave mode damping ratio.
    pub weave_damping_ratio: f64,
    /// Current speed (m/s) — affects stability boundaries.
    pub speed: f64,
}

impl TankSlapperAnalysis {
    /// Create a default analysis for a typical 600 cc sport bike.
    pub fn new() -> Self {
        Self {
            steer_natural_freq: 2.0 * PI * 9.0, // ~9 Hz
            steer_damping_ratio: 0.15,
            weave_natural_freq: 2.0 * PI * 3.0, // ~3 Hz
            weave_damping_ratio: 0.12,
            speed: 30.0,
        }
    }

    /// Return wobble mode damped frequency (Hz).
    pub fn wobble_damped_freq(&self) -> f64 {
        let wd = self.steer_natural_freq
            * (1.0 - self.steer_damping_ratio * self.steer_damping_ratio)
                .max(0.0)
                .sqrt();
        wd / (2.0 * PI)
    }

    /// Return weave mode damped frequency (Hz).
    pub fn weave_damped_freq(&self) -> f64 {
        let wd = self.weave_natural_freq
            * (1.0 - self.weave_damping_ratio * self.weave_damping_ratio)
                .max(0.0)
                .sqrt();
        wd / (2.0 * PI)
    }

    /// Check if the steering system is underdamped (wobble risk).
    pub fn is_wobble_prone(&self) -> bool {
        self.steer_damping_ratio < 0.2
    }

    /// Check if the weave mode is underdamped (weave risk above crossover speed).
    pub fn is_weave_prone(&self, crossover_speed: f64) -> bool {
        self.speed > crossover_speed && self.weave_damping_ratio < 0.1
    }

    /// Compute approximate speed at which weave becomes unstable.
    ///
    /// Simple linearised model: instability speed ∝ ζ * ω_n.
    pub fn weave_instability_speed(&self) -> f64 {
        // Approximate: higher damping and lower natural frequency → more stable
        50.0 * self.weave_damping_ratio / (self.weave_natural_freq / (2.0 * PI))
    }
}

impl Default for TankSlapperAnalysis {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// MotorcycleState
// ---------------------------------------------------------------------------

/// Complete kinematic/dynamic state of the motorcycle at one instant.
#[derive(Debug, Clone)]
pub struct MotorcycleState {
    /// Forward speed (m/s).
    pub speed: f64,
    /// Lean angle from vertical (rad), positive = right.
    pub lean_angle: f64,
    /// Lean rate (rad/s).
    pub lean_rate: f64,
    /// Steering angle (rad), positive = right turn.
    pub steer_angle: f64,
    /// Steering rate (rad/s).
    pub steer_rate: f64,
    /// Yaw rate (rad/s).
    pub yaw_rate: f64,
    /// Front wheel angular velocity (rad/s).
    pub front_wheel_omega: f64,
    /// Rear wheel angular velocity (rad/s).
    pub rear_wheel_omega: f64,
    /// Front suspension compression (m).
    pub front_compression: f64,
    /// Rear suspension compression (m).
    pub rear_compression: f64,
    /// Front tyre normal force (N).
    pub front_normal: f64,
    /// Rear tyre normal force (N).
    pub rear_normal: f64,
    /// Pitch angle (nose up positive, rad).
    pub pitch_angle: f64,
    /// Pitch rate (rad/s).
    pub pitch_rate: f64,
}

impl MotorcycleState {
    /// Create a default upright, stationary state.
    pub fn new() -> Self {
        Self {
            speed: 0.0,
            lean_angle: 0.0,
            lean_rate: 0.0,
            steer_angle: 0.0,
            steer_rate: 0.0,
            yaw_rate: 0.0,
            front_wheel_omega: 0.0,
            rear_wheel_omega: 0.0,
            front_compression: 0.02,
            rear_compression: 0.02,
            front_normal: 0.0,
            rear_normal: 0.0,
            pitch_angle: 0.0,
            pitch_rate: 0.0,
        }
    }
}

impl Default for MotorcycleState {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// MotorcycleParameters
// ---------------------------------------------------------------------------

/// Physical parameters of the motorcycle.
#[derive(Debug, Clone)]
pub struct MotorcycleParameters {
    /// Total mass of bike + rider (kg).
    pub total_mass: f64,
    /// Height of centre of mass (m).
    pub com_height: f64,
    /// Wheelbase (m).
    pub wheelbase: f64,
    /// Distance from front axle to CoM (m).
    pub a: f64,
    /// Front wheel radius (m).
    pub front_wheel_radius: f64,
    /// Rear wheel radius (m).
    pub rear_wheel_radius: f64,
    /// Front wheel moment of inertia (kg·m²).
    pub front_wheel_inertia: f64,
    /// Rear wheel moment of inertia (kg·m²).
    pub rear_wheel_inertia: f64,
    /// Roll (lean) moment of inertia of the whole vehicle (kg·m²).
    pub roll_inertia: f64,
    /// Yaw moment of inertia (kg·m²).
    pub yaw_inertia: f64,
    /// Steering head angle from vertical (rad).
    pub rake_angle: f64,
    /// Fork offset / trail (m).
    pub trail: f64,
}

impl MotorcycleParameters {
    /// Create parameters for a typical 600 cc sport bike.
    pub fn sport_600() -> Self {
        Self {
            total_mass: 200.0,
            com_height: 0.60,
            wheelbase: 1.43,
            a: 0.65, // CoM closer to front
            front_wheel_radius: 0.305,
            rear_wheel_radius: 0.315,
            front_wheel_inertia: 0.6,
            rear_wheel_inertia: 0.8,
            roll_inertia: 25.0,
            yaw_inertia: 20.0,
            rake_angle: deg2rad(24.0),
            trail: 0.096,
        }
    }
}

// ---------------------------------------------------------------------------
// MotorcycleDynamicsSimulator
// ---------------------------------------------------------------------------

/// Full single-track motorcycle dynamics simulator.
///
/// Integrates the equations of motion for lean, yaw, steer, and longitudinal
/// dynamics using a semi-explicit Euler scheme.  Tyre forces use the Pacejka
/// Magic Formula.  Suspension is modelled as decoupled spring-dampers at each axle.
#[derive(Debug, Clone)]
pub struct MotorcycleDynamicsSimulator {
    /// Physical parameters.
    pub params: MotorcycleParameters,
    /// Current state.
    pub state: MotorcycleState,
    /// Front Pacejka tyre.
    pub front_tire: PacejkaMotoTire,
    /// Rear Pacejka tyre.
    pub rear_tire: PacejkaMotoTire,
    /// Telescopic fork front suspension.
    pub front_suspension: TelescopicFork,
    /// Monoshock rear suspension.
    pub rear_suspension: Monoshock,
    /// Aerodynamic fairing.
    pub fairing: AerodynamicFairing,
    /// Lean-angle limiter.
    pub lean_limiter: LeanAngleLimiter,
    /// Rider balance controller.
    pub rider: RiderBalanceController,
    /// Counter-steering model.
    pub counter_steer: CounterSteeringModel,
    /// Wheelie/stoppie detector.
    pub wheel_detector: WheelieStoppieDetector,
    /// Chain drive.
    pub chain: ChainDrive,
    /// Engine torque applied at the rear wheel (N·m).
    pub engine_torque: f64,
    /// Front brake torque (N·m, positive = braking).
    pub front_brake_torque: f64,
    /// Rear brake torque (N·m, positive = braking).
    pub rear_brake_torque: f64,
    /// Accumulated engine shaft angle (rad).
    pub engine_angle: f64,
    /// Accumulated rear wheel angle (rad).
    pub rear_wheel_angle: f64,
}

impl MotorcycleDynamicsSimulator {
    /// Create a new simulator with sport-bike defaults.
    pub fn new() -> Self {
        let params = MotorcycleParameters::sport_600();
        let mut state = MotorcycleState::new();
        // Initialise normal forces from static weight distribution
        let m = params.total_mass;
        state.rear_normal = m * G * params.a / params.wheelbase;
        state.front_normal = m * G * (params.wheelbase - params.a) / params.wheelbase;
        let mut counter_steer = CounterSteeringModel::new();
        counter_steer.speed = state.speed;
        counter_steer.wheelbase = params.wheelbase;
        counter_steer.trail = params.trail;

        Self {
            params,
            state,
            front_tire: PacejkaMotoTire::new(),
            rear_tire: PacejkaMotoTire::new(),
            front_suspension: TelescopicFork::new(),
            rear_suspension: Monoshock::new(),
            fairing: AerodynamicFairing::new_sport_bike(),
            lean_limiter: LeanAngleLimiter::new(),
            rider: RiderBalanceController::new(60.0, 8.0),
            counter_steer,
            wheel_detector: WheelieStoppieDetector::new(),
            chain: ChainDrive::new(15, 42),
            engine_torque: 0.0,
            front_brake_torque: 0.0,
            rear_brake_torque: 0.0,
            engine_angle: 0.0,
            rear_wheel_angle: 0.0,
        }
    }

    /// Set engine throttle torque (N·m).
    pub fn set_engine_torque(&mut self, torque: f64) {
        self.engine_torque = torque.max(0.0);
    }

    /// Set front and rear brake torques (N·m).
    pub fn set_brakes(&mut self, front: f64, rear: f64) {
        self.front_brake_torque = front.max(0.0);
        self.rear_brake_torque = rear.max(0.0);
    }

    /// Update normal forces based on current acceleration and pitch.
    fn update_normal_forces(&mut self, longitudinal_accel: f64) {
        let m = self.params.total_mass;
        let h = self.params.com_height;
        let l = self.params.wheelbase;
        let a = self.params.a;

        // Weight transfer due to acceleration
        let delta_fz = m * longitudinal_accel * h / l;

        let static_rear = m * G * a / l;
        let static_front = m * G * (l - a) / l;

        // Aero downforce split
        let v = self.state.speed;
        let aero_down = self.fairing.downforce(v, self.state.pitch_angle);
        let aero_front = aero_down * 0.45;
        let aero_rear = aero_down * 0.55;

        self.state.front_normal = (static_front - delta_fz + aero_front).max(0.0);
        self.state.rear_normal = (static_rear + delta_fz + aero_rear).max(0.0);
    }

    /// Advance the simulation by one time step `dt` (seconds).
    pub fn step(&mut self, dt: f64) {
        let v = self.state.speed;
        let phi = self.state.lean_angle;
        let _dphi = self.state.lean_rate;
        let delta = self.state.steer_angle;

        // --- Tyre slip ---
        let rr = self.params.rear_wheel_radius;
        let rf = self.params.front_wheel_radius;

        let rear_slip_ratio = if v.abs() > 0.1 {
            (self.state.rear_wheel_omega * rr - v) / v.abs()
        } else {
            0.0
        };
        let front_slip_ratio = if v.abs() > 0.1 {
            (self.state.front_wheel_omega * rf - v) / v.abs()
        } else {
            0.0
        };

        // Lateral slip angles (small angle model)
        let rear_slip_angle = if v.abs() > 0.1 {
            -(self.state.yaw_rate * self.params.a / v)
        } else {
            0.0
        };
        let b = self.params.wheelbase - self.params.a;
        let front_slip_angle = if v.abs() > 0.1 {
            delta - (self.state.yaw_rate * b / v)
        } else {
            0.0
        };

        // Update relaxation filters
        self.front_tire.update_relaxation(front_slip_angle, v, dt);
        self.rear_tire.update_relaxation(rear_slip_angle, v, dt);

        // --- Tyre forces ---
        let (rear_fy, rear_fx, _rear_mz) = self.rear_tire.forces(
            self.rear_tire.filtered_slip_angle,
            rear_slip_ratio,
            self.state.rear_normal,
            phi,
        );
        let (front_fy, front_fx, _front_mz) = self.front_tire.forces(
            self.front_tire.filtered_slip_angle,
            front_slip_ratio,
            self.state.front_normal,
            phi,
        );

        // --- Longitudinal dynamics ---
        let drive_torque = self.chain.transmitted_torque(self.engine_torque, 1.0);
        let drive_force = drive_torque / rr;
        let brake_force = (self.rear_brake_torque / rr) + (self.front_brake_torque / rf);
        let aero_drag = self.fairing.drag_force(v, self.state.pitch_angle);
        let long_force = drive_force + rear_fx + front_fx - brake_force - aero_drag;
        let long_accel = long_force / self.params.total_mass;

        // --- Lean (roll) dynamics ---
        let gravity_lean = self.params.total_mass * G * self.params.com_height * phi.sin();
        let centrifugal_lean =
            -self.params.total_mass * self.state.yaw_rate * v * self.params.com_height * phi.cos();
        let gyro_lean = self
            .counter_steer
            .gyroscopic_lean_moment(self.state.steer_rate);
        let rider_torque = self.rider.update(phi, dt);
        let limiter_torque = self.lean_limiter.restoring_torque(phi);
        let net_lean_torque =
            gravity_lean + centrifugal_lean + gyro_lean + rider_torque + limiter_torque;
        let lean_accel = net_lean_torque / self.params.roll_inertia;

        // --- Yaw dynamics ---
        let net_yaw_torque = (rear_fy * self.params.a) - (front_fy * b);
        let yaw_accel = net_yaw_torque / self.params.yaw_inertia;

        // --- Wheel angular accelerations ---
        let rear_wheel_alpha = (drive_torque - self.rear_brake_torque - rear_fx.abs() * rr)
            / self.params.rear_wheel_inertia;
        let front_wheel_alpha =
            (-self.front_brake_torque - front_fx.abs() * rf) / self.params.front_wheel_inertia;

        // --- Suspension step ---
        let road_vel = 0.0; // flat ground
        self.front_suspension.step(road_vel, dt);
        self.rear_suspension.step(road_vel, dt);

        // --- Chain dynamics ---
        self.engine_angle += dt * self.state.rear_wheel_omega * self.chain.sprocket_ratio();
        self.rear_wheel_angle += dt * self.state.rear_wheel_omega;
        self.chain
            .step(self.engine_angle, self.rear_wheel_angle, dt);

        // --- Update normal forces ---
        self.update_normal_forces(long_accel);

        // --- Wheelie / stoppie detection ---
        self.wheel_detector
            .update(self.state.front_normal, self.state.rear_normal, dt);

        // --- Integrate state (semi-explicit Euler) ---
        self.state.speed += dt * long_accel;
        self.state.speed = self.state.speed.max(0.0);

        self.state.lean_rate += dt * lean_accel;
        self.state.lean_angle += dt * self.state.lean_rate;
        self.state.lean_angle = self.lean_limiter.apply_hard_limit(self.state.lean_angle);

        self.state.yaw_rate += dt * yaw_accel;

        self.state.rear_wheel_omega += dt * rear_wheel_alpha;
        self.state.rear_wheel_omega = self.state.rear_wheel_omega.max(0.0);
        self.state.front_wheel_omega += dt * front_wheel_alpha;
        self.state.front_wheel_omega = self.state.front_wheel_omega.max(0.0);

        self.state.front_compression = self.front_suspension.compression;
        self.state.rear_compression = self.rear_suspension.compression;

        // Update counter-steer model speed
        self.counter_steer.speed = self.state.speed;
    }
}

impl Default for MotorcycleDynamicsSimulator {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// GyroscopicEffect
// ---------------------------------------------------------------------------

/// Gyroscopic effect calculator for rotating wheels.
///
/// Computes the gyroscopic precession torque that arises when the spin
/// axis of a wheel is rotated by a yaw or lean rate.
#[derive(Debug, Clone)]
pub struct GyroscopicEffect {
    /// Wheel moment of inertia about the spin axis (kg·m²).
    pub wheel_inertia: f64,
    /// Wheel angular velocity (rad/s).
    pub omega: f64,
}

impl GyroscopicEffect {
    /// Create a new gyroscopic calculator.
    pub fn new(wheel_inertia: f64, omega: f64) -> Self {
        Self {
            wheel_inertia,
            omega,
        }
    }

    /// Gyroscopic angular momentum magnitude (kg·m²/s).
    pub fn angular_momentum(&self) -> f64 {
        self.wheel_inertia * self.omega
    }

    /// Gyroscopic torque (N·m) when the spin axis precesses at `precession_rate` (rad/s).
    ///
    /// `T = H × ψ̇` where `H` is angular momentum and `ψ̇` is the precession rate.
    pub fn torque(&self, precession_rate: f64) -> f64 {
        self.angular_momentum() * precession_rate
    }

    /// Combined gyroscopic lean torque for both wheels.
    pub fn combined_lean_torque(
        front: &GyroscopicEffect,
        rear: &GyroscopicEffect,
        yaw_rate: f64,
    ) -> f64 {
        (front.angular_momentum() + rear.angular_momentum()) * yaw_rate
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // PacejkaMotoTire tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_pacejka_evaluate_zero() {
        let c = PacejkaMotoCoeffs::default();
        let f = c.evaluate(0.0);
        assert!(
            f.abs() < 1e-10,
            "Pacejka at zero slip should be zero, got {f}"
        );
    }

    #[test]
    fn test_pacejka_evaluate_positive_slip() {
        let c = PacejkaMotoCoeffs::default();
        let f = c.evaluate(0.1);
        assert!(f > 0.0, "Pacejka should be positive for positive slip");
    }

    #[test]
    fn test_pacejka_antisymmetric() {
        let c = PacejkaMotoCoeffs::default();
        let f_pos = c.evaluate(0.2);
        let f_neg = c.evaluate(-0.2);
        assert!(
            (f_pos + f_neg).abs() < 1e-10,
            "Pacejka should be antisymmetric: {f_pos} + {f_neg}"
        );
    }

    #[test]
    fn test_tire_forces_zero_slip() {
        let tire = PacejkaMotoTire::new();
        let (fy, fx, mz) = tire.forces(0.0, 0.0, 1000.0, 0.0);
        assert!(fy.abs() < 1e-6, "Zero slip -> zero lateral force, got {fy}");
        assert!(
            fx.abs() < 1e-6,
            "Zero slip -> zero longitudinal force, got {fx}"
        );
        assert!(
            mz.abs() < 1e-6,
            "Zero lateral force -> zero self-aligning torque, got {mz}"
        );
    }

    #[test]
    fn test_tire_forces_positive_slip_angle() {
        let tire = PacejkaMotoTire::new();
        let (fy, _fx, _mz) = tire.forces(0.1, 0.0, 2000.0, 0.0);
        assert!(
            fy > 0.0,
            "Positive slip angle should give positive lateral force"
        );
    }

    #[test]
    fn test_tire_relaxation_converges() {
        let mut tire = PacejkaMotoTire::new();
        let raw = 0.15;
        for _ in 0..100 {
            tire.update_relaxation(raw, 20.0, 0.01);
        }
        assert!(
            (tire.filtered_slip_angle - raw).abs() < 0.001,
            "Relaxation should converge to raw value"
        );
    }

    // -----------------------------------------------------------------------
    // TelescopicFork tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_fork_force_positive_compression() {
        let fork = TelescopicFork::new();
        let f = fork.force(0.02, 0.0);
        assert!(f > 0.0, "Compressed fork should push back");
    }

    #[test]
    fn test_fork_force_zero_compression() {
        let fork = TelescopicFork::new();
        // Even zero compression gives preload force
        let f = fork.force(0.0, 0.0);
        assert!(f >= 0.0, "Fork force should be non-negative");
    }

    #[test]
    fn test_fork_travel_clamped() {
        let mut fork = TelescopicFork::new();
        fork.step(100.0, 1.0); // very large velocity for 1 second
        assert!(
            fork.compression <= fork.max_travel,
            "Fork travel should not exceed max_travel"
        );
    }

    // -----------------------------------------------------------------------
    // Monoshock tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_monoshock_wheel_rate() {
        let shock = Monoshock::new();
        let wr = shock.wheel_rate();
        let expected = shock.spring_rate * shock.motion_ratio * shock.motion_ratio;
        assert!((wr - expected).abs() < 1e-9);
    }

    #[test]
    fn test_monoshock_wheel_force_positive() {
        let shock = Monoshock::new();
        let f = shock.wheel_force(0.02, 0.0);
        assert!(
            f > 0.0,
            "Compressed monoshock should give positive wheel force"
        );
    }

    // -----------------------------------------------------------------------
    // ChainDrive tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_chain_sprocket_ratio() {
        let chain = ChainDrive::new(15, 42);
        let ratio = chain.sprocket_ratio();
        assert!((ratio - 42.0 / 15.0).abs() < 1e-10);
    }

    #[test]
    fn test_chain_transmitted_torque() {
        let chain = ChainDrive::new(15, 42);
        let t = chain.transmitted_torque(50.0, 2.0);
        let expected = 50.0 * 2.0 * chain.sprocket_ratio();
        assert!((t - expected).abs() < 1e-9);
    }

    #[test]
    fn test_chain_tension() {
        let chain = ChainDrive::new(15, 42);
        let tension = chain.chain_tension(280.0, 0.10);
        assert!((tension - 2800.0).abs() < 1e-9);
    }

    // -----------------------------------------------------------------------
    // AerodynamicFairing tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_fairing_drag_positive() {
        let fairing = AerodynamicFairing::new_sport_bike();
        let drag = fairing.drag_force(30.0, 0.0);
        assert!(drag > 0.0, "Drag should be positive at non-zero speed");
    }

    #[test]
    fn test_fairing_drag_zero_speed() {
        let fairing = AerodynamicFairing::new_sport_bike();
        let drag = fairing.drag_force(0.0, 0.0);
        assert!(drag.abs() < 1e-12, "Zero speed gives zero drag");
    }

    #[test]
    fn test_fairing_drag_scales_v_squared() {
        let fairing = AerodynamicFairing::new_sport_bike();
        let d1 = fairing.drag_force(10.0, 0.0);
        let d2 = fairing.drag_force(20.0, 0.0);
        assert!(
            (d2 / d1 - 4.0).abs() < 0.01,
            "Drag should scale with v², got ratio {}",
            d2 / d1
        );
    }

    #[test]
    fn test_fairing_downforce() {
        let fairing = AerodynamicFairing::new_sport_bike();
        // cl_alpha is negative so positive aoa gives positive downforce
        let df = fairing.downforce(40.0, 0.05);
        assert!(df > 0.0, "Should produce downforce at positive aoa");
    }

    // -----------------------------------------------------------------------
    // LeanAngleLimiter tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_limiter_hard_limit_clamps() {
        let limiter = LeanAngleLimiter::new();
        let extreme = deg2rad(80.0);
        let clamped = limiter.apply_hard_limit(extreme);
        assert!(
            clamped <= limiter.hard_limit + 1e-10,
            "Hard limit should clamp lean angle"
        );
    }

    #[test]
    fn test_limiter_no_warning_below_threshold() {
        let limiter = LeanAngleLimiter::new();
        assert!(!limiter.is_warning(deg2rad(30.0)));
    }

    #[test]
    fn test_limiter_warning_above_threshold() {
        let limiter = LeanAngleLimiter::new();
        assert!(limiter.is_warning(deg2rad(55.0)));
    }

    #[test]
    fn test_limiter_restoring_torque_direction() {
        let limiter = LeanAngleLimiter::new();
        let torque = limiter.restoring_torque(deg2rad(55.0));
        assert!(torque < 0.0, "Restoring torque should oppose positive lean");
    }

    #[test]
    fn test_limiter_no_restoring_below_warn() {
        let limiter = LeanAngleLimiter::new();
        let torque = limiter.restoring_torque(deg2rad(20.0));
        assert!(
            torque.abs() < 1e-10,
            "No restoring torque below warn threshold"
        );
    }

    // -----------------------------------------------------------------------
    // RiderBalanceController tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_rider_p_term_positive_error() {
        let mut ctrl = RiderBalanceController::new(50.0, 5.0);
        ctrl.set_target(0.0);
        // If current lean is negative, error = 0 - (-0.1) = +0.1, torque > 0
        let torque = ctrl.update(-0.1, 0.01);
        assert!(torque > 0.0, "Positive error should give positive torque");
    }

    #[test]
    fn test_rider_torque_clamped() {
        let mut ctrl = RiderBalanceController::new(1000.0, 100.0);
        ctrl.set_target(0.0);
        let torque = ctrl.update(-1.0, 0.01);
        assert!(
            torque <= ctrl.torque_limit,
            "Torque should be clamped to torque_limit"
        );
    }

    #[test]
    fn test_rider_zero_error_zero_torque() {
        let mut ctrl = RiderBalanceController::new(50.0, 5.0);
        ctrl.set_target(0.1);
        // First call sets prev_error; second call with same lean gives d=0
        ctrl.update(0.1, 0.01);
        let torque = ctrl.update(0.1, 0.01);
        assert!(
            torque.abs() < 1e-9,
            "Zero error should give zero torque (steady state)"
        );
    }

    // -----------------------------------------------------------------------
    // CounterSteeringModel tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_counter_steer_gyro_moment_nonzero() {
        let mut cs = CounterSteeringModel::new();
        cs.speed = 20.0;
        let moment = cs.gyroscopic_lean_moment(0.5);
        assert!(
            moment.abs() > 0.0,
            "Non-zero steer rate should give gyroscopic moment"
        );
    }

    #[test]
    fn test_counter_steer_zero_speed() {
        let mut cs = CounterSteeringModel::new();
        cs.speed = 0.0;
        let torque = cs.required_steer_torque(0.5);
        assert_eq!(torque, 0.0, "Zero speed gives zero steering torque");
    }

    #[test]
    fn test_counter_steer_steady_state() {
        let cs = CounterSteeringModel::new();
        let delta = cs.steady_state_steer(50.0, deg2rad(20.0));
        assert!(delta.is_finite(), "Steady-state steer should be finite");
    }

    // -----------------------------------------------------------------------
    // WheelieStoppieDetector tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_detector_initial_no_wheelie() {
        let det = WheelieStoppieDetector::new();
        assert!(!det.is_wheelie());
        assert!(!det.is_stoppie());
    }

    #[test]
    fn test_detector_wheelie_detection() {
        let mut det = WheelieStoppieDetector::new();
        // Front lifts (very low normal force), rear stays down
        det.update(1.0, 1000.0, 0.1);
        assert!(det.is_wheelie(), "Should detect wheelie");
    }

    #[test]
    fn test_detector_stoppie_detection() {
        let mut det = WheelieStoppieDetector::new();
        // Rear lifts, front stays down
        det.update(1000.0, 1.0, 0.1);
        assert!(det.is_stoppie(), "Should detect stoppie");
    }

    #[test]
    fn test_detector_resets_on_contact() {
        let mut det = WheelieStoppieDetector::new();
        det.update(1.0, 1000.0, 0.1); // trigger wheelie
        det.update(500.0, 1000.0, 0.1); // front comes back down
        assert!(
            !det.front_lifted,
            "Front lift should clear when normal force recovers"
        );
    }

    // -----------------------------------------------------------------------
    // TankSlapperAnalysis tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_wobble_freq_positive() {
        let analysis = TankSlapperAnalysis::new();
        let freq = analysis.wobble_damped_freq();
        assert!(freq > 0.0, "Wobble frequency should be positive");
    }

    #[test]
    fn test_weave_freq_positive() {
        let analysis = TankSlapperAnalysis::new();
        let freq = analysis.weave_damped_freq();
        assert!(freq > 0.0, "Weave frequency should be positive");
    }

    #[test]
    fn test_wobble_prone_low_damping() {
        let mut analysis = TankSlapperAnalysis::new();
        analysis.steer_damping_ratio = 0.1;
        assert!(
            analysis.is_wobble_prone(),
            "Low damping should be wobble-prone"
        );
    }

    #[test]
    fn test_wobble_not_prone_high_damping() {
        let mut analysis = TankSlapperAnalysis::new();
        analysis.steer_damping_ratio = 0.5;
        assert!(
            !analysis.is_wobble_prone(),
            "High damping should not be wobble-prone"
        );
    }

    // -----------------------------------------------------------------------
    // MotorcycleDynamicsSimulator tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_simulator_creation_normal_forces() {
        let sim = MotorcycleDynamicsSimulator::new();
        let m = sim.params.total_mass;
        let total = sim.state.front_normal + sim.state.rear_normal;
        assert!(
            (total - m * G).abs() < 1.0,
            "Initial normal forces should sum to weight: expected {}, got {}",
            m * G,
            total
        );
    }

    #[test]
    fn test_simulator_step_no_panic() {
        let mut sim = MotorcycleDynamicsSimulator::new();
        sim.set_engine_torque(80.0);
        for _ in 0..100 {
            sim.step(0.01);
        }
        assert!(sim.state.speed.is_finite(), "Speed should remain finite");
        assert!(
            sim.state.lean_angle.is_finite(),
            "Lean angle should remain finite"
        );
    }

    #[test]
    fn test_simulator_accelerates_with_throttle() {
        let mut sim = MotorcycleDynamicsSimulator::new();
        sim.state.rear_wheel_omega = 100.0;
        sim.state.front_wheel_omega = 100.0;
        sim.state.speed = 20.0;
        sim.set_engine_torque(100.0);
        let initial_speed = sim.state.speed;
        for _ in 0..50 {
            sim.step(0.01);
        }
        assert!(
            sim.state.speed >= initial_speed,
            "Throttle should not decelerate the bike"
        );
    }

    #[test]
    fn test_simulator_lean_limited() {
        let mut sim = MotorcycleDynamicsSimulator::new();
        // Force a large lean
        sim.state.lean_angle = deg2rad(70.0);
        sim.step(0.01);
        assert!(
            sim.state.lean_angle.abs() <= sim.lean_limiter.hard_limit + 1e-9,
            "Lean angle should not exceed hard limit"
        );
    }

    // -----------------------------------------------------------------------
    // GyroscopicEffect tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_gyro_angular_momentum() {
        let g = GyroscopicEffect::new(0.8, 200.0);
        assert!((g.angular_momentum() - 160.0).abs() < 1e-9);
    }

    #[test]
    fn test_gyro_torque() {
        let g = GyroscopicEffect::new(0.8, 200.0);
        let t = g.torque(0.5);
        assert!((t - 80.0).abs() < 1e-9);
    }

    #[test]
    fn test_gyro_combined() {
        let front = GyroscopicEffect::new(0.6, 150.0);
        let rear = GyroscopicEffect::new(0.8, 150.0);
        let t = GyroscopicEffect::combined_lean_torque(&front, &rear, 0.3);
        let expected = (0.6 * 150.0 + 0.8 * 150.0) * 0.3;
        assert!((t - expected).abs() < 1e-9);
    }

    // -----------------------------------------------------------------------
    // Utility function tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_deg2rad_90() {
        assert!((deg2rad(90.0) - PI / 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_rad2deg_pi() {
        assert!((rad2deg(PI) - 180.0).abs() < 1e-10);
    }

    #[test]
    fn test_clamp() {
        assert_eq!(clamp(5.0, 0.0, 3.0), 3.0);
        assert_eq!(clamp(-1.0, 0.0, 3.0), 0.0);
        assert_eq!(clamp(2.0, 0.0, 3.0), 2.0);
    }

    #[test]
    fn test_sign() {
        assert_eq!(sign(3.5), 1.0);
        assert_eq!(sign(-2.1), -1.0);
        assert_eq!(sign(0.0), 0.0);
    }
}
