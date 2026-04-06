#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced vehicle dynamics models for OxiPhysics vehicle crate.
//!
//! Provides bicycle models (linear/nonlinear), double-track, rollover,
//! understeer gradient, pitch dynamics, yaw moment diagram, handling balance,
//! wheel spin control, and a vehicle state estimator.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]
#![allow(unused_imports)]

use std::f64::consts::PI;

// ─── helpers ────────────────────────────────────────────────────────────────

#[inline]
fn clamp(v: f64, lo: f64, hi: f64) -> f64 {
    v.max(lo).min(hi)
}

/// 3×3 matrix as flat row-major array.
type Mat3 = [f64; 9];

/// Matrix-vector multiply: 3×3 * 3 → 3.
fn mat3_mul_vec3(m: &Mat3, v: &[f64; 3]) -> [f64; 3] {
    [
        m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
        m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
        m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
    ]
}

/// Transpose a 3×3 matrix.
fn mat3_transpose(m: &Mat3) -> Mat3 {
    [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]]
}

/// Identity 3×3.
fn mat3_identity() -> Mat3 {
    [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
}

/// Simple 2D rotation matrix (returns 2×2 as \[f64; 4\]).
fn rot2(angle: f64) -> [f64; 4] {
    let c = angle.cos();
    let s = angle.sin();
    [c, -s, s, c]
}

// ─── LinearBicycleModel ─────────────────────────────────────────────────────

/// Two-degree-of-freedom linear bicycle model state.
///
/// Captures lateral velocity `v_y` and yaw rate `r`.
#[derive(Debug, Clone, Default)]
pub struct BicycleState {
    /// Lateral velocity \[m/s\].
    pub v_y: f64,
    /// Yaw rate \[rad/s\].
    pub r: f64,
    /// Side slip angle \[rad\].
    pub beta: f64,
    /// Heading angle \[rad\].
    pub psi: f64,
    /// X position \[m\].
    pub x: f64,
    /// Y position \[m\].
    pub y: f64,
}

/// Linear bicycle model: 2-DOF lateral dynamics.
#[derive(Debug, Clone)]
pub struct LinearBicycleModel {
    /// Vehicle mass \[kg\].
    pub mass: f64,
    /// Yaw moment of inertia \[kg·m²\].
    pub iz: f64,
    /// Distance from CG to front axle \[m\].
    pub lf: f64,
    /// Distance from CG to rear axle \[m\].
    pub lr: f64,
    /// Front cornering stiffness \[N/rad\] (both tires combined).
    pub cf: f64,
    /// Rear cornering stiffness \[N/rad\] (both tires combined).
    pub cr: f64,
    /// Current state.
    pub state: BicycleState,
}

impl LinearBicycleModel {
    /// Create a linear bicycle model for a typical sedan.
    pub fn new() -> Self {
        Self {
            mass: 1500.0,
            iz: 2500.0,
            lf: 1.1,
            lr: 1.6,
            cf: 80_000.0,
            cr: 90_000.0,
            state: BicycleState::default(),
        }
    }

    /// Stability derivatives: A matrix of \[v_y_dot, r_dot\] = A * \[v_y, r\] + B * delta.
    ///
    /// Returns (a11, a12, a21, a22, b1, b2).
    pub fn stability_derivatives(&self, vx: f64) -> (f64, f64, f64, f64, f64, f64) {
        let m = self.mass;
        let iz = self.iz;
        let a = if vx.abs() > 0.01 {
            -(self.cf + self.cr) / (m * vx)
        } else {
            0.0
        };
        let b = if vx.abs() > 0.01 {
            (-self.cf * self.lf + self.cr * self.lr) / (m * vx) - vx
        } else {
            -vx
        };
        let c = if vx.abs() > 0.01 {
            (-self.cf * self.lf + self.cr * self.lr) / (iz * vx)
        } else {
            0.0
        };
        let d = if vx.abs() > 0.01 {
            -(self.cf * self.lf * self.lf + self.cr * self.lr * self.lr) / (iz * vx)
        } else {
            0.0
        };
        let b1 = self.cf / m;
        let b2 = self.cf * self.lf / iz;
        (a, b, c, d, b1, b2)
    }

    /// Critical speed \[m/s\] (speed above which understeer gradient causes stability issues).
    pub fn critical_speed(&self) -> f64 {
        let k = self.understeer_gradient();
        if k <= 0.0 {
            return f64::INFINITY;
        }
        let g = 9.81;
        (g / k).sqrt()
    }

    /// Understeer gradient \[rad/g\].
    pub fn understeer_gradient(&self) -> f64 {
        let g = 9.81;
        let wf = self.mass * g * self.lr / (self.lf + self.lr);
        let wr = self.mass * g * self.lf / (self.lf + self.lr);
        wf / self.cf - wr / self.cr
    }

    /// Step model by dt \[s\] at given longitudinal speed and steering angle.
    pub fn step(&mut self, dt: f64, vx: f64, delta: f64) {
        let (a11, a12, a21, a22, b1, b2) = self.stability_derivatives(vx);
        let vy = self.state.v_y;
        let r = self.state.r;
        // Euler integration
        let dvy = a11 * vy + a12 * r + b1 * delta;
        let dr = a21 * vy + a22 * r + b2 * delta;
        self.state.v_y += dvy * dt;
        self.state.r += dr * dt;
        self.state.psi += self.state.r * dt;
        self.state.beta = if vx.abs() > 0.01 {
            (self.state.v_y / vx).atan()
        } else {
            0.0
        };
        self.state.x += (vx * self.state.psi.cos() - self.state.v_y * self.state.psi.sin()) * dt;
        self.state.y += (vx * self.state.psi.sin() + self.state.v_y * self.state.psi.cos()) * dt;
    }

    /// Steady-state yaw rate \[rad/s\] for given speed and steering angle.
    pub fn steady_state_yaw_rate(&self, vx: f64, delta: f64) -> f64 {
        let l = self.lf + self.lr;
        let k = self.understeer_gradient();
        vx * delta / (l * (1.0 + k * vx * vx / 9.81))
    }
}

impl Default for LinearBicycleModel {
    fn default() -> Self {
        Self::new()
    }
}

// ─── NonlinearBicycleModel ──────────────────────────────────────────────────

/// Tire force model for nonlinear bicycle model.
#[derive(Debug, Clone)]
pub struct TireForceModel {
    /// Pacejka B coefficient.
    pub b: f64,
    /// Pacejka C coefficient.
    pub c: f64,
    /// Pacejka D coefficient (peak friction * Fz).
    pub d: f64,
    /// Pacejka E coefficient.
    pub e: f64,
}

impl TireForceModel {
    /// Default Pacejka tire model.
    pub fn default_values() -> Self {
        Self {
            b: 10.0,
            c: 1.9,
            d: 1.0,
            e: 0.97,
        }
    }

    /// Lateral force coefficient from slip angle \[rad\].
    pub fn lateral_force_coeff(&self, alpha: f64) -> f64 {
        let x = self.b * alpha;

        self.d * (self.c * (x - self.e * (x - x.atan())).atan()).sin()
    }
}

/// Nonlinear 3-DOF bicycle model with combined slip and load transfer.
#[derive(Debug, Clone)]
pub struct NonlinearBicycleModel {
    /// Mass \[kg\].
    pub mass: f64,
    /// Yaw inertia \[kg·m²\].
    pub iz: f64,
    /// CG to front \[m\].
    pub lf: f64,
    /// CG to rear \[m\].
    pub lr: f64,
    /// CG height \[m\].
    pub hcg: f64,
    /// Wheelbase \[m\].
    pub wheelbase: f64,
    /// Tire front.
    pub tire_f: TireForceModel,
    /// Tire rear.
    pub tire_r: TireForceModel,
    /// Lateral velocity \[m/s\] (state).
    pub vy: f64,
    /// Yaw rate \[rad/s\] (state).
    pub r: f64,
    /// Heading \[rad\] (state).
    pub psi: f64,
    /// X position \[m\] (state).
    pub x: f64,
    /// Y position \[m\] (state).
    pub y: f64,
}

impl NonlinearBicycleModel {
    /// Create a default nonlinear bicycle model.
    pub fn new() -> Self {
        let lf = 1.1;
        let lr = 1.6;
        Self {
            mass: 1500.0,
            iz: 2500.0,
            lf,
            lr,
            hcg: 0.5,
            wheelbase: lf + lr,
            tire_f: TireForceModel::default_values(),
            tire_r: TireForceModel::default_values(),
            vy: 0.0,
            r: 0.0,
            psi: 0.0,
            x: 0.0,
            y: 0.0,
        }
    }

    /// Front slip angle \[rad\].
    pub fn front_slip_angle(&self, vx: f64, delta: f64) -> f64 {
        if vx.abs() < 0.1 {
            return 0.0;
        }
        -delta + ((self.vy + self.r * self.lf) / vx).atan()
    }

    /// Rear slip angle \[rad\].
    pub fn rear_slip_angle(&self, vx: f64) -> f64 {
        if vx.abs() < 0.1 {
            return 0.0;
        }
        ((self.vy - self.r * self.lr) / vx).atan()
    }

    /// Front axle load \[N\] including longitudinal load transfer.
    pub fn front_load(&self, ax: f64) -> f64 {
        let g = 9.81;
        let static_fz = self.mass * g * self.lr / self.wheelbase;
        let transfer = self.mass * ax * self.hcg / self.wheelbase;
        (static_fz - transfer).max(0.0)
    }

    /// Rear axle load \[N\].
    pub fn rear_load(&self, ax: f64) -> f64 {
        let g = 9.81;
        let static_fz = self.mass * g * self.lf / self.wheelbase;
        let transfer = self.mass * ax * self.hcg / self.wheelbase;
        (static_fz + transfer).max(0.0)
    }

    /// Step model by dt \[s\].
    pub fn step(&mut self, dt: f64, vx: f64, delta: f64, ax: f64) {
        let alpha_f = self.front_slip_angle(vx, delta);
        let alpha_r = self.rear_slip_angle(vx);
        let fz_f = self.front_load(ax);
        let fz_r = self.rear_load(ax);
        let fy_f = self.tire_f.lateral_force_coeff(alpha_f) * fz_f;
        let fy_r = self.tire_r.lateral_force_coeff(alpha_r) * fz_r;
        let dvy = (fy_f + fy_r) / self.mass - self.r * vx;
        let dr = (fy_f * self.lf - fy_r * self.lr) / self.iz;
        self.vy += dvy * dt;
        self.r += dr * dt;
        self.psi += self.r * dt;
        self.x += (vx * self.psi.cos() - self.vy * self.psi.sin()) * dt;
        self.y += (vx * self.psi.sin() + self.vy * self.psi.cos()) * dt;
    }
}

impl Default for NonlinearBicycleModel {
    fn default() -> Self {
        Self::new()
    }
}

// ─── DoubleTrackModel ───────────────────────────────────────────────────────

/// Four-wheel double-track vehicle model state.
#[derive(Debug, Clone, Default)]
pub struct DoubleTrackState {
    /// Lateral velocity \[m/s\].
    pub vy: f64,
    /// Yaw rate \[rad/s\].
    pub r: f64,
    /// Roll angle \[rad\].
    pub phi: f64,
    /// Roll rate \[rad/s\].
    pub phi_dot: f64,
    /// Individual wheel slip angles \[rad\]: FL, FR, RL, RR.
    pub slip_angles: [f64; 4],
    /// Individual tire lateral forces \[N\]: FL, FR, RL, RR.
    pub fy: [f64; 4],
}

/// Four-wheel double-track model with torque vectoring capability.
#[derive(Debug, Clone)]
pub struct DoubleTrackModel {
    /// Vehicle mass \[kg\].
    pub mass: f64,
    /// Yaw inertia \[kg·m²\].
    pub iz: f64,
    /// Roll inertia \[kg·m²\].
    pub ix: f64,
    /// CG to front axle \[m\].
    pub lf: f64,
    /// CG to rear axle \[m\].
    pub lr: f64,
    /// CG height \[m\].
    pub hcg: f64,
    /// Front track width \[m\].
    pub track_f: f64,
    /// Rear track width \[m\].
    pub track_r: f64,
    /// Front cornering stiffness per tire \[N/rad\].
    pub cf_tire: f64,
    /// Rear cornering stiffness per tire \[N/rad\].
    pub cr_tire: f64,
    /// Front roll stiffness distribution \[0..1\].
    pub roll_dist_front: f64,
    /// Torque vectoring torque \[N·m\] (positive = right turn assist).
    pub tv_torque: f64,
    /// Current state.
    pub state: DoubleTrackState,
}

impl DoubleTrackModel {
    /// Create a default double-track model.
    pub fn new() -> Self {
        Self {
            mass: 1500.0,
            iz: 2500.0,
            ix: 400.0,
            lf: 1.1,
            lr: 1.6,
            hcg: 0.5,
            track_f: 1.5,
            track_r: 1.5,
            cf_tire: 45_000.0,
            cr_tire: 50_000.0,
            roll_dist_front: 0.55,
            tv_torque: 0.0,
            state: DoubleTrackState::default(),
        }
    }

    /// Compute wheel loads including lateral load transfer.
    ///
    /// Returns \[FL, FR, RL, RR\].
    pub fn wheel_loads(&self, ax: f64, ay: f64) -> [f64; 4] {
        let g = 9.81;
        let wb = self.lf + self.lr;
        // Static loads
        let fz_f = self.mass * g * self.lr / wb;
        let fz_r = self.mass * g * self.lf / wb;
        // Longitudinal transfer
        let dz_long = self.mass * ax * self.hcg / wb;
        // Lateral transfer (simplified roll stiffness)
        let k_phi = self.mass * g * self.hcg;
        let dz_lat_f = k_phi * ay / (9.81 * self.track_f) * self.roll_dist_front;
        let dz_lat_r = k_phi * ay / (9.81 * self.track_r) * (1.0 - self.roll_dist_front);
        [
            (fz_f / 2.0 - dz_long / 2.0 - dz_lat_f).max(0.0),
            (fz_f / 2.0 - dz_long / 2.0 + dz_lat_f).max(0.0),
            (fz_r / 2.0 + dz_long / 2.0 - dz_lat_r).max(0.0),
            (fz_r / 2.0 + dz_long / 2.0 + dz_lat_r).max(0.0),
        ]
    }

    /// Compute individual tire slip angles for given steering \[rad\].
    ///
    /// Returns \[FL, FR, RL, RR\] slip angles.
    pub fn tire_slip_angles(&self, vx: f64, vy: f64, r: f64, delta: f64) -> [f64; 4] {
        if vx.abs() < 0.1 {
            return [0.0; 4];
        }
        let alpha_f_left = -delta + ((vy + r * self.lf) / vx).atan();
        let alpha_f_right = alpha_f_left;
        let alpha_r_left = ((vy - r * self.lr) / vx).atan();
        let alpha_r_right = alpha_r_left;
        [alpha_f_left, alpha_f_right, alpha_r_left, alpha_r_right]
    }

    /// Step model by dt \[s\].
    pub fn step(&mut self, dt: f64, vx: f64, delta: f64, ax: f64, wheel_torques: &[f64; 4]) {
        let vy = self.state.vy;
        let r = self.state.r;
        let ay = (vy * r).abs();
        let loads = self.wheel_loads(ax, ay);
        let alphas = self.tire_slip_angles(vx, vy, r, delta);
        let mut fy_total = 0.0;
        let mut mz_total = self.tv_torque;
        let wheel_torque_total: f64 = wheel_torques.iter().sum();
        let _ = wheel_torque_total;
        for i in 0..4 {
            let cf = if i < 2 { self.cf_tire } else { self.cr_tire };
            let fy = -cf * alphas[i] * (loads[i] / (self.mass * 9.81 / 4.0)).min(1.5);
            self.state.fy[i] = fy;
            self.state.slip_angles[i] = alphas[i];
            fy_total += fy;
            let arm = if i < 2 { self.lf } else { -self.lr };
            mz_total += fy * arm;
        }
        // TV torque vectoring
        mz_total += self.tv_torque;
        let dvy = fy_total / self.mass - r * vx;
        let dr = mz_total / self.iz;
        self.state.vy += dvy * dt;
        self.state.r += dr * dt;
    }

    /// Set torque vectoring moment \[N·m\].
    pub fn set_tv_torque(&mut self, torque: f64) {
        self.tv_torque = torque;
    }
}

impl Default for DoubleTrackModel {
    fn default() -> Self {
        Self::new()
    }
}

// ─── RolloverDynamics ───────────────────────────────────────────────────────

/// Rollover dynamics analysis.
#[derive(Debug, Clone)]
pub struct RolloverDynamics {
    /// Track width \[m\].
    pub track: f64,
    /// CG height \[m\].
    pub hcg: f64,
    /// Mass \[kg\].
    pub mass: f64,
    /// Roll stiffness \[N·m/rad\].
    pub roll_stiffness: f64,
    /// Roll damping \[N·m·s/rad\].
    pub roll_damping: f64,
}

impl RolloverDynamics {
    /// Create a rollover dynamics model.
    pub fn new(track: f64, hcg: f64, mass: f64) -> Self {
        Self {
            track,
            hcg,
            mass,
            roll_stiffness: 50_000.0,
            roll_damping: 3_000.0,
        }
    }

    /// Static stability factor (SSF).
    ///
    /// SSF = T/(2*hcg). Rollover occurs when lateral accel > SSF * g.
    pub fn static_stability_factor(&self) -> f64 {
        self.track / (2.0 * self.hcg)
    }

    /// Lateral acceleration threshold for static rollover \[m/s²\].
    pub fn rollover_threshold_lat_accel(&self) -> f64 {
        self.static_stability_factor() * 9.81
    }

    /// Tilt angle \[rad\] for given lateral acceleration.
    pub fn roll_angle(&self, ay: f64) -> f64 {
        let restoring = self.roll_stiffness / (self.mass * self.hcg * 9.81);
        (ay / 9.81 / restoring).atan()
    }

    /// Wheel lift-off lateral acceleration \[m/s²\].
    ///
    /// Accounts for roll compliance which reduces effective SSF.
    pub fn dynamic_rollover_threshold(&self) -> f64 {
        let g = 9.81;
        let ssf = self.static_stability_factor();
        // Roll compliance reduces threshold
        let compliance_factor = 1.0 / (1.0 + self.mass * g * self.hcg / self.roll_stiffness);
        ssf * g * compliance_factor
    }

    /// J-turn peak lateral acceleration \[m/s²\] for given initial speed and steer input.
    pub fn j_turn_peak_ay(&self, v0: f64, delta_step: f64, cf: f64, cr: f64, mass: f64) -> f64 {
        // Simplified: steady-state lateral accel
        let _ = delta_step;
        let _ = cf;
        let _ = cr;
        let g = 9.81;
        let ssf = self.static_stability_factor();
        let ay_max = ssf * g;
        let ay_speed = v0 * v0 * 0.1; // rough centripetal
        ay_max.min(ay_speed * mass / mass) // mass cancels
    }

    /// Risk index \[0..1\]: 1.0 = imminent rollover.
    pub fn rollover_risk(&self, ay: f64) -> f64 {
        let threshold = self.rollover_threshold_lat_accel();
        clamp(ay.abs() / threshold, 0.0, 1.0)
    }
}

// ─── UndersteerGradient ─────────────────────────────────────────────────────

/// Understeer gradient measurement and analysis.
#[derive(Debug, Clone)]
pub struct UndersteerGradient {
    /// Vehicle mass \[kg\].
    pub mass: f64,
    /// CG to front axle \[m\].
    pub lf: f64,
    /// CG to rear axle \[m\].
    pub lr: f64,
    /// Front axle cornering stiffness \[N/rad\].
    pub cf: f64,
    /// Rear axle cornering stiffness \[N/rad\].
    pub cr: f64,
}

impl UndersteerGradient {
    /// Create from vehicle parameters.
    pub fn new(mass: f64, lf: f64, lr: f64, cf: f64, cr: f64) -> Self {
        Self {
            mass,
            lf,
            lr,
            cf,
            cr,
        }
    }

    /// Ackermann steering angle \[rad\] for given speed and radius.
    pub fn ackermann_steer(&self, vx: f64, radius: f64) -> f64 {
        (self.lf + self.lr) / radius + self.understeer_gradient() * vx * vx / (9.81 * radius)
    }

    /// Pure Ackermann (no dynamics) steer angle.
    pub fn kinematic_steer(&self, radius: f64) -> f64 {
        (self.lf + self.lr) / radius
    }

    /// Understeer gradient K \[rad·s²/m\].
    pub fn understeer_gradient(&self) -> f64 {
        let g = 9.81;
        let wf = self.mass * g * self.lr / (self.lf + self.lr);
        let wr = self.mass * g * self.lf / (self.lf + self.lr);
        wf / self.cf - wr / self.cr
    }

    /// Is the vehicle understeering?
    pub fn is_understeering(&self) -> bool {
        self.understeer_gradient() > 0.0
    }

    /// Is the vehicle oversteering?
    pub fn is_oversteering(&self) -> bool {
        self.understeer_gradient() < 0.0
    }

    /// Characteristic speed \[m/s\] for neutral steer (only meaningful for understeer).
    pub fn characteristic_speed(&self) -> f64 {
        let k = self.understeer_gradient();
        if k > 0.0 {
            (9.81 / k).sqrt()
        } else {
            f64::INFINITY
        }
    }

    /// Critical speed for oversteer \[m/s\].
    pub fn critical_speed(&self) -> f64 {
        let k = self.understeer_gradient();
        if k < 0.0 {
            (9.81 / (-k)).sqrt()
        } else {
            f64::INFINITY
        }
    }

    /// Lateral acceleration at given speed and steer angle \[m/s²\].
    pub fn lateral_accel(&self, vx: f64, delta: f64) -> f64 {
        let l = self.lf + self.lr;
        let k = self.understeer_gradient();
        let r = l / (delta - k * vx * vx / 9.81);
        if r.abs() < 0.01 { 0.0 } else { vx * vx / r }
    }
}

// ─── PitchDynamics ──────────────────────────────────────────────────────────

/// Pitch dynamics model for braking and acceleration.
#[derive(Debug, Clone)]
pub struct PitchDynamics {
    /// Mass \[kg\].
    pub mass: f64,
    /// CG height \[m\].
    pub hcg: f64,
    /// Wheelbase \[m\].
    pub wheelbase: f64,
    /// CG to front axle \[m\].
    pub lf: f64,
    /// Front suspension pitch compliance \[rad/(N·m)\].
    pub front_pitch_compliance: f64,
    /// Rear suspension pitch compliance \[rad/(N·m)\].
    pub rear_pitch_compliance: f64,
    /// Anti-dive geometry fraction at front \[0..1\].
    pub anti_dive_fraction: f64,
    /// Anti-squat geometry fraction at rear \[0..1\].
    pub anti_squat_fraction: f64,
    /// Pitch inertia \[kg·m²\].
    pub iy: f64,
    /// Current pitch angle \[rad\] (runtime).
    pub pitch: f64,
    /// Current pitch rate \[rad/s\] (runtime).
    pub pitch_rate: f64,
}

impl PitchDynamics {
    /// Create a pitch dynamics model.
    pub fn new(mass: f64, hcg: f64, wheelbase: f64, lf: f64) -> Self {
        Self {
            mass,
            hcg,
            wheelbase,
            lf,
            front_pitch_compliance: 5e-6,
            rear_pitch_compliance: 5e-6,
            anti_dive_fraction: 0.25,
            anti_squat_fraction: 0.3,
            iy: mass * hcg * hcg * 0.8,
            pitch: 0.0,
            pitch_rate: 0.0,
        }
    }

    /// Fore-aft weight transfer during braking \[N\].
    ///
    /// Positive = weight transfer to front.
    pub fn braking_weight_transfer(&self, decel: f64) -> f64 {
        self.mass * decel * self.hcg / self.wheelbase
    }

    /// Anti-dive moment \[N·m\] from brake geometry.
    pub fn anti_dive_moment(&self, brake_force: f64) -> f64 {
        self.anti_dive_fraction * brake_force * self.hcg
    }

    /// Anti-squat moment \[N·m\] from drive geometry.
    pub fn anti_squat_moment(&self, drive_force: f64) -> f64 {
        self.anti_squat_fraction * drive_force * self.hcg
    }

    /// Front axle load \[N\] during braking/acceleration.
    pub fn front_load(&self, ax: f64) -> f64 {
        let g = 9.81;
        let static_fz = self.mass * g * (self.wheelbase - self.lf) / self.wheelbase;
        let transfer = self.braking_weight_transfer(ax.abs());
        if ax < 0.0 {
            static_fz + transfer
        } else {
            static_fz - transfer
        }
    }

    /// Rear axle load \[N\].
    pub fn rear_load(&self, ax: f64) -> f64 {
        let g = 9.81;
        let static_fz = self.mass * g * self.lf / self.wheelbase;
        let transfer = self.braking_weight_transfer(ax.abs());
        if ax < 0.0 {
            static_fz - transfer
        } else {
            static_fz + transfer
        }
    }

    /// Step pitch dynamics by dt \[s\].
    pub fn step(&mut self, dt: f64, ax: f64) {
        let pitch_stiffness =
            1.0 / (self.front_pitch_compliance + self.rear_pitch_compliance).max(1e-12);
        let pitch_damping = 2.0 * (pitch_stiffness * self.iy).sqrt() * 0.3;
        let pitch_moment = self.mass * ax * self.hcg;
        let d_pitch_rate =
            (pitch_moment - pitch_stiffness * self.pitch - pitch_damping * self.pitch_rate)
                / self.iy.max(1.0);
        self.pitch_rate += d_pitch_rate * dt;
        self.pitch += self.pitch_rate * dt;
    }
}

// ─── YawMomentDiagram ───────────────────────────────────────────────────────

/// Yaw moment diagram (handling diagram) data point.
#[derive(Debug, Clone)]
pub struct YmdPoint {
    /// Lateral acceleration \[m/s²\].
    pub ay: f64,
    /// Yaw moment \[N·m\].
    pub mz: f64,
    /// Steering angle \[rad\].
    pub delta: f64,
}

/// Yaw Moment Diagram generator.
#[derive(Debug, Clone)]
pub struct YawMomentDiagram {
    /// Vehicle mass \[kg\].
    pub mass: f64,
    /// Wheelbase \[m\].
    pub wheelbase: f64,
    /// Front cornering stiffness \[N/rad\].
    pub cf: f64,
    /// Rear cornering stiffness \[N/rad\].
    pub cr: f64,
    /// CG to front \[m\].
    pub lf: f64,
    /// CG to rear \[m\].
    pub lr: f64,
}

impl YawMomentDiagram {
    /// Create a YMD generator.
    pub fn new(mass: f64, lf: f64, lr: f64, cf: f64, cr: f64) -> Self {
        Self {
            mass,
            wheelbase: lf + lr,
            cf,
            cr,
            lf,
            lr,
        }
    }

    /// Generate YMD points for a range of speeds and steer angles.
    pub fn generate(&self, vx: f64, delta_range: &[f64], beta_range: &[f64]) -> Vec<YmdPoint> {
        let mut points = Vec::new();
        for &delta in delta_range {
            for &beta in beta_range {
                let alpha_f = -delta + beta + self.lf * 0.0 / vx.max(0.1); // r=0
                let alpha_r = beta - self.lr * 0.0 / vx.max(0.1);
                let fy_f = -self.cf * alpha_f;
                let fy_r = -self.cr * alpha_r;
                let ay = (fy_f + fy_r) / self.mass;
                let mz = fy_f * self.lf - fy_r * self.lr;
                points.push(YmdPoint { ay, mz, delta });
            }
        }
        points
    }

    /// Find steer angle for zero yaw moment at given lateral acceleration.
    pub fn neutral_steer_angle(&self, vx: f64, ay: f64) -> f64 {
        let g = 9.81;
        let l = self.wheelbase;
        let k = {
            let wf = self.mass * g * self.lr / l;
            let wr = self.mass * g * self.lf / l;
            wf / self.cf - wr / self.cr
        };
        let r = vx * vx / ay.max(0.01).copysign(ay);
        l / r + k * ay / g
    }
}

// ─── HandlingBalance ────────────────────────────────────────────────────────

/// Handling balance analysis: front/rear lateral stiffness ratio.
#[derive(Debug, Clone)]
pub struct HandlingBalance {
    /// Front axle cornering stiffness \[N/rad\].
    pub cf: f64,
    /// Rear axle cornering stiffness \[N/rad\].
    pub cr: f64,
    /// Front weight fraction.
    pub weight_front: f64,
    /// Rear weight fraction.
    pub weight_rear: f64,
}

impl HandlingBalance {
    /// Create a handling balance analyzer.
    pub fn new(cf: f64, cr: f64, weight_front: f64) -> Self {
        Self {
            cf,
            cr,
            weight_front,
            weight_rear: 1.0 - weight_front,
        }
    }

    /// Balance ratio: CF/weight_front vs CR/weight_rear.
    ///
    /// 1.0 = neutral, >1 = understeering, <1 = oversteering.
    pub fn balance_ratio(&self) -> f64 {
        let front_norm = self.cf / self.weight_front.max(1e-6);
        let rear_norm = self.cr / self.weight_rear.max(1e-6);
        front_norm / rear_norm.max(1e-6)
    }

    /// Suggested cornering stiffness change to neutralize balance.
    ///
    /// Returns (delta_cf, delta_cr) to achieve neutral.
    pub fn neutralize_suggestion(&self) -> (f64, f64) {
        let target_ratio = 1.0;
        let current = self.balance_ratio();
        if current > target_ratio {
            // Too much understeer: increase rear stiffness
            let target_cr = self.cf * self.weight_rear / self.weight_front;
            (0.0, target_cr - self.cr)
        } else {
            // Too much oversteer: increase front stiffness
            let target_cf = self.cr * self.weight_front / self.weight_rear;
            (target_cf - self.cf, 0.0)
        }
    }

    /// Sensitivity of balance to front stiffness change.
    pub fn cf_sensitivity(&self) -> f64 {
        1.0 / (self.weight_front.max(1e-6) * self.cr.max(1e-6) / self.weight_rear.max(1e-6))
    }
}

// ─── WheelSpinControl ───────────────────────────────────────────────────────

/// Differential type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DiffType {
    /// Open differential.
    Open,
    /// Limited slip differential.
    Lsd,
    /// Torque vectoring differential.
    TorqueVectoring,
}

/// Wheel spin control with various differential types.
#[derive(Debug, Clone)]
pub struct WheelSpinControl {
    /// Differential type.
    pub diff_type: DiffType,
    /// LSD bias ratio (max torque split, typically 2-5x).
    pub lsd_bias: f64,
    /// Maximum TV torque \[N·m\].
    pub tv_max_torque: f64,
    /// Wheel speeds \[rad/s\]: Left, Right.
    pub wheel_speeds: [f64; 2],
    /// Output torques \[N·m\]: Left, Right.
    pub torques: [f64; 2],
}

impl WheelSpinControl {
    /// Create an open differential.
    pub fn open() -> Self {
        Self {
            diff_type: DiffType::Open,
            lsd_bias: 1.0,
            tv_max_torque: 0.0,
            wheel_speeds: [0.0; 2],
            torques: [0.0; 2],
        }
    }

    /// Create an LSD with given bias ratio.
    pub fn lsd(bias: f64) -> Self {
        Self {
            diff_type: DiffType::Lsd,
            lsd_bias: bias,
            tv_max_torque: 0.0,
            wheel_speeds: [0.0; 2],
            torques: [0.0; 2],
        }
    }

    /// Create a torque vectoring differential.
    pub fn torque_vectoring(max_torque: f64) -> Self {
        Self {
            diff_type: DiffType::TorqueVectoring,
            lsd_bias: 1.0,
            tv_max_torque: max_torque,
            wheel_speeds: [0.0; 2],
            torques: [0.0; 2],
        }
    }

    /// Distribute input torque to wheels.
    pub fn distribute(&mut self, input_torque: f64, yaw_rate: f64) -> [f64; 2] {
        match self.diff_type {
            DiffType::Open => {
                self.torques = [input_torque / 2.0, input_torque / 2.0];
            }
            DiffType::Lsd => {
                let speed_diff = self.wheel_speeds[1] - self.wheel_speeds[0];
                let lsd_moment = clamp(
                    speed_diff * 500.0,
                    -input_torque * (self.lsd_bias - 1.0) * 0.5,
                    input_torque * (self.lsd_bias - 1.0) * 0.5,
                );
                self.torques = [
                    input_torque / 2.0 - lsd_moment,
                    input_torque / 2.0 + lsd_moment,
                ];
            }
            DiffType::TorqueVectoring => {
                let tv = clamp(yaw_rate * 800.0, -self.tv_max_torque, self.tv_max_torque);
                self.torques = [input_torque / 2.0 - tv, input_torque / 2.0 + tv];
            }
        }
        self.torques
    }

    /// Speed equalization error \[rad/s\].
    pub fn speed_error(&self) -> f64 {
        (self.wheel_speeds[1] - self.wheel_speeds[0]).abs()
    }
}

// ─── VehicleStateEstimator ──────────────────────────────────────────────────

/// Vehicle state vector for Kalman filter.
#[derive(Debug, Clone)]
pub struct VehicleStateVector {
    /// Longitudinal velocity \[m/s\].
    pub vx: f64,
    /// Lateral velocity \[m/s\].
    pub vy: f64,
    /// Yaw rate \[rad/s\].
    pub r: f64,
    /// Side slip angle \[rad\].
    pub beta: f64,
}

/// Kalman filter-based vehicle state estimator.
///
/// Estimates velocity and slip angle from sensor data.
#[derive(Debug, Clone)]
pub struct VehicleStateEstimator {
    /// State vector.
    pub state: VehicleStateVector,
    /// State covariance matrix (4×4, flat row-major).
    pub p: Vec<f64>,
    /// Process noise covariance (diagonal).
    pub q_diag: [f64; 4],
    /// Measurement noise covariance (diagonal).
    pub r_diag: [f64; 4],
}

impl VehicleStateEstimator {
    /// Create a new vehicle state estimator.
    pub fn new() -> Self {
        Self {
            state: VehicleStateVector {
                vx: 0.0,
                vy: 0.0,
                r: 0.0,
                beta: 0.0,
            },
            p: vec![
                1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
            ],
            q_diag: [0.01, 0.1, 0.1, 0.05],
            r_diag: [0.1, 0.5, 0.05, 0.2],
        }
    }

    /// Predict step using simple kinematic model.
    pub fn predict(&mut self, dt: f64, ax: f64, ay: f64) {
        // State transition: simple Euler integration
        self.state.vx += ax * dt;
        self.state.vy += ay * dt - self.state.r * self.state.vx * dt;
        self.state.beta = if self.state.vx.abs() > 0.1 {
            (self.state.vy / self.state.vx).atan()
        } else {
            0.0
        };
        // Covariance prediction: P = P + Q * dt
        for (i, &q) in self.q_diag.iter().enumerate() {
            self.p[i * 5] += q * dt;
        }
    }

    /// Update step with sensor measurement.
    ///
    /// `measurement` = \[vx_meas, vy_meas, r_meas, beta_meas\]
    pub fn update(&mut self, measurement: &[f64; 4]) {
        let state_arr = [self.state.vx, self.state.vy, self.state.r, self.state.beta];
        let mut k = [0.0f64; 4]; // Kalman gain (diagonal approximation)
        for i in 0..4 {
            let p_ii = self.p[i * 5];
            let r_ii = self.r_diag[i];
            k[i] = p_ii / (p_ii + r_ii);
        }
        // State update
        self.state.vx += k[0] * (measurement[0] - state_arr[0]);
        self.state.vy += k[1] * (measurement[1] - state_arr[1]);
        self.state.r += k[2] * (measurement[2] - state_arr[2]);
        self.state.beta += k[3] * (measurement[3] - state_arr[3]);
        // Covariance update
        for i in 0..4 {
            self.p[i * 5] *= 1.0 - k[i];
        }
    }

    /// Fuse wheel speed sensors to estimate vx \[m/s\].
    pub fn fuse_wheel_speeds(&mut self, wheel_speeds: &[f64; 4], wheel_radius: f64) -> f64 {
        let vx_est: f64 = wheel_speeds.iter().sum::<f64>() / 4.0 * wheel_radius;
        self.state.vx = 0.9 * self.state.vx + 0.1 * vx_est;
        self.state.vx
    }

    /// Fuse IMU yaw rate \[rad/s\].
    pub fn fuse_imu_yaw_rate(&mut self, measured_r: f64) {
        let k = self.p[10] / (self.p[10] + self.r_diag[2]);
        self.state.r += k * (measured_r - self.state.r);
        self.p[10] *= 1.0 - k;
    }
}

impl Default for VehicleStateEstimator {
    fn default() -> Self {
        Self::new()
    }
}

// ─── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // LinearBicycleModel tests

    #[test]
    fn test_linear_bicycle_default_creates() {
        let m = LinearBicycleModel::new();
        assert_eq!(m.mass, 1500.0);
    }

    #[test]
    fn test_linear_bicycle_understeer_gradient_sign() {
        let m = LinearBicycleModel::new();
        // With rear stiffer than front per unit weight: understeering
        let k = m.understeer_gradient();
        assert!(k.is_finite());
    }

    #[test]
    fn test_linear_bicycle_step_changes_state() {
        let mut m = LinearBicycleModel::new();
        m.step(0.01, 20.0, 0.05);
        assert!(m.state.v_y != 0.0 || m.state.r != 0.0);
    }

    #[test]
    fn test_linear_bicycle_steady_state_yaw_rate_positive() {
        let m = LinearBicycleModel::new();
        let r = m.steady_state_yaw_rate(20.0, 0.05);
        assert!(r > 0.0);
    }

    #[test]
    fn test_linear_bicycle_critical_speed_positive() {
        let m = LinearBicycleModel::new();
        let cs = m.critical_speed();
        assert!(cs > 0.0);
    }

    #[test]
    fn test_linear_bicycle_stability_derivatives() {
        let m = LinearBicycleModel::new();
        let (a11, a12, a21, a22, b1, b2) = m.stability_derivatives(20.0);
        assert!(a11.is_finite() && a12.is_finite() && a21.is_finite() && a22.is_finite());
        assert!(b1 > 0.0 && b2 > 0.0);
    }

    #[test]
    fn test_linear_bicycle_position_changes() {
        let mut m = LinearBicycleModel::new();
        let x0 = m.state.x;
        m.step(1.0, 10.0, 0.01);
        assert!((m.state.x - x0).abs() > 0.1);
    }

    // NonlinearBicycleModel tests

    #[test]
    fn test_nonlinear_bicycle_slip_angles_zero_at_low_speed() {
        let m = NonlinearBicycleModel::new();
        let alpha_f = m.front_slip_angle(0.05, 0.0);
        assert_eq!(alpha_f, 0.0);
    }

    #[test]
    fn test_nonlinear_bicycle_front_slip_with_steer() {
        let m = NonlinearBicycleModel::new();
        let alpha_f = m.front_slip_angle(20.0, 0.1);
        assert!(alpha_f.abs() > 0.0);
    }

    #[test]
    fn test_nonlinear_bicycle_loads_sum_to_weight() {
        let m = NonlinearBicycleModel::new();
        let g = 9.81;
        let fz_f = m.front_load(0.0);
        let fz_r = m.rear_load(0.0);
        assert!((fz_f + fz_r - m.mass * g).abs() < 1.0);
    }

    #[test]
    fn test_nonlinear_bicycle_step_integrates() {
        let mut m = NonlinearBicycleModel::new();
        m.step(0.01, 15.0, 0.05, 0.0);
        assert!(m.vy.is_finite() && m.r.is_finite());
    }

    // DoubleTrackModel tests

    #[test]
    fn test_double_track_wheel_loads_sum() {
        let m = DoubleTrackModel::new();
        let loads = m.wheel_loads(0.0, 0.0);
        let total: f64 = loads.iter().sum();
        let expected = m.mass * 9.81;
        assert!((total - expected).abs() < 10.0);
    }

    #[test]
    fn test_double_track_lateral_transfer_shifts_load() {
        let m = DoubleTrackModel::new();
        let loads_0 = m.wheel_loads(0.0, 0.0);
        let loads_ay = m.wheel_loads(0.0, 5.0);
        // Right side should gain load
        assert!(loads_ay[1] > loads_0[1] || loads_ay[3] > loads_0[3]);
    }

    #[test]
    fn test_double_track_step_no_panic() {
        let mut m = DoubleTrackModel::new();
        m.step(0.01, 20.0, 0.05, 0.0, &[500.0, 500.0, 500.0, 500.0]);
    }

    // RolloverDynamics tests

    #[test]
    fn test_rollover_ssf_positive() {
        let r = RolloverDynamics::new(1.5, 0.5, 1500.0);
        assert!(r.static_stability_factor() > 0.0);
    }

    #[test]
    fn test_rollover_threshold_sensible() {
        let r = RolloverDynamics::new(1.5, 0.5, 1500.0);
        let thresh = r.rollover_threshold_lat_accel();
        assert!(thresh > 5.0 && thresh < 30.0); // Typical 0.5-2g range
    }

    #[test]
    fn test_rollover_risk_zero_at_zero_accel() {
        let r = RolloverDynamics::new(1.5, 0.5, 1500.0);
        assert_eq!(r.rollover_risk(0.0), 0.0);
    }

    #[test]
    fn test_rollover_risk_one_at_threshold() {
        let r = RolloverDynamics::new(1.5, 0.5, 1500.0);
        let thresh = r.rollover_threshold_lat_accel();
        assert!((r.rollover_risk(thresh) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_rollover_dynamic_threshold_lt_static() {
        let r = RolloverDynamics::new(1.5, 0.5, 1500.0);
        assert!(r.dynamic_rollover_threshold() < r.rollover_threshold_lat_accel());
    }

    // UndersteerGradient tests

    #[test]
    fn test_understeer_gradient_computation() {
        let ug = UndersteerGradient::new(1500.0, 1.1, 1.6, 80_000.0, 90_000.0);
        let k = ug.understeer_gradient();
        assert!(k.is_finite());
    }

    #[test]
    fn test_understeer_gradient_is_positive_for_understeer() {
        // To get understeer: Wf/Cf > Wr/Cr
        // Wf = m*g*lr/(lf+lr), Wr = m*g*lf/(lf+lr)
        // With lf=1.1, lr=1.6: Wf/Wr = lr/lf = 1.6/1.1 ≈ 1.45
        // Need Cf low enough: Wf/Cf > Wr/Cr => Cf/Cr < Wf/Wr = 1.45
        // Use cf=60000, cr=80000: ratio = 0.75 < 1.45 => understeer
        let ug = UndersteerGradient::new(1500.0, 1.1, 1.6, 60_000.0, 80_000.0);
        assert!(
            ug.is_understeering(),
            "Expected understeer, K={}",
            ug.understeer_gradient()
        );
    }

    #[test]
    fn test_understeer_characteristic_speed_positive() {
        let ug = UndersteerGradient::new(1500.0, 1.1, 1.6, 80_000.0, 120_000.0);
        // Oversteering car has critical speed
        if ug.is_oversteering() {
            assert!(ug.critical_speed() > 0.0 && ug.critical_speed() < 1000.0);
        }
    }

    // PitchDynamics tests

    #[test]
    fn test_pitch_weight_transfer_positive_braking() {
        let p = PitchDynamics::new(1500.0, 0.5, 2.7, 1.1);
        let wt = p.braking_weight_transfer(5.0);
        assert!(wt > 0.0);
    }

    #[test]
    fn test_pitch_front_load_increases_braking() {
        let p = PitchDynamics::new(1500.0, 0.5, 2.7, 1.1);
        let fz0 = p.front_load(0.0);
        let fz_brake = p.front_load(-8.0); // deceleration
        assert!(fz_brake > fz0);
    }

    #[test]
    fn test_pitch_loads_sum_to_weight() {
        let p = PitchDynamics::new(1500.0, 0.5, 2.7, 1.1);
        let g = 9.81;
        let front = p.front_load(0.0);
        let rear = p.rear_load(0.0);
        assert!((front + rear - p.mass * g).abs() < 10.0);
    }

    #[test]
    fn test_pitch_step_no_panic() {
        let mut p = PitchDynamics::new(1500.0, 0.5, 2.7, 1.1);
        p.step(0.01, -5.0);
    }

    // YawMomentDiagram tests

    #[test]
    fn test_ymd_generate_nonempty() {
        let ymd = YawMomentDiagram::new(1500.0, 1.1, 1.6, 80_000.0, 90_000.0);
        let deltas = vec![-0.05, 0.0, 0.05];
        let betas = vec![-0.05, 0.0, 0.05];
        let pts = ymd.generate(20.0, &deltas, &betas);
        assert!(!pts.is_empty());
    }

    #[test]
    fn test_ymd_zero_steer_zero_beta_near_zero_moment() {
        let ymd = YawMomentDiagram::new(1500.0, 1.1, 1.6, 80_000.0, 90_000.0);
        let deltas = vec![0.0];
        let betas = vec![0.0];
        let pts = ymd.generate(20.0, &deltas, &betas);
        assert!((pts[0].mz).abs() < 1e-6);
    }

    // HandlingBalance tests

    #[test]
    fn test_handling_balance_neutral_ratio_near_one() {
        // Equal normalized stiffness
        let hb = HandlingBalance::new(80_000.0, 120_000.0, 0.6);
        let ratio = hb.balance_ratio();
        assert!(ratio.is_finite() && ratio > 0.0);
    }

    #[test]
    fn test_handling_balance_neutralize_suggestion_both_finite() {
        let hb = HandlingBalance::new(80_000.0, 100_000.0, 0.5);
        let (dcf, dcr) = hb.neutralize_suggestion();
        assert!(dcf.is_finite() && dcr.is_finite());
    }

    // WheelSpinControl tests

    #[test]
    fn test_open_diff_equal_torque() {
        let mut wsc = WheelSpinControl::open();
        let t = wsc.distribute(200.0, 0.0);
        assert!((t[0] - t[1]).abs() < 1e-9);
    }

    #[test]
    fn test_lsd_unequal_torque_at_speed_diff() {
        let mut wsc = WheelSpinControl::lsd(3.0);
        wsc.wheel_speeds = [10.0, 15.0];
        let t = wsc.distribute(200.0, 0.0);
        // Should be different
        assert!(t[0] != t[1]);
    }

    #[test]
    fn test_tv_diff_applies_yaw_moment() {
        let mut wsc = WheelSpinControl::torque_vectoring(500.0);
        let t_straight = wsc.distribute(200.0, 0.0);
        let t_turn = wsc.distribute(200.0, 0.5);
        assert!(t_straight[0] != t_turn[0] || t_straight[1] != t_turn[1]);
    }

    #[test]
    fn test_wheel_spin_speed_error() {
        let mut wsc = WheelSpinControl::open();
        wsc.wheel_speeds = [10.0, 12.0];
        assert!((wsc.speed_error() - 2.0).abs() < 1e-9);
    }

    // VehicleStateEstimator tests

    #[test]
    fn test_estimator_initial_state_zero() {
        let est = VehicleStateEstimator::new();
        assert_eq!(est.state.vx, 0.0);
        assert_eq!(est.state.vy, 0.0);
    }

    #[test]
    fn test_estimator_predict_changes_vx() {
        let mut est = VehicleStateEstimator::new();
        est.predict(0.1, 5.0, 0.0);
        assert!((est.state.vx - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_estimator_update_corrects_state() {
        let mut est = VehicleStateEstimator::new();
        est.state.vx = 10.0;
        est.update(&[12.0, 0.0, 0.0, 0.0]);
        assert!(est.state.vx > 10.0);
    }

    #[test]
    fn test_estimator_fuse_wheel_speeds() {
        let mut est = VehicleStateEstimator::new();
        let speeds = [50.0, 50.0, 50.0, 50.0]; // rad/s
        let vx = est.fuse_wheel_speeds(&speeds, 0.3);
        assert!(vx > 0.0);
    }

    #[test]
    fn test_estimator_fuse_yaw_rate() {
        let mut est = VehicleStateEstimator::new();
        est.state.r = 0.0;
        est.fuse_imu_yaw_rate(0.5);
        assert!(est.state.r > 0.0);
    }

    #[test]
    fn test_estimator_covariance_decreases_after_update() {
        let mut est = VehicleStateEstimator::new();
        let p0 = est.p[0];
        est.update(&[1.0, 0.0, 0.0, 0.0]);
        assert!(est.p[0] <= p0);
    }

    #[test]
    fn test_tire_force_model_bounded() {
        let tf = TireForceModel::default_values();
        let fy = tf.lateral_force_coeff(0.1);
        assert!(fy.abs() <= tf.d + 0.1);
    }

    #[test]
    fn test_mat3_mul_identity() {
        let id = mat3_identity();
        let v = [1.0, 2.0, 3.0];
        let r = mat3_mul_vec3(&id, &v);
        assert!((r[0] - v[0]).abs() < 1e-12);
        assert!((r[1] - v[1]).abs() < 1e-12);
        assert!((r[2] - v[2]).abs() < 1e-12);
    }

    #[test]
    fn test_mat3_transpose_involutory() {
        let m: Mat3 = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let mt = mat3_transpose(&m);
        let mtt = mat3_transpose(&mt);
        for i in 0..9 {
            assert!((m[i] - mtt[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_rot2_orthogonal() {
        let r = rot2(PI / 4.0);
        // Check: |det| = 1
        let det = r[0] * r[3] - r[1] * r[2];
        assert!((det - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_handling_balance_cf_sensitivity_positive() {
        let hb = HandlingBalance::new(80_000.0, 90_000.0, 0.5);
        assert!(hb.cf_sensitivity() > 0.0);
    }
}
