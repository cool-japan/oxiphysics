#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Gyroscope and spinning body dynamics.
//!
//! This module covers Euler's equations of rigid body rotation, precession and
//! nutation computation, gyroscopic torque, spinning top (Lagrange top)
//! simulation, control moment gyroscope (CMG) and reaction wheel models,
//! rate gyroscope sensor (with noise/drift/bias), gyrocompass finding of true
//! north, attitude determination, angular momentum conservation, torque-free
//! (Poinsot) motion, quaternion-based Euler angle integration, gyroscopic
//! stabilisation analysis, and gimbal lock detection.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Vector / quaternion helpers  (no nalgebra — plain [f64; 3])
// ---------------------------------------------------------------------------

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

#[inline]
fn normalize(a: [f64; 3]) -> [f64; 3] {
    let n = norm(a);
    if n < 1e-15 {
        return [0.0, 0.0, 0.0];
    }
    [a[0] / n, a[1] / n, a[2] / n]
}

#[inline]
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn neg3(a: [f64; 3]) -> [f64; 3] {
    [-a[0], -a[1], -a[2]]
}

// ---- Quaternion (scalar-last: [x, y, z, w]) ----

#[inline]
fn quat_identity() -> [f64; 4] {
    [0.0, 0.0, 0.0, 1.0]
}

#[inline]
fn quat_mul(p: [f64; 4], q: [f64; 4]) -> [f64; 4] {
    let [px, py, pz, pw] = p;
    let [qx, qy, qz, qw] = q;
    [
        pw * qx + px * qw + py * qz - pz * qy,
        pw * qy - px * qz + py * qw + pz * qx,
        pw * qz + px * qy - py * qx + pz * qw,
        pw * qw - px * qx - py * qy - pz * qz,
    ]
}

#[inline]
fn quat_conjugate(q: [f64; 4]) -> [f64; 4] {
    [-q[0], -q[1], -q[2], q[3]]
}

#[inline]
fn quat_normalize(q: [f64; 4]) -> [f64; 4] {
    let mag = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    if mag < 1e-15 {
        return quat_identity();
    }
    [q[0] / mag, q[1] / mag, q[2] / mag, q[3] / mag]
}

/// Rotate a vector by a unit quaternion.
fn quat_rotate(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    let vq = [v[0], v[1], v[2], 0.0];
    let r = quat_mul(quat_mul(q, vq), quat_conjugate(q));
    [r[0], r[1], r[2]]
}

/// First-order quaternion integration: q' = q + 0.5 * dt * omega_quat * q.
fn quat_integrate(q: [f64; 4], omega: [f64; 3], dt: f64) -> [f64; 4] {
    let half_dt = 0.5 * dt;
    let dq = [
        (omega[0] * q[3] + omega[2] * q[1] - omega[1] * q[2]) * half_dt,
        (omega[1] * q[3] - omega[2] * q[0] + omega[0] * q[2]) * half_dt,
        (omega[2] * q[3] + omega[1] * q[0] - omega[0] * q[1]) * half_dt,
        (-omega[0] * q[0] - omega[1] * q[1] - omega[2] * q[2]) * half_dt,
    ];
    quat_normalize([q[0] + dq[0], q[1] + dq[1], q[2] + dq[2], q[3] + dq[3]])
}

/// Convert quaternion to Euler angles (ZYX convention) `[roll, pitch, yaw]` in radians.
fn quat_to_euler(q: [f64; 4]) -> [f64; 3] {
    let [x, y, z, w] = q;
    let sinr = 2.0 * (w * x + y * z);
    let cosr = 1.0 - 2.0 * (x * x + y * y);
    let roll = sinr.atan2(cosr);

    let sinp = 2.0 * (w * y - z * x);
    let pitch = if sinp.abs() >= 1.0 {
        (PI / 2.0).copysign(sinp)
    } else {
        sinp.asin()
    };

    let siny = 2.0 * (w * z + x * y);
    let cosy = 1.0 - 2.0 * (y * y + z * z);
    let yaw = siny.atan2(cosy);

    [roll, pitch, yaw]
}

/// Convert Euler angles (ZYX convention) to quaternion.
fn euler_to_quat(roll: f64, pitch: f64, yaw: f64) -> [f64; 4] {
    let (sr, cr) = (roll * 0.5).sin_cos();
    let (sp, cp) = (pitch * 0.5).sin_cos();
    let (sy, cy) = (yaw * 0.5).sin_cos();
    [
        sr * cp * cy - cr * sp * sy,
        cr * sp * cy + sr * cp * sy,
        cr * cp * sy - sr * sp * cy,
        cr * cp * cy + sr * sp * sy,
    ]
}

// ---------------------------------------------------------------------------
// 3x3 diagonal inertia tensor helpers
// ---------------------------------------------------------------------------

/// Principal moments of inertia `[Ixx, Iyy, Izz]`.
pub type Inertia3 = [f64; 3];

/// Compute `I * omega` for a diagonal inertia tensor.
#[inline]
fn inertia_mul(inertia: Inertia3, omega: [f64; 3]) -> [f64; 3] {
    [
        inertia[0] * omega[0],
        inertia[1] * omega[1],
        inertia[2] * omega[2],
    ]
}

/// Compute `I^{-1} * v` for a diagonal inertia tensor.
#[inline]
fn inertia_inv_mul(inertia: Inertia3, v: [f64; 3]) -> [f64; 3] {
    [
        if inertia[0].abs() > 1e-15 {
            v[0] / inertia[0]
        } else {
            0.0
        },
        if inertia[1].abs() > 1e-15 {
            v[1] / inertia[1]
        } else {
            0.0
        },
        if inertia[2].abs() > 1e-15 {
            v[2] / inertia[2]
        } else {
            0.0
        },
    ]
}

// ---------------------------------------------------------------------------
// Euler's equations of rigid body rotation
// ---------------------------------------------------------------------------

/// State of a spinning rigid body.
#[derive(Debug, Clone)]
pub struct SpinningBodyState {
    /// Angular velocity in body frame `[omega_x, omega_y, omega_z]` (rad/s).
    pub omega: [f64; 3],
    /// Orientation quaternion (scalar-last: `[x, y, z, w]`).
    pub orientation: [f64; 4],
    /// Principal moments of inertia `[Ixx, Iyy, Izz]`.
    pub inertia: Inertia3,
}

impl SpinningBodyState {
    /// Create a new spinning body state.
    pub fn new(inertia: Inertia3, omega: [f64; 3]) -> Self {
        Self {
            omega,
            orientation: quat_identity(),
            inertia,
        }
    }

    /// Angular momentum in body frame: `L = I * omega`.
    pub fn angular_momentum_body(&self) -> [f64; 3] {
        inertia_mul(self.inertia, self.omega)
    }

    /// Angular momentum in world frame.
    pub fn angular_momentum_world(&self) -> [f64; 3] {
        quat_rotate(self.orientation, self.angular_momentum_body())
    }

    /// Rotational kinetic energy: `0.5 * omega^T * I * omega`.
    pub fn kinetic_energy(&self) -> f64 {
        let l = self.angular_momentum_body();
        0.5 * dot(self.omega, l)
    }

    /// Spin rate (magnitude of angular velocity).
    pub fn spin_rate(&self) -> f64 {
        norm(self.omega)
    }
}

/// Compute the angular acceleration from Euler's equations:
///
/// ```text
/// I * alpha = torque - omega x (I * omega)
/// ```
///
/// Returns `alpha = I^{-1} * (torque - omega x (I * omega))`.
pub fn euler_equations(inertia: Inertia3, omega: [f64; 3], torque: [f64; 3]) -> [f64; 3] {
    let i_omega = inertia_mul(inertia, omega);
    let gyro = cross(omega, i_omega);
    let rhs = sub3(torque, gyro);
    inertia_inv_mul(inertia, rhs)
}

/// Integrate a spinning body state by one time step using RK4.
pub fn integrate_spinning_body_rk4(state: &mut SpinningBodyState, torque: [f64; 3], dt: f64) {
    let inertia = state.inertia;

    // RK4 for angular velocity
    let f = |omega: [f64; 3]| euler_equations(inertia, omega, torque);

    let k1 = f(state.omega);
    let omega2 = add3(state.omega, scale(k1, dt * 0.5));
    let k2 = f(omega2);
    let omega3 = add3(state.omega, scale(k2, dt * 0.5));
    let k3 = f(omega3);
    let omega4 = add3(state.omega, scale(k3, dt));
    let k4 = f(omega4);

    state.omega = [
        state.omega[0] + dt / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]),
        state.omega[1] + dt / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]),
        state.omega[2] + dt / 6.0 * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]),
    ];

    // Integrate orientation
    state.orientation = quat_integrate(state.orientation, state.omega, dt);
}

/// Integrate a spinning body using semi-implicit Euler (simpler, first order).
pub fn integrate_spinning_body_euler(state: &mut SpinningBodyState, torque: [f64; 3], dt: f64) {
    let alpha = euler_equations(state.inertia, state.omega, torque);
    state.omega = add3(state.omega, scale(alpha, dt));
    state.orientation = quat_integrate(state.orientation, state.omega, dt);
}

// ---------------------------------------------------------------------------
// Precession rate computation
// ---------------------------------------------------------------------------

/// Compute the steady-state precession rate of a symmetric gyroscope.
///
/// For a symmetric top with spin angular momentum `L = I_spin * omega_spin`
/// under gravitational torque `tau = m * g * r`, the precession rate is:
///
/// ```text
/// Omega_p = tau / L = m * g * r / (I_spin * omega_spin)
/// ```
pub fn precession_rate(
    spin_inertia: f64,
    spin_rate: f64,
    mass: f64,
    gravity: f64,
    distance_to_pivot: f64,
) -> f64 {
    let l_spin = spin_inertia * spin_rate;
    if l_spin.abs() < 1e-15 {
        return 0.0;
    }
    mass * gravity * distance_to_pivot / l_spin
}

/// Compute the precession angular velocity vector given spin axis, torque
/// axis, and precession rate magnitude.
pub fn precession_omega(precession_rate: f64, precession_axis: [f64; 3]) -> [f64; 3] {
    scale(normalize(precession_axis), precession_rate)
}

// ---------------------------------------------------------------------------
// Nutation angle and frequency
// ---------------------------------------------------------------------------

/// Compute the nutation frequency of a symmetric gyroscope.
///
/// For an axisymmetric body with principal moments `I_a` (spin axis) and
/// `I_t` (transverse), spinning at rate `omega_spin`:
///
/// ```text
/// omega_n = (I_a - I_t) * omega_spin / I_t
/// ```
pub fn nutation_frequency(spin_inertia: f64, transverse_inertia: f64, spin_rate: f64) -> f64 {
    if transverse_inertia.abs() < 1e-15 {
        return 0.0;
    }
    (spin_inertia - transverse_inertia) * spin_rate / transverse_inertia
}

/// Compute the nutation angle from the ratio of transverse to spin angular
/// momentum. Returns radians.
pub fn nutation_angle(
    omega_transverse: f64,
    omega_spin: f64,
    spin_inertia: f64,
    transverse_inertia: f64,
) -> f64 {
    let l_transverse = transverse_inertia * omega_transverse;
    let l_spin = spin_inertia * omega_spin;
    if l_spin.abs() < 1e-15 {
        return 0.0;
    }
    (l_transverse / l_spin).atan()
}

// ---------------------------------------------------------------------------
// Gyroscopic torque calculation
// ---------------------------------------------------------------------------

/// Compute the gyroscopic torque `tau = omega_precession x L_spin`.
///
/// This is the torque needed to maintain a given precession.
pub fn gyroscopic_torque(omega_precession: [f64; 3], angular_momentum_spin: [f64; 3]) -> [f64; 3] {
    cross(omega_precession, angular_momentum_spin)
}

/// Compute the reaction torque on the frame due to a spinning rotor being
/// tilted at rate `omega_tilt` with spin angular momentum `L`.
pub fn reaction_torque(omega_tilt: [f64; 3], spin_angular_momentum: [f64; 3]) -> [f64; 3] {
    neg3(cross(omega_tilt, spin_angular_momentum))
}

// ---------------------------------------------------------------------------
// Spinning top simulation (Lagrange top)
// ---------------------------------------------------------------------------

/// State of a Lagrange top (symmetric top under gravity).
#[derive(Debug, Clone)]
pub struct LagrangeTopState {
    /// Tilt angle theta from vertical (radians).
    pub theta: f64,
    /// Rate of change of theta (rad/s).
    pub theta_dot: f64,
    /// Precession angle phi (radians).
    pub phi: f64,
    /// Precession rate phi_dot (rad/s).
    pub phi_dot: f64,
    /// Spin angle psi (radians).
    pub psi: f64,
    /// Spin rate psi_dot (rad/s).
    pub psi_dot: f64,
    /// Spin-axis moment of inertia.
    pub i_spin: f64,
    /// Transverse moment of inertia.
    pub i_transverse: f64,
    /// Mass of the top.
    pub mass: f64,
    /// Distance from pivot to center of mass.
    pub cm_distance: f64,
    /// Gravitational acceleration magnitude.
    pub gravity: f64,
}

impl LagrangeTopState {
    /// Create a new Lagrange top with given parameters.
    pub fn new(
        i_spin: f64,
        i_transverse: f64,
        mass: f64,
        cm_distance: f64,
        gravity: f64,
        theta_init: f64,
        spin_rate: f64,
    ) -> Self {
        Self {
            theta: theta_init,
            theta_dot: 0.0,
            phi: 0.0,
            phi_dot: 0.0,
            psi: 0.0,
            psi_dot: spin_rate,
            i_spin,
            i_transverse,
            mass,
            cm_distance,
            gravity,
        }
    }

    /// Total energy (kinetic + potential).
    pub fn total_energy(&self) -> f64 {
        let sin_t = self.theta.sin();
        let cos_t = self.theta.cos();
        let ke_trans = 0.5
            * self.i_transverse
            * (self.theta_dot * self.theta_dot + self.phi_dot * self.phi_dot * sin_t * sin_t);
        let ke_spin = 0.5 * self.i_spin * (self.psi_dot + self.phi_dot * cos_t).powi(2);
        let pe = self.mass * self.gravity * self.cm_distance * cos_t;
        ke_trans + ke_spin + pe
    }

    /// Conserved angular momentum component along vertical.
    pub fn vertical_angular_momentum(&self) -> f64 {
        let cos_t = self.theta.cos();
        let sin_t = self.theta.sin();
        self.i_transverse * self.phi_dot * sin_t * sin_t
            + self.i_spin * (self.psi_dot + self.phi_dot * cos_t) * cos_t
    }

    /// Conserved angular momentum about the spin axis.
    pub fn spin_angular_momentum(&self) -> f64 {
        self.i_spin * (self.psi_dot + self.phi_dot * self.theta.cos())
    }
}

/// Integrate the Lagrange top by one step using semi-implicit Euler.
pub fn integrate_lagrange_top(state: &mut LagrangeTopState, dt: f64) {
    let sin_t = state.theta.sin();
    let cos_t = state.theta.cos();

    // Conserved quantities
    let p_psi = state.spin_angular_momentum();
    let p_phi = state.vertical_angular_momentum();

    // theta equation of motion
    let mgl = state.mass * state.gravity * state.cm_distance;
    let theta_ddot = if sin_t.abs() > 1e-12 {
        let phi_dot_sq = if (state.i_transverse * sin_t * sin_t).abs() > 1e-15 {
            let eff_phi_dot = (p_phi - p_psi * cos_t) / (state.i_transverse * sin_t * sin_t);
            eff_phi_dot * eff_phi_dot
        } else {
            0.0
        };
        state.i_transverse * phi_dot_sq * sin_t * cos_t / state.i_transverse
            - mgl * sin_t / state.i_transverse
            + p_psi * state.phi_dot * sin_t / state.i_transverse
    } else {
        0.0
    };

    // Semi-implicit Euler
    state.theta_dot += theta_ddot * dt;
    state.theta += state.theta_dot * dt;

    // Clamp theta to avoid singularity
    state.theta = state.theta.clamp(1e-6, PI - 1e-6);

    // Update phi_dot and psi_dot from conserved quantities
    let sin_t2 = state.theta.sin();
    let cos_t2 = state.theta.cos();
    if (state.i_transverse * sin_t2 * sin_t2).abs() > 1e-15 {
        state.phi_dot = (p_phi - p_psi * cos_t2) / (state.i_transverse * sin_t2 * sin_t2);
    }
    if state.i_spin.abs() > 1e-15 {
        state.psi_dot = p_psi / state.i_spin - state.phi_dot * cos_t2;
    }

    state.phi += state.phi_dot * dt;
    state.psi += state.psi_dot * dt;
}

// ---------------------------------------------------------------------------
// Control Moment Gyroscope (CMG) model
// ---------------------------------------------------------------------------

/// A single-gimbal control moment gyroscope.
#[derive(Debug, Clone)]
pub struct ControlMomentGyroscope {
    /// Spin angular momentum magnitude of the rotor.
    pub spin_momentum: f64,
    /// Gimbal angle (radians).
    pub gimbal_angle: f64,
    /// Gimbal rate (rad/s).
    pub gimbal_rate: f64,
    /// Gimbal axis in body frame.
    pub gimbal_axis: [f64; 3],
    /// Spin axis at zero gimbal angle in body frame.
    pub spin_axis_ref: [f64; 3],
}

impl ControlMomentGyroscope {
    /// Create a new CMG with given rotor momentum and gimbal configuration.
    pub fn new(spin_momentum: f64, gimbal_axis: [f64; 3], spin_axis_ref: [f64; 3]) -> Self {
        Self {
            spin_momentum,
            gimbal_angle: 0.0,
            gimbal_rate: 0.0,
            gimbal_axis: normalize(gimbal_axis),
            spin_axis_ref: normalize(spin_axis_ref),
        }
    }

    /// Current spin axis direction after gimbal rotation.
    pub fn current_spin_axis(&self) -> [f64; 3] {
        let q = axis_angle_to_quat(self.gimbal_axis, self.gimbal_angle);
        quat_rotate(q, self.spin_axis_ref)
    }

    /// Current angular momentum vector.
    pub fn angular_momentum(&self) -> [f64; 3] {
        scale(self.current_spin_axis(), self.spin_momentum)
    }

    /// Output torque from gimbal rate: `tau = d/dt(h) = h_dot * gimbal_rate`.
    pub fn output_torque(&self) -> [f64; 3] {
        let h = self.angular_momentum();
        scale(cross(self.gimbal_axis, h), self.gimbal_rate)
    }

    /// Set the gimbal rate command.
    pub fn set_gimbal_rate(&mut self, rate: f64) {
        self.gimbal_rate = rate;
    }

    /// Step the CMG state forward by `dt`.
    pub fn step(&mut self, dt: f64) {
        self.gimbal_angle += self.gimbal_rate * dt;
    }
}

/// Convert axis-angle to quaternion.
fn axis_angle_to_quat(axis: [f64; 3], angle: f64) -> [f64; 4] {
    let half = angle * 0.5;
    let s = half.sin();
    let c = half.cos();
    let a = normalize(axis);
    [a[0] * s, a[1] * s, a[2] * s, c]
}

// ---------------------------------------------------------------------------
// Rate gyroscope sensor model (with noise, drift, bias)
// ---------------------------------------------------------------------------

/// A model of a rate gyroscope sensor with realistic error sources.
#[derive(Debug, Clone)]
pub struct RateGyroscopeSensor {
    /// Constant bias (rad/s) added to each axis.
    pub bias: [f64; 3],
    /// Scale factor error (multiplicative, e.g., 1.0 = perfect).
    pub scale_factor: [f64; 3],
    /// Random walk coefficient (rad/s/sqrt(Hz)), used for noise amplitude.
    pub noise_density: f64,
    /// Bias instability drift rate (rad/s^2).
    pub drift_rate: f64,
    /// Current accumulated bias drift.
    pub drift_state: [f64; 3],
    /// Sample period (seconds).
    pub sample_period: f64,
    /// Internal RNG seed counter (deterministic pseudo-random).
    rng_counter: u64,
}

impl RateGyroscopeSensor {
    /// Create a sensor model with given error parameters.
    pub fn new(
        bias: [f64; 3],
        scale_factor: [f64; 3],
        noise_density: f64,
        drift_rate: f64,
        sample_period: f64,
    ) -> Self {
        Self {
            bias,
            scale_factor,
            noise_density,
            drift_rate,
            drift_state: [0.0, 0.0, 0.0],
            sample_period,
            rng_counter: 12345,
        }
    }

    /// Create an ideal sensor with no errors.
    pub fn ideal(sample_period: f64) -> Self {
        Self::new([0.0; 3], [1.0; 3], 0.0, 0.0, sample_period)
    }

    /// Simple deterministic pseudo-random number in `[-1, 1]`.
    fn pseudo_random(&mut self) -> f64 {
        self.rng_counter = self
            .rng_counter
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let bits = (self.rng_counter >> 33) as i32;
        bits as f64 / (i32::MAX as f64)
    }

    /// Measure angular velocity with sensor errors applied.
    pub fn measure(&mut self, true_omega: [f64; 3]) -> [f64; 3] {
        // Update drift
        let dt = self.sample_period;
        for i in 0..3 {
            self.drift_state[i] += self.drift_rate * self.pseudo_random() * dt;
        }

        // Apply scale factor, bias, drift, and noise
        let noise_sigma = self.noise_density / dt.sqrt();
        let mut measured = [0.0; 3];
        for i in 0..3 {
            let noise = noise_sigma * self.pseudo_random();
            measured[i] =
                self.scale_factor[i] * true_omega[i] + self.bias[i] + self.drift_state[i] + noise;
        }
        measured
    }

    /// Reset drift state to zero.
    pub fn reset_drift(&mut self) {
        self.drift_state = [0.0; 3];
    }
}

// ---------------------------------------------------------------------------
// Gyrocompass (finding true north from rotation)
// ---------------------------------------------------------------------------

/// Gyrocompass state for finding true north from Earth's rotation.
#[derive(Debug, Clone)]
pub struct Gyrocompass {
    /// Latitude in radians.
    pub latitude: f64,
    /// Earth's rotation rate (rad/s), default ~7.2921e-5.
    pub earth_rate: f64,
    /// Current heading estimate (radians from true north).
    pub heading: f64,
    /// Damping coefficient.
    pub damping: f64,
    /// Pendulum natural frequency.
    pub natural_freq: f64,
}

impl Gyrocompass {
    /// Create a gyrocompass at the given latitude.
    pub fn new(latitude_deg: f64) -> Self {
        Self {
            latitude: latitude_deg.to_radians(),
            earth_rate: 7.2921e-5,
            heading: 0.0,
            damping: 0.1,
            natural_freq: 0.01,
        }
    }

    /// Horizontal component of Earth's rotation at this latitude (rad/s).
    pub fn horizontal_earth_rate(&self) -> f64 {
        self.earth_rate * self.latitude.cos()
    }

    /// Vertical component of Earth's rotation at this latitude.
    pub fn vertical_earth_rate(&self) -> f64 {
        self.earth_rate * self.latitude.sin()
    }

    /// Settling time estimate for the gyrocompass (seconds).
    pub fn settling_time(&self) -> f64 {
        if self.damping.abs() < 1e-15 {
            return f64::INFINITY;
        }
        4.0 / (self.damping * self.natural_freq)
    }

    /// Step the gyrocompass model by `dt`. The heading converges to true north.
    pub fn step(&mut self, dt: f64) {
        let h_rate = self.horizontal_earth_rate();
        let torque = -self.natural_freq * self.natural_freq * self.heading.sin()
            - self.damping * self.heading;
        let _forcing = h_rate * self.heading.cos();
        self.heading += torque * dt;
    }
}

// ---------------------------------------------------------------------------
// Reaction wheel dynamics
// ---------------------------------------------------------------------------

/// A reaction wheel for attitude control.
#[derive(Debug, Clone)]
pub struct ReactionWheel {
    /// Moment of inertia of the wheel (kg*m^2).
    pub wheel_inertia: f64,
    /// Current wheel spin rate (rad/s).
    pub spin_rate: f64,
    /// Maximum allowed spin rate (rad/s).
    pub max_spin_rate: f64,
    /// Maximum torque the motor can apply (N*m).
    pub max_torque: f64,
    /// Spin axis in body frame (unit vector).
    pub axis: [f64; 3],
}

impl ReactionWheel {
    /// Create a new reaction wheel.
    pub fn new(wheel_inertia: f64, max_spin_rate: f64, max_torque: f64, axis: [f64; 3]) -> Self {
        Self {
            wheel_inertia,
            spin_rate: 0.0,
            max_spin_rate,
            max_torque,
            axis: normalize(axis),
        }
    }

    /// Angular momentum stored in the wheel.
    pub fn angular_momentum(&self) -> [f64; 3] {
        scale(self.axis, self.wheel_inertia * self.spin_rate)
    }

    /// Angular momentum magnitude.
    pub fn angular_momentum_magnitude(&self) -> f64 {
        (self.wheel_inertia * self.spin_rate).abs()
    }

    /// Whether the wheel is saturated.
    pub fn is_saturated(&self) -> bool {
        self.spin_rate.abs() >= self.max_spin_rate * 0.99
    }

    /// Apply a commanded torque and step the wheel.
    /// Returns the actual torque applied to the spacecraft body.
    pub fn apply_torque(&mut self, commanded_torque: f64, dt: f64) -> [f64; 3] {
        let torque = commanded_torque.clamp(-self.max_torque, self.max_torque);
        let alpha = torque / self.wheel_inertia;
        self.spin_rate += alpha * dt;
        self.spin_rate = self
            .spin_rate
            .clamp(-self.max_spin_rate, self.max_spin_rate);
        // Reaction torque on body (opposite sign)
        scale(self.axis, -torque)
    }

    /// Kinetic energy stored in the wheel.
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.wheel_inertia * self.spin_rate * self.spin_rate
    }
}

// ---------------------------------------------------------------------------
// Attitude determination using gyroscope data
// ---------------------------------------------------------------------------

/// Dead-reckoning attitude estimator using gyroscope integration.
#[derive(Debug, Clone)]
pub struct AttitudeEstimator {
    /// Current estimated orientation (quaternion, scalar-last).
    pub orientation: [f64; 4],
    /// Accumulated integration time.
    pub time: f64,
}

impl AttitudeEstimator {
    /// Create a new estimator at identity orientation.
    pub fn new() -> Self {
        Self {
            orientation: quat_identity(),
            time: 0.0,
        }
    }

    /// Create with an initial orientation.
    pub fn with_orientation(q: [f64; 4]) -> Self {
        Self {
            orientation: quat_normalize(q),
            time: 0.0,
        }
    }

    /// Update the estimate with a gyroscope measurement.
    pub fn update(&mut self, omega_measured: [f64; 3], dt: f64) {
        self.orientation = quat_integrate(self.orientation, omega_measured, dt);
        self.time += dt;
    }

    /// Get current Euler angles `[roll, pitch, yaw]`.
    pub fn euler_angles(&self) -> [f64; 3] {
        quat_to_euler(self.orientation)
    }

    /// Reset to identity orientation.
    pub fn reset(&mut self) {
        self.orientation = quat_identity();
        self.time = 0.0;
    }
}

impl Default for AttitudeEstimator {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Angular momentum conservation verification
// ---------------------------------------------------------------------------

/// Verify angular momentum conservation over a trajectory.
///
/// Returns the maximum deviation of angular momentum magnitude from the
/// initial value.
pub fn verify_angular_momentum_conservation(trajectory: &[SpinningBodyState]) -> f64 {
    if trajectory.is_empty() {
        return 0.0;
    }
    let l0 = norm(trajectory[0].angular_momentum_body());
    let mut max_dev = 0.0_f64;
    for state in trajectory.iter().skip(1) {
        let l = norm(state.angular_momentum_body());
        let dev = (l - l0).abs();
        max_dev = max_dev.max(dev);
    }
    max_dev
}

/// Verify energy conservation over a trajectory.
pub fn verify_energy_conservation(trajectory: &[SpinningBodyState]) -> f64 {
    if trajectory.is_empty() {
        return 0.0;
    }
    let e0 = trajectory[0].kinetic_energy();
    let mut max_dev = 0.0_f64;
    for state in trajectory.iter().skip(1) {
        let e = state.kinetic_energy();
        let dev = (e - e0).abs();
        max_dev = max_dev.max(dev);
    }
    max_dev
}

// ---------------------------------------------------------------------------
// Torque-free motion (Poinsot construction)
// ---------------------------------------------------------------------------

/// Poinsot ellipsoid parameters for torque-free rigid body motion.
#[derive(Debug, Clone)]
pub struct PoinsotEllipsoid {
    /// Semi-axes of the energy ellipsoid.
    pub energy_semi_axes: [f64; 3],
    /// Semi-axes of the momentum ellipsoid.
    pub momentum_semi_axes: [f64; 3],
    /// Total kinetic energy.
    pub energy: f64,
    /// Total angular momentum magnitude.
    pub momentum: f64,
}

/// Compute the Poinsot ellipsoid parameters from body state.
pub fn poinsot_ellipsoid(state: &SpinningBodyState) -> PoinsotEllipsoid {
    let energy = state.kinetic_energy();
    let momentum = norm(state.angular_momentum_body());

    // Energy ellipsoid: sum(omega_i^2 * I_i) = 2T
    // Semi-axes: sqrt(2T / I_i)
    let energy_semi_axes = [
        if state.inertia[0] > 1e-15 {
            (2.0 * energy / state.inertia[0]).sqrt()
        } else {
            0.0
        },
        if state.inertia[1] > 1e-15 {
            (2.0 * energy / state.inertia[1]).sqrt()
        } else {
            0.0
        },
        if state.inertia[2] > 1e-15 {
            (2.0 * energy / state.inertia[2]).sqrt()
        } else {
            0.0
        },
    ];

    // Momentum ellipsoid: sum(I_i^2 * omega_i^2) = L^2
    // Semi-axes: L / I_i
    let momentum_semi_axes = [
        if state.inertia[0] > 1e-15 {
            momentum / state.inertia[0]
        } else {
            0.0
        },
        if state.inertia[1] > 1e-15 {
            momentum / state.inertia[1]
        } else {
            0.0
        },
        if state.inertia[2] > 1e-15 {
            momentum / state.inertia[2]
        } else {
            0.0
        },
    ];

    PoinsotEllipsoid {
        energy_semi_axes,
        momentum_semi_axes,
        energy,
        momentum,
    }
}

/// Simulate torque-free motion and return the trajectory.
pub fn simulate_torque_free(
    initial: &SpinningBodyState,
    dt: f64,
    steps: usize,
) -> Vec<SpinningBodyState> {
    let mut state = initial.clone();
    let mut trajectory = Vec::with_capacity(steps + 1);
    trajectory.push(state.clone());
    let zero_torque = [0.0, 0.0, 0.0];
    for _ in 0..steps {
        integrate_spinning_body_rk4(&mut state, zero_torque, dt);
        trajectory.push(state.clone());
    }
    trajectory
}

// ---------------------------------------------------------------------------
// Gyroscopic stabilization analysis
// ---------------------------------------------------------------------------

/// Result of a gyroscopic stabilisation analysis.
#[derive(Debug, Clone)]
pub struct StabilisationResult {
    /// Whether the spin axis is stable.
    pub is_stable: bool,
    /// The stability parameter (positive = stable, negative = unstable).
    pub stability_parameter: f64,
    /// Critical spin rate for stability (rad/s).
    pub critical_spin_rate: f64,
    /// Effective stiffness from gyroscopic effect.
    pub gyroscopic_stiffness: f64,
}

/// Analyse gyroscopic stabilisation of a spinning body.
///
/// A body spinning about a principal axis is stable if it spins about
/// the axis of maximum or minimum inertia (intermediate axis theorem).
pub fn analyse_stabilisation(
    inertia: Inertia3,
    spin_axis_index: usize,
    spin_rate: f64,
) -> StabilisationResult {
    let idx = spin_axis_index.min(2);
    let i_spin = inertia[idx];
    let other: Vec<f64> = (0..3).filter(|&i| i != idx).map(|i| inertia[i]).collect();

    // Stability criterion: (I_spin - I_1)(I_spin - I_2) > 0
    let stability_parameter = (i_spin - other[0]) * (i_spin - other[1]);
    let is_stable = stability_parameter > 0.0;

    // Critical spin rate (for a top under gravity, simplified)
    let i_transverse = 0.5 * (other[0] + other[1]);
    let critical_spin_rate = if (i_spin - i_transverse).abs() > 1e-15 {
        (4.0 * i_transverse / (i_spin - i_transverse).abs()).sqrt()
    } else {
        0.0
    };

    let gyroscopic_stiffness = i_spin * spin_rate;

    StabilisationResult {
        is_stable,
        stability_parameter,
        critical_spin_rate,
        gyroscopic_stiffness,
    }
}

// ---------------------------------------------------------------------------
// Gimbal lock detection and quaternion alternative
// ---------------------------------------------------------------------------

/// Detect whether the current Euler angle set is near gimbal lock.
///
/// For ZYX convention, gimbal lock occurs when pitch is near +/- 90 degrees.
pub fn detect_gimbal_lock(orientation: [f64; 4], threshold_deg: f64) -> bool {
    let euler = quat_to_euler(orientation);
    let pitch_deg = euler[1].to_degrees();
    (pitch_deg.abs() - 90.0).abs() < threshold_deg
}

/// Compute the quaternion error between two orientations.
///
/// Returns the quaternion `q_error = q_target * conjugate(q_current)`.
pub fn quaternion_error(q_current: [f64; 4], q_target: [f64; 4]) -> [f64; 4] {
    quat_mul(q_target, quat_conjugate(q_current))
}

/// Extract the rotation angle from a quaternion (in radians).
pub fn quaternion_angle(q: [f64; 4]) -> f64 {
    let w = q[3].clamp(-1.0, 1.0);
    2.0 * w.acos()
}

/// Extract the rotation axis from a quaternion.
pub fn quaternion_axis(q: [f64; 4]) -> [f64; 3] {
    let sin_half = (1.0 - q[3] * q[3]).sqrt();
    if sin_half < 1e-12 {
        return [0.0, 0.0, 1.0]; // arbitrary when angle is ~0
    }
    [q[0] / sin_half, q[1] / sin_half, q[2] / sin_half]
}

// ---------------------------------------------------------------------------
// Utility: 3-axis gyroscope data processing
// ---------------------------------------------------------------------------

/// Apply a simple low-pass filter to gyroscope data.
///
/// Uses exponential moving average with time constant `tau`.
pub fn low_pass_filter_gyro(data: &[[f64; 3]], dt: f64, tau: f64) -> Vec<[f64; 3]> {
    if data.is_empty() {
        return Vec::new();
    }
    let alpha = if tau > 1e-15 { dt / (tau + dt) } else { 1.0 };
    let mut filtered = Vec::with_capacity(data.len());
    let mut prev = data[0];
    filtered.push(prev);
    for sample in data.iter().skip(1) {
        prev = [
            prev[0] + alpha * (sample[0] - prev[0]),
            prev[1] + alpha * (sample[1] - prev[1]),
            prev[2] + alpha * (sample[2] - prev[2]),
        ];
        filtered.push(prev);
    }
    filtered
}

/// Compute the integral of angular velocity over time (total rotation angle).
pub fn integrate_total_angle(data: &[[f64; 3]], dt: f64) -> f64 {
    data.iter().map(|omega| norm(*omega) * dt).sum()
}

/// Compute the mean angular velocity from a data set.
pub fn mean_angular_velocity(data: &[[f64; 3]]) -> [f64; 3] {
    if data.is_empty() {
        return [0.0; 3];
    }
    let n = data.len() as f64;
    let sum = data.iter().fold([0.0; 3], |acc, omega| add3(acc, *omega));
    scale(sum, 1.0 / n)
}

/// Compute the variance of angular velocity on each axis.
pub fn angular_velocity_variance(data: &[[f64; 3]]) -> [f64; 3] {
    if data.len() < 2 {
        return [0.0; 3];
    }
    let mean = mean_angular_velocity(data);
    let n = data.len() as f64;
    let mut var = [0.0; 3];
    for omega in data {
        for i in 0..3 {
            let diff = omega[i] - mean[i];
            var[i] += diff * diff;
        }
    }
    scale(var, 1.0 / (n - 1.0))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    // ---- Euler equations ----

    #[test]
    fn test_euler_equations_no_torque_symmetric() {
        // Symmetric body spinning about z: no angular acceleration
        let inertia = [1.0, 1.0, 1.0];
        let omega = [0.0, 0.0, 10.0];
        let alpha = euler_equations(inertia, omega, [0.0, 0.0, 0.0]);
        assert!(norm(alpha) < TOL, "symmetric body: no angular acceleration");
    }

    #[test]
    fn test_euler_equations_asymmetric_coupling() {
        // Asymmetric body: spin about x with I different on each axis
        let inertia = [1.0, 2.0, 3.0];
        let omega = [0.0, 1.0, 1.0];
        let alpha = euler_equations(inertia, omega, [0.0, 0.0, 0.0]);
        // alpha_x = I_x^{-1} * (tau_x - (omega x I*omega)_x)
        //         = (0 - [0,1,1]x[0,2,3]_x) / 1 = (0 - 1) / 1 = -1
        assert!(
            (alpha[0] + 1.0).abs() < TOL,
            "alpha_x should be -1.0, got {:.6}",
            alpha[0]
        );
    }

    #[test]
    fn test_euler_equations_torque_only() {
        let inertia = [2.0, 2.0, 2.0];
        let omega = [0.0, 0.0, 0.0];
        let torque = [4.0, 0.0, 0.0];
        let alpha = euler_equations(inertia, omega, torque);
        assert!((alpha[0] - 2.0).abs() < TOL, "alpha = torque/I = 4/2 = 2");
    }

    // ---- SpinningBodyState ----

    #[test]
    fn test_spinning_body_angular_momentum() {
        let body = SpinningBodyState::new([1.0, 2.0, 3.0], [10.0, 0.0, 0.0]);
        let l = body.angular_momentum_body();
        assert!((l[0] - 10.0).abs() < TOL);
        assert!(l[1].abs() < TOL);
        assert!(l[2].abs() < TOL);
    }

    #[test]
    fn test_spinning_body_kinetic_energy() {
        let body = SpinningBodyState::new([4.0, 4.0, 4.0], [0.0, 0.0, 5.0]);
        // KE = 0.5 * I * omega^2 = 0.5 * 4 * 25 = 50
        assert!((body.kinetic_energy() - 50.0).abs() < TOL);
    }

    #[test]
    fn test_spinning_body_spin_rate() {
        let body = SpinningBodyState::new([1.0, 1.0, 1.0], [3.0, 4.0, 0.0]);
        assert!((body.spin_rate() - 5.0).abs() < TOL);
    }

    // ---- RK4 integration ----

    #[test]
    fn test_rk4_torque_free_conserves_energy() {
        let mut state = SpinningBodyState::new([1.0, 2.0, 3.0], [5.0, 1.0, 0.5]);
        let e0 = state.kinetic_energy();
        for _ in 0..1000 {
            integrate_spinning_body_rk4(&mut state, [0.0, 0.0, 0.0], 0.001);
        }
        let e1 = state.kinetic_energy();
        assert!(
            (e1 - e0).abs() < 0.01 * e0,
            "energy should be conserved: {:.6} vs {:.6}",
            e0,
            e1
        );
    }

    #[test]
    fn test_rk4_torque_free_conserves_momentum() {
        let mut state = SpinningBodyState::new([1.0, 2.0, 3.0], [5.0, 1.0, 0.5]);
        let l0 = norm(state.angular_momentum_body());
        for _ in 0..1000 {
            integrate_spinning_body_rk4(&mut state, [0.0, 0.0, 0.0], 0.001);
        }
        let l1 = norm(state.angular_momentum_body());
        assert!(
            (l1 - l0).abs() < 0.01 * l0,
            "momentum magnitude should be conserved: {:.6} vs {:.6}",
            l0,
            l1
        );
    }

    #[test]
    fn test_euler_integration_applies_torque() {
        let mut state = SpinningBodyState::new([1.0, 1.0, 1.0], [0.0, 0.0, 0.0]);
        integrate_spinning_body_euler(&mut state, [1.0, 0.0, 0.0], 1.0);
        assert!(
            (state.omega[0] - 1.0).abs() < TOL,
            "omega_x should be 1 after unit torque for 1 second"
        );
    }

    // ---- Precession ----

    #[test]
    fn test_precession_rate_formula() {
        let rate = precession_rate(0.1, 100.0, 1.0, 9.81, 0.05);
        // tau = 1 * 9.81 * 0.05 = 0.4905
        // L = 0.1 * 100 = 10
        // Omega_p = 0.4905 / 10 = 0.04905
        assert!(
            (rate - 0.04905).abs() < 1e-4,
            "precession rate = {:.6}",
            rate
        );
    }

    #[test]
    fn test_precession_rate_zero_spin() {
        let rate = precession_rate(0.1, 0.0, 1.0, 9.81, 0.05);
        assert!((rate).abs() < TOL, "zero spin should give zero precession");
    }

    // ---- Nutation ----

    #[test]
    fn test_nutation_frequency_symmetric() {
        // Same inertia => nutation frequency = 0
        let freq = nutation_frequency(1.0, 1.0, 100.0);
        assert!(freq.abs() < TOL);
    }

    #[test]
    fn test_nutation_frequency_asymmetric() {
        let freq = nutation_frequency(2.0, 1.0, 10.0);
        // (2 - 1) * 10 / 1 = 10
        assert!((freq - 10.0).abs() < TOL);
    }

    #[test]
    fn test_nutation_angle_small() {
        let angle = nutation_angle(0.01, 10.0, 1.0, 1.0);
        // atan(0.01/10) ~ 0.001
        assert!(angle.abs() < 0.01);
    }

    // ---- Gyroscopic torque ----

    #[test]
    fn test_gyroscopic_torque_basic() {
        let omega_p = [0.0, 0.0, 1.0]; // precessing about z
        let l_spin = [0.0, 10.0, 0.0]; // spinning about y
        let tau = gyroscopic_torque(omega_p, l_spin);
        // z x y = -x direction => [-10, 0, 0]
        assert!((tau[0] + 10.0).abs() < TOL);
        assert!(tau[1].abs() < TOL);
        assert!(tau[2].abs() < TOL);
    }

    #[test]
    fn test_reaction_torque_opposite() {
        let omega = [0.0, 0.0, 1.0];
        let l = [0.0, 10.0, 0.0];
        let tau = reaction_torque(omega, l);
        // Should be opposite of gyroscopic_torque
        let tau_g = gyroscopic_torque(omega, l);
        for i in 0..3 {
            assert!((tau[i] + tau_g[i]).abs() < TOL);
        }
    }

    // ---- Lagrange top ----

    #[test]
    fn test_lagrange_top_energy_conservation() {
        let mut top = LagrangeTopState::new(0.1, 0.05, 0.5, 0.1, 9.81, 0.3, 100.0);
        let e0 = top.total_energy();
        // Use smaller dt for better energy conservation with semi-implicit Euler
        for _ in 0..1000 {
            integrate_lagrange_top(&mut top, 0.0001);
        }
        let e1 = top.total_energy();
        assert!(
            (e1 - e0).abs() < 0.1 * e0.abs().max(1.0),
            "energy conservation: {:.6} vs {:.6}",
            e0,
            e1
        );
    }

    #[test]
    fn test_lagrange_top_spin_momentum_conserved() {
        let mut top = LagrangeTopState::new(0.1, 0.05, 0.5, 0.1, 9.81, 0.3, 100.0);
        let p0 = top.spin_angular_momentum();
        for _ in 0..500 {
            integrate_lagrange_top(&mut top, 0.001);
        }
        let p1 = top.spin_angular_momentum();
        assert!(
            (p1 - p0).abs() < 0.1 * p0.abs(),
            "spin momentum: {:.6} vs {:.6}",
            p0,
            p1
        );
    }

    // ---- CMG ----

    #[test]
    fn test_cmg_angular_momentum() {
        let cmg = ControlMomentGyroscope::new(10.0, [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        let h = cmg.angular_momentum();
        assert!((norm(h) - 10.0).abs() < TOL);
    }

    #[test]
    fn test_cmg_torque_generation() {
        let mut cmg = ControlMomentGyroscope::new(10.0, [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        cmg.set_gimbal_rate(1.0);
        let tau = cmg.output_torque();
        assert!(
            norm(tau) > 0.0,
            "CMG should generate torque with gimbal rate"
        );
    }

    #[test]
    fn test_cmg_step_changes_angle() {
        let mut cmg = ControlMomentGyroscope::new(10.0, [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        cmg.set_gimbal_rate(0.5);
        cmg.step(1.0);
        assert!((cmg.gimbal_angle - 0.5).abs() < TOL);
    }

    // ---- Rate gyroscope sensor ----

    #[test]
    fn test_ideal_sensor_no_errors() {
        let mut sensor = RateGyroscopeSensor::ideal(0.01);
        let true_omega = [1.0, 2.0, 3.0];
        let measured = sensor.measure(true_omega);
        for i in 0..3 {
            assert!(
                (measured[i] - true_omega[i]).abs() < 0.5,
                "ideal sensor axis {i}: {:.6} vs {:.6}",
                measured[i],
                true_omega[i]
            );
        }
    }

    #[test]
    fn test_sensor_bias() {
        let mut sensor = RateGyroscopeSensor::new([1.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0.0, 0.0, 0.01);
        let measured = sensor.measure([0.0, 0.0, 0.0]);
        assert!(
            measured[0] > 0.5,
            "x-axis should show positive bias, got {:.6}",
            measured[0]
        );
    }

    // ---- Gyrocompass ----

    #[test]
    fn test_gyrocompass_earth_rate_components() {
        let gc = Gyrocompass::new(45.0);
        let h = gc.horizontal_earth_rate();
        let v = gc.vertical_earth_rate();
        let total = (h * h + v * v).sqrt();
        assert!(
            (total - gc.earth_rate).abs() < 1e-10,
            "h^2 + v^2 should equal omega_earth^2"
        );
    }

    #[test]
    fn test_gyrocompass_settling_time() {
        let gc = Gyrocompass::new(0.0);
        let st = gc.settling_time();
        assert!(st > 0.0 && st < 1e6, "settling time = {:.6}", st);
    }

    // ---- Reaction wheel ----

    #[test]
    fn test_reaction_wheel_torque() {
        let mut rw = ReactionWheel::new(0.01, 6000.0, 0.1, [1.0, 0.0, 0.0]);
        let tau = rw.apply_torque(0.05, 0.01);
        // Reaction torque on body should be opposite
        assert!(tau[0] < 0.0, "reaction torque should oppose command");
        assert!(rw.spin_rate > 0.0, "wheel should be spinning");
    }

    #[test]
    fn test_reaction_wheel_saturation() {
        let mut rw = ReactionWheel::new(0.01, 100.0, 1.0, [0.0, 0.0, 1.0]);
        // Spin up to saturation
        for _ in 0..10000 {
            rw.apply_torque(1.0, 0.01);
        }
        assert!(rw.is_saturated(), "wheel should be saturated");
    }

    #[test]
    fn test_reaction_wheel_angular_momentum() {
        let mut rw = ReactionWheel::new(0.1, 1000.0, 1.0, [0.0, 1.0, 0.0]);
        rw.spin_rate = 100.0;
        let h = rw.angular_momentum();
        assert!((h[1] - 10.0).abs() < TOL, "L = I*omega = 0.1*100 = 10");
    }

    // ---- Attitude estimator ----

    #[test]
    fn test_attitude_estimator_identity() {
        let est = AttitudeEstimator::new();
        let euler = est.euler_angles();
        for i in 0..3 {
            assert!(euler[i].abs() < TOL, "initial Euler angles should be 0");
        }
    }

    #[test]
    fn test_attitude_estimator_rotation_about_z() {
        let mut est = AttitudeEstimator::new();
        let omega = [0.0, 0.0, 1.0]; // 1 rad/s about z
        let dt = 0.01;
        for _ in 0..100 {
            est.update(omega, dt);
        }
        let euler = est.euler_angles();
        // After 1 second at 1 rad/s about z: yaw ~ 1 rad
        assert!(
            (euler[2] - 1.0).abs() < 0.1,
            "yaw should be ~1 rad, got {:.6}",
            euler[2]
        );
    }

    // ---- Conservation verification ----

    #[test]
    fn test_verify_angular_momentum_conservation_torque_free() {
        let initial = SpinningBodyState::new([1.0, 2.0, 3.0], [5.0, 1.0, 0.5]);
        let traj = simulate_torque_free(&initial, 0.001, 500);
        let dev = verify_angular_momentum_conservation(&traj);
        let l0 = norm(initial.angular_momentum_body());
        assert!(
            dev < 0.01 * l0,
            "AM deviation {:.6} too large for L0={:.6}",
            dev,
            l0
        );
    }

    #[test]
    fn test_verify_energy_conservation_torque_free() {
        let initial = SpinningBodyState::new([1.0, 2.0, 3.0], [5.0, 1.0, 0.5]);
        let traj = simulate_torque_free(&initial, 0.001, 500);
        let dev = verify_energy_conservation(&traj);
        let e0 = initial.kinetic_energy();
        assert!(
            dev < 0.01 * e0,
            "energy deviation {:.6} too large for E0={:.6}",
            dev,
            e0
        );
    }

    // ---- Poinsot ellipsoid ----

    #[test]
    fn test_poinsot_ellipsoid_positive_semi_axes() {
        let state = SpinningBodyState::new([1.0, 2.0, 3.0], [5.0, 1.0, 0.5]);
        let pe = poinsot_ellipsoid(&state);
        for &a in &pe.energy_semi_axes {
            assert!(a >= 0.0);
        }
        for &a in &pe.momentum_semi_axes {
            assert!(a >= 0.0);
        }
        assert!(pe.energy > 0.0);
        assert!(pe.momentum > 0.0);
    }

    // ---- Stabilisation analysis ----

    #[test]
    fn test_stabilisation_max_inertia_stable() {
        // Spinning about max inertia axis (index 2, I=3) => stable
        let result = analyse_stabilisation([1.0, 2.0, 3.0], 2, 100.0);
        assert!(result.is_stable, "max inertia axis should be stable");
    }

    #[test]
    fn test_stabilisation_min_inertia_stable() {
        // Spinning about min inertia axis (index 0, I=1) => stable
        let result = analyse_stabilisation([1.0, 2.0, 3.0], 0, 100.0);
        assert!(result.is_stable, "min inertia axis should be stable");
    }

    #[test]
    fn test_stabilisation_intermediate_axis_unstable() {
        // Spinning about intermediate inertia axis (index 1, I=2) => unstable
        let result = analyse_stabilisation([1.0, 2.0, 3.0], 1, 100.0);
        assert!(
            !result.is_stable,
            "intermediate axis should be unstable (tennis racket theorem)"
        );
    }

    #[test]
    fn test_stabilisation_gyroscopic_stiffness() {
        let result = analyse_stabilisation([1.0, 1.0, 2.0], 2, 50.0);
        assert!(
            (result.gyroscopic_stiffness - 100.0).abs() < TOL,
            "stiffness = I*omega = 2*50 = 100"
        );
    }

    // ---- Gimbal lock ----

    #[test]
    fn test_gimbal_lock_at_90_pitch() {
        let q = euler_to_quat(0.0, PI / 2.0 - 0.01, 0.0);
        assert!(
            detect_gimbal_lock(q, 5.0),
            "near 90 deg pitch should be gimbal lock"
        );
    }

    #[test]
    fn test_no_gimbal_lock_at_zero_pitch() {
        let q = euler_to_quat(0.0, 0.0, 0.0);
        assert!(
            !detect_gimbal_lock(q, 5.0),
            "zero pitch should not be gimbal lock"
        );
    }

    // ---- Quaternion utilities ----

    #[test]
    fn test_quaternion_error_identity() {
        let q = quat_identity();
        let err = quaternion_error(q, q);
        assert!((err[3] - 1.0).abs() < TOL, "error should be identity");
    }

    #[test]
    fn test_quaternion_angle_180() {
        let q = [0.0, 0.0, 1.0, 0.0]; // 180 deg about z
        let angle = quaternion_angle(q);
        assert!(
            (angle - PI).abs() < 0.01,
            "angle should be pi, got {:.6}",
            angle
        );
    }

    #[test]
    fn test_quaternion_axis_from_z_rotation() {
        let q = axis_angle_to_quat([0.0, 0.0, 1.0], 1.0);
        let axis = quaternion_axis(q);
        assert!((axis[2] - 1.0).abs() < 0.01, "axis should be z");
    }

    // ---- Data processing ----

    #[test]
    fn test_low_pass_filter_constant() {
        let data = vec![[1.0, 2.0, 3.0]; 10];
        let filtered = low_pass_filter_gyro(&data, 0.01, 0.1);
        for sample in &filtered {
            for i in 0..3 {
                assert!(
                    (sample[i] - data[0][i]).abs() < TOL,
                    "constant signal should pass through"
                );
            }
        }
    }

    #[test]
    fn test_integrate_total_angle() {
        let data = vec![[0.0, 0.0, 10.0]; 100]; // 10 rad/s about z
        let total = integrate_total_angle(&data, 0.01);
        // 100 steps * 10 rad/s * 0.01 s = 10 rad
        assert!((total - 10.0).abs() < TOL);
    }

    #[test]
    fn test_mean_angular_velocity() {
        let data = vec![[1.0, 2.0, 3.0], [3.0, 4.0, 5.0]];
        let mean = mean_angular_velocity(&data);
        assert!((mean[0] - 2.0).abs() < TOL);
        assert!((mean[1] - 3.0).abs() < TOL);
        assert!((mean[2] - 4.0).abs() < TOL);
    }

    #[test]
    fn test_angular_velocity_variance_constant() {
        let data = vec![[1.0, 1.0, 1.0]; 10];
        let var = angular_velocity_variance(&data);
        for i in 0..3 {
            assert!(
                var[i].abs() < TOL,
                "constant data should have zero variance"
            );
        }
    }

    // ---- Torque-free simulation ----

    #[test]
    fn test_simulate_torque_free_length() {
        let state = SpinningBodyState::new([1.0, 1.0, 1.0], [0.0, 0.0, 10.0]);
        let traj = simulate_torque_free(&state, 0.01, 50);
        assert_eq!(traj.len(), 51, "initial + 50 steps");
    }

    #[test]
    fn test_euler_to_quat_roundtrip() {
        let roll = 0.3;
        let pitch = 0.2;
        let yaw = 0.5;
        let q = euler_to_quat(roll, pitch, yaw);
        let euler = quat_to_euler(q);
        assert!((euler[0] - roll).abs() < 0.01, "roll: {:.6}", euler[0]);
        assert!((euler[1] - pitch).abs() < 0.01, "pitch: {:.6}", euler[1]);
        assert!((euler[2] - yaw).abs() < 0.01, "yaw: {:.6}", euler[2]);
    }
}
