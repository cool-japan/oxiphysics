#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Aerospace maneuver dynamics.
//!
//! Provides 6-DOF aircraft rigid-body dynamics, flight control surfaces,
//! inertial navigation, aerodynamic force/moment computation, maneuver
//! planning (Dubins path, energy-altitude trades), and recovery systems
//! (parachute deployment, landing impact).
//!
//! # Example
//!
//! ```no_run
//! use oxiphysics_rigid::aerospace_maneuver::{AircraftDynamics, AerodynamicForces};
//!
//! let ac = AircraftDynamics::new(8000.0, [4000.0, 3000.0, 5000.0]);
//! assert!(ac.mass > 0.0);
//!
//! let forces = AerodynamicForces::new(0.3, 0.05, 0.02, 1.225);
//! let q = forces.dynamic_pressure(100.0);
//! assert!(q > 0.0);
//! ```

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ── Vector helpers ────────────────────────────────────────────────────────────

/// Cross product of two 3-vectors.
#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Dot product of two 3-vectors.
#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// L2 norm of a 3-vector.
#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Add two 3-vectors.
#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors.
#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector by a scalar.
#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Normalize a 3-vector; returns zero-vector if the norm is tiny.
#[inline]
fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n < 1e-30 {
        [0.0; 3]
    } else {
        [a[0] / n, a[1] / n, a[2] / n]
    }
}

/// Normalize a quaternion `[x, y, z, w]`.
#[inline]
fn quat_normalize(q: [f64; 4]) -> [f64; 4] {
    let n = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    if n < 1e-30 {
        [0.0, 0.0, 0.0, 1.0]
    } else {
        [q[0] / n, q[1] / n, q[2] / n, q[3] / n]
    }
}

/// Rotate a 3-vector by a unit quaternion `[x, y, z, w]`.
#[inline]
fn quat_rotate(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    let [qx, qy, qz, qw] = q;
    let [vx, vy, vz] = v;
    // t = 2 * cross(q.xyz, v)
    let tx = 2.0 * (qy * vz - qz * vy);
    let ty = 2.0 * (qz * vx - qx * vz);
    let tz = 2.0 * (qx * vy - qy * vx);
    [
        vx + qw * tx + qy * tz - qz * ty,
        vy + qw * ty + qz * tx - qx * tz,
        vz + qw * tz + qx * ty - qy * tx,
    ]
}

/// Integrate quaternion attitude with angular velocity `omega` (body frame) over `dt`.
#[inline]
fn quat_integrate(q: [f64; 4], omega: [f64; 3], dt: f64) -> [f64; 4] {
    let [qx, qy, qz, qw] = q;
    let [ox, oy, oz] = omega;
    let half_dt = 0.5 * dt;
    let dqx = half_dt * (qw * ox + qy * oz - qz * oy);
    let dqy = half_dt * (qw * oy + qz * ox - qx * oz);
    let dqz = half_dt * (qw * oz + qx * oy - qy * ox);
    let dqw = half_dt * (-qx * ox - qy * oy - qz * oz);
    quat_normalize([qx + dqx, qy + dqy, qz + dqz, qw + dqw])
}

/// Standard atmosphere: density (kg/m³) at altitude h (m) using ISA model.
pub fn isa_density(h: f64) -> f64 {
    const RHO0: f64 = 1.225; // sea-level density kg/m³
    const H_SCALE: f64 = 8500.0; // approximate scale height m
    RHO0 * (-h / H_SCALE).exp()
}

/// Standard atmosphere: speed of sound (m/s) at altitude h (m).
pub fn isa_speed_of_sound(h: f64) -> f64 {
    // Simple linear lapse up to 11 km, constant above
    let t = if h < 11_000.0 {
        288.15 - 0.0065 * h
    } else {
        216.65
    };
    (1.4 * 287.058 * t).sqrt()
}

// ── AircraftDynamics ──────────────────────────────────────────────────────────

/// Six-degree-of-freedom rigid-body aircraft dynamics.
///
/// State: position (m), velocity (m/s) body-frame, quaternion attitude,
/// angular-velocity (rad/s) body-frame, altitude (m).
#[derive(Debug, Clone)]
pub struct AircraftDynamics {
    /// Aircraft total mass (kg).
    pub mass: f64,
    /// Moment of inertia diagonal \[Ixx, Iyy, Izz\] (kg·m²).
    pub inertia: [f64; 3],
    /// Position in inertial NED frame (m).
    pub position: [f64; 3],
    /// Velocity in body frame (m/s): \[u, v, w\] = \[fore, right, down\].
    pub velocity_body: [f64; 3],
    /// Attitude quaternion \[x, y, z, w\] (body ← NED).
    pub attitude: [f64; 4],
    /// Angular velocity in body frame \[p, q, r\] (rad/s).
    pub omega: [f64; 3],
    /// Accumulated force in body frame (N).
    pub force_body: [f64; 3],
    /// Accumulated torque in body frame (N·m).
    pub torque_body: [f64; 3],
    /// Gravitational acceleration magnitude (m/s²).
    pub gravity: f64,
}

impl AircraftDynamics {
    /// Create a new aircraft with given mass and inertia tensor diagonal.
    pub fn new(mass: f64, inertia: [f64; 3]) -> Self {
        Self {
            mass,
            inertia,
            position: [0.0; 3],
            velocity_body: [0.0; 3],
            attitude: [0.0, 0.0, 0.0, 1.0], // identity
            omega: [0.0; 3],
            force_body: [0.0; 3],
            torque_body: [0.0; 3],
            gravity: 9.80665,
        }
    }

    /// True airspeed (m/s) — magnitude of body-frame velocity.
    pub fn true_airspeed(&self) -> f64 {
        norm3(self.velocity_body)
    }

    /// Angle of attack α (rad) = atan2(w, u).
    pub fn angle_of_attack(&self) -> f64 {
        let [u, _v, w] = self.velocity_body;
        w.atan2(u)
    }

    /// Sideslip angle β (rad) = asin(v / V).
    pub fn sideslip(&self) -> f64 {
        let v_total = self.true_airspeed();
        if v_total < 1e-6 {
            0.0
        } else {
            (self.velocity_body[1] / v_total).clamp(-1.0, 1.0).asin()
        }
    }

    /// Altitude above reference (m) — negative of NED down component.
    pub fn altitude(&self) -> f64 {
        -self.position[2]
    }

    /// Apply external force in body frame.
    pub fn apply_force_body(&mut self, f: [f64; 3]) {
        self.force_body = add3(self.force_body, f);
    }

    /// Apply external torque in body frame.
    pub fn apply_torque_body(&mut self, t: [f64; 3]) {
        self.torque_body = add3(self.torque_body, t);
    }

    /// Clear accumulated forces and torques.
    pub fn clear_forces(&mut self) {
        self.force_body = [0.0; 3];
        self.torque_body = [0.0; 3];
    }

    /// Add gravitational body force (NED gravity → body frame).
    pub fn add_gravity_force(&mut self) {
        // Gravity in NED: [0, 0, m*g]
        let g_ned = [0.0, 0.0, self.mass * self.gravity];
        // Rotate gravity from NED into body frame (conjugate of attitude).
        let q_conj = [
            -self.attitude[0],
            -self.attitude[1],
            -self.attitude[2],
            self.attitude[3],
        ];
        let g_body = quat_rotate(q_conj, g_ned);
        self.force_body = add3(self.force_body, g_body);
    }

    /// Integrate equations of motion for one time step `dt` (s).
    ///
    /// Uses semi-implicit Euler (velocity → position, angular vel → attitude).
    pub fn step(&mut self, dt: f64) {
        let [fx, fy, fz] = self.force_body;
        let [tx, ty, tz] = self.torque_body;
        let [ixx, iyy, izz] = self.inertia;
        let [p, q, r] = self.omega;
        let [u, v, w] = self.velocity_body;

        // Linear acceleration in body frame
        let ax = fx / self.mass + r * v - q * w;
        let ay = fy / self.mass - r * u + p * w;
        let az = fz / self.mass + q * u - p * v;

        // Angular acceleration (Euler's equations)
        let dp = (tx - (izz - iyy) * q * r) / ixx;
        let dq = (ty - (ixx - izz) * p * r) / iyy;
        let dr = (tz - (iyy - ixx) * p * q) / izz;

        // Update velocities
        self.velocity_body = [u + ax * dt, v + ay * dt, w + az * dt];
        self.omega = [p + dp * dt, q + dq * dt, r + dr * dt];

        // Rotate velocity to NED and integrate position
        let vel_ned = quat_rotate(self.attitude, self.velocity_body);
        self.position = add3(self.position, scale3(vel_ned, dt));

        // Integrate attitude
        self.attitude = quat_integrate(self.attitude, self.omega, dt);

        self.clear_forces();
    }

    /// Compute Euler angles \[roll φ, pitch θ, yaw ψ\] (rad) from quaternion.
    pub fn euler_angles(&self) -> [f64; 3] {
        let [qx, qy, qz, qw] = self.attitude;
        let phi = (2.0 * (qw * qx + qy * qz)).atan2(1.0 - 2.0 * (qx * qx + qy * qy));
        let theta = (2.0 * (qw * qy - qz * qx)).clamp(-1.0, 1.0).asin();
        let psi = (2.0 * (qw * qz + qx * qy)).atan2(1.0 - 2.0 * (qy * qy + qz * qz));
        [phi, theta, psi]
    }

    /// Load factor n = lift / weight.
    pub fn load_factor(&self, lift: f64) -> f64 {
        lift / (self.mass * self.gravity)
    }

    /// Kinetic energy (J).
    pub fn kinetic_energy(&self) -> f64 {
        let v2 = dot3(self.velocity_body, self.velocity_body);
        let [ixx, iyy, izz] = self.inertia;
        let [p, q, r] = self.omega;
        0.5 * self.mass * v2 + 0.5 * (ixx * p * p + iyy * q * q + izz * r * r)
    }

    /// Potential energy (J) at current altitude (NED convention: positive up).
    pub fn potential_energy(&self) -> f64 {
        self.mass * self.gravity * self.altitude()
    }

    /// Total mechanical energy (J).
    pub fn total_energy(&self) -> f64 {
        self.kinetic_energy() + self.potential_energy()
    }
}

// ── FlightControlSurface ──────────────────────────────────────────────────────

/// Aerodynamic control surface with first-order actuator dynamics.
///
/// Models elevator, aileron, rudder, or flap with rate limiting and
/// actuator lag.
#[derive(Debug, Clone)]
pub struct FlightControlSurface {
    /// Surface type label.
    pub kind: SurfaceKind,
    /// Current deflection angle (rad).
    pub deflection: f64,
    /// Commanded deflection angle (rad).
    pub command: f64,
    /// Maximum deflection magnitude (rad).
    pub max_deflection: f64,
    /// Maximum slew rate (rad/s).
    pub max_rate: f64,
    /// First-order actuator time constant (s). 0 = ideal.
    pub time_constant: f64,
    /// Control effectiveness coefficient (N·m / rad).
    pub effectiveness: f64,
}

/// Kind of flight control surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceKind {
    /// Elevator — pitch control.
    Elevator,
    /// Aileron — roll control.
    Aileron,
    /// Rudder — yaw control.
    Rudder,
    /// Flap — high-lift device.
    Flap,
    /// Spoiler — drag / roll control.
    Spoiler,
}

impl FlightControlSurface {
    /// Construct a new control surface.
    pub fn new(
        kind: SurfaceKind,
        max_deflection: f64,
        max_rate: f64,
        time_constant: f64,
        effectiveness: f64,
    ) -> Self {
        Self {
            kind,
            deflection: 0.0,
            command: 0.0,
            max_deflection,
            max_rate,
            time_constant,
            effectiveness,
        }
    }

    /// Set a new commanded deflection (will be clamped to limits).
    pub fn set_command(&mut self, cmd: f64) {
        self.command = cmd.clamp(-self.max_deflection, self.max_deflection);
    }

    /// Advance actuator state by `dt` seconds.
    ///
    /// Applies slew-rate limiting and first-order lag.
    pub fn step(&mut self, dt: f64) {
        if self.time_constant < 1e-12 {
            // Ideal actuator — apply rate limit only
            let error = self.command - self.deflection;
            let max_delta = self.max_rate * dt;
            self.deflection += error.clamp(-max_delta, max_delta);
        } else {
            // First-order lag
            let tau = self.time_constant;
            let target_rate = (self.command - self.deflection) / tau;
            let actual_rate = target_rate.clamp(-self.max_rate, self.max_rate);
            self.deflection = (self.deflection + actual_rate * dt)
                .clamp(-self.max_deflection, self.max_deflection);
        }
    }

    /// Moment generated by the surface (N·m).
    pub fn moment(&self) -> f64 {
        self.effectiveness * self.deflection
    }

    /// Lift increment due to flap deflection using simplified model.
    ///
    /// Returns ΔCL per radian of deflection times dynamic pressure times area.
    pub fn delta_cl(&self, cl_delta: f64) -> f64 {
        cl_delta * self.deflection
    }
}

// ── InertialNavigation ────────────────────────────────────────────────────────

/// Inertial measurement unit state and navigation filter.
///
/// Integrates accelerometer and gyro measurements to estimate position,
/// velocity, and attitude; applies a complementary filter for drift correction.
#[derive(Debug, Clone)]
pub struct InertialNavigation {
    /// Estimated position (m) in NED frame.
    pub position: [f64; 3],
    /// Estimated velocity (m/s) in NED frame.
    pub velocity: [f64; 3],
    /// Estimated attitude quaternion \[x, y, z, w\].
    pub attitude: [f64; 4],
    /// Gyroscope bias (rad/s).
    pub gyro_bias: [f64; 3],
    /// Accelerometer bias (m/s²).
    pub accel_bias: [f64; 3],
    /// Complementary filter coefficient (0–1); higher → trust accelerometer more.
    pub alpha: f64,
    /// Accumulated dead-reckoning distance (m).
    pub dead_reckoning_dist: f64,
}

impl InertialNavigation {
    /// Create a new INS at the given initial position and attitude.
    pub fn new(position: [f64; 3], attitude: [f64; 4]) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            attitude: quat_normalize(attitude),
            gyro_bias: [0.0; 3],
            accel_bias: [0.0; 3],
            alpha: 0.02,
            dead_reckoning_dist: 0.0,
        }
    }

    /// Update INS with raw IMU measurements.
    ///
    /// * `accel_body` — accelerometer reading in body frame (m/s², includes gravity).
    /// * `gyro_body` — gyroscope reading in body frame (rad/s).
    /// * `dt` — time step (s).
    pub fn update_imu(&mut self, accel_body: [f64; 3], gyro_body: [f64; 3], dt: f64) {
        // Correct for bias
        let accel = sub3(accel_body, self.accel_bias);
        let gyro = sub3(gyro_body, self.gyro_bias);

        // Integrate attitude
        self.attitude = quat_integrate(self.attitude, gyro, dt);

        // Rotate accelerometer from body to NED; subtract gravity
        let accel_ned = quat_rotate(self.attitude, accel);
        let accel_inertial = [accel_ned[0], accel_ned[1], accel_ned[2] - 9.80665];

        // Integrate velocity and position
        self.velocity = add3(self.velocity, scale3(accel_inertial, dt));
        let delta_pos = scale3(self.velocity, dt);
        self.position = add3(self.position, delta_pos);
        self.dead_reckoning_dist += norm3(delta_pos);
    }

    /// Complementary filter correction using an external attitude reference.
    ///
    /// Blends IMU attitude with `reference_attitude` using coefficient `alpha`.
    pub fn complementary_correction(&mut self, reference_attitude: [f64; 4]) {
        let a = self.alpha;
        let q = self.attitude;
        let r = quat_normalize(reference_attitude);
        // Linear blend (approximate SLERP for small differences)
        let blended = [
            (1.0 - a) * q[0] + a * r[0],
            (1.0 - a) * q[1] + a * r[1],
            (1.0 - a) * q[2] + a * r[2],
            (1.0 - a) * q[3] + a * r[3],
        ];
        self.attitude = quat_normalize(blended);
    }

    /// Apply a GPS position fix to correct accumulated drift.
    pub fn gps_correction(&mut self, gps_position: [f64; 3], weight: f64) {
        let w = weight.clamp(0.0, 1.0);
        for i in 0..3 {
            self.position[i] = (1.0 - w) * self.position[i] + w * gps_position[i];
        }
    }

    /// Current altitude (m) — negative of NED down component.
    pub fn altitude(&self) -> f64 {
        -self.position[2]
    }

    /// Horizontal speed (m/s).
    pub fn horizontal_speed(&self) -> f64 {
        (self.velocity[0] * self.velocity[0] + self.velocity[1] * self.velocity[1]).sqrt()
    }

    /// Euler angles \[roll, pitch, yaw\] (rad).
    pub fn euler_angles(&self) -> [f64; 3] {
        let [qx, qy, qz, qw] = self.attitude;
        let phi = (2.0 * (qw * qx + qy * qz)).atan2(1.0 - 2.0 * (qx * qx + qy * qy));
        let theta = (2.0 * (qw * qy - qz * qx)).clamp(-1.0, 1.0).asin();
        let psi = (2.0 * (qw * qz + qx * qy)).atan2(1.0 - 2.0 * (qy * qy + qz * qz));
        [phi, theta, psi]
    }
}

// ── AerodynamicForces ─────────────────────────────────────────────────────────

/// Aerodynamic force/moment computation using stability derivatives.
///
/// Computes lift, drag, and side force from aerodynamic coefficients,
/// dynamic pressure, and reference geometry.
#[derive(Debug, Clone)]
pub struct AerodynamicForces {
    /// Lift coefficient at zero AoA, CL0.
    pub cl0: f64,
    /// Lift curve slope, CLα (per rad).
    pub cl_alpha: f64,
    /// Zero-lift drag coefficient, CD0.
    pub cd0: f64,
    /// Induced drag factor k (CD = CD0 + k·CL²).
    pub oswald_k: f64,
    /// Side force coefficient per sideslip, CYβ (per rad).
    pub cy_beta: f64,
    /// Pitch moment coefficient, Cm0.
    pub cm0: f64,
    /// Pitch moment slope, Cmα (per rad).
    pub cm_alpha: f64,
    /// Roll moment slope, Clβ (dihedral effect, per rad).
    pub cl_beta: f64,
    /// Yaw moment slope, Cnβ (weathercock stability, per rad).
    pub cn_beta: f64,
    /// Air density (kg/m³).
    pub rho: f64,
    /// Wing reference area (m²).
    pub s_ref: f64,
    /// Mean aerodynamic chord (m).
    pub c_bar: f64,
    /// Wing span (m).
    pub b: f64,
}

impl AerodynamicForces {
    /// Create with basic lift/drag/moment coefficients.
    ///
    /// # Arguments
    /// * `cl_alpha` — lift curve slope (per rad).
    /// * `cd0` — zero-lift drag coefficient.
    /// * `cm_alpha` — pitch-moment slope (per rad).
    /// * `rho` — air density (kg/m³).
    pub fn new(cl_alpha: f64, cd0: f64, cm_alpha: f64, rho: f64) -> Self {
        Self {
            cl0: 0.3,
            cl_alpha,
            cd0,
            oswald_k: 0.04,
            cy_beta: -0.8,
            cm0: 0.0,
            cm_alpha,
            cl_beta: -0.1,
            cn_beta: 0.15,
            rho,
            s_ref: 30.0,
            c_bar: 2.5,
            b: 12.0,
        }
    }

    /// Dynamic pressure q = 0.5 * ρ * V² (Pa).
    pub fn dynamic_pressure(&self, airspeed: f64) -> f64 {
        0.5 * self.rho * airspeed * airspeed
    }

    /// Lift coefficient at angle of attack `alpha` (rad).
    pub fn cl(&self, alpha: f64) -> f64 {
        self.cl0 + self.cl_alpha * alpha
    }

    /// Drag coefficient using parabolic polar.
    pub fn cd(&self, alpha: f64) -> f64 {
        let cl = self.cl(alpha);
        self.cd0 + self.oswald_k * cl * cl
    }

    /// Side-force coefficient at sideslip `beta` (rad).
    pub fn cy(&self, beta: f64) -> f64 {
        self.cy_beta * beta
    }

    /// Lift force (N) in wind axes.
    pub fn lift(&self, airspeed: f64, alpha: f64) -> f64 {
        self.dynamic_pressure(airspeed) * self.s_ref * self.cl(alpha)
    }

    /// Drag force (N) in wind axes.
    pub fn drag(&self, airspeed: f64, alpha: f64) -> f64 {
        self.dynamic_pressure(airspeed) * self.s_ref * self.cd(alpha)
    }

    /// Side force (N).
    pub fn side_force(&self, airspeed: f64, beta: f64) -> f64 {
        self.dynamic_pressure(airspeed) * self.s_ref * self.cy(beta)
    }

    /// Pitching moment (N·m).
    pub fn pitching_moment(&self, airspeed: f64, alpha: f64) -> f64 {
        let q_bar = self.dynamic_pressure(airspeed);
        q_bar * self.s_ref * self.c_bar * (self.cm0 + self.cm_alpha * alpha)
    }

    /// Rolling moment due to dihedral effect (N·m).
    pub fn rolling_moment(&self, airspeed: f64, beta: f64) -> f64 {
        self.dynamic_pressure(airspeed) * self.s_ref * self.b * self.cl_beta * beta
    }

    /// Yawing moment due to weathercock stability (N·m).
    pub fn yawing_moment(&self, airspeed: f64, beta: f64) -> f64 {
        self.dynamic_pressure(airspeed) * self.s_ref * self.b * self.cn_beta * beta
    }

    /// Compute full aerodynamic force vector in wind axes \[Drag, SideForce, Lift\].
    ///
    /// Returns `[Fx_wind, Fy_wind, Fz_wind]` where +Fz_wind = lift.
    pub fn wind_axis_forces(&self, airspeed: f64, alpha: f64, beta: f64) -> [f64; 3] {
        let l = self.lift(airspeed, alpha);
        let d = self.drag(airspeed, alpha);
        let y = self.side_force(airspeed, beta);
        [-d, y, l]
    }

    /// Convert wind-axes forces to body-axes.
    ///
    /// Rotation by α about Y then β about Z.
    pub fn to_body_forces(&self, airspeed: f64, alpha: f64, beta: f64) -> [f64; 3] {
        let [fx_w, fy_w, fz_w] = self.wind_axis_forces(airspeed, alpha, beta);
        let ca = alpha.cos();
        let sa = alpha.sin();
        let cb = beta.cos();
        let sb = beta.sin();
        // Wind → stability (rotate by −α about Y)
        // stability → body (rotate by −β about Z)
        // Combined rotation
        let fx = ca * cb * fx_w - sb * fy_w - sa * cb * fz_w;
        let fy = ca * sb * fx_w + cb * fy_w - sa * sb * fz_w;
        let fz = sa * fx_w + ca * fz_w;
        [fx, fy, fz]
    }

    /// L/D ratio (lift-to-drag).
    pub fn lift_to_drag(&self, alpha: f64) -> f64 {
        let cl = self.cl(alpha);
        let cd = self.cd(alpha);
        if cd.abs() < 1e-30 { 0.0 } else { cl / cd }
    }

    /// Angle of attack for maximum L/D.
    ///
    /// Analytical solution: CLα·α* = sqrt(CD0 / k) − CL0.
    pub fn alpha_best_ld(&self) -> f64 {
        let cl_star = (self.cd0 / self.oswald_k).sqrt();
        (cl_star - self.cl0) / self.cl_alpha.max(1e-6)
    }

    /// Stall AoA (approximate): where CL = CLmax ≈ 1.6.
    pub fn stall_alpha(&self) -> f64 {
        (1.6 - self.cl0) / self.cl_alpha.max(1e-6)
    }
}

// ── ManeuverPlanning ──────────────────────────────────────────────────────────

/// Dubins path segment type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DubinsSegment {
    /// Left (counter-clockwise) arc.
    Left,
    /// Straight.
    Straight,
    /// Right (clockwise) arc.
    Right,
}

/// A three-segment Dubins path.
#[derive(Debug, Clone)]
pub struct DubinsPath {
    /// Sequence of segment types (always length 3).
    pub segments: [DubinsSegment; 3],
    /// Length of each segment (rad for turns, m for straights).
    pub lengths: [f64; 3],
    /// Turning radius (m).
    pub radius: f64,
}

impl DubinsPath {
    /// Total path length (m).
    pub fn total_length(&self) -> f64 {
        self.lengths[0] * self.radius + self.lengths[1] + self.lengths[2] * self.radius
    }
}

/// Maneuver planning: Dubins paths, minimum-time turns, energy-altitude trades.
#[derive(Debug, Clone)]
pub struct ManeuverPlanning {
    /// Minimum turning radius (m).
    pub min_radius: f64,
    /// Maximum bank angle (rad).
    pub max_bank: f64,
    /// Nominal airspeed (m/s).
    pub airspeed: f64,
    /// Gravitational acceleration (m/s²).
    pub gravity: f64,
}

impl ManeuverPlanning {
    /// Create a maneuver planner.
    pub fn new(min_radius: f64, max_bank: f64, airspeed: f64) -> Self {
        Self {
            min_radius,
            max_bank,
            airspeed,
            gravity: 9.80665,
        }
    }

    /// Compute the shortest Dubins RSR path between two configurations.
    ///
    /// Configurations: `(x, y, heading)` in radians.
    pub fn dubins_rsr(&self, start: [f64; 3], end: [f64; 3]) -> DubinsPath {
        let r = self.min_radius;
        let [x1, y1, h1] = start;
        let [x2, y2, h2] = end;

        // Centers of right-turn circles
        let cx1 = x1 + r * (h1 - PI / 2.0).cos();
        let cy1 = y1 + r * (h1 - PI / 2.0).sin();
        let cx2 = x2 + r * (h2 - PI / 2.0).cos();
        let cy2 = y2 + r * (h2 - PI / 2.0).sin();

        let d = ((cx2 - cx1).powi(2) + (cy2 - cy1).powi(2)).sqrt();
        let theta = (cy2 - cy1).atan2(cx2 - cx1);

        let t1 = (theta - h1 + PI / 2.0).rem_euclid(2.0 * PI);
        let t2 = (h2 - theta + PI / 2.0).rem_euclid(2.0 * PI);

        DubinsPath {
            segments: [
                DubinsSegment::Right,
                DubinsSegment::Straight,
                DubinsSegment::Right,
            ],
            lengths: [t1, d, t2],
            radius: r,
        }
    }

    /// Minimum-time level turn: time to achieve a heading change of `delta_psi` (rad).
    ///
    /// Uses maximum bank angle to minimize turn radius.
    pub fn min_time_turn(&self, delta_psi: f64) -> f64 {
        let bank = self.max_bank;
        let turn_rate = self.gravity * bank.tan() / self.airspeed;
        delta_psi.abs() / turn_rate.max(1e-6)
    }

    /// Energy-altitude trade: specific excess power (m/s).
    ///
    /// SEP = (thrust − drag) * V / (mass * g).
    pub fn specific_excess_power(&self, thrust: f64, drag: f64, mass: f64) -> f64 {
        (thrust - drag) * self.airspeed / (mass * self.gravity)
    }

    /// Zoom climb: altitude gain from converting kinetic energy at constant g.
    ///
    /// Δh = (V₂² − V₁²) / (2 * g).
    pub fn zoom_climb_altitude(&self, v_initial: f64, v_final: f64) -> f64 {
        (v_initial * v_initial - v_final * v_final) / (2.0 * self.gravity)
    }

    /// Minimum radius for a level banked turn at given load factor n.
    pub fn min_turn_radius(&self, load_factor: f64) -> f64 {
        let n = load_factor.max(1.0);
        self.airspeed * self.airspeed / (self.gravity * (n * n - 1.0).max(0.0).sqrt())
    }

    /// Corner velocity: airspeed for tightest structural and aerodynamic turn.
    ///
    /// V* = sqrt(2 * n_max * W / (ρ * S * CLmax)).
    pub fn corner_velocity(
        &self,
        n_max: f64,
        weight: f64,
        rho: f64,
        s_ref: f64,
        cl_max: f64,
    ) -> f64 {
        (2.0 * n_max * weight / (rho * s_ref * cl_max)).sqrt()
    }

    /// Time-to-climb from h1 to h2 assuming constant SEP.
    pub fn time_to_climb(&self, h1: f64, h2: f64, sep: f64) -> f64 {
        if sep < 1e-6 {
            f64::INFINITY
        } else {
            (h2 - h1) / sep
        }
    }
}

// ── RecoverySystem ────────────────────────────────────────────────────────────

/// State of the recovery/parachute system.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecoveryState {
    /// Normal flight; no deployment.
    Normal,
    /// Drogue chute deployed.
    DrogueDeployed,
    /// Main chute fully open.
    MainOpen,
    /// Landed.
    Landed,
}

/// Parachute deployment and terminal velocity dynamics.
#[derive(Debug, Clone)]
pub struct RecoverySystem {
    /// System state.
    pub state: RecoveryState,
    /// Payload mass (kg).
    pub mass: f64,
    /// Drogue chute drag area (m²).
    pub cd_drogue: f64,
    /// Main chute drag area (m²).
    pub cd_main: f64,
    /// Effective canopy area in current state (m²).
    pub effective_cd_area: f64,
    /// Current vertical velocity (m/s, positive downward).
    pub velocity: f64,
    /// Altitude (m).
    pub altitude: f64,
    /// Air density at deployment (kg/m³).
    pub rho: f64,
    /// Time since deployment (s).
    pub time_since_deploy: f64,
    /// Inflation time constant (s).
    pub inflation_tau: f64,
    /// Landing impact velocity threshold (m/s).
    pub landing_threshold: f64,
}

impl RecoverySystem {
    /// Create a new recovery system.
    pub fn new(mass: f64, cd_drogue: f64, cd_main: f64) -> Self {
        Self {
            state: RecoveryState::Normal,
            mass,
            cd_drogue,
            cd_main,
            effective_cd_area: 0.0,
            velocity: 0.0,
            altitude: 3000.0,
            rho: 1.1,
            time_since_deploy: 0.0,
            inflation_tau: 1.5,
            landing_threshold: 7.0,
        }
    }

    /// Deploy drogue chute.
    pub fn deploy_drogue(&mut self) {
        self.state = RecoveryState::DrogueDeployed;
        self.time_since_deploy = 0.0;
    }

    /// Deploy main chute (from drogue or directly).
    pub fn deploy_main(&mut self) {
        self.state = RecoveryState::MainOpen;
        self.time_since_deploy = 0.0;
    }

    /// Terminal velocity for fully inflated canopy.
    pub fn terminal_velocity(&self) -> f64 {
        let cda = match self.state {
            RecoveryState::DrogueDeployed => self.cd_drogue,
            RecoveryState::MainOpen => self.cd_main,
            _ => 1e-6,
        };
        ((2.0 * self.mass * 9.80665) / (self.rho * cda)).sqrt()
    }

    /// Advance recovery simulation by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        if self.state == RecoveryState::Landed || self.state == RecoveryState::Normal {
            return;
        }
        self.time_since_deploy += dt;

        // Inflate canopy with exponential approach
        let cd_target = match self.state {
            RecoveryState::DrogueDeployed => self.cd_drogue,
            RecoveryState::MainOpen => self.cd_main,
            _ => 0.0,
        };
        let fill = 1.0 - (-self.time_since_deploy / self.inflation_tau).exp();
        self.effective_cd_area = cd_target * fill;

        // Equations of motion (1-D, positive downward)
        let drag = 0.5 * self.rho * self.effective_cd_area * self.velocity * self.velocity;
        let net_force = self.mass * 9.80665 - drag;
        let accel = net_force / self.mass;
        self.velocity = (self.velocity + accel * dt).max(0.0);
        self.altitude -= self.velocity * dt;

        if self.altitude <= 0.0 {
            self.altitude = 0.0;
            self.state = RecoveryState::Landed;
        }
    }

    /// Landing impact kinetic energy (J).
    pub fn impact_energy(&self) -> f64 {
        0.5 * self.mass * self.velocity * self.velocity
    }

    /// Whether landing is within safe velocity threshold.
    pub fn is_safe_landing(&self) -> bool {
        self.velocity <= self.landing_threshold
    }
}

// ── Atmosphere and Wind Model ─────────────────────────────────────────────────

/// Simple turbulence / wind shear model.
#[derive(Debug, Clone)]
pub struct WindModel {
    /// Mean wind vector in NED (m/s).
    pub mean_wind: [f64; 3],
    /// Turbulence intensity (m/s RMS).
    pub turbulence_intensity: f64,
    /// Wind shear: change in horizontal wind per meter altitude (1/s).
    pub shear_gradient: f64,
}

impl WindModel {
    /// Create a wind model.
    pub fn new(mean_wind: [f64; 3], turbulence_intensity: f64, shear_gradient: f64) -> Self {
        Self {
            mean_wind,
            turbulence_intensity,
            shear_gradient,
        }
    }

    /// Wind velocity at altitude `h` (NED).
    pub fn wind_at_altitude(&self, h: f64) -> [f64; 3] {
        let shear = self.shear_gradient * h;
        [
            self.mean_wind[0] + shear,
            self.mean_wind[1],
            self.mean_wind[2],
        ]
    }

    /// Gust intensity factor (dimensionless) at given altitude.
    pub fn gust_factor(&self, h: f64) -> f64 {
        // Logarithmic wind profile
        if h < 1.0 {
            0.5
        } else {
            (h / 10.0).ln().max(0.1) * 0.3 + 1.0
        }
    }
}

// ── MachNumber utilities ──────────────────────────────────────────────────────

/// Compute Mach number from true airspeed and altitude.
pub fn mach_number(airspeed: f64, altitude: f64) -> f64 {
    let a = isa_speed_of_sound(altitude);
    if a < 1.0 { airspeed } else { airspeed / a }
}

/// Stagnation (total) pressure ratio for subsonic Mach number.
pub fn total_pressure_ratio(mach: f64) -> f64 {
    (1.0 + 0.2 * mach * mach).powf(3.5)
}

/// Prandtl–Glauert compressibility correction to lift coefficient.
pub fn prandtl_glauert(cl: f64, mach: f64) -> f64 {
    let beta = (1.0 - mach * mach).max(0.01).sqrt();
    cl / beta
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-8;

    // ── AircraftDynamics ──────────────────────────────────────────────────

    #[test]
    fn test_aircraft_new() {
        let ac = AircraftDynamics::new(8000.0, [4000.0, 3000.0, 5000.0]);
        assert!((ac.mass - 8000.0).abs() < EPS);
        assert_eq!(ac.inertia, [4000.0, 3000.0, 5000.0]);
    }

    #[test]
    fn test_aircraft_true_airspeed_zero() {
        let ac = AircraftDynamics::new(1000.0, [100.0, 100.0, 100.0]);
        assert!(ac.true_airspeed().abs() < EPS);
    }

    #[test]
    fn test_aircraft_true_airspeed_nonzero() {
        let mut ac = AircraftDynamics::new(1000.0, [100.0, 100.0, 100.0]);
        ac.velocity_body = [100.0, 0.0, 0.0];
        assert!((ac.true_airspeed() - 100.0).abs() < EPS);
    }

    #[test]
    fn test_aircraft_angle_of_attack() {
        let mut ac = AircraftDynamics::new(1000.0, [100.0, 100.0, 100.0]);
        ac.velocity_body = [100.0, 0.0, 0.0];
        assert!(ac.angle_of_attack().abs() < EPS);
    }

    #[test]
    fn test_aircraft_sideslip_zero() {
        let mut ac = AircraftDynamics::new(1000.0, [100.0, 100.0, 100.0]);
        ac.velocity_body = [100.0, 0.0, 0.0];
        assert!(ac.sideslip().abs() < EPS);
    }

    #[test]
    fn test_aircraft_altitude() {
        let mut ac = AircraftDynamics::new(1000.0, [100.0, 100.0, 100.0]);
        ac.position = [0.0, 0.0, -1000.0]; // NED: down = positive, so h=1000
        assert!((ac.altitude() - 1000.0).abs() < EPS);
    }

    #[test]
    fn test_aircraft_kinetic_energy() {
        let mut ac = AircraftDynamics::new(2.0, [1.0, 1.0, 1.0]);
        ac.velocity_body = [3.0, 0.0, 0.0];
        // KE_linear = 0.5 * 2 * 9 = 9
        assert!((ac.kinetic_energy() - 9.0).abs() < 1e-6);
    }

    #[test]
    fn test_aircraft_step_gravity() {
        let mut ac = AircraftDynamics::new(100.0, [100.0, 100.0, 100.0]);
        ac.add_gravity_force();
        ac.step(0.1);
        // With gravity, position should change
        assert!(ac.position[2] != 0.0 || ac.velocity_body[2] != 0.0);
    }

    #[test]
    fn test_aircraft_euler_angles_identity() {
        let ac = AircraftDynamics::new(1000.0, [100.0, 100.0, 100.0]);
        let [phi, theta, psi] = ac.euler_angles();
        assert!(phi.abs() < EPS);
        assert!(theta.abs() < EPS);
        assert!(psi.abs() < EPS);
    }

    #[test]
    fn test_aircraft_load_factor() {
        let ac = AircraftDynamics::new(1000.0, [100.0, 100.0, 100.0]);
        let weight = 1000.0 * 9.80665;
        let n = ac.load_factor(weight);
        assert!((n - 1.0).abs() < 1e-4);
    }

    #[test]
    fn test_aircraft_clear_forces() {
        let mut ac = AircraftDynamics::new(1000.0, [100.0, 100.0, 100.0]);
        ac.apply_force_body([100.0, 200.0, 300.0]);
        ac.clear_forces();
        assert!(norm3(ac.force_body) < EPS);
        assert!(norm3(ac.torque_body) < EPS);
    }

    #[test]
    fn test_aircraft_total_energy() {
        let mut ac = AircraftDynamics::new(1.0, [1.0, 1.0, 1.0]);
        ac.velocity_body = [10.0, 0.0, 0.0];
        ac.position = [0.0, 0.0, -100.0];
        let ke = ac.kinetic_energy();
        let pe = ac.potential_energy();
        let total = ac.total_energy();
        assert!((total - ke - pe).abs() < 1e-8);
    }

    // ── FlightControlSurface ──────────────────────────────────────────────

    #[test]
    fn test_surface_new_zero_deflection() {
        let s = FlightControlSurface::new(SurfaceKind::Elevator, 0.35, 1.0, 0.1, 5000.0);
        assert!(s.deflection.abs() < EPS);
        assert_eq!(s.kind, SurfaceKind::Elevator);
    }

    #[test]
    fn test_surface_ideal_actuator() {
        let mut s = FlightControlSurface::new(SurfaceKind::Aileron, 0.4, 1.0, 0.0, 3000.0);
        s.set_command(0.2);
        s.step(0.1);
        // Rate = 1.0 rad/s, dt=0.1 → can move 0.1 rad; command is 0.2 → moves 0.1
        assert!((s.deflection - 0.1).abs() < 1e-6);
    }

    #[test]
    fn test_surface_lag_actuator() {
        let mut s = FlightControlSurface::new(SurfaceKind::Rudder, 0.35, 10.0, 0.5, 4000.0);
        s.set_command(0.3);
        for _ in 0..1000 {
            s.step(0.01);
        }
        // Should approach command after many steps (10s >> tau=0.5s)
        assert!((s.deflection - 0.3).abs() < 0.01);
    }

    #[test]
    fn test_surface_clamp_command() {
        let mut s = FlightControlSurface::new(SurfaceKind::Flap, 0.5, 1.0, 0.0, 1000.0);
        s.set_command(2.0); // beyond limit
        assert!((s.command - 0.5).abs() < EPS);
    }

    #[test]
    fn test_surface_moment() {
        let mut s = FlightControlSurface::new(SurfaceKind::Elevator, 0.5, 5.0, 0.0, 10_000.0);
        s.deflection = 0.1;
        assert!((s.moment() - 1000.0).abs() < 1e-6);
    }

    #[test]
    fn test_surface_delta_cl() {
        let mut s = FlightControlSurface::new(SurfaceKind::Flap, 0.5, 1.0, 0.0, 1000.0);
        s.deflection = 0.2;
        let dcl = s.delta_cl(3.0); // 3.0 per rad * 0.2 rad = 0.6
        assert!((dcl - 0.6).abs() < EPS);
    }

    // ── InertialNavigation ────────────────────────────────────────────────

    #[test]
    fn test_ins_new() {
        let ins = InertialNavigation::new([0.0, 0.0, -100.0], [0.0, 0.0, 0.0, 1.0]);
        assert!((ins.altitude() - 100.0).abs() < EPS);
    }

    #[test]
    fn test_ins_altitude() {
        let ins = InertialNavigation::new([0.0, 0.0, -500.0], [0.0, 0.0, 0.0, 1.0]);
        assert!((ins.altitude() - 500.0).abs() < EPS);
    }

    #[test]
    fn test_ins_imu_update_position() {
        let mut ins = InertialNavigation::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]);
        // Accelerometer reads gravity in NED frame (z-down): +9.80665 on z axis
        ins.update_imu([0.0, 0.0, 9.80665], [0.0, 0.0, 0.0], 0.1);
        // After gravity subtraction, net accel ≈ 0 → position shouldn't change much
        let dist = norm3(ins.position);
        assert!(
            dist < 0.1,
            "Position should not drift much under gravity compensation: {dist}"
        );
    }

    #[test]
    fn test_ins_euler_angles_identity() {
        let ins = InertialNavigation::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]);
        let [phi, theta, psi] = ins.euler_angles();
        assert!(phi.abs() < EPS);
        assert!(theta.abs() < EPS);
        assert!(psi.abs() < EPS);
    }

    #[test]
    fn test_ins_gps_correction() {
        let mut ins = InertialNavigation::new([100.0, 200.0, -300.0], [0.0, 0.0, 0.0, 1.0]);
        ins.gps_correction([0.0, 0.0, 0.0], 1.0); // full correction to origin
        assert!(norm3(ins.position) < EPS);
    }

    #[test]
    fn test_ins_horizontal_speed() {
        let mut ins = InertialNavigation::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]);
        ins.velocity = [3.0, 4.0, 0.0];
        assert!((ins.horizontal_speed() - 5.0).abs() < EPS);
    }

    #[test]
    fn test_ins_complementary_correction() {
        let mut ins = InertialNavigation::new([0.0, 0.0, 0.0], [0.1, 0.0, 0.0, 0.99499]);
        let ref_q = [0.0, 0.0, 0.0, 1.0];
        for _ in 0..100 {
            ins.complementary_correction(ref_q);
        }
        // After many corrections, attitude should approach reference
        assert!(ins.attitude[0].abs() < 0.05);
    }

    // ── AerodynamicForces ─────────────────────────────────────────────────

    #[test]
    fn test_aero_dynamic_pressure() {
        let aero = AerodynamicForces::new(5.73, 0.025, -0.5, 1.225);
        let q = aero.dynamic_pressure(100.0);
        assert!((q - 0.5 * 1.225 * 10_000.0).abs() < EPS);
    }

    #[test]
    fn test_aero_lift_positive() {
        let aero = AerodynamicForces::new(5.73, 0.025, -0.5, 1.225);
        let l = aero.lift(100.0, 0.1);
        assert!(l > 0.0);
    }

    #[test]
    fn test_aero_drag_positive() {
        let aero = AerodynamicForces::new(5.73, 0.025, -0.5, 1.225);
        let d = aero.drag(100.0, 0.0);
        assert!(d > 0.0);
    }

    #[test]
    fn test_aero_cl_increases_with_alpha() {
        let aero = AerodynamicForces::new(5.73, 0.025, -0.5, 1.225);
        let cl1 = aero.cl(0.0);
        let cl2 = aero.cl(0.2);
        assert!(cl2 > cl1);
    }

    #[test]
    fn test_aero_cd_parabolic_polar() {
        let aero = AerodynamicForces::new(5.73, 0.025, -0.5, 1.225);
        let cd_zero = aero.cd(0.0);
        assert!(cd_zero >= aero.cd0);
    }

    #[test]
    fn test_aero_lift_to_drag_ratio() {
        let aero = AerodynamicForces::new(5.73, 0.025, -0.5, 1.225);
        let ld = aero.lift_to_drag(0.0);
        assert!(ld > 0.0);
    }

    #[test]
    fn test_aero_alpha_best_ld() {
        let aero = AerodynamicForces::new(5.73, 0.025, -0.5, 1.225);
        let alpha_star = aero.alpha_best_ld();
        // Should be a finite real number
        assert!(alpha_star.is_finite());
    }

    #[test]
    fn test_aero_stall_alpha() {
        let aero = AerodynamicForces::new(5.73, 0.025, -0.5, 1.225);
        let alpha_stall = aero.stall_alpha();
        assert!(alpha_stall > 0.0);
    }

    #[test]
    fn test_aero_pitching_moment() {
        let aero = AerodynamicForces::new(5.73, 0.025, -0.5, 1.225);
        let m = aero.pitching_moment(100.0, 0.1);
        assert!(m.is_finite());
    }

    #[test]
    fn test_aero_to_body_forces() {
        let aero = AerodynamicForces::new(5.73, 0.025, -0.5, 1.225);
        let [fx, fy, fz] = aero.to_body_forces(100.0, 0.05, 0.0);
        assert!(fx.is_finite() && fy.is_finite() && fz.is_finite());
    }

    // ── ManeuverPlanning ──────────────────────────────────────────────────

    #[test]
    fn test_maneuver_min_time_turn() {
        let plan = ManeuverPlanning::new(500.0, 60.0_f64.to_radians(), 100.0);
        let t = plan.min_time_turn(PI / 2.0);
        assert!(t > 0.0);
        assert!(t.is_finite());
    }

    #[test]
    fn test_maneuver_sep() {
        let plan = ManeuverPlanning::new(500.0, 60.0_f64.to_radians(), 100.0);
        let sep = plan.specific_excess_power(50_000.0, 20_000.0, 8000.0);
        assert!((sep - 30_000.0 * 100.0 / (8000.0 * 9.80665)).abs() < 1e-4);
    }

    #[test]
    fn test_maneuver_zoom_climb() {
        let plan = ManeuverPlanning::new(500.0, 60.0_f64.to_radians(), 100.0);
        let dh = plan.zoom_climb_altitude(200.0, 0.0);
        let expected = 200.0 * 200.0 / (2.0 * 9.80665);
        assert!((dh - expected).abs() < 1e-4);
    }

    #[test]
    fn test_maneuver_turn_radius() {
        let plan = ManeuverPlanning::new(500.0, 60.0_f64.to_radians(), 100.0);
        let r = plan.min_turn_radius(2.0);
        assert!(r > 0.0);
        assert!(r.is_finite());
    }

    #[test]
    fn test_maneuver_corner_velocity() {
        let plan = ManeuverPlanning::new(500.0, 60.0_f64.to_radians(), 100.0);
        let v = plan.corner_velocity(6.0, 80_000.0, 1.225, 30.0, 1.6);
        assert!(v > 0.0);
    }

    #[test]
    fn test_maneuver_time_to_climb() {
        let plan = ManeuverPlanning::new(500.0, 60.0_f64.to_radians(), 100.0);
        let t = plan.time_to_climb(1000.0, 5000.0, 10.0);
        assert!((t - 400.0).abs() < EPS);
    }

    #[test]
    fn test_maneuver_time_to_climb_zero_sep() {
        let plan = ManeuverPlanning::new(500.0, 60.0_f64.to_radians(), 100.0);
        let t = plan.time_to_climb(0.0, 1000.0, 0.0);
        assert!(t.is_infinite());
    }

    #[test]
    fn test_dubins_rsr_path() {
        let plan = ManeuverPlanning::new(100.0, 60.0_f64.to_radians(), 50.0);
        let path = plan.dubins_rsr([0.0, 0.0, 0.0], [300.0, 0.0, 0.0]);
        assert!(path.total_length() > 0.0);
        assert_eq!(path.segments[1], DubinsSegment::Straight);
    }

    // ── RecoverySystem ────────────────────────────────────────────────────

    #[test]
    fn test_recovery_new() {
        let r = RecoverySystem::new(100.0, 5.0, 50.0);
        assert_eq!(r.state, RecoveryState::Normal);
        assert!((r.mass - 100.0).abs() < EPS);
    }

    #[test]
    fn test_recovery_terminal_velocity_main() {
        let mut r = RecoverySystem::new(100.0, 5.0, 50.0);
        r.deploy_main();
        // Fully inflate
        r.time_since_deploy = 100.0;
        r.effective_cd_area = r.cd_main;
        let vt = r.terminal_velocity();
        let expected = (2.0_f64 * 100.0 * 9.80665 / (1.1 * 50.0)).sqrt();
        assert!((vt - expected).abs() < 0.01);
    }

    #[test]
    fn test_recovery_deploy_drogue() {
        let mut r = RecoverySystem::new(100.0, 5.0, 50.0);
        r.deploy_drogue();
        assert_eq!(r.state, RecoveryState::DrogueDeployed);
    }

    #[test]
    fn test_recovery_deploy_main() {
        let mut r = RecoverySystem::new(100.0, 5.0, 50.0);
        r.deploy_main();
        assert_eq!(r.state, RecoveryState::MainOpen);
    }

    #[test]
    fn test_recovery_step_slows_descent() {
        let mut r = RecoverySystem::new(100.0, 5.0, 50.0);
        r.deploy_main();
        r.velocity = 50.0; // fast initial descent
        r.altitude = 1000.0;
        for _ in 0..1000 {
            r.step(0.1);
            if r.state == RecoveryState::Landed {
                break;
            }
        }
        // Either landed or velocity reduced
        assert!(r.velocity < 50.0 || r.state == RecoveryState::Landed);
    }

    #[test]
    fn test_recovery_impact_energy() {
        let mut r = RecoverySystem::new(50.0, 5.0, 50.0);
        r.velocity = 6.0;
        let ke = r.impact_energy();
        assert!((ke - 0.5 * 50.0 * 36.0).abs() < EPS);
    }

    #[test]
    fn test_recovery_safe_landing() {
        let mut r = RecoverySystem::new(100.0, 5.0, 50.0);
        r.velocity = 5.0;
        assert!(r.is_safe_landing());
        r.velocity = 10.0;
        assert!(!r.is_safe_landing());
    }

    // ── Atmosphere / Mach utilities ───────────────────────────────────────

    #[test]
    fn test_isa_density_sea_level() {
        let rho = isa_density(0.0);
        assert!((rho - 1.225).abs() < EPS);
    }

    #[test]
    fn test_isa_density_decreases_with_altitude() {
        assert!(isa_density(10_000.0) < isa_density(0.0));
    }

    #[test]
    fn test_isa_speed_of_sound_sea_level() {
        let a = isa_speed_of_sound(0.0);
        // Expected ≈ 340.3 m/s
        assert!((a - 340.3).abs() < 0.5, "speed of sound = {a}");
    }

    #[test]
    fn test_mach_number() {
        let m = mach_number(340.3, 0.0);
        assert!((m - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_total_pressure_ratio_sea_level() {
        let pr = total_pressure_ratio(0.0);
        assert!((pr - 1.0).abs() < EPS);
    }

    #[test]
    fn test_prandtl_glauert() {
        let cl_corrected = prandtl_glauert(1.0, 0.0);
        assert!((cl_corrected - 1.0).abs() < EPS);
        let cl_sub = prandtl_glauert(1.0, 0.5);
        assert!(cl_sub > 1.0); // compressibility increases effective CL
    }

    #[test]
    fn test_wind_model_at_altitude() {
        let wm = WindModel::new([10.0, 0.0, 0.0], 1.0, 0.01);
        let w = wm.wind_at_altitude(100.0);
        assert!((w[0] - 11.0).abs() < EPS);
    }

    #[test]
    fn test_wind_model_gust_factor() {
        let wm = WindModel::new([0.0; 3], 0.5, 0.0);
        let g1 = wm.gust_factor(1.0);
        let g100 = wm.gust_factor(100.0);
        assert!(g100 > g1);
    }
}
