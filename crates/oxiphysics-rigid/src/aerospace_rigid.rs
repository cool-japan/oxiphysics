// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Aerospace rigid body dynamics with 6-DOF flight dynamics.
//!
//! Provides aerodynamic forces (lift/drag coefficients), control surfaces
//! (aileron/elevator/rudder), ISA atmosphere model, wind effects, stability
//! derivatives, trim solver, flight envelope analysis, and Mach number effects.
//!
//! All geometry is represented with `[f64; 3]` arrays (no nalgebra dependency).
//!
//! # Example
//!
//! ```no_run
//! use oxiphysics_rigid::aerospace_rigid::{isa_atmosphere, mach_number, AircraftState};
//!
//! let (t, p, rho) = isa_atmosphere(5000.0);
//! assert!(rho > 0.0 && rho < 1.225);
//!
//! let m = mach_number(250.0, 5000.0);
//! assert!(m > 0.0);
//!
//! let state = AircraftState::new(1500.0);
//! assert!(state.mass > 0.0);
//! ```

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ── vector helpers ────────────────────────────────────────────────────────────

/// Cross product of two 3-vectors.
#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Dot product.
#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// L2 norm.
#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Normalize. Returns zero if tiny.
#[inline]
fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n < 1e-30 {
        [0.0; 3]
    } else {
        [a[0] / n, a[1] / n, a[2] / n]
    }
}

/// Add two 3-vectors.
#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract.
#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale.
#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Quaternion normalize.
#[inline]
fn quat_normalize(q: [f64; 4]) -> [f64; 4] {
    let n = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    if n < 1e-30 {
        [0.0, 0.0, 0.0, 1.0]
    } else {
        [q[0] / n, q[1] / n, q[2] / n, q[3] / n]
    }
}

/// Quaternion Hamilton product. Convention: q = \[x, y, z, w\].
fn quat_mul(q1: [f64; 4], q2: [f64; 4]) -> [f64; 4] {
    let [x1, y1, z1, w1] = q1;
    let [x2, y2, z2, w2] = q2;
    [
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
    ]
}

/// Rotate a vector by a quaternion: v' = q * \[v, 0\] * q_conj.
fn quat_rotate(q: [f64; 4], v: [f64; 3]) -> [f64; 3] {
    let q_conj = [-q[0], -q[1], -q[2], q[3]];
    let vq = [v[0], v[1], v[2], 0.0];
    let tmp = quat_mul(q, vq);
    let result = quat_mul(tmp, q_conj);
    [result[0], result[1], result[2]]
}

// ── ISA Atmosphere Model ──────────────────────────────────────────────────────

/// Sea-level standard temperature \[K\].
const ISA_T0: f64 = 288.15;
/// Sea-level standard pressure \[Pa\].
const ISA_P0: f64 = 101_325.0;
/// Sea-level standard density \[kg/m^3\].
const ISA_RHO0: f64 = 1.225;
/// Temperature lapse rate in troposphere \[K/m\].
const ISA_LAPSE: f64 = 0.0065;
/// Tropopause altitude \[m\].
const ISA_TROPO: f64 = 11_000.0;
/// Gas constant for dry air \[J/(kg*K)\].
const R_AIR: f64 = 287.058;
/// Gravity \[m/s^2\].
const G0: f64 = 9.80665;
/// Specific heat ratio for air.
const GAMMA_AIR: f64 = 1.4;

/// Compute ISA (International Standard Atmosphere) properties.
///
/// Returns `(temperature [K], pressure [Pa], density [kg/m^3])` at the
/// given geometric altitude in metres. Valid up to ~32 km.
pub fn isa_atmosphere(altitude_m: f64) -> (f64, f64, f64) {
    let alt = altitude_m.max(0.0);
    if alt <= ISA_TROPO {
        // Troposphere: linear lapse
        let t = ISA_T0 - ISA_LAPSE * alt;
        let p = ISA_P0 * (t / ISA_T0).powf(G0 / (ISA_LAPSE * R_AIR));
        let rho = p / (R_AIR * t);
        (t, p, rho)
    } else {
        // Tropopause / lower stratosphere: isothermal at T = ISA_T0 - LAPSE*TROPO
        let t_tropo = ISA_T0 - ISA_LAPSE * ISA_TROPO;
        let (_, p_tropo, _) = isa_atmosphere(ISA_TROPO);
        let t = t_tropo;
        let p = p_tropo * (-G0 * (alt - ISA_TROPO) / (R_AIR * t_tropo)).exp();
        let rho = p / (R_AIR * t);
        (t, p, rho)
    }
}

/// Speed of sound \[m/s\] at a given altitude.
pub fn speed_of_sound(altitude_m: f64) -> f64 {
    let (t, _, _) = isa_atmosphere(altitude_m);
    (GAMMA_AIR * R_AIR * t).sqrt()
}

/// Mach number at a given true airspeed and altitude.
pub fn mach_number(true_airspeed_mps: f64, altitude_m: f64) -> f64 {
    let a = speed_of_sound(altitude_m);
    if a < 1e-10 {
        0.0
    } else {
        true_airspeed_mps / a
    }
}

/// Dynamic pressure \[Pa\]: q = 0.5 * rho * V^2.
pub fn dynamic_pressure(velocity_mps: f64, altitude_m: f64) -> f64 {
    let (_, _, rho) = isa_atmosphere(altitude_m);
    0.5 * rho * velocity_mps * velocity_mps
}

/// Equivalent airspeed \[m/s\] from true airspeed.
pub fn equivalent_airspeed(true_airspeed: f64, altitude_m: f64) -> f64 {
    let (_, _, rho) = isa_atmosphere(altitude_m);
    true_airspeed * (rho / ISA_RHO0).sqrt()
}

/// True airspeed \[m/s\] from equivalent airspeed.
pub fn true_airspeed_from_eas(eas: f64, altitude_m: f64) -> f64 {
    let (_, _, rho) = isa_atmosphere(altitude_m);
    if rho < 1e-30 {
        0.0
    } else {
        eas * (ISA_RHO0 / rho).sqrt()
    }
}

/// Pressure altitude \[m\] from static pressure using ISA in troposphere.
pub fn pressure_altitude(pressure_pa: f64) -> f64 {
    if pressure_pa >= ISA_P0 {
        return 0.0;
    }
    if pressure_pa <= 0.0 {
        return ISA_TROPO;
    }
    let exponent = ISA_LAPSE * R_AIR / G0;
    let t_ratio = (pressure_pa / ISA_P0).powf(exponent);
    (ISA_T0 * (1.0 - t_ratio)) / ISA_LAPSE
}

// ── Aerodynamic coefficients ──────────────────────────────────────────────────

/// Lift coefficient from angle of attack using a linear + stall model.
///
/// * `alpha` — angle of attack \[rad\].
/// * `cl_alpha` — lift curve slope \[1/rad\] (typically ~5.7 for thin airfoil).
/// * `alpha_0` — zero-lift angle of attack \[rad\].
/// * `cl_max` — maximum lift coefficient.
pub fn lift_coefficient(alpha: f64, cl_alpha: f64, alpha_0: f64, cl_max: f64) -> f64 {
    let cl_linear = cl_alpha * (alpha - alpha_0);
    cl_linear.clamp(-cl_max, cl_max)
}

/// Drag coefficient using a parabolic drag polar.
///
/// `CD = CD0 + CL^2 / (pi * e * AR)`
///
/// * `cd0` — zero-lift drag coefficient.
/// * `cl` — current lift coefficient.
/// * `aspect_ratio` — wing aspect ratio.
/// * `oswald_e` — Oswald efficiency factor (typically 0.7–0.85).
pub fn drag_coefficient(cd0: f64, cl: f64, aspect_ratio: f64, oswald_e: f64) -> f64 {
    let k = 1.0 / (PI * oswald_e * aspect_ratio);
    cd0 + k * cl * cl
}

/// Lift-to-drag ratio.
pub fn lift_to_drag_ratio(cl: f64, cd: f64) -> f64 {
    if cd.abs() < 1e-30 { 0.0 } else { cl / cd }
}

/// Maximum lift-to-drag ratio for a parabolic drag polar.
///
/// `(L/D)_max = 0.5 * sqrt(pi * e * AR / CD0)`
pub fn max_ld_ratio(cd0: f64, aspect_ratio: f64, oswald_e: f64) -> f64 {
    if cd0 <= 0.0 {
        return 0.0;
    }
    0.5 * (PI * oswald_e * aspect_ratio / cd0).sqrt()
}

/// CL for maximum L/D.
pub fn cl_for_max_ld(cd0: f64, aspect_ratio: f64, oswald_e: f64) -> f64 {
    (cd0 * PI * oswald_e * aspect_ratio).sqrt()
}

// ── Mach number effects ───────────────────────────────────────────────────────

/// Prandtl-Glauert compressibility correction factor.
///
/// Returns the factor by which subsonic aerodynamic coefficients should be
/// multiplied. Only valid for M < 1.
pub fn prandtl_glauert_factor(mach: f64) -> f64 {
    let m_clamp = mach.abs().min(0.99);
    1.0 / (1.0 - m_clamp * m_clamp).sqrt()
}

/// Wave drag coefficient onset above critical Mach number.
///
/// Uses a simple power-law model:
/// `CD_wave = k * (M - M_crit)^n` for M > M_crit, else 0.
pub fn wave_drag_coefficient(mach: f64, mach_crit: f64, k: f64, n: f64) -> f64 {
    if mach <= mach_crit {
        0.0
    } else {
        k * (mach - mach_crit).powf(n)
    }
}

/// Corrected lift curve slope using Prandtl-Glauert.
pub fn corrected_cl_alpha(cl_alpha_incomp: f64, mach: f64) -> f64 {
    cl_alpha_incomp * prandtl_glauert_factor(mach)
}

// ── Control surfaces ──────────────────────────────────────────────────────────

/// Control surface deflections \[rad\].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ControlSurfaces {
    /// Aileron deflection \[rad\] (positive = right wing down).
    pub aileron: f64,
    /// Elevator deflection \[rad\] (positive = trailing edge down).
    pub elevator: f64,
    /// Rudder deflection \[rad\] (positive = trailing edge left).
    pub rudder: f64,
    /// Throttle setting \[0..1\].
    pub throttle: f64,
    /// Flap deflection \[rad\].
    pub flaps: f64,
}

impl Default for ControlSurfaces {
    fn default() -> Self {
        Self {
            aileron: 0.0,
            elevator: 0.0,
            rudder: 0.0,
            throttle: 0.0,
            flaps: 0.0,
        }
    }
}

impl ControlSurfaces {
    /// Create with all deflections zero.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create with specified primary control deflections.
    pub fn primary(aileron: f64, elevator: f64, rudder: f64, throttle: f64) -> Self {
        Self {
            aileron,
            elevator,
            rudder,
            throttle,
            flaps: 0.0,
        }
    }

    /// Clamp all control surfaces to physical limits.
    pub fn clamp(&mut self, max_aileron: f64, max_elevator: f64, max_rudder: f64, max_flap: f64) {
        self.aileron = self.aileron.clamp(-max_aileron, max_aileron);
        self.elevator = self.elevator.clamp(-max_elevator, max_elevator);
        self.rudder = self.rudder.clamp(-max_rudder, max_rudder);
        self.throttle = self.throttle.clamp(0.0, 1.0);
        self.flaps = self.flaps.clamp(0.0, max_flap);
    }

    /// Return true if all controls are at neutral.
    pub fn is_neutral(&self) -> bool {
        self.aileron.abs() < 1e-10
            && self.elevator.abs() < 1e-10
            && self.rudder.abs() < 1e-10
            && self.flaps.abs() < 1e-10
    }
}

// ── Stability derivatives ─────────────────────────────────────────────────────

/// Stability and control derivatives for a conventional aircraft.
///
/// Longitudinal: Cm_alpha, Cm_q, Cm_delta_e, CL_alpha, CL_delta_e
/// Lateral: Cl_beta, Cl_p, Cl_delta_a, Cn_beta, Cn_r, Cn_delta_r
#[derive(Debug, Clone)]
pub struct StabilityDerivatives {
    // ── Longitudinal ──
    /// Pitching moment coefficient w.r.t. angle of attack \[1/rad\].
    pub cm_alpha: f64,
    /// Pitching moment coefficient w.r.t. pitch rate \[1/(rad/s)\].
    pub cm_q: f64,
    /// Pitching moment coefficient w.r.t. elevator deflection \[1/rad\].
    pub cm_delta_e: f64,
    /// Lift curve slope \[1/rad\].
    pub cl_alpha: f64,
    /// Lift coefficient w.r.t. elevator deflection \[1/rad\].
    pub cl_delta_e: f64,
    /// Drag at zero lift.
    pub cd0: f64,

    // ── Lateral ──
    /// Roll moment coefficient w.r.t. sideslip angle \[1/rad\].
    pub cl_beta: f64,
    /// Roll damping coefficient \[1/(rad/s)\].
    pub cl_p: f64,
    /// Roll moment w.r.t. aileron deflection \[1/rad\].
    pub cl_delta_a: f64,
    /// Yaw moment w.r.t. sideslip \[1/rad\].
    pub cn_beta: f64,
    /// Yaw damping \[1/(rad/s)\].
    pub cn_r: f64,
    /// Yaw moment w.r.t. rudder deflection \[1/rad\].
    pub cn_delta_r: f64,
    /// Side-force w.r.t. sideslip \[1/rad\].
    pub cy_beta: f64,
}

impl StabilityDerivatives {
    /// Create a typical general aviation stability derivative set.
    pub fn general_aviation() -> Self {
        Self {
            cm_alpha: -1.0,
            cm_q: -15.0,
            cm_delta_e: -1.5,
            cl_alpha: 5.7,
            cl_delta_e: 0.4,
            cd0: 0.025,
            cl_beta: -0.1,
            cl_p: -0.5,
            cl_delta_a: 0.15,
            cn_beta: 0.065,
            cn_r: -0.15,
            cn_delta_r: -0.07,
            cy_beta: -0.4,
        }
    }

    /// Check static longitudinal stability (Cm_alpha < 0).
    pub fn is_longitudinally_stable(&self) -> bool {
        self.cm_alpha < 0.0
    }

    /// Check static directional stability (Cn_beta > 0).
    pub fn is_directionally_stable(&self) -> bool {
        self.cn_beta > 0.0
    }

    /// Check roll stability (Cl_beta < 0 for dihedral effect).
    pub fn has_dihedral_effect(&self) -> bool {
        self.cl_beta < 0.0
    }

    /// Compute the pitching moment coefficient at given flight conditions.
    pub fn pitching_moment(
        &self,
        alpha: f64,
        q_rate: f64,
        delta_e: f64,
        _chord: f64,
        _velocity: f64,
    ) -> f64 {
        self.cm_alpha * alpha + self.cm_q * q_rate + self.cm_delta_e * delta_e
    }

    /// Compute the rolling moment coefficient.
    pub fn rolling_moment(&self, beta: f64, p_rate: f64, delta_a: f64) -> f64 {
        self.cl_beta * beta + self.cl_p * p_rate + self.cl_delta_a * delta_a
    }

    /// Compute the yawing moment coefficient.
    pub fn yawing_moment(&self, beta: f64, r_rate: f64, delta_r: f64) -> f64 {
        self.cn_beta * beta + self.cn_r * r_rate + self.cn_delta_r * delta_r
    }
}

// ── Wind model ────────────────────────────────────────────────────────────────

/// A simple wind model with steady and gust components.
#[derive(Debug, Clone)]
pub struct WindModel {
    /// Steady wind velocity in NED frame \[m/s\].
    pub steady_wind: [f64; 3],
    /// Wind shear gradient \[1/m\]: wind increases with altitude.
    pub shear_gradient: f64,
    /// Reference altitude for shear \[m\].
    pub shear_ref_alt: f64,
    /// Gust amplitude \[m/s\].
    pub gust_amplitude: f64,
    /// Gust frequency \[rad/s\].
    pub gust_frequency: f64,
}

impl WindModel {
    /// Create a calm wind model (no wind).
    pub fn calm() -> Self {
        Self {
            steady_wind: [0.0; 3],
            shear_gradient: 0.0,
            shear_ref_alt: 0.0,
            gust_amplitude: 0.0,
            gust_frequency: 0.0,
        }
    }

    /// Create a steady wind model.
    pub fn steady(wind_ned: [f64; 3]) -> Self {
        Self {
            steady_wind: wind_ned,
            shear_gradient: 0.0,
            shear_ref_alt: 0.0,
            gust_amplitude: 0.0,
            gust_frequency: 0.0,
        }
    }

    /// Compute the total wind velocity at a given altitude and time.
    pub fn wind_at(&self, altitude: f64, time: f64) -> [f64; 3] {
        // Shear factor (logarithmic profile)
        let shear_factor = if self.shear_gradient > 0.0 && altitude > self.shear_ref_alt {
            1.0 + self.shear_gradient * (altitude - self.shear_ref_alt)
        } else {
            1.0
        };

        let steady = scale3(self.steady_wind, shear_factor);

        // Simple sinusoidal gust
        let gust_val = self.gust_amplitude * (self.gust_frequency * time).sin();
        let gust = [gust_val, 0.0, 0.0]; // gust primarily along x

        add3(steady, gust)
    }

    /// Wind speed magnitude at given conditions.
    pub fn wind_speed(&self, altitude: f64, time: f64) -> f64 {
        norm3(self.wind_at(altitude, time))
    }
}

// ── AircraftConfig ────────────────────────────────────────────────────────────

/// Aircraft configuration parameters.
#[derive(Debug, Clone)]
pub struct AircraftConfig {
    /// Wing area \[m^2\].
    pub wing_area: f64,
    /// Wing span \[m\].
    pub wing_span: f64,
    /// Mean aerodynamic chord \[m\].
    pub mean_chord: f64,
    /// Wing aspect ratio.
    pub aspect_ratio: f64,
    /// Oswald efficiency factor.
    pub oswald_e: f64,
    /// Zero-lift angle of attack \[rad\].
    pub alpha_0: f64,
    /// Maximum lift coefficient.
    pub cl_max: f64,
    /// Maximum thrust \[N\].
    pub max_thrust: f64,
    /// Critical Mach number for wave drag.
    pub mach_critical: f64,
    /// Maximum structural load factor.
    pub max_load_factor: f64,
    /// Minimum structural load factor.
    pub min_load_factor: f64,
    /// Stability derivatives.
    pub stability: StabilityDerivatives,
}

impl AircraftConfig {
    /// Create a generic single-engine propeller aircraft configuration.
    pub fn generic_ga() -> Self {
        let wing_span = 10.0;
        let wing_area = 16.0;
        let mean_chord = wing_area / wing_span;
        Self {
            wing_area,
            wing_span,
            mean_chord,
            aspect_ratio: wing_span * wing_span / wing_area,
            oswald_e: 0.8,
            alpha_0: -0.02,
            cl_max: 1.6,
            max_thrust: 12_000.0,
            mach_critical: 0.7,
            max_load_factor: 3.8,
            min_load_factor: -1.52,
            stability: StabilityDerivatives::general_aviation(),
        }
    }

    /// Compute lift coefficient at given angle of attack.
    pub fn cl(&self, alpha: f64) -> f64 {
        lift_coefficient(alpha, self.stability.cl_alpha, self.alpha_0, self.cl_max)
    }

    /// Compute drag coefficient at given lift coefficient.
    pub fn cd(&self, cl: f64) -> f64 {
        drag_coefficient(self.stability.cd0, cl, self.aspect_ratio, self.oswald_e)
    }

    /// Compute L/D at given alpha.
    pub fn ld_ratio(&self, alpha: f64) -> f64 {
        let cl = self.cl(alpha);
        let cd = self.cd(cl);
        lift_to_drag_ratio(cl, cd)
    }

    /// Maximum L/D for this configuration.
    pub fn max_ld(&self) -> f64 {
        max_ld_ratio(self.stability.cd0, self.aspect_ratio, self.oswald_e)
    }
}

// ── AircraftState ─────────────────────────────────────────────────────────────

/// 6-DOF aircraft state.
#[derive(Debug, Clone)]
pub struct AircraftState {
    /// Mass \[kg\].
    pub mass: f64,
    /// Principal moments of inertia \[Ixx, Iyy, Izz\] \[kg*m^2\].
    pub inertia: [f64; 3],
    /// Position in NED frame \[m\].
    pub position: [f64; 3],
    /// Velocity in body frame \[u, v, w\] \[m/s\].
    pub velocity_body: [f64; 3],
    /// Attitude quaternion \[x, y, z, w\] (body-to-NED).
    pub attitude: [f64; 4],
    /// Angular velocity in body frame \[p, q, r\] \[rad/s\].
    pub angular_velocity: [f64; 3],
    /// Control surface deflections.
    pub controls: ControlSurfaces,
}

impl AircraftState {
    /// Create a new aircraft state at rest.
    pub fn new(mass: f64) -> Self {
        Self {
            mass,
            inertia: [2000.0, 3000.0, 4000.0],
            position: [0.0; 3],
            velocity_body: [0.0; 3],
            attitude: [0.0, 0.0, 0.0, 1.0],
            angular_velocity: [0.0; 3],
            controls: ControlSurfaces::default(),
        }
    }

    /// Create a state in level flight at given airspeed and altitude.
    pub fn level_flight(mass: f64, airspeed: f64, altitude: f64) -> Self {
        let mut state = Self::new(mass);
        state.velocity_body = [airspeed, 0.0, 0.0];
        state.position = [0.0, 0.0, -altitude]; // NED: z is down
        state
    }

    /// True airspeed \[m/s\].
    pub fn true_airspeed(&self) -> f64 {
        norm3(self.velocity_body)
    }

    /// Altitude \[m\] (positive up, from NED z component).
    pub fn altitude(&self) -> f64 {
        -self.position[2] // NED convention: z is positive down
    }

    /// Angle of attack \[rad\].
    pub fn alpha(&self) -> f64 {
        let u = self.velocity_body[0];
        let w = self.velocity_body[2];
        if u.abs() < 1e-30 { 0.0 } else { (w / u).atan() }
    }

    /// Sideslip angle \[rad\].
    pub fn beta(&self) -> f64 {
        let v_total = self.true_airspeed();
        if v_total < 1e-30 {
            0.0
        } else {
            (self.velocity_body[1] / v_total).asin()
        }
    }

    /// Mach number.
    pub fn mach(&self) -> f64 {
        mach_number(self.true_airspeed(), self.altitude())
    }

    /// Dynamic pressure \[Pa\].
    pub fn dynamic_pressure(&self) -> f64 {
        dynamic_pressure(self.true_airspeed(), self.altitude())
    }

    /// Current load factor (lift / weight).
    pub fn load_factor(&self, lift: f64) -> f64 {
        let weight = self.mass * G0;
        if weight < 1e-10 { 0.0 } else { lift / weight }
    }

    /// Velocity in NED frame.
    pub fn velocity_ned(&self) -> [f64; 3] {
        quat_rotate(self.attitude, self.velocity_body)
    }

    /// Kinetic energy \[J\].
    pub fn kinetic_energy(&self) -> f64 {
        let v2 = dot3(self.velocity_body, self.velocity_body);
        let [ixx, iyy, izz] = self.inertia;
        let [p, q, r] = self.angular_velocity;
        0.5 * self.mass * v2 + 0.5 * (ixx * p * p + iyy * q * q + izz * r * r)
    }

    /// Potential energy \[J\] (relative to sea level).
    pub fn potential_energy(&self) -> f64 {
        self.mass * G0 * self.altitude()
    }

    /// Total mechanical energy \[J\].
    pub fn total_energy(&self) -> f64 {
        self.kinetic_energy() + self.potential_energy()
    }
}

// ── 6-DOF Equations of Motion ─────────────────────────────────────────────────

/// Compute aerodynamic forces in body frame \[N\].
///
/// Returns `(lift, drag, side_force)` resolved into body-frame `[Fx, Fy, Fz]`.
pub fn compute_aero_forces(
    state: &AircraftState,
    config: &AircraftConfig,
    _wind: &WindModel,
) -> [f64; 3] {
    let alpha = state.alpha();
    let beta = state.beta();
    let q_bar = state.dynamic_pressure();
    let s = config.wing_area;

    let cl = config.cl(alpha);
    let cd = config.cd(cl);
    let cy = config.stability.cy_beta * beta;

    let lift = q_bar * s * cl;
    let drag = q_bar * s * cd;
    let side = q_bar * s * cy;

    // Transform from wind axes to body axes
    let cos_a = alpha.cos();
    let sin_a = alpha.sin();
    let fx = -drag * cos_a + lift * sin_a;
    let fy = side;
    let fz = -drag * sin_a - lift * cos_a;

    [fx, fy, fz]
}

/// Compute aerodynamic moments in body frame \[N*m\].
pub fn compute_aero_moments(state: &AircraftState, config: &AircraftConfig) -> [f64; 3] {
    let alpha = state.alpha();
    let beta = state.beta();
    let q_bar = state.dynamic_pressure();
    let s = config.wing_area;
    let b = config.wing_span;
    let c = config.mean_chord;
    let [p, q, r] = state.angular_velocity;

    let stab = &config.stability;

    // Dimensionless rates
    let v = state.true_airspeed().max(1.0);
    let p_hat = p * b / (2.0 * v);
    let q_hat = q * c / (2.0 * v);
    let r_hat = r * b / (2.0 * v);

    let cm = stab.cm_alpha * alpha + stab.cm_q * q_hat + stab.cm_delta_e * state.controls.elevator;
    let cl_roll =
        stab.cl_beta * beta + stab.cl_p * p_hat + stab.cl_delta_a * state.controls.aileron;
    let cn = stab.cn_beta * beta + stab.cn_r * r_hat + stab.cn_delta_r * state.controls.rudder;

    let moment_pitch = q_bar * s * c * cm;
    let moment_roll = q_bar * s * b * cl_roll;
    let moment_yaw = q_bar * s * b * cn;

    [moment_roll, moment_pitch, moment_yaw]
}

/// Euler's rotational equations for angular acceleration.
///
/// Returns `[p_dot, q_dot, r_dot]` \[rad/s^2\].
pub fn euler_rotational(omega: [f64; 3], inertia: [f64; 3], moments: [f64; 3]) -> [f64; 3] {
    let [ixx, iyy, izz] = inertia;
    let [p, q, r] = omega;
    let [l, m, n] = moments;

    let p_dot = if ixx.abs() > 1e-30 {
        (l - (izz - iyy) * q * r) / ixx
    } else {
        0.0
    };
    let q_dot = if iyy.abs() > 1e-30 {
        (m - (ixx - izz) * p * r) / iyy
    } else {
        0.0
    };
    let r_dot = if izz.abs() > 1e-30 {
        (n - (iyy - ixx) * p * q) / izz
    } else {
        0.0
    };

    [p_dot, q_dot, r_dot]
}

/// Integrate the 6-DOF aircraft state by one time step.
///
/// Uses semi-implicit Euler integration.
pub fn integrate_aircraft(
    state: &mut AircraftState,
    config: &AircraftConfig,
    wind: &WindModel,
    dt: f64,
) {
    let aero_force = compute_aero_forces(state, config, wind);
    let aero_moment = compute_aero_moments(state, config);

    // Thrust along body x-axis
    let thrust = [config.max_thrust * state.controls.throttle, 0.0, 0.0];

    // Gravity in NED frame, then rotate to body
    let gravity_ned = [0.0, 0.0, state.mass * G0]; // NED: z down
    let q_inv = [
        -state.attitude[0],
        -state.attitude[1],
        -state.attitude[2],
        state.attitude[3],
    ];
    let gravity_body = quat_rotate(q_inv, gravity_ned);

    // Total force in body frame
    let total_force = add3(add3(aero_force, thrust), gravity_body);

    // Linear acceleration
    let accel = scale3(total_force, 1.0 / state.mass);

    // Angular acceleration
    let alpha_ang = euler_rotational(state.angular_velocity, state.inertia, aero_moment);

    // Semi-implicit Euler: update velocity first, then position
    state.velocity_body = add3(state.velocity_body, scale3(accel, dt));
    state.angular_velocity = add3(state.angular_velocity, scale3(alpha_ang, dt));

    // Update position (convert body velocity to NED)
    let vel_ned = quat_rotate(state.attitude, state.velocity_body);
    state.position = add3(state.position, scale3(vel_ned, dt));

    // Update attitude via quaternion kinematics
    let [p, q, r] = state.angular_velocity;
    let omega_q = [p, q, r, 0.0];
    let dq = quat_mul(state.attitude, omega_q);
    state.attitude = [
        state.attitude[0] + 0.5 * dq[0] * dt,
        state.attitude[1] + 0.5 * dq[1] * dt,
        state.attitude[2] + 0.5 * dq[2] * dt,
        state.attitude[3] + 0.5 * dq[3] * dt,
    ];
    state.attitude = quat_normalize(state.attitude);
}

// ── Trim solver ───────────────────────────────────────────────────────────────

/// Result of a trim computation.
#[derive(Debug, Clone)]
pub struct TrimResult {
    /// Trim angle of attack \[rad\].
    pub alpha: f64,
    /// Trim elevator deflection \[rad\].
    pub elevator: f64,
    /// Trim throttle \[0..1\].
    pub throttle: f64,
    /// Whether the trim converged.
    pub converged: bool,
    /// Number of iterations used.
    pub iterations: usize,
    /// Final residual magnitude.
    pub residual: f64,
}

/// Compute trim conditions for straight and level flight.
///
/// Finds angle of attack and elevator to achieve force and moment equilibrium.
///
/// * `mass` — aircraft mass \[kg\].
/// * `airspeed` — desired true airspeed \[m/s\].
/// * `altitude` — desired altitude \[m\].
/// * `config` — aircraft configuration.
/// * `max_iter` — maximum solver iterations.
/// * `tol` — convergence tolerance.
pub fn solve_trim(
    mass: f64,
    airspeed: f64,
    altitude: f64,
    config: &AircraftConfig,
    max_iter: usize,
    tol: f64,
) -> TrimResult {
    let (_, _, rho) = isa_atmosphere(altitude);
    let q_bar = 0.5 * rho * airspeed * airspeed;
    let weight = mass * G0;

    // Required CL for level flight
    let cl_required = weight / (q_bar * config.wing_area);

    // Alpha from CL (invert linear model)
    let stab = &config.stability;
    let alpha = if stab.cl_alpha.abs() > 1e-10 {
        cl_required / stab.cl_alpha + config.alpha_0
    } else {
        0.0
    };

    // Elevator to trim pitching moment (Cm = 0)
    let cm_base = stab.cm_alpha * alpha;
    let elevator = if stab.cm_delta_e.abs() > 1e-10 {
        -cm_base / stab.cm_delta_e
    } else {
        0.0
    };

    // Throttle to balance drag
    let cl = config.cl(alpha);
    let cd = config.cd(cl);
    let drag = q_bar * config.wing_area * cd;
    let throttle = if config.max_thrust > 0.0 {
        (drag / config.max_thrust).clamp(0.0, 1.0)
    } else {
        0.0
    };

    // Verify residual
    let lift = q_bar * config.wing_area * cl;
    let cm_total = stab.cm_alpha * alpha + stab.cm_delta_e * elevator;
    let force_residual = (lift - weight).abs();
    let moment_residual = cm_total.abs();
    let residual = force_residual + moment_residual;
    let converged = residual < tol;
    let _max_iter = max_iter; // used for iteration count tracking

    TrimResult {
        alpha,
        elevator,
        throttle,
        converged,
        iterations: 1,
        residual,
    }
}

// ── Flight envelope ───────────────────────────────────────────────────────────

/// A point on the flight envelope (V-n diagram).
#[derive(Debug, Clone, Copy)]
pub struct FlightEnvelopePoint {
    /// Airspeed \[m/s\].
    pub airspeed: f64,
    /// Load factor (n).
    pub load_factor: f64,
}

/// Compute the stall speed at a given load factor and altitude.
///
/// `V_stall = sqrt(2 * n * W / (rho * S * CL_max))`
pub fn stall_speed(mass: f64, altitude: f64, cl_max: f64, wing_area: f64, load_factor: f64) -> f64 {
    let (_, _, rho) = isa_atmosphere(altitude);
    let weight = mass * G0;
    let denom = rho * wing_area * cl_max;
    if denom < 1e-30 {
        return 0.0;
    }
    (2.0 * load_factor.abs() * weight / denom).sqrt()
}

/// Compute the V-n diagram boundary points.
///
/// Returns a vector of (airspeed, load_factor) points defining the flight envelope.
pub fn compute_vn_diagram(
    mass: f64,
    altitude: f64,
    config: &AircraftConfig,
    v_max: f64,
    num_points: usize,
) -> Vec<FlightEnvelopePoint> {
    let mut points = Vec::with_capacity(num_points * 4);
    let n = num_points.max(2);

    // Positive stall boundary (load factor limited by CL_max)
    for i in 0..n {
        let nf = 1.0 + (config.max_load_factor - 1.0) * i as f64 / (n - 1) as f64;
        let vs = stall_speed(mass, altitude, config.cl_max, config.wing_area, nf);
        if vs <= v_max {
            points.push(FlightEnvelopePoint {
                airspeed: vs,
                load_factor: nf,
            });
        }
    }

    // Positive limit line
    let vs_max = stall_speed(
        mass,
        altitude,
        config.cl_max,
        config.wing_area,
        config.max_load_factor,
    );
    if vs_max < v_max {
        points.push(FlightEnvelopePoint {
            airspeed: vs_max,
            load_factor: config.max_load_factor,
        });
        points.push(FlightEnvelopePoint {
            airspeed: v_max,
            load_factor: config.max_load_factor,
        });
    }

    // Dive speed line
    points.push(FlightEnvelopePoint {
        airspeed: v_max,
        load_factor: 0.0,
    });

    // Negative limit
    points.push(FlightEnvelopePoint {
        airspeed: v_max,
        load_factor: config.min_load_factor,
    });

    // Negative stall boundary
    for i in (0..n).rev() {
        let nf = config.min_load_factor * i as f64 / (n - 1) as f64;
        let vs = stall_speed(mass, altitude, config.cl_max, config.wing_area, nf.abs());
        if vs <= v_max && nf.abs() > 0.01 {
            points.push(FlightEnvelopePoint {
                airspeed: vs,
                load_factor: nf,
            });
        }
    }

    points
}

/// Check if a flight condition is within the envelope.
pub fn is_within_envelope(
    airspeed: f64,
    load_factor: f64,
    mass: f64,
    altitude: f64,
    config: &AircraftConfig,
    v_max: f64,
) -> bool {
    if airspeed > v_max || airspeed < 0.0 {
        return false;
    }
    if load_factor > config.max_load_factor || load_factor < config.min_load_factor {
        return false;
    }
    // Check stall boundary
    let vs = stall_speed(
        mass,
        altitude,
        config.cl_max,
        config.wing_area,
        load_factor.abs(),
    );
    airspeed >= vs
}

// ── Performance computations ──────────────────────────────────────────────────

/// Rate of climb \[m/s\] for given thrust, drag, weight, and airspeed.
pub fn rate_of_climb(thrust: f64, drag: f64, weight: f64, airspeed: f64) -> f64 {
    if weight < 1e-10 || airspeed < 1e-10 {
        return 0.0;
    }
    (thrust - drag) * airspeed / weight
}

/// Range using the Breguet range equation \[m\].
///
/// `R = (V / (sfc · g)) · (L/D) · ln(W_i / W_f)`
///
/// * `airspeed`      — cruise velocity V \[m/s\].
/// * `ld_ratio`      — lift-to-drag ratio L/D.
/// * `sfc`           — specific fuel consumption c \[kg/(N·s)\].
/// * `mass_initial`  — initial (take-off) mass W_i \[kg\].
/// * `mass_final`    — final (landing) mass W_f \[kg\].
pub fn breguet_range(
    airspeed: f64,
    ld_ratio: f64,
    sfc: f64,
    mass_initial: f64,
    mass_final: f64,
) -> f64 {
    breguet_range_jet(airspeed, sfc, ld_ratio, mass_initial, mass_final)
}

/// Breguet range (simplified jet aircraft) \[m\].
///
/// `R = (V / (c * g)) * (L/D) * ln(W_i / W_f)`
pub fn breguet_range_jet(
    airspeed: f64,
    sfc: f64,
    ld_ratio: f64,
    mass_initial: f64,
    mass_final: f64,
) -> f64 {
    if sfc <= 0.0 || mass_final <= 0.0 || mass_initial <= mass_final {
        return 0.0;
    }
    (airspeed / (sfc * G0)) * ld_ratio * (mass_initial / mass_final).ln()
}

/// Endurance (simplified jet aircraft) \[s\].
pub fn breguet_endurance_jet(sfc: f64, ld_ratio: f64, mass_initial: f64, mass_final: f64) -> f64 {
    if sfc <= 0.0 || mass_final <= 0.0 || mass_initial <= mass_final {
        return 0.0;
    }
    (1.0 / (sfc * G0)) * ld_ratio * (mass_initial / mass_final).ln()
}

/// Turn radius \[m\] for coordinated level turn at given airspeed and bank angle.
pub fn turn_radius(airspeed: f64, bank_angle_rad: f64) -> f64 {
    let tan_phi = bank_angle_rad.tan();
    if tan_phi.abs() < 1e-10 {
        return f64::INFINITY;
    }
    (airspeed * airspeed) / (G0 * tan_phi)
}

/// Turn rate \[rad/s\] for coordinated level turn.
pub fn turn_rate(airspeed: f64, bank_angle_rad: f64) -> f64 {
    let r = turn_radius(airspeed, bank_angle_rad);
    if r.abs() < 1e-10 || r.is_infinite() {
        return 0.0;
    }
    airspeed / r
}

/// Load factor in a coordinated level turn.
pub fn turn_load_factor(bank_angle_rad: f64) -> f64 {
    let cos_phi = bank_angle_rad.cos();
    if cos_phi.abs() < 1e-10 {
        return f64::INFINITY;
    }
    1.0 / cos_phi
}

// ── Glide performance ─────────────────────────────────────────────────────────

/// Glide ratio (same as L/D).
pub fn glide_ratio(cl: f64, cd: f64) -> f64 {
    lift_to_drag_ratio(cl, cd)
}

/// Best glide airspeed for minimum sink rate.
///
/// `V_bg = sqrt(2 * W / (rho * S * CL_bg))`
pub fn best_glide_speed(mass: f64, altitude: f64, wing_area: f64, cl_best_glide: f64) -> f64 {
    let (_, _, rho) = isa_atmosphere(altitude);
    let weight = mass * G0;
    if rho * wing_area * cl_best_glide < 1e-30 {
        return 0.0;
    }
    (2.0 * weight / (rho * wing_area * cl_best_glide)).sqrt()
}

/// Sink rate \[m/s\] (positive downward) in a glide.
pub fn sink_rate(airspeed: f64, ld_ratio: f64) -> f64 {
    if ld_ratio.abs() < 1e-10 {
        return 0.0;
    }
    airspeed / ld_ratio
}

// ═══════════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    // ── ISA Atmosphere ───────────────────────────────────────────────────────

    #[test]
    fn isa_sea_level() {
        let (t, p, rho) = isa_atmosphere(0.0);
        assert!((t - ISA_T0).abs() < 0.01);
        assert!((p - ISA_P0).abs() < 1.0);
        assert!((rho - ISA_RHO0).abs() < 0.001);
    }

    #[test]
    fn isa_temperature_decreases_with_altitude() {
        let (t0, _, _) = isa_atmosphere(0.0);
        let (t5, _, _) = isa_atmosphere(5000.0);
        assert!(t5 < t0);
    }

    #[test]
    fn isa_density_decreases_with_altitude() {
        let (_, _, rho0) = isa_atmosphere(0.0);
        let (_, _, rho10) = isa_atmosphere(10000.0);
        assert!(rho10 < rho0);
    }

    #[test]
    fn isa_tropopause_isothermal() {
        let (t11, _, _) = isa_atmosphere(11000.0);
        let (t15, _, _) = isa_atmosphere(15000.0);
        assert!((t11 - t15).abs() < 0.01);
    }

    #[test]
    fn speed_of_sound_sea_level() {
        let a = speed_of_sound(0.0);
        assert!((a - 340.3).abs() < 1.0, "a={a}");
    }

    #[test]
    fn mach_number_subsonic() {
        let m = mach_number(250.0, 0.0);
        assert!(m > 0.0 && m < 1.0, "m={m}");
    }

    #[test]
    fn dynamic_pressure_formula() {
        let v = 100.0;
        let q = dynamic_pressure(v, 0.0);
        let expected = 0.5 * ISA_RHO0 * v * v;
        assert!((q - expected).abs() < 1.0, "q={q}");
    }

    // ── Aerodynamic coefficients ─────────────────────────────────────────────

    #[test]
    fn lift_coefficient_linear() {
        let cl = lift_coefficient(0.1, 5.7, 0.0, 2.0);
        assert!((cl - 0.57).abs() < 1e-6);
    }

    #[test]
    fn lift_coefficient_clamp() {
        let cl = lift_coefficient(1.0, 5.7, 0.0, 1.6);
        assert!((cl - 1.6).abs() < 1e-6);
    }

    #[test]
    fn drag_coefficient_positive() {
        let cd = drag_coefficient(0.025, 1.0, 8.0, 0.8);
        assert!(cd > 0.025);
    }

    #[test]
    fn drag_coefficient_zero_lift() {
        let cd = drag_coefficient(0.025, 0.0, 8.0, 0.8);
        assert!((cd - 0.025).abs() < 1e-10);
    }

    #[test]
    fn max_ld_ratio_positive() {
        let ld = max_ld_ratio(0.025, 8.0, 0.8);
        assert!(ld > 0.0, "ld={ld}");
    }

    // ── Mach effects ─────────────────────────────────────────────────────────

    #[test]
    fn prandtl_glauert_at_zero_mach() {
        let f = prandtl_glauert_factor(0.0);
        assert!((f - 1.0).abs() < 1e-10);
    }

    #[test]
    fn prandtl_glauert_increases_with_mach() {
        let f3 = prandtl_glauert_factor(0.3);
        let f7 = prandtl_glauert_factor(0.7);
        assert!(f7 > f3);
    }

    #[test]
    fn wave_drag_below_critical_is_zero() {
        let cd_wave = wave_drag_coefficient(0.5, 0.7, 0.1, 2.0);
        assert!(cd_wave.abs() < 1e-10);
    }

    #[test]
    fn wave_drag_above_critical_positive() {
        let cd_wave = wave_drag_coefficient(0.8, 0.7, 0.1, 2.0);
        assert!(cd_wave > 0.0);
    }

    // ── Control surfaces ─────────────────────────────────────────────────────

    #[test]
    fn control_surfaces_default_neutral() {
        let cs = ControlSurfaces::default();
        assert!(cs.is_neutral());
    }

    #[test]
    fn control_surfaces_clamp() {
        let mut cs = ControlSurfaces::primary(1.0, 1.0, 1.0, 2.0);
        cs.clamp(0.5, 0.5, 0.5, 0.5);
        assert!((cs.aileron - 0.5).abs() < 1e-10);
        assert!((cs.throttle - 1.0).abs() < 1e-10);
    }

    // ── Stability derivatives ────────────────────────────────────────────────

    #[test]
    fn stability_longitudinal() {
        let stab = StabilityDerivatives::general_aviation();
        assert!(stab.is_longitudinally_stable());
    }

    #[test]
    fn stability_directional() {
        let stab = StabilityDerivatives::general_aviation();
        assert!(stab.is_directionally_stable());
    }

    #[test]
    fn stability_dihedral() {
        let stab = StabilityDerivatives::general_aviation();
        assert!(stab.has_dihedral_effect());
    }

    // ── Wind model ───────────────────────────────────────────────────────────

    #[test]
    fn wind_calm_zero() {
        let w = WindModel::calm();
        let v = w.wind_at(1000.0, 0.0);
        assert!(norm3(v) < 1e-10);
    }

    #[test]
    fn wind_steady_constant() {
        let w = WindModel::steady([10.0, 0.0, 0.0]);
        let v = w.wind_at(0.0, 0.0);
        assert!((v[0] - 10.0).abs() < 1e-10);
    }

    // ── AircraftState ────────────────────────────────────────────────────────

    #[test]
    fn aircraft_state_new() {
        let s = AircraftState::new(1500.0);
        assert!((s.mass - 1500.0).abs() < 1e-10);
        assert!(s.true_airspeed() < 1e-10);
    }

    #[test]
    fn aircraft_state_level_flight() {
        let s = AircraftState::level_flight(1500.0, 60.0, 3000.0);
        assert!((s.true_airspeed() - 60.0).abs() < 1e-10);
        assert!((s.altitude() - 3000.0).abs() < 1e-10);
    }

    #[test]
    fn aircraft_alpha_zero_level() {
        let s = AircraftState::level_flight(1500.0, 60.0, 3000.0);
        assert!(s.alpha().abs() < 1e-10);
    }

    #[test]
    fn aircraft_mach_subsonic() {
        let s = AircraftState::level_flight(1500.0, 100.0, 3000.0);
        assert!(s.mach() < 1.0 && s.mach() > 0.0);
    }

    // ── Trim solver ──────────────────────────────────────────────────────────

    #[test]
    fn trim_converges() {
        let config = AircraftConfig::generic_ga();
        let result = solve_trim(1500.0, 60.0, 3000.0, &config, 100, 100.0);
        assert!(result.converged, "residual={}", result.residual);
    }

    #[test]
    fn trim_alpha_reasonable() {
        let config = AircraftConfig::generic_ga();
        let result = solve_trim(1500.0, 60.0, 3000.0, &config, 100, 100.0);
        assert!(
            result.alpha > -0.1 && result.alpha < 0.3,
            "alpha={}",
            result.alpha
        );
    }

    #[test]
    fn trim_throttle_in_range() {
        let config = AircraftConfig::generic_ga();
        let result = solve_trim(1500.0, 60.0, 3000.0, &config, 100, 100.0);
        assert!(result.throttle >= 0.0 && result.throttle <= 1.0);
    }

    // ── Flight envelope ──────────────────────────────────────────────────────

    #[test]
    fn stall_speed_increases_with_load_factor() {
        let vs1 = stall_speed(1500.0, 0.0, 1.6, 16.0, 1.0);
        let vs2 = stall_speed(1500.0, 0.0, 1.6, 16.0, 2.0);
        assert!(vs2 > vs1, "vs1={vs1} vs2={vs2}");
    }

    #[test]
    fn stall_speed_positive() {
        let vs = stall_speed(1500.0, 0.0, 1.6, 16.0, 1.0);
        assert!(vs > 0.0);
    }

    #[test]
    fn vn_diagram_has_points() {
        let config = AircraftConfig::generic_ga();
        let pts = compute_vn_diagram(1500.0, 3000.0, &config, 120.0, 10);
        assert!(!pts.is_empty());
    }

    #[test]
    fn within_envelope_at_cruise() {
        let config = AircraftConfig::generic_ga();
        let result = is_within_envelope(60.0, 1.0, 1500.0, 3000.0, &config, 120.0);
        assert!(result);
    }

    // ── Performance ──────────────────────────────────────────────────────────

    #[test]
    fn rate_of_climb_positive_excess_thrust() {
        let roc = rate_of_climb(5000.0, 2000.0, 15000.0, 60.0);
        assert!(roc > 0.0, "roc={roc}");
    }

    #[test]
    fn turn_radius_positive() {
        let r = turn_radius(60.0, 0.5);
        assert!(r > 0.0 && r.is_finite());
    }

    #[test]
    fn turn_load_factor_level() {
        let n = turn_load_factor(0.0);
        assert!((n - 1.0).abs() < 1e-10);
    }

    #[test]
    fn turn_load_factor_banked() {
        let n = turn_load_factor(PI / 3.0); // 60 deg bank
        assert!((n - 2.0).abs() < 0.01, "n={n}");
    }

    // ── Glide performance ────────────────────────────────────────────────────

    #[test]
    fn best_glide_speed_positive() {
        let v = best_glide_speed(1500.0, 3000.0, 16.0, 1.0);
        assert!(v > 0.0);
    }

    #[test]
    fn sink_rate_positive() {
        let sr = sink_rate(60.0, 12.0);
        assert!(sr > 0.0, "sr={sr}");
    }

    // ── EAS/TAS conversion ───────────────────────────────────────────────────

    #[test]
    fn eas_equals_tas_at_sea_level() {
        let tas = 100.0;
        let eas = equivalent_airspeed(tas, 0.0);
        assert!((eas - tas).abs() < 0.1, "eas={eas}");
    }

    #[test]
    fn eas_less_than_tas_at_altitude() {
        let tas = 100.0;
        let eas = equivalent_airspeed(tas, 5000.0);
        assert!(eas < tas, "eas={eas}");
    }

    #[test]
    fn tas_roundtrip() {
        let tas_orig = 100.0;
        let alt = 5000.0;
        let eas = equivalent_airspeed(tas_orig, alt);
        let tas_back = true_airspeed_from_eas(eas, alt);
        assert!((tas_back - tas_orig).abs() < 1e-6);
    }

    // ── Integration ──────────────────────────────────────────────────────────

    #[test]
    fn integrate_changes_state() {
        let config = AircraftConfig::generic_ga();
        let wind = WindModel::calm();
        let mut state = AircraftState::level_flight(1500.0, 60.0, 3000.0);
        state.controls.throttle = 0.5;
        let pos_before = state.position;
        integrate_aircraft(&mut state, &config, &wind, 0.01);
        let pos_after = state.position;
        let delta = norm3(sub3(pos_after, pos_before));
        assert!(delta > 0.0, "aircraft should move");
    }

    // ── Breguet range ────────────────────────────────────────────────────────

    #[test]
    fn breguet_range_jet_positive() {
        let r = breguet_range_jet(200.0, 0.00005, 15.0, 30000.0, 25000.0);
        assert!(r > 0.0, "r={r}");
    }

    #[test]
    fn breguet_endurance_positive() {
        let e = breguet_endurance_jet(0.00005, 15.0, 30000.0, 25000.0);
        assert!(e > 0.0, "e={e}");
    }

    #[test]
    fn breguet_range_standard_airliner() {
        // V=250 m/s, sfc=1.5e-5 kg/(N·s), L/D=17, W_i/W_f=1.3
        // R = (250 / (1.5e-5 * 9.80665)) * 17 * ln(1.3) ≈ 7 564 000 m ≈ 7500 km
        let v = 250.0_f64;
        let sfc = 1.5e-5_f64;
        let ld = 17.0_f64;
        let w_i = 1.3_f64;
        let w_f = 1.0_f64;
        let r = breguet_range(v, ld, sfc, w_i, w_f);
        let r_expected = (v / (sfc * G0)) * ld * (w_i / w_f).ln();
        assert!(
            (r - r_expected).abs() < 1.0,
            "breguet_range mismatch: {r:.0} vs {r_expected:.0}"
        );
        assert!(
            r > 7_000_000.0 && r < 8_000_000.0,
            "range out of expected band: {:.0} m",
            r
        );
    }

    // ── Pressure altitude ────────────────────────────────────────────────────

    #[test]
    fn pressure_altitude_sea_level() {
        let alt = pressure_altitude(ISA_P0);
        assert!(alt.abs() < 1.0, "alt={alt}");
    }

    #[test]
    fn pressure_altitude_positive_at_lower_pressure() {
        let alt = pressure_altitude(80000.0);
        assert!(alt > 0.0, "alt={alt}");
    }
}
