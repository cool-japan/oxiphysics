// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Aerospace rigid body dynamics with 6-DOF flight dynamics.
//!
//! Provides aerodynamic forces (lift/drag coefficients), control surfaces
//! (aileron/elevator/rudder), ISA atmosphere model, wind effects, stability
//! derivatives, trim solver, flight envelope analysis, and Mach number effects.
//!
//! The aerodynamic model is driven by the *air-relative* velocity (ground
//! velocity minus the wind at the aircraft's altitude and time, see
//! [`AircraftState::air_data`]), and the body-frame translational equations
//! include the `omega x v` transport term, so the ballistic motion of a rotating
//! body is reproduced (to the first-order accuracy of the integrator).
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

use std::f64::consts::PI;

// ── vector helpers ────────────────────────────────────────────────────────────

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

/// Cross product `a x b`.
#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
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

/// Quaternion conjugate (the inverse of a unit quaternion). Convention: q = \[x, y, z, w\].
#[inline]
fn quat_conjugate(q: [f64; 4]) -> [f64; 4] {
    [-q[0], -q[1], -q[2], q[3]]
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
    let q_conj = quat_conjugate(q);
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
    /// Simulation time \[s\], advanced by [`integrate_aircraft`]. It is the
    /// time argument of [`WindModel::wind_at`], so gusts evolve as the
    /// simulation runs.
    pub time: f64,
}

/// Flight condition of the aircraft relative to the *air mass*.
///
/// The aerodynamic forces and moments depend on the velocity of the aircraft
/// with respect to the surrounding air, not with respect to the ground; the two
/// differ by the local wind. Build one with [`AircraftState::air_data`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AirData {
    /// Velocity relative to the air mass, resolved in the body frame
    /// `[u, v, w]` \[m/s\].
    pub velocity_body: [f64; 3],
    /// True airspeed `|v_air|` \[m/s\].
    pub airspeed: f64,
    /// Angle of attack \[rad\].
    pub alpha: f64,
    /// Sideslip angle \[rad\].
    pub beta: f64,
    /// Dynamic pressure `0.5 * rho * |v_air|^2` \[Pa\].
    pub dynamic_pressure: f64,
    /// Mach number of the air-relative flow.
    pub mach: f64,
}

/// Angle of attack and sideslip `(alpha, beta)` \[rad\] of a body-frame flow
/// velocity `[u, v, w]`: `alpha = atan2(w, u)`, `beta = asin(v / |V|)`.
/// Both are zero for a vanishing flow.
fn flow_angles(velocity_body: [f64; 3]) -> (f64, f64) {
    let speed = norm3(velocity_body);
    if speed < 1e-30 {
        return (0.0, 0.0);
    }
    let alpha = velocity_body[2].atan2(velocity_body[0]);
    let beta = (velocity_body[1] / speed).clamp(-1.0, 1.0).asin();
    (alpha, beta)
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
            time: 0.0,
        }
    }

    /// Create a state in level flight at given airspeed and altitude.
    pub fn level_flight(mass: f64, airspeed: f64, altitude: f64) -> Self {
        let mut state = Self::new(mass);
        state.velocity_body = [airspeed, 0.0, 0.0];
        state.position = [0.0, 0.0, -altitude]; // NED: z is down
        state
    }

    /// Speed of the aircraft relative to the ground \[m/s\]; equal to the true
    /// airspeed in calm air. Use [`AircraftState::air_data`] when there is wind.
    pub fn true_airspeed(&self) -> f64 {
        norm3(self.velocity_body)
    }

    /// Altitude \[m\] (positive up, from NED z component).
    pub fn altitude(&self) -> f64 {
        -self.position[2] // NED convention: z is positive down
    }

    /// Angle of attack of the ground-relative velocity \[rad\] (calm air).
    /// Use [`AircraftState::air_data`] when there is wind.
    pub fn alpha(&self) -> f64 {
        flow_angles(self.velocity_body).0
    }

    /// Sideslip angle of the ground-relative velocity \[rad\] (calm air).
    /// Use [`AircraftState::air_data`] when there is wind.
    pub fn beta(&self) -> f64 {
        flow_angles(self.velocity_body).1
    }

    /// Mach number of the ground-relative velocity (calm air).
    pub fn mach(&self) -> f64 {
        mach_number(self.true_airspeed(), self.altitude())
    }

    /// Dynamic pressure of the ground-relative velocity \[Pa\] (calm air).
    pub fn dynamic_pressure(&self) -> f64 {
        dynamic_pressure(self.true_airspeed(), self.altitude())
    }

    /// Wind velocity at the aircraft's current altitude and time, resolved in
    /// the body frame \[m/s\].
    ///
    /// [`WindModel::wind_at`] returns the wind in the NED frame; the attitude
    /// quaternion maps body to NED, so its conjugate maps the wind back.
    pub fn wind_body(&self, wind: &WindModel) -> [f64; 3] {
        let wind_ned = wind.wind_at(self.altitude(), self.time);
        quat_rotate(quat_conjugate(self.attitude), wind_ned)
    }

    /// Flight condition relative to the air mass: the body-frame velocity minus
    /// the body-frame wind, and the airspeed, flow angles, dynamic pressure and
    /// Mach number derived from it.
    ///
    /// An aircraft drifting with the wind (`velocity == wind`) has zero airspeed
    /// and therefore feels no aerodynamic force, whatever its ground speed.
    pub fn air_data(&self, wind: &WindModel) -> AirData {
        let velocity_air = sub3(self.velocity_body, self.wind_body(wind));
        let airspeed = norm3(velocity_air);
        let altitude = self.altitude();
        let (alpha, beta) = flow_angles(velocity_air);
        AirData {
            velocity_body: velocity_air,
            airspeed,
            alpha,
            beta,
            dynamic_pressure: dynamic_pressure(airspeed, altitude),
            mach: mach_number(airspeed, altitude),
        }
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
///
/// The angle of attack, sideslip and dynamic pressure are evaluated from the
/// velocity relative to the air mass (ground velocity minus `wind` at the
/// aircraft's altitude and time, see [`AircraftState::air_data`]).
pub fn compute_aero_forces(
    state: &AircraftState,
    config: &AircraftConfig,
    wind: &WindModel,
) -> [f64; 3] {
    let air = state.air_data(wind);
    let alpha = air.alpha;
    let beta = air.beta;
    let q_bar = air.dynamic_pressure;
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
///
/// Like [`compute_aero_forces`], the flow angles, dynamic pressure and the
/// non-dimensional body rates use the velocity relative to the air mass, so the
/// moments respond to `wind`.
pub fn compute_aero_moments(
    state: &AircraftState,
    config: &AircraftConfig,
    wind: &WindModel,
) -> [f64; 3] {
    let air = state.air_data(wind);
    let alpha = air.alpha;
    let beta = air.beta;
    let q_bar = air.dynamic_pressure;
    let s = config.wing_area;
    let b = config.wing_span;
    let c = config.mean_chord;
    let [p, q, r] = state.angular_velocity;

    let stab = &config.stability;

    // Dimensionless rates
    let v = air.airspeed.max(1.0);
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

/// Translational acceleration of a body-fixed frame \[m/s^2\].
///
/// Newton's second law in a rotating frame, `m (dv/dt + omega x v) = F`, gives
/// `dv/dt = F / m - omega x v`, where `v` and `F` are the body-frame velocity
/// and force and `omega` the body angular velocity. The `omega x v` transport
/// term turns the velocity vector with the body; without it a freely moving
/// rotating body would be accelerated even when no force acts on it.
pub fn body_frame_acceleration(
    force_body: [f64; 3],
    mass: f64,
    angular_velocity: [f64; 3],
    velocity_body: [f64; 3],
) -> [f64; 3] {
    sub3(
        scale3(force_body, 1.0 / mass),
        cross3(angular_velocity, velocity_body),
    )
}

/// Integrate the 6-DOF aircraft state by one time step.
///
/// Uses semi-implicit Euler integration. The aerodynamic forces and moments are
/// evaluated from the velocity relative to the air mass (`wind` sampled at the
/// state's altitude and `state.time`), the translational dynamics include the
/// `omega x v` transport term of the rotating body frame (see
/// [`body_frame_acceleration`]), and `state.time` is advanced by `dt`.
pub fn integrate_aircraft(
    state: &mut AircraftState,
    config: &AircraftConfig,
    wind: &WindModel,
    dt: f64,
) {
    let aero_force = compute_aero_forces(state, config, wind);
    let aero_moment = compute_aero_moments(state, config, wind);

    // Thrust along body x-axis
    let thrust = [config.max_thrust * state.controls.throttle, 0.0, 0.0];

    // Gravity in NED frame, then rotate to body
    let gravity_ned = [0.0, 0.0, state.mass * G0]; // NED: z down
    let gravity_body = quat_rotate(quat_conjugate(state.attitude), gravity_ned);

    // Total force in body frame
    let total_force = add3(add3(aero_force, thrust), gravity_body);

    // Linear acceleration in the rotating body frame (includes -omega x v)
    let accel = body_frame_acceleration(
        total_force,
        state.mass,
        state.angular_velocity,
        state.velocity_body,
    );

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
    state.time += dt;
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

    // ── Wind and body-frame transport (F3957) ────────────────────────────────

    /// Attitude quaternion \[x, y, z, w\] of a pure yaw by `angle` about the NED
    /// z (down) axis; body x then points along the heading.
    fn yaw_quaternion(angle: f64) -> [f64; 4] {
        [0.0, 0.0, (0.5 * angle).sin(), (0.5 * angle).cos()]
    }

    fn close(a: f64, b: f64, rel_tol: f64) -> bool {
        (a - b).abs() <= rel_tol * a.abs().max(b.abs()).max(1e-12)
    }

    fn flying_north(airspeed: f64) -> AircraftState {
        AircraftState::level_flight(1500.0, airspeed, 3000.0)
    }

    #[test]
    fn air_data_in_calm_air_equals_ground_relative_values() {
        let mut s = flying_north(60.0);
        s.velocity_body = [60.0, 4.0, 8.0];
        let air = s.air_data(&WindModel::calm());
        assert!(close(air.airspeed, s.true_airspeed(), 1e-12));
        assert!(close(air.alpha, s.alpha(), 1e-12));
        assert!(close(air.beta, s.beta(), 1e-12));
        assert!(close(air.dynamic_pressure, s.dynamic_pressure(), 1e-12));
        assert!(close(air.mach, s.mach(), 1e-12));
        assert!(norm3(sub3(air.velocity_body, s.velocity_body)) < 1e-12);
        assert!(air.alpha > 0.0 && air.beta > 0.0);
    }

    #[test]
    fn headwind_and_tailwind_shift_airspeed_by_the_wind_speed() {
        let s = flying_north(60.0);
        // Heading north at identity attitude: wind blowing south is a headwind.
        let head = s.air_data(&WindModel::steady([-10.0, 0.0, 0.0]));
        assert!((head.airspeed - 70.0).abs() < 1e-12, "{}", head.airspeed);
        let tail = s.air_data(&WindModel::steady([10.0, 0.0, 0.0]));
        assert!((tail.airspeed - 50.0).abs() < 1e-12, "{}", tail.airspeed);
        // Dynamic pressure follows the air-relative speed, not the ground speed.
        assert!(close(
            head.dynamic_pressure,
            dynamic_pressure(70.0, s.altitude()),
            1e-12
        ));
        assert!(close(
            tail.dynamic_pressure,
            dynamic_pressure(50.0, s.altitude()),
            1e-12
        ));
        // Neither changes the flow direction.
        assert!(head.alpha.abs() < 1e-12 && head.beta.abs() < 1e-12);
        assert!(tail.alpha.abs() < 1e-12 && tail.beta.abs() < 1e-12);
    }

    #[test]
    fn updraft_raises_the_angle_of_attack() {
        let s = flying_north(60.0);
        // NED z points down, so a negative z wind component is an updraft.
        let air = s.air_data(&WindModel::steady([0.0, 0.0, -5.0]));
        assert!(
            (air.alpha - 5.0_f64.atan2(60.0)).abs() < 1e-12,
            "{}",
            air.alpha
        );
        assert!(air.alpha > 0.0);
    }

    #[test]
    fn wind_is_rotated_into_the_body_frame() {
        let mut s = flying_north(60.0);
        s.attitude = yaw_quaternion(PI / 2.0); // heading east
        // Wind blowing east is a tailwind for an eastbound aircraft.
        let tail = WindModel::steady([0.0, 10.0, 0.0]);
        let wb = s.wind_body(&tail);
        assert!((wb[0] - 10.0).abs() < 1e-12 && wb[1].abs() < 1e-12 && wb[2].abs() < 1e-12);
        assert!((s.air_data(&tail).airspeed - 50.0).abs() < 1e-12);
        // Wind blowing north hits the right wing (body y points south): positive
        // sideslip, with the airspeed vector gaining a +y component.
        let from_right = WindModel::steady([10.0, 0.0, 0.0]);
        let air = s.air_data(&from_right);
        assert!((s.wind_body(&from_right)[1] + 10.0).abs() < 1e-12);
        assert!((air.velocity_body[1] - 10.0).abs() < 1e-12);
        assert!(air.beta > 0.0);
        assert!((air.beta - (10.0_f64 / (3700.0_f64).sqrt()).asin()).abs() < 1e-12);
    }

    #[test]
    fn aero_forces_and_moments_depend_only_on_air_relative_velocity() {
        let config = AircraftConfig::generic_ga();
        let mut calm = flying_north(60.0);
        calm.velocity_body = [60.0, 3.0, 8.0];
        calm.angular_velocity = [0.05, -0.04, 0.03];
        calm.controls = ControlSurfaces::primary(0.02, -0.03, 0.01, 0.0);
        // Same air-relative flow, but 10 m/s higher over the ground in a tailwind.
        let mut windy = calm.clone();
        windy.velocity_body = [70.0, 3.0, 8.0];
        let tailwind = WindModel::steady([10.0, 0.0, 0.0]);
        let f_calm = compute_aero_forces(&calm, &config, &WindModel::calm());
        let f_windy = compute_aero_forces(&windy, &config, &tailwind);
        let m_calm = compute_aero_moments(&calm, &config, &WindModel::calm());
        let m_windy = compute_aero_moments(&windy, &config, &tailwind);
        for k in 0..3 {
            assert!(
                close(f_calm[k], f_windy[k], 1e-9),
                "F[{k}]: {f_calm:?} vs {f_windy:?}"
            );
            assert!(
                close(m_calm[k], m_windy[k], 1e-9),
                "M[{k}]: {m_calm:?} vs {m_windy:?}"
            );
        }
        // Ignoring the wind (the old behaviour) would give a different force.
        let f_ignoring_wind = compute_aero_forces(&windy, &config, &WindModel::calm());
        assert!(norm3(sub3(f_ignoring_wind, f_windy)) > 100.0);
    }

    #[test]
    fn drifting_with_the_wind_feels_no_aerodynamic_force() {
        let config = AircraftConfig::generic_ga();
        let mut s = flying_north(0.0);
        s.velocity_body = [10.0, 0.0, 0.0];
        s.controls = ControlSurfaces::primary(0.1, 0.1, 0.1, 0.0);
        let wind = WindModel::steady([10.0, 0.0, 0.0]);
        let air = s.air_data(&wind);
        assert!(air.airspeed < 1e-12);
        assert!(air.alpha == 0.0 && air.beta == 0.0);
        assert!(norm3(compute_aero_forces(&s, &config, &wind)) < 1e-9);
        assert!(norm3(compute_aero_moments(&s, &config, &wind)) < 1e-9);
    }

    #[test]
    fn crosswind_produces_sideslip_and_a_restoring_side_force() {
        let config = AircraftConfig::generic_ga();
        let s = flying_north(60.0);
        // Air moving west: relative to the air the aircraft moves east (+y body).
        let wind = WindModel::steady([0.0, -10.0, 0.0]);
        let air = s.air_data(&wind);
        let beta = (10.0_f64 / (3700.0_f64).sqrt()).asin();
        assert!((air.beta - beta).abs() < 1e-12);
        let force = compute_aero_forces(&s, &config, &wind);
        let expected_side =
            air.dynamic_pressure * config.wing_area * config.stability.cy_beta * beta;
        assert!(
            close(force[1], expected_side, 1e-12),
            "{} vs {expected_side}",
            force[1]
        );
        assert!(force[1] < 0.0, "side force must oppose positive sideslip");
        // Weathercock stability: positive sideslip gives a positive yawing moment.
        let moment = compute_aero_moments(&s, &config, &wind);
        assert!(moment[2] > 0.0, "yawing moment {}", moment[2]);
    }

    #[test]
    fn gust_is_sampled_at_the_state_time() {
        let mut wind = WindModel::calm();
        wind.gust_amplitude = 5.0;
        wind.gust_frequency = 2.0;
        let mut s = flying_north(60.0);
        assert!(
            (s.air_data(&wind).airspeed - 60.0).abs() < 1e-12,
            "no gust at t = 0"
        );
        s.time = PI / 4.0; // sin(2 t) = 1: a 5 m/s tailwind gust
        assert!((s.air_data(&wind).airspeed - 55.0).abs() < 1e-12);
        s.time = 3.0 * PI / 4.0; // sin(2 t) = -1: a 5 m/s headwind gust
        assert!((s.air_data(&wind).airspeed - 65.0).abs() < 1e-12);
    }

    #[test]
    fn integrate_aircraft_advances_the_simulation_time() {
        let config = AircraftConfig::generic_ga();
        let wind = WindModel::calm();
        let mut s = flying_north(60.0);
        assert_eq!(s.time, 0.0);
        for _ in 0..10 {
            integrate_aircraft(&mut s, &config, &wind, 0.01);
        }
        assert!((s.time - 0.1).abs() < 1e-12, "time = {}", s.time);
    }

    #[test]
    fn body_frame_acceleration_contains_the_transport_term() {
        let v = [60.0, -3.0, 8.0];
        let w = [0.1, 0.5, -0.2];
        // Without force the acceleration is -omega x v, perpendicular to v, so
        // the speed (the norm of the body-frame velocity) is stationary.
        let a = body_frame_acceleration([0.0; 3], 1500.0, w, v);
        let wxv = cross3(w, v);
        for k in 0..3 {
            assert!((a[k] + wxv[k]).abs() < 1e-12);
        }
        assert!(dot3(a, v).abs() < 1e-9, "a.v = {}", dot3(a, v));
        // A pure pitch rate q turns the (u, w) components: du/dt = -q w, dw/dt = q u.
        let a = body_frame_acceleration([0.0; 3], 1500.0, [0.0, 0.5, 0.0], [60.0, 0.0, 8.0]);
        assert!((a[0] - (-0.5 * 8.0)).abs() < 1e-12 && (a[2] - (0.5 * 60.0)).abs() < 1e-12);
        // With a force, F/m is added.
        let a = body_frame_acceleration([3000.0, 0.0, -1500.0], 1500.0, w, v);
        assert!((a[0] - (2.0 - wxv[0])).abs() < 1e-12);
        assert!((a[2] - (-1.0 - wxv[2])).abs() < 1e-12);
    }

    #[test]
    fn rotating_body_without_aerodynamics_follows_a_ballistic_path() {
        // With the wing removed there is no aerodynamic force or moment and no
        // thrust, so only gravity acts. Whatever the body does, the NED velocity
        // must change by g*t along "down" and the path must be a parabola.
        // Dropping the omega x v term makes the body-frame velocity rotate with
        // the airframe, so the NED velocity would swing by about q*|v|*t = 30 m/s.
        let mut config = AircraftConfig::generic_ga();
        config.wing_area = 0.0;
        let wind = WindModel::calm();
        let mut s = AircraftState::level_flight(1500.0, 60.0, 3000.0);
        s.velocity_body = [60.0, 0.0, 5.0];
        s.angular_velocity = [0.0, 0.5, 0.0]; // constant pitch rate
        let v0 = s.velocity_ned();
        let p0 = s.position;
        let dt = 1.0e-3;
        let steps = 1000;
        for _ in 0..steps {
            integrate_aircraft(&mut s, &config, &wind, dt);
        }
        let t = dt * steps as f64;
        let v_expected = [v0[0], v0[1], v0[2] + G0 * t];
        let v = s.velocity_ned();
        assert!(
            norm3(sub3(v, v_expected)) < 5.0e-2,
            "v_ned = {v:?}, expected {v_expected:?}"
        );
        let p_expected = [
            p0[0] + v0[0] * t,
            p0[1] + v0[1] * t,
            p0[2] + v0[2] * t + 0.5 * G0 * t * t,
        ];
        assert!(
            norm3(sub3(s.position, p_expected)) < 0.1,
            "p = {:?}, expected {p_expected:?}",
            s.position
        );
        // The airframe really did rotate (about 0.5 rad), and the speed follows
        // the inertial dynamics.
        assert!((s.angular_velocity[1] - 0.5).abs() < 1e-12);
        let speed_expected = norm3(v_expected);
        assert!((norm3(s.velocity_body) - speed_expected).abs() < 5.0e-2);
    }

    #[test]
    fn steady_wind_shifts_the_ground_motion_by_the_wind_vector() {
        // Same air-relative flight, once in calm air at 60 m/s and once in a
        // 10 m/s headwind at 50 m/s over the ground: every aerodynamic input is
        // identical, so the ground velocity must stay 10 m/s lower (the wind
        // vector) throughout, and the ground track falls behind by 10 m/s * t.
        let config = AircraftConfig::generic_ga();
        let dt = 0.005;
        let steps = 400;
        let wind_ned = [-10.0, 0.0, 0.0];
        let fly = |wind: &WindModel, ground_speed: f64| {
            let mut s = flying_north(ground_speed);
            s.controls = ControlSurfaces::primary(0.01, -0.03, 0.0, 0.4);
            for _ in 0..steps {
                integrate_aircraft(&mut s, &config, wind, dt);
            }
            s
        };
        let calm_wind = WindModel::calm();
        let head_wind = WindModel::steady(wind_ned);
        let calm = fly(&calm_wind, 60.0);
        let windy = fly(&head_wind, 50.0);
        let t = dt * steps as f64;
        let velocity_shift = sub3(windy.velocity_ned(), calm.velocity_ned());
        assert!(
            norm3(sub3(velocity_shift, wind_ned)) < 0.1,
            "ground velocity shift {velocity_shift:?}, expected {wind_ned:?}"
        );
        let position_shift = sub3(windy.position, calm.position);
        assert!(
            norm3(sub3(position_shift, scale3(wind_ned, t))) < 0.2,
            "ground track shift {position_shift:?}, expected {:?}",
            scale3(wind_ned, t)
        );
        let airspeed_calm = calm.air_data(&calm_wind).airspeed;
        let airspeed_windy = windy.air_data(&head_wind).airspeed;
        assert!(
            (airspeed_calm - airspeed_windy).abs() < 0.1,
            "airspeeds {airspeed_calm} vs {airspeed_windy}"
        );
    }

    #[test]
    fn headwind_increases_the_drag_deceleration_in_the_integrator() {
        let config = AircraftConfig::generic_ga();
        let dt = 0.01;
        let u_after = |wind: &WindModel| {
            let mut s = flying_north(60.0);
            integrate_aircraft(&mut s, &config, wind, dt);
            s.velocity_body[0]
        };
        let u_calm = u_after(&WindModel::calm());
        let u_headwind = u_after(&WindModel::steady([-20.0, 0.0, 0.0]));
        let u_tailwind = u_after(&WindModel::steady([20.0, 0.0, 0.0]));
        assert!(u_calm < 60.0, "drag must decelerate the aircraft");
        assert!(
            u_headwind < u_calm,
            "headwind: {u_headwind} vs calm {u_calm}"
        );
        assert!(
            u_tailwind > u_calm,
            "tailwind: {u_tailwind} vs calm {u_calm}"
        );
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
