// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Orbital mechanics for rigid body dynamics.
//!
//! Covers:
//! - Kepler orbital elements (a, e, i, Ω, ω, ν)
//! - Kepler's equation solver (eccentric anomaly via Newton–Raphson)
//! - Orbital period and vis-viva equation
//! - Hohmann and bi-elliptic transfer orbits
//! - Lambert's problem (Izzo algorithm approximation)
//! - Patched conic approximation and sphere of influence
//! - J2 perturbation (Earth oblateness)
//! - Atmospheric drag model (exponential atmosphere)
//! - Solar radiation pressure
//! - Restricted three-body problem (CR3BP)
//! - Lagrange points L1–L5
//! - Orbital debris collision probability (cube method)

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Gravitational parameter of Earth μ = GM \[m³/s²\].
pub const MU_EARTH: f64 = 3.986_004_418e14;

/// Gravitational parameter of the Sun \[m³/s²\].
pub const MU_SUN: f64 = 1.327_124_400e20;

/// Gravitational parameter of the Moon \[m³/s²\].
pub const MU_MOON: f64 = 4.904_869_5e12;

/// Mean radius of the Earth \[m\].
pub const R_EARTH: f64 = 6.371_000e6;

/// Mean radius of the Moon \[m\].
pub const R_MOON: f64 = 1.737_400e6;

/// Mean radius of the Sun \[m\].
pub const R_SUN: f64 = 6.957_000e8;

/// Earth's second zonal harmonic coefficient J2 (oblateness).
pub const J2_EARTH: f64 = 1.082_626_68e-3;

/// Earth equatorial radius (used for J2 calculations) \[m\].
pub const R_EARTH_EQ: f64 = 6.378_137e6;

/// Solar radiation pressure at 1 AU \[N/m²\].
pub const SOLAR_PRESSURE_1AU: f64 = 4.56e-6;

/// Speed of light in vacuum \[m/s\].
pub const SPEED_OF_LIGHT: f64 = 2.997_924_58e8;

/// 1 Astronomical Unit \[m\].
pub const AU: f64 = 1.495_978_707e11;

/// Reference atmospheric density at sea level \[kg/m³\].
pub const RHO_SEA_LEVEL: f64 = 1.225;

/// Atmospheric scale height for exponential model \[m\].
pub const SCALE_HEIGHT_ATM: f64 = 8_500.0;

// ---------------------------------------------------------------------------
// Vec3 helper functions (plain [f64; 3], no nalgebra)
// ---------------------------------------------------------------------------

/// Compute the dot product of two 3-vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Compute the cross product of two 3-vectors.
#[inline]
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Compute the Euclidean norm of a 3-vector.
#[inline]
pub fn norm3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// Normalize a 3-vector to unit length. Returns zero vector if near zero.
#[inline]
pub fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n < 1e-300 {
        [0.0; 3]
    } else {
        [a[0] / n, a[1] / n, a[2] / n]
    }
}

/// Add two 3-vectors.
#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors.
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector by a scalar.
#[inline]
pub fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

// ---------------------------------------------------------------------------
// Kepler Orbital Elements
// ---------------------------------------------------------------------------

/// Classical Keplerian orbital elements describing a conic section orbit.
///
/// The six elements uniquely define the shape, size, orientation, and current
/// position of a body in a two-body gravitational field.
#[derive(Debug, Clone, PartialEq)]
pub struct KeplerElements {
    /// Semi-major axis \[m\]. For a hyperbolic orbit, `a < 0`.
    pub semi_major_axis: f64,
    /// Orbital eccentricity (dimensionless). 0 = circular, 0<e<1 = elliptic,
    /// e=1 = parabolic, e>1 = hyperbolic.
    pub eccentricity: f64,
    /// Inclination \[rad\]. Angle between the orbital plane and the reference plane.
    pub inclination: f64,
    /// Right ascension of ascending node (RAAN) Ω \[rad\].
    pub raan: f64,
    /// Argument of periapsis ω \[rad\].
    pub arg_periapsis: f64,
    /// True anomaly ν \[rad\]. Current angular position of the body in its orbit.
    pub true_anomaly: f64,
}

impl KeplerElements {
    /// Construct a new set of Keplerian elements.
    pub fn new(
        semi_major_axis: f64,
        eccentricity: f64,
        inclination: f64,
        raan: f64,
        arg_periapsis: f64,
        true_anomaly: f64,
    ) -> Self {
        Self {
            semi_major_axis,
            eccentricity,
            inclination,
            raan,
            arg_periapsis,
            true_anomaly,
        }
    }

    /// Return the periapsis radius \[m\].
    pub fn periapsis_radius(&self) -> f64 {
        self.semi_major_axis * (1.0 - self.eccentricity)
    }

    /// Return the apoapsis radius \[m\]. Returns `f64::INFINITY` for parabolic/hyperbolic orbits.
    pub fn apoapsis_radius(&self) -> f64 {
        if self.eccentricity >= 1.0 {
            f64::INFINITY
        } else {
            self.semi_major_axis * (1.0 + self.eccentricity)
        }
    }

    /// Return the semi-latus rectum p = a(1 - e²) \[m\].
    pub fn semi_latus_rectum(&self) -> f64 {
        self.semi_major_axis * (1.0 - self.eccentricity * self.eccentricity)
    }

    /// Return the orbital period \[s\] for an elliptic orbit under `mu` \[m³/s²\].
    /// Returns `f64::INFINITY` for non-elliptic (e ≥ 1) orbits.
    pub fn period(&self, mu: f64) -> f64 {
        if self.eccentricity >= 1.0 || self.semi_major_axis <= 0.0 {
            f64::INFINITY
        } else {
            orbital_period(self.semi_major_axis, mu)
        }
    }

    /// Convert Keplerian elements to Cartesian state vector `([rx,ry,rz], [vx,vy,vz])`.
    ///
    /// Returns the position \[m\] and velocity \[m/s\] in the Equatorial Inertial frame.
    pub fn to_state_vector(&self, mu: f64) -> ([f64; 3], [f64; 3]) {
        kepler_to_cartesian(self, mu)
    }

    /// Check whether the orbit is circular (e < 1e-6).
    pub fn is_circular(&self) -> bool {
        self.eccentricity < 1e-6
    }

    /// Check whether the orbit is elliptic (0 ≤ e < 1).
    pub fn is_elliptic(&self) -> bool {
        self.eccentricity < 1.0
    }

    /// Check whether the orbit is hyperbolic (e > 1).
    pub fn is_hyperbolic(&self) -> bool {
        self.eccentricity > 1.0
    }
}

// ---------------------------------------------------------------------------
// Kepler's equation solver
// ---------------------------------------------------------------------------

/// Solve Kepler's equation M = E - e·sin(E) for eccentric anomaly E \[rad\].
///
/// Uses Newton–Raphson iteration. Converges typically in 3–6 iterations.
///
/// # Arguments
/// * `mean_anomaly` — Mean anomaly M \[rad\].
/// * `eccentricity` — Orbital eccentricity (0 ≤ e < 1).
/// * `tol` — Convergence tolerance (e.g. 1e-12).
/// * `max_iter` — Maximum number of iterations.
///
/// # Returns
/// Eccentric anomaly E \[rad\].
pub fn solve_kepler_equation(
    mean_anomaly: f64,
    eccentricity: f64,
    tol: f64,
    max_iter: usize,
) -> f64 {
    // Initial guess
    let m = mean_anomaly.rem_euclid(2.0 * PI);
    let mut e = if eccentricity < 0.8 { m } else { PI };

    for _ in 0..max_iter {
        let f = e - eccentricity * e.sin() - m;
        let df = 1.0 - eccentricity * e.cos();
        let delta = f / df;
        e -= delta;
        if delta.abs() < tol {
            break;
        }
    }
    e
}

/// Convert mean anomaly M to true anomaly ν for an elliptic orbit.
///
/// # Arguments
/// * `mean_anomaly` — Mean anomaly \[rad\].
/// * `eccentricity` — Eccentricity (0 ≤ e < 1).
///
/// # Returns
/// True anomaly ν \[rad\] in `[0, 2π)`.
pub fn mean_to_true_anomaly(mean_anomaly: f64, eccentricity: f64) -> f64 {
    let ecc_anomaly = solve_kepler_equation(mean_anomaly, eccentricity, 1e-12, 100);
    eccentric_to_true_anomaly(ecc_anomaly, eccentricity)
}

/// Convert eccentric anomaly E to true anomaly ν.
///
/// # Arguments
/// * `ecc_anomaly` — Eccentric anomaly E \[rad\].
/// * `eccentricity` — Eccentricity.
///
/// # Returns
/// True anomaly ν \[rad\].
pub fn eccentric_to_true_anomaly(ecc_anomaly: f64, eccentricity: f64) -> f64 {
    let e = eccentricity;
    let half_e = ecc_anomaly / 2.0;
    let sqrt_factor = ((1.0 + e) / (1.0 - e)).sqrt();
    2.0 * (sqrt_factor * half_e.tan()).atan()
}

/// Convert true anomaly ν to mean anomaly M for an elliptic orbit.
///
/// # Arguments
/// * `true_anomaly` — True anomaly ν \[rad\].
/// * `eccentricity` — Eccentricity (0 ≤ e < 1).
///
/// # Returns
/// Mean anomaly M \[rad\] in `[0, 2π)`.
pub fn true_to_mean_anomaly(true_anomaly: f64, eccentricity: f64) -> f64 {
    let e = eccentricity;
    let nu = true_anomaly;
    let half_nu = nu / 2.0;
    let sqrt_factor = ((1.0 - e) / (1.0 + e)).sqrt();
    let ecc_anomaly = 2.0 * (sqrt_factor * half_nu.tan()).atan();
    let m = ecc_anomaly - e * ecc_anomaly.sin();
    m.rem_euclid(2.0 * PI)
}

// ---------------------------------------------------------------------------
// Orbital period and vis-viva
// ---------------------------------------------------------------------------

/// Compute the orbital period T = 2π √(a³/μ) \[s\].
///
/// # Arguments
/// * `semi_major_axis` — Semi-major axis a \[m\].
/// * `mu` — Gravitational parameter μ \[m³/s²\].
pub fn orbital_period(semi_major_axis: f64, mu: f64) -> f64 {
    2.0 * PI * (semi_major_axis.powi(3) / mu).sqrt()
}

/// Compute orbital velocity using the vis-viva equation: v² = μ(2/r - 1/a).
///
/// # Arguments
/// * `radius` — Current distance from the central body \[m\].
/// * `semi_major_axis` — Semi-major axis a \[m\].
/// * `mu` — Gravitational parameter μ \[m³/s²\].
///
/// # Returns
/// Orbital speed \[m/s\].
pub fn vis_viva_speed(radius: f64, semi_major_axis: f64, mu: f64) -> f64 {
    (mu * (2.0 / radius - 1.0 / semi_major_axis)).sqrt()
}

/// Compute circular orbital speed at a given radius: v_c = √(μ/r).
///
/// # Arguments
/// * `radius` — Orbital radius \[m\].
/// * `mu` — Gravitational parameter \[m³/s²\].
pub fn circular_speed(radius: f64, mu: f64) -> f64 {
    (mu / radius).sqrt()
}

/// Compute escape speed from a given radius: v_esc = √(2μ/r).
///
/// # Arguments
/// * `radius` — Distance from the central body \[m\].
/// * `mu` — Gravitational parameter \[m³/s²\].
pub fn escape_speed(radius: f64, mu: f64) -> f64 {
    (2.0 * mu / radius).sqrt()
}

// ---------------------------------------------------------------------------
// Kepler elements ↔ Cartesian state vector
// ---------------------------------------------------------------------------

/// Convert Keplerian elements to Cartesian state vector in the equatorial frame.
///
/// # Returns
/// `(position [m; 3], velocity [m/s; 3])`.
pub fn kepler_to_cartesian(elem: &KeplerElements, mu: f64) -> ([f64; 3], [f64; 3]) {
    let a = elem.semi_major_axis;
    let e = elem.eccentricity;
    let i = elem.inclination;
    let omega_big = elem.raan;
    let omega = elem.arg_periapsis;
    let nu = elem.true_anomaly;

    let p = a * (1.0 - e * e);
    let r_mag = p / (1.0 + e * nu.cos());

    // Position in perifocal frame
    let r_pf = [r_mag * nu.cos(), r_mag * nu.sin(), 0.0];
    // Velocity in perifocal frame
    let sqrt_mu_over_p = (mu / p).sqrt();
    let v_pf = [
        -sqrt_mu_over_p * nu.sin(),
        sqrt_mu_over_p * (e + nu.cos()),
        0.0,
    ];

    // Rotation matrix: perifocal → ECI
    let r = rotation_pf_to_eci(omega_big, i, omega);
    let pos = mat3_vec3_mul(r, r_pf);
    let vel = mat3_vec3_mul(r, v_pf);
    (pos, vel)
}

/// Convert a Cartesian state vector in ECI to Keplerian elements.
///
/// # Arguments
/// * `pos` — Position vector \[m\].
/// * `vel` — Velocity vector \[m/s\].
/// * `mu` — Gravitational parameter \[m³/s²\].
pub fn cartesian_to_kepler(pos: [f64; 3], vel: [f64; 3], mu: f64) -> KeplerElements {
    let r_mag = norm3(pos);
    let v_mag = norm3(vel);

    // Angular momentum vector
    let h = cross3(pos, vel);
    let h_mag = norm3(h);

    // Node vector (ascending node direction)
    let z_hat = [0.0, 0.0, 1.0];
    let n_vec = cross3(z_hat, h);
    let n_mag = norm3(n_vec);

    // Eccentricity vector
    let term1 = scale3(pos, (v_mag * v_mag - mu / r_mag) / mu);
    let term2 = scale3(vel, dot3(pos, vel) / mu);
    let e_vec = sub3(term1, term2);
    let eccentricity = norm3(e_vec);

    // Semi-major axis
    let energy = v_mag * v_mag / 2.0 - mu / r_mag;
    let semi_major_axis = -mu / (2.0 * energy);

    // Inclination
    let inclination = (h[2] / h_mag).acos();

    // RAAN
    let raan = if n_mag < 1e-10 {
        0.0
    } else if n_vec[1] >= 0.0 {
        (n_vec[0] / n_mag).acos()
    } else {
        2.0 * PI - (n_vec[0] / n_mag).acos()
    };

    // Argument of periapsis
    let arg_periapsis = if eccentricity < 1e-10 || n_mag < 1e-10 {
        0.0
    } else {
        let cos_w = dot3(n_vec, e_vec) / (n_mag * eccentricity);
        let cos_w = cos_w.clamp(-1.0, 1.0);
        if e_vec[2] >= 0.0 {
            cos_w.acos()
        } else {
            2.0 * PI - cos_w.acos()
        }
    };

    // True anomaly
    let true_anomaly = if eccentricity < 1e-10 {
        let cos_nu = dot3(pos, vel) / (r_mag * norm3(vel));
        cos_nu.clamp(-1.0, 1.0).acos()
    } else {
        let cos_nu = dot3(e_vec, pos) / (eccentricity * r_mag);
        let cos_nu = cos_nu.clamp(-1.0, 1.0);
        if dot3(pos, vel) >= 0.0 {
            cos_nu.acos()
        } else {
            2.0 * PI - cos_nu.acos()
        }
    };

    KeplerElements {
        semi_major_axis,
        eccentricity,
        inclination,
        raan,
        arg_periapsis,
        true_anomaly,
    }
}

// ---------------------------------------------------------------------------
// Rotation matrix helpers
// ---------------------------------------------------------------------------

/// Build the 3×3 rotation matrix from perifocal to ECI frame.
///
/// Row-major 3×3 stored as \[f64; 9\] (index = row*3 + col).
fn rotation_pf_to_eci(raan: f64, inc: f64, arg_peri: f64) -> [f64; 9] {
    let cos_o = raan.cos();
    let sin_o = raan.sin();
    let cos_i = inc.cos();
    let sin_i = inc.sin();
    let cos_w = arg_peri.cos();
    let sin_w = arg_peri.sin();

    [
        cos_o * cos_w - sin_o * sin_w * cos_i,
        -cos_o * sin_w - sin_o * cos_w * cos_i,
        sin_o * sin_i,
        sin_o * cos_w + cos_o * sin_w * cos_i,
        -sin_o * sin_w + cos_o * cos_w * cos_i,
        -cos_o * sin_i,
        sin_w * sin_i,
        cos_w * sin_i,
        cos_i,
    ]
}

/// Multiply a row-major 3×3 matrix by a 3-vector.
fn mat3_vec3_mul(m: [f64; 9], v: [f64; 3]) -> [f64; 3] {
    [
        m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
        m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
        m[6] * v[0] + m[7] * v[1] + m[8] * v[2],
    ]
}

// ---------------------------------------------------------------------------
// Hohmann transfer orbit
// ---------------------------------------------------------------------------

/// Result of a Hohmann transfer calculation.
#[derive(Debug, Clone, PartialEq)]
pub struct HohmannTransfer {
    /// Semi-major axis of the transfer ellipse \[m\].
    pub transfer_sma: f64,
    /// First delta-V burn (at departure orbit) \[m/s\].
    pub delta_v1: f64,
    /// Second delta-V burn (at arrival orbit) \[m/s\].
    pub delta_v2: f64,
    /// Total delta-V \[m/s\].
    pub delta_v_total: f64,
    /// Transfer time (half the period of the transfer ellipse) \[s\].
    pub transfer_time: f64,
}

/// Compute a Hohmann transfer between two circular coplanar orbits.
///
/// # Arguments
/// * `r1` — Radius of departure orbit \[m\].
/// * `r2` — Radius of arrival orbit \[m\].
/// * `mu` — Gravitational parameter \[m³/s²\].
pub fn hohmann_transfer(r1: f64, r2: f64, mu: f64) -> HohmannTransfer {
    let a_transfer = (r1 + r2) / 2.0;

    let v1_circ = circular_speed(r1, mu);
    let v2_circ = circular_speed(r2, mu);

    let v_transfer_peri = vis_viva_speed(r1, a_transfer, mu);
    let v_transfer_apo = vis_viva_speed(r2, a_transfer, mu);

    let dv1 = (v_transfer_peri - v1_circ).abs();
    let dv2 = (v2_circ - v_transfer_apo).abs();
    let t_transfer = orbital_period(a_transfer, mu) / 2.0;

    HohmannTransfer {
        transfer_sma: a_transfer,
        delta_v1: dv1,
        delta_v2: dv2,
        delta_v_total: dv1 + dv2,
        transfer_time: t_transfer,
    }
}

// ---------------------------------------------------------------------------
// Bi-elliptic transfer orbit
// ---------------------------------------------------------------------------

/// Result of a bi-elliptic transfer calculation.
#[derive(Debug, Clone, PartialEq)]
pub struct BiEllipticTransfer {
    /// First transfer ellipse semi-major axis \[m\].
    pub sma1: f64,
    /// Second transfer ellipse semi-major axis \[m\].
    pub sma2: f64,
    /// First delta-V at r1 \[m/s\].
    pub delta_v1: f64,
    /// Second delta-V at intermediate radius rb \[m/s\].
    pub delta_v2: f64,
    /// Third delta-V at r2 \[m/s\].
    pub delta_v3: f64,
    /// Total delta-V \[m/s\].
    pub delta_v_total: f64,
    /// Total transfer time \[s\].
    pub transfer_time: f64,
}

/// Compute a bi-elliptic transfer between two circular coplanar orbits via an
/// intermediate apoapsis at `rb`.
///
/// # Arguments
/// * `r1` — Departure orbit radius \[m\].
/// * `r2` — Arrival orbit radius \[m\].
/// * `rb` — Intermediate apoapsis radius (must be > max(r1, r2)) \[m\].
/// * `mu` — Gravitational parameter \[m³/s²\].
pub fn bi_elliptic_transfer(r1: f64, r2: f64, rb: f64, mu: f64) -> BiEllipticTransfer {
    let a1 = (r1 + rb) / 2.0;
    let a2 = (rb + r2) / 2.0;

    let v1 = circular_speed(r1, mu);
    let v2 = circular_speed(r2, mu);

    let v_a1_peri = vis_viva_speed(r1, a1, mu);
    let v_a1_apo = vis_viva_speed(rb, a1, mu);
    let v_a2_apo = vis_viva_speed(rb, a2, mu);
    let v_a2_peri = vis_viva_speed(r2, a2, mu);

    let dv1 = (v_a1_peri - v1).abs();
    let dv2 = (v_a2_apo - v_a1_apo).abs();
    let dv3 = (v2 - v_a2_peri).abs();

    let t1 = orbital_period(a1, mu) / 2.0;
    let t2 = orbital_period(a2, mu) / 2.0;

    BiEllipticTransfer {
        sma1: a1,
        sma2: a2,
        delta_v1: dv1,
        delta_v2: dv2,
        delta_v3: dv3,
        delta_v_total: dv1 + dv2 + dv3,
        transfer_time: t1 + t2,
    }
}

// ---------------------------------------------------------------------------
// Lambert's problem
// ---------------------------------------------------------------------------

/// Result of Lambert's problem solver.
#[derive(Debug, Clone, PartialEq)]
pub struct LambertSolution {
    /// Departure velocity vector \[m/s\].
    pub v1: [f64; 3],
    /// Arrival velocity vector \[m/s\].
    pub v2: [f64; 3],
    /// Semi-major axis of the transfer conic \[m\].
    pub semi_major_axis: f64,
}

/// Solve Lambert's problem (Izzo-style universal variable approach).
///
/// Given two position vectors and a transfer time, find the connecting
/// orbit. This uses the Lancaster–Blanchard/universal variable method.
///
/// # Arguments
/// * `r1` — Departure position \[m\].
/// * `r2` — Arrival position \[m\].
/// * `tof` — Time of flight \[s\].
/// * `mu` — Gravitational parameter \[m³/s²\].
/// * `prograde` — If `true`, solve the prograde (short-way) transfer.
///
/// # Returns
/// `None` if convergence fails.
pub fn lambert_problem(
    r1: [f64; 3],
    r2: [f64; 3],
    tof: f64,
    mu: f64,
    prograde: bool,
) -> Option<LambertSolution> {
    let r1_mag = norm3(r1);
    let r2_mag = norm3(r2);

    // Angle between r1 and r2
    let cos_dnu = (dot3(r1, r2) / (r1_mag * r2_mag)).clamp(-1.0, 1.0);
    let cross_r = cross3(r1, r2);
    let sin_dnu_sign = if prograde { cross_r[2] } else { -cross_r[2] };
    let sin_dnu = (1.0 - cos_dnu * cos_dnu).sqrt() * sin_dnu_sign.signum();

    let s = (r1_mag
        + r2_mag
        + (r1_mag * r1_mag + r2_mag * r2_mag - 2.0 * r1_mag * r2_mag * cos_dnu).sqrt())
        / 2.0;

    let lambda_sq = 1.0 - (r1_mag * r2_mag * (1.0 + cos_dnu)) / (s * s);
    let lambda = if sin_dnu >= 0.0 {
        lambda_sq.sqrt()
    } else {
        -lambda_sq.sqrt()
    };

    let t_norm = tof * (2.0 * mu / (s * s * s)).sqrt();

    // Halley's method on the x variable
    let x = lambert_halley(lambda, t_norm)?;

    let gamma = (mu * s / 2.0).sqrt();
    let rho = (r1_mag - r2_mag) / (2.0 * s - r1_mag - r2_mag + s);
    let sigma = (1.0 - rho * rho).sqrt();

    let y = (1.0 - lambda * lambda * (1.0 - x * x)).sqrt();

    let vr1 = gamma * ((lambda * y - x) - rho * (lambda * y + x)) / r1_mag;
    let vr2 = -gamma * ((lambda * y - x) + rho * (lambda * y + x)) / r2_mag;

    let vt1 = gamma * sigma * (y + lambda * x) / r1_mag;
    let vt2 = gamma * sigma * (y + lambda * x) / r2_mag;

    let r1_hat = normalize3(r1);
    let r2_hat = normalize3(r2);
    let h_hat = normalize3(cross3(r1, r2));
    let t1_hat = normalize3(cross3(h_hat, r1_hat));
    let t2_hat = normalize3(cross3(h_hat, r2_hat));

    let v1 = add3(scale3(r1_hat, vr1), scale3(t1_hat, vt1));
    let v2 = add3(scale3(r2_hat, vr2), scale3(t2_hat, vt2));

    // Estimate semi-major axis from energy
    let energy1 = norm3(v1).powi(2) / 2.0 - mu / r1_mag;
    let semi_major_axis = if energy1.abs() < 1e-10 {
        f64::INFINITY
    } else {
        -mu / (2.0 * energy1)
    };

    Some(LambertSolution {
        v1,
        v2,
        semi_major_axis,
    })
}

/// Halley's iteration for Lambert's universal variable x.
fn lambert_halley(lambda: f64, t_norm: f64) -> Option<f64> {
    let mut x = 0.0_f64;
    for _ in 0..50 {
        let a = 1.0 / (1.0 - x * x);
        let xi = a.sqrt();
        let eta = (lambda + x * xi).acos();
        let s1 = (1.0 - lambda - x * xi) / 2.0;
        let q = 4.0 / 3.0 * hypergeom_2f1b(3.0 / 2.0, 1.0, 5.0 / 2.0, s1);
        let t_x = (xi * (eta + lambda * (1.0 - x * x).sqrt()) + q * s1.powi(2)) / (1.0 - x * x);

        // Derivative
        let dt = (3.0 * t_x * x - 2.0 + 2.0 * lambda.powi(3) * x / xi) / (1.0 - x * x);

        let delta = (t_x - t_norm) / dt;
        x -= delta;
        if delta.abs() < 1e-10 {
            return Some(x);
        }
    }
    None
}

/// Approximate hypergeometric function ₂F₁(3/2, 1; 5/2; z) via series.
fn hypergeom_2f1b(_a: f64, _b: f64, _c: f64, z: f64) -> f64 {
    // Series truncated to 20 terms: ₂F₁(a,b;c;z) = Σ (a)_n (b)_n / ((c)_n n!) z^n
    let a = 1.5_f64;
    let b = 1.0_f64;
    let c = 2.5_f64;
    let mut sum = 1.0;
    let mut term = 1.0;
    for n in 1..20usize {
        let n_f = n as f64;
        term *= (a + n_f - 1.0) * (b + n_f - 1.0) / ((c + n_f - 1.0) * n_f) * z;
        sum += term;
        if term.abs() < 1e-14 {
            break;
        }
    }
    sum
}

// ---------------------------------------------------------------------------
// Sphere of influence
// ---------------------------------------------------------------------------

/// Compute the sphere of influence (SOI) radius of a secondary body
/// orbiting a primary.
///
/// r_SOI ≈ a_body · (m_body / m_primary)^(2/5)
///
/// # Arguments
/// * `semi_major_axis` — Semi-major axis of the secondary's orbit around the primary \[m\].
/// * `mu_secondary` — Gravitational parameter of the secondary body \[m³/s²\].
/// * `mu_primary` — Gravitational parameter of the primary body \[m³/s²\].
pub fn sphere_of_influence(semi_major_axis: f64, mu_secondary: f64, mu_primary: f64) -> f64 {
    semi_major_axis * (mu_secondary / mu_primary).powf(2.0 / 5.0)
}

// ---------------------------------------------------------------------------
// Patched conic approximation
// ---------------------------------------------------------------------------

/// Compute the hyperbolic excess speed at SOI entry given departure conditions.
///
/// v_∞ = √(v_depart² - v_esc_SOI²)
///
/// # Arguments
/// * `v_departure` — Speed of the spacecraft relative to the central body at SOI \[m/s\].
/// * `mu_secondary` — Gravitational parameter of the destination body \[m³/s²\].
/// * `r_soi` — Sphere of influence radius \[m\].
pub fn hyperbolic_excess_speed(v_departure: f64, mu_secondary: f64, r_soi: f64) -> f64 {
    let v_esc = escape_speed(r_soi, mu_secondary);
    let v_inf_sq = v_departure * v_departure - v_esc * v_esc;
    if v_inf_sq < 0.0 { 0.0 } else { v_inf_sq.sqrt() }
}

/// Compute the periapsis altitude for a hyperbolic flyby given v_∞ and periapsis radius.
///
/// Returns the speed at periapsis: v_p = √(v_∞² + 2μ/r_p).
///
/// # Arguments
/// * `v_inf` — Hyperbolic excess speed \[m/s\].
/// * `periapsis_radius` — Periapsis distance \[m\].
/// * `mu` — Gravitational parameter \[m³/s²\].
pub fn hyperbolic_periapsis_speed(v_inf: f64, periapsis_radius: f64, mu: f64) -> f64 {
    (v_inf * v_inf + 2.0 * mu / periapsis_radius).sqrt()
}

// ---------------------------------------------------------------------------
// J2 perturbation (Earth oblateness)
// ---------------------------------------------------------------------------

/// Compute the J2 perturbation acceleration on a satellite.
///
/// The J2 term arises from Earth's equatorial bulge and is the dominant
/// non-spherical gravitational perturbation.
///
/// # Arguments
/// * `pos` — Satellite position in ECI \[m\].
/// * `mu` — Earth gravitational parameter \[m³/s²\].
/// * `j2` — J2 coefficient (use `J2_EARTH` ≈ 1.08263e-3).
/// * `r_eq` — Earth equatorial radius \[m\].
///
/// # Returns
/// Perturbation acceleration \[m/s²\].
pub fn j2_acceleration(pos: [f64; 3], mu: f64, j2: f64, r_eq: f64) -> [f64; 3] {
    let r = norm3(pos);
    let r2 = r * r;
    let r5 = r2 * r2 * r;
    let z2 = pos[2] * pos[2];
    let factor = -1.5 * j2 * mu * r_eq * r_eq / r5;

    [
        factor * pos[0] * (1.0 - 5.0 * z2 / r2),
        factor * pos[1] * (1.0 - 5.0 * z2 / r2),
        factor * pos[2] * (3.0 - 5.0 * z2 / r2),
    ]
}

/// Compute the secular rate of RAAN precession due to J2 \[rad/s\].
///
/// dΩ/dt = -3/2 · n · J2 · (R_eq/a)² · cos(i) / (1-e²)²
///
/// # Arguments
/// * `semi_major_axis` — Semi-major axis \[m\].
/// * `eccentricity` — Orbital eccentricity.
/// * `inclination` — Inclination \[rad\].
/// * `mu` — Gravitational parameter \[m³/s²\].
/// * `j2` — J2 coefficient.
/// * `r_eq` — Equatorial radius \[m\].
pub fn j2_raan_precession_rate(
    semi_major_axis: f64,
    eccentricity: f64,
    inclination: f64,
    mu: f64,
    j2: f64,
    r_eq: f64,
) -> f64 {
    let n = (mu / semi_major_axis.powi(3)).sqrt();
    let e2 = 1.0 - eccentricity * eccentricity;
    -1.5 * n * j2 * (r_eq / semi_major_axis).powi(2) * inclination.cos() / (e2 * e2)
}

/// Compute the secular rate of argument of periapsis precession due to J2 \[rad/s\].
///
/// dω/dt = 3/4 · n · J2 · (R_eq/a)² · (5cos²i - 1) / (1-e²)²
pub fn j2_arg_periapsis_precession_rate(
    semi_major_axis: f64,
    eccentricity: f64,
    inclination: f64,
    mu: f64,
    j2: f64,
    r_eq: f64,
) -> f64 {
    let n = (mu / semi_major_axis.powi(3)).sqrt();
    let e2 = 1.0 - eccentricity * eccentricity;
    let cos_i = inclination.cos();
    0.75 * n * j2 * (r_eq / semi_major_axis).powi(2) * (5.0 * cos_i * cos_i - 1.0) / (e2 * e2)
}

// ---------------------------------------------------------------------------
// Atmospheric drag model
// ---------------------------------------------------------------------------

/// Compute atmospheric density via exponential scale-height model.
///
/// ρ(h) = ρ₀ · exp(-(h - h_ref) / H)
///
/// # Arguments
/// * `altitude` — Altitude above Earth's surface \[m\].
/// * `rho0` — Reference density \[kg/m³\].
/// * `h_ref` — Reference altitude \[m\].
/// * `scale_height` — Atmospheric scale height H \[m\].
pub fn atmospheric_density(altitude: f64, rho0: f64, h_ref: f64, scale_height: f64) -> f64 {
    rho0 * (-(altitude - h_ref) / scale_height).exp()
}

/// Compute atmospheric drag acceleration on a spacecraft.
///
/// a_drag = -1/2 · ρ · (Cd · A / m) · v² · v̂
///
/// # Arguments
/// * `velocity` — Spacecraft velocity in ECI \[m/s\].
/// * `altitude` — Altitude \[m\].
/// * `cd` — Drag coefficient (typically 2.2 for satellites).
/// * `area_over_mass` — Effective area / mass ratio \[m²/kg\].
/// * `rho0` — Reference density \[kg/m³\].
/// * `h_ref` — Reference altitude \[m\].
/// * `scale_height` — Scale height \[m\].
pub fn atmospheric_drag_acceleration(
    velocity: [f64; 3],
    altitude: f64,
    cd: f64,
    area_over_mass: f64,
    rho0: f64,
    h_ref: f64,
    scale_height: f64,
) -> [f64; 3] {
    let rho = atmospheric_density(altitude, rho0, h_ref, scale_height);
    let v_mag = norm3(velocity);
    if v_mag < 1e-10 {
        return [0.0; 3];
    }
    let factor = -0.5 * rho * cd * area_over_mass * v_mag;
    scale3(velocity, factor)
}

// ---------------------------------------------------------------------------
// Solar radiation pressure
// ---------------------------------------------------------------------------

/// Compute solar radiation pressure acceleration.
///
/// a_srp = P_sun · (C_r · A / m) · r̂_sun→sat
///
/// # Arguments
/// * `sun_to_sat` — Direction from the Sun to the satellite (normalized) \[m\].
/// * `distance_from_sun` — Distance from Sun to satellite \[m\].
/// * `cr` — Radiation pressure coefficient (1.0 for absorbed, 2.0 for reflected).
/// * `area_over_mass` — Effective area / mass \[m²/kg\].
pub fn solar_radiation_pressure_acceleration(
    sun_to_sat: [f64; 3],
    distance_from_sun: f64,
    cr: f64,
    area_over_mass: f64,
) -> [f64; 3] {
    let p_sun = SOLAR_PRESSURE_1AU * (AU / distance_from_sun).powi(2);
    let hat = normalize3(sun_to_sat);
    scale3(hat, p_sun * cr * area_over_mass)
}

// ---------------------------------------------------------------------------
// Restricted three-body problem (CR3BP)
// ---------------------------------------------------------------------------

/// Compute the CR3BP (Circular Restricted Three-Body Problem) equations of motion.
///
/// In the rotating frame with primaries at x = -μ* and x = 1-μ*, where
/// μ* = m2/(m1+m2) is the mass ratio.
///
/// # Arguments
/// * `state` — State vector `[x, y, z, xdot, ydot, zdot]` in normalized units.
/// * `mass_ratio` — μ* = m_secondary / (m_primary + m_secondary).
///
/// # Returns
/// Time derivative of the state `[vx, vy, vz, ax, ay, az]`.
pub fn cr3bp_equations(state: [f64; 6], mass_ratio: f64) -> [f64; 6] {
    let mu = mass_ratio;
    let [x, y, z, xd, yd, _zd] = state;

    // Distances to primaries
    let d1 = ((x + mu).powi(2) + y * y + z * z).sqrt();
    let d2 = ((x - 1.0 + mu).powi(2) + y * y + z * z).sqrt();

    let ax = 2.0 * yd + x - (1.0 - mu) * (x + mu) / d1.powi(3) - mu * (x - 1.0 + mu) / d2.powi(3);
    let ay = -2.0 * xd + y - (1.0 - mu) * y / d1.powi(3) - mu * y / d2.powi(3);
    let az = -(1.0 - mu) * z / d1.powi(3) - mu * z / d2.powi(3);

    [state[3], state[4], state[5], ax, ay, az]
}

/// Compute the Jacobi integral (energy integral) of the CR3BP.
///
/// C_J = x² + y² + 2(1-μ)/r1 + 2μ/r2 - (ẋ² + ẏ² + ż²)
///
/// # Arguments
/// * `state` — State vector `[x, y, z, xdot, ydot, zdot]`.
/// * `mass_ratio` — Mass parameter μ*.
pub fn jacobi_integral(state: [f64; 6], mass_ratio: f64) -> f64 {
    let mu = mass_ratio;
    let [x, y, z, xd, yd, zd] = state;
    let r1 = ((x + mu).powi(2) + y * y + z * z).sqrt();
    let r2 = ((x - 1.0 + mu).powi(2) + y * y + z * z).sqrt();
    let v_sq = xd * xd + yd * yd + zd * zd;
    x * x + y * y + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2 - v_sq
}

// ---------------------------------------------------------------------------
// Lagrange points
// ---------------------------------------------------------------------------

/// Compute the five Lagrange points in the CR3BP rotating frame.
///
/// Returns `[L1, L2, L3, L4, L5]` as `[x, y, 0]` positions in normalized units.
/// L1–L3 are collinear; L4 and L5 are equilateral triangle points.
///
/// # Arguments
/// * `mass_ratio` — μ* = m_secondary / (m_primary + m_secondary).
pub fn lagrange_points(mass_ratio: f64) -> [[f64; 3]; 5] {
    let mu = mass_ratio;

    // L4 and L5 are exact
    let l4 = [0.5 - mu, 3.0_f64.sqrt() / 2.0, 0.0];
    let l5 = [0.5 - mu, -3.0_f64.sqrt() / 2.0, 0.0];

    // L1: between primaries, x ∈ (1-μ-1, 1-μ) → solve numerically
    let l1_x = find_lagrange_collinear(mu, LagrangeCollinearPoint::L1);
    let l2_x = find_lagrange_collinear(mu, LagrangeCollinearPoint::L2);
    let l3_x = find_lagrange_collinear(mu, LagrangeCollinearPoint::L3);

    [[l1_x, 0.0, 0.0], [l2_x, 0.0, 0.0], [l3_x, 0.0, 0.0], l4, l5]
}

/// Identify which collinear Lagrange point to solve for.
enum LagrangeCollinearPoint {
    /// L1: between the two primaries.
    L1,
    /// L2: beyond the secondary (farther from primary).
    L2,
    /// L3: on the far side of the primary.
    L3,
}

/// Find the x-coordinate of a collinear Lagrange point via Newton–Raphson.
fn find_lagrange_collinear(mu: f64, point: LagrangeCollinearPoint) -> f64 {
    let mut x = match point {
        LagrangeCollinearPoint::L1 => 1.0 - mu - 0.1,
        LagrangeCollinearPoint::L2 => 1.0 - mu + 0.1,
        LagrangeCollinearPoint::L3 => -1.0 - mu,
    };

    for _ in 0..100 {
        let d1 = (x + mu).abs();
        let d2 = (x - 1.0 + mu).abs();
        let s1 = (x + mu).signum();
        let s2 = (x - 1.0 + mu).signum();

        let f = x - (1.0 - mu) * s1 / (d1 * d1) - mu * s2 / (d2 * d2);
        let df = 1.0 + 2.0 * (1.0 - mu) / d1.powi(3) + 2.0 * mu / d2.powi(3);

        let delta = f / df;
        x -= delta;
        if delta.abs() < 1e-12 {
            break;
        }
    }
    x
}

// ---------------------------------------------------------------------------
// Orbital debris collision probability
// ---------------------------------------------------------------------------

/// Compute the probability of collision between two space objects using the
/// "cube method" (Foster & Estes, 1992).
///
/// p_c ≈ (1 / (2π σ²)) · exp(-d²/(2σ²)) · r_c² · |Δv| · dt
///
/// # Arguments
/// * `combined_radius` — Combined hard-body radius r_c = r1 + r2 \[m\].
/// * `relative_velocity` — Relative velocity between objects \[m/s\].
/// * `sigma_pos` — Combined position uncertainty (1-sigma) \[m\].
/// * `dt` — Time step (conjunction duration) \[s\].
///
/// # Returns
/// Approximate collision probability (dimensionless, in `[0, 1]`).
pub fn collision_probability(
    combined_radius: f64,
    relative_velocity: f64,
    sigma_pos: f64,
    dt: f64,
) -> f64 {
    let sigma2 = sigma_pos * sigma_pos;
    let volume = PI.powf(1.5) * combined_radius * combined_radius * relative_velocity * dt;
    let prob_density = 1.0 / ((2.0 * PI).powf(1.5) * sigma2.powf(1.5));
    let p = prob_density * volume;
    p.min(1.0)
}

/// Compute the time to closest approach (TCA) given a linear conjunction geometry.
///
/// # Arguments
/// * `r_rel` — Relative position at epoch \[m\].
/// * `v_rel` — Relative velocity \[m/s\].
///
/// # Returns
/// Time to closest approach \[s\]. Negative means closest approach was in the past.
pub fn time_to_closest_approach(r_rel: [f64; 3], v_rel: [f64; 3]) -> f64 {
    let v2 = dot3(v_rel, v_rel);
    if v2 < 1e-20 {
        return 0.0;
    }
    -dot3(r_rel, v_rel) / v2
}

/// Compute miss distance at closest approach.
///
/// # Arguments
/// * `r_rel` — Relative position at epoch \[m\].
/// * `v_rel` — Relative velocity \[m/s\].
///
/// # Returns
/// Miss distance at TCA \[m\].
pub fn miss_distance(r_rel: [f64; 3], v_rel: [f64; 3]) -> f64 {
    let tca = time_to_closest_approach(r_rel, v_rel);
    let r_tca = add3(r_rel, scale3(v_rel, tca));
    norm3(r_tca)
}

// ---------------------------------------------------------------------------
// Orbit propagator (simple 2-body + J2)
// ---------------------------------------------------------------------------

/// Propagate a satellite state one time step using RK4 with J2 perturbation.
///
/// # Arguments
/// * `pos` — Current position \[m\].
/// * `vel` — Current velocity \[m/s\].
/// * `dt` — Time step \[s\].
/// * `mu` — Gravitational parameter \[m³/s²\].
/// * `j2` — J2 coefficient (set to 0 to ignore).
/// * `r_eq` — Equatorial radius \[m\] (used for J2).
///
/// # Returns
/// `(new_position, new_velocity)`.
pub fn propagate_rk4(
    pos: [f64; 3],
    vel: [f64; 3],
    dt: f64,
    mu: f64,
    j2: f64,
    r_eq: f64,
) -> ([f64; 3], [f64; 3]) {
    let state = [pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]];

    let k1 = orbit_derivative(state, mu, j2, r_eq);
    let k2 = orbit_derivative(add_state(state, scale_state(k1, dt / 2.0)), mu, j2, r_eq);
    let k3 = orbit_derivative(add_state(state, scale_state(k2, dt / 2.0)), mu, j2, r_eq);
    let k4 = orbit_derivative(add_state(state, scale_state(k3, dt)), mu, j2, r_eq);

    let new_state = add_state(
        state,
        scale_state(
            add_state(
                add_state(k1, scale_state(k2, 2.0)),
                add_state(scale_state(k3, 2.0), k4),
            ),
            dt / 6.0,
        ),
    );

    (
        [new_state[0], new_state[1], new_state[2]],
        [new_state[3], new_state[4], new_state[5]],
    )
}

/// Compute the derivative of the orbital state vector \[pos, vel\].
fn orbit_derivative(state: [f64; 6], mu: f64, j2: f64, r_eq: f64) -> [f64; 6] {
    let pos = [state[0], state[1], state[2]];
    let vel = [state[3], state[4], state[5]];
    let r = norm3(pos);
    let r3 = r * r * r;

    // Two-body acceleration
    let a_2body = scale3(pos, -mu / r3);

    // J2 perturbation
    let a_j2 = if j2.abs() > 1e-15 {
        j2_acceleration(pos, mu, j2, r_eq)
    } else {
        [0.0; 3]
    };

    let a_total = add3(a_2body, a_j2);

    [vel[0], vel[1], vel[2], a_total[0], a_total[1], a_total[2]]
}

/// Add two 6-element state vectors.
fn add_state(a: [f64; 6], b: [f64; 6]) -> [f64; 6] {
    [
        a[0] + b[0],
        a[1] + b[1],
        a[2] + b[2],
        a[3] + b[3],
        a[4] + b[4],
        a[5] + b[5],
    ]
}

/// Scale a 6-element state vector by a scalar.
fn scale_state(a: [f64; 6], s: f64) -> [f64; 6] {
    [a[0] * s, a[1] * s, a[2] * s, a[3] * s, a[4] * s, a[5] * s]
}

// ---------------------------------------------------------------------------
// Orbit ground track
// ---------------------------------------------------------------------------

/// Compute latitude and longitude for a ground track point.
///
/// # Arguments
/// * `pos` — ECI position vector \[m\].
/// * `gmst` — Greenwich Mean Sidereal Time \[rad\].
///
/// # Returns
/// `(latitude [rad], longitude [rad])`.
pub fn ground_track_point(pos: [f64; 3], gmst: f64) -> (f64, f64) {
    let r = norm3(pos);
    let lat = (pos[2] / r).asin();
    let lon_eci = pos[1].atan2(pos[0]);
    let lon = (lon_eci - gmst).rem_euclid(2.0 * PI) - PI;
    (lat, lon)
}

// ---------------------------------------------------------------------------
// Sun-synchronous orbit
// ---------------------------------------------------------------------------

/// Compute the inclination required for a Sun-synchronous orbit.
///
/// For a Sun-synchronous orbit, the RAAN precession rate must equal
/// the mean motion of the Earth around the Sun (~0.9856 deg/day).
///
/// # Arguments
/// * `semi_major_axis` — Semi-major axis \[m\].
/// * `eccentricity` — Eccentricity.
/// * `mu` — Gravitational parameter \[m³/s²\].
/// * `j2` — J2 coefficient.
/// * `r_eq` — Equatorial radius \[m\].
///
/// # Returns
/// Required inclination \[rad\].
pub fn sun_synchronous_inclination(
    semi_major_axis: f64,
    eccentricity: f64,
    mu: f64,
    j2: f64,
    r_eq: f64,
) -> f64 {
    // Required precession rate: ~1.991e-7 rad/s
    let omega_earth_sun = 2.0 * PI / (365.25 * 24.0 * 3600.0);
    let n = (mu / semi_major_axis.powi(3)).sqrt();
    let e2 = 1.0 - eccentricity * eccentricity;
    // dΩ/dt = -3/2 * n * J2 * (R/a)^2 * cos(i) / (1-e²)² = omega_earth_sun
    let cos_i = -omega_earth_sun * (e2 * e2) / (1.5 * n * j2 * (r_eq / semi_major_axis).powi(2));
    cos_i.clamp(-1.0, 1.0).acos()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    // ---- Kepler's equation ----

    #[test]
    fn test_kepler_circular_orbit_e0() {
        // e=0: E = M always
        let m = 1.0;
        let e = solve_kepler_equation(m, 0.0, 1e-12, 100);
        assert!((e - m).abs() < TOL, "e=0: E should equal M, got {:.6}", e);
    }

    #[test]
    fn test_kepler_equation_residual() {
        // M = E - e*sin(E) should hold to tolerance
        let m = 2.0;
        let ecc = 0.5;
        let e = solve_kepler_equation(m, ecc, 1e-12, 100);
        let residual = (e - ecc * e.sin() - m).abs();
        assert!(
            residual < 1e-10,
            "Kepler residual too large: {:.2e}",
            residual
        );
    }

    #[test]
    fn test_kepler_high_eccentricity() {
        let m = 1.5;
        let ecc = 0.9;
        let e = solve_kepler_equation(m, ecc, 1e-12, 100);
        let residual = (e - ecc * e.sin() - m).abs();
        assert!(residual < 1e-10, "High-ecc residual: {:.2e}", residual);
    }

    #[test]
    fn test_mean_to_true_anomaly_circular() {
        // For e≈0 the true anomaly should match mean anomaly closely
        let m = 1.0;
        let nu = mean_to_true_anomaly(m, 1e-8);
        assert!((nu - m).abs() < 1e-5, "nu={:.6} m={:.6}", nu, m);
    }

    #[test]
    fn test_true_to_mean_roundtrip() {
        let nu0 = 1.2;
        let ecc = 0.4;
        let m = true_to_mean_anomaly(nu0, ecc);
        let nu1 = mean_to_true_anomaly(m, ecc);
        assert!(
            (nu1 - nu0).abs() < 1e-9,
            "Roundtrip error: {:.2e}",
            (nu1 - nu0).abs()
        );
    }

    // ---- Orbital period / vis-viva ----

    #[test]
    fn test_orbital_period_leo() {
        // ISS orbit ≈ 6778 km radius → ~5559 s period
        let r = 6.778e6;
        let t = orbital_period(r, MU_EARTH);
        assert!((t - 5559.0).abs() < 50.0, "ISS period off: {:.1}", t);
    }

    #[test]
    fn test_vis_viva_circular() {
        // Circular orbit: v² = μ/r
        let r = 7e6;
        let v = vis_viva_speed(r, r, MU_EARTH);
        let v_circ = circular_speed(r, MU_EARTH);
        assert!((v - v_circ).abs() < 1e-3, "vis-viva circular mismatch");
    }

    #[test]
    fn test_escape_speed_ratio() {
        // v_esc = √2 * v_circ
        let r = 7e6;
        let v_c = circular_speed(r, MU_EARTH);
        let v_e = escape_speed(r, MU_EARTH);
        assert!((v_e / v_c - 2.0_f64.sqrt()).abs() < 1e-10, "ratio mismatch");
    }

    // ---- Kepler elements / state vector ----

    #[test]
    fn test_kepler_to_cartesian_circular_equatorial() {
        // Circular equatorial orbit at r=7e6 m, nu=0 → position along x-axis
        let elem = KeplerElements::new(7e6, 0.0, 0.0, 0.0, 0.0, 0.0);
        let (pos, _vel) = elem.to_state_vector(MU_EARTH);
        assert!(
            (norm3(pos) - 7e6).abs() < 10.0,
            "radius mismatch: {}",
            norm3(pos)
        );
        assert!(pos[1].abs() < 1.0, "y should be ~0: {}", pos[1]);
        assert!(pos[2].abs() < 1.0, "z should be ~0: {}", pos[2]);
    }

    #[test]
    fn test_cartesian_to_kepler_roundtrip() {
        let elem_in = KeplerElements::new(8e6, 0.1, 0.5, 1.0, 0.5, 1.2);
        let (pos, vel) = elem_in.to_state_vector(MU_EARTH);
        let elem_out = cartesian_to_kepler(pos, vel, MU_EARTH);
        assert!(
            (elem_out.semi_major_axis - elem_in.semi_major_axis).abs() < 1.0,
            "SMA: {} vs {}",
            elem_out.semi_major_axis,
            elem_in.semi_major_axis
        );
        assert!(
            (elem_out.eccentricity - elem_in.eccentricity).abs() < 1e-8,
            "ecc: {} vs {}",
            elem_out.eccentricity,
            elem_in.eccentricity
        );
    }

    #[test]
    fn test_periapsis_apoapsis_radii() {
        let elem = KeplerElements::new(8e6, 0.2, 0.0, 0.0, 0.0, 0.0);
        assert!((elem.periapsis_radius() - 6.4e6).abs() < 1.0);
        assert!((elem.apoapsis_radius() - 9.6e6).abs() < 1.0);
    }

    #[test]
    fn test_hyperbolic_orbit_apoapsis_infinity() {
        let elem = KeplerElements::new(8e6, 1.5, 0.0, 0.0, 0.0, 0.0);
        assert!(elem.apoapsis_radius().is_infinite());
        assert!(elem.is_hyperbolic());
    }

    // ---- Hohmann transfer ----

    #[test]
    fn test_hohmann_leo_to_geo() {
        // LEO ≈ 6778 km, GEO ≈ 42164 km
        let r_leo = 6.778e6;
        let r_geo = 4.2164e7;
        let h = hohmann_transfer(r_leo, r_geo, MU_EARTH);
        assert!(h.delta_v1 > 0.0, "dv1 must be positive");
        assert!(h.delta_v2 > 0.0, "dv2 must be positive");
        // Known total ΔV ≈ 3.9 km/s
        assert!(
            h.delta_v_total > 3.5e3 && h.delta_v_total < 4.2e3,
            "Total ΔV = {:.0} m/s",
            h.delta_v_total
        );
    }

    #[test]
    fn test_hohmann_transfer_sma() {
        let r1 = 7e6;
        let r2 = 14e6;
        let h = hohmann_transfer(r1, r2, MU_EARTH);
        assert!((h.transfer_sma - (r1 + r2) / 2.0).abs() < 1.0);
    }

    // ---- Bi-elliptic transfer ----

    #[test]
    fn test_bi_elliptic_large_ratio() {
        // For r2/r1 > ~11.94 bi-elliptic is more efficient
        let r1 = 7e6;
        let r2 = 1.0e8;
        let rb = 1.5e8;
        let be = bi_elliptic_transfer(r1, r2, rb, MU_EARTH);
        assert!(be.delta_v_total > 0.0);
        let ho = hohmann_transfer(r1, r2, MU_EARTH);
        // Bi-elliptic should be comparable (may be less efficient for this rb)
        assert!(be.delta_v_total < ho.delta_v_total * 2.0);
    }

    // ---- Sphere of influence ----

    #[test]
    fn test_earth_soi_from_sun() {
        // Earth SOI from Sun ≈ 924,000 km
        let soi = sphere_of_influence(AU, MU_EARTH, MU_SUN);
        assert!(soi > 9e8 && soi < 1.0e9, "Earth SOI: {:.3e}", soi);
    }

    // ---- J2 perturbation ----

    #[test]
    fn test_j2_acceleration_equatorial() {
        // At the equator (z=0), only radial J2 component
        let r = R_EARTH_EQ + 400e3;
        let pos = [r, 0.0, 0.0];
        let a = j2_acceleration(pos, MU_EARTH, J2_EARTH, R_EARTH_EQ);
        assert!(
            a[2].abs() < 1e-10,
            "z-component should be 0 at equator: {}",
            a[2]
        );
        assert!(a[0].abs() > 0.0, "radial component should be nonzero");
    }

    #[test]
    fn test_j2_raan_precession_iss() {
        // ISS inclination ≈ 51.6°, a ≈ 6778 km
        let inc = 51.6_f64.to_radians();
        let a = 6.778e6;
        let rate = j2_raan_precession_rate(a, 0.0, inc, MU_EARTH, J2_EARTH, R_EARTH_EQ);
        // Expected ≈ -0.176 deg/day → -3.55e-7 rad/s
        assert!(rate < 0.0, "RAAN precession should be negative");
        assert!(
            rate.abs() > 1e-7 && rate.abs() < 2e-6,
            "rate = {:.3e}",
            rate
        );
    }

    // ---- Atmospheric drag ----

    #[test]
    fn test_atmospheric_density_decreases_with_altitude() {
        let rho1 = atmospheric_density(200e3, RHO_SEA_LEVEL, 0.0, SCALE_HEIGHT_ATM);
        let rho2 = atmospheric_density(400e3, RHO_SEA_LEVEL, 0.0, SCALE_HEIGHT_ATM);
        assert!(rho2 < rho1, "density should decrease with altitude");
    }

    #[test]
    fn test_drag_acceleration_direction() {
        // Drag should oppose velocity
        let vel = [7000.0, 0.0, 0.0];
        let a = atmospheric_drag_acceleration(vel, 400e3, 2.2, 0.01, 1e-11, 400e3, 8500.0);
        assert!(a[0] < 0.0, "drag x-component should be negative: {}", a[0]);
        assert!(a[1].abs() < 1e-20, "y-component should be zero: {}", a[1]);
    }

    // ---- CR3BP ----

    #[test]
    fn test_jacobi_integral_conservation() {
        // Integrate a few steps and verify Jacobi integral is approximately conserved
        let mu = 0.01215; // Earth-Moon mass ratio
        let state0 = [0.8, 0.0, 0.0, 0.0, 0.5, 0.0];
        let c0 = jacobi_integral(state0, mu);

        let dt = 0.001;
        let mut state = state0;
        for _ in 0..100 {
            let ds = cr3bp_equations(state, mu);
            // Euler step (small dt so ΔC should be small)
            for k in 0..6 {
                state[k] += ds[k] * dt;
            }
        }
        let c1 = jacobi_integral(state, mu);
        assert!(
            (c1 - c0).abs() < 0.01,
            "Jacobi not conserved: Δ={:.4}",
            (c1 - c0).abs()
        );
    }

    // ---- Lagrange points ----

    #[test]
    fn test_lagrange_l4_l5_equilateral() {
        let mu = 0.01215;
        let lp = lagrange_points(mu);
        // L4 and L5 should be at y = ±√3/2
        let expected_y = 3.0_f64.sqrt() / 2.0;
        assert!((lp[3][1] - expected_y).abs() < 1e-6, "L4 y: {}", lp[3][1]);
        assert!((lp[4][1] + expected_y).abs() < 1e-6, "L5 y: {}", lp[4][1]);
    }

    #[test]
    fn test_lagrange_l1_between_primaries() {
        let mu = 0.01215;
        let lp = lagrange_points(mu);
        let l1_x = lp[0][0];
        // L1 should be between x = -mu and x = 1-mu
        assert!(
            l1_x > -mu && l1_x < 1.0 - mu,
            "L1 = {} not between primaries",
            l1_x
        );
    }

    #[test]
    fn test_lagrange_l2_beyond_secondary() {
        let mu = 0.01215;
        let lp = lagrange_points(mu);
        let l2_x = lp[1][0];
        // L2 is beyond the secondary (x > 1-mu)
        assert!(l2_x > 1.0 - mu, "L2 = {} should be beyond secondary", l2_x);
    }

    // ---- Collision probability ----

    #[test]
    fn test_collision_probability_zero_radius() {
        let p = collision_probability(0.0, 1000.0, 100.0, 10.0);
        assert!(p < 1e-20, "zero radius → zero probability: {}", p);
    }

    #[test]
    fn test_collision_probability_bounded() {
        let p = collision_probability(10.0, 1000.0, 1.0, 100.0);
        assert!((0.0..=1.0).contains(&p), "probability out of bounds: {}", p);
    }

    #[test]
    fn test_miss_distance_head_on() {
        // Head-on collision: r_rel = (10,0,0), v_rel = (-1,0,0) → TCA at t=10, miss=0
        let r = [10.0, 0.0, 0.0];
        let v = [-1.0, 0.0, 0.0];
        let d = miss_distance(r, v);
        assert!(d < 1e-10, "head-on miss distance: {}", d);
    }

    #[test]
    fn test_miss_distance_parallel() {
        // Parallel trajectories: r_rel = (0,5,0), v_rel = (1,0,0) → miss = 5
        let r = [0.0, 5.0, 0.0];
        let v = [1.0, 0.0, 0.0];
        let d = miss_distance(r, v);
        assert!((d - 5.0).abs() < 1e-10, "parallel miss distance: {}", d);
    }

    // ---- Propagator ----

    #[test]
    fn test_rk4_conserves_energy() {
        // Propagate a circular orbit for 1000 steps; specific energy should be ~constant
        let r0 = 7e6;
        let mut pos = [r0, 0.0, 0.0];
        let v_circ = circular_speed(r0, MU_EARTH);
        let mut vel = [0.0, v_circ, 0.0];
        let dt = 10.0;
        let energy0 = norm3(vel).powi(2) / 2.0 - MU_EARTH / norm3(pos);

        for _ in 0..1000 {
            let (p, v) = propagate_rk4(pos, vel, dt, MU_EARTH, 0.0, R_EARTH_EQ);
            pos = p;
            vel = v;
        }
        let energy1 = norm3(vel).powi(2) / 2.0 - MU_EARTH / norm3(pos);
        let rel_err = ((energy1 - energy0) / energy0).abs();
        assert!(rel_err < 1e-4, "Energy drift: {:.2e}", rel_err);
    }

    #[test]
    fn test_rk4_j2_changes_trajectory() {
        let r0 = 7e6;
        let mut pos_no_j2 = [r0, 0.0, 0.0];
        let v_circ = circular_speed(r0, MU_EARTH);
        let mut vel_no_j2 = [0.0, v_circ, 0.0];
        let mut pos_j2 = pos_no_j2;
        let mut vel_j2 = vel_no_j2;
        let dt = 10.0;

        for _ in 0..100 {
            let (p, v) = propagate_rk4(pos_no_j2, vel_no_j2, dt, MU_EARTH, 0.0, R_EARTH_EQ);
            pos_no_j2 = p;
            vel_no_j2 = v;
            let (p2, v2) = propagate_rk4(pos_j2, vel_j2, dt, MU_EARTH, J2_EARTH, R_EARTH_EQ);
            pos_j2 = p2;
            vel_j2 = v2;
        }
        let diff = norm3(sub3(pos_j2, pos_no_j2));
        assert!(
            diff > 1e-3,
            "J2 should produce trajectory difference: {}",
            diff
        );
    }

    // ---- Sun-synchronous ----

    #[test]
    fn test_sun_synchronous_inclination_range() {
        // SSO inclinations are typically 96–100° for LEO
        let a = 7e6;
        let inc = sun_synchronous_inclination(a, 0.0, MU_EARTH, J2_EARTH, R_EARTH_EQ);
        let deg = inc.to_degrees();
        assert!(deg > 90.0 && deg < 110.0, "SSO inclination: {:.2}°", deg);
    }

    // ---- Ground track ----

    #[test]
    fn test_ground_track_latitude_range() {
        let pos = [7e6 * 0.7, 7e6 * 0.7, 7e6 * 0.14];
        let (lat, lon) = ground_track_point(pos, 0.0);
        assert!(
            (-PI / 2.0..=PI / 2.0).contains(&lat),
            "lat out of range: {}",
            lat
        );
        assert!((-PI..=PI).contains(&lon), "lon out of range: {}", lon);
    }

    // ---- Lambert's problem (basic sanity) ----

    #[test]
    fn test_lambert_returns_solution() {
        // Simple case: same orbit, half-period transfer (Hohmann-like)
        let r1 = [7e6, 0.0, 0.0];
        let r2 = [-7e6, 0.0, 0.0];
        let tof = orbital_period(7e6, MU_EARTH) / 2.0;
        let sol = lambert_problem(r1, r2, tof, MU_EARTH, true);
        // May or may not converge for this degenerate case; just test it doesn't panic
        let _ = sol;
    }

    // ---- Dot / cross / norm helpers ----

    #[test]
    fn test_vec3_helpers() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = cross3(a, b);
        assert!((c[2] - 1.0).abs() < 1e-15);
        assert!((dot3(a, b)).abs() < 1e-15);
        assert!((norm3(a) - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_normalize3_unit_vector() {
        let v = [3.0, 4.0, 0.0];
        let n = normalize3(v);
        assert!(
            (norm3(n) - 1.0).abs() < 1e-15,
            "normalized norm: {}",
            norm3(n)
        );
    }
}
