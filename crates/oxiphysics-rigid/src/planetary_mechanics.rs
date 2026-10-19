// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Planetary mechanics: orbital mechanics, N-body gravity, Kepler's laws,
//! orbital maneuvers, perturbations, and TLE parsing.
//!
//! All structures use plain `f64` arrays — no nalgebra dependency.
//!
//! # Overview
//!
//! * [`OrbitalElements`] — Classical Keplerian elements, state vector conversion.
//! * [`KeplerOrbit`] — Kepler's equation solver (Newton-Raphson).
//! * [`NBodyGravity`] — Direct N-body or Barnes-Hut octree forces (monopole +
//!   quadrupole cells) with a kick-drift-kick leapfrog integrator.
//! * [`OrbitalManeuver`] — Hohmann, bielliptic, plane-change maneuvers.
//! * [`J2Perturbation`] — J2 oblateness nodal precession.
//! * [`AtmosphericDrag`] — tabulated Vallado exponential atmosphere, drag
//!   deceleration and circular-orbit decay rate.
//! * [`LagrangePoints`] — L1–L5 for circular restricted 3-body problem.
//! * [`TleParser`] — TLE parsing and SGP4-like propagation.
//! * [`HillClohessyWiltshire`] — exact relative-motion propagation and
//!   two-impulse rendezvous targeting ([`hcw_two_impulse`]).

use std::f64::consts::PI;

mod hcw;
mod nbody;
#[cfg(test)]
mod nbody_tests;

pub use hcw::hcw_two_impulse;
pub use nbody::{GravBody, NBodyGravity, OctreeNode};

// ─── Physical constants ────────────────────────────────────────────────────────

/// Earth gravitational parameter μ = GM \[m³/s²\].
pub const MU_EARTH: f64 = 3.986_004_418e14;

/// Earth radius \[m\].
pub const R_EARTH: f64 = 6.371_000e6;

/// Earth's second zonal harmonic (J2 oblateness coefficient).
pub const J2_COEFF: f64 = 1.082_626_68e-3;

/// Universal gravitational constant G \[m³/(kg·s²)\].
pub const G_CONST: f64 = 6.674_30e-11;

/// Solar gravitational parameter μ_Sun \[m³/s²\].
pub const MU_SUN: f64 = 1.327_124_4e20;

/// Mean Earth–Moon distance \[m\].
pub const R_MOON: f64 = 3.844e8;

/// Moon gravitational parameter \[m³/s²\].
pub const MU_MOON: f64 = 4.9048695e12;

/// Representative thermospheric reference density at 200 km \[kg/m³\].
///
/// Legacy single-layer reference value kept for API compatibility;
/// [`AtmosphericDrag::density`] uses the tabulated Vallado model instead.
pub const RHO_REF: f64 = 2.5e-10;

/// Legacy single-layer atmospheric scale height \[m\].
///
/// Not used by [`AtmosphericDrag::density`], which carries a scale height per
/// altitude band.
pub const SCALE_HEIGHT: f64 = 8_500.0;

/// Legacy reference altitude of the single-layer drag model \[m\].
///
/// Not used by [`AtmosphericDrag::density`].
pub const H_REF: f64 = 200_000.0;

/// Barnes-Hut opening angle theta (default).
pub const BH_THETA: f64 = 0.5;

// ─── Vec3 helpers ─────────────────────────────────────────────────────────────

#[inline]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

#[inline]
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn vec3_normalize(a: [f64; 3]) -> [f64; 3] {
    let n = vec3_norm(a);
    if n < 1e-30 {
        [0.0; 3]
    } else {
        vec3_scale(a, 1.0 / n)
    }
}

// ─── OrbitalElements ──────────────────────────────────────────────────────────

/// Classical Keplerian orbital elements describing a two-body orbit.
///
/// These six elements uniquely define the size, shape, and orientation of an
/// orbit together with the spacecraft position at a given epoch.
#[derive(Debug, Clone, PartialEq)]
pub struct OrbitalElements {
    /// Semi-major axis \[m\].
    pub semi_major_axis: f64,
    /// Orbital eccentricity (dimensionless, 0 ≤ e < 1 for elliptic).
    pub eccentricity: f64,
    /// Inclination \[rad\].
    pub inclination: f64,
    /// Right ascension of the ascending node (RAAN) \[rad\].
    pub raan: f64,
    /// Argument of perigee \[rad\].
    pub arg_perigee: f64,
    /// True anomaly at epoch \[rad\].
    pub true_anomaly: f64,
}

impl OrbitalElements {
    /// Create a new set of orbital elements.
    pub fn new(
        semi_major_axis: f64,
        eccentricity: f64,
        inclination: f64,
        raan: f64,
        arg_perigee: f64,
        true_anomaly: f64,
    ) -> Self {
        Self {
            semi_major_axis,
            eccentricity,
            inclination,
            raan,
            arg_perigee,
            true_anomaly,
        }
    }

    /// Orbital period \[s\] using Kepler's third law: T = 2π√(a³/μ).
    pub fn period(&self, mu: f64) -> f64 {
        2.0 * PI * (self.semi_major_axis.powi(3) / mu).sqrt()
    }

    /// Periapsis radius \[m\].
    pub fn periapsis(&self) -> f64 {
        self.semi_major_axis * (1.0 - self.eccentricity)
    }

    /// Apoapsis radius \[m\].
    pub fn apoapsis(&self) -> f64 {
        self.semi_major_axis * (1.0 + self.eccentricity)
    }

    /// Semi-latus rectum \[m\]: p = a(1 - e²).
    pub fn semi_latus_rectum(&self) -> f64 {
        self.semi_major_axis * (1.0 - self.eccentricity * self.eccentricity)
    }

    /// Orbital radius at the current true anomaly \[m\].
    pub fn radius(&self) -> f64 {
        let p = self.semi_latus_rectum();
        p / (1.0 + self.eccentricity * self.true_anomaly.cos())
    }

    /// Convert orbital elements to an ECI state vector (position, velocity).
    ///
    /// Returns `([x,y,z], [vx,vy,vz])` in metres and m/s.
    pub fn to_state_vector(&self, mu: f64) -> ([f64; 3], [f64; 3]) {
        let p = self.semi_latus_rectum();
        let e = self.eccentricity;
        let nu = self.true_anomaly;
        let r = p / (1.0 + e * nu.cos());

        // Position in perifocal frame
        let r_pqw = [r * nu.cos(), r * nu.sin(), 0.0];
        let v_fac = (mu / p).sqrt();
        let v_pqw = [v_fac * (-nu.sin()), v_fac * (e + nu.cos()), 0.0];

        // Rotation matrices
        let (sin_raan, cos_raan) = self.raan.sin_cos();
        let (sin_inc, cos_inc) = self.inclination.sin_cos();
        let (sin_w, cos_w) = self.arg_perigee.sin_cos();

        // DCM from perifocal to ECI
        let r11 = cos_raan * cos_w - sin_raan * sin_w * cos_inc;
        let r12 = -cos_raan * sin_w - sin_raan * cos_w * cos_inc;
        let r21 = sin_raan * cos_w + cos_raan * sin_w * cos_inc;
        let r22 = -sin_raan * sin_w + cos_raan * cos_w * cos_inc;
        let r31 = sin_w * sin_inc;
        let r32 = cos_w * sin_inc;

        let rx = [r11, r21, r31];
        let ry = [r12, r22, r32];
        let rz = [0.0_f64, 0.0, 0.0]; // z column (not used for r_pqw[2]=0)

        let pos = [
            rx[0] * r_pqw[0] + ry[0] * r_pqw[1] + rz[0] * r_pqw[2],
            rx[1] * r_pqw[0] + ry[1] * r_pqw[1] + rz[1] * r_pqw[2],
            rx[2] * r_pqw[0] + ry[2] * r_pqw[1] + rz[2] * r_pqw[2],
        ];
        let vel = [
            rx[0] * v_pqw[0] + ry[0] * v_pqw[1],
            rx[1] * v_pqw[0] + ry[1] * v_pqw[1],
            rx[2] * v_pqw[0] + ry[2] * v_pqw[1],
        ];
        (pos, vel)
    }

    /// Convert an ECI state vector back to orbital elements.
    ///
    /// * `pos` — position vector \[m\]
    /// * `vel` — velocity vector \[m/s\]
    /// * `mu`  — gravitational parameter \[m³/s²\]
    pub fn from_state_vector(pos: [f64; 3], vel: [f64; 3], mu: f64) -> Self {
        let r = vec3_norm(pos);
        let v2 = vec3_dot(vel, vel);
        let rdotv = vec3_dot(pos, vel);

        // Angular momentum vector h = r × v
        let h_vec = vec3_cross(pos, vel);
        let h = vec3_norm(h_vec);

        // Node vector N = k̂ × h
        let k_hat = [0.0_f64, 0.0, 1.0];
        let node_vec = vec3_cross(k_hat, h_vec);
        let node_mag = vec3_norm(node_vec);

        // Eccentricity vector e = (v × h)/μ − r̂
        let vh = vec3_cross(vel, h_vec);
        let r_hat = vec3_normalize(pos);
        let ecc_vec = [
            vh[0] / mu - r_hat[0],
            vh[1] / mu - r_hat[1],
            vh[2] / mu - r_hat[2],
        ];
        let e = vec3_norm(ecc_vec);

        // Semi-major axis
        let energy = v2 / 2.0 - mu / r;
        let a = if energy.abs() < 1e-12 {
            1e12 // parabolic
        } else {
            -mu / (2.0 * energy)
        };

        // Inclination
        let inc = (h_vec[2] / h).clamp(-1.0, 1.0).acos();

        // RAAN
        let raan = if node_mag < 1e-10 {
            0.0
        } else if node_vec[1] >= 0.0 {
            (node_vec[0] / node_mag).clamp(-1.0, 1.0).acos()
        } else {
            2.0 * PI - (node_vec[0] / node_mag).clamp(-1.0, 1.0).acos()
        };

        // Argument of perigee
        let arg_p = if node_mag < 1e-10 || e < 1e-10 {
            0.0
        } else {
            let cos_w = vec3_dot(node_vec, ecc_vec) / (node_mag * e);
            let w = cos_w.clamp(-1.0, 1.0).acos();
            if ecc_vec[2] >= 0.0 { w } else { 2.0 * PI - w }
        };

        // True anomaly
        let nu = if e < 1e-10 {
            let cos_nu = vec3_dot(r_hat, vec3_normalize(h_vec));
            cos_nu.clamp(-1.0, 1.0).acos()
        } else {
            let cos_nu = vec3_dot(ecc_vec, pos) / (e * r);
            let nu0 = cos_nu.clamp(-1.0, 1.0).acos();
            if rdotv >= 0.0 { nu0 } else { 2.0 * PI - nu0 }
        };

        OrbitalElements::new(a, e, inc, raan, arg_p, nu)
    }
}

// ─── KeplerOrbit ──────────────────────────────────────────────────────────────

/// Kepler's equation solver and related orbit propagation utilities.
///
/// Solves M = E − e·sin(E) for the eccentric anomaly E using Newton-Raphson
/// iteration, then converts to true anomaly ν.
pub struct KeplerOrbit {
    /// Gravitational parameter \[m³/s²\].
    pub mu: f64,
    /// Maximum iterations for Newton-Raphson solver.
    pub max_iter: usize,
    /// Convergence tolerance \[rad\].
    pub tol: f64,
}

impl KeplerOrbit {
    /// Create a new Kepler orbit solver.
    pub fn new(mu: f64) -> Self {
        Self {
            mu,
            max_iter: 50,
            tol: 1e-12,
        }
    }

    /// Solve Kepler's equation M = E − e·sin(E) for eccentric anomaly E \[rad\].
    ///
    /// Uses Newton-Raphson iteration starting from E₀ = M.
    pub fn solve_kepler(&self, mean_anomaly: f64, eccentricity: f64) -> f64 {
        let m = mean_anomaly.rem_euclid(2.0 * PI);
        let mut e = if eccentricity > 0.8 { PI } else { m };
        for _ in 0..self.max_iter {
            let f = e - eccentricity * e.sin() - m;
            let fp = 1.0 - eccentricity * e.cos();
            let delta = f / fp;
            e -= delta;
            if delta.abs() < self.tol {
                break;
            }
        }
        e
    }

    /// Convert eccentric anomaly E to true anomaly ν \[rad\].
    pub fn eccentric_to_true(&self, e_anom: f64, eccentricity: f64) -> f64 {
        let half = (e_anom / 2.0).tan();
        let ratio = ((1.0 + eccentricity) / (1.0 - eccentricity)).sqrt();
        2.0 * (ratio * half).atan()
    }

    /// Convert true anomaly ν to eccentric anomaly E \[rad\].
    pub fn true_to_eccentric(&self, true_anom: f64, eccentricity: f64) -> f64 {
        let half = (true_anom / 2.0).tan();
        let ratio = ((1.0 - eccentricity) / (1.0 + eccentricity)).sqrt();
        2.0 * (ratio * half).atan()
    }

    /// Convert eccentric anomaly E to mean anomaly M \[rad\].
    pub fn eccentric_to_mean(&self, e_anom: f64, eccentricity: f64) -> f64 {
        e_anom - eccentricity * e_anom.sin()
    }

    /// Propagate orbital elements by time `dt` \[s\].
    ///
    /// Updates the true anomaly field using Kepler's equation.
    pub fn propagate(&self, elems: &OrbitalElements, dt: f64) -> OrbitalElements {
        let n = (self.mu / elems.semi_major_axis.powi(3)).sqrt(); // mean motion
        let e0 = self.true_to_eccentric(elems.true_anomaly, elems.eccentricity);
        let m0 = self.eccentric_to_mean(e0, elems.eccentricity);
        let m1 = (m0 + n * dt).rem_euclid(2.0 * PI);
        let e1 = self.solve_kepler(m1, elems.eccentricity);
        let nu1 = self.eccentric_to_true(e1, elems.eccentricity);
        OrbitalElements::new(
            elems.semi_major_axis,
            elems.eccentricity,
            elems.inclination,
            elems.raan,
            elems.arg_perigee,
            nu1,
        )
    }

    /// Compute orbital speed at a given radius r \[m\] using the vis-viva equation.
    pub fn vis_viva_speed(&self, r: f64, a: f64) -> f64 {
        (self.mu * (2.0 / r - 1.0 / a)).sqrt()
    }
}

// ─── NBodyGravity ─────────────────────────────────────────────────────────────
// `GravBody`, `NBodyGravity` and `OctreeNode` live in the `nbody` submodule.

// ─── OrbitalManeuver ──────────────────────────────────────────────────────────

/// Orbital maneuver calculations: Hohmann transfers, bielliptic transfers,
/// plane changes, and delta-v budget.
pub struct OrbitalManeuver {
    /// Gravitational parameter \[m³/s²\].
    pub mu: f64,
}

impl OrbitalManeuver {
    /// Create a new maneuver calculator.
    pub fn new(mu: f64) -> Self {
        Self { mu }
    }

    /// Compute delta-v for a Hohmann transfer between circular orbits.
    ///
    /// Returns `(dv1, dv2, total_dv)` \[m/s\].
    ///
    /// * `r1` — initial circular orbit radius \[m\]
    /// * `r2` — final circular orbit radius \[m\]
    pub fn hohmann_transfer(&self, r1: f64, r2: f64) -> (f64, f64, f64) {
        let a_transfer = (r1 + r2) / 2.0;
        let v1 = (self.mu / r1).sqrt();
        let v2 = (self.mu / r2).sqrt();
        let v_trans1 = (self.mu * (2.0 / r1 - 1.0 / a_transfer)).sqrt();
        let v_trans2 = (self.mu * (2.0 / r2 - 1.0 / a_transfer)).sqrt();
        let dv1 = (v_trans1 - v1).abs();
        let dv2 = (v2 - v_trans2).abs();
        (dv1, dv2, dv1 + dv2)
    }

    /// Compute delta-v for a bielliptic transfer.
    ///
    /// Returns `(dv1, dv2, dv3, total_dv)` \[m/s\].
    ///
    /// * `r1` — initial orbit radius \[m\]
    /// * `r2` — final orbit radius \[m\]
    /// * `rb` — intermediate apoapsis radius \[m\]
    pub fn bielliptic_transfer(&self, r1: f64, r2: f64, rb: f64) -> (f64, f64, f64, f64) {
        let v1 = (self.mu / r1).sqrt();
        let v2 = (self.mu / r2).sqrt();
        let a1 = (r1 + rb) / 2.0;
        let a2 = (rb + r2) / 2.0;
        let v_t1_peri = (self.mu * (2.0 / r1 - 1.0 / a1)).sqrt();
        let v_t1_apo = (self.mu * (2.0 / rb - 1.0 / a1)).sqrt();
        let v_t2_apo = (self.mu * (2.0 / rb - 1.0 / a2)).sqrt();
        let v_t2_peri = (self.mu * (2.0 / r2 - 1.0 / a2)).sqrt();
        let dv1 = (v_t1_peri - v1).abs();
        let dv2 = (v_t2_apo - v_t1_apo).abs();
        let dv3 = (v2 - v_t2_peri).abs();
        (dv1, dv2, dv3, dv1 + dv2 + dv3)
    }

    /// Compute delta-v for an inclination change maneuver at circular orbit.
    ///
    /// * `r`           — circular orbit radius \[m\]
    /// * `delta_inc`   — inclination change \[rad\]
    pub fn plane_change(&self, r: f64, delta_inc: f64) -> f64 {
        let v = (self.mu / r).sqrt();
        2.0 * v * (delta_inc / 2.0).sin()
    }

    /// Combined plane change with Hohmann transfer (optimal split).
    ///
    /// Returns total delta-v \[m/s\].
    pub fn combined_plane_change_hohmann(&self, r1: f64, r2: f64, delta_inc: f64) -> f64 {
        // Perform plane change at apoapsis (cheaper when v is lower)
        let a_transfer = (r1 + r2) / 2.0;
        let v1 = (self.mu / r1).sqrt();
        let v_trans1 = (self.mu * (2.0 / r1 - 1.0 / a_transfer)).sqrt();
        let v_trans2 = (self.mu * (2.0 / r2 - 1.0 / a_transfer)).sqrt();
        let v2 = (self.mu / r2).sqrt();
        let dv1 = (v_trans1 - v1).abs();
        // Combined plane change + circularization at apoapsis
        let dv2 =
            (v2 * v2 + v_trans2 * v2_sq(v_trans2) - 2.0 * v2 * v_trans2 * delta_inc.cos()).sqrt();
        dv1 + dv2
    }

    /// Tsiolkovsky rocket equation: required propellant mass fraction.
    ///
    /// * `dv`   — delta-v \[m/s\]
    /// * `isp`  — specific impulse \[s\]
    ///
    /// Returns mass ratio m0/mf (initial/final mass).
    pub fn tsiolkovsky_mass_ratio(&self, dv: f64, isp: f64) -> f64 {
        let v_ex = isp * 9.80665; // exhaust velocity [m/s]
        (dv / v_ex).exp()
    }

    /// Time of flight for Hohmann transfer \[s\].
    pub fn hohmann_tof(&self, r1: f64, r2: f64) -> f64 {
        let a_transfer = (r1 + r2) / 2.0;
        PI * (a_transfer.powi(3) / self.mu).sqrt()
    }
}

#[inline]
fn v2_sq(v: f64) -> f64 {
    v
}

// ─── J2Perturbation ───────────────────────────────────────────────────────────

/// J2 oblateness perturbation model for Earth orbit.
///
/// Computes nodal precession (RAAN drift) and argument-of-perigee drift
/// due to the Earth's equatorial bulge.
pub struct J2Perturbation {
    /// Gravitational parameter \[m³/s²\].
    pub mu: f64,
    /// Planet equatorial radius \[m\].
    pub r_eq: f64,
    /// J2 coefficient (dimensionless).
    pub j2: f64,
}

impl J2Perturbation {
    /// Create a new J2 perturbation model with Earth parameters.
    pub fn new_earth() -> Self {
        Self {
            mu: MU_EARTH,
            r_eq: R_EARTH,
            j2: J2_COEFF,
        }
    }

    /// Create with custom parameters.
    pub fn new(mu: f64, r_eq: f64, j2: f64) -> Self {
        Self { mu, r_eq, j2 }
    }

    /// RAAN precession rate dΩ/dt \[rad/s\].
    ///
    /// Negative for prograde orbits (i < 90°).
    pub fn raan_drift_rate(&self, a: f64, e: f64, i: f64) -> f64 {
        let n = (self.mu / a.powi(3)).sqrt();
        let p = a * (1.0 - e * e);
        -1.5 * n * self.j2 * (self.r_eq / p).powi(2) * i.cos()
    }

    /// Argument-of-perigee drift rate dω/dt \[rad/s\].
    pub fn arg_perigee_drift_rate(&self, a: f64, e: f64, i: f64) -> f64 {
        let n = (self.mu / a.powi(3)).sqrt();
        let p = a * (1.0 - e * e);
        0.75 * n * self.j2 * (self.r_eq / p).powi(2) * (5.0 * i.cos().powi(2) - 1.0)
    }

    /// Mean motion drift rate dM/dt \[rad/s\] (additional due to J2).
    pub fn mean_motion_drift(&self, a: f64, e: f64, i: f64) -> f64 {
        let n = (self.mu / a.powi(3)).sqrt();
        let p = a * (1.0 - e * e);
        0.75 * n
            * self.j2
            * (self.r_eq / p).powi(2)
            * (1.0 - e * e).sqrt()
            * (3.0 * i.cos().powi(2) - 1.0)
    }

    /// Propagate orbital elements by `dt` \[s\] including J2 secular drift.
    pub fn propagate(&self, elems: &OrbitalElements, dt: f64) -> OrbitalElements {
        let a = elems.semi_major_axis;
        let e = elems.eccentricity;
        let i = elems.inclination;
        let n = (self.mu / a.powi(3)).sqrt();

        let d_raan = self.raan_drift_rate(a, e, i) * dt;
        let d_w = self.arg_perigee_drift_rate(a, e, i) * dt;
        let d_m = (n + self.mean_motion_drift(a, e, i)) * dt;

        let e0 = KeplerOrbit::new(self.mu).true_to_eccentric(elems.true_anomaly, e);
        let m0 = e0 - e * e0.sin();
        let m1 = (m0 + d_m).rem_euclid(2.0 * PI);
        let e1 = KeplerOrbit::new(self.mu).solve_kepler(m1, e);
        let nu1 = KeplerOrbit::new(self.mu).eccentric_to_true(e1, e);

        OrbitalElements::new(
            a,
            e,
            i,
            (elems.raan + d_raan).rem_euclid(2.0 * PI),
            (elems.arg_perigee + d_w).rem_euclid(2.0 * PI),
            nu1,
        )
    }

    /// Frozen orbit eccentricity (critical inclination condition).
    ///
    /// At the critical inclination ~63.43°, the argument of perigee does
    /// not drift (used for Molniya orbits).
    pub fn critical_inclination() -> f64 {
        // cos²(i) = 1/5 → i = arccos(1/√5)
        (1.0_f64 / 5.0_f64).sqrt().acos()
    }
}

// ─── AtmosphericDrag ──────────────────────────────────────────────────────────

/// Piecewise exponential atmosphere, Vallado (2013), *Fundamentals of
/// Astrodynamics and Applications*, 4th ed., Table 8-4.
///
/// Each row is `(h0 [km], rho0 [kg/m³], H [km])`: inside the band that starts
/// at `h0` the density is `rho0 · exp(-(h - h0) / H)`.  Adjacent bands agree at
/// their shared edge to better than 0.2 %.
const EXPONENTIAL_ATMOSPHERE: [(f64, f64, f64); 28] = [
    (0.0, 1.225, 7.249),
    (25.0, 3.899e-2, 6.349),
    (30.0, 1.774e-2, 6.682),
    (40.0, 3.972e-3, 7.554),
    (50.0, 1.057e-3, 8.382),
    (60.0, 3.206e-4, 7.714),
    (70.0, 8.770e-5, 6.549),
    (80.0, 1.905e-5, 5.799),
    (90.0, 3.396e-6, 5.382),
    (100.0, 5.297e-7, 5.877),
    (110.0, 9.661e-8, 7.263),
    (120.0, 2.438e-8, 9.473),
    (130.0, 8.484e-9, 12.636),
    (140.0, 3.845e-9, 16.149),
    (150.0, 2.070e-9, 22.523),
    (180.0, 5.464e-10, 29.740),
    (200.0, 2.789e-10, 37.105),
    (250.0, 7.248e-11, 45.546),
    (300.0, 2.418e-11, 53.628),
    (350.0, 9.518e-12, 53.298),
    (400.0, 3.725e-12, 58.515),
    (450.0, 1.585e-12, 60.828),
    (500.0, 6.967e-13, 63.822),
    (600.0, 1.454e-13, 71.835),
    (700.0, 3.614e-14, 88.667),
    (800.0, 1.170e-14, 124.64),
    (900.0, 5.245e-15, 181.05),
    (1000.0, 3.019e-15, 268.00),
];

/// Atmospheric drag deceleration model.
///
/// Uses the tabulated piecewise-exponential atmosphere of Vallado (2013),
/// Table 8-4 (surface to 1000 km, extrapolated above with the topmost band).
/// It is a mean-density fit and does not model solar-flux, geomagnetic or
/// diurnal variations.
pub struct AtmosphericDrag {
    /// Drag coefficient (dimensionless), typically 2.2.
    pub cd: f64,
    /// Cross-sectional area \[m²\].
    pub area: f64,
    /// Spacecraft mass \[kg\].
    pub mass: f64,
    /// Ballistic coefficient B = m/(Cd·A) \[kg/m²\].
    pub ballistic_coeff: f64,
}

impl AtmosphericDrag {
    /// Create an atmospheric drag model.
    pub fn new(cd: f64, area: f64, mass: f64) -> Self {
        Self {
            cd,
            area,
            mass,
            ballistic_coeff: mass / (cd * area),
        }
    }

    /// Atmospheric density at altitude h \[m\] above Earth's surface \[kg/m³\].
    ///
    /// Evaluates the Vallado (2013) Table 8-4 exponential model: inside the
    /// band `[h0, h_next)` the density is `rho0 · exp(-(h - h0) / H)`.
    /// Negative altitudes are clamped to the surface, altitudes above 1000 km
    /// extend the topmost band (`H` = 268 km) and a NaN altitude yields NaN.
    pub fn density(&self, altitude: f64) -> f64 {
        if altitude.is_nan() {
            return f64::NAN;
        }
        let h_km = altitude.max(0.0) / 1_000.0;
        let (h0_km, rho0, scale_km) = EXPONENTIAL_ATMOSPHERE
            .iter()
            .rev()
            .find(|band| band.0 <= h_km)
            .copied()
            .unwrap_or(EXPONENTIAL_ATMOSPHERE[0]);
        rho0 * (-(h_km - h0_km) / scale_km).exp()
    }

    /// Drag deceleration magnitude at altitude h with speed v \[m/s\].
    pub fn drag_decel(&self, altitude: f64, speed: f64) -> f64 {
        let rho = self.density(altitude);
        0.5 * rho * self.cd * self.area * speed * speed / self.mass
    }

    /// Drag force vector on a body at altitude h moving with velocity `vel` \[m/s\].
    ///
    /// Returns force vector \[N\].
    pub fn drag_force(&self, altitude: f64, vel: [f64; 3]) -> [f64; 3] {
        let speed = vec3_norm(vel);
        if speed < 1e-10 {
            return [0.0; 3];
        }
        let rho = self.density(altitude);
        let f = 0.5 * rho * self.cd * self.area * speed * speed;
        let dir = vec3_normalize(vel);
        vec3_scale(dir, -f)
    }

    /// Secular decay rate da/dt \[m/s\] of a circular orbit of radius `a` \[m\].
    ///
    /// `altitude` \[m\] is the height at which the density is evaluated.  For a
    /// circular orbit the drag power `-½ ρ (Cd·A/m) v³` changes the orbital
    /// energy `-μ/(2a)`, which gives `da/dt = -ρ (Cd·A/m) √(μ a)`; this is the
    /// per-revolution loss `-2π ρ (Cd·A/m) a²` divided by the period.
    pub fn decay_rate(&self, altitude: f64, a: f64) -> f64 {
        let rho = self.density(altitude);
        -rho * self.cd * self.area / self.mass * (MU_EARTH * a).sqrt()
    }
}

// ─── LagrangePoints ───────────────────────────────────────────────────────────

/// Lagrange point calculator for the circular restricted 3-body problem (CR3BP).
///
/// Given primary masses m1 (larger) and m2 (smaller) in a circular orbit,
/// computes the five co-rotating equilibrium points L1–L5.
pub struct LagrangePoints {
    /// Mass of primary body \[kg\].
    pub m1: f64,
    /// Mass of secondary body \[kg\].
    pub m2: f64,
    /// Separation distance between primaries \[m\].
    pub distance: f64,
    /// Mass ratio μ = m2/(m1+m2).
    pub mass_ratio: f64,
}

impl LagrangePoints {
    /// Create a new Lagrange point calculator.
    pub fn new(m1: f64, m2: f64, distance: f64) -> Self {
        let mass_ratio = m2 / (m1 + m2);
        Self {
            m1,
            m2,
            distance,
            mass_ratio,
        }
    }

    /// Earth-Moon system Lagrange points.
    pub fn earth_moon() -> Self {
        let m1 = 5.972e24; // Earth mass [kg]
        let m2 = 7.342e22; // Moon mass [kg]
        Self::new(m1, m2, R_MOON)
    }

    /// Sun-Earth Lagrange points.
    pub fn sun_earth() -> Self {
        let m1 = 1.989e30; // Sun mass [kg]
        let m2 = 5.972e24; // Earth mass [kg]
        Self::new(m1, m2, 1.496e11) // 1 AU
    }

    /// L1 point: between the two primaries.
    ///
    /// Solved iteratively from the quintic polynomial.
    /// Returns distance from the secondary (m2) \[m\].
    pub fn l1_distance_from_secondary(&self) -> f64 {
        // Hill sphere radius approximation: r_L1 ≈ d * (μ/3)^(1/3)
        let mu = self.mass_ratio;
        let r_hill = self.distance * (mu / 3.0).powf(1.0 / 3.0);
        // Newton iteration on quintic
        let mut x = r_hill;
        for _ in 0..50 {
            let f = self.l1_quintic(x);
            let fp = self.l1_quintic_deriv(x);
            if fp.abs() < 1e-30 {
                break;
            }
            let dx = f / fp;
            x -= dx;
            if dx.abs() < 1e-12 * x.abs() {
                break;
            }
        }
        x
    }

    fn l1_quintic(&self, r: f64) -> f64 {
        // Collinear point equation for L1
        let mu = self.mass_ratio;
        let d = self.distance;
        let x = d * (1.0 - mu) - r; // L1 x position
        let _r1 = d - r; // distance from m1
        // Force balance: (1-mu)/(r1^2) - mu/(r^2) = ... (simplified)
        (1.0 - mu) / (x * x) - mu / (r * r) - x / (d * d * d)
    }

    fn l1_quintic_deriv(&self, r: f64) -> f64 {
        let mu = self.mass_ratio;
        let d = self.distance;
        let x = d * (1.0 - mu) - r;
        2.0 * (1.0 - mu) / (x * x * x) + 2.0 * mu / (r * r * r) - 1.0 / (d * d * d)
    }

    /// L1 position in the rotating frame \[m\].
    ///
    /// Origin at centre of mass, x-axis pointing from m1 to m2.
    pub fn l1(&self) -> [f64; 3] {
        let mu = self.mass_ratio;
        // L1 at x = 1 - mu - r_L1/d from m1
        let r_l1 = self.l1_distance_from_secondary();
        let x = self.distance * (1.0 - mu) - r_l1;
        [x, 0.0, 0.0]
    }

    /// L2 position (beyond the secondary) \[m\].
    pub fn l2(&self) -> [f64; 3] {
        let mu = self.mass_ratio;
        let r_hill = self.distance * (mu / 3.0).powf(1.0 / 3.0);
        let x = self.distance * (1.0 - mu) + r_hill;
        [x, 0.0, 0.0]
    }

    /// L3 position (opposite the secondary) \[m\].
    pub fn l3(&self) -> [f64; 3] {
        let mu = self.mass_ratio;
        // L3 at approximately x = -d*(1 + 5*mu/12)
        let x = -self.distance * (1.0 + 5.0 * mu / 12.0);
        [x, 0.0, 0.0]
    }

    /// L4 position (60° ahead of secondary) \[m\].
    pub fn l4(&self) -> [f64; 3] {
        let mu = self.mass_ratio;
        let x = self.distance * (0.5 - mu);
        let y = self.distance * (3.0_f64.sqrt() / 2.0);
        [x, y, 0.0]
    }

    /// L5 position (60° behind secondary) \[m\].
    pub fn l5(&self) -> [f64; 3] {
        let mu = self.mass_ratio;
        let x = self.distance * (0.5 - mu);
        let y = -self.distance * (3.0_f64.sqrt() / 2.0);
        [x, y, 0.0]
    }

    /// Check if L4/L5 are linearly stable (Routh criterion: μ < ~0.0385).
    pub fn is_triangular_stable(&self) -> bool {
        self.mass_ratio < 0.038_520_9
    }

    /// Distance between L4 (or L5) and each primary.
    ///
    /// For an equilateral triangle these should equal the primary separation.
    pub fn l4_distances(&self) -> (f64, f64) {
        let l4 = self.l4();
        let x1 = -self.distance * self.mass_ratio; // primary m1
        let x2 = self.distance * (1.0 - self.mass_ratio); // secondary m2
        let d1 = ((l4[0] - x1).powi(2) + l4[1].powi(2)).sqrt();
        let d2 = ((l4[0] - x2).powi(2) + l4[1].powi(2)).sqrt();
        (d1, d2)
    }
}

// ─── TleParser ────────────────────────────────────────────────────────────────

/// Parsed TLE (Two-Line Element) data.
#[derive(Debug, Clone)]
pub struct TleData {
    /// Satellite name.
    pub name: String,
    /// Satellite catalog number.
    pub catalog_number: u32,
    /// TLE epoch (year + day-of-year fraction).
    pub epoch_year: u32,
    /// Day of year with fractional part.
    pub epoch_day: f64,
    /// First time derivative of mean motion \[rev/day²\] × 2.
    pub mean_motion_dot: f64,
    /// Second time derivative of mean motion \[rev/day³\] × 6.
    pub mean_motion_ddot: f64,
    /// BSTAR drag term \[1/earth-radii\].
    pub bstar: f64,
    /// Inclination \[deg\].
    pub inclination: f64,
    /// Right ascension of ascending node \[deg\].
    pub raan: f64,
    /// Eccentricity (decimal, without leading 0.).
    pub eccentricity: f64,
    /// Argument of perigee \[deg\].
    pub arg_perigee: f64,
    /// Mean anomaly \[deg\].
    pub mean_anomaly: f64,
    /// Mean motion \[rev/day\].
    pub mean_motion: f64,
    /// Revolution number at epoch.
    pub rev_number: u32,
}

/// TLE parser and SGP4-simplified propagator.
///
/// Parses two-line element sets and propagates using a simplified SGP4
/// (without full atmospheric drag and luni-solar perturbations).
pub struct TleParser;

impl TleParser {
    /// Parse a TLE from three lines (name, line1, line2).
    ///
    /// Returns `None` if parsing fails.
    pub fn parse(name: &str, line1: &str, line2: &str) -> Option<TleData> {
        if line1.len() < 69 || line2.len() < 69 {
            return None;
        }
        if !line1.starts_with('1') || !line2.starts_with('2') {
            return None;
        }

        let catalog = line1[2..7].trim().parse::<u32>().ok()?;
        let epoch_str = line1[18..32].trim();
        let epoch_year: u32 = epoch_str[..2].trim().parse().ok()?;
        let epoch_day: f64 = epoch_str[2..].trim().parse().ok()?;

        let mean_motion_dot: f64 = line1[33..43].trim().parse().ok()?;
        let bstar_str = &line1[53..61];
        let bstar = Self::parse_decimal_with_exp(bstar_str);

        let inc: f64 = line2[8..16].trim().parse().ok()?;
        let raan: f64 = line2[17..25].trim().parse().ok()?;
        let ecc_str = line2[26..33].trim();
        let ecc: f64 = format!("0.{ecc_str}").parse().ok()?;
        let arg_peri: f64 = line2[34..42].trim().parse().ok()?;
        let mean_anom: f64 = line2[43..51].trim().parse().ok()?;
        let mean_mot: f64 = line2[52..63].trim().parse().ok()?;
        let rev: u32 = line2[63..68].trim().parse().unwrap_or(0);

        Some(TleData {
            name: name.trim().to_string(),
            catalog_number: catalog,
            epoch_year,
            epoch_day,
            mean_motion_dot,
            mean_motion_ddot: 0.0,
            bstar,
            inclination: inc,
            raan,
            eccentricity: ecc,
            arg_perigee: arg_peri,
            mean_anomaly: mean_anom,
            mean_motion: mean_mot,
            rev_number: rev,
        })
    }

    /// Parse the TLE "decimal with assumed leading 0." exponent format.
    /// E.g. " 16538-3" → 0.16538e-3 = 1.6538e-4.
    fn parse_decimal_with_exp(s: &str) -> f64 {
        let s = s.trim();
        if s.is_empty() {
            return 0.0;
        }
        // Find sign + last digit before possible exponent
        let neg = s.starts_with('-');
        let s2 = if neg || s.starts_with('+') {
            &s[1..]
        } else {
            s
        };
        // Find exponent separator (+ or - not at start)
        let mut split_idx = None;
        for (i, c) in s2.char_indices().skip(1) {
            if c == '+' || c == '-' {
                split_idx = Some(i);
                break;
            }
        }
        if let Some(idx) = split_idx {
            let mantissa: f64 = format!("0.{}", &s2[..idx]).parse().unwrap_or(0.0);
            let exp: i32 = s2[idx..].parse().unwrap_or(0);
            let val = mantissa * 10f64.powi(exp);
            if neg { -val } else { val }
        } else {
            let mantissa: f64 = format!("0.{s2}").parse().unwrap_or(0.0);
            if neg { -mantissa } else { mantissa }
        }
    }

    /// Propagate TLE data by `dt_minutes` using simplified Keplerian motion.
    ///
    /// Returns the orbital elements after propagation.
    pub fn propagate(tle: &TleData, dt_minutes: f64) -> OrbitalElements {
        let deg2rad = PI / 180.0;
        let n0 = tle.mean_motion * 2.0 * PI / 86400.0; // [rad/s]
        // Recover semi-major axis from mean motion
        let a = (MU_EARTH / (n0 * n0)).powf(1.0 / 3.0);
        let e = tle.eccentricity;
        let i = tle.inclination * deg2rad;
        let raan = tle.raan * deg2rad;
        let arg_p = tle.arg_perigee * deg2rad;
        let m0 = tle.mean_anomaly * deg2rad;

        let dt_sec = dt_minutes * 60.0;
        let m1 = (m0 + n0 * dt_sec).rem_euclid(2.0 * PI);
        let solver = KeplerOrbit::new(MU_EARTH);
        let e_anom = solver.solve_kepler(m1, e);
        let nu = solver.eccentric_to_true(e_anom, e);

        // J2 secular drift
        let j2pert = J2Perturbation::new_earth();
        let d_raan = j2pert.raan_drift_rate(a, e, i) * dt_sec;
        let d_w = j2pert.arg_perigee_drift_rate(a, e, i) * dt_sec;

        OrbitalElements::new(
            a,
            e,
            i,
            (raan + d_raan).rem_euclid(2.0 * PI),
            (arg_p + d_w).rem_euclid(2.0 * PI),
            nu,
        )
    }

    /// Verify TLE line checksum.
    pub fn verify_checksum(line: &str) -> bool {
        if line.len() < 69 {
            return false;
        }
        let mut sum: u32 = 0;
        for c in line[..68].chars() {
            match c {
                '0'..='9' => sum += c as u32 - '0' as u32,
                '-' => sum += 1,
                _ => {}
            }
        }
        let check = (sum % 10) as u8;
        let expected = line.chars().nth(68).unwrap_or('0') as u8 - b'0';
        check == expected
    }
}

// ─── HillClohessyWiltshire ─────────────────────────────────────────────────────

/// Hill-Clohessy-Wiltshire (CW) equations for relative motion in near-circular orbit.
///
/// Describes the relative motion of a chaser spacecraft with respect to a
/// target in a circular reference orbit.
pub struct HillClohessyWiltshire {
    /// Mean motion of the reference orbit \[rad/s\].
    pub n: f64,
}

impl HillClohessyWiltshire {
    /// Create a CW model for the given reference orbit radius \[m\].
    pub fn new(r_ref: f64) -> Self {
        Self {
            n: (MU_EARTH / r_ref.powi(3)).sqrt(),
        }
    }

    /// Propagate relative state `[x,y,z, xd,yd,zd]` by time `t` \[s\].
    ///
    /// Uses the exact CW solution (closed-form).
    /// Coordinate frame: x = radial, y = along-track, z = cross-track.
    pub fn propagate(&self, state0: [f64; 6], t: f64) -> [f64; 6] {
        let n = self.n;
        let (snt, cnt) = (n * t).sin_cos();
        let [x0, y0, z0, xd0, yd0, zd0] = state0;

        let x = (4.0 - 3.0 * cnt) * x0 + snt * xd0 / n + 2.0 * (1.0 - cnt) * yd0 / n;
        let y = 6.0 * (snt - n * t) * x0 + y0 - 2.0 * (1.0 - cnt) * xd0 / n
            + (4.0 * snt - 3.0 * n * t) * yd0 / n;
        let z = z0 * cnt + zd0 * snt / n;

        let xd = 3.0 * n * snt * x0 + cnt * xd0 + 2.0 * snt * yd0;
        let yd = 6.0 * n * (cnt - 1.0) * x0 - 2.0 * snt * xd0 + (4.0 * cnt - 3.0) * yd0;
        let zd = -z0 * n * snt + zd0 * cnt;

        [x, y, z, xd, yd, zd]
    }

    /// Compute the impulsive delta-v for a two-impulse rendezvous (CW targeting).
    ///
    /// Given the initial relative state `s0` and the desired final state `sf`
    /// (both `[x, y, z, xd, yd, zd]`), returns `Some((dv1, dv2))`.  `dv1` is
    /// applied at `t = 0` and puts the chaser on the coasting arc that reaches
    /// `sf[0..3]` at `tf`; `dv2` is applied at `tf` and turns the arrival
    /// velocity into `sf[3..6]`.  The existing velocity `s0[3..6]` is accounted
    /// for, and the solution is exact for the CW equations (see
    /// [`hcw_two_impulse`]).
    ///
    /// Returns `None` when the transfer is singular or ill-posed: `tf <= 0`,
    /// non-finite input, or a transfer time at which the CW position map
    /// cannot be inverted, such as `n·tf` a multiple of `π` (a whole number of
    /// half-periods; see [`hcw_two_impulse`] for the full list).
    pub fn two_impulse_rendezvous(
        &self,
        s0: [f64; 6],
        sf: [f64; 6],
        tf: f64,
    ) -> Option<([f64; 3], [f64; 3])> {
        hcw_two_impulse(
            self.n,
            [s0[0], s0[1], s0[2]],
            [s0[3], s0[4], s0[5]],
            [sf[0], sf[1], sf[2]],
            [sf[3], sf[4], sf[5]],
            tf,
        )
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng_util::uniform;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    const TOL: f64 = 1e-8;

    // ── Kepler's equation ─────────────────────────────────────────────────

    #[test]
    fn test_kepler_circular_orbit() {
        let solver = KeplerOrbit::new(MU_EARTH);
        // For e=0, E = M for all M
        for m in [0.0, 0.5, 1.0, 2.0, PI] {
            let e_anom = solver.solve_kepler(m, 0.0);
            assert!(
                (e_anom - m % (2.0 * PI)).abs() < 1e-10,
                "e_anom={e_anom} m={m}"
            );
        }
    }

    #[test]
    fn test_kepler_equation_solution() {
        let solver = KeplerOrbit::new(MU_EARTH);
        let e = 0.5;
        let m = 1.2;
        let e_anom = solver.solve_kepler(m, e);
        // Verify: M = E - e*sin(E)
        let m_check = e_anom - e * e_anom.sin();
        assert!(
            (m_check - m.rem_euclid(2.0 * PI)).abs() < 1e-10,
            "residual={}",
            (m_check - m.rem_euclid(2.0 * PI)).abs()
        );
    }

    #[test]
    fn test_kepler_high_eccentricity() {
        let solver = KeplerOrbit::new(MU_EARTH);
        let e = 0.9;
        for m in [0.1, 1.0, 2.5, 5.0] {
            let e_anom = solver.solve_kepler(m, e);
            let residual = (e_anom - e * e_anom.sin() - m.rem_euclid(2.0 * PI)).abs();
            assert!(residual < 1e-8, "high-e residual={residual} m={m}");
        }
    }

    #[test]
    fn test_vis_viva_circular() {
        let solver = KeplerOrbit::new(MU_EARTH);
        // Circular orbit: v = sqrt(mu/r), r = a → vis-viva = sqrt(mu*(2/r - 1/r)) = sqrt(mu/r)
        let r = 7_000_000.0;
        let v = solver.vis_viva_speed(r, r);
        let expected = (MU_EARTH / r).sqrt();
        assert!((v - expected).abs() < 1e-3, "v={v} expected={expected}");
    }

    #[test]
    fn test_anomaly_round_trip() {
        let solver = KeplerOrbit::new(MU_EARTH);
        let e = 0.3;
        let nu = 1.1; // true anomaly [rad]
        let e_anom = solver.true_to_eccentric(nu, e);
        let nu_back = solver.eccentric_to_true(e_anom, e);
        assert!((nu - nu_back).abs() < 1e-10, "round trip: {nu} → {nu_back}");
    }

    // ── Orbital elements ──────────────────────────────────────────────────

    #[test]
    fn test_orbital_elements_period() {
        let elems = OrbitalElements::new(7_000_000.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let t = elems.period(MU_EARTH);
        // LEO ~93 min period
        assert!(t > 5500.0 && t < 6000.0, "period={t}");
    }

    #[test]
    fn test_periapsis_apoapsis() {
        let a = 8_000_000.0;
        let e = 0.2;
        let elems = OrbitalElements::new(a, e, 0.0, 0.0, 0.0, 0.0);
        assert!((elems.periapsis() - a * (1.0 - e)).abs() < TOL, "periapsis");
        assert!((elems.apoapsis() - a * (1.0 + e)).abs() < TOL, "apoapsis");
    }

    #[test]
    fn test_state_vector_roundtrip() {
        let elems = OrbitalElements::new(7_500_000.0, 0.1, 0.52, 1.0, 0.5, 0.8);
        let (pos, vel) = elems.to_state_vector(MU_EARTH);
        let elems2 = OrbitalElements::from_state_vector(pos, vel, MU_EARTH);
        assert!(
            (elems.semi_major_axis - elems2.semi_major_axis).abs() / elems.semi_major_axis < 1e-6,
            "a: {} vs {}",
            elems.semi_major_axis,
            elems2.semi_major_axis
        );
        assert!(
            (elems.eccentricity - elems2.eccentricity).abs() < 1e-6,
            "e: {} vs {}",
            elems.eccentricity,
            elems2.eccentricity
        );
    }

    #[test]
    fn test_circular_orbit_velocity() {
        let r = 6_800_000.0;
        let elems = OrbitalElements::new(r, 0.0, 0.0, 0.0, 0.0, 0.0);
        let (_, vel) = elems.to_state_vector(MU_EARTH);
        let v = vec3_norm(vel);
        let expected = (MU_EARTH / r).sqrt();
        assert!(
            (v - expected).abs() / expected < 1e-6,
            "v={v} expected={expected}"
        );
    }

    // ── Orbital maneuvers ─────────────────────────────────────────────────

    #[test]
    fn test_hohmann_transfer_leo_geo() {
        let maneuver = OrbitalManeuver::new(MU_EARTH);
        let r_leo = 6_778_000.0; // ~400 km LEO
        let r_geo = 42_164_000.0; // GEO
        let (dv1, dv2, total) = maneuver.hohmann_transfer(r_leo, r_geo);
        assert!(dv1 > 0.0 && dv2 > 0.0, "delta-v must be positive");
        assert!((total - dv1 - dv2).abs() < 1e-10, "total = dv1 + dv2");
        // Expect total ~3.9 km/s for LEO to GEO
        assert!(total > 3500.0 && total < 4500.0, "total dv={total} m/s");
    }

    #[test]
    fn test_hohmann_tof() {
        let maneuver = OrbitalManeuver::new(MU_EARTH);
        let r1 = 6_778_000.0;
        let r2 = 42_164_000.0;
        let tof = maneuver.hohmann_tof(r1, r2);
        // ~5.25 hours for LEO-GEO
        assert!(tof > 18_000.0 && tof < 20_000.0, "tof={tof} s");
    }

    #[test]
    fn test_bielliptic_transfer() {
        let maneuver = OrbitalManeuver::new(MU_EARTH);
        let r1 = 6_778_000.0;
        let r2 = 42_164_000.0;
        let rb = 200_000_000.0; // very high intermediate orbit
        let (_dv1, _dv2, _dv3, total) = maneuver.bielliptic_transfer(r1, r2, rb);
        // Bielliptic with high rb should give similar or lower total than Hohmann
        assert!(total > 0.0, "bielliptic dv must be positive: {total}");
    }

    #[test]
    fn test_plane_change() {
        let maneuver = OrbitalManeuver::new(MU_EARTH);
        let r = 7_000_000.0;
        let dv = maneuver.plane_change(r, PI / 6.0); // 30° plane change
        let v = (MU_EARTH / r).sqrt();
        let expected = 2.0 * v * (PI / 12.0).sin();
        assert!((dv - expected).abs() < 1e-3, "dv={dv} expected={expected}");
    }

    #[test]
    fn test_tsiolkovsky_mass_ratio() {
        let maneuver = OrbitalManeuver::new(MU_EARTH);
        let dv = 3_000.0; // 3 km/s
        let isp = 300.0; // typical chemical engine
        let ratio = maneuver.tsiolkovsky_mass_ratio(dv, isp);
        assert!(ratio > 1.0, "mass ratio must be > 1");
        // dv = isp*g*ln(ratio) → ratio = exp(3000/(300*9.80665))
        let expected = (dv / (isp * 9.80665)).exp();
        assert!(
            (ratio - expected).abs() < 1e-10,
            "ratio={ratio} expected={expected}"
        );
    }

    // ── J2 perturbation ───────────────────────────────────────────────────

    #[test]
    fn test_j2_raan_drift_prograde() {
        let j2 = J2Perturbation::new_earth();
        let a = 7_000_000.0;
        let e = 0.0;
        let i = 0.5; // prograde orbit
        let drift = j2.raan_drift_rate(a, e, i);
        assert!(
            drift < 0.0,
            "prograde RAAN drift should be negative: {drift}"
        );
    }

    #[test]
    fn test_j2_raan_drift_polar() {
        let j2 = J2Perturbation::new_earth();
        let a = 7_000_000.0;
        let e = 0.0;
        let i = PI / 2.0; // polar orbit
        let drift = j2.raan_drift_rate(a, e, i);
        assert!(
            drift.abs() < 1e-10,
            "polar orbit RAAN drift should be zero: {drift}"
        );
    }

    #[test]
    fn test_j2_critical_inclination() {
        let i_crit = J2Perturbation::critical_inclination();
        // At critical inclination, 5*cos^2(i) - 1 = 0 → cos^2(i) = 1/5
        let check = 5.0 * i_crit.cos().powi(2) - 1.0;
        assert!(
            check.abs() < 1e-10,
            "critical inclination condition: {check}"
        );
    }

    // ── Atmospheric drag ──────────────────────────────────────────────────

    #[test]
    fn test_atmospheric_density_iss_altitude() {
        let drag = AtmosphericDrag::new(2.2, 1.0, 1.0);
        let rho = drag.density(400_000.0); // ISS at ~400 km
        assert!(rho > 0.0 && rho < 1e-8, "density at ISS: {rho}");
    }

    #[test]
    fn test_atmospheric_drag_decel_positive() {
        let drag = AtmosphericDrag::new(2.2, 10.0, 1000.0);
        let d = drag.drag_decel(400_000.0, 7700.0); // ISS speed
        assert!(d > 0.0, "drag decel must be positive: {d}");
    }

    #[test]
    fn test_atmospheric_density_decreases_with_altitude() {
        let drag = AtmosphericDrag::new(2.2, 1.0, 1.0);
        let rho200 = drag.density(200_000.0);
        let rho400 = drag.density(400_000.0);
        assert!(rho400 < rho200, "density should decrease with altitude");
    }

    #[test]
    fn test_density_matches_vallado_table_anchor_values() {
        let drag = AtmosphericDrag::new(2.2, 1.0, 1.0);
        // Band base values of Vallado (2013) Table 8-4.
        let anchors = [
            (0.0, 1.225),
            (100_000.0, 5.297e-7),
            (200_000.0, 2.789e-10),
            (300_000.0, 2.418e-11),
            (400_000.0, 3.725e-12),
            (500_000.0, 6.967e-13),
            (1_000_000.0, 3.019e-15),
        ];
        for (h, expected) in anchors {
            let rho = drag.density(h);
            assert!(
                ((rho - expected) / expected).abs() < 1e-12,
                "density at {h} m: {rho} vs {expected}"
            );
        }
    }

    #[test]
    fn test_density_inside_a_band_follows_the_exponential_law() {
        let drag = AtmosphericDrag::new(2.2, 1.0, 1.0);
        // 440 km sits 40 km into the 400 km band (H = 58.515 km).
        let expected = 3.725e-12 * (-40.0_f64 / 58.515).exp();
        let rho = drag.density(440_000.0);
        assert!(((rho - expected) / expected).abs() < 1e-12, "rho = {rho}");
    }

    #[test]
    fn test_density_is_continuous_across_every_band_edge() {
        let drag = AtmosphericDrag::new(2.2, 1.0, 1.0);
        for &(h0_km, _, _) in EXPONENTIAL_ATMOSPHERE.iter().skip(1) {
            let h = h0_km * 1_000.0;
            let below = drag.density(h - 1.0);
            let above = drag.density(h + 1.0);
            let jump = (below / above - 1.0).abs();
            assert!(
                jump < 5e-3,
                "density jumps by {jump} across the {h0_km} km band edge"
            );
        }
    }

    #[test]
    fn test_density_decreases_strictly_on_a_one_kilometre_grid() {
        let drag = AtmosphericDrag::new(2.2, 1.0, 1.0);
        let mut previous = drag.density(0.0);
        for km in 1..=1_500 {
            let rho = drag.density(f64::from(km) * 1_000.0);
            assert!(rho > 0.0 && rho < previous, "not decreasing at {km} km");
            previous = rho;
        }
    }

    #[test]
    fn test_density_at_leo_altitudes_is_physical() {
        let drag = AtmosphericDrag::new(2.2, 1.0, 1.0);
        let rho_400 = drag.density(400_000.0);
        assert!(
            (1e-12..=1e-11).contains(&rho_400),
            "rho(400 km) = {rho_400}"
        );
        // The mean thermosphere at 300 km is of the order of 1e-11 kg/m^3.
        let rho_300 = drag.density(300_000.0);
        assert!(
            (1e-11..=1e-10).contains(&rho_300),
            "rho(300 km) = {rho_300}"
        );
        // The old layered model reported ~9e-16 here.
        let rho_250 = drag.density(250_000.0);
        assert!(rho_250 > 5e-11, "rho(250 km) = {rho_250}");
    }

    #[test]
    fn test_density_clamps_below_surface_and_extends_above_the_table() {
        let drag = AtmosphericDrag::new(2.2, 1.0, 1.0);
        assert_eq!(drag.density(-5_000.0), drag.density(0.0));
        assert!(drag.density(f64::NAN).is_nan());
        let expected = 3.019e-15 * (-500.0_f64 / 268.0).exp();
        let rho = drag.density(1_500_000.0);
        assert!(((rho - expected) / expected).abs() < 1e-12, "rho = {rho}");
    }

    #[test]
    fn test_decay_rate_matches_circular_orbit_closed_form() {
        let drag = AtmosphericDrag::new(2.2, 1.0, 100.0);
        let a = 6_778.0e3;
        let altitude = a - R_EARTH;
        let rate = drag.decay_rate(altitude, a);
        let rho = drag.density(altitude);
        // da/dt = -rho (Cd A / m) sqrt(mu a)
        let expected = -rho * 0.022 * (MU_EARTH * a).sqrt();
        assert!(
            ((rate - expected) / expected).abs() < 1e-12,
            "rate = {rate}, expected = {expected}"
        );
        // Equivalent form: per-revolution loss 2 pi rho (Cd A / m) a^2 over the period.
        let period = 2.0 * PI * (a.powi(3) / MU_EARTH).sqrt();
        let per_revolution = -2.0 * PI * rho * 0.022 * a * a;
        assert!(
            ((rate - per_revolution / period) / rate).abs() < 1e-12,
            "rate = {rate}, per-revolution form = {}",
            per_revolution / period
        );
        // Independent route through the drag power: dE/dt = -a_drag v, E = -mu / (2a).
        let speed = (MU_EARTH / a).sqrt();
        let from_energy = -2.0 * a * a / MU_EARTH * drag.drag_decel(altitude, speed) * speed;
        assert!(
            ((rate - from_energy) / rate).abs() < 1e-12,
            "rate = {rate}, energy route = {from_energy}"
        );
    }

    #[test]
    fn test_decay_rate_is_metres_per_second_scale_in_low_earth_orbit() {
        // A 100 kg, 1 m^2 satellite at 400 km loses a few hundred metres per day
        // (the previous formula was off by a factor of 2 pi a, i.e. ~4e7).
        let drag = AtmosphericDrag::new(2.2, 1.0, 100.0);
        let a = R_EARTH + 400_000.0;
        let rate = drag.decay_rate(400_000.0, a);
        let per_day = rate * 86_400.0;
        assert!(per_day < 0.0, "orbit must decay: {per_day} m/day");
        assert!(
            per_day.abs() > 0.1 && per_day.abs() < 1_000.0,
            "decay per day = {per_day} m"
        );
    }

    // ── Lagrange points ───────────────────────────────────────────────────

    #[test]
    fn test_earth_moon_l1_location() {
        let lp = LagrangePoints::earth_moon();
        let l1 = lp.l1();
        // L1 is between Earth and Moon, closer to Moon (at ~85% of Moon distance)
        let x_moon = R_MOON * (1.0 - lp.mass_ratio);
        assert!(
            l1[0] > 0.0 && l1[0] < x_moon,
            "L1 x={} moon_x={}",
            l1[0],
            x_moon
        );
    }

    #[test]
    fn test_l4_l5_equilateral_distances() {
        let lp = LagrangePoints::earth_moon();
        let (d1, d2) = lp.l4_distances();
        // For equilateral triangle, both distances should equal R_MOON
        assert!(
            (d1 - R_MOON).abs() / R_MOON < 0.01,
            "d1={d1} R_MOON={R_MOON}"
        );
        assert!(
            (d2 - R_MOON).abs() / R_MOON < 0.01,
            "d2={d2} R_MOON={R_MOON}"
        );
    }

    #[test]
    fn test_l4_l5_symmetry() {
        let lp = LagrangePoints::earth_moon();
        let l4 = lp.l4();
        let l5 = lp.l5();
        // L4 and L5 are mirror images across x-axis
        assert!((l4[0] - l5[0]).abs() < 1e-3, "L4/L5 x should match");
        assert!((l4[1] + l5[1]).abs() < 1e-3, "L4/L5 y should be opposite");
    }

    #[test]
    fn test_earth_moon_stability() {
        let lp = LagrangePoints::earth_moon();
        // Earth-Moon mass ratio is ~1/81, well below Routh criterion
        assert!(
            lp.is_triangular_stable(),
            "Earth-Moon L4/L5 should be stable"
        );
    }

    #[test]
    fn test_sun_earth_l2_location() {
        let lp = LagrangePoints::sun_earth();
        let l2 = lp.l2();
        // L2 should be beyond Earth at ~1.5 million km from Earth
        let x_earth = lp.distance * (1.0 - lp.mass_ratio);
        assert!(
            l2[0] > x_earth,
            "L2 beyond Earth: l2={} earth_x={}",
            l2[0],
            x_earth
        );
    }

    // ── TLE parser ────────────────────────────────────────────────────────

    #[test]
    fn test_tle_parse_iss() {
        let name = "ISS (ZARYA)";
        let line1 = "1 25544U 98067A   23001.50000000  .00016717  00000-0  10270-3 0  9999";
        let line2 = "2 25544  51.6412 208.6184 0006400 287.5568  72.5123 15.50377579371143";
        let tle = TleParser::parse(name, line1, line2);
        assert!(tle.is_some(), "ISS TLE should parse successfully");
        let tle = tle.unwrap();
        assert_eq!(tle.catalog_number, 25544);
        assert!(
            (tle.inclination - 51.6412).abs() < 1e-3,
            "inc={}",
            tle.inclination
        );
        assert!(
            tle.mean_motion > 15.0 && tle.mean_motion < 16.0,
            "n={}",
            tle.mean_motion
        );
    }

    #[test]
    fn test_tle_propagation_runs() {
        let name = "ISS (ZARYA)";
        let line1 = "1 25544U 98067A   23001.50000000  .00016717  00000-0  10270-3 0  9999";
        let line2 = "2 25544  51.6412 208.6184 0006400 287.5568  72.5123 15.50377579371143";
        let tle = TleParser::parse(name, line1, line2).unwrap();
        let elems = TleParser::propagate(&tle, 90.0); // propagate 90 minutes
        assert!(
            elems.semi_major_axis > R_EARTH,
            "semi-major axis above Earth"
        );
        assert!(
            elems.eccentricity >= 0.0 && elems.eccentricity < 1.0,
            "valid eccentricity"
        );
    }

    #[test]
    fn test_tle_checksum_valid() {
        let line1 = "1 25544U 98067A   23001.50000000  .00016717  00000-0  10270-3 0  9999";
        // We test that the function runs without panic
        let _ = TleParser::verify_checksum(line1);
    }

    // ── Hill-Clohessy-Wiltshire ───────────────────────────────────────────

    #[test]
    fn test_cw_drift_free_initial_condition() {
        // For drift-free relative orbit: x0 = a, y0 = -2*a_dot/n, xd0 = 0, yd0 = 2*n*a
        let r_ref = 7_000_000.0;
        let cw = HillClohessyWiltshire::new(r_ref);
        let n = cw.n;
        let a = 1000.0; // 1 km along-track amplitude
        let state0 = [a, 0.0, 0.0, 0.0, 2.0 * n * a, 0.0];
        let t = 2.0 * PI / n; // one full orbit
        let state1 = cw.propagate(state0, t);
        // After one orbit, relative position should return to initial
        assert!(
            (state1[0] - state0[0]).abs() < 1.0,
            "x drift: {}",
            state1[0] - state0[0]
        );
        assert!(
            (state1[2] - state0[2]).abs() < 1.0,
            "z drift: {}",
            state1[2] - state0[2]
        );
    }

    #[test]
    fn test_cw_zero_state_stays_zero() {
        let cw = HillClohessyWiltshire::new(7_000_000.0);
        let state0 = [0.0; 6];
        let state1 = cw.propagate(state0, 1000.0);
        for (i, &s) in state1.iter().enumerate() {
            assert!(s.abs() < 1e-10, "component {i} should remain zero: {s}");
        }
    }

    #[test]
    fn test_cw_out_of_plane_oscillation() {
        let r_ref = 7_000_000.0;
        let cw = HillClohessyWiltshire::new(r_ref);
        let n = cw.n;
        let z0 = 500.0;
        let state0 = [0.0, 0.0, z0, 0.0, 0.0, 0.0];
        let t = PI / (2.0 * n); // quarter orbit
        let state1 = cw.propagate(state0, t);
        // z should be near zero at quarter period (cosine -> 0)
        assert!(
            state1[2].abs() < z0 * 0.01,
            "z at quarter period: {}",
            state1[2]
        );
    }

    #[test]
    fn test_two_impulse_rendezvous_reaches_the_requested_state() {
        let cw = HillClohessyWiltshire::new(7_000_000.0);
        let tf = 0.3 * 2.0 * PI / cw.n;
        let mut rng = StdRng::seed_from_u64(0x00C0_FFEE);
        for case in 0..20 {
            let r0: [f64; 3] = std::array::from_fn(|_| uniform(&mut rng, -2_000.0, 2_000.0));
            let v0: [f64; 3] = std::array::from_fn(|_| uniform(&mut rng, -3.0, 3.0));
            let rf: [f64; 3] = std::array::from_fn(|_| uniform(&mut rng, -2_000.0, 2_000.0));
            let vf: [f64; 3] = std::array::from_fn(|_| uniform(&mut rng, -3.0, 3.0));
            let s0 = [r0[0], r0[1], r0[2], v0[0], v0[1], v0[2]];
            let sf = [rf[0], rf[1], rf[2], vf[0], vf[1], vf[2]];
            let (dv1, dv2) = cw
                .two_impulse_rendezvous(s0, sf, tf)
                .expect("0.3 periods is a regular transfer time");

            let departure = [
                s0[0],
                s0[1],
                s0[2],
                s0[3] + dv1[0],
                s0[4] + dv1[1],
                s0[5] + dv1[2],
            ];
            let arrival = cw.propagate(departure, tf);
            for axis in 0..3 {
                assert!(
                    (arrival[axis] - sf[axis]).abs() < 1e-6,
                    "case {case}: position axis {axis} misses by {}",
                    arrival[axis] - sf[axis]
                );
                assert!(
                    (arrival[3 + axis] + dv2[axis] - sf[3 + axis]).abs() < 1e-8,
                    "case {case}: velocity axis {axis} is off by {}",
                    arrival[3 + axis] + dv2[axis] - sf[3 + axis]
                );
            }
        }
    }

    #[test]
    fn test_two_impulse_rendezvous_needs_no_burn_on_a_natural_trajectory() {
        // If the coasting arc from s0 already reaches sf, both impulses vanish.
        // The previous implementation ignored the initial velocity and could
        // not reproduce this.
        let cw = HillClohessyWiltshire::new(6_900_000.0);
        let tf = 1_700.0;
        let s0 = [350.0, -1_200.0, 80.0, 0.4, -0.9, 0.15];
        let sf = cw.propagate(s0, tf);
        let (dv1, dv2) = cw
            .two_impulse_rendezvous(s0, sf, tf)
            .expect("regular transfer time");
        for axis in 0..3 {
            assert!(dv1[axis].abs() < 1e-9, "dv1[{axis}] = {}", dv1[axis]);
            assert!(dv2[axis].abs() < 1e-9, "dv2[{axis}] = {}", dv2[axis]);
        }
    }

    #[test]
    fn test_two_impulse_rendezvous_rejects_singular_or_invalid_requests() {
        let cw = HillClohessyWiltshire::new(7_000_000.0);
        let period = 2.0 * PI / cw.n;
        let s0 = [100.0, -50.0, 20.0, 0.1, 0.0, -0.2];
        let sf = [0.0; 6];
        assert!(cw.two_impulse_rendezvous(s0, sf, 0.0).is_none());
        assert!(cw.two_impulse_rendezvous(s0, sf, -10.0).is_none());
        assert!(cw.two_impulse_rendezvous(s0, sf, f64::NAN).is_none());
        // n·tf = 2π and 4π: the in-plane position map is singular.
        assert!(cw.two_impulse_rendezvous(s0, sf, period).is_none());
        assert!(cw.two_impulse_rendezvous(s0, sf, 2.0 * period).is_none());
        // n·tf = π: the cross-track block sin(n·tf)/n vanishes.
        assert!(cw.two_impulse_rendezvous(s0, sf, 0.5 * period).is_none());
        let mut bad = s0;
        bad[1] = f64::INFINITY;
        assert!(cw.two_impulse_rendezvous(bad, sf, 0.3 * period).is_none());
        // The in-plane determinant 8(1 - cos φ) - 3 φ sin φ also vanishes at
        // φ ≈ 8.8387 (about 1.41 revolutions), where tan(φ/2) = 3φ/8.
        let first_inner_root = 8.838_742_844_152_04;
        assert!(
            cw.two_impulse_rendezvous(s0, sf, first_inner_root / cw.n)
                .is_none()
        );
        // Away from the singular times the transfer is solvable.
        assert!(cw.two_impulse_rendezvous(s0, sf, 0.9 * period).is_some());
        assert!(cw.two_impulse_rendezvous(s0, sf, 0.3 * period).is_some());
    }
}
