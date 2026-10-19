// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SPH-based sediment transport simulation module.
//!
//! This module provides a complete SPH framework for modelling sediment
//! transport phenomena in fluid–sediment coupled systems. It focuses on
//! particle-level algorithms that integrate directly with SPH fluid solvers.
//!
//! ## Features
//!
//! - **Bed load** via Shields criterion and Meyer-Peter-Muller formula
//! - **Suspended load** advection–diffusion with SPH discretisation
//! - **Erosion / deposition** exchange between bed and fluid phases
//! - **Morphological evolution** of the bed surface (Exner equation)
//! - **Sediment concentration** field tracked per SPH particle
//! - **Settling velocity** models (Stokes, Ferguson-Church, Richardson-Zaki)
//! - **Turbulent diffusion** of suspended sediment via SPH Laplacian
//! - **Scour modelling** around structures with local amplification
//! - **Multi-fraction** transport with hiding/exposure corrections
//! - **Cohesive sediment** flocculation and hindered settling

use std::collections::HashMap;
use std::f64::consts::PI;

use crate::neighbor::SpatialHash3D;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Gravitational acceleration magnitude (m/s^2).
const G: f64 = 9.81;

/// Density of quartz sediment grains (kg/m^3).
const RHO_S_DEFAULT: f64 = 2650.0;

/// Density of fresh water (kg/m^3).
const RHO_F_DEFAULT: f64 = 1000.0;

/// Kinematic viscosity of water at 20 C (m^2/s).
const NU_WATER: f64 = 1.0e-6;

/// Von Karman constant for turbulent boundary layers.
const KAPPA: f64 = 0.41;

/// Default critical Shields parameter for uniform sand.
const THETA_CR_DEFAULT: f64 = 0.047;

/// Meyer-Peter-Muller coefficient.
const MPM_COEFF: f64 = 8.0;

// ---------------------------------------------------------------------------
// Vector helpers (no nalgebra – use [f64; 3])
// ---------------------------------------------------------------------------

/// Dot product of two 3-vectors.
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean length of a 3-vector.
fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

/// Scale a 3-vector by scalar.
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Add two 3-vectors.
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract b from a.
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

// ---------------------------------------------------------------------------
// SPH kernel (cubic spline)
// ---------------------------------------------------------------------------

/// Gradient magnitude of the cubic-spline kernel in 3-D.
fn cubic_spline_dw(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma = 1.0 / (PI * h * h * h * h);
    if q < 1.0e-14 {
        0.0
    } else if q < 1.0 {
        sigma * (-3.0 * q + 2.25 * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        sigma * (-0.75 * t * t)
    } else {
        0.0
    }
}

/// Kernel gradient vector from particle i towards particle j.
fn grad_w(xi: [f64; 3], xj: [f64; 3], h: f64) -> [f64; 3] {
    let rij = sub3(xi, xj);
    let r = len3(rij);
    if r < 1.0e-30 {
        return [0.0; 3];
    }
    let dw = cubic_spline_dw(r, h);
    scale3(rij, dw / r)
}

// ---------------------------------------------------------------------------
// Settling velocity models
// ---------------------------------------------------------------------------

/// Stokes settling velocity for a spherical grain.
///
/// Valid for low Reynolds number (Re_p < 1).
///
/// # Arguments
/// * `d` - Grain diameter (m)
/// * `rho_s` - Sediment density (kg/m^3)
/// * `rho_f` - Fluid density (kg/m^3)
/// * `nu` - Kinematic viscosity (m^2/s)
pub fn settling_velocity_stokes(d: f64, rho_s: f64, rho_f: f64, nu: f64) -> f64 {
    let s = rho_s / rho_f;
    (s - 1.0) * G * d * d / (18.0 * nu)
}

/// Ferguson-Church settling velocity, valid across all grain Reynolds numbers.
///
/// Uses the universal drag relation of Ferguson & Church (2004):
///
/// `w_s = R g d^2 / (C1 nu + sqrt(0.75 C2 R g d^3))`,  `R = rho_s / rho_f - 1`,
/// with `C1 = 18` and `C2 = 1` for natural grains.  The formula itself lives in
/// [`crate::sediment_transport::settling_velocity_ferguson_church`], which this
/// function delegates to so that there is a single implementation.  It reduces
/// to Stokes' law for fine grains and to a `sqrt(d)` law for coarse ones.
///
/// A grain lighter than the fluid (`rho_s < rho_f`) rises: the result is then the
/// negative rise velocity, matching the sign convention of
/// [`settling_velocity_stokes`].  Non-positive diameters, viscosities or fluid
/// densities give `0.0`.
///
/// # Arguments
/// * `d` - Grain diameter (m)
/// * `rho_s` - Sediment density (kg/m^3)
/// * `rho_f` - Fluid density (kg/m^3)
/// * `nu` - Kinematic viscosity (m^2/s)
pub fn settling_velocity_ferguson_church(d: f64, rho_s: f64, rho_f: f64, nu: f64) -> f64 {
    if rho_f <= 0.0 {
        return 0.0;
    }
    if rho_s >= rho_f {
        crate::sediment_transport::settling_velocity_ferguson_church(d, rho_s, rho_f, nu)
    } else {
        // The drag law is symmetric in the submerged weight: mirror the density
        // contrast and flip the sign to obtain the rise velocity.
        -crate::sediment_transport::settling_velocity_ferguson_church(
            d,
            2.0 * rho_f - rho_s,
            rho_f,
            nu,
        )
    }
}

/// Richardson-Zaki hindered settling for concentrated suspensions.
///
/// `w = w_s0 * (1 - c)^n` where n depends on particle Reynolds number.
///
/// # Arguments
/// * `w_s0` - Single-particle settling velocity (m/s)
/// * `c` - Volumetric concentration (0..1)
/// * `n` - Richardson-Zaki exponent (typically 4.65 for Re_p < 0.2)
pub fn settling_velocity_hindered(w_s0: f64, c: f64, n: f64) -> f64 {
    if c >= 1.0 {
        return 0.0;
    }
    let c_clamp = c.clamp(0.0, 0.99);
    w_s0 * (1.0 - c_clamp).powf(n)
}

/// Dimensionless particle diameter D*.
///
/// `D* = d * ((s-1)*g / nu^2)^(1/3)`
pub fn dimensionless_diameter(d: f64, rho_s: f64, rho_f: f64, nu: f64) -> f64 {
    let s = rho_s / rho_f;
    d * ((s - 1.0) * G / (nu * nu)).powf(1.0 / 3.0)
}

// ---------------------------------------------------------------------------
// Shields criterion and bed-load initiation
// ---------------------------------------------------------------------------

/// Compute the Shields parameter (dimensionless bed shear stress).
///
/// `theta = tau_b / ((rho_s - rho_f) * g * d)`
pub fn shields_parameter(tau_b: f64, rho_s: f64, rho_f: f64, d: f64) -> f64 {
    let denom = (rho_s - rho_f) * G * d;
    if denom.abs() < 1.0e-30 {
        return 0.0;
    }
    tau_b / denom
}

/// Critical Shields parameter from the Soulsby-Whitehouse (1997) curve.
///
/// Applicable across a wide range of D*.
pub fn critical_shields_soulsby_whitehouse(d_star: f64) -> f64 {
    if d_star < 1.0e-10 {
        return 0.3;
    }
    0.3 / (1.0 + 1.2 * d_star) + 0.055 * (1.0 - (-0.02 * d_star).exp())
}

/// Check whether bed-load transport is initiated.
///
/// Returns `true` if the Shields parameter exceeds the critical value.
pub fn is_motion_initiated(theta: f64, theta_cr: f64) -> bool {
    theta > theta_cr
}

/// Bed shear stress from flow velocity using a quadratic friction law.
///
/// `tau_b = rho_f * Cf * |u|^2`
pub fn bed_shear_stress(rho_f: f64, cf: f64, u_mag: f64) -> f64 {
    rho_f * cf * u_mag * u_mag
}

/// Friction coefficient from Manning's roughness and depth.
///
/// `Cf = g * n^2 / h^(1/3)`
pub fn friction_coeff_manning(n_manning: f64, h: f64) -> f64 {
    if h < 1.0e-10 {
        return 0.0;
    }
    G * n_manning * n_manning / h.powf(1.0 / 3.0)
}

/// Friction velocity from bed shear stress.
///
/// `u* = sqrt(tau_b / rho_f)`
pub fn friction_velocity(tau_b: f64, rho_f: f64) -> f64 {
    if tau_b < 0.0 || rho_f < 1.0e-30 {
        return 0.0;
    }
    (tau_b / rho_f).sqrt()
}

// ---------------------------------------------------------------------------
// Bed-load transport formulae
// ---------------------------------------------------------------------------

/// Meyer-Peter-Muller bed-load transport rate (dimensionless).
///
/// `phi = 8 * (theta - theta_cr)^1.5`  for theta > theta_cr, else 0.
pub fn bedload_mpm(theta: f64, theta_cr: f64) -> f64 {
    if theta <= theta_cr {
        return 0.0;
    }
    MPM_COEFF * (theta - theta_cr).powf(1.5)
}

/// Convert dimensionless transport rate to volumetric transport rate.
///
/// `q_b = phi * sqrt((s-1) * g * d^3)`
pub fn bedload_volumetric(phi: f64, d: f64, rho_s: f64, rho_f: f64) -> f64 {
    let s = rho_s / rho_f;
    phi * ((s - 1.0) * G * d * d * d).sqrt()
}

/// Van Rijn (1984) bed-load transport formula (dimensionless).
///
/// Uses transport stage parameter T = (theta - theta_cr) / theta_cr.
pub fn bedload_van_rijn(theta: f64, theta_cr: f64, d_star: f64) -> f64 {
    if theta <= theta_cr || theta_cr < 1.0e-30 || d_star < 1.0e-30 {
        return 0.0;
    }
    let t_param = (theta - theta_cr) / theta_cr;
    0.053 * t_param.powf(2.1) / d_star.powf(0.3)
}

/// Einstein bed-load function (probability-based).
///
/// `phi = K * theta^1.5 * exp(-b / theta)` with empirical constants.
pub fn bedload_einstein(theta: f64) -> f64 {
    let k = 40.0;
    let b = 0.35;
    if theta < 1.0e-30 {
        return 0.0;
    }
    k * theta.powf(1.5) * (-b / theta).exp()
}

/// Slope correction factor for bed-load on inclined beds.
///
/// Adjusts critical Shields parameter for longitudinal slope angle.
pub fn slope_correction(theta_cr0: f64, bed_slope_angle: f64, repose_angle: f64) -> f64 {
    if repose_angle.abs() < 1.0e-10 {
        return theta_cr0;
    }
    let ratio = bed_slope_angle.sin() / repose_angle.sin();
    theta_cr0 * (1.0 - ratio).max(0.0)
}

// ---------------------------------------------------------------------------
// Suspended-load transport
// ---------------------------------------------------------------------------

/// Rouse number for suspended sediment.
///
/// `Z = w_s / (kappa * u*)`
pub fn rouse_number(w_s: f64, u_star: f64) -> f64 {
    if u_star.abs() < 1.0e-30 {
        return f64::MAX;
    }
    w_s / (KAPPA * u_star)
}

/// Rouse concentration profile.
///
/// Returns the sediment concentration at elevation z above the bed,
/// given reference concentration c_a at height a, and total depth h.
pub fn rouse_profile(z: f64, h: f64, a: f64, c_a: f64, rouse_z: f64) -> f64 {
    if z <= a || z >= h || a <= 0.0 || h <= a {
        return c_a;
    }
    let factor = ((h - z) / z) * (a / (h - a));
    if factor <= 0.0 {
        return 0.0;
    }
    c_a * factor.powf(rouse_z)
}

/// Depth-averaged suspended-sediment concentration by trapezoidal integration
/// of the Rouse profile.
pub fn rouse_depth_averaged(h: f64, a: f64, c_a: f64, rouse_z: f64, n_steps: usize) -> f64 {
    if h <= a || n_steps == 0 {
        return c_a;
    }
    let dz = (h - a) / n_steps as f64;
    let mut sum = 0.0;
    for i in 0..=n_steps {
        let z = a + i as f64 * dz;
        let c = rouse_profile(z, h, a, c_a, rouse_z);
        let w = if i == 0 || i == n_steps { 0.5 } else { 1.0 };
        sum += w * c;
    }
    sum * dz / (h - a)
}

/// Reference concentration at the bed from van Rijn (1984).
///
/// `c_a = 0.015 * (d/a) * T^1.5 / D*^0.3`
pub fn reference_concentration_van_rijn(d: f64, a: f64, t_param: f64, d_star: f64) -> f64 {
    if a < 1.0e-30 || d_star < 1.0e-30 {
        return 0.0;
    }
    0.015 * (d / a) * t_param.powf(1.5) / d_star.powf(0.3)
}

// ---------------------------------------------------------------------------
// Erosion and deposition
// ---------------------------------------------------------------------------

/// Partheniades erosion rate for cohesive sediment.
///
/// `E = M * (tau_b / tau_ce - 1)`  when tau_b > tau_ce.
pub fn erosion_partheniades(tau_b: f64, tau_ce: f64, m: f64) -> f64 {
    if tau_b <= tau_ce || tau_ce < 1.0e-30 {
        return 0.0;
    }
    m * (tau_b / tau_ce - 1.0)
}

/// Krone deposition rate for cohesive sediment.
///
/// `D = w_s * c * (1 - tau_b / tau_cd)`  when tau_b < tau_cd.
pub fn deposition_krone(w_s: f64, c: f64, tau_b: f64, tau_cd: f64) -> f64 {
    if tau_cd < 1.0e-30 || tau_b >= tau_cd {
        return 0.0;
    }
    w_s * c * (1.0 - tau_b / tau_cd)
}

/// Non-cohesive erosion rate (pick-up function, van Rijn style).
///
/// `E = alpha * w_s * T^1.5` where T is the transport-stage parameter.
pub fn erosion_noncohesive(w_s: f64, t_param: f64, alpha: f64) -> f64 {
    if t_param <= 0.0 {
        return 0.0;
    }
    alpha * w_s * t_param.powf(1.5)
}

/// Net vertical sediment flux (positive = erosion, negative = deposition).
pub fn net_sediment_flux(erosion: f64, deposition: f64) -> f64 {
    erosion - deposition
}

// ---------------------------------------------------------------------------
// Morphological evolution (Exner equation)
// ---------------------------------------------------------------------------

/// Bed porosity correction factor `1 / (1 - p)`.
pub fn porosity_factor(porosity: f64) -> f64 {
    if porosity >= 1.0 {
        return 1.0e10;
    }
    1.0 / (1.0 - porosity)
}

/// Exner equation: rate of change of bed elevation.
///
/// `dz_b/dt = -1/(1-p) * div(q_b) + (D - E)/(1-p)`
///
/// Here we take `div_qb` as a pre-computed divergence of bed-load flux
/// and net exchange = deposition - erosion.
pub fn exner_bed_change(div_qb: f64, deposition: f64, erosion: f64, porosity: f64) -> f64 {
    let pf = porosity_factor(porosity);
    pf * (-div_qb + deposition - erosion)
}

/// Bed elevation after one time step (explicit Euler).
pub fn bed_elevation_step(z_old: f64, dz_dt: f64, dt: f64) -> f64 {
    z_old + dz_dt * dt
}

/// Morphological acceleration factor.
///
/// Speeds up bed evolution relative to hydrodynamics.
pub fn morfac_scale(dz_dt: f64, morfac: f64) -> f64 {
    dz_dt * morfac
}

// ---------------------------------------------------------------------------
// Turbulent diffusion in SPH
// ---------------------------------------------------------------------------

/// Turbulent sediment diffusivity from eddy viscosity.
///
/// `epsilon_s = beta * epsilon_t` where beta ~ 1.0 for fine sediment.
pub fn sediment_diffusivity(epsilon_t: f64, beta: f64) -> f64 {
    beta * epsilon_t
}

/// Parabolic eddy viscosity profile.
///
/// `nu_t = kappa * u* * z * (1 - z/h)`
pub fn eddy_viscosity_parabolic(u_star: f64, z: f64, h: f64) -> f64 {
    if h < 1.0e-30 || z < 0.0 || z > h {
        return 0.0;
    }
    KAPPA * u_star * z * (1.0 - z / h)
}

/// SPH Laplacian approximation of diffusion for concentration field.
///
/// Returns the contribution to `dc_i/dt` from neighbour j.
/// Uses the Brookshaw (2003) formulation:
///   `(eps_i + eps_j) * (c_i - c_j) * (r_ij . grad_W) / (|r_ij|^2 + eta^2)`
///
/// with `r_ij = x_i - x_j` and `grad_W = grad_i W_ij`.  Because the kernel
/// decreases with distance, `r_ij . grad_W = r W'(r) <= 0`, so a particle that
/// is richer in sediment than its neighbour (`c_i > c_j`) loses it (negative
/// rate) and a poorer one gains it: diffusion flattens concentration peaks.
///
/// Multiply by `m_j / rho_j` and sum over all neighbours.
pub fn sph_diffusion_pair(
    c_i: f64,
    c_j: f64,
    eps_i: f64,
    eps_j: f64,
    xi: [f64; 3],
    xj: [f64; 3],
    h: f64,
    eta_sq: f64,
) -> f64 {
    let rij = sub3(xi, xj);
    let r2 = dot3(rij, rij);
    let gw = grad_w(xi, xj, h);
    let rdotgw = dot3(rij, gw);
    let eps_avg = eps_i + eps_j;
    eps_avg * (c_i - c_j) * rdotgw / (r2 + eta_sq)
}

/// Compute the full diffusion rate for particle i given all neighbours.
///
/// `dc_i/dt = sum_j (m_j/rho_j) * diffusion_pair(...)`
pub fn sph_diffusion_rate(
    _i: usize,
    positions: &[[f64; 3]],
    concentrations: &[f64],
    diffusivities: &[f64],
    masses: &[f64],
    densities: &[f64],
    neighbours: &[usize],
    h: f64,
    i_idx: usize,
) -> f64 {
    let eta_sq = 0.01 * h * h;
    let mut dc_dt = 0.0;
    for &j in neighbours {
        if j == i_idx {
            continue;
        }
        let pair = sph_diffusion_pair(
            concentrations[i_idx],
            concentrations[j],
            diffusivities[i_idx],
            diffusivities[j],
            positions[i_idx],
            positions[j],
            h,
            eta_sq,
        );
        dc_dt += (masses[j] / densities[j]) * pair;
    }
    dc_dt
}

// ---------------------------------------------------------------------------
// Sediment concentration tracking
// ---------------------------------------------------------------------------

/// A single SPH particle carrying sediment concentration.
#[derive(Debug, Clone)]
pub struct SedimentSphParticle {
    /// Position (m).
    pub pos: [f64; 3],
    /// Velocity (m/s).
    pub vel: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Density (kg/m^3).
    pub density: f64,
    /// Volumetric sediment concentration (0..1).
    pub concentration: f64,
    /// Turbulent diffusivity (m^2/s).
    pub diffusivity: f64,
    /// Smoothing length (m).
    pub h: f64,
    /// Flag: is this a bed particle?
    pub is_bed: bool,
}

impl SedimentSphParticle {
    /// Create a new fluid sediment particle.
    pub fn new_fluid(pos: [f64; 3], vel: [f64; 3], mass: f64, density: f64, h: f64) -> Self {
        Self {
            pos,
            vel,
            mass,
            density,
            concentration: 0.0,
            diffusivity: 1.0e-4,
            h,
            is_bed: false,
        }
    }

    /// Create a new bed particle.
    pub fn new_bed(pos: [f64; 3], mass: f64, density: f64, h: f64) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            mass,
            density,
            concentration: 0.6, // typical packed-bed concentration
            diffusivity: 0.0,
            h,
            is_bed: true,
        }
    }

    /// Speed magnitude.
    pub fn speed(&self) -> f64 {
        len3(self.vel)
    }

    /// Mixture density considering sediment.
    pub fn mixture_density(&self, rho_s: f64, rho_f: f64) -> f64 {
        self.concentration * rho_s + (1.0 - self.concentration) * rho_f
    }
}

// ---------------------------------------------------------------------------
// Scour modelling
// ---------------------------------------------------------------------------

/// Configuration for local scour around structures.
#[derive(Debug, Clone)]
pub struct ScourConfig {
    /// Amplification factor for bed shear stress near the structure.
    pub amplification: f64,
    /// Structure centre position.
    pub centre: [f64; 3],
    /// Structure characteristic dimension (e.g. pile diameter).
    pub diameter: f64,
    /// Influence radius (multiples of diameter).
    pub influence_radii: f64,
}

impl ScourConfig {
    /// Create a default scour config for a cylindrical pile.
    pub fn cylinder(centre: [f64; 3], diameter: f64) -> Self {
        Self {
            amplification: 2.5, // typical horseshoe vortex amplification
            centre,
            diameter,
            influence_radii: 3.0,
        }
    }

    /// Compute the local amplification factor at position x.
    ///
    /// Returns 1.0 outside the influence zone, linearly interpolates
    /// to `self.amplification` at the structure surface.
    pub fn local_amplification(&self, x: [f64; 3]) -> f64 {
        let dx = [x[0] - self.centre[0], x[1] - self.centre[1], 0.0];
        let r = len3(dx);
        let r_struct = self.diameter / 2.0;
        let r_inf = self.influence_radii * self.diameter;
        if r >= r_inf {
            1.0
        } else if r <= r_struct {
            self.amplification
        } else {
            let t = (r - r_struct) / (r_inf - r_struct);
            self.amplification + t * (1.0 - self.amplification)
        }
    }

    /// Equilibrium scour depth estimate (Melville & Coleman 2000, simplified).
    ///
    /// `d_se / D = 2.4 * K` where K is a shape/flow factor (K=1 for circular pile).
    pub fn equilibrium_depth(&self, k_factor: f64) -> f64 {
        2.4 * k_factor * self.diameter
    }

    /// Time-dependent scour depth using exponential approach.
    ///
    /// `d_s(t) = d_se * (1 - exp(-t / T_s))`
    pub fn scour_depth_time(&self, d_se: f64, t: f64, t_scale: f64) -> f64 {
        if t_scale < 1.0e-30 {
            return d_se;
        }
        d_se * (1.0 - (-t / t_scale).exp())
    }
}

// ---------------------------------------------------------------------------
// Multi-fraction transport
// ---------------------------------------------------------------------------

/// Grain fraction descriptor.
#[derive(Debug, Clone)]
pub struct GrainFraction {
    /// Representative diameter (m).
    pub diameter: f64,
    /// Density (kg/m^3).
    pub rho_s: f64,
    /// Fraction by mass in the bed (0..1).
    pub fraction: f64,
}

/// Hiding-exposure correction factor (Egiazaroff, 1965).
///
/// Modifies the critical Shields parameter for fraction i in a mixture.
///
/// `xi = (log10(19) / log10(19 * d_i / d_m))^2`
pub fn hiding_exposure_egiazaroff(d_i: f64, d_m: f64) -> f64 {
    if d_m < 1.0e-30 || d_i < 1.0e-30 {
        return 1.0;
    }
    let log19 = 19.0_f64.log10();
    let denom = (19.0 * d_i / d_m).log10();
    if denom.abs() < 1.0e-10 {
        return 1.0;
    }
    (log19 / denom).powi(2)
}

/// Compute bed-load transport for each fraction in a mixture.
///
/// Uses Meyer-Peter-Muller with hiding-exposure correction.
///
/// Returns a vector of dimensionless transport rates, one per fraction.
pub fn multi_fraction_bedload(
    fractions: &[GrainFraction],
    tau_b: f64,
    rho_f: f64,
    theta_cr0: f64,
) -> Vec<f64> {
    let d_m = mean_diameter(fractions);
    fractions
        .iter()
        .map(|f| {
            let xi = hiding_exposure_egiazaroff(f.diameter, d_m);
            let theta_cr_i = theta_cr0 * xi;
            let theta_i = shields_parameter(tau_b, f.rho_s, rho_f, f.diameter);
            let phi_i = bedload_mpm(theta_i, theta_cr_i);
            f.fraction * phi_i
        })
        .collect()
}

/// Mean diameter of a grain mixture.
pub fn mean_diameter(fractions: &[GrainFraction]) -> f64 {
    let total_frac: f64 = fractions.iter().map(|f| f.fraction).sum();
    if total_frac < 1.0e-30 {
        return 0.0;
    }
    fractions
        .iter()
        .map(|f| f.fraction * f.diameter)
        .sum::<f64>()
        / total_frac
}

/// Sorting coefficient (geometric std deviation) of a mixture.
pub fn sorting_coefficient(fractions: &[GrainFraction]) -> f64 {
    let d_m = mean_diameter(fractions);
    if d_m < 1.0e-30 || fractions.is_empty() {
        return 0.0;
    }
    let total_frac: f64 = fractions.iter().map(|f| f.fraction).sum();
    let var = fractions
        .iter()
        .map(|f| f.fraction * (f.diameter.ln() - d_m.ln()).powi(2))
        .sum::<f64>()
        / total_frac;
    var.sqrt().exp()
}

// ---------------------------------------------------------------------------
// Cohesive sediment / flocculation
// ---------------------------------------------------------------------------

/// Flocculation model parameters.
#[derive(Debug, Clone)]
pub struct FlocculationParams {
    /// Primary particle diameter (m).
    pub d_primary: f64,
    /// Maximum floc diameter (m).
    pub d_max: f64,
    /// Flocculation rate constant (1/s).
    pub k_floc: f64,
    /// Breakup rate constant (1/s).
    pub k_break: f64,
    /// Fractal dimension of flocs (typically 1.8-2.2).
    pub fractal_dim: f64,
}

impl FlocculationParams {
    /// Create default parameters for estuarine mud.
    pub fn estuarine_mud() -> Self {
        Self {
            d_primary: 4.0e-6,
            d_max: 2.0e-3,
            k_floc: 0.01,
            k_break: 0.005,
            fractal_dim: 2.0,
        }
    }

    /// Equilibrium floc diameter for given shear rate.
    ///
    /// `d_eq = d_max / (1 + (G / G_ref)^q)` where G is shear rate.
    pub fn equilibrium_floc_diameter(&self, shear_rate: f64) -> f64 {
        let g_ref = 10.0; // reference shear rate (1/s)
        let q = 0.5;
        self.d_max / (1.0 + (shear_rate / g_ref).powf(q))
    }

    /// Rate of change of floc diameter.
    ///
    /// `dd/dt = k_floc * (d_eq - d)` when d < d_eq (aggregation)
    /// `dd/dt = -k_break * (d - d_eq)` when d > d_eq (breakup)
    pub fn floc_growth_rate(&self, d_current: f64, d_eq: f64) -> f64 {
        if d_current < d_eq {
            self.k_floc * (d_eq - d_current)
        } else {
            -self.k_break * (d_current - d_eq)
        }
    }

    /// Effective density of a floc (less dense than primary grain).
    ///
    /// Uses fractal dimension approach:
    /// `rho_floc = rho_f + (rho_s - rho_f) * (d_p / d_floc)^(3 - D_f)`
    pub fn floc_density(&self, d_floc: f64, rho_s: f64, rho_f: f64) -> f64 {
        if d_floc < self.d_primary {
            return rho_s;
        }
        let ratio = self.d_primary / d_floc;
        rho_f + (rho_s - rho_f) * ratio.powf(3.0 - self.fractal_dim)
    }
}

// ---------------------------------------------------------------------------
// SPH sediment transport system
// ---------------------------------------------------------------------------

/// Configuration for the SPH sediment transport solver.
#[derive(Debug, Clone)]
pub struct SedimentSphConfig {
    /// Fluid density (kg/m^3).
    pub rho_f: f64,
    /// Default sediment grain density (kg/m^3).
    pub rho_s: f64,
    /// Representative grain diameter (m).
    pub d50: f64,
    /// Bed porosity (0..1).
    pub porosity: f64,
    /// Morphological acceleration factor.
    pub morfac: f64,
    /// Critical Shields parameter.
    pub theta_cr: f64,
    /// Manning roughness coefficient.
    pub manning_n: f64,
    /// Turbulent Schmidt number (epsilon_s / epsilon_t).
    pub schmidt_number: f64,
    /// Time step for morphological update (s).
    pub dt_morph: f64,
    /// Enable suspended load.
    pub enable_suspended: bool,
    /// Enable bed-load.
    pub enable_bedload: bool,
    /// Critical bed shear stress for deposition `tau_cd` (Pa) in the Krone law
    /// `D = w_s c (1 - tau_b / tau_cd)`.
    pub tau_cd: f64,
    /// Pick-up coefficient `alpha` of the non-cohesive erosion law
    /// `E = alpha w_s T^1.5` (`T` = transport-stage parameter).
    pub erosion_coeff: f64,
    /// Horizontal radius (m) within which neighbouring bed columns are coupled
    /// to estimate the bed-load flux divergence by least squares.  Choose about
    /// 1.5 column spacings so that the four face and four diagonal neighbours
    /// of a regular grid are included.
    pub column_radius: f64,
    /// Height above the bed (m) up to which fluid particles are sampled for the
    /// near-bed velocity and sediment concentration of a column.
    pub near_bed_height: f64,
}

impl Default for SedimentSphConfig {
    fn default() -> Self {
        Self {
            rho_f: RHO_F_DEFAULT,
            rho_s: RHO_S_DEFAULT,
            d50: 0.5e-3,
            porosity: 0.4,
            morfac: 1.0,
            theta_cr: THETA_CR_DEFAULT,
            manning_n: 0.025,
            schmidt_number: 1.0,
            dt_morph: 0.01,
            enable_suspended: true,
            enable_bedload: true,
            tau_cd: 0.1,
            erosion_coeff: 0.01,
            column_radius: 0.15,
            near_bed_height: 0.1,
        }
    }
}

/// State of a bed column for morphological tracking.
#[derive(Debug, Clone)]
pub struct BedColumn {
    /// Horizontal position \[x, y\].
    pub xy: [f64; 2],
    /// Current bed elevation (m).
    pub z_bed: f64,
    /// Cumulative erosion (m, positive means material removed).
    pub cumulative_erosion: f64,
    /// Cumulative deposition (m, positive means material added).
    pub cumulative_deposition: f64,
    /// Active layer thickness (m).
    pub active_layer: f64,
}

impl BedColumn {
    /// Create a new bed column at a given position and elevation.
    pub fn new(x: f64, y: f64, z_bed: f64) -> Self {
        Self {
            xy: [x, y],
            z_bed,
            cumulative_erosion: 0.0,
            cumulative_deposition: 0.0,
            active_layer: 0.01, // 1 cm default active layer
        }
    }

    /// Update bed elevation.
    pub fn update(&mut self, dz: f64) {
        if dz > 0.0 {
            self.cumulative_deposition += dz;
        } else {
            self.cumulative_erosion += dz.abs();
        }
        self.z_bed += dz;
    }
}

/// Local transport state of one bed column, evaluated in the first pass of
/// [`SedimentSphSolver::morphological_step`].
#[derive(Debug, Clone, Copy)]
struct ColumnTransport {
    /// Whether fluid was found near the bed above the column.
    active: bool,
    /// Bed-load flux vector `q_b` in the horizontal plane (m^2/s).
    bedload_flux: [f64; 2],
    /// Deposition rate `D` (m/s of bed elevation, before the porosity factor).
    deposition: f64,
    /// Erosion rate `E` (m/s of bed elevation, before the porosity factor).
    erosion: f64,
}

impl ColumnTransport {
    /// A column without fluid above it.
    fn inactive() -> Self {
        Self {
            active: false,
            bedload_flux: [0.0; 2],
            deposition: 0.0,
            erosion: 0.0,
        }
    }
}

/// Horizontal bucket grid over the fluid particles.
///
/// A bed column only sees fluid particles whose horizontal distance is below
/// their support radius `2 h_p`.  Bucketing the particles by `(x, y)` cell with
/// a cell size equal to the largest support radius lets each column visit the
/// 3 x 3 buckets around it instead of every particle.
struct FluidBuckets {
    cell_size: f64,
    cells: HashMap<(i64, i64), Vec<usize>>,
}

impl FluidBuckets {
    /// Bucket the fluid (non-bed) particles.
    fn build(particles: &[SedimentSphParticle]) -> Self {
        let support = particles
            .iter()
            .filter(|p| !p.is_bed)
            .map(|p| 2.0 * p.h)
            .filter(|s| s.is_finite())
            .fold(0.0_f64, f64::max);
        let cell_size = if support > 0.0 { support } else { 1.0 };
        let mut cells: HashMap<(i64, i64), Vec<usize>> = HashMap::new();
        for (idx, p) in particles.iter().enumerate() {
            if !p.is_bed {
                cells
                    .entry(Self::key(p.pos[0], p.pos[1], cell_size))
                    .or_default()
                    .push(idx);
            }
        }
        Self { cell_size, cells }
    }

    /// Bucket coordinates of the horizontal position `(x, y)`.
    fn key(x: f64, y: f64, cell_size: f64) -> (i64, i64) {
        (
            (x / cell_size).floor() as i64,
            (y / cell_size).floor() as i64,
        )
    }

    /// Indices of the fluid particles in the 3 x 3 buckets around `(x, y)`.
    fn candidates(&self, x: f64, y: f64) -> impl Iterator<Item = usize> + '_ {
        let (kx, ky) = Self::key(x, y, self.cell_size);
        (-1_i64..=1).flat_map(move |ox| {
            (-1_i64..=1).flat_map(move |oy| {
                let key = (kx.saturating_add(ox), ky.saturating_add(oy));
                self.cells.get(&key).into_iter().flatten().copied()
            })
        })
    }
}

/// Solve `M g = rhs` for the symmetric positive semi-definite 2 x 2 matrix
/// `M = [[a, b], [b, c]]` with the Moore-Penrose pseudo-inverse.
///
/// A (nearly) rank-deficient `M`, such as the moment matrix of columns on a
/// straight line, yields the minimum-norm solution, i.e. no variation across
/// the line.
fn solve_symmetric_pseudo_inverse(a: f64, b: f64, c: f64, rhs: [f64; 2]) -> [f64; 2] {
    let half_trace = 0.5 * (a + c);
    let radius = (0.25 * (a - c) * (a - c) + b * b).sqrt();
    let lambda_max = half_trace + radius;
    let lambda_min = half_trace - radius;
    if lambda_max.is_nan() || lambda_max <= 0.0 {
        return [0.0; 2];
    }
    // Unit eigenvector of the dominant eigenvalue.
    let v = if b != 0.0 {
        let (vx, vy) = (lambda_max - c, b);
        let norm = vx.hypot(vy);
        [vx / norm, vy / norm]
    } else if a >= c {
        [1.0, 0.0]
    } else {
        [0.0, 1.0]
    };
    let w = [-v[1], v[0]];
    let along_v = (v[0] * rhs[0] + v[1] * rhs[1]) / lambda_max;
    let mut g = [along_v * v[0], along_v * v[1]];
    if lambda_min > 1.0e-9 * lambda_max {
        let along_w = (w[0] * rhs[0] + w[1] * rhs[1]) / lambda_min;
        g[0] += along_w * w[0];
        g[1] += along_w * w[1];
    }
    g
}

/// Least-squares estimate of the divergence `div(q)` of a horizontal flux field.
///
/// Fits `q(x) ~ q_c + J (x - x_c)` to the neighbouring `(position, flux)`
/// samples in the least-squares sense around the `centre` sample and returns
/// `tr(J) = dq_x/dx + dq_y/dy`.  The fit is exact for linear fields; with
/// neighbours on a single line only the variation along that line is resolved,
/// and without neighbours the result is `0.0`.
fn least_squares_divergence(
    centre: ([f64; 2], [f64; 2]),
    neighbours: &[([f64; 2], [f64; 2])],
) -> f64 {
    let (x_c, q_c) = centre;
    let (mut m_xx, mut m_xy, mut m_yy) = (0.0_f64, 0.0_f64, 0.0_f64);
    let mut rhs_x = [0.0_f64; 2];
    let mut rhs_y = [0.0_f64; 2];
    for &(x, q) in neighbours {
        let d = [x[0] - x_c[0], x[1] - x_c[1]];
        let dq = [q[0] - q_c[0], q[1] - q_c[1]];
        m_xx += d[0] * d[0];
        m_xy += d[0] * d[1];
        m_yy += d[1] * d[1];
        rhs_x[0] += dq[0] * d[0];
        rhs_x[1] += dq[0] * d[1];
        rhs_y[0] += dq[1] * d[0];
        rhs_y[1] += dq[1] * d[1];
    }
    let grad_qx = solve_symmetric_pseudo_inverse(m_xx, m_xy, m_yy, rhs_x);
    let grad_qy = solve_symmetric_pseudo_inverse(m_xx, m_xy, m_yy, rhs_y);
    grad_qx[0] + grad_qy[1]
}

/// SPH sediment transport solver.
///
/// Manages the coupling between the SPH fluid and the movable bed.
#[derive(Debug, Clone)]
pub struct SedimentSphSolver {
    /// Configuration.
    pub config: SedimentSphConfig,
    /// Fluid particles.
    pub particles: Vec<SedimentSphParticle>,
    /// Bed columns for morphological tracking.
    pub bed_columns: Vec<BedColumn>,
    /// Current simulation time (s).
    pub time: f64,
}

impl SedimentSphSolver {
    /// Create a new solver with given configuration.
    pub fn new(config: SedimentSphConfig) -> Self {
        Self {
            config,
            particles: Vec::new(),
            bed_columns: Vec::new(),
            time: 0.0,
        }
    }

    /// Add a fluid particle to the system.
    pub fn add_fluid_particle(&mut self, pos: [f64; 3], vel: [f64; 3], mass: f64, h: f64) {
        self.particles.push(SedimentSphParticle::new_fluid(
            pos,
            vel,
            mass,
            self.config.rho_f,
            h,
        ));
    }

    /// Add a bed particle.
    pub fn add_bed_particle(&mut self, pos: [f64; 3], mass: f64, h: f64) {
        self.particles.push(SedimentSphParticle::new_bed(
            pos,
            mass,
            self.config.rho_s,
            h,
        ));
    }

    /// Add a bed column for morphological tracking.
    pub fn add_bed_column(&mut self, x: f64, y: f64, z_bed: f64) {
        self.bed_columns.push(BedColumn::new(x, y, z_bed));
    }

    /// Count fluid (non-bed) particles.
    pub fn fluid_count(&self) -> usize {
        self.particles.iter().filter(|p| !p.is_bed).count()
    }

    /// Count bed particles.
    pub fn bed_count(&self) -> usize {
        self.particles.iter().filter(|p| p.is_bed).count()
    }

    /// Compute bed shear stress at a given near-bed velocity magnitude and depth.
    pub fn compute_bed_shear(&self, u_near_bed: f64, depth: f64) -> f64 {
        let cf = friction_coeff_manning(self.config.manning_n, depth);
        bed_shear_stress(self.config.rho_f, cf, u_near_bed)
    }

    /// Compute bed-load transport rate (dimensionless) at a point.
    pub fn compute_bedload_rate(&self, tau_b: f64) -> f64 {
        let theta = shields_parameter(tau_b, self.config.rho_s, self.config.rho_f, self.config.d50);
        bedload_mpm(theta, self.config.theta_cr)
    }

    /// Single morphological time step for all bed columns.
    ///
    /// Advances every bed column with the Exner equation
    ///
    /// `(1 - p) dz_b/dt = -div(q_b) + D - E`
    ///
    /// (see [`exner_bed_change`]) in two passes, so that the update of one
    /// column never influences the fluxes seen by another within the step:
    ///
    /// 1. **Local transport** (per column).  The fluid particles whose
    ///    horizontal distance to the column is below their support radius
    ///    `2 h_p` and that lie above the bed give the water depth (height of
    ///    the highest such particle above the bed); those within
    ///    [`near_bed_height`](SedimentSphConfig::near_bed_height) give the
    ///    near-bed horizontal velocity `u` and sediment concentration `c` as
    ///    plain averages over the *same* set.  The Manning friction at that
    ///    depth gives the bed shear stress `tau_b`, then the Shields parameter
    ///    `theta`.  The bed-load flux vector is
    ///    `q_b = q_b(theta) u / |u|` (Meyer-Peter-Muller), the deposition rate
    ///    `D` follows the Krone law with
    ///    [`tau_cd`](SedimentSphConfig::tau_cd) and the erosion rate `E` the
    ///    van Rijn pick-up law with
    ///    [`erosion_coeff`](SedimentSphConfig::erosion_coeff), both using the
    ///    Ferguson-Church settling velocity of the grain.
    /// 2. **Flux divergence and bed update.**  `div(q_b)` of every column is
    ///    estimated from the fluxes of the columns within
    ///    [`column_radius`](SedimentSphConfig::column_radius) by a local
    ///    least-squares linear fit, and the bed changes by
    ///    `exner_bed_change(div q_b, D, E, porosity) * morfac * dt`.
    ///
    /// A uniform flow over a flat bed has a uniform flux, hence no divergence,
    /// and does not change the bed by bed-load.  Columns without fluid above
    /// them are left untouched and do not take part in the flux stencils.
    pub fn morphological_step(&mut self, dt: f64) {
        let transport = self.column_transport();
        let divergence = self.bedload_flux_divergence(&transport);
        let cfg = &self.config;
        for ((col, local), div_qb) in self
            .bed_columns
            .iter_mut()
            .zip(transport.iter())
            .zip(divergence.iter())
        {
            if !local.active {
                continue;
            }
            let dz_dt = exner_bed_change(*div_qb, local.deposition, local.erosion, cfg.porosity);
            col.update(morfac_scale(dz_dt, cfg.morfac) * dt);
        }
        self.time += dt;
    }

    /// Pass 1 of [`morphological_step`](Self::morphological_step): the local
    /// bed-load flux and the erosion / deposition rates of every bed column.
    fn column_transport(&self) -> Vec<ColumnTransport> {
        let cfg = &self.config;
        let buckets = FluidBuckets::build(&self.particles);
        let w_s = settling_velocity_ferguson_church(cfg.d50, cfg.rho_s, cfg.rho_f, NU_WATER);
        self.bed_columns
            .iter()
            .map(|col| {
                let mut u_sum = [0.0_f64; 2];
                let mut c_sum = 0.0_f64;
                let mut count = 0.0_f64;
                let mut depth = 0.0_f64;
                for idx in buckets.candidates(col.xy[0], col.xy[1]) {
                    let p = &self.particles[idx];
                    let dx = p.pos[0] - col.xy[0];
                    let dy = p.pos[1] - col.xy[1];
                    let height = p.pos[2] - col.z_bed;
                    if !(dx * dx + dy * dy < (2.0 * p.h).powi(2) && height > 0.0) {
                        continue;
                    }
                    depth = depth.max(height);
                    if height < cfg.near_bed_height {
                        u_sum[0] += p.vel[0];
                        u_sum[1] += p.vel[1];
                        c_sum += p.concentration;
                        count += 1.0;
                    }
                }
                if count < 1.0 {
                    return ColumnTransport::inactive();
                }
                let u_avg = [u_sum[0] / count, u_sum[1] / count];
                let speed = u_avg[0].hypot(u_avg[1]);
                let c_avg = c_sum / count;

                let cf = friction_coeff_manning(cfg.manning_n, depth);
                let tau_b = bed_shear_stress(cfg.rho_f, cf, speed);
                let theta = shields_parameter(tau_b, cfg.rho_s, cfg.rho_f, cfg.d50);

                let mut bedload_flux = [0.0_f64; 2];
                if cfg.enable_bedload && speed > 0.0 {
                    let phi = bedload_mpm(theta, cfg.theta_cr);
                    let q_mag = bedload_volumetric(phi, cfg.d50, cfg.rho_s, cfg.rho_f);
                    bedload_flux = [q_mag * u_avg[0] / speed, q_mag * u_avg[1] / speed];
                }

                let (deposition, erosion) = if cfg.enable_suspended {
                    let deposition = deposition_krone(w_s, c_avg, tau_b, cfg.tau_cd);
                    let erosion = if cfg.theta_cr > 0.0 && theta > cfg.theta_cr {
                        let stage = (theta - cfg.theta_cr) / cfg.theta_cr;
                        erosion_noncohesive(w_s, stage, cfg.erosion_coeff)
                    } else {
                        0.0
                    };
                    (deposition, erosion)
                } else {
                    (0.0, 0.0)
                };
                ColumnTransport {
                    active: true,
                    bedload_flux,
                    deposition,
                    erosion,
                }
            })
            .collect()
    }

    /// Pass 2 of [`morphological_step`](Self::morphological_step): the
    /// divergence of the bed-load flux at every column from a least-squares
    /// fit over the active columns within `column_radius`.
    fn bedload_flux_divergence(&self, transport: &[ColumnTransport]) -> Vec<f64> {
        let radius = self.config.column_radius;
        let n = self.bed_columns.len();
        if !(radius > 0.0 && radius.is_finite()) {
            return vec![0.0; n];
        }
        let anchors: Vec<[f64; 3]> = self
            .bed_columns
            .iter()
            .map(|col| [col.xy[0], col.xy[1], 0.0])
            .collect();
        let neighbour_lists = SpatialHash3D::all_neighbors(&anchors, radius);
        (0..n)
            .map(|i| {
                if !transport[i].active {
                    return 0.0;
                }
                let stencil: Vec<([f64; 2], [f64; 2])> = neighbour_lists[i]
                    .iter()
                    .filter(|&&j| transport[j].active)
                    .map(|&j| (self.bed_columns[j].xy, transport[j].bedload_flux))
                    .collect();
                least_squares_divergence(
                    (self.bed_columns[i].xy, transport[i].bedload_flux),
                    &stencil,
                )
            })
            .collect()
    }

    /// Total sediment volume eroded across all bed columns (m).
    pub fn total_erosion(&self) -> f64 {
        self.bed_columns.iter().map(|c| c.cumulative_erosion).sum()
    }

    /// Total sediment volume deposited across all bed columns (m).
    pub fn total_deposition(&self) -> f64 {
        self.bed_columns
            .iter()
            .map(|c| c.cumulative_deposition)
            .sum()
    }
}

// ---------------------------------------------------------------------------
// Exner equation grid (flux divergence form)
// ---------------------------------------------------------------------------

/// A regular 2-D Cartesian grid for solving the Exner bed-change equation.
///
/// The Exner equation in flux-divergence form is:
/// ```text
/// (1 - λ_p) ∂η/∂t + ∇·q_b = 0
/// ```
/// where `η` is the bed elevation, `λ_p` is the bed porosity, and `q_b` is
/// the bed-load flux vector.  Forward finite differences are used.
#[derive(Debug, Clone)]
pub struct ExnerGrid {
    /// Number of cells in the x-direction.
    pub grid_nx: usize,
    /// Number of cells in the y-direction.
    pub grid_ny: usize,
    /// Cell size in x (m).
    pub dx: f64,
    /// Cell size in y (m).
    pub dy: f64,
    /// Bed porosity λ_p (0..1).
    pub porosity: f64,
    /// Bed elevation η (m), stored in row-major order `[j * grid_nx + i]`.
    pub bed_elevation: Vec<f64>,
    /// x-component of bed-load flux q_bx at cell centres, row-major.
    pub flux_x: Vec<f64>,
    /// y-component of bed-load flux q_by at cell centres, row-major.
    pub flux_y: Vec<f64>,
}

impl ExnerGrid {
    /// Create a new flat Exner grid with zero fluxes.
    ///
    /// * `grid_nx`, `grid_ny` — number of cells in each direction.
    /// * `dx`, `dy` — cell dimensions (m).
    /// * `porosity` — bed porosity (typical ~0.4).
    /// * `initial_elevation` — uniform initial bed elevation (m).
    pub fn new(
        grid_nx: usize,
        grid_ny: usize,
        dx: f64,
        dy: f64,
        porosity: f64,
        initial_elevation: f64,
    ) -> Self {
        let n = grid_nx * grid_ny;
        Self {
            grid_nx,
            grid_ny,
            dx,
            dy,
            porosity,
            bed_elevation: vec![initial_elevation; n],
            flux_x: vec![0.0_f64; n],
            flux_y: vec![0.0_f64; n],
        }
    }

    /// Advance bed elevation by `dt` seconds using the Exner equation.
    ///
    /// Uses forward finite differences for the flux divergence.  At domain
    /// boundaries a backward difference is used to maintain a one-sided stencil.
    ///
    /// `∂η/∂t = -(∇·q_b) / (1 - λ_p)`
    pub fn compute_bed_change(&mut self, dt: f64) {
        let nx = self.grid_nx;
        let ny = self.grid_ny;
        let dx = self.dx;
        let dy = self.dy;
        let lambda_p = 1.0 - self.porosity;
        let lambda_p_safe = if lambda_p.abs() > 1e-30 {
            lambda_p
        } else {
            1.0
        };

        let mut deta = vec![0.0_f64; nx * ny];

        for j in 0..ny {
            for i in 0..nx {
                // Forward difference in x
                let dqbx_dx = if i + 1 < nx {
                    (self.flux_x[j * nx + i + 1] - self.flux_x[j * nx + i]) / dx
                } else {
                    // Backward difference at x-boundary
                    let i_prev = i.saturating_sub(1);
                    (self.flux_x[j * nx + i] - self.flux_x[j * nx + i_prev]) / dx
                };

                // Forward difference in y
                let dqby_dy = if j + 1 < ny {
                    (self.flux_y[(j + 1) * nx + i] - self.flux_y[j * nx + i]) / dy
                } else {
                    // Backward difference at y-boundary
                    let j_prev = j.saturating_sub(1);
                    (self.flux_y[j * nx + i] - self.flux_y[j_prev * nx + i]) / dy
                };

                let div_qb = dqbx_dx + dqby_dy;
                deta[j * nx + i] = -dt * div_qb / lambda_p_safe;
            }
        }

        for (eta, d) in self.bed_elevation.iter_mut().zip(deta.iter()) {
            *eta += d;
        }
    }
}

// ---------------------------------------------------------------------------
// Avalanche / angle of repose
// ---------------------------------------------------------------------------

/// Angle-of-repose avalanche correction.
///
/// If the local bed slope exceeds the angle of repose, redistribute
/// sediment downslope. Returns the correction dz for the uphill cell.
pub fn avalanche_correction(
    z_up: f64,
    z_down: f64,
    dx: f64,
    repose_angle: f64,
    relaxation: f64,
) -> f64 {
    let slope = (z_up - z_down) / dx;
    let slope_limit = repose_angle.tan();
    if slope <= slope_limit {
        return 0.0;
    }
    // Redistribute excess slope
    let excess = (slope - slope_limit) * dx;
    -relaxation * excess * 0.5
}

/// Apply avalanche correction to a 1-D array of bed elevations.
pub fn apply_avalanche_1d(bed: &mut [f64], dx: f64, repose_angle: f64, relaxation: f64) {
    let n = bed.len();
    if n < 2 {
        return;
    }
    let old = bed.to_vec();
    for i in 1..n {
        let dz = avalanche_correction(old[i - 1], old[i], dx, repose_angle, relaxation);
        bed[i - 1] += dz;
        bed[i] -= dz;
    }
    for i in (0..n - 1).rev() {
        let dz = avalanche_correction(old[i + 1], old[i], dx, repose_angle, relaxation);
        bed[i + 1] += dz;
        bed[i] -= dz;
    }
}

// ---------------------------------------------------------------------------
// Bed-form geometry
// ---------------------------------------------------------------------------

/// Ripple height estimate from van Rijn (1984).
///
/// `delta = d * 0.11 * (D*/20)^{-0.5} * (1 - exp(-T/2))`
pub fn ripple_height_van_rijn(d: f64, d_star: f64, t_param: f64) -> f64 {
    if d_star < 1.0e-10 || t_param < 0.0 {
        return 0.0;
    }
    d * 0.11 * (d_star / 20.0).powf(-0.5) * (1.0 - (-t_param / 2.0).exp())
}

/// Ripple wavelength (typically ~7 * ripple height).
pub fn ripple_wavelength(ripple_h: f64) -> f64 {
    7.0 * ripple_h
}

/// Dune height estimate (Yalin, 1964).
///
/// `H_d / h = 1/6 * (1 - Fr^2)` for `Fr < 1`.
pub fn dune_height_yalin(h: f64, froude: f64) -> f64 {
    if froude >= 1.0 {
        return 0.0;
    }
    h / 6.0 * (1.0 - froude * froude)
}

// ---------------------------------------------------------------------------
// Transport mode classification
// ---------------------------------------------------------------------------

/// Transport regime based on Rouse number.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TransportMode {
    /// Bed-load dominant (Z > 2.5).
    BedLoad,
    /// Suspended load with significant bed-load (1.2 < Z < 2.5).
    Mixed,
    /// Suspended load dominant (Z < 1.2).
    SuspendedLoad,
    /// Wash load (Z < 0.8).
    WashLoad,
}

/// Classify sediment transport mode from Rouse number.
pub fn classify_transport(rouse_z: f64) -> TransportMode {
    if rouse_z > 2.5 {
        TransportMode::BedLoad
    } else if rouse_z > 1.2 {
        TransportMode::Mixed
    } else if rouse_z > 0.8 {
        TransportMode::SuspendedLoad
    } else {
        TransportMode::WashLoad
    }
}

// ---------------------------------------------------------------------------
// Longshore transport (CERC formula)
// ---------------------------------------------------------------------------

/// CERC longshore sediment transport rate.
///
/// `Q_ls = K / (16 * (s-1) * (1-p)) * sqrt(g / gamma_b) * H_b^{5/2} * sin(2*alpha_b)`
///
/// where `H_b` is breaking wave height and `alpha_b` is breaking angle.
pub fn cerc_longshore_transport(
    k_cerc: f64,
    h_b: f64,
    alpha_b: f64,
    rho_s: f64,
    rho_f: f64,
    porosity: f64,
    gamma_b: f64,
) -> f64 {
    let s = rho_s / rho_f;
    let denom = 16.0 * (s - 1.0) * (1.0 - porosity);
    if denom.abs() < 1.0e-30 || gamma_b < 1.0e-30 {
        return 0.0;
    }
    k_cerc / denom * (G / gamma_b).sqrt() * h_b.powf(2.5) * (2.0 * alpha_b).sin()
}

// ---------------------------------------------------------------------------
// SPH-specific bed interaction kernel
// ---------------------------------------------------------------------------

/// Compute the bed interaction force on a fluid particle near the bed.
///
/// Models the drag/friction between a fluid particle and the fixed bed
/// using a penalty-like approach.
pub fn bed_interaction_force(
    particle_pos: [f64; 3],
    particle_vel: [f64; 3],
    bed_z: f64,
    bed_normal: [f64; 3],
    stiffness: f64,
    damping: f64,
) -> [f64; 3] {
    let penetration = bed_z - particle_pos[2];
    if penetration <= 0.0 {
        return [0.0; 3];
    }
    // Normal repulsion
    let f_normal = scale3(bed_normal, stiffness * penetration);
    // Tangential damping (friction)
    let vn = dot3(particle_vel, bed_normal);
    let f_damp = scale3(bed_normal, -damping * vn);
    add3(f_normal, f_damp)
}

/// Compute buoyancy-corrected weight for a sediment grain in fluid.
///
/// `F_b = (rho_s - rho_f) * V_grain * g`
pub fn submerged_weight(d: f64, rho_s: f64, rho_f: f64) -> f64 {
    let volume = PI / 6.0 * d * d * d;
    (rho_s - rho_f) * volume * G
}

// ---------------------------------------------------------------------------
// Concentration advection (SPH)
// ---------------------------------------------------------------------------

/// Advect sediment concentration forward in time for one particle.
///
/// Uses explicit Euler: `c_new = c_old + dt * (dc/dt_diffusion + dc/dt_settling)`
pub fn advect_concentration(c_old: f64, dc_dt_diff: f64, _w_s: f64, _dcdz: f64, dt: f64) -> f64 {
    let dc = dc_dt_diff * dt;
    let c_new = c_old + dc;
    c_new.clamp(0.0, 1.0)
}

/// Compute total suspended sediment load (integral of c over all fluid particles).
pub fn total_suspended_load(particles: &[SedimentSphParticle]) -> f64 {
    particles
        .iter()
        .filter(|p| !p.is_bed)
        .map(|p| p.concentration * p.mass / p.density)
        .sum()
}

/// Average concentration across all fluid particles.
pub fn average_concentration(particles: &[SedimentSphParticle]) -> f64 {
    let fluid: Vec<&SedimentSphParticle> = particles.iter().filter(|p| !p.is_bed).collect();
    if fluid.is_empty() {
        return 0.0;
    }
    let sum: f64 = fluid.iter().map(|p| p.concentration).sum();
    sum / fluid.len() as f64
}

// ---------------------------------------------------------------------------
// Density current model
// ---------------------------------------------------------------------------

/// Turbidity current front velocity (Ellison & Turner, 1959).
///
/// `U_f = Fr_d * sqrt(g' * h_c)` where g' is reduced gravity.
pub fn turbidity_front_velocity(g_prime: f64, h_current: f64, fr_d: f64) -> f64 {
    if g_prime < 0.0 || h_current < 0.0 {
        return 0.0;
    }
    fr_d * (g_prime * h_current).sqrt()
}

/// Reduced gravity from concentration.
///
/// `g' = g * (rho_s - rho_f) * c / rho_f`
pub fn reduced_gravity(c: f64, rho_s: f64, rho_f: f64) -> f64 {
    G * (rho_s - rho_f) * c / rho_f
}

/// Entrainment coefficient for turbidity current (Parker et al., 1986).
///
/// `E_w = 0.075 / sqrt(1 + 718 * Ri^2.4)` where Ri is Richardson number.
pub fn entrainment_coefficient(ri: f64) -> f64 {
    0.075 / (1.0 + 718.0 * ri.abs().powf(2.4)).sqrt()
}

// ---------------------------------------------------------------------------
// Bed armouring
// ---------------------------------------------------------------------------

/// Armour ratio: ratio of surface D50 to subsurface D50.
///
/// Values > 1 indicate a coarsened (armoured) surface.
pub fn armour_ratio(d50_surface: f64, d50_subsurface: f64) -> f64 {
    if d50_subsurface < 1.0e-30 {
        return 1.0;
    }
    d50_surface / d50_subsurface
}

/// Check if the bed is armoured (surface coarser than substrate).
pub fn is_armoured(d50_surface: f64, d50_subsurface: f64, threshold: f64) -> bool {
    armour_ratio(d50_surface, d50_subsurface) > threshold
}

// ---------------------------------------------------------------------------
// Utility: transport capacity
// ---------------------------------------------------------------------------

/// Sediment transport capacity (Engelund-Hansen total load).
///
/// `q_t = 0.05 * rho_f * u^5 / (g^2 * d50 * Delta^2 * Cf^{-1/2})`
pub fn transport_capacity_eh(u: f64, d50: f64, rho_s: f64, rho_f: f64, cf: f64) -> f64 {
    let delta = (rho_s - rho_f) / rho_f;
    if delta.abs() < 1.0e-30 || d50 < 1.0e-30 || cf < 1.0e-30 {
        return 0.0;
    }
    let g2 = G * G;
    0.05 * rho_f * u.powi(5) / (g2 * d50 * delta * delta) * cf.sqrt()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests;
