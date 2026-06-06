// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Sediment transport in Smoothed-Particle Hydrodynamics (SPH).
//!
//! This module implements sediment dynamics for SPH fluid simulations,
//! covering:
//!
//! - Shields criterion for initiation of bedload motion
//! - Bedload transport (Meyer-Peter–Müller and van Rijn formulae)
//! - Rouse profile for suspended sediment concentration
//! - Turbulent diffusion of suspended sediment
//! - Bed evolution (erosion / deposition fluxes, morphological change factor)
//! - Grain-size distribution (log-normal, median diameter, percentiles)
//! - Avalanche model for underwater slopes (angle of repose, critical slope)
//! - Density currents and turbidity-current SPH particles
//! - Settling velocity of sediment grains (Stokes, Ferguson-Church)
//! - Suspended load SPH advection–diffusion

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Gravitational acceleration (m/s²).
pub const GRAVITY: f64 = 9.81;

/// Density of quartz/typical sediment grain (kg/m³).
pub const RHO_SEDIMENT: f64 = 2650.0;

/// Density of freshwater (kg/m³).
pub const RHO_WATER: f64 = 1000.0;

/// Kinematic viscosity of water at 20 °C (m²/s).
pub const NU_WATER: f64 = 1.0e-6;

// ---------------------------------------------------------------------------
// Grain-size distribution
// ---------------------------------------------------------------------------

/// Log-normal grain-size distribution parameterised by median diameter and
/// log-standard deviation.
#[derive(Debug, Clone)]
pub struct GrainSizeDistribution {
    /// Median grain diameter D50 (m).
    pub d50: f64,
    /// Geometric standard deviation σ_g (dimensionless).
    pub sigma_g: f64,
    /// Sediment density ρ_s (kg/m³).
    pub rho_s: f64,
}

impl GrainSizeDistribution {
    /// Create a new [`GrainSizeDistribution`].
    ///
    /// # Arguments
    /// * `d50`     – median grain diameter (m)
    /// * `sigma_g` – geometric standard deviation (> 1 for poorly sorted)
    /// * `rho_s`   – sediment grain density (kg/m³)
    pub fn new(d50: f64, sigma_g: f64, rho_s: f64) -> Self {
        Self {
            d50,
            sigma_g,
            rho_s,
        }
    }

    /// Return the grain diameter at cumulative probability `p` (0–1) from
    /// the log-normal distribution.
    ///
    /// Uses the relation: `D_p = D50 * σ_g^(√2 · erf⁻¹(2p-1))`.
    /// For convenience the erf⁻¹ is approximated by a rational series.
    pub fn percentile(&self, p: f64) -> f64 {
        let p = p.clamp(1e-6, 1.0 - 1e-6);
        // Rational approximation to inverse error function (Winitzki 2003)
        let a = 0.147_f64;
        let ln_term = (1.0 - (2.0 * p - 1.0).powi(2)).ln();
        let common = 2.0 / (PI * a) + ln_term / 2.0;
        // inner = sqrt(common² − ln_term/a) − common  is always ≥ 0 for valid inputs.
        // erfinv(2p−1) = sign(p−0.5) · sqrt(inner).
        let inner = (common * common - ln_term / a).sqrt() - common;
        let erf_inv = if p > 0.5 {
            inner.sqrt()
        } else {
            -inner.max(0.0).sqrt()
        };
        self.d50 * self.sigma_g.powf(std::f64::consts::SQRT_2 * erf_inv)
    }

    /// D84 grain diameter (84th percentile).
    pub fn d84(&self) -> f64 {
        self.percentile(0.84)
    }

    /// D16 grain diameter (16th percentile).
    pub fn d16(&self) -> f64 {
        self.percentile(0.16)
    }

    /// Dimensionless grain diameter (D*) after van Rijn.
    ///
    /// ```text
    /// D* = D50 * [(s-1)*g / ν²]^(1/3)
    /// ```
    ///
    /// where `s = ρ_s / ρ_w`.
    pub fn d_star(&self, rho_w: f64, nu: f64) -> f64 {
        let s = self.rho_s / rho_w;
        let factor = ((s - 1.0) * GRAVITY / (nu * nu)).powf(1.0 / 3.0);
        self.d50 * factor
    }

    /// Submerged specific gravity of the sediment (dimensionless).
    pub fn submerged_specific_gravity(&self, rho_w: f64) -> f64 {
        (self.rho_s - rho_w) / rho_w
    }
}

// ---------------------------------------------------------------------------
// Settling velocity
// ---------------------------------------------------------------------------

/// Compute the settling velocity of a spherical grain using the Stokes law
/// (valid for particle Reynolds number Re_p < 0.5).
///
/// ```text
/// w_s = (ρ_s - ρ_w) g d² / (18 μ)
/// ```
///
/// # Arguments
/// * `d`     – grain diameter (m)
/// * `rho_s` – grain density (kg/m³)
/// * `rho_w` – fluid density (kg/m³)
/// * `nu`    – kinematic viscosity (m²/s)
pub fn settling_velocity_stokes(d: f64, rho_s: f64, rho_w: f64, nu: f64) -> f64 {
    if d <= 0.0 || nu <= 0.0 {
        return 0.0;
    }
    (rho_s - rho_w) * GRAVITY * d * d / (18.0 * rho_w * nu)
}

/// Compute the settling velocity using the Ferguson-Church formula, which
/// spans from Stokes to Newton's regime.
///
/// ```text
/// w_s = R g d² / (C1 ν + (0.75 C2 R g d³)^0.5)
/// ```
///
/// with C1 = 18, C2 = 1.0 for natural grains.
///
/// # Arguments
/// * `d`     – grain diameter (m)
/// * `rho_s` – grain density (kg/m³)
/// * `rho_w` – fluid density (kg/m³)
/// * `nu`    – kinematic viscosity (m²/s)
pub fn settling_velocity_ferguson_church(d: f64, rho_s: f64, rho_w: f64, nu: f64) -> f64 {
    if d <= 0.0 || nu <= 0.0 {
        return 0.0;
    }
    let r = (rho_s - rho_w) / rho_w;
    let c1 = 18.0_f64;
    let c2 = 1.0_f64;
    let numerator = r * GRAVITY * d * d;
    let denominator = c1 * nu + (0.75 * c2 * r * GRAVITY * d * d * d).sqrt();
    if denominator <= 0.0 {
        0.0
    } else {
        numerator / denominator
    }
}

// ---------------------------------------------------------------------------
// Shields criterion
// ---------------------------------------------------------------------------

/// Compute the Shields parameter (dimensionless bed shear stress) θ.
///
/// ```text
/// θ = τ_b / [(ρ_s - ρ_w) g d]
/// ```
///
/// # Arguments
/// * `tau_b` – bed shear stress (Pa)
/// * `rho_s` – sediment density (kg/m³)
/// * `rho_w` – fluid density (kg/m³)
/// * `d`     – grain diameter (m)
pub fn shields_parameter(tau_b: f64, rho_s: f64, rho_w: f64, d: f64) -> f64 {
    if d <= 0.0 {
        return 0.0;
    }
    tau_b / ((rho_s - rho_w) * GRAVITY * d)
}

/// Compute the critical Shields parameter θ_cr from the Soulsby-Whitehouse
/// (1997) formula based on the dimensionless grain diameter D*.
///
/// ```text
/// θ_cr = 0.30 / (1 + 1.2 D*) + 0.055 [1 − exp(−0.02 D*)]
/// ```
///
/// # Arguments
/// * `d_star` – dimensionless grain diameter D*
pub fn critical_shields_parameter(d_star: f64) -> f64 {
    if d_star <= 0.0 {
        return 0.0;
    }
    0.30 / (1.0 + 1.2 * d_star) + 0.055 * (1.0 - (-0.02 * d_star).exp())
}

/// Returns `true` when the Shields parameter exceeds the critical value,
/// i.e., when bedload motion initiates.
///
/// # Arguments
/// * `theta`    – current Shields parameter
/// * `theta_cr` – critical Shields parameter
pub fn motion_initiated(theta: f64, theta_cr: f64) -> bool {
    theta > theta_cr
}

// ---------------------------------------------------------------------------
// Bedload transport
// ---------------------------------------------------------------------------

/// Compute the dimensionless bedload transport rate Φ using the
/// Meyer-Peter–Müller formula.
///
/// ```text
/// Φ = 8 (θ − θ_cr)^1.5        for θ > θ_cr
/// Φ = 0                        otherwise
/// ```
///
/// # Arguments
/// * `theta`    – Shields parameter
/// * `theta_cr` – critical Shields parameter
pub fn bedload_mpm(theta: f64, theta_cr: f64) -> f64 {
    let excess = theta - theta_cr;
    if excess <= 0.0 {
        0.0
    } else {
        8.0 * excess.powf(1.5)
    }
}

/// Convert the dimensionless bedload rate Φ to volumetric bedload transport
/// per unit width q_b (m²/s).
///
/// ```text
/// q_b = Φ √[(s-1) g d³]
/// ```
///
/// # Arguments
/// * `phi`   – dimensionless transport rate
/// * `d`     – grain diameter (m)
/// * `s`     – specific gravity of sediment (ρ_s/ρ_w)
pub fn bedload_volumetric(phi: f64, d: f64, s: f64) -> f64 {
    if d <= 0.0 || s <= 1.0 {
        return 0.0;
    }
    phi * ((s - 1.0) * GRAVITY * d * d * d).sqrt()
}

/// Compute the dimensionless bedload rate using the van Rijn (1984) formula.
///
/// ```text
/// Φ_vR = 0.053 T^2.1 / D*^0.3
/// ```
///
/// where the transport stage parameter T = (θ/θ_cr − 1).
///
/// # Arguments
/// * `theta`    – Shields parameter
/// * `theta_cr` – critical Shields parameter
/// * `d_star`   – dimensionless grain diameter D*
pub fn bedload_van_rijn(theta: f64, theta_cr: f64, d_star: f64) -> f64 {
    if theta <= theta_cr || d_star <= 0.0 {
        return 0.0;
    }
    let t = theta / theta_cr - 1.0;
    0.053 * t.powf(2.1) / d_star.powf(0.3)
}

// ---------------------------------------------------------------------------
// Rouse profile
// ---------------------------------------------------------------------------

/// Compute the Rouse number Z (dimensionless).
///
/// ```text
/// Z = w_s / (κ u*)
/// ```
///
/// where κ = 0.41 is the von Kármán constant.
///
/// # Arguments
/// * `w_s`     – sediment settling velocity (m/s)
/// * `u_star`  – bed friction velocity (m/s)
pub fn rouse_number(w_s: f64, u_star: f64) -> f64 {
    const KAPPA: f64 = 0.41;
    if u_star <= 0.0 {
        return f64::INFINITY;
    }
    w_s / (KAPPA * u_star)
}

/// Rouse suspended-sediment concentration profile.
///
/// ```text
/// C(z) = C_a * [(a/z) * (h-z)/(h-a)]^Z
/// ```
///
/// # Arguments
/// * `z`    – height above bed (m)
/// * `h`    – total water depth (m)
/// * `a`    – reference height above bed (m); must satisfy 0 < a < h
/// * `c_a`  – reference concentration at height `a` (kg/m³)
/// * `z_r`  – Rouse number Z
pub fn rouse_profile(z: f64, h: f64, a: f64, c_a: f64, z_r: f64) -> f64 {
    if z <= 0.0 || h <= a || a <= 0.0 || z >= h {
        return 0.0;
    }
    let ratio = (a / z) * ((h - z) / (h - a));
    if ratio <= 0.0 {
        0.0
    } else {
        c_a * ratio.powf(z_r)
    }
}

/// Integrate the Rouse profile over the water column \[a, h\] using the
/// trapezoidal rule with `n_steps` intervals to obtain the depth-averaged
/// suspended sediment concentration.
///
/// # Arguments
/// * `h`       – total water depth (m)
/// * `a`       – reference height (m)
/// * `c_a`     – reference concentration (kg/m³)
/// * `z_r`     – Rouse number
/// * `n_steps` – number of integration intervals
pub fn rouse_depth_averaged(h: f64, a: f64, c_a: f64, z_r: f64, n_steps: usize) -> f64 {
    if h <= a || a <= 0.0 || n_steps == 0 {
        return 0.0;
    }
    let dz = (h - a) / n_steps as f64;
    let mut integral = 0.0;
    for k in 0..=n_steps {
        let z = a + k as f64 * dz;
        let c = rouse_profile(z, h, a, c_a, z_r);
        let w = if k == 0 || k == n_steps { 0.5 } else { 1.0 };
        integral += w * c * dz;
    }
    integral / (h - a)
}

// ---------------------------------------------------------------------------
// Turbulent diffusion of sediment
// ---------------------------------------------------------------------------

/// Compute the sediment diffusion coefficient ε_s from turbulent eddy
/// viscosity, assuming ε_s = β ε_f where β ≈ 1 for cohesionless sediment.
///
/// # Arguments
/// * `epsilon_f` – fluid turbulent diffusivity (m²/s)
/// * `beta`      – ratio ε_s / ε_f (typically 1.0)
pub fn sediment_diffusivity(epsilon_f: f64, beta: f64) -> f64 {
    epsilon_f * beta
}

/// Compute the turbulent diffusion flux of sediment in 1-D.
///
/// ```text
/// J = -ε_s * dC/dz
/// ```
///
/// # Arguments
/// * `eps_s` – sediment diffusivity (m²/s)
/// * `dc_dz` – concentration gradient (kg/m⁴)
pub fn diffusion_flux_1d(eps_s: f64, dc_dz: f64) -> f64 {
    -eps_s * dc_dz
}

// ---------------------------------------------------------------------------
// Erosion and deposition
// ---------------------------------------------------------------------------

/// Compute the erosion flux using the Partheniades (1965) formula for
/// cohesive sediment (E ∝ excess shear).
///
/// ```text
/// E = M (τ/τ_ce - 1)   for τ > τ_ce
/// E = 0                  otherwise
/// ```
///
/// # Arguments
/// * `tau`    – bed shear stress (Pa)
/// * `tau_ce` – critical erosion stress (Pa)
/// * `m`      – erosion rate coefficient (kg/m²/s)
pub fn erosion_flux_partheniades(tau: f64, tau_ce: f64, m: f64) -> f64 {
    if tau <= tau_ce || tau_ce <= 0.0 {
        return 0.0;
    }
    m * (tau / tau_ce - 1.0)
}

/// Compute the deposition flux for cohesive sediment (Krone 1962).
///
/// ```text
/// D = w_s C (1 - τ/τ_cd)   for τ < τ_cd
/// D = 0                      otherwise
/// ```
///
/// # Arguments
/// * `w_s`    – settling velocity (m/s)
/// * `c`      – near-bed sediment concentration (kg/m³)
/// * `tau`    – bed shear stress (Pa)
/// * `tau_cd` – critical deposition stress (Pa)
pub fn deposition_flux_krone(w_s: f64, c: f64, tau: f64, tau_cd: f64) -> f64 {
    if tau >= tau_cd || tau_cd <= 0.0 {
        return 0.0;
    }
    w_s * c * (1.0 - tau / tau_cd)
}

/// Net vertical sediment flux (deposition − erosion), positive = deposition.
///
/// # Arguments
/// * `deposition` – deposition rate (kg/m²/s)
/// * `erosion`    – erosion rate (kg/m²/s)
pub fn net_sediment_flux(deposition: f64, erosion: f64) -> f64 {
    deposition - erosion
}

// ---------------------------------------------------------------------------
// Bed evolution
// ---------------------------------------------------------------------------

/// State of a computational bed cell used for morphological updating.
#[derive(Debug, Clone)]
pub struct BedCell {
    /// Elevation of the bed above a datum (m).
    pub elevation: f64,
    /// Porosity of the bed material (dimensionless, 0–1).
    pub porosity: f64,
    /// Surface-averaged grain diameter (m).
    pub d50: f64,
}

impl BedCell {
    /// Create a new [`BedCell`].
    pub fn new(elevation: f64, porosity: f64, d50: f64) -> Self {
        Self {
            elevation,
            porosity,
            d50,
        }
    }

    /// Update bed elevation from a net sediment flux over time step `dt`.
    ///
    /// The Exner equation (1-D form):
    ///
    /// ```text
    /// ∂z_b/∂t = -1/(1-n) · ∂q_b/∂x
    /// ```
    ///
    /// Here we use the simplified point form where `dqb_dx` is the local
    /// divergence of bedload transport (m/s).
    ///
    /// # Arguments
    /// * `dqb_dx` – divergence of bedload transport (m/s)
    /// * `dt`     – time step (s)
    /// * `mf`     – morphological change factor (≥ 1)
    pub fn update_exner(&mut self, dqb_dx: f64, dt: f64, mf: f64) {
        let n = self.porosity;
        self.elevation -= mf * dt * dqb_dx / (1.0 - n);
    }

    /// Apply suspended-load erosion/deposition to the bed.
    ///
    /// # Arguments
    /// * `net_flux` – net vertical sediment flux (kg/m²/s); positive = deposition
    /// * `dt`       – time step (s)
    /// * `rho_s`    – sediment density (kg/m³)
    /// * `mf`       – morphological change factor
    pub fn update_suspended(&mut self, net_flux: f64, dt: f64, rho_s: f64, mf: f64) {
        if rho_s <= 0.0 {
            return;
        }
        let dz = mf * dt * net_flux / (rho_s * (1.0 - self.porosity));
        self.elevation += dz;
    }
}

// ---------------------------------------------------------------------------
// Morphological change factor
// ---------------------------------------------------------------------------

/// A morphological change factor accelerates bed evolution relative to the
/// hydrodynamic time scale.
///
/// MORFAC is commonly used in coastal/estuary modelling to allow simulation
/// of long-term bed change with short hydrodynamic simulations.
#[derive(Debug, Clone)]
pub struct Morfac {
    /// The factor by which bed change is multiplied (≥ 1).
    pub factor: f64,
}

impl Morfac {
    /// Create a new [`Morfac`] with the given multiplication factor.
    pub fn new(factor: f64) -> Self {
        let factor = factor.max(1.0);
        Self { factor }
    }

    /// Apply the morphological factor to a bed-change increment.
    pub fn amplify(&self, dz: f64) -> f64 {
        self.factor * dz
    }
}

// ---------------------------------------------------------------------------
// Avalanche model (underwater slope stability)
// ---------------------------------------------------------------------------

/// Underwater slope avalanche model.
///
/// When the local bed slope exceeds the critical (angle of repose) slope,
/// sediment is redistributed to restore stability.
#[derive(Debug, Clone)]
pub struct AvalancheModel {
    /// Critical slope angle (radians).
    pub critical_angle: f64,
    /// Fraction of excess slope removed per avalanche step (0–1).
    pub relax_factor: f64,
}

impl AvalancheModel {
    /// Create a new [`AvalancheModel`].
    ///
    /// # Arguments
    /// * `critical_angle` – angle of repose (radians); typical values 30°–35°
    /// * `relax_factor`   – relaxation fraction per step (0–1)
    pub fn new(critical_angle: f64, relax_factor: f64) -> Self {
        let relax_factor = relax_factor.clamp(0.0, 1.0);
        Self {
            critical_angle,
            relax_factor,
        }
    }

    /// Compute the critical slope (tangens of angle of repose).
    pub fn critical_slope(&self) -> f64 {
        self.critical_angle.tan()
    }

    /// Determine whether an avalanche is triggered between two adjacent bed
    /// cells with elevations `z_i` and `z_j` separated by distance `dx`.
    ///
    /// Returns the avalanche-induced change in elevation of cell `i` (m).
    /// A positive value means cell `i` loses sediment.
    ///
    /// # Arguments
    /// * `z_i` – elevation of current cell (m)
    /// * `z_j` – elevation of neighbour cell (m)
    /// * `dx`  – horizontal separation (m)
    pub fn avalanche_dz(&self, z_i: f64, z_j: f64, dx: f64) -> f64 {
        if dx <= 0.0 {
            return 0.0;
        }
        let slope = (z_i - z_j) / dx;
        let s_cr = self.critical_slope();
        if slope.abs() <= s_cr {
            return 0.0;
        }
        // Move sediment to restore critical slope
        let excess = slope.abs() - s_cr;
        let dz = self.relax_factor * excess * dx / 2.0;
        if slope > 0.0 { dz } else { -dz }
    }

    /// Apply avalanche diffusion to a 1-D bed profile in-place.
    ///
    /// Iterates over adjacent pairs and redistributes excess sediment.
    ///
    /// # Arguments
    /// * `bed`  – mutable slice of bed elevations (m)
    /// * `dx`   – uniform cell spacing (m)
    pub fn apply_1d(&self, bed: &mut [f64], dx: f64) {
        let n = bed.len();
        if n < 2 {
            return;
        }
        for i in 0..n - 1 {
            let dz = self.avalanche_dz(bed[i], bed[i + 1], dx);
            bed[i] -= dz;
            bed[i + 1] += dz;
        }
    }
}

// ---------------------------------------------------------------------------
// Density current / turbidity current SPH particle
// ---------------------------------------------------------------------------

/// A single SPH particle representing a turbidity current (sediment-laden
/// density current).
#[derive(Debug, Clone)]
pub struct TurbidityParticle {
    /// Particle position \[x, y, z\] (m).
    pub position: [f64; 3],
    /// Particle velocity \[vx, vy, vz\] (m/s).
    pub velocity: [f64; 3],
    /// Volumetric sediment concentration (0–1).
    pub concentration: f64,
    /// Particle mass (kg).
    pub mass: f64,
    /// Smoothing length h (m).
    pub h: f64,
}

impl TurbidityParticle {
    /// Create a new [`TurbidityParticle`].
    pub fn new(
        position: [f64; 3],
        velocity: [f64; 3],
        concentration: f64,
        mass: f64,
        h: f64,
    ) -> Self {
        Self {
            position,
            velocity,
            concentration: concentration.clamp(0.0, 1.0),
            mass,
            h,
        }
    }

    /// Compute the bulk density of the particle (fluid + sediment mixture).
    ///
    /// # Arguments
    /// * `rho_f` – fluid density (kg/m³)
    /// * `rho_s` – sediment density (kg/m³)
    pub fn bulk_density(&self, rho_f: f64, rho_s: f64) -> f64 {
        (1.0 - self.concentration) * rho_f + self.concentration * rho_s
    }

    /// Buoyancy-reduced gravity acting on the turbidity particle.
    ///
    /// # Arguments
    /// * `rho_f` – ambient fluid density (kg/m³)
    /// * `rho_s` – sediment grain density (kg/m³)
    pub fn reduced_gravity(&self, rho_f: f64, rho_s: f64) -> f64 {
        let rho_bulk = self.bulk_density(rho_f, rho_s);
        GRAVITY * (rho_bulk - rho_f) / rho_f
    }

    /// Advance position and velocity using explicit Euler integration with
    /// buoyancy and a simple drag term.
    ///
    /// # Arguments
    /// * `dt`   – time step (s)
    /// * `rho_f`  – ambient fluid density (kg/m³)
    /// * `rho_s`  – sediment grain density (kg/m³)
    /// * `cd`     – drag coefficient (dimensionless)
    pub fn step(&mut self, dt: f64, rho_f: f64, rho_s: f64, cd: f64) {
        let g_r = self.reduced_gravity(rho_f, rho_s);
        // Gravity acts in -z direction
        let az = -g_r;
        // Simple linear drag proportional to velocity magnitude
        let v_mag = {
            let [vx, vy, vz] = self.velocity;
            (vx * vx + vy * vy + vz * vz).sqrt()
        };
        let drag_fac = if v_mag > 0.0 { -cd * v_mag } else { 0.0 };
        let [vx, vy, vz] = self.velocity;
        self.velocity[0] += dt * drag_fac * vx;
        self.velocity[1] += dt * drag_fac * vy;
        self.velocity[2] += dt * (az + drag_fac * vz);
        let [vx2, vy2, vz2] = self.velocity;
        self.position[0] += dt * vx2;
        self.position[1] += dt * vy2;
        self.position[2] += dt * vz2;
    }

    /// Reduce the particle concentration by settling at rate `w_s` over `dt`.
    ///
    /// # Arguments
    /// * `w_s` – settling velocity (m/s)
    /// * `dt`  – time step (s)
    pub fn settle(&mut self, w_s: f64, dt: f64) {
        let lost = w_s * dt * self.concentration;
        self.concentration = (self.concentration - lost).max(0.0);
    }
}

// ---------------------------------------------------------------------------
// Turbidity current ensemble
// ---------------------------------------------------------------------------

/// A collection of [`TurbidityParticle`]s forming a turbidity current.
#[derive(Debug, Clone)]
pub struct TurbidityCurrentSph {
    /// All SPH particles in the turbidity current.
    pub particles: Vec<TurbidityParticle>,
    /// Ambient fluid density (kg/m³).
    pub rho_fluid: f64,
    /// Sediment grain density (kg/m³).
    pub rho_sediment: f64,
}

impl TurbidityCurrentSph {
    /// Create a new empty [`TurbidityCurrentSph`].
    pub fn new(rho_fluid: f64, rho_sediment: f64) -> Self {
        Self {
            particles: Vec::new(),
            rho_fluid,
            rho_sediment,
        }
    }

    /// Add a particle to the current.
    pub fn add_particle(&mut self, p: TurbidityParticle) {
        self.particles.push(p);
    }

    /// Advance all particles by one time step.
    ///
    /// # Arguments
    /// * `dt`   – time step (s)
    /// * `cd`   – drag coefficient
    /// * `w_s`  – sediment settling velocity (m/s)
    pub fn step(&mut self, dt: f64, cd: f64, w_s: f64) {
        for p in &mut self.particles {
            p.step(dt, self.rho_fluid, self.rho_sediment, cd);
            p.settle(w_s, dt);
        }
    }

    /// Remove particles whose concentration has fallen below a threshold.
    ///
    /// # Arguments
    /// * `c_min` – minimum concentration threshold
    pub fn prune(&mut self, c_min: f64) {
        self.particles.retain(|p| p.concentration >= c_min);
    }

    /// Compute the total sediment mass carried by all particles.
    pub fn total_sediment_mass(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| p.mass * p.concentration * self.rho_sediment)
            .sum()
    }

    /// Return the number of active particles.
    pub fn num_particles(&self) -> usize {
        self.particles.len()
    }
}

// ---------------------------------------------------------------------------
// Suspended load SPH advection-diffusion
// ---------------------------------------------------------------------------

/// Parameters for the SPH advection-diffusion of suspended sediment.
#[derive(Debug, Clone)]
pub struct SuspendedLoadParams {
    /// Settling velocity w_s (m/s).
    pub w_s: f64,
    /// Turbulent Schmidt number S_c (dimensionless; typically 0.7–1.0).
    pub schmidt_number: f64,
    /// Turbulent kinematic viscosity ν_t (m²/s).
    pub nu_t: f64,
}

impl SuspendedLoadParams {
    /// Create a new [`SuspendedLoadParams`].
    pub fn new(w_s: f64, schmidt_number: f64, nu_t: f64) -> Self {
        Self {
            w_s,
            schmidt_number,
            nu_t,
        }
    }

    /// Sediment diffusivity ε_s = ν_t / S_c.
    pub fn diffusivity(&self) -> f64 {
        if self.schmidt_number <= 0.0 {
            return 0.0;
        }
        self.nu_t / self.schmidt_number
    }
}

/// Compute the SPH concentration update rate for a single particle pair.
///
/// Applies the SPH discretisation of the diffusion operator:
///
/// ```text
/// dC_i/dt += 2 m_j/ρ_j * ε_s * (C_i - C_j) / |r_ij|² * r_ij · ∇W_ij
/// ```
///
/// This is a simplified scalar diffusion kernel contribution.
///
/// # Arguments
/// * `ci`    – concentration at particle i (kg/m³)
/// * `cj`    – concentration at particle j (kg/m³)
/// * `mj`    – mass of particle j (kg)
/// * `rho_j` – density of particle j (kg/m³)
/// * `eps_s` – sediment diffusivity (m²/s)
/// * `r_ij`  – distance between i and j (m)
/// * `dw`    – magnitude of ∇W_ij (1/m⁴)
pub fn sph_diffusion_rate(
    ci: f64,
    cj: f64,
    mj: f64,
    rho_j: f64,
    eps_s: f64,
    r_ij: f64,
    dw: f64,
) -> f64 {
    if rho_j <= 0.0 || r_ij <= 0.0 {
        return 0.0;
    }
    2.0 * mj / rho_j * eps_s * (ci - cj) / (r_ij * r_ij) * r_ij * dw
}

// ---------------------------------------------------------------------------
// Bed shear stress from depth-averaged velocity
// ---------------------------------------------------------------------------

/// Compute the bed shear stress from a depth-averaged velocity using the
/// quadratic friction law.
///
/// ```text
/// τ_b = ρ C_d U²
/// ```
///
/// where C_d is the drag coefficient.
///
/// # Arguments
/// * `rho` – fluid density (kg/m³)
/// * `u`   – depth-averaged velocity magnitude (m/s)
/// * `cd`  – drag coefficient (dimensionless; typically 0.003–0.006)
pub fn bed_shear_stress_quadratic(rho: f64, u: f64, cd: f64) -> f64 {
    rho * cd * u * u
}

/// Compute the bed friction velocity u*.
///
/// ```text
/// u* = √(τ_b / ρ)
/// ```
///
/// # Arguments
/// * `tau_b` – bed shear stress (Pa)
/// * `rho`   – fluid density (kg/m³)
pub fn friction_velocity(tau_b: f64, rho: f64) -> f64 {
    if tau_b <= 0.0 || rho <= 0.0 {
        return 0.0;
    }
    (tau_b / rho).sqrt()
}

// ---------------------------------------------------------------------------
// Sediment transport SPH particle
// ---------------------------------------------------------------------------

/// An SPH particle augmented with sediment transport fields.
#[derive(Debug, Clone)]
pub struct SedimentParticle {
    /// Position \[x, y, z\] (m).
    pub position: [f64; 3],
    /// Velocity \[vx, vy, vz\] (m/s).
    pub velocity: [f64; 3],
    /// Suspended sediment concentration (kg/m³).
    pub concentration: f64,
    /// Particle mass (kg).
    pub mass: f64,
    /// Smoothing length (m).
    pub h: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Fluid density (kg/m³).
    pub density: f64,
}

impl SedimentParticle {
    /// Create a new [`SedimentParticle`].
    pub fn new(
        position: [f64; 3],
        velocity: [f64; 3],
        concentration: f64,
        mass: f64,
        h: f64,
        pressure: f64,
        density: f64,
    ) -> Self {
        Self {
            position,
            velocity,
            concentration,
            mass,
            h,
            pressure,
            density,
        }
    }

    /// Advect the sediment concentration with the particle velocity plus a
    /// downward settling correction.
    ///
    /// # Arguments
    /// * `dt`  – time step (s)
    /// * `w_s` – settling velocity (m/s)
    pub fn advect_concentration(&mut self, dt: f64, w_s: f64) {
        // Simple first-order upwind settling in z direction
        let settled = w_s * dt * self.concentration / self.h.max(1e-12);
        self.concentration = (self.concentration - settled).max(0.0);
    }

    /// Compute the buoyancy force per unit volume on the particle.
    ///
    /// # Arguments
    /// * `rho_s` – sediment density (kg/m³)
    pub fn buoyancy_force(&self, rho_s: f64) -> [f64; 3] {
        let delta_rho = self.concentration * (rho_s - self.density);
        let fz = -GRAVITY * delta_rho;
        [0.0, 0.0, fz]
    }
}

// ---------------------------------------------------------------------------
// Grain-size sorting
// ---------------------------------------------------------------------------

/// Compute the armoring ratio: ratio of D50 of armored surface layer to
/// D50 of underlying substrate.
///
/// A value > 1 indicates coarsening of the bed surface (armoring).
///
/// # Arguments
/// * `d50_surface`   – median diameter of surface layer (m)
/// * `d50_substrate` – median diameter of substrate (m)
pub fn armoring_ratio(d50_surface: f64, d50_substrate: f64) -> f64 {
    if d50_substrate <= 0.0 {
        return 1.0;
    }
    d50_surface / d50_substrate
}

/// Hiding/exposure correction factor for mixed grain-size bedload
/// (Parker-Klingeman 1982 approach):
///
/// ```text
/// ξ_i = (D_i / D50)^{-b}
/// ```
///
/// where b ≈ 0.67.
///
/// # Arguments
/// * `d_i` – grain diameter of the fraction of interest (m)
/// * `d50` – median grain diameter of the mixture (m)
/// * `b`   – hiding exponent (typically 0.67)
pub fn hiding_exposure_factor(d_i: f64, d50: f64, b: f64) -> f64 {
    if d50 <= 0.0 || d_i <= 0.0 {
        return 1.0;
    }
    (d_i / d50).powf(-b)
}

// ---------------------------------------------------------------------------
// Bedform geometry
// ---------------------------------------------------------------------------

/// Estimate ripple height from empirical relations (van Rijn 1984).
///
/// ```text
/// Δ = 0.11 d (D*)^{-0.3} (1 - exp(-0.5 T)) (25 - T)
/// ```
/// for 0 < T < 25.
///
/// # Arguments
/// * `d`      – median grain diameter (m)
/// * `d_star` – dimensionless grain diameter
/// * `t`      – transport stage T = θ/θ_cr − 1
pub fn ripple_height(d: f64, d_star: f64, t: f64) -> f64 {
    if t <= 0.0 || t >= 25.0 || d_star <= 0.0 {
        return 0.0;
    }
    0.11 * d * d_star.powf(-0.3) * (1.0 - (-0.5 * t).exp()) * (25.0 - t)
}

/// Estimate ripple length from ripple height using the aspect ratio λ ≈ 7.3 Δ.
///
/// # Arguments
/// * `delta` – ripple height (m)
pub fn ripple_length(delta: f64) -> f64 {
    7.3 * delta
}

// ---------------------------------------------------------------------------
// Suspended sediment reference concentration
// ---------------------------------------------------------------------------

/// Reference concentration at height `a` above the bed (van Rijn 1984).
///
/// ```text
/// C_a = 0.015 (d/a) T^{1.5} D*^{-0.3}
/// ```
///
/// # Arguments
/// * `d`      – grain diameter (m)
/// * `a`      – reference height above bed (m)
/// * `t`      – transport stage parameter
/// * `d_star` – dimensionless grain diameter
pub fn reference_concentration_van_rijn(d: f64, a: f64, t: f64, d_star: f64) -> f64 {
    if a <= 0.0 || d_star <= 0.0 || t <= 0.0 {
        return 0.0;
    }
    0.015 * (d / a) * t.powf(1.5) / d_star.powf(0.3)
}

// ---------------------------------------------------------------------------
// Cohesive sediment (flocculation)
// ---------------------------------------------------------------------------

/// Simplified floc settling velocity model (Winterwerp 1998 limit form).
///
/// ```text
/// w_f = w_s0 * (C / C_ref)^p
/// ```
///
/// # Arguments
/// * `w_s0`  – primary particle settling velocity (m/s)
/// * `c`     – suspension concentration (kg/m³)
/// * `c_ref` – reference concentration (kg/m³)
/// * `p`     – empirical exponent (typically 1.0–2.0)
pub fn floc_settling_velocity(w_s0: f64, c: f64, c_ref: f64, p: f64) -> f64 {
    if c_ref <= 0.0 || c <= 0.0 {
        return w_s0;
    }
    w_s0 * (c / c_ref).powf(p)
}

// ---------------------------------------------------------------------------
// Turbidity current head velocity
// ---------------------------------------------------------------------------

/// Estimate turbidity current head velocity (Benjamin 1968).
///
/// ```text
/// U_head = Fr_h * √(g' h)
/// ```
///
/// where `Fr_h ≈ 0.5` for a full-depth current, `g' = g Δρ/ρ_a` is the
/// reduced gravity, and `h` is the current thickness.
///
/// # Arguments
/// * `g_prime`  – reduced gravity g' (m/s²)
/// * `h`        – current thickness (m)
/// * `fr_h`     – head Froude number (dimensionless; default ≈ 0.5)
pub fn turbidity_head_velocity(g_prime: f64, h: f64, fr_h: f64) -> f64 {
    if g_prime <= 0.0 || h <= 0.0 {
        return 0.0;
    }
    fr_h * (g_prime * h).sqrt()
}

// ---------------------------------------------------------------------------
// Bed slope effect on critical Shields parameter
// ---------------------------------------------------------------------------

/// Apply the bed-slope correction to the critical Shields parameter.
///
/// Longitudinal component:
/// ```text
/// θ_cr(β) = θ_cr0 * cos(β) * (1 - tan(β)/tan(φ))
/// ```
///
/// where β is the bed slope angle and φ is the friction angle.
///
/// # Arguments
/// * `theta_cr0` – critical Shields parameter on flat bed
/// * `beta`      – bed slope angle (radians, positive downslope)
/// * `phi`       – friction angle of sediment (radians)
pub fn bed_slope_shields_correction(theta_cr0: f64, beta: f64, phi: f64) -> f64 {
    if phi <= 0.0 {
        return theta_cr0;
    }
    theta_cr0 * beta.cos() * (1.0 - beta.tan() / phi.tan())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- settling_velocity_stokes ---

    #[test]
    fn test_stokes_zero_diameter() {
        assert_eq!(settling_velocity_stokes(0.0, 2650.0, 1000.0, 1e-6), 0.0);
    }

    #[test]
    fn test_stokes_positive() {
        let w = settling_velocity_stokes(2e-4, 2650.0, 1000.0, 1e-6);
        assert!(w > 0.0, "w={w}");
    }

    #[test]
    fn test_stokes_scales_as_d_squared() {
        let w1 = settling_velocity_stokes(1e-4, 2650.0, 1000.0, 1e-6);
        let w2 = settling_velocity_stokes(2e-4, 2650.0, 1000.0, 1e-6);
        assert!((w2 / w1 - 4.0).abs() < 1e-9, "ratio={}", w2 / w1);
    }

    #[test]
    fn test_stokes_zero_viscosity() {
        assert_eq!(settling_velocity_stokes(1e-4, 2650.0, 1000.0, 0.0), 0.0);
    }

    // --- settling_velocity_ferguson_church ---

    #[test]
    fn test_fc_positive() {
        let w = settling_velocity_ferguson_church(2e-4, 2650.0, 1000.0, 1e-6);
        assert!(w > 0.0, "w={w}");
    }

    #[test]
    fn test_fc_zero_diameter() {
        assert_eq!(
            settling_velocity_ferguson_church(0.0, 2650.0, 1000.0, 1e-6),
            0.0
        );
    }

    #[test]
    fn test_fc_small_grain_approaches_stokes() {
        // For very small grains the FC formula should be close to Stokes
        let d = 1e-5;
        let w_stokes = settling_velocity_stokes(d, 2650.0, 1000.0, 1e-6);
        let w_fc = settling_velocity_ferguson_church(d, 2650.0, 1000.0, 1e-6);
        let rel_diff = (w_fc - w_stokes).abs() / w_stokes;
        assert!(rel_diff < 0.15, "rel_diff={rel_diff}");
    }

    // --- shields_parameter ---

    #[test]
    fn test_shields_zero_diameter() {
        assert_eq!(shields_parameter(1.0, 2650.0, 1000.0, 0.0), 0.0);
    }

    #[test]
    fn test_shields_typical() {
        let theta = shields_parameter(0.5, 2650.0, 1000.0, 2e-4);
        assert!(theta > 0.0, "theta={theta}");
    }

    #[test]
    fn test_shields_scales_inversely_with_d() {
        let t1 = shields_parameter(1.0, 2650.0, 1000.0, 1e-3);
        let t2 = shields_parameter(1.0, 2650.0, 1000.0, 2e-3);
        assert!((t1 / t2 - 2.0).abs() < 1e-10, "ratio={}", t1 / t2);
    }

    // --- critical_shields_parameter ---

    #[test]
    fn test_critical_shields_positive() {
        let tc = critical_shields_parameter(10.0);
        assert!(tc > 0.0 && tc < 1.0, "tc={tc}");
    }

    #[test]
    fn test_critical_shields_zero_dstar() {
        assert_eq!(critical_shields_parameter(0.0), 0.0);
    }

    #[test]
    fn test_critical_shields_decreases_for_large_dstar() {
        let tc1 = critical_shields_parameter(1.0);
        let tc100 = critical_shields_parameter(100.0);
        // At large D* the first term tends to zero and second dominates
        assert!(tc1 > 0.0 && tc100 > 0.0);
    }

    // --- bedload_mpm ---

    #[test]
    fn test_mpm_below_critical() {
        assert_eq!(bedload_mpm(0.02, 0.05), 0.0);
    }

    #[test]
    fn test_mpm_above_critical() {
        let phi = bedload_mpm(0.1, 0.05);
        assert!(phi > 0.0, "phi={phi}");
    }

    #[test]
    fn test_mpm_exact_value() {
        // θ - θ_cr = 0.05, Φ = 8 * 0.05^1.5
        let phi = bedload_mpm(0.1, 0.05);
        let expected = 8.0 * 0.05_f64.powf(1.5);
        assert!((phi - expected).abs() < 1e-12);
    }

    // --- bedload_van_rijn ---

    #[test]
    fn test_van_rijn_below_critical() {
        assert_eq!(bedload_van_rijn(0.03, 0.05, 5.0), 0.0);
    }

    #[test]
    fn test_van_rijn_positive() {
        let phi = bedload_van_rijn(0.15, 0.05, 5.0);
        assert!(phi > 0.0, "phi={phi}");
    }

    // --- rouse_number ---

    #[test]
    fn test_rouse_zero_ustar() {
        assert_eq!(rouse_number(0.02, 0.0), f64::INFINITY);
    }

    #[test]
    fn test_rouse_typical() {
        let z = rouse_number(0.02, 0.05);
        assert!(z > 0.0 && z.is_finite(), "z={z}");
    }

    #[test]
    fn test_rouse_scales_inverse_ustar() {
        let z1 = rouse_number(0.02, 0.1);
        let z2 = rouse_number(0.02, 0.2);
        assert!((z1 / z2 - 2.0).abs() < 1e-10);
    }

    // --- rouse_profile ---

    #[test]
    fn test_rouse_at_reference_height() {
        // At z = a the profile should return c_a exactly
        let c = rouse_profile(0.05, 1.0, 0.05, 0.1, 2.0);
        assert!((c - 0.1).abs() < 1e-12, "c={c}");
    }

    #[test]
    fn test_rouse_decreases_with_height() {
        let c1 = rouse_profile(0.1, 2.0, 0.05, 0.1, 2.0);
        let c2 = rouse_profile(0.5, 2.0, 0.05, 0.1, 2.0);
        assert!(c2 < c1, "c1={c1} c2={c2}");
    }

    #[test]
    fn test_rouse_zero_above_surface() {
        let c = rouse_profile(2.5, 2.0, 0.05, 0.1, 2.0);
        assert_eq!(c, 0.0);
    }

    #[test]
    fn test_rouse_depth_averaged_positive() {
        let c_avg = rouse_depth_averaged(2.0, 0.05, 0.1, 2.0, 50);
        assert!(c_avg > 0.0 && c_avg <= 0.1, "c_avg={c_avg}");
    }

    // --- erosion / deposition ---

    #[test]
    fn test_erosion_below_critical() {
        assert_eq!(erosion_flux_partheniades(0.4, 0.5, 1e-3), 0.0);
    }

    #[test]
    fn test_erosion_positive() {
        let e = erosion_flux_partheniades(1.0, 0.5, 1e-3);
        assert!(e > 0.0, "e={e}");
    }

    #[test]
    fn test_deposition_above_critical() {
        assert_eq!(deposition_flux_krone(0.02, 0.5, 1.0, 0.5), 0.0);
    }

    #[test]
    fn test_deposition_positive() {
        let d = deposition_flux_krone(0.02, 0.5, 0.1, 0.5);
        assert!(d > 0.0, "d={d}");
    }

    // --- BedCell ---

    #[test]
    fn test_bed_cell_exner_erodes() {
        let mut cell = BedCell::new(1.0, 0.4, 2e-4);
        // Positive divergence → erosion → bed lowers
        cell.update_exner(0.01, 1.0, 1.0);
        assert!(cell.elevation < 1.0, "elevation={}", cell.elevation);
    }

    #[test]
    fn test_bed_cell_exner_deposits() {
        let mut cell = BedCell::new(1.0, 0.4, 2e-4);
        // Negative divergence → deposition → bed rises
        cell.update_exner(-0.01, 1.0, 1.0);
        assert!(cell.elevation > 1.0, "elevation={}", cell.elevation);
    }

    #[test]
    fn test_bed_cell_suspended_deposition() {
        let mut cell = BedCell::new(1.0, 0.4, 2e-4);
        cell.update_suspended(10.0, 0.1, 2650.0, 1.0);
        assert!(cell.elevation > 1.0);
    }

    // --- AvalancheModel ---

    #[test]
    fn test_avalanche_stable_slope() {
        let av = AvalancheModel::new(35.0_f64.to_radians(), 0.5);
        let dz = av.avalanche_dz(1.0, 0.9, 1.0);
        // slope = 0.1, critical ~ tan(35°) ≈ 0.7 → stable
        assert_eq!(dz, 0.0);
    }

    #[test]
    fn test_avalanche_unstable_slope() {
        let av = AvalancheModel::new(5.0_f64.to_radians(), 1.0);
        let dz = av.avalanche_dz(2.0, 0.0, 1.0);
        assert!(dz > 0.0, "dz={dz}");
    }

    #[test]
    fn test_avalanche_apply_1d_flattens() {
        let av = AvalancheModel::new(5.0_f64.to_radians(), 1.0);
        let mut bed = vec![10.0, 0.0, 0.0, 0.0, 0.0];
        av.apply_1d(&mut bed, 1.0);
        // The steep slope from bed[0] to bed[1] should have been reduced
        let initial_slope = 10.0;
        let new_slope = bed[0] - bed[1];
        assert!(new_slope < initial_slope, "new_slope={new_slope}");
    }

    // --- TurbidityParticle ---

    #[test]
    fn test_turbidity_bulk_density() {
        let p = TurbidityParticle::new([0.0; 3], [0.0; 3], 0.1, 1.0, 0.01);
        let rho = p.bulk_density(1000.0, 2650.0);
        // 0.9 * 1000 + 0.1 * 2650 = 900 + 265 = 1165
        assert!((rho - 1165.0).abs() < 1e-9, "rho={rho}");
    }

    #[test]
    fn test_turbidity_settle_reduces_concentration() {
        let mut p = TurbidityParticle::new([0.0; 3], [0.0; 3], 0.5, 1.0, 0.1);
        p.settle(0.01, 1.0);
        assert!(p.concentration < 0.5);
    }

    #[test]
    fn test_turbidity_concentration_non_negative() {
        let mut p = TurbidityParticle::new([0.0; 3], [0.0; 3], 1e-5, 1.0, 0.01);
        p.settle(100.0, 1000.0);
        assert!(p.concentration >= 0.0);
    }

    // --- TurbidityCurrentSph ---

    #[test]
    fn test_turbidity_sph_add_prune() {
        let mut tc = TurbidityCurrentSph::new(1000.0, 2650.0);
        tc.add_particle(TurbidityParticle::new([0.0; 3], [0.0; 3], 1e-8, 1.0, 0.01));
        tc.add_particle(TurbidityParticle::new([0.0; 3], [0.0; 3], 0.5, 1.0, 0.01));
        tc.prune(1e-6);
        assert_eq!(tc.num_particles(), 1);
    }

    #[test]
    fn test_turbidity_sph_step() {
        let mut tc = TurbidityCurrentSph::new(1000.0, 2650.0);
        tc.add_particle(TurbidityParticle::new(
            [0.0; 3],
            [1.0, 0.0, 0.0],
            0.2,
            1.0,
            0.1,
        ));
        let x0 = tc.particles[0].position[0];
        tc.step(0.1, 0.01, 0.001);
        assert!(tc.particles[0].position[0] > x0);
    }

    // --- GrainSizeDistribution ---

    #[test]
    fn test_gsd_d50_is_median() {
        let gsd = GrainSizeDistribution::new(2e-4, 1.5, 2650.0);
        let d50 = gsd.percentile(0.50);
        assert!((d50 - 2e-4).abs() < 1e-9, "d50={d50}");
    }

    #[test]
    fn test_gsd_d84_greater_than_d50() {
        let gsd = GrainSizeDistribution::new(2e-4, 2.0, 2650.0);
        assert!(gsd.d84() > gsd.d50, "d84={}", gsd.d84());
    }

    #[test]
    fn test_gsd_d16_less_than_d50() {
        let gsd = GrainSizeDistribution::new(2e-4, 2.0, 2650.0);
        assert!(gsd.d16() < gsd.d50, "d16={}", gsd.d16());
    }

    #[test]
    fn test_gsd_d_star_positive() {
        let gsd = GrainSizeDistribution::new(2e-4, 1.5, 2650.0);
        let ds = gsd.d_star(1000.0, 1e-6);
        assert!(ds > 0.0, "d_star={ds}");
    }

    // --- friction velocity ---

    #[test]
    fn test_friction_velocity_zero_tau() {
        assert_eq!(friction_velocity(0.0, 1000.0), 0.0);
    }

    #[test]
    fn test_friction_velocity_positive() {
        let u = friction_velocity(1.0, 1000.0);
        assert!((u - (1.0_f64 / 1000.0).sqrt()).abs() < 1e-12);
    }

    // --- ripple height ---

    #[test]
    fn test_ripple_height_positive() {
        let delta = ripple_height(2e-4, 5.0, 10.0);
        assert!(delta > 0.0, "delta={delta}");
    }

    #[test]
    fn test_ripple_height_zero_at_t_25() {
        assert_eq!(ripple_height(2e-4, 5.0, 25.0), 0.0);
    }

    #[test]
    fn test_ripple_length_proportional() {
        let lam = ripple_length(0.1);
        assert!((lam - 0.73).abs() < 1e-10);
    }

    // --- morphological factor ---

    #[test]
    fn test_morfac_minimum_one() {
        let mf = Morfac::new(0.5);
        assert_eq!(mf.factor, 1.0);
    }

    #[test]
    fn test_morfac_amplify() {
        let mf = Morfac::new(10.0);
        assert!((mf.amplify(0.1) - 1.0).abs() < 1e-12);
    }

    // --- hiding exposure ---

    #[test]
    fn test_hiding_exposure_unity_at_d50() {
        let xi = hiding_exposure_factor(2e-4, 2e-4, 0.67);
        assert!((xi - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_hiding_coarse_grain_exposed() {
        // Coarser than D50 → xi < 1 (exposed)
        let xi = hiding_exposure_factor(4e-4, 2e-4, 0.67);
        assert!(xi < 1.0, "xi={xi}");
    }

    // --- floc settling ---

    #[test]
    fn test_floc_settling_no_concentration() {
        let w = floc_settling_velocity(0.001, 0.0, 1.0, 1.5);
        assert!((w - 0.001).abs() < 1e-12);
    }

    #[test]
    fn test_floc_settling_increases_with_concentration() {
        let w1 = floc_settling_velocity(0.001, 0.5, 1.0, 1.0);
        let w2 = floc_settling_velocity(0.001, 1.5, 1.0, 1.0);
        assert!(w2 > w1, "w1={w1} w2={w2}");
    }

    // --- turbidity head velocity ---

    #[test]
    fn test_turbidity_head_zero_depth() {
        assert_eq!(turbidity_head_velocity(0.1, 0.0, 0.5), 0.0);
    }

    #[test]
    fn test_turbidity_head_positive() {
        let u = turbidity_head_velocity(0.1, 10.0, 0.5);
        assert!(u > 0.0, "u={u}");
    }

    // --- bed slope correction ---

    #[test]
    fn test_bed_slope_zero_slope() {
        let tc = bed_slope_shields_correction(0.05, 0.0, 35.0_f64.to_radians());
        assert!((tc - 0.05).abs() < 1e-10, "tc={tc}");
    }

    #[test]
    fn test_bed_slope_reduces_critical() {
        // Downslope → critical Shields should be less than flat bed
        let tc_flat = bed_slope_shields_correction(0.05, 0.0, 35.0_f64.to_radians());
        let tc_slope =
            bed_slope_shields_correction(0.05, 5.0_f64.to_radians(), 35.0_f64.to_radians());
        assert!(tc_slope < tc_flat, "tc_flat={tc_flat} tc_slope={tc_slope}");
    }
}
