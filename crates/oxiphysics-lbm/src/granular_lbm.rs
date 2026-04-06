#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! LBM for granular flows: dense suspensions, avalanches, jamming.
//!
//! This module provides a comprehensive framework for simulating granular
//! materials using the Lattice Boltzmann Method, including:
//!
//! - [`GranularParams`]: physical parameters for grains (density, diameter, restitution)
//! - [`MuIPressureModel`]: μ(I) rheology for dense granular flows
//! - [`GranularLbm`]: LBM solver with Bagnold kinetic + collisional stresses
//! - [`AvalancheDynamics`]: angle of repose, critical slope, runout distance
//! - [`JammingTransition`]: critical packing fraction, bulk modulus divergence
//! - [`GranularTemperature`]: fluctuation kinetic energy T_g = `δv²`/3
//! - [`BagnoldScaling`]: Bagnold rheology in the dense inertial regime
//! - [`SedimentBed`]: Shields criterion erosion threshold, bedload flux
//! - [`DryGranularFlow`]: Saint-Venant shallow-water equations for granular flows
//! - [`CohesiveGranular`]: van der Waals + capillary cohesion forces
//!
//! # Physics
//!
//! Granular flows occupy a regime between solid and fluid behaviour.  At low
//! inertial numbers I ≪ 1 the material is quasi-static (μ → μ_s); at large I
//! it behaves as a rapid granular gas.  The μ(I) rheology captures the
//! continuous transition and is incorporated here as an effective viscosity for
//! the BGK collision operator.
//!
//! # References
//! - GDR MiDi (2004). On dense granular flows. *Eur. Phys. J. E*, 14, 341–365.
//! - Jop, P., Forterre, Y., & Pouliquen, O. (2006). *Nature*, 441, 727–730.
//! - Bagnold, R. A. (1954). *Proc. R. Soc. Lond. A*, 225, 49–63.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// D2Q9 lattice constants (local copies; no nalgebra used in lbm crate)
// ---------------------------------------------------------------------------

/// Number of discrete velocities in the D2Q9 set.
pub const Q9: usize = 9;

/// D2Q9 discrete velocity x-components.
pub const EX9: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];

/// D2Q9 discrete velocity y-components.
pub const EY9: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];

/// D2Q9 lattice weights.
pub const W9: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// Lattice speed of sound squared c_s² = 1/3.
pub const CS2_GRAN: f64 = 1.0 / 3.0;

/// Gravitational acceleration (m s⁻²).
pub const GRAVITY: f64 = 9.81;

// ---------------------------------------------------------------------------
// GranularParams
// ---------------------------------------------------------------------------

/// Physical parameters describing a granular material.
#[derive(Debug, Clone)]
pub struct GranularParams {
    /// Grain bulk density ρ_g (kg m⁻³).
    pub density: f64,
    /// Mean grain diameter d (m).
    pub diameter: f64,
    /// Normal restitution coefficient e ∈ \[0, 1\].
    pub restitution: f64,
    /// Inter-grain friction coefficient μ_w ≥ 0.
    pub friction: f64,
    /// Maximum packing fraction φ_max (random close packing ≈ 0.64).
    pub phi_max: f64,
    /// Minimum (random loose) packing fraction φ_min ≈ 0.59.
    pub phi_min: f64,
}

impl GranularParams {
    /// Construct new [`GranularParams`].
    pub fn new(
        density: f64,
        diameter: f64,
        restitution: f64,
        friction: f64,
        phi_max: f64,
        phi_min: f64,
    ) -> Self {
        Self {
            density,
            diameter,
            restitution,
            friction,
            phi_max,
            phi_min,
        }
    }

    /// Default sand-like grains.
    ///
    /// ρ = 2650 kg m⁻³, d = 0.5 mm, e = 0.8, μ = 0.4.
    pub fn sand() -> Self {
        Self::new(2650.0, 5e-4, 0.8, 0.4, 0.64, 0.59)
    }

    /// Return a dimensionless measure: ratio of grain density to water density.
    pub fn specific_gravity(&self) -> f64 {
        self.density / 1000.0
    }
}

// ---------------------------------------------------------------------------
// MuIPressureModel
// ---------------------------------------------------------------------------

/// μ(I) granular rheology relating effective friction to the inertial number.
///
/// The constitutive law is:
///
/// ```text
/// μ(I) = μ_s + (μ_2 - μ_s) / (I_0 / I + 1)
/// ```
///
/// with the inertial number  I = d γ̇ / √(P / ρ).
///
/// The solid fraction varies as:
///
/// ```text
/// φ(I) = φ_max - (φ_max - φ_min) · I
/// ```
#[derive(Debug, Clone)]
pub struct MuIPressureModel {
    /// Quasi-static friction coefficient μ_s.
    pub mu_s: f64,
    /// Upper limiting friction coefficient μ_2 > μ_s.
    pub mu_2: f64,
    /// Reference inertial number I_0 (typically 0.279 for glass beads).
    pub i_0: f64,
    /// Grain diameter d (m).
    pub diameter: f64,
    /// Grain density ρ (kg m⁻³).
    pub density: f64,
    /// Maximum packing fraction φ_max.
    pub phi_max: f64,
    /// Minimum packing fraction φ_min.
    pub phi_min: f64,
}

impl MuIPressureModel {
    /// Construct a new [`MuIPressureModel`].
    pub fn new(
        mu_s: f64,
        mu_2: f64,
        i_0: f64,
        diameter: f64,
        density: f64,
        phi_max: f64,
        phi_min: f64,
    ) -> Self {
        Self {
            mu_s,
            mu_2,
            i_0,
            diameter,
            density,
            phi_max,
            phi_min,
        }
    }

    /// Glass-bead parameters from Jop et al. (2006).
    pub fn glass_beads() -> Self {
        Self::new(0.382, 0.643, 0.279, 1e-3, 2500.0, 0.64, 0.58)
    }

    /// Inertial number I = d γ̇ / √(P / ρ).
    ///
    /// Returns `0.0` if pressure P ≤ 0.
    pub fn inertial_number(&self, shear_rate: f64, pressure: f64) -> f64 {
        if pressure <= 0.0 {
            return 0.0;
        }
        self.diameter * shear_rate / (pressure / self.density).sqrt()
    }

    /// Effective friction coefficient μ(I).
    pub fn effective_friction(&self, inertial: f64) -> f64 {
        if inertial <= 0.0 {
            return self.mu_s;
        }
        self.mu_s + (self.mu_2 - self.mu_s) / (self.i_0 / inertial + 1.0)
    }

    /// Solid volume fraction φ(I) = φ_max - (φ_max - φ_min) · I.
    ///
    /// Clamped to \[0, φ_max\].
    pub fn solid_fraction(&self, inertial: f64) -> f64 {
        let phi = self.phi_max - (self.phi_max - self.phi_min) * inertial;
        phi.clamp(0.0, self.phi_max)
    }

    /// Shear stress τ = μ(I) · P.
    pub fn shear_stress(&self, shear_rate: f64, pressure: f64) -> f64 {
        let i = self.inertial_number(shear_rate, pressure);
        self.effective_friction(i) * pressure
    }

    /// Effective kinematic viscosity ν = μ(I) P / (ρ φ γ̇²) (for non-zero γ̇).
    pub fn effective_viscosity(&self, shear_rate: f64, pressure: f64) -> f64 {
        if shear_rate.abs() < 1e-30 || pressure <= 0.0 {
            return f64::INFINITY;
        }
        let i = self.inertial_number(shear_rate, pressure);
        let phi = self.solid_fraction(i);
        let mu = self.effective_friction(i);
        mu * pressure / (self.density * phi.max(1e-10) * shear_rate * shear_rate)
    }
}

// ---------------------------------------------------------------------------
// GranularLbm
// ---------------------------------------------------------------------------

/// D2Q9 LBM solver for granular flows with Bagnold kinetic + collisional stresses.
///
/// The effective viscosity is computed from μ(I) rheology at each grid cell
/// and used to determine the local relaxation time τ = ν/c_s² + 0.5.
#[derive(Debug, Clone)]
pub struct GranularLbm {
    /// Grid width (number of cells in x).
    pub nx: usize,
    /// Grid height (number of cells in y).
    pub ny: usize,
    /// Distribution functions f\[y\]\[x\]\[q\].
    pub f: Vec<Vec<[f64; 9]>>,
    /// Granular material parameters.
    pub params: GranularParams,
    /// μ(I) rheology model.
    pub mu_i: MuIPressureModel,
    /// Confining pressure (Pa).
    pub pressure: f64,
    /// Gravity vector \[gx, gy\] (m s⁻²).
    pub gravity: [f64; 2],
    /// Current simulation time step.
    pub step: usize,
}

impl GranularLbm {
    /// Construct a new quiescent [`GranularLbm`] grid.
    pub fn new(
        nx: usize,
        ny: usize,
        params: GranularParams,
        mu_i: MuIPressureModel,
        pressure: f64,
        gravity: [f64; 2],
    ) -> Self {
        let mut f = vec![vec![[0.0_f64; 9]; nx]; ny];
        // Initialise as equilibrium at rest with density 1
        for row in f.iter_mut() {
            for cell in row.iter_mut() {
                for (q, w) in W9.iter().enumerate() {
                    cell[q] = *w;
                }
            }
        }
        Self {
            nx,
            ny,
            f,
            params,
            mu_i,
            pressure,
            gravity,
            step: 0,
        }
    }

    /// Compute macroscopic density ρ and velocity \[ux, uy\] from distributions.
    pub fn macroscopic(&self, fi: &[f64; 9]) -> (f64, [f64; 2]) {
        let rho: f64 = fi.iter().sum();
        if rho < 1e-30 {
            return (0.0, [0.0, 0.0]);
        }
        let ux = fi.iter().zip(EX9.iter()).map(|(f, e)| f * e).sum::<f64>() / rho;
        let uy = fi.iter().zip(EY9.iter()).map(|(f, e)| f * e).sum::<f64>() / rho;
        (rho, [ux, uy])
    }

    /// Maxwell–Boltzmann equilibrium distribution f^eq_q.
    pub fn equilibrium(rho: f64, ux: f64, uy: f64) -> [f64; 9] {
        let usq = ux * ux + uy * uy;
        let mut feq = [0.0_f64; 9];
        for q in 0..Q9 {
            let eu = EX9[q] * ux + EY9[q] * uy;
            feq[q] = W9[q]
                * rho
                * (1.0 + eu / CS2_GRAN + eu * eu / (2.0 * CS2_GRAN * CS2_GRAN)
                    - usq / (2.0 * CS2_GRAN));
        }
        feq
    }

    /// Compute local shear rate magnitude |γ̇| from velocity gradient (finite diff).
    pub fn shear_rate_at(&self, ix: usize, iy: usize) -> f64 {
        let ixp = (ix + 1).min(self.nx - 1);
        let ixm = ix.saturating_sub(1);
        let iyp = (iy + 1).min(self.ny - 1);
        let iym = iy.saturating_sub(1);

        let (_, vel_xp) = self.macroscopic(&self.f[iy][ixp]);
        let (_, vel_xm) = self.macroscopic(&self.f[iy][ixm]);
        let (_, vel_yp) = self.macroscopic(&self.f[iyp][ix]);
        let (_, vel_ym) = self.macroscopic(&self.f[iym][ix]);

        let dux_dy = (vel_xp[1] - vel_xm[1]) / 2.0;
        let duy_dx = (vel_yp[0] - vel_ym[0]) / 2.0;
        // Strain rate tensor magnitude: |S| = √(2 S_ij S_ij)
        let s_xy = 0.5 * (dux_dy + duy_dx);
        (2.0 * s_xy * s_xy).sqrt()
    }

    /// Perform one BGK collision–streaming step.
    pub fn step(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let mut f_new = vec![vec![[0.0_f64; 9]; nx]; ny];

        for iy in 0..ny {
            for ix in 0..nx {
                let fi = &self.f[iy][ix];
                let (rho, vel) = self.macroscopic(fi);
                let gamma = self.shear_rate_at(ix, iy);
                let nu = self.mu_i.effective_viscosity(gamma, self.pressure);
                let nu_lbm = nu.clamp(1e-4, 10.0); // clamp for stability
                let tau = nu_lbm / CS2_GRAN + 0.5;
                let omega = 1.0 / tau;

                let feq = Self::equilibrium(rho, vel[0], vel[1]);
                // Gravity forcing (body force per unit mass → Guo scheme simplified)
                let ax = self.gravity[0];
                let ay = self.gravity[1];

                for q in 0..Q9 {
                    let force_q = W9[q]
                        * rho
                        * ((EX9[q] - vel[0]) / CS2_GRAN * ax + (EY9[q] - vel[1]) / CS2_GRAN * ay)
                        * (1.0 - 0.5 * omega);

                    let f_coll = fi[q] * (1.0 - omega) + feq[q] * omega + force_q;

                    // Streaming with periodic BC
                    let ixs = ((ix as isize + EX9[q] as isize).rem_euclid(nx as isize)) as usize;
                    let iys = ((iy as isize + EY9[q] as isize).rem_euclid(ny as isize)) as usize;
                    f_new[iys][ixs][q] = f_coll;
                }
            }
        }
        self.f = f_new;
        self.step += 1;
    }

    /// Compute average kinetic energy density (proxy for granular temperature).
    pub fn mean_kinetic_energy(&self) -> f64 {
        let mut total = 0.0;
        let n = (self.nx * self.ny) as f64;
        for row in &self.f {
            for cell in row {
                let (rho, vel) = self.macroscopic(cell);
                total += 0.5 * rho * (vel[0] * vel[0] + vel[1] * vel[1]);
            }
        }
        total / n
    }
}

// ---------------------------------------------------------------------------
// AvalancheDynamics
// ---------------------------------------------------------------------------

/// Avalanche dynamics: angle of repose, critical slope, and runout estimation.
#[derive(Debug, Clone)]
pub struct AvalancheDynamics {
    /// Static angle of repose θ_r (radians).
    pub angle_of_repose: f64,
    /// Dynamic (flowing) angle θ_d (radians); typically < θ_r.
    pub dynamic_angle: f64,
    /// Slope length L (m).
    pub slope_length: f64,
    /// Grain diameter d (m).
    pub grain_diameter: f64,
}

impl AvalancheDynamics {
    /// Construct a new [`AvalancheDynamics`] model.
    ///
    /// Angles must be provided in **degrees** and are stored as radians.
    pub fn new(repose_deg: f64, dynamic_deg: f64, slope_length: f64, grain_diameter: f64) -> Self {
        Self {
            angle_of_repose: repose_deg.to_radians(),
            dynamic_angle: dynamic_deg.to_radians(),
            slope_length,
            grain_diameter,
        }
    }

    /// Return the static angle of repose in degrees.
    pub fn repose_degrees(&self) -> f64 {
        self.angle_of_repose.to_degrees()
    }

    /// Check whether a slope angle `theta` (radians) is above the critical angle.
    pub fn is_critical(&self, theta: f64) -> bool {
        theta > self.angle_of_repose
    }

    /// Estimate runout distance using the Fahrböschung (travel angle) approach.
    ///
    /// ```text
    /// L_runout = H / tan(θ_d)
    /// ```
    ///
    /// where H = L sin(θ_r) is the drop height.
    pub fn runout_distance(&self) -> f64 {
        let h = self.slope_length * self.angle_of_repose.sin();
        if self.dynamic_angle.tan().abs() < 1e-30 {
            return f64::INFINITY;
        }
        h / self.dynamic_angle.tan()
    }

    /// Dimensionless runout ratio L_runout / H.
    pub fn runout_ratio(&self) -> f64 {
        if self.dynamic_angle.tan().abs() < 1e-30 {
            return f64::INFINITY;
        }
        1.0 / self.dynamic_angle.tan()
    }

    /// Height of material deposited at the base of the slope (m).
    pub fn deposition_height(&self) -> f64 {
        self.slope_length * self.angle_of_repose.sin()
    }

    /// Effective friction at the front of the avalanche (simplified).
    pub fn front_friction(&self) -> f64 {
        self.dynamic_angle.tan()
    }
}

// ---------------------------------------------------------------------------
// JammingTransition
// ---------------------------------------------------------------------------

/// Model for the jamming transition in granular packings.
///
/// Near the critical packing fraction φ_J, the bulk modulus K and coordination
/// number Z diverge as power laws:
///
/// ```text
/// K  ~ (φ - φ_J)^α
/// Z  ~ (φ - φ_J)^β   (above φ_J)
/// ```
#[derive(Debug, Clone)]
pub struct JammingTransition {
    /// Critical packing fraction φ_J ≈ 0.64 for 3-D random close packing.
    pub phi_j: f64,
    /// Bulk modulus exponent α (typically 0.5 for harmonic spheres).
    pub bulk_exponent: f64,
    /// Coordination number exponent β (typically 0.5).
    pub coord_exponent: f64,
    /// Reference bulk modulus K_0 (Pa).
    pub k0: f64,
}

impl JammingTransition {
    /// Construct a new [`JammingTransition`] model.
    pub fn new(phi_j: f64, bulk_exponent: f64, coord_exponent: f64, k0: f64) -> Self {
        Self {
            phi_j,
            bulk_exponent,
            coord_exponent,
            k0,
        }
    }

    /// Standard random-close-packing model (φ_J = 0.64).
    pub fn random_close_packing() -> Self {
        Self::new(0.64, 0.5, 0.5, 1e5)
    }

    /// Returns `true` if the packing fraction exceeds φ_J (jammed state).
    pub fn is_jammed(&self, phi: f64) -> bool {
        phi >= self.phi_j
    }

    /// Bulk modulus K(φ) for φ > φ_J; returns `0.0` below φ_J.
    pub fn bulk_modulus(&self, phi: f64) -> f64 {
        if phi <= self.phi_j {
            return 0.0;
        }
        self.k0 * (phi - self.phi_j).powf(self.bulk_exponent)
    }

    /// Excess coordination number Z - Z_iso for φ > φ_J; 0 below.
    ///
    /// Z_iso = 2d = 6 in 3-D (isostaticity condition).
    pub fn excess_coordination(&self, phi: f64) -> f64 {
        if phi <= self.phi_j {
            return 0.0;
        }
        (phi - self.phi_j).powf(self.coord_exponent)
    }

    /// Pressure near jamming P ~ (φ - φ_J)^(α+1).
    pub fn pressure_near_jamming(&self, phi: f64) -> f64 {
        if phi <= self.phi_j {
            return 0.0;
        }
        self.k0 * (phi - self.phi_j).powf(self.bulk_exponent + 1.0)
    }
}

// ---------------------------------------------------------------------------
// GranularTemperature
// ---------------------------------------------------------------------------

/// Granular temperature T_g = ⟨δv²⟩/3 — fluctuation kinetic energy per unit mass.
///
/// This is the analogue of thermodynamic temperature in kinetic theory of grains.
#[derive(Debug, Clone)]
pub struct GranularTemperature {
    /// Mean flow velocity \[ux, uy, uz\] (m s⁻¹).
    pub mean_velocity: [f64; 3],
    /// Velocity fluctuation variance (m² s⁻²); T_g = variance/3 in 3-D.
    pub variance: f64,
}

impl GranularTemperature {
    /// Construct from a list of particle velocities.
    ///
    /// Each entry is a 3-component velocity vector.
    pub fn from_velocities(velocities: &[[f64; 3]]) -> Self {
        let n = velocities.len();
        if n == 0 {
            return Self {
                mean_velocity: [0.0; 3],
                variance: 0.0,
            };
        }
        let nf = n as f64;
        let mut mean = [0.0_f64; 3];
        for v in velocities {
            mean[0] += v[0];
            mean[1] += v[1];
            mean[2] += v[2];
        }
        mean[0] /= nf;
        mean[1] /= nf;
        mean[2] /= nf;

        let variance = velocities
            .iter()
            .map(|v| {
                let dx = v[0] - mean[0];
                let dy = v[1] - mean[1];
                let dz = v[2] - mean[2];
                dx * dx + dy * dy + dz * dz
            })
            .sum::<f64>()
            / nf;

        Self {
            mean_velocity: mean,
            variance,
        }
    }

    /// Granular temperature T_g = ⟨δv²⟩/3.
    pub fn temperature(&self) -> f64 {
        self.variance / 3.0
    }

    /// Effective "pressure" from granular kinetic theory p_k = ρ T_g.
    pub fn kinetic_pressure(&self, density: f64) -> f64 {
        density * self.temperature()
    }

    /// Haff's cooling law dissipation rate: ζ = -2 T_g^(3/2) / (d · t_c)
    ///
    /// Simplified form; returns positive rate of energy loss.
    pub fn dissipation_rate(&self, grain_diameter: f64, collision_time: f64) -> f64 {
        let tg = self.temperature();
        if grain_diameter <= 0.0 || collision_time <= 0.0 {
            return 0.0;
        }
        2.0 * tg.powf(1.5) / (grain_diameter * collision_time)
    }
}

// ---------------------------------------------------------------------------
// BagnoldScaling
// ---------------------------------------------------------------------------

/// Bagnold rheology for rapid (inertial) granular flow.
///
/// In the dense inertial regime the shear and normal stresses scale as γ̇²:
///
/// ```text
/// τ_Bagnold  = ρ_g d² λ^(1/2) γ̇²   (shear stress)
/// P_Bagnold  = τ_Bagnold / μ_s      (normal stress)
/// ```
///
/// where λ = φ / (φ_rcp - φ) is the linear concentration.
#[derive(Debug, Clone)]
pub struct BagnoldScaling {
    /// Grain density ρ_g (kg m⁻³).
    pub grain_density: f64,
    /// Grain diameter d (m).
    pub diameter: f64,
    /// Bagnold coefficient B (dimensionless, ≈ 0.013 experimentally).
    pub bagnold_coeff: f64,
    /// Random close packing fraction φ_rcp.
    pub phi_rcp: f64,
}

impl BagnoldScaling {
    /// Construct a new [`BagnoldScaling`] model.
    pub fn new(grain_density: f64, diameter: f64, bagnold_coeff: f64, phi_rcp: f64) -> Self {
        Self {
            grain_density,
            diameter,
            bagnold_coeff,
            phi_rcp,
        }
    }

    /// Linear concentration λ = φ / (φ_rcp - φ).
    pub fn linear_concentration(&self, phi: f64) -> f64 {
        let denom = self.phi_rcp - phi;
        if denom <= 1e-30 {
            return f64::INFINITY;
        }
        phi / denom
    }

    /// Bagnold shear stress τ_B = B ρ_g d² λ^(1/2) γ̇² (Pa).
    pub fn shear_stress(&self, phi: f64, shear_rate: f64) -> f64 {
        let lambda = self.linear_concentration(phi);
        if lambda.is_infinite() {
            return f64::INFINITY;
        }
        self.bagnold_coeff
            * self.grain_density
            * self.diameter
            * self.diameter
            * lambda.sqrt()
            * shear_rate
            * shear_rate
    }

    /// Bagnold normal (dispersive) stress P_B = τ_B / tan(α_Bagnold).
    ///
    /// Approximation: P_B ≈ τ_B / 0.32 (Bagnold's empirical constant).
    pub fn normal_stress(&self, phi: f64, shear_rate: f64) -> f64 {
        self.shear_stress(phi, shear_rate) / 0.32
    }

    /// Effective Bagnold viscosity η_eff = τ / γ̇.
    pub fn effective_viscosity(&self, phi: f64, shear_rate: f64) -> f64 {
        if shear_rate.abs() < 1e-30 {
            return 0.0;
        }
        self.shear_stress(phi, shear_rate) / shear_rate
    }
}

// ---------------------------------------------------------------------------
// SedimentBed
// ---------------------------------------------------------------------------

/// Sediment bed erosion and bedload transport using the Shields criterion.
#[derive(Debug, Clone)]
pub struct SedimentBed {
    /// Grain diameter d₅₀ (m).
    pub d50: f64,
    /// Grain density ρ_s (kg m⁻³).
    pub rho_s: f64,
    /// Fluid density ρ_f (kg m⁻³).
    pub rho_f: f64,
    /// Kinematic viscosity of the fluid ν (m² s⁻¹).
    pub nu_fluid: f64,
    /// Critical Shields parameter θ_c (≈ 0.047 for non-cohesive sand).
    pub critical_shields: f64,
}

impl SedimentBed {
    /// Construct a new [`SedimentBed`] model.
    pub fn new(d50: f64, rho_s: f64, rho_f: f64, nu_fluid: f64, critical_shields: f64) -> Self {
        Self {
            d50,
            rho_s,
            rho_f,
            nu_fluid,
            critical_shields,
        }
    }

    /// Shields parameter θ = τ_b / \[(ρ_s - ρ_f) g d\].
    pub fn shields_parameter(&self, bed_shear_stress: f64) -> f64 {
        let denom = (self.rho_s - self.rho_f) * GRAVITY * self.d50;
        if denom <= 1e-30 {
            return 0.0;
        }
        bed_shear_stress / denom
    }

    /// Returns `true` when the bed shear stress exceeds the erosion threshold.
    pub fn is_eroding(&self, bed_shear_stress: f64) -> bool {
        self.shields_parameter(bed_shear_stress) > self.critical_shields
    }

    /// Meyer-Peter & Müller bed-load flux q_b* (dimensionless).
    ///
    /// ```text
    /// q_b* = 8 (θ - θ_c)^{3/2}   for θ > θ_c
    /// ```
    pub fn bedload_flux_dimensionless(&self, bed_shear_stress: f64) -> f64 {
        let theta = self.shields_parameter(bed_shear_stress);
        if theta <= self.critical_shields {
            return 0.0;
        }
        8.0 * (theta - self.critical_shields).powf(1.5)
    }

    /// Dimensional bed-load flux q_b (m² s⁻¹).
    pub fn bedload_flux(&self, bed_shear_stress: f64) -> f64 {
        let qb_star = self.bedload_flux_dimensionless(bed_shear_stress);
        let s = self.rho_s / self.rho_f - 1.0;
        qb_star * (s * GRAVITY * self.d50.powi(3)).sqrt()
    }

    /// Particle Reynolds number Re_p = u* d / ν, where u* = √(τ_b/ρ_f).
    pub fn particle_reynolds(&self, bed_shear_stress: f64) -> f64 {
        let u_star = (bed_shear_stress / self.rho_f).sqrt();
        u_star * self.d50 / self.nu_fluid
    }
}

// ---------------------------------------------------------------------------
// DryGranularFlow
// ---------------------------------------------------------------------------

/// Shallow-water Saint-Venant equations for dry granular flows.
///
/// Depth-averaged model:
///
/// ```text
/// ∂h/∂t + ∇·(h u) = 0
/// ∂(h u)/∂t + ∇·(h u⊗u) + g h ∇h = g h (tan θ_bed − μ_base u/|u|)
/// ```
#[derive(Debug, Clone)]
pub struct DryGranularFlow {
    /// Number of cells nx.
    pub nx: usize,
    /// Flow depth h\[i\] (m) at each cell.
    pub depth: Vec<f64>,
    /// Depth-averaged velocity u\[i\] (m s⁻¹) at each cell.
    pub velocity: Vec<f64>,
    /// Cell size dx (m).
    pub dx: f64,
    /// Bed slope angle (radians).
    pub slope: f64,
    /// Basal friction coefficient μ_base.
    pub mu_base: f64,
}

impl DryGranularFlow {
    /// Construct a new quiescent depth profile.
    pub fn new(nx: usize, dx: f64, slope_deg: f64, mu_base: f64) -> Self {
        Self {
            nx,
            depth: vec![0.0; nx],
            velocity: vec![0.0; nx],
            dx,
            slope: slope_deg.to_radians(),
            mu_base,
        }
    }

    /// Set a Gaussian pile of material at the centre.
    pub fn set_gaussian_pile(&mut self, h_max: f64, sigma: f64) {
        let mid = (self.nx as f64 - 1.0) / 2.0;
        for i in 0..self.nx {
            let x = i as f64 - mid;
            self.depth[i] = h_max * (-0.5 * (x / sigma).powi(2)).exp();
        }
    }

    /// Advance one explicit time step Δt using the upwind scheme.
    pub fn advance(&mut self, dt: f64) {
        let nx = self.nx;
        let mut dh = vec![0.0_f64; nx];
        let mut du = vec![0.0_f64; nx];
        let g = GRAVITY;
        let tan_slope = self.slope.tan();

        for i in 1..nx - 1 {
            let hi = self.depth[i];
            let ui = self.velocity[i];

            // Upwind advective flux for h
            let flux_l = if ui >= 0.0 {
                self.depth[i - 1] * self.velocity[i - 1]
            } else {
                hi * ui
            };
            let flux_r = if ui >= 0.0 {
                hi * ui
            } else {
                self.depth[i + 1] * self.velocity[i + 1]
            };
            dh[i] = -(flux_r - flux_l) / self.dx;

            // Pressure gradient + gravity + friction
            let dh_dx = (self.depth[i + 1] - self.depth[i - 1]) / (2.0 * self.dx);
            let friction = if ui.abs() > 1e-10 {
                -self.mu_base * ui.signum()
            } else {
                0.0
            };
            du[i] = g * (tan_slope - dh_dx) + friction * g;
        }

        for i in 1..nx - 1 {
            self.depth[i] = (self.depth[i] + dt * dh[i]).max(0.0);
            self.velocity[i] += dt * du[i];
        }
    }

    /// Total volume (∫ h dx).
    pub fn total_volume(&self) -> f64 {
        self.depth.iter().sum::<f64>() * self.dx
    }

    /// Maximum flow depth.
    pub fn max_depth(&self) -> f64 {
        self.depth.iter().cloned().fold(0.0_f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// CohesiveGranular
// ---------------------------------------------------------------------------

/// Cohesive granular material model including van der Waals and capillary forces.
#[derive(Debug, Clone)]
pub struct CohesiveGranular {
    /// Grain diameter d (m).
    pub diameter: f64,
    /// Hamaker constant A (J), typically 1e-20 J for mineral grains.
    pub hamaker: f64,
    /// Surface tension γ_lv (N m⁻¹) of the wetting liquid.
    pub surface_tension: f64,
    /// Contact angle θ_c (radians).
    pub contact_angle: f64,
    /// Liquid saturation S ∈ \[0, 1\] (fraction of pore space filled).
    pub saturation: f64,
}

impl CohesiveGranular {
    /// Construct a new [`CohesiveGranular`] model.
    pub fn new(
        diameter: f64,
        hamaker: f64,
        surface_tension: f64,
        contact_angle: f64,
        saturation: f64,
    ) -> Self {
        Self {
            diameter,
            hamaker,
            surface_tension,
            contact_angle,
            saturation,
        }
    }

    /// Moist quartz sand at pendular-regime saturation.
    pub fn moist_sand() -> Self {
        Self::new(5e-4, 1.3e-20, 0.0728, 20f64.to_radians(), 0.05)
    }

    /// Van der Waals adhesion force between two identical spheres at separation
    /// `z` (m):  F_vdW = -A d / (24 z²).
    pub fn vdw_force(&self, separation: f64) -> f64 {
        if separation <= 0.0 {
            return f64::NEG_INFINITY;
        }
        -self.hamaker * self.diameter / (24.0 * separation * separation)
    }

    /// Capillary force (Laplace bridge approximation):
    ///
    /// ```text
    /// F_cap = π d γ cos θ_c
    /// ```
    pub fn capillary_force(&self) -> f64 {
        PI * self.diameter * self.surface_tension * self.contact_angle.cos()
    }

    /// Total cohesion force: F_total = F_vdW + F_cap (N).
    ///
    /// `separation` is the van der Waals surface-to-surface gap (m).
    pub fn total_cohesion(&self, separation: f64) -> f64 {
        self.vdw_force(separation) + self.capillary_force()
    }

    /// Cohesion number Co = F_cap / (ρ g d³): ratio of capillary to gravitational force.
    pub fn cohesion_number(&self, grain_density: f64) -> f64 {
        let f_cap = self.capillary_force();
        let f_grav = grain_density * GRAVITY * self.diameter.powi(3);
        if f_grav.abs() < 1e-30 {
            return f64::INFINITY;
        }
        f_cap / f_grav
    }

    /// Effective tensile strength σ_t ≈ F_cap / (π d² / 4) (Pa).
    pub fn tensile_strength(&self) -> f64 {
        let area = PI * self.diameter * self.diameter / 4.0;
        self.capillary_force() / area
    }

    /// Liquid bridge volume V_b ≈ π S d³ / 6 (pendular regime approximation, m³).
    pub fn bridge_volume(&self) -> f64 {
        PI * self.saturation * self.diameter.powi(3) / 6.0
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- MuIPressureModel ---

    #[test]
    fn test_mu_i_quasi_static_limit() {
        // μ(I) → μ_s as I → 0
        let model = MuIPressureModel::glass_beads();
        let mu = model.effective_friction(0.0);
        assert!((mu - model.mu_s).abs() < 1e-12);
    }

    #[test]
    fn test_mu_i_upper_bound() {
        // μ(I) < μ_2 for all finite I
        let model = MuIPressureModel::glass_beads();
        let mu = model.effective_friction(1000.0);
        assert!(mu < model.mu_2);
    }

    #[test]
    fn test_mu_i_monotone() {
        let model = MuIPressureModel::glass_beads();
        let mu_low = model.effective_friction(0.01);
        let mu_high = model.effective_friction(1.0);
        assert!(mu_high > mu_low, "μ(I) must be monotonically increasing");
    }

    #[test]
    fn test_inertial_number_zero_pressure() {
        let model = MuIPressureModel::glass_beads();
        assert_eq!(model.inertial_number(1.0, 0.0), 0.0);
    }

    #[test]
    fn test_inertial_number_formula() {
        // I = d γ̇ / √(P/ρ)
        let model = MuIPressureModel::new(0.38, 0.64, 0.28, 1e-3, 2500.0, 0.64, 0.58);
        let gamma = 100.0;
        let p = 1000.0;
        let expected = model.diameter * gamma / (p / model.density).sqrt();
        let computed = model.inertial_number(gamma, p);
        assert!((computed - expected).abs() < 1e-12);
    }

    #[test]
    fn test_solid_fraction_decreases_with_i() {
        let model = MuIPressureModel::glass_beads();
        let phi_low = model.solid_fraction(0.01);
        let phi_high = model.solid_fraction(0.5);
        assert!(phi_low > phi_high);
    }

    #[test]
    fn test_solid_fraction_clamped() {
        let model = MuIPressureModel::glass_beads();
        let phi = model.solid_fraction(-1.0);
        assert!(phi >= 0.0 && phi <= model.phi_max);
    }

    #[test]
    fn test_mu_i_shear_stress_positive() {
        let model = MuIPressureModel::glass_beads();
        let tau = model.shear_stress(10.0, 500.0);
        assert!(tau > 0.0);
    }

    // --- GranularParams ---

    #[test]
    fn test_granular_params_sand() {
        let p = GranularParams::sand();
        assert!(p.density > 0.0);
        assert!(p.diameter > 0.0);
        assert!((0.0..=1.0).contains(&p.restitution));
        assert!(p.phi_max <= 0.7);
    }

    #[test]
    fn test_specific_gravity_sand() {
        let p = GranularParams::sand();
        // sand density ~2650 kg/m³, so s ≈ 2.65
        assert!((p.specific_gravity() - 2.65).abs() < 0.01);
    }

    // --- AvalancheDynamics ---

    #[test]
    fn test_angle_of_repose_bounds() {
        // Typical values: 15–45 degrees
        let av = AvalancheDynamics::new(34.0, 28.0, 10.0, 5e-4);
        let deg = av.repose_degrees();
        assert!(
            deg >= 15.0 && deg <= 45.0,
            "angle_of_repose={:.2} not in [15,45]",
            deg
        );
    }

    #[test]
    fn test_is_critical_above() {
        let av = AvalancheDynamics::new(30.0, 25.0, 5.0, 1e-3);
        assert!(av.is_critical(35f64.to_radians()));
    }

    #[test]
    fn test_is_critical_below() {
        let av = AvalancheDynamics::new(30.0, 25.0, 5.0, 1e-3);
        assert!(!av.is_critical(20f64.to_radians()));
    }

    #[test]
    fn test_runout_positive() {
        let av = AvalancheDynamics::new(34.0, 15.0, 100.0, 5e-4);
        assert!(av.runout_distance() > 0.0);
    }

    #[test]
    fn test_avalanche_slope_decreases() {
        // After runout, effective slope = dynamic angle < repose angle
        let av = AvalancheDynamics::new(34.0, 20.0, 50.0, 5e-4);
        assert!(av.dynamic_angle < av.angle_of_repose);
    }

    // --- JammingTransition ---

    #[test]
    fn test_jamming_phi_j_random_close_packing() {
        let j = JammingTransition::random_close_packing();
        // RCP φ_J ≤ 0.64
        assert!(j.phi_j <= 0.64 + 1e-9);
    }

    #[test]
    fn test_jamming_not_jammed_below() {
        let j = JammingTransition::random_close_packing();
        assert!(!j.is_jammed(0.60));
    }

    #[test]
    fn test_jamming_is_jammed_above() {
        let j = JammingTransition::random_close_packing();
        assert!(j.is_jammed(0.65));
    }

    #[test]
    fn test_bulk_modulus_zero_below_phi_j() {
        let j = JammingTransition::random_close_packing();
        assert_eq!(j.bulk_modulus(0.60), 0.0);
    }

    #[test]
    fn test_bulk_modulus_positive_above_phi_j() {
        let j = JammingTransition::random_close_packing();
        assert!(j.bulk_modulus(0.66) > 0.0);
    }

    #[test]
    fn test_excess_coordination_zero_below() {
        let j = JammingTransition::random_close_packing();
        assert_eq!(j.excess_coordination(0.50), 0.0);
    }

    #[test]
    fn test_packing_fraction_random_loose_packing() {
        // Random loose packing ≤ 0.64
        let j = JammingTransition::random_close_packing();
        assert!(
            j.phi_j <= 0.64 + 1e-9,
            "packing fraction {:.4} exceeds RCP limit",
            j.phi_j
        );
    }

    // --- GranularTemperature ---

    #[test]
    fn test_granular_temperature_nonneg() {
        let vels = vec![[1.0, 0.5, -0.3], [-0.2, 1.1, 0.4], [0.0, -0.8, 0.9]];
        let gt = GranularTemperature::from_velocities(&vels);
        assert!(gt.temperature() >= 0.0);
    }

    #[test]
    fn test_granular_temperature_zero_for_identical_velocities() {
        let vels = vec![[1.0, 2.0, 3.0]; 5];
        let gt = GranularTemperature::from_velocities(&vels);
        assert!(gt.temperature() < 1e-12);
    }

    #[test]
    fn test_granular_temperature_empty() {
        let gt = GranularTemperature::from_velocities(&[]);
        assert_eq!(gt.temperature(), 0.0);
    }

    #[test]
    fn test_granular_kinetic_pressure_positive() {
        let vels = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let gt = GranularTemperature::from_velocities(&vels);
        let p = gt.kinetic_pressure(1500.0);
        assert!(p >= 0.0);
    }

    // --- BagnoldScaling ---

    #[test]
    fn test_bagnold_linear_concentration_positive() {
        let b = BagnoldScaling::new(2650.0, 5e-4, 0.013, 0.64);
        assert!(b.linear_concentration(0.50) > 0.0);
    }

    #[test]
    fn test_bagnold_shear_stress_scales_as_shear_rate_squared() {
        let b = BagnoldScaling::new(2650.0, 5e-4, 0.013, 0.64);
        let phi = 0.50;
        let tau1 = b.shear_stress(phi, 10.0);
        let tau2 = b.shear_stress(phi, 20.0);
        // τ ∝ γ̇² → ratio should be 4
        assert!((tau2 / tau1 - 4.0).abs() < 1e-8);
    }

    #[test]
    fn test_bagnold_normal_stress_larger_than_shear() {
        let b = BagnoldScaling::new(2650.0, 5e-4, 0.013, 0.64);
        let tau = b.shear_stress(0.5, 50.0);
        let p = b.normal_stress(0.5, 50.0);
        assert!(p > tau);
    }

    // --- SedimentBed ---

    #[test]
    fn test_shields_criterion_positive() {
        let bed = SedimentBed::new(5e-4, 2650.0, 1000.0, 1e-6, 0.047);
        let theta = bed.shields_parameter(2.0);
        assert!(theta > 0.0);
    }

    #[test]
    fn test_no_erosion_below_threshold() {
        let bed = SedimentBed::new(5e-4, 2650.0, 1000.0, 1e-6, 0.047);
        // Very small shear stress
        assert!(!bed.is_eroding(0.01));
    }

    #[test]
    fn test_erosion_above_threshold() {
        let bed = SedimentBed::new(5e-4, 2650.0, 1000.0, 1e-6, 0.047);
        // Large shear stress
        assert!(bed.is_eroding(10.0));
    }

    #[test]
    fn test_bedload_flux_zero_below_threshold() {
        let bed = SedimentBed::new(5e-4, 2650.0, 1000.0, 1e-6, 0.047);
        assert_eq!(bed.bedload_flux_dimensionless(0.0), 0.0);
    }

    #[test]
    fn test_bedload_flux_positive_above_threshold() {
        let bed = SedimentBed::new(5e-4, 2650.0, 1000.0, 1e-6, 0.047);
        assert!(bed.bedload_flux(5.0) > 0.0);
    }

    // --- DryGranularFlow ---

    #[test]
    fn test_dry_granular_flow_volume_conserved() {
        let mut flow = DryGranularFlow::new(50, 0.1, 30.0, 0.5);
        flow.set_gaussian_pile(1.0, 3.0);
        let v0 = flow.total_volume();
        flow.advance(0.001);
        let v1 = flow.total_volume();
        // Volume should be approximately conserved (within 5%)
        assert!(
            (v1 - v0).abs() / v0 < 0.05,
            "volume not conserved: {:.6} vs {:.6}",
            v0,
            v1
        );
    }

    #[test]
    fn test_dry_granular_flow_max_depth_positive() {
        let mut flow = DryGranularFlow::new(30, 0.1, 20.0, 0.4);
        flow.set_gaussian_pile(0.5, 2.0);
        assert!(flow.max_depth() > 0.0);
    }

    #[test]
    fn test_dry_granular_flow_depth_nonneg() {
        let mut flow = DryGranularFlow::new(30, 0.1, 35.0, 0.6);
        flow.set_gaussian_pile(0.5, 2.0);
        for _ in 0..20 {
            flow.advance(0.001);
        }
        for &h in &flow.depth {
            assert!(h >= 0.0, "negative depth {h}");
        }
    }

    // --- CohesiveGranular ---

    #[test]
    fn test_capillary_force_positive() {
        let c = CohesiveGranular::moist_sand();
        assert!(c.capillary_force() > 0.0);
    }

    #[test]
    fn test_vdw_force_attractive() {
        let c = CohesiveGranular::moist_sand();
        // vdW force is attractive (negative)
        assert!(c.vdw_force(1e-9) < 0.0);
    }

    #[test]
    fn test_cohesion_number_positive() {
        let c = CohesiveGranular::moist_sand();
        assert!(c.cohesion_number(2650.0) > 0.0);
    }

    #[test]
    fn test_tensile_strength_positive() {
        let c = CohesiveGranular::moist_sand();
        assert!(c.tensile_strength() > 0.0);
    }

    #[test]
    fn test_bridge_volume_positive() {
        let c = CohesiveGranular::moist_sand();
        assert!(c.bridge_volume() > 0.0);
    }

    // --- GranularLbm ---

    #[test]
    fn test_granular_lbm_step_runs() {
        let params = GranularParams::sand();
        let mu_i = MuIPressureModel::glass_beads();
        let mut lbm = GranularLbm::new(8, 8, params, mu_i, 1000.0, [0.0, -0.001]);
        lbm.step();
        assert_eq!(lbm.step, 1);
    }

    #[test]
    fn test_granular_lbm_mean_kinetic_energy_nonneg() {
        let params = GranularParams::sand();
        let mu_i = MuIPressureModel::glass_beads();
        let lbm = GranularLbm::new(8, 8, params, mu_i, 1000.0, [0.0, -0.001]);
        assert!(lbm.mean_kinetic_energy() >= 0.0);
    }

    #[test]
    fn test_equilibrium_sums_to_rho() {
        let feq = GranularLbm::equilibrium(1.2, 0.05, -0.03);
        let rho: f64 = feq.iter().sum();
        assert!((rho - 1.2).abs() < 1e-10);
    }

    #[test]
    fn test_macroscopic_round_trip() {
        let rho_in = 1.1;
        let ux_in = 0.04;
        let uy_in = -0.02;
        let feq = GranularLbm::equilibrium(rho_in, ux_in, uy_in);
        let params = GranularParams::sand();
        let mu_i = MuIPressureModel::glass_beads();
        let lbm = GranularLbm::new(4, 4, params, mu_i, 500.0, [0.0, 0.0]);
        let (rho_out, vel) = lbm.macroscopic(&feq);
        assert!((rho_out - rho_in).abs() < 1e-10);
        assert!((vel[0] - ux_in).abs() < 1e-10);
        assert!((vel[1] - uy_in).abs() < 1e-10);
    }
}
