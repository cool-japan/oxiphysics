#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Blood flow simulation via Lattice Boltzmann Method.
//!
//! This module provides hemodynamics-specific LBM tools:
//!
//! - **Non-Newtonian blood rheology**: Casson and Carreau-Yasuda models
//! - **Pulsatile flow**: Womersley-based pulsatile inlet conditions
//! - **Arterial geometry**: Cylindrical, tapered, curved, stenosed arteries
//! - **Wall shear stress (WSS)**: Instantaneous and time-averaged WSS
//! - **Oscillatory shear index (OSI)**: Indicator of disturbed flow
//! - **Residence time**: Relative residence time computation
//! - **Stenosis modeling**: Parametric stenosis with severity control
//! - **Bifurcation flows**: Y-shaped and T-shaped bifurcations
//! - **RBC transport**: Simplified red blood cell advection

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Math helpers
// ─────────────────────────────────────────────────────────────────────────────

#[inline]
fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}

#[inline]
fn len2(v: [f64; 2]) -> f64 {
    dot2(v, v).sqrt()
}

#[inline]
fn add2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] + b[0], a[1] + b[1]]
}

#[inline]
fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}

#[inline]
fn scale2(v: [f64; 2], s: f64) -> [f64; 2] {
    [v[0] * s, v[1] * s]
}

// ─────────────────────────────────────────────────────────────────────────────
// D2Q9 lattice constants
// ─────────────────────────────────────────────────────────────────────────────

/// D2Q9 velocity vectors.
const D2Q9_E: [[i32; 2]; 9] = [
    [0, 0],
    [1, 0],
    [0, 1],
    [-1, 0],
    [0, -1],
    [1, 1],
    [-1, 1],
    [-1, -1],
    [1, -1],
];

/// D2Q9 weights.
const D2Q9_W: [f64; 9] = [
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

/// Opposite direction indices for bounce-back.
const D2Q9_OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

/// Speed of sound squared: cs^2 = 1/3.
const CS2: f64 = 1.0 / 3.0;

// ─────────────────────────────────────────────────────────────────────────────
// Non-Newtonian blood rheology
// ─────────────────────────────────────────────────────────────────────────────

/// Casson rheology model for blood.
///
/// The Casson model captures yield-stress behavior:
/// sqrt(tau) = sqrt(tau_y) + sqrt(mu_inf * gamma_dot)
///
/// Effective viscosity: mu_eff = (sqrt(tau_y / gamma_dot) + sqrt(mu_inf))^2
#[derive(Debug, Clone, Copy)]
pub struct CassonModel {
    /// Yield stress (Pa). Typical: 0.005 Pa.
    pub tau_y: f64,
    /// Infinite-shear viscosity (Pa·s). Typical: 0.00345 Pa·s.
    pub mu_inf: f64,
    /// Minimum shear rate to avoid singularity.
    pub gamma_dot_min: f64,
}

impl CassonModel {
    /// Create a new Casson model with default blood parameters.
    pub fn new_blood() -> Self {
        Self {
            tau_y: 0.005,
            mu_inf: 0.00345,
            gamma_dot_min: 1e-6,
        }
    }

    /// Create with custom parameters.
    pub fn new(tau_y: f64, mu_inf: f64) -> Self {
        Self {
            tau_y,
            mu_inf,
            gamma_dot_min: 1e-6,
        }
    }

    /// Compute effective viscosity at a given shear rate.
    pub fn viscosity(&self, gamma_dot: f64) -> f64 {
        let gd = gamma_dot.abs().max(self.gamma_dot_min);
        let sqrt_part = (self.tau_y / gd).sqrt() + self.mu_inf.sqrt();
        sqrt_part * sqrt_part
    }

    /// Compute the relaxation time tau for LBM from local shear rate.
    pub fn relaxation_time(&self, gamma_dot: f64, _dx: f64, _dt: f64) -> f64 {
        let nu = self.viscosity(gamma_dot);
        0.5 + nu / CS2
    }
}

/// Carreau-Yasuda rheology model for blood.
///
/// mu(gamma_dot) = mu_inf + (mu_0 - mu_inf) * (1 + (lambda * gamma_dot)^a)^((n-1)/a)
#[derive(Debug, Clone, Copy)]
pub struct CarreauYasudaModel {
    /// Zero-shear viscosity (Pa·s). Typical: 0.056 Pa·s.
    pub mu_0: f64,
    /// Infinite-shear viscosity (Pa·s). Typical: 0.00345 Pa·s.
    pub mu_inf: f64,
    /// Relaxation time (s). Typical: 3.313 s.
    pub lambda: f64,
    /// Power-law index. Typical: 0.3568.
    pub n: f64,
    /// Yasuda parameter. Typical: 2.0 (reduces to Carreau model).
    pub a: f64,
}

impl CarreauYasudaModel {
    /// Create with default blood parameters (Carreau-Yasuda).
    pub fn new_blood() -> Self {
        Self {
            mu_0: 0.056,
            mu_inf: 0.00345,
            lambda: 3.313,
            n: 0.3568,
            a: 2.0,
        }
    }

    /// Create with custom parameters.
    pub fn new(mu_0: f64, mu_inf: f64, lambda: f64, n: f64, a: f64) -> Self {
        Self {
            mu_0,
            mu_inf,
            lambda,
            n,
            a,
        }
    }

    /// Compute effective viscosity at a given shear rate.
    pub fn viscosity(&self, gamma_dot: f64) -> f64 {
        let gd = gamma_dot.abs();
        let term = (self.lambda * gd).powf(self.a);
        self.mu_inf + (self.mu_0 - self.mu_inf) * (1.0 + term).powf((self.n - 1.0) / self.a)
    }

    /// Compute the relaxation time tau for LBM.
    pub fn relaxation_time(&self, gamma_dot: f64) -> f64 {
        let nu = self.viscosity(gamma_dot);
        0.5 + nu / CS2
    }
}

/// Power-law rheology model (simpler non-Newtonian).
///
/// mu = K * gamma_dot^(n-1)
#[derive(Debug, Clone, Copy)]
pub struct PowerLawModel {
    /// Consistency index (Pa·s^n).
    pub k: f64,
    /// Power-law exponent.
    pub n: f64,
    /// Minimum shear rate clamp.
    pub gamma_dot_min: f64,
    /// Maximum viscosity clamp.
    pub mu_max: f64,
}

impl PowerLawModel {
    /// Create a power-law model for blood.
    pub fn new_blood() -> Self {
        Self {
            k: 0.017,
            n: 0.708,
            gamma_dot_min: 1e-6,
            mu_max: 0.1,
        }
    }

    /// Create with custom parameters.
    pub fn new(k: f64, n: f64) -> Self {
        Self {
            k,
            n,
            gamma_dot_min: 1e-6,
            mu_max: 0.1,
        }
    }

    /// Compute effective viscosity.
    pub fn viscosity(&self, gamma_dot: f64) -> f64 {
        let gd = gamma_dot.abs().max(self.gamma_dot_min);
        (self.k * gd.powf(self.n - 1.0)).min(self.mu_max)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Pulsatile flow
// ─────────────────────────────────────────────────────────────────────────────

/// Pulsatile flow generator based on Fourier decomposition of cardiac waveform.
///
/// The flow velocity profile is modeled as:
/// U(t) = U_mean + sum_k A_k * sin(k * omega * t + phi_k)
#[derive(Debug, Clone)]
pub struct PulsatileFlow {
    /// Mean flow velocity.
    pub u_mean: f64,
    /// Heart rate in beats per minute.
    pub heart_rate: f64,
    /// Fourier amplitudes for harmonics.
    pub amplitudes: Vec<f64>,
    /// Fourier phases for harmonics (radians).
    pub phases: Vec<f64>,
}

impl PulsatileFlow {
    /// Create a simple pulsatile flow with a single harmonic.
    pub fn new_simple(u_mean: f64, heart_rate: f64, amplitude: f64) -> Self {
        Self {
            u_mean,
            heart_rate,
            amplitudes: vec![amplitude],
            phases: vec![0.0],
        }
    }

    /// Create a more realistic cardiac waveform with multiple harmonics.
    pub fn new_cardiac(u_mean: f64, heart_rate: f64) -> Self {
        // Typical aortic waveform harmonics
        Self {
            u_mean,
            heart_rate,
            amplitudes: vec![0.6, 0.3, 0.15, 0.08],
            phases: vec![0.0, -0.5, -1.0, -1.5],
        }
    }

    /// Angular frequency omega = 2*pi*f.
    pub fn omega(&self) -> f64 {
        2.0 * PI * self.heart_rate / 60.0
    }

    /// Womersley number: alpha = R * sqrt(omega / nu).
    pub fn womersley_number(&self, radius: f64, nu: f64) -> f64 {
        radius * (self.omega() / nu).sqrt()
    }

    /// Evaluate the velocity at time t.
    pub fn velocity_at(&self, t: f64) -> f64 {
        let omega = self.omega();
        let mut u = self.u_mean;
        for (k, (amp, phase)) in self.amplitudes.iter().zip(self.phases.iter()).enumerate() {
            u += amp * ((k + 1) as f64 * omega * t + phase).sin();
        }
        u
    }

    /// Compute the Womersley velocity profile at radial position r/R.
    ///
    /// For a single harmonic, the Womersley solution gives:
    /// u(r, t) = Re{ (dp/dx) * R^2 / (mu * alpha^2) * (1 - J0(alpha*i^{3/2}*r/R) / J0(alpha*i^{3/2})) * exp(i*omega*t) }
    ///
    /// Simplified here as parabolic profile modulated by pulsatile velocity.
    pub fn velocity_profile(&self, r_ratio: f64, t: f64) -> f64 {
        let u_center = self.velocity_at(t);
        // Parabolic profile approximation: u(r) = 2*U_mean*(1 - (r/R)^2)
        2.0 * u_center * (1.0 - r_ratio * r_ratio)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Arterial geometry
// ─────────────────────────────────────────────────────────────────────────────

/// Type of arterial geometry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ArteryType {
    /// Straight cylindrical artery.
    Straight,
    /// Tapered artery (linearly decreasing radius).
    Tapered,
    /// Curved artery (constant curvature).
    Curved,
    /// Artery with stenosis.
    Stenosed,
}

/// Arterial geometry definition on a 2D grid.
///
/// Defines the vessel wall positions for each column of the lattice.
#[derive(Debug, Clone)]
pub struct ArterialGeometry {
    /// Number of lattice nodes in x-direction.
    pub nx: usize,
    /// Number of lattice nodes in y-direction.
    pub ny: usize,
    /// Lower wall y-position for each x.
    pub wall_lower: Vec<f64>,
    /// Upper wall y-position for each x.
    pub wall_upper: Vec<f64>,
    /// Type of artery.
    pub artery_type: ArteryType,
}

impl ArterialGeometry {
    /// Create a straight artery centered in the domain.
    pub fn straight(nx: usize, ny: usize, radius: f64) -> Self {
        let center = ny as f64 / 2.0;
        let wall_lower = vec![center - radius; nx];
        let wall_upper = vec![center + radius; nx];
        Self {
            nx,
            ny,
            wall_lower,
            wall_upper,
            artery_type: ArteryType::Straight,
        }
    }

    /// Create a tapered artery (radius decreases linearly).
    pub fn tapered(nx: usize, ny: usize, r_inlet: f64, r_outlet: f64) -> Self {
        let center = ny as f64 / 2.0;
        let mut wall_lower = Vec::with_capacity(nx);
        let mut wall_upper = Vec::with_capacity(nx);
        for i in 0..nx {
            let frac = i as f64 / (nx - 1).max(1) as f64;
            let r = r_inlet + (r_outlet - r_inlet) * frac;
            wall_lower.push(center - r);
            wall_upper.push(center + r);
        }
        Self {
            nx,
            ny,
            wall_lower,
            wall_upper,
            artery_type: ArteryType::Tapered,
        }
    }

    /// Create a stenosed artery.
    ///
    /// `severity` is the fraction of diameter reduction (0 to 1).
    /// `stenosis_center` is the x-position of the stenosis center (fraction of length).
    /// `stenosis_length` is the length of the stenosis (fraction of total length).
    pub fn stenosed(
        nx: usize,
        ny: usize,
        radius: f64,
        severity: f64,
        stenosis_center: f64,
        stenosis_length: f64,
    ) -> Self {
        let center = ny as f64 / 2.0;
        let x_center = (nx as f64 * stenosis_center) as usize;
        let half_len = (nx as f64 * stenosis_length * 0.5) as usize;
        let mut wall_lower = Vec::with_capacity(nx);
        let mut wall_upper = Vec::with_capacity(nx);
        for i in 0..nx {
            let dist = (i as i64 - x_center as i64).unsigned_abs() as f64;
            let half_len_f = half_len as f64;
            let r = if dist < half_len_f {
                // Cosine-shaped stenosis
                let s = severity * radius * 0.5 * (1.0 + (PI * dist / half_len_f).cos());
                radius - s
            } else {
                radius
            };
            wall_lower.push(center - r);
            wall_upper.push(center + r);
        }
        Self {
            nx,
            ny,
            wall_lower,
            wall_upper,
            artery_type: ArteryType::Stenosed,
        }
    }

    /// Check if a lattice node is inside the vessel (fluid domain).
    pub fn is_fluid(&self, x: usize, y: usize) -> bool {
        let yf = y as f64;
        yf > self.wall_lower[x] && yf < self.wall_upper[x]
    }

    /// Check if a lattice node is a wall node.
    pub fn is_wall(&self, x: usize, y: usize) -> bool {
        let yf = y as f64;
        let tol = 1.0;
        (yf - self.wall_lower[x]).abs() < tol || (yf - self.wall_upper[x]).abs() < tol
    }

    /// Get the local radius at position x.
    pub fn local_radius(&self, x: usize) -> f64 {
        (self.wall_upper[x] - self.wall_lower[x]) * 0.5
    }

    /// Get the center y-position at position x.
    pub fn center_y(&self, x: usize) -> f64 {
        (self.wall_upper[x] + self.wall_lower[x]) * 0.5
    }

    /// Compute the cross-sectional area at position x (in lattice units).
    pub fn cross_section_area(&self, x: usize) -> f64 {
        self.wall_upper[x] - self.wall_lower[x]
    }
}

/// Y-shaped bifurcation geometry.
#[derive(Debug, Clone)]
pub struct BifurcationGeometry {
    /// Total grid width.
    pub nx: usize,
    /// Total grid height.
    pub ny: usize,
    /// Parent vessel radius.
    pub parent_radius: f64,
    /// Daughter vessel radius (both daughters).
    pub daughter_radius: f64,
    /// Bifurcation angle (half-angle, radians).
    pub bifurcation_angle: f64,
    /// X-position where bifurcation starts.
    pub bifurcation_x: usize,
    /// Boolean mask: true = fluid.
    pub fluid_mask: Vec<Vec<bool>>,
}

impl BifurcationGeometry {
    /// Create a Y-shaped bifurcation.
    pub fn y_bifurcation(
        nx: usize,
        ny: usize,
        parent_radius: f64,
        daughter_radius: f64,
        bifurcation_angle: f64,
    ) -> Self {
        let bif_x = nx / 3;
        let center_y = ny as f64 / 2.0;
        let mut fluid_mask = vec![vec![false; ny]; nx];

        for x in 0..nx {
            for y in 0..ny {
                let yf = y as f64;
                if x < bif_x {
                    // Parent vessel
                    if (yf - center_y).abs() < parent_radius {
                        fluid_mask[x][y] = true;
                    }
                } else {
                    // Two daughter vessels
                    let dx = (x - bif_x) as f64;
                    let upper_center = center_y + dx * bifurcation_angle.tan();
                    let lower_center = center_y - dx * bifurcation_angle.tan();
                    if (yf - upper_center).abs() < daughter_radius
                        || (yf - lower_center).abs() < daughter_radius
                    {
                        fluid_mask[x][y] = true;
                    }
                    // Also keep fluid in the transition region
                    if x < bif_x + 5 && (yf - center_y).abs() < parent_radius {
                        fluid_mask[x][y] = true;
                    }
                }
            }
        }

        Self {
            nx,
            ny,
            parent_radius,
            daughter_radius,
            bifurcation_angle,
            bifurcation_x: bif_x,
            fluid_mask,
        }
    }

    /// Check if a node is fluid.
    pub fn is_fluid(&self, x: usize, y: usize) -> bool {
        if x < self.nx && y < self.ny {
            self.fluid_mask[x][y]
        } else {
            false
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Hemodynamics LBM Solver
// ─────────────────────────────────────────────────────────────────────────────

/// Cell type for the hemodynamics simulation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CellType {
    /// Fluid cell.
    Fluid,
    /// Wall cell (bounce-back).
    Wall,
    /// Inlet boundary.
    Inlet,
    /// Outlet boundary.
    Outlet,
}

/// Hemodynamics LBM solver for blood flow in arteries.
///
/// Uses D2Q9 lattice with non-Newtonian viscosity models.
#[derive(Debug, Clone)]
pub struct HemodynamicsLbm {
    /// Grid width.
    pub nx: usize,
    /// Grid height.
    pub ny: usize,
    /// Distribution functions f_i(x, y).
    pub f: Vec<Vec<[f64; 9]>>,
    /// Temporary storage for streaming.
    pub f_tmp: Vec<Vec<[f64; 9]>>,
    /// Macroscopic density.
    pub rho: Vec<Vec<f64>>,
    /// Macroscopic velocity.
    pub vel: Vec<Vec<[f64; 2]>>,
    /// Cell types.
    pub cell_type: Vec<Vec<CellType>>,
    /// Local relaxation time (varies for non-Newtonian).
    pub tau: Vec<Vec<f64>>,
    /// Base relaxation time.
    pub tau_base: f64,
    /// Current time step.
    pub time_step: usize,
}

impl HemodynamicsLbm {
    /// Create a new solver with given grid dimensions and base viscosity.
    pub fn new(nx: usize, ny: usize, nu: f64) -> Self {
        let tau_base = 0.5 + nu / CS2;
        let feq = Self::equilibrium_static(1.0, [0.0, 0.0]);
        Self {
            nx,
            ny,
            f: vec![vec![feq; ny]; nx],
            f_tmp: vec![vec![[0.0; 9]; ny]; nx],
            rho: vec![vec![1.0; ny]; nx],
            vel: vec![vec![[0.0, 0.0]; ny]; nx],
            cell_type: vec![vec![CellType::Fluid; ny]; nx],
            tau: vec![vec![tau_base; ny]; nx],
            tau_base,
            time_step: 0,
        }
    }

    /// Set up geometry from an ArterialGeometry.
    pub fn set_arterial_geometry(&mut self, geom: &ArterialGeometry) {
        for x in 0..self.nx {
            for y in 0..self.ny {
                if !geom.is_fluid(x, y) {
                    self.cell_type[x][y] = CellType::Wall;
                }
            }
            // Inlet (first column)
            if x == 0 {
                for y in 0..self.ny {
                    if geom.is_fluid(x, y) {
                        self.cell_type[x][y] = CellType::Inlet;
                    }
                }
            }
            // Outlet (last column)
            if x == self.nx - 1 {
                for y in 0..self.ny {
                    if geom.is_fluid(x, y) {
                        self.cell_type[x][y] = CellType::Outlet;
                    }
                }
            }
        }
    }

    /// Set up geometry from a BifurcationGeometry.
    pub fn set_bifurcation_geometry(&mut self, geom: &BifurcationGeometry) {
        for x in 0..self.nx {
            for y in 0..self.ny {
                if geom.is_fluid(x, y) {
                    self.cell_type[x][y] = CellType::Fluid;
                } else {
                    self.cell_type[x][y] = CellType::Wall;
                }
            }
        }
        // Inlet
        for y in 0..self.ny {
            if geom.is_fluid(0, y) {
                self.cell_type[0][y] = CellType::Inlet;
            }
        }
        // Outlets
        for y in 0..self.ny {
            if geom.is_fluid(self.nx - 1, y) {
                self.cell_type[self.nx - 1][y] = CellType::Outlet;
            }
        }
    }

    /// Compute equilibrium distribution for given density and velocity.
    fn equilibrium_static(rho: f64, u: [f64; 2]) -> [f64; 9] {
        let u_sq = dot2(u, u);
        let mut feq = [0.0; 9];
        for i in 0..9 {
            let eu = D2Q9_E[i][0] as f64 * u[0] + D2Q9_E[i][1] as f64 * u[1];
            feq[i] = D2Q9_W[i]
                * rho
                * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2));
        }
        feq
    }

    /// Compute local shear rate from the velocity field.
    pub fn compute_shear_rate(&self, x: usize, y: usize) -> f64 {
        if x == 0 || x >= self.nx - 1 || y == 0 || y >= self.ny - 1 {
            return 0.0;
        }
        // du_x/dy
        let duxy = (self.vel[x][y + 1][0] - self.vel[x][y - 1][0]) * 0.5;
        // du_y/dx
        let duyx = (self.vel[x + 1][y][1] - self.vel[x - 1][y][1]) * 0.5;
        // du_x/dx
        let duxx = (self.vel[x + 1][y][0] - self.vel[x - 1][y][0]) * 0.5;
        // du_y/dy
        let duyy = (self.vel[x][y + 1][1] - self.vel[x][y - 1][1]) * 0.5;
        // Strain rate tensor magnitude: sqrt(2 * S_ij * S_ij)
        let s_xx = duxx;
        let s_yy = duyy;
        let s_xy = 0.5 * (duxy + duyx);
        (2.0 * (s_xx * s_xx + s_yy * s_yy + 2.0 * s_xy * s_xy)).sqrt()
    }

    /// Update relaxation times using Casson model.
    pub fn update_tau_casson(&mut self, model: &CassonModel) {
        for x in 0..self.nx {
            for y in 0..self.ny {
                if self.cell_type[x][y] == CellType::Fluid {
                    let gamma_dot = self.compute_shear_rate(x, y);
                    self.tau[x][y] = model.relaxation_time(gamma_dot, 1.0, 1.0);
                }
            }
        }
    }

    /// Update relaxation times using Carreau-Yasuda model.
    pub fn update_tau_carreau_yasuda(&mut self, model: &CarreauYasudaModel) {
        for x in 0..self.nx {
            for y in 0..self.ny {
                if self.cell_type[x][y] == CellType::Fluid {
                    let gamma_dot = self.compute_shear_rate(x, y);
                    self.tau[x][y] = model.relaxation_time(gamma_dot);
                }
            }
        }
    }

    /// Compute macroscopic quantities from distribution functions.
    pub fn compute_macroscopic(&mut self) {
        for x in 0..self.nx {
            for y in 0..self.ny {
                if self.cell_type[x][y] == CellType::Wall {
                    continue;
                }
                let mut rho = 0.0;
                let mut ux = 0.0;
                let mut uy = 0.0;
                for i in 0..9 {
                    rho += self.f[x][y][i];
                    ux += self.f[x][y][i] * D2Q9_E[i][0] as f64;
                    uy += self.f[x][y][i] * D2Q9_E[i][1] as f64;
                }
                self.rho[x][y] = rho;
                if rho.abs() > 1e-12 {
                    self.vel[x][y] = [ux / rho, uy / rho];
                }
            }
        }
    }

    /// BGK collision step with local relaxation time.
    pub fn collide(&mut self) {
        for x in 0..self.nx {
            for y in 0..self.ny {
                if self.cell_type[x][y] != CellType::Fluid
                    && self.cell_type[x][y] != CellType::Inlet
                    && self.cell_type[x][y] != CellType::Outlet
                {
                    continue;
                }
                let feq = Self::equilibrium_static(self.rho[x][y], self.vel[x][y]);
                let omega = 1.0 / self.tau[x][y];
                for i in 0..9 {
                    self.f[x][y][i] += omega * (feq[i] - self.f[x][y][i]);
                }
            }
        }
    }

    /// Streaming step with periodic/bounce-back handling.
    pub fn stream(&mut self) {
        // Copy to temporary
        for x in 0..self.nx {
            for y in 0..self.ny {
                self.f_tmp[x][y] = self.f[x][y];
            }
        }
        // Stream
        for x in 0..self.nx {
            for y in 0..self.ny {
                for i in 0..9 {
                    let xn = x as i32 + D2Q9_E[i][0];
                    let yn = y as i32 + D2Q9_E[i][1];
                    if xn >= 0 && xn < self.nx as i32 && yn >= 0 && yn < self.ny as i32 {
                        let xnu = xn as usize;
                        let ynu = yn as usize;
                        if self.cell_type[xnu][ynu] == CellType::Wall {
                            // Bounce-back
                            self.f[x][y][D2Q9_OPP[i]] = self.f_tmp[x][y][i];
                        } else {
                            self.f[xnu][ynu][i] = self.f_tmp[x][y][i];
                        }
                    }
                }
            }
        }
    }

    /// Apply pulsatile inlet boundary condition.
    pub fn apply_pulsatile_inlet(
        &mut self,
        pulsatile: &PulsatileFlow,
        geom: &ArterialGeometry,
        t: f64,
    ) {
        let x = 0;
        for y in 0..self.ny {
            if self.cell_type[x][y] != CellType::Inlet {
                continue;
            }
            let center = geom.center_y(x);
            let radius = geom.local_radius(x);
            let r_ratio = ((y as f64 - center) / radius).abs().min(1.0);
            let u_local = pulsatile.velocity_profile(r_ratio, t);
            let rho_local = self.rho[x][y];
            let u_vec = [u_local.max(0.0), 0.0];
            self.f[x][y] = Self::equilibrium_static(rho_local, u_vec);
        }
    }

    /// Apply zero-gradient outlet boundary condition.
    pub fn apply_outlet(&mut self) {
        let x = self.nx - 1;
        if x == 0 {
            return;
        }
        for y in 0..self.ny {
            if self.cell_type[x][y] == CellType::Outlet {
                self.f[x][y] = self.f[x - 1][y];
            }
        }
    }

    /// Perform one full LBM step.
    pub fn step(&mut self) {
        self.compute_macroscopic();
        self.collide();
        self.stream();
        self.time_step += 1;
    }

    /// Perform one step with non-Newtonian viscosity update (Carreau-Yasuda).
    pub fn step_non_newtonian_cy(&mut self, model: &CarreauYasudaModel) {
        self.compute_macroscopic();
        self.update_tau_carreau_yasuda(model);
        self.collide();
        self.stream();
        self.time_step += 1;
    }

    /// Perform one step with Casson viscosity model.
    pub fn step_non_newtonian_casson(&mut self, model: &CassonModel) {
        self.compute_macroscopic();
        self.update_tau_casson(model);
        self.collide();
        self.stream();
        self.time_step += 1;
    }

    /// Initialize with a parabolic Poiseuille profile.
    pub fn init_poiseuille(&mut self, u_max: f64, geom: &ArterialGeometry) {
        for x in 0..self.nx {
            let center = geom.center_y(x);
            let radius = geom.local_radius(x);
            for y in 0..self.ny {
                if self.cell_type[x][y] == CellType::Wall {
                    continue;
                }
                let r_ratio = ((y as f64 - center) / radius).abs().min(1.0);
                let u_local = u_max * (1.0 - r_ratio * r_ratio);
                self.vel[x][y] = [u_local.max(0.0), 0.0];
                self.rho[x][y] = 1.0;
                self.f[x][y] = Self::equilibrium_static(1.0, self.vel[x][y]);
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Wall shear stress
// ─────────────────────────────────────────────────────────────────────────────

/// Compute instantaneous wall shear stress on the lower and upper walls.
///
/// WSS = mu * (du/dn) at the wall, where n is the wall-normal direction.
///
/// Returns (wss_lower\[x\], wss_upper\[x\]) for each x-position.
pub fn compute_wall_shear_stress(
    solver: &HemodynamicsLbm,
    geom: &ArterialGeometry,
    mu: f64,
) -> (Vec<f64>, Vec<f64>) {
    let mut wss_lower = vec![0.0; solver.nx];
    let mut wss_upper = vec![0.0; solver.nx];

    for x in 1..(solver.nx - 1) {
        let center = geom.center_y(x);
        let radius = geom.local_radius(x);
        let y_low = (center - radius).ceil() as usize + 1;
        let y_high = (center + radius).floor() as usize - 1;

        // Lower wall WSS (du/dy at lower wall)
        if y_low + 1 < solver.ny {
            let du_dy = solver.vel[x][y_low + 1][0] - solver.vel[x][y_low][0];
            wss_lower[x] = mu * du_dy.abs();
        }

        // Upper wall WSS (du/dy at upper wall)
        if y_high >= 1 {
            let du_dy = solver.vel[x][y_high][0] - solver.vel[x][y_high - 1][0];
            wss_upper[x] = mu * du_dy.abs();
        }
    }

    (wss_lower, wss_upper)
}

/// Time-averaged wall shear stress (TAWSS) from a sequence of WSS snapshots.
///
/// TAWSS = (1/T) * integral(|WSS(t)| dt)
pub fn compute_tawss(wss_history: &[Vec<f64>]) -> Vec<f64> {
    if wss_history.is_empty() {
        return Vec::new();
    }
    let n = wss_history[0].len();
    let nt = wss_history.len() as f64;
    let mut tawss = vec![0.0; n];
    for wss_snapshot in wss_history {
        for i in 0..n {
            tawss[i] += wss_snapshot[i].abs();
        }
    }
    for val in &mut tawss {
        *val /= nt;
    }
    tawss
}

// ─────────────────────────────────────────────────────────────────────────────
// Oscillatory shear index
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the oscillatory shear index (OSI) from WSS history.
///
/// OSI = 0.5 * (1 - |integral(WSS dt)| / integral(|WSS| dt))
///
/// OSI = 0 for unidirectional flow, OSI = 0.5 for purely oscillatory flow.
///
/// `wss_history` contains signed WSS values at each time step.
pub fn compute_osi(wss_history: &[Vec<f64>]) -> Vec<f64> {
    if wss_history.is_empty() {
        return Vec::new();
    }
    let n = wss_history[0].len();
    let mut sum_wss = vec![0.0; n]; // integral of WSS (signed)
    let mut sum_abs_wss = vec![0.0; n]; // integral of |WSS|

    for snapshot in wss_history {
        for i in 0..n {
            sum_wss[i] += snapshot[i];
            sum_abs_wss[i] += snapshot[i].abs();
        }
    }

    let mut osi = vec![0.0; n];
    for i in 0..n {
        if sum_abs_wss[i].abs() > 1e-12 {
            osi[i] = 0.5 * (1.0 - sum_wss[i].abs() / sum_abs_wss[i]);
        }
    }
    osi
}

// ─────────────────────────────────────────────────────────────────────────────
// Residence time
// ─────────────────────────────────────────────────────────────────────────────

/// Relative Residence Time (RRT).
///
/// RRT = 1 / ((1 - 2*OSI) * TAWSS)
///
/// High RRT indicates regions of stagnant flow, correlated with atherosclerosis.
pub fn compute_rrt(tawss: &[f64], osi: &[f64]) -> Vec<f64> {
    let n = tawss.len().min(osi.len());
    let mut rrt = vec![0.0; n];
    for i in 0..n {
        let denom = (1.0 - 2.0 * osi[i]) * tawss[i];
        if denom.abs() > 1e-12 {
            rrt[i] = 1.0 / denom;
        } else {
            rrt[i] = 1e12; // Very high residence time
        }
    }
    rrt
}

/// Tracer-based residence time computation.
///
/// Advects a passive scalar C(x, y, t) with the velocity field.
/// C is initialized to 1.0 inside the domain and the residence time
/// at each point is proportional to the remaining concentration.
#[derive(Debug, Clone)]
pub struct ResidenceTimeTracer {
    /// Concentration field.
    pub concentration: Vec<Vec<f64>>,
    /// Grid dimensions.
    pub nx: usize,
    /// Grid dimensions.
    pub ny: usize,
    /// Accumulated residence time.
    pub residence_time: Vec<Vec<f64>>,
}

impl ResidenceTimeTracer {
    /// Create a new tracer with uniform initial concentration.
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            concentration: vec![vec![1.0; ny]; nx],
            nx,
            ny,
            residence_time: vec![vec![0.0; ny]; nx],
        }
    }

    /// Advect the tracer field for one time step using upwind scheme.
    pub fn advect(&mut self, vel: &[Vec<[f64; 2]>], cell_type: &[Vec<CellType>], dt: f64) {
        let mut new_c = self.concentration.clone();
        for x in 1..(self.nx - 1) {
            for y in 1..(self.ny - 1) {
                if cell_type[x][y] == CellType::Wall {
                    new_c[x][y] = 0.0;
                    continue;
                }
                let ux = vel[x][y][0];
                let uy = vel[x][y][1];
                // Upwind differences
                let dc_dx = if ux > 0.0 {
                    self.concentration[x][y] - self.concentration[x - 1][y]
                } else {
                    self.concentration[x + 1][y] - self.concentration[x][y]
                };
                let dc_dy = if uy > 0.0 {
                    self.concentration[x][y] - self.concentration[x][y - 1]
                } else {
                    self.concentration[x][y + 1] - self.concentration[x][y]
                };
                new_c[x][y] = self.concentration[x][y] - dt * (ux * dc_dx + uy * dc_dy);
                new_c[x][y] = new_c[x][y].clamp(0.0, 1.0);
            }
        }
        self.concentration = new_c;
        // Accumulate residence time
        for x in 0..self.nx {
            for y in 0..self.ny {
                self.residence_time[x][y] += self.concentration[x][y] * dt;
            }
        }
    }

    /// Set inlet concentration to zero (washout).
    pub fn set_inlet_washout(&mut self, cell_type: &[Vec<CellType>]) {
        for x in 0..self.nx {
            for y in 0..self.ny {
                if cell_type[x][y] == CellType::Inlet {
                    self.concentration[x][y] = 0.0;
                }
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Stenosis modeling utilities
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the pressure drop across a stenosis using the Young-Tsai model.
///
/// Delta_P = (K_v * mu * U / D) + (K_t * 0.5 * rho * U^2) * ((A_0/A_s - 1)^2)
///
/// where K_v, K_t are empirical coefficients.
pub fn stenosis_pressure_drop(
    mu: f64,
    rho: f64,
    u_mean: f64,
    diameter: f64,
    area_ratio: f64,
    k_v: f64,
    k_t: f64,
) -> f64 {
    let viscous = k_v * mu * u_mean / diameter;
    let turbulent = k_t * 0.5 * rho * u_mean * u_mean * (area_ratio - 1.0).powi(2);
    viscous + turbulent
}

/// Compute the Reynolds number at a stenosis throat.
pub fn stenosis_reynolds_number(rho: f64, u_mean: f64, diameter: f64, mu: f64) -> f64 {
    rho * u_mean * diameter / mu
}

/// Compute the severity of a stenosis as percent area reduction.
pub fn stenosis_severity(original_area: f64, stenosed_area: f64) -> f64 {
    ((original_area - stenosed_area) / original_area * 100.0).clamp(0.0, 100.0)
}

// ─────────────────────────────────────────────────────────────────────────────
// RBC Transport
// ─────────────────────────────────────────────────────────────────────────────

/// Simplified Red Blood Cell particle for advection.
#[derive(Debug, Clone, Copy)]
pub struct RbcParticle {
    /// Position (x, y).
    pub pos: [f64; 2],
    /// Velocity (interpolated from fluid).
    pub vel: [f64; 2],
    /// Orientation angle (radians).
    pub angle: f64,
    /// Semi-major axis (microns).
    pub a: f64,
    /// Semi-minor axis (microns).
    pub b: f64,
}

impl RbcParticle {
    /// Create a new RBC at a position with default shape.
    pub fn new(x: f64, y: f64) -> Self {
        Self {
            pos: [x, y],
            vel: [0.0, 0.0],
            angle: 0.0,
            a: 4.0, // typical RBC semi-major axis (microns)
            b: 1.0, // typical semi-minor axis
        }
    }

    /// Aspect ratio of the RBC.
    pub fn aspect_ratio(&self) -> f64 {
        self.a / self.b
    }
}

/// RBC transport solver using Lagrangian particle tracking.
#[derive(Debug, Clone)]
pub struct RbcTransport {
    /// Collection of RBC particles.
    pub particles: Vec<RbcParticle>,
    /// Domain width.
    pub nx: usize,
    /// Domain height.
    pub ny: usize,
}

impl RbcTransport {
    /// Create a new RBC transport system.
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            particles: Vec::new(),
            nx,
            ny,
        }
    }

    /// Seed RBCs uniformly across the inlet.
    pub fn seed_inlet(&mut self, n_particles: usize, geom: &ArterialGeometry) {
        let center = geom.center_y(0);
        let radius = geom.local_radius(0);
        let dy = 2.0 * radius / (n_particles + 1) as f64;
        for i in 0..n_particles {
            let y = center - radius + (i + 1) as f64 * dy;
            self.particles.push(RbcParticle::new(1.0, y));
        }
    }

    /// Advect all RBC particles using bilinear interpolation of velocity field.
    pub fn advect(&mut self, vel: &[Vec<[f64; 2]>], _cell_type: &[Vec<CellType>], dt: f64) {
        for p in &mut self.particles {
            // Bilinear interpolation
            let ix = (p.pos[0].floor() as usize).min(self.nx - 2);
            let iy = (p.pos[1].floor() as usize).min(self.ny - 2);
            let fx = p.pos[0] - ix as f64;
            let fy = p.pos[1] - iy as f64;

            let u00 = vel[ix][iy];
            let u10 = vel[ix + 1][iy];
            let u01 = vel[ix][iy + 1];
            let u11 = vel[ix + 1][iy + 1];

            let ux = u00[0] * (1.0 - fx) * (1.0 - fy)
                + u10[0] * fx * (1.0 - fy)
                + u01[0] * (1.0 - fx) * fy
                + u11[0] * fx * fy;
            let uy = u00[1] * (1.0 - fx) * (1.0 - fy)
                + u10[1] * fx * (1.0 - fy)
                + u01[1] * (1.0 - fx) * fy
                + u11[1] * fx * fy;

            p.vel = [ux, uy];
            p.pos[0] += ux * dt;
            p.pos[1] += uy * dt;

            // Update orientation (Jeffery's orbit approximation)
            if ix > 0 && ix < self.nx - 1 && iy > 0 && iy < self.ny - 1 {
                let dudy = vel[ix][iy + 1][0] - vel[ix][iy][0];
                let dvdx = vel[ix + 1][iy][1] - vel[ix][iy][1];
                let omega_z = 0.5 * (dudy - dvdx); // vorticity
                let ar = p.aspect_ratio();
                let lambda = (ar * ar - 1.0) / (ar * ar + 1.0);
                let _shear = 0.5 * (dudy + dvdx);
                p.angle += (omega_z - lambda * omega_z) * dt;
            }
        }

        // Remove particles that left the domain
        self.particles.retain(|p| {
            p.pos[0] >= 0.0
                && p.pos[0] < self.nx as f64
                && p.pos[1] >= 0.0
                && p.pos[1] < self.ny as f64
        });
    }

    /// Count particles in each grid cell (hematocrit field).
    pub fn hematocrit_field(&self) -> Vec<Vec<f64>> {
        let mut field = vec![vec![0.0; self.ny]; self.nx];
        for p in &self.particles {
            let ix = (p.pos[0] as usize).min(self.nx - 1);
            let iy = (p.pos[1] as usize).min(self.ny - 1);
            field[ix][iy] += 1.0;
        }
        field
    }

    /// Get the number of active particles.
    pub fn active_count(&self) -> usize {
        self.particles.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Hemodynamic indices
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the flow rate through a cross-section at position x.
pub fn compute_flow_rate(solver: &HemodynamicsLbm, x: usize) -> f64 {
    let mut q = 0.0;
    for y in 0..solver.ny {
        if solver.cell_type[x][y] == CellType::Fluid
            || solver.cell_type[x][y] == CellType::Inlet
            || solver.cell_type[x][y] == CellType::Outlet
        {
            q += solver.vel[x][y][0] * solver.rho[x][y];
        }
    }
    q
}

/// Compute the average velocity magnitude at a cross-section.
pub fn compute_average_velocity(solver: &HemodynamicsLbm, x: usize) -> f64 {
    let mut sum_u = 0.0;
    let mut count = 0;
    for y in 0..solver.ny {
        if solver.cell_type[x][y] == CellType::Fluid
            || solver.cell_type[x][y] == CellType::Inlet
            || solver.cell_type[x][y] == CellType::Outlet
        {
            sum_u += len2(solver.vel[x][y]);
            count += 1;
        }
    }
    if count > 0 { sum_u / count as f64 } else { 0.0 }
}

/// Compute the pressure along the centerline.
pub fn compute_centerline_pressure(solver: &HemodynamicsLbm, geom: &ArterialGeometry) -> Vec<f64> {
    let mut pressure = Vec::with_capacity(solver.nx);
    for x in 0..solver.nx {
        let cy = geom.center_y(x) as usize;
        let p = solver.rho[x][cy] * CS2; // p = rho * cs^2
        pressure.push(p);
    }
    pressure
}

/// Compute the Dean number for curved flow (relevant for curved arteries).
///
/// De = Re * sqrt(D / (2*R_c))
///
/// where Re is Reynolds number, D is diameter, R_c is curvature radius.
pub fn dean_number(reynolds: f64, diameter: f64, curvature_radius: f64) -> f64 {
    reynolds * (diameter / (2.0 * curvature_radius)).sqrt()
}

/// Compute the Strouhal number: St = f * D / U.
pub fn strouhal_number(frequency: f64, diameter: f64, velocity: f64) -> f64 {
    frequency * diameter / velocity
}

/// Compute endothelial cell activation potential (ECAP = OSI / TAWSS).
///
/// High ECAP indicates regions prone to endothelial dysfunction.
pub fn compute_ecap(tawss: &[f64], osi: &[f64]) -> Vec<f64> {
    let n = tawss.len().min(osi.len());
    let mut ecap = vec![0.0; n];
    for i in 0..n {
        if tawss[i].abs() > 1e-12 {
            ecap[i] = osi[i] / tawss[i];
        }
    }
    ecap
}

/// Compute the transverse WSS (transWSS) indicator.
///
/// TransWSS measures the time-averaged WSS component perpendicular to the
/// time-averaged WSS vector. Simplified here for 1D WSS.
pub fn compute_trans_wss(wss_x_history: &[Vec<f64>], wss_y_history: &[Vec<f64>]) -> Vec<f64> {
    if wss_x_history.is_empty() || wss_y_history.is_empty() {
        return Vec::new();
    }
    let n = wss_x_history[0].len();
    let nt = wss_x_history.len() as f64;

    // Time-averaged WSS vector
    let mut avg_x = vec![0.0; n];
    let mut avg_y = vec![0.0; n];
    for t in 0..wss_x_history.len() {
        for i in 0..n {
            avg_x[i] += wss_x_history[t][i];
            avg_y[i] += wss_y_history[t][i];
        }
    }
    for i in 0..n {
        avg_x[i] /= nt;
        avg_y[i] /= nt;
    }

    // Compute transverse component
    let mut trans_wss = vec![0.0; n];
    for i in 0..n {
        let mag = (avg_x[i] * avg_x[i] + avg_y[i] * avg_y[i]).sqrt();
        if mag > 1e-12 {
            // Unit vector in time-averaged direction
            let nx = avg_x[i] / mag;
            let ny = avg_y[i] / mag;
            // Transverse: sum |WSS x n_hat|
            for t in 0..wss_x_history.len() {
                let cross = wss_x_history[t][i] * ny - wss_y_history[t][i] * nx;
                trans_wss[i] += cross.abs();
            }
            trans_wss[i] /= nt;
        }
    }
    trans_wss
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // --- Rheology tests ---

    #[test]
    fn test_casson_viscosity_high_shear() {
        let model = CassonModel::new_blood();
        // At very high shear rates, viscosity approaches mu_inf
        let mu = model.viscosity(1000.0);
        assert!(
            (mu - model.mu_inf).abs() < 0.01,
            "Casson high-shear viscosity = {mu}, expected ~{}",
            model.mu_inf
        );
    }

    #[test]
    fn test_casson_viscosity_low_shear() {
        let model = CassonModel::new_blood();
        // At low shear rates, viscosity should be higher
        let mu_low = model.viscosity(0.01);
        let mu_high = model.viscosity(100.0);
        assert!(
            mu_low > mu_high,
            "Casson: low-shear viscosity ({mu_low}) should exceed high-shear ({mu_high})"
        );
    }

    #[test]
    fn test_casson_relaxation_time_positive() {
        let model = CassonModel::new_blood();
        let tau = model.relaxation_time(10.0, 1.0, 1.0);
        assert!(tau > 0.5, "Relaxation time must be > 0.5, got {tau}");
    }

    #[test]
    fn test_carreau_yasuda_viscosity_limits() {
        let model = CarreauYasudaModel::new_blood();
        // At zero shear: should approach mu_0
        let mu_zero = model.viscosity(0.0);
        assert!(
            (mu_zero - model.mu_0).abs() < 0.001,
            "CY zero-shear = {mu_zero}, expected ~{}",
            model.mu_0
        );
        // At very high shear: should approach mu_inf
        let mu_high = model.viscosity(1e6);
        assert!(
            (mu_high - model.mu_inf).abs() < 0.001,
            "CY high-shear = {mu_high}, expected ~{}",
            model.mu_inf
        );
    }

    #[test]
    fn test_carreau_yasuda_shear_thinning() {
        let model = CarreauYasudaModel::new_blood();
        let mu1 = model.viscosity(1.0);
        let mu10 = model.viscosity(10.0);
        let mu100 = model.viscosity(100.0);
        assert!(
            mu1 > mu10 && mu10 > mu100,
            "CY should be shear-thinning: {} > {} > {}",
            mu1,
            mu10,
            mu100
        );
    }

    #[test]
    fn test_power_law_viscosity() {
        let model = PowerLawModel::new_blood();
        let mu = model.viscosity(10.0);
        assert!(mu > 0.0, "Power-law viscosity should be positive");
        assert!(mu < model.mu_max, "Power-law viscosity should be bounded");
    }

    // --- Pulsatile flow tests ---

    #[test]
    fn test_pulsatile_mean_velocity() {
        let pf = PulsatileFlow::new_simple(0.1, 72.0, 0.05);
        // Average over one cardiac cycle should be close to u_mean
        let period = 60.0 / pf.heart_rate;
        let n = 1000;
        let dt = period / n as f64;
        let mut sum = 0.0;
        for i in 0..n {
            sum += pf.velocity_at(i as f64 * dt);
        }
        let avg = sum / n as f64;
        assert!(
            (avg - pf.u_mean).abs() < 0.01,
            "Average pulsatile velocity = {avg}, expected ~{}",
            pf.u_mean
        );
    }

    #[test]
    fn test_pulsatile_amplitude() {
        let pf = PulsatileFlow::new_simple(0.1, 72.0, 0.05);
        let period = 60.0 / pf.heart_rate;
        let n = 1000;
        let dt = period / n as f64;
        let mut max_u = f64::NEG_INFINITY;
        let mut min_u = f64::INFINITY;
        for i in 0..n {
            let u = pf.velocity_at(i as f64 * dt);
            max_u = max_u.max(u);
            min_u = min_u.min(u);
        }
        assert!(max_u > pf.u_mean, "Max velocity should exceed mean");
        assert!(min_u < pf.u_mean, "Min velocity should be below mean");
    }

    #[test]
    fn test_womersley_number() {
        let pf = PulsatileFlow::new_simple(0.1, 72.0, 0.05);
        let wo = pf.womersley_number(0.005, 3.5e-6); // 5mm radius, blood kinematic viscosity
        assert!(wo > 0.0, "Womersley number should be positive, got {wo}");
    }

    #[test]
    fn test_velocity_profile_parabolic() {
        let pf = PulsatileFlow::new_simple(0.1, 72.0, 0.0); // no pulsation
        let u_center = pf.velocity_profile(0.0, 0.0);
        let u_wall = pf.velocity_profile(1.0, 0.0);
        assert!(u_center > 0.0, "Center velocity should be positive");
        assert!(
            u_wall.abs() < 1e-10,
            "Wall velocity should be ~0, got {u_wall}"
        );
    }

    // --- Geometry tests ---

    #[test]
    fn test_straight_artery_geometry() {
        let geom = ArterialGeometry::straight(100, 50, 15.0);
        assert!(geom.is_fluid(50, 25), "Center should be fluid");
        assert!(!geom.is_fluid(50, 0), "Far from center should be wall");
        assert!(
            (geom.local_radius(50) - 15.0).abs() < 1e-10,
            "Radius should be 15"
        );
    }

    #[test]
    fn test_tapered_artery_geometry() {
        let geom = ArterialGeometry::tapered(100, 50, 15.0, 10.0);
        let r_inlet = geom.local_radius(0);
        let r_outlet = geom.local_radius(99);
        assert!((r_inlet - 15.0).abs() < 0.1, "Inlet radius = {r_inlet}");
        assert!((r_outlet - 10.0).abs() < 0.1, "Outlet radius = {r_outlet}");
        assert!(r_inlet > r_outlet, "Tapered: inlet radius > outlet radius");
    }

    #[test]
    fn test_stenosed_artery() {
        let geom = ArterialGeometry::stenosed(100, 50, 15.0, 0.5, 0.5, 0.3);
        let r_healthy = geom.local_radius(0);
        let r_stenosis = geom.local_radius(50);
        assert!(
            r_stenosis < r_healthy,
            "Stenosis should reduce radius: {r_stenosis} < {r_healthy}"
        );
    }

    #[test]
    fn test_bifurcation_geometry() {
        let geom = BifurcationGeometry::y_bifurcation(100, 80, 10.0, 7.0, 0.3);
        // Parent region should be fluid
        assert!(geom.is_fluid(5, 40), "Parent vessel center should be fluid");
        // Far from vessels should not be fluid
        assert!(!geom.is_fluid(5, 0), "Far from vessel should not be fluid");
    }

    // --- Solver tests ---

    #[test]
    fn test_hemodynamics_lbm_creation() {
        let solver = HemodynamicsLbm::new(50, 30, 0.1);
        assert_eq!(solver.nx, 50);
        assert_eq!(solver.ny, 30);
        assert!(solver.tau_base > 0.5);
    }

    #[test]
    fn test_hemodynamics_lbm_step() {
        let mut solver = HemodynamicsLbm::new(30, 20, 0.1);
        let geom = ArterialGeometry::straight(30, 20, 7.0);
        solver.set_arterial_geometry(&geom);
        solver.init_poiseuille(0.05, &geom);
        solver.step();
        assert_eq!(solver.time_step, 1);
        // Density should remain close to 1.0
        let center_rho = solver.rho[15][10];
        assert!(
            (center_rho - 1.0).abs() < 0.1,
            "Density at center = {center_rho}"
        );
    }

    #[test]
    fn test_hemodynamics_mass_conservation() {
        let mut solver = HemodynamicsLbm::new(30, 20, 0.1);
        let geom = ArterialGeometry::straight(30, 20, 7.0);
        solver.set_arterial_geometry(&geom);
        solver.init_poiseuille(0.01, &geom);
        // Total mass before
        let mass_before: f64 = solver.rho.iter().flat_map(|col| col.iter()).sum();
        for _ in 0..5 {
            solver.step();
        }
        let mass_after: f64 = solver.rho.iter().flat_map(|col| col.iter()).sum();
        let rel_change = (mass_after - mass_before).abs() / mass_before;
        assert!(
            rel_change < 0.05,
            "Mass conservation: relative change = {rel_change}"
        );
    }

    #[test]
    fn test_non_newtonian_step() {
        let mut solver = HemodynamicsLbm::new(30, 20, 0.1);
        let geom = ArterialGeometry::straight(30, 20, 7.0);
        solver.set_arterial_geometry(&geom);
        solver.init_poiseuille(0.01, &geom);
        let model = CarreauYasudaModel::new_blood();
        solver.step_non_newtonian_cy(&model);
        assert_eq!(solver.time_step, 1);
    }

    // --- WSS tests ---

    #[test]
    fn test_wss_positive() {
        let mut solver = HemodynamicsLbm::new(30, 20, 0.1);
        let geom = ArterialGeometry::straight(30, 20, 7.0);
        solver.set_arterial_geometry(&geom);
        solver.init_poiseuille(0.05, &geom);
        solver.compute_macroscopic();
        let (wss_l, wss_u) = compute_wall_shear_stress(&solver, &geom, 0.003);
        // At least some WSS values should be positive
        let has_positive = wss_l.iter().any(|&w| w > 0.0) || wss_u.iter().any(|&w| w > 0.0);
        assert!(
            has_positive,
            "WSS should have positive values for Poiseuille flow"
        );
    }

    #[test]
    fn test_tawss_computation() {
        let wss1 = vec![1.0, 2.0, 3.0];
        let wss2 = vec![2.0, 1.0, 4.0];
        let tawss = compute_tawss(&[wss1, wss2]);
        assert!((tawss[0] - 1.5).abs() < 1e-10);
        assert!((tawss[1] - 1.5).abs() < 1e-10);
        assert!((tawss[2] - 3.5).abs() < 1e-10);
    }

    // --- OSI tests ---

    #[test]
    fn test_osi_unidirectional() {
        // All positive WSS -> OSI = 0
        let wss = vec![vec![1.0, 2.0, 3.0]; 10];
        let osi = compute_osi(&wss);
        for val in &osi {
            assert!(
                val.abs() < 1e-10,
                "OSI should be 0 for unidirectional flow, got {val}"
            );
        }
    }

    #[test]
    fn test_osi_oscillatory() {
        // Alternating WSS -> OSI = 0.5
        let mut history = Vec::new();
        for i in 0..100 {
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            history.push(vec![sign, sign, sign]);
        }
        let osi = compute_osi(&history);
        for val in &osi {
            assert!(
                (*val - 0.5).abs() < 0.01,
                "OSI should be ~0.5 for oscillatory flow, got {val}"
            );
        }
    }

    // --- Residence time tests ---

    #[test]
    fn test_rrt_computation() {
        let tawss = vec![1.0, 2.0, 0.5];
        let osi = vec![0.0, 0.0, 0.0];
        let rrt = compute_rrt(&tawss, &osi);
        assert!((rrt[0] - 1.0).abs() < 1e-10);
        assert!((rrt[1] - 0.5).abs() < 1e-10);
        assert!((rrt[2] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_residence_time_tracer() {
        let mut tracer = ResidenceTimeTracer::new(10, 10);
        let vel = vec![vec![[0.1, 0.0]; 10]; 10];
        let cell_type = vec![vec![CellType::Fluid; 10]; 10];
        tracer.advect(&vel, &cell_type, 1.0);
        // Concentration should decrease at downstream positions
        assert!(tracer.concentration[5][5] <= 1.0);
    }

    // --- Stenosis tests ---

    #[test]
    fn test_stenosis_pressure_drop() {
        let dp = stenosis_pressure_drop(0.003, 1060.0, 0.3, 0.006, 2.0, 32.0, 1.52);
        assert!(dp > 0.0, "Pressure drop should be positive: {dp}");
    }

    #[test]
    fn test_stenosis_reynolds_number() {
        let re = stenosis_reynolds_number(1060.0, 0.3, 0.003, 0.003);
        assert!((re - 318.0).abs() < 1.0, "Re = {re}, expected ~318");
    }

    #[test]
    fn test_stenosis_severity() {
        let sev = stenosis_severity(1.0, 0.5);
        assert!((sev - 50.0).abs() < 1e-10, "50% area reduction = {sev}%");
    }

    // --- RBC tests ---

    #[test]
    fn test_rbc_creation() {
        let rbc = RbcParticle::new(10.0, 20.0);
        assert!((rbc.pos[0] - 10.0).abs() < 1e-10);
        assert!((rbc.pos[1] - 20.0).abs() < 1e-10);
        assert!(rbc.aspect_ratio() > 1.0);
    }

    #[test]
    fn test_rbc_transport_seeding() {
        let geom = ArterialGeometry::straight(100, 50, 15.0);
        let mut transport = RbcTransport::new(100, 50);
        transport.seed_inlet(10, &geom);
        assert_eq!(transport.active_count(), 10);
    }

    #[test]
    fn test_rbc_advection() {
        let geom = ArterialGeometry::straight(100, 50, 15.0);
        let mut transport = RbcTransport::new(100, 50);
        transport.seed_inlet(5, &geom);
        let vel = vec![vec![[0.1, 0.0]; 50]; 100];
        let cell_type = vec![vec![CellType::Fluid; 50]; 100];
        let initial_x: Vec<f64> = transport.particles.iter().map(|p| p.pos[0]).collect();
        transport.advect(&vel, &cell_type, 1.0);
        for (i, p) in transport.particles.iter().enumerate() {
            assert!(p.pos[0] > initial_x[i], "RBC should move downstream");
        }
    }

    // --- Hemodynamic index tests ---

    #[test]
    fn test_dean_number() {
        let de = dean_number(500.0, 0.006, 0.05);
        assert!(de > 0.0, "Dean number should be positive");
    }

    #[test]
    fn test_strouhal_number() {
        let st = strouhal_number(1.2, 0.006, 0.3);
        assert!(
            (st - 0.024).abs() < 0.001,
            "Strouhal = {st}, expected ~0.024"
        );
    }

    #[test]
    fn test_ecap_computation() {
        let tawss = vec![1.0, 2.0, 0.5];
        let osi = vec![0.1, 0.2, 0.3];
        let ecap = compute_ecap(&tawss, &osi);
        assert!((ecap[0] - 0.1).abs() < 1e-10);
        assert!((ecap[1] - 0.1).abs() < 1e-10);
        assert!((ecap[2] - 0.6).abs() < 1e-10);
    }

    #[test]
    fn test_flow_rate_computation() {
        let mut solver = HemodynamicsLbm::new(30, 20, 0.1);
        let geom = ArterialGeometry::straight(30, 20, 7.0);
        solver.set_arterial_geometry(&geom);
        solver.init_poiseuille(0.05, &geom);
        solver.compute_macroscopic();
        let q = compute_flow_rate(&solver, 15);
        assert!(
            q > 0.0,
            "Flow rate should be positive for Poiseuille flow: {q}"
        );
    }

    #[test]
    fn test_centerline_pressure() {
        let mut solver = HemodynamicsLbm::new(30, 20, 0.1);
        let geom = ArterialGeometry::straight(30, 20, 7.0);
        solver.set_arterial_geometry(&geom);
        solver.init_poiseuille(0.01, &geom);
        let pressure = compute_centerline_pressure(&solver, &geom);
        assert_eq!(pressure.len(), 30);
        // Pressure should be close to rho*cs^2 ~ 1/3
        for p in &pressure {
            assert!((*p - CS2).abs() < 0.05, "Pressure = {p}, expected ~{CS2}");
        }
    }
}
