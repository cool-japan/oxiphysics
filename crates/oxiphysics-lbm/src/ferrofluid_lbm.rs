#![allow(clippy::needless_range_loop)]
#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ferrofluid Lattice Boltzmann Method (FF-LBM).
//!
//! This module implements a D2Q9 LBM solver for magnetic fluid dynamics:
//!
//! - [`FerrofluidParams`]: physical parameters (viscosity, permeability, etc.)
//! - [`MagneticField`]: external field with gradient tensor
//! - [`FerrofluidCell`]: per-cell state (density, velocity, magnetization, order parameter)
//! - [`FerrofluidLbm`]: BGK-LBM solver with Kelvin body force
//! - [`FerrofluidVortex`]: rotating-field driven vortex flow
//!
//! # Physics
//!
//! Ferrofluids are stable colloidal suspensions of single-domain magnetic
//! nanoparticles in a carrier liquid.  The dominant body force is the
//! Kelvin force  **f** = μ₀ (**M**·∇)**H**, which is incorporated as an
//! external forcing term in the BGK collision operator via the Guo scheme.
//!
//! Magnetization is computed with the Langevin model corrected for
//! chain-like particle aggregates that form along the field direction,
//! causing the well-known magnetoviscous effect.
//!
//! # References
//! - Rosensweig, R.E. (1985). *Ferrohydrodynamics*. Cambridge Univ. Press.
//! - Odenbach, S. (2002). *Magnetoviscous Effects in Ferrofluids*. Springer.
//! - Guo, Z., Zheng, C., & Shi, B. (2002). Discrete lattice effects on the
//!   forcing term in the lattice Boltzmann method. *Phys. Rev. E*, 65, 046308.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// D2Q9 lattice constants
// ---------------------------------------------------------------------------

/// Number of discrete velocities in D2Q9.
pub const D2Q9_Q: usize = 9;

/// D2Q9 discrete velocity x-components.
pub const D2Q9_EX: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];

/// D2Q9 discrete velocity y-components.
pub const D2Q9_EY: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];

/// D2Q9 lattice weights.
pub const D2Q9_W: [f64; 9] = [
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

/// Lattice speed of sound squared (c_s² = 1/3 in lattice units).
pub const CS2: f64 = 1.0 / 3.0;

/// Vacuum magnetic permeability μ₀ (H m⁻¹).
pub const MU_0: f64 = 1.256_637_061_4e-6;

// ---------------------------------------------------------------------------
// Parameter structs
// ---------------------------------------------------------------------------

/// Physical and numerical parameters for a ferrofluid simulation.
#[derive(Debug, Clone)]
pub struct FerrofluidParams {
    /// Dynamic viscosity of the carrier fluid η₀ (Pa·s).
    pub viscosity: f64,
    /// Relative magnetic permeability μ_r of the ferrofluid.
    pub magnetic_permeability: f64,
    /// Saturation magnetization M_s (A m⁻¹).
    pub magnetization_saturation: f64,
    /// Volume fraction φ of magnetic nanoparticles (0–1).
    pub particle_volume_fraction: f64,
    /// Mean chain length L (number of particles per aggregate).
    pub chain_length: f64,
    /// Magnetoviscous coupling constant λ (dimensionless).
    pub coupling_constant: f64,
}

impl Default for FerrofluidParams {
    fn default() -> Self {
        Self {
            viscosity: 1.0e-3,
            magnetic_permeability: 2.5,
            magnetization_saturation: 3.2e4,
            particle_volume_fraction: 0.05,
            chain_length: 3.0,
            coupling_constant: 1.2,
        }
    }
}

/// External magnetic field at a point in space.
#[derive(Debug, Clone, Default)]
pub struct MagneticField {
    /// Field vector **H** or **B** (A m⁻¹ or T depending on context).
    pub field: [f64; 3],
    /// Field gradient tensor ∂H_i/∂x_j  (row-major, 3×3).
    pub gradient: [[f64; 3]; 3],
}

/// State of a single LBM cell in the ferrofluid simulation.
#[derive(Debug, Clone)]
pub struct FerrofluidCell {
    /// Macroscopic mass density ρ (kg m⁻³ in physical units).
    pub density: f64,
    /// Macroscopic velocity **u** (m s⁻¹).
    pub velocity: [f64; 3],
    /// Local magnetization vector **M** (A m⁻¹).
    pub magnetization: [f64; 3],
    /// Nematic chain order parameter S ∈ \[0, 1\].
    pub chain_order_parameter: f64,
}

impl Default for FerrofluidCell {
    fn default() -> Self {
        Self {
            density: 1.0,
            velocity: [0.0; 3],
            magnetization: [0.0; 3],
            chain_order_parameter: 0.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Core physics functions
// ---------------------------------------------------------------------------

/// Kelvin body force **f** = μ₀ (**M**·∇)**H**.
///
/// # Arguments
/// * `magnetization` - magnetization vector **M** (A m⁻¹)
/// * `field_grad`    - gradient tensor ∂H_j/∂x_i (row = component, col = space)
///
/// # Returns
/// Force density vector (N m⁻³).
pub fn kelvin_force(magnetization: &[f64; 3], field_grad: &[[f64; 3]; 3]) -> [f64; 3] {
    let mut f = [0.0f64; 3];
    for i in 0..3 {
        for j in 0..3 {
            f[i] += magnetization[j] * field_grad[i][j];
        }
        f[i] *= MU_0;
    }
    f
}

/// Langevin magnetization model.
///
/// Computes the equilibrium magnetization magnitude for a given effective
/// field amplitude using:
///   M = M_s · L(α)  where  L(α) = coth(α) − 1/α  and  α = m·μ₀·H / (k_B T).
///
/// # Arguments
/// * `h_eff` - effective field magnitude H (A m⁻¹)
/// * `ms`    - saturation magnetization M_s (A m⁻¹)
/// * `alpha` - Langevin parameter  α = μ₀ m H / (k_B T)  (dimensionless)
///
/// # Returns
/// Magnetization magnitude (A m⁻¹).
pub fn langevin_magnetization(h_eff: f64, ms: f64, alpha: f64) -> f64 {
    let _ = h_eff; // h_eff is encoded in alpha by the caller
    if alpha.abs() < 1.0e-6 {
        return ms * alpha / 3.0; // linear regime L(α) ≈ α/3
    }
    let coth = 1.0 / alpha.tanh();
    ms * (coth - 1.0 / alpha)
}

/// Magnetic pressure (Maxwell stress normal component).
///
/// p_mag = (μ₀ / 2) · |**M**|² (scalar isotropic part relevant at interfaces).
///
/// # Arguments
/// * `m` - magnetization vector **M**
/// * `b` - magnetic flux density **B** (T)
///
/// # Returns
/// Magnetic pressure (Pa).
pub fn magnetic_pressure(m: &[f64; 3], b: &[f64; 3]) -> f64 {
    let _ = b;
    let m2: f64 = m.iter().map(|x| x * x).sum();
    0.5 * MU_0 * m2
}

/// Magnetic energy density  u_m = −μ₀ **M**·**H** / 2 (J m⁻³).
///
/// # Arguments
/// * `m` - magnetization **M** (A m⁻¹)
/// * `h` - field **H** (A m⁻¹)
///
/// # Returns
/// Energy density (J m⁻³) — negative because energy decreases in increasing field.
pub fn magnetic_energy_density(m: &[f64; 3], h: &[f64; 3]) -> f64 {
    let dot: f64 = m.iter().zip(h.iter()).map(|(a, b)| a * b).sum();
    -0.5 * MU_0 * dot
}

/// Rosensweig instability threshold.
///
/// The normal-field (Rosensweig) instability occurs when
///   M² > (2/μ₀) · (ρ g / k_c)^(1/2)  (simplified flat-interface criterion).
///
/// Returns `true` when the supplied magnetization magnitude exceeds the
/// threshold, indicating that the flat interface is unstable.
///
/// # Arguments
/// * `m_mag`   - magnetization magnitude (A m⁻¹)
/// * `density` - fluid density ρ (kg m⁻³)
/// * `gravity` - gravitational acceleration g (m s⁻²)
/// * `surface_tension` - surface tension γ (N m⁻¹)
///
/// # Returns
/// `true` when the Rosensweig instability threshold is exceeded.
pub fn rosensweig_instability(
    m_mag: f64,
    density: f64,
    gravity: f64,
    surface_tension: f64,
) -> bool {
    // Critical wave-number k_c = (ρ g / γ)^(1/2)
    let k_c = (density * gravity / surface_tension).sqrt();
    let threshold_sq = 2.0 * (density * gravity / k_c) / MU_0;
    m_mag * m_mag > threshold_sq
}

/// Chain-formation order parameter from field alignment.
///
/// A simple phenomenological model: S = tanh(λ · φ · L · α / 3)
/// where α is the Langevin parameter and λ is the coupling constant.
///
/// # Arguments
/// * `alpha`            - Langevin parameter α
/// * `phi`              - particle volume fraction
/// * `chain_len`        - mean chain length L
/// * `coupling`         - magnetoviscous coupling λ
///
/// # Returns
/// Order parameter S ∈ \[0, 1).
pub fn chain_order_parameter(alpha: f64, phi: f64, chain_len: f64, coupling: f64) -> f64 {
    (coupling * phi * chain_len * alpha / 3.0).tanh()
}

/// Magnetoviscous effective viscosity enhancement.
///
/// The rotational hindrance of particle chains increases the effective
/// viscosity:  η_eff = η₀ · (1 + 1.5 φ L S)
///
/// # Arguments
/// * `eta0`     - carrier viscosity η₀ (Pa·s)
/// * `phi`      - volume fraction
/// * `chain_len`- mean chain length
/// * `order`    - chain order parameter S
///
/// # Returns
/// Effective dynamic viscosity (Pa·s).
pub fn magnetoviscous_viscosity(eta0: f64, phi: f64, chain_len: f64, order: f64) -> f64 {
    eta0 * (1.0 + 1.5 * phi * chain_len * order)
}

// ---------------------------------------------------------------------------
// LBM helpers
// ---------------------------------------------------------------------------

/// Compute the BGK equilibrium distribution for D2Q9.
///
/// # Arguments
/// * `rho` - density
/// * `ux`  - x-velocity
/// * `uy`  - y-velocity
///
/// # Returns
/// Array of 9 equilibrium populations.
pub fn equilibrium_d2q9(rho: f64, ux: f64, uy: f64) -> [f64; 9] {
    let u2 = ux * ux + uy * uy;
    let mut feq = [0.0f64; 9];
    for i in 0..9 {
        let eu = D2Q9_EX[i] * ux + D2Q9_EY[i] * uy;
        feq[i] =
            D2Q9_W[i] * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
    }
    feq
}

/// Guo forcing term for a single population.
///
/// Δf_i = w_i (1 − 1/(2τ)) · \[(e_i − u)/c_s² + (e_i·u)e_i/c_s⁴\] · F
///
/// # Arguments
/// * `i`    - population index (0..9)
/// * `ux`   - macroscopic x-velocity
/// * `uy`   - macroscopic y-velocity
/// * `fx`   - body force x-component
/// * `fy`   - body force y-component
/// * `tau`  - BGK relaxation time
///
/// # Returns
/// Forcing correction Δf_i.
pub fn guo_forcing(i: usize, ux: f64, uy: f64, fx: f64, fy: f64, tau: f64) -> f64 {
    let ex = D2Q9_EX[i];
    let ey = D2Q9_EY[i];
    let eu = ex * ux + ey * uy;
    let ef = ex * fx + ey * fy;
    let uf = ux * fx + uy * fy;
    D2Q9_W[i] * (1.0 - 1.0 / (2.0 * tau)) * ((ef - uf) / CS2 + eu * ef / (CS2 * CS2))
}

// ---------------------------------------------------------------------------
// Main solver
// ---------------------------------------------------------------------------

/// Ferrofluid LBM solver on a 2-D D2Q9 lattice.
///
/// Solves the incompressible Navier-Stokes equations with the Kelvin body
/// force using the BGK collision operator and Guo forcing scheme.
pub struct FerrofluidLbm {
    /// Number of cells in x-direction.
    pub nx: usize,
    /// Number of cells in y-direction.
    pub ny: usize,
    /// Grid spacing Δx (m).
    pub dx: f64,
    /// Time step Δt (s).
    pub dt: f64,
    /// Physical parameters.
    pub params: FerrofluidParams,
    /// Distribution functions f\[cell * Q + pop\].
    f: Vec<f64>,
    /// Streaming buffer (same layout as f).
    f_tmp: Vec<f64>,
    /// Macroscopic cell states.
    pub cells: Vec<FerrofluidCell>,
    /// Applied magnetic field (uniform over domain for now).
    magnetic_field: MagneticField,
    /// BGK relaxation time τ.
    tau: f64,
}

impl FerrofluidLbm {
    /// Create a new ferrofluid LBM solver.
    ///
    /// # Arguments
    /// * `nx`     - grid width (cells)
    /// * `ny`     - grid height (cells)
    /// * `params` - ferrofluid physical parameters
    pub fn new(nx: usize, ny: usize, params: FerrofluidParams) -> Self {
        let n = nx * ny;
        // τ from kinematic viscosity: ν = c_s²(τ − 0.5) Δt  (lattice units Δx=Δt=1)
        // Using lattice units: τ = ν/c_s² + 0.5, ν ≡ params.viscosity (lattice units)
        let tau = params.viscosity / CS2 + 0.5;
        Self {
            nx,
            ny,
            dx: 1.0,
            dt: 1.0,
            params: params.clone(),
            f: vec![0.0; n * D2Q9_Q],
            f_tmp: vec![0.0; n * D2Q9_Q],
            cells: vec![FerrofluidCell::default(); n],
            magnetic_field: MagneticField::default(),
            tau,
        }
    }

    /// Return the flat cell index for grid position (x, y).
    #[inline]
    pub fn idx(&self, x: usize, y: usize) -> usize {
        y * self.nx + x
    }

    /// Initialize the domain with uniform density and velocity.
    ///
    /// # Arguments
    /// * `density`  - initial density ρ₀
    /// * `velocity` - initial velocity \[ux, uy\] (lattice units)
    pub fn initialize_uniform(&mut self, density: f64, velocity: [f64; 2]) -> &mut Self {
        let feq = equilibrium_d2q9(density, velocity[0], velocity[1]);
        for cell in self.cells.iter_mut() {
            cell.density = density;
            cell.velocity = [velocity[0], velocity[1], 0.0];
        }
        for c in 0..self.nx * self.ny {
            for i in 0..D2Q9_Q {
                self.f[c * D2Q9_Q + i] = feq[i];
            }
        }
        self
    }

    /// Set the applied magnetic field for the next step(s).
    ///
    /// # Arguments
    /// * `field` - magnetic field description (field vector + gradient tensor)
    pub fn apply_magnetic_field(&mut self, field: &MagneticField) -> &mut Self {
        self.magnetic_field = field.clone();
        // Update cell magnetization using Langevin model
        let ms = self.params.magnetization_saturation;
        let phi = self.params.particle_volume_fraction;
        let l = self.params.chain_length;
        let coupling = self.params.coupling_constant;
        let h_mag: f64 = self
            .magnetic_field
            .field
            .iter()
            .map(|v| v * v)
            .sum::<f64>()
            .sqrt();
        // alpha = μ₀ m H / k_B T  approximated as coupling * h_mag / ms
        let alpha = if ms > 1e-30 {
            coupling * h_mag / ms
        } else {
            0.0
        };
        let m_val = langevin_magnetization(h_mag, ms, alpha);
        let h_dir: [f64; 3] = if h_mag > 1e-30 {
            [
                self.magnetic_field.field[0] / h_mag,
                self.magnetic_field.field[1] / h_mag,
                self.magnetic_field.field[2] / h_mag,
            ]
        } else {
            [0.0; 3]
        };
        let order = chain_order_parameter(alpha, phi, l, coupling);
        for cell in self.cells.iter_mut() {
            cell.magnetization = [m_val * h_dir[0], m_val * h_dir[1], m_val * h_dir[2]];
            cell.chain_order_parameter = order;
        }
        self
    }

    /// BGK collision step with Kelvin body force (Guo scheme).
    pub fn collision_step(&mut self) -> &mut Self {
        let tau = self.tau;
        let omega = 1.0 / tau;
        let grad = self.magnetic_field.gradient;

        for c in 0..self.nx * self.ny {
            let rho = self.cells[c].density;
            let ux = self.cells[c].velocity[0];
            let uy = self.cells[c].velocity[1];
            let mag = self.cells[c].magnetization;

            // Kelvin body force (use only x,y components for 2-D)
            let kf = kelvin_force(&mag, &grad);
            let fx = kf[0];
            let fy = kf[1];

            let feq = equilibrium_d2q9(rho, ux, uy);
            for i in 0..D2Q9_Q {
                let fi = self.f[c * D2Q9_Q + i];
                let forcing = guo_forcing(i, ux, uy, fx, fy, tau);
                self.f[c * D2Q9_Q + i] = fi - omega * (fi - feq[i]) + forcing;
            }
        }
        self
    }

    /// Streaming step with periodic boundary conditions.
    pub fn streaming_step(&mut self) -> &mut Self {
        let nx = self.nx;
        let ny = self.ny;
        for y in 0..ny {
            for x in 0..nx {
                let c = y * nx + x;
                for i in 0..D2Q9_Q {
                    let ex = D2Q9_EX[i] as isize;
                    let ey = D2Q9_EY[i] as isize;
                    let xd = ((x as isize + ex).rem_euclid(nx as isize)) as usize;
                    let yd = ((y as isize + ey).rem_euclid(ny as isize)) as usize;
                    let cd = yd * nx + xd;
                    self.f_tmp[cd * D2Q9_Q + i] = self.f[c * D2Q9_Q + i];
                }
            }
        }
        std::mem::swap(&mut self.f, &mut self.f_tmp);

        // Update macroscopic quantities
        let phi = self.params.particle_volume_fraction;
        let l = self.params.chain_length;
        let coupling = self.params.coupling_constant;
        let h_mag: f64 = self
            .magnetic_field
            .field
            .iter()
            .map(|v| v * v)
            .sum::<f64>()
            .sqrt();
        let ms = self.params.magnetization_saturation;
        let alpha = if ms > 1e-30 {
            coupling * h_mag / ms
        } else {
            0.0
        };

        for c in 0..nx * ny {
            let mut rho = 0.0f64;
            let mut ux = 0.0f64;
            let mut uy = 0.0f64;
            for i in 0..D2Q9_Q {
                let fi = self.f[c * D2Q9_Q + i];
                rho += fi;
                ux += D2Q9_EX[i] * fi;
                uy += D2Q9_EY[i] * fi;
            }
            if rho > 1e-30 {
                ux /= rho;
                uy /= rho;
            }
            self.cells[c].density = rho;
            self.cells[c].velocity[0] = ux;
            self.cells[c].velocity[1] = uy;
            self.cells[c].chain_order_parameter = chain_order_parameter(alpha, phi, l, coupling);
        }
        self
    }

    /// Perform one full LBM step (collision then streaming).
    pub fn step(&mut self) -> &mut Self {
        self.collision_step();
        self.streaming_step();
        self
    }

    /// Return the average kinetic energy density over the domain.
    pub fn kinetic_energy_density(&self) -> f64 {
        let mut ke = 0.0;
        for cell in &self.cells {
            let u2 = cell.velocity[0].powi(2) + cell.velocity[1].powi(2);
            ke += 0.5 * cell.density * u2;
        }
        ke / (self.nx * self.ny) as f64
    }

    /// Return the average magnetic energy density over the domain.
    pub fn avg_magnetic_energy(&self) -> f64 {
        let h = &self.magnetic_field.field;
        let mut total = 0.0;
        for cell in &self.cells {
            total += magnetic_energy_density(&cell.magnetization, h).abs();
        }
        total / (self.nx * self.ny) as f64
    }
}

// ---------------------------------------------------------------------------
// FerrofluidVortex — rotating-field driven flow
// ---------------------------------------------------------------------------

/// Rotating magnetic field driven vortex model.
///
/// A uniform magnetic field rotating at angular frequency ω induces a
/// body-torque on the magnetic particles, which through viscous coupling
/// drives a swirling flow.  The induced angular velocity of the fluid
/// is estimated from the balance of magnetic torque and viscous drag.
#[derive(Debug, Clone)]
pub struct FerrofluidVortex {
    /// Rotation frequency ω (rad s⁻¹).
    pub omega: f64,
    /// Peak field amplitude H₀ (A m⁻¹).
    pub field_amplitude: f64,
    /// Ferrofluid physical parameters.
    pub params: FerrofluidParams,
    /// Current phase φ (rad).
    pub phase: f64,
}

impl FerrofluidVortex {
    /// Create a new vortex driver.
    ///
    /// # Arguments
    /// * `omega`           - field rotation frequency (rad s⁻¹)
    /// * `field_amplitude` - peak field H₀ (A m⁻¹)
    /// * `params`          - ferrofluid parameters
    pub fn new(omega: f64, field_amplitude: f64, params: FerrofluidParams) -> Self {
        Self {
            omega,
            field_amplitude,
            params,
            phase: 0.0,
        }
    }

    /// Return the instantaneous applied field at time t.
    ///
    /// **H**(t) = H₀ (cos(ωt), sin(ωt), 0)
    pub fn field_at_time(&self, t: f64) -> [f64; 3] {
        let angle = self.omega * t + self.phase;
        [
            self.field_amplitude * angle.cos(),
            self.field_amplitude * angle.sin(),
            0.0,
        ]
    }

    /// Estimate the steady-state azimuthal fluid velocity induced by the
    /// rotating field at radius r from the centre (simplified model).
    ///
    /// u_θ ≈ (3 φ η_mag ω r) / (2 η_eff)
    pub fn azimuthal_velocity(&self, r: f64) -> f64 {
        let phi = self.params.particle_volume_fraction;
        let eta0 = self.params.viscosity;
        let h0 = self.field_amplitude;
        let ms = self.params.magnetization_saturation;
        let coupling = self.params.coupling_constant;
        let alpha = if ms > 1e-30 { coupling * h0 / ms } else { 0.0 };
        let order = chain_order_parameter(alpha, phi, self.params.chain_length, coupling);
        let eta_eff = magnetoviscous_viscosity(eta0, phi, self.params.chain_length, order);
        let eta_mag = 1.5 * phi * eta0 * order;
        3.0 * phi * eta_mag * self.omega * r / (2.0 * eta_eff)
    }

    /// Advance the phase by one time step.
    pub fn advance(&mut self, dt: f64) {
        self.phase += self.omega * dt;
        self.phase = self.phase.rem_euclid(2.0 * PI);
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- physical function tests -------------------------------------------

    #[test]
    fn test_langevin_linear_regime() {
        // For small α, L(α) ≈ α/3
        let ms = 1.0e4;
        let alpha = 1.0e-4;
        let m = langevin_magnetization(0.0, ms, alpha);
        let expected = ms * alpha / 3.0;
        assert!((m - expected).abs() < 1e-6 * expected.abs().max(1.0));
    }

    #[test]
    fn test_langevin_saturation_limit() {
        // For large α, L(α) → 1 − 1/α, so M approaches M_s from below
        let ms = 3.2e4;
        let alpha = 1000.0;
        let m = langevin_magnetization(ms, ms, alpha);
        // L(1000) = 1 - 1/1000 = 0.999, so m = ms * 0.999
        assert!(m > 0.998 * ms);
        assert!(m <= ms);
    }

    #[test]
    fn test_langevin_zero_alpha() {
        let m = langevin_magnetization(0.0, 1.0e4, 0.0);
        assert!(m.abs() < 1.0e-6);
    }

    #[test]
    fn test_kelvin_force_zero_gradient() {
        let m = [1.0e4, 0.0, 0.0];
        let grad = [[0.0; 3]; 3];
        let f = kelvin_force(&m, &grad);
        assert_eq!(f, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_kelvin_force_uniform_gradient() {
        // M = M_x x̂, ∂H_x/∂x = g → F_x = μ₀ M_x g
        let m_x = 1.0e4;
        let g = 100.0;
        let m = [m_x, 0.0, 0.0];
        let mut grad = [[0.0f64; 3]; 3];
        grad[0][0] = g; // ∂H_x/∂x
        let f = kelvin_force(&m, &grad);
        let expected_fx = MU_0 * m_x * g;
        assert!((f[0] - expected_fx).abs() < 1e-20 * expected_fx.abs().max(1.0));
    }

    #[test]
    fn test_magnetic_pressure_positive() {
        let m = [1.0e4, 0.0, 0.0];
        let b = [0.0; 3];
        let p = magnetic_pressure(&m, &b);
        assert!(p > 0.0);
    }

    #[test]
    fn test_magnetic_pressure_scale() {
        let m = [1.0e4, 0.0, 0.0];
        let b = [0.0; 3];
        let p = magnetic_pressure(&m, &b);
        let expected = 0.5 * MU_0 * 1.0e8;
        assert!((p - expected).abs() < 1e-10 * expected);
    }

    #[test]
    fn test_magnetic_energy_density_negative() {
        let m = [1.0e4, 0.0, 0.0];
        let h = [1.0e4, 0.0, 0.0];
        let u = magnetic_energy_density(&m, &h);
        assert!(u < 0.0);
    }

    #[test]
    fn test_magnetic_energy_density_formula() {
        let m = [2.0, 0.0, 0.0];
        let h = [3.0, 0.0, 0.0];
        let u = magnetic_energy_density(&m, &h);
        let expected = -0.5 * MU_0 * 6.0;
        assert!((u - expected).abs() < 1e-20);
    }

    #[test]
    fn test_rosensweig_below_threshold() {
        // Very weak field → no instability
        let instab = rosensweig_instability(1.0, 1000.0, 9.81, 0.03);
        assert!(!instab);
    }

    #[test]
    fn test_rosensweig_above_threshold() {
        // Huge magnetization → instability
        let instab = rosensweig_instability(1.0e10, 1000.0, 9.81, 0.03);
        assert!(instab);
    }

    #[test]
    fn test_chain_order_zero_field() {
        let s = chain_order_parameter(0.0, 0.05, 3.0, 1.2);
        assert!(s.abs() < 1e-12);
    }

    #[test]
    fn test_chain_order_bounded() {
        // For large α the order parameter saturates to (at most) 1.0
        let s = chain_order_parameter(100.0, 0.1, 5.0, 2.0);
        assert!(s >= 0.0 && s <= 1.0);
    }

    #[test]
    fn test_magnetoviscous_no_chains() {
        let eta = magnetoviscous_viscosity(1e-3, 0.05, 1.0, 0.0);
        assert!((eta - 1e-3).abs() < 1e-15);
    }

    #[test]
    fn test_magnetoviscous_enhancement() {
        let eta0 = 1e-3;
        let eta = magnetoviscous_viscosity(eta0, 0.1, 3.0, 1.0);
        assert!(eta > eta0);
    }

    // ---- D2Q9 equilibrium tests --------------------------------------------

    #[test]
    fn test_equilibrium_mass_conservation() {
        let feq = equilibrium_d2q9(1.0, 0.1, 0.05);
        let rho: f64 = feq.iter().sum();
        assert!((rho - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_equilibrium_momentum_conservation() {
        let ux = 0.1;
        let uy = -0.05;
        let feq = equilibrium_d2q9(1.0, ux, uy);
        let jx: f64 = feq.iter().zip(D2Q9_EX.iter()).map(|(f, e)| f * e).sum();
        let jy: f64 = feq.iter().zip(D2Q9_EY.iter()).map(|(f, e)| f * e).sum();
        assert!((jx - ux).abs() < 1e-12);
        assert!((jy - uy).abs() < 1e-12);
    }

    #[test]
    fn test_equilibrium_at_rest() {
        let feq = equilibrium_d2q9(1.0, 0.0, 0.0);
        // At rest: populations equal weights times rho
        for i in 0..9 {
            assert!(
                (feq[i] - D2Q9_W[i]).abs() < 1e-14,
                "pop {i}: {:.15}",
                feq[i]
            );
        }
    }

    // ---- solver integration tests ------------------------------------------

    #[test]
    fn test_lbm_new() {
        let params = FerrofluidParams::default();
        let lbm = FerrofluidLbm::new(8, 8, params);
        assert_eq!(lbm.nx, 8);
        assert_eq!(lbm.ny, 8);
    }

    #[test]
    fn test_lbm_initialize_uniform() {
        let params = FerrofluidParams::default();
        let mut lbm = FerrofluidLbm::new(4, 4, params);
        lbm.initialize_uniform(1.0, [0.0, 0.0]);
        for cell in &lbm.cells {
            assert!((cell.density - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_lbm_step_density_conservation() {
        let params = FerrofluidParams::default();
        let mut lbm = FerrofluidLbm::new(8, 8, params);
        lbm.initialize_uniform(1.0, [0.01, 0.0]);
        let rho_before: f64 = lbm.cells.iter().map(|c| c.density).sum();
        lbm.step();
        let rho_after: f64 = lbm.cells.iter().map(|c| c.density).sum();
        assert!((rho_before - rho_after).abs() < 1e-8 * rho_before);
    }

    #[test]
    fn test_apply_magnetic_field_sets_magnetization() {
        let params = FerrofluidParams {
            magnetization_saturation: 1.0e4,
            coupling_constant: 1.0,
            ..Default::default()
        };
        let mut lbm = FerrofluidLbm::new(4, 4, params);
        lbm.initialize_uniform(1.0, [0.0, 0.0]);
        let field = MagneticField {
            field: [1.0e4, 0.0, 0.0],
            gradient: [[0.0; 3]; 3],
        };
        lbm.apply_magnetic_field(&field);
        // All cells should have nonzero magnetization in x
        for cell in &lbm.cells {
            assert!(cell.magnetization[0] > 0.0);
        }
    }

    #[test]
    fn test_kinetic_energy_at_rest() {
        let params = FerrofluidParams::default();
        let mut lbm = FerrofluidLbm::new(4, 4, params);
        lbm.initialize_uniform(1.0, [0.0, 0.0]);
        let ke = lbm.kinetic_energy_density();
        assert!(ke.abs() < 1e-12);
    }

    #[test]
    fn test_kinetic_energy_nonzero_velocity() {
        let params = FerrofluidParams::default();
        let mut lbm = FerrofluidLbm::new(4, 4, params);
        lbm.initialize_uniform(1.0, [0.1, 0.0]);
        let ke = lbm.kinetic_energy_density();
        assert!(ke > 0.0);
    }

    #[test]
    fn test_multiple_steps_stability() {
        let params = FerrofluidParams::default();
        let mut lbm = FerrofluidLbm::new(8, 8, params);
        lbm.initialize_uniform(1.0, [0.01, 0.005]);
        for _ in 0..20 {
            lbm.step();
        }
        // Check no NaN / Inf
        for cell in &lbm.cells {
            assert!(cell.density.is_finite());
            assert!(cell.velocity[0].is_finite());
            assert!(cell.velocity[1].is_finite());
        }
    }

    // ---- FerrofluidVortex tests --------------------------------------------

    #[test]
    fn test_vortex_field_amplitude() {
        let params = FerrofluidParams::default();
        let vortex = FerrofluidVortex::new(10.0, 1.0e4, params);
        let h = vortex.field_at_time(0.0);
        assert!((h[0] - 1.0e4).abs() < 1e-8);
        assert!(h[1].abs() < 1e-8);
    }

    #[test]
    fn test_vortex_field_orthogonal_at_quarter_period() {
        let params = FerrofluidParams::default();
        let omega = 2.0 * PI;
        let vortex = FerrofluidVortex::new(omega, 1.0e4, params);
        let t_quarter = 0.25; // T/4
        let h = vortex.field_at_time(t_quarter);
        assert!(h[0].abs() < 1e-8);
        assert!((h[1] - 1.0e4).abs() < 1e-8);
    }

    #[test]
    fn test_vortex_azimuthal_velocity_positive() {
        let params = FerrofluidParams::default();
        let vortex = FerrofluidVortex::new(100.0, 1.0e5, params);
        let u = vortex.azimuthal_velocity(0.01);
        assert!(u >= 0.0);
    }

    #[test]
    fn test_vortex_advance_phase() {
        let params = FerrofluidParams::default();
        let mut vortex = FerrofluidVortex::new(1.0, 1.0e4, params);
        vortex.advance(PI);
        assert!((vortex.phase - PI).abs() < 1e-12);
    }

    #[test]
    fn test_guo_forcing_zero_force() {
        let delta = guo_forcing(0, 0.0, 0.0, 0.0, 0.0, 1.0);
        assert!(delta.abs() < 1e-15);
    }

    #[test]
    fn test_avg_magnetic_energy_with_field() {
        let params = FerrofluidParams::default();
        let mut lbm = FerrofluidLbm::new(4, 4, params);
        lbm.initialize_uniform(1.0, [0.0, 0.0]);
        let field = MagneticField {
            field: [1.0e4, 0.0, 0.0],
            gradient: [[0.0; 3]; 3],
        };
        lbm.apply_magnetic_field(&field);
        let e = lbm.avg_magnetic_energy();
        assert!(e > 0.0);
    }
}
