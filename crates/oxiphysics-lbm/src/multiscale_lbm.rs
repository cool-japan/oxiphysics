// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Multiscale and hybrid LBM-continuum methods.
//!
//! This module provides methods that bridge different length and time scales:
//!
//! - [`HybridLbmNavier`]: Hybrid LBM / Navier-Stokes domain decomposition
//!   (D1Q3 BGK populations on the LBM side)
//! - [`AdaptiveLbm`]: Multi-level D1Q3 advection-diffusion LBM with subcycled
//!   refinement levels and Dupuis-Chopard transfer operators
//! - [`MesoscaleLbm`]: Nano-to-meso-to-macro three-scale coupling on a D1Q3
//!   hydrodynamic lattice
//! - [`LbmMolecularDynamicsCoupling`]: LBM / MD continuum-atomistic coupling
//! - [`ChapmanEnskogExpansion`]: Chapman-Enskog perturbative expansion
//! - [`multiscale_diffusion_coefficient`]: Combined diffusion coefficient
//!
//! The three lattice solvers share the D1Q3 kernel (`w = {2/3, 1/6, 1/6}`,
//! `c = {0, +1, -1}`, `cs² = 1/3`): populations are relaxed with BGK collision
//! towards a genuine equilibrium, streamed one node per step and reduced to
//! macroscopic fields by moments, so mass (and, where the model conserves it,
//! momentum) is preserved to round-off.
//!
//! # Physics Background
//!
//! Multiscale methods connect atomistic (MD), mesoscopic (LBM), and
//! macroscopic (Navier-Stokes/Euler) descriptions. The Schwarz alternating
//! method and flux-matching boundary conditions allow seamless information
//! transfer across scale interfaces.
//!
//! ## Chapman-Enskog Expansion
//! The Chapman-Enskog expansion recovers the Navier-Stokes equations from
//! the Boltzmann equation at low Knudsen number. At first order one obtains
//! the Newtonian viscous stress; at second order (Burnett) corrections arise.
//!
//! ## Irving-Kirkwood Flux
//! The Irving-Kirkwood procedure maps discrete MD particle data onto
//! continuum fields using a kernel localization function, preserving mass,
//! momentum, and energy conservation.
//!
//! # References
//! - Succi, S. (2001). *The Lattice Boltzmann Equation*. Oxford University Press.
//! - E, W. & Engquist, B. (2003). The heterogeneous multiscale methods.
//!   *Commun. Math. Sci.*, 1(1), 87–132.
//! - Chapman, S. & Cowling, T. G. (1970). *The Mathematical Theory of Non-Uniform
//!   Gases* (3rd ed.). Cambridge University Press.

mod adaptive;
mod d1q3;
mod hybrid;
mod mesoscale;
#[cfg(test)]
mod regression_tests;

pub use adaptive::AdaptiveLbm;
pub use hybrid::HybridLbmNavier;
pub use mesoscale::MesoscaleLbm;

use std::f64::consts::PI;

// ============================================================================
// Helper utilities
// ============================================================================

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

#[inline]
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}

// ============================================================================
// LbmMolecularDynamicsCoupling
// ============================================================================

/// Coupling between LBM (continuum mesoscale) and Molecular Dynamics (atomistic).
///
/// The Irving-Kirkwood procedure maps MD particle data onto LBM grid nodes
/// using a Gaussian kernel, preserving continuum conservation laws.
pub struct LbmMolecularDynamicsCoupling {
    /// LBM velocity field at grid nodes (x-component).
    pub lbm_vel: Vec<f64>,
    /// MD particle positions.
    pub md_positions: Vec<[f64; 3]>,
    /// MD particle velocities.
    pub md_velocities: Vec<[f64; 3]>,
    /// LBM grid spacing (m).
    pub dx: f64,
    /// Kernel width for coarse-graining (m).
    pub kernel_width: f64,
    /// MD particle masses (kg).
    pub md_masses: Vec<f64>,
    /// MD particle forces (N).
    pub md_forces: Vec<[f64; 3]>,
}

impl LbmMolecularDynamicsCoupling {
    /// Create a new LBM-MD coupling.
    pub fn new(n_lbm: usize, n_md: usize, dx: f64, kernel_width: f64) -> Self {
        Self {
            lbm_vel: vec![0.0; n_lbm],
            md_positions: vec![[0.0; 3]; n_md],
            md_velocities: vec![[0.0; 3]; n_md],
            dx,
            kernel_width,
            md_masses: vec![1.0; n_md],
            md_forces: vec![[0.0; 3]; n_md],
        }
    }

    /// Number of LBM grid nodes.
    pub fn n_lbm(&self) -> usize {
        self.lbm_vel.len()
    }

    /// Gaussian kernel value: φ(r) = exp(-r²/(2σ²)) / (σ√(2π)).
    fn gaussian_kernel(&self, r: f64) -> f64 {
        let sigma = self.kernel_width;
        if sigma < 1e-300 {
            return 0.0;
        }
        (-r * r / (2.0 * sigma * sigma)).exp() / (sigma * (2.0 * PI).sqrt())
    }

    /// Irving-Kirkwood stress flux: projects MD momentum onto LBM grid nodes
    /// via the Gaussian kernel and returns the coarse-grained momentum density.
    pub fn irving_kirkwood_flux(&self) -> Vec<f64> {
        let n = self.n_lbm();
        let mut flux = vec![0.0; n];
        for (idx, pos) in self.md_positions.iter().enumerate() {
            let mass = self.md_masses[idx];
            let vx = self.md_velocities[idx][0];
            for (gi, f) in flux.iter_mut().enumerate() {
                let x_grid = gi as f64 * self.dx;
                let r = (pos[0] - x_grid).abs();
                *f += mass * vx * self.gaussian_kernel(r);
            }
        }
        flux
    }

    /// Density projection: coarse-grain MD particle masses onto LBM grid.
    pub fn density_projection(&self) -> Vec<f64> {
        let n = self.n_lbm();
        let mut rho = vec![0.0; n];
        for (idx, pos) in self.md_positions.iter().enumerate() {
            let mass = self.md_masses[idx];
            for (gi, r_out) in rho.iter_mut().enumerate() {
                let x_grid = gi as f64 * self.dx;
                let r = (pos[0] - x_grid).abs();
                *r_out += mass * self.gaussian_kernel(r);
            }
        }
        rho
    }

    /// Velocity projection: coarse-grain MD velocities onto LBM grid (x-component).
    pub fn velocity_projection(&self) -> Vec<f64> {
        let rho = self.density_projection();
        let momentum = self.irving_kirkwood_flux();
        rho.iter()
            .zip(momentum.iter())
            .map(|(&r, &m)| if r.abs() > 1e-300 { m / r } else { 0.0 })
            .collect()
    }

    /// Force decomposition: decomposes total MD forces into conservative
    /// and dissipative components and returns the conservative part.
    ///
    /// Conservative force on particle i: F_c = Σ_{j≠i} (-∇V_ij).
    /// Here we return the magnitude of the total force as a surrogate.
    pub fn force_decomposition(&self) -> Vec<f64> {
        self.md_forces.iter().map(|f| len3(*f)).collect()
    }
}

// ============================================================================
// ChapmanEnskogExpansion
// ============================================================================

/// Chapman-Enskog perturbative expansion connecting Boltzmann to Navier-Stokes.
///
/// The expansion parameter is the Knudsen number Kn = λ/L, where λ is the
/// mean free path and L is a macroscopic length scale.
pub struct ChapmanEnskogExpansion {
    /// Knudsen number Kn = λ/L.
    pub knudsen: f64,
    /// Dynamic viscosity μ (Pa·s).
    pub viscosity: f64,
    /// Thermal conductivity κ (W/(m·K)).
    pub thermal_conductivity: f64,
    /// Specific heat at constant volume c_v (J/(kg·K)).
    pub cv: f64,
    /// Temperature T (K).
    pub temperature: f64,
    /// Reference density ρ₀ (kg/m³).
    pub density: f64,
}

impl ChapmanEnskogExpansion {
    /// Create a new Chapman-Enskog expansion object.
    pub fn new(
        knudsen: f64,
        viscosity: f64,
        thermal_conductivity: f64,
        cv: f64,
        temperature: f64,
        density: f64,
    ) -> Self {
        Self {
            knudsen,
            viscosity,
            thermal_conductivity,
            cv,
            temperature,
            density,
        }
    }

    /// First-order (Navier-Stokes) viscous stress tensor component.
    ///
    /// σ_{xx}^(1) = -2μ · ∂u/∂x  (Newtonian constitutive relation).
    ///
    /// # Arguments
    /// * `du_dx` – velocity gradient ∂u/∂x (s⁻¹)
    pub fn first_order_stress(&self, du_dx: f64) -> f64 {
        -2.0 * self.viscosity * du_dx
    }

    /// Second-order (Burnett) correction to the stress tensor.
    ///
    /// σ_{xx}^(2) ≈ Kn² · μ² / ρ · ∂²u/∂x²
    ///
    /// # Arguments
    /// * `d2u_dx2` – second velocity gradient ∂²u/∂x² (m⁻¹ s⁻¹)
    pub fn second_order_stress(&self, d2u_dx2: f64) -> f64 {
        if self.density.abs() < 1e-300 {
            return 0.0;
        }
        self.knudsen * self.knudsen * self.viscosity * self.viscosity / self.density * d2u_dx2
    }

    /// Fourier heat flux at first order.
    ///
    /// q = -κ · ∂T/∂x
    ///
    /// # Arguments
    /// * `dt_dx` – temperature gradient ∂T/∂x (K/m)
    pub fn heat_flux(&self, dt_dx: f64) -> f64 {
        -self.thermal_conductivity * dt_dx
    }

    /// Navier-Stokes limit: effective bulk viscosity ζ in the NS equations.
    ///
    /// For a monatomic ideal gas ζ = 0; for polyatomic gases ζ ≠ 0.
    /// Returns μ (shear viscosity) as the leading-order effective viscosity.
    pub fn navier_stokes_limit(&self) -> f64 {
        self.viscosity
    }

    /// Burnett coefficients (ω₁, ω₂) for the second-order stress corrections.
    ///
    /// For a hard-sphere gas the Burnett coefficients are approximately:
    /// ω₁ = 3/2 and ω₂ = 1 (Kogan 1969).
    pub fn burnett_coefficients(&self) -> (f64, f64) {
        // Standard Burnett coefficients for Maxwell molecules
        let omega1 = 3.0 / 2.0;
        let omega2 = 1.0 + self.knudsen * 0.1; // perturbative Kn correction
        (omega1, omega2)
    }
}

// ============================================================================
// multiscale_diffusion_coefficient
// ============================================================================

/// Combined multiscale diffusion coefficient blending LBM and molecular values.
///
/// Uses a geometric interpolation:
///
/// D_eff = D_lbm^(1-α) · D_molecular^α
///
/// where α = `scale_factor` ∈ \[0, 1\] weights the molecular contribution.
///
/// # Arguments
/// * `lbm_d`        – LBM diffusion coefficient D_lbm (m²/s)
/// * `molecular_d`  – Molecular diffusion coefficient D_mol (m²/s)
/// * `scale_factor` – Blending weight α ∈ \[0, 1\] (0 = pure LBM, 1 = pure MD)
pub fn multiscale_diffusion_coefficient(lbm_d: f64, molecular_d: f64, scale_factor: f64) -> f64 {
    let a = clamp(scale_factor, 0.0, 1.0);
    if lbm_d <= 0.0 || molecular_d <= 0.0 {
        return 0.0;
    }
    lbm_d.powf(1.0 - a) * molecular_d.powf(a)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ---- LbmMolecularDynamicsCoupling ----------------------------------------

    #[test]
    fn test_lbm_md_density_projection_positive() {
        let mut c = LbmMolecularDynamicsCoupling::new(10, 5, 0.1, 0.15);
        // Place particles in the domain
        for i in 0..5 {
            c.md_positions[i] = [i as f64 * 0.1, 0.0, 0.0];
        }
        let rho = c.density_projection();
        let total: f64 = rho.iter().sum();
        assert!(total > 0.0);
    }

    #[test]
    fn test_lbm_md_velocity_projection_finite() {
        let mut c = LbmMolecularDynamicsCoupling::new(8, 3, 0.1, 0.2);
        c.md_positions[0] = [0.3, 0.0, 0.0];
        c.md_velocities[0] = [1.0, 0.0, 0.0];
        let vel = c.velocity_projection();
        for v in vel {
            assert!(v.is_finite());
        }
    }

    #[test]
    fn test_lbm_md_irving_kirkwood_length() {
        let c = LbmMolecularDynamicsCoupling::new(12, 4, 0.05, 0.1);
        let flux = c.irving_kirkwood_flux();
        assert_eq!(flux.len(), 12);
    }

    #[test]
    fn test_lbm_md_force_decomposition_non_negative() {
        let mut c = LbmMolecularDynamicsCoupling::new(4, 3, 0.1, 0.2);
        c.md_forces[0] = [1.0, 2.0, -3.0];
        let fd = c.force_decomposition();
        for f in fd {
            assert!(f >= 0.0);
        }
    }

    // ---- ChapmanEnskogExpansion ---------------------------------------------

    #[test]
    fn test_ce_first_order_stress_sign() {
        let ce = ChapmanEnskogExpansion::new(0.01, 1.8e-5, 0.025, 718.0, 300.0, 1.2);
        // Positive du_dx → negative stress (compression → positive, extension → negative)
        let stress = ce.first_order_stress(1000.0);
        assert!(stress < 0.0);
    }

    #[test]
    fn test_ce_first_order_stress_zero_gradient() {
        let ce = ChapmanEnskogExpansion::new(0.01, 1.8e-5, 0.025, 718.0, 300.0, 1.2);
        assert_eq!(ce.first_order_stress(0.0), 0.0);
    }

    #[test]
    fn test_ce_second_order_stress_scales_with_kn_squared() {
        let ce1 = ChapmanEnskogExpansion::new(0.1, 1e-3, 0.025, 718.0, 300.0, 1.2);
        let ce2 = ChapmanEnskogExpansion::new(0.2, 1e-3, 0.025, 718.0, 300.0, 1.2);
        let s1 = ce1.second_order_stress(1.0);
        let s2 = ce2.second_order_stress(1.0);
        // Kn doubled → stress quadrupled
        assert!((s2 / s1 - 4.0).abs() < 1e-6);
    }

    #[test]
    fn test_ce_heat_flux_fourier_law() {
        let ce = ChapmanEnskogExpansion::new(0.01, 1.8e-5, 0.025, 718.0, 300.0, 1.2);
        let q = ce.heat_flux(-10.0); // negative gradient → positive flux
        assert!(
            q > 0.0,
            "Positive temperature gradient direction → positive heat flux"
        );
    }

    #[test]
    fn test_ce_navier_stokes_limit_equals_viscosity() {
        let mu = 2.5e-4;
        let ce = ChapmanEnskogExpansion::new(0.001, mu, 0.025, 718.0, 300.0, 1.2);
        assert!((ce.navier_stokes_limit() - mu).abs() < 1e-15);
    }

    #[test]
    fn test_ce_burnett_coefficients_positive() {
        let ce = ChapmanEnskogExpansion::new(0.05, 1e-5, 0.025, 718.0, 300.0, 1.0);
        let (w1, w2) = ce.burnett_coefficients();
        assert!(w1 > 0.0 && w2 > 0.0);
    }

    // ---- multiscale_diffusion_coefficient -----------------------------------

    #[test]
    fn test_mdc_pure_lbm() {
        let d = multiscale_diffusion_coefficient(1e-5, 1e-9, 0.0);
        assert!((d - 1e-5).abs() < 1e-20);
    }

    #[test]
    fn test_mdc_pure_molecular() {
        let d = multiscale_diffusion_coefficient(1e-5, 1e-9, 1.0);
        assert!((d - 1e-9).abs() < 1e-24);
    }

    #[test]
    fn test_mdc_midpoint_geometric_mean() {
        let lbm = 1e-4;
        let mol = 1e-8;
        let d = multiscale_diffusion_coefficient(lbm, mol, 0.5);
        let expected = (lbm * mol).sqrt();
        assert!((d - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_mdc_zero_lbm_returns_zero() {
        let d = multiscale_diffusion_coefficient(0.0, 1e-9, 0.5);
        assert_eq!(d, 0.0);
    }

    #[test]
    fn test_mdc_between_limits() {
        let lbm = 1e-4;
        let mol = 1e-8;
        let d = multiscale_diffusion_coefficient(lbm, mol, 0.3);
        assert!(d <= lbm && d >= mol);
    }
}
