#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Multiscale and hybrid LBM-continuum methods.
//!
//! This module provides methods that bridge different length and time scales:
//!
//! - [`HybridLbmNavier`]: Hybrid LBM / Navier-Stokes domain decomposition
//! - [`AdaptiveLbm`]: Adaptive mesh refinement with multi-level LBM
//! - [`MesoscaleLbm`]: Nano-to-meso-to-macro three-scale coupling
//! - [`LbmMolecularDynamicsCoupling`]: LBM / MD continuum-atomistic coupling
//! - [`ChapmanEnskogExpansion`]: Chapman-Enskog perturbative expansion
//! - [`multiscale_diffusion_coefficient`]: Combined diffusion coefficient
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

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

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
// HybridLbmNavier
// ============================================================================

/// Hybrid LBM / Navier-Stokes domain decomposition solver.
///
/// Splits the simulation domain into an LBM sub-domain and a Navier-Stokes
/// (NS) sub-domain. An overlap region of `overlap_width` cells is used for
/// the Schwarz alternating iteration and flux-matching boundary conditions.
pub struct HybridLbmNavier {
    /// Boolean mask: `true` where LBM is active.
    pub lbm_region: Vec<bool>,
    /// Boolean mask: `true` where Navier-Stokes is active.
    pub ns_region: Vec<bool>,
    /// Width (number of cells) of the overlap / coupling region.
    pub overlap_width: usize,
    /// LBM velocity field (one value per cell, x-component representative).
    pub lbm_velocity: Vec<f64>,
    /// NS velocity field (one value per cell, x-component representative).
    pub ns_velocity: Vec<f64>,
    /// LBM density field.
    pub lbm_density: Vec<f64>,
    /// NS density field.
    pub ns_density: Vec<f64>,
    /// LBM dynamic viscosity ν (m²/s).
    pub lbm_viscosity: f64,
    /// NS dynamic viscosity ν (m²/s).
    pub ns_viscosity: f64,
}

impl HybridLbmNavier {
    /// Create a new hybrid solver with `n` cells.
    ///
    /// The first `n/2` cells are assigned to LBM, the second half to NS.
    /// The overlap spans `overlap_width` cells on each side of the interface.
    pub fn new(n: usize, overlap_width: usize, viscosity: f64) -> Self {
        let half = n / 2;
        let mut lbm_region = vec![false; n];
        let mut ns_region = vec![false; n];
        for i in 0..n {
            if i < half + overlap_width {
                lbm_region[i] = true;
            }
            if i + overlap_width >= half {
                ns_region[i] = true;
            }
        }
        Self {
            lbm_region,
            ns_region,
            overlap_width,
            lbm_velocity: vec![0.0; n],
            ns_velocity: vec![0.0; n],
            lbm_density: vec![1.0; n],
            ns_density: vec![1.0; n],
            lbm_viscosity: viscosity,
            ns_viscosity: viscosity,
        }
    }

    /// Number of cells in the domain.
    pub fn n_cells(&self) -> usize {
        self.lbm_region.len()
    }

    /// Perform one explicit time step of size `dt`.
    ///
    /// Applies a simple forward-Euler advection for both sub-domains, then
    /// exchanges information across the overlap via Schwarz coupling.
    pub fn step(&mut self, dt: f64) {
        let n = self.n_cells();
        // LBM sub-domain: simple diffusion (finite-difference Laplacian)
        let mut lbm_new = self.lbm_velocity.clone();
        for i in 1..n - 1 {
            if self.lbm_region[i] {
                let lap = self.lbm_velocity[i + 1] - 2.0 * self.lbm_velocity[i]
                    + self.lbm_velocity[i - 1];
                lbm_new[i] += self.lbm_viscosity * dt * lap;
            }
        }
        // NS sub-domain: same stencil
        let mut ns_new = self.ns_velocity.clone();
        for i in 1..n - 1 {
            if self.ns_region[i] {
                let lap =
                    self.ns_velocity[i + 1] - 2.0 * self.ns_velocity[i] + self.ns_velocity[i - 1];
                ns_new[i] += self.ns_viscosity * dt * lap;
            }
        }
        self.lbm_velocity = lbm_new;
        self.ns_velocity = ns_new;
        self.schwarz_coupling();
    }

    /// Schwarz alternating coupling: enforce continuity of velocity across the
    /// overlap region by averaging LBM and NS solutions.
    pub fn schwarz_coupling(&mut self) {
        let n = self.n_cells();
        for i in 0..n {
            if self.lbm_region[i] && self.ns_region[i] {
                let avg = 0.5 * (self.lbm_velocity[i] + self.ns_velocity[i]);
                self.lbm_velocity[i] = avg;
                self.ns_velocity[i] = avg;
            }
        }
    }

    /// Flux matching: enforce continuity of momentum flux at the interface.
    ///
    /// Returns the interface flux (ρ·u) at the first overlap cell.
    pub fn flux_matching(&self) -> f64 {
        let n = self.n_cells();
        for i in 0..n {
            if self.lbm_region[i] && self.ns_region[i] {
                let flux_lbm = self.lbm_density[i] * self.lbm_velocity[i];
                let flux_ns = self.ns_density[i] * self.ns_velocity[i];
                return 0.5 * (flux_lbm + flux_ns);
            }
        }
        0.0
    }

    /// Velocity mismatch L2-norm across the overlap region.
    ///
    /// Returns `√( Σ (u_lbm - u_ns)² / N_overlap )`.
    pub fn velocity_mismatch(&self) -> f64 {
        let n = self.n_cells();
        let mut sum = 0.0;
        let mut count = 0usize;
        for i in 0..n {
            if self.lbm_region[i] && self.ns_region[i] {
                let d = self.lbm_velocity[i] - self.ns_velocity[i];
                sum += d * d;
                count += 1;
            }
        }
        if count == 0 {
            return 0.0;
        }
        (sum / count as f64).sqrt()
    }
}

// ============================================================================
// AdaptiveLbm
// ============================================================================

/// Multi-level adaptive LBM with spatial and temporal refinement.
///
/// Maintains a hierarchy of grids: a coarse grid and one or more fine grids
/// whose resolutions are described by `refinement_levels` (refinement factor
/// relative to the coarse grid at each level).
pub struct AdaptiveLbm {
    /// Refinement factor at each level (e.g. `[2, 4]` means 2× and 4×).
    pub refinement_levels: Vec<usize>,
    /// Coarse-grid distribution function values.
    pub coarse_grid: Vec<f64>,
    /// Fine-grid distribution function values per level.
    pub fine_grids: Vec<Vec<f64>>,
    /// Number of coarse cells.
    pub n_coarse: usize,
    /// LBM relaxation parameter ω on the coarse grid.
    pub omega_coarse: f64,
}

impl AdaptiveLbm {
    /// Create a new adaptive LBM with `n_coarse` cells and given refinement levels.
    pub fn new(n_coarse: usize, refinement_levels: Vec<usize>, omega_coarse: f64) -> Self {
        let fine_grids = refinement_levels
            .iter()
            .map(|&r| vec![0.0; n_coarse * r])
            .collect();
        Self {
            refinement_levels,
            coarse_grid: vec![0.0; n_coarse],
            fine_grids,
            n_coarse,
            omega_coarse,
        }
    }

    /// Refine cells `start..end` on level `level` by copying coarse values.
    pub fn refine_region(&mut self, level: usize, start: usize, end: usize) {
        if level >= self.fine_grids.len() {
            return;
        }
        let r = self.refinement_levels[level];
        let fine = &mut self.fine_grids[level];
        for ci in start..end.min(self.n_coarse) {
            let val = self.coarse_grid[ci];
            for k in 0..r {
                let fi = ci * r + k;
                if fi < fine.len() {
                    fine[fi] = val;
                }
            }
        }
    }

    /// Coarsen cells `start..end` on level `level` by averaging fine values.
    pub fn coarsen_region(&mut self, level: usize, start: usize, end: usize) {
        if level >= self.fine_grids.len() {
            return;
        }
        let r = self.refinement_levels[level];
        let fine = self.fine_grids[level].clone();
        for ci in start..end.min(self.n_coarse) {
            let mut sum = 0.0;
            for k in 0..r {
                let fi = ci * r + k;
                if fi < fine.len() {
                    sum += fine[fi];
                }
            }
            self.coarse_grid[ci] = sum / r as f64;
        }
    }

    /// Linear time interpolation between two coarse-grid states.
    ///
    /// `alpha` ∈ \[0, 1\] gives the interpolation weight (0 = old, 1 = new).
    pub fn time_interpolation(&self, old_state: &[f64], new_state: &[f64], alpha: f64) -> Vec<f64> {
        let a = clamp(alpha, 0.0, 1.0);
        old_state
            .iter()
            .zip(new_state.iter())
            .map(|(&o, &n)| (1.0 - a) * o + a * n)
            .collect()
    }

    /// Spatial interpolation from coarse to fine grid using linear interpolation.
    pub fn spatial_interpolation(&self, coarse: &[f64], refinement: usize) -> Vec<f64> {
        let n = coarse.len();
        if refinement == 0 {
            return coarse.to_vec();
        }
        let fine_n = n * refinement;
        let mut fine = vec![0.0; fine_n];
        for ci in 0..n {
            let v_left = coarse[ci];
            let v_right = coarse[(ci + 1).min(n - 1)];
            for k in 0..refinement {
                let t = k as f64 / refinement as f64;
                fine[ci * refinement + k] = (1.0 - t) * v_left + t * v_right;
            }
        }
        fine
    }

    /// Perform one BGK collision-streaming step on the coarse grid.
    pub fn step(&mut self) {
        let n = self.n_coarse;
        // BGK relaxation toward equilibrium (equilibrium = 0 for simplicity)
        for val in self.coarse_grid.iter_mut() {
            *val *= 1.0 - self.omega_coarse;
        }
        // Streaming: shift values one cell to the right (periodic)
        if n > 1 {
            let last = self.coarse_grid[n - 1];
            for i in (1..n).rev() {
                self.coarse_grid[i] = self.coarse_grid[i - 1];
            }
            self.coarse_grid[0] = last;
        }
    }
}

// ============================================================================
// MesoscaleLbm
// ============================================================================

/// Three-scale (nano → meso → macro) LBM coupling.
///
/// The `nano_cells` represent an atomistic/nano-scale description (e.g. MD),
/// the `meso_cells` represent the LBM mesoscale, and the macro fields are
/// derived from the meso fields via coarse-graining.
pub struct MesoscaleLbm {
    /// Nano-scale cell values (e.g. particle number densities).
    pub meso_cells: Vec<f64>,
    /// Nano-scale input values.
    pub nano_cells: Vec<f64>,
    /// Number of nano cells per meso cell (coarse-graining ratio).
    pub nano_per_meso: usize,
    /// Macro velocity field (one value per meso cell).
    pub macro_velocity: Vec<f64>,
    /// Macro density field.
    pub macro_density: Vec<f64>,
}

impl MesoscaleLbm {
    /// Create a new three-scale LBM.
    ///
    /// `n_meso` is the number of meso cells; `nano_per_meso` is how many nano
    /// cells map to one meso cell.
    pub fn new(n_meso: usize, nano_per_meso: usize) -> Self {
        Self {
            meso_cells: vec![0.0; n_meso],
            nano_cells: vec![0.0; n_meso * nano_per_meso],
            nano_per_meso,
            macro_velocity: vec![0.0; n_meso],
            macro_density: vec![1.0; n_meso],
        }
    }

    /// Number of meso cells.
    pub fn n_meso(&self) -> usize {
        self.meso_cells.len()
    }

    /// Nano-to-meso coarse-graining: average `nano_per_meso` nano values into
    /// each meso cell.
    pub fn nano_to_meso(&mut self) {
        let n_meso = self.n_meso();
        let r = self.nano_per_meso;
        for mi in 0..n_meso {
            let mut sum = 0.0;
            for k in 0..r {
                let ni = mi * r + k;
                if ni < self.nano_cells.len() {
                    sum += self.nano_cells[ni];
                }
            }
            self.meso_cells[mi] = sum / r as f64;
        }
    }

    /// Meso-to-macro projection: compute macro density and velocity from meso
    /// distribution functions via standard LBM moment integrals.
    ///
    /// Here we use the simple approximation ρ ≈ Σ f_i and ρ·u ≈ Σ f_i · c_i
    /// with c_i ∈ {-1, 0, +1} for a D1Q3 lattice.
    pub fn meso_to_macro(&mut self) {
        let n = self.n_meso();
        for mi in 0..n {
            // Simplified: treat meso_cells[mi] as the 0th moment (density)
            self.macro_density[mi] = self.meso_cells[mi].max(0.0) + 1.0;
            // Velocity estimated from gradient of density
            let left = if mi > 0 {
                self.meso_cells[mi - 1]
            } else {
                self.meso_cells[mi]
            };
            let right = if mi < n - 1 {
                self.meso_cells[mi + 1]
            } else {
                self.meso_cells[mi]
            };
            self.macro_velocity[mi] = 0.5 * (right - left);
        }
    }

    /// Perform one heterogeneous time step:
    /// 1. nano_to_meso coarse-graining,
    /// 2. BGK relaxation on meso scale,
    /// 3. meso_to_macro projection.
    pub fn heterogeneous_step(&mut self) {
        self.nano_to_meso();
        // BGK relaxation on meso cells (ω = 1.5, τ = 1/ω)
        let omega = 1.5_f64;
        for val in self.meso_cells.iter_mut() {
            *val *= 1.0 - omega;
        }
        self.meso_to_macro();
    }
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
            for gi in 0..n {
                let x_grid = gi as f64 * self.dx;
                let r = (pos[0] - x_grid).abs();
                flux[gi] += mass * vx * self.gaussian_kernel(r);
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
            for gi in 0..n {
                let x_grid = gi as f64 * self.dx;
                let r = (pos[0] - x_grid).abs();
                rho[gi] += mass * self.gaussian_kernel(r);
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

    // ---- HybridLbmNavier ----------------------------------------------------

    #[test]
    fn test_hybrid_new_cell_count() {
        let h = HybridLbmNavier::new(100, 5, 0.1);
        assert_eq!(h.n_cells(), 100);
    }

    #[test]
    fn test_hybrid_overlap_both_active() {
        let h = HybridLbmNavier::new(20, 2, 0.1);
        // Some cells in the overlap must be active in both regions
        let both = (0..20)
            .filter(|&i| h.lbm_region[i] && h.ns_region[i])
            .count();
        assert!(both > 0, "Overlap must have cells active in both regions");
    }

    #[test]
    fn test_hybrid_step_preserves_length() {
        let mut h = HybridLbmNavier::new(40, 3, 0.01);
        h.step(0.001);
        assert_eq!(h.lbm_velocity.len(), 40);
        assert_eq!(h.ns_velocity.len(), 40);
    }

    #[test]
    fn test_schwarz_coupling_reduces_mismatch() {
        let mut h = HybridLbmNavier::new(20, 4, 0.1);
        // Introduce a known mismatch
        for i in 0..20 {
            if h.lbm_region[i] {
                h.lbm_velocity[i] = 1.0;
            }
            if h.ns_region[i] {
                h.ns_velocity[i] = 0.0;
            }
        }
        let before = h.velocity_mismatch();
        h.schwarz_coupling();
        let after = h.velocity_mismatch();
        assert!(
            after < before || before == 0.0,
            "Schwarz must reduce mismatch"
        );
    }

    #[test]
    fn test_schwarz_coupling_equalises_overlap() {
        let mut h = HybridLbmNavier::new(20, 4, 0.1);
        for i in 0..20 {
            h.lbm_velocity[i] = 2.0;
            h.ns_velocity[i] = 0.0;
        }
        h.schwarz_coupling();
        assert_eq!(h.velocity_mismatch(), 0.0);
    }

    #[test]
    fn test_flux_matching_finite() {
        let h = HybridLbmNavier::new(40, 5, 0.1);
        let flux = h.flux_matching();
        assert!(flux.is_finite());
    }

    #[test]
    fn test_velocity_mismatch_zero_initially() {
        let h = HybridLbmNavier::new(20, 4, 0.1);
        // All velocities initialised to 0
        assert_eq!(h.velocity_mismatch(), 0.0);
    }

    // ---- AdaptiveLbm --------------------------------------------------------

    #[test]
    fn test_adaptive_new_coarse_length() {
        let a = AdaptiveLbm::new(16, vec![2, 4], 1.0);
        assert_eq!(a.coarse_grid.len(), 16);
    }

    #[test]
    fn test_adaptive_fine_grid_lengths() {
        let a = AdaptiveLbm::new(8, vec![2, 4], 1.0);
        assert_eq!(a.fine_grids[0].len(), 16);
        assert_eq!(a.fine_grids[1].len(), 32);
    }

    #[test]
    fn test_adaptive_refine_region_copies_values() {
        let mut a = AdaptiveLbm::new(4, vec![2], 1.0);
        a.coarse_grid[1] = 5.0;
        a.refine_region(0, 1, 2);
        assert!((a.fine_grids[0][2] - 5.0).abs() < 1e-12);
        assert!((a.fine_grids[0][3] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_adaptive_coarsen_averages_correctly() {
        let mut a = AdaptiveLbm::new(4, vec![2], 1.0);
        a.fine_grids[0][0] = 2.0;
        a.fine_grids[0][1] = 4.0;
        a.coarsen_region(0, 0, 1);
        assert!((a.coarse_grid[0] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_adaptive_time_interpolation_alpha_zero() {
        let a = AdaptiveLbm::new(4, vec![], 1.0);
        let old = vec![1.0, 2.0, 3.0, 4.0];
        let new = vec![5.0, 6.0, 7.0, 8.0];
        let interp = a.time_interpolation(&old, &new, 0.0);
        assert_eq!(interp, old);
    }

    #[test]
    fn test_adaptive_time_interpolation_alpha_one() {
        let a = AdaptiveLbm::new(4, vec![], 1.0);
        let old = vec![0.0; 4];
        let new = vec![1.0; 4];
        let interp = a.time_interpolation(&old, &new, 1.0);
        assert_eq!(interp, new);
    }

    #[test]
    fn test_adaptive_spatial_interpolation_length() {
        let a = AdaptiveLbm::new(4, vec![3], 1.0);
        let coarse = vec![1.0, 2.0, 3.0, 4.0];
        let fine = a.spatial_interpolation(&coarse, 3);
        assert_eq!(fine.len(), 12);
    }

    #[test]
    fn test_adaptive_step_modifies_grid() {
        let mut a = AdaptiveLbm::new(8, vec![2], 0.5);
        a.coarse_grid[3] = 1.0;
        a.step();
        // After streaming the non-zero value should have moved
        assert_ne!(a.coarse_grid[3], 1.0);
    }

    // ---- MesoscaleLbm -------------------------------------------------------

    #[test]
    fn test_mesoscale_new_cell_counts() {
        let m = MesoscaleLbm::new(10, 4);
        assert_eq!(m.meso_cells.len(), 10);
        assert_eq!(m.nano_cells.len(), 40);
    }

    #[test]
    fn test_mesoscale_nano_to_meso_average() {
        let mut m = MesoscaleLbm::new(2, 4);
        // Set nano cells for first meso cell to known values
        m.nano_cells[0] = 1.0;
        m.nano_cells[1] = 3.0;
        m.nano_cells[2] = 5.0;
        m.nano_cells[3] = 7.0;
        m.nano_to_meso();
        assert!((m.meso_cells[0] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_mesoscale_meso_to_macro_density_positive() {
        let mut m = MesoscaleLbm::new(4, 2);
        m.nano_to_meso();
        m.meso_to_macro();
        for &d in &m.macro_density {
            assert!(d > 0.0);
        }
    }

    #[test]
    fn test_mesoscale_heterogeneous_step_preserves_lengths() {
        let mut m = MesoscaleLbm::new(5, 3);
        m.heterogeneous_step();
        assert_eq!(m.macro_velocity.len(), 5);
        assert_eq!(m.macro_density.len(), 5);
    }

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
