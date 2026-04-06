#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Poromechanics and Biot's consolidation theory.
//!
//! Implements coupled solid-fluid mechanics for saturated porous media:
//!
//! - [`BiotCoeff`]: Biot effective stress coefficient and modulus
//! - [`PorousElement`]: Finite element for porous media with coupling matrix
//! - [`ConsolidationSolver`]: Time-stepping solver for Biot consolidation
//! - [`TerzaghiConsolidation`]: 1-D Terzaghi consolidation (analytical)
//! - [`DarcyFlow`]: Darcy seepage velocity computation
//! - [`terzaghi_u`]: Settlement degree of consolidation series expansion
//! - [`consolidation_time`]: Real time from dimensionless time factor
//! - [`effective_stress`]: Terzaghi effective stress principle

#![allow(dead_code)]

use std::f64::consts::PI;

// ============================================================================
// BiotCoeff
// ============================================================================

/// Biot effective-stress coefficient and Biot modulus for a saturated porous medium.
///
/// The Biot coefficient α relates total stress changes to pore-pressure changes:
/// σ′ = σ - α p
///
/// The Biot modulus M controls the compressibility of the fluid-solid mixture.
pub struct BiotCoeff {
    /// Biot effective-stress coefficient α (dimensionless, 0 < α ≤ 1).
    pub alpha: f64,
    /// Biot modulus M (Pa) – inverse of the storage coefficient.
    pub modulus: f64,
}

impl BiotCoeff {
    /// Create a `BiotCoeff` with the given α and M values.
    pub fn new(alpha: f64, modulus: f64) -> Self {
        Self { alpha, modulus }
    }

    /// Compute Biot parameters from drained bulk modulus, fluid bulk modulus and porosity.
    ///
    /// Uses the Biot-Gassmann relations:
    /// α = 1 - K_dry / K_s
    /// 1/M = φ/K_f + (α - φ)/K_s
    ///
    /// # Arguments
    /// * `k_dry` – drained (skeleton) bulk modulus \[Pa\]
    /// * `k_s`   – solid grain bulk modulus \[Pa\]
    /// * `k_f`   – pore fluid bulk modulus \[Pa\]
    /// * `phi`   – porosity (0 < φ < 1)
    pub fn compute_from_bulk(k_dry: f64, k_s: f64, k_f: f64, phi: f64) -> Self {
        let alpha = 1.0 - k_dry / k_s;
        let inv_m = phi / k_f + (alpha - phi) / k_s;
        let modulus = if inv_m.abs() < 1e-300 {
            f64::INFINITY
        } else {
            1.0 / inv_m
        };
        Self { alpha, modulus }
    }

    /// Storage coefficient S = α² / M + φ / K_f (Pa⁻¹).
    pub fn storage_coefficient(&self, phi: f64, k_f: f64) -> f64 {
        self.alpha * self.alpha / self.modulus + phi / k_f
    }

    /// Undrained bulk modulus K_u = K_dry + α² M.
    pub fn undrained_bulk(k_dry: f64, alpha: f64, modulus: f64) -> f64 {
        k_dry + alpha * alpha * modulus
    }

    /// Skempton's B coefficient B = α M / K_u.
    pub fn skempton_b(k_dry: f64, alpha: f64, modulus: f64) -> f64 {
        let k_u = Self::undrained_bulk(k_dry, alpha, modulus);
        if k_u.abs() < 1e-300 {
            0.0
        } else {
            alpha * modulus / k_u
        }
    }
}

// ============================================================================
// PorousElement
// ============================================================================

/// Finite element representing a porous solid-fluid element.
///
/// Stores connectivity, permeability tensor, porosity, and solid stiffness
/// for assembly of the coupled u-p system.
pub struct PorousElement {
    /// Global node IDs for this element.
    pub node_ids: Vec<usize>,
    /// Permeability tensor k \[m²\] stored in row-major 3×3 order.
    pub permeability: [f64; 9],
    /// Porosity φ (dimensionless).
    pub porosity: f64,
    /// Drained Young's modulus E \[Pa\].
    pub solid_stiffness: f64,
}

impl PorousElement {
    /// Create a new `PorousElement`.
    ///
    /// # Arguments
    /// * `node_ids`       – global DOF node indices
    /// * `permeability`   – 3×3 permeability tensor in row-major order \[m²\]
    /// * `porosity`       – porosity φ
    /// * `solid_stiffness`– drained Young's modulus E \[Pa\]
    pub fn new(
        node_ids: Vec<usize>,
        permeability: [f64; 9],
        porosity: f64,
        solid_stiffness: f64,
    ) -> Self {
        Self {
            node_ids,
            permeability,
            porosity,
            solid_stiffness,
        }
    }

    /// Compute a simplified coupling matrix L = α B^T m N_p dV.
    ///
    /// Returns a flat vector of length `3 * n_nodes` representing coupling
    /// contributions for a single Gauss point with unit volume and Biot α = 1.
    pub fn compute_coupling_matrix(&self) -> Vec<f64> {
        let n = self.node_ids.len();
        // Each displacement node has 3 DOFs; each pressure node has 1 DOF.
        // Simplified: return B^T m shape for unit Jacobian, α=1.
        let mut coupling = vec![0.0_f64; 3 * n];
        let n_f64 = n as f64;
        for i in 0..n {
            // Uniform shape function gradient contribution (1/n per direction)
            coupling[3 * i] = 1.0 / n_f64;
            coupling[3 * i + 1] = 1.0 / n_f64;
            coupling[3 * i + 2] = 1.0 / n_f64;
        }
        coupling
    }

    /// Permeability in the x-direction k_xx \[m²\].
    pub fn k_xx(&self) -> f64 {
        self.permeability[0]
    }

    /// Isotropic permeability scalar (mean of diagonal entries).
    pub fn isotropic_permeability(&self) -> f64 {
        (self.permeability[0] + self.permeability[4] + self.permeability[8]) / 3.0
    }
}

// ============================================================================
// ConsolidationSolver
// ============================================================================

/// Time-domain solver for Biot's consolidation equations.
///
/// Solves the coupled u-p system:
/// \[ K_uu  K_up \] { u̇ }   \[ 0  0   \] { u }   { f_u }
/// \[ K_up^T  0  \] { ṗ } + \[ 0  H_pp\] { p } = { f_p }
///
/// using a first-order implicit (backward Euler) time integration.
pub struct ConsolidationSolver {
    /// Displacement DOF vector.
    pub u_dofs: Vec<f64>,
    /// Pore-pressure DOF vector.
    pub p_dofs: Vec<f64>,
    /// Mechanical stiffness matrix K_uu (row-major dense).
    pub k_uu: Vec<f64>,
    /// Coupling matrix K_up (n_u × n_p, row-major dense).
    pub k_up: Vec<f64>,
    /// Compressibility matrix K_pp (diagonal, stored as vector).
    pub k_pp: Vec<f64>,
    /// Permeability / flow matrix H_pp (diagonal, stored as vector).
    pub h_pp: Vec<f64>,
    /// Current simulation time \[s\].
    pub time: f64,
}

impl ConsolidationSolver {
    /// Create a new `ConsolidationSolver`.
    ///
    /// # Arguments
    /// * `n_u`  – number of displacement DOFs
    /// * `n_p`  – number of pressure DOFs
    pub fn new(n_u: usize, n_p: usize) -> Self {
        Self {
            u_dofs: vec![0.0; n_u],
            p_dofs: vec![0.0; n_p],
            k_uu: vec![0.0; n_u * n_u],
            k_up: vec![0.0; n_u * n_p],
            k_pp: vec![0.0; n_p],
            h_pp: vec![0.0; n_p],
            time: 0.0,
        }
    }

    /// Perform one backward-Euler consolidation time step of size `dt` \[s\].
    ///
    /// For this simplified implementation the pressure update uses the diagonal
    /// compressibility / permeability matrices only (no off-diagonal coupling),
    /// which is exact for single-element test problems and provides a skeleton
    /// for full assembly in production use.
    pub fn step(&mut self, dt: f64) {
        let n_p = self.p_dofs.len();
        for i in 0..n_p {
            // Implicit update: (K_pp + dt * H_pp) p_{n+1} = K_pp * p_n
            let lhs = self.k_pp[i] + dt * self.h_pp[i];
            if lhs.abs() > 1e-300 {
                self.p_dofs[i] = self.k_pp[i] * self.p_dofs[i] / lhs;
            }
        }
        self.time += dt;
    }

    /// Number of displacement DOFs.
    pub fn n_u(&self) -> usize {
        self.u_dofs.len()
    }

    /// Number of pressure DOFs.
    pub fn n_p(&self) -> usize {
        self.p_dofs.len()
    }
}

// ============================================================================
// TerzaghiConsolidation
// ============================================================================

/// Analytical 1-D Terzaghi consolidation model.
///
/// Describes the time-dependent settlement of a saturated clay layer of
/// thickness H drained at the top (single drainage) or both surfaces (double
/// drainage) under a suddenly applied load.
pub struct TerzaghiConsolidation {
    /// Drainage height H \[m\] (half-thickness for double drainage).
    pub drainage_height: f64,
    /// Coefficient of consolidation c_v \[m²/s\].
    pub cv: f64,
}

impl TerzaghiConsolidation {
    /// Create a new `TerzaghiConsolidation`.
    ///
    /// # Arguments
    /// * `drainage_height` – drainage path length H \[m\]
    /// * `cv`              – coefficient of consolidation c_v \[m²/s\]
    pub fn new(drainage_height: f64, cv: f64) -> Self {
        Self {
            drainage_height,
            cv,
        }
    }

    /// Dimensionless time factor T_v = c_v t / H².
    pub fn time_factor(&self, time: f64) -> f64 {
        self.cv * time / (self.drainage_height * self.drainage_height)
    }

    /// Average degree of consolidation U (0–1) using Terzaghi series.
    ///
    /// Uses `terzaghi_u` with T_v = c_v t / H².
    pub fn degree_of_consolidation(&self, time: f64) -> f64 {
        let tv = self.time_factor(time);
        terzaghi_u(tv)
    }

    /// Settlement at time `time` given final settlement `s_final` \[m\].
    pub fn settlement(&self, time: f64, s_final: f64) -> f64 {
        self.degree_of_consolidation(time) * s_final
    }

    /// Time \[s\] required to reach degree of consolidation `u_target` (0–1).
    ///
    /// Uses Casagrande's approximate formula.
    pub fn time_to_consolidation(&self, u_target: f64) -> f64 {
        // Approximate: T_v ≈ π/4 * U² for U ≤ 0.6, else use log approximation
        let tv = if u_target <= 0.6 {
            PI / 4.0 * u_target * u_target
        } else {
            -0.9332 * (1.0 - u_target).ln() - 0.0851
        };
        tv * self.drainage_height * self.drainage_height / self.cv
    }
}

// ============================================================================
// DarcyFlow
// ============================================================================

/// Darcy seepage flow model.
///
/// Computes the specific discharge (Darcy velocity) q = -k/μ * ∇h
/// for laminar flow through a porous medium.
pub struct DarcyFlow {
    /// Intrinsic permeability k \[m²\].
    pub permeability: f64,
    /// Dynamic viscosity μ \[Pa·s\].
    pub viscosity: f64,
    /// Hydraulic gradient ∇h \[m/m\] (dimensionless).
    pub hydraulic_gradient: [f64; 3],
}

impl DarcyFlow {
    /// Create a new `DarcyFlow`.
    pub fn new(permeability: f64, viscosity: f64, hydraulic_gradient: [f64; 3]) -> Self {
        Self {
            permeability,
            viscosity,
            hydraulic_gradient,
        }
    }

    /// Compute the Darcy velocity q = -k/μ * ∇h \[m/s\].
    pub fn darcy_velocity(&self) -> [f64; 3] {
        let factor = -self.permeability / self.viscosity;
        [
            factor * self.hydraulic_gradient[0],
            factor * self.hydraulic_gradient[1],
            factor * self.hydraulic_gradient[2],
        ]
    }

    /// Magnitude of the Darcy velocity |q| \[m/s\].
    pub fn darcy_speed(&self) -> f64 {
        let q = self.darcy_velocity();
        (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt()
    }

    /// Hydraulic conductivity K = k ρ g / μ \[m/s\].
    ///
    /// # Arguments
    /// * `density` – fluid density ρ \[kg/m³\]
    /// * `g`       – gravitational acceleration \[m/s²\]
    pub fn hydraulic_conductivity(&self, density: f64, g: f64) -> f64 {
        self.permeability * density * g / self.viscosity
    }

    /// Reynolds number for Darcy flow Re = ρ |q| d / μ.
    ///
    /// # Arguments
    /// * `density`          – fluid density ρ \[kg/m³\]
    /// * `particle_diameter`– mean grain diameter d \[m\]
    pub fn reynolds_number(&self, density: f64, particle_diameter: f64) -> f64 {
        density * self.darcy_speed() * particle_diameter / self.viscosity
    }
}

// ============================================================================
// Free functions
// ============================================================================

/// Average degree of consolidation U using Terzaghi Fourier series.
///
/// U(T_v) = 1 - Σ_{m=0}^{∞} (8 / (2m+1)² π²) exp(-(2m+1)² π² T_v / 4)
///
/// The series is truncated after 50 terms for practical accuracy.
pub fn terzaghi_u(tv: f64) -> f64 {
    if tv <= 0.0 {
        return 0.0;
    }
    let mut sum = 0.0_f64;
    for m in 0_u32..50 {
        let two_m1 = (2 * m + 1) as f64;
        let coeff = 8.0 / (two_m1 * two_m1 * PI * PI);
        let exponent = -two_m1 * two_m1 * PI * PI * tv / 4.0;
        sum += coeff * exponent.exp();
    }
    (1.0 - sum).clamp(0.0, 1.0)
}

/// Real time t \[s\] from dimensionless time factor T_v, c_v \[m²/s\] and H \[m\].
///
/// t = T_v * H² / c_v
pub fn consolidation_time(tv: f64, cv: f64, h: f64) -> f64 {
    tv * h * h / cv
}

/// Effective stress σ′ = σ - p (Terzaghi effective stress principle) \[Pa\].
///
/// # Arguments
/// * `total`         – total stress σ \[Pa\]
/// * `pore_pressure` – pore fluid pressure p \[Pa\]
pub fn effective_stress(total: f64, pore_pressure: f64) -> f64 {
    total - pore_pressure
}

/// Consolidation coefficient c_v from permeability, compressibility and fluid properties.
///
/// c_v = k / (μ m_v)  where m_v is the coefficient of volume compressibility.
///
/// # Arguments
/// * `k`   – hydraulic conductivity k \[m/s\] (K = k / (ρ g) style input in m/s)
/// * `mu`  – dynamic viscosity \[Pa·s\]
/// * `m_v` – coefficient of volume compressibility \[1/Pa\]
pub fn consolidation_coefficient(k: f64, mu: f64, m_v: f64) -> f64 {
    if mu.abs() < 1e-300 || m_v.abs() < 1e-300 {
        return 0.0;
    }
    k / (mu * m_v)
}

/// Pore pressure dissipation ratio at depth z in a layer of height H at time factor T_v.
///
/// u(z, T_v) / u_0 = Σ_{m=0}^{∞} (4 / (2m+1) π) sin((2m+1) π z / (2H)) exp(...)
pub fn pore_pressure_ratio(z: f64, h: f64, tv: f64) -> f64 {
    if h < 1e-300 {
        return 0.0;
    }
    let mut sum = 0.0_f64;
    for m in 0_u32..50 {
        let two_m1 = (2 * m + 1) as f64;
        let coeff = 4.0 / (two_m1 * PI);
        let spatial = (two_m1 * PI * z / (2.0 * h)).sin();
        let temporal = (-two_m1 * two_m1 * PI * PI * tv / 4.0).exp();
        sum += coeff * spatial * temporal;
    }
    sum.clamp(0.0, 1.0)
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // --- BiotCoeff ---

    #[test]
    fn test_biot_new() {
        let b = BiotCoeff::new(0.8, 1e9);
        assert!((b.alpha - 0.8).abs() < 1e-12);
        assert!((b.modulus - 1e9).abs() < 1.0);
    }

    #[test]
    fn test_biot_compute_from_bulk_incompressible_solid() {
        // If K_dry ≈ 0 (very soft skeleton) and K_s >> K_f, α → 1
        let b = BiotCoeff::compute_from_bulk(0.0, 1e12, 2.2e9, 0.4);
        assert!((b.alpha - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_biot_compute_from_bulk_typical_sandstone() {
        // Typical sandstone: K_dry=10 GPa, K_s=36 GPa, K_f=2.2 GPa, phi=0.2
        let b = BiotCoeff::compute_from_bulk(10e9, 36e9, 2.2e9, 0.2);
        assert!(b.alpha > 0.0 && b.alpha < 1.0);
        assert!(b.modulus > 0.0);
    }

    #[test]
    fn test_biot_alpha_range() {
        let b = BiotCoeff::compute_from_bulk(5e9, 36e9, 2.2e9, 0.3);
        assert!(b.alpha >= 0.0 && b.alpha <= 1.0);
    }

    #[test]
    fn test_biot_storage_coefficient() {
        let b = BiotCoeff::new(1.0, 1e9);
        let s = b.storage_coefficient(0.4, 2.2e9);
        assert!(s > 0.0);
    }

    #[test]
    fn test_biot_undrained_bulk() {
        let k_u = BiotCoeff::undrained_bulk(10e9, 0.8, 5e9);
        // K_u > K_dry
        assert!(k_u > 10e9);
    }

    #[test]
    fn test_biot_skempton_b() {
        let b_val = BiotCoeff::skempton_b(10e9, 0.8, 5e9);
        assert!(b_val > 0.0 && b_val <= 1.0);
    }

    // --- PorousElement ---

    #[test]
    fn test_porous_element_new() {
        let k = [1e-12, 0.0, 0.0, 0.0, 1e-12, 0.0, 0.0, 0.0, 1e-12];
        let elem = PorousElement::new(vec![0, 1, 2, 3], k, 0.3, 1e7);
        assert_eq!(elem.node_ids.len(), 4);
        assert!((elem.porosity - 0.3).abs() < 1e-12);
    }

    #[test]
    fn test_porous_element_coupling_matrix_length() {
        let k = [1e-12; 9];
        let elem = PorousElement::new(vec![0, 1, 2, 3], k, 0.3, 1e7);
        let coupling = elem.compute_coupling_matrix();
        assert_eq!(coupling.len(), 12); // 3 * 4
    }

    #[test]
    fn test_porous_element_coupling_matrix_values() {
        let k = [1e-12; 9];
        let elem = PorousElement::new(vec![0, 1, 2, 3], k, 0.3, 1e7);
        let coupling = elem.compute_coupling_matrix();
        // Each entry should be 1/4 = 0.25
        for &v in &coupling {
            assert!((v - 0.25).abs() < 1e-12);
        }
    }

    #[test]
    fn test_porous_element_k_xx() {
        let mut k = [0.0_f64; 9];
        k[0] = 2e-12;
        let elem = PorousElement::new(vec![0], k, 0.2, 1e8);
        assert!((elem.k_xx() - 2e-12).abs() < 1e-300);
    }

    #[test]
    fn test_porous_element_isotropic_permeability() {
        let k = [1e-12, 0.0, 0.0, 0.0, 2e-12, 0.0, 0.0, 0.0, 3e-12];
        let elem = PorousElement::new(vec![0], k, 0.2, 1e8);
        assert!((elem.isotropic_permeability() - 2e-12).abs() < 1e-24);
    }

    // --- ConsolidationSolver ---

    #[test]
    fn test_consolidation_solver_new() {
        let solver = ConsolidationSolver::new(6, 2);
        assert_eq!(solver.n_u(), 6);
        assert_eq!(solver.n_p(), 2);
        assert!((solver.time - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_consolidation_solver_step_time() {
        let mut solver = ConsolidationSolver::new(3, 1);
        solver.step(0.1);
        assert!((solver.time - 0.1).abs() < 1e-12);
        solver.step(0.1);
        assert!((solver.time - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_consolidation_solver_pressure_decay() {
        let mut solver = ConsolidationSolver::new(3, 1);
        solver.p_dofs[0] = 100.0;
        solver.k_pp[0] = 1.0;
        solver.h_pp[0] = 1.0;
        // After one step: p = k_pp * p0 / (k_pp + dt * h_pp) = 100/(1+1) = 50
        solver.step(1.0);
        assert!((solver.p_dofs[0] - 50.0).abs() < 1e-10);
    }

    #[test]
    fn test_consolidation_solver_pressure_zero_kpp() {
        let mut solver = ConsolidationSolver::new(3, 1);
        solver.p_dofs[0] = 100.0;
        solver.k_pp[0] = 0.0;
        solver.h_pp[0] = 1.0;
        // lhs = 0 + dt * 1 = dt; k_pp = 0 → p = 0
        solver.step(0.5);
        assert!((solver.p_dofs[0]).abs() < 1e-12);
    }

    // --- TerzaghiConsolidation ---

    #[test]
    fn test_terzaghi_time_factor_zero() {
        let tc = TerzaghiConsolidation::new(1.0, 1e-6);
        assert!((tc.time_factor(0.0)).abs() < 1e-12);
    }

    #[test]
    fn test_terzaghi_time_factor() {
        let cv = 1e-6;
        let h = 2.0;
        let tc = TerzaghiConsolidation::new(h, cv);
        let t = 1e5;
        let tv = cv * t / (h * h);
        assert!((tc.time_factor(t) - tv).abs() < 1e-20);
    }

    #[test]
    fn test_terzaghi_degree_zero_time() {
        let tc = TerzaghiConsolidation::new(1.0, 1e-6);
        assert!((tc.degree_of_consolidation(0.0)).abs() < 1e-10);
    }

    #[test]
    fn test_terzaghi_degree_large_time() {
        let tc = TerzaghiConsolidation::new(1.0, 1e-6);
        // t → ∞ should give U → 1
        let u = tc.degree_of_consolidation(1e20);
        assert!((u - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_terzaghi_degree_monotonic() {
        let tc = TerzaghiConsolidation::new(1.0, 1e-7);
        let u1 = tc.degree_of_consolidation(1e4);
        let u2 = tc.degree_of_consolidation(2e4);
        assert!(u2 >= u1);
    }

    #[test]
    fn test_terzaghi_settlement() {
        let tc = TerzaghiConsolidation::new(1.0, 1e-6);
        let s = tc.settlement(1e20, 0.1);
        assert!((s - 0.1).abs() < 1e-5);
    }

    #[test]
    fn test_terzaghi_time_to_consolidation_50pct() {
        let cv = 1e-6;
        let h = 1.0;
        let tc = TerzaghiConsolidation::new(h, cv);
        let t50 = tc.time_to_consolidation(0.5);
        // Verify degree at computed time is near 0.5
        let u_check = tc.degree_of_consolidation(t50);
        assert!((u_check - 0.5).abs() < 0.05);
    }

    // --- DarcyFlow ---

    #[test]
    fn test_darcy_velocity_1d() {
        let flow = DarcyFlow::new(1e-12, 1e-3, [1.0, 0.0, 0.0]);
        let v = flow.darcy_velocity();
        // q = -k/mu * grad = -1e-12 / 1e-3 = -1e-9 m/s
        assert!((v[0] - (-1e-9)).abs() < 1e-21);
        assert!((v[1]).abs() < 1e-30);
        assert!((v[2]).abs() < 1e-30);
    }

    #[test]
    fn test_darcy_velocity_zero_gradient() {
        let flow = DarcyFlow::new(1e-12, 1e-3, [0.0, 0.0, 0.0]);
        let v = flow.darcy_velocity();
        assert!(v[0].abs() < 1e-30 && v[1].abs() < 1e-30 && v[2].abs() < 1e-30);
    }

    #[test]
    fn test_darcy_speed() {
        let flow = DarcyFlow::new(1e-12, 1e-3, [1.0, 0.0, 0.0]);
        let speed = flow.darcy_speed();
        assert!((speed - 1e-9).abs() < 1e-21);
    }

    #[test]
    fn test_darcy_hydraulic_conductivity() {
        // k=1e-12 m², rho=1000, g=9.81, mu=1e-3
        let flow = DarcyFlow::new(1e-12, 1e-3, [0.0; 3]);
        let k_cond = flow.hydraulic_conductivity(1000.0, 9.81);
        let expected = 1e-12 * 1000.0 * 9.81 / 1e-3;
        assert!((k_cond - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_darcy_reynolds_number() {
        let flow = DarcyFlow::new(1e-12, 1e-3, [1.0, 0.0, 0.0]);
        let re = flow.reynolds_number(1000.0, 1e-3);
        assert!(re >= 0.0);
    }

    // --- Free functions ---

    #[test]
    fn test_terzaghi_u_zero() {
        assert!((terzaghi_u(0.0)).abs() < 1e-12);
    }

    #[test]
    fn test_terzaghi_u_large() {
        assert!((terzaghi_u(100.0) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_terzaghi_u_half() {
        // At T_v ≈ 0.197, U ≈ 0.5
        let u = terzaghi_u(0.197);
        assert!((u - 0.5).abs() < 0.02);
    }

    #[test]
    fn test_terzaghi_u_monotonic() {
        let u1 = terzaghi_u(0.1);
        let u2 = terzaghi_u(0.3);
        let u3 = terzaghi_u(1.0);
        assert!(u1 <= u2 && u2 <= u3);
    }

    #[test]
    fn test_consolidation_time_roundtrip() {
        let cv = 2e-8;
        let h = 3.0;
        let tv = 0.5;
        let t = consolidation_time(tv, cv, h);
        let tv_back = cv * t / (h * h);
        assert!((tv_back - tv).abs() < 1e-12);
    }

    #[test]
    fn test_effective_stress_positive() {
        let sigma_eff = effective_stress(100e3, 30e3);
        assert!((sigma_eff - 70e3).abs() < 1e-8);
    }

    #[test]
    fn test_effective_stress_zero_pore() {
        assert!((effective_stress(50e3, 0.0) - 50e3).abs() < 1e-10);
    }

    #[test]
    fn test_effective_stress_negative() {
        // Tensile effective stress
        let s = effective_stress(10e3, 50e3);
        assert!((s - (-40e3)).abs() < 1e-8);
    }

    #[test]
    fn test_consolidation_coefficient_basic() {
        let cc = consolidation_coefficient(1e-9, 1e-3, 1e-7);
        assert!(cc > 0.0);
    }

    #[test]
    fn test_consolidation_coefficient_zero_viscosity() {
        let cc = consolidation_coefficient(1e-9, 0.0, 1e-7);
        assert!((cc).abs() < 1e-12);
    }

    #[test]
    fn test_pore_pressure_ratio_surface() {
        // At z=0 (drainage surface) ratio → 0
        let r = pore_pressure_ratio(0.0, 1.0, 1.0);
        assert!(r.abs() < 0.01);
    }

    #[test]
    fn test_pore_pressure_ratio_zero_time() {
        // At T_v=0, ratio should be ~1 at mid-layer
        let r = pore_pressure_ratio(1.0, 1.0, 0.0);
        assert!(r > 0.9);
    }

    #[test]
    fn test_pore_pressure_ratio_range() {
        for i in 0..10 {
            let tv = i as f64 * 0.1;
            let r = pore_pressure_ratio(0.5, 1.0, tv);
            assert!(r >= 0.0 && r <= 1.0);
        }
    }
}
