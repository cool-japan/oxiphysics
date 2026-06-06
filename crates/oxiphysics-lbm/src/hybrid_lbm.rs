// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Hybrid LBM–Navier-Stokes coupling.
//!
//! This module provides infrastructure for coupling LBM and continuum NS
//! solvers within a single simulation domain.  It covers:
//!
//! - Region classification ([`HybridRegion`], [`HybridMesh`])
//! - Coupling parameters ([`HybridCouplingParams`])
//! - LBM-to-NS and NS-to-LBM field conversion helpers
//! - Interface interpolation and domain decomposition
//! - Conservation-check utilities

// ============================================================================
// D2Q9 velocity set (shared constant)
// ============================================================================

/// Standard D2Q9 discrete velocities.
pub const D2Q9_VELOCITIES: [[i32; 2]; 9] = [
    [0, 0],   // 0 – rest
    [1, 0],   // 1
    [0, 1],   // 2
    [-1, 0],  // 3
    [0, -1],  // 4
    [1, 1],   // 5
    [-1, 1],  // 6
    [-1, -1], // 7
    [1, -1],  // 8
];

/// D2Q9 equilibrium weights.
pub const D2Q9_WEIGHTS: [f64; 9] = [
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

// ============================================================================
// Structs
// ============================================================================

/// A rectangular region of the hybrid mesh classified as LBM or NS domain.
#[derive(Debug, Clone)]
pub struct HybridRegion {
    /// `true` if this region is solved with LBM, `false` for Navier-Stokes.
    pub is_lbm: bool,
    /// Axis-aligned bounding box `[[xmin, xmax\], [ymin, ymax]]`.
    pub bounds: [[f64; 2]; 2],
}

impl HybridRegion {
    /// Create a new `HybridRegion`.
    pub fn new(is_lbm: bool, bounds: [[f64; 2]; 2]) -> Self {
        Self { is_lbm, bounds }
    }

    /// Check whether a point (x, y) lies inside this region (inclusive).
    pub fn contains(&self, x: f64, y: f64) -> bool {
        x >= self.bounds[0][0]
            && x <= self.bounds[0][1]
            && y >= self.bounds[1][0]
            && y <= self.bounds[1][1]
    }
}

/// A hybrid mesh composed of multiple [`HybridRegion`]s.
#[derive(Debug, Clone)]
pub struct HybridMesh {
    /// Ordered list of regions that partition the domain.
    pub regions: Vec<HybridRegion>,
    /// Nodes on the LBM/NS interface, each given as `[x, y]`.
    pub interface_nodes: Vec<[f64; 2]>,
}

impl HybridMesh {
    /// Create an empty `HybridMesh`.
    pub fn new() -> Self {
        Self {
            regions: Vec::new(),
            interface_nodes: Vec::new(),
        }
    }

    /// Add a region to the mesh.
    pub fn add_region(&mut self, region: HybridRegion) {
        self.regions.push(region);
    }

    /// Detect and record interface nodes between adjacent LBM and NS regions.
    ///
    /// A node in `candidate_nodes` is added to `interface_nodes` if it lies on
    /// the boundary between an LBM region and an NS region.
    pub fn detect_interface(&mut self, candidate_nodes: &[[f64; 2]]) {
        self.interface_nodes.clear();
        for &node in candidate_nodes {
            let in_lbm = self
                .regions
                .iter()
                .any(|r| r.is_lbm && r.contains(node[0], node[1]));
            let in_ns = self
                .regions
                .iter()
                .any(|r| !r.is_lbm && r.contains(node[0], node[1]));
            if in_lbm && in_ns {
                self.interface_nodes.push(node);
            }
        }
    }
}

impl Default for HybridMesh {
    fn default() -> Self {
        Self::new()
    }
}

/// Parameters controlling the LBM–NS overlap coupling.
#[derive(Debug, Clone)]
pub struct HybridCouplingParams {
    /// Width of the overlap (buffer) zone between LBM and NS domains.
    pub overlap_width: f64,
    /// Polynomial order used for interface interpolation (1 = linear).
    pub interpolation_order: usize,
}

impl HybridCouplingParams {
    /// Create coupling parameters.
    pub fn new(overlap_width: f64, interpolation_order: usize) -> Self {
        Self {
            overlap_width,
            interpolation_order,
        }
    }
}

// ============================================================================
// LBM → NS conversions
// ============================================================================

/// Compute the macroscopic velocity \[ux, uy\] from D2Q9 distributions.
///
/// u_α = Σ_i f_i * e_i_α / rho,  where rho = Σ_i f_i.
///
/// Returns `[0.0, 0.0]` when rho ≤ 0.
pub fn lbm_to_ns_velocity(f: &[f64], d2q9_vels: &[[i32; 2]; 9]) -> [f64; 2] {
    let rho: f64 = f.iter().sum();
    if rho <= 0.0 {
        return [0.0, 0.0];
    }
    let mut u = [0.0_f64; 2];
    for (i, &fi) in f.iter().enumerate().take(9) {
        u[0] += fi * d2q9_vels[i][0] as f64;
        u[1] += fi * d2q9_vels[i][1] as f64;
    }
    [u[0] / rho, u[1] / rho]
}

/// Compute macroscopic density ρ = Σ_i f_i from D2Q9 distributions.
pub fn lbm_to_ns_density(f: &[f64]) -> f64 {
    f.iter().sum()
}

/// Reconstruct D2Q9 equilibrium distribution f_eq^α.
///
/// f_eq^α = w_α * ρ * (1 + (e·u)/cs² + (e·u)²/(2cs⁴) − |u|²/(2cs²))
///
/// with cs² = 1/3.
pub fn ns_to_lbm_equilibrium(rho: f64, u: [f64; 2], alpha: usize) -> f64 {
    if alpha >= 9 {
        return 0.0;
    }
    let cs2 = 1.0 / 3.0;
    let w = D2Q9_WEIGHTS[alpha];
    let ex = D2Q9_VELOCITIES[alpha][0] as f64;
    let ey = D2Q9_VELOCITIES[alpha][1] as f64;
    let eu = ex * u[0] + ey * u[1];
    let u2 = u[0] * u[0] + u[1] * u[1];
    w * rho * (1.0 + eu / cs2 + eu * eu / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2))
}

// ============================================================================
// Interface utilities
// ============================================================================

/// Linear interpolation between two values at parameter t ∈ \[0, 1\].
///
/// result = val_a * (1 − t) + val_b * t
pub fn interface_interpolation_linear(val_a: f64, val_b: f64, t: f64) -> f64 {
    val_a * (1.0 - t) + val_b * t
}

/// Split a 1-D domain of `n` cells into an LBM portion and an NS portion.
///
/// Returns `(n_lbm, n_ns)` where `n_lbm` is clamped to `[0, n]` and
/// `n_ns = n − n_lbm`.
pub fn domain_decompose_1d(n: usize, n_lbm: usize) -> (usize, usize) {
    let lbm = n_lbm.min(n);
    (lbm, n - lbm)
}

// ============================================================================
// Conservation checks
// ============================================================================

/// Return `true` when the relative kinetic-energy mismatch is within `tol`.
///
/// |lbm_ke − ns_ke| / max(|lbm_ke|, |ns_ke|, ε) ≤ tol
pub fn energy_conservation_check(lbm_ke: f64, ns_ke: f64, tol: f64) -> bool {
    let denom = lbm_ke.abs().max(ns_ke.abs()).max(1e-30);
    (lbm_ke - ns_ke).abs() / denom <= tol
}

/// Compute the mass flux through a 1-D interface (discrete sum).
///
/// flux = Σ_i ρ_i * (u_i · n) * dx
///
/// `density` and `velocity` are parallel slices of the interface cells.
/// `normal` is the outward unit normal of the interface.
pub fn mass_flux_interface(density: &[f64], velocity: &[f64], normal: [f64; 2], dx: f64) -> f64 {
    density
        .iter()
        .zip(velocity.iter())
        .map(|(&rho, &u)| {
            // project scalar velocity onto normal (treat velocity as x-component)
            rho * u * (normal[0] + normal[1]) * dx
        })
        .sum()
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // -- HybridRegion ---------------------------------------------------------

    #[test]
    fn test_region_contains_inside() {
        let r = HybridRegion::new(true, [[0.0, 1.0], [0.0, 1.0]]);
        assert!(r.contains(0.5, 0.5));
    }

    #[test]
    fn test_region_contains_boundary() {
        let r = HybridRegion::new(true, [[0.0, 1.0], [0.0, 1.0]]);
        assert!(r.contains(0.0, 0.0));
        assert!(r.contains(1.0, 1.0));
    }

    #[test]
    fn test_region_outside() {
        let r = HybridRegion::new(false, [[0.0, 1.0], [0.0, 1.0]]);
        assert!(!r.contains(2.0, 0.5));
    }

    #[test]
    fn test_region_is_lbm_flag() {
        let r = HybridRegion::new(true, [[0.0, 1.0], [0.0, 1.0]]);
        assert!(r.is_lbm);
    }

    // -- HybridMesh -----------------------------------------------------------

    #[test]
    fn test_mesh_new_empty() {
        let m = HybridMesh::new();
        assert!(m.regions.is_empty());
        assert!(m.interface_nodes.is_empty());
    }

    #[test]
    fn test_mesh_add_region() {
        let mut m = HybridMesh::new();
        m.add_region(HybridRegion::new(true, [[0.0, 1.0], [0.0, 1.0]]));
        m.add_region(HybridRegion::new(false, [[1.0, 2.0], [0.0, 1.0]]));
        assert_eq!(m.regions.len(), 2);
    }

    #[test]
    fn test_mesh_detect_interface() {
        let mut m = HybridMesh::new();
        // LBM: [0,1] x [0,1], NS: [0,1] x [0,1] (full overlap)
        m.add_region(HybridRegion::new(true, [[0.0, 1.0], [0.0, 1.0]]));
        m.add_region(HybridRegion::new(false, [[0.0, 1.0], [0.0, 1.0]]));
        let candidates = vec![[0.5, 0.5], [2.0, 2.0]];
        m.detect_interface(&candidates);
        assert_eq!(m.interface_nodes.len(), 1);
        assert_eq!(m.interface_nodes[0], [0.5, 0.5]);
    }

    #[test]
    fn test_mesh_default() {
        let m = HybridMesh::default();
        assert!(m.regions.is_empty());
    }

    // -- HybridCouplingParams -------------------------------------------------

    #[test]
    fn test_coupling_params() {
        let p = HybridCouplingParams::new(0.1, 1);
        assert!((p.overlap_width - 0.1).abs() < 1e-12);
        assert_eq!(p.interpolation_order, 1);
    }

    // -- lbm_to_ns_velocity ---------------------------------------------------

    #[test]
    fn test_velocity_rest() {
        let f = D2Q9_WEIGHTS;
        let u = lbm_to_ns_velocity(&f, &D2Q9_VELOCITIES);
        assert!(u[0].abs() < 1e-12);
        assert!(u[1].abs() < 1e-12);
    }

    #[test]
    fn test_velocity_zero_density() {
        let f = [0.0_f64; 9];
        let u = lbm_to_ns_velocity(&f, &D2Q9_VELOCITIES);
        assert_eq!(u, [0.0, 0.0]);
    }

    #[test]
    fn test_velocity_nonzero() {
        // Shift mass into direction 1 (+x)
        let mut f = [0.0_f64; 9];
        f[0] = 0.5;
        f[1] = 0.5;
        let u = lbm_to_ns_velocity(&f, &D2Q9_VELOCITIES);
        assert!(u[0] > 0.0);
    }

    // -- lbm_to_ns_density ----------------------------------------------------

    #[test]
    fn test_density_rest() {
        let rho = lbm_to_ns_density(&D2Q9_WEIGHTS);
        assert!((rho - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_density_zero() {
        assert_eq!(lbm_to_ns_density(&[0.0_f64; 9]), 0.0);
    }

    // -- ns_to_lbm_equilibrium ------------------------------------------------

    #[test]
    fn test_equilibrium_sum_to_rho() {
        let rho = 1.2;
        let u = [0.1, 0.05];
        let sum: f64 = (0..9).map(|a| ns_to_lbm_equilibrium(rho, u, a)).sum();
        assert!((sum - rho).abs() < 1e-10, "sum={sum}, rho={rho}");
    }

    #[test]
    fn test_equilibrium_rest_matches_weights() {
        let rho = 1.0;
        let u = [0.0, 0.0];
        for (a, &w) in D2Q9_WEIGHTS.iter().enumerate() {
            let feq = ns_to_lbm_equilibrium(rho, u, a);
            assert!((feq - w).abs() < 1e-12, "a={a}");
        }
    }

    #[test]
    fn test_equilibrium_out_of_range() {
        assert_eq!(ns_to_lbm_equilibrium(1.0, [0.0, 0.0], 9), 0.0);
    }

    // -- interface_interpolation_linear ---------------------------------------

    #[test]
    fn test_interp_t0() {
        assert!((interface_interpolation_linear(3.0, 7.0, 0.0) - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_interp_t1() {
        assert!((interface_interpolation_linear(3.0, 7.0, 1.0) - 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_interp_midpoint() {
        assert!((interface_interpolation_linear(0.0, 10.0, 0.5) - 5.0).abs() < 1e-12);
    }

    // -- domain_decompose_1d --------------------------------------------------

    #[test]
    fn test_decompose_half() {
        let (l, n) = domain_decompose_1d(10, 5);
        assert_eq!(l, 5);
        assert_eq!(n, 5);
    }

    #[test]
    fn test_decompose_all_lbm() {
        let (l, n) = domain_decompose_1d(10, 10);
        assert_eq!(l, 10);
        assert_eq!(n, 0);
    }

    #[test]
    fn test_decompose_overflow_clamped() {
        let (l, n) = domain_decompose_1d(10, 20);
        assert_eq!(l, 10);
        assert_eq!(n, 0);
    }

    #[test]
    fn test_decompose_zero_lbm() {
        let (l, n) = domain_decompose_1d(10, 0);
        assert_eq!(l, 0);
        assert_eq!(n, 10);
    }

    // -- energy_conservation_check -------------------------------------------

    #[test]
    fn test_energy_equal() {
        assert!(energy_conservation_check(1.0, 1.0, 1e-6));
    }

    #[test]
    fn test_energy_small_diff() {
        assert!(energy_conservation_check(1.0, 1.0001, 1e-2));
    }

    #[test]
    fn test_energy_large_diff() {
        assert!(!energy_conservation_check(1.0, 2.0, 1e-6));
    }

    #[test]
    fn test_energy_both_zero() {
        assert!(energy_conservation_check(0.0, 0.0, 1e-6));
    }

    // -- mass_flux_interface --------------------------------------------------

    #[test]
    fn test_mass_flux_zero_density() {
        let rho = vec![0.0, 0.0];
        let u = vec![1.0, 1.0];
        assert_eq!(mass_flux_interface(&rho, &u, [1.0, 0.0], 0.1), 0.0);
    }

    #[test]
    fn test_mass_flux_positive() {
        let rho = vec![1.0];
        let u = vec![1.0];
        let flux = mass_flux_interface(&rho, &u, [1.0, 0.0], 0.1);
        assert!(flux > 0.0);
    }

    #[test]
    fn test_mass_flux_symmetry() {
        let rho = vec![1.0, 1.0];
        let u = vec![1.0, -1.0];
        let flux = mass_flux_interface(&rho, &u, [1.0, 0.0], 1.0);
        // Velocities cancel
        assert!(flux.abs() < 1e-12);
    }
}
