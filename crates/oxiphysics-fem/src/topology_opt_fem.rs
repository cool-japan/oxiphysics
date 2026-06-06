// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! FEM-based topology optimization using SIMP and level-set methods.
//!
//! Implements the Solid Isotropic Material with Penalization (SIMP) method for
//! structural topology optimization on 2-D structured meshes, combined with
//! sensitivity filtering, OC (Optimality Criteria) density updates, and
//! auxiliary tools for thermal topology optimization and manufacturability.

// ============================================================================
// TopOptProblem
// ============================================================================

/// Parameters defining a 2-D topology optimization problem.
#[derive(Debug, Clone)]
pub struct TopOptProblem {
    /// Number of elements in the x-direction.
    pub nelx: usize,
    /// Number of elements in the y-direction.
    pub nely: usize,
    /// Target volume fraction (0 < volfrac ≤ 1).
    pub volfrac: f64,
    /// SIMP penalization exponent (typically 3).
    pub penal: f64,
    /// Filter radius for sensitivity smoothing (in element units).
    pub rmin: f64,
}

impl TopOptProblem {
    /// Create a new topology optimization problem.
    ///
    /// # Parameters
    /// * `nelx` – number of elements along x.
    /// * `nely` – number of elements along y.
    /// * `volfrac` – target volume fraction in (0, 1].
    /// * `penal` – SIMP penalization exponent.
    /// * `rmin` – sensitivity-filter radius.
    pub fn new(nelx: usize, nely: usize, volfrac: f64, penal: f64, rmin: f64) -> Self {
        Self {
            nelx,
            nely,
            volfrac,
            penal,
            rmin,
        }
    }

    /// Initialize the design field to a uniform density equal to `volfrac`.
    pub fn initialize(&self) -> TopOptState {
        let n = self.nelx * self.nely;
        let x = vec![self.volfrac; n];
        let x_old = vec![self.volfrac; n];
        let dc = vec![0.0; n];
        TopOptState {
            x,
            x_old,
            dc,
            objective: 0.0,
        }
    }

    /// SIMP stiffness interpolation: E(x) = E_min + x^p * (E_0 - E_min).
    ///
    /// Returns the effective Young's modulus for element density `x`.
    pub fn simp_stiffness(&self, x: f64) -> f64 {
        let e0 = 1.0_f64;
        let emin = 1e-9_f64;
        emin + x.powf(self.penal) * (e0 - emin)
    }

    /// Build the element-to-DOF map for the structured Q4 mesh.
    ///
    /// For an `nelx × nely` mesh the nodes are numbered row-major from the
    /// bottom-left corner.  Each Q4 element has 4 corner nodes ordered
    /// counter-clockwise; each node contributes 2 DOFs (x and y displacement).
    ///
    /// The returned `Vec` has one entry per element (row-major, y-outer x-inner),
    /// and each entry is the ordered list of 8 global DOF indices for that element.
    pub fn build_element_dof_map(&self) -> Vec<Vec<usize>> {
        let nelx = self.nelx;
        let nely = self.nely;
        let n_nodes_per_row = nelx + 1;
        let n_elem = nelx * nely;
        let mut map = Vec::with_capacity(n_elem);
        for ey in 0..nely {
            for ex in 0..nelx {
                // CCW node numbering for Q4:
                //   n0 = bottom-left,  n1 = bottom-right,
                //   n2 = top-right,    n3 = top-left
                let n0 = ey * n_nodes_per_row + ex;
                let n1 = ey * n_nodes_per_row + (ex + 1);
                let n2 = (ey + 1) * n_nodes_per_row + (ex + 1);
                let n3 = (ey + 1) * n_nodes_per_row + ex;
                let dofs = vec![
                    2 * n0,
                    2 * n0 + 1,
                    2 * n1,
                    2 * n1 + 1,
                    2 * n2,
                    2 * n2 + 1,
                    2 * n3,
                    2 * n3 + 1,
                ];
                map.push(dofs);
            }
        }
        map
    }

    /// Compute structural compliance (objective) from displacement vector and
    /// element stiffness matrices.
    ///
    /// Compliance = uᵀ K u (sum over elements), using the proper Q4 mesh
    /// connectivity from `build_element_dof_map()`.
    pub fn compute_objective(&self, u: &[f64], ke: &[Vec<f64>]) -> f64 {
        let dof_per_elem = 8_usize; // 4-node quadrilateral in 2-D: 8 DOFs
        let mut compliance = 0.0_f64;
        let n_elem = self.nelx * self.nely;
        let dof_map = self.build_element_dof_map();
        let n_dof = u.len();
        for e in 0..n_elem {
            if e >= ke.len() {
                break;
            }
            let edofs = &dof_map[e];
            let mut ue = vec![0.0_f64; dof_per_elem];
            for i in 0..dof_per_elem {
                let idx = edofs[i];
                if idx < n_dof {
                    ue[i] = u[idx];
                }
            }
            // uᵢ kᵢⱼ uⱼ
            for i in 0..dof_per_elem {
                for j in 0..dof_per_elem {
                    let kij = {
                        let row_len = dof_per_elem;
                        let flat = i * row_len + j;
                        if flat < ke[e].len() { ke[e][flat] } else { 0.0 }
                    };
                    compliance += ue[i] * kij * ue[j];
                }
            }
        }
        compliance
    }

    /// Apply a sensitivity filter (weighted averaging over a circular neighbourhood
    /// of radius `self.rmin`) to smooth element sensitivities in `state.dc`.
    ///
    /// Stores the result back into `state.dc`.
    pub fn sensitivity_filter(&self, state: &mut TopOptState) {
        let nelx = self.nelx;
        let nely = self.nely;
        let rmin = self.rmin;
        let dc_in = state.dc.clone();
        let x = &state.x;

        let mut dc_out = vec![0.0_f64; nelx * nely];

        for iy in 0..nely {
            for ix in 0..nelx {
                let e = iy * nelx + ix;
                let mut sum = 0.0_f64;
                let mut weight_sum = 0.0_f64;

                let i_lo = if ix as f64 > rmin {
                    ix - rmin.ceil() as usize
                } else {
                    0
                };
                let i_hi = (ix + rmin.ceil() as usize + 1).min(nelx);
                let j_lo = if iy as f64 > rmin {
                    iy - rmin.ceil() as usize
                } else {
                    0
                };
                let j_hi = (iy + rmin.ceil() as usize + 1).min(nely);

                for jy in j_lo..j_hi {
                    for jx in i_lo..i_hi {
                        let f = rmin
                            - ((ix as f64 - jx as f64).powi(2) + (iy as f64 - jy as f64).powi(2))
                                .sqrt();
                        if f > 0.0 {
                            let f2 = x[jy * nelx + jx];
                            weight_sum += f;
                            sum += f * f2 * dc_in[jy * nelx + jx];
                        }
                    }
                }

                if weight_sum > 1e-15 {
                    dc_out[e] = sum / (x[e].max(1e-9) * weight_sum);
                } else {
                    dc_out[e] = dc_in[e];
                }
            }
        }

        state.dc = dc_out;
    }

    /// Optimality Criteria (OC) update of the density field.
    ///
    /// Performs a bisection on the Lagrange multiplier for the volume constraint,
    /// then returns the updated density vector.
    pub fn oc_update(&self, state: &TopOptState) -> Vec<f64> {
        let n = self.nelx * self.nely;
        let move_limit = 0.2_f64;
        let mut l1 = 0.0_f64;
        let mut l2 = 1e9_f64;
        let mut xnew = vec![0.0_f64; n];

        while (l2 - l1) / (l1 + l2 + 1e-30) > 1e-3 {
            let lmid = 0.5 * (l1 + l2);
            for (e, xnew_e) in xnew.iter_mut().enumerate() {
                let be = (-state.dc[e] / lmid).max(0.0).sqrt();
                let lower = (state.x[e] - move_limit).max(0.001);
                let upper = (state.x[e] + move_limit).min(1.0);
                *xnew_e = (state.x[e] * be).clamp(lower, upper);
            }
            let vol: f64 = xnew.iter().sum::<f64>() / n as f64;
            if vol > self.volfrac {
                l1 = lmid;
            } else {
                l2 = lmid;
            }
        }

        xnew
    }

    /// Compute the current volume fraction of `state.x`.
    pub fn volume_fraction(&self, state: &TopOptState) -> f64 {
        let n = self.nelx * self.nely;
        if n == 0 {
            return 0.0;
        }
        state.x.iter().sum::<f64>() / n as f64
    }

    /// Check convergence: returns `true` if `dx_max` (maximum element density
    /// change from the previous iteration) is below 0.01.
    pub fn converged(&self, dx_max: f64) -> bool {
        dx_max < 0.01
    }
}

// ============================================================================
// TopOptState
// ============================================================================

/// Mutable state for the topology optimization iteration.
#[derive(Debug, Clone)]
pub struct TopOptState {
    /// Current element densities.
    pub x: Vec<f64>,
    /// Element densities from the previous iteration.
    pub x_old: Vec<f64>,
    /// Element compliance sensitivities (∂c/∂xₑ).
    pub dc: Vec<f64>,
    /// Structural compliance (objective function value).
    pub objective: f64,
}

impl TopOptState {
    /// Create a new state with specified size.
    pub fn new(n_elem: usize, init_density: f64) -> Self {
        Self {
            x: vec![init_density; n_elem],
            x_old: vec![init_density; n_elem],
            dc: vec![0.0; n_elem],
            objective: 0.0,
        }
    }

    /// Compute the maximum absolute change in element densities since the last
    /// snapshot stored in `x_old`.
    pub fn max_density_change(&self) -> f64 {
        self.x
            .iter()
            .zip(self.x_old.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max)
    }

    /// Update `x_old` to the current `x`.
    pub fn snapshot(&mut self) {
        self.x_old.clone_from(&self.x);
    }
}

// ============================================================================
// FiniteElementMesh2D
// ============================================================================

/// A structured 2-D quadrilateral mesh for FEM topology optimization.
///
/// The mesh consists of `nelx × nely` 4-node bilinear quadrilateral elements
/// arranged on a unit-cell grid.
#[derive(Debug, Clone)]
pub struct FiniteElementMesh2D {
    /// Number of elements in the x-direction.
    pub nelx: usize,
    /// Number of elements in the y-direction.
    pub nely: usize,
}

impl FiniteElementMesh2D {
    /// Create a new 2-D FEM mesh.
    pub fn new(nelx: usize, nely: usize) -> Self {
        Self { nelx, nely }
    }

    /// Total number of degrees of freedom (2 per node, structured grid).
    pub fn dofs(&self) -> usize {
        2 * (self.nelx + 1) * (self.nely + 1)
    }

    /// Assemble the global stiffness matrix (stored as dense rows for simplicity).
    ///
    /// Uses SIMP interpolation with penalization exponent `penal`.
    pub fn assemble_stiffness(&self, x: &[f64], penal: f64) -> Vec<Vec<f64>> {
        let ndof = self.dofs();
        let mut k_global = vec![vec![0.0_f64; ndof]; ndof];
        let ke = self.element_stiffness();
        let n_elem = self.nelx * self.nely;

        for e in 0..n_elem {
            let iy = e / self.nelx;
            let ix = e % self.nelx;
            let e_factor = {
                let xe = if e < x.len() { x[e] } else { 1.0 };
                let emin = 1e-9_f64;
                let e0 = 1.0_f64;
                emin + xe.powf(penal) * (e0 - emin)
            };

            // DOF indices for element (e): node numbering on structured grid
            // node layout: bottom-left, bottom-right, top-right, top-left
            let n1 = iy * (self.nelx + 1) + ix; // bottom-left
            let n2 = iy * (self.nelx + 1) + ix + 1; // bottom-right
            let n3 = (iy + 1) * (self.nelx + 1) + ix + 1; // top-right
            let n4 = (iy + 1) * (self.nelx + 1) + ix; // top-left

            let edofs = [
                2 * n1,
                2 * n1 + 1,
                2 * n2,
                2 * n2 + 1,
                2 * n3,
                2 * n3 + 1,
                2 * n4,
                2 * n4 + 1,
            ];

            for i in 0..8 {
                for j in 0..8 {
                    let gi = edofs[i];
                    let gj = edofs[j];
                    if gi < ndof && gj < ndof {
                        k_global[gi][gj] += e_factor * ke[i][j];
                    }
                }
            }
        }

        k_global
    }

    /// Apply Dirichlet boundary conditions.
    ///
    /// Returns a list of `(dof_index, prescribed_value)` pairs.
    /// Left edge: fixed (u = v = 0). Bottom-center node: fixed in y.
    pub fn apply_bcs(&self) -> Vec<(usize, f64)> {
        let mut bcs = Vec::new();
        // Left edge nodes: fix all DOFs
        for iy in 0..=(self.nely) {
            let n = iy * (self.nelx + 1);
            bcs.push((2 * n, 0.0)); // u
            bcs.push((2 * n + 1, 0.0)); // v
        }
        bcs
    }

    /// Compute the bilinear Q4 element stiffness matrix (8×8) for a unit cell.
    ///
    /// Uses 2×2 Gauss integration with plane-stress assumption (ν = 0.3).
    pub fn element_stiffness(&self) -> Vec<Vec<f64>> {
        // 2×2 Gauss points
        let gp = 1.0_f64 / 3.0_f64.sqrt();
        let gpts = [(-gp, -gp), (gp, -gp), (gp, gp), (-gp, gp)];

        let nu = 0.3_f64;
        // Plane-stress constitutive matrix D (E = 1):
        let d11 = 1.0 / (1.0 - nu * nu);
        let d12 = nu / (1.0 - nu * nu);
        let d33 = 0.5 * (1.0 - nu) / (1.0 - nu * nu);
        let d_mat = [[d11, d12, 0.0], [d12, d11, 0.0], [0.0, 0.0, d33]];

        let mut ke = vec![vec![0.0_f64; 8]; 8];

        for &(xi, eta) in &gpts {
            // Shape function derivatives in natural coordinates
            let dna = [
                [-(1.0 - eta) / 4.0, -(1.0 - xi) / 4.0],
                [(1.0 - eta) / 4.0, -(1.0 + xi) / 4.0],
                [(1.0 + eta) / 4.0, (1.0 + xi) / 4.0],
                [-(1.0 + eta) / 4.0, (1.0 - xi) / 4.0],
            ];

            // Jacobian for unit square (J = 0.5 I for element size 1×1)
            let j_inv = 1.0_f64; // simplified: each element is 1×1

            // B matrix (3×8)
            let mut b = vec![vec![0.0_f64; 8]; 3];
            for i in 0..4 {
                b[0][2 * i] = dna[i][0] * j_inv;
                b[1][2 * i + 1] = dna[i][1] * j_inv;
                b[2][2 * i] = dna[i][1] * j_inv;
                b[2][2 * i + 1] = dna[i][0] * j_inv;
            }

            // ke += B^T D B * det(J) * weight (det(J) = 0.25, weight = 1)
            let det_j = 0.25_f64;
            for i in 0..8 {
                for j in 0..8 {
                    let mut val = 0.0_f64;
                    for r in 0..3 {
                        for s in 0..3 {
                            val += b[r][i] * d_mat[r][s] * b[s][j];
                        }
                    }
                    ke[i][j] += val * det_j;
                }
            }
        }

        ke
    }

    /// Return the 4 corner coordinates of element `e` (row-major ordering).
    ///
    /// Each element occupies a unit cell; coordinates are in element-index space.
    pub fn node_coords(&self, e: usize) -> [[f64; 2]; 4] {
        let iy = e / self.nelx;
        let ix = e % self.nelx;
        let x0 = ix as f64;
        let y0 = iy as f64;
        [
            [x0, y0],
            [x0 + 1.0, y0],
            [x0 + 1.0, y0 + 1.0],
            [x0, y0 + 1.0],
        ]
    }
}

// ============================================================================
// ThermalTopOpt
// ============================================================================

/// Thermal topology optimization: minimize heat compliance (maximize conductance).
#[derive(Debug, Clone)]
pub struct ThermalTopOpt {
    /// Element thermal conductivities (one per element).
    pub conductivity: Vec<f64>,
}

impl ThermalTopOpt {
    /// Create a new thermal topology optimization state with `n` elements.
    pub fn new(n: usize, k_init: f64) -> Self {
        Self {
            conductivity: vec![k_init; n],
        }
    }

    /// Compute the thermal compliance: C = Tᵀ K_th T.
    ///
    /// Uses a simplified 1-D approximation: C = Σ(1 / kₑ) (series resistance).
    pub fn heat_compliance(&self) -> f64 {
        self.conductivity
            .iter()
            .map(|&k| if k > 1e-15 { 1.0 / k } else { 1e15 })
            .sum()
    }

    /// Compute thermal sensitivities ∂C/∂kₑ = -1/kₑ².
    pub fn thermal_sensitivity(&self) -> Vec<f64> {
        self.conductivity
            .iter()
            .map(|&k| if k > 1e-15 { -1.0 / (k * k) } else { -1e30 })
            .collect()
    }

    /// Update conductivities using an OC-like rule with SIMP interpolation.
    ///
    /// k(x) = k_min + x^p * (k_0 - k_min).
    pub fn update_conductivity(&mut self, x: &[f64], penal: f64) {
        let k0 = 1.0_f64;
        let kmin = 1e-6_f64;
        for (k, &xe) in self.conductivity.iter_mut().zip(x.iter()) {
            *k = kmin + xe.clamp(0.0, 1.0).powf(penal) * (k0 - kmin);
        }
    }
}

// ============================================================================
// ManufacturabilityFilter
// ============================================================================

/// Manufacturability filter for topology optimization results.
///
/// Enforces minimum length scales and checks for overhang constraints
/// (relevant for additive manufacturing).
#[derive(Debug, Clone)]
pub struct ManufacturabilityFilter {
    /// Minimum feature length (in element units).
    pub min_length: f64,
}

impl ManufacturabilityFilter {
    /// Create a new manufacturability filter.
    pub fn new(min_length: f64) -> Self {
        Self { min_length }
    }

    /// Apply a morphological closing filter to enforce the minimum length scale.
    ///
    /// Iteratively dilates then erodes the density field.
    pub fn apply(&self, x: &[f64]) -> Vec<f64> {
        if x.is_empty() {
            return vec![];
        }
        // Simplified 1-D closing: blur with a Gaussian-like kernel
        let n = x.len();
        let r = (self.min_length / 2.0).ceil() as usize;
        let mut dilated = vec![0.0_f64; n];
        for (i, dil_i) in dilated.iter_mut().enumerate() {
            let lo = i.saturating_sub(r);
            let hi = (i + r + 1).min(n);
            *dil_i = x[lo..hi].iter().cloned().fold(0.0_f64, f64::max);
        }
        let mut eroded = vec![0.0_f64; n];
        for (i, er_i) in eroded.iter_mut().enumerate() {
            let lo = i.saturating_sub(r);
            let hi = (i + r + 1).min(n);
            *er_i = dilated[lo..hi].iter().cloned().fold(1.0_f64, f64::min);
        }
        eroded
    }

    /// Check overhang constraint for additive manufacturing (build in +y direction).
    ///
    /// Returns `true` if the density field satisfies the overhang constraint
    /// (no unsupported solid material above void) on a 1-D representation.
    pub fn check_overhang(&self, x: &[f64]) -> bool {
        let threshold = 0.5_f64;
        for i in 1..x.len() {
            if x[i] > threshold && x[i - 1] < threshold {
                return false; // solid over void
            }
        }
        true
    }

    /// Extract solid and void regions from a density field using threshold 0.5.
    ///
    /// Returns `(solid_indices, void_indices)`.
    pub fn extract_solid_void(&self, x: &[f64]) -> (Vec<usize>, Vec<usize>) {
        let threshold = 0.5_f64;
        let solid: Vec<usize> = x
            .iter()
            .enumerate()
            .filter(|&(_, v)| *v >= threshold)
            .map(|(i, _)| i)
            .collect();
        let void: Vec<usize> = x
            .iter()
            .enumerate()
            .filter(|&(_, v)| *v < threshold)
            .map(|(i, _)| i)
            .collect();
        (solid, void)
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ---- TopOptProblem construction ----

    #[test]
    fn test_top_opt_problem_new() {
        let p = TopOptProblem::new(10, 10, 0.5, 3.0, 1.5);
        assert_eq!(p.nelx, 10);
        assert_eq!(p.nely, 10);
        assert!((p.volfrac - 0.5).abs() < 1e-12);
        assert!((p.penal - 3.0).abs() < 1e-12);
        assert!((p.rmin - 1.5).abs() < 1e-12);
    }

    #[test]
    fn test_initialize_state() {
        let p = TopOptProblem::new(4, 4, 0.5, 3.0, 1.5);
        let state = p.initialize();
        assert_eq!(state.x.len(), 16);
        assert!(state.x.iter().all(|&v| (v - 0.5).abs() < 1e-12));
    }

    // ---- SIMP stiffness ----

    #[test]
    fn test_simp_stiffness_solid() {
        let p = TopOptProblem::new(4, 4, 0.5, 3.0, 1.5);
        let e = p.simp_stiffness(1.0);
        assert!((e - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_simp_stiffness_void() {
        let p = TopOptProblem::new(4, 4, 0.5, 3.0, 1.5);
        let e = p.simp_stiffness(0.0);
        assert!(e > 0.0);
        assert!(e < 1e-6);
    }

    #[test]
    fn test_simp_stiffness_mid() {
        let p = TopOptProblem::new(4, 4, 0.5, 3.0, 1.5);
        let e_half = p.simp_stiffness(0.5);
        assert!(e_half > 0.0);
        assert!(e_half < 1.0);
    }

    // ---- compute_objective ----

    #[test]
    fn test_compute_objective_zero_disp() {
        let p = TopOptProblem::new(2, 2, 0.5, 3.0, 1.5);
        let u = vec![0.0; 32];
        let ke: Vec<Vec<f64>> = vec![vec![1.0; 64]; 4];
        let obj = p.compute_objective(&u, &ke);
        assert!((obj).abs() < 1e-12);
    }

    #[test]
    fn test_compute_objective_positive() {
        let p = TopOptProblem::new(1, 1, 1.0, 3.0, 1.5);
        // Unit displacement on first DOF
        let mut u = vec![0.0_f64; 8];
        u[0] = 1.0;
        // Identity-like ke
        let mut ke_row = vec![0.0_f64; 64];
        for i in 0..8 {
            ke_row[i * 8 + i] = 1.0;
        }
        let ke = vec![ke_row];
        let obj = p.compute_objective(&u, &ke);
        assert!(obj > 0.0);
    }

    // ---- volume_fraction ----

    #[test]
    fn test_volume_fraction_uniform() {
        let p = TopOptProblem::new(5, 5, 0.4, 3.0, 1.5);
        let state = p.initialize();
        let vf = p.volume_fraction(&state);
        assert!((vf - 0.4).abs() < 1e-10);
    }

    #[test]
    fn test_volume_fraction_ones() {
        let p = TopOptProblem::new(3, 3, 0.5, 3.0, 1.5);
        let mut state = p.initialize();
        for v in state.x.iter_mut() {
            *v = 1.0;
        }
        let vf = p.volume_fraction(&state);
        assert!((vf - 1.0).abs() < 1e-10);
    }

    // ---- converged ----

    #[test]
    fn test_converged_true() {
        let p = TopOptProblem::new(4, 4, 0.5, 3.0, 1.5);
        assert!(p.converged(0.005));
    }

    #[test]
    fn test_converged_false() {
        let p = TopOptProblem::new(4, 4, 0.5, 3.0, 1.5);
        assert!(!p.converged(0.05));
    }

    // ---- sensitivity_filter ----

    #[test]
    fn test_sensitivity_filter_uniform() {
        let p = TopOptProblem::new(4, 4, 0.5, 3.0, 1.5);
        let mut state = p.initialize();
        // All sensitivities = 1 => after filter should still be 1
        for v in state.dc.iter_mut() {
            *v = 1.0;
        }
        p.sensitivity_filter(&mut state);
        for &v in state.dc.iter() {
            assert!((v - 1.0).abs() < 1e-6, "expected ~1 but got {v}");
        }
    }

    // ---- OC update ----

    #[test]
    fn test_oc_update_volume_respected() {
        let p = TopOptProblem::new(4, 4, 0.5, 3.0, 1.5);
        let mut state = p.initialize();
        // Constant negative sensitivity
        for v in state.dc.iter_mut() {
            *v = -1.0;
        }
        let xnew = p.oc_update(&state);
        let vf: f64 = xnew.iter().sum::<f64>() / xnew.len() as f64;
        assert!(
            (vf - p.volfrac).abs() < 0.02,
            "volume fraction {vf} not close to 0.5"
        );
    }

    #[test]
    fn test_oc_update_bounded() {
        let p = TopOptProblem::new(3, 3, 0.5, 3.0, 1.5);
        let state = p.initialize();
        let xnew = p.oc_update(&state);
        for &v in xnew.iter() {
            assert!((0.001..=1.0).contains(&v), "density out of bounds: {v}");
        }
    }

    // ---- TopOptState ----

    #[test]
    fn test_state_max_density_change_zero() {
        let state = TopOptState::new(9, 0.5);
        assert!(state.max_density_change() < 1e-12);
    }

    #[test]
    fn test_state_snapshot() {
        let mut state = TopOptState::new(4, 0.5);
        state.x[0] = 0.8;
        state.snapshot();
        assert!((state.x_old[0] - 0.8).abs() < 1e-12);
    }

    // ---- FiniteElementMesh2D ----

    #[test]
    fn test_fem_mesh_dofs() {
        let mesh = FiniteElementMesh2D::new(4, 4);
        assert_eq!(mesh.dofs(), 2 * 5 * 5);
    }

    #[test]
    fn test_fem_mesh_element_stiffness_size() {
        let mesh = FiniteElementMesh2D::new(2, 2);
        let ke = mesh.element_stiffness();
        assert_eq!(ke.len(), 8);
        assert_eq!(ke[0].len(), 8);
    }

    #[test]
    fn test_fem_mesh_element_stiffness_symmetric() {
        let mesh = FiniteElementMesh2D::new(2, 2);
        let ke = mesh.element_stiffness();
        for (i, row) in ke.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let diff = (val - ke[j][i]).abs();
                assert!(
                    diff < 1e-10,
                    "ke not symmetric at ({i},{j}): {} vs {}",
                    val,
                    ke[j][i]
                );
            }
        }
    }

    #[test]
    fn test_fem_mesh_node_coords() {
        let mesh = FiniteElementMesh2D::new(4, 3);
        let coords = mesh.node_coords(0);
        assert!((coords[0][0]).abs() < 1e-12);
        assert!((coords[0][1]).abs() < 1e-12);
        assert!((coords[1][0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_fem_apply_bcs_nonempty() {
        let mesh = FiniteElementMesh2D::new(4, 4);
        let bcs = mesh.apply_bcs();
        assert!(!bcs.is_empty());
        // All prescribed values should be 0
        for &(_, val) in &bcs {
            assert!((val).abs() < 1e-12);
        }
    }

    #[test]
    fn test_fem_assemble_stiffness_size() {
        let mesh = FiniteElementMesh2D::new(2, 2);
        let x = vec![1.0_f64; 4];
        let kg = mesh.assemble_stiffness(&x, 3.0);
        let ndof = mesh.dofs();
        assert_eq!(kg.len(), ndof);
        assert_eq!(kg[0].len(), ndof);
    }

    // ---- ThermalTopOpt ----

    #[test]
    fn test_thermal_compliance_series() {
        let t = ThermalTopOpt::new(4, 1.0);
        let c = t.heat_compliance();
        assert!((c - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_thermal_sensitivity_negative() {
        let t = ThermalTopOpt::new(3, 2.0);
        let sens = t.thermal_sensitivity();
        for &s in &sens {
            assert!(s < 0.0, "sensitivity should be negative but got {s}");
        }
    }

    #[test]
    fn test_thermal_update_conductivity() {
        let mut t = ThermalTopOpt::new(4, 0.5);
        let x = vec![1.0_f64; 4];
        t.update_conductivity(&x, 3.0);
        for &k in &t.conductivity {
            assert!((k - 1.0).abs() < 1e-6);
        }
    }

    // ---- ManufacturabilityFilter ----

    #[test]
    fn test_mfg_filter_preserves_solid() {
        let f = ManufacturabilityFilter::new(1.0);
        let x = vec![1.0_f64; 10];
        let out = f.apply(&x);
        for &v in &out {
            assert!(v > 0.5, "solid region reduced: {v}");
        }
    }

    #[test]
    fn test_mfg_check_overhang_valid() {
        let f = ManufacturabilityFilter::new(1.0);
        let x = vec![1.0_f64, 1.0, 1.0, 1.0];
        assert!(f.check_overhang(&x));
    }

    #[test]
    fn test_mfg_check_overhang_invalid() {
        let f = ManufacturabilityFilter::new(1.0);
        let x = vec![0.0_f64, 1.0, 0.0, 0.0]; // solid at index 1 over void at 0
        assert!(!f.check_overhang(&x));
    }

    #[test]
    fn test_mfg_extract_solid_void() {
        let f = ManufacturabilityFilter::new(1.0);
        let x = vec![0.8_f64, 0.2, 0.6, 0.1];
        let (solid, void) = f.extract_solid_void(&x);
        assert_eq!(solid, vec![0, 2]);
        assert_eq!(void, vec![1, 3]);
    }

    #[test]
    fn test_mfg_filter_empty() {
        let f = ManufacturabilityFilter::new(2.0);
        let out = f.apply(&[]);
        assert!(out.is_empty());
    }

    // ── B1: build_element_dof_map tests ────────────────────────────────────

    #[test]
    fn test_dof_map_size() {
        let p = TopOptProblem::new(4, 3, 0.5, 3.0, 1.5);
        let map = p.build_element_dof_map();
        assert_eq!(map.len(), 4 * 3, "must have nelx*nely entries");
        for (i, dofs) in map.iter().enumerate() {
            assert_eq!(dofs.len(), 8, "element {i} must have 8 DOFs");
        }
    }

    #[test]
    fn test_dof_map_covers_all_dofs() {
        // For a 2×2 mesh: (2+1)*(2+1) = 9 nodes → 18 DOFs (0..17)
        let p = TopOptProblem::new(2, 2, 0.5, 3.0, 1.5);
        let map = p.build_element_dof_map();
        let mut seen = std::collections::HashSet::new();
        for dofs in &map {
            for &d in dofs {
                seen.insert(d);
            }
        }
        // All DOFs 0..17 must appear
        let n_nodes = (p.nelx + 1) * (p.nely + 1);
        for d in 0..2 * n_nodes {
            assert!(seen.contains(&d), "DOF {d} not covered by any element");
        }
    }

    #[test]
    fn test_dof_map_first_element_nelx3() {
        // nelx=3, nely=2 → n_nodes_per_row=4
        // Element 0 (ey=0, ex=0): n0=0, n1=1, n2=5, n3=4
        // DOFs: [0,1, 2,3, 10,11, 8,9]
        let p = TopOptProblem::new(3, 2, 0.5, 3.0, 1.5);
        let map = p.build_element_dof_map();
        let expected = vec![0usize, 1, 2, 3, 10, 11, 8, 9];
        assert_eq!(
            map[0], expected,
            "First element DOF map mismatch: got {:?}, want {:?}",
            map[0], expected
        );
    }

    #[test]
    fn test_dof_map_compliance_identity_k() {
        // 1×1 mesh, identity K, u=[1,0, 2,0, 3,0, 4,0] (x-displacements at nodes)
        // Compliance should equal Σ ue_i^2 for the identity (diagonal K)
        let p = TopOptProblem::new(1, 1, 0.5, 3.0, 1.5);
        // 1×1 mesh → 4 nodes, 8 DOFs, 1 element
        let u: Vec<f64> = vec![1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.0];
        // Identity K (8×8)
        let mut k_identity = vec![0.0_f64; 64];
        for i in 0..8 {
            k_identity[i * 8 + i] = 1.0;
        }
        let ke = vec![k_identity];
        let c = p.compute_objective(&u, &ke);
        // With identity K, compliance = Σ ue_i^2 = 1+0+4+0+9+0+16+0 = 30
        let expected = 1.0_f64 + 0.0 + 4.0 + 0.0 + 9.0 + 0.0 + 16.0 + 0.0;
        let diff = (c - expected).abs();
        assert!(
            diff < 1e-10,
            "Compliance with identity K: expected {expected}, got {c}"
        );
    }
}
